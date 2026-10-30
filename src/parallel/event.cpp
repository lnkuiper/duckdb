#include "duckdb/parallel/event.hpp"
#include "duckdb/common/assert.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include "duckdb/execution/executor.hpp"

namespace duckdb {

Event::Event(Executor &executor_p)
    : executor(executor_p), finished_tasks(0), total_tasks(0), finished_dependencies(0), total_dependencies(0),
      finished(false) {
}

OptionalEventDependency::OptionalEventDependency(const shared_ptr<Event> &dependent_p) : dependent(dependent_p) {
}

void OptionalEventDependency::Complete() {
	if (completed.exchange(true)) {
		return;
	}
	auto parent = dependent.lock();
	if (parent) {
		parent->CompleteDependency(PipelineDependencyType::OPTIONAL_DEPENDENCY);
	}
}

void Event::CompleteDependency(PipelineDependencyType type) {
	if (type == PipelineDependencyType::REQUIRED && optional_dependencies &&
	    ++finished_required_dependencies == total_required_dependencies) {
		auto &release =
		    CanSkipOptionalDependencies() ? optional_dependencies->prerequisites : optional_dependencies->dependents;
		// Keep this required dependency outstanding while releasing edges that can reenter scheduling.
		for (auto &dependency : release) {
			dependency->Complete();
		}
	}
	idx_t current_finished = ++finished_dependencies;
	D_ASSERT(current_finished <= total_dependencies);
	if (current_finished == total_dependencies) {
		D_ASSERT(total_tasks == 0);
		Schedule();
		if (total_tasks == 0 && AutoFinishWithoutTasks()) {
			Finish();
		}
	}
}

void Event::SetOptionalDependencies(unique_ptr<OptionalEventDependencies> dependencies) {
	D_ASSERT(total_required_dependencies > 0);
	D_ASSERT(!optional_dependencies);
	optional_dependencies = std::move(dependencies);
}

void Event::Finish() {
	D_ASSERT(!finished);
	FinishEvent();
	finished = true;
	// finished processing the pipeline, now we can schedule pipelines that depend on this pipeline
	for (auto &parent_entry : parents) {
		if (parent_entry.optional_dependency) {
			parent_entry.optional_dependency->Complete();
			continue;
		}
		auto parent = parent_entry.event.lock();
		if (!parent) { // LCOV_EXCL_START
			continue;
		} // LCOV_EXCL_STOP
		// mark a dependency as completed for each of the parents
		parent->CompleteDependency();
	}
	FinalizeFinish();
}

shared_ptr<OptionalEventDependency> Event::AddDependency(Event &event, PipelineDependencyType type) {
	total_dependencies++;
	shared_ptr<OptionalEventDependency> optional_dependency;
	if (type == PipelineDependencyType::OPTIONAL_DEPENDENCY) {
		optional_dependency = make_shared_ptr<OptionalEventDependency>(shared_from_this());
	} else {
		total_required_dependencies++;
	}
	event.parents.push_back({weak_ptr<Event>(shared_from_this()), optional_dependency});
#ifdef DEBUG
	if (type == PipelineDependencyType::REQUIRED) {
		event.parents_raw.push_back(*this);
	}
#endif
	return optional_dependency;
}

const vector<reference<Event>> &Event::GetParentsVerification() const {
	return parents_raw;
}

void Event::FinishTask() {
	D_ASSERT(finished_tasks.load() < total_tasks.load());
	idx_t current_tasks = total_tasks;
	idx_t current_finished = ++finished_tasks;
	D_ASSERT(current_finished <= current_tasks);
	if (current_finished == current_tasks) {
		Finish();
	}
}

ClientContext &Event::GetClientContext() {
	return executor.context;
}

void Event::InsertEvent(shared_ptr<Event> replacement_event) {
	replacement_event->parents = std::move(parents);
#ifdef DEBUG
	replacement_event->parents_raw = std::move(parents_raw);
#endif
	replacement_event->AddDependency(*this);
	executor.AddEvent(std::move(replacement_event));
}

void Event::SetTasks(vector<shared_ptr<Task>> tasks) {
	auto &ts = TaskScheduler::GetScheduler(executor.context);
	D_ASSERT(total_tasks == 0);
	D_ASSERT(!tasks.empty());
	this->total_tasks = tasks.size();
	ts.ScheduleTasks(executor.GetToken(), tasks);
}

} // namespace duckdb
