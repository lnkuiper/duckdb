//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/parallel/event.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/execution/executor.hpp"
#include "duckdb/common/atomic.hpp"
#include "duckdb/common/common.hpp"
#include "duckdb/common/vector.hpp"

namespace duckdb {
class Executor;
class Task;
class Event;

class OptionalEventDependency {
public:
	explicit OptionalEventDependency(const shared_ptr<Event> &dependent);
	void Complete();

private:
	weak_ptr<Event> dependent;
	atomic<bool> completed {false};
};

struct EventDependent {
	weak_ptr<Event> event;
	shared_ptr<OptionalEventDependency> optional_dependency;
};

struct OptionalEventDependencies {
	vector<shared_ptr<OptionalEventDependency>> prerequisites;
	vector<shared_ptr<OptionalEventDependency>> dependents;
};

class Event : public enable_shared_from_this<Event> {
public:
	explicit Event(Executor &executor);
	virtual ~Event() = default;

public:
	virtual void Schedule() = 0;
	//! Whether an event with no scheduled tasks should automatically finish after scheduling
	virtual bool AutoFinishWithoutTasks() const {
		return true;
	}
	//! Called right after the event is finished
	virtual void FinishEvent() {
	}
	//! Called after the event is entirely finished
	virtual void FinalizeFinish() {
	}
	//! Decide whether to bypass the group's prerequisites or release its optional dependents.
	virtual bool CanSkipOptionalDependencies() {
		return false;
	}

	void FinishTask();
	void Finish();

	shared_ptr<OptionalEventDependency> AddDependency(Event &event,
	                                                  PipelineDependencyType type = PipelineDependencyType::REQUIRED);
	bool HasDependencies() const {
		return total_dependencies != 0;
	}
	bool HasTasks() const {
		return total_tasks != 0;
	}
	const vector<reference<Event>> &GetParentsVerification() const;

	void CompleteDependency(PipelineDependencyType type = PipelineDependencyType::REQUIRED);
	void SetOptionalDependencies(unique_ptr<OptionalEventDependencies> dependencies);

	void SetTasks(vector<shared_ptr<Task>> tasks);

	void InsertEvent(shared_ptr<Event> replacement_event);

	bool IsFinished() const {
		return finished;
	}

	virtual void PrintPipeline() {
	}

	ClientContext &GetClientContext();

	template <class TARGET>
	TARGET &Cast() {
		DynamicCastCheck<TARGET>(this);
		return reinterpret_cast<TARGET &>(*this);
	}
	template <class TARGET>
	const TARGET &Cast() const {
		DynamicCastCheck<TARGET>(this);
		return reinterpret_cast<const TARGET &>(*this);
	}

protected:
	Executor &executor;
	//! The current threads working on the event
	atomic<idx_t> finished_tasks;
	//! The maximum amount of threads that can work on the event
	atomic<idx_t> total_tasks;

	//! The event starts after every dependency has finished or been released.
	atomic<idx_t> finished_dependencies;
	idx_t total_dependencies;
	atomic<idx_t> finished_required_dependencies {0};
	idx_t total_required_dependencies = 0;
	unique_ptr<OptionalEventDependencies> optional_dependencies;

	//! The events that depend on this event to run
	vector<EventDependent> parents;
	//! Required edges must remain acyclic even when optional edges can be released.
	vector<reference<Event>> parents_raw;

	//! Whether or not the event is finished executing
	atomic<bool> finished;
};

} // namespace duckdb
