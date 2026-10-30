#include "catch.hpp"
#include "duckdb.hpp"
#include "duckdb/parallel/event.hpp"
#include "duckdb/parallel/pipeline_schedule.hpp"

#include <thread>

using namespace duckdb;

TEST_CASE("Runtime filter ordering respects transitive and sibling dependencies", "[runtime_filter_schedule]") {
	PipelineSchedule schedule;
	for (idx_t i = 0; i < 7; i++) {
		schedule.stages.emplace_back(PipelineScheduleStageType::EXECUTE, nullptr);
	}
	for (idx_t i = 1; i <= 4; i++) {
		schedule.stages[i].dependencies.push_back(i - 1);
	}
	// A sibling build must execute before this producer prepares its finalization.
	schedule.stages[2].dependencies.push_back(6);
	REQUIRE(schedule.TryAddDependency(5, 4));
	REQUIRE(schedule.DependsOn(5, 0));
	REQUIRE(schedule.TryAddDependency(5, 2));
	REQUIRE(schedule.stages[5].dependencies.size() == 1);
	REQUIRE_FALSE(schedule.TryAddDependency(6, 4));
	REQUIRE(schedule.stages[6].dependencies.empty());
	REQUIRE_FALSE(schedule.TryAddDependency(4, 4));
	REQUIRE_FALSE(schedule.HasCycle());
}

namespace {

class DependencyTestEvent : public Event {
public:
	explicit DependencyTestEvent(Executor &executor) : Event(executor) {
	}

	bool bypass = false;
	bool auto_finish = true;
	atomic<idx_t> decisions {0};
	atomic<idx_t> schedules {0};

	bool CanSkipOptionalDependencies() override {
		decisions++;
		return bypass;
	}

	void Schedule() override {
		schedules++;
	}

	bool AutoFinishWithoutTasks() const override {
		return auto_finish;
	}

	void Start() {
		D_ASSERT(!HasDependencies());
		Schedule();
		if (auto_finish) {
			Finish();
		}
	}
};

} // namespace

TEST_CASE("Required scheduling waits survive optional path release", "[runtime_filter_schedule]") {
	DuckDB db(nullptr);
	Connection con(db);
	Executor executor(*con.context);
	executor.Reset();
	PipelineSchedule schedule;
	vector<shared_ptr<DependencyTestEvent>> events;
	for (idx_t i = 0; i < 3; i++) {
		schedule.stages.emplace_back(PipelineScheduleStageType::EXECUTE, nullptr);
		events.push_back(make_shared_ptr<DependencyTestEvent>(executor));
	}
	schedule.stages[0].dependencies.emplace_back(1);
	schedule.stages[1].dependencies.emplace_back(2, PipelineDependencyType::OPTIONAL_DEPENDENCY);
	REQUIRE(schedule.TryAddDependency(0, 2));
	shared_ptr<OptionalEventDependency> optional;
	for (idx_t i = 0; i < schedule.stages.size(); i++) {
		for (auto &dependency : schedule.stages[i].dependencies) {
			auto edge = events[i]->AddDependency(*events[dependency.index], dependency.type);
			if (edge) {
				optional = std::move(edge);
			}
		}
	}
	optional->Complete();
	REQUIRE(events[1]->IsFinished());
	REQUIRE_FALSE(events[0]->IsFinished());
	events[2]->Start();
	REQUIRE(events[0]->IsFinished());
	REQUIRE(events[0]->schedules == 1);
}

TEST_CASE("Optional admission rejects completion and finalization targets", "[runtime_filter_schedule]") {
	PipelineSchedule schedule;
	for (auto type : {PipelineScheduleStageType::PREPARE_FINISH, PipelineScheduleStageType::FINISH,
	                  PipelineScheduleStageType::COMPLETE, PipelineScheduleStageType::INITIALIZE}) {
		schedule.stages.emplace_back(type, nullptr);
	}
	schedule.stages[0].dependencies.emplace_back(3);
	schedule.stages[1].dependencies.emplace_back(0);
	schedule.stages[2].dependencies.emplace_back(1);
	schedule.stages[2].dependency_group = 0;
	REQUIRE_FALSE(schedule.TryAddOptionalDependency(2, 2));
	REQUIRE_FALSE(schedule.TryAddOptionalDependency(1, 2));
	REQUIRE_FALSE(schedule.TryAddOptionalDependency(0, 2));
	REQUIRE_FALSE(schedule.HasCycle());
}

TEST_CASE("Optional event dependencies remain enforced by default", "[runtime_filter_schedule]") {
	DuckDB db(nullptr);
	Connection con(db);
	Executor executor(*con.context);
	executor.Reset();
	auto predecessor = make_shared_ptr<DependencyTestEvent>(executor);
	auto dependent = make_shared_ptr<DependencyTestEvent>(executor);
	dependent->AddDependency(*predecessor, PipelineDependencyType::OPTIONAL_DEPENDENCY);
	REQUIRE(dependent->schedules == 0);
	predecessor->Start();
	REQUIRE(dependent->IsFinished());
	REQUIRE(dependent->schedules == 1);
	REQUIRE(dependent->decisions == 0);
}

TEST_CASE("Optional bypass waits for every required input and tolerates late completion", "[runtime_filter_schedule]") {
	DuckDB db(nullptr);
	Connection con(db);
	Executor executor(*con.context);
	executor.Reset();
	auto first = make_shared_ptr<DependencyTestEvent>(executor);
	auto second = make_shared_ptr<DependencyTestEvent>(executor);
	auto optional = make_shared_ptr<DependencyTestEvent>(executor);
	auto dependent = make_shared_ptr<DependencyTestEvent>(executor);
	dependent->bypass = true;
	dependent->AddDependency(*first);
	dependent->AddDependency(*second);
	auto group = make_uniq<OptionalEventDependencies>();
	group->prerequisites.push_back(dependent->AddDependency(*optional, PipelineDependencyType::OPTIONAL_DEPENDENCY));
	dependent->SetOptionalDependencies(std::move(group));
	first->Start();
	REQUIRE(dependent->decisions == 0);
	second->Start();
	REQUIRE(dependent->IsFinished());
	optional->Start();
	REQUIRE(dependent->decisions == 1);
	REQUIRE(dependent->schedules == 1);
}

TEST_CASE("One dependency decision resolves both barriers or releases an optional dependent",
          "[runtime_filter_schedule]") {
	DuckDB db(nullptr);
	Connection con(db);
	Executor executor(*con.context);
	executor.Reset();
	auto input = make_shared_ptr<DependencyTestEvent>(executor);
	auto prepare = make_shared_ptr<DependencyTestEvent>(executor);
	auto finish = make_shared_ptr<DependencyTestEvent>(executor);
	auto publication = make_shared_ptr<DependencyTestEvent>(executor);
	auto consumer = make_shared_ptr<DependencyTestEvent>(executor);
	auto sibling_input = make_shared_ptr<DependencyTestEvent>(executor);
	auto sibling_prepare = make_shared_ptr<DependencyTestEvent>(executor);
	SECTION("Bypass both prerequisites") {
		prepare->bypass = true;
	}
	SECTION("Release dependent before waiting for prerequisites") {
		prepare->bypass = false;
	}
	prepare->AddDependency(*input);
	finish->AddDependency(*prepare);
	publication->AddDependency(*finish);
	sibling_input->AddDependency(*consumer);
	sibling_prepare->AddDependency(*sibling_input);
	auto group = make_uniq<OptionalEventDependencies>();
	group->prerequisites.push_back(prepare->AddDependency(*sibling_input, PipelineDependencyType::OPTIONAL_DEPENDENCY));
	group->prerequisites.push_back(
	    finish->AddDependency(*sibling_prepare, PipelineDependencyType::OPTIONAL_DEPENDENCY));
	group->dependents.push_back(consumer->AddDependency(*publication, PipelineDependencyType::OPTIONAL_DEPENDENCY));
	prepare->SetOptionalDependencies(std::move(group));
	input->Start();
	REQUIRE(prepare->IsFinished());
	REQUIRE(finish->IsFinished());
	REQUIRE(publication->IsFinished());
	REQUIRE(consumer->IsFinished());
	REQUIRE(sibling_prepare->IsFinished());
	REQUIRE(prepare->decisions == 1);
	REQUIRE(prepare->schedules == 1);
	REQUIRE(finish->decisions == 0);
	REQUIRE(finish->schedules == 1);
	REQUIRE(consumer->schedules == 1);
}

TEST_CASE("Optional completion and bypass race schedules an event once", "[runtime_filter_schedule]") {
	DuckDB db(nullptr);
	Connection con(db);
	Executor executor(*con.context);
	executor.Reset();
	for (idx_t iteration = 0; iteration < 64; iteration++) {
		auto first = make_shared_ptr<DependencyTestEvent>(executor);
		auto second = make_shared_ptr<DependencyTestEvent>(executor);
		auto optional = make_shared_ptr<DependencyTestEvent>(executor);
		auto dependent = make_shared_ptr<DependencyTestEvent>(executor);
		dependent->bypass = true;
		dependent->AddDependency(*first);
		dependent->AddDependency(*second);
		auto group = make_uniq<OptionalEventDependencies>();
		auto edge = dependent->AddDependency(*optional, PipelineDependencyType::OPTIONAL_DEPENDENCY);
		group->prerequisites.push_back(edge);
		dependent->SetOptionalDependencies(std::move(group));
		std::thread t1([&] { first->Start(); });
		std::thread t2([&] { second->Start(); });
		std::thread t3([&] { optional->Start(); });
		edge->Complete();
		t1.join();
		t2.join();
		t3.join();
		REQUIRE(dependent->IsFinished());
		REQUIRE(dependent->decisions == 1);
		REQUIRE(dependent->schedules == 1);
	}
}

TEST_CASE("Inserted finalization events retain optional edge completion", "[runtime_filter_schedule]") {
	DuckDB db(nullptr);
	Connection con(db);
	Executor executor(*con.context);
	executor.Reset();
	auto producer = make_shared_ptr<DependencyTestEvent>(executor);
	auto replacement = make_shared_ptr<DependencyTestEvent>(executor);
	auto consumer = make_shared_ptr<DependencyTestEvent>(executor);
	producer->auto_finish = false;
	replacement->auto_finish = false;
	auto edge = consumer->AddDependency(*producer, PipelineDependencyType::OPTIONAL_DEPENDENCY);
	producer->InsertEvent(replacement);
	producer->Start();
	producer->Finish();
	REQUIRE(replacement->schedules == 1);
	REQUIRE(consumer->schedules == 0);
	SECTION("Wait for finalization") {
		replacement->Finish();
	}
	SECTION("Release dependent while finalization is running") {
		edge->Complete();
		REQUIRE(consumer->IsFinished());
		replacement->Finish();
	}
	REQUIRE(consumer->IsFinished());
	REQUIRE(consumer->schedules == 1);
}

TEST_CASE("Optional edges do not retain an abandoned consumer", "[runtime_filter_schedule]") {
	DuckDB db(nullptr);
	Connection con(db);
	Executor executor(*con.context);
	executor.Reset();
	auto producer = make_shared_ptr<DependencyTestEvent>(executor);
	auto consumer = make_shared_ptr<DependencyTestEvent>(executor);
	auto edge = consumer->AddDependency(*producer, PipelineDependencyType::OPTIONAL_DEPENDENCY);
	weak_ptr<Event> weak_consumer = consumer;
	consumer.reset();
	REQUIRE(weak_consumer.expired());
	edge->Complete();
	producer->Start();
	REQUIRE(producer->IsFinished());
}
