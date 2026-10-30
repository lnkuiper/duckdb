#include "duckdb/parallel/pipeline_schedule.hpp"

#include "duckdb/execution/operator/scan/physical_table_scan.hpp"
#include "duckdb/parallel/meta_pipeline.hpp"

namespace duckdb {

struct PipelineScheduleStageStack {
	PipelineScheduleStageStack(idx_t initialize_p, idx_t execute_p, idx_t prepare_finish_p, idx_t finish_p,
	                           idx_t complete_p)
	    : initialize(initialize_p), execute(execute_p), prepare_finish(prepare_finish_p), finish(finish_p),
	      complete(complete_p) {
	}

	idx_t initialize;
	idx_t execute;
	idx_t prepare_finish;
	idx_t finish;
	idx_t complete;
};

using pipeline_schedule_stage_map_t = reference_map_t<Pipeline, PipelineScheduleStageStack>;

PipelineScheduleStage::PipelineScheduleStage(PipelineScheduleStageType type_p, shared_ptr<Pipeline> pipeline_p)
    : type(type_p), pipeline(std::move(pipeline_p)) {
}

static idx_t AddStage(PipelineSchedule &result, PipelineScheduleStageType type, const shared_ptr<Pipeline> &pipeline) {
	result.stages.emplace_back(type, pipeline);
	return result.stages.size() - 1;
}

static void AddDependency(PipelineSchedule &result, idx_t dependent, idx_t dependency,
                          PipelineDependencyType type = PipelineDependencyType::REQUIRED) {
	result.stages[dependent].dependencies.emplace_back(dependency, type);
}

static void AddExternalInputDependency(PipelineSchedule &result, idx_t dependent, idx_t dependency,
                                       Pipeline &consumer) {
	auto dependency_index = result.stages[dependent].dependencies.size();
	result.stages[dependent].dependencies.push_back(dependency);
	result.external_input_dependencies.emplace_back(dependent, dependency_index, consumer);
}

static optional_ptr<Pipeline> GetExternalInputConsumer(const PipelineSchedule &schedule, idx_t dependent,
                                                       idx_t dependency_index) {
	for (auto &dependency : schedule.external_input_dependencies) {
		if (dependency.dependent == dependent && dependency.dependency_index == dependency_index) {
			return dependency.consumer.get();
		}
	}
	return nullptr;
}

enum class CycleVisitState : uint8_t { UNVISITED, VISITING, VISITED };

vector<PipelineScheduleEdge> PipelineSchedule::GetCycle() const {
	vector<CycleVisitState> states(stages.size(), CycleVisitState::UNVISITED);
	for (idx_t stage_idx = 0; stage_idx < stages.size(); stage_idx++) {
		if (states[stage_idx] != CycleVisitState::UNVISITED) {
			continue;
		}
		vector<pair<idx_t, idx_t>> stack;
		stack.emplace_back(stage_idx, 0);
		states[stage_idx] = CycleVisitState::VISITING;
		while (!stack.empty()) {
			auto &entry = stack.back();
			auto &dependencies = stages[entry.first].dependencies;
			if (entry.second == dependencies.size()) {
				states[entry.first] = CycleVisitState::VISITED;
				stack.pop_back();
				continue;
			}
			auto dependency = dependencies[entry.second++].index;
			D_ASSERT(dependency < stages.size());
			if (states[dependency] == CycleVisitState::UNVISITED) {
				states[dependency] = CycleVisitState::VISITING;
				stack.emplace_back(dependency, 0);
				continue;
			}
			if (states[dependency] != CycleVisitState::VISITING) {
				continue;
			}
			vector<PipelineScheduleEdge> result;
			idx_t cycle_start = 0;
			while (stack[cycle_start].first != dependency) {
				cycle_start++;
				D_ASSERT(cycle_start < stack.size());
			}
			for (idx_t path_idx = cycle_start; path_idx + 1 < stack.size(); path_idx++) {
				auto dependent = stack[path_idx].first;
				auto dependency_idx = stack[path_idx].second - 1;
				result.emplace_back(dependent, stages[dependent].dependencies[dependency_idx].index,
				                    GetExternalInputConsumer(*this, dependent, dependency_idx));
			}
			auto dependent = stack.back().first;
			auto dependency_idx = stack.back().second - 1;
			result.emplace_back(dependent, dependency, GetExternalInputConsumer(*this, dependent, dependency_idx));
			return result;
		}
	}
	return {};
}

bool PipelineSchedule::HasCycle() const {
	return !GetCycle().empty();
}

bool PipelineSchedule::DependsOn(idx_t dependent, idx_t dependency, PipelineDependencyType type) const {
	D_ASSERT(dependent < stages.size() && dependency < stages.size());
	vector<bool> visited(stages.size(), false);
	vector<idx_t> pending {dependent};
	while (!pending.empty()) {
		auto current = pending.back();
		pending.pop_back();
		if (current == dependency) {
			return true;
		}
		if (visited[current]) {
			continue;
		}
		visited[current] = true;
		for (auto &child : stages[current].dependencies) {
			if (type != PipelineDependencyType::REQUIRED || child.type == PipelineDependencyType::REQUIRED) {
				pending.push_back(child.index);
			}
		}
	}
	return false;
}

bool PipelineSchedule::TryAddDependency(idx_t dependent, idx_t dependency) {
	if ((!decision_inputs.empty() && decision_inputs[dependent]) || dependent == dependency ||
	    DependsOn(dependency, dependent)) {
		return false;
	}
	if (!DependsOn(dependent, dependency, PipelineDependencyType::REQUIRED)) {
		stages[dependent].dependencies.push_back(dependency);
	}
	return true;
}

bool PipelineSchedule::TryAddOptionalDependency(idx_t dependent, idx_t dependency) {
	auto group = stages[dependency].dependency_group;
	if (dependent == dependency || stages[dependent].type != PipelineScheduleStageType::EXECUTE || !group.IsValid() ||
	    stages[dependency].type != PipelineScheduleStageType::COMPLETE ||
	    (!decision_inputs.empty() && decision_inputs[dependent])) {
		return false;
	}
	auto &decision = stages[group.GetIndex()];
	vector<idx_t> pending;
	for (auto &input : decision.dependencies) {
		if (input.type == PipelineDependencyType::REQUIRED) {
			pending.push_back(input.index);
		}
	}
	vector<bool> inputs(stages.size(), false);
	while (!pending.empty()) {
		auto index = pending.back();
		pending.pop_back();
		if (index == dependent) {
			return false;
		}
		if (inputs[index]) {
			continue;
		}
		inputs[index] = true;
		auto &stage = stages[index];
		if (stage.type != PipelineScheduleStageType::INITIALIZE && stage.type != PipelineScheduleStageType::EXECUTE) {
			return false;
		}
		if (stage.type == PipelineScheduleStageType::EXECUTE) {
			auto source = stage.pipeline->GetSource();
			// Blocking and streaming inputs can wait for other stages before reaching the decision.
			if (stage.pipeline->IsExternalInput() || !source || source->type != PhysicalOperatorType::TABLE_SCAN ||
			    stage.pipeline->GetSink() != decision.pipeline->GetSink()) {
				return false;
			}
		}
		for (auto &input : stage.dependencies) {
			pending.push_back(input.index);
		}
	}
	decision_inputs.resize(stages.size(), false);
	for (idx_t i = 0; i < inputs.size(); i++) {
		decision_inputs[i] = decision_inputs[i] || inputs[i];
	}
	stages[dependent].dependencies.emplace_back(dependency, PipelineDependencyType::OPTIONAL_DEPENDENCY);
	return true;
}

static PipelineScheduleStageStack AddBasePipeline(PipelineSchedule &result, const shared_ptr<Pipeline> &pipeline) {
	auto initialize = AddStage(result, PipelineScheduleStageType::INITIALIZE, pipeline);
	auto execute = AddStage(result, PipelineScheduleStageType::EXECUTE, pipeline);
	auto prepare_finish = AddStage(result, PipelineScheduleStageType::PREPARE_FINISH, pipeline);
	auto finish = AddStage(result, PipelineScheduleStageType::FINISH, pipeline);
	auto complete = AddStage(result, PipelineScheduleStageType::COMPLETE, pipeline);
	// Dependencies: initialize -> execute -> prepare finish -> finish -> complete.
	AddDependency(result, execute, initialize);
	AddDependency(result, prepare_finish, execute);
	AddDependency(result, finish, prepare_finish);
	AddDependency(result, complete, finish);
	for (auto index : {prepare_finish, finish, complete}) {
		result.stages[index].dependency_group = prepare_finish;
	}
	return PipelineScheduleStageStack(initialize, execute, prepare_finish, finish, complete);
}

static bool RequiresInitializeOnSchedule(Pipeline &pipeline) {
	auto source = pipeline.GetSource();
	if (!source || source->type != PhysicalOperatorType::TABLE_SCAN) {
		return false;
	}
	auto &table_function = source->Cast<PhysicalTableScan>();
	return table_function.function.global_initialization == TableFunctionInitialization::INITIALIZE_ON_SCHEDULE;
}

static void AddMetaPipeline(PipelineSchedule &result, MetaPipeline &meta_pipeline,
                            pipeline_schedule_stage_map_t &stage_map) {
	vector<shared_ptr<Pipeline>> pipelines;
	meta_pipeline.GetPipelines(pipelines, false);
	D_ASSERT(!pipelines.empty());

	auto base_stack = AddBasePipeline(result, pipelines[0]);
	stage_map.emplace(*pipelines[0], base_stack);

	for (idx_t pipeline_idx = 1; pipeline_idx < pipelines.size(); pipeline_idx++) {
		auto &pipeline = pipelines[pipeline_idx];
		auto execute = AddStage(result, PipelineScheduleStageType::EXECUTE, pipeline);
		auto finish_group = meta_pipeline.GetFinishGroup(*pipeline);
		if (finish_group) {
			// This pipeline is part of a finish group.
			auto group_entry = stage_map.find(*finish_group);
			D_ASSERT(group_entry != stage_map.end());
			auto &group_stack = group_entry->second;
			// Dependencies: base finish -> pipeline execute -> group prepare finish.
			AddDependency(result, execute, base_stack.finish);
			AddDependency(result, group_stack.prepare_finish, execute);
			stage_map.emplace(*pipeline,
			                  PipelineScheduleStageStack(base_stack.initialize, execute, group_stack.prepare_finish,
			                                             group_stack.finish, base_stack.complete));
		} else if (meta_pipeline.HasFinishEvent(*pipeline)) {
			// This pipeline has its own finish event despite writing to the shared sink, so Finalize runs twice.
			auto prepare_finish = AddStage(result, PipelineScheduleStageType::PREPARE_FINISH, pipeline);
			auto finish = AddStage(result, PipelineScheduleStageType::FINISH, pipeline);
			// Dependencies:
			// base finish -> pipeline execute -> pipeline prepare finish -> pipeline finish -> base complete.
			AddDependency(result, execute, base_stack.finish);
			AddDependency(result, prepare_finish, execute);
			AddDependency(result, finish, prepare_finish);
			AddDependency(result, base_stack.complete, finish);
			result.stages[prepare_finish].dependency_group = prepare_finish;
			result.stages[finish].dependency_group = prepare_finish;
			// Completion now covers more than one independently finalized group.
			result.stages[base_stack.complete].dependency_group.SetInvalid();
			stage_map.emplace(*pipeline, PipelineScheduleStageStack(base_stack.initialize, execute, prepare_finish,
			                                                        finish, base_stack.complete));
		} else {
			// No additional finish event: execute after base initialization and join the base finish chain.
			AddDependency(result, execute, base_stack.initialize);
			AddDependency(result, base_stack.prepare_finish, execute);
			stage_map.emplace(*pipeline,
			                  PipelineScheduleStageStack(base_stack.initialize, execute, base_stack.prepare_finish,
			                                             base_stack.finish, base_stack.complete));
		}
	}

	for (auto &pipeline : pipelines) {
		if (RequiresInitializeOnSchedule(*pipeline)) {
			// Certain functions must be eagerly initialized during scheduling. Record them here so the schedule
			// consumer initializes the function before it starts scheduling events.
			result.initialize_on_schedule_pipelines.push_back(*pipeline);
		}
	}
}

unique_ptr<PipelineSchedule> BuildPipelineSchedule(const vector<shared_ptr<MetaPipeline>> &meta_pipelines,
                                                   PipelineScheduleMode mode) {
	auto result = make_uniq<PipelineSchedule>();
	pipeline_schedule_stage_map_t stage_map;
	const bool allow_missing_meta_pipelines = mode == PipelineScheduleMode::ACTIVE_SUBSET;

	for (auto &meta_pipeline : meta_pipelines) {
		AddMetaPipeline(*result, *meta_pipeline, stage_map);
	}

	for (auto &entry : stage_map) {
		auto &pipeline = entry.first.get();
		// A pipeline dependency is complete only after its sink has finalized.
		for (auto &dependency : pipeline.GetDependencies()) {
			auto dep = dependency.lock();
			D_ASSERT(dep);
			auto dep_entry = stage_map.find(*dep);
			if (dep_entry != stage_map.end()) {
				AddDependency(*result, entry.second.execute, dep_entry->second.complete);
			}
		}
		// A dataflow dependency only needs the producer's source and sink state to be initialized.
		for (auto &dependency : pipeline.GetDataflowDependencies()) {
			auto dep = dependency.lock();
			D_ASSERT(dep);
			auto dep_entry = stage_map.find(*dep);
			if (dep_entry != stage_map.end()) {
				AddDependency(*result, entry.second.execute, dep_entry->second.initialize);
			}
		}
	}

	// Dependencies within a MetaPipeline order their execute stages directly.
	for (auto &entry : stage_map) {
		for (auto &intra_dependency : entry.first.get().GetIntraDependencies()) {
			auto dep = intra_dependency.pipeline.lock();
			D_ASSERT(dep);
			auto dependency_entry = stage_map.find(*dep);
			if (dependency_entry == stage_map.end()) {
				if (!allow_missing_meta_pipelines) {
					throw InternalException("Missing dependency pipeline in meta-pipeline dependency");
				}
				continue;
			}
			AddDependency(*result, entry.second.execute, dependency_entry->second.execute, intra_dependency.type);
		}
	}

	// External producers execute in place of the consumer's source. They wait until the consumer state and its
	// blocking inputs are ready, while independent scheduling constraints remain attached to the consumer itself.
	// The consumer execution and finalization stages remain open until every producer has finished pushing.
	for (auto &entry : stage_map) {
		auto &consumer = entry.first.get();
		if (!consumer.IsExternalInput()) {
			continue;
		}
		for (auto &producer_ref : consumer.GetExternalInputProducers()) {
			auto producer = producer_ref.lock();
			D_ASSERT(producer);
			auto producer_entry = stage_map.find(*producer);
			if (producer_entry == stage_map.end()) {
				if (!allow_missing_meta_pipelines) {
					throw InternalException("Missing external input producer pipeline");
				}
				continue;
			}
			AddExternalInputDependency(*result, producer_entry->second.execute, entry.second.initialize, consumer);
			for (auto &dependency : consumer.GetDependencies()) {
				auto dep = dependency.lock();
				D_ASSERT(dep);
				auto dep_entry = stage_map.find(*dep);
				if (dep_entry != stage_map.end()) {
					AddExternalInputDependency(*result, producer_entry->second.execute, dep_entry->second.complete,
					                           consumer);
				}
			}
			for (auto &dependency : consumer.GetDataflowDependencies()) {
				auto dep = dependency.lock();
				D_ASSERT(dep);
				auto dep_entry = stage_map.find(*dep);
				if (dep_entry != stage_map.end()) {
					AddExternalInputDependency(*result, producer_entry->second.execute, dep_entry->second.initialize,
					                           consumer);
				}
			}
			AddExternalInputDependency(*result, entry.second.execute, producer_entry->second.execute, consumer);
			AddExternalInputDependency(*result, entry.second.prepare_finish, producer_entry->second.execute, consumer);
		}
	}

	// These dependencies make coordinated sibling join builds run in this order:
	// 1. all child pipelines execute through Combine
	// 2. all child pipelines run PrepareFinalize
	// 3. all child pipelines run Finalize
	// Operators communicate their memory usage through the TemporaryMemoryManager in PrepareFinalize. When the
	// child pipelines Finalize, all required memory is known and the manager can make an informed decision.
	for (auto &meta_pipeline : meta_pipelines) {
		vector<shared_ptr<MetaPipeline>> children;
		meta_pipeline->GetMetaPipelines(children, false, true);
		for (auto &child1 : children) {
			if (child1->Type() != MetaPipelineType::JOIN_BUILD) {
				continue;
			}
			auto child1_entry = stage_map.find(*child1->GetBasePipeline());
			if (child1_entry == stage_map.end()) {
				if (!allow_missing_meta_pipelines) {
					throw InternalException("Missing sibling join-build pipeline");
				}
				continue;
			}
			for (auto &child2 : children) {
				if (child2->Type() != MetaPipelineType::JOIN_BUILD || RefersToSameObject(*child1, *child2)) {
					continue;
				}
				if (!RefersToSameObject(*child1->GetParent(), *child2->GetParent())) {
					continue;
				}
				auto child2_entry = stage_map.find(*child2->GetBasePipeline());
				if (child2_entry == stage_map.end()) {
					if (!allow_missing_meta_pipelines) {
						throw InternalException("Missing sibling join-build dependency pipeline");
					}
					continue;
				}
				// Include every input feeding the sibling's prepare event.
				for (auto &input : result->stages[child2_entry->second.prepare_finish].dependencies) {
					if (input.type == PipelineDependencyType::REQUIRED) {
						AddDependency(*result, child1_entry->second.prepare_finish, input.index,
						              PipelineDependencyType::OPTIONAL_DEPENDENCY);
					}
				}
				// All children must finish PrepareFinalize before any child enters Finalize.
				AddDependency(*result, child1_entry->second.finish, child2_entry->second.prepare_finish,
				              PipelineDependencyType::OPTIONAL_DEPENDENCY);
			}
		}
	}

	return result;
}

} // namespace duckdb
