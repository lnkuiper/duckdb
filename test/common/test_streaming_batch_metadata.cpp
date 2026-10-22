#include "catch.hpp"
#include "test_helpers.hpp"
#include "duckdb/common/vector/vector_writer.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/execution/operator/filter/physical_filter.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/bound_reference_expression.hpp"
#include "duckdb/parallel/pipeline_executor.hpp"

using namespace duckdb;

namespace {

struct MetadataSourceGlobal : GlobalSourceState {
	explicit MetadataSourceGlobal(BatchMetadataInfo requested_p) : requested(requested_p) {
	}
	BatchMetadataInfo requested;
	atomic<idx_t> next_batch {0};
	idx_t MaxThreads() override {
		return 4;
	}
};

struct MetadataSourceLocal : LocalSourceState {
	idx_t row = 0;
	idx_t end = 0;
	idx_t first = 0;
	idx_t batch = 0;
};

class MetadataSource : public PhysicalOperator {
public:
	MetadataSource(PhysicalPlan &plan, bool supported_p, bool completion_p = false)
	    : PhysicalOperator(plan, PhysicalOperatorType::DUMMY_SCAN, {LogicalType::BIGINT}, 7), supported(supported_p),
	      completion(completion_p), batch_limit(completion_p ? 3 : 2) {
	}

	bool supported;
	bool completion;
	atomic<idx_t> batch_limit;

	BatchMetadataInfo SourceMetadata() const override {
		auto result = supported ? BatchMetadataInfo::RowPositions() : BatchMetadataInfo();
		if (completion) {
			result.Merge(BatchMetadataInfo::Completion());
		}
		return result;
	}
	bool IsSource() const override {
		return true;
	}
	bool ParallelSource() const override {
		return true;
	}
	duckdb::unique_ptr<GlobalSourceState> GetGlobalSourceState(ClientContext &,
	                                                           const OperatorPartitionInfo &info) const override {
		REQUIRE_FALSE(info.RequiresPartitioning());
		return make_uniq<MetadataSourceGlobal>(info.metadata);
	}
	duckdb::unique_ptr<LocalSourceState> GetLocalSourceState(ExecutionContext &, GlobalSourceState &) const override {
		return make_uniq<MetadataSourceLocal>();
	}
	SourceResultType GetDataInternal(ExecutionContext &, DataChunk &chunk, OperatorSourceInput &input) const override {
		auto &local = input.local_state.Cast<MetadataSourceLocal>();
		auto &global = input.global_state.Cast<MetadataSourceGlobal>();
		if (local.row == local.end) {
			auto batch = global.next_batch.fetch_add(1);
			if (batch >= batch_limit) {
				return SourceResultType::FINISHED;
			}
			local.batch = batch;
			local.first = batch == 2 ? 5 : batch * 5;
			local.row = local.first;
			if (completion && batch == 1) {
				REQUIRE(global.requested.completion);
				input.metadata.batch.batch_index = batch;
				input.metadata.batch.first_row_index = 5;
				input.metadata.batch.row_count = 0;
				input.metadata.offset_in_batch = 0;
				input.batch_index_state = SourceBatchIndexState::ADVANCED;
				local.end = local.row;
				return SourceResultType::BLOCKED;
			}
			local.end = MinValue<idx_t>(local.row + 5, 7);
		}
		auto batch = local.batch;
		auto first = local.first;
		auto end = MinValue<idx_t>(first + 5, 7);
		auto count = MinValue<idx_t>(2, end - local.row);
		if (global.requested.Any()) {
			input.metadata.batch.batch_index = batch;
			input.metadata.batch.row_count = end - first;
			input.metadata.batch.first_row_index = first;
			input.metadata.offset_in_batch = local.row - first;
		}
		{
			auto writer = FlatVector::Writer<int64_t>(chunk.data[0], count);
			for (idx_t i = 0; i < count; i++) {
				writer.WriteValue(NumericCast<int64_t>(local.row++));
			}
		}
		chunk.CheckCardinality(count);
		return SourceResultType::HAVE_MORE_OUTPUT;
	}
};

struct MetadataSplitState : OperatorState {
	idx_t offset = 0;
	bool SupportsReuse() const override {
		return true;
	}
	void Reset() override {
		offset = 0;
	}
};

struct MetadataCacheState : CachingOperatorState {
	bool block_flush = true;
	bool SupportsReuse() const override {
		return true;
	}
	void Reset() override {
		ResetCachingState();
		block_flush = true;
	}
};

class MetadataCache : public CachingPhysicalOperator {
public:
	explicit MetadataCache(PhysicalPlan &plan)
	    : CachingPhysicalOperator(plan, PhysicalOperatorType::INVALID, {LogicalType::BIGINT}, 7) {
	}
	bool ParallelOperator() const override {
		return true;
	}
	BatchMetadataInfo OutputMetadata(const BatchMetadataInfo &input) const override {
		return input;
	}
	duckdb::unique_ptr<OperatorState> GetOperatorState(ExecutionContext &) const override {
		return make_uniq<MetadataCacheState>();
	}
	OperatorBatchResultType FlushBatch(ExecutionContext &context, OperatorChunk &output,
	                                   OperatorBatchInput &input) const override {
		auto &state = input.local_state.Cast<MetadataCacheState>();
		if (state.block_flush) {
			state.block_flush = false;
			return OperatorBatchResultType::BLOCKED;
		}
		state.block_flush = true;
		return CachingPhysicalOperator::FlushBatch(context, output, input);
	}
	OperatorResultType ExecuteInternal(ExecutionContext &, OperatorChunk &input, OperatorChunk &output,
	                                   GlobalOperatorState &, OperatorState &) const override {
		output.Reference(input);
		return OperatorResultType::NEED_MORE_INPUT;
	}
};

class MetadataSplit : public PhysicalOperator {
public:
	explicit MetadataSplit(PhysicalPlan &plan)
	    : PhysicalOperator(plan, PhysicalOperatorType::INVALID, {LogicalType::BIGINT}, 7) {
	}
	bool ParallelOperator() const override {
		return true;
	}
	BatchMetadataInfo OutputMetadata(const BatchMetadataInfo &input) const override {
		return input;
	}
	duckdb::unique_ptr<OperatorState> GetOperatorState(ExecutionContext &) const override {
		return make_uniq<MetadataSplitState>();
	}
	OperatorResultType Execute(ExecutionContext &, OperatorChunk &input, OperatorChunk &output, GlobalOperatorState &,
	                           OperatorState &state_p) const override {
		auto &state = state_p.Cast<MetadataSplitState>();
		output.chunk.Slice(input.chunk, state.offset, state.offset + 1);
		output.metadata = input.metadata;
		if (output.metadata.offset_in_batch.IsValid()) {
			output.metadata.offset_in_batch = output.metadata.offset_in_batch.GetIndex() + state.offset;
		}
		if (++state.offset < input.chunk.size()) {
			return OperatorResultType::HAVE_MORE_OUTPUT;
		}
		state.offset = 0;
		return OperatorResultType::NEED_MORE_INPUT;
	}
};

//! Consumes metadata but deliberately makes no guarantee for its output.
class MetadataBarrier : public PhysicalOperator {
public:
	explicit MetadataBarrier(PhysicalPlan &plan)
	    : PhysicalOperator(plan, PhysicalOperatorType::INVALID, {LogicalType::BIGINT}, 7) {
	}
	bool ParallelOperator() const override {
		return true;
	}
	BatchMetadataInfo RequestedMetadata() const override {
		return BatchMetadataInfo::RowPositions();
	}
	OperatorResultType Execute(ExecutionContext &, OperatorChunk &input, OperatorChunk &output, GlobalOperatorState &,
	                           OperatorState &) const override {
		REQUIRE(input.metadata.offset_in_batch.IsValid());
		output.Reference(input);
		return OperatorResultType::NEED_MORE_INPUT;
	}
};

struct MetadataObservation : GlobalOperatorState {
	mutex lock;
	duckdb::vector<ChunkMetadata> metadata;
	duckdb::vector<idx_t> values;
};

class MetadataObserver : public PhysicalOperator {
public:
	explicit MetadataObserver(PhysicalPlan &plan)
	    : PhysicalOperator(plan, PhysicalOperatorType::INVALID, {LogicalType::BIGINT}, 7) {
	}
	bool ParallelOperator() const override {
		return true;
	}
	BatchMetadataInfo RequestedMetadata() const override {
		return BatchMetadataInfo::RowPositions();
	}
	BatchMetadataInfo OutputMetadata(const BatchMetadataInfo &input) const override {
		return input;
	}
	duckdb::unique_ptr<GlobalOperatorState> GetGlobalOperatorState(ClientContext &) const override {
		return make_uniq<MetadataObservation>();
	}
	OperatorResultType Execute(ExecutionContext &, OperatorChunk &input, OperatorChunk &output,
	                           GlobalOperatorState &global_p, OperatorState &) const override {
		auto &global = global_p.Cast<MetadataObservation>();
		duckdb::lock_guard<mutex> guard(global.lock);
		global.metadata.push_back(input.metadata);
		global.values.push_back(input.chunk.GetValue(0, 0).GetValue<idx_t>());
		output.Reference(input);
		return OperatorResultType::NEED_MORE_INPUT;
	}
};

struct MetadataSinkState : LocalSinkState {
	bool block = true;
	ChunkMetadata blocked_metadata;
};

class MetadataSink : public PhysicalOperator {
public:
	MetadataSink(PhysicalPlan &plan, bool require_metadata_p)
	    : PhysicalOperator(plan, PhysicalOperatorType::INVALID, {LogicalType::BIGINT}, 7),
	      require_metadata(require_metadata_p) {
	}
	bool require_metadata;

	OperatorPartitionInfo RequiredPartitionInfo() const override {
		auto result = OperatorPartitionInfo::NoPartitionInfo();
		if (require_metadata) {
			result.metadata = BatchMetadataInfo::RowPositions();
		}
		return result;
	}
	bool IsSink() const override {
		return true;
	}
	bool ParallelSink() const override {
		return true;
	}
	duckdb::unique_ptr<LocalSinkState> GetLocalSinkState(ExecutionContext &) const override {
		return make_uniq<MetadataSinkState>();
	}
	SinkResultType Sink(ExecutionContext &, DataChunk &chunk, OperatorSinkInput &input) const override {
		REQUIRE_FALSE(input.local_state.partition_info.batch_index.IsValid());
		auto &state = input.local_state.Cast<MetadataSinkState>();
		if (require_metadata) {
			REQUIRE(input.metadata);
			auto &metadata = *input.metadata;
			metadata.Verify(BatchMetadataInfo::RowPositions(), chunk.size());
			REQUIRE(metadata.batch.first_row_index.GetIndex() + metadata.offset_in_batch.GetIndex() ==
			        chunk.GetValue(0, 0).GetValue<idx_t>());
			if (state.block) {
				state.blocked_metadata = metadata;
			} else {
				auto &blocked = state.blocked_metadata;
				REQUIRE(metadata.batch.batch_index == blocked.batch.batch_index);
				REQUIRE(metadata.batch.stream_id == blocked.batch.stream_id);
				REQUIRE(metadata.batch.row_count == blocked.batch.row_count);
				REQUIRE(metadata.batch.first_row_index == blocked.batch.first_row_index);
				REQUIRE(metadata.offset_in_batch == blocked.offset_in_batch);
			}
		}
		state.block = !state.block;
		return state.block ? SinkResultType::NEED_MORE_INPUT : SinkResultType::BLOCKED;
	}
};

} // namespace

TEST_CASE("Streaming operators request chunk metadata independently of sink ordering", "[streaming_batch_metadata]") {
	bool supported = true;
	bool barrier = false;
	bool require_sink_metadata = false;
	bool cache = false;
	SECTION("Split chunks retain positions through blocked sinks and executor reset") {
	}
	SECTION("Caches retain separate batch descriptors without completion requests") {
		cache = true;
	}
	SECTION("Unsupported sources leave optional metadata unavailable") {
		supported = false;
	}
	SECTION("An operator's default output contract invalidates even copied metadata") {
		barrier = true;
	}
	SECTION("Mandatory sink metadata survives blocked retries") {
		require_sink_metadata = true;
	}

	DuckDB db;
	Connection con(db);
#ifndef DUCKDB_NO_THREADS
	REQUIRE_NO_FAIL(con.Query("SET threads=4"));
#endif
	auto pending = con.Submit("SELECT 1");
	REQUIRE_NO_FAIL(*pending);
	auto &context = *con.context;
	PhysicalPlan plan(Allocator::Get(context));
	auto &source = plan.Make<MetadataSource>(supported);
	auto &split = plan.Make<MetadataSplit>();
	auto &observer = plan.Make<MetadataObserver>();
	auto &sink = plan.Make<MetadataSink>(require_sink_metadata);
	auto pipeline = make_shared_ptr<Pipeline>(context.GetExecutor());
	PipelineBuildState state;
	state.SetPipelineSource(*pipeline, source);
	state.SetPipelineSink(*pipeline, &sink, 0);
	// Ready reverses the construction order.
	state.AddPipelineOperator(*pipeline, observer);
	if (barrier) {
		state.AddPipelineOperator(*pipeline, plan.Make<MetadataBarrier>());
	}
	state.AddPipelineOperator(*pipeline, split);
	if (cache) {
		state.AddPipelineOperator(*pipeline, plan.Make<MetadataCache>());
	}
	pipeline->Ready();
	REQUIRE(pipeline->GetMetadataPlan().source_partition_info.metadata.row_positions == supported);
	REQUIRE_FALSE(pipeline->GetMetadataPlan().sink_partition_info.RequiresPartitioning());
	pipeline->Reset();
#ifndef DUCKDB_NO_THREADS
	REQUIRE(pipeline->GetMaxThreads() == 4);
#endif
	PipelineExecutor executor(context, *pipeline);
	const duckdb::vector<idx_t> expected {0, 1, 2, 3, 4, 5, 6};
	for (idx_t run = 0; run < 2; run++) {
		if (run > 0) {
			pipeline->ResetForReschedule(true);
		}
		executor.PrepareForExecution();
		idx_t calls = 0;
		while (executor.Execute(1) != PipelineExecuteResult::FINISHED) {
			REQUIRE(++calls < 100);
		}
		auto &observed = observer.op_state->Cast<MetadataObservation>();
		REQUIRE(observed.values == expected);
		REQUIRE(observed.metadata.size() == 7);
		for (idx_t row = 0; row < 7; row++) {
			auto &metadata = observed.metadata[row];
			if (supported && !barrier) {
				REQUIRE(metadata.batch.stream_id.GetIndex() == pipeline->GetStreamId());
				REQUIRE(metadata.batch.batch_index.GetIndex() == row / 5 + 1);
				REQUIRE(metadata.batch.row_count.GetIndex() == (row < 5 ? 5 : 2));
				REQUIRE(metadata.batch.first_row_index.GetIndex() + metadata.offset_in_batch.GetIndex() == row);
			} else {
				REQUIRE_FALSE(metadata.batch.batch_index.IsValid());
				REQUIRE_FALSE(metadata.offset_in_batch.IsValid());
			}
		}
	}
}

namespace {

struct BatchCompletionResult : GlobalSinkState {
	duckdb::vector<idx_t> rows;
	duckdb::vector<BatchMetadata> completed;
	duckdb::vector<idx_t> admitted;
};

struct BatchCompletionState : LocalSinkState {
	bool block_rows = true;
	bool block_completion = true;
	bool block_next = true;
	BatchMetadata pending_completion;
};

class BatchCompletionSink : public PhysicalOperator {
public:
	explicit BatchCompletionSink(PhysicalPlan &plan)
	    : PhysicalOperator(plan, PhysicalOperatorType::INVALID, {LogicalType::BIGINT}, 7) {
	}
	bool IsSink() const override {
		return true;
	}
	bool ParallelSink() const override {
		return true;
	}
	OperatorPartitionInfo RequiredPartitionInfo() const override {
		auto result = OperatorPartitionInfo::NoPartitionInfo();
		result.metadata = BatchMetadataInfo::Completion();
		return result;
	}
	BatchMetadataInfo RequestedMetadata() const override {
		return BatchMetadataInfo::RowPositions();
	}
	duckdb::unique_ptr<GlobalSinkState> GetGlobalSinkState(ClientContext &) const override {
		return make_uniq<BatchCompletionResult>();
	}
	duckdb::unique_ptr<LocalSinkState> GetLocalSinkState(ExecutionContext &) const override {
		return make_uniq<BatchCompletionState>();
	}
	SinkResultType Sink(ExecutionContext &, DataChunk &chunk, OperatorSinkInput &input) const override {
		auto &local = input.local_state.Cast<BatchCompletionState>();
		auto &global = input.global_state.Cast<BatchCompletionResult>();
		REQUIRE(input.metadata);
		auto &metadata = *input.metadata;
		REQUIRE(metadata.batch.batch_index == local.partition_info.batch_index);
		for (auto &batch : global.completed) {
			REQUIRE(batch.batch_index != metadata.batch.batch_index);
		}
		if (metadata.offset_in_batch.IsValid()) {
			REQUIRE(metadata.batch.first_row_index.GetIndex() + metadata.offset_in_batch.GetIndex() ==
			        chunk.GetValue(0, 0).GetValue<idx_t>());
		}
		if (local.block_rows) {
			local.block_rows = false;
			return SinkResultType::BLOCKED;
		}
		local.block_rows = true;
		for (idx_t row = 0; row < chunk.size(); row++) {
			global.rows.push_back(chunk.GetValue(0, row).GetValue<idx_t>());
		}
		return SinkResultType::NEED_MORE_INPUT;
	}
	SinkNextBatchType CompleteBatch(ExecutionContext &, OperatorSinkBatchInput &input) const override {
		auto &local = input.local_state.Cast<BatchCompletionState>();
		auto &global = input.global_state.Cast<BatchCompletionResult>();
		REQUIRE(input.batch.batch_index == local.partition_info.batch_index);
		REQUIRE(local.partition_info.min_batch_index.GetIndex() <= input.batch.batch_index.GetIndex());
		REQUIRE(input.batch.batch_index.GetIndex() == global.completed.size() + 1);
		if (local.block_completion) {
			local.pending_completion = input.batch;
			local.block_completion = false;
			return SinkNextBatchType::BLOCKED;
		}
		REQUIRE(input.batch.batch_index == local.pending_completion.batch_index);
		REQUIRE(input.batch.row_count == local.pending_completion.row_count);
		REQUIRE(input.batch.first_row_index == local.pending_completion.first_row_index);
		local.block_completion = true;
		global.completed.push_back(input.batch);
		return SinkNextBatchType::READY;
	}
	SinkNextBatchType NextBatch(ExecutionContext &, OperatorSinkNextBatchInput &input) const override {
		auto &local = input.local_state.Cast<BatchCompletionState>();
		if (local.block_next) {
			local.block_next = false;
			return SinkNextBatchType::BLOCKED;
		}
		local.block_next = true;
		input.global_state.Cast<BatchCompletionResult>().admitted.push_back(
		    local.partition_info.batch_index.GetIndex());
		return SinkNextBatchType::READY;
	}
};

} // namespace

TEST_CASE("Batch completion drains caches before acknowledgement and minimum release", "[streaming_batch_metadata]") {
	bool filter = false;
	bool cache_enabled = true;
	SECTION("Dense batches, an explicit empty advance, split cached output and blocked callbacks") {
	}
	SECTION("Entirely filtered batches still complete with unknown cardinality") {
		filter = true;
	}
	SECTION("Completion is independent of cache allocation") {
		cache_enabled = false;
	}

	DuckDB db;
	Connection con(db);
	if (!cache_enabled) {
		REQUIRE_NO_FAIL(con.Query("SET enable_caching_operators=false"));
	}
	auto pending = con.Submit("SELECT 1");
	REQUIRE_NO_FAIL(*pending);
	auto &context = *con.context;
	PhysicalPlan plan(Allocator::Get(context));
	auto &source = plan.Make<MetadataSource>(true, true);
	auto &sink = plan.Make<BatchCompletionSink>();
	auto pipeline = make_shared_ptr<Pipeline>(context.GetExecutor());
	PipelineBuildState state;
	state.SetPipelineSource(*pipeline, source);
	state.SetPipelineSink(*pipeline, &sink, 0);
	state.AddPipelineOperator(*pipeline, plan.Make<MetadataSplit>());
	if (filter) {
		duckdb::vector<duckdb::unique_ptr<Expression>> predicates;
		predicates.push_back(BoundComparisonExpression::Create(
		    ExpressionType::COMPARE_GREATERTHANOREQUALTO, make_uniq<BoundReferenceExpression>(LogicalType::BIGINT, 0),
		    make_uniq<BoundConstantExpression>(Value::BIGINT(5))));
		state.AddPipelineOperator(
		    *pipeline,
		    plan.Make<PhysicalFilter>(duckdb::vector<LogicalType> {LogicalType::BIGINT}, std::move(predicates), 7));
	}
	state.AddPipelineOperator(*pipeline, plan.Make<MetadataCache>());
	pipeline->Ready();
	REQUIRE(pipeline->GetMetadataPlan().source_partition_info.metadata.completion);
	REQUIRE_FALSE(pipeline->GetMetadataPlan().source_partition_info.RequiresPartitioning());
	pipeline->Reset();
	PipelineExecutor executor(context, *pipeline);
	for (idx_t run = 0; run < 2; run++) {
		if (run > 0) {
			pipeline->ResetForReschedule(true);
		}
		executor.PrepareForExecution();
		idx_t calls = 0;
		while (executor.Execute(1) != PipelineExecuteResult::FINISHED) {
			REQUIRE(++calls < 200);
		}
		auto &result = sink.sink_state->Cast<BatchCompletionResult>();
		auto expected = filter ? duckdb::vector<idx_t> {5, 6} : duckdb::vector<idx_t> {0, 1, 2, 3, 4, 5, 6};
		REQUIRE(result.rows == expected);
		REQUIRE(result.completed.size() == 3);
		for (idx_t batch = 0; batch < 3; batch++) {
			auto &metadata = result.completed[batch];
			REQUIRE(metadata.batch_index.GetIndex() == batch + 1);
			REQUIRE(metadata.stream_id.GetIndex() == pipeline->GetStreamId());
			if (filter) {
				REQUIRE_FALSE(metadata.row_count.IsValid());
				REQUIRE_FALSE(metadata.first_row_index.IsValid());
			} else {
				REQUIRE(metadata.row_count.GetIndex() == (batch == 0 ? 5 : batch == 1 ? 0 : 2));
				REQUIRE(metadata.first_row_index.GetIndex() == (batch == 0 ? 0 : 5));
			}
		}
	}
}

TEST_CASE("Completion-only workers release their reservation even without an input batch",
          "[streaming_batch_metadata]") {
	DuckDB db;
	Connection con(db);
	auto pending = con.Submit("SELECT 1");
	REQUIRE_NO_FAIL(*pending);
	auto &context = *con.context;
	PhysicalPlan plan(Allocator::Get(context));
	auto &source = plan.Make<MetadataSource>(true, true);
	auto &sink = plan.Make<BatchCompletionSink>();
	auto pipeline = make_shared_ptr<Pipeline>(context.GetExecutor());
	PipelineBuildState state;
	state.SetPipelineSource(*pipeline, source);
	state.SetPipelineSink(*pipeline, &sink, 0);
	state.AddPipelineOperator(*pipeline, plan.Make<MetadataCache>());
	pipeline->Ready();
	pipeline->Reset();
	PipelineExecutor active(context, *pipeline);
	PipelineExecutor idle(context, *pipeline);
	auto &result = sink.sink_state->Cast<BatchCompletionResult>();
	idx_t calls = 0;
	while (result.admitted.empty()) {
		REQUIRE(active.Execute(1) != PipelineExecuteResult::FINISHED);
		REQUIRE(++calls < 20);
	}
	REQUIRE(result.admitted[0] == 1);
	// No further batches are available, but the first worker still owns the first one.
	source.Cast<MetadataSource>().batch_limit = 1;
	while (idle.Execute(1) != PipelineExecuteResult::FINISHED) {
		REQUIRE(++calls < 50);
	}
	REQUIRE(result.completed.empty());
	REQUIRE(pipeline->UpdateBatchIndex(1, 1) == 1);
	while (active.Execute(1) != PipelineExecuteResult::FINISHED) {
		REQUIRE(++calls < 100);
	}
	const duckdb::vector<idx_t> expected {0, 1, 2, 3, 4};
	REQUIRE(result.rows == expected);
	REQUIRE(result.completed.size() == 1);
	auto terminal = PipelineBuildState::BATCH_INCREMENT - 1;
	REQUIRE(pipeline->UpdateBatchIndex(terminal, terminal) == terminal);
}
