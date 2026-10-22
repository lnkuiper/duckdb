#include "catch.hpp"
#include "test_helpers.hpp"
#include "duckdb/common/vector/vector_writer.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
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
};

class MetadataSource : public PhysicalOperator {
public:
	MetadataSource(PhysicalPlan &plan, bool supported_p)
	    : PhysicalOperator(plan, PhysicalOperatorType::DUMMY_SCAN, {LogicalType::BIGINT}, 7), supported(supported_p) {
	}

	bool supported;

	BatchMetadataInfo SourceMetadata() const override {
		return supported ? BatchMetadataInfo::RowPositions() : BatchMetadataInfo();
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
			if (batch >= 2) {
				return SourceResultType::FINISHED;
			}
			local.row = batch * 5;
			local.end = MinValue<idx_t>(local.row + 5, 7);
		}
		auto batch = local.row / 5;
		auto first = batch * 5;
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
	SECTION("Split chunks retain positions through blocked sinks and executor reset") {
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
