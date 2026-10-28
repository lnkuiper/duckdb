#include "catch.hpp"
#include "test_helpers.hpp"
#include "duckdb/common/local_file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/execution/operator/helper/physical_limit.hpp"
#include "duckdb/execution/operator/set/physical_cte.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/parallel/pipeline_broadcast_exchange.hpp"
#include "duckdb/parallel/pipeline_event.hpp"
#include "duckdb/parallel/pipeline_executor.hpp"
#include "duckdb/storage/buffer_manager.hpp"

#include <chrono>
#include <condition_variable>
#include <thread>

namespace duckdb {

#ifndef DUCKDB_NO_THREADS
namespace {
class ManualPipelineEvent : public PipelineEvent {
public:
	explicit ManualPipelineEvent(shared_ptr<Pipeline> pipeline) : PipelineEvent(std::move(pipeline)) {
		total_tasks = 2;
	}
};

class BroadcastReadSignal : public InterruptDoneSignalState {
public:
	bool IsSignalled() {
		annotated_lock_guard<annotated_mutex> guard(lock);
		return done;
	}
};

class GatedBroadcastFileSystem : public LocalFileSystem {
public:
	explicit GatedBroadcastFileSystem(string prefix_p) : prefix(std::move(prefix_p)) {
	}

	string GetName() const override {
		return "GatedBroadcastFileSystem";
	}

	bool CanHandleFile(const string &path) override {
		return StringUtil::StartsWith(path, prefix);
	}

	void Read(FileHandle &handle, void *buffer, int64_t nr_bytes, idx_t location) override {
		{
			unique_lock<mutex> guard(lock);
			if (armed) {
				entered = true;
				cv.notify_all();
				cv.wait(guard, [&] { return !armed; });
			}
		}
		LocalFileSystem::Read(handle, buffer, nr_bytes, location);
	}

	void Arm() {
		lock_guard<mutex> guard(lock);
		armed = true;
	}

	bool WaitForRead() {
		unique_lock<mutex> guard(lock);
		return cv.wait_for(guard, std::chrono::seconds(10), [&] { return entered; });
	}

	void Release() {
		lock_guard<mutex> guard(lock);
		armed = false;
		cv.notify_all();
	}

private:
	string prefix;
	mutex lock;
	std::condition_variable cv;
	bool armed = false;
	bool entered = false;
};
} // namespace

TEST_CASE("A completed LIMIT worker does not unregister another broadcast reader", "[broadcast_spool]") {
	DuckDB db;
	Connection con(db);
	auto &context = *con.context;
	REQUIRE_NO_FAIL(con.Query("SET threads=2"));
	auto pending = con.Submit("SELECT 1");
	REQUIRE_NO_FAIL(*pending);

	vector<LogicalType> types {LogicalType::BIGINT};
	auto exchange = make_shared_ptr<PipelineBroadcastExchange>(
	    context, types, PipelineBroadcastExchangeCompletionMode::RUN_TO_COMPLETION,
	    OrderPreservationType::INSERTION_ORDER, true, PipelineBroadcastExchangeBufferMode::BUFFER_ALL);
	auto consumer = exchange->RegisterConsumer();
	exchange->SelectBufferedConsumer(consumer, PipelineBroadcastExchangeScanMode::BATCH);
	exchange->Reset();
	auto producer = exchange->GetLocalState(context, 0);
	SourcePartitionInfo partition;
	InterruptState interrupt;
	DataChunk chunk;
	chunk.Initialize(context, types);
	for (idx_t i = 1; i <= 5; i++) {
		partition.batch_index = i <= 2 ? 1 : i - 1;
		partition.min_batch_index = partition.batch_index;
		if (i != 2) {
			REQUIRE(exchange->NextBatch(*producer, partition, interrupt) == SinkNextBatchType::READY);
		}
		chunk.SetChildCardinality(1);
		chunk.data[0].SetValue(0, Value::BIGINT(NumericCast<int64_t>(i)));
		REQUIRE(exchange->Push(chunk, *producer, partition, interrupt) == SinkResultType::NEED_MORE_INPUT);
	}
	exchange->Finish();

	// Keep an earlier batch assigned while a later worker reaches its local limit.
	auto reader = exchange->GetScanState();
	optional_idx batch;
	SourceBatchIndexState batch_state;
	chunk.Reset();
	REQUIRE(exchange->Scan(consumer, chunk, *reader, batch, batch_state, interrupt) ==
	        SourceResultType::HAVE_MORE_OUTPUT);
	REQUIRE(chunk.GetValue(0, 0) == Value::BIGINT(1));

	PhysicalPlan plan(Allocator::Get(context));
	auto &source = plan.Make<PhysicalCTEConsumerSource>(types, 5, TableIndex(0), exchange, consumer);
	auto &limit = plan.Make<PhysicalLimit>(types, BoundLimitNode::ConstantValue(1), BoundLimitNode(), 1);
	auto pipeline = make_shared_ptr<Pipeline>(context.GetExecutor());
	PipelineBuildState build_state;
	build_state.SetPipelineSource(*pipeline, source);
	build_state.SetPipelineSink(*pipeline, &limit, 0);
	pipeline->Ready();
	pipeline->Reset();
	auto event = make_shared_ptr<ManualPipelineEvent>(pipeline);
	auto first = make_shared_ptr<PipelineTask>(*pipeline, event);
	auto second = make_shared_ptr<PipelineTask>(*pipeline, event);
	REQUIRE(first->ExecuteTask(TaskExecutionMode::PROCESS_ALL) == TaskExecutionResult::TASK_FINISHED);
	REQUIRE(!event->IsFinished());

	chunk.Reset();
	REQUIRE(exchange->Scan(consumer, chunk, *reader, batch, batch_state, interrupt) ==
	        SourceResultType::HAVE_MORE_OUTPUT);
	REQUIRE(chunk.GetValue(0, 0) == Value::BIGINT(2));
	exchange->FinishScan(consumer, *reader);
	REQUIRE(second->ExecuteTask(TaskExecutionMode::PROCESS_ALL) == TaskExecutionResult::TASK_FINISHED);
	REQUIRE(event->IsFinished());

	reader = exchange->GetScanState();
	chunk.Reset();
	REQUIRE(exchange->Scan(consumer, chunk, *reader, batch, batch_state, interrupt) == SourceResultType::FINISHED);
}

TEST_CASE("Local broadcast completion wakes batch claimers while the producer is active", "[broadcast_spool]") {
	DuckDB db;
	Connection con(db);
	vector<LogicalType> types {LogicalType::BIGINT};
	PipelineBroadcastExchange exchange(*con.context, types,
	                                   PipelineBroadcastExchangeCompletionMode::STOP_WHEN_UNCONSUMED,
	                                   OrderPreservationType::INSERTION_ORDER, true);
	auto consumer = exchange.RegisterConsumer();
	exchange.SelectBufferedConsumer(consumer, PipelineBroadcastExchangeScanMode::BATCH);
	exchange.Reset();
	auto waiting = exchange.GetScanState();
	auto finishing = exchange.GetScanState();
	auto signal = make_shared_ptr<BroadcastReadSignal>();
	InterruptState interrupt {weak_ptr<InterruptDoneSignalState>(signal)};
	DataChunk chunk;
	chunk.Initialize(*con.context, types);
	optional_idx batch;
	SourceBatchIndexState batch_state;
	REQUIRE(exchange.Scan(consumer, chunk, *waiting, batch, batch_state, interrupt) == SourceResultType::BLOCKED);
	REQUIRE(!signal->IsSignalled());
	exchange.FinishScan(consumer, *finishing);
	REQUIRE(signal->IsSignalled());
	REQUIRE(exchange.Scan(consumer, chunk, *waiting, batch, batch_state, interrupt) == SourceResultType::FINISHED);
}

TEST_CASE("Broadcast spool reads survive local completion and are discarded on cancellation", "[broadcast_spool]") {
	const bool cancel = GENERATE(false, true);
	DuckDB db;
	Connection con(db);
	REQUIRE_NO_FAIL(con.Query("SET threads=1"));
	REQUIRE_NO_FAIL(con.Query("SET memory_limit='32MB'"));
	auto path = TestCreatePath("broadcast_completion_spill");
	auto fs = make_uniq<GatedBroadcastFileSystem>(path);
	auto &gate = *fs;
	FileSystem::GetFileSystem(*con.context).RegisterSubSystem(std::move(fs));
	REQUIRE_NO_FAIL(con.Query("SET temp_directory=" + Value(path).ToSQLString()));
	auto source = con.Query("SELECT repeat('abcdefghijklmnopqrstuvwxyz', 16) || i::VARCHAR "
	                        "FROM range(997) t(i)");
	REQUIRE_NO_FAIL(*source);
	auto input = source->Fetch();
	REQUIRE(input);
	auto types = input->GetTypes();
	PipelineBroadcastExchange exchange(*con.context, types, PipelineBroadcastExchangeCompletionMode::RUN_TO_COMPLETION,
	                                   OrderPreservationType::INSERTION_ORDER, true,
	                                   PipelineBroadcastExchangeBufferMode::BUFFER_ALL);
	auto consumer = exchange.RegisterConsumer();
	exchange.SelectBufferedConsumer(consumer, PipelineBroadcastExchangeScanMode::BATCH);
	exchange.Reset();
	auto local = exchange.GetLocalState(*con.context, 0);
	SourcePartitionInfo partition;
	InterruptState interrupt;
	const idx_t chunk_count = (128000 + input->size() - 1) / input->size();
	for (idx_t i = 0; i < chunk_count; i++) {
		partition.batch_index = i + 1;
		partition.min_batch_index = i + 1;
		REQUIRE(exchange.NextBatch(*local, partition, interrupt) == SinkNextBatchType::READY);
		REQUIRE(exchange.Push(*input, *local, partition, interrupt) == SinkResultType::NEED_MORE_INPUT);
	}
	REQUIRE(BufferManager::GetBufferManager(*con.context).GetUsedSwap() > 0);
	exchange.Finish();

	auto reader = exchange.GetScanState();
	auto finishing = exchange.GetScanState();
	DataChunk output;
	output.Initialize(*con.context, types);
	atomic<bool> stopped {false};
	string error;
	gate.Arm();
	std::thread task([&] {
		try {
			optional_idx batch;
			SourceBatchIndexState batch_state;
			do {
				if (exchange.Scan(consumer, output, *reader, batch, batch_state, interrupt) ==
				    SourceResultType::FINISHED) {
					break;
				}
			} while (!stopped);
		} catch (std::exception &ex) {
			error = ex.what();
		}
	});
	const bool read_reserved = gate.WaitForRead();
	stopped = true;
	if (cancel) {
		exchange.Cancel();
	} else {
		exchange.FinishScan(consumer, *finishing);
	}
	gate.Release();
	task.join();
	INFO(error);
	REQUIRE(read_reserved);
	REQUIRE(error.empty());
	if (cancel) {
		REQUIRE(output.size() == 0);
	} else {
		REQUIRE(output.size() == input->size());
		REQUIRE(output.GetValue(0, 0) == input->GetValue(0, 0));
		REQUIRE(output.GetValue(0, output.size() - 1) == input->GetValue(0, input->size() - 1));
		exchange.FinishScan(consumer, *reader);
	}
	exchange.Cancel();
}

TEST_CASE("Completed broadcast spools support concurrent spilling scans", "[broadcast_spool]") {
	DuckDB db;
	Connection con(db);
	REQUIRE_NO_FAIL(con.Query("SET threads=1"));
	REQUIRE_NO_FAIL(con.Query("SET memory_limit='64MB'"));
	REQUIRE_NO_FAIL(con.Query("SET temp_directory=" + Value(TestCreatePath("broadcast_spool_spill")).ToSQLString()));
	auto source =
	    con.Query("SELECT i, CASE WHEN i % 7 != 0 THEN repeat('abcdefghijklmnopqrstuvwxyz', 16) || i::VARCHAR "
	              "END s, [s, NULL, s] l FROM range(1, 998) t(i)");
	REQUIRE_NO_FAIL(*source);
	auto input = source->Fetch();
	REQUIRE(input);
	auto types = input->GetTypes();
	auto &allocator = BufferAllocator::Get(*con.context);
	PipelineBroadcastExchange exchange(*con.context, types, PipelineBroadcastExchangeCompletionMode::RUN_TO_COMPLETION,
	                                   OrderPreservationType::NO_ORDER, false,
	                                   PipelineBroadcastExchangeBufferMode::BUFFER_ALL);
	const idx_t reader_count = 4;
	const idx_t chunk_count = (64000 + input->size() - 1) / input->size();
	for (idx_t i = 0; i < reader_count; i++) {
		exchange.SelectBufferedConsumer(exchange.RegisterConsumer(), PipelineBroadcastExchangeScanMode::CHUNK);
	}

	// Reuse the exchange, retaining each previous scan's pins until its output has been checked.
	for (idx_t epoch = 0; epoch < 2; epoch++) {
		exchange.Reset();
		auto local = exchange.GetLocalState(*con.context, 0);
		SourcePartitionInfo partition;
		InterruptState interrupt;
		for (idx_t i = 0; i < chunk_count; i++) {
			REQUIRE(exchange.Push(*input, *local, partition, interrupt) == SinkResultType::NEED_MORE_INPUT);
		}
		REQUIRE(BufferManager::GetBufferManager(*con.context).GetUsedSwap() > 0);
		exchange.Finish();

		vector<shared_ptr<PipelineBroadcastExchangeScanState>> states;
		vector<unique_ptr<DataChunk>> outputs;
		vector<string> errors(reader_count);
		vector<idx_t> counts(reader_count, 0);
		vector<std::thread> readers;
		for (idx_t i = 0; i < reader_count; i++) {
			states.push_back(exchange.GetScanState());
			auto output = make_uniq<DataChunk>();
			output->Initialize(allocator, types);
			outputs.push_back(std::move(output));
		}
		for (idx_t i = 0; i < reader_count; i++) {
			readers.emplace_back([&, i] {
				try {
					auto &output = *outputs[i];
					optional_idx batch;
					SourceBatchIndexState batch_state;
					while (exchange.Scan(i, output, *states[i], batch, batch_state, interrupt) !=
					       SourceResultType::FINISHED) {
						for (idx_t row = 0; row < output.size(); row++) {
							auto input_row = counts[i]++ % input->size();
							for (idx_t col = 0; col < types.size(); col++) {
								if (!Value::NotDistinctFrom(output.GetValue(col, row),
								                            input->GetValue(col, input_row))) {
									throw InternalException("Incorrect broadcast spool value");
								}
							}
						}
					}
				} catch (std::exception &ex) {
					errors[i] = ex.what();
				}
			});
		}
		for (auto &reader : readers) {
			reader.join();
		}
		for (idx_t i = 0; i < reader_count; i++) {
			INFO(errors[i]);
			REQUIRE(errors[i].empty());
			REQUIRE(counts[i] == chunk_count * input->size());
		}
		// All consumers have retired: returned vectors still borrow their local reader's pins.
		exchange.Cancel();
		for (auto &output : outputs) {
			REQUIRE(output->size() == input->size());
			for (idx_t col = 0; col < types.size(); col++) {
				REQUIRE(Value::NotDistinctFrom(output->GetValue(col, 0), input->GetValue(col, 0)));
				REQUIRE(Value::NotDistinctFrom(output->GetValue(col, input->size() - 1),
				                               input->GetValue(col, input->size() - 1)));
			}
		}
	}
}
#endif

} // namespace duckdb
