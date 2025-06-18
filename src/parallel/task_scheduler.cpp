#include "duckdb/parallel/task_scheduler.hpp"

#include "duckdb/common/chrono.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"

#ifndef DUCKDB_NO_THREADS
#include "concurrentqueue.h"
#include "duckdb/common/thread.hpp"
#include "lightweightsemaphore.h"

#include <thread>
#else
#include <queue>
#endif

#if defined(_WIN32)
#include <windows.h>
#elif defined(__GNUC__)
#include <sched.h>
#include <unistd.h>
#endif

namespace duckdb {

struct SchedulerThread {
#ifndef DUCKDB_NO_THREADS
	explicit SchedulerThread(unique_ptr<thread> thread_p) : internal_thread(std::move(thread_p)) {
	}

	unique_ptr<thread> internal_thread;
#endif
};

#ifndef DUCKDB_NO_THREADS
typedef duckdb_moodycamel::ProducerToken producer_token_t;
typedef duckdb_moodycamel::ConcurrentQueue<shared_ptr<Task>> concurrent_queue_t;
typedef duckdb_moodycamel::LightweightSemaphore lightweight_semaphore_t;

struct PartitionedQueue {
	PartitionedQueue() : queue(), semaphore(), is_stealing(false) {
	}
	concurrent_queue_t queue;
	lightweight_semaphore_t semaphore;
	atomic<bool> is_stealing;
};

struct ConcurrentQueue {
private:
	typedef std::make_signed<std::size_t>::type ssize_t;

public:
	producer_token_t CreateProducerToken();

	void SetThreads(const idx_t &num_threads);

	void Enqueue(ProducerToken &token, shared_ptr<Task> task);
	void EnqueueBulk(ProducerToken &token, vector<shared_ptr<Task>> &tasks);
	bool DequeueFromProducer(ProducerToken &token, shared_ptr<Task> &task);
	bool TryDequeue(shared_ptr<Task> &task);

	idx_t GetNumberOfTasks() const;
	idx_t GetProducerCount() const;
	idx_t GetTaskCountForProducer(ProducerToken &token) const;

	bool Wait();
	bool Wait(int64_t timeout_usecs);
	void Signal(const idx_t &n = 1);

	bool WaitForTask(const bool &allocator_background_threads, const idx_t &allocator_flush_threshold,
	                 shared_ptr<Task> &task, const optional_idx &thread_idx);

private:
	atomic<idx_t> num_threads;

	PartitionedQueue global;

	//! To reduce contention, threads do single access within their group,
	//! and the partition only does bulk access with the globals
	static constexpr idx_t PARTITIONS = 16;
	static constexpr idx_t THREADS_PER_PARTITION = 8;

	array<PartitionedQueue, PARTITIONS> partitions;
};

struct QueueProducerToken {
	explicit QueueProducerToken(ConcurrentQueue &queue) : queue_token(queue.CreateProducerToken()) {
	}

	producer_token_t queue_token;
};

producer_token_t ConcurrentQueue::CreateProducerToken() {
	return producer_token_t(global.queue);
}

void ConcurrentQueue::SetThreads(const idx_t &num_threads_p) {
	num_threads = num_threads_p;

	// Move group tasks to the global queue again
	// Note that we have lost producer information
	// This is unfortunate, but only happens when the number of threads changes,
	// which should be infrequent WHILE queries are being executed
	vector<shared_ptr<Task>> tasks(THREADS_PER_PARTITION);
	for (auto &partition : partitions) {
		const auto num_tasks = partition.queue.try_dequeue_bulk(tasks.begin(), THREADS_PER_PARTITION);
		global.queue.enqueue_bulk(std::make_move_iterator(tasks.begin()), num_tasks);
	}
}

void ConcurrentQueue::Enqueue(ProducerToken &token, shared_ptr<Task> task) {
	vector<shared_ptr<Task>> tasks {std::move(task)};
	EnqueueBulk(token, tasks);
}

void ConcurrentQueue::EnqueueBulk(ProducerToken &token, vector<shared_ptr<Task>> &tasks) {
	lock_guard<mutex> producer_lock(token.producer_lock);
	for (const auto &task : tasks) {
		task->token = token;
	}
	if (global.queue.enqueue_bulk(token.token->queue_token, std::make_move_iterator(tasks.begin()), tasks.size())) {
		Signal(tasks.size());
	} else {
		throw InternalException("Could not schedule task!");
	}
}

bool ConcurrentQueue::DequeueFromProducer(ProducerToken &token, shared_ptr<Task> &task) {
	lock_guard<mutex> producer_lock(token.producer_lock);
	return global.queue.try_dequeue_from_producer(token.token->queue_token, task);
}

bool ConcurrentQueue::TryDequeue(shared_ptr<Task> &task) {
	return global.queue.try_dequeue(task);
}

bool ConcurrentQueue::WaitForTask(const bool &allocator_background_threads, const idx_t &allocator_flush_threshold,
                                  shared_ptr<Task> &task, const optional_idx &thread_idx) {
	// Initial wait time of 0.5s (in mu s) before flushing
	static constexpr int64_t INITIAL_WAIT = 500000;

	const auto num_threads_local = num_threads.load();
	auto &partition =
	    thread_idx.IsValid() ? partitions[(thread_idx.GetIndex() / THREADS_PER_PARTITION) % PARTITIONS] : global;

	// Wait for the semaphore, doing some memory management while waiting
	if (!Allocator::SupportsFlush()) {
		// allocator can't flush, start an untimed wait
		partition.semaphore.wait();
	} else if (!partition.semaphore.wait(INITIAL_WAIT)) {
		// allocator can flush, we flush this threads outstanding allocations after it was idle for 0.5s
		Allocator::ThreadFlush(allocator_background_threads, allocator_flush_threshold, num_threads_local);
		const auto decay_delay = Allocator::DecayDelay();
		if (!decay_delay.IsValid()) {
			partition.semaphore.wait(); // no decay delay specified - untimed wait
		} else {
			const auto remaining_wait = UnsafeNumericCast<int64_t>(decay_delay.GetIndex()) * 1000000 - INITIAL_WAIT;
			if (!partition.semaphore.wait(remaining_wait)) {
				// in total, the thread was idle for the entire decay delay (note: seconds converted to mus)
				// mark it as idle and start an untimed wait
				Allocator::ThreadIdle();
				partition.semaphore.wait();
			}
		}
	}

	if (thread_idx.IsValid()) {
		// We have a thread index, try to get a task from the partition
		auto expected = false;
		if (partition.is_stealing.compare_exchange_strong(expected, true, std::memory_order_acquire)) {
			// Value was false before exchanging, this thread is the only thread in here
			bool got_task = partition.queue.try_dequeue(task);
			if (!got_task) {
				// No tasks left in the partition queue, steal up to THREADS_PER_PARTITION from global
				const auto num_threads_aligned = AlignValueFloor<idx_t, THREADS_PER_PARTITION>(num_threads_local);
				auto num_tasks = THREADS_PER_PARTITION;
				if (num_threads_local != num_threads_aligned && thread_idx.GetIndex() >= num_threads_aligned) {
					// Round down if we're the last partition
					num_tasks = thread_idx.GetIndex() - num_threads_aligned + 1;
				}

				// Try to get this number of tasks from the global semaphore/queue
				num_tasks = NumericCast<idx_t>(global.semaphore.tryWaitMany(NumericCast<ssize_t>(num_tasks)));
				vector<shared_ptr<Task>> tasks(num_tasks);
				num_tasks = global.queue.try_dequeue_bulk(tasks.begin(), num_tasks);

				// Get a task for this thread (if any)
				if (num_tasks != 0) {
					task = std::move(tasks[--num_tasks]);
					got_task = true;
				}

				// Bulk enqueue in the partition and mark that we're done stealing so other threads can grab a task
				partition.queue.enqueue_bulk(std::make_move_iterator(tasks.begin()), num_tasks);
			}

			// We are done stealing
			partition.is_stealing.store(false, std::memory_order_release);

			// If we got a task, immediately return
			if (got_task) {
				return true;
			}
		} else {
			// Another thread is currently stealing, wait
			while (partition.is_stealing.load(std::memory_order_acquire)) {
				TaskScheduler::YieldThread();
			}

			// Try to get one from this partition
			if (partition.queue.try_dequeue(task)) {
				return true;
			}
		}
	}

	// We got here if we had a thread index but did not get a task,
	// or if we did not have a thread index. In both cases we try to get a task from the global queue
	return global.queue.try_dequeue(task);
}

idx_t ConcurrentQueue::GetNumberOfTasks() const {
	return global.queue.size_approx();
}

idx_t ConcurrentQueue::GetProducerCount() const {
	return global.queue.size_producers_approx();
}

idx_t ConcurrentQueue::GetTaskCountForProducer(ProducerToken &token) const {
	lock_guard<mutex> producer_lock(token.producer_lock);
	return global.queue.size_producer_approx(token.token->queue_token);
}

bool ConcurrentQueue::Wait() {
	return global.semaphore.wait();
}

bool ConcurrentQueue::Wait(int64_t timeout_usecs) {
	return global.semaphore.wait(timeout_usecs);
}

void ConcurrentQueue::Signal(const idx_t &n) {
	const auto count = NumericCast<ssize_t>(n);

	// Signal the global semaphore with the appropriate amount
	global.semaphore.signal(count);

	// Give each group a single dequeue signal
	for (auto &partition : partitions) {
		partition.semaphore.signal(THREADS_PER_PARTITION);
	}
}

#else
struct ConcurrentQueue {
	reference_map_t<QueueProducerToken, std::queue<shared_ptr<Task>>> q;
	mutex qlock;

	void Enqueue(ProducerToken &token, shared_ptr<Task> task);
	void EnqueueBulk(ProducerToken &token, vector<shared_ptr<Task>> &tasks);
	bool DequeueFromProducer(ProducerToken &token, shared_ptr<Task> &task);
};

void ConcurrentQueue::Enqueue(ProducerToken &token, shared_ptr<Task> task) {
	lock_guard<mutex> lock(qlock);
	task->token = token;
	q[std::ref(*token.token)].push(std::move(task));
}

void ConcurrentQueue::EnqueueBulk(ProducerToken &token, vector<shared_ptr<Task>> &tasks) {
	lock_guard<mutex> lock(qlock);
	for (auto &task : tasks) {
		task->token = token;
		q[std::ref(*token.token)].push(std::move(task));
	}
}

bool ConcurrentQueue::DequeueFromProducer(ProducerToken &token, shared_ptr<Task> &task) {
	lock_guard<mutex> lock(qlock);
	D_ASSERT(!q.empty());

	const auto it = q.find(std::ref(*token.token));
	if (it == q.end() || it->second.empty()) {
		return false;
	}

	task = std::move(it->second.front());
	it->second.pop();

	return true;
}

struct QueueProducerToken {
	explicit QueueProducerToken(ConcurrentQueue &queue) : queue(&queue) {
	}

	~QueueProducerToken() {
		lock_guard<mutex> lock(queue->qlock);
		queue->q.erase(*this);
	}

private:
	ConcurrentQueue *queue;
};
#endif

ProducerToken::ProducerToken(TaskScheduler &scheduler, unique_ptr<QueueProducerToken> token)
    : scheduler(scheduler), token(std::move(token)) {
}

ProducerToken::~ProducerToken() {
}

TaskScheduler::TaskScheduler(DatabaseInstance &db)
    : db(db), queue(make_uniq<ConcurrentQueue>()),
      allocator_flush_threshold(db.config.options.allocator_flush_threshold),
      allocator_background_threads(db.config.options.allocator_background_threads), requested_thread_count(0),
      current_thread_count(1) {
	SetAllocatorBackgroundThreads(db.config.options.allocator_background_threads);
}

TaskScheduler::~TaskScheduler() {
#ifndef DUCKDB_NO_THREADS
	try {
		RelaunchThreadsInternal(0);
	} catch (...) {
		// nothing we can do in the destructor if this fails
	}
#endif
}

TaskScheduler &TaskScheduler::GetScheduler(ClientContext &context) {
	return TaskScheduler::GetScheduler(DatabaseInstance::GetDatabase(context));
}

TaskScheduler &TaskScheduler::GetScheduler(DatabaseInstance &db) {
	return db.GetScheduler();
}

unique_ptr<ProducerToken> TaskScheduler::CreateProducer() {
	auto token = make_uniq<QueueProducerToken>(*queue);
	return make_uniq<ProducerToken>(*this, std::move(token));
}

void TaskScheduler::ScheduleTask(ProducerToken &token, shared_ptr<Task> task) {
	// Enqueue a task for the given producer token and signal any sleeping threads
	queue->Enqueue(token, std::move(task));
}

void TaskScheduler::ScheduleManyTasks(ProducerToken &producer, vector<shared_ptr<Task>> &tasks) {
	queue->EnqueueBulk(producer, tasks);
}

bool TaskScheduler::GetTaskFromProducer(ProducerToken &token, shared_ptr<Task> &task) {
	return queue->DequeueFromProducer(token, task);
}

void TaskScheduler::ExecuteForever(atomic<bool> *marker, const optional_idx &thread_idx) {
#ifndef DUCKDB_NO_THREADS
	auto &config = DBConfig::GetConfig(db);
	shared_ptr<Task> task;
	// loop until the marker is set to false
	while (*marker) {
		if (queue->WaitForTask(allocator_background_threads, allocator_flush_threshold, task, thread_idx)) {
			auto process_mode = config.options.scheduler_process_partial ? TaskExecutionMode::PROCESS_PARTIAL
			                                                             : TaskExecutionMode::PROCESS_ALL;
			auto execute_result = task->Execute(process_mode);

			switch (execute_result) {
			case TaskExecutionResult::TASK_FINISHED:
			case TaskExecutionResult::TASK_ERROR:
				task.reset();
				break;
			case TaskExecutionResult::TASK_NOT_FINISHED: {
				// task is not finished - reschedule immediately
				auto &token = *task->token;
				queue->Enqueue(token, std::move(task));
				break;
			}
			case TaskExecutionResult::TASK_BLOCKED:
				task->Deschedule();
				task.reset();
				break;
			}
		}
	}
	// this thread will exit, flush all of its outstanding allocations
	if (Allocator::SupportsFlush()) {
		Allocator::ThreadFlush(allocator_background_threads, 0, NumericCast<idx_t>(requested_thread_count.load()));
		Allocator::ThreadIdle();
	}
#else
	throw NotImplementedException("DuckDB was compiled without threads! Background thread loop is not allowed.");
#endif
}

idx_t TaskScheduler::ExecuteTasks(atomic<bool> *marker, idx_t max_tasks) {
#ifndef DUCKDB_NO_THREADS
	idx_t completed_tasks = 0;
	// loop until the marker is set to false
	while (*marker && completed_tasks < max_tasks) {
		shared_ptr<Task> task;
		if (!queue->TryDequeue(task)) {
			return completed_tasks;
		}
		auto execute_result = task->Execute(TaskExecutionMode::PROCESS_ALL);

		switch (execute_result) {
		case TaskExecutionResult::TASK_FINISHED:
		case TaskExecutionResult::TASK_ERROR:
			task.reset();
			completed_tasks++;
			break;
		case TaskExecutionResult::TASK_NOT_FINISHED:
			throw InternalException("Task should not return TASK_NOT_FINISHED in PROCESS_ALL mode");
		case TaskExecutionResult::TASK_BLOCKED:
			task->Deschedule();
			task.reset();
			break;
		}
	}
	return completed_tasks;
#else
	throw NotImplementedException("DuckDB was compiled without threads! Background thread loop is not allowed.");
#endif
}

void TaskScheduler::ExecuteTasks(idx_t max_tasks) {
#ifndef DUCKDB_NO_THREADS
	shared_ptr<Task> task;
	for (idx_t i = 0; i < max_tasks; i++) {
		queue->Wait(TASK_TIMEOUT_USECS);
		if (!queue->TryDequeue(task)) {
			return;
		}
		try {
			auto execute_result = task->Execute(TaskExecutionMode::PROCESS_ALL);
			switch (execute_result) {
			case TaskExecutionResult::TASK_FINISHED:
			case TaskExecutionResult::TASK_ERROR:
				task.reset();
				break;
			case TaskExecutionResult::TASK_NOT_FINISHED:
				throw InternalException("Task should not return TASK_NOT_FINISHED in PROCESS_ALL mode");
			case TaskExecutionResult::TASK_BLOCKED:
				task->Deschedule();
				task.reset();
				break;
			}
		} catch (...) {
			return;
		}
	}
#else
	throw NotImplementedException("DuckDB was compiled without threads! Background thread loop is not allowed.");
#endif
}

#ifndef DUCKDB_NO_THREADS
static void ThreadExecuteTasks(TaskScheduler *scheduler, atomic<bool> *marker, const idx_t &thread_idx) {
	scheduler->ExecuteForever(marker, thread_idx);
}
#endif

int32_t TaskScheduler::NumberOfThreads() {
	return current_thread_count.load();
}

idx_t TaskScheduler::GetNumberOfTasks() const {
#ifndef DUCKDB_NO_THREADS
	return queue->GetNumberOfTasks();
#else
	idx_t task_count = 0;
	for (auto &producer : queue->q) {
		task_count += producer.second.size();
	}
	return task_count;
#endif
}

idx_t TaskScheduler::GetProducerCount() const {
#ifndef DUCKDB_NO_THREADS
	return queue->GetProducerCount();
#else
	return queue->q.size();
#endif
}

idx_t TaskScheduler::GetTaskCountForProducer(ProducerToken &token) const {
#ifndef DUCKDB_NO_THREADS
	return queue->GetTaskCountForProducer(token);
#else
	const auto it = queue->q.find(std::ref(*token.token));
	if (it == queue->q.end()) {
		return 0;
	}
	return it->second.size();
#endif
}

void TaskScheduler::SetThreads(idx_t total_threads, idx_t external_threads) {
	if (total_threads == 0) {
		throw SyntaxException("Number of threads must be positive!");
	}
#ifndef DUCKDB_NO_THREADS
	if (total_threads < external_threads) {
		throw SyntaxException("Number of threads can't be smaller than number of external threads!");
	}
#else
	if (total_threads != external_threads) {
		throw NotImplementedException(
		    "DuckDB was compiled without threads! Setting total_threads != external_threads is not allowed.");
	}
#endif
	requested_thread_count = NumericCast<int32_t>(total_threads - external_threads);
}

void TaskScheduler::SetAllocatorFlushTreshold(idx_t threshold) {
	allocator_flush_threshold = threshold;
}

void TaskScheduler::SetAllocatorBackgroundThreads(bool enable) {
	allocator_background_threads = enable;
	Allocator::SetBackgroundThreads(enable);
}

void TaskScheduler::Signal(idx_t n) {
#ifndef DUCKDB_NO_THREADS
	queue->Signal(n);
#endif
}

void TaskScheduler::YieldThread() {
#ifndef DUCKDB_NO_THREADS
	std::this_thread::yield();
#endif
}

idx_t TaskScheduler::GetEstimatedCPUId() {
#if defined(EMSCRIPTEN)
	// FIXME: Wasm + multithreads can likely be implemented as
	//   return return (idx_t)std::hash<std::thread::id>()(std::this_thread::get_id());
	return 0;
#else
	// this code comes from jemalloc
#if defined(_WIN32)
	return (idx_t)GetCurrentProcessorNumber();
#elif defined(_GNU_SOURCE)
	auto cpu = sched_getcpu();
	if (cpu < 0) {
#ifndef DUCKDB_NO_THREADS
		// fallback to thread id
		return (idx_t)std::hash<std::thread::id>()(std::this_thread::get_id());
#else

		return 0;
#endif
	}
	return (idx_t)cpu;
#elif defined(__aarch64__) && defined(__APPLE__)
	/* Other oses most likely use tpidr_el0 instead */
	uintptr_t c;
	asm volatile("mrs %x0, tpidrro_el0" : "=r"(c)::"memory");
	return (idx_t)(c & (1 << 3) - 1);
#else
#ifndef DUCKDB_NO_THREADS
	// fallback to thread id
	return (idx_t)std::hash<std::thread::id>()(std::this_thread::get_id());
#else
	return 0;
#endif
#endif
#endif
}

void TaskScheduler::RelaunchThreads() {
	lock_guard<mutex> t(thread_lock);
	auto n = requested_thread_count.load();
	RelaunchThreadsInternal(n);
}

void TaskScheduler::RelaunchThreadsInternal(int32_t n) {
#ifndef DUCKDB_NO_THREADS
	auto &config = DBConfig::GetConfig(db);
	auto new_thread_count = NumericCast<idx_t>(n);
	if (threads.size() == new_thread_count) {
		current_thread_count = NumericCast<int32_t>(threads.size() + config.options.external_threads);
		return;
	}
	if (threads.size() > new_thread_count) {
		// we are reducing the number of threads: clear all threads first
		for (idx_t i = 0; i < threads.size(); i++) {
			*markers[i] = false;
		}
		Signal(threads.size());
		// now join the threads to ensure they are fully stopped before erasing them
		for (idx_t i = 0; i < threads.size(); i++) {
			threads[i]->internal_thread->join();
		}
		// erase the threads/markers
		threads.clear();
		markers.clear();
	}
	if (threads.size() < new_thread_count) {
		// we are increasing the number of threads: launch them and run tasks on them
		idx_t create_new_threads = new_thread_count - threads.size();
		for (idx_t i = 0; i < create_new_threads; i++) {
			// launch a thread and assign it a cancellation marker
			auto marker = unique_ptr<atomic<bool>>(new atomic<bool>(true));
			unique_ptr<thread> worker_thread;
			try {
				worker_thread = make_uniq<thread>(ThreadExecuteTasks, this, marker.get(), threads.size());
			} catch (std::exception &ex) {
				// thread constructor failed - this can happen when the system has too many threads allocated
				// in this case we cannot allocate more threads - stop launching them
				break;
			}
			auto thread_wrapper = make_uniq<SchedulerThread>(std::move(worker_thread));

			threads.push_back(std::move(thread_wrapper));
			markers.push_back(std::move(marker));
		}
	}
	current_thread_count = NumericCast<int32_t>(threads.size() + config.options.external_threads);
	if (Allocator::SupportsFlush()) {
		Allocator::FlushAll();
	}
	// notify the task queue with the new number of threads
	// this causes it to move all group tasks to the global queue again
	queue->SetThreads(threads.size());
#endif
}

} // namespace duckdb
