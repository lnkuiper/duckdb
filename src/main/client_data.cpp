#include "duckdb/main/client_data.hpp"
#include "duckdb/common/multi_file/multi_file_list.hpp"

#include "duckdb/catalog/catalog_search_path.hpp"
#include "duckdb/common/constants.hpp"
#include "duckdb/common/opener_file_system.hpp"
#include "duckdb/common/random_engine.hpp"
#include "duckdb/common/serializer/buffered_file_writer.hpp"
#include "duckdb/main/attached_database.hpp"
#include "duckdb/main/client_config.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/client_context_state.hpp"
#include "duckdb/main/context.hpp"
#include "duckdb/main/client_context_file_opener.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/database_manager.hpp"
#include "duckdb/main/query_profiler.hpp"

namespace duckdb {

ClientData::ClientData(ClientContext &context) : catalog_search_path(make_uniq<CatalogSearchPath>(context)) {
	auto &db = DatabaseInstance::GetDatabase(context);

	default_database = DatabaseManager::Get(db).GetDefaultDatabase();
	temporary_objects = make_shared_ptr<AttachedDatabase>(db, AttachedDatabaseType::TEMP_DATABASE);
	temporary_objects->oid = DatabaseManager::Get(db).NextOid();
	random_engine = make_uniq<RandomEngine>();

	temporary_objects->Initialize();
}

ClientData::~ClientData() {
	// This needs to be destroyed first because PreparedStatementData holds a PhysicalPlan,
	// and the PhysicalPlan holds objects that might rely on the client file opener/file system/buffer manager
	prepared_statements.clear();
}

ClientData &ClientData::Get(ClientContext &context) {
	return *context.client_data;
}

const ClientData &ClientData::Get(const ClientContext &context) {
	return *context.client_data;
}

RandomEngine &RandomEngine::Get(ClientContext &context) {
	return *ClientData::Get(context).random_engine;
}

namespace {
struct RandomQueryState : public ClientContextState {
	mutex lock;
	optional_idx query_id;
	uint64_t seed = 0;
	RandomEngine worker_seeds {0};
};
} // namespace

uint64_t RandomEngine::GetQuerySeed(ClientContext &context, bool per_worker) {
	auto &random_engine = Get(context);
	if (!context.transaction.HasActiveTransaction()) {
		lock_guard<mutex> guard(random_engine.lock);
		return random_engine.NextSeed();
	}
	auto state = context.GetContext().GetRegisteredState().GetOrCreate<RandomQueryState>("random_query_seed");
	lock_guard<mutex> guard(state->lock);
	auto query_id = context.transaction.GetActiveQuery();
	if (!state->query_id.IsValid() || state->query_id.GetIndex() != query_id) {
		lock_guard<mutex> random_guard(random_engine.lock);
		state->seed = random_engine.NextSeed();
		// Worker allocation must not advance the session RNG or reuse the scan streams.
		state->worker_seeds.SetSeed(~state->seed);
		state->query_id = query_id;
	}
	return per_worker ? state->worker_seeds.NextRandomInteger64() : state->seed;
}

} // namespace duckdb
