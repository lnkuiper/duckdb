#include "duckdb/reoptimizer/reoptimizer.hpp"

#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/optimizer/optimizer.hpp"

using namespace std;
using namespace duckdb;

ReOptimizer::ReOptimizer(Binder &binder, ClientContext &context) : context(context), binder(binder) {
}
