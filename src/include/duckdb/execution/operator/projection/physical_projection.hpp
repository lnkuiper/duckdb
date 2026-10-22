//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/execution/operator/projection/physical_projection.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include "duckdb/execution/physical_operator.hpp"
#include "duckdb/planner/expression.hpp"

namespace duckdb {

class PhysicalProjection : public PhysicalOperator {
public:
	OperatorInitializationInfo GetInitializationInfo() const override;
	BatchMetadataInfo RequestedMetadata() const override;
	static constexpr const PhysicalOperatorType TYPE = PhysicalOperatorType::PROJECTION;

public:
	PhysicalProjection(PhysicalPlan &physical_plan, vector<LogicalType> types,
	                   vector<unique_ptr<Expression>> select_list, idx_t estimated_cardinality);

	vector<unique_ptr<Expression>> select_list;

public:
	BatchMetadataInfo OutputMetadata(const BatchMetadataInfo &input) const override {
		return input;
	}
	void InitializeOutputChunk(Allocator &allocator, DataChunk &chunk) const override;
	unique_ptr<OperatorState> GetOperatorState(ExecutionContext &context) const override;
	unique_ptr<GlobalOperatorState> GetGlobalOperatorState(ClientContext &context) const override {
		return GetGlobalOperatorState(context, OperatorStateInitInput());
	}
	unique_ptr<GlobalOperatorState> GetGlobalOperatorState(ClientContext &context,
	                                                       const OperatorStateInitInput &input) const override;
	OperatorResultType Execute(ExecutionContext &context, OperatorChunk &input, OperatorChunk &chunk,
	                           GlobalOperatorState &gstate, OperatorState &state) const override;

	bool ParallelOperator() const override;
	PipelineExternalInputSupport GetExternalInputSupport() const override {
		return PipelineExternalInputSupport::SUPPORTED;
	}

	InsertionOrderPreservingMap<string> ParamsToString() const override;
};

} // namespace duckdb
