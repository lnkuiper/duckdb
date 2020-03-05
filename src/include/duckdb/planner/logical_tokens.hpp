//===----------------------------------------------------------------------===//
//                         DuckDB
//
// duckdb/planner/logical_tokens.hpp
//
//
//===----------------------------------------------------------------------===//

#pragma once

namespace duckdb {

class LogicalOperator;

class LogicalAggregate;
class LogicalAnyJoin;
class LogicalChunkGet;
class LogicalComparisonJoin;
class LogicalCreate;
class LogicalCreateTable;
class LogicalCreateIndex;
class LogicalCrossProduct;
class LogicalDelete;
class LogicalDelimGet;
class LogicalDelimJoin;
class LogicalDistinct;
class LogicalEmptyResult;
class LogicalExpressionGet;
class LogicalFilter;
class LogicalGet;
class LogicalIndexScan;
class LogicalJoin;
class LogicalLimit;
class LogicalOrder;
class LogicalTopN;
class LogicalProjection;
class LogicalInsert;
class LogicalCopyFromFile;
class LogicalCopyToFile;
class LogicalExplain;
class LogicalSetOperation;
class LogicalUpdate;
class LogicalTableFunction;
class LogicalPrepare;
class LogicalPruneColumns;
class LogicalWindow;
class LogicalExecute;
class LogicalSimple;
class LogicalRecursiveCTE;
class LogicalCTERef;

} // namespace duckdb
