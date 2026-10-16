"""Emit private execution helpers for proven fixed zero/one-child schedules."""


def hooks(rules):
    result = {}
    for _, initialize, finalize in rules.values():
        result[initialize] = "void"
        result[finalize] = "arena_ptr<TransformResultValue>"
    return result


SINGLE_CHILD_STATE = """
namespace {

struct SingleChildTransformState {
	SingleChildTransformState(ParseResult &parse_result, const TransformFrameOps &info,
	                          arena_ptr<TransformResultValue> child_result)
	    : parse_result(parse_result), info(info), child_slots(child_result ? 1 : 0) {
		child_results[0] = std::move(child_result);
	}

	void ReserveChildSlots(idx_t count) {
		D_ASSERT(count <= 1);
		child_slots = count;
	}
	void PushChild(TransformInput input, idx_t slot) {
		D_ASSERT(slot < child_slots);
		D_ASSERT(!child);
		child.emplace(input);
	}
	template <class T>
	T TakeResult(idx_t slot) {
		return TakeTransformChildResult<T>({child_results, child_slots}, info, slot);
	}

	ParseResult &parse_result;
	const TransformFrameOps &info;
	arena_ptr<TransformResultValue> child_results[1];
	optional<TransformInput> child;
	idx_t child_slots;
};

} // namespace
"""


def single_child_source_declarations(rules):
    lines = [
        SINGLE_CHILD_STATE,
        "namespace {",
        "struct GeneratedSingleChildTransform {",
    ]
    for name, result in hooks(rules).items():
        lines.append(
            f"template <class PROCESS> static {result} {name}("
            "PEGTransformer &transformer, ParseResult &parse_result, PROCESS &process);"
        )
    for rule in rules:
        lines.append(
            f"static TransformStep Execute{rule}(PEGTransformer &transformer, "
            "ParseResult &parse_result, arena_ptr<TransformResultValue> child_result);"
        )
    lines.extend(["};", "} // namespace", ""])
    return lines


def single_child_source_definitions(rules):
    lines = []
    for name, result in hooks(rules).items():
        forward = "" if result == "void" else "return "
        lines.extend(
            [
                f"{result} PEGTransformerFactory::{name}(PEGTransformer &transformer, GeneratedTransformProcess &process) {{",
                f"    {forward}GeneratedSingleChildTransform::{name}(transformer, process.parse_result, process);",
                "}",
                "",
            ]
        )
    for rule, (ops, initialize, finalize) in rules.items():
        lines.append(
            f"""
TransformStep GeneratedSingleChildTransform::Execute{rule}(PEGTransformer &transformer, ParseResult &parse_result,
                                                         arena_ptr<TransformResultValue> child_result) {{
    SingleChildTransformState process(parse_result, {ops}, std::move(child_result));
    if (!process.child_results[0]) {{
        {initialize}(transformer, parse_result, process);
        if (process.child) {{
            return TransformStep::Child(*process.child);
        }}
    }}
    auto result = {finalize}(transformer, parse_result, process);
    D_ASSERT(!process.child);
    return TransformStep::Complete(std::move(result));
}}
"""
        )
    return lines
