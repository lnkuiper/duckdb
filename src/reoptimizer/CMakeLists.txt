add_library_unity(duckdb_reoptimizer
                  OBJECT
                  reoptimizer.cpp)
set(ALL_OBJECT_FILES
    ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:duckdb_reoptimizer>
    PARENT_SCOPE)
