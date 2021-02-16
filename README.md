# Re-Optimization in DuckDB
Master thesis project on implementing re-optimization in DuckDB, evaluated on the Join Order Benchmark, inspired by the paper by Perron et al. "How I Learned to Stop Worrying and Love Re-Optimization" which did the same in PostgreSQL.

There are significant differences between PostgreSQL and DuckDB that made this an interesting piece of research.
Many re-optimization schemes were evaluated.
All of them yielded an improved run-time for the longest-running queries.

Most of my work can be found in `src/reoptimizer/reoptimizer.cpp` and the corresponding header file, save for some modifications to LogicalOperators and Query Profiling).

Last upstream merge commit id: `c65bab0d7234c62807f6e28df71ad68c25929b00`

## Usage
The project is built with CMake, using the `make` command. The IMDB data is downloaded using `make imdb`.

Re-optimization is enabled by setting `ClientContext::enable_reoptimizer` to `true`. Different re-optimization algorithms are available by modifying line `reoptimizer.cpp:41`. The re-optimization scheme by Perron et al. can be used by changing `client_context.cpp:192` to `plan = reoptimizer.SimulatedReOptimize(move(plan), query);`.
