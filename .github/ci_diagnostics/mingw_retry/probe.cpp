#include <chrono>
#include <cstdio>
#include <thread>
#ifdef _WIN32
#include <windows.h>
#endif

static void Probe() {
    using clock = std::chrono::steady_clock;
    constexpr auto wait = std::chrono::milliseconds(10);
    auto shortest = clock::duration::max();
    unsigned early = 0;
    for (unsigned i = 0; i < 100; ++i) {
        const auto start = clock::now();
        std::this_thread::sleep_for(wait);
        const auto elapsed = clock::now() - start;
        if (elapsed < shortest) {
            shortest = elapsed;
        }
        early += elapsed < wait;
    }
    std::printf("shortest wait: %.9f ms; waits below 10 ms: %u/100\n",
                std::chrono::duration<double, std::milli>(shortest).count(), early);
}

int main() {
#ifdef __VERSION__
    std::printf("compiler: %s\n", __VERSION__);
#endif
#ifdef _GLIBCXX_USE_NANOSLEEP
    std::puts("libstdc++ uses nanosleep");
#endif
#ifdef _GLIBCXX_USE_CLOCK_MONOTONIC
    std::puts("libstdc++ uses CLOCK_MONOTONIC");
#endif
    std::puts("default timer period:");
    std::thread worker(Probe);
    worker.join();
#ifdef _WIN32
    const auto module = LoadLibraryA("winmm.dll");
    if (module) {
        using PeriodFunction = UINT (WINAPI *)(UINT);
        auto begin_period = reinterpret_cast<PeriodFunction>(GetProcAddress(module, "timeBeginPeriod"));
        auto end_period = reinterpret_cast<PeriodFunction>(GetProcAddress(module, "timeEndPeriod"));
        if (begin_period && end_period && begin_period(1) == 0) {
            std::puts("1 ms timer period (standalone probe only):");
            std::thread high_resolution_worker(Probe);
            high_resolution_worker.join();
            end_period(1);
        }
        FreeLibrary(module);
    }
#endif
}
