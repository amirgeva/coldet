#include <chrono>
#include "utils.h"


uint64_t get_tick_count()
{
    using clk = std::chrono::steady_clock;
    static auto start = clk::now();
    auto now = clk::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
}
