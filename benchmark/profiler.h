#pragma once

#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>

class ProfilerStatistics
{
    struct Sample{
        uint64_t count;
        uint64_t duration;
        Sample& operator+=(const Sample& other)
        {
            count += other.count;
            duration += other.duration;
            return *this;
        }
    };
    std::unordered_map<std::string, std::vector<Sample>> m_Statistics;
public:
    static ProfilerStatistics& instance()
    {
        static ProfilerStatistics instance;
        return instance;
    }

    void add(const char* name, uint64_t duration, size_t count=1)
    {
        auto it = m_Statistics.find(name);
        if (it == m_Statistics.end())
        {
            auto& v = m_Statistics[name];
            v.reserve(100000);
            v.push_back(Sample{count,duration});
        }
        else it->second.push_back(Sample{count,duration});
    }

    void print_statistics() const
    {
        for (const auto& [name, times] : m_Statistics)
        {
            Sample total{0,0};
            for (const Sample& s : times)
            {
                total += s;
            }
            double average = static_cast<double>(total.duration) / total.count;
            printf("Profiler: %s, calls: %zu, average time: %.2f us\n", name.c_str(), times.size(), average);
        }
    }
};

class Profiler
{
public:
    Profiler(const char* name, size_t count = 1)
        : name_(name)
        , count_(count)
        , start_time_(std::chrono::high_resolution_clock::now())
    {}

    ~Profiler()
    {
        auto end_time =  std::chrono::high_resolution_clock::now();
        uint64_t duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time_).count();
        ProfilerStatistics::instance().add(name_, duration, count_);
    }
private:
    const char* name_;
    size_t count_;
    std::chrono::high_resolution_clock::time_point start_time_;
};


#define PROFILE_SCOPE(name) Profiler profiler_##__LINE__(name)
#define PROFILE_SCOPE_N(name, n) Profiler profiler_##__LINE__(name, n)