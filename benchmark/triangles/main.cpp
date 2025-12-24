#include "../profiler.h"
#include <mytritri.h>

using namespace COLDET;

extern "C" {
int tri_tri_intersect(const float V0[3],const float V1[3],const float V2[3],
                      const float U0[3],const float U1[3],const float U2[3]);
};

int main(int argc, char* argv[])
{
    std::vector<Triangle> triangles;
    const int num_triangles = 1000;
    for (int i = 0; i < num_triangles; ++i)
    {
        triangles.emplace_back(
            Vector3D(static_cast<float>(rand() % 1000), static_cast<float>(rand() % 1000), static_cast<float>(rand() % 1000)),
            Vector3D(static_cast<float>(rand() % 1000), static_cast<float>(rand() % 1000), static_cast<float>(rand() % 1000)),
            Vector3D(static_cast<float>(rand() % 1000), static_cast<float>(rand() % 1000), static_cast<float>(rand() % 1000))
        );
    }
    {
        size_t test_count = (num_triangles * (num_triangles - 1)) / 2;
        PROFILE_SCOPE_N("tri_tri_intersect benchmark", test_count);
        int intersection_count = 0;
        for (size_t i = 0; i < triangles.size(); ++i)
        {
            for (size_t j = i + 1; j < triangles.size(); ++j)
            {
                const Triangle& t1 = triangles[i];
                const Triangle& t2 = triangles[j];
                int res = tri_tri_intersect(t1.v1.get(), t1.v2.get(), t1.v3.get(), t2.v1.get(), t2.v2.get(), t2.v3.get());
                if (res)
                {
                    ++intersection_count;
                }
            }
        }
    }
    ProfilerStatistics::instance().print_statistics();
    return 0;
}