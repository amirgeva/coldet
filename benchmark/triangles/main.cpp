#include "../profiler.h"
#include <mytritri.h>
#include <tritri.h>

using namespace COLDET;

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
    const size_t test_count = (num_triangles * (num_triangles - 1)) / 2;
    int plain_count = 0, shape_count = 0;
    {
        PROFILE_SCOPE_N("tri_tri_intersect", test_count);
        for (size_t i = 0; i < triangles.size(); ++i)
            for (size_t j = i + 1; j < triangles.size(); ++j)
            {
                const Triangle& t1 = triangles[i];
                const Triangle& t2 = triangles[j];
                if (tri_tri_intersect(t1.v1.get(), t1.v2.get(), t1.v3.get(), t2.v1.get(), t2.v2.get(), t2.v3.get()))
                    ++plain_count;
            }
    }
    {
        // with the shapes the triangles computed when they were built
        PROFILE_SCOPE_N("tri_tri_intersect with shapes", test_count);
        for (size_t i = 0; i < triangles.size(); ++i)
            for (size_t j = i + 1; j < triangles.size(); ++j)
            {
                const Triangle& t1 = triangles[i];
                const Triangle& t2 = triangles[j];
                if (tri_tri_intersect(t1.v1.get(), t1.v2.get(), t1.v3.get(), t1.shape,
                                      t2.v1.get(), t2.v2.get(), t2.v3.get(), t2.shape))
                    ++shape_count;
            }
    }
    printf("%d intersecting pairs (%d with shapes)\n", plain_count, shape_count);
    ProfilerStatistics::instance().print_statistics();
    return 0;
}