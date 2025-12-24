#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <gtest/gtest.h>
#include <obb.h>

using namespace obb;


struct OBBTestResult {
    bool contains_all = true;
    float max_violation = 0.0f; // max amount outside (in local units)
    float obb_volume = 0.0f;
    float aabb_volume = 0.0f;
    float volume_ratio_obb_over_aabb = 0.0f; // <= 1 is "better" (smaller than AABB)
};

OBBTestResult test_obb_compactness(
    const OBB& obb,
    const std::vector<Vector3f>& points,
    float eps = 1e-5f)
{
    OBBTestResult r;
    if (points.empty()) return r;

    const Eigen::Matrix3f R = obb.rotation;
    const Vector3f c = obb.center;
    const Vector3f e = obb.half_size;

    // OBB volume
    r.obb_volume = 8.0f * e.x() * e.y() * e.z();

    // AABB of points in world
    Vector3f minw(std::numeric_limits<float>::infinity(),
                  std::numeric_limits<float>::infinity(),
                  std::numeric_limits<float>::infinity());
    Vector3f maxw(-std::numeric_limits<float>::infinity(),
                  -std::numeric_limits<float>::infinity(),
                  -std::numeric_limits<float>::infinity());

    for (const auto& p : points) {
        // containment
        Vector3f q = R.transpose() * (p - c); // local point
        Vector3f d = q.cwiseAbs() - e;        // positive components are outside

        float viol = std::max({d.x(), d.y(), d.z(), 0.0f});
        r.max_violation = std::max(r.max_violation, viol);
        if (viol > eps) r.contains_all = false;

        // AABB
        minw = minw.cwiseMin(p);
        maxw = maxw.cwiseMax(p);
    }

    Vector3f aabb_size = (maxw - minw).cwiseMax(Vector3f::Zero());
    r.aabb_volume = aabb_size.x() * aabb_size.y() * aabb_size.z();

    if (r.aabb_volume > 0.0f)
        r.volume_ratio_obb_over_aabb = r.obb_volume / r.aabb_volume;
    else
        r.volume_ratio_obb_over_aabb = 0.0f; // degenerate case

    return r;
}


TEST(OBB, build)
{
    std::vector<Vector3f> points;
    float s = sqrt(2);
    Matrix3f R; R << s,0,s, 0,1,0, -s,0,s;
                     
    for(size_t i=0; i<1000; ++i)
    {
        Vector3f p;
        p.setRandom();
        p.z()*=5;
        points.push_back(R * p);
    }
    OBB box = compute_obb(points);
    {
        std::ofstream f("pts.csv");
        for(int sx=-1; sx<=1; sx+=2)
        for(int sy=-1; sy<=1; sy+=2)
        for(int sz=-1; sz<=1; sz+=2)
        {
            Vector3f corner = box.center
                            + sx * box.half_size.x() * box.rotation.col(0)
                            + sy * box.half_size.y() * box.rotation.col(1)
                            + sz * box.half_size.z() * box.rotation.col(2);
            f << corner.x() << "," << corner.y() << "," << corner.z() << "\n";
        }
        for (const auto& p : points) {
            f << p.x() << "," << p.y() << "," << p.z() << "\n";
        }
    }
    EXPECT_TRUE(box.half_size.x() > 0);
    EXPECT_TRUE(box.half_size.y() > 0);
    EXPECT_TRUE(box.half_size.z() > 0);
    OBBTestResult res = test_obb_compactness(box, points);
    EXPECT_TRUE(res.contains_all);
    EXPECT_LT(res.max_violation, 1e-5f);
    EXPECT_LT(res.volume_ratio_obb_over_aabb, 0.5f); // OBB should be significantly more compact than AABB
}

static inline Matrix3f x_rotation(float angle_rad)
{
    Matrix3f R;
    float c = cos(angle_rad);
    float s = sin(angle_rad);
    R << 1,0,0,
         0,c,-s,
         0,s,c;
    return R;
}

static inline Matrix3f y_rotation(float angle_rad)
{
    Matrix3f R;
    float c = cos(angle_rad);
    float s = sin(angle_rad);
    R << c,0,s,
         0,1,0,
        -s,0,c;
    return R;
}

static inline Matrix3f z_rotation(float angle_rad)
{
    Matrix3f R;
    float c = cos(angle_rad);
    float s = sin(angle_rad);
    R << c,-s,0,
         s, c,0,
         0, 0,1;
    return R;
}

static inline Matrix3f random_rotation()
{
    static std::mt19937 rng(12345);
    static std::uniform_real_distribution<float> dist(0.0f, 2.0f * 3.14159265f);

    float ax = dist(rng);
    float ay = dist(rng);
    float az = dist(rng);

    return z_rotation(az) * y_rotation(ay) * x_rotation(ax);
}

struct range
{
    float min=-std::numeric_limits<float>::infinity();
    float max=std::numeric_limits<float>::infinity();

    bool valid() const { return min <= max; }

    range& operator&= (const range& other)
    {
        min = std::max(min, other.min);
        max = std::min(max, other.max);
        return *this;
    }

    range operator& (const range& other) const
    {
        range r = *this;
        return r &= other;
    }
};

static bool intersects(const Vector3f& half_size, const Ray& ray, float max_t)
{
    range total;
    for(Eigen::Index i=0;i<3;++i)
    {
        if (fabs(ray.direction[i]) > 1e-6f)
        {
            float t1 = (-half_size[i] - ray.origin[i]) / ray.direction[i];
            float t2 = ( half_size[i] - ray.origin[i]) / ray.direction[i];
            if (t1 > t2) std::swap(t1, t2);
            range r;
            r.min = t1;
            r.max = t2;
            total &= r;
        }
    }
    return total.valid() && total.min;
}

static inline bool brute_force_intersection(const OBB& a, const OBB& b)
{
    OBB b_in_a = b;
    // Transform b into a's local space
    b_in_a.rotation = a.rotation.transpose() * b.rotation;
    b_in_a.center = a.rotation.transpose() * (b.center - a.center);
    // Now a is axis-aligned at origin
    Vector3f a_corners[8],b_corners[8];
    {
        const Eigen::Matrix3f R = b_in_a.rotation;
        const Vector3f c = b_in_a.center;
        const Vector3f e = b_in_a.half_size;
        int idx = 0;
        for (int sx : {-1, +1})
        for (int sy : {-1, +1})
        for (int sz : {-1, +1})
        {
            a_corners[idx] = Vector3f(
                sx * a.half_size.x(),
                sy * a.half_size.y(),
                sz * a.half_size.z());
            b_corners[idx++] = c
                            + sx * e.x() * R.col(0)
                            + sy * e.y() * R.col(1)
                            + sz * e.z() * R.col(2);
        }
    }
    for(int i=0; i<8; ++i)
    {
        if (contains(a, b_corners[i]))
            return true;
        if (contains(b_in_a, a_corners[i]))
            return true;
    }
    for(int i=0; i<8; ++i)
    {
        Vector3f origin=b_corners[i];
        for(int j=i+1;j<8;++j)
        {
            Vector3f direction=(b_corners[j] - origin);
            Ray ray{origin, direction.normalized()};
            if (intersects(a.half_size, ray, direction.norm()))
                return true;
        }
    }
    return false;
}

TEST(OBB, intersection)
{
    size_t total_intersections=0;
    for(size_t iter=0;iter<1000;++iter)
    {
        OBB a,b;
        a.center.setRandom(); a.center *= 10.0f;
        b.center.setRandom(); b.center *= 10.0f;
        a.half_size.setRandom(); a.half_size = a.half_size.cwiseAbs() + Vector3f::Constant(0.5f);
        b.half_size.setRandom(); b.half_size = b.half_size.cwiseAbs() + Vector3f::Constant(0.5f);
        a.rotation = random_rotation();
        b.rotation = random_rotation();
        bool bf = brute_force_intersection(a,b);
        bool opt = intersects(a,b);
        if (opt!=bf)
        {
            std::cout << "Failed intersection test on iteration " << iter << std::endl;
            std::vector<Vector3f> corners;
            obb_corners_world(a, corners);
            obb_corners_world(b, corners);
            std::ofstream f("/mnt/c/tmp/pts.csv");
            for (const auto& p : corners)
            {
                f << p.x() << "," << p.y() << "," << p.z() << "\n";
            }
        }
        total_intersections += opt ? 1 : 0;
        ASSERT_EQ(bf,opt);
    }
    std::cout << "Total intersections: " << total_intersections << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
