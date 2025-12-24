#pragma once

#include <vector>
#include <Eigen/Dense>

namespace obb
{

using Vector3f = Eigen::Matrix<float, 3, 1>;
using Matrix3f = Eigen::Matrix<float, 3, 3>;

struct OBB
{
    Matrix3f rotation;
    Vector3f center;
    Vector3f half_size;
};

struct Ray
{
    Vector3f origin;
    Vector3f direction; // should be normalized
};

// Compute a tight bounding oriented box (OBB) for a set of 3D points
OBB compute_obb(const std::vector<Vector3f>& points);
OBB compute_obb(const std::vector<OBB>& obbs);
bool contains(const OBB& obb, const Vector3f& point);

bool intersects(const OBB& a, const OBB& b);
void obb_corners_world(const OBB& b, std::vector<Vector3f>& out);

}