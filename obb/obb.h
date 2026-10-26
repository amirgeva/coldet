#pragma once

#include <vector>
#include <Eigen/Dense>

namespace obb
{

using Vector3f = Eigen::Matrix<float, 3, 1>;
using Matrix3f = Eigen::Matrix<float, 3, 3>;

struct Transform
{
	Transform(const Matrix3f& R = Matrix3f::Identity(),
			  const Vector3f& T = Vector3f::Zero())
		: rotation(R)
		, center(T)
	{}

	Matrix3f rotation;
	Vector3f center;

	Vector3f apply(const Vector3f& point) const
	{
		return rotation * point + center;
	}
};

struct OBB
{
	Transform transform;
	Vector3f half_size;
};

struct Ray
{
	Vector3f origin;
	Vector3f direction; // should be normalized
};

// Compute a tight bounding oriented box (OBB) for a set of 3D points
OBB compute_obb(const std::vector<Vector3f>& points, const std::vector<uint32_t>& indices);
OBB compute_obb(const std::vector<Vector3f>& points);
OBB compute_obb(const std::vector<OBB>& obbs);
bool contains(const OBB& obb, const Vector3f& point, Vector3f& error);
bool contains_local(const OBB& obb, const Vector3f& local_point);

bool intersects(const OBB& a, const OBB& b);
bool intersects(const OBB& a, const Transform& Ta, const OBB& b, const Transform& Tb);
void obb_corners_self(const OBB& b, std::vector<Vector3f>& out);
void obb_corners_world(const OBB& b, std::vector<Vector3f>& out);
Vector3f to_local(const OBB& obb, const Vector3f& point, bool absolute);

}