#include <vector>
#include <numeric>
#include <limits>
#include <algorithm>
#include "obb.h"

namespace obb
{

static inline void sort_eigen_desc(const Vector3f& evals, int idx[3]) {
	idx[0] = 0; idx[1] = 1; idx[2] = 2;
	std::sort(idx, idx + 3, [&](int a, int b){ return evals[a] > evals[b]; });
}

OBB compute_obb(const std::vector<Vector3f>& points)
{
	std::vector<uint32_t> indices(points.size());
	std::iota(indices.begin(), indices.end(), 0);
	return compute_obb(points, indices);
}

OBB compute_obb(const std::vector<Vector3f>& points, const std::vector<uint32_t>& indices)
{
	OBB obb;
	obb.transform.rotation.setIdentity();
	obb.transform.center.setZero();
	obb.half_size.setZero();

	if (indices.empty())
		return obb;

	// 1) Centroid
	Vector3f mean = Vector3f::Zero();
	for (const auto& i : indices) mean += points[i];
	mean /= static_cast<float>(indices.size());

	if (indices.size() == 1) {
		obb.transform.center = mean;
		return obb;
	}

	// 2) Covariance
	Eigen::Matrix3f cov = Eigen::Matrix3f::Zero();
	for (const auto& i : indices) {
		Vector3f d = points[i] - mean;
		cov += d * d.transpose();
	}
	cov /= static_cast<float>(indices.size());

	// 3) PCA axes (eigenvectors of symmetric covariance)
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es(cov);
	Eigen::Matrix3f evecs = es.eigenvectors(); // columns correspond to ascending eigenvalues
	Vector3f evals = es.eigenvalues();

	// Sort eigenvectors by descending eigenvalue (largest variance first)
	int order[3];
	sort_eigen_desc(evals, order);

	Eigen::Matrix3f R;
	R.col(0) = evecs.col(order[0]).normalized();
	R.col(1) = evecs.col(order[1]).normalized();
	R.col(2) = evecs.col(order[2]).normalized();

	// Ensure right-handed rotation (determinant +1)
	if (R.determinant() < 0.0f) {
		R.col(2) *= -1.0f;
	}

	// 4) Project points into OBB local space and find tight AABB there
	Vector3f minv(std::numeric_limits<float>::infinity(),
				  std::numeric_limits<float>::infinity(),
				  std::numeric_limits<float>::infinity());
	Vector3f maxv(-std::numeric_limits<float>::infinity(),
				  -std::numeric_limits<float>::infinity(),
				  -std::numeric_limits<float>::infinity());

	for (const auto& i : indices) {
		Vector3f q = R.transpose() * (points[i] - mean); // local coords about centroid
		minv = minv.cwiseMin(q);
		maxv = maxv.cwiseMax(q);
	}

	Vector3f center_local = 0.5f * (minv + maxv);
	static const Vector3f margin = 1.0e-3 * Vector3f(1, 1, 1);
	obb.half_size = 0.5f * (maxv - minv) + margin;

	// 5) Back to world
	Vector3f center_world = mean + R * center_local;

	obb.transform.rotation = R;
	obb.transform.center = center_world;

	return obb;
}

void obb_corners_self(const OBB& b, std::vector<Vector3f>& out)
{
	const Vector3f e = b.half_size;
	for (int sx : {-1, +1})
	for (int sy : {-1, +1})
	for (int sz : {-1, +1})
	{
		Vector3f local(sx * e.x(), sy * e.y(), sz * e.z());
		out.push_back(local);
	}
}

void obb_corners_world(const OBB& b, std::vector<Vector3f>& out)
{
	const Eigen::Matrix3f R = b.transform.rotation;
	const Vector3f c = b.transform.center;
	const Vector3f e = b.half_size;
	for (int sx : {-1, +1})
	for (int sy : {-1, +1})
	for (int sz : {-1, +1})
	{
		Vector3f local(sx * e.x(), sy * e.y(), sz * e.z());
		out.push_back(c + R * local);
	}
}


OBB compute_obb(const std::vector<OBB>& obbs)
{
	OBB out;
	out.transform.rotation.setIdentity();
	out.transform.center.setZero();
	out.half_size.setZero();

	if (obbs.empty())
		return out;

	std::vector<Vector3f> pts;
	pts.reserve(obbs.size() * 8);

	for (const auto& b : obbs)
	{
		obb_corners_world(b, pts); // existing helper
	}

	return compute_obb(pts); // existing helper (PCA + tight extents)
}

bool contains(const OBB& obb, const Vector3f& point, Vector3f& error)
{
	const Vector3f local = to_local(obb, point, true);
	error = local - obb.half_size;
	return error.maxCoeff() <= 0.0f;
}

bool contains_local(const OBB& obb, const Vector3f& local)
{
	return (std::abs(local.x()) <= obb.half_size.x() &&
			std::abs(local.y()) <= obb.half_size.y() &&
			std::abs(local.z()) <= obb.half_size.z());
}

bool intersects(const OBB& a, const Transform& Ta, const OBB& b, const Transform& Tb)
{
	// Create temporary transformed OBBs
	OBB transformed_a;
	transformed_a.transform.center = Ta.apply(a.transform.center);
	transformed_a.transform.rotation = Ta.rotation * a.transform.rotation;
	transformed_a.half_size = a.half_size;

	OBB transformed_b;
	transformed_b.transform.center = Tb.apply(b.transform.center);
	transformed_b.transform.rotation = Tb.rotation * b.transform.rotation;
	transformed_b.half_size = b.half_size;

	// Use the existing OBB-OBB intersection function
	return intersects(transformed_a, transformed_b);
}

bool intersects(const OBB& A, const OBB& B)
{
	// OBB-OBB intersection via SAT (15 axes).
	// Assumes A.transform.rotation and B.transform.rotation are orthonormal (pure rotations).

	const float EPS = 1e-6f;

	const Vector3f a = A.half_size;
	const Vector3f b = B.half_size;

	// Compute rotation matrix expressing B in A’s frame
	const Matrix3f R = A.transform.rotation.transpose() * B.transform.rotation;

	// Compute translation vector t from A to B in A’s frame
	const Vector3f t = A.transform.rotation.transpose() * (B.transform.center - A.transform.center);

	// Compute common subexpressions. Add epsilon to counteract arithmetic errors
	Matrix3f AbsR = R.cwiseAbs();
	AbsR.array() += EPS;

	float ra, rb;

	// Test axes L = A0, A1, A2
	for (int i = 0; i < 3; ++i) {
		ra = a[i];
		rb = b[0] * AbsR(i, 0) + b[1] * AbsR(i, 1) + b[2] * AbsR(i, 2);
		if (std::abs(t[i]) > ra + rb) return false;
	}

	// Test axes L = B0, B1, B2
	for (int j = 0; j < 3; ++j) {
		ra = a[0] * AbsR(0, j) + a[1] * AbsR(1, j) + a[2] * AbsR(2, j);
		rb = b[j];
		const float tj = std::abs(t[0] * R(0, j) + t[1] * R(1, j) + t[2] * R(2, j));
		if (tj > ra + rb) return false;
	}

	// Test axis L = Ai x Bj (9 tests)
	for (int i = 0; i < 3; ++i) {
		const int i1 = (i + 1) % 3;
		const int i2 = (i + 2) % 3;

		for (int j = 0; j < 3; ++j) {
			const int j1 = (j + 1) % 3;
			const int j2 = (j + 2) % 3;

			// |t · (Ai x Bj)| = | t_i2 * R(i1,j) - t_i1 * R(i2,j) |
			const float tproj = std::abs(t[i2] * R(i1, j) - t[i1] * R(i2, j));

			ra = a[i1] * AbsR(i2, j) + a[i2] * AbsR(i1, j);
			rb = b[j1] * AbsR(i, j2) + b[j2] * AbsR(i, j1);

			if (tproj > ra + rb) return false;
		}
	}

	return true; // no separating axis found
}

Vector3f to_local(const OBB& obb, const Vector3f& point, bool absolute)
{
	Vector3f local = obb.transform.rotation.transpose() * (point - obb.transform.center);
	if (absolute) {
		local.x() = std::abs(local.x());
		local.y() = std::abs(local.y());
		local.z() = std::abs(local.z());
	}
	return local;
}


} // namespace obb