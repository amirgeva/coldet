#include <vector>
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
    OBB obb;
    obb.rotation.setIdentity();
    obb.center.setZero();
    obb.half_size.setZero();

    if (points.empty())
        return obb;

    // 1) Centroid
    Vector3f mean = Vector3f::Zero();
    for (const auto& p : points) mean += p;
    mean /= static_cast<float>(points.size());

    if (points.size() == 1) {
        obb.center = mean;
        return obb;
    }

    // 2) Covariance
    Eigen::Matrix3f cov = Eigen::Matrix3f::Zero();
    for (const auto& p : points) {
        Vector3f d = p - mean;
        cov += d * d.transpose();
    }
    cov /= static_cast<float>(points.size());

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

    for (const auto& p : points) {
        Vector3f q = R.transpose() * (p - mean); // local coords about centroid
        minv = minv.cwiseMin(q);
        maxv = maxv.cwiseMax(q);
    }

    Vector3f center_local = 0.5f * (minv + maxv);
    obb.half_size = 0.5f * (maxv - minv);

    // 5) Back to world
    Vector3f center_world = mean + R * center_local;

    obb.rotation = R;
    obb.center = center_world;

    return obb;
}

void obb_corners_world(const OBB& b, std::vector<Vector3f>& out)
{
    const Eigen::Matrix3f R = b.rotation;
    const Vector3f c = b.center;
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
    out.rotation.setIdentity();
    out.center.setZero();
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

bool contains(const OBB& obb, const Vector3f& point)
{
    const Eigen::Matrix3f R = obb.rotation;
    const Vector3f c = obb.center;
    const Vector3f d = point - c;
    const Vector3f local = R.transpose() * d;
    return (std::abs(local.x()) <= obb.half_size.x() &&
            std::abs(local.y()) <= obb.half_size.y() &&
            std::abs(local.z()) <= obb.half_size.z());
}

bool intersects(const OBB& A, const OBB& B)
{
    // OBB-OBB intersection via SAT (15 axes).
    // Assumes A.rotation and B.rotation are orthonormal (pure rotations).

    const float EPS = 1e-6f;

    const Vector3f a = A.half_size;
    const Vector3f b = B.half_size;

    // Compute rotation matrix expressing B in A’s frame
    const Matrix3f R = A.rotation.transpose() * B.rotation;

    // Compute translation vector t from A to B in A’s frame
    const Vector3f t = A.rotation.transpose() * (B.center - A.center);

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



} // namespace obb