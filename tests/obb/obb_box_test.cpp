// Tests of the OBB primitives: fitting, containment and box / box intersection
#include <gtest/gtest.h>
#include "test_helpers.h"

namespace
{

// Largest distance of a point outside the box (0 when all are inside)
float max_violation(const OBB& box, const std::vector<Vector3f>& points)
{
	float worst = 0;
	for (const auto& p : points)
	{
		Vector3f d = to_local(box, p, true) - box.half_size;
		worst = std::max(worst, d.maxCoeff());
	}
	return worst;
}

OBB make_box(const Matrix3f& R, const Vector3f& center, const Vector3f& half)
{
	OBB b;
	b.transform = Transform{ R, center };
	b.half_size = half;
	return b;
}

} // namespace

TEST(obb_fit, SinglePoint)
{
	OBB box = compute_obb(std::vector<Vector3f>{ Vector3f(1, 2, 3) });
	EXPECT_EQ(box.half_size, Vector3f::Zero());
	EXPECT_EQ(box.transform.center, Vector3f(1, 2, 3));
}

TEST(obb_fit, EmptyInput)
{
	OBB box = compute_obb(std::vector<Vector3f>{});
	EXPECT_EQ(box.half_size, Vector3f::Zero());
}

TEST(obb_fit, CollinearPoints)
{
	const Vector3f dir = Vector3f(1, 2, -1).normalized(), base(3, -1, 2);
	std::vector<Vector3f> points;
	for (int i = 0; i <= 10; ++i) points.push_back(base + dir * float(i));
	OBB box = compute_obb(points);
	std::vector<float> h(box.half_size.data(), box.half_size.data() + 3);
	std::sort(h.begin(), h.end());
	EXPECT_NEAR(h[2], 5.0f, 1e-5f) << "half the segment's length";
	EXPECT_NEAR(h[1], 0.0f, 1e-5f);
	EXPECT_NEAR(h[0], 0.0f, 1e-5f);
	EXPECT_LT(max_violation(box, points), 1e-4f);
	EXPECT_NEAR(box.transform.rotation.determinant(), 1.0f, 1e-5f);
}

TEST(obb_fit, CoplanarPoints)
{
	std::mt19937 rng(1);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	const Matrix3f R = random_rotation();
	std::vector<Vector3f> points;
	for (int i = 0; i < 200; ++i) points.push_back(R * Vector3f(unit(rng) * 3, unit(rng), 0) + Vector3f(5, 5, 5));
	OBB box = compute_obb(points);
	EXPECT_NEAR(box.half_size.minCoeff(), 0.0f, 1e-5f) << "flat box";
	EXPECT_LT(max_violation(box, points), 1e-4f);
}

// The box always contains its points, whatever their distribution
TEST(obb_fit, ContainsRandomClouds)
{
	std::mt19937 rng(2);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	for (int k = 0; k < 200; ++k)
	{
		const Matrix3f R = random_rotation();
		const Vector3f scale(std::exp(unit(rng) * 3), std::exp(unit(rng) * 3), std::exp(unit(rng) * 3));
		const Vector3f center = Vector3f(unit(rng), unit(rng), unit(rng)) * 100.0f;
		std::vector<Vector3f> points;
		const int n = 2 + k % 50;
		for (int i = 0; i < n; ++i)
			points.push_back(R * Vector3f(unit(rng), unit(rng), unit(rng)).cwiseProduct(scale) + center);
		OBB box = compute_obb(points);
		EXPECT_LT(max_violation(box, points), 1e-4f * (1.0f + center.norm())) << "cloud " << k;
		EXPECT_NEAR(box.transform.rotation.determinant(), 1.0f, 1e-4f);
		EXPECT_TRUE((box.transform.rotation.transpose() * box.transform.rotation).isIdentity(1e-4f));
	}
}

// Fitting a subset by index is the same as fitting a copy of the subset
TEST(obb_fit, IndexedSubset)
{
	std::mt19937 rng(3);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	std::vector<Vector3f> points;
	for (int i = 0; i < 100; ++i) points.emplace_back(unit(rng), unit(rng) * 3, unit(rng) * 0.5f);
	std::vector<uint32_t> indices;
	std::vector<Vector3f> subset;
	for (uint32_t i = 0; i < 100; i += 3) { indices.push_back(i); subset.push_back(points[i]); }
	OBB a = compute_obb(points, indices), b = compute_obb(subset);
	EXPECT_TRUE(a.half_size.isApprox(b.half_size, 1e-5f));
	EXPECT_TRUE(a.transform.center.isApprox(b.transform.center, 1e-5f));
}

// A box fit to other boxes contains all their corners
TEST(obb_fit, FromBoxes)
{
	std::mt19937 rng(4);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	std::vector<OBB> boxes;
	std::vector<Vector3f> corners;
	for (int i = 0; i < 10; ++i)
	{
		boxes.push_back(make_box(random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 5.0f,
								 Vector3f(unit(rng), unit(rng), unit(rng)).cwiseAbs() + Vector3f::Constant(0.1f)));
		obb_corners_world(boxes.back(), corners);
	}
	OBB all = compute_obb(boxes);
	EXPECT_LT(max_violation(all, corners), 1e-4f);
}

TEST(obb_fit, ContainmentHelpers)
{
	OBB box = make_box(z_rotation(0.5f), Vector3f(1, 2, 3), Vector3f(1, 2, 0.5f));
	Vector3f error;
	EXPECT_TRUE(contains(box, box.transform.apply(Vector3f(0.9f, -1.9f, 0.4f)), error));
	EXPECT_FALSE(contains(box, box.transform.apply(Vector3f(1.1f, 0, 0)), error));
	EXPECT_NEAR(error.x(), 0.1f, 1e-5f) << "error is the distance outside, per local axis";
	EXPECT_TRUE(contains_local(box, Vector3f(-1, 2, -0.5f))) << "corners are inside";
	EXPECT_FALSE(contains_local(box, Vector3f(0, 0, 0.6f)));
	EXPECT_TRUE(to_local(box, box.transform.apply(Vector3f(0.25f, -0.5f, 0.125f)), false).isApprox(Vector3f(0.25f, -0.5f, 0.125f), 1e-5f));
	std::vector<Vector3f> self, world;
	obb_corners_self(box, self);
	obb_corners_world(box, world);
	ASSERT_EQ(self.size(), 8u);
	ASSERT_EQ(world.size(), 8u);
	for (int i = 0; i < 8; ++i)
		EXPECT_TRUE(world[i].isApprox(box.transform.apply(self[i]), 1e-5f));
}

// With diagnostics on, every node is verified while the tree is built
TEST(obb_tree, BuildWithDiagnostics)
{
	set_diagnostics(true);
	for (size_t leaf : { 1, 4, 16 })
	{
		OBBTree tree;
		for (const auto& t : torus(2.0f, 0.5f, 0.02f)) tree.add_triangle(t);
		EXPECT_NO_THROW(tree.build(leaf)) << "leaf size " << leaf;
		EXPECT_TRUE(tree.verify());
	}
	set_diagnostics(false);
}

// Box / box intersection with parallel axes: every axis-aligned rotation, at contact and gaps
TEST(obb_intersect, AxisAlignedRotations)
{
	const Vector3f half(1, 1, 1);
	struct Case { const char* name; Vector3f offset; bool hit; };
	const Case cases[] = {
		{ "overlap", { 1.5f, 0.25f, 0 }, true },
		{ "face contact", { 2, 0.5f, 0.25f }, true },
		{ "edge contact", { 2, 2, 0 }, true },
		{ "corner contact", { 2, 2, 2 }, true },
		{ "face gap", { 2.125f, 0, 0 }, false },
		{ "edge gap", { 2.125f, 2.125f, 0 }, false },
		{ "corner gap", { 2.125f, 2.125f, 2.125f }, false },
	};
	for (const Matrix3f& Ra : axis_rotations())
		for (const Matrix3f& Rb : axis_rotations())
			for (const auto& c : cases)
			{
				OBB a = make_box(Ra, Vector3f::Zero(), half), b = make_box(Rb, c.offset, half);
				EXPECT_EQ(intersects(a, b), c.hit) << c.name;
				EXPECT_EQ(intersects(b, a), c.hit) << c.name << " (swapped)";
			}
}

// Boxes without thickness (single flat triangles) in contact, near and far from the origin
TEST(obb_intersect, FlatBoxes)
{
	for (float distance : { 0.0f, 1e3f, 1e4f })
		for (int k = 0; k < 50; ++k)
		{
			const Matrix3f R = random_rotation();
			const Vector3f base = R * Vector3f(1, -2, 3).normalized() * distance;
			// two flat boxes in the same plane, overlapping
			OBB a = make_box(R, base, Vector3f(1, 1, 0));
			OBB b = make_box(R * z_rotation(0.3f), base + R * Vector3f(1.5f, 0, 0), Vector3f(1, 1, 0));
			EXPECT_TRUE(intersects(a, b)) << "coplanar flat boxes, distance " << distance;
			// stacked: a flat box on top of a solid one
			OBB c = make_box(R, base + R * Vector3f(0, 0, -1), Vector3f(1, 1, 1));
			EXPECT_TRUE(intersects(a, c)) << "flat box resting on a box, distance " << distance;
			// lifted clearly off the plane
			OBB d = make_box(R, base + R * Vector3f(0.5f, 0, 0.05f), Vector3f(1, 1, 0));
			EXPECT_FALSE(intersects(a, d)) << "parallel flat boxes 0.05 apart, distance " << distance;
		}
}

// The transformed variant is the same as transforming both boxes first
TEST(obb_intersect, TransformedVariant)
{
	std::mt19937 rng(9);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	size_t hits = 0;
	for (int k = 0; k < 2000; ++k)
	{
		OBB a = make_box(random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)), Vector3f(unit(rng), unit(rng), unit(rng)).cwiseAbs() + Vector3f::Constant(0.1f));
		OBB b = make_box(random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)), Vector3f(unit(rng), unit(rng), unit(rng)).cwiseAbs() + Vector3f::Constant(0.1f));
		Transform Ta{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 1.5f };
		Transform Tb{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 1.5f };
		OBB wa = make_box(Ta.rotation * a.transform.rotation, Ta.apply(a.transform.center), a.half_size);
		OBB wb = make_box(Tb.rotation * b.transform.rotation, Tb.apply(b.transform.center), b.half_size);
		bool expected = intersects(wa, wb);
		EXPECT_EQ(intersects(a, Ta, b, Tb), expected);
		hits += expected ? 1 : 0;
	}
	EXPECT_GT(hits, 100u);
}
