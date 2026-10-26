// Structured scenes run through COLDET and the OBB tree, checked against brute force
#include <gtest/gtest.h>
#include "test_helpers.h"

namespace
{

const float PI = 3.14159265f;

Transform translate(float x, float y, float z) { return Transform{ Matrix3f::Identity(), Vector3f(x, y, z) }; }

// Applies T after the transform U: (T*U)(p) = T(U(p))
Transform compose(const Transform& T, const Transform& U)
{
	return Transform{ T.rotation * U.rotation, T.rotation * U.center + T.center };
}

bool has_pair(const PairList& pairs, size_t i, size_t j)
{
	return std::binary_search(pairs.begin(), pairs.end(), TrianglePair(i, j));
}

} // namespace

// Face to face contact: a box resting on a tessellated floor (coplanar triangles)
TEST(scene, BoxOnFloor)
{
	const auto floor = grid_mesh(8, 5.0f);
	const auto box = box_mesh(Vector3f(1, 1, 1));
	struct Placement { const char* name; Transform world; };
	const Placement worlds[] = {
		{ "identity", Transform{} },
		{ "rotated 90 about z", Transform{ z_rotation(PI / 2), Vector3f::Zero() } },
		{ "random rotation", Transform{ random_rotation(), Vector3f(0.3f, -0.2f, 0.1f) } },
		{ "far from origin", Transform{ random_rotation(), Vector3f(1000, -2000, 500) } },
	};
	for (const auto& w : worlds)
		for (float spin : { 0.0f, PI / 4, 0.3f })
		{
			SCOPED_TRACE(std::string(w.name) + ", box spin " + std::to_string(spin));
			const Transform box_on_floor{ z_rotation(spin), Vector3f(0.25f, 0.5f, 1.0f) };
			PairList resting = check_scene(floor, w.world, box, compose(w.world, box_on_floor));
			EXPECT_FALSE(resting.empty()) << "resting box must touch the floor";
			if (w.world.center.norm() < 10.0f)
			{
				// Near the origin, the placement is precise: COLDET (which works in the floor's
				// frame) must find the contacts of the exact geometry, apart from borderline ones
				PairList exact = brute_force_all(floor, transformed(box, box_on_floor));
				std::sort(exact.begin(), exact.end());
				PairList cd = coldet_pairs(floor, w.world, box, compose(w.world, box_on_floor));
				EXPECT_EQ(compare_pairs(cd, exact, floor, Transform{}, box, box_on_floor).second, 0u)
					<< "COLDET found " << cd.size() << " of " << exact.size() << " contacts";
			}

			const Transform lifted{ z_rotation(spin), Vector3f(0.25f, 0.5f, 1.01f) };
			EXPECT_TRUE(check_scene(floor, w.world, box, compose(w.world, lifted)).empty()) << "lifted box";

			const Transform sunk{ z_rotation(spin), Vector3f(0.25f, 0.5f, 0.5f) };
			EXPECT_FALSE(check_scene(floor, w.world, box, compose(w.world, sunk)).empty()) << "sunk box";
		}
}

// Parallel axes: every axis-aligned relative rotation, with overlap, face, edge and
// corner contact, and a gap.  This is the degenerate case for separating axis tests.
TEST(scene, AxisAlignedRotations)
{
	const auto box = box_mesh(Vector3f(1, 1, 1));
	struct Offset { const char* name; Vector3f offset; bool hit; };
	const Offset offsets[] = {
		{ "overlap", { 1.5f, 0.25f, 0 }, true },
		{ "face contact", { 2, 0.5f, 0.25f }, true },
		{ "edge contact", { 2, 2, 0.5f }, true },
		{ "corner contact", { 2, 2, 2 }, true },
		{ "face gap", { 2.125f, 0, 0 }, false },
		{ "edge gap", { 2.125f, 2.125f, 0 }, false },
	};
	for (const Matrix3f& R : axis_rotations())
		for (const auto& o : offsets)
		{
			SCOPED_TRACE(std::string(o.name));
			PairList bf = check_scene(box, Transform{}, box, Transform{ R, o.offset });
			EXPECT_EQ(!bf.empty(), o.hit) << "rotation\n" << R;
		}
}

// Pure translations of a sphere through a box, on a grid of positions
TEST(scene, PureTranslations)
{
	const auto box = box_mesh(Vector3f(1, 0.5f, 2));
	const auto sphere = sphere_mesh(0.75f, 8, 16);
	size_t hits = 0;
	for (float x = -2.5f; x <= 2.5f; x += 0.5f)
		for (float y = -1.5f; y <= 1.5f; y += 0.75f)
			for (float z = -3.0f; z <= 3.0f; z += 1.0f)
			{
				SCOPED_TRACE("translation " + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z));
				hits += check_scene(box, Transform{}, sphere, translate(x, y, z)).empty() ? 0 : 1;
			}
	EXPECT_GT(hits, 0u);
}

// A model against itself at the same placement: every triangle meets itself and its neighbours
TEST(scene, SelfTest)
{
	const auto t = torus(2.0f, 0.5f, 0.05f);
	for (const Transform& T : { Transform{}, Transform{ random_rotation(), Vector3f(3, -1, 2) } })
	{
		PairList bf = check_scene(t, T, t, T);
		for (size_t i = 0; i < t.size(); ++i)
			EXPECT_TRUE(has_pair(bf, i, i)) << "triangle " << i << " against itself";
	}
}

// A tiny model against a large detailed one: one tree stays at its root while the other
// descends to its leaves
TEST(scene, TinyVsHuge)
{
	const float R = 100.0f;
	const auto big = sphere_mesh(R, 40, 80);
	const auto tiny = box_mesh(Vector3f(0.05f, 0.05f, 0.05f));
	const Matrix3f rot = random_rotation();
	for (const Vector3f& dir : { Vector3f(0, 0, 1), Vector3f(1, 2, 3).normalized(), Vector3f(-1, 0.5f, 0).normalized() })
	{
		// the sphere's facets lie up to R*(1-cos(pi/40)) ~ 0.3 inside the true radius
		struct Case { const char* name; float radius; bool hit; };
		for (const Case& c : { Case{ "outside", R + 1.0f, false }, Case{ "inside", R - 2.0f, false } })
		{
			SCOPED_TRACE(c.name);
			EXPECT_EQ(!check_scene(big, Transform{}, tiny, Transform{ rot, dir * c.radius }).empty(), c.hit);
		}
		// Straddle the surface: scan along the direction until the tiny box crosses it
		bool any = false;
		for (float r = R - 0.4f; r <= R + 0.1f; r += 0.025f)
			any |= !check_scene(big, Transform{}, tiny, Transform{ rot, dir * r }).empty();
		EXPECT_TRUE(any) << "tiny box crossing the surface";
	}
}

// Triangles of very different sizes in the same models
TEST(scene, MixedTriangleSizes)
{
	std::mt19937 rng(11);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f), log_size(-3.0f, 0.0f);
	auto soup = [&](size_t n) {
		std::vector<Triangle> res;
		for (size_t i = 0; i < n; ++i)
		{
			float s = std::pow(10.0f, log_size(rng));
			Vector3f c(unit(rng) * 3, unit(rng) * 3, unit(rng) * 3);
			auto v = [&]() { return Vector3f(c + Vector3f(unit(rng), unit(rng), unit(rng)) * s); };
			res.push_back(Triangle{ v(), v(), v() });
		}
		return res;
	};
	const auto a = soup(800), b = soup(800);
	size_t pairs = 0;
	for (int i = 0; i < 10; ++i)
		pairs += check_scene(a, Transform{ random_rotation(), Vector3f::Zero() }, b,
							 Transform{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) }).size();
	EXPECT_GT(pairs, 0u);
}

// Models far from the origin, where float precision is coarse
TEST(scene, FarFromOrigin)
{
	const auto t = torus(2.0f, 0.5f, 0.04f);
	std::mt19937 rng(3);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	size_t hits = 0;
	for (float distance : { 1e3f, 1e4f })
		for (int i = 0; i < 5; ++i)
		{
			const Vector3f base = Vector3f(1, -1, 1).normalized() * distance;
			const Transform Ta{ random_rotation(), base };
			const Transform Tb{ random_rotation(), base + Vector3f(unit(rng), unit(rng), unit(rng)) * 2.0f };
			hits += check_scene(t, Ta, t, Tb).empty() ? 0 : 1;
		}
	EXPECT_GT(hits, 0u);
}

// Both libraries test surfaces only: a model fully inside another does not collide
TEST(scene, ContainmentIsNotACollision)
{
	const auto big = box_mesh(Vector3f(5, 5, 5));
	const auto small = box_mesh(Vector3f(1, 1, 1));
	EXPECT_TRUE(check_scene(big, Transform{}, small, Transform{ random_rotation(), Vector3f(0.5f, -1, 1) }).empty());
	const auto sphere = sphere_mesh(1.0f, 8, 16);
	EXPECT_TRUE(check_scene(big, Transform{}, sphere, translate(1, 2, -1)).empty());
}

// Models with a single triangle or very few
TEST(scene, FewTriangles)
{
	const std::vector<Triangle> one{ Triangle{ {0, 0, 0}, {2, 0, 0}, {0, 2, 0} } };
	const std::vector<Triangle> vertical{ Triangle{ {0.5f, 0.25f, -1}, {0.5f, 0.25f, 1}, {0.5f, 1.25f, 0.5f} } };
	EXPECT_EQ(check_scene(one, Transform{}, vertical, Transform{}).size(), 1u);
	EXPECT_TRUE(check_scene(one, Transform{}, vertical, translate(0, 0, 3)).empty());

	const auto t = torus(2.0f, 0.5f, 0.05f);
	const Transform through{ x_rotation(PI / 2), Vector3f(2, 0, 0) }; // the triangle cuts the tube
	EXPECT_FALSE(check_scene(one, through, t, Transform{}).empty());
	EXPECT_FALSE(check_scene(t, Transform{}, one, through).empty());

	const auto two = box_mesh(Vector3f(1, 1, 1));
	const std::vector<Triangle> face(two.begin(), two.begin() + 2);
	EXPECT_FALSE(check_scene(face, Transform{}, two, translate(0.5f, 0.5f, 0)).empty());
}

// Models without triangles never collide
TEST(scene, EmptyModels)
{
	const std::vector<Triangle> empty;
	const auto box = box_mesh(Vector3f(1, 1, 1));
	EXPECT_TRUE(check_scene(empty, Transform{}, box, Transform{}).empty());
	EXPECT_TRUE(check_scene(box, Transform{}, empty, Transform{}).empty());
	EXPECT_TRUE(check_scene(empty, Transform{}, empty, Transform{}).empty());
}

// A lopsided model: a dense cluster of triangles plus one far away triangle
TEST(scene, ClusteredModel)
{
	std::mt19937 rng(8);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	std::vector<Triangle> cluster;
	for (int i = 0; i < 500; ++i)
	{
		auto v = [&]() { return Vector3f(Vector3f(unit(rng), unit(rng), unit(rng)) * 0.2f); };
		cluster.push_back(Triangle{ v(), v(), v() });
	}
	cluster.push_back(Triangle{ {50, 0, 0}, {52, 0, 0}, {50, 2, 0} });
	const std::vector<Triangle> probe{ Triangle{ {50.5f, 0.25f, -1}, {50.5f, 0.25f, 1}, {50.5f, 1.25f, 0.5f} } };
	PairList bf = check_scene(cluster, Transform{}, probe, Transform{});
	ASSERT_EQ(bf.size(), 1u);
	EXPECT_EQ(bf[0], TrianglePair(500, 0));
	EXPECT_FALSE(check_scene(cluster, Transform{}, box_mesh(Vector3f(0.1f, 0.1f, 0.1f)), Transform{}).empty());
}

// Degenerate triangles (segments and points) mixed into a mesh
TEST(scene, DegenerateTriangles)
{
	auto t = torus(2.0f, 0.5f, 0.05f);
	for (size_t i = 0; i < t.size(); i += 7)
		t[i].c = t[i].b;						// segment
	for (size_t i = 3; i < t.size(); i += 11)
		t[i].b = t[i].c = t[i].a;				// point
	for (int k = 0; k < 5; ++k)
		check_scene(t, Transform{}, t, Transform{ random_rotation(), Vector3f(0.5f, 0.25f, 0) });
}

// The same triangle twice in one model: both copies are reported
TEST(scene, DuplicateTriangles)
{
	auto box = box_mesh(Vector3f(1, 1, 1));
	box.push_back(box[0]);
	const std::vector<Triangle> probe{ Triangle{ {-0.5f, -0.5f, -2}, {0.5f, 0, 0}, {-0.5f, 0.5f, -2} } };
	PairList bf = check_scene(box, Transform{}, probe, Transform{});
	EXPECT_TRUE(has_pair(bf, 0, 0));
	EXPECT_TRUE(has_pair(bf, box.size() - 1, 0));
}
