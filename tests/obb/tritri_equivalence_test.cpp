// The C++ triangle test (src/tritri.cpp) against the frozen C original it was ported from
// (tests/reference/tritri_reference.c): results must be identical, for every entry point.
#include <gtest/gtest.h>
#include "test_helpers.h"

extern "C" int reference_tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3],
										   const float U0[3], const float U1[3], const float U2[3]);

namespace
{

struct Counts { size_t pairs = 0, hits = 0, mismatches = 0; std::string first; };

// Compares every entry point, in both argument orders, against the reference
void compare(const Triangle& x, const Triangle& y, Counts& c)
{
	const TriangleShape sx(x.a.data(), x.b.data(), x.c.data()), sy(y.a.data(), y.b.data(), y.c.data());
	for (int order = 0; order < 2; ++order)
	{
		const Triangle& v = order ? y : x;
		const Triangle& u = order ? x : y;
		const TriangleShape& sv = order ? sy : sx;
		const TriangleShape& su = order ? sx : sy;
		const bool expected = reference_tri_tri_intersect(v.a.data(), v.b.data(), v.c.data(), u.a.data(), u.b.data(), u.c.data()) != 0;
		const bool plain = tri_tri_intersect(v.a.data(), v.b.data(), v.c.data(), u.a.data(), u.b.data(), u.c.data());
		const bool both = tri_tri_intersect(v.a.data(), v.b.data(), v.c.data(), sv, u.a.data(), u.b.data(), u.c.data(), su);
		const bool first = tri_tri_intersect(v.a.data(), v.b.data(), v.c.data(), sv, u.a.data(), u.b.data(), u.c.data());
		++c.pairs;
		c.hits += expected ? 1 : 0;
		if (plain != expected || both != expected || first != expected)
		{
			if (c.mismatches++ == 0)
			{
				char buf[512];
				std::snprintf(buf, sizeof(buf), "reference %d, plain %d, shapes %d, first shape %d for "
					"{%.9g,%.9g,%.9g} {%.9g,%.9g,%.9g} {%.9g,%.9g,%.9g} / {%.9g,%.9g,%.9g} {%.9g,%.9g,%.9g} {%.9g,%.9g,%.9g}",
					int(expected), int(plain), int(both), int(first),
					v.a.x(), v.a.y(), v.a.z(), v.b.x(), v.b.y(), v.b.z(), v.c.x(), v.c.y(), v.c.z(),
					u.a.x(), u.a.y(), u.a.z(), u.b.x(), u.b.y(), u.b.z(), u.c.x(), u.c.y(), u.c.z());
				c.first = buf;
			}
		}
	}
}

void report(const char* name, const Counts& c)
{
	EXPECT_EQ(c.mismatches, 0u) << name << ": first mismatch: " << c.first;
	EXPECT_GT(c.hits, 0u) << name;
	EXPECT_LT(c.hits, c.pairs) << name;
	std::cout << name << ": " << c.pairs << " tests, " << c.hits << " intersecting" << std::endl;
}

Triangle place(const Triangle& t, const Matrix3f& R, const Vector3f& off, float scale)
{
	auto f = [&](const Vector3f& p) -> Vector3f { return (R * p + off) * scale; };
	return Triangle{ f(t.a), f(t.b), f(t.c) };
}

} // namespace

// General position, at several scales and distances from the origin
TEST(tritri_equivalence, RandomPairs)
{
	std::mt19937 rng(1001);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	Counts c;
	for (float scale : { 1e-3f, 1.0f, 1e3f })
		for (float distance : { 0.0f, 10.0f, 1000.0f })
			for (int i = 0; i < 30000; ++i)
			{
				const Vector3f center = Vector3f(unit(rng), unit(rng), unit(rng)) * distance;
				auto v = [&]() -> Vector3f { return (center + Vector3f(unit(rng), unit(rng), unit(rng))) * scale; };
				compare(Triangle{ v(), v(), v() }, Triangle{ v(), v(), v() }, c);
			}
	report("random pairs", c);
}

// Pairs whose vertices straddle the coplanarity tolerance: some vertices are zeroed and
// others not, which exercises every branch of the interval computation
TEST(tritri_equivalence, NearCoplanar)
{
	std::mt19937 rng(1002);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f), k(0.0f, 3.0f);
	Counts c;
	for (float distance : { 0.0f, 100.0f, 1000.0f })
		for (int i = 0; i < 60000; ++i)
		{
			// in the z=0 plane, then each vertex lifted by 0..3 times the size tolerance
			const float tol = 1e-6f * 2.0f;
			auto p = [&](bool lift) { return Vector3f(unit(rng), unit(rng), lift ? k(rng) * tol * (unit(rng) < 0 ? -1.0f : 1.0f) : 0.0f); };
			Triangle a{ p(false), p(false), p(false) }, b{ p(true), p(true), p(true) };
			const Matrix3f R = random_rotation();
			const Vector3f off = Vector3f(unit(rng), unit(rng), unit(rng)) * distance;
			compare(place(a, R, off, 1.0f), place(b, R, off, 1.0f), c);
		}
	report("near coplanar", c);
}

// Exactly coplanar pairs on a coarse grid, which makes collinear and touching edges common
TEST(tritri_equivalence, Coplanar)
{
	std::mt19937 rng(1003);
	std::uniform_int_distribution<int> q(-4, 4);
	const auto rotations = axis_rotations();
	Counts c;
	for (int i = 0; i < 100000; ++i)
	{
		auto p = [&]() { return Vector3f(float(q(rng)) / 4.0f, float(q(rng)) / 4.0f, 0.0f); };
		const Matrix3f& R = rotations[i % rotations.size()];
		Triangle a = place(Triangle{ p(), p(), p() }, R, Vector3f::Zero(), 1.0f);
		Triangle b = place(Triangle{ p(), p(), p() }, R, Vector3f::Zero(), 1.0f);
		compare(a, b, c);
		if (i % 4 == 0) // the same, rotated arbitrarily
		{
			const Matrix3f Q = random_rotation();
			compare(place(a, Q, Vector3f::Zero(), 1.0f), place(b, Q, Vector3f::Zero(), 1.0f), c);
		}
	}
	report("coplanar", c);
}

// Degenerate triangles: segments, points and slivers, against all kinds of triangles
TEST(tritri_equivalence, Degenerate)
{
	std::mt19937 rng(1004);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	std::uniform_int_distribution<int> kind(0, 4);
	Counts c;
	for (int i = 0; i < 100000; ++i)
	{
		auto v = [&]() { return Vector3f(unit(rng), unit(rng), unit(rng)); };
		auto make = [&](int k) {
			Triangle t{ v(), v(), v() };
			switch (k)
			{
			case 1: t.c = t.b; break;								// segment
			case 2: t.b = t.c = t.a; break;							// point
			case 3: t.c = t.a + (t.b - t.a) * 0.375f; break;		// collinear
			case 4: t.c = t.a + (t.b - t.a) * 0.5f + Vector3f(1e-7f, -1e-7f, 1e-7f); break; // sliver
			default: break;
			}
			return t;
		};
		compare(make(kind(rng)), make(kind(rng)), c);
	}
	report("degenerate", c);
}

// Mesh neighbours (shared edges and vertices), and two tori placed close to each other:
// the pairs the trees actually test, with many borderline contacts
TEST(tritri_equivalence, MeshPairs)
{
	const auto t = torus(2.0f, 0.5f, 0.04f);
	Counts c;
	for (size_t i = 0; i < t.size(); ++i)
		for (size_t j = i; j < std::min(t.size(), i + 64); ++j)
			compare(t[i], t[j], c);
	std::mt19937 rng(1005);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	for (int k = 0; k < 4; ++k)
	{
		const auto a = transformed(t, Transform{ random_rotation(), Vector3f::Zero() });
		const auto b = transformed(t, Transform{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) });
		for (const auto& x : a)
			for (const auto& y : b)
				compare(x, y, c);
	}
	report("mesh pairs", c);
}
