// Tests of tri_tri_intersect against known answers and an independent reference
#include <array>
#include <functional>
#include <cmath>
#include <sstream>
#include <gtest/gtest.h>
#include "test_helpers.h"

namespace
{

struct TriCase
{
	const char* name;
	Triangle t2;	// tested against the base triangle (0,0,0) (2,0,0) (0,2,0)
	bool expected;
	bool robust;	// clear margins: also valid under arbitrary rotations
};

const Triangle BASE{ {0, 0, 0}, {2, 0, 0}, {0, 2, 0} };

// All coordinates are multiples of 1/8, so they stay exact under the axis
// permutations, sign flips, power of two scales and offsets applied below.
const std::vector<TriCase> CASES = {
	// General position
	{ "crossing interior",            { {0.5f, 0.25f, -1}, {0.5f, 0.25f, 1}, {0.5f, 1.25f, 0.5f} }, true, true },
	{ "crosses plane, outside",       { {3.5f, 0.25f, -1}, {3.5f, 0.25f, 1}, {3.5f, 1.25f, 0.5f} }, false, true },
	{ "parallel above",               { {0, 0, 0.5f}, {2, 0, 0.5f}, {0, 2, 0.5f} }, false, true },
	{ "pierces near hypotenuse",      { {1.5f, 0.25f, -1}, {1.5f, 0.25f, 1}, {1.5f, 1.25f, 0.5f} }, true, true },
	{ "misses beyond hypotenuse",     { {1.5f, 0.625f, -1}, {1.5f, 0.625f, 1}, {1.5f, 1.625f, 0.5f} }, false, true },
	{ "big triangle through small",   { {-8, -8, -1}, {8, -8, -1}, {0, 8, 0.875f} }, true, true },
	// Touching
	{ "vertex on face",               { {0.5f, 0.5f, 0}, {1, 0.5f, 1}, {0.5f, 1, 1} }, true, false },
	{ "vertex just above face",       { {0.5f, 0.5f, 0.125f}, {1, 0.5f, 1}, {0.5f, 1, 1} }, false, true },
	{ "edge touches edge",            { {1, 0, -1}, {1, 0, 1}, {1, -1, 0} }, true, false },
	{ "edge just misses edge",        { {1, -0.125f, -1}, {1, -0.125f, 1}, {1, -1.125f, 0} }, false, true },
	{ "edge lies on face",            { {0.25f, 0.25f, 0}, {1, 0.25f, 0}, {0.5f, 0.5f, 1} }, true, false },
	{ "vertex touches vertex",        { {2, 0, 0}, {3, 0, 1}, {3, 1, 1} }, true, false },
	// Mesh neighbours
	{ "shared edge, folded",          { {2, 0, 0}, {0, 2, 0}, {2, 2, 1} }, true, false },
	{ "shared edge, coplanar",        { {2, 0, 0}, {0, 2, 0}, {2, 2, 0} }, true, false },
	{ "shared vertex",                { {0, 0, 0}, {-1, 0, 1}, {0, -1, 1} }, true, false },
	// Coplanar
	{ "coplanar overlapping",         { {0.5f, 0.5f, 0}, {3, 0.5f, 0}, {0.5f, 3, 0} }, true, true },
	{ "coplanar contained",           { {0.25f, 0.25f, 0}, {0.75f, 0.25f, 0}, {0.25f, 0.75f, 0} }, true, true },
	{ "coplanar containing",          { {-1, -1, 0}, {5, -1, 0}, {-1, 5, 0} }, true, true },
	{ "coplanar disjoint",            { {3, 3, 0}, {4, 3, 0}, {3, 4, 0} }, false, true },
	{ "coplanar gap at hypotenuse",   { {1.125f, 1.125f, 0}, {2, 1.125f, 0}, {1.125f, 2, 0} }, false, true },
	{ "coplanar vertex contact",      { {2, 0, 0}, {3, 0, 0}, {3, 1, 0} }, true, false },
	{ "coplanar collinear edges",     { {1, 0, 0}, {3, 0, 0}, {2, -1, 0} }, true, false },
	{ "coplanar star (no vertex in)", { {0.5f, -0.5f, 0}, {1, 2, 0}, {-0.5f, 0.75f, 0} }, true, true },
	// Near coplanar: distances within / beyond the relative tolerance
	{ "near coplanar, within eps",    { {0.5f, 0.5f, 1e-7f}, {3, 0.5f, 1e-7f}, {0.5f, 3, 1e-7f} }, true, false },
	{ "near coplanar, beyond eps",    { {0.5f, 0.5f, 0.015625f}, {3, 0.5f, 0.015625f}, {0.5f, 3, 0.015625f} }, false, false },
	// Degenerate second triangle
	{ "segment piercing face",        { {0.5f, 0.5f, -1}, {0.5f, 0.5f, 1}, {0.5f, 0.5f, 0.25f} }, true, false },
	{ "segment missing face",         { {3, 3, -1}, {3, 3, 1}, {3, 3, 0.25f} }, false, false },
	{ "segment beside face, shadow on it", { {0.5f, 0.5f, 1}, {3, 3, -1}, {3, 3, -1} }, false, false },
	{ "collinear points beside face", { {0.5f, 0.5f, 1}, {1.75f, 1.75f, 0}, {3, 3, -1} }, false, false },
	{ "segment tip on face",          { {0.5f, 0.5f, 0}, {1, 1, 1}, {1, 1, 1} }, true, false },
	{ "segment along an edge",        { {0.5f, 0, 0}, {1.5f, 0, 0}, {1.5f, 0, 0} }, true, false },
	{ "point on face",                { {0.5f, 0.5f, 0}, {0.5f, 0.5f, 0}, {0.5f, 0.5f, 0} }, true, false },
	{ "point off face",               { {0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f} }, false, false },
	{ "segment across face, in plane",{ {-1, 0.5f, 0}, {3, 0.5f, 0}, {1, 0.5f, 0} }, true, false },
	{ "segment in plane, outside",    { {-1, 3, 0}, {3, 3, 0}, {1, 3, 0} }, false, false },
};

bool tri_tri(const Triangle& a, const Triangle& b)
{
	return tri_tri_intersect(a.a.data(), a.b.data(), a.c.data(), b.a.data(), b.b.data(), b.c.data()) != 0;
}

Triangle map(const Triangle& t, const std::function<Vector3f(const Vector3f&)>& f)
{
	return Triangle{ f(t.a), f(t.b), f(t.c) };
}

// The 6 vertex orders of a triangle
std::array<Triangle, 6> orders(const Triangle& t)
{
	return { Triangle{ t.a, t.b, t.c }, Triangle{ t.b, t.c, t.a }, Triangle{ t.c, t.a, t.b },
			 Triangle{ t.a, t.c, t.b }, Triangle{ t.c, t.b, t.a }, Triangle{ t.b, t.a, t.c } };
}

// Checks one placement under all vertex orders and both argument orders.
// Returns a description of the first failure, or an empty string.
std::string check_all_orders(const Triangle& t1, const Triangle& t2, bool expected)
{
	for (const auto& a : orders(t1))
		for (const auto& b : orders(t2))
		{
			if (tri_tri(a, b) != expected) return "tri_tri(base, t2)";
			if (tri_tri(b, a) != expected) return "tri_tri(t2, base)";
		}
	return "";
}

// 48 signed axis permutations (rotations and reflections), exact in floating point
std::vector<Matrix3f> signed_permutations()
{
	std::vector<Matrix3f> res;
	int perm[6][3] = { {0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0} };
	for (auto& p : perm)
		for (int signs = 0; signs < 8; ++signs)
		{
			Matrix3f R = Matrix3f::Zero();
			for (int i = 0; i < 3; ++i)
				R(i, p[i]) = (signs >> i) & 1 ? -1.0f : 1.0f;
			res.push_back(R);
		}
	return res;
}

// gtest prints parameters by name instead of their bytes
void PrintTo(const TriCase& c, std::ostream* os) { *os << c.name; }

std::vector<TriCase> robust_cases()
{
	std::vector<TriCase> res;
	for (const auto& c : CASES) if (c.robust) res.push_back(c);
	return res;
}

std::string case_name(const ::testing::TestParamInfo<TriCase>& info)
{
	std::string n = info.param.name;
	for (char& ch : n) if (!isalnum((unsigned char)ch)) ch = '_';
	return n;
}

class TriTriCase : public ::testing::TestWithParam<TriCase> {};
// Touching and degenerate cases only hold under exact symmetries, not arbitrary rotations
class TriTriRobustCase : public ::testing::TestWithParam<TriCase> {};

} // namespace

// Every case under exact symmetries: 48 axis permutations / reflections,
// power of two scales from 2^-10 to 2^10, and offsets from the origin.
TEST_P(TriTriCase, ExactSymmetries)
{
	const TriCase& c = GetParam();
	const Vector3f offsets[] = { {0, 0, 0}, {64, -32, 16}, {-1024, 512, 2048} };
	size_t failures = 0, total = 0;
	std::string first;
	for (const Matrix3f& R : signed_permutations())
		for (int exp : { -10, -5, 0, 5, 10 })
			for (const Vector3f& off : offsets)
			{
				// keep a 1/8 resolution exact: |coordinate / scale| must stay below 2^20
				const float scale = std::ldexp(1.0f, exp);
				auto f = [&](const Vector3f& p) { return Vector3f(R * p * scale + off * scale); };
				std::string err = check_all_orders(map(BASE, f), map(c.t2, f), c.expected);
				++total;
				if (!err.empty() && failures++ == 0)
				{
					std::ostringstream os;
					os << err << " with scale 2^" << exp << ", offset " << off.transpose() << ", R=\n" << R;
					first = os.str();
				}
			}
	EXPECT_EQ(failures, 0u) << c.name << ": expected " << c.expected << ", failed " << failures << "/" << total
							<< " placements. First: " << first;
}

// Cases with clear margins under arbitrary rotations and offsets
TEST_P(TriTriRobustCase, RandomRotations)
{
	const TriCase& c = GetParam();
	std::mt19937 rng(77);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	for (float distance : { 0.0f, 10.0f, 100.0f, 1000.0f })
		for (float scale : { 1e-3f, 1.0f, 1e3f })
		{
			size_t failures = 0;
			for (int i = 0; i < 200; ++i)
			{
				Matrix3f R = random_rotation();
				Vector3f off = Vector3f(unit(rng), unit(rng), unit(rng)) * distance;
				auto f = [&](const Vector3f& p) { return Vector3f((R * p + off) * scale); };
				if (!check_all_orders(map(BASE, f), map(c.t2, f), c.expected).empty())
					++failures;
			}
			EXPECT_EQ(failures, 0u) << c.name << ": distance from origin " << distance
									<< " triangle sizes, scale " << scale;
		}
}

INSTANTIATE_TEST_SUITE_P(tritri, TriTriCase, ::testing::ValuesIn(CASES), case_name);
INSTANTIATE_TEST_SUITE_P(tritri, TriTriRobustCase, ::testing::ValuesIn(robust_cases()), case_name);

namespace
{

using V3d = Eigen::Vector3d;
struct TriD { V3d v[3]; };

// Independent reference: separating axis test in double precision, over all 17 candidate
// axes (2 normals, 9 edge-edge cross products, 6 in-plane edge normals for coplanar cases).
// Returns the largest separation found (> 0: disjoint, <= 0: intersecting or touching),
// normalized by the triangles' size.
double sat_separation(const TriD& a, const TriD& b, bool in_plane_axes_only = false)
{
	auto edges = [](const TriD& t) { return std::array<V3d, 3>{ t.v[1] - t.v[0], t.v[2] - t.v[1], t.v[0] - t.v[2] }; };
	auto ea = edges(a), eb = edges(b);
	V3d na = ea[0].cross(ea[1]), nb = eb[0].cross(eb[1]);
	// For coplanar triangles, the normal and edge-edge axes all give a separation of 0
	std::vector<V3d> axes;
	if (!in_plane_axes_only)
	{
		axes = { na, nb };
		for (auto& x : ea) for (auto& y : eb) axes.push_back(x.cross(y));
	}
	for (auto& x : ea) axes.push_back(na.cross(x));
	for (auto& y : eb) axes.push_back(nb.cross(y));
	double size = 0;
	for (auto& e : ea) size = std::max(size, e.norm());
	for (auto& e : eb) size = std::max(size, e.norm());
	double best = -1e300;
	for (V3d axis : axes)
	{
		double len = axis.norm();
		if (len < 1e-9 * size * size) continue; // degenerate axis
		axis /= len;
		double amin = 1e300, amax = -1e300, bmin = 1e300, bmax = -1e300;
		for (auto& p : a.v) { double d = axis.dot(p); amin = std::min(amin, d); amax = std::max(amax, d); }
		for (auto& p : b.v) { double d = axis.dot(p); bmin = std::min(bmin, d); bmax = std::max(bmax, d); }
		best = std::max(best, std::max(bmin - amax, amin - bmax));
	}
	return best / size;
}

} // namespace

// Random triangle pairs at several scales and distances from the origin, compared to the
// double precision reference.  Pairs within 1e-4 (relative) of touching are skipped.
TEST(tritri, MatchesReference)
{
	std::mt19937 rng(99);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	size_t compared = 0, hits = 0;
	for (float scale : { 1e-3f, 1.0f, 1e3f })
		for (float distance : { 0.0f, 100.0f })
		{
			size_t failures = 0;
			for (int i = 0; i < 20000; ++i)
			{
				Vector3f center = Vector3f(unit(rng), unit(rng), unit(rng)) * distance * scale;
				auto vertex = [&]() { return Vector3f(center + Vector3f(unit(rng), unit(rng), unit(rng)) * scale); };
				Triangle a{ vertex(), vertex(), vertex() }, b{ vertex(), vertex(), vertex() };
				TriD da{ { a.a.cast<double>(), a.b.cast<double>(), a.c.cast<double>() } };
				TriD db{ { b.a.cast<double>(), b.b.cast<double>(), b.c.cast<double>() } };
				double sep = sat_separation(da, db);
				if (std::abs(sep) < 1e-4) continue;
				bool expected = sep <= 0;
				++compared;
				hits += expected ? 1 : 0;
				if (tri_tri(a, b) != expected) ++failures;
			}
			EXPECT_EQ(failures, 0u) << "scale " << scale << ", distance " << distance;
		}
	std::cout << "Compared " << compared << " random pairs, " << hits << " intersecting" << std::endl;
}

// Regression: coplanar triangles 0.75 apart, with an edge of each on one line (a floor
// grid diagonal lined up with a box face's diagonal, after a rotation).  The parallel edges'
// crossing test used to pass on rounding noise.
TEST(tritri, CoplanarCollinearEdges)
{
	const Triangle a{ {-1.89438319f, -0.434346676f, -1.07463372f}, {-1.33999681f, -1.09678233f, -1.9781462f}, {-0.242805243f, -0.979608953f, -1.39082944f} };
	const Triangle b{ {-0.471508384f, 0.150591969f, 0.407180786f}, {2.17101622f, -0.721827686f, -0.098732233f}, {0.41550979f, -0.909305036f, -1.03843927f} };
	EXPECT_EQ(check_all_orders(a, b, false), "") << "collinear edges, triangles apart";
	const Triangle base{ {0, 0, 0}, {2, 0, 0}, {0, 2, 0} };
	EXPECT_EQ(check_all_orders(base, Triangle{ {2.5f, 0, 0}, {3, 0, 0}, {2.5f, -1, 0} }, false), "") << "collinear edges with a gap";
	EXPECT_EQ(check_all_orders(base, Triangle{ {2, 0, 0}, {3, 0, 0}, {2.5f, -1, 0} }, true), "") << "collinear edges meeting at a vertex";
	EXPECT_EQ(check_all_orders(base, Triangle{ {0.5f, 0, 0}, {1.5f, 0, 0}, {1, -1, 0} }, true), "") << "collinear edges, one within the other";
}

// Coplanar random pairs: exercises the coplanar code path against the reference
TEST(tritri, CoplanarMatchesReference)
{
	std::mt19937 rng(5);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	size_t compared = 0, hits = 0, failures = 0;
	const auto rotations = axis_rotations();
	for (int i = 0; i < 20000; ++i)
	{
		// exact in-plane coordinates (multiples of 1/64), then an exact axis rotation
		auto q = [&]() { return std::round(unit(rng) * 64.0f) / 64.0f; };
		Triangle a{ {q(), q(), 0}, {q(), q(), 0}, {q(), q(), 0} };
		Triangle b{ {q(), q(), 0}, {q(), q(), 0}, {q(), q(), 0} };
		const Matrix3f& R = rotations[i % rotations.size()];
		a = map(a, [&](const Vector3f& p) { return Vector3f(R * p); });
		b = map(b, [&](const Vector3f& p) { return Vector3f(R * p); });
		TriD da{ { a.a.cast<double>(), a.b.cast<double>(), a.c.cast<double>() } };
		TriD db{ { b.a.cast<double>(), b.b.cast<double>(), b.c.cast<double>() } };
		auto area = [](const TriD& t) { return (t.v[1] - t.v[0]).cross(t.v[2] - t.v[0]).norm(); };
		if (area(da) < 1e-3 || area(db) < 1e-3) continue; // skip slivers
		double sep = sat_separation(da, db, true);
		if (std::abs(sep) < 1e-4) continue;
		bool expected = sep <= 0;
		++compared;
		hits += expected ? 1 : 0;
		if (tri_tri(a, b) != expected) ++failures;
	}
	EXPECT_EQ(failures, 0u);
	std::cout << "Compared " << compared << " coplanar pairs, " << hits << " intersecting" << std::endl;
}
