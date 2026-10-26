// Helpers shared by the OBB / COLDET test files
#pragma once

#include <algorithm>
#include <random>
#include <vector>
#include <obbtree.h>
#include <coldet.h>
#include <tritri.h>
#include <gtest/gtest.h>

using namespace obb;
using COLDET::tri_tri_intersect;
using COLDET::TriangleShape;



inline Matrix3f x_rotation(float angle_rad)
{
	Matrix3f R;
	float c = cos(angle_rad);
	float s = sin(angle_rad);
	R << 1, 0, 0,
		0, c, -s,
		0, s, c;
	return R;
}

inline Matrix3f y_rotation(float angle_rad)
{
	Matrix3f R;
	float c = cos(angle_rad);
	float s = sin(angle_rad);
	R << c, 0, s,
		0, 1, 0,
		-s, 0, c;
	return R;
}

inline Matrix3f z_rotation(float angle_rad)
{
	Matrix3f R;
	float c = cos(angle_rad);
	float s = sin(angle_rad);
	R << c, -s, 0,
		s, c, 0,
		0, 0, 1;
	return R;
}

inline Matrix3f random_rotation()
{
	static std::mt19937 rng(12345);
	static std::uniform_real_distribution<float> dist(0.0f, 2.0f * 3.14159265f);

	float ax = dist(rng);
	float ay = dist(rng);
	float az = dist(rng);

	return z_rotation(az) * y_rotation(ay) * x_rotation(ax);
}

using PairList = std::vector<obb::TrianglePair>;

// Triangles of size ~2*tri_scale, scattered in a 20x20x20 cube
inline std::vector<Triangle> random_soup(size_t n, float tri_scale, std::mt19937& rng)
{
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	// the return type must be explicit: an Eigen expression would refer to a destroyed temporary
	auto random_vector = [&](float scale) -> Vector3f { return Vector3f(unit(rng), unit(rng), unit(rng)) * scale; };
	std::vector<Triangle> soup(n);
	for (auto& tri : soup)
	{
		Transform T{ random_rotation(), random_vector(10.0f) };
		tri.a = T.apply(random_vector(tri_scale));
		tri.b = T.apply(random_vector(tri_scale));
		tri.c = T.apply(random_vector(tri_scale));
	}
	return soup;
}

// Torus around the z axis (ring radius R, tube radius r), with triangle edges
// about triangle_ratio times the torus diameter
inline std::vector<Triangle> torus(float R, float r, float triangle_ratio)
{
	const double PI = 3.14159265358979;
	const double edge = 2.0 * (R + r) * triangle_ratio;
	const size_t major = size_t(2.0 * PI * R / edge);
	const size_t minor = size_t(2.0 * PI * r / edge);
	auto P = [&](size_t i, size_t j) {
		double u = 2.0 * PI * double(i % major) / double(major);
		double v = 2.0 * PI * double(j % minor) / double(minor);
		double w = R + r * cos(v);
		return Vector3f(float(w * cos(u)), float(w * sin(u)), float(r * sin(v)));
	};
	std::vector<Triangle> res;
	res.reserve(2 * major * minor);
	for (size_t i = 0; i < major; ++i)
		for (size_t j = 0; j < minor; ++j)
		{
			res.push_back(Triangle{ P(i, j), P(i + 1, j), P(i + 1, j + 1) });
			res.push_back(Triangle{ P(i, j), P(i + 1, j + 1), P(i, j + 1) });
		}
	return res;
}

inline CollisionModel3D* make_coldet_model(const std::vector<Triangle>& soup)
{
	CollisionModel3D* model = newCollisionModel3D(false);
	model->setTriangleNumber(int(soup.size()));
	for (const auto& tri : soup)
		model->addTriangle(tri.a.data(), tri.b.data(), tri.c.data());
	model->finalize();
	return model;
}

inline OBBTree make_obb_tree(const std::vector<Triangle>& soup, size_t max_triangles_per_leaf = 4)
{
	OBBTree tree;
	for (const auto& tri : soup)
		tree.add_triangle(tri);
	tree.build(max_triangles_per_leaf);
	return tree;
}

// COLDET uses row vectors (p' = p * M), so the rotation is stored transposed
inline void to_coldet_matrix(const Transform& T, float m[16])
{
	for (int r = 0; r < 3; ++r)
	{
		for (int c = 0; c < 3; ++c)
			m[r * 4 + c] = T.rotation(c, r);
		m[r * 4 + 3] = 0.0f;
		m[12 + r] = T.center[r];
	}
	m[15] = 1.0f;
}

inline std::vector<Triangle> transformed(const std::vector<Triangle>& soup, const Transform& T)
{
	std::vector<Triangle> res(soup.size());
	for (size_t i = 0; i < soup.size(); ++i)
		res[i] = Triangle{ T.apply(soup[i].a), T.apply(soup[i].b), T.apply(soup[i].c) };
	return res;
}

inline PairList brute_force_all(const std::vector<Triangle>& a, const std::vector<Triangle>& b)
{
	PairList pairs;
	std::vector<TriangleShape> sa, sb;
	for (const auto& t : a) sa.emplace_back(t.a.data(), t.b.data(), t.c.data());
	for (const auto& t : b) sb.emplace_back(t.a.data(), t.b.data(), t.c.data());
	for (size_t i = 0; i < a.size(); ++i)
		for (size_t j = 0; j < b.size(); ++j)
			if (tri_tri_intersect(a[i].a.data(), a[i].b.data(), a[i].c.data(), sa[i],
								  b[j].a.data(), b[j].b.data(), b[j].c.data(), sb[j]))
				pairs.emplace_back(i, j);
	return pairs;
}

// Separation of two triangles in double precision, by the separating axis test over all 17
// candidate axes: > 0 is the gap between disjoint triangles, <= 0 minus the smallest overlap
// depth.  Coplanar triangles have a separation of 0 along their normal.
inline double separation(const Triangle& ta, const Triangle& tb)
{
	using V3d = Eigen::Vector3d;
	const V3d a[3] = { ta.a.cast<double>(), ta.b.cast<double>(), ta.c.cast<double>() };
	const V3d b[3] = { tb.a.cast<double>(), tb.b.cast<double>(), tb.c.cast<double>() };
	const V3d ea[3] = { a[1] - a[0], a[2] - a[1], a[0] - a[2] }, eb[3] = { b[1] - b[0], b[2] - b[1], b[0] - b[2] };
	const V3d na = ea[0].cross(ea[1]), nb = eb[0].cross(eb[1]);
	std::vector<V3d> axes = { na, nb };
	for (auto& x : ea) for (auto& y : eb) axes.push_back(x.cross(y));
	for (auto& x : ea) axes.push_back(na.cross(x));
	for (auto& y : eb) axes.push_back(nb.cross(y));
	double best = -1e300;
	for (V3d axis : axes)
	{
		if (axis.squaredNorm() < 1e-30) continue;
		axis.normalize();
		double amin = 1e300, amax = -1e300, bmin = 1e300, bmax = -1e300;
		for (auto& p : a) { double d = axis.dot(p); amin = std::min(amin, d); amax = std::max(amax, d); }
		for (auto& p : b) { double d = axis.dot(p); bmin = std::min(bmin, d); bmax = std::max(bmax, d); }
		best = std::max(best, std::max(bmin - amax, amin - bmax));
	}
	return best;
}

// COLDET tests triangles in model A's frame, the OBB tree in world space, so rounding can
// flip pairs that barely touch.  Returns how many pairs differ between the two sorted lists,
// and how many of those are not borderline.  A pair is borderline when:
//  - its result changes when both triangles are shrunk or grown by 1% about their centroids
//    (far from the origin, the triangle test's tolerance grows with the coordinates, so a
//    smaller perturbation may not leave it), or
//  - the triangles touch to within 'precision' (gap or overlap depth): with a precision of a
//    few ulps of the coordinates, contacts the float input can't resolve, like faces at rest
//    on each other far from the origin.
inline std::pair<size_t, size_t> compare_pairs(const PairList& p1, const PairList& p2,
											   const std::vector<Triangle>& soup_a, const Transform& Ta,
											   const std::vector<Triangle>& soup_b, const Transform& Tb,
											   double precision = 0.0)
{
	auto scaled = [](const Triangle& t, const Transform& T, float f) {
		Vector3f a = T.apply(t.a), b = T.apply(t.b), c = T.apply(t.c), m = (a + b + c) / 3.0f;
		return Triangle{ m + (a - m) * f, m + (b - m) * f, m + (c - m) * f };
	};
	auto hits = [&](const obb::TrianglePair& pr, float f) {
		Triangle x = scaled(soup_a[pr.first], Ta, f), y = scaled(soup_b[pr.second], Tb, f);
		return tri_tri_intersect(x.a.data(), x.b.data(), x.c.data(), y.a.data(), y.b.data(), y.c.data()) != 0;
	};
	PairList diff;
	std::set_symmetric_difference(p1.begin(), p1.end(), p2.begin(), p2.end(), std::back_inserter(diff));
	size_t unexplained = 0;
	for (const auto& pr : diff)
	{
		if (!hits(pr, 0.99f) && hits(pr, 1.01f))
			continue;
		if (precision > 0.0 && std::abs(separation(scaled(soup_a[pr.first], Ta, 1.0f), scaled(soup_b[pr.second], Tb, 1.0f))) <= precision)
			continue;
		++unexplained;
	}
	return { diff.size(), unexplained };
}


// ---- Mesh builders ----

// Closed box centered at the origin, 12 triangles with outward winding
inline std::vector<Triangle> box_mesh(const Vector3f& half)
{
	auto P = [&](int sx, int sy, int sz) { return Vector3f(sx * half.x(), sy * half.y(), sz * half.z()); };
	std::vector<Triangle> res;
	auto quad = [&](const Vector3f& a, const Vector3f& b, const Vector3f& c, const Vector3f& d) {
		res.push_back(Triangle{ a, b, c });
		res.push_back(Triangle{ a, c, d });
	};
	quad(P(-1, -1, -1), P(-1, 1, -1), P(1, 1, -1), P(1, -1, -1)); // -z
	quad(P(-1, -1, 1), P(1, -1, 1), P(1, 1, 1), P(-1, 1, 1));     // +z
	quad(P(-1, -1, -1), P(1, -1, -1), P(1, -1, 1), P(-1, -1, 1)); // -y
	quad(P(-1, 1, -1), P(-1, 1, 1), P(1, 1, 1), P(1, 1, -1));     // +y
	quad(P(-1, -1, -1), P(-1, -1, 1), P(-1, 1, 1), P(-1, 1, -1)); // -x
	quad(P(1, -1, -1), P(1, 1, -1), P(1, 1, 1), P(1, -1, 1));     // +x
	return res;
}

// Flat n x n grid of quads in the z=0 plane, covering [-half, half]^2
inline std::vector<Triangle> grid_mesh(size_t n, float half)
{
	std::vector<Triangle> res;
	auto P = [&](size_t i, size_t j) {
		return Vector3f(-half + 2.0f * half * float(i) / float(n), -half + 2.0f * half * float(j) / float(n), 0.0f);
	};
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
		{
			res.push_back(Triangle{ P(i, j), P(i + 1, j), P(i + 1, j + 1) });
			res.push_back(Triangle{ P(i, j), P(i + 1, j + 1), P(i, j + 1) });
		}
	return res;
}

// UV sphere centered at the origin
inline std::vector<Triangle> sphere_mesh(float radius, size_t rings, size_t segments)
{
	const double PI = 3.14159265358979;
	auto P = [&](size_t r, size_t s) {
		double theta = PI * double(r) / double(rings);
		double phi = 2.0 * PI * double(s % segments) / double(segments);
		return Vector3f(float(radius * sin(theta) * cos(phi)), float(radius * sin(theta) * sin(phi)), float(radius * cos(theta)));
	};
	std::vector<Triangle> res;
	for (size_t r = 0; r < rings; ++r)
		for (size_t s = 0; s < segments; ++s)
		{
			if (r > 0) res.push_back(Triangle{ P(r, s), P(r + 1, s), P(r, s + 1) });
			if (r + 1 < rings) res.push_back(Triangle{ P(r, s + 1), P(r + 1, s), P(r + 1, s + 1) });
		}
	return res;
}

inline std::vector<Triangle> concat(std::vector<Triangle> a, const std::vector<Triangle>& b)
{
	a.insert(a.end(), b.begin(), b.end());
	return a;
}

// All 24 rotations that map the coordinate axes onto themselves
inline std::vector<Matrix3f> axis_rotations()
{
	std::vector<Matrix3f> res;
	int perm[6][3] = { {0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0} };
	for (auto& p : perm)
		for (int signs = 0; signs < 8; ++signs)
		{
			Matrix3f R = Matrix3f::Zero();
			for (int i = 0; i < 3; ++i)
				R(i, p[i]) = (signs >> i) & 1 ? -1.0f : 1.0f;
			if (R.determinant() > 0.0f)
				res.push_back(R);
		}
	return res;
}

// ---- Scene checker ----

// All colliding pairs found by COLDET, sorted.  Also checks that its first-hit query agrees
// and reports one of them.
inline PairList coldet_pairs(const std::vector<Triangle>& a, const Transform& Ta,
							 const std::vector<Triangle>& b, const Transform& Tb, bool is_static = false)
{
	CollisionModel3D* ca = newCollisionModel3D(is_static);
	CollisionModel3D* cb = newCollisionModel3D(is_static);
	for (const auto& t : a) ca->addTriangle(t.a.data(), t.b.data(), t.c.data());
	for (const auto& t : b) cb->addTriangle(t.a.data(), t.b.data(), t.c.data());
	ca->finalize();
	cb->finalize();
	float ma[16], mb[16];
	to_coldet_matrix(Ta, ma);
	to_coldet_matrix(Tb, mb);
	ca->setTransform(ma);
	cb->setTransform(mb);
	std::vector<std::pair<int, int>> raw;
	ca->collisionAll(cb, raw);
	PairList cd(raw.begin(), raw.end());
	std::sort(cd.begin(), cd.end());
	bool hit = ca->collision(cb);
	EXPECT_EQ(hit, !cd.empty());
	if (hit)
	{
		int i1, i2;
		ca->getCollidingTriangles(i1, i2);
		EXPECT_TRUE(std::binary_search(cd.begin(), cd.end(), TrianglePair(i1, i2)));
	}
	delete ca;
	delete cb;
	return cd;
}

// Runs both libraries on one placement of two models, and checks them against brute force:
//  - the OBB tree (world space, like brute force) must find exactly the brute force pairs,
//    for several leaf sizes
//  - COLDET must find the same pairs, apart from borderline ones (see compare_pairs)
//  - each first-hit query must agree with its find-all result, and report a real pair
// Returns the brute force pairs, for scene specific expectations.
inline PairList check_scene(const std::vector<Triangle>& a, const Transform& Ta,
							const std::vector<Triangle>& b, const Transform& Tb)
{
	PairList bf = brute_force_all(transformed(a, Ta), transformed(b, Tb));
	std::sort(bf.begin(), bf.end());

	for (size_t leaf : { 1, 4 })
	{
		SCOPED_TRACE("OBB tree, leaf size " + std::to_string(leaf));
		OBBTree ta = make_obb_tree(a, leaf), tb = make_obb_tree(b, leaf);
		EXPECT_TRUE(ta.verify());
		EXPECT_TRUE(tb.verify());
		ta.set_dynamic_transform(Ta);
		tb.set_dynamic_transform(Tb);
		IntersectionResult all = ta.find_all_intersections(tb);
		std::sort(all.pairs.begin(), all.pairs.end());
		EXPECT_EQ(all.pairs, bf);
		IntersectionResult first = ta.check_intersection(tb);
		EXPECT_EQ(first.hit, !bf.empty());
		if (first.hit)
			EXPECT_TRUE(std::binary_search(bf.begin(), bf.end(), TrianglePair(first.triangle_index1, first.triangle_index2)));
		if (leaf == 1)
		{
			// the tree's own brute force query
			IntersectionResult own = ta.bf_check_intersection(tb);
			EXPECT_EQ(own.hit, !bf.empty());
			if (own.hit)
				EXPECT_TRUE(std::binary_search(bf.begin(), bf.end(), TrianglePair(own.triangle_index1, own.triangle_index2)));
		}
	}

	// COLDET works in model A's frame, from float placement matrices: far from the origin,
	// those place the models a few ulps of their translation off, so allow contacts within that
	float coord = 0.0f;
	for (const auto& t : transformed(a, Ta)) coord = std::max({ coord, t.a.cwiseAbs().maxCoeff(), t.b.cwiseAbs().maxCoeff(), t.c.cwiseAbs().maxCoeff() });
	for (const auto& t : transformed(b, Tb)) coord = std::max({ coord, t.a.cwiseAbs().maxCoeff(), t.b.cwiseAbs().maxCoeff(), t.c.cwiseAbs().maxCoeff() });
	const double precision = 64.0 * std::numeric_limits<float>::epsilon() * coord;
	for (bool is_static : { false, true })
	{
		SCOPED_TRACE(is_static ? "COLDET static model" : "COLDET dynamic model");
		PairList cd = coldet_pairs(a, Ta, b, Tb, is_static);
		EXPECT_EQ(compare_pairs(cd, bf, a, Ta, b, Tb, precision).second, 0u);
	}
	return bf;
}
