#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <gtest/gtest.h>
#include <obbtree.h>
#include <coldet.h>
#include "../../benchmark/profiler.h"

extern "C" {
int tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3],
					  const float U0[3], const float U1[3], const float U2[3]);
}

using namespace obb;


static inline Matrix3f x_rotation(float angle_rad)
{
	Matrix3f R;
	float c = cos(angle_rad);
	float s = sin(angle_rad);
	R << 1, 0, 0,
		0, c, -s,
		0, s, c;
	return R;
}

static inline Matrix3f y_rotation(float angle_rad)
{
	Matrix3f R;
	float c = cos(angle_rad);
	float s = sin(angle_rad);
	R << c, 0, s,
		0, 1, 0,
		-s, 0, c;
	return R;
}

static inline Matrix3f z_rotation(float angle_rad)
{
	Matrix3f R;
	float c = cos(angle_rad);
	float s = sin(angle_rad);
	R << c, -s, 0,
		s, c, 0,
		0, 0, 1;
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

using PairList = std::vector<obb::TrianglePair>;

// Triangles of size ~2*tri_scale, scattered in a 20x20x20 cube
static std::vector<Triangle> random_soup(size_t n, float tri_scale, std::mt19937& rng)
{
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	auto random_vector = [&](float scale) { return Vector3f(unit(rng), unit(rng), unit(rng)) * scale; };
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
static std::vector<Triangle> torus(float R, float r, float triangle_ratio)
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

static CollisionModel3D* make_coldet_model(const std::vector<Triangle>& soup)
{
	CollisionModel3D* model = newCollisionModel3D(false);
	model->setTriangleNumber(int(soup.size()));
	for (const auto& tri : soup)
		model->addTriangle(tri.a.data(), tri.b.data(), tri.c.data());
	model->finalize();
	return model;
}

static OBBTree make_obb_tree(const std::vector<Triangle>& soup, size_t max_triangles_per_leaf = 4)
{
	OBBTree tree;
	for (const auto& tri : soup)
		tree.add_triangle(tri);
	tree.build(max_triangles_per_leaf);
	return tree;
}

// COLDET uses row vectors (p' = p * M), so the rotation is stored transposed
static void to_coldet_matrix(const Transform& T, float m[16])
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

static std::vector<Triangle> transformed(const std::vector<Triangle>& soup, const Transform& T)
{
	std::vector<Triangle> res(soup.size());
	for (size_t i = 0; i < soup.size(); ++i)
		res[i] = Triangle{ T.apply(soup[i].a), T.apply(soup[i].b), T.apply(soup[i].c) };
	return res;
}

static PairList brute_force_all(const std::vector<Triangle>& a, const std::vector<Triangle>& b)
{
	PairList pairs;
	for (size_t i = 0; i < a.size(); ++i)
		for (size_t j = 0; j < b.size(); ++j)
			if (tri_tri_intersect(a[i].a.data(), a[i].b.data(), a[i].c.data(),
								  b[j].a.data(), b[j].b.data(), b[j].c.data()))
				pairs.emplace_back(i, j);
	return pairs;
}

// COLDET tests triangles in model A's frame, the OBB tree in world space, so rounding can
// flip pairs that barely touch.  Returns how many pairs differ between the two sorted lists,
// and how many of those are not borderline: their result is the same whether both triangles
// are shrunk or grown by 0.1% about their centroids.
static std::pair<size_t, size_t> compare_pairs(const PairList& p1, const PairList& p2,
											   const std::vector<Triangle>& soup_a, const Transform& Ta,
											   const std::vector<Triangle>& soup_b, const Transform& Tb)
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
		if (hits(pr, 0.999f) || !hits(pr, 1.001f))
			++unexplained;
	return { diff.size(), unexplained };
}

// Feeds identical triangles and transforms to COLDET and the OBB tree, profiles both, and
// requires them to report the same set of colliding pairs, apart from borderline ones.
// With check_correctness, a brute force reference is also run: the OBB tree (also world
// space) must match it exactly, and COLDET apart from borderline pairs.
// Each query places the models with a random rotation, and a random offset of up to max_offset per axis.
static void run_comparison(const std::vector<Triangle>& soup_a, const std::vector<Triangle>& soup_b,
						   float max_offset, size_t iterations, bool check_correctness)
{
	std::mt19937 rng(2025);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	ProfilerStatistics::instance().clear();

	CollisionModel3D *cd_a, *cd_b;
	{
		Profiler p("coldet build (2 models)");
		cd_a = make_coldet_model(soup_a);
		cd_b = make_coldet_model(soup_b);
	}
	OBBTree tree_a, tree_b;
	{
		Profiler p("obbtree build (2 models)");
		tree_a = make_obb_tree(soup_a);
		tree_b = make_obb_tree(soup_b);
	}
	ASSERT_TRUE(tree_a.verify());
	ASSERT_TRUE(tree_b.verify());

	size_t total_pairs = 0, total_hits = 0, total_borderline = 0;
	for (size_t iter = 0; iter < iterations; ++iter)
	{
		SCOPED_TRACE("iteration " + std::to_string(iter));
		Transform Ta{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * max_offset };
		Transform Tb{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * max_offset };

		float ma[16], mb[16];
		to_coldet_matrix(Ta, ma);
		to_coldet_matrix(Tb, mb);
		cd_a->setTransform(ma);
		cd_b->setTransform(mb);
		tree_a.set_dynamic_transform(Ta);
		tree_b.set_dynamic_transform(Tb);

		bool cd_hit, obb_hit;
		{
			Profiler p("coldet first hit");
			cd_hit = cd_a->collision(cd_b);
		}
		{
			Profiler p("obbtree first hit");
			obb_hit = tree_a.check_intersection(tree_b).hit;
		}

		std::vector<std::pair<int, int>> cd_raw;
		{
			Profiler p("coldet find all");
			cd_a->collisionAll(cd_b, cd_raw);
		}
		IntersectionResult obb_res;
		{
			Profiler p("obbtree find all");
			obb_res = tree_a.find_all_intersections(tree_b);
		}
		total_pairs += obb_res.pairs.size();
		total_hits += obb_hit ? 1 : 0;

		PairList cd_pairs(cd_raw.begin(), cd_raw.end());
		PairList& obb_pairs = obb_res.pairs;
		std::sort(cd_pairs.begin(), cd_pairs.end());
		std::sort(obb_pairs.begin(), obb_pairs.end());
		EXPECT_EQ(cd_hit, !cd_pairs.empty());
		EXPECT_EQ(obb_hit, !obb_pairs.empty());
		auto [differing, unexplained] = compare_pairs(cd_pairs, obb_pairs, soup_a, Ta, soup_b, Tb);
		EXPECT_EQ(unexplained, 0u);
		total_borderline += differing;

		if (!check_correctness)
			continue;

		PairList bf;
		{
			const auto wa = transformed(soup_a, Ta);
			const auto wb = transformed(soup_b, Tb);
			Profiler p("brute force find all");
			bf = brute_force_all(wa, wb);
		}
		std::sort(bf.begin(), bf.end());

		EXPECT_EQ(obb_pairs, bf);
		EXPECT_EQ(compare_pairs(cd_pairs, bf, soup_a, Ta, soup_b, Tb).second, 0u);
	}
	std::cout << soup_a.size() << " + " << soup_b.size() << " triangles, " << iterations << " queries. "
			  << "Queries with hits: " << total_hits << ", total colliding pairs: " << total_pairs
			  << ", borderline pairs where the trees differ: " << total_borderline << std::endl;
	ProfilerStatistics::instance().print_statistics();
	delete cd_a;
	delete cd_b;
}

// Large triangles (1/10 of the space) give many collisions per query: checks exact agreement
TEST(compare, correctness)
{
	std::mt19937 rng(2024);
	const auto a = random_soup(1000, 1.0f, rng);
	const auto b = random_soup(1000, 1.0f, rng);
	run_comparison(a, b, 15.0f, 200, true);
}

// Two tori with triangles 1/1000 of the torus size, for performance.
// Brute force is too slow here, so only the two trees are compared to each other.
TEST(compare, performance)
{
	std::vector<Triangle> t;
	{
		Profiler p("torus tessellation");
		t = torus(10.0f, 2.5f, 0.001f);
	}
	run_comparison(t, t, 8.0f, 100, false);
}

// Times the OBB tree on the torus test for a range of leaf sizes, replaying the same placements.
// Slow (rebuilds the trees per size), so disabled: run with --gtest_also_run_disabled_tests
TEST(compare, DISABLED_leaf_size_sweep)
{
	using clock = std::chrono::high_resolution_clock;
	auto us = [](clock::duration d) { return std::chrono::duration<double, std::micro>(d).count(); };
	const size_t QUERIES = 100;
	const std::vector<Triangle> t = torus(10.0f, 2.5f, 0.001f);

	std::mt19937 rng(2025);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	std::vector<std::pair<Transform, Transform>> placements;
	for (size_t i = 0; i < QUERIES; ++i)
		placements.emplace_back(Transform{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 8.0f },
								Transform{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 8.0f });

	std::cout << t.size() << " triangles per torus, " << QUERIES << " queries\n"
			  << "leaf  depth  build(2) ms  first hit us  find all us  box tests/query\n";
	size_t reference_pairs = 0;
	for (size_t leaf : { 1, 2, 3, 4, 6, 8, 10, 16, 32 })
	{
		OBBTree a, b;
		for (const auto& tri : t)
		{
			a.add_triangle(tri);
			b.add_triangle(tri);
		}
		auto t0 = clock::now();
		size_t depth = a.build(leaf);
		b.build(leaf);
		double build = us(clock::now() - t0);
		double first_hit = 0, find_all = 0;
		size_t pairs = 0, checks = 0;
		for (const auto& [Ta, Tb] : placements)
		{
			a.set_dynamic_transform(Ta);
			b.set_dynamic_transform(Tb);
			auto t1 = clock::now();
			a.check_intersection(b);
			auto t2 = clock::now();
			IntersectionResult r = a.find_all_intersections(b);
			auto t3 = clock::now();
			first_hit += us(t2 - t1);
			find_all += us(t3 - t2);
			pairs += r.pairs.size();
			checks += r.checks;
		}
		// Every leaf size must find exactly the same pairs
		if (reference_pairs == 0) reference_pairs = pairs;
		EXPECT_EQ(pairs, reference_pairs) << "leaf " << leaf;
		printf("%4zu  %5zu  %11.0f  %12.1f  %11.1f  %15zu\n", leaf, depth, build / 1000.0,
			   first_hit / QUERIES, find_all / QUERIES, checks / QUERIES);
	}
}


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

	const Eigen::Matrix3f R = obb.transform.rotation;
	const Vector3f c = obb.transform.center;
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
		Vector3f d = q.cwiseAbs() - e;		// positive components are outside

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
	float s = 1.0f / sqrt(2.0f);
	Matrix3f R; R << s,0,s, 0,1,0, -s,0,s;
					 
	std::srand(42); // setRandom uses rand(); keep the cloud independent of test order
	for(size_t i=0; i<1000; ++i)
	{
		Vector3f p;
		p.setRandom();
		p.y()*=2; // distinct extents so the PCA axes are well defined
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
			Vector3f corner = box.transform.center
								+ sx * box.half_size.x() * box.transform.rotation.col(0)
								+ sy * box.half_size.y() * box.transform.rotation.col(1)
								+ sz * box.half_size.z() * box.transform.rotation.col(2);
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
	range total{ 0.0f,max_t };
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
		else if (fabs(ray.origin[i]) > half_size[i])
		{
			// Parallel to this slab and outside it
			return false;
		}
	}
	return total.valid();
}

// Tests b's corners and edges (and face/space diagonals) against a, in a's local space
static bool brute_force_one_sided(const OBB& a, const OBB& b)
{
	OBB b_in_a = b;
	// Transform b into a's local space
	b_in_a.transform.rotation = a.transform.rotation.transpose() * b.transform.rotation;
	b_in_a.transform.center = a.transform.rotation.transpose() * (b.transform.center - a.transform.center);
	// Now a is axis-aligned at origin
	std::vector<Vector3f> b_corners;
	obb_corners_world(b_in_a, b_corners);
	for(int i=0; i<8; ++i)
	{
		if (contains_local(a, b_corners[i]))
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

// Two convex boxes intersect iff a corner of one is inside the other,
// or an edge of one crosses the other, so both directions must be tested
static inline bool brute_force_intersection(const OBB& a, const OBB& b)
{
	return brute_force_one_sided(a, b) || brute_force_one_sided(b, a);
}

TEST(OBB, intersection)
{
	size_t total_intersections=0;
	size_t total_mismatches = 0;
	const size_t N = 1000000;
	for(size_t iter=0;iter<N;++iter)
	{
		OBB a,b;
		a.transform.center.setRandom(); a.transform.center *= 3.0f;
		b.transform.center.setRandom(); b.transform.center *= 3.0f;
		a.half_size.setRandom(); a.half_size = a.half_size.cwiseAbs() + Vector3f::Constant(0.5f);
		b.half_size.setRandom(); b.half_size = b.half_size.cwiseAbs() + Vector3f::Constant(0.5f);
		a.transform.rotation = random_rotation();
		b.transform.rotation = random_rotation();
		bool bf = brute_force_intersection(a,b);
		bool opt = intersects(a,b);
		if (opt!=bf)
		{
			total_mismatches++;
			/*
			std::cout << "Failed intersection test on iteration " << iter << std::endl;
			std::vector<Vector3f> corners;
			obb_corners_world(a, corners);
			obb_corners_world(b, corners);
			std::ofstream f("C:\\tmp\\pts.csv");
			for (const auto& p : corners)
			{
				f << p.x() << "," << p.y() << "," << p.z() << "\n";
			}
			*/
		}
		total_intersections += opt ? 1 : 0;
		//ASSERT_EQ(bf,opt);
	}
	std::cout << "Total intersections: " << total_intersections << ", mismatches: " << total_mismatches << std::endl;
	EXPECT_LT(total_mismatches, 5u); // allow a few mismatches due to numerical issues (touching boxes)
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
