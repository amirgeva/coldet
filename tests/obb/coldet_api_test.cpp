// Tests of the rest of the COLDET API: ray and sphere queries, collision details,
// transforms, error handling, timeouts and the utility functions
#include <chrono>
#include <cmath>
#include <gtest/gtest.h>
#include "test_helpers.h"

namespace
{

using V3d = Eigen::Vector3d;

struct Model
{
	CollisionModel3D* cm;
	explicit Model(const std::vector<Triangle>& soup, bool is_static = false, const Transform& T = Transform{})
		: cm(newCollisionModel3D(is_static))
	{
		for (const auto& t : soup) cm->addTriangle(t.a.data(), t.b.data(), t.c.data());
		cm->finalize();
		set(T);
	}
	~Model() { delete cm; }
	void set(const Transform& T)
	{
		float m[16];
		to_coldet_matrix(T, m);
		cm->setTransform(m);
	}
	CollisionModel3D* operator->() { return cm; }
};

// Ray / triangle in double precision (Moller-Trumbore).  Returns the ray parameter, or -1 for
// no hit; 'margin' is set to how close the hit is to the triangle's border, relative to its size
double ray_triangle(const V3d& O, const V3d& D, const Triangle& t, double& margin)
{
	const V3d a = t.a.cast<double>(), b = t.b.cast<double>(), c = t.c.cast<double>();
	const V3d e1 = b - a, e2 = c - a, p = D.cross(e2);
	const double det = e1.dot(p);
	margin = 1e300;
	if (std::abs(det) < 1e-12 * e1.norm() * e2.norm()) return -1;
	const V3d s = O - a;
	const double u = s.dot(p) / det;
	const V3d q = s.cross(e1);
	const double v = D.dot(q) / det;
	const double t_hit = e2.dot(q) / det;
	margin = std::min({ u, v, 1 - u - v });
	if (margin < 0 || t_hit < 0) return -1;
	return t_hit;
}

// Closest point on a triangle to p, in double precision (Ericson, Real-Time Collision Detection 5.1.5)
V3d closest_point_on_triangle(const V3d& p, const V3d& a, const V3d& b, const V3d& c)
{
	const V3d ab = b - a, ac = c - a, ap = p - a;
	const double d1 = ab.dot(ap), d2 = ac.dot(ap);
	if (d1 <= 0 && d2 <= 0) return a;
	const V3d bp = p - b;
	const double d3 = ab.dot(bp), d4 = ac.dot(bp);
	if (d3 >= 0 && d4 <= d3) return b;
	const double vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) return a + ab * (d1 / (d1 - d3));
	const V3d cp = p - c;
	const double d5 = ab.dot(cp), d6 = ac.dot(cp);
	if (d6 >= 0 && d5 <= d6) return c;
	const double vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) return a + ac * (d2 / (d2 - d6));
	const double va = d3 * d6 - d5 * d4;
	if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
	const double denom = 1 / (va + vb + vc);
	return a + ab * (vb * denom) + ac * (vc * denom);
}

double distance_to_triangle(const Vector3f& p, const Triangle& t)
{
	const V3d q = p.cast<double>();
	return (q - closest_point_on_triangle(q, t.a.cast<double>(), t.b.cast<double>(), t.c.cast<double>())).norm();
}

Vector3f point_of(const float* p) { return Vector3f(p[0], p[1], p[2]); }

Triangle triangle_of(const float* t) { return Triangle{ point_of(t), point_of(t + 3), point_of(t + 6) }; }

float max_diff(const Triangle& x, const Triangle& y)
{
	return std::max({ (x.a - y.a).norm(), (x.b - y.b).norm(), (x.c - y.c).norm() });
}

} // namespace

// ---- Ray queries ----

TEST(coldet_ray, KnownAnswers)
{
	Model box(box_mesh(Vector3f(1, 1, 1)));
	const float dir[3] = { 1, 0, 0 };
	float p[3];
	int t1, t2;

	const float outside[3] = { -5, 0.25f, 0.3f };
	ASSERT_TRUE(box->rayCollision(outside, dir, true));
	box->getCollisionPoint(p, false);
	EXPECT_NEAR(p[0], -1, 1e-5f);
	EXPECT_NEAR(p[1], 0.25f, 1e-5f);
	EXPECT_NEAR(p[2], 0.3f, 1e-5f);
	box->getCollidingTriangles(t1, t2);
	EXPECT_TRUE(t1 == 8 || t1 == 9) << "closest hit is on the -x face, got triangle " << t1;
	EXPECT_TRUE(box->rayCollision(outside, dir, false)) << "first hit mode";

	const float above[3] = { -5, 2, 0 };
	EXPECT_FALSE(box->rayCollision(above, dir, true)) << "ray passing above the box";
	const float away[3] = { -1, 0, 0 }, outward[3] = { -5, 0, 0 };
	EXPECT_FALSE(box->rayCollision(outward, away, true)) << "ray pointing away from the box";

	EXPECT_FALSE(box->rayCollision(outside, dir, true, 0.0f, 3.0f)) << "segment ending before the box";
	EXPECT_TRUE(box->rayCollision(outside, dir, true, 0.0f, 5.0f)) << "segment reaching the box";
	EXPECT_TRUE(box->rayCollision(outside, dir, true, 2.0f, 5.0f)) << "segment starting at 2, reaching the box";
	EXPECT_FALSE(box->rayCollision(outside, dir, true, 6.5f, 20.0f)) << "segment starting beyond the far face";

	const float inside[3] = { 0, 0.1f, 0.2f };
	ASSERT_TRUE(box->rayCollision(inside, dir, true)) << "ray from inside the box";
	box->getCollisionPoint(p, false);
	EXPECT_NEAR(p[0], 1, 1e-5f);
}

TEST(coldet_ray, TransformedModel)
{
	const Transform T{ z_rotation(0.5f) * x_rotation(0.25f), Vector3f(10, -3, 4) };
	Model box(box_mesh(Vector3f(1, 1, 1)), false, T);
	Model static_box(box_mesh(Vector3f(1, 1, 1)), true, T);
	// a ray towards the box center, in world space
	const Vector3f O = T.apply(Vector3f(-5, 0.2f, 0.1f));
	const Vector3f D = T.rotation * Vector3f(1, 0, 0);
	for (Model* m : { &box, &static_box })
	{
		ASSERT_TRUE((*m)->rayCollision(O.data(), D.data(), true));
		float p[3];
		(*m)->getCollisionPoint(p, false);
		EXPECT_LT((point_of(p) - T.apply(Vector3f(-1, 0.2f, 0.1f))).norm(), 1e-4f);
		(*m)->getCollisionPoint(p, true);
		EXPECT_LT((point_of(p) - Vector3f(-1, 0.2f, 0.1f)).norm(), 1e-4f) << "model space point";
	}
}

// Random rays against a random soup, compared to a double precision reference
TEST(coldet_ray, MatchesReference)
{
	std::mt19937 rng(21);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	const auto soup = random_soup(300, 1.0f, rng);
	const Transform T{ random_rotation(), Vector3f(2, -1, 3) };
	Model model(soup, false, T);
	const auto world = transformed(soup, T);
	size_t compared = 0, hits = 0;
	for (int i = 0; i < 3000; ++i)
	{
		Vector3f O = Vector3f(unit(rng), unit(rng), unit(rng)) * 15.0f + T.center;
		Vector3f D = Vector3f(unit(rng), unit(rng), unit(rng)).normalized();
		double best = 1e300, closest_margin = 1e300;
		size_t best_tri = 0;
		bool borderline = false;
		for (size_t k = 0; k < world.size(); ++k)
		{
			double margin;
			double t = ray_triangle(O.cast<double>(), D.cast<double>(), world[k], margin);
			if (std::abs(margin) < 1e-4) borderline = true;
			if (t >= 0 && t < 1e-3) borderline = true;
			if (t >= 0 && t < best) { best = t; best_tri = k; }
			if (t >= 0) closest_margin = std::min(closest_margin, t);
		}
		if (borderline) continue;
		const bool expected = best < 1e300;
		++compared;
		hits += expected ? 1 : 0;
		EXPECT_EQ(model->rayCollision(O.data(), D.data(), false), expected) << "first hit, ray " << i;
		ASSERT_EQ(model->rayCollision(O.data(), D.data(), true), expected) << "closest, ray " << i;
		if (expected)
		{
			float p[3];
			model->getCollisionPoint(p, false);
			const Vector3f ref = O + D * float(best);
			EXPECT_LT((point_of(p) - ref).norm(), 1e-3f) << "closest point, ray " << i;
			int t1, t2;
			model->getCollidingTriangles(t1, t2);
			double margin;
			EXPECT_NEAR(ray_triangle(O.cast<double>(), D.cast<double>(), world[t1], margin), best, 1e-3)
				<< "reported triangle " << t1 << " vs closest " << best_tri;
		}
	}
	EXPECT_GT(hits, 100u);
	std::cout << "Compared " << compared << " rays, " << hits << " hitting" << std::endl;
}

// ---- Sphere queries ----

TEST(coldet_sphere, KnownAnswers)
{
	Model box(box_mesh(Vector3f(1, 1, 1)));
	struct Case { const char* name; Vector3f center; float radius; bool hit; };
	const Case cases[] = {
		{ "near face, reaching", { -1.5f, 0.2f, 0.1f }, 0.6f, true },
		{ "near face, short", { -1.5f, 0.2f, 0.1f }, 0.4f, false },
		{ "near edge, reaching", { 1.3f, 1.3f, 0 }, 0.5f, true },		// distance 0.424
		{ "near edge, short", { 1.3f, 1.3f, 0 }, 0.4f, false },
		{ "near corner, reaching", { 1.3f, 1.3f, 1.3f }, 0.6f, true },	// distance 0.520
		{ "near corner, short", { 1.3f, 1.3f, 1.3f }, 0.5f, false },
		{ "inside, not touching", { 0, 0, 0 }, 0.5f, false },			// surfaces only
		{ "inside, touching", { 0.8f, 0, 0 }, 0.5f, true },
		{ "containing the box", { 0, 0, 0 }, 3.0f, true },				// a solid ball overlaps every triangle
	};
	for (const auto& c : cases)
	{
		EXPECT_EQ(box->sphereCollision(c.center.data(), c.radius), c.hit) << c.name;
		if (c.hit)
		{
			float p[3];
			box->getCollisionPoint(p, false);
			EXPECT_LE((point_of(p) - c.center).norm(), c.radius + 1e-5f) << c.name << ": point within the sphere";
			int t1, t2;
			box->getCollidingTriangles(t1, t2);
			EXPECT_LT(distance_to_triangle(point_of(p), box_mesh(Vector3f(1, 1, 1))[t1]), 1e-5) << c.name << ": point on the triangle";
		}
	}
}

// Random spheres against a random soup, compared to a double precision reference
TEST(coldet_sphere, MatchesReference)
{
	std::mt19937 rng(31);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f), radius(0.05f, 1.5f);
	const auto soup = random_soup(300, 1.0f, rng);
	const Transform T{ random_rotation(), Vector3f(-4, 2, 1) };
	for (bool is_static : { false, true })
	{
		Model model(soup, is_static, T);
		const auto world = transformed(soup, T);
		size_t compared = 0, hits = 0;
		for (int i = 0; i < 2000; ++i)
		{
			const Vector3f C = Vector3f(unit(rng), unit(rng), unit(rng)) * 12.0f + T.center;
			const float r = radius(rng);
			double nearest = 1e300;
			for (const auto& t : world) nearest = std::min(nearest, distance_to_triangle(C, t));
			if (std::abs(nearest - r) < 1e-4) continue;
			const bool expected = nearest < r;
			++compared;
			hits += expected ? 1 : 0;
			ASSERT_EQ(model->sphereCollision(C.data(), r), expected) << "sphere " << i << ", static " << is_static;
			if (expected)
			{
				float p[3];
				model->getCollisionPoint(p, false);
				EXPECT_LE((point_of(p) - C).norm(), r + 1e-4f);
				int t1, t2;
				model->getCollidingTriangles(t1, t2);
				EXPECT_LT(distance_to_triangle(point_of(p), world[t1]), 1e-4);
			}
		}
		EXPECT_GT(hits, 100u);
		std::cout << "Compared " << compared << " spheres, " << hits << " hitting" << std::endl;
	}
}

// ---- Model collision details ----

// After a hit, the reported triangles and collision point must describe a real contact
TEST(coldet_collision, Details)
{
	std::mt19937 rng(41);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	const auto sa = random_soup(400, 1.0f, rng), sb = random_soup(400, 1.0f, rng);
	Model a(sa), b(sb);
	size_t hits = 0;
	for (int i = 0; i < 100; ++i)
	{
		const Transform Ta{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 5.0f };
		const Transform Tb{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 5.0f };
		a.set(Ta);
		b.set(Tb);
		if (!a->collision(b.cm)) continue;
		++hits;
		int i1, i2;
		a->getCollidingTriangles(i1, i2);
		const Triangle wa = transformed({ sa[i1] }, Ta)[0], wb = transformed({ sb[i2] }, Tb)[0];
		EXPECT_TRUE(tri_tri_intersect(wa.a.data(), wa.b.data(), wa.c.data(), wb.a.data(), wb.b.data(), wb.c.data()));

		float t1[9], t2[9];
		a->getCollidingTriangles(t1, t2, true);
		EXPECT_EQ(max_diff(triangle_of(t1), sa[i1]), 0.0f) << "model space: this model's triangle";
		const Transform inv_a{ Ta.rotation.transpose(), -(Ta.rotation.transpose() * Ta.center) };
		EXPECT_LT(max_diff(triangle_of(t2), transformed({ wb }, inv_a)[0]), 1e-4f) << "model space: other triangle, in this model's space";
		a->getCollidingTriangles(t1, t2, false);
		EXPECT_LT(max_diff(triangle_of(t1), wa), 1e-4f) << "world space: this model's triangle";
		EXPECT_LT(max_diff(triangle_of(t2), wb), 1e-4f) << "world space: other triangle";

		float p[3];
		a->getCollisionPoint(p, false);
		EXPECT_LT(distance_to_triangle(point_of(p), wa), 1e-3) << "collision point on this model's triangle";
		EXPECT_LT(distance_to_triangle(point_of(p), wb), 1e-3) << "collision point on the other triangle";
		a->getCollisionPoint(p, true);
		EXPECT_LT(distance_to_triangle(point_of(p), sa[i1]), 1e-3) << "model space collision point";
	}
	EXPECT_GT(hits, 20u);
}

// other_transform overrides the other model's own transform
TEST(coldet_collision, OtherTransformOverride)
{
	std::mt19937 rng(51);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	const auto sa = random_soup(300, 1.0f, rng), sb = random_soup(300, 1.0f, rng);
	Model a(sa), b(sb), b2(sb);
	size_t hits = 0;
	for (int i = 0; i < 50; ++i)
	{
		const Transform Ta{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 10.0f };
		const Transform Tb{ random_rotation(), Vector3f(unit(rng), unit(rng), unit(rng)) * 10.0f };
		a.set(Ta);
		b.set(Transform{ Matrix3f::Identity(), Vector3f(1000, 0, 0) }); // overridden below
		b2.set(Tb);
		float mb[16];
		to_coldet_matrix(Tb, mb);
		const bool expected = a->collision(b2.cm);
		EXPECT_EQ(a->collision(b.cm, -1, 0, mb), expected);
		std::vector<std::pair<int, int>> p1, p2;
		a->collisionAll(b2.cm, p1);
		a->collisionAll(b.cm, p2, mb);
		EXPECT_EQ(p1, p2);
		hits += expected ? 1 : 0;
	}
	EXPECT_GT(hits, 0u);
}

// A model tested against itself, using other_transform for the second placement
TEST(coldet_collision, SelfWithOtherTransform)
{
	const auto t = torus(2.0f, 0.5f, 0.05f);
	Model m(t);
	float mb[16];
	to_coldet_matrix(Transform{ x_rotation(1.5707963f), Vector3f(2, 0, 0) }, mb); // a linked ring
	std::vector<std::pair<int, int>> pairs;
	EXPECT_FALSE(m->collisionAll(m.cm, pairs, mb)) << "linked rings don't touch";
	to_coldet_matrix(Transform{ Matrix3f::Identity(), Vector3f(1, 0, 0) }, mb);
	EXPECT_TRUE(m->collisionAll(m.cm, pairs, mb)) << "rings in one plane, cutting through each other";
	EXPECT_FALSE(pairs.empty());
}

// ---- Model construction ----

TEST(coldet_model, Inconsistencies)
{
	const float v0[3] = { 0, 0, 0 }, v1[3] = { 1, 0, 0 }, v2[3] = { 0, 1, 0 };
	CollisionModel3D* a = newCollisionModel3D(false);
	CollisionModel3D* b = newCollisionModel3D(false);
	a->addTriangle(v0, v1, v2);
	b->addTriangle(v0, v1, v2);
	b->finalize();
	EXPECT_THROW(a->collision(b), Inconsistency) << "collision before finalize";
	EXPECT_THROW(b->collision(a), Inconsistency) << "collision with a model that is not finalized";
	a->finalize();
	EXPECT_NO_THROW(a->collision(b));
	EXPECT_THROW(a->finalize(), Inconsistency) << "finalize twice";
	EXPECT_THROW(a->addTriangle(v0, v1, v2), Inconsistency) << "add triangle after finalize";
	EXPECT_THROW(a->addTriangle(0, 0, 0, 1, 0, 0, 0, 1, 0), Inconsistency) << "add triangle (coordinates) after finalize";
	delete a;
	delete b;
}

TEST(coldet_model, RadiusAndTriangleNumber)
{
	std::mt19937 rng(61);
	const auto soup = random_soup(200, 1.0f, rng);
	float expected = 0;
	for (const auto& t : soup) expected = std::max({ expected, t.a.norm(), t.b.norm(), t.c.norm() });
	CollisionModel3D* m = newCollisionModel3D(false);
	m->setTriangleNumber(int(soup.size()));
	for (const auto& t : soup) m->addTriangle(t.a.x(), t.a.y(), t.a.z(), t.b.x(), t.b.y(), t.b.z(), t.c.x(), t.c.y(), t.c.z());
	m->finalize();
	EXPECT_NEAR(m->getRadius(), expected, 1e-5f);
	// both addTriangle forms and setTriangleNumber build the same model
	Model reference(soup);
	const Transform T{ random_rotation(), Vector3f(1, 2, 0) };
	float mt[16];
	to_coldet_matrix(T, mt);
	m->setTransform(mt);
	reference.set(T);
	Model probe(random_soup(200, 1.0f, rng));
	std::vector<std::pair<int, int>> p1, p2;
	m->collisionAll(probe.cm, p1);
	reference->collisionAll(probe.cm, p2);
	std::sort(p1.begin(), p1.end());
	std::sort(p2.begin(), p2.end());
	EXPECT_EQ(p1, p2);
	EXPECT_FALSE(p1.empty());
	delete m;
}

// ---- Timeout ----

// A query that takes longer than MaxProcessingTime throws TimeoutExpired
TEST(coldet_model, Timeout)
{
	// Two wavy surfaces 0.01 apart: the large boxes overlap, so the search descends deep into
	// both trees before rejecting, but the surfaces never touch
	auto surface = [](size_t n, float lift) {
		std::vector<Triangle> res;
		auto P = [&](size_t i, size_t j) {
			float x = 10.0f * float(i) / float(n), y = 10.0f * float(j) / float(n);
			return Vector3f(x, y, 0.2f * std::sin(x) * std::cos(y) + lift);
		};
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < n; ++j)
			{
				res.push_back(Triangle{ P(i, j), P(i + 1, j), P(i + 1, j + 1) });
				res.push_back(Triangle{ P(i, j), P(i + 1, j + 1), P(i, j + 1) });
			}
		return res;
	};
	Model a(surface(600, 0.0f)), b(surface(600, 0.01f));
	std::vector<std::pair<int, int>> pairs;
	auto t0 = std::chrono::steady_clock::now();
	a->collisionAll(b.cm, pairs);
	double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
	ASSERT_TRUE(pairs.empty()) << "the surfaces must not touch";
	if (ms < 10.0)
		GTEST_SKIP() << "the query takes only " << ms << " ms, too fast to time out reliably";
	EXPECT_THROW(a->collision(b.cm, -1, 1), TimeoutExpired) << "query takes " << ms << " ms";
	EXPECT_FALSE(a->collision(b.cm, -1, 60000)) << "a generous limit does not time out";
}

// ---- Utility functions ----

TEST(coldet_utils, SphereRayCollision)
{
	const float center[3] = { 5, 0, 0 };
	float p[3];
	const float origin[3] = { 0, 0, 0 };
	const float toward[3] = { 2, 0, 0 }; // not normalized
	ASSERT_TRUE(SphereRayCollision(center, 1, origin, toward, p)) << "ray toward the sphere";
	EXPECT_NEAR(p[0], 4, 1e-5f);
	EXPECT_NEAR(p[1], 0, 1e-5f);

	const float miss[3] = { 1, 1, 0 };
	EXPECT_FALSE(SphereRayCollision(center, 1, origin, miss, p)) << "ray passing the sphere";

	const float away[3] = { -1, 0, 0 };
	EXPECT_FALSE(SphereRayCollision(center, 1, origin, away, p)) << "sphere behind the ray";

	const float inside[3] = { 5, 0, 0 };
	ASSERT_TRUE(SphereRayCollision(center, 1, inside, toward, p)) << "ray from inside the sphere";
	EXPECT_NEAR(p[0], 6, 1e-5f) << "exits through the far side";
}

TEST(coldet_utils, SphereSphereCollision)
{
	float c1[3] = { 0, 0, 0 }, c2[3] = { 3, 4, 0 }, p[3];
	ASSERT_TRUE(SphereSphereCollision(c1, 2, c2, 4, p)) << "overlapping spheres (distance 5)";
	// the point is on the first sphere, toward the second
	EXPECT_NEAR(p[0], 1.2f, 1e-5f);
	EXPECT_NEAR(p[1], 1.6f, 1e-5f);
	EXPECT_NEAR(p[2], 0, 1e-5f);
	EXPECT_FALSE(SphereSphereCollision(c1, 2, c2, 2.5f, p)) << "separate spheres";
	float c3[3] = { 0, 0, 0 };
	EXPECT_TRUE(SphereSphereCollision(c1, 1, c3, 1, p)) << "concentric spheres";
}
