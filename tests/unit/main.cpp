// Tests of the multiple object systems (sweep and prune, and spheres) against brute force
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <set>
#include <vector>
#include <gtest/gtest.h>
#include <coldet.h>
#include <cdmath3d.h>
#include <multiobject.h>
#include <sweep_prune.h>
#include <tritri.h>

using namespace COLDET;

namespace
{

using PairSet = std::set<std::pair<int, int>>;

std::pair<int, int> ordered(int a, int b) { return a < b ? std::make_pair(a, b) : std::make_pair(b, a); }

// ---- Sweep and prune boxes ----

struct Cube { float pos[3]; float radius; bool valid = false; };

bool overlap(const Cube& a, const Cube& b)
{
	for (int i = 0; i < 3; ++i)
		if (std::abs(a.pos[i] - b.pos[i]) >= a.radius + b.radius) return false;
	return true;
}

PairSet brute_force(const std::vector<Cube>& cubes)
{
	PairSet res;
	for (int i = 0; i < int(cubes.size()); ++i)
		for (int j = i + 1; j < int(cubes.size()); ++j)
			if (cubes[i].valid && cubes[j].valid && overlap(cubes[i], cubes[j]))
				res.insert({ i, j });
	return res;
}

PairSet as_set(const SNP::intr_set& s)
{
	PairSet res;
	for (const auto& i : s) res.insert(ordered(i.id1, i.id2));
	return res;
}

} // namespace

// One axis: after every add, move and remove, the endpoints must stay sorted, the indices
// must point at them, and the overlap matrix must match brute force
TEST(sweep_prune, AxisMatchesBruteForce)
{
	const int N = 40;
	std::mt19937 rng(2);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f), length(0.1f, 2.0f);
	std::uniform_int_distribution<int> pick(0, N - 1), op(0, 9);
	SNP::SweepPruneAxe axe;
	axe.initialize(N);
	std::vector<float> lo(N), hi(N);
	std::vector<bool> valid(N, false);
	auto check = [&](int step, const char* what) {
		const auto& eps = axe.get_end_points();
		for (size_t i = 1; i < eps.size(); ++i)
			ASSERT_LE(eps[i - 1].pos, eps[i].pos) << what << " at step " << step << ": endpoints not sorted at " << i;
		for (int id = 0; id < N; ++id)
		{
			if (!valid[id]) continue;
			const auto& ind = axe.get_indices()[id];
			ASSERT_EQ(eps[ind.start].id, id) << what << " at step " << step;
			ASSERT_TRUE(eps[ind.start].start);
			ASSERT_EQ(eps[ind.stop].id, id);
			ASSERT_FALSE(eps[ind.stop].start);
			ASSERT_NEAR(eps[ind.start].pos, lo[id], 1e-4f);
			ASSERT_NEAR(eps[ind.stop].pos, hi[id], 1e-4f);
			for (int other = 0; other < N; ++other)
			{
				if (!valid[other] || other == id) continue;
				const bool expected = lo[id] < hi[other] && lo[other] < hi[id];
				ASSERT_EQ(axe.get_matrix_row(id).test(other), expected)
					<< what << " at step " << step << ": " << id << " [" << lo[id] << "," << hi[id] << "] vs "
					<< other << " [" << lo[other] << "," << hi[other] << "]";
			}
		}
	};
	for (int step = 0; step < 5000; ++step)
	{
		const int id = pick(rng);
		if (!valid[id])
		{
			lo[id] = unit(rng) * 5.0f;
			hi[id] = lo[id] + length(rng);
			ASSERT_TRUE(axe.add_interval(SNP::ushort(id), lo[id], hi[id]));
			valid[id] = true;
			check(step, "add");
		}
		else if (op(rng) == 0)
		{
			axe.remove_interval(id);
			valid[id] = false;
			check(step, "remove");
		}
		else
		{
			const float d = unit(rng) * 2.0f;
			axe.move_interval(id, d);
			lo[id] += d;
			hi[id] += d;
			check(step, "move");
		}
	}
}

// Add, move and remove boxes; the incrementally maintained intersections must always match
TEST(sweep_prune, MatchesBruteForce)
{
	const int N = 120;
	std::mt19937 rng(1);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f), radius(0.2f, 1.5f);
	std::uniform_int_distribution<int> pick(0, N - 1), op(0, 9);
	SNP::SweepPrune sp(N);
	std::vector<Cube> cubes(N);
	auto add = [&](int id) {
		Cube& c = cubes[id];
		for (float& p : c.pos) p = unit(rng) * 8.0f;
		c.radius = radius(rng);
		c.valid = true;
		sp.add_box(id, c.pos, c.radius);
	};
	for (int id = 0; id < N; ++id) add(id);
	ASSERT_EQ(as_set(sp.get_intersections()), brute_force(cubes)) << "after adding";
	for (int step = 0; step < 2000; ++step)
	{
		const int id = pick(rng);
		const int o = op(rng);
		if (!cubes[id].valid) add(id);
		else if (o == 0)
		{
			sp.remove_box(id);
			cubes[id].valid = false;
		}
		else
		{
			float delta[3];
			for (int i = 0; i < 3; ++i) { delta[i] = unit(rng) * 2.0f; cubes[id].pos[i] += delta[i]; }
			sp.move_box(id, delta);
		}
		ASSERT_EQ(as_set(sp.get_intersections()), brute_force(cubes)) << "step " << step << ", operation " << o;
	}
}

namespace
{

// ---- Multiple object systems ----

float rand_float(std::mt19937& rng) { return std::uniform_real_distribution<float>(-0.5f, 0.5f)(rng); }

CollisionModel3D* random_model(std::mt19937& rng)
{
	CollisionModel3D* m = newCollisionModel3D(false);
	for (int i = 0; i < 32; ++i)
	{
		float v[9];
		for (float& x : v) x = rand_float(rng);
		m->addTriangle(&v[0], &v[3], &v[6]);
	}
	m->finalize();
	return m;
}

Matrix3D translation(const float* pos)
{
	Matrix3D m = Matrix3D::Identity;
	m.m[3][0] = pos[0];
	m.m[3][1] = pos[1];
	m.m[3][2] = pos[2];
	return m;
}

// A placement the system reads through TransformUpdater
class Placement : public TransformUpdater
{
public:
	Matrix3D matrix = Matrix3D::Identity;
	const float* update() const override { return reinterpret_cast<const float*>(&matrix); }
};

enum class Kind { Updater, Position, Sphere };

struct Object
{
	Kind kind;
	CollisionModel3D* model = nullptr;	// owned by the scene, shared by both systems
	Placement placement;
	float pos[3];
	float radius = 0;
	int id = -1;						// same in both systems
	bool valid = false;
};

// The same objects in both systems, and a brute force reference
class Scene
{
public:
	std::vector<Object> objects;
	MultiObjectSystem* sweep = newSweepPruneSystem(1024);
	MultiObjectSystem* spheres = newSpheresSystem();
	std::mt19937 rng{ 7 };

	~Scene()
	{
		delete sweep;
		delete spheres;
		for (auto& o : objects) delete o.model;
	}

	void random_position(float* p)
	{
		std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
		for (int i = 0; i < 3; ++i) p[i] = unit(rng) * 6.0f;
	}

	void add(Object& o)
	{
		random_position(o.pos);
		o.placement.matrix = translation(o.pos);
		int a = -1, b = -1;
		switch (o.kind)
		{
		case Kind::Updater:
			a = sweep->add_object(o.model, &o.placement);
			b = spheres->add_object(o.model, &o.placement);
			break;
		case Kind::Position:
			a = sweep->add_object(o.model, o.pos);
			b = spheres->add_object(o.model, o.pos);
			break;
		case Kind::Sphere:
			a = sweep->add_sphere_object(o.radius, o.pos);
			b = spheres->add_sphere_object(o.radius, o.pos);
			break;
		}
		ASSERT_EQ(a, b) << "both systems allocate the same ids";
		o.id = a;
		o.valid = true;
	}

	void move(Object& o)
	{
		std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
		for (int i = 0; i < 3; ++i) o.pos[i] += unit(rng);
		o.placement.matrix = translation(o.pos);
		sweep->move_object(o.id, o.pos);
		spheres->move_object(o.id, o.pos);
	}

	void remove(Object& o)
	{
		sweep->remove_object(o.id);
		spheres->remove_object(o.id);
		o.valid = false;
	}

	Scene()
	{
		std::uniform_real_distribution<float> radius(0.3f, 1.2f);
		objects.resize(90);
		for (size_t i = 0; i < objects.size(); ++i)
		{
			Object& o = objects[i];
			o.kind = i < 50 ? Kind::Updater : (i < 70 ? Kind::Position : Kind::Sphere);
			if (o.kind != Kind::Sphere)
			{
				o.model = random_model(rng);
				o.radius = o.model->getRadius();
			}
			else o.radius = radius(rng);
		}
		for (auto& o : objects) add(o);
	}

	const Object& by_id(int id) const
	{
		for (const auto& o : objects) if (o.valid && o.id == id) return o;
		throw std::runtime_error("no object with id " + std::to_string(id));
	}

	// Pairs whose bounding spheres overlap, or that really collide (exact)
	PairSet brute_force(bool exact)
	{
		PairSet res;
		for (size_t i = 0; i < objects.size(); ++i)
			for (size_t j = i + 1; j < objects.size(); ++j)
			{
				Object &a = objects[i], &b = objects[j];
				if (!a.valid || !b.valid) continue;
				float d2 = 0;
				for (int k = 0; k < 3; ++k) d2 += (a.pos[k] - b.pos[k]) * (a.pos[k] - b.pos[k]);
				float sum = a.radius + b.radius;
				if (d2 >= sum * sum) continue;
				bool hit = true;
				if (exact)
				{
					if (a.model && b.model)
					{
						a.model->setTransform(translation(a.pos).m[0]);
						b.model->setTransform(translation(b.pos).m[0]);
						hit = a.model->collision(b.model);
					}
					else if (a.model || b.model)
					{
						Object& m = a.model ? a : b;
						Object& s = a.model ? b : a;
						m.model->setTransform(translation(m.pos).m[0]);
						hit = m.model->sphereCollision(s.pos, s.radius);
					}
				}
				if (hit) res.insert(ordered(a.id, b.id));
			}
		return res;
	}

	static PairSet results(MultiObjectSystem* sys, bool exact)
	{
		PairSet res;
		int n = sys->check_collisions(exact);
		for (int i = 0; i < n; ++i)
		{
			CollisionDetails d;
			EXPECT_TRUE(sys->get_collision(i, d));
			res.insert(ordered(d.id1, d.id2));
		}
		return res;
	}
};

// Sweep and prune reports the pairs whose cubes around the bounding spheres overlap
PairSet cube_pairs(const Scene& s)
{
	PairSet res;
	for (size_t i = 0; i < s.objects.size(); ++i)
		for (size_t j = i + 1; j < s.objects.size(); ++j)
		{
			const Object &a = s.objects[i], &b = s.objects[j];
			if (!a.valid || !b.valid) continue;
			bool hit = true;
			for (int k = 0; k < 3; ++k)
				if (std::abs(a.pos[k] - b.pos[k]) >= a.radius + b.radius) hit = false;
			if (hit) res.insert(ordered(a.id, b.id));
		}
	return res;
}

} // namespace

// Both systems against brute force, while objects move, leave and come back
TEST(multi_object, MatchesBruteForce)
{
	Scene s;
	std::uniform_int_distribution<size_t> pick(0, s.objects.size() - 1);
	std::uniform_int_distribution<int> op(0, 9);
	size_t exact_hits = 0;
	for (int step = 0; step < 300; ++step)
	{
		SCOPED_TRACE("step " + std::to_string(step));
		Object& o = s.objects[pick(s.rng)];
		if (!o.valid) s.add(o);
		else if (op(s.rng) == 0) s.remove(o);
		else s.move(o);

		const PairSet exact = s.brute_force(true);
		EXPECT_EQ(Scene::results(s.spheres, true), exact) << "sphere system, exact";
		EXPECT_EQ(Scene::results(s.sweep, true), exact) << "sweep and prune, exact";
		EXPECT_EQ(Scene::results(s.spheres, false), s.brute_force(false)) << "sphere system, proximity";
		EXPECT_EQ(Scene::results(s.sweep, false), cube_pairs(s)) << "sweep and prune, proximity (cubes)";
		exact_hits += exact.size();
		if (::testing::Test::HasFailure()) break;
	}
	EXPECT_GT(exact_hits, 0u);
}

// The details of exact collisions describe real contacts
TEST(multi_object, CollisionDetails)
{
	Scene s;
	size_t checked = 0;
	for (MultiObjectSystem* sys : { s.sweep, s.spheres })
	{
		int n = sys->check_collisions(true);
		for (int i = 0; i < n; ++i)
		{
			CollisionDetails d;
			ASSERT_TRUE(sys->get_collision(i, d));
			const Object &a = s.by_id(d.id1), &b = s.by_id(d.id2);
			if (a.model && b.model)
			{
				EXPECT_TRUE(tri_tri_intersect(d.t1, d.t1 + 3, d.t1 + 6, d.t2, d.t2 + 3, d.t2 + 6))
					<< "the reported world space triangles intersect";
				++checked;
			}
			// the collision point is near both objects
			for (const Object* o : { &a, &b })
			{
				float d2 = 0;
				for (int k = 0; k < 3; ++k) d2 += (d.point[k] - o->pos[k]) * (d.point[k] - o->pos[k]);
				EXPECT_LE(std::sqrt(d2), o->radius + 1e-3f);
			}
		}
		CollisionDetails d;
		EXPECT_FALSE(sys->get_collision(-1, d));
		EXPECT_FALSE(sys->get_collision(n, d));
	}
	EXPECT_GT(checked, 0u);
}

TEST(multi_object, CollisionModels)
{
	Scene s;
	for (MultiObjectSystem* sys : { s.sweep, s.spheres })
	{
		for (const auto& o : s.objects)
			EXPECT_EQ(sys->get_collision_model(o.id), o.model) << "id " << o.id;
		EXPECT_EQ(sys->get_collision_model(-1), nullptr);
		EXPECT_EQ(sys->get_collision_model(100000), nullptr);
		EXPECT_EQ(sys->add_object(static_cast<CollisionModel3D*>(nullptr), s.objects[0].pos), -1);
	}
	// removed objects have no model, and their ids are reused
	Object& o = s.objects[3];
	const int id = o.id;
	s.remove(o);
	EXPECT_EQ(s.sweep->get_collision_model(id), nullptr);
	EXPECT_EQ(s.spheres->get_collision_model(id), nullptr);
	s.add(o);
	EXPECT_EQ(o.id, id) << "the free slot is reused";
	EXPECT_EQ(s.sweep->get_collision_model(id), o.model);
}

// Removed objects no longer collide
TEST(multi_object, RemovedObjectsDontCollide)
{
	Scene s;
	for (auto& o : s.objects) s.remove(o);
	EXPECT_EQ(s.sweep->check_collisions(false), 0);
	EXPECT_EQ(s.spheres->check_collisions(false), 0);
	EXPECT_EQ(s.sweep->check_collisions(true), 0);
	EXPECT_EQ(s.spheres->check_collisions(true), 0);
}
