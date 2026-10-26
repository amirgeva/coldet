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

// Sweep and prune boxes are closed: touching boxes overlap
bool overlap(const Cube& a, const Cube& b)
{
	for (int i = 0; i < 3; ++i)
		if (std::abs(a.pos[i] - b.pos[i]) > a.radius + b.radius) return false;
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

namespace
{

// On a grid, values snap to multiples of 0.5, so equal endpoints are common
float snap(float x, bool grid) { return grid ? std::round(x * 2.0f) / 2.0f : x; }

// One axis: after every add, move and remove, the endpoints must stay sorted, the indices
// must point at them, and the overlap matrix must match brute force.  Intervals are closed:
// touching intervals overlap.
void axis_random_ops(bool grid)
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
				const bool expected = lo[id] <= hi[other] && lo[other] <= hi[id];
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
			lo[id] = snap(unit(rng) * 5.0f, grid);
			hi[id] = lo[id] + snap(length(rng), grid); // on the grid, possibly of zero length
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
			const float d = snap(unit(rng) * 2.0f, grid);
			axe.move_interval(id, d);
			lo[id] += d;
			hi[id] += d;
			check(step, "move");
		}
	}
}

// Add, move and remove boxes; the incrementally maintained intersections must always match
// space: positions in [-space, space] per axis.  move: moves of up to move per axis, or
// with teleport, moves to anywhere in the space (jumping past many endpoints at once)
void box_random_ops(bool grid, float space = 8.0f, float move = 2.0f, bool teleport = false)
{
	const int N = 120;
	std::mt19937 rng(1);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f), radius(0.2f, 1.5f);
	std::uniform_int_distribution<int> pick(0, N - 1), op(0, 9);
	SNP::SweepPrune sp(N);
	std::vector<Cube> cubes(N);
	size_t touching = 0;
	auto add = [&](int id) {
		Cube& c = cubes[id];
		for (float& p : c.pos) p = snap(unit(rng) * space, grid);
		c.radius = snap(radius(rng), grid);
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
			for (int i = 0; i < 3; ++i)
			{
				delta[i] = teleport ? snap(unit(rng) * space, grid) - cubes[id].pos[i] : snap(unit(rng) * move, grid);
				cubes[id].pos[i] += delta[i];
			}
			sp.move_box(id, delta);
		}
		ASSERT_EQ(as_set(sp.get_intersections()), brute_force(cubes)) << "step " << step << ", operation " << o;
		for (int other = 0; other < N; ++other)
			if (other != id && cubes[id].valid && cubes[other].valid && overlap(cubes[id], cubes[other]))
				for (int i = 0; i < 3; ++i)
					if (std::abs(cubes[id].pos[i] - cubes[other].pos[i]) == cubes[id].radius + cubes[other].radius)
						++touching;
	}
	if (grid)
		EXPECT_GT(touching, 20u) << "the grid must produce touching boxes";
}

} // namespace

TEST(sweep_prune, AxisMatchesBruteForce) { axis_random_ops(false); }
// Equal endpoints everywhere, and zero length intervals
TEST(sweep_prune, AxisTiesMatchBruteForce) { axis_random_ops(true); }
TEST(sweep_prune, MatchesBruteForce) { box_random_ops(false); }
// Boxes on a grid: touching faces, edges and corners, and equal positions
TEST(sweep_prune, TiesMatchBruteForce) { box_random_ops(true); }
// Moves that jump across the whole space, past many endpoints at once
TEST(sweep_prune, LargeMoves) { box_random_ops(false, 8.0f, 0.0f, true); box_random_ops(true, 8.0f, 0.0f, true); }
// Nearly every box overlaps nearly every other
TEST(sweep_prune, DenseScene) { box_random_ops(false, 1.0f, 0.5f); box_random_ops(true, 1.0f, 0.5f); }

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

// A random rotation, as the upper 3x3 of a row vector matrix
Matrix3D random_rotation(std::mt19937& rng)
{
	std::normal_distribution<float> n;
	float q[4], len = 0;
	for (float& x : q) { x = n(rng); len += x * x; }
	len = std::sqrt(len);
	for (float& x : q) x /= len;
	const float w = q[0], x = q[1], y = q[2], z = q[3];
	const float r[3][3] = {
		{ 1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y) },
		{ 2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x) },
		{ 2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y) } };
	Matrix3D m = Matrix3D::Identity;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			m.m[i][j] = r[i][j];
	return m;
}

enum class Kind { Updater, Position, Sphere };

struct Object
{
	Kind kind;
	CollisionModel3D* model = nullptr;	// owned by the scene, shared by both systems
	Placement placement;				// updater objects: rotation and position
	float pos[3];
	float radius = 0;
	int id = -1;						// same in both systems
	bool valid = false;

	// where the object's model is: what the systems must use for exact checks
	Matrix3D matrix() const { return kind == Kind::Updater ? placement.matrix : translation(pos); }
};

struct SceneOptions
{
	int updater_objects = 50, position_objects = 20, sphere_objects = 20;
	size_t shared_models = 0;	// 0: a model per object, otherwise objects share a pool of this many
	bool rotate = false;		// updater objects get random rotations
	float space = 6.0f;			// positions in [-space, space] per axis
	float step = 1.0f;			// moves of up to step per axis
	bool grid = false;			// positions, moves and sphere radii on a 0.25 grid (equal endpoints)
	bool updater_moves = false;	// updater objects move only through their matrix, not move_object
	bool teleport = false;		// moves go anywhere in the space, instead of by up to step
	int point_spheres = 0;		// this many of the sphere objects have radius 0
	int empty_models = 0;		// this many of the mesh objects have a model without triangles
	int capacity = 1024;		// sweep and prune's max_objects
	unsigned seed = 7;
};

// The same objects in both systems, and a brute force reference
class Scene
{
public:
	SceneOptions opt;
	std::vector<Object> objects;
	std::vector<CollisionModel3D*> models;
	MultiObjectSystem* sweep;
	MultiObjectSystem* spheres = newSpheresSystem();
	std::mt19937 rng;

	explicit Scene(const SceneOptions& o = SceneOptions())
		: opt(o), sweep(newSweepPruneSystem(o.capacity)), rng(o.seed)
	{
		std::uniform_real_distribution<float> radius(0.3f, 1.2f);
		const int n = opt.updater_objects + opt.position_objects + opt.sphere_objects;
		const size_t meshes = size_t(opt.updater_objects + opt.position_objects);
		const size_t pool = opt.shared_models ? std::min(opt.shared_models, meshes) : meshes;
		for (size_t i = 0; i < pool; ++i)
			models.push_back(random_model(rng));
		CollisionModel3D* empty = nullptr;
		if (opt.empty_models > 0)
		{
			empty = newCollisionModel3D(false);
			empty->finalize();
			models.push_back(empty);
		}
		objects.resize(n);
		for (int i = 0; i < n; ++i)
		{
			Object& obj = objects[i];
			obj.kind = i < opt.updater_objects ? Kind::Updater
					 : (i < opt.updater_objects + opt.position_objects ? Kind::Position : Kind::Sphere);
			if (obj.kind != Kind::Sphere)
			{
				obj.model = i < opt.empty_models ? empty : models[i % pool];
				obj.radius = obj.model->getRadius();
			}
			else if (i >= n - opt.point_spheres)
				obj.radius = 0.0f;
			else obj.radius = opt.grid ? std::max(0.25f, std::round(radius(rng) * 4.0f) / 4.0f) : radius(rng);
		}
		for (auto& obj : objects) add(obj);
	}

	~Scene()
	{
		delete sweep;
		delete spheres;
		for (auto* m : models) delete m;
	}

	void random_position(float* p)
	{
		std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
		for (int i = 0; i < 3; ++i) p[i] = grid_snap(unit(rng) * opt.space);
	}

	float grid_snap(float x) const { return opt.grid ? std::round(x * 4.0f) / 4.0f : x; }

	void set_placement(Object& o)
	{
		for (int i = 0; i < 3; ++i) o.placement.matrix.m[3][i] = o.pos[i];
	}

	void add(Object& o)
	{
		random_position(o.pos);
		o.placement.matrix = opt.rotate ? random_rotation(rng) : Matrix3D::Identity;
		set_placement(o);
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
		for (int i = 0; i < 3; ++i)
			o.pos[i] = opt.teleport ? grid_snap(unit(rng) * opt.space) : o.pos[i] + grid_snap(unit(rng) * opt.step);
		set_placement(o);
		if (o.kind == Kind::Updater && opt.updater_moves) return;
		sweep->move_object(o.id, o.pos);
		spheres->move_object(o.id, o.pos);
	}

	void remove(Object& o)
	{
		sweep->remove_object(o.id);
		spheres->remove_object(o.id);
		o.valid = false;
	}

	const Object& by_id(int id) const
	{
		for (const auto& o : objects) if (o.valid && o.id == id) return o;
		throw std::runtime_error("no object with id " + std::to_string(id));
	}

	// Pairs whose bounding spheres overlap, or that really collide (exact).  The second
	// model's placement is passed as other_transform, so a shared model is placed correctly.
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
						Matrix3D ma = a.matrix(), mb = b.matrix();
						a.model->setTransform(ma.m[0]);
						hit = a.model->collision(b.model, -1, 0, mb.m[0]);
					}
					else if (a.model || b.model)
					{
						Object& m = a.model ? a : b;
						Object& s = a.model ? b : a;
						m.model->setTransform(m.matrix().m[0]);
						hit = m.model->sphereCollision(s.pos, s.radius);
					}
				}
				if (hit) res.insert(ordered(a.id, b.id));
			}
		return res;
	}

	// Sweep and prune reports the pairs whose (closed) cubes around the bounding spheres overlap
	PairSet cube_pairs() const
	{
		PairSet res;
		for (size_t i = 0; i < objects.size(); ++i)
			for (size_t j = i + 1; j < objects.size(); ++j)
			{
				const Object &a = objects[i], &b = objects[j];
				if (!a.valid || !b.valid) continue;
				bool hit = true;
				for (int k = 0; k < 3; ++k)
					if (std::abs(a.pos[k] - b.pos[k]) > a.radius + b.radius) hit = false;
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

	// Both systems against brute force, exact and proximity.  Returns the exact hits.
	size_t check(const std::string& where)
	{
		const PairSet exact = brute_force(true);
		EXPECT_EQ(results(spheres, true), exact) << where << ": sphere system, exact";
		EXPECT_EQ(results(sweep, true), exact) << where << ": sweep and prune, exact";
		EXPECT_EQ(results(spheres, false), brute_force(false)) << where << ": sphere system, proximity";
		EXPECT_EQ(results(sweep, false), cube_pairs()) << where << ": sweep and prune, proximity (cubes)";
		return exact.size();
	}

	// Random moves, removals and returns, checking after each one.  Returns the exact hits.
	size_t run(int steps)
	{
		std::uniform_int_distribution<size_t> pick(0, objects.size() - 1);
		std::uniform_int_distribution<int> op(0, 9);
		size_t hits = 0;
		for (int step = 0; step < steps; ++step)
		{
			Object& o = objects[pick(rng)];
			if (!o.valid) add(o);
			else if (op(rng) == 0) remove(o);
			else move(o);
			hits += check("step " + std::to_string(step));
			if (::testing::Test::HasFailure()) break;
		}
		return hits;
	}
};

} // namespace

// Both systems against brute force, while objects move, leave and come back
TEST(multi_object, MatchesBruteForce)
{
	Scene s;
	EXPECT_GT(s.run(300), 0u);
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

// Objects sharing their collision model (instances of one mesh): each must be tested at
// its own placement
TEST(multi_object, SharedModels)
{
	SceneOptions opt;
	opt.shared_models = 3;
	Scene s(opt);
	EXPECT_GT(s.run(300), 0u);
}

// Objects on a grid: sweep and prune sees many equal endpoints
TEST(multi_object, GridPositions)
{
	SceneOptions opt;
	opt.grid = true;
	opt.sphere_objects = 60; // with radii on the grid too, so their cubes touch
	Scene s(opt);
	EXPECT_GT(s.run(300), 0u);
}

// Sweep and prune holds at most max_objects objects: adding more fails (-1) without
// disturbing the others, and a freed slot can be used again
TEST(multi_object, Capacity)
{
	const int capacity = 8;
	MultiObjectSystem* sys = newSweepPruneSystem(capacity);
	std::mt19937 rng(3);
	std::uniform_real_distribution<float> unit(-1.0f, 1.0f);
	struct Sphere { float pos[3]; float radius; int id; };
	std::vector<Sphere> spheres;
	auto add = [&]() {
		Sphere s;
		for (float& p : s.pos) p = unit(rng) * 2.0f;
		s.radius = 0.75f;
		s.id = sys->add_sphere_object(s.radius, s.pos);
		return s;
	};
	auto check = [&](const char* where) {
		PairSet expected;
		for (size_t i = 0; i < spheres.size(); ++i)
			for (size_t j = i + 1; j < spheres.size(); ++j)
			{
				float d2 = 0;
				for (int k = 0; k < 3; ++k) d2 += (spheres[i].pos[k] - spheres[j].pos[k]) * (spheres[i].pos[k] - spheres[j].pos[k]);
				float sum = spheres[i].radius + spheres[j].radius;
				if (d2 < sum * sum) expected.insert(ordered(spheres[i].id, spheres[j].id));
			}
		EXPECT_EQ(Scene::results(sys, true), expected) << where;
	};
	for (int i = 0; i < capacity; ++i)
	{
		spheres.push_back(add());
		EXPECT_EQ(spheres.back().id, i);
	}
	check("full");
	for (int i = 0; i < 3; ++i)
		EXPECT_EQ(add().id, -1) << "no room for more objects";
	check("after failed adds");
	sys->remove_object(spheres[2].id);
	spheres.erase(spheres.begin() + 2);
	Sphere s = add();
	EXPECT_EQ(s.id, 2) << "the freed slot is used";
	spheres.push_back(s);
	check("after reusing a slot");
	EXPECT_EQ(add().id, -1);
	delete sys;
}

// Models placed by their updater with arbitrary rotations: the bounding radius is measured
// from the model's origin, so it holds for any rotation, and exact checks use the rotation
TEST(multi_object, RotatedModels)
{
	SceneOptions opt;
	opt.rotate = true;
	Scene s(opt);
	// the rotations must matter: without them, the exact results differ
	{
		const PairSet rotated = s.brute_force(true);
		std::vector<Matrix3D> saved;
		for (auto& o : s.objects)
		{
			saved.push_back(o.placement.matrix);
			o.placement.matrix = translation(o.pos);
		}
		const PairSet unrotated = s.brute_force(true);
		for (size_t i = 0; i < s.objects.size(); ++i) s.objects[i].placement.matrix = saved[i];
		EXPECT_NE(rotated, unrotated) << "rotations don't change the outcome, so the test can't tell";
	}
	EXPECT_GT(s.run(300), 0u);
	SceneOptions shared = opt;
	shared.shared_models = 2;
	shared.seed = 8;
	Scene t(shared);
	EXPECT_GT(t.run(300), 0u) << "shared, rotated models";
}

// Objects placed by an updater are moved by changing what it returns, without move_object:
// the systems must follow the updater
TEST(multi_object, UpdaterDrivenMoves)
{
	SceneOptions opt;
	opt.updater_moves = true;
	opt.rotate = true;
	Scene s(opt);
	EXPECT_GT(s.run(300), 0u);
}

// Moves that jump across the whole scene
TEST(multi_object, LargeMoves)
{
	SceneOptions opt;
	opt.teleport = true;
	Scene s(opt);
	EXPECT_GT(s.run(300), 0u);
}

// A crowded scene: most bounding volumes overlap
TEST(multi_object, DenseScene)
{
	SceneOptions opt;
	opt.updater_objects = 20;
	opt.position_objects = 10;
	opt.sphere_objects = 120;
	opt.space = 1.5f;
	opt.step = 0.5f;
	opt.grid = true;
	Scene s(opt);
	EXPECT_GT(s.run(100), 1000u);
}

// Point objects (radius 0) and models without triangles
TEST(multi_object, PointsAndEmptyModels)
{
	SceneOptions opt;
	opt.point_spheres = 10;
	opt.empty_models = 5;
	opt.space = 3.0f;
	Scene s(opt);
	// they must take part: some proximity pairs involve them
	size_t involved = 0;
	for (const auto& pr : s.brute_force(false))
	{
		const Object &a = s.by_id(pr.first), &b = s.by_id(pr.second);
		for (const Object* o : { &a, &b })
			if (o->radius == 0.0f) ++involved;
	}
	EXPECT_GT(involved, 0u);
	EXPECT_GT(s.run(300), 0u);
}

// Repeated checks, redundant and invalid operations: results stay consistent with brute force
TEST(multi_object, RepeatedAndInvalidOperations)
{
	Scene s;
	for (MultiObjectSystem* sys : { s.sweep, s.spheres })
	{
		const PairSet first = Scene::results(sys, true);
		EXPECT_EQ(Scene::results(sys, true), first) << "checking again without changes";
		EXPECT_EQ(Scene::results(sys, false), Scene::results(sys, false));
	}
	s.check("initial");

	Object& o = s.objects[5];
	const int id = o.id;
	s.remove(o);
	s.check("removed");
	s.sweep->remove_object(id);
	s.spheres->remove_object(id);
	s.check("removed twice");
	const float far[3] = { 100, 100, 100 };
	s.sweep->move_object(id, far);
	s.spheres->move_object(id, far);
	s.check("moved a removed object");

	for (int bad : { -1, 100000 })
	{
		s.sweep->remove_object(bad);
		s.spheres->remove_object(bad);
		s.sweep->move_object(bad, far);
		s.spheres->move_object(bad, far);
	}
	s.check("invalid ids");

	Object& p = s.objects[7];
	s.sweep->move_object(p.id, p.pos);
	s.spheres->move_object(p.id, p.pos);
	s.check("moved to where it is");

	s.add(o);
	EXPECT_EQ(o.id, id);
	s.check("added back");
}
