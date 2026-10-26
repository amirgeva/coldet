// Sweep and prune against the all pairs sphere system, as the number of objects grows.
// Each frame moves 10% of the objects, then checks for collisions.  The scene's
// density is constant, so each object has about the same number of neighbours at any size.
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <vector>
#include "../profiler.h"
#include <multiobject.h>

using namespace COLDET;

namespace
{

using PairSet = std::set<std::pair<int, int>>;

PairSet results(MultiObjectSystem* sys, bool exact)
{
	PairSet res;
	int n = sys->check_collisions(exact);
	for (int i = 0; i < n; ++i)
	{
		CollisionDetails d;
		sys->get_collision(i, d);
		res.insert(d.id1 < d.id2 ? std::make_pair(d.id1, d.id2) : std::make_pair(d.id2, d.id1));
	}
	return res;
}

} // namespace

// Optional argument: the size of the moves (default 0.5, against bounding radii of 0.3..1)
int main(int argc, char* argv[])
{
	const int FRAMES = 100;
	const float step = argc > 1 ? float(std::atof(argv[1])) : 0.5f;
	bool all_equal = true;
	for (int n : { 100, 300, 1000, 3000, 10000 })
	{
		std::mt19937 rng(n);
		std::uniform_real_distribution<float> unit(-1.0f, 1.0f), radius(0.3f, 1.0f);
		// constant density: about 1 object per 8 volume units
		const float space = std::cbrt(float(n) * 8.0f) / 2.0f;
		MultiObjectSystem* sweep = newSweepPruneSystem(n);
		MultiObjectSystem* spheres = newSpheresSystem();
		std::vector<std::vector<float>> pos(n, std::vector<float>(3));
		for (int i = 0; i < n; ++i)
		{
			for (float& p : pos[i]) p = unit(rng) * space;
			const float r = radius(rng);
			sweep->add_sphere_object(r, pos[i].data());
			spheres->add_sphere_object(r, pos[i].data());
		}
		const std::string sweep_name = "sweep and prune, " + std::to_string(n) + " objects, frame";
		const std::string spheres_name = "spheres, " + std::to_string(n) + " objects, frame";
		size_t contacts = 0;
		for (int frame = 0; frame < FRAMES; ++frame)
		{
			std::uniform_int_distribution<int> pick(0, n - 1);
			std::vector<int> moved;
			for (int k = 0; k < n / 10; ++k)
			{
				const int id = pick(rng);
				for (float& p : pos[id]) p += unit(rng) * step;
				moved.push_back(id);
			}
			int a, b;
			{
				Profiler p(sweep_name.c_str());
				for (int id : moved) sweep->move_object(id, pos[id].data());
				a = sweep->check_collisions(true);
			}
			{
				Profiler p(spheres_name.c_str());
				for (int id : moved) spheres->move_object(id, pos[id].data());
				b = spheres->check_collisions(true);
			}
			contacts += a;
			// both systems must report the same contacts (each is tested against brute force)
			if (a != b || results(sweep, true) != results(spheres, true))
			{
				std::printf("MISMATCH: %d objects, frame %d: %d vs %d contacts\n", n, frame, a, b);
				all_equal = false;
			}
		}
		std::printf("%d objects: %.1f contacts per frame\n", n, double(contacts) / FRAMES);
		delete sweep;
		delete spheres;
	}
	ProfilerStatistics::instance().print_statistics();
	return all_equal ? 0 : 1;
}
