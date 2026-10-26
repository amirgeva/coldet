# ColDet

A C++ library for 3D collision detection between triangle meshes.

ColDet answers "do these two models touch, and where?" for arbitrary triangle soups placed
anywhere in space with rotations and translations. It also handles rays, spheres, and scenes
with many moving objects. It started in 2000 as a SourceForge project. This repository is
a modernized version, with a CMake build, an extensive test suite and a number of robustness
fixes. It also includes an experimental oriented bounding box (OBB) tree as an alternative
to the original hierarchy.

- **Mesh against mesh**: first hit, or every colliding triangle pair, with the triangles and
  a point of contact.
- **Ray and sphere queries** against a mesh: first or closest hit along a ray or segment.
- **Multiple object scenes**: a sweep and prune broad phase (or an all pairs sphere test),
  followed by exact triangle tests.
- **Robust at any scale**: scale independent tolerances, coplanar and degenerate triangles,
  and models far from the origin are all handled and tested.

## Contents

- [Building](#building)
- [Using the library](#using-the-library)
  - [Collision models](#collision-models)
  - [Transforms](#transforms)
  - [Ray and sphere queries](#ray-and-sphere-queries)
  - [Scenes with many objects](#scenes-with-many-objects)
  - [The triangle intersection test](#the-triangle-intersection-test)
  - [The OBB tree (experimental)](#the-obb-tree-experimental)
- [Behaviour and tolerances](#behaviour-and-tolerances)
- [Tests](#tests)
- [Benchmarks](#benchmarks)
- [Performance](#performance)
- [Repository layout](#repository-layout)
- [Known limitations](#known-limitations)
- [License](#license)

## Building

Requirements:

- A C++17 compiler (developed with Visual Studio 2022)
- CMake 3.20 or later
- [Eigen 3](https://eigen.tuxfamily.org) (used by the OBB library)
- [GoogleTest](https://github.com/google/googletest), for the tests
- TBB is detected if present, but nothing uses it at the moment

With [vcpkg](https://vcpkg.io):

```sh
vcpkg install eigen3 gtest
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
```

This builds:

| Target | What |
|---|---|
| `coldet` | The collision detection library (static) |
| `obb` | The OBB tree library (static, depends on `coldet` and Eigen) |
| `test.obb`, `coldet_unit_tests` | The tests (see [Tests](#tests)) |
| `bm.triangles`, `bm.multiobject` | Benchmarks |

Turn off the tests with `-DBUILD_TESTING=OFF`. GoogleTest is then not needed.

## Using the library

All public headers are in `src/`. The model API (`coldet.h`) is in the global namespace. The
multiple object systems and the triangle test are in namespace `COLDET`.

### Collision models

A `CollisionModel3D` holds the triangles of one mesh, in the mesh's own coordinates. Add the
triangles, call `finalize()` to build the bounding box hierarchy, then place the model with
a transform and test it against other models.

```cpp
#include <coldet.h>

CollisionModel3D* a = newCollisionModel3D();      // dynamic model
CollisionModel3D* b = newCollisionModel3D(true);  // static: its transform rarely changes

const float v0[3] = { 0, 0, 0 }, v1[3] = { 1, 0, 0 }, v2[3] = { 0, 1, 0 };
a->addTriangle(v0, v1, v2);                       // ... all of a's triangles
b->addTriangle(0, 0, -1, 0, 0, 1, 1, 1, 0);       // or with coordinates
a->finalize();
b->finalize();

a->setTransform(ma);                              // see Transforms below
b->setTransform(mb);

if (a->collision(b))
{
    int ta, tb;
    a->getCollidingTriangles(ta, tb);             // triangle indices, in a and in b
    float point[3];
    a->getCollisionPoint(point, false);           // false: in world space
}

std::vector<std::pair<int, int>> pairs;
a->collisionAll(b, pairs);                        // every colliding pair (a's index, b's index)
```

Notes:

- Triangles are numbered in the order they were added.
- `collision()` stops at the first colliding pair. `collisionAll()` finds all of them.
- `collision(other, AccuracyDepth, MaxProcessingTime, other_transform)`:
  - `MaxProcessingTime` (in milliseconds) throws `TimeoutExpired` when exceeded.
  - `other_transform` places the other model for this test only, without changing its own
    transform. This is how to test a model against itself, or against a model shared by
    several objects.
  - `AccuracyDepth` is not supported.
- `getCollidingTriangles(t1, t2, ModelSpace)` returns the coordinates of the colliding
  triangles. With `ModelSpace` true they are in this model's space (the other triangle is
  transformed into it). Otherwise they are in world space.
- Calling `collision()` before `finalize()`, or adding triangles after it, throws
  `Inconsistency`.
- A model without triangles never collides.
- Collisions are between *surfaces*. A model entirely inside another one, without touching its
  triangles, does not collide.

### Transforms

`setTransform()` takes a 4x4 matrix of 16 floats in the *row vector* convention: a point
`p` is transformed as `p' = p * M`. So the upper 3x3 holds the rotation (transposed,
compared to the column vector convention), and the translation is in elements 12, 13 and 14:

```cpp
// rotation R (column vector convention) and translation t
float m[16] = { R00, R10, R20, 0,
                R01, R11, R21, 0,
                R02, R12, R22, 0,
                t0,  t1,  t2,  1 };
```

Transforms must be rigid: rotation and translation only, no scaling. Relative transforms
between models are computed in double precision, so models far from the origin keep their
precision relative to each other.

A static model (`newCollisionModel3D(true)`) caches the inverse of its transform. Use it for
models that move rarely.

### Ray and sphere queries

```cpp
const float origin[3] = { 0.2f, 0.2f, 5 }, dir[3] = { 0, 0, -1 };
if (a->rayCollision(origin, dir, true))          // true: the closest hit
    a->getCollisionPoint(p, false);

a->rayCollision(origin, dir, false, 0.0f, 3.0f); // a segment, from 0 to 3 along dir

const float center[3] = { 0.2f, 0.2f, 0.1f };
bool touching = a->sphereCollision(center, 0.5f);
```

- Rays and spheres are given in world space.
- The sphere is solid: a sphere that contains triangles collides with them.
- After a hit, `getCollidingTriangles()` and `getCollisionPoint()` describe it. Only the
  first triangle is valid.

There are also two utility functions:
- `SphereRayCollision(center, radius, origin, direction, point)` gives the first point where
  a ray meets a sphere. From inside the sphere, that's the exit point.
- `SphereSphereCollision(c1, r1, c2, r2, point)` tests two overlapping spheres. The point is
  on the first sphere, toward the second.

### Scenes with many objects

A `MultiObjectSystem` tracks many objects and finds the pairs that collide. Two
implementations share one interface:

- `newSweepPruneSystem(max_objects)`: an incremental sweep and prune broad phase, over the
  cubes around each object's bounding sphere. It holds up to `max_objects` objects. This is
  the one to use for anything but small scenes.
- `newSpheresSystem()`: tests all pairs of bounding spheres. Simple, with no size limit.

```cpp
#include <multiobject.h>
using namespace COLDET;

MultiObjectSystem* world = newSweepPruneSystem(1024);

Placement player;                                   // a TransformUpdater, see below
int p = world->add_object(mesh, &player);           // placed by an updater
const float crate_pos[3] = { 3, 0, 0 };
int c = world->add_object(mesh, crate_pos);         // placed at a position, same mesh
const float ball_pos[3] = { 0, 2, 0 };
int s = world->add_sphere_object(0.5f, ball_pos);   // no mesh: a sphere

player.matrix[12] = 2.5f;                           // moves the player
const float new_crate_pos[3] = { 3, 0.5f, 0 };
world->move_object(c, new_crate_pos);               // moves an object added by position

int n = world->check_collisions(true);              // true: exact, at triangle level
for (int i = 0; i < n; ++i)
{
    CollisionDetails d;
    world->get_collision(i, d);                     // d.id1, d.id2, d.point, d.t1, d.t2
}
```

A `TransformUpdater` gives the system an object's current transform:

```cpp
class Placement : public TransformUpdater
{
public:
    float matrix[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    const float* update() const override { return matrix; }
};
```

Behaviour:

- **Updaters**: an object added with an updater follows it. Each `check_collisions()` reads
  the updater's current transform, so `move_object()` isn't needed for those objects.
  Objects added at a position move with `move_object()`.
- **Shared meshes**: objects may share one `CollisionModel3D`, for example instances of one
  mesh. Each object is tested at its own placement. The system sets model transforms while
  checking, so don't rely on a shared model's transform afterwards.
- **Proximity versus exact**: `check_collisions(false)` returns the pairs whose bounding
  volumes overlap. Sweep and prune uses the cubes around the bounding spheres, and touching
  cubes count as overlapping. The sphere system uses the spheres themselves. With `true`,
  each candidate pair is tested exactly: mesh against mesh, mesh against sphere, or sphere
  against sphere. Both systems then give the same results.
- **Bounding radius**: an object's bounding radius is its model's `getRadius()`, the distance
  of the farthest vertex from the model's origin. It holds for any rotation.
- **Ids**: ids are reused after `remove_object()`. Adding an object to a full sweep and prune
  system returns -1. Invalid ids passed to `move_object()` and `remove_object()` are ignored.
- **Ownership**: models are not owned by the system.

### The triangle intersection test

The exact test under everything is available directly (`tritri.h`):

```cpp
#include <tritri.h>

bool hit = COLDET::tri_tri_intersect(a0, a1, a2, b0, b1, b2);

// testing the same triangles many times: compute their shapes once
COLDET::TriangleShape sa(a0, a1, a2), sb(b0, b1, b2);
hit = COLDET::tri_tri_intersect(a0, a1, a2, sa, b0, b1, b2, sb);
```

It is based on Tomas Möller's "A Fast Triangle-Triangle Intersection Test" (1997). It extends
the original with scale independent tolerances, robust coplanar handling, and degenerate
triangles (see [Behaviour and tolerances](#behaviour-and-tolerances)).

A `TriangleShape` holds the parts of the test that depend only on one triangle. Computing it
once saves about 15% per test. It is valid only for the exact vertices it was computed from, so
a transformed triangle needs its own.

### The OBB tree (experimental)

`obb/` contains an alternative hierarchy of PCA fitted oriented bounding boxes (with Eigen
types, in namespace `obb`). It works on the same kind of triangle soups, and tests the leaves
with the same triangle test.

```cpp
#include <obbtree.h>
using namespace obb;

OBBTree a, b;
a.add_triangle(Triangle{ Vector3f(0, 0, 0), Vector3f(1, 0, 0), Vector3f(0, 1, 0) }); // ...
a.build(4);                                         // up to 4 triangles per leaf
b.build(4);

a.set_dynamic_transform(Transform{ R, t });         // world placement: rotation, translation
b.set_dynamic_transform(Transform{ R2, t2 });
IntersectionResult first = a.check_intersection(b);        // first.hit, triangle_index1/2
IntersectionResult all = a.find_all_intersections(b);      // all.pairs
```

Compared with the original hierarchy, in the benchmarks below:

- It is faster on real surfaces, whose deep nodes are thin curved patches that oriented boxes
  fit tightly.
- It is slower to build.
- It is slightly slower on unstructured random triangle soups.

A leaf size of 4 was best in the leaf size sweep. Below that, the tree doesn't get any finer,
because the split rule stops splitting once all triangles fall on one side.

## Behaviour and tolerances

Collision results are reported for triangles that intersect *or touch*. Touching is decided
within a small tolerance:

- **Scale independent**: a vertex counts as on a triangle's plane within `1e-6` times the
  triangle's size. Models in millimetres and in kilometres behave the same.
- **Far from the origin**: the tolerance never drops below a few float ulps of the
  coordinates. Coplanar faces of models placed thousands of units from the origin are still
  found in contact. Gaps below float resolution there may go either way.
- **Coplanar triangles**: these are tested in 2D, with parallel and collinear edges handled
  explicitly.
- **Degenerate triangles**: segments, points and slivers are tested as segments, not projected.
- **Flat triangles in contact**: the box tests carry small epsilons, so boxes without
  thickness (single flat triangles) still meet under any rotation.

Pairs that barely touch can come out differently in different coordinate frames. COLDET
tests in one model's frame, the OBB tree in world space. The tests treat such pairs as
borderline (see below).

## Tests

```sh
cd build
ctest -C Release
```

There are about 120 tests, discovered from two executables:

**`test.obb`** (in `tests/obb/`):

| Suite | What |
|---|---|
| `tritri.*` | The triangle test against known answers, and against an independent double precision separating axis reference. The 35 hand-built cases (crossing, touching, mesh neighbours, coplanar, near coplanar, degenerate) run under 48 exact axis symmetries, scales from 2⁻¹⁰ to 2¹⁰, offsets and every vertex and argument order. Cases with clear margins also run under random rotations. |
| `tritri_equivalence.*` | The C++ triangle test against the frozen original C implementation (`tests/reference/tritri_reference.c`). Results must be bit for bit identical over about 29M pairs. |
| `scene.*` | Structured scenes through both libraries against brute force:<br>- a box resting on a floor (rotated, far from the origin)<br>- all 24 axis-aligned rotations at contact and at gaps<br>- translations, self tests, tiny against huge, mixed triangle sizes<br>- containment, empty, single triangle and clustered models<br>- degenerate and duplicate triangles |
| `coldet_*` | The rest of the model API:<br>- ray and sphere queries against double precision references<br>- collision details and points, `other_transform`, static models<br>- errors, `getRadius`, timeouts<br>- the sphere utility functions |
| `obb_*` | OBB fitting (including degenerate point sets), containment, box tests with parallel axes and flat boxes, tree diagnostics. |
| `compare.*` | COLDET against the OBB tree. A 1000-triangle random soup must match brute force exactly. Two 3.16M-triangle tori are timed, and both trees must agree. |

**`coldet_unit_tests`** (in `tests/unit/`):

| Suite | What |
|---|---|
| `sweep_prune.*` | One axis and whole boxes against brute force after every add, move and remove. Covers ties on a grid, zero length intervals, large moves and dense scenes. |
| `multi_object.*` | Both systems against brute force (exact and proximity) while objects move, leave and return. Covers shared meshes, rotations, updater driven moves, grids, capacity, points, empty models and invalid operations. |

Brute force comparisons are exact where both sides work in the same frame. Where they don't, a
differing pair is accepted only if it is borderline: its result changes when both
triangles shrink or grow by 1%, or the triangles touch within the input's float precision.

Disabled tests run with `--gtest_also_run_disabled_tests`:
- `compare.DISABLED_leaf_size_sweep` times the OBB tree at leaf sizes 1 to 32 on the torus
  test (about a minute).

> **Debug builds on Windows with vcpkg**: the build copies vcpkg's *release* `gtest.dll` next
> to the Debug test executables, which then crash at startup. To run Debug tests, copy
> `gtest.dll` and `gtest_main.dll` from `<vcpkg>/installed/x64-windows/debug/bin` next to them.

## Benchmarks

- `bm.triangles`: the triangle test on 500k pairs, with and without precomputed shapes.
- `bm.multiobject [move]`: sweep and prune against the sphere system, from 100 to 10,000
  objects at constant density. Each frame moves 10% of the objects by up to `move` per axis
  (default 0.5), then checks for collisions. Both systems must report the same contacts.
- `test.obb --gtest_filter=compare.*`: COLDET against the OBB tree, on a random soup and on two
  tori.

## Performance

Release build, single thread.

Two tori of 3.16M triangles each, triangles 1/1000 of the torus size, 100 random placements:

| | COLDET | OBB tree (leaf 4) |
|---|---|---|
| First hit | 67 µs | 13.5 µs |
| All pairs | 13.8 ms | 7.0 ms |
| Build (both models) | 2.5 s | 7.3 s |

Multiple objects, per frame, with moves of ±0.05 / ±0.5:

| Objects | Sphere system | Sweep and prune |
|---|---|---|
| 100 | 4.5 µs | 1.7 / 5.4 µs |
| 1,000 | 386 µs | 33 / 125 µs |
| 10,000 | 38 ms | 1.9 / 9.6 ms |

Sweep and prune costs grow with the number of endpoints each move passes. In a uniformly
filled 3D scene, endpoints per unit length on each axis grow as N^(2/3), so its cost is
superlinear, and small moves (coherent motion) are what keep it low.

## Repository layout

```
src/                 the coldet library
  coldet.h           public model API (CollisionModel3D, rays, spheres)
  multiobject.h      public multiple object API (MultiObjectSystem)
  tritri.h/.cpp      the triangle / triangle test
  box*.cpp/h         the bounding box hierarchy
  coldet*.cpp        model implementation
  sweep_prune.h      the sweep and prune broad phase
  multi_impl.h       the multiple object systems
obb/                 the experimental OBB tree library
tests/obb/           most tests (test.obb)
tests/unit/          sweep and prune and multiple object tests (coldet_unit_tests)
tests/reference/     the original C triangle test, frozen as a reference
benchmark/           benchmarks and the profiler they use
```

## Known limitations

- **Sweep and prune scaling**: endpoint passes grow with scene density, and each axis keeps an
  N×N bit matrix of overlaps (12.5 MB per axis at 10,000 objects). A spatial grid of
  sweeps, or a hash of overlapping pairs, would scale better.
- **Sweep and prune size**: endpoint ids are 16 bit, so `newSweepPruneSystem` supports up to
  65,536 objects. This limit is not tested.
- **`AccuracyDepth`** is not supported.
- **OBB tree build**: it is slower than the original hierarchy, and its split rule creates
  leaves larger than requested when triangles straddle the split plane.
- **Surfaces only**: models inside other models don't collide (see
  [Collision models](#collision-models)).

## License

GNU Lesser General Public License, version 2.1. See [LICENSE](LICENSE).

Copyright (C) 2000-2013 Amir Geva.
