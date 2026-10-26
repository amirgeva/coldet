/** \file tritri.h
    Triangle / triangle intersection test.

    Based on Tomas Moller, "A Fast Triangle-Triangle Intersection Test",
    Journal of Graphics Tools, 2(2), 1997, with:
    - scale independent tolerances, also covering coordinates far from the origin
    - robust handling of coplanar triangles, including parallel and collinear edges
    - degenerate triangles (segments and points), tested as segments
*/
#pragma once

namespace COLDET
{

/** What the test derives from a triangle's vertices alone.  Computing it once and reusing
    it for every test of the triangle saves that work per pair.  It is only valid for the
    exact vertices it was computed from: a transformed triangle needs its own. */
struct TriangleShape
{
  TriangleShape() = default;
  TriangleShape(const float v0[3], const float v1[3], const float v2[3]);

  /** |N|^2, for the normal N = (v1-v0) x (v2-v0) */
  float normal_sq = 0.0f;
  /** EPSILON^2 times the longer of |v1-v0|^2 and |v2-v0|^2: the size part of the tolerance */
  float size_eps_sq = 0.0f;
  /** The triangle is a segment or a point (or a sliver below EPSILON) */
  bool degenerate = false;
};

/** True if the triangles intersect or touch */
bool tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3],
                       const float U0[3], const float U1[3], const float U2[3]);

/** The same, with both triangles' shapes computed in advance */
bool tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3], const TriangleShape& v,
                       const float U0[3], const float U1[3], const float U2[3], const TriangleShape& u);

/** The same, with the first triangle's shape computed in advance */
bool tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3], const TriangleShape& v,
                       const float U0[3], const float U1[3], const float U2[3]);

} // namespace COLDET
