/*  Triangle / triangle intersection test.
    Based on Tomas Moller, "A Fast Triangle-Triangle Intersection Test",
    Journal of Graphics Tools, 2(2), 1997.

    Note: every expression keeps the operand order of the original C implementation, so
    that results are identical to it, bit for bit.  Keep it that way when editing:
    floating point operations are not associative.
*/
#include <cmath>
#include <utility>
#include "tritri.h"

namespace COLDET
{

namespace
{

/* Tolerances.  The signed distances to a triangle's plane are scaled by |N| (twice the
   triangle's area), so the tolerance is scaled to match: a vertex is treated as lying on
   the plane when it is within the larger of EPSILON times the triangle's longest edge,
   and COORD_EPSILON times the coordinate magnitude of the triangles (estimated from their
   first vertices).  The first keeps the test independent of the model's scale, the second
   covers the rounding error of coordinates far from the origin, which would otherwise
   break up coplanar triangles. */
constexpr double EPSILON = 0.000001;
constexpr double COORD_EPSILON = 0.0000005; // a few float ulps
constexpr float EPSILON_SQ = float(EPSILON * EPSILON);
constexpr float COORD_EPSILON_SQ = float(COORD_EPSILON * COORD_EPSILON);

struct Vec
{
  float v[3];
  operator const float*() const { return v; }
};

inline Vec sub(const float* a, const float* b) { return { { a[0] - b[0], a[1] - b[1], a[2] - b[2] } }; }

inline Vec cross(const float* a, const float* b)
{
  return { { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] } };
}

inline float dot(const float* a, const float* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

/* Squared tolerance for the signed distances to the plane of a triangle.
   coord_sq is the squared coordinate magnitude of the two triangles */
inline float plane_epsilon_sq(const TriangleShape& s, float coord_sq)
{
  float lmax = s.size_eps_sq;
  const float c = COORD_EPSILON_SQ * coord_sq;
  if (c > lmax) lmax = c;
  return s.normal_sq * lmax;
}

/* ---- Coplanar triangles ---- */

/* 2D test of segment V0+t*(Ax,Ay) against U0U1, in the projection plane (i0,i1).
   Based on Franklin Antonio's gem: "Faster Line Segment Intersection", Graphics Gems III,
   pp. 199-202.  For (nearly) parallel segments the determinants f, d, e are all rounding
   noise, so those are handled apart: parallel lines don't cross, collinear segments are
   compared by their extents along the line. */
bool edge_edge_test(const float* V0, const float* U0, const float* U1, float Ax, float Ay, int i0, int i1)
{
  const float Bx = U0[i0] - U1[i0];
  const float By = U0[i1] - U1[i1];
  const float Cx = V0[i0] - U0[i0];
  const float Cy = V0[i1] - U0[i1];
  const float f = Ay * Bx - Ax * By;
  const float d = By * Cx - Bx * Cy;
  const float la = Ax * Ax + Ay * Ay;
  const float lb = Bx * Bx + By * By;
  if (f * f <= EPSILON_SQ * la * lb)
  {
    // parallel: collinear when V0 is on the line of U0U1
    const float lc = Cx * Cx + Cy * Cy;
    if (d * d > EPSILON_SQ * lb * (lc > la ? lc : la)) return false;
    if (la == 0.0f) return false;
    // overlap of U0U1 with V0V1, as parameters along V0V1 (0..1)
    const float aa = 1.0f / la;
    float t0 = -(Cx * Ax + Cy * Ay) * aa;   // U0
    float t1 = t0 - (Bx * Ax + By * Ay) * aa; // U1 = U0 - B
    if (t0 > t1) std::swap(t0, t1);
    return t1 >= 0.0f && t0 <= 1.0f;
  }
  if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f))
  {
    const float e = Ax * Cy - Ay * Cx;
    if (f > 0)
    {
      if (e >= 0 && e <= f) return true;
    }
    else
    {
      if (e <= 0 && e >= f) return true;
    }
  }
  return false;
}

/* edge V0V1 against the three edges of U */
bool edge_against_tri_edges(const float* V0, const float* V1, const float* U0, const float* U1, const float* U2,
                            int i0, int i1)
{
  const float Ax = V1[i0] - V0[i0];
  const float Ay = V1[i1] - V0[i1];
  return edge_edge_test(V0, U0, U1, Ax, Ay, i0, i1) ||
         edge_edge_test(V0, U1, U2, Ax, Ay, i0, i1) ||
         edge_edge_test(V0, U2, U0, Ax, Ay, i0, i1);
}

/* is V0 strictly inside triangle U, in the projection plane (i0,i1) */
bool point_in_tri(const float* V0, const float* U0, const float* U1, const float* U2, int i0, int i1)
{
  float a = U1[i1] - U0[i1];
  float b = -(U1[i0] - U0[i0]);
  float c = -a * U0[i0] - b * U0[i1];
  const float d0 = a * V0[i0] + b * V0[i1] + c;

  a = U2[i1] - U1[i1];
  b = -(U2[i0] - U1[i0]);
  c = -a * U1[i0] - b * U1[i1];
  const float d1 = a * V0[i0] + b * V0[i1] + c;

  a = U0[i1] - U2[i1];
  b = -(U0[i0] - U2[i0]);
  c = -a * U2[i0] - b * U2[i1];
  const float d2 = a * V0[i0] + b * V0[i1] + c;
  return d0 * d1 > 0.0 && d0 * d2 > 0.0;
}

bool coplanar_tri_tri(const float* N, const float* V0, const float* V1, const float* V2,
                      const float* U0, const float* U1, const float* U2)
{
  // project onto the axis-aligned plane that maximizes the area of the triangles
  const float A[3] = { std::fabs(N[0]), std::fabs(N[1]), std::fabs(N[2]) };
  int i0, i1;
  if (A[0] > A[1])
  {
    if (A[0] > A[2]) { i0 = 1; i1 = 2; } // A[0] is greatest
    else             { i0 = 0; i1 = 1; } // A[2] is greatest
  }
  else
  {
    if (A[2] > A[1]) { i0 = 0; i1 = 1; } // A[2] is greatest
    else             { i0 = 0; i1 = 2; } // A[1] is greatest
  }

  // test all edges of triangle 1 against the edges of triangle 2
  if (edge_against_tri_edges(V0, V1, U0, U1, U2, i0, i1)) return true;
  if (edge_against_tri_edges(V1, V2, U0, U1, U2, i0, i1)) return true;
  if (edge_against_tri_edges(V2, V0, U0, U1, U2, i0, i1)) return true;

  // finally, test if one triangle is totally contained in the other
  return point_in_tri(V0, U0, U1, U2, i0, i1) || point_in_tri(U0, V0, V1, V2, i0, i1);
}

/* ---- Degenerate triangles (segments and points) ----
   The general test relies on both triangles' planes.  A degenerate triangle has no
   plane: its edges are tested as segments against the other triangle instead. */

/* inclusive test of a point X in the plane of triangle U (normal N) */
bool point_in_tri_3d(const float* X, const float* U0, const float* U1, const float* U2,
                     const float* N, float normal_sq)
{
  const float* U[4] = { U0, U1, U2, U0 };
  for (int i = 0; i < 3; ++i)
  {
    const Vec E = sub(U[i + 1], U[i]);
    const Vec W = sub(X, U[i]);
    const Vec C = cross(E, W);
    const float d = dot(C, N);
    if (d < 0.0f && d * d > EPSILON_SQ * dot(C, C) * normal_sq) return false;
  }
  return true;
}

/* segment PQ against the non degenerate triangle U with normal N.
   eps is the squared tolerance of the plane distances (see plane_epsilon_sq) */
bool segment_tri(const float* P, const float* Q, const float* U0, const float* U1, const float* U2,
                 const float* N, float normal_sq, float eps)
{
  float dp = dot(N, sub(P, U0));
  float dq = dot(N, sub(Q, U0));
  if (dp * dp < eps) dp = 0.0f;
  if (dq * dq < eps) dq = 0.0f;
  if (dp * dq > 0.0f) return false; // both ends on the same side
  if (dp == 0.0f && dq == 0.0f)     // the segment lies in the plane
  {
    // an end in the triangle (inclusive, unlike the coplanar test), or crossing its edges
    if (point_in_tri_3d(P, U0, U1, U2, N, normal_sq) || point_in_tri_3d(Q, U0, U1, U2, N, normal_sq)) return true;
    return coplanar_tri_tri(N, P, Q, Q, U0, U1, U2);
  }
  const float s = dp / (dp - dq);
  const Vec X = { { P[0] + s * (Q[0] - P[0]), P[1] + s * (Q[1] - P[1]), P[2] + s * (Q[2] - P[2]) } };
  return point_in_tri_3d(X, U0, U1, U2, N, normal_sq);
}

inline float clamp01(float x) { return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x); }

/* squared distance between segments P1Q1 and P2Q2 (Ericson, Real-Time Collision Detection 5.1.9) */
float segment_segment_dist_sq(const float* P1, const float* Q1, const float* P2, const float* Q2)
{
  const float tiny = 1e-30f;
  const Vec d1 = sub(Q1, P1), d2 = sub(Q2, P2), r = sub(P1, P2);
  const float a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
  float s, t;
  if (a <= tiny && e <= tiny) return dot(r, r);
  if (a <= tiny) { s = 0.0f; t = clamp01(f / e); }
  else
  {
    const float c = dot(d1, r);
    if (e <= tiny) { t = 0.0f; s = clamp01(-c / a); }
    else
    {
      const float b = dot(d1, d2);
      const float denom = a * e - b * b;
      s = denom > 0.0f ? clamp01((b * f - c * e) / denom) : 0.0f;
      t = (b * s + f) / e;
      if (t < 0.0f) { t = 0.0f; s = clamp01(-c / a); }
      else if (t > 1.0f) { t = 1.0f; s = clamp01((b - c) / a); }
    }
  }
  const Vec c1 = { { P1[0] + d1[0] * s, P1[1] + d1[1] * s, P1[2] + d1[2] * s } };
  const Vec c2 = { { P2[0] + d2[0] * t, P2[1] + d2[1] * t, P2[2] + d2[2] * t } };
  const Vec w = sub(c1, c2);
  return dot(w, w);
}

/* the longest edge of a triangle, as a segment.  Returns its squared length */
float longest_edge(const float* const T[3], const float*& P, const float*& Q)
{
  float best = -1.0f;
  for (int i = 0; i < 3; ++i)
  {
    const float l = dot(sub(T[(i + 1) % 3], T[i]), sub(T[(i + 1) % 3], T[i]));
    if (l > best) { best = l; P = T[i]; Q = T[(i + 1) % 3]; }
  }
  return best;
}

/* D is degenerate, W (with shape w) may or may not be */
bool degenerate_tri_tri(const float* D0, const float* D1, const float* D2,
                        const float* W0, const float* W1, const float* W2, const TriangleShape& w,
                        float coord_sq)
{
  if (w.degenerate)
  {
    // both are segments (or points): compare their distance to the tolerance
    const float* TD[3] = { D0, D1, D2 };
    const float* TW[3] = { W0, W1, W2 };
    const float *P1, *Q1, *P2, *Q2;
    const float l1 = longest_edge(TD, P1, Q1);
    const float l2 = longest_edge(TW, P2, Q2);
    float tol_sq = EPSILON_SQ * (l1 > l2 ? l1 : l2);
    if (COORD_EPSILON_SQ * coord_sq > tol_sq) tol_sq = COORD_EPSILON_SQ * coord_sq;
    return segment_segment_dist_sq(P1, Q1, P2, Q2) <= tol_sq;
  }
  const Vec N = cross(sub(W1, W0), sub(W2, W0));
  const float eps = plane_epsilon_sq(w, coord_sq);
  return segment_tri(D0, D1, W0, W1, W2, N, w.normal_sq, eps) ||
         segment_tri(D1, D2, W0, W1, W2, N, w.normal_sq, eps) ||
         segment_tri(D2, D0, W0, W1, W2, N, w.normal_sq, eps);
}

/* ---- General case ---- */

/* the two ends of a triangle's interval on the intersection line */
inline void isect(float vv0, float vv1, float vv2, float d0, float d1, float d2, float& isect0, float& isect1)
{
  isect0 = vv0 + (vv1 - vv0) * d0 / (d0 - d1);
  isect1 = vv0 + (vv2 - vv0) * d0 / (d0 - d2);
}

/* A triangle's interval on the intersection line, from its vertices' projections p0..p2
   and signed distances d0..d2 to the other triangle's plane.  False if coplanar. */
bool compute_intervals(float p0, float p1, float p2, float d0, float d1, float d2, float d0d1, float d0d2,
                       float& isect0, float& isect1)
{
  if (d0d1 > 0.0f)
    isect(p2, p0, p1, d2, d0, d1, isect0, isect1); // d0, d1 on the same side, d2 on the other or on the plane
  else if (d0d2 > 0.0f)
    isect(p1, p0, p2, d1, d0, d2, isect0, isect1);
  else if (d1 * d2 > 0.0f || d0 != 0.0f)
    isect(p0, p1, p2, d0, d1, d2, isect0, isect1);
  else if (d1 != 0.0f)
    isect(p1, p0, p2, d1, d0, d2, isect0, isect1);
  else if (d2 != 0.0f)
    isect(p2, p0, p1, d2, d0, d1, isect0, isect1);
  else
    return false; // coplanar
  return true;
}

/* the test, with u computed only when needed (when V's plane doesn't reject U), if not given */
bool intersect(const float* V0, const float* V1, const float* V2, const TriangleShape& v,
               const float* U0, const float* U1, const float* U2, const TriangleShape* u_given)
{
  const Vec N1 = cross(sub(V1, V0), sub(V2, V0)); // plane of V: N1.(X-V0)=0

  float coord_sq = dot(V0, V0);
  if (dot(U0, U0) > coord_sq) coord_sq = dot(U0, U0);
  TriangleShape u_computed;
  if (v.degenerate)
  {
    const TriangleShape& u = u_given ? *u_given : (u_computed = TriangleShape(U0, U1, U2));
    return degenerate_tri_tri(V0, V1, V2, U0, U1, U2, u, coord_sq);
  }

  // signed distances of U to the plane of V.  Measuring relative to V0 avoids cancellation
  // when coordinates are large
  float du0 = dot(N1, sub(U0, V0));
  float du1 = dot(N1, sub(U1, V0));
  float du2 = dot(N1, sub(U2, V0));
  float eps = plane_epsilon_sq(v, coord_sq);
  if (du0 * du0 < eps) du0 = 0.0f;
  if (du1 * du1 < eps) du1 = 0.0f;
  if (du2 * du2 < eps) du2 = 0.0f;
  const float du0du1 = du0 * du1;
  const float du0du2 = du0 * du2;
  if (du0du1 > 0.0f && du0du2 > 0.0f) return false; // all on the same side, none on the plane

  const TriangleShape& u = u_given ? *u_given : (u_computed = TriangleShape(U0, U1, U2));
  const Vec N2 = cross(sub(U1, U0), sub(U2, U0)); // plane of U: N2.(X-U0)=0
  if (u.degenerate)
    return degenerate_tri_tri(U0, U1, U2, V0, V1, V2, v, coord_sq);

  // signed distances of V to the plane of U
  float dv0 = dot(N2, sub(V0, U0));
  float dv1 = dot(N2, sub(V1, U0));
  float dv2 = dot(N2, sub(V2, U0));
  eps = plane_epsilon_sq(u, coord_sq);
  if (dv0 * dv0 < eps) dv0 = 0.0f;
  if (dv1 * dv1 < eps) dv1 = 0.0f;
  if (dv2 * dv2 < eps) dv2 = 0.0f;
  const float dv0dv1 = dv0 * dv1;
  const float dv0dv2 = dv0 * dv2;
  if (dv0dv1 > 0.0f && dv0dv2 > 0.0f) return false;

  // a nearly degenerate triangle has a poor normal: project a coplanar test with the other one
  const float* NC = v.normal_sq >= u.normal_sq ? N1 : N2;

  // direction of the intersection line, and the index of its largest component
  const Vec D = cross(N1, N2);
  float max = std::fabs(D[0]);
  int index = 0;
  const float b = std::fabs(D[1]);
  const float c = std::fabs(D[2]);
  if (b > max) { max = b; index = 1; }
  if (c > max) { max = c; index = 2; }

  // the simplified projection onto the line, and each triangle's interval on it
  float isect1[2], isect2[2];
  if (!compute_intervals(V0[index], V1[index], V2[index], dv0, dv1, dv2, dv0dv1, dv0dv2, isect1[0], isect1[1]))
    return coplanar_tri_tri(NC, V0, V1, V2, U0, U1, U2);
  if (!compute_intervals(U0[index], U1[index], U2[index], du0, du1, du2, du0du1, du0du2, isect2[0], isect2[1]))
    return coplanar_tri_tri(NC, V0, V1, V2, U0, U1, U2);

  if (isect1[0] > isect1[1]) std::swap(isect1[0], isect1[1]);
  if (isect2[0] > isect2[1]) std::swap(isect2[0], isect2[1]);
  return !(isect1[1] < isect2[0] || isect2[1] < isect1[0]);
}

} // namespace

TriangleShape::TriangleShape(const float v0[3], const float v1[3], const float v2[3])
{
  const Vec E1 = sub(v1, v0), E2 = sub(v2, v0);
  const Vec N = cross(E1, E2);
  const float l1 = dot(E1, E1), l2 = dot(E2, E2);
  normal_sq = dot(N, N);
  // the longer of two edges: the third is at most their sum
  size_eps_sq = (l2 > l1 ? l2 : l1) * EPSILON_SQ;
  // the sine of the angle at v0 is below EPSILON (or an edge has no length)
  degenerate = normal_sq <= EPSILON_SQ * l1 * l2;
}

bool tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3],
                       const float U0[3], const float U1[3], const float U2[3])
{
  return intersect(V0, V1, V2, TriangleShape(V0, V1, V2), U0, U1, U2, nullptr);
}

bool tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3], const TriangleShape& v,
                       const float U0[3], const float U1[3], const float U2[3], const TriangleShape& u)
{
  return intersect(V0, V1, V2, v, U0, U1, U2, &u);
}

bool tri_tri_intersect(const float V0[3], const float V1[3], const float V2[3], const TriangleShape& v,
                       const float U0[3], const float U1[3], const float U2[3])
{
  return intersect(V0, V1, V2, v, U0, U1, U2, nullptr);
}

} // namespace COLDET
