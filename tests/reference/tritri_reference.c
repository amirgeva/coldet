/* Frozen copy of the original C implementation (src/tritri.c before its port to C++),
   kept only as the reference the C++ version (src/tritri.cpp) must match, bit for bit.
   Public names are prefixed with reference_. */

/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * int tri_tri_intersect(float V0[3],float V1[3],float V2[3],
 *                         float U0[3],float U1[3],float U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 */

#include <math.h>


/* if USE_EPSILON_TEST is true then we do a check: 
         if |dv|<eps then dv=0.0;
   else no check is done (which is less robust)
   The signed distances are scaled by |N| (twice the area of the triangle
   defining the plane), so eps is scaled to match: a vertex is treated as
   lying on the plane when it is within the larger of EPSILON times the
   triangle's longest edge, and COORD_EPSILON times the coordinate magnitude
   of the triangles (estimated from their first vertices).  The first keeps
   the test independent of the model's scale, the second covers the rounding
   error of coordinates far from the origin, which would otherwise break up
   coplanar triangles.
*/
#define USE_EPSILON_TEST TRUE  
#define EPSILON 0.000001
#define COORD_EPSILON 0.0000005 /* a few float ulps */


/* some macros */
#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2)          \
            dest[0]=v1[0]-v2[0]; \
            dest[1]=v1[1]-v2[1]; \
            dest[2]=v1[2]-v2[2]; 

/* sort so that a<=b */
#define SORT(a,b)       \
             if(a>b)    \
             {          \
               float c; \
               c=a;     \
               a=b;     \
               b=c;     \
             }

#define ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1) \
              isect0=VV0+(VV1-VV0)*D0/(D0-D1);    \
              isect1=VV0+(VV2-VV0)*D0/(D0-D2);


#define COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1) \
  if(D0D1>0.0f)                                         \
  {                                                     \
    /* here we know that D0D2<=0.0 */                   \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else if(D0D2>0.0f)                                    \
  {                                                     \
    /* here we know that d0d1<=0.0 */                   \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D1*D2>0.0f || D0!=0.0f)                       \
  {                                                     \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
    ISECT(VV0,VV1,VV2,D0,D1,D2,isect0,isect1);          \
  }                                                     \
  else if(D1!=0.0f)                                     \
  {                                                     \
    ISECT(VV1,VV0,VV2,D1,D0,D2,isect0,isect1);          \
  }                                                     \
  else if(D2!=0.0f)                                     \
  {                                                     \
    ISECT(VV2,VV0,VV1,D2,D0,D1,isect0,isect1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    /* triangles are coplanar */                        \
    return reference_coplanar_tri_tri(NC,V0,V1,V2,U0,U1,U2);      \
  }



/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */ 
/* 2D test of segment V0+t*(Ax,Ay) against U0U1, in the projection plane (i0,i1).
   For (nearly) parallel segments the determinants f, d, e are all rounding noise, so
   those are handled apart: parallel lines don't cross, collinear segments are
   compared by their extents along the line. */
static int edge_edge_test(const float V0[3],const float U0[3],const float U1[3],
                          float Ax,float Ay,short i0,short i1)
{
  float Bx,By,Cx,Cy,e,d,f,la,lb,lc,t0,t1,aa;
  Bx=U0[i0]-U1[i0];
  By=U0[i1]-U1[i1];
  Cx=V0[i0]-U0[i0];
  Cy=V0[i1]-U0[i1];
  f=Ay*Bx-Ax*By;
  d=By*Cx-Bx*Cy;
  la=Ax*Ax+Ay*Ay;
  lb=Bx*Bx+By*By;
  if(f*f<=(float)(EPSILON*EPSILON)*la*lb)
  {
    /* parallel: collinear when V0 is on the line of U0U1 */
    lc=Cx*Cx+Cy*Cy;
    if(d*d>(float)(EPSILON*EPSILON)*lb*(lc>la ? lc : la)) return 0;
    if(la==0.0f) return 0;
    /* overlap of U0U1 with V0V1, as parameters along V0V1 (0..1) */
    aa=1.0f/la;
    t0=-(Cx*Ax+Cy*Ay)*aa;                 /* U0 */
    t1=t0-(Bx*Ax+By*Ay)*aa;               /* U1 = U0 - B */
    if(t0>t1) { float t=t0; t0=t1; t1=t; }
    return t1>=0.0f && t0<=1.0f;
  }
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))
  {
    e=Ax*Cy-Ay*Cx;
    if(f>0)
    {
      if(e>=0 && e<=f) return 1;
    }
    else
    {
      if(e<=0 && e>=f) return 1;
    }
  }
  return 0;
}

#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  if(edge_edge_test(V0,U0,U1,Ax,Ay,i0,i1)) return 1;

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  float Ax,Ay;                                 \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  float a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}

int reference_coplanar_tri_tri(const float N[3],const float V0[3],const float V1[3],const float V2[3],
                     const float U0[3],const float U1[3],const float U2[3])
{
   float A[3];
   short i0,i1;
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
#pragma warning( disable : 4244 )
   A[0]=fabs(N[0]);
   A[1]=fabs(N[1]);
   A[2]=fabs(N[2]);
#pragma warning( default : 4244 )
   if(A[0]>A[1])
   {
      if(A[0]>A[2])  
      {
          i0=1;      /* A[0] is greatest */
          i1=2;
      }
      else
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
   }
   else   /* A[0]<=A[1] */
   {
      if(A[2]>A[1])
      {
          i0=0;      /* A[2] is greatest */
          i1=1;                                           
      }
      else
      {
          i0=0;      /* A[1] is greatest */
          i1=2;
      }
    }               
                
    /* test all edges of triangle 1 against the edges of triangle 2 */
    EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);
                
    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    POINT_IN_TRI(V0,U0,U1,U2);
    POINT_IN_TRI(U0,V0,V1,V2);

    return 0;
}


/* squared tolerance for the signed distances to the plane with normal N=E1xE2.
   coord_sq is the squared coordinate magnitude of the two triangles */
static float plane_epsilon_sq(const float N[3],const float E1[3],const float E2[3],float coord_sq)
{
  float l1,l2,lmax,c;
  /* the longer of two edges: the third is at most their sum */
  l1=DOT(E1,E1);
  l2=DOT(E2,E2);
  lmax=l1;
  if(l2>lmax) lmax=l2;
  lmax*=(float)(EPSILON*EPSILON);
  c=(float)(COORD_EPSILON*COORD_EPSILON)*coord_sq;
  if(c>lmax) lmax=c;
  return DOT(N,N)*lmax;
}

/* ---- degenerate triangles (segments and points) ----
   The general test relies on both triangles' planes.  A degenerate triangle has no
   plane: its edges are tested as segments against the other triangle instead. */

/* the sine of the angle at the first vertex is below EPSILON (or an edge has no length) */
static int is_degenerate(const float N[3],const float E1[3],const float E2[3])
{
  return DOT(N,N) <= (float)(EPSILON*EPSILON)*DOT(E1,E1)*DOT(E2,E2);
}

/* inclusive test of a point X in the plane of triangle U (normal N) */
static int point_in_tri_3d(const float X[3],const float U0[3],const float U1[3],const float U2[3],
                           const float N[3])
{
  const float* U[4];
  float E[3],W[3],C[3],d;
  int i;
  U[0]=U0; U[1]=U1; U[2]=U2; U[3]=U0;
  for(i=0;i<3;++i)
  {
    SUB(E,U[i+1],U[i]);
    SUB(W,X,U[i]);
    CROSS(C,E,W);
    d=DOT(C,N);
    if(d<0.0f && d*d>(float)(EPSILON*EPSILON)*DOT(C,C)*DOT(N,N)) return 0;
  }
  return 1;
}

/* segment PQ against the non degenerate triangle U with normal N.
   eps is the squared tolerance of the plane distances (see plane_epsilon_sq) */
static int segment_tri(const float P[3],const float Q[3],
                       const float U0[3],const float U1[3],const float U2[3],
                       const float N[3],float eps)
{
  float T[3],X[3],dp,dq,s;
  SUB(T,P,U0); dp=DOT(N,T);
  SUB(T,Q,U0); dq=DOT(N,T);
  if(dp*dp<eps) dp=0.0f;
  if(dq*dq<eps) dq=0.0f;
  if(dp*dq>0.0f) return 0;            /* both ends on the same side */
  if(dp==0.0f && dq==0.0f)            /* the segment lies in the plane */
  {
    /* an end in the triangle (inclusive, unlike the coplanar test), or crossing its edges */
    if(point_in_tri_3d(P,U0,U1,U2,N) || point_in_tri_3d(Q,U0,U1,U2,N)) return 1;
    return reference_coplanar_tri_tri(N,P,Q,Q,U0,U1,U2);
  }
  s=dp/(dp-dq);
  X[0]=P[0]+s*(Q[0]-P[0]);
  X[1]=P[1]+s*(Q[1]-P[1]);
  X[2]=P[2]+s*(Q[2]-P[2]);
  return point_in_tri_3d(X,U0,U1,U2,N);
}

static float clamp01(float x) { return x<0.0f ? 0.0f : (x>1.0f ? 1.0f : x); }

/* squared distance between segments P1Q1 and P2Q2 (Ericson, Real-Time Collision Detection 5.1.9) */
static float segment_segment_dist_sq(const float P1[3],const float Q1[3],const float P2[3],const float Q2[3])
{
  float d1[3],d2[3],r[3],c1[3],c2[3],w[3];
  float a,e,f,b,c,denom,s,t;
  const float tiny=1e-30f;
  SUB(d1,Q1,P1); SUB(d2,Q2,P2); SUB(r,P1,P2);
  a=DOT(d1,d1); e=DOT(d2,d2); f=DOT(d2,r);
  if(a<=tiny && e<=tiny) return DOT(r,r);
  if(a<=tiny) { s=0.0f; t=clamp01(f/e); }
  else
  {
    c=DOT(d1,r);
    if(e<=tiny) { t=0.0f; s=clamp01(-c/a); }
    else
    {
      b=DOT(d1,d2);
      denom=a*e-b*b;
      s = denom>0.0f ? clamp01((b*f-c*e)/denom) : 0.0f;
      t=(b*s+f)/e;
      if(t<0.0f) { t=0.0f; s=clamp01(-c/a); }
      else if(t>1.0f) { t=1.0f; s=clamp01((b-c)/a); }
    }
  }
  c1[0]=P1[0]+d1[0]*s; c1[1]=P1[1]+d1[1]*s; c1[2]=P1[2]+d1[2]*s;
  c2[0]=P2[0]+d2[0]*t; c2[1]=P2[1]+d2[1]*t; c2[2]=P2[2]+d2[2]*t;
  SUB(w,c1,c2);
  return DOT(w,w);
}

/* the longest edge of a triangle, as a segment */
static float longest_edge(const float* T[3],const float** P,const float** Q)
{
  float E[3],l,best=-1.0f;
  int i;
  for(i=0;i<3;++i)
  {
    SUB(E,T[(i+1)%3],T[i]);
    l=DOT(E,E);
    if(l>best) { best=l; *P=T[i]; *Q=T[(i+1)%3]; }
  }
  return best;
}

/* D is degenerate, W may or may not be */
static int degenerate_tri_tri(const float D0[3],const float D1[3],const float D2[3],
                              const float W0[3],const float W1[3],const float W2[3],
                              float coord_sq)
{
  float E1[3],E2[3],N[3],eps;
  SUB(E1,W1,W0);
  SUB(E2,W2,W0);
  CROSS(N,E1,E2);
  if(is_degenerate(N,E1,E2))
  {
    /* both are segments (or points): compare their distance to the tolerance */
    const float *TD[3],*TW[3],*P1,*Q1,*P2,*Q2;
    float l1,l2,tol_sq;
    TD[0]=D0; TD[1]=D1; TD[2]=D2;
    TW[0]=W0; TW[1]=W1; TW[2]=W2;
    l1=longest_edge(TD,&P1,&Q1);
    l2=longest_edge(TW,&P2,&Q2);
    tol_sq=(float)(EPSILON*EPSILON)*(l1>l2 ? l1 : l2);
    if((float)(COORD_EPSILON*COORD_EPSILON)*coord_sq>tol_sq) tol_sq=(float)(COORD_EPSILON*COORD_EPSILON)*coord_sq;
    return segment_segment_dist_sq(P1,Q1,P2,Q2)<=tol_sq;
  }
  eps=plane_epsilon_sq(N,E1,E2,coord_sq);
  return segment_tri(D0,D1,W0,W1,W2,N,eps) ||
         segment_tri(D1,D2,W0,W1,W2,N,eps) ||
         segment_tri(D2,D0,W0,W1,W2,N,eps);
}

int reference_tri_tri_intersect(const float V0[3],const float V1[3],const float V2[3],
                      const float U0[3],const float U1[3],const float U2[3])
{
  float E1[3],E2[3],T[3];
  float N1[3],N2[3],eps; /* squared tolerance */
  const float* NC; /* normal used for the coplanar test */
  float coord_sq;
  float du0,du1,du2,dv0,dv1,dv2;
  float D[3];
  float isect1[2], isect2[2];
  float du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  float vp0,vp1,vp2;
  float up0,up1,up2;
  float b,c,max;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  /* plane equation 1: N1.(X-V0)=0 */

  coord_sq=DOT(V0,V0);
  if(DOT(U0,U0)>coord_sq) coord_sq=DOT(U0,U0);
  if(is_degenerate(N1,E1,E2))
    return degenerate_tri_tri(V0,V1,V2,U0,U1,U2,coord_sq);

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane.
     Measuring relative to V0 avoids cancellation when coordinates are large */
  SUB(T,U0,V0); du0=DOT(N1,T);
  SUB(T,U1,V0); du1=DOT(N1,T);
  SUB(T,U2,V0); du2=DOT(N1,T);

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  eps=plane_epsilon_sq(N1,E1,E2,coord_sq);
  if(du0*du0<eps) du0=0.0;
  if(du1*du1<eps) du1=0.0;
  if(du2*du2<eps) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  /* plane equation 2: N2.(X-U0)=0 */
  if(is_degenerate(N2,E1,E2))
    return degenerate_tri_tri(U0,U1,U2,V0,V1,V2,coord_sq);

  /* put V0,V1,V2 into plane equation 2 */
  SUB(T,V0,U0); dv0=DOT(N2,T);
  SUB(T,V1,U0); dv1=DOT(N2,T);
  SUB(T,V2,U0); dv2=DOT(N2,T);

#if USE_EPSILON_TEST==TRUE
  eps=plane_epsilon_sq(N2,E1,E2,coord_sq);
  if(dv0*dv0<eps) dv0=0.0;
  if(dv1*dv1<eps) dv1=0.0;
  if(dv2*dv2<eps) dv2=0.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;
        
  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* a degenerate triangle has no normal: project the coplanar test with the other one */
  NC = DOT(N1,N1)>=DOT(N2,N2) ? N1 : N2;

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
#pragma warning( disable : 4244 )
  max=fabs(D[0]);
  index=0;
  b=fabs(D[1]);
  c=fabs(D[2]);
#pragma warning( default : 4244 )
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

        /* this is the simplified projection onto L*/
        vp0=V0[index];
        vp1=V1[index];
        vp2=V2[index];

        up0=U0[index];
        up1=U1[index];
        up2=U2[index];

  /* compute interval for triangle 1 */
  COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,isect1[0],isect1[1]);

  /* compute interval for triangle 2 */
  COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,isect2[0],isect2[1]);

  SORT(isect1[0],isect1[1]);
  SORT(isect2[0],isect2[1]);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
  return 1;
}

