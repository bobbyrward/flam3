/*
    FLAM3 - cosmic recursive fractal flames
    Copyright (C) 1992-2006  Scott Draves <source@flam3.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/

#define _GNU_SOURCE

#include "private.h"
#include "img.h"
#include "config.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <errno.h>

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/mach_error.h>
#define flam3_os "OSX"
#else
#ifdef WIN32
#define WINVER 0x0500
#include <windows.h>
#define flam3_os "WIN"
#else
#define flam3_os "LNX"
#endif
#endif



char *flam3_version() {

  if (strcmp(SVN_REV, "exported"))
    return flam3_os "-" VERSION "." SVN_REV;
  return VERSION;
}


#define SUB_BATCH_SIZE     10000
#define CHOOSE_XFORM_GRAIN 10000

#define random_distrib(v) ((v)[random()%vlen(v)])

#define badvalue(x) (((x)!=(x))||((x)>1e10)||((x)<-1e10))

extern void sincos(double x, double *s, double *c);

typedef struct {

   double tx,ty; /* Starting coordinates */

   double precalc_atan, precalc_sina;  /* Precalculated, if needed */
   double precalc_cosa, precalc_sqrt;
   double precalc_sumsq,precalc_atanyx;

   flam3_xform *xform; /* For the important values */

   /* Output Coords */

   double p0, p1;

   /* Pointer to the isaac RNG state */
   randctx *rc;

} flam3_iter_helper;



/* Variation functions */
static void var0_linear(void *, double);
static void var1_sinusoidal(void *, double);
static void var2_spherical(void *, double);
static void var3_swirl(void *, double);
static void var4_horseshoe(void *, double);
static void var5_polar(void *, double);
static void var6_handkerchief(void *, double);
static void var7_heart(void *, double);
static void var8_disc(void *, double);
static void var9_spiral(void *, double);
static void var10_hyperbolic(void *, double);
static void var11_diamond(void *, double);
static void var12_ex(void *, double);
static void var13_julia(void *, double);
static void var14_bent(void *, double);
static void var15_waves(void *, double);
static void var16_fisheye(void *, double);
static void var17_popcorn(void *, double);
static void var18_exponential(void *, double);
static void var19_power(void *, double);
static void var20_cosine(void *, double);
static void var21_rings(void *, double);
static void var22_fan(void *, double);
static void var23_blob(void *, double);
static void var24_pdj(void *, double);
static void var25_fan2(void *, double);
static void var26_rings2(void *, double);
static void var27_eyefish(void *, double);
static void var28_bubble(void *, double);
static void var29_cylinder(void *, double);
static void var30_perspective(void *, double);
static void var31_noise(void *, double);
static void var32_juliaN_generic(void *, double);
static void var33_juliaScope_generic(void *, double);
static void var34_blur(void *, double);
static void var35_gaussian(void *, double);
static void var36_radial_blur(void *, double);
static void var37_pie(void *, double);
static void var38_ngon(void *, double);
/*static void var39_image(void *, double);*/
static void var39_curl(void *, double);
static void var40_rectangles(void *, double);
static void var41_arch(void *helper, double weight);
static void var42_tangent(void *helper, double weight);
static void var43_square(void *helper, double weight);
static void var44_rays(void *helper, double weight);
static void var45_blade(void *helper, double weight);
static void var46_secant2(void *helper, double weight);
static void var47_twintrian(void *helper, double weight);
static void var48_cross(void *helper, double weight);
static void var49_disc2(void *helper, double weight);
static void var50_supershape(void *helper, double weight);
static void var51_flower(void *helper, double weight);
static void var52_conic(void *helper, double weight);
static void var53_parabola(void *helper, double weight);
static void var54_bent2(void *helper, double weight);
static void var55_bipolar(void *helper, double weight);
static void var56_boarders(void *helper, double weight);
static void var57_butterfly(void *helper, double weight);
static void var58_cell(void *helper, double weight);
static void var59_cpow(void *helper, double weight);
static void var60_curve(void *helper, double weight);
static void var61_edisc(void *helper, double weight);
static void var62_elliptic(void *helper, double weight);
static void var63_escher(void *helper, double weight);
static void var64_foci(void *helper, double weight);
static void var65_lazysusan(void *helper, double weight);
static void var66_loonie(void *helper, double weight);
static void var67_pre_blur(void *helper, double weight);
static void var68_modulus(void *helper, double weight);
static void var69_oscope(void *helper, double weight);
static void var70_polar2(void *helper, double weight);
static void var71_popcorn2(void *helper, double weight);
static void var72_scry(void *helper, double weight);
static void var73_separation(void *helper, double weight);
static void var74_split(void *helper, double weight);
static void var75_splits(void *helper, double weight);
static void var76_stripes(void *helper, double weight);
static void var77_wedge(void *helper, double weight);
static void var78_wedge_julia(void *helper, double weight);
static void var79_wedge_sph(void *helper, double weight);
static void var80_whorl(void *helper, double weight);
static void var81_waves2(void *helper, double weight);

/* Precalculation functions */
static void perspective_precalc(flam3_xform *xf);
static void juliaN_precalc(flam3_xform *xf);
static void juliaScope_precalc(flam3_xform *xf);
static void radial_blur_precalc(flam3_xform *xf);
static void waves_precalc(flam3_xform *xf);
static void disc2_precalc(flam3_xform *xf);
static void supershape_precalc(flam3_xform *xf);
static void wedgeJulia_precalc(flam3_xform *xf);

static int id_matrix(double s[3][2]);
static double flam3_atof(char *nstr);
static int flam3_atoi(char *nstr);
static void copy_matrix(double to[3][2], double from[3][2]);
static void convert_linear_to_polar(flam3_genome *cp, int ncps, int xfi, int cflag, double cxang[4][2], double cxmag[4][2], double cxtrn[4][2]);
static void interp_and_convert_back(double *c, int ncps, int xfi, double cxang[4][2], double cxmag[4][2], double cxtrn[4][2],double store_array[3][2]);
void prepare_xform_fn_ptrs(flam3_genome *, randctx *);
static void initialize_xforms(flam3_genome *thiscp, int start_here);
static int parse_flame_element(xmlNode *);
static int parse_xform_xml(xmlNode *chld_node,flam3_xform *this_xform, int *num_xaos, flam3_chaos_entry *xaos, int numstd, int motionxf);
static int apply_xform(flam3_genome *cp, int fn, double *p, double *q, randctx *rc);
static double adjust_percentage(double in);
static void xform_precalc(flam3_genome *cp, int xi);

void decode64( char *instring, char *outstring );
void encode64( FILE *infile, FILE *outfile, int linesize);
void decodeblock64 (unsigned char in[4], unsigned char out[3]);
void encodeblock64 (unsigned char in[3], unsigned char out[4], int len);
void b64decode(char* instr, char *outstr);

/*
 * VARIATION FUNCTIONS
 * must be of the form void (void *, double)
 */
static void var0_linear (void *helper, double weight) {
   /* linear */
   /* nx = tx;
      ny = ty;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   f->p0 += weight * f->tx;
   f->p1 += weight * f->ty;
}

static void var1_sinusoidal (void *helper, double weight) {
   /* sinusoidal */
   /* nx = sin(tx);
      ny = sin(ty);
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   f->p0 += weight * sin(f->tx);
   f->p1 += weight * sin(f->ty);
}

static void var2_spherical (void *helper, double weight) {
   /* spherical */
   /* double r2 = tx * tx + ty * ty + 1e-6;
      nx = tx / r2;
      ny = ty / r2;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r2 = weight / ( f->precalc_sumsq + EPS);

   f->p0 += r2 * f->tx;
   f->p1 += r2 * f->ty;
}

static void var3_swirl (void *helper, double weight) {
   /* swirl */
   /* double r2 = tx * tx + ty * ty;    /k here is fun
      double c1 = sin(r2);
      double c2 = cos(r2);
      nx = c1 * tx - c2 * ty;
      ny = c2 * tx + c1 * ty;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r2 = f->precalc_sumsq;
   double c1,c2;
   double nx,ny;
   
   sincos(r2,&c1,&c2);
//   double c1 = sin(r2);
//   double c2 = cos(r2);
   nx = c1 * f->tx - c2 * f->ty;
   ny = c2 * f->tx + c1 * f->ty;

   f->p0 += weight * nx;
   f->p1 += weight * ny;
}

static void var4_horseshoe (void *helper, double weight) {
   /* horseshoe */
   /* a = atan2(tx, ty);
      c1 = sin(a);
      c2 = cos(a);
      nx = c1 * tx - c2 * ty;
      ny = c2 * tx + c1 * ty;
      p[0] += v * nx;
      p[1] += v * ny;  */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   double r = weight / (f->precalc_sqrt + EPS);

   f->p0 += (f->tx - f->ty) * (f->tx + f->ty) * r;
   f->p1 += 2.0 * f->tx * f->ty * r;
}

static void var5_polar (void *helper, double weight) {
   /* polar */
   /* nx = atan2(tx, ty) / M_PI;
      ny = sqrt(tx * tx + ty * ty) - 1.0;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double nx = f->precalc_atan * M_1_PI;
   double ny = f->precalc_sqrt - 1.0;

   f->p0 += weight * nx;
   f->p1 += weight * ny;
}

static void var6_handkerchief (void *helper, double weight) {
   /* folded handkerchief */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      p[0] += v * sin(a+r) * r;
      p[1] += v * cos(a-r) * r; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double a = f->precalc_atan;
   double r = f->precalc_sqrt;

   f->p0 += weight * r * sin(a+r);
   f->p1 += weight * r * cos(a-r);
}

static void var7_heart (void *helper, double weight) {
   /* heart */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      a *= r;
      p[0] += v * sin(a) * r;
      p[1] += v * cos(a) * -r; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double a = f->precalc_sqrt * f->precalc_atan;
   double ca,sa;
   double r = weight * f->precalc_sqrt;
   
   sincos(a,&sa,&ca);

   f->p0 += r * sa;
   f->p1 += (-r) * ca;
}

static void var8_disc (void *helper, double weight) {
   /* disc */
   /* nx = tx * M_PI;
      ny = ty * M_PI;
      a = atan2(nx, ny);
      r = sqrt(nx*nx + ny*ny);
      p[0] += v * sin(r) * a / M_PI;
      p[1] += v * cos(r) * a / M_PI; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double a = f->precalc_atan * M_1_PI;
   double r = M_PI * f->precalc_sqrt;
   double sr,cr;
   sincos(r,&sr,&cr);

   f->p0 += weight * sr * a;
   f->p1 += weight * cr * a;
}

static void var9_spiral (void *helper, double weight) {
   /* spiral */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty) + 1e-6;
      p[0] += v * (cos(a) + sin(r)) / r;
      p[1] += v * (sin(a) - cos(r)) / r; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = f->precalc_sqrt + EPS;
   double r1 = weight/r;
   double sr,cr;
   sincos(r,&sr,&cr);

   f->p0 += r1 * (f->precalc_cosa + sr);
   f->p1 += r1 * (f->precalc_sina - cr);
}

static void var10_hyperbolic (void *helper, double weight) {
   /* hyperbolic */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty) + 1e-6;
      p[0] += v * sin(a) / r;
      p[1] += v * cos(a) * r; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = f->precalc_sqrt + EPS;

   f->p0 += weight * f->precalc_sina / r;
   f->p1 += weight * f->precalc_cosa * r;
}

static void var11_diamond (void *helper, double weight) {
   /* diamond */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      p[0] += v * sin(a) * cos(r);
      p[1] += v * cos(a) * sin(r); */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = f->precalc_sqrt;
   double sr,cr;
   sincos(r,&sr,&cr);

   f->p0 += weight * f->precalc_sina * cr;
   f->p1 += weight * f->precalc_cosa * sr;
}

static void var12_ex (void *helper, double weight) {
   /* ex */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      n0 = sin(a+r);
      n1 = cos(a-r);
      m0 = n0 * n0 * n0 * r;
      m1 = n1 * n1 * n1 * r;
      p[0] += v * (m0 + m1);
      p[1] += v * (m0 - m1); */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double a = f->precalc_atan;
   double r = f->precalc_sqrt;

   double n0 = sin(a+r);
   double n1 = cos(a-r);

   double m0 = n0 * n0 * n0 * r;
   double m1 = n1 * n1 * n1 * r;

   f->p0 += weight * (m0 + m1);
   f->p1 += weight * (m0 - m1);
}

static void var13_julia (void *helper, double weight) {
   /* julia */
   /* a = atan2(tx, ty)/2.0;
      if (flam3_random_bit()) a += M_PI;
      r = pow(tx*tx + ty*ty, 0.25);
      nx = r * cos(a);
      ny = r * sin(a);
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r;
   double a = 0.5 * f->precalc_atan;
   double sa,ca;

   if (flam3_random_isaac_bit(f->rc)) //(flam3_random_bit())
      a += M_PI;

   r = weight * sqrt(f->precalc_sqrt);
   
   sincos(a,&sa,&ca);

   f->p0 += r * ca;
   f->p1 += r * sa;
}

static void var14_bent (void *helper, double weight) {
   /* bent */
   /* nx = tx;
      ny = ty;
      if (nx < 0.0) nx = nx * 2.0;
      if (ny < 0.0) ny = ny / 2.0;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double nx = f->tx;
   double ny = f->ty;

   if (nx < 0.0)
      nx = nx * 2.0;
   if (ny < 0.0)
      ny = ny / 2.0;

   f->p0 += weight * nx;
   f->p1 += weight * ny;
}

static void var15_waves (void *helper, double weight) {
   /* waves */
   /* dx = coef[2][0];
      dy = coef[2][1];
      nx = tx + coef[1][0]*sin(ty/((dx*dx)+EPS));
      ny = ty + coef[1][1]*sin(tx/((dy*dy)+EPS));
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double c10 = f->xform->c[1][0];
   double c11 = f->xform->c[1][1];

   double nx = f->tx + c10 * sin( f->ty * f->xform->waves_dx2 );
   double ny = f->ty + c11 * sin( f->tx * f->xform->waves_dy2 );

   f->p0 += weight * nx;
   f->p1 += weight * ny;
}

static void var16_fisheye (void *helper, double weight) {
   /* fisheye */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      r = 2 * r / (r + 1);
      nx = r * cos(a);
      ny = r * sin(a);
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = f->precalc_sqrt;

   r = 2 * weight / (r+1);

   f->p0 += r * f->ty;
   f->p1 += r * f->tx;
}

static void var17_popcorn (void *helper, double weight) {
   /* popcorn */
   /* dx = tan(3*ty);
      dy = tan(3*tx);
      nx = tx + coef[2][0] * sin(dx);
      ny = ty + coef[2][1] * sin(dy);
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double dx = tan(3*f->ty);
   double dy = tan(3*f->tx);

   double nx = f->tx + f->xform->c[2][0] * sin(dx);
   double ny = f->ty + f->xform->c[2][1] * sin(dy);

   f->p0 += weight * nx;
   f->p1 += weight * ny;
}

static void var18_exponential (void *helper, double weight) {
   /* exponential */
   /* dx = exp(tx-1.0);
      dy = M_PI * ty;
      nx = cos(dy) * dx;
      ny = sin(dy) * dx;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double dx = weight * exp(f->tx - 1.0);
   double dy = M_PI * f->ty;
   double sdy,cdy;
   
   sincos(dy,&sdy,&cdy);
   

   f->p0 += dx * cdy;
   f->p1 += dx * sdy;
}

static void var19_power (void *helper, double weight) {
   /* power */
   /* a = atan2(tx, ty);
      sa = sin(a);
      r = sqrt(tx*tx + ty*ty);
      r = pow(r, sa);
      nx = r * precalc_cosa;
      ny = r * sa;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = weight * pow(f->precalc_sqrt, f->precalc_sina);

   f->p0 += r * f->precalc_cosa;
   f->p1 += r * f->precalc_sina;
}

static void var20_cosine (void *helper, double weight) {
   /* cosine */
   /* nx = cos(tx * M_PI) * cosh(ty);
      ny = -sin(tx * M_PI) * sinh(ty);
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double a = f->tx * M_PI;
   double sa,ca;
   double nx,ny;
   
   sincos(a,&sa,&ca);
   nx =  ca * cosh(f->ty);
   ny = -sa * sinh(f->ty);

   f->p0 += weight * nx;
   f->p1 += weight * ny;
}

static void var21_rings (void *helper, double weight) {
   /* rings */
   /* dx = coef[2][0];
      dx = dx * dx + EPS;
      r = sqrt(tx*tx + ty*ty);
      r = fmod(r + dx, 2*dx) - dx + r*(1-dx);
      a = atan2(tx, ty);
      nx = cos(a) * r;
      ny = sin(a) * r;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double dx = f->xform->c[2][0] * f->xform->c[2][0] + EPS;
   double r = f->precalc_sqrt;
   r = weight * (fmod(r+dx, 2*dx) - dx + r * (1 - dx));

   f->p0 += r * f->precalc_cosa;
   f->p1 += r * f->precalc_sina;
}

static void var22_fan (void *helper, double weight) {
   /* fan */
   /* dx = coef[2][0];
      dy = coef[2][1];
      dx = M_PI * (dx * dx + EPS);
      dx2 = dx/2;
      a = atan(tx,ty);
      r = sqrt(tx*tx + ty*ty);
      a += (fmod(a+dy, dx) > dx2) ? -dx2 : dx2;
      nx = cos(a) * r;
      ny = sin(a) * r;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double dx = M_PI * (f->xform->c[2][0] * f->xform->c[2][0] + EPS);
   double dy = f->xform->c[2][1];
   double dx2 = 0.5 * dx;

   double a = f->precalc_atan;
   double r = weight * f->precalc_sqrt;
   double sa,ca;

   a += (fmod(a+dy,dx) > dx2) ? -dx2 : dx2;
   sincos(a,&sa,&ca);

   f->p0 += r * ca;
   f->p1 += r * sa;
}

static void var23_blob (void *helper, double weight) {
   /* blob */
   /* a = atan2(tx, ty);
      r = sqrt(tx*tx + ty*ty);
      r = r * (bloblow + (blobhigh-bloblow) * (0.5 + 0.5 * sin(blobwaves * a)));
      nx = sin(a) * r;
      ny = cos(a) * r;

      p[0] += v * nx;
      p[1] += v * ny; */


   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = f->precalc_sqrt;
   double a = f->precalc_atan;
   double bdiff = f->xform->blob_high - f->xform->blob_low;

   r = r * (f->xform->blob_low +
            bdiff * (0.5 + 0.5 * sin(f->xform->blob_waves * a)));

   f->p0 += weight * f->precalc_sina * r;
   f->p1 += weight * f->precalc_cosa * r;
}

static void var24_pdj (void *helper, double weight) {
   /* pdj */
   /* nx1 = cos(pdjb * tx);
      nx2 = sin(pdjc * tx);
      ny1 = sin(pdja * ty);
      ny2 = cos(pdjd * ty);

      p[0] += v * (ny1 - nx1);
      p[1] += v * (nx2 - ny2); */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double nx1 = cos(f->xform->pdj_b * f->tx);
   double nx2 = sin(f->xform->pdj_c * f->tx);
   double ny1 = sin(f->xform->pdj_a * f->ty);
   double ny2 = cos(f->xform->pdj_d * f->ty);

   f->p0 += weight * (ny1 - nx1);
   f->p1 += weight * (nx2 - ny2);
}

static void var25_fan2 (void *helper, double weight) {
   /* fan2 */
   /* a = precalc_atan;
      r = precalc_sqrt;

      dy = fan2y;
      dx = M_PI * (fan2x * fan2x + EPS);
      dx2 = dx / 2.0;

      t = a + dy - dx * (int)((a + dy)/dx);

      if (t > dx2)
         a = a - dx2;
      else
         a = a + dx2;

      nx = sin(a) * r;
      ny = cos(a) * r;

      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   double dy = f->xform->fan2_y;
   double dx = M_PI * (f->xform->fan2_x * f->xform->fan2_x + EPS);
   double dx2 = 0.5 * dx;
   double a = f->precalc_atan;
   double sa,ca;
   double r = weight * f->precalc_sqrt;

   double t = a + dy - dx * (int)((a + dy)/dx);

   if (t>dx2)
      a = a-dx2;
   else
      a = a+dx2;
      
   sincos(a,&sa,&ca);

   f->p0 += r * sa;
   f->p1 += r * ca;
}

static void var26_rings2 (void *helper, double weight) {
   /* rings2 */
   /* r = precalc_sqrt;
      dx = rings2val * rings2val + EPS;
      r += dx - 2.0*dx*(int)((r + dx)/(2.0 * dx)) - dx + r * (1.0-dx);
      nx = precalc_sina * r;
      ny = precalc_cosa * r;
      p[0] += v * nx;
      p[1] += v * ny; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = f->precalc_sqrt;
   double dx = f->xform->rings2_val * f->xform->rings2_val + EPS;

   r += -2.0*dx*(int)((r+dx)/(2.0*dx)) + r * (1.0-dx);

   f->p0 += weight * f->precalc_sina * r;
   f->p1 += weight * f->precalc_cosa * r;
}

static void var27_eyefish (void *helper, double weight) {
   /* eyefish */
   /* r = 2.0 * v / (precalc_sqrt + 1.0);
      p[0] += r*tx;
      p[1] += r*ty; */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = (weight * 2.0) / (f->precalc_sqrt + 1.0);

   f->p0 += r * f->tx;
   f->p1 += r * f->ty;
}

static void var28_bubble (void *helper, double weight) {
   /* bubble */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double r = weight / (0.25 * (f->precalc_sumsq) + 1);

  f->p0 += r * f->tx;
  f->p1 += r * f->ty;
}

static void var29_cylinder (void *helper, double weight) {
   /* cylinder (01/06) */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   f->p0 += weight * sin(f->tx);
   f->p1 += weight * f->ty;
}

static void var30_perspective (void *helper, double weight) {
   /* perspective (01/06) */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double t = 1.0 / (f->xform->perspective_dist - f->ty * f->xform->persp_vsin + EPS);

   f->p0 += weight * f->xform->perspective_dist * f->tx * t;
   f->p1 += weight * f->xform->persp_vfcos * f->ty * t;
}

static void var31_noise (void *helper, double weight) {
   /* noise (03/06) */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double tmpr, sinr, cosr, r;

   tmpr = flam3_random_isaac_01(f->rc) * 2 * M_PI;
   sincos(tmpr,&sinr,&cosr);

   r = weight * flam3_random_isaac_01(f->rc);

   f->p0 += f->tx * r * cosr;
   f->p1 += f->ty * r * sinr;
}

static void var32_juliaN_generic (void *helper, double weight) {
   /* juliaN (03/06) */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   int t_rnd = trunc((f->xform->juliaN_rN)*flam3_random_isaac_01(f->rc));
   
   double tmpr = (f->precalc_atanyx + 2 * M_PI * t_rnd) / f->xform->juliaN_power;

   double r = weight * pow(f->precalc_sumsq, f->xform->juliaN_cn);
   double sina, cosa;
   sincos(tmpr,&sina,&cosa);

   f->p0 += r * cosa;
   f->p1 += r * sina;
}

static void var33_juliaScope_generic (void *helper, double weight) {
   /* juliaScope (03/06) */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   int t_rnd = trunc((f->xform->juliaScope_rN) * flam3_random_isaac_01(f->rc));

   double tmpr, r;
   double sina, cosa;

   if ((t_rnd & 1) == 0)
      tmpr = (2 * M_PI * t_rnd + f->precalc_atanyx) / f->xform->juliaScope_power;
   else
      tmpr = (2 * M_PI * t_rnd - f->precalc_atanyx) / f->xform->juliaScope_power;

   sincos(tmpr,&sina,&cosa);

   r = weight * pow(f->precalc_sumsq, f->xform->juliaScope_cn);

   f->p0 += r * cosa;
   f->p1 += r * sina;
}

static void var34_blur (void *helper, double weight) {
   /* blur (03/06) */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double tmpr, sinr, cosr, r;

   tmpr = flam3_random_isaac_01(f->rc) * 2 * M_PI;
   sincos(tmpr,&sinr,&cosr);

   r = weight * flam3_random_isaac_01(f->rc);

   f->p0 += r * cosr;
   f->p1 += r * sinr;
}

static void var35_gaussian (void *helper, double weight) {
   /* gaussian (09/06) */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double ang, r, sina, cosa;

   ang = flam3_random_isaac_01(f->rc) * 2 * M_PI;
   sincos(ang,&sina,&cosa);

   r = weight * ( flam3_random_isaac_01(f->rc) + flam3_random_isaac_01(f->rc)
                   + flam3_random_isaac_01(f->rc) + flam3_random_isaac_01(f->rc) - 2.0 );

   f->p0 += r * cosa;
   f->p1 += r * sina;
}

static void var36_radial_blur (void *helper, double weight) {
   /* radial blur (09/06) */
   /* removed random storage 6/07 */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double rndG, ra, rz, tmpa, sa, ca;

   /* Get pseudo-gaussian */
   rndG = weight * (flam3_random_isaac_01(f->rc) + flam3_random_isaac_01(f->rc)
                   + flam3_random_isaac_01(f->rc) + flam3_random_isaac_01(f->rc) - 2.0);

   /* Calculate angle & zoom */
   ra = f->precalc_sqrt;
   tmpa = f->precalc_atanyx + f->xform->radialBlur_spinvar*rndG;
   sincos(tmpa,&sa,&ca);
   rz = f->xform->radialBlur_zoomvar * rndG - 1;

   f->p0 += ra * ca + rz * f->tx;
   f->p1 += ra * sa + rz * f->ty;
}

static void var37_pie(void *helper, double weight) {
   /* pie by Joel Faber (June 2006) */
   flam3_iter_helper *f = (flam3_iter_helper *) helper;
   double a, r, sa, ca;
   int sl;

   sl = (int) (flam3_random_isaac_01(f->rc) * f->xform->pie_slices + 0.5);
   a = f->xform->pie_rotation +
       2.0 * M_PI * (sl + flam3_random_isaac_01(f->rc) * f->xform->pie_thickness) / f->xform->pie_slices;
   r = weight * flam3_random_isaac_01(f->rc);
   sincos(a,&sa,&ca);

   f->p0 += r * ca;
   f->p1 += r * sa;
}

static void var38_ngon(void *helper, double weight) {
   /* ngon by Joel Faber (09/06) */
   flam3_iter_helper *f = (flam3_iter_helper *) helper;
   double r_factor,theta,phi,b, amp;

   r_factor = pow(f->precalc_sumsq, f->xform->ngon_power/2.0);

   theta = f->precalc_atanyx;
   b = 2*M_PI/f->xform->ngon_sides;

   phi = theta - (b*floor(theta/b));
   if (phi > b/2)
      phi -= b;

   amp = f->xform->ngon_corners * (1.0 / (cos(phi) + EPS) - 1.0) + f->xform->ngon_circle;
   amp /= (r_factor + EPS);

   f->p0 += weight * f->tx * amp;
   f->p1 += weight * f->ty * amp;
}

static void var39_curl(void *helper, double weight)
{
    flam3_iter_helper *f = (flam3_iter_helper *) helper;

    double re = 1.0 + f->xform->curl_c1 * f->tx + f->xform->curl_c2 * (f->tx * f->tx - f->ty * f->ty);
    double im = f->xform->curl_c1 * f->ty + 2.0 * f->xform->curl_c2 * f->tx * f->ty;

    double r = weight / (re*re + im*im);

    f->p0 += (f->tx * re + f->ty * im) * r;
    f->p1 += (f->ty * re - f->tx * im) * r;
}

static void var40_rectangles(void *helper, double weight)
{
    flam3_iter_helper *f = (flam3_iter_helper *) helper;

    if (f->xform->rectangles_x==0)
       f->p0 += weight * f->tx;
    else
       f->p0 += weight * ((2 * floor(f->tx / f->xform->rectangles_x) + 1) * f->xform->rectangles_x - f->tx);

    if (f->xform->rectangles_y==0)
       f->p1 += weight * f->ty;
    else
       f->p1 += weight * ((2 * floor(f->ty / f->xform->rectangles_y) + 1) * f->xform->rectangles_y - f->ty);

}

static void var41_arch(void *helper, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Arch;
   var
     sinr, cosr: double;
   begin
     SinCos(random * vars[29]*pi, sinr, cosr);
     FPx := FPx + sinr*vars[29];
     FPy := FPy + sqr(sinr)/cosr*vars[29];
   end;
   */
   
   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   double ang = flam3_random_isaac_01(f->rc) * weight * M_PI;
   double sinr,cosr;
   sincos(ang,&sinr,&cosr);

   f->p0 += weight * sinr;
   f->p1 += weight * (sinr*sinr)/cosr;

}

static void var42_tangent(void *helper, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Tangent;
   begin
     FPx := FPx + vars[30] * (sin(FTx)/cos(FTy));
     FPy := FPy + vars[30] * (sin(FTy)/cos(FTy));
   end;
   */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   f->p0 += weight * sin(f->tx)/cos(f->ty);
   f->p1 += weight * tan(f->ty);

}

static void var43_square(void *helper, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.SquareBlur;
   begin
     FPx := FPx + vars[31] * (random - 0.5);
     FPy := FPy + vars[31] * (random - 0.5);
   end;
   */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   f->p0 += weight * (flam3_random_isaac_01(f->rc) - 0.5);
   f->p1 += weight * (flam3_random_isaac_01(f->rc) - 0.5);

}

static void var44_rays(void *helper, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Rays;
   var
     r, sinr, cosr, tgr: double;
   begin
     SinCos(random * vars[32]*pi, sinr, cosr);
     r := vars[32] / (sqr(FTx) + sqr(FTy) + EPS);
     tgr := sinr/cosr;
     FPx := FPx + tgr * (cos(FTx)*vars[32]) * r;
     FPy := FPy + tgr * (sin(FTy)*vars[32]) * r;
   end;
   */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   double ang = weight * flam3_random_isaac_01(f->rc) * M_PI;
   double r = weight / (f->precalc_sumsq + EPS);
   double tanr = weight * tan(ang) * r;


   f->p0 += tanr * cos(f->tx);
   f->p1 += tanr * sin(f->ty);

}

static void var45_blade(void *helper, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Blade;
   var
     r, sinr, cosr: double;
   begin
     r := sqrt(sqr(FTx) + sqr(FTy))*vars[33];
     SinCos(r*random, sinr, cosr);
     FPx := FPx + vars[33] * FTx * (cosr + sinr);
     FPy := FPy + vars[33] * FTx * (cosr - sinr);
   end;
   */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   double r = flam3_random_isaac_01(f->rc) * weight * f->precalc_sqrt;
   double sinr,cosr;
   
   sincos(r,&sinr,&cosr);

   f->p0 += weight * f->tx * (cosr + sinr);
   f->p1 += weight * f->tx * (cosr - sinr);

}

static void var46_secant2(void *helper, double weight)
{
   /* Intended as a 'fixed' version of secant */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   double r = weight * f->precalc_sqrt;
   double cr = cos(r);
   double icr = 1.0/cr;

   f->p0 += weight * f->tx;
   
   if (cr<0)
      f->p1 += weight*(icr + 1);
   else
      f->p1 += weight*(icr - 1);
}

static void var47_twintrian(void *helper, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.TwinTrian;
   var
     r, diff, sinr, cosr: double;
   begin
     r := sqrt(sqr(FTx) + sqr(FTy))*vars[35];
     SinCos(r*random, sinr, cosr);
     diff := Math.Log10(sinr*sinr)+cosr;
     FPx := FPx + vars[35] * FTx * diff;
     FPy := FPy + vars[35] * FTx * (diff - (sinr*pi));
   end;
   */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   double r = flam3_random_isaac_01(f->rc) * weight * f->precalc_sqrt;
   double sinr,cosr,diff;
   
   sincos(r,&sinr,&cosr);
   diff = log10(sinr*sinr)+cosr;
   
   if (badvalue(diff))
      diff = -30.0;      

   f->p0 += weight * f->tx * diff;
   f->p1 += weight * f->tx * (diff - sinr*M_PI);

}

static void var48_cross(void *helper, double weight)
{
   /* Z+ variation Jan 07
   procedure TXForm.Cross;
   var
     r: double;
   begin
     r := vars[36]*sqrt(1/(sqr(sqr(FTx)-sqr(FTy))+EPS));
     FPx := FPx + FTx * r;
     FPy := FPy + FTy * r;
   end;
   */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   double s = f->tx*f->tx - f->ty*f->ty;
   double r = weight * sqrt(1.0 / (s*s+EPS));

   f->p0 += f->tx * r;
   f->p1 += f->ty * r;

}

static void var49_disc2(void *helper, double weight)
{
   /* Z+ variation Jan 07
   c := vvar/PI;
   k := rot*PI;
     sinadd := Sin(add);
     cosadd := Cos(add);
   cosadd := cosadd - 1;
   if (add > 2*PI) then begin
     cosadd := cosadd * (1 + add - 2*PI);
     sinadd := sinadd * (1 + add - 2*PI)
   end
   else if (add < -2*PI) then begin
     cosadd := cosadd * (1 + add + 2*PI);
     sinadd := sinadd * (1 + add + 2*PI)
   end
   end;
   procedure TVariationDisc2.CalcFunction;
   var
     r, sinr, cosr: extended;
   begin
     SinCos(k * (FTx^+FTy^), sinr, cosr);   //rot*PI
     r := c * arctan2(FTx^, FTy^); //vvar/PI
     FPx^ := FPx^ + (sinr + cosadd) * r;
     FPy^ := FPy^ + (cosr + sinadd) * r;
   */

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   double r,t,sinr, cosr;

   t = f->xform->disc2_timespi * (f->tx + f->ty);
   sincos(t,&sinr,&cosr);
   r = weight * f->precalc_atan / M_PI;

   f->p0 += (sinr + f->xform->disc2_cosadd) * r;
   f->p1 += (cosr + f->xform->disc2_sinadd) * r;

}

static void var50_supershape(void *helper, double weight) {

   flam3_iter_helper *f = (flam3_iter_helper *) helper;

   double theta;
   double t1,t2,r;
   double st,ct;
   double myrnd;

   theta = f->xform->supershape_pm_4 * f->precalc_atanyx + M_PI_4;
   
   sincos(theta,&st,&ct);

   t1 = fabs(ct);
   t1 = pow(t1,f->xform->supershape_n2);

   t2 = fabs(st);
   t2 = pow(t2,f->xform->supershape_n3);
   
   myrnd = f->xform->supershape_rnd;

   r = weight * ( (myrnd*flam3_random_isaac_01(f->rc) + (1.0-myrnd)*f->precalc_sqrt) - f->xform->supershape_holes) 
      * pow(t1+t2,f->xform->supershape_pneg1_n1) / f->precalc_sqrt;

   f->p0 += r * f->tx;
   f->p1 += r * f->ty;
}

static void var51_flower(void *helper, double weight) {
    /* cyberxaos, 4/2007 */
    /*   theta := arctan2(FTy^, FTx^);
         r := (random-holes)*cos(petals*theta);
         FPx^ := FPx^ + vvar*r*cos(theta);
         FPy^ := FPy^ + vvar*r*sin(theta);*/
         
    flam3_iter_helper *f = (flam3_iter_helper *)helper;
    double theta = f->precalc_atanyx;
    double r = weight * (flam3_random_isaac_01(f->rc) - f->xform->flower_holes) * 
                    cos(f->xform->flower_petals*theta) / f->precalc_sqrt;

    f->p0 += r * f->tx;
    f->p1 += r * f->ty;
}
    
static void var52_conic(void *helper, double weight) {
    /* cyberxaos, 4/2007 */
    /*   theta := arctan2(FTy^, FTx^);
         r :=  (random - holes)*((eccentricity)/(1+eccentricity*cos(theta)));
         FPx^ := FPx^ + vvar*r*cos(theta);
         FPy^ := FPy^ + vvar*r*sin(theta); */
         
    flam3_iter_helper *f = (flam3_iter_helper *)helper;
    double ct = f->tx / f->precalc_sqrt;
    double r = weight * (flam3_random_isaac_01(f->rc) - f->xform->conic_holes) * 
                    f->xform->conic_eccen / (1 + f->xform->conic_eccen*ct) / f->precalc_sqrt;

    f->p0 += r * f->tx;
    f->p1 += r * f->ty;
}

static void var53_parabola(void *helper, double weight) {
    /* cyberxaos, 4/2007 */
    /*   r := sqrt(sqr(FTx^) + sqr(FTy^));
         FPx^ := FPx^ + parabola_height*vvar*sin(r)*sin(r)*random;  
         FPy^ := FPy^ + parabola_width*vvar*cos(r)*random; */
         
    flam3_iter_helper *f = (flam3_iter_helper *)helper;
    double r = f->precalc_sqrt;
    double sr,cr;
    
    sincos(r,&sr,&cr);
    
    f->p0 += f->xform->parabola_height * weight * sr*sr * flam3_random_isaac_01(f->rc);
    f->p1 += f->xform->parabola_width * weight * cr * flam3_random_isaac_01(f->rc);
    
}      

static void var54_bent2 (void *helper, double weight) {

   /* Bent2 in the Apophysis Plugin Pack */   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;   
   double nx = f->tx;
   double ny = f->ty;

   if (nx < 0.0)
      nx = nx * f->xform->bent2_x;
   if (ny < 0.0)
      ny = ny * f->xform->bent2_y;

   f->p0 += weight * nx;
   f->p1 += weight * ny;
}

static void var55_bipolar (void *helper, double weight) {

   /* Bipolar in the Apophysis Plugin Pack */   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;   
   double x2y2 = f->precalc_sumsq;
   double t = x2y2+1;
   double x2 = 2*f->tx;
   double ps = -M_PI_2 * f->xform->bipolar_shift;
   double y = 0.5 * atan2(2.0 * f->ty, x2y2 - 1.0) + ps;
   
   if (y > M_PI_2)
       y = -M_PI_2 + fmod(y + M_PI_2, M_PI);
   else if (y < -M_PI_2)
       y = M_PI_2 - fmod(M_PI_2 - y, M_PI);

   f->p0 += weight * 0.25 * M_2_PI * log ( (t+x2) / (t-x2) );
   f->p1 += weight * M_2_PI * y;
}

static void var56_boarders (void *helper, double weight) {

   /* Boarders in the Apophysis Plugin Pack */   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;   
   double roundX, roundY, offsetX, offsetY;
    
   roundX = rint(f->tx);
   roundY = rint(f->ty);
   offsetX = f->tx - roundX;
   offsetY = f->ty - roundY;
    
   if (flam3_random_isaac_01(f->rc) >= 0.75) {
      f->p0 += weight*(offsetX*0.5 + roundX);
      f->p1 += weight*(offsetY*0.5 + roundY);
   } else {
      
      if (fabs(offsetX) >= fabs(offsetY)) {
         
         if (offsetX >= 0.0) {
            f->p0 += weight*(offsetX*0.5 + roundX + 0.25);
            f->p1 += weight*(offsetY*0.5 + roundY + 0.25 * offsetY / offsetX);
         } else {
            f->p0 += weight*(offsetX*0.5 + roundX - 0.25);
            f->p1 += weight*(offsetY*0.5 + roundY - 0.25 * offsetY / offsetX);  
         }
         
      } else {
         
         if (offsetY >= 0.0) {
            f->p1 += weight*(offsetY*0.5 + roundY + 0.25);
            f->p0 += weight*(offsetX*0.5 + roundX + offsetX/offsetY*0.25);
         } else {
            f->p1 += weight*(offsetY*0.5 + roundY - 0.25);
            f->p0 += weight*(offsetX*0.5 + roundX - offsetX/offsetY*0.25);
         }
      }
   }
}

static void var57_butterfly (void *helper, double weight) {

   /* Butterfly in the Apophysis Plugin Pack */   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;   
   
   /* wx is weight*4/sqrt(3*pi) */
   double wx = weight*1.3029400317411197908970256609023;
   
   double y2 = f->ty*2.0;
   double r = wx*sqrt(fabs(f->ty * f->tx)/(EPS + f->tx*f->tx + y2*y2));
   
   f->p0 += r * f->tx;
   f->p1 += r * y2;
   
}

static void var58_cell (void *helper, double weight) {

   /* Cell in the Apophysis Plugin Pack */   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;   

   double inv_cell_size = 1.0/f->xform->cell_size;
    
   /* calculate input cell */
   int x = floor(f->tx*inv_cell_size);
   int y = floor(f->ty*inv_cell_size);

   /* Offset from cell origin */
   double dx = f->tx - x*f->xform->cell_size;
   double dy = f->ty - y*f->xform->cell_size;
   
   /* interleave cells */
   if (y >= 0) {
      if (x >= 0) {
         y *= 2;
         x *= 2;
      } else {
         y *= 2;
         x = -(2*x+1);
      }
   } else {
      if (x >= 0) {
         y = -(2*y+1);
         x *= 2;
      } else {
         y = -(2*y+1);
         x = -(2*x+1);
      }
   }
   
   f->p0 += weight * (dx + x*f->xform->cell_size);
   f->p1 -= weight * (dy + y*f->xform->cell_size);
   
}

static void var59_cpow (void *helper, double weight) {

   /* Cpow in the Apophysis Plugin Pack */   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;   
   
   double a = f->precalc_atanyx;
   double lnr = 0.5 * log(f->precalc_sumsq);
   double va = 2.0 * M_PI / f->xform->cpow_power;
   double vc = f->xform->cpow_r / f->xform->cpow_power;
   double vd = f->xform->cpow_i / f->xform->cpow_power;
   double ang = vc*a + vd*lnr + va*floor(f->xform->cpow_power*flam3_random_isaac_01(f->rc));
   double sa,ca;
   
   double m = weight * exp(vc * lnr - vd * a);
   
   sincos(ang,&sa,&ca);
   
   f->p0 += m * ca;
   f->p1 += m * sa;
   
}

static void var60_curve (void *helper, double weight) {

   /* Curve in the Apophysis Plugin Pack */   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;   
   double pc_xlen = f->xform->curve_xlength*f->xform->curve_xlength;
   double pc_ylen = f->xform->curve_ylength*f->xform->curve_ylength;
   
   if (pc_xlen<1E-20) pc_xlen = 1E-20;
   
   if (pc_ylen<1E-20) pc_ylen = 1E-20;

   f->p0 += weight * (f->tx + f->xform->curve_xamp * exp(-f->ty*f->ty/pc_xlen));
   f->p1 += weight * (f->ty + f->xform->curve_yamp * exp(-f->tx*f->tx/pc_ylen));
      
}

static void var61_edisc (void *helper, double weight) {

   /* Edisc in the Apophysis Plugin Pack */   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double tmp = f->precalc_sumsq + 1.0;
   double tmp2 = 2.0 * f->tx;
   double r1 = sqrt(tmp+tmp2);
   double r2 = sqrt(tmp-tmp2);
   double xmax = (r1+r2) * 0.5;
   double a1 = log(xmax + sqrt(xmax - 1.0));
   double a2 = -acos(f->tx/xmax);
   double w = weight / 11.57034632;
   double snv,csv,snhu,cshu;
   
   sincos(a1,&snv,&csv);
   
   snhu = sinh(a2);
   cshu = cosh(a2);
   
   if (f->ty > 0.0) snv = -snv;
   
   f->p0 += w * cshu * csv;
   f->p1 += w * snhu * snv;
   
}

static void var62_elliptic (void *helper, double weight) {

   /* Elliptic in the Apophysis Plugin Pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   double tmp = f->precalc_sumsq + 1.0;
   double x2 = 2.0 * f->tx;
   double xmax = 0.5 * (sqrt(tmp+x2) + sqrt(tmp-x2));
   double a = f->tx / xmax;
   double b = 1.0 - a*a;
   double ssx = xmax - 1.0;
   double w = weight / M_PI_2;
   
   if (b<0)
      b = 0;
   else
      b = sqrt(b);
      
   if (ssx<0)
      ssx = 0;
   else
      ssx = sqrt(ssx);
      
   f->p0 += w * atan2(a,b);
   
   if (f->ty > 0)
      f->p1 += w * log(xmax + ssx);
   else
      f->p1 -= w * log(xmax + ssx);
      
}

static void var63_escher (void *helper, double weight) {

   /* Escher in the Apophysis Plugin Pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double seb,ceb;
   double vc,vd;
   double m,n;
   double sn,cn;

   double a = f->precalc_atanyx;
   double lnr = 0.5 * log(f->precalc_sumsq);

   sincos(f->xform->escher_beta,&seb,&ceb);
   
   vc = 0.5 * (1.0 + ceb);
   vd = 0.5 * seb;

   m = weight * exp(vc*lnr - vd*a);
   n = vc*a + vd*lnr;
   
   sincos(n,&sn,&cn);
   
   f->p0 += m * cn;
   f->p1 += m * sn;
      
}

static void var64_foci (void *helper, double weight) {

   /* Foci in the Apophysis Plugin Pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double expx = exp(f->tx) * 0.5;
   double expnx = 0.25 / expx;
   double sn,cn,tmp;
   
   sincos(f->ty,&sn,&cn);
   tmp = weight/(expx + expnx - cn);
   
   f->p0 += tmp * (expx - expnx);
   f->p1 += tmp * sn;
      
}

static void var65_lazysusan (void *helper, double weight) {

   /* Lazysusan in the Apophysis Plugin Pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double x = f->tx - f->xform->lazysusan_x;
   double y = f->ty + f->xform->lazysusan_y;
   double r = sqrt(x*x + y*y);
   double sina, cosa;
   
   if (r<weight) {
      double a = atan2(y,x) + f->xform->lazysusan_spin +
                 f->xform->lazysusan_twist*(weight-r);
      sincos(a,&sina,&cosa);
      r = weight * r;
      
      f->p0 += r*cosa + f->xform->lazysusan_x;
      f->p1 += r*sina - f->xform->lazysusan_y;
   } else {
      
      r = weight * (1.0 + f->xform->lazysusan_space / r);
      
      f->p0 += r*x + f->xform->lazysusan_x;
      f->p1 += r*y - f->xform->lazysusan_y;
   
   }
      
}

static void var66_loonie (void *helper, double weight) {

   /* Loonie in the Apophysis Plugin Pack */

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double r2 = f->precalc_sumsq;
   double w2 = weight*weight;
   
   if (r2 < w2) {
      double r = weight * sqrt(w2/r2 - 1.0);
      f->p0 += r * f->tx;
      f->p1 += r * f->ty;
   } else {
      f->p0 += weight * f->tx;
      f->p1 += weight * f->ty;
   }
         
}

static void var67_pre_blur (void *helper, double weight) {

   /* pre-xform: PreBlur (Apo 2.08) */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   /* Get pseudo-gaussian */
   double rndG = weight * (flam3_random_isaac_01(f->rc) + flam3_random_isaac_01(f->rc)
                   + flam3_random_isaac_01(f->rc) + flam3_random_isaac_01(f->rc) - 2.0);
   double rndA = flam3_random_isaac_01(f->rc) * 2.0 * M_PI;
   double sinA,cosA;
   
   sincos(rndA,&sinA,&cosA);
   
   /* Note: original coordinate changed */
   f->tx += rndG * cosA;
   f->ty += rndG * sinA;
         
}

static void var68_modulus (void *helper, double weight) {

   /* Modulus in the Apophysis Plugin Pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double xr = 2*f->xform->modulus_x;
   double yr = 2*f->xform->modulus_y;
   
   if (f->tx > f->xform->modulus_x)
      f->p0 += weight * (-f->xform->modulus_x + fmod(f->tx + f->xform->modulus_x, xr));
   else if (f->tx < -f->xform->modulus_x)
      f->p0 += weight * ( f->xform->modulus_x - fmod(f->xform->modulus_x - f->tx, xr));
   else
      f->p0 += weight * f->tx;
      
   if (f->ty > f->xform->modulus_y)
      f->p1 += weight * (-f->xform->modulus_y + fmod(f->ty + f->xform->modulus_y, yr));
   else if (f->ty < -f->xform->modulus_y)
      f->p1 += weight * ( f->xform->modulus_y - fmod(f->xform->modulus_y - f->ty, yr));
   else
      f->p1 += weight * f->ty;
         
}

static void var69_oscope (void *helper, double weight) {

   /* oscilloscope from the apophysis plugin pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double tpf = 2 * M_PI * f->xform->oscope_frequency;
   double t;
   
   if (f->xform->oscope_damping == 0.0)
      t = f->xform->oscope_amplitude * cos(tpf*f->tx) + f->xform->oscope_separation;
   else {
      t = f->xform->oscope_amplitude * exp(-fabs(f->tx)*f->xform->oscope_damping)
          * cos(tpf*f->tx) + f->xform->oscope_separation;
   }
   
   if (fabs(f->ty) <= t) {
      f->p0 += weight*f->tx;
      f->p1 -= weight*f->ty;
   } else {
      f->p0 += weight*f->tx;
      f->p1 += weight*f->ty;
   } 
}

static void var70_polar2 (void *helper, double weight) {

   /* polar2 from the apophysis plugin pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double p2v = weight / M_PI;
   
   f->p0 += p2v * f->precalc_atan;
   f->p1 += p2v/2.0 * log(f->precalc_sumsq);
}

static void var71_popcorn2 (void *helper, double weight) {

   /* popcorn2 from the apophysis plugin pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   f->p0 += weight * ( f->tx + f->xform->popcorn2_x * sin(tan(f->ty*f->xform->popcorn2_c)));
   f->p1 += weight * ( f->ty + f->xform->popcorn2_y * sin(tan(f->tx*f->xform->popcorn2_c)));

}

static void var72_scry (void *helper, double weight) {

   /* scry from the apophysis plugin pack */
   /* note that scry does not multiply by weight, but as the */
   /* values still approach 0 as the weight approaches 0, it */
   /* should be ok                                           */ 

   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   double t = f->precalc_sumsq;
   double r = 1.0 / (f->precalc_sqrt * (t + 1.0/(weight+EPS)));
   
   f->p0 += f->tx * r;
   f->p1 += f->ty * r;

}

static void var73_separation (void *helper, double weight) {

   /* separation from the apophysis plugin pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   double sx2 = f->xform->separation_x * f->xform->separation_x;
   double sy2 = f->xform->separation_y * f->xform->separation_y;
   
   if (f->tx > 0.0)
      f->p0 += weight * (sqrt(f->tx*f->tx + sx2)- f->tx*f->xform->separation_xinside);
   else
      f->p0 -= weight * (sqrt(f->tx*f->tx + sx2)+ f->tx*f->xform->separation_xinside);
   
   if (f->ty > 0.0)
      f->p1 += weight * (sqrt(f->ty*f->ty + sy2)- f->ty*f->xform->separation_yinside);
   else
      f->p1 -= weight * (sqrt(f->ty*f->ty + sy2)+ f->ty*f->xform->separation_yinside);
   
}

static void var74_split (void *helper, double weight) {
   
   /* Split from apo plugins pack */
   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   if (cos(f->tx*f->xform->split_xsize*M_PI) >= 0)
      f->p1 += weight*f->ty;
   else
      f->p1 -= weight*f->ty;
      
   if (cos(f->ty*f->xform->split_ysize*M_PI) >= 0)
      f->p0 += weight * f->tx;
   else
      f->p0 -= weight * f->tx;

}

static void var75_splits (void *helper, double weight) {
   
   /* Splits from apo plugins pack */
   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   if (f->tx >= 0)
      f->p0 += weight*(f->tx+f->xform->splits_x);
   else
      f->p0 += weight*(f->tx-f->xform->splits_x);
      
   if (f->ty >= 0)
      f->p1 += weight*(f->ty+f->xform->splits_y);
   else
      f->p1 += weight*(f->ty-f->xform->splits_y);

}

static void var76_stripes (void *helper, double weight) {
   
   /* Stripes from apo plugins pack */
   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   double roundx,offsetx;
   
   roundx = floor(f->tx + 0.5);
   offsetx = f->tx - roundx;

   f->p0 += weight * (offsetx*(1.0-f->xform->stripes_space)+roundx);
   f->p1 += weight * (f->ty + offsetx*offsetx*f->xform->stripes_warp);

}

static void var77_wedge (void *helper, double weight) {
   
   /* Wedge from apo plugins pack */
   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   double r = f->precalc_sqrt;
   double a = f->precalc_atanyx + f->xform->wedge_swirl * r;
   double c = floor( (f->xform->wedge_count * a + M_PI)*M_1_PI*0.5);
   
   double comp_fac = 1 - f->xform->wedge_angle*f->xform->wedge_count*M_1_PI*0.5;
   double sa, ca;
   
   a = a * comp_fac + c * f->xform->wedge_angle;
   
   sincos(a,&sa,&ca);

   r = weight * (r + f->xform->wedge_hole);
   
   f->p0 += r*ca;
   f->p1 += r*sa;

}

static void var78_wedge_julia (void *helper, double weight) {
   /* wedge_julia from apo plugin pack */
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   double r = weight * pow(f->precalc_sumsq, f->xform->wedgeJulia_cn);
   int t_rnd = (int)((f->xform->wedgeJulia_rN)*flam3_random_isaac_01(f->rc));
   double a = (f->precalc_atanyx + 2 * M_PI * t_rnd) / f->xform->wedge_julia_power;
   double c = floor( (f->xform->wedge_julia_count * a + M_PI)*M_1_PI*0.5 );
   double sa,ca;
   
   a = a * f->xform->wedgeJulia_cf + c * f->xform->wedge_julia_angle;
   
   sincos(a,&sa,&ca);

   f->p0 += r * ca;
   f->p1 += r * sa;
}

static void var79_wedge_sph (void *helper, double weight) {
   
   /* Wedge_sph from apo plugins pack */
   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   double r = 1.0/(f->precalc_sqrt+EPS);
   double a = f->precalc_atanyx + f->xform->wedge_sph_swirl * r;
   double c = floor( (f->xform->wedge_sph_count * a + M_PI)*M_1_PI*0.5);
   
   double comp_fac = 1 - f->xform->wedge_sph_angle*f->xform->wedge_sph_count*M_1_PI*0.5;
   double sa, ca;
   
   a = a * comp_fac + c * f->xform->wedge_sph_angle;

   sincos(a,&sa,&ca);   
   r = weight * (r + f->xform->wedge_sph_hole);
   
   f->p0 += r*ca;
   f->p1 += r*sa;

}

static void var80_whorl (void *helper, double weight) {
   
   /* whorl from apo plugins pack */
   
   /*
    * !!! Note !!!
    * This code uses the variation weight in a non-standard fashion, and
    * it may change or even be removed in future versions of flam3.
    */

   flam3_iter_helper *f = (flam3_iter_helper *)helper;

   double r = f->precalc_sqrt;
   double a,sa,ca;

   if (r<weight)
      a = f->precalc_atanyx + f->xform->whorl_inside/(weight-r);
   else
      a = f->precalc_atanyx + f->xform->whorl_outside/(weight-r);
   
   sincos(a,&sa,&ca);
   
   f->p0 += weight*r*ca;
   f->p1 += weight*r*sa;

}

static void var81_waves2 (void *helper, double weight) {
   
   /* waves2 from Joel F */
   
   flam3_iter_helper *f = (flam3_iter_helper *)helper;
   
   f->p0 += weight*(f->tx + f->xform->waves2_scalex)*sin(f->ty * f->xform->waves2_freqx);
   f->p1 += weight*(f->ty + f->xform->waves2_scaley)*sin(f->tx * f->xform->waves2_freqy);

}

/* Precalc functions */

static void perspective_precalc(flam3_xform *xf) {
   double ang = xf->perspective_angle * M_PI / 2.0;
   xf->persp_vsin = sin(ang);
   xf->persp_vfcos = xf->perspective_dist * cos(ang);
}

static void juliaN_precalc(flam3_xform *xf) {
   xf->juliaN_rN = fabs(xf->juliaN_power);
   xf->juliaN_cn = xf->juliaN_dist / (double)xf->juliaN_power / 2.0;
}

static void wedgeJulia_precalc(flam3_xform *xf) {
   xf->wedgeJulia_cf = 1.0 - xf->wedge_julia_angle * xf->wedge_julia_count * M_1_PI * 0.5;
   xf->wedgeJulia_rN = fabs(xf->wedge_julia_power);
   xf->wedgeJulia_cn = xf->wedge_julia_dist / xf->wedge_julia_power / 2.0;
}

static void juliaScope_precalc(flam3_xform *xf) {
   xf->juliaScope_rN = fabs(xf->juliaScope_power);
   xf->juliaScope_cn = xf->juliaScope_dist / (double)xf->juliaScope_power / 2.0;
}

static void radial_blur_precalc(flam3_xform *xf) {
   sincos(xf->radialBlur_angle * M_PI / 2.0,
             &xf->radialBlur_spinvar, &xf->radialBlur_zoomvar);
//   xf->radialBlur_spinvar = sin(xf->radialBlur_angle * M_PI / 2);
//   xf->radialBlur_zoomvar = cos(xf->radialBlur_angle * M_PI / 2);
}

static void waves_precalc(flam3_xform *xf) {
   double dx = xf->c[2][0];
   double dy = xf->c[2][1];

   xf->waves_dx2 = 1.0/(dx * dx + EPS);
   xf->waves_dy2 = 1.0/(dy * dy + EPS);
}

static void disc2_precalc(flam3_xform *xf) {
   double add = xf->disc2_twist;
   double k;

   xf->disc2_timespi = xf->disc2_rot * M_PI;

   sincos(add,&xf->disc2_sinadd,&xf->disc2_cosadd);
   xf->disc2_cosadd -= 1;

   if (add > 2 * M_PI) {
      k = (1 + add - 2*M_PI);
      xf->disc2_cosadd *= k;
      xf->disc2_sinadd *= k;
   }

   if (add < -2 * M_PI) {
      k = (1 + add + 2*M_PI);
      xf->disc2_cosadd *= k;
      xf->disc2_sinadd *= k;
   }

}

static void supershape_precalc(flam3_xform *xf) {
   xf->supershape_pm_4 = xf->supershape_m / 4.0;
   xf->supershape_pneg1_n1 = -1.0 / xf->supershape_n1;
}

void prepare_xform_fn_ptrs(flam3_genome *cp, randctx *rc) {

   double d;
   int i,j,totnum;

   /* Loop over valid xforms */
   for (i = 0; i < cp->num_xforms; i++) {
      d = cp->xform[i].density;
      if (d < 0.0) {
         fprintf(stderr, "xform %d weight must be non-negative, not %g.\n",i,d);
         exit(1);
      }

      if (i != cp->final_xform_index && d == 0.0)
         continue;

      totnum = 0;

      cp->xform[i].vis_adjusted = adjust_percentage(cp->xform[i].visibility);

      cp->xform[i].precalc_angles_flag=0;
      cp->xform[i].precalc_atan_xy_flag=0;
      cp->xform[i].precalc_atan_yx_flag=0;
      cp->xform[i].has_preblur=0;

      for (j = 0; j < flam3_nvariations; j++) {

         if (cp->xform[i].var[j]!=0) {

            cp->xform[i].varFunc[totnum] = j;
            cp->xform[i].active_var_weights[totnum] = cp->xform[i].var[j];

            if (j==VAR_POLAR) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_HANDKERCHIEF) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_HEART) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_DISC) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_SPIRAL) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_HYPERBOLIC) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_DIAMOND) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_EX) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_JULIA) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_POWER) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_RINGS) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_FAN) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_BLOB) {
               cp->xform[i].precalc_atan_xy_flag=1;
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_FAN2) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_RINGS2) {
               cp->xform[i].precalc_angles_flag=1;
            } else if (j==VAR_JULIAN) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_JULIASCOPE) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_RADIAL_BLUR) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_NGON) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_DISC2) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_SUPER_SHAPE) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_FLOWER) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_CONIC) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_CPOW) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_ESCHER) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_PRE_BLUR) {
               cp->xform[i].has_preblur=cp->xform[i].var[j];
            } else if (j==VAR_POLAR2) {
               cp->xform[i].precalc_atan_xy_flag=1;
            } else if (j==VAR_WEDGE) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_WEDGE_JULIA) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_WEDGE_SPH) {
               cp->xform[i].precalc_atan_yx_flag=1;
            } else if (j==VAR_WHORL) {
               cp->xform[i].precalc_atan_yx_flag=1;
            }
            
            totnum++;
         }
      }

      cp->xform[i].num_active_vars = totnum;

   }
}

FILE *cout = NULL;
void write_color(double c) {
    if (NULL == cout) {
   cout = fopen("cout.txt", "w");
    }
    fprintf(cout, "%.30f\n", c);
}

unsigned short * flam3_create_xform_distrib(flam3_genome *cp) {

   /* Xform distrib is created in this function             */   
   int numrows;
   int dist_row,i;
   unsigned short *xform_distrib;
   
   numrows = cp->num_xforms - (cp->final_xform_index>=0) + 1;
   xform_distrib = calloc(numrows*CHOOSE_XFORM_GRAIN,sizeof(unsigned short));
   
   /* First, set up the first row of the xform_distrib (raw weights) */
   flam3_create_chaos_distrib(cp, -1, xform_distrib);
   
   /* Check for non-unity chaos */
   cp->chaos_enable = 1 - flam3_check_unity_chaos(cp);
   
   if (cp->chaos_enable) {
   
      /* Now set up a row for each of the xforms */
      dist_row = 0;
      for (i=0;i<cp->num_xforms;i++) {
      
         if (cp->final_xform_index == i)
            continue;
         else
            dist_row++;
         
         flam3_create_chaos_distrib(cp, i, &(xform_distrib[CHOOSE_XFORM_GRAIN*(dist_row)]));
      }
   }
   
   return(xform_distrib);
}

int flam3_check_unity_chaos(flam3_genome *cp) {

   int i,j;
   int num_std;
   int unity=1;
   num_std = cp->num_xforms - (cp->final_xform_index >= 0);
   
   for (i=0;i<num_std;i++) {
      for (j=0;j<num_std;j++) {
         if (cp->chaos[i][j] < 1.0)
            unity=0;
      }
   }
   
   return(unity);
}

void flam3_create_chaos_distrib(flam3_genome *cp, int xi, unsigned short *xform_distrib) {

   /* Xform distrib is a preallocated array of CHOOSE_XFORM_GRAIN chars */
   /* address of array is passed in, contents are modified              */
   double t,r,dr;
   int i,j;
   int num_std;
   
   //fprintf(stdout,"storing at %ld\n",xform_distrib);
   
   num_std = cp->num_xforms - (cp->final_xform_index >= 0);
   
   dr = 0.0;
   for (i = 0; i < num_std; i++) {
      double d = cp->xform[i].density;
                     
      if (xi>=0)
         d *= cp->chaos[xi][i];
         
      //fprintf(stdout,"%f ",d);
      if (d < 0.0) {
         fprintf(stderr, "xform weight must be non-negative, not %g.\n", d);
         exit(1);
      }         
      
      dr += d;
   }
   
   //fprintf(stdout,"dr=%f\n",dr);
   
   if (dr == 0.0) {
      fprintf(stderr, "cannot iterate empty flame.\n");
      exit(1);
   }
   
   dr = dr / CHOOSE_XFORM_GRAIN;

   j = 0;
   t = cp->xform[0].density;
   if (xi>=0)
     t *= cp->chaos[xi][0];
   r = 0.0;
   for (i = 0; i < CHOOSE_XFORM_GRAIN; i++) {
      while (r >= t) {
         j++;

         if (xi>=0)
            t += cp->xform[j].density*cp->chaos[xi][j];
         else
            t += cp->xform[j].density;
         
      }
      //fprintf(stdout,"%d ",j);
      xform_distrib[i] = j;
      r += dr;
   }
   //fprintf(stdout,"\n---\n");
}
/*
 * run the function system described by CP forward N generations.  store
 * the N resulting 4-vectors in SAMPLES.  the initial point is passed in
 * SAMPLES[0..3].  ignore the first FUSE iterations.
 */


int flam3_iterate(flam3_genome *cp, int n, int fuse,  double *samples, unsigned short *xform_distrib, randctx *rc) {
   int i;
   double p[4], q[4];
   int consec = 0;
   int badvals = 0;
   int lastxf=0;
   int fn;
   
   p[0] = samples[0];
   p[1] = samples[1];
   p[2] = samples[2];
   p[3] = samples[3];

   /* Perform precalculations */   
   for (i=0;i<cp->num_xforms;i++)
      xform_precalc(cp,i);

   for (i = -4*fuse; i < 4*n; i+=4) {
      if (cp->chaos_enable)
         fn = xform_distrib[ lastxf*CHOOSE_XFORM_GRAIN + (((unsigned)irand(rc)) % CHOOSE_XFORM_GRAIN)];
      else
         fn = xform_distrib[ ((unsigned)irand(rc)) % CHOOSE_XFORM_GRAIN ];
      
      
      if (1) {
         if (apply_xform(cp, fn, p, q, rc)>0) {
            //fprintf(stdout,"bad result from xform %d\n",fn);
            consec ++;
            badvals ++;
            if (consec<5) {
               p[0] = q[0];
               p[1] = q[1];
               p[2] = q[2];
               p[3] = q[3];
               i -= 4;
               continue;
            } else
               consec = 0;
         } else
            consec = 0;
      } else {
         apply_xform(cp, fn, p, q, rc);
      }

      /* Store the last used transform */
      lastxf = fn+1;

//      fprintf(stderr,"lastxf (post)=%d\n",lastxf);

      p[0] = q[0];
      p[1] = q[1];
      p[2] = q[2];
      p[3] = q[3];

      if (cp->final_xform_enable == 1) {
         apply_xform(cp, cp->final_xform_index, p, q, rc);
      }

      /* if fuse over, store it */
      if (i >= 0) {
         samples[i] = q[0];
         samples[i+1] = q[1];
         samples[i+2] = q[2];
         samples[i+3] = q[3];
      }
   }
   
   return(badvals);
}

static int apply_xform(flam3_genome *cp, int fn, double *p, double *q, randctx *rc)
{
   flam3_iter_helper f;
   int var_n;
   double next_color,s,s1;

   f.rc = rc;

   s = cp->xform[fn].color_speed;
   s1 = 0.5 - 0.5 * s;

   next_color = (p[2] + cp->xform[fn].color) * s1 + s * p[2];
   q[2] = next_color;
   q[3] = cp->xform[fn].vis_adjusted;
//   if (cp->xform[fn].vis_adjusted<1.0) {
//      q[3] = (flam3_random_isaac_01(rc) < cp->xform[fn].vis_adjusted);
//   } else {
//      q[3] = 1.0;
//   }

   f.tx = cp->xform[fn].c[0][0] * p[0] + cp->xform[fn].c[1][0] * p[1] + cp->xform[fn].c[2][0];
   f.ty = cp->xform[fn].c[0][1] * p[0] + cp->xform[fn].c[1][1] * p[1] + cp->xform[fn].c[2][1];

   /* Pre-xforms go here, and modify the f.tx and f.ty values */
   if (cp->xform[fn].has_preblur!=0.0)
      var67_pre_blur(&f, cp->xform[fn].has_preblur);

   /* Always calculate sumsq and sqrt */
   f.precalc_sumsq = f.tx*f.tx + f.ty*f.ty;
   f.precalc_sqrt = sqrt(f.precalc_sumsq);

   /* Check to see if we can precalculate any parts */
   /* Precalculate atanxy, sin, cos */
   if (cp->xform[fn].precalc_atan_xy_flag > 0) {
      f.precalc_atan = atan2(f.tx,f.ty);
   }
   
   if (cp->xform[fn].precalc_angles_flag > 0) {
      f.precalc_sina = f.tx / f.precalc_sqrt;
      f.precalc_cosa = f.ty / f.precalc_sqrt;
   }

   /* Precalc atanyx */
   if (cp->xform[fn].precalc_atan_yx_flag > 0) {
      f.precalc_atanyx = atan2(f.ty,f.tx);
   }

   f.p0 = 0.0;
   f.p1 = 0.0;
   f.xform = &(cp->xform[fn]);

   
   for (var_n=0; var_n < cp->xform[fn].num_active_vars; var_n++) {
//      (*cp->xform[fn].varFunc[var_n])(&f, cp->xform[fn].active_var_weights[var_n]);
//   }
//

      switch (cp->xform[fn].varFunc[var_n]) {
 
         case (VAR_LINEAR):
            var0_linear(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_SINUSOIDAL):
                var1_sinusoidal(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_SPHERICAL):
                var2_spherical(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_SWIRL):
                var3_swirl(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_HORSESHOE):
                var4_horseshoe(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_POLAR): 
                var5_polar(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_HANDKERCHIEF):
                var6_handkerchief(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_HEART):
                var7_heart(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_DISC):
                var8_disc(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_SPIRAL):
                var9_spiral(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_HYPERBOLIC):
                var10_hyperbolic(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_DIAMOND):
                var11_diamond(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_EX):
                var12_ex(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_JULIA): 
                var13_julia(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_BENT):
                var14_bent(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_WAVES):
                var15_waves(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_FISHEYE): 
                var16_fisheye(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_POPCORN):
                var17_popcorn(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_EXPONENTIAL):
                var18_exponential(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_POWER): 
                var19_power(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_COSINE):
                var20_cosine(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_RINGS):
                var21_rings(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_FAN):
                var22_fan(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_BLOB):
                var23_blob(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_PDJ):
                var24_pdj(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_FAN2):
                var25_fan2(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_RINGS2): 
                var26_rings2(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_EYEFISH): 
                var27_eyefish(&f, cp->xform[fn].active_var_weights[var_n]); break;             
         case (VAR_BUBBLE):
                var28_bubble(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_CYLINDER):
                var29_cylinder(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_PERSPECTIVE):
                var30_perspective(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_NOISE):
                var31_noise(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_JULIAN): 
                var32_juliaN_generic(&f, cp->xform[fn].active_var_weights[var_n]); break;            
         case (VAR_JULIASCOPE):
                var33_juliaScope_generic(&f, cp->xform[fn].active_var_weights[var_n]);break;
         case (VAR_BLUR):
                var34_blur(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_GAUSSIAN_BLUR):
                var35_gaussian(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_RADIAL_BLUR):
                var36_radial_blur(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_PIE):
                var37_pie(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_NGON):
                var38_ngon(&f, cp->xform[fn].active_var_weights[var_n]); break;          
         case (VAR_CURL):
                var39_curl(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_RECTANGLES):
                var40_rectangles(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_ARCH):
                var41_arch(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_TANGENT):
                var42_tangent(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_SQUARE):
                var43_square(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_RAYS):
                var44_rays(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_BLADE): 
                var45_blade(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_SECANT2): 
                var46_secant2(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_TWINTRIAN): 
                var47_twintrian(&f, cp->xform[fn].active_var_weights[var_n]); break;               
         case (VAR_CROSS):
                var48_cross(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_DISC2):
                var49_disc2(&f, cp->xform[fn].active_var_weights[var_n]); break;            
         case (VAR_SUPER_SHAPE):
                var50_supershape(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_FLOWER):
                var51_flower(&f, cp->xform[fn].active_var_weights[var_n]); break;            
         case (VAR_CONIC):
                var52_conic(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_PARABOLA): 
                var53_parabola(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_BENT2): 
                var54_bent2(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_BIPOLAR): 
                var55_bipolar(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_BOARDERS): 
                var56_boarders(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_BUTTERFLY): 
                var57_butterfly(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_CELL): 
                var58_cell(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_CPOW): 
                var59_cpow(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_CURVE): 
                var60_curve(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_EDISC): 
                var61_edisc(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_ELLIPTIC): 
                var62_elliptic(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_ESCHER): 
                var63_escher(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_FOCI): 
                var64_foci(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_LAZYSUSAN): 
                var65_lazysusan(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_LOONIE): 
                var66_loonie(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_MODULUS): 
                var68_modulus(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_OSCILLOSCOPE): 
                var69_oscope(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_POLAR2): 
                var70_polar2(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_POPCORN2): 
                var71_popcorn2(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_SCRY): 
                var72_scry(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_SEPARATION): 
                var73_separation(&f, cp->xform[fn].active_var_weights[var_n]); break;              
         case (VAR_SPLIT):
                var74_split(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_SPLITS):
                var75_splits(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_STRIPES):
                var76_stripes(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_WEDGE):
                var77_wedge(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_WEDGE_JULIA):
                var78_wedge_julia(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_WEDGE_SPH):
                var79_wedge_sph(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_WHORL):
                var80_whorl(&f, cp->xform[fn].active_var_weights[var_n]); break;
         case (VAR_WAVES2):
                var81_waves2(&f, cp->xform[fn].active_var_weights[var_n]); break;
      }
//      if (badvalue(f.p0) || badvalue(f.p1)) {
//         fprintf(stderr,"%d\n",cp->xform[fn].varFunc[var_n]);
//         break;
//      }
   }
   /* apply the post transform */
   if (!id_matrix(cp->xform[fn].post)) {
      q[0] = cp->xform[fn].post[0][0] * f.p0 + cp->xform[fn].post[1][0] * f.p1 + cp->xform[fn].post[2][0];
      q[1] = cp->xform[fn].post[0][1] * f.p0 + cp->xform[fn].post[1][1] * f.p1 + cp->xform[fn].post[2][1];
   } else {
      q[0] = f.p0;
      q[1] = f.p1;
   }

   if (badvalue(q[0]) || badvalue(q[1])) {
      q[0] = flam3_random_isaac_11(rc);
      q[1] = flam3_random_isaac_11(rc);
      return(1);
   } else
      return(0);

}

void flam3_xform_preview(flam3_genome *cp, int xi, double range, int numvals, int depth, double *result, randctx *rc) {

   /* We will evaluate the 'xi'th xform 'depth' times, over the following values:           */
   /* x in [-range : range], y in [-range : range], with 2* (2*numvals+1)^2 values returned */ 
   double p[4];
   double incr;
   int outi;
   int xx,yy,dd;
   
   /* Prepare the function pointers */
   prepare_xform_fn_ptrs(cp,rc);
   
   /* Calculate increment */
   incr = range / (double)numvals;
   
   /* Perform precalculations */
   xform_precalc(cp,xi);
   
   outi=0;
   
   /* Loop over the grid */
   for (xx=-numvals;xx<=numvals;xx++) {
      for (yy=-numvals;yy<=numvals;yy++) {
      
         /* Calculate the input coordinates */
         p[0] = (double)xx * incr;
         p[1] = (double)yy * incr;
         
         /* Loop over the depth */
         for (dd=0;dd<depth;dd++)
            apply_xform(cp, xi, p, p, rc);
         
         result[outi] = p[0];
         result[outi+1] = p[1];
         
         outi += 2;
      }
   }
}         
 
static void xform_precalc(flam3_genome *cp, int xi) {

   perspective_precalc(&(cp->xform[xi]));
   juliaN_precalc(&(cp->xform[xi]));
   juliaScope_precalc(&(cp->xform[xi]));
   radial_blur_precalc(&(cp->xform[xi]));
   waves_precalc(&(cp->xform[xi]));
   disc2_precalc(&(cp->xform[xi]));
   supershape_precalc(&(cp->xform[xi]));
   wedgeJulia_precalc(&(cp->xform[xi]));   
}   

static double adjust_percentage(double in) {

   double out;

   if (in==0.0)
      out = 0.0;
   else
      out = pow(10.0, -log(1.0/in)/log(2) );

   return(out);
}

/* correlation dimension, after clint sprott.
   computes slope of the correlation sum at a size scale
   the order of 2% the size of the attractor or the camera. */
double flam3_dimension(flam3_genome *cp, int ntries, int clip_to_camera) {
  double fd;
  double *hist;
  double bmin[2];
  double bmax[2];
  double d2max;
  int lp;
  long int default_isaac_seed = (long int)time(0);
  randctx rc;
  int i, n1=0, n2=0, got, nclipped;

  /* Set up the isaac rng */
  for (lp = 0; lp < RANDSIZ; lp++)
     rc.randrsl[lp] = default_isaac_seed;

  irandinit(&rc,1);

  if (ntries < 2) ntries = 3000*1000;

  if (clip_to_camera) {
    double scale, ppux, corner0, corner1;
    scale = pow(2.0, cp->zoom);
    ppux = cp->pixels_per_unit * scale;
    corner0 = cp->center[0] - cp->width / ppux / 2.0;
    corner1 = cp->center[1] - cp->height / ppux / 2.0;
    bmin[0] = corner0;
    bmin[1] = corner1;
    bmax[0] = corner0 + cp->width  / ppux;
    bmax[1] = corner1 + cp->height / ppux;
  } else {
    flam3_estimate_bounding_box(cp, 0.0, 0, bmin, bmax, &rc);
  }

  d2max =
    (bmax[0] - bmin[0]) * (bmax[0] - bmin[0]) +
    (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]);

  //  fprintf(stderr, "d2max=%g %g %g %g %g\n", d2max,
  //  bmin[0], bmin[1], bmax[0], bmax[1]);

  hist = malloc(2 * ntries * sizeof(double));

  got = 0;
  nclipped = 0;
  while (got < 2*ntries) {
    double subb[4*SUB_BATCH_SIZE];
    int i4, clipped;
    unsigned short *xform_distrib;
    subb[0] = flam3_random_isaac_11(&rc);
    subb[1] = flam3_random_isaac_11(&rc);
    subb[2] = 0.0;
    subb[3] = 0.0;
    prepare_xform_fn_ptrs(cp,&rc);
    xform_distrib = flam3_create_xform_distrib(cp);
    flam3_iterate(cp, SUB_BATCH_SIZE, 20, subb, xform_distrib, &rc);
    free(xform_distrib);
    i4 = 0;
    for (i = 0; i < SUB_BATCH_SIZE; i++) {
      if (got == 2*ntries) break;
      clipped = clip_to_camera &&
   ((subb[i4] < bmin[0]) ||
    (subb[i4+1] < bmin[1]) ||
    (subb[i4] > bmax[0]) ||
    (subb[i4+1] > bmax[1]));
      if (!clipped) {
   hist[got] = subb[i4];
   hist[got+1] = subb[i4+1];
   got += 2;
      } else {
   nclipped++;
   if (nclipped > 10 * ntries) {
       fprintf(stderr, "warning: too much clipping, "
          "flam3_dimension giving up.\n");
       return sqrt(-1.0);
   }
      }
      i4 += 4;
    }
  }
  if (0)
    fprintf(stderr, "cliprate=%g\n", nclipped/(ntries+(double)nclipped));

  for (i = 0; i < ntries; i++) {
    int ri;
    double dx, dy, d2;
    double tx, ty;

    tx = hist[2*i];
    ty = hist[2*i+1];

    do {
      ri = 2 * (random() % ntries);
    } while (ri == i);

    dx = hist[ri] - tx;
    dy = hist[ri+1] - ty;
    d2 = dx*dx + dy*dy;
    if (d2 < 0.004 * d2max) n2++;
    if (d2 < 0.00004 * d2max) n1++;
  }

  fd = 0.434294 * log(n2 / (n1 - 0.5));

  if (0)
    fprintf(stderr, "n1=%d n2=%d\n", n1, n2);

  free(hist);
  return fd;
}

void flam3_colorhist(flam3_genome *cp, int num_batches, double *hist) {

  int lp,plp;
  int mycolor;
  long int default_isaac_seed = (long int)time(0);
  randctx rc;
  unsigned short *xform_distrib;
  double sub_batch[4*SUB_BATCH_SIZE];

  /* Set up the isaac rng */
  for (lp = 0; lp < RANDSIZ; lp++)
     rc.randrsl[lp] = default_isaac_seed;

  irandinit(&rc,1);
  
  memset(hist,0,256*sizeof(double));
  
  for (lp=0;lp<num_batches;lp++) {
  
    sub_batch[0] = flam3_random_isaac_11(&rc);
    sub_batch[1] = flam3_random_isaac_11(&rc);
    sub_batch[2] = 0;
    sub_batch[3] = 0;

    // get into the attractor
    prepare_xform_fn_ptrs(cp,&rc);
    xform_distrib = flam3_create_xform_distrib(cp);
    flam3_iterate(cp, SUB_BATCH_SIZE, 20, sub_batch, xform_distrib, &rc);
    free(xform_distrib);
    
    // histogram the colors in the sub_batch array
    for (plp=0;plp<4*SUB_BATCH_SIZE;plp+=4) {
      mycolor = (int)(sub_batch[plp+2]*CMAP_SIZE);
      if (mycolor<0) mycolor=0;
      if (mycolor>CMAP_SIZE_M1) mycolor=CMAP_SIZE_M1;
      
      hist[mycolor] += 1;
    }
  }
} 
  


double flam3_lyapunov(flam3_genome *cp, int ntries) {
  double p[4];
  double x, y;
  double xn, yn;
  double xn2, yn2;
  double dx, dy, r;
  double eps = 1e-5;
  int i;
  double sum = 0.0;
  unsigned short *xform_distrib;

  int lp;
  long int default_isaac_seed = (long int)time(0);
  randctx rc;

  /* Set up the isaac rng */
  for (lp = 0; lp < RANDSIZ; lp++)
     rc.randrsl[lp] = default_isaac_seed;

  irandinit(&rc,1);


  if (ntries < 1) ntries = 10000;

  for (i = 0; i < ntries; i++) {
    x = flam3_random_isaac_11(&rc);
    y = flam3_random_isaac_11(&rc);

    p[0] = x;
    p[1] = y;
    p[2] = 0.0;
    p[3] = 0.0;

    // get into the attractor
    prepare_xform_fn_ptrs(cp,&rc);
    xform_distrib = flam3_create_xform_distrib(cp);
    flam3_iterate(cp, 1, 20+(random()%10), p, xform_distrib, &rc);
    free(xform_distrib);

    x = p[0];
    y = p[1];

    // take one deterministic step
    srandom(i);

    prepare_xform_fn_ptrs(cp,&rc);
    xform_distrib = flam3_create_xform_distrib(cp);
    flam3_iterate(cp, 1, 0, p, xform_distrib, &rc);
    free(xform_distrib);

    xn = p[0];
    yn = p[1];

    do {
      dx = flam3_random_isaac_11(&rc);
      dy = flam3_random_isaac_11(&rc);
      r = sqrt(dx * dx + dy * dy);
    } while (r == 0.0);
    dx /= r;
    dy /= r;

    dx *= eps;
    dy *= eps;

    p[0] = x + dx;
    p[1] = y + dy;
    p[2] = 0.0;

    // take the same step but with eps
    srandom(i);
    prepare_xform_fn_ptrs(cp,&rc);
    xform_distrib = flam3_create_xform_distrib(cp);
    flam3_iterate(cp, 1, 0, p, xform_distrib, &rc);
    free(xform_distrib);

    xn2 = p[0];
    yn2 = p[1];

    r = sqrt((xn-xn2)*(xn-xn2) + (yn-yn2)*(yn-yn2));

    sum += log(r/eps);
  }
  return sum/(log(2.0)*ntries);
}

/* args must be non-overlapping */
static void mult_matrix(double s1[2][2], double s2[2][2], double d[2][2]) {
   d[0][0] = s1[0][0] * s2[0][0] + s1[1][0] * s2[0][1];
   d[1][0] = s1[0][0] * s2[1][0] + s1[1][0] * s2[1][1];
   d[0][1] = s1[0][1] * s2[0][0] + s1[1][1] * s2[0][1];
   d[1][1] = s1[0][1] * s2[1][0] + s1[1][1] * s2[1][1];
}

/* BY is angle in degrees */
void flam3_rotate(flam3_genome *cp, double by, int interpolation_type) {
   int i;
   for (i = 0; i < cp->num_xforms; i++) {
      double r[2][2];
      double T[2][2];
      double U[2][2];
      double dtheta = by * 2.0 * M_PI / 360.0;

      /* Don't rotate xforms with > 0 animate values */
      if (cp->xform[i].animate > 0.0)
         continue;

      if (cp->xform[i].padding == 1) {
         if (interpolation_type == flam3_inttype_compat) {
            /* gen 202 era flam3 did not rotate padded xforms */
            continue;
         } else if (interpolation_type == flam3_inttype_older) {
            /* not sure if 198 era flam3 rotated padded xforms */
            continue;
         } else if (interpolation_type == flam3_inttype_linear) {
            /* don't rotate for prettier symsings */
            continue;
         } else if (interpolation_type == flam3_inttype_log) {
            /* Current flam3: what do we prefer? */
            //continue;
         }
      }

      /* Do NOT rotate final xforms */
      if (cp->final_xform_enable==1 && cp->final_xform_index==i)
         continue;

      r[1][1] = r[0][0] = cos(dtheta);
      r[0][1] = sin(dtheta);
      r[1][0] = -r[0][1];
      T[0][0] = cp->xform[i].c[0][0];
      T[1][0] = cp->xform[i].c[1][0];
      T[0][1] = cp->xform[i].c[0][1];
      T[1][1] = cp->xform[i].c[1][1];
      mult_matrix(r, T, U);
      cp->xform[i].c[0][0] = U[0][0];
      cp->xform[i].c[1][0] = U[1][0];
      cp->xform[i].c[0][1] = U[0][1];
      cp->xform[i].c[1][1] = U[1][1];
   }
}

static double det_matrix(double s[2][2]) {
   return s[0][0] * s[1][1] - s[0][1] * s[1][0];
}

static int id_matrix(double s[3][2]) {
  return
    (s[0][0] == 1.0) &&
    (s[0][1] == 0.0) &&
    (s[1][0] == 0.0) &&
    (s[1][1] == 1.0) &&
    (s[2][0] == 0.0) &&
    (s[2][1] == 0.0);
}

static void copy_matrix(double to[3][2], double from[3][2]) {

    to[0][0] = from[0][0];
    to[0][1] = from[0][1];
    to[1][0] = from[1][0];
    to[1][1] = from[1][1];
    to[2][0] = from[2][0];
    to[2][1] = from[2][1];
}


static void clear_matrix(double m[3][2]) {
   m[0][0] = 0.0;
   m[0][1] = 0.0;
   m[1][0] = 0.0;
   m[1][1] = 0.0;
   m[2][0] = 0.0;
   m[2][1] = 0.0;
}

/* element-wise linear */
static void sum_matrix(double s, double m1[3][2], double m2[3][2]) {

   m2[0][0] += s * m1[0][0];
   m2[0][1] += s * m1[0][1];
   m2[1][0] += s * m1[1][0];
   m2[1][1] += s * m1[1][1];
   m2[2][0] += s * m1[2][0];
   m2[2][1] += s * m1[2][1];

}

void print_chaos(flam3_genome *cp) {
   int numstd = cp->num_xforms - (cp->final_xform_index >= 0);
   int i,j;
   
   fprintf(stderr,"-----\n");
   for (i=0;i<numstd;i++) {
      fprintf(stderr,"chaos %d: ",i);
      for (j=0;j<numstd;j++) {
         fprintf(stderr,"%f ",cp->chaos[i][j]);
      }
      fprintf(stderr,"\n");
   }
}


static void interpolate_cmap(flam3_palette cmap, double blend,
              int index0, double hue0, int index1, double hue1) {
                 
   flam3_palette p0,p1;
   int i, j;

   flam3_get_palette(index0, p0, hue0);
   flam3_get_palette(index1, p1, hue1);

   for (i = 0; i < 256; i++) {
      double t[5], s[5];
    
      rgb2hsv(p0[i].color, s);
      rgb2hsv(p1[i].color, t);
      
      s[3] = p0[i].color[3];
      t[3] = p1[i].color[3];
      
      s[4] = p0[i].index;
      t[4] = p1[i].index;
    
      for (j = 0; j < 5; j++)
         t[j] = ((1.0-blend) * s[j]) + (blend * t[j]);
         
      hsv2rgb(t, cmap[i].color);
      cmap[i].color[3] = t[3];
      cmap[i].index = t[4];      
   }
}

/*   tTime2 = tTime * tTime
     tTime3 = tTime2 * tTime
tpFinal = (((-tp1 + 3 * tp2 - 3 * tp3 + tp4) * tTime3)
                + ((2 * tp1 - 5 * tp2 + 4 * tp3 - tp4) * tTime2)
                + ((-tp1 + tp3) * tTime)
                + (2 * tp2))
                / 2
*/
static void interpolate_catmull_rom(flam3_genome cps[], double t, flam3_genome *result) {
    double t2 = t * t;
    double t3 = t2 * t;
    double cmc[4];

    cmc[0] = (2*t2 - t - t3) / 2;
    cmc[1] = (3*t3 - 5*t2 + 2) / 2;
    cmc[2] = (4*t2 - 3*t3 + t) / 2;
    cmc[3] = (t3 - t2) / 2;

    flam3_interpolate_n(result, 4, cps, cmc);
}


#define INTERP(x)  do { result->x = 0.0; \
   for (k = 0; k < ncp; k++) result->x += c[k] * cpi[k].x; } while(0)

#define INTERI(x)  do { double tt = 0.0; \
   for (k = 0; k < ncp; k++) tt += c[k] * cpi[k].x; \
   result->x = (int)rint(tt); } while(0)


/* all cpi and result must be aligned (have the same number of xforms,
   and have final xform in the same slot) */
void flam3_interpolate_n(flam3_genome *result, int ncp,
          flam3_genome *cpi, double *c) {
    int i, j, k, numstd;

   if (flam3_palette_interpolation_hsv == cpi[0].palette_interpolation) {
   
      for (i = 0; i < 256; i++) {
         double t[3], s[4];
         s[0] = s[1] = s[2] = s[3] = s[4] = 0.0;
         
         for (k = 0; k < ncp; k++) {
            rgb2hsv(cpi[k].palette[i].color, t);
            for (j = 0; j < 3; j++)
               s[j] += c[k] * t[j];
            
            s[3] += c[k] * cpi[k].palette[i].color[3];
            s[4] += c[k] * cpi[k].palette[i].index;
            
         }
       
         hsv2rgb(s, result->palette[i].color);
         result->palette[i].color[3] = s[3];
         result->palette[i].index = s[4];
       
         for (j = 0; j < 4; j++) {
           if (result->palette[i].color[j] < 0.0)
              result->palette[i].color[j] = 0.0;
           if (result->palette[i].color[j] > 1.0)
              result->palette[i].color[j] = 1.0;
         }
         
         if (result->palette[i].index < 0.0)
            result->palette[i].index = 0.0;
         if (result->palette[i].index > 255.0)
            result->palette[i].index = 255.0;
      }
   } else {
      /* Sweep - not the best option for float indices */
      for (i = 0; i < 256; i++) {
         j = (i < (256 * c[0])) ? 0 : 1;
         result->palette[i] = cpi[j].palette[i];
      }
   }

   result->palette_index = flam3_palette_random;
   result->symmetry = 0;
   result->spatial_filter_select = cpi[0].spatial_filter_select;
   result->temporal_filter_type = cpi[0].temporal_filter_type;
   result->palette_mode = cpi[0].palette_mode;

   result->interpolation_type = cpi[0].interpolation_type;
   INTERP(brightness);
   INTERP(contrast);
   INTERP(highlight_power);
   INTERP(gamma);
   INTERP(vibrancy);
   INTERP(hue_rotation);
   INTERI(width);
   INTERI(height);
   INTERI(spatial_oversample);
   INTERP(center[0]);
   INTERP(center[1]);
   INTERP(rot_center[0]);
   INTERP(rot_center[1]);
   INTERP(background[0]);
   INTERP(background[1]);
   INTERP(background[2]);
   INTERP(pixels_per_unit);
   INTERP(spatial_filter_radius);
   INTERP(temporal_filter_exp);
   INTERP(temporal_filter_width);
   INTERP(sample_density);
   INTERP(zoom);
   INTERP(rotate);
   INTERI(nbatches);
   INTERI(ntemporal_samples);
   INTERP(estimator);
   INTERP(estimator_minimum);
   INTERP(estimator_curve);
   INTERP(gam_lin_thresh);
   
   /* Interpolate the chaos array */
   numstd = cpi[0].num_xforms - (cpi[0].final_xform_index >= 0);
   for (i=0;i<numstd;i++) {
      for (j=0;j<numstd;j++) {
         INTERP(chaos[i][j]);
         if (result->chaos[i][j]<0) result->chaos[i][j]=0;
         if (result->chaos[i][j]>1) result->chaos[i][j]=1.0;
      }
   }

   for (i = 0; i < cpi[0].num_xforms; i++) {
      double td;
      int all_id;
      INTERP(xform[i].density);
      td = result->xform[i].density;
      result->xform[i].density = (td < 0.0) ? 0.0 : td;
      INTERP(xform[i].color);
      if (result->xform[i].color<0) result->xform[i].color=0;
      if (result->xform[i].color>1) result->xform[i].color=1;
      
      INTERP(xform[i].visibility);      
      INTERP(xform[i].color_speed);
      INTERP(xform[i].animate);
      INTERP(xform[i].blob_low);
      INTERP(xform[i].blob_high);
      INTERP(xform[i].blob_waves);
      INTERP(xform[i].pdj_a);
      INTERP(xform[i].pdj_b);
      INTERP(xform[i].pdj_c);
      INTERP(xform[i].pdj_d);
      INTERP(xform[i].fan2_x);
      INTERP(xform[i].fan2_y);
      INTERP(xform[i].rings2_val);
      INTERP(xform[i].perspective_angle);
      INTERP(xform[i].perspective_dist);
      INTERP(xform[i].juliaN_power);
      INTERP(xform[i].juliaN_dist);
      INTERP(xform[i].juliaScope_power);
      INTERP(xform[i].juliaScope_dist);
      INTERP(xform[i].radialBlur_angle);
      INTERP(xform[i].pie_slices);
      INTERP(xform[i].pie_rotation);
      INTERP(xform[i].pie_thickness);
      INTERP(xform[i].ngon_sides);
      INTERP(xform[i].ngon_power);
      INTERP(xform[i].ngon_circle);
      INTERP(xform[i].ngon_corners);
      INTERP(xform[i].curl_c1);
      INTERP(xform[i].curl_c2);
      INTERP(xform[i].rectangles_x);
      INTERP(xform[i].rectangles_y);
      INTERP(xform[i].amw_amp);
      INTERP(xform[i].disc2_rot);
      INTERP(xform[i].disc2_twist);
      INTERP(xform[i].supershape_rnd);
      INTERP(xform[i].supershape_m);
      INTERP(xform[i].supershape_n1);
      INTERP(xform[i].supershape_n2);
      INTERP(xform[i].supershape_n3);
      INTERP(xform[i].supershape_holes);
      INTERP(xform[i].flower_petals);
      INTERP(xform[i].flower_holes);
      INTERP(xform[i].conic_eccen);
      INTERP(xform[i].conic_holes);
      INTERP(xform[i].parabola_height);
      INTERP(xform[i].parabola_width);
      INTERP(xform[i].bent2_x);
      INTERP(xform[i].bent2_y);
      INTERP(xform[i].bipolar_shift);
      INTERP(xform[i].cell_size);
      INTERP(xform[i].cpow_r);
      INTERP(xform[i].cpow_i);
      INTERP(xform[i].cpow_power);
      INTERP(xform[i].curve_xamp);
      INTERP(xform[i].curve_yamp);
      INTERP(xform[i].curve_xlength);
      INTERP(xform[i].curve_ylength);
      INTERP(xform[i].escher_beta);
      INTERP(xform[i].lazysusan_x);
      INTERP(xform[i].lazysusan_y);
      INTERP(xform[i].lazysusan_twist);
      INTERP(xform[i].lazysusan_space);
      INTERP(xform[i].lazysusan_spin);
      INTERP(xform[i].modulus_x);
      INTERP(xform[i].modulus_y);
      INTERP(xform[i].oscope_separation);
      INTERP(xform[i].oscope_frequency);
      INTERP(xform[i].oscope_amplitude);
      INTERP(xform[i].oscope_damping);
      INTERP(xform[i].popcorn2_x);
      INTERP(xform[i].popcorn2_y);
      INTERP(xform[i].popcorn2_c);
      INTERP(xform[i].separation_x);
      INTERP(xform[i].separation_xinside);
      INTERP(xform[i].separation_y);
      INTERP(xform[i].separation_yinside);
      INTERP(xform[i].split_xsize);
      INTERP(xform[i].split_ysize);
      INTERP(xform[i].splits_x);
      INTERP(xform[i].splits_y);
      INTERP(xform[i].stripes_space);
      INTERP(xform[i].stripes_warp);
      INTERP(xform[i].wedge_angle);
      INTERP(xform[i].wedge_hole);
      INTERP(xform[i].wedge_count);
      INTERP(xform[i].wedge_swirl);
      INTERP(xform[i].wedge_julia_angle);
      INTERP(xform[i].wedge_julia_count);
      INTERP(xform[i].wedge_julia_power);
      INTERP(xform[i].wedge_julia_dist);
      INTERP(xform[i].wedge_sph_angle);
      INTERP(xform[i].wedge_sph_hole);
      INTERP(xform[i].wedge_sph_count);
      INTERP(xform[i].wedge_sph_swirl);
      INTERP(xform[i].whorl_inside);
      INTERP(xform[i].whorl_outside);
      INTERP(xform[i].waves2_scalex);
      INTERP(xform[i].waves2_scaley);
      INTERP(xform[i].waves2_freqx);
      INTERP(xform[i].waves2_freqy);

      for (j = 0; j < flam3_nvariations; j++)
         INTERP(xform[i].var[j]);

      if (flam3_inttype_log == cpi[0].interpolation_type) {
         int col;
         double cxmag[4][2];  // XXX why only 4? should be ncp
         double cxang[4][2];
         double cxtrn[4][2];

         /* affine part */
         clear_matrix(result->xform[i].c);
         convert_linear_to_polar(cpi,ncp,i,0,cxang,cxmag,cxtrn);
         interp_and_convert_back(c, ncp, i, cxang, cxmag, cxtrn,result->xform[i].c);

         /* post part */
         all_id = 1;
         for (k=0; k<ncp; k++)
            all_id &= id_matrix(cpi[k].xform[i].post);
         
         clear_matrix(result->xform[i].post);
         if (all_id) {
            result->xform[i].post[0][0] = 1.0;
            result->xform[i].post[1][1] = 1.0;
         } else {
            convert_linear_to_polar(cpi,ncp,i,1,cxang,cxmag,cxtrn);
            interp_and_convert_back(c, ncp, i, cxang, cxmag, cxtrn,result->xform[i].post);
         }
         

         if (0) {
            fprintf(stderr,"original coefs\n");
            for (k=0;k<3;k++) {
               for (col=0;col<2;col++) {
                  fprintf(stderr,"%f ",cpi[0].xform[i].c[k][col]);
               }
            }
            fprintf(stderr,"\n");
            
            fprintf(stderr,"new coefs\n");
            for (k=0;k<3;k++) {
               for (col=0;col<2;col++) {
                  fprintf(stderr,"%f ",result->xform[i].c[k][col]);
               }
            }
            fprintf(stderr,"\n");

         }
         
      } else {

         /* Interpolate c matrix & post */
         clear_matrix(result->xform[i].c);
         clear_matrix(result->xform[i].post);
         all_id = 1;
         for (k = 0; k < ncp; k++) {
            sum_matrix(c[k], cpi[k].xform[i].c, result->xform[i].c);
            sum_matrix(c[k], cpi[k].xform[i].post, result->xform[i].post);

            all_id &= id_matrix(cpi[k].xform[i].post);

         }
         if (all_id) {
            clear_matrix(result->xform[i].post);
            result->xform[i].post[0][0] = 1.0;
            result->xform[i].post[1][1] = 1.0;
         }
      }

   }
}

#define APPMOT(x)  do { addto->x += mot[i].x * motion_funcs(func,freq*blend); } while (0);

double motion_funcs(int funcnum, double timeval) {

   /* motion funcs should be cyclic, and equal to 0 at integral time values */
   /* abs peak values should be not be greater than 1                       */
   if (funcnum==MOTION_SIN) {
      return (sin(2.0*M_PI*timeval));
   } else if (funcnum==MOTION_TRIANGLE) {
      double fr = fmod(timeval,1.0);
      
      if (fr<0) fr+= 1.0;
      
      if (fr<=.25)
         fr = 4.0 * fr;
      else if (fr<=.75)
         fr = -4.0 * fr + 2.0;
      else
         fr = 4.0 * fr - 4.0;
     
      return(fr);
   } else if (funcnum==MOTION_COS) {
      return(1.0-cos(2.0*M_PI*timeval));
   }
   
}

void apply_motion_parameters(flam3_xform *xf, flam3_xform *addto, double blend) {

   int i,j,k;
   int freq;
   int func;
   flam3_xform* mot;
   
   mot = xf->motion;

   /* Loop over the motion elements and add their contribution to the original vals */
   for (i=0; i<xf->num_motion; i++) {   
   
      freq = mot->motion_freq;
      func = mot->motion_func;
      
      APPMOT(density); /* Must ensure > 0 after all is applied */
      APPMOT(color); /* Must ensure [0,1] after all is applied */
      
      APPMOT(visibility);      
      APPMOT(color_speed);
      APPMOT(animate);
      APPMOT(blob_low);
      APPMOT(blob_high);
      APPMOT(blob_waves);
      APPMOT(pdj_a);
      APPMOT(pdj_b);
      APPMOT(pdj_c);
      APPMOT(pdj_d);
      APPMOT(fan2_x);
      APPMOT(fan2_y);
      APPMOT(rings2_val);
      APPMOT(perspective_angle);
      APPMOT(perspective_dist);
      APPMOT(juliaN_power);
      APPMOT(juliaN_dist);
      APPMOT(juliaScope_power);
      APPMOT(juliaScope_dist);
      APPMOT(radialBlur_angle);
      APPMOT(pie_slices);
      APPMOT(pie_rotation);
      APPMOT(pie_thickness);
      APPMOT(ngon_sides);
      APPMOT(ngon_power);
      APPMOT(ngon_circle);
      APPMOT(ngon_corners);
      APPMOT(curl_c1);
      APPMOT(curl_c2);
      APPMOT(rectangles_x);
      APPMOT(rectangles_y);
      APPMOT(amw_amp);
      APPMOT(disc2_rot);
      APPMOT(disc2_twist);
      APPMOT(supershape_rnd);
      APPMOT(supershape_m);
      APPMOT(supershape_n1);
      APPMOT(supershape_n2);
      APPMOT(supershape_n3);
      APPMOT(supershape_holes);
      APPMOT(flower_petals);
      APPMOT(flower_holes);
      APPMOT(conic_eccen);
      APPMOT(conic_holes);
      APPMOT(parabola_height);
      APPMOT(parabola_width);
      APPMOT(bent2_x);
      APPMOT(bent2_y);
      APPMOT(bipolar_shift);
      APPMOT(cell_size);
      APPMOT(cpow_r);
      APPMOT(cpow_i);
      APPMOT(cpow_power);
      APPMOT(curve_xamp);
      APPMOT(curve_yamp);
      APPMOT(curve_xlength);
      APPMOT(curve_ylength);
      APPMOT(escher_beta);
      APPMOT(lazysusan_x);
      APPMOT(lazysusan_y);
      APPMOT(lazysusan_twist);
      APPMOT(lazysusan_space);
      APPMOT(lazysusan_spin);
      APPMOT(modulus_x);
      APPMOT(modulus_y);
      APPMOT(oscope_separation);
      APPMOT(oscope_frequency);
      APPMOT(oscope_amplitude);
      APPMOT(oscope_damping);
      APPMOT(popcorn2_x);
      APPMOT(popcorn2_y);
      APPMOT(popcorn2_c);
      APPMOT(separation_x);
      APPMOT(separation_xinside);
      APPMOT(separation_y);
      APPMOT(separation_yinside);
      APPMOT(split_xsize);
      APPMOT(split_ysize);
      APPMOT(splits_x);
      APPMOT(splits_y);
      APPMOT(stripes_space);
      APPMOT(stripes_warp);
      APPMOT(wedge_angle);
      APPMOT(wedge_hole);
      APPMOT(wedge_count);
      APPMOT(wedge_swirl);
      APPMOT(wedge_julia_angle);
      APPMOT(wedge_julia_count);
      APPMOT(wedge_julia_power);
      APPMOT(wedge_julia_dist);
      APPMOT(wedge_sph_angle);
      APPMOT(wedge_sph_hole);
      APPMOT(wedge_sph_count);
      APPMOT(wedge_sph_swirl);
      APPMOT(whorl_inside);
      APPMOT(whorl_outside);
      APPMOT(waves2_scalex);
      APPMOT(waves2_scaley);
      APPMOT(waves2_freqx);
      APPMOT(waves2_freqy);

      for (j = 0; j < flam3_nvariations; j++)
         APPMOT(var[j]);
         
      for (j=0; j<3; j++) {
         for (k=0; k<2; k++) {
            APPMOT(c[j][k]);
            APPMOT(post[j][k]);
         }
      }
         
   }
   
   /* Make sure certain params are within reasonable bounds */
   if (addto->color<0) addto->color=0;
   if (addto->color>1) addto->color=1;
   if (addto->density<0) addto->density=0;
   
}
      
      
   

void establish_asymmetric_refangles(flam3_genome *cp, int ncps) {

   int k, xfi, col;
   
   double cxang[4][2],d,c1[2];

   for (xfi=0; xfi<cp[0].num_xforms; xfi++) {
   
     /* Final xforms don't rotate regardless of their symmetry */
     if (cp[0].final_xform_enable==1 && xfi==cp[0].final_xform_index)
        continue;

     for (k=0; k<ncps;k++) {

          /* Establish the angle for each component */
          /* Should potentially functionalize */
          for (col=0;col<2;col++) {
          
               c1[0] = cp[k].xform[xfi].c[col][0];
               c1[1] = cp[k].xform[xfi].c[col][1];
               
               cxang[k][col] = atan2(c1[1],c1[0]);
          }
     }
      
     for (k=1; k<ncps; k++) {
     
          for (col=0;col<2;col++) {

               int sym0,sym1;
          int padsymflag;

               d = cxang[k][col]-cxang[k-1][col];

               /* Adjust to avoid the -pi/pi discontinuity */
               if (d > M_PI+EPS)
               cxang[k][col] -= 2*M_PI;
               else if (d < -(M_PI-EPS) )
               cxang[k][col] += 2*M_PI;

               /* If this is an asymmetric case, store the NON-symmetric angle    */
               /* Check them pairwise and store the reference angle in the second */
               /* to avoid overwriting if asymmetric on both sides                */
               /* Depending on the interpolation type, treat padded xforms as symmetric */
               //padsymflag = (cp[k-1].interpolation_type==flam3_inttype_log);
               padsymflag = 0;
          
               sym0 = (cp[k-1].xform[xfi].animate>0 || (cp[k-1].xform[xfi].padding==1 && padsymflag));
               sym1 = (cp[k].xform[xfi].animate>0 || (cp[k].xform[xfi].padding==1 && padsymflag));

               if ( sym1 && !sym0 )
                  cp[k].xform[xfi].wind[col] = cxang[k-1][col] + 2*M_PI;
               else if ( sym0 && !sym1 )
                  cp[k].xform[xfi].wind[col] = cxang[k][col] + 2*M_PI;

          }
     }
   }
}



static void convert_linear_to_polar(flam3_genome *cp, int ncps, int xfi, int cflag, double cxang[4][2], double cxmag[4][2], double cxtrn[4][2]) {

   double c1[2],d,t,refang;
   int col,k;
   int zlm[2];

   for (k=0; k<ncps;k++) {

      /* Establish the angles and magnitudes for each component */
      /* Keep translation linear */
      zlm[0]=zlm[1]=0;
      for (col=0;col<2;col++) {
      
         if (cflag==0) {
            c1[0] = cp[k].xform[xfi].c[col][0];
            c1[1] = cp[k].xform[xfi].c[col][1];
            t = cp[k].xform[xfi].c[2][col];            
         } else {
            c1[0] = cp[k].xform[xfi].post[col][0];
            c1[1] = cp[k].xform[xfi].post[col][1];
            t = cp[k].xform[xfi].post[2][col];
         }
         
         cxang[k][col] = atan2(c1[1],c1[0]);
         cxmag[k][col] = sqrt(c1[0]*c1[0] + c1[1]*c1[1]);
         
         if (cxmag[k][col]== 0.0)
            zlm[col]=1;
         
         cxtrn[k][col] = t;
      }
      
      if (zlm[0]==1 && zlm[1]==0)
         cxang[k][0] = cxang[k][1];
      else if (zlm[0]==0 && zlm[1]==1)
         cxang[k][1] = cxang[k][0];
      
   }
   
   /* Make sure the rotation is the shorter direction around the circle */
   /* by adjusting each angle in succession, and rotate clockwise if 180 degrees */
   {
      for (col=0; col<2; col++) {
         for (k=1;k<ncps;k++) {

            /* Adjust angles differently if we have an asymmetric case */   
            if (cp[k].xform[xfi].wind[col]>0 && cflag==0) {

               /* Adjust the angles to make sure that it's within wind:wind+2pi */
               refang = cp[k].xform[xfi].wind[col] - 2*M_PI;

               /* Make sure both angles are within [refang refang+2*pi] */
               while(cxang[k-1][col] < refang)
                    cxang[k-1][col] += 2*M_PI;
               
               while(cxang[k-1][col] > refang + 2*M_PI)
                    cxang[k-1][col] -= 2*M_PI;
                    
               while(cxang[k][col] < refang)
                    cxang[k][col] += 2*M_PI;
               
               while(cxang[k][col] > refang + 2*M_PI)
                    cxang[k][col] -= 2*M_PI;

            } else {

          /* Normal way of adjusting angles */
               d = cxang[k][col]-cxang[k-1][col];
         
          /* Adjust to avoid the -pi/pi discontinuity */
          if (d > M_PI+EPS)
                  cxang[k][col] -= 2*M_PI;
          else if (d < -(M_PI-EPS) ) /* Forces clockwise rotation at 180 */
             cxang[k][col] += 2*M_PI;
       }

            
         }
      }
   }
}

static void interp_and_convert_back(double *c, int ncps, int xfi, double cxang[4][2], double cxmag[4][2], double cxtrn[4][2],double store_array[3][2]) {

   int i,col;
   
   double accang[2],accmag[2];
   double expmag;
   int accmode[2];
   
   accang[0] = 0.0;
   accang[1] = 0.0;
   accmag[0] = 0.0;
   accmag[1] = 0.0;

   accmode[0]=accmode[1]=0;
   
   /* accumulation mode defaults to logarithmic, but in special */
   /* cases we want to switch to linear accumulation            */
   for (col=0; col<2; col++) {
      for (i=0; i<ncps; i++) {
         if (log(cxmag[i][col])<-10)
            accmode[col]=1; // Mode set to linear interp
      }
   }
   
   for (i=0; i<ncps; i++) {
      for (col=0; col<2; col++) {
      
         accang[col] += c[i] * cxang[i][col];
         
         if (accmode[col]==0)
            accmag[col] += c[i] * log(cxmag[i][col]);
         else 
            accmag[col] += c[i] * (cxmag[i][col]);
            
         /* translation is ready to go */
         store_array[2][col] += c[i] * cxtrn[i][col];
      }
   }
   
   /* Convert the angle back to rectangular */
   for (col=0;col<2;col++) {
      if (accmode[col]==0)
          expmag = exp(accmag[col]);
     else
          expmag = accmag[col];
      
      store_array[col][0] = expmag * cos(accang[col]);
      store_array[col][1] = expmag * sin(accang[col]);
   }
   
}
   
void flam3_align(flam3_genome *dst, flam3_genome *src, int nsrc) {
   int i, tfx, tnx, max_nx = 0, max_fx = 0;
   int already_aligned=1;
   int xf,j;
   int ii,fnd;
   double normed;
   
   max_nx = src[0].num_xforms - (src[0].final_xform_index >= 0);
   max_fx = src[0].final_xform_enable;
   
   for (i = 1; i < nsrc; i++) {
      tnx = src[i].num_xforms - (src[i].final_xform_index >= 0);
      if (max_nx != tnx) {
         already_aligned = 0;
         if (tnx > max_nx) max_nx = tnx;
      }
      
      tfx = src[i].final_xform_enable;
      if (max_fx != tfx) {
         already_aligned = 0;
         max_fx |= tfx;
      }
   }

   /* Pad the cps to equal xforms */
   for (i = 0; i < nsrc; i++) {
      flam3_copyx(&dst[i], &src[i], max_nx, max_fx);
   }
      
   /* Check to see if there's a parametric variation present in one xform   */
   /* but not in an aligned xform.  If this is the case, use the parameters */
   /* from the xform with the variation as the defaults for the blank one.  */

   /* Skip if this genome is compatibility mode */
   if (dst[i].interpolation_type == flam3_inttype_compat ||
       dst[i].interpolation_type == flam3_inttype_older)
      return;
   
   /* All genomes will have the same number of xforms at this point */
   /* num = max_nx + max_fx */
   for (i = 0; i<nsrc; i++) {

       for (xf = 0; xf<max_nx+max_fx; xf++) {
                  
          /* Loop over the variations to see which of them are set to 0 */
          /* Note that there are no parametric variations < 23 */
          for (j = 23; j < flam3_nvariations; j++) {
         
              if (dst[i].xform[xf].var[j]==0) {
            
                 if (i>0) {
                              
                    /* Check to see if the prior genome's xform is populated */
                    if (dst[i-1].xform[xf].var[j] != 0) {
                  
                       /* Copy the prior genome's parameters and continue */
                       flam3_copy_params(&(dst[i].xform[xf]), &(dst[i-1].xform[xf]), j);
                       continue;
                    }

                 } else if (i<nsrc-1) {

                    /* Check to see if the next genome's xform is populated */
                    if (dst[i+1].xform[xf].var[j] != 0) {
                  
                       /* Copy the next genome's parameters and continue */
                       flam3_copy_params(&(dst[i].xform[xf]), &(dst[i+1].xform[xf]), j);
                       continue;
                    }
                 }
              }
          } /* variations */

          if (dst[i].xform[xf].padding == 1 && !already_aligned) {
         
             /* This is a new xform.  Let's see if we can choose a better 'identity' xform. */
             /* Check the neighbors to see if any of these variations are used: */
             /* rings2, fan2, blob, perspective, julian, juliascope, ngon, curl, super_shape, split */
             /* If so, we can use a better starting point for these */
            
             /* Remove linear from the list */
             dst[i].xform[xf].var[0] = 0.0;
            
             /* Look through all of the 'companion' xforms to see if we get a match on any of these */
             fnd=0;

             /* Only do the next substitution for log interpolation */
             if ( (i==0 && dst[i].interpolation_type == flam3_inttype_log)
                  || (i>0 && dst[i-1].interpolation_type==flam3_inttype_log) ) {

             for (ii=-1; ii<=1; ii+=2) {

                /* Skip if out of bounds */
                if (i+ii<0 || i+ii>=nsrc)
                   continue;
                  
                /* Skip if this is also padding */
                if (dst[i+ii].xform[xf].padding==1)
                   continue;

                /* Spherical / Ngon (trumps all others due to holes)       */
                /* Interpolate these against a 180 degree rotated identity */
                /* with weight -1.                                         */
                /* Added JULIAN/JULIASCOPE to get rid of black wedges      */
                if (dst[i+ii].xform[xf].var[VAR_SPHERICAL]>0 ||
                      dst[i+ii].xform[xf].var[VAR_NGON]>0 || 
                      dst[i+ii].xform[xf].var[VAR_JULIAN]>0 || 
                      dst[i+ii].xform[xf].var[VAR_JULIASCOPE]>0 ||
                      dst[i+ii].xform[xf].var[VAR_POLAR]>0 ||
                      dst[i+ii].xform[xf].var[VAR_WEDGE_SPH]>0 ||
                      dst[i+ii].xform[xf].var[VAR_WEDGE_JULIA]>0) {
                 
                   dst[i].xform[xf].var[VAR_LINEAR] = -1.0;
                   /* Set the coefs appropriately */
                   dst[i].xform[xf].c[0][0] = -1.0;
                   dst[i].xform[xf].c[0][1] = 0.0;
                   dst[i].xform[xf].c[1][0] = 0.0;
                   dst[i].xform[xf].c[1][1] = -1.0;
                   dst[i].xform[xf].c[2][0] = 0.0;
                   dst[i].xform[xf].c[2][1] = 0.0;               
                   fnd=-1;
                }
             }

             }

             if (fnd==0) {

                for (ii=-1; ii<=1; ii+=2) {

                   /* Skip if out of bounds */
                   if (i+ii<0 || i+ii>=nsrc)
                      continue;
                     
                   /* Skip if also padding */
                   if (dst[i+ii].xform[xf].padding==1)
                      continue;

                   /* Rectangles */
                   if (dst[i+ii].xform[xf].var[VAR_RECTANGLES]>0) {
                      dst[i].xform[xf].var[VAR_RECTANGLES] = 1.0;
                      dst[i].xform[xf].rectangles_x = 0.0;
                      dst[i].xform[xf].rectangles_y = 0.0;
                      fnd++;
                   }

                   /* Rings 2 */
                   if (dst[i+ii].xform[xf].var[VAR_RINGS2]>0) {
                      dst[i].xform[xf].var[VAR_RINGS2] = 1.0;
                      dst[i].xform[xf].rings2_val = 0.0;
                      fnd++;
                   }
                  
                   /* Fan 2 */
                   if (dst[i+ii].xform[xf].var[VAR_FAN2]>0) {
                      dst[i].xform[xf].var[VAR_FAN2] = 1.0;
                      dst[i].xform[xf].fan2_x = 0.0;
                      dst[i].xform[xf].fan2_y = 0.0;
                      fnd++;
                   }
               
                   /* Blob */
                   if (dst[i+ii].xform[xf].var[VAR_BLOB]>0) {
                      dst[i].xform[xf].var[VAR_BLOB] = 1.0;
                      dst[i].xform[xf].blob_low = 1.0;
                      dst[i].xform[xf].blob_high = 1.0;
                      dst[i].xform[xf].blob_waves = 1.0;
                      fnd++;
                   }
               
                   /* Perspective */
                   if (dst[i+ii].xform[xf].var[VAR_PERSPECTIVE]>0) {
                      dst[i].xform[xf].var[VAR_PERSPECTIVE] = 1.0;
                      dst[i].xform[xf].perspective_angle = 0.0;
                      /* Keep the perspective distance as-is */
                      fnd++;
                   }
               
                   /* Curl */
                   if (dst[i+ii].xform[xf].var[VAR_CURL]>0) {
                      dst[i].xform[xf].var[VAR_CURL] = 1.0;
                      dst[i].xform[xf].curl_c1 = 0.0;
                      dst[i].xform[xf].curl_c2 = 0.0;
                      fnd++;
                   }

                   /* Super-Shape */
                   if (dst[i+ii].xform[xf].var[VAR_SUPER_SHAPE]>0) {
                      dst[i].xform[xf].var[VAR_SUPER_SHAPE] = 1.0;
                      /* Keep supershape_m the same */
                      dst[i].xform[xf].supershape_n1 = 2.0;
                      dst[i].xform[xf].supershape_n2 = 2.0;
                      dst[i].xform[xf].supershape_n3 = 2.0;
                      dst[i].xform[xf].supershape_rnd = 0.0;
                      dst[i].xform[xf].supershape_holes = 0.0;
                      fnd++;
                   }
                }
             }

             /* If we didn't have any matches with those, */
             /* try the affine ones, fan and rings        */
             if (fnd==0) {
            
                for (ii=-1; ii<=1; ii+=2) {

                   /* Skip if out of bounds */
                   if (i+ii<0 || i+ii>=nsrc)
                      continue;                  

                   /* Skip if also a padding xform */
                   if (dst[i+ii].xform[xf].padding==1)
                      continue;
                     
                   /* Fan */
                   if (dst[i+ii].xform[xf].var[VAR_FAN]>0) {
                      dst[i].xform[xf].var[VAR_FAN] = 1.0;
                      fnd++;
                   }

                   /* Rings */
                   if (dst[i+ii].xform[xf].var[VAR_RINGS]>0) {
                      dst[i].xform[xf].var[VAR_RINGS] = 1.0;
                      fnd++;
                   }

                }
               
                if (fnd>0) {
                   /* Set the coefs appropriately */
                   dst[i].xform[xf].c[0][0] = 0.0;
                   dst[i].xform[xf].c[0][1] = 1.0;
                   dst[i].xform[xf].c[1][0] = 1.0;
                   dst[i].xform[xf].c[1][1] = 0.0;
                   dst[i].xform[xf].c[2][0] = 0.0;
                   dst[i].xform[xf].c[2][1] = 0.0;               
                }
             }
                                          
             /* If we still have no matches, switch back to linear */
             if (fnd==0)

                dst[i].xform[xf].var[VAR_LINEAR] = 1.0;

             else if (fnd>0) {

                /* Otherwise, go through and normalize the weights. */
                normed = 0.0;
                for (j = 0; j < flam3_nvariations; j++)
                   normed += dst[i].xform[xf].var[j];
                  
                for (j = 0; j < flam3_nvariations; j++)
                   dst[i].xform[xf].var[j] /= normed;

             }         
          }
       } /* xforms */
   } /* genomes */
                              
}


/*
 * create a control point that interpolates between the control points
 * passed in CPS.  CPS must be sorted by time.
 */
void flam3_interpolate(flam3_genome cps[], int ncps,
             double time, flam3_genome *result) {
   int i1, i2;
   double c[2];
   flam3_genome cpi[4];
   int cpsi,xfi;
   int numstd,n;

   if (1 == ncps) {
      flam3_copy(result, &(cps[0]));
      return;
   }
      
   if (cps[0].time >= time) {
      i1 = 0;
      i2 = 1;
   } else if (cps[ncps - 1].time <= time) {
      i1 = ncps - 2;
      i2 = ncps - 1;
   } else {
      i1 = 0;
      while (cps[i1].time < time)
         i1++;

      i1--;
      i2 = i1 + 1;

   }

   c[0] = (cps[i2].time - time) / (cps[i2].time - cps[i1].time);
   c[1] = 1.0 - c[0];

   memset(cpi, 0, 4*sizeof(flam3_genome));
   
   /* To interpolate the xforms, we will make copies of the source cps  */
   /* and ensure that they both have the same number before progressing */
   if (flam3_interpolation_linear == cps[i1].interpolation) {
       flam3_align(&cpi[0], &cps[i1], 2);
     
   } else {
       if (0 == i1) {
      fprintf(stderr, "error: cannot use smooth interpolation on first segment.\n");
      exit(1);
       }
       if (ncps-1 == i2) {
      fprintf(stderr, "error: cannot use smooth interpolation on last segment.\n");
      exit(1);
       }
       flam3_align(&cpi[0], &cps[i1-1], 4);
   }
   
   /* Clear the destination cp */
   clear_cp(result, 1);
      
   if (cpi[0].final_xform_index >= 0) {
      flam3_add_xforms(result, cpi[0].num_xforms-1, 0, 0);
      flam3_add_xforms(result, 1, 0, 1);
   } else
      flam3_add_xforms(result, cpi[0].num_xforms, 0, 0);
   

   result->time = time;
   result->interpolation = flam3_interpolation_linear;
   result->interpolation_type = cpi[0].interpolation_type;
   result->palette_interpolation = flam3_palette_interpolation_hsv;

   if (flam3_interpolation_linear == cps[i1].interpolation) {
       flam3_interpolate_n(result, 2, cpi, c);
   } else {
       interpolate_catmull_rom(cpi, c[1], result);
       clear_cp(&(cpi[2]),0);
       clear_cp(&(cpi[3]),0);
//       free(cpi[2].xform);
//       free(cpi[3].xform);
   }
   
   clear_cp(&(cpi[0]),0);
   clear_cp(&(cpi[1]),0);
//   free(cpi[0].xform);
//   free(cpi[1].xform);

#if 0
   flam3_print(stdout, result, NULL);
   for (i = 0; i < sizeof(result->xform[0].post); i++) {
       printf("%d ", ((unsigned char *)result->xform[0].post)[i]);
   }
   printf("\n");
#endif
}

static int compare_xforms(const void *av, const void *bv) {
   flam3_xform *a = (flam3_xform *) av;
   flam3_xform *b = (flam3_xform *) bv;
   double aa[2][2];
   double bb[2][2];
   double ad, bd;

   aa[0][0] = a->c[0][0];
   aa[0][1] = a->c[0][1];
   aa[1][0] = a->c[1][0];
   aa[1][1] = a->c[1][1];
   bb[0][0] = b->c[0][0];
   bb[0][1] = b->c[0][1];
   bb[1][0] = b->c[1][0];
   bb[1][1] = b->c[1][1];
   ad = det_matrix(aa);
   bd = det_matrix(bb);

   if (a->color_speed < b->color_speed) return 1;
   if (a->color_speed > b->color_speed) return -1;
   if (a->color_speed) {
      if (ad < 0) return -1;
      if (bd < 0) return 1;
      ad = atan2(a->c[0][0], a->c[0][1]);
      bd = atan2(b->c[0][0], b->c[0][1]);
   }

   if (ad < bd) return -1;
   if (ad > bd) return 1;
   return 0;
}


static void initialize_xforms(flam3_genome *thiscp, int start_here) {

   int i,j;

   for (i = start_here ; i < thiscp->num_xforms ; i++) {
       thiscp->xform[i].padding = 0;
       thiscp->xform[i].density = 0.0;
       thiscp->xform[i].color_speed = 0.0;
       thiscp->xform[i].animate = 0.0;
       thiscp->xform[i].color = i&1;
       thiscp->xform[i].visibility = 1.0;
       thiscp->xform[i].var[0] = 1.0;
       thiscp->xform[i].motion_freq = 0;
       thiscp->xform[i].motion_func = 0;
       thiscp->xform[i].num_motion = 0;
       thiscp->xform[i].motion = NULL;
       for (j = 1; j < flam3_nvariations; j++)
          thiscp->xform[i].var[j] = 0.0;
       thiscp->xform[i].c[0][0] = 1.0;
       thiscp->xform[i].c[0][1] = 0.0;
       thiscp->xform[i].c[1][0] = 0.0;
       thiscp->xform[i].c[1][1] = 1.0;
       thiscp->xform[i].c[2][0] = 0.0;
       thiscp->xform[i].c[2][1] = 0.0;
       thiscp->xform[i].post[0][0] = 1.0;
       thiscp->xform[i].post[0][1] = 0.0;
       thiscp->xform[i].post[1][0] = 0.0;
       thiscp->xform[i].post[1][1] = 1.0;
       thiscp->xform[i].post[2][0] = 0.0;
       thiscp->xform[i].post[2][1] = 0.0;
       thiscp->xform[i].wind[0] = 0.0;
       thiscp->xform[i].wind[1] = 0.0;
       thiscp->xform[i].blob_low = 0.0;
       thiscp->xform[i].blob_high = 1.0;
       thiscp->xform[i].blob_waves = 1.0;
       thiscp->xform[i].pdj_a = 0.0;
       thiscp->xform[i].pdj_b = 0.0;
       thiscp->xform[i].pdj_c = 0.0;
       thiscp->xform[i].pdj_d = 0.0;
       thiscp->xform[i].fan2_x = 0.0;
       thiscp->xform[i].fan2_y = 0.0;
       thiscp->xform[i].rings2_val = 0.0;
       thiscp->xform[i].perspective_angle = 0.0;
       thiscp->xform[i].perspective_dist = 0.0;
       thiscp->xform[i].persp_vsin = 0.0;
       thiscp->xform[i].persp_vfcos = 0.0;
       thiscp->xform[i].radialBlur_angle = 0.0;
       thiscp->xform[i].disc2_rot = 0.0;
       thiscp->xform[i].disc2_twist = 0.0;
       thiscp->xform[i].disc2_sinadd = 0.0;
       thiscp->xform[i].disc2_cosadd = 0.0;
       thiscp->xform[i].disc2_timespi = 0.0;
       thiscp->xform[i].flower_petals = 0.0;
       thiscp->xform[i].flower_holes = 0.0;
       thiscp->xform[i].parabola_height = 0.0;
       thiscp->xform[i].parabola_width = 0.0;
       thiscp->xform[i].bent2_x = 1.0;
       thiscp->xform[i].bent2_y = 1.0;
       thiscp->xform[i].bipolar_shift = 0.0;
       thiscp->xform[i].cell_size = 1.0;
       thiscp->xform[i].cpow_r = 1.0;
       thiscp->xform[i].cpow_i = 0.0;
       thiscp->xform[i].cpow_power = 1.0;
       thiscp->xform[i].curve_xamp = 0.0;
       thiscp->xform[i].curve_yamp = 0.0;
       thiscp->xform[i].curve_xlength = 1.0;
       thiscp->xform[i].curve_ylength = 1.0;
       thiscp->xform[i].escher_beta = 0.0;
       thiscp->xform[i].lazysusan_space = 0.0;
       thiscp->xform[i].lazysusan_twist = 0.0;
       thiscp->xform[i].lazysusan_spin = 0.0;
       thiscp->xform[i].lazysusan_x = 0.0;
       thiscp->xform[i].lazysusan_y = 0.0;
       thiscp->xform[i].modulus_x = 0.0;
       thiscp->xform[i].modulus_y = 0.0;
       thiscp->xform[i].oscope_separation = 1.0;
       thiscp->xform[i].oscope_frequency = M_PI;
       thiscp->xform[i].oscope_amplitude = 1.0;
       thiscp->xform[i].oscope_damping = 0.0;
       thiscp->xform[i].popcorn2_c = 0.0;
       thiscp->xform[i].popcorn2_x = 0.0;
       thiscp->xform[i].popcorn2_y = 0.0;
       thiscp->xform[i].separation_x = 0.0;
       thiscp->xform[i].separation_xinside = 0.0;
       thiscp->xform[i].separation_y = 0.0;
       thiscp->xform[i].separation_yinside = 0.0;
       thiscp->xform[i].split_xsize = 0.0;
       thiscp->xform[i].split_ysize = 0.0;
       thiscp->xform[i].splits_x = 0.0;
       thiscp->xform[i].splits_y = 0.0;
       thiscp->xform[i].stripes_space = 0.0;
       thiscp->xform[i].stripes_warp = 0.0;
       thiscp->xform[i].wedge_angle = 0.0;
       thiscp->xform[i].wedge_hole = 0.0;
       thiscp->xform[i].wedge_count = 1.0;
       thiscp->xform[i].wedge_swirl = 0.0;
       thiscp->xform[i].wedge_sph_angle = 0.0;
       thiscp->xform[i].wedge_sph_hole = 0.0;
       thiscp->xform[i].wedge_sph_count = 1.0;
       thiscp->xform[i].wedge_sph_swirl = 0.0;

       thiscp->xform[i].wedge_julia_power = 1.0;
       thiscp->xform[i].wedge_julia_dist = 0.0;
       thiscp->xform[i].wedge_julia_count = 1.0;
       thiscp->xform[i].wedge_julia_angle = 0.0;
       thiscp->xform[i].wedgeJulia_cf = 0.0;
       thiscp->xform[i].wedgeJulia_cn = 0.5;
       thiscp->xform[i].wedgeJulia_rN = 1.0;
       thiscp->xform[i].whorl_inside = 0.0;
       thiscp->xform[i].whorl_outside = 0.0;
       
       thiscp->xform[i].waves2_scalex = 0.0;       
       thiscp->xform[i].waves2_scaley = 0.0;       
       thiscp->xform[i].waves2_freqx = 0.0;       
       thiscp->xform[i].waves2_freqy = 0.0;       
       
       thiscp->xform[i].juliaN_power = 1.0;
       thiscp->xform[i].juliaN_dist = 1.0;
       thiscp->xform[i].juliaN_rN = 1.0;
       thiscp->xform[i].juliaN_cn = 0.5;
       thiscp->xform[i].juliaScope_power = 1.0;
       thiscp->xform[i].juliaScope_dist = 1.0;
       thiscp->xform[i].juliaScope_rN = 1.0;
       thiscp->xform[i].juliaScope_cn = 0.5;
       thiscp->xform[i].radialBlur_spinvar = 0.0;
       thiscp->xform[i].radialBlur_zoomvar = 1.0;
       thiscp->xform[i].pie_slices = 6.0;
       thiscp->xform[i].pie_rotation = 0.0;
       thiscp->xform[i].pie_thickness = 0.5;
       thiscp->xform[i].ngon_sides = 5;
       thiscp->xform[i].ngon_power = 3;
       thiscp->xform[i].ngon_circle = 1;
       thiscp->xform[i].ngon_corners = 2;
       thiscp->xform[i].curl_c1 = 1.0;
       thiscp->xform[i].curl_c2 = 0.0;
       thiscp->xform[i].rectangles_x = 1.0;
       thiscp->xform[i].rectangles_y = 1.0;
       thiscp->xform[i].amw_amp = 1.0;
       thiscp->xform[i].supershape_rnd = 0.0;
       thiscp->xform[i].supershape_m = 0.0;
       thiscp->xform[i].supershape_n1 = 1.0;
       thiscp->xform[i].supershape_n2 = 1.0;
       thiscp->xform[i].supershape_n3 = 1.0;
       thiscp->xform[i].supershape_holes = 0.0;
       thiscp->xform[i].conic_eccen = 1.0;
       thiscp->xform[i].conic_holes = 0.0;


   }
}

void flam3_copy_params(flam3_xform *dest, flam3_xform *src, int varn) {

   /* We only want to copy param var coefs for this one */
   if (varn==VAR_BLOB) {
      /* Blob */
      dest->blob_low = src->blob_low;
      dest->blob_high = src->blob_high;
      dest->blob_waves = src->blob_waves;
   } else if (varn==VAR_PDJ) {
      /* PDJ */
      dest->pdj_a = src->pdj_a;
      dest->pdj_b = src->pdj_b;
      dest->pdj_c = src->pdj_c;
      dest->pdj_d = src->pdj_d;
   } else if (varn==VAR_FAN2) {
      /* Fan2 */
      dest->fan2_x = src->fan2_x;
      dest->fan2_y = src->fan2_y;
   } else if (varn==VAR_RINGS2) {
      /* Rings2 */
      dest->rings2_val = src->rings2_val;
   } else if (varn==VAR_PERSPECTIVE) {
      /* Perspective */
      dest->perspective_angle = src->perspective_angle;
      dest->perspective_dist = src->perspective_dist;
      dest->persp_vsin = src->persp_vsin;
      dest->persp_vfcos = src->persp_vfcos;
   } else if (varn==VAR_JULIAN) {
      /* Julia_N */
      dest->juliaN_power = src->juliaN_power;
      dest->juliaN_dist = src->juliaN_dist;
      dest->juliaN_rN = src->juliaN_rN;
      dest->juliaN_cn = src->juliaN_cn;
   } else if (varn==VAR_JULIASCOPE) {
      /* Julia_Scope */
      dest->juliaScope_power = src->juliaScope_power;
      dest->juliaScope_dist = src->juliaScope_dist;
      dest->juliaScope_rN = src->juliaScope_rN;
      dest->juliaScope_cn = src->juliaScope_cn;
   } else if (varn==VAR_RADIAL_BLUR) {
      /* Radial Blur */
      dest->radialBlur_angle = src->radialBlur_angle;
   } else if (varn==VAR_PIE) {
      /* Pie */
      dest->pie_slices = src->pie_slices;
      dest->pie_rotation = src->pie_rotation;
      dest->pie_thickness = src->pie_thickness;
   } else if (varn==VAR_NGON) {
      /* Ngon */
      dest->ngon_sides = src->ngon_sides;
      dest->ngon_power = src->ngon_power;
      dest->ngon_corners = src->ngon_corners;
      dest->ngon_circle = src->ngon_circle;
   } else if (varn==VAR_CURL) {
      /* Curl */
      dest->curl_c1 = src->curl_c1;
      dest->curl_c2 = src->curl_c2;
   } else if (varn==VAR_RECTANGLES) {
      /* Rect */
      dest->rectangles_x = src->rectangles_x;
      dest->rectangles_y = src->rectangles_y;
   } else if (varn==VAR_DISC2) {
      /* Disc2 */
      dest->disc2_rot = src->disc2_rot;
      dest->disc2_twist = src->disc2_twist;
   } else if (varn==VAR_SUPER_SHAPE) {
      /* Supershape */
      dest->supershape_rnd = src->supershape_rnd;
      dest->supershape_m = src->supershape_m;
      dest->supershape_n1 = src->supershape_n1;
      dest->supershape_n2 = src->supershape_n2;
      dest->supershape_n3 = src->supershape_n3;
      dest->supershape_holes = src->supershape_holes;
   } else if (varn==VAR_FLOWER) {
      /* Flower */
      dest->flower_petals = src->flower_petals;
      dest->flower_petals = src->flower_petals;
   } else if (varn==VAR_CONIC) {
      /* Conic */
      dest->conic_eccen = src->conic_eccen;
      dest->conic_holes = src->conic_holes;
   } else if (varn==VAR_PARABOLA) {
      /* Parabola */
      dest->parabola_height = src->parabola_height;
      dest->parabola_width = src->parabola_width;
   } else if (varn==VAR_BENT2) {
      /* Bent2 */
      dest->bent2_x = src->bent2_x;
      dest->bent2_y = src->bent2_y;
   } else if (varn==VAR_BIPOLAR) {
      /* Bipolar */
      dest->bipolar_shift = src->bipolar_shift;
   } else if (varn==VAR_CELL) {
      /* Cell */
      dest->cell_size = src->cell_size;
   } else if (varn==VAR_CPOW) {
      /* Cpow */
      dest->cpow_i = src->cpow_i;
      dest->cpow_r = src->cpow_r;
      dest->cpow_power = src->cpow_power;
   } else if (varn==VAR_CURVE) {
      /* Curve */
      dest->curve_xamp = src->curve_xamp;
      dest->curve_yamp = src->curve_yamp;
      dest->curve_xlength = src->curve_xlength;
      dest->curve_ylength = src->curve_ylength;
   } else if (varn==VAR_ESCHER) {
      /* Escher */
      dest->escher_beta = src->escher_beta;
   } else if (varn==VAR_LAZYSUSAN) {
      /* Lazysusan */
      dest->lazysusan_x = src->lazysusan_x;
      dest->lazysusan_y = src->lazysusan_y;
      dest->lazysusan_spin = src->lazysusan_spin;
      dest->lazysusan_space = src->lazysusan_space;
      dest->lazysusan_twist = src->lazysusan_twist;
   } else if (varn==VAR_MODULUS) {
      /* Modulus */
      dest->modulus_x = src->modulus_x;
      dest->modulus_y = src->modulus_y;
   } else if (varn==VAR_OSCILLOSCOPE) {
      /* Oscope */
      dest->oscope_separation = src->oscope_separation;
      dest->oscope_frequency = src->oscope_frequency;
      dest->oscope_amplitude = src->oscope_amplitude;
      dest->oscope_damping = src->oscope_damping;
   } else if (varn==VAR_POPCORN2) {
      /* Popcorn2 */
      dest->popcorn2_x = src->popcorn2_x;
      dest->popcorn2_y = src->popcorn2_y;
      dest->popcorn2_c = src->popcorn2_c;
   } else if (varn==VAR_SEPARATION) {
      /* Separation */
      dest->separation_x = src->separation_x;
      dest->separation_y = src->separation_y;
      dest->separation_xinside = src->separation_xinside;
      dest->separation_yinside = src->separation_yinside;
   } else if (varn==VAR_SPLIT) {
      /* Split */
      dest->split_xsize = src->split_xsize;
      dest->split_ysize = src->split_ysize;
   } else if (varn==VAR_SPLITS) {
      /* Splits */
      dest->splits_x = src->splits_x;
      dest->splits_y = src->splits_y;
   } else if (varn==VAR_STRIPES) {
      /* Stripes */
      dest->stripes_space = src->stripes_space;
      dest->stripes_warp = src->stripes_warp;
   } else if (varn==VAR_WEDGE) {
      /* Wedge */
      dest->wedge_angle = src->wedge_angle;
      dest->wedge_hole = src->wedge_hole;
      dest->wedge_count = src->wedge_count;
      dest->wedge_swirl = src->wedge_swirl;
   } else if (varn==VAR_WEDGE_JULIA) {
      /* Wedge_Julia */
      dest->wedge_julia_angle = src->wedge_julia_angle;
      dest->wedge_julia_count = src->wedge_julia_count;
      dest->wedge_julia_power = src->wedge_julia_power;
      dest->wedge_julia_dist = src->wedge_julia_dist;
      dest->wedgeJulia_cf = src->wedgeJulia_cf;
      dest->wedgeJulia_cn = src->wedgeJulia_cn;
      dest->wedgeJulia_rN = src->wedgeJulia_rN;
   } else if (varn==VAR_WEDGE_SPH) {
      /* Wedge_sph */
      dest->wedge_sph_angle = src->wedge_sph_angle;
      dest->wedge_sph_hole = src->wedge_sph_hole;
      dest->wedge_sph_count = src->wedge_sph_count;
      dest->wedge_sph_swirl = src->wedge_sph_swirl;      
   } else if (varn==VAR_WHORL) {
      /* whorl */
      dest->whorl_inside = src->whorl_inside;
      dest->whorl_outside = src->whorl_outside;
   } else if (varn==VAR_WAVES2) {
      /* waves2 */
      dest->waves2_scalex = src->waves2_scalex;
      dest->waves2_scaley = src->waves2_scaley;
      dest->waves2_freqx = src->waves2_freqx;
      dest->waves2_freqy = src->waves2_freqy;
   }
}

/* Motion support functions */
void flam3_add_motion_element(flam3_xform *xf) {

   /* Add one to the xform's count of motion elements */
   xf->num_motion++;
   
   /* Reallocate the motion storage to include the empty space */
   xf->motion = (struct xform *)realloc(xf->motion, xf->num_motion * sizeof(struct xform));
   
   /* Initialize the motion element */
   /* In this case, all elements should be set to 0 */
   memset( &(xf->motion[xf->num_motion-1]), 0, sizeof(struct xform));
   
}   
   

/* Xform support functions */
void flam3_add_xforms(flam3_genome *thiscp, int num_to_add, int interp_padding, int final_flag) {

   int i,j;
   int old_num = thiscp->num_xforms;
   int oldstd,numstd;
   flam3_xform tmp;
   
   oldstd = thiscp->num_xforms - (thiscp->final_xform_index >= 0);

//   if (thiscp->num_xforms > 0)
      thiscp->xform = (flam3_xform *)realloc(thiscp->xform, (thiscp->num_xforms + num_to_add) * sizeof(flam3_xform));
//   else
//      thiscp->xform = (flam3_xform *)malloc(num_to_add * sizeof(flam3_xform));

   thiscp->num_xforms += num_to_add;

   /* Initialize all the new xforms */
   initialize_xforms(thiscp, old_num);

   /* Set the padding flag for the new xforms */
   if (interp_padding) {
      for (i = old_num ; i < thiscp->num_xforms ; i++)
         thiscp->xform[i].padding=1;
   }
   
   /* If the final xform is not the last xform in the list, make it so */
   if (thiscp->final_xform_index >= 0 && thiscp->final_xform_index != thiscp->num_xforms-1) {
      tmp = thiscp->xform[thiscp->final_xform_index];
      for (i=thiscp->final_xform_index; i < thiscp->num_xforms-1; i++)
         thiscp->xform[i] = thiscp->xform[i+1];
      
      thiscp->final_xform_index = thiscp->num_xforms-1;
      thiscp->xform[thiscp->final_xform_index] = tmp;
   }
   
   if (final_flag) {
      /* Set the final xform index */
      thiscp->final_xform_enable = 1;
      thiscp->final_xform_index = thiscp->num_xforms-1;
   } else {
      /* Handle the chaos array */
      numstd = thiscp->num_xforms - (thiscp->final_xform_index>=0);
      
      /* Pad existing rows */
      for (i=0;i<oldstd;i++) {
         thiscp->chaos[i] = realloc(thiscp->chaos[i], numstd * sizeof(double));
         for (j=oldstd; j<numstd; j++)
            thiscp->chaos[i][j] = 1.0;
      }
      
      /* Add new rows */
      thiscp->chaos = realloc(thiscp->chaos,numstd * sizeof(double *));
      for (i=oldstd; i<numstd; i++) {
         thiscp->chaos[i] = malloc(numstd * sizeof(double));
         for (j=0;j<numstd;j++)
            thiscp->chaos[i][j] = 1.0;
      }
   }
}

void flam3_delete_xform(flam3_genome *thiscp, int idx_to_delete) {

   int i,j;
   int nth_std=-1;
   int num_std = thiscp->num_xforms - (thiscp->final_xform_index >= 0);

   /* Handle the final xform index */
   if (thiscp->final_xform_index != idx_to_delete) {
      /* We're going to delete the nth std xform. */

      /* Delete the nth_std row of the chaos array */
      free(thiscp->chaos[idx_to_delete]);

      /* Shift the pointers down one */
      for (i=idx_to_delete+1;i<num_std;i++)
         thiscp->chaos[i-1] = thiscp->chaos[i];
      
      /* Realloc the pointer array */
      thiscp->chaos = realloc(thiscp->chaos,(num_std-1)*sizeof(double *));
      num_std--;
      
      /* Loop over all of the rows and remove the nth_std element from them */
      for (i=0;i<num_std;i++) {
         for (j=idx_to_delete+1;j<num_std+1;j++) {
            thiscp->chaos[i][j-1] = thiscp->chaos[i][j];
         }
         /* Realloc the vector to have one less element */
         thiscp->chaos[i] = realloc(thiscp->chaos[i],num_std*sizeof(double));
         
      }
   }      
      
   
   if (thiscp->final_xform_index == idx_to_delete) {
      thiscp->final_xform_index = -1;
      thiscp->final_xform_enable = 0;
   } else if (thiscp->final_xform_index > idx_to_delete) {
      thiscp->final_xform_index--;
   }

   /* Move all of the xforms down one */
   for (i=idx_to_delete; i<thiscp->num_xforms-1; i++)
      thiscp->xform[i] = thiscp->xform[i+1];

   thiscp->num_xforms--;

   /* Reduce the memory storage by one xform */
   thiscp->xform = (flam3_xform *)realloc(thiscp->xform, sizeof(flam3_xform) * thiscp->num_xforms);
   
}



/* Copy one control point to another */
void flam3_copy(flam3_genome *dest, flam3_genome *src) {
   
   int i;
   int numstd;

   /* If there are any xforms in dest before the copy, clean them up */
   clear_cp(dest, 1);

   /* Copy main contents of genome */
   memcpy(dest, src, sizeof(flam3_genome));

   /* Only the pointer to the xform was copied, not the actual xforms. */
   /* We need to create new xform memory storage for this new cp       */
   /* This goes for chaos, too.                                        */
   dest->num_xforms = 0;
   dest->final_xform_index = -1;
   dest->xform = NULL;
   dest->chaos = NULL;

   /* Add the standard xforms first */
   numstd = src->num_xforms-(src->final_xform_index>=0);
   flam3_add_xforms(dest, numstd, 0, 0);
   for (i=0;i<numstd;i++)
      dest->xform[i] = src->xform[i];
      
   /* Add the final x if it's present */
   if (src->final_xform_index>=0) {
      flam3_add_xforms(dest, 1, 0, 1);
      dest->xform[dest->final_xform_index] = src->xform[src->final_xform_index];
   }
   //memcpy(dest->xform, src->xform, dest->num_xforms * sizeof(flam3_xform));
   
   /* Also, only the pointer to the chaos array was copied.
    * We have to take care of that as well.                 */   
   for (i=0;i<numstd;i++)
      memcpy(dest->chaos[i],src->chaos[i], numstd * sizeof(double));
         
}

void flam3_copyx(flam3_genome *dest, flam3_genome *src, int dest_std_xforms, int dest_final_xform) {

   int i,j,numsrcstd;
   
   /* If there are any xforms in dest before the copy, clean them up */
   clear_cp(dest, 1);

   /* Copy main contents of genome */
   memcpy(dest, src, sizeof(flam3_genome));

   /* Only the pointer to the xform was copied, not the actual xforms. */
   /* We need to create new xform memory storage for this new cp       */
   /* This goes for chaos, too.                                        */
   dest->num_xforms = 0;
   dest->xform = NULL;
   dest->chaos = NULL;
   dest->final_xform_index = -1;

   /* Add the padded standard xform list */
   /* Set the pad to 1 for these */
   flam3_add_xforms(dest, dest_std_xforms, 1, 0);

   numsrcstd = src->num_xforms - (src->final_xform_index >= 0);

   for(i=0;i<numsrcstd;i++) {

      /* When we copy the old xform, the pad is set to 0 */
      dest->xform[i] = src->xform[i];      
      /* Copy the initial chaos from the src - the rest are already 1 */
      memcpy(dest->chaos[i], src->chaos[i], numsrcstd*sizeof(double));
      
   }   
   
   /* Add the final xform if necessary */
   if (dest_final_xform > 0) {
      flam3_add_xforms(dest, dest_final_xform, 1, 1);

      if (src->final_xform_enable > 0) {
         dest->xform[dest->num_xforms-1] = src->xform[src->final_xform_index];
      } else {
         /* Interpolated-against final xforms need animate & color_speed set to 1.0 */
         dest->xform[dest->num_xforms-1].animate=1.0;
         dest->xform[dest->num_xforms-1].color_speed=1.0;
      }

   } else {
      dest->final_xform_index = -1;
      dest->final_xform_enable = 0;
   }

}


static flam3_genome xml_current_cp;
static flam3_genome *xml_all_cp;
/*static flam3_image_store *im_store;*/
static int xml_all_ncps;
static int flam3_conversion_failed;
/*static int xml_num_images;*/

static int flam3_atoi(char *nstr) {

    /* Note that this is NOT thread-safe, but simplifies things significantly. */
    int res;
    char *endp;

    /* Reset errno */
    errno=0;

    /* Convert the string using strtol */
    res = strtol(nstr, &endp, 10);

    /* Check errno & return string */
    if (endp!=nstr+strlen(nstr))
       flam3_conversion_failed = 1;
    if (errno)
       flam3_conversion_failed = 1;

    return(res);
}

static double flam3_atof(char *nstr) {

    /* Note that this is NOT thread-safe, but simplifies things significantly. */
    double res;
    char *endp;

    /* Reset errno */
    errno=0;

    /* Convert the string using strtod */
    res = strtod(nstr, &endp);
    
    /* Check errno & return string */
    if (endp!=nstr+strlen(nstr))
       flam3_conversion_failed = 1;
    if (errno)
       flam3_conversion_failed = 1;

    return(res);
}

void clear_cp(flam3_genome *cp, int default_flag) {
    cp->palette_index = flam3_palette_random;
    cp->center[0] = 0.0;
    cp->center[1] = 0.0;
    cp->rot_center[0] = 0.0;
    cp->rot_center[1] = 0.0;
    cp->gamma = 4.0;
    cp->vibrancy = 1.0;
    cp->contrast = 1.0;
    cp->brightness = 4.0;
    cp->symmetry = 0;
    cp->hue_rotation = 0.0;
    cp->rotate = 0.0;
    cp->edits = NULL;
    cp->pixels_per_unit = 50;
    cp->interpolation = flam3_interpolation_linear;
    cp->palette_interpolation = flam3_palette_interpolation_hsv;

    cp->genome_index = 0;
    memset(cp->parent_fname,0,flam3_parent_fn_len);

    if (default_flag==flam3_defaults_on) {
       /* If defaults are on, set to reasonable values */
       cp->highlight_power = -1.0;
       cp->background[0] = 0.0;
       cp->background[1] = 0.0;
       cp->background[2] = 0.0;
       cp->width = 100;
       cp->height = 100;
       cp->spatial_oversample = 1;
       cp->spatial_filter_radius = 0.5;
       cp->zoom = 0.0;
       cp->sample_density = 1;
       /* Density estimation stuff defaulting to ON */
       cp->estimator = 9.0;
       cp->estimator_minimum = 0.0;
       cp->estimator_curve = 0.4;
       cp->gam_lin_thresh = 0.01;
//       cp->motion_exp = 0.0;
       cp->nbatches = 1;
       cp->ntemporal_samples = 1000;
       cp->spatial_filter_select = flam3_gaussian_kernel;
       cp->interpolation_type = flam3_inttype_log;
       cp->temporal_filter_type = flam3_temporal_box;
       cp->temporal_filter_width = 1.0;
       cp->temporal_filter_exp = 0.0;
       cp->palette_mode = flam3_palette_mode_step;

    } else {
       /* Defaults are off, so set to UN-reasonable values. */
       cp->highlight_power = -1.0;
       cp->background[0] = -1.0;
       cp->background[1] = -1.0;
       cp->background[2] = -1.0;
       cp->zoom = 999999999;
       cp->spatial_oversample = -1;
       cp->spatial_filter_radius = -1;
       cp->nbatches = -1;
       cp->ntemporal_samples = -1;
       cp->width = -1;
       cp->height = -1;
       cp->sample_density = -1;
       cp->estimator = -1;
       cp->estimator_minimum = -1;
       cp->estimator_curve = -1;
       cp->gam_lin_thresh = -1;
//       cp->motion_exp = -999;
       cp->nbatches = 0;
       cp->ntemporal_samples = 0;
       cp->spatial_filter_select = -1;
       cp->interpolation_type = -1;
       cp->temporal_filter_type = -1;
       cp->temporal_filter_width = -1;
       cp->temporal_filter_exp = -999;
       cp->palette_mode = -1;
    }

    if (cp->xform != NULL && cp->num_xforms > 0) {
        int i;
        int ns = cp->num_xforms - (cp->final_xform_index>=0);
        
       for (i=0;i<ns;i++) {
          free(cp->chaos[i]);
       }
       free(cp->chaos);
       cp->chaos=NULL;
              
       free(cp->xform);
       cp->xform=NULL;
       
       cp->num_xforms = 0;
    }

    cp->final_xform_enable = 0;
    cp->final_xform_index = -1;

}

static double flam3_spatial_support[flam3_num_spatialfilters] = {

   1.5, /* gaussian */
   1.0, /* hermite */
   0.5, /* box */
   1.0, /* triangle */
   1.5, /* bell */
   2.0, /* b spline */
   2.0, /* mitchell */
   1.0, /* blackman */
   2.0, /* catrom */
   1.0, /* hanning */
   1.0, /* hamming */
   3.0, /* lanczos3 */
   2.0, /* lanczos2 */
   1.5  /* quadratic */
};

double flam3_spatial_filter(int knum, double x) {

   if (knum==0)
      return flam3_gaussian_filter(x);
   else if (knum==1)
      return flam3_hermite_filter(x);
   else if (knum==2)
      return flam3_box_filter(x);
   else if (knum==3)
      return flam3_triangle_filter(x);
   else if (knum==4)
      return flam3_bell_filter(x);
   else if (knum==5)
      return flam3_b_spline_filter(x);
   else if (knum==6)
      return flam3_mitchell_filter(x);
   else if (knum==7)
      return flam3_sinc(x)*flam3_blackman_filter(x);
   else if (knum==8)
      return flam3_catrom_filter(x);
   else if (knum==9)
      return flam3_sinc(x)*flam3_hanning_filter(x);
   else if (knum==10)
      return flam3_sinc(x)*flam3_hamming_filter(x);
   else if (knum==11)
      return flam3_lanczos3_filter(x)*flam3_sinc(x/3.0);   
   else if (knum==12)
      return flam3_lanczos2_filter(x)*flam3_sinc(x/2.0);
   else if (knum==13)
      return flam3_quadratic_filter(x);
   else {
      fprintf(stderr,"Unknown filter kernel %d!\n",knum);
      exit(1);
   }
}
      
char *flam3_variation_names[1+flam3_nvariations] = {
  "linear",
  "sinusoidal",
  "spherical",
  "swirl",
  "horseshoe",
  "polar",
  "handkerchief",
  "heart",
  "disc",
  "spiral",
  "hyperbolic",
  "diamond",
  "ex",
  "julia",
  "bent",
  "waves",
  "fisheye",
  "popcorn",
  "exponential",
  "power",
  "cosine",
  "rings",
  "fan",
  "blob",
  "pdj",
  "fan2",
  "rings2",
  "eyefish",
  "bubble",
  "cylinder",
  "perspective",
  "noise",
  "julian",
  "juliascope",
  "blur",
  "gaussian_blur",
  "radial_blur",
  "pie",
  "ngon",
  "curl",
  "rectangles",
  "arch",
  "tangent",
  "square",
  "rays",
  "blade",
  "secant2",
  "twintrian",
  "cross",
  "disc2",
  "super_shape",
  "flower",
  "conic",
  "parabola",
  "bent2",
  "bipolar",
  "boarders",
  "butterfly",
  "cell",
  "cpow",
  "curve",
  "edisc",
  "elliptic",
  "escher",
  "foci",
  "lazysusan",
  "loonie",
  "pre_blur",
  "modulus",
  "oscilloscope",
  "polar2",
  "popcorn2",
  "scry",
  "separation",
  "split",
  "splits",
  "stripes",
  "wedge",
  "wedge_julia",
  "wedge_sph",
  "whorl",
  "waves2",
  0
};


static int var2n(const char *s) {
  int i;
  for (i = 0; i < flam3_nvariations; i++)
    if (!strcmp(s, flam3_variation_names[i])) return i;
  return flam3_variation_none;
}

static int parse_flame_element(xmlNode *flame_node) {
   flam3_genome *cp = &xml_current_cp;
   xmlNode *chld_node, *motion_node;
   xmlNodePtr edit_node;
   xmlAttrPtr att_ptr, cur_att;
   int solo_xform=-1;
   char *att_str;
   int num_std_xforms=-1;
   char *cpy;
   char tmps[2];
   int i;
   int ix = 0;
   flam3_xform tmpcpy;
   flam3_chaos_entry *xaos=NULL;
   int num_xaos=0;

   /* Reset the conversion error flag */
   /* NOT threadsafe                  */
   flam3_conversion_failed=0;

   /* Store this flame element in the current cp */

   /* The top level element is a flame element. */
   /* Read the attributes of it and store them. */
   att_ptr = flame_node->properties;

   if (att_ptr==NULL) {
      fprintf(stderr, "Error : <flame> element has no attributes.\n");
      return(1);
   }

   memset(cp->flame_name,0,flam3_name_len+1);

   for (cur_att = att_ptr; cur_att; cur_att = cur_att->next) {

       att_str = (char *) xmlGetProp(flame_node,cur_att->name);
       
      /* Compare attribute names */
      if (!xmlStrcmp(cur_att->name, (const xmlChar *)"time")) {
         cp->time = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"interpolation")) {
     if (!strcmp("linear", att_str)) {
         cp->interpolation = flam3_interpolation_linear;
     } else if  (!strcmp("smooth", att_str)) {
         cp->interpolation = flam3_interpolation_smooth;
     } else {
         fprintf(stderr, "warning: unrecognized interpolation type %s.\n", att_str);
     }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"palette_interpolation")) {
     if (!strcmp("hsv", att_str)) {
         cp->palette_interpolation = flam3_palette_interpolation_hsv;
     } else if  (!strcmp("sweep", att_str)) {
         cp->palette_interpolation = flam3_palette_interpolation_sweep;
     } else {
         fprintf(stderr, "warning: unrecognized palette interpolation type %s.\n", att_str);
     }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"interpolation_space") ||
                 !xmlStrcmp(cur_att->name, (const xmlChar *)"interpolation_type")) {
      
         if (!strcmp("linear", att_str))
            cp->interpolation_type = flam3_inttype_linear;
         else if (!strcmp("log", att_str))
            cp->interpolation_type = flam3_inttype_log;
         else if (!strcmp("old", att_str))
            cp->interpolation_type = flam3_inttype_compat;
         else if (!strcmp("older", att_str))
            cp->interpolation_type = flam3_inttype_older;
         else
            fprintf(stderr,"warning: unrecognized interpolation_type %s.\n",att_str);
     
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"name")) {
         strncpy(cp->flame_name, att_str, flam3_name_len);
         i = (int)strlen(cp->flame_name)-1;
         while(i-->0) {
            if (isspace(cp->flame_name[i]))
               cp->flame_name[i] = '_';
         }

      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"palette")) {
         cp->palette_index = flam3_atoi(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"size")) {
         if (sscanf(att_str, "%d %d%1s", &cp->width, &cp->height, tmps) != 2) {
            fprintf(stderr,"error: invalid size attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"center")) {
         if (sscanf(att_str, "%lf %lf%1s", &cp->center[0], &cp->center[1], tmps) != 2) {
            fprintf(stderr,"error: invalid center attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
         cp->rot_center[0] = cp->center[0];
         cp->rot_center[1] = cp->center[1];
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"scale")) {
         cp->pixels_per_unit = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rotate")) {
         cp->rotate = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"zoom")) {
         cp->zoom = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oversample")) {
         cp->spatial_oversample = flam3_atoi(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"supersample")) {
         cp->spatial_oversample = flam3_atoi(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"filter")) {
         cp->spatial_filter_radius = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"filter_shape")) {
         if (!strcmp("gaussian", att_str))
            cp->spatial_filter_select = flam3_gaussian_kernel;
         else if (!strcmp("hermite", att_str))
            cp->spatial_filter_select = flam3_hermite_kernel;
         else if (!strcmp("box", att_str))
            cp->spatial_filter_select = flam3_box_kernel;
         else if (!strcmp("triangle", att_str))
            cp->spatial_filter_select = flam3_triangle_kernel;
         else if (!strcmp("bell", att_str))
            cp->spatial_filter_select = flam3_bell_kernel;
         else if (!strcmp("bspline", att_str))
            cp->spatial_filter_select = flam3_b_spline_kernel;
         else if (!strcmp("mitchell", att_str))
            cp->spatial_filter_select = flam3_mitchell_kernel;
         else if (!strcmp("blackman", att_str))
            cp->spatial_filter_select = flam3_blackman_kernel;
         else if (!strcmp("catrom", att_str))
            cp->spatial_filter_select = flam3_catrom_kernel;
         else if (!strcmp("hanning", att_str))
            cp->spatial_filter_select = flam3_hanning_kernel;
         else if (!strcmp("hamming", att_str))
            cp->spatial_filter_select = flam3_hamming_kernel;
         else if (!strcmp("lanczos3", att_str))
            cp->spatial_filter_select = flam3_lanczos3_kernel;
         else if (!strcmp("lanczos2", att_str))
            cp->spatial_filter_select = flam3_lanczos2_kernel;
         else if (!strcmp("quadratic", att_str))
            cp->spatial_filter_select = flam3_quadratic_kernel;
         else
            fprintf(stderr, "warning: unrecognized kernel shape %s.  Using gaussian.\n", att_str);

      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"temporal_filter_type")) {
         if (!strcmp("box", att_str))
            cp->temporal_filter_type = flam3_temporal_box;
         else if (!strcmp("gaussian", att_str))
            cp->temporal_filter_type = flam3_temporal_gaussian;
         else if (!strcmp("exp",att_str))
            cp->temporal_filter_type = flam3_temporal_exp;
         else
            fprintf(stderr, "warning: unrecognized temporal filter %s.  Using box.\n",att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"temporal_filter_width")) {
         cp->temporal_filter_width = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"temporal_filter_exp")) {
         cp->temporal_filter_exp = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"palette_mode")) {
         if (!strcmp("step", att_str))
            cp->palette_mode = flam3_palette_mode_step;
         else if (!strcmp("linear", att_str))
            cp->palette_mode = flam3_palette_mode_linear;
         else
            fprintf(stderr,"warning: unrecognized palette mode %s.  Using step.\n",att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"quality")) {
         cp->sample_density = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"passes")) {
         cp->nbatches = flam3_atoi(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"temporal_samples")) {
         cp->ntemporal_samples = flam3_atoi(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"background")) {
         if (sscanf(att_str, "%lf %lf %lf%1s", &cp->background[0], &cp->background[1], &cp->background[2], tmps) != 3) {
            fprintf(stderr,"error: invalid background attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"brightness")) {
         cp->brightness = flam3_atof(att_str);
/*      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"contrast")) {
         cp->contrast = flam3_atof(att_str);*/
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"gamma")) {
         cp->gamma = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"highlight_power")) {
         cp->highlight_power = atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"vibrancy")) {
         cp->vibrancy = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"hue")) {
         cp->hue_rotation = fmod(flam3_atof(att_str), 1.0);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"estimator_radius")) {
         cp->estimator = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"estimator_minimum")) {
         cp->estimator_minimum = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"estimator_curve")) {
         cp->estimator_curve = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"gamma_threshold")) {
         cp->gam_lin_thresh = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"soloxform")) {
         solo_xform = flam3_atof(att_str);
      }

      xmlFree(att_str);

   }

   /* Finished with flame attributes.  Now look at children of flame element. */
   for (chld_node=flame_node->children; chld_node; chld_node = chld_node->next) {

      /* Is this a color node? */
      if (!xmlStrcmp(chld_node->name, (const xmlChar *)"color")) {
         double index = -1;
         double r=0.0,g=0.0,b=0.0;

         /* Loop through the attributes of the color element */
         att_ptr = chld_node->properties;

         if (att_ptr==NULL) {
            fprintf(stderr,"Error:  No attributes for color element.\n");
            return(1);
         }

         for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

            att_str = (char *) xmlGetProp(chld_node,cur_att->name);

            if (!xmlStrcmp(cur_att->name, (const xmlChar *)"index")) {
               index = flam3_atof(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rgb")) {
               if (sscanf(att_str, "%lf %lf %lf%1s", &r, &g, &b, tmps) != 3) {
                  fprintf(stderr,"error: invalid rgb attribute '%s'\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
            } else {
               fprintf(stderr,"Error:  Unknown color attribute '%s'\n",cur_att->name);
               xmlFree(att_str);
               return(1);
            }

            xmlFree(att_str);
         }

         if (index >= 0 && index <= 255.0) {
            cp->palette[ix].color[0] = r / 255.0;
            cp->palette[ix].color[1] = g / 255.0;
            cp->palette[ix].color[2] = b / 255.0;
            cp->palette[ix].color[3] = 1.0;
            cp->palette[ix].index = index;
            ix++;
         } else {
            fprintf(stderr,"Error:  Color element with bad/missing index attribute (%d)\n",index);
            return(1);
         }
      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"colors")) {

         int count;

         /* Loop through the attributes of the colors element */
         att_ptr = chld_node->properties;

         if (att_ptr==NULL) {
            fprintf(stderr,"Error: No attributes for colors element.\n");
            return(1);
         }

         for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

            att_str = (char *) xmlGetProp(chld_node,cur_att->name);

            if (!xmlStrcmp(cur_att->name, (const xmlChar *)"count")) {
               count = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"data")) {
               if (flam3_parse_hexformat_colors(att_str, cp, count, 4) > 0) {
                  fprintf(stderr,"error parsing hexformatted colors\n");
                  xmlFree(att_str);
                  return(1);
               }
            } else {
               fprintf(stderr,"Error:  Unknown color attribute '%s'\n",cur_att->name);
               xmlFree(att_str);
               return(1);
            }

            xmlFree(att_str);
         }


      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"palette")) {

         /* This could be either the old form of palette or the new form */
         /* Make sure BOTH are not specified, otherwise either are ok    */
         int numcolors=0;
         int numbytes=0;
         int old_format=0;
         int new_format=0;
         int index0, index1;
         double hue0, hue1;
         double blend = 0.5;
         index0 = index1 = flam3_palette_random;
         hue0 = hue1 = 0.0;

         /* Loop through the attributes of the palette element */
         att_ptr = chld_node->properties;

         if (att_ptr==NULL) {
            fprintf(stderr,"Error:  No attributes for palette element.\n");
            return(1);
         }

         for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

            att_str = (char *) xmlGetProp(chld_node,cur_att->name);

            if (!xmlStrcmp(cur_att->name, (const xmlChar *)"index0")) {
               old_format++;
               index0 = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"index1")) {
               old_format++;
               index1 = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"hue0")) {
               old_format++;
               hue0 = flam3_atof(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"hue1")) {
               old_format++;
               hue1 = flam3_atof(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"blend")) {
               old_format++;
               blend = flam3_atof(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"count")) {
               new_format++;
               numcolors = flam3_atoi(att_str);
            } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"format")) {
               new_format++;
               if (!strcmp(att_str,"RGB"))
                  numbytes=3;
               else if (!strcmp(att_str,"RGBA"))
                  numbytes=4;
               else {
                  fprintf(stderr,"Error: Unrecognized palette format string (%s)\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
            } else {
               fprintf(stderr,"Error:  Unknown palette attribute '%s'\n",cur_att->name);
               xmlFree(att_str);
               return(1);
            }

            xmlFree(att_str);
         }

         /* Old or new format? */
         if (new_format>0 && old_format>0) {
            fprintf(stderr,"Error: mixing of old and new palette tag syntax not allowed.\n");
            return(1);
         }

         if (old_format>0)
            interpolate_cmap(cp->palette, blend, index0, hue0, index1, hue1);
         else {

            char *pal_str;

            /* Read formatted string from contents of tag */

            pal_str = (char *) xmlNodeGetContent(chld_node);
            /* !!! check hexform parse for errors */

            if (flam3_parse_hexformat_colors(pal_str, cp, numcolors, numbytes) > 0) {
               fprintf(stderr,"error reading hexformatted colors\n");
               xmlFree(pal_str);
               return(1);
            }

            xmlFree(pal_str);
         }
      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"symmetry")) {

         int kind=0;
         int bef,aft;

         /* Loop through the attributes of the symmetry element */
         att_ptr = chld_node->properties;

         if (att_ptr==NULL) {
            fprintf(stderr,"Error:  No attributes for symmetry element.\n");
            exit(1);
         }

         for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

            att_str = (char *) xmlGetProp(chld_node,cur_att->name);

            if (!xmlStrcmp(cur_att->name, (const xmlChar *)"kind")) {
               kind = flam3_atoi(att_str);
            } else {
               fprintf(stderr,"Error:  Unknown symmetry attribute '%s'\n",cur_att->name);
               xmlFree(att_str);
               return(1);
            }

            xmlFree(att_str);
         }

         bef = cp->num_xforms;
         flam3_add_symmetry(cp,kind);
         aft = cp->num_xforms;
         num_std_xforms += (aft-bef);
         

      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"xform") ||
                  !xmlStrcmp(chld_node->name, (const xmlChar *)"finalxform")) {

         int xf = cp->num_xforms;
         
         if (!xmlStrcmp(chld_node->name, (const xmlChar *)"finalxform")) {

            if (cp->final_xform_index >=0) {
               fprintf(stderr,"Error:  Cannot specify more than one final xform.\n");
               return(1);
            }

            flam3_add_xforms(cp, 1, 0, 1);
            cp->xform[xf].var[0]=0.0;
            cp->final_xform_index = xf;
            /* Now, if present, the xform enable defaults to on */
            cp->final_xform_enable = 1;
            
         } else {

            /* Add one to the counter */
            flam3_add_xforms(cp, 1, 0, 0);
            cp->xform[xf].var[0]=0.0;
            num_std_xforms++;
                     
         }
         
         if (parse_xform_xml(chld_node, &(cp->xform[xf]), &num_xaos, xaos, num_std_xforms, 0) != 0)
            return(1);
         
         if (cp->final_xform_index == xf && cp->xform[xf].density != 0.0) {
            fprintf(stderr,"Error: Final xforms should not have weight specified.\n");
            return(1);
         }
         
         /* Check for non-zero motion_* params */
         if (cp->xform[xf].motion_freq != 0 || cp->xform[xf].motion_func != 0) {
            fprintf(stderr,"Error: Motion parameters should not be specified in xforms.\n");
            return(1);
         }
         
         
         /* Motion Language:  Check the xform element for children - should be named 'motion'. */
         for (motion_node=chld_node->children; motion_node; motion_node = motion_node->next) {
                  
            if (!xmlStrcmp(motion_node->name, (const xmlChar *)"motion")) {
            
               int nm = cp->xform[xf].num_motion;

               /* Add motion element to xform */
               flam3_add_motion_element( &cp->xform[xf] );            
               
               /* Read motion xml */
               if (parse_xform_xml(motion_node, &(cp->xform[xf].motion[nm]), NULL, NULL, 0, 1) != 0)
                  return(1);
               
            }
            
         }
            

      } else if (!xmlStrcmp(chld_node->name, (const xmlChar *)"edit")) {

         /* Create a new XML document with this edit node as the root node */
         cp->edits = xmlNewDoc( (const xmlChar *)"1.0");
         edit_node = xmlCopyNode( chld_node, 1 );
         xmlDocSetRootElement(cp->edits, edit_node);

      }
   } /* Done parsing flame element. */
   
   num_std_xforms++;

   for (i=0;i<num_std_xforms;i++) {
      
      /* Adjust visibility with solo xform setting */
      if (solo_xform>=0 && i!=solo_xform)
         cp->xform[i].visibility = 0.0;

   }
   
   /* Set the chaos array entries with the values in the xaos list */
   for (i=0;i<num_xaos;i++)
      cp->chaos[xaos[i].from][xaos[i].to] = xaos[i].scalar;
      
   free(xaos);
   
   /* If there is a final xform in this cp, move it to the end of the list */   
   if (cp->final_xform_index >=0 && cp->final_xform_index != (cp->num_xforms-1)) {
      /* Make a copy of the final xform */
      tmpcpy = cp->xform[cp->final_xform_index];
      
      /* Move each other xform up one */
      for (i=cp->final_xform_index+1;i<cp->num_xforms;i++)
         cp->xform[i-1] = cp->xform[i];
         
      /* Put the final at the end */
      cp->xform[cp->num_xforms-1] = tmpcpy;
      
      cp->final_xform_index = cp->num_xforms - 1;
   }
      
   //print_chaos(cp);
   
   /* Check for bad parse */
   if (flam3_conversion_failed) {
      fprintf(stderr,"error: parsing a double or int attribute's value.\n");
      return(1);
   }
   
   return(0);

}

static int parse_xform_xml(xmlNode *chld_node,flam3_xform *this_xform, int *num_xaos, flam3_chaos_entry *xaos, int numstd, int motionxf) {

   xmlAttrPtr att_ptr, cur_att;
   char *att_str, *cpy;
   char tmps[2];
   int j,k;

   /* Loop through the attributes of the xform element */
   att_ptr = chld_node->properties;

   if (att_ptr==NULL) {
      fprintf(stderr,"Error: No attributes for element.\n");
      return(1);
   }

   for (cur_att=att_ptr; cur_att; cur_att = cur_att->next) {

      att_str = (char *) xmlGetProp(chld_node,cur_att->name);

      cpy = att_str;
      if (!xmlStrcmp(cur_att->name, (const xmlChar *)"weight")) {
         this_xform->density = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"symmetry")) {
         /* Deprecated.  Set both color_speed and animate to this value. */
         this_xform->color_speed = flam3_atof(att_str);
         this_xform->animate = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"color_speed")) {
         this_xform->color_speed = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"animate")) {
         this_xform->animate = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"motion_freq")) {
         this_xform->motion_freq = flam3_atoi(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"motion_func")) {
         if (!strcmp("sin", att_str)) {
            this_xform->motion_func = MOTION_SIN;
         } else if (!strcmp("triangle",att_str)) {
            this_xform->motion_func = MOTION_TRIANGLE;
         } else if (!strcmp("cos",att_str)) {
            this_xform->motion_func = MOTION_COS;
         } else {
            fprintf(stderr,"Error: unknown motion function '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }

      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"color")) {
         double tmpc1;
         this_xform->color = 0.0;
         /* Try two coords first */
         if (sscanf(att_str, "%lf %lf%1s", &this_xform->color, &tmpc1, tmps) != 2) {
            /* Try one color */
            if (sscanf(att_str, "%lf%1s", &this_xform->color,tmps) != 1) {
               fprintf(stderr,"Error: malformed xform color attribute '%s'\n",att_str);
               xmlFree(att_str);
               return(1);
            }
         }            
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"var1")) {
         for (j=0; j < flam3_nvariations; j++) {
            this_xform->var[j] = 0.0;
         }
         j = flam3_atoi(att_str);

         if (j < 0 || j >= flam3_nvariations) {
            fprintf(stderr,"Error:  Bad variation (%d)\n",j);
            xmlFree(att_str);
            return(1);
         }

         this_xform->var[j] = 1.0;
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"var")) {
         for (j=0; j < flam3_nvariations; j++) {
            char *cpy2;
            errno=0;
            this_xform->var[j] = strtod(cpy, &cpy2);
            if (errno != 0 || cpy==cpy2) {
               fprintf(stderr,"error: bad value in var attribute '%s'\n",att_str);
               xmlFree(att_str);
               return(1);
            }
            cpy=cpy2;
         }

         if (cpy != att_str+strlen(att_str)) {
            fprintf(stderr,"error: extra chars at the end of var attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"chaos")) {      
         /* Chaos scalars */
         
         char *tok;
         double scal;
         int toi=0;
         
         if (motionxf==1) {
            fprintf(stderr,"error: motion element cannot have a chaos attribute.\n");
            xmlFree(att_str);
            return(1);
         }
         
         /* The att string contains at least one value, delimited by a space */
         tok = strtok(cpy," ");
         while (tok!=NULL) {
            scal = flam3_atof(tok);
            
            /* Skip 1.0 entries */
            if (scal==1.0) {
               toi++;
               tok = strtok(NULL," ");
               continue;
            }
            
            /* Realloc the xaos list */
            xaos = realloc(xaos,(*num_xaos+1) * sizeof(flam3_chaos_entry));

            /* Populate the xaos list */
            xaos[*num_xaos].from = numstd;
            xaos[*num_xaos].to = toi;
            xaos[*num_xaos].scalar = scal;
            toi++;
            *num_xaos++;
                              
            /* Get the next token */
            tok = strtok(NULL," ");
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"plotmode")) {

         if (motionxf==1) {
            fprintf(stderr,"error: motion element cannot have a plotmode attribute.\n");
            xmlFree(att_str);
            return(1);
         }
         
         if (!strcmp("off", att_str))
            this_xform->visibility = 0.0;
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"visibility")) {
         this_xform->visibility = flam3_atof(att_str);
         
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"coefs")) {
         for (k=0; k<3; k++) {
            for (j=0; j<2; j++) {
               char *cpy2;
               errno = 0;
               this_xform->c[k][j] = strtod(cpy, &cpy2);
               if (errno != 0 || cpy==cpy2) {
                  fprintf(stderr,"error: bad value in coefs attribute '%s'\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
               cpy=cpy2;
            }
         }
         if (cpy != att_str+strlen(att_str)) {
            fprintf(stderr,"error: extra chars at the end of coefs attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"post")) {
         for (k = 0; k < 3; k++) {
            for (j = 0; j < 2; j++) {
               char *cpy2;
               errno = 0;
               this_xform->post[k][j] = strtod(cpy, &cpy2);
               if (errno != 0 || cpy==cpy2) {
                  fprintf(stderr,"error: bad value in post attribute '%s'\n",att_str);
                  xmlFree(att_str);
                  return(1);
               }
               cpy=cpy2;
            }
         }
         if (cpy != att_str+strlen(att_str)) {
            fprintf(stderr,"error: extra chars at end of post attribute '%s'\n",att_str);
            xmlFree(att_str);
            return(1);
         }
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"blob_low")) {
         this_xform->blob_low = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"blob_high")) {
         this_xform->blob_high = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"blob_waves")) {
         this_xform->blob_waves = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pdj_a")) {
         this_xform->pdj_a = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pdj_b")) {
         this_xform->pdj_b = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pdj_c")) {
         this_xform->pdj_c = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pdj_d")) {
         this_xform->pdj_d = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"fan2_x")) {
         this_xform->fan2_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"fan2_y")) {
         this_xform->fan2_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rings2_val")) {
         this_xform->rings2_val = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"perspective_angle")) {
         this_xform->perspective_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"perspective_dist")) {
         this_xform->perspective_dist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"julian_power")) {
         this_xform->juliaN_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"julian_dist")) {
         this_xform->juliaN_dist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"juliascope_power")) {
         this_xform->juliaScope_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"juliascope_dist")) {
         this_xform->juliaScope_dist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"radial_blur_angle")) {
         this_xform->radialBlur_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pie_slices")) {
         this_xform->pie_slices = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pie_rotation")) {
         this_xform->pie_rotation = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"pie_thickness")) {
         this_xform->pie_thickness = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"ngon_sides")) {
         this_xform->ngon_sides = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"ngon_power")) {
         this_xform->ngon_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"ngon_circle")) {
         this_xform->ngon_circle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"ngon_corners")) {
         this_xform->ngon_corners = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curl_c1")) {
         this_xform->curl_c1 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curl_c2")) {
         this_xform->curl_c2 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rectangles_x")) {
         this_xform->rectangles_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"rectangles_y")) {
         this_xform->rectangles_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"amw_amp")) {
         this_xform->amw_amp = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"disc2_rot")) {
         this_xform->disc2_rot = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"disc2_twist")) {
         this_xform->disc2_twist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_rnd")) {
         this_xform->supershape_rnd = flam3_atof(att_str);
         /* Limit to [0,1] */
         if (this_xform->supershape_rnd<0)
            this_xform->supershape_rnd=0;
         else if (this_xform->supershape_rnd>1)
            this_xform->supershape_rnd=1;
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_m")) {
         this_xform->supershape_m = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_n1")) {
         this_xform->supershape_n1 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_n2")) {
         this_xform->supershape_n2 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_n3")) {
         this_xform->supershape_n3 = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"super_shape_holes")) {
         this_xform->supershape_holes = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"flower_petals")) {
         this_xform->flower_petals = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"flower_holes")) {
         this_xform->flower_holes = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"conic_eccentricity")) {
         this_xform->conic_eccen = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"conic_holes")) {
         this_xform->conic_holes = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"parabola_height")) {
         this_xform->parabola_height = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"parabola_width")) {
         this_xform->parabola_width = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bent2_x")) {
         this_xform->bent2_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bent2_y")) {
         this_xform->bent2_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"bipolar_shift")) {
         this_xform->bipolar_shift = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"cell_size")) {
         this_xform->cell_size = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"cpow_i")) {
         this_xform->cpow_i = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"cpow_r")) {
         this_xform->cpow_r = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"cpow_power")) {
         this_xform->cpow_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curve_xamp")) {
         this_xform->curve_xamp = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curve_yamp")) {
         this_xform->curve_yamp = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curve_xlength")) {
         this_xform->curve_xlength = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"curve_ylength")) {
         this_xform->curve_ylength = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"escher_beta")) {
         this_xform->escher_beta = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_x")) {
         this_xform->lazysusan_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_y")) {
         this_xform->lazysusan_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_spin")) {
         this_xform->lazysusan_spin = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_space")) {
         this_xform->lazysusan_space = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"lazysusan_twist")) {
         this_xform->lazysusan_twist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"modulus_x")) {
         this_xform->modulus_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"modulus_y")) {
         this_xform->modulus_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscope_separation")) {
         this_xform->oscope_separation = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscope_frequency")) {
         this_xform->oscope_frequency = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscope_amplitude")) {
         this_xform->oscope_amplitude = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"oscope_damping")) {
         this_xform->oscope_damping = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"popcorn2_x")) {
         this_xform->popcorn2_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"popcorn2_y")) {
         this_xform->popcorn2_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"popcorn2_c")) {
         this_xform->popcorn2_c = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"separation_x")) {
         this_xform->separation_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"separation_xinside")) {
         this_xform->separation_xinside = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"separation_y")) {
         this_xform->separation_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"separation_yinside")) {
         this_xform->separation_yinside = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"split_xsize")) {
         this_xform->split_xsize = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"split_ysize")) {
         this_xform->split_ysize = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"splits_x")) {
         this_xform->splits_x = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"splits_y")) {
         this_xform->splits_y = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"stripes_space")) {
         this_xform->stripes_space = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"stripes_warp")) {
         this_xform->stripes_warp = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_angle")) {
         this_xform->wedge_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_hole")) {
         this_xform->wedge_hole = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_count")) {
         this_xform->wedge_count = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_swirl")) {
         this_xform->wedge_swirl = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_julia_angle")) {
         this_xform->wedge_julia_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_julia_count")) {
         this_xform->wedge_julia_count = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_julia_power")) {
         this_xform->wedge_julia_power = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_julia_dist")) {
         this_xform->wedge_julia_dist = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_sph_angle")) {
         this_xform->wedge_sph_angle = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_sph_hole")) {
         this_xform->wedge_sph_hole = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_sph_count")) {
         this_xform->wedge_sph_count = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"wedge_sph_swirl")) {
         this_xform->wedge_sph_swirl = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"whorl_inside")) {
         this_xform->whorl_inside = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"whorl_outside")) {
         this_xform->whorl_outside = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"waves2_scalex")) {
         this_xform->waves2_scalex = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"waves2_scaley")) {
         this_xform->waves2_scaley = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"waves2_freqx")) {
         this_xform->waves2_freqx = flam3_atof(att_str);
      } else if (!xmlStrcmp(cur_att->name, (const xmlChar *)"waves2_freqy")) {
         this_xform->waves2_freqy = flam3_atof(att_str);
      } else {
         int v = var2n((char *) cur_att->name);
         if (v != flam3_variation_none)
            this_xform->var[v] = flam3_atof(att_str);
         else
            fprintf(stderr,"Warning: unrecognized variation %s.  Ignoring.\n",(char *)cur_att->name);
      }


      xmlFree(att_str);
   }
   return(0);
}

int flam3_count_nthreads(void) {
   int nthreads;
#ifndef HAVE_LIBPTHREAD
   return(1);
#endif

#ifdef WIN32
   SYSTEM_INFO sysInfo;
   GetSystemInfo(&sysInfo);
   nthreads = sysInfo.dwNumberOfProcessors;
#else
#ifdef __APPLE__
   kern_return_t    kr;
   host_name_port_t   host;
   unsigned int     size;
   struct host_basic_info   hi;

   host = mach_host_self();
   size = sizeof(hi)/sizeof(int);
   kr = host_info(host, HOST_BASIC_INFO, (host_info_t)&hi, &size);
   if (kr != KERN_SUCCESS) {
       mach_error("host_info():", kr);
       exit(EXIT_FAILURE);
   }
   nthreads = hi.avail_cpus;
#else 
#ifndef _SC_NPROCESSORS_ONLN
   char line[MAXBUF];
   FILE *f = fopen("/proc/cpuinfo", "r");
   if (NULL == f) goto def;
   nthreads = 0;
   while (fgets(line, MAXBUF, f)) {
      if (!strncmp("processor\t:", line, 11))
         nthreads++;
   }
   fclose(f);
   if (nthreads < 1) goto def;
   return (nthreads);
def:
   fprintf(stderr, "could not read /proc/cpuinfo, using one render thread.\n");
   nthreads = 1;
#else
   nthreads = sysconf(_SC_NPROCESSORS_ONLN);
   if (nthreads < 1) nthreads = 1;
#endif
#endif
#endif
   return (nthreads);
}

static void scan_for_flame_nodes(xmlNode *cur_node, char *parent_file, int default_flag) {

   xmlNode *this_node = NULL;
   size_t f3_storage; //,im_storage;
   int pfe_success;

   /* Loop over this level of elements */
   for (this_node=cur_node; this_node; this_node = this_node->next) {

      /* Check to see if this element is a <flame> element */
      if (this_node->type == XML_ELEMENT_NODE && !xmlStrcmp(this_node->name, (const xmlChar *)"flame")) {

         /* This is a flame element.  Parse it. */
         clear_cp(&xml_current_cp, default_flag);

         pfe_success = parse_flame_element(this_node);
         
         if (pfe_success>0) {
            fprintf(stderr,"error parsing flame element\n");
            xml_all_cp = NULL; /* leaks memory but terminates */
            xml_all_ncps = 0;
            return;
         }

         /* Copy this cp into the array */
         f3_storage = (1+xml_all_ncps)*sizeof(flam3_genome);
         xml_all_cp = realloc(xml_all_cp, f3_storage);
         /* Clear out the realloc'd memory */
         memset(&(xml_all_cp[xml_all_ncps]),0,sizeof(flam3_genome));

         if (xml_current_cp.palette_index != flam3_palette_random) {
            flam3_get_palette(xml_current_cp.palette_index, xml_current_cp.palette,
               xml_current_cp.hue_rotation);
         }

         xml_current_cp.genome_index = xml_all_ncps;
         memset(xml_current_cp.parent_fname, 0, flam3_parent_fn_len);
         strncpy(xml_current_cp.parent_fname,parent_file,flam3_parent_fn_len-1);

         flam3_copy(&(xml_all_cp[xml_all_ncps]), &xml_current_cp);
         xml_all_ncps ++;

/*      } else if (this_node->type == XML_ELEMENT_NODE && !xmlStrcmp(this_node->name, (const xmlChar *)"image")) {

         im_storage = (1+xml_num_images)*sizeof(flam3_image_store);
         im_store = realloc(im_store, im_storage);
         parse_image_element(this_node);
         xml_num_images++;

      }*/
      } else {
         /* Check all of the children of this element */
         scan_for_flame_nodes(this_node->children, parent_file, default_flag);
      }
   }
}

flam3_genome *flam3_parse_xml2(char *xmldata, char *xmlfilename, int default_flag, int *ncps) {

   xmlDocPtr doc; /* Parsed XML document tree */
   xmlNode *rootnode;
   char *bn;
   int i;

   /* Parse XML string into internal document */
   /* Forbid network access during read       */
   doc = xmlReadMemory(xmldata, (int)strlen(xmldata), xmlfilename, NULL, XML_PARSE_NONET);

   /* Check for errors */
   if (doc==NULL) {
      fprintf(stderr, "Failed to parse %s\n", xmlfilename);
      return NULL;
   }

   /* What is the root node of the document? */
   rootnode = xmlDocGetRootElement(doc);

   /* Scan for <flame> nodes, starting with this node */
   xml_all_cp = NULL;
   xml_all_ncps = 0;
   memset(&xml_current_cp, 0, sizeof(flam3_genome));
   /*xml_num_images = 0;*/

   bn = basename(xmlfilename);
   scan_for_flame_nodes(rootnode, bn, default_flag);

   xmlFreeDoc(doc);

   *ncps = xml_all_ncps;
   
   /* Check to see if the first control point or the second-to-last */
   /* control point has interpolation="smooth".  This is invalid    */
   /* and should be reset to linear (with a warning).               */
   if (*ncps>=1) {
      if (xml_all_cp[0].interpolation == flam3_interpolation_smooth) {
         fprintf(stderr,"Warning: smooth interpolation cannot be used for first segment.\n"
                        "         switching to linear.\n");
         xml_all_cp[0].interpolation = flam3_interpolation_linear;
      }
   }
   
   if (*ncps>=2) {
      if (xml_all_cp[*ncps-2].interpolation == flam3_interpolation_smooth) {
         fprintf(stderr,"Warning: smooth interpolation cannot be used for last segment.\n"
                        "         switching to linear.\n");
         xml_all_cp[*ncps-2].interpolation = flam3_interpolation_linear;
      }
   }
   
   /* Finally, ensure that consecutive 'rotate' parameters never exceed */
   /* a difference of more than 180 degrees (+/-) for interpolation.    */
   /* An adjustment of +/- 360 degrees is made until this is true.      */
   if (*ncps>1) {
   
      for (i=1;i<*ncps;i++) {

         /* Only do this adjustment if we're not in compat mode */
         if (flam3_inttype_compat != xml_all_cp[i-1].interpolation_type
        && flam3_inttype_older != xml_all_cp[i-1].interpolation_type) {
      
            while (xml_all_cp[i].rotate < xml_all_cp[i-1].rotate-180)
               xml_all_cp[i].rotate += 360;
            
            while (xml_all_cp[i].rotate > xml_all_cp[i-1].rotate+180)
               xml_all_cp[i].rotate -= 360;
         }
      }
   }
   
   return xml_all_cp;
}

flam3_genome * flam3_parse_from_file(FILE *f, char *fname, int default_flag, int *ncps) {
   int i, c, slen = 5000;
   char *s;
   flam3_genome *ret;

   /* Incrementally read XML file into a string */
   s = malloc(slen);
   i = 0;
   do {
      c = getc(f);
      if (EOF == c)
         break;
      s[i++] = c;
      if (i == slen-1) {
         slen *= 2;
         s = realloc(s, slen);
      }
   } while (1);

   /* Null-terminate the read XML data */
   s[i] = 0;

   /* Parse the XML string */
   if (fname)
      ret = flam3_parse_xml2(s, fname, default_flag, ncps);
   else
      ret = flam3_parse_xml2(s, "stdin", default_flag, ncps);
      
   free(s);
   
   return(ret);

}

static void flam3_edit_print(FILE *f, xmlNodePtr editNode, int tabs, int formatting) {

   char *tab_string = "   ";
   int ti,strl;
   xmlAttrPtr att_ptr=NULL,cur_att=NULL;
   xmlNodePtr chld_ptr=NULL, cur_chld=NULL;
   int edit_or_sheep = 0, indent_printed = 0;
   char *ai;
   int tablim = argi("print_edit_depth",0);

   char *att_str,*cont_str,*cpy_string;

   if (tablim>0 && tabs>tablim)
   return;

   /* If this node is an XML_ELEMENT_NODE, print it and it's attributes */
   if (editNode->type==XML_ELEMENT_NODE) {

      /* Print the node at the tab specified */
      if (formatting) {
         for (ti=0;ti<tabs;ti++)
            fprintf(f,"%s",tab_string);
      }

      fprintf(f,"<%s",editNode->name);

      /* This can either be an edit node or a sheep node */
      /* If it's an edit node, add one to the tab        */
      if (!xmlStrcmp(editNode->name, (const xmlChar *)"edit")) {
         edit_or_sheep = 1;
         tabs ++;
      } else if (!xmlStrcmp(editNode->name, (const xmlChar *)"sheep"))
         edit_or_sheep = 2;
      else
         edit_or_sheep = 0;


      /* Print the attributes */
      att_ptr = editNode->properties;

      for (cur_att = att_ptr; cur_att; cur_att = cur_att->next) {

         att_str = (char *) xmlGetProp(editNode,cur_att->name);
         fprintf(f," %s=\"%s\"",cur_att->name,att_str);
         xmlFree(att_str);
      }

      /* Does this node have children? */
      if (!editNode->children || (tablim>0 && tabs>tablim)) {
         /* Close the tag and subtract the tab */
         fprintf(f,"/>");
         if (formatting)
            fprintf(f,"\n");
         tabs--;
      } else {

         /* Close the tag */
         fprintf(f,">");

         if (formatting)
            fprintf(f,"\n");

         /* Loop through the children and print them */
         chld_ptr = editNode->children;

         indent_printed = 0;

         for (cur_chld=chld_ptr; cur_chld; cur_chld = cur_chld->next) {

            /* If child is an element, indent first and then print it. */
            if (cur_chld->type==XML_ELEMENT_NODE &&
               (!xmlStrcmp(cur_chld->name, (const xmlChar *)"edit") ||
      (!xmlStrcmp(cur_chld->name, (const xmlChar *)"sheep")))) {

               if (indent_printed) {
                  indent_printed = 0;
                  fprintf(f,"\n");
               }

               flam3_edit_print(f, cur_chld, tabs, 1);

            } else {

               /* Child is a text node.  We don't want to indent more than once. */
               if (xmlIsBlankNode(cur_chld))
                  continue;

               if (indent_printed==0 && formatting==1) {
                  for (ti=0;ti<tabs;ti++)
                     fprintf(f,"%s",tab_string);
                  indent_printed = 1;
               }

               /* Print nodes without formatting. */
               flam3_edit_print(f, cur_chld, tabs, 0);

            }
         }

         if (indent_printed && formatting)
            fprintf(f,"\n");

         /* Tab out. */
         tabs --;
         if (formatting) {
            for (ti=0;ti<tabs;ti++)
               fprintf(f,"%s",tab_string);
         }

         /* Close the tag */
         fprintf(f,"</%s>",editNode->name);

         if (formatting) {
            fprintf(f,"\n");
         }
      }

   } else if (editNode->type==XML_TEXT_NODE) {

      /* Print text node */
      cont_str = (char *) xmlNodeGetContent(editNode);
      cpy_string = &(cont_str[0]);
      while (isspace(*cpy_string))
         cpy_string++;

      strl = (int)strlen(cont_str)-1;

      while (isspace(cont_str[strl]))
         strl--;

      cont_str[strl+1] = 0;

      fprintf(f,"%s",cpy_string);

   }
}

void flam3_apply_template(flam3_genome *cp, flam3_genome *templ) {

   /* Check for invalid values - only replace those with valid ones */
   if (templ->background[0] >= 0)
      cp->background[0] = templ->background[0];
   if (templ->background[1] >= 0)
      cp->background[1] = templ->background[1];
   if (templ->background[1] >= 0)
      cp->background[2] = templ->background[2];
   if (templ->zoom < 999999998)
      cp->zoom = templ->zoom;
   if (templ->spatial_oversample > 0)
      cp->spatial_oversample = templ->spatial_oversample;
   if (templ->spatial_filter_radius >= 0)
      cp->spatial_filter_radius = templ->spatial_filter_radius;
   if (templ->sample_density > 0)
      cp->sample_density = templ->sample_density;
   if (templ->nbatches > 0)
      cp->nbatches = templ->nbatches;
   if (templ->ntemporal_samples > 0)
      cp->ntemporal_samples = templ->ntemporal_samples;
   if (templ->width > 0) {
      /* preserving scale should be an option */
      cp->pixels_per_unit = cp->pixels_per_unit * templ->width / cp->width;
      cp->width = templ->width;
   }
   if (templ->height > 0)
      cp->height = templ->height;
   if (templ->estimator >= 0)
      cp->estimator = templ->estimator;
   if (templ->estimator_minimum >= 0)
      cp->estimator_minimum = templ->estimator_minimum;
   if (templ->estimator_curve >= 0)
      cp->estimator_curve = templ->estimator_curve;
   if (templ->gam_lin_thresh >= 0)
      cp->gam_lin_thresh = templ->gam_lin_thresh;
//   if (templ->motion_exp != -999)
//      cp->motion_exp = templ->motion_exp;
   if (templ->nbatches>0)
      cp->nbatches = templ->nbatches;
   if (templ->ntemporal_samples>0)
      cp->ntemporal_samples = templ->ntemporal_samples;
   if (templ->spatial_filter_select>0)
      cp->spatial_filter_select = templ->spatial_filter_select;
   if (templ->interpolation >= 0)
      cp->interpolation = templ->interpolation;
   if (templ->interpolation_type >= 0)
      cp->interpolation_type = templ->interpolation_type;
   if (templ->temporal_filter_type >= 0)
      cp->temporal_filter_type = templ->temporal_filter_type;
   if (templ->temporal_filter_width > 0)
      cp->temporal_filter_width = templ->temporal_filter_width;
   if (templ->temporal_filter_exp > -900)
      cp->temporal_filter_exp = templ->temporal_filter_exp;
   if (templ->highlight_power >=0)
      cp->highlight_power = templ->highlight_power;
   if (templ->palette_mode >= 0)
      cp->palette_mode = templ->palette_mode;

}

char *flam3_print_to_string(flam3_genome *cp) {

   FILE *tmpflame;
   long stringbytes;
   char *genome_string;
   
   tmpflame = tmpfile();
   flam3_print(tmpflame,cp,NULL,flam3_dont_print_edits);
   stringbytes = ftell(tmpflame);
   fseek(tmpflame,0L, SEEK_SET);
   genome_string = (char *)calloc(stringbytes+1,1);
   if (stringbytes != fread(genome_string, 1, stringbytes, tmpflame)) {
       perror("FLAM3: reading string from temp file");
   }
   fclose(tmpflame);
   
   return(genome_string);
}
   

void flam3_print(FILE *f, flam3_genome *cp, char *extra_attributes, int print_edits) {
   int i, j,numstd;

   fprintf(f, "<flame time=\"%g\"", cp->time);
   
   if (cp->flame_name[0]!=0)
      fprintf(f, " name=\"%s\"",cp->flame_name);
   
   fprintf(f, " size=\"%d %d\"", cp->width, cp->height);
   fprintf(f, " center=\"%g %g\"", cp->center[0], cp->center[1]);
   fprintf(f, " scale=\"%g\"", cp->pixels_per_unit);

   if (cp->zoom != 0.0)
      fprintf(f, " zoom=\"%g\"", cp->zoom);

   fprintf(f, " rotate=\"%g\"", cp->rotate);
   fprintf(f, " supersample=\"%d\"", cp->spatial_oversample);
   fprintf(f, " filter=\"%g\"", cp->spatial_filter_radius);

   /* Need to print the correct kernel to use */
   if (cp->spatial_filter_select == flam3_gaussian_kernel)
      fprintf(f, " filter_shape=\"gaussian\"");
   else if (cp->spatial_filter_select == flam3_hermite_kernel)
      fprintf(f, " filter_shape=\"hermite\"");
   else if (cp->spatial_filter_select == flam3_box_kernel)
      fprintf(f, " filter_shape=\"box\"");
   else if (cp->spatial_filter_select == flam3_triangle_kernel)
      fprintf(f, " filter_shape=\"triangle\"");
   else if (cp->spatial_filter_select == flam3_bell_kernel)
      fprintf(f, " filter_shape=\"bell\"");
   else if (cp->spatial_filter_select == flam3_b_spline_kernel)
      fprintf(f, " filter_shape=\"bspline\"");
   else if (cp->spatial_filter_select == flam3_mitchell_kernel)
      fprintf(f, " filter_shape=\"mitchell\"");
   else if (cp->spatial_filter_select == flam3_blackman_kernel)
      fprintf(f, " filter_shape=\"blackman\"");
   else if (cp->spatial_filter_select == flam3_catrom_kernel)
      fprintf(f, " filter_shape=\"catrom\"");
   else if (cp->spatial_filter_select == flam3_hanning_kernel)
      fprintf(f, " filter_shape=\"hanning\"");
   else if (cp->spatial_filter_select == flam3_hamming_kernel)
      fprintf(f, " filter_shape=\"hamming\"");
   else if (cp->spatial_filter_select == flam3_lanczos3_kernel)
      fprintf(f, " filter_shape=\"lanczos3\"");
   else if (cp->spatial_filter_select == flam3_lanczos2_kernel)
      fprintf(f, " filter_shape=\"lanczos2\"");
   else if (cp->spatial_filter_select == flam3_quadratic_kernel)
      fprintf(f, " filter_shape=\"quadratic\"");

   if (cp->temporal_filter_type == flam3_temporal_box)
      fprintf(f, " temporal_filter_type=\"box\"");
   else if (cp->temporal_filter_type == flam3_temporal_gaussian)
      fprintf(f, " temporal_filter_type=\"gaussian\"");
   else if (cp->temporal_filter_type == flam3_temporal_exp)
      fprintf(f, " temporal_filter_type=\"exp\" temporal_filter_exp=\"%g\"",cp->temporal_filter_exp);

   fprintf(f, " temporal_filter_width=\"%g\"",cp->temporal_filter_width);
   


   fprintf(f, " quality=\"%g\"", cp->sample_density);
   fprintf(f, " passes=\"%d\"", cp->nbatches);
   fprintf(f, " temporal_samples=\"%d\"", cp->ntemporal_samples);
   fprintf(f, " background=\"%g %g %g\"",
      cp->background[0], cp->background[1], cp->background[2]);
   fprintf(f, " brightness=\"%g\"", cp->brightness);
   fprintf(f, " gamma=\"%g\"", cp->gamma);
   fprintf(f, " highlight_power=\"%g\"", cp->highlight_power);
   fprintf(f, " vibrancy=\"%g\"", cp->vibrancy);
   fprintf(f, " estimator_radius=\"%g\" estimator_minimum=\"%g\" estimator_curve=\"%g\"",
      cp->estimator, cp->estimator_minimum, cp->estimator_curve);
   fprintf(f, " gamma_threshold=\"%g\"", cp->gam_lin_thresh);
   
   if (flam3_palette_mode_step == cp->palette_mode)
      fprintf(f, " palette_mode=\"step\"");
   else if (flam3_palette_mode_linear == cp->palette_mode)
      fprintf(f, " palette_mode=\"linear\"");

   if (flam3_interpolation_linear != cp->interpolation)
       fprintf(f, " interpolation=\"smooth\"");
       
   if (flam3_inttype_linear == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"linear\"");
   else if (flam3_inttype_log == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"log\"");
   else if (flam3_inttype_compat == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"old\"");
   else if (flam3_inttype_older == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"older\"");


   if (flam3_palette_interpolation_hsv != cp->palette_interpolation)
       fprintf(f, " palette_interpolation=\"sweep\"");

//   if (cp->motion_exp != 0.0)
//      fprintf(f, " motion_exponent=\"%g\"", cp->motion_exp);


   if (extra_attributes)
      fprintf(f, " %s", extra_attributes);

   fprintf(f, ">\n");

   if (cp->symmetry)
      fprintf(f, "   <symmetry kind=\"%d\"/>\n", cp->symmetry);
      
   //print_chaos(cp);
   
   numstd=-1;   
   for (i = 0; i < cp->num_xforms; i++) {
      int blob_var=0,pdj_var=0,fan2_var=0,rings2_var=0,perspective_var=0;
      int juliaN_var=0,juliaScope_var=0,radialBlur_var=0,pie_var=0,disc2_var=0;
      int ngon_var=0,curl_var=0,rectangles_var=0,supershape_var=0;
      int flower_var=0,conic_var=0,parabola_var=0,bent2_var=0,bipolar_var=0;
      int cell_var=0,cpow_var=0,curve_var=0,escher_var=0,lazys_var=0;
      int modulus_var=0,oscope_var=0,popcorn2_var=0,separation_var=0;
      int split_var=0,splits_var=0,stripes_var=0,wedge_var=0,wedgeJ_var=0;
      int wedgeS_var=0,whorl_var=0,waves2_var=0;
      if ( !(cp->symmetry &&  cp->xform[i].color_speed == 1.0 && cp->xform[i].animate == 1.0)) {

         if (i==cp->final_xform_index )
            fprintf(f, "   <finalxform color=\"%g\" ", cp->xform[i].color);
         else {
            numstd++;
            fprintf(f, "   <xform weight=\"%g\" color=\"%g\" ", cp->xform[i].density, cp->xform[i].color);
         }
         
         fprintf(f, "color_speed=\"%g\" ", cp->xform[i].color_speed);
         
         if (i!=cp->final_xform_index)
            fprintf(f, "animate=\"%g\" ", cp->xform[i].animate);
            
//         fprintf(f, "symmetry=\"%g\" ", cp->xform[i].symmetry);

         for (j = 0; j < flam3_nvariations; j++) {
            double v = cp->xform[i].var[j];
            if (0.0 != v) {
               fprintf(f, "%s=\"%g\" ", flam3_variation_names[j], v);
               if (j==VAR_BLOB)
                  blob_var=1;
               else if (j==VAR_PDJ)
                  pdj_var=1;
               else if (j==VAR_FAN2)
                  fan2_var=1;
               else if (j==VAR_RINGS2)
                  rings2_var=1;
               else if (j==VAR_PERSPECTIVE)
                  perspective_var=1;
               else if (j==VAR_JULIAN)
                  juliaN_var=1;
               else if (j==VAR_JULIASCOPE)
                  juliaScope_var=1;
               else if (j==VAR_RADIAL_BLUR)
                  radialBlur_var=1;
               else if (j==VAR_PIE)
                  pie_var=1;
               else if (j==VAR_NGON)
                  ngon_var=1;
               else if (j==VAR_CURL)
                  curl_var=1;
               else if (j==VAR_RECTANGLES)
                  rectangles_var=1;
               else if (j==VAR_DISC2)
                  disc2_var=1;
               else if (j==VAR_SUPER_SHAPE)
                  supershape_var=1;
               else if (j==VAR_FLOWER)
                  flower_var=1;
               else if (j==VAR_CONIC)
                  conic_var=1;
               else if (j==VAR_PARABOLA)
                  parabola_var=1;
               else if (j==VAR_BENT2)
                  bent2_var=1;
               else if (j==VAR_BIPOLAR)
                  bipolar_var=1;
               else if (j==VAR_CELL)
                  cell_var=1;
               else if (j==VAR_CPOW)
                  cpow_var=1;
               else if (j==VAR_CURVE)
                  curve_var=1;
               else if (j==VAR_ESCHER)
                  escher_var=1;
               else if (j==VAR_LAZYSUSAN)
                  lazys_var=1;
               else if (j==VAR_MODULUS)
                  modulus_var=1;
               else if (j==VAR_OSCILLOSCOPE)
                  oscope_var=1;
               else if (j==VAR_POPCORN2)
                  popcorn2_var=1;
               else if (j==VAR_SPLIT)
                  split_var=1;
               else if (j==VAR_SPLITS)
                  splits_var=1;
               else if (j==VAR_STRIPES)
                  stripes_var=1;
               else if (j==VAR_WEDGE)
                  wedge_var=1;
               else if (j==VAR_WEDGE_JULIA)
                  wedgeJ_var=1;
               else if (j==VAR_WEDGE_SPH)
                  wedgeS_var=1;
               else if (j==VAR_WHORL)
                  whorl_var=1;
               else if (j==VAR_WAVES2)
                  waves2_var=1;
            }
         }

         if (blob_var==1) {
            fprintf(f, "blob_low=\"%g\" ", cp->xform[i].blob_low);
            fprintf(f, "blob_high=\"%g\" ", cp->xform[i].blob_high);
            fprintf(f, "blob_waves=\"%g\" ", cp->xform[i].blob_waves);
         }

         if (pdj_var==1) {
            fprintf(f, "pdj_a=\"%g\" ", cp->xform[i].pdj_a);
            fprintf(f, "pdj_b=\"%g\" ", cp->xform[i].pdj_b);
            fprintf(f, "pdj_c=\"%g\" ", cp->xform[i].pdj_c);
            fprintf(f, "pdj_d=\"%g\" ", cp->xform[i].pdj_d);
         }

         if (fan2_var==1) {
            fprintf(f, "fan2_x=\"%g\" ", cp->xform[i].fan2_x);
            fprintf(f, "fan2_y=\"%g\" ", cp->xform[i].fan2_y);
         }

         if (rings2_var==1) {
            fprintf(f, "rings2_val=\"%g\" ", cp->xform[i].rings2_val);
         }

         if (perspective_var==1) {
            fprintf(f, "perspective_angle=\"%g\" ", cp->xform[i].perspective_angle);
            fprintf(f, "perspective_dist=\"%g\" ", cp->xform[i].perspective_dist);
         }

         if (juliaN_var==1) {
            fprintf(f, "julian_power=\"%g\" ", cp->xform[i].juliaN_power);
            fprintf(f, "julian_dist=\"%g\" ", cp->xform[i].juliaN_dist);
         }

         if (juliaScope_var==1) {
            fprintf(f, "juliascope_power=\"%g\" ", cp->xform[i].juliaScope_power);
            fprintf(f, "juliascope_dist=\"%g\" ", cp->xform[i].juliaScope_dist);
         }

         if (radialBlur_var==1) {
            fprintf(f, "radial_blur_angle=\"%g\" ", cp->xform[i].radialBlur_angle);
         }

         if (pie_var==1) {
            fprintf(f, "pie_slices=\"%g\" ", cp->xform[i].pie_slices);
            fprintf(f, "pie_rotation=\"%g\" ", cp->xform[i].pie_rotation);
            fprintf(f, "pie_thickness=\"%g\" ", cp->xform[i].pie_thickness);
         }

         if (ngon_var==1) {
            fprintf(f, "ngon_sides=\"%g\" ", cp->xform[i].ngon_sides);
            fprintf(f, "ngon_power=\"%g\" ", cp->xform[i].ngon_power);
            fprintf(f, "ngon_corners=\"%g\" ", cp->xform[i].ngon_corners);
            fprintf(f, "ngon_circle=\"%g\" ", cp->xform[i].ngon_circle);
         }

         if (curl_var==1) {
            fprintf(f, "curl_c1=\"%g\" ", cp->xform[i].curl_c1);
            fprintf(f, "curl_c2=\"%g\" ", cp->xform[i].curl_c2);
         }

         if (rectangles_var==1) {
            fprintf(f, "rectangles_x=\"%g\" ", cp->xform[i].rectangles_x);
            fprintf(f, "rectangles_y=\"%g\" ", cp->xform[i].rectangles_y);
         }

         if (disc2_var==1) {
            fprintf(f, "disc2_rot=\"%g\" ", cp->xform[i].disc2_rot);
            fprintf(f, "disc2_twist=\"%g\" ", cp->xform[i].disc2_twist);
         }

         if (supershape_var==1) {
            fprintf(f, "super_shape_rnd=\"%g\" ", cp->xform[i].supershape_rnd);
            fprintf(f, "super_shape_m=\"%g\" ", cp->xform[i].supershape_m);
            fprintf(f, "super_shape_n1=\"%g\" ", cp->xform[i].supershape_n1);
            fprintf(f, "super_shape_n2=\"%g\" ", cp->xform[i].supershape_n2);
            fprintf(f, "super_shape_n3=\"%g\" ", cp->xform[i].supershape_n3);
            fprintf(f, "super_shape_holes=\"%g\" ", cp->xform[i].supershape_holes);
         }

         if (flower_var==1) {
            fprintf(f, "flower_petals=\"%g\" ", cp->xform[i].flower_petals);
            fprintf(f, "flower_holes=\"%g\" ", cp->xform[i].flower_holes);
         }

         if (conic_var==1) {
            fprintf(f, "conic_eccentricity=\"%g\" ", cp->xform[i].conic_eccen);
            fprintf(f, "conic_holes=\"%g\" ", cp->xform[i].conic_holes);
         }

         if (parabola_var==1) {
            fprintf(f, "parabola_height=\"%g\" ", cp->xform[i].parabola_height);
            fprintf(f, "parabola_width=\"%g\" ", cp->xform[i].parabola_width);
         }
         
         if (bent2_var==1) {
            fprintf(f, "bent2_x=\"%g\" ", cp->xform[i].bent2_x);
            fprintf(f, "bent2_y=\"%g\" ", cp->xform[i].bent2_y);
         }
         
         if (bipolar_var==1) {
            fprintf(f, "bipolar_shift=\"%g\" ", cp->xform[i].bipolar_shift);
         }

         if (cell_var==1) {
            fprintf(f, "cell_size=\"%g\" ", cp->xform[i].cell_size);
         }

         if (cpow_var==1) {
            fprintf(f, "cpow_i=\"%g\" ", cp->xform[i].cpow_i);
            fprintf(f, "cpow_r=\"%g\" ", cp->xform[i].cpow_r);
            fprintf(f, "cpow_power=\"%g\" ", cp->xform[i].cpow_power);
         }

         if (curve_var==1) {
            fprintf(f, "curve_xamp=\"%g\" ", cp->xform[i].curve_xamp);
            fprintf(f, "curve_yamp=\"%g\" ", cp->xform[i].curve_yamp);
            fprintf(f, "curve_xlength=\"%g\" ", cp->xform[i].curve_xlength);
            fprintf(f, "curve_ylength=\"%g\" ", cp->xform[i].curve_ylength);
         }

         if (escher_var==1) {
            fprintf(f, "escher_beta=\"%g\" ", cp->xform[i].escher_beta);
         }

         if (lazys_var==1) {
            fprintf(f, "lazysusan_x=\"%g\" ", cp->xform[i].lazysusan_x);
            fprintf(f, "lazysusan_y=\"%g\" ", cp->xform[i].lazysusan_y);
            fprintf(f, "lazysusan_spin=\"%g\" ", cp->xform[i].lazysusan_spin);
            fprintf(f, "lazysusan_space=\"%g\" ", cp->xform[i].lazysusan_space);
            fprintf(f, "lazysusan_twist=\"%g\" ", cp->xform[i].lazysusan_twist);
         }

         if (modulus_var==1) {
            fprintf(f, "modulus_x=\"%g\" ", cp->xform[i].modulus_x);
            fprintf(f, "modulus_y=\"%g\" ", cp->xform[i].modulus_y);
         }

         if (oscope_var==1) {
            fprintf(f, "oscope_separation=\"%g\" ", cp->xform[i].oscope_separation);
            fprintf(f, "oscope_frequency=\"%g\" ", cp->xform[i].oscope_frequency);
            fprintf(f, "oscope_amplitude=\"%g\" ", cp->xform[i].oscope_amplitude);
            fprintf(f, "oscope_damping=\"%g\" ", cp->xform[i].oscope_damping);
         }

         if (popcorn2_var==1) {
            fprintf(f, "popcorn2_x=\"%g\" ", cp->xform[i].popcorn2_x);
            fprintf(f, "popcorn2_y=\"%g\" ", cp->xform[i].popcorn2_y);
            fprintf(f, "popcorn2_c=\"%g\" ", cp->xform[i].popcorn2_c);
         }

         if (separation_var==1) {
            fprintf(f, "separation_x=\"%g\" ", cp->xform[i].separation_x);
            fprintf(f, "separation_y=\"%g\" ", cp->xform[i].separation_y);
            fprintf(f, "separation_xinside=\"%g\" ", cp->xform[i].separation_xinside);
            fprintf(f, "separation_yinside=\"%g\" ", cp->xform[i].separation_yinside);
         }

         if (split_var==1) {
            fprintf(f, "split_xsize=\"%g\" ", cp->xform[i].split_xsize);
            fprintf(f, "split_ysize=\"%g\" ", cp->xform[i].split_ysize);
         }

         if (splits_var==1) {
            fprintf(f, "splits_x=\"%g\" ", cp->xform[i].splits_x);
            fprintf(f, "splits_y=\"%g\" ", cp->xform[i].splits_y);
         }

         if (stripes_var==1) {
            fprintf(f, "stripes_space=\"%g\" ", cp->xform[i].stripes_space);
            fprintf(f, "stripes_warp=\"%g\" ", cp->xform[i].stripes_warp);
         }

         if (wedge_var==1) {
            fprintf(f, "wedge_angle=\"%g\" ", cp->xform[i].wedge_angle);
            fprintf(f, "wedge_hole=\"%g\" ", cp->xform[i].wedge_hole);
            fprintf(f, "wedge_count=\"%g\" ", cp->xform[i].wedge_count);
            fprintf(f, "wedge_swirl=\"%g\" ", cp->xform[i].wedge_swirl);            
         }

         if (wedgeJ_var==1) {
            fprintf(f, "wedge_julia_angle=\"%g\" ", cp->xform[i].wedge_julia_angle);
            fprintf(f, "wedge_julia_count=\"%g\" ", cp->xform[i].wedge_julia_count);
            fprintf(f, "wedge_julia_power=\"%g\" ", cp->xform[i].wedge_julia_power);
            fprintf(f, "wedge_julia_dist=\"%g\" ", cp->xform[i].wedge_julia_dist);            
         }

         if (wedgeS_var==1) {
            fprintf(f, "wedge_sph_angle=\"%g\" ", cp->xform[i].wedge_sph_angle);
            fprintf(f, "wedge_sph_hole=\"%g\" ", cp->xform[i].wedge_sph_hole);
            fprintf(f, "wedge_sph_count=\"%g\" ", cp->xform[i].wedge_sph_count);
            fprintf(f, "wedge_sph_swirl=\"%g\" ", cp->xform[i].wedge_sph_swirl);            
         }

         if (whorl_var==1) {
            fprintf(f, "whorl_inside=\"%g\" ", cp->xform[i].whorl_inside);
            fprintf(f, "whorl_outside=\"%g\" ", cp->xform[i].whorl_outside);
         }

         if (waves2_var==1) {
            fprintf(f, "waves2_scalex=\"%g\" ", cp->xform[i].waves2_scalex);
            fprintf(f, "waves2_scaley=\"%g\" ", cp->xform[i].waves2_scaley);
            fprintf(f, "waves2_freqx=\"%g\" ", cp->xform[i].waves2_freqx);
            fprintf(f, "waves2_freqy=\"%g\" ", cp->xform[i].waves2_freqy);
         }

         fprintf(f, "coefs=\"");
         for (j = 0; j < 3; j++) {
            if (j) fprintf(f, " ");
            fprintf(f, "%g %g", cp->xform[i].c[j][0], cp->xform[i].c[j][1]);
         }
         fprintf(f, "\"");

         if (!id_matrix(cp->xform[i].post)) {
            fprintf(f, " post=\"");
            for (j = 0; j < 3; j++) {
               if (j) fprintf(f, " ");
               fprintf(f, "%g %g", cp->xform[i].post[j][0], cp->xform[i].post[j][1]);
            }
            fprintf(f, "\"");
         }
         
         if (i!=cp->final_xform_index ) {
         
            /* Print out the chaos row for this xform */
            int numcols = cp->num_xforms - (cp->final_xform_index >= 0);
            while (numcols > 0 && cp->chaos[i][numcols-1]==1.0)
               numcols--;
            
            if (numcols>0) {   
               fprintf(f, " chaos=\"");
               for (j=0;j<numcols;j++)
                  fprintf(f, "%g ",cp->chaos[i][j]);
               fprintf(f, "\"");
            }
            
            if (cp->xform[i].visibility<1.0)
               fprintf(f, " visibility=\"%g\"",cp->xform[i].visibility);
                
         }

         fprintf(f, "/>\n");

      }
   }

   for (i = 0; i < 256; i++) {
      double r, g, b;
      r = (cp->palette[i].color[0] * 255.0);
      g = (cp->palette[i].color[1] * 255.0);
      b = (cp->palette[i].color[2] * 255.0);
      if (getenv("intpalette"))
         fprintf(f, "   <color index=\"%d\" rgb=\"%d %d %d\"/>\n", i, (int)rint(r), (int)rint(g), (int)rint(b));
      else {
#ifdef USE_FLOAT_INDICES
         fprintf(f, "   <color index=\"%.10g\" rgb=\"%.6g %.6g %.6g\"/>\n", cp->palette[i].index, r, g, b);
#else
         fprintf(f, "   <color index=\"%d\" rgb=\"%.6g %.6g %.6g\"/>\n", i, r, g, b);
#endif
      }
   }

   if (cp->edits != NULL && print_edits==flam3_print_edits) {

      /* We need a custom script for printing these */
      /* and it needs to be recursive               */
      xmlNodePtr elem_node = xmlDocGetRootElement(cp->edits);
      flam3_edit_print(f,elem_node, 1, 1);
   }
   fprintf(f, "</flame>\n");

}

/* returns a uniform variable from 0 to 1 */
double flam3_random01() {
   return (random() & 0xfffffff) / (double) 0xfffffff;
}

double flam3_random11() {
   return ((random() & 0xfffffff) - 0x7ffffff) / (double) 0x7ffffff;
}

/* This function must be called prior to rendering a frame */
void flam3_init_frame(flam3_frame *f) {

   char *ai;
   char *isaac_seed = args("isaac_seed",NULL);
   long int default_isaac_seed = (long int)time(0);  

   /* Clear out the isaac state */
   memset(f->rc.randrsl, 0, RANDSIZ*sizeof(ub4));

   /* Set the isaac seed */
   if (NULL == isaac_seed) {
      int lp;
      /* No isaac seed specified.  Use the system time to initialize. */
      for (lp = 0; lp < RANDSIZ; lp++)
         f->rc.randrsl[lp] = default_isaac_seed;
   } else {
      /* Use the specified string */
      strncpy((char *)&f->rc.randrsl,(const char *)isaac_seed, RANDSIZ*sizeof(ub4));
   }

   /* Initialize the random number generator */
   irandinit(&f->rc,1);
}

/* returns uniform variable from ISAAC rng */
double flam3_random_isaac_01(randctx *ct) {
   return ((int)irand(ct) & 0xfffffff) / (double) 0xfffffff;
}

double flam3_random_isaac_11(randctx *ct) {
   return (((int)irand(ct) & 0xfffffff) - 0x7ffffff) / (double) 0x7ffffff;
}

int flam3_random_bit() {
  /* might not be threadsafe */
  static int n = 0;
  static int l;
  if (0 == n) {
    l = random();
    n = 20;
  } else {
    l = l >> 1;
    n--;
  }
  return l & 1;
}

int flam3_random_isaac_bit(randctx *ct) {
   int tmp = irand(ct);
   return tmp & 1;
}

int flam3_parse_hexformat_colors(char *colstr, flam3_genome *cp, int numcolors, int chan) {

   int c_idx=0;
   int col_count=0;
   int r,g,b;
   int sscanf_ret;
   char tmps[2];
   
   /* Strip whitespace prior to first color */
   while (isspace( (int)colstr[c_idx]))
      c_idx++;

   do {

      /* Parse an RGB triplet at a time... */
      if (chan==3)
         sscanf_ret = sscanf(&(colstr[c_idx]),"%2x%2x%2x",&r,&g,&b);
      else
         sscanf_ret = sscanf(&(colstr[c_idx]),"00%2x%2x%2x",&r,&g,&b);

      if (sscanf_ret != 3) {
         fprintf(stderr, "Error:  Problem reading hexadecimal color data.\n");
         return(1);
      }

      c_idx += 2*chan;

      while (isspace( (int)colstr[c_idx]))
         c_idx++;

      cp->palette[col_count].color[0] = r / 255.0;
      cp->palette[col_count].color[1] = g / 255.0;
      cp->palette[col_count].color[2] = b / 255.0;
      cp->palette[col_count].color[3] = 1.0;
      cp->palette[col_count].index = col_count;

      col_count++;

   } while (col_count<numcolors);
   
   if (sscanf(&(colstr[c_idx]),"%1s",tmps)>0) {
      fprintf(stderr,"error: extra data at end of hex color data '%s'\n",&(colstr[c_idx]));
      return(1);
   }
   
   return(0);
}

/* sum of entries of vector to 1 */
static int normalize_vector(double *v, int n) {
   double t = 0.0;
   int i;
   for (i = 0; i < n; i++)
      t += v[i];
   if (0.0 == t) return 1;
   t = 1.0 / t;
   for (i = 0; i < n; i++)
      v[i] *= t;
   return 0;
}



static double round6(double x) {
  x *= 1e6;
  if (x < 0) x -= 1.0;
  return 1e-6*(int)(x+0.5);
}

/* sym=2 or more means rotational
   sym=1 means identity, ie no symmetry
   sym=0 means pick a random symmetry (maybe none)
   sym=-1 means bilateral (reflection)
   sym=-2 or less means rotational and reflective
*/
void flam3_add_symmetry(flam3_genome *cp, int sym) {
   int i, j, k;
   double a;
   int result = 0;
   flam3_xform tmp;
   int orig_xf = cp->num_xforms - (cp->final_xform_index >= 0);

   if (0 == sym) {
      static int sym_distrib[] = {
         -4, -3,
         -2, -2, -2,
         -1, -1, -1,
         2, 2, 2,
         3, 3,
         4, 4,
      };
      if (random()&1) {
         sym = random_distrib(sym_distrib);
      } else if (random()&31) {
         sym = (random()%13)-6;
      } else {
         sym = (random()%51)-25;
      }
   }

   if (1 == sym || 0 == sym) return;

   cp->symmetry = sym;

   if (sym < 0) {

      i = cp->num_xforms;
      flam3_add_xforms(cp,1,0,0);

      cp->xform[i].density = 1.0;
      cp->xform[i].color_speed = 1.0;
      cp->xform[i].animate = 1.0;
      cp->xform[i].var[0] = 1.0;
      for (j = 1; j < flam3_nvariations; j++)
         cp->xform[i].var[j] = 0;
      cp->xform[i].color = 1.0;
      cp->xform[i].c[0][0] = -1.0;
      cp->xform[i].c[0][1] = 0.0;
      cp->xform[i].c[1][0] = 0.0;
      cp->xform[i].c[1][1] = 1.0;
      cp->xform[i].c[2][0] = 0.0;
      cp->xform[i].c[2][1] = 0.0;

      result++;
      sym = -sym;
   }

   a = 2*M_PI/sym;

   for (k = 1; k < sym; k++) {

      i = cp->num_xforms;
      flam3_add_xforms(cp, 1, 0,0);

      cp->xform[i].density = 1.0;
      cp->xform[i].color_speed = 1.0;
      cp->xform[i].animate = 1.0;
      cp->xform[i].var[0] = 1.0;
      for (j = 1; j < flam3_nvariations; j++)
         cp->xform[i].var[j] = 0;
      cp->xform[i].color = (sym<3) ? 0.0 : ((k-1.0)/(sym-2.0));
      cp->xform[i].c[0][0] = round6(cos(k*a));
      cp->xform[i].c[0][1] = round6(sin(k*a));
      cp->xform[i].c[1][0] = round6(-cp->xform[i].c[0][1]);
      cp->xform[i].c[1][1] = cp->xform[i].c[0][0];
      cp->xform[i].c[2][0] = 0.0;
      cp->xform[i].c[2][1] = 0.0;

      result++;
   }   
   
   qsort((char *) &cp->xform[cp->num_xforms-result], result,
      sizeof(flam3_xform), compare_xforms);
      
}

static int random_var() {
  return random() % flam3_nvariations;
}

static int random_varn(int n) {
   return random() % n;
}

void flam3_random(flam3_genome *cp, int *ivars, int ivars_n, int sym, int spec_xforms) {

   int i, j, nxforms, var, samed, multid, samepost, postid, addfinal=0;
   int finum = -1;
   int n;
   double sum;

   static int xform_distrib[] = {
     2, 2, 2, 2,
     3, 3, 3, 3,
     4, 4, 4,
     5, 5,
     6
   };

   clear_cp(cp,flam3_defaults_on);

   cp->hue_rotation = (random()&7) ? 0.0 : flam3_random01();
   cp->palette_index = flam3_get_palette(flam3_palette_random, cp->palette, cp->hue_rotation);
   cp->time = 0.0;
   cp->interpolation = flam3_interpolation_linear;
   cp->palette_interpolation = flam3_palette_interpolation_hsv;

   /* Choose the number of xforms */
   if (spec_xforms>0) {
      nxforms = spec_xforms;
      flam3_add_xforms(cp,nxforms,0,0);
   } else {
      nxforms = random_distrib(xform_distrib);
      flam3_add_xforms(cp,nxforms,0,0);
      /* Add a final xform 15% of the time */
      addfinal = flam3_random01() < 0.15;
      if (addfinal) {
         flam3_add_xforms(cp,1,0,1);
         nxforms = nxforms + addfinal;
         finum = nxforms-1;
      }
   }   

   /* If first input variation is 'flam3_variation_random' */
   /* choose one to use or decide to use multiple    */
   if (flam3_variation_random == ivars[0]) {
      if (flam3_random_bit()) {
         var = random_var();
      } else {
         var = flam3_variation_random;
      }
   } else {
      var = flam3_variation_random_fromspecified;
   }

   samed = flam3_random_bit();
   multid = flam3_random_bit();
   postid = flam3_random01() < 0.6;
   samepost = flam3_random_bit();

   /* Loop over xforms */
   for (i = 0; i < nxforms; i++) {
      int j, k;
      cp->xform[i].density = 1.0 / nxforms;
      cp->xform[i].color = i&1;
      cp->xform[i].color_speed = 0.0;
      cp->xform[i].animate = 0.0;
      for (j = 0; j < 3; j++) {
         for (k = 0; k < 2; k++) {
            cp->xform[i].c[j][k] = flam3_random11();
            cp->xform[i].post[j][k] = (double)(k==j);
         }
      }
      
      if ( i != finum ) {
      
         if (!postid) {

            for (j = 0; j < 3; j++)
            for (k = 0; k < 2; k++) {
               if (samepost || (i==0))
                  cp->xform[i].post[j][k] = flam3_random11();
               else
                  cp->xform[i].post[j][k] = cp->xform[0].post[j][k];
            }
         }

         /* Clear all variation coefs */
         for (j = 0; j < flam3_nvariations; j++)
            cp->xform[i].var[j] = 0.0;

         if (flam3_variation_random != var && 
               flam3_variation_random_fromspecified != var) {

            /* Use only one variation specified for all xforms */
            cp->xform[i].var[var] = 1.0;

         } else if (multid && flam3_variation_random == var) {

           /* Choose a random var for this xform */
             cp->xform[i].var[random_var()] = 1.0;

         } else {

            if (samed && i > 0) {

               /* Copy the same variations from the previous xform */
               for (j = 0; j < flam3_nvariations; j++) {
                  cp->xform[i].var[j] = cp->xform[i-1].var[j];
                  flam3_copy_params(&(cp->xform[i]),&(cp->xform[i-1]),j);
               }

            } else {

               /* Choose a random number of vars to use, at least 2 */
               /* but less than flam3_nvariations.Probability leans */
               /* towards fewer variations.                         */
               n = 2;
               while ((flam3_random_bit()) && (n<flam3_nvariations))
                  n++;
               
               /* Randomly choose n variations, and change their weights. */
               /* A var can be selected more than once, further reducing  */
               /* the probability that multiple vars are used.            */
               for (j = 0; j < n; j++) {
                  if (flam3_variation_random_fromspecified != var)
                     cp->xform[i].var[random_var()] = flam3_random01();
                  else
                     cp->xform[i].var[ivars[random_varn(ivars_n)]] = flam3_random01();
               }

               /* Normalize weights to 1.0 total. */
               sum = 0.0;
               for (j = 0; j < flam3_nvariations; j++)
                  sum += cp->xform[i].var[j];
               if (sum == 0.0)
                  cp->xform[i].var[random_var()] = 1.0;
               else {
                  for (j = 0; j < flam3_nvariations; j++)
                     cp->xform[i].var[j] /= sum;
               }
            }
         }
      } else {
         /* Handle final xform randomness. */
         n = 1;
         if (flam3_random_bit()) n++;
         
         /* Randomly choose n variations, and change their weights. */
         /* A var can be selected more than once, further reducing  */
         /* the probability that multiple vars are used.            */
         for (j = 0; j < n; j++) {
            if (flam3_variation_random_fromspecified != var)
               cp->xform[i].var[random_var()] = flam3_random01();
            else
               cp->xform[i].var[ivars[random_varn(ivars_n)]] = flam3_random01();
         }

         /* Normalize weights to 1.0 total. */
         sum = 0.0;
         for (j = 0; j < flam3_nvariations; j++)
            sum += cp->xform[i].var[j];
         if (sum == 0.0)
            cp->xform[i].var[random_var()] = 1.0;
         else {
            for (j = 0; j < flam3_nvariations; j++)
               cp->xform[i].var[j] /= sum;
         }
      }  

      
      if (cp->xform[i].var[VAR_WAVES] > 0) {
         waves_precalc(&(cp->xform[i]));
      }

      /* Generate random params for parametric variations, if selected. */
      if (cp->xform[i].var[VAR_BLOB] > 0) {
         /* Create random params for blob */
         cp->xform[i].blob_low = 0.2 + 0.5 * flam3_random01();
         cp->xform[i].blob_high = 0.8 + 0.4 * flam3_random01();
         cp->xform[i].blob_waves = (int)(2 + 5 * flam3_random01());
      }

      if (cp->xform[i].var[VAR_PDJ] > 0) {
         /* Create random params for PDJ */
         cp->xform[i].pdj_a = 3.0 * flam3_random11();
         cp->xform[i].pdj_b = 3.0 * flam3_random11();
         cp->xform[i].pdj_c = 3.0 * flam3_random11();
         cp->xform[i].pdj_d = 3.0 * flam3_random11();
      }

      if (cp->xform[i].var[VAR_FAN2] > 0) {
         /* Create random params for fan2 */
         cp->xform[i].fan2_x = flam3_random11();
         cp->xform[i].fan2_y = flam3_random11();
      }

      if (cp->xform[i].var[VAR_RINGS2] > 0) {
         /* Create random params for rings2 */
         cp->xform[i].rings2_val = 2*flam3_random01();
      }

      if (cp->xform[i].var[VAR_PERSPECTIVE] > 0) {

         /* Create random params for perspective */
         cp->xform[i].perspective_angle = flam3_random01();
         cp->xform[i].perspective_dist = 2*flam3_random01() + 1.0;

      }

      if (cp->xform[i].var[VAR_JULIAN] > 0) {

         /* Create random params for juliaN */
         cp->xform[i].juliaN_power = (int)(5*flam3_random01() + 2);
         cp->xform[i].juliaN_dist = 1.0;

      }

      if (cp->xform[i].var[VAR_JULIASCOPE] > 0) {

         /* Create random params for juliaScope */
         cp->xform[i].juliaScope_power = (int)(5*flam3_random01() + 2);
         cp->xform[i].juliaScope_dist = 1.0;

      }

      if (cp->xform[i].var[VAR_RADIAL_BLUR] > 0) {

         /* Create random params for radialBlur */
         cp->xform[i].radialBlur_angle = (2 * flam3_random01() - 1);

      }

      if (cp->xform[i].var[VAR_PIE] > 0) {
         /* Create random params for pie */
         cp->xform[i].pie_slices = (int) 10.0*flam3_random01();
         cp->xform[i].pie_thickness = flam3_random01();
         cp->xform[i].pie_rotation = 2.0 * M_PI * flam3_random11();
      }

      if (cp->xform[i].var[VAR_NGON] > 0) {
         /* Create random params for ngon */
         cp->xform[i].ngon_sides = (int) flam3_random01()* 10 + 3;
         cp->xform[i].ngon_power = 3*flam3_random01() + 1;
         cp->xform[i].ngon_circle = 3*flam3_random01();
         cp->xform[i].ngon_corners = 2*flam3_random01()*cp->xform[i].ngon_circle;
      }

      if (cp->xform[i].var[VAR_CURL] > 0) {
         /* Create random params for curl */
         cp->xform[i].curl_c1 = flam3_random01();
         cp->xform[i].curl_c2 = flam3_random01();
      }

      if (cp->xform[i].var[VAR_RECTANGLES] > 0) {
         /* Create random params for rectangles */
         cp->xform[i].rectangles_x = flam3_random01();
         cp->xform[i].rectangles_y = flam3_random01();
      }

      if (cp->xform[i].var[VAR_DISC2] > 0) {
      /* Create random params for disc2 */
      cp->xform[i].disc2_rot = 0.5 * flam3_random01();
      cp->xform[i].disc2_twist = 0.5 * flam3_random01();

      }

      if (cp->xform[i].var[VAR_SUPER_SHAPE] > 0) {
         /* Create random params for supershape */
         cp->xform[i].supershape_rnd = flam3_random01();
         cp->xform[i].supershape_m = (int) flam3_random01()*6;
         cp->xform[i].supershape_n1 = flam3_random01()*40;
         cp->xform[i].supershape_n2 = flam3_random01()*20;
         cp->xform[i].supershape_n3 = cp->xform[i].supershape_n2;
         cp->xform[i].supershape_holes = 0.0;
      }

      if (cp->xform[i].var[VAR_FLOWER] > 0) {
         /* Create random params for flower */
         cp->xform[i].flower_petals = 4 * flam3_random01();
         cp->xform[i].flower_holes = flam3_random01();
      }

      if (cp->xform[i].var[VAR_CONIC] > 0) {
         /* Create random params for conic */
         cp->xform[i].conic_eccen = flam3_random01();
         cp->xform[i].conic_holes = flam3_random01();
      }

      if (cp->xform[i].var[VAR_PARABOLA] > 0) {
         /* Create random params for parabola */
         cp->xform[i].parabola_height = 0.5 + flam3_random01();
         cp->xform[i].parabola_width = 0.5 + flam3_random01();
      }

      if (cp->xform[i].var[VAR_BENT2] > 0) {
         /* Create random params for bent2 */
         cp->xform[i].bent2_x = 3*(-0.5 + flam3_random01());
         cp->xform[i].bent2_y = 3*(-0.5 + flam3_random01());
      }

      if (cp->xform[i].var[VAR_BIPOLAR] > 0) {
         /* Create random params for bipolar */
         cp->xform[i].bipolar_shift = 2.0 * flam3_random01() - 1;
      }

      if (cp->xform[i].var[VAR_CELL] > 0) {
         /* Create random params for cell */
         cp->xform[i].cell_size = 2.0 * flam3_random01() + 0.5;
      }

      if (cp->xform[i].var[VAR_CPOW] > 0) {
         /* Create random params for cpow */
         cp->xform[i].cpow_r = 3.0 * flam3_random01();
         cp->xform[i].cpow_i = flam3_random01() - 0.5;
         cp->xform[i].cpow_power = (int)(5.0 * flam3_random01());
      }

      if (cp->xform[i].var[VAR_CURVE] > 0) {
         /* Create random params for curve */
         cp->xform[i].curve_xamp = 5 * (flam3_random01()-.5);
         cp->xform[i].curve_yamp = 4 * (flam3_random01()-.5);
         cp->xform[i].curve_xlength = 2 * (flam3_random01()+.5);
         cp->xform[i].curve_ylength = 2 * (flam3_random01()+.5);
      }

      if (cp->xform[i].var[VAR_ESCHER] > 0) {
         /* Create random params for escher */
         cp->xform[i].escher_beta = M_PI * flam3_random11();
      }

      if (cp->xform[i].var[VAR_LAZYSUSAN] > 0) {
         /* Create random params for lazysusan */
         cp->xform[i].lazysusan_x = 2.0*flam3_random11();
         cp->xform[i].lazysusan_y = 2.0*flam3_random11();
         cp->xform[i].lazysusan_spin = M_PI*flam3_random11();
         cp->xform[i].lazysusan_space = 2.0*flam3_random11();
         cp->xform[i].lazysusan_twist = 2.0*flam3_random11();
      }

      if (cp->xform[i].var[VAR_MODULUS] > 0) {
         /* Create random params for modulus */
         cp->xform[i].modulus_x = flam3_random11();
         cp->xform[i].modulus_y = flam3_random11();
      }

      if (cp->xform[i].var[VAR_OSCILLOSCOPE] > 0) {
         /* Create random params for oscope */
         cp->xform[i].oscope_separation = 1.0 + flam3_random11();
         cp->xform[i].oscope_frequency = M_PI * flam3_random11();
         cp->xform[i].oscope_amplitude = 1.0 + 2 * flam3_random01();
         cp->xform[i].oscope_damping = flam3_random01();
      }

      if (cp->xform[i].var[VAR_POPCORN2] > 0) {
         /* Create random params for popcorn2 */
         cp->xform[i].popcorn2_x = 0.2 * flam3_random01();
         cp->xform[i].popcorn2_y = 0.2 * flam3_random01();
         cp->xform[i].popcorn2_c = 5 * flam3_random01();
      }

      if (cp->xform[i].var[VAR_SEPARATION] > 0) {
         /* Create random params for separation */
         cp->xform[i].separation_x = 1 + flam3_random11();
         cp->xform[i].separation_y = 1 + flam3_random11();
         cp->xform[i].separation_xinside = flam3_random11();
         cp->xform[i].separation_yinside = flam3_random11();
      }

      if (cp->xform[i].var[VAR_SPLIT] > 0) {
         /* Create random params for split */
         cp->xform[i].split_xsize = flam3_random11();
         cp->xform[i].split_ysize = flam3_random11();
      }

      if (cp->xform[i].var[VAR_SPLITS] > 0) {
         /* Create random params for splits */
         cp->xform[i].splits_x = flam3_random11();
         cp->xform[i].splits_y = flam3_random11();
      }

      if (cp->xform[i].var[VAR_STRIPES] > 0) {
         /* Create random params for stripes */
         cp->xform[i].stripes_space = flam3_random01();
         cp->xform[i].stripes_warp = 5*flam3_random01();
      }

      if (cp->xform[i].var[VAR_WEDGE] > 0) {
         /* Create random params for wedge */
         cp->xform[i].wedge_angle = M_PI*flam3_random01();
         cp->xform[i].wedge_hole = 0.5*flam3_random11();
         cp->xform[i].wedge_count = floor(5*flam3_random01())+1;
         cp->xform[i].wedge_swirl = flam3_random01();
      }

      if (cp->xform[i].var[VAR_WEDGE_JULIA] > 0) {

         /* Create random params for wedge_julia */
         cp->xform[i].wedge_julia_power = (int)(5*flam3_random01() + 2);
         cp->xform[i].wedge_julia_dist = 1.0;
         cp->xform[i].wedge_julia_count = (int)(3*flam3_random01() + 1);
         cp->xform[i].wedge_julia_angle = M_PI * flam3_random01();

      }

      if (cp->xform[i].var[VAR_WEDGE_SPH] > 0) {
         /* Create random params for wedge_sph */
         cp->xform[i].wedge_sph_angle = M_PI*flam3_random01();
         cp->xform[i].wedge_sph_hole = 0.5*flam3_random11();
         cp->xform[i].wedge_sph_count = floor(5*flam3_random01())+1;
         cp->xform[i].wedge_sph_swirl = flam3_random01();
      }

      if (cp->xform[i].var[VAR_WHORL] > 0) {
         /* Create random params for whorl */
         cp->xform[i].whorl_inside = flam3_random01();
         cp->xform[i].whorl_outside = flam3_random01();
      }

      if (cp->xform[i].var[VAR_WAVES2] > 0) {
         /* Create random params for waves2 */
         cp->xform[i].waves2_scalex = 0.5 + flam3_random01();
         cp->xform[i].waves2_scaley = 0.5 + flam3_random01();
         cp->xform[i].waves2_freqx = 4 * flam3_random01();
         cp->xform[i].waves2_freqy = 4 * flam3_random01();
      }         

   }

   /* Randomly add symmetry (but not if we've already added a final xform) */
   if (sym || (!(random()%4) && !addfinal))
      flam3_add_symmetry(cp, sym);
   else
      cp->symmetry = 0;

   //qsort((char *) cp->xform, (cp->num_xforms-addfinal), sizeof(flam3_xform), compare_xforms);


}


static int sort_by_x(const void *av, const void *bv) {
    double *a = (double *) av;
    double *b = (double *) bv;
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return 1;
    return 0;
}

static int sort_by_y(const void *av, const void *bv) {
    double *a = (double *) av;
    double *b = (double *) bv;
    if (a[1] < b[1]) return -1;
    if (a[1] > b[1]) return 1;
    return 0;
}


/* Memory helper functions because 

    Python on Windows uses the MSVCR71.dll version of the C Runtime and 
    mingw uses the MSVCRT.dll version. */

void *flam3_malloc(size_t size) {

   return (malloc(size));
   
}

void flam3_free(void *ptr) {

   free(ptr);
   
}

/*
 * find a 2d bounding box that does not enclose eps of the fractal density
 * in each compass direction.
 */
void flam3_estimate_bounding_box(flam3_genome *cp, double eps, int nsamples,
             double *bmin, double *bmax, randctx *rc) {
   int i;
   int low_target, high_target;
   double min[2], max[2];
   double *points;
   unsigned short *xform_distrib;

   if (nsamples <= 0) nsamples = 10000;
   low_target = (int)(nsamples * eps);
   high_target = nsamples - low_target;

   points = (double *) malloc(sizeof(double) * 4 * nsamples);
   points[0] = flam3_random_isaac_11(rc);
   points[1] = flam3_random_isaac_11(rc);
   points[2] = 0.0;
   points[3] = 0.0;

   prepare_xform_fn_ptrs(cp,rc);
   xform_distrib = flam3_create_xform_distrib(cp);
   flam3_iterate(cp, nsamples, 20, points, xform_distrib, rc);
   free(xform_distrib);

   min[0] = min[1] =  1e10;
   max[0] = max[1] = -1e10;

   for (i = 0; i < nsamples; i++) {
      double *p = &points[4*i];
      if (p[0] < min[0]) min[0] = p[0];
      if (p[1] < min[1]) min[1] = p[1];
      if (p[0] > max[0]) max[0] = p[0];
      if (p[1] > max[1]) max[1] = p[1];
   }

   if (low_target == 0) {
      bmin[0] = min[0];
      bmin[1] = min[1];
      bmax[0] = max[0];
      bmax[1] = max[1];
      free(points);
      return;
   }

   qsort(points, nsamples, sizeof(double) * 4, sort_by_x);
   bmin[0] = points[4 * low_target];
   bmax[0] = points[4 * high_target];

   qsort(points, nsamples, sizeof(double) * 4, sort_by_y);
   bmin[1] = points[4 * low_target + 1];
   bmax[1] = points[4 * high_target + 1];
   free(points);
}





typedef double bucket_double[5];
typedef double abucket_double[4];
typedef unsigned int bucket_int[5];
typedef unsigned int abucket_int[4];
typedef unsigned short bucket_short[5];
typedef unsigned short abucket_short[4];
typedef float bucket_float[5];
typedef float abucket_float[4];

#ifdef HAVE_GCC_64BIT_ATOMIC_OPS
static inline void
double_atomic_add(double *dest, double delta)
{
   uint64_t *int_ptr = (uint64_t *)dest;
   union {
      double dblval;
      uint64_t intval;
   } old_val, new_val;
   int success;

   do {
      old_val.dblval = *dest;
      new_val.dblval = old_val.dblval + delta;
      success = __sync_bool_compare_and_swap(
         int_ptr, old_val.intval, new_val.intval);
   } while (!success);
}
#endif /* HAVE_GCC_64BIT_ATOMIC_OPS */

#ifdef HAVE_GCC_ATOMIC_OPS
static inline void
float_atomic_add(float *dest, float delta)
{
   uint32_t *int_ptr = (uint32_t *)dest;
   union {
      float fltval;
      uint32_t intval;
   } old_val, new_val;
   int success;

   do {
      old_val.fltval = *dest;
      new_val.fltval = old_val.fltval + delta;
      success = __sync_bool_compare_and_swap(
         int_ptr, old_val.intval, new_val.intval);
   } while (!success);
}

static inline void
uint_atomic_add(unsigned int *dest, unsigned int delta)
{
   unsigned int old_val, new_val;
   int success;

   do {
      old_val = *dest;
      if (UINT_MAX - old_val > delta)
         new_val = old_val + delta;
      else
         new_val = UINT_MAX;
      success = __sync_bool_compare_and_swap(
         dest, old_val, new_val);
   } while (!success);
}

static inline void
ushort_atomic_add(unsigned short *dest, unsigned short delta)
{
   unsigned short old_val, new_val;
   int success;

   do {
      old_val = *dest;
      if (USHRT_MAX - old_val > delta)
         new_val = old_val + delta;
      else
         new_val = USHRT_MAX;
      success = __sync_bool_compare_and_swap(
         dest, old_val, new_val);
   } while (!success);
}
#endif /* HAVE_GCC_ATOMIC_OPS */

/* 64-bit datatypes */
#define B_ACCUM_T double
#define A_ACCUM_T double
#define bucket bucket_double
#define abucket abucket_double
#define abump_no_overflow(dest, delta) do {dest += delta;} while (0)
#define add_c_to_accum(acc,i,ii,j,jj,wid,hgt,c) do { \
   if ( (j) + (jj) >=0 && (j) + (jj) < (hgt) && (i) + (ii) >=0 && (i) + (ii) < (wid)) { \
   abucket *a = (acc) + ( (i) + (ii) ) + ( (j) + (jj) ) * (wid); \
   abump_no_overflow(a[0][0],(c)[0]); \
   abump_no_overflow(a[0][1],(c)[1]); \
   abump_no_overflow(a[0][2],(c)[2]); \
   abump_no_overflow(a[0][3],(c)[3]); \
   } \
} while (0)
/* single-threaded */
#define USE_LOCKS
#define bump_no_overflow(dest, delta)  do {dest += delta;} while (0)
#define render_rectangle render_rectangle_double
#define iter_thread iter_thread_double
#include "rect.c"
#ifdef HAVE_GCC_64BIT_ATOMIC_OPS
/* multi-threaded */
#undef USE_LOCKS
#undef bump_no_overflow
#undef render_rectangle
#undef iter_thread
#define bump_no_overflow(dest, delta)  double_atomic_add(&dest, delta)
#define render_rectangle render_rectangle_double_mt
#define iter_thread iter_thread_double_mt
#include "rect.c"
#else /* !HAVE_GCC_64BIT_ATOMIC_OPS */
#define render_rectangle_double_mt render_rectangle_double
#endif /* HAVE_GCC_64BIT_ATOMIC_OPS */
#undef render_rectangle
#undef iter_thread
#undef add_c_to_accum
#undef A_ACCUM_T
#undef B_ACCUM_T
#undef bucket
#undef abucket
#undef bump_no_overflow
#undef abump_no_overflow

/* 32-bit datatypes */
#define B_ACCUM_T unsigned int
#define A_ACCUM_T unsigned int
#define bucket bucket_int
#define abucket abucket_int
#define abump_no_overflow(dest, delta) do { \
   if (UINT_MAX - dest > delta) dest += delta; else dest = UINT_MAX; \
} while (0)
#define add_c_to_accum(acc,i,ii,j,jj,wid,hgt,c) do { \
   if ( (j) + (jj) >=0 && (j) + (jj) < (hgt) && (i) + (ii) >=0 && (i) + (ii) < (wid)) { \
   abucket *a = (acc) + ( (i) + (ii) ) + ( (j) + (jj) ) * (wid); \
   abump_no_overflow(a[0][0],(c)[0]); \
   abump_no_overflow(a[0][1],(c)[1]); \
   abump_no_overflow(a[0][2],(c)[2]); \
   abump_no_overflow(a[0][3],(c)[3]); \
   } \
} while (0)
/* single-threaded */
#define USE_LOCKS
#define bump_no_overflow(dest, delta) do { \
   if (UINT_MAX - dest > delta) dest += delta; else dest = UINT_MAX; \
} while (0)
#define render_rectangle render_rectangle_int
#define iter_thread iter_thread_int
#include "rect.c"
#ifdef HAVE_GCC_ATOMIC_OPS
/* multi-threaded */
#undef USE_LOCKS
#undef bump_no_overflow
#undef render_rectangle
#undef iter_thread
#define bump_no_overflow(dest, delta)  uint_atomic_add(&dest, delta)
#define render_rectangle render_rectangle_int_mt
#define iter_thread iter_thread_int_mt
#include "rect.c"
#else /* !HAVE_GCC_ATOMIC_OPS */
#define render_rectangle_int_mt render_rectangle_int
#endif /* HAVE_GCC_ATOMIC_OPS */
#undef iter_thread
#undef render_rectangle
#undef add_c_to_accum
#undef A_ACCUM_T
#undef B_ACCUM_T
#undef bucket
#undef abucket
#undef bump_no_overflow
#undef abump_no_overflow

/* experimental 32-bit datatypes (called 33) */
#define B_ACCUM_T unsigned int
#define A_ACCUM_T float
#define bucket bucket_int
#define abucket abucket_float
#define abump_no_overflow(dest, delta) do {dest += delta;} while (0)
#define add_c_to_accum(acc,i,ii,j,jj,wid,hgt,c) do { \
   if ( (j) + (jj) >=0 && (j) + (jj) < (hgt) && (i) + (ii) >=0 && (i) + (ii) < (wid)) { \
   abucket *a = (acc) + ( (i) + (ii) ) + ( (j) + (jj) ) * (wid); \
   abump_no_overflow(a[0][0],(c)[0]); \
   abump_no_overflow(a[0][1],(c)[1]); \
   abump_no_overflow(a[0][2],(c)[2]); \
   abump_no_overflow(a[0][3],(c)[3]); \
   } \
} while (0)
/* single-threaded */
#define USE_LOCKS
#define bump_no_overflow(dest, delta) do { \
   if (UINT_MAX - dest > delta) dest += delta; else dest = UINT_MAX; \
} while (0)
#define render_rectangle render_rectangle_float
#define iter_thread iter_thread_float
#include "rect.c"
#ifdef HAVE_GCC_ATOMIC_OPS
/* multi-threaded */
#undef USE_LOCKS
#undef bump_no_overflow
#undef render_rectangle
#undef iter_thread
#define bump_no_overflow(dest, delta)  uint_atomic_add(&dest, delta)
#define render_rectangle render_rectangle_float_mt
#define iter_thread iter_thread_float_mt
#include "rect.c"
#else /* !HAVE_GCC_ATOMIC_OPS */
#define render_rectangle_float_mt render_rectangle_float
#endif /* HAVE_GCC_ATOMIC_OPS */
#undef iter_thread
#undef render_rectangle
#undef add_c_to_accum
#undef A_ACCUM_T
#undef B_ACCUM_T
#undef bucket
#undef abucket
#undef bump_no_overflow
#undef abump_no_overflow


/* 16-bit datatypes */
#define B_ACCUM_T unsigned short
#define A_ACCUM_T unsigned short
#define bucket bucket_short
#define abucket abucket_short
#define MAXBUCKET (1<<14)
#define abump_no_overflow(dest, delta) do { \
   if (USHRT_MAX - dest > delta) dest += delta; else dest = USHRT_MAX; \
} while (0)
#define add_c_to_accum(acc,i,ii,j,jj,wid,hgt,c) do { \
   if ( (j) + (jj) >=0 && (j) + (jj) < (hgt) && (i) + (ii) >=0 && (i) + (ii) < (wid)) { \
   abucket *a = (acc) + ( (i) + (ii) ) + ( (j) + (jj) ) * (wid); \
   abump_no_overflow(a[0][0],(c)[0]); \
   abump_no_overflow(a[0][1],(c)[1]); \
   abump_no_overflow(a[0][2],(c)[2]); \
   abump_no_overflow(a[0][3],(c)[3]); \
   } \
} while (0)
/* single-threaded */
#define USE_LOCKS
#define bump_no_overflow(dest, delta) do { \
   if (USHRT_MAX - dest > delta) dest += delta; else dest = USHRT_MAX; \
} while (0)
#define render_rectangle render_rectangle_short
#define iter_thread iter_thread_short
#include "rect.c"
#ifdef HAVE_GCC_ATOMIC_OPS
/* multi-threaded */
#undef USE_LOCKS
#undef bump_no_overflow
#undef render_rectangle
#undef iter_thread
#define bump_no_overflow(dest, delta)  ushort_atomic_add(&dest, delta)
#define render_rectangle render_rectangle_short_mt
#define iter_thread iter_thread_short_mt
#include "rect.c"
#else /* !HAVE_GCC_ATOMIC_OPS */
#define render_rectangle_short_mt render_rectangle_short
#endif /* HAVE_GCC_ATOMIC_OPS */
#undef iter_thread
#undef render_rectangle
#undef add_c_to_accum
#undef A_ACCUM_T
#undef B_ACCUM_T
#undef bucket
#undef abucket
#undef bump_no_overflow
#undef abump_no_overflow

double flam3_render_memory_required(flam3_frame *spec)
{
  flam3_genome *cps = spec->genomes;
  int real_bits = spec->bits;

  if (33 == real_bits) real_bits = 32;

  /* note 4 channels * 2 buffers cancels out 8 bits per byte */
  /* does not yet include memory for density estimation filter */

  return
    (double) cps[0].spatial_oversample * cps[0].spatial_oversample *
    (double) cps[0].width * cps[0].height * real_bits;
}

void bits_error(flam3_frame *spec) {
      fprintf(stderr, "flam3: bits must be 16, 32, 33, or 64 not %d.\n",
         spec->bits);
      exit(1);
}

void flam3_render(flam3_frame *spec, void *out,
        int out_width, int field, int nchan, int trans,
        stat_struct *stats) {
  if (spec->nthreads == 1) {
    /* single-threaded */
    switch (spec->bits) {
    case 16:
      render_rectangle_short(spec, out, out_width, field, nchan, trans, stats);
      break;
    case 32:
      render_rectangle_int(spec, out, out_width, field, nchan, trans, stats);
      break;
    case 33:
      render_rectangle_float(spec, out, out_width, field, nchan, trans, stats);
      break;
    case 64:
      render_rectangle_double(spec, out, out_width, field, nchan, trans, stats);
      break;
    default:
      bits_error(spec);
      break;
    }
  } else {
    /* multi-threaded */
    switch (spec->bits) {
    case 16:
      render_rectangle_short_mt(spec, out, out_width, field, nchan, trans, stats);
      break;
    case 32:
      render_rectangle_int_mt(spec, out, out_width, field, nchan, trans, stats);
      break;
    case 33:
      render_rectangle_float_mt(spec, out, out_width, field, nchan, trans, stats);
      break;
    case 64:
      render_rectangle_double_mt(spec, out, out_width, field, nchan, trans, stats);
      break;
    default:
      bits_error(spec);
      break;
    }
  }
}

/*
 *   filter function definitions
 * from Graphics Gems III code
 * and ImageMagick resize.c
 */

double flam3_hermite_filter(double t) {
   /* f(t) = 2|t|^3 - 3|t|^2 + 1, -1 <= t <= 1 */
   if(t < 0.0) t = -t;
   if(t < 1.0) return((2.0 * t - 3.0) * t * t + 1.0);
   return(0.0);
}

double flam3_box_filter(double t) {
   if((t > -0.5) && (t <= 0.5)) return(1.0);
   return(0.0);
}

double flam3_triangle_filter(double t) {
   if(t < 0.0) t = -t;
   if(t < 1.0) return(1.0 - t);
   return(0.0);
}

double flam3_bell_filter(double t) {
   /* box (*) box (*) box */
   if(t < 0) t = -t;
   if(t < .5) return(.75 - (t * t));
   if(t < 1.5) {
      t = (t - 1.5);
      return(.5 * (t * t));
   }
   return(0.0);
}

double flam3_b_spline_filter(double t) {

   /* box (*) box (*) box (*) box */
   double tt;

   if(t < 0) t = -t;
   if(t < 1) {
      tt = t * t;
      return((.5 * tt * t) - tt + (2.0 / 3.0));
   } else if(t < 2) {
      t = 2 - t;
      return((1.0 / 6.0) * (t * t * t));
   }
   return(0.0);
}

double flam3_sinc(double x) {
   x *= M_PI;
   if(x != 0) return(sin(x) / x);
   return(1.0);
}

double flam3_blackman_filter(double x) {
  return(0.42+0.5*cos(M_PI*x)+0.08*cos(2*M_PI*x));
}

double flam3_catrom_filter(double x) {
  if (x < -2.0)
    return(0.0);
  if (x < -1.0)
    return(0.5*(4.0+x*(8.0+x*(5.0+x))));
  if (x < 0.0)
    return(0.5*(2.0+x*x*(-5.0-3.0*x)));
  if (x < 1.0)
    return(0.5*(2.0+x*x*(-5.0+3.0*x)));
  if (x < 2.0)
    return(0.5*(4.0+x*(-8.0+x*(5.0-x))));
  return(0.0);
}

double flam3_mitchell_filter(double t) {
   double tt;

   tt = t * t;
   if(t < 0) t = -t;
   if(t < 1.0) {
      t = (((12.0 - 9.0 * flam3_mitchell_b - 6.0 * flam3_mitchell_c) * (t * tt))
         + ((-18.0 + 12.0 * flam3_mitchell_b + 6.0 * flam3_mitchell_c) * tt)
         + (6.0 - 2 * flam3_mitchell_b));
      return(t / 6.0);
   } else if(t < 2.0) {
      t = (((-1.0 * flam3_mitchell_b - 6.0 * flam3_mitchell_c) * (t * tt))
         + ((6.0 * flam3_mitchell_b + 30.0 * flam3_mitchell_c) * tt)
         + ((-12.0 * flam3_mitchell_b - 48.0 * flam3_mitchell_c) * t)
         + (8.0 * flam3_mitchell_b + 24 * flam3_mitchell_c));
      return(t / 6.0);
   }
   return(0.0);
}

double flam3_hanning_filter(double x) {
  return(0.5+0.5*cos(M_PI*x));
}

double flam3_hamming_filter(double x) {
  return(0.54+0.46*cos(M_PI*x));
}

double flam3_lanczos3_filter(double t) {
   if(t < 0) t = -t;
   if(t < 3.0) return(flam3_sinc(t) * flam3_sinc(t/3.0));
   return(0.0);
}

double flam3_lanczos2_filter(double t) {
   if(t < 0) t = -t;
   if(t < 2.0) return(flam3_sinc(t) * flam3_sinc(t/2.0));
   return(0.0);
}

double flam3_gaussian_filter(double x) {
  return(exp((-2.0*x*x))*sqrt(2.0/M_PI));
}

double flam3_quadratic_filter(double x) {
  if (x < -1.5)
    return(0.0);
  if (x < -0.5)
    return(0.5*(x+1.5)*(x+1.5));
  if (x < 0.5)
    return(0.75-x*x);
  if (x < 1.5)
    return(0.5*(x-1.5)*(x-1.5));
  return(0.0);
}

void b64decode(char* instr, char *outstr)
{
    char *cur, *start;
    int d, dlast, phase;
    char c;
    static int table[256] = {
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* 00-0F */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* 10-1F */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,62,-1,-1,-1,63,  /* 20-2F */
        52,53,54,55,56,57,58,59,60,61,-1,-1,-1,-1,-1,-1,  /* 30-3F */
        -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,  /* 40-4F */
        15,16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,-1,  /* 50-5F */
        -1,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,  /* 60-6F */
        41,42,43,44,45,46,47,48,49,50,51,-1,-1,-1,-1,-1,  /* 70-7F */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* 80-8F */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* 90-9F */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* A0-AF */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* B0-BF */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* C0-CF */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* D0-DF */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  /* E0-EF */
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1   /* F0-FF */
    };

    d = dlast = phase = 0;
    start = instr;
    for (cur = instr; *cur != '\0'; ++cur )
    {
   // jer: this is my bit that treats line endings as physical breaks
   if(*cur == '\n' || *cur == '\r'){phase = dlast = 0; continue;}
        d = table[(int)*cur];
        if(d != -1)
        {
            switch(phase)
            {
            case 0:
                ++phase;
                break;
            case 1:
                c = ((dlast << 2) | ((d & 0x30) >> 4));
                *outstr++ = c;
                ++phase;
                break;
            case 2:
                c = (((dlast & 0xf) << 4) | ((d & 0x3c) >> 2));
                *outstr++ = c;
                ++phase;
                break;
            case 3:
                c = (((dlast & 0x03 ) << 6) | d);
                *outstr++ = c;
                phase = 0;
                break;
            }
            dlast = d;
        }
    }
}

/* Helper functions for After Effects Plugin */
/* Function to calculate the size of a 'flattened' genome (required by AE API) */
size_t flam3_size_flattened_genome(flam3_genome *cp) {

   size_t flatsize;
   
   flatsize = sizeof(flam3_genome);
   flatsize += cp->num_xforms * sizeof(flam3_xform);
   
   return(flatsize);
}

/* Function to flatten the contents of a genome into a buffer */
void flam3_flatten_genome(flam3_genome *cp, void *buf) {

   int i;
   char *bufoff;

   /* Copy genome first */
   memcpy(buf, (const void *)cp, sizeof(flam3_genome));
   
   /* Copy the xforms */
   bufoff = (char *)buf + sizeof(flam3_genome);
   for (i=0; i<cp->num_xforms; i++) {
      memcpy(bufoff, (const void *)(&cp->xform[i]), sizeof(flam3_xform));
      bufoff += sizeof(flam3_xform);
   }
}

/* Function to unflatten a genome buffer */
void flam3_unflatten_genome(void *buf, flam3_genome *cp) {

   int i;
   char *bufoff;
   
   /* Copy out the genome */
   memcpy((void *)cp, (const void *)buf, sizeof(flam3_genome));
   
   /* Allocate space for the xforms */
   cp->xform = (flam3_xform *)malloc(cp->num_xforms * sizeof(flam3_xform));
   
   /* Initialize the xforms (good habit to be in) */
   initialize_xforms(cp, 0);
   
   /* Copy out the xforms from the buffer */
   bufoff = (char *)buf + sizeof(flam3_genome);
   for (i=0; i<cp->num_xforms; i++) {
      memcpy(bufoff, (const void *)(&cp->xform[i]), sizeof(flam3_xform));
      bufoff += sizeof(flam3_xform);
   }
}

void flam3_srandom() {
   unsigned int seed;
   char *s = getenv("seed");

   if (s)
      seed = atoi(s);
   else
      seed = time(0) + getpid();

   srandom(seed);
}

double flam3_calc_alpha(double density, double gamma, double linrange) {

   double dnorm = density;
   double funcval = pow(linrange, gamma);
   double frac,alpha;
   
   if (dnorm>0) {
      if (dnorm < linrange) {
         frac = dnorm/linrange;
         alpha = (1.0-frac) * dnorm * (funcval / linrange) + frac * pow(dnorm,gamma);
      } else
         alpha = pow(dnorm,gamma);
   } else
      alpha = 0;
      
   return(alpha);
}

void flam3_calc_newrgb(double *cbuf, double ls, double highpow, double *newrgb) {

   int rgbi;
   double newls,lsratio;
   double newhsv[3];
   double a, maxa=-1.0, maxc;
   
   /* Identify the most saturated channel */
   for (rgbi=0;rgbi<3;rgbi++) {
      a = ls * (cbuf[rgbi]/PREFILTER_WHITE);
      if (a>maxa) {
         maxa = a;
         maxc = cbuf[rgbi]/PREFILTER_WHITE;
      }
   }
   
   /* If a channel is saturated and we have a non-negative highlight power */
   /* modify the color to prevent hue shift                                */
   if (maxa>255 && highpow>=0.0) {
      newls = 255.0/maxc;
      lsratio = pow(newls/ls,highpow);

      /* Calculate the max-value color (ranged 0 - 1) */
      for (rgbi=0;rgbi<3;rgbi++)
         newrgb[rgbi] = newls*(cbuf[rgbi]/PREFILTER_WHITE)/255.0;

      /* Reduce saturation by the lsratio */
      rgb2hsv(newrgb,newhsv);
      newhsv[1] *= lsratio;
      hsv2rgb(newhsv,newrgb);

      for (rgbi=0;rgbi<3;rgbi++)
         newrgb[rgbi] *= 255.0;
      
   } else {
      for (rgbi=0;rgbi<3;rgbi++)
         newrgb[rgbi] = ls*(cbuf[rgbi]/PREFILTER_WHITE);
   }
}
