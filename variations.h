/*
    FLAM3 - cosmic recursive fractal flames
    Copyright (C) 1992-2009 Spotworks LLC

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef variations_included
#define variations_included

#include "private.h"



/* Variation functions */
void var0_linear(void *, double);
void var1_sinusoidal(void *, double);
void var2_spherical(void *, double);
void var3_swirl(void *, double);
void var4_horseshoe(void *, double);
void var5_polar(void *, double);
void var6_handkerchief(void *, double);
void var7_heart(void *, double);
void var8_disc(void *, double);
void var9_spiral(void *, double);
void var10_hyperbolic(void *, double);
void var11_diamond(void *, double);
void var12_ex(void *, double);
void var13_julia(void *, double);
void var14_bent(void *, double);
void var15_waves(void *, double);
void var16_fisheye(void *, double);
void var17_popcorn(void *, double);
void var18_exponential(void *, double);
void var19_power(void *, double);
void var20_cosine(void *, double);
void var21_rings(void *, double);
void var22_fan(void *, double);
void var23_blob(void *, double);
void var24_pdj(void *, double);
void var25_fan2(void *, double);
void var26_rings2(void *, double);
void var27_eyefish(void *, double);
void var28_bubble(void *, double);
void var29_cylinder(void *, double);
void var30_perspective(void *, double);
void var31_noise(void *, double);
void var32_juliaN_generic(void *, double);
void var33_juliaScope_generic(void *, double);
void var34_blur(void *, double);
void var35_gaussian(void *, double);
void var36_radial_blur(void *, double);
void var37_pie(void *, double);
void var38_ngon(void *, double);
void var39_curl(void *, double);
void var40_rectangles(void *, double);
void var41_arch(void *helper, double weight);
void var42_tangent(void *helper, double weight);
void var43_square(void *helper, double weight);
void var44_rays(void *helper, double weight);
void var45_blade(void *helper, double weight);
void var46_secant2(void *helper, double weight);
void var47_twintrian(void *helper, double weight);
void var48_cross(void *helper, double weight);
void var49_disc2(void *helper, double weight);
void var50_supershape(void *helper, double weight);
void var51_flower(void *helper, double weight);
void var52_conic(void *helper, double weight);
void var53_parabola(void *helper, double weight);
void var54_bent2(void *helper, double weight);
void var55_bipolar(void *helper, double weight);
void var56_boarders(void *helper, double weight);
void var57_butterfly(void *helper, double weight);
void var58_cell(void *helper, double weight);
void var59_cpow(void *helper, double weight);
void var60_curve(void *helper, double weight);
void var61_edisc(void *helper, double weight);
void var62_elliptic(void *helper, double weight);
void var63_escher(void *helper, double weight);
void var64_foci(void *helper, double weight);
void var65_lazysusan(void *helper, double weight);
void var66_loonie(void *helper, double weight);
void var67_pre_blur(void *helper, double weight);
void var68_modulus(void *helper, double weight);
void var69_oscope(void *helper, double weight);
void var70_polar2(void *helper, double weight);
void var71_popcorn2(void *helper, double weight);
void var72_scry(void *helper, double weight);
void var73_separation(void *helper, double weight);
void var74_split(void *helper, double weight);
void var75_splits(void *helper, double weight);
void var76_stripes(void *helper, double weight);
void var77_wedge(void *helper, double weight);
void var78_wedge_julia(void *helper, double weight);
void var79_wedge_sph(void *helper, double weight);
void var80_whorl(void *helper, double weight);
void var81_waves2(void *helper, double weight);

/* Precalculation functions */
void perspective_precalc(flam3_xform *xf);
void juliaN_precalc(flam3_xform *xf);
void juliaScope_precalc(flam3_xform *xf);
void radial_blur_precalc(flam3_xform *xf);
void waves_precalc(flam3_xform *xf);
void disc2_precalc(flam3_xform *xf);
void supershape_precalc(flam3_xform *xf);
void wedgeJulia_precalc(flam3_xform *xf);

void xform_precalc(flam3_genome *cp, int xi);
void prepare_xform_fn_ptrs(flam3_genome *, randctx *);

int apply_xform(flam3_genome *cp, int fn, double *p, double *q, randctx *rc);
void initialize_xforms(flam3_genome *thiscp, int start_here);
#endif
