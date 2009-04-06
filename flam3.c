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

#define _GNU_SOURCE

#include "private.h"
#include "img.h"
#include "config.h"
#include "variations.h"
#include "interpolation.h"
#include "parser.h"
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
#define DE_THRESH   50

#define random_distrib(v) ((v)[random()%vlen(v)])

int flam3_create_spatial_filter(flam3_frame *spec, int field, double **filter);
double flam3_create_temporal_filter(int numsteps, int filter_type, double filter_exp, double filter_width,
                                    double **temporal_filter, double **temporal_deltas);
flam3_de_helper flam3_create_de_filters(double max_rad, double min_rad, double curve, int ss);
                                    
                                    
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
      
      if (apply_xform(cp, fn, p, q, rc)>0) {
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

      /* Store the last used transform */
      lastxf = fn+1;

      p[0] = q[0];
      p[1] = q[1];
      p[2] = q[2];
      p[3] = q[3];

      if (cp->final_xform_enable == 1)
         apply_xform(cp, cp->final_xform_index, p, q, rc);

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

#define APPMOT(x)  do { addto->x += mot[i].x * motion_funcs(func,freq*blend); } while (0);

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
   }
   
   clear_cp(&(cpi[0]),0);
   clear_cp(&(cpi[1]),0);

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
   
/* Motion support functions */
void flam3_delete_motion_elements(flam3_xform *xf) {

   /* Free the motion elements */
   if (xf->num_motion>0) {
      free(xf->motion);
      xf->num_motion = 0;
   }
   
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
   
   int i,ii,j;
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
   for (i=0;i<numstd;i++) {
      dest->xform[i] = src->xform[i];
      dest->xform[i].num_motion=0;
      dest->xform[i].motion=NULL;
      if (src->xform[i].num_motion>0) {
         for (j=0;j<src->xform[i].num_motion;j++)
            flam3_add_motion_element(&dest->xform[i]);
         memcpy(dest->xform[i].motion,src->xform[i].motion,src->xform[i].num_motion*sizeof(flam3_xform));
      }
      
   }
      
   /* Add the final x if it's present */
   if (src->final_xform_index>=0) {
      i = src->final_xform_index;
      flam3_add_xforms(dest, 1, 0, 1);
      ii = dest->final_xform_index;
      dest->xform[ii] = src->xform[i];
      dest->xform[ii].num_motion = 0;
      dest->xform[ii].motion=NULL;
      
      if (src->xform[i].num_motion>0) {
         for (j=0;j<src->xform[i].num_motion;j++)
            flam3_add_motion_element(&dest->xform[ii]);
         memcpy(dest->xform[ii].motion,src->xform[i].motion,src->xform[i].num_motion*sizeof(flam3_xform));
      }
      
   }
   
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
      dest->xform[i].num_motion=0;
      dest->xform[i].motion=NULL;
      if (src->xform[i].num_motion>0) {
         for (j=0;j<src->xform[i].num_motion;j++)
            flam3_add_motion_element(&dest->xform[i]);
         memcpy(dest->xform[i].motion,src->xform[i].motion,src->xform[i].num_motion*sizeof(flam3_xform));
      }
      /* Copy the initial chaos from the src - the rest are already 1 */
      memcpy(dest->chaos[i], src->chaos[i], numsrcstd*sizeof(double));
      
   }   
   
   /* Add the final xform if necessary */
   if (dest_final_xform > 0) {
      flam3_add_xforms(dest, dest_final_xform, 1, 1);

      if (src->final_xform_enable > 0) {
         i = src->final_xform_index;
         dest->xform[dest->num_xforms-1] = src->xform[i];
         dest->xform[dest->num_xforms-1].num_motion = 0;
         dest->xform[dest->num_xforms-1].motion=NULL;
         if (src->xform[i].num_motion>0) {
            for (j=0;j<src->xform[i].num_motion;j++)
               flam3_add_motion_element(&dest->xform[dest->num_xforms-1]);
            memcpy(dest->xform[dest->num_xforms-1].motion,src->xform[i].motion,src->xform[i].num_motion*sizeof(flam3_xform));
         }
         
      } else {
         /* Interpolated-against final xforms need animate & color_speed set to 1.0 */
         dest->xform[dest->num_xforms-1].num_motion = 0;
         dest->xform[dest->num_xforms-1].motion=NULL;
         dest->xform[dest->num_xforms-1].animate=1.0;
         dest->xform[dest->num_xforms-1].color_speed=1.0;
      }

   } else {
      dest->final_xform_index = -1;
      dest->final_xform_enable = 0;
   }

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
    cp->edits = NULL; //BAD
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

      for (i=0;i<cp->num_xforms;i++)
         flam3_delete_motion_elements(&cp->xform[i]);
                       
       free(cp->xform);
       cp->xform=NULL;
       
       cp->num_xforms = 0;
    }
    
//    if (NULL != cp->edits)
//       xmlFreeDoc(cp->edits);
    

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

flam3_genome *flam3_parse_xml2(char *xmldata, char *xmlfilename, int default_flag, int *ncps) {

   xmlDocPtr doc; /* Parsed XML document tree */
   xmlNode *rootnode;
   char *bn;
   int i;
   int loc_all_ncps=0;
   flam3_genome *loc_all_cp=NULL;
   

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
   bn = basename(xmlfilename);
   
   /* Have to use &loc_all_cp since the memory gets allocated in scan_for_flame_nodes */
   scan_for_flame_nodes(rootnode, bn, default_flag,&loc_all_cp,&loc_all_ncps);
   
   xmlFreeDoc(doc);

   *ncps = loc_all_ncps;
   
   /* Check to see if the first control point or the second-to-last */
   /* control point has interpolation="smooth".  This is invalid    */
   /* and should be reset to linear (with a warning).               */
   if (loc_all_ncps>=1) {
      if (loc_all_cp[0].interpolation == flam3_interpolation_smooth) {
         fprintf(stderr,"Warning: smooth interpolation cannot be used for first segment.\n"
                        "         switching to linear.\n");
         loc_all_cp[0].interpolation = flam3_interpolation_linear;
      }
   }
   
   if (loc_all_ncps>=2) {
      if (loc_all_cp[(loc_all_ncps)-2].interpolation == flam3_interpolation_smooth) {
         fprintf(stderr,"Warning: smooth interpolation cannot be used for last segment.\n"
                        "         switching to linear.\n");
         loc_all_cp[loc_all_ncps-2].interpolation = flam3_interpolation_linear;
      }
   }
   
   /* Finally, ensure that consecutive 'rotate' parameters never exceed */
   /* a difference of more than 180 degrees (+/-) for interpolation.    */
   /* An adjustment of +/- 360 degrees is made until this is true.      */
   if (*ncps>1) {
   
      for (i=1;i<*ncps;i++) {

         /* Only do this adjustment if we're not in compat mode */
         if (flam3_inttype_compat != loc_all_cp[i-1].interpolation_type
        && flam3_inttype_older != loc_all_cp[i-1].interpolation_type) {
      
            while (loc_all_cp[i].rotate < loc_all_cp[i-1].rotate-180)
               loc_all_cp[i].rotate += 360;
            
            while (loc_all_cp[i].rotate > loc_all_cp[i-1].rotate+180)
               loc_all_cp[i].rotate -= 360;
         }
      }
   }
   
   xmlCleanupParser();
   
   return loc_all_cp;
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

   if (extra_attributes)
      fprintf(f, " %s", extra_attributes);

   fprintf(f, ">\n");

   if (cp->symmetry)
      fprintf(f, "   <symmetry kind=\"%d\"/>\n", cp->symmetry);
   
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

int flam3_create_spatial_filter(flam3_frame *spec, int field, double **filter) {

   int sf_kernel = spec->genomes[0].spatial_filter_select;
   int supersample = spec->genomes[0].spatial_oversample;
   double sf_radius = spec->genomes[0].spatial_filter_radius;
   double aspect_ratio = spec->pixel_aspect_ratio;   
   double sf_supp = flam3_spatial_support[sf_kernel];
   
   double fw = 2.0 * sf_supp * supersample * sf_radius / aspect_ratio;
   double adjust, ii, jj;
   
   int fwidth = ((int) fw) + 1;
   int i,j;
   
   
   /* Make sure the filter kernel has same parity as oversample */
   if ((fwidth ^ supersample) & 1)
      fwidth++;

   /* Calculate the coordinate scaling factor for the kernel values */
   if (fw > 0.0)
      adjust = sf_supp * fwidth / fw;
   else
      adjust = 1.0;

   /* Calling function MUST FREE THE RETURNED KERNEL, lest ye leak memory */
   (*filter) = (double *)malloc(sizeof(double) * fwidth * fwidth);

   /* fill in the coefs */
   for (i = 0; i < fwidth; i++)
      for (j = 0; j < fwidth; j++) {
      
         /* Calculate the function inputs for the kernel function */
         ii = ((2.0 * i + 1.0) / (double)fwidth - 1.0)*adjust;
         jj = ((2.0 * j + 1.0) / (double)fwidth - 1.0)*adjust;

         /* Scale for scanlines */
         if (field) jj *= 2.0;

         /* Adjust for aspect ratio */
         jj /= aspect_ratio;

         (*filter)[i + j * fwidth] = 
               flam3_spatial_filter(sf_kernel,ii) * flam3_spatial_filter(sf_kernel,jj);
      }


   if (normalize_vector((*filter), fwidth * fwidth)) {
      fprintf(stderr, "Spatial filter value is too small: %g.  Terminating.\n",sf_radius);
      exit(1);
   }   
   
   return (fwidth);
}

double flam3_create_temporal_filter(int numsteps, int filter_type, double filter_exp, double filter_width,
                                    double **temporal_filter, double **temporal_deltas) {

   double maxfilt = 0.0;
   double sumfilt = 0.0;
   double slpx,halfsteps;
   double *deltas, *filter;
   
   int i;

   /* Allocate memory - this must be freed in the calling routine! */   
   deltas = (double *)malloc(numsteps*sizeof(double));
   filter = (double *)malloc(numsteps*sizeof(double));
   
   /* Deal with only one step */
   if (numsteps==1) {
      deltas[0] = 0;
      filter[0] = 1.0;
      *temporal_deltas = deltas;
      *temporal_filter = filter;
      return(1.0);
   }
      

   /* Define the temporal deltas */   
   for (i = 0; i < numsteps; i++)
      deltas[i] = ((double)i /(double)(numsteps - 1) - 1.0)*filter_width;
      
   /* Define the filter coefs */
   if (flam3_temporal_exp == filter_type) {

      for (i=0; i < numsteps; i++) {

         if (filter_exp>=0)
            slpx = ((double)i+1.0)/numsteps;
         else
            slpx = (double)(numsteps - i)/numsteps;

         /* Scale the color based on these values */
         filter[i] = pow(slpx,fabs(filter_exp));
         
         /* Keep the max */
         if (filter[i]>maxfilt)
            maxfilt = filter[i];
      }

   } else if (flam3_temporal_gaussian == filter_type) {

      halfsteps = numsteps/2.0;
      for (i=0; i < numsteps; i++) {
      
         /* Gaussian */
         filter[i] = flam3_spatial_filter(flam3_gaussian_kernel,
                           flam3_spatial_support[flam3_gaussian_kernel]*fabs(i - halfsteps)/halfsteps);
         /* Keep the max */
         if (filter[i]>maxfilt)
            maxfilt = filter[i];
      }
      
   } else { // (flam3_temporal_box)

      for (i=0; i < numsteps; i++)
         filter[i] = 1.0;
         
	   maxfilt = 1.0;
	   
   }

   /* Adjust the filter so that the max is 1.0, and */
   /* calculate the K2 scaling factor  */
   for (i=0;i<numsteps;i++) {
      filter[i] /= maxfilt;
      sumfilt += filter[i];
   }
         
   sumfilt /= numsteps;
   
   *temporal_deltas = deltas;
   *temporal_filter = filter;
   
   return(sumfilt);
}                                     

flam3_de_helper flam3_create_de_filters(double max_rad, double min_rad, double curve, int ss) {

   flam3_de_helper de;
   double comp_max_radius, comp_min_radius;
   double num_de_filters_d;
   int num_de_filters,de_max_ind;
   int de_row_size, de_half_size, de_kernel_index;
   int filtloop;
   int keep_thresh=100;

   if (curve <= 0.0) {
      fprintf(stderr,"estimator curve must be > 0\n");
      exit(1);
   }

   if (max_rad < min_rad) {
      fprintf(stderr,"estimator must be larger than estimator_minimum.\n");
      fprintf(stderr,"(%f > %f) ? \n",max_rad,min_rad);
      exit(1);
   }

   /* We should scale the filter width by the oversample          */
   /* The '+1' comes from the assumed distance to the first pixel */
   comp_max_radius = max_rad*ss + 1;
   comp_min_radius = min_rad*ss + 1;

   /* Calculate how many filter kernels we need based on the decay function */
   /*                                                                       */
   /*    num filters = (de_max_width / de_min_width)^(1/estimator_curve)    */
   /*                                                                       */
   num_de_filters_d = pow( comp_max_radius/comp_min_radius, 1.0/curve );
   num_de_filters = ceil(num_de_filters_d);
         
   /* Condense the smaller kernels to save space */
   if (num_de_filters>keep_thresh) { 
      de_max_ind = ceil(DE_THRESH + pow(num_de_filters-DE_THRESH,curve))+1;
      de.max_filtered_counts = pow( (double)(de_max_ind-DE_THRESH), 1.0/curve) + DE_THRESH;
   } else {
      de_max_ind = num_de_filters;
      de.max_filtered_counts = de_max_ind;
   }

   /* Allocate the memory for these filters */
   /* and the hit/width lookup vector       */
   de_row_size = 2*ceil(comp_max_radius)-1;
   de_half_size = (de_row_size-1)/2;
   de.kernel_size = (de_half_size+1)*(2+de_half_size)/2;

   de.filter_coefs = (double *) malloc (de_max_ind * de.kernel_size * sizeof(double));
   de.filter_widths = (double *) malloc (de_max_ind * sizeof(double));

   /* Generate the filter coefficients */
   de.max_filter_index = 0;
   for (filtloop=0;filtloop<de_max_ind;filtloop++) {

      double de_filt_sum=0.0, de_filt_d;
      double de_filt_h;
      int dej,dek;
      double adjloop;
      double coef;
      int filter_coef_idx;

      /* Calculate the filter width for this number of hits in a bin */
      if (filtloop<keep_thresh)
         de_filt_h = (comp_max_radius / pow(filtloop+1,curve));
      else {
         adjloop = pow(filtloop-keep_thresh,(1.0/curve)) + keep_thresh;
         de_filt_h = (comp_max_radius / pow(adjloop+1,curve));
      }

      /* Once we've reached the min radius, don't populate any more */
      if (de_filt_h <= comp_min_radius) {
         de_filt_h = comp_min_radius;
         de.max_filter_index = filtloop;
      }

      de.filter_widths[filtloop] = de_filt_h;

      /* Calculate norm of kernel separately (easier) */
      for (dej=-de_half_size; dej<=de_half_size; dej++) {
         for (dek=-de_half_size; dek<=de_half_size; dek++) {
            
            de_filt_d = sqrt( (double)(dej*dej+dek*dek) ) / de_filt_h;

            /* Only populate the coefs within this radius */
            if (de_filt_d <= 1.0) {

               /* Gaussian */
               de_filt_sum += flam3_spatial_filter(flam3_gaussian_kernel,
                        flam3_spatial_support[flam3_gaussian_kernel]*de_filt_d);

               /* Epanichnikov */
//             de_filt_sum += (1.0 - (de_filt_d * de_filt_d));
            }
         }
      }

      filter_coef_idx = filtloop*de.kernel_size;

      /* Calculate the unique entries of the kernel */
      for (dej=0; dej<=de_half_size; dej++) {
         for (dek=0; dek<=dej; dek++) {
            de_filt_d = sqrt( (double)(dej*dej+dek*dek) ) / de_filt_h;

            /* Only populate the coefs within this radius */
            if (de_filt_d>1.0)
               de.filter_coefs[filter_coef_idx] = 0.0;
            else {

               /* Gaussian */
               de.filter_coefs[filter_coef_idx] = flam3_spatial_filter(flam3_gaussian_kernel,
                        flam3_spatial_support[flam3_gaussian_kernel]*de_filt_d)/de_filt_sum; 
                            
               /* Epanichnikov */
//             de_filter_coefs[filter_coef_idx] = (1.0 - (de_filt_d * de_filt_d))/de_filt_sum;
            }
                  
            filter_coef_idx ++;
         }
      }

      if (de.max_filter_index>0)
         break;
   }

   if (de.max_filter_index==0)
      de.max_filter_index = de_max_ind-1;

   
   return(de);
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
