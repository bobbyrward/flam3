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

#include "interpolation.h"

double adjust_percentage(double in) {

   if (in==0.0)
      return(0.0);
   else
      return(pow(10.0, -log(1.0/in)/log(2)));

}

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

double det_matrix(double s[2][2]) {
   return s[0][0] * s[1][1] - s[0][1] * s[1][0];
}

int id_matrix(double s[3][2]) {
   return
      (s[0][0] == 1.0) &&
      (s[0][1] == 0.0) &&
      (s[1][0] == 0.0) &&
      (s[1][1] == 1.0) &&
      (s[2][0] == 0.0) &&
      (s[2][1] == 0.0);
}

void copy_matrix(double to[3][2], double from[3][2]) {

   to[0][0] = from[0][0];
   to[0][1] = from[0][1];
   to[1][0] = from[1][0];
   to[1][1] = from[1][1];
   to[2][0] = from[2][0];
   to[2][1] = from[2][1];
}


void clear_matrix(double m[3][2]) {
   m[0][0] = 0.0;
   m[0][1] = 0.0;
   m[1][0] = 0.0;
   m[1][1] = 0.0;
   m[2][0] = 0.0;
   m[2][1] = 0.0;
}

void sum_matrix(double s, double m1[3][2], double m2[3][2]) {

   m2[0][0] += s * m1[0][0];
   m2[0][1] += s * m1[0][1];
   m2[1][0] += s * m1[1][0];
   m2[1][1] += s * m1[1][1];
   m2[2][0] += s * m1[2][0];
   m2[2][1] += s * m1[2][1];
}

void mult_matrix(double s1[2][2], double s2[2][2], double d[2][2]) {
   d[0][0] = s1[0][0] * s2[0][0] + s1[1][0] * s2[0][1];
   d[1][0] = s1[0][0] * s2[1][0] + s1[1][0] * s2[1][1];
   d[0][1] = s1[0][1] * s2[0][0] + s1[1][1] * s2[0][1];
   d[1][1] = s1[0][1] * s2[1][0] + s1[1][1] * s2[1][1];
}

int compare_xforms(const void *av, const void *bv) {
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

void interpolate_cmap(flam3_palette cmap, double blend,
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

void interp_and_convert_back(double *c, int ncps, int xfi, double cxang[4][2], 
                             double cxmag[4][2], double cxtrn[4][2],double store_array[3][2]) {

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

void convert_linear_to_polar(flam3_genome *cp, int ncps, int xfi, int cflag, 
                             double cxang[4][2], double cxmag[4][2], double cxtrn[4][2]) {

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

void interpolate_catmull_rom(flam3_genome cps[], double t, flam3_genome *result) {
   double t2 = t * t;
   double t3 = t2 * t;
   double cmc[4];

   cmc[0] = (2*t2 - t - t3) / 2;
   cmc[1] = (3*t3 - 5*t2 + 2) / 2;
   cmc[2] = (4*t2 - 3*t3 + t) / 2;
   cmc[3] = (t3 - t2) / 2;

   flam3_interpolate_n(result, 4, cps, cmc);
}


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

   /* Interpolate each xform */
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
