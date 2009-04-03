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

/* this file is included into flam3.c once for each buffer bit-width */

/*
 * for batch
 *   generate de filters
 *   for temporal_sample_batch
 *     interpolate
 *     compute colormap
 *     for subbatch
 *       compute samples
 *       buckets += cmap[samples]
 *   accum += time_filter[temporal_sample_batch] * log[buckets] * de_filter
 * image = filter(accum)
 */


/* allow this many iterations for settling into attractor */
#define FUSE 15

#define PREFILTER_WHITE 255
#define WHITE_LEVEL 255


#if defined(HAVE_LIBPTHREAD) && defined(USE_LOCKS)
  /* mutex for bucket accumulator */
  pthread_mutex_t bucket_mutex;
#endif

static void iter_thread(void *fth) {
   double sub_batch;
   int j;
   flam3_thread_helper *fthp = (flam3_thread_helper *)fth;
   flam3_iter_constants *ficp = fthp->fic;
   /* Timer information needs to be persistent - only updated by 1 thread*/
   static time_t progress_timer;
   static time_t progress_timer_history[64];
   static double progress_history[64];
   static int progress_history_mark = 0;
   double eta = 0.0;
   flam3_palette cmap_scaled;
   
   if (fthp->timer_initialize) {
   	progress_timer = 0;
   	memset(progress_timer_history,0,64*sizeof(time_t));
   	memset(progress_history,0,64*sizeof(double));
   	progress_history_mark = 0;
   }
   
   for (sub_batch = 0; sub_batch < ficp->batch_size; sub_batch+=SUB_BATCH_SIZE) {
      int sub_batch_size, badcount;
      time_t newt = time(NULL);
      /* sub_batch is double so this is sketchy */
      sub_batch_size = (sub_batch + SUB_BATCH_SIZE > ficp->batch_size) ?
	(ficp->batch_size - sub_batch) : SUB_BATCH_SIZE;

      if (fthp->first_thread && newt != progress_timer) {
         double percent = 100.0 *
             ((((sub_batch / (double) ficp->batch_size) + ficp->temporal_sample_num)
             / ficp->ntemporal_samples) + ficp->batch_num)/ficp->nbatches;
         int old_mark = 0;
         int ticker;

	 if (ficp->spec->verbose)
	     fprintf(stderr, "\rchaos: %5.1f%%", percent);
         progress_timer = newt;
         if (progress_timer_history[progress_history_mark] &&
                progress_history[progress_history_mark] < percent)
            old_mark = progress_history_mark;

         if (percent > 0) {
            eta = (100 - percent) * (progress_timer - progress_timer_history[old_mark])
                  / (percent - progress_history[old_mark]);

	    if (ficp->spec->verbose) {
		ticker = (progress_timer&1)?':':'.';
		if (eta < 1000)
		    ticker = ':';
		if (eta > 100)
		    fprintf(stderr, "  ETA%c %.1f minutes", ticker, eta / 60);
		else
		    fprintf(stderr, "  ETA%c %ld seconds ", ticker, (long) ceil(eta));
		fprintf(stderr, "              \r");
		fflush(stderr);
	    }
         }

         progress_timer_history[progress_history_mark] = progress_timer;
         progress_history[progress_history_mark] = percent;
         progress_history_mark = (progress_history_mark + 1) % 64;
      }

      if (ficp->spec->progress) {
	if (fthp->first_thread) {
	  if ((*ficp->spec->progress)(ficp->spec->progress_parameter,
				      sub_batch/(double)ficp->batch_size, 0, eta)) {
	    ficp->aborted = 1;
            #ifdef HAVE_LIBPTHREAD
	      pthread_exit((void *)0);
            #else
              return;
            #endif
	  }
	} else {
          #ifdef HAVE_LIBPTHREAD
	    if (ficp->aborted) pthread_exit((void *)0);
          #else
            if (ficp->aborted) return;
          #endif
	}
      }

      /* Seed iterations */
      fthp->iter_storage[0] = flam3_random_isaac_11(&(fthp->rc));
      fthp->iter_storage[1] = flam3_random_isaac_11(&(fthp->rc));
      fthp->iter_storage[2] = flam3_random_isaac_01(&(fthp->rc));
      fthp->iter_storage[3] = flam3_random_isaac_01(&(fthp->rc));

      /* Execute iterations */
      badcount = flam3_iterate(&(fthp->cp), sub_batch_size, FUSE, fthp->iter_storage, ficp->xform_distrib, &(fthp->rc));

      #if defined(HAVE_LIBPTHREAD) && defined(USE_LOCKS)
        /* Lock mutex for access to accumulator */
        pthread_mutex_lock(&bucket_mutex);
      #endif

      /* Add the badcount to the counter */
      ficp->badvals += badcount;

      /* Put them in the bucket accumulator */
      for (j = 0; j < sub_batch_size*4; j+=4) {
         double p0, p1, p00, p11;
         double dbl_index0,dbl_frac;
         double interpcolor[4];
         int k, ci, color_index0, color_index1;
         double *p = &(fthp->iter_storage[j]);
         bucket *b;

         if (fthp->cp.rotate != 0.0) {
            p00 = p[0] - fthp->cp.rot_center[0];
            p11 = p[1] - fthp->cp.rot_center[1];
            p0 = p00 * ficp->rot[0][0] + p11 * ficp->rot[0][1] + fthp->cp.rot_center[0];
            p1 = p00 * ficp->rot[1][0] + p11 * ficp->rot[1][1] + fthp->cp.rot_center[1];
         } else {
            p0 = p[0];
            p1 = p[1];
         }

         if (p0 >= ficp->bounds[0] && p1 >= ficp->bounds[1] && p0 <= ficp->bounds[2] && p1 <= ficp->bounds[3]) {

            double logvis=1.0;
            bucket *buckets = (bucket *)(ficp->buckets);

            /* Skip if invisible */
            if (p[3]==0)
               continue;
            else
               logvis = p[3];
            
            b = buckets + (int)(ficp->ws0 * p0 - ficp->wb0s0) +
                ficp->width * (int)(ficp->hs1 * p1 - ficp->hb1s1);

#ifdef USE_FLOAT_INDICES
            color_index0 = 0;
            
            //fprintf(stdout,"%.16f\n",p[2]*256.0);
            
            while(color_index0 < CMAP_SIZE_M1) {
            	if (ficp->dmap[color_index0+1].index > p[2])
            	   break;
            	else
            	   color_index0++;
            	   //fprintf(stderr,"- - %f\n",cmap_scaled[color_index0].index);
            }
            
            //fprintf(stderr,"p[2],ci0 = %f, %d\n",p[2],color_index0);

            if (p[3]==1.0) {
               bump_no_overflow(b[0][0], ficp->dmap[color_index0].color[0]);
               bump_no_overflow(b[0][1], ficp->dmap[color_index0].color[1]);
               bump_no_overflow(b[0][2], ficp->dmap[color_index0].color[2]);
               bump_no_overflow(b[0][3], ficp->dmap[color_index0].color[3]);
               bump_no_overflow(b[0][4], 255.0);
            } else {
               bump_no_overflow(b[0][0], logvis*ficp->dmap[color_index0].color[0]);
               bump_no_overflow(b[0][1], logvis*ficp->dmap[color_index0].color[1]);
               bump_no_overflow(b[0][2], logvis*ficp->dmap[color_index0].color[2]);
               bump_no_overflow(b[0][3], logvis*ficp->dmap[color_index0].color[3]);
               bump_no_overflow(b[0][4], 255.0);
#else
            dbl_index0 = p[2] * CMAP_SIZE;
            color_index0 = (int) (dbl_index0);
            
            if (flam3_palette_mode_linear == fthp->cp.palette_mode) {
               if (color_index0 < 0) {
                  color_index0 = 0;
                  dbl_frac = 0;
               } else if (color_index0 >= CMAP_SIZE_M1) {
                  color_index0 = CMAP_SIZE_M1-1;
                  dbl_frac = 1.0;
               } else {
                  /* interpolate between color_index0 and color_index0+1 */
                  dbl_frac = dbl_index0 - (double)color_index0;
               }
                        
               for (ci=0;ci<4;ci++) {
                  interpcolor[ci] = ficp->dmap[color_index0].color[ci] * (1.0-dbl_frac) + 
                                    ficp->dmap[color_index0+1].color[ci] * dbl_frac;               
               }
               
            } else { /* Palette mode step */
            
               if (color_index0 < 0) {
                  color_index0 = 0;
               } else if (color_index0 >= CMAP_SIZE_M1) {
                  color_index0 = CMAP_SIZE_M1;
               }
                        
               for (ci=0;ci<4;ci++)
                  interpcolor[ci] = ficp->dmap[color_index0].color[ci];               
            }

            if (p[3]==1.0) {
               bump_no_overflow(b[0][0], interpcolor[0]);
               bump_no_overflow(b[0][1], interpcolor[1]);
               bump_no_overflow(b[0][2], interpcolor[2]);
               bump_no_overflow(b[0][3], interpcolor[3]);
               bump_no_overflow(b[0][4], 255.0);
            //fprintf(stderr,"%10.8f %10.8f %10.8f %10.8f\n",b[0][0],b[0][1],b[0][2],b[0][3]);
            } else {
               bump_no_overflow(b[0][0], logvis*interpcolor[0]);
               bump_no_overflow(b[0][1], logvis*interpcolor[1]);
               bump_no_overflow(b[0][2], logvis*interpcolor[2]);
               bump_no_overflow(b[0][3], logvis*interpcolor[3]);
               bump_no_overflow(b[0][4], 255.0);
            }
#endif

         }
      }
      
      #if defined(HAVE_LIBPTHREAD) && defined(USE_LOCKS)
        /* Release mutex */
        pthread_mutex_unlock(&bucket_mutex);
      #endif

   }
   #ifdef HAVE_LIBPTHREAD
     pthread_exit((void *)0);
   #endif
}

static void render_rectangle(flam3_frame *spec, void *out,
			     int out_width, int field, int nchan,
			     int transp, stat_struct *stats) {
   long nbuckets;
   int i, j, k, batch_num, temporal_sample_num;
   double nsamples, batch_size;
   bucket  *buckets;
   abucket *accumulate;
   double *points;
   double *filter, *temporal_filter, *temporal_deltas, *batch_filter;
   double ppux=0, ppuy=0;
   int image_width, image_height;    /* size of the image to produce */
   int filter_width;
   int bytes_per_channel = spec->bytes_per_channel;
   int oversample = spec->genomes[0].spatial_oversample;
   double highpow = spec->genomes[0].highlight_power;
   int nbatches = spec->genomes[0].nbatches;
   /* ntemporal_samples ( & nbatches?) should be computed inside the batch
      loop the same place the DE parms are.  problem is that batches now
      depends on times.  fix by making all batches happen at the same time */
   int ntemporal_samples = spec->genomes[0].ntemporal_samples;
   flam3_palette dmap;
   int gutter_width;
   double vibrancy = 0.0;
   double gamma = 0.0;
   double background[3];
   int vib_gam_n = 0;
   double keep_thresh=100.0;
   time_t progress_timer = 0, progress_began=0;
   int verbose = spec->verbose;
   char *fname = getenv("image");
   int gnm_idx,max_gnm_de_fw,de_offset;
   flam3_genome cp;
   unsigned short *xform_distrib;
   int num_std_xf;
   char *ai;
   flam3_iter_constants fic;
   flam3_thread_helper *fth;
#ifdef HAVE_LIBPTHREAD
   pthread_attr_t pt_attr;
   pthread_t *myThreads=NULL;
#endif
   int thread_status;
   int thi;
   time_t tstart,tend;
   
   double sumfilt;

   tstart = time(NULL);

   fic.badvals = 0;
   fic.aborted = 0;

   stats->num_iters = 0;

   memset(&cp,0, sizeof(flam3_genome));

   if (nbatches < 1) {
       fprintf(stderr, "nbatches must be positive," " not %d.\n", nbatches);
       exit(1);
   }

   if (oversample < 1) {
       fprintf(stderr, "oversample must be positive," " not %d.\n", oversample);
       exit(1);
   }

   /* Initialize the thread helper structures */
   fth = (flam3_thread_helper *)calloc(spec->nthreads,sizeof(flam3_thread_helper));
   for (i=0;i<spec->nthreads;i++) {
      fth[i].cp.final_xform_index=-1;
   }
   image_width = spec->genomes[0].width;
   if (field) {
      image_height = spec->genomes[0].height / 2;
      if (field == flam3_field_odd)
    out += nchan * bytes_per_channel * out_width;
      out_width *= 2;
   } else
      image_height = spec->genomes[0].height;


   /* Spatial Filter kernel creation */


   if (1) {
   int sf_sel = spec->genomes[0].spatial_filter_select;
   double sf_supp = flam3_spatial_support[sf_sel];
   double fw =  (2.0 * sf_supp * oversample *
         spec->genomes[0].spatial_filter_radius /
         spec->pixel_aspect_ratio);
   double adjust;

   filter_width = ((int) fw) + 1;
   /* make sure it has same parity as oversample */
   if ((filter_width ^ oversample) & 1)
      filter_width++;

   if (fw > 0.0)
      adjust = sf_supp * filter_width / fw;
   else
      adjust = 1.0;

#if 0
   fprintf(stderr, "fw = %g filter_width = %d adjust=%g\n",
         fw, filter_width, adjust);
#endif

   filter = (double *) malloc(sizeof(double) * filter_width * filter_width);
   /* fill in the coefs */
   for (i = 0; i < filter_width; i++)
      for (j = 0; j < filter_width; j++) {
         double ii = ((2.0 * i + 1.0) / filter_width - 1.0)*adjust;
         double jj = ((2.0 * j + 1.0) / filter_width - 1.0)*adjust;

         if (field) jj *= 2.0;

         jj /= spec->pixel_aspect_ratio;

         filter[i + j * filter_width] = 
            flam3_spatial_filter(sf_sel,ii) * flam3_spatial_filter(sf_sel,jj);
      }


   if (normalize_vector(filter, filter_width * filter_width)) {
      fprintf(stderr, "spatial filter value is too small: %g.\n",
         spec->genomes[0].spatial_filter_radius);
      exit(1);
   }
#if 0
      printf("vvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
      for (j = 0; j < filter_width; j++) {
    for (i = 0; i < filter_width; i++) {
      printf(" %5d", (int)(10000 * filter[i + j * filter_width]));
    }
    printf("\n");
      }
      printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
      fflush(stdout);
#endif
   }

   temporal_filter = (double *) malloc(sizeof(double) * nbatches * ntemporal_samples);
   temporal_deltas = (double *) malloc(sizeof(double) * nbatches * ntemporal_samples);
   batch_filter = (double *) malloc(sizeof(double) * nbatches);
   
   for (i=0;i<nbatches;i++)
      batch_filter[i] = 1.0/(double)nbatches;

   /* batch filter */
   for (i=0; i<nbatches; i++)
      batch_filter[i]=1.0/nbatches;

   if (nbatches*ntemporal_samples > 1) {

      double nn = nbatches*ntemporal_samples;
      double maxfilt = 0.0;
      sumfilt=0.0;

      /* Define temporal deltas from center time */
      for (i = 0; i < nn; i++)
            temporal_deltas[i] = (2.0 * ((double) i / (nn - 1)) - 1.0)
                                  * (spec->genomes[0].temporal_filter_width/2.0);

      /* fill in the coefs */
      if (flam3_temporal_exp == spec->genomes[0].temporal_filter_type) {
         for (i = 0; i < nn; i++) {
            double slpx;

            if (spec->genomes[0].temporal_filter_exp>=0)
               slpx = ((double)i+1.0)/nn;
            else
               slpx = (double)(nn - i)/nn;

            /* Scale the color based on these values */
            temporal_filter[i] = pow(slpx,fabs(spec->genomes[0].temporal_filter_exp));
            if (temporal_filter[i]>maxfilt) maxfilt = temporal_filter[i];
         }

      } else if (flam3_temporal_gaussian == spec->genomes[0].temporal_filter_type) {

         double nn2 = nn/2.0;
         for (i = 0; i < nn; i++) {
            temporal_filter[i] = flam3_spatial_filter(flam3_gaussian_kernel,
                            flam3_spatial_support[flam3_gaussian_kernel]*fabs(i - nn2)/nn2);
            if (temporal_filter[i]>maxfilt) maxfilt = temporal_filter[i];
         }

      } else { // (flam3_temporal_box == spec->genomes[0].temporal_filter_type)
         for (i=0; i<nn; i++) {
            temporal_filter[i] = 1.0;
	      }
	      maxfilt = 1.0;
      }

      /* Adjust the filter so that the max is 1.0, and */
      /* calculate the K2 scaling factor  */
      for (i=0;i<nn;i++) {
         temporal_filter[i] /= maxfilt;
         sumfilt += temporal_filter[i];
      }
         
      sumfilt = sumfilt / nn;

   } else {
      sumfilt = 1.0;
      temporal_filter[0] = 1.0;
      temporal_deltas[0] = 0.0;
   }

   /* the number of additional rows of buckets we put at the edge so
      that the filter doesn't go off the edge */

   gutter_width = (filter_width - oversample) / 2;

    /* Check the size of the density estimation filter. */
    /* If the 'radius' of the density estimation filter is greater than the */
    /* gutter width, we have to pad with more.  Otherwise, we can use the same value. */
   max_gnm_de_fw=0;
   for (gnm_idx = 0; gnm_idx < spec->ngenomes; gnm_idx++) {

      int this_width = (int)ceil(spec->genomes[gnm_idx].estimator) * oversample;
      if (this_width > max_gnm_de_fw)
         max_gnm_de_fw = this_width;
   }

   /* Add OS-1 for the averaging at the edges, if it's > 0 already */
   if (max_gnm_de_fw>0)
      max_gnm_de_fw += (oversample-1);

   /* max_gnm_de_fw is now the number of pixels of additional gutter      */
   /* necessary to appropriately perform the density estimation filtering */
   /* Check to see if it's greater than the gutter_width                  */

   if (max_gnm_de_fw > gutter_width) {
      de_offset = max_gnm_de_fw - gutter_width;
      gutter_width = max_gnm_de_fw;
   } else
      de_offset = 0;

   fic.height = oversample * image_height + 2 * gutter_width;
   fic.width  = oversample * image_width  + 2 * gutter_width;

   nbuckets = (long)fic.width * (long)fic.height;
   if (1) {

     char *last_block = NULL;
     size_t memory_rqd = (sizeof(bucket) * nbuckets +
             sizeof(abucket) * nbuckets +
             4 * sizeof(double) * (size_t)SUB_BATCH_SIZE * spec->nthreads);
     last_block = (char *) malloc(memory_rqd);
     if (NULL == last_block) {
	 fprintf(stderr, "render_rectangle: cannot malloc %g bytes.\n", (double)memory_rqd);
       fprintf(stderr, "render_rectangle: h=%d w=%d nb=%ld.\n", fic.width, fic.height, nbuckets);
       exit(1);
     }
     /* else fprintf(stderr, "render_rectangle: mallocked %dMb.\n", Mb); */

     buckets = (bucket *) last_block;
     accumulate = (abucket *) (last_block + sizeof(bucket) * nbuckets);
     points = (double *)  (last_block + (sizeof(bucket) + sizeof(abucket)) * nbuckets);
   }

   if (verbose) {
     fprintf(stderr, "chaos: ");
     progress_began = time(NULL);
   }

   background[0] = background[1] = background[2] = 0.0;
   memset((char *) accumulate, 0, sizeof(abucket) * nbuckets);

   for (batch_num = 0; batch_num < nbatches; batch_num++) {
      double de_time;
      double sample_density=0.0;
      int de_row_size, de_kernel_index=0, de_half_size;
      int de_cutoff_val=0, de_count_limit;
      double *de_filter_coefs=NULL,*de_filter_widths=NULL;
      double num_de_filters_d;
      int num_de_filters=0,filtloop,de_max_ind;
      double comp_max_radius,comp_min_radius;
      double k1, area, k2;

      de_time = spec->time + temporal_deltas[batch_num*ntemporal_samples];

      memset((char *) buckets, 0, sizeof(bucket) * nbuckets);

      /* interpolate and get a control point                      */
      /* ONLY FOR DENSITY FILTER WIDTH PURPOSES                   */
      /* additional interpolation will be done in the temporal_sample loop */
      flam3_interpolate(spec->genomes, spec->ngenomes, de_time, &cp);

      /* if instructed to by the genome, create the density estimation */
      /* filter kernels.  Check boundary conditions as well.           */
      if (cp.estimator < 0.0 || cp.estimator_minimum < 0.0) {
         fprintf(stderr,"density estimator filter widths must be >= 0\n");
         exit(1);
      }

      if (spec->bits <= 32) {
     if (cp.estimator > 0.0) {
         fprintf(stderr, "warning: density estimation disabled with %d bit buffers.\n", spec->bits);
         cp.estimator = 0.0;
     }
      }

      if (cp.estimator > 0.0) {


    if (cp.estimator_curve <= 0.0) {
       fprintf(stderr,"estimator curve must be > 0\n");
       exit(1);
    }

         if (cp.estimator < cp.estimator_minimum) {
            fprintf(stderr,"estimator must be larger than estimator_minimum.\n");
            fprintf(stderr,"(%f > %f) ? \n",cp.estimator,cp.estimator_minimum);
            exit(1);
         }

         /* We should scale the filter width by the oversample          */
         /* The '+1' comes from the assumed distance to the first pixel */
         comp_max_radius = cp.estimator*oversample + 1;
         comp_min_radius = cp.estimator_minimum*oversample + 1;

         /* Calculate how many filter kernels we need based on the decay function */
         /*                                                                       */
         /*    num filters = (de_max_width / de_min_width)^(1/estimator_curve)    */
         /*                                                                       */
         num_de_filters_d = pow( comp_max_radius/comp_min_radius, (1.0/cp.estimator_curve) );
         num_de_filters = ceil(num_de_filters_d);
         
         /* Condense the smaller kernels to save space */
         if (num_de_filters>keep_thresh) { 
            de_max_ind = ceil(keep_thresh + pow(num_de_filters-keep_thresh,cp.estimator_curve))+1;
            de_count_limit = pow( (double)(de_max_ind-100), 1.0/cp.estimator_curve) + 100;
         } else {
            de_max_ind = num_de_filters;
            de_count_limit = de_max_ind;
         }

         /* Allocate the memory for these filters */
         /* and the hit/width lookup vector       */
         de_row_size = 2*ceil(comp_max_radius)-1;
         de_half_size = (de_row_size-1)/2;
         de_kernel_index = (de_half_size+1)*(2+de_half_size)/2;

         de_filter_coefs = (double *) malloc (de_max_ind * de_kernel_index * sizeof(double));
         de_filter_widths = (double *) malloc (de_max_ind * sizeof(double));

         /* Generate the filter coefficients */
         de_cutoff_val = 0;
         for (filtloop=0;filtloop<de_max_ind;filtloop++) {

            double de_filt_sum=0.0, de_filt_d;
            double de_filt_h;
            double dej,dek,adjloop;
            int filter_coef_idx;

            if (filtloop<keep_thresh)
               de_filt_h = (comp_max_radius / pow(filtloop+1,cp.estimator_curve));
            else {
	       adjloop = pow(filtloop-keep_thresh,(1/cp.estimator_curve))+keep_thresh;
               de_filt_h = (comp_max_radius / pow(adjloop+1,cp.estimator_curve));
            }

            if (de_filt_h <= comp_min_radius) {
               de_filt_h = comp_min_radius;
               de_cutoff_val = filtloop;
            }

            de_filter_widths[filtloop] = de_filt_h;

            /* Calculate norm of kernel separately (easier) */
            for (dej=-de_half_size; dej<=de_half_size; dej++) {
               for (dek=-de_half_size; dek<=de_half_size; dek++) {
                  de_filt_d = sqrt( (double)(dej*dej+dek*dek) ) / de_filt_h;

                  if (de_filt_d <= 1.0) {

                     /* Gaussian */
                     de_filt_sum += flam3_spatial_filter(flam3_gaussian_kernel,
                        flam3_spatial_support[flam3_gaussian_kernel]*de_filt_d);

//                     /* Epanichnikov */
//                     de_filt_sum += (1.0 - (de_filt_d * de_filt_d));
                  }
               }
            }

            filter_coef_idx = filtloop*de_kernel_index;

            /* Calculate the unique entries of the kernel */
            for (dej=0; dej<=de_half_size; dej++) {
               for (dek=0; dek<=dej; dek++) {
                  de_filt_d = sqrt( (double)(dej*dej+dek*dek) ) / de_filt_h;

                  if (de_filt_d>1.0)
                     de_filter_coefs[filter_coef_idx] = 0.0;
                  else {

                     /* Gaussian */
                     de_filter_coefs[filter_coef_idx] = flam3_spatial_filter(flam3_gaussian_kernel,
                            flam3_spatial_support[flam3_gaussian_kernel]*de_filt_d)/de_filt_sum; 
                            
//                     de_filter_coefs[filter_coef_idx] = adjust_percentage(de_filter_coefs[filter_coef_idx]);

//                     /* Epanichnikov */
//                     de_filter_coefs[filter_coef_idx] = (1.0 - (de_filt_d * de_filt_d))/de_filt_sum;
                  }
                  
                  filter_coef_idx ++;
               }
            }

            if (de_cutoff_val>0)
               break;
         }

         if (de_cutoff_val==0)
            de_cutoff_val = num_de_filters-1;

      }
      
#if 0
      printf("DE Filters: %d\n",de_max_ind);
      if (1) { int dej,dek1,dek2,fci;
      fci=0;
//      for (dej=0; dej<de_max_ind;dej++) {
      for (dej=0; dej<1;dej++) {
         printf("de_filter_widths[%d] = %f\n",dej,de_filter_widths[dej]);
         for (dek1=0;dek1<=de_half_size;dek1++) {
            for (dek2=0;dek2<=dek1;dek2++) {
               printf("%f ",de_filter_coefs[fci++]);
            }
            printf("\n");
         }
         printf("--------------------------\n");
      }
      
      }
#endif
      for (temporal_sample_num = 0; temporal_sample_num < ntemporal_samples; temporal_sample_num++) {
         double temporal_sample_time;
         double color_scalar = temporal_filter[batch_num*ntemporal_samples + temporal_sample_num];

         temporal_sample_time = spec->time +
        temporal_deltas[batch_num*ntemporal_samples + temporal_sample_num];

         /* Interpolate and get a control point */
         flam3_interpolate(spec->genomes, spec->ngenomes, temporal_sample_time, &cp);

         prepare_xform_fn_ptrs(&cp, &spec->rc);

         /* compute the colormap entries.                             */
         /* the input colormap is 256 long with entries from 0 to 1.0 */
         for (j = 0; j < CMAP_SIZE; j++) {
            dmap[j].index = cp.palette[(j * 256) / CMAP_SIZE].index / 256.0;
            for (k = 0; k < 4; k++)
               dmap[j].color[k] = (cp.palette[(j * 256) / CMAP_SIZE].color[k] * WHITE_LEVEL) * color_scalar;
         }

         /* compute camera */
         if (1) {
            double t0, t1, shift=0.0, corner0, corner1;
            double scale;

            if (cp.sample_density <= 0.0) {
              fprintf(stderr,
                 "sample density (quality) must be greater than zero,"
                 " not %g.\n", cp.sample_density);
              exit(1);
            }

            scale = pow(2.0, cp.zoom);
            sample_density = cp.sample_density * scale * scale;

            ppux = cp.pixels_per_unit * scale;
            ppuy = field ? (ppux / 2.0) : ppux;
            ppux /=  spec->pixel_aspect_ratio;
            switch (field) {
               case flam3_field_both: shift =  0.0; break;
               case flam3_field_even: shift = -0.5; break;
               case flam3_field_odd:  shift =  0.5; break;
            }
            shift = shift / ppux;
            t0 = (double) gutter_width / (oversample * ppux);
            t1 = (double) gutter_width / (oversample * ppuy);
            corner0 = cp.center[0] - image_width / ppux / 2.0;
            corner1 = cp.center[1] - image_height / ppuy / 2.0;
            fic.bounds[0] = corner0 - t0;
            fic.bounds[1] = corner1 - t1 + shift;
            fic.bounds[2] = corner0 + image_width  / ppux + t0;
            fic.bounds[3] = corner1 + image_height / ppuy + t1 + shift;
            fic.size[0] = 1.0 / (fic.bounds[2] - fic.bounds[0]);
            fic.size[1] = 1.0 / (fic.bounds[3] - fic.bounds[1]);
            fic.rot[0][0] = cos(cp.rotate * 2 * M_PI / 360.0);
            fic.rot[0][1] = -sin(cp.rotate * 2 * M_PI / 360.0);
            fic.rot[1][0] = -fic.rot[0][1];
            fic.rot[1][1] = fic.rot[0][0];
            fic.ws0 = fic.width * fic.size[0];
            fic.wb0s0 = fic.ws0 * fic.bounds[0];
            fic.hs1 = fic.height * fic.size[1];
            fic.hb1s1 = fic.hs1 * fic.bounds[1];

         }

//         nsamples = sample_density * nbuckets / (oversample * oversample);
         nsamples = sample_density * image_width * image_height;
#if 0
         fprintf(stderr, "sample_density=%g nsamples=%g nbuckets=%ld time=%g\n",
       sample_density, nsamples, nbuckets, temporal_sample_time);
#endif

         batch_size = nsamples / (nbatches * ntemporal_samples);

         stats->num_iters += batch_size;
         
         /* Set up the xform_distrib array */
         xform_distrib = flam3_create_xform_distrib(&cp);
         
         /* Fill in the iter constants */
         fic.xform_distrib = xform_distrib;
         fic.spec = spec;
         fic.batch_size = batch_size / (double)spec->nthreads;
         fic.temporal_sample_num = temporal_sample_num;
         fic.ntemporal_samples = ntemporal_samples;
         fic.batch_num = batch_num;
         fic.nbatches = nbatches;

         fic.fname_specified = 0;
         fic.dmap = (flam3_palette_entry *)dmap;
         fic.color_scalar = color_scalar;
         fic.buckets = (void *)buckets;

         /* Initialize the thread helper structures */
         for (thi = 0; thi < spec->nthreads; thi++)
         {
            int rk;
            /* Create a new isaac state for this thread */
            for (rk = 0; rk < RANDSIZ; rk++)
               fth[thi].rc.randrsl[rk] = irand(&spec->rc);

            irandinit(&(fth[thi].rc),1);

            if (0==thi) {

               fth[thi].first_thread=1;
               if (0==batch_num && 0==temporal_sample_num)
               	fth[thi].timer_initialize = 1;
               else
               	fth[thi].timer_initialize = 0;
               	
            } else {
               fth[thi].first_thread=0;
	         	fth[thi].timer_initialize = 0;
            }

            fth[thi].iter_storage = &(points[thi*SUB_BATCH_SIZE*4]);
            fth[thi].fic = &fic;
            flam3_copy(&(fth[thi].cp),&cp);

         }

#ifdef HAVE_LIBPTHREAD
         /* Let's make some threads */
         myThreads = (pthread_t *)malloc(spec->nthreads * sizeof(pthread_t));

         #if defined(USE_LOCKS)
         pthread_mutex_init(&bucket_mutex, NULL);
         #endif

         pthread_attr_init(&pt_attr);
         pthread_attr_setdetachstate(&pt_attr,PTHREAD_CREATE_JOINABLE);

         for (thi=0; thi <spec->nthreads; thi ++)
            pthread_create(&myThreads[thi], &pt_attr, (void *)iter_thread, (void *)(&(fth[thi])));

         pthread_attr_destroy(&pt_attr);

         /* Wait for them to return */
         for (thi=0; thi < spec->nthreads; thi++)
            pthread_join(myThreads[thi], (void **)&thread_status);

         #if defined(USE_LOCKS)
         pthread_mutex_destroy(&bucket_mutex);
         #endif
         
         free(myThreads);
#else
         for (thi=0; thi < spec->nthreads; thi++)
            iter_thread( (void *)(&(fth[thi])) );
#endif
         
         /* Free the xform_distrib array */
         free(xform_distrib);
             
	 if (fic.aborted) {
	   if (verbose) fprintf(stderr, "\naborted!\n");
	   goto done;
	 }

         vibrancy += cp.vibrancy;
         gamma += cp.gamma;
         background[0] += cp.background[0];
         background[1] += cp.background[1];
         background[2] += cp.background[2];
         vib_gam_n++;

      }

      k1 =(cp.contrast * cp.brightness *
      PREFILTER_WHITE * 268.0 * batch_filter[batch_num]) / 256;
      area = image_width * image_height / (ppux * ppuy);
      k2 = (oversample * oversample * nbatches) /
   (cp.contrast * area * WHITE_LEVEL * sample_density * sumfilt);
#if 0
      printf("iw=%d,ih=%d,ppux=%f,ppuy=%f\n",image_width,image_height,ppux,ppuy);
      printf("contrast=%f, brightness=%f, PREFILTER=%d, temporal_filter=%f\n",
        cp.contrast, cp.brightness, PREFILTER_WHITE, temporal_filter[batch_num]);
      printf("oversample=%d, nbatches=%d, area = %f, WHITE_LEVEL=%d, sample_density=%f\n",
        oversample, nbatches, area, WHITE_LEVEL, sample_density);
      printf("k1=%f,k2=%15.12f\n",k1,k2);
#endif

      if (num_de_filters == 0) {

         for (j = 0; j < fic.height; j++) {
            for (i = 0; i < fic.width; i++) {

               abucket *a = accumulate + i + j * fic.width;
               bucket *b = buckets + i + j * fic.width;
               double c[4], ls;

               c[0] = (double) b[0][0];
               c[1] = (double) b[0][1];
               c[2] = (double) b[0][2];
               c[3] = (double) b[0][3];
               if (0.0 == c[3])
                  continue;

               ls = (k1 * log(1.0 + c[3] * k2))/c[3];
               c[0] *= ls;
               c[1] *= ls;
               c[2] *= ls;
               c[3] *= ls;

               abump_no_overflow(a[0][0], c[0]);
               abump_no_overflow(a[0][1], c[1]);
               abump_no_overflow(a[0][2], c[2]);
               abump_no_overflow(a[0][3], c[3]);
            }
         }
      } else {
      
         int ss = (int)floor(oversample / 2.0);
         int scf = !((int)oversample & 1);
         double scfact = pow(oversample/(oversample+1.0), 2.0);
         
         /* Density estimation code */
         /* Remember, we already padded with an extra pixel at the beginning */
         for (j = oversample-1; j < fic.height-(oversample-1); j++) {
            for (i = oversample-1; i < fic.width-(oversample-1); i++) {

               int ii,jj;
               double f_select=0.0;
               int f_select_int,f_coef_idx;
               int arr_filt_width;
               bucket *b;
               double c[4],ls;

               b = buckets + i + j*fic.width;

               /* Don't do anything if there's no hits here */
               if (b[0][4] == 0 || b[0][3] == 0)
                  continue;

               /* Count density in ssxss area   */
               /* Scale if OS>1 for equal iters */
               for (ii=-ss; ii<=ss; ii++) {
                  for (jj=-ss; jj<=ss; jj++) {
                     b = buckets + (i + ii) + (j + jj)*fic.width;
                     f_select += b[0][4]/255.0;
                  }
               }
               
               if (scf)
                  f_select *= scfact;
                  
               if (f_select > de_count_limit)
                  f_select_int = de_cutoff_val;                  
               else if (f_select<=keep_thresh)
                  f_select_int = (int)ceil(f_select)-1;
               else
                  f_select_int = (int)keep_thresh +
                     (int)floor(pow(f_select-keep_thresh,cp.estimator_curve));

               /* If the filter selected below the min specified clamp it to the min */
               if (f_select_int >= de_cutoff_val)
                  f_select_int = de_cutoff_val;

               /* We only have to calculate the values for ~1/8 of the square */
               f_coef_idx = f_select_int*de_kernel_index;

               arr_filt_width = (int)ceil(de_filter_widths[f_select_int])-1;

               b = buckets + i + j*fic.width;

               for (jj=0; jj<=arr_filt_width; jj++) {
                  for (ii=0; ii<=jj; ii++) {

                     /* Skip if coef is 0 */
                     if (de_filter_coefs[f_coef_idx]==0.0) {
                        f_coef_idx++;
                        continue;
                     }
                     
                     c[0] = (double) b[0][0];
                     c[1] = (double) b[0][1];
                     c[2] = (double) b[0][2];
                     c[3] = (double) b[0][3];

                     ls = de_filter_coefs[f_coef_idx]*(k1 * log(1.0 + c[3] * k2))/c[3];

                     c[0] *= ls;
                     c[1] *= ls;
                     c[2] *= ls;
                     c[3] *= ls;

                     if (jj==0 && ii==0) {
                        add_c_to_accum(accumulate,i,ii,j,jj,fic.width,fic.height,c);
                     }
                     else if (ii==0) {
                        add_c_to_accum(accumulate,i,jj,j,0,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,-jj,j,0,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,0,j,jj,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,0,j,-jj,fic.width,fic.height,c);
                     } else if (jj==ii) {
                        add_c_to_accum(accumulate,i,ii,j,jj,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,-ii,j,jj,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,ii,j,-jj,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,-ii,j,-jj,fic.width,fic.height,c);
                     } else {
                        add_c_to_accum(accumulate,i,ii,j,jj,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,-ii,j,jj,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,ii,j,-jj,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,-ii,j,-jj,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,jj,j,ii,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,-jj,j,ii,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,jj,j,-ii,fic.width,fic.height,c);
                        add_c_to_accum(accumulate,i,-jj,j,-ii,fic.width,fic.height,c);
                     }

                     f_coef_idx++;

                  }
               }

	       if (0) {
		 abucket *a = accumulate + i + j * fic.width;
		 a[0][0] = f_select_int;
	       }
            }
       if (verbose && time(NULL) != progress_timer) {
      progress_timer = time(NULL);
      fprintf(stderr, "\rdensity estimation: %d/%d          ", j, fic.height);
      fflush(stderr);
       }


      if (fic.spec->progress) {
	if ((*fic.spec->progress)(fic.spec->progress_parameter,
				  j/(double)fic.height, 1, 0.0)) {
	   if (verbose) fprintf(stderr, "\naborted!\n");
	  goto done;
	  }
      }

         }
         
      } /* End density estimation loop */


      /* If allocated, free the de filter memory for the next batch */
      if (num_de_filters > 0) {
         free(de_filter_coefs);
         free(de_filter_widths);
      }

   }

   if (verbose) {
     fprintf(stderr, "\rchaos: 100.0%%  took: %ld seconds   \n", time(NULL) - progress_began);
     fprintf(stderr, "filtering...");
   }
   

   /* filter the accumulation buffer down into the image */
   if (1) {
      int x, y;
      double t[4],newrgb[3],newhsv[3];
      double g = 1.0 / (gamma / vib_gam_n);
      double tmp,a;
      double alpha,ls;
      int rgbi;

      double linrange = cp.gam_lin_thresh;
      double funcval = pow(linrange,g);
      double frac;
      
      int eclip = argi("earlyclip",0);

      vibrancy /= vib_gam_n;
      background[0] /= vib_gam_n/256.0;
      background[1] /= vib_gam_n/256.0;
      background[2] /= vib_gam_n/256.0;
      
      /* If we're in the early clip mode, perform this first step to  */
      /* apply the gamma correction and clipping before the spat filt */
      
      if (eclip) {

         for (j = 0; j < fic.height; j++) {
            for (i = 0; i < fic.width; i++) {

               abucket *ac = accumulate + i + j*fic.width;
               
               if (ac[0][3]<=0) {
                  alpha = 0.0;
                  ls = 0.0;
               } else {
                  tmp=ac[0][3]/PREFILTER_WHITE;
                  alpha = flam3_calc_alpha(tmp,g,linrange);
                  ls = vibrancy * 256.0 * alpha / tmp;
                  if (alpha<0.0) alpha = 0.0;
                  if (alpha>1.0) alpha = 1.0;
               }
            
               t[0] = (double)ac[0][0];
               t[1] = (double)ac[0][1];
               t[2] = (double)ac[0][2];
               t[3] = (double)ac[0][3];
            
               flam3_calc_newrgb(t, ls, highpow, newrgb);
                  
               for (rgbi=0;rgbi<3;rgbi++) {
                  a = newrgb[rgbi];
                  a += (1.0-vibrancy) * 256.0 * pow( t[rgbi] / PREFILTER_WHITE, g);
                  if (nchan<=3 || transp==0)
                     a += ((1.0 - alpha) * background[rgbi]);
                  else {
                     if (alpha>0)
                        a /= alpha;
                     else
                        a = 0;
                  }

                  /* Clamp here to ensure proper filter functionality */
                  if (a>255) a = 255;
                  if (a<0) a = 0;
               
                  /* Replace values in accumulation buffer with these new ones */
                  ac[0][rgbi] = a;
               }

               ac[0][3] = alpha;

            }
         }
      }

      /* Apply the spatial filter */
      y = de_offset;
      for (j = 0; j < image_height; j++) {
         x = de_offset;
         for (i = 0; i < image_width; i++) {
            int ii, jj,rgbi;
            void *p;
            unsigned short *p16;
            unsigned char *p8;
            t[0] = t[1] = t[2] = t[3] = 0.0;
            for (ii = 0; ii < filter_width; ii++) {
               for (jj = 0; jj < filter_width; jj++) {
                  double k = filter[ii + jj * filter_width];
                  abucket *ac = accumulate + x+ii + (y+jj)*fic.width;
                  

                  t[0] += k * ac[0][0];
                  t[1] += k * ac[0][1];
                  t[2] += k * ac[0][2];
                  t[3] += k * ac[0][3];


               }
            }

            p = out + nchan * bytes_per_channel * (i + j * out_width);
            p8 = (unsigned char *)p;
            p16 = (unsigned short *)p;
            
            /* The old way, spatial filter first and then clip after gamma */
            if (!eclip) {
            
               tmp=t[3]/PREFILTER_WHITE;
               
               if (t[3]<=0) {
                  alpha = 0.0;
                  ls = 0.0;
               } else { 
                  alpha = flam3_calc_alpha(tmp,g,linrange);
                  ls = vibrancy * 256.0 * alpha / tmp;
                  if (alpha<0.0) alpha = 0.0;
                  if (alpha>1.0) alpha = 1.0;
               }
              
               flam3_calc_newrgb(t, ls, highpow, newrgb);

               for (rgbi=0;rgbi<3;rgbi++) {
                  a = newrgb[rgbi];
                  a += (1.0-vibrancy) * 256.0 * pow( t[rgbi] / PREFILTER_WHITE, g);
                  if (nchan<=3 || transp==0)
                     a += ((1.0 - alpha) * background[rgbi]);
                  else {
                     if (alpha>0)
                        a /= alpha;
                     else
                        a = 0;
                  }

                  /* Clamp here to ensure proper filter functionality */
                  if (a>255) a = 255;
                  if (a<0) a = 0;
               
                  /* Replace values in accumulation buffer with these new ones */
                  t[rgbi] = a;
               }
               t[3] = alpha;
            }

            for (rgbi=0;rgbi<3;rgbi++) {

               a = t[rgbi];

               if (a > 255)
                  a = 255;
               if (a < 0)
                  a = 0;
               
               if (2==bytes_per_channel) {
                  a *= 256.0; /* Scales to [0-65280] */
                  p16[rgbi] = (unsigned short) a;
               } else {
                  p8[rgbi] = (unsigned char) a;
               }
            }


            if (t[3]>1)
               t[3]=1;
            if (t[3]<0)
               t[3]=0;

            /* alpha */
            if (nchan>3) {
               if (transp==1) {
                  if (2==bytes_per_channel)
                     p16[3] = (unsigned short) (t[3] * 65535);
                  else
                     p8[3] = (unsigned char) (t[3] * 255);
               } else {
                  if (2==bytes_per_channel)
                     p16[3] = 65535;
                  else
                     p8[3] = 255;
               }
            }

            x += oversample;
         }
         y += oversample;
      }
   }

 done:

   stats->badvals = fic.badvals;

   free(temporal_filter);
   free(temporal_deltas);
   free(batch_filter);
   free(filter);
   free(buckets);
//   free(accumulate);
//   free(points);
   /* We have to clear the cps in fth first */
   for (thi = 0; thi < spec->nthreads; thi++) {
      clear_cp(&(fth[thi].cp),0);
   }   
   free(fth);
   clear_cp(&cp,0);

   if (getenv("insert_palette")) {
     int ph = 100;
     if (ph >= image_height) ph = image_height;
     /* insert the palette into the image */
     for (j = 0; j < ph; j++) {
       for (i = 0; i < image_width; i++) {
	 unsigned char *p = out + nchan * (i + j * out_width);
	 p[0] = (unsigned char)dmap[i * 256 / image_width].color[0];
	 p[1] = (unsigned char)dmap[i * 256 / image_width].color[1];
	 p[2] = (unsigned char)dmap[i * 256 / image_width].color[2];
       }
     }
   }

   tend = time(NULL);
   stats->render_seconds = (int)(tend-tstart);

}
