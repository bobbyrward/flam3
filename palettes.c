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

#include "private.h"
#include "palettes.h"

lib_palette *the_palettes = NULL;
int npalettes;

static void parse_palettes(xmlNode *node) {
   xmlAttrPtr attr;
   char *val;
   lib_palette *pal;
    
   while (node) {
      if (node->type == XML_ELEMENT_NODE && !xmlStrcmp(node->name, (const xmlChar *)"palette")) {
         attr = node->properties;
	      pal = &the_palettes[npalettes];
	      memset(pal, 0, sizeof(lib_palette));

         while (attr) {
            val = (char *) xmlGetProp(node, attr->name);
            if (!xmlStrcmp(attr->name, (const xmlChar *)"data")) {
               int count = 256;
               int c_idx = 0;
               int r,g,b;
               int col_count = 0;
               int sscanf_ret;
               
               c_idx=0;
               col_count = 0;
               
               do {
                  sscanf_ret = sscanf((char *)&(val[c_idx]),"00%2x%2x%2x",&r,&g,&b);
                  if (sscanf_ret != 3) {
                     printf("Error:  Problem reading hexadecimal color data.\n");
                     exit(1);
                  }
			         c_idx += 8;
                  while (isspace( (int)val[c_idx]))
                     c_idx++;

                  pal->colors[col_count][0] = r;
                  pal->colors[col_count][1] = g;
                  pal->colors[col_count][2] = b;
                  
                  col_count++;
               } while (col_count<count);
            } else if (!xmlStrcmp(attr->name, (const xmlChar *)"number")) {
               pal->number = atoi(val);
            } else if (!xmlStrcmp(attr->name, (const xmlChar *)"name")) {
               strncpy(pal->name, val, flam3_name_len);
               pal->name[flam3_name_len-1] = 0;
            }
            
            xmlFree(val);
            attr = attr->next;
         }
         
         npalettes++;
         the_palettes = realloc(the_palettes, (1 + npalettes) * sizeof(lib_palette));
      } else
         parse_palettes(node->children);
         
	   node = node->next;
   }
}

static int init_palettes(char *filename) {
    FILE *fp;
    xmlDocPtr doc;
    xmlNode *rootnode;
   int i, c, slen = 5000;
   char *s;

    fp = fopen(filename, "rb");
    if (NULL == fp) {
	fprintf(stderr, "flam3: could not open palette file ");
	perror(filename);
	exit(1);
    }
   
   /* Incrementally read XML file into a string */
   s = malloc(slen);
   i = 0;
   do {
      c = getc(fp);
      if (EOF == c) {
	  if (ferror(fp)) {
	      perror(filename);
	      exit(1);
	  }
         break;
      }
      s[i++] = c;
      if (i == slen-1) {
         slen *= 2;
         s = realloc(s, slen);
      }
   } while (1);

   fclose(fp);
   s[i] = 0;
   
   doc = xmlReadMemory(s, (int)strlen(s), filename, NULL, XML_PARSE_NONET);
   if (NULL == doc) {
       fprintf(stderr, "error parsing %s (%s).\n", filename, s);
       exit(1);
   }
   rootnode = xmlDocGetRootElement(doc);
   the_palettes = malloc(sizeof(lib_palette));
   npalettes = 0;
   parse_palettes(rootnode);
   xmlFreeDoc(doc);
   
   free(s);
   xmlCleanupParser();
   return(1);
}

int flam3_get_palette(int n, flam3_palette c, double hue_rotation) {
   int cmap_len = 256;
   int idx, i, j;

   if (NULL == the_palettes) {
   
      char *d = getenv("flam3_palettes");
      init_palettes(d ? d : (PACKAGE_DATA_DIR "/flam3-palettes.xml"));
   }
   
   if (flam3_palette_random == n)
      n = the_palettes[random()%npalettes].number;

   for (idx = 0; idx < npalettes; idx++) {
      
      if (n == the_palettes[idx].number) {
      	/* Loop over elements of cmap */
	      for (i = 0; i < cmap_len; i++) {
            int ii = (i * 256) / cmap_len;
            double rgb[3], hsv[3];
            
            /* Colors are in 0-1 space */
            for (j = 0; j < 3; j++)
               rgb[j] = the_palettes[idx].colors[ii][j] / 255.0;
	       
	         rgb2hsv(rgb, hsv);
	         hsv[0] += hue_rotation * 6.0;
	         hsv2rgb(hsv, rgb);
	         
	         c[i].index = i;
	       
	         for (j = 0; j < 3; j++)
		         c[i].color[j] = rgb[j];
		         
		      c[i].color[3] = 1.0;
         }
	   
	      return n;
      }
   }
   
   fprintf(stderr, "warning: palette number %d not found, using white.\n", n);

   for (i = 0; i < cmap_len; i++) {
   	c[i].index = i;
      for (j = 0; j < 4; j++)
         c[i].color[j] = 1.0;
   }

   return flam3_palette_random;
}

/* rgb 0 - 1,
   h 0 - 6, s 0 - 1, v 0 - 1 */
void rgb2hsv(rgb, hsv)
   double *rgb; double *hsv;
 {
  double rd, gd, bd, h, s, v, max, min, del, rc, gc, bc;

  rd = rgb[0];
  gd = rgb[1];
  bd = rgb[2];

  /* compute maximum of rd,gd,bd */
  if (rd>=gd) { if (rd>=bd) max = rd;  else max = bd; }
         else { if (gd>=bd) max = gd;  else max = bd; }

  /* compute minimum of rd,gd,bd */
  if (rd<=gd) { if (rd<=bd) min = rd;  else min = bd; }
         else { if (gd<=bd) min = gd;  else min = bd; }

  del = max - min;
  v = max;
  if (max != 0.0) s = (del) / max;
             else s = 0.0;

  h = 0;
  if (s != 0.0) {
    rc = (max - rd) / del;
    gc = (max - gd) / del;
    bc = (max - bd) / del;

    if      (rd==max) h = bc - gc;
    else if (gd==max) h = 2 + rc - bc;
    else if (bd==max) h = 4 + gc - rc;

    if (h<0) h += 6;
  }

  hsv[0] = h;
  hsv[1] = s;
  hsv[2] = v;
}


/* h 0 - 6, s 0 - 1, v 0 - 1
   rgb 0 - 1 */
void hsv2rgb(hsv, rgb)
   double *hsv;
   double *rgb;
{
   double h = hsv[0], s = hsv[1], v = hsv[2];
  int    j;
  double rd, gd, bd;
  double f, p, q, t;

   while (h >= 6.0) h = h - 6.0;
   while (h <  0.0) h = h + 6.0;
   j = (int) floor(h);
   f = h - j;
   p = v * (1-s);
   q = v * (1 - (s*f));
   t = v * (1 - (s*(1 - f)));
   
   switch (j) {
    case 0:  rd = v;  gd = t;  bd = p;  break;
    case 1:  rd = q;  gd = v;  bd = p;  break;
    case 2:  rd = p;  gd = v;  bd = t;  break;
    case 3:  rd = p;  gd = q;  bd = v;  break;
    case 4:  rd = t;  gd = p;  bd = v;  break;
    case 5:  rd = v;  gd = p;  bd = q;  break;
    default: rd = v;  gd = t;  bd = p;  break;
   }

   rgb[0] = rd;
   rgb[1] = gd;
   rgb[2] = bd;
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
   double a, maxa=-1.0, maxc=0;
   
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


#if 0
#if 0
void
cmap2image(n, f)
   int n;
   FILE *f;
{
   int i;

   fprintf(f, "P3\n16 16\n# cmap %d\n255\n", n);
   for (i = 0; i < 256; i++)
      fprintf(f, "%d %d %d ",
	      the_cmaps[n][i][0],
	      the_cmaps[n][i][1],
	      the_cmaps[n][i][2]);
   fprintf(f, "\n");
}
   
void
cmap2image2(n, f)
   int n;
   FILE *f;
{
   int i, j;

   fprintf(f, "P3\n50 256\n# cmap %d\n255\n", n);
   for (i = 0; i < 256; i++)
       for (j = 0; j < 50; j++)
	   fprintf(f, "%d %d %d ",
		   the_cmaps[n][i][0],
		   the_cmaps[n][i][1],
		   the_cmaps[n][i][2]);
   fprintf(f, "\n");
}
   
#endif
main(argc, argv)
   int argc;
   char **argv;
{
#if 0
  int i;
  FILE *f;

  /*  printf("%d\n", vlen(the_cmaps));
      exit(0); */
  for (i = 0 ; i < vlen(the_cmaps); i++) {
    char fname[100];
    sprintf(fname, "cmap%04d.ppm", i);
    f = fopen(fname, "w");
    if (NULL == f) {
      perror(fname);
      exit(1);
    }
    cmap2image2(i, f);
    fclose(f);
  }
#elsif 0
   int n = atoi(argv[1]);
   if (n >= 0 && n < vlen(the_cmaps))
       cmap2image2(n, stdout);
   else {
       fprintf(stderr, "bad cmap index: %d.\n", n);
       exit(1);
   }
#else
   int i, j;
   init_palettes("pal.xml");

   printf("<palettes>\n");
   for (i = 0; i < npalettes; i++) {
       printf("<palette number=\"%d\" data=\"", the_palettes[i].number);
       for (j = 0; j < 256; j++) {
	   printf("00%02x%02x%02x",
		  the_palettes[i].colors[j][0],
		  the_palettes[i].colors[j][1],
		  the_palettes[i].colors[j][2]);
	   if ((j&7)==7) printf("\n");
       }
       printf("\"/>\n");
   }
   printf("</palettes>\n");

#endif
}
#endif
