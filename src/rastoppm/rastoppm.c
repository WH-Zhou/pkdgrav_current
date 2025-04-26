/* rastoppm.c -- converts Sun rasterfiles into ppm format */

/* DCR 98-01-27 */

#include <stdio.h>
#include <stdlib.h>
#include "rasterfile.h"

#ifndef _SYS_BSD_TYPES_H
#ifndef _SYS_TYPES_H
typedef unsigned char u_char;
typedef unsigned long u_long;
#endif
#endif

#define RLE_ESC 128 /* escape character for RLE encoding */

#define PUT(c)\
{\
	i += (num_colors == 0 ? 8 : 1);\
	if (!odd_byte || i % (w + 1)) {\
		if (num_colors == 0)\
			bits_to_bytes(c,fpo);\
		else {\
			(void) fputc(r[c],fpo);\
			(void) fputc(g[c],fpo);\
			(void) fputc(b[c],fpo);\
		}\
	}\
}

/*** END OF PREAMBLE ***/

static void
bits_to_bytes(int c,FILE *fp)
{
	int i,m;

	for (i=0,m=0x80;i<8;i++,m>>=1) {
		if ((c & m) == m) { /* why does on mean off?... */
			(void) fputc(0,fp);
			(void) fputc(0,fp);
			(void) fputc(0,fp);
			}
		else {
			(void) fputc(255,fp);
			(void) fputc(255,fp);
			(void) fputc(255,fp);
			}
		}
	}

static int
read_sun_long(FILE *fp,int *l)
{
  /*
   * Taken by permission from read_sun_long() in xvsunras.c from the
   * xv 3.10a distribution, Copyright (C) 1994 by John Bradley.
   *
   */

  int c0,c1,c2,c3;

  c0 = fgetc(fp);
  c1 = fgetc(fp);
  c2 = fgetc(fp);
  c3 = fgetc(fp);

  *l = (int) ((((u_long) c0 & 0xff) << 24) |
       (((u_long) c1 & 0xff) << 16) |
       (((u_long) c2 & 0xff) <<  8) |
	   (((u_long) c3 & 0xff)));

  if (ferror(fp) || feof(fp))
    return EOF;

  return 0;
}

static int
convert(FILE *fpi,FILE *fpo)
{
	struct rasterfile header;
	u_char *r=NULL,*g=NULL,*b=NULL;
	int i,c,w,h,num_colors,odd_byte,ras_length;

	/* Read header, taking byte ordering into account */

	(void) read_sun_long(fpi,&header.ras_magic);
	(void) read_sun_long(fpi,&header.ras_width);
	(void) read_sun_long(fpi,&header.ras_height);
	(void) read_sun_long(fpi,&header.ras_depth);
	(void) read_sun_long(fpi,&header.ras_length);
	(void) read_sun_long(fpi,&header.ras_type);
	(void) read_sun_long(fpi,&header.ras_maptype);
	if (read_sun_long(fpi,&header.ras_maplength) != 0) {
		(void) fprintf(stderr,"Error reading ras header.\n");
		return 1;
		}
	if (header.ras_magic != RAS_MAGIC) {
		(void) fprintf(stderr,"Input not rasterfile.\n");
		return 1;
		}
	if (header.ras_depth != 1 && header.ras_depth != 8) {
		(void) fprintf(stderr,"Rasterfile must have 1- or 8-bit depth.\n");
		return 1;
		}
	switch (header.ras_type) {
	case RT_OLD:
	case RT_STANDARD:
	case RT_BYTE_ENCODED:
		break;
	default:
		(void) fprintf(stderr,"Unsupported ras encoding.\n");
		return 1;
		}
	num_colors = header.ras_maplength/3;
	switch (header.ras_maptype) {
	case RMT_EQUAL_RGB:
		if (num_colors <= 0) {
			(void) fprintf(stderr,"No colors found!\n");
			return 1;
			}
		break;
	case RMT_NONE:
		if (num_colors != 0) {
			(void) fprintf(stderr,"Expected zero-length colormap\n");
			return 1;
			}
		break;
	default:
		(void) fprintf(stderr,"Unsupported colormap.\n");
		return 1;
		}

	if (num_colors > 0) {
		r = (u_char *) malloc((size_t) num_colors);
		g = (u_char *) malloc((size_t) num_colors);
		b = (u_char *) malloc((size_t) num_colors);

		(void) fread((void *) r,1,(size_t) num_colors,fpi);
		(void) fread((void *) g,1,(size_t) num_colors,fpi);
		(void) fread((void *) b,1,(size_t) num_colors,fpi);
		}

	w = header.ras_width;
	h = header.ras_height;

	odd_byte = w % 2; /* ras image lines rounded to 16 bits */
	ras_length = (w + odd_byte)*h;
	if (header.ras_type == RT_STANDARD && header.ras_length != ras_length) {
		(void) fprintf(stderr,"Inconsistent ras data length.");
		return 1;
		}
	if (ras_length <= 0) {
		(void) fprintf(stderr,"No ras data found!");
		return 1;
		}

	(void) fprintf(fpo,"P6\n%d %d\n255\n",w,h); /* ppm header */

	if (header.ras_type == RT_BYTE_ENCODED) { /* RLE raster */
		int nc;
		i = 0;
		while (i < ras_length) {
			if ((c = fgetc(fpi)) < 0) {
				(void) fprintf(stderr,"Corrupt rasterfile.\n");
				return 1;
				}
			if (c == RLE_ESC) {
				if ((nc = fgetc(fpi)) < 0) {
					(void) fprintf(stderr,"Corrupt rasterfile.\n");
					return 1;
					}
				if (nc == 0) PUT(RLE_ESC) /* note: no semi-colon */
				else if (nc == 1) {
					if ((c = fgetc(fpi)) != RLE_ESC) {
						(void) fprintf(stderr,"Corrupt RLE ras file.\n");
						return 1;
						}
					PUT(RLE_ESC)
					PUT(RLE_ESC)
					}
				else {
					c = fgetc(fpi);
					if (c < 0 || (num_colors > 0 && c >= num_colors)) {
						(void) fprintf(stderr,"Corrupt rasterfile.\n");
						return 1;
						}
					for (;nc>=0;--nc) PUT(c)
					}
				}
			else {
				if (num_colors > 0 && c >= num_colors) {
					(void) fprintf(stderr,"Corrupt rasterfile.\n");
					return 1;
					}
				PUT(c)
				}
			}
		}
	else {/* standard raster */
		int j;
		for (i=0;i<h;++i) {
			for (j=0;j<w;++j) {
				c = fgetc(fpi);
				if (c < 0 || (num_colors > 0 && c >= num_colors)) {
					(void) fprintf(stderr,"Corrupt rasterfile.\n");
					return 1;
					}
				if (num_colors == 0)
					bits_to_bytes(c,fpo);
				else {
					(void) fputc(r[c],fpo);
					(void) fputc(g[c],fpo);
					(void) fputc(b[c],fpo);
					}
				}
			if (odd_byte)
				(void) fgetc(fpi);
			}
		}

	if (num_colors > 0) {
		free((void *) b);
		free((void *) g);
		free((void *) r);
		}

	return 0;
	}

int
main(int argc,char *argv[])
{
	FILE *fpi=stdin,*fpo=stdout;

	if (argc > 3) {
		(void) fprintf(stderr,"Usage: %s [ ras_file [ ppm_file ] ]\n",argv[0]);
		return 1;
		}

	if (argc > 1 && (fpi = fopen(argv[1],"r")) == NULL) {
		(void) fprintf(stderr,"Unable to open \"%s\" for reading\n",argv[1]);
		return 1;
		}

	if (argc == 3 && (fpo = fopen(argv[2],"w")) == NULL) {
		(void) fprintf(stderr,"Unable to open \"%s\" for writing\n",argv[2]);
		(void) fclose(fpi);
		return 1;
		}

	return convert(fpi,fpo);
	}

/* rastoppm.c */
