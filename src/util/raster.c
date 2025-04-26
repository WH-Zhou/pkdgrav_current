/*
 * raster.c -- DCR 94-08-30
 * ========================
 * Routines for handling raster data.
 *
 * Global functions: MakeColormap(), AllocImage(), DrawPoint(), DrawVector(),
 *   DumpRaster(), FreeImage().
 *
 */

#define USE_RLE /*DEBUG PBMplus doesn't like my RLE method; xv doesn't mind*/

/*** Include files ***/

#include <string.h>
#include <mem_util.h>
#include <rasterfile.h>
#include <scalar.h>
#include "raster.h"

/*** Local functions ***/

static void set_point(IMAGE_T *, int, int, COLOR);
static void clip(IMAGE_T *, int *, int *);
static void write_sun_long(unsigned long int, FILE *);
static void write_rle_data(IMAGE_T *, FILE *);

/*** End of preamble ***/

void MakeColormap(COLORMAP_T *colormap)
{
  /* Constructs colormap with simple colors plus gray scale */

  int i;
  COLOR red[NUM_COLORS], green[NUM_COLORS], blue[NUM_COLORS];
  double gray_scale;

  /* Set descriptive fields */

  colormap->type = RMT_EQUAL_RGB;
  colormap->size = NUM_COLORS;

  /* Allocate memory for colors */

  colormap->r = (COLOR *) Array(NUM_COLORS, sizeof(COLOR));
  colormap->g = (COLOR *) Array(NUM_COLORS, sizeof(COLOR));
  colormap->b = (COLOR *) Array(NUM_COLORS, sizeof(COLOR));

  /* Assign bright colours, including black */

  red[BLACK]   = 0;   green[BLACK]   = 0;   blue[BLACK]   = 0;
  red[WHITE]   = 255; green[WHITE]   = 255; blue[WHITE]   = 255;
  red[RED]     = 255; green[RED]     = 0;   blue[RED]     = 0;
  red[GREEN]   = 0;   green[GREEN]   = 255; blue[GREEN]   = 0;
  red[BLUE]    = 0;   green[BLUE]    = 0;   blue[BLUE]    = 255;
  red[YELLOW]  = 255; green[YELLOW]  = 255; blue[YELLOW]  = 0;
  red[MAGENTA] = 255; green[MAGENTA] = 0;   blue[MAGENTA] = 255;
  red[CYAN]    = 0;   green[CYAN]    = 255; blue[CYAN]    = 255;
  red[GOLD]    = 255; green[GOLD]    = 215; blue[GOLD]    = 0;
  red[PINK]    = 255; green[PINK]    = 192; blue[PINK]    = 203;
  red[ORANGE]  = 255; green[ORANGE]  = 165; blue[ORANGE]  = 0;
  red[KHAKI]   = 240; green[KHAKI]   = 230; blue[KHAKI]   = 140;
  red[VIOLET]  = 238; green[VIOLET]  = 130; blue[VIOLET]  = 238;
  red[MAROON]  = 176; green[MAROON]  = 48;  blue[MAROON]  = 96;
  red[AQUA]    = 127; green[AQUA]    = 255; blue[AQUA]    = 212;
  red[NAVY]    =   0; green[NAVY]    = 0;   blue[NAVY]    = 128;

  /* Assign gray scale */

  gray_scale = (double) LAST_GRAY / (LAST_GRAY - FIRST_GRAY);

  for (i = FIRST_GRAY; i <= LAST_GRAY; i++)
	  red[i] = green[i] = blue[i] = (COLOR) (gray_scale * (i - FIRST_GRAY));

  /* Copy data to colormap */

  for (i = 0; i < NUM_COLORS; i++) {
    colormap->r[i] = red[i];
    colormap->g[i] = green[i];
    colormap->b[i] = blue[i];
  }
}

void AllocImage(IMAGE_T **image, int w, int h)
{
  /* Allocates memory for "*image", width "w", height "h" */

  int i;

  *image = (IMAGE_T *) Elem(sizeof(IMAGE_T));

  (*image)->w = w;
  (*image)->h = h;

  /* Note order of rows, columns in following */

  (*image)->data = (COLOR **) Grid((size_t) h, (size_t) w, sizeof(COLOR));

  for (i = 0; i < h; i++)
	memset((*image)->data[i], 0, (size_t) w);
}

void DrawPoint(IMAGE_T *image, int x, int y, COLOR c)
{
  /* Draws a point (color "c") at position "x", "y" in "image" */

  if (!image || x < 0 || x >= image->w || y < 0 || y >= image->h)
    return;

  set_point(image, x, y, c);
}

void DrawVector(IMAGE_T *image, int x1, int y1, int x2, int y2, COLOR c)
{
  /*
   * Draws a line (color "c") between position "x1", "y1" and "x2", "y2"
   * in "image". Lines extending beyond image dimension are clipped.
   *
   */

  int x, y, dx, dy, x_min, y_min, x_max, y_max, sx;
  float m, b, ya0, yb0, ya, yb;

  dx = x2 - x1;
  dy = y2 - y1;

  x_min = MIN(x1, x2);
  y_min = MIN(y1, y2);
  x_max = MAX(x1, x2);
  y_max = MAX(y1, y2);

  if (x_max < 0 || x_min >= image->w || y_max < 0 || y_min >= image->h)
    return;

  clip(image, &x_min, &y_min);
  clip(image, &x_max, &y_max);

  if (dx == 0 && dy == 0)
    DrawPoint(image, x1, y1, c);
  else if (dx == 0)
    for (y = y_min; y <= y_max; y++)
      set_point(image, x1, y, c);
  else if (dy == 0)
    for (x = x_min; x <= x_max; x++)
      set_point(image, x, y1, c);
  else {
	m = (float) dy / (float) dx;
    b = (float) y1 - m * (float) x1;
    sx = SGN(dx);
    
    for (x = x_min; x <= x_max; x++) {
	  ya0 = m * (float) x + b;
      yb0 = ya0 + m * (float) sx;
      ya = MIN(ya0, yb0);
      yb = MAX(ya0, yb0);
      for (y = (int) ya; y <= (int) yb; y++)
		if (y >= y_min && y <= y_max)
		  set_point(image, x, y, c);
    }
  }
}

void DumpRaster(COLORMAP_T *colormap, IMAGE_T *image, const char *filename)
{
  /* Dumps "colormap" and "image" to rasterfile "filename" */

  FILE *fp;

  /* Open rasterfile */

  if ((fp = fopen(filename, "w")) == NULL) {
    perror(filename);
    return;
  }

  /* Dump header, ensuring integers are output in correct (Sun) byte order */

  write_sun_long(RAS_MAGIC, fp);
  write_sun_long((size_t) image->w, fp);
  write_sun_long((size_t) image->h, fp);
  write_sun_long(8, fp);
  write_sun_long((size_t) ((image->w + (image->w % 2)) * image->h), fp);
#ifdef USE_RLE
  write_sun_long(RT_BYTE_ENCODED, fp);
#else
  write_sun_long(RT_STANDARD, fp);
#endif
  write_sun_long((size_t) colormap->type, fp);
  write_sun_long(colormap->size * 3, fp);

  /* Dump colormap */

  (void) fwrite(colormap->r, sizeof(COLOR), colormap->size, fp);
  (void) fwrite(colormap->g, sizeof(COLOR), colormap->size, fp);
  (void) fwrite(colormap->b, sizeof(COLOR), colormap->size, fp);

  /* Dump raster in run-length-encoded format */

#ifdef USE_RLE
  write_rle_data(image, fp);
#else
  {
    int i;

    for (i = 0; i < image->h; i++)
      (void) fwrite(image->data[i], sizeof(COLOR), image->w, fp);
  }
#endif

  /* Close rasterfile */

  if (fclose(fp) == EOF)
    perror(filename);
}

void FreeImage(IMAGE_T *image)
{
  /* Deallocates memory associated with "image" */

  FreeGrid((size_t) image->h, (void **) image->data);
  FreeElem(image);
}

static void set_point(IMAGE_T *image, int x, int y, COLOR c)
{
  /* Assigns value "c" to data element at "x", "y" in "image" */

  image->data[y][x] = c; /* Note order of rows and columns */
}

static void clip(IMAGE_T *image, int *x, int *y)
{
  /* Ensures "*x", "*y" lie in "image" */

  *x = MAX(*x, 0);
  *x = MIN(*x, image->w - 1);

  *y = MAX(*y, 0);
  *y = MIN(*y, image->h - 1);
}

static void write_sun_long(unsigned long int l, FILE *fp)
{
  /* Writes a Sun-ordered unsigned long "l" to "fp". Taken from xv 3.01 */

  char c;

  c = (char) ((l >> 24) & 0xff);
  (void) putc(c, fp);
  c = (char) ((l >> 16) & 0xff);
  (void) putc(c, fp);
  c = (char) ((l >> 8) & 0xff);
  (void) putc(c, fp);
  c = (char) (l & 0xff);
  (void) putc(c, fp);
}

#define RLE_ESC 128 /* Escape character for RLE encoding */

/*DEBUG I/O error checks needed*/
static void write_rle_data(IMAGE_T *image, FILE *fp)
{
  /* Writes "image" data to "fp" in run-length-encoded (RLE) format */

  int i, j, odd;
  unsigned char c, cc=0, nc;

  odd = image->w % 2;

  i = j = 0;

  c = image->data[0][0];

  while (1) {
    nc = 0;
    while (1) {
      j++;
      if (j == image->w) {
	i++;
	j = 0;
	if (odd)
	  (void) putc(0, fp); /* Image lines rounded to multiple of 16 bits */
      }
      if (i == image->h)
	break;
      cc = image->data[i][j];
      if (cc != c || nc == 255)
	break;
      ++nc;
    }
    if (nc == 0) {
      if (c == RLE_ESC) {
	(void) putc(RLE_ESC, fp);
	(void) putc(0, fp);
      }
      else
	(void) putc(c, fp);
    }
    else if (nc == 1 && c != RLE_ESC) {
      (void) putc(c, fp);
      (void) putc(c, fp);
    }
    else {
      (void) putc(RLE_ESC, fp);
      (void) putc(nc, fp);
      (void) putc(c, fp);
    }
    if (i == image->h)
      break;
    c = cc;
  }
}

#undef RLE_ESC

/* raster.c */
