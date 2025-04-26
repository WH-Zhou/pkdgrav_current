/*
 * colors.h -- DCR 94-09-29
 * ========================
 *
 */

#ifndef COLORS_H

/* Color type definition */

typedef unsigned char COLOR;

/* Other definitions */

#define NUM_COLORS        256 /* Use all available colors */
#define NUM_BRIGHT_COLORS 15  /* Excludes black */

#define BLACK   ((COLOR) 0) /* Black should ALWAYS be defined as 0 */
#define WHITE   ((COLOR) 1) /* White should ALWAYS be defined as 1 */
#define RED     ((COLOR) 2) /* Remaining "bright" colors */
#define GREEN   ((COLOR) 3)
#define BLUE    ((COLOR) 4)
#define YELLOW  ((COLOR) 5)
#define MAGENTA ((COLOR) 6)
#define CYAN    ((COLOR) 7)
#define GOLD	((COLOR) 8)
#define PINK    ((COLOR) 9)
#define ORANGE  ((COLOR) 10)
#define KHAKI	((COLOR) 11)
#define VIOLET	((COLOR) 12)
#define MAROON	((COLOR) 13)
#define AQUA	((COLOR) 14)
#define NAVY	((COLOR) 15)

/* Fill rest of colormap with gray scale */

#define FIRST_GRAY ((COLOR) 16)  /* Same as BLACK */
#define LAST_GRAY  ((COLOR) 255) /* Same as WHITE */
#define GRAY       ((COLOR) 136)

#define COLORS_H
#endif /* COLORS_H */

/* colors.h */
