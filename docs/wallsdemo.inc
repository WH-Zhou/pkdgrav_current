/*
** povray.inc -- DCR 94-10-19
** ==========================
** Definitions for use with POV-Ray scene files.
*/

#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"
#include "stones.inc"

/* definitions */

#declare Bright = texture {pigment {color White} finish {ambient 1.0}}

#declare ParticleSolid = texture {pigment {color White} finish {phong 0.2 diffuse 1.0 ambient 0.02}}
#declare ParticleFancy = texture {T_Stone5 pigment {color LightSteelBlue} finish {ambient 0.5}}
#declare ParticleBright = texture {Bright}

#declare StyleParticle = texture {ParticleSolid}

#declare WallPlain = texture {pigment {transmit 1.0} finish {ambient 0.02}}
#declare WallAgate = texture {pigment {agate}}
#declare WallBrick = texture {pigment {brick brick_size 0.03 scale <1.5,1,1> mortar 0.003}}

/* defaults */

#declare StylePlane		= texture {WallBrick}
#declare StyleTriangle		= texture {WallPlain}
#declare StyleRectangle		= texture {WallAgate}
#declare StyleDisk		= texture {WallAgate}
#declare StyleCylinderInfinite	= texture {WallPlain}
#declare StyleCylinderFinite	= texture {WallAgate}
#declare StyleShell		= texture {WallAgate}

/* default background color */

background {color Black}

/* povray.inc */
