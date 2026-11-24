#ifndef WRAPPER_LATIN_COVER_H_INCLUDED
#define WRAPPER_LATIN_COVER_H_INCLUDED

#define latin_cover(a,b) FUNCNAME_LATINCOVER(C_DTPI(a,b))
#define latin_cover_2d(a,b,c) FUNCNAME_LATINCOVER2D(C_DT2PI(a,b,c))
#define latin_cover_3d(a,b,c,d) FUNCNAME_LATINCOVER3D(C_DT3PI(a,b,c,d))

__MATHSUITE __JBURKARDT void * FUNCNAME_LATINCOVER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LATINCOVER2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LATINCOVER3D(void *);

#endif
