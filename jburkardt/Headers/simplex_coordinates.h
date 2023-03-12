#ifndef WRAPPER_SIMPLEX_COORDINATES_H_INCLUDED
#define WRAPPER_SIMPLEX_COORDINATES_H_INCLUDED

#define simplex_coordinates1(a) FUNCNAME_SIMPLEXCOORDINATES1(C_SUSHRT(a))
#define simplex_coordinates2(a) FUNCNAME_SIMPLEXCOORDINATES2(C_SUSHRT(a))
#define simplex_volume(a,b) R_DBL(FUNCNAME_SIMPLEXVOLUME(C_DTPIT(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEXVOLUME(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEXCOORDINATES1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEXCOORDINATES2(void *);

#endif
