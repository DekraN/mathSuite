#ifndef WRAPPER_POLYGON_TRIANGULATE_H_INCLUDED
#define WRAPPER_POLYGON_TRIANGULATE_H_INCLUDED

#define between(a,b,c,d,e,f) R_UCHR(FUNCNAME_BETWEEN(C_PDBL6(a,b,c,d,e,f)))
#define collinear(a,b,c,d,e,f) R_UCHR(FUNCNAME_COLLINEAR(C_PDBL6(a,b,c,d,e,f)))
#define diagonal(a,b,c,d,e,f,g) R_UCHR(FUNCNAME_DIAGONAL(C_3DT2PI2PIT(a,b,c,d,e,f,g)))
#define diagonalie(a,b,c,d,e,f) R_UCHR(FUNCNAME_DIAGONALIE(C_DTPIT2DTPIPIT(a,e,b,c,d,f)))
#define in_cone(a,b,c,d,e,f,g) R_UCHR(FUNCNAME_INCONE(C_3DT2PI2PIT(a,b,c,d,e,f,g)))
#define intersect(a,b,c,d,e,f,g,h) R_UCHR(FUNCNAME_INTERSECT(C_PDBL8(a,b,c,d,e,f,g,h)))
#define intersect_prop(a,b,c,d,e,f,g,h) R_UCHR(FUNCNAME_INTERSECTPROP(C_PDBL8(a,b,c,d,e,f,g,h)))
#define triangle_area_ptriang(a,b,c,d,e,f) R_DBL(FUNCNAME_TRIANGLEAREAPTRIANG(C_PDBL6(a,b,c,d,e,f)))


__MATHSUITE __JBURKARDT void * FUNCNAME_DIAGONAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIAGONALIE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_INCONE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETWEEN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_COLLINEAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEAREAPTRIANG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_INTERSECT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_INTERSECTPROP(void *);

#endif
