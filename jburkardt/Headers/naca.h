#ifndef WRAPPER_NACA_H_INCLUDED
#define WRAPPER_NACA_H_INCLUDED

#define naca4_cambered(a,b,c,d,e,f,g,h,i,j) FUNCNAME_NACA4CAMBERED(C_4ITDT5PIT(a,b,c,d,e,f,g,h,i,j))
#define naca4_symmetric(a,b,c,d) FUNCNAME_NACA4SYMMETRIC(C_DT2ITPIT(c,a,b,d))

__MATHSUITE __JBURKARDT void * FUNCNAME_NACA4CAMBERED(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NACA4SYMMETRIC(void *);

#endif
