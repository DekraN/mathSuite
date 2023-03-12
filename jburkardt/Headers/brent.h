#ifndef WRAPPER_BRENT_H_INCLUDED
#define WRAPPER_BRENT_H_INCLUDED

#define glomin(a,b,c,d,e,f,g,h,i) R_DBL(FUNCNAME_GLOMIN(C_7ITFITPIT(a,b,c,d,e,f,g,h,i)))
#define local_min(a,b,c,d,e,f) R_DBL(FUNCNAME_LOCALMIN(C_4PITFITPIT(a,b,c,d,e,f)))
#define r8_epsilon() R_DBL(FUNCNAME_R8EPSILON(NULL))
#define zero(a,b,c,d,e) R_DBL(FUNCNAME_ZERO(C_4PITFIT(a,b,c,d,e)))

__MATHSUITE __JBURKARDT void * FUNCNAME_GLOMIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOCALMIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ZERO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8EPSILON(void *);

#endif
