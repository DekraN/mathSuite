#ifndef WRAPPER_POISSON_OPENMP_H_INCLUDED
#define WRAPPER_POISSON_OPENMP_H_INCLUDED

#define sweep(a,b,c,d,e,f,g,h,i) FUNCNAME_SWEEP(C_2DT2ITPPIT2DT2PPIT(a,b,c,d,e,f,g,h,i))
#define u_exact(a,b) R_DBL(FUNCNAME_UEXACT(C_PDBL2(a,b)))
#define uxxyy_exact(a,b) R_DBL(FUNCNAME_UXXYYEXACT(C_PDBL2(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_SWEEP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_UEXACT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_UXXYYEXACT(void *);

#endif
