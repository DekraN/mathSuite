#ifndef WRAPPER_MULTIGRID_POISSON_1D_H_INCLUDED
#define WRAPPER_MULTIGRID_POISSON_1D_H_INCLUDED

#define ctof(a,b,c,d) FUNCNAME_CTOF(C_2DT2PIT(a,c,b,d))
#define ftoc(a,b,c,d,e,f) FUNCNAME_FTOC(C_2DT4PIT(a,d,b,c,e,f))
#define gauss_seidel(a,b,c,d) FUNCNAME_GAUSSSEIDEL(C_DT3PIT(a,b,c,d))
#define monogrid_poisson_1d(a,b,c,d,e,f,g,h,i) FUNCNAME_MONOGRIDPOISSON1D(C_DT4IT2FITPDTPIT(a,b,c,d,e,f,g,h,i))
#define multigrid_poisson_1d(a,b,c,d,e,f,g,h,i) FUNCNAME_MULTIGRIDPOISSON1D(C_DT4IT2FITPDTPIT(a,b,c,d,e,f,g,h,i))

__MATHSUITE __JBURKARDT void * FUNCNAME_CTOF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FTOC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAUSSSEIDEL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MONOGRIDPOISSON1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MULTIGRIDPOISSON1D(void *);


#endif
