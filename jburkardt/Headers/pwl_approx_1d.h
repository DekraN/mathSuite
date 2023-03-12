#ifndef WRAPPER_PWL_APPROX_1D_H_INCLUDED
#define WRAPPER_PWL_APPROX_1D_H_INCLUDED

#define pwl_approx_1d(a,b,c,d,e) FUNCNAME_PWLAPPROX1D(C_2DT3PIT(a,d,b,c,e))
#define pwl_approx_1d_matrix(a,b,c,d,e) FUNCNAME_PWLAPPROX1DMATRIX(C_2DT3PIT(a,d,b,c,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_PWLAPPROX1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PWLAPPROX1DMATRIX(void *);

#endif
