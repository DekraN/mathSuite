#ifndef WRAPPER_VANDERMONDE_APPROX_1D_H_INCLUDED
#define WRAPPER_VANDERMONDE_APPROX_1D_H_INCLUDED

#define vandermonde_approx_1d_coef(a,b,c,d) FUNCNAME_VANDERMONDEAPPROX1DCOEF(C_2DT2PIT(a,b,c,d))
#define vandermonde_approx_1d_matrix(a,b,c) FUNCNAME_VANDERMONDEAPPROX1DMATRIX(C_2DTPIT(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_VANDERMONDEAPPROX1DCOEF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VANDERMONDEAPPROX1DMATRIX(void *);

#endif