#ifndef WRAPPER_LAGUERRE_POLYNOMIAL_H_INCLUDED
#define WRAPPER_LAGUERRE_POLYNOMIAL_H_INCLUDED

#define l_polynomial(a,b) R_DBL(FUNCNAME_LPOLYNOMIAL(C_DTIT(a,b)))
#define l_polynomial_coefficients(a,b) R_UCHR(FUNCNAME_LPOLYNOMIALCOEFFICIENTS(C_DTPIT(a,b)))
#define lm_polynomial(a,b,c) R_DBL(FUNCNAME_LMPOLYNOMIAL(C_2DTIT(a,b,c)))
#define lm_polynomial_coefficients(a,b,c) R_UCHR(FUNCNAME_LMPOLYNOMIALCOEFFICIENTS(C_2DTPIT(a,b,c)))
#define lf_function(a,b,c) R_DBL(FUNCNAME_LFFUNCTION(C_2ITDT(b,c,a)))
#define l_quadrature_rule(a,b,c) R_UCHR(FUNCNAME_LQUADRATURERULE(C_DT2PIT(a,b,c)))
#define lf_quadrature_rule(a,b,c,d) R_UCHR(FUNCNAME_LFQUADRATURERULE(C_DTIT2PIT(a,b,c,d)))
#define lm_quadrature_rule(a,b,c,d) R_UCHR(FUNCNAME_LMQUADRATURERULE(C_2DT2PIT(a,b,c,d)))
#define l_integral(a) R_INT(FUNCNAME_LINTEGRAL(C_SUSHRT(a)))
#define lf_integral(a,b) R_DBL(FUNCNAME_LFINTEGRAL(C_DTIT(a,b)))
#define lm_integral(a,b) R_DBL(FUNCNAME_LMINTEGRAL(C_PUSHRT2(a,b)))
#define l_polynomial_zeros(a,b) R_UCHR(FUNCNAME_LPOLYNOMIALZEROS(C_DTPIT(a,b)))
#define lf_function_zeros(a,b,c) R_UCHR(FUNCNAME_LFFUNCTIONZEROS(C_DTPITIT(a,c,b)))
#define lm_polynomial_zeros(a,b,c) R_UCHR(FUNCNAME_LMPOLYNOMIALZEROS(C_2DTPIT(a,b,c)))
#define lm_polynomial_values(a,b,c,d,e) FUNCNAME_LMPOLYNOMIALVALUES(C_3PDT2PIT(a,b,c,d,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_LPOLYNOMIALCOEFFICIENTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LMPOLYNOMIALCOEFFICIENTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LQUADRATURERULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LFQUADRATURERULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LMQUADRATURERULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LPOLYNOMIALZEROS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LMPOLYNOMIALZEROS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LFFUNCTIONZEROS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LPOLYNOMIAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LFINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LMINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LMPOLYNOMIAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LFFUNCTION(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LMPOLYNOMIALVALUES(void *);

#endif
