#ifndef WRAPPER_HERMITE_POLYNOMIAL_H_INCLUDED
#define WRAPPER_HERMITE_POLYNOMIAL_H_INCLUDED

#define h_integral(a) R_DBL(FUNCNAME_HINTEGRAL(C_SUSHRT(a)))
#define h_polynomial_coefficients(a,b) R_UCHR(FUNCNAME_HPOLYNOMIALCOEFFICIENTS(C_DTPIT(a,b)))
#define h_polynomial_value(a,b) R_DBL(FUNCNAME_HPOLYNOMIALVALUE(C_DTIT(a,b)))
#define h_polynomial_zeros(a,b) R_UCHR(FUNCNAME_HPOLYNOMIALZEROS(C_DTPIT(a,b)))
#define h_quadrature_rule(a,b,c) R_UCHR(FUNCNAME_HQUADRATURERULE(C_DT2PIT(a,b,c)))
#define he_double_product_integral(a) R_DBL(FUNCNAME_HEDOUBLEPRODUCTINTEGRAL(a))
#define he_integral(a) R_DBL(FUNCNAME_HEINTEGRAL(C_SUSHRT(a)))
#define he_polynomial_coefficients(a,b) R_UCHR(FUNCNAME_HEPOLYNOMIALCOEFFICIENTS(C_DTPIT(a,b)))
#define he_polynomial_value(a,b) R_DBL(FUNCNAME_HEPOLYNOMIALVALUE(C_DTIT(a,b)))
#define he_polynomial_zeros(a,b) R_UCHR(FUNCNAME_HEPOLYNOMIALZEROS(C_DTPIT(a,b)))
#define he_quadrature_rule(a,b,c) R_UCHR(FUNCNAME_HEQUADRATURERULE(C_DT2PIT(a,b,c)))
#define he_triple_product_integral(a) R_DBL(FUNCNAME_HETRIPLEPRODUCTINTEGRAL(a))
#define hen_polynomial_value(a,b) R_DBL(FUNCNAME_HENPOLYNOMIALVALUE(C_DTIT(a,b)))
#define hf_function_value(a,b) R_DBL(FUNCNAME_HFFUNCTIONVALUE(C_DTIT(a,b)))
#define hf_quadrature_rule(a,b,c) R_UCHR(FUNCNAME_HFQUADRATURERULE(C_DT2PIT(a,b,c)))
#define hn_polynomial_value(a,b) R_DBL(FUNCNAME_HNPOLYNOMIALVALUE(C_DTIT(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_HPOLYNOMIALCOEFFICIENTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HPOLYNOMIALZEROS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HQUADRATURERULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEDOUBLEPRODUCTINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEPOLYNOMIALCOEFFICIENTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEPOLYNOMIALZEROS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEQUADRATURERULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HETRIPLEPRODUCTINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HFQUADRATURERULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HENPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HFFUNCTIONVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HNPOLYNOMIALVALUE(void *);

#endif
