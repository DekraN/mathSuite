#ifndef WRAPPER_LINE_INTEGRALS_H_INCLUDED
#define WRAPPER_LINE_INTEGRALS_H_INCLUDED

#define line01_length() R_DBL(FUNCNAME_LINE01LENGTH(NULL))
#define line01_monomial_integral(a) R_DBL(FUNCNAME_LINE01MONOMIALINTEGRAL(C_SUSHRT(a)))
#define line01_sample(a,b) FUNCNAME_LINE01SAMPLE(C_DTPI(a,b))
#define monomial_value_1d(a,b,c) FUNCNAME_MONOMIALVALUE1D(C_2DTPIT(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_LINE01SAMPLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MONOMIALVALUE1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINE01LENGTH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINE01MONOMIALINTEGRAL(void *);

#endif