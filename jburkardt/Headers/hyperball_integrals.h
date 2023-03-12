#ifndef WRAPPER_HYPERBALL_INTEGRALS_H_INCLUDED
#define WRAPPER_HYPERBALL_INTEGRALS_H_INCLUDED

#define r8mat_normal_01_new(a,b,c) FUNCNAME_R8MATNORMAL01NEW(C_2DTPI(a,b,c))
#define hyperball01_monomial_integral(a,b) R_DBL(FUNCNAME_HYPERBALL01MONOMIALINTEGRAL(C_DTPI(a,b)))
#define hyperball01_sample(a,b,c) FUNCNAME_HYPERBALL01SAMPLE(C_2DTPI(a,b,c))
#define hyperball01_volume(a) R_DBL(FUNCNAME_HYPERBALL01VOLUME(C_SUSHRT(a)))

__MATHSUITE __JBURKARDT void * FUNCNAME_HYPERBALL01MONOMIALINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATNORMAL01NEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HYPERBALL01SAMPLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HYPERBALL01VOLUME(void *);

#endif
