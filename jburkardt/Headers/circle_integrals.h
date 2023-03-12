#ifndef WRAPPER_CIRCLE_INTEGRALS_H_INCLUDED
#define WRAPPER_CIRCLE_INTEGRALS_H_INCLUDED

#define circle01_length() R_DBL(FUNCNAME_CIRCLE01LENGTH(NULL))
#define circle01_monomial_integral(a) R_DBL(FUNCNAME_CIRCLE01MONOMIALINTEGRAL(a))
#define circle01_sample(a,b,c) FUNCNAME_CIRCLE01SAMPLE(C_DTPITPI(a,c,b))
#define r8vec_uniform_01(a,b,c) FUNCNAME_R8VECUNIFORM01(C_DTPITPI(a,c,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLE01MONOMIALINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECUNIFORM01(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLE01SAMPLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLE01LENGTH(void *);

#endif
