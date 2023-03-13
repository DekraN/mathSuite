#ifndef WRAPPER_DISK_INTEGRALS_H_INCLUDED
#define WRAPPER_DISK_INTEGRALS_H_INCLUDED

#define disk01_area() R_DBL(FUNCNAME_DISK01AREA(NULL))
#define disk01_monomial_integral(a) R_DBL(FUNCNAME_DISK01MONOMIALINTEGRAL(C_SUSHRT(a)))
#define disk01_sample(a,b) FUNCNAME_DISK01SAMPLE(C_DTPI(a,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_DISK01MONOMIALINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DISK01SAMPLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DISK01AREA(void *);

#endif