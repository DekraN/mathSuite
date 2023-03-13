#ifndef WRAPPER_QUADMOM_H_INCLUDED
#define WRAPPER_QUADMOM_H_INCLUDED

#define moment_method(a,b,c) FUNCNAME_MOMENTMETHOD(C_DT3PIT(a,b,c))
#define moments_laguerre(a) FUNCNAME_MOMENTSLAGUERRE(C_SUSHRT(a))
#define moments_legendre(a,b,c) FUNCNAME_MOMENTSLEGENDRE(C_2ITDT(b,c,a))
#define moments_normal(a,b,c) FUNCNAME_MOMENTSNORMAL(C_2ITDT(b,c,a))
#define moments_truncated_normal_ab(a,b,c,d,e) FUNCNAME_MOMENTSTRUNCATEDNORMALAB(C_DT4IT(a,b,c,d,e))
#define moments_truncated_normal_a(a,b,c,d) FUNCNAME_MOMENTSTRUNCATEDNORMALA(C_DT3IT(a,b,c,d))
#define moments_truncated_normal_b(a,b,c,d) FUNCNAME_MOMENTSTRUNCATEDNORMALB(C_DT3IT(a,b,c,d))
#define normal_01_cdf(a) R_DBL(FUNCNAME_NORMAL01CDF(C_SDBL(a)))
#define normal_01_pdf(a) R_DBL(FUNCNAME_NORMAL01PDF(C_SDBL(a)))
#define truncated_normal_ab_moment(a,b,c,d,e) R_DBL(FUNCNAME_TRUNCATEDNORMALABMOMENT(C_DT4IT(a,b,c,d,e)))
#define truncated_normal_a_moment(a,b,c,d) R_DBL(FUNCNAME_TRUNCATEDNORMALAMOMENT(C_DT3IT(a,b,c,d)))
#define truncated_normal_b_moment(a,b,c,d) R_DBL(FUNCNAME_TRUNCATEDNORMALBMOMENT(C_DT3IT(a,b,c,d)))
#define r8mat_cholesky_factor_upper(a,b,c) FUNCNAME_R8MATCHOLESKYFATTORUPPER(C_DTPITPB(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTMETHOD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTSLAGUERRE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTSLEGENDRE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTSNORMAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTSTRUNCATEDNORMALAB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTSTRUNCATEDNORMALA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTSTRUNCATEDNORMALB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATCHOLESKYFATTORUPPER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NORMAL01CDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NORMAL01PDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRUNCATEDNORMALAMOMENT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRUNCATEDNORMALABMOMENT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRUNCATEDNORMALBMOMENT(void *);

#endif