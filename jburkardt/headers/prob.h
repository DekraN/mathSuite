#ifndef WRAPPER_PROB_H_INCLUDED
#define WRAPPER_PROB_H_INCLUDED

#define zipf_cdf(a,b) R_DBL(FUNCNAME_ZIPFCDF(C_DTIT(a,b)))
#define zipf_pdf(a,b) R_DBL(FUNCNAME_ZIPFPDF(C_DTIT(a,b)))
#define buffon_laplace_pdf(a,b,c) R_DBL(FUNCNAME_BUFFONLAPLACEPDF(C_PDBL3(a,b,c)))
#define laplace_cdf(a,b,c) R_DBL(FUNCNAME_LAPLACECDF(C_PDBL3(a,b,c)))
#define laplace_cdf_inv(a,b,c) R_DBL(FUNCNAME_LAPLACECDFINV(C_PDBL3(a,b,c)))
#define laplace_pdf(a,b,c) R_DBL(FUNCNAME_LAPLACEPDF(C_PDBL3(a,b,c)))
#define log_normal_cdf(a,b,c) R_DBL(FUNCNAME_LOGNORMALCDF(C_PDBL3(a,b,c)))
#define normal_cdf(a,b,c) R_DBL(FUNCNAME_NORMALCDF(C_PDBL3(a,b,c)))
#define log_normal_pdf(a,b,c) R_DBL(FUNCNAME_LOGNORMALPDF(C_PDBL3(a,b,c)))
#define log_series_pdf(a,b) R_DBL(FUNCNAME_LOGSERIESPDF(C_DTIT(a,b)))
#define genlogistic_cdf(a,b,c,d) R_DBL(FUNCNAME_GENLOGISTICCDF(C_PDBL4(a,b,c,d)))
#define logistic_cdf(a,b,c) R_DBL(FUNCNAME_LOGISTICCDF(C_PDBL3(a,b,c)))
#define binomial_coef(a,b) R_USHRT(FUNCNAME_BINOMIALCOEF(C_PUSHRT2(a,b)))
#define negative_binomial_cdf(a,b,c) R_DBL(FUNCNAME_NEGATIVEBINOMIALCDF(C_2DTIT(a,b,c)))
#define negative_binomial_pdf(a,b,c) R_DBL(FUNCNAME_NEGATIVEBINOMIALPDF(C_2DTIT(a,b,c)))
#define poisson_pdf(a,b) R_DBL(FUNCNAME_POISSONPDF(C_DTIT(a,b)))
#define rayleigh_cdf(a,b) R_DBL(FUNCNAME_RAYLEIGHCDF(C_PDBL2(a,b)))
#define rayleigh_pdf(a,b) R_DBL(FUNCNAME_RAYLEIGHPDF(C_PDBL2(a,b)))
#define von_mises_pdf(a,b,c) R_DBL(FUNCNAME_VONMISESPDF(C_PDBL3(a,b,c)))
#define bessel_i0(a) R_DBL(FUNCNAME_BESSELI0(C_SDBL(a)))
#define weibull_cdf(a,b,c,d) R_DBL(FUNCNAME_WEIBULLCDF(C_PDBL4(a,b,c,d)))
#define weibull_pdf(a,b,c,d) R_DBL(FUNCNAME_WEIBULLPDF(C_PDBL4(a,b,c,d)))

__MATHSUITE __JBURKARDT void * FUNCNAME_BESSELI0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ZIPFCDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ZIPFPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOGSERIESPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BINOMIALCOEF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POISSONPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RAYLEIGHCDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RAYLEIGHPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BUFFONLAPLACEPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAPLACECDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAPLACECDFINV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAPLACEPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NORMALCDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOGNORMALCDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOGNORMALPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOGISTICCDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NEGATIVEBINOMIALCDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NEGATIVEBINOMIALPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VONMISESPDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GENLOGISTICCDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WEIBULLCDF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WEIBULLPDF(void *);

#endif
