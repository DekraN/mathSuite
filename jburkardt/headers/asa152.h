#ifndef WRAPPER_ASA152_H_INCLUDED
#define WRAPPER_ASA152_H_INCLUDED

#define alnfac(a) R_DBL(FUNCNAME_ALNFAC(C_SUSHRT(a)))
#define chyper(a,b,c,d) R_DBL(FUNCNAME_CHYPER(C_B3DT(a,b,c,d)))
#define hypergeometric_cdf_values(a,b,c,d,e,f) FUNCNAME_HYPERGEOMETRICCDFVALUES(C_5PDTPIT(a,b,c,d,e,f))
#define hypergeometric_pdf_values(a,b,c,d,e,f) FUNCNAME_HYPERGEOMETRICPDFVALUES(C_5PDTPIT(a,b,c,d,e,f))

__MATHSUITE __JBURKARDT void * FUNCNAME_ALNFAC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHYPER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HYPERGEOMETRICCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HYPERGEOMETRICPDFVALUES(void *);

#endif
