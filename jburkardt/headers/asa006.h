#ifndef WRAPPER_ASA006_H_INCLUDED
#define WRAPPER_ASA006_H_INCLUDED

#define cholesky(a,b,c,d,e) R_UCHAR(FUNCNAME_CHOLESKY(C_PIT2IPITPDT(a,b,c,d,e)))
#define subchl(a,b,c,d,e,f,g) R_UCHAR(FUNCNAME_SUBCHL(C_PITDTPDTPITPDTDTPIT(a,b,c,d,e,f,g)))

__MATHSUITE __JBURKARDT void * FUNCNAME_CHOLESKY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBCHL(void *);

#endif
