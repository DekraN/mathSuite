#ifndef WRAPPER_ASA266_H_INCLUDED
#define WRAPPER_ASA266_H_INCLUDED

#define alogam(a) R_DBL(FUNCNAME_ALOGAM(C_SDBL(a)))
#define gamain(a,b) R_DBL(FUNCNAME_GAMAIN(C_PDBL2(a,b)))
#define r8col_mean(a,b,c) FUNCNAME_R8COLMEAN(C_2DTPIT(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLMEAN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ALOGAM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAMAIN(void *);

#endif
