#ifndef WRAPPER_ASA005_H_INCLUDED
#define WRAPPER_ASA005_H_INCLUDED

#define alnorm(a,b) R_DBL(FUNCNAME_ALNORM(C_ITB(a,b)))
#define prncst(a,b,c) R_DBL(FUNCNAME_PRNCST(C_2ITDT(a,c,b)))
#define tfn(a,b) R_DBL(FUNCNAME_TFN(C_PDBL2(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_ALNORM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TFN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PRNCST(void *);

#endif
