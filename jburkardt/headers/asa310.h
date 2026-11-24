#ifndef WRAPPER_ASA310_H_INCLUDED
#define WRAPPER_ASA310_H_INCLUDED

#define ncbeta(a,b,c,d,e) R_DBL(FUNCNAME_NCBETA(C_PDBL5(a,b,c,d,e)))

__MATHSUITE __JBURKARDT void * FUNCNAME_NCBETA(void *);

#endif
