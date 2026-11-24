#ifndef WRAPPER_ASA172_H_INCLUDED
#define WRAPPER_ASA172_H_INCLUDED

#define revers(a,b) FUNCNAME_REVERS(C_DTPIT(a,b))
#define simdo(a,b,c,d,e,f) R_UCHR(FUNCNAME_SIMDO(C_2BDTPITPDTPIT(a,b,c,d,e,f)))

__MATHSUITE __JBURKARDT void * FUNCNAME_REVERS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMDO(void *);

#endif
