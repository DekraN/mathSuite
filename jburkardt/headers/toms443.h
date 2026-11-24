#ifndef WRAPPER_TOMS443_H_INCLUDED
#define WRAPPER_TOMS443_H_INCLUDED

#define wew_a(a,b) R_DBL(FUNCNAME_WEWA(C_ITPIT(a,b)))
#define wew_b(a,b) R_DBL(FUNCNAME_WEWB(C_ITPIT(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_WEWA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WEWB(void *);

#endif
