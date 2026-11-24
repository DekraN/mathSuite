#ifndef WRAPPER_ASA091_H_INCLUDED
#define WRAPPER_ASA091_H_INCLUDED

#define gammad(a,b) R_DBL(FUNCNAME_GAMMAD(C_PDBL2(a,b)))
#define ppchi2(a,b,c) R_DBL(FUNCNAME_PPCHI2(C_PDBL3(a,b,c)))

__MATHSUITE  void * FUNCNAME_GAMMAD(void *);
__MATHSUITE  void * FUNCNAME_PPCHI2(void *);

#endif
