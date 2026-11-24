#ifndef WRAPPER_CPV_H_INCLUDED
#define WRAPPER_CPV_H_INCLUDED

#define cpv(a,b,c,d) R_DBL(FUNCNAME_CPV(C_FIT2ITDT(a,b,c,d)))

__MATHSUITE __JBURKARDT void * FUNCNAME_CPV(void *);

#endif
