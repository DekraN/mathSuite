#ifndef WRAPPER_FILON_H_INCLUDED
#define WRAPPER_FILON_H_INCLUDED

#define filon_tab_cos(a,b,c,d,e) R_DBL(FUNCNAME_FILONTABCOS(C_DT3ITPIT(a,c,d,e,b)))
#define filon_tab_sin(a,b,c,d,e) R_DBL(FUNCNAME_FILONTABSIN(C_DT3ITPIT(a,c,d,e,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_FILONTABCOS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FILONTABSIN(void *);

#endif
