#ifndef WRAPPER_WALSH_H_INCLUDED
#define WRAPPER_WALSH_H_INCLUDED

#define ffwt(a,b) FUNCNAME_FFWT(C_DTPIT(a,b))
#define fwt(a,b) FUNCNAME_FWT(C_DTPIT(a,b))
#define walsh(a,b) FUNCNAME_WALSH(C_DTPIT(a,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_FFWT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FWT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WALSH(void *);

#endif
