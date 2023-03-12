#ifndef WRAPPER_C_CALLS_F77_H_INCLUDED
#define WRAPPER_C_CALLS_F77_H_INCLUDED

#define kronrod(a,b,c,d,e) FUNCNAME_KRONROD(C_DTIT3PIT(a,b,c,d,e))
#define kronrod_adjust(a,b,c,d,e,f) FUNCNAME_KRONRODADJUST(C_2ITDT3PIT(a,b,c,d,e,f))

__MATHSUITE __JBURKARDT void * FUNCNAME_KRONROD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_KRONRODADJUST(void *);

#endif
