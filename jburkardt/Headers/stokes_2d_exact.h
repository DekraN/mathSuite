#ifndef WRAPPER_STOKES_2D_EXACT_H_INCLUDED
#define WRAPPER_STOKES_2D_EXACT_H_INCLUDED

#define resid_stokes1(a,b,c,d,e,f) FUNCNAME_RESIDSTOKES1(C_DT5PIT(a,b,c,d,e,f))
#define resid_stokes2(a,b,c,d,e,f) FUNCNAME_RESIDSTOKES2(C_DT5PIT(a,b,c,d,e,f))
#define resid_stokes3(a,b,c,d,e,f) FUNCNAME_RESIDSTOKES3(C_DT5PIT(a,b,c,d,e,f))
#define rhs_stokes1(a,b,c,d,e,f) FUNCNAME_RHSSTOKES1(C_DT5PIT(a,b,c,d,e,f))
#define rhs_stokes2(a,b,c,d,e,f) FUNCNAME_RHSSTOKES2(C_DT5PIT(a,b,c,d,e,f))
#define rhs_stokes3(a,b,c,d,e,f) FUNCNAME_RHSSTOKES3(C_DT5PIT(a,b,c,d,e,f))

__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDSTOKES1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDSTOKES2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDSTOKES3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RHSSTOKES1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RHSSTOKES2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RHSSTOKES3(void *);

#endif

