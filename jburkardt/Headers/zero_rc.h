#ifndef WRAPPER_ZERO_RC_H_INCLUDED
#define WRAPPER_ZERO_RC_H_INCLUDED

#define r8mat_fs(a,b,c,d) FUNCNAME_R8MATFS(C_2DT2PIT(a,b,c,d))
#define root_rc(a,b,c,d,e) R_DBL(FUNCNAME_ROOTRC(C_2IT3PIT(a,b,c,d,e)))
#define roots_rc(a,b,c,d,e,f) FUNCNAME_ROOTSRC(C_DT5PIT(a,b,c,d,e,f))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATFS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROOTRC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROOTSRC(void *);

#endif
