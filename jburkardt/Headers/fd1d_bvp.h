#ifndef WRAPPER_FD1D_BVP_H_INCLUDED
#define WRAPPER_FD1D_BVP_H_INCLUDED

#define fd1d_bvp(a,b,c,d,e,f) FUNCNAME_FD1DBVP(C_DT4FITPIT(a,b,c,d,e,f))
#define r83np_fs(a,b,c) FUNCNAME_R83NPFS(C_DT2PIT(a,b,c))
#define r8vec_even(a,b,c,d) FUNCNAME_R8VECEVEN(C_DT4FITPIT(a,b,c,d))


__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECEVEN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R83NPFS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FD1DBVP(void *);

#endif
