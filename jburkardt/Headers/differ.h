#ifndef WRAPPER_DIFFER_H_INCLUDED
#define WRAPPER_DIFFER_H_INCLUDED

#define differ_backward(a,b,c,d,e) FUNCNAME_DIFFERBACKWARD(C_IT2DT2PIT(a,b,c,d,e))
#define differ_central(a,b,c,d,e) FUNCNAME_DIFFERCENTRAL(C_IT2DT2PIT(a,b,c,d,e))
#define differ_forward(a,b,c,d,e) FUNCNAME_DIFFERFORWARD(C_IT2DT2PIT(a,b,c,d,e))
#define differ_inverse(a,b) FUNCNAME_DIFFERINVERSE(C_DTPIT(a,b))
#define differ_matrix(a,b) FUNCNAME_DIFFERMATRIX(C_DTPIT(a,b))
#define differ_solve(a,b,c) FUNCNAME_DIFFERSOLVE(C_2DTPIT(a,c,b))
#define differ_stencil(a,b,c,d,e) FUNCNAME_DIFFERSTENCIL(C_IT2DT2PIT(a,b,c,d,e))
#define r8_factorial(a) R_DBL(FUNCNAME_R8FACTORIAL(C_SUSHRT(a)))
#define r8mat_fs_new(a,b,c) FUNCNAME_R8MATFSNEW(C_DT2PIT(a,b,c))
#define r8mat_sub_new(a,b,c,d) FUNCNAME_R8MATSUBNEW(C_2DT2PIT(a,b,c,d))
#define r8vm_sl(a,b,c,d,e,f) FUNCNAME_R8VMSL(C_DT2PITDTPITPS(a,b,c,d,e,f))
#define r8vm_sl_new(a,b,c,d,e) FUNCNAME_R8VMSLNEW(C_DT2PITDTPS(a,b,c,d,e))
#define r8mat_mm_new(a,b,c,d,e) FUNCNAME_R8MATMMNEW(C_3DT2PIT(a,b,c,d,e))


__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFERBACKWARD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFERCENTRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFERFORWARD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFERSTENCIL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_INVERSEERROR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VMSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFERINVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFERMATRIX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFERSOLVE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATFSNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATSUBNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VMSLNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMMNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8FACTORIAL(void *);

#endif
