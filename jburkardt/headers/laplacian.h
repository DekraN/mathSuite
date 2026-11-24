#ifndef WRAPPER_LAPLACIAN_H_INCLUDED
#define WRAPPER_LAPLACIAN_H_INCLUDED

#define cholesky_upper_error(a,b,c) R_DBL(FUNCNAME_CHOLESKYUPPERERROR(C_DT2PIT(a,b,c)))
#define r8mat_mtm_new(a,b,c,d,e) FUNCNAME_R8MATMTMNEW(C_3DT2PIT(a,b,c,d,e))
#define eigen_error(a,b,c,d,e) R_DBL(FUNCNAME_EIGENERROR(C_2DT3PIT(a,b,c,d,e)))
#define l1dd_apply(a,b,c) FUNCNAME_L1DDAPPLY(C_DTPITIT(a,c,b))
#define l1dd_cholesky(a,b) FUNCNAME_L1DDCHOLESKY(C_DTIT(a,b))
#define l1dd_eigen(a,b,c,d) FUNCNAME_L1DDEIGEN(C_DTIT2PIT(a,b,c,d))
#define l1dd(a,b) FUNCNAME_L1DD(C_DTIT(a,b))
#define l1dd_inverse(a,b) FUNCNAME_L1DDINVERSE(C_DTIT(a,b))
#define l1dd_lu(a,b,c,d) FUNCNAME_L1DDLU(C_DTIT2PIT(a,b,c,d))
#define l1dn_apply(a,b,c) FUNCNAME_L1DNAPPLY(C_DTPITIT(a,c,b))
#define l1dn_cholesky(a,b) FUNCNAME_L1DNCHOLESKY(C_DTIT(a,b))
#define l1dn_eigen(a,b,c,d) FUNCNAME_L1DNEIGEN(C_DTIT2PIT(a,b,c,d))
#define l1dn(a,b) FUNCNAME_L1DN(C_DTIT(a,b))
#define l1dn_inverse(a,b) FUNCNAME_L1DNINVERSE(C_DTIT(a,b))
#define l1dn_lu(a,b,c,d) FUNCNAME_L1DNLU(C_DTIT2PIT(a,b,c,d))
#define l1nd_apply(a,b,c) FUNCNAME_L1NDAPPLY(C_DTPITIT(a,c,b))
#define l1nd_cholesky(a,b) FUNCNAME_L1NDCHOLESKY(C_DTIT(a,b))
#define l1nd_eigen(a,b,c,d) FUNCNAME_L1NDEIGEN(C_DTIT2PIT(a,b,c,d))
#define l1nd(a,b) FUNCNAME_L1ND(C_DTIT(a,b))
#define l1nd_inverse(a,b) FUNCNAME_L1NDINVERSE(C_DTIT(a,b))
#define l1nd_lu(a,b,c,d) FUNCNAME_L1NDLU(C_DTIT2PIT(a,b,c,d))
#define l1nn_apply(a,b,c) FUNCNAME_L1NNAPPLY(C_DTPITIT(a,c,b))
#define l1nn_cholesky(a,b) FUNCNAME_L1NNCHOLESKY(C_DTIT(a,b))
#define l1nn_eigen(a,b,c,d) FUNCNAME_L1NNEIGEN(C_DTIT2PIT(a,b,c,d))
#define l1nn(a,b) FUNCNAME_L1NN(C_DTIT(a,b))
#define l1nn_lu(a,b,c,d) FUNCNAME_L1NNLU(C_DTIT2PIT(a,b,c,d))
#define l1pp_apply(a,b,c) FUNCNAME_L1PPAPPLY(C_DTPITIT(a,c,b))
#define l1pp_cholesky(a,b) FUNCNAME_L1PPCHOLESKY(C_DTIT(a,b))
#define l1pp_eigen(a,b,c,d) FUNCNAME_L1PPEIGEN(C_DTIT2PIT(a,b,c,d))
#define l1pp(a,b) FUNCNAME_L1PP(C_DTIT(a,b))
#define l1pp_lu(a,b,c,d) FUNCNAME_L1PPLU(C_DTIT2PIT(a,b,c,d))

__MATHSUITE __JBURKARDT void * FUNCNAME_CHOLESKYUPPERERROR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EIGENERROR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DDEIGEN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DDLU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DNEIGEN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DNLU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NDEIGEN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NDLU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NNEIGEN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NNLU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1PPEIGEN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1PPLU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMTMNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DDAPPLY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DDCHOLESKY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DDINVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DNAPPLY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DNCHOLESKY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1DNINVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NDAPPLY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NDCHOLESKY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1ND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NDINVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NNAPPLY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NNCHOLESKY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1NN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1PPAPPLY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1PPCHOLESKY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_L1PP(void *);

#endif
