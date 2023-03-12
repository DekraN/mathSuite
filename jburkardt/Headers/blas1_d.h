#ifndef WRAPPER_BLAS1_D_H_INCLUDED
#define WRAPPER_BLAS1_D_H_INCLUDED

#define daxpy(a,b,c,d,e,f) FUNCNAME_DAXPY(C_DTITPITDTPITDT(a,b,c,d,e,f))
#define ddot(a,b,c,d,e) R_DBL(FUNCNAME_DDOT(C_3DT2PIT(a,c,e,b,d)))
#define dnrm2(a,b,c) R_DBL(FUNCNAME_DNRM2(C_DTPITI(a,b,c)))
#define drot(a,b,c,d,e,f,g) FUNCNAME_DROT(C_DTPITIPITI2IT(a,b,c,d,e,f,g))
#define drotg(a,b,c,d) FUNCNAME_DROTG(C_PPDBL4(a,b,c,d))
#define dscal(a,b,c,d) FUNCNAME_DSCAL(C_DTITPITI(a,b,c,d))
#define dswap(a,b,c,d,e) FUNCNAME_DSWAP(C_DTPITIPITI(a,b,c,d,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_DAXPY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DDOT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DNRM2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DROT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DROTG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSCAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSWAP(void *);

#endif
