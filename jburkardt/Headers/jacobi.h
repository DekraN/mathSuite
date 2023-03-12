#ifndef WRAPPER_JACOBI_H_INCLUDED
#define WRAPPER_JACOBI_H_INCLUDED

#define dif2(a,b) FUNCNAME_DIF2(C_PUSHRT2(a,b))
#define r8mat_residual_norm(a,b,c,d,e) R_DBL(FUNCNAME_R8MATRESIDUALNORM(C_2DT3PIT(a,b,c,d,e)))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATRESIDUALNORM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIF2(void *);

#endif
