#ifndef WRAPPER_STOCHASTIC_RK_H_INCLUDED
#define WRAPPER_STOCHASTIC_RK_H_INCLUDED

#define rk1_ti_step(a,b,c,d,e,f,g) R_DBL(FUNCNAME_RK1TISTEP(C_4IT2FITPI(a,b,c,d,e,f,g)))
#define rk2_ti_step(a,b,c,d,e,f,g) R_DBL(FUNCNAME_RK2TISTEP(C_4IT2FITPI(a,b,c,d,e,f,g)))
#define rk3_ti_step(a,b,c,d,e,f,g) R_DBL(FUNCNAME_RK3TISTEP(C_4IT2FITPI(a,b,c,d,e,f,g)))
#define rk4_ti_step(a,b,c,d,e,f,g) R_DBL(FUNCNAME_RK4TISTEP(C_4IT2FITPI(a,b,c,d,e,f,g)))

__MATHSUITE __JBURKARDT void * FUNCNAME_RK1TISTEP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RK2TISTEP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RK3TISTEP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RK4TISTEP(void *);

#endif
