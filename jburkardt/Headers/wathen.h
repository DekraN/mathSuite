#ifndef WRAPPER_WATHEN_H_INCLUDED
#define WRAPPER_WATHEN_H_INCLUDED

#define mv_st(a,b,c,d,e,f,g) FUNCNAME_MVST(C_3DT2PI2PIT(a,b,c,d,e,f,g))
#define mv_gb(a,b,c,d,e,f) FUNCNAME_MVGB(C_4DT2PIT(a,b,c,d,e,f))
#define mv_ge(a,b,c,d) FUNCNAME_MVGE(C_2DT2PIT(a,b,c,d))
#define cg_st(a,b,c,d,e,f,g) FUNCNAME_CGST(C_2DT2PI3PIT(a,b,c,d,e,f,g))
#define wathen_bandwidth(a,b,c,d,e) FUNCNAME_WATHENBANDWIDTH(C_2DT3PDT(a,b,c,d,e))
#define wathen_gb(a,b,c,d) FUNCNAME_WATHENGB(C_3DTPI(a,b,c,d))
#define wathen_ge(a,b,c,d) FUNCNAME_WATHENGE(C_3DTPI(a,b,c,d))
#define wathen_order(a,b) R_USHRT(FUNCNAME_WATHENORDER(C_PUSHRT2(a,b)))
#define wathen_st(a,b,c,d,e,f) FUNCNAME_WATHENST(C_3DTPI2PDT(a,b,c,d,e,f))
#define wathen_st_size(a,b) R_USHRT(FUNCNAME_WATHENSTSIZE(C_PUSHRT2(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_CGST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WATHENBANDWIDTH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MVST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MVGB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MVGE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WATHENGB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WATHENGE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WATHENST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WATHENORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WATHENSTSIZE(void *);

#endif
