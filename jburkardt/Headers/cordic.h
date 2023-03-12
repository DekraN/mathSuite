#ifndef WRAPPER_CORDIC_H_INCLUDED
#define WRAPPER_CORDIC_H_INCLUDED

#define angle_shift(a,b) R_DBL(FUNCNAME_ANGLESHIFT(C_PDBL2(a,b)))
#define acos_cordic(a,b) R_DBL(FUNCNAME_ACOSCORDIC(C_PDBL2(a,b)))
#define asin_cordic(a,b) R_DBL(FUNCNAME_ASINCORDIC(C_PDBL2(a,b)))
#define atan_cordic(a,b,c) R_DBL(FUNCNAME_ATANCORDIC(C_2ITDT(a,b,c)))
#define cbrt_cordic(a,b) R_DBL(FUNCNAME_CBRTCORDIC(C_DTIT(b,a)))
#define exp_cordic(a,b) R_DBL(FUNCNAME_EXPCORDIC(C_DTIT(b,a)))
#define ln_cordic(a,b) R_DBL(FUNCNAME_LNCORDIC(C_DTIT(b,a)))
#define sqrt_cordic(a,b) R_DBL(FUNCNAME_SQRTCORDIC(C_DTIT(b,a)))
#define tan_cordic(a,b) R_DBL(FUNCNAME_TANCORDIC(C_DTIT(b,a)))
#define acos_values(a,b,c) FUNCNAME_ACOSVALUES(C_PDT2PIT(a,b,c))
#define arcsin_values(a,b,c) FUNCNAME_ARCSINVALUES(C_PDT2PIT(a,b,c))
#define arctan_values(a,b,c) FUNCNAME_ARCTANVALUES(C_PDT2PIT(a,b,c))
#define cbrt_values(a,b,c) FUNCNAME_CBRTVALUES(C_PDT2PIT(a,b,c))
#define cos_values(a,b,c) FUNCNAME_COSVALUES(C_PDT2PIT(a,b,c))
#define exp_values(a,b,c) FUNCNAME_EXPVALUES(C_PDT2PIT(a,b,c))
#define ln_values(a,b,c) FUNCNAME_LNVALUES(C_PDT2PIT(a,b,c))
#define sin_values(a,b,c) FUNCNAME_SINVALUES(C_PDT2PIT(a,b,c))
#define sqrt_values(a,b,c) FUNCNAME_SQRTVALUES(C_PDT2PIT(a,b,c))
#define tan_values(a,b,c) FUNCNAME_TANVALUES(C_PDT2PIT(a,b,c))


__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLESHIFT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ACOSCORDIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ASINCORDIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CBRTCORDIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EXPCORDIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LNCORDIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SQRTCORDIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TANCORDIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ATANCORDIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ARCSINVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ARCTANVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ACOSVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CBRTVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_COSVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EXPVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LNVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SINVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SQRTVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TANVALUES(void *);

#endif
