#ifndef WRAPPER_MATRIX_EXPONENTIAL_H_INCLUDED
#define WRAPPER_MATRIX_EXPONENTIAL_H_INCLUDED

#define r8mat_norm_li(a,b,c) R_DBL(FUNCNAME_R8MATNORMLI(C_2DTPIT(a,b,c)))
#define r8mat_fss_new(a,b,c,d) FUNCNAME_R8MATFSSNEW(C_2DT2PIT(a,c,b,d))
#define c8mat_copy(a,b,c,d) FUNCNAME_C8MATCOPYNEW(C_DTPCXDTPCX(a,c,b,d))
#define c8mat_fss(a,b,c,d) FUNCNAME_C8MATFSS(C_DTPCXDTPCX(a,b,c,d))
#define c8mat_fss_new(a,b,c,d) FUNCNAME_C8MATFSSNEW(C_DTPCXDTPCX(a,b,c,d))
#define r8mat_significant(a,b,c,d) R_UCHR(FUNCNAME_R8MATSIGNIFICANT(C_2DT2PIT(a,b,c,d)))
#define r8mat_minvm(a,b,c,d,e) FUNCNAME_R8MATMINVM(C_2DT3PIT(a,b,c,d,e))
#define r8mat_add(a,b,c,d,e,f,g) FUNCNAME_R8MATADD(C_2ITDTPITDT2PIT(c,e,a,d,b,f,g))
#define r8mat_scale(a,b,c,d) FUNCNAME_R8MATSCALE(C_2DTPITIT(a,b,d,c))
#define c8mat_minvm(a,b,c,d,e) FUNCNAME_C8MATMINVM(C_2DT3PCX(a,b,c,d,e))
#define c8mat_mm(a,b,c,d,e,f) FUNCNAME_C8MATMM(C_3DT3PCX(a,b,c,d,e,f))
#define c8mat_add_r8(a,b,c,d,e,f,g) FUNCNAME_C8MATADDR8(C_2DTITPCXIT2PCX(a,b,c,d,e,f,g))
#define c8mat_identity_new(a) FUNCNAME_C8MATIDENTITYNEW(C_SUSHRT(a))
#define c8mat_scale_r8(a,b,c,d) FUNCNAME_C8MATSCALER8(C_2DTITPCX(a,b,c,d))
#define c8mat_norm_li(a,b,c) R_DBL(FUNCNAME_C8MATNORMLI(C_2DTPCX(a,b,c)))
#define c8mat_copy_new(a,b,c) FUNCNAME_C8MATCOPYNEW(C_2DTPCX(a,b,c))
#define c8mat_expm1(a,b) FUNCNAME_C8MATEXPM1(C_DTPCX(a,b))
#define r8mat_expm1(a,b) FUNCNAME_R8MATEXPM1(C_DTPIT(a,b))
#define r8mat_expm2(a,b) FUNCNAME_R8MATEXPM2(C_DTPIT(a,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATNORMLI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATSIGNIFICANT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMINVM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATMINVM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATMM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATADD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATSCALE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATFSSNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATIDENTITYNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATCOPY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATCOPYNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATFSS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATFSSNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATADDR8(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATSCALER8(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATNORMLI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATEXPM1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATEXPM2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8MATEXPM1(void *);

#endif
