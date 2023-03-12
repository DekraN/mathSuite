#ifndef WRAPPER_LINPACK_D_H_INCLUDED
#define WRAPPER_LINPACK_D_H_INCLUDED

#define dcopy(a,b,c,d,e) FUNCNAME_DCOPY(C_3DT2PIT(a,c,e,d,b))
#define dasum(a,b,c) R_DBL(FUNCNAME_DASUM(C_2DTPIT(a,c,b)))
#define idamax(a,b,c) R_INT(FUNCNAME_IDAMAX(C_2DTPIT(a,c,b)))
#define dchdc(a,b,c,d,e,f) R_USHRT(FUNCNAME_DCHDC(C_PIT2DTPITPDTDT(a,b,c,d,e,f)))
#define dchdd(a,b,c,d,e,f,g,h,i,j,k) R_SHRT(FUNCNAME_DCHDD(C_PIT2DT2PIT2DT4PIT(a,b,c,d,e,f,g,h,i,j,k)))
#define dchex(a,b,c,d,e,f,g,h,i,j,k) FUNCNAME_DCHEX(C_PIT4DTPIT2DT2PITDT(a,b,c,d,e,f,g,h,i,j,k))
#define dchud(a,b,c,d,e,f,g,h,i,j,k) FUNCNAME_DCHUD(C_PIT2DT2PIT2DT4PIT(a,b,c,d,e,f,g,h,i,j,k))
#define dgbco(a,b,c,d,e,f,g) R_DBL(FUNCNAME_DGBCO(C_PIT4DTPDTPIT(a,b,c,d,e,f,g)))
#define dgbdi(a,b,c,d,e,f,g) FUNCNAME_DGBDI(C_PIT4DTPIPIT(a,b,c,d,e,f,g))
#define dgbfa(a,b,c,d,e,f) R_USHRT(FUNCNAME_DGBFA(C_PIT4DTPDT(a,b,c,d,e,f)))
#define dgbsl(a,b,c,d,e,f,g,h) FUNCNAME_DGBSL(C_PIT4DTPDTPITDT(a,b,c,d,e,f,g,h))
#define dgeco(a,b,c,d,e) R_DBL(FUNCNAME_DGECO(C_PIT2DTPDTPIT(a,b,c,d,e)))
#define dgedi(a,b,c,d,e,f,g,h) FUNCNAME_DGEDI(C_PIT2DTPI2PITDT(a,b,c,d,e,f,g,h))
#define dgefa(a,b,c,d) R_USHRT(FUNCNAME_DGEFA(C_PIT2DTPDT(a,b,c,d)))
#define dgesl(a,b,c,d,e,f) FUNCNAME_DGESL(C_PIT2DTPDTPITDT(a,b,c,e,d,f))
#define dgtsl(a,b,c,d,e) R_USHRT(FUNCNAME_DGTSL(C_DT4PIT(a,b,c,d,e)))
#define dpbco(a,b,c,d,e) R_DBL(FUNCNAME_DPBCO(C_3DT2PIT(b,c,d,e,a)))
#define dpbdi(a,b,c,d,e) FUNCNAME_DPBDI(C_3DT2PIT(b,c,d,e,a)))
#define dpbfa(a,b,c,d) R_USHRT(FUNCNAME_DPBFA(C_3DTPIT(b,c,d,a)))
#define dpbsl(a,b,c,d,e) FUNCNAME_DPBSL(C_3DT2PIT(b,c,d,e,a)))
#define dpoco(a,b,c,d) R_DBL(FUNCNAME_DPOCO(C_2DT2PIT(b,c,a,d)))
#define dpodi(a,b,c,d,e) FUNCNAME_DPODI(C_PIT2DTPITDT(b,c,e,d,a))
#define dpofa(a,b,c) R_USHRT(FUNCNAME_DPOFA(C_2DTPIT(b,c,a)))
#define dposl(a,b,c,d) FUNCNAME_DPOSL(C_2DT2PIT(b,c,a,d))
#define dppco(a,b,c) R_DBL(FUNCNAME_DPPCO(C_DT2PIT(b,c,a)))
#define dppdi(a,b,c,d) FUNCNAME_DPPDI(C_2DT2PIT(b,d,a,c))
#define dppfa(a,b) R_USHRT(FUNCNAME_DPPFA(C_DTPIT(b,a)))
#define dppsl(a,b,c) FUNCNAME_DPPSL(C_DT2PIT(b,c,a)))
#define dptsl(a,b,c,d) FUNCNAME_DPTSL(C_DT3PIT(a,b,c,d))
#define dqrdc(a,b,c,d,e,f,g,h) FUNCNAME_DQRDC(C_PIT3DTPITPIPITDT(a,b,c,d,e,f,g,h))
#define dqrsl(a,b,c,d,e,f,g,h,i,j,k,l) R_USHRT(FUNCNAME_DQRSL(C_PIT3DT7PITDT(a,b,c,d,e,f,g,h,i,j,k,l)))
#define dsico(a,b,c,d,e,f,g) R_DBL(FUNCNAME_DSICO(C_PIT2DTPDTPIT(a,b,c,d,e,f,g)))
#define dsidi(a,b,c,d,e,f,g,h) FUNCNAME_DSIDI(C_DT2PITDT2PIDTPIT(b,e,g,c,d,f,h,a))
#define dsifa(a,b,c,d) R_USHRT(FUNCNAME_DSIFA(C_PIT2DTPDT(a,b,c,d)))
#define dsisl(a,b,c,d,e) FUNCNAME_DSISL(C_DTPITDTPIPIT(b,a,c,d,e))
#define dspco(a,b,c,d) R_DBL(FUNCNAME_DSPCO(C_DT2PITPDT(b,a,d,c)))
#define dspdi(a,b,c,d,e,f,g) FUNCNAME_DSPDI(C_2DT2PI3PIT(b,g,c,e,d,f,a))
#define dspfa(a,b,c) R_USHRT(FUNCNAME_DSPFA(C_DTPITPDT(b,a,c)))
#define dspsl(a,b,c,d) FUNCNAME_DSPSL(C_PITDTPIPIT(a,b,c,d))
#define dsvdc(a,b,c,d,e,f,g,h,i,j,k,l) R_USHRT(FUNCNAME_DSVDC(C_PIT3DT3PITDTPITDTPITDT(a,b,c,d,e,f,g,h,i,j,k,l)))
#define dtrco(a,b,c,d,e) R_DBL(FUNCNAME_DTRCO(C_PIT2DTPITDT(b,c,e,d,a))
#define dtrsl(a,b,c,d,e) R_USHRT(FUNCNAME_DTRSL(C_PIT2DTPITDT(b,c,e,d,a))

__MATHSUITE __JBURKARDT void * FUNCNAME_DCOPY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DASUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_IDAMAX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DCHDC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DCHDD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DCHEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DCHUD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGBCO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGBDI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGBFA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGBSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGECO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGEDI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGEFA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGESL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGTSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPBCO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPBDI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPBFA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPBSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPOCO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPODI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPOFA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPOSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPPCO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPPDI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPPFA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPPSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPTSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DQRDC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DQRSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSICO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSIDI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSIFA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSISL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSPCO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSPDI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSPFA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSPSL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSVDC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DTRCO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DTRDI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DTRSL(void *);

#endif
