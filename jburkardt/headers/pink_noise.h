#ifndef WRAPPER_PINK_NOISE_H_INCLUDED
#define WRAPPER_PINK_NOISE_H_INCLUDED

#define cdelay2(a,b) FUNCNAME_CDELAY2(C_IPI(a,b))
#define corr(a,b,c) FUNCNAME_CORR(C_2DTPIT(a,c,b))
#define cross_corr(a,b,c,d) FUNCNAME_CROSSCORR(C_2DT2PIT(a,d,b,c))
#define ran1f(a,b,c) R_DBL(FUNCNAME_RAN1F(C_DTPITPI(a,b,c)))
#define ranh(a,b,c) R_DBL(FUNCNAME_RANH(C_DTPITPI(a,b,c)))
#define wrap2(a,b) FUNCNAME_WRAP2(C_IPI(a,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_CDELAY2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RAN1F(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RANH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WRAP2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CORR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CROSSCORR(void *);

#endif
