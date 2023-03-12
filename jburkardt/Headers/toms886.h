#ifndef WRAPPER_TOMS886_H_INCLUDED
#define WRAPPER_TOMS886_H_INCLUDED

#define dgemm(a,b,c,d,e,f,g,h,i,j,k,l,m) FUNCNAME_DGEMM(C_2C3DTITPITDTPITDTITPITDT(a,b,c,d,e,f,g,h,i,j,k,l,m))
#define cheb(a,b,c) FUNCNAME_CHEB(C_DTPITIT(a,c,b))
#define franke(a,b) R_DBL(FUNCNAME_FRANKE(C_PDBL2(a,b)))
#define padua2(a,b,c,d,e,f,g,h,i) FUNCNAME_PADUA2(C_3DT6PIT(a,b,c,d,e,f,g,h,i))
#define pd2val(a,b,c,d,e) R_DBL(FUNCNAME_PD2VAL(C_DT2ITDTPIT(a,d,e,b,c)))
#define pdpts(a,b,c,d,e) FUNCNAME_PDPTS(C_DT3PITPI(a,b,c,d,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_DGEMM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PADUA2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PDPTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FRANKE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PD2VAL(void *);


#endif
