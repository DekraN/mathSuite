#ifndef WRAPPER_TOMS655_H_INCLUDED
#define WRAPPER_TOMS655_H_INCLUDED

#define scmm(a,b,c,d,e,f) FUNCNAME_SCMM(C_2DT4IT(a,b,c,d,e,f))
#define cawiq(a,b,c,d,e,f,g,h,i,j,k) FUNCNAME_CAWIQ(C_DTPITPIDTPI2DT2PITDTIT(a,b,c,d,e,f,g,h,i,j,k))
#define cgqfs(a,b,c,d,e,f,g) FUNCNAME_CGQFS(C_2DT2ITDT2PIT(a,b,c,d,e,f,g))
#define chkqf(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) FUNCNAME_CHKQF(C_2PITPI2DTPI4DT2ITDT2IT(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o))
#define chkqfs(a,b,c,d,e,f,g,h,i,j,k,l,m,n) FUNCNAME_CHKQFS(C_2PITPI2DTPIDTPIT3DT2ITDT(a,b,c,d,e,f,g,h,i,j,k,l,m,n))
#define ciqf(a,b,c,d,e,f,g,h,i,j,k,l) FUNCNAME_CIQF(C_DTPITPIDTPI2DT4ITDT(a,b,c,d,e,f,g,h,i,j,k,l))
#define ciqfs(a,b,c,d,e,f,g,h,i,j) FUNCNAME_CIQFS(C_DTPITPIDTPI2DT2ITDT(a,b,c,d,e,f,g,h,i,j))
#define cliqf(a,b,c,d,e,f,g,h) FUNCNAME_CLIQF(C_DTPITDT4ITDT(a,b,c,d,e,f,g,h))
#define cliqfs(a,b,c,d,e,f) FUNCNAME_CLIQFS(C_DTPITDT4ITDT(a,b,c,d,e,f))
#define cwiqd(a,b,c,d,e,f,g,h,i) FUNCNAME_CWIQD(C_3DTITPITDT3PIT(a,b,c,d,e,f,g,h,i))
#define sct(a,b,c,d,e) FUNCNAME_SCT(C_DT2ITDTPIT(c,d,e,a,b))
#define wm(a,b,c,d) FUNCNAME_WM(C_DT2ITDT(a,c,d,b))
#define wtfn(a,b,c,d,e) FUNCNAME_WTFN(C_DT2ITDTPIT(b,d,e,c,a))

__MATHSUITE __JBURKARDT void * FUNCNAME_CGQFS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHKQF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHKQFS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SCMM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CAWIQ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIQF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIQFS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CLIQF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CLIQFS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CWIQD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SCT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WTFN(void *);

__MATHSUITE __JBURKARDT  ityp   eiqf ( const register dim_typ nt, ityp [static nt], int [static nt], ityp [], const register dim_typ, int [static nt],dim_typ, ityp  ( ityp, dim_typ ) );
__MATHSUITE __JBURKARDT  ityp   eiqfs ( const register dim_typ nt, ityp [static nt], ityp [static nt], ityp ( ityp, dim_typ ) );

#endif
