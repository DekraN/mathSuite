#ifndef WRAPPER_BESSEL_H_INCLUDED
#define WRAPPER_BESSEL_H_INCLUDED

#define bessj0(a) R_DBL(FUNCNAME_BESSJ0(C_SDBL(a)))
#define bessj(a,b) R_DBL(FUNCNAME_BESSJ(C_DTIT(a,b)))
#define bessj1(a) R_DBL(FUNCNAME_BESSJ1(C_SDBL(a)))
#define bessy0(a) R_DBL(FUNCNAME_BESSY0(C_SDBL(a)))
#define bessy1(a) R_DBL(FUNCNAME_BESSY1(C_SDBL(a)))
#define bessy(a,b) R_DBL(FUNCNAME_BESSY(C_DTIT(a,b)))
#define bessi0(a) R_DBL(FUNCNAME_BESSI0(C_SDBL(a)))
#define bessi1(a) R_DBL(FUNCNAME_BESSI1(C_SDBL(a)))
#define bessi(a,b) R_DBL(FUNCNAME_BESSI(C_DTIT(a,b)))
#define bessk0(a) R_DBL(FUNCNAME_BESSK0(C_SDBL(a)))
#define bessk1(a) R_DBL(FUNCNAME_BESSK1(C_SDBL(a)))
#define bessk(a,b) R_DBL(FUNCNAME_BESSK(C_DTIT(a,b)))


__MATHSUITE __JBURKARDT void * FUNCNAME_BESSJ0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSJ1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSY0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSY1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSI0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSI1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSK0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSK1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSJ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSK(void *);


#endif
