#ifndef WRAPPER_CHEBYSHEV_SERIES_H_INCLUDED
#define WRAPPER_CHEBYSHEV_SERIES_H_INCLUDED

#define echebser0(a,b,c) R_DBL(FUNCNAME_ECHEBSER0(C_DTPITIT(c,b,a)))
#define echebser1(a,b,c,d) R_DBL(FUNCNAME_ECHEBSER1(C_DTIT2PIT(c,a,b,d)))
#define echebser2(a,b,c,d,e) R_DBL(FUNCNAME_ECHEBSER2(C_DTIT3PIT(c,a,d,e,b)))
#define echebser3(a,b,c,d,e,f) R_DBL(FUNCNAME_ECHEBSER3(C_DT2PITIT2PIT(c,d,e,a,f,b)))
#define echebser4(a,b,c,d,e,f,g) R_DBL(FUNCNAME_ECHEBSER4(C_ITPITDT4PIT(a,b,c,d,e,f,g)))
#define evenchebser0(a,b,c) R_DBL(FUNCNAME_EVENCHEBSER0(C_DTPITIT(c,b,a)))
#define evenchebser1(a,b,c,d) R_DBL(FUNCNAME_EVENCHEBSER1(C_DTIT2PIT(c,a,b,d)))
#define evenchebser2(a,b,c,d,e) R_DBL(FUNCNAME_EVENCHEBSER2(C_DTIT3PIT(c,a,d,e,b)))
#define oddchebser0(a,b,c) R_DBL(FUNCNAME_ODDCHEBSER0(C_DTPITIT(c,b,a)))
#define oddchebser1(a,b,c,d) R_DBL(FUNCNAME_ODDCHEBSER1(C_DTIT2PIT(c,a,b,d)))
#define oddchebser2(a,b,c,d,e) R_DBL(FUNCNAME_ODDCHEBSER2(C_DTIT3PIT(c,a,d,e,b)))
#define a_to_i4(a) R_INT(FUNCNAME_ATOI4(C_SCHR(a)))

__MATHSUITE __JBURKARDT void * FUNCNAME_ECHEBSER0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ECHEBSER1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ECHEBSER2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ECHEBSER3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ECHEBSER4(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EVENCHEBSER0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EVENCHEBSER1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EVENCHEBSER2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ODDCHEBSER0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ODDCHEBSER1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ODDCHEBSER2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ATOI4(void *);

#endif
