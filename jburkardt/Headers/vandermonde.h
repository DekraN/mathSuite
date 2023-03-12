#ifndef WRAPPER_VANDERMONDE_H_INCLUDED
#define WRAPPER_VANDERMONDE_H_INCLUDED

#define bivand1(a,b,c) FUNCNAME_BIVAND1(C_DT2PIT(a,b,c))
#define bivand2(a,b,c) FUNCNAME_BIVAND2(C_DT2PIT(a,b,c))
#define dvand(a,b,c) FUNCNAME_DVAND(C_DT2PIT(a,b,c))
#define dvandprg(a,b,c,d,e,f) FUNCNAME_DVANDPRG(C_DT5PIT(a,b,c,d,e,f))
#define pvand(a,b,c) FUNCNAME_PVAND(C_DT2PIT(a,b,c))
#define pvandprg(a,b,c,d,e,f) FUNCNAME_PVANDPRG(C_DT5PIT(a,b,c,d,e,f))

__MATHSUITE __JBURKARDT void * FUNCNAME_DVANDPRG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PVANDPRG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BIVAND1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BIVAND2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DVAND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PVAND(void *);

#endif
