#ifndef WRAPPER_HAAR_H_INCLUDED
#define WRAPPER_HAAR_H_INCLUDED

#define haar_1d(a,b) FUNCNAME_HAAR1D(C_DTPIT(a,b))
#define haar_1d_inverse(a,b) FUNCNAME_HAAR1DINVERSE(C_DTPIT(a,b))
#define haar_2d(a,b,c) FUNCNAME_HAAR2D(C_2DTPIT(a,b,c))
#define haar_2d_inverse(a,b,c) FUNCNAME_HAAR2DINVERSE(C_2DTPIT(a,b,c))
#define dif_deriv(a,b,c,d,e,f) FUNCNAME_DIFDERIV(C_DT2PITPDT2PIT(a,b,c,d,e,f))
#define r8vec_copy_new(a,b) FUNCNAME_R8VECCOPYNEW(C_DTPIT(a,b))
#define r8vec_ones_new(a) FUNCNAME_R8VECONESNEW(a)

__MATHSUITE __JBURKARDT void * FUNCNAME_HAAR1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HAAR1DINVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HAAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HAAR2DINVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFDERIV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCOPYNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECONESNEW(void *);

#endif
