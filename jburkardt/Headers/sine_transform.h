#ifndef WRAPPER_SINE_TRANSFORM_H_INCLUDED
#define WRAPPER_SINE_TRANSFORM_H_INCLUDED

#define sine_transform_data(a,b) FUNCNAME_SINETRANSFORMDATA(C_DTPIT(a,b))
#define sine_transform_function(a,b,c,d) FUNCNAME_SINETRANSFORMFUNCTION(C_DT2ITFIT(a,b,c,d))
#define sine_transform_interpolant(a,b,c,d,e,f,g,h) FUNCNAME_SINETRANSFORMINTERPOLANT(C_2DT4IT2PIT(a,g,b,c,d,e,f,h))

__MATHSUITE __JBURKARDT void * FUNCNAME_SINETRANSFORMDATA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SINETRANSFORMINTERPOLANT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SINETRANSFORMFUNCTION(void *);

#endif

