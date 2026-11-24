#ifndef WRAPPER_KRONROD_H_INCLUDED
#define WRAPPER_KRONROD_H_INCLUDED

#define abwe1(a,b,c,d,e,f,g,h) FUNCNAME_ABWE1(C_2DT2ITI3PIT(a,b,c,d,e,f,g,h))
#define abwe2(a,b,c,d,e,f,g,h,i) FUNCNAME_ABWE2(C_2DT2ITI4PIT(a,b,c,d,e,f,g,h,i))

__MATHSUITE __JBURKARDT void * FUNCNAME_ABWE1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ABWE2(void *);

#endif
