#ifndef WRAPPER_PWL_INTERP_2D_H_INCLUDED
#define WRAPPER_PWL_INTERP_2D_H_INCLUDED

#define r8vec_bracket5(a,b,c) R_INT(FUNCNAME_R8VECBRACKET5(C_DTPITIT(a,b,c)))
#define pwl_interp_2d(a,b,c,d,e,f,g,h) FUNCNAME_PWLINTERP2D(C_2DT3PITDT2PIT(a,b,c,d,e,f,g,h))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECBRACKET5(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PWLINTERP2D(void *);

#endif
