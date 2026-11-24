#ifndef WRAPPER_PWL_INTERP_1D_H_INCLUDED
#define WRAPPER_PWL_INTERP_1D_H_INCLUDED

#define pwl_interp_1d(a,b,c,d,e) FUNCNAME_PWLINTERP1D(C_2DT3PIT(a,d,b,c,e))
#define pwl_basis_1d(a,b,c,d) FUNCNAME_PWLBASIS1D(C_2DT2PIT(a,c,b,d))
#define pwl_value_1d(a,b,c,d,e) FUNCNAME_PWLVALUE1D(C_2DT3PIT(a,d,b,c,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_PWLINTERP1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PWLBASIS1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PWLVALUE1D(void *);

#endif
