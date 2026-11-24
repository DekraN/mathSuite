#ifndef WRAPPER_HERMITE_H_INCLUDED
#define WRAPPER_HERMITE_H_INCLUDED

#define dif_shift_x(a,b,c,d) FUNCNAME_DIFSHIFTX(C_DTIT2PIT(a,d,b,c))
#define dif_shift_zero(a,b,c) FUNCNAME_DIFSHIFTZERO(C_DT2PIT(a,b,c))
#define dif_to_r8poly(a,b,c,d) FUNCNAME_DIFTOR8POLY(C_DT3PIT(a,b,c,d))
#define dif_vals(a,b,c,d,e) FUNCNAME_DIFVALS(C_DTIT2PIT(a,d,b,c))
#define hermite_basis_0(a,b,c,d) R_DBL(FUNCNAME_HERMITEBASIS0(C_2DTPITIT(a,c,b,d)))
#define hermite_basis_1(a,b,c,d) R_DBL(FUNCNAME_HERMITEBASIS1(C_2DTPITIT(a,c,b,d)))
#define hermite_interpolant(a,b,c,d,e,f,g,h) R_DBL(FUNCNAME_HERMITEINTERPOLANT(C_DT7PIT(a,b,c,d,e,f,g,h)))
#define hermite_interpolant_rule(a,b,c,d) FUNCNAME_HERMITEINTERPOLANTRULE(C_DT2ITPIT(a,b,c,d))
#define hermite_interpolant_value(a,b,c,d,e,f,g,h,i) FUNCNAME_HERMITEINTERPOLANTVALUE(C_DT4PITDT3PIT(a,b,c,d,e,f,g,h,i))
#define r8poly_ant_val(a,b,c) R_DBL(FUNCNAME_R8POLYANTVAL(C_DTPITIT(a,b,c)))
#define r8poly_degree(a,b) R_DBL(FUNCNAME_R8POLYDEGREE(C_DTPIT(a,b)))
#define r8vec_chebyshev_new(a,b,c) FUNCNAME_R8VECCHEBYSHEVNEW(C_2ITDT(b,c,a))

__MATHSUITE __JBURKARDT void * FUNCNAME_DIFDERIV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFSHIFTX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFSHIFTZERO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFTOR8POLY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HERMITEBASIS0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HERMITEBASIS1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HERMITEINTERPOLANT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HERMITEINTERPOLANTVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYANTVAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYDEGREE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCOPY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFVALS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HERMITEINTERPOLANTRULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCHEBYSHEVNEW(void *);

#endif
