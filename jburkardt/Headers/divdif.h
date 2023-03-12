#ifndef WRAPPER_DIVDIF_H_INCLUDED
#define WRAPPER_DIVDIF_H_INCLUDED

#define cheby_t_zero(a) FUNCNAME_CHEBYTZERO(C_SUSHRT(a))
#define cheby_u_zero(a) FUNCNAME_CHEBYUZERO(C_SUSHRT(a))
#define data_to_dif(a,b,c,d) FUNCNAME_DATATODIF(C_DT3PIT(a,b,c,d))
#define data_to_r8poly(a,b,c,d) FUNCNAME_DATATOR8POLY(C_DT3PIT(a,b,c,d))
#define dif_antideriv(a,b,c,d,e,f) FUNCNAME_DIFANTIDERIV(C_DT2PITPDT2PIT(a,b,c,d,e,f))
#define dif_append(a,b,c,d,e,f,g,h) FUNCNAME_DIFAPPEND(C_DT2PIT2ITPDT2PIT(a,b,c,d,e,f,g,h))
#define dif_basis(a,b,c) FUNCNAME_DIFBASIS(C_DT2PIT(a,b,c))
#define dif_basis_deriv(a,b,c,d) FUNCNAME_DIFBASISDERIV(C_DT3PIT(a,b,c,d))
#define dif_basis_derivk(a,b,c,d,e) FUNCNAME_DIFBASISDERIVK(C_2DT3PIT(a,c,b,d,e))
#define dif_basis_i(a,b,c,d) FUNCNAME_DIFBASISI(C_2DT2PIT(a,b,c,d))
#define dif_deriv_table(a,b,c,d,e) FUNCNAME_DIFDERIVTABLE(C_DT4PIT(a,b,c,d,e))
#define dif_derivk_table(a,b,c,d,e,f) FUNCNAME_DIFDERIVKTABLE(C_DT2PITI2PIT(a,b,c,d,e,f))
#define dif_val(a,b,c,d) R_DBL(FUNCNAME_DIFVAL(C_DTIT2PIT(a,d,b,c)))
#define lagrange_rule(a,b) FUNCNAME_LAGRANGERULE(C_DTPIT(a,b))
#define lagrange_sum(a,b,c,d,e) R_DBL(FUNCNAME_LAGRANGESUM(C_DTIT3PIT(a,e,b,c,d)))
#define lagrange_val(a,b,c,d) R_DBL(FUNCNAME_LAGRANGEVAL(C_DTIT2PIT(a,d,b,c)))
#define nc_rule(a,b,c,d,e) FUNCNAME_NCRULE(C_DT2IT2PIT(a,b,c,d,e))
#define ncc_rule(a,b,c) FUNCNAME_NCCRULE(C_DT2PIT(a,b,c))
#define nco_rule(a,b,c) FUNCNAME_NCORULE(C_DT2PIT(a,b,c))
#define r8_swap(a,b) FUNCNAME_R8SWAP(C_PPDBL2(a,b))
#define r8poly_ant_cof(a,b,c) FUNCNAME_R8POLYANTCOF(C_DT2PIT(a,b,c))
#define r8poly_basis(a,b,c) FUNCNAME_R8POLYBASIS(C_DT2PIT(a,b,c))
#define r8poly_basis_1(a,b,c,d) FUNCNAME_R8POLYBASIS1(C_2DT2PIT(a,b,c,d))
#define r8poly_der_cof(a,b,c) FUNCNAME_R8POLYDERCOF(C_DT2PIT(a,b,c))
#define r8poly_der_val(a,b,c) R_DBL(FUNCNAME_R8POLYDERVAL(C_DTPITIT(a,b,c)))
#define r8poly_order(a,b) R_USHRT(FUNCNAME_R8POLYORDER(C_DTPIT(a,b)))
#define r8poly_shift(a,b,c,d) FUNCNAME_R8POLYSHIFT(C_DT2ITPIT(c,a,b,d))
#define r8poly_val_horner(a,b,c) R_DBL(FUNCNAME_R8POLYVALHORNER(C_DTPITIT(a,b,c)))
#define r8vec_indicator(a,b) FUNCNAME_R8VECINDICATOR(C_DTPIT(a,b))
#define roots_to_dif(a,b,c,d,e) FUNCNAME_ROOTSTODIF(C_DTPITPDT2PIT(a,b,c,d,e))
#define roots_to_r8poly(a,b,c,d) FUNCNAME_ROOTSTOR8POLY(C_PITDTPIPIT(b,a,c,d))





__MATHSUITE __JBURKARDT void * FUNCNAME_DATATODIF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DATATOR8POLY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFANTIDERIV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFAPPEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFBASIS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFBASISDERIV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFBASISDERIVK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFBASISI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFDERIVTABLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFDERIVKTABLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFVAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGRANGESUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGRANGEVAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NCRULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NCCRULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NCORULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SWAP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYANTCOF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYBASIS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYBASIS1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYDERCOF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYDERVAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYSHIFT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYVALHORNER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDICATOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROOTSTODIF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROOTSTOR8POLY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBYTZERO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBYUZERO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGRANGERULE(void *);

#endif
