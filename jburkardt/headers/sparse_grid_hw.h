#ifndef WRAPPER_SPARSE_GRID_HW_H_INCLUDED
#define WRAPPER_SPARSE_GRID_HW_H_INCLUDED

#define r8cvv_offset(a,b) FUNCNAME_R8CVVOFFSET(C_DTPI(a,b))
#define r8cvv_rget_new(a,b,c,d,e) FUNCNAME_R8CVVRGETNEW(C_2DTPIIPIT(a,c,d,e,b))
#define r8cvv_rset(a,b,c,d,e,f) FUNCNAME_R8CVVRSET(C_DTPITDTPIIPIT(a,b,c,d,e,f))
#define cce_order(a) R_USHRT(FUNCNAME_CCEORDER(C_SUSHRT(a)))
#define ccl_order(a) R_USHRT(FUNCNAME_CCLORDER(C_SUSHRT(a)))
#define ccs_order(a) R_USHRT(FUNCNAME_CCSORDER(C_SUSHRT(a)))
#define fn_integral(a) R_DBL(FUNCNAME_FNINTEGRAL(C_SUSHRT(a)))
#define fn_value(a,b,c) FUNCNAME_FNVALUE(C_2DTPIT(a,b,c))
#define fu_integral(a) R_DBL(FUNCNAME_FUINTEGRAL(C_SUSHRT(a)))
#define fu_value(a,b,c) FUNCNAME_FUVALUE(C_2DTPIT(a,b,c))
#define get_seq(a,b,c) FUNCNAME_GETSEQ(C_PUSHRT3(a,b,c))
#define gqn(a,b,c) FUNCNAME_GQN(C_DT2PIT(a,b,c))
#define gqn_order(a) R_USHRT(FUNCNAME_GQNORDER(C_SUSHRT(a)))
#define gqn2_order(a) R_USHRT(FUNCNAME_GQN2ORDER(C_SUSHRT(a)))
#define gqu(a,b,c) FUNCNAME_GQU(C_DT2PIT(a,b,c))
#define gqu_order(a) R_USHRT(FUNCNAME_GQUORDER(C_SUSHRT(a)))
#define kpn(a,b,c) FUNCNAME_KPN(C_DT2PIT(a,b,c))
#define kpu(a,b,c) FUNCNAME_KPU(C_DT2PIT(a,b,c))
#define kpu_order(a) R_USHRT(FUNCNAME_KPUORDER(C_SUSHRT(a)))
#define num_seq(a,b) R_USHRT(FUNCNAME_NUMSEQ(C_PINT2(a,b)))
#define rule_adjust(a,b,c,d,e,f,g) FUNCNAME_RULEADJUST(C_4ITDT2PIT(a,b,c,d,e,f,g))
#define rule_sort(a,b,c,d) FUNCNAME_RULESORT(C_2DT2PIT(a,b,c,d))
#define symmetric_sparse_size(a,b,c,d) R_USHRT(FUNCNAME_SYMMETRICSPARSESIZE(C_2DTPITIT(a,b,c,d)))
#define tensor_product(a,b,c,d,e,f,g,h) FUNCNAME_TENSORPRODUCT(C_DTPIT2DTPI3PIT(a,d,c,f,b,e,g,h))
#define tensor_product_cell(a,b,c,d,e,f,g,h,i) FUNCNAME_TENSORPRODUCTCELL(C_DT2PITDT2PIDT2PIT(a,b,c,d,e,f,g,h,i))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8CVVRSET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GQN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GQU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_KPN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_KPU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RULEADJUST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RULESORT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SYMMETRICSPARSESIZE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TENSORPRODUCT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TENSORPRODUCTCELL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CVVOFFSET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CVVRGETNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FNVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FUVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GETSEQ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CCEORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CCLORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CCSORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FNINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FUINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GQNORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GQN2ORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GQUORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_KPNORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_KPUORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NUMSEQ(void *);

#endif
