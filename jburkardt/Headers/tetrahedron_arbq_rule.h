#ifndef WRAPPER_TETRAHEDRON_ARBQ_RULE_H_INCLUDED
#define WRAPPER_TETRAHEDRON_ARBQ_RULE_H_INCLUDED

#define r8mat_row_copy(a,b,c,d,e) FUNCNAME_R8MATROWCOPY(C_3DT2PIT(a,b,c,d,e))
#define ref_to_koorn(a) FUNCNAME_REFTOKOORN(a)
#define tetrahedron_rule01(a,b,c) FUNCNAME_TETRAHEDRONRULE01(C_DT2PIT(a,b,c))
#define tetrahedron_rule02(a,b,c) FUNCNAME_TETRAHEDRONRULE02(C_DT2PIT(a,b,c))
#define tetrahedron_rule03(a,b,c) FUNCNAME_TETRAHEDRONRULE03(C_DT2PIT(a,b,c))
#define tetrahedron_rule04(a,b,c) FUNCNAME_TETRAHEDRONRULE04(C_DT2PIT(a,b,c))
#define tetrahedron_rule05(a,b,c) FUNCNAME_TETRAHEDRONRULE05(C_DT2PIT(a,b,c))
#define tetrahedron_rule06(a,b,c) FUNCNAME_TETRAHEDRONRULE06(C_DT2PIT(a,b,c))
#define tetrahedron_rule07(a,b,c) FUNCNAME_TETRAHEDRONRULE07(C_DT2PIT(a,b,c))
#define tetrahedron_rule08(a,b,c) FUNCNAME_TETRAHEDRONRULE08(C_DT2PIT(a,b,c))
#define tetrahedron_rule09(a,b,c) FUNCNAME_TETRAHEDRONRULE09(C_DT2PIT(a,b,c))
#define tetrahedron_rule10(a,b,c) FUNCNAME_TETRAHEDRONRULE10(C_DT2PIT(a,b,c))
#define tetrahedron_rule11(a,b,c) FUNCNAME_TETRAHEDRONRULE11(C_DT2PIT(a,b,c))
#define tetrahedron_rule12(a,b,c) FUNCNAME_TETRAHEDRONRULE12(C_DT2PIT(a,b,c))
#define tetrahedron_rule13(a,b,c) FUNCNAME_TETRAHEDRONRULE13(C_DT2PIT(a,b,c))
#define tetrahedron_rule14(a,b,c) FUNCNAME_TETRAHEDRONRULE14(C_DT2PIT(a,b,c))
#define tetrahedron_rule15(a,b,c) FUNCNAME_TETRAHEDRONRULE15(C_DT2PIT(a,b,c))
#define tetrahedron_arbq(a,b,c,d) FUNCNAME_TETRAHEDRONARBQ(C_2DT2PIT(a,b,c,d))
#define tetrahedron_arbq_size(a) R_USHRT(FUNCNAME_TETRAHEDRONARBQSIZE(C_SUSHRT(a)))
#define tetrahedron_ref(a,b,c,d) R_USHRT(FUNCNAME_TETRAHEDRONREF(C_PPDBL4(a,b,c,d)))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATROWCOPY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE01(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE02(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE03(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE04(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE05(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE06(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE07(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE08(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE09(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE10(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE11(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE12(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE13(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE14(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRULE15(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONARBQ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONREF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_REFTOKOORN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONARBQSIZE(void *);

#endif