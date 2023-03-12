#ifndef WRAPPER_COMBO_H_INCLUDED
#define WRAPPER_COMBO_H_INCLUDED

#define bell_values(a,b,c) FUNCNAME_BELLVALUES(C_PPUSHRT3(a,b,c))
#define i4vec_search_binary_a(a,b,c) R_INT(FUNCNAME_I4VECSEARCHBINARYA(C_DTPII(a,b,c)))
#define i4vec_search_binary_d(a,b,c) R_INT(FUNCNAME_I4VECSEARCHBINARYD(C_DTPII(a,b,c)))
#define i4vec_sort_insert_a(a,b) FUNCNAME_I4VECSORTINSERTA(C_DTPI(a,b))
#define i4vec_sort_insert_d(a,b) FUNCNAME_I4VECSORTINSERTD(C_DTPI(a,b))
#define i4vec_backtrack(a,b,c,d,e,f,g,h) FUNCNAME_I4VECBACKTRACK(C_2DT6PI(a,b,c,d,e,f,g,h))
#define i4vec_sum(a,b) R_INT(FUNCNAME_I4VECSUM(C_DTPI(a,b)))
#define knapsack_01(a,b,c) FUNCNAME_KNAPSACK01(C_2DTPI(a,c,b))
#define perm_check(a,b) R_UCHR(FUNCNAME_PERMCHECK(C_DTPI(a,b)))
#define perm_inv(a,b) FUNCNAME_PERMINV(C_DTPI(a,b))
#define perm_lex_rank(a,b) R_USHRT(FUNCNAME_PERMLEXRANK(C_DTPI(a,b)))
#define perm_lex_unrank(a,b) FUNCNAME_PERMLEXUNRANK(C_PUSHRT2(a,b))



__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECSEARCHBINARYA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECSEARCHBINARYD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECSORTINSERTA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECSORTINSERTD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECBACKTRACK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECSUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERMCHECK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERMINV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERMLEXRANK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_KNAPSACK01(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERMLEXUNRANK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BELLVALUES(void *);


#endif
