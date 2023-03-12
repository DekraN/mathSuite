#ifndef WRAPPER_SUBSET_SUM_H_INCLUDED
#define WRAPPER_SUBSET_SUM_H_INCLUDED

#define gen_laguerre_recur(a,b,c,d,e,f,g,h) FUNCNAME_GENLAGUERRERECUR(C_PITDT2IT4PIT(a,e,d,f,b,c,g,h))
#define hermite_recur(a,b,c,d,e) FUNCNAME_HERMITERECUR(C_DTIT3PIT(e,d,a,b,c))
#define jacobi_recur(a,b,c,d,e,f,g,h,i) FUNCNAME_JACOBIRECUR(C_2ITDT2PITIT3PIT(f,g,e,a,b,d,c,h,i))
#define laguerre_recur(a,b,c,d,e,f,g) FUNCNAME_LAGUERRERECUR(C_ITPITDT4PIT(d,a,e,f,g,b,c))
#define gen_laguerre_root(a,b,c,d,e,f,g) FUNCNAME_GENLAGUERREROOT(C_ITPITDT4PIT(c,a,b,d,e,f,g))
#define hermite_root(a,b,c,d) FUNCNAME_HERMITEROOT(C_DT3PIT(b,c,d,a))
#define jacobi_root(a,b,c,d,e,f,g,h) FUNCNAME_JACOBIROOT(C_PITDT2IT4PIT(a,b,c,d,e,f,g,h))
#define laguerre_root(a,b,c,d,e,f) FUNCNAME_LAGUERREROOT(C_DT5PIT(b,a,c,d,e,f))
#define perm0_cycle(a,b,c,d,e) FUNCNAME_PERM0CYCLE(C_DTPIPSPIB(a,b,c,d,e))
#define i4_sqrt_cf(a,b,c,d) FUNCNAME_I4SQRTCF(C_2DT2PI(a,b,c,d))
#define ubvec_to_ui4(a,b) R_UINT(FUNCNAME_UBVECTOUI4(C_DTPI(a,b)))
#define ubvec_xor(a,b,c,d) FUNCNAME_UBVECXOR(C_DT3PI(a,b,c,d))
#define ui4_to_ubvec(a,b,c) FUNCNAME_UI4TOUBVEC(C_DTPII(b,c,a))
#define i4_sqrt(a,b,c) FUNCNAME_I4SQRT(C_DT2PI(a,b,c))
#define triang(a,b,c) FUNCNAME_TRIANG(C_DT2PI(a,b,c))
#define perm0_inverse(a,b) FUNCNAME_PERM0INVERSE(C_DTPI(a,b))
#define i4mat_2perm0(a,b,c,d,e) FUNCNAME_I4MAT2PERM0(C_2DT3PI(a,b,c,d,e))
#define i4_bset(a,b) R_INT(FUNCNAME_I4BSET(C_IDT(a,b)))
#define i4_btest(a,b) R_INT(FUNCNAME_I4BTEST(C_IDT(a,b)))
#define index_rank0(a,b,c) R_SHRT(FUNCNAME_INDEXRANK0(C_2DTPI(a,b,c)))
#define index_unrank0(a,b,c,d) FUNCNAME_INDEXUNRANK0(C_3DTPI(a,b,c,d))
#define perm0_lex_next(a,b,c) FUNCNAME_PERM0LEXNEXT(C_DTPIPB(a,b,c))
#define perm0_free(a,b,c,d) FUNCNAME_PERM0FREE(C_2DT2PI(a,c,b,d))
#define ksub_random(a,b,c,d) FUNCNAME_KSUBRANDOM(C_2DT2PI(a,b,c,d))
#define subset_sum_count(a,b,c,d,e) R_USHRT(FUNCNAME_SUBSETSUMCOUNT(C_2DTPI2DT(a,c,b,d,e)))
#define subset_sum_find(a,b,c,d,e,f) FUNCNAME_SUBSETSUMFIND(C_2DTPI2DTPI(a,b,c,d,e,f))
#define subset_sum_table_to_list_length(a,b) R_USHRT(FUNCNAME_SUBSETSUMTABLETOLISTLENGTH(C_DTPI(a,b)))
#define subset_sum_table(a,b,c) FUNCNAME_SUBSETSUMTABLE(C_2DTPI(a,b,c))
#define subset_sum_table_to_list(a,b,c) FUNCNAME_SUBSETSUMTABLETOLIST(C_2DTPI(a,c,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_GENLAGUERRERECUR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HERMITERECUR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_JACOBIRECUR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGUERRERECUR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GENLAGUERREROOT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HERMITEROOT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_JACOBIROOT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGUERREROOT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERM0CYCLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4SQRTCF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_UBVECTOUI4(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_UBVECXOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_UI4TOUBVEC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4SQRT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4BSET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4BTEST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4MAT2PERM0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_INDEXRANK0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_INDEXUNRANK0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERM0LEXNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERM0FREE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_KSUBRANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETSUMCOUNT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETSUMTABLETOLISTLENGTH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERM0INVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETSUMFIND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETSUMTABLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETSUMTABLETOLIST(void *);

#endif
