#ifndef WRAPPER_R8LIB_H_INCLUDED
#define WRAPPER_R8LIB_H_INCLUDED

#define r8_factorial2(a) R_DBL(FUNCNAME_R8FACTORIAL2(C_SUSHRT(a)))
#define r8_fractional(a) R_DBL(FUNCNAME_R8FRACTIONAL(C_SDBL(a)))
#define r8_in_01(a) R_UCHR(FUNCNAME_R8IN01(C_SDBL(a)))
#define r8_is_int(a) R_UCHR(FUNCNAME_R8ISINT(C_SDBL(a)))
#define r8_mant(a,b,c,d) FUNCNAME_R8MANT(C_PITPIPITPI(a,b,c,d))
#define r8_nint(a) R_INT(FUNCNAME_R8NINT(C_SDBL(a)))
#define r8_normal_01(a) R_DBL(FUNCNAME_R8NORMAL01(a))
#define r8_mop(a) R_DBL(FUNCNAME_R8MOP(C_SINT(a)))
#define r8_normal_ab(a,b,c) R_DBL(FUNCNAME_R8NORMALAB(C_2ITPI(a,b,c)))
#define r8_power_fast(a,b,c) R_DBL(FUNCNAME_R8POWERFAST(C_ITIPI(a,b,c)))
#define r8_pythag(a,b) R_DBL(FUNCNAME_R8PYTHAG(C_PDBL2(a,b)))
#define r8_reverse_bytes(a) R_DBL(FUNCNAME_R8REVERSEBYTES(C_SDBL(a)))
#define r8_round_i4(a) R_INT(FUNCNAME_R8ROUNDI4(C_SDBL(a)))
#define r8_round2(a,b) R_DBL(FUNCNAME_R8ROUND2(C_DTIT(a,b)))
#define r8_roundb(a,b,c) R_DBL(FUNCNAME_R8ROUNDB(C_IDTIT(a,b,c)))
#define r8_roundx(a,b) R_DBL(FUNCNAME_R8ROUNDX(C_DTIT(a,b)))
#define r8_sign_opposite_strict(a,b) R_DBL(FUNCNAME_R8SIGNOPPOSITESTRICT(C_PDBL2(a,b)))
#define r8_tiny() R_DBL(FUNCNAME_R8TINY(NULL))
#define r8_to_dhms(a,b,c,d,e) FUNCNAME_R8TODHMS(C_IT4PI(a,b,c,d,e))
#define r8_to_i4(a,b,c,d,e) R_INT(FUNCNAME_R8TOI4(C_3IT2I(a,b,c,d,e)))
#define r8_to_r8_discrete(a,b,c,d) R_DBL(FUNCNAME_R8TOR8DISCRETE(C_3ITI(a,b,c,d)))
#define r8_uniform_ab(a,b,c) R_DBL(FUNCNAME_R8UNIFORMAB(C_2ITPI(a,b,c)))
#define r8_unswap3(a,b,c) FUNCNAME_R8UNSWAP3(C_PPDBL3(a,b,c))
#define r8_walsh_1d(a,b) R_DBL(FUNCNAME_R8WALSH1D(C_ITI(a,b)))
#define r82vec_order_type(a,b) R_SHRT(FUNCNAME_R82VECORDERTYPE(C_DTPIT(a,b)))
#define r82vec_permute(a,b,c,d) FUNCNAME_R82VECPERMUTE(C_DTPIIPIT(a,b,c,d))
#define r82vec_sort_heap_index_a(a,b,c) FUNCNAME_R82VECSORTHEAPINDEXA(C_DTPITI(a,c,b))
#define r8block_zero_new(a,b,c) FUNCNAME_R8BLOCKZERONEW(C_PPUSHRT3(a,b,c))
#define r8col_compare(a,b,c,d,e) R_SHRT(FUNCNAME_R8COLCOMPARE(C_2DTPIT2DT(a,b,c,d,e)))
#define r8col_duplicates(a,b,c,d) FUNCNAME_R8COLDUPLICATES(C_3DTPI(a,b,c,d))
#define r8col_find(a,b,c,d) R_SHRT(FUNCNAME_R8COLFIND(C_2DT2PIT(a,b,c,d)))
#define r8col_first_index(a,b,c,d) FUNCNAME_R8COLFIRSTINDEX(C_2DTPITIT(a,b,c,d))
#define r8col_insert(a,b,c,d,e) R_SHRT(FUNCNAME_R8COLINSERT(C_3DT2PIT(a,b,c,d,e)))
#define r8col_max(a,b,c) FUNCNAME_R8COLMAX(C_2DTPIT(a,b,c))
#define r8col_max_index(a,b,c) FUNCNAME_R8COLMAXINDEX(C_2DTPIT(a,b,c))
#define r8col_max_one(a,b,c) FUNCNAME_R8COLMAXONE(C_2DTPIT(a,b,c))
#define r8col_min(a,b,c) FUNCNAME_R8COLMIN(C_2DTPIT(a,b,c))
#define r8col_min_index(a,b,c) FUNCNAME_R8COLMININDEX(C_2DTPIT(a,b,c))
#define r8col_part_quick_a(a,b,c,d,e) FUNCNAME_R8COLPARTQUICKA(C_2DTPIT2PDT(a,b,c,d,e))
#define r8col_permute(a,b,c,d,e) FUNCNAME_R8COLPERMUTE(C_2DTPIIPIT(a,b,c,d,e))
#define r8col_sort_heap_a(a,b,c) FUNCNAME_R8COLSORTHEAPA(C_2DTPIT(a,b,c))
#define r8col_sort_heap_index_a(a,b,c,d) FUNCNAME_R8COLSORTHEAPINDEXA(C_2DTIPIT(a,b,c,d))
#define r8col_sort_quick_a(a,b,c,d) FUNCNAME_R8COLSORTQUICKA(C_2DTIPIT(a,b,c,d))
#define r8col_sorted_tol_undex(a,b,c,d,e,f,g) FUNCNAME_R8COLSORTEDTOLUNDEX(C_2DTPITDTIT2PI(a,b,c,d,e,f,g))
#define r8col_sorted_tol_unique(a,b,c,d) R_USHRT(FUNCNAME_R8COLSORTEDTOLUNIQUE(C_2DTPITIT(a,b,c,d)))
#define r8col_sorted_tol_unique_count(a,b,c,d) R_USHRT(FUNCNAME_R8COLSORTEDTOLUNIQUECOUNT(C_2DTPITIT(a,b,c,d)))
#define r8col_sorted_undex(a,b,c,d,e,f) FUNCNAME_R8COLSORTEDUNDEX(C_DTIPDTDT2PI(a,b,c,d,e,f))
#define r8col_sorted_unique(a,b,c,d) R_USHRT(FUNCNAME_R8COLSORTEDUNIQUE(C_2DTPIT(a,b,c,d)))
#define r8col_sorted_unique_count(a,b,c,d) R_USHRT(FUNCNAME_R8COLSORTEDUNIQUECOUNT(C_2DTPIT(a,b,c,d)))
#define r8col_sortr_a(a,b,c,d) FUNCNAME_R8COLSORTRA(C_2DTIPIT(a,b,d,c))
#define r8col_sum(a,b,c) FUNCNAME_R8COLSUM(C_2DTPIT(a,b,c))
#define r8col_swap(a,b,c,d,e) FUNCNAME_R8COLSWAP(C_2DTPIT2DT(a,b,c,d,e))
#define r8col_to_r8vec(a,b,c) FUNCNAME_R8COLTOR8VEC(C_2DTPIT(a,b,c))
#define r8col_tol_undex(a,b,c,d,e,f,g) FUNCNAME_R8COLTOLUNDEX(C_2DTPITDTIT2PI(a,b,c,d,e,f,g))
#define r8mat_copy(a,b,c,d) FUNCNAME_R8MATCOPY(C_2DT2PIT(a,b,c,d))
#define r8mat_delete(a,b,c) FUNCNAME_R8MATDELETE(C_PPIT2DT(a,b,c))
#define r8mat_det(a,b) R_DBL(FUNCNAME_R8MATDET(C_DTPIT(a,b)))
#define r8mat_det_2d(a) R_DBL(FUNCNAME_R8MATDET2D(a))
#define r8mat_det_3d(a) R_DBL(FUNCNAME_R8MATDET3D(a))
#define r8mat_det_4d(a) R_DBL(FUNCNAME_R8MATDET4D(a))
#define r8mat_det_5d(a) R_DBL(FUNCNAME_R8MATDET5D(a))
#define r8mat_flip_cols(a,b,c) FUNCNAME_R8MATFLIPCOLS(C_2DTPIT(a,b,c))
#define r8mat_flip_rows(a,b,c) FUNCNAME_R8MATFLIPROWS(C_2DTPIT(a,b,c))
#define r8mat_identity_new(a) FUNCNAME_R8MATIDENTITYNEW(C_SUSHRT(a))
#define r8mat_in_01(a,b,c) R_UCHR(FUNCNAME_R8MATIN01(C_2DTPIT(a,b,c)))
#define r8mat_indicator_new(a,b) FUNCNAME_R8MATINDICATORNEW(C_PUSHRT2(a,b))
#define r8mat_inverse_2d(a) FUNCNAME_R8MATINVERSE2D(a)
#define r8mat_inverse_3d(a) FUNCNAME_R8MATINVERSE3D(a)
#define r8mat_inverse_4d(a) FUNCNAME_R8MATINVERSE4D(a)
// #define r8mat_inverse_5d(a) FUNCNAME_R8MATINVERSE5D(a)
#define r8mat_l_inverse(a,b) FUNCNAME_R8MATLINVERSE(C_DTPIT(a,b))
#define r8mat_l_solve(a,b,c) FUNCNAME_R8MATLSOLVE(C_DT2PIT(a,b,c))
#define r8mat_lt_solve(a,b,c) FUNCNAME_R8MATLTSOLVE(C_DT2PIT(a,b,c))
#define r8mat_max(a,b,c) R_DBL(FUNCNAME_R8MATMAX(C_2DTPIT(a,b,c)))
#define r8mat_min(a,b,c) R_DBL(FUNCNAME_R8MATMIN(C_2DTPIT(a,b,c)))
#define r8mat_mm(a,b,c,d,e,f) FUNCNAME_R8MATMM(C_3DT3PIT(a,b,c,d,e,f))
#define r8mat_mv(a,b,c,d,e) FUNCNAME_R8MATMV(C_2DT3PIT(a,b,c,d,e))
#define r8mat_new(a,b) FUNCNAME_R8MATNEW(C_PUSHRT2(a,b))
#define r8mat_norm_eis(a,b,c) R_DBL(FUNCNAME_R8MATNORMEIS(C_2DTPIT(a,b,c)))
#define r8mat_norm_fro(a,b,c) R_DBL(FUNCNAME_R8MATNORMFRO(C_2DTPIT(a,b,c)))
#define r8mat_mxm(a,b,c,d,e,f) FUNCNAME_R8MATMXM(C_3DT3PIT(a,b,c,d,e,f))
#define r8mat_nullspace(a,b,c,d) FUNCNAME_R8MATNULLSPACE(C_3DTPIT(a,b,d,c))
#define r8mat_nullspace_size(a,b,c) R_USHRT(FUNCNAME_R8MATNULLSPACESIZE(C_2DTPIT(a,b,c)))
#define r8mat_ref(a,b,c) FUNCNAME_R8MATREF(C_2DTPIT(a,b,c))
#define r8mat_rref(a,b,c) FUNCNAME_R8MATRREF(C_2DTPIT(a,b,c))
#define r8mat_symm_eigen(a,b,c) FUNCNAME_R8MATSYMMEIGEN(C_DT2PIT(a,b,c))
#define r8mat_symm_jacobi(a,b) FUNCNAME_R8MATSYMMJACOBI(C_DTPIT(a,b))
#define r8mat_to_r8plu(a,b,c,d) R_USHRT(FUNCNAME_R8MATTOR8PLU(C_PITDTPIPIT(b,a,c,d)))
#define r8mat_trace(a,b) R_DBL(FUNCNAME_R8MATTRACE(C_DTPIT(a,b)))
#define r8mat_transpose(a,b,c) FUNCNAME_R8MATTRANSPOSE(C_2DTPIT(a,b,c))
#define r8mat_transpose_in_place(a,b) FUNCNAME_R8MATTRANSPOSEINPLACE(C_DTPIT(a,b))
#define r8mat_zero(a,b,c) FUNCNAME_R8MATZERO(C_2DTPIT(a,b,c))
#define r8plu_det(a,b,c) R_DBL(FUNCNAME_R8PLUDET(C_DTPITPI(a,c,b)))
#define r8pp_delete(a,b,c) FUNCNAME_R8PPDELETE(C_PPIT2DT(a,b,c))
#define r8pp_new(a,b) FUNCNAME_R8PPNEW(C_PUSHRT2(a,b))
#define r8r8_compare(a,b,c,d) R_SHRT(FUNCNAME_R8R8COMPARE(C_PDBL4(a,b,c,d)))
#define r8r8r8_compare(a,b,c,d,e,f) R_SHRT(FUNCNAME_R8R8R8COMPARE(C_PDBL6(a,b,c,d,e,f)))
#define r8row_max(a,b,c) FUNCNAME_R8ROWMAX(C_2DTPIT(a,b,c))
#define r8row_mean(a,b,c) FUNCNAME_R8ROWMEAN(C_2DTPIT(a,b,c))
#define r8row_min(a,b,c) FUNCNAME_R8ROWMIN(C_2DTPIT(a,b,c))
#define r8row_sum(a,b,c) FUNCNAME_R8ROWSUM(C_2DTPIT(a,b,c))
#define r8row_swap(a,b,c,d,e) FUNCNAME_R8ROWSWAP(C_2DTPIT2DT(a,b,c,d,e))
#define r8row_to_r8vec(a,b,c) FUNCNAME_R8ROWTOR8VEC(C_2DTPIT(a,b,c))
#define r8vec_01_to_ab(a,b,c,d) FUNCNAME_R8VEC01TOAB(C_DT2ITPIT(a,c,d,b))
#define r8vec_ab_to_01(a,b) FUNCNAME_R8VECABTO01(C_DTPIT(a,b))
#define r8vec_ab_to_cd(a,b,c,d) FUNCNAME_R8VECABTOCD(C_DT2ITPIT(a,c,d,b))
#define r8vec_amax(a,b) R_DBL(FUNCNAME_R8VECAMAX(C_DTPIT(a,b)))
#define r8vec_amax_index(a,b) R_USHRT(FUNCNAME_R8VECAMAXINDEX(C_DTPIT(a,b)))
#define r8vec_amin(a,b) R_DBL(FUNCNAME_R8VECAMIN(C_DTPIT(a,b)))
#define r8vec_amin_index(a,b) R_USHRT(FUNCNAME_R8VECAMININDEX(C_DTPIT(a,b)))
#define r8vec_any_normal(a,b) FUNCNAME_R8VECANYNORMAL(C_DTPIT(a,b))
#define r8vec_circular_variance(a,b) R_DBL(FUNCNAME_R8VECCIRCULARVARIANCE(C_DTPIT(a,b)))
#define r8vec_compare(a,b,c) R_SHRT(FUNCNAME_R8VECCOMPARE(C_DT2PIT(a,b,c)))
#define r8vec_convolve_circ(a,b,c) FUNCNAME_R8VECCONVOLVECIRC(C_DT2PIT(a,b,c))
#define r8vec_correlation(a,b,c) R_DBL(FUNCNAME_R8VECCORRELATION(C_DT2PIT(a,b,c)))
#define r8vec_covar(a,b,c) R_DBL(FUNCNAME_R8VECCOVAR(C_DT2PIT(a,b,c)))
#define r8vec_cross_product_2d(a,b) R_DBL(FUNCNAME_R8VECCROSSPRODUCT2D(C_PPDBL2(a,b)))
#define r8vec_cross_product_affine_2d(a,b,c) R_DBL(FUNCNAME_R8VECCROSSPRODUCTAFFINE2D(C_PPDBL3(a,b,c)))
#define r8vec_cross_product_3d(a,b) FUNCNAME_R8VECCROSSPRODUCT3D(C_PPDBL2(a,b))
#define r8vec_cross_product_affine_3d(a,b,c) R_DBL(FUNCNAME_R8VECCROSSPRODUCTAFFINE3D(C_PPDBL3(a,b,c)))
#define r8vec_dif(a,b) FUNCNAME_R8VECDIF(C_DTIT(a,b))
#define r8vec_diff_norm(a,b,c) R_DBL(FUNCNAME_R8VECDIFFNORM(C_DT2PIT(a,b,c)))
#define r8vec_diff_norm_l1(a,b,c) R_DBL(FUNCNAME_R8VECDIFFNORML1(C_DT2PIT(a,b,c)))
#define r8vec_diff_norm_l2(a,b,c) R_DBL(FUNCNAME_R8VECDIFFNORML2(C_DT2PIT(a,b,c)))
#define r8vec_diff_norm_li(a,b,c) R_DBL(FUNCNAME_R8VECDIFFNORMLI(C_DT2PIT(a,b,c)))
#define r8vec_diff_norm_squared(a,b,c) R_DBL(FUNCNAME_R8VECDIFFNORMSQUARED(C_DT2PIT(a,b,c)))
#define r8vec_distance(a,b,c) R_DBL(FUNCNAME_R8VECDISTANCE(C_DT2PIT(a,b,c)))
#define r8vec_distinct(a,b) R_UCHR(FUNCNAME_R8VECDISTINCT(C_DTPIT(a,b)))
#define r8vec_divide(a,b,c) FUNCNAME_R8VECDIVIDE(C_DTPITIT(a,b,c))
#define r8vec_dot_product(a,b,c) R_DBL(FUNCNAME_R8VECDOTPRODUCT(C_DT2PIT(a,b,c)))
#define r8vec_dot_product_affine(a,b,c,d) R_DBL(FUNCNAME_R8VECDOTPRODUCTAFFINE(C_DT2PIT(a,b,c,d)))
#define r8vec_eq(a,b,c) R_UCHR(FUNCNAME_R8VECEQ(C_DT2PIT(a,b,c)))
#define r8vec_even_select(a,b,c,d) R_DBL(FUNCNAME_R8VECEVENSELECT(C_DT2ITDT(a,b,c,d)))
#define r8vec_even2(a,b,c,d,e,f) FUNCNAME_R8VECEVEN2(C_DTPIDTPITPIPIT(a,b,c,d,e,f))
#define r8vec_even3(a,b,c,d) FUNCNAME_R8VECEVEN3(C_2DT2PIT(a,b,c,d))
#define r8vec_expand_linear(a,b,c) FUNCNAME_R8VECEXPANDLINEAR(C_2DTPIT(a,c,b))
#define r8vec_first_index(a,b,c) FUNCNAME_R8VECFIRSTINDEX(C_DTPITIT(a,b,c))
#define r8vec_frac(a,b,c) R_DBL(FUNCNAME_R8VECFRAC(C_2DTPIT(a,c,b)))
#define r8vec_fraction(a,b) FUNCNAME_R8VECFRACTION(C_DTPIT(a,b))
#define r8vec_gt(a,b,c) R_UCHR(FUNCNAME_R8VECGT(C_DT2PIT(a,b,c)))
#define r8vec_heap_a(a,b) FUNCNAME_R8VECHEAPA(C_DTPIT(a,b))
#define r8vec_heap_d(a,b) FUNCNAME_R8VECHEAPD(C_DTPIT(a,b))
#define i4vec_zero(a,b) FUNCNAME_I4VECZERO(C_DTPI(a,b))
#define r8vec_histogram(a,b,c,d,e) FUNCNAME_R8VECHISTOGRAM(C_DT2ITDTPIT(a,c,d,e,b))
#define r8vec_house_column(a,b,c) FUNCNAME_R8VECHOUSECOLUMN(C_2DTPIT(a,c,b))
#define r8vec_i4vec_dot_product(a,b,c) R_DBL(FUNCNAME_R8VECI4VECDOTPRODUCT(C_DTPITPI(a,b,c)))
#define r8vec_in_01(a,b) R_UCHR(FUNCNAME_R8VECIN01(C_DTPIT(a,b)))
#define r8vec_index_delete_all(a,b,c,d,e,f,g) FUNCNAME_R8VECINDEXDELETEALL(C_DTPITPIITPDTPITPI(a,b,c,d,e,f,g))
#define r8vec_index_delete_dupes(a,b,c,d,e,f) FUNCNAME_R8VECINDEXDELETEDUPES(C_DTPITPIPDTPITPI(a,b,c,d,e,f))
#define r8vec_index_delete_one(a,b,c,d,e,f,g) FUNCNAME_R8VECINDEXDELETEONE(C_DTPITPIITPDTPITPI(a,b,c,d,e,f,g))
#define r8vec_index_insert(a,b,c,d) FUNCNAME_R8VECINDEXINSERT(C_PDTPITPIIT(a,b,c,d))
#define r8vec_index_insert_unique(a,b,c,d) FUNCNAME_R8VECINDEXINSERTUNIQUE(C_PDTPITPIIT(a,b,c,d))
#define r8vec_index_order(a,b,c) FUNCNAME_R8VECINDEXORDER(C_DTPITPI(a,b,c))
#define r8vec_index_search(a,b,c,d,e,f,g) FUNCNAME_R8VECINDEXSEARCH(C_DTPITPIIT3PDT(a,b,c,d,e,f,g))
#define r8vec_index_sort_unique(a,b,c,d,e) FUNCNAME_R8VECINDEXSORTUNIQUE(C_DT2PITPDTPI(a,b,d,c,e))
#define r8vec_index_sorted_range(a,b,c,d,e,f,g) FUNCNAME_R8VECINDEXSORTEDRANGE(C_DTPITPI2IT2PI(a,b,c,d,e,f,g))
#define r8vec_indexed_heap_d(a,b,c) FUNCNAME_R8VECINDEXEDHEAPD(C_DTPITPI(a,b,c))
#define r8vec_indexed_heap_d_extract(a,b,c) R_USHRT(FUNCNAME_R8VECINDEXEDHEAPDEXTRACT(C_PDTPITPI(a,b,c)))
#define r8vec_indexed_heap_d_insert(a,b,c,d) FUNCNAME_R8VECINDEXEDHEAPDINSERT(C_DTPITPDTPI(d,b,a,c))
#define r8vec_indexed_heap_d_max(a,b,c) R_USHRT(FUNCNAME_R8VECINDEXEDHEAPDMAX(C_DTPITPI(a,b,c)))
#define r8vec_indicator0(a,b) FUNCNAME_R8VECINDICATOR0(C_DTPIT(a,b))
#define r8vec_indicator0_new(a) FUNCNAME_R8VECINDICATOR0NEW(C_SUSHRT(a))
#define r8vec_indicator1(a,b) FUNCNAME_R8VECINDICATOR1(C_DTPIT(a,b))
#define r8vec_indicator1_new(a) FUNCNAME_R8VECINDICATOR1NEW(C_SUSHRT(a))
#define r8vec_insert(a,b,c,d) FUNCNAME_R8VECINSERT(C_2DTPITIT(a,c,b,d))
#define r8vec_is_int(a,b) R_UCHR(FUNCNAME_R8VECISINT(C_DTPIT(a,b)))
#define r8vec_is_nonnegative(a,b) R_UCHR(FUNCNAME_R8VECISNONNEGATIVE(C_DTPIT(a,b)))
#define r8vec_is_zero(a,b) R_UCHR(FUNCNAME_R8VECISZERO(C_DTPIT(a,b)))
#define r8vec_lt(a,b,c) R_UCHR(FUNCNAME_R8VECLT(C_DT2PIT(a,b,c)))
#define r8vec_max_index(a,b) R_SHRT(FUNCNAME_R8VECMAXINDEX(C_DTPIT(a,b)))
#define r8vec_min_index(a,b) R_SHRT(FUNCNAME_R8VECMININDEX(C_DTPIT(a,b)))
#define r8vec_min_pos(a,b) R_DBL(FUNCNAME_R8VECMINPOS(C_DTPIT(a,b)))
#define r8vec_mirror_next(a,b) R_INT(FUNCNAME_R8VECMIRRORNEXT(C_DTPIT(a,b)))
#define r8vec_negative_strict(a,b) R_UCHR(FUNCNAME_R8VECNEGATIVESTRICT(C_DTPIT(a,b)))
#define r8vec_nint(a,b) FUNCNAME_R8VECNINT(C_DTPIT(a,b))
#define r8vec_norm(a,b) R_DBL(FUNCNAME_R8VECNORM(C_DTPIT(a,b)))
#define r8vec_norm_affine(a,b,c) R_DBL(FUNCNAME_R8VECNORMAFFINE(C_DT2PIT(a,b,c)))
#define r8vec_norm_l1(a,b) R_DBL(FUNCNAME_R8VECNORML1(C_DTPIT(a,b)))
#define r8vec_norm_l2(a,b) R_DBL(FUNCNAME_R8VECNORML2(C_DTPIT(a,b)))
#define r8vec_norm_li(a,b) R_DBL(FUNCNAME_R8VECNORMLI(C_DTPIT(a,b)))
#define r8vec_norm_lp(a,b,c) R_DBL(FUNCNAME_R8VECNORMLP(C_DTPITIT(a,b,c)))
#define r8vec_normalize(a,b) FUNCNAME_R8VECNORMALIZE(C_DTPIT(a,b))
#define r8vec_normalize_l1(a,b) FUNCNAME_R8VECNORMALIZEL1(C_DTPIT(a,b))
#define r8vec_normsq(a,b) R_DBL(FUNCNAME_R8VECNORMSQ(C_DTPIT(a,b)))
#define r8vec_normsq_affine(a,b,c) R_DBL(FUNCNAME_R8VECNORMSQAFFINE(C_DT2PIT(a,b,c)))
#define r8vec_order_type(a,b) R_SHRT(FUNCNAME_R8VECORDERTYPE(C_DTPIT(a,b)))
#define r8vec_part_quick_a(a,b,c,d) FUNCNAME_R8VECPARTQUICKA(C_DTPIT2PDT(a,b,c,d))
#define r8vec_permute(a,b,c,d) FUNCNAME_R8VECPERMUTE(C_DTPIIPIT(a,b,c,d))
#define r8vec_permute_cyclic(a,b,c) FUNCNAME_R8VECPERMUTECYCLIC(C_2DTPIT(a,b,c))
#define r8vec_permute_uniform(a,b,c,d) FUNCNAME_R8VECPERMUTEUNIFORM(C_DTPITPI(a,b,c,d))
#define r8vec_polarize(a,b,c,d,e) FUNCNAME_R8VECPOLARIZE(C_DT4PIT(a,b,c,d,e))
#define r8vec_positive_strict(a,b) R_UCHR(FUNCNAME_R8VECPOSITIVESTRICT(C_DTPIT(a,b)))
#define r8vec_range(a,b,c,d,e,f,g) FUNCNAME_R8VECRANGE(C_DTPIT2IT3PIT(a,b,c,d,e,f,g))
#define r8vec_range_2(a,b,c,d) FUNCNAME_R8VECRANGE2(C_DT3PIT(a,b,c,d))
#define r8vec_reverse(a,b) FUNCNAME_R8VECREVERSE(C_DTPIT(a,b))
#define r8vec_rotate(a,b,c) FUNCNAME_R8VECROTATE(C_2DTPIT(a,c,b))
#define r8vec_scalar_triple_product(a,b,c) R_DBL(FUNCNAME_R8VECSCALARTRIPLEPRODUCT(C_PPDBL3(a,b,c)))
#define r8vec_search_binary_a(a,b,c) R_INT(FUNCNAME_R8VECSEARCHBINARYA(C_DTPITIT(a,b,c)))
#define r8vec_shift(a,b,c) FUNCNAME_R8VECSHIFT(C_2DTPIT(a,b,c))
#define r8vec_shift_circular(a,b,c) FUNCNAME_R8VECSHIFTCIRCULAR(C_2DTPIT(a,b,c))
#define r8vec_sort_bubble_a(a,b) FUNCNAME_R8VECSORTBUBBLEA(C_DTPIT(a,b))
#define r8vec_sort_bubble_d(a,b) FUNCNAME_R8VECSORTBUBBLED(C_DTPIT(a,b))
#define r8vec_sort_heap_a(a,b) FUNCNAME_R8VECSORTHEAPA(C_DTPIT(a,b))
#define r8vec_sort_heap_d(a,b) FUNCNAME_R8VECSORTHEAPD(C_DTPIT(a,b))
#define r8vec_sort_heap_index_a(a,b,c) FUNCNAME_R8VECSORTHEAPINDEXA(C_DTPITPI(a,b,c))
#define r8vec_sort_heap_index_a_new(a,b) FUNCNAME_R8VECSORTHEAPINDEXANEW(C_DTPIT(a,b))
#define r8vec_sort_heap_index_d(a,b,c) FUNCNAME_R8VECSORTHEAPINDEXD(C_DTPITPI(a,b,c))
#define r8vec_sort_heap_index_d_new(a,b) FUNCNAME_R8VECSORTHEAPINDEXDNEW(C_DTPIT(a,b))
#define r8vec_sort_heap_mask_a(a,b,c,d) FUNCNAME_R8VECSORTHEAPMASKA(C_2DTPIPIT(a,c,d,b))
#define r8vec_sort_insert_a(a,b) FUNCNAME_R8VECSORTINSERTA(C_DTPIT(a,b))
#define r8vec_sort_insert_index_a(a,b) FUNCNAME_R8VECSORTINSERTINDEXA(C_DTPIT(a,b))
#define r8vec_sort_quick_a(a,b) FUNCNAME_R8VECSORTQUICKA(C_DTPIT(a,b))
#define r8vec_sort_shell_a(a,b) FUNCNAME_R8VECSORTSHELLA(C_DTPIT(a,b))
#define r8vec_sorted_merge_a(a,b,c,d,e) FUNCNAME_R8VECSORTEDMERGEA(C_PIT2DTPDTPIT(b,c,a,e,d))
#define r8vec_sorted_nearest(a,b,c) R_SHRT(FUNCNAME_R8VECSORTEDNEAREST(C_DTPITIT(a,b,c)))
#define r8vec_sorted_range(a,b,c,d,e,f) FUNCNAME_R8VECSORTEDRANGE(C_DTPIT2PI2IT(a,b,e,f,c,d))
#define r8vec_sorted_split(a,b,c,d,e) FUNCNAME_R8VECSORTEDSPLIT(C_DTPITIT2PI(a,b,c,d,e))
#define r8vec_sorted_undex(a,b,c,d,e,f) FUNCNAME_R8VECSORTEDUNDEX(C_DTPITDTIT2PI(a,b,c,d,e,f))
#define r8vec_sorted_unique(a,b,c,d) FUNCNAME_R8VECSORTEDUNIQUE(C_DTPITITPDT(a,b,c,d))
#define r8vec_sorted_unique_count(a,b,c) R_USHRT(FUNCNAME_R8VECSORTEDUNIQUECOUNT(C_DTPITIT(a,b,c)))
#define r8vec_sorted_unique_hist(a,b,c,d,e,f,g) FUNCNAME_R8VECSORTEDUNIQUEHIST(C_DTPITITDTPDTPITPI(a,b,c,d,e,f,g))
#define r8vec_split(a,b,c) R_USHRT(FUNCNAME_R8VECSPLIT(C_DTPITIT(a,b,c)))
#define r8vec_stutter(a,b,c,d) FUNCNAME_R8VECSTUTTER(C_2DT2PIT(a,c,b,d))
#define r8vec_stutter_new(a,b,c) FUNCNAME_R8VECSTUTTERNEW(C_2DTPIT(a,c,b))
#define r8vec_swap(a,b,c) FUNCNAME_R8VECSWAP(C_DT2PIT(a,b,c))
#define r8vec_undex(a,b,c,d,e,f) FUNCNAME_R8VECUNDEX(C_DTPITDTIT2PI(a,b,c,d,e,f))
#define r8vec_uniform_ab(a,b,c,d,e) FUNCNAME_R8VECUNIFORMAB(C_DT2ITPIPIT(a,b,c,d,e))
#define r8vec_unique_count(a,b,c) R_USHRT(FUNCNAME_R8VECUNIQUECOUNT(C_DTPITIT(a,b,c)))
#define r8vec_unique_index(a,b,c) FUNCNAME_R8VECUNIQUEINDEX(C_DTPITIT(a,b,c))
#define r8vec_vector_triple_product(a,b,c) FUNCNAME_R8VECVECTORTRIPLEPRODUCT(C_PPDBL3(a,b,c))
#define r8vec_zero(a,b) FUNCNAME_R8VECZERO(C_DTPIT(a,b))
#define r8vec2_compare(a,b,c,d,e) R_SHRT(FUNCNAME_R8VEC2COMPARE(C_3DT2PIT(a,d,e,b,c)))
#define r8vec2_sort_a(a,b,c) FUNCNAME_R8VEC2SORTA(C_DT2PIT(a,b,c))
#define r8vec2_sort_d(a,b,c) FUNCNAME_R8VEC2SORTD(C_DT2PIT(a,b,c))
#define r8vec2_sort_heap_index_a(a,b,c) FUNCNAME_R8VEC2SORTHEAPINDEXA(C_DT2PIT(a,b,c))
#define r8vec2_sorted_unique(a,b,c,d) FUNCNAME_R8VEC2SORTEDUNIQUE(C_DT2PITPDT(a,b,c,d))
#define r8vec2_sorted_unique_index(a,b,c,d,e) FUNCNAME_R8VEC2SORTEDUNIQUEINDEX(C_DT2PITPDTPI(a,b,c,d,e))
#define r8vec2_sum_max_index(a,b,c) R_SHRT(FUNCNAME_R8VEC2SUMMAXINDEX(C_DT2PIT(a,b,c)))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8R8COMPARE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8R8R8COMPARE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8WALSH1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8TOR8DISCRETE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8TOI4(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8NINT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8FACTORIAL2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MANT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8NORMALAB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POWERFAST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8TODHMS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8UNIFORMAB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R82VECORDERTYPE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R82VECPERMUTE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLCOMPARE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLFIND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLINSERT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLMAXONE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLPARTQUICKA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLPERMUTE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTHEAPA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTQUICKA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTEDTOLUNDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTEDTOLUNIQUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTEDTOLUNIQUECOUNT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTEDUNDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTEDUNIQUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTEDUNIQUECOUNT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTRA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSWAP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLTOLUNDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATCOPY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATDELETE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATDET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATDET2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATDET3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATDET4D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATDET5D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATFLIPCOLS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATFLIPROWS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATIN01(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMAX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATNORMEIS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATNORMFRO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMXM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATNULLSPACESIZE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATREF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATRREF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATSYMMJACOBI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATTOR8PLU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATTRACE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATTRANSPOSEINPLACE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATZERO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8PLUDET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8PPDELETE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROWSWAP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VEC01TOAB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECABTO01(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECAMAX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECAMAXINDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECAMIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECAMININDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCIRCULARVARIANCE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCOMPARE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCORRELATION(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCOVAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCROSSPRODUCT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCROSSPRODUCTAFFINE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDIFFNORM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDIFFNORML1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDIFFNORML2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDIFFNORMLI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDIFFNORMSQUARED(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDISTANCE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDISTINCT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDIVIDE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDOTPRODUCT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDOTPRODUCTAFFINE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECEQ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECEVEN2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECEVEN3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECFRAC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECGT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECHEAPA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECHEAPD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECZERO(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECI4VECDOTPRODUCT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECIN01(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXDELETEALL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXDELETEDUPES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXDELETEONE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXINSERT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXINSERTUNIQUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXSEARCH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXSORTUNIQUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXSORTEDRANGE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXEDHEAPD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXEDHEAPDEXTRACT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDEXEDHEAPDINSERT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECINDEXEDHEAPDMAX(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECINDICATOR0(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECINDICATOR1(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECINSERT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECISINT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECISNONNEGATIVE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECISZERO(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECLT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECMAXINDEX(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECMININDEX(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECMINPOS(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECMIRRORNEXT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNEGATIVESTRICT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORM(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORMAFFINE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORML1(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORML2(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORMLI(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORMLP(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORMALIZE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORMALIZEL1(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORMSQ(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECNORMSQAFFINE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECORDERTYPE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECPARTQUICKA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECPERMUTE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECPERMUTECYCLIC(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECPERMUTEUNIFORM(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECPOLARIZE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECPOSITIVESTRICT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECRANGE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECRANGE2(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECREVERSE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECROTATE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSCALARTRIPLEPRODUCT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSEARCHBINARYA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSHIFT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSHIFTCIRCULAR(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTBUBBLEA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTBUBBLED(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTHEAPA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTHEAPD(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTHEAPINDEXA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTHEAPINDEXD(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTINSERTA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTQUICKA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTSHELLA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTEDNEAREST(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTEDRANGE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTEDSPLIT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTEDUNDEX(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTEDUNIQUECOUNT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSORTEDUNIQUEHIST(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSPLIT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSTUTTER(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECSWAP(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECUNDEX(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECUNIFORMAB(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECUNIQUECOUNT(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VECZERO(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VEC2COMPARE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VEC2SORTA(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VEC2SORTD(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VEC2SORTEDUNIQUE(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VEC2SORTEDUNIQUEINDEX(void *);
__MATHSUITE __JBURKARDT void *  FUNCNAME_R8VEC2SUMMAXINDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R82CHEBY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R82VECSORTHEAPINDEXA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8BLOCKZERONEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLDUPLICATES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLFIRSTINDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLMAX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLMAXINDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLMIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLMININDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSORTHEAPINDEXA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLSUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COLTOR8VEC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATIDENTITYNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATINDICATORNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATINVERSE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATINVERSE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATINVERSE4D(void *);
// __MATHSUITE __JBURKARDT void * FUNCNAME_R8MATINVERSE5D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATLINVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATLSOLVE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATLTSOLVE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATNULLSPACE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATSYMMEIGEN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATTRANSPOSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8PPNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROWMAX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROWMEAN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROWMIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROWSUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROWTOR8VEC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECABTOCD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECANYNORMAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCONVOLVECIRC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCROSSPRODUCT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECCROSSPRODUCTAFFINE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECDIF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECEXPANDLINEAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECFIRSTINDEX(void *);	
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECFRACTION(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECHISTOGRAM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECHOUSECOLUMN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDICATOR0NEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECINDICATOR1NEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECNINT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSORTHEAPINDEXANEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSORTHEAPINDEXDNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSORTHEAPMASKA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSORTINSERTINDEXA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSORTEDMERGEA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSORTEDUNIQUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSTUTTERNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECUNIQUEINDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECVECTORTRIPLEPRODUCT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VEC2SORTHEAPINDEXA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8TINY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8FRACTIONAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8IN01(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ISINT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8NORMAL01(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MOP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8REVERSEBYTES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROUNDI4(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8PYTHAG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROUND2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROUNDX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SIGNOPPOSITESTRICT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ROUNDB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECEVENSELECT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8UNSWAP3(void *);

#endif
