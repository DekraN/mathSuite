#ifndef WRAPPER_JBURKARDT_H_INCLUDED
#define WRAPPER_JBURKARDT_H_INCLUDED

#define r8poly_value(a,b,c) R_DBL(FUNCNAME_R8POLYVALUE(C_DTPITIT(a,b,c)))
#define p_polynomial_prime2(a,b) R_DBL(FUNCNAME_PPOLYNOMIALPRIME2(C_DTIT(a,b)))
#define p_polynomial_value(a,b) R_DBL(FUNCNAME_PPOLYNOMIALVALUE(C_DTIT(a,b)))
#define pm_polynomial_value(a,b,c) R_DBL(FUNCNAME_PMPOLYNOMIALVALUE(C_2DTIT(a,b,c)))
#define pmn_polynomial_value(a,b,c) R_DBL(FUNCNAME_PMNPOLYNOMIALVALUE(C_2DTIT(a,b,c)))
#define pmns_polynomial_value(a,b,c) R_DBL(FUNCNAME_PMNSPOLYNOMIALVALUE(C_2DTIT(a,b,c)))
#define pn_polynomial_value(a,b) R_DBL(FUNCNAME_PNPOLYNOMIALVALUE(C_DTIT(a,b)))
#define pn_polynomial_coefficients(a,b) R_UCHR(FUNCNAME_PNPOLYNOMIALCOEFFICIENTS(C_DTPIT(a,b)))
#define r8vec_normal_01_new(a,b) FUNCNAME_R8VECNORMAL01NEW(C_DTPI(a,b))
#define i4vec_uniform_ab_new(a,b,c,d) FUNCNAME_I4VECUNIFORMABNEW(C_3DTPI(a,b,c,d))
#define gauss(a,b,c,d,e) FUNCNAME_GAUSS(C_DT4PI(a,b,c,d,e))
#define r_jacobi(a,b,c,d,e) FUNCNAME_RJACOBI(C_DT2IT2PIT(a,b,c,d,e))
#define diaphony_compute(a,b,c) R_DBL(FUNCNAME_DIAPHONYCOMPUTE(C_2DTPIT(a,b,c)))
#define r8_modp(a,b) R_DBL(FUNCNAME_R8MODP(C_PDBL2(a,b)))
#define inverse_error(a,b,c) R_DBL(FUNCNAME_INVERSEERROR(C_DT2PIT(a,b,c)))
#define r8po_sl2(a,b,c) FUNCNAME_R8POSL2(C_DT2PIT(a,b,c))
#define chebyshev1_exactness(a,b,c,d) FUNCNAME_CHEBYSHEV1EXACTNESS(C_2DT2PIT(a,d,b,c))
#define chebyshev1_integral(a) R_DBL(FUNCNAME_CHEBYSHEV1INTEGRAL(C_SUSHRT(a)))
#define chebyshev2_exactness(a,b,c,d) FUNCNAME_CHEBYSHEV2EXACTNESS(C_2DT2PIT(a,d,b,c))
#define chebyshev2_integral(a) R_DBL(FUNCNAME_CHEBYSHEV2INTEGRAL(C_SUSHRT(a)))
#define r8_cosd(a) R_DBL(FUNCNAME_R8COSD(C_SDBL(a)))
#define cotd(a) R_DBL(FUNCNAME_COTD(C_SDBL(a)))
#define r8_cscd(a) R_DBL(FUNCNAME_R8CSCD(C_SDBL(a)))
#define r8_secd(a) R_DBL(FUNCNAME_R8SECD(C_SDBL(a)))
#define r8_sind(a) R_DBL(FUNCNAME_R8SIND(C_SDBL(a)))
#define r8_tand(a) R_DBL(FUNCNAME_R8TAND(C_SDBL(a)))
#define c4_cos(a,b) R_COMPLEX(FUNCNAME_C4COS(C_COMPLEX(a,b)))
#define c4_sin(a,b) R_COMPLEX(FUNCNAME_C4SIN(C_COMPLEX(a,b)))
#define r8_rand(a) R_DBL(FUNCNAME_R8RAND(C_SDBL(a)))
#define r8_randgs(a,b) R_DBL(FUNCNAME_R8RANDGS(C_PDBL2(a,b)))
#define r8_admp(a,b,c) FUNCNAME_R8ADMP(C_IT2PIT(a,b,c))
#define r8_ai(a) R_DBL(FUNCNAME_R8AI(C_SDBL(a)))
#define r8_aid(a) R_DBL(FUNCNAME_R8AID(C_SDBL(a)))
#define r8_aide(a) R_DBL(FUNCNAME_R8AIDE(C_SDBL(a)))
#define r8_aie(a) R_DBL(FUNCNAME_R8AIE(C_SDBL(a)))
#define r8_aimp(a,b,c) FUNCNAME_R8AIMP(C_IT2PIT(a,b,c))
#define r8_aint(a) R_DBL(FUNCNAME_R8AINT(C_SDBL(a)))
#define r8_asinh(a) R_DBL(FUNCNAME_R8ASINH(C_SDBL(a)))
#define r8_b1mp(a,b,c) FUNCNAME_R8B1MP(C_IT2PIT(a,b,c))
#define r8_bi(a) R_DBL(FUNCNAME_R8BI(C_SDBL(a)))
#define r8_bid(a) R_DBL(FUNCNAME_R8BID(C_SDBL(a)))
#define r8_bide(a) R_DBL(FUNCNAME_R8BIDE(C_SDBL(a)))
#define r8_bie(a) R_DBL(FUNCNAME_R8BIE(C_SDBL(a)))
#define r8_binom(a,b) R_DBL(FUNCNAME_R8BINOM(C_PUSHRT2(a,b)))
#define r8_chi(a) R_DBL(FUNCNAME_R8CHI(C_SDBL(a)))
#define r8_chu(a,b,c) R_DBL(FUNCNAME_R8CHU(C_PDBL3(a,b,c)))
#define r8_chu_scaled(a,b,c) R_DBL(FUNCNAME_R8CHUSCALED(C_PDBL3(a,b,c)))
#define r8_ci(a) R_DBL(FUNCNAME_R8CI(C_SDBL(a)))
#define r8_cin(a) R_DBL(FUNCNAME_R8CIN(C_SDBL(a)))
#define r8_cinh(a) R_DBL(FUNCNAME_R8CINH(C_SDBL(a)))
#define r8_cos_deg(a) R_DBL(FUNCNAME_R8COSDEG(C_SDBL(a)))
#define r8_dawson(a) R_DBL(FUNCNAME_R8DAWSON(C_SDBL(a)))
#define r8_e1(a) R_DBL(FUNCNAME_R8E1(C_SDBL(a)))
#define r8_ei(a) R_DBL(FUNCNAME_R8EI(C_SDBL(a)))
#define r8_erf(a) R_DBL(FUNCNAME_R8ERF(C_SDBL(a)))
#define r8_erfc(a) R_DBL(FUNCNAME_R8ERFC(C_SDBL(a)))
#define r8_exprel(a) R_DBL(FUNCNAME_R8EXPREL(C_SDBL(a)))
#define r8_fac(a) R_DBL(FUNCNAME_R8FAC(C_SUSHRT(a)))
#define r8_gami(a,b) R_DBL(FUNCNAME_R8GAMI(C_PDBL2(a,b)))
#define r8_gamic(a,b) R_DBL(FUNCNAME_R8GAMIC(C_PDBL2(a,b)))
#define r8_gamit(a,b) R_DBL(FUNCNAME_R8GAMIT(C_PDBL2(a,b)))
#define r8_gamr(a) R_DBL(FUNCNAME_R8GAMR(C_SDBL(a)))
#define r8_gmic(a,b,c) R_DBL(FUNCNAME_R8GMIC(C_PDBL3(a,b,c)))
#define r8_gmit(a,b,c,d,e) R_DBL(FUNCNAME_R8GMIT(C_PDBL5(a,b,c,d,e)))
#define r8_int(a) R_DBL(FUNCNAME_R8INT(C_SDBL(a)))
#define r8_lbeta(a,b) R_DBL(FUNCNAME_R8LBETA(C_PDBL2(a,b)))
#define r8_lgams(a,b,c) FUNCNAME_R8LGAMS(C_IT2PIT(a,b,c))
#define r8_lgic(a,b,c) R_DBL(FUNCNAME_R8LGIC(C_PDBL3(a,b,c)))
#define r8_lgit(a,b,c) R_DBL(FUNCNAME_R8LGIT(C_PDBL3(a,b,c)))
#define r8_li(a) R_DBL(FUNCNAME_R8LI(C_SDBL(a)))
#define r8_lngam(a) R_DBL(FUNCNAME_R8LNGAM(C_SDBL(a)))
#define r8_lnrel(a) R_DBL(FUNCNAME_R8LNREL(C_SDBL(a)))
#define r8_mod(a,b) R_DBL(FUNCNAME_R8MOD(C_PDBL2(a,b)))
#define r8_pak(a,b) R_DBL(FUNCNAME_R8PAK(C_DTIT(b,a)))
#define r8_poch(a,b) R_DBL(FUNCNAME_R8POCH(C_PDBL2(a,b)))
#define comp_next(a,b,c,d,e,f) FUNCNAME_COMPNEXT(C_2DT4PI(a,b,c,d,e,f))
#define r8_poch1(a,b) R_DBL(FUNCNAME_R8POCH1(C_PDBL2(a,b)))
#define r8_ren() R_DBL(FUNCNAME_R8REN())
#define r8_shi(a) R_DBL(FUNCNAME_R8SHI(C_SDBL(a)))
#define r8_si(a) R_DBL(FUNCNAME_R8SI(C_SDBL(a)))
#define r8_sifg(a,b,c) FUNCNAME_R8SIFG(C_IT2PIT(a,b,c))
#define r8_sin_deg(a,b,c) FUNCNAME_R8SINDEG(C_IT2PIT(a,b,c))
#define r8_spence(a) R_DBL(FUNCNAME_R8SPENCE(C_SDBL(a)))
#define r8_sqrt(a) R_DBL(FUNCNAME_R8SQRT(C_SDBL(a)))
#define r8_tan(a) R_DBL(FUNCNAME_R8TAN(C_SDBL(a)))
#define r8_tanh(a) R_DBL(FUNCNAME_R8TANH(C_SDBL(a)))
#define r8_upak(a,b,c) FUNCNAME_R8UPAK(C_ITPITPDT(a,b,c))
#define r8vec_uniform_ab_new(a,b,c,d) FUNCNAME_R8VECUNIFORMABNEW(C_DT2ITPI(a,b,c,d))
#define r8vec_sorted_nearest0(a,b,c) R_SHRT(FUNCNAME_R8VECSORTEDNEAREST0(C_DTPITIT(a,b,c)))
#define jacobi1(a,b,c,d) FUNCNAME_JACOBI1(C_DT3PIT(a,b,c,d))
#define lagrange_basis_function_1d(a,b,c,d) R_DBL(FUNCNAME_LAGRANGEBASISFUNCTION1D(C_2DTPITIT(a,c,b,d)))
#define lagrange_interp_2d(a,b,c,d,e,f,g,h) FUNCNAME_LAGRANGEINTERP2D(C_2DT3PITDT2PIT(a,b,c,d,e,f,g,h))
#define laguerre_monomial_quadrature(a,b,c,d,e) R_DBL(FUNCNAME_LAGUERREMONOMIALQUADRATURE(C_3DT2PIT(a,b,c,d,e)))
#define lu_error(a,b,c,d) R_DBL(FUNCNAME_LUERROR(C_DT3PIT(a,b,c,d)))
#define r8vec_linspace2(a,b,c,d) FUNCNAME_R8VECLINSPACE2(C_DT2ITPIT(a,b,c,d))
#define r8_machar(a,b,c,d,e,f,g,h,i,j,k,l,m) FUNCNAME_R8MACHAR(C_9PLI4PIT(a,b,c,d,e,f,g,h,i,j,k,l,m))
#define explode(a,b,c) R_USHRT(FUNCNAME_EXPLODE(C_2ITDT(a,b,c)))
#define grid_2d(a,b,c,d,e,f,g,h) FUNCNAME_GRID2D(C_2DT4IT2PIT(a,d,b,c,e,f,g,h))
#define polygon_triangulate(a,b,c) FUNCNAME_POLYGONTRIANGULATE(C_DT2PIT(a,b,c))
#define polygon_integral_xx(a,b) R_DBL(FUNCNAME_POLYGONINTEGRALXX(C_DTPIT(a,b)))
#define polygon_integral_y(a,b) R_DBL(FUNCNAME_POLYGONINTEGRALY(C_DTPIT(a,b)))
#define polygon_integral_yy(a,b) R_DBL(FUNCNAME_POLYGONINTEGRALYY(C_DTPIT(a,b)))
#define diaedg(a,b,c,d,e,f,g,h) R_SHRT(FUNCNAME_DIAEDG(C_PDBL8(a,b,c,d,e,f,g,h)))
#define lrline(a,b,c,d,e,f,g) R_SHRT(FUNCNAME_LRLINE(C_PDBL7(a,b,c,d,e,f,g)))
#define pwl_interp_2d_scattered_value(a,b,c,d,e,f,g,h) FUNCNAME_PWLINTERP2DSCATTEREDVALUE(C_DT2PITDT2PIDTPIT(a,b,c,d,e,f,g,h))
#define r8tris2(a,b,c,d,e) R_UCHR(FUNCNAME_R8TRIS2(C_DTPIT3PI(a,b,c,d,e)))
#define swapec(a,b,c,d,e,f,g,h,i,j) R_SHRT(FUNCNAME_SWAPEC(C_DT3PIDTPITDT3PI(a,b,c,d,e,f,g,h,i,j)))
#define vbedg(a,b,c,d,e,f,g,h,i,j,k) FUNCNAME_VBEDG(C_2ITDTPITDT6PI(a,b,c,d,e,f,g,h,i,j,k))
#define r8vec_swtb(a,b,c) FUNCNAME_R8VECSWTB(C_DT2PIT(a,b,c))
#define r8vec_swtf(a,b,c,d) FUNCNAME_R8VECSWTF(C_DT3PIT(a,b,c,d))
#define sparse_grid_cfn_size(a,b) R_USHRT(FUNCNAME_SPARSEGRIDCFNSIZE(C_PUSHRT2(a,b)))
#define square_monomial(a,b,c) R_DBL(FUNCNAME_SQUAREMONOMIAL(C_PI2PIT(c,a,b)))
#define square_rule(a,b,c,d,e) FUNCNAME_SQUARERULE(C_2PITPI2PIT(a,b,c,d,e))
#define square_volume(a,b) R_DBL(FUNCNAME_SQUAREVOLUME(C_PPDBL2(a,b)))
#define angle_shift_deg(a,b) R_DBL(FUNCNAME_ANGLESHIFTDEG(C_PDBL2(a,b)))
#define angle_to_rgb(a) FUNCNAME_ANGLETORGB(C_SDBL(a))
#define axis_limits(a,b,c,d,e,f,g) FUNCNAME_AXISLIMITS(C_2ITDT3PITPDT(a,b,c,d,e,f,g))
#define bar_check(a) R_UCHR(FUNCNAME_BARCHECK(a))
#define bar_code(a) FUNCNAME_BARCODE(a)
#define bar_digit_code_left(a) FUNCNAME_BARDIGITCODELEFT(C_SUSHRT(a))
#define bar_digit_code_right(a) FUNCNAME_BARDIGITCODERIGHT(C_SUSHRT(a))
#define bmi_english(a,b,c) R_DBL(FUNCNAME_BMIENGLISH(C_PDBL3(a,b,c)))
#define bmi_metric(a,b) R_DBL(FUNCNAME_BMIMETRIC(C_PDBL2(a,b)))
#define euler_constant() R_DBL(FUNCNAME_EULERCONSTANT())
#define feet_to_meters(a) R_DBL(FUNCNAME_FEETTOMETERS(C_SDBL(a)))
#define gauss_sum(a,b,c,d,e,f) R_DBL(FUNCNAME_GAUSSSUM(C_2DT4PIT(a,b,c,d,e,f)))
#define grid1(a,b,c,d) FUNCNAME_GRID1(C_2DT2PIT(a,b,c,d))
#define grid1n(a,b,c,d,e) FUNCNAME_GRID1N(C_2DT2PIT(a,b,c,d,e))
#define grid2(a,b,c,d,e,f) FUNCNAME_GRID2(C_4DT2PIT(a,b,c,d,e,f))
#define grid2n(a,b,c,d,e,f) FUNCNAME_GRID2N(C_4DT2PIT(a,b,c,d,e,f))
#define grid3(a,b,c,d,e,f) FUNCNAME_GRID3(C_3DT3PIT(a,b,c,d,e,f))
#define grid3n(a,b,c,d,e,f,g,h) FUNCNAME_GRID3N(C_5DT3PIT(a,b,c,d,e,f,g,h))
#define grid4(a,b,c,d,e,f,g,h,i,j) FUNCNAME_GRID4(C_7DT3PIT(a,b,c,d,e,f,g,h,i,j))
#define grid4n(a,b,c,d,e,f,g,h,i,j,k,l) FUNCNAME_GRID4N(C_9PLI4PIT(a,b,c,d,e,f,g,h,i,j,k,l))
#define i2_reverse_bytes(a) R_SHRT(FUNCNAME_I2REVERSEBYTES(C_SSHRT(a)))
#define pounds_to_kilograms(a) R_DBL(FUNCNAME_POUNDSTOKILOGRAMS(C_SDBL(a)))
#define versine_pulse(a,b,c,d,e) R_DBL(FUNCNAME_VERSINEPULSE(C_PDBL5(a,b,c,d,e)))
#define asm_enum(a) R_INT(FUNCNAME_ASMENUM(C_SINT(a)))
#define asm_triangle(a,b) FUNCNAME_ASMTRIANGLE(C_DTPI(a,b))
#define bell(a,b) FUNCNAME_BELL(C_DTPI(a,b))
#define change_greedy(a,b,c,d,e) FUNCNAME_CHANGEGREEDY(C_2DT3PI(a,b,c,d,e))
#define change_next(a,b,c,d,e,f) FUNCNAME_CHANGENEXT(C_2DT3PIPB(a,b,c,d,e,f))
#define chinese_check(a,b) R_INT(FUNCNAME_CHINESECHECK(C_DTPI(a,b)))
#define congruence(a,b,c,d) R_INT(FUNCNAME_CONGRUENCE(C_3IPI(a,b,c,d)))
#define count_pose_random(a,b,c) FUNCNAME_COUNTPOSERANDOM(C_PPINT3(a,b,c))
#define debruijn(a,b,c) FUNCNAME_DEBRUIJN(C_2DTPI(a,b,c))
#define derange_enum(a) R_INT(FUNCNAME_DERANGEENUM(C_SUSHRT(a)))
#define derange_enum2(a,b) FUNCNAME_DERANGEENUM2(C_DTPI(a,b))
#define derange_enum3(a) R_USHRT(FUNCNAME_DERANGEENUM3(C_SUSHRT(a)))
#define derange0_back_candidate(a,b,c,d,e,f) FUNCNAME_DERANGE0BACKCANDIDATE(C_2DT4PI(a,c,b,d,e,f))
#define derange0_back_next(a,b,c) FUNCNAME_DERANGE0BACKNEXT(C_DTPIPB(a,b,c))
#define derange0_check(a,b) R_UCHR(FUNCNAME_DERANGE0CHECK(C_DTPI(a,b)))
#define derange0_weed_next(a,b,c,d,e) FUNCNAME_DERANGE0WEEDNEXT(C_DTPIPB2PDT(a,b,c,d,e))
#define digraph_arc_euler(a,b,c,d,e,f) FUNCNAME_DIGRAPHARCEULER(C_2DT4PI(a,b,c,d,e,f))
#define diophantine(a,b,c,d,e,f) FUNCNAME_DIOPHANTINE(C_2DT4PI(a,b,c,d,e,f))
#define diophantine_solution_minimize(a,b,c,d) FUNCNAME_DIOPHANTINESOLUTIONMINIMIZE(C_2I2PI(a,b,c,d))
#define equiv_next(a,b,c,d,e) FUNCNAME_EQUIVNEXT(C_DTPDT2PIPB(a,b,c,d,e))
#define equiv_next2(a,b,c) FUNCNAME_EQUIVNEXT2(C_DTPIPB(c,b,a))
#define equiv_random(a,b,c,d) FUNCNAME_EQUIVRANDOM(C_DTPIPDTPI(a,b,c,d))
#define euler_row(a,b) FUNCNAME_EULERROW(C_DTPI(a,b))
#define frobenius_number_order2(a,b) R_INT(FUNCNAME_FROBENIUSNUMBERORDER2(C_PINT2(a,b)))
#define gray_next(a,b,c,d) FUNCNAME_GRAYNEXT(C_DTPIPDTPI(a,b,c,d))
#define gray_rank(a) R_USHRT(FUNCNAME_GRAYRANK(C_SINT(a)))
#define gray_rank2(a) R_USHRT(FUNCNAME_GRAYRANK2(C_SINT(a)))
#define gray_unrank(a) R_INT(FUNCNAME_GRAYUNRANK(C_SUSHRT(a)))
#define gray_unrank2(a) R_INT(FUNCNAME_GRAYUNRANK2(C_SUSHRT(a)))
#define i4_partition_random(a,b,c,d,e,f) FUNCNAME_I4PARTITIONRANDOM(C_DT4PIPDT(a,b,c,d,e,f))
#define josephus(a,b,c) R_USHRT(FUNCNAME_JOSEPHUS(C_PUSHRT3(a,b,c)))
#define matrix_product_opt(a,b,c,d) FUNCNAME_MATRIXPRODUCTOPT(C_DTPIPDTPI(a,b,c,d))
#define moebius_matrix(a,b,c) FUNCNAME_MOEBIUSMATRIX(C_DT2PI(a,b,c))
#define nim_sum(a,b) R_UINT(FUNCNAME_NIMSUM(C_PUINT2(a,b)))
#define padovan(a,b) FUNCNAME_PADOVAN(C_DTPI(a,b))
#define pell_basic(a,b,c) FUNCNAME_PELLBASIC(C_DT2PI(a,b,c))
#define pell_next(a,b,c,d,e,f,g) FUNCNAME_PELLNEXT(C_DT4I2PI(a,b,c,d,e,f,g))
#define perrin(a,b) FUNCNAME_PERRIN(C_DTPI(a,b))
#define pord_check(a,b) R_UCHR(FUNCNAME_PORDCHECK(C_DTPI(a,b)))
#define power_series1(a,b,c,d) FUNCNAME_POWERSERIES1(C_DTIT2PIT(a,b,c,d))
#define power_series2(a,b,c) FUNCNAME_POWERSERIES2(C_DT2PIT(a,b,c))
#define power_series3(a,b,c,d) FUNCNAME_POWERSERIES3(C_DT3PIT(a,b,c,d))
#define power_series4(a,b,c,d) FUNCNAME_POWERSERIES4(C_DT3PIT(a,b,c,d))
#define r8mat_permanent(a,b) R_DBL(FUNCNAME_R8MATPERMANENT(C_DTPI(a,b)))
#define r8vec_backtrack(a,b,c,d,e,f,g,h) FUNCNAME_R8VECBACKTRACK(C_2DTPITPI2PDTPITPI(a,b,c,d,e,f,g,h))
#define rat_farey(a,b,c,d,e) FUNCNAME_RATFAREY(C_2DTPDT2PI(a,b,c,d,e))
#define rat_farey2(a,b,c) FUNCNAME_RATFAREY2(C_DT2PI(a,b,c))
#define schroeder(a,b) FUNCNAME_SCHROEDER(C_DTPI(a,b))
#define subset_gray_rank(a,b) R_USHRT(FUNCNAME_SUBSETGRAYRANK(C_DTPI(a,b)))
#define subset_gray_unrank(a,b,c) FUNCNAME_SUBSETGRAYUNRANK(C_2DTPI(a,b,c))
#define subset_lex_next(a,b,c,d,e) FUNCNAME_SUBSETLEXNEXT(C_3DT2PI(a,b,c,d,e))
#define subset_random(a,b,c) FUNCNAME_SUBSETRANDOM(C_DT2PI(a,b,c))
#define thue_binary_next(a,b) FUNCNAME_THUEBINARYNEXT(C_PDTPI(a,b))
#define thue_ternary_next(a,b) FUNCNAME_THUETERNARYNEXT(C_PDTPI(a,b))
#define vector_constrained_next(a,b,c,d,e,f) FUNCNAME_VECTORCONSTRAINEDNEXT(C_DT4PIPB(a,b,c,d,e,f))
#define vector_constrained_next2(a,b,c,d,e,f) FUNCNAME_VECTORCONSTRAINEDNEXT2(C_DT4PIPB(a,b,c,d,e,f))
#define vector_constrained_next3(a,b,c,d,e,f) FUNCNAME_VECTORCONSTRAINEDNEXT3(C_DT3PIPITPB(a,b,c,d,e,f))
#define vector_constrained_next4(a,b,c,d,e,f) FUNCNAME_VECTORCONSTRAINEDNEXT4(C_DTPIT2PIPITITPB(a,b,c,d,e,f))
#define vector_constrained_next5(a,b,c,d,e,f) FUNCNAME_VECTORCONSTRAINEDNEXT5(C_DTPI2DTPIPB(a,b,c,d,e,f))
#define vector_constrained_next6(a,b,c,d,e,f,g,h) FUNCNAME_VECTORCONSTRAINEDNEXT6(C_DTPIT3PI2ITPB(a,b,c,d,e,f,g,h))
#define vector_constrained_next7(a,b,c,d,e,f,g) FUNCNAME_VECTORCONSTRAINEDNEXT7(C_DTPIT2PI2ITPB(a,b,c,d,e,f,g))
#define vector_next(a,b,c,d,e) FUNCNAME_VECTORNEXT(C_DT3PIPB(a,b,c,d,e))
#define ytb_enum(a) R_USHRT(FUNCNAME_YTBENUM(C_SUSHRT(a)))
#define ytb_next(a,b,c,d) FUNCNAME_YTBNEXT(C_DT2PIPB(a,b,c,d))
#define ytb_random(a,b,c,d) FUNCNAME_YTBRANDOM(C_DT3PI(a,b,c,d))
#define simple_f(a,b,c) FUNCNAME_SIMPLEF(C_IT2PIT(a,b,c))
#define cbt_traverse(a) FUNCNAME_CBTTRAVERSE(C_SUSHRT(a))
#define triangle_order6_physical_to_reference(a,b,c,d) FUNCNAME_TRIANGLEORDER6PHYSICALTOREFERENCE(C_DT3PIT(b,c,d,a))
#define triangle_order6_reference_to_physical(a,b,c,d) FUNCNAME_TRIANGLEORDER6REFERENCETOPHYSICAL(C_DT3PIT(b,c,d,a))
#define vand1(a,b) FUNCNAME_VAND1(C_DTPIT(a,b))
#define vandermonde_value_1d(a,b,c,d) FUNCNAME_VANDERMONDEVALUE1D(C_2DT2PIT(a,c,b,d))
#define cg_gb(a,b,c,d,e,f) FUNCNAME_CGGB(C_3DT3PIT(a,b,c,d,e,f))
#define cg_ge(a,b,c,d) FUNCNAME_CGGE(C_DT3PIT(a,b,c,d))
#define rnorm(a,b,c) R_DBL(FUNCNAME_RNORM(C_PI2PIT(a,b,c)))
#define wshrt(a,b,c,d,e) FUNCNAME_WSHRT(C_DTPITDTPIPIT(c,a,d,e,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_FROBENIUSNUMBERORDER2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EXPLODE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MACHAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C4COS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C4SIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POLYVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PNPOLYNOMIALCOEFFICIENTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAUSS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RJACOBI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIAPHONYCOMPUTE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_INVERSEERROR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBYSHEV1EXACTNESS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBYSHEV2EXACTNESS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_COMPNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECCOMPONENTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSORTEDNEAREST0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGRANGEBASISFUNCTION1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGUERREMONOMIALQUADRATURE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LUERROR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECLINSPACE2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONINTEGRALXX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONINTEGRALXY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONINTEGRALY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONINTEGRALYY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8TRIS2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SWAPEC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VBEDG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSWTF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SQUAREMONOMIAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SQUARERULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BARCHECK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAUSSSUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ASMTRIANGLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BELL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHANGEGREEDY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHANGENEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHINESECHECK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_COUNTPOSERANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DEBRUIJN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DERANGEENUM2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DERANGE0BACKCANDIDATE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DERANGE0BACKNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DERANGE0CHECK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DERANGE0WEEDNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIGRAPHARCEULER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EQUIVNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EQUIVNEXT2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EQUIVRANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EULERROW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAYNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4PARTITIONRANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MATRIXPRODUCTOPT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOEBIUSMATRIX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PADOVAN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PERRIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PORDCHECK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POWERSERIES1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POWERSERIES2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POWERSERIES3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POWERSERIES4(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATPERMANENT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECBACKTRACK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RATFAREY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RATFAREY2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SCHROEDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETGRAYRANK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETGRAYUNRANK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETLEXNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUBSETRANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_THUEBINARYNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_THUETERNARYNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORCONSTRAINEDNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORCONSTRAINEDNEXT2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORCONSTRAINEDNEXT3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORCONSTRAINEDNEXT4(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORCONSTRAINEDNEXT5(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORCONSTRAINEDNEXT6(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORCONSTRAINEDNEXT7(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_YTBNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_YTBRANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEORDER6PHYSICALTOREFERENCE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEORDER6REFERENCETOPHYSICAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUAEQUAD0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CGGB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CGGE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FFT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FFT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECNORMAL01NEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4VECUNIFORMABNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POSL2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECUNIFORMABNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_JACOBI1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGRANGEINTERP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONTRIANGULATE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PWLINTERP2DSCATTEREDVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSWTB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLETORGB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BARCODE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BARDIGITCODELEFT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BARDIGITCODERIGHT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID1N(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID2N(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID3N(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID4(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRID4N(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VAND1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VANDERMONDEVALUE1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBYSHEV1INTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBYSHEV2INTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COSD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_COTD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CSCD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SECD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SIND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8TAND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8AI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8AID(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POCH1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8B1MP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8REN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8AIDE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8AIE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8AINT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ASINH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8BI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8BID(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8BIDE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8BIE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CHI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CINH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8COSDEG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8DAWSON(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8E1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8EI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ERF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ERFC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8EXPREL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8FAC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8GAMR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8INT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8LI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8LNGAM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8LNREL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8RAND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SHI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SINDEG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SPENCE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SQRT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8TAN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8TANH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SQUAREVOLUME(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FEETTOMETERS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I2REVERSEBYTES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POUNDSTOKILOGRAMS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ASMENUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DERANGEENUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DERANGEENUM3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAYRANK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAYRANK2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAYUNRANK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAYUNRANK2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_YTBENUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PNPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PPOLYNOMIALPRIME2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MODP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8BINOM(void *);	
__MATHSUITE __JBURKARDT void * FUNCNAME_R8GAMI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8GAMIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8GAMIT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8LBETA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MOD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8PAK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8POCH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8RANDGS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PMPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PMNPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PMNSPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CHU(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8CHUSCALED(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8GMIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8LGIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8LGIT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8GMIT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8ADMP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8AIMP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8LGAMS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8SIFG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8UPAK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIAEDG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LRLINE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPARSEGRIDCFNSIZE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLESHIFTDEG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_AXISLIMITS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BMIENGLISH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BMIMETRIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EULERCONSTANT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VERSINEPULSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CONGRUENCE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIOPHANTINE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIOPHANTINESOLUTIONMINIMIZE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_JOSEPHUS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NIMSUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PELLBASIC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PELLNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CBTTRAVERSE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RNORM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WSHRT(void *);

#endif
