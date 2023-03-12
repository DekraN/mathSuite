#ifndef WRAPPER_VALUES_H_INCLUDED
#define WRAPPER_VALUES_H_INCLUDED

#define student_noncentral_cdf_values(a,b,c,d,e) FUNCNAME_STUDENTNONCENTRALCDFVALUES(C_PDTPI3PIT(a,b,c,d,e))
#define normal_01_cdf_values(a,b,c) FUNCNAME_NORMAL01CDFVALUES(C_PDT2PIT(a,b,c))
#define gamma_values(a,b,c) FUNCNAME_GAMMAVALUES(C_PDT2PIT(a,b,c))
#define i4_fall_values(a,b,c,d) FUNCNAME_I4FALLVALUES(C_PPUSHRT4(a,b,c,d))
#define bessel_j0_values(a,b,c) FUNCNAME_BESSELJ0VALUES(C_PDT2PIT(a,b,c))
#define erfc_values(a,b,c) FUNCNAME_ERFCVALUES(C_PDT2PIT(a,b,c))
#define h_polynomial_values(a,b,c,d) FUNCNAME_HPOLYNOMIALVALUES(C_2PDT2PIT(a,b,c,d))
#define he_polynomial_values(a,b,c,d) FUNCNAME_HEPOLYNOMIALVALUES(C_2PDT2PIT(a,b,c,d))
#define hf_function_values(a,b,c,d) FUNCNAME_HFFUNCTIONVALUES(C_2PDT2PIT(a,b,c,d))
#define hep_values(a,b,c,d) FUNCNAME_HEPVALUES(C_2PDT2PIT(a,b,c,d))
#define hypersphere_01_area_values(a,b,c) FUNCNAME_HYPERSPHERE01AREAVALUES(C_2PDTPIT(a,b,c))
#define hypersphere_01_volume_values(a,b,c) FUNCNAME_HYPERSPHERE01VOLUMEVALUES(C_2PDTPIT(a,b,c))
#define i4_factorial2_values(a,b,c) FUNCNAME_I4FACTORIAL2VALUES(C_2PDTPU(a,b,c))
#define i4_rise_values(a,b,c,d) FUNCNAME_I4RISEVALUES(C_3PDTPU(a,b,c,d))
#define j_polynomial_values(a,b,c,d,e,f) FUNCNAME_JPOLYNOMIALVALUES(C_2PDT4PIT(a,b,c,d,e,f))
#define l_polynomial_values(a,b,c,d) FUNCNAME_LPOLYNOMIALVALUES(C_2PDT2PIT(a,b,c,d))
#define lf_function_values(a,b,c,d,e) FUNCNAME_LFFUNCTIONVALUES(C_PDT2PITPDTPIT(a,c,d,b,e))
#define pmn_polynomial_values(a,b,c,d,e) FUNCNAME_PMNPOLYNOMIALVALUES(C_3PDT2PIT(a,b,c,d,e))
#define pmns_polynomial_values(a,b,c,d,e) FUNCNAME_PMNSPOLYNOMIALVALUES(C_3PDT2PIT(a,b,c,d,e))
#define lobatto_polynomial_values(a,b,c,d) FUNCNAME_LOBATTOPOLYNOMIALVALUES(C_2PDT2PIT(a,b,c,d))
#define bernoulli_number_values(a,b,c) FUNCNAME_BERNOULLINUMBERVALUES(C_2PDTPIT(a,b,c))
#define bernstein_poly_values(a,b,c,d,e) FUNCNAME_BERNSTEINPOLYVALUES(C_3PDT2PIT(a,b,c,d,e))
#define beta_values(a,b,c,d) FUNCNAME_BETAVALUES(C_PDT3PIT(a,b,c,d))
#define catalan_values(a,b,c) FUNCNAME_CATALANVALUES(C_PPUSHRT3(a,b,c))
#define cheby_t_poly_values(a,b,c,d) FUNCNAME_CHEBYTPOLYVALUES(C_2PDT2PIT(a,b,c,d))
#define cheby_u_poly_values(a,b,c,d) FUNCNAME_CHEBYUPOLYVALUES(C_2PDT2PIT(a,b,c,d))
#define collatz_count_values(a,b,c) FUNCNAME_COLLATZCOUNTVALUES(C_PPUSHRT3(a,b,c))
#define cos_power_int_values(a,b,c,d,e) FUNCNAME_COSPOWERINTVALUES(C_PDT2PITPDTPIT(a,b,c,d,e))
#define euler_number_values(a,b,c) FUNCNAME_EULERNUMBERVALUES(C_2PDTPI(a,b,c))
#define gegenbauer_poly_values(a,b,c,d,e) FUNCNAME_GEGENBAUERPOLYVALUES(C_PDT2PITPDTPIT(a,c,d,b,e))
#define gud_values(a,b,c) FUNCNAME_GUDVALUES(C_PDT2PIT(a,b,c))
#define hyper_2f1_values(a,b,c,d,e,f) FUNCNAME_HYPER2F1VALUES(C_4PITPDTPIT(b,c,d,e,a,f))
#define lerch_values(a,b,c,d,e) FUNCNAME_LERCHVALUES(C_PDT2PITPDTPIT(a,b,d,c,e))
#define mertens_values(a,b,c) FUNCNAME_MERTENSVALUES(C_2PDTPS(a,b,c))
#define moebius_values(a,b,c) FUNCNAME_MOEBIUSVALUES(C_2PDTPS(a,b,c))
#define omega_values(a,b,c) FUNCNAME_OMEGAVALUES(C_2PDTPI(a,c,b))
#define partition_distinct_count_values(a,b,c) FUNCNAME_PARTITIONDISTINCTCOUNTVALUES(C_PPUSHRT3(a,b,c))
#define phi_values(a,b,c) FUNCNAME_PHIVALUES(C_PPUSHRT3(a,b,c))
#define r8_factorial_log_values(a,b,c) FUNCNAME_R8FACTORIALLOGVALUES(C_2PDTPIT(a,b,c))
#define sigma_values(a,b,c) FUNCNAME_SIGMAVALUES(C_PPUSHRT3(a,b,c))
#define sin_power_int_values(a,b,c,d,e) FUNCNAME_SINPOWERINTVALUES(C_PDT2PITPDTPIT(a,b,c,d,e))
#define spherical_harmonic_values(a,b,c,d,e,f,g) FUNCNAME_SPHERICALHARMONICVALUES(C_2PDTPS4PIT(a,b,c,d,e,f,g))
#define tau_values(a,b,c) FUNCNAME_TAUVALUES(C_PPUSHRT3(a,b,c))
#define zeta_values(a,b,c) FUNCNAME_ZETAVALUES(C_2PDTPIT(a,b,c))
#define bessel_i0_values(a,b,c) FUNCNAME_BESSELI0VALUES(C_PDT2PIT(a,b,c))
#define bessel_i1_values(a,b,c) FUNCNAME_BESSELI1VALUES(C_PDT2PIT(a,b,c))
#define bessel_ix_values(a,b,c,d) FUNCNAME_BESSELIXVALUES(C_PDT3PIT(a,b,c,d))
#define cauchy_cdf_values(a,b,c,d,e) FUNCNAME_CAUCHYCDFVALUES(C_PDT4PIT(a,b,c,d,e))
#define exponential_cdf_values(a,b,c,d) FUNCNAME_EXPONENTIALCDFVALUES(C_PDT3PIT(a,b,c,d))
#define extreme_values_cdf_values(a,b,c,d,e) FUNCNAME_EXTREMEVALUESCDFVALUES(C_PDT4PIT(a,b,c,d,e))
#define geometric_cdf_values(a,b,c,d) FUNCNAME_GEOMETRICCDFVALUES(C_2PDT2PIT(a,b,c,d))
#define laplace_cdf_values(a,b,c,d,e) FUNCNAME_LAPLACECDFVALUES(C_PDT4PIT(a,b,c,d,e))
#define log_normal_cdf_values(a,b,c,d,e) FUNCNAME_LOGNORMALCDFVALUES(C_PDT4PIT(a,b,c,d,e))
#define log_series_cdf_values(a,b,c,d) FUNCNAME_LOGSERIESCDFVALUES(C_PDTPITPDTPIT(a,c,b,d))
#define logistic_cdf_values(a,b,c,d,e) FUNCNAME_LOGISTICCDFVALUES(C_PDT4PIT(a,b,c,d,e))
#define rayleigh_cdf_values(a,b,c,d) FUNCNAME_RAYLEIGHCDFVALUES(C_PDT3PIT(a,b,c,d))
#define von_mises_cdf_values(a,b,c,d,e) FUNCNAME_VONMISESCDFVALUES(C_PDT4PIT(a,b,c,d,e))
#define weibull_cdf_values(a,b,c,d,e) FUNCNAME_WEIBULLCDFVALUES(C_PDT4PIT(a,b,c,d,e))
#define frobenius_number_order2_values(a,b,c,d) FUNCNAME_FROBENIUSNUMBERORDER2VALUES(C_PPUSHRT4(a,b,c,d))
#define r8_fall_values(a,b,c,d) FUNCNAME_R8FALLVALUES(C_PDTPITPDTPIT(a,c,b,d))
#define r8_rise_values(a,b,c,d) FUNCNAME_R8RISEVALUES(C_PDTPITPDTPIT(a,c,b,d))
#define lambert_w_values(a,b,c) FUNCNAME_LAMBERTWVALUES(C_PDT2PIT(a,b,c))
#define r8_factorial2_values(a,b,c) FUNCNAME_R8FACTORIAL2VALUES(C_2PDTPIT(a,b,c))
#define bivariate_normal_cdf_values(a,b,c,d,e) FUNCNAME_BIVARIATENORMALCDFVALUES(C_PDT4PIT(a,b,c,d,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_NORMAL01CDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAMMAVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSELJ0VALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ERFCVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HYPERSPHERE01AREAVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HYPERSPHERE01VOLUMEVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4FACTORIAL2VALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BERNOULLINUMBERVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CATALANVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_COLLATZCOUNTVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EULERNUMBERVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GUDVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MERTENSVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOEBIUSVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_OMEGAVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARTITIONDISTINCTCOUNTVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PHIVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8FACTORIALLOGVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIGMAVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TAUVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ZETAVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSELI0VALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSELI1VALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAMBERTWVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8FACTORIAL2VALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4FALLVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HPOLYNOMIALVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEPOLYNOMIALVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HFFUNCTIONVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEPVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4RISEVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOBATTOPOLYNOMIALVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_AGMVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETAVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBYTPOLYVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBYUPOLYVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BESSELIXVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EXPONENTIALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GEOMETRICCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOGSERIESCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RAYLEIGHCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FROBENIUSNUMBERORDER2VALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8FALLVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8RISEVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_STUDENTNONCENTRALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PMNPOLYNOMIALVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PMNSPOLYNOMIALVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BERNSTEINPOLYVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_COSPOWERINTVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GEGENBAUERPOLYVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LERCHVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SINPOWERINTVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CAUCHYCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EXTREMEVALUESCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAPLACECDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOGNORMALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOGISTICCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VONMISESCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_WEIBULLCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BIVARIATENORMALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LFFUNCTIONVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_JPOLYNOMIALVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LPOLYNOMIALVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HYPER2F1VALUES(void *);

#endif
