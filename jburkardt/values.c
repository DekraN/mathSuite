#ifndef __DISABLE_BESSEL

#include "../../dutils.h"


__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_GAMMAVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	gamma_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSELJ0VALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	bessel_j0_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_ERFCVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	erfc_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_HYPERSPHERE01AREAVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	hypersphere_01_area_values(&ref, &ref2, &ref3);
	mpfr_set_d(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_HYPERSPHERE01VOLUMEVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	hypersphere_01_volume_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_I4FACTORIAL2VALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	unsigned ref3 = 0;
	i4_factorial2_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BERNOULLINUMBERVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	bernoulli_number_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_CATALANVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	catalan_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_COLLATZCOUNTVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	collatz_count_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_EULERNUMBERVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	int ref3 = 0;
	euler_number_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_si(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_GUDVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	gud_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_MERTENSVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	short ref3 = 0;
	mertens_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
    mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_si(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_MOEBIUSVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	short ref3 = 0;
	moebius_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_si(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_OMEGAVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	int ref2 = 0;
	dim_typ ref3 = 0;
	omega_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_si(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_PARTITIONDISTINCTCOUNTVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	partition_distinct_count_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_PHIVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	phi_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_R8FACTORIALLOGVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	r8_factorial_log_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_SIGMAVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	sigma_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_TAUVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	tau_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_ZETAVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	zeta_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSELI0VALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	bessel_i0_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSELI1VALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	bessel_i1_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_LAMBERTWVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	lambert_w_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_R8FACTORIAL2VALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	r8_factorial2_values(&ref, &ref2, &ref3);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_I4FALLVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	dim_typ ref4 = 0;
	i4_fall_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_HPOLYNOMIALVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	h_polynomial_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_HEPOLYNOMIALVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	he_polynomial_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_HFFUNCTIONVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	hf_function_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_HEPVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	hep_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_I4RISEVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	unsigned ref4 = 0;
	i4_rise_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_LOBATTOPOLYNOMIALVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	lobatto_polynomial_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_AGMVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	agm_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BETAVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	beta_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_CHEBYTPOLYVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	cheby_t_poly_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_CHEBYUPOLYVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	cheby_u_poly_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSELIXVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	bessel_ix_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_EXPONENTIALCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	exponential_cdf_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_GEOMETRICCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	geometric_cdf_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_LOGSERIESCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	dim_typ ref3 = 0;
	ityp ref4 = 0.00;
	log_series_cdf_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_RAYLEIGHCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	rayleigh_cdf_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_FROBENIUSNUMBERORDER2VALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	dim_typ ref4 = 0;
	frobenius_number_order2_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_R8FALLVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	dim_typ ref3 = 0;
	ityp ref4 = 0.00;
	r8_fall_values(&ref, &ref2, &ref3, &ref4);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_R8RISEVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	dim_typ ref3 = 0;
	ityp ref4 = 0.00;
	r8_rise_values(&ref, &ref2, &ref3, &ref4),
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_STUDENTNONCENTRALCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	int ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	student_noncentral_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_si(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_PMNPOLYNOMIALVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	pmn_polynomial_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_PMNSPOLYNOMIALVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	pmns_polynomial_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BERNSTEINPOLYVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	dim_typ ref3 = 0;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	bernstein_poly_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_COSPOWERINTVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	dim_typ ref4 = 0;
	ityp ref5 = 0.00;
	cos_power_int_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_GEGENBAUERPOLYVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	gegenbauer_poly_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_LERCHVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	dim_typ ref3 = 0;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	lerch_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_SINPOWERINTVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	dim_typ ref4 = 0;
	ityp ref5 = 0.00;
	sin_power_int_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_CAUCHYCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	cauchy_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_EXTREMEVALUESCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	extreme_values_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_LAPLACECDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	laplace_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_LOGNORMALCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	log_normal_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_LOGISTICCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	logistic_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_VONMISESCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	von_mises_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_WEIBULLCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	weibull_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BIVARIATENORMALCDFVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	bivariate_normal_cdf_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_LFFUNCTIONVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	lf_function_values(&ref, &ref2, &ref3, &ref4, &ref5);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_JPOLYNOMIALVALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	dim_typ ref2 = 0;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	ityp ref6 = 0.00;
	j_polynomial_values(&ref, &ref2, &ref3, &ref4, &ref5, &ref6);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[5]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_HYPER2F1VALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dim_typ ref = 0;
	ityp ref2 = 0.00;
	ityp ref3 = 0.00;
	ityp ref4 = 0.00;
	ityp ref5 = 0.00;
	ityp ref6 = 0.00;
	hyper_2f1_values(&ref, &ref2, &ref3, &ref4, &ref5, &ref6);
	mpfr_set_ui(val, true, MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[1]), ref2, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[2]), ref3, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[3]), ref4, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[4]), ref5, MPFR_RNDN);
	mpfr_set_d(*(nodes->data.function.refs[5]), ref5, MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

#endif
