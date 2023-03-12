#ifndef __DISABLE_EXPRTRIG

#include "../dutils.h"

__MATHSUITE void * FUNCNAME_SIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);
    	
	mpfr_sin(val, args[0], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_sinh(val, args[0], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_csc(val, args[0], MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
        
    mpfr_csch(val, args[0], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ASIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && !(TRIGONOMETRIC_DOMAIN(args[0])))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    
    mpfr_asin(val, args[0], MPFR_RNDN);
    if(isSett(BOOLS_DEGREESENTERING))
		mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ASINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_asinh(val, args[0], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ACSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_cmp_d(args[0], -1) > 0 && mpfr_cmp_d(args[0], 1) < 0)
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    }
    mpfr_acsc(val, args[0]);
    if(isSett(BOOLS_DEGREESENTERING))
		mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ACSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_acsch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);
	
        mpfr_cos(val, args[0], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_cosh(val, args[0], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);
    	
    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_sec(val, args[0], MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_sech(val, args[0], MPFR_RNDN);       
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ACOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && !(TRIGONOMETRIC_DOMAIN(args[0])))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL; 
    
    mpfr_acos(val, args[0], MPFR_RNDN);   
    if(isSett(BOOLS_DEGREESENTERING))
		mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ACOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_cmp_d(args[0], 1) < 0)
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    
    mpfr_acosh(val, args[0], MPFR_RNDN);   
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ASEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_cmp_d(args[0], -1) > 0 && mpfr_cmp_d(args[0], 1) < 0)
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    mpfr_asec(val, args[0]);
	if(isSett(BOOLS_DEGREESENTERING))
		mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ASECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && (mpfLZero(args[0]) || mpfr_cmp_d(args[0], 1) >= 0))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    mpfr_asech(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_TAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
		mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_tan(val, args[0], MPFR_RNDN); 
	mpfr_clears(mpftmp, mpftmp2, NULL);  
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_TANH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_tanh(val, args[0], MPFR_RNDN);        
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);
    	
    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_cot(val, args[0], MPFR_RNDN); 
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COTH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    
    mpfr_coth(val, args[0], MPFR_RNDN); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ATAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_atan(val, args[0], MPFR_RNDN); 
    if(isSett(BOOLS_DEGREESENTERING))
		mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ATANH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && !(TRIGONOMETRIC_DOMAIN(args[0])))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    
    mpfr_atanh(val, args[0], MPFR_RNDN); 
    if(isSett(BOOLS_DEGREESENTERING))	
		mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ACOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_acot(val, args[0]);
    if(isSett(BOOLS_DEGREESENTERING))
		mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ACOTH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && TRIGONOMETRIC_DOMAIN(args[0]))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    mpfr_acoth(val, args[0]);
    if(isSett(BOOLS_DEGREESENTERING))
		mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_hsin(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hsinh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);
    	
    mpfr_qsin(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qsinh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_hcos(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hcosh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_qcos(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qcosh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HSEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))	
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
	}
    mpfr_hsec(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HSECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hsech(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QSEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qsec(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QSECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qsech(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
		mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }			
	mpfr_clears(mpftmp, mpftmp2, NULL); 
    mpfr_hcsc(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hcsch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
	{
		mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qcsc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qcsch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HTAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_htan(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HTANH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_htanh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QTAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);
    	
    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qtan(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QTANH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qtanh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
	{
		mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hcot(val, args[0]);
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCOTH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hcoth(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
	{
		mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qcot(val, args[0]);
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCOTH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qcoth(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_VSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_vsin(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_VSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_vsinh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
		mpfr_rad(args[0], args[0]);

    mpfr_cvsin(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_cvsinh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_VCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{

	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_vcos(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_VCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_vcosh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{

	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_cvcos(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_cvcosh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

        mpfr_hvsin(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hvsinh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_hcvsin(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_hcvsinh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_qvsin(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qvsinh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	 
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_qcvsin(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qcvsinh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_hvcos(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hvcosh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_hcvcos(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hcvcosh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_qvcos(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qvcosh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_qcvcos(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qcvcosh(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ESEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
   
    if(isSett(BOOLS_DEGREESENTERING))	
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_esec(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ESECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_esech(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ECSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
	{
		mpfr_clears(args[0], mpftmp, mpftmp2, NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_ecsc(val, args[0]);
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ECSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_ecsch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HESEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
		mpfr_rad(args[0], args[0]);
		
    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hesec(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HESECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hesech(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HECSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
   
    if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
	{
		mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hecsc(val, args[0]);
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HECSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hecsch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QESEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING)) 	
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qesec(val, args[0]);
   	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QESECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qesech(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QECSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
	{
		mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qecsc(val, args[0]);
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QECSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qecsch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SINC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_sinc(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SINCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_sinch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HSINC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_hsinc(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HSINCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hsinch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QSINC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_qsinc(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QSINCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qsinch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_cosc(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_cosch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCOSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hcosc(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCOSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hcosch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCOSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qcosc(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCOSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qcosch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SECC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_secc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SECCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_secch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HSECC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hsecc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HSECCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hsecch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QSECC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{

	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qsecc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QSECCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qsecch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSCC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_cscc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSCCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_cscch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCSCC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 	
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hcscc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCSCCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_hcscch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCSCC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);


    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qcscc(val, args[0]);        	
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCSCCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qcscch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_TANC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);


    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_tanc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_TANCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_tanch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HTANC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_htanc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HTANCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_htanch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QTANC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);
    	
    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qtanc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QTANCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_qtanch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COTC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_cotc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COTCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_cotch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCOTC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(isSett(BOOLS_DEGREESENTERING)) 
    	mpfr_rad(args[0], args[0]);


    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hcotc(val, args[0]);
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HCOTCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_hcotch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCOTC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    
	if(isSett(BOOLS_DEGREESENTERING))
    	mpfr_rad(args[0], args[0]);

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && !trigonometric_domain(args[0], mpftmp, mpftmp2))
    {
    	mpfr_clears(args[0], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qcotc(val, args[0]);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QCOTCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]))
    {
    	mpfr_clear(args[0]);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    mpfr_qcotch(val, args[0]);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ATAN2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_atan2(val, args[0], args[1], MPFR_RNDN);
    if(isSett(BOOLS_DEGREESENTERING)) mpfr_deg(val, val);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = csin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
			 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = csinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = ccos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = ccosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
			
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CTAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ctan(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
		
	mpfr_clears(mpftmp, mpftmp2, NULL); 	
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CTANH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = ctanh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ccsc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
		
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ccsch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
			
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = csec(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = csech(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ccot(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCOTH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ccoth(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = chsin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chsinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}
	
    const register double complex result = cqsin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqsinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = chcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chcosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = cqcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqcosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHSEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chsec(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHSECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chsech(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQSEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqsec(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQSECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqsech(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chcsc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chcsch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqcsc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqcsch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
			
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHTAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chtan(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHTANH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chtanh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQTAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqtan(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQTANH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqtanh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
			
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chcot(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCOTH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chcoth(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqcot(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
	mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCOTH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[0]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqcoth(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CPXVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = cpxvsin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CPXVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cpxvsinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
		
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = ccvsin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = ccvsinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CPXVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = cpxvcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

    if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CPXVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	const register double complex result = cpxvcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
	mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
  
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = ccvcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = ccvcosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

        const register double complex result = chvsin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
        mpfr_set_d(val, creal(result), MPFR_RNDN);

		if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
			mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
		}
		else
			mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chvsinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = chcvsin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chcvsinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
  
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}
	
    const register double complex result = cqvsin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqvsinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCVSIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}
	
        const register double complex result = cqcvsin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
        mpfr_set_d(val, creal(result), MPFR_RNDN);

		if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
			mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
		}
		else
			mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCVSINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqcvsinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = chvcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chvcosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

        const register double complex result = chcvcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
        mpfr_set_d(val, creal(result), MPFR_RNDN);

		if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
			mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
		}
		else
			mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chcvcosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}
	
        const register double complex result = cqvcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
        mpfr_set_d(val, creal(result), MPFR_RNDN);

        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
			mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
		}
		else
			mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqvcosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
			
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCVCOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = cqcvcos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCVCOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqcvcosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CESEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cesec(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
        
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CESECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cesech(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CECSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cecsc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CECSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cecsch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHESEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chesec(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
        
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHESECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chesech(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHECSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = checsc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHECSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = checsch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQESEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqesec(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL);       
        
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQESECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqesech(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQECSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqecsc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQECSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqecsch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
	
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSINC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = csinc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSINCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = csinch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
		
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHSINC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    const register double complex result = chsinc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHSINCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chsinch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQSINC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}
	
    const register double complex result = cqsinc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQSINCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqsinch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCOSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}
	
    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ccosc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCOSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = ccosch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCOSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chcosc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCOSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chcosch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCOSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqcosc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCOSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqcosch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSECC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
           
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}
	
    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = csecc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CSECCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = csecch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHSECC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL);  
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chsecc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHSECCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chsecch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQSECC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{

    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqsecc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQSECCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqsecch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
		
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCSCC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ccscc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
	mpfr_clears(mpftmp, mpftmp2, NULL);   
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCSCCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = ccscch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCSCC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chcscc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCSCCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chcscch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCSCC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqcscc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCSCCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqcscch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CTANC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ctanc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CTANCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = ctanch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHTANC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL);  
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chtanc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHTANCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = chtanch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQTANC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_inits(mpftmp, mpftmp2, NULL); 
    mpfr_const_pi(mpftmp, MPFR_RNDN);
    mpfr_div_ui(mpftmp, mpftmp, 2, MPFR_RNDN);
	mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (mpfr_zero_p(args[0]) || !trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqtanc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQTANCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cqtanch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCOTC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ccotc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
   	mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CCOTCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = ccotch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCOTC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chcotc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CHCOTCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
   	 	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = chcotch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCOTC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		mpfr_rad(args[0], args[0]);
		mpfr_rad(args[1], args[1]);
	}

    mpfr_t mpftmp, mpftmp2;
    mpfr_init(mpftmp2);
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    mpfr_const_pi(mpftmp2, MPFR_RNDN);
    if(dcheck && mpfr_zero_p(args[1]) && (!trigonometric_domain(args[0], mpftmp, mpftmp2)))
    {
    	mpfr_clears(args[0], args[1], mpftmp, mpftmp2, NULL); 
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqcotc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
    mpfr_clears(mpftmp, mpftmp2, NULL); 

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CQCOTCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{

    if(dcheck && mpfr_zero_p(args[0]) && mpfr_zero_p(args[1]))
    {
    	mpfr_clears(args[0], args[1], NULL);
    	return (*err = EXPR_ERROR_DIVBYZERO), NULL;
    }
    const register double complex result = cqcotch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CASIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[1]) && !(TRIGONOMETRIC_DOMAIN(args[0])))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    const register double complex result = casin(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
    if(isSett(BOOLS_DEGREESENTERING)) 
		mpfr_deg(val, val);


	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CASINH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = casinh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CACOS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[1]) && !(TRIGONOMETRIC_DOMAIN(args[0])))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    	
    const register double complex result = cacos(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
    if(isSett(BOOLS_DEGREESENTERING)) 
		mpfr_deg(val, val);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
   
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CACOSH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[1]) && mpfr_cmp_d(args[0], 1) <0)
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    const register double complex result = cacosh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CATAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = catan(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
    if(isSett(BOOLS_DEGREESENTERING)) 
		mpfr_deg(val, val);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CATANH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[1]) && !(TRIGONOMETRIC_DOMAIN(args[0])))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    	
    const register double complex result = catanh(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CACSC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[1]) && mpfr_cmp_d(args[0], -1) > 0 && mpfr_cmp_d(args[0], 1) < 0)
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    	
    const register double complex result = cacsc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
    if(isSett(BOOLS_DEGREESENTERING)) 
		mpfr_deg(val, val);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CACSCH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cacsch(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CASEC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[1]) && mpfr_cmp_d(args[0], -1) > 0 && mpfr_cmp_d(args[0], 1) < 0)
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    	
    const register double complex result = casec(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
    if(isSett(BOOLS_DEGREESENTERING)) 
		mpfr_deg(val, val);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CASECH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[1]) && (mpfLZero(args[0]) || mpfr_cmp_d(args[0], 1) >= 0))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    const register double complex result = casech(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CACOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cacot(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
    if(isSett(BOOLS_DEGREESENTERING)) 
		mpfr_deg(val, val);

	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CACOTH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(dcheck && mpfr_zero_p(args[1]) && TRIGONOMETRIC_DOMAIN(args[0]))
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    	
    const register double complex result = cacoth(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
    mpfr_set_d(val, creal(result), MPFR_RNDN);
    
	if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set_d(*(nodes->data.function.refs[IMAG_PART]), cimag(result), MPFR_RNDN);
	}
	else
		mpfr_set_d(*(nodes->data.function.refs[REAL_PART]), cimag(result), MPFR_RNDN);
 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

#endif
