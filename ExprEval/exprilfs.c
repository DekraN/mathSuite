#ifndef __DISABLE_EXPRILFS

#include "../dutils.h"

__MATHSUITE void * FUNCNAME_ATLTUAE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	answer_to_life_the_universe_and_everything(val);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EXIT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    exit(mpfr_get_si(args[0], MPFR_RNDN));	
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MSS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    access(mss) = false;
    msprintf(COLOR_USER, "\nScripting Mode has been correctly disabled.\n\n");
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ABS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_abs(val, args[0], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MOD(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_fmod(val, args[0], args[1], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_IPART(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_modf(val, args[1], args[0], MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FPART(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_modf(val, args[1], args[0], MPFR_RNDN); 
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MIN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
     
    
    if((*err = exprEvalNode(obj, nodes->data.function.nodes, 0, val, EXPREVAL_NALLOC)))
    	return err;
    	
    mpfr_init(v);
    
	for(dim_typ pos = 1; pos < nodes->data.function.nodecount; ++pos)
    {
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
			mpfr_clear(v);
			return err;
        }
        else if(mpfr_less_p(v, val))
            mpfr_set(val, v, MPFR_RNDN);
    }
	
	mpfr_clear(v);
	// _MPFR_MIN(val, nodes->data.function.nodecount, args);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MAX(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
     
    
    if((*err = exprEvalNode(obj, nodes->data.function.nodes, 0, val, EXPREVAL_NALLOC)))
    	return err;
    	
    mpfr_init(v);
    
	for(dim_typ pos = 1; pos < nodes->data.function.nodecount; ++pos)
    {
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
			mpfr_clear(v);
			return err;
        }
        else if(mpfr_greater_p(v, val))
            mpfr_set(val, v, MPFR_RNDN);
    }
	
	mpfr_clear(v);
	// _MPFR_MIN(val, nodes->data.function.nodecount, args);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_BITCOUNTER(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
   	mpfr_set_si(val, mpfr_zero_p(args[2]) ? countbits(mpfr_get_si(args[0], MPFR_RNDN)) : ucountbits(mpfr_get_si(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MGET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	equalSpecialMatrix(SPECIALMATRIX_CURRENTMATRIX, SPECIALMATRIX_LMPMATRIX);
	mpfr_set_ui(val, true, MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MLET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(mpfr_cmp_ui(args[0], 0) >= 0 && mpfr_cmp_ui(args[1], 0) >= 0 && mpfr_cmp_ui(args[1], MAX_SPECIALMATRICES+accessCurrentSession()->MLSystem.itemsno) < 0 && mpfr_cmp_ui(args[1], MAX_SPECIALMATRICES+accessCurrentSession()->MLSystem.itemsno) < 0)
	{
		equalSpecialMatrix(mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN));
		mpfr_set_ui(val, true, MPFR_RNDN);
	}
	else
	{
		mpfr_set_ui(val, false, MPFR_RNDN);
		return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
	}

 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MDEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    msyprintf(COLOR_SYSTEM, "\nIt has been deleted the entire Matrix List.\n\n");
	delAllMatrixFromMatrixList();
	mpfr_set_ui(val, true, MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MOFF(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_add_ui(val, args[1], MAX_SPECIALMATRICES, MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SOL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	randomizedSelect(val, &args[1], 0, nodes->data.function.nodecount-2, mpfr_get_ui(args[0], MPFR_RNDN));
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_VERSION(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_d(val, strtod(PROGRAM_VERSION, NULL), MPFR_RNDN);
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EXITCHAR(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_si(args[1], 'B') || mpfr_cmp_si(args[1], 'A'))
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[1], MPFR_RNDN);
        access(curLayout)->exit_char = mpfr_get_si(val, MPFR_RNDN);
    }
    else
    	mpfr_set_si(val, access(curLayout)->exit_char, MPFR_RNDN);

    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_PREC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_PRECISION) < 0 || mpfr_cmp_ui(args[0], MAX_PRECISION) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->precision = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
    	mpfr_set_ui(val, access(curLayout)->precision, MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_STABFACT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_STABFACT) < 0 || mpfr_cmp_ui(args[0], MAX_STABFACT) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    	mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->stabilizer_factor = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
        mpfr_set_ui(val, access(curLayout)->stabilizer_factor, MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_BLOCKSIZE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_BLOCKSIZE) < 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->block_size = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
        mpfr_set_ui(val, access(curLayout)->block_size, MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MINOSMMDIM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_OSMM_DIM) < 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->min_osmm_dim = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
    	mpfr_set_ui(val, access(curLayout)->min_osmm_dim, MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MINSTRASSENDIM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_STRASSEN_DIM) < 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->min_strassen_dim = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
        mpfr_set_ui(val, access(curLayout)->min_strassen_dim, MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MINSRNUMBER(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_MEMOIZABLE_INDEX+1) < 0 || mpfr_cmp_ui(args[0], MIN_STIRLING_NUMBER) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->min_stirling_number = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
   		mpfr_set_ui(val, access(curLayout)->min_stirling_number, MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ALGEBRA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_ALGEBRA) < 0 || mpfr_cmp_ui(args[0], MAX_ALGEBRA) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->algebra = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
    	mpfr_set_ui(val, access(curLayout)->algebra, MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_OUTLIERCONST(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_d(args[0], MIN_OUTLIER_CONSTANT) < 0 || mpfr_cmp_d(args[0], MAX_OUTLIER_CONSTANT) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->outlier_constant = mpfr_get_flt(val, MPFR_RNDN);
    }
    else
    	mpfr_set_flt(val, access(curLayout)->outlier_constant, MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RSEED(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_RANDOMSEED) < 0 || mpfr_cmp_ui(args[0], MAX_RANDOMSEED) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(random_seed) = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
    	mpfr_set_ui(val, access(random_seed), MPFR_RNDN);
         	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MMIFIBO(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_MEMOIZABLE_INDEX) < 0 || mpfr_cmp_ui(args[0], MAX_FIBONACCI_MEMOIZABLE_INDEX) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
    	mpfr_set_ui(val, access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI], MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MMIEVENDFACT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_MEMOIZABLE_INDEX) < 0 || mpfr_cmp_ui(args[0], MAX_EVEN_DOUBLEFACTORIAL_MEMOIZABLE_INDEX) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] = mpfr_get_ui(val, MPFR_RNDN);
    }
    else
    	mpfr_set_ui(val, access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL], MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MMIODDDFACT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(nodes->data.function.nodecount)
    {
        if(mpfr_cmp_ui(args[0], MIN_MEMOIZABLE_INDEX) < 0 || mpfr_cmp_ui(args[0], MAX_ODD_DOUBLEFACTORIAL_MEMOIZABLE_INDEX) > 0)
            return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
        mpfr_set(val, args[0], MPFR_RNDN);
        access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] = mpfr_get_ui(val, MPFR_RNDN);

    }
    else
		mpfr_set_ui(val, access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL], MPFR_RNDN);
        	
    return (*err = EXPR_ERROR_NOERROR), NULL;
} 

__MATHSUITE void * FUNCNAME_BINSUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    char * const c = binaryAlgSum(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), BAS_SUM);
    
    if(c)
    {
		mpfr_set_si(val, strtod(c, NULL), MPFR_RNDN);
		free(c);
	}
	else
		return (*err = EXPR_ERROR_BADEXPR), NULL;
		
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_BINSUB(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    char * const c = binaryAlgSum(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), BAS_SUB);
    
    if(c)
    {
        mpfr_set_si(val, c ? strtod(c, NULL) : false, MPFR_RNDN);
        free(c);
    }
    else
        return (*err = EXPR_ERROR_BADEXPR), NULL;

    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COMPLEMENT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(mpfr_zero_p(args[1]))
        mpfr_set_d(val, ~((int64_t)mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
    else
    {
        char * const c = binNumComp(mpfr_get_d(args[0], MPFR_RNDN));
        if(c)
        {
            mpfr_set_si(val, strtod(c, NULL), MPFR_RNDN);
            free(c);
        }
        else
            return (*err = EXPR_ERROR_BADEXPR), NULL;
    }
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_POW(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_pow(val, args[0], args[1], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SQRT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_sqrt(val, args[0], MPFR_RNDN); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CBRT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_cbrt(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ROOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE tmp;
    mpfr_init(tmp);
    mpfr_pow_si(tmp, args[1], -1, MPFR_RNDN);
    mpfr_pow(val, args[0], tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HARRIS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    harris(val, args[0], nodes->data.function.nodecount-1, &args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_HARRIS2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	harris2(val, args[0], &args[1]); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_THARRIS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_ui(val, t_harris(args[0], args[1], nodes->data.function.nodecount-2, &args[2]), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_THARRIS2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_ui(val, t_harris2(args[0], args[1], &args[2]), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_R8NORMAL01VALUES(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	int ref = 0;
	mpfr_set_d(val, r8_normal_01(&ref), MPFR_RNDN);
	mpfr_set_ui(*(nodes->data.function.refs[0]), ref, MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

typedef struct
{
	int a; 
	int b;
	ityp c;
} testfunc_ARGS;

static __MATHSUITE ityp testfunc(void * data)
{
	testfunc_ARGS * testfunc_data = data;
	int a = testfunc_data->a;
	int b = testfunc_data->b;
	ityp c = testfunc_data->c;
	return a+b;
}

__MATHSUITE void * FUNCNAME_TEST(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    //printf("\n\nTEST OF MPFR\n\n");
		       // body of the test
	mpfr_t x, y, z, k;
	mpfr_init_set_d(x, 2.00, MPFR_RNDN); 
	mpfr_init_set_d(y, 500.00, MPFR_RNDN);      
	mpfr_init_set_d(k, 20.54, MPFR_RNDN);
	mpfr_init(z);
	
	mpfr_pow(z, x, y, MPFR_RNDN);  
	mpfr_prec_t pc = mpfr_get_prec(z);
	
	mpfr_printf("Z value: %Rf, Z prec: %Pu bits", z, pc);
	/* When the program is about to exit, do ... */
	// mpfr_set_d(val, testfunc(&((testfunc_ARGS){mpfr_get_si(x, MPFR_RNDN), mpfr_get_si(y, MPFR_RNDN)})), MPFR_RNDN);
	mpfr_set_d(val, testfunc(&((void * [3]){(void *)mpfr_get_si(x, MPFR_RNDN), (void *)mpfr_get_si(y, MPFR_RNDN), (void *) mpfr_get_si(k, MPFR_RNDN)})), MPFR_RNDN);
	mpfr_clear (x);
	mpfr_clear (y);
	mpfr_clear (z);
	mpfr_clear (k);
	mpfr_free_cache ();           /* free the cache for constants like pi */
	
	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_log10(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_log2(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_POW10(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_d(val, 10.00, MPFR_RNDN);
    mpfr_pow(val, val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_log(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EXP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_exp(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EXPC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_expc(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EXP10(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_exp10(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EXP10C(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_exp10c(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EXP2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_exp2(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EXP2C(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_exp2c(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}



__MATHSUITE void * FUNCNAME_LOGN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE l1, l2; 
	
	mpfr_inits(l1, l2, NULL);
    errno = 0;
    mpfr_log(l1, args[0], MPFR_RNDN);
    
    if(errno)
    {
    	mpfr_clear(l1);
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    }
    
    mpfr_log(l2, args[1], MPFR_RNDN); 
    
    if(errno)
    {
    	mpfr_clears(l1, l2, NULL);
    	return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
    }

    if(mpfr_zero_p(l2))
        {
#if(EXPR_ERROR_LEVEL >= EXPR_ERROR_LEVEL_CHECK)
        return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
#else
		mpfr_set_d(val, 0.00, MPFR_RNDN);
        return (*err = EXPR_ERROR_NOERROR), NULL;
#endif
        }

	mpfr_div(val, l1, l2, MPFR_RNDN);
	mpfr_clears(l1, l2, NULL);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOGC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_log10c(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LNC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_logc(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG2C(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_log2c(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG1P(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_log1p(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG1PC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_log1pc(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG101P(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_log101p(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG101PC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_log101pc(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG21P(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_log21p(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LOG21PC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_log21pc(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFREXPM1(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_expm1(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFREINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_eint(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRLI2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_li2(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRGAMMA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_gamma(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRLNGAMMA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_lngamma(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRDIGAMMA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_digamma(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRZETA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_zeta(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRERF(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_erf(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRERFC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_erfc(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRJ0(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_j0(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRJ1(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_j1(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRJN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_jn(val, mpfr_get_ui(args[0], MPFR_RNDN), args[1], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRY0(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_y0(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRY1(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_y1(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRYN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_yn(val, mpfr_get_ui(args[0], MPFR_RNDN), args[1], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRFMA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_fma(val, args[1], args[2], args[3], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRFMS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_fms(val, args[1], args[2], args[3], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRAGM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_agm(val, args[0], args[1], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRHYPOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_hypot(val, args[0], args[1], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MPFRAI(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_ai(val, args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CELFAH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_fah(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FAHCEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_fah(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CELKEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_kel(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_KELCEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_kel(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CELRANK(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_rank(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RANKCEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_rank(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CELREA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_rea(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_REACEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_rea(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CELNEW(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_new(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_NEWCEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_new(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CELDEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_del(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_DELCEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_del(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CELROM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_rom(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ROMCEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	cel_rom(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FAHKEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	fah_kel(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_KELFAH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	fah_kel(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FAHRANK(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	fah_rank(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RANKFAH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	fah_rank(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FAHREA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	fah_rea(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_REAFAH(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	fah_rea(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_REARANK(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	rea_rank(val, args[0], DIRECT_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RANKREA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	rea_rank(val, args[0], INVERSE_CONVERSION);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CEXP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cexp(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CEXPC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cexpc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CEXP10(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cexp10(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CEXP10C(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cexp10c(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CEXP2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cexp2(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CEXP2C(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = cexp2c(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CPOW(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	const register double complex result = cpow(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I, mpfr_get_d(args[2], MPFR_RNDN)+mpfr_get_d(args[3], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CROOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	const register double complex result = crootnX(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I, mpfr_get_d(args[2], MPFR_RNDN)+mpfr_get_d(args[3], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CSQRT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = csqrt(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CCBRT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = ccbrt(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOGN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	const register double complex result = clogbN(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I, mpfr_get_d(args[2], MPFR_RNDN)+mpfr_get_d(args[3], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLNC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clogc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog10(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOGC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog10c(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog2(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG2C(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog2c(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG1P(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog1p(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG1PC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog1pc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG101P(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog101p(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG101PC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog101pc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG21P(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog21p(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CLOG21PC(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register double complex result = clog21pc(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I);
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

__MATHSUITE void * FUNCNAME_CARG(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_d(val, isSett(BOOLS_DEGREESENTERING) ? degd(carg(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I)) : carg(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CABS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_d(val, cabs(mpfr_get_d(args[0], MPFR_RNDN)+mpfr_get_d(args[1], MPFR_RNDN)*I), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QABS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	_qabs(val, args);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_OABS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	_oabs(val, args);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SABS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	_sabs(val, args);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CEIL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_ceil(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FLOOR(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_floor(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SECONDGRADEQSOLVER(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{   	
   	EXPRTYPE root[2];
   	mpfr_inits(root[0], root[1], NULL); 
   	
    if(!_secondGradeEquationSolver(args, root))
       return (*err = EXPR_ERROR_NOERROR), NULL;

	mpfr_set(*(nodes->data.function.refs[ROOT_X1]), root[ROOT_X1], MPFR_RNDN);
	mpfr_set(val, root[ROOT_X1], MPFR_RNDN);
	mpfr_set(*(nodes->data.function.refs[ROOT_X2]), root[ROOT_X2], MPFR_RNDN);
	mpfr_clears(root[ROOT_X1], root[ROOT_X2], NULL);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COMPLEXADD(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{   
    EXPRTYPE complexRes[MAX_COMPLEX_UNITS];
    mpfr_inits(complexRes[REAL_PART], complexRes[IMAG_PART], NULL); 
    _complexAdd(args, complexRes);
    mpfr_set(val, complexRes[REAL_PART], MPFR_RNDN);
    
    if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set(*(nodes->data.function.refs[IMAG_PART]), complexRes[IMAG_PART], MPFR_RNDN);
	}
	else
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), complexRes[IMAG_PART], MPFR_RNDN);
			
	mpfr_clears(complexRes[REAL_PART], complexRes[IMAG_PART], NULL);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COMPLEXMUL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE complexRes[MAX_COMPLEX_UNITS];
    mpfr_inits(complexRes[REAL_PART], complexRes[IMAG_PART], NULL); 
    _complexMul(args, complexRes);
    mpfr_set(val, complexRes[REAL_PART], MPFR_RNDN);
    
    if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
	{
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), val, MPFR_RNDN);
		mpfr_set(*(nodes->data.function.refs[IMAG_PART]), complexRes[IMAG_PART], MPFR_RNDN);
	}
	else
		mpfr_set(*(nodes->data.function.refs[REAL_PART]), complexRes[IMAG_PART], MPFR_RNDN);	
				
    mpfr_clears(complexRes[REAL_PART], complexRes[IMAG_PART], NULL);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QUATERNIONSADD(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE quaternionsRes[MAX_QUATERNIONS_UNITS];
    
	mpfr_inits(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_KPART], NULL);
    _quaternionsAdd(args, quaternionsRes);
	mpfr_set(*(nodes->data.function.refs[REAL_PART]), quaternionsRes[QUATERNIONS_REALPART], MPFR_RNDN);
	mpfr_set(val, quaternionsRes[QUATERNIONS_REALPART], MPFR_RNDN);
    mpfr_set(*(nodes->data.function.refs[QUATERNIONS_IPART]), quaternionsRes[QUATERNIONS_IPART], MPFR_RNDN);
    mpfr_set(*(nodes->data.function.refs[QUATERNIONS_JPART]), quaternionsRes[QUATERNIONS_JPART], MPFR_RNDN);
    mpfr_set(*(nodes->data.function.refs[QUATERNIONS_KPART]), quaternionsRes[QUATERNIONS_KPART], MPFR_RNDN);
    mpfr_clears(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_KPART], NULL);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_QUATERNIONSMUL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE quaternionsRes[MAX_QUATERNIONS_UNITS];
	
	mpfr_inits(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_KPART], NULL);
    _quaternionsMul(args, quaternionsRes);
	mpfr_set(*(nodes->data.function.refs[REAL_PART]), quaternionsRes[QUATERNIONS_REALPART], MPFR_RNDN);
	mpfr_set(val, quaternionsRes[QUATERNIONS_REALPART], MPFR_RNDN);
    mpfr_set(*(nodes->data.function.refs[QUATERNIONS_IPART]), quaternionsRes[QUATERNIONS_IPART], MPFR_RNDN);
    mpfr_set(*(nodes->data.function.refs[QUATERNIONS_JPART]), quaternionsRes[QUATERNIONS_JPART], MPFR_RNDN);
    mpfr_set(*(nodes->data.function.refs[QUATERNIONS_KPART]), quaternionsRes[QUATERNIONS_KPART], MPFR_RNDN);
    mpfr_clears(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_KPART], NULL);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RAND(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    long a;

    /* Perform random routine directly */
    a = ((long)(mpfr_get_si(*(nodes->data.function.refs[0]), MPFR_RNDN))) * 214013L + 2531011L;
    mpfr_set_d(*(nodes->data.function.refs[0]), (ityp)a, MPFR_RNDN);
    mpfr_set_d(val, (ityp)((a >> 16) & 0x7FFF) / (ityp)(32768), MPFR_RNDN); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RANDOM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register ityp diff = mpfr_get_d(args[1], MPFR_RNDN) - mpfr_get_d(args[0], MPFR_RNDN);
    /* Perform random routine directly */
    const register long a = ((long)(mpfr_get_si(*(nodes->data.function.refs[0]), MPFR_RNDN))) * 214013L + 2531011L;
    mpfr_set_d(*(nodes->data.function.refs[0]), (ityp)a, MPFR_RNDN);
    const register ityp rval = (ityp)((a >> 16) & 0x7FFF) / (ityp)(32767);
	mpfr_add_d(val, args[0], rval*diff, MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RANDOMIZE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    static int curcall = 0;
    mpfr_set_d(*(nodes->data.function.refs[0]), (ityp)((clock() + 1024 + (++curcall) * time(NULL))), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_DEG(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_deg(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RAD(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_rad(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RECTTOPOLR(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE tmp;
    EXPRTYPE tmp2;
    EXPRTYPE tmp3; 
    mpfr_inits(tmp, tmp2, tmp3, NULL); 
    mpfr_pow_ui(tmp, args[0], 2, MPFR_RNDN);
    mpfr_pow_ui(tmp2, args[1], 2, MPFR_RNDN);
    mpfr_add(tmp3, tmp, tmp2, MPFR_RNDN);
    mpfr_sqrt(val, tmp3, MPFR_RNDN);
    mpfr_clears(tmp, tmp2, tmp3, NULL); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_RECTTOPOLA(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE tmp;
    mpfr_init(tmp);
    mpfr_atan2(tmp, args[1], args[0], MPFR_RNDN);
    
    if(mpfr_cmp_d(tmp, 0.00))
    {
    	mpfr_const_pi(val, MPFR_RNDN);
    	mpfr_mul_ui(val, val, 2, MPFR_RNDN);
    }
    else
    	mpfr_set(val, tmp, MPFR_RNDN);
    	
    mpfr_clear(tmp);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_POLTORECTX(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE tmp;
    mpfr_init(tmp);
    mpfr_cos(tmp, args[1], MPFR_RNDN);
    mpfr_mul(val, args[0], tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_POLTORECTY(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE tmp;
    mpfr_init(tmp);
    mpfr_sin(tmp, args[1], MPFR_RNDN);
    mpfr_mul(val, args[0], tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_BASECHANGE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_si(val, changeBase(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN)), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_NPRIMENUMBER(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	N_prime_Number(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_PRIMORIAL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_primr(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FIRSTNPRIMENUMBERSSUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_fpnsum(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FIBONACCIAL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_fibnc(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_LCM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	math_lcm(val, args[0], args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_GCD(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	math_GCD(val, args[0], args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_GEOMETRICSUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_gsum(val, args[0], args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FACTORIAL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_fac_ui(val, mpfr_get_ui(args[0], MPFR_RNDN), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_DOUBLEFACTORIAL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_dfact(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_STIRLING(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	stirling(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FIBONACCI(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_fibo(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_PERMUTATIONS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	perm(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_PERMUTATIONSREP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    // Overheaded but it is necessary if Domain Check is enabled.
    mpfr_t tmp, tmp2;
    mpfr_init(tmp);
    mpfr_init_set_ui(tmp2, nodes->data.function.nodecount-1, MPFR_RNDN);
    mpfr_summation(tmp, tmp2, SUMMATION_SUM, args);
    
    if(dcheck && mpfr_cmp(tmp, args[0]))
    	return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
    
	perm_rep(val, args[0], tmp2, &args[1]); 
	mpfr_clears(tmp, tmp2, NULL);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_KPERMUTATIONS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	kperm(val, args[0], args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_KPERMUTATIONSREP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	kperm_rep(val, args[0], args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COMBINATIONS(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	comb(val, args[0], args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COMBINATIONSREP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	comb_rep(val, args[0], args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ARMONICSUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_asum(val, args[0]); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_GENARMONICSUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_gasum(val, args[0], args[1]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FIBONACCISUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_fsum(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FACTORIALSUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_fasum(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_DOUBLEFACTORIALSUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_sfasum(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FIRSTNNUMBERSSUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_fnnsum(val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
	
	mpfr_init(v);
	mpfr_set_ui(val, 0, MPFR_RNDN);

    for(dim_typ pos = 0; pos < nodes->data.function.nodecount; ++pos)
    	if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
			mpfr_clear(v);
			return err;
        }
        else
            mpfr_add(val, val, v, MPFR_RNDN);
            
	mpfr_clear(v);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_PRODUCT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
	
	mpfr_init(v);
	mpfr_set_ui(val, 1, MPFR_RNDN);

    for(dim_typ pos = 0; pos < nodes->data.function.nodecount; ++pos)
    	if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
			mpfr_clear(v);
			return err;
        }
        else
            mpfr_mul(val, val, v, MPFR_RNDN);
            
	mpfr_clear(v);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MODE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
	math_mode(val, tmp, args); 
	mpfr_clear(tmp);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_VARIANCE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
    math_variance(val, tmp, args);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_VARIANCE2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
    math_variance2(val, tmp, args);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COVARIANCE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    if(!(nodes->data.function.nodecount & 1))
    	return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
    
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount >> 1, MPFR_RNDN); 
	math_covariance(val, tmp, args, &args[nodes->data.function.nodecount >> 1]);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_COVARIANCE2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(!(nodes->data.function.nodecount & 1))
    	return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
    
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount >> 1, MPFR_RNDN); 
	math_covariance2(val, tmp, args, &args[nodes->data.function.nodecount >> 1]);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_STDDEV(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
	math_stddev(val, tmp, args);
 	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_STDDEV2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
	math_stddev2(val, tmp, args);
 	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_STDCOD(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(!(nodes->data.function.nodecount & 1))
    	return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;

    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
	math_stdcod(val, tmp, args, &args[nodes->data.function.nodecount >> 1]);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_STDCOD2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(!(nodes->data.function.nodecount & 1))
    	return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
    
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount >> 1, MPFR_RNDN); 
	math_stdcod(val, tmp, args, &args[nodes->data.function.nodecount >> 1]);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_PEARSON(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(!(nodes->data.function.nodecount & 1))
    	return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
    
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount >> 1, MPFR_RNDN); 
	math_pearson(val, tmp, args, &args[nodes->data.function.nodecount >> 1]);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_PEARSON2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(!(nodes->data.function.nodecount & 1))
    	return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
    
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, nodes->data.function.nodecount >> 1, MPFR_RNDN); 
	math_pearson2(val, tmp, args, &args[nodes->data.function.nodecount >> 1]);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_OUTLIER(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	const register dim_typ dim = nodes->data.function.nodecount-1;
	
	if(mpfr_get_ui(args[1], MPFR_RNDN) >= dim)
		return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
	
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, dim, MPFR_RNDN); 
    mpfr_set_ui(val, math_outlier(tmp, args[0] , &args[1]), MPFR_RNDN);
    mpfr_clear(tmp);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_OUTLIER2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	const register dim_typ dim = nodes->data.function.nodecount-2;
	
	if(mpfr_get_ui(args[1], MPFR_RNDN) >= dim)
		return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
	
    mpfr_t tmp;
	mpfr_init_set_ui(tmp, dim, MPFR_RNDN); 
    mpfr_set_ui(val, math_outlier2(tmp, args[0], mpfr_get_flt(args[1], MPFR_RNDN), &args[2]), MPFR_RNDN);
    mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MAP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	const register dim_typ funcID = mpfr_get_ui(args[1], MPFR_RNDN);
   
    if(mpfr_cmp_ui(args[1], funcID) || funcID < FID_SIN || funcID >= MAX_FIDS)
    {
		mpfr_set_ui(val, false, MPFR_RNDN);
		return (*err = EXPR_ERROR_OUTOFRANGE), NULL;
	}
	
    ext_math.functions[funcID](val, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_GEOMEAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
  	mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
	math_geomean(val, tmp, args);
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ARMEAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
  	mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
	math_armean(val, tmp, args);
	mpfr_clear(tmp);  
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_POWMEAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
  	mpfr_init_set_ui(tmp, nodes->data.function.nodecount-1, MPFR_RNDN); 
	math_powmean(val, tmp, args[0], &args[1]);
	mpfr_clear(tmp);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CENTRALVALUE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_t tmp;
  	mpfr_init_set_ui(tmp, nodes->data.function.nodecount, MPFR_RNDN); 
	math_scale(val, tmp, args); 
	mpfr_clear(tmp); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FIRSTQUARTILE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    math_first_quartile(val, nodes->data.function.nodecount, args);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MEDIAN(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	math_median(val, nodes->data.function.nodecount, args);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_THIRDQUARTILE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	math_third_quartile(val, nodes->data.function.nodecount, args);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_IF(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v, d1;
	*err = exprEvalNode(obj, nodes->data.function.nodes, 0, v, EXPREVAL_ALLOC);
	
    if(!err)
    {
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, 1+mpfr_zero_p(args[0]), d1, EXPREVAL_ALLOC)))
        {
			mpfr_clears(v, d1, NULL);
            return err;
        }
        mpfr_set(val, d1, MPFR_RNDN);
        mpfr_clears(v, d1, NULL);
    }
    else
    {
		mpfr_clear(v);
        return err;
    }
    
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SELECT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v, d1;
	*err = exprEvalNode(obj, nodes->data.function.nodes, 0, v, EXPREVAL_ALLOC);
	
    if(!err)
    {
	    if(mpfLZero(args[0]))
        {
	        if((*err = exprEvalNode(obj, nodes->data.function.nodes, 1, d1, EXPREVAL_ALLOC)))
	        {
	        	mpfr_clears(v, d1, NULL);
	            return err;
	        }
        }
    	else if(mpfr_zero_p(args[0]))
        {
	        if((*err = exprEvalNode(obj, nodes->data.function.nodes, 2, d1, EXPREVAL_ALLOC)))
	        {
	        	mpfr_clears(v, d1, NULL);
	            return err;
	        }
        }
    	else
        {
	        if(nodes->data.function.nodecount == 3)
            {
	       	    if((*err = exprEvalNode(obj, nodes->data.function.nodes, 2, d1, EXPREVAL_ALLOC)))
	            {
		        	mpfr_clears(v, d1, NULL);
		            return err;
	        	}
            }
        	else
            {
	            if((*err = exprEvalNode(obj, nodes->data.function.nodes, 3, d1, EXPREVAL_ALLOC)))
	            {
		        	mpfr_clears(v, d1, NULL);
		            return err;
	        	}
            }
        }
    	mpfr_set(val, args[1], MPFR_RNDN);
    	mpfr_clears(v, d1, NULL);
	}
    else
    {
		mpfr_clear(v);
        return err;
    }
    
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_EQUAL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_si(val, !mpfr_cmp(args[0], args[1]), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_ABOVE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_ui(val, mpfr_cmp(args[0],args[1]) > 0, MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_BELOW(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_ui(val, mpfr_cmp(args[0],args[1]) < 0, MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_AVG(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	
	mpfr_init(v);
	mpfr_set_ui(val, 0, MPFR_RNDN);

    for(dim_typ pos = 0; pos < nodes->data.function.nodecount; ++pos)
    {
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
			mpfr_clear(v);
			return err;
        }
        mpfr_add(val, val, v, MPFR_RNDN);
    }

	mpfr_clear(v);
    mpfr_div_ui(val, val, (nodes->data.function.nodecount), MPFR_RNDN);   
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CLIP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set(val, mpfr_cmp(args[0], args[1]) < 0 ? args[1] : mpfr_cmp(args[0], args[2]) > 0 ? args[2] : args[0], MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_CLAMP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE tmp, tmp2, tmp3;
    mpfr_inits(tmp, tmp2, tmp3, NULL);  
    mpfr_sub(tmp2, args[0], args[1], MPFR_RNDN);
    mpfr_sub(tmp3, args[2], args[1], MPFR_RNDN);
    mpfr_fmod(tmp, tmp2, tmp3, MPFR_RNDN);
    mpfr_add(val, tmp, args[1+mpfLZero(tmp)], MPFR_RNDN);
	mpfr_clears(tmp, tmp2, tmp3, NULL); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_PNTCHANGE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    EXPRTYPE odiff, ndiff, perc;
    
    mpfr_inits(odiff, ndiff, NULL);
    mpfr_sub(odiff, args[1], args[0], MPFR_RNDN);
    mpfr_sub(ndiff, args[3], args[2], MPFR_RNDN);

    if(mpfr_zero_p(odiff))
    {
        mpfr_set(val, args[0], MPFR_RNDN);
		mpfr_clears(odiff, ndiff, NULL); 
        return (*err = EXPR_ERROR_NOERROR), NULL;
    }

	mpfr_init(perc);
	mpfr_sub(perc, args[4], args[0], MPFR_RNDN); 
	mpfr_div(perc, perc, odiff, MPFR_RNDN);
	mpfr_mul(val, perc, ndiff, MPFR_RNDN);
	mpfr_add(val, val, args[2], MPFR_RNDN);
	mpfr_clears(odiff, ndiff, perc, NULL); 
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_POLY(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	eval(val, &args[1], nodes->data.function.nodecount-1, args[0]);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_AND(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_ui(val, (!mpfr_zero_p(args[0])) && (!mpfr_zero_p(args[1])), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_OR(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    mpfr_set_ui(val, (!mpfr_zero_p(args[0])) || (!mpfr_zero_p(args[1])), MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_NOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{

    mpfr_set_ui(val, mpfr_zero_p(args[0]), MPFR_RNDN);	
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_FOR(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	dim_typ pos;
    EXPRTYPE v, test;
    
    mpfr_inits(v, test, NULL);
    *err = exprEvalNode(obj, nodes->data.function.nodes, 0, v, EXPREVAL_NALLOC);

    if(!err)
        *err = exprEvalNode(obj, nodes->data.function.nodes, 1, test, EXPREVAL_NALLOC);
        
    if(!err)
    {
		while(!mpfr_zero_p(test))
	    {
		    for(pos = 3; pos < nodes->data.function.nodecount; ++pos)
		        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, val, EXPREVAL_NALLOC)))
		        {
	    			mpfr_clears(v, test, NULL);
		            return err;
		        }

		    if((*err = exprEvalNode(obj, nodes->data.function.nodes, 2, v, EXPREVAL_NALLOC)))
		    {
	    		mpfr_clears(v, test, NULL);
		        return err;
		    }

		    if((*err = exprEvalNode(obj, nodes->data.function.nodes, 1, test, EXPREVAL_NALLOC)))
		    {
	    		mpfr_clears(v, test, NULL);
		        return err;
		    }
	    }
	    mpfr_clears(v, test, NULL);
    }
    else
    {
 		mpfr_clears(v, test, NULL);
        return err;
    }
    
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MANY(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	// nop
	// mpfr_set_ui(val, true, MPFR_RNDN);
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

#endif
