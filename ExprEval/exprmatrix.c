#ifndef __DISABLE_EXPRMATRIX

#include "../dutils.h"

__MATHSUITE void * FUNCNAME_MATRIXDET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register dim_typ d1ex = sqrt(nodes->data.function.nodecount);
    const register dim_typ dim = d1ex*d1ex;

    if(dim == nodes->data.function.nodecount)
        det(val, args, d1ex, NULL);
    else
        return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
        
 	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MATRIXNORM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register dim_typ d1ex = sqrt(nodes->data.function.nodecount-1);
    const register dim_typ dim = d1ex*d1ex;

    if(dim == nodes->data.function.nodecount-1)
		norm(val, &args[1], d1ex, mpfr_get_ui(args[0], MPFR_RNDN));
    else
        return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MATRIXTRACE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register dim_typ d1ex = sqrt(nodes->data.function.nodecount);
    const register dim_typ dim = d1ex*d1ex;

    if(dim == nodes->data.function.nodecount)    
		_matrixTrace(val, args, d1ex);
    else
        return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
        
    return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MATRIXRANK(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register dim_typ d2ex = (nodes->data.function.nodecount-1) / mpfr_get_ui(args[0], MPFR_RNDN);
    const register dim_typ d1ex = (nodes->data.function.nodecount-1) / d2ex;
    const register dim_typ dim = d1ex*d2ex;

    if(d1ex == mpfr_get_ui(args[0], MPFR_RNDN))
        {

		    dim_typ2 dex =
		    {
		    	d1ex,
		    	d2ex
		    };
		
		    mpfr_set_ui(val, rank(&args[1], dex), MPFR_RNDN);
        
        }
    else 
        return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_MATRIXILLCHK(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register dim_typ d1ex = sqrt(nodes->data.function.nodecount);
    const register dim_typ dim = d1ex*d1ex;

    if(dim == nodes->data.function.nodecount)
        {

            
		dim_typ2 dims =
		{
			d1ex,
			d1ex
		};
		
        mpfr_t *matrix;

		if(!equalMatrix(&matrix, args, dims, EQUALMATRIX_REALLOC))
    		return (*err = EXPR_ERROR_MEMORY), NULL;
    	
		EXPRTYPE norms1;
		
        mpfr_init(norms1);
		norms(norms1, args, d1ex);
		
		if(!invertMatrix(matrix, args, d1ex))
		{
			matrixFree(&matrix, dims);
            return (*err = EXPR_ERROR_UNKNOWN), NULL;
        }
		
		norms(val, args, d1ex);
		mpfr_mul(val, val, norms1, MPFR_RNDN);
        matrixFree(&matrix, dims);
        mpfr_clear(norms1);
        }
    else
        return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL;
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE void * FUNCNAME_SCALARPROD(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
	dim_typ i, pos;
	
    const register dim_typ d1ex = nodes->data.function.nodecount;
    const register dim_typ dim = d1ex>>1; 
    
    if((dim<<1) == d1ex)
    {
        EXPRTYPE vector[d1ex];
        mpfr_set_ui(val, 0, MPFR_RNDN);
        mpfr_init(v);

        for(pos = 0; pos<d1ex; ++pos )
        {
            if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
            {
            	mpfr_clear(v);
            	return err;
            }
            else if(pos < dim)
        		mpfr_init_set(vector[pos], v, MPFR_RNDN);
            else
            {
				mpfr_mul(v, v, vector[pos-dim], MPFR_RNDN);
				mpfr_add(val, val, v, MPFR_RNDN);
			}
        }
		mpfr_clear(v);
    }
    else
		return (*err = EXPR_ERROR_BADNUMBERARGUMENTS), NULL; 
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

#endif
