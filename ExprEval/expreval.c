/*
    File: expreval.c
    Auth: Brian Allen Vanderburg II
    Date: Wednesday, April 30, 2003
    Desc: Evaluation routines for the ExprEval library

    This file is part of ExprEval.
*/

/* Includes */
#include "exprincl.h"

#include "exprpriv.h"

/* Defines for error checking */
#include <errno.h>

#if(EXPR_ERROR_LEVEL >= EXPR_ERROR_LEVEL_CHECK)
#define EXPR_RESET_ERR() errno = 0
#define EXPR_CHECK_ERR() if(errno) err = EXPR_ERROR_OUTOFRANGE
#else
#define EXPR_RESET_ERR()
#define EXPR_CHECK_ERR()
#endif


/* This routine will evaluate an expression */
EXPRERRTYPE exprEval(exprObj *obj, EXPRTYPE val)
    {
    EXPRTYPE dummy;
    
    if(val == NULL)
        val = dummy;

    /* Make sure it was parsed successfully */
        /* Do NOT reset the break count.  Let is accumulate
           between calls until breaker function is called */
        return (!obj->parsedbad && obj->parsedgood && obj->headnode) ? exprEvalNode(obj, obj->headnode, 0, val, EXPREVAL_ALLOC) : EXPR_ERROR_BADEXPR;
    }

/* Evaluate a node */
EXPRERRTYPE exprEvalNode(exprObj *obj, exprNode *nodes, int curnode, EXPRTYPE val, const bool alloc)
    {
    uint64_t i;
    EXPRERRTYPE err;
    int pos;
    EXPRTYPE d1, d2;

	if(alloc)
		mpfr_init(val);
	
    if(obj == NULL || nodes == NULL)
        return EXPR_ERROR_NULLPOINTER;

    /* Update n to point to correct node */
    nodes += curnode;

    /* Check breaker count */
    if(obj->breakcur-- <= 0)
        {
        /* Reset count before returning */
        obj->breakcur = obj->breakcount;
            
        if(exprGetBreakResult(obj))
            return EXPR_ERROR_BREAK;
        }

    switch(nodes->type)
        {
        case EXPR_NODETYPE_MULTI:
            {
            /* Multi for multiple expressions in one string */
            for(pos = 0; pos < nodes->data.oper.nodecount; ++pos)
                {
                if(pos)
                	mpfr_clear(val);
                if((err = exprEvalNode(obj, nodes->data.oper.nodes, pos, val, !alloc)))
                    return err;
                }
            break;
            }

        case EXPR_NODETYPE_ADD:
            {
            /* Addition */
            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 0, d1, EXPREVAL_ALLOC)))
	            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 1, d2, EXPREVAL_ALLOC)))
	            {
	            	mpfr_add(val, d1, d2, MPFR_RNDN);
					mpfr_clears(d1, d2, NULL);       
	        	}
	            else
	            {
	  				mpfr_clears(d1, d2, NULL);      
	                return err;
	        	}	

            break;
            }

        case EXPR_NODETYPE_SUBTRACT:
            {
            /* Subtraction */
            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 0, d1, EXPREVAL_ALLOC)))
	            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 1, d2, EXPREVAL_ALLOC)))
	        	{
	            	mpfr_sub(val, d1, d2, MPFR_RNDN);
					mpfr_clears(d1, d2, NULL); 
	        	}
	            else
	            {
					mpfr_clears(d1, d2, NULL); 
	                return err;
	        	}

            break;
            }

        case EXPR_NODETYPE_MULTIPLY:
            {
            /* Multiplication */
            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 0, d1, EXPREVAL_ALLOC)))
	            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 1, d2, EXPREVAL_ALLOC)))
	            {
	            	mpfr_mul(val, d1, d2, MPFR_RNDN);
					mpfr_clears(d1, d2, NULL); 
	        	}
	            else
	            {
					mpfr_clears(d1, d2, NULL); 
	                return err;
	        	}

            break;
            }

        case EXPR_NODETYPE_DIVIDE:
            {
            /* Division */
            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 0, d1, EXPREVAL_ALLOC)))
	            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 1, d2, EXPREVAL_ALLOC)))
	                {
	                if(mpfr_zero_p(d2) == 0)
	                {
		            	mpfr_div(val, d1, d2, MPFR_RNDN);
						mpfr_clears(d1, d2, NULL); 
		        	}
	                else
	                    {
						mpfr_clears(d1, d2, NULL); 
	#if(EXPR_ERROR_LEVEL >= EXPR_ERROR_LEVEL_CHECK)
	                    return EXPR_ERROR_DIVBYZERO;
	#else
	                    mpfr_set_d(val, 0.00, MPFR_RNDN);
	                    return EXPR_ERROR_NOERROR;
	#endif
	                    }
	                }
	            else
	            {
	            	mpfr_clears(d1, d2, NULL); 
	                return err;
	            }

            break;
            }

        case EXPR_NODETYPE_EXPONENT:
            {
            /* Exponent */
            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 0, d1, EXPREVAL_ALLOC)))
	            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 1, d2, EXPREVAL_ALLOC)))
	                {
	                EXPR_RESET_ERR();
	                mpfr_pow(val, d1, d2, MPFR_RNDN);
					mpfr_clears(d1, d2, NULL); 
	                EXPR_CHECK_ERR();
	                }
	            else
	            {
					mpfr_clears(d1, d2, NULL); 
	                return err;
	      	  	}		

            break;
            }

        case EXPR_NODETYPE_NEGATE:
            {
            /* Negative value */
            if(!(err = exprEvalNode(obj, nodes->data.oper.nodes, 0, d1, EXPREVAL_ALLOC)))
            {
				mpfr_neg(val, d1, MPFR_RNDN);
				mpfr_clear(d1);
        	}
            else
            {
				mpfr_clear(d1);
                return err;
        	}
           
            break;
            }


        case EXPR_NODETYPE_VALUE:
            {
            /* Directly access the value */
            mpfr_set(val, nodes->data.value.value, MPFR_RNDN);
            break;
            }

        case EXPR_NODETYPE_VARIABLE:
            {
            /* Directly access the variable or constant */
            mpfr_set(val, *(nodes->data.variable.vaddr), MPFR_RNDN);
            break;
            }

        case EXPR_NODETYPE_ASSIGN:
            {
            /* Evaluate assignment subnode */
            if(!(err = exprEvalNode(obj, nodes->data.assign.node, 0, val, !alloc)))
                /* Directly assign the variable */
                mpfr_set(*(nodes->data.assign.vaddr), val, MPFR_RNDN);
            else
                return err;
            
            break;
            }

        default: // case EXPR_NODETYPE_FUNCTION:
            {
            	dim_typ j;
            	const bool isEager = nodes->data.function.mp_ctrl->type == CTRL_EVALTYPE && nodes->data.function.mp_ctrl->ReturnCTRLSystem.evalType == EAGER_EVALUATION;
            	EXPRTYPE args[isEager*nodes->data.function.nodecount];
            	
            	if(isEager || nodes->data.function.mp_ctrl->type != CTRL_EVALTYPE)
            		for(i=0; i<nodes->data.function.nodecount; ++i)	
						if((err = exprEvalNode(obj, nodes->data.function.nodes, i, args[i], EXPREVAL_ALLOC)))
						{
							for(j=0; j<=i; ++j)	
								mpfr_clear(args[j]);
							return err;
						}		
            	
            	EXPR_RESET_ERR();
                /* Call the correct function */
            	if(nodes->data.function.rt)
            		nodes->data.function.rt((*(nodes->data.function.fptr))(val, obj, nodes, &err, args), args, nodes, val);
            	else
					(*(nodes->data.function.fptr))(val, obj, nodes, &err, args);
                EXPR_CHECK_ERR();
                    
	            if(isEager)
					for(i=0; i<nodes->data.function.nodecount; ++i)
						mpfr_clear(args[i]);
						
				return err;
            }

		
        // default:
            /* Unknown node type */
            //return EXPR_ERROR_UNKNOWN;
        
        }

    return EXPR_ERROR_NOERROR;
    }



