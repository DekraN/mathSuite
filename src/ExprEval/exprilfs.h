/*
    File: exprilfs.h
    Auth: Brian Allen Vanderburg II
    Date: Tuesday, February 28, 2006
    Desc: Inline Function Solvers for exprEvalNode

    This file is part of ExprEval.
*/

/*
    This is here to help prevent expreval.c from getting
    too crowded.

    Provided variables:
    obj: expression object point
    nodes: function node with paramters
    d1, d2: variables
    err: error
    val: value pointer for resuld
    pos: integer

    Also EXPR_RESET_ERR() and EXPR_CHECK_ERR()

    The chunks below are included inside a statement that looks like this:

    switch(nodes->data.function.type)
        {
        #include "exprilfs.h"

        default:
            {
            return EXPR_ERROR_UNKNOWN;
            }
        }
*/

/* exit */
case EXPR_NODEFUNC_EXIT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);

    if(!err)
        safeExit(*val);
    else
        return err;

    break;;
    }

/* mss */
case EXPR_NODEFUNC_MSS:
    access(mss) = false;
    printf2(COLOR_USER, "\nScripting Mode has been correctly disabled.\n\n");
    break;

/* abs */
case EXPR_NODEFUNC_ABS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        *val = (-1+((d1 >= 0)<<1))*d1;
    else
        return err;

    break;
    }

/* mod */
case EXPR_NODEFUNC_MOD:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = fmod(d1, d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* ipart */
case EXPR_NODEFUNC_IPART:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        modf(d1, val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* fpart */
case EXPR_NODEFUNC_FPART:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = modf(d1, &d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* min */
case EXPR_NODEFUNC_MIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                {
                if(d2 < d1)
                    d1 = d2;
                }
            else
                return err;
            }
        }
    else
        return err;

    *val = d1;

    break;
    }

/* max */
case EXPR_NODEFUNC_MAX:
    {
    int pos;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                {
                if(d2 > d1)
                    d1 = d2;
                }
            else
                return err;
            }
        }
    else
        return err;

    *val = d1;

    break;
    }

/* bcnt */
case EXPR_NODEFUNC_BITCOUNTER:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        *val = d2 ? ucountbits(d1) : countbits(d1);
    else
        return err;

    break;
    }

/* vers */
case EXPR_NODEFUNC_VERSION:
    {
    *val = strtod(PROG__VERSION, NULL);
    break;
    }
    
/* ec */
case EXPR_NODEFUNC_EXITCHAR:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) == 'A' || (*val) == 'B')
                return err;
            access(curLayout)->exit_char = (*val);
            }
        }
        else
            (*val) = access(curLayout)->exit_char;
    break;
    }

/* prec */
case EXPR_NODEFUNC_PREC:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_PRECISION || (*val) > MAX_PRECISION)
                return err;
            access(curLayout)->precision = *val;
            }
        }
    else
        (*val) = access(curLayout)->precision;

    break;
    }

/* stabfact */
case EXPR_NODEFUNC_STABFACT:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_STABFACT || (*val) > MAX_STABFACT)
                return err;
            access(curLayout)->stabilizer_factor = *val;
            }
        }
        else
            (*val) = access(curLayout)->stabilizer_factor;
    break;
    }

/* msrn */
case EXPR_NODEFUNC_MINSRNUMBER:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_MEMOIZABLE_INDEX+1 || (*val) > MIN_STIRLING_NUMBER)
                return err;
            access(curLayout)->min_stirling_number = *val;
            }
        }
        else
            (*val) = access(curLayout)->min_stirling_number;
    break;
    }

/* alg */
case EXPR_NODEFUNC_ALGEBRA:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_ALGEBRA || (*val) > MAX_ALGEBRA)
                return err;
            access(curLayout)->algebra = (*val);
            }
        }
        else
            (*val) = access(curLayout)->algebra;
    break;
    }
    
/* oc */
case EXPR_NODEFUNC_OUTLIERCONST:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_OUTLIER_CONSTANT || (*val) > MAX_OUTLIER_CONSTANT)
                return err;
            access(curLayout)->outlier_constant = (*val);
            }
        }
        else
            (*val) = access(curLayout)->outlier_constant;
    break;
    }

/* rseed */
case EXPR_NODEFUNC_RSEED:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_RANDOMSEED || (*val) > MAX_RANDOMSEED)
                return err;
            access(random_seed) = *val;
            }
        }
        else
            (*val) = access(random_seed);
    break;
    }

/* mmi_fibo */
case EXPR_NODEFUNC_MMIFIBO:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_MEMOIZABLE_INDEX || (*val) > MAX_FIBONACCI_MEMOIZABLE_INDEX)
                return err;
            access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = *val;
            }
        }
        else
            (*val) = access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI];
    break;
    }

/* mmi_fact */
case EXPR_NODEFUNC_MMIFACT:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_MEMOIZABLE_INDEX || (*val) > MAX_FATTORIALE_MEMOIZABLE_INDEX)
                return err;
            access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = *val;
            }
        }
        else
            (*val) = access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE];
    break;
    }

/* mmi_esfact */
case EXPR_NODEFUNC_MMIEVENSFACT:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_MEMOIZABLE_INDEX || (*val) > MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX)
                return err;
            access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = *val;
            }
        }
        else
            (*val) = access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE];
    break;
    }

/* mmi_osfact */
case EXPR_NODEFUNC_MMIODDSFACT:
    {
    if(nodes->data.function.nodecount)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 0, val);
        if(!err)
            {
            if((*val) < MIN_MEMOIZABLE_INDEX || (*val) > MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX)
                return err;
            access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] = *val;
            }
        }
        else
            (*val) = access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE];
    break;
    }


/* bsum */
case EXPR_NODEFUNC_BINSUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        char *c = NULL;
        if((c = binaryAlgSum(d1, d2, BAS_SUM)))
            {
            *val = strtod(c, NULL);
            free(c);
            }
        else
            return err;
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* bsub */
case EXPR_NODEFUNC_BINSUB:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        char *c = NULL;
        if((c = binaryAlgSum(d1, d2, BAS_SUB)))
            {
            *val = strtod(c, NULL);
            free(c);
            }
        else
            return err;
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* comp */
case EXPR_NODEFUNC_COMP:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    d2 = 0.00;

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        if(!d2)
            *val = ~((int64_t)d1);
        else
            {
            char *c = NULL;
            if((c = binNumComp(d1)))
                {
                *val = strtod(c, NULL);
                free(c);
                }
            else
                return err;
            }
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* pow */
case EXPR_NODEFUNC_POW:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = pow(d1, d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sqrt */
case EXPR_NODEFUNC_SQRT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sqrt(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cbrt */
case EXPR_NODEFUNC_CBRT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cbrt(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* root */
case EXPR_NODEFUNC_ROOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = pow(d1, (1/d2));
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sin */
case EXPR_NODEFUNC_SIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sinh */
case EXPR_NODEFUNC_SINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* csc */
case EXPR_NODEFUNC_CSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = csc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* csch */
case EXPR_NODEFUNC_CSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = csch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* asin */
case EXPR_NODEFUNC_ASIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = asin(d1);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* asinh */
case EXPR_NODEFUNC_ASINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = asinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* acsc */
case EXPR_NODEFUNC_ACSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = acsc(d1);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* acsch */
case EXPR_NODEFUNC_ACSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = acsch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cos */
case EXPR_NODEFUNC_COS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cosh */
case EXPR_NODEFUNC_COSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sec */
case EXPR_NODEFUNC_SEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sec(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sech */
case EXPR_NODEFUNC_SECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sech(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* acos */
case EXPR_NODEFUNC_ACOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = acos(d1);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* acosh */
case EXPR_NODEFUNC_ACOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = acosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* asec */
case EXPR_NODEFUNC_ASEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = asec(d1);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* asech */
case EXPR_NODEFUNC_ASECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = asech(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* tan */
case EXPR_NODEFUNC_TAN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = tan(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* tanh */
case EXPR_NODEFUNC_TANH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = tanh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cot */
case EXPR_NODEFUNC_COT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cot(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* coth */
case EXPR_NODEFUNC_COTH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = coth(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* atan */
case EXPR_NODEFUNC_ATAN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = atan(d1);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* atanh */
case EXPR_NODEFUNC_ATANH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = atanh(d1);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* acot */
case EXPR_NODEFUNC_ACOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = acot(d1);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* acoth */
case EXPR_NODEFUNC_ACOTH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = acoth(d1);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hsin */
case EXPR_NODEFUNC_HSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hsinh */
case EXPR_NODEFUNC_HSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hsinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qsin */
case EXPR_NODEFUNC_QSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qsinh */
case EXPR_NODEFUNC_QSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qsinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcos */
case EXPR_NODEFUNC_HCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcosh */
case EXPR_NODEFUNC_HCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcos */
case EXPR_NODEFUNC_QCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcosh */
case EXPR_NODEFUNC_QCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hsec */
case EXPR_NODEFUNC_HSEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hsec(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hsech */
case EXPR_NODEFUNC_HSECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hsech(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qsec */
case EXPR_NODEFUNC_QSEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qsec(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qsech */
case EXPR_NODEFUNC_QSECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qsech(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcsc */
case EXPR_NODEFUNC_HCSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcsc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcsch */
case EXPR_NODEFUNC_HCSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcsch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcsc */
case EXPR_NODEFUNC_QCSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcsc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcsch */
case EXPR_NODEFUNC_QCSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcsch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* htan */
case EXPR_NODEFUNC_HTAN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = htan(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* htanh */
case EXPR_NODEFUNC_HTANH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = htanh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qtan */
case EXPR_NODEFUNC_QTAN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qtan(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qtanh */
case EXPR_NODEFUNC_QTANH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qtanh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcot */
case EXPR_NODEFUNC_HCOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcot(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcoth */
case EXPR_NODEFUNC_HCOTH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcoth(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcot */
case EXPR_NODEFUNC_QCOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcot(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcoth */
case EXPR_NODEFUNC_QCOTH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcoth(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* vsin */
case EXPR_NODEFUNC_VSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = vsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* vsinh */
case EXPR_NODEFUNC_VSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = vsinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cvsin */
case EXPR_NODEFUNC_CVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cvsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cvsinh */
case EXPR_NODEFUNC_CVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cvsinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* vcos */
case EXPR_NODEFUNC_VCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = vcos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* vcosh */
case EXPR_NODEFUNC_VCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = vcosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cvcos */
case EXPR_NODEFUNC_CVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cvcos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cvcosh */
case EXPR_NODEFUNC_CVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cvcosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hvsin */
case EXPR_NODEFUNC_HVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hvsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hvsinh */
case EXPR_NODEFUNC_HVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hvsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcvsin */
case EXPR_NODEFUNC_HCVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcvsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcvsinh */
case EXPR_NODEFUNC_HCVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcvsinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qvsin */
case EXPR_NODEFUNC_QVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qvsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qvsinh */
case EXPR_NODEFUNC_QVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qvsinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcvsin */
case EXPR_NODEFUNC_QCVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcvsin(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcvsinh */
case EXPR_NODEFUNC_QCVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcvsinh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hvcos */
case EXPR_NODEFUNC_HVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hvcos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hvcosh */
case EXPR_NODEFUNC_HVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hvcosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcvcos */
case EXPR_NODEFUNC_HCVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcvcos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcvcosh */
case EXPR_NODEFUNC_HCVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcvcosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qvcos */
case EXPR_NODEFUNC_QVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qvcos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qvcosh */
case EXPR_NODEFUNC_QVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qvcosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcvcos */
case EXPR_NODEFUNC_QCVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);


    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcvcos(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcvcosh */
case EXPR_NODEFUNC_QCVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcvcosh(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* esec */
case EXPR_NODEFUNC_ESEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = esec(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* esech */
case EXPR_NODEFUNC_ESECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = esec(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* ecsc */
case EXPR_NODEFUNC_ECSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = ecsc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* ecsch */
case EXPR_NODEFUNC_ECSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = ecsch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hesec */
case EXPR_NODEFUNC_HESEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hesec(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hesech */
case EXPR_NODEFUNC_HESECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hesech(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hecsc */
case EXPR_NODEFUNC_HECSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hecsc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hecsch */
case EXPR_NODEFUNC_HECSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hecsch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qesec */
case EXPR_NODEFUNC_QESEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qesec(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qesech */
case EXPR_NODEFUNC_QESECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qesech(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qecsc */
case EXPR_NODEFUNC_QECSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qecsc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qecsch */
case EXPR_NODEFUNC_QECSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qecsch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sinc */
case EXPR_NODEFUNC_SINC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sinc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sinch */
case EXPR_NODEFUNC_SINCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sinch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hsinc */
case EXPR_NODEFUNC_HSINC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hsinc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hsinch */
case EXPR_NODEFUNC_HSINCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hsinch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qsinc */
case EXPR_NODEFUNC_QSINC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qsinc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qsinch */
case EXPR_NODEFUNC_QSINCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qsinch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cosc */
case EXPR_NODEFUNC_COSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cosc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cosch */
case EXPR_NODEFUNC_COSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cosch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcosc */
case EXPR_NODEFUNC_HCOSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcosc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcosch */
case EXPR_NODEFUNC_HCOSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcosch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcosc */
case EXPR_NODEFUNC_QCOSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcosc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcosch */
case EXPR_NODEFUNC_QCOSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcosch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* secc */
case EXPR_NODEFUNC_SECC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = secc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* secch */
case EXPR_NODEFUNC_SECCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = secch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hsecc */
case EXPR_NODEFUNC_HSECC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hsecc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hsecch */
case EXPR_NODEFUNC_HSECCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hsecch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qsecc */
case EXPR_NODEFUNC_QSECC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qsecc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qsecch */
case EXPR_NODEFUNC_QSECCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qsecch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cscc */
case EXPR_NODEFUNC_CSCC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cscc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cscch */
case EXPR_NODEFUNC_CSCCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cscch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcscc */
case EXPR_NODEFUNC_HCSCC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcscc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcscch */
case EXPR_NODEFUNC_HCSCCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcscch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcscc */
case EXPR_NODEFUNC_QCSCC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcscc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcscch */
case EXPR_NODEFUNC_QCSCCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcscch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* tanc */
case EXPR_NODEFUNC_TANC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = tanc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* tanch */
case EXPR_NODEFUNC_TANCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = tanch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* htanc */
case EXPR_NODEFUNC_HTANC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = htanc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* htanch */
case EXPR_NODEFUNC_HTANCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = htanch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qtanc */
case EXPR_NODEFUNC_QTANC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qtanc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qtanch */
case EXPR_NODEFUNC_QTANCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qtanch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cotc */
case EXPR_NODEFUNC_COTC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cotc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cotch */
case EXPR_NODEFUNC_COTCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cotch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcotc */
case EXPR_NODEFUNC_HCOTC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcotc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* hcotch */
case EXPR_NODEFUNC_HCOTCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = hcotch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcotc */
case EXPR_NODEFUNC_QCOTC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
    if(isSett(BOOLS_DEGREESENTERING)) d1 = rad(d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcotc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* qcotch */
case EXPR_NODEFUNC_QCOTCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = qcotch(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* atan2 */
case EXPR_NODEFUNC_ATAN2:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = atan2(d1, d2);
        if(isSett(BOOLS_DEGREESENTERING)) *val = deg(*val);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* csin */
case EXPR_NODEFUNC_CSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csin(d1+d2*I);
        *val = creal(result);
		if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* csinh */
case EXPR_NODEFUNC_CSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccos */
case EXPR_NODEFUNC_CCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccosh */
case EXPR_NODEFUNC_CCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ctan */
case EXPR_NODEFUNC_CTAN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ctan(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ctanh */
case EXPR_NODEFUNC_CTANH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ctanh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccsc */
case EXPR_NODEFUNC_CCSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccsc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccsch */
case EXPR_NODEFUNC_CCSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccsch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* csec */
case EXPR_NODEFUNC_CSEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csec(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* csech */
case EXPR_NODEFUNC_CSECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csech(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccot */
case EXPR_NODEFUNC_CCOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccot(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccoth */
case EXPR_NODEFUNC_CCOTH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccoth(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chsin */
case EXPR_NODEFUNC_CHSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chsin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }    
    
/* chsinh */
case EXPR_NODEFUNC_CHSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chsinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqsin */
case EXPR_NODEFUNC_CQSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqsin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqsinh */
case EXPR_NODEFUNC_CQSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqsinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcos */
case EXPR_NODEFUNC_CHCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcosh */
case EXPR_NODEFUNC_CHCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
    
/* cqcos */
case EXPR_NODEFUNC_CQCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcosh */
case EXPR_NODEFUNC_CQCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chsec */
case EXPR_NODEFUNC_CHSEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chsec(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chsech */
case EXPR_NODEFUNC_CHSECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chsech(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqsec */
case EXPR_NODEFUNC_CQSEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqsec(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqsech */
case EXPR_NODEFUNC_CQSECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqsech(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* chcsc */
case EXPR_NODEFUNC_CHCSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcsc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcsch */
case EXPR_NODEFUNC_CHCSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcsch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcsc */
case EXPR_NODEFUNC_CQCSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcsc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcsch */
case EXPR_NODEFUNC_CQCSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcsch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chtan */
case EXPR_NODEFUNC_CHTAN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chtan(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chtanh */
case EXPR_NODEFUNC_CHTANH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chtanh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqtan */
case EXPR_NODEFUNC_CQTAN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqtan(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqtanh */
case EXPR_NODEFUNC_CQTANH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqtanh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcot */
case EXPR_NODEFUNC_CHCOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcot(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcoth */
case EXPR_NODEFUNC_CHCOTH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcoth(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcot */
case EXPR_NODEFUNC_CQCOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcot(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcoth */
case EXPR_NODEFUNC_CQCOTH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcoth(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cpxvsin */
case EXPR_NODEFUNC_CPXVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cpxvsin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cpxvsinh */
case EXPR_NODEFUNC_CPXVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cpxvsinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccvsin */
case EXPR_NODEFUNC_CCVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccvsin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccvsinh */
case EXPR_NODEFUNC_CCVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccvsinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cpxvcos */
case EXPR_NODEFUNC_CPXVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cpxvcos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cpxvcosh */
case EXPR_NODEFUNC_CPXVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cpxvcosh(d1+d2*I);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccvcos */
case EXPR_NODEFUNC_CCVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccvcos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccvcosh */
case EXPR_NODEFUNC_CCVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccvcosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chvsin */
case EXPR_NODEFUNC_CHVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chvsin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chvsinh */
case EXPR_NODEFUNC_CHVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chvsinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcvsin */
case EXPR_NODEFUNC_CHCVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcvsin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcvsinh */
case EXPR_NODEFUNC_CHCVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcvsinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqvsin */
case EXPR_NODEFUNC_CQVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqvsin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqvsinh */
case EXPR_NODEFUNC_CQVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqvsinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcvsin */
case EXPR_NODEFUNC_CQCVSIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcvsin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcvsinh */
case EXPR_NODEFUNC_CQCVSINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcvsinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chvcos */
case EXPR_NODEFUNC_CHVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chvcos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chvcosh */
case EXPR_NODEFUNC_CHVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chvcosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcvcos */
case EXPR_NODEFUNC_CHCVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcvcos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcvcosh */
case EXPR_NODEFUNC_CHCVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcvcosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqvcos */
case EXPR_NODEFUNC_CQVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqvcos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqvcosh */
case EXPR_NODEFUNC_CQVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqvcosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcvcos */
case EXPR_NODEFUNC_CQCVCOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcvcos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcvcosh */
case EXPR_NODEFUNC_CQCVCOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcvcosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cesec */
case EXPR_NODEFUNC_CESEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cesec(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cesech */
case EXPR_NODEFUNC_CESECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cesech(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cecsc */
case EXPR_NODEFUNC_CECSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cecsc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cecsch */
case EXPR_NODEFUNC_CECSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cecsch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chesec */
case EXPR_NODEFUNC_CHESEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chesec(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chesech */
case EXPR_NODEFUNC_CHESECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chesech(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* checsc */
case EXPR_NODEFUNC_CHECSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = checsc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* checsch */
case EXPR_NODEFUNC_CHECSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = checsch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqesec */
case EXPR_NODEFUNC_CQESEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqesec(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqesech */
case EXPR_NODEFUNC_CQESECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqesech(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqecsc */
case EXPR_NODEFUNC_CQECSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqecsc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqecsch */
case EXPR_NODEFUNC_CQECSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqecsch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* csinc */
case EXPR_NODEFUNC_CSINC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csinc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* csinch */
case EXPR_NODEFUNC_CSINCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csinch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chsinc */
case EXPR_NODEFUNC_CHSINC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chsinc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chsinch */
case EXPR_NODEFUNC_CHSINCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chsinch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqsinc */
case EXPR_NODEFUNC_CQSINC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqsinc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqsinch */
case EXPR_NODEFUNC_CQSINCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqsinch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccosc */
case EXPR_NODEFUNC_CCOSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccosc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccosch */
case EXPR_NODEFUNC_CCOSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccosch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcosc */
case EXPR_NODEFUNC_CHCOSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcosc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* chcosch */
case EXPR_NODEFUNC_CHCOSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcosch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcosc */
case EXPR_NODEFUNC_CQCOSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcosc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* cqcosch */
case EXPR_NODEFUNC_CQCOSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcosch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* csecc */
case EXPR_NODEFUNC_CSECC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csecc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* csecch */
case EXPR_NODEFUNC_CSECCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csecch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chsecc */
case EXPR_NODEFUNC_CHSECC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chsecc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* chsecch */
case EXPR_NODEFUNC_CHSECCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chsecch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqsecc */
case EXPR_NODEFUNC_CQSECC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqsecc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* cqsecch */
case EXPR_NODEFUNC_CQSECCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqsecch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccscc */
case EXPR_NODEFUNC_CCSCC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccscc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* ccscch */
case EXPR_NODEFUNC_CCSCCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccscch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcscc */
case EXPR_NODEFUNC_CHCSCC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcscc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* chcscch */
case EXPR_NODEFUNC_CHCSCCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcscch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcscc */
case EXPR_NODEFUNC_CQCSCC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcscc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* cqcscch */
case EXPR_NODEFUNC_CQCSCCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcscch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ctanc */
case EXPR_NODEFUNC_CTANC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ctanc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* ctanch */
case EXPR_NODEFUNC_CTANCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ctanch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chtanc */
case EXPR_NODEFUNC_CHTANC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chtanc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* chtanch */
case EXPR_NODEFUNC_CHTANCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chtanch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqtanc */
case EXPR_NODEFUNC_CQTANC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqtanc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* cqtanch */
case EXPR_NODEFUNC_CQTANCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqtanch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccotc */
case EXPR_NODEFUNC_CCOTC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccotc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* ccotch */
case EXPR_NODEFUNC_CCOTCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccotch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* chcotc */
case EXPR_NODEFUNC_CHCOTC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcotc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* chcotch */
case EXPR_NODEFUNC_CHCOTCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = chcotch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cqcotc */
case EXPR_NODEFUNC_CQCOTC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcotc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* cqcotch */
case EXPR_NODEFUNC_CQCOTCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cqcotch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* casin */
case EXPR_NODEFUNC_CASIN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = casin(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* casinh */
case EXPR_NODEFUNC_CASINH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = casinh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cacos */
case EXPR_NODEFUNC_CACOS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cacos(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* cacosh */
case EXPR_NODEFUNC_CACOSH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cacosh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* catan */
case EXPR_NODEFUNC_CATAN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = catan(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* catanh */
case EXPR_NODEFUNC_CATANH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = catanh(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cacsc */
case EXPR_NODEFUNC_CACSC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cacsc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* cacsch */
case EXPR_NODEFUNC_CACSCH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cacsch(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* casec */
case EXPR_NODEFUNC_CASEC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = casec(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* casech */
case EXPR_NODEFUNC_CASECH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = casech(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cacot */
case EXPR_NODEFUNC_CACOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        
    if(isSett(BOOLS_DEGREESENTERING))
	{
		d1 = rad(d1);
		d2 = rad(d2);
	}

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cacot(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
	}
	
/* cacoth */
case EXPR_NODEFUNC_CACOTH:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cacoth(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* det */
case EXPR_NODEFUNC_MATRIXDET:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    const register dim_typ d1ex = (dim_typ)d1;
    const register dim_typ dim = d1ex*d1ex;

    if((!err) && nodes->data.function.nodecount >= dim+1)
        {
        ityp *matrix = NULL;
        
        if(!matrixAlloc(&matrix, (dim_typ2){d1ex,d1ex}))
            return err;

        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                matrix[pos-1] = d2;
            else
                return err;
            }

        bool flag = false;

        *val = det(matrix, d1ex, &flag);
        matrixFree(&matrix);
        }
    else
        return err;

    break;
    }


/* norm */
case EXPR_NODEFUNC_MATRIXNORM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(err)
        return err;

    err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    const register dim_typ d1ex = (dim_typ)d2;
    const register dim_typ dim = d1ex*d1ex;

    if((!err) && nodes->data.function.nodecount >= dim+1)
        {
        ityp d3;
        ityp *matrix = NULL;
        
        if(!matrixAlloc(&matrix, (dim_typ2){d1ex,d1ex}))
            return err;

        for(pos = 1; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d3);
            if(!err)
                matrix[pos-MAX_DIMENSIONS] = d3;
            else
                return err;
            }

        *val = norm(matrix, d1ex, d1);
        matrixFree(&matrix);
        }
    else
        return err;

    break;
    }

/* trace */
case EXPR_NODEFUNC_MATRIXTRACE:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    const register dim_typ d1ex = (dim_typ)d1;
    const register dim_typ dim = d1ex*d1ex;

    if((!err) && nodes->data.function.nodecount >= dim+1)
        {
        	
        ityp *matrix = NULL;
        
        if(!matrixAlloc(&matrix, (dim_typ2){d1ex,d1ex}))
            return err;

        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                matrix[pos-1] = d2;
            else
                return err;
            }

        *val = _matrixTrace(matrix, d1ex);
        matrixFree(&matrix);
        }
    else
        return err;

    break;
    }

/* rank */
case EXPR_NODEFUNC_MATRIXRANK:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(err)
        return err;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d2);

    const dim_typ dex[MAX_DIMENSIONS] =
    {
        (dim_typ)d1,
        (dim_typ)d2
    };

    const register dim_typ dim = dex[ROWS]*dex[COLUMNS];

    if((!err) && nodes->data.function.nodecount >= dim+2)
        {
        ityp d3;
        ityp *matrix = NULL;
        
        if(!matrixAlloc(&matrix, dex))
            return err;

        for(pos = 1; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d3);
            if(!err)
                matrix[pos-MAX_DIMENSIONS] = d3;
            else
                return err;
            }

        *val = rank(matrix, dex);
        matrixFree(&matrix);
        }
    else
        return err;

    break;
    }

/* illchk */
case EXPR_NODEFUNC_MATRIXILLCHK:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    const register dim_typ d1ex = (dim_typ)d1;
    const register dim_typ dim = d1ex*d1ex;
    
    ityp *matrix = NULL;
    
    if(!matrixAlloc(&matrix, (dim_typ2){d1ex, d1ex}))
    	return err;
            
    if((!err) && nodes->data.function.nodecount >= dim+1)
        {

        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                matrix[pos-1] = d2;
            else
                return err;
            }

        const register ityp norms1 = norms(matrix, d1ex);
        const dim_typ d1ex_per2 = d1ex<<1;

        matrix = realloc(matrix, sizeof(ityp)*d1ex_per2*d1ex_per2);
        errMem(matrix, err);
        
        if(!invertMatrix(matrix, d1ex))
            return err;

        *val = norms1*norms(matrix, d1ex);
        matrixFree(&matrix);
        }
    else
        return err;

    break;
    }

/* sprod */
case EXPR_NODEFUNC_SCALARPROD:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    const register dim_typ d1ex = (dim_typ)d1;
    const register dim_typ dim = d1ex<<1;

    if((!err) && nodes->data.function.nodecount >= dim+1)
        {
        ityp vector[d1ex];

        *val = 0.00;

        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                if(pos <= d1ex)
                    vector[pos-1] = d2;
                else
                    *val += vector[pos-d1ex-1];
            else
                return err;
            }
        }
    else
        return err;

    break;
    }


/* log */
case EXPR_NODEFUNC_LOG:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log10(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* log2 */
case EXPR_NODEFUNC_LOG2:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log2(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* pow10 */
case EXPR_NODEFUNC_POW10:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = pow(10.0, d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* ln */
case EXPR_NODEFUNC_LN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* exp */
case EXPR_NODEFUNC_EXP:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = exp(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* expc */
case EXPR_NODEFUNC_EXPC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = expc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* exp10 */
case EXPR_NODEFUNC_EXP10:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = exp10(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* exp2 */
case EXPR_NODEFUNC_EXP2:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = exp2(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* exp2 */
case EXPR_NODEFUNC_EXP2C:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = exp2c(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

case EXPR_NODEFUNC_LOGN:
    {
    EXPRTYPE l1, l2;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        l1 = log(d1);
        EXPR_CHECK_ERR();
        l2 = log(d2);
        EXPR_CHECK_ERR();


        if(l2 == 0.0)
            {
#if(EXPR_ERROR_LEVEL >= EXPR_ERROR_LEVEL_CHECK)
            return EXPR_ERROR_OUTOFRANGE;
#else
            *val = 0.0;
            return EXPR_ERROR_NOERROR;
#endif
            }

        *val = l1 / l2;
        }
    else
        return err;

    break;
    }

/* logc */
case EXPR_NODEFUNC_LOGC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log10c(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* lnc */
case EXPR_NODEFUNC_LNC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = logc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* log2c */
case EXPR_NODEFUNC_LOG2C:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log2c(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* log1p */
case EXPR_NODEFUNC_LOG1P:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log1p(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* log1pc */
case EXPR_NODEFUNC_LOG1PC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log1pc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* log101p */
case EXPR_NODEFUNC_LOG101P:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log101p(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* log101pc */
case EXPR_NODEFUNC_LOG101PC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log101pc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* log21p */
case EXPR_NODEFUNC_LOG21P:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log21p(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* log21pc */
case EXPR_NODEFUNC_LOG21PC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = log21pc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cexp */
case EXPR_NODEFUNC_CEXP:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cexp(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cexpc */
case EXPR_NODEFUNC_CEXPC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cexpc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cexp10 */
case EXPR_NODEFUNC_CEXP10:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cexp10(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cexp10c */
case EXPR_NODEFUNC_CEXP10C:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cexp10c(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cexp2 */
case EXPR_NODEFUNC_CEXP2:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cexp2(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cexp2c */
case EXPR_NODEFUNC_CEXP2C:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = cexp2c(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cpow */
case EXPR_NODEFUNC_CPOW:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
    	{
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        if(!err)
	    	{
	    	ityp d3, d4;
	        err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS+REAL_PART, &d3);
        	if(!err)
		    	{
		        err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS+IMAG_PART, &d4);
			    if(!err)
			        {
			        EXPR_RESET_ERR();
					const register double complex result = cpow(d1+d2*I, d3+d4*I);
					*val = creal(result);
					if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
					{
						*(nodes->data.function.refs[REAL_PART]) = *val;
						*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
					}
					else
						*(nodes->data.function.refs[REAL_PART]) = cimag(result);
					
			        EXPR_CHECK_ERR();
			        }
			    else
			    	return err;
				}
			else
				return err;
			}
		else
			return err;
		}
    else
        return err;

    break;
    }
    
/* croot */
case EXPR_NODEFUNC_CROOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
    	{
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        if(!err)
	    	{
	    	ityp d3, d4;
	        err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS+REAL_PART, &d3);
        	if(!err)
		    	{
		        err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS+IMAG_PART, &d4);
			    if(!err)
			        {
			        EXPR_RESET_ERR();
					const register double complex result = crootnX(d1+d2*I, d3+d4*I);
					*val = creal(result);
					if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
					{
						*(nodes->data.function.refs[REAL_PART]) = *val;
						*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
					}
					else
						*(nodes->data.function.refs[REAL_PART]) = cimag(result);
					
			        EXPR_CHECK_ERR();
			        }
			    else
			    	return err;
				}
			else
				return err;
			}
		else
			return err;
		}
    else
        return err;

    break;
    }
    
/* csqrt */
case EXPR_NODEFUNC_CSQRT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = csqrt(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* ccbrt  */
case EXPR_NODEFUNC_CCBRT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = ccbrt(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clogn */
case EXPR_NODEFUNC_CLOGN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
    	{
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        if(!err)
	    	{
	    	ityp d3, d4;
	        err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS+REAL_PART, &d3);
        	if(!err)
		    	{
		        err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS+IMAG_PART, &d4);
			    if(!err)
			        {
			        EXPR_RESET_ERR();
					const register double complex result = clogbN(d1+d2*I, d3+d4*I);
					*val = creal(result);
					if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
					{
						*(nodes->data.function.refs[REAL_PART]) = *val;
						*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
					}
					else
						*(nodes->data.function.refs[REAL_PART]) = cimag(result);
					
			        EXPR_CHECK_ERR();
			        }
			    else
			    	return err;
				}
			else
				return err;
			}
		else
			return err;
		}
    else
        return err;
        
    break;
    }
    
/* cln */
case EXPR_NODEFUNC_CLN:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clnc */
case EXPR_NODEFUNC_CLNC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clogc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog */
case EXPR_NODEFUNC_CLOG:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog10(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clogc */
case EXPR_NODEFUNC_CLOGC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog10c(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog2 */
case EXPR_NODEFUNC_CLOG2:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog2(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog2c */
case EXPR_NODEFUNC_CLOG2C:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog2c(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog1p */
case EXPR_NODEFUNC_CLOG1P:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog1p(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog1pc */
case EXPR_NODEFUNC_CLOG1PC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog1pc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog101p */
case EXPR_NODEFUNC_CLOG101P:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog101p(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog101pc */
case EXPR_NODEFUNC_CLOG101PC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog101pc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog21p */
case EXPR_NODEFUNC_CLOG21P:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog21p(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* clog21pc */
case EXPR_NODEFUNC_CLOG21PC:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        const register double complex result = clog21pc(d1+d2*I);
        *val = creal(result);
        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
		{
			*(nodes->data.function.refs[REAL_PART]) = *val;
			*(nodes->data.function.refs[IMAG_PART]) = cimag(result);
		}
		else
			*(nodes->data.function.refs[REAL_PART]) = cimag(result);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* carg */
case EXPR_NODEFUNC_CARG:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = isSett(BOOLS_DEGREESENTERING) ? deg(carg(d1+d2*I)) : carg(d1+d2*I);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* cabs */
case EXPR_NODEFUNC_CABS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = cabs(d1+d2*I);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
    
/* qabs */
case EXPR_NODEFUNC_QABS:
    {
    ityp cpx[MAX_QUATERNIONS_UNITS];
    err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_REALPART, cpx);

    if(!err)
    	{
        err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_IPART, &cpx[QUATERNIONS_IPART]);
        if(!err)
	    	{
	    	ityp d3, d4;
	        err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_JPART, &cpx[QUATERNIONS_JPART]);
        	if(!err)
		    	{
		        err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_KPART, &cpx[QUATERNIONS_KPART]);
			    if(!err)
			        {
			        EXPR_RESET_ERR();
					*val = _qabs(cpx);
			        EXPR_CHECK_ERR();
			        }
			    else
			    	return err;
				}
			else
				return err;
			}
		else
			return err;
		}
    else
        return err;

    break;
    }
    
/* oabs */
case EXPR_NODEFUNC_OABS:
    {
    ityp cpx[MAX_OCTONIONS_UNITS];
    err = exprEvalNode(obj, nodes->data.function.nodes, OCTONIONS_REALPART, cpx);

    if(!err)
    	{
        err = exprEvalNode(obj, nodes->data.function.nodes, OCTONIONS_E1PART, &cpx[OCTONIONS_E1PART]);
        if(!err)
	    	{
	    	ityp d3, d4;
	        err = exprEvalNode(obj, nodes->data.function.nodes, OCTONIONS_E2PART, &cpx[OCTONIONS_E2PART]);
        	if(!err)
		    	{
		        err = exprEvalNode(obj, nodes->data.function.nodes, OCTONIONS_E3PART, &cpx[OCTONIONS_E3PART]);
		        if(!err)
			    	{
			        err = exprEvalNode(obj, nodes->data.function.nodes, OCTONIONS_E4PART, &cpx[OCTONIONS_E4PART]);
			        if(!err)
				    	{
				        err = exprEvalNode(obj, nodes->data.function.nodes, OCTONIONS_E5PART, &cpx[OCTONIONS_E5PART]);
				        if(!err)
					    	{
					        err = exprEvalNode(obj, nodes->data.function.nodes, OCTONIONS_E6PART, &cpx[OCTONIONS_E6PART]);
					        if(!err)
						    	{
						        err = exprEvalNode(obj, nodes->data.function.nodes, OCTONIONS_E7PART, &cpx[OCTONIONS_E7PART]);
							    if(!err)
							        {
							        EXPR_RESET_ERR();
									*val = _oabs(cpx);
							        EXPR_CHECK_ERR();
							        }
							    else
							    	return err;
								}	
							else
								return err;
							}
						else
							return err;
						}
					else
						return err;
					}
			    else
			    	return err;
				}
			else
				return err;
			}
		else
			return err;
		}
    else
        return err;

    break;
    }
    
/* sabs */
case EXPR_NODEFUNC_SABS:
    {
    ityp cpx[MAX_SEDENIONS_UNITS];
    err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_REALPART, cpx);

    if(!err)
    	{
        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E1PART, &cpx[SEDENIONS_E1PART]);
        if(!err)
	    	{
	    	ityp d3, d4;
	        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E2PART, &cpx[SEDENIONS_E2PART]);
        	if(!err)
		    	{
		        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E3PART, &cpx[SEDENIONS_E3PART]);
		        if(!err)
			    	{
			        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E4PART, &cpx[SEDENIONS_E4PART]);
			        if(!err)
				    	{
				        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E5PART, &cpx[SEDENIONS_E5PART]);
				        if(!err)
					    	{
					        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E6PART, &cpx[SEDENIONS_E6PART]);
					        if(!err)
						    	{
						        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E7PART, &cpx[SEDENIONS_E7PART]);
						        if(!err)
							    	{
							        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E8PART, &cpx[SEDENIONS_E8PART]);
							        if(!err)
								    	{
								        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E9PART, &cpx[SEDENIONS_E9PART]);
								        if(!err)
									    	{
									        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E10PART, &cpx[SEDENIONS_E10PART]);
									        if(!err)
										    	{
										        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E11PART, &cpx[SEDENIONS_E11PART]);
										        if(!err)
											    	{
											        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E12PART, &cpx[SEDENIONS_E12PART]);
											        if(!err)
												    	{
												        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E13PART, &cpx[SEDENIONS_E13PART]);
												        if(!err)
													    	{
													        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E14PART, &cpx[SEDENIONS_E14PART]);
													        if(!err)
														    	{
														        err = exprEvalNode(obj, nodes->data.function.nodes, SEDENIONS_E15PART, &cpx[SEDENIONS_E15PART]);
															    if(!err)
															        {
															        EXPR_RESET_ERR();
																	*val = _sabs(cpx);
															        EXPR_CHECK_ERR();
															        }
															    else
															    	return err;
																}
															else
																return err;
															}
														else
															return err;
														}
													else
														return err;
													}	
												else
													return err;
												}
											else
												return err;
											}
										else
											return err;
										}
									else
										return err;
									}
							    else
							    	return err;
								}
							else
								return err;
							}
						else
							return err;
						}
					else
						return err;
					}
			    else
			    	return err;
				}
			else
				return err;
			}
		else
			return err;
		}
    else
        return err;

    break;
    } 

/* ceil */
case EXPR_NODEFUNC_CEIL:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        *val = ceil(d1);
        }
    else
        return err;

    break;
    }

/* floor */
case EXPR_NODEFUNC_FLOOR:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        *val = floor(d1);
        }
    else
        return err;

    break;
    }

/* sgeqsolver */
case EXPR_NODEFUNC_SGEQSOLVER:
    {
    EXPRTYPE coeff_c;

    err = exprEvalNode(obj, nodes->data.function.nodes, COEFF_A, &d1);

    if(!err)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, COEFF_B, &d2);
        if(!err)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, COEFF_C, &coeff_c);
            if(!err)
                {


                ityp *abc = NULL;

                if(!matrixAlloc(&abc, (dim_typ2){1, MAX_ABSTRACT_DIMENSIONS}))
                    return err;

                *(abc) = d1;
                *(abc + COEFF_B) = d2;
               	*(abc + COEFF_C) = coeff_c;
               	
               	ityp root[MAX_DIMENSIONS];
               	
                if(!_secondGradeEquationSolver(abc, root))
                   return err;

                *val = *(nodes->data.function.refs[ROOT_X1]) = root[ROOT_X1];
                *(nodes->data.function.refs[ROOT_X2]) = root[ROOT_X2];
                }
                else
                    return err;
            }
            else
                return err;
        }

    break;
    }

/* csum */
case EXPR_NODEFUNC_COMPLEXADD:
    {
    EXPRTYPE secondNum[MAX_COMPLEX_UNITS];

    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        if(!err)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS, secondNum);
            if(!err)
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS+IMAG_PART, &secondNum[IMAG_PART]);

                if(!err)
                    {
                    ityp *cpx = NULL;

                    if(!matrixAlloc(&cpx, (dim_typ2){MAX_COMPLEX_UNITS, MAX_COMPLEX_UNITS}))
                        return err;

                    *(cpx) = d1;
                    *(cpx + IMAG_PART) = d2;
                    *(cpx + MAX_COMPLEX_UNITS) = secondNum[REAL_PART];
                    *(cpx + MAX_COMPLEX_UNITS + IMAG_PART) = secondNum[IMAG_PART];
                    
                    ityp complexRes[MAX_COMPLEX_UNITS];


                    _complexAdd(cpx, complexRes);
                    matrixFree(&cpx);
                    
                    *val = complexRes[REAL_PART];
			        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
					{
						*(nodes->data.function.refs[REAL_PART]) = *val;
						*(nodes->data.function.refs[IMAG_PART]) = complexRes[IMAG_PART];
					}
					else
						*(nodes->data.function.refs[REAL_PART]) = complexRes[IMAG_PART];
						
                    }
                    else
                        return err;
                }
                else
                    return err;
            }
            else
                return err;
        }
        else
        	return err;

    break;
    }

/* cprod */
case EXPR_NODEFUNC_COMPLEXMUL:
    {
    EXPRTYPE secondNum[MAX_COMPLEX_UNITS];

    err = exprEvalNode(obj, nodes->data.function.nodes, REAL_PART, &d1);

    if(!err)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, IMAG_PART, &d2);
        if(!err)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS, secondNum);
            if(!err)
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, MAX_COMPLEX_UNITS+IMAG_PART, &secondNum[IMAG_PART]);

                if(!err)
                    {
                    ityp *cpx = NULL;

                    if(!matrixAlloc(&cpx, (dim_typ2){MAX_COMPLEX_UNITS, MAX_COMPLEX_UNITS}))
                        return err;

                    *(cpx) = d1;
                    *(cpx + IMAG_PART) = d2;
                    *(cpx + MAX_COMPLEX_UNITS) = secondNum[REAL_PART];
                    *(cpx + MAX_COMPLEX_UNITS + IMAG_PART) = secondNum[IMAG_PART];
                    
                    ityp complexRes[MAX_COMPLEX_UNITS];

                    _complexMul(cpx, complexRes);
                    matrixFree(&cpx);
                    
                    *val = complexRes[REAL_PART];
			        if(nodes->data.function.refcount == MAX_COMPLEX_UNITS)
					{
						*(nodes->data.function.refs[REAL_PART]) = *val;
						*(nodes->data.function.refs[IMAG_PART]) = complexRes[IMAG_PART];
					}
					else
						*(nodes->data.function.refs[REAL_PART]) = complexRes[IMAG_PART];
                    }
                    else
                        return err;
                }
                else
                    return err;
            }
            else
                return err;
        }
        else
        	return err;

    break;
    }
    
/* qsum */
case EXPR_NODEFUNC_QUATERNIONSADD:
    {
    ityp d3, d4;
    EXPRTYPE secondNum[MAX_QUATERNIONS_UNITS];

    err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_REALPART, &d1);

    if(!err)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_IPART, &d2);
        if(!err)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_JPART, &d3);
            if(!err)
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_KPART, &d4);
				if(!err)
                	{
               		err = exprEvalNode(obj, nodes->data.function.nodes, MAX_QUATERNIONS_UNITS, secondNum);
               		if(!err)
                		{
                		err = exprEvalNode(obj, nodes->data.function.nodes, MAX_QUATERNIONS_UNITS+QUATERNIONS_IPART, &secondNum[QUATERNIONS_IPART]);
                		if(!err)
               				{
                			err = exprEvalNode(obj, nodes->data.function.nodes, MAX_QUATERNIONS_UNITS+QUATERNIONS_JPART, &secondNum[QUATERNIONS_JPART]);
                			if(!err)
               					{
                				err = exprEvalNode(obj, nodes->data.function.nodes, MAX_QUATERNIONS_UNITS+QUATERNIONS_KPART, &secondNum[QUATERNIONS_KPART]);
				                if(!err)
				                    {
				                    ityp *quaternions = NULL;
				
				                    if(!matrixAlloc(&quaternions, (dim_typ2){MAX_QUATERNIONS_UNITS, MAX_QUATERNIONS_UNITS}))
				                        return err;
				
				                    *(quaternions) = d1;
				                    *(quaternions + QUATERNIONS_IPART) = d2;
				                    *(quaternions + QUATERNIONS_JPART) = d3;
				                    *(quaternions + QUATERNIONS_KPART) = d4;
				                    *(quaternions + MAX_QUATERNIONS_UNITS) = secondNum[QUATERNIONS_REALPART];
				                    *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART) = secondNum[QUATERNIONS_IPART];
				                    *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART) = secondNum[QUATERNIONS_JPART];
				                    *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART) = secondNum[QUATERNIONS_KPART];
				                    
				                    ityp quaternionsRes[MAX_QUATERNIONS_UNITS];
				
				
				                    _quaternionsAdd(quaternions, quaternionsRes);
				                    matrixFree(&quaternions);
				
									*val = *(nodes->data.function.refs[REAL_PART]) = quaternionsRes[QUATERNIONS_REALPART];
				                    *(nodes->data.function.refs[QUATERNIONS_IPART]) = quaternionsRes[QUATERNIONS_IPART];
				                    *(nodes->data.function.refs[QUATERNIONS_JPART]) = quaternionsRes[QUATERNIONS_JPART];
				                    *(nodes->data.function.refs[QUATERNIONS_KPART]) = quaternionsRes[QUATERNIONS_KPART];

				                    }
				                    else
				                        return err;
				                }
				                    else
				                        return err;
				       		}
				                else
				                	return err;
				        }
				        	else
				                return err;
				    }
				    	else
				        	return err;
            	}
                else
                    return err;
            }
            else
                return err;
        }
		else
			return err;

    break;
    }
    
/* qprod */
case EXPR_NODEFUNC_QUATERNIONSMUL:
    {
    ityp d3, d4;
    EXPRTYPE secondNum[MAX_QUATERNIONS_UNITS];

    err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_REALPART, &d1);

    if(!err)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_IPART, &d2);
        if(!err)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_JPART, &d3);
            if(!err)
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, QUATERNIONS_KPART, &d4);
				if(!err)
                	{
               		err = exprEvalNode(obj, nodes->data.function.nodes, MAX_QUATERNIONS_UNITS, secondNum);
               		if(!err)
                		{
                		err = exprEvalNode(obj, nodes->data.function.nodes, MAX_QUATERNIONS_UNITS+QUATERNIONS_IPART, &secondNum[QUATERNIONS_IPART]);
                		if(!err)
               				{
                			err = exprEvalNode(obj, nodes->data.function.nodes, MAX_QUATERNIONS_UNITS+QUATERNIONS_JPART, &secondNum[QUATERNIONS_JPART]);
                			if(!err)
               					{
                				err = exprEvalNode(obj, nodes->data.function.nodes, MAX_QUATERNIONS_UNITS+QUATERNIONS_KPART, &secondNum[QUATERNIONS_KPART]);
				                if(!err)
				                    {
				                    ityp *quaternions = NULL;
				
				                    if(!matrixAlloc(&quaternions, (dim_typ2){MAX_QUATERNIONS_UNITS, MAX_QUATERNIONS_UNITS}))
				                        return err;
				
				                    *(quaternions) = d1;
				                    *(quaternions + QUATERNIONS_IPART) = d2;
				                    *(quaternions + QUATERNIONS_JPART) = d3;
				                    *(quaternions + QUATERNIONS_KPART) = d4;
				                    *(quaternions + MAX_QUATERNIONS_UNITS) = secondNum[QUATERNIONS_REALPART];
				                    *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART) = secondNum[QUATERNIONS_IPART];
				                    *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART) = secondNum[QUATERNIONS_JPART];
				                    *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART) = secondNum[QUATERNIONS_KPART];
				                    
				                    ityp quaternionsRes[MAX_QUATERNIONS_UNITS];
				
				
				                    _quaternionsAdd(quaternions, quaternionsRes);
				                    matrixFree(&quaternions);
				
									*val = *(nodes->data.function.refs[REAL_PART]) = quaternionsRes[QUATERNIONS_REALPART];
				                    *(nodes->data.function.refs[QUATERNIONS_IPART]) = quaternionsRes[QUATERNIONS_IPART];
				                    *(nodes->data.function.refs[QUATERNIONS_JPART]) = quaternionsRes[QUATERNIONS_JPART];
				                    *(nodes->data.function.refs[QUATERNIONS_KPART]) = quaternionsRes[QUATERNIONS_KPART];

				                    }
				                    else
				                        return err;
				                }
				                    else
				                        return err;
				       		}
				                else
				                	return err;
				        }
				        	else
				                return err;
				    }
				    	else
				        	return err;
            	}
                else
                    return err;
            }
            else
                return err;
        }
		else
			return err;

    break;
    }

/* rand */
case EXPR_NODEFUNC_RAND:
    {
    long a;

    /* Perform random routine directly */
    a = ((long)(*(nodes->data.function.refs[0]))) * 214013L + 2531011L;
    *(nodes->data.function.refs[0]) = (EXPRTYPE)a;

    *val =  (EXPRTYPE)((a >> 16) & 0x7FFF) / (EXPRTYPE)(32768);
    break;
    }

/* random */
case EXPR_NODEFUNC_RANDOM:
    {
    EXPRTYPE diff, rval;
    long a;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        diff = d2 - d1;

        /* Perform random routine directly */
        a = ((long)(*(nodes->data.function.refs[0]))) * 214013L + 2531011L;
        *(nodes->data.function.refs[0]) = (EXPRTYPE)a;

        rval = (EXPRTYPE)((a >> 16) & 0x7FFF) / (EXPRTYPE)(32767);

        *val = (rval * diff) + d1;
        }
    else
        return err;

    break;
    }

/* randomize */
case EXPR_NODEFUNC_RANDOMIZE:
    {
    static int curcall = 0;

    ++curcall;

    *(nodes->data.function.refs[0]) = (EXPRTYPE)((clock() + 1024 + curcall) * time(NULL));

    break;
    }

/* deg */
case EXPR_NODEFUNC_DEG:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        *val = deg(d1);
        }
    else
        return err;

    break;
    }

/* rad */
case EXPR_NODEFUNC_RAD:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        *val = rad(d1);
        }
    else
        return err;

    break;
    }

/* recttopolr */
case EXPR_NODEFUNC_RECTTOPOLR:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sqrt((d1 * d1) + (d2 * d2));
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* recttopola */
case EXPR_NODEFUNC_RECTTOPOLA:
    {
    EXPRTYPE tmp;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        tmp = atan2(d2, d1);
        EXPR_CHECK_ERR();

        if(tmp < 0.0)
            *val = tmp = (2.0 * M_PI);
        else
            *val = tmp;
        }
    else
        return err;

    break;
    }

/* poltorectx */
case EXPR_NODEFUNC_POLTORECTX:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = d1 * cos(d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* poltorecty */
case EXPR_NODEFUNC_POLTORECTY:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = d1 * sin(d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* cbase */
case EXPR_NODEFUNC_CBASE:
    {
    EXPRTYPE v;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &v);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 2, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = changeBase((int64_t)v, (sel_typ)d1, (sel_typ)d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* npnum */
case EXPR_NODEFUNC_NPNUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = N_prime_Number(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* primr */
case EXPR_NODEFUNC_PRIMORIAL:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = primr(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* fpnsum */
case EXPR_NODEFUNC_FPNSUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = fpnsum(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* fibnc */
case EXPR_NODEFUNC_FIBONACCIAL:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = fibnc(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* mcm */
case EXPR_NODEFUNC_LCM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = math_mcm((uint64_t)d1, (uint64_t)d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* MCD */
case EXPR_NODEFUNC_GCD:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = math_MCD((uint64_t)d1, (uint64_t)d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* gsum */
case EXPR_NODEFUNC_GSUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = gsum(d1, (int64_t)d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* fact */
case EXPR_NODEFUNC_FACT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = fact(d1); // setting ERROR condition
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sfact */
case EXPR_NODEFUNC_SFACT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sfact(d1); // setting ERROR condition
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* strlng */
case EXPR_NODEFUNC_STIRLING:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = stirling(d1); // setting ERROR condition
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* fibo */
case EXPR_NODEFUNC_FIBO:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = fibo(d1); // setting ERROR condition
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* perm */
case EXPR_NODEFUNC_PERMS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = perm(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* permrep */
case EXPR_NODEFUNC_PERMSREP:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        const register dim_typ dim = nodes->data.function.nodecount-1;
        ityp args[dim];
        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos-1] = d2;
            else
                return err;
            }
        // Overheaded but it is necessary if Domain Check is enabled.
        if(dcheck && summation(dim, SUMMATION_SUM, args) != d1)
        	return(err = EXPR_ERROR_SYNTAX);
        *val = perm_rep(d1, dim, args);
        }
    else
        return err;

    break;
    }
    
/* kperm */
case EXPR_NODEFUNC_KPERMS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = kperm(d1, d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* kpermrep */
case EXPR_NODEFUNC_KPERMSREP:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = kperm_rep(d1, d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* comb */
case EXPR_NODEFUNC_COMBS:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = comb(d1, d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }
    
/* combrep */
case EXPR_NODEFUNC_COMBSREP:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = comb_rep(d1, d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }


/* asum */
case EXPR_NODEFUNC_ASUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = asum(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* gasum */
case EXPR_NODEFUNC_GASUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = gasum((uint64_t)d1, (uint64_t)d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* fsum */
case EXPR_NODEFUNC_FSUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = fsum(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* fasum */
case EXPR_NODEFUNC_FASUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = fasum(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sfasum */
case EXPR_NODEFUNC_SFASUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = sfasum(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* fnnsum */
case EXPR_NODEFUNC_FNNSUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = fnnsum(d1);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* sum */
case EXPR_NODEFUNC_SUM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = summation(nodes->data.function.nodecount, INVERSE_OPS, args);
        }
    else
        return err;

    break;
    }

/* product */
case EXPR_NODEFUNC_PRODUCT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = productory(nodes->data.function.nodecount, INVERSE_OPS, args);
        }
    else
        return err;

    break;
    }

/* media */
case EXPR_NODEFUNC_MEDIA:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_media(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }
    
/* var */
case EXPR_NODEFUNC_VARIANCE:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_variance(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }
    
/* cov */
case EXPR_NODEFUNC_COVARIANCE:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        if(nodes->data.function.nodecount%2)
        	return(err = EXPR_ERROR_SYNTAX);
        const register dim_typ dim = nodes->data.function.nodecount<<1;
        ityp args[MAX_DIMENSIONS][dim];
        args[FIRST_VECTOR][0] = d1;
        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[FIRST_VECTOR][pos] = d2;
            else
                return err;
            }
        args[SECOND_VECTOR][0] = d1;
        for( ; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[SECOND_VECTOR][pos-dim] = d2;
            else
                return err;
            }
        *val = math_covariance(dim, args[FIRST_VECTOR], args[SECOND_VECTOR]);
        }
    else
        return err;

    break;
    }
    
/* stddev */
case EXPR_NODEFUNC_STDDEV:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_stddev(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }
    
/* outlier */
case EXPR_NODEFUNC_OUTLIER:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        const register dim_typ dim = nodes->data.function.nodecount-1;
        ityp args[dim];
        const ityp outlier_idx = d1;
        if(dim <= outlier_idx)
        	return(err = EXPR_ERROR_SYNTAX);
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos-1] = d2;
            else
                return err;
            }
        *val = math_outlier(dim, outlier_idx , args);
        }
    else
        return err;

    break;
    }
    
/* outlier2 */
case EXPR_NODEFUNC_OUTLIER2:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        const register dim_typ dim = nodes->data.function.nodecount-MAX_DIMENSIONS;
        ityp args[dim];
        const ityp outlier_idx = d1;
        if(dim <= outlier_idx)
        	return(err = EXPR_ERROR_SYNTAX);
        if(!(err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2)))
    		{

	       	const ityp outlier_constant = d2;
	        ityp d3 = 0.00;
	        for(pos = 1; ++pos < nodes->data.function.nodecount; )
	            {
	            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d3);
	            if(!err)
	                args[pos-MAX_DIMENSIONS] = d3;
	            else
	                return err;
	            }
	        *val = math_outlier2(dim, outlier_idx, outlier_constant, args);
	    	}
	    else
	    	return err;
        }
    else
        return err;

    break;
    }
    
/* map */
case EXPR_NODEFUNC_MAP:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);
	dim_typ funcID;
    if(!err)
        {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(err || (!err && d2 != (funcID = (dim_typ)d2) || funcID < FID_SIN || funcID > LAST_FID))
                return(err = EXPR_ERROR_SYNTAX);
     	   *val = ext_math.functions[funcID](d1);
        }
    else
        return err;

    break;
    }


/* geomedia */
case EXPR_NODEFUNC_GEOMEDIA:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_geomedia(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }

/* armedia */
case EXPR_NODEFUNC_ARMEDIA:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_armedia(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }

/* powmedia */
case EXPR_NODEFUNC_POWMEDIA:
    {
    EXPRTYPE v;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &v);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 1; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos-1] = d2;
            else
                return err;
            }
        *val = math_powmedia(nodes->data.function.nodecount, v, args);
        }
    else
        return err;

    break;
    }

/* cval */
case EXPR_NODEFUNC_CVAL:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_scale(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }
    
/* q1 */
case EXPR_NODEFUNC_FIRSTQUARTILE:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_first_quartile(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }

/* mediana */
case EXPR_NODEFUNC_MEDIANA:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_mediana(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }
    
/* q3 */
case EXPR_NODEFUNC_THIRDQUARTILE:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
            else
                return err;
            }
        *val = math_third_quartile(nodes->data.function.nodecount, args);
        }
    else
        return err;

    break;
    }

/* if */
case EXPR_NODEFUNC_IF:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        if(d1 != 0.0)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, 1, val);
            if(err)
                return err;
            }
        else
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, 2, val);
            if(err)
                return err;
            }
        }
    else
        return err;

    break;
    }

/* select */
case EXPR_NODEFUNC_SELECT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        if(d1 < 0.0)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, 1, val);
            if(err)
                return err;
            }
        else if(d1 == 0.0)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, 2, val);
            if(err)
                return err;
            }
        else
            {
            if(nodes->data.function.nodecount == 3)
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, 2, val);
                if(err)
                    return err;
                }
            else
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, 3, val);
                if(err)
                    return err;
                }
            }
        }
    else
        return err;

    break;
    }

/* equal */
case EXPR_NODEFUNC_EQUAL:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        *val = (d1 == d2) ? 1.0 : 0.0;
        }
    else
        return err;

    break;
    }

/* above */
case EXPR_NODEFUNC_ABOVE:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        *val = (d1 > d2) ? 1.0 : 0.0;
        }
    else
        return err;

    break;
    }

/* below */
case EXPR_NODEFUNC_BELOW:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        *val = (d1 < d2) ? 1.0 : 0.0;
        }
    else
        return err;

    break;
    }

/* avg */
case EXPR_NODEFUNC_AVG:
    {
    d2 = 0.0;

    for(pos = 0; pos < nodes->data.function.nodecount; ++pos)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d1);
        if(!err)
            {
            d2 += d1;
            }
        else
            return err;
        }

    *val = d2 / (EXPRTYPE)(nodes->data.function.nodecount);

    break;
    }

/* clip */
case EXPR_NODEFUNC_CLIP:
    {
    EXPRTYPE v;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &v);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 2, &d2);

    if(!err)
        {
        if(v < d1)
            *val = d1;
        else if(v > d2)
            *val = d2;
        else
            *val = v;
        }
    else
        return err;

    break;
    }

/* clamp */
case EXPR_NODEFUNC_CLAMP:
    {
    EXPRTYPE v, tmp;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &v);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 2, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        tmp = fmod(v - d1, d2 - d1);
        EXPR_CHECK_ERR();

        if(tmp < 0.0)
            *val = tmp * d2;
        else
            *val = tmp + d1;
        }
    else
        return err;

    break;
    }

/* pntchange */
case EXPR_NODEFUNC_PNTCHANGE:
    {
    EXPRTYPE n1, n2, pnt;
    EXPRTYPE odiff, ndiff, perc;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 2, &n1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 3, &n2);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 4, &pnt);

    if(!err)
        {
        odiff = d2 - d1;
        ndiff = n2 - n1;

        if(odiff == 0.0)
            {
            *val = d1;
            return EXPR_ERROR_NOERROR;
            }

        perc = (pnt - d1) / odiff;

        *val = n1 + (perc * ndiff);
        }
    else
        return err;

    break;
    }

/* poly */
case EXPR_NODEFUNC_POLY:
    {
    EXPRTYPE total, curpow;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        curpow = (EXPRTYPE)(nodes->data.function.nodecount) - 2.0;
        total = 0.0;

        for(pos = 0; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(err)
                return err;

            EXPR_RESET_ERR();
            total = total + (d2 * pow(d1, curpow));
            EXPR_CHECK_ERR();

            curpow = curpow - 1.0;
            }
        }
    else
        return err;

    *val = total;
    break;
    }

/* and */
case EXPR_NODEFUNC_AND:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        if(d1 == 0.0 || d2 == 0.0)
            *val = 0.0;
        else
            *val = 1.0;
        }
    else
        return err;

    break;
    }

/* or */
case EXPR_NODEFUNC_OR:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        if(d1 != 0.0 || d2 != 0.0)
            *val = 1.0;
        else
            *val = 0.0;
        }
    else
        return err;

    break;
    }

/* not */
case EXPR_NODEFUNC_NOT:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        if(d1 != 0.0)
            *val = 0.0;
        else
            *val = 1.0;
        }
    else
        return err;

    break;
    }

/* for */
case EXPR_NODEFUNC_FOR:
    {
    int pos;
    EXPRTYPE test;

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &test);

    if(!err)
        {
        while(test != 0.00)
            {
            for(pos = 2; ++pos < nodes->data.function.nodecount; )
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, pos, val);
                if(err)
                    return err;
                }

            err = exprEvalNode(obj, nodes->data.function.nodes, 2, &d1);
            if(err)
                return err;

            err = exprEvalNode(obj, nodes->data.function.nodes, 1, &test);
            if(err)
                return err;
            }
        }
    else
        return err;

    break;
    }

/* many */
case EXPR_NODEFUNC_MANY:
    {
    for(pos = 0; pos < nodes->data.function.nodecount; ++pos)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, pos, val);
        if(err)
            return err;
        }

    break;
    }

