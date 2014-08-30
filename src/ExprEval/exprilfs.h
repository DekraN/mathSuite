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

    break;
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
        *val = (-1+2*(d1 >= 0))*d1;
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
            if((*val) < MIN_MEMOIZABLE_INDEX+1 || (*val) > MIN_STIRLINGREQUIRES_NUMBER)
                return err;
            access(curLayout)->min_stirlingrequires_number = *val;
            }
        }
        else
            (*val) = access(curLayout)->min_stirlingrequires_number;
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
            _changeAlgebraDims(*val);
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
            access(curLayout)->max_memoizable_indeces[FUNCTION_FIBONACCI] = *val;
            }
        }
        else
            (*val) = access(curLayout)->max_memoizable_indeces[FUNCTION_FIBONACCI];
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
            access(curLayout)->max_memoizable_indeces[FUNCTION_FATTORIALE] = *val;
            }
        }
        else
            (*val) = access(curLayout)->max_memoizable_indeces[FUNCTION_FATTORIALE];
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
            access(curLayout)->max_memoizable_indeces[FUNCTION_EVEN_SFATTORIALE] = *val;
            }
        }
        else
            (*val) = access(curLayout)->max_memoizable_indeces[FUNCTION_EVEN_SFATTORIALE];
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
            access(curLayout)->max_memoizable_indeces[FUNCTION_ODD_SFATTORIALE] = *val;
            }
        }
        else
            (*val) = access(curLayout)->max_memoizable_indeces[FUNCTION_ODD_SFATTORIALE];
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
case EXPR_NODEFUNC_HHCVSIN:
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
case EXPR_NODEFUNC_HHCVSINH:
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

/* det */
case EXPR_NODEFUNC_MATRIXDET:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    const dim_typ d1ex = (dim_typ)d1;
    const dim_typ dim = d1ex*d1ex;

    if((!err) && nodes->data.function.nodecount >= dim+1)
        {
        ityp vector[dim];

        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                vector[pos-1] = d2;
            else
                return err;
            }

        ityp **matrix = NULL;

        if(!matrixAlloc(&matrix, (dim_typ2){d1ex,d1ex}))
            return err;

        vectorToMatrix((dim_typ2){d1ex,d1ex}, vector, matrix);

        bool flag = false;

        *val = det(matrix, d1ex, &flag);
        matrixFree(&matrix, d1ex);
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

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d2);

    const dim_typ d1ex = (dim_typ)d2;
    const dim_typ dim = d1ex*d1ex;

    if((!err) && nodes->data.function.nodecount >= dim+1)
        {
        ityp d3;
        ityp vector[dim];

        for(pos = 1; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d3);
            if(!err)
                vector[pos-2] = d3;
            else
                return err;
            }

        ityp **matrix = NULL;

        if(!matrixAlloc(&matrix, (dim_typ2){d1ex,d1ex}))
            return err;

        vectorToMatrix((dim_typ2){d1ex,d1ex}, vector, matrix);

        *val = norm(matrix, d1ex, d1);
        matrixFree(&matrix, d1ex);
        }
    else
        return err;

    break;
    }

/* trace */
case EXPR_NODEFUNC_MATRIXTRACE:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    const dim_typ d1ex = (dim_typ)d1;
    const dim_typ dim = d1ex*d1ex;

    if((!err) && nodes->data.function.nodecount >= dim+1)
        {
        ityp vector[dim];

        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                vector[pos-1] = d2;
            else
                return err;
            }

        ityp **matrix = NULL;

        if(!matrixAlloc(&matrix, (dim_typ2){d1ex,d1ex}))
            return err;

        vectorToMatrix((dim_typ2){d1ex,d1ex}, vector, matrix);

        *val = trace(matrix, d1ex);
        matrixFree(&matrix, d1ex);
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

    const dim_typ dim = dex[RAWS]*dex[COLUMNS];

    if((!err) && nodes->data.function.nodecount >= dim+2)
        {
        ityp d3;
        ityp vector[dim];

        for(pos = 1; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d3);
            if(!err)
                vector[pos-2] = d3;
            else
                return err;
            }

        ityp **matrix = NULL;

        if(!matrixAlloc(&matrix, dex))
            return err;

        vectorToMatrix(dex, vector, matrix);

        *val = rank(matrix, dex);
        matrixFree(&matrix, dex[RAWS]);
        }
    else
        return err;

    break;
    }

/* illchk */
case EXPR_NODEFUNC_MATRIXILLCHK:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    const dim_typ d1ex = (dim_typ)d1;
    const dim_typ dim = d1ex*d1ex;

    if((!err) && nodes->data.function.nodecount >= dim+1)
        {
        ityp vector[dim];

        for(pos = 0; ++pos < dim; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                vector[pos-1] = d2;
            else
                return err;
            }

        ityp **matrix = NULL;

        if(!matrixAlloc(&matrix, (dim_typ2){d1ex, d1ex}))
            return err;

        vectorToMatrix((dim_typ2){d1ex, d1ex}, vector, matrix);


        const ityp norms1 = norms(matrix, d1ex);
        const dim_typ d1ex_per2 = d1ex<<1;

        matrix = realloc(matrix, sizeof(ityp *)*d1ex_per2);
        errMem(matrix, err);

        dim_typ i;

        for(i=0; i<d1ex; ++i)
        {
            matrix[i] = realloc(matrix[i], sizeof(ityp)*d1ex_per2);
            errMem(matrix[i], err);
        }

        if(!invertMatrix(matrix, d1ex))
            return err;

        *val = norms1*norms(matrix, d1ex);
        matrixFree(&matrix, d1ex);
        }
    else
        return err;

    break;
    }

/* sprod */
case EXPR_NODEFUNC_SCALARPROD:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    const dim_typ d1ex = (dim_typ)d1;
    const dim_typ dim = d1ex<<1;

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

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);
        if(!err)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, 1, &coeff_c);
            if(!err)
                {


                ityp **abc = NULL;

                if(!matrixAlloc(&abc, (dim_typ2){1, MAX_ABSTRACT_DIMENSIONS}))
                    return err;

                abc[0][COEFF_A] = d1;
                abc[0][COEFF_B] = d2;
                abc[0][COEFF_C] = coeff_c;

                // ityp root[MAX_DIMENSIONS];

                // if(!_secondGradeEquationSolver(abc, root))
                   //  return err;

                if(!_secondGradeEquationSolver(abc[0], *(nodes->data.function.refs)))
                   return err;

                *val = *(nodes->data.function.refs[ROOT_X1]);

                // *(nodes->data.function.refs[ROOT_X1]) = root[ROOT_X1];
                // *(nodes->data.function.refs[ROOT_X2]) = root[ROOT_X2];
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
case EXPR_NODEFUNC_COMPLEXSUM:
    {
    EXPRTYPE secondNum[MAX_DIMENSIONS];

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);
        if(!err)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, 1, secondNum);
            if(!err)
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, 1, &secondNum[IMAG_PART]);

                if(!err)
                    {
                    ityp **complex = NULL;

                    if(!matrixAlloc(&complex, (dim_typ2){MAX_DIMENSIONS, MAX_DIMENSIONS}))
                        return err;

                    complex[XRAW][REAL_PART] = d1;
                    complex[XRAW][IMAG_PART] = d2;
                    complex[YRAW][REAL_PART] = secondNum[REAL_PART];
                    complex[YRAW][IMAG_PART] = secondNum[IMAG_PART];

                    // ityp complexRes[MAX_DIMENSIONS];

                    _complexSum(complex, *(nodes->data.function.refs)); // _complexSum(complex, complexRes);
                    matrixFree(&complex, MAX_DIMENSIONS);

                    *val = *(nodes->data.function.refs[IMAG_PART]);

                    // *(nodes->data.function.refs[REAL_PART]) = complexRes[REAL_PART];
                    // *(nodes->data.function.refs[IMAG_PART]) = complexRes[IMAG_PART];
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

    break;
    }

/* cprod */
case EXPR_NODEFUNC_COMPLEXPROD:
    {
    EXPRTYPE secondNum[MAX_DIMENSIONS];

    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);
        if(!err)
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, 1, secondNum);
            if(!err)
                {
                err = exprEvalNode(obj, nodes->data.function.nodes, 1, &secondNum[IMAG_PART]);

                if(!err)
                    {
                    ityp **complex = NULL;

                    if(!matrixAlloc(&complex, (dim_typ2){MAX_DIMENSIONS, MAX_DIMENSIONS}))
                        return err;

                    complex[XRAW][REAL_PART] = d1;
                    complex[XRAW][IMAG_PART] = d2;
                    complex[YRAW][REAL_PART] = secondNum[REAL_PART];
                    complex[YRAW][IMAG_PART] = secondNum[IMAG_PART];

                    // ityp complexRes[MAX_DIMENSIONS];

                    _complexProd(complex, *(nodes->data.function.refs)); // _complexProd(complex, complexRes);
                    matrixFree(&complex, MAX_DIMENSIONS);

                    *val = *(nodes->data.function.refs[IMAG_PART]);
                    // *(nodes->data.function.refs[REAL_PART]) = complexRes[REAL_PART];
                    // *(nodes->data.function.refs[IMAG_PART]) = complexRes[IMAG_PART];
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
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

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
case EXPR_NODEFUNC_PERM:
    {
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

    if(!err)
        {
        EXPR_RESET_ERR();
        *val = perm(d1, d2);
        EXPR_CHECK_ERR();
        }
    else
        return err;

    break;
    }

/* comb */
case EXPR_NODEFUNC_COMB:
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
        {
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d1);
        }

    if(!err)
        {
        ityp args[nodes->data.function.nodecount];
        args[0] = d1;
        for(pos = 1; ++pos < nodes->data.function.nodecount; )
            {
            err = exprEvalNode(obj, nodes->data.function.nodes, pos, &d2);
            if(!err)
                args[pos] = d2;
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

/* day */
case EXPR_NODEFUNC_DAY:
    {
    EXPRTYPE d3;
    err = exprEvalNode(obj, nodes->data.function.nodes, 0, &d1);

    if(!err)
        {
            err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);
            if(!err)
            {
                err = exprEvalNode(obj, nodes->data.function.nodes, 2, &d3);
                *val = getDayNumber(d1, d2, d3);
            }
            else
                return err;
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
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

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
        err = exprEvalNode(obj, nodes->data.function.nodes, 1, &d2);

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

