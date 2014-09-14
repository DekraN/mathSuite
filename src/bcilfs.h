/*

BASECALC support header
This has been created in order
to prevent programs.c getting too crowded.
*/


case BCALC_EXIT:
    safeExit(a);
    break;

case BCALC_ADDIZIONE:
    c = (a + b);
    break;

case BCALC_SOTTRAZIONE:
    c = (a - b);
    break;

case BCALC_MOLTIPLICAZIONE:
    c = (a * b);
    break;

case BCALC_DIVISIONE:
    if(!b)
    {
        if(dcheck)
            printErr(33, "a/0 Operation (Division per ZERO) isn't allowed");
        else
            if(!a) printf2(COLOR_SYSTEM, "\nRESULT: Indeterminate Form 0/0.\n");
            else printf2(COLOR_SYSTEM, "\nRESULT: %cinf.\n", a > 0 ? '+' : '-');
        return;
    }
    c = (a / b);
    break;

case BCALC_RESTO:
    if(dcheck && (!b))
    {
        printErr(33, "You are trying to get Rest from a 'per-zero DIVISION");
        return;
    }
    c = ((int)a % (int)b);
    break;

case BCALC_ADDIZIONEBINARIA:
{
    char *c = NULL;
    if((c = binaryAlgSum(a, b, BAS_SUM)))
        free(c);
    return;
}

case BCALC_SOTTRAZIONEBINARIA:
{
    char *c = NULL;
    if((c = binaryAlgSum(a, b, BAS_SUB)))
       free(c);
    return;
}

case BCALC_COMPLEMENTO:
{
    char *c = NULL;
    if((c = binNumComp(a)))
        free(c);
    return;
}

case BCALC_ELEVAMENTOAPOTENZA:
    if(dcheck && (!a) && (!b))
    {
        printf2(COLOR_SYSTEM, "\nRESULT: Indeterminate Form 0^0.\n");
        return;
    }
    c = pow(a, b);
    break;

case BCALC_EXPANDEXPC:
    if(a)
    {
        if(dcheck && !b)
        {
            printErr(33, "expc function isn't defined in 0.\n(DOMAIN: R\\{x| x=0} )");
            return;
        }
        c = expc(b);
    }
    else
        c = exp(b);
    break;

case BCALC_EXP10ANDEXP10C:
    if(a)
    {
        if(dcheck && !b)
        {
            printErr(33, "exp10c function isn't defined in 0.\n(DOMAIN: R\\{x| x=0} )");
            return;
        }
        c = exp10c(b);
    }
    else
        c = exp10(b);
    break;


case BCALC_EXP2ANDEXP2C:
    if(a)
    {
        if(dcheck && !b)
        {
            printErr(33, "exp2c function isn't defined in 0.\n(DOMAIN: R\\{x| x=0} )");
            return;
        }
        c = exp2c(b);
    }
    else
        c = exp2(b);
    break;

case BCALC_RADICENESIMA:
{
    if(dcheck && a < 0 && !(((int)b)%2))
    {
        printErr(33, "Invalid N Root Function Argument.\n(DOMAIN: [0,+inf] )");
        return;
    }
    c = rootnX(a, b);
    // c = pow(a, (1/b));
    break;
}

case BCALC_RADICEQUADRATA:
    if(dcheck && a < 0)
    {
        printErr(33, "Invalid Square Root Function Argument.\n(DOMAIN: [0,+inf] )");
        return;
    }
    c = sqrt(a);
    break;

case BCALC_RADICECUBICA:
    c = cbrt(a);
    break;

case BCALC_LOGARITMO:

    if(dcheck && b <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlog function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]0,+inf[ )");
        return;
    }

    c = (a ? log10(b) : log(b));
    break;

case BCALC_LOGARITMO2:

    if(dcheck && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlog2 function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]0,+inf[ )");
        return;
    }

    c = log2(a);
    break;

case BCALC_LOGARITMOBN:

    if(dcheck && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlogBN functions isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]0,+inf[ )");
        return;
    }

    c = logbN(a, b);
    break;

case BCALC_LOGARITMOC:

    if(dcheck && b <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlogc function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]0,+inf[ )");
        return;
    }

    c = (a ? log10c(b) : logc(b));
    break;

case BCALC_LOGARITMO2C:

    if(dcheck && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlog2c function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]0,+inf] )");
        return;
    }

    c = log2c(a);
    break;

case BCALC_LOGARITMO1P:

    if(dcheck && a <= -1)
    {
        printf2(COLOR_SYSTEM, "\nlog1p function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]-1,+inf[ )");
        return;
    }

    c = log1p(a);
    break;

case BCALC_LOGARITMO1PC:

    if(dcheck && (a <= -1 || !a))
    {
        printf2(COLOR_SYSTEM, "\nlog1pc function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]-1,+inf]\\0 )");
        return;
    }

    c = log1pc(a);
    break;
    
case BCALC_LOGARITMO101P:

    if(dcheck && a <= -1)
    {
        printf2(COLOR_SYSTEM, "\nlog101p function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]-1,+inf[ )");
        return;
    }

    c = log101p(a);
    break;

case BCALC_LOGARITMO101PC:

    if(dcheck && (a <= -1 || !a))
    {
        printf2(COLOR_SYSTEM, "\nlog101pc function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]-1,+inf]\\0 )");
        return;
    }

    c = log101pc(a);
    break;
    
case BCALC_LOGARITMO21P:

    if(dcheck && a <= -1)
    {
        printf2(COLOR_SYSTEM, "\nlog21p function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]-1,+inf[ )");
        return;
    }

    c = log21p(a);
    break;

case BCALC_LOGARITMO21PC:

    if(dcheck && (a <= -1 || !a))
    {
        printf2(COLOR_SYSTEM, "\nlog21pc function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ".\n");
        printErr(33, "(DOMAIN: ]-1,+inf]\\0 )");
        return;
    }

    c = log21pc(a);
    break;

case BCALC_CEXP:
	viewComplexResult(CRV_DONTDGCHECK, "cexp", cexp(a+b*I));
	return;
case BCALC_CEXPC:
	if(dcheck && (!a) && !b)
        printErr(33, "cexpc function isn't defined in (0,0).\n(DOMAIN: C\\(0,0))");
	else
		viewComplexResult(CRV_DONTDGCHECK, "cexpc", cexpc(a+b*I));
	return;
	
case BCALC_CEXP10:	
	viewComplexResult(CRV_DONTDGCHECK, "cexp10", cexp10(a+b*I));
	return;
	
case BCALC_CEXP10C:
	if(dcheck && (!a) && !b)
        printErr(33, "cexp10c function isn't defined in (0,0).\n(DOMAIN: C\\(0,0))");
	else
		viewComplexResult(CRV_DONTDGCHECK, "cexp10c", cexp10c(a+b*I));
	return;
	
case BCALC_CEXP2:
	viewComplexResult(CRV_DONTDGCHECK, "cexp2", cexp2(a+b*I));
	return;
	
case BCALC_CEXP2C:
	if(dcheck && (!a) && !b)
        printErr(33, "cexp2c function isn't defined in (0,0).\n(DOMAIN: C\\(0,0))");
	else
	viewComplexResult(CRV_DONTDGCHECK, "cexp2c", cexp2c(a+b*I));
	return;
	
case BCALC_CSQRT:
	if(dcheck && (!b) && a < 0)
        printErr(33, "Invalid Complex Square Root Function Argument.\n(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) < 0} )");
	else
		viewComplexResult(CRV_DONTDGCHECK, "cqsrt", csqrt(a+b*I));
	return;
	
case BCALC_CCBRT:
	viewComplexResult(CRV_DONTDGCHECK, "ccbrt", ccbrt(a+b*I));
	return;
	
case BCALC_CLOG:
	if(dcheck && (!b) && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nclog function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= 0} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog", clog(a+b*I));
	return;
	
case BCALC_CLOGC:
	if(dcheck && (!b) && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nclogc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= 0} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clogc", clogc(a+b*I));
	return;	
	
case BCALC_CLOG10:
	if(dcheck && (!b) && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nclog10 function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= 0} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog10", clog10(a+b*I));
	return;
	
case BCALC_CLOG10C:
	if(dcheck && (!b) && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nclog10c function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= 0} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog10c", clog10c(a+b*I));
	return;
	
case BCALC_CLOG2:
	if(dcheck && (!b) && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nclog2 function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= 0} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog2", clog2(a+b*I));
	return;
	
case BCALC_CLOG2C:
	if(dcheck && (!b) && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nclog2c function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= 0} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog2c", clog2c(a+b*I));
	return;
	
case BCALC_CLOG1P:
	if(dcheck && (!b) && a <= -1)
    {
        printf2(COLOR_SYSTEM, "\nclog1p function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= -1} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog1p", clog1p(a+b*I));
	return;
	
case BCALC_CLOG1PC:
	if(dcheck && (!b) && (a <= -1 || !a))
    {
        printf2(COLOR_SYSTEM, "\nclog1pc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) <= -1 V Re(x) = 0)} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog1pc", clog1pc(a+b*I));
	return;
	
case BCALC_CLOG101P:
	if(dcheck && (!b) && a <= -1)
    {
        printf2(COLOR_SYSTEM, "\nclog101p function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= -1} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog101p", clog101p(a+b*I));
	return;
	
case BCALC_CLOG101PC:
	if(dcheck && (!b) && (a <= -1 || !a))
    {
        printf2(COLOR_SYSTEM, "\nclog101pc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) <= -1 V Re(x) = 0)} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog101pc", clog101pc(a+b*I));
	return;
	
case BCALC_CLOG21P:
	if(dcheck && (!b) && a <= -1)
    {
        printf2(COLOR_SYSTEM, "\nclog21p function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) <= -1} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog21p", clog21p(a+b*I));
	return;
	
case BCALC_CLOG21PC:
	if(dcheck && (!b) && (a <= -1 || !a))
    {
        printf2(COLOR_SYSTEM, "\nclog21pc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) <= -1 V Re(x) = 0)} )");
    }
    else
		viewComplexResult(CRV_DONTDGCHECK, "clog21pc", clog21pc(a+b*I));
	return;
	
case BCALC_CARG:
	c = carg(a+b*I);
	break;
	
case BCALC_CABS:
	c = cabs(a+b*I);
	break;
	
case BCALC_BITCOUNTER:
    c = countbits(a);
    break;

case BCALC_UBITCOUNTER:
    c = ucountbits(a);
    break;

case BCALC_VERSION:
    c = strtod(PROG__VERSION, NULL);
    break;
    
case BCALC_EXITCHAR:
	if(a)
	{
		if(a < 0 || a == 'A' || a == 'B')
		{
			printErr(33, "Invalid inserted Exit Char");
			return;
		}
		c = a;
		access(curLayout)->exit_char = a;
	}
	else
		c = access(curLayout)->exit_char;
	break;

case BCALC_PREC:

    if(a)
    {
        if(a < MIN_PRECISION || a > MAX_PRECISION)
        {
            printErr(33, "Invalid inserted Precision Value.\nMust be a non-negative integer between %hu and %hu", MIN_PRECISION, MAX_PRECISION);
            return;
        }
        c = access(curLayout)->precision = a;
    }
    else
        c = access(curLayout)->precision;
    break;


case BCALC_SFACT:
	
    if(a)
    {
        if(a < MIN_STABFACT || a > MAX_STABFACT)
        {
            printErr(33, "Invalid Inserted Stabilizer Factor Value.\nMust be a non-negative integer between %hu and %hu", MIN_STABFACT, MAX_STABFACT);
            return;
        }
        c = access(curLayout)->stabilizer_factor = a;
    }
    else
        c = access(curLayout)->stabilizer_factor;
    break;
    
case BCALC_MINSRNUMBER:
	
	if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MIN_STIRLING_NUMBER)
        {
            printErr(33, "Invalid Inserted Stabilizer Factor Value.\nMust be a non-negative integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MIN_STIRLING_NUMBER);
            return;
        }
        c = access(curLayout)->min_stirling_number = a;
    }
    else
        c = access(curLayout)->min_stirling_number;
    break;

case BCALC_ALGEBRA:
    if(a)
    {
    if(a < MIN_ALGEBRA || a > MAX_ALGEBRA)
        {
            printErr(33, "Invalid Inserted Algebra Identifier.\nMust be a non-negative integer between %hu and %hu", MIN_ALGEBRA, MAX_ALGEBRA);
            return;
        }
        c = a;
        access(curLayout)->algebra = a;
    }
    else
        c = access(curLayout)->algebra;
    break;
    
case BCALC_OUTLIERCONST:
	if(a)
    {
    if(a < MIN_OUTLIER_CONSTANT || a > MAX_OUTLIER_CONSTANT)
        {
            printErr(33, "Invalid Inserted Outlier Constant.\nMust be a non-negative float between %.*f and %.*f", DEFAULT_PRECISION, MIN_OUTLIER_CONSTANT, DEFAULT_PRECISION, MAX_OUTLIER_CONSTANT);
            return;
        }
        c = a;
        access(curLayout)->outlier_constant = a;
    }
    else
        c = access(curLayout)->outlier_constant;
    break;
	

case BCALC_RSEED:

    if(a)
    {
        if(a < MIN_RANDOMSEED || a > MAX_RANDOMSEED)
        {
            printErr(33, "Invalid inserted Random Seed Value.\nMust be a non-negative integer between %hu and %hu", MIN_RANDOMSEED, MAX_RANDOMSEED);
            return;
        }
        c = access(random_seed) = a;
    }
    else
        c = access(random_seed);
    break;

case BCALC_MMIFIBO:

    if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MAX_FIBONACCI_MEMOIZABLE_INDEX)
        {
            printErr(33, "Invalid inserted %s MIM Value.\nMust be a non-negative integer between %hu and %hu", suite_c.memoizers_names[FUNCTION_FIBONACCI], MIN_MEMOIZABLE_INDEX+1, MAX_FIBONACCI_MEMOIZABLE_INDEX);;
            return;
        }
        c = access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = a;
    }
    else
        c = access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI];
    break;

case BCALC_MMIFACT:

    if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MAX_FATTORIALE_MEMOIZABLE_INDEX)
        {
            printErr(33, "Invalid inserted %s MIM Value.\nMust be a non-negative integer between %hu and %hu", suite_c.memoizers_names[FUNCTION_FATTORIALE], MIN_MEMOIZABLE_INDEX+1, MAX_FATTORIALE_MEMOIZABLE_INDEX);
            return;
        }
        c = access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = a;
    }
    else
        c = access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE];
    break;
    
case BCALC_MMIEVENSFACT:

    if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX)
        {
            printErr(33, "Invalid inserted %s MIM Value.\nMust be a non-negative integer between %hu and %hu", suite_c.memoizers_names[FUNCTION_EVEN_SFATTORIALE], MIN_MEMOIZABLE_INDEX+1, MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX);
            return;
        }
        c = access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = a;
    }
    else
        c = access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE];
    break;
    
case BCALC_MMIODDSFACT:

    if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX)
        {
            printErr(33, "Invalid inserted %s MIM Value.\nMust be a non-negative integer between %hu and %hu", suite_c.memoizers_names[FUNCTION_ODD_SFATTORIALE], MIN_MEMOIZABLE_INDEX+1, MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX);
            return;
        }
        c = access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] = a;
    }
    else
        c = access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE];
    break;


case BCALC_TRASFORMAANGOLI:
    c = (a ? deg(b) : rad(b));
    break;

case BCALC_SINANDSINH:
    c = (a ? sinh(b) : sin(b));
    break;

case BCALC_COSANDCOSH:
    c = (a  ? cosh(b) : cos(b));
    break;

case BCALC_TANANDTANH:
    if(a) c = tanh(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ntan function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x| x = PI/2 + k*PI} )");
            return;
        }
        c = tan(b);
    }
    break;

case BCALC_CSCANDCSCH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\ncsch function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: ]-inf,0[U]0,+inf[ )");
            return;
        }
        c = csch(b);
    }
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ncsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = csc(b);
    }
    break;

case BCALC_SECANDSECH:
    if(a) c = sech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nsec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = sec(b);
    }
    break;

case BCALC_COTANDCOTH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\ncoth function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: ]-inf,0[U]0,+inf[ )");
            return;
        }
        c = coth(b);
    }
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ncot function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = cot(b);
    }
    break;

case BCALC_HSINANDHSINH:
    c = (a ? hsinh(b) : hsin(b));
    break;

case BCALC_QSINANDQSINH:
    c = (a ? qsinh(b) : qsin(b));
    break;

case BCALC_HCOSANDHCOSH:
    c = (a ? hcosh(b) : hcos(b));
    break;

case BCALC_QCOSANDQCOSH:
    c = (a ? qcosh(b) : qcos(b));
    break;

case BCALC_HSECANDHSECH:
    if(a) c = hsech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhsec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = hsec(b);
    }
    break;

case BCALC_QSECANDQSECH:
    if(a) c = qsech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqsec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = qsec(b);
    }
    break;

case BCALC_HCSCANDHCSCH:
    if(a) c = hcsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhcsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hcsc(b);
    }
    break;

case BCALC_QCSCANDQCSC:
    if(a) c = qcsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqcsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qcsc(b);
    }
    break;

case BCALC_HTANANDHTANH:
    if(a) c = htanh(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhtan function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x| x = PI/2 + k*PI} )");
            return;
        }
        c = htan(b);
    }
    break;

case BCALC_QTANANDQTANH:
    if(a) c = qtanh(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqtan function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x| x = PI/2 + k*PI} )");
            return;
        }
        c = qtan(b);
    }
    break;

case BCALC_HCOTANDHCOTH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\nhcoth function isn't defined in: 0.\n");
            printErr(33, "(DOMAIN: ]-inf,0[U]0,+inf[ )");
            return;
        }
        c = hcoth(b);
    }
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhcot function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hcot(b);
    }
    break;

case BCALC_QCOTANDQCOTH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\nqcoth function isn't defined in: 0.\n");
            printErr(33, "(DOMAIN: ]-inf,0[U]0,+inf[ )");
            return;
        }
        c = qcoth(b);
    }
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqcot function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qcot(b);
    }
    break;

case BCALC_VSINANDVSINH:
    c = (a ? vsinh(b) : vsin(b));
    break;

case BCALC_CVSINANDCVSINH:
    c = (a ? cvsinh(b) : cvsin(b));
    break;

case BCALC_VCOSANDVCOSH:
    c = (a ? vcosh(b) : vcos(b));
    break;

case BCALC_CVCOSANDCVCOSH:
    c = (a ? cvcosh(b) : cvcos(b));
    break;

case BCALC_HVSINANDHVSINH:
    c = (a ? hvsinh(b) : hvsin(b));
    break;

case BCALC_HCVSINANDHCVSINH:
    c = (a ? hcvsinh(b) : hcvsin(b));
    break;

case BCALC_QVSINANDQVSINH:
    c = (a ? qvsinh(b) : qvsin(b));
    break;

case BCALC_QCVSINANDQCVSINH:
    c = (a ? qcvsinh(b) : qcvsin(b));
    break;

case BCALC_HVCOSANDHVCOSH:
    c = (a ? hvcosh(b) : hvcos(b));
    break;

case BCALC_HCVCOSANDHCVCOSH:
    c = (a ? hcvcosh(b) : hcvcos(b));
    break;

case BCALC_QVCOSANDQVCOSH:
    c = (a ? qvcosh(b) : qvcos(b));
    break;

case BCALC_QCVCOSANDQCVCOSH:
    c = (a ? qcvcosh(b) : qcvcos(b));
    break;

case BCALC_ESECANDESECH:
    if(a) c = esech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nesec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = esec(b);
    }
    break;

case BCALC_ECSCANDECSCH:
    if(a) c = ecsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\necsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = ecsc(b);
    }
    break;

case BCALC_HESECANDHESECH:
    if(a) c = hesech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhesec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = hesec(b);
    }
    break;

case BCALC_HECSCANDHECSCH:
    if(a) c = hecsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhecsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hecsc(b);
    }
    break;

case BCALC_QESECANDQESECH:
    if(a) c = qesech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqesec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = qesec(b);
    }
    break;

case BCALC_QECSCANDQECSCH:
    if(a) c = qecsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqecsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qecsc(b);
    }
    break;

case BCALC_SINCANDSINCH:
    c = (a ? sinch(b) : sinc(b));
    break;

case BCALC_HSINCANDHSINCH:
    c = (a ? hsinch(b) : hsinc(b));
    break;

case BCALC_QSINCANDQSINCH:
    c = (a ? qsinch(b) : qsinc(b));
    break;

case BCALC_COSCANDCOSCH:
    if(a) c = cosch(b);
    else
    {
        if(dcheck && !b)
        {
            printf2(COLOR_SYSTEM, "\ncosc function isn't defined in: 0.\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0} )");
            return;
        }
        c = cosc(b);
    }
    break;

case BCALC_HCOSCANDHCOSCH:
    if(a) c = hcosch(b);
    else
    {
        if(dcheck && !b)
        {
            printf2(COLOR_SYSTEM, "\nhcosc function isn't defined in: 0.\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0} )");
            return;
        }
        c = hcosc(b);
    }
    break;

case BCALC_QCOSCANDQCOSCH:
    if(a) c = qcosch(b);
    else
    {
        if(dcheck && !b)
        {
            printf2(COLOR_SYSTEM, "\nqcosc function isn't defined in: 0.\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0} )");
            return;
        }
        c = qcosc(b);
    }
    break;

case BCALC_SECCANDSECCH:
    if(a) c = secch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nsecc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = secc(b);
    }
    break;

case BCALC_HSECCANDHSECCH:
    if(a) c = hsecch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhsecc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = hsecc(b);
    }
    break;

case BCALC_QSECCANDQSECCH:
    if(a) c = qsecch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqsecc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = qsecc(b);
    }
    break;

case BCALC_CSCCANDCSCCH:
    if(a) c = cscch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ncscc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = cscc(b);
    }
    break;

case BCALC_HCSCCANDHCSCCH:
    if(a) c = hcscch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhcscc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hcscc(b);
    }
    break;

case BCALC_QCSCCANDQCSCCH:
    if(a) c = qcscch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqcscc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qcscc(b);
    }
    break;

case BCALC_TANCANDTANCH:
    if(a) c = tanch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ntanc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = tanc(b);
    }
    break;

case BCALC_HTANCANDHTANCH:
    if(a) c = htanch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhtanc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = htanc(b);
    }
    break;

case BCALC_QTANCANDQTANCH:
    if(a) c = qtanch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqtanc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = qtanc(b);
    }
    break;

case BCALC_COTCANDCOTCH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\ncotch function isn't defined in: 0.\n");
            printErr(33, "(DOMAIN: ]-inf,0[U]0,+inf[ )");
            return;
        }
        c = cotch(b);
    }
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ncotc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = cotc(b);
    }
    break;

case BCALC_HCOTCANDHCOTCH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\nhcotch function isn't defined in: 0.\n");
            printErr(33, "(DOMAIN: ]-inf,0[U]0,+inf[ )");
            return;
        }
        c = hcotch(b);
    }
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhcotc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hcotc(b);
    }
    break;

case BCALC_QCOTCANDQCOTCH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\nqcotch function isn't defined in: 0.\n");
            printErr(33, "(DOMAIN: ]-inf,0[U]0,+inf[ )");
            return;
        }
        c = qcotch(b);
    }
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqcotc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qcotc(b);
    }
    break;


case BCALC_ASINANDASINH:
    if(a == 0)
    {
        if(dcheck && !(TRIGONOMETRIC_DOMAIN(b)))
        {
            printf2(COLOR_SYSTEM, "\nasin function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: [-1,1] )");
            return;
        }
        c = asin(b);
    }
    else c = asinh(b);
    break;

case BCALC_ACOSANDACOSH:
    if(a)
    {

        if(dcheck && b < 1)
        {
            printf2(COLOR_SYSTEM, "\nacosh function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: [1,+inf] )");
            return;
        }
        c = acosh(b);
    }
    else
    {
        if(dcheck && !(TRIGONOMETRIC_DOMAIN(b)))
        {
            printf2(COLOR_SYSTEM, "\nacos function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: [-1,1] )");
            return;
        }
        c = acos(b);
    }
    break;

case BCALC_ATANANDATANH:
    if(a)
    {
        if(dcheck && !(TRIGONOMETRIC_DOMAIN(b)))
        {
            printf2(COLOR_SYSTEM, "\natanh function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: [-1,1] )");
            return;
        }
        c = atanh(b);
    }
    else c = atan(b);
    break;

case BCALC_ATAN2:
    if(dcheck && (!a) && (!b))
    {
        printf2(COLOR_SYSTEM, "\natan2 function isn't defined in: (0,0). ");
        printErr(33, "(DOMAIN: (R^2)\\{P=(y,x)|P=(0,0)} )");
        return;
    }
    c = atan2(a, b);
    break;

case BCALC_ACSCANDACSCH:
    if(a) c = acsch(b);
    else
    {
        if(dcheck && b > -1 && b < 1)
        {
            printf2(COLOR_SYSTEM, "\nacsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: ]-inf,-1]U[1,+inf[ )");
            return;
        }
        c = acsc(b);
    }
    break;

case BCALC_ASECANDASECH:
    if(a)
    {
        if(dcheck && (b < 0 || b >= 1))
        {
            printf2(COLOR_SYSTEM, "\nasech function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: [0, 1) )");
            return;
        }
        c = asech(b);

    }
    else
    {
        if(dcheck && b > -1 && b < 1)
        {
            printf2(COLOR_SYSTEM, "\nasec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: ]-inf,-1]U[1,+inf[ )");
            return;
        }
        c = asec(b);
    }
    break;

case BCALC_ACOTANDACOTH:
    if(a)
    {
        if(dcheck && TRIGONOMETRIC_DOMAIN(b))
        {
            printf2(COLOR_SYSTEM, "\nacoth function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ".\n");
            printErr(33, "(DOMAIN: ]-inf,-1[U]1,+inf[ )");
            return;
        }
        c = acoth(b);
    }
    else c = acot(b);
    break;
    
case BCALC_CSIN:
	viewComplexResult(CRV_DODGCHECK, "csin", csin(a+b*I));
	return;
case BCALC_CSINH:
	viewComplexResult(CRV_DODGCHECK, "csinh", csinh(a+b*I));
	return;
case BCALC_CCOS:	
	viewComplexResult(CRV_DODGCHECK, "ccos", ccos(a+b*I));
	return;
case BCALC_CCOSH:
	viewComplexResult(CRV_DODGCHECK, "ccosh", ccosh(a+b*I));
	return;	
case BCALC_CTAN:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nctan function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ctan", ctan(a+b*I));
	return;
case BCALC_CTANH:
	viewComplexResult(CRV_DODGCHECK, "ctanh", ctanh(a+b*I));
	return;
case BCALC_CCSC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nccsc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ccsc", ccsc(a+b*I));
	return;
case BCALC_CCSCH:
	if(dcheck && (!a) && (!b))
    {
        printf2(COLOR_SYSTEM, "\ncsch function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ccsch", ccsch(a+b*I));
	return;
case BCALC_CSEC:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncsec function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x|Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "csec", csec(a+b*I));
	return;
case BCALC_CSECH:
	viewComplexResult(CRV_DODGCHECK, "csech", csech(a+b*I));
	return;
case BCALC_CCOT:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nccot function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
        return;
    }
	viewComplexResult(CRV_DODGCHECK, "ccot", ccot(a+b*I));
	return;
case BCALC_CCOTH:
	if(dcheck && (!a) && (!b))
    {
        printf2(COLOR_SYSTEM, "\nccoth function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ccoth", ccoth(a+b*I));
	return;
case BCALC_CHSIN:
	viewComplexResult(CRV_DODGCHECK, "chsin", chsin(a+b*I));
	return;
case BCALC_CHSINH:
	viewComplexResult(CRV_DODGCHECK, "chsinh", chsinh(a+b*I));
	return;
case BCALC_CQSIN:
	viewComplexResult(CRV_DODGCHECK, "cqsin", cqsin(a+b*I));
	return;
case BCALC_CQSINH:
	viewComplexResult(CRV_DODGCHECK, "cqsinh", cqsinh(a+b*I));
	return;
case BCALC_CHCOS:
	viewComplexResult(CRV_DODGCHECK, "chcos", chcos(a+b*I));
	return;
case BCALC_CHCOSH:
	viewComplexResult(CRV_DODGCHECK, "chcosh", chcosh(a+b*I));
	return;
case BCALC_CQCOS:
	viewComplexResult(CRV_DODGCHECK, "cqcos", cqcos(a+b*I));
	return;	
case BCALC_CQCOSH:
	viewComplexResult(CRV_DODGCHECK, "cqcosh", cqcosh(a+b*I));
	return;
case BCALC_CHSEC:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchsec function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
	}
	else
		viewComplexResult(CRV_DODGCHECK, "chsec", chsec(a+b*I));
	return;
case BCALC_CHSECH:
	viewComplexResult(CRV_DODGCHECK, "chsech", chsech(a+b*I));
	return;
case BCALC_CQSEC:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqsec function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqsec", cqsec(a+b*I));
	return;
case BCALC_CQSECH:
	viewComplexResult(CRV_DODGCHECK, "cqsech", cqsech(a+b*I));
	return;
case BCALC_CHCSC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchcsc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chcsc", chcsc(a+b*I));
	return;
case BCALC_CHCSCH:
	viewComplexResult(CRV_DODGCHECK, "chcsch", chcsch(a+b*I));
	return;
case BCALC_CQCSC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqcsc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqcsc", cqcsc(a+b*I));
	return;
case BCALC_CQCSCH:
	viewComplexResult(CRV_DODGCHECK, "cqcsch", cqcsch(a+b*I));
	return;
case BCALC_CHTAN:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchtan function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chtan", chtan(a+b*I));
	return;
case BCALC_CHTANH:
	viewComplexResult(CRV_DODGCHECK, "chtanh", chtanh(a+b*I));
	return;
case BCALC_CQTAN:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqtan function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqtan", cqtan(a+b*I));
	return;
case BCALC_CQTANH:
	viewComplexResult(CRV_DODGCHECK, "cqtanh", cqtanh(a+b*I));
	return;
case BCALC_CHCOT:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchcot function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chcot", chcot(a+b*I));
	return;
case BCALC_CHCOTH:
	if(dcheck && (!a) && (!b))
    {
        printf2(COLOR_SYSTEM, "\nchcoth function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chcoth", chcoth(a+b*I));
	return;
case BCALC_CQCOT:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqcot function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqcot", cqcot(a+b*I));
	return;
case BCALC_CQCOTH:
	if(dcheck && (!a) && (!b))
    {
        printf2(COLOR_SYSTEM, "\ncqcoth function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqcoth", cqcoth(a+b*I));
	return;
case BCALC_CVSIN:
	viewComplexResult(CRV_DODGCHECK, "cpxvsin", cpxvsin(a+b*I));
	return;
case BCALC_CVSINH:
	viewComplexResult(CRV_DODGCHECK, "cpxvsinh", cpxvsinh(a+b*I));
	return;
case BCALC_CCVSIN:
	viewComplexResult(CRV_DODGCHECK, "ccvsin", ccvsin(a+b*I));
	return;
case BCALC_CCVSINH:
	viewComplexResult(CRV_DODGCHECK, "ccvsinh", ccvsinh(a+b*I));
	return;
case BCALC_CVCOS:
	viewComplexResult(CRV_DODGCHECK, "cpxvcos", cpxvcos(a+b*I));
	return;
case BCALC_CVCOSH:
	viewComplexResult(CRV_DODGCHECK, "cpxvcosh", cpxvcosh(a+b*I));
	return;
case BCALC_CCVCOS:
	viewComplexResult(CRV_DODGCHECK, "ccvcos", ccvcos(a+b*I));
	return;
case BCALC_CCVCOSH:
	viewComplexResult(CRV_DODGCHECK, "ccvcosh", ccvcosh(a+b*I));
	return;
case BCALC_CHVSIN:
	viewComplexResult(CRV_DODGCHECK, "chvsin", chvsin(a+b*I));
	return;
case BCALC_CHVSINH:
	viewComplexResult(CRV_DODGCHECK, "chvsinh", chvsinh(a+b*I));
	return;
case BCALC_CHCVSIN:
	viewComplexResult(CRV_DODGCHECK, "chcvsin", chcvsin(a+b*I));
	return;
case BCALC_CHCVSINH:
	viewComplexResult(CRV_DODGCHECK, "chcvsinh", chcvsinh(a+b*I));
	return;
case BCALC_CQVSIN:	
	viewComplexResult(CRV_DODGCHECK, "cqvsin", cqvsin(a+b*I));
	return;
case BCALC_CQVSINH:
	viewComplexResult(CRV_DODGCHECK, "cqvsinh", cqvsinh(a+b*I));
	return;
case BCALC_CQCVSIN:
	viewComplexResult(CRV_DODGCHECK, "cqcvsin", cqcvsin(a+b*I));
	return;
case BCALC_CQCVSINH:
	viewComplexResult(CRV_DODGCHECK, "cqcvsinh", cqcvsinh(a+b*I));
	return;
case BCALC_CHVCOS:	
	viewComplexResult(CRV_DODGCHECK, "chcvcos", chcvcos(a+b*I));
	return;
case BCALC_CHVCOSH:
	viewComplexResult(CRV_DODGCHECK, "chcvcosh", chcvcosh(a+b*I));
	return;
case BCALC_CHCVCOS:
	viewComplexResult(CRV_DODGCHECK, "chcvcos", chcvcos(a+b*I));
	return;
case BCALC_CHCVCOSH:
	viewComplexResult(CRV_DODGCHECK, "chcvcosh", chcvcosh(a+b*I));
	return;
case BCALC_CQVCOS:
	viewComplexResult(CRV_DODGCHECK, "cqvcos", cqvcos(a+b*I));
	return;
case BCALC_CQVCOSH:
	viewComplexResult(CRV_DODGCHECK, "cqvcosh", cqvcosh(a+b*I));
	return;
case BCALC_CQCVCOS:
	viewComplexResult(CRV_DODGCHECK, "cqcvcos", cqcvcos(a+b*I));
	return;
case BCALC_CQCVCOSH:
	viewComplexResult(CRV_DODGCHECK, "cqcvcosh", cqcvcosh(a+b*I));
	return;
case BCALC_CESEC:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncesec function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cesec", cesec(a+b*I));
	return;
case BCALC_CESECH:
	viewComplexResult(CRV_DODGCHECK, "cesech", cesech(a+b*I));
	return;
case BCALC_CECSC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncecsc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cecsc", cecsc(a+b*I));
	return;
case BCALC_CECSCH:
	viewComplexResult(CRV_DODGCHECK, "cecsch", cecsch(a+b*I));
	return;
case BCALC_CHESEC:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchesec function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chesec", chesec(a+b*I));
	return;
case BCALC_CHESECH:	
	viewComplexResult(CRV_DODGCHECK, "chesech", chesech(a+b*I));
	return;
case BCALC_CHECSC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchecsc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: R\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
	
	}
	else
		viewComplexResult(CRV_DODGCHECK, "checsc", checsc(a+b*I));
	return;
case BCALC_CHECSCH:
	viewComplexResult(CRV_DODGCHECK, "checsch", checsch(a+b*I));
	return;
case BCALC_CQESEC:
	if(dcheck && (!b) && (!trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqesec function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = PI/2 + k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqesec", cqesec(a+b*I));
	return;
case BCALC_CQESECH:
	viewComplexResult(CRV_DODGCHECK, "cqesech", cqesech(a+b*I));
	return;
case BCALC_CQECSC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqecsc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqecsc", cqecsc(a+b*I));
	return;
case BCALC_CQECSCH:
	viewComplexResult(CRV_DODGCHECK, "cqecsch", cqecsch(a+b*I));
	return;
case BCALC_CSINC:
	viewComplexResult(CRV_DODGCHECK, "csinc", csinc(a+b*I));
	return;
case BCALC_CSINCH:
	viewComplexResult(CRV_DODGCHECK, "csinch", csinch(a+b*I));
	return;
case BCALC_CHSINC:
	viewComplexResult(CRV_DODGCHECK, "chsinc", chsinc(a+b*I));
	return;
case BCALC_CHSINCH:
	viewComplexResult(CRV_DODGCHECK, "chsinch", chsinch(a+b*I));
	return;
case BCALC_CQSINC:
	viewComplexResult(CRV_DODGCHECK, "cqsinc", cqsinc(a+b*I));
	return;
case BCALC_CQSINCH:
	viewComplexResult(CRV_DODGCHECK, "cqsinch", cqsinch(a+b*I));
	return;
case BCALC_CCOSC:
	if(dcheck && (!a) && !b)
    {
        printf2(COLOR_SYSTEM, "\nccosc function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ccosc", ccosc(a+b*I));
	return;
case BCALC_CCOSCH:
	viewComplexResult(CRV_DODGCHECK, "ccosch", ccosch(a+b*I));
	return;
case BCALC_CHCOSC:	
	if(dcheck && (!a) && !b)
    {
        printf2(COLOR_SYSTEM, "\nchcosc function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chcosc", chcosc(a+b*I));
	return;
case BCALC_CHCOSCH:
	viewComplexResult(CRV_DODGCHECK, "chcosch", chcosch(a+b*I));
	return;
case BCALC_CQCOSC:
	if(dcheck && (!a) && !b)
    {
        printf2(COLOR_SYSTEM, "\ncqcosc function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
	else	
		viewComplexResult(CRV_DODGCHECK, "cqcosc", cqcosc(a+b*I));
	return;
case BCALC_CQCOSCH:
	viewComplexResult(CRV_DODGCHECK, "cqcosch", cqcosch(a+b*I));
	return;
case BCALC_CSECC:	
	if(dcheck && (!b) && ((!a) || !trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncsecc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) = 0 V Re(x) = PI/2 + k*PI)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "csecc", csecc(a+b*I));
	return;
case BCALC_CSECCH:
	viewComplexResult(CRV_DODGCHECK, "csecch", csecch(a+b*I));
	return;
case BCALC_CHSECC:
	if(dcheck && (!b) && ((!a) || !trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchsecc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) = 0 V Re(x) = PI/2 + k*PI)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chsecc", chsecc(a+b*I));
	return;
case BCALC_CHSECCH:
	viewComplexResult(CRV_DODGCHECK, "chsecch", chsecch(a+b*I));
	return;
case BCALC_CQSECC:	
	if(dcheck && (!b) && ((!a) || !trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqsecc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) = 0 V Re(x) = PI/2 + k*PI)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqsecc", cqsecc(a+b*I));
	return;
case BCALC_CQSECCH:
	viewComplexResult(CRV_DODGCHECK, "cqsecch", cqsecch(a+b*I));
	return;
case BCALC_CCSCC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nccscc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ccscc", ccscc(a+b*I));
	return;
case BCALC_CCSCCH:
	viewComplexResult(CRV_DODGCHECK, "ccscc", ccscch(a+b*I));
	return;
case BCALC_CHCSCC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchcscc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else	
		viewComplexResult(CRV_DODGCHECK, "chcscc", chcscc(a+b*I));
	return;
case BCALC_CHCSCCH:
	viewComplexResult(CRV_DODGCHECK, "chcscch", chcscch(a+b*I));
	return;
case BCALC_CQCSCC:	
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqcscc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqcscc", cqcscc(a+b*I));
	return;
case BCALC_CQCSCCH:
	viewComplexResult(CRV_DODGCHECK, "cqcscch", cqcscch(a+b*I));
	return;
case BCALC_CTANC:
	if(dcheck && (!b) && ((!a) || !trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nctanc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) = 0 V Re(x) = PI/2 + k*PI)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ctanc", ctanc(a+b*I));
	return;
case BCALC_CTANCH:
	viewComplexResult(CRV_DODGCHECK, "ctanch", ctanch(a+b*I));
	return;
case BCALC_CHTANC:
	if(dcheck && (!b) && ((!a) || !trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchtanc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) = 0 V Re(x) = PI/2 + k*PI)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chtanc", chtanc(a+b*I));
	return;
case BCALC_CHTANCH:
	viewComplexResult(CRV_DODGCHECK, "chtanch", chtanch(a+b*I));
	return;
case BCALC_CQTANC:
	if(dcheck && (!b) && ((!a) || !trigonometric_domain(a, M_PI_2, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqtanc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) = 0 V Re(x) = PI/2 + k*PI)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqtanc", cqtanc(a+b*I));
	return;
case BCALC_CQTANCH:
	viewComplexResult(CRV_DODGCHECK, "cqtanch", cqtanch(a+b*I));
	return;
case BCALC_CCOTC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nccotc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ccotc", ccotc(a+b*I));
	return;	
case BCALC_CCOTCH:
	if(dcheck && (!a) && !b)
    {
        printf2(COLOR_SYSTEM, "\nccotch function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "ccotch", ccotch(a+b*I));
	return;
case BCALC_CHCOTC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\nchcotc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chcotc", chcotc(a+b*I));
	return;
case BCALC_CHCOTCH:
	if(dcheck && (!a) && !b)
    {
        printf2(COLOR_SYSTEM, "\nchcotch function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "chcotch", chcotch(a+b*I));
	return;
case BCALC_CQCOTC:
	if(dcheck && (!b) && (!trigonometric_domain(a, 0, M_PI)))
    {
        printf2(COLOR_SYSTEM, "\ncqcotc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) = k*PI} )");
    }
    else	
		viewComplexResult(CRV_DODGCHECK, "cqcotc", cqcotc(a+b*I));
	return;
case BCALC_CQCOTCH:
	if(dcheck && (!a) && !b)
    {
        printf2(COLOR_SYSTEM, "\ncqcotch function isn't defined in: (0,0).\n");
        printErr(33, "(DOMAIN: C\\(0,0))");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cqcotch", cqcotch(a+b*I));
	return;
case BCALC_CASIN:
	if(dcheck && (!b) && !(TRIGONOMETRIC_DOMAIN(a)))
    {
        printf2(COLOR_SYSTEM, "\ncasin function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) < -1 V Re(x) > 1)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "casin", casin(a+b*I));
	return;
case BCALC_CASINH:
	viewComplexResult(CRV_DODGCHECK, "casinh", casinh(a+b*I));
	return;
case BCALC_CACOS:
	if(dcheck && (!b) && !(TRIGONOMETRIC_DOMAIN(a)))
    {
        printf2(COLOR_SYSTEM, "\ncacos function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
    	printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) < -1 V Re(x) > 1)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cacos", cacos(a+b*I));
	return;
case BCALC_CACOSH:
	if(dcheck && (!b) && a < 1)
    {
        printf2(COLOR_SYSTEM, "\ncacosh function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) < 1} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cacosh", cacosh(a+b*I));
	return;
case BCALC_CATAN:
	viewComplexResult(CRV_DODGCHECK, "catan", catan(a+b*I));
	return;	
case BCALC_CATANH:
	if(dcheck && (!b) && !(TRIGONOMETRIC_DOMAIN(a)))
    {
        printf2(COLOR_SYSTEM, "\ncatanh function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
    	printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) < -1 V Re(x) > 1)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "catanh", catanh(a+b*I));
	return;
case BCALC_CACSC:	
	if(dcheck && (!b) && a > -1 && a < 1)
    {
        printf2(COLOR_SYSTEM, "\ncacsc function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) > -1 ^ Re(x) < 1} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cacsc", cacsc(a+b*I));
	return;
case BCALC_CACSCH:
	viewComplexResult(CRV_DODGCHECK, "cacsch", cacsch(a+b*I));
	return;
case BCALC_CASEC:	
	if(dcheck && (!b) && a > -1 && a < 1)
    {
        printf2(COLOR_SYSTEM, "\ncasec function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) > -1 ^ Re(x) < 1} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "casec", casec(a+b*I));
	return;
case BCALC_CASECH:
	if(dcheck && (!b) && (a < 0 || a >= 1))
    {
        printf2(COLOR_SYSTEM, "\ncasech function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ (Re(x) < 0 V Re(x) >= 1)} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "casech", casech(a+b*I));
	return;
case BCALC_CACOT:
	viewComplexResult(CRV_DODGCHECK, "cacot", cacot(a+b*I));
	return;
case BCALC_CACOTH:
	if(dcheck && (!b) && TRIGONOMETRIC_DOMAIN(a))
    {
        printf2(COLOR_SYSTEM, "\ncacoth function isn't defined in: (");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ",0).\n");
        printErr(33, "(DOMAIN: C\\{x| Im(x) = 0 ^ Re(x) >= -1 ^ Re(x) <= -1} )");
    }
    else
		viewComplexResult(CRV_DODGCHECK, "cacoth", cacoth(a+b*I));
	return;

case BCALC_MCD:
    if(dcheck && (a < 1 || b < 1))
    {
        printErr(33, "Invalid inserted Value");
        return;
    }
    c = math_MCD((uint64_t)a, (uint64_t)b);
    break;

case BCALC_MCM:
    if(dcheck && (a < 1 || b < 1))
    {
        printErr(33, "Invalid inserted Value");
        return;
    }
    c = math_mcm((uint64_t)a, (uint64_t)b);
    break;

case BCALC_APPROSSIMAZIONE:
    c = a ? ceil(b) : floor(b);
    break;

case BCALC_SOMMASUCCESSIONEGEOMETRICA:
{
    int64_t exponent;
    exponent = (int64_t) b;
    if(a == 1)
    {
        printErr(33, "Invalid inserted Value");
        return;
    }

    c = gsum(a, exponent);

    break;
}

case BCALC_SOMMASUCCESSIONEARMONICAGEN:
    case BCALC_SOMMASUCCESSIONEARMONICA:
    {
        if(a <= 1)
        {
            printErr(33, "Invalid inserted Value");
            return;
        }
        uint64_t numero;
        uint64_t exponent;

        numero = (uint64_t) a;
        exponent = oprID == BCALC_SOMMASUCCESSIONEARMONICA ? 1 : (uint64_t) b;

        c = gasum(numero, exponent);

        break;
    }

case BCALC_SOMMASUCCESSIONEFIBONACCI:
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value");
        return;
    }

    c = fsum(a);
    break;

case BCALC_SOMMASUCCESSIONEFATTORIALE:
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value");
        return;
    }

    c = fasum(a);
    break;

case BCALC_SOMMASUCCESSIONESEMIFATTORIALE:
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value");
        return;
    }

    c = sfasum(a);
    break;

case BCALC_SOMMATORIA:
{
    static ityp sum = 0.00;
    if(b)
    	c = (sum += a*(INVERSE_OPS?(-1):(1)));
    else
    {
    	c = sum;
    	sum = 0.00;
    }
    break;
}


case BCALC_PRODUTTORIA:
{
    static ityp product = 1.00;
    if(b)
    	c = (product *= mpow2(a, (INVERSE_OPS?(-1):(1))));
    else
    {
    	c = product;
    	product = 1.00;
    }
    break;
}

case BCALC_MEDIA:
{
    static ityp media = 0;
    static uint64_t accumulate = 0;

    if(b)
    {
        media += a;
        ++ accumulate;
        viewInsertedValue;
    }

    c = media/accumulate;
    media = accumulate = 0;

    break;
}

case BCALC_VARIANZA:
{
    static ityp *vector = NULL;
    static uint64_t accumulate = 0;

    if(b)
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        if(!(accumulate))
            vector = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));
        else if(lazy_exec)
            vector = realloc(vector, sizeof(ityp)*(accumulate+1));
        else if(!((accumulate+1) % access(curLayout)->stabilizer_factor))
            vector = realloc(vector, sizeof(ityp)*(accumulate+access(curLayout)->stabilizer_factor));

        errMem(vector, VSPACE);
        vector[accumulate++] = a;
        viewInsertedValue;
    }

    c = math_variance(accumulate, vector);
    free(vector);
    vector = NULL;
    accumulate = 0;

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;
}

case BCALC_STDDEV:
{
    static ityp *vector = NULL;
    static uint64_t accumulate = 0;

    if(b)
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        if(!(accumulate))
            vector = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));
        else if(lazy_exec)
            vector = realloc(vector, sizeof(ityp)*(accumulate+1));
        else if(!((accumulate+1) % access(curLayout)->stabilizer_factor))
            vector = realloc(vector, sizeof(ityp)*(accumulate+access(curLayout)->stabilizer_factor));

        errMem(vector, VSPACE);
        vector[accumulate++] = a;
        viewInsertedValue;
    }

    c = math_stddev(accumulate, vector);
    free(vector);
    vector = NULL;
    accumulate = 0;

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;
}

case BCALC_OUTLIER:
{
    static ityp *vector = NULL;
    static uint64_t accumulate = 0;
    static int64_t outlier_idx = -1;

    if(b == outlier_idx || !accumulate)
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS
        
        if(!accumulate)
        	outlier_idx = b;

        if(!(accumulate))
            vector = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));
        else if(lazy_exec)
            vector = realloc(vector, sizeof(ityp)*(accumulate+1));
        else if(!((accumulate+1) % access(curLayout)->stabilizer_factor))
            vector = realloc(vector, sizeof(ityp)*(accumulate+access(curLayout)->stabilizer_factor));

        errMem(vector, VSPACE);
        
        vector[accumulate++] = a;
        viewInsertedValue;
    }
    
    if(accumulate <= outlier_idx)
    {
    	printErr(33, "You haven't entered at least %llu elements!", outlier_idx);
    	return;
    }
        
    PRINT2N();
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, vector[outlier_idx]);
    printf2(COLOR_USER, " is%s an OUTLIER for the Entered Vector.", math_outlier(accumulate, outlier_idx, vector) ? NULL_CHAR : "n't");
    PRINTN();
    
    free(vector);
    vector = NULL;

    outlier_idx = -1,
    accumulate = 0;

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}

case BCALC_MAP:
{
    static ityp *vector = NULL;
    static uint64_t accumulate = 0;
    static dim_typ funcID = MAX_FIDS;

    if(b == funcID || !accumulate)
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS
        
        if(!accumulate && (b != (funcID = (dim_typ)b) || funcID < FID_SIN || funcID > LAST_FID))
        {
        	printErr(33, "You've entered an Invalid Function ID. Must be an integer between %hu and %hu", FID_SIN, LAST_FID);
        	return;
        }
        

        if(!(accumulate))
            vector = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));
        else if(lazy_exec)
            vector = realloc(vector, sizeof(ityp)*(accumulate+1));
        else if(!((accumulate+1) % access(curLayout)->stabilizer_factor))
            vector = realloc(vector, sizeof(ityp)*(accumulate+access(curLayout)->stabilizer_factor));

        errMem(vector, VSPACE);
        
        vector[accumulate++] = a;
        viewInsertedValue;
    }
        
    PRINT2N();
    SHOWPAUSEMESSAGE();
    PRINTL();
    printf2(COLOR_USER, "The inserted Vector: \n");
	PRINTL();
	
	uint64_t i;
	
    ityp (* const map_func)(ityp) = ext_math.functions[funcID];
	for(i=0; i<accumulate; ++i)
	{
		printf2(COLOR_USER, "- %llu: ", i);
		printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, vector[i]);
		printf2(COLOR_USER, ";\n");
		vector[i] = map_func(vector[i]);
		if(catchPause()) return;
	}
	
	PRINTL();
	printf2(COLOR_USER, "mapped with the %s function results in:\n", ext_math.funcnames[funcID]);
	PRINT2N();
    SHOWPAUSEMESSAGE();
	PRINTL();
	
	for(i=0; i<accumulate; ++i)
	{
		printf2(COLOR_USER, "- %llu: ", i);
		printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, vector[i]);
		printf2(COLOR_USER, ";\n");
		if(catchPause()) return;
	}
	
    free(vector);
    vector = NULL;

	funcID = MAX_FIDS;
    accumulate = 0;

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}


case BCALC_MEDIAGEOMETRICA:
{
    static ityp media = 0;
    static ityp accumulate = 0;

    if(!(a))
    {
        c = pow(media, 1/accumulate);
        media = accumulate = 0;
        break;
    }

    media *= a;
    ++ accumulate;
    viewInsertedValue;

}

case BCALC_MEDIAARMONICA:
{
    static ityp media = 0;
    static uint64_t accumulate = 0;

    if(!(a))
    {
        c = accumulate/media;
        media = accumulate = 0;
        break;
    }

    media += 1/a;
    ++ accumulate;
    viewInsertedValue;
}

case BCALC_MEDIAPOTENZA:
{
    static ityp media = 1;
    static ityp accumulate = 0;

    if(!(a) && !(b))
    {
        c = pow(media, 1/b);
        media = 1;
        accumulate = 0;
        break;
    }

    if(dcheck && !(a))
    {
        printErr(33, "Invalid inserted Value. 0 hasn't weight in this Operation");
        return;
    }

    if(dcheck && !(b))
    {
        printErr(33, "Invalid inserted Value.\nRoot must be a non-zero integer");
        return;
    }

    media += pow(a,b);
    ++ accumulate;
    viewInsertedValue;
}

case BCALC_VALORECENTRALE:
{

    static ityp *vector = NULL;
    static uint64_t accumulate = 0;

    if(a)
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        if(!(accumulate))
            vector = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));
        else if(lazy_exec)
            vector = realloc(vector, sizeof(ityp)*(accumulate+1));
        else if(!((accumulate+1) % access(curLayout)->stabilizer_factor))
            vector = realloc(vector, sizeof(ityp)*(accumulate+access(curLayout)->stabilizer_factor));

        errMem(vector, VSPACE);
        vector[accumulate++] = a;
        viewInsertedValue;
    }

    c = math_scale(accumulate, vector);
    free(vector);
    vector = NULL;
    accumulate = 0;

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;
}

case BCALC_PRIMOQUARTILE:
{

    static ityp *vector;
    static uint64_t accumulate = 0;

    if(a)
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        if(!(accumulate))
            vector = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));
        else if(lazy_exec)
            vector = realloc(vector, sizeof(ityp)*(accumulate+1));
        else if(!((accumulate+1) % access(curLayout)->stabilizer_factor))
            vector = realloc(vector, sizeof(ityp)*(accumulate+access(curLayout)->stabilizer_factor));

        errMem(vector, VSPACE);
        vector[accumulate++] = a;
        viewInsertedValue;
    }

    c = math_first_quartile(accumulate, vector);
    accumulate = 0;
    free(vector);
    vector = NULL;

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;

}

case BCALC_MEDIANA:
{

    static ityp *vector;
    static uint64_t accumulate = 0;

    if(a)
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        if(!(accumulate))
            vector = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));
        else if(lazy_exec)
            vector = realloc(vector, sizeof(ityp)*(accumulate+1));
        else if(!((accumulate+1) % access(curLayout)->stabilizer_factor))
            vector = realloc(vector, sizeof(ityp)*(accumulate+access(curLayout)->stabilizer_factor));

        errMem(vector, VSPACE);
        vector[accumulate++] = a;
        viewInsertedValue;
    }

    c = math_mediana(accumulate, vector);
    accumulate = 0;
    free(vector);
    vector = NULL;

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;

}

case BCALC_TERZOQUARTILE:
{

    static ityp *vector;
    static uint64_t accumulate = 0;

    if(a)
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        if(!(accumulate))
            vector = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));
        else if(lazy_exec)
            vector = realloc(vector, sizeof(ityp)*(accumulate+1));
        else if(!((accumulate+1) % access(curLayout)->stabilizer_factor))
            vector = realloc(vector, sizeof(ityp)*(accumulate+access(curLayout)->stabilizer_factor));

        errMem(vector, VSPACE);
        vector[accumulate++] = a;
        viewInsertedValue;
    }

    c = math_third_quartile(accumulate, vector);
    accumulate = 0;
    free(vector);
    vector = NULL;

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;

}

case BCALC_SOMMAPRIMINNUMERI:
    if(a < 0)
    {
        printErr(1, "Invalid inserted Value");
        return;
    }
    printf2(COLOR_USER, "\nFirst %llu Natural Numbers SUM RESULT is: %llu.\n\n", (uint64_t) a, fnnsum(a));
    return;

case BCALC_FATTORIALE:
{
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value\nMust be a non-negative and non-zero integer");
        return;
    }

    uint64_t numero;
    numero = (uint64_t) a;

    if(b)
    {
        uint64_t i;

        printf2(COLOR_SYSTEM, "\nFirst %llu FACTORIAL Succession Numbers are:\n\n", numero);
        SHOWPAUSEMESSAGE();
        PRINTL();

        for(i=0; i<=numero; ++i)
        {
            printf2(COLOR_USER, "%llu: ", i+1);
            printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, fact(i));
            printf2(COLOR_USER, ";\n"); // tab[i] = fibo(i);
            if(catchPause()) return;
        }

        PRINTL();

        return;
    }

    c = fact(a);
    break;

}

case BCALC_SEMIFATTORIALE:
{
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value\nMust be a non-negative and non-zero integer");
        return;
    }

    uint64_t numero;
    numero = (uint64_t) a;

    if(b)
    {
        uint64_t i;

        printf2(COLOR_SYSTEM, "\nFirst %llu SEMI-FACTORIAL Succession Numbers are:\n\n", numero);
        SHOWPAUSEMESSAGE();
        PRINTL();

        for(i=0; i<=numero; ++i)
        {
            printf2(COLOR_USER, "%llu: ", i+1);
            printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, sfact(i));
            printf2(COLOR_USER, ";\n"); // tab[i] = fibo(i);
            if(catchPause()) return;
        }

        PRINTL();

        return;
    }

    c = sfact(a);
    break;

}

case BCALC_STIRLING:
{
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value\nMust be a non-negative and non-zero integer");
        return;
    }

    uint64_t numero;
    numero = (uint64_t) a;

    if(b)
    {
        uint64_t i;

        printf2(COLOR_SYSTEM, "\nFirst %llu STIRLING's APPROXIMATION Succession Numbers are:\n\n", numero);
        SHOWPAUSEMESSAGE();
        PRINTL();

        for(i=0; i<=numero; ++i)
        {
            printf2(COLOR_USER, "%llu: ", i+1);
            printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, stirling(i));
            printf2(COLOR_USER, ";\n"); // tab[i] = fibo(i);
            if(catchPause()) return;
        }

        PRINTL();

        return;
    }

    c = stirling(a);
    break;

}

case BCALC_FIBONACCI:
{
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value.\nMust be a non-negative and non-zero integer");
        return;
    }

    uint64_t numero;
    numero = (uint64_t) a;

    if(b)
    {
        uint64_t i;

        printf2(COLOR_SYSTEM, "\nFirst %llu FIBONACCI Succession Numbers are:\n\n", numero+1);
        SHOWPAUSEMESSAGE();
        PRINTL();

        for(i=0; i<numero; ++i)
        {
            printf2(COLOR_USER, "%llu: ", i+1);
            printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, fibo(i));
            printf2(COLOR_USER, ";\n"); // tab[i] = fibo(i);
            if(catchPause()) return;
        }

        PRINTL();
        return;
    }

    c = fibo(a);

    break;

}

case BCALC_PERMUTATIONS:
    c = perm(a);
    break;
    
case BCALC_PERMUTATIONSREP:
{
	static ityp denominator = 1.00;
	static ityp sum = 0.00;
	static uint64_t accumulate = 0;

    if(b)
    {
        denominator *= fact(b);
        sum += b;
        ++ accumulate;
        viewInsertedValue;
    }
    
    if(dcheck && sum != a)
    {
    	printErr(33, "You've inserted invalid multiplicities, because their sum isn't %.*f", access(curLayout)->precision, a);
		denominator = 1.00;
		sum = 0.00;
		return;
	}	
    
    c = a/denominator;
    denominator = 1.00,
    sum = accumulate = 0;

    break;
}
    
case BCALC_KPERMUTATIONS:
	c = kperm(a, b);
    
case BCALC_KPERMUTATIONSREP:
	c = kperm_rep(a, b);

case BCALC_COMBINATIONS:
    c = comb(a, b);
    break;
    
case BCALC_COMBINATIONSREP:
	c = comb_rep(a, b);

case BCALC_CAMBIAMENTODIBASE:

    if(dcheck && a > access(curLayout)->max_changebase_binary_convnum)
    {
        printErr(33, "Invalid inserted Value.\nMust be an integer between 0 and %d", access(curLayout)->max_changebase_binary_convnum);
        return;
    }
    if(dcheck && (b < access(curLayout)->basecalc_minbase || b > access(curLayout)->basecalc_maxbase))
    {
        printErr(33, "Invalid inserted Value.\nBase must be an integer between %d and %d",
                access(curLayout)->basecalc_minbase, access(curLayout)->basecalc_maxbase);
        return;
    }

    if(dcheck && b > 10)
        printf2(COLOR_SYSTEM, "\nWARNING: It is highly recommended not to handle\nOUTPUT Conversion Basis bigger than 10.\n\n");
    c = changeBase((int64_t)a, DEFAULT_BASE, b);
    break;

case BCALC_GENERATOREMATRICIRANDOM:
{

    if(dcheck && (a < 1 || b < 1 || a > access(curLayout)->matrix_max_rows || b > access(curLayout)->matrix_max_columns))
    {
        printErr(33, "Invalid [ROWS COLUMNS] format.\nYou have to insert non-negative and non-zero ROWS and COLUMNS.\n\
and must be respectively less than: %hu and %hu", access(curLayout)->matrix_max_rows, access(curLayout)->matrix_max_columns);
        return;
    }

    dim_typ2 dim =
    {
        (dim_typ) a,
        (dim_typ) b
    };
    ityp *matrix = NULL;

    if(!matrixAlloc(&matrix, dim))
        return;

    randomMatrix(matrix, dim);
    matrixFree(&matrix);
    return;
}

case BCALC_PASCALTRIANGLE:
{
    uint64_t numero;
    numero = (uint64_t) a;
    if(a < access(curLayout)->pascal_triangle_min_rows || a > access(curLayout)->pascal_triangle_max_rows)
    {
        printErr(33, "Invalid inserted Value.\nMust be an integer between %u and %u", access(curLayout)->pascal_triangle_min_rows, access(curLayout)->pascal_triangle_max_rows);
        return;
    }

    getPascalTriangle(numero+1);
    return;
}

case BCALC_NUMERIROMANI:
{
    const dim_typ numero = (dim_typ) a;

    if(dcheck && (a < access(curLayout)->min_roman_number || a > access(curLayout)->max_roman_number))
    {
        printErr(33, "Invalid inserted Value.\nMust be an integer between %hu and %hu", access(curLayout)->min_roman_number, access(curLayout)->max_roman_number);
        return;
    }

    char str[ROMAN_NUMBER_STRING] = NULL_CHAR;
    getRomanNumber(numero, str);
    printf2(COLOR_USER, "\nInserted Value: %hu ROMAN NUMBER Conversion is:\n%s.\n\n", numero, str);
    return;
}

case BCALC_PRIMINNUMERIPRIMI:
{
    if(a > b || a < 2 || b < 2)
    {
        printErr(33, "Invalid inserted Value\nBoth Interval Params must be >= 3.\n");
        return;
    }

    uint64_t startpoint;
    uint64_t endpoint;

    startpoint = (uint64_t) a;
    endpoint = (uint64_t) b;
    printf2(COLOR_SYSTEM, "\nFIRST PRIME NUMBERS between %llu and %llu are:\n\n", startpoint, endpoint);
    prime_N_Number(startpoint, endpoint);
    return;
}

case BCALC_NESIMONUMEROPRIMO:
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value.\nMust be a non-negative and non-zero integer");
        return;
    }
    c = N_prime_Number(a); // numero);
    break;


case BCALC_PRIMORIALE:
    c = primr(a);
    break;

case BCALC_SOMMAPRIMINNUMERIPRIMI:
    c = fpnsum(a);
    break;

case BCALC_FIBONACCIALE:
    c = fibnc(a);
    break;

case BCALC_INFORMAZIONI:
    printOpersIdentifiers();
    return;
