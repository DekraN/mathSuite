/*

BASECALC support header
This has been created in order
to prevent programs.c getting too crowded.
*/


case BASECALC_EXIT:
    safeExit(a);
    break;

case BASECALC_ADDIZIONE:
    c = (a + b);
    break;

case BASECALC_SOTTRAZIONE:
    c = (a - b);
    break;

case BASECALC_MOLTIPLICAZIONE:
    c = (a * b);
    break;

case BASECALC_DIVISIONE:
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

case BASECALC_RESTO:
    if(dcheck && (!b))
    {
        printErr(33, "You are trying to get Rest from a 'per-zero DIVISION");
        return;
    }
    c = ((int)a % (int)b);
    break;

case BASECALC_ADDIZIONEBINARIA:
{
    char *c = NULL;
    if((c = binaryAlgSum(a, b, BAS_SUM)))
        free(c);
    return;
}

case BASECALC_SOTTRAZIONEBINARIA:
{
    char *c = NULL;
    if((c = binaryAlgSum(a, b, BAS_SUB)))
       free(c);
    return;
}

case BASECALC_COMPLEMENTO:
{
    char *c = NULL;
    if((c = binNumComp(a)))
        free(c);
    return;
}

case BASECALC_ELEVAMENTOAPOTENZA:
    if(dcheck && (!a) && (!b))
    {
        printf2(COLOR_SYSTEM, "\nRESULT: Indeterminate Form 0^0.\n");
        return;
    }
    c = pow(a, b);
    break;

case BASECALC_EXPANDEXPC:
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

case BASECALC_EXP10ANDEXP10C:
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


case BASECALC_EXP2ANDEXP2C:
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

case BASECALC_RADICENESIMA:
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

case BASECALC_RADICEQUADRATA:
    if(dcheck && a < 0)
    {
        printErr(33, "Invalid Quad Root Function Argument.\n(DOMAIN: [0,+inf] )");
        return;
    }
    c = sqrt(a);
    break;

case BASECALC_RADICECUBICA:
    c = cbrt(a);
    break;

case BASECALC_LOGARITMO:

    if(dcheck && b <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlog function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
        printf2(COLOR_SYSTEM, ". ");
        printErr(33, "(DOMAIN: ]0,+inf[ )");
        return;
    }

    c = (a ? log10(b) : log(b));
    break;

case BASECALC_LOGARITMO2:

    if(dcheck && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlog2 function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ". ");
        printErr(33, "(DOMAIN: ]0,+inf[ )");
        return;
    }

    c = log2(a);
    break;

case BASECALC_LOGARITMOBN:

    if(dcheck && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlogBN functions isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ". ");
        printErr(33, "(DOMAIN: ]0,+inf[ )");
        return;
    }

    c = logbN(a, b);
    break;

case BASECALC_LOGARITMOC:

    if(dcheck && b <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlogc function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
        printf2(COLOR_SYSTEM, ". ");
        printErr(33, "(DOMAIN: ]0,+inf[ )");
        return;
    }

    c = (a ? log10c(b) : logc(b));
    break;

case BASECALC_LOGARITMO2C:

    if(dcheck && a <= 0)
    {
        printf2(COLOR_SYSTEM, "\nlog2c function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ". ");
        printErr(33, "(DOMAIN: ]0,+inf] )");
        return;
    }

    c = log2c(a);
    break;

case BASECALC_LOGARITMO1P:

    if(dcheck && a <= -1)
    {
        printf2(COLOR_SYSTEM, "\nlog1p function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ". ");
        printErr(33, "(DOMAIN: ]-1,+inf[ )");
        return;
    }

    c = log1p(a);
    break;

case BASECALC_LOGARITMO1PC:

    if(dcheck && a <= -1)
    {
        printf2(COLOR_SYSTEM, "\nlog2c function isn't defined in: ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, a);
        printf2(COLOR_SYSTEM, ". ");
        printErr(33, "(DOMAIN: ]-1,+inf] )");
        return;
    }

    c = log1pc(a);
    break;

case BASECALC_BITCOUNTER:
    c = countbits(a);
    break;

case BASECALC_UBITCOUNTER:
    c = ucountbits(a);
    break;

case BASECALC_VERSION:
    c = strtod(PROG__VERSION, NULL);
    break;
    
case BASECALC_EXITCHAR:
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

case BASECALC_PREC:

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


case BASECALC_SFACT:
	
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
    
case BASECALC_MINSRNUMBER:
	
	if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MIN_STIRLINGREQUIRES_NUMBER)
        {
            printErr(33, "Invalid Inserted Stabilizer Factor Value.\nMust be a non-negative integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MIN_STIRLINGREQUIRES_NUMBER);
            return;
        }
        c = access(curLayout)->min_stirlingrequires_number = a;
    }
    else
        c = access(curLayout)->min_stirlingrequires_number;
    break;

case BASECALC_ALGEBRA:
    if(a)
    {
    if(a < MIN_ALGEBRA || a > MAX_ALGEBRA)
        {
            printErr(33, "Invalid Inserted Algebra Identifier.\nMust be a non-negative integer between %hu and %hu", MIN_ALGEBRA, MAX_ALGEBRA);
            return;
        }
        c = a;
        _changeAlgebraDims(a);
    }
    else
        c = access(curLayout)->algebra;
    break;
    
case BASECALC_OUTLIERCONST:
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
	

case BASECALC_RSEED:

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

case BASECALC_MMIFIBO:

    if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MAX_FIBONACCI_MEMOIZABLE_INDEX)
        {
            printErr(33, "Invalid inserted %s MIM Value.\nMust be a non-negative integer between %hu and %hu", suite_c.memoizers_names[FUNCTION_FIBONACCI], MIN_MEMOIZABLE_INDEX+1, MAX_FIBONACCI_MEMOIZABLE_INDEX);;
            return;
        }
        c = access(curLayout)->max_memoizable_indeces[FUNCTION_FIBONACCI] = a;
    }
    else
        c = access(curLayout)->max_memoizable_indeces[FUNCTION_FIBONACCI];
    break;

case BASECALC_MMIFACT:

    if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MAX_FATTORIALE_MEMOIZABLE_INDEX)
        {
            printErr(33, "Invalid inserted %s MIM Value.\nMust be a non-negative integer between %hu and %hu", suite_c.memoizers_names[FUNCTION_FATTORIALE], MIN_MEMOIZABLE_INDEX+1, MAX_FATTORIALE_MEMOIZABLE_INDEX);
            return;
        }
        c = access(curLayout)->max_memoizable_indeces[FUNCTION_FATTORIALE] = a;
    }
    else
        c = access(curLayout)->max_memoizable_indeces[FUNCTION_FATTORIALE];
    break;
    
case BASECALC_MMIEVENSFACT:

    if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX)
        {
            printErr(33, "Invalid inserted %s MIM Value.\nMust be a non-negative integer between %hu and %hu", suite_c.memoizers_names[FUNCTION_EVEN_SFATTORIALE], MIN_MEMOIZABLE_INDEX+1, MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX);
            return;
        }
        c = access(curLayout)->max_memoizable_indeces[FUNCTION_EVEN_SFATTORIALE] = a;
    }
    else
        c = access(curLayout)->max_memoizable_indeces[FUNCTION_EVEN_SFATTORIALE];
    break;
    
case BASECALC_MMIODDSFACT:

    if(a)
    {
        if(a < MIN_MEMOIZABLE_INDEX+1 || a > MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX)
        {
            printErr(33, "Invalid inserted %s MIM Value.\nMust be a non-negative integer between %hu and %hu", suite_c.memoizers_names[FUNCTION_ODD_SFATTORIALE], MIN_MEMOIZABLE_INDEX+1, MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX);
            return;
        }
        c = access(curLayout)->max_memoizable_indeces[FUNCTION_ODD_SFATTORIALE] = a;
    }
    else
        c = access(curLayout)->max_memoizable_indeces[FUNCTION_ODD_SFATTORIALE];
    break;


case BASECALC_TRASFORMAANGOLI:
    c = (a ? deg(b) : rad(b));
    break;

case BASECALC_SINANDSINH:
    c = (a ? sinh(b) : sin(b));
    break;

case BASECALC_COSANDCOSH:
    c = (a  ? cosh(b) : cos(b));
    break;

case BASECALC_TANANDTANH:
    if(a) c = tanh(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ntan function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x| x = PI/2 + k*PI} )");
            return;
        }
        c = tan(b);
    }
    break;

case BASECALC_CSCANDCSCH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\ncsch function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = csc(b);
    }
    break;

case BASECALC_SECANDSECH:
    if(a) c = sech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nsec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = sec(b);
    }
    break;

case BASECALC_COTANDCOTH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\ncoth function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = cot(b);
    }
    break;

case BASECALC_HSINANDHSINH:
    c = (a ? hsinh(b) : hsin(b));
    break;

case BASECALC_QSINANDQSINH:
    c = (a ? qsinh(b) : qsin(b));
    break;

case BASECALC_HCOSANDHCOSH:
    c = (a ? hcosh(b) : hcos(b));
    break;

case BASECALC_QCOSANDQCOSH:
    c = (a ? qcosh(b) : qcos(b));
    break;

case BASECALC_HSECANDHSECH:
    if(a) c = hsech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhsec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = hsec(b);
    }
    break;

case BASECALC_QSECANDQSECH:
    if(a) c = qsech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqsec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = qsec(b);
    }
    break;

case BASECALC_HCSCANDHCSCH:
    if(a) c = hcsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhcsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hcsc(b);
    }
    break;

case BASECALC_QCSCANDQCSC:
    if(a) c = qcsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqcsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qcsc(b);
    }
    break;

case BASECALC_HTANANDHTANH:
    if(a) c = htanh(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhtan function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x| x = PI/2 + k*PI} )");
            return;
        }
        c = htan(b);
    }
    break;

case BASECALC_QTANANDQTANH:
    if(a) c = qtanh(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqtan function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x| x = PI/2 + k*PI} )");
            return;
        }
        c = qtan(b);
    }
    break;

case BASECALC_HCOTANDHCOTH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\nhcoth function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hcot(b);
    }
    break;

case BASECALC_QCOTANDQCOTH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\nqcoth function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qcot(b);
    }
    break;

case BASECALC_VSINANDVSINH:
    c = (a ? vsinh(b) : vsin(b));
    break;

case BASECALC_CVSINANDCVSINH:
    c = (a ? cvsinh(b) : cvsin(b));
    break;

case BASECALC_VCOSANDVCOSH:
    c = (a ? vcosh(b) : vcos(b));
    break;

case BASECALC_CVCOSANDCVCOSH:
    c = (a ? cvcosh(b) : cvcos(b));
    break;

case BASECALC_HVSINANDHVSINH:
    c = (a ? hvsinh(b) : hvsin(b));
    break;

case BASECALC_HCVSINANDHCVSINH:
    c = (a ? hcvsinh(b) : hcvsin(b));
    break;

case BASECALC_QVSINANDQVSINH:
    c = (a ? qvsinh(b) : qvsin(b));
    break;

case BASECALC_QCVSINANDQCVSINH:
    c = (a ? qcvsinh(b) : qcvsin(b));
    break;

case BASECALC_HVCOSANDHVCOSH:
    c = (a ? hvcosh(b) : hvcos(b));
    break;

case BASECALC_HCVCOSANDHCVCOSH:
    c = (a ? hcvcosh(b) : hcvcos(b));
    break;

case BASECALC_QVCOSANDQVCOSH:
    c = (a ? qvcosh(b) : qvcos(b));
    break;

case BASECALC_QCVCOSANDQCVCOSH:
    c = (a ? qcvcosh(b) : qcvcos(b));
    break;

case BASECALC_ESECANDESECH:
    if(a) c = esech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nesec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = esec(b);
    }
    break;

case BASECALC_ECSCANDECSCH:
    if(a) c = ecsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\necsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = ecsc(b);
    }
    break;

case BASECALC_HESECANDHESECH:
    if(a) c = hesech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhesec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = hesec(b);
    }
    break;

case BASECALC_HECSCANDHECSCH:
    if(a) c = hecsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhecsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hecsc(b);
    }
    break;

case BASECALC_QESECANDQESECH:
    if(a) c = qesech(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqesec function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = PI/2 + k*PI} )");
            return;
        }
        c = qesec(b);
    }
    break;

case BASECALC_QECSCANDQECSCH:
    if(a) c = qecsch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqecsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qecsc(b);
    }
    break;

case BASECALC_SINCANDSINCH:
    c = (a ? sinch(b) : sinc(b));
    break;

case BASECALC_HSINCANDHSINCH:
    c = (a ? hsinch(b) : hsinc(b));
    break;

case BASECALC_QSINCANDQSINCH:
    c = (a ? qsinch(b) : qsinc(b));
    break;

case BASECALC_COSCANDCOSCH:
    if(a) c = cosch(b);
    else
    {
        if(dcheck && !b)
        {
            printf2(COLOR_SYSTEM, "\ncosc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0} )");
            return;
        }
        c = cosc(b);
    }
    break;

case BASECALC_HCOSCANDHCOSCH:
    if(a) c = hcosch(b);
    else
    {
        if(dcheck && !b)
        {
            printf2(COLOR_SYSTEM, "\nhcosc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0} )");
            return;
        }
        c = hcosc(b);
    }
    break;

case BASECALC_QCOSCANDQCOSCH:
    if(a) c = qcosch(b);
    else
    {
        if(dcheck && !b)
        {
            printf2(COLOR_SYSTEM, "\nqcosc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0} )");
            return;
        }
        c = qcosc(b);
    }
    break;

case BASECALC_SECCANDSECCH:
    if(a) c = secch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nsecc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = secc(b);
    }
    break;

case BASECALC_HSECCANDHSECCH:
    if(a) c = hsecch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhsecc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = hsecc(b);
    }
    break;

case BASECALC_QSECCANDQSECCH:
    if(a) c = qsecch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqsecc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = qsecc(b);
    }
    break;

case BASECALC_CSCCANDCSCCH:
    if(a) c = cscch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ncscc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = cscc(b);
    }
    break;

case BASECALC_HCSCCANDHCSCCH:
    if(a) c = hcscch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhcscc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hcscc(b);
    }
    break;

case BASECALC_QCSCCANDQCSCCH:
    if(a) c = qcscch(b);
    else
    {
        if(dcheck && (!trigonometric_domain(b, 0, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqcscc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qcscc(b);
    }
    break;

case BASECALC_TANCANDTANCH:
    if(a) c = tanch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\ntanc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = tanc(b);
    }
    break;

case BASECALC_HTANCANDHTANCH:
    if(a) c = htanch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nhtanc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = htanc(b);
    }
    break;

case BASECALC_QTANCANDQTANCH:
    if(a) c = qtanch(b);
    else
    {
        if(dcheck && ((!b) || !trigonometric_domain(b, M_PI_2, M_PI)))
        {
            printf2(COLOR_SYSTEM, "\nqtanc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = 0 V x = PI/2 + k*PI} )");
            return;
        }
        c = qtanc(b);
    }
    break;

case BASECALC_COTCANDCOTCH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\ncotch function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = cotc(b);
    }
    break;

case BASECALC_HCOTCANDHCOTCH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\nhcotch function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = hcotc(b);
    }
    break;

case BASECALC_QCOTCANDQCOTCH:
    if(a)
    {
        if(dcheck && (!b))
        {
            printf2(COLOR_SYSTEM, "\nqcotch function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: R\\{x|x = k*PI} )");
            return;
        }
        c = qcotc(b);
    }
    break;


case BASECALC_ASINANDASINH:
    if(a == 0)
    {
        if(dcheck && !(TRIGONOMETRIC_DOMAIN(b)))
        {
            printf2(COLOR_SYSTEM, "\nasin function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: [-1,1] )");
            return;
        }
        c = asin(b);
    }
    else c = asinh(b);
    break;

case BASECALC_ACOSANDACOSH:
    if(a)
    {

        if(dcheck && b < 1)
        {
            printf2(COLOR_SYSTEM, "\nacosh function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: [-1,1] )");
            return;
        }
        c = acos(b);
    }
    break;

case BASECALC_ATANANDATANH:
    if(a)
    {
        if(dcheck && !(TRIGONOMETRIC_DOMAIN(b)))
        {
            printf2(COLOR_SYSTEM, "\natanh function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: [-1,1] )");
            return;
        }
        c = atanh(b);
    }
    else c = atan(b);
    break;

case BASECALC_ATAN2:
    if(dcheck && (!a) && (!b))
    {
        printf2(COLOR_SYSTEM, "\natan2 function isn't defined in: (0,0). ");
        printErr(33, "(DOMAIN: (R^2)\\{P=(y,x)|P=(0,0)} )");
        return;
    }
    c = atan2(a, b);
    break;

case BASECALC_ACSCANDACSCH:
    if(a) c = acsch(b);
    else
    {
        if(dcheck && b > -1 && b < 1)
        {
            printf2(COLOR_SYSTEM, "\nacsc function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: ]-inf,-1]U[1,+inf[ )");
            return;
        }
        c = acsc(b);
    }
    break;

case BASECALC_ASECANDASECH:
    if(a)
    {
        if(dcheck && (b < 0 || b >= 1))
        {
            printf2(COLOR_SYSTEM, "\nasech function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
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
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: ]-inf,-1]U[1,+inf[ )");
            return;
        }
        c = asec(b);
    }
    break;

case BASECALC_ACOTANDACOTH:
    if(a)
    {
        if(dcheck && TRIGONOMETRIC_DOMAIN(b))
        {
            printf2(COLOR_SYSTEM, "\nacoth function isn't defined in: ");
            printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, b);
            printf2(COLOR_SYSTEM, ". ");
            printErr(33, "(DOMAIN: ]-inf,-1[U]1,+inf[ )");
            return;
        }
        c = acoth(b);
    }
    else c = acot(b);
    break;

case BASECALC_MCD:
    if(dcheck && (a < 1 || b < 1))
    {
        printErr(33, "Invalid inserted Value");
        return;
    }
    c = math_MCD((uint64_t)a, (uint64_t)b);
    break;

case BASECALC_MCM:
    if(dcheck && (a < 1 || b < 1))
    {
        printErr(33, "Invalid inserted Value");
        return;
    }
    c = math_mcm((uint64_t)a, (uint64_t)b);
    break;

case BASECALC_APPROSSIMAZIONE:
    c = a ? ceil(b) : floor(b);
    break;

case BASECALC_SOMMASUCCESSIONEGEOMETRICA:
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

case BASECALC_SOMMASUCCESSIONEARMONICAGEN:
    case BASECALC_SOMMASUCCESSIONEARMONICA:
    {
        if(a <= 1)
        {
            printErr(33, "Invalid inserted Value");
            return;
        }
        uint64_t numero;
        uint64_t exponent;

        numero = (uint64_t) a;
        exponent = oprID == BASECALC_SOMMASUCCESSIONEARMONICA ? 1 : (uint64_t) b;

        c = gasum(numero, exponent);

        break;
    }

case BASECALC_SOMMASUCCESSIONEFIBONACCI:
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value");
        return;
    }

    c = fsum(a);
    break;

case BASECALC_SOMMASUCCESSIONEFATTORIALE:
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value");
        return;
    }

    c = fasum(a);
    break;

case BASECALC_SOMMASUCCESSIONESEMIFATTORIALE:
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value");
        return;
    }

    c = sfasum(a);
    break;

case BASECALC_SOMMATORIA:
{
    static ityp sum = 0;
    c = (sum += (b ? a*(INVERSE_OPS?(-1):(1)) : 0));
    // c = sum = (sum+(a*(INVERSE_OPS ? (-1):(1)))) : 0);
    break;
}


case BASECALC_PRODUTTORIA:
{
    static ityp product = 1;
    c = (product *= (b ? mpow2(a, (INVERSE_OPS?(-1):(1))) : 1));
    break;
}

case BASECALC_MEDIA:
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

case BASECALC_VARIANZA:
{
    static ityp *vector = NULL;
    static uint64_t accumulate = 0;

    if(b)
    {

        #if WINOS
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

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;
}

case BASECALC_STDDEV:
{
    static ityp *vector = NULL;
    static uint64_t accumulate = 0;

    if(b)
    {

        #if WINOS
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

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;
}

case BASECALC_OUTLIER:
{
    static ityp *vector = NULL;
    static uint64_t accumulate = 0;
    static int64_t outlier_idx = -1;

    if(b == outlier_idx || !accumulate)
    {

        #if WINOS
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

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}

case BASECALC_MAP:
{
    static ityp *vector = NULL;
    static uint64_t accumulate = 0;
    static dim_typ funcID = MAX_FIDS;

    if(b == funcID || !accumulate)
    {

        #if WINOS
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

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}


case BASECALC_MEDIAGEOMETRICA:
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

case BASECALC_MEDIAARMONICA:
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

case BASECALC_MEDIAPOTENZA:
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

case BASECALC_VALORECENTRALE:
{

    static ityp *vector = NULL;
    static uint64_t accumulate = 0;

    if(a)
    {

        #if WINOS
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

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;
}

case BASECALC_PRIMOQUARTILE:
{

    static ityp *vector;
    static uint64_t accumulate = 0;

    if(a)
    {

        #if WINOS
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

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;

}

case BASECALC_MEDIANA:
{

    static ityp *vector;
    static uint64_t accumulate = 0;

    if(a)
    {

        #if WINOS
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

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;

}

case BASECALC_TERZOQUARTILE:
{

    static ityp *vector;
    static uint64_t accumulate = 0;

    if(a)
    {

        #if WINOS
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

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    break;

}

case BASECALC_SOMMAPRIMINNUMERI:
    if(a < 0)
    {
        printErr(1, "Invalid inserted Value");
        return;
    }
    printf2(COLOR_USER, "\nFirst %llu Natural Numbers SUM RESULT is: %llu.\n\n", (uint64_t) a, fnnsum(a));
    return;

case BASECALC_FATTORIALE:
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

case BASECALC_SEMIFATTORIALE:
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

case BASECALC_STIRLING:
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

case BASECALC_FIBONACCI:
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
case BASECALC_PERMUTATION:
    c = perm(a, b);
    break;

case BASECALC_COMBINATION:
    c = comb(a, b);
    break;

case BASECALC_CAMBIAMENTODIBASE:

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

case BASECALC_GENERATOREMATRICIRANDOM:
{

    if(dcheck && (a < 1 || b < 1 || a > access(curLayout)->matrix_max_raws || b > access(curLayout)->matrix_max_columns))
    {
        printErr(33, "Invalid [RAWS COLUMNS] format.\nYou have to insert non-negative and non-zero RAWS and COLUMNS.\n\
and must be respectively less than: %hu and %hu", access(curLayout)->matrix_max_raws, access(curLayout)->matrix_max_columns);
        return;
    }

    dim_typ2 dim =
    {
        (dim_typ) a,
        (dim_typ) b
    };
    ityp **matrix = NULL;

    if(!matrixAlloc(&matrix, dim))
        return;

    randomMatrix(matrix, dim);
    matrixFree(&matrix, dim[RAWS]);
    return;
}

case BASECALC_PASCALTRIANGLE:
{
    uint64_t numero;
    numero = (uint64_t) a;
    if(a < access(curLayout)->pascal_triangle_min_raws || a > access(curLayout)->pascal_triangle_max_raws)
    {
        printErr(33, "Invalid inserted Value.\nMust be an integer between %u and %u", access(curLayout)->pascal_triangle_min_raws, access(curLayout)->pascal_triangle_max_raws);
        return;
    }

    getPascalTriangle(numero+1);
    return;
}

case BASECALC_NUMERIROMANI:
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

case BASECALC_PRIMINNUMERIPRIMI:
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

case BASECALC_NESIMONUMEROPRIMO:
    if(a < 1)
    {
        printErr(33, "Invalid inserted Value.\nMust be a non-negative and non-zero integer");
        return;
    }
    c = N_prime_Number(a); // numero);
    break;


case BASECALC_PRIMORIALE:
    c = primr(a);
    break;

case BASECALC_SOMMAPRIMINNUMERIPRIMI:
    c = fpnsum(a);
    break;

case BASECALC_FIBONACCIALE:
    c = fibnc(a);
    break;

case BASECALC_INFORMAZIONI:
{
    printOpersIdentifiers();
    return;
}
