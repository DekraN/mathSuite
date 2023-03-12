// adv_calc.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"


#ifndef __DISABLE_CRYPTOGRAPHICHASH
__MSSHELL_WRAPPER_ static void  cryptographicHashFunctions(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __apnt cryptographicHashFunctions(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_CRYPTOGRAPHICHASH_PROGS, cryptographic_hash_prog, adv_calc[ADVCALC_CRYPTOGRAPHICHASH].name, BY_CHARS); // MAX_MATMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return;
}
#endif

__MSSHELL_WRAPPER_ static void  secondGradeEquationSolver(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  complexAdd(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  complexMul(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  getDate(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  polynomEval(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  getDate(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  absoluteOrientation(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  routhTable(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  juryTable(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  newtonDifferenceTables(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  lagrangeInterpolation(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  funcIntegration(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  straightLineFitting(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  parabolicCurveFitting(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  linearSystemsSolver(const sel_typ argc, char ** argv);


sprog adv_calc[MAX_ADVCALC_PROGS] =
{
	[ADVCALC_CRYPTOGRAPHICHASH] =
    {
    	CMD_CRYPTOGRAPHICHASH,
        NAME_CRYPTOGRAPHICHASH,
        USAGE_CRYPTOGRAPHICHASH,
        #ifndef __DISABLE_DATABASE
        LEVEL_CRYPTOGRAPHICHASH,
        #endif
        cryptographicHashFunctions,
        ARGC_CRYPTOGRAPHICHASH,
        AUTOMATIC,
        FATHER
    },
    [ADVCALC_SECONDGRADEEQUATIONSOLVER] =
    {
    	CMD_SECONDGRADEQSOLVER,
        NAME_SECONDGRADEQSOLVER,
        USAGE_SECONDGRADEQSOLVER,
        #ifndef __DISABLE_DATABASE
        LEVEL_SECONDGRADEQSOLVER,
        #endif
        secondGradeEquationSolver,
        ARGC_SECONDGRADEQSOLVER,
        BY_USER,
        CHILD
    },
    [ADVCALC_COMPLEXNUMBERSSUM] =
    {
    	CMD_COMPLEXADD,
        NAME_COMPLEXADD,
        USAGE_COMPLEXADD,
        #ifndef __DISABLE_DATABASE
        LEVEL_COMPLEXADD,
        #endif
        complexAdd,
        ARGC_COMPLEXADD,
        BY_USER,
        CHILD
    },
    [ADVCALC_COMPLEXNUMBERSPROD] =
    {
    	CMD_COMPLEXMUL,
        NAME_COMPLEXMUL,
        USAGE_COMPLEXMUL,
        #ifndef __DISABLE_DATABASE
        LEVEL_COMPLEXMUL,
        #endif
        complexMul,
        ARGC_COMPLEXMUL,
        BY_USER,
        CHILD
    },
    [ADVCALC_POLYNOMEVALUATOR] =
    {
    	CMD_POLYNOMEVALUATOR,
    	NAME_POLYNOMEVALUATOR,
    	USAGE_POLYNOMEVALUATOR,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_POLYNOMEVALUATOR,
    	#endif
    	polynomEval,
    	ARGC_POLYNOMEVALUATOR,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_POLYNOMDEVALUATOR] =
    {
    	CMD_POLYNOMDEVALUATOR,
    	NAME_POLYNOMDEVALUATOR,
    	USAGE_POLYNOMDEVALUATOR,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_POLYNOMDEVALUATOR,
    	#endif
    	polynomEval,
    	ARGC_POLYNOMDEVALUATOR,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_GETFORMATTEDDATE] =
    {
    	CMD_GETDATE,
        NAME_GETDATE,
        USAGE_GETDATE,
        #ifndef __DISABLE_DATABASE
        LEVEL_GETDATE,
        #endif
        getDate,
        ARGC_GETDATE,
        BY_USER,
        CHILD
    },
    [ADVCALC_ABSOLUTEORIENTATION] =
    {
    	CMD_ABSOLUTEORIENTATION,
    	NAME_ABSOLUTEORIENTATION,
    	USAGE_ABSOLUTEORIENTATION,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_ABSOLUTEORIENTATION,
    	#endif
    	absoluteOrientation,
    	ARGC_ABSOLUTEORIENTATION,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_ROUTHTABLE] =
    {
    	CMD_ROUTHTABLE,
    	NAME_ROUTHTABLE,
    	USAGE_ROUTHTABLE,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_ROUTHTABLE,
    	#endif
    	routhTable,
    	ARGC_ROUTHTABLE,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_JURYTABLE] =
    {
    	CMD_JURYTABLE,
    	NAME_JURYTABLE,
    	USAGE_JURYTABLE,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_JURYTABLE,
    	#endif
    	juryTable,
    	ARGC_JURYTABLE,
    	BY_USER,
    	CHILD
    },
    [ADVCALC_NEWTONDIFFTABLES] =
    {
    	CMD_NEWTONDIFFTABLES,
        NAME_NEWTONDIFFTABLES,
        USAGE_NEWTONDIFFTABLES,
        #ifndef __DISABLE_DATABASE
        LEVEL_NEWTONDIFFTABLES,
        #endif
        newtonDifferenceTables,
        ARGC_NEWTONDIFFTABLES,
        BY_USER,
        CHILD
    },
    [ADVCALC_LAGRANGEINTERPOLATION] =
    {
    	CMD_LAGRANGEINTERPOLATION,
        NAME_LAGRANGEINTERPOLATION,
        USAGE_LAGRANGEINTERPOLATION,
        #ifndef __DISABLE_DATABASE
        LEVEL_LAGRANGEINTERPOLATION,
        #endif
        lagrangeInterpolation,
        ARGC_LAGRANGEINTERPOLATION,
        BY_USER,
        CHILD
    },
    [ADVCALC_FUNCTIONINTEGRATION] =
    {
    	CMD_FID,
        NAME_FID,
        USAGE_FID,
        #ifndef __DISABLE_DATABASE
        LEVEL_FID,
        #endif
        funcIntegration,
        ARGC_FID,
        BY_USER,
        CHILD
    },
    [ADVCALC_STRAIGHTLINEFITTING] =
    {
    	CMD_STRAIGHTLINEFITTING,
        NAME_STRAIGHTLINEFITTING,
        USAGE_STRAIGHTLINEFITTING,
        #ifndef __DISABLE_DATABASE
        LEVEL_STRAIGHTLINEFITTING,
        #endif
        straightLineFitting,
        ARGC_STRAIGHTLINEFITTING,
        BY_USER,
        CHILD
    },
    [ADVCALC_PARABOLICCURVEFITTING] =
    {
    	CMD_PARABOLICCURVEFITTING,
        NAME_PARABOLICCURVEFITTING,
        USAGE_PARABOLICCURVEFITTING,
        #ifndef __DISABLE_DATABASE
        LEVEL_PARABOLICCURVEFITTING,
        #endif
        parabolicCurveFitting,
        ARGC_PARABOLICCURVEFITTING,
        BY_USER,
        CHILD
    },
    [ADVCALC_LINEARSYSTEMSSOLVER] =
    {
    	CMD_LINEARSYSTEMSSOLVER,
        NAME_LINEARSYSTEMSSOLVER,
        USAGE_LINEARSYSTEMSSOLVER,
        #ifndef __DISABLE_DATABASE
        LEVEL_LINEARSYSTEMSSOLVER,
        #endif
        linearSystemsSolver,
        ARGC_LINEARSYSTEMSSOLVER,
        BY_USER,
        CHILD
    }
};

__MSSHELL_WRAPPER_ static void  secondGradeEquationSolver(const sel_typ argc, char ** argv)
{
    mpfr_t *abc;

    if(argc)
    {
        dim_typ dim[2];
        if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &abc, dim, &dim[COLUMNS])) || dim[ROWS] != 1 || dim[COLUMNS] != 3)
        {
            matrixFree(&abc, dim);
            printUsage(&adv_calc[ADVCALC_SECONDGRADEEQUATIONSOLVER]);
            return;
        }
    }
    else
    {

        msprintf(COLOR_CREDITS, "\nEnter COEFFICIENTS: 'a', 'b' e 'c' of SECOND GRADE Equation:\n\"a*(x^2) + b*x + c\"");
        msprintf(COLOR_CREDITS, "\nby inserting related inline [1 x 3] Matrix.\n\n");
        PRINTHOWTOBACKMESSAGE();

        if(!insertNMMatrix(&abc, (dim_typ2){1, 3}))
            return;
    }

    mpfr_t root[2];
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

    if(difftime)
    	gettimeofday(&tvBegin, NULL);

    if(_secondGradeEquationSolver(abc, root))
    {
    	char buf[MAX_BUFSIZ];
        mpfr_sprintf(buf, "\n1st ROOT = %Rf;\n2nd ROOT = %Rf.\n\n", root[ROOT_X1], root[ROOT_X2]);
        msprintf(COLOR_USER, buf);
    }

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    matrixFree(&abc, ((dim_typ2){1, 3}));
    return;
}

__MSSHELL_WRAPPER_ static void  complexAdd(const sel_typ argc, char ** argv)
{
    mpfr_t *cpx;
    const sel_typ algebra_units = (!access(curLayout)->algebra) ? MAX_COMPLEX_UNITS : exp2(access(curLayout)->algebra);

    if(argc)
    {
        dim_typ dim[2];
		if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &cpx, dim, &dim[COLUMNS])) || dim[ROWS] != 2 || dim[COLUMNS] != algebra_units)
        {
            matrixFree(&cpx, dim);
            printUsage(&adv_calc[ADVCALC_COMPLEXNUMBERSSUM]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter 2x%hu MATRIX whose ROWS contains respectively\nREAL PART and IMAGINARY PART%s of both two Operands.\n\n", algebra_units, algebra_units > MAX_COMPLEX_UNITS ? "s":NULL_CHAR);

        if(!insertNMMatrix(&cpx, (dim_typ2){2, algebra_units}))
            return;
    }

    mpfr_t complexRes[algebra_units];

    static void (* const complexAddFunc[_MAX_ALGEBRA][2])(mpfr_t *restrict, mpfr_t [static MAX_SEDENIONS_UNITS]) =
	{
		{
			_complexAdd,
			_complexSub
		},
		{
			_quaternionsAdd,
			_quaternionsSub
		},
		{
			_octonionsAdd,
			_octonionsSub
		},
		{
			_sedenionsAdd,
			_sedenionsSub
		}
	};

	dim_typ i;
	ityp atime;
	char buf[MAX_BUFSIZ];
	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

	if(difftime)
		gettimeofday(&tvBegin, NULL);

	complexAddFunc[((access(curLayout)->algebra == ALGEBRA_COMPLEXNUMBERS || !access(curLayout)->algebra) ? ALGEBRA_COMPLEXNUMBERS : access(curLayout)->algebra)-1][INVERSE_OPS](cpx, complexRes);

	if(difftime)
		atime = getDiffTime(&tvBegin);

	msprintf(COLOR_USER, "\nRESULT of Operation: (");
    PRINT2N();

    for(i=0; i<algebra_units; ++i)
    {
        mpfr_sprintf(buf, OUTPUT_CONVERSION_FORMAT, *(cpx + i));
        msprintf(COLOR_USER, "%s%s%s", buf, suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") ":" +");
    }

    static char oprchar[2] = "+-";
    printf("%c\n", oprchar[INVERSE_OPS]);

    for(i=0; i<algebra_units; ++i)
    {
        mpfr_sprintf(buf, OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*SECOND_NUMBER) + i));
        msprintf(COLOR_USER, "%s%s%s", buf, suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") is = to:\n":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        mpfr_sprintf(buf, OUTPUT_CONVERSION_FORMAT, complexRes[i]);
        msprintf(COLOR_USER, "%s%s%s", buf, suite_c.algebra_imaginary_units_names[algebra_units][i], i ==algebra_units-1 ? ";\n\n":" + ");
    }

    matrixFree(&cpx, ((dim_typ2){MAX_COMPLEX_UNITS, algebra_units}));

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, atime);
	}

    return;
}

__MSSHELL_WRAPPER_ static void  complexMul(const sel_typ argc, char ** argv)
{
    mpfr_t *cpx;
    const fsel_typ algebra_units = (!access(curLayout)->algebra) ? 2 : exp2(access(curLayout)->algebra);

    if(argc)
    {
        dim_typ dim[2];
		if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &cpx, dim, &dim[COLUMNS])) || dim[ROWS] != 2 || dim[COLUMNS] != 2)
        {
            matrixFree(&cpx, dim);
            printUsage(&adv_calc[ADVCALC_COMPLEXNUMBERSPROD]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter 2x%hu MATRIX whose ROWS contains respectively\nREAL PART and IMAGINARY PART%s of both two Operands.\n\n", algebra_units, algebra_units > MAX_COMPLEX_UNITS ? "s":NULL_CHAR);

        if(!insertNMMatrix(&cpx, (dim_typ2){2, algebra_units}))
            return;
    }

    mpfr_t complexRes[algebra_units];

    static void (* const complexMulFunc[_MAX_ALGEBRA][2])(mpfr_t *restrict, mpfr_t [static MAX_SEDENIONS_UNITS]) =
	{
		{
			_complexMul,
			_complexDiv
		},
		{
			_quaternionsMul,
			_quaternionsDiv
		},
		{
			_octonionsMul,
			_octonionsDiv
		},
		{
			_sedenionsMul,
			_sedenionsDiv
		}
	};

	dim_typ i;
	ityp atime; 
	char buf[MAX_BUFSIZ];
	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

	if(difftime)
		gettimeofday(&tvBegin, NULL);

	complexMulFunc[((access(curLayout)->algebra == ALGEBRA_COMPLEXNUMBERS || !access(curLayout)->algebra) ? ALGEBRA_COMPLEXNUMBERS : access(curLayout)->algebra)-1][INVERSE_OPS](cpx, complexRes);
    if(difftime)
    	atime = getDiffTime(&tvBegin);

	msprintf(COLOR_USER, "\nRESULT of Operation: (");
    PRINT2N();

    for(i=0; i<algebra_units; ++i)
    {
        mpfr_sprintf(buf, OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*FIRST_NUMBER) + i));
        msprintf(COLOR_USER, "%s%s%s", buf, suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") ":" +");
    }

    static char oprchar[2] = "*/";
    printf("%c\n", oprchar[INVERSE_OPS]);

    for(i=0; i<algebra_units; ++i)
    {
        mpfr_sprintf(buf, OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*SECOND_NUMBER) + i));
        msprintf(COLOR_USER, "%s%s%s", buf, suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") is = to:\n":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        mpfr_sprintf(buf, OUTPUT_CONVERSION_FORMAT, complexRes[i]);
        msprintf(COLOR_USER, "%s%s%s", buf, suite_c.algebra_imaginary_units_names[algebra_units][i], i ==algebra_units-1 ? ";\n\n":" + ");
    }

    matrixFree(&cpx, ((dim_typ2){MAX_COMPLEX_UNITS, algebra_units}));

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, atime);
	}
    return;
}

__MSSHELL_WRAPPER_ static void  polynomEval(const sel_typ argc, char ** argv)
{
	mpfr_t *table;
	dim_typ dim[2];

	if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &table, dim, &dim[COLUMNS])))
        {
            matrixFree(&table, dim);
            printUsage(&adv_calc[ADVCALC_POLYNOMEVALUATOR]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter the Polynom n-dimensioned Row-Matrix.\n\n");
        if(!insertMatrix(table, dim[ROWS], dim[COLUMNS], false))
            return;
    }

    mpfr_t scal;
    PRINTHOWTOBACKMESSAGE();
    PRINTN();

    if(argc > 1)
    {
        if(!parse(argv[1], &scal))
        {
        	mpfr_clear(scal);
        	matrixFree(&table, dim);
            return;
        }
    }
    else
	{
		while(requires(scal, NULL, "Enter a double floating-point Scalar Number.\n", "Inserted Scalar", PARSER_NOSETTINGS) && isNullVal(scal))
        {
        	mpfr_clear(scal);
            CLEARBUFFER();
            if(exitHandleCheck)
			{
				matrixFree(&table, dim);
				return;
			}
		}
	}

    // dim_typ times = 1;
    struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	dim_typ times = __pmode__ == ADVCALC_POLYNOMDEVALUATOR;
	//	const bool dermode = __pmode__ == ADVCALC_POLYNOMDEVALUATOR;

	if(difftime)
		gettimeofday(&tvBegin, NULL);

	mpfr_t result, tmp;
	char der_order[MINMIN_STRING] = NULL_CHAR;
	
	mpfr_init(tmp);

	if(times)
	{
		if(argc > 2)
		{
            if(!parse(argv[1], &tmp))
            {
            	mpfr_clear(tmp);
            	matrixFree(&table, dim);
                return;
            }
            times = mpfr_get_ui(tmp, MPFR_RNDN);
	    }
	    else
		{

			while(requires(tmp, NULL, "Enter a non-zero positive integer representative of the Derivative Order.\n", "Inserted Derivative Order:", PARSER_NOSETTINGS) && isNullVal(tmp))
	        {
	        	mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
					return;
	            printErr(5, "Invalid inserted Value.\nMust be a non-zero integer between 1 and %hu", dim[COLUMNS]);

			}
		}

		mpfr_init(result);
	    for(dim_typ i=0; i<times; ++i)
			deval(result, table, dim[COLUMNS], scal);
		msprintf(COLOR_USER, "The Derivative of the inserted POLYNOM is the POLYNOM: \n");
		printMatrix(stdout, &table, dim);
		sprintf(der_order, "%hu-Derivative ", times);
	}
	else
	{
		mpfr_init(result);
        eval(result, table, dim[COLUMNS], scal);
    }

	char buf[MAX_BUFSIZ];
	mpfr_sprintf(buf, "\nInserted POLYNOM %sEvaluated in: %Rf RESULT is: %Rf.\n\n", der_order, result);
	msprintf(COLOR_USER, buf);

	if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

	mpfr_clears(scal, result, tmp, NULL); 
    matrixFree(&table, dim);
	return;
}

__MSSHELL_WRAPPER_ static void  getDate(const sel_typ argc, char ** argv)
{
    mpfr_t tmp, mpftmp, mpftmp2;
    ityp tmp2, tmp3;

    sel_typ dd, mm;
    sel_typ yyyy;

    if(argc)
    {
        if(argc != 3)
        {
            printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
            return;
        }
		if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (dd = mpfr_get_ui(tmp, MPFR_RNDN))))
        {
        	mpfr_clear(tmp);
            printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
            return;
        }
        if(argv[1][strlen(argv[1])-1] != TERMINATING_CHAR)
            strcat(argv[1], TERMINATING_STRING);
        if((!parse(argv[1], &tmp)) || mpfr_cmp_ui(tmp, (mm = mpfr_get_ui(tmp, MPFR_RNDN))))
        {
        	mpfr_clear(tmp);
            printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
            return;
        }
        if(argv[2][strlen(argv[2])-1] != TERMINATING_CHAR)
            strcat(argv[2], TERMINATING_STRING);
        if((!parse(argv[2], &tmp)) || mpfr_cmp_ui(tmp, (yyyy = mpfr_get_ui(tmp, MPFR_RNDN))))
        {
        	mpfr_clear(tmp);
            printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
            return;
        }

    }
    else
    {
        char seperator[SIGN_STRING] = NULL_CHAR;

        strcpy(seperator, "]\n[");

        msprintf(COLOR_CREDITS, "\nEnter Numeric DATE as expected format:\n");
        msprintf(COLOR_CREDITS, "[DD%sMM%sYYYY]\n", seperator, seperator);
        PRINTHOWTOBACKMESSAGE();

        while((requires(tmp, NULL, NULL_CHAR, "Inserted DAY is", PARSER_SHOWRESULT) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (dd = mpfr_get_ui(tmp, MPFR_RNDN))) || dd < 1 || dd > MAX_MONTH_DAYS) ||
     (requires(mpftmp ,NULL, NULL_CHAR, "Inserted MONTH is", PARSER_SHOWRESULT) || isNullVal(mpftmp) || mpfr_cmp_ui(mpftmp, (mm = mpfr_get_ui(mpftmp, MPFR_RNDN))) || mm < 1 || mm > MAX_MONTHS) ||
     (requires(mpftmp2, NULL, NULL_CHAR, "Inserted YEAR is", PARSER_SHOWRESULT) || isNullVal(mpftmp2) || mpfr_cmp_ui(mpftmp2, (yyyy = mpfr_get_ui(mpftmp2, MPFR_RNDN))) || yyyy < 1))
        {
			mpfr_clears(tmp, mpftmp, mpftmp2, NULL); 
            CLEARBUFFER();
            if(exitHandleCheck)
				return;
            printErr(33, "Invalid [DD MM YYYY] format.\nYou have to insert DAY, MONTH and YEAR as non-negative and non-zero integers.\n\
and the first two Params must respectively be less than %hu and %hu", MAX_MONTH_DAYS, MAX_MONTHS);
        }
    }

    msprintf(COLOR_USER, "\nFormatted DATE correspondent to the Inserted Numeric ONE:\n%hu/%hu/%llu is: %s %hu %s %llu.\n\n",
           tmp, mpftmp, mpftmp2, getDayName(getDayNumber(dd, mm, yyyy)), dd, getMonthName(mm), yyyy);
	mpfr_clears(tmp, mpftmp, mpftmp2, NULL); 
    return;
}

__MSSHELL_WRAPPER_ static void  absoluteOrientation(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix = NULL;
    mpfr_t *matrix2 = NULL;
    dim_typ dim[2];
    dim_typ dim2[2];

    if(argc)
    {
    	if(argc == 2)
    	{
    		if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            	strcat(argv[0], TERMINATING_STRING);
	        if(!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS]))
	        {
	        	matrixFree(&matrix, dim);
	        	printUsage(&adv_calc[ADVCALC_ABSOLUTEORIENTATION]);
	            return;
	    	}
	    	if(argv[1][strlen(argv[1])-1] != TERMINATING_CHAR)
            	strcat(argv[1], TERMINATING_STRING);
	        if((!matrixToken(argv[1], &matrix, dim2, &dim2[COLUMNS])) || dim[ROWS] != dim2[ROWS] || dim[COLUMNS] != dim2[COLUMNS])
	        {
	        	matrixFree(&matrix, dim);
	        	matrixFree(&matrix2, dim2);
	        	printUsage(&adv_calc[ADVCALC_ABSOLUTEORIENTATION]);
	            return;   
	    	}
		}
		else
		{
			printUsage(&adv_calc[ADVCALC_ABSOLUTEORIENTATION]);
	        return;
		}
    }
    else
	{
		
		if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], false))
			return;
		
		if(!insertNMMatrix(&matrix2, dim))
		{
			matrixFree(&matrix, dim);
			return;
		}
	}
	
	mpfr_t *result;
	//const register dim_typ rdim = dim[COLUMNS];
    msprintf(COLOR_SYSTEM, "\nThe TRANSFORMATION MATRIX R of Inserted Point Matrices is:\n\n");
    
    if(!matrixAlloc(&result, (dim_typ2){dim[ROWS], dim[ROWS]}))
	{
		matrixFree(&matrix, dim);
		matrixFree(&matrix2, dim2);
        return;
    }

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	if(absor(dim, matrix, matrix2, &result) == ABSOR_ERROR)
	{
		matrixFree(&matrix, dim);
		matrixFree(&matrix2, dim2);
		matrixFree(&result, ((dim_typ2){dim[ROWS], dim[ROWS]}));
		printErr(33, "Something got wrong during Absolute Orientation Transformation");
		return;
	}
    
    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}
    
    printMatrix(stdout, &result, ((dim_typ2){dim[ROWS], dim[ROWS]}));
    matrixFree(&matrix, dim);
    matrixFree(&matrix2, dim2);
    matrixFree(&result, ((dim_typ2){dim[ROWS], dim[ROWS]}));
    return;
}

__MSSHELL_WRAPPER_ static void  routhTable(const sel_typ argc, char ** argv)
{
	mpfr_t *table;
	dim_typ dim[2];

	if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &table, dim, &dim[COLUMNS])) || dim[COLUMNS] <= 2)
        {
            matrixFree(&table, dim);
            printUsage(&adv_calc[ADVCALC_ROUTHTABLE]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter the Polynom n>2 dimensioned Row-Matrix.\n\n");
        if(!insertMatrix(table, dim[ROWS], dim[COLUMNS], false))
            return;
        if(dim[COLUMNS] <= 2)
        {
        	matrixFree(&table, dim);
            printUsage(&adv_calc[ADVCALC_ROUTHTABLE]);
            return;
        }
    }

    short permanences;
    struct timeval tvBegin;
    fsel_typ nullrow = 0; // could not be null-row-ed the first row
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

	if(difftime)
		gettimeofday(&tvBegin, NULL);

    if((permanences = _routhTable(&table, dim[COLUMNS], &nullrow)) == ROUTHTABLE_ALLOC_ERROR)
    	printErr(12, "Routh Table Evaluator Dynamic Memory Allocation Problem");
    else
    {
    	msprintf(COLOR_USER, "The ROUTH TABLE of the inserted Polynom Matrix is: \n");
    	printMatrix(stdout, &table, ((dim_typ2){dim[COLUMNS], ((dim_typ)((dim[COLUMNS]*0.5) + 1))}));
    	msprintf(COLOR_USER, "PERMANENCES: %hu, VARIATIONS: %hu", permanences, dim[COLUMNS]-1-permanences);
    	if(nullrow)
    		msprintf(COLOR_USER, "\nIt has been used the AUXILIARY POLYNOM's Derivative on the %huth NULL ROW.\n\n", nullrow+1);
    	if(difftime)
	    {
	    	PRINTL();
	        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}
    }

    matrixFree(&table, dim);
	return;
}

__MSSHELL_WRAPPER_ static void  juryTable(const sel_typ argc, char ** argv)
{
	mpfr_t *table;
	dim_typ dim[2];

	if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &table, dim, &dim[COLUMNS])) || dim[COLUMNS] <= 2)
        {
            matrixFree(&table, dim);
            printUsage(&adv_calc[ADVCALC_JURYTABLE]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter the Polynom n>2 dimensioned Row-Matrix.\n\n");
        if(!insertMatrix(table, dim[ROWS], dim[COLUMNS], false))
            return;
        if(dim[COLUMNS] <= 2)
        {
        	matrixFree(&table, dim);
            printUsage(&adv_calc[ADVCALC_JURYTABLE]);
            return;
        }
    }

    struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

	if(difftime)
		gettimeofday(&tvBegin, NULL);

	sel_typ result;

    if((result = _juryTable(&table, dim[COLUMNS])) == JURYTABLE_ALLOC_ERROR)
    	printErr(12, "Jury Table Evaluator Dynamic Memory Allocation Problem");
    else
    {
    	msprintf(COLOR_USER, "\nThe JURY TABLE of the inserted Polynom Matrix is: \n");
    	printMatrix(stdout, &table, ((dim_typ2){((dim[COLUMNS]-1)<<1)-3,dim[COLUMNS]}));
    	if(difftime)
	    {
	    	PRINTL();
	        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}
		msprintf(COLOR_USER, "JURY Criterion is: %s.\n", result == JURYTABLE_SATISFIED?"SATISFIED":"NOT SATISFIED");
	}

    matrixFree(&table, dim);
	return;
}


__MSSHELL_WRAPPER_ static void  newtonDifferenceTables(const sel_typ argc, char ** argv)
{

    mpfr_t x[access(curLayout)->max_newton_difftables_dim];
    mpfr_t y[access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim];

    dim_typ n;
    mpfr_t tmp;

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (n = mpfr_get_ui(tmp, MPFR_RNDN))) || n < access(curLayout)->min_newton_difftables_dim || n > access(curLayout)->max_newton_difftables_dim)
        {
        	mpfr_clear(tmp);
            printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter Difference Table DIMENSION.\n");
        PRINTHOWTOBACKMESSAGE();
        
        while(requires(tmp, NULL, NULL_CHAR, "Inserted Difference Table DIMENSION is:", PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (n = mpfr_get_ui(tmp, MPFR_RNDN))) || n < access(curLayout)->min_newton_difftables_dim || n > access(curLayout)->max_newton_difftables_dim)
        {
        	mpfr_clear(tmp);
            CLEARBUFFER();
            if(exitHandleCheck)
				return;
            printErr(5, "Invalid inserted Value.\nMust be a non-negative integer between %hu and %hu", access(curLayout)->min_newton_difftables_dim, access(curLayout)->max_newton_difftables_dim);
        }
    }

    dim_typ i, j;

    if(argc > 1)
        if(argc == n+1)
        {
            char *token = NULL;
			for(i=0; i<n; ++i)
            {

                if((token = strtok(argv[i+1], ",")))
                {
                	if(token[strlen(token)-1] != TERMINATING_CHAR)
            			strcat(token, TERMINATING_STRING);
                    if(!parse(token, &x[i]))
                    {
                    	for(j=0; j<i; ++j)
                    		mpfr_clear(x[j]);
                        printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                        return;
                    }
                }
                else
                {
                    printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                    return;
                }

                if((token = strtok(NULL, ",")))
                {
                	if(token[strlen(token)-1] != TERMINATING_CHAR)
            			strcat(token, TERMINATING_STRING);
                    if(!parse(token, &y[i][0]))
                    {
                    	for(dim_typ j=0; j<i; ++j)
                    		mpfr_clear(y[j][0]);
                        printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                        return;
                    }
                }
                else
                {
                    printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                    return;
                }
            }
        }
        else
        {
            printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
            return;
        }
    else
        for(i=0; i<n; ++i)
        {
            msprintf(COLOR_CREDITS, "Enter couple No. %hu as expected format:\n[X]\n[Y]\n", i);
            /// scanf("%f %f",&x[i],&y[i][0]); // PAY STRICT ATTENTION TO THIS HANDLE

            while(requires(x[i], NULL, NULL_CHAR, "Inserted X is:", PARSER_SHOWRESULT) || isNullVal(x[i]) || requires(y[i][0], NULL, NULL_CHAR, "Inserted Y is:", PARSER_SHOWRESULT) || isNullVal(y[i][0])) 
            {
				mpfr_clear(tmp);
                CLEARBUFFER();
                if(exitHandleCheck)
					return;
                printErr(5, "Invalid inserted Value");
            }
        }

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

	if(difftime)
		gettimeofday(&tvBegin, NULL);

    newtonDifferenceTable(n, y, FORWARD_DIFFTAB);
    showNewtonDifferenceTable(n, x, y, FORWARD_DIFFTAB);
    newtonDifferenceTable(n, y, BACKWARD_DIFFTAB);
    showNewtonDifferenceTable(n, x, y, BACKWARD_DIFFTAB);
    mpfr_clear(tmp);
    
    for(i=0; i<access(curLayout)->max_newton_difftables_dim; ++i)
    {
    	mpfr_clear(x[i]);
    	for(j=0; j<access(curLayout)->max_newton_difftables_dim; ++j)
    		mpfr_clear(y[i][j]);
    }

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    return;
}

__MSSHELL_WRAPPER_ static void  lagrangeInterpolation(const sel_typ argc, char ** argv)
{
    mpfr_t tmp;
    dim_typ dim;

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
        	strcat(argv[0], TERMINATING_STRING);
        if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (dim = mpfr_get_ui(tmp, MPFR_RNDN))) || dim < 1 || dim > USHRT_MAX)
        {
        	mpfr_clear(tmp);
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter Data DIMENSION.\n\n");
        while(requires(tmp, NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT) && isNullVal(tmp) || mpfr_cmp_ui(tmp, (dim = mpfr_get_ui(tmp, MPFR_RNDN))) || dim < 1 || dim > USHRT_MAX)
        {
        	mpfr_clear(tmp);
            CLEARBUFFER();
            if(exitHandleCheck)
				return;
            printErr(33, "Invalid inserted Data DIMENSION.\nMust be an integer between 1 and %z", USHRT_MAX);
        }
    }

    mpfr_t *xy = NULL;
    dim_typ i, j;

    if(argc > 1)
    {
        dim_typ rc[2];
		if(argv[1][strlen(argv[1])-1] != TERMINATING_CHAR)
        	strcat(argv[1], TERMINATING_STRING);
        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != 2 || rc[COLUMNS] != dim)
        {
            matrixFree(&xy, rc);
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter related Matrix filled with Data you want Interpolation to Process,\n");
        msprintf(COLOR_CREDITS, "by putting on each rows the %hu X and Y Values.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){2, dim}))
            return;
    }

    mpfr_t xp;

    if(argc > 2)
	{
		if(argv[2][strlen(argv[2])-1] != TERMINATING_CHAR)
        	strcat(argv[2], TERMINATING_STRING);
        if(!parse(argv[2], &xp))
        {
        	matrixFree(&xy, ((dim_typ2){2, dim}));
        	mpfr_clear(xp);
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
	}
    else
    {
        msprintf(COLOR_CREDITS, "Enter X VALUE to find Y one.\n");
        while(requires(xp, NULL, NULL_CHAR, "X VALUE | Y = F(X) is:", PARSER_SHOWRESULT) || isNullVal(xp))
        {
        	mpfr_clear(xp);
            CLEARBUFFER();
            if(exitHandleCheck)
            {
            	matrixFree(&xy, ((dim_typ2){2, dim}));
                return;
            }
            printErr(5, "Invalid inserted Value");
        }
    }

    mpfr_t yp;
    mpfr_t dr, nr;

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	mpfr_inits(nr, dr, NULL); 
	mpfr_init_set_ui(yp, 0, MPFR_RNDN);

	for(i = 0; i < dim; ++i)
    {
		mpfr_set_ui(dr, 1, MPFR_RNDN);
		mpfr_set(nr, dr, MPFR_RNDN);
        for(j = 0; j<dim; ++j)
            if(i!=j)
            {
            	mpfr_sub(tmp, xp, *(xy + (dim*XROW) + j), MPFR_RNDN);
            	mpfr_mul(nr, nr, tmp, MPFR_RNDN);
            	mpfr_sub(tmp, *(xy + (dim*XROW) + i), *(xy + (dim*XROW) + j), MPFR_RNDN);
            	mpfr_mul(dr, dr, tmp, MPFR_RNDN);
            }

		mpfr_mul(tmp, dr, *(xy + (dim*YROW) + i), MPFR_RNDN);
		mpfr_div(tmp, nr, tmp, MPFR_RNDN);
		mpfr_add(yp, yp, tmp, MPFR_RNDN);
    }

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

	char buf[MAX_BUFSIZ];
    matrixFree(&xy, ((dim_typ2){2, dim}));
    mpfr_sprintf(buf, "\nRequested Y VALUE is: %Rf.\n\n", yp);
    msprintf(COLOR_USER, buf);
    mpfr_clears(tmp, nr, dr, xp, yp, NULL); 
    return;
}

__MSSHELL_WRAPPER_ static void  funcIntegration(const sel_typ argc, char ** argv)
{
    mpfr_t x0, xn;
    mpfr_t h, s;

    dim_typ funcID;
    dim_typ j;

    funcID = selectListItem(MAX_FIDS, MAX_FIDS > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
                            "Select desired Function you want to Integrate", ext_math.funcnames);

    if(funcID == MAX_FIDS) return;

    mpfr_t tmp;
    sel_typ mode;

    if(argc)
    {
        if((mode = strtod(argv[0], NULL)) == 3) return;
        if(mode < 0 || mode > 3)
        {
            printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nSelect Defined Integration Calculus Mode:\n");
        msprintf(COLOR_CREDITS, "- A for Simpsons' 1/3 method,\n- B for Simpsons' 3/8 method,\n- C per Trapezoidal Method;\n");
        msprintf(COLOR_CREDITS, "- D to go Back...\n\n");

        do if((mode = toupper(getch())) == 'D') return;
        while(mode < 'A' && mode > 'D');

        mode -= 'A';
    }

    dim_typ i, n;

    if(argc > 1)
        if(argc == 4)
        {
        	if(argv[1][strlen(argv[1])-1] != TERMINATING_CHAR)
        		strcat(argv[1], TERMINATING_STRING);
            if(!parse(argv[1], &x0))
            {
            	mpfr_clear(x0);
                printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                return;
            }
			if(argv[2][strlen(argv[2])-1] != TERMINATING_CHAR)
        		strcat(argv[2], TERMINATING_STRING);
            if(!parse(argv[2], &xn))
            {
            	mpfr_clears(x0, xn, NULL); 
                printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                return;
            }
			if(argv[3][strlen(argv[3])-1] != TERMINATING_CHAR)
        		strcat(argv[3], TERMINATING_STRING);
            if((!parse(argv[3], &tmp)) || mpfr_cmp_ui(tmp, (n = mpfr_get_ui(tmp, MPFR_RNDN))))
            {
            	mpfr_clears(x0, xn, tmp, NULL); 
                printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
                return;
            }
        }
        else
        {
            printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
            return;
        }
    else
    {
        char seperator[SIGN_STRING] = NULL_CHAR;

        strcpy(seperator, "]\n[");
        msprintf(COLOR_CREDITS, "\nEnter INTEGRATION Extremes and Intervals NUMBER as expected format:\n");
        msprintf(COLOR_CREDITS, "[x0%sxN%sNo]\n", seperator, seperator);

        while(requires(x0, NULL, NULL_CHAR, "First INTEGRATION Extreme is:", PARSER_SHOWRESULT) || isNullVal(x0) ||
				requires(xn, NULL, NULL_CHAR, "Second INTEGRATION Extreme is:", PARSER_SHOWRESULT) || isNullVal(xn) ||
				requires(tmp, NULL, NULL_CHAR, "Number of inserted INTERVALS is:", PARSER_SHOWRESULT) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (n = mpfr_get_ui(tmp, MPFR_RNDN))) || n < 1 || n > INT_MAX)
        {
			mpfr_clears(x0, xn, tmp, NULL); 
            CLEARBUFFER();
            if(exitHandleCheck)
				return;
            printErr(33, "Invalid inserted Intervals NUMBER.\nMust be an integer between 1 and %z", INT_MAX);
        }
    }
    
    mpfr_sub(h, xn, x0, MPFR_RNDN);
    mpfr_div_ui(h, h, n, MPFR_RNDN);

 	mpfr_t result, tmp2;
    void (* const y)(mpfr_t, register mpfr_t) = ext_math.functions[funcID];
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

    if(difftime)
    	gettimeofday(&tvBegin, NULL);
    	
    mpfr_inits(result, tmp2, NULL); 

    switch(mode)
    {
        case SIMPSON1DIV8_RULE:
        	mpfr_add(s, x0, h, MPFR_RNDN); 
        	y(s, s);
        	mpfr_mul_ui(s, s, 4, MPFR_RNDN);
        	y(tmp, xn);
        	mpfr_add(s, s, tmp, MPFR_RNDN);
        	y(tmp, x0);
        	mpfr_add(s, s, tmp, MPFR_RNDN);

            for(i = 3; i<=n-1; i+=2)
        	{
				mpfr_mul_ui(tmp2, h, i-1, MPFR_RNDN);
				mpfr_add(tmp2, tmp2, x0, MPFR_RNDN);
				y(tmp2, tmp2);
				mpfr_mul_ui(tmp2, tmp2, 2, MPFR_RNDN);
				mpfr_mul_ui(tmp, h, i, MPFR_RNDN);
				mpfr_add(tmp, tmp, x0, MPFR_RNDN);
				y(tmp, tmp);
				mpfr_mul_ui(tmp, tmp, 4, MPFR_RNDN);
				mpfr_add(s, tmp2, tmp, MPFR_RNDN);
			}

			mpfr_div_ui(result, h, 3, MPFR_RNDN);
			mpfr_mul(result, result, s, MPFR_RNDN);
            break;

        case SIMPSON3DIV8_RULE:
        {
            bool flag;

			y(s, xn);
			y(tmp, x0);
			mpfr_add(s, s, tmp, MPFR_RNDN);

            for(i = 1; i<=n-1;++i)
            {
                for(j=1;j<=n-1;++j)
                    if((flag = i == 3*j))
                        break;
            	mpfr_mul_ui(tmp, h, i, MPFR_RNDN);
            	mpfr_add(tmp, tmp, x0, MPFR_RNDN);
            	mpfr_mul_ui(tmp, tmp, 3-flag, MPFR_RNDN);
				mpfr_add(s, s, tmp, MPFR_RNDN);
            }

			mpfr_div_ui(result, h, 8, MPFR_RNDN);
			mpfr_mul_ui(result, result, 3, MPFR_RNDN);
			mpfr_mul(result, result, s, MPFR_RNDN);
            break;
        }

        case TRAPEZOIDAL_RULE:
        	y(s, xn);
        	y(tmp, x0);
			mpfr_add(s, s, tmp, MPFR_RNDN);

            for(i = 0; ++i < n; )
            {
            	mpfr_mul_ui(tmp, h, i, MPFR_RNDN);
            	mpfr_add(tmp, tmp, x0, MPFR_RNDN);
            	y(tmp, tmp);
            	mpfr_mul_ui(tmp, tmp, 2, MPFR_RNDN);
            	mpfr_add(s, s, tmp, MPFR_RNDN);
            }
		
			mpfr_div_ui(result, h, 2, MPFR_RNDN);
			mpfr_mul(result, result, s, MPFR_RNDN);
            break;
    }

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

	char buf[MAX_BUFSIZ];
    mpfr_sprintf(buf, "%s(x) Function Integral Value calculated between: x0 = %Rf and xN = %Rf,\nwith %hu Intervals NUMBER is: %Rf\n\n", ext_math.funcnames[funcID], x0, xn, n, result);
    msprintf(COLOR_USER, buf);
   	mpfr_clears(result, x0, xn, tmp, tmp2, NULL); 
    return;
}

__MSSHELL_WRAPPER_ static void  straightLineFitting(const sel_typ argc, char ** argv)
{
    mpfr_t tmp;
    dim_typ dim;

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
        	strcat(argv[0], TERMINATING_STRING);
        if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (dim = mpfr_get_ui(tmp, MPFR_RNDN))))
        {
        	mpfr_clear(tmp);
            printUsage(&adv_calc[ADVCALC_STRAIGHTLINEFITTING]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter Data DIMENSION.\n\n");
        while(requires(tmp, NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (dim = mpfr_get_ui(tmp, MPFR_RNDN))) || dim < 1 || dim > USHRT_MAX)
        {
        	mpfr_clear(tmp);
            CLEARBUFFER();
            if(exitHandleCheck)
				return;
            printErr(33, "Invalid inserted Data DIMENSION.\nMust be an integer between 1 and %zu", USHRT_MAX);
        }
    }

    mpfr_t *xy;
    dim_typ i;

    // we must seek for the BACKTRACKING FEATURE

    if(argc > 1)
    {
        dim_typ rc[2];
		if(argv[1][strlen(argv[1])-1] != TERMINATING_CHAR)
        	strcat(argv[1], TERMINATING_STRING);
        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != 2 || rc[COLUMNS] != dim)
        {
            matrixFree(&xy, rc);
       	 	mpfr_clear(tmp);
            printUsage(&adv_calc[ADVCALC_STRAIGHTLINEFITTING]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter the Matrix filled with EQUATIONS DATA,\n");
        msprintf(COLOR_CREDITS, "by putting on each ROWS the %hu X and Y VALUES.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){2, dim}))
        {
        	mpfr_clear(tmp);
            return;
        }
    }

    mpfr_t sum_x, sum_xy, sum_x2, sum_y;
    mpfr_init_set_ui(sum_x, 0, MPFR_RNDN);
    mpfr_init_set_ui(sum_xy, 0, MPFR_RNDN);
    mpfr_init_set_ui(sum_x2, 0, MPFR_RNDN);
    mpfr_init_set_ui(sum_y, 0, MPFR_RNDN);

    const register dim_typ cache[2] =
    {
    	dim*XROW,
    	dim*YROW
    };

    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

    if(difftime)
    	gettimeofday(&tvBegin, NULL);

	for(i = 0; i < dim; ++i)
    {
		mpfr_add(sum_x, sum_x, *(xy + cache[XROW] + i), MPFR_RNDN); 
		mpfr_add(sum_y, sum_y, *(xy + cache[YROW] + i), MPFR_RNDN);
		mpfr_mul(tmp, *(xy + cache[XROW] + i), *(xy + cache[YROW] + i), MPFR_RNDN);
		mpfr_add(sum_xy, sum_xy, tmp, MPFR_RNDN);
		mpfr_pow_ui(tmp, *(xy + cache[XROW] + i), 2, MPFR_RNDN);
		mpfr_add(sum_x2, sum_x2, tmp, MPFR_RNDN);
    }

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}
	
	mpfr_t b, a;
	char buf[MAX_BUFSIZ];
	
    matrixFree(&xy, ((dim_typ2){2, dim}));
    mpfr_inits(b, a, NULL); 
    mpfr_mul_ui(b, sum_xy, dim, MPFR_RNDN);
    mpfr_mul(tmp, sum_x, sum_y, MPFR_RNDN);
    mpfr_sub(b, b, tmp, MPFR_RNDN);
    mpfr_pow_ui(tmp, sum_x, 2, MPFR_RNDN);
    mpfr_mul_ui(a, sum_x2, dim, MPFR_RNDN);
    mpfr_sub(a, a, tmp, MPFR_RNDN);
	mpfr_div(b, b, a, MPFR_RNDN);
	mpfr_mul(a, b, sum_x, MPFR_RNDN); 
	mpfr_div_ui(a, a, dim, MPFR_RNDN);
    mpfr_sprintf(buf, "\na = %Rf; b = %Rf.\nThe Equation is %Rf + %sX.\n\n", a, b, a, b);
    msprintf(COLOR_USER, buf);
	mpfr_clears(b, a, tmp, sum_x, sum_xy, sum_x2, sum_y, NULL); 
    return;
}

__MSSHELL_WRAPPER_ static void  parabolicCurveFitting(const sel_typ argc, char ** argv)
{
    mpfr_t tmp;
    dim_typ dim;

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
        	strcat(argv[0], TERMINATING_STRING);
        if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (dim = mpfr_get_ui(tmp, MPFR_RNDN))))
        {
        	mpfr_clear(tmp);
            printUsage(&adv_calc[ADVCALC_PARABOLICCURVEFITTING]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter Data DIMENSION.\n\n");
        while(requires(tmp, NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (dim = mpfr_get_ui(tmp, MPFR_RNDN))) || dim < 1 || dim > USHRT_MAX)
        {
    		mpfr_clear(tmp);
            CLEARBUFFER();
            if(exitHandleCheck)
				return; 
            printErr(33, "Invalid inserted Data DIMENSIONM.\nMust be an integer between 1 and %z", USHRT_MAX);
        }
    }

    mpfr_t *xy;
    // we must seek for the BACKTRACKING FEATURE

    if(argc > 1)
    {
        dim_typ rc[2];
		if(argv[1][strlen(argv[1])-1] != TERMINATING_CHAR)
        	strcat(argv[1], TERMINATING_STRING);
        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != 2 || rc[COLUMNS] != dim)
        {
        	mpfr_clear(tmp); 
            matrixFree(&xy, rc);
            printUsage(&adv_calc[ADVCALC_PARABOLICCURVEFITTING]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter the Matrix filled with EQUATIONS DATA,\n");
        msprintf(COLOR_CREDITS, "by putting on each Rows the %hu X and Y VALUES.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){2, dim}))
        {
        	mpfr_clear(tmp);
            return;
        }
    }

    mpfr_t *matrix;

    if(!matrixAlloc(&matrix, (dim_typ2){3, 4}))
    {
    	mpfr_clear(tmp);
        matrixFree(&xy, ((dim_typ2){2, dim}));
        return;
    }

	mpfr_set_ui(*matrix, dim, MPFR_RNDN);
    dim_typ i;

    const register dim_typ cache[2] =
    {
    	dim*XROW,
    	dim*YROW
    };

    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

    if(difftime)
    	gettimeofday(&tvBegin, NULL);

    #pragma omp parallel for
	for(i = 0; i < dim; ++i)
    {
    	mpfr_add(*(matrix + 4), *(matrix + 4), *(xy + cache[XROW] + i), MPFR_RNDN); 
    	mpfr_set(*(matrix + 1), *(matrix + 4), MPFR_RNDN);
		mpfr_add(*(matrix + 3), *(matrix + 3), *(xy + cache[YROW] + i), MPFR_RNDN); 
		mpfr_pow_ui(tmp, *(xy + cache[XROW] + i), 2, MPFR_RNDN);
		mpfr_add(*(matrix + 8), *(matrix + 8), tmp, MPFR_RNDN);
		mpfr_set(*(matrix + 2), *(matrix + 8), MPFR_RNDN); 
		mpfr_pow_ui(tmp, *(xy + cache[XROW] + i), 3, MPFR_RNDN);
		mpfr_add(*(matrix + 9), *(matrix + 9), tmp, MPFR_RNDN);
		mpfr_set(*(matrix + 6), *(matrix + 9), MPFR_RNDN); 
		mpfr_pow_ui(tmp, *(xy + cache[XROW] + i), 4, MPFR_RNDN);
		mpfr_add(*(matrix + 10), *(matrix + 10), tmp, MPFR_RNDN);
		mpfr_mul(tmp, *(xy + cache[XROW] + i), *(xy + cache[YROW] * i), MPFR_RNDN);
		mpfr_add(*(matrix + 7), *(matrix + 7), tmp, MPFR_RNDN); 
        mpfr_pow_ui(tmp, *(xy + cache[XROW] + i), 2, MPFR_RNDN);
        mpfr_mul(tmp, tmp, *(xy + cache[YROW]*i), MPFR_RNDN);
        mpfr_add(*(matrix + 11), *(matrix + 11), tmp, MPFR_RNDN); 
    }

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    dim_typ j, k;
    mpfr_t ratio;
    char buf[MAX_BUFSIZ];
    
    mpfr_init(ratio);
    matrixFree(&xy, ((dim_typ2){2, dim}));

    #pragma omp parallel for
	for(i = 0; i < 3; ++i)
        #pragma omp parallel for
		for(j = 0; j < 3; ++j)
            if(i!=j)
            {
            	mpfr_div(ratio, *(matrix + (j<<2) + i), *(matrix + (i<<2) + i), MPFR_RNDN); 
                for(k = 0; k < 4; ++k)
                {
					mpfr_mul(ratio, ratio, *(matrix + (i<<2) + k), MPFR_RNDN); 
					mpfr_sub(*(matrix + (j<<2) + k), *(matrix + (j<<2) + k), ratio, MPFR_RNDN); 
				}
            }

    #pragma omp parallel for
	for(i = 0; i < 3; ++i)
		for(j = 0; j < 4; ++j)
			mpfr_div(*(matrix + (i<<2) + j), *(matrix + (i<<2) + j), *(matrix + (i<<2) + i), MPFR_RNDN); 

    msprintf(COLOR_USER, "\nPARABOLIC CURVE Fitting is:\n");
    PRINTL();

    for(i = 0; i < 3; ++i)
    {
        mpfr_sprintf(buf, "%c => %Rf;\n", i+97, *(matrix + (i<<2) + 3));
        msprintf(COLOR_USER, buf);
    }

	mpfr_clears(tmp, ratio, NULL); 
    matrixFree(&matrix, ((dim_typ2){3, 4}));
    PRINTL();
    PRINTN();
    return;
}

__MSSHELL_WRAPPER_ static void  linearSystemsSolver(const sel_typ argc, char ** argv)
{
    dim_typ dim[2];
    mpfr_t *matrix;

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
        	strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[COLUMNS] != dim[ROWS]+1)
        {
            matrixFree(&matrix, dim);
            printUsage(&adv_calc[ADVCALC_LINEARSYSTEMSSOLVER]);
            return;
        }
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter Complete Matrix correspondent to the Linear System you want the Program to solve.\n\n");
        if((!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], false)) || dim[COLUMNS] != dim[ROWS]+1)
        {
            if(dim[COLUMNS] != dim[ROWS]+1)
                printErr(33, "You have to insert an [R X R+1] Matrix");
            return;
        }
    }

	char buf[MAX_BUFSIZ];
    dim_typ i, j, k;
    mpfr_t a, b;

    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIME);

    if(difftime)
    	gettimeofday(&tvBegin, NULL);
    	
    mpfr_inits(a, b, NULL); 

    for(i = 0; i < dim[ROWS]; ++i)
        for(j = 0; j < dim[ROWS]; ++j)
            if(i != j)
            {
            	mpfr_set(a, *(matrix + (dim[COLUMNS]*j) + i), MPFR_RNDN); 
            	mpfr_set(b, *(matrix + (dim[COLUMNS]*i) + i), MPFR_RNDN);
                for(k = 0; k < dim[COLUMNS]; ++k)
                {
                	mpfr_div(a, a, b, MPFR_RNDN);
					mpfr_mul(a, a, *(matrix + (dim[COLUMNS]*i) + k), MPFR_RNDN);
					mpfr_sub(*(matrix + (dim[COLUMNS]*j) + k), *(matrix + (dim[COLUMNS]*j) + k), a, MPFR_RNDN); 
                }
            }

	#pragma omp parallel for
    for(i = 0; i < dim[ROWS]; ++i)
        for(j = 0; j < dim[COLUMNS]; ++j)
        	mpfr_div(*(matrix + (dim[COLUMNS]*i) + j), *(matrix + (dim[COLUMNS]*i) + j), *(matrix + (dim[COLUMNS]*i) + i), MPFR_RNDN); 

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    msprintf(COLOR_USER, "Simultaneous Solutions of the given Linear System are:\n\n");
    PRINTL();

    for(i = 0; i < dim[ROWS] ; ++i)
    {
        mpfr_sprintf(buf, "%c => %Rf;\n", i+97, *(matrix + (dim[COLUMNS]*i) + dim[ROWS]));
        msprintf(COLOR_USER, buf);
    }

	mpfr_clears(a, b, NULL); 
    matrixFree(&matrix, dim);
    PRINTL();
    PRINTN();
    return;
}
