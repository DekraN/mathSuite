// adv_calc.c 10/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"


__MSSHELL_WRAPPER_ static void _MS__private secondGradeEquationSolver(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private complexAdd(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private complexMul(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private getDate(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private simplexMethod(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private newtonDifferenceTables(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private lagrangeInterpolation(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private greatestEigenValue(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private funcIntegration(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private straightLineFitting(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private parabolicCurveFitting(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private linearSystemsSolver(const sel_typ argc, char ** argv);


sprog adv_calc[MAX_ADVCALC_PROGS] =
{
    [ADVCALC_SECONDGRADEEQUATIONSOLVER] =
    {
        "Second Grade Equations Solver",
        CMD_SECONDGRADEQSOLVER,
        USAGE_SECONDGRADEQSOLVER,
        secondGradeEquationSolver,
        BY_USER,
        CHILD
    },
    [ADVCALC_COMPLEXNUMBERSSUM] =
    {
        "Complex and HyperComplex Numbers Addition",
        CMD_COMPLEXADD,
        USAGE_COMPLEXADD,
        complexAdd,
        BY_USER,
        CHILD
    },
    [ADVCALC_COMPLEXNUMBERSPROD] =
    {
        "Complex and HyperComplex Numbers Multiplication",
        CMD_COMPLEXMUL,
        USAGE_COMPLEXMUL,
        complexMul,
        BY_USER,
        CHILD
    },
    [ADVCALC_SIMPLEXMETHOD] =
    {
        "Non-Dual Simplex Method",
        CMD_SIMPLEXMETHOD,
        USAGE_SIMPLEXMETHOD,
        simplexMethod,
        BY_USER,
        CHILD
    },
    [ADVCALC_NEWTONDIFFTABLES] =
    {
        "Newton Difference Tables",
        CMD_NEWTONDIFFTABLES,
        USAGE_NEWTONDIFFTABLES,
        newtonDifferenceTables,
        BY_USER,
        CHILD
    },
    [ADVCALC_LAGRANGEINTERPOLATION] =
    {
        "Lagrange Unequal Interpolation",
        CMD_LAGRANGEINTERPOLATION,
        USAGE_LAGRANGEINTERPOLATION,
        lagrangeInterpolation,
        BY_USER,
        CHILD
    },
    [ADVCALC_GREATESTEIGENVALUE] =
    {
        "Greatest Eigen Value",
        CMD_GREATESTEIGENVALUE,
        USAGE_GREATESTEIGENVALUE,
        greatestEigenValue,
        BY_USER,
        CHILD
    },
    [ADVCALC_FUNCTIONINTEGRATION] =
    {
        "Function Integration",
        CMD_FID,
        USAGE_FID,
        funcIntegration,
        BY_USER,
        CHILD
    },
    [ADVCALC_STRAIGHTLINEFITTING] =
    {
        "Straight Line Fitting",
        CMD_STRAIGHTLINEFITTING,
        USAGE_STRAIGHTLINEFITTING,
        straightLineFitting,
        BY_USER,
        CHILD
    },
    [ADVCALC_PARABOLICCURVEFITTING] =
    {
        "Parabolic Curve Fitting",
        CMD_PARABOLICCURVEFITTING,
        USAGE_PARABOLICCURVEFITTING,
        parabolicCurveFitting,
        BY_USER,
        CHILD
    },
    [ADVCALC_LINEARSYSTEMSSOLVER] =
    {
        "Linear Systems Solver",
        CMD_LINEARSYSTEMSSOLVER,
        USAGE_LINEARSYSTEMSSOLVER,
        linearSystemsSolver,
        BY_USER,
        CHILD
    }
};

#define PARSING_ADVCALC_ALLOWED isSett(BOOLS_ADVCALCPARSING)

__MSSHELL_WRAPPER_ static void _MS__private secondGradeEquationSolver(const sel_typ argc, char ** argv)
{
    ityp *abc = NULL;

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &abc, dim, &dim[COLUMNS])) || dim[ROWS] != 1 || dim[COLUMNS] != MAX_ABSTRACT_DIMENSIONS)
        {
            matrixFree(&abc);
            printUsage(&adv_calc[ADVCALC_SECONDGRADEEQUATIONSOLVER]);
            return;
        }
    }
    else
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter COEFFICIENTS: 'a', 'b' e 'c' of SECOND GRADE Equation:\n\"a*(x^2) + b*x + c\"");
        printf2(COLOR_CREDITS, "\nby inserting related inline [1 x 3] Matrix.\n\n");

       if(PARSING_ADVCALC_ALLOWED)
            PRINTHOWTOBACKMESSAGE();

        if(!insertNMMatrix(&abc, (dim_typ2){1, MAX_ABSTRACT_DIMENSIONS}))
            return;
    }

    ityp root[MAX_DIMENSIONS];
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
    
    if(difftime)
    	gettimeofday(&tvBegin, NULL);

    if(_secondGradeEquationSolver(abc, root))
    {
        printf2(COLOR_USER, "\n1st ROOT = ");
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, root[ROOT_X1]);
        printf2(COLOR_USER, ";\n2nd ROOT = ");
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, root[ROOT_X2]);
        printf2(COLOR_USER, ".\n\n");
    }
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
        PRINTL();
	}
    
    matrixFree(&abc);

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private complexAdd(const sel_typ argc, char ** argv)
{
    ityp *cpx = NULL;
    const sel_typ algebra_units = (!access(curLayout)->algebra) ? MAX_COMPLEX_UNITS : exp2(access(curLayout)->algebra);

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &cpx, dim, &dim[COLUMNS])) || dim[ROWS] != MAX_DIMENSIONS || dim[COLUMNS] != algebra_units)
        {
            matrixFree(&cpx);
            printUsage(&adv_calc[ADVCALC_COMPLEXNUMBERSSUM]);
            return;
        }
    }
    else
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif

        printf2(COLOR_CREDITS, "\nEnter 2x%hu MATRIX whose ROWS contains respectively\nREAL PART and IMAGINARY PART%s of both two Operands.\n\n", algebra_units, algebra_units > MAX_COMPLEX_UNITS ? "s":NULL_CHAR);

        if(!insertNMMatrix(&cpx, (dim_typ2){MAX_DIMENSIONS, algebra_units}))
            return;
    }

    ityp complexRes[algebra_units];
    
    static void (* const complexAddFunc[_MAX_ALGEBRA][MAX_DIMENSIONS])(ityp *restrict, ityp [static MAX_SEDENIONS_UNITS]) =
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
	
	struct timeval tvBegin;
	ityp atime;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	complexAddFunc[((access(curLayout)->algebra == ALGEBRA_COMPLEXNUMBERS || !access(curLayout)->algebra) ? ALGEBRA_COMPLEXNUMBERS : access(curLayout)->algebra)-1][INVERSE_OPS](cpx, complexRes);
    
	if(difftime)
		atime = getDiffTime(&tvBegin);
    
	printf2(COLOR_USER, "\nRESULT of Operation: (");
	
    dim_typ i;

    PRINT2N();

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, *(cpx + i));
        printf2(COLOR_USER, "%s%s ", suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") ":" +");
    }
    
    static char oprchar[MAX_DIMENSIONS] = "+-";
    printf("%c\n", oprchar[INVERSE_OPS]);

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*SECOND_NUMBER) + i));
        printf2(COLOR_USER, "%s%s ", suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") is = to:\n":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, complexRes[i]);
        printf2(COLOR_USER, "%s%s", suite_c.algebra_imaginary_units_names[algebra_units][i], i ==algebra_units-1 ? ";\n\n":" + ");
    }

    matrixFree(&cpx);
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, atime);
        PRINTL();
	}

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINSO

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private complexMul(const sel_typ argc, char ** argv)
{
    ityp *cpx = NULL;
    const fsel_typ algebra_units = (!access(curLayout)->algebra) ? MAX_DIMENSIONS : exp2(access(curLayout)->algebra);

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &cpx, dim, &dim[COLUMNS])) || dim[ROWS] != MAX_DIMENSIONS || dim[COLUMNS] != MAX_DIMENSIONS)
        {
            matrixFree(&cpx);
            printUsage(&adv_calc[ADVCALC_COMPLEXNUMBERSPROD]);
            return;
        }
    }
    else
    {
        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter 2x%hu MATRIX whose ROWS contains respectively\nREAL PART and IMAGINARY PART%s of both two Operands.\n\n", algebra_units, algebra_units > MAX_COMPLEX_UNITS ? "s":NULL_CHAR);

        if(!insertNMMatrix(&cpx, (dim_typ2){MAX_DIMENSIONS, algebra_units}))
            return;
    }

    ityp complexRes[algebra_units];
    
    static void (* const complexMulFunc[_MAX_ALGEBRA][MAX_DIMENSIONS])(ityp *restrict, ityp [static MAX_SEDENIONS_UNITS]) =
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
	
	struct timeval tvBegin;
	ityp atime;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
	
	complexMulFunc[((access(curLayout)->algebra == ALGEBRA_COMPLEXNUMBERS || !access(curLayout)->algebra) ? ALGEBRA_COMPLEXNUMBERS : access(curLayout)->algebra)-1][INVERSE_OPS](cpx, complexRes);
    if(difftime)
    	atime = getDiffTime(&tvBegin);
    	
	printf2(COLOR_USER, "\nRESULT of Operation: (");

    dim_typ i;
    
    PRINT2N();

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*FIRST_NUMBER) + i));
        printf2(COLOR_USER, "%s%s ", suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") ":" +");
    }
    
    static char oprchar[MAX_DIMENSIONS] = "*/";
    printf("%c\n", oprchar[INVERSE_OPS]);

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, *(cpx + (algebra_units*SECOND_NUMBER) + i));
        printf2(COLOR_USER, "%s%s ", suite_c.algebra_imaginary_units_names[algebra_units][i], i==algebra_units-1 ? ") is = to:\n":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, complexRes[i]);
        printf2(COLOR_USER, "%s%s", suite_c.algebra_imaginary_units_names[algebra_units][i], i ==algebra_units-1 ? ";\n\n":" + ");
    }

    matrixFree(&cpx);
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, atime);
        PRINTL();
	}

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private simplexMethod(const sel_typ argc, char ** argv)
{

    sel_typ mode;
    ityp *tableau = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {
        if((mode = strtod(argv[0], NULL)) == MAX_DIMENSIONS) return;
        if(mode < MIN_PROBLEM || mode > MAX_DIMENSIONS)
        {
            printUsage(&adv_calc[ADVCALC_SIMPLEXMETHOD]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nSelect Simplex Method's Problem Type:\n");
        printf2(COLOR_CREDITS, "- A for min Problem,\n- B for max Problem;\n");
        printf2(COLOR_CREDITS, "- %c to go Back...\n\n", access(curLayout)->exit_char);

        do if((mode = toupper(getch())) == access(curLayout)->exit_char) return;
        while(mode < 'A' && mode > access(curLayout)->exit_char);

        mode -= 'A';
    }

    if(argc > 1)
    {
        if(!matrixToken(argv[1], &tableau, dim, &dim[COLUMNS]))
        {
            matrixFree(&tableau);
            printUsage(&adv_calc[ADVCALC_SIMPLEXMETHOD]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "Enter the Constraints' Coefficients Matrix, and use the last row\nto insert the Target Function Coefficients.\n");
        printf2(COLOR_CREDITS, "NOTE: In the last row you have to insert an element\nbefore exiting Insert Process, in order to align Matrix Dimensions.\n\n");
        if(!insertMatrix(tableau, dim[ROWS], dim[COLUMNS], false))
            return;
    }

    dim_typ i;
    dim_typ dimc[MAX_DIMENSIONS];
    const dim_typ dimrows_minus1 = dim[ROWS]-1;

    ityp *constraint_types = NULL;

    if(argc > 2)
    {
        if((!matrixToken(argv[2], &constraint_types, dimc, &dimc[COLUMNS])) || dimc[ROWS] != 1 || dimc[COLUMNS] != dim[ROWS])
        {
            matrixFree(&tableau);
            printUsage(&adv_calc[ADVCALC_SIMPLEXMETHOD]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nEnter Constraints Types: 0 for <=, non-zero element for >=.\n");
        if(!insertNMMatrix(&constraint_types, (dim_typ2){1,dimrows_minus1}))
        {
            matrixFree(&tableau);
            return;
        }
    }

    ityp *bfs = NULL;
    const dim_typ bfsdims[MAX_DIMENSIONS] =
    {
        1,
        dim[COLUMNS]-1
    };

    if(!matrixAlloc(&bfs, bfsdims))
    {
        matrixFree(&tableau);
        return;
    }

    sel_typ exit_state;

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
    if((exit_state = _simplexMethod(&tableau, &bfs, dim, constraint_types, mode)) == SIMPLEXMETHOD_INFBFS_ERROR)
        printErr(33, "This Problem has a Solution whose limit is Infinite");
    else if(exit_state == SIMPLEXMETHOD_FARBFS_ERROR)
        printErr(33, "No convergence after %hu iterations!", access(curLayout)->max_simplex_iterations);
    else if(exit_state == SIMPLEXMETHOD_ALLOC_ERROR)
        printErr(12, "Simplex Method Heap Dynamic Memory Allocation Problem");
    else
    {
        printf2(COLOR_USER, "\nRelaxed Problem BFS with Artificial Variables is: ");
        printMatrix(stdout, bfs, (dim_typ2){1,dim[ROWS]+dim[COLUMNS]-2});
    }
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
        PRINTL();
	}

    matrixFree(&tableau);
    matrixFree(&constraint_types);
    matrixFree(&bfs);

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}


__MSSHELL_WRAPPER_ static void _MS__private newtonDifferenceTables(const sel_typ argc, char ** argv)
{

    ityp x[access(curLayout)->max_newton_difftables_dim];
    ityp y[access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim];

    dim_typ n;
    ityp tmp;

    if(argc)
    {
        if(PARSING_ADVCALC_ALLOWED)
        {
            if((!parse(argv[0], &tmp)) || tmp != (n = (dim_typ)tmp) || n < access(curLayout)->min_newton_difftables_dim || n > access(curLayout)->max_newton_difftables_dim)
            {
                printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                return;
            }
        }
        else if((tmp = strtod(argv[0], NULL)) != (n = (dim_typ)tmp) || n < access(curLayout)->min_newton_difftables_dim || n > access(curLayout)->max_newton_difftables_dim)
        {
            printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nEnter Difference Table DIMENSION.\n");

        if(PARSING_ADVCALC_ALLOWED)
            PRINTHOWTOBACKMESSAGE();
        while((PARSING_ADVCALC_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted Difference Table DIMENSION is:", PARSER_NOSETTINGS)))) :
                (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (n = (dim_typ)tmp) || n < access(curLayout)->min_newton_difftables_dim || n > access(curLayout)->max_newton_difftables_dim)
        {
            CLEARBUFFER();
            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
            if(exitHandleCheck) return;
            printErr(5, "Invalid inserted Value.\nMust be a non-negative integer between %hu and %hu", access(curLayout)->min_newton_difftables_dim, access(curLayout)->max_newton_difftables_dim);
        }
    }

    dim_typ i;

    if(argc > 1)
        if(argc == n+1)
        {
            char *token = NULL;
			for(i=0; i<n; ++i)
            {

                if((token = strtok(argv[i+1], ",")))
                {
                    if(PARSING_ADVCALC_ALLOWED)
                    {
                        if(!parse(token, &x[i]))
                        {
                            printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                            return;
                        }
                    }
                    else
                        x[i] = strtod(token, NULL);
                }
                else
                {
                    printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                    return;
                }

                if((token = strtok(NULL, ",")))
                {
                    if(PARSING_ADVCALC_ALLOWED)
                    {
                        if(!parse(token, &y[i][0]))
                        {
                            printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                            return;
                        }
                    }
                    else
                        y[i][0] = strtod(token, NULL);
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
            printf2(COLOR_CREDITS, "Enter couple No. %hu as expected format:\n[X%sY]\n", i, PARSING_ADVCALC_ALLOWED ? "]\n[" : " ");
            /// scanf("%f %f",&x[i],&y[i][0]); // PAY STRICT ATTENTION TO THIS HANDLE

            while((PARSING_ADVCALC_ALLOWED ? (isNullVal((x[i] = requires(NULL, NULL_CHAR, "Inserted X is:", PARSER_SHOWRESULT)))) ||
                (isNullVal((y[i][0] = requires(NULL, NULL_CHAR, "Inserted Y is:", PARSER_SHOWRESULT)))) :
                !scanf2(2, "%lf %lf", &x[i], &y[i][0])))
            {
                CLEARBUFFER();
                if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
                if(exitHandleCheck) return;
                printErr(5, "Invalid inserted Value");
            }
        }

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
    newtonDifferenceTable(n, y, FORWARD_DIFFTAB);
    showNewtonDifferenceTable(n, x, y, FORWARD_DIFFTAB);
    newtonDifferenceTable(n, y, BACKWARD_DIFFTAB);
    showNewtonDifferenceTable(n, x, y, BACKWARD_DIFFTAB);
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
        PRINTL();
	}

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private lagrangeInterpolation(const sel_typ argc, char ** argv)
{
    ityp tmp;
    dim_typ dim;

    if(argc)
    {
        if(PARSING_ADVCALC_ALLOWED)
        {
            if((!parse(argv[0], &tmp)) || tmp != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
            {
                printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
                return;
            }
        }
        else if((tmp = strtod(argv[0], NULL)) != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
        {
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nEnter Data DIMENSION.\n\n");

        while((PARSING_ADVCALC_ALLOWED ? isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT))) :
            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
        {
            CLEARBUFFER();
            if(access(exitHandle) == EXITHANDLE_GETCMD) continue; // Highly experimental
            if(exitHandleCheck) return;
            printErr(33, "Invalid inserted Data DIMENSION.\nMust be an integer between 1 and %z", USHRT_MAX);
        }
    }

    ityp *xy = NULL;
    dim_typ i, j;

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy);
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
    }
    else
    {
        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter related Matrix filled with Data you want Interpolation to Process,\n");
        printf2(COLOR_CREDITS, "by putting on each rows the %hu X and Y Values.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){MAX_DIMENSIONS, dim}))
            return;
    }

    ityp xp;

    if(argc > 2)
        if(PARSING_ADVCALC_ALLOWED)
        {
            if(!parse(argv[2], &xp))
            {
                printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
                #ifdef WINOS
                    SetExitButtonState(ENABLED);
                #endif // WINOS
                return;
            }
        }
        else
            xp = strtod(argv[2], NULL);
    else
    {
        printf2(COLOR_CREDITS, "Enter X VALUE to find Y one.\n");
        while((PARSING_ADVCALC_ALLOWED ? (isNullVal((xp = requires(NULL, NULL_CHAR, "X VALUE | Y = F(X) is:", PARSER_SHOWRESULT)))) :
            !scanf2(1, INPUT_CONVERSION_FORMAT, &xp)))
        {
            CLEARBUFFER();
            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
            if(exitHandleCheck)
            {
                #ifdef WINOS
                    SetExitButtonState(ENABLED);
                #endif // WINOS
                return;
            }
            printErr(5, "Invalid inserted Value");
        }
    }

    ityp yp = 0;
    ityp dr, nr;

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	for(i = 0; i < dim; ++i)
    {

        dr = nr = 1.00;
        for(j = 0; j<dim; ++j)
            if(i!=j)
            {
                nr *= xp - *(xy + (dim*XROW) + j);
                dr *= *(xy + (dim*XROW) + i) - *(xy + (dim*XROW) + j);
            }

        yp += nr/dr*(*(xy + (dim*YROW) + i));
    }
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
        PRINTL();
	}

    matrixFree(&xy);

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    printf2(COLOR_USER, "\nRequested Y VALUE is: ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, yp);
    printf2(COLOR_USER, ".\n\n");

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private greatestEigenValue(const sel_typ argc, char ** argv)
{

    ityp **matrix1 = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    matrix1 = malloc(sizeof(ityp*));
    errMem(matrix1, VSPACE);

    if(argc)
    {
        if((!matrixToken(argv[0], matrix1, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
            matrixFree(matrix1);
            free(matrix1);
            printUsage(&adv_calc[ADVCALC_GREATESTEIGENVALUE]);
            return;
        }
    }
    else if(!enterMatrix(matrix1, dim, &dim[COLUMNS], true, true))
    {
        free(matrix1);
        return;
    }

    ityp **matrix2 = NULL;
    matrix2 = malloc(sizeof(ityp*));
    errMem(matrix2, (matrixFree(matrix1), free(matrix1)));

    if(!matrixAlloc(matrix2, (dim_typ2){dim[ROWS], 1}))
    {
        matrixFree(matrix1);
        free(matrix1);
        free(matrix2);
        #ifdef WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    ityp eigenValue;

    dim_typ i;

    #pragma omp parallel for
	for(i=0; i<dim[ROWS]; ++i)
        *((*matrix2) + i) = 1.00;

    ityp **result = NULL;
    
    result = malloc(sizeof(ityp*));
    errMem(result, (matrixFree(matrix1), matrixFree(matrix2), free(matrix1), free(matrix2)));

    if(!matrixAlloc(result, (dim_typ2){dim[ROWS], 1}))
    {
        matrixFree(matrix1);
        matrixFree(matrix2);
        free(matrix1);
        free(matrix2);
        free(result);
        #ifdef WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    ityp resultVector[dim[ROWS]];
	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);

    for( ;; )
    {
        _matrixMultiplication(matrix1, matrix2, result, (dim_typ3){dim[ROWS], dim[ROWS], 1});
        // matrixToVector((*result), (dim_typ2){dim[ROWS], 1}, resultVector, MATRIX_TO_VECTOR);
        eigenValue = MAX(3, resultVector);

        #pragma omp parallel for
		for(i=0; i<dim[ROWS]; ++i)
            *((*result) + i) /= eigenValue;

        if(isEqualMatrix((*matrix2), (*result), (dim_typ2){dim[ROWS], 1}))
            break;
	
		#pragma omp parallel for
        for(i=0; i<dim[ROWS]; ++i)
            *((*matrix2) + i) = *((*result) + i);
    }
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    printf2(COLOR_USER, "\nGreatest EIGEN Value is: ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, eigenValue);
    printf2(COLOR_USER, ".\nAn EIGEN Vector is:\n\n");

    printMatrix(stdout, (*result), (dim_typ2){dim[ROWS],1});

    matrixFree(matrix1);
    matrixFree(matrix2);
    matrixFree(result);

    free(matrix1);
    free(matrix2);
    free(result);

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private funcIntegration(const sel_typ argc, char ** argv)
{
    ityp x0, xn;
    ityp h, s;

    dim_typ funcID;
    dim_typ j;

    funcID = selectListItem(MAX_FIDS, MAX_FIDS > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
                            "Select desired Function you want to Integrate", ext_math.funcnames);

    if(funcID == MAX_FIDS) return;

    ityp tmp;


    sel_typ mode;



    if(argc)
    {
        if((mode = strtod(argv[0], NULL)) == MAX_ABSTRACT_DIMENSIONS) return;
        if(mode < 0 || mode > MAX_ABSTRACT_DIMENSIONS)
        {
            printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nSelect Defined Integration Calculus Mode:\n");
        printf2(COLOR_CREDITS, "- A for Simpsons' 1/3 method,\n- B for Simpsons' 3/8 method,\n- C per Trapezoidal Method;\n");
        printf2(COLOR_CREDITS, "- D to go Back...\n\n");

        do if((mode = toupper(getch())) == 'D') return;
        while(mode < 'A' && mode > 'D');

        mode -= 'A';
    }

    dim_typ i, n;

    if(argc > 1)
        if(argc == 4)
        {
            if(PARSING_ADVCALC_ALLOWED)
            {
                if(!parse(argv[1], &x0))
                {
                    printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                    return;
                }

                if(!parse(argv[2], &xn))
                {
                    printUsage(&adv_calc[ADVCALC_NEWTONDIFFTABLES]);
                    return;
                }

                if((!parse(argv[3], &tmp)) || tmp != (n = (dim_typ)tmp))
                {
                    printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
                    return;
                }
            }
            else
            {
                x0 = strtod(argv[1], NULL);
                xn = strtod(argv[2], NULL);
                if((tmp = strtod(argv[3], NULL)) != (n = (dim_typ)tmp))
                {
                    printUsage(&adv_calc[ADVCALC_FUNCTIONINTEGRATION]);
                    return;
                }
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

        strcpy(seperator, PARSING_ADVCALC_ALLOWED ? "]\n[" : BLANK_STRING);
        printf2(COLOR_CREDITS, "\nEnter INTEGRATION Extremes and Intervals NUMBER as expected format:\n");
        printf2(COLOR_CREDITS, "[x0%sxN%sNo]\n", seperator, seperator);

        while((PARSING_ADVCALC_ALLOWED ? (isNullVal((x0 = requires(NULL, NULL_CHAR, "First INTEGRATION Extreme is:", PARSER_SHOWRESULT)))) ||
            (isNullVal((xn = requires(NULL, NULL_CHAR, "Second INTEGRATION Extreme is:", PARSER_SHOWRESULT)))) || (isNullVal((tmp = requires(NULL, NULL_CHAR, "Il Numero di INTERVALLI inserito e'", PARSER_SHOWRESULT)))) :
            (!scanf2(3, "%lf %lf %lf", &x0, &xn, &tmp))) || tmp != (n = (dim_typ)tmp) || n < 1 || n > INT_MAX)
        {
            CLEARBUFFER();
            if(access(exitHandle) == EXITHANDLE_GETCMD) continue; // Highly experimental
            if(exitHandleCheck) return;
            printErr(33, "Invalid inserted Intervals NUMBER.\nMust be an integer between 1 and %z", INT_MAX);
        }
    }

    ityp result = 0.00;

    h = (xn - x0) / n;

    ityp (* const y)(register ityp) = ext_math.functions[funcID];
    
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
    
    if(difftime)
    	gettimeofday(&tvBegin, NULL);

    switch(mode)
    {
        case SIMPSON1DIV8_RULE:
            s = y(x0)+y(xn)+4*y(x0+h);
            for(i = 3; i<=n-1; i+=2)
                s += 4*y(x0+i*h) + 2*y(x0+(i-1)*h);

            result = (h/3)*s;
            break;

        case SIMPSON3DIV8_RULE:
        {
            bool flag;

            s = y(x0)+y(xn);

            for(i = 1; i<=n-1;++i)
            {
                for(j=1;j<=n-1;++j)
                    if((flag = i == 3*j))
                        break;
                s += flag ? 2*y(x0+i*h) : 3*y(x0+i*h);
            }

            result = (3*h/8)*s;
            break;
        }

        case TRAPEZOIDAL_RULE:
            s = y(x0) + y(xn);
            for(i = 0; ++i < n; )
                s += 2*y(x0+i*h);
            result = (h*0.5)*s;
            break;
    }
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    printf2(COLOR_USER, "%s(x) Function Integral Value calculated between: x0 = ", ext_math.funcnames[funcID]);
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, x0);
    printf2(COLOR_USER, " and xN = ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, xn);
    printf2(COLOR_USER, ",\nwith %hu Intervals NUMBER is: %6.*f\n\n", n, access(curLayout)->precision, result);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private straightLineFitting(const sel_typ argc, char ** argv)
{
    ityp tmp;
    dim_typ dim;

    if(argc)
    {
        if(PARSING_ADVCALC_ALLOWED)
        {
            if((!parse(argv[0], &tmp)) || tmp != (dim = (dim_typ)tmp))
            {
                printUsage(&adv_calc[ADVCALC_STRAIGHTLINEFITTING]);
                return;
            }
        }
        else if((tmp = strtod(argv[0], NULL)) != (dim = (dim_typ)tmp))
        {
            printUsage(&adv_calc[ADVCALC_STRAIGHTLINEFITTING]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nEnter Data DIMENSION.\n\n");
        while((PARSING_ADVCALC_ALLOWED ? isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT))) :
            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
        {
            CLEARBUFFER();
            if(access(exitHandle) == EXITHANDLE_GETCMD) continue; // Highly experimental
            if(exitHandleCheck) return;
            printErr(33, "Invalid inserted Data DIMENSION.\nMust be an integer between 1 and %zu", USHRT_MAX);
        }
    }

    ityp *xy = NULL;

    dim_typ i;

    // we must seek for the BACKTRACKING FEATURE

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy);
            printUsage(&adv_calc[ADVCALC_STRAIGHTLINEFITTING]);
            return;
        }
    }
    else
    {

        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter the Matrix filled with EQUATIONS DATA,\n");
        printf2(COLOR_CREDITS, "by putting on each ROWS the %hu X and Y VALUES.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){MAX_DIMENSIONS, dim}))
            return;
    }

    ityp sum_x, sum_xy, sum_x2, sum_y;

    sum_x = sum_xy = sum_x2 = sum_y = 0.00;
    
    const register dim_typ cache[MAX_DIMENSIONS] =
    {
    	dim*XROW,
    	dim*YROW
    };
    
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
    
    if(difftime)
    	gettimeofday(&tvBegin, NULL);
    	
	for(i = 0; i < dim; ++i)
    {

        sum_x += *(xy + cache[XROW] + i);
        sum_y += *(xy + cache[YROW] + i);
        sum_xy += *(xy + cache[XROW] + i) * *(xy + cache[YROW] + i);
        sum_x2 += pow(*(xy + cache[XROW] + i), 2);
    }
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    matrixFree(&xy);

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    const register ityp b = (dim*sum_xy - sum_x*sum_y)/(dim*sum_x2 - pow(sum_x,2));
    const register ityp a = (sum_y - b*sum_x)/dim;

    printf2(COLOR_USER, "\na = ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, a);
    printf2(COLOR_USER, "; b = ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, b);
    printf2(COLOR_USER, ".\nThe Equation is ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, a);
    printf2(COLOR_USER, " + ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, b);
    printf2(COLOR_USER, "X.\n\n");

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private parabolicCurveFitting(const sel_typ argc, char ** argv)
{
    ityp tmp;
    dim_typ dim;

    if(argc)
    {
        if(PARSING_ADVCALC_ALLOWED)
        {
            if((!parse(argv[0], &tmp)) || tmp != (dim = (dim_typ)tmp))
            {
                printUsage(&adv_calc[ADVCALC_PARABOLICCURVEFITTING]);
                return;
            }
        }
        else if((tmp = strtod(argv[0], NULL)) != (dim = (dim_typ)tmp))
        {
            printUsage(&adv_calc[ADVCALC_PARABOLICCURVEFITTING]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nEnter Data DIMENSION.\n\n");
        while((PARSING_ADVCALC_ALLOWED ? isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted Data DIMENSION is:", PARSER_SHOWRESULT))) :
            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (dim = (dim_typ)tmp) || dim < 1 || dim > USHRT_MAX)
        {
            CLEARBUFFER();
            if(access(exitHandle) == EXITHANDLE_GETCMD) continue; // Highly experimental
            if(exitHandleCheck) return;
            printErr(33, "Invalid inserted Data DIMENSIONM.\nMust be an integer between 1 and %z", USHRT_MAX);
        }
    }

    ityp *xy = NULL;


    // we must seek for the BACKTRACKING FEATURE

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[ROWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy);
            printUsage(&adv_calc[ADVCALC_PARABOLICCURVEFITTING]);
            return;
        }
    }
    else
    {
        #ifdef WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter the Matrix filled with EQUATIONS DATA,\n");
        printf2(COLOR_CREDITS, "by putting on each Rows the %hu X and Y VALUES.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){MAX_DIMENSIONS, dim}))
            return;
    }

    ityp *matrix = NULL;

    if(!matrixAlloc(&matrix, (dim_typ2){MAX_ABSTRACT_DIMENSIONS, 4}))
    {
        matrixFree(&xy);
        #ifdef WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }
    
    *(matrix) = dim;

    dim_typ i;
    
    const register dim_typ cache[MAX_DIMENSIONS] =
    {
    	dim*XROW,
    	dim*YROW
    };
    
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
    
    if(difftime)
    	gettimeofday(&tvBegin, NULL);

    #pragma omp parallel for
	for(i = 0; i < dim; ++i)
    {
        *(matrix + 1) = (*(matrix + 4) += *(xy + cache[XROW] + i));
        *(matrix + 3) += *(xy + cache[YROW] + i);
        *(matrix + 2) = (*(matrix + 8) += pow(*(xy + cache[XROW] + i), 2)); 
        *(matrix + 6) = (*(matrix + 9) += pow(*(xy + cache[XROW] + i), 3));
        *(matrix + 10) += pow(*(xy + cache[XROW] + i), 4); 
        *(matrix + 7) += (*(xy + cache[XROW] + i) * *(xy + cache[YROW] * i));
        *(matrix + 11) += (pow(*(xy + cache[XROW] + i), 2) * *(xy + cache[YROW]*i));
    }
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    matrixFree(&xy);

    dim_typ j;
    dim_typ k;

    #pragma omp parallel for
	for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
        #pragma omp parallel for
		for(j = 0; j < MAX_ABSTRACT_DIMENSIONS; ++j)
            if(i!=j)
            {
                const ityp ratio = *(matrix + (j<<MAX_DIMENSIONS) + i)/ *(matrix + (i<<MAX_DIMENSIONS) + i);
                for(k = 0; k < 4; ++k)
                    *(matrix + (j<<MAX_DIMENSIONS) + k) -= ratio * *(matrix + (i<<MAX_DIMENSIONS) + k);
            }

    #pragma omp parallel for
	for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
    {
        const ityp a = *(matrix + (i<<MAX_DIMENSIONS) + i);
		for(j = 0; j < 4; ++j)
            *(matrix + (i<<MAX_DIMENSIONS) + j) /= a;
    }

    printf2(COLOR_USER, "\nPARABOLIC CURVE Fitting is:\n");
    PRINTL();

    for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
    {
        printf2(COLOR_USER, "%c => ", i+97);
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, *(matrix + (i<<MAX_DIMENSIONS) + 3));
        printf2(COLOR_USER, ";\n");
    }

    matrixFree(&matrix);

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    PRINTL();
    PRINTN();

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private linearSystemsSolver(const sel_typ argc, char ** argv)
{
    dim_typ dim[MAX_DIMENSIONS];
    ityp *matrix = NULL;

    if(argc)
    {

        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[COLUMNS] != dim[ROWS]+1)
        {
            matrixFree(&matrix);
            printUsage(&adv_calc[ADVCALC_LINEARSYSTEMSSOLVER]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nEnter Complete Matrix correspondent to the Linear System you want the Program to solve.\n\n");
        if((!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], false)) || dim[COLUMNS] != dim[ROWS]+1)
        {
            if(dim[COLUMNS] != dim[ROWS]+1)
                printErr(33, "You have to insert an [R X R+1] Matrix");
            return;
        }
    }

    dim_typ i, j, k;
    ityp a, b;
    
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIMEADVCALC);
    
    if(difftime)
    	gettimeofday(&tvBegin, NULL);

    for(i = 0; i < dim[ROWS]; ++i)
        for(j = 0; j < dim[ROWS]; ++j)
            if(i != j)
            {
                a = *(matrix + (dim[COLUMNS]*j) + i);
                b = *(matrix + (dim[COLUMNS]*i) + i);
                for(k = 0; k < dim[COLUMNS]; ++k)
                    *(matrix + (dim[COLUMNS]*j) + k) -= (a/b) * *(matrix + (dim[COLUMNS]*i) + k);
            }

	#pragma omp parallel for
    for(i = 0; i < dim[ROWS]; ++i)
    {
        const ityp c = *(matrix + (dim[COLUMNS]*i) + i);
        for(j = 0; j < dim[COLUMNS]; ++j)
            *(matrix + (dim[COLUMNS]*i) + j) /= c;
    }
    
    if(difftime)
    {
    	PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    printf2(COLOR_USER, "Simultaneous Solutions of the given Linear System are:\n\n");
    PRINTL();

    for(i = 0; i < dim[ROWS] ; ++i)
    {
        printf2(COLOR_USER, "%c => ", i+97);
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, *(matrix + (dim[COLUMNS]*i) + dim[ROWS]));
        printf2(COLOR_USER, ";\n");
    }

    matrixFree(&matrix);

    #ifdef WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    PRINTL();
    PRINTN();
    return;
}
