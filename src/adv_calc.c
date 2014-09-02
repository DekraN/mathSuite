// adv_calc.c 29/08/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"


__MSSHELL_WRAPPER_ static void _MS__private secondGradeEquationSolver(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private complexSum(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private complexProd(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private getDate(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private simplexMethod(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private newtonDifferenceTables(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private lagrangeInterpolation(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private greatestEigenValue(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private funcIntegration(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private straightLineFitting(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private parabolicCurveFitting(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private linearSystemsSolver(const register sel_typ argc, char ** argv);


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
        "Complex and HyperComplex Numbers Sum",
        CMD_COMPLEXSUM,
        USAGE_COMPLEXSUM,
        complexSum,
        BY_USER,
        CHILD
    },
    [ADVCALC_COMPLEXNUMBERSPROD] =
    {
        "Complex and HyperComplex Numbers Product",
        CMD_COMPLEXPROD,
        USAGE_COMPLEXPROD,
        complexProd,
        BY_USER,
        CHILD
    },
    [ADVCALC_GETFORMATTEDDATE] =
    {
        "Formatted Date by Numeric One",
        CMD_GETDATE,
        USAGE_GETDATE,
        getDate,
        BY_USER,
        CHILD,
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

__MSSHELL_WRAPPER_ static void _MS__private secondGradeEquationSolver(const register sel_typ argc, char ** argv)
{
    ityp **abc = NULL;

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &abc, dim, &dim[COLUMNS])) || dim[RAWS] != 1 || dim[COLUMNS] != MAX_ABSTRACT_DIMENSIONS)
        {
            matrixFree(&abc, dim[RAWS]);
            printUsage(&adv_calc[ADVCALC_SECONDGRADEEQUATIONSOLVER]);
            return;
        }
    }
    else
    {

        #if WINOS
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

    if(_secondGradeEquationSolver(abc[0], root))
    {
        printf2(COLOR_USER, "\n1st ROOT = ");
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, root[ROOT_X1]);
        printf2(COLOR_USER, ";\n2nd ROOT = ");
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, root[ROOT_X2]);
        printf2(COLOR_USER, ".\n\n");
    }

    matrixFree(&abc, 1);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private complexSum(const register sel_typ argc, char ** argv)
{
    ityp **complex = NULL;
    const fsel_typ algebra_units = (!access(curLayout)->algebra) ? MAX_DIMENSIONS : exp2(access(curLayout)->algebra);

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &complex, dim, &dim[COLUMNS])) || dim[RAWS] != MAX_DIMENSIONS || dim[COLUMNS] != algebra_units)
        {
            matrixFree(&complex, dim[RAWS]);
            printUsage(&adv_calc[ADVCALC_COMPLEXNUMBERSSUM]);
            return;
        }
    }
    else
    {

        #if WINOS
            SetExitButtonState(DISABLED);
        #endif

        printf2(COLOR_CREDITS, "\nEnter 2x%hu MATRIX whose RAWS contains respectively\nREAL PART and IMAGINARY PART%s of both two Operands.\n\n", algebra_units, algebra_units > MAX_COMPLEX_UNITS ? "s":NULL_CHAR);

        if(!insertNMMatrix(&complex, (dim_typ2){MAX_DIMENSIONS, algebra_units}))
            return;
    }

    ityp complexRes[algebra_units];

    switch(access(curLayout)->algebra)
    {
        case ALGEBRA_REALNUMBERS:
        case ALGEBRA_COMPLEXNUMBERS:
            _complexSum(complex, complexRes);
            break;
        case ALGEBRA_QUATERNIONS:
            _quaternionsSum(complex, complexRes);
            break;
        case ALGEBRA_OCTONIONS:
            _octonionsSum(complex, complexRes);
            break;
        case ALGEBRA_SEDENIONS:
            _sedenionsSum(complex, complexRes);
            break;
    }

    printf2(COLOR_USER, "\nRESULT of Operation: (");

    dim_typ i;

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, complex[FIRST_NUMBER][i]);
        printf2(COLOR_USER, "%s%s ", suite_c.algebra_imaginary_units_names[i], i==algebra_units-1 ? ")":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, complex[SECOND_NUMBER][i]);
        printf2(COLOR_USER, "%s%s ", suite_c.algebra_imaginary_units_names[i], i==algebra_units-1 ? ") is = to: ":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, complexRes[i]);
        printf2(COLOR_USER, i==algebra_units-1 ? ";\n\n":" + ");
    }

    matrixFree(&complex, MAX_DIMENSIONS);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINSO

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private complexProd(const register sel_typ argc, char ** argv)
{
    ityp **complex = NULL;
    const fsel_typ algebra_units = (!access(curLayout)->algebra) ? MAX_DIMENSIONS : exp2(access(curLayout)->algebra);

    if(argc)
    {
        dim_typ dim[MAX_DIMENSIONS];

        if((!matrixToken(argv[0], &complex, dim, &dim[COLUMNS])) || dim[RAWS] != MAX_DIMENSIONS || dim[COLUMNS] != MAX_DIMENSIONS)
        {
            matrixFree(&complex, dim[RAWS]);
            printUsage(&adv_calc[ADVCALC_COMPLEXNUMBERSPROD]);
            return;
        }
    }
    else
    {
        #if WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter 2x%hu MATRIX whose RAWS contains respectively\nREAL PART and IMAGINARY PART%s of both two Operands.\n\n", algebra_units, algebra_units > MAX_COMPLEX_UNITS ? "s":NULL_CHAR);

        if(!insertNMMatrix(&complex, (dim_typ2){MAX_DIMENSIONS, algebra_units}))
            return;
    }

    ityp complexRes[algebra_units];

    switch(access(curLayout)->algebra)
    {
        case ALGEBRA_REALNUMBERS:
        case ALGEBRA_COMPLEXNUMBERS:
            _complexProd(complex, complexRes);
            break;
        case ALGEBRA_QUATERNIONS:
            _quaternionsProduct(complex, complexRes);
            break;
        case ALGEBRA_OCTONIONS:
            _octonionsProduct(complex, complexRes);
            break;
        case ALGEBRA_SEDENIONS:
            _sedenionsProduct(complex, complexRes);
            break;
    }


    printf2(COLOR_USER, "\nRESULT of Operation: (");

    dim_typ i;

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, complex[FIRST_NUMBER][i]);
        printf2(COLOR_USER, "%s%s ", suite_c.algebra_imaginary_units_names[i], i==algebra_units-1 ? ")":" +");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, complex[SECOND_NUMBER][i]);
        printf2(COLOR_USER, "%s%s ", suite_c.algebra_imaginary_units_names[i], i==algebra_units-1 ? ") is = to: ":" *");
    }

    for(i=0; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, complexRes[i]);
        printf2(COLOR_USER, i==algebra_units-1 ? ";\n\n":" * ");
    }

    matrixFree(&complex, MAX_DIMENSIONS);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private getDate(const register sel_typ argc, char ** argv)
{
    ityp tmp;
    ityp tmp2, tmp3;

    ityp dd, mm;
    ityp yyyy;

    if(argc)
    {
        if(argc != MAX_ABSTRACT_DIMENSIONS)
        {
            printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
            return;
        }

        if(PARSING_ADVCALC_ALLOWED)
        {
            if((!parse(argv[0], &tmp)) || tmp != (dd = (sel_typ)tmp))
            {
                printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
                return;
            }
            if((!parse(argv[1], &tmp)) || tmp != (mm = (sel_typ)tmp))
            {
                printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
                return;
            }
            if((!parse(argv[2], &tmp)) || tmp != (yyyy = (uint64_t)tmp))
            {
                printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
                return;
            }
        }
        else
        {
            if((tmp = strtod(argv[0], NULL)) != (dd = (sel_typ)tmp))
            {
                printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
                return;
            }
            if((tmp = strtod(argv[1], NULL)) != (mm = (sel_typ)tmp))
            {
                printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
                return;
            }
            if((tmp = strtod(argv[2], NULL)) != (yyyy = (uint64_t)tmp))
            {
                printUsage(&adv_calc[ADVCALC_GETFORMATTEDDATE]);
                return;
            }
        }

    }
    else
    {
        char seperator[SIGN_STRING] = NULL_CHAR;

        strcpy(seperator, PARSING_ADVCALC_ALLOWED ? "]\n[" : BLANK_STRING);

        printf2(COLOR_CREDITS, "\nEnter Numeric DATE as expected format:\n");
        printf2(COLOR_CREDITS, "[DD%sMM%sYYYY]\n", seperator, seperator);

        if(PARSING_ADVCALC_ALLOWED)
            PRINTHOWTOBACKMESSAGE();

        while(PARSING_ADVCALC_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted DAY is", PARSER_SHOWRESULT))) || tmp != (dd = (sel_typ)tmp) || dd < 1 || dd > MAX_MONTH_DAYS) ||
        (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted MONTH is", PARSER_SHOWRESULT))) || tmp != (mm = (sel_typ)tmp) || mm < 1 || mm > MAX_MONTHS) ||
        (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted YEAR is", PARSER_SHOWRESULT))) || tmp != (yyyy = (uint64_t)tmp) || yyyy < 1)      :
        (!scanf2(3, "%lf %lf %lf", &tmp, &tmp2, &tmp3)) || tmp != (dd = (sel_typ)tmp) || tmp2 != (mm = (sel_typ)tmp2) || tmp3 != (yyyy = (uint64_t)tmp3)
          || dd < 1 || mm < 1 || yyyy < 1 || dd > MAX_MONTH_DAYS || mm > MAX_MONTHS)
        {
            CLEARBUFFER();
            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
            if(exitHandleCheck) return;
            printErr(33, "Invalid [DD MM YYYY] format.\nYou have to insert DAY, MONTH and YEAR as non-negative and non-zero integers.\n\
and the first two Params must respectively be less than %hu and %hu", MAX_MONTH_DAYS, MAX_MONTHS);
        }
    }

    printf2(COLOR_USER, "\nFormatted DATE correspondent to the Inserted Numeric ONE:\n%hu/%hu/%llu is: %s %hu %s %llu.\n\n",
           dd, mm, yyyy, getDayName(getDayNumber(dd, mm, yyyy)), dd, getMonthName(mm), yyyy);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private simplexMethod(const register sel_typ argc, char ** argv)
{

    sel_typ mode;
    ityp **tableau = NULL;
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
            matrixFree(&tableau, dim[RAWS]);
            printUsage(&adv_calc[ADVCALC_SIMPLEXMETHOD]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "Enter the Constraints' Coefficients Matrix, and use the last raw\nto insert the Target Function Coefficients.\n");
        printf2(COLOR_CREDITS, "NOTE: In the last raw you have to insert an element\nbefore exiting Insert Process, in order to align Matrix Dimensions.\n\n");
        if(!insertMatrix(tableau, dim[RAWS], dim[COLUMNS], false))
            return;
    }

    dim_typ i;
    dim_typ dimc[MAX_DIMENSIONS];
    const dim_typ dimraws_minus1 = dim[RAWS]-1;

    bool constraint_types[dimraws_minus1];
    ityp **constraint_types_matrix = NULL;

    if(argc > 2)
    {
        if((!matrixToken(argv[2], &constraint_types_matrix, dimc, &dimc[COLUMNS])) || dimc[RAWS] != 1 || dimc[COLUMNS] != dim[RAWS])
        {
            matrixFree(&tableau, dim[RAWS]);
            printUsage(&adv_calc[ADVCALC_SIMPLEXMETHOD]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nEnter Constraints Types: 0 for <=, non-zero element for >=.\n");
        if(!insertNMMatrix(&constraint_types_matrix, (dim_typ2){1,dimraws_minus1}))
        {
            matrixFree(&tableau, dim[RAWS]);
            return;
        }
    }

    #pragma omp parallel for
	for(i=0; i<dimraws_minus1; ++i)
        constraint_types[i] = constraint_types_matrix[RAWS][i] != 0.00;

    matrixFree(&constraint_types_matrix, 1);

    ityp **bfs = NULL;
    const dim_typ bfsdims[MAX_DIMENSIONS] =
    {
        1,
        dim[COLUMNS]-1
    };

    if(!matrixAlloc(&bfs, bfsdims))
    {
        matrixFree(&tableau, dim[RAWS]);
        return;
    }

    sel_typ exit_state;

    if((exit_state = _simplexMethod(&tableau, &bfs, dim, constraint_types, mode)) == SIMPLEXMETHOD_INFBFS_ERROR)
        printErr(33, "This Problem has a Solution whose limit is Infinite");
    else if(exit_state == SIMPLEXMETHOD_FARBFS_ERROR)
        printErr(33, "No convergence after %hu iterations!", access(curLayout)->max_simplex_iterations);
    else if(exit_state == SIMPLEXMETHOD_ALLOC_ERROR)
        printErr(12, "Simplex Method Heap Dynamic Memory Allocation Problem");
    else
    {
        printf2(COLOR_USER, "\nRelaxed Problem BFS with Artificial Variables is: ");
        printMatrix(stdout, bfs, (dim_typ2){1,dim[RAWS]+dim[COLUMNS]-2});
    }

    matrixFree(&tableau, dim[RAWS]);
    matrixFree(&bfs, 1);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}


__MSSHELL_WRAPPER_ static void _MS__private newtonDifferenceTables(const register sel_typ argc, char ** argv)
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

    newtonDifferenceTable(n, y, FORWARD_DIFFTAB);
    showNewtonDifferenceTable(n, x, y, FORWARD_DIFFTAB);
    newtonDifferenceTable(n, y, BACKWARD_DIFFTAB);
    showNewtonDifferenceTable(n, x, y, BACKWARD_DIFFTAB);

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private lagrangeInterpolation(const register sel_typ argc, char ** argv)
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

    ityp **xy = NULL;
    dim_typ i, j;

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[RAWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy, rc[RAWS]);
            printUsage(&adv_calc[ADVCALC_LAGRANGEINTERPOLATION]);
            return;
        }
    }
    else
    {
        #if WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter related Matrix filled with Data you want Interpolation to Process,\n");
        printf2(COLOR_CREDITS, "by putting on each raws the %hu X and Y Values.\n\n", dim);

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
                #if WINOS
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
                #if WINOS
                    SetExitButtonState(ENABLED);
                #endif // WINOS
                return;
            }
            printErr(5, "Invalid inserted Value");
        }
    }

    ityp yp = 0;
    ityp dr, nr;

	for(i = 0; i < dim; ++i)
    {

        dr = nr = 1.00;
        for(j = 0; j<dim; ++j)
            if(i!=j)
            {
                nr *= xp - xy[XRAW][j];
                dr *= xy[XRAW][i] - xy[XRAW][j];
            }

        yp += nr/dr*xy[YRAW][i];
    }

    matrixFree(&xy, MAX_DIMENSIONS);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    printf2(COLOR_USER, "\nRequested Y VALUE is: ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, yp);
    printf2(COLOR_USER, ".\n\n");

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private greatestEigenValue(const register sel_typ argc, char ** argv)
{

    ityp ***matrix1 = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    matrix1 = malloc(sizeof(ityp**));
    errMem(matrix1, VSPACE);

    if(argc)
    {
        if((!matrixToken(argv[0], matrix1, dim, &dim[COLUMNS])) || dim[RAWS] != dim[COLUMNS])
        {
            matrixFree(matrix1, dim[RAWS]);
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

    ityp ***matrix2 = NULL;
    matrix2 = malloc(sizeof(ityp**));
    errMem(matrix2, (matrixFree(matrix1, dim[RAWS]), free(matrix1)));

    if(!matrixAlloc(matrix2, (dim_typ2){dim[RAWS], 1}))
    {
        matrixFree(matrix1, dim[RAWS]);
        free(matrix1);
        free(matrix2);
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    ityp eigenValue;

    dim_typ i;

    #pragma omp parallel for
	for(i=0; i<dim[RAWS]; ++i)
        (*matrix2)[i][0] = 1.00;

    ityp ***result = NULL;

    result = malloc(sizeof(ityp**));
    errMem(result, (matrixFree(matrix1, dim[RAWS]), matrixFree(matrix2, dim[RAWS]), free(matrix1), free(matrix2)));

    if(!matrixAlloc(result, (dim_typ2){dim[RAWS], 1}))
    {
        matrixFree(matrix1, dim[RAWS]);
        matrixFree(matrix2, dim[RAWS]);
        free(matrix1);
        free(matrix2);
        free(result);
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    ityp resultVector[dim[RAWS]];


    for( ;; )
    {
        _matrixProduct(matrix1, matrix2, result, (dim_typ3){dim[RAWS], dim[RAWS], 1});
        matrixToVector((*result), (dim_typ2){dim[RAWS], 1}, resultVector, MATRIX_TO_VECTOR);
        eigenValue = MAX(3, resultVector);

        #pragma omp parallel for
		for(i=0; i<dim[RAWS]; ++i)
            (*result)[i][0] /= eigenValue;

        if(isEqualMatrix((*matrix2), (*result), (dim_typ2){dim[RAWS], 1}))
            break;
	
		#pragma omp parallel for
        for(i=0; i<dim[RAWS]; ++i)
            (*matrix2)[i][0] = (*result)[i][0];
    }

    printf2(COLOR_USER, "\nGreatest EIGEN Value is: ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, eigenValue);
    printf2(COLOR_USER, ".\nAn EIGEN Vector is:\n\n");

    printMatrix(stdout, (*result), (dim_typ2){dim[RAWS],1});

    matrixFree(matrix1, dim[RAWS]);
    matrixFree(matrix2, dim[RAWS]);
    matrixFree(result, dim[RAWS]);

    free(matrix1);
    free(matrix2);
    free(result);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private funcIntegration(const register sel_typ argc, char ** argv)
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

    printf2(COLOR_USER, "%s(x) Function Integral Value calculated between: x0 = ", ext_math.funcnames[funcID]);
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, x0);
    printf2(COLOR_USER, " and xN = ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, xn);
    printf2(COLOR_USER, ",\nwith %hu Intervals NUMBER is: %6.*f\n\n", n, access(curLayout)->precision, result);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private straightLineFitting(const register sel_typ argc, char ** argv)
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

    ityp **xy = NULL;

    dim_typ i;

    // we must seek for the BACKTRACKING FEATURE

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[RAWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy, rc[RAWS]);
            printUsage(&adv_calc[ADVCALC_STRAIGHTLINEFITTING]);
            return;
        }
    }
    else
    {

        #if WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter the Matrix filled with EQUATIONS DATA,\n");
        printf2(COLOR_CREDITS, "by putting on each RAWS the %hu X and Y VALUES.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){MAX_DIMENSIONS, dim}))
            return;
    }

    ityp sum_x, sum_xy, sum_x2, sum_y;

    sum_x = sum_xy = sum_x2 = sum_y = 0.00;

	for(i = 0; i < dim; ++i)
    {

        sum_x += xy[XRAW][i];
        sum_y += xy[YRAW][i];
        sum_xy += xy[XRAW][i] * xy[YRAW][i];
        sum_x2 += pow(xy[XRAW][i], 2);
    }

    matrixFree(&xy, MAX_DIMENSIONS);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    ityp a, b;

    b = (dim*sum_xy - sum_x*sum_y)/(dim*sum_x2 - pow(sum_x,2));
    a = (sum_y - b*sum_x)/dim;

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

__MSSHELL_WRAPPER_ static void _MS__private parabolicCurveFitting(const register sel_typ argc, char ** argv)
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

    ityp **xy = NULL;


    // we must seek for the BACKTRACKING FEATURE

    if(argc > 1)
    {
        dim_typ rc[MAX_DIMENSIONS];

        if((!matrixToken(argv[1], &xy, rc, &rc[COLUMNS])) || rc[RAWS] != MAX_DIMENSIONS || rc[COLUMNS] != dim)
        {
            matrixFree(&xy, rc[RAWS]);
            printUsage(&adv_calc[ADVCALC_PARABOLICCURVEFITTING]);
            return;
        }
    }
    else
    {
        #if WINOS
            SetExitButtonState(DISABLED);
        #endif // WINOS

        printf2(COLOR_CREDITS, "\nEnter the Matrix filled with EQUATIONS DATA,\n");
        printf2(COLOR_CREDITS, "by putting on each Raws the %hu X and Y VALUES.\n\n", dim);

        if(!insertNMMatrix(&xy, (dim_typ2){MAX_DIMENSIONS, dim}))
            return;
    }

    ityp **matrix = NULL;

    if(!matrixAlloc(&matrix, (dim_typ2){MAX_ABSTRACT_DIMENSIONS, 4}))
    {
        matrixFree(&xy, MAX_DIMENSIONS);
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    matrix[0][0] = dim;

    dim_typ i;

    #pragma omp parallel for
	for(i = 0; i < dim; ++i)
    {
        matrix[0][1] = (matrix[1][0] += xy[XRAW][i]);
        matrix[0][3] += xy[YRAW][i];
        matrix[0][2] = (matrix[2][0] += pow(xy[XRAW][i], 2));
        matrix[1][2] = (matrix[2][1] += pow(xy[XRAW][i], 3));
        matrix[2][2] += pow(xy[XRAW][i], 4);
        matrix[1][3] += xy[XRAW][i]*xy[YRAW][i];
        matrix[2][3] += pow(xy[XRAW][i], 2) * xy[YRAW][i];
    }

    matrixFree(&xy, MAX_DIMENSIONS);

    dim_typ j;
    dim_typ k;

    #pragma omp parallel for
	for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
        #pragma omp parallel for
		for(j = 0; j < MAX_ABSTRACT_DIMENSIONS; ++j)
            if(i!=j)
            {
                const ityp ratio = matrix[j][i]/matrix[i][i];
                #pragma omp parallel for num_threads(4)
                for(k = 0; k < 4; ++k)
                    matrix[j][k] -= ratio * matrix[i][k];
            }

    #pragma omp parallel for
	for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
    {
        const ityp a = matrix[i][i];
        #pragma omp parallel for num_threads(4)
		for(j = 0; j < 4; ++j)
            matrix[i][j] /= a;
    }

    printf2(COLOR_USER, "\nPARABOLIC CURVE Fitting is:\n");
    PRINTL();

    for(i = 0; i < MAX_ABSTRACT_DIMENSIONS; ++i)
    {
        printf2(COLOR_USER, "%c => ", i+97);
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, matrix[i][3]);
        printf2(COLOR_USER, ";\n");
    }

    matrixFree(&matrix, 4);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    PRINTL();
    PRINTN();

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private linearSystemsSolver(const register sel_typ argc, char ** argv)
{
    dim_typ dim[MAX_DIMENSIONS];
    ityp **matrix = NULL;

    if(argc)
    {

        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[COLUMNS] != dim[RAWS]+1)
        {
            matrixFree(&matrix, dim[RAWS]);
            printUsage(&adv_calc[ADVCALC_LINEARSYSTEMSSOLVER]);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "\nEnter Complete Matrix correspondent to the Linear System you want the Program to solve.\n\n");
        if((!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], false)) || dim[COLUMNS] != dim[RAWS]+1)
        {
            if(dim[COLUMNS] != dim[RAWS]+1)
                printErr(33, "You have to insert an [R X R+1] Matrix");
            return;
        }
    }

    dim_typ i, j, k;
    ityp a, b;

    for(i = 0; i < dim[RAWS]; ++i)
        for(j = 0; j < dim[RAWS]; ++j)
            if(i != j)
            {
                a = matrix[j][i];
                b = matrix[i][i];
                for(k = 0; k < dim[COLUMNS]; ++k)
                    matrix[j][k] -= (a/b) * matrix[i][k];
            }

	#pragma omp parallel for
    for(i = 0; i < dim[RAWS]; ++i)
    {
        const ityp c = matrix[i][i];
        for(j = 0; j < dim[COLUMNS]; ++j)
            matrix[i][j] /= c;
    }

    printf2(COLOR_USER, "Simultaneous Solutions of the given Linear System are:\n\n");
    PRINTL();

    for(i = 0; i < dim[RAWS] ; ++i)
    {
        printf2(COLOR_USER, "%c => ", i+97);
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, matrix[i][dim[RAWS]]);
        printf2(COLOR_USER, ";\n");
    }

    matrixFree(&matrix, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    PRINTL();
    PRINTN();
    return;
}
