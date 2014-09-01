// algebra.c 01/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"

#define CHECK_INVERSE_OPERATIONS() if(isSett(BOOLS_AUTOTURNBACK))access(curLayout)->bools &= ~suite_c.bools[BOOLS_INVERSEOPERATIONS].bmask

#ifdef ALLOW_MATMANAGER
__MSSHELL_WRAPPER_ static void _MS__private __apnt matrixManager(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __apnt matrixManager(const register sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_MATMANAGER_PROGS,
                        mat_manager, alg_operations[ALGEBRA_MATRICESMANAGER].name,
                        #if MAX_MATMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif
                        );
    return;
}
#endif

__MSSHELL_WRAPPER_ static void _MS__private matrixSort(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixNorm(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixDet(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixTrace(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixRank(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixSVD(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixInv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixTranspose(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixSum(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixProduct(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixKProduct(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private perScalarProduct(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private illConditionChecking(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private matrixFattLU(const register sel_typ argc, char ** argv);

sprog alg_operations[MAX_ALGEBRA_OPERATIONS] =
{
    #ifdef ALLOW_MATMANAGER
    [ALGEBRA_MATRICESMANAGER] =
    {
        "Matrices Manager",
        CMD_MATRIXMANAGER,
        USAGE_MATRIXMANAGER,
        matrixManager,
        AUTOMATIC,
        FATHER
    },
    #endif
    [ALGEBRA_MATRIXSORT] =
    {
        "Matrix Sort",
        CMD_MATRIXSORT,
        USAGE_MATRIXSORT,
        matrixSort,
        BY_USER,
        CHILD
    },
    [ALGEBRA_NORMCALCULATOR] =
    {
        "Square Matrix Norm Calculator",
        CMD_MATRIXNORM,
        USAGE_MATRIXNORM,
        matrixNorm,
        BY_USER,
        CHILD
    },
    [ALGEBRA_DETERMINANTCALCULATOR] =
    {
        "Square Matrix Determinant Calculator",
        CMD_MATRIXDET,
        USAGE_MATRIXDET,
        matrixDet,
        BY_USER,
        CHILD
    },
    [ALGEBRA_TRACECALCULATOR] =
    {
        "Square Matrix Trace Calculator",
        CMD_MATRIXTRACE,
        USAGE_MATRIXTRACE,
        matrixTrace,
        BY_USER,
        CHILD
    },
    [ALGEBRA_RANKCALCULATOR] =
    {
        "Matrix Rank Calculator",
        CMD_MATRIXRANK,
        USAGE_MATRIXRANK,
        matrixRank,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRIXSVD] =
    {
        "Matrix SVD",
        CMD_MATRIXSVD,
        USAGE_MATRIXSVD,
        matrixSVD,
        BY_USER,
        CHILD
    },
    [ALGEBRA_INVERSEMATRIX] =
    {
        "Square Matrix Inversion",
        CMD_MATRIXINV,
        USAGE_MATRIXINV,
        matrixInv,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRIXTRANSPOSE] =
    {
        "Matrix Transposition",
        CMD_MATRIXTRANSPOSE,
        USAGE_MATRIXTRANSPOSE,
        matrixTranspose,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRICESSUM] =
    {
        "Matrices Sum",
        CMD_MATRIXSUM,
        USAGE_MATRIXSUM,
        matrixSum,
        BY_USER,
        CHILD
    },
    [ALGEBRA_TENSORSSUM] =
    {
        "Tensors Sum",
        CMD_TENSORSUM,
        USAGE_TENSORSUM,
        matrixSum,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRICESPRODUCT] =
    {
        "Matrices Product",
        CMD_MATRIXPRODUCT,
        USAGE_MATRIXPRODUCT,
        matrixProduct,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRICESKPRODUCT] =
    {
        "Matrices Kronecker Product",
        CMD_MATRIXKPRODUCT,
        USAGE_MATRIXKPRODUCT,
        matrixKProduct,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRIXPOWER] =
    {
        "Square Matrix N Power",
        CMD_MATRIXPOWER,
        USAGE_MATRIXPOWER,
        matrixProduct,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRIXKPOWER] =
    {
        "Matrix Kronecker N Power",
        CMD_MATRIXKPOWER,
        USAGE_MATRIXKPOWER,
        matrixKProduct,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRIXPERVECTOR] =
    {
        "Matrix per Vector Product",
        CMD_MATRIXPERVECTOR,
        USAGE_MATRIXPERVECTOR,
        matrixProduct,
        BY_USER,
        CHILD
    },
    [ALGEBRA_SCALARPRODUCT] =
    {
        "Vectors Scalar Product",
        CMD_SCALARPRODUCT,
        USAGE_SCALARPRODUCT,
        matrixProduct,
        BY_USER,
        CHILD
    },
    [ALGEBRA_PERSCALARPRODUCT] =
    {
        "Matrix per Scalar Product",
        CMD_PERSCALARPRODUCT,
        USAGE_PERSCALARPRODUCT,
        perScalarProduct,
        BY_USER,
        CHILD
    },
    [ALGEBRA_SCALARDIVISIONMATRIX] =
    {
        "Scalar Division Matrix",
        CMD_SCALARDIVISIONMATRIX,
        USAGE_SCALARDIVISIONMATRIX,
        perScalarProduct,
        BY_USER,
        CHILD
    },
    [ALGEBRA_ILLCONDITIONCHECKING] =
    {
        "Square Matrix Ill Condition Checking",
        CMD_ILLCONDITIONCHECKING,
        USAGE_ILLCONDITIONCHECKING,
        illConditionChecking,
        BY_USER,
        CHILD
    },
    [ALGEBRA_MATRIXLUFACTORIZATION] =
    {
        "Square Matrix LU-Decomposition",
        CMD_MATRIXFATTLU,
        USAGE_MATRIXFATTLU,
        matrixFattLU,
        BY_USER,
        CHILD
    }
};

#define PARSING_MATRIX_ALLOWED isSett(BOOLS_MATRIXPARSING)

__MSSHELL_WRAPPER_ static void _MS__private matrixSort(const register sel_typ argc, char ** argv)
{

    sel_typ tmp;
    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {
        if(argc == MAX_DIMENSIONS)
        {

            if((tmp = strtod(argv[0], NULL)) == MAX_ABSTRACT_DIMENSIONS) return;
            if(tmp < 0 || tmp > MAX_ABSTRACT_DIMENSIONS)
            {
                printUsage(&alg_operations[ALGEBRA_MATRIXSORT]);
                return;
            }

            if(!tmp)
            {
                if(!matrixToken(argv[1], &matrix, dim, &dim[COLUMNS]))
                    return;
            }
            else
            {
                if(!insertDims(&dim[RAWS], &dim[COLUMNS]))
                    return;

                if(!matrixAlloc(&matrix, dim))
                    return;

                if(!randomMatrix(matrix, dim))
                {
                    matrixFree(&matrix, dim[RAWS]);
                    return;
                }
            }

        }
        else
        {
            printUsage(&alg_operations[ALGEBRA_MATRIXSORT]);
            return;
        }
    }
    else
    {

        for( ; ; )
        {

            CLEARBUFFER();

            printf2(COLOR_CREDITS, "\nSelect Matrix Filling Mode.\n");
            printf2(COLOR_CREDITS, "A for User Matrix Inserting, B for Random Matrix Generation;\n");
            printf2(COLOR_CREDITS, "- %c to go Back..\n\n", access(curLayout)->exit_char);

            PRINTL();

            do if((tmp = toupper(getch())) == access(curLayout)->exit_char) return;
            while(tmp < 'A' && tmp > access(curLayout)->exit_char);

            if(tmp == 'A')
            {
                if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], false))
                    continue;
            }

            if(tmp == 'B')
            {
                if(!insertDims(&dim[RAWS], &dim[COLUMNS]))
                    continue;
                if(!matrixAlloc(&matrix, dim))
                    continue;
                if(!randomMatrix(matrix, dim))
                {
                    matrixFree(&matrix, dim[RAWS]);
                    continue;
                }
            }

            break;
        }
        tmp -= 'A';
    }

    CLEARBUFFER();

    // PRINTL();


    PRINTL();

    ityp vector[dim[RAWS]*dim[COLUMNS]];

    matrixToVector(matrix, dim, vector, MATRIX_TO_VECTOR);


    qsort(vector,dim[RAWS]*dim[COLUMNS], sizeof(ityp), cmpfunc);

    vectorToMatrix(dim, vector, matrix);

    printf2(COLOR_USER, "\nMatrix has been sort with Ascending Order:\n\n");
    PRINTL();
    printMatrix(stdout, matrix, dim);
    matrixFree(&matrix, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixNorm(const register sel_typ argc, char ** argv)
{
    sel_typ tmp;

    if(argc)
    {
        if((tmp = strtod(argv[0], NULL)) == MAX_ABSTRACT_DIMENSIONS) return;
        if(tmp < 0 || tmp > MAX_ABSTRACT_DIMENSIONS)
        {
            printUsage(&alg_operations[ALGEBRA_NORMCALCULATOR]);
            return;
        }
    }
    else
    {
        printf("\nSelect Norm Calculation Mode:\n");
        printf("- A for P-Norm, B for 1-Norm, C for Inf-Norm;\n");
        printf("- D to go Back...\n\n");

        PRINTL();

        do if((tmp = toupper(getch())) == 'D') return;
        while(tmp < 'A' && tmp > 'D');

        tmp -= 'A';
    }

    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc > 1)
    {
        if((!matrixToken(argv[1], &matrix, dim, &dim[COLUMNS])) || dim[RAWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim[RAWS]);
            printUsage(&alg_operations[ALGEBRA_NORMCALCULATOR]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], true))
        return;

    const register ityp nrm = tmp ? norm(matrix, dim[RAWS], tmp-1) : norms(matrix, dim[RAWS]);

    matrixFree(&matrix, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        *(access(exprVars)->e_ANS) = nrm;

    printf2(COLOR_USER, "\nInserted Quad Matrix NORM is: ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, nrm);
    printf2(COLOR_USER, ".\n\n");

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixDet(const register sel_typ argc, char ** argv)
{
    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {

        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[RAWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim[RAWS]);
            printUsage(&alg_operations[ALGEBRA_DETERMINANTCALCULATOR]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], true))
        return;

    bool flag = false;
    const register ityp dt = det(matrix, dim[RAWS], &flag);

    matrixFree(&matrix, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        *(access(exprVars)->e_ANS) = dt;

    printf2(COLOR_USER, "\nInserted Quad Matrix %sDETERMINANT is: ", flag ? "ABS-":NULL_CHAR);
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, dt);
    printf2(COLOR_USER, ".\n\n");

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixTrace(const register sel_typ argc, char ** argv)
{
    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];


    if(argc)
    {
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[RAWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim[RAWS]);
            printUsage(&alg_operations[ALGEBRA_TRACECALCULATOR]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], true))
        return;

    const register ityp trc = trace(matrix, dim[RAWS]);

    matrixFree(&matrix, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        *(access(exprVars)->e_ANS) = trc;

    printf2(COLOR_USER, "\nInserted Quad Matrix TRACE is: ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, trc);
    printf2(COLOR_USER, ".\n\n");

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixRank(const register sel_typ argc, char ** argv)
{
    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {
        if(!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS]))
            return;
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], false))
        return;
    const register dim_typ rk = rank(matrix, dim);

    matrixFree(&matrix, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        *(access(exprVars)->e_ANS) = rk;

    printf2(COLOR_USER, "\nInserted Matrix RANK is: %hu.\n\n", rk);

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixSVD(const register sel_typ argc, char ** argv)
{
    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {
        if(!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS]))
            return;
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], false))
        return;

    ityp * S = NULL;
    ityp ** V = NULL;

    const register dim_typ columns = dim[COLUMNS];

    S = malloc(sizeof(ityp)*columns);
    errMem(S, VSPACE);

    if(!matrixAlloc(&V, (dim_typ2){columns, columns}))
    {
        free(S);
        matrixFree(&V, columns);
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }


    dsvd(matrix, dim, S, V);
    matrixFree(&V, columns);
    matrixFree(&matrix, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    printf2(COLOR_USER, "\nSINGULAR VALUES Vector is:\n");

    SHOWPAUSEMESSAGE();

    for(uint64_t i=0; i<columns; ++i)
    {
        printf2(COLOR_USER, "- %llu: ", i);
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, S[i]);
        PRINTN();
        if(catchPause()) return;
    }

    free(S);

    PRINTL();
    PRINT2N();

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixInv(const register sel_typ argc, char ** argv)
{
    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[RAWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim[RAWS]);
            printUsage(&alg_operations[ALGEBRA_INVERSEMATRIX]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], true))
        return;

    const dim_typ dimraws_per2 = dim[RAWS]<<1;

    matrix = realloc(matrix, sizeof(ityp *)*dimraws_per2);
    errMem(matrix, (matrixFree(&matrix, dimraws_per2)));

    dim_typ i;

    for(i=0; i<dim[RAWS]; ++i)
    {
        matrix[i] = realloc(matrix[i], sizeof(ityp)*dimraws_per2);
        errMem(matrix[i], (matrixFree(&matrix, dimraws_per2), free_foreach2(&matrix, dimraws_per2)));
    }


    if(!invertMatrix(matrix, dim[RAWS]))
    {
        printErr(1, "You cannot invert SINGULAR Matrices");
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    printf2(COLOR_SYSTEM, "\nINVERSE MATRIX of inserted Quad Matrix is:\n\n");

    dim_typ j;

    // VISUALIZZAZIONE PARTICOLARE MATRICE
    for(i=0; i<dim[RAWS]; ++i)
    {
        for(j = dim[RAWS]; j < dimraws_per2; ++j)
        {
            printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, matrix[i][j]);
            printf2(COLOR_USER, j == dimraws_per2-1 ? ";" : ",");
            PRINTSPACE();
        }
        PRINTN();
    }

    matrixFree(&matrix, dimraws_per2);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixTranspose(const register sel_typ argc, char ** argv)
{
    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {
        if(!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS]))
            return;
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], false))
        return;

    printf2(COLOR_SYSTEM, "\nTRANSPOSED MATRIX of Inserted Matrix is:\n\n");

    ityp **matrix2 = NULL;

    if(!matrixAlloc(&matrix2, (dim_typ2){dim[COLUMNS], dim[RAWS]}))
    {
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    transpose(matrix, matrix2, dim); // dim[RAWS], dim[COLUMNS]);
    printMatrix(stdout, matrix2, (dim_typ2){dim[COLUMNS], dim[RAWS]});

    matrixFree(&matrix, dim[RAWS]);
    matrixFree(&matrix2, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixSum(const register sel_typ argc, char ** argv)
{
    dim_typ i, j, k, l;
    dim_typ dim[MAX_DIMENSIONS];
    dim_typ cdim[MAX_DIMENSIONS];
    uint64_t tdim=1;

    ityp ****matrix1 = NULL;
    ityp ****matrix2 = NULL;
    ityp ****matrix_sum = NULL;
    const bool complex_entries = access(curLayout)->algebra != 0;
    const bool tensor_mode = __pmode__ == ALGEBRA_TENSORSSUM;
    const sel_typ algebra_units = exp2(access(curLayout)->algebra);

    ityp tmp;

    if(tensor_mode)
    {
        tmp = 0.00;
        if(argc)
        {
            if(PARSING_MATRIX_ALLOWED)
            {
                if((!parse(argv[0], &tmp)) || tmp != (tdim = (uint64_t)tmp) || tdim < 1)
                {
                    printUsage(&alg_operations[__pmode__]);
                    return;
                }
            }
            else if((tmp = strtod(argv[0], NULL)) != (tdim=(uint64_t)tmp) || tdim < 1)
            {
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {

            printf2(COLOR_CREDITS, "\nEnter Tensor Dimension.\n");

            if(PARSING_SYSTEM_ALLOWED)
                PRINTHOWTOBACKMESSAGE();

            PRINT2N();

            while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
                    (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (tdim = (uint64_t)tmp) || tdim < 1)
            {
                CLEARBUFFER();
                if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
                if(isNullVal(tmp) && exitHandleCheck) return;
                printErr(1, "Invalid inserted Value.\nTensor Dimension must be a non-negative and >= 1 integer");
            }

        }

    }

    char matrix_type[MINMIN_STRING];
    strcpy(matrix_type, INVERSE_OPS ? "Difference":"Sum");

    matrix1 = malloc(sizeof(ityp***)*tdim);
    errMem(matrix1, VSPACE);

    for(i=0; i<tdim; ++i)
    {
        matrix1[i] = malloc(sizeof(ityp**)*algebra_units);
        if(checkErrMem(matrix1[i]))
        {
        	#pragma omp parallel for
            for(j=0; j<i; ++j)
                free(matrix1[j]);
            free(matrix1);
            return;
        }
    }

    sel_typ argv_pos = tensor_mode;

    for(k=0; k<tdim; ++k, ++argv_pos)
    {
        if(argv_pos)
            printf2(COLOR_CREDITS, "Performing the %hu Matrices Component %s Operation.\n", k+1, matrix_type);

        if(argc > argv_pos)
        {
            if(!matrixToken(argv[argv_pos], matrix1[k], dim, &dim[COLUMNS]))
            {
            	if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
		       		#pragma omp parallel for num_threads(tdim)
		            for(i=0; i<tdim; ++i)
		                free(matrix1[i]);
		        else
		        	#pragma omp parallel for
		            for(i=0; i<tdim; ++i)
		                free(matrix1[i]);
                free(matrix1);
                return;
            }
            #if WINOS
                SetExitButtonState(DISABLED);
            #endif // WINOS
        }
        else if(!enterMatrix(matrix1[k], dim, &dim[COLUMNS], false, !complex_entries))
        {
        	if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	        	#pragma omp parallel for num_threads(tdim)
	            for(i=0; i<tdim; ++i)
	                free(matrix1[i]);
	        else
	        	#pragma omp parallel for
	            for(i=0; i<tdim; ++i)
	                free(matrix1[i]);
            free(matrix1);
            return;
        }

        // Entering eventually IMAGINARY PART MATRIX of COMPLEX MATRIX1
        for(i=1; i<algebra_units; ++i, ++argv_pos)
        {
            matrix1[k][i] = NULL;
            printf2(COLOR_CREDITS, "Enter IMAGINARY %s PART of the First Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printf2(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[RAWS], dim[COLUMNS]);
            if(argc > argv_pos)
            {
                if((!matrixToken(argv[argv_pos], &matrix1[k][i], cdim, &cdim[COLUMNS])) || cdim[RAWS] != dim[RAWS] || cdim[COLUMNS] != dim[COLUMNS])
                {
                	#pragma omp parallel for
                    for(j=0; j<=i; ++j)
                        matrixFree(&matrix1[k][j], cdim[RAWS]);
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
	                    for(j=0; j<algebra_units; ++j)
	                        matrixFree(&matrix1[k][j], cdim[RAWS]);
                    else
                    	#pragma omp parallel for
                    	for(j=0; j<algebra_units; ++j)
	                        matrixFree(&matrix1[k][j], cdim[RAWS]);
	                    
	                if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
                    	#pragma omp parallel for num_threads(tdim)
	                    for(j=0; j<tdim; ++j)
	                        free(matrix1[j]);
	                else
	                	#pragma omp parallel for
	                    for(j=0; j<tdim; ++j)
	                        free(matrix1[j]);
	                    
                    free(matrix1);
                    printUsage(&alg_operations[__pmode__]);
                    return;
                }
            }
            else
            {
                if(!insertNMMatrix(&matrix1[k][i], dim))
                {
                	
                	#pragma omp parallel for
                    for(j=0; j<i; ++j)
                        matrixFree(&matrix1[k][j], dim[RAWS]);
                    if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
	                    for(j=0; j<algebra_units; ++j)
	                        matrixFree(&matrix1[k][j], dim[RAWS]);
	            	else
	            		#pragma omp parallel for
	                    for(j=0; j<algebra_units; ++j)
	                        matrixFree(&matrix1[k][j], dim[RAWS]);
	                
	                if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(tdim)
						for(j=0; j<tdim; ++j)
	                        free(matrix1[j]);
	                else
	                	#pragma omp parallel for
						for(j=0; j<tdim; ++j)
	                        free(matrix1[j]);

                    free(matrix1);
                    return;
                }
            }
        }

        if(complex_entries)
        {
            printf2(COLOR_USER, "\n\nFirst Complex Matrix correctly inserted:\n\n");
            printMatrix(stdout, matrix1[k][0], dim);
            for(i=1; i<algebra_units; ++i)
            {
                printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
                printMatrix(stdout, matrix1[k][i], dim);
            }
        }

        if(!k)
        {
            matrix2 = malloc(sizeof(ityp***)*tdim);
            if(checkErrMem(matrix2))
            {
                if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	            	#pragma omp parallel for num_threads(algebra_units)
	                for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[k][j], dim[RAWS]);
	            else
	            	#pragma omp parallel for
	                for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[k][j], dim[RAWS]);
	            
	        	if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
				    #pragma omp parallel for num_threads(tdim)
					for(i=0; i<tdim; ++i)
	                    free(matrix1[i]);
	            else
	            	#pragma omp parallel for
					for(i=0; i<tdim; ++i)
	                    free(matrix1[i]);
                free(matrix1);
            }

			for(i=0; i<tdim; ++i)
            {
                matrix2[i] = malloc(sizeof(ityp**)*algebra_units);
                if(checkErrMem(matrix2[i]))
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(j=0; j<algebra_units; ++j)
	                        matrixFree(&matrix1[k][j], dim[RAWS]);
	            	else
	            		#pragma omp parallel for
						for(j=0; j<algebra_units; ++j)
	                        matrixFree(&matrix1[k][j], dim[RAWS]);
	                
	                if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(tdim)
						for(j=0; j<tdim; ++j)
	                        free(matrix1[j]);
	                else
	                	#pragma omp parallel for
						for(j=0; j<tdim; ++j)
	                        free(matrix1[j]);
	                        
                    #pragma omp parallel for
				    for(j=0; j<i; ++j)
                        free(matrix2[j]);
                    free(matrix1);
                    free(matrix2);
                    return;
                }
            }
        }

        if(argc > ++argv_pos)
        {
            dim_typ rc[MAX_DIMENSIONS];

            if((!matrixToken(argv[argv_pos], matrix2[k], rc, &rc[COLUMNS])) || rc[RAWS] != dim[RAWS] || rc[COLUMNS] != dim[COLUMNS])
            {
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	            	#pragma omp parallel for num_threads(algebra_units)
	                for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[k][i], dim[RAWS]);
	            else
	            	#pragma omp parallel for
	                for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[k][i], dim[RAWS]);
	                    
                matrixFree(matrix2[k], rc[RAWS]);
                
                if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(tdim)
					for(i=0; i<tdim; ++i)
	                {
	                    free(matrix1[i]);
	                    free(matrix2[i]);
	                }
	            else
	            	#pragma omp parallel for
					for(i=0; i<tdim; ++i)
	                {
	                    free(matrix1[i]);
	                    free(matrix2[i]);
	                }
	                
                free(matrix1);
                free(matrix2);
                #if WINOS
                    SetExitButtonState(ENABLED);
                #endif // WINOS
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {

            printf2(COLOR_CREDITS, "Enter second [%hu x %hu] Matrix Elements.\n", dim[RAWS], dim[COLUMNS]);
            printf2(COLOR_CREDITS, "Elements Inserting Carriage will be automatically redirected on second scansion (for RAWS Scanning).\n\n");

            if(!insertNMMatrix(matrix2[k], dim)) // dim[RAWS], dim[COLUMNS]))
            {
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
				    for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[k][i], dim[RAWS]);
	            else
	            	#pragma omp parallel for
				    for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[k][i], dim[RAWS]);
	            
	            if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(tdim)
					for(i=0; i<tdim; ++i)
	                {
	                    free(matrix1[i]);
	                    free(matrix2[i]);
	                }
	            else
		            #pragma omp parallel for
						for(i=0; i<tdim; ++i)
		                {
		                    free(matrix1[i]);
		                    free(matrix2[i]);
		                }
		                
                free(matrix1);
                free(matrix2);
                return;
            }
        }

        for(i=1; i<algebra_units; ++i, ++argv_pos)
        {
            matrix2[k][i] = NULL;
            printf2(COLOR_CREDITS, "Enter IMAGINARY %s PART of the Second Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printf2(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[RAWS], dim[COLUMNS]);
            if(argc > argv_pos)
            {

                if((!matrixToken(argv[argv_pos], &matrix2[k][i], cdim, &cdim[COLUMNS])) || cdim[RAWS] != dim[RAWS] || cdim[COLUMNS] != dim[COLUMNS])
                {
                    #pragma omp parallel for
					for(j=0; j<=i; ++j)
                    {
                        matrixFree(&matrix1[k][j], dim[RAWS]);
                        matrixFree(&matrix2[k][j], cdim[RAWS]);
                    }
                    if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(j=0; j<algebra_units; ++j)
	                    {
	                        matrixFree(&matrix1[k][j], dim[RAWS]);
	                        matrixFree(&matrix2[k][j], cdim[RAWS]);
	                    }
	                else
	                	#pragma omp parallel for
						for(j=0; j<algebra_units; ++j)
	                    {
	                        matrixFree(&matrix1[k][j], dim[RAWS]);
	                        matrixFree(&matrix2[k][j], cdim[RAWS]);
	                    }
	                    
	                if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(tdim)
					    for(j=0; j<tdim; ++j)
	                    {
	                        free(matrix1[j]);
	                        free(matrix2[j]);
	                    }
	                else
		                #pragma omp parallel for
						    for(j=0; j<tdim; ++j)
		                    {
		                        free(matrix1[j]);
		                        free(matrix2[j]);
		                    }
                    free(matrix1);
                    free(matrix2);
                    printUsage(&alg_operations[__pmode__]);
                    return;
                }

            }
            else
            {
                if(!insertNMMatrix(&matrix2[k][i], (dim_typ2){dim[RAWS], dim[COLUMNS]}))
                {
                	#pragma omp parallel for
                    for(j=0; j<i; ++j)
                        matrixFree(&matrix1[k][j], dim[RAWS]);
                    if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(j=0; j<algebra_units; ++j)
	                        matrixFree(&matrix1[k][j], dim[RAWS]);
	                else
	                	#pragma omp parallel for
						for(j=0; j<algebra_units; ++j)
	                        matrixFree(&matrix1[k][j], dim[RAWS]);
	                        
                    matrixFree(matrix2[k], dim[RAWS]);
                	if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(tdim)
						for(j=0; j<tdim; ++j)
	                    {
	                        free(matrix1[j]);
	                        free(matrix2[j]);
	                    }
	                else
	                	#pragma omp parallel for
						for(j=0; j<tdim; ++j)
	                    {
	                        free(matrix1[j]);
	                        free(matrix2[j]);
	                    }
                    free(matrix1);
                    free(matrix2);
                    return;
                }
            }
        }

        printf2(COLOR_USER, "\nSecond Matrix correctly inserted:\n\n");
        printMatrix(stdout, matrix2[k][0], dim); // dim[RAWS], dim[COLUMNS]);

        for(i=1; i<algebra_units; ++i)
        {
            printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printMatrix(stdout, matrix2[k][i], dim);
        }

        printf2(COLOR_SYSTEM, "\nMatrix %s [%hu x %hu] between First and Second Matrix is:\n", matrix_type, dim[RAWS], dim[COLUMNS]);

        if(!k)
        {
            matrix_sum = malloc(sizeof(ityp***)*tdim);
            if(checkErrMem(matrix_sum))
            {
                #pragma omp parallel for num_threads(tdim)
				for(i=k; i<tdim; ++i)
                    #pragma omp parallel for num_threads(algebra_units)
					for(j=0; j<algebra_units; ++j)
                    {
                        matrixFree(&matrix1[i][j], dim[RAWS]);
                        matrixFree(&matrix2[i][j], dim[RAWS]);
                    }
                if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(tdim)
					for(i=0; i<tdim; ++i)
	                {
	                    free(matrix1[i]);
	                    free(matrix2[i]);
	                }
	            else
	            	#pragma omp parallel for
					for(i=0; i<tdim; ++i)
	                {
	                    free(matrix1[i]);
	                    free(matrix2[i]);
	                }
                free(matrix1);
                free(matrix2);
            }

			for(i=0; i<tdim; ++i)
            {
                matrix_sum[i] = malloc(sizeof(ityp**)*algebra_units);
                if(checkErrMem(matrix_sum[i]))
                {
                    #pragma omp parallel for num_threads(tdim)
					for(l=k; l<tdim; ++l)
                        #pragma omp parallel for num_threads(algebra_units)
						for(j=0; j<algebra_units; ++j)
                        {
                            matrixFree(&matrix1[l][j], dim[RAWS]);
                            matrixFree(&matrix2[l][j], dim[RAWS]);
                        }

					if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(tdim)
						for(j=0; j<tdim; ++j)
	                    {
	                        free(matrix1[j]);
	                        free(matrix2[j]);
	                    }
	                else
	                	#pragma omp parallel for
						for(j=0; j<tdim; ++j)
	                    {
	                        free(matrix1[j]);
	                        free(matrix2[j]);
	                    }
                    #pragma omp parallel for
					for(j=0; j<i; ++j)
                        free(matrix_sum[j]);
                    free(matrix1);
                    free(matrix2);
                    free(matrix_sum);
                    return;
                }
            }
        }

        bool mxsumErr = false;

		for(i=0; i<algebra_units; ++i)
            if((mxsumErr = !matrixAlloc(&matrix_sum[k][i], dim)))
                break;

        if(mxsumErr)
        {
            #pragma omp parallel for num_threads(tdim)
			for(i=k; i<tdim; ++i)
                #pragma omp parallel for num_threads(algebra_units)
				for(j=0; j<algebra_units; ++j)
                {
                    matrixFree(&matrix1[i][j], dim[RAWS]);
                    matrixFree(&matrix2[i][j], dim[RAWS]);
                }
            if(tdim <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && tdim > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	            #pragma omp parallel for num_threads(tdim)
				for(i=0; i<tdim; ++i)
	            {
	                free(matrix1[i]);
	                free(matrix2[i]);
	                free(matrix_sum[i]);
	            }
	        else
	        	#pragma omp parallel for
				for(i=0; i<tdim; ++i)
	            {
	                free(matrix1[i]);
	                free(matrix2[i]);
	                free(matrix_sum[i]);
	            }
            free(matrix1);
            free(matrix2);
            free(matrix_sum);
            #if WINOS
                SetExitButtonState(ENABLED);
            #endif // WINOS
            return;
        }

        access(matrixSumFunc)(matrix1[k], matrix2[k], matrix_sum[k], dim);

        CHECK_INVERSE_OPERATIONS();

        PRINTL();
        printMatrix(stdout, matrix_sum[k][0], dim);
        for(i=1; i<algebra_units; ++i)
        {
            printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printMatrix(stdout, matrix_sum[k][i], dim);
            matrixFree(&matrix1[k][i], dim[RAWS]);
            matrixFree(&matrix2[k][i], dim[RAWS]);
            matrixFree(&matrix_sum[k][i], dim[RAWS]);
        }
        // Freeing Matrix

        matrixFree(matrix1[k], dim[RAWS]);
        matrixFree(matrix2[k], dim[RAWS]);
        matrixFree(matrix_sum[k], dim[RAWS]);

        free(matrix1[k]);
        free(matrix2[k]);
        free(matrix_sum[k]);
    }

    free(matrix1);
    free(matrix2);
    free(matrix_sum);


    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}


__MSSHELL_WRAPPER_ static void _MS__private matrixProduct(const register sel_typ argc, char ** argv)
{
    dim_typ dim[MAX_ABSTRACT_DIMENSIONS];
    dim_typ cdim[MAX_DIMENSIONS];
    dim_typ ii, ij;
    char tmp;


    ityp ***matrix1 = NULL;
    const bool assert_m = __pmode__ == ALGEBRA_MATRIXPOWER;
    const bool complex_entries = access(curLayout)->algebra != 0;
    const sel_typ algebra_units = exp2(access(curLayout)->algebra);

    matrix1 = malloc(sizeof(ityp**)*algebra_units);
    errMem(matrix1, VSPACE);

    // Richiesta di inserimento valori Prima Matrice.

    if(argc)
    {
        if((!matrixToken(argv[0], matrix1, dim, &dim[COLUMNS])) || (assert_m && dim[RAWS] != dim[COLUMNS]))
        {
            matrixFree(matrix1, dim[RAWS]);
            free(matrix1);
            printUsage(&alg_operations[__pmode__]);
            return;
        }
    }
    else if(!enterMatrix(matrix1, dim, &dim[COLUMNS], assert_m, !complex_entries))
    {
        free(matrix1);
        return;
    }

    sel_typ argv_pos = 0;


    // Entering eventually IMAGINARY PART MATRIX of COMPLEX MATRIX1
    for(ii=1; ii<algebra_units; ++ii, ++argv_pos)
    {
        matrix1[ii] = NULL;
        printf2(COLOR_CREDITS, "Enter IMAGINARY %s PART of the First Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
        printf2(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[RAWS], dim[COLUMNS]);
        if(argc > argv_pos)
        {
            if((!matrixToken(argv[argv_pos], &matrix1[ii], cdim, &cdim[COLUMNS])) || cdim[RAWS] != dim[RAWS] || cdim[COLUMNS] != dim[COLUMNS])
            {
                #pragma omp parallel for
				for(ij=0; ij<=ii; ++ij)
                    matrixFree(&matrix1[ij], cdim[RAWS]);
                free(matrix1);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {
            if(!insertNMMatrix(&matrix1[ii], dim))
            {
                #pragma omp parallel for
				for(ij=0; ij<ii; ++ij)
                    matrixFree(&matrix1[ij], dim[RAWS]);
                free(matrix1);
                return;
            }
        }
    }

    if(complex_entries)
    {
        printf2(COLOR_USER, "\n\nFirst Complex Matrix correctly inserted:\n\n");
        printMatrix(stdout, (*matrix1), dim);
        for(ii=1; ii<algebra_units; ++ii)
        {
            printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
            printMatrix(stdout, matrix1[ii], dim);
        }
    }


    // Richiesta di inserimento valori Seconda Matrice.

    tmp = 1;

    uint64_t tmp2 = 0;
    ityp ***matrix2 = NULL;

    matrix2 = malloc(sizeof(ityp**)*algebra_units);
    errMem(matrix2, (free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE), free(matrix1)));

    ++ argv_pos;

    if(assert_m)
    {
        ityp tmp3 = 0.00;

        if(argc > argv_pos)
        {
            if(PARSING_MATRIX_ALLOWED)
            {
                if((!parse(argv[argv_pos], &tmp3)) || tmp3 != (tmp2 = (uint64_t)tmp3) || tmp2 < 1)
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
	                else
	                	#pragma omp parallel for
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
	                        
                    free(matrix1);
                    free(matrix2);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif // WINOS
                    printUsage(&alg_operations[__pmode__]);
                    return;
                }
            }
            else if((tmp3 = strtod(argv[argv_pos], NULL)) != (tmp2 = (uint64_t)tmp3) || tmp2 < 1)
            {
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(ii=0; ii<algebra_units; ++ii)
	                    matrixFree(&matrix1[ii], dim[RAWS]);
	            else
	            	#pragma omp parallel for
					for(ii=0; ii<algebra_units; ++ii)
	                    matrixFree(&matrix1[ii], dim[RAWS]);
                free(matrix1);
                free(matrix2);
                #if WINOS
                    SetExitButtonState(ENABLED);
                #endif // WINOS
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {
            printf2(COLOR_CREDITS, "\nEnter N Power to which you want to raise the Inserted Quad Matrix.\n");

            if(PARSING_SYSTEM_ALLOWED)
                PRINTHOWTOBACKMESSAGE();

            PRINT2N();

            while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp3 = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
                (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp3))) || tmp3 != (tmp2 = (uint64_t)tmp3) || tmp2 < 1)
            {
                CLEARBUFFER();
                if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
                if(isNullVal(tmp3) && exitHandleCheck)
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
	                else
	                	#pragma omp parallel for
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
                    free(matrix1);
                    free(matrix2);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif // WINOS
                    return;
                }
                printErr(1, "Invalid inserted Value.\nExponent must be a non-negative and >= 1 integer");
            }

            CLEARBUFFER();

            if(tmp2 != 1)
            {
                if(!matrixAlloc(matrix2, dim))
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
	                else
	                	#pragma omp parallel for 
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
                    free(matrix1);
                    free(matrix2);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif
                    return;
                }
                equalMatrix(matrix2, (*matrix1), dim);
                if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(ii=1; ii<algebra_units; ++ii)
	                {
	                    (*matrix2)[ii] = NULL;
	                    equalMatrix(&matrix2[ii], matrix1[ii], dim);
	                }
	            else
		            #pragma omp parallel for
						for(ii=1; ii<algebra_units; ++ii)
		                {
		                    (*matrix2)[ii] = NULL;
		                    equalMatrix(&matrix2[ii], matrix1[ii], dim);
		                }
	            	
            }

            dim[COLUMNS2] = dim[RAWS];
        }
    }
    else
    {
        const bool assert = (__pmode__ == ALGEBRA_SCALARPRODUCT || __pmode__ == ALGEBRA_MATRIXPERVECTOR);

        dim_typ i, j;
        dim_typ start_col_index;
        volatile sel_typ back_tracking;

        start_col_index = 0;

        if(argc > argv_pos)
        {
            dim_typ raw2;
            if((!matrixToken(argv[argv_pos], matrix2, &raw2, &dim[COLUMNS2])) || raw2 != dim[COLUMNS])
            {
                matrixFree(matrix2, raw2);
                if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(ii=0; ii<algebra_units; ++ii)
	                    matrixFree(&matrix1[ii], dim[RAWS]);
	            else
	            	#pragma omp parallel for
					for(ii=0; ii<algebra_units; ++ii)
	                    matrixFree(&matrix1[ii], dim[RAWS]);
                free(matrix1);
                free(matrix2);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }// Overhead here isn't inevitable because the order is important
        else
            if(isSett(BOOLS_INSERTMODE))
            {
                printf2(COLOR_CREDITS, "Enter second Matrix Elements.\n");
                printf2(COLOR_CREDITS, "And when you reach desired columns number, %s.\n\n", isSett(BOOLS_MATRIXPARSING) ? "press ENTER":"enter an alphanumeric value");

                (*matrix2) = malloc(MAX_DIMENSIONS*sizeof(ityp *));
                errMem((*matrix2), (free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE), free(matrix1)));

                (*matrix2)[0] = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));

                for(dim[COLUMNS2] = 0; (tmp > 0 || tmp == -2) && !(assert); ++ dim[COLUMNS2])
                {
                    const dim_typ analog_columns = dim[COLUMNS2]+1;
                    if(lazy_exec)
                        (*matrix2)[0] = realloc((*matrix2)[0], sizeof(ityp)*analog_columns);
                    else if(!(analog_columns % STABILIZER_FACTOR))
                        (*matrix2)[0] = realloc((*matrix2)[0], sizeof(ityp)*(dim[COLUMNS2] + access(curLayout)->stabilizer_factor));

                    errMem((*matrix2)[0], (free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE), free(matrix1)));

                    if((tmp = insertElement((*matrix2), (dim_typ2){0, dim[COLUMNS2]}, false)) == -1 && tmp != -2 && getItemsListNo(MATRICES) != STARTING_MATNO)
                    {
                        if(access(curMatrix)->dim[RAWS] != dim[COLUMNS])
                        {
                            printErr(33, "You cannot use Current Matrix because\nit doesn't have %hu Raws right as First Matrix", dim[COLUMNS]);
                            dim[COLUMNS2] --;
                            tmp = 1;
                        }
                        else
                        {
                            equalMatrix(matrix2, access(curMatrix)->matrix, access(curMatrix)->dim);
                            printf2(COLOR_USER, "\nYou're correctly using Current Matrix.\n\n");
                            dim[COLUMNS2] = access(curMatrix)->dim[COLUMNS] -1;
                        }
                    }

                    if(!checkBackTracking(tmp, &dim[COLUMNS2]))
                    {
                    	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                        #pragma omp parallel for num_threads(algebra_units)
							for(ii=0; ii<algebra_units; ++ii)
	                            matrixFree(&matrix1[ii], dim[RAWS]);
	                    else
	                    	#pragma omp parallel for
							for(ii=0; ii<algebra_units; ++ii)
	                            matrixFree(&matrix1[ii], dim[RAWS]);
                        matrixFree(matrix2, dim[COLUMNS]);
                        free(matrix1);
                        free(matrix2);
                        #if WINOS
                            SetExitButtonState(ENABLED);
                        #endif // WINOS
                        return;
                    }
                }

                if(tmp != -1)
                {
                    dim[COLUMNS2] += 1 -(2*!(assert));

                    if(!(assert))
                    {
                        (*matrix2)[0] = realloc((*matrix2)[0], sizeof(ityp)*dim[COLUMNS2]); // TRUNKING MATRIX
                        errMem((*matrix2)[0], (free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE), free(matrix1)));
                    }

                    (*matrix2) = realloc((*matrix2), sizeof(ityp *)*(dim[COLUMNS]+!(lazy_exec)+!(assert)));
                    errMem((*matrix2), (free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE), free(matrix1)));

                    for(i = !(assert); tmp != -1 && i<dim[COLUMNS]; ++i)
                    {
                        (*matrix2)[i] = malloc(sizeof(ityp)*dim[COLUMNS2]);
                        errMem((*matrix2)[i], ((free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE), free(matrix1))););
                        for(j = start_col_index; tmp != -1 && j<dim[COLUMNS2]; ++j)
                        {
                            while((tmp = insertElement((*matrix2), (dim_typ2){i, j}, false)) != 1 && tmp != -2 && !(char_insert))
                                if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                                    if(access(curMatrix)->dim[RAWS] != dim[COLUMNS])
                                        printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Raws right as First Matrix", dim[COLUMNS]);
                                    else
                                    {
                                        equalMatrix(matrix2, access(curMatrix)->matrix, access(curMatrix)->dim);
                                        printf2(COLOR_USER, "\nYou're correctly using Current Matrix.\n\n");
                                        dim[COLUMNS2] = access(curMatrix)->dim[COLUMNS] -1;
                                        break;
                                    }
                                else
                                    printErr(1, "You cannot enter characters");

                            if((back_tracking = checkBackTracking2(tmp, &i, &j, &start_col_index, dim[COLUMNS2])) == 1)
                                break;
                            else if(back_tracking == 2)
                            {
                            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                                #pragma omp parallel for num_threads(algebra_units)
									for(ii=0; ii<algebra_units; ++ii)
	                                    matrixFree(&matrix1[ii], dim[RAWS]);
	                            else
	                            	 #pragma omp parallel for
									for(ii=0; ii<algebra_units; ++ii)
	                                    matrixFree(&matrix1[ii], dim[RAWS]);
                                matrixFree(matrix2, dim[COLUMNS]);
                                free(matrix1);
                                free(matrix2);
                                #if WINOS
                                    SetExitButtonState(ENABLED);
                                #endif // WINOS
                                return;
                            }
                        }

                    }
                }
            }
            else
            {
                if(!insertDim(&dim[COLUMNS2], COLUMNS))
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
	                else
	                	#pragma omp parallel for
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
                    free(matrix1);
                    free(matrix2);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif // WINOS
                    return;
                }

                CLEARBUFFER();

                if(!matrixAlloc(matrix2, (dim_typ2){dim[COLUMNS], dim[COLUMNS2]}))
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
	                else
	                	 #pragma omp parallel for
						for(ii=0; ii<algebra_units; ++ii)
	                        matrixFree(&matrix1[ii], dim[RAWS]);
                    free(matrix1);
                    free(matrix2);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif // WINOS
                    return;
                }

                for(i=0; tmp != -1 && i<dim[COLUMNS2] && !(assert); ++i)
                {
                    while((tmp = insertElement((*matrix2), (dim_typ2){0, i}, false)) != 1 && tmp != -2 && !(char_insert))
                        if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                            if(access(curMatrix)->dim[RAWS] != dim[COLUMNS])
                                printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Columns Right as First Matrix", dim[COLUMNS]);
                            else
                            {
                                equalMatrix(matrix2, access(curMatrix)->matrix, access(curMatrix)->dim);
                                printf2(COLOR_USER, "\nYou're correctly using Current Matrix.\n\n");
                                break;
                            }
                        else
                            printErr(1, "You cannot enter characters");

                    if(!checkBackTracking(tmp, &i))
                    {
                    	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                        #pragma omp parallel for num_threads(algebra_units)
							for(ii=0; ii<algebra_units; ++ii)
	                            matrixFree(&matrix1[ii], dim[RAWS]);
	                	else
	                		#pragma omp parallel for
							for(ii=0; ii<algebra_units; ++ii)
	                            matrixFree(&matrix1[ii], dim[RAWS]);
                        matrixFree(matrix2, dim[COLUMNS]);
                        free(matrix1);
                        free(matrix2);
                        #if WINOS
                            SetExitButtonState(ENABLED);
                        #endif
                        return;
                    }

                }

                for(i = !(assert); tmp != -1 && i<dim[COLUMNS]; ++i)
                    for(j=start_col_index; tmp != -1 && j<dim[COLUMNS2]; ++j)
                    {
                        while((tmp = insertElement((*matrix2), (dim_typ2){i, j}, false)) != 1 && tmp != -2 && !(char_insert))
                            if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                                if(access(curMatrix)->dim[RAWS] != dim[COLUMNS])
                                    printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Columns Right as First Matrix", dim[COLUMNS]);
                                else
                                {
                                    equalMatrix(matrix2, access(curMatrix)->matrix, access(curMatrix)->dim);
                                    printf2(COLOR_USER, "\nYou're correctly using Current Matrix.\n\n");
                                    break;
                                }
                            else
                                printErr(1, "You cannot enter characters");

                        if((back_tracking = checkBackTracking2(tmp, &i, &j, &start_col_index, dim[COLUMNS2])) == 1)
                            break;
                        else if(back_tracking == 2)
                        {
                        	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                            #pragma omp parallel for num_threads(algebra_units)
								for(ii=0; ii<algebra_units; ++ii)
	                                matrixFree(&matrix1[ii], dim[RAWS]);
	                        else
	                        	#pragma omp parallel for
								for(ii=0; ii<algebra_units; ++ii)
	                                matrixFree(&matrix1[ii], dim[RAWS]);
                            matrixFree(matrix2, dim[COLUMNS]);
                            free(matrix1);
                            free(matrix2);
                            #if WINOS
                                SetExitButtonState(ENABLED);
                            #endif // WINOS
                            return;
                        }

                    }

            }

        // Entering eventually IMAGINARY PART MATRIX of COMPLEX MATRIX2
        for(ii=1; ii<algebra_units; ++ii, ++argv_pos)
        {
            matrix2[ii] = NULL;
            printf2(COLOR_CREDITS, "Enter IMAGINARY %s PART of the Second Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printf2(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[COLUMNS], dim[COLUMNS2]);

            if(argc > argv_pos)
            {
                if((!matrixToken(argv[argv_pos], &matrix2[ii], cdim, &cdim[COLUMNS])) || cdim[RAWS] != dim[COLUMNS] || cdim[COLUMNS] != dim[COLUMNS2])
                {
                    #pragma omp parallel for
					for(ij=0; ij<=ii; ++ij)
                    {
                        matrixFree(&matrix1[ij], dim[RAWS]);
                        matrixFree(&matrix2[ij], cdim[RAWS]);
                    }
                    free(matrix1);
                    free(matrix2);
                    printUsage(&alg_operations[__pmode__]);
                    return;
                }
            }
            else
            {
                if(!insertNMMatrix(&matrix2[ii], (dim_typ2){dim[COLUMNS], dim[COLUMNS2]}))
                {
                    matrixFree(matrix1, dim[RAWS]);
                    matrixFree(matrix2, dim[COLUMNS]);
                    #pragma omp parallel for
					for(ij=1; ij<ii; ++ij)
                        matrixFree(&matrix1[ij], dim[RAWS]);
                    free(matrix1);
                    free(matrix2);
                    return;
                }
            }
        }

        printf2(COLOR_USER, "\n\nSecond %s Matrix correctly inserted:\n\n", complex_entries ? "Complex":NULL_CHAR);
        printMatrix(stdout, (*matrix2) , (dim_typ2){i, dim[COLUMNS2]});
        for(ii=1; ii<algebra_units; ++ii)
        {
            printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
            printMatrix(stdout, matrix2[ii], (dim_typ2){i, dim[COLUMNS2]});
        }

    }

    PRINT2N();

    if(assert_m)
        printf2(COLOR_SYSTEM, "First Matrix [%hu x %hu] raisen to the Power of: %llu is the [%hu x %hu] Matrix:", dim[RAWS], dim[RAWS], tmp2, dim[RAWS], dim[RAWS]);
    else
        printf2(COLOR_SYSTEM, "Matrix %s [%hu x %hu] between First and Second Matrix is:", INVERSE_OPS ? "Quotient":"Product", dim[RAWS], dim[COLUMNS2]);

    PRINT2N();

    if(__pmode__ == ALGEBRA_SCALARPRODUCT)
        printf2(COLOR_SYSTEM, "Whose unique element is the Scalar Product between two inserted Vectors is:");

    PRINTN();

    dim_typ x;
    ityp ***matrix_product = NULL;
    const sel_typ idx = assert_m ? RAWS : COLUMNS2;


    matrix_product = malloc(sizeof(ityp**)*algebra_units);
    errMem(matrix_product, (free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE),
                            free_foreach(matrix2, algebra_units, dim[idx], tmp2 == 1),
                            free(matrix1),
                            free(matrix2)));

    bool mxprodErr = false;

    if(tmp2 != 1)
		for(ii=0; ii<algebra_units; ++ii)
            if((mxprodErr = !matrixAlloc(&matrix_product[ii], (dim_typ2){dim[RAWS], dim[COLUMNS2]})))
                break;

    x = 1;

    if(tmp2 != 1)
    {
        if(mxprodErr && (!tmp2))
        {
        	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	            #pragma omp parallel for num_threads(algebra_units)
				for(ii=0; ii<algebra_units; ++ii)
	            {
	                matrixFree(&matrix_product[ii], dim[RAWS]);
	                matrixFree(&matrix1[ii], dim[RAWS]);
	            }
	        else
	            #pragma omp parallel for
				for(ii=0; ii<algebra_units; ++ii)
	            {
	                matrixFree(&matrix_product[ii], dim[RAWS]);
	                matrixFree(&matrix1[ii], dim[RAWS]);
	            }
            free(matrix1);
            free(matrix2);
            free(matrix_product);
            #if WINOS
                SetExitButtonState(ENABLED);
            #endif // WINOS
            return;
        }
        do
        {

            access(matrixProdFunc)(matrix1, matrix2, matrix_product, dim);

            if(assert_m)
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(ii=0; ii<algebra_units; ++ii)
	                    equalMatrix(&matrix2[ii], matrix_product[ii], (dim_typ2){dim[RAWS], dim[COLUMNS2]});
	            else
	            	#pragma omp parallel for
					for(ii=0; ii<algebra_units; ++ii)
	                    equalMatrix(&matrix2[ii], matrix_product[ii], (dim_typ2){dim[RAWS], dim[COLUMNS2]});
        }
        while(++x < tmp2);
    }
    else
        if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	        #pragma omp parallel for num_threads(algebra_units)
			for(ii=0; ii<algebra_units; ++ii)
	            matrix_product[ii] = matrix1[ii];
		else
			#pragma omp parallel for
			for(ii=0; ii<algebra_units; ++ii)
	            matrix_product[ii] = matrix1[ii];

    CHECK_INVERSE_OPERATIONS();

    PRINTL();
    printMatrix(stdout, (*matrix_product), (dim_typ2){dim[RAWS], dim[COLUMNS2]});
    for(ii=1; ii<algebra_units; ++ii)
    {
        printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
        printMatrix(stdout, matrix_product[ii], (dim_typ2){dim[RAWS], dim[COLUMNS2]});
    }


    // Overhead here isn't inevitable because the order is important
    if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	    #pragma omp parallel for num_threads(algebra_units)
		for(ii=0; ii<algebra_units; ++ii)
	        matrixFree(&matrix_product[ii], dim[RAWS]);
	else
		#pragma omp parallel for
		for(ii=0; ii<algebra_units; ++ii)
	        matrixFree(&matrix_product[ii], dim[RAWS]);

    if(tmp2 != 1)
        if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	        #pragma omp parallel for num_threads(algebra_units)
			for(ii=0; ii<algebra_units; ++ii)
	            matrixFree(&matrix2[ii], dim[idx]);
	    else
	    	#pragma omp parallel for
			for(ii=0; ii<algebra_units; ++ii)
	            matrixFree(&matrix2[ii], dim[idx]);

    if(!tmp2)
        if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	        #pragma omp parallel for num_threads(algebra_units)
			for(ii=0; ii<algebra_units; ++ii)
	            matrixFree(&matrix1[ii], dim[RAWS]);
	    else
	    	#pragma omp parallel for
			for(ii=0; ii<algebra_units; ++ii)
	            matrixFree(&matrix1[ii], dim[RAWS]);

    free(matrix1);
    free(matrix2);
    free(matrix_product);


    #if WINOS
        SetExitButtonState(ENABLED);
    #endif

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixKProduct(const register sel_typ argc, char ** argv)
{

    dim_typ dim[MAX_DIMENSIONS];
    dim_typ cdim[MAX_DIMENSIONS];
    dim_typ i, j;

    ityp ***matrix1 = NULL;
    const bool complex_entries = access(curLayout)->algebra != 0;
    const sel_typ algebra_units = exp2(access(curLayout)->algebra);

    matrix1 = malloc(sizeof(ityp**)*algebra_units);
    errMem(matrix1, VSPACE);

    if(argc)
    {
        if(!matrixToken(argv[0], matrix1, dim, &dim[COLUMNS]))
        {
            free(matrix1);
            return;
        }
    }
    else if(!enterMatrix(matrix1, dim, &dim[COLUMNS], false, !complex_entries))
    {
        free(matrix1);
        return;
    }

    sel_typ argv_pos = 0;

    // Entering eventually IMAGINARY PART MATRIX of COMPLEX MATRIX1
    for(i=1; i<algebra_units; ++i, ++argv_pos)
    {
        matrix1[i] = NULL;
        printf2(COLOR_CREDITS, "Enter IMAGINARY %s PART of the First Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
        printf2(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[RAWS], dim[COLUMNS]);
        if(argc > argv_pos)
        {
            if((!matrixToken(argv[argv_pos], &matrix1[i], cdim, &cdim[COLUMNS])) || cdim[RAWS] != dim[RAWS] || cdim[COLUMNS] != dim[COLUMNS])
            {
                #pragma omp parallel for
				for(j=0; j<=i; ++j)
                    matrixFree(&matrix1[j], cdim[RAWS]);
                free(matrix1);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {
            if(!insertNMMatrix(&matrix1[i], dim))
            {
                matrixFree(matrix1, dim[RAWS]);
                #pragma omp parallel for
				for(j=1; j<i; ++j)
                    matrixFree(&matrix1[j], dim[RAWS]);
                free(matrix1);
                return;
            }
        }
    }

    if(complex_entries)
    {
        printf2(COLOR_USER, "\n\nFirst Complex Matrix correctly inserted:\n\n");
        printMatrix(stdout, (*matrix1), dim);
        for(i=1; i<algebra_units; ++i)
        {
            printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printMatrix(stdout, matrix1[i], dim);
        }
    }

    ityp ***matrix2 = NULL;
    dim_typ dim2[MAX_DIMENSIONS];
    const bool assert_m = __pmode__ == ALGEBRA_MATRIXKPOWER;

    matrix2 = malloc(sizeof(ityp**)*algebra_units);
    errMem(matrix2, (free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE), free(matrix1)));

    ++ argv_pos;
    uint64_t tmp = 0;

    if(assert_m)
    {
        ityp tmp2 = 0.00;

        if(argc > argv_pos)
        {
            if(PARSING_MATRIX_ALLOWED)
            {
                if((!parse(argv[argv_pos], &tmp2)) || tmp2 != (tmp = (uint64_t)tmp2) || tmp < 1)
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(i=0; i<algebra_units; ++i)
	                        matrixFree(&matrix1[i], dim[RAWS]);
	                else
	                	#pragma omp parallel for
						for(i=0; i<algebra_units; ++i)
	                        matrixFree(&matrix1[i], dim[RAWS]);
                    free(matrix1);
                    free(matrix2);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif // WINOS
                    printUsage(&alg_operations[__pmode__]);
                    return;
                }
            }
            else if((tmp2 = strtod(argv[argv_pos], NULL)) != (tmp = (uint64_t)tmp2) || tmp < 1)
            {
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[i], dim[RAWS]);
	            else
	            	#pragma omp parallel
					for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[i], dim[RAWS]);
                free(matrix1);
                free(matrix2);
                #if WINOS
                    SetExitButtonState(ENABLED);
                #endif // WINOS
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {
            printf2(COLOR_CREDITS, "\nEnter N Power to which you want to raise the Inserted Matrix.\n");

            if(PARSING_SYSTEM_ALLOWED)
                PRINTHOWTOBACKMESSAGE();

            PRINT2N();

            while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp2 = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
                (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp2))) || tmp2 != (tmp = (uint64_t)tmp2) || tmp < 1)
            {
                CLEARBUFFER();
                if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
                if(isNullVal(tmp2) && exitHandleCheck)
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(i=0; i<algebra_units; ++i)
	                        matrixFree(&matrix1[i], dim[RAWS]);
	                else
	                	#pragma omp parallel for
						for(i=0; i<algebra_units; ++i)
	                        matrixFree(&matrix1[i], dim[RAWS]);
                    free(matrix1);
                    free(matrix2);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif // WINOS
                    return;
                }
                printErr(1, "Invalid inserted Value.\nExponent must be a non-negative and >= 1 integer");
            }

            CLEARBUFFER();

            if(tmp != 1)
            {
                if(!matrixAlloc(matrix2, dim))
                {
                	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                    #pragma omp parallel for num_threads(algebra_units)
						for(i=0; i<algebra_units; ++i)
	                        matrixFree(&matrix1[i], dim[RAWS]);
	                else
	                	#pragma omp parallel for
						for(i=0; i<algebra_units; ++i)
	                        matrixFree(&matrix1[i], dim[RAWS]);
                    free(matrix1);
                    free(matrix2);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif
                    return;
                }
                if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(i=0; i<algebra_units; ++i)
	                {
	                    (*matrix2)[i] = NULL;
	                    equalMatrix(&matrix2[i], matrix1[i], dim);
	                }
	            else
	            	#pragma omp parallel for
					for(i=0; i<algebra_units; ++i)
	                {
	                    (*matrix2)[i] = NULL;
	                    equalMatrix(&matrix2[i], matrix1[i], dim);
	                }
            }
            dim2[RAWS] = dim[RAWS];
            dim2[COLUMNS] = dim[COLUMNS];
        }
    }
    else
    {
        if(argc > argv_pos)
        {

            if((!matrixToken(argv[argv_pos], matrix2, dim2, &dim2[COLUMNS])))
            {
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[i], dim[RAWS]);
	            else
	            	#pragma omp parallel for
					for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[i], dim[RAWS]);
                matrixFree(matrix2, dim2[RAWS]);
                free(matrix1);
                free(matrix2);
                #if WINOS
                    SetExitButtonState(ENABLED);
                #endif // WINOS
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {

            if(!enterMatrix(matrix2, dim2, &dim2[COLUMNS], false, !complex_entries))
            {
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[i], dim[RAWS]);
	            else
	            	#pragma omp parallel for
					for(i=0; i<algebra_units; ++i)
	                    matrixFree(&matrix1[i], dim[RAWS]);
                free(matrix1);
                free(matrix2);
                return;
            }
        }

        for(i=1; i<algebra_units; ++i, ++argv_pos)
        {
            matrix2[i] = NULL;
            printf2(COLOR_CREDITS, "Enter IMAGINARY %s PART of the Second Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printf2(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim2[RAWS], dim2[COLUMNS]);
            if(argc > argv_pos)
            {
                if((!matrixToken(argv[argv_pos], &matrix2[i], cdim, &cdim[COLUMNS])) || cdim[RAWS] != dim2[RAWS] || cdim[COLUMNS] != dim2[COLUMNS])
                {
                    #pragma omp parallel for
					for(j=0; j<=i; ++j)
                    {
                        matrixFree(&matrix1[j], dim[RAWS]);
                        matrixFree(&matrix2[j], cdim[RAWS]);
                    }
                    free(matrix1);
                    free(matrix2);
                    printUsage(&alg_operations[__pmode__]);
                    return;
                }
            }
            else
            {
                if(!insertNMMatrix(&matrix2[i], dim2))
                {
                    #pragma omp parallel for
					for(j=0; j<i; ++j)
                    {
                        matrixFree(&matrix1[j], dim[RAWS]);
                        matrixFree(&matrix2[j], dim[RAWS]);
                    }
                    free(matrix1);
                    free(matrix2);
                    return;
                }
            }
        }

        if(complex_entries)
        {
            printf2(COLOR_USER, "\nSecond Complex Matrix correctly inserted:\n\n");
            printMatrix(stdout, (*matrix2), dim2);
            for(i=1; i<algebra_units; ++i)
            {
                printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
                printMatrix(stdout, matrix2[i], dim2);
            }
        }
    }


    const dim_typ dim3[MAX_DIMENSIONS] =
    {
        assert_m ? powi(dim[RAWS], tmp) : dim[RAWS]*dim2[RAWS],
        assert_m ? powi(dim[COLUMNS], tmp) : dim[COLUMNS]*dim2[COLUMNS]
    };

    PRINT2N();

    if(assert_m)
        printf2(COLOR_SYSTEM, "First Matrix [%hu x %hu] raisen to the Kronecker Power of: %llu is the [%hu x %hu] Matrix:", dim[RAWS], dim[COLUMNS], tmp, dim3[RAWS], dim3[COLUMNS]);
    else
        printf2(COLOR_SYSTEM, "Matrix Kronecker %s [%hu x %hu] between First and Second Matrix is:", INVERSE_OPS ? "Quotient":"Product", dim3[RAWS], dim3[COLUMNS]);

    PRINT2N();

    dim_typ x;
    ityp ***matrix_product = NULL;

    matrix_product = malloc(sizeof(ityp**)*algebra_units);
    errMem(matrix_product, (free_foreach(matrix1, algebra_units, dim[RAWS], NORMAL_MODE),
                            free_foreach(matrix2, algebra_units, dim2[RAWS], tmp == 1),
                            free(matrix1),
                            free(matrix2)));

    bool mxprodErr = false;

    if(tmp != 1)
		for(i=0; i<algebra_units; ++i)
            if((mxprodErr = !matrixAlloc(&matrix_product[i], dim3)))
                break;

    const dim_typ dims[MAX_DIMENSIONS][MAX_DIMENSIONS] =
    {
        {
            dim[RAWS],
            dim[COLUMNS]
        },
        {
            dim2[RAWS],
            dim2[COLUMNS]
        }
    };

    x = 1;


    if(tmp != 1)
    {
        if(mxprodErr && !tmp)
        {
        	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	            #pragma omp parallel for num_threads(algebra_units)
			    for(i=0; i<algebra_units; ++i)
	            {
	                matrixFree(&matrix_product[i], dim[RAWS]);
	                matrixFree(&matrix1[i], dim[RAWS]);
	            }
	        else
	        	#pragma omp parallel for
			    for(i=0; i<algebra_units; ++i)
	            {
	                matrixFree(&matrix_product[i], dim[RAWS]);
	                matrixFree(&matrix1[i], dim[RAWS]);
	            }
            free(matrix1);
            free(matrix2);
            free(matrix_product);
            #if WINOS
                SetExitButtonState(ENABLED);
            #endif // WINOS
            return;
        }
        do
        {

            access(matrixKProdFunc)(matrix1, matrix2, matrix_product, dims);

            if(assert_m)
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(i=0; i<algebra_units; ++i)
	                    equalMatrix(&matrix2[i], matrix_product[i], dim);
	            else
	            	#pragma omp parallel for
					for(i=0; i<algebra_units; ++i)
	                    equalMatrix(&matrix2[i], matrix_product[i], dim);
        }
        while(++x < tmp);
    }
    else
        if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	        #pragma omp parallel for num_threads(algebra_units)
			for(i=0; i<algebra_units; ++i)
	            matrix_product[i] = matrix1[i];
	    else
	    	#pragma omp parallel for
			for(i=0; i<algebra_units; ++i)
	            matrix_product[i] = matrix1[i];


    CHECK_INVERSE_OPERATIONS();
    PRINTL();
    printMatrix(stdout, (*matrix_product), dim3);

    /// UNTESTED!!!


    for(i=1; i<algebra_units; ++i)
    {
        printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
        printMatrix(stdout, matrix_product[i], dim3);
    }
    // Freeing Matrix (Real)

    // Overhead here isn't inevitable because the order is important

	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	    #pragma omp parallel for num_threads(algebra_units)
		for(i=0; i<algebra_units; ++i)
	        matrixFree(&matrix_product[i], dim3[RAWS]);
	else
		#pragma omp parallel for
		for(i=0; i<algebra_units; ++i)
	        matrixFree(&matrix_product[i], dim3[RAWS]);

    if(tmp != 1)
    	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	        #pragma omp parallel for num_threads(algebra_units)
			for(i=0; i<algebra_units; ++i)
	            matrixFree(&matrix2[i], dim2[RAWS]);
	    else
	    	#pragma omp parallel for
			for(i=0; i<algebra_units; ++i)
	            matrixFree(&matrix2[i], dim2[RAWS]);

    if(!tmp)
    	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	        #pragma omp parallel for num_threads(algebra_units)
			for(i=0; i<algebra_units; ++i)
	            matrixFree(&matrix1[i], dim[RAWS]);
	    else	
	    	#pragma omp parallel for
			for(i=0; i<algebra_units; ++i)
	            matrixFree(&matrix1[i], dim[RAWS]);

    free(matrix1);
    free(matrix2);
    free(matrix_product);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS


    return;
}


__MSSHELL_WRAPPER_ static void _MS__private perScalarProduct(const register sel_typ argc, char ** argv)
{
    ityp scal;
    dim_typ ii, ij;
    ityp ***matrix = NULL;

    const bool complex_entries = access(curLayout)->algebra != 0;
    const sel_typ algebra_units = exp2(access(curLayout)->algebra);
    dim_typ dim[MAX_DIMENSIONS];
    dim_typ cdim[MAX_DIMENSIONS];

    matrix = malloc(sizeof(ityp)*algebra_units);
    errMem(matrix, VSPACE);

    if(argc)
    {
        if(!matrixToken(argv[0], matrix, dim, &dim[COLUMNS]))
            return;
    }
    else if(!enterMatrix(matrix, dim, &dim[COLUMNS], false, !complex_entries))
        return;

    for(ii=1; ii<algebra_units; ++ii)
    {
        matrix[ii] = NULL;
        printf2(COLOR_CREDITS, "Enter IMAGINARY %s PART of the Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
        printf2(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[RAWS], dim[COLUMNS]);
        if(argc > 1)
        {
            if((!matrixToken(argv[1], &matrix[ii], cdim, &cdim[COLUMNS])) || cdim[RAWS] != dim[RAWS] || cdim[COLUMNS] != dim[COLUMNS])
            {  
				#pragma omp parallel for
				for(ij=0; ij<=ii; ++ij)
                    matrixFree(&matrix[ij], dim[RAWS]);
                printUsage(&alg_operations[ALGEBRA_PERSCALARPRODUCT]);
                return;
            }
        }
        else
        {
            if(!insertNMMatrix(&matrix[ii], dim))
            {
                #pragma omp parallel for
				for(ij=0; ij<ii; ++ij)
                    matrixFree(&matrix[ij], dim[RAWS]);
                return;
            }
        }
    }

    if(complex_entries)
    {
        printf2(COLOR_USER, "\nMatrix correctly inserted:\n\n");
        printMatrix(stdout, (*matrix), dim);
        for(ii=1; ii<algebra_units; ++ii)
        {
            printf2(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
            printMatrix(stdout, matrix[ii], dim);
        }
    }

    printf2(COLOR_CREDITS, "Enter a double floating-point Scalar Number.\n");

    if(PARSING_SYSTEM_ALLOWED)
        PRINTHOWTOBACKMESSAGE();

    PRINTN();

    const sel_typ argv_pos = 1+algebra_units;

    if(argc > argv_pos)
    {
        if(PARSING_SYSTEM_ALLOWED)
        {
            if(!parse(argv[argv_pos], &scal))
            {
            	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	                #pragma omp parallel for num_threads(algebra_units)
					for(ii=0; ii<algebra_units; ++ii)
	                    matrixFree(&matrix[ii], dim[RAWS]);
	            else
	            	#pragma omp parallel for
					for(ii=0; ii<algebra_units; ++ii)
	                    matrixFree(&matrix[ii], dim[RAWS]);
                return;
            }
        }
        else
            scal = strtod(argv[argv_pos], NULL);
    }
    else
        while(PARSING_SYSTEM_ALLOWED ? ((!(scal = requires(NULL, NULL_CHAR, "Inserted Scalar", PARSER_SHOWRESULT))) && dcheck && INVERSE_OPS) :
        scanf(INPUT_CONVERSION_FORMAT, &scal) != 1 || (dcheck && INVERSE_OPS && !scal))
        printErr(22, "Invalid inserted Value");

    if(isNullVal(scal))
    {
    	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
	        #pragma omp parallel for num_threads(algebra_units)
			for(ii=0; ii<algebra_units; ++ii)
	            matrixFree(&matrix[ii], dim[RAWS]);
	    else
	    	#pragma omp parallel for
			for(ii=0; ii<algebra_units; ++ii)
	            matrixFree(&matrix[ii], dim[RAWS]);
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    CLEARBUFFER();

    if(__pmode__ == ALGEBRA_SCALARDIVISIONMATRIX)
    {
        printf2(COLOR_SYSTEM, "Scalar ");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, scal);
        printf2(COLOR_SYSTEM, " has been correctly divided by the Matrix:\n\n");
    }
    else
    {
        printf2(COLOR_SYSTEM, "Matrix has been correctly %s by the Scalar ", INVERSE_OPS ? "divided" : "multiplied");
        printf2(COLOR_SYSTEM, OUTPUT_CONVERSION_FORMAT, scal);
        printf2(COLOR_SYSTEM, ":\n\n");
    }

    PRINTL();

    if(INVERSE_OPS && __pmode__ == ALGEBRA_PERSCALARPRODUCT)
        scal = (1/scal);

    dim_typ i, j, k;
    ityp (* const mul_func)(register ityp, register ityp) = __pmode__ == ALGEBRA_PERSCALARPRODUCT ? math_mul : math_div;

	if(algebra_units <= MIN_EXTENSIVE_MULTITHREADING_CORESNO || (isSett(BOOLS_EXTENSIVEMULTITHREADING) && algebra_units > MIN_EXTENSIVE_MULTITHREADING_CORESNO))
		#pragma omp parallel for num_threads(algebra_units)
	    for(k=0; k<algebra_units; ++k)
	    	#pragma omp parallel for
	        for(i=0; i<dim[RAWS]; ++i)
	        	#pragma omp parallel for
	            for(j=0; j<dim[COLUMNS]; ++j)
	                matrix[k][i][j] = mul_func(scal, matrix[k][i][j]);
	else
		#pragma omp parallel for
	    for(k=0; k<algebra_units; ++k)
	    	#pragma omp parallel for
	        for(i=0; i<dim[RAWS]; ++i)
	        	#pragma omp parallel for
	            for(j=0; j<dim[COLUMNS]; ++j)
	                matrix[k][i][j] = mul_func(scal, matrix[k][i][j]);

    CHECK_INVERSE_OPERATIONS();

    printMatrix(stdout, (*matrix), dim);
    matrixFree(matrix, dim[RAWS]);
    for(ii=1; ii<algebra_units; ++i)
    {
        printf("\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
        printMatrix(stdout, matrix[ii], dim);
        matrixFree(&matrix[ii], dim[RAWS]);
    }

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS


    return;
}

__MSSHELL_WRAPPER_ static void _MS__private illConditionChecking(const register sel_typ argc, char ** argv)
{
    ityp **matrix = NULL;
    dim_typ dim[MAX_DIMENSIONS];

    if(argc)
    {
        if(!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS]))
            return;

        if(dim[RAWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim[RAWS]);
            printUsage(&alg_operations[ALGEBRA_ILLCONDITIONCHECKING]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], true))
        return;

    const ityp norms1 = norms(matrix, dim[RAWS]);
    const dim_typ dimraws_per2 = dim[RAWS]<<1;

    matrix = realloc(matrix, sizeof(ityp *)*dimraws_per2);
    errMem(matrix, VSPACE);

    dim_typ i;

	for(i=0; i<dim[RAWS]; ++i)
    {
        matrix[i] = realloc(matrix[i], sizeof(ityp)*dimraws_per2);
        errMem(matrix[i], VSPACE);
    }

    if(!invertMatrix(matrix, dim[RAWS]))
    {
        printf("You cannot invert SINGULAR Matrices", 1);
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    printf2(COLOR_USER, "\nInserted Quad Matrix ILL CONDITION CHECKING is: ");
    printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, norms1*norms(matrix,dim[RAWS]));
    printf2(COLOR_USER, ".\n\n");

    matrixFree(&matrix, dimraws_per2);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private matrixFattLU(const register sel_typ argc, char ** argv)
{
    dim_typ dim[MAX_DIMENSIONS];
    ityp **matrix;
    ityp **L;
    ityp **U;

    matrix = L = U = NULL;

    if(argc)
    {
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[RAWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim[RAWS]);
            printUsage(&alg_operations[ALGEBRA_MATRIXLUFACTORIZATION]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[RAWS], dim[COLUMNS], true))
        return;

    bool assert[MAX_DIMENSIONS];

    if(((assert[LOWER_TRIANGULAR] = !matrixAlloc(&L, dim))) || ((assert[UPPER_TRIANGULAR] = !matrixAlloc(&U, dim))))
    {
        matrixFree(&matrix, dim[RAWS]);
        if(assert[LOWER_TRIANGULAR] && !assert[UPPER_TRIANGULAR])
            matrixFree(&L, dim[RAWS]);
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return;
    }

    if(!FattLU(dim[RAWS], matrix, L, U))
        printErr(1, "You cannot LU-decompose SINGULAR Matrices");
    else
    {
        printf2(COLOR_USER, "\nInserted Matrix has been correctly LU-decomposed\ninto Triangular Matrices respectively: L as Lower and U as Upper:\n\n");
        printMatrix(stdout, L, dim);
        printMatrix(stdout, U, dim);
    }

    matrixFree(&matrix, dim[RAWS]);
    matrixFree(&L, dim[RAWS]);
    matrixFree(&U, dim[RAWS]);

    #if WINOS
        SetExitButtonState(ENABLED);
    #endif // WINOS

    return;
}
