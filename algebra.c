// algebra.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"

#ifndef __DISABLE_MATMANAGER
__MSSHELL_WRAPPER_ static void  __apnt matrixManager(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __apnt matrixManager(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_MATMANAGER_PROGS, mat_manager, alg_operations[ALGOPS_MATRICESMANAGER].name, BY_CHARS);
    return;
}
#endif

#ifndef __DISABLE_SPECIALMATRICES
__MSSHELL_WRAPPER_ static void  __apnt specialMatrices(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __apnt specialMatrices(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(SPECMAT_LETSPECIALMATRIX, specmat_prog, alg_operations[ALGOPS_SPECIALMATRICES].name, BY_CHARS);
    return;
}
#endif

__MSSHELL_WRAPPER_ static void  matrixSort(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixEigenValues(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixNorm(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixDet(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixTrace(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixRank(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixSVD(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixInv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixCoFactor(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixTranspose(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixAdd(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixMultiplication(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixKProduct(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  perScalarMultiplication(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  illConditionChecking(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  matrixFattLU(const sel_typ argc, char ** argv);

sprog alg_operations[MAX_ALGEBRA_OPERATIONS] =
{
    #ifndef __DISABLE_MATMANAGER
    [ALGOPS_MATRICESMANAGER] =
    {
    	CMD_MATRIXMANAGER,
        NAME_MATRIXMANAGER,
        USAGE_MATRIXMANAGER,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXMANAGER,
        #endif
        matrixManager,
        ARGC_MATRIXMANAGER,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifndef __DISABLE_SPECIALMATRICES
    [ALGOPS_SPECIALMATRICES] =
    {
    	CMD_SPECIALMATRICES,
    	NAME_SPECIALMATRICES,
    	USAGE_SPECIALMATRICES,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_SPECIALMATRICES,
    	#endif
    	specialMatrices,
    	ARGC_SPECIALMATRICES,
    	AUTOMATIC,
    	FATHER   	
    },
    #endif
    [ALGOPS_MATRIXSORT] =
    {
    	CMD_MATRIXSORT,
        NAME_MATRIXSORT,
        USAGE_MATRIXSORT,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXSORT,
        #endif
        matrixSort,
        ARGC_MATRIXSORT,
        BY_USER,
        CHILD
    },
    [ALGOPS_MATRIXEIGVALUES] =
    {
    	CMD_MATRIXEIGVALUES,
    	NAME_MATRIXEIGVALUES,
    	USAGE_MATRIXEIGVALUES,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_MATRIXEIGVALUES,
    	#endif
    	matrixEigenValues,
    	ARGC_MATRIXEIGVALUES,
    	BY_USER,
    	CHILD
    },
    [ALGOPS_NORMCALCULATOR] =
    {
    	CMD_MATRIXNORM,
        NAME_MATRIXNORM,
        USAGE_MATRIXNORM,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXNORM,
        #endif
        matrixNorm,
        ARGC_MATRIXNORM,
        BY_USER,
        CHILD
    },
    [ALGOPS_DETERMINANTCALCULATOR] =
    {
    	CMD_MATRIXDET,
        NAME_MATRIXDET,
        USAGE_MATRIXDET,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXDET,
        #endif
        matrixDet,
        ARGC_MATRIXDET,
        BY_USER,
        CHILD
    },
    [ALGOPS_TRACECALCULATOR] =
    {
    	CMD_MATRIXTRACE,
        NAME_MATRIXTRACE,
        USAGE_MATRIXTRACE,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXTRACE,
        #endif
        matrixTrace,
        ARGC_MATRIXTRACE,
        BY_USER,
        CHILD
    },
    [ALGOPS_RANKCALCULATOR] =
    {
    	CMD_MATRIXRANK,
        NAME_MATRIXRANK,
        USAGE_MATRIXRANK,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXRANK,
        #endif
        matrixRank,
        ARGC_MATRIXRANK,
        BY_USER,
        CHILD
    },
    [ALGOPS_MATRIXSVD] =
    {
    	CMD_MATRIXSVD,
        NAME_MATRIXSVD,
        USAGE_MATRIXSVD,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXSVD,
        #endif
        matrixSVD,
        ARGC_MATRIXSVD,
        BY_USER,
        CHILD
    },
    [ALGOPS_INVERSEMATRIX] =
    {
    	CMD_MATRIXINV,
        NAME_MATRIXINV,
        USAGE_MATRIXINV,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXINV,
        #endif
        matrixInv,
        ARGC_MATRIXINV,
        BY_USER,
        CHILD
    },
    [ALGOPS_FASTINVERSEMATRIX] =
    {
    	CMD_FASTMATRIXINV,
    	NAME_FASTMATRIXINV,
    	USAGE_FASTMATRIXINV,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_FASTMATRIXINV,
    	#endif
    	matrixInv,
    	ARGC_FASTMATRIXINV,
    	BY_USER,
    	CHILD
    },
    [ALGOPS_MATRIXCOFACTOR] =
    {
    	CMD_MATRIXCOFACTOR,
    	NAME_MATRIXCOFACTOR,
    	USAGE_MATRIXCOFACTOR,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_MATRIXCOFACTOR,
    	#endif
    	matrixCoFactor,
    	ARGC_MATRIXCOFACTOR,
    	BY_USER,
    	CHILD
    },
    [ALGOPS_MATRIXADJOINT] =
    {
    	CMD_MATRIXADJOINT,
    	NAME_MATRIXADJOINT,
    	USAGE_MATRIXADJOINT,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_MATRIXADJOINT,
    	#endif
    	matrixCoFactor,
    	ARGC_MATRIXADJOINT,
    	BY_USER,
    	CHILD
    },
    [ALGOPS_MATRIXTRANSPOSE] =
    {
    	CMD_MATRIXTRANSPOSE,
        NAME_MATRIXTRANSPOSE,
        USAGE_MATRIXTRANSPOSE,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXTRANSPOSE,
        #endif
        matrixTranspose,
        ARGC_MATRIXTRANSPOSE,
        BY_USER,
        CHILD
    },
    [ALGOPS_MATRIXADD] =
    {
    	CMD_MATRIXADD,
        NAME_MATRIXADD,
        USAGE_MATRIXADD,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXADD,
        #endif
        matrixAdd,
        ARGC_MATRIXADD,
        BY_USER,
        CHILD
    },
    [ALGOPS_TENSORSSUM] =
    {
    	CMD_TENSORADD,
        NAME_TENSORADD,
        USAGE_TENSORADD,
        #ifndef __DISABLE_DATABASE
        LEVEL_TENSORADD,
        #endif
        matrixAdd,
        ARGC_TENSORADD,
        BY_USER,
        CHILD
    },
    [ALGOPS_MATRIXMULTIPLICATION] =
    {
    	CMD_MATRIXMULTIPLICATION,
        NAME_MATRIXMULTIPLICATION,
        USAGE_MATRIXMULTIPLICATION,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXMULTIPLICATION,
        #endif
        matrixMultiplication,
        ARGC_MATRIXMULTIPLICATION,
        BY_USER,
        CHILD
    },
    [ALGOPS_KRONPRODUCT] =
    {
    	CMD_MATRIXKPRODUCT,
        NAME_MATRIXKPRODUCT,
        USAGE_MATRIXKPRODUCT,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXKPRODUCT,
        #endif
        matrixKProduct,
        ARGC_MATRIXKPRODUCT,
        BY_USER,
        CHILD
    },
    [ALGOPS_MATRIXPOWER] =
    {
    	CMD_MATRIXPOWER,
        NAME_MATRIXPOWER,
        USAGE_MATRIXPOWER,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXPOWER,
        #endif
        matrixMultiplication,
        ARGC_MATRIXPOWER,
        BY_USER,
        CHILD
    },
    [ALGOPS_MATRIXKPOWER] =
    {
    	CMD_MATRIXKPOWER,
        NAME_MATRIXKPOWER,
        USAGE_MATRIXKPOWER,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXKPOWER,
        #endif
        matrixKProduct,
        ARGC_MATRIXKPOWER,
        BY_USER,
        CHILD
    },
    [ALGOPS_MATRIXPERVECTOR] =
    {
    	CMD_MATRIXPERVECTOR,
        NAME_MATRIXPERVECTOR,
        USAGE_MATRIXPERVECTOR,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXPERVECTOR,
        #endif
        matrixMultiplication,
        ARGC_MATRIXPERVECTOR,
        BY_USER,
        CHILD
    },
    [ALGOPS_DOTPRODUCT] =
    {
    	CMD_DOTPRODUCT,
        NAME_DOTPRODUCT,
        USAGE_DOTPRODUCT,
        #ifndef __DISABLE_DATABASE
        LEVEL_DOTPRODUCT,
        #endif
        matrixMultiplication,
        ARGC_DOTPRODUCT,
        BY_USER,
        CHILD
    },
    [ALGOPS_PERSCALARMULTIPLICATION] =
    {
    	CMD_PERSCALARMULTIPLICATION,
        NAME_PERSCALARMULTIPLICATION,
        USAGE_PERSCALARMULTIPLICATION,
        #ifndef __DISABLE_DATABASE
        LEVEL_PERSCALARMULTIPLICATION,
        #endif
        perScalarMultiplication,
        ARGC_PERSCALARMULTIPLICATION,
        BY_USER,
        CHILD
    },
    [ALGOPS_SCALARDIVISIONMATRIX] =
    {
    	CMD_SCALARDIVISIONMATRIX,
        NAME_SCALARDIVISIONMATRIX,
        USAGE_SCALARDIVISIONMATRIX,
        #ifndef __DISABLE_DATABASE
    	LEVEL_SCALARDIVISIONMATRIX,
    	#endif
        perScalarMultiplication,
        ARGC_SCALARDIVISIONMATRIX,
        BY_USER,
        CHILD
    },
    [ALGOPS_ILLCONDITIONCHECKING] =
    {
    	CMD_ILLCONDITIONCHECKING,
        NAME_ILLCONDITIONCHECKING,
        USAGE_ILLCONDITIONCHECKING,
        #ifndef __DISABLE_DATABASE
        LEVEL_ILLCONDITIONCHECKING,
        #endif
        illConditionChecking,
        ARGC_ILLCONDITIONCHECKING,
        BY_USER,
        CHILD
    },
    [ALGOPS_MATRIXLUFACTORIZATION] =
    {
    	CMD_MATRIXFATTLU,
        NAME_MATRIXFATTLU,
        USAGE_MATRIXFATTLU,
        #ifndef __DISABLE_DATABASE
        LEVEL_MATRIXFATTLU,
        #endif
        matrixFattLU,
        ARGC_MATRIXFATTLU,
        BY_USER,
        CHILD
    }
};

__MSSHELL_WRAPPER_ static void  matrixSort(const sel_typ argc, char ** argv)
{

    sel_typ tmp;
    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
        if(argc == 2)
        {

            if((tmp = strtod(argv[0], NULL)) == 3) return;
            if(tmp < 0 || tmp > 3)
            {
                printUsage(&alg_operations[ALGOPS_MATRIXSORT]);
                return;
            }

            if(!tmp)
            {
            	if(argv[1][strlen(argv[1])-1] != TERMINATING_CHAR)
            		strcat(argv[1], TERMINATING_STRING);
                if(!matrixToken(argv[1], &matrix, dim, &dim[COLUMNS]))
                {
                	matrixFree(&matrix, dim);
                	printUsage(&alg_operations[ALGOPS_MATRIXSORT]);
                    return;
            	}
            }
            else
            {
                if(!insertDims(&dim[ROWS], &dim[COLUMNS]))
                    return;

                if(!matrixAlloc(&matrix, dim))
                    return;

                if(!randomMatrix(matrix, dim))
                {
                    matrixFree(&matrix, dim);
                    return;
                }
            }

        }
        else
        {
            printUsage(&alg_operations[ALGOPS_MATRIXSORT]);
            return;
        }
    }
    else
    {

        for( ; ; )
        {

            CLEARBUFFER();

            msprintf(COLOR_CREDITS, "\nSelect Matrix Filling Mode.\n");
            msprintf(COLOR_CREDITS, "A for User Matrix Inserting, B for Random Matrix Generation;\n");
            msprintf(COLOR_CREDITS, "- %c to go Back..\n\n", access(curLayout)->exit_char);

            PRINTL();

            do if((tmp = toupper(getch())) == access(curLayout)->exit_char) return;
            while(tmp < 'A' && tmp > access(curLayout)->exit_char);

            if(tmp == 'A')
            {
                if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], false))
                    continue;
            }

            if(tmp == 'B')
            {
                if(!insertDims(&dim[ROWS], &dim[COLUMNS]))
                    continue;
                if(!matrixAlloc(&matrix, dim))
                    continue;
                if(!randomMatrix(matrix, dim))
                {
                    matrixFree(&matrix, dim);
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
    
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
    
    if(difftime)
    	gettimeofday(&tvBegin, NULL);
    qsort(matrix,dim[ROWS]*dim[COLUMNS], sizeof(mpfr_t), mpfr_cmpfunc);
    
    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    msprintf(COLOR_USER, "\nMatrix has been sort with Ascending Order:\n\n");
    PRINTL();
    printMatrix(stdout, &matrix, dim);
    matrixFree(&matrix, dim);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixEigenValues(const sel_typ argc, char ** argv)
{

    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
        	strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim);
            printUsage(&adv_calc[ALGOPS_MATRIXEIGVALUES]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], true))
        return;
        
    mpfr_t *eigenValues;
    dim_typ2 dim2 =
    {
    	1,
    	dim[ROWS]
    };
    
    if(!matrixAlloc(&eigenValues, dim2))
    {
    	matrixFree(&matrix, dim);
    	return;
    }
    
    mpfr_t *eigenVectors;
    
    if(!matrixAlloc(&eigenVectors, dim))
    {
    	matrixFree(&matrix, dim);
    	matrixFree(&eigenValues, dim2);
    	return;
    }
    
    sel_typ exit_state;
    
    struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
    
    if((exit_state = _matrixEigenValues(matrix, eigenValues, eigenVectors, dim[ROWS])) == EIGVALUES_INFEVS_ERROR)
    	printErr(33, "No convergence after %hu iterations! Probably Complex EigenValues", access(curLayout)->max_eigvalues_iterations);
    else if(exit_state == EIGVALUES_ALLOC_ERROR)
    	printErr(12, "Eigen Vectors Processing Heap Dynamic Memory Allocation Problem");
    else
    {
    	msprintf(COLOR_USER, "\nMatrix EigenValues are: ");
    	printMatrix(stdout, &eigenValues, dim2);
    	msprintf(COLOR_USER, "\nand its EigenVectors Matrix is: ");
    	printMatrix(stdout, &eigenVectors, dim);
    	
    	if(difftime)
		{
			PRINTL();
		    msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}
    }
    
    matrixFree(&matrix, dim);
    matrixFree(&eigenValues, dim2);
    matrixFree(&eigenVectors, dim);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixNorm(const sel_typ argc, char ** argv)
{
    sel_typ tmp;

    if(argc)
    {
        if((tmp = strtod(argv[0], NULL)) == 3) return;
        if(tmp < 0 || tmp > 3)
        {
            printUsage(&alg_operations[ALGOPS_NORMCALCULATOR]);
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

    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc > 1)
    {
    	if(argv[1][strlen(argv[1])-1] != TERMINATING_CHAR)
            strcat(argv[1], TERMINATING_STRING);
        if((!matrixToken(argv[1], &matrix, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim);
            printUsage(&alg_operations[ALGOPS_NORMCALCULATOR]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], true))
        return;

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	mpfr_t nrm;
	mpfr_init(nrm);
		
	if(tmp)
		norm(nrm, matrix, dim[ROWS], tmp-1);
	else
		norms(nrm, matrix, dim[ROWS]);
    
    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    matrixFree(&matrix, dim);

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        mpfr_set(*(access(exprVars)->e_ANS), nrm, MPFR_RNDN);

	char buf[MAX_BUFSIZ];
    mpfr_sprintf(buf, "\nInserted Square Matrix NORM is: %Rf.\n\n", nrm);
    msprintf(COLOR_USER, buf);
    mpfr_clear(nrm);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixDet(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
		if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim);
            printUsage(&alg_operations[ALGOPS_DETERMINANTCALCULATOR]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], true))
        return;


    bool flag = false;
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIME); 
    
    if(difftime)
    	gettimeofday(&tvBegin, NULL);
    	
	mpfr_t dt;
	
	mpfr_init(dt);
    det(dt, matrix, dim[ROWS], &flag);
    
    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    matrixFree(&matrix, dim);

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        mpfr_set(*(access(exprVars)->e_ANS), dt, MPFR_RNDN);
	
	char buf[MAX_BUFSIZ];
    mpfr_sprintf(buf, "\nInserted Square Matrix %sDETERMINANT is: %Rf.\n\n", flag ? "ABS-":NULL_CHAR, dt);
    msprintf(COLOR_USER, buf);
	mpfr_clear(dt);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixTrace(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix;
    dim_typ dim[2];


    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
        	strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim);
            printUsage(&alg_operations[ALGOPS_TRACECALCULATOR]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], true))
        return;

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
		
	mpfr_t trc;
	
	mpfr_init(trc);
	_matrixTrace(trc, matrix, dim[ROWS]);
    
    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    matrixFree(&matrix, dim);

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        mpfr_set(*(access(exprVars)->e_ANS), trc, MPFR_RNDN); 

	char buf[MAX_BUFSIZ];
    mpfr_sprintf(buf, "\nInserted Square Matrix TRACE is: %Rf.\n\n", trc);
    msprintf(COLOR_USER, buf);
	mpfr_clear(trc);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixRank(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
        	strcat(argv[0], TERMINATING_STRING);
        if(!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS]))
        {
        	matrixFree(&matrix, dim);
        	printUsage(&alg_operations[ALGOPS_RANKCALCULATOR]);
            return;
    	}
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], false))
        return;
        
    struct timeval tvBegin;
    const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
    
    if(difftime)
    	gettimeofday(&tvBegin, NULL);
    
	const register dim_typ rk = rank(matrix, dim);
    
    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    matrixFree(&matrix, dim);

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        mpfr_set_ui(*(access(exprVars)->e_ANS), rk, MPFR_RNDN);

    msprintf(COLOR_USER, "\nInserted Matrix RANK is: %hu.\n\n", rk);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixSVD(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if(!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS]))
    	{
    		matrixFree(&matrix, dim);
    		printUsage(&alg_operations[ALGOPS_MATRIXSVD]);
            return;
   		}
    }	
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], false))
        return;
        
    mpfr_t tmp, tmp2;
	mpfr_t *svd_vec;
	mpfr_t *Lm;
	mpfr_t *Rm;
	
	const register dim_typ secure_dim = MAX(dim[ROWS], dim[COLUMNS]);
    
    if(!matrixAlloc(&Lm, (dim_typ2){secure_dim, secure_dim}))
    {
    	matrixFree(&matrix, dim);
        return;
    }
    
    if(!matrixAlloc(&Rm, (dim_typ2){secure_dim, secure_dim}))
    {
    	matrixFree(&matrix, dim);
    	matrixFree(&Lm, ((dim_typ2){secure_dim, secure_dim}));
        return;
    }
    
    // const register dim_typ maxv = dim[dim[ROWS] >= dim[COLUMNS]];
    
    if(!matrixAlloc(&svd_vec, (dim_typ2){secure_dim,1}))
    {
    	matrixFree(&matrix, dim);
    	matrixFree(&Lm, ((dim_typ2){secure_dim, secure_dim}));
    	matrixFree(&Rm, ((dim_typ2){secure_dim, secure_dim}));
        return;
    }
    
	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	mpfr_init_set_d(tmp, SVD_EPSILON, MPFR_RNDN);
	mpfr_init_set_d(tmp2, SVD_TOLERANCE, MPFR_RNDN);
	
	svd(dim, SVD_WITHU, SVD_WITHV, tmp, tmp2, matrix, svd_vec, Lm, Rm);
    // dsvd(matrix, dim, S, V);
    
    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}
    
    matrixFree(&matrix, dim);
    mpfr_clears(tmp, tmp2, NULL);
    msprintf(COLOR_USER, "\nSINGULAR VALUES Vector is:\n");
    
    printMatrix(stdout, &svd_vec, ((dim_typ2){secure_dim,1})); 
	// printMatrix(stdout, S, (dim_typ2){dim[COLUMNS],1});
    
    printMatrix(stdout, &Lm, ((dim_typ2){dim[ROWS], dim[ROWS]}));
    printMatrix(stdout, &Rm, ((dim_typ2){dim[COLUMNS], dim[COLUMNS]}));
    
    SHOWPAUSEMESSAGE();
    // matrixFree(&S, ((dim_typ2){dim[COLUMNS],1}));
    // matrixFree(&V, ((dim_typ2){dim[COLUMNS],dim[COLUMNS]}));
    
    matrixFree(&svd_vec, ((dim_typ2){secure_dim,1}));
	matrixFree(&Lm, ((dim_typ2){secure_dim, secure_dim}));
	matrixFree(&Rm, ((dim_typ2){secure_dim, secure_dim}));
    
    PRINTL();
    PRINT2N();
    return;
}

__MSSHELL_WRAPPER_ static void  matrixInv(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim);
            printUsage(&alg_operations[ALGOPS_INVERSEMATRIX]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], true))
        return;

    mpfr_t *matrix2;

    if(!matrixAlloc(&matrix2, dim))
        return;

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	const register sel_typ exitVal = __pmode__ == ALGOPS_FASTINVERSEMATRIX ? invertMatrixFast(matrix, matrix2, dim[ROWS]) : invertMatrix(matrix, matrix2, dim[ROWS]);
		
    if(exitVal == INVERTMATRIX_SINGULAR)
        printErr(1, "You cannot invert SINGULAR Matrices");
    else if(exitVal == INVERTMATRIX_ALLOCERROR)
    	printErr(12, "Matrix Inverse Calculator Dynamic Memory Allocation Problem");
    else
    {
    	if(difftime)
	    {
	    	PRINTL();
	        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}
	
	    msprintf(COLOR_SYSTEM, "\nINVERSE MATRIX of inserted Quad Matrix is:\n\n");
	    printMatrix(stdout, __pmode__ == ALGOPS_FASTINVERSEMATRIX ? &matrix2 : &matrix, dim);
    }
    
    // printMatrix(stdout, matrix2, dim);
    matrixFree(&matrix, dim);
    matrixFree(&matrix2, dim);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixCoFactor(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim);
            printUsage(&alg_operations[ALGOPS_INVERSEMATRIX]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], true))
        return;
        
    mpfr_t *matrix2;

    if(!matrixAlloc(&matrix2, dim))
        return;

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	static bool (* const cof_funcs[2])(mpfr_t *restrict, mpfr_t *restrict, dim_typ) =
	{
		CoFactor,
		adjoint
	};
	
	const bool which_prog = __pmode__-ALGOPS_MATRIXCOFACTOR;
	bool (* const cof_func)(mpfr_t *restrict, mpfr_t *restrict, dim_typ) = cof_funcs[which_prog];
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	if(!cof_func(matrix, matrix2, dim[ROWS]))
	{
		printErr(12, "CoFactor Matrix Evaluating Process Heap Dynamic Memory Allocation Problem");
		matrixFree(&matrix, dim);
		matrixFree(&matrix2, dim);
		return;
	}

    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}
	
	msprintf(COLOR_SYSTEM, "\%s MATRIX of Inserted Matrix is:\n\n", which_prog?"ADJOINT":"COFACTOR");
    printMatrix(stdout, which_prog?&matrix:&matrix2, ((dim_typ2){dim[ROWS], dim[ROWS]}));

    matrixFree(&matrix, dim);
    matrixFree(&matrix2, dim);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixTranspose(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
        	strcat(argv[0], TERMINATING_STRING);
        if(!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS]))
        {
        	matrixFree(&matrix, dim);
        	printUsage(&alg_operations[ALGOPS_MATRIXTRANSPOSE]);
            return;
    	}
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], false))
        return;

    msprintf(COLOR_SYSTEM, "\nTRANSPOSED MATRIX of Inserted Matrix is:\n\n");

    mpfr_t *matrix2;
    
    if(!matrixAlloc(&matrix2, (dim_typ2){dim[COLUMNS], dim[ROWS]}))
    {
    	matrixFree(&matrix, dim);
        return;
    }

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
    transpose(matrix, matrix2, dim);
    
    if(difftime)
    {
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}
    
    printMatrix(stdout, &matrix2, ((dim_typ2){dim[COLUMNS], dim[ROWS]}));

    matrixFree(&matrix, dim);
    matrixFree(&matrix2, ((dim_typ2){dim[COLUMNS], dim[ROWS]}));
    return;
}

__MSSHELL_WRAPPER_ static void  matrixAdd(const sel_typ argc, char ** argv)
{
    dim_typ i, j, k, l;
    dim_typ dim[2];
    dim_typ cdim[2];
    uint64_t tdim=1;
    mpfr_t tmp;
    mpfr_t ***matrix1 = NULL;
    mpfr_t ***matrix2 = NULL;
    mpfr_t ***matrix_sum = NULL;
    const bool complex_entries = access(curLayout)->algebra != 0;
    const bool tensor_mode = __pmode__ == ALGOPS_TENSORSSUM;
    const sel_typ algebra_units = exp2(access(curLayout)->algebra);


    if(tensor_mode)
    {
        if(argc)
        {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (tdim = mpfr_get_ui(tmp, MPFR_RNDN))) || tdim < 1)
            {
            	mpfr_clear(tmp);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {

            msprintf(COLOR_CREDITS, "\nEnter Tensor Dimension.\n");
            PRINTHOWTOBACKMESSAGE();
            PRINT2N();
			
            while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (tdim = mpfr_get_ui(tmp, MPFR_RNDN))) || tdim < 1)
            {
				mpfr_clear(tmp);
                CLEARBUFFER();
                if(isNullVal(tmp) && exitHandleCheck)
					return;
                printErr(1, "Invalid inserted Value.\nTensor Dimension must be a non-negative and >= 1 integer");
            }

        }

    }

    char matrix_type[MINMIN_STRING];
    
    strcpy(matrix_type, INVERSE_OPS ? "Difference":"Sum");
    matrix1 = malloc(sizeof(mpfr_t**)*tdim);
    errMem(matrix1, VSPACE);
    mpfr_clear(tmp);

    for(i=0; i<tdim; ++i)
    {
        matrix1[i] = malloc(sizeof(mpfr_t*)*algebra_units);
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
            msprintf(COLOR_CREDITS, "Performing the %hu Matrices Component %s Operation.\n", k+1, matrix_type);

        if(argc > argv_pos)
        {
        	if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            	strcat(argv[argv_pos], TERMINATING_STRING);
            if(!matrixToken(argv[argv_pos], matrix1[k], dim, &dim[COLUMNS]))
            {
	       		#pragma omp parallel for num_threads(tdim)
	            for(i=0; i<tdim; ++i)
	                free(matrix1[i]);
                free(matrix1);
                return;
            }
        }
        else if(!enterMatrix(matrix1[k], dim, &dim[COLUMNS], false, !complex_entries))
        {
        	#pragma omp parallel for num_threads(tdim)
            for(i=0; i<tdim; ++i)
                free(matrix1[i]);
            free(matrix1);
            return;
        }

        // Entering eventually IMAGINARY PART MATRIX of COMPLEX MATRIX1
        for(i=1; i<algebra_units; ++i, ++argv_pos)
        {
            matrix1[k][i] = NULL;
            msprintf(COLOR_CREDITS, "Enter IMAGINARY %s PART of the First Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            msprintf(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[ROWS], dim[COLUMNS]);
            if(argc > argv_pos)
            {
            	if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            		strcat(argv[argv_pos], TERMINATING_STRING);
                if((!matrixToken(argv[argv_pos], &matrix1[k][i], cdim, &cdim[COLUMNS])) || cdim[ROWS] != dim[ROWS] || cdim[COLUMNS] != dim[COLUMNS])
                {
                        
	                #pragma omp parallel for num_threads(algebra_units)
	                for(j=0; j<algebra_units; ++j)
	                    matrixFree(&matrix1[k][j], cdim);
	                    
                	#pragma omp parallel for num_threads(tdim)
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
                    #pragma omp parallel for num_threads(algebra_units)
                    for(j=0; j<algebra_units; ++j)
                        matrixFree(&matrix1[k][j], dim);

                    #pragma omp parallel for num_threads(tdim)
					for(j=0; j<tdim; ++j)
                        free(matrix1[j]);

                    free(matrix1);
                    return;
                }
            }
        }

        if(complex_entries)
        {
            msprintf(COLOR_USER, "\n\nFirst Complex Matrix correctly inserted:\n\n");
            printMatrix(stdout, &matrix1[k][0], dim);
            for(i=1; i<algebra_units; ++i)
            {
                msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
                printMatrix(stdout, &matrix1[k][i], dim);
            }
        }

        if(!k)
        {
            matrix2 = malloc(sizeof(mpfr_t**)*tdim);
            if(checkErrMem(matrix2))
            {
            	#pragma omp parallel for num_threads(algebra_units)
                for(i=0; i<algebra_units; ++i)
                    matrixFree(&matrix1[k][j], dim);
	            
			    #pragma omp parallel for num_threads(tdim)
				for(i=0; i<tdim; ++i)
                    free(matrix1[i]);
                    
                free(matrix1);
                return;
            }

			for(i=0; i<tdim; ++i)
            {
                matrix2[i] = malloc(sizeof(mpfr_t*)*algebra_units);
                if(checkErrMem(matrix2[i]))
                {
                    #pragma omp parallel for num_threads(algebra_units)
					for(j=0; j<algebra_units; ++j)
                        matrixFree(&matrix1[k][j], dim);
	                
                    #pragma omp parallel for num_threads(tdim)
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
            dim_typ rc[2];
			if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            	strcat(argv[argv_pos], TERMINATING_STRING);
            if((!matrixToken(argv[argv_pos], matrix2[k], rc, &rc[COLUMNS])) || rc[ROWS] != dim[ROWS] || rc[COLUMNS] != dim[COLUMNS])
            {
            	#pragma omp parallel for num_threads(algebra_units)
                for(i=0; i<algebra_units; ++i)
                    matrixFree(&matrix1[k][i], dim);
	                    
                matrixFree(matrix2[k], rc);
                
                #pragma omp parallel for num_threads(tdim)
				for(i=0; i<tdim; ++i)
                {
                    free(matrix1[i]);
                    free(matrix2[i]);
                }
	                
                free(matrix1);
                free(matrix2);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {

            msprintf(COLOR_CREDITS, "Enter second [%hu x %hu] Matrix Elements.\n", dim[ROWS], dim[COLUMNS]);
            msprintf(COLOR_CREDITS, "Elements Inserting Carriage will be automatically redirected on second scansion (for ROWS Scanning).\n\n");

            if(!insertNMMatrix(matrix2[k], dim)) // dim[ROWS], dim[COLUMNS]))
            {
                #pragma omp parallel for num_threads(algebra_units)
			    for(i=0; i<algebra_units; ++i)
                    matrixFree(&matrix1[k][i], dim);
	            
                #pragma omp parallel for num_threads(tdim)
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
            msprintf(COLOR_CREDITS, "Enter IMAGINARY %s PART of the Second Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            msprintf(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[ROWS], dim[COLUMNS]);
            if(argc > argv_pos)
            {
				if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            		strcat(argv[argv_pos], TERMINATING_STRING);
                if((!matrixToken(argv[argv_pos], &matrix2[k][i], cdim, &cdim[COLUMNS])) || cdim[ROWS] != dim[ROWS] || cdim[COLUMNS] != dim[COLUMNS])
                {
                    #pragma omp parallel for
					for(j=0; j<=i; ++j)
                    {
                        matrixFree(&matrix1[k][j], dim);
                        matrixFree(&matrix2[k][j], dim); 
                    }
                    
                    #pragma omp parallel for num_threads(algebra_units)
					for(j=0; j<algebra_units; ++j)
                    {
                        matrixFree(&matrix1[k][j], dim);
                        matrixFree(&matrix2[k][j], dim);
                    }
	                    
                    #pragma omp parallel for num_threads(tdim)
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
                if(!insertNMMatrix(&matrix2[k][i], (dim_typ2){dim[ROWS], dim[COLUMNS]}))
                {
                	#pragma omp parallel for
                    for(j=0; j<i; ++j)
                        matrixFree(&matrix1[k][j], dim);
                        
                    #pragma omp parallel for num_threads(algebra_units)
					for(j=0; j<algebra_units; ++j)
                        matrixFree(&matrix1[k][j], dim);
	                        
                    matrixFree(matrix2[k], dim);
                    
                    #pragma omp parallel for num_threads(tdim)
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

        msprintf(COLOR_USER, "\nSecond Matrix correctly inserted:\n\n");
        printMatrix(stdout, &matrix2[k][0], dim); // dim[ROWS], dim[COLUMNS]);

        for(i=1; i<algebra_units; ++i)
        {
            msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printMatrix(stdout, &matrix2[k][i], dim);
        }

        msprintf(COLOR_SYSTEM, "\nMatrix %s [%hu x %hu] between First and Second Matrix is:\n", matrix_type, dim[ROWS], dim[COLUMNS]);

        if(!k)
        {
            matrix_sum = malloc(sizeof(mpfr_t**)*tdim);
            if(checkErrMem(matrix_sum))
            {
                #pragma omp parallel for num_threads(tdim)
				for(i=k; i<tdim; ++i)
                    #pragma omp parallel for num_threads(algebra_units)
					for(j=0; j<algebra_units; ++j)
                    {
                        matrixFree(&matrix1[i][j], dim);
                        matrixFree(&matrix2[i][j], dim);
                    }
                    
                #pragma omp parallel for num_threads(tdim)
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
                matrix_sum[i] = malloc(sizeof(mpfr_t*)*algebra_units);
                if(checkErrMem(matrix_sum[i]))
                {
                    #pragma omp parallel for num_threads(tdim)
					for(l=k; l<tdim; ++l)
                        #pragma omp parallel for num_threads(algebra_units)
						for(j=0; j<algebra_units; ++j)
                        {
                            matrixFree(&matrix1[l][j], dim);
                            matrixFree(&matrix2[l][j], dim);
                        }

                    #pragma omp parallel for num_threads(tdim)
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
                    matrixFree(&matrix1[i][j], dim);
                    matrixFree(&matrix2[i][j], dim);
                }
                
            #pragma omp parallel for num_threads(tdim)
			for(i=0; i<tdim; ++i)
            {
                free(matrix1[i]);
                free(matrix2[i]);
                free(matrix_sum[i]);
            }
            
            free(matrix1);
            free(matrix2);
            free(matrix_sum);
            return;
        }
        
        static void (* const matrixAddFuncs[_MAX_ALGEBRA][2])(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]) =
   		{
   			{
	    		_matrixAdd,
	    		_matrixSub
	    	},
	    	{
	    		_matrixCAdd,
	    		_matrixCSub
	    	},
	    	{
	    		_matrixQAdd,
	    		_matrixQSub 
	    	},
	    	{
	    		_matrixOAdd,
	    		_matrixOSub
	    	},
	    	{
	    		_matrixSAdd,
	    		_matrixSSub
	    	}
    	}; 
        
        struct timeval tvBegin;
        const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
        
        if(difftime)
        	gettimeofday(&tvBegin, NULL);
        	
        matrixAddFuncs[access(curLayout)->algebra][INVERSE_OPS](matrix1[k], matrix2[k], matrix_sum[k], dim, (dim_typ3){dim[COLUMNS],dim[COLUMNS],dim[COLUMNS]});
        
        if(difftime)
    	{
	    	PRINTL();
	        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}

        PRINTL();
        printMatrix(stdout, &matrix_sum[k][0], dim);
        for(i=1; i<algebra_units; ++i)
        {
            msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printMatrix(stdout, &matrix_sum[k][i], dim);
            matrixFree(&matrix1[k][i], dim);
            matrixFree(&matrix2[k][i], dim);
            matrixFree(&matrix_sum[k][i], dim);
        }
        // Freeing Matrix

        matrixFree(matrix1[k], dim);
        matrixFree(matrix2[k], dim);
        matrixFree(matrix_sum[k], dim);

        free(matrix1[k]);
        free(matrix2[k]);
        free(matrix_sum[k]);
    }

    free(matrix1);
    free(matrix2);
    free(matrix_sum);
    return;
}

static unsigned int nextPowerOf2(unsigned int n)
{
    unsigned int p = 1;
    if (n && !(n & (n - 1)))
        return n;
 
    while (p < n)
        p <<= 1;
    return p;
}


__MSSHELL_WRAPPER_ static void  matrixMultiplication(const sel_typ argc, char ** argv)
{
    dim_typ dim[3];
    dim_typ cdim[2];
    dim_typ ii, ij;
    char tmp;

    mpfr_t **matrix1;
    const bool assert_m = __pmode__ == ALGOPS_MATRIXPOWER;
    const bool complex_entries = access(curLayout)->algebra != 0;
    const sel_typ algebra_units = exp2(access(curLayout)->algebra);

    matrix1 = malloc(sizeof(mpfr_t*)*algebra_units);
    errMem(matrix1, VSPACE);

    // Richiesta di inserimento valori Prima Matrice.

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], matrix1, dim, &dim[COLUMNS])) || (assert_m && dim[ROWS] != dim[COLUMNS]))
        {
            matrixFree(matrix1, dim);
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
        msprintf(COLOR_CREDITS, "Enter IMAGINARY %s PART of the First Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
        msprintf(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[ROWS], dim[COLUMNS]);
        if(argc > argv_pos)
        {
        	if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            	strcat(argv[argv_pos], TERMINATING_STRING);
            if((!matrixToken(argv[argv_pos], &matrix1[ii], cdim, &cdim[COLUMNS])) || cdim[ROWS] != dim[ROWS] || cdim[COLUMNS] != dim[COLUMNS])
            {
                #pragma omp parallel for
				for(ij=0; ij<=ii; ++ij)
                    matrixFree(&matrix1[ij], dim);
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
                    matrixFree(&matrix1[ij], dim);
                free(matrix1);
                return;
            }
        }
    }

    if(complex_entries)
    {
        msprintf(COLOR_USER, "\n\nFirst Complex Matrix correctly inserted:\n\n");
        printMatrix(stdout, matrix1, dim);
        for(ii=1; ii<algebra_units; ++ii)
        {
            msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
            printMatrix(stdout, &matrix1[ii], dim);
        }
    }
    
    // Richiesta di inserimento valori Seconda Matrice.

    tmp = 1;

    uint64_t tmp2 = 0;
    mpfr_t **matrix2 = malloc(sizeof(mpfr_t*)*algebra_units);
    errMem(matrix2, (free_foreach(matrix1, algebra_units, (dim_typ2){dim[ROWS], dim[COLUMNS]}, NORMAL_MODE), free(matrix1)));

    ++ argv_pos;

    if(assert_m)
    {
    	mpfr_t tmp3;
        if(argc > argv_pos)
        {
        	if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            	strcat(argv[argv_pos], TERMINATING_STRING);
            if((!parse(argv[argv_pos], &tmp3)) || mpfr_cmp_ui(tmp3, (tmp2 = mpfr_get_ui(tmp3, MPFR_RNDN))) || tmp2 < 1)
            {
                #pragma omp parallel for num_threads(algebra_units)
				for(ii=0; ii<algebra_units; ++ii)
                    matrixFree(&matrix1[ii], dim);
                free(matrix1);
                free(matrix2);
                mpfr_clear(tmp3);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {
            msprintf(COLOR_CREDITS, "\nEnter N Power to which you want to raise the Inserted Square Matrix.\n");
            PRINTHOWTOBACKMESSAGE();
            PRINT2N();

            while(requires(tmp3, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp3) || mpfr_cmp_ui(tmp3, (tmp2 = mpfr_get_ui(tmp3, MPFR_RNDN))) || tmp2 < 1)
            {
                CLEARBUFFER();
                if(isNullVal(tmp3) && exitHandleCheck)
                {
                    #pragma omp parallel for num_threads(algebra_units)
					for(ii=0; ii<algebra_units; ++ii)
                        matrixFree(&matrix1[ii], dim);
                    mpfr_clear(tmp3);
                    free(matrix1);
                    free(matrix2);
                    return;
                }
            	mpfr_clear(tmp3);
                printErr(1, "Invalid inserted Value.\nExponent must be a non-negative and >= 1 integer");
            }

            CLEARBUFFER();

            if(tmp2 != 1)
            {
				for(ii=0; ii<algebra_units; ++ii)
                    if(!equalMatrix(&matrix2[ii], matrix1[ii], dim, EQUALMATRIX_REALLOC))
                    {
                    	#pragma omp parallel for num_threads(algebra_units)
                    	for(ij=0; ij<algebra_units; ++ij)
                        	matrixFree(&matrix1[ij], dim);
                    	#pragma omp parallel for
						for(ij=0; ij<ii; ++ij)
	                        matrixFree(&matrix2[ij], ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
	                    mpfr_clear(tmp3); 
	                    free(matrix1);
	                    free(matrix2);
	                    return;
                    }
	            	
            }
            dim[COLUMNS2] = dim[ROWS];
        }
    }
    else
    {

        dim_typ i, j;
        dim_typ start_col_index=0;
        volatile sel_typ back_tracking;
        const bool assert = (__pmode__ == ALGOPS_DOTPRODUCT || __pmode__ == ALGOPS_MATRIXPERVECTOR);
        if(argc > argv_pos)
        {
            if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            	strcat(argv[argv_pos], TERMINATING_STRING);
            if((!matrixToken(argv[argv_pos], matrix2, &i, &dim[COLUMNS2])) || i != dim[COLUMNS])
            {
                matrixFree(matrix2, ((dim_typ2){i, dim[COLUMNS2]}));
                #pragma omp parallel for num_threads(algebra_units)
				for(ii=0; ii<algebra_units; ++ii)
                    matrixFree(&matrix1[ii], dim);
                free(matrix1);
                free(matrix2);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }// Overhead here isn't inevitable because the order is important
        else
        {
            msprintf(COLOR_CREDITS, "Enter second Matrix Elements.\n");
            msprintf(COLOR_CREDITS, "And when you reach desired columns number, press ENTER.\n\n");
			(*matrix2) = malloc(sizeof(mpfr_t)<<1);

            for(dim[COLUMNS2] = 0; (tmp > 0 || tmp == -2) && !(assert); ++ dim[COLUMNS2])
            {
                const dim_typ analog_columns = dim[COLUMNS2]+1;
                
          (*matrix2) = realloc((*matrix2), (sizeof(mpfr_t)<<1)*analog_columns);
                errMem((*matrix2), (free_foreach(matrix1, algebra_units, (dim_typ2){dim[ROWS], dim[COLUMNS]}, NORMAL_MODE), free(matrix1)));

                if((tmp = insertElement((*matrix2), (dim_typ2){0, dim[COLUMNS2]}, analog_columns, false, MPFR_INIT_ELEMENTS)) == -1 && tmp != -2 && getItemsListNo(MATRICES) != STARTING_MATNO)
                {
                    if(access(curMatrix)->dim[ROWS] != dim[COLUMNS])
                    {
                        printErr(33, "You cannot use Current Matrix because\nit doesn't have %hu Rows right as First Matrix", dim[COLUMNS]);
                        -- dim[COLUMNS2];
                        tmp = 1;
                    }
                    else
                    {
                    	matrixFree(matrix2, ((dim_typ2){1,analog_columns}));
                        if(!equalMatrix(matrix2, access(curMatrix)->matrix, access(curMatrix)->dim, EQUALMATRIX_REALLOC))
                        {
                        	#pragma omp parallel for num_threads(algebra_units)
							for(ii=0; ii<algebra_units; ++ii)
			                    matrixFree(&matrix1[ii], dim);
			                free(matrix1);
			                free(matrix2);
                        	return;
                        }
                        msprintf(COLOR_USER, "\nYou're correctly using Current Matrix.\n\n");
                        dim[COLUMNS2] = access(curMatrix)->dim[COLUMNS] -1;
                    }
                }

                if(!checkBackTracking(tmp, &dim[COLUMNS2]))
                {
                    #pragma omp parallel for num_threads(algebra_units)
					for(ii=0; ii<algebra_units; ++ii)
                        matrixFree(&matrix1[ii], dim);
                    matrixFree(matrix2, ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
                    free(matrix1);
                    free(matrix2);
                    return;
                }
            }

            if(tmp != -1)
            {
                dim[COLUMNS2] += 1 -((!(assert))<<1);

                for(i = !(assert); tmp != -1 && i<dim[COLUMNS]; ++i)
                {
                	(*matrix2) = realloc((*matrix2),(sizeof(mpfr_t)<<1)*((!assert)?dim[COLUMNS2]:1)*(dim[COLUMNS]+!(assert)));
                    errMem((*matrix2), ((free_foreach(matrix1, algebra_units, (dim_typ2){dim[ROWS], dim[COLUMNS]}, NORMAL_MODE), free(matrix1))););
                    for(j = start_col_index; tmp != -1 && j<dim[COLUMNS2]; ++j)
                    {
                        while((tmp = insertElement((*matrix2), (dim_typ2){i, j}, dim[COLUMNS2], false, MPFR_INIT_ELEMENTS)) != 1 && tmp != -2)
                            if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                                if(access(curMatrix)->dim[ROWS] != dim[COLUMNS])
                                    printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Rows right as First Matrix", dim[COLUMNS]);
                                else
                                {
                                	matrixFree(matrix2, ((dim_typ2){1,dim[COLUMNS2]}));
                                    if(!equalMatrix(matrix2, access(curMatrix)->matrix, access(curMatrix)->dim, EQUALMATRIX_REALLOC))
                                    {
                                    	#pragma omp parallel for num_threads(algebra_units)
                                    	for(ii=0; ii<algebra_units; ++ii)
			                                matrixFree(&matrix1[ii], dim);
			                            free(matrix1);
			                            free(matrix2);
                                    }
                                    msprintf(COLOR_USER, "\nYou're correctly using Current Matrix.\n\n");
                                    dim[COLUMNS2] = access(curMatrix)->dim[COLUMNS] -1;
                                    break;
                                }
                            else
                                printErr(1, "You cannot enter characters");

                        if((back_tracking = checkBackTracking2(tmp, &i, &j, &start_col_index, dim[COLUMNS2])) == 1)
                            break;
                        else if(back_tracking == 2)
                        {
                            #pragma omp parallel for num_threads(algebra_units)
							for(ii=0; ii<algebra_units; ++ii)
                                matrixFree(&matrix1[ii], dim);
                            matrixFree(matrix2, ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
                            free(matrix1);
                            free(matrix2);
                            return;
                        }
                    }

                }
            }
        }

        // Entering eventually IMAGINARY PART MATRIX of COMPLEX MATRIX2
        for(ii=1; ii<algebra_units; ++ii, ++argv_pos)
        {
            matrix2[ii] = NULL;
            msprintf(COLOR_CREDITS, "Enter IMAGINARY %s PART of the Second Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
            msprintf(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[COLUMNS], dim[COLUMNS2]);

            if(argc > argv_pos)
            {
            	if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            		strcat(argv[argv_pos], TERMINATING_STRING);
                if((!matrixToken(argv[argv_pos], &matrix2[ii], cdim, &cdim[COLUMNS])) || cdim[ROWS] != dim[COLUMNS] || cdim[COLUMNS] != dim[COLUMNS2])
                {
                	#pragma omp parallel for num_threads(algebra_units)
					for(ij=0; ij<algebra_units; ++ij)
                        matrixFree(&matrix1[ij], dim);
                    #pragma omp parallel for
					for(ij=0; ij<=ii; ++ij)
                        matrixFree(&matrix2[ij], ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
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
                    matrixFree(matrix1, dim);
                    matrixFree(matrix2, ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
                    #pragma omp parallel for
					for(ij=1; ij<ii; ++ij)
                        matrixFree(&matrix1[ij], dim);
                    free(matrix1);
                    free(matrix2);
                    return;
                }
                
                msprintf(COLOR_USER, "\n\nSecond %s Matrix correctly inserted:\n\n", complex_entries ? "Complex":NULL_CHAR);
		        printMatrix(stdout, matrix2, ((dim_typ2){i, dim[COLUMNS2]}));
		        for(ii=1; ii<algebra_units; ++ii)
		        {
		            msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
		            printMatrix(stdout, &matrix2[ii], ((dim_typ2){i, dim[COLUMNS2]}));
        		}
            }
    	}

    }

    PRINT2N();

    if(assert_m)
        msprintf(COLOR_SYSTEM, "First Matrix [%hu x %hu] raised to the Power of: %llu is the [%hu x %hu] Matrix:", dim[ROWS], dim[ROWS], tmp2, dim[ROWS], dim[ROWS]);
    else
        msprintf(COLOR_SYSTEM, "Matrix Multiplication [%hu x %hu] between First and Second Matrix is:", dim[ROWS], dim[COLUMNS2]);

    PRINT2N();

    if(__pmode__ == ALGOPS_DOTPRODUCT)
        msprintf(COLOR_SYSTEM, "Whose unique element is the Dot Product between two inserted Vectors is:");

    PRINTN();

    dim_typ x = 1;
    mpfr_t **matrix_product = malloc(sizeof(mpfr_t*)*algebra_units);
    errMem(matrix_product, (free_foreach(matrix1, algebra_units, (dim_typ2){dim[ROWS], dim[COLUMNS]}, NORMAL_MODE),
                            free_foreach(matrix2, algebra_units, (dim_typ2){dim[COLUMNS], dim[COLUMNS2]}, tmp2 == 1),
                            free(matrix1),
                            free(matrix2)));
    const sel_typ idx = assert_m ? ROWS : COLUMNS2;
    bool hasPadding = false;
    bool mxprodErr = false;
    
    static void (* const matrixMulFuncs[_MAX_ALGEBRA])(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]) =
    {
    	_matrixMultiplication,
        _matrixCMultiplication,
        _matrixQMultiplication,
        _matrixOMultiplication,
        _matrixSMultiplication
	};
	
	static void (* const matrixAddFuncs[_MAX_ALGEBRA])(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]) =
    {
    	_matrixAdd,
        _matrixCAdd,
        _matrixQAdd,
        _matrixOAdd,
        _matrixSAdd
	};
	
	static void (* const matrixSubFuncs[_MAX_ALGEBRA])(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]) =
    {
    	_matrixSub,
        _matrixCSub,
        _matrixQSub,
        _matrixOSub,
        _matrixSSub
	};
	
	dim_typ2 olddim =
	{
		dim[ROWS],
		dim[COLUMNS2]
	};
	
	mpfr_t logan, logan2, logan3;
	mpfr_init_set_ui(logan, dim[ROWS], MPFR_RNDN);
	mpfr_init_set_ui(logan2, dim[COLUMNS], MPFR_RNDN);
	mpfr_init_set_ui(logan3, dim[COLUMNS2], MPFR_RNDN);

	mpfr_log2(logan, logan, MPFR_RNDN);
	mpfr_log2(logan2, logan2, MPFR_RNDN);
	mpfr_log2(logan3, logan3, MPFR_RNDN);
		
	if(dim[ROWS] > access(curLayout)->min_strassen_dim && dim[COLUMNS] > access(curLayout)->min_strassen_dim && dim[COLUMNS2] > access(curLayout)->min_strassen_dim && (mpfr_cmp_ui(logan, mpfr_get_ui(logan, MPFR_RNDN)) || mpfr_cmp_ui(logan2, mpfr_get_ui(logan2, MPFR_RNDN)) || mpfr_cmp_ui(logan3, mpfr_get_ui(logan3, MPFR_RNDN))))
	{
		hasPadding = true;
		const register dim_typ np2 = nextPowerOf2(MAX(dim[COLUMNS2], MAX(nextPowerOf2(dim[ROWS]), nextPowerOf2(dim[COLUMNS]))));
		mpfr_t **tmpmat = malloc(sizeof(mpfr_t*)*algebra_units);
		mpfr_t **tmpmat2 = malloc(sizeof(mpfr_t*)*algebra_units);
		dim_typ i, j;
		#pragma omp parallel for num_threads(algebra_units) private(i)
		for(ii=0; ii<algebra_units; ++ii)
		{	
			// #pragma omp parallel sections private(i)
			{
				// #pragma omp section
				{
					/*
					if(!matrixAlloc(&tmpmat[ii], (dim_typ2){np2,np2}))
					{
						// Overhead here isn't inevitable because the order is important
					    #pragma omp parallel for num_threads(algebra_units)
						for(ii=0; ii<algebra_units; ++ii)
					        matrixFree(&matrix_product[ii], ((dim_typ2){dim[ROWS], dim[COLUMNS2]}));
					
					    if(tmp2 != 1)
					        #pragma omp parallel for num_threads(algebra_units)
							for(ii=0; ii<algebra_units; ++ii)
					            matrixFree(&matrix2[ii], ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
					
					    if(!tmp2)
					        #pragma omp parallel for num_threads(algebra_units)
							for(ii=0; ii<algebra_units; ++ii)
					            matrixFree(&matrix1[ii], dim);
					    return;
					}
					*/
					
					(void) matrixAlloc(&tmpmat[ii], (dim_typ2){np2,np2});
					for(i=0; i<dim[ROWS]; ++i)
						for(j=0; j<dim[COLUMNS]; ++j)
							mpfr_set(*(tmpmat[ii] + np2*i + j), *(matrix1[ii] + dim[COLUMNS]*i + j), MPFR_RNDN);
					matrixFree(&matrix1[ii], dim);
					
				}
				// #pragma omp section
				{
					/*
					if(!matrixAlloc(&tmpmat2[ii], (dim_typ2){np2,np2}))
					{
						// Overhead here isn't inevitable because the order is important
					    #pragma omp parallel for num_threads(algebra_units)
						for(ii=0; ii<algebra_units; ++ii)
					        matrixFree(&matrix_product[ii], ((dim_typ2){dim[ROWS], dim[COLUMNS2]}));
					
					    if(tmp2 != 1)
					        #pragma omp parallel for num_threads(algebra_units)
							for(ii=0; ii<algebra_units; ++ii)
					            matrixFree(&matrix2[ii], ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
					
					    if(!tmp2)
					        #pragma omp parallel for num_threads(algebra_units)
							for(ii=0; ii<algebra_units; ++ii)
					            matrixFree(&matrix1[ii], dim);
					    return;
					}
					*/
					(void) matrixAlloc(&tmpmat2[ii], (dim_typ2){np2,np2});
					for(i=0; i<dim[COLUMNS]; ++i)
						for(j=0; j<dim[COLUMNS2]; ++j)
							mpfr_set(*(tmpmat2[ii] + np2*i + j), *(matrix2[ii] + dim[COLUMNS2]*i + j), MPFR_RNDN);
					matrixFree(&matrix2[ii], ((dim_typ2){dim[COLUMNS],dim[COLUMNS2]}));

				}
			}
		}
		free(matrix1);
		free(matrix2);
		matrix1 = tmpmat;
		matrix2 = tmpmat2;
		dim[ROWS] = dim[COLUMNS] = dim[COLUMNS2] = np2; 
	}
	
	mpfr_clears(logan, logan2, logan3, NULL);
	
	if(tmp2 != 1)
	{
		for(ii=0; ii<algebra_units; ++ii)
            if((mxprodErr = !matrixAlloc(&matrix_product[ii], (dim_typ2){dim[ROWS], dim[COLUMNS2]})))
                break;
                
        if(mxprodErr && (!tmp2))
        {
        	#pragma omp parallel for num_threads(algebra_units)
			for(ij=0; ij<ii; ++ij)
	            matrixFree(&matrix2[ij], ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
            #pragma omp parallel for num_threads(algebra_units)
			for(ii=0; ii<algebra_units; ++ii)
            {
                matrixFree(&matrix_product[ii], ((dim_typ2){dim[ROWS], dim[COLUMNS2]}));
                matrixFree(&matrix1[ii], dim);
            }
            free(matrix1);
            free(matrix2);
            free(matrix_product);
            return;
        }
        struct timeval tvBegin;
        const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
        
        dim_typ i;
        
        if(difftime)
        	gettimeofday(&tvBegin, NULL);
        do
        {
        	
			mmult_fast(dim[ROWS], (dim_typ3){dim[COLUMNS],dim[COLUMNS2],dim[COLUMNS2]},algebra_units, matrix1, matrix2, matrix_product, matrixMulFuncs[access(curLayout)->algebra], matrixAddFuncs[access(curLayout)->algebra], matrixSubFuncs[access(curLayout)->algebra]);
        	/// matrixMulFuncs[access(curLayout)->algebra](matrix1, matrix2, matrix_product, dim);
			
			if(assert_m)
                //#pragma omp parallel for num_threads(algebra_units)
				for(ii=0; ii<algebra_units; ++ii)
				{
                 (void) equalMatrix(&matrix2[ii], matrix_product[ii], (dim_typ2){dim[ROWS], dim[COLUMNS2]}, EQUALMATRIX_NOREALLOC);
                	if(x != tmp2-1)
                    	for(i=0; i<dim[ROWS]*dim[COLUMNS]; ++i)	
                    		mpfr_set_ui(matrix_product[ii][i], 0, MPFR_RNDN); 
            	}
        }
        while(++x < tmp2);
        
        if(difftime)
    	{
	    	PRINTL();
	        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}
		
    }
    else
        #pragma omp parallel for num_threads(algebra_units)
		for(ii=0; ii<algebra_units; ++ii)
            matrix_product[ii] = matrix1[ii];

    PRINTL();
    
    dim_typ matrix_product_dim[2];
	
	if(hasPadding)
	{
		dim_typ i, j; 
		// #pragma omp parallel for num_threads(algebra_units)
		for(ii=0; ii<algebra_units; ++ii)
		{
			for(i=0; i<dim[ROWS]; ++i)
				for(j=0; j<dim[COLUMNS]; ++j)
					mpfr_set(*(matrix_product[ii] + olddim[COLUMNS]*i + j), *(matrix_product[ii] + dim[COLUMNS2]*i + j), MPFR_RNDN); 
			
			if(!(matrix_product[ii] = realloc(matrix_product[ii], sizeof(mpfr_t)*olddim[ROWS]*olddim[COLUMNS])))
			{
				// Overhead here isn't inevitable because the order is important
			    #pragma omp parallel for num_threads(algebra_units)
				for(ii=0; ii<algebra_units; ++ii)
			        matrixFree(&matrix_product[ii], ((dim_typ2){dim[ROWS], dim[COLUMNS2]}));
			
			    if(tmp2 != 1)
			        #pragma omp parallel for num_threads(algebra_units)
					for(ii=0; ii<algebra_units; ++ii)
			            matrixFree(&matrix2[ii], ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
			
			    if(!tmp2)
			        #pragma omp parallel for num_threads(algebra_units)
					for(ii=0; ii<algebra_units; ++ii)
			            matrixFree(&matrix1[ii], dim);
			            
			    free(matrix1);
			    free(matrix2);
			    free(matrix_product);
			    return;
			}
		}
		matrix_product_dim[ROWS] = olddim[ROWS];
		matrix_product_dim[COLUMNS] = olddim[COLUMNS];

	}
	else
	{
		matrix_product_dim[ROWS] = dim[ROWS];
		matrix_product_dim[COLUMNS] = dim[COLUMNS2];
	}

	printMatrix(stdout, matrix_product, matrix_product_dim);
	
	for(ii=1; ii<algebra_units; ++ii)
    {
        msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
        printMatrix(stdout, &matrix_product[ii], matrix_product_dim);
    }

    // Overhead here isn't inevitable because the order is important
    #pragma omp parallel for num_threads(algebra_units)
	for(ii=0; ii<algebra_units; ++ii)
        matrixFree(&matrix_product[ii], matrix_product_dim);

    if(tmp2 != 1)
        #pragma omp parallel for num_threads(algebra_units)
		for(ii=0; ii<algebra_units; ++ii)
		{
            matrixFree(&matrix2[ii], ((dim_typ2){dim[COLUMNS], dim[COLUMNS2]}));
            matrixFree(&matrix1[ii], dim);
        }
            
    free(matrix1);
    free(matrix2);
    free(matrix_product);
    return;
}

__MSSHELL_WRAPPER_ static void  matrixKProduct(const sel_typ argc, char ** argv)
{

    dim_typ dim[2];
    dim_typ cdim[2];
    dim_typ i, j;
    const bool complex_entries = access(curLayout)->algebra != 0;
    const sel_typ algebra_units = exp2(access(curLayout)->algebra);
    mpfr_t **matrix1 = malloc(sizeof(mpfr_t*)*algebra_units);
    errMem(matrix1, VSPACE);

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
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
        msprintf(COLOR_CREDITS, "Enter IMAGINARY %s PART of the First Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
        msprintf(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[ROWS], dim[COLUMNS]);
        if(argc > argv_pos)
        {
        	if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            	strcat(argv[argv_pos], TERMINATING_STRING);
            if((!matrixToken(argv[argv_pos], &matrix1[i], cdim, &cdim[COLUMNS])) || cdim[ROWS] != dim[ROWS] || cdim[COLUMNS] != dim[COLUMNS])
            {
                #pragma omp parallel for
				for(j=0; j<=i; ++j)
                    matrixFree(&matrix1[j], dim);
                free(matrix1);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {
            if(!insertNMMatrix(&matrix1[i], dim))
            {
                matrixFree(matrix1, dim);
                #pragma omp parallel for
				for(j=1; j<i; ++j)
                    matrixFree(&matrix1[j], dim);
                free(matrix1);
                return;
            }
        }
    }

    if(complex_entries)
    {
        msprintf(COLOR_USER, "\n\nFirst Complex Matrix correctly inserted:\n\n");
        printMatrix(stdout, matrix1, dim);
        for(i=1; i<algebra_units; ++i)
        {
            msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            printMatrix(stdout, &matrix1[i], dim);
        }
    }

    mpfr_t **matrix2 = malloc(sizeof(mpfr_t*)*algebra_units);
    errMem(matrix2, (free_foreach(matrix1, algebra_units, (dim_typ2){dim[ROWS], dim[COLUMNS]}, NORMAL_MODE), free(matrix1)));
    dim_typ dim2[2];
    const bool assert_m = __pmode__ == ALGOPS_MATRIXKPOWER;

    ++ argv_pos;
    uint64_t tmp = 0;

    if(assert_m)
    {
        mpfr_t tmp2;
        if(argc > argv_pos)
        {
        	if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            	strcat(argv[argv_pos], TERMINATING_STRING);
            if((!parse(argv[argv_pos], &tmp2)) || mpfr_cmp_ui(tmp2, (tmp = mpfr_get_ui(tmp2, MPFR_RNDN))) || tmp < 1)
            {
                #pragma omp parallel for num_threads(algebra_units)
				for(i=0; i<algebra_units; ++i)
                    matrixFree(&matrix1[i], dim);
                mpfr_clear(tmp2);
                free(matrix1);
                free(matrix2);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {
            msprintf(COLOR_CREDITS, "\nEnter N Power to which you want to raise the Inserted Matrix.\n");
            PRINTHOWTOBACKMESSAGE();
            PRINT2N();
			
            while(requires(tmp2, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp2) || mpfr_cmp_ui(tmp2, (tmp = mpfr_get_ui(tmp2, MPFR_RNDN))) || tmp < 1)
            {
                CLEARBUFFER();
                if(isNullVal(tmp2) && exitHandleCheck)
                {
                    #pragma omp parallel for num_threads(algebra_units)
					for(i=0; i<algebra_units; ++i)
                        matrixFree(&matrix1[i], dim);
                    mpfr_clear(tmp2);
                    free(matrix1);
                    free(matrix2);
                    return;
                }
                mpfr_clear(tmp2);
                printErr(1, "Invalid inserted Value.\nExponent must be a non-negative and >= 1 integer");
            }

            CLEARBUFFER();

            if(tmp != 1)
            {
            	const dim_typ tmpdim[2] =
            	{
            		powi(dim[ROWS], tmp),
            		powi(dim[COLUMNS], tmp)
            	};
            	
            	for(i=0; i<algebra_units; ++i) 
            	{
	            	if(!matrixAlloc(&matrix2[i], tmpdim))
	                {
	                    #pragma omp parallel for num_threads(algebra_units)
						for(i=0; i<algebra_units; ++i)
	                        matrixFree(&matrix1[i], dim);
	                    #pragma omp parallel for num_threads(i)
						for(j=0; j<i; ++j)
                      	 	matrixFree(&matrix2[j], dim);
	                    mpfr_clear(tmp2);
	                    free(matrix1);
	                    free(matrix2);
	                    return;
	                }
	             (void) equalMatrix(&matrix2[i], matrix1[i], dim, EQUALMATRIX_NOREALLOC);
	            }
            	
            	/*
                if(!equalMatrix(matrix2, (*matrix1), dim, EQUALMATRIX_REALLOC))
                {
                    #pragma omp parallel for num_threads(algebra_units)
					for(i=0; i<algebra_units; ++i)
                        matrixFree(&matrix1[i], dim);
                    mpfr_clear(tmp2);
                    free(matrix1);
                    free(matrix2);
                    return;
                }
                
                #pragma omp parallel for num_threads(algebra_units)
				for(i=1; i<algebra_units; ++i)
                    if(!equalMatrix(&matrix2[i], matrix1[i], dim, EQUALMATRIX_REALLOC))
                    {
                    	#pragma omp parallel for num_threads(algebra_units)
						for(i=0; i<algebra_units; ++i)
                      	 	matrixFree(&matrix1[i], dim);
                      	#pragma omp parallel for num_threads(i)
						for(j=0; j<i; ++j)
                      	 	matrixFree(&matrix2[j], dim);
	                    mpfr_clear(tmp2);
	                    free(matrix1);
	                    free(matrix2);
                    }
                */

            }
            dim2[ROWS] = dim[ROWS];
            dim2[COLUMNS] = dim[COLUMNS];
        }
    }
    else
    {
        if(argc > argv_pos)
        {
			if(argv[argv_pos][strlen(argv[argv_pos])-1] != TERMINATING_CHAR)
            	strcat(argv[argv_pos], TERMINATING_STRING);
            if((!matrixToken(argv[argv_pos], matrix2, dim2, &dim2[COLUMNS])))
            {
                #pragma omp parallel for num_threads(algebra_units)
				for(i=0; i<algebra_units; ++i)
                    matrixFree(&matrix1[i], dim);
                matrixFree(matrix2, dim2);
                free(matrix1);
                free(matrix2);
                printUsage(&alg_operations[__pmode__]);
                return;
            }
        }
        else
        {

            if(!enterMatrix(matrix2, dim2, &dim2[COLUMNS], false, !complex_entries))
            {
                #pragma omp parallel for num_threads(algebra_units)
				for(i=0; i<algebra_units; ++i)
                    matrixFree(&matrix1[i], dim);
                free(matrix1);
                free(matrix2);
                return;
            }
        }

        for(i=1; i<algebra_units; ++i, ++argv_pos)
        {
            matrix2[i] = NULL;
            msprintf(COLOR_CREDITS, "Enter IMAGINARY %s PART of the Second Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
            msprintf(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim2[ROWS], dim2[COLUMNS]);
            if(argc > argv_pos)
            {
                if((!matrixToken(argv[argv_pos], &matrix2[i], cdim, &cdim[COLUMNS])) || cdim[ROWS] != dim2[ROWS] || cdim[COLUMNS] != dim2[COLUMNS])
                {
                    #pragma omp parallel for
					for(j=0; j<=i; ++j)
                    {
                        matrixFree(&matrix1[j], dim);
                        matrixFree(&matrix2[j], dim2);
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
                        matrixFree(&matrix1[j], dim);
                        matrixFree(&matrix2[j], dim2);
                    }
                    free(matrix1);
                    free(matrix2);
                    return;
                }
            }
        }

        if(complex_entries)
        {
            msprintf(COLOR_USER, "\nSecond Complex Matrix correctly inserted:\n\n");
            printMatrix(stdout, matrix2, dim2);
            for(i=1; i<algebra_units; ++i)
            {
                msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
                printMatrix(stdout, &matrix2[i], dim2);
            }
        }
    }

    const dim_typ dim3[2] =
    {
        assert_m ? powi(dim[ROWS], tmp) : dim[ROWS]*dim2[ROWS],
        assert_m ? powi(dim[COLUMNS], tmp) : dim[COLUMNS]*dim2[COLUMNS]
    };

    PRINT2N();

    if(assert_m)
        msprintf(COLOR_SYSTEM, "First Matrix [%hu x %hu] raised to the Kronecker Power of: %llu is the [%hu x %hu] Matrix:", dim[ROWS], dim[COLUMNS], tmp, dim3[ROWS], dim3[COLUMNS]);
    else
        msprintf(COLOR_SYSTEM, "Matrix Kronecker Product [%hu x %hu] between First and Second Matrix is:", dim3[ROWS], dim3[COLUMNS]);

    PRINT2N();

    dim_typ x = 1;
    mpfr_t **matrix_product = malloc(sizeof(mpfr_t*)*algebra_units);
    errMem(matrix_product, (free_foreach(matrix1, algebra_units, (dim_typ2){dim[ROWS], dim[COLUMNS]}, NORMAL_MODE),
                            free_foreach(matrix2, algebra_units, (dim_typ2){dim[COLUMNS], dim[COLUMNS2]}, tmp == 1),
                            free(matrix1),
                            free(matrix2)));

    bool mxprodErr = false;

    if(tmp != 1)
		for(i=0; i<algebra_units; ++i)
            if((mxprodErr = !matrixAlloc(&matrix_product[i], dim3)))
                break;

    dim_typ dims[2][2];
    
    dims[FIRST_MATRIX][ROWS] = dim[ROWS];
    dims[FIRST_MATRIX][COLUMNS] = dim[COLUMNS];
    
    
    static void (* const matrixKProdFuncs[_MAX_ALGEBRA])(mpfr_t **, mpfr_t **, mpfr_t **, register dim_typ [static 2][2]) =
    {
    	_matrixKProduct,
        _matrixKCProduct,
        _matrixKQProduct,
        _matrixKOProduct,
        _matrixKSProduct
    };
    
    if(tmp != 1)
    {
        if(mxprodErr && !tmp)
        {
            #pragma omp parallel for num_threads(algebra_units)
		    for(i=0; i<algebra_units; ++i)
            {
                matrixFree(&matrix_product[i], dim3);
                matrixFree(&matrix1[i], dim);
            }
            free(matrix1);
            free(matrix2);
            free(matrix_product);
            return;
        }
        struct timeval tvBegin;
        const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
        
        if(difftime) 
        	gettimeofday(&tvBegin, NULL); 
  		do
        {
        	
			dims[SECOND_MATRIX][ROWS] = powi(dim2[ROWS], x);
			dims[SECOND_MATRIX][COLUMNS] = powi(dim2[COLUMNS], x);
            matrixKProdFuncs[access(curLayout)->algebra](matrix1, matrix2, matrix_product, dims);
            
            if(assert_m)
				for(i=0; i<algebra_units; ++i)
                 (void) equalMatrix(&matrix2[i], matrix_product[i], (dim_typ2){dims[SECOND_MATRIX][ROWS]*x, dims[SECOND_MATRIX][COLUMNS]*x}, EQUALMATRIX_NOREALLOC);
        }
        while(++x < tmp);
        
        if(difftime)
    	{
	    	PRINTL();
	        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}
		
    }
    else
        #pragma omp parallel for num_threads(algebra_units)
		for(i=0; i<algebra_units; ++i)
            matrix_product[i] = matrix1[i];

    PRINTL();
   	printMatrix(stdout, matrix_product, dim3);

    for(i=1; i<algebra_units; ++i)
    {
        msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][i]);
        printMatrix(stdout, &matrix_product[i], dim3);
    }
    // Freeing Matrix (Real)
    // Overhead here isn't inevitable because the order is important

    #pragma omp parallel for num_threads(algebra_units)
	for(i=0; i<algebra_units; ++i)
        matrixFree(&matrix_product[i], dim3);

    if(tmp != 1)
        #pragma omp parallel for num_threads(algebra_units)
		for(i=0; i<algebra_units; ++i)
		{
            matrixFree(&matrix2[i], assert_m ? dim3 : dim2);
            matrixFree(&matrix1[i], dim);
        }

    free(matrix1);
    free(matrix2);
    free(matrix_product);
    return;
}


__MSSHELL_WRAPPER_ static void  perScalarMultiplication(const sel_typ argc, char ** argv)
{
    mpfr_t scal;
    dim_typ ii, ij;
    mpfr_t **matrix;

    const bool complex_entries = access(curLayout)->algebra != 0;
    const sel_typ algebra_units = exp2(access(curLayout)->algebra);
    dim_typ dim[2];
    dim_typ cdim[2];

    matrix = malloc(sizeof(mpfr_t)*algebra_units);
    errMem(matrix, VSPACE);

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if(!matrixToken(argv[0], matrix, dim, &dim[COLUMNS]))
        {
        	matrixFree(matrix, dim);
            free(matrix);
        	printUsage(&alg_operations[ALGOPS_PERSCALARMULTIPLICATION]);
            return;
    	}
    }
    else if(!enterMatrix(matrix, dim, &dim[COLUMNS], false, !complex_entries))
        return;

    for(ii=1; ii<algebra_units; ++ii)
    {
        matrix[ii] = NULL;
        msprintf(COLOR_CREDITS, "Enter IMAGINARY %s PART of the Complex Matrix\n", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
        msprintf(COLOR_CREDITS, "by entering another [%hu x %hu] Matrix.\n\n", dim[ROWS], dim[COLUMNS]);
        if(argc > 1)
        {
            if((!matrixToken(argv[1], &matrix[ii], cdim, &cdim[COLUMNS])) || cdim[ROWS] != dim[ROWS] || cdim[COLUMNS] != dim[COLUMNS])
            {  
				#pragma omp parallel for
				for(ij=0; ij<=ii; ++ij)
                    matrixFree(&matrix[ij], dim);
                printUsage(&alg_operations[ALGOPS_PERSCALARMULTIPLICATION]);
                return;
            }
        }
        else
        {
            if(!insertNMMatrix(&matrix[ii], dim))
            {
                #pragma omp parallel for
				for(ij=0; ij<ii; ++ij)
                    matrixFree(&matrix[ij], dim);
                return;
            }
        }
    }

    if(complex_entries)
    {
        msprintf(COLOR_USER, "\nMatrix correctly inserted:\n\n");
        printMatrix(stdout, matrix, dim);
        for(ii=1; ii<algebra_units; ++ii)
        {
            msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
            printMatrix(stdout, &matrix[ii], dim);
        }
    }

    PRINTHOWTOBACKMESSAGE();
    PRINTN();

    const sel_typ argv_pos = 1+algebra_units;

    if(argc > argv_pos)
    {
        if(!parse(argv[argv_pos], &scal))
        {
        	mpfr_clear(scal);
            #pragma omp parallel for num_threads(algebra_units)
			for(ii=0; ii<algebra_units; ++ii)
                matrixFree(&matrix[ii], dim);
            return;
        }
    }
    else
        while(requires(scal, NULL, "Enter a double floating-point Scalar Number.\n", "Inserted Scalar", PARSER_SHOWRESULT) || (dcheck && INVERSE_OPS && !scal))
        {
        	printErr(22, "Invalid inserted Value");
        	mpfr_clear(scal);
        }

    if(isNullVal(scal))
    {
        #pragma omp parallel for num_threads(algebra_units)
		for(ii=0; ii<algebra_units; ++ii)
            matrixFree(&matrix[ii], dim);
        mpfr_clear(scal); 
        return;
    }
    CLEARBUFFER();
    char buf[MAX_BUFSIZ];
    
    if(__pmode__ == ALGOPS_SCALARDIVISIONMATRIX)
        mpfr_sprintf(buf, "Scalar %Rf has been correctly divided by the Matrix:\n\n", scal);
    else
        mpfr_sprintf(buf, "Matrix has been correctly %s by the Scalar %Rf.\n\n", INVERSE_OPS ? "divided" : "multiplied", scal);
    
    msprintf(COLOR_SYSTEM, buf);

    PRINTL();

    if(INVERSE_OPS && __pmode__ == ALGOPS_PERSCALARMULTIPLICATION)
        mpfr_pow_si(scal, scal, -1, MPFR_RNDN);
	
    dim_typ i, j, k;
    // int (* mul_func)(mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t) = __pmode__ == ALGOPS_PERSCALARMULTIPLICATION ? mpfr_mul : mpfr_div;

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	const register dim_typ rows = dim[ROWS];
	const register dim_typ columns = dim[COLUMNS];
	
	// I cannot use mul_func pointer due to a stupid nonsense warning
		
	
	if(__pmode__ == ALGOPS_PERSCALARMULTIPLICATION)
		#pragma omp parallel for num_threads(algebra_units)
	    for(k=0; k<algebra_units; ++k)
	    	#pragma omp parallel for
	        for(i=0; i<rows; ++i)
	        	#pragma omp parallel for
	            for(j=0; j<columns; ++j)
	             (void) mpfr_mul(matrix[k][dim[COLUMNS]*i+j], scal, matrix[k][dim[COLUMNS]*i+j], MPFR_RNDN);
    else
	    #pragma omp parallel for num_threads(algebra_units)
	    for(k=0; k<algebra_units; ++k)
	    	#pragma omp parallel for
	        for(i=0; i<rows; ++i)
	        	#pragma omp parallel for
	            for(j=0; j<columns; ++j)
	             (void) mpfr_div(matrix[k][dim[COLUMNS]*i+j], scal, matrix[k][dim[COLUMNS]*i+j], MPFR_RNDN);
	                
	if(difftime)
	{
    	PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	}

    printMatrix(stdout, matrix, dim);
    mpfr_clear(scal);
    matrixFree(matrix, dim);
    for(ii=1; ii<algebra_units; ++i)
    {
        msprintf(COLOR_USER, "\n+ %s*", suite_c.algebra_imaginary_units_names[access(curLayout)->algebra][ii]);
        printMatrix(stdout, &matrix[ii], dim);
        matrixFree(&matrix[ii], dim);
    }
    return;
}

__MSSHELL_WRAPPER_ static void  illConditionChecking(const sel_typ argc, char ** argv)
{
    mpfr_t *matrix;
    dim_typ dim[2];

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
        	matrixFree(&matrix, dim);
        	printUsage(&alg_operations[ALGOPS_ILLCONDITIONCHECKING]);
            return;
    	}
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], true))
        return;
        
    mpfr_t *matrix2;

	if(!matrixAlloc(&matrix2, dim))
    	return;

	dim_typ i;
	mpfr_t norms1, norms2; 
	mpfr_init(norms1);
	norms(norms1, matrix2, dim[ROWS]);
    const dim_typ dimrows_per2 = dim[ROWS]<<1;
	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
    
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
	const register sel_typ exitVal = invertMatrix(matrix, matrix2, dim[ROWS]);
		
    if(exitVal == INVERTMATRIX_SINGULAR)
        printErr(1, "You cannot test SINGULAR Matrices");
    else if(exitVal == INVERTMATRIX_ALLOCERROR)
    	printErr(12, "Matrix Inverse Calculator Dynamic Memory Allocation Problem");
    else
    {
    	if(difftime)
		{
	    	PRINTL();
	        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}
    	char buf[MAX_BUFSIZ];
		mpfr_init(norms2);
		norms(norms2, matrix, dim[ROWS]);
		// norms(norms2, matrix2, dim[ROWS]);
		mpfr_mul(norms2, norms1, norms2, MPFR_RNDN); 
	    mpfr_sprintf(buf, "\nInserted Square Matrix ILL CONDITION CHECKING is: %Rf.\n\n", norms2);
	    msprintf(COLOR_USER, buf);
	    mpfr_clear(norms2);
    }
	
    matrixFree(&matrix, dim);
    matrixFree(&matrix2, dim);
    mpfr_clear(norms1); 
    return;
}

__MSSHELL_WRAPPER_ static void  matrixFattLU(const sel_typ argc, char ** argv)
{
    dim_typ dim[2];
    mpfr_t *matrix;
    mpfr_t *L, *U;

    if(argc)
    {
    	if(argv[0][strlen(argv[0])-1] != TERMINATING_CHAR)
            strcat(argv[0], TERMINATING_STRING);
        if((!matrixToken(argv[0], &matrix, dim, &dim[COLUMNS])) || dim[ROWS] != dim[COLUMNS])
        {
            matrixFree(&matrix, dim);
            printUsage(&alg_operations[ALGOPS_MATRIXLUFACTORIZATION]);
            return;
        }
    }
    else if(!insertMatrix(matrix, dim[ROWS], dim[COLUMNS], true))
        return;

    bool assert[2];

    if(((assert[LOWER_TRIANGULAR] = !matrixAlloc(&L, dim))) || ((assert[UPPER_TRIANGULAR] = !matrixAlloc(&U, dim))))
    {
        matrixFree(&matrix, dim);
        if(assert[LOWER_TRIANGULAR] && !assert[UPPER_TRIANGULAR])
            matrixFree(&L, dim);
        return;
    }

	struct timeval tvBegin;
	const bool difftime = isSett(BOOLS_SHOWDIFFTIME);
	
	if(difftime)
		gettimeofday(&tvBegin, NULL);
		
    if(!FattLU(dim[ROWS], matrix, L, U))
        printErr(1, "You cannot LU-decompose SINGULAR Matrices");
    else
    {
    	if(difftime)
		{
			PRINTL();
		    msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
		}
        msprintf(COLOR_USER, "\nInserted Matrix has been correctly LU-decomposed\ninto Triangular Matrices respectively: L as Lower and U as Upper:\n\n");
        printMatrix(stdout, &L, dim);
        printMatrix(stdout, &U, dim);
    }

    matrixFree(&matrix, dim);
    matrixFree(&L, dim);
    matrixFree(&U, dim);
    return;
}
