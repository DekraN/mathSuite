#include "dutils.h"

#ifndef __DISABLE_SPECIALMATRICES

__MSSHELL_WRAPPER_ static void  printLmpMatrix(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  printMatrixList(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  printCurrentMLMatrix(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  setCurrentMLMatrix(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  getCurrentMatrix(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  letSpecialMatrix(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  delMatrixList(const sel_typ argc, char ** argv);


sprog specmat_prog[MAX_SPECMAT_PROGS] =
{
	[SPECMAT_PRINTLMPMATRIX] =
    {
    	CMD_PRINTLMPMATRIX,
    	NAME_PRINTLMPMATRIX,
    	USAGE_PRINTLMPMATRIX,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_PRINTLMPMATRIX,
    	#endif
    	printLmpMatrix,
    	ARGC_PRINTLMPMATRIX,
    	AUTOMATIC,
    	CHILD
    },
    [SPECMAT_PRINTMATRIXLIST] =
    {
    	CMD_PRINTMATRIXLIST,
    	NAME_PRINTMATRIXLIST,
    	USAGE_PRINTMATRIXLIST,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_PRINTMATRIXLIST,
    	#endif
    	printMatrixList,
    	ARGC_PRINTMATRIXLIST,
    	AUTOMATIC,
    	CHILD
    },
    [SPECMAT_PRINTCURRENTMLMATRIX] =
    {
    	CMD_PRINTCURRENTMLMATRIX,
    	NAME_PRINTCURRENTMLMATRIX,
    	USAGE_PRINTCURRENTMLMATRIX,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_PRINTCURRENTMLMATRIX,
    	#endif
    	printCurrentMLMatrix,
    	ARGC_PRINTCURRENTMLMATRIX,
    	AUTOMATIC,
    	CHILD
    },
    [SPECMAT_SETCURRENTMLMATRIX] =
    {
    	CMD_SETCURRENTMLMATRIX,
    	NAME_SETCURRENTMLMATRIX,
    	USAGE_SETCURRENTMLMATRIX,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_SETCURRENTMLMATRIX,
    	#endif
    	setCurrentMLMatrix,
    	ARGC_SETCURRENTMLMATRIX,
    	AUTOMATIC,
    	CHILD
    },
    [SPECMAT_GETCURRENTMATRIX] =
    {
    	CMD_GETCURRENTMATRIX,
    	NAME_GETCURRENTMATRIX,
    	USAGE_GETCURRENTMATRIX,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_GETCURRENTMATRIX,
    	#endif
    	getCurrentMatrix,
    	ARGC_GETCURRENTMATRIX,
    	AUTOMATIC,
    	CHILD
    },
    [SPECMAT_DELMATRIXLIST] =
    {
    	CMD_DELMATRIXLIST,
    	NAME_DELMATRIXLIST,
    	USAGE_DELMATRIXLIST,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_DELMATRIXLIST,
    	#endif
    	delMatrixList,
    	ARGC_DELMATRIXLIST,
    	AUTOMATIC,
    	CHILD
    },
    [SPECMAT_LETSPECIALMATRIX] =
    {
    	CMD_LETSPECIALMATRIX,
    	NAME_LETSPECIALMATRIX,
    	USAGE_LETSPECIALMATRIX,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_LETSPECIALMATRIX,
    	#endif
    	letSpecialMatrix,
    	ARGC_LETSPECIALMATRIX,
    	AUTOMATIC,
    	CHILD
    }
};

__MSSHELL_WRAPPER_ static void  printLmpMatrix(const sel_typ argc, char ** argv)
{
	matrixObj * const currentLmpMatrix = accessCurrentLmpMatrix;
	if(currentLmpMatrix && currentLmpMatrix->matrix)	
		printMatrix(stdout, &(currentLmpMatrix->matrix), currentLmpMatrix->dim);
	else
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "No Last Matrix Printed Available");
	return;
}

__MSSHELL_WRAPPER_ static void  printMatrixList(const sel_typ argc, char ** argv)
{
	
	matrixObj * currentMLMatrix = NULL;
	
	if(accessCurrentSession()->MLSystem.list && accessCurrentSession()->MLSystem.itemsno)
		for(dim_typ i=0; i<accessCurrentSession()->MLSystem.itemsno; ++i)
			if((currentMLMatrix = getMatrixFromMatrixList(i)) && currentMLMatrix->matrix)
				printMatrix(stdout, &(currentMLMatrix->matrix), currentMLMatrix->dim);
	else
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "No MatrixList Available");
	return;
}

__MSSHELL_WRAPPER_ static void  printCurrentMLMatrix(const sel_typ argc, char ** argv)
{
	matrixObj * const currentMLMatrix = getMatrixFromMatrixList(accessCurrentSession()->MLSystem.cur_item);
	if(currentMLMatrix && currentMLMatrix->matrix)	
		printMatrix(stdout, &(currentMLMatrix->matrix), currentMLMatrix->dim);
	else
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "No MatrixList Available");
	return;
}


__MSSHELL_WRAPPER_ static void  setCurrentMLMatrix(const sel_typ argc, char ** argv)
{
	session * const currentSession = accessCurrentSession();
	PRINTN();

	if(!currentSession->MLSystem.itemsno)
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "No MatrixList Available");
	else
	{
		if(argv && argv[0])
	    {
	    	mpfr_t tmp;
	    	if(!parse(argv[0], &tmp))
		    {
		    	mpfr_clear(tmp);
		        printErr(1, "Parse Error on "CMD_SETCURRENTMLMATRIX" command.");
		        return;
		    }
		
		    if(mpfr_cmp_ui(tmp, 0) < 0 || mpfr_cmp_ui(tmp, accessCurrentSession()->MLSystem.itemsno) >= 0)
		    {
		    	mpfr_clear(tmp);
		        printErr(33, CMD_SETCURRENTMLMATRIX" accepts only integers between %hu and %hu", 0, currentSession->MLSystem.itemsno-1);
		        printUsage(&specmat_prog[SPECMAT_SETCURRENTMLMATRIX]);
				return;
		    }
		    
		    currentSession->MLSystem.cur_item = mpfr_get_si(tmp, MPFR_RNDN);
		    mpfr_clear(tmp);
	    	msprintf(COLOR_USER, "\nThe Current MatrixList Matrix is: %hu.\n", currentSession->MLSystem.cur_item);
	    }
	    else
	    {
	    	printErr(33, CMD_SETCURRENTMLMATRIX" accepts only integers between %hu and %hu", 0, currentSession->MLSystem.itemsno-1);
		    printUsage(&specmat_prog[SPECMAT_SETCURRENTMLMATRIX]);
			return;
	    }
	}
    
	return;
}

__MSSHELL_WRAPPER_ static void  getCurrentMatrix(const sel_typ argc, char ** argv)
{
	equalSpecialMatrix(SPECIALMATRIX_CURRENTMATRIX, SPECIALMATRIX_LMPMATRIX);
	return;
}

__MSSHELL_WRAPPER_ static void  letSpecialMatrix(const sel_typ argc, char ** argv)
{
	session * const currentSession = accessCurrentSession();
	PRINTN();
	
	if(argc > 1)
    {
    	
    	char destString[MAX_BUFSIZ] = NULL_CHAR;
    	char sourceString[MAX_BUFSIZ] = NULL_CHAR;
    	strcpy(destString, argv[0]);
    	strcpy(sourceString, argv[1]);

    	mpfr_t dest, source; 
    	
    	if(!parse(destString, &dest))
	    {
	    	mpfr_clear(dest);
	        printErr(1, "Parse Error on "CMD_LETSPECIALMATRIX" command.");
	        return;
	    }
	
	    if(mpfr_cmp_ui(dest, 0) < 0 || mpfr_cmp_ui(dest, MAX_SPECIALMATRICES+accessCurrentSession()->MLSystem.itemsno) >= 0)
	    {
	    	mpfr_clear(dest);
	        printErr(33, CMD_LETSPECIALMATRIX" accepts only integers between 0 and %hu", MAX_SPECIALMATRICES+currentSession->MLSystem.itemsno-1);
	        printUsage(&specmat_prog[SPECMAT_LETSPECIALMATRIX]);
			return;
	    }
		    
	    if(!parse(sourceString, &source))
	    {
	    	mpfr_clear(dest);
	    	// mpfr_clear(source);
	        printErr(1, "Parse Error on "CMD_LETSPECIALMATRIX" command.");
	        return;
	    }
	    
	    if(mpfr_cmp_ui(source, 0) < 0 || mpfr_cmp_ui(source, MAX_SPECIALMATRICES+accessCurrentSession()->MLSystem.itemsno) >= 0)
	    {
	    	mpfr_clear(dest);
	    	mpfr_clear(source);
	        printErr(33, CMD_LETSPECIALMATRIX" accepts only integers between 0 and %hu", MAX_SPECIALMATRICES+currentSession->MLSystem.itemsno-1);
	        printUsage(&specmat_prog[SPECMAT_LETSPECIALMATRIX]);
			return;
	    }
	    
	    equalSpecialMatrix(mpfr_get_ui(dest, MPFR_RNDN), mpfr_get_ui(source, MPFR_RNDN));
		mpfr_clears(dest, source, NULL);
    }
    else
    {
    	printErr(33, CMD_LETSPECIALMATRIX" accepts only integers between 0 and %hu", MAX_SPECIALMATRICES+currentSession->MLSystem.itemsno-1);
	    printUsage(&specmat_prog[SPECMAT_LETSPECIALMATRIX]);
		return;
    }
    
	return;
}

__MSSHELL_WRAPPER_ static void  delMatrixList(const sel_typ argc, char ** argv)
{
	delAllMatrixFromMatrixList();
	msyprintf(COLOR_SYSTEM, "\nIt has been deleted the entire Matrix List.\n\n");
	return;
}


#endif
