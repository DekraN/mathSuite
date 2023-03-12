// mat_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifndef __DISABLE_MATMANAGER

// FUNCTIONS DECLARATIONS

// static bool extractMat(ityp ***matrix, dim_typ *righe, dim_typ *colonne, char item_path[MAX_PATH_LENGTH]);

__MSSHELL_WRAPPER_ static void  __lmp_prog setCurMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void  __lmp_prog createMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog viewMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog printMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog updMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog updAllMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog delMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog editMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog renMat(const sel_typ argc, char ** argv);
#ifndef __DISABLE_DATABASE
	__MSSHELL_WRAPPER_ static void  __lmp_prog persistMat(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void  __lmp_prog retrieveMat(const sel_typ argc, char ** argv);
#endif

sprog mat_manager[MAX_MATMANAGER_PROGS] =
{
    [MATRICES_SETCURRENT] =
    {
    	CMD_SETCURMAT,
        NAME_SETCURMAT,
        USAGE_SETCURMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_SETCURMAT,
        #endif
        setCurMat,
        ARGC_SETCURMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_OPEN] =
    {
    	CMD_OPENMAT,
        NAME_OPENMAT,
        USAGE_OPENMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_OPENMAT,
        #endif
        createMat,
        ARGC_OPENMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_CREATE] =
    {
    	CMD_CREATEMAT,
        NAME_CREATEMAT,
        USAGE_CREATEMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_CREATEMAT,
        #endif
        createMat,
        ARGC_CREATEMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_READ] =
    {
    	CMD_READMAT,
        NAME_READMAT,
        USAGE_READMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_READMAT,
        #endif
        viewMat,
        ARGC_READMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_PRINT] =
    {
    	CMD_PRTMAT,
        NAME_PRTMAT,
        USAGE_PRTMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_PRTMAT,
        #endif
        printMat,
        ARGC_PRTMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_UPDATE] =
    {
    	CMD_UPDMAT,
        NAME_UPDMAT,
        USAGE_UPDMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDMAT,
        #endif
        updMat,
        ARGC_UPDMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_UPDATEALL] =
    {
    	CMD_UPDALLMAT,
        NAME_UPDALLMAT,
        USAGE_UPDALLMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDALLMAT,
        #endif
        updAllMat,
        ARGC_UPDALLMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_SAVE] =
    {
    	CMD_SAVEMAT,
        NAME_SAVEMAT,
        USAGE_SAVEMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_SAVEMAT,
        #endif
        createMat,
        ARGC_SAVEMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_DELETE] =
    {
    	CMD_DELMAT,
        NAME_DELMAT,
        USAGE_DELMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELMAT,
        #endif
        delMat,
        ARGC_DELMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_DELETEPHYSICAL] =
    {
    	CMD_DELMAT2,
        NAME_DELMAT2,
        USAGE_DELMAT2,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELMAT2,
        #endif
        delMat,
        ARGC_DELMAT2,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_DELETEALL] =
    {
    	CMD_DELALLMATS,
        NAME_DELALLMATS,
        USAGE_DELALLMATS,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELALLMATS,
        #endif
        delMat,
        ARGC_DELALLMATS,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_DELETEALLPHYSICALS] =
    {
    	CMD_DELALLMATS2,
        NAME_DELALLMATS2,
        USAGE_DELALLMATS2,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELALLMATS2,
        #endif
        delMat,
        ARGC_DELALLMATS2,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_EDIT] =
    {
    	CMD_EDITMAT,
        NAME_EDITMAT,
        USAGE_EDITMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_EDITMAT,
        #endif
        editMat,
        ARGC_EDITMAT,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_RENAME] =
    {
    	CMD_RENMAT,
        NAME_RENMAT,
        USAGE_RENMAT,
        #ifndef __DISABLE_DATABASE
        LEVEL_RENMAT,
        #endif
        renMat,
        ARGC_RENMAT,
        AUTOMATIC,
        CHILD
    }
    #ifndef __DISABLE_DATABASE
    , [MATRICES_PERSIST] =
    {
    	CMD_PERSISTMAT,
    	NAME_PERSISTMAT,
    	USAGE_PERSISTMAT,
    	LEVEL_PERSISTMAT,
    	persistMat,
    	ARGC_PERSISTMAT,
    	AUTOMATIC,
    	CHILD
    },
    [MATRICES_RETRIEVE] =
    {
    	CMD_RETRIEVEMAT,
    	NAME_RETRIEVEMAT,
    	USAGE_RETRIEVEMAT,
    	LEVEL_RETRIEVEMAT,
    	retrieveMat,
    	ARGC_RETRIEVEMAT,
    	AUTOMATIC,
    	CHILD
    }
    #endif
};

// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static inline void  __lmp_prog setCurMat(const sel_typ argc, char ** argv)
{
    setCurItem(argc ? getItemID(argv[0], &mat_manager[MATRICES_SETCURRENT], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static inline void  __lmp_prog createMat(const sel_typ argc, char ** argv)
{
    createItem(argc ? argv[0] : NULL, MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog viewMat(const sel_typ argc, char ** argv)
{
    viewItem(argc ? getItemID(argv[0], &mat_manager[MATRICES_READ], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog printMat(const sel_typ argc, char ** argv)
{
    printListItem(argc ? getItemID(argv[0], &mat_manager[MATRICES_PRINT], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog updMat(const sel_typ argc, char ** argv)
{
    updItem(MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog updAllMat(const sel_typ argc, char ** argv)
{
    updAll(MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog delMat(const sel_typ argc, char ** argv)
{
    delItem(argc ? getItemID(argv[0], &mat_manager[__pmode__], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static void  __lmp_prog editMat(const sel_typ argc, char ** argv)
{
    dim_typ which_mat;

    if((which_mat = argc ? getItemID(argv[0], &mat_manager[MATRICES_EDIT], MATRICES) : itemSelect(MATRICES)) == NULL_ITEM(MATRICES))
        return;

    // const dim_typ old_mat = suite.lists[MATRICES].cur_item;
    // old_mat = suite.lists[MATRICES].cur_item;
    /// passToItem(which_mat, MATRICES, true);

    matrixObj * const tmp = ((matrixObj *)(listNo(which_mat, MATRICES)->data));

    if(tmp->matrix)
    {
        EXPRTYPE *matrix;
        dim_typ olddim[2];

        olddim[ROWS] = tmp->dim[ROWS];
        olddim[COLUMNS] = tmp->dim[COLUMNS];

        if(argc > 1)
        {
            if(!matrixToken(argv[1], &matrix, &(tmp->dim[ROWS]), &(tmp->dim[COLUMNS])))
            {
                matrixFree(&matrix, ((dim_typ2){tmp->dim[ROWS], tmp->dim[COLUMNS]}));
                tmp->dim[ROWS] = olddim[ROWS];
                tmp->dim[COLUMNS] = olddim[COLUMNS];
                printUsage(&mat_manager[MATRICES_EDIT]);
                return;
            }
        }
        else if(!insertMatrix(matrix, tmp->dim[ROWS], tmp->dim[COLUMNS], false))
        {
            tmp->dim[ROWS] = olddim[ROWS];
            tmp->dim[COLUMNS] = olddim[COLUMNS];
            return;
        }

		matrixFree(&(tmp->matrix), olddim);
        if(!equalMatrix(&(tmp->matrix), matrix, tmp->dim, EQUALMATRIX_REALLOC))
            return;
    }
    else if(!insertMatrix(tmp->matrix, tmp->dim[ROWS], tmp->dim[COLUMNS], false))
        return;

    if(saveItem(which_mat, MATRICES))
    {
        msprintf(COLOR_USER, "\nSelected Matrix:\n");
        msprintf(COLOR_CREDITS, listNo(which_mat, MATRICES)->path);
        msprintf(COLOR_USER, "\nhas been successfully modified.\n\n");
    }

    /// passToItem(old_mat, MATRICES, false);

    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog renMat(const sel_typ argc, char ** argv)
{
    renItem(argc > 1 ? argv[1] : NULL, argc ? getItemID(argv[0], &mat_manager[MATRICES_RENAME], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

#ifndef __DISABLE_DATABASE
	__MSSHELL_WRAPPER_ static void  __lmp_prog persistMat(const sel_typ argc, char ** argv)
	{
	    dim_typ which_mat;
	
	    if((which_mat = argc ? getItemID(argv[0], &mat_manager[MATRICES_EDIT], MATRICES) : itemSelect(MATRICES)) == NULL_ITEM(MATRICES))
	        return;
	
	    matrixObj * const tmp = ((matrixObj *)(listNo(which_mat, MATRICES)->data));
		            
		if(!access(curLayout)->database.con)
			dbEstablishConnection(0, NULL);
		
		if(access(curLayout)->database.con)
		{
			char matname[MAX_PATH_LENGTH];
			strcpy(matname, listNo(which_mat, MATRICES)->path);
			(void) strtok(matname, EXTENSION_DOT);
			deleteMat(matname, DELETEMAT_ROW);
			_printMatrix(NULL, &(tmp->matrix), tmp->dim, matname); // which_item);
			msprintf(COLOR_USER, "\nIt has been stored in DB the %s:\n%s."DEFAULT_MATRIX_FILE_EXTENSION, suite_c.listsnames[MATRICES], matname);
		}
	
	    return;
	}
	
	__MSSHELL_WRAPPER_ static void  __lmp_prog retrieveMat(const sel_typ argc, char ** argv)
	{
		dim_typ which_mat;
	
	    if((which_mat = argc ? getItemID(argv[0], &mat_manager[MATRICES_EDIT], MATRICES) : itemSelect(MATRICES)) == NULL_ITEM(MATRICES))
	        return;
	        
	    matrixObj * const tmp = malloc(sizeof(matrixObj));
	    errMem(tmp, VSPACE);
	        
		if(!access(curLayout)->database.con)
	    	dbEstablishConnection(0, NULL); // (char * [4]){NULL, NULL, NULL, NULL});
	
		if(access(curLayout)->database.con)
		{
			char mat[MAX_BUFSIZ] = NULL_CHAR;
			char matname [MAX_PATH_LENGTH];
			strcpy(matname, listNo(which_mat, MATRICES)->path);
			(void) strtok(matname, EXTENSION_DOT);
			if((!readMat(matname, mat)) || !matrixToken(mat, &(tmp->matrix), &tmp->dim[ROWS], &tmp->dim[COLUMNS]))
			{
				printErr(ERROR_INPUTOUTPUT, "An error occurred during parsing from the DB the %s:\n%s."DEFAULT_MATRIX_FILE_EXTENSION, suite_c.listsnames[MATRICES], matname);
				free(tmp);
				return;	
			}
			msprintf(COLOR_USER, "\nIt has been parsed from the DB the %s:\n%s."DEFAULT_MATRIX_FILE_EXTENSION, suite_c.listsnames[MATRICES], matname);
		}
	
		listNo(which_mat, MATRICES)->data = tmp;
	    return;
	}
#endif

#endif
