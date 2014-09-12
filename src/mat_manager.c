// mat_manager.c 10/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifdef ALLOW_MATMANAGER

// FUNCTIONS DECLARATIONS

// static bool extractMat(ityp ***matrix, dim_typ *righe, dim_typ *colonne, char item_path[MAX_PATH_LENGTH]);

__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog setCurMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __lmp_prog createMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog viewMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog printMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog updMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog updAllMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog delMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog editMat(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog renMat(const sel_typ argc, char ** argv);

sprog mat_manager[MAX_MATMANAGER_PROGS] =
{
    [MATRICES_SETCURRENT] =
    {
        "Select Current Matrix",
        CMD_SETCURMAT,
        USAGE_SETCURMAT,
        setCurMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_OPEN] =
    {
        "Open Matrix ("DEFAULT_MATRIX_FILE_EXTENSION")",
        CMD_OPENMAT,
        USAGE_OPENMAT,
        createMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_CREATE] =
    {
        "Create new Matrix ("DEFAULT_MATRIX_FILE_EXTENSION")",
        CMD_CREATEMAT,
        USAGE_CREATEMAT,
        createMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_READ] =
    {
        "View Matrix",
        CMD_READMAT,
        USAGE_READMAT,
        viewMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_PRINT] =
    {
        "Print Matrix",
        CMD_PRTMAT,
        USAGE_PRTMAT,
        printMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_UPDATE] =
    {
        "Save Matrix",
        CMD_UPDMAT,
        USAGE_UPDMAT,
        updMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_UPDATEALL] =
    {
        "Save all Matrices",
        CMD_UPDALLMAT,
        USAGE_UPDALLMAT,
        updAllMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_SAVE] =
    {
        "Save Matrix as ("DEFAULT_MATRIX_FILE_EXTENSION")",
        CMD_SAVEMAT,
        USAGE_SAVEMAT,
        createMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_DELETE] =
    {
        "Delete Matrix",
        CMD_DELMAT,
        USAGE_DELMAT,
        delMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_DELETEPHYSICAL] =
    {
        "Delete Physical Matrix",
        CMD_DELMAT2,
        USAGE_DELMAT2,
        delMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_DELETEALL] =
    {
        "Delete all Matrices",
        CMD_DELALLMATS,
        USAGE_DELALLMATS,
        delMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_DELETEALLPHYSICALS] =
    {
        "Delete all Physical Matrices",
        CMD_DELALLMATS2,
        USAGE_DELALLMATS2,
        delMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_EDIT] =
    {
        "Edit Matrix",
        CMD_EDITMAT,
        USAGE_EDITMAT,
        editMat,
        AUTOMATIC,
        CHILD
    },
    [MATRICES_RENAME] =
    {
        "Rename Matrix ("DEFAULT_MATRIX_FILE_EXTENSION")",
        CMD_RENMAT,
        USAGE_RENMAT,
        renMat,
        AUTOMATIC,
        CHILD
    }
};

// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog setCurMat(const sel_typ argc, char ** argv)
{
    setCurItem(argc ? getItemID(argv[0], &mat_manager[MATRICES_SETCURRENT], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static inline void _MS__private __lmp_prog createMat(const sel_typ argc, char ** argv)
{
    createItem(argc ? argv[0] : NULL, MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog viewMat(const sel_typ argc, char ** argv)
{
    viewItem(argc ? getItemID(argv[0], &mat_manager[MATRICES_READ], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog printMat(const sel_typ argc, char ** argv)
{
    printListItem(argc ? getItemID(argv[0], &mat_manager[MATRICES_PRINT], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog updMat(const sel_typ argc, char ** argv)
{
    updItem(MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog updAllMat(const sel_typ argc, char ** argv)
{
    updAll(MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog delMat(const sel_typ argc, char ** argv)
{
    delItem(argc ? getItemID(argv[0], &mat_manager[__pmode__], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog editMat(const sel_typ argc, char ** argv)
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
        ityp *matrix = NULL;
        dim_typ olddim[MAX_DIMENSIONS];

        olddim[ROWS] = tmp->dim[ROWS];
        olddim[COLUMNS] = tmp->dim[COLUMNS];

        if(argc > 1)
        {
            if(!matrixToken(argv[1], &matrix, &(tmp->dim[ROWS]), &(tmp->dim[COLUMNS])))
            {
                matrixFree(&matrix);
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

        if(!equalMatrix(&(tmp->matrix), matrix, tmp->dim))
        {
            matrixFree(&matrix);
            #if WINOS
                SetExitButtonState(ENABLED);
            #endif
            return;
        }
    }
    else if(!insertMatrix(tmp->matrix, tmp->dim[ROWS], tmp->dim[COLUMNS], false))
        return;

    if(saveItem(which_mat, MATRICES))
    {
        printf2(COLOR_USER, "\nSelected Matrix:\n");
        printf2(COLOR_CREDITS, listNo(which_mat, MATRICES)->path);
        printf2(COLOR_USER, "\nhas been successfully modified.\n\n");
    }

    /// passToItem(old_mat, MATRICES, false);

    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog renMat(const sel_typ argc, char ** argv)
{
    renItem(argc > 1 ? argv[1] : NULL, argc ? getItemID(argv[0], &mat_manager[MATRICES_RENAME], MATRICES) : getItemsListNo(MATRICES), MATRICES);
    return;
}

#endif
