// lists_manager.c 23/08/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#if((defined ALLOW_VARLISTMANAGER) || (defined ALLOW_MATMANAGER) || (defined ALLOW_LOGSMANAGER) || (defined ALLOW_LAYOUTSMANAGER))

// FUNCTIONS DECLARATIONS

__MSSHELL_WRAPPER_ __MSNATIVE_ static void _MS__private __system printItem(const dim_typ, sel_typ);
__MSSHELL_WRAPPER_ __MSNATIVE_ static void _MS__private __system showCurItem(sel_typ);
__MSSHELL_WRAPPER_ __MSNATIVE_ static bool _MS__private __system listDeleteProc(const char [static MAX_PATH_LENGTH], sel_typ);

// FUNZIONI LISTE

__MSSHELL_WRAPPER_ __MSNATIVE_ static void _MS__private __system showCurItem(sel_typ mode)
{
    char str[MINMIN_STRING] = NULL_CHAR;
    strcpy(str, suite_c.listsnames[mode]);
    toupper_s(str);
    printf2(COLOR_CREDITS, "CURRENT %s:\n%s\n", str, listNo(access(lists)[mode].cur_item, mode)->path);
    return;
}

__MSNATIVE_ inline sel_typ __system __export checkItemTypeByExtension(const char extension[static MAX_EXTENSION_LENGTH])
{

    if(!strcmp(extension, DEFAULT_VARLIST_FILE_EXTENSION))
        return ENVS;

    if(!strcmp(extension, DEFAULT_MATRIX_FILE_EXTENSION))
        return MATRICES;

    if(!strcmp(extension, DEFAULT_LOG_FILE_EXTENSION))
        return LOGS;

    if(!strcmp(extension, DEFAULT_LAYOUT_FILE_EXTENSION))
        return LAYOUTS;

    return MAX_LISTS;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ bool _MS__private __system __export listInsertProc(const char item_path[static MAX_PATH_LENGTH], sel_typ mode)
{
    const bool assert = getItemsListNo(mode) == STARTING_ITEMSNO;

    if(assert)
        initList(mode);

    if(!file_exists(item_path))
    {
        printErr(2, "An error during:\n%s\n%s File opening process might have occurred", item_path, suite_c.listsnames[mode]);
        return false;
    }

    if(!listInsert(item_path, mode)) // if(!headInsert(suite.lists[mode].items_list, str))
    {
        printErr(11, "Failed to add:\n%s\nto the %s Items List", item_path, suite_c.listsnames[mode]);
        return false;
    }

    createItemData(access(lists)[mode].itemsno-1, mode);

    if(assert)
        refreshItem(STARTING_ITEMSNO, mode); // binding BOOLS_AUTOSETCURITEM in this case...

    return true;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ static bool _MS__private __system listDeleteProc(const char item_path[static MAX_PATH_LENGTH], sel_typ mode)
{
    // this will free the envitem SPECIFIED VARLIST

    const dim_typ which_item = searchItem(item_path, mode);
    void * data = listNo(which_item, mode)->data;

    switch(mode)
    {
        case ENVS:
            exprValListFree(((exprType *)data)->var_list);
            break;

        case MATRICES:
        {
            matrixObj * const tmp = ((matrixObj*)data);
            if(tmp->matrix)
                matrixFree(&(tmp->matrix), tmp->dim[RAWS]);
            break;
        }

        case LOGS:
            free(((logObj *)data)->buffer);
            break;

        case LAYOUTS:
            // do something
            break;

        default: return false;
    }

    free(data);

    if(!listDelete(which_item, mode))
    {
        printErr(11, "Failed to delete:\n%s\n%s File", item_path, suite_c.listsnames[mode]);
        return false;
    }

    printf2(COLOR_USER, "\nIt has been correctly deleted the %s File:", suite_c.listsnames[mode]);
    printf2(COLOR_CREDITS, "\n%s\n\n", item_path);
    return true;
}

__MSNATIVE_ dim_typ __system __export searchItem(const char item_path[static MAX_PATH_LENGTH], sel_typ mode)
{
    dim_typ i, j;
    nodelist * pntr[MAX_DIMENSIONS];

    pntr[NEXT_LISTNODE] = access(lists)[mode].items_list[NEXT_LISTNODE];
    pntr[PREV_LISTNODE] = access(lists)[mode].items_list[PREV_LISTNODE];


    for(i=0, j=getItemsListNo(mode)-1; (i < getItemsListNo(mode)) && (j >= 0); ++i, --j)
    {
        if((!strcmp(pntr[NEXT_LISTNODE]->path, item_path)))
            return i;

        if((!strcmp(pntr[PREV_LISTNODE]->path, item_path)))
            return j;

        pntr[NEXT_LISTNODE] = pntr[NEXT_LISTNODE]->ref[NEXT_LISTNODE];
        pntr[PREV_LISTNODE] = pntr[PREV_LISTNODE]->ref[PREV_LISTNODE];

    }

    return NULL_ITEM(mode);
}

__MSNATIVE_ inline nodelist * __system __export listNo(dim_typ which_item, sel_typ mode)
{
    const dim_typ itemno = getItemsListNo(mode);
    const bool assert = which_item > (itemno*0.5);
    const sel_typ arm_index = 1-(2*assert);
    nodelist * pntr = NULL;

    dim_typ i;

    for(i=assert*(itemno-1), pntr = access(lists)[mode].items_list[assert]; (i*arm_index) < (arm_index*which_item); pntr = pntr->ref[assert], i += arm_index);

    return pntr;
};

__MSNATIVE_ bool __system __export listInsert(const char item_path[static MAX_PATH_LENGTH], sel_typ mode)
{
    nodelist *var, *temp;

    var = (nodelist *) malloc(sizeof(nodelist));
    errMem(var, false);

    strcpy(var->path, item_path);

    if(!access(lists)[mode].items_list[NEXT_LISTNODE])
    {
        access(lists)[mode].items_list[NEXT_LISTNODE]=var;
        access(lists)[mode].items_list[NEXT_LISTNODE]->ref[PREV_LISTNODE]=NULL;
        access(lists)[mode].items_list[NEXT_LISTNODE]->ref[NEXT_LISTNODE]=NULL;
        access(lists)[mode].items_list[PREV_LISTNODE] = access(lists)[mode].items_list[NEXT_LISTNODE];
    }
    else
    {
        access(lists)[mode].items_list[PREV_LISTNODE] = access(lists)[mode].items_list[NEXT_LISTNODE];

        while(access(lists)[mode].items_list[PREV_LISTNODE])
        {
            temp=access(lists)[mode].items_list[PREV_LISTNODE];
            access(lists)[mode].items_list[PREV_LISTNODE]=access(lists)[mode].items_list[PREV_LISTNODE]->ref[NEXT_LISTNODE];
        }

        access(lists)[mode].items_list[PREV_LISTNODE]=var;
        temp->ref[NEXT_LISTNODE]=access(lists)[mode].items_list[PREV_LISTNODE];
        access(lists)[mode].items_list[PREV_LISTNODE]->ref[PREV_LISTNODE]=temp;
        access(lists)[mode].items_list[PREV_LISTNODE]->ref[NEXT_LISTNODE]=NULL;
    }

    ++ access(lists)[mode].itemsno;

    return true;
}

/// WORKS PERFECTLY
/// !!!
__MSNATIVE_ bool __system __export listDelete(dim_typ which_item, sel_typ mode)
{
    nodelist * const link = listNo(which_item, mode);
    nodelist * const prev = link->ref[PREV_LISTNODE];
    nodelist * const next = link->ref[NEXT_LISTNODE];

    if(prev)
        if(next)
        {
            prev->ref[NEXT_LISTNODE] = next;
            next->ref[PREV_LISTNODE] = prev;
        }
        else
        {
            prev->ref[NEXT_LISTNODE] = NULL;
            access(lists)[mode].items_list[PREV_LISTNODE] = prev;
        }
    else
        if (next)
        {
            next->ref[PREV_LISTNODE] = NULL;
            access(lists)[mode].items_list[NEXT_LISTNODE] = next;
        }
        else
            access(lists)[mode].items_list[PREV_LISTNODE] = access(lists)[mode].items_list[NEXT_LISTNODE] = NULL;

    free(link);
    access(lists)[mode].itemsno --;
    return true;
}

// FUNZIONI DI CRUD RELATIVE ALLE VARIE ISTANZE
// DELLE LISTE DI VARI OGGETTI

__MSNATIVE_ dim_typ __system __export itemSelect(sel_typ mode)
{

    char access_point[MIN_STRING];

    switch(mode)
    {
        #ifdef ALLOW_VARLISTMANAGER
        case ENVS:
            strcpy(access_point, change_settings[SETTINGS_ENVSMANAGER].name);
            break;
        #endif
        #ifdef ALLOW_MATMANAGER
        case MATRICES:
            strcpy(access_point, alg_operations[ALGEBRA_MATRICESMANAGER].name);
            break;
        #endif
        #ifdef ALLOW_LOGSMANAGER
        case LOGS:
            strcpy(access_point, change_settings[SETTINGS_LOGSMANAGER].name);
            break;
        #endif
        #ifdef ALLOW_LAYOUTSMANAGER
        case LAYOUTS:
            strcpy(access_point, change_settings[SETTINGS_LAYOUTSMANAGER].name);
            break;
        #endif

        default: return NULL_ITEM(mode);
    }

    PRINTN();

    dim_typ which_item;

    which_item = NULL_ITEM(mode); // getItemsListNo(mode);

    do
    {
        if(isSett(BOOLS_ITEMSSELECTBYPATH))
        {
            char string[MAX_PATH_LENGTH];
            printf2(COLOR_CREDITS, "Enter the desired %s File Path.\n", suite_c.listsnames[mode]);
            printf2(COLOR_CREDITS, "or insert %c to exit SubProgram: %s.\n", SCANFEXIT_CHAR, access_point);
            PRINTL();
            while(scanf("%s", string) != 1 || string[0] == SCANFEXIT_CHAR || (which_item = searchItem(string, mode)) == NULL_ITEM(mode))
            {
                CLEARBUFFER();
                if(string[0] == SCANFEXIT_CHAR) return NULL_ITEM(mode);
                printErr(2, "Inserted Path:\n%s\nrefers to a non-existent %s Item Type", string, suite_c.listsnames[mode]);
            }

            CLEARBUFFER();
        }
        else
        {
            dim_typ i;
            printf("Select the %s Item from the List.\n\n", suite_c.listsnames[mode]);
            // PRINT2N();
            PRINTL();

            if(getItemsListNo(mode) > MAX_CASEINSENSITIVE_CHARS_ALPHABET)
            {
                for(i=0; i<getItemsListNo(mode); ++i)
                    printf("- %hu: %s;\n", i, listNo(i, mode)->path);
                printf("- %hu: Back to %s SubProgram.\n\n", i, access_point);
                PRINTL();


                // dim_typ which_item;
                ityp tmp;

                while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
                        scanf(INPUT_CONVERSION_FORMAT, &tmp) != 1) || tmp != (which_item = (dim_typ)tmp) || tmp < 0 || tmp > getItemsListNo(mode))
                    printErr(1, "Invalid Program Mode");
            }
            else
            {
                for(i=0; i<getItemsListNo(mode); ++i)
                    printf("- %c: %s;\n", i+'A', listNo(i, mode)->path);
                printf("- %c: Back to %s SubProgram.\n\n", i+'A', access_point);
                PRINTL();

                sel_typ tmp;

                do
                    tmp = toupper(getch());
                while(tmp < 'A' || tmp > getItemsListNo(mode)+'A');

                which_item = tmp-'A';
            }

            if(which_item == getItemsListNo(mode)) return NULL_ITEM(mode);


        }

        if(listNo(which_item, mode) == NULL)
            printErr(14, "Failed to select %s Item Type", suite_c.listsnames[mode]);
    }
    while(listNo(which_item, mode) == NULL); // Control CYCLE SELECTION

    return which_item;
}

__MSNATIVE_ inline void __system __export refreshItem(dim_typ which_item, sel_typ mode)
{
    switch(mode)
    {
        case ENVS:
            access(exprVars) = ((exprType *)(listNo(which_item, ENVS)->data));
            break;

        case MATRICES:
            access(curMatrix) = ((matrixObj *)(listNo(which_item, MATRICES)->data));
            /// setCurrentMatrix();
            break;

        case LOGS:
            access(curLog) = ((logObj *)(listNo(which_item, LOGS)->data));
            break;

        case LAYOUTS:
            access(curLayout) = ((layoutObj *)(listNo(which_item, LAYOUTS)->data));
            _changeAlgebraDims(access(curLayout)->algebra);
            break;
    }
    return;
}


__MSNATIVE_ bool __system __export saveItem(dim_typ which_item, sel_typ mode)
{
    if(getItemsListNo(mode) == STARTING_ITEMSNO)
    {
        printErr(1, "There aren't any %s Items currently loaded in memory", suite_c.listsnames[mode]);
        return false;
    }

    if(mode == LAYOUTS)
        getProgramSettings(which_item);
    else
    {
        FILE *fp;

        if((fp = checkForFHErrors(listNo(which_item, mode)->path, "w")) == NULL)
            return false;

        switch(mode)
        {
            case ENVS:
                getVarListEx(fp, which_item);
                break;
            case MATRICES:
            {
                if(!access(curMatrix)->matrix)
                {
                    printErr(11, "CURRENT MATRIX CURRENTLY NOT ALLOCATED");
                    return false;
                }
                matrixObj * const tmp = ((matrixObj *)(listNo(which_item, MATRICES)->data));
                printMatrix(fp, tmp->matrix, tmp->dim); // tmp->dim[RAWS], tmp->dim[COLUMNS]);
                // printfMatrix();
                break;
        	}
            case LOGS:
                logWrite(fp, ((logObj *)(listNo(which_item, LOGS)->data))); /// which_item);
                break;

            default: return false;
        }

        fclose(fp);
    }

    return true;
}

__MSNATIVE_ inline void __system __export passToItem(dim_typ which_item, sel_typ mode, bool savecurrent)
{

    if(isSett(BOOLS_ITEMSAUTOSAVING) && savecurrent)
        saveItem(access(lists)[mode].cur_item, mode);


    // freeExprEvalLists();
    access(lists)[mode].cur_item = which_item;
    updInfo();


    // suite.exprVars.curenv = which_env;
    showCurItem(mode);


    // AT REFRESH ITEM PHASE, the Item will be redirected.
    refreshItem(which_item, mode);

    return;
}

__MSNATIVE_ void __system __export setCurItem(const dim_typ itemID, sel_typ mode)
{
    dim_typ which_item;

    if((which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID) == NULL_ITEM(mode))
        return;

    printf2(COLOR_CREDITS, "\nIt has been selected the %s Item:\n", suite_c.listsnames[mode]);
    printf2(COLOR_CREDITS, listNo(which_item, mode)->path);
    PRINT2N();

    passToItem(which_item, mode, true);
    return;
}

__MSNATIVE_ void __system __export createItemData(dim_typ which_item, sel_typ mode)
{
    switch(mode)
    {
        case ENVS:
            refreshExprEvalVarList(which_item);
            break;

        case MATRICES:
            setCurrentMatrix(which_item);
            break;

        case LOGS:
        {
            logObj * const tmp = malloc(sizeof(logObj));
            errMem(tmp, VSPACE);
            tmp->buffer = malloc(sizeof(char)*DEFAULT_BUFSIZE);
            errMem(tmp->buffer, VSPACE);
            strcpy(tmp->buffer, NULL_CHAR); // initializing log buffer
            tmp->buflen = DEFAULT_BUFSIZE;
            listNo(which_item, LOGS)->data = tmp;
            break;
        }

        case LAYOUTS:
            setProgramSettings(which_item);
            break;

    }

    return;
}

__MSNATIVE_ __WINCALL void __system __export createItem(const char *string, sel_typ mode)
{
    dim_typ which_item;
    char name[MAX_PATH_LENGTH] = NULL_CHAR;
    const bool assert = (__pmode__ == ENVS_OPEN || __pmode__ == MATRICES_OPEN || __pmode__ == LOGS_OPEN);
    // bool mustcreatefile;

    char iname[MINMIN_STRING] = NULL_CHAR;
    strcpy(iname, suite_c.listsnames[mode]);

    // nodelist * item;
    #if WINOS
    // mustcreatefile = false;
    if(isnSett(BOOLS_ITEMSSELECTBYPATH))
    {
        bool wHandler;

        switch(mode)
        {
            case ENVS:
            {
                wHandler = windowsFileHandler(name,
                "MathSuite Vars Lists (*."DEFAULT_VARLIST_FILE_EXTENSION")\0*."DEFAULT_VARLIST_FILE_EXTENSION"\0Text Documents (*.txt)\0*.txt\0All Files (*.*)\0*.*\0",
                DEFAULT_VARLIST_FILE_EXTENSION, assert);
                break;
            }
            case MATRICES:
            {
                wHandler = windowsFileHandler(name,
                "MathSuite Matrices (*."DEFAULT_MATRIX_FILE_EXTENSION")\0*."DEFAULT_MATRIX_FILE_EXTENSION"\0Text Documents (*.txt)\0*.txt\0All Files (*.*)\0*.*\0",
                DEFAULT_MATRIX_FILE_EXTENSION, assert);
                break;
            }
            case LOGS:
            {
                wHandler = windowsFileHandler(name,
                "File LOG (*."DEFAULT_LOG_FILE_EXTENSION")\0*."DEFAULT_LOG_FILE_EXTENSION"\0Text Documents (*.txt)\0*.txt\0All Files (*.*)\0*.*\0",
                DEFAULT_LOG_FILE_EXTENSION, assert);
                break;
            }
            case LAYOUTS:
            {
                wHandler = windowsFileHandler(name,
                "Settings Configuration (*."DEFAULT_LAYOUT_FILE_EXTENSION")\0*."DEFAULT_LAYOUT_FILE_EXTENSION"\0Text Documents (*.txt)\0*.txt\0All Files (*.*)\0*.*\0",
                DEFAULT_LAYOUT_FILE_EXTENSION, assert);
                break;
            }
            default:
                wHandler = false;
                break;
        }

        if(wHandler)
        {
            if(searchItem(name, mode) != NULL_ITEM(mode))
            {
                printErr(17,  "Inserted Path:\n%s\nrefers to an already existent %s Item", name, iname);
                return;
            }
        }
        else
        {
            printErr(14, "Failed to select %s File", iname);
            return;
        }
    }
    else
    #endif
    {
        if(string)
        {
            if((which_item = searchItem(string, mode)) != NULL_ITEM(mode) || (assert && !file_exists(string)))
            {
                if(assert && !file_exists(name))
                    printErr(2, "Inserted Path:\n%s\nrefers to a non-existent File", name);
                if(which_item != NULL_ITEM(mode))
                    printErr(17,  "Inserted Path:\n%s\nrefers to an already existent %s Item", name, iname);
                return;
            }
            strcpy(name, string);
        }
        else
        {
            printf2(COLOR_CREDITS, "Enter the Path of the desired %s Item.\n", iname);
            printf2(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
            PRINTL();

            while(scanf("%s", name) != 1 || name[0] == SCANFEXIT_CHAR || (which_item = searchItem(name, mode)) != NULL_ITEM(mode) || (assert && !file_exists(name)));
            {
                CLEARBUFFER();
                if(name[0] == SCANFEXIT_CHAR) return;
                if(assert && !file_exists(name))
                {
                    printErr(2, "Inserted Path:\n%s\nrefers to a non-existent File", name);
                    return;
                }
                if(which_item != NULL_ITEM(mode))
                {
                    printErr(17,  "Inserted Path:\n%s\nrefers to an already existent %s Item", name, iname);
                    return;
                }
                // mustcreatefile = true;
            }
        }

        CLEARBUFFER();
    }

    const bool zero_starting = getItemsListNo(mode) == STARTING_ITEMSNO;

    if(zero_starting)
        initList(mode);

    if(!listInsert(name, mode))
    {
        printErr(11, "Failed to add:\n%s\nto the %s Items List", name, suite_c.listsnames[mode]);
        return;
    }

    which_item = searchItem(name, mode);


    printf2(COLOR_SYSTEM, "\nID %s: %hu\n", iname, which_item);


    if(__pmode__ == ENVS_SAVE || __pmode__ == MATRICES_SAVE || __pmode__ == LOGS_SAVE || __pmode__ == LAYOUTS_SAVE)
    {
        if(isSett(BOOLS_ITEMSAUTOSAVING))
            saveItem(access(lists)[mode].cur_item, mode);

        access(lists)[mode].cur_item = which_item; // BOOLS_AUTOSETCURITEM IS HARDCODED IN THIS CASE.
        updInfo();


        // suite.exprVars.curenv = which_env;
        createItemData(which_item, mode);
        // suite.exprVars.var_list = ((exprValList *)(listNo(suite.lists[mode].items_list, which_item, mode)->data));
        refreshItem(which_item, mode);
        /// suite.exprVars = ((exprType *)(listNo(suite.lists[ENVS].items_list, which_item, mode)->data));
        updItem(mode);
    }
    else
    {
        if((__pmode__ == ENVS_CREATE || __pmode__ == MATRICES_CREATE || __pmode__ == LOGS_CREATE || __pmode__ == LAYOUTS_CREATE) && !writeFile(name))
            return;

        createItemData(which_item, mode);

        if(isSett(BOOLS_AUTOSETCURITEM))
            passToItem(which_item, mode, !zero_starting);
    }

    printf2(COLOR_CREDITS, "\nIt has been correctly %s the %s Item Type:\n%s\n\n", assert ? "opened" : "saved", iname, name);

    return;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ static void _MS__private __system printItem(const dim_typ itemID, sel_typ mode)
{
    dim_typ which_item;

    if((which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID) == NULL_ITEM(mode))
        return;

        PRINTN();

        switch(mode)
        {
            case ENVS:
                getVarListEx(stdout, which_item);
                break;

            case MATRICES:
            {
                printf2(COLOR_SYSTEM, "\nSelected Matrix:\n\n");
                matrixObj * const tmp = ((matrixObj *)(listNo(which_item, MATRICES)->data));
                printMatrix(stdout, tmp->matrix, tmp->dim);
                break;
        	}

            case LOGS:
                logWrite(stdout, ((logObj *)(listNo(which_item, LOGS)->data)));
                break;

            case LAYOUTS:
                viewProgramSettings(which_item);
                break;
            default: return;
        }

    return;
}

__MSNATIVE_ void __system __export viewItem(const dim_typ itemID, sel_typ mode)
{
    if(access(lists)[mode].itemsno == ITEMS_LASTONESTANDING)
    {
        if(mode == ENVS)
            getVarList(stdout);
        else if(mode == MATRICES)
        {
            if(!access(curMatrix)->matrix)
                printErr(11, "CURRENT MATRIX CURRENTLY NOT ALLOCATED");
            else
                printMatrix(stdout, access(curMatrix)->matrix, access(curMatrix)->dim);
        }
        else if(mode == LOGS)
            logWrite(stdout, ((logObj *)(listNo(access(lists)[LOGS].cur_item, LOGS)->data)));
        else if(mode == LAYOUTS)
            viewProgramSettings(access(lists)[LAYOUTS].cur_item);
    }
    else
        printItem(itemID, mode);


    return;
}

__MSNATIVE_ void __system __export printListItem(const dim_typ itemID, sel_typ mode)
{
    dim_typ which_item;

    if((which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID) == NULL_ITEM(mode))
        return;

    char str[MAX_PATH_LENGTH];

    strcpy(str, listNo(which_item, mode)->path);

    if(printFile(str))
        printf2(COLOR_USER, "It has been correctly printed the %s Item Type:\n%s.\n\n", suite_c.listsnames[mode], str);

    return;
}

__MSNATIVE_ inline void __system __export updItem(sel_typ mode)
{
    if(saveItem(access(lists)[mode].cur_item, mode))
    {
        printf2(COLOR_USER, "Current %s Item has been correctly saved in:\n", suite_c.listsnames[mode]);
        printf2(COLOR_CREDITS, listNo(access(lists)[mode].cur_item, mode)->path);
        PRINTN();
    }
    return;
}

__MSNATIVE_ inline void __system __export updAll(sel_typ mode)
{
    dim_typ i;
    dim_typ j = 0;
    const dim_typ itemsno = getItemsListNo(mode);

    if(getItemsListNo(mode) != STARTING_ITEMSNO)
        for(i=0; i<itemsno; ++i)
            if(saveItem(i, mode))
                ++ j;

    printf2(COLOR_USER, "%hu elements has been correctly saved.\n\n", j);
    return;
}

__MSNATIVE_ void __system __export delItem(const dim_typ itemID, sel_typ mode)
{

    if(mode == LAYOUTS && getItemsListNo(LAYOUTS) == ITEMS_LASTONESTANDING)
    {
        printErr(1, "You cannot delete the Last Layout Item Remaining");
        return;
    }

    dim_typ which_item;

    which_item = MAIN_ITEM;

    const bool assert = __pmode__ == ENVS_DELETEALL || __pmode__ == ENVS_DELETEALLPHYSICALS ||
                        __pmode__ == MATRICES_DELETEALL  || __pmode__ == MATRICES_DELETEALLPHYSICALS  ||
                        __pmode__ == LOGS_DELETEALL || __pmode__ == LOGS_DELETEALLPHYSICALS  ;

    if(assert && getItemsListNo(mode) == STARTING_ITEMSNO)
    {
        printErr(1, "Actually there are no Items of Type: %s to delete", suite_c.listsnames[mode]);
        return;
    }

    if((!assert) && (which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID) == NULL_ITEM(mode))
        return;

    char name[MAX_PATH_LENGTH];
    strcpy(name, (listNo(which_item, mode)->path));

    if(which_item == access(lists)[mode].cur_item)
    {
        if(isSett(BOOLS_ITEMSAUTOSAVING))
            saveItem(which_item, mode);

        const bool one_last = getItemsListNo(mode) != ITEMS_LASTONESTANDING;

        if(!listDeleteProc(name, mode))
            return;

        if(isSett(BOOLS_AUTOSETCURITEM) && one_last)
        {
            access(lists)[mode].cur_item -= which_item == getItemsListNo(mode);
            refreshItem(access(lists)[mode].cur_item, mode);
            // showCurItem(mode);
        }

    }
    else
    {
        if(!listDeleteProc(name, mode))
            return;

        if(isSett(BOOLS_AUTOSETCURITEM) && which_item < access(lists)[mode].cur_item)
            access(lists)[mode].cur_item --;
            // showCurItem(mode);
        // suite.exprVars.curenv --;
    }

    showCurItem(mode);

    if(__pmode__ == ENVS_DELETEPHYSICAL || __pmode__ == ENVS_DELETEALLPHYSICALS || __pmode__ == MATRICES_DELETEPHYSICAL   ||
       __pmode__ == MATRICES_DELETEALLPHYSICALS || __pmode__ == LOGS_DELETEPHYSICAL || __pmode__ == LOGS_DELETEALLPHYSICALS||
       __pmode__ == LAYOUTS_DELETEPHYSICAL || __pmode__ == LAYOUTS_DELETEALLPHYSICALS)
    {
        int err;
        if((err = remove(name)))
            printErr(2, "An error during:\n%s\n%s File deleting process might have occurred", name, suite_c.listsnames[mode]);
    }

    if(assert && getItemsListNo(mode) != STARTING_ITEMSNO)
        delItem(NULL_ITEM(mode), mode);

    return;
}

// CMDNAME Manager getItemID DEBUGGER/DISPATCHER
__MSNATIVE_ dim_typ __system __export getItemID(char *string, sprog * const prog, sel_typ mode)
{
    ityp tmp;
    dim_typ which_item;

    tmp = 0.00;

    if(PARSING_SYSTEM_ALLOWED)
    {
        if((!parse(string, &tmp)) || tmp != (which_item = (dim_typ)tmp) || which_item < 0 || which_item > getItemsListNo(mode)-1)
        {
            printUsage(prog);
            return NULL_ITEM(mode);
        }
    }
    else if((tmp = strtoul(string, NULL, sizeof(which_item))) != (which_item = (dim_typ)tmp) || which_item < 0 || which_item > getItemsListNo(mode)-1)
    {
        printUsage(prog);
        return NULL_ITEM(mode);
    }

    return which_item;
}

__MSNATIVE_ void __system __export relItem(const dim_typ itemID, bool mode)
{
    dim_typ which_item;

    if((which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID) == NULL_ITEM(mode))
        return;

    printf2(COLOR_USER, "\nIt has been correctly reloaded the %s Item Type:\n", suite_c.listsnames[mode]);

    char name[MAX_PATH_LENGTH];
    strcpy(name, listNo(which_item, mode)->path);

    if(mode == LOGS)
        _flushLogBuf(((logObj *)(listNo(which_item, LOGS)->data)));

    printf2(COLOR_CREDITS, name);
    PRINT2N();

    writeFile(name);

    return;
}

__MSNATIVE_ void __system __export renItem(const char *string, const dim_typ itemID, sel_typ mode)
{
    dim_typ which_item;
    bool assert;

    char iname[MINMIN_STRING];
    strcpy(iname, suite_c.listsnames[mode]);

    size_t len;
    char name[MAX_PATH_LENGTH];
    char newname[MAX_PATH_LENGTH];

    len = 0;
    assert = false;

    if(string && itemID != NULL_ITEM(mode))
    {
        which_item = itemID;
        strcpy(name, listNo(which_item, mode)->path);
        if(strrchr(string, SCANFEXIT_CHAR) != NULL || (len = strlen(string)) > MAX_PATH_LENGTH || len < 2)
        strcpy(newname, string);
    }
    else
        do
        {

            CLEARBUFFER();

            if((which_item = itemSelect(mode)) == NULL_ITEM(mode))
                return;

            strcpy(name, (listNo(which_item, mode)->path));

            printf2(COLOR_CREDITS, "\nEnter the newname you wish to give to the selected %s Item.\n", iname);
            printf2(COLOR_CREDITS, "or insert "SCANFEXIT_STRING" to go Back.\n\n");

            CLEARBUFFER();


            while(scanf("%s", newname) != 1 || (assert = (newname[0] == SCANFEXIT_CHAR)) ||
                strrchr(newname, SCANFEXIT_CHAR) != NULL || (len = strlen(newname)) > MAX_PATH_LENGTH || len < 2)
            {
                CLEARBUFFER();
                if(assert) break;
                printErr(20+(18*(len>MAX_PATH_LENGTH)), "Invalid inserted Name: %s", newname);
            }

        }
        while(assert);

    CLEARBUFFER();

    switch(mode)
    {
        case ENVS:
            strcat(newname, EXTENSION_DOT NULL_CHAR DEFAULT_VARLIST_FILE_EXTENSION);
            break;
        case MATRICES:
            strcat(newname, EXTENSION_DOT NULL_CHAR DEFAULT_MATRIX_FILE_EXTENSION);
            break;
        case LOGS:
            strcat(newname, EXTENSION_DOT NULL_CHAR DEFAULT_LOG_FILE_EXTENSION);
            break;
        case LAYOUTS:
            strcat(newname, EXTENSION_DOT NULL_CHAR DEFAULT_LAYOUT_FILE_EXTENSION);
            break;
    }

    if(frename(name, newname))
    {
        strcpy(listNo(which_item, mode)->path, newname);
        printf2(COLOR_USER, "%s Item Type has been correctly renamed:", iname);
        printf2(COLOR_CREDITS, "\n%s", name);
        printf2(COLOR_USER, "\nin:\n");
        printf2(COLOR_CREDITS, newname);
        printf2(COLOR_USER, ".\n\n");
    }

    return;
}

#endif
