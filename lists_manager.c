// lists_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#if((!defined __DISABLE_VARLISTMANAGER) || (!defined __DISABLE_MATMANAGER) || (!defined __DISABLE_LOGSMANAGER) || (!defined __DISABLE_LAYOUTSMANAGER))

// FUNCTIONS DECLARATIONS

__MSSHELL_WRAPPER_ __MATHSUITE static void   printItem(const dim_typ, sel_typ);
__MSSHELL_WRAPPER_ __MATHSUITE static void   showCurItem(sel_typ);
__MSSHELL_WRAPPER_ __MATHSUITE static bool   listDeleteProc(const char [static MAX_PATH_LENGTH], sel_typ);

// FUNZIONI LISTE

__MSSHELL_WRAPPER_ __MATHSUITE static void   showCurItem(sel_typ mode)
{
    char str[MINMIN_STRING] = NULL_CHAR;
    strcpy(str, suite_c.listsnames[mode]);
    toupper_s(str);
    msprintf(COLOR_CREDITS, "CURRENT %s:\n%s\n", str, listNo(access(lists)[mode].cur_item, mode)->path);
    return;
}

__MATHSUITE inline sel_typ   checkItemTypeByExtension(const char extension[static MAX_EXTENSION_LENGTH])
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

__MSSHELL_WRAPPER_ __MATHSUITE bool    listInsertProc(const char item_path[static MAX_PATH_LENGTH], sel_typ mode)
{
    const bool assert = getItemsListNo(mode) == STARTING_ITEMSNO;
	
    if(assert)
        initList(mode);

    if(!file_exists(item_path))
    {
        printErr(2, "An error during:\n%s\n%s File opening process might have occurred", item_path, suite_c.listsnames[mode]);
        return false;
    }
    
    if(searchItem(item_path, mode) != NULL_ITEM(mode))
    {
		printErr(17,  "Inserted Path:\n%s\nrefers to an already existing %s Item", item_path, suite_c.listsnames[mode]);    
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

__MSSHELL_WRAPPER_ __MATHSUITE static bool   listDeleteProc(const char item_path[static MAX_PATH_LENGTH], sel_typ mode)
{
    // this will free the envitem SPECIFIED VARLIST

    const dim_typ which_item = searchItem(item_path, mode);
    void * data = listNo(which_item, mode)->data;

	switch(mode)
	{
		case ENVS:
			(void) strtok(listNo(which_item, mode)->path, EXTENSION_DOT);
			deleteVarList(listNo(which_item, mode)->path);
			exprValListFree(((exprType *)data)->var_list);
			break; 
		case MATRICES:
		{
			matrixObj * const tmp = ((matrixObj*)data);
		    if(tmp->matrix)
		        matrixFree(&(tmp->matrix), tmp->dim);
		 (void) strtok(listNo(which_item, mode)->path, EXTENSION_DOT);
		    deleteMat(listNo(which_item, mode)->path, DELETEMAT_MAT);
		    break;
		}
		case LOGS:
		    free(((logObj *)data)->buffer);
		    break;
		case LAYOUTS:
			// do nothing
			break;
	}

    free(data);

    if(!listDelete(which_item, mode))
    {
        printErr(11, "Failed to delete:\n%s\n%s File", item_path, suite_c.listsnames[mode]);
        return false;
    }

    msprintf(COLOR_USER, "\nIt has been correctly deleted the %s File:", suite_c.listsnames[mode]);
    msprintf(COLOR_CREDITS, "\n%s\n\n", item_path);
    return true;
}

__MATHSUITE dim_typ   searchItem(const char item_path[static MAX_PATH_LENGTH], sel_typ mode)
{
    dim_typ i, j;
    nodelist * pntr[2];

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

__MATHSUITE inline nodelist *   listNo(dim_typ which_item, sel_typ mode)
{
	dim_typ i;
    const dim_typ itemno = getItemsListNo(mode);
    const bool assert = which_item > (itemno*0.5);
    const sel_typ arm_index = 1-(assert<<1);
    nodelist * pntr = NULL;

    for(i=assert*(itemno-1), pntr = access(lists)[mode].items_list[assert]; (i*arm_index) < (arm_index*which_item); pntr = pntr->ref[assert], i += arm_index);
    return pntr;
};

__MATHSUITE bool   listInsert(const char item_path[static MAX_PATH_LENGTH], sel_typ mode)
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
__MATHSUITE bool   listDelete(dim_typ which_item, sel_typ mode)
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

__MATHSUITE dim_typ   itemSelect(sel_typ mode)
{

    //char access_point[MIN_STRING];
    
    static const char *access_points_names[INFO_STRING] =
    {
    	#ifndef __DISABLE_VARLISTMANAGER
			change_settings[SETTINGS_ENVSMANAGER].name,
		#endif
		#ifndef __DISABLE_MATMANAGER
			alg_operations[ALGOPS_MATRICESMANAGER].name,
		#endif
		#ifndef __DISABLE_LOGSMANAGER
			change_settings[SETTINGS_LOGSMANAGER].name,
		#endif
		#ifndef __DISABLE_LAYOUTSMANAGER
			change_settings[SETTINGS_LAYOUTSMANAGER].name
		#endif
    };
    
	dim_typ which_item = NULL_ITEM(mode); // getItemsListNo(mode);
    const char * const access_point = access_points_names[mode];
    
    PRINTN();

    do
    {
        if(isSett(BOOLS_ITEMSSELECTBYPATH))
        {
            char string[MAX_PATH_LENGTH];
            msprintf(COLOR_CREDITS, "Enter the desired %s File Path.\n", suite_c.listsnames[mode]);
            msprintf(COLOR_CREDITS, "or insert %c to exit SubProgram: %s.\n", SCANFEXIT_CHAR, access_point);
            PRINTL();
            while(scanf("%s", string) != 1 || string[0] == SCANFEXIT_CHAR || (which_item = searchItem(string, mode)) == NULL_ITEM(mode))
            {
                CLEARBUFFER();
                if(string[0] == SCANFEXIT_CHAR) return NULL_ITEM(mode);
                printErr(2, "Inserted Path:\n%s\nrefers to a non-existing %s Item Type", string, suite_c.listsnames[mode]);
            }

            CLEARBUFFER();
        }
        else
        {
            dim_typ i;
            printf("Select the %s Item from the List.\n\n", suite_c.listsnames[mode]);
            PRINTL();

            if(getItemsListNo(mode) > MAX_CASEINSENSITIVE_CHARS_ALPHABET)
            {
                for(i=0; i<getItemsListNo(mode); ++i)
                    printf("- %hu: %s;\n", i, listNo(i, mode)->path);
                printf("- %hu: Back to %s SubProgram.\n\n", i, access_point);
                PRINTL();


                mpfr_t tmp;
                while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (which_item = mpfr_get_ui(tmp, MPFR_RNDN))) || which_item < 0 || which_item > getItemsListNo(mode))
                {
                	mpfr_clear(tmp);
					printErr(1, "Invalid Program Mode");
				}
                mpfr_clear(tmp); 
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

__MATHSUITE inline void   refreshItem(dim_typ which_item, sel_typ mode)
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
            break;
    }
    return;
}


__MATHSUITE bool   saveItem(dim_typ which_item, sel_typ mode)
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
				#ifndef __DISABLE_DATABASE
					if(!access(curLayout)->database.con)
	        			dbEstablishConnection(0, NULL);
				
					if(access(curLayout)->database.con)
						getVarListEx(NULL, which_item);
				#endif
				break;
			case MATRICES:
			{
				if(!access(curMatrix)->matrix)
	            {
	                printErr(11, "CURRENT MATRIX CURRENTLY NOT ALLOCATED");
	                return false;
	            }
	            matrixObj * const tmp = ((matrixObj *)(listNo(which_item, MATRICES)->data));
	            printMatrix(fp, &(tmp->matrix), tmp->dim); // tmp->dim[ROWS], tmp->dim[COLUMNS]);
	            
	            
	            #ifndef __DISABLE_DATABASE
		            if(!access(curLayout)->database.con)
	        			dbEstablishConnection(0, NULL);
				
					if(access(curLayout)->database.con)
					{
						char matname[MAX_PATH_LENGTH];
		        		strcpy(matname, listNo(which_item, MATRICES)->path);
		        		(void) strtok(matname, EXTENSION_DOT);
						deleteMat(matname, DELETEMAT_ROW);
						_printMatrix(NULL, &(tmp->matrix), tmp->dim, matname); // which_item);
					}
				#endif
				
				break;
			}
			case LOGS:
				logWrite(fp, ((logObj *)(listNo(which_item, LOGS)->data))); /// which_item);
				break;
		}
        fclose(fp);
    }

    return true;
}

__MATHSUITE inline void   passToItem(dim_typ which_item, sel_typ mode, const bool options) // bool savecurrent)
{
    if(isSett(BOOLS_ITEMSAUTOSAVING) && (options & PTI_SAVECURRENT)) // savecurrent)
        saveItem(access(lists)[mode].cur_item, mode);
        
    accessCurrentSession()->lastItem[mode] = access(lists)[mode].cur_item = which_item;
    
    if(options & PTI_UPDINFO)
    	updInfo(NULL);  
    	
    if(options & PTI_SHOWITEM)
    	showCurItem(mode);
    	
    refreshItem(which_item, mode);
    return;
}

__MATHSUITE void   setCurItem(const dim_typ itemID, sel_typ mode)
{
    const dim_typ which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID;

    if(which_item == NULL_ITEM(mode))
        return;

    msprintf(COLOR_CREDITS, "\nIt has been selected the %s Item:\n", suite_c.listsnames[mode]);
    msprintf(COLOR_CREDITS, listNo(which_item, mode)->path);
    PRINT2N();
    passToItem(which_item, mode, PTI_SAVECURRENT | PTI_SHOWITEM | PTI_UPDINFO);
    return;
}

__MATHSUITE void   createItemData(dim_typ which_item, sel_typ mode)
{
	static void (* const func_data[MAX_LISTS])(dim_typ) =
	{
		refreshExprEvalVarList,
		setCurrentMatrix,
		setCurrentLog,
		setProgramSettings
	};
	
	func_data[mode](which_item);
    return;
}

__MATHSUITE __WINCALL void   createItem(const char *string, sel_typ mode)
{
    dim_typ which_item;
    char name[MAX_PATH_LENGTH] = NULL_CHAR;
    const bool assert = (__pmode__ == ENVS_OPEN || __pmode__ == MATRICES_OPEN || __pmode__ == LOGS_OPEN);
   	const char *const iname = suite_c.listsnames[mode];

    // nodelist * item;
    #ifdef WINOS
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
                printErr(17,  "Inserted Path:\n%s\nrefers to an already existing %s Item", name, iname);
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
                    printErr(2, "Inserted Path:\n%s\nrefers to a non-existing File", name);
                if(which_item != NULL_ITEM(mode))
                    printErr(17,  "Inserted Path:\n%s\nrefers to an already existing %s Item", name, iname);
                return;
            }
            strcpy(name, string);
        }
        else
        {
            msprintf(COLOR_CREDITS, "Enter the Path of the desired %s Item.\n", iname);
            msprintf(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
            PRINTL();

            while(scanf("%s", name) != 1 || name[0] == SCANFEXIT_CHAR || (which_item = searchItem(name, mode)) != NULL_ITEM(mode) || (assert && !file_exists(name)))
            {
                CLEARBUFFER();
                if(name[0] == SCANFEXIT_CHAR) return;
                if(assert && !file_exists(name))
                {
                    printErr(2, "Inserted Path:\n%s\nrefers to a non-existing File", name);
                    return;
                }
                if(which_item != NULL_ITEM(mode))
                {
                    printErr(17,  "Inserted Path:\n%s\nrefers to an already existing %s Item", name, iname);
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


    if(__pmode__ == ENVS_SAVE || __pmode__ == MATRICES_SAVE || __pmode__ == LOGS_SAVE || __pmode__ == LAYOUTS_SAVE)
    {
        if(isSett(BOOLS_ITEMSAUTOSAVING))
            saveItem(access(lists)[mode].cur_item, mode);

        access(lists)[mode].cur_item = which_item; // BOOLS_AUTOSETCURITEM IS HARDCODED IN THIS CASE.
        createItemData(which_item, mode);
        refreshItem(which_item, mode);
        updItem(mode);
    }
    else
    {
        if((__pmode__ == ENVS_CREATE || __pmode__ == MATRICES_CREATE || __pmode__ == LOGS_CREATE || __pmode__ == LAYOUTS_CREATE) && !writeFile(name))
            return;

        createItemData(which_item, mode);

        if(isSett(BOOLS_AUTOSETCURITEM))
            passToItem(which_item, mode, (PTI_SAVECURRENT*(!zero_starting)) | PTI_SHOWITEM | PTI_UPDINFO);
    }

    updInfo(NULL);  
    msprintf(COLOR_CREDITS, "\nIt has been correctly %s the %s Item Type:\n%s\n\n", assert ? "opened" : "saved", iname, name);
    return;
}

__MSSHELL_WRAPPER_ __MATHSUITE static void   printItem(const dim_typ itemID, sel_typ mode)
{
    const dim_typ which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID;

    if(which_item == NULL_ITEM(mode))
        return;

        PRINTN();

        switch(mode)
        {
            case ENVS:
                getVarListEx(stdout, which_item);
                break;

            case MATRICES:
            {
                msprintf(COLOR_SYSTEM, "\nSelected Matrix:\n\n");
                matrixObj * const tmp = ((matrixObj *)(listNo(which_item, MATRICES)->data));
                printMatrix(stdout, &(tmp->matrix), tmp->dim);
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

__MATHSUITE void   viewItem(const dim_typ itemID, sel_typ mode)
{
    if(access(lists)[mode].itemsno == ITEMS_LASTONESTANDING)
    {
    	switch(mode)
    	{
    		case ENVS:
    			getVarList(stdout);
    			break;
    		case MATRICES:
    			if(!access(curMatrix)->matrix)
                	printErr(11, "CURRENT MATRIX CURRENTLY NOT ALLOCATED");
	            else
	                printMatrix(stdout, &(access(curMatrix)->matrix), access(curMatrix)->dim);
	            break;
	        case LOGS:
	            logWrite(stdout, ((logObj *)(listNo(access(lists)[LOGS].cur_item, LOGS)->data)));
	            break;
	        case LAYOUTS:
	        	viewProgramSettings(access(lists)[LAYOUTS].cur_item);
	        	break;
    	} 
    }
    else
        printItem(itemID, mode);


    return;
}

__MATHSUITE void   printListItem(const dim_typ itemID, sel_typ mode)
{
    dim_typ which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID;

    if(which_item == NULL_ITEM(mode))
        return;

    const char * const str = listNo(which_item, mode)->path;

    if(printFile(str))
        msprintf(COLOR_USER, "It has been correctly printed the %s Item Type:\n%s.\n\n", suite_c.listsnames[mode], str);

    return;
}

__MATHSUITE inline void   updItem(sel_typ mode)
{
    if(saveItem(access(lists)[mode].cur_item, mode))
    {
        msprintf(COLOR_USER, "Current %s Item has been correctly saved in:\n", suite_c.listsnames[mode]);
        msprintf(COLOR_CREDITS, listNo(access(lists)[mode].cur_item, mode)->path);
        PRINTN();
    }
    return;
}

__MATHSUITE inline void   updAll(sel_typ mode)
{
    dim_typ i, j = 0;
    const dim_typ itemsno = getItemsListNo(mode);

    if(getItemsListNo(mode) != STARTING_ITEMSNO)
        for(i=0; i<itemsno; ++i)
            if(saveItem(i, mode))
                ++ j;

	const register int backupSocket = access(lastSessionSocket);
	access(lastSessionSocket) = INVALID_SOCKET;
    msprintf(COLOR_USER, "%hu elements has been correctly saved.\n\n", j);
    access(lastSessionSocket) = backupSocket;
    return;
}

__MATHSUITE void   delItem(const dim_typ itemID, sel_typ mode)
{

    if(mode == LAYOUTS && getItemsListNo(LAYOUTS) == ITEMS_LASTONESTANDING)
    {
        printErr(1, "You cannot delete the Last Layout Item Remaining");
        return;
    }

    dim_typ which_item = MAIN_ITEM;
    const bool assert = __pmode__ == ENVS_DELETEALL || __pmode__ == ENVS_DELETEALLPHYSICALS ||
                        __pmode__ == MATRICES_DELETEALL  || __pmode__ == MATRICES_DELETEALLPHYSICALS ||
                        __pmode__ == LOGS_DELETEALL || __pmode__ == LOGS_DELETEALLPHYSICALS  ;

    if(assert && getItemsListNo(mode) == STARTING_ITEMSNO)
    {
        printErr(1, "Actually there are no Items of Type: %s to delete", suite_c.listsnames[mode]);
        return;
    }

    if((!assert) && (which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID) == NULL_ITEM(mode))
        return;

    const char * const name = listNo(which_item, mode)->path;

    if(which_item == access(lists)[mode].cur_item)
    {
        if(isSett(BOOLS_ITEMSAUTOSAVING))
            saveItem(which_item, mode);

        const bool one_last = getItemsListNo(mode) != ITEMS_LASTONESTANDING;
        
        if(__pmode__ == ENVS_DELETEPHYSICAL || __pmode__ == ENVS_DELETEALLPHYSICALS || __pmode__ == MATRICES_DELETEPHYSICAL   ||
       __pmode__ == MATRICES_DELETEALLPHYSICALS || __pmode__ == LAYOUTS_DELETEPHYSICAL || __pmode__ == LAYOUTS_DELETEALLPHYSICALS)
		{
		    int err;
		    
		    if((err = remove(name)))
		        printErr(2, "An error during:\n%s\n%s File deleting process might have occurred", name, suite_c.listsnames[mode]);
		}


        if(!listDeleteProc(name, mode))
            return;

        if(isSett(BOOLS_AUTOSETCURITEM) && one_last)
        {
            access(lists)[mode].cur_item -= which_item == getItemsListNo(mode);
            refreshItem(access(lists)[mode].cur_item, mode);
        }

    }
    else
    {
    	if(__pmode__ == ENVS_DELETEPHYSICAL || __pmode__ == ENVS_DELETEALLPHYSICALS || __pmode__ == MATRICES_DELETEPHYSICAL   ||
       __pmode__ == MATRICES_DELETEALLPHYSICALS || __pmode__ == LAYOUTS_DELETEPHYSICAL || __pmode__ == LAYOUTS_DELETEALLPHYSICALS)
	    {
	        int err;
	        
	        if((err = remove(name)))
	            printErr(2, "An error during:\n%s\n%s File deleting process might have occurred", name, suite_c.listsnames[mode]);
	    }


        if(!listDeleteProc(name, mode))
            return;

        if(isSett(BOOLS_AUTOSETCURITEM) && which_item < access(lists)[mode].cur_item)
            access(lists)[mode].cur_item --;
    }

    showCurItem(mode);
    
    if(assert && getItemsListNo(mode) != STARTING_ITEMSNO)
        delItem(NULL_ITEM(mode), mode);

    return;
}

// CMDNAME Manager getItemID DEBUGGER/DISPATCHER
__MATHSUITE dim_typ   getItemID(char *string, sprog * const prog, sel_typ mode)
{
    mpfr_t tmp; 
    dim_typ which_item;
    
    if((!parse(string, &tmp)) || mpfr_cmp_ui(tmp, (which_item = mpfr_get_ui(tmp, MPFR_RNDN))) || which_item < 0 || which_item > getItemsListNo(mode)-1)
    {
    	mpfr_clear(tmp);
        printUsage(prog);
        return NULL_ITEM(mode);
    }
    
   	mpfr_clear(tmp);
    return which_item;
}

__MATHSUITE void   relItem(const dim_typ itemID, bool mode)
{
    const dim_typ which_item = itemID == NULL_ITEM(mode) ? itemSelect(mode) : itemID;

    if(which_item == NULL_ITEM(mode))
        return;

    msprintf(COLOR_USER, "\nIt has been correctly reloaded the %s Item Type:\n", suite_c.listsnames[mode]);

    const char * const name = listNo(which_item, mode)->path;
    
    if(mode == LOGS)
        _flushLogBuf(((logObj *)(listNo(which_item, LOGS)->data)));

    msprintf(COLOR_CREDITS, name);
    PRINT2N();
    writeFile(name);
    return;
}

__MATHSUITE void   renItem(const char *string, const dim_typ itemID, sel_typ mode)
{
	size_t len = 0;
    dim_typ which_item;
    bool assert = false;
    char newname[MAX_PATH_LENGTH];
    const char * const iname = suite_c.listsnames[mode];
    char * name;
    
    if(string && itemID != NULL_ITEM(mode))
    {
        which_item = itemID;
        name = listNo(which_item, mode)->path;
        if(strrchr(string, SCANFEXIT_CHAR) != NULL || (len = strlen(string)) > MAX_PATH_LENGTH || len < 2)
        strcpy(newname, string);
    }
    else
        do
        {

            CLEARBUFFER();

            if((which_item = itemSelect(mode)) == NULL_ITEM(mode))
                return;

            name = listNo(which_item, mode)->path;
            msprintf(COLOR_CREDITS, "\nEnter the newname you wish to give to the selected %s Item.\n", iname);
            msprintf(COLOR_CREDITS, "or insert "SCANFEXIT_STRING" to go Back.\n\n");
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
    
    static const char init_string[MAX_LISTS][MAX_EXTENSION_LENGTH] =
    {
    	DEFAULT_VARLIST_FILE_EXTENSION,
    	DEFAULT_MATRIX_FILE_EXTENSION,
    	DEFAULT_LOG_FILE_EXTENSION,
    	DEFAULT_LAYOUT_FILE_EXTENSION
    };
    
    strcat(newname, EXTENSION_DOT NULL_CHAR);
    strcat(newname, init_string[mode]);

    if(frename(name, newname))
    {
        strcpy(listNo(which_item, mode)->path, newname);
        msprintf(COLOR_USER, "%s Item Type has been correctly renamed:", iname);
        msprintf(COLOR_CREDITS, "\n%s", name);
        msprintf(COLOR_USER, "\nin:\n");
        msprintf(COLOR_CREDITS, newname);
        msprintf(COLOR_USER, ".\n\n");
    }

    return;
}

#endif
