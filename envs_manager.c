// envs_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifndef __DISABLE_VARLISTMANAGER

// FUNCTIONS DECLARATIONS


__MSSHELL_WRAPPER_ static void  __lmp_prog setCurEnv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void  __lmp_prog createEnv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog viewEnv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog printEnv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog updEnv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog updAllEnvs(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog delEnv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog relEnv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog syncEnv(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog renEnv(const sel_typ argc, char ** argv);
#ifndef __DISABLE_DATABASE
	__MSSHELL_WRAPPER_ static void  __lmp_prog persistEnv(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void  __lmp_prog retrieveEnv(const sel_typ argc, char ** argv);
#endif

sprog envs_manager[MAX_ENVSMANAGER_PROGS] =
{
    [ENVS_SETCURRENT] =
    {
    	CMD_SETCURENV,
        NAME_SETCURENV,
        USAGE_SETCURENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_SETCURENV,
        #endif
        setCurEnv,
        ARGC_SETCURENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_OPEN] =
    {
    	CMD_OPENENV,
        NAME_OPENENV,
        USAGE_OPENENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_OPENENV,
        #endif
        createEnv,
        ARGC_OPENENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_CREATE] =
    {
    	CMD_CREATEENV,
        NAME_CREATEENV,
        USAGE_CREATEENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_CREATEENV,
        #endif
        createEnv,
        ARGC_CREATEENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_READ] =
    {
    	CMD_READENV,
        NAME_READENV,
        USAGE_READENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_READENV,
        #endif
        viewEnv,
        ARGC_READENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_PRINT] =
    {
    	CMD_PRTENV,
        NAME_PRTENV,
        USAGE_PRTENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_PRTENV,
        #endif
        printEnv,
        ARGC_PRTENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_UPDATE] =
    {
    	CMD_UPDENV,
        NAME_UPDENV,
        USAGE_UPDENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDENV,
        #endif
        updEnv,
        ARGC_UPDENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_UPDATEALL] =
    {
    	CMD_UPDALLENVS,
        NAME_UPDALLENVS,
        USAGE_UPDALLENVS,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDALLENVS,
        #endif
        updAllEnvs,
        ARGC_UPDALLENVS,
        AUTOMATIC,
        CHILD
    },
    [ENVS_SAVE] =
    {
    	CMD_SAVEENV,
        NAME_SAVEENV,
        USAGE_SAVEENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_SAVEENV,
        #endif
        createEnv,
        ARGC_SAVEENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_DELETE] =
    {
    	CMD_DELENV,
        NAME_DELENV,
        USAGE_DELENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELENV,
        #endif
        delEnv,
        ARGC_DELENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_DELETEPHYSICAL] =
    {
    	CMD_DELENV2,
        NAME_DELENV2,
        USAGE_DELENV2,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELENV2,
        #endif
        delEnv,
        ARGC_DELENV2,
        AUTOMATIC,
        CHILD
    },
    [ENVS_DELETEALL] =
    {
    	CMD_DELALLENVS,
        NAME_DELALLENVS,
        USAGE_DELALLENVS,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELALLENVS,
        #endif
        delEnv,
        ARGC_DELALLENVS,
        AUTOMATIC,
        CHILD
    },
    [ENVS_DELETEALLPHYSICALS] =
    {
    	CMD_DELALLENVS2,
        NAME_DELALLENVS2,
        USAGE_DELALLENVS2,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELALLENVS2,
        #endif
        delEnv,
        ARGC_DELALLENVS2,
        AUTOMATIC,
        CHILD
    },
    [ENVS_RELOAD] =
    {
    	CMD_RESETENV,
        NAME_RESETENV,
        USAGE_RESETENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_RESETENV,
        #endif
        relEnv,
        ARGC_RESETENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_SYNC] =
    {
    	CMD_SYNCENV,
        NAME_SYNCENV,
        USAGE_SYNCENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_SYNCENV,
        #endif
        syncEnv,
        ARGC_SYNCENV,
        AUTOMATIC,
        CHILD
    },
    [ENVS_RENAME] =
    {
    	CMD_RENENV,
        NAME_RENENV,
        USAGE_RENENV,
        #ifndef __DISABLE_DATABASE
        LEVEL_RENENV,
        #endif
        renEnv,
        ARGC_RENENV,
        AUTOMATIC,
        CHILD
    }
    #ifndef __DISABLE_DATABASE
    , [ENVS_PERSIST] =
    {
    	CMD_PERSISTENV,
    	NAME_PERSISTENV,
    	USAGE_PERSISTENV,
    	LEVEL_PERSISTENV,
    	persistEnv,
    	ARGC_PERSISTENV,
    	AUTOMATIC,
    	CHILD 
    },
    [ENVS_RETRIEVE] =
	{
		CMD_RETRIEVEENV,
		NAME_RETRIEVEENV,
		USAGE_RETRIEVEENV,
		LEVEL_RETRIEVEENV,
		retrieveEnv,
		ARGC_RETRIEVEENV,
		AUTOMATIC,
		CHILD
	}
	#endif
};

// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static inline void  __lmp_prog setCurEnv(const sel_typ argc, char ** argv)
{
    setCurItem(argc ? getItemID(argv[0], &envs_manager[ENVS_SETCURRENT], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static inline void  __lmp_prog createEnv(const sel_typ argc, char ** argv)
{
    createItem(argc ? argv[0] : NULL, ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog viewEnv(const sel_typ argc, char ** argv)
{
    viewItem(argc ? getItemID(argv[0], &envs_manager[ENVS_READ], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog printEnv(const sel_typ argc, char ** argv)
{
    printListItem(argc ? getItemID(argv[0], &envs_manager[ENVS_PRINT], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog updEnv(const sel_typ argc, char ** argv)
{
    updItem(ENVS);
    return;
}
__MSSHELL_WRAPPER_ static inline void  __lmp_prog updAllEnvs(const sel_typ argc, char ** argv)
{
    updAll(ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog delEnv(const sel_typ argc, char ** argv)
{
    delItem(argc ? getItemID(argv[0], &envs_manager[__pmode__], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog relEnv(const sel_typ argc, char ** argv)
{
    relItem(argc ? getItemID(argv[0], &envs_manager[ENVS_RELOAD], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog syncEnv(const sel_typ argc, char ** argv)
{
    dim_typ which_env;

    if((which_env = argc ? getItemID(argv[0], &envs_manager[ENVS_SYNC], ENVS) : itemSelect(ENVS)) == NULL_ENV)
       return;

	FILE *fp;
    char name[MAX_PATH_LENGTH];
    strcpy(name, listNo(which_env, ENVS)->path);

    if((fp = checkForFHErrors(name, "a")) == NULL)
        return;

    getVarList(fp);
    fclose(fp);
    msprintf(COLOR_USER, "\nThe Env:\n");
    msprintf(COLOR_CREDITS, name);
    msprintf(COLOR_USER, "\nhas been correctly synchronized with the Env:\n");
    msprintf(COLOR_CREDITS, listNo(access(lists)[ENVS].cur_item, ENVS)->path);
    PRINT2N();

    if(isSett(BOOLS_AUTOSETCURITEM))
        passToItem(which_env, ENVS, PTI_SAVECURRENT | PTI_SHOWITEM | PTI_UPDINFO);

    return;
}

#ifndef __DISABLE_DATABASE
	__MSSHELL_WRAPPER_ static inline void  __lmp_prog renEnv(const sel_typ argc, char ** argv)
	{
	    renItem(argc > 1 ? argv[1] : NULL, argc ? getItemID(argv[0], &envs_manager[ENVS_RENAME], ENVS) : getItemsListNo(ENVS), ENVS);
	    return;
	}
	
	__MSSHELL_WRAPPER_ static inline void  __lmp_prog persistEnv(const sel_typ argc, char ** argv)
	{
	    dim_typ which_env;
	
	    if((which_env = argc ? getItemID(argv[0], &envs_manager[ENVS_SYNC], ENVS) : itemSelect(ENVS)) == NULL_ENV)
	       return;
	       
	    if(!access(curLayout)->database.con)
			dbEstablishConnection(0, NULL);
	
		if(access(curLayout)->database.con)
			getVarListEx(NULL, which_env);
			
		msprintf(COLOR_USER, "\nIt has been stored in DB the %s:\n%s", suite_c.listsnames[ENVS], listNo(which_env, ENVS)->path);
	    return;
	}
	
	__MSSHELL_WRAPPER_ static void  __lmp_prog retrieveEnv(const sel_typ argc, char ** argv)
	{
		
		dim_typ which_env;
	
	    if((which_env = argc ? getItemID(argv[0], &envs_manager[ENVS_SYNC], ENVS) : itemSelect(ENVS)) == NULL_ENV)
	       return;
	       
		if(!access(curLayout)->database.con)
	        dbEstablishConnection(0, NULL); // (char * [4]){NULL, NULL, NULL, NULL});
	    
		if(access(curLayout)->database.con)
		{
			char tab[MAX_PATH_LENGTH];
			strcpy(tab, listNo(which_env, ENVS)->path);
			(void) strtok(tab, EXTENSION_DOT);
			dim_typ varlistno;
			dim_typ correctness = 0;
			char var[MAX_VARS][VAR_MAXNAMELENGTH] = { NULL_CHAR };
			char val[MAX_VARS][VAR_MAXVALLENGTH] = { NULL_CHAR };
			readVarList(tab, &varlistno, var, val);
				
			for(dim_typ i=0; i<varlistno; ++i)
			{
				sprintf(var[i], "%s=%s;", var[i], val[i]);
				correctness += _parse(var[i], NULL, ((exprValList *)listNo(which_env, ENVS)->data));
			}
			if(correctness == varlistno)
				msprintf(COLOR_USER, "\nIt has been parsed from the DB the %s:\n%s."DEFAULT_VARLIST_FILE_EXTENSION, suite_c.listsnames[ENVS], tab);
			else
				printErr(ERROR_INPUTOUTPUT, "An error occurred during parsing from the DB the %s:\n%s."DEFAULT_VARLIST_FILE_EXTENSION, suite_c.listsnames[ENVS], tab);
		}
	
	    return;
	}
#endif

#endif
