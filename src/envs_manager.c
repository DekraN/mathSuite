// envs_manager.c 10/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifdef ALLOW_VARLISTMANAGER

// FUNCTIONS DECLARATIONS


__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog setCurEnv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __lmp_prog createEnv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog viewEnv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog printEnv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog updEnv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog updAllEnvs(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog delEnv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog relEnv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog syncEnv(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog renEnv(const register sel_typ argc, char ** argv);

sprog envs_manager[MAX_ENVSMANAGER_PROGS] =
{
    [ENVS_SETCURRENT] =
    {
        "Select Current Env",
        CMD_SETCURENV,
        USAGE_SETCURENV,
        setCurEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_OPEN] =
    {
        "Open Env ("DEFAULT_VARLIST_FILE_EXTENSION")",
        CMD_OPENENV,
        USAGE_OPENENV,
        createEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_CREATE] =
    {
        "Create new Env ("DEFAULT_VARLIST_FILE_EXTENSION")",
        CMD_CREATEENV,
        USAGE_CREATEENV,
        createEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_READ] =
    {
        "View Env",
        CMD_READENV,
        USAGE_READENV,
        viewEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_PRINT] =
    {
        "Print Env",
        CMD_PRTENV,
        USAGE_PRTENV,
        printEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_UPDATE] =
    {
        "Save Env",
        CMD_UPDENV,
        USAGE_UPDENV,
        updEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_UPDATEALL] =
    {
        "Save all Envs",
        CMD_UPDALLENVS,
        USAGE_UPDALLENVS,
        updAllEnvs,
        AUTOMATIC,
        CHILD
    },
    [ENVS_SAVE] =
    {
        "Save Env as ("DEFAULT_VARLIST_FILE_EXTENSION")",
        CMD_SAVEENV,
        USAGE_SAVEENV,
        createEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_DELETE] =
    {
        "Delete Env",
        CMD_DELENV,
        USAGE_DELENV,
        delEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_DELETEPHYSICAL] =
    {
        "Delete Physical Env",
        CMD_DELENV2,
        USAGE_DELENV2,
        delEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_DELETEALL] =
    {
        "Delete all Envs",
        CMD_DELALLENVS,
        USAGE_DELALLENVS,
        delEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_DELETEALLPHYSICALS] =
    {
        "Delete all Physical Envs",
        CMD_DELALLENVS2,
        USAGE_DELALLENVS2,
        delEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_RELOAD] =
    {
        "Reset Env",
        CMD_RESETENV,
        USAGE_RESETENV,
        relEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_SYNC] =
    {
        "Synchronize Env",
        CMD_SYNCENV,
        USAGE_SYNCENV,
        syncEnv,
        AUTOMATIC,
        CHILD
    },
    [ENVS_RENAME] =
    {
        "Rename Env ("DEFAULT_VARLIST_FILE_EXTENSION")",
        CMD_RENENV,
        USAGE_RENENV,
        renEnv,
        AUTOMATIC,
        CHILD
    }
};

// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog setCurEnv(const register sel_typ argc, char ** argv)
{
    setCurItem(argc ? getItemID(argv[0], &envs_manager[ENVS_SETCURRENT], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static inline void _MS__private __lmp_prog createEnv(const register sel_typ argc, char ** argv)
{
    createItem(argc ? argv[0] : NULL, ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog viewEnv(const register sel_typ argc, char ** argv)
{
    viewItem(argc ? getItemID(argv[0], &envs_manager[ENVS_READ], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog printEnv(const register sel_typ argc, char ** argv)
{
    printListItem(argc ? getItemID(argv[0], &envs_manager[ENVS_PRINT], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog updEnv(const register sel_typ argc, char ** argv)
{
    updItem(ENVS);
    return;
}
__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog updAllEnvs(const register sel_typ argc, char ** argv)
{
    updAll(ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog delEnv(const register sel_typ argc, char ** argv)
{
    delItem(argc ? getItemID(argv[0], &envs_manager[__pmode__], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog relEnv(const register sel_typ argc, char ** argv)
{
    relItem(argc ? getItemID(argv[0], &envs_manager[ENVS_RELOAD], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog syncEnv(const register sel_typ argc, char ** argv)
{
    dim_typ which_env;

    if((which_env = argc ? getItemID(argv[0], &envs_manager[ENVS_SYNC], ENVS) : itemSelect(ENVS)) == NULL_ENV)
       return;


    char name[MAX_PATH_LENGTH];
    strcpy(name, listNo(which_env, ENVS)->path);

    FILE *fp;
    if((fp = checkForFHErrors(name, "a")) == NULL)
        return;

    getVarList(fp);
    fclose(fp);

    printf2(COLOR_USER, "\nThe Env:\n");

    printf2(COLOR_CREDITS, name);
    printf2(COLOR_USER, "\nhas been correctly synchronized with the Env:\n");
    printf2(COLOR_CREDITS, listNo(access(lists)[ENVS].cur_item, ENVS)->path);


    PRINT2N();

    if(isSett(BOOLS_AUTOSETCURITEM))
        passToItem(which_env, ENVS, true);

    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog renEnv(const register sel_typ argc, char ** argv)
{
    renItem(argc > 1 ? argv[1] : NULL, argc ? getItemID(argv[0], &envs_manager[ENVS_RENAME], ENVS) : getItemsListNo(ENVS), ENVS);
    return;
}

#endif
