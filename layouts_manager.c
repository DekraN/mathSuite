// layouts_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifndef __DISABLE_LAYOUTSMANAGER

// FUNCTIONS DECLARATIONS

__MSSHELL_WRAPPER_ static void   __lmp_prog setCurLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void   __lmp_prog createLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   __lmp_prog viewLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   __lmp_prog printLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   __lmp_prog updLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   __lmp_prog updAllLayouts(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   __lmp_prog delLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   __lmp_prog relLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   __lmp_prog renLayout(const sel_typ argc, char ** argv);

sprog layouts_manager[MAX_LAYOUTSMANAGER_PROGS] =
{
    [LAYOUTS_SETCURRENT] =
    {
    	CMD_SETCURLAYOUT,
        NAME_SETCURLAYOUT,
        USAGE_SETCURLAYOUT,
        #ifndef __DISABLE_DATABASE
    	LEVEL_SETCURLAYOUT,
    	#endif
        setCurLayout,
        ARGC_SETCURLAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_OPEN] =
    {
    	CMD_OPENLAYOUT,
        NAME_OPENLAYOUT,
        USAGE_OPENLAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_OPENLAYOUT,
        #endif
        createLayout,
        ARGC_OPENLAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_CREATE] =
    {
    	CMD_CREATELAYOUT,
        NAME_CREATELAYOUT,
        USAGE_CREATELAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_CREATELAYOUT,
        #endif
        createLayout,
        ARGC_CREATELAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_READ] =
    {
    	CMD_READLAYOUT,
        NAME_READLAYOUT,
        USAGE_READLAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_READLAYOUT,
        #endif
        viewLayout,
        ARGC_READLAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_PRINT] =
    {
    	CMD_PRTLAYOUT,
        NAME_PRTLAYOUT,
        USAGE_PRTLAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_PRTLAYOUT,
        #endif
        printLayout,
        ARGC_PRTLAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_UPDATE] =
    {
    	CMD_UPDLAYOUT,
        NAME_UPDLAYOUT,
        USAGE_UPDLAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDLAYOUT,
        #endif
        updLayout,
        ARGC_UPDLAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_UPDATEALL] =
    {
    	CMD_UPDALLLAYOUTS,
        NAME_UPDALLLAYOUTS,
        USAGE_UPDALLLAYOUTS,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDALLLAYOUTS,
        #endif
        updAllLayouts,
        ARGC_UPDALLLAYOUTS,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_SAVE] =
    {
    	CMD_SAVELAYOUT,
        NAME_SAVELAYOUT,
        USAGE_SAVELAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_SAVELAYOUT,
        #endif
        createLayout,
        ARGC_SAVELAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_DELETE] =
    {
    	CMD_DELLAYOUT,
        NAME_DELLAYOUT,
        USAGE_DELLAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELLAYOUT,
        #endif
        delLayout,
        ARGC_DELLAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_DELETEPHYSICAL] =
    {
    	CMD_DELLAYOUT2,
        NAME_DELLAYOUT2,
        USAGE_DELLAYOUT2,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELLAYOUT2,
        #endif
        delLayout,
        ARGC_DELLAYOUT2,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_DELETEALL] =
    {
    	CMD_DELALLLAYOUTS,
        NAME_DELALLLAYOUTS,
        USAGE_DELALLLAYOUTS,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELALLLAYOUTS,
        #endif
        delLayout,
        ARGC_DELALLLAYOUTS,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_DELETEALLPHYSICALS] =
    {
    	CMD_DELALLLAYOUTS2,
        NAME_DELALLLAYOUTS2,
        USAGE_DELALLLAYOUTS2,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELALLLAYOUTS2,
        #endif
        delLayout,
        ARGC_DELALLLAYOUTS2,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_RELOAD] =
    {
    	CMD_RESETLAYOUT,
        NAME_RESETLAYOUT,
        USAGE_RESETLAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_RESETLAYOUT,
        #endif
        relLayout,
        ARGC_RESETLAYOUT,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_RENAME] =
    {
    	CMD_RENLAYOUT,
        NAME_RENLAYOUT,
        USAGE_RENLAYOUT,
        #ifndef __DISABLE_DATABASE
        LEVEL_RENLAYOUT,
        #endif
        renLayout,
        ARGC_RENLAYOUT,
        AUTOMATIC,
        CHILD
    }
};

// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static inline void   __lmp_prog setCurLayout(const sel_typ argc, char ** argv)
{
    setCurItem(argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_SETCURRENT], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static inline void   __lmp_prog createLayout(const sel_typ argc, char ** argv)
{
    createItem(argc ? argv[0] : NULL, LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void   __lmp_prog viewLayout(const sel_typ argc, char ** argv)
{
    viewItem(argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_READ], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void   __lmp_prog printLayout(const sel_typ argc, char ** argv)
{
    printListItem(argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_PRINT], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void   __lmp_prog updLayout(const sel_typ argc, char ** argv)
{
    updItem(LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void   __lmp_prog updAllLayouts(const sel_typ argc, char ** argv)
{
    updAll(LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void   __lmp_prog delLayout(const sel_typ argc, char ** argv)
{
    delItem(argc ? getItemID(argv[0], &layouts_manager[__pmode__], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void   __lmp_prog relLayout(const sel_typ argc, char ** argv)
{
    relItem(argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_RELOAD], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void   __lmp_prog renLayout(const sel_typ argc, char ** argv)
{
    renItem(argc > 1 ? argv[1] : NULL, argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_RENAME], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

#endif
