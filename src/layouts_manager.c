// layouts_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifdef ALLOW_LAYOUTSMANAGER

// FUNCTIONS DECLARATIONS

__MSSHELL_WRAPPER_ static void _MS__private __system __lmp_prog setCurLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __system __lmp_prog createLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system __lmp_prog viewLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system __lmp_prog printLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system __lmp_prog updLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system __lmp_prog updAllLayouts(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system __lmp_prog delLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system __lmp_prog relLayout(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system __lmp_prog renLayout(const sel_typ argc, char ** argv);

sprog layouts_manager[MAX_LAYOUTSMANAGER_PROGS] =
{
    [LAYOUTS_SETCURRENT] =
    {
        "Select Current Layout",
        CMD_SETCURLAYOUT,
        USAGE_SETCURLAYOUT,
        setCurLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_OPEN] =
    {
        "Open Layout ("DEFAULT_LAYOUT_FILE_EXTENSION")",
        CMD_OPENLAYOUT,
        USAGE_OPENLAYOUT,
        createLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_CREATE] =
    {
        "Create new Layout ("DEFAULT_LAYOUT_FILE_EXTENSION")",
        CMD_CREATELAYOUT,
        USAGE_CREATELAYOUT,
        createLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_READ] =
    {
        "View Layout",
        CMD_READLAYOUT,
        USAGE_READLAYOUT,
        viewLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_PRINT] =
    {
        "Print Layout",
        CMD_PRTLAYOUT,
        USAGE_PRTLAYOUT,
        printLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_UPDATE] =
    {
        "Save Layout",
        CMD_UPDLAYOUT,
        USAGE_UPDLAYOUT,
        updLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_UPDATEALL] =
    {
        "Save all Layouts",
        CMD_UPDALLLAYOUTS,
        USAGE_UPDALLLAYOUTS,
        updAllLayouts,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_SAVE] =
    {
        "Save Layout as ("DEFAULT_LAYOUT_FILE_EXTENSION")",
        CMD_SAVELAYOUT,
        USAGE_SAVELAYOUT,
        createLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_DELETE] =
    {
        "Delete Layout",
        CMD_DELLAYOUT,
        USAGE_DELLAYOUT,
        delLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_DELETEPHYSICAL] =
    {
        "Delete Physical Layout",
        CMD_DELLAYOUT2,
        USAGE_DELLAYOUT2,
        delLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_DELETEALL] =
    {
        "Delete all Layouts",
        CMD_DELALLLAYOUTS,
        USAGE_DELALLLAYOUTS,
        delLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_DELETEALLPHYSICALS] =
    {
        "Delete all Physical Layouts",
        CMD_DELALLLAYOUTS2,
        USAGE_DELALLLAYOUTS2,
        delLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_RELOAD] =
    {
        "Reset Layout",
        CMD_RESETLAYOUT,
        USAGE_RESETLAYOUT,
        relLayout,
        AUTOMATIC,
        CHILD
    },
    [LAYOUTS_RENAME] =
    {
        "Rename Layout ("DEFAULT_LAYOUT_FILE_EXTENSION")",
        CMD_RENLAYOUT,
        USAGE_RENLAYOUT,
        renLayout,
        AUTOMATIC,
        CHILD
    }
};

// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static inline void _MS__private __system __lmp_prog setCurLayout(const sel_typ argc, char ** argv)
{
    setCurItem(argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_SETCURRENT], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static inline void _MS__private __system __lmp_prog createLayout(const sel_typ argc, char ** argv)
{
    createItem(argc ? argv[0] : NULL, LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __system __lmp_prog viewLayout(const sel_typ argc, char ** argv)
{
    viewItem(argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_READ], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __system __lmp_prog printLayout(const sel_typ argc, char ** argv)
{
    printListItem(argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_PRINT], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __system __lmp_prog updLayout(const sel_typ argc, char ** argv)
{
    updItem(LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __system __lmp_prog updAllLayouts(const sel_typ argc, char ** argv)
{
    updAll(LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __system __lmp_prog delLayout(const sel_typ argc, char ** argv)
{
    delItem(argc ? getItemID(argv[0], &layouts_manager[__pmode__], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __system __lmp_prog relLayout(const sel_typ argc, char ** argv)
{
    relItem(argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_RELOAD], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __system __lmp_prog renLayout(const sel_typ argc, char ** argv)
{
    renItem(argc > 1 ? argv[1] : NULL, argc ? getItemID(argv[0], &layouts_manager[LAYOUTS_RENAME], LAYOUTS) : getItemsListNo(LAYOUTS), LAYOUTS);
    return;
}

#endif
