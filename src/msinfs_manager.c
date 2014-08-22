// msinfs_manager.c 20/08/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifdef ALLOW_MSINFSMANAGER

__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __system msinfLoader(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __system msinfCreate(const register sel_typ argc, char ** argv);

sprog msinfs_manager[MAX_MSINFSMANAGER_PROGS] =
{
    [MSINFS_LOADER] =
    {
        "Load Startup File ("DEFAULT_PATHLIST_FILE_EXTENSION")",
        CMD_MSINFLOADER,
        USAGE_MSINFLOADER,
        msinfLoader,
        AUTOMATIC,
        CHILD
    },
    [MSINFS_CREATE] =
    {
        "Create new Startup File ("DEFAULT_PATHLIST_FILE_EXTENSION")",
        CMD_MSINFCREATE,
        USAGE_MSINFCREATE,
        msinfCreate,
        AUTOMATIC,
        CHILD
    }
};

__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __system msinfLoader(const register sel_typ argc, char ** argv)
{
    char path[MAX_PATH_LENGTH];
    #if WINOS
        if(isnSett(BOOLS_ITEMSSELECTBYPATH))
        {
            const bool wHandler = windowsFileHandler(path,  "MathSuite Informations (*."DEFAULT_PATHLIST_FILE_EXTENSION")\0*."DEFAULT_PATHLIST_FILE_EXTENSION"\0Text Document (*.txt)\0*.txt\0Generic DAT Files (*.DAT)\0*.DAT\0All Files (*.*)\0*.*\0",
                                    DEFAULT_PATHLIST_FILE_EXTENSION, true);
            if(wHandler)
            {
                if(_msinfLoader(path))
                    sprint("%s\nFile has been correctly loaded.\n\n", path);
                else
                    printErr(2, "An error during:\n%s\nFile opening process might have occurred", path);
            }
            else
                printErr(14, "Failed to select "DEFAULT_PATHLIST_FILE_EXTENSION" File");
        }
        else
    #endif
        {
            if(argc)
            {
                if(!file_exists(argv[0]))
                {
                    printErr(2, "Inserted Path:\n%s\nrefers to non-existent File", argv[0]);
                    return;
                }
                strcpy(path, argv[0]);
            }
            else
            {
                bool assert;

                printf2(COLOR_CREDITS, "Enter "DEFAULT_PATHLIST_FILE_EXTENSION" Path of the File you wish to load.\n");
                printf2(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
                PRINTL();

                while(scanf("%s", path) != 1 || path[0] == SCANFEXIT_CHAR || (assert = !file_exists(path)))
                {
                    CLEARBUFFER();
                    if(path[0] == SCANFEXIT_CHAR) return;
                    if(assert)
                    {
                        printErr(2, "Inserted Path:\n%s\nrefers to non-existent File", path);
                        return;
                    }
                // mustcreatefile = true;
                }
            }
            if(_msinfLoader(path))
                sprint("%s\nFile has been correctly loaded.\n\n", path);
            else
                printErr(2, "An error during:\n%s\nFile opening process might have occurred", path);
        }

    return;
}

__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __system msinfCreate(const register sel_typ argc, char ** argv)
{
    char path[MAX_PATH_LENGTH];
    #if WINOS
        if(isnSett(BOOLS_ITEMSSELECTBYPATH))
        {
            const bool wHandler = windowsFileHandler(path,  "MathSuite Informations (*."DEFAULT_PATHLIST_FILE_EXTENSION")\0*."DEFAULT_PATHLIST_FILE_EXTENSION"\0Text Document (*.txt)\0*.txt\0Generic DAT Files (*.DAT)\0*.DAT\0All Files (*.*)\0*.*\0",
                                    DEFAULT_PATHLIST_FILE_EXTENSION, false);
            if(wHandler)
            {
                if(_msinfCreate(path))
                    sprint("%s\nFile has been correctly loaded.\n\n", path);
                else
                    printErr(9, "An error during:\n%s\nFile opening process might have occurred", path);
            }
            else
                printErr(14, "Failed to select "DEFAULT_PATHLIST_FILE_EXTENSION" File");
        }
        else
    #endif
        {
            if(argc)
            {
                if(file_exists(argv[0]))
                {
                    printErr(1, "Inserted Path:\n%s\nrefers to an already existent File", argv[0]);
                    return;
                }
                strcpy(path, argv[0]);
            }
            else
            {
                bool assert;

                printf2(COLOR_CREDITS, "Enter "DEFAULT_PATHLIST_FILE_EXTENSION" Path of the File you wish to load.\n");
                printf2(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
                PRINTL();

                while(scanf("%s", path) != 1 || path[0] == SCANFEXIT_CHAR || (assert = file_exists(path)))
                {
                    CLEARBUFFER();
                    if(path[0] == SCANFEXIT_CHAR) return;
                    if(assert)
                    {
                        printErr(1, "Inserted Path:\n%s\nrefers to an already existent File", path);
                        return;
                    }
                // mustcreatefile = true;
                }
            }
            if(_msinfCreate(path))
                sprint("%s\nFile has been correctly saved.\n\n", path);
            else
                printErr(9, "An error during:\n%s\nFile opening process might have occurred.\n\n", path);
        }
    return;
}

#endif // ALLOW_MSINFSMANAGER
