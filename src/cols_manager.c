// cols_manager.c 04/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"

#if WINOS
#ifdef ALLOW_COLSMANAGER

__MSSHELL_WRAPPER_ static void _MS__private __system changeColors(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __system colFileLoader(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system backupColFile(const register sel_typ argc, char ** argv);

sprog cols_manager[MAX_COLSMANAGER_PROGS] =
{
    [COLORS_CHANGE] =
    {
        "Edit Program Colors",
        CMD_CHANGECOLORS,
        USAGE_CHANGECOLORS,
        changeColors,
        AUTOMATIC,
        CHILD
    },
    [COLORS_FILESLOADER] =
    {
        "Load Colors Settings File ("DEFAULT_COLORS_FILE_EXTENSION")",
        CMD_COLORSFILESLOADER,
        USAGE_COLORSFILESLOADER,
        colFileLoader,
        AUTOMATIC,
        CHILD
    },
    [COLORS_BACKUPFILES] =
    {
        "Save Colors Settings File",
        CMD_BACKUPCOLORSFILES,
        USAGE_BACKUPCOLORSFILES,
        backupColFile,
        AUTOMATIC,
        CHILD
    }
};

__MSSHELL_WRAPPER_ static void _MS__private __system changeColors(const register sel_typ argc, char ** argv)
{
    dim_typ i, j;
    dim_typ old_col = INIT_COLOR;

    if(argc)
    {
        if(argc > 1)
        {
            ityp tmp;
            if(PARSING_SYSTEM_ALLOWED)
            {
                if((!parse(argv[0], &tmp)) || tmp != (i = (dim_typ)tmp) || i < 0 || i > MAX_COLOR_TYPES)
                {
                    printErr(1, "Inserted Value refers to a bad Color ID");
                    printUsage(&cols_manager[COLORS_CHANGE]);
                    return;
                }

                old_col = access(colors)[i];

                if((!parse(argv[1], &tmp)) || tmp != (j = (dim_typ)tmp) || j < 0 || j > MAX_COLORS)
                {
                    printErr(1, "Inserted Value refers to a bad Color ID");
                    printUsage(&cols_manager[COLORS_CHANGE]);
                    return;
                }
            }
            else
            {
                if((tmp = strtod(argv[0], NULL)) != (i = (dim_typ)tmp) || (tmp = strtod(argv[1], NULL)) != (j = (dim_typ)tmp) || i < 0 || j < 0 || i > MAX_COLOR_TYPES || j > MAX_COLORS)
                {
                    printErr(1, "Inserted Value refers to a bad Color ID");
                    printUsage(&change_settings[COLORS_CHANGE]);
                    return;
                }
            }
        }
        else
        {
            printUsage(&change_settings[COLORS_CHANGE]);
            return;
        }
    }
    else
    {
        if((i = selectListItem(MAX_COLOR_TYPES, MAX_COLOR_TYPES > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
            "Select the Color Type you wish to change", suite_c.colors_types_names)) == MAX_COLOR_TYPES) return;

        old_col = access(colors)[i];

        if((j = selectListItem(MAX_COLORS, MAX_COLORS > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
            "Select the Color to set to the Color Type you've selected to change", suite_c.colors_names)) == MAX_COLORS) return;

    }

    access(colors)[i] = j;
    printf2(COLOR_USER, "%s has been correctly changed from: %s to: %s.\n\n", suite_c.colors_types_names[i], suite_c.colors_names[old_col], suite_c.colors_names[j]);
    return;

}

__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __system colFileLoader(const register sel_typ argc, char ** argv)
{
    char path[MAX_PATH_LENGTH];
    if(isnSett(BOOLS_ITEMSSELECTBYPATH))
    {
        const bool wHandler = windowsFileHandler(path,  "Settings Configuration (*."DEFAULT_COLORS_FILE_EXTENSION")\0*."DEFAULT_COLORS_FILE_EXTENSION"\0Text Documents (*.txt)\0*.txt\0All Files (*.*)\0*.*\0",
                                DEFAULT_COLORS_FILE_EXTENSION, true);
        if(wHandler)
        {
            _colFileLoader(path);
            sprint("%s\nFile has been correctly loaded.\n\n", path);
        }
        else
            printErr(14, "Failed to select Colors Settings File");
    }
    else
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

            printf2(COLOR_CREDITS, "Enter the Path of the "DEFAULT_LAYOUT_FILE_EXTENSION" File you wish to load.\n");
            printf2(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
            PRINTL();

            while(scanf("%s", path) != 1 || path[0] == SCANFEXIT_CHAR || (assert = !file_exists(path)))
            {
                CLEARBUFFER();
                if(path[0] == SCANFEXIT_CHAR) return;
                if(assert)
                {
                    printErr(2, "Inserted Path:\n%s\nrefers to non-existent File", argv[0]);
                    return;
                }
            // mustcreatefile = true;
            }
        }
        _colFileLoader(path);
        sprint("%s\nFile has been correctly loaded.\n\n", path);
    }
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system backupColFile(const register sel_typ argc, char ** argv)
{
    _backupColFile();
    printf2(COLOR_USER, "%s\nColors Settings File has been correctly saved.\n\n", access(colors_path));
    return;
}

#endif
#endif
