// cols_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"

#ifndef __DISABLE_COLSMANAGER

__MSSHELL_WRAPPER_ static void   changeColors(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void   colFileLoader(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   backupColFile(const sel_typ argc, char ** argv);

sprog cols_manager[MAX_COLSMANAGER_PROGS] =
{
    [COLORS_CHANGE] =
    {
    	CMD_CHANGECOLORS,
        NAME_CHANGECOLORS,
        USAGE_CHANGECOLORS,
        #ifndef __DISABLE_DATABASE
        LEVEL_CHANGECOLORS,
        #endif
        changeColors,
        ARGC_CHANGECOLORS,
        AUTOMATIC,
        CHILD
    },
    [COLORS_FILESLOADER] =
    {
    	CMD_COLORSFILESLOADER,
        NAME_COLORSFILESLOADER,
        USAGE_COLORSFILESLOADER,
        #ifndef __DISABLE_DATABASE
        LEVEL_COLORSFILESLOADER,
        #endif
        colFileLoader,
        ARGC_COLORSFILESLOADER,
        AUTOMATIC,
        CHILD
    },
    [COLORS_BACKUPFILES] =
    {
    	CMD_BACKUPCOLORSFILES,
        NAME_BACKUPCOLORSFILES,
        USAGE_BACKUPCOLORSFILES,
        LEVEL_BACKUPCOLORSFILES,
        backupColFile,
        ARGC_BACKUPCOLORSFILES,
        AUTOMATIC,
        CHILD
    }
};

__MSSHELL_WRAPPER_ static void   changeColors(const sel_typ argc, char ** argv)
{
    dim_typ i, j;
    dim_typ old_col = INIT_COLOR;

    if(argc)
    {
        if(argc > 1)
        {
            mpfr_t tmp;
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (i = mpfr_get_ui(tmp, MPFR_RNDN))) || i < 0 || i > MAX_COLOR_TYPES)
            {
            	mpfr_clear(tmp);
                printErr(1, "Inserted Value refers to a bad Color ID");
                printUsage(&cols_manager[COLORS_CHANGE]);
                return;
            }

            old_col = access(colors)[i];

            if((!parse(argv[1], &tmp)) || mpfr_cmp_ui(tmp, (j = mpfr_get_ui(tmp, MPFR_RNDN))) || j < 0 || j > MAX_COLORS)
            {
            	mpfr_clear(tmp);
                printErr(1, "Inserted Value refers to a bad Color ID");
                printUsage(&cols_manager[COLORS_CHANGE]);
                return;
            }
            
            mpfr_clear(tmp);
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
    msprintf(COLOR_USER, "%s has been correctly changed from: %s to: %s.\n\n", suite_c.colors_types_names[i], suite_c.colors_names[old_col], suite_c.colors_names[j]);
    return;

}

__MSSHELL_WRAPPER_ __WINCALL static void   colFileLoader(const sel_typ argc, char ** argv)
{
    char path[MAX_PATH_LENGTH];
    #ifdef WINOS
		if(isnSett(BOOLS_ITEMSSELECTBYPATH))
		{
		    const bool wHandler = windowsFileHandler(path,  "Settings Configuration (*."DEFAULT_COLORS_FILE_EXTENSION")\0*."DEFAULT_COLORS_FILE_EXTENSION"\0Text Documents (*.txt)\0*.txt\0All Files (*.*)\0*.*\0",
		                            DEFAULT_COLORS_FILE_EXTENSION, true);
		    if(wHandler)
		    {
		        _colFileLoader(path);
		        msyprintf(COLOR_SYSTEM, "%s\nFile has been correctly loaded.\n\n", path);
		    }
		    else
		        printErr(14, "Failed to select Colors Settings File");
		}
		else
	#endif
    {
        if(argc)
        {
            if(!file_exists(argv[0]))
            {
                printErr(2, "Inserted Path:\n%s\nrefers to non-existing File", argv[0]);
                return;
            }
            strcpy(path, argv[0]);
        }
        else
        {
            bool assert;

            msprintf(COLOR_CREDITS, "Enter the Path of the "DEFAULT_LAYOUT_FILE_EXTENSION" File you wish to load.\n");
            msprintf(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
            PRINTL();

            while(scanf("%s", path) != 1 || path[0] == SCANFEXIT_CHAR || (assert = !file_exists(path)))
            {
                CLEARBUFFER();
                if(path[0] == SCANFEXIT_CHAR) return;
                if(assert)
                {
                    printErr(2, "Inserted Path:\n%s\nrefers to non-existing File", argv[0]);
                    return;
                }
            // mustcreatefile = true;
            }
        }
        _colFileLoader(path);
        msyprintf(COLOR_SYSTEM, "%s\nFile has been correctly loaded.\n\n", path);
    }
    return;
}

__MSSHELL_WRAPPER_ static void   backupColFile(const sel_typ argc, char ** argv)
{
    _backupColFile();
    msprintf(COLOR_USER, "%s\nColors Settings File has been correctly saved.\n\n", access(colors_path));
    return;
}

#endif
