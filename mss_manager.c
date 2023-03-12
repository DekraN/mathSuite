// mss_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifndef __DISABLE_MSSMANAGER

__MSSHELL_WRAPPER_ static void   handleCmdLine(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   execScriptFiles(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   showUsage(const sel_typ argc, char ** argv);

sprog mss_manager[MAX_MSSMANAGER_PROGS] =
{
    [MSSMANAGER_COMMANDLINE] =
    {
    	CMD_CMDLINE,
        NAME_CMDLINE,
        USAGE_CMDLINE,
        #ifndef __DISABLE_DATABASE
        LEVEL_CMDLINE,
        #endif
        handleCmdLine,
        ARGC_CMDLINE,
        BY_USER, /// AUTOMATIC
        CHILD
    },
    [MSSMANAGER_EXECSCRIPTFILES] =
    {
    	CMD_EXECMSSFILES,
        NAME_EXECMSSFILES,
        USAGE_EXECMSSFILES,
        #ifndef __DISABLE_DATABASE
        LEVEL_EXECMSSFILES,
        #endif
        execScriptFiles,
        ARGC_EXECMSSFILES,
        AUTOMATIC,
        CHILD
    },
    [MSSMANAGER_SHOWUSAGES] =
    {
    	CMD_SHOWUSAGES,
        NAME_SHOWUSAGES,
        USAGE_SHOWUSAGES,
        #ifndef __DISABLE_DATABASE
        LEVEL_SHOWUSAGES,
        #endif
        showUsage,
        ARGC_SHOWUSAGES,
        AUTOMATIC,
        CHILD
    }
};

__MSSHELL_WRAPPER_ static void   handleCmdLine(const sel_typ argc, char ** argv)
{
	char str[INFO_STRING];
    msprintf(COLOR_CREDITS, "Enter one Programs Macro-List Command.\n\n");
 (void) gets(str);

    if(!str[0])
        return;
    	
	_handleCmdLine(str);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static void   execScriptFiles(const sel_typ argc, char ** argv)
{
    char path[MAX_PATH_LENGTH];
    #ifdef WINOS
        if(isnSett(BOOLS_ITEMSSELECTBYPATH))
        {
            const bool wHandler = windowsFileHandler(path,  "MathSuite ScriptFiles (*."DEFAULT_SCRIPTFILES_EXTENSION")\0*."DEFAULT_SCRIPTFILES_EXTENSION"\0Text Documents (*.txt)\0*.txt\0All Files (*.*)\0*.*\0",
                                    DEFAULT_SCRIPTFILES_EXTENSION, true);
            if(wHandler)
            {
                if(_execScriptFiles(path))
                    msyprintf(COLOR_SYSTEM, "\nScriptfile:\n%s.\nhas been correctly executed.\n\n", path);
                else
                printErr(2, "An error during:\n%s\nScriptfile opening process might have occurred", path);
            }
            else
                printErr(14, "Failed to select "DEFAULT_SCRIPTFILES_EXTENSION" File");
        }
        else
    #endif
        {
            if(argc)
            {
                if(!file_exists(argv[0]))
                {
                    printErr(2, "Inserted Path:\n%s\nrefers to non-existing Scriptfile", argv[0]);
                    return;
                }
                strcpy(path, argv[0]);
            }
            else
            {
                bool assert;

                msprintf(COLOR_CREDITS, "Enter Path of the "DEFAULT_SCRIPTFILES_EXTENSION" you wish to load.\n");
                msprintf(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
                PRINTL();

                while(scanf("%s", path) != 1 || path[0] == SCANFEXIT_CHAR || (assert = !file_exists(path)))
                {
                    CLEARBUFFER();
                    if(path[0] == SCANFEXIT_CHAR) return;
                    if(assert)
                    {
                        printErr(2, "Inserted Path:\n%s\nrefers to non-existing Scriptfile", path);
                        return;
                    }
                // mustcreatefile = true;
                }
            }
            if(_execScriptFiles(path))
                msyprintf(COLOR_SYSTEM, "\nScriptfile:\n%s\nhas been correctly executed.\n\n", path);
            else
                printErr(2, "An error during:\n%s\nScriptfile opening process might have occurred", path);
        }
    return;
}

__MSSHELL_WRAPPER_ static void   showUsage(const sel_typ argc, char ** argv)
{
    if(argc)
    {
        sprog * const prog = searchProgram(argv[0]);
        
        if(prog)
            _showUsage(prog);
        else
            printErr(1, "SubProgram hasn't been found in Program Macro-List");
            
        return;
    }

    msprintf(COLOR_USER, "\n\n*** Program MACRO-LIST ***\n");
    SHOWPAUSEMESSAGE();
    PRINTL();

    dim_typ i;

    for(i=0; i<MAX_PROGRAMMI; ++i)
    {
        _showUsage(&main_menu[i]);
        if(catchPause()) return;
    }

    for(i=0; i<MAX_ADVCALC_PROGS; ++i)
    {
        _showUsage(&adv_calc[i]);
        if(catchPause()) return;
    }

    for(i=0; i<MAX_ALGEBRA_OPERATIONS; ++i)
    {
        _showUsage(&alg_operations[i]);
        if(catchPause()) return;
    }

	#ifndef __DISABLE_VARLISTMANAGER
	    for(i=0; i<MAX_ENVSMANAGER_PROGS; ++i)
	    {
	        _showUsage(&envs_manager[i]);
	        if(catchPause()) return;
	    }
	#endif
	
	#ifndef __DISABLE_MATMANAGER
	    for(i=0; i<MAX_MATMANAGER_PROGS; ++i)
	    {
	        _showUsage(&mat_manager[i]);
	        if(catchPause()) return;
	    }
	#endif
    
    #ifndef __DISABLE_LOGSMANAGER
	    for(i=0; i<MAX_LOGSMANAGER_PROGS; ++i)
	    {
	        _showUsage(&logs_manager[i]);
	        if(catchPause()) return;
	    }
	#endif

	#ifndef __DISABLE_SYSLOGMANAGER
	    for(i=0; i<MAX_SYSLOGMANAGER_PROGS; ++i)
	    {
	        _showUsage(&syslog_manager[i]);
	        if(catchPause()) return;
	    }
	#endif

	#ifndef __DISABLE_LAYOUTSMANAGER
	    for(i=0; i<MAX_LAYOUTSMANAGER_PROGS; ++i)
	    {
	        _showUsage(&layouts_manager[i]);
	        if(catchPause()) return;
	    }
	#endif

    for(i=0; i<MAX_SETTINGS; ++i)
    {
        _showUsage(&change_settings[i]);
        if(catchPause()) return;
    }

	#ifndef __DISABLE_COLSMANAGER
	    for(i=0; i<MAX_COLSMANAGER_PROGS; ++i)
	    {
	        _showUsage(&cols_manager[i]);
	        if(catchPause()) return;
	    }
	#endif

	#ifndef __DISABLE_LFSMANAGER
	    for(i=0; i<MAX_LFSMANAGER_PROGS; ++i)
	    {
	        _showUsage(&lfs_manager[i]);
	        if(catchPause()) return;
	    }
    #endif
    
    for(i=0; i<MAX_MSSMANAGER_PROGS; ++i)
    {
    	_showUsage(&mss_manager[i]);
    	if(catchPause()) return;
    }

    return;
}

#endif
