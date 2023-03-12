// syslog_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifndef __DISABLE_SYSLOGMANAGER


// FUNCTIONS DECLARATIONS

__MSSHELL_WRAPPER_ static void   createLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   viewLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   printLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   editLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   updLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   relLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   selLogBufLen(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   flushLogBuf(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void   renLog(const sel_typ argc, char ** argv);


sprog syslog_manager[MAX_SYSLOGMANAGER_PROGS] =
{
    [SYSLOG_OPEN] =
    {
    	CMD_OPENSYSLOG,
        NAME_OPENSYSLOG,
        USAGE_OPENSYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_OPENSYSLOG,
        #endif
        createLog,
        ARGC_OPENSYSLOG,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_CREATE] =
    {
    	CMD_CREATESYSLOG,
        NAME_CREATESYSLOG,
        USAGE_CREATESYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_CREATESYSLOG,
        #endif
        createLog,
        ARGC_CREATESYSLOG,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_READ] =
    {
    	CMD_READSYSLOG,
        NAME_READSYSLOG,
        USAGE_READSYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_READSYSLOG,
        #endif
        viewLog,
        ARGC_READSYSLOG,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_PRINT] =
    {
    	CMD_PRTSYSLOG,
        NAME_PRTSYSLOG,
        USAGE_PRTSYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_PRTSYSLOG,
        #endif
        printLog,
        ARGC_PRTSYSLOG,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_EDIT] =
    {
    	CMD_EDITSYSLOG,
        NAME_EDITSYSLOG,
        USAGE_EDITSYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_EDITSYSLOG,
        #endif
        editLog,
        ARGC_EDITSYSLOG,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_UPDATE] =
    {
    	CMD_UPDSYSLOG,
        NAME_UPDSYSLOG,
        USAGE_UPDSYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDSYSLOG,
        #endif
        updLog,
        ARGC_UPDSYSLOG,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_SAVE] =
    {
    	CMD_SAVESYSLOG,
        NAME_SAVESYSLOG,
        USAGE_SAVESYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_SAVESYSLOG,
        #endif
        createLog,
        ARGC_SAVESYSLOG,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_RELOAD] =
    {
    	CMD_RESETSYSLOG,
        NAME_RESETSYSLOG,
        USAGE_RESETSYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_RESETSYSLOG,
        #endif
        relLog,
        ARGC_RESETSYSLOG,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_SETBUFLEN] =
    {
    	CMD_SETSYSLOGBUFLEN,
        NAME_SETSYSLOGBUFLEN,
        USAGE_SETSYSLOGBUFLEN,
        #ifndef __DISABLE_DATABASE
        LEVEL_SETSYSLOGBUFLEN,
        #endif
        selLogBufLen,
        ARGC_SETSYSLOGBUFLEN,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_EMPTYBUFFER] =
    {
    	CMD_EMPTYSYSLOGBUF,
        NAME_EMPTYSYSLOGBUF,
        USAGE_EMPTYSYSLOGBUF,
        #ifndef __DISABLE_DATABASE
        LEVEL_EMPTYSYSLOGBUF,
        #endif
        flushLogBuf,
        ARGC_EMPTYSYSLOGBUF,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_RENAME] =
    {
    	CMD_RENSYSLOG,
        NAME_RENSYSLOG,
        USAGE_RENSYSLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_RENSYSLOG,
        #endif
        renLog,
        ARGC_RENSYSLOG,
        AUTOMATIC,
        CHILD
    }
};


// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static void   createLog(const sel_typ argc, char ** argv)
{
    char name[MAX_PATH_LENGTH] = NULL_CHAR;
    const bool assert = __pmode__ == SYSLOG_OPEN;
    #ifdef WINOS
    if(isnSett(BOOLS_ITEMSSELECTBYPATH))
    {
        if(!windowsFileHandler(name,
                                "File Log (*."DEFAULT_SYSLOG_FILE_EXTENSION")\0*."DEFAULT_SYSLOG_FILE_EXTENSION"\0Documento di testo (*.txt)\0*.txt\0Tutti i Files (*.*)\0*.*\0",
                                DEFAULT_SYSLOG_FILE_EXTENSION, assert))
        // else
        {
            printErr(14, "Failed System Log Selection");
            return;
        }
    }
    else
    #endif
    {
        if(argc)
        {
            if(!file_exists(argv[0]))
            {
                printErr(2, "Inserted Path refers to non-existing System Log:\n%s", argv[0]);
                printUsage(&syslog_manager[__pmode__]);
                return;
            }
            strcpy(name, argv[0]);
        }
        else
        {
            // char string[MAX_PATH_LENGTH];
            msprintf(COLOR_CREDITS, "Enter desired System Log Path.\n");
            msprintf(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
            PRINTL();
            while(scanf("%s", name) != 1 || name[0] == SCANFEXIT_CHAR || (assert && !file_exists(name)))
            {
                CLEARBUFFER();
                if(name[0] == SCANFEXIT_CHAR) return;
                if(assert && !file_exists(name))
                {
                    printErr(2, "Inserted Path refers to non-existing System Log:\n%s", name);
                    return;
                }
            }
            CLEARBUFFER();
        }
    }

    FILE *fp = NULL;

    if((fp = checkForFHErrors(name, "a")) == NULL)
        return;

    static bool once_executed = false;

    if(!once_executed)
    {
        access(sysLog) = malloc(sizeof(logObj));
        errMem(access(sysLog), VSPACE);
    	access(sysLog)->buffer = calloc(DEFAULT_BUFSIZE, sizeof(char));
    	errMem(access(sysLog)->buffer, VSPACE);
    	strcpy(access(sysLog)->buffer, NULL_CHAR); // initializing log buffer
    	access(sysLog)->buflen = DEFAULT_BUFSIZE;
        once_executed = true;
    }

    strcpy(access(sysLogPath), name);
	msyprintf(COLOR_SYSTEM, "System Log has been %s in the following Path:\n%s\n\n", assert ? "opened" : "saved", name);
    fclose(fp);

    return;
}

__MSSHELL_WRAPPER_ static void   viewLog(const sel_typ argc, char ** argv)
{
    logPrint(access(sysLog));
    return;
}

__MSSHELL_WRAPPER_ static void   printLog(const sel_typ argc, char ** argv)
{
    printFile(access(sysLogPath));
    return;
}

static void editLog(const sel_typ argc, char ** argv)
{
    _editLog(access(sysLogPath));
    return;
}

__MSSHELL_WRAPPER_ static void   updLog(const sel_typ argc, char ** argv)
{
    FILE *fp;

    if((fp = checkForFHErrors(access(sysLogPath), "a")) == NULL)
        return;

    uint64_t i;
    const size_t slen = strlen(access(sysLog)->buffer);

    for(i=0; i<slen; ++i)
        putc(access(sysLog)->buffer[i], fp);

    msyprintf(COLOR_SYSTEM, "System Log:\n%s\nhas been correctly saved on Disk.\n\n", access(sysLogPath));

    fclose(fp);

    _flushLogBuf(access(sysLog));

    return;
}

__MSSHELL_WRAPPER_ static void   relLog(const sel_typ argc, char ** argv)
{
    writeFile(access(sysLogPath));
    return;
}

__MSSHELL_WRAPPER_ static void   selLogBufLen(const sel_typ argc, char ** argv)
{
    size_t buflen = 0;

    if(argc)
    {
	    mpfr_t tmp;
	    if(!parse(argv[0], &tmp))
	    {
	        printUsage(&syslog_manager[SYSLOG_SETBUFLEN]);
	        return;
	    }
	    buflen = mpfr_get_ui(tmp, MPFR_RNDN);
    }
    setLogBufLen(access(sysLogPath), access(sysLog), buflen);
    return;
}

__MSSHELL_WRAPPER_ static void   flushLogBuf(const sel_typ argc, char ** argv)
{
    msyprintf(COLOR_SYSTEM, "\nSystem Log:\n%s\nBuffer has been properly flushed.\n\n", access(sysLogPath));
    _flushLogBuf(access(sysLog));
    return;
}

__MSSHELL_WRAPPER_ static void   renLog(const sel_typ argc, char ** argv)
{
	size_t len;
    char name[MAX_PATH_LENGTH];
    char newname[MAX_PATH_LENGTH];

    CLEARBUFFER();
    strcpy(name, access(sysLogPath));
    CLEARBUFFER();
    len = 0;

    if(argc)
    {
        if(strrchr(argv[0], SCANFEXIT_CHAR) != NULL || (len = strlen(argv[0])) > MAX_PATH_LENGTH || len < 2)
        {
            printErr(20+(18*(len>MAX_PATH_LENGTH)), "Invalid inserted Name: %s", argv[0]);
            printUsage(&syslog_manager[SYSLOG_RENAME]);
            return;
        }
        strcpy(newname, argv[0]);
    }
    else
    {
        msprintf(COLOR_CREDITS, "\nEnter System Log newname.\n");
        msprintf(COLOR_CREDITS, "Enter %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);

        while(scanf("%s", newname) != 1 || newname[0] == SCANFEXIT_CHAR ||
            strrchr(newname, SCANFEXIT_CHAR) != NULL || (len = strlen(newname)) > MAX_PATH_LENGTH || len < 2)
        {
            CLEARBUFFER();
            if(newname[0] == SCANFEXIT_CHAR) return;
            printErr(20+(18*(len>MAX_PATH_LENGTH)), "Invalid inserted name: %s", newname);
        }
    }

    strcat(newname, ".LOG");

    CLEARBUFFER();


    if(frename(name, newname))
    {
        strcpy(access(sysLogPath), newname);
        msyprintf(COLOR_SYSTEM, "System Log:\n%s\nhas been correctly renamed to:\n%s.\n\n", name, newname);
    }

    return;
}

#endif
