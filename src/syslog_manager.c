// syslog_manager.c 20/08/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifdef ALLOW_SYSLOGMANAGER


// FUNCTIONS DECLARATIONS

__MSSHELL_WRAPPER_ static void _MS__private __system createLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system viewLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system printLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system editLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system updLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system relLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system selLogBufLen(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system flushLogBuf(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __system renLog(const register sel_typ argc, char ** argv);


sprog syslog_manager[MAX_SYSLOGMANAGER_PROGS] =
{
    [SYSLOG_OPEN] =
    {
        "Open System Log ("DEFAULT_SYSLOG_FILE_EXTENSION")",
        CMD_OPENSYSLOG,
        USAGE_OPENSYSLOG,
        createLog,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_CREATE] =
    {
        "Create new System Log ("DEFAULT_SYSLOG_FILE_EXTENSION")",
        CMD_CREATESYSLOG,
        USAGE_CREATESYSLOG,
        createLog,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_READ] =
    {
        "View System Log",
        CMD_READSYSLOG,
        USAGE_READSYSLOG,
        viewLog,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_PRINT] =
    {
        "Print System Log",
        CMD_PRTSYSLOG,
        USAGE_PRTSYSLOG,
        printLog,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_EDIT] =
    {
        "Edit System Log ("DEFAULT_SYSLOG_FILE_EXTENSION")",
        CMD_EDITSYSLOG,
        USAGE_EDITSYSLOG,
        editLog,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_UPDATE] =
    {
        "Save System Log",
        CMD_UPDSYSLOG,
        USAGE_UPDSYSLOG,
        updLog,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_SAVE] =
    {
        "Save System Log as ("DEFAULT_SYSLOG_FILE_EXTENSION")",
        CMD_SAVESYSLOG,
        USAGE_SAVESYSLOG,
        createLog,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_RELOAD] =
    {
        "Reset System Log",
        CMD_RESETSYSLOG,
        USAGE_RESETSYSLOG,
        relLog,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_SETBUFLEN] =
    {
        "Select System Log Buffer Length",
        CMD_SETSYSLOGBUFLEN,
        USAGE_SETSYSLOGBUFLEN,
        selLogBufLen,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_EMPTYBUFFER] =
    {
        "Flush System Log Buffer",
        CMD_EMPTYSYSLOGBUF,
        USAGE_EMPTYSYSLOGBUF,
        flushLogBuf,
        AUTOMATIC,
        CHILD
    },
    [SYSLOG_RENAME] =
    {
        "Rename System Log ("DEFAULT_SYSLOG_FILE_EXTENSION")",
        CMD_RENSYSLOG,
        USAGE_RENSYSLOG,
        renLog,
        AUTOMATIC,
        CHILD
    }
};


// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static void _MS__private __system createLog(const register sel_typ argc, char ** argv)
{
    char name[MAX_PATH_LENGTH] = NULL_CHAR;
    const bool assert = __pmode__ == SYSLOG_OPEN;
    #if WINOS
    if(isnSett(BOOLS_ITEMSSELECTBYPATH))
    {
        if(!windowsFileHandler(name,
                                "File Log (*."DEFAULT_SYSLOG_FILE_EXTENSION")\0*."DEFAULT_SYSLOG_FILE_EXTENSION"\0Documento di testo (*.txt)\0*.txt\0File DAT Generici (*.DAT)\0*.DAT\0Tutti i Files (*.*)\0*.*\0",
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
                printErr(2, "Inserted Path refers to non-existent System Log:\n%s", argv[0]);
                printUsage(&syslog_manager[__pmode__]);
                return;
            }
            strcpy(name, argv[0]);
        }
        else
        {
            // char string[MAX_PATH_LENGTH];
            printf2(COLOR_CREDITS, "Enter desired System Log Path.\n");
            printf2(COLOR_CREDITS, "or insert %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);
            PRINTL();
            while(scanf("%s", name) != 1 || name[0] == SCANFEXIT_CHAR || (assert && !file_exists(name)));
            {
                CLEARBUFFER();
                if(name[0] == SCANFEXIT_CHAR) return;
                if(assert && !file_exists(name))
                {
                    printErr(2, "Inserted Path refers to non-existent System Log:\n%s", name);
                    return;
                }
            }
            CLEARBUFFER();
        }
    }

    FILE *fp = NULL;

    if((fp = checkForFHErrors(name, "w")) == NULL)
        return;

    static bool once_executed = false;

    if(once_executed == false)
    {
        access(sysLog) = malloc(sizeof(logObj));
        errMem(access(sysLog), VSPACE);
        once_executed = true;
    }

    sprint("System Log has been %s in the following Path:\n%s\n\n", assert ? "opened" : "saved", name);
    strcpy(access(sysLogPath), name);

    fclose(fp);

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system viewLog(const register sel_typ argc, char ** argv)
{
    logPrint(access(sysLog));
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system printLog(const register sel_typ argc, char ** argv)
{
    printFile(access(sysLogPath));
    return;
}

static void editLog(const register sel_typ argc, char ** argv)
{
    _editLog(access(sysLogPath));
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system updLog(const register sel_typ argc, char ** argv)
{
    FILE *fp;

    if((fp = checkForFHErrors(access(sysLogPath), "a")) == NULL)
        return;

    uint64_t i;
    const size_t slen = strlen(access(sysLog)->buffer);

    for(i=0; i<slen; ++i)
        putc(access(sysLog)->buffer[i], fp);

    sprint("System Log:\n%s\nhas been correctly saved on Disk.\n\n", access(sysLogPath));

    fclose(fp);

    _flushLogBuf(access(sysLog));

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system relLog(const register sel_typ argc, char ** argv)
{
    writeFile(access(sysLogPath));
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system selLogBufLen(const register sel_typ argc, char ** argv)
{
    size_t buflen = 0;

    if(argc)
    {
        if(PARSING_SYSTEM_ALLOWED)
        {
            ityp tmp;
            if(!parse(argv[0], &tmp))
            {
                printUsage(&syslog_manager[SYSLOG_SETBUFLEN]);
                return;
            }
            buflen = tmp;
        }
        else
            buflen = strtoul(argv[0], NULL, sizeof(access(sysLog)->buflen));
    }
    setLogBufLen(access(sysLogPath), access(sysLog), buflen);
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system flushLogBuf(const register sel_typ argc, char ** argv)
{
    sprint("\nSystem Log:\n%s\nBuffer has been properly flushed.\n\n", access(sysLogPath));
    _flushLogBuf(access(sysLog));
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __system renLog(const register sel_typ argc, char ** argv)
{
    char name[MAX_PATH_LENGTH];
    char newname[MAX_PATH_LENGTH];

    CLEARBUFFER();

    strcpy(name, access(sysLogPath));

    size_t len;

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
        printf2(COLOR_CREDITS, "\nEnter System Log newname.\n");
        printf2(COLOR_CREDITS, "Enter %c to exit SubProgram.\n\n", SCANFEXIT_CHAR);

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
        sprint("System Log:\n%s\nhas been correctly renamed to:\n%s.\n\n", name, newname);
    }

    return;
}

#endif
