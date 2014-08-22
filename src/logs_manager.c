// logs_manager.c 20/08/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifdef ALLOW_LOGSMANAGER

// FUNCTIONS DECLARATIONS

__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog setCurLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __lmp_prog createLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog viewLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog printLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog updLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog updAllLogs(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog delLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog relLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __lmp_prog editLog(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog selLogBufLen(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog flushLogBuf(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog renLog(const register sel_typ argc, char ** argv);

sprog logs_manager[MAX_LOGSMANAGER_PROGS] =
{
    [LOGS_SETCURRENT] =
    {
        "Select Current Log",
        CMD_SETCURLOG,
        USAGE_SETCURLOG,
        setCurLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_OPEN] =
    {
        "Open Log ("DEFAULT_LOG_FILE_EXTENSION")",
        CMD_OPENLOG,
        USAGE_OPENLOG,
        createLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_CREATE] =
    {
        "Create new Log ("DEFAULT_LOG_FILE_EXTENSION")",
        CMD_CREATELOG,
        USAGE_CREATELOG,
        createLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_READ] =
    {
        "View Log",
        CMD_READLOG,
        USAGE_READLOG,
        viewLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_PRINT] =
    {
        "Print Log",
        CMD_PRTLOG,
        USAGE_PRTLOG,
        printLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_UPDATE] =
    {
        "Save Log",
        CMD_UPDLOG,
        USAGE_UPDLOG,
        updLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_UPDATEALL] =
    {
        "Save all Logs",
        CMD_UPDALLLOGS,
        USAGE_UPDALLLOGS,
        updAllLogs,
        AUTOMATIC,
        CHILD
    },
    [LOGS_SAVE] =
    {
        "Save Log as ("DEFAULT_LOG_FILE_EXTENSION")",
        CMD_SAVELOG,
        USAGE_SAVELOG,
        createLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_DELETE] =
    {
        "Delete Log",
        CMD_DELLOG,
        USAGE_DELLOG,
        delLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_DELETEPHYSICAL] =
    {
        "Delete Physical Log",
        CMD_DELLOG2,
        USAGE_DELLOG2,
        delLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_DELETEALL] =
    {
        "Delete all Logs",
        CMD_DELALLLOGS,
        USAGE_DELALLLOGS,
        delLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_DELETEALLPHYSICALS] =
    {
        "Delete all Physical Logs",
        CMD_DELALLLOGS2,
        USAGE_DELALLLOGS2,
        delLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_RELOAD] =
    {
        "Reset Log",
        CMD_RESETLOG,
        USAGE_RESETLOG,
        relLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_EDIT] =
    {
        "Edit Log ("DEFAULT_LOG_FILE_EXTENSION")",
        CMD_EDITLOG,
        USAGE_EDITLOG,
        editLog,
        AUTOMATIC,
        CHILD
    },
    [LOGS_SETBUFLEN] =
    {
        "Select Log Buffer Length",
        CMD_SETLOGBUFLEN,
        USAGE_SETLOGBUFLEN,
        selLogBufLen,
        AUTOMATIC,
        CHILD
    },
    [LOGS_EMPTYBUFFER] =
    {
        "Flush Log Buffer",
        CMD_EMPTYLOGBUF,
        USAGE_EMPTYLOGBUF,
        flushLogBuf,
        AUTOMATIC,
        CHILD
    },
    [LOGS_RENAME] =
    {
        "Rename Log ("DEFAULT_LOG_FILE_EXTENSION")",
        CMD_RENLOG,
        USAGE_RENLOG,
        renLog,
        AUTOMATIC,
        CHILD
    }
};

// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog setCurLog(const register sel_typ argc, char ** argv)
{
    setCurItem(argc ? getItemID(argv[0], &logs_manager[LOGS_SETCURRENT], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static inline void _MS__private __lmp_prog createLog(const register sel_typ argc, char ** argv)
{
    createItem(argc ? argv[0] : NULL, LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog viewLog(const register sel_typ argc, char ** argv)
{
    viewItem(argc ? getItemID(argv[0], &logs_manager[LOGS_READ], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog printLog(const register sel_typ argc, char ** argv)
{
    printListItem(argc ? getItemID(argv[0], &logs_manager[LOGS_PRINT], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog updLog(const register sel_typ argc, char ** argv)
{
    updItem(LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog updAllLogs(const register sel_typ argc, char ** argv)
{
    updAll(LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog delLog(const register sel_typ argc, char ** argv)
{
    delItem(argc ? getItemID(argv[0], &logs_manager[__pmode__], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog relLog(const register sel_typ argc, char ** argv)
{
    relItem(argc ? getItemID(argv[0], &logs_manager[LOGS_RELOAD], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static void _MS__private __lmp_prog editLog(const register sel_typ argc, char ** argv)
{
    dim_typ which_log;

    if((which_log = argc ? getItemID(argv[0], &logs_manager[LOGS_EDIT], LOGS) : itemSelect(LOGS)) == NULL_ITEM(LOGS))
        return;

    char name[MAX_PATH_LENGTH] = NULL_CHAR;
    strcpy(name, listNo(which_log, LOGS)->path);
    _editLog(name);

    return;
}

__MSSHELL_WRAPPER_ void _MS__private __lmp_prog setLogBufLen(const char path[static MAX_PATH_LENGTH], logObj * const which_log, const size_t buflen)
{
    const size_t old_buflen = which_log->buflen;

    ityp tmp;
    tmp = 0.00;

    if(buflen)
    {
        if(buflen != (which_log->buflen = (size_t)tmp) || which_log->buflen < MIN_BUFLEN || which_log->buflen > MAX_BUFLEN)
        {
            printErr(33, "Invalid inserted Dimension Value.\nMust be an integer between %zu and %zu", MIN_BUFLEN, MAX_BUFLEN);
            return;
        }
    }
    else
    {
        printf2(COLOR_CREDITS, "Enter the Length of the Selected Log Buffer:\n%s.\n\n", path);

        if(PARSING_SYSTEM_ALLOWED)
            PRINTHOWTOBACKMESSAGE();

        while((PARSING_SYSTEM_ALLOWED ? ((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)) == NULL_VAL) :
        (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (which_log->buflen = (size_t)tmp) || which_log->buflen < MIN_BUFLEN || which_log->buflen > MAX_BUFLEN)
        {
            CLEARBUFFER();

            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
            if(exitHandleCheck)
            {
                which_log->buflen = old_buflen;
                return;
            }
            printErr(33, "Invalid inserted Dimension Value.\nMust be an integer between %zu and %zu", MIN_BUFLEN, MAX_BUFLEN);
        }
    }

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        *(access(exprVars)->e_ANS) = which_log->buflen;

    printf2(COLOR_USER, "\nIt has been correctly selected %zu as the Selected Log Buffer Length:\n", which_log->buflen);
    printf2(COLOR_CREDITS, path);
    printf2(COLOR_USER, ".\n\n");
    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog selLogBufLen(const register sel_typ argc, char ** argv)
{
    dim_typ which_log;

    if((which_log = argc ? getItemID(argv[0], &logs_manager[LOGS_SETBUFLEN], LOGS) : itemSelect(LOGS)) == NULL_ITEM(LOGS))
        return;

    size_t buflen = 0;

    if(argc > 1)
    {
        ityp tmp;
        if(PARSING_SYSTEM_ALLOWED)
        {
            if(!parse(argv[1], &tmp))
            {
                printUsage(&logs_manager[LOGS_SETBUFLEN]);
                return;
            }
            buflen = tmp;
        }
        else
            buflen = strtoul(argv[1], NULL, sizeof(buflen));
    }

    nodelist * const node = listNo(which_log, LOGS);

    setLogBufLen(node->path, ((logObj *)(node->data)), buflen);

    return;
}

__MSSHELL_WRAPPER_ static void _MS__private __lmp_prog flushLogBuf(const register sel_typ argc, char ** argv)
{
    dim_typ which_log;

    if((which_log = argc ? getItemID(argv[0], &logs_manager[LOGS_EMPTYBUFFER], LOGS) : itemSelect(LOGS)) == NULL_ITEM(LOGS))
        return;

    nodelist * const tmp = listNo(which_log, LOGS);

    printf2(COLOR_USER, "\nIt has been properly emptied the Buffer of the Log:\n");
    printf2(COLOR_CREDITS, tmp->path);
    printf2(COLOR_USER, ".\n\n");
    _flushLogBuf(((logObj *)(tmp->data)));
    return;
}

__MSSHELL_WRAPPER_ static inline void _MS__private __lmp_prog renLog(const register sel_typ argc, char ** argv)
{
    renItem(argc > 1 ? argv[1] : NULL, argc ? getItemID(argv[0], &logs_manager[LOGS_RENAME], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}
#endif
