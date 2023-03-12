// logs_manager.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h"
#ifndef __DISABLE_LOGSMANAGER

// FUNCTIONS DECLARATIONS

__MSSHELL_WRAPPER_ static void  __lmp_prog setCurLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void  __lmp_prog createLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog viewLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog printLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog updLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog updAllLogs(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog delLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog relLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __WINCALL static void  __lmp_prog editLog(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog selLogBufLen(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog flushLogBuf(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  __lmp_prog renLog(const sel_typ argc, char ** argv);

sprog logs_manager[MAX_LOGSMANAGER_PROGS] =
{
    [LOGS_SETCURRENT] =
    {
    	CMD_SETCURLOG,
        NAME_SETCURLOG,
        USAGE_SETCURLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_SETCURLOG,
        #endif
        setCurLog,
        ARGC_SETCURLOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_OPEN] =
    {
    	CMD_OPENLOG,
        NAME_OPENLOG,
        USAGE_OPENLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_OPENLOG,
        #endif
        createLog,
        ARGC_OPENLOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_CREATE] =
    {
    	CMD_CREATELOG,
        NAME_CREATELOG,
        USAGE_CREATELOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_CREATELOG,
        #endif
        createLog,
        ARGC_CREATELOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_READ] =
    {
    	CMD_READLOG,
        NAME_READLOG,
        USAGE_READLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_READLOG,
        #endif
        viewLog,
        ARGC_READLOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_PRINT] =
    {
    	CMD_PRTLOG,
        NAME_PRTLOG,
        USAGE_PRTLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_PRTLOG,
        #endif
        printLog,
        ARGC_PRTLOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_UPDATE] =
    {
    	CMD_UPDLOG,
        NAME_UPDLOG,
        USAGE_UPDLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDLOG,
        #endif
        updLog,
        ARGC_UPDLOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_UPDATEALL] =
    {
    	CMD_UPDALLLOGS,
        NAME_UPDALLLOGS,
        USAGE_UPDALLLOGS,
        #ifndef __DISABLE_DATABASE
        LEVEL_UPDALLLOGS,
        #endif
        updAllLogs,
        ARGC_UPDALLLOGS,
        AUTOMATIC,
        CHILD
    },
    [LOGS_SAVE] =
    {
    	CMD_SAVELOG,
        NAME_SAVELOG,
        USAGE_SAVELOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_SAVELOG,
        #endif
        createLog,
        ARGC_SAVELOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_DELETE] =
    {
    	CMD_DELLOG,
        NAME_DELLOG,
        USAGE_DELLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELLOG,
        #endif
        delLog,
        ARGC_DELLOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_DELETEPHYSICAL] =
    {
    	CMD_DELLOG2,
        NAME_DELLOG2,
        USAGE_DELLOG2,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELLOG2,
        #endif
        delLog,
        ARGC_DELLOG2,
        AUTOMATIC,
        CHILD
    },
    [LOGS_DELETEALL] =
    {
    	CMD_DELALLLOGS,
        NAME_DELALLLOGS,
        USAGE_DELALLLOGS,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELALLLOGS,
        #endif
        delLog,
        ARGC_DELALLLOGS,
        AUTOMATIC,
        CHILD
    },
    [LOGS_DELETEALLPHYSICALS] =
    {
    	CMD_DELALLLOGS2,
        NAME_DELALLLOGS2,
        USAGE_DELALLLOGS2,
        #ifndef __DISABLE_DATABASE
        LEVEL_DELALLLOGS2,
        #endif
        delLog,
        ARGC_DELALLLOGS2,
        AUTOMATIC,
        CHILD
    },
    [LOGS_RELOAD] =
    {
    	CMD_RESETLOG,
        NAME_RESETLOG,
        USAGE_RESETLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_RESETLOG,
        #endif
        relLog,
        ARGC_RESETLOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_EDIT] =
    {
    	CMD_EDITLOG,
        NAME_EDITLOG,
        USAGE_EDITLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_EDITLOG,
        #endif
        editLog,
        ARGC_EDITLOG,
        AUTOMATIC,
        CHILD
    },
    [LOGS_SETBUFLEN] =
    {
    	CMD_SETLOGBUFLEN,
        NAME_SETLOGBUFLEN,
        USAGE_SETLOGBUFLEN,
        #ifndef __DISABLE_DATABASE
        LEVEL_SETLOGBUFLEN,
        #endif
        selLogBufLen,
        ARGC_SETLOGBUFLEN,
        AUTOMATIC,
        CHILD
    },
    [LOGS_EMPTYBUFFER] =
    {
    	CMD_EMPTYLOGBUF,
        NAME_EMPTYLOGBUF,
        USAGE_EMPTYLOGBUF,
        #ifndef __DISABLE_DATABASE
        LEVEL_EMPTYLOGBUF,
        #endif
        flushLogBuf,
        ARGC_EMPTYLOGBUF,
        AUTOMATIC,
        CHILD
    },
    [LOGS_RENAME] =
    {
    	CMD_RENLOG,
        NAME_RENLOG,
        USAGE_RENLOG,
        #ifndef __DISABLE_DATABASE
        LEVEL_RENLOG,
        #endif
        renLog,
        ARGC_RENLOG,
        AUTOMATIC,
        CHILD
    }
};

// FUNCTIONS DEFINITIONS

__MSSHELL_WRAPPER_ static inline void  __lmp_prog setCurLog(const sel_typ argc, char ** argv)
{
    setCurItem(argc ? getItemID(argv[0], &logs_manager[LOGS_SETCURRENT], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static inline void  __lmp_prog createLog(const sel_typ argc, char ** argv)
{
    createItem(argc ? argv[0] : NULL, LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog viewLog(const sel_typ argc, char ** argv)
{
    viewItem(argc ? getItemID(argv[0], &logs_manager[LOGS_READ], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog printLog(const sel_typ argc, char ** argv)
{
    printListItem(argc ? getItemID(argv[0], &logs_manager[LOGS_PRINT], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog updLog(const sel_typ argc, char ** argv)
{
    updItem(LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog updAllLogs(const sel_typ argc, char ** argv)
{
    updAll(LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog delLog(const sel_typ argc, char ** argv)
{
    delItem(argc ? getItemID(argv[0], &logs_manager[__pmode__], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog relLog(const sel_typ argc, char ** argv)
{
    relItem(argc ? getItemID(argv[0], &logs_manager[LOGS_RELOAD], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}

__MSSHELL_WRAPPER_ __WINCALL static void  __lmp_prog editLog(const sel_typ argc, char ** argv)
{
    dim_typ which_log;

    if((which_log = argc ? getItemID(argv[0], &logs_manager[LOGS_EDIT], LOGS) : itemSelect(LOGS)) == NULL_ITEM(LOGS))
        return;

    char name[MAX_PATH_LENGTH] = NULL_CHAR;
    strcpy(name, listNo(which_log, LOGS)->path);
    _editLog(name);

    return;
}

__MSSHELL_WRAPPER_ void  __lmp_prog setLogBufLen(const char path[static MAX_PATH_LENGTH], logObj * const which_log, const size_t buflen)
{
    const size_t old_buflen = which_log->buflen;

    mpfr_t tmp;

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
        msprintf(COLOR_CREDITS, "Enter the Length of the Selected Log Buffer:\n%s.\n\n", path);
        PRINTHOWTOBACKMESSAGE();
        
        mpfr_t tmp;

        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (which_log->buflen = mpfr_get_ui(tmp, MPFR_RNDN))) || which_log->buflen < MIN_BUFLEN || which_log->buflen > MAX_BUFLEN)
        {
            CLEARBUFFER();

            if(exitHandleCheck)
            {
                which_log->buflen = old_buflen;
                mpfr_clear(tmp);
                return;
            }
            printErr(33, "Invalid inserted Dimension Value.\nMust be an integer between %zu and %zu", MIN_BUFLEN, MAX_BUFLEN);
        }
    }

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
        mpfr_set_ui(*(access(exprVars)->e_ANS), which_log->buflen, MPFR_RNDN); 

    msprintf(COLOR_USER, "\nIt has been correctly selected %zu as the Selected Log Buffer Length:\n", which_log->buflen);
    msprintf(COLOR_CREDITS, path);
    msprintf(COLOR_USER, ".\n\n");
    mpfr_clear(tmp);
    return;
}

__MSSHELL_WRAPPER_ static void  __lmp_prog selLogBufLen(const sel_typ argc, char ** argv)
{
    dim_typ which_log;

    if((which_log = argc ? getItemID(argv[0], &logs_manager[LOGS_SETBUFLEN], LOGS) : itemSelect(LOGS)) == NULL_ITEM(LOGS))
        return;

    size_t buflen = 0;

    if(argc > 1)
    {
        mpfr_t tmp;
        if(!parse(argv[1], &tmp))
        {
            printUsage(&logs_manager[LOGS_SETBUFLEN]);
            return;
        }
        buflen = mpfr_get_ui(tmp, MPFR_RNDN);
    }

    nodelist * const node = listNo(which_log, LOGS);
    setLogBufLen(node->path, ((logObj *)(node->data)), buflen);
    return;
}

__MSSHELL_WRAPPER_ static void  __lmp_prog flushLogBuf(const sel_typ argc, char ** argv)
{
    dim_typ which_log;

    if((which_log = argc ? getItemID(argv[0], &logs_manager[LOGS_EMPTYBUFFER], LOGS) : itemSelect(LOGS)) == NULL_ITEM(LOGS))
        return;

    nodelist * const tmp = listNo(which_log, LOGS);

    msprintf(COLOR_USER, "\nIt has been properly emptied the Buffer of the Log:\n");
    msprintf(COLOR_CREDITS, tmp->path);
    msprintf(COLOR_USER, ".\n\n");
    _flushLogBuf(((logObj *)(tmp->data)));
    return;
}

__MSSHELL_WRAPPER_ static inline void  __lmp_prog renLog(const sel_typ argc, char ** argv)
{
    renItem(argc > 1 ? argv[1] : NULL, argc ? getItemID(argv[0], &logs_manager[LOGS_RENAME], LOGS) : getItemsListNo(LOGS), LOGS);
    return;
}
#endif
