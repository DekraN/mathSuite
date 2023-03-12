#include "dutils.h"

#ifndef __DISABLE_SYSTEM

__MSSHELL_WRAPPER_ static void scriptingMode(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void exitProgram(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void doCommit(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void execMutant(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void doCheckOut(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void cleanMutantFiles(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void doPurge(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void purgeDeep(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void purgeFast(const sel_typ argc, char ** argv);

sprog system_prog[MAX_SYSTEM_PROGS] =
{
	[SYSTEM_SCRIPTINGMODE] =
	{
		CMD_SCRIPTINGMODE,
		NAME_SCRIPTINGMODE,
		USAGE_SCRIPTINGMODE,
		#ifndef __DISABLE_DATABASE
		LEVEL_SCRIPTINGMODE,
		#endif
		scriptingMode,
		ARGC_SCRIPTINGMODE,
		AUTOMATIC,
		CHILD
	},
	[SYSTEM_EXIT] =
	{
		CMD_EXIT,
		NAME_EXIT,
		USAGE_EXIT,
		#ifndef __DISABLE_DATABASE
		LEVEL_EXIT,
		#endif
		exitProgram,
		ARGC_EXIT,
		AUTOMATIC,
		CHILD
	}, 
	[SYSTEM_MUTANT] =
	{
		CMD_MUTANT,
		NAME_MUTANT,
		USAGE_MUTANT,
		#ifndef __DISABLE_DATABASE
		LEVEL_MUTANT,
		#endif
		execMutant,
		ARGC_MUTANT,
		AUTOMATIC,
		CHILD
	},
	[SYSTEM_COMMIT] =
	{
		CMD_COMMIT,
		NAME_COMMIT,
		USAGE_COMMIT,
		#ifndef __DISABLE_DATABASE
		LEVEL_COMMIT,
		#endif
		doCommit,
		ARGC_COMMIT,
		AUTOMATIC,
		CHILD
	},
	[SYSTEM_CHECKOUT] =
	{
		CMD_CHECKOUT,
		NAME_CHECKOUT,
		USAGE_CHECKOUT,
		#ifndef __DISABLE_DATABASE
		LEVEL_CHECKOUT,
		#endif
		doCheckOut,
		ARGC_CHECKOUT,
		AUTOMATIC,
		CHILD
	},
	[SYSTEM_CLEAN] =
	{
		CMD_CLEAN,
		NAME_CLEAN,
		USAGE_CLEAN,
		#ifndef __DISABLE_DATABASE
		LEVEL_CLEAN,
		#endif
		cleanMutantFiles,
		ARGC_MUTANT,
		AUTOMATIC,
		CHILD
	},
	[SYSTEM_PURGE] =
	{
		CMD_PURGE,
		NAME_PURGE,
		USAGE_PURGE,
		#ifndef __DISABLE_DATABASE
		LEVEL_PURGE,
		#endif
		doPurge,
		ARGC_PURGE,
		AUTOMATIC,
		CHILD
	},
	[SYSTEM_PURGEDEEP] =
	{
		CMD_PURGEDEEP,
		NAME_PURGEDEEP,
		USAGE_PURGEDEEP,
		#ifndef __DISABLE_DATABASE
		LEVEL_PURGEDEEP,
		#endif
		purgeDeep,
		ARGC_PURGEDEEP,
		AUTOMATIC,
		CHILD
	},
	[SYSTEM_PURGEFAST] =
	{
		CMD_PURGEFAST,
		NAME_PURGEFAST,
		USAGE_PURGEFAST,
		#ifndef __DISABLE_DATABASE
		LEVEL_PURGEFAST,
		#endif
		purgeFast,
		ARGC_PURGEFAST,
		AUTOMATIC,
		CHILD
	}
};

__MSSHELL_WRAPPER_ static void scriptingMode(const sel_typ argc, char ** argv)
{
	
	access(mss) = true;
    msprintf(COLOR_USER, "\nScripting Mode has been enabled.\n\n");
	return;
}


__MSSHELL_WRAPPER_ static void exitProgram(const sel_typ argc, char ** argv)
{
	PRINTN();
    register int excode = EXIT_SUCCESS;
    
    if(argc)
    {
    	mpfr_t tmp;
    	if(!parse(argv[0], &tmp))
	    {
	    	mpfr_clear(tmp);
	        printErr(1, "Parse Error on "CMD_EXIT" command.");
	        return;
	    }
	
	    if(mpfr_cmp_si(tmp, INT_MIN) < 0 || mpfr_cmp_si(tmp, INT_MAX) > 0)
	    {
	    	mpfr_clear(tmp);
	        printErr(33, CMD_EXIT" accepts only integers between %d and %d", INT_MIN, INT_MAX);
	        printUsage(&system_prog[SYSTEM_EXIT]);
	        return;
	    }
	    excode = mpfr_get_si(tmp, MPFR_RNDN);
	    mpfr_clear(tmp);
    }
    
    exit(excode);
	return;
}

__MSSHELL_WRAPPER_ static void execMutant(const sel_typ argc, char ** argv)
{
	msprintf(COLOR_USER, MUTANT_WELCOME);
	system(MUTANT_FILENAME);
	if(isSett(BOOLS_EXITAFTERMUTCODEXEC))
		exit(EXIT_MUTANTCODE);
	return;
}

__MSSHELL_WRAPPER_ static void doCommit(const sel_typ argc, char ** argv)
{
	system(_COMMIT_CMD);
	msyprintf(COLOR_SYSTEM, "Successfully performed \""CMD_COMMIT"\" command.");
	return;
}


__MSSHELL_WRAPPER_ static void doCheckOut(const sel_typ argc, char ** argv)
{
	if(!file_exists(MUTANT_DIR"\\"MUTANTOBJ_ARCHIVE))
	{
		printErr(2, MUTANTOBJ_ARCHIVE" doesn't exist");
		return;
	}
	
	system(_CHECKOUT_CMD);
	msyprintf(COLOR_SYSTEM, "Successfully performed \""CMD_CHECKOUT"\" command.");
	return;
}

__MSSHELL_WRAPPER_ static void cleanMutantFiles(const sel_typ argc, char ** argv)
{
	system(MAKE_EXECUTABLE" clean");
	msyprintf(COLOR_SYSTEM, "Successfully performed \""CMD_CLEAN"\" command.");
	return;
}

__MSSHELL_WRAPPER_ static void doPurge(const sel_typ argc, char ** argv)
{
	removeCriticalFiles(RCF_EXPREVAL | RCF_SKIP);
	msyprintf(COLOR_SYSTEM, "Successfully performed \""CMD_PURGE"\" command.");
	return;
}

__MSSHELL_WRAPPER_ static void purgeDeep(const sel_typ argc, char ** argv)
{
	removeCriticalFiles(RCF_DEEP | RCF_SKIP);
	msyprintf(COLOR_SYSTEM, "Successfully performed \""CMD_PURGEDEEP"\" command.");
	return;
}

__MSSHELL_WRAPPER_ static void purgeFast(const sel_typ argc, char ** argv)
{
	removeCriticalFiles(RCF_SKIP);
	msyprintf(COLOR_SYSTEM, "Successfully performed \""CMD_PURGEFAST"\" command.");
	return;
}


#endif
