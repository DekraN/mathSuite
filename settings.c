// settings.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#ifndef FREEZES_SETTINGSMANAGER

#include "dutils.h"

// Settings Function Declarations and Definitions
//

#ifndef __DISABLE_SYSTEM
	__MSSHELL_WRAPPER_ static void   __apnt sysManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   __apnt sysManager(const sel_typ argc, char ** argv)
	{
	
	    operationsGroupMenu(MAX_SYSTEM_PROGS, system_prog, change_settings[SETTINGS_SYSMANAGER].name, BY_CHARS); // MAX_ENVSMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
	    return;
	}
#endif

#ifndef __DISABLE_VARLISTMANAGER
	__MSSHELL_WRAPPER_ static void   __apnt envsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   __apnt envsManager(const sel_typ argc, char ** argv)
	{
	
	    operationsGroupMenu(MAX_ENVSMANAGER_PROGS, envs_manager, change_settings[SETTINGS_ENVSMANAGER].name, BY_CHARS); // MAX_ENVSMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
	    return;
	}
#endif

#ifndef __DISABLE_LOGSMANAGER
	__MSSHELL_WRAPPER_ static void   __apnt logsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   __apnt logsManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_LOGSMANAGER_PROGS, logs_manager, change_settings[SETTINGS_LOGSMANAGER].name, BY_CHARS); // MAX_LOGSMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
	    return;
	}
#endif

#ifndef __DISABLE_SYSLOGMANAGER
	__MSSHELL_WRAPPER_ static void   __apnt sysLogManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   __apnt sysLogManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_SYSLOGMANAGER_PROGS, syslog_manager, change_settings[SETTINGS_SYSLOGMANAGER].name, BY_CHARS); // MAX_SYSLOGMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
	    return;
	}
#endif

#ifndef __DISABLE_LAYOUTSMANAGER
	__MSSHELL_WRAPPER_ static void   __apnt layoutsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   __apnt layoutsManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_LAYOUTSMANAGER_PROGS, layouts_manager, change_settings[SETTINGS_LAYOUTSMANAGER].name, BY_CHARS); // MAX_LAYOUTSMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
	    return;
	}
#endif

#ifndef __DISABLE_COLSMANAGER
	__MSSHELL_WRAPPER_ static void   __apnt colorsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   __apnt colorsManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_COLSMANAGER_PROGS, cols_manager, change_settings[SETTINGS_COLORSMANAGER].name, BY_CHARS); // MAX_COLSMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
	    return;
	}
#endif
	
#ifndef __DISABLE_LFSMANAGER
	__MSSHELL_WRAPPER_ static void   __apnt lfsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   __apnt lfsManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_LFSMANAGER_PROGS, lfs_manager, change_settings[SETTINGS_LFSMANAGER].name, BY_CHARS); // MAX_LFSMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
	    return;
	}
#endif

#ifndef __DISABLE_PRECEDIT
	__MSSHELL_WRAPPER_ static void   changePrecision(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changePrecision(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_precision = access(curLayout)->precision;
	    mpfr_t tmp;
	
	    if(argc)
	    {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->precision = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->precision < MIN_PRECISION || access(curLayout)->precision > MAX_PRECISION)
            {
            	mpfr_clear(tmp);
            	access(curLayout)->precision = old_precision;
                printErr(33, "Invalid inserted Precision Value.\nMust be an integer between %hu and %hu", MIN_PRECISION, MAX_PRECISION);
                printUsage(&change_settings[SETTINGS_CHANGEPRECISION]);
                return;
            }
	    }
	    else
	    {
	        msprintf(COLOR_CREDITS, "\nEnter a Value as Program Precision: between %hu and %hu.\n", MIN_PRECISION, MAX_PRECISION);
	        PRINTHOWTOBACKMESSAGE();
	        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (access(curLayout)->precision = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->precision < MIN_PRECISION || access(curLayout)->precision > MAX_PRECISION)
	        {
	    		mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
	            {
	                access(curLayout)->precision = old_precision;
	                return;
	            }
	            printErr(33, "Invalid inserted Precision Value.\nMust be an integer between %hu and %hu", MIN_PRECISION, MAX_PRECISION);
	        }
	    }
	
		mpfr_clear(tmp);
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        mpfr_set_ui(*(access(exprVars)->e_ANS), access(curLayout)->precision, MPFR_RNDN);
	
	    msyprintf(COLOR_SYSTEM, "Program Precision Value has been correctly changed from: %hu to: %hu.\n", old_precision, access(curLayout)->precision);
	    return;
	}
#endif

#ifndef __DISABLE_STABFACTEDIT
	__MSSHELL_WRAPPER_ static void   changeStabilizerFactor(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeStabilizerFactor(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_stabilizer_factor = access(curLayout)->stabilizer_factor;
	    mpfr_t tmp;
	
	    if(argc)
	    {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->stabilizer_factor = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->stabilizer_factor < MIN_STABFACT || access(curLayout)->stabilizer_factor > MAX_STABFACT)
            {
            	mpfr_clear(tmp);
            	access(curLayout)->stabilizer_factor = old_stabilizer_factor;
                printErr(33, "Invalid inserted Stabilizer Factor Value.\nMust be an integer between %hu and %hu", MIN_STABFACT, MAX_STABFACT);
                printUsage(&change_settings[SETTINGS_CHANGESTABILIZERFACTOR]);
                return;
            }
	    }
	    else
	    {
	        msprintf(COLOR_CREDITS, "\nEnter a Value as Program Stabilizer Factor: a non-negative and non-zero integer.\n");
	        PRINTHOWTOBACKMESSAGE();
	        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (access(curLayout)->stabilizer_factor = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->stabilizer_factor < MIN_STABFACT || access(curLayout)->stabilizer_factor > MAX_STABFACT)
	        {
	        	mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
	            {
	                access(curLayout)->stabilizer_factor = old_stabilizer_factor;
	                return;
	            }
	            printErr(33, "Invalid inserted Stabilizer Factor Value.\nMust be an integer between %hu and %hu", MIN_STABFACT, MAX_STABFACT);
	        }
	    }
	    
	    mpfr_clear(tmp);
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        mpfr_set_ui(*(access(exprVars)->e_ANS), access(curLayout)->stabilizer_factor, MPFR_RNDN);
	
	    msyprintf(COLOR_SYSTEM, "Program Stabilizer Factor Value has been correctly changed from: %hu to: %hu.\n", old_stabilizer_factor, access(curLayout)->stabilizer_factor);
	    return;
	}
#endif

#ifndef __DISABLE_BLOCKSIZEEDIT
	__MSSHELL_WRAPPER_ static void   changeBlockSize(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeBlockSize(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_block_size= access(curLayout)->block_size;
	    mpfr_t tmp;
	
	    if(argc)
	    {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->block_size = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->block_size < MIN_BLOCKSIZE)
            {
            	mpfr_clear(tmp);
            	access(curLayout)->block_size = old_block_size;
                printErr(33, "Invalid inserted Block Size Value.\nMust be an integer greater or equal than %hu", MIN_BLOCKSIZE);
                printUsage(&change_settings[SETTINGS_CHANGEBLOCKSIZE]);
                return;
            }
	    }
	    else
	    {
	        msprintf(COLOR_CREDITS, "\nEnter a Value as Program Block Size: a non-negative and non-zero integer.\n");
	        PRINTHOWTOBACKMESSAGE();
	        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (access(curLayout)->block_size = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->block_size < MIN_BLOCKSIZE)
	        {
	        	mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
	            {
	                access(curLayout)->block_size = old_block_size;
	                return;
	            }
	            printErr(33, "Invalid inserted Block Size Value.\nMust be an integer greater or equal than %hu", MIN_BLOCKSIZE);
	        }
	    }
	
		mpfr_clear(tmp);
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        mpfr_set_ui(*(access(exprVars)->e_ANS), access(curLayout)->block_size, MPFR_RNDN);
	
	    msyprintf(COLOR_SYSTEM, "Program Block Size Value has been correctly changed from: %hu to: %hu.\n", old_block_size, access(curLayout)->block_size);
	    return;
	}
#endif

#ifndef __DISABLE_MINOSMMDIMEDIT
	__MSSHELL_WRAPPER_ static void   changeMinOSMMDim(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeMinOSMMDim(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_min_osmm_dim= access(curLayout)->min_osmm_dim;
	    mpfr_t tmp;
	
	    if(argc)
	    {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->min_osmm_dim = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->min_osmm_dim < MIN_OSMM_DIM)
            {
            	mpfr_clear(tmp);
            	access(curLayout)->min_osmm_dim = old_min_osmm_dim;
                printErr(33, "Invalid inserted Min OSMM Dimension Value.\nMust be an integer greater or equal than %hu", MIN_OSMM_DIM);
                printUsage(&change_settings[SETTINGS_CHANGEMINOSMMDIM]);
                return;
            }
	    }
	    else
	    {
	        msprintf(COLOR_CREDITS, "\nEnter a Value as Min OSMM Dimension: a non-negative and non-zero integer.\n");
	        PRINTHOWTOBACKMESSAGE();
	        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (access(curLayout)->min_osmm_dim = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->min_osmm_dim < MIN_OSMM_DIM)
	        {
	        	mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
	            {
	                access(curLayout)->min_osmm_dim = old_min_osmm_dim;
	                return;
	            }
	            printErr(33, "Invalid inserted Min OSMM Dimension Value.\nMust be an integer greater or equal than %hu", MIN_OSMM_DIM);
	        }
	    }
	    
		mpfr_clear(tmp);
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        mpfr_set_ui(*(access(exprVars)->e_ANS), access(curLayout)->min_osmm_dim, MPFR_RNDN);
	
	    msyprintf(COLOR_SYSTEM, "Program Min OSMM Dimension Value has been correctly changed from: %hu to: %hu.\n", old_min_osmm_dim, access(curLayout)->min_osmm_dim);
	    return;
	}
#endif

#ifndef __DISABLE_MINSTRASSENDIMEDIT
	__MSSHELL_WRAPPER_ static void   changeMinStrassenDim(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeMinStrassenDim(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_min_strassen_dim= access(curLayout)->min_strassen_dim;
		mpfr_t tmp;
	
	    if(argc)
	    {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->min_strassen_dim = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->min_strassen_dim < MIN_STRASSEN_DIM)
            {
            	mpfr_clear(tmp);
            	access(curLayout)->min_strassen_dim = old_min_strassen_dim;
                printErr(33, "Invalid inserted Min Strassen Dimension Value.\nMust be an integer greater or equal than %hu", MIN_STRASSEN_DIM);
                printUsage(&change_settings[SETTINGS_CHANGEMINSTRASSENDIM]);
                return;
            }
	    }
	    else
	    {
	        msprintf(COLOR_CREDITS, "\nEnter a Value as Min Strassen Dimension: a non-negative and non-zero integer.\n");
	        PRINTHOWTOBACKMESSAGE();

	        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (access(curLayout)->min_strassen_dim = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->min_strassen_dim < MIN_STRASSEN_DIM)
	        {
	        	mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
	            {
	                access(curLayout)->min_strassen_dim = old_min_strassen_dim;
	                return;
	            }
	            printErr(33, "Invalid inserted Min Strassen Dimension Value.\nMust be an integer greater or equal than %hu", MIN_STRASSEN_DIM);
	        }
		}
		
		mpfr_clear(tmp);
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        mpfr_set_ui(*(access(exprVars)->e_ANS), access(curLayout)->min_strassen_dim, MPFR_RNDN);
	
	    msyprintf(COLOR_SYSTEM, "Program Min Strassen Dimension Value has been correctly changed from: %hu to: %hu.\n", old_min_strassen_dim, access(curLayout)->min_strassen_dim);
	    return;
	}
#endif

#ifndef __DISABLE_MINSRNUMBEREDIT
	__MSSHELL_WRAPPER_ static void   changeMinSRNumber(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeMinSRNumber(const sel_typ argc, char ** argv)
	{
		const dim_typ old_minsr_number = access(curLayout)->min_stirling_number;
	    mpfr_t tmp;
	
	    if(argc)
	    {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->min_stirling_number = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->min_stirling_number < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->min_stirling_number > MIN_STIRLING_NUMBER)
            {
            	mpfr_clear(tmp);
            	access(curLayout)->min_stirling_number = old_minsr_number;
                printErr(33, "Invalid inserted Min Stirling Number.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MIN_STIRLING_NUMBER);
                printUsage(&change_settings[SETTINGS_CHANGEMINSRNUMBER]);
                return;
            }
	    }
	    else
	    {
	        msprintf(COLOR_CREDITS, "\nEnter a Value as Min Stirling Number: a non-negative and non-zero integer.\n");
	        PRINTHOWTOBACKMESSAGE();
	        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (access(curLayout)->min_stirling_number = mpfr_get_ui(tmp, MPFR_RNDN))) || access(curLayout)->min_stirling_number < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->min_stirling_number > MIN_STIRLING_NUMBER)
	        {
	        	mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
	            {
	                access(curLayout)->min_stirling_number = old_minsr_number;
	                return;
	            }
	            printErr(33, "Invalid inserted Min Stirling Number.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MIN_STIRLING_NUMBER);
	        }
	    }
	    
		mpfr_clear(tmp);
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        mpfr_set_ui(*(access(exprVars)->e_ANS), access(curLayout)->min_stirling_number, MPFR_RNDN);
	
	    msyprintf(COLOR_SYSTEM, "Program Min Stirling Number has been correctly changed from: %hu to: %hu.\n", old_minsr_number, access(curLayout)->min_stirling_number);
		return;
	}

#endif

#ifndef __DISABLE_ALGEBRAEDIT
	__MSSHELL_WRAPPER_ static void   changeAlgebra(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeAlgebra(const sel_typ argc, char ** argv)
	{
	
	    dim_typ i;
	   	 
	    if(argc)
	    {
	        mpfr_t tmp;
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (i = mpfr_get_ui(tmp, MPFR_RNDN))) || i < MIN_ALGEBRA || i > MAX_ALGEBRA)
            {
            	mpfr_clear(tmp);
                printErr(1, "Invalid inserted Value: not correspondent to any Algebra Identifier");
                printUsage(&change_settings[SETTINGS_CHANGEALGEBRA]);
                return;
            }
        	mpfr_clear(tmp);
	    }
	    else if((i = selectListItem(_MAX_ALGEBRA, _MAX_ALGEBRA > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
	            "Select Algebra Identifier in which to perform Algebra Operations", suite_c.algebra_elements_names)) == _MAX_ALGEBRA) return;

	    access(curLayout)->algebra = i;
	    msyprintf(COLOR_SYSTEM, "%s Algebra has been correctly selected.\n\n", suite_c.algebra_elements_names[i]);
	    return;
	}
#endif

#ifndef __DISABLE_OUTLIERCONSTEDIT
	__MSSHELL_WRAPPER_ static void   changeOutlierConst(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeOutlierConst(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_outlier_constant = access(curLayout)->outlier_constant;
		mpfr_t tmp;
		 
	    if(argc)
	    {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_d(tmp, MIN_OUTLIER_CONSTANT) < 0 || mpfr_cmp_d(tmp, MAX_OUTLIER_CONSTANT) > 0)
            {
            	mpfr_clear(tmp);
            	access(curLayout)->outlier_constant = old_outlier_constant;
                printErr(33, "Invalid inserted Outlier Constant.\nMust be a float between %.*f and %.*f", DEFAULT_PRECISION, MIN_OUTLIER_CONSTANT, DEFAULT_PRECISION, MAX_OUTLIER_CONSTANT);
                printUsage(&change_settings[SETTINGS_CHANGEOUTLIERCONST]);
                return;
            }
            access(curLayout)->outlier_constant = mpfr_get_flt(tmp, MPFR_RNDN); 
	    }
	    else
	    {
	        msprintf(COLOR_CREDITS, "\nEnter a Value as Program Outlier Constant: a non-negative and non-zero float.\n");
	        PRINTHOWTOBACKMESSAGE();
	        
	        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_d(tmp, MIN_OUTLIER_CONSTANT) < 0 || mpfr_cmp_d(tmp, MAX_OUTLIER_CONSTANT) > 0)
	        {
	        	mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
	            {
	                access(curLayout)->outlier_constant = old_outlier_constant;
	                return;
	            }
	            printErr(33, "Invalid inserted Outlier Constant.\nMust be a float between %.*f and %.*f", DEFAULT_PRECISION, MIN_OUTLIER_CONSTANT, DEFAULT_PRECISION, MAX_OUTLIER_CONSTANT);
	        }
	    	access(curLayout)->outlier_constant = mpfr_get_flt(tmp, MPFR_RNDN); 
	    }
		
		mpfr_clear(tmp);
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        mpfr_set_ui(*(access(exprVars)->e_ANS), access(curLayout)->outlier_constant, MPFR_RNDN); 
	
	    msyprintf(COLOR_SYSTEM, "Program Outlier Constant has been correctly changed from: %.*f to: %.*f.\n", DEFAULT_PRECISION, old_outlier_constant, DEFAULT_PRECISION, access(curLayout)->outlier_constant);
	    return;
	}
#endif

#ifndef __DISABLE_EXITCHAREDIT
	__MSSHELL_WRAPPER_ static void   changeExitChar(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeExitChar(const sel_typ argc, char ** argv)
	{
	    const char old_exit_char = access(curLayout)->exit_char;
	
	    if(argc)
	    {
	        if(argv[0][0] == 'A' || argv[0][0] == 'B')
	        {
	            printErr(5, "Invalid inserted Char");
	            printUsage(&change_settings[SETTINGS_CHANGEEXITCHAR]);
	            return;
	    	}
	        access(curLayout)->exit_char = argv[0][0];
	    }
	    else
	    {
	        printf("\nDefine a Character for Back to MAIN MENU function.\n\n");
	        while(scanf(" %c", &access(curLayout)->exit_char) != 1 || access(curLayout)->exit_char == 'A' || access(curLayout)->exit_char == 'B')
	            printErr(5, "Invalid inserted Char");
	    }
	
	    CLEARBUFFER();
	    msyprintf(COLOR_SYSTEM, "Back to MAIN MENU Character has been correctly changed from: %c to: %c.\n", old_exit_char, access(curLayout)->exit_char);
	    return;
	}
#endif

#ifndef __DISABLE_RANDOMSEEDEDIT
	__MSSHELL_WRAPPER_ static void   changeRandomSeed(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeRandomSeed(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_random_seed = access(random_seed);
	    mpfr_t tmp;
	
	    if(argc)
	    {
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (access(random_seed) = mpfr_get_ui(tmp, MPFR_RNDN))) || access(random_seed) < MIN_RANDOMSEED || access(random_seed) > MAX_RANDOMSEED)
            {
            	mpfr_clear(tmp);
            	access(random_seed) = old_random_seed;
                printErr(5, "Invalid inserted Value.\nMust be an integer between %hu and %hu", MIN_RANDOMSEED, MAX_RANDOMSEED); // strettamente positivo e minore di %u", UINT_MAX);
                printUsage(&change_settings[SETTINGS_CHANGERANDOMSEED]);
                return;
            }
	    }
	    else
	    {
	        msprintf(COLOR_CREDITS, "\nEnter RandomSeed Value\n");
	        PRINTHOWTOBACKMESSAGE();

	        while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (access(random_seed) = mpfr_get_ui(tmp, MPFR_RNDN))) || access(random_seed) < MIN_RANDOMSEED || access(random_seed) > MAX_RANDOMSEED)
	        {
	        	mpfr_clear(tmp);
	            CLEARBUFFER();
	            if(exitHandleCheck)
	            {
	                access(random_seed) = old_random_seed;
	                return;
	            }
	            printErr(5, "Invalid inserted Value.\nMust be an integer between %hu and %hu", MIN_RANDOMSEED, MAX_RANDOMSEED); // , UINT_MAX);
	        }
	    }
	
		mpfr_clear(tmp);
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        mpfr_set_ui(*(access(exprVars)->e_ANS), access(random_seed), MPFR_RNDN);
	
	    srand(access(random_seed));
	    msyprintf(COLOR_SYSTEM, "Program RandomSeed Value has been correctly changed from: %hu to: %hu.\n", old_random_seed, access(random_seed));
	
	    return;
	}
#endif

#ifndef __DISABLE_BOOLVARSEDIT
	__MSSHELL_WRAPPER_ static void   changeBoolValues(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeBoolValues(const sel_typ argc, char ** argv)
	{
	    dim_typ i;
	    if(argc)
	    {
	        mpfr_t tmp;
            if((!parse(argv[0], &tmp)) || mpfr_cmp_ui(tmp, (i = mpfr_get_ui(tmp, MPFR_RNDN))) || i < 0 || i > MAX_BOOL_SETTINGS)
            {
            	mpfr_clear(tmp);
                printErr(1, "Invalid inserted Value: not correspondent to any Boolean Settings ID");
                printUsage(&change_settings[SETTINGS_CHANGEBOOLVALUES]);
                return;
            }
	    }
	    else if((i = selectListItem(MAX_BOOL_SETTINGS, MAX_BOOL_SETTINGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
	            "Select Bool Settings you wish to change", suite_c.bools_names)) == MAX_BOOL_SETTINGS) return;
	
	    access(curLayout)->bools ^= suite_c.bools[i].bmask;
	    msyprintf(COLOR_SYSTEM, "%s Boolean Settings\nhas been correctly %s.\n\n", suite_c.bools_names[i], isSett(i) ? "ENABLED":"DISABLED");
	    return;
	}
#endif

#ifndef __DISABLE_MMIEDIT
	__MSSHELL_WRAPPER_ static void   changeMaxMemoizableIndices(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   changeMaxMemoizableIndices(const sel_typ argc, char ** argv)
	{
	
	    const dim_typ old_max_memoizable_index[MAX_MEMOIZABLE_FUNCTIONS] =
	    {
	        access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI],
	        //access(curLayout)->max_memoizable_indices[FUNCTION_FACTORIAL],
	        access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL],
	        access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL]
	    };
	
		dim_typ i;
		mpfr_t tmp;
	
	    if(argc)
	    {
	        if(argc > LAST_MEMOIZABLE_FUNCTION)
	        {
                if((!parse(argv[FUNCTION_FIBONACCI], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = mpfr_get_ui(tmp, MPFR_RNDN))) ||
                access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] > MAX_FIBONACCI_MEMOIZABLE_INDEX)
                {
                	mpfr_clear(tmp);
                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
                    printErr(33, "Invalid inserted MIM_FIBO Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_FIBONACCI_MEMOIZABLE_INDEX);
                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
                    return;
                }

				/*
                if((!parse(argv[FUNCTION_FACTORIAL], &tmp)) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_FACTORIAL] = (dim_typ)tmp) ||
                access(curLayout)->max_memoizable_indices[FUNCTION_FACTORIAL] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FACTORIAL] > MAX_FACTORIAL_MEMOIZABLE_INDEX)
                {
                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
                	access(curLayout)->max_memoizable_indices[FUNCTION_FACTORIAL] = old_max_memoizable_index[FUNCTION_FACTORIAL];
                    printErr(33, "Invalid inserted MIM_FACT Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_FACTORIAL_MEMOIZABLE_INDEX);
                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
                    return;
                }
                */
                
                if((!parse(argv[FUNCTION_EVEN_DOUBLEFACTORIAL], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] = mpfr_get_ui(tmp, MPFR_RNDN))) ||
                access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] > MAX_EVEN_DOUBLEFACTORIAL_MEMOIZABLE_INDEX)
                {
                	mpfr_clear(tmp);
                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
                	// access(curLayout)->max_memoizable_indices[FUNCTION_FACTORIAL] = old_max_memoizable_index[FUNCTION_FACTORIAL];
                	access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] = old_max_memoizable_index[FUNCTION_EVEN_DOUBLEFACTORIAL];
                    printErr(33, "Invalid inserted MIM_SFACT_EVEN Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_EVEN_DOUBLEFACTORIAL_MEMOIZABLE_INDEX);
                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
                    return;
                }
                
                if((!parse(argv[FUNCTION_ODD_DOUBLEFACTORIAL], &tmp)) || mpfr_cmp_ui(tmp, (access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] = mpfr_get_ui(tmp, MPFR_RNDN))) ||
                access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] > MAX_ODD_DOUBLEFACTORIAL_MEMOIZABLE_INDEX)
                {
                	mpfr_clear(tmp);
                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
                	// access(curLayout)->max_memoizable_indices[FUNCTION_FACTORIAL] = old_max_memoizable_index[FUNCTION_FACTORIAL];
                	access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] = old_max_memoizable_index[FUNCTION_EVEN_DOUBLEFACTORIAL];
                	access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] = old_max_memoizable_index[FUNCTION_ODD_DOUBLEFACTORIAL];
                    printErr(33, "Invalid inserted MIM_SFACT_ODD Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_ODD_DOUBLEFACTORIAL_MEMOIZABLE_INDEX);
                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
                    return;
                }
	        }
	        else
	        {
	            printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	            return;
	        }
	    }
	    else
	    {
	
			char seperator[SIGN_STRING];
			
			strcpy(seperator, "]\n[");
	        msprintf(COLOR_CREDITS, "\nEnter the three MAX MEMOIZABLE INDICES of the functions respectively\n:%s, %s and %s, as expected format:\n", suite_c.memoizers_names[FUNCTION_FIBONACCI], suite_c.memoizers_names[FUNCTION_EVEN_DOUBLEFACTORIAL], suite_c.memoizers_names[FUNCTION_ODD_DOUBLEFACTORIAL]);
	        msprintf(COLOR_CREDITS, "[MIM_FIBO%sMIM_FACT%sMIM_EVEN_SFACT%sMIM_ODD_SFACT]\n", seperator, seperator, seperator);
	
	        PRINTHOWTOBACKMESSAGE();
	        
	        mpfr_t tmp2, tmp3;
			
	        while((requires(tmp, NULL, NULL_CHAR, "Inserted MAX Fibonacci MEMOIZABLE INDEX is:", PARSER_SHOWRESULT) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = mpfr_get_ui(tmp, MPFR_RNDN))) ||
	            access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] > MAX_MEMOIZABLE_INDEX) ||
			 (requires(tmp2, NULL, NULL_CHAR, "Inserted MAX Even Double Factorial MEMOIZABLE INDEX is:", PARSER_SHOWRESULT) || isNullVal(tmp2) || mpfr_cmp_ui(tmp2, (access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] = mpfr_get_ui(tmp2, MPFR_RNDN))) ||
	            access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] > MAX_MEMOIZABLE_INDEX) ||
				(requires(tmp3, NULL, NULL_CHAR, "Inserted MAX Odd Double Factorial MEMOIZABLE INDEX is:", PARSER_SHOWRESULT) || isNullVal(tmp3) || mpfr_cmp_ui(tmp3, (access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] = mpfr_get_ui(tmp3, MPFR_RNDN))) ||
	            access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] > MAX_MEMOIZABLE_INDEX) || access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] < MIN_MEMOIZABLE_INDEX+1 || 
				 access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] > MAX_MEMOIZABLE_INDEX ||
	              access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL] > MAX_MEMOIZABLE_INDEX || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL] > MAX_MEMOIZABLE_INDEX)
	        {
	        	mpfr_clears(tmp, tmp2, tmp3, NULL); 
	            CLEARBUFFER();
	            if(exitHandleCheck) // if(tmp3[ROWS] == NULL_VAL || tmp3[COLUMNS] == NULL_VAL)
	            {
	            	#pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
	            	for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
	            		access(curLayout)->max_memoizable_indices[i] = old_max_memoizable_index[i];
	                return;
	            }
	            printErr(33, "Invalid [MIM_FIBO MIM_FACT MIM_EVENSFACT MIM_ODDSFACT] format.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_MEMOIZABLE_INDEX);
	        }
	        mpfr_clears(tmp2, tmp3, NULL); 
	    }
	    
	    for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
			msyprintf(COLOR_SYSTEM, "MAX %s MEMOIZABLE INDEX Value has been changed from: %hu to: %hu.\n", old_max_memoizable_index[i], access(curLayout)->max_memoizable_indices[i]);
	
		mpfr_clear(tmp);
	    return;
	}
	
	__MSSHELL_WRAPPER_ static void   emptyMemoizersBuffers(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   emptyMemoizersBuffers(const sel_typ argc, char ** argv)
	{
	
	    sel_typ tmp;
	
	    if(argc)
	    {
	        mpfr_t tmp2;
            if((!parse(argv[0], &tmp2)) || mpfr_cmp_ui(tmp2, (tmp = mpfr_get_ui(tmp2, MPFR_RNDN))) || (tmp != FUNCTION_FIBONACCI && tmp != FUNCTION_EVEN_DOUBLEFACTORIAL && tmp != FUNCTION_ODD_DOUBLEFACTORIAL && tmp != MAX_MEMOIZABLE_FUNCTIONS))
            {
            	mpfr_clear(tmp2);
                printUsage(&change_settings[SETTINGS_EMPTYMEMOIZERSBUFFERS]);
                return;
            }
	    }
	    else if((tmp = selectListItem(_MAX_MEMOIZABLE_FUNCTIONS, _MAX_MEMOIZABLE_FUNCTIONS > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
	            "Select a System Math Environment Memoizer to empty", suite_c.memoizers_names)) == _MAX_MEMOIZABLE_FUNCTIONS) return;
	
		if(tmp == MAX_MEMOIZABLE_FUNCTIONS)
			flushAllMemoizersBuffers();
		else
			_flushMemoizersBuffers(tmp);
	
	    return;
	}

#endif

#ifndef __DISABLE_BUFFERSEDIT
	__MSSHELL_WRAPPER_ static void   emptyBuffers(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   emptyBuffers(const sel_typ argc, char ** argv)
	{
	
	    sel_typ tmp;
	    
	    const char buffers_names[_MAX_STD_BUFFERS][INFO_STRING] =
		{
			"stdin",
			"stdout",
			"stderr",
			"Every"
		};
	
	    if(argc)
	    {
	        mpfr_t tmp2;
            if((!parse(argv[0], &tmp2)) || mpfr_cmp_ui(tmp2, (tmp = mpfr_get_ui(tmp2, MPFR_RNDN))) || (tmp != 0 && tmp != 1))
            {
            	mpfr_clear(tmp2);
                printUsage(&change_settings[SETTINGS_EMPTYBUFFERS]);
                return;
            }
            mpfr_clear(tmp2);
	    }
	    else if((tmp = selectListItem(_MAX_STD_BUFFERS, _MAX_STD_BUFFERS > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
	            "Select a System Buffer to empty", buffers_names)) == _MAX_STD_BUFFERS) return;
		
		FILE * buffers[MAX_STD_BUFFERS] =
		{
			stdin,
			stdout,
			stderr
		};
	
		if(tmp == MAX_STD_BUFFERS)
		{
			fflush(stdin);
	    	fflush(stdout);
	    	fflush(stderr);
		}
		else
	    	fflush(buffers[tmp]);
	
	    msyprintf(COLOR_SYSTEM, "\n%s BUFFER has been properly flushed.\n", buffers_names[tmp]);
	    return;
	}
#endif

#ifndef __DISABLE_SETDEFAULTS
	__MSSHELL_WRAPPER_ static void   setDefaultValues(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void   setDefaultValues(const sel_typ argc, char ** argv)
	{
	    setDefaults();
	    msprintf(COLOR_USER, "Default Settings has been correctly loaded:\n");
	    viewProgramSettings(access(lists)[LAYOUTS].cur_item);
	    flushAllMemoizersBuffers();
	    PRINTL();
	    return;
	}
#endif

sprog change_settings[MAX_SETTINGS] =
{
	#ifndef __DISABLE_SYSTEM
    [SETTINGS_SYSMANAGER] =
    {
    	CMD_SYSMANAGER,
        NAME_SYSMANAGER,
        USAGE_SYSMANAGER,
        #ifndef __DISABLE_DATABASE
        LEVEL_SYSMANAGER,
        #endif
        sysManager,
        ARGC_SYSMANAGER,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifndef __DISABLE_VARLISTMANAGER
    [SETTINGS_ENVSMANAGER] =
    {
    	CMD_ENVSMANAGER,
        NAME_ENVSMANAGER,
        USAGE_ENVSMANAGER,
        #ifndef __DISABLE_DATABASE
        LEVEL_ENVSMANAGER,
        #endif
        envsManager,
        ARGC_ENVSMANAGER,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifndef __DISABLE_LOGSMANAGER
    [SETTINGS_LOGSMANAGER] =
    {
    	CMD_LOGSMANAGER,
        NAME_LOGSMANAGER,
        USAGE_LOGSMANAGER,
        #ifndef __DISABLE_DATABASE
        LEVEL_LOGSMANAGER,
        #endif
        logsManager,
        ARGC_LOGSMANAGER,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifndef __DISABLE_SYSLOGMANAGER
    [SETTINGS_SYSLOGMANAGER] =
    {
    	CMD_SYSLOGMANAGER,
        NAME_SYSLOGMANAGER,
        USAGE_SYSLOGMANAGER,
        #ifndef __DISABLE_DATABASE
        LEVEL_SYSLOGMANAGER,
        #endif
        sysLogManager,
        ARGC_SYSLOGMANAGER,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifndef __DISABLE_LAYOUTSMANAGER
    [SETTINGS_LAYOUTSMANAGER] =
    {
    	CMD_LAYOUTSMANAGER,
        NAME_LAYOUTSMANAGER,
        USAGE_LAYOUTSMANAGER,
        #ifndef __DISABLE_DATABASE
        LEVEL_LAYOUTSMANAGER,
        #endif
        layoutsManager,
        ARGC_LAYOUTSMANAGER,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifndef __DISABLE_COLSMANAGER
    [SETTINGS_COLORSMANAGER] =
    {
    	CMD_COLORSMANAGER,
        NAME_COLORSMANAGER,
        USAGE_COLORSMANAGER,
        #ifndef __DISABLE_DATABASE
        LEVEL_COLORSMANAGER,
        #endif
        colorsManager,
        ARGC_COLORSMANAGER,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifndef __DISABLE_LFSMANAGER
    [SETTINGS_LFSMANAGER] =
    {
    	CMD_LFSMANAGER,
        NAME_LFSMANAGER,
        USAGE_LFSMANAGER,
        #ifndef __DISABLE_DATABASE
        LEVEL_LFSMANAGER,
        #endif
        lfsManager,
        ARGC_LFSMANAGER,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_PRECEDIT
    [SETTINGS_CHANGEPRECISION] =
    {
    	CMD_CHANGEPRECISION,
        NAME_CHANGEPRECISION,
        USAGE_CHANGEPRECISION,
        #ifndef __DISABLE_DATABASE
        LEVEL_CHANGEPRECISION,
        #endif
        changePrecision,
        ARGC_CHANGEPRECISION,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_STABFACTEDIT
    [SETTINGS_CHANGESTABILIZERFACTOR] =
    {
    	CMD_CHANGESTABFACT,
        NAME_CHANGESTABFACT,
        USAGE_CHANGESTABFACT,
        #ifndef __DISABLE_DATABASE
        LEVEL_CHANGESTABFACT,
        #endif
        changeStabilizerFactor,
        ARGC_CHANGESTABFACT,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_BLOCKSIZEEDIT
    {
    	CMD_CHANGEBLOCKSIZE,
    	NAME_CHANGEBLOCKSIZE,
    	USAGE_CHANGEBLOCKSIZE,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_CHANGEBLOCKSIZE,
    	#endif
    	changeBlockSize,
    	ARGC_CHANGEBLOCKSIZE,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifndef __DISABLE_MINOSMMDIMEDIT
    {
    	CMD_CHANGEMINOSMMDIM,
    	NAME_CHANGEMINOSMMDIM,
    	USAGE_CHANGEMINOSMMDIM,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_CHANGEMINOSMMDIM,
    	#endif
    	changeMinOSMMDim,
    	ARGC_CHANGEMINOSMMDIM,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifndef __DISABLE_MINSTRASSENDIMEDIT
    {
    	CMD_CHANGEMINSTRASSENDIM,
    	NAME_CHANGEMINSTRASSENDIM,
    	USAGE_CHANGEMINSTRASSENDIM,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_CHANGEMINSTRASSENDIM,
    	#endif
    	changeMinStrassenDim,
    	ARGC_CHANGEMINSTRASSENDIM,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifndef __DISABLE_MINSRNUMBEREDIT
    {
    	CMD_CHANGEMINSRNUMBER,
    	NAME_CHANGEMINSRNUMBER,
    	USAGE_CHANGEMINSRNUMBER,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_CHANGEMINSRNUMBER,
    	#endif
    	changeMinSRNumber,
    	ARGC_CHANGEMINSRNUMBER,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifndef __DISABLE_ALGEBRAEDIT
    [SETTINGS_CHANGEALGEBRA] =
    {
    	CMD_CHANGEALGEBRA,
        NAME_CHANGEALGEBRA,
        USAGE_CHANGEALGEBRA,
        #ifndef __DISABLE_DATABASE
        LEVEL_CHANGEALGEBRA,
        #endif
        changeAlgebra,
        ARGC_CHANGEALGEBRA,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_OUTLIERCONSTEDIT
    {
    	CMD_CHANGEOUTLIERCONST,
    	NAME_CHANGEOUTLIERCONST,
    	USAGE_CHANGEOUTLIERCONST,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_CHANGEOUTLIERCONST,
    	#endif
    	changeOutlierConst,
    	ARGC_CHANGEOUTLIERCONST,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifndef __DISABLE_EXITCHAREDIT
    [SETTINGS_CHANGEEXITCHAR] =
    {
    	CMD_CHANGEEXITCHAR,
        NAME_CHANGEEXITCHAR,
        USAGE_CHANGEEXITCHAR,
        #ifndef __DISABLE_DATABASE
        LEVEL_CHANGEEXITCHAR,
        #endif
        changeExitChar,
        ARGC_CHANGEEXITCHAR,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_RANDOMSEEDEDIT
    [SETTINGS_CHANGERANDOMSEED] =
    {
    	CMD_CHANGERANDOMSEED,
        NAME_CHANGERANDOMSEED,
        USAGE_CHANGERANDOMSEED,
        #ifndef __DISABLE_DATABASE
        LEVEL_CHANGERANDOMSEED,
        #endif
        changeRandomSeed,
        ARGC_CHANGERANDOMSEED,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_BOOLVARSEDIT
    [SETTINGS_CHANGEBOOLVALUES] =
    {
    	CMD_CHANGEBOOLVALUES,
        NAME_CHANGEBOOLVALUES,
        USAGE_CHANGEBOOLVALUES,
        #ifndef __DISABLE_DATABASE
        LEVEL_CHANGEBOOLVALUES,
        #endif
        changeBoolValues,
        ARGC_CHANGEBOOLVALUES,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_MMIEDIT
    [SETTINGS_CHANGEMAXMEMOIZABLEINDICES] =
    {
    	CMD_CHANGEMAXMEMIDX,
        NAME_CHANGEMAXMEMIDX,
        USAGE_CHANGEMAXMEMIDX,
        #ifndef __DISABLE_DATABASE
        LEVEL_CHANGEMAXMEMIDX,
        #endif
        changeMaxMemoizableIndices,
        ARGC_CHANGEMAXMEMIDX,
        AUTOMATIC,
        CHILD
    },
    [SETTINGS_EMPTYMEMOIZERSBUFFERS] =
    {
    	CMD_EMPTYMEMOIZERSBUFFERS,
        NAME_EMPTYMEMOIZERSBUFFERS,
        USAGE_EMPTYMEMOIZERSBUFFERS,
        #ifndef __DISABLE_DATABASE
        LEVEL_EMPTYMEMOIZERSBUFFERS,
        #endif
        emptyMemoizersBuffers,
        ARGC_EMPTYMEMOIZERSBUFFERS,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_BUFFERSEDIT
    [SETTINGS_EMPTYBUFFERS] =
    {
    	CMD_EMPTYBUFFERS,
        NAME_EMPTYBUFFERS,
        USAGE_EMPTYBUFFERS,
        #ifndef __DISABLE_DATABASE
        LEVEL_EMPTYBUFFERS,
        #endif
        emptyBuffers,
        ARGC_EMPTYBUFFERS,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifndef __DISABLE_SETDEFAULTS
    [SETTINGS_SETDEFAULTS] =
    {
    	CMD_SETDEFAULTS,
        NAME_SETDEFAULTS,
        USAGE_SETDEFAULTS,
        #ifndef __DISABLE_DATABASE
        LEVEL_SETDEFAULTS,
        #endif
        setDefaultValues,
        ARGC_SETDEFAULTS,
        AUTOMATIC,
        CHILD
    }
    #endif
};

#endif
