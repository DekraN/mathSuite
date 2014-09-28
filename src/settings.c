// settings.c 16/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#ifndef FREEZES_SETTINGSMANAGER

#include "dutils.h"

// Settings Function Declarations and Definitions
//

#ifdef ALLOW_VARLISTMANAGER
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt envsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt envsManager(const sel_typ argc, char ** argv)
	{
	
	    operationsGroupMenu(MAX_ENVSMANAGER_PROGS,
	                        envs_manager, change_settings[SETTINGS_ENVSMANAGER].name,
	                        #if MAX_ENVSMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
	                            BY_NUMBERS
	                        #else
	                            BY_CHARS
	                        #endif
	                        );
	    return;
	}
#endif

#ifdef ALLOW_LOGSMANAGER
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt logsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt logsManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_LOGSMANAGER_PROGS,
	                        logs_manager, change_settings[SETTINGS_LOGSMANAGER].name,
	                        #if MAX_LOGSMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
	                            BY_NUMBERS
	                        #else
	                            BY_CHARS
	                        #endif
	                        );
	    return;
	}
#endif

#ifdef ALLOW_SYSLOGMANAGER
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt sysLogManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt sysLogManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_SYSLOGMANAGER_PROGS,
	                        syslog_manager, change_settings[SETTINGS_SYSLOGMANAGER].name,
	                        #if MAX_SYSLOGMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
	                            BY_NUMBERS
	                        #else
	                            BY_CHARS
	                        #endif
	                        );
	    return;
	}
#endif

#ifdef ALLOW_LAYOUTSMANAGER
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt layoutsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt layoutsManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_LAYOUTSMANAGER_PROGS,
	                        layouts_manager, change_settings[SETTINGS_LAYOUTSMANAGER].name,
	                        #if MAX_LAYOUTSMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
	                            BY_NUMBERS
	                        #else
	                            BY_CHARS
	                        #endif
	                        );
	    return;
	}
#endif

#ifdef ALLOW_COLSMANAGER
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt colorsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt colorsManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_COLSMANAGER_PROGS,
	                        cols_manager, change_settings[SETTINGS_COLORSMANAGER].name,
	                        #if MAX_COLSMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
	                            BY_NUMBERS
	                        #else
	                            BY_CHARS
	                        #endif
	                        );
	    return;
	}
#endif
	
#ifdef ALLOW_LFSMANAGER
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt lfsManager(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system __apnt lfsManager(const sel_typ argc, char ** argv)
	{
	    operationsGroupMenu(MAX_LFSMANAGER_PROGS,
	                    lfs_manager, change_settings[SETTINGS_LFSMANAGER].name,
	                    #if MAX_COLSMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
	                        BY_NUMBERS
	                    #else
	                        BY_CHARS
	                    #endif
	                    );
	    return;
	}
#endif

#ifdef ALLOW_PRECEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changePrecision(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changePrecision(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_precision = access(curLayout)->precision;
	    ityp tmp;
	
	    // uint64_t tmp2;
	
	    if(argc)
	    {
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	
	            if((!parse(argv[0], &tmp)) || tmp != (access(curLayout)->precision = (fsel_typ)tmp) || access(curLayout)->precision < MIN_PRECISION || access(curLayout)->precision > MAX_PRECISION)
	            {
	            	access(curLayout)->precision = old_precision;
	                printErr(33, "Invalid inserted Precision Value.\nMust be an integer between %hu and %hu", MIN_PRECISION, MAX_PRECISION);
	                printUsage(&change_settings[SETTINGS_CHANGEPRECISION]);
	                return;
	            }
	        }
	        else if((tmp = strtod(argv[0], NULL)) != (access(curLayout)->precision = (fsel_typ)tmp) || access(curLayout)->precision < MIN_PRECISION || access(curLayout)->precision > MAX_PRECISION)
	        {
	        	access(curLayout)->precision = old_precision;
	            printErr(33, "Invalid inserted Precision Value.\nMust be an integer between %hu and %hu", MIN_PRECISION, MAX_PRECISION);
	            printUsage(&change_settings[SETTINGS_CHANGEPRECISION]);
	            return;
	        }
	    }
	    else
	    {
	        printf2(COLOR_CREDITS, "\nEnter a Value as Program Precision: between %hu and %hu.\n", MIN_PRECISION, MAX_PRECISION);
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	        while((PARSING_SYSTEM_ALLOWED ? ((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)) == NULL_VAL) :
	            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (access(curLayout)->precision = (fsel_typ)tmp) || access(curLayout)->precision < MIN_PRECISION || access(curLayout)->precision > MAX_PRECISION)
	        {
	            CLEARBUFFER();
	
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck)
	            {
	                access(curLayout)->precision = old_precision;
	                return;
	            }
	            printErr(33, "Invalid inserted Precision Value.\nMust be an integer between %hu and %hu", MIN_PRECISION, MAX_PRECISION);
	        }
	    }
	
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        *(access(exprVars)->e_ANS) = access(curLayout)->precision;
	
	    sprint("Program Precision Value has been correctly changed from: %hu to: %hu.\n", old_precision, access(curLayout)->precision);
	    return;
	}
#endif

#ifdef ALLOW_STABFACTEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeStabilizerFactor(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeStabilizerFactor(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_stabilizer_factor = access(curLayout)->stabilizer_factor;
	    ityp tmp = 0.00;
	
	    if(argc)
	    {
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || tmp != (access(curLayout)->stabilizer_factor = (fsel_typ)tmp) || access(curLayout)->stabilizer_factor < MIN_STABFACT || access(curLayout)->stabilizer_factor > MAX_STABFACT)
	            {
	            	access(curLayout)->stabilizer_factor = old_stabilizer_factor;
	                printErr(33, "Invalid inserted Stabilizer Factor Value.\nMust be an integer between %hu and %hu", MIN_STABFACT, MAX_STABFACT);
	                printUsage(&change_settings[SETTINGS_CHANGESTABILIZERFACTOR]);
	                return;
	            }
	        }
	        else if((tmp = strtod(argv[0], NULL)) != (access(curLayout)->stabilizer_factor = (fsel_typ)tmp) || access(curLayout)->stabilizer_factor < MIN_STABFACT || access(curLayout)->stabilizer_factor > MAX_STABFACT)
	        {
	        	access(curLayout)->stabilizer_factor = old_stabilizer_factor;
	            printErr(33, "Invalid inserted Stabilizer Factor Value.\nMust be an integer between %hu and %hu", MIN_STABFACT, MAX_STABFACT);
	            printUsage(&change_settings[SETTINGS_CHANGESTABILIZERFACTOR]);
	            return;
	        }
	    }
	    else
	    {
	        printf2(COLOR_CREDITS, "\nEnter a Value as Program Stabilizer Factor: a non-negative and non-zero integer.\n");
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	        while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
	            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (access(curLayout)->stabilizer_factor = (fsel_typ)tmp) || access(curLayout)->stabilizer_factor < MIN_STABFACT || access(curLayout)->stabilizer_factor > MAX_STABFACT)
	        {
	            CLEARBUFFER();
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck)
	            {
	                access(curLayout)->stabilizer_factor = old_stabilizer_factor;
	                return;
	            }
	            printErr(33, "Invalid inserted Stabilizer Factor Value.\nMust be an integer between %hu and %hu", MIN_STABFACT, MAX_STABFACT);
	        }
	    }
	
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        *(access(exprVars)->e_ANS) = access(curLayout)->stabilizer_factor;
	
	    sprint("Program Stabilizer Factor Value has been correctly changed from: %hu to: %hu.\n", old_stabilizer_factor, access(curLayout)->stabilizer_factor);
	    return;
	}
#endif

#ifdef ALLOW_BLOCKSIZEEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeBlockSize(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeBlockSize(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_block_size= access(curLayout)->block_size;
	    ityp tmp = 0.00;
	
	    if(argc)
	    {
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || tmp != (access(curLayout)->block_size = (fsel_typ)tmp) || access(curLayout)->block_size < MIN_BLOCKSIZE)
	            {
	            	access(curLayout)->block_size = old_block_size;
	                printErr(33, "Invalid inserted Block Size Value.\nMust be an integer greater or equal than %hu", MIN_BLOCKSIZE);
	                printUsage(&change_settings[SETTINGS_CHANGEBLOCKSIZE]);
	                return;
	            }
	        }
	        else if((tmp = strtod(argv[0], NULL)) != (access(curLayout)->block_size = (fsel_typ)tmp) || access(curLayout)->block_size < MIN_BLOCKSIZE)
	        {
	        	access(curLayout)->block_size = old_block_size;
	            printErr(33, "Invalid inserted Block Size Value.\nMust be an integer greater or equal than %hu", MIN_BLOCKSIZE);
	            printUsage(&change_settings[SETTINGS_CHANGEBLOCKSIZE]);
	            return;
	        }
	    }
	    else
	    {
	        printf2(COLOR_CREDITS, "\nEnter a Value as Program Block Size: a non-negative and non-zero integer.\n");
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	        while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
	            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (access(curLayout)->block_size = (fsel_typ)tmp) || access(curLayout)->block_size < MIN_BLOCKSIZE)
	        {
	            CLEARBUFFER();
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck)
	            {
	                access(curLayout)->block_size = old_block_size;
	                return;
	            }
	            printErr(33, "Invalid inserted Block Size Value.\nMust be an integer greater or equal than %hu", MIN_BLOCKSIZE);
	        }
	    }
	
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        *(access(exprVars)->e_ANS) = access(curLayout)->block_size;
	
	    sprint("Program Block Size Value has been correctly changed from: %hu to: %hu.\n", old_block_size, access(curLayout)->block_size);
	    return;
	}
#endif

#ifdef ALLOW_MINOSMMDIMEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeMinOSMMDim(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeMinOSMMDim(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_min_osmm_dim= access(curLayout)->min_osmm_dim;
	    ityp tmp = 0.00;
	
	    if(argc)
	    {
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || tmp != (access(curLayout)->min_osmm_dim = (fsel_typ)tmp) || access(curLayout)->min_osmm_dim < MIN_OSMM_DIM)
	            {
	            	access(curLayout)->min_osmm_dim = old_min_osmm_dim;
	                printErr(33, "Invalid inserted Min OSMM Dimension Value.\nMust be an integer greater or equal than %hu", MIN_OSMM_DIM);
	                printUsage(&change_settings[SETTINGS_CHANGEMINOSMMDIM]);
	                return;
	            }
	        }
	        else if((tmp = strtod(argv[0], NULL)) != (access(curLayout)->min_osmm_dim = (fsel_typ)tmp) || access(curLayout)->min_osmm_dim < MIN_OSMM_DIM)
	        {
	        	access(curLayout)->min_osmm_dim = old_min_osmm_dim;
	            printErr(33, "Invalid inserted Min OSMM Dimension Value.\nMust be an integer greater or equal than %hu", MIN_OSMM_DIM);
	            printUsage(&change_settings[SETTINGS_CHANGEMINOSMMDIM]);
	            return;
	        }
	    }
	    else
	    {
	        printf2(COLOR_CREDITS, "\nEnter a Value as Min OSMM Dimension: a non-negative and non-zero integer.\n");
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	        while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
	            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (access(curLayout)->min_osmm_dim = (fsel_typ)tmp) || access(curLayout)->min_osmm_dim < MIN_OSMM_DIM)
	        {
	            CLEARBUFFER();
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck)
	            {
	                access(curLayout)->min_osmm_dim = old_min_osmm_dim;
	                return;
	            }
	            printErr(33, "Invalid inserted Min OSMM Dimension Value.\nMust be an integer greater or equal than %hu", MIN_OSMM_DIM);
	        }
	    }
	
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        *(access(exprVars)->e_ANS) = access(curLayout)->min_osmm_dim;
	
	    sprint("Program Min OSMM Dimension Value has been correctly changed from: %hu to: %hu.\n", old_min_osmm_dim, access(curLayout)->min_osmm_dim);
	    return;
	}
#endif

#ifdef ALLOW_MINSTRASSENDIMEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeMinStrassenDim(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeMinStrassenDim(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_min_strassen_dim= access(curLayout)->min_strassen_dim;
	    ityp tmp = 0.00;
	
	    if(argc)
	    {
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || tmp != (access(curLayout)->min_strassen_dim = (fsel_typ)tmp) || access(curLayout)->min_strassen_dim < MIN_STRASSEN_DIM)
	            {
	            	access(curLayout)->min_strassen_dim = old_min_strassen_dim;
	                printErr(33, "Invalid inserted Min Strassen Dimension Value.\nMust be an integer greater or equal than %hu", MIN_STRASSEN_DIM);
	                printUsage(&change_settings[SETTINGS_CHANGEMINSTRASSENDIM]);
	                return;
	            }
	        }
	        else if((tmp = strtod(argv[0], NULL)) != (access(curLayout)->min_strassen_dim = (fsel_typ)tmp) || access(curLayout)->min_strassen_dim < MIN_STRASSEN_DIM)
	        {
	        	access(curLayout)->min_strassen_dim = old_min_strassen_dim;
	            printErr(33, "Invalid inserted Min Strassen Dimension Value.\nMust be an integer greater or equal than %hu", MIN_STRASSEN_DIM);
	            printUsage(&change_settings[SETTINGS_CHANGEMINSTRASSENDIM]);
	            return;
	        }
	    }
	    else
	    {
	        printf2(COLOR_CREDITS, "\nEnter a Value as Min Strassen Dimension: a non-negative and non-zero integer.\n");
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	        while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
	            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (access(curLayout)->min_strassen_dim = (fsel_typ)tmp) || access(curLayout)->min_strassen_dim < MIN_STRASSEN_DIM)
	        {
	            CLEARBUFFER();
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck)
	            {
	                access(curLayout)->min_strassen_dim = old_min_strassen_dim;
	                return;
	            }
	            printErr(33, "Invalid inserted Min Strassen Dimension Value.\nMust be an integer greater or equal than %hu", MIN_STRASSEN_DIM);
	        }
	    }
	
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        *(access(exprVars)->e_ANS) = access(curLayout)->min_strassen_dim;
	
	    sprint("Program Min Strassen Dimension Value has been correctly changed from: %hu to: %hu.\n", old_min_strassen_dim, access(curLayout)->min_strassen_dim);
	    return;
	}
#endif

#ifdef ALLOW_MINSRNUMBEREDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeMinSRNumber(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeMinSRNumber(const sel_typ argc, char ** argv)
	{
		const dim_typ old_minsr_number = access(curLayout)->min_stirling_number;
	    ityp tmp;
	
	    if(argc)
	    {
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || tmp != (access(curLayout)->min_stirling_number = (dim_typ)tmp) || access(curLayout)->min_stirling_number < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->min_stirling_number > MIN_STIRLING_NUMBER)
	            {
	            	access(curLayout)->min_stirling_number = old_minsr_number;
	                printErr(33, "Invalid inserted Min Stirling Number.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MIN_STIRLING_NUMBER);
	                printUsage(&change_settings[SETTINGS_CHANGEMINSRNUMBER]);
	                return;
	            }
	        }
	        else if((tmp = strtod(argv[0], NULL)) != (access(curLayout)->min_stirling_number = (dim_typ)tmp) || access(curLayout)->min_stirling_number < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->min_stirling_number > MIN_STIRLING_NUMBER)
	        {
	        	access(curLayout)->min_stirling_number = old_minsr_number;
	            printErr(33, "Invalid inserted Min Stirling Number.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MIN_STIRLING_NUMBER);
	            printUsage(&change_settings[SETTINGS_CHANGEMINSRNUMBER]);
	            return;
	        }
	    }
	    else
	    {
	        printf2(COLOR_CREDITS, "\nEnter a Value as Min Stirling Number: a non-negative and non-zero integer.\n");
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	        while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
	            (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (access(curLayout)->min_stirling_number = (dim_typ)tmp) || access(curLayout)->min_stirling_number < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->min_stirling_number > MIN_STIRLING_NUMBER)
	        {
	            CLEARBUFFER();
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck)
	            {
	                access(curLayout)->min_stirling_number = old_minsr_number;
	                return;
	            }
	            printErr(33, "Invalid inserted Min Stirling Number.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MIN_STIRLING_NUMBER);
	        }
	    }
	
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        *(access(exprVars)->e_ANS) = access(curLayout)->min_stirling_number;
	
	    sprint("Program Min Stirling Number has been correctly changed from: %hu to: %hu.\n", old_minsr_number, access(curLayout)->min_stirling_number);
		return;
	}

#endif

#ifdef ALLOW_ALGEBRAEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeAlgebra(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeAlgebra(const sel_typ argc, char ** argv)
	{
	
	    dim_typ i;
	    
	    if(argc)
	    {
	        ityp tmp = 0.00;
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || tmp != (i = (dim_typ)tmp) || i < MIN_ALGEBRA || i > MAX_ALGEBRA)
	            {
	                printErr(1, "Invalid inserted Value: not correspondent to any Algebra Identifier");
	                printUsage(&change_settings[SETTINGS_CHANGEALGEBRA]);
	                return;
	            }
	        }
	        else if((tmp = strtod(argv[0], NULL)) != (i = (dim_typ)tmp) || i < MIN_ALGEBRA || i > MAX_ALGEBRA)
	        {
	            printErr(1, "Invalid inserted Value: not correspondent to any Algebra Identifier");
	            printUsage(&change_settings[SETTINGS_CHANGEALGEBRA]);
	            return;
	        }
	    }
	    else if((i = selectListItem(_MAX_ALGEBRA, _MAX_ALGEBRA > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
	            "Select Algebra Identifier in which to perform Algebra Operations", suite_c.algebra_elements_names)) == _MAX_ALGEBRA) return;
	
	    access(curLayout)->algebra = i;
	    sprint("%s Algebra has been correctly selected.\n\n", suite_c.algebra_elements_names[i]);
	    return;
	}
#endif

#ifdef ALLOW_OUTLIERCONSTEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeOutlierConst(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeOutlierConst(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_outlier_constant = access(curLayout)->outlier_constant;
		ityp tmp = 0.00; 
	    if(argc)
	    {
	    	ityp tmp = 0.00;
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || (access(curLayout)->outlier_constant=tmp) < MIN_OUTLIER_CONSTANT || access(curLayout)->outlier_constant > MAX_OUTLIER_CONSTANT)
	            {
	            	access(curLayout)->outlier_constant = old_outlier_constant;
	                printErr(33, "Invalid inserted Outlier Constant.\nMust be a float between %.*f and %.*f", DEFAULT_PRECISION, MIN_OUTLIER_CONSTANT, DEFAULT_PRECISION, MAX_OUTLIER_CONSTANT);
	                printUsage(&change_settings[SETTINGS_CHANGEOUTLIERCONST]);
	                return;
	            }
	        }
	        else if((access(curLayout)->outlier_constant = strtof(argv[0], NULL)) < MIN_STABFACT || access(curLayout)->outlier_constant > MAX_STABFACT)
	        {
	        	access(curLayout)->outlier_constant = old_outlier_constant;
	            printErr(33, "Invalid inserted Outlier Constant.\nMust be a float between %.*f and %.*f", DEFAULT_PRECISION, MIN_OUTLIER_CONSTANT, MAX_OUTLIER_CONSTANT);
	            printUsage(&change_settings[SETTINGS_CHANGEOUTLIERCONST]);
	            return;
	        }
	    }
	    else
	    {
	        printf2(COLOR_CREDITS, "\nEnter a Value as Program Outlier Constant: a non-negative and non-zero float.\n");
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	        while((PARSING_SYSTEM_ALLOWED ? (isNullVal((access(curLayout)->outlier_constant = (float)requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
	            (!scanf2(1, INPUT_CONVERSION_FORMAT, &access(curLayout)->outlier_constant))) || access(curLayout)->outlier_constant < MIN_OUTLIER_CONSTANT || access(curLayout)->outlier_constant > MAX_OUTLIER_CONSTANT)
	        {
	            CLEARBUFFER();
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck)
	            {
	                access(curLayout)->outlier_constant = old_outlier_constant;
	                return;
	            }
	            printErr(33, "Invalid inserted Outlier Constant.\nMust be a float between %.*f and %.*f", DEFAULT_PRECISION, MIN_OUTLIER_CONSTANT, DEFAULT_PRECISION, MAX_OUTLIER_CONSTANT);
	        }
	    }
	
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        *(access(exprVars)->e_ANS) = access(curLayout)->outlier_constant;
	
	    sprint("Program Outlier Constant has been correctly changed from: %.*f to: %.*f.\n", DEFAULT_PRECISION, old_outlier_constant, DEFAULT_PRECISION, access(curLayout)->outlier_constant);
	    return;
	}
#endif

#ifdef ALLOW_EXITCHAREDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeExitChar(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeExitChar(const sel_typ argc, char ** argv)
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
	    sprint("Back to MAIN MENU Character has been correctly changed from: %c to: %c.\n", old_exit_char, access(curLayout)->exit_char);
	    return;
	}
#endif

#ifdef ALLOW_RANDOMSEEDEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeRandomSeed(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeRandomSeed(const sel_typ argc, char ** argv)
	{
	    const fsel_typ old_random_seed = access(random_seed);
	    ityp tmp;
	
	    if(argc)
	    {
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || tmp != (access(random_seed) = (fsel_typ)tmp) || access(random_seed) < MIN_RANDOMSEED || access(random_seed) > MAX_RANDOMSEED)
	            {
	            	access(random_seed) = old_random_seed;
	                printErr(5, "Invalid inserted Value.\nMust be an integer between %hu and %hu", MIN_RANDOMSEED, MAX_RANDOMSEED); // strettamente positivo e minore di %u", UINT_MAX);
	                printUsage(&change_settings[SETTINGS_CHANGERANDOMSEED]);
	                return;
	            }
	        }
	        else if((tmp = strtoul(argv[0], NULL, sizeof(access(random_seed)))) != (access(random_seed) = (fsel_typ)tmp) || access(random_seed) < MIN_RANDOMSEED || access(random_seed) > MAX_RANDOMSEED)
	        {
	        	access(random_seed) = old_random_seed;
	            printErr(5, "Invalid inserted Value.\nMust be an integer between %hu and %hu", MIN_RANDOMSEED, MAX_RANDOMSEED); // , UINT_MAX);
	            printUsage(&change_settings[SETTINGS_CHANGERANDOMSEED]);
	            return;
	        }
	    }
	    else
	    {
	        printf2(COLOR_CREDITS, "\nEnter RandomSeed Value\n");
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	
	        while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
	                (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (access(random_seed) = (fsel_typ)tmp) || access(random_seed) < MIN_RANDOMSEED || access(random_seed) > MAX_RANDOMSEED)
	        {
	            CLEARBUFFER();
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck)
	            {
	                access(random_seed) = old_random_seed;
	                return;
	            }
	            printErr(5, "Invalid inserted Value.\nMust be an integer between %hu and %hu", MIN_RANDOMSEED, MAX_RANDOMSEED); // , UINT_MAX);
	        }
	    }
	
	    CLEARBUFFER();
	
	    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS))
	        *(access(exprVars)->e_ANS) = access(random_seed);
	
	    srand(access(random_seed));
	    sprint("Program RandomSeed Value has been correctly changed from: %hu to: %hu.\n", old_random_seed, access(random_seed));
	
	    return;
	}
#endif

#ifdef ALLOW_BOOLVARSEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeBoolValues(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeBoolValues(const sel_typ argc, char ** argv)
	{
	    dim_typ i;
	    if(argc)
	    {
	        ityp tmp;
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	            if((!parse(argv[0], &tmp)) || tmp != (i = (dim_typ)tmp) || i < 0 || i > MAX_BOOL_SETTINGS)
	            {
	                printErr(1, "Invalid inserted Value: not correspondent to any Boolean Settings ID");
	                printUsage(&change_settings[SETTINGS_CHANGEBOOLVALUES]);
	                return;
	            }
	        }
	        else if((tmp = strtod(argv[0], NULL)) != (i = (dim_typ)tmp) || i < 0 || i > MAX_BOOL_SETTINGS)
	        {
	            printErr(1, "Invalid inserted Value: not correspondent to any Boolean Settings ID");
	            printUsage(&change_settings[SETTINGS_CHANGEBOOLVALUES]);
	            return;
	        }
	    }
	    else if((i = selectListItem(MAX_BOOL_SETTINGS, MAX_BOOL_SETTINGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET,
	            "Select Bool Settings you wish to change", suite_c.bools_names)) == MAX_BOOL_SETTINGS) return;
	
	    access(curLayout)->bools ^= suite_c.bools[i].bmask;
	    sprint("%s Boolean Settings\nhas been correctly %s.\n\n", suite_c.bools_names[i], isSett(i) ? "ENABLED":"DISABLED");
	
	    return;
	}
#endif

#ifdef ALLOW_MMIEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system changeMaxMemoizableIndices(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system changeMaxMemoizableIndices(const sel_typ argc, char ** argv)
	{
	
	    const dim_typ old_max_memoizable_index[MAX_MEMOIZABLE_FUNCTIONS] =
	    {
	        access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI],
	        access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE],
	        access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE],
	        access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE]
	    };
	
		dim_typ i;
		ityp tmp = 0.00;
	
	    if(argc)
	    {
	        if(argc > LAST_MEMOIZABLE_FUNCTION)
	        {
	            if(PARSING_SYSTEM_ALLOWED)
	            {
	                if((!parse(argv[FUNCTION_FIBONACCI], &tmp)) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = (dim_typ)tmp) ||
	                access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] > MAX_FIBONACCI_MEMOIZABLE_INDEX)
	                {
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
	                    printErr(33, "Invalid inserted MIM_FIBO Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_FIBONACCI_MEMOIZABLE_INDEX);
	                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	                    return;
	                }
	
	                if((!parse(argv[FUNCTION_FATTORIALE], &tmp)) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = (dim_typ)tmp) ||
	                access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] > MAX_FATTORIALE_MEMOIZABLE_INDEX)
	                {
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = old_max_memoizable_index[FUNCTION_FATTORIALE];
	                    printErr(33, "Invalid inserted MIM_FACT Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_FATTORIALE_MEMOIZABLE_INDEX);
	                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	                    return;
	                }
	                
	                if((!parse(argv[FUNCTION_EVEN_SFATTORIALE], &tmp)) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = (dim_typ)tmp) ||
	                access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] > MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX)
	                {
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = old_max_memoizable_index[FUNCTION_FATTORIALE];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = old_max_memoizable_index[FUNCTION_EVEN_SFATTORIALE];
	                    printErr(33, "Invalid inserted MIM_SFACT_EVEN Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX);
	                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	                    return;
	                }
	                
	                if((!parse(argv[FUNCTION_ODD_SFATTORIALE], &tmp)) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] = (dim_typ)tmp) ||
	                access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] > MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX)
	                {
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = old_max_memoizable_index[FUNCTION_FATTORIALE];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = old_max_memoizable_index[FUNCTION_EVEN_SFATTORIALE];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] = old_max_memoizable_index[FUNCTION_ODD_SFATTORIALE];
	                    printErr(33, "Invalid inserted MIM_SFACT_ODD Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX);
	                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	                    return;
	                }
	                
	                
	            }
	            else
	            {
	                if((tmp = strtod(argv[FUNCTION_FIBONACCI], NULL)) != (access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = (dim_typ)tmp) ||
	                access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] > MAX_FIBONACCI_MEMOIZABLE_INDEX)
	                {
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
	                    printErr(33, "Invalid inserted MIM_FIBO Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_FIBONACCI_MEMOIZABLE_INDEX);
	                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	                    return;
	                }
	
	                if((tmp = strtod(argv[FUNCTION_FATTORIALE], NULL)) != (access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = (dim_typ)tmp) ||
	                access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] > MAX_FATTORIALE_MEMOIZABLE_INDEX)
	                {
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = old_max_memoizable_index[FUNCTION_FATTORIALE];
	                    printErr(33, "Invalid inserted MIM_FACT Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_FATTORIALE_MEMOIZABLE_INDEX);
	                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	                    return;
	                }
	                
	                if((tmp = strtod(argv[FUNCTION_EVEN_SFATTORIALE], NULL)) != (access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = (dim_typ)tmp) ||
	                access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] > MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX)
	                {
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = old_max_memoizable_index[FUNCTION_FATTORIALE];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = old_max_memoizable_index[FUNCTION_EVEN_SFATTORIALE];
	                    printErr(33, "Invalid inserted MIM_SFACT_EVEN Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX);
	                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	                    return;
	                }
	                
	                if((tmp = strtod(argv[FUNCTION_ODD_SFATTORIALE], NULL)) != (access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] = (dim_typ)tmp) ||
	                access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] > MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX)
	                {
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = old_max_memoizable_index[FUNCTION_FIBONACCI];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = old_max_memoizable_index[FUNCTION_FATTORIALE];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = old_max_memoizable_index[FUNCTION_EVEN_SFATTORIALE];
	                	access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] = old_max_memoizable_index[FUNCTION_ODD_SFATTORIALE];
	                    printErr(33, "Invalid inserted MIM_SFACT_EVEN Value.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX);
	                    printUsage(&change_settings[SETTINGS_CHANGEMAXMEMOIZABLEINDICES]);
	                    return;
	                }
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
			strcpy(seperator, PARSING_SYSTEM_ALLOWED ? "]\n[" : BLANK_STRING);
	        printf2(COLOR_CREDITS, "\nEnter the three MAX MEMOIZABLE INDICES of the functions respectively\n:%s, %s, %s and %s, as expected format:\n", suite_c.memoizers_names[FUNCTION_FIBONACCI], suite_c.memoizers_names[FUNCTION_FATTORIALE], suite_c.memoizers_names[FUNCTION_EVEN_SFATTORIALE], suite_c.memoizers_names[FUNCTION_ODD_SFATTORIALE]);
	        printf2(COLOR_CREDITS, "[MIM_FIBO%sMIM_FACT%sMIM_EVEN_SFACT%sMIM_ODD_SFACT]\n", seperator, seperator, seperator);
	
	        if(PARSING_SYSTEM_ALLOWED)
	            PRINTHOWTOBACKMESSAGE();
	            
	        ityp tmp2, tmp3, tmp4;
	        
	        tmp2 = tmp3 = tmp4 = 0.00;
	
	        while(PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted MAX Fibonacci MEMOIZABLE INDEX is:", PARSER_SHOWRESULT))) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = (dim_typ)tmp) ||
	            access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] > MAX_MEMOIZABLE_INDEX) ||
	            (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted MAX Factorial MEMOIZABLE INDEX is:", PARSER_SHOWRESULT))) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = (dim_typ)tmp) ||
	            access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] > MAX_MEMOIZABLE_INDEX) ||
			    (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted MAX Even SemiFactorial MEMOIZABLE INDEX is:", PARSER_SHOWRESULT))) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = (dim_typ)tmp) ||
	            access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] > MAX_MEMOIZABLE_INDEX) ||
				(isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted MAX Odd SemiFactorial MEMOIZABLE INDEX is:", PARSER_SHOWRESULT))) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] = (dim_typ)tmp) ||
	            access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] > MAX_MEMOIZABLE_INDEX) :
	            (!scanf2(2, "%lf %lf %lf %lf", &tmp, &tmp2, &tmp3, &tmp4)) || tmp != (access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] = (dim_typ)tmp) || tmp2 != (access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] = (dim_typ)tmp2) || tmp3 != (access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] = (dim_typ)tmp3) || tmp4 != (access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] = (dim_typ)tmp4) ||
	              access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] < MIN_MEMOIZABLE_INDEX+1 || access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI] > MAX_MEMOIZABLE_INDEX ||
	              access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE] > MAX_MEMOIZABLE_INDEX || access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE] > MAX_MEMOIZABLE_INDEX || access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE] > MAX_MEMOIZABLE_INDEX)
	        {
	            CLEARBUFFER();
	            if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
	            if(exitHandleCheck) // if(tmp3[ROWS] == NULL_VAL || tmp3[COLUMNS] == NULL_VAL)
	            {
	            	#pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
	            	for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
	            		access(curLayout)->max_memoizable_indices[i] = old_max_memoizable_index[i];
	                return;
	            }
	            printErr(33, "Invalid [MIM_FIBO MIM_FACT MIM_EVENSFACT MIM_ODDSFACT] format.\nMust be an integer between %hu and %hu", MIN_MEMOIZABLE_INDEX+1, MAX_MEMOIZABLE_INDEX);
	        }
	    }
	    
	    for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
			sprint("MAX %s MEMOIZABLE INDEX Value has been changed from: %hu to: %hu.\n", old_max_memoizable_index[i], access(curLayout)->max_memoizable_indices[i]);
	
	    return;
	}
	
	__MSSHELL_WRAPPER_ static void _MS__private __system emptyMemoizersBuffers(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system emptyMemoizersBuffers(const sel_typ argc, char ** argv)
	{
	
	    sel_typ tmp;
	
	    if(argc)
	    {
	        ityp tmp2 = 0.00;
	
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	
	            if((!parse(argv[0], &tmp2)) || tmp2 != (tmp = (sel_typ)tmp2) || (tmp != FUNCTION_FIBONACCI && tmp != FUNCTION_FATTORIALE && tmp != FUNCTION_EVEN_SFATTORIALE && tmp != FUNCTION_ODD_SFATTORIALE && tmp != MAX_MEMOIZABLE_FUNCTIONS))
	            {
	                printUsage(&change_settings[SETTINGS_EMPTYMEMOIZERSBUFFERS]);
	                return;
	            }
	        }
	        else if((tmp2 = strtod(argv[0], NULL)) != (tmp = (sel_typ)tmp2) || (tmp != FUNCTION_FIBONACCI && tmp != FUNCTION_FATTORIALE && tmp != FUNCTION_EVEN_SFATTORIALE && tmp != FUNCTION_ODD_SFATTORIALE && tmp != MAX_MEMOIZABLE_FUNCTIONS))
	        {
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

#ifdef ALLOW_BUFFERSEDIT
	__MSSHELL_WRAPPER_ static void _MS__private __system emptyBuffers(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system emptyBuffers(const sel_typ argc, char ** argv)
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
	        ityp tmp2 = 0.00;
	
	        if(PARSING_SYSTEM_ALLOWED)
	        {
	
	            if((!parse(argv[0], &tmp2)) || tmp2 != (tmp = (sel_typ)tmp2) || (tmp != 0 && tmp != 1))
	            {
	                printUsage(&change_settings[SETTINGS_EMPTYBUFFERS]);
	                return;
	            }
	        }
	        else if((tmp2 == strtod(argv[0], NULL)) != (tmp = (sel_typ)tmp2) || (tmp != 0 && tmp != 1))
	        {
	            printUsage(&change_settings[SETTINGS_EMPTYBUFFERS]);
	            return;
	        }
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
	
	    sprint("\n%s BUFFER has been properly flushed.\n", buffers_names[tmp]);
	    return;
	}
#endif

#ifdef ALLOW_SETDEFAULTS
	__MSSHELL_WRAPPER_ static void _MS__private __system setDefaultValues(const sel_typ argc, char ** argv);
	__MSSHELL_WRAPPER_ static void _MS__private __system setDefaultValues(const sel_typ argc, char ** argv)
	{
	    setDefaults();
	    printf2(COLOR_USER, "Default Settings has been correctly loaded:\n");
	    viewProgramSettings(access(lists)[LAYOUTS].cur_item);
	    flushAllMemoizersBuffers();
	    PRINTL();
	    return;
	}
#endif

sprog change_settings[MAX_SETTINGS] =
{
    #ifdef ALLOW_VARLISTMANAGER
    [SETTINGS_ENVSMANAGER] =
    {
        "Variables Environments Manager",
        CMD_ENVSMANAGER,
        USAGE_ENVSMANAGER,
        envsManager,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifdef ALLOW_LOGSMANAGER
    [SETTINGS_LOGSMANAGER] =
    {
        "Logs Manager",
        CMD_LOGSMANAGER,
        USAGE_LOGSMANAGER,
        logsManager,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifdef ALLOW_SYSLOGMANAGER
    [SETTINGS_SYSLOGMANAGER] =
    {
        "System Log Manager",
        CMD_SYSLOGMANAGER,
        USAGE_SYSLOGMANAGER,
        sysLogManager,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifdef ALLOW_LAYOUTSMANAGER
    [SETTINGS_LAYOUTSMANAGER] =
    {
        "Layouts Manager",
        CMD_LAYOUTSMANAGER,
        USAGE_LAYOUTSMANAGER,
        layoutsManager,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifdef ALLOW_COLSMANAGER
    [SETTINGS_COLORSMANAGER] =
    {
        "Program Colors Manager",
        CMD_COLORSMANAGER,
        USAGE_COLORSMANAGER,
        colorsManager,
        AUTOMATIC,
        FATHER
    },
    #endif
    #ifdef ALLOW_LFSMANAGER
    [SETTINGS_LFSMANAGER] =
    {
        "Startup Files Manager",
        CMD_LFSMANAGER,
        USAGE_LFSMANAGER,
        lfsManager,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_PRECEDIT
    [SETTINGS_CHANGEPRECISION] =
    {
        "Change Program Precision",
        CMD_CHANGEPRECISION,
        USAGE_CHANGEPRECISION,
        changePrecision,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_STABFACTEDIT
    [SETTINGS_CHANGESTABILIZERFACTOR] =
    {
        "Change Stabilizer Factor",
        CMD_CHANGESTABFACT,
        USAGE_CHANGESTABFACT,
        changeStabilizerFactor,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_BLOCKSIZEEDIT
    {
    	"Change Block Size",
    	CMD_CHANGEBLOCKSIZE,
    	USAGE_CHANGEBLOCKSIZE,
    	changeBlockSize,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifdef ALLOW_MINOSMMDIMEDIT
    {
    	"Change Min OSMM Dimension",
    	CMD_CHANGEMINOSMMDIM,
    	USAGE_CHANGEMINOSMMDIM,
    	changeMinOSMMDim,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifdef ALLOW_MINSTRASSENDIMEDIT
    {
    	"Change Min Strassen Dimension",
    	CMD_CHANGEMINSTRASSENDIM,
    	USAGE_CHANGEMINSTRASSENDIM,
    	changeMinStrassenDim,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifdef ALLOW_MINSRNUMBEREDIT
    {
    	"Change Min Stirling Number",
    	CMD_CHANGEMINSRNUMBER,
    	USAGE_CHANGEMINSRNUMBER,
    	changeMinSRNumber,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifdef ALLOW_ALGEBRAEDIT
    [SETTINGS_CHANGEALGEBRA] =
    {
        "Change Algebra",
        CMD_CHANGEALGEBRA,
        USAGE_CHANGEALGEBRA,
        changeAlgebra,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_OUTLIERCONSTEDIT
    {
    	"Change Outlier Constant",
    	CMD_CHANGEOUTLIERCONST,
    	USAGE_CHANGEOUTLIERCONST,
    	changeOutlierConst,
    	AUTOMATIC,
    	CHILD
    },
    #endif
    #ifdef ALLOW_EXITCHAREDIT
    [SETTINGS_CHANGEEXITCHAR] =
    {
        "Change Exit Char",
        CMD_CHANGEEXITCHAR,
        USAGE_CHANGEEXITCHAR,
        changeExitChar,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_RANDOMSEEDEDIT
    [SETTINGS_CHANGERANDOMSEED] =
    {
        "Change Random Seed",
        CMD_CHANGERANDOMSEED,
        USAGE_CHANGERANDOMSEED,
        changeRandomSeed,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_BOOLVARSEDIT
    [SETTINGS_CHANGEBOOLVALUES] =
    {
        "Enable-Disable Bool Settings",
        CMD_CHANGEBOOLVALUES,
        USAGE_CHANGEBOOLVALUES,
        changeBoolValues,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_MMIEDIT
    [SETTINGS_CHANGEMAXMEMOIZABLEINDICES] =
    {
        "Change Max Memoizable Indices",
        CMD_CHANGEMAXMEMIDX,
        USAGE_CHANGEMAXMEMIDX,
        changeMaxMemoizableIndices,
        AUTOMATIC,
        CHILD
    },
    [SETTINGS_EMPTYMEMOIZERSBUFFERS] =
    {
        "Flush Memoizers Buffers",
        CMD_EMPTYMEMOIZERSBUFFERS,
        USAGE_EMPTYMEMOIZERSBUFFERS,
        emptyMemoizersBuffers,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_BUFFERSEDIT
    [SETTINGS_EMPTYBUFFERS] =
    {
        "Flush Buffers",
        CMD_EMPTYBUFFERS,
        USAGE_EMPTYBUFFERS,
        emptyBuffers,
        AUTOMATIC,
        CHILD
    },
    #endif
    #ifdef ALLOW_SETDEFAULTS
    [SETTINGS_SETDEFAULTS] =
    {
        "Restore Default Settings",
        CMD_SETDEFAULTS,
        USAGE_SETDEFAULTS,
        setDefaultValues,
        AUTOMATIC,
        CHILD
    }
    #endif
};

#endif
