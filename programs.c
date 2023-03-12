// programs.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h" // DA RENDERE VISIBILE       SIA AL COMPILATORE CHE AL LINKER
// #include "ExprEval/exprincl.h" // In order To use Redefined MATH_ constants

__MSSHELL_WRAPPER_ void __apnt multiUser(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MULTIUSER_LOGIN, multiuser_prog, main_menu[MAIN_MULTIUSER].name, BY_CHARS); // MAX_ADVCALC_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return ;
}

__MSSHELL_WRAPPER_ void __apnt advancedCalculator(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_ADVCALC_PROGS, adv_calc, main_menu[MAIN_ADVANCEDCALCULATOR].name, BY_CHARS); // MAX_ADVCALC_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return ;
}

__MSSHELL_WRAPPER_ void basicCalculator(const sel_typ argc, char ** argv)
{
	mpfr_t tmp;
 (void) requires(tmp, argc ? argv[0] : NULL, argc ? NULL : "Enter an Expression", "Result is", PARSER_SHOWRESULT | PARSER_SHOWVARLIST | PARSER_SHOWDIFFTIME | PARSER_SAVERESULT);
    mpfr_clear(tmp);
	return ;
}

__MSSHELL_WRAPPER_ void __apnt mssManager(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_MSSMANAGER_PROGS, mss_manager, main_menu[MAIN_MSSMANAGER].name, BY_CHARS); // MAX_MSSMANAGER_PROGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return ;
}

__MSSHELL_WRAPPER_ __MATHSUITE void    operationsGroupMenu(dim_typ dim, sprog programs[static dim], const char access_point[static INFO_STRING], bool list)
{
    for( ; ; )
    {
		dim_typ i;
        dim_typ tmp;
        
        printf("\nSelect desired Program:\n");
        PRINTL();
        CLEARBUFFER();

        if(list == BY_NUMBERS)
        {
            for(i=0; i<dim; ++i)
                printf("- %hu: %s;\n", i, programs[i].name);

            printf("- %hu: Clear SCREEN;\n", i);
            printf("- %hu: PROGRAM Informations;\n", i+1);
            printf("- %hu: Exit from PROGRAM.\n\n", i+2);
            PRINTL();


            mpfr_t tmp2;
			while(requires(tmp2, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp2) || mpfr_cmp_ui(tmp2, (tmp = mpfr_get_ui(tmp2, MPFR_RNDN))) || tmp < 0 || tmp > dim+2)
            {
            	mpfr_clear(tmp2);
				printErr(1, "Invalid PROGRAM Mode");
			}
            mpfr_clear(tmp2);

        }
        else
        {
            // printf(programs[MAIN_ALGEBRAOPERATIONS].name);
            for(i=0; i<dim; ++i)
                printf("- %c: %s;\n", i+'A', programs[i].name);
            printf("- %c: Clear SCREEN;\n", i+'A');
            printf("- %c: PROGRAM Informations;\n", i+1+'A');
            printf("- %c: Exit from PROGRAM.\n\n", i+2+'A');
            PRINTL();

            sel_typ tmp2;

            do
                tmp2 = toupper(getch());
            while(tmp2 < 'A' || tmp2 > dim+2+'A');

            // tmp = tmp2-65;
            // tmp = (toupper(tmp2)-65);

            tmp = tmp2-'A';
        }

        __pmode__ = tmp;
        CLEARBUFFER();
        PRINT2N();


        if(tmp == dim)
        {
        	clearScreen();
            msprintf(COLOR_USER, "\nSCREEN has been correctly cleaned.\n\n");
            PRINTL();
            operationsGroupMenu(dim, programs, access_point, list);
            return;
        }

        if(tmp == dim+1)
        {
            if(!strcmp(access_point, NULL_CHAR))
                progInfo(WITH_DESCRIPTION);
            else
            {
                char low_name[FILENAME_MAX] = NULL_CHAR;
                char tmpname[MAX_PATH_LENGTH] = DESCRIPTIONS_FOLDER;
                strfnm(access_point, low_name);
                strcat(tmpname, low_name);
                PRINT2N();
                msprintf(COLOR_CREDITS, "\t\t       %s\n", access_point);

                if(!readFile(tmpname))
                    printErr(2, "Non-existing DESCRIPTION File in Directory "DESCRIPTIONS_FOLDER":\n%s", low_name);

            }

            continue;
        }

        if(tmp == dim+2) break;


        sprog prog_chosen = programs[tmp];

        // PRINTING PROGRAM NAME
        char str[INFO_STRING];

        strcpy(str, prog_chosen.name);
        toupper_s(str);
        msprintf(COLOR_CREDITS, str);
        PRINTN();

        bool rep_check;

		if(accessCurrentUser->level < prog_chosen.level)
			printErr(ERROR_PERMISSIONDENIED, "You haven't permission to execute: \"%s\" program\nRequired Level: %s", prog_chosen.name, suite_c.usertypes_names[prog_chosen.level]);
		else
	        do
	        {
	            // printf("\n_____________________________________________________\n\n");
	
	            struct timeval tvBegin;
	            const bool asrt = isSett(BOOLS_SHOWEXECTIME);
	
	            CLEARBUFFER();
	
	            if(asrt)
	            	gettimeofday(&tvBegin, NULL);
	
	
	            prog_chosen.program_function(0, NULL); // RICHIAMA LA FUNZIONE O METODO DEL PROGRAMMA SELEZIONATO
	            // avendo a disposizione l'indirizzo della funzione corrispondente al subprogram scelto.
	            if(asrt)
	                printf("\nExecution Average Time: %.*f.\n\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	
				/// CRITICAL POINT
	            // refreshExprEvalLists();
	
	            if((rep_check = (!prog_chosen.isFather) && isSett(BOOLS_PROGREPEATCHECK) && (!prog_chosen.automatic)))
	            {
	                PRINTL();
	                msprintf(COLOR_CREDITS, "Press any key to repeat\nor press %c to go Back to Main Menu.\n", access(curLayout)->exit_char);
	            }
	            CLEARBUFFER();
	        }
	        while(rep_check && getch() != access(curLayout)->exit_char);
    }
    
    CLEARBUFFER();
    PRINTL();

    // #undef prog_chosen

    return;
}

// New Families Access Points
//
__MSSHELL_WRAPPER_ void __apnt algebraOperations(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_ALGEBRA_OPERATIONS, alg_operations, main_menu[MAIN_ALGEBRAOPERATIONS].name, BY_NUMBERS); // MAX_ALGEBRA_OPERATIONS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return;
}

__MSSHELL_WRAPPER_ void __apnt changeProgramSettings(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_SETTINGS, change_settings, main_menu[MAIN_CHANGESETTINGS].name, BY_CHARS); // MAX_SETTINGS+MAX_OMNIPRESENT_ELEMENTS<MAX_CASEINSENSITIVE_CHARS_ALPHABET);
    return;
}

__MSSHELL_WRAPPER_ __MATHSUITE void  progInfo(sel_typ skip)
{
	PRINTL();
	msprintf(COLOR_CREDITS, "\t        "PROGRAM_NAME" V"PROGRAM_VERSION" by ");
    msprintf(access(colors)[MEMBER_COLORAUTHOR], PROGRAM_AUTHOR);
    msprintf(COLOR_CREDITS, " and ");
    msprintf(access(colors)[MEMBER_COLORAUTHOR], PROGRAM_COAUTHOR);
    msprintf(COLOR_PURPLE, "\n\t     Powerful Calculus Environment and Matrices Handling Engine.\n\n");
    PRINTL();
    
    if(skip)
    {
    	/*
        FILE *fp = NULL;

        if((fp = fopen(DESCRIPTIONS_FOLDER"mathSuite."DEFAULT_HELP_FILE_EXTENSION"", "r")) == NULL)
            printErr(2, "Unable to open File containing\nProgram Description");
        else
        {
        */
        msprintf(COLOR_CREDITS, "LAST UPDATE DATE: "PROGRAM_LASTUPDATEDATE"\n");
		msprintf(COLOR_USER, "\nThis software is licensed under Creative Commons CC BY-SA 2.0\n");
		msprintf(COLOR_USER, "contact me for more informations: marco_chiarelli@yahoo.it\n");
		msprintf(COLOR_USER, "or marcochiarelli.nextgenlab@gmail.com. For further informations about\n");
		msprintf(COLOR_USER, "ExprEval or MSCenv mathematical set of functions, see basicCalculator.txt\n");
		msprintf(COLOR_USER, "or visit the project page at: https://sourceforge.net/projects/mathsuite/\n");
		
		/*
        if(!readFile(DESCRIPTIONS_FOLDER"mathSuite."DEFAULT_HELP_FILE_EXTENSION""))
            printErr(2, "Non-existing DESCRIPTION File in Directory "DESCRIPTIONS_FOLDER":\nsuite."DEFAULT_HELP_FILE_EXTENSION);
        */

        PRINTL();
        // }
    }

    return;
}

 inline void   freeExprEvalLists()
{
    exprValListFree(access(const_list));
    exprFuncListFree(&access(func_list));
    return;
}

 void   refreshExprEvalVarList(dim_typ which_env)
{
    EXPRERRTYPE err;
    jmp_buf jumper;

    ityp diff = 0.00;
    exprType * const tmp = malloc(sizeof(exprType));

    errMem(tmp, VSPACE);

    tmp->var_list = INIT_VARLIST;
    ///    tmp->const_list = INIT_CONSTLIST;
    tmp->e_ANS = INIT_ANS;
    mpfr_init_set_d(tmp->global_var, DEFAULT_GLOBALVAL, MPFR_RNDN);

    /* Create varlist */
    err = exprValListCreate(&(tmp->var_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        msyprintf(COLOR_SYSTEM, "Var List Creation Error\n");
        longjmp(jumper, err);
    }

    /* Init variable list */
    // err = exprValListAddAddress((*vlist), "global", &(suite.exprVars.global_var));
    err = exprValListAddAddress(tmp->var_list, "global", &(tmp->global_var));
    if(err != EXPR_ERROR_NOERROR)
    {
        msyprintf(COLOR_SYSTEM, "Var List Init Error\n");
        longjmp(jumper, err);
    }

    // err = exprValListAdd((*vlist), DEFAULT_ENVS_ANSVALNAME, 0.0);
    mpfr_t mpftmp;
    mpfr_init_set_ui(mpftmp, 0, MPFR_RNDN);
    err = exprValListAdd(tmp->var_list, DEFAULT_ENVS_ANSVALNAME, mpftmp);
    if(err != EXPR_ERROR_NOERROR)
    {
        msyprintf(COLOR_SYSTEM, "Error adding variable \'%s\'\n", DEFAULT_ENVS_ANSVALNAME);
        mpfr_clear(mpftmp); 
        longjmp(jumper, err);
    }
    
    mpfr_clear(mpftmp); 
    // exprValListGetAddress((*vlist), DEFAULT_ENVS_ANSVALNAME, &(suite.exprVars.e_res));
    exprValListGetAddress(tmp->var_list, DEFAULT_ENVS_ANSVALNAME, &(tmp->e_ANS));
    if(tmp->e_ANS == NULL)
    {
        msyprintf(COLOR_SYSTEM, "Unable to get address of \'%s\'\n", DEFAULT_ENVS_ANSVALNAME);
        longjmp(jumper, EXPR_ERROR_UNKNOWN);
    }

    FILE *fp = NULL;
    const bool assert = getItemsListNo(ENVS) != STARTING_ENVSNO;
    const bool asrt = isSett(BOOLS_SHOWDIFFTIME);

    if(assert && (fp = checkForFHErrors(listNo(which_env, ENVS)->path, "r")) == NULL)
        return;

    exprObj *e = INIT_OBJLIST;
    int start, end;
    struct timeval tvBegin;

    if(assert)
    {
        char str[MAX_BUFSIZ];
        
    	while(fgets(str, MAX_FILE_LINES, fp) != NULL)
        {
            err = exprCreate(&e, access(func_list), tmp->var_list, access(const_list), NULL, 0);
            if(err != EXPR_ERROR_NOERROR)
            {
                msyprintf(COLOR_SYSTEM, "Expr Creation Error.\n");
                exprFree(e);
                continue;
            }

            err = exprParse(e, str);
            if(err != EXPR_ERROR_NOERROR)
            {
                exprGetErrorPosition(e, &start, &end);
                msyprintf(COLOR_SYSTEM, "Parse Error (%d,%d).\n", start, end);
                exprFree(e);
                continue;
            }

            if(asrt)
            	gettimeofday(&tvBegin, NULL);

            mpfr_t val;
       	 	err = exprEval(e, val);  
				      
            if(err != EXPR_ERROR_NOERROR)
            {
        		msyprintf(COLOR_ERROR, "Eval Error:\n%s.\n", suite_c.exprErrorString[err == EXPR_ERROR_UNKNOWN ? 0 : err+1  ]);
                exprFree(e);
                continue;
            }

            if(asrt)
                diff += getDiffTime(&tvBegin);

            exprFree(e);
        	mpfr_clear(val);
        }
        
        fclose(fp);
        
        #ifndef __DISABLE_DATABASE
        if(!access(curLayout)->database.con)
        	dbEstablishConnection(0, NULL); // (char * [4]){NULL, NULL, NULL, NULL});
    
    	if(access(curLayout)->database.con)
		{
			char tab[MAX_PATH_LENGTH];
			strcpy(tab, listNo(which_env, ENVS)->path);
			(void) strtok(tab, EXTENSION_DOT);
			dim_typ varlistno;
			char var[MAX_VARS][VAR_MAXNAMELENGTH] = { NULL_CHAR };
			char val[MAX_VARS][VAR_MAXVALLENGTH] = { NULL_CHAR };
			readVarList(tab, &varlistno, var, val);
			
			for(dim_typ i=0; i<varlistno; ++i)
			{
				sprintf(var[i], "%s=%s;", var[i], val[i]);
				err = exprCreate(&e, access(func_list), tmp->var_list, access(const_list), NULL, 0);
	            if(err != EXPR_ERROR_NOERROR)
	            {
	                msyprintf(COLOR_SYSTEM, "Expr Creation Error.\n");
	                exprFree(e);
	                continue;
	            }
	
	            err = exprParse(e, var[i]);
	            if(err != EXPR_ERROR_NOERROR)
	            {
	                exprGetErrorPosition(e, &start, &end);
	                msyprintf(COLOR_SYSTEM, "Parse Error (%d,%d).\n", start, end);
	                exprFree(e);
	                continue;
	            }
	
	            if(asrt)
	            	gettimeofday(&tvBegin, NULL);
	
	            mpfr_t val;
	       	 	err = exprEval(e, val);  
					      
	            if(err != EXPR_ERROR_NOERROR)
	            {
        			msyprintf(COLOR_ERROR, "Eval Error:\n%s.\n", suite_c.exprErrorString[err == EXPR_ERROR_UNKNOWN ? 0 : err+1  ]);
	                exprFree(e);
	                continue;
	            }
	
	            if(asrt)
	                diff += getDiffTime(&tvBegin);
	
	            exprFree(e);
	        	mpfr_clear(val);
			}
		} 
    	#endif
		
        if(asrt)
        {
            PRINTL();
            msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, diff);
            PRINTL();
        }
    }

    listNo(which_env, ENVS)->data = tmp;
    return;
}

 void   refreshExprEvalLists()
{
    freeExprEvalLists();
    access(func_list) = INIT_FUNCLIST;
    access(const_list) = INIT_CONSTLIST;
    /* Set error buffer */

    EXPRERRTYPE err;
    jmp_buf jumper;

    err = setjmp(jumper);

    if(err)
        {

        if(err != -1)
            msyprintf(COLOR_SYSTEM, "Error: %d\n", err);

        return;
        }

    exprFuncListCreate(&access(func_list));

    /* Init funclist */
    err = exprFuncListInit(access(func_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        msyprintf(COLOR_SYSTEM, "Error initializing internal functions\n");
        longjmp(jumper, err);
    }

    /* Create constlist */
    err = exprValListCreate(&access(const_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        msyprintf(COLOR_SYSTEM, "Const List Creation Error\n");
        longjmp(jumper, err);
    }

    /* Init constlist */
    err = exprValListInit(access(const_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        msyprintf(COLOR_SYSTEM, "Error initializing internal constants\n");
        longjmp(jumper, err);
    }

    return;
}

__MATHSUITE inline void   setCurrentMatrix(dim_typ which_mat)
{
    if(getItemsListNo(MATRICES) != STARTING_MATNO && !extractMat(which_mat))
    {
        matrixObj * const tmp = ((matrixObj *)(listNo(which_mat, MATRICES)->data));
        matrixFree(&(tmp->matrix), tmp->dim);
    }

    return;
}

__MATHSUITE inline void   setCurrentLog(dim_typ which_log)
{
	logObj * const tmp = malloc(sizeof(logObj));
    errMem(tmp, VSPACE);
    tmp->buffer = calloc(DEFAULT_BUFSIZE, sizeof(char));
    errMem(tmp->buffer, VSPACE);
    strcpy(tmp->buffer, NULL_CHAR); // initializing log buffer
    errMem(tmp->buffer, VSPACE);
	tmp->buflen = DEFAULT_BUFSIZE;
    listNo(which_log, LOGS)->data = tmp;
    return;
}

__MSSHELL_WRAPPER_ __MATHSUITE void  setDefaults()
{

    static bool once_executed = false;
    access(mode) = PROGRAM_BUSY;

    if(once_executed)
        resetProgramSettings(access(curLayout), listNo(access(lists)[LAYOUTS].cur_item, LAYOUTS)->path);
    else
    	once_executed = true;
    	
    randomize();
    access(random_seed) = starting_random_seed;

    refreshExprEvalLists();
    
    matrixObj * const currentLmpMatrix = accessCurrentLmpMatrix;

    if(currentLmpMatrix)
    {
        if(currentLmpMatrix->matrix)
            matrixFree(&(currentLmpMatrix->matrix), currentLmpMatrix->dim); 
        free(currentLmpMatrix);
    }
    
    resetLmpMatrix();
    return;
}
