// programs.c 16/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h" // DA RENDERE VISIBILE       SIA AL COMPILATORE CHE AL LINKER
// #include "ExprEval/exprincl.h" // In order To use Redefined MATH_ constants


__MSSHELL_WRAPPER_ void basicCalculator(const sel_typ argc, char ** argv)
{
    (void) requires(argc ? argv[0] : NULL, argc ? NULL : "Enter an Expression", "Result is", PARSER_SHOWRESULT | PARSER_SHOWVARLIST | PARSER_SHOWDIFFTIME | PARSER_SAVERESULT);
    return ;
}

__MSSHELL_WRAPPER_ void __apnt calcolatoreAvanzato(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_ADVCALC_PROGS, adv_calc,
                        main_menu[MAIN_ADVANCEDCALCULATOR].name,
                        #if MAX_ADVCALC_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif // MAX_SETTINGS
                        );
    return ;
}

__MSSHELL_WRAPPER_ void __apnt mssManager(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_MSSMANAGER_PROGS, mss_manager,
                        main_menu[MAIN_MSSMANAGER].name,
                        #if MAX_MSSMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif // MAX_SETTINGS
                        );
    return ;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ void _MS__private __system __export operationsGroupMenu(dim_typ dim, sprog programs[static dim], const char access_point[static INFO_STRING], bool list)
{
    for( ;; )
    {

        printf("\nSelect desired Program:\n");
        PRINTL();
        // PRINTN();

        dim_typ i;
        dim_typ tmp;

        // tmp = PROGRAM_BUSY;

        CLEARBUFFER();

        if(list == BY_NUMBERS)
        {
            for(i=0; i<dim; ++i)
                printf("- %hu: %s;\n", i, programs[i].name);

            printf("- %hu: Clear SCREEN;\n", i);
            printf("- %hu: PROGRAM Informations;\n", i+1);
            printf("- %hu: Exit from PROGRAM.\n\n", i+2);
            PRINTL();


            ityp tmp2;

            while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp2 = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
                    scanf(INPUT_CONVERSION_FORMAT, &tmp2) != 1) || tmp2 != (tmp = (dim_typ)tmp2) || tmp < 0 || tmp > dim+2)
                printErr(1, "Invalid PROGRAM Mode");

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
            pulisciSchermo;
            printf2(COLOR_USER, "\nSCREEN has been correctly cleaned.\n\n");
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
                printf2(COLOR_CREDITS, "\t\t       %s\n", access_point);

                if(!readFile(tmpname))
                    printErr(2, "Non-existent DESCRIPTION File in Directory "DESCRIPTIONS_FOLDER":\n%s", low_name);

            }

            continue;
        }

        if(tmp == dim+2) break;



        sprog prog_chosen = programs[tmp];


        // PRINTING PROGRAM NAME
        char str[INFO_STRING];

        strcpy(str, prog_chosen.name);
        toupper_s(str);
        printf2(COLOR_CREDITS, str);
        PRINTN();

        if(isSett(BOOLS_SHOWDESCRIPTION))
        {
            if(prog_chosen.program_function == basicCalculator)
            {
                char low_name[FILENAME_MAX];
                strfnm(prog_chosen.name, low_name);

                if(file_exists(low_name))
                {
                    #ifdef WINOS
                        sprintf(str, "notepad %s", low_name);
                        (void) system(str);
                    #else
                        readFile(low_name);
                    #endif
                }
                else
                {
                    printErr(2, "Unable to open File:\n%s", low_name);

                    printf("\nDo you want to Reach Expression Help Page\n"); // : %s\n", EXPREVAL_TMPL);
                    printf("for more informations about EXPREVAL basilar inline functions?\n");
                    printf("[Y (Yes) / N (Not)]\n");

                    char selection;

                    do selection = toupper(getch());
                    while(selection != 'Y' && selection != 'N');

                    if(selection == 'Y')
                    {
                        sprintf(str, "START %s", EXPREVAL_TMPL);
                        (void) system(str);
                    }

                    PRINTN();
                    CLEARBUFFER();

                }
        	}
            else
            {
                char low_name[FILENAME_MAX];
                char tmpname[MAX_PATH_LENGTH] = DESCRIPTIONS_FOLDER;
                strfnm(prog_chosen.name, low_name);
                strcat(tmpname, low_name);
                printf2(COLOR_CREDITS, "\t\t       %s\n", prog_chosen.name);

                if(!readFile(tmpname))
                    printErr(2, "Non-existent DESCRIPTION File in Directory "DESCRIPTIONS_FOLDER":\n%s", low_name);

                printf2(COLOR_SYSTEM, ".\nCMDNAME: %s ; USAGE: %s .\n", prog_chosen.cmdname, prog_chosen.usage);
                PRINTL();
            }
        }

        bool rep_check;

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

            if(isSett(BOOLS_RESETLISTS)) refreshExprEvalLists();

            if((rep_check = (!prog_chosen.isFather) && isSett(BOOLS_PROGREPEATCHECK) && (!prog_chosen.automatic)))
            {
                PRINTL();
                printf2(COLOR_CREDITS, "Press any key to repeat\nor press %c to go Back to Main Menu.\n", access(curLayout)->exit_char);
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
    operationsGroupMenu(MAX_ALGEBRA_OPERATIONS,
                        alg_operations, main_menu[MAIN_ALGEBRAOPERATIONS].name,
                        #if MAX_ALGEBRA_OPERATIONS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif // MAX_ALGEBRA_OPERATIONS
                        );
    return;
}

__MSSHELL_WRAPPER_ void __apnt changeProgramSettings(const sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_SETTINGS, change_settings,
                        main_menu[MAIN_CHANGESETTINGS].name,
                        #if MAX_SETTINGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif // MAX_SETTINGS
                        );
    return;
}

//

__MSSHELL_WRAPPER_ __MSNATIVE_ void __system progInfo(sel_typ skip)
{
    if(skip)
    {
        FILE *fp = NULL;

        if((fp = fopen(DESCRIPTIONS_FOLDER"suite."DEFAULT_HELP_FILE_EXTENSION"", "r")) == NULL)
            printErr(2, "Unable to open File containing\nProgram Description");
        else
        {
        	
            printf2(COLOR_CREDITS, "\t          "PROG__NAME" V"PROG__VERSION" by ");
            printf2(access(colors)[MEMBER_COLORAUTHOR], PROG__AUTHOR);
            printf2(COLOR_CREDITS, ".\n\n\t   //________________________________________________________\\\\\
\n\nLAST UPDATE DATE: "PROG__LASTUPDATEDATE"\n");

            if(!readFile(DESCRIPTIONS_FOLDER"suite."DEFAULT_HELP_FILE_EXTENSION""))
                printErr(2, "Non-existent DESCRIPTION File in Directory "DESCRIPTIONS_FOLDER":\nsuite."DEFAULT_HELP_FILE_EXTENSION);

            PRINTL();
        }
    }
    else
    {
        printf2(COLOR_CREDITS, PROG__NAME" V"PROG__VERSION" by ");
        printf2(access(colors)[MEMBER_COLORAUTHOR], PROG__AUTHOR);
        printf2(COLOR_CREDITS, ".\n");
        PRINTL();
    }

    return;
}

__MSUTIL_ inline void __system __export freeExprEvalLists()
{
    exprValListFree(access(const_list));
    exprFuncListFree(access(func_list));
    return;
}

__MSUTIL_ void __system __export refreshExprEvalVarList(dim_typ which_env)
{
    int err;
    jmp_buf jumper;

    ityp diff = 0.00;
    exprType * const tmp = malloc(sizeof(exprType));

    errMem(tmp, VSPACE);

    tmp->var_list = INIT_VARLIST;
    ///    tmp->const_list = INIT_CONSTLIST;
    tmp->e_ANS = INIT_ANS;
    tmp->global_var = DEFAULT_GLOBALVAL;

    /* Create varlist */
    err = exprValListCreate(&(tmp->var_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Var List Creation Error\n");
        longjmp(jumper, err);
    }

    /* Init variable list */
    // err = exprValListAddAddress((*vlist), "global", &(suite.exprVars.global_var));
    err = exprValListAddAddress(tmp->var_list, "global", &(tmp->global_var));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Var List Init Error\n");
        longjmp(jumper, err);
    }

    // err = exprValListAdd((*vlist), DEFAULT_ENVS_ANSVALNAME, 0.0);
    err = exprValListAdd(tmp->var_list, DEFAULT_ENVS_ANSVALNAME, 0.00);
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Error adding variable \'%s\'\n", DEFAULT_ENVS_ANSVALNAME);
        longjmp(jumper, err);
    }

    // exprValListGetAddress((*vlist), DEFAULT_ENVS_ANSVALNAME, &(suite.exprVars.e_res));
    exprValListGetAddress(tmp->var_list, DEFAULT_ENVS_ANSVALNAME, &(tmp->e_ANS));
    if(tmp->e_ANS == NULL)
    {
        sprint("Unable to get address of \'%s\'\n", DEFAULT_ENVS_ANSVALNAME);
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
        // char c;
        char str[MAX_BUFSIZ]; // char str[MIN_STRING];

        // strcpy(str, NULL_CHAR);

        while(fgets(str, MAX_FILE_LINES, fp) != NULL)
        {

            err = exprCreate(&e, access(func_list), tmp->var_list, access(const_list), NULL, 0);
            if(err != EXPR_ERROR_NOERROR)
            {
                sprint("Expr Creation Error.\n");
                exprFree(e);
                continue;
            }


            err = exprParse(e, str);
            if(err != EXPR_ERROR_NOERROR)
            {
                exprGetErrorPosition(e, &start, &end);
                sprint("Parse Error (%d,%d).\n", start, end);
                exprFree(e);
                continue;
            }

            if(asrt)
            	gettimeofday(&tvBegin, NULL);

            ityp val;
            err = exprEval(e, &val);

            if(err != EXPR_ERROR_NOERROR)
            {
                sprint("Eval Error: %d.\n", err);
                exprFree(e);
                continue;
            }

            if(asrt)
                diff += getDiffTime(&tvBegin);

            exprFree(e);
        }

        fclose(fp);

        if(asrt)
        {
            PRINTL();
            printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, diff);
            PRINTL();
        }
    }

    listNo(which_env, ENVS)->data = tmp;

    return;
}

__MSUTIL_ void __system __export refreshExprEvalLists()
{
    freeExprEvalLists();


    access(func_list) = INIT_FUNCLIST;
    access(const_list) = INIT_CONSTLIST;
    /* Set error buffer */

    int err;
    jmp_buf jumper;

    err = setjmp(jumper);

    if(err)
        {
        /* Free stuff */

        if(access(func_list))
            exprFuncListFree(access(func_list));

        if(err != -1)
            sprint("Error: %d\n", err);

        return;
        }

    err = exprFuncListCreate(&(access(func_list)));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Func List Creation Error\n");
        longjmp(jumper, 1);
    }

    /* Init funclist */
    err = exprFuncListInit(access(func_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Error initializing internal functions\n");
        longjmp(jumper, err);
    }

    /* Create constlist */
    err = exprValListCreate(&access(const_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Const List Creation Error\n");
        longjmp(jumper, err);
    }

    /* Init constlist */
    err = exprValListInit(access(const_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Error initializing internal constants\n");
        longjmp(jumper, err);
    }

    return;
}

__MSNATIVE_ inline void __system __export setCurrentMatrix(dim_typ which_mat)
{
    if(getItemsListNo(MATRICES) != STARTING_MATNO && !extractMat(which_mat))
    {
        matrixObj * const tmp = ((matrixObj *)(listNo(which_mat, MATRICES)->data));
        matrixFree(&(tmp->matrix));
    }

    return;
}

__MSNATIVE_ inline void __system __export setCurrentLog(dim_typ which_log)
{
	logObj * const tmp = malloc(sizeof(logObj));
    errMem(tmp, VSPACE);
    tmp->buffer = malloc(sizeof(char)*DEFAULT_BUFSIZE);
    errMem(tmp->buffer, VSPACE);
    strcpy(tmp->buffer, NULL_CHAR); // initializing log buffer
	tmp->buflen = DEFAULT_BUFSIZE;
    listNo(which_log, LOGS)->data = tmp;
    return;
}


__MSSHELL_WRAPPER_ __MSNATIVE_ void __system setDefaults()
{

    static bool once_executed = false;
    access(mode) = PROGRAM_BUSY;

    if(once_executed)
        resetProgramSettings(access(curLayout), listNo(access(lists)[LAYOUTS].cur_item, LAYOUTS)->path);
    else
    	once_executed = true;
    	
    randomize;
    access(random_seed) = starting_random_seed;

    dim_typ i;
    refreshExprEvalLists();

    if(access(lmpMatrix))
    {
        if(access(lmpMatrix)->matrix)
            matrixFree(&(access(lmpMatrix)->matrix));
        free(access(lmpMatrix));
    }
    
    access(lmpMatrix) = malloc(sizeof(matrixObj));
    errMem(access(lmpMatrix), VSPACE);

	resetLmpMatrix();

    return;
}
