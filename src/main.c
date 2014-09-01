/*!////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
//!////////////////////////////////////////////////////////////////////////////////////////////////////////////// ///
/*!________________________________________________________________________________________________________________*/
///                         mathSuite v5.75 --- by Marco Chiarelli aka DekraN aka Wesker013                  	  ///
///        							LAST UPDATE: 15:00 of 01/09/2014 by Myself. 								  ///
/// 	This program is protected by Creative Commons CC BY-SA 2.0 License. For more informations contact me. 	  ///
///                     You can contact me at: marco_chiarelli@yahoo.it or on the secundary mail:                 ///
/// marcochiarelli.nextgenlab@gmail.com in order to report a bug or simply for sending me an advice that could be ///
///                        useful or could improve the speed or optimize my MSCenv System.                        ///
//!-------------------------------------------------------------------------------------------------------------- ///
//!______________________________________________________________________________________________________________ ///
//!-------------------------------------------------------------------------------------------------------------- ///
///                contact me at: marco_chiarelli@yahoo.it or marcochiarelli.nextgenlab@gmail.com                 ///
///        I'll be glad to fix your scripts or simply to take away your doubts about the program                  ///
//!-------------------------------------------------------------------------------------------------------------- ///
//!-------------------------------------------------------------------------------------------------------------- ///
/// Thanks to giggikr: http://forum.html.it/forum/showthread/t-1374455.html for his own function, cambiabase,     ///
///                         which I renamed, modified and adapted to this program. Thanks to:                     ///
/// http://elite.polito.it/files/courses/12BHD/progr/Esercizi-C-v2_01.pdf for some of their scripts.              ///
/// Thanks to Bibek Subedi, for his invertMatrix function, which I renamed, modified and adapted to this program. ///
/// Link Source: http://programming-technique.blogspot.it/2011/09/numerical-methods-inverse-of-nxn-matrix.html    ///
/// Massive thanks to Brian Allen Vanderburg II for his fabulous C parser and inline functions solver, EXPREVAL,  ///
/// which elegantly gave in theory infinite functionalities and potential to my program. That's the project link  ///
/// with Online Documentation: http://expreval.sourceforge.net/ Thanks to: http://www.cprogramming.com/tips/ and  ///
///             http://stackoverflow.com/questions/599365/what-is-your-favorite-c-programming-trick               ///
///                       http://stackoverflow.com/questions/132241/hidden-features-of-c                          ///
///    that are some websites in which I found a lot of useful C tips and tricks, and they were an important      ///
///  checkpoint for resources retrieving in order to speed-up and optimize my program. Still greatly thanks to    ///
///    Bibek Subedi for his website: http://www.programming-techniques.com/ which put in front of my eyes a new   ///
/// world of C programming. I also recently renewed the program code by improving a lot of his C tricks and tips. ///
/// For example, the upper-triangular Matrixes conversion, which was useful to enhance some functions like det(), ///
///      							sgeqsolver ExprEval inline command, etc.									  ///
///                 Greatly thanks to Daniel Veillard for his fabulous XML Parser, LIBXML2.                       ///
/// 		Thanks to Francesco Palma for reporting me some bugs, and finally, massive thanks to my				  ///
/// Informatic Fundaments Teacher, Mario Alessandro Bochicchio, which gave me a lot of C advices and some general ///
///            tricks and tips, that enlarged my professional informatic horizonts. That's all...                 ///
/*!________________________________________________________________________________________________________________*/
//!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/


#include "dutils.h" // DA RENDERE VISIBILE SIA AL COMPILATORE CHE AL LINKER

#ifdef STACKALLOC
    struct program suite =
    {
        SYSTEM_LOG,
        INIT_MATRIXSUMFUNC,
        INIT_MATRIXPRODFUNC,
        INIT_MATRIXKPRODFUNC,
        INIT_EXPRTYPE,
        INIT_MATRIXOBJ,
        INIT_MATRIXOBJ,
        INIT_LOGOBJ,
        INIT_LOGOBJ,
        INIT_LAYOUTOBJ,
        {
            {
                NULL,
                0
            },
            {
                NULL,
                0
            },
            {
            	NULL,
            	0
            },
            {
            	NULL,
            	0
            }
        },
        {
            {
                {
                    INIT_VARLIST,
                    INIT_VARLIST
                },
                MAIN_ENV,
                STARTING_ENVSNO
            },
            {
                {
                    INIT_MATLIST,
                    INIT_MATLIST
                },
                MAIN_MAT,
                STARTING_MATNO
            },
            {
                {
                    INIT_LOGSLIST,
                    INIT_LOGSLIST
                },
                MAIN_LOG,
                STARTING_LOGSNO
            },
            {
                {
                    INIT_LAYOUTLIST,
                    INIT_LAYOUTLIST
                },
                MAIN_LAYOUT,
                STARTING_LAYOUTNO
            }
        },
        INIT_FUNCLIST,
        INIT_CONSTLIST,
        DEFAULT_COLORS_PATH,
        INITIALIZING_RANDOM_SEED,
        {
            INIT_COLOR,
            INIT_COLOR,
            INIT_COLOR,
            INIT_COLOR,
            INIT_COLOR,
            INIT_COLOR
        },
        PROGRAM_BUSY,
        INVALID_EXITHANDLE,
        false,
        false
    };
#else
    struct program * suite;
#endif

const struct prog_constants suite_c =
{
    {
        "Varlist",
        "Matrix",
        "Log",
        "Layout"
    },
    #ifdef XMLCALL
	    {
	    	STRING_FALSE,
	    	STRING_TRUE
	    },
	#endif
    {
        "Show Programs Descriptions",
        "Set Item as Current at Creation",
        "Items Selection By Path",
        "Lists Items Autosaving",
        #if WINOS
            "System Log Security Check",
        #endif
        "System INPUT Parser",
        "Base Calculator Parser",
        "Advanced Calculator Parser",
        "Matrices Items Inserting Parser",
        "Reset Consts List and Functions List",
        "Save last Results",
        "Show Varlist after Operations",
        "Show Average Time after Operations",
        "Show Average Time after Programs Executions",
        "Show Date and Time in Textline",
        "Domains Checks",
        "Matrices Chars Insert",
        "Matrices Automatic Inserting",
        "Lazy Execution",
        "Inverse Operations",
        "Inverse Operations Automatic Deactivator",
        "Degrees Entering on Trigonometric Operations",
        "Request Programs Repetitions",
        "Extensive MultiThreading"
    },
    {
        {
            BITMASK_SHOWDESCRIPTION,
            DEFAULT_SHOWDESCRIPTION
        },
        {
            BITMASK_AUTOSETCURITEM,
            DEFAULT_AUTOSETCURITEM
        },
        {
            BITMASK_ITEMSSELECTBYPATH,
            DEFAULT_ITEMSSELECTBYPATH
        },
        {
            BITMASK_ITEMSAUTOSAVING,
            DEFAULT_ITEMSAUTOSAVING
        },
        #if WINOS
        {
            BITMASK_SYSLOGSECURITYCHECK,
            DEFAULT_SYSLOGSECURITYCHECK
        },
        #endif // WINOS
        {
            BITMASK_SYSTEMPARSING,
            DEFAULT_SYSTEMPARSING
        },
        {
            BITMASK_BASECALCPARSING,
            DEFAULT_BASECALCPARSING
        },
        {
            BITMASK_ADVCALCPARSING,
            DEFAULT_ADVCALCPARSING
        },
        {
            BITMASK_MATRIXPARSING,
            DEFAULT_MATRIXPARSING
        },
        {
            BITMASK_RESETLISTS,
            DEFAULT_RESETLISTS
        },
        {
            BITMASK_SAVERESULTS,
            DEFAULT_SAVERESULTS
        },
        {
            BITMASK_SHOWVARLIST,
            DEFAULT_SHOWVARLIST
        },
        {
            BITMASK_SHOWDIFFTIME,
            DEFAULT_SHOWDIFFTIME
        },
        {
            BITMASK_SHOWEXECTIME,
            DEFAULT_SHOWEXECTIME
        },
        {
            BITMASK_PRINTTIME,
            DEFAULT_PRINTTIME
        },
        {
            BITMASK_DOMAINCHECK,
            DEFAULT_DOMAINCHECK
        },
        {
            BITMASK_CHARINSERT,
            DEFAULT_CHARINSERT
        },
        {
            BITMASK_INSERTMODE,
            DEFAULT_INSERTMODE
        },
        {
            BITMASK_LAZYEXECUTION,
            DEFAULT_LAZYEXECUTION
        },
        {
            BITMASK_INVERSEOPERATIONS,
            DEFAULT_INVERSEOPERATIONS
        },
        {
            BITMASK_AUTOTURNBACK,
            DEFAULT_AUTOTURNBACK
        },
        {
            BITMASK_DEGREESENTERING,
            DEFAULT_DEGREESENTERING
        },
        {
            BITMASK_PROGREPEATCHECK,
            DEFAULT_PROGREPEATCHECK
        },
        {
        	BITMASK_EXTENSIVEMULTITHREADING,
        	DEFAULT_EXTENSIVEMULTITHREADING
        }
    },
    {
    	MAX_FIBONACCI_MEMOIZABLE_INDEX,
    	MAX_FATTORIALE_MEMOIZABLE_INDEX,
    	MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX,
    	MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX
    },
    {
    	"Fibonacci",
    	"Factorial",
		"Even SemiFactorial",
		"Odd SemiFactorial",
		"Every"
    },
    {
        "Real Numbers",
        "Complex Numbers",
        "Quaternions",
        "Octonions",
        "Sedenions"
    },
    {
        {
            REAL_UNIT_NAME
        },
        {
            REAL_UNIT_NAME,
            "i"
        },
        {
            REAL_UNIT_NAME,
            "i",
            "j",
            "k"
        },
        {
            REAL_UNIT_NAME,
            "e1",
            "e2",
            "e3",
            "e4",
            "e5",
            "e6",
            "e7"
        },
        {
            REAL_UNIT_NAME,
            "e1",
            "e2",
            "e3",
            "e4",
            "e5",
            "e6",
            "e7",
            "e8",
            "e9",
            "e10",
            "e11",
            "e12",
            "e13",
            "e14",
            "e15"
        }
    },
    {
        _matrixSum,
        _matrixCSum,
        _matrixQSum,
        _matrixOSum,
        _matrixSSum
    },
    {
        _matrixProduct,
        _matrixCProduct,
        _matrixQProduct,
        _matrixOProduct,
        _matrixSProduct
    },
    {
        _matrixKProduct,
        _matrixKCProduct,
        _matrixKQProduct,
        _matrixKOProduct,
        _matrixKSProduct
    },
    {
        "Default                COLOR",
        "Errors                 COLOR",
        "Credits and Directives COLOR",
        "User Executions        COLOR",
        "System                 COLOR",
        "Signatures/Decorations COLOR"
    },
    {
        "BLACK",
        "BLUE",
        "GREEN",
        "AZURE",
        "DARK RED",
        "PURPLE",
        "YELLOW",
        "WHITE",
        "GREY",
        "LUM BLUE",
        "LUM GREEN",
        "CRYSTAL",
        "LUM RED",
        "LUM PURPLE",
        "LUM YELLOW",
        "WHITEPLUS",
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        "LUM WHITE",
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        "LUM AZURE",
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING,
        UNKNOWN_STRING
    }
};

sprog main_menu[MAX_PROGRAMMI] =
{
    [MAIN_BASECALCULATOR] =
    {
        "Base Calculator",
        CMD_BASECALC,
        USAGE_BASECALC,
        calcolatoreDiBase,
        BY_USER,
        CHILD
    },
    [MAIN_ADVANCEDCALCULATOR] =
    {
        "Advanced Calculator",
        CMD_ADVCALC,
        USAGE_ADVCALC,
        calcolatoreAvanzato,
        AUTOMATIC,
        FATHER
    },
    [MAIN_ALGEBRAOPERATIONS] =
    {
        "Linear Algebra Operations",
        CMD_LINEARALGEBRAOPERATIONS,
        USAGE_LINEARALGEBRAOPERATIONS,
        algebraOperations,
        AUTOMATIC,
        FATHER
    },
    #ifdef ALLOW_MSSMANAGER
    [MAIN_MSSMANAGER] =
    {
        "PROGRAM Scripts Manager",
        CMD_MSSMANAGER,
        USAGE_MSSMANAGER,
        mssManager,
        AUTOMATIC,
        FATHER
    },
    #endif
    [MAIN_CHANGESETTINGS] =
    {
        "Edit PROGRAM Settings",
        CMD_CHANGESETTINGS,
        USAGE_CHANGESETTINGS,
        changeProgramSettings,
        AUTOMATIC,
        FATHER
    }
};


int main(int argc, char **argv)
{
    #ifndef STACKALLOC
        suite = malloc(sizeof(struct program));
        errMem(suite, HEAPALLOC_ERROR);

        // INITIALIZING suite vars
        strcpy(access(sysLogPath), SYSTEM_LOG);

        access(matrixSumFunc) = INIT_MATRIXSUMFUNC;
        access(matrixProdFunc) = INIT_MATRIXPRODFUNC;
        access(matrixKProdFunc) = INIT_MATRIXKPRODFUNC;

        access(exprVars) = INIT_EXPRTYPE;
        access(curMatrix) = INIT_MATRIXOBJ;
        access(lmpMatrix) = INIT_MATRIXOBJ;
        access(curLog) = INIT_LOGOBJ;
        access(sysLog) = INIT_LOGOBJ;
        access(curLayout) = INIT_LAYOUTOBJ;

        #pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
        for(dim_typ i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
        {
        	access(sysMem)[i].memoizer = NULL;
        	access(sysMem)[i].current_max_index = 0;
        }

        access(lists)[ENVS].items_list[NEXT_LISTNODE] = INIT_VARLIST;
        access(lists)[ENVS].items_list[PREV_LISTNODE] = INIT_VARLIST;
        access(lists)[ENVS].cur_item = MAIN_ENV;
        access(lists)[ENVS].itemsno = STARTING_ENVSNO;

        access(lists)[MATRICES].items_list[NEXT_LISTNODE] = INIT_MATLIST;
        access(lists)[MATRICES].items_list[PREV_LISTNODE] = INIT_MATLIST;
        access(lists)[MATRICES].cur_item = MAIN_MAT;
        access(lists)[MATRICES].itemsno = STARTING_MATNO;

        access(lists)[LOGS].items_list[NEXT_LISTNODE] = INIT_LOGSLIST;
        access(lists)[LOGS].items_list[PREV_LISTNODE] = INIT_LOGSLIST;
        access(lists)[LOGS].cur_item = MAIN_LOG;
        access(lists)[LOGS].itemsno = STARTING_LOGSNO;

        access(lists)[LAYOUTS].items_list[NEXT_LISTNODE] = INIT_LAYOUTLIST;
        access(lists)[LAYOUTS].items_list[PREV_LISTNODE] = INIT_LAYOUTLIST;
        access(lists)[LAYOUTS].cur_item = MAIN_LAYOUT;
        access(lists)[LAYOUTS].itemsno = STARTING_LAYOUTNO;

        access(func_list) = INIT_FUNCLIST;
        access(const_list) = INIT_CONSTLIST;

        strcpy(access(colors_path), COLORS_PATH);

        access(random_seed) = INITIALIZING_RANDOM_SEED;

        __pmode__ = PROGRAM_BUSY;

        access(exitHandle) = INVALID_EXITHANDLE;
        access(mss) = access(sigresult) = false;

    #endif

	xmlDocPtr doc = NULL;       /* document pointer */
	xmlNodePtr root_node = NULL, node = NULL, node1 = NULL;/* node pointers */
	xmlDtdPtr dtd = NULL;       /* DTD pointer */

    if(!file_exists(access(colors_path)))
    {

    	xmlInitParser();
    	LIBXML_TEST_VERSION;

    	doc = xmlNewDoc(BAD_CAST "1.0");
    	root_node = xmlNewNode(NULL, BAD_CAST "colors ");

    	xmlDocSetRootElement(doc, root_node);

    	char tmp_string[INFO_STRING] = NULL_CHAR;
		sprintf(tmp_string, "%hu", DEFAULT_COLOR);
		xmlNewChild(root_node, NULL, BAD_CAST "defaultColor", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", _COLOR_ERROR);
		xmlNewChild(root_node, NULL, BAD_CAST "errorsColor", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", _COLOR_CREDITS);
		xmlNewChild(root_node, NULL, BAD_CAST "creditsColor", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", _COLOR_USER);
		xmlNewChild(root_node, NULL, BAD_CAST "userColor", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", _COLOR_SYSTEM);
		xmlNewChild(root_node, NULL, BAD_CAST "systemColor", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", _COLOR_AUTHOR);
		xmlNewChild(root_node, NULL, BAD_CAST "authorColor", BAD_CAST tmp_string);

		xmlSaveFormatFileEnc(access(colors_path), doc, XML_ENCODING, 1);
	    /// xmlFreeNode(root_node);
	    xmlFreeDoc(doc);
	    xmlCleanupParser();
    }

	#if WINOS
		_colFileLoader(access(colors_path));
	    if(!ShowWindow(GetConsoleWindowNT(), SW_MAXIMIZE))
	       printErr(22, "ShowWindow SW_MAXIMIZE failed with error: %lu", GetLastError());
	    SetDefaultColor();
	#endif

    /// SetExitButtonState(DISABLED); // testing handler

    setDefaults();

    char path[MAX_PATH_LENGTH];
    char * extension = NULL;

    if(argv[1] && strlen(argv[1]) > 1 && ((extension = strrchr(argv[1], '.')+1) == NULL || (!strcmp(extension, DEFAULT_PATHLIST_FILE_EXTENSION))))
        strcpy(path, argv[1]);
    else
        sprintf(path, "%s.%s", DEFAULT_MAIN_PATH, DEFAULT_PATHLIST_FILE_EXTENSION);

    dim_typ i;

    const sel_typ mode = extension ? checkItemTypeByExtension(extension) : MAX_LISTS;

    _lfLoader(path);

    if((mode != LAYOUTS && !getItemsListNo(LAYOUTS)) || (mode == LAYOUTS && !file_exists(argv[1])))
    {

    	xmlInitParser();
    	LIBXML_TEST_VERSION;

    	doc = xmlNewDoc(BAD_CAST "1.0");
    	root_node = xmlNewNode(NULL, BAD_CAST "settings");

    	xmlDocSetRootElement(doc, root_node);

    	node = xmlNewNode(NULL, BAD_CAST "generalSettings");

    	xmlNewChild(node, NULL, BAD_CAST "exitChar", BAD_CAST EXIT_CHAR);
    	char tmp_string[INFO_STRING] = NULL_CHAR;
    	sprintf(tmp_string, "%hu", DEFAULT_PRECISION);
    	xmlNewChild(node, NULL, BAD_CAST "programPrecision", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", DEFAULT_STABILIZER_FACTOR);
    	xmlNewChild(node, NULL, BAD_CAST "stabilizerFactor", BAD_CAST tmp_string);
    	sprintf(tmp_string, "%hu", MIN_STIRLINGREQUIRES_NUMBER);
   	    xmlNewChild(node, NULL, BAD_CAST "minStirlingRequiresNumber", BAD_CAST tmp_string);
   	    sprintf(tmp_string, "%hu", DEFAULT_ALGEBRA);
   	    xmlNewChild(node, NULL, BAD_CAST "algebra", BAD_CAST tmp_string);
   	    sprintf(tmp_string, "%.*f", DEFAULT_PRECISION, OUTLIER_CONSTANT);
   	    xmlNewChild(node, NULL, BAD_CAST "outlierConstant", BAD_CAST tmp_string);

		xmlAddChild(root_node, node);
		xmlFreeNode(node);

		node = xmlNewNode(NULL, BAD_CAST "matricesOptions");

		sprintf(tmp_string, "%hu", MAX_RIGHE);
		xmlNewChild(node, NULL, BAD_CAST "maxRaws", BAD_CAST tmp_string);
		sprintf(tmp_string,"%hu", MAX_COLONNE);
		xmlNewChild(node, NULL, BAD_CAST "maxColumns", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", MAX_DSVD_ITERATIONS);
		xmlNewChild(node, NULL, BAD_CAST "maxDSVDIterations", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", MAX_SIMPLEXMETHOD_ITERATIONS);
		xmlNewChild(node, NULL, BAD_CAST "maxSimplexIterations", BAD_CAST tmp_string);

		xmlAddChild(root_node, node);
		xmlFreeNode(node);

		node = xmlNewNode(NULL, BAD_CAST "memoizerOptions");

		#pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
		for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
		{
			char str[2*INFO_STRING] = NULL_CHAR;
			char strboolized[INFO_STRING] = NULL_CHAR;
			strboolize(suite_c.memoizers_names[i], strboolized);
			strboolized[0] = toupper(strboolized[0]);
			sprintf(str, "/settings/memoizerOptions/max%sMemoizableIndex", strboolized);
			xmlNewChild(node, NULL, BAD_CAST str, BAD_CAST suite_c.bools_identifiers[suite_c.max_memoizable_indices[i]]);
		}

		xmlAddChild(root_node, node);
		xmlFreeNode(node);

		node = xmlNewNode(NULL, BAD_CAST "baseConversions");

		sprintf(tmp_string, "%hu", BASECALC_CAMBIAMENTODIBASE_MINBASE);
		xmlNewChild(node, NULL, BAD_CAST "minBase", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", BASECALC_CAMBIAMENTODIBASE_MAXBASE);
		xmlNewChild(node, NULL, BAD_CAST "maxBase", BAD_CAST tmp_string);
		sprintf(tmp_string, "%d", MAX_MINBASE_CONVERTIBLE_NUM);
		xmlNewChild(node, NULL, BAD_CAST "maxChangebaseBinaryConvnum", BAD_CAST tmp_string);

		xmlAddChild(root_node, node);
		xmlFreeNode(node);

		node = xmlNewNode(NULL, BAD_CAST "newtonDifferenceTables");

		sprintf(tmp_string, "%hu", MIN_NEWTON_DIFFTABLES_DIM);
		xmlNewChild(node, NULL, BAD_CAST "minDimension", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", MAX_NEWTON_DIFFTABLES_DIM);
		xmlNewChild(node, NULL, BAD_CAST "maxDimension", BAD_CAST tmp_string);

		xmlAddChild(root_node, node);
		xmlFreeNode(node);

		node = xmlNewNode(NULL, BAD_CAST "romanNumbers");

		sprintf(tmp_string, "%hu", MIN_PROCESSABLE_ROMAN_NUMBER);
		xmlNewChild(node, NULL, BAD_CAST "minProcessableNumber", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", MAX_PROCESSABLE_ROMAN_NUMBER);
		xmlNewChild(node, NULL, BAD_CAST "maxProcessableNumber", BAD_CAST tmp_string);

		xmlAddChild(root_node, node);
		xmlFreeNode(node);

		node = xmlNewNode(NULL, BAD_CAST "pascalsTriangle");

		sprintf(tmp_string, "%hu", MIN_PASCALTRIANGLE_RAWS);
		xmlNewChild(node, NULL, BAD_CAST "minRaws", BAD_CAST tmp_string);
		sprintf(tmp_string, "%hu", MAX_PASCALTRIANGLE_RAWS);
		xmlNewChild(node, NULL, BAD_CAST "maxRaws", BAD_CAST tmp_string);

		xmlAddChild(root_node, node);
		xmlFreeNode(node);

		node = xmlNewNode(NULL, BAD_CAST "booleanKeys");

		#pragma omp parallel for
        for(i=0; i<MAX_BOOL_SETTINGS; ++i)
        {
        	char name[4*MIN_STRING] = NULL_CHAR;
			char strboolized[2*MIN_STRING] = NULL_CHAR;
			strboolize(suite_c.bools_names[i], strboolized);
        	sprintf(name, "/settings/booleanKeys/%s", strboolized);
        	xmlNewChild(node, NULL, BAD_CAST name, BAD_CAST suite_c.bools_identifiers[suite_c.bools[i].default_val]);
        }

		xmlAddChild(root_node, node);
		xmlFreeNode(node);

	    xmlSaveFormatFileEnc(SETTINGS_FILENAME, doc, XML_ENCODING, 1);
	    /// xmlFreeNode(root_node);
	    xmlFreeDoc(doc);
	    xmlCleanupParser();

	    listInsertProc(SETTINGS_FILENAME, LAYOUTS);
    }

    if(mode != MAX_LISTS)
        listInsertProc(argv[1], mode);

    if(extension && !strcmp(extension, DEFAULT_SCRIPTFILES_EXTENSION))
        _execScriptFiles(argv[1]);
    else if(argc > 2)
        _handleCmdLine(argc-2, &argv[2]);

    progInfo(WITHOUT_DESCRIPTION);

    updInfo();
    operationsGroupMenu(MAX_PROGRAMMI, main_menu, NULL_CHAR, BY_CHARS);

    prepareToExit();
    system("PAUSE"); // if working on Windows Environment...

    return NOERROR_EXIT;
}

