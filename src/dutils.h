// dutils.h
// 01/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be included
exclusively by main.c, geometry.c, programs.c, algebra.c and settings.c project files of my suite program!
I do not assume any responsibilities about the use with any other code-scripts.
*/

// HEADER Guard
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <errno.h>
#include <conio.h>
#include <setjmp.h>

#include <windows.h>
#include <errors.h>
#include <signal.h>

#include <omp.h>

#define WINOS ((defined WIN32) || (defined WINDOWS))

#define false 0
#define true 1


// PROJECT HEADERS
// #include "ext_math.h"
#include "defaults.h"
#include "ExprEval\exprincl.h"

#define XMLCALL

#ifdef XMLCALL
	#include "libxml/encoding.h"
	#include "libxml/xmlwriter.h"
	#include "libxml/xpath.h"
	
	#define XML_ENCODING "UTF-8"
	
	#define SETTINGS_FILENAME "./settings.xml"
	
	#define STRING_FALSE "false"
	#define STRING_TRUE "true"
	
	#define XML_FILENAMES_LENGTH MINMIN_STRING
	#define MAX_XML_FIELDSTRINGS SIGN_STRING
#endif
// #include "ExprEval/expreval.h"

// #include "ExprEval/expreval.h"


#ifdef __cplusplus
extern "C" {
#endif


#define __MSSHELL_WRAPPER_
#define __MSNATIVE_
#define __MSUTIL_


#define __WINCALL32
#define __WINCALL __WINCALL32

#define __export // introduced with the passing to Orwell Dev-C++ IDE
#define __apnt
#define __lmp_prog
#define __system
#define _MS__private

#define PROG__NAME "mathSuite"
#define PROG__VERSION "5.75"
#define PROG__AUTHOR "Marco Chiarelli"
#define PROG__LASTUPDATEDATE "01/09/2014"


// INITIALIZING EXPREVAL DEFAULT CONSTANTS
#define INIT_OBJLIST NULL
#define INIT_FUNCLIST NULL
#define INIT_CONSTLIST NULL
#define INIT_VARLIST NULL
#define INIT_ANS NULL

#define VSPACE

#define NOERROR_EXIT 0
#define HEAPALLOC_ERROR -1
#define MSINF_FAULT_ERROR -2


/// MINMAX powerful macros
/// designed to enhance
/// statistic media functions,
/// especially the mediana one.

#define MIN_MODE false
#define MAX_MODE true

#define MIN(a,b) MINMAX(a,b,MIN_MODE,NULL)
#define MAX(a,b) MINMAX(a,b,MAX_MODE,NULL)


/// #define INI_NAME "./settings.xml"

/// #define WINOS true /// ((defined(WIN32)) || (defined(WINDOWS)))

#define XML_NAME "./settings.xml"

#if WINOS
    #define MAX_PATH_LENGTH MAX_PATH
    #define pulisciSchermo system("cls")
#else
    #include <unistd.h>
    #define MAX_PATH_LENGTH PATH_MAX
    #define DEFAULT_LINUX_SPOOLFOLDER "/var/spool"
    #define pulisciSchermo system("clear")
#endif

#define getItemsListNo(mode) access(lists)[mode].itemsno
#define initList(mode) access(lists)[mode].items_list[PREV_LISTNODE] = access(lists)[mode].items_list[NEXT_LISTNODE] = (nodelist *)NULL;

#define DESCRIPTIONS_FOLDER DEFAULT_DESCRIPTIONS_FOLDER


#define INVALID_ITEM NULL
#define STARTING_ITEMSNO 0
#define ITEMS_LASTONESTANDING STARTING_ITEMSNO+1
#define NULL_ITEM(mode) getItemsListNo(mode)

#define NULL_ENV NULL_ITEM(ENVS)
#define NULL_MAT NULL_ITEM(MATRICES)
#define NULL_LOG NULL_ITEM(LOGS)

#define STARTING_ENVSNO STARTING_ITEMSNO
#define STARTING_MATNO  STARTING_ITEMSNO
#define STARTING_LOGSNO STARTING_ITEMSNO
#define STARTING_LAYOUTNO STARTING_ITEMSNO



#define MAIN_ITEM DEFAULT_MAIN_ITEM

#define MAIN_ENV MAIN_ITEM
#define MAIN_MAT MAIN_ITEM
#define MAIN_LOG MAIN_ITEM
#define MAIN_LAYOUT MAIN_ITEM



#define INIT_MATRIXSUMFUNC (void (*)(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS]))NULL
#define INIT_MATRIXPRODFUNC (void (*)(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_ABSTRACT_DIMENSIONS]))NULL
#define INIT_MATRIXKPRODFUNC (void (*)(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS]))NULL



// #define MAIN_ENV DEFAULT_MAIN_ENV
// #define MAIN_MAT DEFAULT_MAIN_MAT
// #define MAIN_LOG DEFAULT_MAIN_LOG


#define INIT_MATLIST NULL
#define INIT_LOGSLIST NULL
#define INIT_CURRENTMATRIX NULL
#define INIT_LASTMATRIXPRINTED NULL
#define INIT_LAYOUTLIST NULL
#define INIT_DIM 0

//
#define INIT_CM_RAWS INIT_DIM
#define INIT_CM_COLUMNS INIT_DIM
#define INIT_LMP_RAWS INIT_DIM
#define INIT_LMP_COLUMNS INIT_DIM
//


#define STARTING_DIM INIT_DIM+1
#define STARTING_CM_RAWS STARTING_DIM
#define STARTING_CM_COLUMNS STARTING_DIM
#define STARTING_LMP_RAWS STARTING_DIM
#define STARTING_LMP_COLUMNS STARTING_DIM
// #define MAT_LASTONESTANDING STARTING_MATNO+1
// #define LOGS_LASTONESTANDING STARTING_LOGNO+1

#define INIT_EXPRTYPE (exprType *)NULL
#define INIT_MATRIXOBJ (matrixObj *)NULL
#define INIT_LOGOBJ (logObj *)NULL
#define INIT_LAYOUTOBJ (layoutObj *)NULL

enum
{
    ERROR_FORBIDDENOPERATION = 1,
    ERROR_NOSUCHFILEORDIR,
    ERROR_NOSUCHPROCESS,
    ERROR_INTERRUPTEDFUNCCALL,
    ERROR_INPUTOUTPUT,
    ERROR_NOSUCHDEVICEADDRESS,
    ERROR_TOOLONGARGLIST,
    ERROR_EXECFORMAT,
    ERROR_BADFILEDESCRIPTOR,
    ERROR_NOCHILDPROCESS,
    ERROR_RESOURCETEMPORARYUNAVAILABLE,
    ERROR_NOTENOUGHSPACE,
    ERROR_PERMISSIONDENIED,
    ERROR_BADADDRESS,
    UNKNOWN_ERROR,
    ERROR_RESOURCEDEVICE,
    ERROR_FILEEXISTS,
    ERROR_IMPROPERLINK,
    ERROR_NOSUCHDEVICE,
    ERROR_NOTADIRECTORY,
    ERROR_ISASHIFT2,
    ERROR_INVALIDARGUMENT,
    ERROR_TOOMANYOPENFILESINSYSTEM,
    ERROR_TOOMANYOPENFILES,
    ERROR_INAPPROPRIATEIOCONTROLOPERATIONS,
    ERROR_FILETOOLARGE = ERROR_INAPPROPRIATEIOCONTROLOPERATIONS+2,
    ERROR_NOSPACELEFTONDEVICE,
    ERROR_INVALIDSEEK,
    ERROR_READONLYFILESYSTEM,
    ERROR_TOOMANYLINKS,
    ERROR_BROKENPIPE,
    ERROR_DOMAIN,
    ERROR_TOOLARGERESULT,
    ERROR_RESORUCEDEADLOCKAVOIDED = ERROR_TOOLARGERESULT+2,
    ERROR_TOOLONGFILENAME,
    ERROR_NOLOCKSAVAILABLE,
    ERROR_NOTIMPLEMENTEDFUNCTION,
    ERROR_DIRECTORYNOTEMPTY,
    ERROR_ILLEGALBYTESEQUENCE
};


#define MAX_BUFSIZ 2048


#define BUFLEN_FACTOR 4 // 5 // in order to increase Maximum Buffer Capacity (Length)

#define MIN_BUFLEN 1
#define MAX_BUFLEN (MAX_BUFSIZ * BUFLEN_FACTOR)

#define DEFAULT_BUFSIZE MAX_BUFSIZ

#define MINMINMIN_STRING 10
#define MINMIN_STRING 20
#define MIN_STRING 50
#define MAX_STRING 256
#define MAXX_STRING (2*MAX_STRING)
#define MAXMAX_STRING (4*MAX_STRING)
#define MMAX_STRING (5*MAX_STRING)
#define KSTRING 1000

#define SIGN_STRING MINMINMIN_STRING
//
#define INFO_STRING MIN_STRING

#define MAX_DESCRIPTION_STRING 3000

#define PATHCONTAINS_STRING (MAX_PATH_LENGTH+(2*MIN_STRING))
#define ERR_STRING PATHCONTAINS_STRING

#define MAX_FILE_LINES MAXX_STRING /// 255 /// INFO_STRING /// 100

#define MAX_IDENTIFIER_LENGTH MIN_STRING

#define MAX_EXTENSION_LENGTH SIGN_STRING

#define NULL_CHAR ""
#define BLANK_STRING " "
#define BLANK_CHAR ' '

#define UNKNOWN_STRING "(unknown)"

#define UNDERSCORE_SEPERATOR_STRING "_"

#define FILENAMES_SEPERATOR UNDERSCORE_SEPERATOR_STRING


// #define RAD_RATIO (M_PI/180)
// #define RAD_RATIO_FACTOR 0.333333333

#define INT64_TMAX ULLONG_MAX
#define NULL_VAL MAX_VAL // INT64_TMAX // starting_random_seed
// #define NULL_VAL2 MAX_VAL-(RAD_RATIO)
// #define NULL_VAL3 MAX_VAL-(RAD_RATIO_FACTOR*RAD_RATIO)

// #define isNullVal(x) (x == NULL_VAL || x == NULL_VAL2 || x == NULL_VAL3)
// #define isNullVal3(x) (x == NULL_VAL || x == NULL_VAL3)

#define isNullVal(x) (x == NULL_VAL)
#define exitHandleCheck (access(exitHandle) == EXITHANDLE_SETCMD || access(exitHandle) == EXITHANDLE_BACKCMD || !access(exitHandle))

#define MATRIXGET_COMMAND DEFAULT_MATRIXGET_COMMAND
#define MATRIXSET_COMMAND DEFAULT_MATRIXSET_COMMAND
#define MATRIXBACK_COMMAND DEFAULT_MATRIXBACK_COMMAND

#define SYSTEM_LOG DEFAULT_SYSTEM_LOG

// #define SECURITY_CHECK ENABLED

#define TERMINATING_CHAR ';'
#define TERMINATING_STRING ";"

#define SCANFEXIT_CHAR '.'
#define SCANFEXIT_STRING "."

#define EXTENSION_DOT "."

#define MATRIXES_SEPERATOR_STRING BLANK_STRING
#define MATRIXES_SEPERATOR_CHAR BLANK_CHAR

#define MIN_EXTENSIVE_MULTITHREADING_CORESNO DEFAULT_MIN_EXTENSIVE_MULTITHREADING_CORESNO

// Simplex METHOD MACROS

enum
{
    SIMPLEXMETHOD_FOUNDBFS = 0,
    SIMPLEXMETHOD_ALLOC_ERROR,
    SIMPLEXMETHOD_INFBFS_ERROR,
    SIMPLEXMETHOD_FARBFS_ERROR
};

#define MIN_PROBLEM false
#define MAX_PROBLEM true

#define MAX_SIMPLEXMETHOD_ITERATIONS DEFAULT_MAX_SIMPLEXMETHOD_ITERATIONS
//

#define MAX_DSVD_ITERATIONS DEFAULT_MAX_DSVD_ITERATIONS


#define OUTLIER_CONSTANT DEFAULT_OUTLIER_CONSTANT

#define MIN_OUTLIER_CONSTANT DEFAULT_MIN_OUTLIER_CONSTANT
#define MAX_OUTLIER_CONSTANT DEFAULT_MAX_OUTLIER_CONSTANT

#define MIDDLE_ELEMENT false
#define INITIAL_ELEMENT true

#define BY_START false
#define BY_GLOBALVAL true

#if WINOS
	#define SetDefaultColor() SetColor(COLOR_DEFAULT)
#endif

#define ROMAN_NUMBER_STRING 2*MIN_STRING

#define ROMAN_NUMBER_MAPSTRING 5

#define MIN_PROCESSABLE_ROMAN_NUMBER 1
#define MAX_PROCESSABLE_ROMAN_NUMBER 3000

#define MIN_PASCALTRIANGLE_RAWS 1
#define MAX_PASCALTRIANGLE_RAWS 50

#define EXIT_MESSAGE DEFAULT_EXIT_MESSAGE

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


// DEFAULT-CUSTOM TYPES TYPEDEFS


typedef long long int64_t;
typedef unsigned long long uint64_t;

//

typedef unsigned char sel_typ; // DEFINING SELECTION TYPE
typedef unsigned short fsel_typ; // DEFINING FORMATTABLE SELECTION TYPE
typedef sel_typ bool; // DEFINING ASSERT TYPE (like bool vars)
// but unlike typedefining a { false, true } enumeration, this system
// combinated to false,true macro definition occupies less memory.

// typedef fsel_typ dim_typ;

#define COMPARE_STRING_FACTOR MAX_STRING

//
#if (MAX_RIGHE > 0 && MAX_COLONNE > 0)
#if (MAX_RIGHE < COMPARE_STRING_FACTOR && MAX_COLONNE < COMPARE_STRING_FACTOR)
typedef fsel_typ dim_typ;
#else
typedef unsigned dim_typ;
#endif
#else
    #error Invalid definition of MAX_RIGHE or MAX_COLONNE.
#endif
//

typedef tipo_predefinito ityp;

// #include "ExprEval/expreval.h"

#define EXPREVAL_TMPL ".\ExprEval\exprtmpl.html" /// "http://expreval.sourceforge.net/exprtmpl.html"

// INCLUDING EXPREVAL STUFFS BY DIRECTLY
/* Define type of data to use */
// typedef double EXPRTYPE;
typedef ityp EXPRTYPE;

/* Defines for various things */

/* Max id size */
#define EXPR_MAXIDENTSIZE 255

/* Error values */
enum
    {
    EXPR_ERROR_UNKNOWN = -1, /* Unknown error */
    EXPR_ERROR_NOERROR = 0, /* No Error */
    EXPR_ERROR_MEMORY, /* Memory allocation failed */
    EXPR_ERROR_NULLPOINTER, /* Null pointer passed to function */
    EXPR_ERROR_NOTFOUND, /* Item not found in a list */
    EXPR_ERROR_UNMATCHEDCOMMENT, /* Unmatched comment tags */
    EXPR_ERROR_INVALIDCHAR, /* Invalid characters in expression */
    EXPR_ERROR_ALREADYEXISTS, /* An item already called create */
    EXPR_ERROR_ALREADYPARSEDBAD, /* Expression parsed already, but unsuccessfully. call free or clear */
    EXPR_ERROR_ALREADYPARSEDGOOD, /* Expression parsed already, successfully, call free or clear */
    EXPR_ERROR_EMPTYEXPR, /* Empty expression string passed to parse */
    EXPR_ERROR_UNMATCHEDPAREN, /* Unmatched parenthesis */
    EXPR_ERROR_SYNTAX, /* Syntax error in expression */
    EXPR_ERROR_MISSINGSEMICOLON, /* Missing semicolon at end of expression */
    EXPR_ERROR_BADIDENTIFIER, /* Identifier was to big or not formed right */
    EXPR_ERROR_NOSUCHFUNCTION, /* Function does not exist in function list */
    EXPR_ERROR_BADNUMBERARGUMENTS, /* Bad number of arguments in a function call */
    EXPR_ERROR_BADEXPR, /* This is a bad expression to evaluate. It has not been parsed or has unsuccessfully */
    EXPR_ERROR_UNABLETOASSIGN, /* Unable to do an assignment, maybe no variable list */
    EXPR_ERROR_DIVBYZERO, /* Attempted a division by zero */
    EXPR_ERROR_NOVARLIST, /* No variable list found but one is needed */
    EXPR_ERROR_BREAK, /* Expression was broken by break function */
    EXPR_ERROR_CONSTANTASSIGN, /* Assignment to a constant */
    EXPR_ERROR_REFCONSTANT, /* Constant used as a reference parameter */
    EXPR_ERROR_OUTOFRANGE, /* A bad value was passed to a function */

    EXPR_ERROR_USER /* Custom errors should be larger than this */
    };

/* Macros */

/* Forward declarations */
typedef struct _exprNode exprNode;
typedef struct _exprFuncList exprFuncList;
typedef struct _exprValList exprValList;
typedef struct _exprObj exprObj;

/* Function types */
typedef int (*exprFuncType)(exprObj *obj, exprNode *nodes, int nodecount, EXPRTYPE **refs, int refcount, EXPRTYPE *val);
typedef int (*exprBreakFuncType)(exprObj *obj);



/* Functions */

/* Version information function */
void exprGetVersion(int *major, int *minor);

/* Functions for function lists */
int exprFuncListCreate(exprFuncList **flist);
int exprFuncListAdd(exprFuncList *flist, char *name, exprFuncType ptr, int min, int max, int refmin, int refmax);
int exprFuncListFree(exprFuncList *flist);
int exprFuncListClear(exprFuncList *flist);
int exprFuncListInit(exprFuncList *flist);

/* Functions for value lists */
int exprValListCreate(exprValList **vlist);
int exprValListAdd(exprValList *vlist, char *name, EXPRTYPE val);
int exprValListSet(exprValList *vlist, char *name, EXPRTYPE val);
int exprValListGet(exprValList *vlist, char *name, EXPRTYPE *val);
int exprValListAddAddress(exprValList *vlist, char *name, EXPRTYPE *addr);
int exprValListGetAddress(exprValList *vlist, char *name, EXPRTYPE **addr);
void *exprValListGetNext(exprValList *vlist, char **name, EXPRTYPE *value, EXPRTYPE** addr, void *cookie);
int exprValListFree(exprValList *vlist);
int exprValListClear(exprValList *vlist);
int exprValListInit(exprValList *vlist);

/* Functions for expression objects */
int exprCreate(exprObj **obj, exprFuncList *flist, exprValList *vlist, exprValList *clist,
    exprBreakFuncType breaker, void *userdata);
int exprFree(exprObj *obj);
int exprClear(exprObj *obj);
int exprParse(exprObj *obj, char *expr);
int exprEval(exprObj *obj, EXPRTYPE *val);
int exprEvalNode(exprObj *obj, exprNode *nodes, int curnode, EXPRTYPE *val);
exprFuncList *exprGetFuncList(exprObj *obj);
exprValList *exprGetVarList(exprObj *obj);
exprValList *exprGetConstList(exprObj *obj);
exprBreakFuncType exprGetBreakFunc(exprObj *obj);
int exprGetBreakResult(exprObj *obj);
void* exprGetUserData(exprObj *obj);
void exprSetUserData(exprObj *obj, void *userdata);
void exprSetBreakCount(exprObj *obj, int count);
void exprGetErrorPosition(exprObj *obj, int *start, int *end);

/* Other useful routines */
int exprValidIdent(char *name);
//

// Removed because It caused overhead in Stirling's function
// (There is another sqrt function which is assumed to be multiplied in the Formula
// #define STIRLING_CONSTANT sqrt(2*M_PI)

#define MAX_MEMOIZABLE_INDEX 100

#define MAX_FIBONACCI_MEMOIZABLE_INDEX DEFAULT_MAX_FIBONACCI_MEMOIZABLE_INDEX
#define MAX_FATTORIALE_MEMOIZABLE_INDEX DEFAULT_MAX_FATTORIALE_MEMOIZABLE_INDEX
#define MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX DEFAULT_MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX
#define MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX DEFAULT_MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX


#if (MAX_FIBONACCI_MEMOIZABLE_INDEX > MAX_MEMOIZABLE_INDEX || MAX_FATTORIALE_MEMOIZABLE_INDEX > MAX_MEMOIZABLE_INDEX || MAX_FPNUM_MEMOIZABLE_INDEX > MAX_MEMOIZABLE_INDEX)
    #error Invalid definition of one of MAX FUNCTIONS MEMOIZABLE INDICES.
#endif // MAX_FIBONACCI_MEMOIZABLE_INDEX

#define MIN_STIRLINGREQUIRES_NUMBER DEFAULT_MIN_STIRLINGREQUIRES_NUMBER

#define INPUT false
#define OUTPUT true

#define OVERFLOW_ERROR "\nINFO: An OVERFLOW ERROR might have occurred during (%hu,%hu) Elaborating Matrix Element...\n"

#define CALLOC_MODE false
#define MALLOC_MODE true

#define FIRST_MATRIX 0
#define SECOND_MATRIX 1

#define NORMAL_MODE false
#define NOP_MODE true

//
#define insertMatrix(x,y,z,k) enterMatrix(&x,&y,&z,k,true)
#define matrixFree(x,y) _matrixFree(x,y,NORMAL_MODE)

/// #define matrixAlloc(x,y) _matrixAlloc(x,(dim_typ2)y,MALLOC_MODE)
//

// FORMAT MACROS
//


// COLOR Enumerations
enum
{
    COLOR_BLACK = 0,
    COLOR_BLUE,
    COLOR_GREEN,
    COLOR_AZURE,
    COLOR_DARKRED,
    COLOR_PURPLE,
    COLOR_YELLOW,
    COLOR_WHITE,
    COLOR_GREY,
    COLOR_SMARTBLUE,
    COLOR_LUMGREEN,
    COLOR_CRYSTAL,
    COLOR_SMARTRED,
    COLOR_SMARTPURPLE,
    COLOR_SMARTYELLOW,
    COLOR_WHITEPLUS
};

#define COLOR_SMARTWHITE 39
#define COLOR_LIGHTBLUE 43

#define _COLOR_ERROR 60 // COLOR_SMARTRED
#define _COLOR_CREDITS 57 // COLOR_LIGHTBLUE // COLOR_CRYSTAL // COLOR_SMARTBLUE
#define _COLOR_USER 58 // COLOR_LUMGREEN
#define _COLOR_SYSTEM 60 // COLOR_SMARTYELLOW
#define _COLOR_AUTHOR 54 // COLOR_SMARTPURPLE

enum
{
    MEMBER_DEFAULTCOLOR = 0,
    MEMBER_COLORERROR,
    MEMBER_COLORCREDITS,
    MEMBER_COLORUSER,
    MEMBER_COLORSYSTEM,
    MEMBER_COLORAUTHOR
};

#define INIT_COLOR COLOR_WHITE
#define LAST_MEMBER_COLOR MEMBER_COLORAUTHOR
#define MAX_COLOR_TYPES LAST_MEMBER_COLOR+1
#define COLORS_PATH DEFAULT_COLORS_PATH
#define MAX_COLORS 64

#define COLOR_DEFAULT access(colors)[MEMBER_DEFAULTCOLOR]
#define COLOR_ERROR access(colors)[MEMBER_COLORERROR]
#define COLOR_CREDITS access(colors)[MEMBER_COLORCREDITS]
#define COLOR_USER access(colors)[MEMBER_COLORUSER]
#define COLOR_SYSTEM access(colors)[MEMBER_COLORSYSTEM]
#define COLOR_AUTHOR access(colors)[MEMBER_COLORAUTHOR]

#define _sprint(x,y)                                                             \
{                                                                                \
    va_list ap;                                                                  \
    va_start(ap,x);                                                              \
    char str[MAX_BUFSIZ+INFO_STRING];                                            \
    const bool cond = access(sysLog) != INIT_LOGOBJ;                             \
	SetColor(y);                                                             	 \
    vsprintf(str, x, ap);                                                        \
    prependTimeToString(str, cond);                                              \
    printf(str);                                                  				 \
    SetDefaultColor();                                                           \
    if(cond)                                                                     \
        logCheck(access(sysLog), str, access(sysLogPath));                       \
    va_end(ap);                                                                  \
}																				 \

#define PRINTSPACE() printf(" ");
#define PRINTT() printf("\t");
#define PRINTL() printf("________________________________________________________________________________\n\n") // enhanced formatting

#define PRINTN() printf("\n");
#define PRINT2N() printf("\n\n");
#define PRINT5N() printf("\n\n\n\n\n");
#define PRINT10N() printf("\n\n\n\n\n\n\n\n\n\n");

#define SHOWPAUSEMESSAGE() PRINTL(), printf2(COLOR_SYSTEM, "Press CTRL + C to stop Program Action.\n\n");
#define PRINTHOWTOBACKMESSAGE() printf2(COLOR_SYSTEM, "Press ENTER to exit SubProgram.\n\n");


#define __pmode__ access(mode) // MACRO PROGRAM MODE

#define DISABLED false
#define ENABLED true

#define BY_NUMBERS false
#define BY_CHARS true

#define EXIT_CMD "exit"
#define _MSS_CMD DEFAULT_SCRIPTFILES_EXTENSION"_init"
#define MSS_CMD DEFAULT_SCRIPTFILES_EXTENSION"_exit"

#define INITIALIZING_RANDOM_SEED 0 // Invalid value for Random Seed.
// Used just to initializate the program random seed value.

#define DEFAULT_RANDOM_SEED (unsigned)time(NULL)
volatile tipo_random_seed starting_random_seed;

#define MIN_RANDOMSEED 1
#define MAX_RANDOMSEED RAND_MAX

#define MIN_STABFACT 1
#define MAX_STABFACT DEFAULT_MAX_STABFACT // 20

enum
{
    ALGEBRA_REALNUMBERS,
    ALGEBRA_COMPLEXNUMBERS,
    ALGEBRA_QUATERNIONS,
    ALGEBRA_OCTONIONS,
    ALGEBRA_SEDENIONS
};

#define LAST_ALGEBRA ALGEBRA_SEDENIONS

#define MIN_ALGEBRA ALGEBRA_REALNUMBERS
#define MAX_ALGEBRA LAST_ALGEBRA

#define _MAX_ALGEBRA MAX_ALGEBRA+1

#define REAL_UNIT_NAME "real"

#define MIN_MEMOIZABLE_INDEX 5


#define STABILIZER_FACTOR access(curLayout)->stabilizer_factor


// DEFINIZIONE MACRO SOTTO-PROGRAMMI
//

// Valore di inizializzazione del metadato modalità
// della variabile strutturata suite, definita più sotto
#define PROGRAM_BUSY -1
// o più semplicemente indica che l'utente è in fase di scelta del subprogram.


// Enumerazione ID Sotto-Programmi
// per ogni Famiglia
//

// MAIN PROGRAMS ENUMERATIONS
enum
{
    MAIN_BASECALCULATOR = 0,
    MAIN_ADVANCEDCALCULATOR,
    MAIN_ALGEBRAOPERATIONS,
    #ifdef ALLOW_MSSMANAGER
        MAIN_MSSMANAGER,
    #endif // ALLOW_MSSMANAGER
    MAIN_CHANGESETTINGS
};

// MACRO PROGRAMMI
// Immettere la macro dell'ultimo programma definito
#define LAST_PROGRAM MAIN_CHANGESETTINGS // Cambiare ogni qualvola si aggiunge
// un sottoprogramma alla relativa enumerazione

//
#define MAX_PROGRAMMI LAST_PROGRAM+1

// ADVANCED CALCULATOR MACRO AND ENUMERATIONS
//

enum
{
    MSSMANAGER_COMMANDLINE = 0,
    MSSMANAGER_EXECSCRIPTFILES,
    MSSMANAGER_SHOWUSAGES,
};

#define MAX_ARGS 5
#define LAST_MSSMANAGER_PROG MSSMANAGER_SHOWUSAGES
#define MAX_MSSMANAGER_PROGS LAST_MSSMANAGER_PROG+1

enum
{
    ADVCALC_SECONDGRADEEQUATIONSOLVER = 0,
    ADVCALC_COMPLEXNUMBERSSUM,
    ADVCALC_COMPLEXNUMBERSPROD,
    ADVCALC_GETFORMATTEDDATE,
    ADVCALC_SIMPLEXMETHOD,
    ADVCALC_NEWTONDIFFTABLES,
    ADVCALC_LAGRANGEINTERPOLATION,
    ADVCALC_GREATESTEIGENVALUE,
    ADVCALC_FUNCTIONINTEGRATION,
    ADVCALC_STRAIGHTLINEFITTING,
    ADVCALC_PARABOLICCURVEFITTING,
    ADVCALC_LINEARSYSTEMSSOLVER
};

#define LAST_ADVCALC_PROG ADVCALC_LINEARSYSTEMSSOLVER
#define MAX_ADVCALC_PROGS LAST_ADVCALC_PROG+1

// LINEAR ALGEBRA MACRO AND ENUMERATIONS
//
enum
{
    ALGEBRA_MATRICESMANAGER = 0,
    ALGEBRA_MATRIXSORT,
    ALGEBRA_NORMCALCULATOR,
    ALGEBRA_DETERMINANTCALCULATOR,
    ALGEBRA_TRACECALCULATOR,
    ALGEBRA_RANKCALCULATOR,
    ALGEBRA_MATRIXSVD,
    ALGEBRA_INVERSEMATRIX,
    ALGEBRA_MATRIXTRANSPOSE,
    ALGEBRA_MATRICESSUM,
    ALGEBRA_TENSORSSUM,
    ALGEBRA_MATRICESPRODUCT,
    ALGEBRA_MATRICESKPRODUCT,
    ALGEBRA_MATRIXPOWER,
    ALGEBRA_MATRIXKPOWER,
    ALGEBRA_MATRIXPERVECTOR,
    ALGEBRA_VECTORPERVECTOR,
    ALGEBRA_SCALARPRODUCT,
    ALGEBRA_PERSCALARPRODUCT,
    ALGEBRA_SCALARDIVISIONMATRIX,
    ALGEBRA_ILLCONDITIONCHECKING,
    ALGEBRA_MATRIXLUFACTORIZATION
};


#define LAST_ALGEBRA_PROGRAM ALGEBRA_MATRIXLUFACTORIZATION
#define MAX_ALGEBRA_OPERATIONS LAST_ALGEBRA_PROGRAM+1


//DEFINING ENUMERATIONS OF MAIN_CHANGESETTINGS SUBPROGRAM

/*
IMPORTANTE!!! Lasciare sempre per ultimo l'enumerazione
SETTINGS_SETDEFAULTS, violando tale regola si potrebbe interferire
con il sistema pre-compiling di freezing di alcune impostazioni
dal settings.xml
*/
//
enum
{
    #ifdef ALLOW_VARLISTMANAGER
        SETTINGS_ENVSMANAGER,
    #endif
    #ifdef ALLOW_LOGSMANAGER
        SETTINGS_LOGSMANAGER,
    #endif
    #ifdef ALLOW_SYSLOGMANAGER
        SETTINGS_SYSLOGMANAGER,
    #endif
    #ifdef ALLOW_LAYOUTSMANAGER
    SETTINGS_LAYOUTSMANAGER,
    #endif
    #ifdef ALLOW_COLSMANAGER
    SETTINGS_COLORSMANAGER,
    #endif
    #ifdef ALLOW_LFSMANAGER
        SETTINGS_LFSMANAGER,
    #endif // ALLOW_LFSMANAGER
    #ifdef ALLOW_PRECEDIT
        SETTINGS_CHANGEPRECISION,
    #endif
    #ifdef ALLOW_STABFACTEDIT
        SETTINGS_CHANGESTABILIZERFACTOR,
    #endif
    #ifdef ALLOW_MINSRNUMBEREDIT
    	SETTINGS_CHANGEMINSRNUMBER,
    #endif
    #ifdef ALLOW_ALGEBRAEDIT
        SETTINGS_CHANGEALGEBRA,
    #endif
    #ifdef ALLOW_OUTLIERCONSTEDIT
    	SETTINGS_CHANGEOUTLIERCONST,
    #endif
    #ifdef ALLOW_EXITCHAREDIT
        SETTINGS_CHANGEEXITCHAR,
    #endif
    #ifdef ALLOW_RANDOMSEEDEDIT
        SETTINGS_CHANGERANDOMSEED,
    #endif
    #ifdef ALLOW_BOOLVARSEDIT
        SETTINGS_CHANGEBOOLVALUES,
    #endif
    #ifdef ALLOW_OPERIDSEDIT
        SETTINGS_CHANGEOPERATIONSIDENTIFIERS,
    #endif
    #ifdef ALLOW_MMIEDIT
        SETTINGS_CHANGEMAXMEMOIZABLEINDICES,
        SETTINGS_EMPTYMEMOIZERSBUFFERS,
    #endif
    #ifdef ALLOW_BUFFERSEDIT
        SETTINGS_EMPTYBUFFERS,
    #endif
    #ifdef ALLOW_SETDEFAULTS
        SETTINGS_SETDEFAULTS
    #endif
};

#define LAST_SETTINGS_PROGRAM SETTINGS_SETDEFAULTS
#define MAX_SETTINGS LAST_SETTINGS_PROGRAM+1


enum
{
    BOOLS_SHOWDESCRIPTION = 0,
    BOOLS_AUTOSETCURITEM,
    BOOLS_ITEMSSELECTBYPATH,
    BOOLS_ITEMSAUTOSAVING,
    #if WINOS
        BOOLS_SYSLOGSECURITYCHECK,
    #endif
    BOOLS_SYSTEMPARSING,
    BOOLS_BASECALCPARSING,
    BOOLS_ADVCALCPARSING,
    BOOLS_MATRIXPARSING,
    BOOLS_RESETLISTS,
    BOOLS_SAVERESULTS,
    BOOLS_SHOWVARLIST,
    BOOLS_SHOWDIFFTIME,
    BOOLS_SHOWEXECTIME,
    BOOLS_PRINTTIME,
    BOOLS_DOMAINCHECK,
    BOOLS_CHARINSERT,
    BOOLS_INSERTMODE,
    BOOLS_LAZYEXECUTION,
    BOOLS_INVERSEOPERATIONS,
    BOOLS_AUTOTURNBACK,
    BOOLS_DEGREESENTERING,
    BOOLS_PROGREPEATCHECK,
    BOOLS_EXTENSIVEMULTITHREADING
};

#define BITMASK_SHOWDESCRIPTION 1
#define BITMASK_AUTOSETCURITEM 2
#define BITMASK_ITEMSSELECTBYPATH 4
#define BITMASK_ITEMSAUTOSAVING 8
#define BITMASK_SYSLOGSECURITYCHECK 16
#define BITMASK_SYSTEMPARSING 32
#define BITMASK_BASECALCPARSING 64
#define BITMASK_ADVCALCPARSING 128
#define BITMASK_MATRIXPARSING 256
#define BITMASK_RESETLISTS 512
#define BITMASK_SAVERESULTS 1024
#define BITMASK_SHOWVARLIST 2048
#define BITMASK_SHOWDIFFTIME 4096
#define BITMASK_SHOWEXECTIME 8192
#define BITMASK_PRINTTIME 16384
#define BITMASK_DOMAINCHECK 32768
#define BITMASK_CHARINSERT 65536
#define BITMASK_INSERTMODE 131072
#define BITMASK_LAZYEXECUTION 262144
#define BITMASK_INVERSEOPERATIONS 524288
#define BITMASK_AUTOTURNBACK 1048576
#define BITMASK_DEGREESENTERING 2097152
#define BITMASK_PROGREPEATCHECK 4194304
#define BITMASK_EXTENSIVEMULTITHREADING 8388608

#define LAST_BOOL_SETTINGS BOOLS_EXTENSIVEMULTITHREADING
#define MAX_BOOL_SETTINGS LAST_BOOL_SETTINGS+1

#define isSett(x) ((access(curLayout)->bools & suite_c.bools[x].bmask) == suite_c.bools[x].bmask)
#define isnSett(x) (!(access(curLayout)->bools & suite_c.bools[x].bmask))

// ENVS MANAGER ENUMERATIONS

enum
{
    ENVS_SETCURRENT = 0,
    ENVS_OPEN,
    ENVS_CREATE,     //!C
    ENVS_READ,       //!R
    ENVS_PRINT,        //!R
    ENVS_UPDATE,        //!U
    ENVS_UPDATEALL,    //!U
    ENVS_SAVE,       //!U
    ENVS_DELETE,        //!D
    ENVS_DELETEPHYSICAL,       //!D
    ENVS_DELETEALL,    //!D
    ENVS_DELETEALLPHYSICALS,   //!D
    ENVS_RELOAD,        //!D
    ENVS_SYNC,
    ENVS_RENAME
};

#define LAST_ENVSMANAGER_PROG ENVS_RENAME
#define MAX_ENVSMANAGER_PROGS LAST_ENVSMANAGER_PROG+1

enum
{
    MATRICES_SETCURRENT = 0,
    MATRICES_OPEN,
    MATRICES_CREATE,     //!C
    MATRICES_READ,       //!R
    MATRICES_PRINT,        //!R
    MATRICES_UPDATE,        //!U
    MATRICES_UPDATEALL,     //!U
    MATRICES_SAVE,       //!U
    MATRICES_DELETE,        //!D
    MATRICES_DELETEPHYSICAL,       //!D
    MATRICES_DELETEALL,     //!D
    MATRICES_DELETEALLPHYSICALS,    //!D
    MATRICES_EDIT,       //!U
    MATRICES_RENAME,        //!U
};

#define LAST_MATMANAGER_PROG MATRICES_RENAME
#define MAX_MATMANAGER_PROGS LAST_MATMANAGER_PROG+1

enum
{
    LOGS_SETCURRENT = 0,
    LOGS_OPEN,       //!C
    LOGS_CREATE,     //!C
    LOGS_READ,       //!R
    LOGS_PRINT,        //!R
    LOGS_UPDATE,        //!U
    LOGS_UPDATEALL,    //!U
    LOGS_SAVE,       //!U
    LOGS_DELETE,        //!D
    LOGS_DELETEPHYSICAL,       //!D
    LOGS_DELETEALL,    //!D
    LOGS_DELETEALLPHYSICALS,   //!D
    LOGS_RELOAD,        //!D
    LOGS_EDIT,       //!U
    LOGS_SETBUFLEN,  //!U
    LOGS_EMPTYBUFFER,   //!U
    LOGS_RENAME         //!U
};

#define LAST_LOGSMANAGER_PROG LOGS_RENAME
#define MAX_LOGSMANAGER_PROGS LAST_LOGSMANAGER_PROG+1

enum
{
    SYSLOG_OPEN = 0,    //!C
    SYSLOG_CREATE,      //!C
    SYSLOG_READ,        //!R
    SYSLOG_PRINT,         //!R
    SYSLOG_EDIT,        //!U
    SYSLOG_UPDATE,         //!U
    SYSLOG_SAVE,        //!U
    SYSLOG_RELOAD,         //!D
    SYSLOG_SETBUFLEN,   //!U
    SYSLOG_EMPTYBUFFER,    //!U
    SYSLOG_RENAME          //!U
};

#define LAST_SYSLOGMANAGER_PROG SYSLOG_RENAME
#define MAX_SYSLOGMANAGER_PROGS LAST_SYSLOGMANAGER_PROG+1

enum
{
    LAYOUTS_SETCURRENT = 0,
    LAYOUTS_OPEN,
    LAYOUTS_CREATE,     //!C
    LAYOUTS_READ,       //!R
    LAYOUTS_PRINT,        //!R
    LAYOUTS_UPDATE,        //!U
    LAYOUTS_UPDATEALL,    //!U
    LAYOUTS_SAVE,       //!U
    LAYOUTS_DELETE,        //!D
    LAYOUTS_DELETEPHYSICAL,       //!D
    LAYOUTS_DELETEALL,    //!D
    LAYOUTS_DELETEALLPHYSICALS,   //!D
    LAYOUTS_RELOAD,        //!D
    LAYOUTS_RENAME
};

#define LAST_LAYOUTSMANAGER_PROG LAYOUTS_RENAME
#define MAX_LAYOUTSMANAGER_PROGS LAST_LAYOUTSMANAGER_PROG+1

enum
{
    COLORS_CHANGE = 0,
    COLORS_FILESLOADER,
    COLORS_BACKUPFILES
};

#define LAST_COLSMANAGER_PROG COLORS_BACKUPFILES
#define MAX_COLSMANAGER_PROGS LAST_COLSMANAGER_PROG+1

enum
{
    LFS_LOADER = 0,
    LFS_CREATE
};

#define LAST_LFSMANAGER_PROG LFS_CREATE
#define MAX_LFSMANAGER_PROGS LAST_LFSMANAGER_PROG+1

#define MAX_CASEINSENSITIVE_CHARS_ALPHABET 26

// Macro per alcune variabili booleane del programma
#define dcheck isSett(BOOLS_DOMAINCHECK)
#define char_insert isSett(BOOLS_CHARINSERT)
#define lazy_exec isSett(BOOLS_LAZYEXECUTION)
#define PARSING_SYSTEM_ALLOWED isSett(BOOLS_SYSTEMPARSING)


#define PROGRAM_INFORMATIONS MAX_PROGRAMMI
#define EXIT_FROM_PROGRAM MAX_PROGRAMMI+1


// VALORE DI INIZIALIZZAZIONE CARATTERE DI USCITA DAL PROGRAMMA
#define INITIALIZING_EXIT_CHAR '£'
#define INITIALIZING_DEFAULT_COLOR COLOR_WHITE

#define MAIN_COLOR DEFAULT_COLOR

// ENUMERAZIONE OPERAZIONI

#define OPERATION_ID 0
#define ALIAS_ID 1


enum
{
    BASECALC_EXIT = 0,
    BASECALC_ADDIZIONE,
    BASECALC_SOTTRAZIONE,
    BASECALC_MOLTIPLICAZIONE,
    BASECALC_DIVISIONE,
    BASECALC_RESTO,
    BASECALC_ADDIZIONEBINARIA,
    BASECALC_SOTTRAZIONEBINARIA,
    BASECALC_COMPLEMENTO,
    BASECALC_ELEVAMENTOAPOTENZA,
    BASECALC_EXPANDEXPC,
    BASECALC_EXP10ANDEXP10C,
    BASECALC_EXP2ANDEXP2C,
    BASECALC_RADICENESIMA,
    BASECALC_LOGARITMO,
    BASECALC_LOGARITMO2,
    BASECALC_LOGARITMOBN,
    BASECALC_LOGARITMOC,
    BASECALC_LOGARITMO2C,
    BASECALC_LOGARITMO1P,
    BASECALC_LOGARITMO1PC,
    BASECALC_BITCOUNTER,
    BASECALC_UBITCOUNTER,
    BASECALC_VERSION,
    BASECALC_EXITCHAR,
    BASECALC_PREC,
    BASECALC_SFACT,
    BASECALC_MINSRNUMBER,
    BASECALC_ALGEBRA,
    BASECALC_OUTLIERCONST,
    BASECALC_RSEED,
    BASECALC_MMIFIBO,
    BASECALC_MMIFACT,
    BASECALC_MMIEVENSFACT,
    BASECALC_MMIODDSFACT,
    BASECALC_TRASFORMAANGOLI,
    BASECALC_SINANDSINH,
    BASECALC_COSANDCOSH,
    BASECALC_TANANDTANH,
    BASECALC_CSCANDCSCH,
    BASECALC_SECANDSECH,
    BASECALC_COTANDCOTH,
    BASECALC_HSINANDHSINH,
    BASECALC_QSINANDQSINH,
    BASECALC_HCOSANDHCOSH,
    BASECALC_QCOSANDQCOSH,
    BASECALC_HSECANDHSECH,
    BASECALC_QSECANDQSECH,
    BASECALC_HCSCANDHCSCH,
    BASECALC_QCSCANDQCSC,
    BASECALC_HTANANDHTANH,
    BASECALC_QTANANDQTANH,
    BASECALC_HCOTANDHCOTH,
    BASECALC_QCOTANDQCOTH,
    BASECALC_VSINANDVSINH,
    BASECALC_CVSINANDCVSINH,
    BASECALC_VCOSANDVCOSH,
    BASECALC_CVCOSANDCVCOSH,
    BASECALC_HVSINANDHVSINH,
    BASECALC_HCVSINANDHCVSINH,
    BASECALC_QVSINANDQVSINH,
    BASECALC_QCVSINANDQCVSINH,
    BASECALC_HVCOSANDHVCOSH,
    BASECALC_HCVCOSANDHCVCOSH,
    BASECALC_QVCOSANDQVCOSH,
    BASECALC_QCVCOSANDQCVCOSH,
    BASECALC_ESECANDESECH,
    BASECALC_ECSCANDECSCH,
    BASECALC_HESECANDHESECH,
    BASECALC_HECSCANDHECSCH,
    BASECALC_QESECANDQESECH,
    BASECALC_QECSCANDQECSCH,
    BASECALC_SINCANDSINCH,
    BASECALC_HSINCANDHSINCH,
    BASECALC_QSINCANDQSINCH,
    BASECALC_COSCANDCOSCH,
    BASECALC_HCOSCANDHCOSCH,
    BASECALC_QCOSCANDQCOSCH,
    BASECALC_SECCANDSECCH,
    BASECALC_HSECCANDHSECCH,
    BASECALC_QSECCANDQSECCH,
    BASECALC_CSCCANDCSCCH,
    BASECALC_HCSCCANDHCSCCH,
    BASECALC_QCSCCANDQCSCCH,
    BASECALC_TANCANDTANCH,
    BASECALC_HTANCANDHTANCH,
    BASECALC_QTANCANDQTANCH,
    BASECALC_COTCANDCOTCH,
    BASECALC_HCOTCANDHCOTCH,
    BASECALC_QCOTCANDQCOTCH,
    BASECALC_ASINANDASINH,
    BASECALC_ACOSANDACOSH,
    BASECALC_ATANANDATANH,
    BASECALC_ATAN2,
    BASECALC_ACSCANDACSCH,
    BASECALC_ASECANDASECH,
    BASECALC_ACOTANDACOTH,
    BASECALC_MCM,
    BASECALC_MCD,
    BASECALC_APPROSSIMAZIONE,
    BASECALC_SOMMASUCCESSIONEGEOMETRICA,
    BASECALC_SOMMASUCCESSIONEARMONICA,
    BASECALC_SOMMASUCCESSIONEARMONICAGEN,
    BASECALC_SOMMASUCCESSIONEFIBONACCI,
    BASECALC_SOMMASUCCESSIONEFATTORIALE,
    BASECALC_SOMMASUCCESSIONESEMIFATTORIALE,
    BASECALC_SOMMATORIA,
    BASECALC_PRODUTTORIA,
    BASECALC_MEDIA,
    BASECALC_VARIANZA,
    BASECALC_STDDEV,
    BASECALC_OUTLIER,
    BASECALC_MAP,
    BASECALC_MEDIAGEOMETRICA,
    BASECALC_MEDIAARMONICA,
    BASECALC_MEDIAPOTENZA,
    BASECALC_VALORECENTRALE,
    BASECALC_PRIMOQUARTILE,
    BASECALC_MEDIANA,
    BASECALC_TERZOQUARTILE,
    BASECALC_SOMMAPRIMINNUMERI,
    BASECALC_RADICEQUADRATA,
    BASECALC_RADICECUBICA,
    BASECALC_FATTORIALE,
    BASECALC_SEMIFATTORIALE,
    BASECALC_STIRLING,
    BASECALC_FIBONACCI,
    BASECALC_PERMUTATIONS,
    BASECALC_PERMUTATIONSREP,
    BASECALC_KPERMUTATIONS,
    BASECALC_KPERMUTATIONSREP,
    BASECALC_COMBINATIONS,
    BASECALC_COMBINATIONSREP,
    BASECALC_PASCALTRIANGLE,
    BASECALC_NUMERIROMANI,
    BASECALC_PRIMINNUMERIPRIMI,
    BASECALC_NESIMONUMEROPRIMO,
    BASECALC_PRIMORIALE,
    BASECALC_SOMMAPRIMINNUMERIPRIMI,
    BASECALC_FIBONACCIALE,
    BASECALC_CAMBIAMENTODIBASE,
    BASECALC_GENERATOREMATRICIRANDOM,
    BASECALC_INFORMAZIONI
};


#define LAST_OPERATION BASECALC_INFORMAZIONI
// cambiare ogni qualvolta si aggiunge un'operazione
// all'enumerazione operazioni e quindi all'array di operazioni di default default_operazioni
//

#define MAX_OPERATIONS LAST_OPERATION+1


enum
{
    FID_SIN = 0,
    FID_COS,
    FID_SINH,
    FID_COSH,
    FID_SEC,
    FID_SECH,
    FID_CSC,
    FID_CSCH,
    FID_ASIN,
    FID_ASINH,
    FID_ACOS,
    FID_ACOSH,
    FID_ASEC,
    FID_ASECH,
    FID_ACSC,
    FID_ACSCH,
    FID_TAN,
    FID_TANH,
    FID_ATAN,
    FID_ATANH,
    FID_COT,
    FID_COTH,
    FID_ACOT,
    FID_ACOTH,
    FID_HSIN,
    FID_HSINH,
    FID_QSIN,
    FID_QSINH,
    FID_HCOS,
    FID_HCOSH,
    FID_QCOS,
    FID_QCOSH,
    FID_HSEC,
    FID_HSECH,
    FID_QSEC,
    FID_QSECH,
    FID_HCSC,
    FID_HCSCH,
    FID_QCSC,
    FID_QCSCH,
    FID_HTAN,
    FID_HTANH,
    FID_QTAN,
    FID_QTANH,
    FID_HCOT,
    FID_HCOTH,
    FID_QCOT,
    FID_QCOTH,
    FID_VSIN,
    FID_VSINH,
    FID_CVSIN,
    FID_CVSINH,
    FID_VCOS,
    FID_VCOSH,
    FID_CVCOS,
    FID_CVCOSH,
    FID_HVSIN,
    FID_HVSINH,
    FID_HCVSIN,
    FID_HCVSINH,
    FID_QVSIN,
    FID_QVSINH,
    FID_QCVSIN,
    FID_QCVSINH,
    FID_HVCOS,
    FID_HVCOSH,
    FID_HCVCOS,
    FID_HCVCOSH,
    FID_QVCOS,
    FID_QVCOSH,
    FID_QCVCOS,
    FID_QCVCOSH,
    FID_ESEC,
    FID_ESECH,
    FID_ECSC,
    FID_ECSCH,
    FID_HESEC,
    FID_HESECH,
    FID_HECSC,
    FID_HECSCH,
    FID_QESEC,
    FID_QESECH,
    FID_QECSC,
    FID_QECSCH,
    FID_SINC,
    FID_SINCH,
    FID_HSINC,
    FID_HSINCH,
    FID_QSINC,
    FID_QSINCH,
    FID_COSC,
    FID_COSCH,
    FID_HCOSC,
    FID_HCOSCH,
    FID_QCOSC,
    FID_QCOSCH,
    FID_SECC,
    FID_SECCH,
    FID_HSECC,
    FID_HSECCH,
    FID_QSECC,
    FID_QSECCH,
    FID_CSCC,
    FID_CSCCH,
    FID_HCSCC,
    FID_HCSCCH,
    FID_QCSCC,
    FID_QCSCCH,
    FID_TANC,
    FID_TANCH,
    FID_HTANC,
    FID_HTANCH,
    FID_QTANC,
    FID_QTANCH,
    FID_COTC,
    FID_COTCH,
    FID_HCOTC,
    FID_HCOTCH,
    FID_QCOTC,
    FID_QCOTCH,
    FID_LOG,
    FID_LOG10,
    FID_LOG2,
    FID_LOGC,
    FID_LOG10C,
    FID_LOG2C,
    FID_LOG1P,
    FID_LOG1PC,
    FID_EXP,
    FID_EXPC,
    FID_EXP10,
    FID_EXP10C,
    FID_EXP2,
    FID_EXP2C,
    FID_ASUM,
    FID_FIBO,
    FID_FACT,
    FID_SFACT,
    FID_NPNUM,
    FID_PRIMR,
    FID_FPNSUM,
    FID_FIBNC,
    FID_FSUM,
    FID_FASUM,
    FID_SFASUM,
    FID_FNNSUM,
    FID_FLOOR,
    FID_CEIL,
    FID_DEG,
    FID_RAD
};

#define LAST_FID FID_RAD
#define MAX_FIDS LAST_FID+1


// MACRO PER FAR VEDERE O MENO LA DESCRIZIONE
// NELLE VARIE FUNZIONI DI FORMATTAZIONE INFORMAZIONI
// RELATIVE AL PROGRAMMA.
#define WITH_DESCRIPTION true
#define WITHOUT_DESCRIPTION false


// DEFINING MAX_TYPES
// and DOMAIN_NULL MACROS
// in order to identify
// the operations relative domains.

// #define MAX_TYPES 15
// #define DOMAIN_NULL MAX_TYPES+1

// DEFINING binary_function
// or unary_function operations INVERSE_OPSs
#define unary_function false
#define binary_function true

// calcolatoreDiBase() GENERIC STUFFS...
#define BASECALC_CAMBIAMENTODIBASE_MINBASE 2
#define BASECALC_CAMBIAMENTODIBASE_MAXBASE 10

#define MAX_MINBASE_CONVERTIBLE_NUM 262000

// #define MAX_TYPES 14
// #define DOMAIN_NULL MAX_TYPES+1

#define LOWER_TRIANGULAR 0
#define UPPER_TRIANGULAR 1


enum
{
    SIMPSON1DIV8_RULE = 0,
    SIMPSON3DIV8_RULE,
    TRAPEZOIDAL_RULE
};


// DEFINING FAMILIES STRUCTURES MACROS
//
#define BY_USER false
#define AUTOMATIC true

//
// Macro to assert whether a subprogram is Child or Father
#define CHILD false
#define FATHER true
//

// BASIC CHECK FOR DMA (Dynamic Memory Allocation)
#define errMem(pntr, rval) if(checkErrMem(pntr)) return rval;

#define Elements_in(arrayname) (sizeof(arrayname)/sizeof(arrayname *))
// The macro upwriten is useless because when using dynamic mallocation,
// size_t sizes indexing of matrix are lost.


/*
Dedicati alla modalità di funzionamento
della funzione matrixToVector, per decidere
in che senso deve essere svolta
0 -> Normale, 1 -> Viceversa
*/

#define MATRIX_TO_VECTOR false
#define VECTOR_TO_MATRIX true

/*
Dedicati alla modalità di stampa della matrice
dell'omonima funzione. Come suggeriscono le
stesse macro, passando 0 si stampa una matrice
di valori in virgola mobile, altrimenti di interi.
*/

enum
{
    RAWS = 0,
    COLUMNS,
    COLUMNS2
};
/*
Dedicati alla modalita' booleana di funzionamento
del sottoprogramma procedurale binaryAlgSum.
BAS_SUM indica una somma, BAS_SUB indica una differenza,
sempre da operare ovviamente tramite operazioni strettamente
associate al tipo int, ovvero normale Somma Algerica
*/

#define BAS_SUM true
#define BAS_SUB false



// #define DIM_BITFIELD (sizeof(dim_typ)-1)*8
#define DIM_BITFIELD 8
#define DATE_BITFIELD 5
#define DOMAIN_BITFIELD 4
#define BOOL_BITFIELD 1
#define BITMASK_BITFIELD 24

#define MAX_DAYSWEEK 7
#define MAX_MONTH_DAYS 31
#define MAX_MONTHS 12
#define GETDAYNAME_MAP_LENGTH 12
#define DAYSWEEK_LENGTH SIGN_STRING
#define MONTHS_LENGTH SIGN_STRING

#define MIN_NEWTON_DIFFTABLES_DIM 1
#define MAX_NEWTON_DIFFTABLES_DIM 15


// USED TO SHOW BACKWARD && FORWARD
// NEWTON's DIFFERENCE TABLES
#define BACKWARD_DIFFTAB false
#define FORWARD_DIFFTAB true

enum
{
    EXITHANDLE_EXIT = 0,
    EXITHANDLE_GETCMD,
    EXITHANDLE_SETCMD,
    EXITHANDLE_BACKCMD,
    INVALID_EXITHANDLE
};

// #define INVALID_EXITHANDLE 3

#define PARSER_NOSETTINGS 0

#define PARSER_SHOWRESULT 1
#define PARSER_SHOWVARLIST 2
#define PARSER_SHOWDIFFTIME 4
#define PARSER_SAVERESULT 8


#define LEFT_OPR 0
#define RIGHT_OPR 1

#define MAX_RIGHE_PER_COLONNE (int)(floor(sqrt(access(curLayout)->matrix_max_raws*access(curLayout)->matrix_max_columns))) // DIMENSIONE MASSIMA MATRICI QUADRATE nXn

#define MAX_DIMENSIONS 2
#define MAX_ABSTRACT_DIMENSIONS 3

#define SUMMATION_SUM false
#define SUMMATION_SUB true

#define PRODUCTORY_MUL false
#define PRODUCTORY_DIV true


// used for CURVEs FITTING or INTERPOLATIONS SubPrograms...
//
#define XRAW 0
#define YRAW 1
//

#define FIRST_NUMBER XRAW
#define SECOND_NUMBER YRAW

#define FIRST_VECTOR XRAW
#define SECOND_VECTOR YRAW

#define FIRST_QUARTILE XRAW
#define THIRD_QUARTILE YRAW

#define FIRST_QUARTILE_CONSTANT 0.25
#define SECOND_QUARTILE_CONSTANT 0.50
#define THIRD_QUARTILE_CONSTANT 0.75

// COMPLEX Numbers MACROS...
//
#define REAL_PART XRAW
#define IMAG_PART YRAW
//
#define MAX_COMPLEX_UNITS MAX_DIMENSIONS


// IperCOMPLEX Numbers MACROS
//
enum
{
    QUATERNIONS_REALPART = 0,
    QUATERNIONS_IPART,
    QUATERNIONS_JPART,
    QUATERNIONS_KPART
};

#define LAST_QUATERNIONS_UNIT QUATERNIONS_KPART
#define MAX_QUATERNIONS_UNITS LAST_QUATERNIONS_UNIT+1

enum
{
    OCTONIONS_REALPART = 0,
    OCTONIONS_E1PART,
    OCTONIONS_E2PART,
    OCTONIONS_E3PART,
    OCTONIONS_E4PART,
    OCTONIONS_E5PART,
    OCTONIONS_E6PART,
    OCTONIONS_E7PART
};

#define LAST_OCTONIONS_UNIT OCTONIONS_E7PART
#define MAX_OCTONIONS_UNITS LAST_OCTONIONS_UNIT+1

enum
{
    SEDENIONS_REALPART = 0,
    SEDENIONS_E1PART,
    SEDENIONS_E2PART,
    SEDENIONS_E3PART,
    SEDENIONS_E4PART,
    SEDENIONS_E5PART,
    SEDENIONS_E6PART,
    SEDENIONS_E7PART,
    SEDENIONS_E8PART,
    SEDENIONS_E9PART,
    SEDENIONS_E10PART,
    SEDENIONS_E11PART,
    SEDENIONS_E12PART,
    SEDENIONS_E13PART,
    SEDENIONS_E14PART,
    SEDENIONS_E15PART
};

#define LAST_SEDENIONS_UNIT SEDENIONS_E15PART
#define MAX_SEDENIONS_UNITS LAST_SEDENIONS_UNIT+1

// SecondGradeEquationSolver() STUFFS...
//
#define ROOT_X1 0
#define ROOT_X2 1

enum
{
    COEFF_A = 0,
    COEFF_B,
    COEFF_C
};
//

// Memoizer System STUFFS...
enum
{
	FUNCTION_FIBONACCI = 0,
	FUNCTION_FATTORIALE,
	FUNCTION_EVEN_SFATTORIALE,
	FUNCTION_ODD_SFATTORIALE
};

#define LAST_MEMOIZABLE_FUNCTION FUNCTION_ODD_SFATTORIALE
#define MAX_MEMOIZABLE_FUNCTIONS LAST_MEMOIZABLE_FUNCTION+1

#define _MAX_MEMOIZABLE_FUNCTIONS MAX_MEMOIZABLE_FUNCTIONS+1

#define MAX_STD_BUFFERS MAX_ABSTRACT_DIMENSIONS
#define _MAX_STD_BUFFERS MAX_STD_BUFFERS+1

#define INVALIDRETURNVALUE_NPRIMENUMBER 0.00
#define INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS -1.00

#define INVALIDRETURNVALUE_FIBONACCI INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS
#define INVALIDRETURNVALUE_FATTORIALE INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS
#define INVALIDRETURNVALUE_SFATTORIALE INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS

#define INVALIDRETURNVALUE_EVEN_SFATTORIALE INVALIDRETURNVALUE_SFATTORIALE
#define INVALIDRETURNVALUE_ODD_SFATTORIALE INVALIDRETURNVALUE_SFATTORIALE


// PULIZIA DEL BUFFER
#define CLEARBUFFER() fflush(stdin) // cleanBuffer()

#define TRIGONOMETRIC_DOMAIN(x) (x>=-1&&x<=1)

#define random(x) ((rand()%x)+1)
#define randomize srand(starting_random_seed=DEFAULT_RANDOM_SEED)
#define random2(x,y) randomize,random(y)
#define enhanced_random(x,y) srand((unsigned)x),random(y)


// PROJECT TYPEDEFS

#define MAX_ALIAS DEFAULT_MAX_ALIAS

typedef const dim_typ dim_typ2[MAX_DIMENSIONS];
typedef const dim_typ dim_typ3[MAX_ABSTRACT_DIMENSIONS];

typedef struct
{
    const char name[2*INFO_STRING];
    const char cmdname[SIGN_STRING];
    const char usage[INFO_STRING];
    void (*program_function)(const register sel_typ, char **);
    const bool automatic: BOOL_BITFIELD;
    const bool isFather: BOOL_BITFIELD;
} sprog;

struct operations
{
    char name[2*INFO_STRING];
    char description[2*INFO_STRING];
    char identifiers[MAX_ALIAS][MAX_IDENTIFIER_LENGTH];
    sel_typ domA: DOMAIN_BITFIELD;
    sel_typ domB: DOMAIN_BITFIELD;
    volatile bool bin_or_unary: BOOL_BITFIELD;
};
//

#define NEXT_LISTNODE false
#define PREV_LISTNODE true

typedef struct _nodelist
{
    char path[MAX_PATH_LENGTH];
    void *data; // generic pointer to "data" data.
    struct _nodelist * (ref[MAX_DIMENSIONS]);
} nodelist;

typedef struct
{
    ityp **matrix;
    dim_typ dim[MAX_DIMENSIONS];
} matrixObj;

typedef struct
{
    char *buffer;
    size_t buflen;
} logObj;

typedef struct
{
    dim_typ * memoizer;
    dim_typ current_max_index: DIM_BITFIELD;
} memObj;

typedef struct
{
	char exit_char;
    fsel_typ precision: DIM_BITFIELD;
    fsel_typ stabilizer_factor: DIM_BITFIELD;
    fsel_typ algebra: DIM_BITFIELD;
    float outlier_constant;
    dim_typ matrix_max_raws;
    dim_typ matrix_max_columns;
    dim_typ max_dsvd_iterations;
    dim_typ max_simplex_iterations;
    dim_typ max_memoizable_indeces[MAX_MEMOIZABLE_FUNCTIONS];
    dim_typ min_stirlingrequires_number: DIM_BITFIELD;
    dim_typ basecalc_minbase: DIM_BITFIELD;
    dim_typ basecalc_maxbase: DIM_BITFIELD;
    int max_changebase_binary_convnum;
    dim_typ min_newton_difftables_dim: DIM_BITFIELD;
    dim_typ max_newton_difftables_dim: DIM_BITFIELD;
    dim_typ min_roman_number: DIM_BITFIELD;
    dim_typ max_roman_number;
    dim_typ pascal_triangle_min_raws: DIM_BITFIELD;
    dim_typ pascal_triangle_max_raws: DIM_BITFIELD;
    volatile unsigned bools: BITMASK_BITFIELD;
} layoutObj;

typedef struct
{
    exprValList *var_list;
    EXPRTYPE *e_ANS;
    EXPRTYPE global_var;
} exprType;

enum
{
    ENVS = 0,
    MATRICES,
    LOGS,
    LAYOUTS
};

#define LAST_PROGRAM_ITEMTYPE LAYOUTS
#define MAX_LISTS LAST_PROGRAM_ITEMTYPE+1

// ext_math.C STRUCTS DEFINITIONS

/*
DEFINING MAX_TYPES
and DOMAIN_NULL MACROS
in order to identify
the operations relative domains.
*/
//
enum
{
    DOMAIN_BOOL = 0,
    DOMAIN_UCHAR,
    DOMAIN_SCHAR,
    DOMAIN_INT,
    DOMAIN_SHRT,
    DOMAIN_LONG,
    DOMAIN_LLONG,
    DOMAIN_UINT,
    DOMAIN_USHRT,
    DOMAIN_ULONG,
    DOMAIN_ULLONG,
    DOMAIN_FLT,
    DOMAIN_DBL,
    DOMAIN_LDBL
};

#define LAST_TYPEDOMAIN DOMAIN_LDBL
#define MAX_TYPES LAST_TYPEDOMAIN+1
#define DOMAIN_NULL MAX_TYPES
//

typedef struct
{
    const long double min;
    const long double max;
} TypeRange;


struct ext_math_type
{
    const TypeRange types_range[MAX_TYPES];
    const char funcnames[MAX_FIDS][MIN_STRING];
    ityp (* const functions[MAX_FIDS])(ityp);
    const char days_week_names[MAX_DAYSWEEK][DAYSWEEK_LENGTH];
    const char months_names[MAX_MONTHS][MONTHS_LENGTH];
    const sel_typ getdaynum_map[GETDAYNAME_MAP_LENGTH];
    const char romn_thousand_map[MAX_ABSTRACT_DIMENSIONS][ROMAN_NUMBER_MAPSTRING];
    const char romn_map[MAX_ABSTRACT_DIMENSIONS][9][ROMAN_NUMBER_MAPSTRING];
};
//

struct prog_constants
{
    char listsnames[MAX_LISTS][SIGN_STRING];
    #ifdef XMLCALL
    	const char bools_identifiers[MAX_DIMENSIONS][SIGN_STRING];
    #endif
    const char bools_names[MAX_BOOL_SETTINGS][INFO_STRING];
    struct
    {
        const unsigned bmask;
        const bool default_val;
    } const bools[MAX_BOOL_SETTINGS];
    const dim_typ max_memoizable_indices[MAX_MEMOIZABLE_FUNCTIONS];
	const char memoizers_names[_MAX_MEMOIZABLE_FUNCTIONS][INFO_STRING];
    const char algebra_elements_names[_MAX_ALGEBRA][INFO_STRING];
    const char algebra_imaginary_units_names[_MAX_ALGEBRA][MAX_SEDENIONS_UNITS][SIGN_STRING];
    void (* matrixSumFuncs[_MAX_ALGEBRA])(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
    void (* matrixProdFuncs[_MAX_ALGEBRA])(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_ABSTRACT_DIMENSIONS]);
    void (* matrixKProdFuncs[_MAX_ALGEBRA])(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS]);
    const char colors_types_names[MAX_COLOR_TYPES][INFO_STRING];
    const char colors_names[MAX_COLORS][INFO_STRING];
};

struct program
{
    char sysLogPath[MAX_PATH_LENGTH];
    void (* matrixSumFunc)(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
    void (* matrixProdFunc)(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_ABSTRACT_DIMENSIONS]);
    void (* matrixKProdFunc)(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS]);
    exprType *exprVars;
    matrixObj *curMatrix;
    matrixObj *lmpMatrix;
    logObj *curLog;
    logObj *sysLog;
    layoutObj *curLayout;
    memObj sysMem[MAX_MEMOIZABLE_FUNCTIONS];
    struct
    {
        nodelist * (items_list[MAX_DIMENSIONS]);
        dim_typ cur_item: DIM_BITFIELD;
        dim_typ itemsno: DIM_BITFIELD;
    } lists[MAX_LISTS];
    exprFuncList *func_list;
    exprValList *const_list;
    char colors_path[MAX_PATH_LENGTH];
    fsel_typ random_seed;
    sel_typ colors[MAX_COLOR_TYPES];
    volatile sel_typ mode;
    volatile sel_typ exitHandle;
    volatile bool mss: BOOL_BITFIELD;
    volatile bool sigresult: BOOL_BITFIELD;
};

#ifdef STACKALLOC
    #define access(x) suite.x
    extern struct program suite;
#else
    #define access(x) suite->x
    extern struct program * suite;
#endif


#define INVERSE_OPS isSett(BOOLS_INVERSEOPERATIONS)
#define getVarList(fp) getVarListEx(fp, access(lists)[ENVS].cur_item)


// extern nodelist *env_list;

extern struct operations operazioni[MAX_OPERATIONS];

_CRTIMP extern int errno;

/// extern struct program suite;

extern sprog main_menu[MAX_PROGRAMMI],
                    adv_calc[MAX_ADVCALC_PROGS],
                    alg_operations[MAX_ALGEBRA_OPERATIONS];

#ifdef ALLOW_MSSMANAGER
    extern sprog mss_manager[MAX_MSSMANAGER_PROGS];
#endif // ALLOW_MSSMANAGER

#ifndef FREEZE_SETTINGSMANAGER
    extern sprog change_settings[MAX_SETTINGS];
#endif

#ifdef ALLOW_VARLISTMANAGER
    extern sprog envs_manager[MAX_ENVSMANAGER_PROGS];
#endif

#ifdef ALLOW_MATMANAGER
    extern sprog mat_manager[MAX_MATMANAGER_PROGS];
#endif

#ifdef ALLOW_LOGSMANAGER
    extern sprog logs_manager[MAX_LOGSMANAGER_PROGS];
#endif

#ifdef ALLOW_SYSLOGMANAGER
    extern sprog syslog_manager[MAX_SYSLOGMANAGER_PROGS];
#endif

#ifdef ALLOW_LAYOUTSMANAGER
    extern sprog layouts_manager[MAX_LAYOUTSMANAGER_PROGS];
#endif

#ifdef ALLOW_COLSMANAGER
    extern sprog cols_manager[MAX_COLSMANAGER_PROGS];
#endif

#ifdef ALLOW_LFSMANAGER
    extern sprog lfs_manager[MAX_LFSMANAGER_PROGS];
#endif

extern const struct ext_math_type ext_math;
extern const struct prog_constants suite_c;


/// FUNCTIONS DECLARATIONS
/// geometry.c
__MSUTIL_ void toupper_s(char *);
__MSUTIL_ void tolower_s(char *);
__MSNATIVE_ void strundsc(const char *, char []);
__MSNATIVE_ void strboolize(const char *, char []);
__MSNATIVE_ void strfnm(const char *, char [static MAX_PATH_LENGTH]);
__MSUTIL_ int __export countbits(long);
__MSUTIL_ int __export ucountbits(unsigned long);
__MSUTIL_ char __export *strrev(char *);
__MSUTIL_ char __export *replace(char const * const, char const * const, char const * const);
__MSUTIL_ bool __export __system file_exists(const char *);
__MSNATIVE_ bool __system readFile(const char [static MAX_PATH_LENGTH]);
__MSNATIVE_ bool __system printFile(const char [static MAX_PATH_LENGTH]);
__MSNATIVE_ bool __system writeFile(const char [static MAX_PATH_LENGTH]);
__MSNATIVE_ FILE * __system checkForFHErrors(const char [static MAX_PATH_LENGTH], char [static 1]);
__MSNATIVE_ bool __system frename(const char [static MAX_PATH_LENGTH], const char [static MAX_PATH_LENGTH]);
__MSUTIL_ bool matrixUTConv(ityp **, dim_typ);
__MSNATIVE_ const ityp sarrus(ityp **);
__MSNATIVE_ ityp carlucci(ityp **);
__MSNATIVE_ ityp checkStdMat(ityp **, dim_typ);
__MSNATIVE_ __MSUTIL_ ityp det(ityp **, dim_typ, bool *);
__MSNATIVE_ ityp trace(ityp **, dim_typ);
__MSSHELL_WRAPPER_ __MSNATIVE_ bool randomMatrix(ityp **, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ void transpose(ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]);
__MSUTIL_ bool __export FattLU(dim_typ, ityp **, ityp **, ityp **);
__MSUTIL_ bool __export invertMatrix(ityp **, dim_typ);
__MSUTIL_ bool __export dsvd(ityp **, const register dim_typ [static MAX_DIMENSIONS], ityp *, ityp **);
__MSNATIVE_ dim_typ __export rank(ityp **, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ void __export matrixToVector(ityp **, const register dim_typ dim[static MAX_DIMENSIONS], ityp [static dim[RAWS]*dim[COLUMNS]], bool);
__MSNATIVE_ void __export vectorToMatrix(const register dim_typ dim[static MAX_DIMENSIONS], ityp [static dim[RAWS]*dim[COLUMNS]], ityp **);
__MSNATIVE_ void __system _flushLogBuf(logObj * const);
__MSNATIVE_ __WINCALL void __system _editLog(const char [static MAX_PATH_LENGTH]);
__MSNATIVE_ void __system logCheck(logObj * const , const char *, const char [static MAX_PATH_LENGTH]);
__MSNATIVE_ void __system prependTimeToString(char *, const bool);
__MSNATIVE_ void __system printErr(const int, const char *, ...);
__MSNATIVE_ void __system sprint(const char *, ...);
__MSNATIVE_ void __system fprintf2(FILE *, const char *, ...);
__MSNATIVE_ void __system printf2(const sel_typ, const char *, ...);
__MSNATIVE_ bool __system scanf2(sel_typ, const char *, ...);
__MSNATIVE_ bool __system __export parse(char [], ityp *);
__MSNATIVE_ _MS__private void __system printMatrix(FILE *, ityp **, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ bool __system __export extractMat(dim_typ);
__MSNATIVE_ bool __system __export matrixToken(const char [], ityp ***, dim_typ *, dim_typ *);
__MSSHELL_WRAPPER_ __MSNATIVE_ const sprog * const __system searchProgram(const char [static SIGN_STRING]);
__MSUTIL_ int __export cmpfunc(const void *, const void *);
__MSUTIL_ bool __export min_cmpfunc(const register ityp, const register ityp);
__MSUTIL_ bool __export max_cmpfunc(const register ityp, const register ityp);
__MSNATIVE_ void __system _changeAlgebraDims(const dim_typ);
__MSNATIVE_ void __system _showUsage(const sprog * const);
__MSNATIVE_ void __system printUsage(const sprog * const);
__MSNATIVE_ void __system prepareToExit(void);
__MSNATIVE_ void __system safeExit(const volatile int);
__MSNATIVE_ void __system _handleCmdLine(const register sel_typ, char **);
__MSNATIVE_ bool __system _execScriptFiles(const char [static MAX_PATH_LENGTH]);
__MSNATIVE_ bool __system _lfLoader(const char [static MAX_PATH_LENGTH]);
__MSNATIVE_ bool __system _lfCreate(const char [static MAX_PATH_LENGTH]);

#ifdef XMLCALL
	__MSUTIL_ XMLCALL xmlDoc * __system __export xmlInit(const char [static XML_FILENAMES_LENGTH], xmlXPathContext **);
	__MSUTIL_ XMLCALL void __system __export xmlExit(const char [static XML_FILENAMES_LENGTH], xmlDoc **, xmlXPathObject **, xmlXPathContext **);
	__MSUTIL_ XMLCALL void __system __export xmlWriteInt(xmlXPathObject **, xmlXPathContext *, const char *, const int);
	__MSUTIL_ XMLCALL void __system __export xmlWriteBool(xmlXPathObject **, xmlXPathContext *, const char *, const bool);
	__MSUTIL_ XMLCALL void __system __export xmlWriteFloat(xmlXPathObject **, xmlXPathContext *, const char *, const float);
	__MSUTIL_ XMLCALL void __system __export xmlWriteString(xmlXPathObject **, xmlXPathContext *, const char *, const char [static MAX_XML_FIELDSTRINGS]);
	__MSUTIL_ XMLCALL int __system __export xmlGetInt(xmlXPathObject **, xmlXPathContext *, const char *);
	__MSUTIL_ XMLCALL bool __system __export xmlGetBool(xmlXPathObject **, xmlXPathContext *, const char *);
	__MSUTIL_ XMLCALL float __system __export xmlGetFloat(xmlXPathObject **, xmlXPathContext *, const char *);
	__MSUTIL_ XMLCALL void __system __export xmlGetString(xmlXPathObject **, xmlXPathContext *, const char *, char [static MAX_XML_FIELDSTRINGS]);
#endif

#if WINOS
    __MSUTIL_ __WINCALL BOOL WINAPI __system __export SetExitButtonState(const bool);
    __MSUTIL_ const char * const __system __export getFilename(const char [static MAX_PATH_LENGTH]);
    __MSUTIL_ __WINCALL void __system updInfo(void);
    __MSUTIL_ __WINCALL HWND WINAPI __system __export GetConsoleWindowNT();
    __MSUTIL_ __WINCALL bool __system __export windowsFileHandler(char *, const char *, const char [static MAX_EXTENSION_LENGTH], bool);
#endif

__MSUTIL_ void __system updInfo(void);
__MSUTIL_ void __system __export SetColor(const sel_typ);
__MSNATIVE_ XMLCALL void __system _backupColFile(void);
__MSNATIVE_ XMLCALL void __system getProgramSettings(dim_typ);
__MSNATIVE_ XMLCALL void __system _colFileLoader(const char [static MAX_PATH_LENGTH]);

__MSNATIVE_ ityp __export MINMAX(const register dim_typ dim, const ityp [static dim], const bool, dim_typ *);
__MSNATIVE_ bool __system __export isDomainForbidden(ityp, bool);
__MSUTIL_ void _MS__private __system __export free_foreach(ityp ***, const dim_typ, const dim_typ, bool mode);
__MSUTIL_ void _MS__private __system __export free_foreach2(ityp ***, const dim_typ);
__MSUTIL_ void _MS__private __system __export free_foreach3(ityp ***, const dim_typ);
__MSUTIL_ bool __system __export checkErrMem(const void *);
__MSNATIVE_ bool __system __export matrixAlloc(ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ void __system __export _matrixFree(ityp ***, dim_typ, bool mode);
__MSNATIVE_ bool __system __export equalMatrix(ityp ***, ityp **, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ _MS__private void __system resetLmpMatrix(void);
__MSNATIVE_ void __system __export _flushMemoizersBuffers(sel_typ);
__MSNATIVE_ void __system __export flushAllMemoizersBuffers(void);
__MSNATIVE_ dim_typ __system __export selectListItem(dim_typ dim, bool, const char *, const char [static dim][MIN_STRING]);
__MSNATIVE_ void __system viewProgramSettings(dim_typ);
__MSNATIVE_ XMLCALL void __system resetProgramSettings(layoutObj * const, const char [static MAX_PATH_LENGTH]);
__MSNATIVE_ void __system setProgramSettings(dim_typ);
__MSNATIVE_ bool __system __export catchPause();
__MSNATIVE_ void __system logPrint(logObj * const);
__MSNATIVE_ void __system logWrite(FILE *, logObj * const);
__MSNATIVE_ void __system getVarListEx(FILE *, dim_typ);
__MSNATIVE_ ityp __system __export requires(const char *, const char *, const char *, const unsigned);
__MSNATIVE_ bool __system insertDims(dim_typ *, dim_typ *);
__MSNATIVE_ bool __system insertDim(dim_typ *, bool);
__MSNATIVE_ void __system sigproc(void);
__MSNATIVE_ void __system sigexit(void);
__MSSHELL_WRAPPER_ __MSNATIVE_ void showNewtonDifferenceTable(dim_typ n, ityp [static n], ityp [access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim], bool);
__MSNATIVE_ bool __system insertNMMatrix(ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ volatile char __system insertElement(ityp **, const register dim_typ [static MAX_DIMENSIONS], bool);
__MSNATIVE_ volatile bool __system checkBackTracking(volatile char, dim_typ *);
__MSNATIVE_ volatile __system sel_typ checkBackTracking2(volatile char, dim_typ *, dim_typ *, dim_typ *, dim_typ);
__MSNATIVE_ bool __system __export enterMatrix(ityp ***, dim_typ *, dim_typ *, bool, bool);


/// programs.c
__MSSHELL_WRAPPER_ void __apnt changeProgramSettings(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void calcolatoreDiBase(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void __apnt calcolatoreAvanzato(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void __apnt algebraOperations(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void __apnt mssManager(const register sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __MSNATIVE_ void _MS__private __system __export operationsGroupMenu(dim_typ dim, sprog [static dim], const char [static INFO_STRING], bool);
__MSSHELL_WRAPPER_ __MSNATIVE_ void __system progInfo(sel_typ);
__MSNATIVE_ bool _MS__private __system doesExistOperIdentifier(const char [static MAX_IDENTIFIER_LENGTH], sel_typ [static MAX_DIMENSIONS]);
__MSUTIL_ void __system __export freeExprEvalLists();
__MSUTIL_ void __system __export refreshExprEvalVarList(dim_typ);
__MSUTIL_ void __system __export refreshExprEvalLists();
__MSNATIVE_ void __system __export setCurrentMatrix(dim_typ);
__MSSHELL_WRAPPER_ __MSNATIVE_ void __system setDefaults();



/// ext_math.c

__MSNATIVE_ char * const __system __export binaryAlgSum(ityp, ityp, bool);
__MSNATIVE_ char * const __system __export binNumComp(ityp);
__MSNATIVE_ void __export _complexSum(ityp **, ityp [static MAX_DIMENSIONS]);
__MSNATIVE_ void __export _complexProd(ityp **, ityp [static MAX_DIMENSIONS]);
__MSNATIVE_ void __export _quaternionsSum(ityp **, ityp [static MAX_QUATERNIONS_UNITS]);
__MSNATIVE_ void __export _quaternionsProduct(ityp **, ityp [static MAX_QUATERNIONS_UNITS]);
__MSNATIVE_ void __export _octonionsEMTSum(ityp **, ityp [static MAX_OCTONIONS_UNITS]);
__MSNATIVE_ void __export _octonionsMTSum(ityp **, ityp [static MAX_OCTONIONS_UNITS]);
__MSNATIVE_ void __export _octonionsSum(ityp **, ityp [static MAX_OCTONIONS_UNITS]);
__MSNATIVE_ void __export _octonionsEMTProduct(ityp **, ityp [static MAX_OCTONIONS_UNITS]);
__MSNATIVE_ void __export _octonionsMTProduct(ityp **, ityp [static MAX_OCTONIONS_UNITS]);
__MSNATIVE_ void __export _octonionsProduct(ityp **, ityp [static MAX_OCTONIONS_UNITS]);
__MSNATIVE_ void __export _sedenionsEMTSum(ityp **, ityp [static MAX_SEDENIONS_UNITS]);
__MSNATIVE_ void __export _sedenionsMTSum(ityp **, ityp [static MAX_SEDENIONS_UNITS]);
__MSNATIVE_ void __export _sedenionsSum(ityp **, ityp [static MAX_SEDENIONS_UNITS]);
__MSNATIVE_ void __export _sedenionsEMTProduct(ityp **, ityp [static MAX_SEDENIONS_UNITS]);
__MSNATIVE_ void __export _sedenionsMTProduct(ityp **, ityp [static MAX_SEDENIONS_UNITS]);
__MSNATIVE_ void __export _sedenionsProduct(ityp **, ityp [static MAX_SEDENIONS_UNITS]);
__MSNATIVE_ bool __export _secondGradeEquationSolver(ityp *, ityp [static MAX_DIMENSIONS]);
__MSUTIL_ const char * const __export getDayName(sel_typ);
__MSUTIL_ const char * const __export getMonthName(sel_typ);
__MSUTIL_ const sel_typ __export getDayNumber(sel_typ, sel_typ, uint64_t);
__MSUTIL_ void __export getRomanNumber(dim_typ, char [static ROMAN_NUMBER_STRING]);
__MSUTIL_ int64_t _MS__private __export powi(register int64_t, register int64_t);
__MSUTIL_ ityp __export mpow(ityp, int64_t);
__MSUTIL_ ityp __export mpow2(ityp, int64_t);
__MSNATIVE_ __MSUTIL_ void __export getPascalTriangle(uint64_t);
__MSUTIL_ bool __system __export isPrimeHISP(register uint64_t);
__MSUTIL_ bool __system __export isPrimeLIFP(register uint64_t);
__MSNATIVE_ void __system __export prime_N_Number(register uint64_t, register uint64_t);
__MSNATIVE_ ityp __system __export N_prime_Number(register ityp);
__MSNATIVE_ ityp __export fibnc(register ityp);
__MSNATIVE_ ityp __export primr(register ityp);
__MSNATIVE_ ityp __export fpnsum(register ityp);
__MSNATIVE_ bool __system __export trigonometric_domain(register ityp, register ityp, register ityp);
__MSNATIVE_ ityp __export stirling(register ityp);
__MSNATIVE_ ityp __export fibo(register ityp);
__MSNATIVE_ ityp __export fact(register ityp);
__MSUTIL_ ityp __export sfact(register ityp);
__MSNATIVE_ ityp __export sfact_even(register ityp);
__MSNATIVE_ ityp __export sfact_odd(register ityp);
__MSUTIL_ ityp __export perm(register ityp);
__MSUTIL_ ityp __export perm_rep(register ityp, register uint64_t dim, ityp [static dim]);
__MSUTIL_ ityp __export kperm(register ityp, register ityp);
__MSUTIL_ ityp __export kperm_rep(register ityp, register ityp);
__MSUTIL_ ityp __export comb(register ityp, register ityp);
__MSUTIL_ ityp __export comb_rep(register ityp, register ityp);
__MSNATIVE_ ityp __export gsum(register ityp, register int64_t);
__MSNATIVE_ ityp __export gasum(register uint64_t, register uint64_t);
__MSNATIVE_ ityp __export asum(register ityp);
__MSNATIVE_ ityp __export fsum(register ityp);
__MSNATIVE_ ityp __export fasum(register ityp);
__MSNATIVE_ ityp __export sfasum(register ityp);
__MSNATIVE_ ityp __export fnnsum(register ityp);
__MSNATIVE_ ityp __export summation(uint64_t dim, bool, ityp [static dim]);
__MSNATIVE_ ityp __export productory(uint64_t dim, bool, ityp [static dim]);
__MSNATIVE_ ityp __export math_media(uint64_t dim, ityp[static dim]);
__MSNATIVE_ ityp __export math_variance(uint64_t dim, ityp [static dim]);
__MSNATIVE_ ityp __export math_covariance(uint64_t dim, ityp [static dim], ityp [static dim]);
__MSNATIVE_ ityp __export math_stddev(uint64_t dim, ityp [static dim]);
__MSNATIVE_ bool __export math_outlier(uint64_t dim, uint64_t, ityp [static dim]);
__MSNATIVE_ bool __export math_outlier2(uint64_t dim, uint64_t, float, ityp vector[static dim]);
__MSNATIVE_ ityp __export math_geomedia(uint64_t dim, ityp[static dim]);
__MSNATIVE_ ityp __export math_armedia(uint64_t dim, ityp [static dim]);
__MSNATIVE_ ityp __export math_powmedia(uint64_t dim, uint64_t, ityp [static dim]);
__MSNATIVE_ ityp __export math_scale(uint64_t dim, ityp [static dim]);
__MSNATIVE_ ityp __export math_first_quartile(uint64_t dim, ityp [static dim]);
__MSNATIVE_ ityp __export math_mediana(uint64_t dim, ityp [static dim]);
__MSNATIVE_ ityp __export math_third_quartile(uint64_t dim, ityp [static dim]);
__MSUTIL_ int64_t __system __export changeBase(register int, sel_typ, sel_typ);
__MSUTIL_ uint64_t __export math_MCD(register uint64_t, register uint64_t);
__MSUTIL_ uint64_t __export math_mcm(register uint64_t, register uint64_t);
__MSNATIVE_ ityp __export exp10(register ityp);
__MSNATIVE_ ityp __export expc(register ityp);
__MSNATIVE_ ityp __export exp10c(register ityp);
__MSNATIVE_ ityp __export exp2c(register ityp);
__MSUTIL_ ityp __export logbN(register ityp, register ityp);
__MSNATIVE_ ityp __export logc(register ityp);
__MSNATIVE_ ityp __export log10c(register ityp);
__MSNATIVE_ ityp __export log2c(register ityp);
__MSNATIVE_ ityp __export log1pc(register ityp);
__MSUTIL_ ityp __export rootnX(register ityp, register ityp);
__MSNATIVE_ ityp __export deg(register ityp);
__MSNATIVE_ ityp __export rad(register ityp);
__MSNATIVE_ ityp __export csc(register ityp);
__MSNATIVE_ ityp __export sec(register ityp);
__MSNATIVE_ ityp __export cot(register ityp);
__MSNATIVE_ ityp __export csch(register ityp);
__MSNATIVE_ ityp __export sech(register ityp);
__MSNATIVE_ ityp __export coth(register ityp);
__MSNATIVE_ ityp __export acsc(register ityp);
__MSNATIVE_ ityp __export asec(register ityp);
__MSNATIVE_ ityp __export acot(register ityp);
__MSNATIVE_ ityp __export acsch(register ityp);
__MSNATIVE_ ityp __export asech(register ityp);
__MSNATIVE_ ityp __export acoth(register ityp);
__MSNATIVE_ ityp __export hsin(register ityp);
__MSNATIVE_ ityp __export hsinh(register ityp);
__MSNATIVE_ ityp __export qsin(register ityp);
__MSNATIVE_ ityp __export qsinh(register ityp);
__MSNATIVE_ ityp __export hcos(register ityp);
__MSNATIVE_ ityp __export hcosh(register ityp);
__MSNATIVE_ ityp __export qcos(register ityp);
__MSNATIVE_ ityp __export qcosh(register ityp);
__MSNATIVE_ ityp __export hcsc(register ityp);
__MSNATIVE_ ityp __export hcsch(register ityp);
__MSNATIVE_ ityp __export qcsc(register ityp);
__MSNATIVE_ ityp __export qcsch(register ityp);
__MSNATIVE_ ityp __export hsec(register ityp);
__MSNATIVE_ ityp __export hsech(register ityp);
__MSNATIVE_ ityp __export qsec(register ityp);
__MSNATIVE_ ityp __export qsech(register ityp);
__MSNATIVE_ ityp __export htan(register ityp);
__MSNATIVE_ ityp __export htanh(register ityp);
__MSNATIVE_ ityp __export qtan(register ityp);
__MSNATIVE_ ityp __export qtanh(register ityp);
__MSNATIVE_ ityp __export hcot(register ityp);
__MSNATIVE_ ityp __export hcoth(register ityp);
__MSNATIVE_ ityp __export qcot(register ityp);
__MSNATIVE_ ityp __export qcoth(register ityp);
__MSNATIVE_ ityp __export vsin(register ityp);
__MSNATIVE_ ityp __export cvsin(register ityp);
__MSNATIVE_ ityp __export vcos(register ityp);
__MSNATIVE_ ityp __export cvcos(register ityp);
__MSNATIVE_ ityp __export hvsin(register ityp);
__MSNATIVE_ ityp __export hcvsin(register ityp);
__MSNATIVE_ ityp __export hvcos(register ityp);
__MSNATIVE_ ityp __export hcvcos(register ityp);
__MSNATIVE_ ityp __export qvsin(register ityp);
__MSNATIVE_ ityp __export qcvsin(register ityp);
__MSNATIVE_ ityp __export qvcos(register ityp);
__MSNATIVE_ ityp __export qcvcos(register ityp);
__MSNATIVE_ ityp __export vsinh(register ityp);
__MSNATIVE_ ityp __export cvsinh(register ityp);
__MSNATIVE_ ityp __export vcosh(register ityp);
__MSNATIVE_ ityp __export cvcosh(register ityp);
__MSNATIVE_ ityp __export hvsinh(register ityp);
__MSNATIVE_ ityp __export hcvsinh(register ityp);
__MSNATIVE_ ityp __export hvcosh(register ityp);
__MSNATIVE_ ityp __export hcvcosh(register ityp);
__MSNATIVE_ ityp __export qvsinh(register ityp);
__MSNATIVE_ ityp __export qcvsinh(register ityp);
__MSNATIVE_ ityp __export qvcosh(register ityp);
__MSNATIVE_ ityp __export qcvcosh(register ityp);
__MSNATIVE_ ityp __export esec(register ityp);
__MSNATIVE_ ityp __export ecsc(register ityp);
__MSNATIVE_ ityp __export esech(register ityp);
__MSNATIVE_ ityp __export ecsch(register ityp);
__MSNATIVE_ ityp __export hesec(register ityp);
__MSNATIVE_ ityp __export hecsc(register ityp);
__MSNATIVE_ ityp __export hesech(register ityp);
__MSNATIVE_ ityp __export hecsch(register ityp);
__MSNATIVE_ ityp __export qesec(register ityp);
__MSNATIVE_ ityp __export qecsc(register ityp);
__MSNATIVE_ ityp __export qesech(register ityp);
__MSNATIVE_ ityp __export qecsch(register ityp);
__MSNATIVE_ ityp __export sinc(register ityp);
__MSNATIVE_ ityp __export sinch(register ityp);
__MSNATIVE_ ityp __export hsinc(register ityp);
__MSNATIVE_ ityp __export hsinch(register ityp);
__MSNATIVE_ ityp __export qsinc(register ityp);
__MSNATIVE_ ityp __export qsinch(register ityp);
__MSNATIVE_ ityp __export cosc(register ityp);
__MSNATIVE_ ityp __export cosch(register ityp);
__MSNATIVE_ ityp __export hcosc(register ityp);
__MSNATIVE_ ityp __export hcosch(register ityp);
__MSNATIVE_ ityp __export qcosc(register ityp);
__MSNATIVE_ ityp __export qcosch(register ityp);
__MSNATIVE_ ityp __export secc(register ityp);
__MSNATIVE_ ityp __export secch(register ityp);
__MSNATIVE_ ityp __export hsecc(register ityp);
__MSNATIVE_ ityp __export hsecch(register ityp);
__MSNATIVE_ ityp __export qsecc(register ityp);
__MSNATIVE_ ityp __export qsecch(register ityp);
__MSNATIVE_ ityp __export cscc(register ityp);
__MSNATIVE_ ityp __export cscch(register ityp);
__MSNATIVE_ ityp __export hcscc(register ityp);
__MSNATIVE_ ityp __export hcscch(register ityp);
__MSNATIVE_ ityp __export qcscc(register ityp);
__MSNATIVE_ ityp __export qcscch(register ityp);
__MSNATIVE_ ityp __export tanc(register ityp);
__MSNATIVE_ ityp __export tanch(register ityp);
__MSNATIVE_ ityp __export htanc(register ityp);
__MSNATIVE_ ityp __export htanch(register ityp);
__MSNATIVE_ ityp __export qtanc(register ityp);
__MSNATIVE_ ityp __export qtanch(register ityp);
__MSNATIVE_ ityp __export cotc(register ityp);
__MSNATIVE_ ityp __export cotch(register ityp);
__MSNATIVE_ ityp __export hcotc(register ityp);
__MSNATIVE_ ityp __export hcotch(register ityp);
__MSNATIVE_ ityp __export qcotc(register ityp);
__MSNATIVE_ ityp __export qcotch(register ityp);

__MSNATIVE_ bool __system __export isEqualMatrix(ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ __MSUTIL_ ityp __system __export norms(ityp **, dim_typ);
__MSNATIVE_ __MSUTIL_ ityp __system __export norm(ityp **, dim_typ, bool);
__MSNATIVE_ __MSUTIL_ void __export newtonDifferenceTable(dim_typ, ityp [access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim], bool);
__MSNATIVE_ ityp __system __export math_sum(register ityp, register ityp);
__MSNATIVE_ ityp __system __export math_sub(register ityp, register ityp);
__MSNATIVE_ ityp __system __export math_mul(register ityp, register ityp);
__MSNATIVE_ ityp __system __export math_div(register ityp, register ityp);
__MSNATIVE_ sel_typ _MS__private __system __export _simplexMethod(ityp ***, ityp ***, const register dim_typ dim [static MAX_DIMENSIONS], const bool [static dim[RAWS]-1], bool);
__MSNATIVE_ _MS__private void __system __export _matrixOEMTSum(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixOEMTSum(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixSEMTSum(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixSMTSum(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixOEMTProduct(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixOMTProduct(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixSEMTProduct(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixSMTProduct(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixKOEMTProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ, const register dim_typ, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixKOMTProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ, const register dim_typ, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixKSEMTProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ, const register dim_typ, const register dim_typ, const register dim_typ);
__MSNATIVE_ _MS__private void __system __export _matrixKSMTProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ, const register dim_typ, const register dim_typ, const register dim_typ);
__MSNATIVE_ void _MS__private __system __export _matrixSum(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixCSum(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixQSum(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixOSum(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixSSum(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_ABSTRACT_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixCProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_ABSTRACT_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixQProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_ABSTRACT_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixOProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_ABSTRACT_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixSProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_ABSTRACT_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixKProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixKCProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixKQProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixKOProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS]);
__MSNATIVE_ void _MS__private __system __export _matrixKSProduct(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS]);

/// lists_manager.c
__MSNATIVE_ sel_typ __system __export checkItemTypeByExtension(const char [static MAX_EXTENSION_LENGTH]);
__MSSHELL_WRAPPER_ __MSNATIVE_ bool _MS__private __system __export listInsertProc(const char [static MAX_PATH_LENGTH], sel_typ);
__MSNATIVE_ dim_typ __system __export searchItem(const char [static MAX_PATH_LENGTH], sel_typ);
__MSNATIVE_ nodelist * __system __export listNo(dim_typ, sel_typ);
__MSNATIVE_ bool __system __export listInsert(const char [static MAX_PATH_LENGTH], sel_typ);
__MSNATIVE_ bool __system __export listDelete(dim_typ, sel_typ);
__MSNATIVE_ dim_typ __system __export itemSelect(sel_typ);
__MSNATIVE_ void __system __export refreshItem(dim_typ, sel_typ);
__MSNATIVE_ bool __system __export saveItem(dim_typ, sel_typ);
__MSNATIVE_ void __system __export passToItem(dim_typ, sel_typ, bool);
__MSNATIVE_ void __system __export setCurItem(const dim_typ, sel_typ);
__MSNATIVE_ void __system __export createItemData(dim_typ, sel_typ);
__MSNATIVE_ __WINCALL void __system __export createItem(const char *, sel_typ);
__MSNATIVE_ void __system __export viewItem(const dim_typ, sel_typ);
__MSNATIVE_ void __system __export printListItem(const dim_typ, sel_typ);
__MSNATIVE_ void __system __export updItem(sel_typ);
__MSNATIVE_ void __system __export updAll(sel_typ);
__MSNATIVE_ void __system __export delItem(const dim_typ, sel_typ);
__MSNATIVE_ dim_typ __system __export getItemID(char *, sprog * const, sel_typ);
__MSNATIVE_ void __system __export relItem(const dim_typ, bool);
__MSNATIVE_ void __system __export renItem(const char *, const dim_typ, sel_typ);

/// logs_manager.c
__MSSHELL_WRAPPER_ void _MS__private __lmp_prog setLogBufLen(const char [static MAX_PATH_LENGTH], logObj * const, const size_t);

#ifdef __cplusplus
}
#endif // __cplusplus
