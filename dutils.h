// dutils.h
// 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be included
exclusively by main.c, mengine.c, programs.c, algebra.c and settings.c project files of my suite program!
I do not assume any responsibilities about the use with any other code-scripts.
*/

// HEADER Guard
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <errno.h>
#include <setjmp.h>
#include <search.h>
#include <dirent.h>

#include <sys/time.h>

#include <signal.h>

#include <omp.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

// #include <openssl/md2.h>
#include <openssl/md4.h>
#include <openssl/md5.h>
#include <openssl/mdc2.h>
#include <openssl/ripemd.h>
#include <openssl/whrlpool.h>
#include <openssl/sha.h>


#define false 0
#define true 1


// PROJECT HEADERS
// #include "ext_math.h"
#include "defaults.h"
#include "ExprEval/exprincl.h"
#include "hashmap/hashmap.h"

#define XMLCALL

#ifdef XMLCALL
	#include "libxml/encoding.h"
	#include "libxml/xmlwriter.h"
	#include "libxml/xpath.h"

	#define XML_ENCODING "UTF-8"

	#define SETTINGS_FILENAME "./settings.xml"

	#define XML_FILENAMES_LENGTH MINMIN_STRING
	#define MAX_XML_FIELDSTRINGS SIGN_STRING
#endif
// #include "ExprEval/expreval.h"

// #include "ExprEval/expreval.h"


#ifdef __cplusplus
extern "C" {
#endif


#define __MSSHELL_WRAPPER_
#define __MATHSUITE
#define __JBURKARDT
#define __WINCALL32
#define __WINCALL __WINCALL32
#define __apnt
#define __lmp_prog

#define PROGRAM_NAME "mathSuite"
#define PROGRAM_VERSION "10"
#define PROGRAM_AUTHOR "Marco Chiarelli"
#define PROGRAM_COAUTHOR "John Burkardt"
#define PROGRAM_LASTUPDATEDATE "05/05/2016"


// INITIALIZING EXPREVAL DEFAULT CONSTANTS
#define INIT_OBJLIST NULL
#define INIT_FUNCLIST NULL
#define INIT_CONSTLIST NULL
#define INIT_VARLIST NULL
#define INIT_ANS NULL

#define MAX_VARS 100

#define VSPACE


enum
{
	HASHMAP_ERROR = -2,
	HEAPALLOC_ERROR,
	NOERROR_EXIT
}; 


/// MINMAX powerful macros
/// designed to enhance
/// statistic media functions,
/// especially the median one.

#define MIN_MODE false
#define MAX_MODE true

#define MPFR_MIN(a,b) mpfr_less_p(a,b) ? a : b
#define MPFR_MAX(a,b) mpfr_greater_p(a,b) ? a : b

#define MIN(a,b) b<a?b:a
#define MAX(a,b) b<a?a:b

#define _MIN(a,b,c) MINMAX(a, NULL, b,c,NULL,NULL)
#define _MAX(a,b,c) MINMAX(NULL, a, b,c,NULL,NULL)

#define _MPFR_MIN(a,b,c) MPFR_MINMAX(a,NULL,b,c,NULL,NULL)
#define _MPFR_MAX(a,b,c) MPFR_MINMAX(NULL,a,b,c,NULL,NULL)

#define i4_huge 2147483647
#define r8_huge MAX_VAL

#define PREC (5e-16f)
#define ISZERO(x) (fabs(x)<PREC)

#define INVALID_SOCKET -1

#ifndef __DISABLE_SERVER
	#define SERVER_ADDRESS INADDR_LOOPBACK
	#define SERVER_LENGTH INET6_ADDRSTRLEN
	#define SERVER_PORT 49152
	#define SERVER_FAMILY AF_UNSPEC
	#define SERVER_BACKLOG 10
	#define PORT_STRLEN 6
	#define SERVER_BUFSIZE 1024
	#define WELLKNOWN_PORT 49152
	#define WELLKNOWN_SERVER INADDR_ANY
	#define DEFAULT_RESOLUTION false
	#define DEFAULT_LISTENINGSOCKETS 1
	#define MAX_LISTENINGSOCKETS 10
	#define MAX_ACCEPTABLECONNECTIONS 15
	#define MAX_SOCKETLENGTH 5
#endif

#define DEFAULT_DBSERVER "localhost"
#define DEFAULT_DBUSERNAME "root"
#define DEFAULT_DBPASSWORD "root"
#define DEFAULT_DATABASE "mathsuite"

#define QUERY_LENGTH 500
#define USER_MAXNAMELENGTH 45
#define USER_MAXSURNAMELENGTH 45
// #define USER_IDENTIFIERLENGTH 254
#define USER_MAXUSERNAMELENGTH 45
#define USER_MAXEMAILLENGTH 254
#define USER_MAXIDENTIFIERLENGTH MAX(USER_MAXUSERNAMELENGTH, USER_MAXEMAILLENGTH)

#define IDENTIFIER_USERNAME false
#define IDENTIFIER_EMAIL true

#define DEFAULT_IDENTIFIER IDENTIFIER_USERNAME

#define USER_MAXPASSWORDLENGTH 1000
#define USER_MAXFIELDLENGTH USER_MAXPASSWORDLENGTH
#define VAR_MAXVALLENGTH 1000
#define VAR_MAXNAMELENGTH 45
#define VAR_MAXFIELDLENGTH VAR_MAXVALLENGTH
#define MAT_MAXNAMELENGTH 255
#define TAB_MAXNAMELENGTH 255

#define INVALID_ID -1
#define INVALID_VARLISTID INVALID_ID
#define INVALID_MATRIXID INVALID_ID

#define DELETEMAT_ROW true
#define DELETEMAT_MAT false

enum
{
	INVERTMATRIX_SUCCESS = 0,
	INVERTMATRIX_ALLOCERROR,
	INVERTMATRIX_SINGULAR
};

/// #define INI_NAME "./settings.xml"

/// #define WINOS true /// ((defined(WIN32)) || (defined(WINDOWS)))

#define XML_NAME "./settings.xml"

#ifdef WINOS
	// #define WIN32_LEAN_AND_MEAN
	#include <winsock2.h>
	// #include <windows.h>
	#include <ws2tcpip.h>
	#include <errors.h>
	#include <conio.h>
    #define MAX_PATH_LENGTH MAX_PATH
    #define clearScreen clrscr
    #define pulisciSchermo (void) system("cls")
    #define removedir RemoveDirectory
#else
    #include <termios.h>
    #include <unistd.h>
    #define MAX_PATH_LENGTH 260
    #define AI_NUMERICHOST AI_NUMERICSERV
    #define DEFAULT_LINUX_SPOOLFOLDER "/var/spool"
    #define clearScreen write(1,"\E[H\E[2J",7); //  printf("\033[H\033[J")
    #define pulisciSchermo (void) system("clear")
    #define closesocket close
    #define removedir rmdir
#endif

#define EXIT_HASHMAPERROR 1

#define accessCurrentLmpMatrix accessCurrentSession()->MLSystem.lmpMatrix
#define accessCurrentMatrixList accessCurrentSession()->MLSystem.list
#define popCurrentMatrixList accessCurrentMatrixList->matrix
#define accessCurrentUser accessCurrentSession()->user





// #include <my_global.h>
#include <mysql.h>

#define DESCRIPTIONS_FOLDER DEFAULT_DESCRIPTIONS_FOLDER

#define getItemsListNo(mode) access(lists)[mode].itemsno
#define initList(mode) access(lists)[mode].items_list[PREV_LISTNODE] = access(lists)[mode].items_list[NEXT_LISTNODE] = (nodelist *)NULL

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


#define INIT_MATLIST NULL
#define INIT_LOGSLIST NULL
#define INIT_CURRENTMATRIX NULL
#define INIT_LASTMATRIXPRINTED NULL
#define INIT_LAYOUTLIST NULL
#define INIT_DIM 0

//
#define INIT_CM_ROWS INIT_DIM
#define INIT_CM_COLUMNS INIT_DIM
#define INIT_LMP_ROWS INIT_DIM
#define INIT_LMP_COLUMNS INIT_DIM
//


#define STARTING_DIM INIT_DIM+1
#define STARTING_CM_ROWS STARTING_DIM
#define STARTING_CM_COLUMNS STARTING_DIM
#define STARTING_LMP_ROWS STARTING_DIM
#define STARTING_LMP_COLUMNS STARTING_DIM

#define INIT_EXPRTYPE (exprType *)NULL
#define INIT_MATRIXOBJ (matrixObj *)NULL
#define INIT_LOGOBJ (logObj *)NULL
#define INIT_LAYOUTOBJ (layoutObj *)NULL

#define VECTORTYPE_COLUMN false
#define VECTORTYPE_ROW true

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

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#define MAX_BUFSIZ 2048
#define BUFLEN_FACTOR 4 // 5 // in order to increase Maximum Buffer Capacity (Length);
#define MIN_BUFLEN 1
#define MAX_BUFLEN (MAX_BUFSIZ * BUFLEN_FACTOR)
#define DEFAULT_BUFSIZE MAX_BUFSIZ

#define MINMINMIN_STRING 10
#define MINMIN_STRING 20
#define MIN_STRING 100
#define MAX_STRING 256
#define MAXX_STRING (MAX_STRING<<1)
#define MAXMAX_STRING (MAX_STRING<<2)
#define MMAX_STRING (5*MAX_STRING)
#define KSTRING 1000

#define SIGN_STRING MINMINMIN_STRING
//
#define INFO_STRING MIN_STRING
#define DINFO_STRING INFO_STRING<<1
#define PATHCONTAINS_STRING (MAX_PATH_LENGTH+(MIN_STRING<<1))
#define ERR_STRING PATHCONTAINS_STRING
#define MAX_FILE_LINES MAXX_STRING /// 255 /// INFO_STRING /// 100
#define MAX_IDENTIFIER_LENGTH MIN_STRING
#define MAX_EXTENSION_LENGTH SIGN_STRING

#define NULL_CHAR ""
#define BLANK_STRING " "
#define BLANK_CHAR ' '
#define UNKNOWN_STRING "(unknown)"
#define NONE_STRING "(none)"
#define UNDERSCORE_SEPERATOR_STRING "_"
#define FILENAMES_SEPERATOR UNDERSCORE_SEPERATOR_STRING
#define OPENING_BRACKET '{'
#define CLOSING_BRACKET '}'
#define OPENING_BRACKETSTRING "{"
#define CLOSING_BRACKETSTRING "}"
#define ESCAPE_CHARACTER '\\'

#define DBSERVER_LENGTH MIN_STRING
#define DBUSERNAME_LENGTH MIN_STRING
#define DBPASSWORD_LENGTH MIN_STRING
#define DATABASE_LENGTH MIN_STRING

#define INT64_TMAX ULLONG_MAX
#define NULL_VAL MAX_VAL

#define isNullVal(x) mpfr_cmp_d(x, NULL_VAL) == 0
#define exitHandleCheck (access(exitHandle) == EXITHANDLE_SETCMD || access(exitHandle) == EXITHANDLE_BACKCMD || !access(exitHandle))

#define MATRIXSET_COMMAND DEFAULT_MATRIXSET_COMMAND
#define MATRIXLET_COMMAND DEFAULT_MATRIXLET_COMMAND
#define MATRIXBACK_COMMAND DEFAULT_MATRIXBACK_COMMAND

#define TERMINATING_CHAR ';'
#define TERMINATING_STRING ";"
#define SCANFEXIT_CHAR '.'
#define SCANFEXIT_STRING "."
#define EXTENSION_DOT "."
#define EXTENSION_DOT_CHAR '.'
#define MATRIXES_SEPERATOR_STRING BLANK_STRING
#define MATRIXES_SEPERATOR_CHAR BLANK_CHAR

#define EQUAL_STRING "="
#define ASTERISK_STRING "*"
#define ASTERISK_CHAR '*'

#define SEARCHPROGRAM_FREESTRING BLANK_STRING
#define SEARCHPROGRAM_FREECHAR BLANK_CHAR

#define UPDINFO_UNSTABLESTRING BLANK_STRING
#define UPDINFO_UNSTABLECHAR BLANK_CHAR


#define MIN_EXTENSIVE_MULTITHREADING_CORESNO DEFAULT_MIN_EXTENSIVE_MULTITHREADING_CORESNO

// Bellman_Ford MACROS

#define BELLMANFORD_NEGATIVEWEIGHTCYCLUS false
#define BELLMANFORD_SUCCESS true

// Integer Bisection METHOD MACROS

#define INTBISECTION_NOSIGNCHANGES false
#define INTBISECTION_SUCCESS true

// Cholesky Decomposition MACROS

enum
{
	CHOLESKY_SUCCESS = 0,
	CHOLESKY_INVALIDMATRIXDIM,
	CHOLESKY_INVALIDSTORAGEDIM,
	CHOLESKY_INVALIDMATRIX
};



// Nelmin MACROS

enum
{
	NELMIN_INVALIDREQMIN = 0,
	NELMIN_INVALIDDIMENSION,
	NELMIN_INVALIDKONVGE,
	NELMIN_ALLOC_ERROR,
	NELMIN_SUCCESS,
	NELMIN_KCOUNTEXCEEDED
};

// Invmod MACROS

enum
{
	INVMOD_SUCCESS = 0,
	INVMOD_INVALIDELEMENTS,
	INVMOD_NONZEROELMSINMIXEDMODULUSPOSITIONS,
	INVMOD_NOTINVERTIBLEMATRIX,
};


// Ball 01 Monomial Integral
#define BALL01MONINT_INVALIDRETURNVALUE MAX_VAL

#define LEGENDRE_SET_INVALIDRETURNVALUE MAX_VAL

#define DLANOR_INVALIDRETURNVALUE MAX_VAL
#define DSTREM_INVALIDRETURNVALUE MAX_VAL

#define IMTQLX_ITERATIONS_ERROR false
#define IMTQLX_SUCCESS true

#define GAMAIN_INVALIDRETURNVALUE MAX_VAL
#define BETAIN_INVALIDRETURNVALUE MAX_VAL
#define BETANC_INVALIDRETURNVALUE MAX_VAL
#define NCBETA_INVALIDRETURNVALUE MAX_VAL
#define TNC_INVALIDRETURNVALUE MAX_VAL
#define GAMMAD_INVALIDRETURNVALUE MAX_VAL
#define GAMMDS_INVALIDRETURNVALUE MAX_VAL
#define ALNGAM_INVALIDRETURNVALUE MAX_VAL
#define R4NORMAL01CDFINVERSE_INVALIDRETURNVALUE MAX_VAL
#define R8NORMAL01CDFINVERSE_INVALIDRETURNVALUE MAX_VAL
#define PPCHI2_INVALIDRETURNVALUE MAX_VAL
#define XINBTA_INVALIDRETURNVALUE MAX_VAL
#define TRIGAMMA_INVALIDRETURNVALUE MAX_VAL
#define CHYPER_INVALIDRETURNVALUE MAX_VAL

#define RCONT_INVALIDPARAMS false
#define RCONT_SUCCESS true

// Simdo MACROS

enum
{
	SIMDO_SUCCESS = 0,
	SIMDO_JSUMEXCEEDEDIPROD,
	SIMDO_IVECILLEGALCOMPO
};

// Simplex Lattice Point MACROS

enum
{
	SIMPLEXLP_SUCCESS = 0,
	SIMPLEXLP_INVALIDPARAMS,
	SIMPLEXLP_ALTEREDDATA
};

// Temperature Conversion Functions MACROS
#define INVERSE_CONVERSION false
#define DIRECT_CONVERSION true

// Legendre Polynomial MACROS

#define LEGENDRE_DOMAIN_ERROR false
#define LEGENDRE_SUCCESS true

#define LEGENDRE_INVALIDRETURNVALUE MAX_VAL

// Bernstein Polynomial MACROS

#define BERNSTEIN_DOMAIN_ERROR false
#define BERNSTEIN_SUCCESS true

#define BERNSTEIN_INVALIDRETURNVALUE MAX_VAL

// Chebyshev Polynomial MACROS

#define CHEBYSHEV_DOMAIN_ERROR false
#define CHEBYSHEV_SUCCESS true

#define CHEBYSHEV_INVALIDRETURNVALUE MAX_VAL

// Hermite Polynomial MACROS

#define HERMITE_DOMAIN_ERROR false
#define HERMITE_SUCCESS true

#define HERMITE_INVALIDRETURNVALUE MAX_VAL

// Jacobi Polynomial MACROS

#define JACOBI_INVALIDRETURNVALUE MAX_VAL

// Laguerre Polynomial MACROS

#define LAGUERRE_DOMAIN_ERROR false
#define LAGUERRE_SUCCESS true

#define LAGUERRE_INVALIDRETURNVALUE MAX_VAL

// Bessel Polynomial MACROS

#define BESSEL_ACC 40.0
#define BESSEL_BIGNO 1.0e10
#define BESSEL_BIGNI 1.0e-10

#define BESSEL_INVALIDRETURNVALUE MAX_VAL


// FFT (Fast Fourier Transform); MACROS

#define BACKWARD_FFTRANSFORMATION false
#define FORWARD_FFTRANSFORMATION true

#define FFT_ALLOC_ERROR false
#define FFT_SUCCESS true

// Routh Criterion METHOD MACROS

#define ROUTHTABLE_EPSILON 0.0000000001
#define ROUTHTABLE_ALLOC_ERROR -1

// Jury Criterion METHOD MACROS

enum
{
	JURYTABLE_ALLOC_ERROR = 0,
	JURYTABLE_NOTSATISFIED,
	JURYTABLE_SATISFIED
};

//

#define USEFUL_EPSILON r8_epsilon()
#define USEFUL_TOLERANCE 100.00 * USEFUL_EPSILON

//

// EigenValues and EigenVectors MACROS

enum
{
	EIGVALUES_FOUNDEVS = 0,
	EIGVALUES_ALLOC_ERROR,
	EIGVALUES_INFEVS_ERROR
};

#define EIGENVALUES_PREC 0.00000000001

#define MAX_EIGVALUES_ITERATIONS DEFAULT_MAX_EIGVALUES_ITERATIONS

// SVD Macro
#define SVD_WITHU true
#define SVD_WITHV SVD_WITHU

#define SVD_WITHOUTU false
#define SVD_WITHOUTV SVD_WITHOUTU

#define SVD_TOLERANCE USEFUL_TOLERANCE
#define SVD_EPSILON USEFUL_EPSILON

// DSVD Max Iterations MACRO
#define MAX_DSVD_ITERATIONS DEFAULT_MAX_DSVD_ITERATIONS


#define OUTLIER_CONSTANT DEFAULT_OUTLIER_CONSTANT

#define MIN_OUTLIER_CONSTANT DEFAULT_MIN_OUTLIER_CONSTANT
#define MAX_OUTLIER_CONSTANT DEFAULT_MAX_OUTLIER_CONSTANT

#define MIDDLE_ELEMENT false
#define INITIAL_ELEMENT true

#define BY_START false
#define BY_GLOBALVAL true

#define SetDefaultColor() SetColor(COLOR_DEFAULT)

#define ROMAN_NUMBER_STRING MIN_STRING<<1
#define ROMAN_NUMBER_MAPSTRING 5
#define MIN_PROCESSABLE_ROMAN_NUMBER 1
#define MAX_PROCESSABLE_ROMAN_NUMBER 3000
#define MIN_PASCALTRIANGLE_ROWS 1
#define MAX_PASCALTRIANGLE_ROWS 50

#define EXIT_MESSAGE DEFAULT_EXIT_MESSAGE

#define ADDCONSTANTS
#define ADDPHYSICSCONSTANTS
#define ADDINFORMATICCONSTANTS
#define ADDTIMECONSTANTS
#define ADDFIDCONSTANTS


// MapperParser C Stuffs

#define itypCheckLen(a,b) a == _accessContainer(b).size
#define boolCheckLen(a,b) a == _accessContainer(b).size
#define dimtypCheckLen(a,b) a == _accessContainer(b).size
#define shortCheckLen(a,b) a == _accessContainer(b).size
#define intCheckLen(a,b) a == _accessContainer(b).size


#define vcheck(a,b) if(a){mpfr_set_si(val, true, MPFR_RNDN),b;}else{mpfr_set_si(val, false, MPFR_RNDN);}
#define check(a,b) mpfr_set_d(val, a ? b : MAX_VAL, MPFR_RNDN);

#define checkFID(a) mpfr_sgn(a) > 0 && mpfr_cmp_d(a, MAX_BASEFIDS) < 0
#define getFID(a) ext_math.base_functions[mpfr_get_ui(a, MPFR_RNDN)]

enum
{
	EAGER_EVALUATION = 0,
	LAZY_EVALUATION,
	INVALID_EVALUATION
};

// DEFAULT-CUSTOM TYPES TYPEDEFS


typedef long long int64_t;
typedef unsigned long long uint64_t;

//

typedef unsigned char sel_typ; // DEFINING SELECTION TYPE
typedef unsigned short fsel_typ; // DEFINING FORMATTABLE SELECTION TYPE
typedef sel_typ bool;


#define COMPARE_STRING_FACTOR MAX_STRING

//
#if (MAX_ROWS > 0 && MAX_COLUMNS > 0)
#if (MAX_ROWS < COMPARE_STRING_FACTOR && MAX_COLUMNS < COMPARE_STRING_FACTOR)
typedef fsel_typ dim_typ;
#else
typedef unsigned dim_typ;
#endif
#else
    #error Invalid definition of MAX_ROWS or MAX_COLUMNS.
#endif
//

typedef tipo_predefinito ityp;

#define EXPREVAL_TMPL "./ExprEval/exprtmpl.html" /// "http://expreval.sourceforge.net/exprtmpl.html"


// MPFR MACROS

#define mpfGZero(a) mpfr_sgn(a) > 0
#define mpfLZero(a) mpfr_sgn(a) < 0
#define mpfGeZero(a) mpfGZero(a) > 0 || mpfr_zero_p(a)
#define mpfLeZero(a) mpfLZero(a) < 0 || mpfr_zero_p(a)

enum
{
	CTRL_EVALTYPE = 0,
	CTRL_DIMCHECK,
	CTRL_OPRCHECK,
	CTRL_PATTERNCHECK,
};

typedef struct
{
	union 
	{
		struct
		{
			uint64_t (*ctrl_function)( mpfr_t [], sel_typ []);
			sel_typ * opr;
		} ReturnCTRLSystemPatterned;
		sel_typ opr;
		dim_typ dim;
		bool evalType;
	} ReturnCTRLSystem;
	sel_typ type;
} mpCTRL;

// 

// INCLUDING EXPREVAL STUFFS BY DIRECTLY
/* Define type of data to use */
// typedef double EXPRTYPE;
typedef mpfr_t EXPRTYPE;
typedef sel_typ EXPRERRTYPE;

/* Defines for various things */

/* Max id size */
#define EXPR_MAXIDENTSIZE 255

/* Error values */
enum
    {
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

    EXPR_ERROR_USER, /* Custom errors should be larger than this */
    
    
    EXPR_ERROR_UNKNOWN = 255, /* Unknown error */
	MAX_EXPR_ERROR = EXPR_ERROR_UNKNOWN
    };

/* Macros */

/* Forward declarations */
typedef struct _exprNode exprNode;
typedef map_t exprFuncList;
typedef struct _exprValList exprValList;
typedef struct _exprObj exprObj;

/* Function types */
// typedef int (*exprFuncType)(exprObj *obj, exprNode *nodes, int nodecount, EXPRTYPE **refs, int refcount, EXPRTYPE val);
typedef void * (*exprRealFuncType)(void *);
typedef void * (*exprFuncType)(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
typedef void (*exprReturnerType)( void * data, EXPRTYPE args[], exprNode *nodes, mpfr_t val);
typedef uint64_t (*exprReturnCtrlType)( EXPRTYPE args[], sel_typ opr[]);

typedef int (*exprBreakFuncType)(exprObj *obj);

/* Functions */

/* Version information function */
void exprGetVersion(int *major, int *minor);

/* Functions for function lists */
void exprFuncListCreate(exprFuncList *flist);
EXPRERRTYPE exprFuncListAdd(exprFuncList flist, char *name, exprReturnerType rt, exprFuncType ptr, exprRealFuncType rptr, mpCTRL * mp_ctrl, int MIN, int MAX, int refmin, int refmax);
void exprFuncListFree(exprFuncList *flist);
EXPRERRTYPE exprFuncListInit(exprFuncList flist);

/* Functions for value lists */
EXPRERRTYPE exprValListCreate(exprValList **vlist);
EXPRERRTYPE exprValListAdd(exprValList *vlist, char *name, EXPRTYPE val);
EXPRERRTYPE exprValListSet(exprValList *vlist, char *name, EXPRTYPE val);
EXPRERRTYPE exprValListGet(exprValList *vlist, char *name, EXPRTYPE val);
EXPRERRTYPE exprValListAddAddress(exprValList *vlist, char *name, EXPRTYPE *addr);
EXPRERRTYPE exprValListGetAddress(exprValList *vlist, char *name, EXPRTYPE **addr);
void *exprValListGetNext(exprValList *vlist, char **name, EXPRTYPE value, EXPRTYPE** addr, void *cookie);
EXPRERRTYPE exprValListFree(exprValList *vlist);
EXPRERRTYPE exprValListClear(exprValList *vlist);
int exprValListInit(exprValList *vlist);

/* Functions for expression objects */
EXPRERRTYPE exprCreate(exprObj **obj, exprFuncList flist, exprValList *vlist, exprValList *clist,
    exprBreakFuncType breaker, void *userdata);
EXPRERRTYPE exprFree(exprObj *obj);
EXPRERRTYPE exprClear(exprObj *obj);
EXPRERRTYPE exprParse(exprObj *obj, char *expr);
EXPRERRTYPE exprEval(exprObj *obj, EXPRTYPE val);
EXPRERRTYPE exprEvalNode(exprObj *obj, exprNode *nodes, int curnode, EXPRTYPE val, const bool alloc);
exprFuncList exprGetFuncList(exprObj *obj);
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

#include "ExprEval/exprpriv.h"

#define EXPREVAL_NALLOC false
#define EXPREVAL_ALLOC true

// Removed because It caused overhead in Stirling's function
// (There is another sqrt function which is assumed to be multiplied in the Formula
// #define STIRLING_CONSTANT sqrt(2*M_PI);

#define MAX_MEMOIZABLE_INDEX 100

#define MAX_FIBONACCI_MEMOIZABLE_INDEX DEFAULT_MAX_FIBONACCI_MEMOIZABLE_INDEX
#define MAX_FACTORIAL_MEMOIZABLE_INDEX DEFAULT_MAX_FACTORIAL_MEMOIZABLE_INDEX
#define MAX_EVEN_DOUBLEFACTORIAL_MEMOIZABLE_INDEX DEFAULT_MAX_EVEN_DOUBLEFACTORIAL_MEMOIZABLE_INDEX
#define MAX_ODD_DOUBLEFACTORIAL_MEMOIZABLE_INDEX DEFAULT_MAX_ODD_DOUBLEFACTORIAL_MEMOIZABLE_INDEX


#if (MAX_FIBONACCI_MEMOIZABLE_INDEX > MAX_MEMOIZABLE_INDEX || MAX_FACTORIAL_MEMOIZABLE_INDEX > MAX_MEMOIZABLE_INDEX || MAX_FPNUM_MEMOIZABLE_INDEX > MAX_MEMOIZABLE_INDEX)
    #error Invalid definition of one of MAX FUNCTIONS MEMOIZABLE INDICES.
#endif // MAX_FIBONACCI_MEMOIZABLE_INDEX

#define MIN_STIRLING_NUMBER DEFAULT_MIN_STIRLING_NUMBER

#define INPUT false
#define OUTPUT true

#define OVERFLOW_ERROR "\nINFO: An OVERFLOW ERROR might have occurred during (%hu,%hu) Elaborating Matrix Element...\n"

#define CALLOC_MODE false
#define MALLOC_MODE true

enum
{
	FIRST_MATRIX = 0,
	SECOND_MATRIX,
	THIRD_MATRIX,
	MATRIX_PRODUCT = THIRD_MATRIX,
	MATRIX_SUM = THIRD_MATRIX,
	MATRIX_SUB = THIRD_MATRIX,
	MAX_MATRICES
};

#define MIN_BLOCKSIZE 1
#define _BLOCK_SIZE DEFAULT_BLOCKSIZE
#define MIN_OSMM_DIM DEFAULT_MINOSMMDIM
#define MIN_STRASSEN_DIM DEFAULT_MINSTRASSENDIM

#define STARTING_MINOSMMDIM MIN_OSMM_DIM + 2
#define STARTING_MINSTRASSENDIM MIN_STRASSEN_DIM

#define NORMAL_MODE false
#define NOP_MODE true

#ifndef WINOS
    #define MIN(a,b) (((a)<(b))?(a):(b))
#endif

//
#define insertMatrix(x,y,z,k) enterMatrix(&x,&y,&z,k,true)
#define matrixFree(x,dim) _matrixFree(x,dim,NORMAL_MODE)
#define printMatrix(a,b,c) _printMatrix(a,b,c,NULL)

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
    COLOR_WHITEPLUS,
    COLOR_SMARTWHITE = 49,
    COLOR_LIGHTBLUE = 43
};

#ifdef WINOS
	#define _COLOR_ERROR 60 // COLOR_SMARTRED
	#define _COLOR_CREDITS 57 // COLOR_LIGHTBLUE // COLOR_CRYSTAL // COLOR_SMARTBLUE
	#define _COLOR_USER 58 // COLOR_LUMGREEN
	#define _COLOR_SYSTEM 60 // COLOR_SMARTYELLOW
	#define _COLOR_AUTHOR 54 // COLOR_SMARTPURPLE
	#define MAX_COLORS 64
#else
	#define _COLOR_ERROR 60 // COLOR_SMARTRED
    #define _COLOR_CREDITS 57 // COLOR_LIGHTBLUE // COLOR_CRYSTAL // COLOR_SMARTBLUE
    #define _COLOR_USER 58 // COLOR_LUMGREEN
    #define _COLOR_SYSTEM 59 // COLOR_SMARTYELLOW
    #define _COLOR_AUTHOR 60 // COLOR_SMARTPURPLE
    #define LAST_COLOR COLOR_WHITEPLUS
    #define MAX_COLORS LAST_COLOR+1
#endif

enum
{
    MEMBER_DEFAULTCOLOR = 0,
    MEMBER_COLORERROR,
    MEMBER_COLORCREDITS,
    MEMBER_COLORUSER,
    MEMBER_COLORSYSTEM,
    MEMBER_COLORAUTHOR,
    MAX_COLOR_TYPES
};

#define INIT_COLOR COLOR_WHITE
#define COLORS_PATH DEFAULT_COLORS_PATH

#define COLOR_DEFAULT access(colors)[MEMBER_DEFAULTCOLOR]
#define COLOR_ERROR access(colors)[MEMBER_COLORERROR]
#define COLOR_CREDITS access(colors)[MEMBER_COLORCREDITS]
#define COLOR_USER access(colors)[MEMBER_COLORUSER]
#define COLOR_SYSTEM access(colors)[MEMBER_COLORSYSTEM]
#define COLOR_AUTHOR access(colors)[MEMBER_COLORAUTHOR]

#define PRINTSPACE() msprintf(COLOR_DEFAULT, " ")
#define PRINTT() msprintf(COLOR_DEFAULT, "\t")
#define PRINTL() msprintf(43, "\n________________________________________________________________________________\n\n") // enhanced formatting

#define PRINTN() msprintf(COLOR_DEFAULT, "\n")
#define PRINT2N() msprintf(COLOR_DEFAULT, "\n\n")
#define PRINT5N() msprintf(COLOR_DEFAULT, "\n\n\n\n\n")
#define PRINT10N() msprintf(COLOR_DEFAULT, "\n\n\n\n\n\n\n\n\n\n")

#define SHOWPAUSEMESSAGE() PRINTL(), msprintf(COLOR_SYSTEM, "Press CTRL + C to stop Program Action.\n\n")
#define PRINTHOWTOBACKMESSAGE() msprintf(COLOR_SYSTEM, "Press ENTER to exit SubProgram.\n\n")


#define __pmode__ access(mode) // MACRO PROGRAM MODE

#define DISABLED false
#define ENABLED true

#define BY_NUMBERS false
#define BY_CHARS true

#define MPFR_INIT_ELEMENTS ENABLED
#define MPFR_SET_ELEMENTS DISABLED

// Cryptographic Stuffs

#define MAX_HEADERS 10
// #define HEXHASH_LEN 500
#define HEXHASH_LEN SHA512_HASHREPLENGTH
#define EXIT_MUTANTCODE 1
#define PERM_OWNER 2
#define MSS_COMPILER "gcc"
#define ECHO_CMD "echo"
#define ERROR_CMD "err"
#define PRINT_CMD "print"
#define SYSPRINT_CMD "sprint"
//#define DEVMODE_CMD "faa82036538b8f367fcf7bfd4c63b789"
//#define OWNERMODE_CMD "4ca1b49056a9cb5a1ed7af6aaed39915"
#define DEVMODE_CMD "71bc535a84a94e2c95a7a0d7de5b9c2b2feb6f688ebad343543b8efc79d76c6116fe0b0dce79ef69c15d1d374c030ee3f5a798518567baf5a3615782d82fd198"
#define OWNERMODE_CMD "592bf0c5506db50b5241500fb57e538e21c3165c991691069e33b07882ad9f69abeb30e57e3746c07bff58daf6178aa247f49feba651d980c987f4a92440f0a4"
#define DISABLEFEATURES_CMD "disable"
#define DISABLEFEATURESFAST_CMD "disable-fast"
#define DISABLEDEEP_CMD "disable-deep"
#define MUTABLE_CMD "mutable"
#define INLINE_CMD "inline"
#define MAKE_CMD "make"
#define REBUILD_CMD "rebuild"
#define PHONY_FILE ".phony"
#define MUTANT_DIR "mutant"
#define MUTANTSRC_DIR "mutantSrc"
#define MUTANTSRC_ARCHIVE "mutantSrc.zip"
#define MUTANTOBJ_ARCHIVE "mutantObj.zip"
#define MUTANT_REMOVEDIRSFASTCMD "cd "MUTANT_DIR" && rm -f *.c && rm -f *.h && rm -f mathSuite.ico && rm -f mathSuite_private.rc"
#define MUTANT_CREATING "\nCreating a mathSuite Mutant Version...\n"
#define MUTANT_WELCOME "This is a mutant mathSuite version"
#define _MSS_CMD DEFAULT_SCRIPTFILES_EXTENSION"_init"
#define MSS_CMD DEFAULT_SCRIPTFILES_EXTENSION"_exit"

#define MSS_TMPFILENAMSRC "mutant/tmp.c"

#ifdef WINOS
	#define MUTANT_FILENAME "mathSuiteMutant.exe"
	#define MAKE_EXECUTABLE "cd "MUTANT_DIR" && mingw32-make.exe"
	#define EXTRACTSRC_CMD "powershell Add-Type -A System.IO.Compression.FileSystem ; [IO.Compression.ZipFile]::ExtractToDirectory('"MUTANT_DIR"\\"MUTANTSRC_ARCHIVE"', '"MUTANT_DIR"')"
	#define _COMMIT_CMD "powershell Add-Type -A System.IO.Compression.FileSystem ; rm "MUTANT_DIR"\\"MUTANTOBJ_ARCHIVE" ; [IO.Compression.ZipFile]::CreateFromDirectory('"MUTANT_DIR"\\obj', '"MUTANT_DIR"\\"MUTANTOBJ_ARCHIVE"')"
	#define _CHECKOUT_CMD "powershell Add-Type -A System.IO.Compression.FileSystem ; rm -r "MUTANT_DIR"\\obj ; [IO.Compression.ZipFile]::ExtractToDirectory('"MUTANT_DIR"\\"MUTANTOBJ_ARCHIVE"', '"MUTANT_DIR"\\obj')"
	#define MUTANT_REMOVEDIRSCMD "rm -f *.c && rm -f *.h && rm -f mathSuite.ico && rm -f mathSuite_private.rc && powershell rm -r ExprEval ; rm -r jburkardt ; rm -r libxml"
	#define MSS_TMPFILENAMEXE "mutant\\tmp.exe"
#else
	#define MUTANT_FILENAME "mathSuiteMutant"
	#define MAKE_EXECUTABLE "cd "MUTANT_DIR" && make"
	#define EXTRACTSRC_CMD "cd "MUTANT_DIR" && unzip "MUTANTSRC_ARCHIVE" -d "MUTANTSRC_DIR
	#define _COMMIT_CMD "cd "MUTANT_DIR && zip -r "MUTANTOBJ_ARCHIVE"
	#define _CHECKOUT_cmd "cd "MUTANT_DIR" && rm -rd obj && unzip -r "MUTANTOBJ_ARCHIVE" -d obj"
	#define MUTANT_REMOVEDIRSCMD "rm -f *.c && rm -f *.h && rm -f mathSuite.ico && rm -f mathSuite_private.rc && rm -rd ExprEval && rm -rd jburkardt && rm -rd libxml"
	#define MSS_TMPFILENAMEXE "mutant/tmp"
#endif

#define INITIALIZING_RANDOM_SEED 0 // Invalid value for Random Seed.
// Used just to initializate the program random seed value.

#define DEFAULT_RANDOM_SEED (unsigned)time(NULL)
tipo_random_seed starting_random_seed;

#define MIN_RANDOMSEED 1
#define MAX_RANDOMSEED RAND_MAX

#define MIN_STABFACT 1
#define MAX_STABFACT DEFAULT_MAX_STABFACT // 20

enum
{
    ALGEBRA_REALNUMBERS,
    MIN_ALGEBRA = ALGEBRA_REALNUMBERS,
    ALGEBRA_COMPLEXNUMBERS,
    ALGEBRA_QUATERNIONS,
    ALGEBRA_OCTONIONS,
    ALGEBRA_SEDENIONS,
    MAX_ALGEBRA = ALGEBRA_SEDENIONS,
    _MAX_ALGEBRA
};

#define REAL_UNIT_NAME NULL_CHAR // "real"

#define _cabs(a,b) ___cabs(a, b, ALGEBRA_COMPLEXNUMBERS)
#define _qabs(a,b) ___cabs(a, b, ALGEBRA_QUATERNIONS)
#define _oabs(a,b) ___cabs(a, b, ALGEBRA_OCTONIONS)
#define _sabs(a,b) ___cabs(a, b, ALGEBRA_SEDENIONS)

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
    MAIN_BCALCULATOR = 0,
    MAIN_ADVANCEDCALCULATOR,
    MAIN_ALGEBRAOPERATIONS,
    #ifndef __DISABLE_MULTIUSER
    MAIN_MULTIUSER,
	#endif
    #ifndef __DISABLE_MSSMANAGER
        MAIN_MSSMANAGER,
    #endif // ALLOW_MSSMANAGER
    MAIN_CHANGESETTINGS,
    MAX_PROGRAMMI
};

enum
{
	#ifndef __DISABLE_MULTIUSER
	#ifndef __DISABLE_SERVER
  		MULTIUSER_SERVERMODE = 0,
  	#else
  		MULTIUSER_SERVERMODE = -1,
  	#endif
  	#ifndef __DISABLE_DATABASE
	  	MULTIUSER_DATABASE,
		MULTIUSER_LOGIN,
		MULTIUSER_REGISTER,
		MULTIUSER_EDITUSER,
		MULTIUSER_REMPERMIT,
		PROG_REMPERMIT,
	#endif
	#endif
	MAX_MULTIUSER_PROGS
};

// ADVANCED CALCULATOR MACRO AND ENUMERATIONS
//

enum
{
	#ifndef __DISABLE_MSSMANAGER
	    MSSMANAGER_COMMANDLINE = 0,
	    MSSMANAGER_EXECSCRIPTFILES,
	    MSSMANAGER_SHOWUSAGES,
    #endif
    MAX_MSSMANAGER_PROGS
};

#define MAX_ARGS 10

enum
{
	ADVCALC_CRYPTOGRAPHICHASH = 0,
    ADVCALC_SECONDGRADEEQUATIONSOLVER,
    ADVCALC_COMPLEXNUMBERSSUM,
    ADVCALC_COMPLEXNUMBERSPROD,
    ADVCALC_POLYNOMEVALUATOR,
    ADVCALC_POLYNOMDEVALUATOR,
    ADVCALC_GETFORMATTEDDATE,
    ADVCALC_ABSOLUTEORIENTATION,
    ADVCALC_ROUTHTABLE,
    ADVCALC_JURYTABLE,
    ADVCALC_NEWTONDIFFTABLES,
    ADVCALC_LAGRANGEINTERPOLATION,
    ADVCALC_FUNCTIONINTEGRATION,
    ADVCALC_STRAIGHTLINEFITTING,
    ADVCALC_PARABOLICCURVEFITTING,
    ADVCALC_LINEARSYSTEMSSOLVER,
    MAX_ADVCALC_PROGS
};


#define ABSOR_ERROR true
#define ABSOR_RDIM 0
#define ABSOR_DIM 1

// Cryptographic Programs

enum
{
	#ifndef __DISABLE_CRYPTOGRAPHICHASH
	CRYPTOGRAPHICHASH_MDC2,
	// CRYPTOGRAPHICHASH_MD2,
	CRYPTOGRAPHICHASH_MD4,
	CRYPTOGRAPHICHASH_MD5,
	CRYPTOGRAPHICHASH_RIPEMD160,
	CRYPTOGRAPHICHASH_WHIRLPOOL,
	CRYPTOGRAPHICHASH_SHA1,
	CRYPTOGRAPHICHASH_SHA224,
	CRYPTOGRAPHICHASH_SHA256,
	CRYPTOGRAPHICHASH_SHA384,
	CRYPTOGRAPHICHASH_SHA512,
	#endif
	MAX_CRYPTOGRAPHICHASH_PROGS
};

// SPECIAL MATRICES ENUMERATIONS
//

enum
{
	#ifndef __DISABLE_SPECIALMATRICES
		SPECMAT_PRINTLMPMATRIX,
		SPECMAT_PRINTMATRIXLIST,
		SPECMAT_PRINTCURRENTMLMATRIX,
		SPECMAT_SETCURRENTMLMATRIX,
		SPECMAT_GETCURRENTMATRIX,
		SPECMAT_DELMATRIXLIST,
		SPECMAT_LETSPECIALMATRIX,
	#endif
	MAX_SPECMAT_PROGS
};

// LINEAR ALGEBRA MACRO AND ENUMERATIONS
//
enum
{
	#ifndef __DISABLE_MATMANAGER
    	ALGOPS_MATRICESMANAGER = 0,
    #endif
    #ifndef __DISABLE_SPECIALMATRICES
    	ALGOPS_SPECIALMATRICES,
    #endif
    ALGOPS_MATRIXSORT,
    ALGOPS_MATRIXEIGVALUES,
    ALGOPS_NORMCALCULATOR,
    ALGOPS_DETERMINANTCALCULATOR,
    ALGOPS_TRACECALCULATOR,
    ALGOPS_RANKCALCULATOR,
    ALGOPS_MATRIXSVD,
    ALGOPS_INVERSEMATRIX,
    ALGOPS_FASTINVERSEMATRIX,
    ALGOPS_MATRIXCOFACTOR,
    ALGOPS_MATRIXADJOINT,
    ALGOPS_MATRIXTRANSPOSE,
    ALGOPS_MATRIXADD,
    ALGOPS_TENSORSSUM,
    ALGOPS_MATRIXMULTIPLICATION,
    ALGOPS_KRONPRODUCT,
    ALGOPS_MATRIXPOWER,
    ALGOPS_MATRIXKPOWER,
    ALGOPS_MATRIXPERVECTOR,
    ALGOPS_DOTPRODUCT,
    ALGOPS_PERSCALARMULTIPLICATION,
    ALGOPS_SCALARDIVISIONMATRIX,
    ALGOPS_ILLCONDITIONCHECKING,
    ALGOPS_MATRIXLUFACTORIZATION,
    MAX_ALGEBRA_OPERATIONS
};


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
	#ifndef __DISABLE_SYSTEM
		SETTINGS_SYSMANAGER,
	#endif 
    #ifndef __DISABLE_VARLISTMANAGER
        SETTINGS_ENVSMANAGER,
    #endif
    #ifndef __DISABLE_LOGSMANAGER
        SETTINGS_LOGSMANAGER,
    #endif
    #ifndef __DISABLE_SYSLOGMANAGER
        SETTINGS_SYSLOGMANAGER,
    #endif
    #ifndef __DISABLE_LAYOUTSMANAGER
    SETTINGS_LAYOUTSMANAGER,
    #endif
    #ifndef __DISABLE_COLSMANAGER
    SETTINGS_COLORSMANAGER,
    #endif
    #ifndef __DISABLE_LFSMANAGER
        SETTINGS_LFSMANAGER,
    #endif // ALLOW_LFSMANAGER
    #ifndef __DISABLE_PRECEDIT
        SETTINGS_CHANGEPRECISION,
    #endif
    #ifndef __DISABLE_STABFACTEDIT
        SETTINGS_CHANGESTABILIZERFACTOR,
    #endif
    #ifndef __DISABLE_BLOCKSIZEEDIT
    	SETTINGS_CHANGEBLOCKSIZE,
    #endif
    #ifndef __DISABLE_MINOSMMDIMEDIT
    	SETTINGS_CHANGEMINOSMMDIM,
    #endif
    #ifndef __DISABLE_MINSTRASSENDIMEDIT
    	SETTINGS_CHANGEMINSTRASSENDIM,
    #endif
    #ifndef __DISABLE_MINSRNUMBEREDIT
    	SETTINGS_CHANGEMINSRNUMBER,
    #endif
    #ifndef __DISABLE_ALGEBRAEDIT
        SETTINGS_CHANGEALGEBRA,
    #endif
    #ifndef __DISABLE_OUTLIERCONSTEDIT
    	SETTINGS_CHANGEOUTLIERCONST,
    #endif
    #ifndef __DISABLE_EXITCHAREDIT
        SETTINGS_CHANGEEXITCHAR,
    #endif
    #ifndef __DISABLE_RANDOMSEEDEDIT
        SETTINGS_CHANGERANDOMSEED,
    #endif
    #ifndef __DISABLE_BOOLVARSEDIT
        SETTINGS_CHANGEBOOLVALUES,
    #endif
    #ifndef __DISABLE_MMIEDIT
        SETTINGS_CHANGEMAXMEMOIZABLEINDICES,
        SETTINGS_EMPTYMEMOIZERSBUFFERS,
    #endif
    #ifndef __DISABLE_BUFFERSEDIT
        SETTINGS_EMPTYBUFFERS,
    #endif
    #ifndef __DISABLE_SETDEFAULTS
        SETTINGS_SETDEFAULTS,
    #endif
    MAX_SETTINGS
};


enum
{
    BOOLS_AUTOSETCURITEM = 0,
    BOOLS_ITEMSSELECTBYPATH,
    BOOLS_ITEMSAUTOSAVING,
    BOOLS_SAVERESULTS,
    BOOLS_SHOWVARLIST,
    BOOLS_SHOWDIFFTIME,
    BOOLS_RUNMUTCODEAFTERCOMP,
    BOOLS_EXITAFTERMUTCODEXEC,
    BOOLS_SHOWEXECTIME,
    BOOLS_PRINTROWSLABELS,
    BOOLS_DOMAINCHECK,
    BOOLS_LAZYEXECUTION,
    BOOLS_INVERSEOPERATIONS,
    BOOLS_DEGREESENTERING,
    BOOLS_PROGREPEATCHECK,
    BOOLS_ARRAYSAUTORESET,
    BOOLS_ADJUSTARRAYSINDEX,
    BITMASK_BITFIELD = BOOLS_ADJUSTARRAYSINDEX,
    MAX_BOOL_SETTINGS
};

#define BITMASK_AUTOSETCURITEM 1
#define BITMASK_ITEMSSELECTBYPATH 2
#define BITMASK_ITEMSAUTOSAVING 4
#define BITMASK_SAVERESULTS 8
#define BITMASK_SHOWVARLIST 16
#define BITMASK_SHOWDIFFTIME 32
#define BITMASK_RUNMUTCODEAFTERCOMP 64
#define BITMASK_EXITAFTERMUTCODEXEC 128
#define BITMASK_SHOWEXECTIME 256
#define BITMASK_PRINTROWSLABELS 512
#define BITMASK_DOMAINCHECK 1024
#define BITMASK_LAZYEXECUTION 2048
#define BITMASK_INVERSEOPERATIONS 4096
#define BITMASK_DEGREESENTERING 8192
#define BITMASK_PROGREPEATCHECK 16384
#define BITMASK_ARRAYSAUTORESET 32768
#define BITMASK_ADJUSTARRAYSINDEX 65536

#define MAX_DAYSWEEK 7
#define MAX_MONTH_DAYS 31
#define MAX_MONTHS 12
#define GETDAYNAME_MAP_LENGTH 12
#define DAYSWEEK_LENGTH SIGN_STRING
#define MONTHS_LENGTH SIGN_STRING

#define isSett(x) ((access(curLayout)->bools & suite_c.bools[x].bmask) == suite_c.bools[x].bmask)
#define isnSett(x) (!(access(curLayout)->bools & suite_c.bools[x].bmask))

// ENVS MANAGER STUFFS

#define parse(a,b) _parse(a,b,NULL)

enum 
{
	#ifndef __DISABLE_SYSTEM
	SYSTEM_SCRIPTINGMODE,
	SYSTEM_EXIT,
	SYSTEM_MUTANT,
	SYSTEM_COMMIT,
	SYSTEM_CHECKOUT,
	SYSTEM_CLEAN,
	SYSTEM_PURGE,
	SYSTEM_PURGEDEEP,
	SYSTEM_PURGEFAST,
	#endif
	MAX_SYSTEM_PROGS
};

enum
{
	#ifndef __DISABLE_VARLISTMANAGER
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
    ENVS_RENAME,
    #ifndef __DISABLE_DATABASE
	    ENVS_PERSIST,
	    ENVS_RETRIEVE,
    #endif
    #endif
    MAX_ENVSMANAGER_PROGS
};

enum
{
	#ifndef __DISABLE_MATMANAGER
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
    MATRICES_RENAME,        //!U,
    #ifndef __DISABLE_DATABASE
	    MATRICES_PERSIST,
	    MATRICES_RETRIEVE,
    #endif
    #endif
    MAX_MATMANAGER_PROGS
};

enum
{
	#ifndef __DISABLE_LOGSMANAGER
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
    LOGS_RENAME,         //!U
    #endif
    MAX_LOGSMANAGER_PROGS
};

enum
{
	#ifndef __DISABLE_SYSLOGMANAGER
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
    SYSLOG_RENAME,          //!U
    #endif
    MAX_SYSLOGMANAGER_PROGS
};
enum
{
	#ifndef __DISABLE_LAYOUTSMANAGER
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
    LAYOUTS_RENAME,
    #endif
    MAX_LAYOUTSMANAGER_PROGS
};

enum
{
	#ifndef __DISABLE_COLSMANAGER
    COLORS_CHANGE = 0,
    COLORS_FILESLOADER,
    COLORS_BACKUPFILES,
    #endif
    MAX_COLSMANAGER_PROGS
};

enum
{
	#ifndef __DISABLE_LFSMANAGER
    LFS_LOADER = 0,
    LFS_CREATE,
    #endif
    MAX_LFSMANAGER_PROGS
};

#define MAX_CASEINSENSITIVE_CHARS_ALPHABET 26
#define _MAX_OMNIPRESENT_ELEMENTS 3
#define MAX_OMNIPRESENT_ELEMENTS _MAX_OMNIPRESENT_ELEMENTS-1

// Macro per alcune variabili booleane del programma
#define dcheck isSett(BOOLS_DOMAINCHECK)
#define lazy_exec isSett(BOOLS_LAZYEXECUTION)


#define PROGRAM_INFORMATIONS MAX_PROGRAMMI
#define EXIT_FROM_PROGRAM MAX_PROGRAMMI+1


// VALORE DI INIZIALIZZAZIONE CARATTERE DI USCITA DAL PROGRAMMA
#define INITIALIZING_EXIT_CHAR '£'
#define INITIALIZING_DEFAULT_COLOR COLOR_WHITE

#define MAIN_COLOR DEFAULT_COLOR

#define COMPAREDOUBLE_PRECISION 0.001

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
    FID_RAD,
    MAX_FIDS
};

enum
{
    BASEFID_SIN = 0,
    BASEFID_COS,
    BASEFID_SINH,
    BASEFID_COSH,
    BASEFID_SEC,
    BASEFID_SECH,
    BASEFID_CSC,
    BASEFID_CSCH,
    BASEFID_ASIN,
    BASEFID_ASINH,
    BASEFID_ACOS,
    BASEFID_ACOSH,
    BASEFID_ASEC,
    BASEFID_ASECH,
    BASEFID_ACSC,
    BASEFID_ACSCH,
    BASEFID_TAN,
    BASEFID_TANH,
    BASEFID_ATAN,
    BASEFID_ATANH,
    BASEFID_COT,
    BASEFID_COTH,
    BASEFID_ACOT,
    BASEFID_ACOTH,
    BASEFID_HSIN,
    BASEFID_HSINH,
    BASEFID_QSIN,
    BASEFID_QSINH,
    BASEFID_HCOS,
    BASEFID_HCOSH,
    BASEFID_QCOS,
    BASEFID_QCOSH,
    BASEFID_HSEC,
    BASEFID_HSECH,
    BASEFID_QSEC,
    BASEFID_QSECH,
    BASEFID_HCSC,
    BASEFID_HCSCH,
    BASEFID_QCSC,
    BASEFID_QCSCH,
    BASEFID_HTAN,
    BASEFID_HTANH,
    BASEFID_QTAN,
    BASEFID_QTANH,
    BASEFID_HCOT,
    BASEFID_HCOTH,
    BASEFID_QCOT,
    BASEFID_QCOTH,
    BASEFID_VSIN,
    BASEFID_VSINH,
    BASEFID_CVSIN,
    BASEFID_CVSINH,
    BASEFID_VCOS,
    BASEFID_VCOSH,
    BASEFID_CVCOS,
    BASEFID_CVCOSH,
    BASEFID_HVSIN,
    BASEFID_HVSINH,
    BASEFID_HCVSIN,
    BASEFID_HCVSINH,
    BASEFID_QVSIN,
    BASEFID_QVSINH,
    BASEFID_QCVSIN,
    BASEFID_QCVSINH,
    BASEFID_HVCOS,
    BASEFID_HVCOSH,
    BASEFID_HCVCOS,
    BASEFID_HCVCOSH,
    BASEFID_QVCOS,
    BASEFID_QVCOSH,
    BASEFID_QCVCOS,
    BASEFID_QCVCOSH,
    BASEFID_ESEC,
    BASEFID_ESECH,
    BASEFID_ECSC,
    BASEFID_ECSCH,
    BASEFID_HESEC,
    BASEFID_HESECH,
    BASEFID_HECSC,
    BASEFID_HECSCH,
    BASEFID_QESEC,
    BASEFID_QESECH,
    BASEFID_QECSC,
    BASEFID_QECSCH,
    BASEFID_SINC,
    BASEFID_SINCH,
    BASEFID_HSINC,
    BASEFID_HSINCH,
    BASEFID_QSINC,
    BASEFID_QSINCH,
    BASEFID_COSC,
    BASEFID_COSCH,
    BASEFID_HCOSC,
    BASEFID_HCOSCH,
    BASEFID_QCOSC,
    BASEFID_QCOSCH,
    BASEFID_SECC,
    BASEFID_SECCH,
    BASEFID_HSECC,
    BASEFID_HSECCH,
    BASEFID_QSECC,
    BASEFID_QSECCH,
    BASEFID_CSCC,
    BASEFID_CSCCH,
    BASEFID_HCSCC,
    BASEFID_HCSCCH,
    BASEFID_QCSCC,
    BASEFID_QCSCCH,
    BASEFID_TANC,
    BASEFID_TANCH,
    BASEFID_HTANC,
    BASEFID_HTANCH,
    BASEFID_QTANC,
    BASEFID_QTANCH,
    BASEFID_COTC,
    BASEFID_COTCH,
    BASEFID_HCOTC,
    BASEFID_HCOTCH,
    BASEFID_QCOTC,
    BASEFID_QCOTCH,
    BASEFID_LOG,
    BASEFID_LOG10,
    BASEFID_LOG2,
    BASEFID_LOGC,
    BASEFID_LOG10C,
    BASEFID_LOG2C,
    BASEFID_LOG1P,
    BASEFID_LOG1PC,
    BASEFID_EXP,
    BASEFID_EXPC,
    BASEFID_EXP10,
    BASEFID_EXP10C,
    BASEFID_EXP2,
    BASEFID_EXP2C,
    BASEFID_FLOOR,
    BASEFID_CEIL,
    BASEFID_DEG,
    BASEFID_RAD,
    MAX_BASEFIDS
};


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

// basicCalculator() GENERIC STUFFS...
#define BCALC_BASECHANGE_MINBASE 2
#define BCALC_BASECHANGE_MAXBASE 10

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

// BASIC CHECK FOR DMA (Dynamic Memory Allocation);
#define errMem(pntr, rval) if(checkErrMem(pntr)) return rval
#define errMemEx(pntr, rval) if(checkErrMem(pntr)){ result = rval; return &result; }

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
    ROWS = 0,
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


// #define DIM_BITFIELD (sizeof(dim_typ)-1)*<<3
#define DIM_BITFIELD 8
#define DATE_BITFIELD 5
#define DOMAIN_BITFIELD 4
#define BOOL_BITFIELD 1

#define MIN_NEWTON_DIFFTABLES_DIM 1
#define MAX_NEWTON_DIFFTABLES_DIM 15


// USED TO SHOW BACKWARD && FORWARD
// NEWTON's DIFFERENCE TABLES
#define BACKWARD_DIFFTAB false
#define FORWARD_DIFFTAB true

enum
{
    EXITHANDLE_EXIT = 0,
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

#define RCF_NONE 0
#define RCF_EXPREVAL 1
#define RCF_DEEP 2
#define RCF_SKIP 4

#define PTI_NONE 0
#define PTI_SAVECURRENT 1
#define PTI_SHOWITEM 2
#define PTI_UPDINFO 4

#define EQUALMATRIX_NOREALLOC false
#define EQUALMATRIX_REALLOC true

enum
{
	SPECIALMATRIX_CURRENTMATRIX = 0,
	SPECIALMATRIX_LMPMATRIX,
	SPECIALMATRIX_MATRIXLIST,
	MAX_SPECIALMATRICES
};

#define LEFT_OPR 0
#define RIGHT_OPR 1

#define MAX_MATRIXLISTITEMS 255
#define MAX_RIGHE_PER_COLONNE (int)(floor(sqrt(access(curLayout)->matrix_max_rows*access(curLayout)->matrix_max_columns))) // DIMENSIONE MASSIMA MATRICI QUADRATE nXn

#define SUMMATION_SUM false
#define SUMMATION_SUB true

#define r8vec_sum(a,b) summation(a,SUMMATION_SUM,b)

#define PRODUCTORY_MUL false
#define PRODUCTORY_DIV true

#define r8vec_product(a,b) productory(a,PRODUCTORY_MUL,b) 


// used for CURVEs FITTING or INTERPOLATIONS SubPrograms...
//
#define XROW 0
#define YROW 1
#define ZROW 2 // useful for Hermite k tensors' indices
//

#define F_OBSERVATIONS XROW
#define F_VARIABLES YROW
#define F_MAXCLUSTERS ZROW

#define F_SDIMENSIONS XROW
#define F_NPOINTS YROW

#define FIRST_NUMBER XROW
#define SECOND_NUMBER YROW

#define FIRST_VECTOR XROW
#define SECOND_VECTOR YROW

#define FIRST_QUARTILE XROW
#define THIRD_QUARTILE YROW

#define FIRST_QUARTILE_CONSTANT 0.25
#define SECOND_QUARTILE_CONSTANT 0.50
#define THIRD_QUARTILE_CONSTANT 0.75

// COMPLEX Numbers MACROS...
//
#define REAL_PART XROW
#define IMAG_PART YROW
//
#define MAX_COMPLEX_UNITS 2


// IperCOMPLEX Numbers MACROS
//
enum
{
    QUATERNIONS_REALPART = 0,
    QUATERNIONS_IPART,
    QUATERNIONS_JPART,
    QUATERNIONS_KPART,
    MAX_QUATERNIONS_UNITS
};

enum
{
    OCTONIONS_REALPART = 0,
    OCTONIONS_E1PART,
    OCTONIONS_E2PART,
    OCTONIONS_E3PART,
    OCTONIONS_E4PART,
    OCTONIONS_E5PART,
    OCTONIONS_E6PART,
    OCTONIONS_E7PART,
    MAX_OCTONIONS_UNITS
};

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
    SEDENIONS_E15PART,
    MAX_SEDENIONS_UNITS
};

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
	//FUNCTION_FACTORIAL,
	FUNCTION_EVEN_DOUBLEFACTORIAL,
	FUNCTION_ODD_DOUBLEFACTORIAL
};

#define LAST_MEMOIZABLE_FUNCTION FUNCTION_ODD_DOUBLEFACTORIAL
#define MAX_MEMOIZABLE_FUNCTIONS LAST_MEMOIZABLE_FUNCTION+1

#define _MAX_MEMOIZABLE_FUNCTIONS MAX_MEMOIZABLE_FUNCTIONS+1

#define MAX_STD_BUFFERS 3
#define _MAX_STD_BUFFERS MAX_STD_BUFFERS+1

#define INVALIDRETURNVALUE_NPRIMENUMBER 0.00
#define INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS -1.00

#define INVALIDRETURNVALUE_FIBONACCI INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS
#define INVALIDRETURNVALUE_FACTORIAL INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS
#define INVALIDRETURNVALUE_DOUBLEFACTORIAL INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS

#define INVALIDRETURNVALUE_EVEN_DOUBLEFACTORIAL INVALIDRETURNVALUE_DOUBLEFACTORIAL
#define INVALIDRETURNVALUE_ODD_DOUBLEFACTORIAL INVALIDRETURNVALUE_DOUBLEFACTORIAL

#define INVALIDRETURNVALUE_MEMFUNCS INVALIDRETURNVALUE_MEMOIZABLEFUNCTIONS

// PULIZIA DEL BUFFER
#define CLEARBUFFER() fflush(stdin) // cleanBuffer()

#define TRIGONOMETRIC_DOMAIN(x) (mpfr_cmp_si(x,-1) >= 0 && mpfr_cmp_si(x,1) <= 0)

#define random(x) ((rand()%x)+1)
#define randomize() srand(starting_random_seed=DEFAULT_RANDOM_SEED)
#define random2(x,y) randomize(),random(y)
#define enhanced_random(x,y) srand((unsigned)x),random(y)


// PROJECT TYPEDEFS

#define MAX_ALIAS DEFAULT_MAX_ALIAS

typedef const dim_typ dim_typ2[2];
typedef const dim_typ dim_typ3[3];

/// CRYPTOGRAPHIC STUFFS

__MATHSUITE void hexhash(const register sel_typ, const char *, char *);

#define MDC2_HASHREPLENGTH 32
// #define MD2_HASHREPLENGTH 32
#define MD4_HASHREPLENGTH 32
#define MD5_HASHREPLENGTH MD4_HASHREPLENGTH
#define RIPEMD160_HASHREPLENGTH 40
#define WHIRLPOOL_HASHREPLENGTH 128
#define SHA1_HASHREPLENGTH 40
#define SHA224_HASHREPLENGTH 56
#define SHA256_HASHREPLENGTH 64
#define SHA384_HASHREPLENGTH 96
#define SHA512_HASHREPLENGTH 128

enum
{
	HEXHASH_MDC2,
	// HEXHASH_MD2,
	HEXHASH_MD4,
	HEXHASH_MD5,
	HEXHASH_RIPEMD160,
	HEXHASH_WHIRLPOOL,
	HEXHASH_SHA1,
	HEXHASH_SHA224,
	HEXHASH_SHA256,
	HEXHASH_SHA384,
	HEXHASH_SHA512,
	MAX_HEXHASHES
};

#define hexhashMDC2(a,b) hexhash(HEXHASH_MDC2, a, b)
// #define hexhashMD2(a,b) hexhash(HEXHASH_MD2, a, b)
#define hexhashMD4(a,b) hexhash(HEXHASH_MD4, a, b)
#define hexhashMD5(a,b) hexhash(HEXHASH_MD5, a, b)
#define hexhashRIPEMD160(a,b) hexhash(HEXHASH_RIPEMD160, a, b)
#define hexhashWHIRLPOOL(a,b) hexhash(HEXHASH_WHIRLPOOL, a, b)
#define hexhashSHA1(a,b) hexhash(HEXHASH_SHA1, a, b)
#define hexhashSHA224(a,b) hexhash(HEXHASH_SHA224, a, b)
#define hexhashSHA256(a,b) hexhash(HEXHASH_SHA256, a, b)
#define hexhashSHA384(a,b) hexhash(HEXHASH_SHA384, a, b)
#define hexhashSHA512(a,b) hexhash(HEXHASH_SHA512, a, b)

// END CRYPTOGRAPHIC STUFFS

enum
{
	TYPID_FLOAT,
	TYPID_BOOL,
	TYPID_USHRT,
	TYPID_SHRT,
	TYPID_INT
};

typedef struct
{
	char cmdname[SIGN_STRING];
    const char name[DINFO_STRING];
    const char usage[INFO_STRING];
    #ifndef __DISABLE_DATABASE
    	const sel_typ level;
    #endif
    void (*program_function)(const sel_typ, char **);
    const dim_typ argc: DIM_BITFIELD;
    const bool automatic: BOOL_BITFIELD;
    const bool isFather: BOOL_BITFIELD;
} sprog;

#define NEXT_LISTNODE false
#define PREV_LISTNODE true

typedef struct _nodelist
{
    char path[MAX_PATH_LENGTH];
    void *data; // generic pointer to "data" data.
    struct _nodelist * (ref[2]);
} nodelist;

typedef struct
{
    mpfr_t *matrix;
    dim_typ dim[2];
} matrixObj;

typedef struct _mNodelist
{
	matrixObj * matrix;
	struct _mNodelist * next;
} mNodelist; 

typedef struct
{
    char *buffer;
    size_t buflen;
} logObj;

typedef struct
{
    mpfr_t * memoizer;
    dim_typ current_max_index: DIM_BITFIELD;
} memObj;

#ifndef __DISABLE_DATABASE
typedef struct
{
	int iduser;
	char name[USER_MAXNAMELENGTH];
	char surname[USER_MAXSURNAMELENGTH];
	char username[USER_MAXUSERNAMELENGTH];
	char email[USER_MAXEMAILLENGTH];
	char password[USER_MAXPASSWORDLENGTH];
	sel_typ level;
} userObj;
#endif

typedef struct
{
	#ifndef __DISABLE_MULTIUSER
	#ifndef __DISABLE_SERVER
	struct
	{
		char server[SERVER_LENGTH];
		int port;
		sel_typ family;
		bool resolution;
		dim_typ nListeningSockets: DIM_BITFIELD;
		dim_typ nAcceptableConnections: DIM_BITFIELD;
		int bufferLength;
	} server;
	#endif
	#ifndef __DISABLE_DATABASE
	struct
	{
		MYSQL * con;
		char server[DBSERVER_LENGTH];
		char username[DBUSERNAME_LENGTH];
		char password[DBPASSWORD_LENGTH];
		char database[DATABASE_LENGTH];
		dim_typ maxUsers: DIM_BITFIELD;
		fsel_typ defaultLevel: DIM_BITFIELD;
		bool identifierMode: BOOL_BITFIELD;
	} database;
	#endif
	#endif
	char exit_char;
    fsel_typ precision: DIM_BITFIELD;
    fsel_typ stabilizer_factor: DIM_BITFIELD;
    fsel_typ algebra: DIM_BITFIELD;
    float outlier_constant;
    dim_typ max_matrix_list_items: DIM_BITFIELD;
    dim_typ matrix_max_rows;
    dim_typ matrix_max_columns;
    dim_typ block_size;
    dim_typ min_osmm_dim;
    dim_typ min_strassen_dim;
    dim_typ max_eigvalues_iterations;
    dim_typ max_dsvd_iterations;
    dim_typ max_memoizable_indices[MAX_MEMOIZABLE_FUNCTIONS];
    dim_typ min_stirling_number: DIM_BITFIELD;
    dim_typ basecalc_minbase: DIM_BITFIELD;
    dim_typ basecalc_maxbase: DIM_BITFIELD;
    int max_changebase_binary_convnum;
    dim_typ min_newton_difftables_dim: DIM_BITFIELD;
    dim_typ max_newton_difftables_dim: DIM_BITFIELD;
    dim_typ min_roman_number: DIM_BITFIELD;
    dim_typ max_roman_number;
    dim_typ pascal_triangle_min_rows: DIM_BITFIELD;
    dim_typ pascal_triangle_max_rows: DIM_BITFIELD;
    unsigned bools: BITMASK_BITFIELD;
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
    LAYOUTS,
    LAST_PROGRAM_ITEMTYPE = LAYOUTS,
    MAX_LISTS
};

// ext_math.C STRUCTS DEFINITIONS


struct ext_math_type
{
    const char funcnames[MAX_FIDS][MIN_STRING];
    void (* const functions[MAX_FIDS])(mpfr_t, register mpfr_t);
    ityp (* const base_functions[MAX_BASEFIDS])(const register ityp);
    const char days_week_names[MAX_DAYSWEEK][DAYSWEEK_LENGTH];
    const char months_names[MAX_MONTHS][MONTHS_LENGTH];
    const sel_typ getdaynum_map[GETDAYNAME_MAP_LENGTH];
    const char romn_thousand_map[3][ROMAN_NUMBER_MAPSTRING];
    const char romn_map[3][9][ROMAN_NUMBER_MAPSTRING];
};
//

#define INVALID_ARGC 255
#define ONEARG_ARGC 254

enum
{
	LEVEL_GUEST = 0,
	LEVEL_USER,
	LEVEL_ADMIN,
	LEVEL_DEVELOPER,
	LEVEL_OWNER,
	MAX_USERTYPES
};

struct prog_constants
{
    const char listsnames[MAX_LISTS][SIGN_STRING];
	const char usertypes_names[MAX_USERTYPES][MINMIN_STRING];
	const char identifier_names[2][SIGN_STRING];
    const char bools_names[MAX_BOOL_SETTINGS][MIN_STRING];
    const char exprErrorString[EXPR_ERROR_USER+1][MIN_STRING];
    struct
    {
        const unsigned bmask;
        const bool default_val;
    } const bools[MAX_BOOL_SETTINGS];
    const dim_typ max_memoizable_indices[MAX_MEMOIZABLE_FUNCTIONS];
	const char memoizers_names[_MAX_MEMOIZABLE_FUNCTIONS][INFO_STRING];
    const char algebra_elements_names[_MAX_ALGEBRA][INFO_STRING];
    const char algebra_imaginary_units_names[_MAX_ALGEBRA][MAX_SEDENIONS_UNITS][SIGN_STRING];
    const char colors_types_names[MAX_COLOR_TYPES][INFO_STRING];
    const char colors_names[MAX_COLORS][INFO_STRING];
};

#define MAX_PRMTYPES 5
#define MAX_PRMCONTAINERSIZE 5000
#define accessContainer(a) access(PRMSystem).container[(uint64_t)a]
#define _accessContainer(a) accessContainer(mpfr_get_ui(a, MPFR_RNDN))

#define MAX_CONTAINERSBUFFERLEN 30
#define MAX_FLOATBUFFER MAX_CONTAINERSBUFFERLEN
#define MAX_BOOLBUFFER MAX_CONTAINERSBUFFERLEN
#define MAX_USHRTBUFFER MAX_CONTAINERSBUFFERLEN
#define MAX_SHRTBUFFER MAX_CONTAINERSBUFFERLEN
#define MAX_INTBUFFER MAX_CONTAINERSBUFFERLEN

#define MAX_USERS 255
#define MAX_SESSIONS MAX_USERS

#define MASTERSESSION_NAME "MASTER"
#define MASTERSESSION_SURNAME "SESSION"
#define MASTERSESSION_USERNAME "root"

#define DEFAULT_LEVEL LEVEL_USER

#define WRITEUSER_INSERT true
#define WRITEUSER_UPDATE false

#define USERNAME_ADMIN "root"
#define USERNAME_DEVELOPER "developer"
#define USERNAME_OWNER "owner"

#define INVALID_UID -1
#define MASTER_SESSION 0
#define CLIENTSERVER_BUFSIZE SERVER_BUFSIZE*10

#ifdef __DISABLE_MULTIUSER
	#ifndef __DISABLE_SERVER
		#define __DISABLE_SERVER
	#endif
	#ifndef __DISABLE_DATABASE
		#define __DISABLE_DATABASE
	#endif
#else
#if(defined(__DISABLE_SERVER) && defined(__DISABLE_DATABASE))
	#define __DISABLE_MULTIUSER
#endif
#endif

#include "ExprEval/jburkardt/structs.h"


typedef struct
{
	struct
	{
		mNodelist * list;
		matrixObj * lmpMatrix;
		dim_typ cur_item: DIM_BITFIELD;
    	dim_typ itemsno: DIM_BITFIELD;
	} MLSystem;
	dim_typ lastItem[2];
	userObj * user;
} session; 

struct program
{
	struct
	{
		struct
		{
			void *pnt;
			uint64_t size;
			bool _volatile;
		} container[MAX_PRMCONTAINERSIZE*MAX_PRMTYPES];
		uint64_t currentIndex;
		bool vectorType;
	} PRMSystem;
	#ifndef __DISABLE_DATABASE
		map_t sessions;
	#endif
	struct
    {
        nodelist * (items_list[2]);
        dim_typ cur_item: DIM_BITFIELD;
        dim_typ itemsno: DIM_BITFIELD;
    } lists[MAX_LISTS];
	exprType *exprVars;
    matrixObj *curMatrix;
    // matrixObj *lmpMatrix;
    logObj *curLog;
    logObj *sysLog;
    layoutObj *curLayout;
    memObj sysMem[MAX_MEMOIZABLE_FUNCTIONS];
    exprFuncList func_list;
    exprValList *const_list;
    char sysLogPath[MAX_PATH_LENGTH];
    char colors_path[MAX_PATH_LENGTH];
    #ifndef __DISABLE_SERVER
		int lastSessionSocket;
	#endif
    fsel_typ random_seed;
    sel_typ colors[MAX_COLOR_TYPES];
    sel_typ mode;
    volatile sel_typ exitHandle;
    volatile bool mss: BOOL_BITFIELD;
    volatile bool sigresult: BOOL_BITFIELD;
    #ifndef __DISABLE_SERVER
    	volatile bool server_mode: BOOL_BITFIELD;
    	// char serverBuffer[CLIENTSERVER_BUFSIZE];
    	int serverBufferLength;
    	char *serverBuffer;
    #endif
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

#ifdef WINOS
	_CRTIMP extern int errno;
#else
	extern int errno;
#endif

/// extern struct program suite;

extern sprog main_menu[MAX_PROGRAMMI],
                    adv_calc[MAX_ADVCALC_PROGS],
                    alg_operations[MAX_ALGEBRA_OPERATIONS];
                    
#ifndef __DISABLE_SPECIALMATRICES
	extern sprog specmat_prog[MAX_SPECMAT_PROGS];
#endif
                    
#ifndef __DISABLE_MULTIUSER
	extern sprog multiuser_prog[MAX_MULTIUSER_PROGS];
#endif

#ifndef __DISABLE_CRYPTOGRAPHICHASH
	extern sprog cryptographic_hash_prog[MAX_CRYPTOGRAPHICHASH_PROGS];
#endif

#ifndef __DISABLE_MSSMANAGER
    extern sprog mss_manager[MAX_MSSMANAGER_PROGS];
#endif // ALLOW_MSSMANAGER

#ifndef FREEZE_SETTINGSMANAGER
    extern sprog change_settings[MAX_SETTINGS];
#endif

#ifndef __DISABLE_SYSTEM
	extern sprog system_prog[MAX_SYSTEM_PROGS];
#endif

#ifndef __DISABLE_VARLISTMANAGER
    extern sprog envs_manager[MAX_ENVSMANAGER_PROGS];
#endif

#ifndef __DISABLE_MATMANAGER
    extern sprog mat_manager[MAX_MATMANAGER_PROGS];
#endif

#ifndef __DISABLE_LOGSMANAGER
    extern sprog logs_manager[MAX_LOGSMANAGER_PROGS];
#endif

#ifndef __DISABLE_SYSLOGMANAGER
    extern sprog syslog_manager[MAX_SYSLOGMANAGER_PROGS];
#endif

#ifndef __DISABLE_LAYOUTSMANAGER
    extern sprog layouts_manager[MAX_LAYOUTSMANAGER_PROGS];
#endif

#ifndef __DISABLE_COLSMANAGER
    extern sprog cols_manager[MAX_COLSMANAGER_PROGS];
#endif

#ifndef __DISABLE_LFSMANAGER
    extern sprog lfs_manager[MAX_LFSMANAGER_PROGS];
#endif

extern const struct ext_math_type ext_math;
extern const struct prog_constants suite_c;


/// FUNCTIONS DECLARATIONS
/// mengine.c
 void toupper_s(char *);
 void tolower_s(char *);
 void asteriskize(char *);
__MATHSUITE void strundsc(const char *, char []);
__MATHSUITE void strboolize(const char *, char []);
__MATHSUITE void strfnm(const char *, char [static MAX_PATH_LENGTH]);
 int  countbits(long);
 int  ucountbits(unsigned long);
 char  *strrev(char *);
 char  *replace(char const * const, char const * const, char const * const);
 bool   file_exists(const char *);
__MATHSUITE bool  readFile(const char [static MAX_PATH_LENGTH]);
__MATHSUITE bool  printFile(const char [static MAX_PATH_LENGTH]);
__MATHSUITE bool  writeFile(const char [static MAX_PATH_LENGTH]);
__MATHSUITE FILE *  checkForFHErrors(const char [static MAX_PATH_LENGTH], char [static 1]);
__MATHSUITE bool  frename(const char [static MAX_PATH_LENGTH], const char [static MAX_PATH_LENGTH]);
 bool matrixUTConv(mpfr_t *restrict, dim_typ, sel_typ *);
 __MATHSUITE void det_3d(mpfr_t, mpfr_t *restrict);
 __MATHSUITE void det_4d(mpfr_t, mpfr_t *restrict);
 __MATHSUITE void det_5d(mpfr_t, mpfr_t *restrict);
__MATHSUITE  void det(mpfr_t, mpfr_t *restrict, dim_typ, bool *);
__MATHSUITE void _matrixTrace(mpfr_t, mpfr_t *restrict, dim_typ);
__MSSHELL_WRAPPER_ __MATHSUITE bool randomMatrix(mpfr_t *, const register dim_typ [static 2]);
__MATHSUITE void transpose(mpfr_t *restrict, mpfr_t *restrict, const register dim_typ [static 2]);
 bool  FattLU(dim_typ, mpfr_t *restrict, mpfr_t *restrict, mpfr_t *);
  __MATHSUITE sel_typ  invertMatrix(mpfr_t *restrict, mpfr_t * restrict, dim_typ);
   __MATHSUITE sel_typ  invertMatrixFast(mpfr_t *restrict, mpfr_t *, dim_typ);
// __MATHSUITE bool  invertMatrix(mpfr_t *restrict, mpfr_t *, dim_typ);
 bool  CoFactor(mpfr_t *restrict, mpfr_t *restrict, dim_typ);
 bool   adjoint(mpfr_t *restrict, mpfr_t *restrict, dim_typ);
 dim_typ  svd(const register dim_typ [static 2],const bool,const bool,mpfr_t eps,mpfr_t tol,mpfr_t *restrict,mpfr_t *restrict,mpfr_t *restrict, mpfr_t *restrict);
 bool  dsvd(mpfr_t *restrict, const register dim_typ [static 2], mpfr_t *restrict, mpfr_t *restrict);
 dim_typ  absor(const register dim_typ[static 2], mpfr_t *restrict, mpfr_t *restrict, mpfr_t **);
__MATHSUITE dim_typ  rank(mpfr_t *restrict, const register dim_typ [static 2]);
__MATHSUITE void  _flushLogBuf(logObj * const);
__MATHSUITE __WINCALL void  _editLog(const char [static MAX_PATH_LENGTH]);
__MATHSUITE void  logCheck(logObj * const , const char *, const char [static MAX_PATH_LENGTH]);
__MATHSUITE void  prependTimeToString(char *);
__MATHSUITE void  printErr(const int, const char *, ...);
__MATHSUITE void msyprintf(const sel_typ col, const char *format, ...);
__MATHSUITE void  fprintf2(FILE *, const char *, ...);
__MATHSUITE void  msprintf(const sel_typ, const char *, ...);
__MATHSUITE bool  scanf2(sel_typ, const char *, ...);
__MATHSUITE bool   _parse(char [], mpfr_t *, exprValList *);
__MATHSUITE  void  _printMatrix(FILE *, mpfr_t **, const register dim_typ [static 2], const char []);
__MATHSUITE  void  printItypMatrix(ityp *, const register dim_typ [static 2]);
__MATHSUITE  void  printBoolMatrix(bool *, const register dim_typ [static 2]);
__MATHSUITE  void  printUShortMatrix(dim_typ *, const register dim_typ [static 2]);
__MATHSUITE  void  printShortMatrix(short  *, const register dim_typ [static 2]);
__MATHSUITE  void  printIntMatrix(int  *, const register dim_typ [static 2]);
__MATHSUITE bool   extractMat(dim_typ);
__MATHSUITE bool   matrixToken(const char [], mpfr_t **, dim_typ *, dim_typ *);

__MATHSUITE void  delAll();

__MATHSUITE void  newFloat(const uint64_t, mpfr_t, const bool, ityp *);
__MATHSUITE void  newBool(const uint64_t, mpfr_t, const bool, bool *);
__MATHSUITE void  newUShort(const uint64_t, mpfr_t, const bool, dim_typ *);
__MATHSUITE void  newShort(const uint64_t, mpfr_t, const bool, short *);
__MATHSUITE void  newInt(const uint64_t, mpfr_t, const bool, int *);

__MATHSUITE void  del(const uint64_t);

__MATHSUITE void  aCheckFloat(const uint64_t, mpfr_t, ityp *); 
__MATHSUITE void  aCheckBool(const uint64_t, mpfr_t, bool *);
__MATHSUITE void  aCheckUShort(const uint64_t, mpfr_t, dim_typ *);
__MATHSUITE void  aCheckShort(const uint64_t, mpfr_t, short *);
__MATHSUITE void  aCheckInt(const uint64_t, mpfr_t, int *);


__MATHSUITE void  truncateVector(const register dim_typ dim, mpfr_t [static dim], ityp [static dim]);
__MATHSUITE ityp *  accessFloatContainerPointer(const mpfr_t);
__MATHSUITE bool *  accessBoolContainerPointer(const mpfr_t);
__MATHSUITE dim_typ *  accessUShortContainerPointer(const mpfr_t);
__MATHSUITE short *  accessShortContainerPointer(const mpfr_t);
__MATHSUITE int *  accessIntContainerPointer(const mpfr_t);

__MSSHELL_WRAPPER_ __MATHSUITE sprog * const  searchProgram(char [static SIGN_STRING]);
 int  cmpfunc(const void *, const void *);
 int mpfr_cmpfunc(const void * a, const void * b);
__MATHSUITE void  _showUsage(sprog * const);
__MATHSUITE void  printUsage(sprog * const);
__MATHSUITE void  prepareToExit(void);
__MATHSUITE void resetLmpMatrix();
// Network functions
#ifndef __DISABLE_SERVER
	__MATHSUITE bool getAddressInfo(char *string, struct sockaddr * addr, socklen_t len);
	__MATHSUITE void closeUnsetSocket(fd_set *, int *);
	__MATHSUITE sel_typ getFamily(const register sel_typ family);
#endif
__MATHSUITE void removeDirEnt(const char [static MAX_PATH_LENGTH], const register dim_typ extNo, char [static extNo][MAX_EXTENSION_LENGTH]);
__MATHSUITE void removeCriticalFiles(const bool);
__MATHSUITE void  _handleCmdLine(char *);
__MATHSUITE bool  _execScriptFiles(const char [static MAX_PATH_LENGTH]);
__MATHSUITE bool  _lfLoader(const char [static MAX_PATH_LENGTH]);
__MATHSUITE bool  _lfCreate(const char [static MAX_PATH_LENGTH]);

#ifdef XMLCALL
	 XMLCALL xmlDoc *   xmlInit(const char [static XML_FILENAMES_LENGTH], xmlXPathContext **);
	 XMLCALL void   xmlExit(const char [static XML_FILENAMES_LENGTH], xmlDoc **, xmlXPathObject **, xmlXPathContext **);
	 XMLCALL bool   xmlWriteInt(xmlXPathObject **, xmlXPathContext *, const char *, const int);
	 XMLCALL bool   xmlWriteBool(xmlXPathObject **, xmlXPathContext *, const char *, const bool);
	 XMLCALL bool   xmlWriteFloat(xmlXPathObject **, xmlXPathContext *, const char *, const float);
	 XMLCALL bool   xmlWriteString(xmlXPathObject **, xmlXPathContext *, const char *, const char [static MAX_XML_FIELDSTRINGS]);
	 XMLCALL bool   xmlGetInt(xmlXPathObject **, xmlXPathContext *, const char *, int *);
	 XMLCALL bool   xmlGetBool(xmlXPathObject **, xmlXPathContext *, const char *, bool *);
	 XMLCALL bool   xmlGetFloat(xmlXPathObject **, xmlXPathContext *, const char *, float *);
	 XMLCALL bool   xmlGetString(xmlXPathObject **, xmlXPathContext *, const char *, char [static MAX_XML_FIELDSTRINGS]);
	__MATHSUITE XMLCALL void  _backupColFile(void);
	__MATHSUITE XMLCALL void  _colFileLoader(const char [static MAX_PATH_LENGTH]);
	__MATHSUITE XMLCALL void  getProgramSettings(dim_typ);
	__MATHSUITE XMLCALL void  resetProgramSettings(layoutObj * const, const char [static MAX_PATH_LENGTH]);
#endif

#ifdef WINOS
	__WINCALL void clrscr();
     __WINCALL HWND WINAPI   GetConsoleWindowNT();
     __WINCALL bool   windowsFileHandler(char *, const char *, const char [static MAX_EXTENSION_LENGTH], bool);
    int inet_pton(int af, const char *src, void *dst);
	const char *inet_ntop(int af, const void *src, char *dst, int);
#else
	 __MATHSUITE   int getch( );
#endif

 const char * const   getFilename(const char path[static MAX_PATH_LENGTH]);
 void  updInfo(const char [static SERVER_LENGTH+5]);

 void   SetColor(const sel_typ);

// __MATHSUITE ityp  MINMAX(const register dim_typ dim, const ityp [static dim], const bool, dim_typ *);
__MATHSUITE void MINMAX(ityp *, ityp *, unsigned dim, ityp [static dim], dim_typ *, dim_typ *);
__MATHSUITE void  MPFR_MINMAX(mpfr_t, mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t, mpfr_t);
__MATHSUITE  ityp   getDiffTime(struct timeval * tvBegin);
__MATHSUITE bool   isDomainForbidden(mpfr_t, bool);
 void    free_foreach(mpfr_t **, const dim_typ, const register dim_typ [static 2], bool mode);
 bool   checkErrMem(const void *);
__MATHSUITE bool   matrixAlloc(mpfr_t **, const register dim_typ [static 2]);
__MATHSUITE void   _matrixFree(mpfr_t **, const register dim_typ [static 2], bool);
__MATHSUITE bool   equalMatrix(mpfr_t **, mpfr_t *, const register dim_typ [static 2], const bool);
__MATHSUITE  void  resetLmpMatrix(void);
__MATHSUITE void   _flushMemoizersBuffers(sel_typ);
__MATHSUITE void   flushAllMemoizersBuffers(void);
__MATHSUITE dim_typ   selectListItem(dim_typ dim, bool, const char *, const char [static dim][MIN_STRING]);

#ifndef __DISABLE_DATABASE
	__MATHSUITE bool writeVar(const char *, const char *, mpfr_t);
	__MATHSUITE bool readVar(const char *, const char *, mpfr_t);
	__MATHSUITE int findVarList(const char *);
	__MATHSUITE bool readVarList(const char *, dim_typ *, char [][VAR_MAXNAMELENGTH], char [][VAR_MAXVALLENGTH]);
	__MATHSUITE bool deleteVarList(const char *);
	
	__MATHSUITE bool writeMat(const char *, const char *);
	__MATHSUITE int findMat(const char *); 
	__MATHSUITE bool readMat(const char *, char *);
	__MATHSUITE bool deleteMat(const char *, const bool);
	
	// __MATHSUITE bool readUser(const register int, char [static 6][USER_MAXFIELDLENGTH]);
	__MATHSUITE bool readUser(const register int, userObj **);
	__MATHSUITE bool readUserFromUserName(const char [static USER_MAXIDENTIFIERLENGTH], userObj **);
	__MATHSUITE bool readUserFromEmail(const char [static USER_MAXIDENTIFIERLENGTH], userObj **);
	__MATHSUITE bool writeUser(userObj *restrict, const bool);
	
	__MATHSUITE matrixObj * getMatrixFromMatrixList(const register dim_typ);
	__MATHSUITE bool pushMatrixList(matrixObj *);
	__MATHSUITE void delAllMatrixFromMatrixList(void);	
	__MATHSUITE void delTailMatrixList(void);
	__MATHSUITE session * accessCurrentSession();
	__MATHSUITE void dbEstablishConnection(const sel_typ, char *[static 4]);
	
#endif

__MATHSUITE void  viewProgramSettings(dim_typ);
__MATHSUITE void  setProgramSettings(dim_typ);
__MATHSUITE bool   catchPause();
__MATHSUITE void  logPrint(logObj * const);
__MATHSUITE void  logWrite(FILE *, logObj * const);
__MATHSUITE void  getVarListEx(FILE *, dim_typ);
__MATHSUITE const bool   requires(mpfr_t, const char *, const char *, const char *, const bool);
__MATHSUITE bool  insertDims(dim_typ *, dim_typ *);
__MATHSUITE bool  insertDim(dim_typ *, bool);
__MATHSUITE void  equalSpecialMatrix(int, int);
__MATHSUITE void  sigexit(int);
__MSSHELL_WRAPPER_ __MATHSUITE void showNewtonDifferenceTable(dim_typ n, mpfr_t [static n], mpfr_t [access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim], bool);
__MATHSUITE bool  insertNMMatrix(mpfr_t **, const register dim_typ [static 2]);
__MATHSUITE volatile char  insertElement(mpfr_t *restrict, const register dim_typ [static 2], const register dim_typ, const bool, const bool);
__MATHSUITE volatile bool  checkBackTracking(volatile char, dim_typ *);
__MATHSUITE volatile  sel_typ checkBackTracking2(volatile char, dim_typ *, dim_typ *, dim_typ *, dim_typ);
__MATHSUITE bool   enterMatrix(mpfr_t **, dim_typ *, dim_typ *, bool, bool);


/// programs.c
__MSSHELL_WRAPPER_ void __apnt changeProgramSettings(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void basicCalculator(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void __apnt advancedCalculator(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void __apnt algebraOperations(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void __apnt multiUser(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ void __apnt mssManager(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ __MATHSUITE void    operationsGroupMenu(dim_typ dim, sprog [static dim], const char [static INFO_STRING], bool);
__MSSHELL_WRAPPER_ __MATHSUITE void  progInfo(sel_typ);
 void   freeExprEvalLists();
 void   refreshExprEvalVarList(dim_typ);
 void   refreshExprEvalLists();
__MATHSUITE void   setCurrentMatrix(dim_typ);
__MATHSUITE void   setCurrentLog(dim_typ);
__MSSHELL_WRAPPER_ __MATHSUITE void  setDefaults();

/// ext_math.c

__MATHSUITE char * const   binaryAlgSum(ityp, ityp, bool);
__MATHSUITE char * const   binNumComp(ityp);
__MATHSUITE void  ___cabs(mpfr_t, mpfr_t *restrict, const register sel_typ);
__MATHSUITE void  _complexAdd(mpfr_t *restrict, mpfr_t [static 2]);
__MATHSUITE void  _complexSub(mpfr_t *restrict, mpfr_t [static 2]);
__MATHSUITE void  _complexMul(mpfr_t *restrict, mpfr_t [static 2]);
__MATHSUITE void  _complexDiv(mpfr_t *restrict, mpfr_t [static 2]);
__MATHSUITE void  _quaternionsAdd(mpfr_t *restrict, mpfr_t [static MAX_QUATERNIONS_UNITS]);
__MATHSUITE void  _quaternionsSub(mpfr_t *restrict, mpfr_t [static MAX_QUATERNIONS_UNITS]);
__MATHSUITE void  _quaternionsMul(mpfr_t *restrict, mpfr_t [static MAX_QUATERNIONS_UNITS]);
__MATHSUITE void  _quaternionsDiv(mpfr_t *restrict, mpfr_t [static MAX_QUATERNIONS_UNITS]);
__MATHSUITE void  _octonionsAdd(mpfr_t *restrict, mpfr_t [static MAX_OCTONIONS_UNITS]);
__MATHSUITE void  _octonionsSub(mpfr_t *restrict, mpfr_t [static MAX_OCTONIONS_UNITS]);
__MATHSUITE void  _octonionsMul(mpfr_t *restrict, mpfr_t [static MAX_OCTONIONS_UNITS]);
__MATHSUITE void  _octonionsDiv(mpfr_t *restrict, mpfr_t [static MAX_OCTONIONS_UNITS]);
__MATHSUITE void  _sedenionsAdd(mpfr_t *restrict, mpfr_t [static MAX_SEDENIONS_UNITS]);
__MATHSUITE void  _sedenionsSub(mpfr_t *restrict, mpfr_t [static MAX_SEDENIONS_UNITS]);
__MATHSUITE void  _sedenionsMul(mpfr_t *restrict, mpfr_t [static MAX_SEDENIONS_UNITS]);
__MATHSUITE void  _sedenionsDiv(mpfr_t *restrict, mpfr_t [static MAX_SEDENIONS_UNITS]);
__MATHSUITE bool  _secondGradeEquationSolver(mpfr_t *restrict , mpfr_t [static 2]);
 const char * const  getDayName(sel_typ);
 const char * const  getMonthName(sel_typ);
 const sel_typ  getDayNumber(sel_typ, sel_typ, uint64_t);
 void  getRomanNumber(dim_typ, char [static ROMAN_NUMBER_STRING]);
 void  mpow(mpfr_t, mpfr_t, mpfr_t);
 void  mpow2(mpfr_t, mpfr_t, mpfr_t);
__MATHSUITE  void  getPascalTriangle(mpfr_t);
 bool   isPrimeHISP(register mpfr_t);
 bool   isPrimeLIFP(register mpfr_t);
__MATHSUITE void   prime_N_Number(register mpfr_t, register mpfr_t);
__MATHSUITE void   N_prime_Number(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_fibnc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_primr(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_fpnsum(mpfr_t, register mpfr_t);
__MATHSUITE bool   trigonometric_domain(register mpfr_t, register mpfr_t, register mpfr_t);
__MATHSUITE void  stirling(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_fibo(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_fact(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_dfact(mpfr_t, register mpfr_t);
__MATHSUITE void  sfact_even(mpfr_t, register mpfr_t);
__MATHSUITE void  sfact_odd(mpfr_t, register mpfr_t);
 void  perm(mpfr_t, register mpfr_t);
 void  perm_rep(mpfr_t, register mpfr_t, register mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
 void  kperm(mpfr_t, register mpfr_t, register mpfr_t);
 void  kperm_rep(mpfr_t, register mpfr_t, register mpfr_t);
 void  comb(mpfr_t, register mpfr_t, register mpfr_t);
 void  comb_rep(mpfr_t, register mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_gsum(mpfr_t, register mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_gasum(mpfr_t, register mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_asum(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_fsum(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_fasum(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_sfasum(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_fnnsum(mpfr_t, register mpfr_t);
__MATHSUITE void   answer_to_life_the_universe_and_everything(mpfr_t);
__MATHSUITE ityp  summation(uint64_t dim, bool, ityp [static dim]);
__MATHSUITE ityp  productory(uint64_t dim, bool, ityp [static dim]);
__MATHSUITE void  mpfr_summation(mpfr_t, mpfr_t dim, bool, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  mpfr_product(mpfr_t, mpfr_t dim, bool, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE ityp  mbase_mean(const register dim_typ dim, ityp [static dim]);
__MATHSUITE void  math_mean(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_mode(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_variance(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_variance2(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_covariance(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_covariance2(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_stdcod(mpfr_t rop, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_stdcod2(mpfr_t rop, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_stddev(mpfr_t rop, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_stddev2(mpfr_t rop, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_pearson(mpfr_t rop, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_pearson2(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE bool  math_outlier(mpfr_t dim, mpfr_t, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE bool  math_outlier2(mpfr_t dim, mpfr_t, ityp, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_geomean(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_armean(mpfr_t, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_powmean(mpfr_t, mpfr_t dim, mpfr_t power, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_scale(mpfr_t rop, mpfr_t dim, mpfr_t [static mpfr_get_ui(dim, MPFR_RNDN)]);
__MATHSUITE void  math_first_quartile(mpfr_t rop, const register dim_typ dim, mpfr_t [static dim]);
__MATHSUITE void  math_median(mpfr_t rop, const register dim_typ dim, mpfr_t [static dim]);
__MATHSUITE void  math_third_quartile(mpfr_t rop, const register dim_typ dim, mpfr_t [static dim]);

__MATHSUITE dim_typ randomizedSelect(mpfr_t, mpfr_t [], dim_typ, dim_typ, dim_typ);
__MATHSUITE dim_typ partition(mpfr_t [], dim_typ, dim_typ);
__MATHSUITE dim_typ randomizedPartition(mpfr_t [], dim_typ, dim_typ);

 int64_t   powi(register int64_t, register int64_t); // calculates x^y
 int64_t   changeBase(register int, sel_typ, sel_typ);
 void  math_GCD(mpfr_t, mpfr_t, mpfr_t); // Euclide's Algorythm
 void  math_lcm(mpfr_t, mpfr_t, mpfr_t);



__MATHSUITE ityp  exp10(register ityp);
__MATHSUITE ityp  expc(register ityp);
__MATHSUITE ityp  exp10c(register ityp);
__MATHSUITE ityp  exp2c(register ityp);
 ityp  logbN(register ityp, register ityp);
__MATHSUITE ityp  logc(register ityp);
__MATHSUITE ityp  log10c(register ityp);
__MATHSUITE ityp  log2c(register ityp);
__MATHSUITE ityp  log1pc(register ityp);
__MATHSUITE ityp  log101p(register ityp);
__MATHSUITE ityp  log101pc(register ityp);
__MATHSUITE ityp  log21p(register ityp);
__MATHSUITE ityp  log21pc(register ityp);
__MATHSUITE ityp  deg(register ityp);
__MATHSUITE ityp  rad(register ityp);
__MATHSUITE ityp  csc(register ityp);
__MATHSUITE ityp  sec(register ityp);
__MATHSUITE ityp  cot(register ityp);
__MATHSUITE ityp  csch(register ityp);
__MATHSUITE ityp  sech(register ityp);
__MATHSUITE ityp  coth(register ityp);
__MATHSUITE ityp  acsc(register ityp);
__MATHSUITE ityp  asec(register ityp);
__MATHSUITE ityp  acot(register ityp);
__MATHSUITE ityp  acsch(register ityp);
__MATHSUITE ityp  asech(register ityp);
__MATHSUITE ityp  acoth(register ityp);
__MATHSUITE ityp  hsin(register ityp);
__MATHSUITE ityp  hsinh(register ityp);
__MATHSUITE ityp  qsin(register ityp);
__MATHSUITE ityp  qsinh(register ityp);
__MATHSUITE ityp  hcos(register ityp);
__MATHSUITE ityp  hcosh(register ityp);
__MATHSUITE ityp  qcos(register ityp);
__MATHSUITE ityp  qcosh(register ityp);
__MATHSUITE ityp  hcsc(register ityp);
__MATHSUITE ityp  hcsch(register ityp);
__MATHSUITE ityp  qcsc(register ityp);
__MATHSUITE ityp  qcsch(register ityp);
__MATHSUITE ityp  hsec(register ityp);
__MATHSUITE ityp  hsech(register ityp);
__MATHSUITE ityp  qsec(register ityp);
__MATHSUITE ityp  qsech(register ityp);
__MATHSUITE ityp  htan(register ityp);
__MATHSUITE ityp  htanh(register ityp);
__MATHSUITE ityp  qtan(register ityp);
__MATHSUITE ityp  qtanh(register ityp);
__MATHSUITE ityp  hcot(register ityp);
__MATHSUITE ityp  hcoth(register ityp);
__MATHSUITE ityp  qcot(register ityp);
__MATHSUITE ityp  qcoth(register ityp);
__MATHSUITE ityp  vsin(register ityp);
__MATHSUITE ityp  cvsin(register ityp);
__MATHSUITE ityp  vcos(register ityp);
__MATHSUITE ityp  cvcos(register ityp);
__MATHSUITE ityp  hvsin(register ityp);
__MATHSUITE ityp  hcvsin(register ityp);
__MATHSUITE ityp  hvcos(register ityp);
__MATHSUITE ityp  hcvcos(register ityp);
__MATHSUITE ityp  qvsin(register ityp);
__MATHSUITE ityp  qcvsin(register ityp);
__MATHSUITE ityp  qvcos(register ityp);
__MATHSUITE ityp  qcvcos(register ityp);
__MATHSUITE ityp  vsinh(register ityp);
__MATHSUITE ityp  cvsinh(register ityp);
__MATHSUITE ityp  vcosh(register ityp);
__MATHSUITE ityp  cvcosh(register ityp);
__MATHSUITE ityp  hvsinh(register ityp);
__MATHSUITE ityp  hcvsinh(register ityp);
__MATHSUITE ityp  hvcosh(register ityp);
__MATHSUITE ityp  hcvcosh(register ityp);
__MATHSUITE ityp  qvsinh(register ityp);
__MATHSUITE ityp  qcvsinh(register ityp);
__MATHSUITE ityp  qvcosh(register ityp);
__MATHSUITE ityp  qcvcosh(register ityp);
__MATHSUITE ityp  esec(register ityp);
__MATHSUITE ityp  ecsc(register ityp);
__MATHSUITE ityp  esech(register ityp);
__MATHSUITE ityp  ecsch(register ityp);
__MATHSUITE ityp  hesec(register ityp);
__MATHSUITE ityp  hecsc(register ityp);
__MATHSUITE ityp  hesech(register ityp);
__MATHSUITE ityp  hecsch(register ityp);
__MATHSUITE ityp  qesec(register ityp);
__MATHSUITE ityp  qecsc(register ityp);
__MATHSUITE ityp  qesech(register ityp);
__MATHSUITE ityp  qecsch(register ityp);
__MATHSUITE ityp  sinc(register ityp);
__MATHSUITE ityp  sinch(register ityp);
__MATHSUITE ityp  hsinc(register ityp);
__MATHSUITE ityp  hsinch(register ityp);
__MATHSUITE ityp  qsinc(register ityp);
__MATHSUITE ityp  qsinch(register ityp);
__MATHSUITE ityp  cosc(register ityp);
__MATHSUITE ityp  cosch(register ityp);
__MATHSUITE ityp  hcosc(register ityp);
__MATHSUITE ityp  hcosch(register ityp);
__MATHSUITE ityp  qcosc(register ityp);
__MATHSUITE ityp  qcosch(register ityp);
__MATHSUITE ityp  secc(register ityp);
__MATHSUITE ityp  secch(register ityp);
__MATHSUITE ityp  hsecc(register ityp);
__MATHSUITE ityp  hsecch(register ityp);
__MATHSUITE ityp  qsecc(register ityp);
__MATHSUITE ityp  qsecch(register ityp);
__MATHSUITE ityp  cscc(register ityp);
__MATHSUITE ityp  cscch(register ityp);
__MATHSUITE ityp  hcscc(register ityp);
__MATHSUITE ityp  hcscch(register ityp);
__MATHSUITE ityp  qcscc(register ityp);
__MATHSUITE ityp  qcscch(register ityp);
__MATHSUITE ityp  tanc(register ityp);
__MATHSUITE ityp  tanch(register ityp);
__MATHSUITE ityp  htanc(register ityp);
__MATHSUITE ityp  htanch(register ityp);
__MATHSUITE ityp  qtanc(register ityp);
__MATHSUITE ityp  qtanch(register ityp);
__MATHSUITE ityp  cotc(register ityp);
__MATHSUITE ityp  cotch(register ityp);
__MATHSUITE ityp  hcotc(register ityp);
__MATHSUITE ityp  hcotch(register ityp);
__MATHSUITE ityp  qcotc(register ityp);
__MATHSUITE ityp  qcotch(register ityp);


__MATHSUITE void  mpfr_expc(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_exp10c(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_exp2c(mpfr_t, mpfr_t);
 void  mpfr_logbN(mpfr_t, mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_logc(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_log10c(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_log2c(mpfr_t, mpfr_t);
__MATHSUITE void  mswrap_log1p(mpfr_t, mpfr_t);
__MATHSUITE void  mswrap_floor(mpfr_t, mpfr_t);
__MATHSUITE void  mswrap_ceil(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_log1pc(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_log101p(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_log101pc(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_log21p(mpfr_t, mpfr_t);
__MATHSUITE void  mpfr_log21pc(mpfr_t, mpfr_t);
 void  rootnX(mpfr_t, mpfr_t, mpfr_t);

__MATHSUITE double complex  cexpc(register double complex);
__MATHSUITE double complex  cexp10(register double complex);
__MATHSUITE double complex  cexp10c(register double complex);
__MATHSUITE double complex  cexp2(register double complex);
__MATHSUITE double complex  cexp2c(register double complex);
 double complex  clogbN(register double complex, register double complex);
__MATHSUITE double complex  clogc(register double complex);
__MATHSUITE double complex  clog10(register double complex);
__MATHSUITE double complex  clog10c(register double complex);
__MATHSUITE double complex  clog2(register double complex);
__MATHSUITE double complex  clog2c(register double complex);
__MATHSUITE double complex  clog1p(register double complex);
__MATHSUITE double complex  clog1pc(register double complex);
__MATHSUITE double complex  clog101p(register double complex);
__MATHSUITE double complex  clog101pc(register double complex);
__MATHSUITE double complex  clog21p(register double complex);
__MATHSUITE double complex  clog21pc(register double complex);
 double complex  ccbrt(register double complex);
 double complex  crootnX(register double complex, register double complex);

__MATHSUITE void  cel_fah(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  cel_kel(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  cel_rank(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  cel_rea(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  cel_new(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  cel_del(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  cel_rom(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  fah_kel(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  fah_rank(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  fah_rea(mpfr_t, mpfr_t, const bool);
__MATHSUITE void  rea_rank(mpfr_t, mpfr_t, const bool);

__MATHSUITE void  mswrap_sin(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_cos(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_sinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_cosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_sec(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_sech(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_csc(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_csch(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_asin(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_asinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_acos(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_acosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_tan(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_tanh(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_atan(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_atanh(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_cot(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_coth(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_log(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_log10(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_log2(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_exp(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_exp10(mpfr_t, register mpfr_t);
__MATHSUITE void  mswrap_exp2(mpfr_t, register mpfr_t);
__MATHSUITE ityp  cotan(const register ityp);
__MATHSUITE ityp  degd(const register ityp);
__MATHSUITE void  mpfr_deg(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_rad(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_acsc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_asec(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_acot(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_acsch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_asech(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_acoth(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hsin(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hsinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qsin(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qsinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcos(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcos(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcsc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcsch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcsc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcsch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hsec(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hsech(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qsec(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qsech(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_htan(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_htanh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qtan(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qtanh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcot(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcoth(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcot(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcoth(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_vsin(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cvsin(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_vcos(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cvcos(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hvsin(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcvsin(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hvcos(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcvcos(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qvsin(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcvsin(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qvcos(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcvcos(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_vsinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cvsinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_vcosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cvcosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hvsinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcvsinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hvcosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcvcosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qvsinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcvsinh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qvcosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcvcosh(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_esec(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_ecsc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_esech(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_ecsch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hesec(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hecsc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hesech(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hecsch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qesec(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qecsc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qesech(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qecsch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_sinc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_sinch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hsinc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hsinch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qsinc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qsinch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cosc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cosch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcosc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcosch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcosc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcosch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_secc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_secch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hsecc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hsecch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qsecc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qsecch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cscc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cscch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcscc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcscch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcscc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcscch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_tanc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_tanch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_htanc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_htanch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qtanc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qtanch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cotc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_cotch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcotc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_hcotch(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcotc(mpfr_t, register mpfr_t);
__MATHSUITE void  mpfr_qcotch(mpfr_t, register mpfr_t);

__MATHSUITE double complex  ccsc(register double complex);
__MATHSUITE double complex  csec(register double complex);
__MATHSUITE double complex  ccot(register double complex);
__MATHSUITE double complex  ccsch(register double complex);
__MATHSUITE double complex  csech(register double complex);
__MATHSUITE double complex  ccoth(register double complex);
__MATHSUITE double complex  cacsc(register double complex);
__MATHSUITE double complex  casec(register double complex);
__MATHSUITE double complex  cacot(register double complex);
__MATHSUITE double complex  cacsch(register double complex);
__MATHSUITE double complex  casech(register double complex);
__MATHSUITE double complex  cacoth(register double complex);
__MATHSUITE double complex  chsin(register double complex);
__MATHSUITE double complex  chsinh(register double complex);
__MATHSUITE double complex  cqsin(register double complex);
__MATHSUITE double complex  cqsinh(register double complex);
__MATHSUITE double complex  chcos(register double complex);
__MATHSUITE double complex  chcosh(register double complex);
__MATHSUITE double complex  cqcos(register double complex);
__MATHSUITE double complex  cqcosh(register double complex);
__MATHSUITE double complex  chcsc(register double complex);
__MATHSUITE double complex  chcsch(register double complex);
__MATHSUITE double complex  cqcsc(register double complex);
__MATHSUITE double complex  cqcsch(register double complex);
__MATHSUITE double complex  chsec(register double complex);
__MATHSUITE double complex  chsech(register double complex);
__MATHSUITE double complex  cqsec(register double complex);
__MATHSUITE double complex  cqsech(register double complex);
__MATHSUITE double complex  chtan(register double complex);
__MATHSUITE double complex  chtanh(register double complex);
__MATHSUITE double complex  cqtan(register double complex);
__MATHSUITE double complex  cqtanh(register double complex);
__MATHSUITE double complex  chcot(register double complex);
__MATHSUITE double complex  chcoth(register double complex);
__MATHSUITE double complex  cqcot(register double complex);
__MATHSUITE double complex  cqcoth(register double complex);
__MATHSUITE double complex  cpxvsin(register double complex);
__MATHSUITE double complex  ccvsin(register double complex);
__MATHSUITE double complex  cpxvcos(register double complex);
__MATHSUITE double complex  ccvcos(register double complex);
__MATHSUITE double complex  chvsin(register double complex);
__MATHSUITE double complex  chcvsin(register double complex);
__MATHSUITE double complex  chvcos(register double complex);
__MATHSUITE double complex  chcvcos(register double complex);
__MATHSUITE double complex  cqvsin(register double complex);
__MATHSUITE double complex  cqcvsin(register double complex);
__MATHSUITE double complex  cqvcos(register double complex);
__MATHSUITE double complex  cqcvcos(register double complex);
__MATHSUITE double complex  cpxvsinh(register double complex);
__MATHSUITE double complex  ccvsinh(register double complex);
__MATHSUITE double complex  cpxvcosh(register double complex);
__MATHSUITE double complex  ccvcosh(register double complex);
__MATHSUITE double complex  chvsinh(register double complex);
__MATHSUITE double complex  chcvsinh(register double complex);
__MATHSUITE double complex  chvcosh(register double complex);
__MATHSUITE double complex  chcvcosh(register double complex);
__MATHSUITE double complex  cqvsinh(register double complex);
__MATHSUITE double complex  cqcvsinh(register double complex);
__MATHSUITE double complex  cqvcosh(register double complex);
__MATHSUITE double complex  cqcvcosh(register double complex);
__MATHSUITE double complex  cesec(register double complex);
__MATHSUITE double complex  cecsc(register double complex);
__MATHSUITE double complex  cesech(register double complex);
__MATHSUITE double complex  cecsch(register double complex);
__MATHSUITE double complex  chesec(register double complex);
__MATHSUITE double complex  checsc(register double complex);
__MATHSUITE double complex  chesech(register double complex);
__MATHSUITE double complex  checsch(register double complex);
__MATHSUITE double complex  cqesec(register double complex);
__MATHSUITE double complex  cqecsc(register double complex);
__MATHSUITE double complex  cqesech(register double complex);
__MATHSUITE double complex  cqecsch(register double complex);
__MATHSUITE double complex  csinc(register double complex);
__MATHSUITE double complex  csinch(register double complex);
__MATHSUITE double complex  chsinc(register double complex);
__MATHSUITE double complex  chsinch(register double complex);
__MATHSUITE double complex  cqsinc(register double complex);
__MATHSUITE double complex  cqsinch(register double complex);
__MATHSUITE double complex  ccosc(register double complex);
__MATHSUITE double complex  ccosch(register double complex);
__MATHSUITE double complex  chcosc(register double complex);
__MATHSUITE double complex  chcosch(register double complex);
__MATHSUITE double complex  cqcosc(register double complex);
__MATHSUITE double complex  cqcosch(register double complex);
__MATHSUITE double complex  csecc(register double complex);
__MATHSUITE double complex  csecch(register double complex);
__MATHSUITE double complex  chsecc(register double complex);
__MATHSUITE double complex  chsecch(register double complex);
__MATHSUITE double complex  cqsecc(register double complex);
__MATHSUITE double complex  cqsecch(register double complex);
__MATHSUITE double complex  ccscc(register double complex);
__MATHSUITE double complex  ccscch(register double complex);
__MATHSUITE double complex  chcscc(register double complex);
__MATHSUITE double complex  chcscch(register double complex);
__MATHSUITE double complex  cqcscc(register double complex);
__MATHSUITE double complex  cqcscch(register double complex);
__MATHSUITE double complex  ctanc(register double complex);
__MATHSUITE double complex  ctanch(register double complex);
__MATHSUITE double complex  chtanc(register double complex);
__MATHSUITE double complex  chtanch(register double complex);
__MATHSUITE double complex  cqtanc(register double complex);
__MATHSUITE double complex  cqtanch(register double complex);
__MATHSUITE double complex  ccotc(register double complex);
__MATHSUITE double complex  ccotch(register double complex);
__MATHSUITE double complex  chcotc(register double complex);
__MATHSUITE double complex  chcotch(register double complex);
__MATHSUITE double complex  cqcotc(register double complex);
__MATHSUITE double complex  cqcotch(register double complex);
__MATHSUITE bool   isEqualMatrix(mpfr_t *, mpfr_t *, const register dim_typ [static 2]);
__MATHSUITE  void   norms(mpfr_t, mpfr_t *, dim_typ);
__MATHSUITE  void   norm(mpfr_t, mpfr_t *, dim_typ, bool);
__MATHSUITE  void  newtonDifferenceTable(dim_typ, mpfr_t [access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim], bool);

#include "jburkardt/headers/jburkardt.h"
#include "jburkardt/headers/alpert_rule.h"
#include "jburkardt/headers/asa005.h"
#include "jburkardt/headers/asa006.h"
#include "jburkardt/headers/asa032.h"
#include "jburkardt/headers/asa047.h"
#include "jburkardt/headers/asa058.h"
#include "jburkardt/headers/asa063.h"
#include "jburkardt/headers/asa066.h"
#include "jburkardt/headers/asa076.h"
#include "jburkardt/headers/asa091.h"
#include "jburkardt/headers/asa103.h"
#include "jburkardt/headers/asa109.h"
#include "jburkardt/headers/asa111.h"
#include "jburkardt/headers/asa121.h"
#include "jburkardt/headers/asa144.h"
#include "jburkardt/headers/asa147.h"
#include "jburkardt/headers/asa152.h"
#include "jburkardt/headers/asa172.h"
#include "jburkardt/headers/asa183.h"
#include "jburkardt/headers/asa226.h"
#include "jburkardt/headers/asa243.h"
#include "jburkardt/headers/asa245.h"
#include "jburkardt/headers/asa266.h"
#include "jburkardt/headers/asa310.h"
#include "jburkardt/headers/asa314.h"
#include "jburkardt/headers/ball_integrals.h"
#include "jburkardt/headers/bernstein_polynomial.h"
#include "jburkardt/headers/bessel.h"
#include "jburkardt/headers/blas1_d.h"
#include "jburkardt/headers/brent.h"
#include "jburkardt/headers/brownian_motion_simulation.h"
#include "jburkardt/headers/burgers_solution.h"
#include "jburkardt/headers/c_calls_f77.h"
#include "jburkardt/headers/ccn_rule.h"
#include "jburkardt/headers/cpv.h"
#include "jburkardt/headers/cdflib.h"
#include "jburkardt/headers/chebyshev_polynomial.h"
#include "jburkardt/headers/chebyshev_series.h"
#include "jburkardt/headers/circle_arc_grid.h"
#include "jburkardt/headers/circle_integrals.h"
#include "jburkardt/headers/circle_rule.h"
#include "jburkardt/headers/circle_segment.h"
#include "jburkardt/headers/clenshaw_curtis_rule.h"
#include "jburkardt/headers/colored_noise.h"
#include "jburkardt/headers/combination_lock.h"
#include "jburkardt/headers/combo.h"
#include "jburkardt/headers/cordic.h"
#include "jburkardt/headers/correlation.h"
#include "jburkardt/headers/cosine_transform.h"
#include "jburkardt/headers/cube_arbq_rule.h"
#include "jburkardt/headers/cube_exactness.h"
#include "jburkardt/headers/cube_felippa_rule.h"
#include "jburkardt/headers/cycle_brent.h"
#include "jburkardt/headers/cycle_floyd.h"
#include "jburkardt/headers/differ.h"
#include "jburkardt/headers/disk_integrals.h"
#include "jburkardt/headers/divdif.h"
#include "jburkardt/headers/duel_simulation.h"
#include "jburkardt/headers/ellipse_grid.h"
#include "jburkardt/headers/ellipse_monte_carlo.h"
#include "jburkardt/headers/ellipsoid_monte_carlo.h"
#include "jburkardt/headers/exactness.h"
#include "jburkardt/headers/fd1d_bvp.h"
#include "jburkardt/headers/fem_basis.h"
#include "jburkardt/headers/fem1d_bvp_linear.h"
#include "jburkardt/headers/fem1d_bvp_quadratic.h"
#include "jburkardt/headers/fem1d_heat_steady.h"
#include "jburkardt/headers/fem1d_lagrange.h"
#include "jburkardt/headers/fem1d_pack.h"
#include "jburkardt/headers/fem1d_project.h"
#include "jburkardt/headers/fem2d_bvp_linear.h"
#include "jburkardt/headers/fem2d_bvp_quadratic.h"
#include "jburkardt/headers/fem2d_serene.h"
#include "jburkardt/headers/filon.h"
#include "jburkardt/headers/geometry.h"
#include "jburkardt/headers/haar.h"
#include "jburkardt/headers/hermite.h"
#include "jburkardt/headers/hermite_polynomial.h"
#include "jburkardt/headers/hermite_product_polynomial.h"
#include "jburkardt/headers/hermite_rule.h"
#include "jburkardt/headers/high_card_simulation.h"
#include "jburkardt/headers/hyperball_integrals.h"
#include "jburkardt/headers/hypercube_grid.h"
#include "jburkardt/headers/hypersphere_integrals.h"
#include "jburkardt/headers/hypersphere_properties.h"
#include "jburkardt/headers/i4lib.h"
#include "jburkardt/headers/image_components.h"
#include "jburkardt/headers/image_denoise.h"
#include "jburkardt/headers/interp.h"
#include "jburkardt/headers/isbn.h"
#include "jburkardt/headers/ising_2d_simulation.h"
#include "jburkardt/headers/jacobi.h"
#include "jburkardt/headers/jacobi_eigenvalue.h"
#include "jburkardt/headers/jacobi_polynomial.h"
#include "jburkardt/headers/knapsack_01.h"
#include "jburkardt/headers/kronrod.h"
#include "jburkardt/headers/lagrange_approx_1d.h"
#include "jburkardt/headers/lagrange_interp_1d.h"
#include "jburkardt/headers/lagrange_interp_nd.h"
#include "jburkardt/headers/laguerre_polynomial.h"
#include "jburkardt/headers/laplacian.h"
#include "jburkardt/headers/latin_cover.h"
#include "jburkardt/headers/latin_random.h"
#include "jburkardt/headers/latinize.h"
#include "jburkardt/headers/lebesgue.h"
#include "jburkardt/headers/legendre_polynomial.h"
#include "jburkardt/headers/life_serial.h"
#include "jburkardt/headers/line_cvt_lloyd.h"
#include "jburkardt/headers/line_fekete_rule.h"
#include "jburkardt/headers/line_grid.h"
#include "jburkardt/headers/line_integrals.h"
#include "jburkardt/headers/line_ncc_rule.h"
#include "jburkardt/headers/line_nco_rule.h"
#include "jburkardt/headers/linpack_d.h"
#include "jburkardt/headers/llsq.h"
#include "jburkardt/headers/lobatto_polynomial.h"
#include "jburkardt/headers/lorenz_ode.h"
#include "jburkardt/headers/matrix_exponential.h"
#include "jburkardt/headers/md.h"
#include "jburkardt/headers/mgs.h"
#include "jburkardt/headers/monomial.h"
#include "jburkardt/headers/multigrid_poisson_1d.h"
#include "jburkardt/headers/naca.h"
#include "jburkardt/headers/navier_stokes_2d_exact.h"
#include "jburkardt/headers/nearest_interp_1d.h"
#include "jburkardt/headers/newton_interp_1d.h"
#include "jburkardt/headers/nintlib.h"
#include "jburkardt/headers/normal.h"
#include "jburkardt/headers/ode.h"
#include "jburkardt/headers/padua.h"
#include "jburkardt/headers/partition_problem.h"
#include "jburkardt/headers/patterson_rule.h"
#include "jburkardt/headers/pce_ode_hermite.h"
#include "jburkardt/headers/pdflib.h"
#include "jburkardt/headers/piecewise_linear_product_integral.h"
#include "jburkardt/headers/pink_noise.h"
#include "jburkardt/headers/poisson_openmp.h"
#include "jburkardt/headers/polpak.h"
#include "jburkardt/headers/polygon_grid.h"
#include "jburkardt/headers/polygon_integrals.h"
#include "jburkardt/headers/polygon_triangulate.h"
#include "jburkardt/headers/power_method.h"
#include "jburkardt/headers/power_rule.h"
#include "jburkardt/headers/prob.h"
#include "jburkardt/headers/pwl_approx_1d.h"
#include "jburkardt/headers/pwl_interp_1d.h"
#include "jburkardt/headers/pwl_interp_2d.h"
#include "jburkardt/headers/pyramid_exactness.h"
#include "jburkardt/headers/pyramid_grid.h"
#include "jburkardt/headers/quadmom.h"
#include "jburkardt/headers/r8lib.h"
#include "jburkardt/headers/rbf_interp_1d.h"
#include "jburkardt/headers/rbf_interp_2d.h"
#include "jburkardt/headers/rbf_interp_nd.h"
#include "jburkardt/headers/rk4.h"
#include "jburkardt/headers/rkf45.h"
#include "jburkardt/headers/satisfy.h"
#include "jburkardt/headers/sftpack.h"
#include "jburkardt/headers/sgmga.h"
#include "jburkardt/headers/shepard_interp_1d.h"
#include "jburkardt/headers/shepard_interp_2d.h"
#include "jburkardt/headers/shepard_interp_nd.h"
#include "jburkardt/headers/simplex_coordinates.h"
#include "jburkardt/headers/simplex_gm_rule.h"
#include "jburkardt/headers/simplex_grid.h"
#include "jburkardt/headers/sine_transform.h"
#include "jburkardt/headers/snakes_and_ladders.h"
#include "jburkardt/headers/sparse_grid_cc.h"
#include "jburkardt/headers/sparse_grid_hw.h"
#include "jburkardt/headers/sparse_interp_nd.h"
#include "jburkardt/headers/sphere_fibonacci_grid.h"
#include "jburkardt/headers/sphere_grid.h"
#include "jburkardt/headers/sphere_integrals.h"
#include "jburkardt/headers/sphere_lebedev_rule.h"
#include "jburkardt/headers/sphere_llq_grid.h"
#include "jburkardt/headers/sphere_llt_grid.h"
#include "jburkardt/headers/sphere_quad.h"
#include "jburkardt/headers/sphere_stereograph.h"
#include "jburkardt/headers/square_arbq_rule.h"
#include "jburkardt/headers/square_grid.h"
#include "jburkardt/headers/square_integrals.h"
#include "jburkardt/headers/square_symq_rule.h"
#include "jburkardt/headers/stochastic_diffusion.h"
#include "jburkardt/headers/stochastic_rk.h"
#include "jburkardt/headers/stokes_2d_exact.h"
#include "jburkardt/headers/subset_sum.h"
#include "jburkardt/headers/tet_mesh.h"
#include "jburkardt/headers/tetrahedron_arbq_rule.h"
#include "jburkardt/headers/tetrahedron_exactness.h"
#include "jburkardt/headers/tetrahedron_grid.h"
#include "jburkardt/headers/tetrahedron_integrals.h"
#include "jburkardt/headers/tetrahedron_ncc_rule.h"
#include "jburkardt/headers/simple_rkf45.h"
#include "jburkardt/headers/toeplitz_cholesky.h"
#include "jburkardt/headers/toms097.h"
#include "jburkardt/headers/toms179.h"
#include "jburkardt/headers/toms322.h"
#include "jburkardt/headers/toms443.h"
#include "jburkardt/headers/toms655.h"
#include "jburkardt/headers/toms743.h"
#include "jburkardt/headers/toms886.h"
#include "jburkardt/headers/treepack.h"
#include "jburkardt/headers/triangle_exactness.h"
#include "jburkardt/headers/triangle_fekete_rule.h"
#include "jburkardt/headers/triangle_felippa_rule.h"
#include "jburkardt/headers/triangle_grid.h"
#include "jburkardt/headers/triangle_integrals.h"
#include "jburkardt/headers/triangle_ncc_rule.h"
#include "jburkardt/headers/triangle_nco_rule.h"
#include "jburkardt/headers/triangle_symq_rule.h"
#include "jburkardt/headers/triangulation.h"
#include "jburkardt/headers/triangulation_node_to_element.h"
#include "jburkardt/headers/triangulation_svg.h"
#include "jburkardt/headers/truncated_normal.h"
#include "jburkardt/headers/unicycle.h"
#include "jburkardt/headers/upc.h"
#include "jburkardt/headers/values.h"
#include "jburkardt/headers/vandermonde.h"
#include "jburkardt/headers/vandermonde_approx_1d.h"
#include "jburkardt/headers/vandermonde_approx_2d.h"
#include "jburkardt/headers/vandermonde_interp_2d.h"
#include "jburkardt/headers/walsh.h"
#include "jburkardt/headers/wathen.h"
#include "jburkardt/headers/wavelet.h"
#include "jburkardt/headers/wedge_exactness.h"
#include "jburkardt/headers/wedge_felippa_rule.h"
#include "jburkardt/headers/weekday.h"
#include "jburkardt/headers/wishart.h"
#include "jburkardt/headers/wtime.h"
#include "jburkardt/headers/zero_rc.h"

// MathSuite ExprEval Wrappers
#include "ExprEval/exprcnts.h"
#include "ExprEval/exprilfs.h"
#include "ExprEval/exprmatrix.h"
#include "ExprEval/exprtrig.h"


__MATHSUITE  bool   FFT(ityp **, register int, const bool);
__MATHSUITE  bool   FFT2D(ityp **, register int[static 2], const bool);
__MATHSUITE void   harris(mpfr_t, const register mpfr_t, const register dim_typ, mpfr_t *restrict);
__MATHSUITE void   harris2(mpfr_t, const register mpfr_t, mpfr_t *restrict);
__MATHSUITE bool   t_harris(const register mpfr_t, const register mpfr_t, const register dim_typ, mpfr_t *restrict);
__MATHSUITE bool   t_harris2(const register mpfr_t, const register mpfr_t, mpfr_t *restrict);
__MATHSUITE void   eval(mpfr_t, mpfr_t *restrict, const register dim_typ, const mpfr_t);
__MATHSUITE void   deval(mpfr_t, mpfr_t *restrict, const register dim_typ, const mpfr_t);
__MATHSUITE short    _routhTable(mpfr_t **, const register dim_typ, fsel_typ *);
__MATHSUITE sel_typ    _juryTable(mpfr_t **, const register dim_typ);
__MATHSUITE  sel_typ    _matrixEigenValues(mpfr_t *restrict, mpfr_t *restrict, mpfr_t *restrict, const register dim_typ);
__MATHSUITE void    _matrixAdd(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixCAdd(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixQAdd(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixOAdd(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixSAdd(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixSub(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixCSub(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixQSub(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixOSub(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixSSub(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]);

__MATHSUITE  void   mmult_fast(const register dim_typ, const register dim_typ [static MAX_MATRICES], const register sel_typ, mpfr_t **, mpfr_t **, mpfr_t **,
	void (* const )(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]),
	void (* const )(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]),
	void (* const )(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]));

__MATHSUITE  void   squareOSMM(void (*)(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]), mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ, const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixMultiplication(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixCMultiplication(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixQMultiplication(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixOMultiplication(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixSMultiplication(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]);
__MATHSUITE void    _matrixKProduct(mpfr_t **, mpfr_t **, mpfr_t **, register dim_typ [static 2][2]);
__MATHSUITE void    _matrixKCProduct(mpfr_t **, mpfr_t **, mpfr_t **, register dim_typ [static 2][2]);
__MATHSUITE void    _matrixKQProduct(mpfr_t **, mpfr_t **, mpfr_t **, register dim_typ [static 2][2]);
__MATHSUITE void    _matrixKOProduct(mpfr_t **, mpfr_t **, mpfr_t **, register dim_typ [static 2][2]);
__MATHSUITE void    _matrixKSProduct(mpfr_t **, mpfr_t **, mpfr_t **, register dim_typ [static 2][2]);

/// lists_manager.c
__MATHSUITE sel_typ   checkItemTypeByExtension(const char [static MAX_EXTENSION_LENGTH]);
__MSSHELL_WRAPPER_ __MATHSUITE bool    listInsertProc(const char [static MAX_PATH_LENGTH], sel_typ);
__MATHSUITE dim_typ   searchItem(const char [static MAX_PATH_LENGTH], sel_typ);
__MATHSUITE nodelist *   listNo(dim_typ, sel_typ);
__MATHSUITE bool   listInsert(const char [static MAX_PATH_LENGTH], sel_typ);
__MATHSUITE bool   listDelete(dim_typ, sel_typ);
__MATHSUITE dim_typ   itemSelect(sel_typ);
__MATHSUITE void   refreshItem(dim_typ, sel_typ);
__MATHSUITE bool   saveItem(dim_typ, sel_typ);
__MATHSUITE void passToItem(dim_typ, sel_typ, const bool);
__MATHSUITE void   setCurItem(const dim_typ, sel_typ);
__MATHSUITE void   createItemData(dim_typ, sel_typ);
__MATHSUITE __WINCALL void   createItem(const char *, sel_typ);
__MATHSUITE void   viewItem(const dim_typ, sel_typ);
__MATHSUITE void   printListItem(const dim_typ, sel_typ);
__MATHSUITE void   updItem(sel_typ);
__MATHSUITE void   updAll(sel_typ);
__MATHSUITE void   delItem(const dim_typ, sel_typ);
__MATHSUITE dim_typ   getItemID(char *, sprog * const, sel_typ);
__MATHSUITE void   relItem(const dim_typ, bool);
__MATHSUITE void   renItem(const char *, const dim_typ, sel_typ);

/// logs_manager.c
__MSSHELL_WRAPPER_ void  __lmp_prog setLogBufLen(const char [static MAX_PATH_LENGTH], logObj * const, const size_t);


#ifdef __cplusplus
}
#endif // __cplusplus
