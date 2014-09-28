// Project Default Settings Configuration
// WARNING!!! Assure to know what you're doing by managing this vars!!

// START


// Directives and Initializing some elements...
// Absolutely Don't touch this!!

// HEADER Guard
// #pragma once

// just indefine it if you want to compile this on UNIX Archs
#define WINOS

#define STACKALLOC

#define DEFAULT_AUTHOR "Wesker"
#define DEFAULT_SECNDAUTHOR "DekraN"

#define DEFAULT_MATRIXGET_COMMAND "get"
#define DEFAULT_MATRIXSET_COMMAND "set"
#define DEFAULT_MATRIXBACK_COMMAND "back"

#define DEFAULT_SYSTEM_LOG "syslog.LOG"

#define DEFAULT_DESCRIPTIONS_FOLDER "descriptions/"

#define DEFAULT_HELP_FILE_EXTENSION "txt"

#define DEFAULT_SCRIPTFILES_EXTENSION "mss"
#define DEFAULT_VARLIST_FILE_EXTENSION "vlf"
#define DEFAULT_MATRIX_FILE_EXTENSION "mf"
#define DEFAULT_LOG_FILE_EXTENSION "LOG"
#define DEFAULT_SYSLOG_FILE_EXTENSION DEFAULT_LOG_FILE_EXTENSION
#define DEFAULT_LAYOUT_FILE_EXTENSION "xml"
#define DEFAULT_COLORS_FILE_EXTENSION DEFAULT_LAYOUT_FILE_EXTENSION
#define DEFAULT_PATHLIST_FILE_EXTENSION "lf"
#define DEFAULT_MAIN_PATH "autorun"
#define DEFAULT_ENVS_ANSVALNAME "ANS"
#define DEFAULT_MAIN_ITEM 0

#define TEST_MATRIX_NAME "tmat"

#define DEFAULT_EXIT_MESSAGE "\n\n\nThank you for using this program.\n"


#define DEFAULT_RANDOM_SEED (unsigned)time(NULL)

#ifdef WINOS
	#define DEFAULT_COLOR 58 // COLOR_SMARTWHITE
#else
	#define DEFAULT_COLOR 10
#endif

#define DEFAULT_PRECISION 2
#define SHOWTIME_PRECISION 4
#define DEFAULT_STABILIZER_FACTOR 10
#define DEFAULT_ALGEBRA 0
#define DEFAULT_BASE 10

#define DEFAULT_MAX_SIMPLEXMETHOD_ITERATIONS 500 // 3000
#define DEFAULT_MAX_DSVD_ITERATIONS 30000
#define DEFAULT_MAX_EIGVALUES_ITERATIONS 50000

#define DEFAULT_GLOBALVAL 27.200

#define DEFAULT_MAX_FIBONACCI_MEMOIZABLE_INDEX 25
#define DEFAULT_MAX_FATTORIALE_MEMOIZABLE_INDEX 40
#define DEFAULT_MAX_SFATTORIALE_MEMOIZABLE_INDEX 80


#define DEFAULT_MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX DEFAULT_MAX_SFATTORIALE_MEMOIZABLE_INDEX
#define DEFAULT_MAX_ODD_SFATTORIALE_MEMOIZABLE_INDEX DEFAULT_MAX_EVEN_SFATTORIALE_MEMOIZABLE_INDEX-1


#define DEFAULT_MIN_OUTLIER_CONSTANT 1.00
#define DEFAULT_MAX_OUTLIER_CONSTANT 10.00

#define DEFAULT_OUTLIER_CONSTANT 1.50

#define DEFAULT_MIN_STIRLING_NUMBER DEFAULT_MAX_FATTORIALE_MEMOIZABLE_INDEX

#define DEFAULT_BLOCKSIZE 41
#define DEFAULT_MINOSMMDIM 2
#define DEFAULT_MINSTRASSENDIM 16

#define DEFAULT_MAX_ALIAS 5

#define DEFAULT_MAX_STABFACT 25

// DEFINIZIONE MACRO VARIABILI BOOLEANE
#define DEFAULT_SHOWDESCRIPTION false
#define DEFAULT_AUTOSETCURITEM true
#define DEFAULT_ITEMSSELECTBYPATH false
#define DEFAULT_ITEMSAUTOSAVING true

#ifdef WINOS
    #define DEFAULT_SYSLOGSECURITYCHECK false
#endif // WINOS

#define DEFAULT_SYSTEMPARSING true
#define DEFAULT_ADVCALCPARSING true
#define DEFAULT_MATRIXPARSING true
#define DEFAULT_RESETLISTS false
#define DEFAULT_SAVERESULTS true
#define DEFAULT_SHOWVARLIST true
#define DEFAULT_SHOWDIFFTIME false
#define DEFAULT_SHOWDIFFTIMEADVCALC true
#define DEFAULT_SHOWDIFFTIMEALGOPS true
#define DEFAULT_SHOWEXECTIME false
#define DEFAULT_PRINTTIME false
#define DEFAULT_PRINTROWSLABELS true
#define DEFAULT_DOMAINCHECK true
#define DEFAULT_CHARINSERT false
#define DEFAULT_INSERTMODE true
#define DEFAULT_LAZYEXECUTION true
#define DEFAULT_INVERSEOPERATIONS false
#define DEFAULT_AUTOTURNBACK true
#define DEFAULT_DEGREESENTERING false
#define DEFAULT_PROGREPEATCHECK true
#define DEFAULT_STRASSENOPTIMIZATION true

#define EXIT_CHAR 'C'

/*
Cambia il tipo predefinito di elementi delle matrici.
ATTENZIONE: Ammessi tutti i tipi purch� siano formattabili tramite conversion format,
ma le modifiche laddove le variabili del suddetto tipo vengano coinvolte in formattazioni
per funzioni INPUT o OUTPUT, sono assolutamente richieste per il corretto funzionamento
del programma e devono essere fatte manualmente agendo in maniera diretta
sullo script in questione
*/
//


#define tipo_predefinito double
#define MIN_VAL -DBL_MAX
#define MAX_VAL DBL_MAX
#define tipo_random_seed time_t
#define TYPE_DOMAIN(x) (x < MIN_VAL || x > MAX_VAL)
#define DOMAIN_DEFAULT DOMAIN_DBL


#define DEFAULT_COLORS_PATH "./colors.xml"


#define MIN_PRECISION 0
#define MAX_PRECISION DBL_DIG


#define OUTPUT_CONVERSION_FORMAT "%.*f",access(curLayout)->precision
#define INPUT_CONVERSION_FORMAT "%lf"


// Maximum Definitions Dimensions
//
#define MAX_RIGHE 100
#define MAX_COLONNE 100

// Conditional Compilations STUFFS...
///
/// #define FREEZE_SETTINGSMANAGER
#define ALLOW_MSSMANAGER
#define ALLOW_VARLISTMANAGER
#define ALLOW_LOGSMANAGER
#define ALLOW_SYSLOGMANAGER
#define ALLOW_MATMANAGER
#define ALLOW_LAYOUTSMANAGER
#define ALLOW_LFSMANAGER

#define ALLOW_PRECEDIT
#define ALLOW_STABFACTEDIT
#define ALLOW_BLOCKSIZEEDIT
#define ALLOW_MINOSMMDIMEDIT
#define ALLOW_MINSTRASSENDIMEDIT
#define ALLOW_OUTLIERCONSTEDIT
#define ALLOW_MINSRNUMBEREDIT
#define ALLOW_ALGEBRAEDIT
#define ALLOW_OUTLIERCONSTEDIT
#define ALLOW_BOOLVARSEDIT
#define ALLOW_MMIEDIT
#define ALLOW_EXITCHAREDIT
#define ALLOW_COLSMANAGER
#define ALLOW_RANDOMSEEDEDIT
#define ALLOW_BUFFERSEDIT
#define ALLOW_SETDEFAULTS
///
// END

// Only for experts.
// START

// Default SubPrograms CmdNames

#define CMD_BCALC "calc"
#define CMD_ADVCALC "aclc"
#define CMD_LINEARALGEBRAOPERATIONS "laop"
#define CMD_CHANGESETTINGS "settings"
#define CMD_MSSMANAGER DEFAULT_SCRIPTFILES_EXTENSION
#define USAGE_BCALC "[VALUE] [STRING] [VALUE] or [[EXPR]]"
#define USAGE_ADVCALC NULL_CHAR
#define USAGE_LINEARALGEBRAOPERATIONS NULL_CHAR
#define USAGE_CHANGESETTINGS NULL_CHAR
#define USAGE_MSSMANAGER NULL_CHAR

#define CMD_SECONDGRADEQSOLVER "secgreqs"
#define CMD_COMPLEXADD "cadd"
#define CMD_COMPLEXMUL "cmul"
#define CMD_SIMPLEXMETHOD "sxmeth"
#define CMD_NEWTONDIFFTABLES "dftables"
#define CMD_LAGRANGEINTERPOLATION "intrpl"
#define CMD_FID "fint"
#define CMD_STRAIGHTLINEFITTING "slinef"
#define CMD_PARABOLICCURVEFITTING "pcurvef"
#define CMD_LINEARSYSTEMSSOLVER "meqslvr"
#define USAGE_SECONDGRADEQSOLVER "[1X3MATRIX]"
#define USAGE_COMPLEXADD "[2X2MATRIX]"
#define USAGE_COMPLEXMUL "[2X2MATRIX]"
#define USAGE_SIMPLEXMETHOD "[nXmMATRIX] [n-1BOOLMATRIX]"
#define USAGE_NEWTONDIFFTABLES "[EXPR] [[EXPR]-(EXPR EXPR)COUPLES]"
#define USAGE_LAGRANGEINTERPOLATION "[EXPR] [2X[EXPR]MATRIX] [EXPR]"
#define USAGE_FID "[STRING] [EXPR] [EXPR] [EXPR] EXPR"
#define USAGE_STRAIGHTLINEFITTING "[EXPR] [2X(EXPR)MATRIX]"
#define USAGE_PARABOLICCURVEFITTING "[EXPR] [2X(EXPR)MATRIX]"
#define USAGE_LINEARSYSTEMSSOLVER "[MATRIX]"

// Remov the comment tag below to Add FIDs (Functions Identifiers) to the program
/// #define ADDFIDCONSTANTS

#ifdef ALLOW_MSSMANAGER
    #define CMD_CMDLINE "cmd"
    #define CMD_EXECMSSFILES "exec"
    #define CMD_SHOWUSAGES "info"
    #define USAGE_CMDLINE NULL_CHAR
    #define USAGE_EXECMSSFILES "[STRING]"
    #define USAGE_SHOWUSAGES "[STRING]"
#endif // ALLOW_MSSMANAGER

// linear algebra operations cmd
#define CMD_MATRIXSORT "sort"
#define CMD_MATRIXEIGVALUES "evs"
#define CMD_MATRIXNORM "norm"
#define CMD_MATRIXDET "det"
#define CMD_MATRIXTRACE "trace"
#define CMD_MATRIXRANK "rank"
#define CMD_MATRIXSVD "svd"
#define CMD_MATRIXINV "inv"
#define CMD_MATRIXTRANSPOSE "tran"
#define CMD_MATRIXADD "add"
#define CMD_TENSORADD "tadd"
#define CMD_MATRIXMULTIPLICATION "mul"
#define CMD_MATRIXKPRODUCT "kprod"
// #define CMD_MATRIXKPRODUCT "kron"
#define CMD_MATRIXPOWER "pow"
#define CMD_MATRIXKPOWER "kpow"
#define CMD_MATRIXPERVECTOR "mpv"
#define CMD_DOTPRODUCT "dot"
#define CMD_PERSCALARMULTIPLICATION "psmul"
#define CMD_SCALARDIVISIONMATRIX "sdmat"
#define CMD_ILLCONDITIONCHECKING "ichk"
#define CMD_MATRIXFATTLU "flu"
#define USAGE_MATRIXSORT "[VALUE] [MATRIX]"
#define USAGE_MATRIXEIGVALUES "[MATRIX]"
#define USAGE_MATRIXNORM "[MATRIX]"
#define USAGE_MATRIXDET "[MATRIX]"
#define USAGE_MATRIXTRACE "[MATRIX]"
#define USAGE_MATRIXRANK "[MATRIX]"
#define USAGE_MATRIXSVD "[MATRIX]"
#define USAGE_MATRIXINV "[MATRIX]"
#define USAGE_MATRIXTRANSPOSE "[MATRIX]"
#define USAGE_MATRIXADD "[nXmMATRIX] [nXmMATRIX]"
#define USAGE_TENSORADD "[lx[nXmMATRIX]] [lx[nXmMATRIX]]"
#define USAGE_MATRIXMULTIPLICATION "[iXjMATRIX] [jXkMATRIX]"
#define USAGE_MATRIXKPRODUCT "[iXjMATRIX] [kXlMATRIX]"
#define USAGE_MATRIXPOWER "[MATRIX] [EXPR]"
#define USAGE_MATRIXKPOWER "[MATRIX] [EXPR]"
#define USAGE_MATRIXPERVECTOR "[nXmMATRIX] [mX1MATRIX]"
#define USAGE_DOTPRODUCT "[1XnMATRIX] [nX1MATRIX]"
#define USAGE_PERSCALARMULTIPLICATION "[MATRIX] [EXPR]"
#define USAGE_SCALARDIVISIONMATRIX "[MATRIX] [EXPR]"
#define USAGE_ILLCONDITIONCHECKING "[MATRIX]"
#define USAGE_MATRIXFATTLU "[MATRIX]"

#ifdef ALLOW_PRECEDIT
    #define CMD_CHANGEPRECISION "prec"
    #define USAGE_CHANGEPRECISION "[EXPR]"
#endif
#ifdef ALLOW_STABFACTEDIT
    #define CMD_CHANGESTABFACT "stabfact"
    #define USAGE_CHANGESTABFACT "[EXPR]"
#endif
#ifdef ALLOW_BLOCKSIZEEDIT
	#define CMD_CHANGEBLOCKSIZE "bsiz"
	#define USAGE_CHANGEBLOCKSIZE "[EXPR]"
#endif
#ifdef ALLOW_MINOSMMDIMEDIT
	#define CMD_CHANGEMINOSMMDIM "mosmm"
	#define USAGE_CHANGEMINOSMMDIM "[EXPR]"
#endif
#ifdef ALLOW_MINSTRASSENDIMEDIT
	#define CMD_CHANGEMINSTRASSENDIM "msd"
	#define USAGE_CHANGEMINSTRASSENDIM "[EXPR]"
#endif
#ifdef ALLOW_MINSRNUMBEREDIT
	#define CMD_CHANGEMINSRNUMBER "msrn"
	#define USAGE_CHANGEMINSRNUMBER "[EXPR]"
#endif
#ifdef ALLOW_ALGEBRAEDIT
    #define CMD_CHANGEALGEBRA "alg"
    #define USAGE_CHANGEALGEBRA "[EXPR]"
#endif
#ifdef ALLOW_OUTLIERCONSTEDIT
	#define CMD_CHANGEOUTLIERCONST "oc"
	#define USAGE_CHANGEOUTLIERCONST "[EXPR]"
#endif
#ifdef ALLOW_EXITCHAREDIT
    #define CMD_CHANGEEXITCHAR "echar"
    #define USAGE_CHANGEEXITCHAR "[CHAR]"
#endif
#ifdef ALLOW_RANDOMSEEDEDIT
    #define CMD_CHANGERANDOMSEED "rseed"
    #define USAGE_CHANGERANDOMSEED "[EXPR]"
#endif
#ifdef ALLOW_BOOLVARSEDIT
    #define CMD_CHANGEBOOLVALUES "bools"
    #define USAGE_CHANGEBOOLVALUES "[EXPR]"
#endif
#ifdef ALLOW_MMIEDIT
    #define CMD_CHANGEMAXMEMIDX "mmidx"
    #define CMD_EMPTYMEMOIZERSBUFFERS "embuf"
    #define USAGE_CHANGEMAXMEMIDX "[EXPR] [EXPR] [EXPR] [EXPR]"
    #define USAGE_EMPTYMEMOIZERSBUFFERS "[VALUE]"
#endif
#ifdef ALLOW_BUFFERSEDIT
    #define CMD_EMPTYBUFFERS "ebufs"
    #define USAGE_EMPTYBUFFERS "[VALUE]"
#endif
#ifdef ALLOW_SETDEFAULTS
    #define CMD_SETDEFAULTS "def"
    #define USAGE_SETDEFAULTS NULL_CHAR
#endif

#ifdef ALLOW_VARLISTMANAGER
    #define CMD_ENVSMANAGER "envs"
    #define USAGE_ENVSMANAGER NULL_CHAR
    #define CMD_SETCURENV "cenv"
    #define CMD_OPENENV "eopen"
    #define CMD_CREATEENV "ecreate"
    #define CMD_READENV "eread"
    #define CMD_PRTENV "eprt"
    #define CMD_UPDENV "eupd"
    #define CMD_UPDALLENVS "eupdall"
    #define CMD_SAVEENV "esave"
    #define CMD_DELENV "edel"
    #define CMD_DELENV2 "edel2"
    #define CMD_DELALLENVS "edelall"
    #define CMD_DELALLENVS2 "edelall2"
    #define CMD_RESETENV "erel"
    #define CMD_RENENV "eren"
    #define CMD_SYNCENV "sync"
    #define USAGE_SETCURENV "[VALUE]"
    #define USAGE_OPENENV "[STRING]"
    #define USAGE_CREATEENV "[STRING]"
    #define USAGE_READENV "[VALUE]"
    #define USAGE_PRTENV "[VALUE]"
    #define USAGE_UPDENV NULL_CHAR
    #define USAGE_UPDALLENVS NULL_CHAR
    #define USAGE_SAVEENV "[STRING]"
    #define USAGE_DELENV "[VALUE]"
    #define USAGE_DELENV2 "[VALUE]"
    #define USAGE_DELALLENVS NULL_CHAR
    #define USAGE_DELALLENVS2 NULL_CHAR
    #define USAGE_RESETENV "[VALUE]"
    #define USAGE_RENENV "[VALUE] [STRING]"
    #define USAGE_SYNCENV "[VALUE]"
#endif // modifica_ambienti_variabili

#ifdef ALLOW_MATMANAGER
    #define CMD_MATRIXMANAGER "mat"
    #define USAGE_MATRIXMANAGER NULL_CHAR
    #define CMD_SETCURMAT "cmat"
    #define CMD_OPENMAT "mopen"
    #define CMD_CREATEMAT "mcreate"
    #define CMD_READMAT "mread"
    #define CMD_PRTMAT "mprt"
    #define CMD_UPDMAT "mupd"
    #define CMD_UPDALLMAT "mupdall"
    #define CMD_SAVEMAT "msave"
    #define CMD_DELMAT "mdel"
    #define CMD_DELMAT2 "mdel2"
    #define CMD_DELALLMATS "mdelall"
    #define CMD_DELALLMATS2 "mdelall2"
    #define CMD_RESETMAT "mrel"
    #define CMD_RENMAT "mren"
    #define CMD_EDITMAT "medit"
    #define USAGE_SETCURMAT "[VALUE]"
    #define USAGE_OPENMAT "[STRING]"
    #define USAGE_CREATEMAT "[STRING]"
    #define USAGE_READMAT "[VALUE]"
    #define USAGE_PRTMAT "[VALUE]"
    #define USAGE_UPDMAT NULL_CHAR
    #define USAGE_UPDALLMAT NULL_CHAR
    #define USAGE_SAVEMAT "[STRING]"
    #define USAGE_DELMAT "[VALUE]"
    #define USAGE_DELMAT2 "[VALUE]"
    #define USAGE_DELALLMATS NULL_CHAR
    #define USAGE_DELALLMATS2 NULL_CHAR
    #define USAGE_RESETMAT "[VALUE]"
    #define USAGE_RENMAT "[VALUE] [STRING]"
    #define USAGE_EDITMAT "[VALUE] [MATRIX]"
#endif // modifica_matrici

#ifdef ALLOW_LOGSMANAGER
    #define CMD_LOGSMANAGER "logs"
    #define USAGE_LOGSMANAGER NULL_CHAR
    #define CMD_SETCURLOG "clog"
    #define CMD_OPENLOG "lopen"
    #define CMD_CREATELOG "lcreate"
    #define CMD_READLOG "lread"
    #define CMD_PRTLOG "lprt"
    #define CMD_UPDLOG "lupd"
    #define CMD_UPDALLLOGS "lupdall"
    #define CMD_SAVELOG "lsave"
    #define CMD_DELLOG "ldel"
    #define CMD_DELLOG2 "ldel2"
    #define CMD_DELALLLOGS "ldelall"
    #define CMD_DELALLLOGS2 "ldelall2"
    #define CMD_RESETLOG "lrel"
    #define CMD_RENLOG "lren"
    #define CMD_EDITLOG "ledit"
    #define CMD_SETLOGBUFLEN "slblen"
    #define CMD_EMPTYLOGBUF "elbuf"
    #define USAGE_SETCURLOG "[VALUE]"
    #define USAGE_OPENLOG "[STRING]"
    #define USAGE_CREATELOG "[STRING]"
    #define USAGE_READLOG "[VALUE]"
    #define USAGE_PRTLOG "[VALUE]"
    #define USAGE_UPDLOG NULL_CHAR
    #define USAGE_UPDALLLOGS NULL_CHAR
    #define USAGE_SAVELOG "[STRING]"
    #define USAGE_DELLOG "[VALUE]"
    #define USAGE_DELLOG2 "[VALUE]"
    #define USAGE_DELALLLOGS NULL_CHAR
    #define USAGE_DELALLLOGS2 NULL_CHAR
    #define USAGE_RESETLOG "[VALUE]"
    #define USAGE_RENLOG "[VALUE] [STRING]"
    #define USAGE_EDITLOG "[VALUE]"
    #define USAGE_SETLOGBUFLEN "[VALUE] [EXPR]"
    #define USAGE_EMPTYLOGBUF "[VALUE]"
#endif // modifica_log

#ifdef ALLOW_SYSLOGMANAGER
    #define CMD_SYSLOGMANAGER "slog"
    #define USAGE_SYSLOGMANAGER NULL_CHAR
    #define CMD_OPENSYSLOG "sopen"
    #define CMD_CREATESYSLOG "screate"
    #define CMD_READSYSLOG "sread"
    #define CMD_PRTSYSLOG "sprt"
    #define CMD_EDITSYSLOG "sedit"
    #define CMD_UPDSYSLOG "supd"
    #define CMD_SAVESYSLOG "ssave"
    #define CMD_RESETSYSLOG "srel"
    #define CMD_SETSYSLOGBUFLEN "sslblen"
    #define CMD_EMPTYSYSLOGBUF "eslbuf"
    #define CMD_RENSYSLOG "sren"
    #define USAGE_OPENSYSLOG "[STRING]"
    #define USAGE_CREATESYSLOG "[STRING]"
    #define USAGE_READSYSLOG NULL_CHAR
    #define USAGE_PRTSYSLOG NULL_CHAR
    #define USAGE_EDITSYSLOG NULL_CHAR
    #define USAGE_UPDSYSLOG NULL_CHAR
    #define USAGE_SAVESYSLOG "[STRING]"
    #define USAGE_RESETSYSLOG NULL_CHAR
    #define USAGE_SETSYSLOGBUFLEN "[EXPR]"
    #define USAGE_EMPTYSYSLOGBUF NULL_CHAR
    #define USAGE_RENSYSLOG "[STRING]"
#endif // modifica_syslog

#ifdef ALLOW_LAYOUTSMANAGER
    #define CMD_LAYOUTSMANAGER "layts"
    #define USAGE_LAYOUTSMANAGER NULL_CHAR
    #define CMD_SETCURLAYOUT "clyt"
    #define CMD_OPENLAYOUT "lyopen"
    #define CMD_CREATELAYOUT "lycreate"
    #define CMD_READLAYOUT "lyread"
    #define CMD_PRTLAYOUT "lyprt"
    #define CMD_UPDLAYOUT "lyupd"
    #define CMD_UPDALLLAYOUTS "lyupdall"
    #define CMD_SAVELAYOUT "lysave"
    #define CMD_DELLAYOUT "lydel"
    #define CMD_DELLAYOUT2 "lydel2"
    #define CMD_DELALLLAYOUTS "lydelall"
    #define CMD_DELALLLAYOUTS2 "lydelall2"
    #define CMD_RESETLAYOUT "lyrel"
    #define CMD_RENLAYOUT "lyren"
    #define USAGE_SETCURLAYOUT "[VALUE]"
    #define USAGE_OPENLAYOUT "[STRING]"
    #define USAGE_CREATELAYOUT "[STRING]"
    #define USAGE_READLAYOUT "[VALUE]"
    #define USAGE_PRTLAYOUT "[VALUE]"
    #define USAGE_UPDLAYOUT NULL_CHAR
    #define USAGE_UPDALLLAYOUTS NULL_CHAR
    #define USAGE_SAVELAYOUT "[STRING]"
    #define USAGE_DELLAYOUT "[VALUE]"
    #define USAGE_DELLAYOUT2 "[VALUE]"
    #define USAGE_DELALLLAYOUTS NULL_CHAR
    #define USAGE_DELALLLAYOUTS2 NULL_CHAR
    #define USAGE_RESETLAYOUT "[VALUE]"
    #define USAGE_RENLAYOUT "[VALUE] [STRING]"
#endif // modifica_layouts

#ifdef ALLOW_COLSMANAGER
    #define CMD_COLORSMANAGER "cols"
    #define USAGE_COLORSMANAGER NULL_CHAR
    #define CMD_CHANGECOLORS "ccol"
    #define CMD_COLORSFILESLOADER "cload"
    #define CMD_BACKUPCOLORSFILES "cupd"
    #define USAGE_CHANGECOLORS "[VALUE] [VALUE]"
    #define USAGE_COLORSFILESLOADER "[STRING]"
    #define USAGE_BACKUPCOLORSFILES NULL_CHAR
#endif // modifica_colori

#ifdef ALLOW_LFSMANAGER
    #define CMD_LFSMANAGER DEFAULT_PATHLIST_FILE_EXTENSION
    #define USAGE_LFSMANAGER NULL_CHAR
    #define CMD_lfLoader "msload"
    #define CMD_LF_CREATE "mscreate"
    #define USAGE_lfLoader "[STRING]"
    #define USAGE_LFCREATE "[STRING]"
#endif

// Default Operations Identifiers

// Si raccomanda di racchiudere il simbolo tra apici, come esempi...

#define IDENTIFIER_ADDIZIONEBINARIA "badd"
#define IDENTIFIER_SOTTRAZIONEBINARIA "bsub"
#define IDENTIFIER_COMPLEMENT "comp"
#define IDENTIFIER_ELEVAMENTOAPOTENZA "pow"

#define IDENTIFIER_EXPANDEXPC "exp"
#define IDENTIFIER_EXP10ANDEXP10C "exp10"
#define IDENTIFIER_EXP2ANDEXP2C "exp2"
#define IDENTIFIER_RADICENESIMA "root"
#define IDENTIFIER_LOGARITMON "ln"
#define IDENTIFIER_LOGARITMO "log"
#define IDENTIFIER_LOGARITMO2 "log2"
#define IDENTIFIER_LOGARITMOBN "logn"
#define IDENTIFIER_LOGARITMOC "logc"
#define IDENTIFIER_LOGARITMO2C "log2c"
#define IDENTIFIER_LOGARITMO1P "log1p"
#define IDENTIFIER_LOGARITMO1PC "log1pc"
#define IDENTIFIER_LOGARITMO101P "log101p"
#define IDENTIFIER_LOGARITMO101PC "log101pc"
#define IDENTIFIER_LOGARITMO21P "log21p"
#define IDENTIFIER_LOGARITMO21PC "log21pc"


#define IDENTIFIER_CEXP "cexp"
#define IDENTIFIER_CEXP10 "cexp10"
#define IDENTIFIER_CEXP2 "cexp2"
#define IDENTIFIER_CPOW "cpow"
#define IDENTIFIER_CROOT "croot"
#define IDENTIFIER_CLOGN "clogn"
#define IDENTIFIER_CLN "cln"
#define IDENTIFIER_CLOG "clog"
#define IDENTIFIER_CLOG10 "clog10"
#define IDENTIFIER_CLOG2 "clog2"
#define IDENTIFIER_CLOG1P "clog1p"
#define IDENTIFIER_CLOG101P "clog101p"
#define IDENTIFIER_CLOG21P "clog21p"
#define IDENTIFIER_CARG "carg"
#define IDENTIFIER_CABS "cabs"
#define IDENTIFIER_QABS "qabs"
#define IDENTIFIER_OABS "oabs"
#define IDENTIFIER_SABS "sabs"


#define IDENTIFIER_BITCOUNTER "bcnt"
#define IDENTIFIER_UBITCOUNTER "ubcnt"
#define IDENTIFIER_VERSION "vers"
#define IDENTIFIER_PREC "prec"
#define IDENTIFIER_EXITCHAR "ec"
#define IDENTIFIER_STABFACT "sf"
#define IDENTIFIER_BLOCKSIZE "bsiz"
#define IDENTIFIER_MINOSMMDIM "mosmm"
#define IDENTIFIER_MINSTRASSENDIM "mstr"
#define IDENTIFIER_MINSRNUMBER "msrn"
#define IDENTIFIER_ALGEBRA "alg"
#define IDENTIFIER_OUTLIERCONST "oc"
#define IDENTIFIER_RSEED "rseed"
#define IDENTIFIER_MMIFIBO "mmi_fibo"
#define IDENTIFIER_MMIFACT "mmi_fact"
#define IDENTIFIER_MMIEVENSFACT "mmi_esfact"
#define IDENTIFIER_MMIODDSFACT "mmi_osfact"
#define IDENTIFIER_TRASFORMAANGOLI "angles"
#define IDENTIFIER_SINANDSINH "sin"
#define IDENTIFIER_COSANDCOSH "cos"
#define IDENTIFIER_TANANDTANH "tan"
#define IDENTIFIER_CSCANDCSCH "csc"
#define IDENTIFIER_SECANDSECH "sec"
#define IDENTIFIER_COTANDCOTH "cot"
#define IDENTIFIER_ASINANDASINH "asin"
#define IDENTIFIER_ACOSANDACOSH "acos"
#define IDENTIFIER_ATANANDATANH "atan"
#define IDENTIFIER_ATAN2 "atan2"
#define IDENTIFIER_ACSCANDACSCH "acsc"
#define IDENTIFIER_ASECANDASECH "asec"
#define IDENTIFIER_ACOTANDACOTH "acot"
#define IDENTIFIER_HSINANDHSINH "hsin"
#define IDENTIFIER_QSINANDQSINH "qsin"
#define IDENTIFIER_HCOSANDHCOSH "hcos"
#define IDENTIFIER_QCOSANDQCOSH "qcos"
#define IDENTIFIER_HSECANDHSECH "hsec"
#define IDENTIFIER_QSECANDQSECH "qsec"
#define IDENTIFIER_HCSCANDHCSCH "hcsc"
#define IDENTIFIER_QCSCANDQCSCH "qcsc"
#define IDENTIFIER_HTANANDHTANH "htan"
#define IDENTIFIER_QTANANDQTANH "qtan"
#define IDENTIFIER_HCOTANDHCOTH "hcot"
#define IDENTIFIER_QCOTANDQCOTH "qcot"
#define IDENTIFIER_VSINANDVSINH "vsin"
#define IDENTIFIER_CVSINANDCVSINH "cvsin"
#define IDENTIFIER_VCOSANDVCOSH "vcos"
#define IDENTIFIER_CVCOSANDCVCOSH "cvcos"
#define IDENTIFIER_HVSINANDHVSINH "hvsin"
#define IDENTIFIER_HCVSINANDHCVSINH "hcvsin"
#define IDENTIFIER_HVCOSANDHVCOSH "hvcos"
#define IDENTIFIER_HCVCOSANDHCVCOSH "hcvcos"
#define IDENTIFIER_QVSINANDQVSINH "qvsin"
#define IDENTIFIER_QCVSINANDQCVSINH "qcvsin"
#define IDENTIFIER_QVCOSANDQVCOSH "qvcos"
#define IDENTIFIER_QCVCOSANDQCVCOSH "qcvcos"
#define IDENTIFIER_ESECANDESECH "esec"
#define IDENTIFIER_ECSCANDECSCH "ecsc"
#define IDENTIFIER_HESECANDHESECH "hesec"
#define IDENTIFIER_HECSCANDHECSCH "hecsc"
#define IDENTIFIER_QESECANDQESECH "qesec"
#define IDENTIFIER_QECSCANDQECSCH "qecsc"
#define IDENTIFIER_SINCANDSINCH "sinc"
#define IDENTIFIER_HSINCANDHSINCH "hsinc"
#define IDENTIFIER_QSINCANDQSINCH "qsinc"
#define IDENTIFIER_COSCANDCOSCH "cosc"
#define IDENTIFIER_HCOSCANDHCOSCH "hcosc"
#define IDENTIFIER_QCOSCANDQCOSCH "qcosc"
#define IDENTIFIER_SECCANDSECCH "secc"
#define IDENTIFIER_HSECCANDHSECCH "hsecc"
#define IDENTIFIER_QSECCANDQSECCH "qsecc"
#define IDENTIFIER_CSCCANDCSCCH "cscc"
#define IDENTIFIER_HCSCCANDHCSCCH "hcscc"
#define IDENTIFIER_QCSCCANDQCSCCH "qcscc"
#define IDENTIFIER_TANCANDTANCH "tanc"
#define IDENTIFIER_HTANCANDHTANCH "htanc"
#define IDENTIFIER_QTANCANDQTANCH "qtanc"
#define IDENTIFIER_COTCANDCOTCH "cotc"
#define IDENTIFIER_HCOTCANDHCOTCH "hcotc"
#define IDENTIFIER_QCOTCANDQCOTCH "qcotc"

#define IDENTIFIER_CSIN "csin"
#define IDENTIFIER_CCOS "ccos"
#define IDENTIFIER_CTAN "ctan"
#define IDENTIFIER_CCSC "ccsc"
#define IDENTIFIER_CSEC "csec"
#define IDENTIFIER_CCOT "ccot"
#define IDENTIFIER_CASIN "casin"
#define IDENTIFIER_CACOS "cacos"
#define IDENTIFIER_CATAN "catan"
#define IDENTIFIER_CACSC "cacsc"
#define IDENTIFIER_CASEC "casec"
#define IDENTIFIER_CACOT "cacot"
#define IDENTIFIER_CHSIN "chsin"
#define IDENTIFIER_CQSIN "cqsin"
#define IDENTIFIER_CHCOS "chcos"
#define IDENTIFIER_CQCOS "cqcos"
#define IDENTIFIER_CHSEC "chsec"
#define IDENTIFIER_CQSEC "cqsec"
#define IDENTIFIER_CHCSC "chcsc"
#define IDENTIFIER_CQCSC "cqcsc"
#define IDENTIFIER_CHTAN "chtan"
#define IDENTIFIER_CQTAN "cqtan"
#define IDENTIFIER_CHCOT "chcot"
#define IDENTIFIER_CQCOT "cqcot"
#define IDENTIFIER_CVSIN "cpxvsin"
#define IDENTIFIER_CCVSIN "ccvsin"
#define IDENTIFIER_CVCOS "cpxvcos"
#define IDENTIFIER_CCVCOS "ccvcos"
#define IDENTIFIER_CHVSIN "chvsin"
#define IDENTIFIER_CHCVSIN "chcvsin"
#define IDENTIFIER_CHVCOS "chvcos"
#define IDENTIFIER_CHCVCOS "chcvcos"
#define IDENTIFIER_CQVSIN "cqvsin"
#define IDENTIFIER_CQCVSIN "cqcvsin"
#define IDENTIFIER_CQVCOS "cqvcos"
#define IDENTIFIER_CQCVCOS "cqcvcos"
#define IDENTIFIER_CESEC "cesec"
#define IDENTIFIER_CECSC "cecsc"
#define IDENTIFIER_CHESEC "chesec"
#define IDENTIFIER_CHECSC "checsc"
#define IDENTIFIER_CQESEC "cqesec"
#define IDENTIFIER_CQECSC "cqecsc"
#define IDENTIFIER_CSINC "csinc"
#define IDENTIFIER_CHSINC "chsinc"
#define IDENTIFIER_CQSINC "cqsinc"
#define IDENTIFIER_CCOSC "ccosc"
#define IDENTIFIER_CHCOSC "chcosc"
#define IDENTIFIER_CQCOSC "cqcosc"
#define IDENTIFIER_CSECC "csecc"
#define IDENTIFIER_CHSECC "chsecc"
#define IDENTIFIER_CQSECC "cqsecc"
#define IDENTIFIER_CCSCC "ccscc"
#define IDENTIFIER_CHCSCC "chcscc"
#define IDENTIFIER_CQCSCC "cqcscc"
#define IDENTIFIER_CTANC "ctanc"
#define IDENTIFIER_CHTANC "chtanc"
#define IDENTIFIER_CQTANC "cqtanc"
#define IDENTIFIER_CCOTC "ccotc"
#define IDENTIFIER_CHCOTC "chcotc"
#define IDENTIFIER_CQCOTC "cqcotc"

#define IDENTIFIER_MATRIXNORM "norm"
#define IDENTIFIER_MATRIXTRACE "trace"
#define IDENTIFIER_MATRIXDET "det"
#define IDENTIFIER_MATRIXRANK "rank"
#define IDENTIFIER_MATRIXILLCHK "illchk"
#define IDENTIFIER_SCALARPROD "sprod"
#define IDENTIFIER_MCD "GCD"
#define IDENTIFIER_MCM "lcm"
#define IDENTIFIER_APPROSSIMAZIONE "approx"
#define IDENTIFIER_RISOLUTOREEQUAZIONISECONDOGRADO "sgeqs"
#define IDENTIFIER_COMPLEXADD "cadd"
#define IDENTIFIER_COMPLEXMUL "cmul"
#define IDENTIFIER_QUATERNIONSADD "qadd"
#define IDENTIFIER_QUATERNIONSMUL "qmul"
#define IDENTIFIER_SOMMASUCCESSIONEGEOMETRICA "gsum"
#define IDENTIFIER_SOMMASUCCESSIONEARMONICA "asum"
#define IDENTIFIER_SOMMASUCCESSIONEARMONICAGEN "gasum"
#define IDENTIFIER_SOMMASUCCESSIONEFIBONACCI "fsum"
#define IDENTIFIER_SOMMASUCCESSIONEFATTORIALE "fasum"
#define IDENTIFIER_SOMMASUCCESSIONESEMIFATTORIALE "sfasum"
#define IDENTIFIER_SOMMATORIA "sum"
#define IDENTIFIER_PRODUTTORIA "product"
#define IDENTIFIER_MODE "mode"
#define IDENTIFIER_MEDIA "media"
#define IDENTIFIER_VARIANCE "var"
#define IDENTIFIER_COVARIANCE "cov"
#define IDENTIFIER_STDDEV "stddev"
#define IDENTIFIER_OUTLIER "otlr"
#define IDENTIFIER_OUTLIER2 "otlr2"
#define IDENTIFIER_MAP "map"
#define IDENTIFIER_MEDIAGEOMETRICA "geomedia"
#define IDENTIFIER_MEDIAARMONICA "armedia"
#define IDENTIFIER_MEDIAPOTENZA "powmedia"
#define IDENTIFIER_VALORECENTRALE "cval"
#define IDENTIFIER_FIRSTQUARTILE "q1"
#define IDENTIFIER_MEDIANA "mediana"
#define IDENTIFIER_THIRDQUARTILE "q3"
#define IDENTIFIER_SOMMAPRIMINNUMERI "fnnsum"
#define IDENTIFIER_RADICEQUADRATA "sqrt"
#define IDENTIFIER_RADICECUBICA "cbrt"
#define IDENTIFIER_CSQRT "csqrt"
#define IDENTIFIER_CCBRT "ccbrt"
#define OPERATOR_FATTORIALE "!"
#define OPERATOR_SEMIFATTORIALE "!!"
#define OPERATOR_FIBONACCI "\""
#define IDENTIFIER_FATTORIALE "fact"
#define IDENTIFIER_SEMIFATTORIALE "sfact"
#define IDENTIFIER_STIRLING "stlrng"
#define IDENTIFIER_FIBONACCI "fibo"
#define IDENTIFIER_ASEMIFATTORIALE "asfact"
#define IDENTIFIER_AFIBONACCI "afibo"
#define IDENTIFIER_PERMUTATIONS "perm"
#define IDENTIFIER_PERMUTATIONSREP "permrep"
#define IDENTIFIER_KPERMUTATIONS "kperm"
#define IDENTIFIER_KPERMUTATIONSREP "kpermrep"
#define IDENTIFIER_COMBINATIONS "comb"
#define IDENTIFIER_COMBINATIONSREP "combrep"
#define IDENTIFIER_PASCALTRIANGLE "ptrn"
#define IDENTIFIER_NUMERIROMANI "romn"
#define IDENTIFIER_PRIMINNUMERIPRIMI "fpnum"
#define IDENTIFIER_NESIMONUMEROPRIMO "npnum"
#define OPERATOR_PRIMORIALE "#"
#define IDENTIFIER_PRIMORIALE "primr"
#define IDENTIFIER_SOMMAPRIMINNUMERIPRIMI "fpnsum"
#define IDENTIFIER_FIBONACCIALE "fibnc"
#define IDENTIFIER_CAMBIAMENTODIBASE "cbase"
#define IDENTIFIER_GENERATOREMATRICIRANDOM "randmat"
#define IDENTIFIER_NOTANOPERATION "info" // used To Show the BASECALC Informations
//

// END

// v1.5 14/05/2013
// for any problem contact me at: marco_chiarelli@yahoo.it
