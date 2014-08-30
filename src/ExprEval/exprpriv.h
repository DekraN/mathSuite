/*
    File: exprpriv.h
    Auth: Brian Allen Vanderburg II
    Date: Tuesday, February 28, 2006
    Desc: Private include file for ExprEval library

    This file is part of ExprEval.
*/


/* Include once */
#ifndef __BAVII_EXPRPRIV_H
#define __BAVII_EXPRPRIV_H

/* Need some definitions, NULL, etc */
#include <stddef.h>

/* Include config and main expreval header */
#include "expreval.h"
#include "exprconf.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Version number
*/
#define EXPR_VERSIONMAJOR 2
#define EXPR_VERSIONMINOR 7

/* Node types */
enum
    {
    EXPR_NODETYPE_UNKNOWN = 0,
    EXPR_NODETYPE_MULTI,
    EXPR_NODETYPE_ADD,
    EXPR_NODETYPE_SUBTRACT,
    EXPR_NODETYPE_MULTIPLY,
    EXPR_NODETYPE_DIVIDE,
    EXPR_NODETYPE_EXPONENT,
    EXPR_NODETYPE_NEGATE,
    EXPR_NODETYPE_VALUE,
    EXPR_NODETYPE_VARIABLE,
    EXPR_NODETYPE_ASSIGN,
    EXPR_NODETYPE_FUNCTION
    };

/* Functions can be evaluated directly in EXPREVAL.  If fptr
   is NULL, type is used to determine what the function is */
enum
    {
    EXPR_NODEFUNC_UNKNOWN = 0,
    EXPR_NODEFUNC_EXIT,
    EXPR_NODEFUNC_MSS,
    EXPR_NODEFUNC_ABS,
    EXPR_NODEFUNC_MOD,
    EXPR_NODEFUNC_IPART,
    EXPR_NODEFUNC_FPART,
    EXPR_NODEFUNC_MIN,
    EXPR_NODEFUNC_MAX,
    //
    EXPR_NODEFUNC_BITCOUNTER,

    EXPR_NODEFUNC_VERSION,
    EXPR_NODEFUNC_EXITCHAR,
    EXPR_NODEFUNC_PREC,
    EXPR_NODEFUNC_STABFACT,
   	EXPR_NODEFUNC_MINSRNUMBER,
    EXPR_NODEFUNC_ALGEBRA,
    EXPR_NODEFUNC_OUTLIERCONST,
    EXPR_NODEFUNC_RSEED,
    EXPR_NODEFUNC_MMIFIBO,
    EXPR_NODEFUNC_MMIFACT,
   	EXPR_NODEFUNC_MMIEVENSFACT,
    EXPR_NODEFUNC_MMIODDSFACT,
    //
    EXPR_NODEFUNC_BINSUM,
    EXPR_NODEFUNC_BINSUB,
    EXPR_NODEFUNC_COMP,
    //
    EXPR_NODEFUNC_POW,
    EXPR_NODEFUNC_SQRT,
    //
    EXPR_NODEFUNC_CBRT,
    EXPR_NODEFUNC_ROOT,
    //
    EXPR_NODEFUNC_SIN,
    EXPR_NODEFUNC_SINH,
    //
    EXPR_NODEFUNC_CSC,
    EXPR_NODEFUNC_CSCH,
    EXPR_NODEFUNC_ASIN,
    EXPR_NODEFUNC_ASINH,
    EXPR_NODEFUNC_ACSC,
    EXPR_NODEFUNC_ACSCH,
    EXPR_NODEFUNC_COS,
    EXPR_NODEFUNC_COSH,
    EXPR_NODEFUNC_SEC,
    EXPR_NODEFUNC_SECH,
    EXPR_NODEFUNC_ACOS,
    EXPR_NODEFUNC_ACOSH,
    EXPR_NODEFUNC_ASEC,
    EXPR_NODEFUNC_ASECH,
    EXPR_NODEFUNC_TAN,
    EXPR_NODEFUNC_TANH,
    EXPR_NODEFUNC_COT,
    EXPR_NODEFUNC_COTH,
    EXPR_NODEFUNC_ATAN,
    EXPR_NODEFUNC_ATANH,
    EXPR_NODEFUNC_ACOT,
    EXPR_NODEFUNC_ACOTH,
    EXPR_NODEFUNC_HSIN,
    EXPR_NODEFUNC_HSINH,
    EXPR_NODEFUNC_QSIN,
    EXPR_NODEFUNC_QSINH,
    EXPR_NODEFUNC_HCOS,
    EXPR_NODEFUNC_HCOSH,
    EXPR_NODEFUNC_QCOS,
    EXPR_NODEFUNC_QCOSH,
    EXPR_NODEFUNC_HSEC,
    EXPR_NODEFUNC_HSECH,
    EXPR_NODEFUNC_QSEC,
    EXPR_NODEFUNC_QSECH,
    EXPR_NODEFUNC_HCSC,
    EXPR_NODEFUNC_HCSCH,
    EXPR_NODEFUNC_QCSC,
    EXPR_NODEFUNC_QCSCH,
    EXPR_NODEFUNC_HTAN,
    EXPR_NODEFUNC_HTANH,
    EXPR_NODEFUNC_QTAN,
    EXPR_NODEFUNC_QTANH,
    EXPR_NODEFUNC_HCOT,
    EXPR_NODEFUNC_HCOTH,
    EXPR_NODEFUNC_QCOT,
    EXPR_NODEFUNC_QCOTH,
    EXPR_NODEFUNC_VSIN,
    EXPR_NODEFUNC_VSINH,
    EXPR_NODEFUNC_CVSIN,
    EXPR_NODEFUNC_CVSINH,
    EXPR_NODEFUNC_VCOS,
    EXPR_NODEFUNC_VCOSH,
    EXPR_NODEFUNC_CVCOS,
    EXPR_NODEFUNC_CVCOSH,
    EXPR_NODEFUNC_HVSIN,
    EXPR_NODEFUNC_HVSINH,
    EXPR_NODEFUNC_HCVSIN,
    EXPR_NODEFUNC_HCVSINH,
    EXPR_NODEFUNC_QVSIN,
    EXPR_NODEFUNC_QVSINH,
    EXPR_NODEFUNC_HHCVSIN,
    EXPR_NODEFUNC_HHCVSINH,
    EXPR_NODEFUNC_HVCOS,
    EXPR_NODEFUNC_HVCOSH,
    EXPR_NODEFUNC_HCVCOS,
    EXPR_NODEFUNC_HCVCOSH,
    EXPR_NODEFUNC_QVCOS,
    EXPR_NODEFUNC_QVCOSH,
    EXPR_NODEFUNC_QCVCOS,
    EXPR_NODEFUNC_QCVCOSH,
    EXPR_NODEFUNC_ESEC,
    EXPR_NODEFUNC_ESECH,
    EXPR_NODEFUNC_ECSC,
    EXPR_NODEFUNC_ECSCH,
    EXPR_NODEFUNC_HESEC,
    EXPR_NODEFUNC_HESECH,
    EXPR_NODEFUNC_HECSC,
    EXPR_NODEFUNC_HECSCH,
    EXPR_NODEFUNC_QESEC,
    EXPR_NODEFUNC_QESECH,
    EXPR_NODEFUNC_QECSC,
    EXPR_NODEFUNC_QECSCH,
    EXPR_NODEFUNC_SINC,
    EXPR_NODEFUNC_SINCH,
    EXPR_NODEFUNC_HSINC,
    EXPR_NODEFUNC_HSINCH,
    EXPR_NODEFUNC_QSINC,
    EXPR_NODEFUNC_QSINCH,
    EXPR_NODEFUNC_COSC,
    EXPR_NODEFUNC_COSCH,
    EXPR_NODEFUNC_HCOSC,
    EXPR_NODEFUNC_HCOSCH,
    EXPR_NODEFUNC_QCOSC,
    EXPR_NODEFUNC_QCOSCH,
    EXPR_NODEFUNC_SECC,
    EXPR_NODEFUNC_SECCH,
    EXPR_NODEFUNC_HSECC,
    EXPR_NODEFUNC_HSECCH,
    EXPR_NODEFUNC_QSECC,
    EXPR_NODEFUNC_QSECCH,
    EXPR_NODEFUNC_CSCC,
    EXPR_NODEFUNC_CSCCH,
    EXPR_NODEFUNC_HCSCC,
    EXPR_NODEFUNC_HCSCCH,
    EXPR_NODEFUNC_QCSCC,
    EXPR_NODEFUNC_QCSCCH,
    EXPR_NODEFUNC_TANC,
    EXPR_NODEFUNC_TANCH,
    EXPR_NODEFUNC_HTANC,
    EXPR_NODEFUNC_HTANCH,
    EXPR_NODEFUNC_QTANC,
    EXPR_NODEFUNC_QTANCH,
    EXPR_NODEFUNC_COTC,
    EXPR_NODEFUNC_COTCH,
    EXPR_NODEFUNC_HCOTC,
    EXPR_NODEFUNC_HCOTCH,
    EXPR_NODEFUNC_QCOTC,
    EXPR_NODEFUNC_QCOTCH,
    //
    EXPR_NODEFUNC_ATAN2,
    //
    EXPR_NODEFUNC_MATRIXDET,
    EXPR_NODEFUNC_MATRIXNORM,
    EXPR_NODEFUNC_MATRIXTRACE,
    EXPR_NODEFUNC_MATRIXRANK,
    EXPR_NODEFUNC_MATRIXILLCHK,
    EXPR_NODEFUNC_SCALARPROD,
    //
    EXPR_NODEFUNC_LOG,
    EXPR_NODEFUNC_LOG2,
    EXPR_NODEFUNC_POW10,
    EXPR_NODEFUNC_LN,
    EXPR_NODEFUNC_EXP,
    EXPR_NODEFUNC_EXPC,
    EXPR_NODEFUNC_EXP10,
    EXPR_NODEFUNC_EXP10C,
    EXPR_NODEFUNC_EXP2,
    EXPR_NODEFUNC_EXP2C,
    EXPR_NODEFUNC_LOGN,
    EXPR_NODEFUNC_LOGC,
    EXPR_NODEFUNC_LNC,
    EXPR_NODEFUNC_LOG2C,
    EXPR_NODEFUNC_LOG1P,
    EXPR_NODEFUNC_LOG1PC,
    EXPR_NODEFUNC_CEIL,
    EXPR_NODEFUNC_FLOOR,
    EXPR_NODEFUNC_SGEQSOLVER,
    EXPR_NODEFUNC_COMPLEXSUM,
    EXPR_NODEFUNC_COMPLEXPROD,
    EXPR_NODEFUNC_RAND,
    EXPR_NODEFUNC_RANDOM,
    EXPR_NODEFUNC_RANDOMIZE,
    EXPR_NODEFUNC_DEG,
    EXPR_NODEFUNC_RAD,
    EXPR_NODEFUNC_RECTTOPOLR,
    EXPR_NODEFUNC_RECTTOPOLA,
    EXPR_NODEFUNC_POLTORECTX,
    EXPR_NODEFUNC_POLTORECTY,
    EXPR_NODEFUNC_CBASE,
    EXPR_NODEFUNC_NPNUM,
    EXPR_NODEFUNC_PRIMORIAL,
    EXPR_NODEFUNC_FPNSUM,
    EXPR_NODEFUNC_FIBONACCIAL,
    EXPR_NODEFUNC_LCM,
    EXPR_NODEFUNC_GCD,
    EXPR_NODEFUNC_FACT,
    EXPR_NODEFUNC_SFACT,
    EXPR_NODEFUNC_STIRLING,
    EXPR_NODEFUNC_FIBO,
    EXPR_NODEFUNC_PERM,
    EXPR_NODEFUNC_COMB,
    EXPR_NODEFUNC_GSUM,
    EXPR_NODEFUNC_ASUM,
    EXPR_NODEFUNC_GASUM,
    EXPR_NODEFUNC_FSUM,
    EXPR_NODEFUNC_FASUM,
    EXPR_NODEFUNC_SFASUM,
    EXPR_NODEFUNC_FNNSUM,
    //
    EXPR_NODEFUNC_SUM,
    EXPR_NODEFUNC_PRODUCT,
    EXPR_NODEFUNC_MEDIA,
    EXPR_NODEFUNC_VARIANCE,
    EXPR_NODEFUNC_COVARIANCE,
    EXPR_NODEFUNC_STDDEV,
    EXPR_NODEFUNC_OUTLIER,
    EXPR_NODEFUNC_OUTLIER2,
    EXPR_NODEFUNC_MAP,
    EXPR_NODEFUNC_GEOMEDIA,
    EXPR_NODEFUNC_ARMEDIA,
    EXPR_NODEFUNC_POWMEDIA,
    EXPR_NODEFUNC_CVAL,
    EXPR_NODEFUNC_FIRSTQUARTILE,
    EXPR_NODEFUNC_MEDIANA,
    EXPR_NODEFUNC_THIRDQUARTILE,
    EXPR_NODEFUNC_DAY,
    //
    EXPR_NODEFUNC_IF,
    EXPR_NODEFUNC_SELECT,
    EXPR_NODEFUNC_EQUAL,
    EXPR_NODEFUNC_ABOVE,
    EXPR_NODEFUNC_BELOW,
    EXPR_NODEFUNC_AVG,
    EXPR_NODEFUNC_CLIP,
    EXPR_NODEFUNC_CLAMP,
    EXPR_NODEFUNC_PNTCHANGE,
    EXPR_NODEFUNC_POLY,
    EXPR_NODEFUNC_AND,
    EXPR_NODEFUNC_OR,
    EXPR_NODEFUNC_NOT,
    EXPR_NODEFUNC_FOR,
    EXPR_NODEFUNC_MANY
    };

/* Forward declarations */
typedef struct _exprFunc exprFunc;
typedef struct _exprVal exprVal;

/* Expression object */
struct _exprObj
    {
    struct _exprFuncList *flist; /* Functions */
    struct _exprValList *vlist; /* Variables */
    struct _exprValList *clist; /* Constants */
    struct _exprNode *headnode; /* Head parsed node */

    exprBreakFuncType breakerfunc; /* Break function type */

    void *userdata; /* User data, can be any 32 bit value */
    int parsedgood; /* non-zero if successfully parsed */
    int parsedbad; /* non-zero if parsed but unsuccessful */
    int breakcount; /* how often to check the breaker function */
    int breakcur; /* do we check the breaker function yet */
    int starterr; /* start position of an error */
    int enderr; /* end position of an error */
    };

/* Object for a function */
struct _exprFunc
    {
    char *fname; /* Name of the function */
    exprFuncType fptr; /* Function pointer */
    int min, max; /* Min and max args for the function. */
    int refmin, refmax; /* Min and max ref. variables for the function */
    int type; /* Function node type.  exprEvalNOde solves the function */

    struct _exprFunc *next; /* For linked list */
    };

/* Function list object */
struct _exprFuncList
    {
    struct _exprFunc *head;
    };

/* Object for values */
struct _exprVal
    {
    char *vname; /* Name of the value */
    EXPRTYPE vval; /* Value of the value */
    EXPRTYPE *vptr; /* Pointer to a value.  Used only if not NULL */

    struct _exprVal *next; /* For linked list */
    };

/* Value list */
struct _exprValList
    {
    struct _exprVal *head;
    };

/* Expression node type */
struct _exprNode
    {
    int type; /* Node type */

    union _data /* Union of info for various types */
        {
        struct _oper
            {
            struct _exprNode *nodes; /* Operation arguments */
            int nodecount; /* Number of arguments */
            } oper;

        struct _variable
            {
            EXPRTYPE *vaddr; /* Used if EXPR_FAST_VAR_ACCESS defined */
            } variable;

        struct _value
            {
            EXPRTYPE value; /* Value if type is value */
            } value;

        struct _assign /* Assignment struct */
            {
            EXPRTYPE *vaddr; /* Used if EXPR_FAST_VAR_ACCESS defined */
            struct _exprNode *node; /* Node to evaluate */
            } assign;

        struct _function
            {
            exprFuncType fptr; /* Function pointer */
            struct _exprNode *nodes; /* Array of argument nodes */
            int nodecount; /* Number of argument nodes */
            EXPRTYPE **refs; /* Reference variables */
            int refcount; /* Number of variable references (not a reference counter) */
            int type; /* Type of function for exprEvalNode if fptr is NULL */
            } function;
        } data;
    };



/* Functions for function lists */
int exprFuncListAddType(exprFuncList *flist, char *name, int type, int min, int max, int refmin, int refmax);
int exprFuncListGet(exprFuncList *flist, char *name, exprFuncType *ptr, int *type, int *min, int *max, int *refmin, int *refmax);


#ifdef __cplusplus
}
#endif

#endif /* __BAVII_EXPRPRIV_H */

