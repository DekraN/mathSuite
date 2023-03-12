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
    exprReturnerType rt; /* Return Handler or Returner */
    exprFuncType fptr; /* Function pointer */
    exprRealFuncType rfptr;
    int MIN, MAX; /* Min and MAX args for the function. */
    int refmin, refmax; /* Min and MAX ref. variables for the function */
    // int type; /* Function node type.  exprEvalNOde solves the function */
	mpCTRL * mp_ctrl; /* Control Node for Returning Type*/
    // struct _exprFunc *next; /* For linked list */
    };

/* Function list object */
struct _exprFuncList
    {
    map_t hashmap;
    // struct _exprFunc *head;
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
            exprReturnerType rt;
            exprFuncType fptr; /* Function pointer */
            exprRealFuncType rfptr; 
            struct _exprNode *nodes; /* Array of argument nodes */
            int nodecount; /* Number of argument nodes */
            EXPRTYPE **refs; /* Reference variables */
            int refcount; /* Number of variable references (not a reference counter) */
            // int type; /* Type of function for exprEvalNode if fptr is NULL */
            mpCTRL * mp_ctrl;
            } function;
        } data;
    };



/* Functions for function lists */
// EXPRERRTYPE exprFuncListAddType(exprFuncList *flist, char *name, EXPRERRTYPE (* function)(EXPRTYPE*, exprObj *, exprNode * nodes, EXPRTYPE []), const bool evalType, int type, int MIN, int MAX, int refmin, int refmax) 
EXPRERRTYPE exprFuncListGet(exprFuncList flist, char *name, exprReturnerType *rt, exprFuncType *ptr, exprRealFuncType *rptr, mpCTRL ** mp_ctrl, int *MIN, int *MAX, int *refmin, int *refmax);


#ifdef __cplusplus
}
#endif

#endif /* __BAVII_EXPRPRIV_H */

