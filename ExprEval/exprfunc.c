/*
    File: exprfunc.c
    Auth: Brian Allen Vanderburg II
    Date: Thursday, April 24, 2003
    Desc: Expression function list routines

    This file is part of ExprEval.
*/



/* Includes */
#include "exprincl.h"

#include "exprpriv.h"
#include "exprmem.h"

/* Internal functions */
static exprFunc *exprCreateFunc(char *name, exprReturnerType rt, exprFuncType ptr, exprRealFuncType rptr, mpCTRL * mp_ctrl, int MIN, int MAX, int refmin, int refmax);
static void exprFuncListFreeData(exprFuncList flist);


/* This function creates the function list, */
void exprFuncListCreate(exprFuncList *flist)
    {
	(*flist) = hashmap_new();
    return;
    }

/* Add a function to the list */
EXPRERRTYPE exprFuncListAdd(exprFuncList flist, char *name, exprReturnerType rt, exprFuncType ptr, exprRealFuncType rptr, mpCTRL * mp_ctrl, int MIN, int MAX, int refmin, int refmax)
    {
    exprFunc *tmp;
    int result;

    if(flist == NULL)
        return EXPR_ERROR_NULLPOINTER;

    /* Make sure the name is valid */
    if(!exprValidIdent(name))
        return EXPR_ERROR_BADIDENTIFIER;

    /* Fix values only if none are negative (negative values mean no limit) */

    /* if both are neg, no MIN or MAX number of args */
    /* if MIN is neg, MAX pos, no MIN number of args but a maximum */
    /* if MIN is pos, MAX neg, there is a MIN number of args, but no MAX */
    /* if both pos, then a MIN and MAX limit.  We swap to make sure it works
       right. I.E.  Min of 3 and MAX of 2 would make function unusable */
    if(MIN >= 0 && MAX >= 0)
        {
        if(MIN > MAX)
            {
            result = MIN;
            MIN = MAX;
            MAX = result;
            }
        }

    if(refmin >= 0 && refmax >= 0)
        {
        if(refmin > refmax)
            {
            result = refmin;
            refmin = MAX;
            refmax = result;
            }
        }
        
    /* It did not exist, so add it at the head */
    tmp = exprCreateFunc(name, rt, ptr, rptr, mp_ctrl, MIN, MAX, refmin, refmax);
        
    if(tmp == NULL)
        return EXPR_ERROR_MEMORY;
        
    hashmap_put(flist, name, tmp);
    return EXPR_ERROR_NOERROR;
    }

/* Get the function from a list along with it's MIN an MAX data */
EXPRERRTYPE exprFuncListGet(exprFuncList flist, char *name, exprReturnerType *rt, exprFuncType *ptr, exprRealFuncType *rptr, mpCTRL ** mp_ctrl, int *MIN, int *MAX, int *refmin, int *refmax)
    {
    exprFunc *cur;

    if(flist == NULL)
        return EXPR_ERROR_NULLPOINTER;

    if(name == NULL || name[0] == '\0')
        return EXPR_ERROR_NOTFOUND;

    /* Search for the item */
    if(hashmap_get(flist, name, (void **) &cur) == MAP_OK)	
        {
        /* We found it. */
        *rt = cur->rt;
        *ptr = cur->fptr;
        *rptr = cur->rfptr;
        *MIN = cur->MIN;
        *MAX = cur->MAX;
        *refmin = cur->refmin;
        *refmax = cur->refmax;
        // *type = cur->type;
        *mp_ctrl = cur->mp_ctrl;

        /* return now */
        return EXPR_ERROR_NOERROR;
        }
        
    /* If we got here, we did not find the item in the list */
    return EXPR_ERROR_NOTFOUND;
    }

/* This routine will free the function list */
void exprFuncListFree(exprFuncList *flist)
    {
    exprFunc *cur = NULL;
    
    if(!(*flist))
    	return;
	
	exprFuncListFreeData((*flist));
    hashmap_free((*flist));
    
 	(*flist) = NULL;
    return;
    }
    
static inline int exprFreeHashMap(void *a, void *b)
{
	free(b);
	return MAP_MISSING;
}

/* This routine will free any child nodes, and then free itself */
void exprFuncListFreeData(exprFuncList flist)
    {
		(void) hashmap_iterate(flist, exprFreeHashMap, NULL);
    }

/* This routine will create the function object */
exprFunc *exprCreateFunc(char *name, exprReturnerType rt, exprFuncType ptr, exprRealFuncType rptr, mpCTRL * mp_ctrl, int MIN, int MAX, int refmin, int refmax)
    {
    exprFunc *tmp;
    char *vtmp;

    /* We already checked the name in exprFuncListAdd */

    /* Create it */
    tmp = exprAllocMem(sizeof(exprFunc));
    if(tmp == NULL)
        return NULL;

    /* Allocate space for the name */
    vtmp = exprAllocMem(strlen(name) + 1);

    if(vtmp == NULL)
        {
        exprFreeMem(tmp);
        return NULL;
        }

    /* Copy the data over */
    strcpy(vtmp, name);
    tmp->fname = vtmp;
    tmp->rt = rt;
    tmp->fptr = ptr;
    tmp->rfptr = rptr;
    tmp->MIN = MIN;
    tmp->MAX = MAX;
    tmp->refmin = refmin;
    tmp->refmax = refmax;
    // tmp->type = type;
    tmp->mp_ctrl = mp_ctrl;
    return tmp;
    }
