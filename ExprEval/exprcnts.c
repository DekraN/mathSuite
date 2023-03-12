#ifndef __DISABLE_EXPRCNTS

#include "../dutils.h"

__MATHSUITE __JBURKARDT void * FUNCNAME_EFLOAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newFloat(mpfr_get_ui(args[0], MPFR_RNDN), val, false, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EBOOL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newBool(mpfr_get_ui(args[0], MPFR_RNDN), val, false, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EUSHORT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newUShort(mpfr_get_ui(args[0], MPFR_RNDN), val, false, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_ESHORT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newShort(mpfr_get_ui(args[0], MPFR_RNDN), val, false, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newInt(mpfr_get_ui(args[0], MPFR_RNDN), val, false, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_ENEW(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	switch(mpfr_get_ui(args[0], MPFR_RNDN))
	{
		case TYPID_FLOAT:
			newFloat(mpfr_get_ui(args[1], MPFR_RNDN), val, false, NULL);
			break;
		case TYPID_BOOL:
			newBool(mpfr_get_ui(args[1], MPFR_RNDN), val, false, NULL);
			break;
		case TYPID_USHRT:
			newUShort(mpfr_get_ui(args[1], MPFR_RNDN), val, false, NULL);
			break;
		case TYPID_SHRT:
			newShort(mpfr_get_ui(args[1], MPFR_RNDN), val, false, NULL);
			break;
		case TYPID_INT:
			newInt(mpfr_get_ui(args[1], MPFR_RNDN), val, false, NULL);
			break;
		default:
			mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	}
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_FLOAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
    const register uint64_t dim = nodes->data.function.nodecount;
    ityp * tmp = calloc(dim, sizeof(ityp));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_d(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newFloat(dim, val, false, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_BOOL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
    const register uint64_t dim = nodes->data.function.nodecount;
    bool * tmp = calloc(dim, sizeof(bool));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_ui(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newBool(dim, val, false, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_USHORT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
    const register uint64_t dim = nodes->data.function.nodecount;
    dim_typ * tmp = calloc(dim, sizeof(dim_typ));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_ui(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newUShort(dim, val, false, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SHORT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
    const register uint64_t dim = nodes->data.function.nodecount;
    short * tmp = calloc(dim, sizeof(short));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_si(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newShort(dim, val, false, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_INT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
    const register uint64_t dim = nodes->data.function.nodecount;
    int * tmp = calloc(dim, sizeof(int));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_si(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newInt(dim, val, false, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_NEW(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	dim_typ pos;
	EXPRTYPE v, d1;
	*err = exprEvalNode(obj, nodes->data.function.nodes, 0, v, EXPREVAL_ALLOC);
    	
    const register uint64_t dim = nodes->data.function.nodecount-1;

    if(!err)
        {
        	mpfr_init(d1);
			switch(mpfr_get_ui(v, MPFR_RNDN))
			{
				case TYPID_FLOAT:
				{
					ityp * tmp = calloc(dim, sizeof(ityp));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_d(d1, MPFR_RNDN);
				    newFloat(dim, val, false, tmp);   
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				 
					break;
				}
				case TYPID_BOOL:
				{
					bool * tmp = calloc(dim, sizeof(bool));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_si(d1, MPFR_RNDN);
				    newBool(dim, val, false, tmp);
				    if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				
					break;
				}
				case TYPID_USHRT:
				{
					dim_typ * tmp = calloc(dim, sizeof(dim_typ));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_ui(d1, MPFR_RNDN);
					newUShort(dim, val, false, tmp);
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				
					break;
				}
				case TYPID_SHRT:
				{
					short * tmp = calloc(dim, sizeof(short));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_si(d1, MPFR_RNDN);
					newShort(dim, val, false, tmp);
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				
					break;
				}
				case TYPID_INT:
				{
					int * tmp = calloc(dim, sizeof(int));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_si(d1, MPFR_RNDN);
					newInt(dim, val, false, tmp);
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				
					break;
				}
				default:
					mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			}
				
			mpfr_clears(v, d1, NULL);
		}
		else
		{
			mpfr_clear(v); 
			return err;
		}
    
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EVFLOAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newFloat(mpfr_get_ui(args[0], MPFR_RNDN), val, true, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EVBOOL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newBool(mpfr_get_ui(args[0], MPFR_RNDN), val, true, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EVUSHORT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newUShort(mpfr_get_ui(args[0], MPFR_RNDN), val, true, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EVSHORT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newShort(mpfr_get_ui(args[0], MPFR_RNDN), val, true, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EVINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	newInt(mpfr_get_ui(args[0], MPFR_RNDN), val, true, NULL);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_EVNEW(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	switch(mpfr_get_ui(args[0], MPFR_RNDN))
	{
		case TYPID_FLOAT:
			newFloat(mpfr_get_ui(args[1], MPFR_RNDN), val, true, NULL);
			break;
		case TYPID_BOOL:
			newBool(mpfr_get_ui(args[1], MPFR_RNDN), val, true, NULL);
			break;
		case TYPID_USHRT:
			newUShort(mpfr_get_ui(args[1], MPFR_RNDN), val, true, NULL);
			break;
		case TYPID_SHRT:
			newShort(mpfr_get_ui(args[1], MPFR_RNDN), val, true, NULL);
			break;
		case TYPID_INT:
			newInt(mpfr_get_ui(args[1], MPFR_RNDN), val, true, NULL);
			break;
		default:
			mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	}
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_VFLOAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
    const register uint64_t dim = nodes->data.function.nodecount;
    ityp * tmp = calloc(dim, sizeof(ityp));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_d(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newFloat(dim, val, true, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_VBOOL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
    const register uint64_t dim = nodes->data.function.nodecount;
    bool * tmp = calloc(dim, sizeof(bool));
    	
    mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_ui(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newBool(dim, val, true, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_VUSHORT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
    const register uint64_t dim = nodes->data.function.nodecount;
    dim_typ * tmp = calloc(dim, sizeof(dim_typ));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_ui(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newUShort(dim, val, true, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_VSHORT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
    const register uint64_t dim = nodes->data.function.nodecount;
    short * tmp = calloc(dim, sizeof(short));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
		    tmp[pos] = mpfr_get_si(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newShort(dim, val, true, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_VINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	EXPRTYPE v;
	 
    const register uint64_t dim = nodes->data.function.nodecount;
    int * tmp = calloc(dim, sizeof(int));

	mpfr_init(v);
	
	for(dim_typ pos = 0; pos<dim; ++pos )
        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, v, EXPREVAL_NALLOC)))
        {
        	mpfr_clear(v);
        	return err;
        }
        else
            tmp[pos] = mpfr_get_si(v, MPFR_RNDN);
	
	mpfr_clear(v);
	newInt(dim, val, true, tmp);
	if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
		free(tmp);
		
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_VNEW(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	dim_typ pos;
	EXPRTYPE v, d1;
	*err = exprEvalNode(obj, nodes->data.function.nodes, 0, v, EXPREVAL_ALLOC);
    	
    const register uint64_t dim = nodes->data.function.nodecount-1;

    if(!err)
        {
        	mpfr_init(d1);
			switch(mpfr_get_ui(args[0], MPFR_RNDN))
			{
				case TYPID_FLOAT:
				{
					ityp * tmp = calloc(dim, sizeof(ityp));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_d(d1, MPFR_RNDN);
					newFloat(dim, val, true, NULL);
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				
					break;
				}
				case TYPID_BOOL:
				{
					bool * tmp = calloc(dim, sizeof(bool));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_si(d1, MPFR_RNDN);
					newBool(dim, val, true, tmp);
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				
					break;
				}
				case TYPID_USHRT:
				{
					dim_typ * tmp = calloc(dim, sizeof(dim_typ));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_ui(d1, MPFR_RNDN); 
					newUShort(dim, val, true, tmp);
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
						
					break;
				}
				case TYPID_SHRT:
				{
					short * tmp = calloc(dim, sizeof(short));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_si(d1, MPFR_RNDN);
					newShort(dim, val, true, tmp);
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				
					break;
				}
				case TYPID_INT:
				{
					int * tmp = calloc(dim, sizeof(int));
					for(pos = 1; pos<=dim; ++pos )
				        if((*err = exprEvalNode(obj, nodes->data.function.nodes, pos, d1, EXPREVAL_NALLOC)))
				        {
				        	mpfr_clears(v, d1, NULL);
				        	return err;
				        }
				        else
				            tmp[pos] = mpfr_get_si(d1, MPFR_RNDN);
					newInt(dim, val, true, tmp);
					if(!mpfr_cmp_ui(val, MAX_PRMCONTAINERSIZE))
						free(tmp);
				
					break;
				}
				default:
					mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			}
		
			mpfr_clears(v, d1, NULL);
		}
    
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_VOLATILE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_ui(val, _accessContainer(args[0]).pnt == NULL ? MAX_PRMCONTAINERSIZE : (_accessContainer(args[0])._volatile = !_accessContainer(args[0])._volatile), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SETARRAYSINDEX(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_ui(val, mpfLZero(args[0]) || mpfr_cmp_d(args[0], MAX_PRMCONTAINERSIZE) >= 0 ? MAX_PRMCONTAINERSIZE : (access(PRMSystem).currentIndex = mpfr_get_ui(args[0], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_FSET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{	
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((ityp*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_d(args[2], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_BSET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((bool*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_ui(args[2], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_USET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((dim_typ*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_ui(args[2], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SSET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{	
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((short*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_si(args[2], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_ISET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{	
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((int*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_si(args[2], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	switch(mpfr_get_ui(args[0], MPFR_RNDN))
	{
		case TYPID_FLOAT:
		    mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[2], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((ityp*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_d(args[3], MPFR_RNDN)), MPFR_RNDN);
			break;
		case TYPID_BOOL:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[2], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((bool*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_ui(args[3], MPFR_RNDN)), MPFR_RNDN);
			break;
		case TYPID_USHRT:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[2], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((dim_typ*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_ui(args[3], MPFR_RNDN)), MPFR_RNDN);
			break;
		case TYPID_SHRT:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[2], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((short*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_si(args[3], MPFR_RNDN)), MPFR_RNDN);
			break;
		case TYPID_INT:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[2], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : (*((int*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)) = mpfr_get_si(args[3], MPFR_RNDN)), MPFR_RNDN);
			break;
		default:
			mpfr_set_d(val, MAX_VAL, MPFR_RNDN);
	}
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_FGET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((ityp*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_BGET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((bool*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_UGET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((dim_typ*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SGET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((short*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_IGET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_d(val, _accessContainer(args[0]).pnt == NULL || mpfLZero(args[1]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((int*)(_accessContainer(args[0]).pnt) + mpfr_get_ui(args[1], MPFR_RNDN)), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_GET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	switch(mpfr_get_ui(args[0], MPFR_RNDN))
	{
		case TYPID_FLOAT:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((ityp*)(_accessContainer(args[1]).pnt) + mpfr_get_ui(args[2], MPFR_RNDN)), MPFR_RNDN);
			break;
		case TYPID_BOOL:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((bool*)(_accessContainer(args[1]).pnt) + mpfr_get_ui(args[2], MPFR_RNDN)), MPFR_RNDN);
			break;
		case TYPID_USHRT:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((dim_typ*)(_accessContainer(args[1]).pnt) + mpfr_get_ui(args[2], MPFR_RNDN)), MPFR_RNDN);
			break;
		case TYPID_SHRT:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((short*)(_accessContainer(args[1]).pnt) + mpfr_get_ui(args[2], MPFR_RNDN)), MPFR_RNDN);
			break;
		case TYPID_INT:
			mpfr_set_d(val, _accessContainer(args[1]).pnt == NULL || mpfLZero(args[2]) || mpfr_cmp_d(args[1], _accessContainer(args[0]).size) >= 0 ? MAX_VAL : *((int*)(_accessContainer(args[1]).pnt) + mpfr_get_ui(args[2], MPFR_RNDN)), MPFR_RNDN);
			break;
		default:
			mpfr_set_d(val, MAX_VAL, MPFR_RNDN);
	}
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_DEL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
	{
		mpfr_set(val, args[0], MPFR_RNDN);
		del(mpfr_get_ui(args[0], MPFR_RNDN));
	}
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_DELALL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	delAll();
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_FPRINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printItypMatrix(accessFloatContainerPointer(args[0]),
			(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[0]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[0]).size)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_BPRINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printBoolMatrix(accessBoolContainerPointer(args[0]),
			(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[0]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[0]).size)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_UPRINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printUShortMatrix( accessUShortContainerPointer(args[0]),
			(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[0]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[0]).size)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SPRINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printShortMatrix( accessShortContainerPointer(args[0]),
			(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[0]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[0]).size)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_IPRINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printIntMatrix( accessIntContainerPointer(args[0]),
			(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[0]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[0]).size)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_PRINT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	switch(mpfr_get_ui(args[0], MPFR_RNDN))
	{
		case TYPID_FLOAT:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printItypMatrix( accessFloatContainerPointer(args[1]),
					(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[1]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[1]).size)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		case TYPID_BOOL:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printBoolMatrix( accessBoolContainerPointer(args[1]),
					(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[1]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[1]).size)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		case TYPID_USHRT:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printUShortMatrix( accessUShortContainerPointer(args[1]),
					(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[1]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[1]).size)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		case TYPID_SHRT:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printShortMatrix( accessShortContainerPointer(args[1]),
					(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[1]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[1]).size)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		case TYPID_INT:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printIntMatrix( accessIntContainerPointer(args[1]),
					(dim_typ2){(uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_ROW ? 1 : _accessContainer(args[1]).size), (uint64_t) (access(PRMSystem).vectorType == VECTORTYPE_COLUMN ? 1 : _accessContainer(args[1]).size)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		default:
			mpfr_set_d(val, MAX_VAL, MPFR_RNDN);
	}

	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_FPMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printItypMatrix( accessFloatContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_BPMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printBoolMatrix( accessBoolContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_UPMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printUShortMatrix( accessUShortContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SPMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printShortMatrix( accessShortContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_IPMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	if(_accessContainer(args[0]).pnt == NULL)
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	else
		printIntMatrix( accessIntContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_PMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	switch(mpfr_get_ui(args[0], MPFR_RNDN))
	{
		case TYPID_FLOAT:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printItypMatrix( accessFloatContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		case TYPID_BOOL:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printBoolMatrix( accessBoolContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		case TYPID_USHRT:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printUShortMatrix( accessUShortContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		case TYPID_SHRT:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printShortMatrix( accessShortContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		case TYPID_INT:
			if(_accessContainer(args[1]).pnt == NULL)
				mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
			else
				printIntMatrix( accessIntContainerPointer(args[0]), (dim_typ2){_accessContainer(args[0]).size, mpfr_get_ui(args[1], MPFR_RNDN)});
			mpfr_set(val, args[0], MPFR_RNDN);
			break;
		default:
			mpfr_set_d(val, MAX_VAL, MPFR_RNDN);
	}
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_VTYPE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	mpfr_set_ui(val, (access(PRMSystem).vectorType = !access(PRMSystem).vectorType), MPFR_RNDN);
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_FMATPRM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{	
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t PRM_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
	{
		mpfr_set_ui(val, true, MPFR_RNDN);
        uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				*((ityp*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = mpfr_get_d(mtx[idx], MPFR_RNDN); 
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_BMATPRM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t PRM_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
        uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				*((bool*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = (bool) mpfr_get_ui(mtx[idx], MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_UMATPRM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t PRM_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
        uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				*((dim_typ*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = (dim_typ) mpfr_get_ui(mtx[idx], MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SMATPRM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t PRM_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
        uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				*((short*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = (short) mpfr_get_si(mtx[idx], MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_IMATPRM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t PRM_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
        uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				*((int*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = (int) mpfr_get_si(mtx[idx], MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_MATPRM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[1], MPFR_RNDN);  	
    const register uint64_t PRM_mid = mpfr_get_ui(args[2], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
        uint64_t i, j;
        register uint64_t idx;
        EXPRTYPE *mtx = mtobj->matrix;
		switch(mpfr_get_ui(args[0], MPFR_RNDN))
		{
			case TYPID_FLOAT:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						*((ityp*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = mpfr_get_d(mtx[idx], MPFR_RNDN); 
				break;
			case TYPID_BOOL:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						*((bool*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = (bool) mpfr_get_ui(mtx[idx], MPFR_RNDN);
				break;
			case TYPID_USHRT:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						*((dim_typ*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = (dim_typ) mpfr_get_ui(mtx[idx], MPFR_RNDN); 	
				break;
			case TYPID_SHRT:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						*((short*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = (short) mpfr_get_ui(mtx[idx], MPFR_RNDN);
				break;
			case TYPID_INT:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						*((int*)(accessContainer(PRM_mid).pnt) + (idx = i*dim[COLUMNS] + j)) = (int) mpfr_get_ui(mtx[idx], MPFR_RNDN); 
				break;
			default:
				mpfr_set_ui(val, false, MPFR_RNDN);
		}
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_FPRMMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t PRM_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
		uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				mpfr_set_d(mtx[(idx = i*dim[COLUMNS] + j)], *((ityp*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_BPRMMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t PRM_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
		uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				mpfr_set_ui(mtx[(idx = i*dim[COLUMNS] + j)], *((bool*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_UPRMMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{	
    const register uint64_t PRM_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
		uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				mpfr_set_ui(mtx[(idx = i*dim[COLUMNS] + j)], *((dim_typ*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_SPRMMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t PRM_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
		uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				mpfr_set_si(mtx[(idx = i*dim[COLUMNS] + j)], *((short*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_IPRMMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t PRM_mid = mpfr_get_ui(args[0], MPFR_RNDN);
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[1], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
		uint64_t i, j;
        register uint64_t idx;
		EXPRTYPE *mtx = mtobj->matrix;
		for(i=0; i<dim[ROWS]; ++i)
			for(j=0; j<dim[COLUMNS]; ++j)
				mpfr_set_si(mtx[(idx = i*dim[COLUMNS] + j)], *((int*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

__MATHSUITE __JBURKARDT void * FUNCNAME_PRMMAT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
    const register uint64_t PRM_mid = mpfr_get_ui(args[1], MPFR_RNDN);	
    const register uint64_t mathsuite_mid = mpfr_get_ui(args[2], MPFR_RNDN);
    matrixObj *mtobj = ((matrixObj *)(listNo(mathsuite_mid, MATRICES)->data));
    dim_typ2 dim =
    {
    	mtobj->dim[ROWS],
    	mtobj->dim[COLUMNS]
    };
    
    if(mathsuite_mid < 0 || mathsuite_mid >= getItemsListNo(MATRICES) || PRM_mid < 0 || PRM_mid > access(PRMSystem).currentIndex || accessContainer(PRM_mid).size != dim[ROWS]*dim[COLUMNS])
		mpfr_set_ui(val, false, MPFR_RNDN);
	else
    {
		mpfr_set_ui(val, true, MPFR_RNDN);
        uint64_t i, j;
        register uint64_t idx;
        EXPRTYPE *mtx = mtobj->matrix;
		switch(mpfr_get_ui(args[0], MPFR_RNDN))
		{
			case TYPID_FLOAT:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						mpfr_set_d(mtx[(idx = i*dim[COLUMNS] + j)], *((ityp*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
				break;
			case TYPID_BOOL:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						mpfr_set_ui(mtx[(idx = i*dim[COLUMNS] + j)],  *((bool*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN); 
				break;
			case TYPID_USHRT:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						mpfr_set_ui(mtx[(idx = i*dim[COLUMNS] + j)], *((dim_typ*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
				break;
			case TYPID_SHRT:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						mpfr_set_si(mtx[(idx = i*dim[COLUMNS] + j)], *((short*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
				break;
			case TYPID_INT:
				for(i=0; i<dim[ROWS]; ++i)
					for(j=0; j<dim[COLUMNS]; ++j)
						mpfr_set_si(mtx[(idx = i*dim[COLUMNS] + j)], *((int*)(accessContainer(PRM_mid).pnt) + idx), MPFR_RNDN);
				break;
			default:
				mpfr_set_ui(val, false, MPFR_RNDN);
		}
	}
	
	return (*err = EXPR_ERROR_NOERROR), NULL;
}

#endif
