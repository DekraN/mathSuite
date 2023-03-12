#ifndef __DISABLE_BLAS1D

#include "../../dutils.h"

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_DAXPY(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    daxpy(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN),  accessFloatContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN),  accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN));
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_DDOT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, ddot(mpfr_get_d(args[0], MPFR_RNDN),  accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN),  accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_DNRM2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, dnrm2(mpfr_get_d(args[0], MPFR_RNDN),  accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_DROT(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    drot(mpfr_get_d(args[0], MPFR_RNDN),  accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN),  accessFloatContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN));
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_DROTG(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{ 
	drotg( accessFloatContainerPointer(args[0]),  accessFloatContainerPointer(args[1]),  accessFloatContainerPointer(args[2]),  accessFloatContainerPointer(args[3]));
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_DSCAL(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{	
    dscal(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN),  accessFloatContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN)); 
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_DSWAP(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	dswap(mpfr_get_d(args[0], MPFR_RNDN),  accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN),  accessFloatContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN));
 	return EXPR_ERROR_NOERROR;
}

#endif
