#ifndef __DISABLE_FEM2DSERENE

#include "../../dutils.h"

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BASISSERENE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	aCheckFloat(8, itypCheckLen(8, args[6]) &&itypCheckLen(8, args[7]), val, basis_serene(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN),  accessFloatContainerPointer(args[6]),  accessFloatContainerPointer(args[7])));
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BASISDXSERENE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{	
	aCheckFloat(8, itypCheckLen(8, args[6]) &&itypCheckLen(8, args[7]), val, basis_dx_serene(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN),  accessFloatContainerPointer(args[6]),  accessFloatContainerPointer(args[7])));
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BASISDYSERENE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
	aCheckFloat(8, itypCheckLen(8, args[6]) &&itypCheckLen(8, args[7]), val, basis_dy_serene(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN),  accessFloatContainerPointer(args[6]),  accessFloatContainerPointer(args[7])));
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_FEM2DBVPSERENENODENUM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{	
    mpfr_set_d(val, fem2d_bvp_serene_node_num(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)), MPFR_RNDN);
	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_NOT1D(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, not1d(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)), MPFR_RNDN);
	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_NOT1(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, not1(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_NOT2DX(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, not2dy(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN)), MPFR_RNDN);
	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_NOT2DY(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, not2dy(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN)), MPFR_RNDN);
	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_NOT2(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{	
    mpfr_set_d(val, not2(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN)), MPFR_RNDN);
	return EXPR_ERROR_NOERROR;
}

#endif
