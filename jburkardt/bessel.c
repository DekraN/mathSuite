#ifndef __DISABLE_BESSEL

#include "../../dutils.h"


__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSJ0(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessj0(mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSJ1(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessj1(mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSY0(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessy0(mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSY1(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessy1(mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSI0(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessi0(mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSI1(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessi1(mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSK0(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessk0(mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSK1(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessk1(mpfr_get_d(args[0], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSJ(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessj(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSY(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessy(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSI(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessi(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}

__MATHSUITE __JBURKARDT EXPRERRTYPE FUNCNAME_BESSK(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRTYPE args[])
{
    mpfr_set_d(val, bessk(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)), MPFR_RNDN);
 	return EXPR_ERROR_NOERROR;
}


#endif
