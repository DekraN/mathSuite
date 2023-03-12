#include "dutils.h" 

__MATHSUITE __JBURKARDT inline uint64_t RC_MUL( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(args[opr[0]], MPFR_RNDN)*mpfr_get_ui(args[opr[1]], MPFR_RNDN);
}

__MATHSUITE __JBURKARDT inline uint64_t RC_ADD( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(args[opr[0]], MPFR_RNDN)+mpfr_get_ui(args[opr[1]], MPFR_RNDN);
}

__MATHSUITE __JBURKARDT inline uint64_t RC_ADDADDCONST( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(args[opr[0]], MPFR_RNDN)+mpfr_get_ui(args[opr[1]], MPFR_RNDN)+opr[2];
}

__MATHSUITE __JBURKARDT inline uint64_t RC_3MUL( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(args[opr[0]], MPFR_RNDN)*mpfr_get_ui(args[opr[1]], MPFR_RNDN)*mpfr_get_ui(args[opr[2]], MPFR_RNDN);
}

__MATHSUITE __JBURKARDT inline uint64_t RC_4MUL( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(args[opr[0]], MPFR_RNDN)*mpfr_get_ui(args[opr[1]], MPFR_RNDN)*mpfr_get_ui(args[opr[2]], MPFR_RNDN)*mpfr_get_ui(args[opr[3]], MPFR_RNDN);
}

__MATHSUITE __JBURKARDT inline uint64_t RC_TRIANGLENUMADD( EXPRTYPE args[], sel_typ opr[])
{
	return triangle_num(mpfr_get_ui(args[opr[0]], MPFR_RNDN)+opr[1]);
}

__MATHSUITE __JBURKARDT inline uint64_t RC_MULTRIANGLENUMADD( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(args[opr[0]], MPFR_RNDN)*triangle_num(mpfr_get_ui(args[opr[1]], MPFR_RNDN)+opr[2]);
}

__MATHSUITE __JBURKARDT inline uint64_t RC_SCALADDCONSTMUL( EXPRTYPE args[], sel_typ opr[])
{
	return ((opr[0]*mpfr_get_ui(args[opr[1]], MPFR_RNDN))+opr[2])*mpfr_get_ui(args[opr[3]], MPFR_RNDN);
}

__MATHSUITE __JBURKARDT inline uint64_t RC_ADDCONSTSCAL( EXPRTYPE args[], sel_typ opr[])
{
	return ((opr[0]+mpfr_get_ui(args[opr[1]], MPFR_RNDN))*opr[2]);
}

__MATHSUITE __JBURKARDT inline uint64_t RC_SCAL( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(args[opr[0]], MPFR_RNDN)*opr[1];
}

__MATHSUITE __JBURKARDT inline uint64_t RC_PPMULSCAL( EXPRTYPE args[], sel_typ opr[])
{
	return ((mpfr_get_ui(args[opr[0]], MPFR_RNDN)+opr[1])*(mpfr_get_ui(args[opr[2]], MPFR_RNDN)+opr[3]))*opr[4];
}

__MATHSUITE __JBURKARDT inline uint64_t RC_PPMULDIV( EXPRTYPE args[], sel_typ opr[])
{
	return ((mpfr_get_ui(args[opr[0]], MPFR_RNDN)+opr[1])*(mpfr_get_ui(args[opr[2]], MPFR_RNDN)+opr[3]))/opr[4];
}

__MATHSUITE __JBURKARDT inline uint64_t RC_PPPMUL2MULDIV( EXPRTYPE args[], sel_typ opr[])
{
	return ((mpfr_get_ui(args[opr[0]], MPFR_RNDN)+opr[1])*(mpfr_get_ui(args[opr[2]], MPFR_RNDN)+opr[3]))*mpfr_get_ui(args[opr[4]], MPFR_RNDN)*mpfr_get_ui(args[opr[5]], MPFR_RNDN)/opr[6];
}

__MATHSUITE __JBURKARDT inline uint64_t RC_2ADD2CONSTDIV( EXPRTYPE args[], sel_typ opr[])
{
	return ((mpfr_get_ui(args[opr[0]], MPFR_RNDN)+mpfr_get_ui(args[opr[1]], MPFR_RNDN)+opr[2])*(mpfr_get_ui(args[opr[3]], MPFR_RNDN)+mpfr_get_ui(args[opr[4]], MPFR_RNDN)+opr[5]))/opr[6];
}

__MATHSUITE __JBURKARDT inline uint64_t RC_PPMULADD( EXPRTYPE args[], sel_typ opr[])
{
	return ((mpfr_get_ui(args[opr[0]], MPFR_RNDN)+opr[2])*(mpfr_get_ui(args[opr[1]], MPFR_RNDN)+opr[3]))+opr[4];
}

__MATHSUITE __JBURKARDT inline uint64_t RC_ADDCONST( EXPRTYPE args[], sel_typ opr[])
{
	return (mpfr_get_ui(args[opr[0]], MPFR_RNDN))+opr[1];
}

/*
__MATHSUITE __JBURKARDT inline uint64_t RC_DEREF( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(*(nodes->data.function.refs[opr[0]]), MPFR_RNDN);
}
*/

/*
__MATHSUITE __JBURKARDT inline uint64_t RC_DEREFSCAL( EXPRTYPE args[], sel_typ opr[])
{
	return mpfr_get_ui(*(nodes->data.function.refs[opr[0]]), MPFR_RNDN)*opr[1];
}
*/

__MATHSUITE __JBURKARDT inline uint64_t RC_WAVELET( EXPRTYPE args[], sel_typ opr[])
{
	return pow(2,mpfr_get_ui(args[opr[0]], MPFR_RNDN)) + mpfr_get_ui(args[opr[1]], MPFR_RNDN) + (pow(2,mpfr_get_ui(args[opr[2]], MPFR_RNDN))-1)*mpfr_get_ui(args[opr[3]], MPFR_RNDN) - ((pow(2,mpfr_get_ui(args[opr[4]], MPFR_RNDN))-1)*2);
}

__MATHSUITE __JBURKARDT inline void RE_PINT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckInt((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PUINT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckInt((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PSHRT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckShort((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PUSHRT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckUShort((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PCHR( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckUShort((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PUCHR( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckBool((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PLNG( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckInt((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PULNG( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckInt((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PLLNG( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckInt((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PULLNG( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckInt((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PFLT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckFloat((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PDBL( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckFloat((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_PLDBL( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	aCheckFloat((nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.ctrl_function)(args, nodes->data.function.mp_ctrl->ReturnCTRLSystem.ReturnCTRLSystemPatterned.opr), val, data);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_INT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_si(val, R_INT(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_UINT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_ui(val, R_UINT(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_SHRT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_si(val, R_SHRT(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_USHRT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_ui(val, R_USHRT(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_CHR( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_si(val, R_CHR(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_UCHR( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_ui(val, R_UCHR(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_LNG( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_si(val, R_LNG(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_ULNG( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_ui(val, R_ULNG(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_LLNG( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_si(val, R_LLNG(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_ULLNG( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_ui(val, R_ULLNG(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_FLT( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_flt(val, R_FLT(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_DBL( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_d(val, R_DBL(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_LDBL( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	mpfr_set_ld(val, R_LDBL(data), MPFR_RNDN);
	return;
}

__MATHSUITE __JBURKARDT inline void RE_CPLX( void * data, EXPRTYPE args[], exprNode * nodes, mpfr_t val)
{
	//	unimplemented
	return;
}

__MATHSUITE __JBURKARDT inline void * H_CPLX (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return NULL; // unimplemented
}

__MATHSUITE __JBURKARDT inline void * H_2DT3PCX (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return NULL; // unimplemented
}

__MATHSUITE __JBURKARDT inline void * H_3DT3PCX (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return NULL; // unimplemented
}

__MATHSUITE __JBURKARDT inline void * H_PINT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessIntContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessIntContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessShortContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessUShortContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessUShortContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessBoolContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessIntContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessIntContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessIntContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessIntContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessFloatContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessFloatContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(accessFloatContainerPointer(args[0]));
}

__MATHSUITE __JBURKARDT inline void * H_PINT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT2(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PINT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT3(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PINT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT4(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PINT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT5(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PINT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT6(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PINT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT7(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PINT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT8(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PINT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT9(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PINT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PINT10(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN), mpfr_get_si(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPINT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPINT10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT2(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT3(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT4(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT5(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT6(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT7(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT8(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT9(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUINT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUINT10(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN), mpfr_get_si(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUINT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUINT10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT2(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT3(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT4(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT5(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT6(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT7(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT8(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT9(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PSHRT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PSHRT10(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN), mpfr_get_si(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT2(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT3(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1]), accessShortContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT4(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessShortContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT5(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessShortContainerPointer(args[3]), accessShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT6(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessShortContainerPointer(args[3]), accessShortContainerPointer(args[4]), accessShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT7(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessShortContainerPointer(args[3]), accessShortContainerPointer(args[4]), accessShortContainerPointer(args[5]), accessShortContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT8(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessShortContainerPointer(args[3]), accessShortContainerPointer(args[4]), accessShortContainerPointer(args[5]), accessShortContainerPointer(args[6]), accessShortContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT9(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessShortContainerPointer(args[3]), accessShortContainerPointer(args[4]), accessShortContainerPointer(args[5]), accessShortContainerPointer(args[6]), accessShortContainerPointer(args[7]), accessShortContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPSHRT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPSHRT10(accessShortContainerPointer(args[0]), accessShortContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessShortContainerPointer(args[3]), accessShortContainerPointer(args[4]), accessShortContainerPointer(args[5]), accessShortContainerPointer(args[6]), accessShortContainerPointer(args[7]), accessShortContainerPointer(args[8]), accessShortContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT2(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT3(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT4(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT5(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT6(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT7(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT8(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT9(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUSHRT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUSHRT10(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN), mpfr_get_ui(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUSHRT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUSHRT10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR2(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR3(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR4(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR5(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR6(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR7(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR8(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR9(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PCHR10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCHR10(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN), mpfr_get_si(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR2(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR3(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1]), accessBoolContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR4(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1]), accessBoolContainerPointer(args[2]), accessBoolContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR5(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1]), accessBoolContainerPointer(args[2]), accessBoolContainerPointer(args[3]), accessBoolContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR6(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1]), accessBoolContainerPointer(args[2]), accessBoolContainerPointer(args[3]), accessBoolContainerPointer(args[4]), accessBoolContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR7(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1]), accessBoolContainerPointer(args[2]), accessBoolContainerPointer(args[3]), accessBoolContainerPointer(args[4]), accessBoolContainerPointer(args[5]), accessBoolContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR8(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1]), accessBoolContainerPointer(args[2]), accessBoolContainerPointer(args[3]), accessBoolContainerPointer(args[4]), accessBoolContainerPointer(args[5]), accessBoolContainerPointer(args[6]), accessBoolContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR9(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1]), accessBoolContainerPointer(args[2]), accessBoolContainerPointer(args[3]), accessBoolContainerPointer(args[4]), accessBoolContainerPointer(args[5]), accessBoolContainerPointer(args[6]), accessBoolContainerPointer(args[7]), accessBoolContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPCHR10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPCHR10(accessBoolContainerPointer(args[0]), accessBoolContainerPointer(args[1]), accessBoolContainerPointer(args[2]), accessBoolContainerPointer(args[3]), accessBoolContainerPointer(args[4]), accessBoolContainerPointer(args[5]), accessBoolContainerPointer(args[6]), accessBoolContainerPointer(args[7]), accessBoolContainerPointer(args[8]), accessBoolContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR2(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR3(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR4(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR5(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR6(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR7(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR8(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR9(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PUCHR10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PUCHR10(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN), mpfr_get_ui(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPUCHR10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPUCHR10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG2(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG3(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG4(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG5(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG6(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG7(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG8(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG9(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLNG10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLNG10(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN), mpfr_get_si(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLNG10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLNG10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG2(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG3(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG4(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG5(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG6(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG7(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG8(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG9(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULNG10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULNG10(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN), mpfr_get_ui(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULNG10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULNG10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG2(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG3(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG4(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG5(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG6(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG7(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG8(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG9(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLLNG10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLLNG10(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_si(args[5], MPFR_RNDN), mpfr_get_si(args[6], MPFR_RNDN), mpfr_get_si(args[7], MPFR_RNDN), mpfr_get_si(args[8], MPFR_RNDN), mpfr_get_si(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLLNG10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLLNG10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG2(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG3(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG4(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG5(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG6(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG7(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG8(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG9(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PULLNG10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PULLNG10(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN), mpfr_get_ui(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPULLNG10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPULLNG10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT2(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT3(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), mpfr_get_flt(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT4(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), mpfr_get_flt(args[2], MPFR_RNDN), mpfr_get_flt(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT5(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), mpfr_get_flt(args[2], MPFR_RNDN), mpfr_get_flt(args[3], MPFR_RNDN), mpfr_get_flt(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT6(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), mpfr_get_flt(args[2], MPFR_RNDN), mpfr_get_flt(args[3], MPFR_RNDN), mpfr_get_flt(args[4], MPFR_RNDN), mpfr_get_flt(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT7(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), mpfr_get_flt(args[2], MPFR_RNDN), mpfr_get_flt(args[3], MPFR_RNDN), mpfr_get_flt(args[4], MPFR_RNDN), mpfr_get_flt(args[5], MPFR_RNDN), mpfr_get_flt(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT8(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), mpfr_get_flt(args[2], MPFR_RNDN), mpfr_get_flt(args[3], MPFR_RNDN), mpfr_get_flt(args[4], MPFR_RNDN), mpfr_get_flt(args[5], MPFR_RNDN), mpfr_get_flt(args[6], MPFR_RNDN), mpfr_get_flt(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT9(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), mpfr_get_flt(args[2], MPFR_RNDN), mpfr_get_flt(args[3], MPFR_RNDN), mpfr_get_flt(args[4], MPFR_RNDN), mpfr_get_flt(args[5], MPFR_RNDN), mpfr_get_flt(args[6], MPFR_RNDN), mpfr_get_flt(args[7], MPFR_RNDN), mpfr_get_flt(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PFLT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PFLT10(mpfr_get_flt(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), mpfr_get_flt(args[2], MPFR_RNDN), mpfr_get_flt(args[3], MPFR_RNDN), mpfr_get_flt(args[4], MPFR_RNDN), mpfr_get_flt(args[5], MPFR_RNDN), mpfr_get_flt(args[6], MPFR_RNDN), mpfr_get_flt(args[7], MPFR_RNDN), mpfr_get_flt(args[8], MPFR_RNDN), mpfr_get_flt(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT2(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT3(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT4(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT5(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT6(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT7(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT8(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT9(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPFLT10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPFLT10(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL2(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL3(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL4(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL5(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL6(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL7(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL8(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL9(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL10(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN), mpfr_get_d(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL11 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL11(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN), mpfr_get_d(args[9], MPFR_RNDN), mpfr_get_d(args[10], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PDBL12 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDBL12(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN), mpfr_get_d(args[9], MPFR_RNDN), mpfr_get_d(args[10], MPFR_RNDN), mpfr_get_d(args[11], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL2(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL3(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL4(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL5(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL6(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL7(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL8(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL9(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPDBL10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPDBL10(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL2(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL3(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL4(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL5(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL6(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL7(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL8(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL9(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PLDBL10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PLDBL10(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN), mpfr_get_d(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL2 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL2(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL3 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL3(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL4 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL4(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL5 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL5(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL6 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL6(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL7 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL7(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL8 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL8(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL9 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL9(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PPLDBL10 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPLDBL10(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}


__MATHSUITE __JBURKARDT inline void * H_SINT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SINT(mpfr_get_si(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SUINT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SUINT(mpfr_get_ui(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SSHRT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SSHRT(mpfr_get_si(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SUSHRT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SUSHRT(mpfr_get_ui(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SCHR (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SCHR(mpfr_get_si(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SUCHR (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SUCHR(mpfr_get_ui(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SLNG(mpfr_get_si(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SULNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SULNG(mpfr_get_ui(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SLLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SLLNG(mpfr_get_si(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SULLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SULLNG(mpfr_get_ui(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SFLT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SFLT(mpfr_get_flt(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SDBL (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SDBL(mpfr_get_d(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_SLDBL (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SLDBL(mpfr_get_si(args[0], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_ITB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITB(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT2IPITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT2IPITPDT(accessFloatContainerPointer(args[0]), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PITDTPDTPITPDTDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PITDTPDTPITPDTDTPIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessUShortContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT4PIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2PDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PDT2PIT(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT5PITPDT2DT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT5PITPDT2DT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessUShortContainerPointer(args[6]), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_4PITPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4PITPDTPIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT3PIT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT2PIT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_4PDTPITPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4PDTPITPB(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessBoolContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_B3DT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_B3DT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_5PDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_5PDTPIT(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_2BDTPITPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2BDTPITPDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PIT2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PIT2PDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PIT2PDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PIT2PDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2I2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2I2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DT3ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3ITPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_3PDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3PDT2PIT(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTITPITDTPITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTITPITDTPITDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3DT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITIPITI2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITIPITI2IT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), mpfr_get_si(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTITPITI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTITPITI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITIPITI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITIPITI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_7ITFITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_7ITFITPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), getFID(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_4PITFITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4PITFITPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), getFID(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_4PITFIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4PITFIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), getFID(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2ITPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2ITPITPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITITDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITITDTPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTIT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTIT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDT3PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2IT2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT6PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT6PIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessUShortContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT6PITPSPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT6PITPSPIT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessShortContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT4PITPSPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT4PITPSPIT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessShortContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT5PITPSPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT5PITPSPIT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessShortContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PS2PIT2PUL (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PS2PIT2PUL(accessShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_IPS2PIT2PUL7PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IPS2PIT2PUL7PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessShortContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10]), accessFloatContainerPointer(args[11]), accessFloatContainerPointer(args[12])));
}

__MATHSUITE __JBURKARDT inline void * H_IPS4PIT2PUL4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IPS2PIT2PUL7PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessShortContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10]), accessFloatContainerPointer(args[11]), accessFloatContainerPointer(args[12])));
}

__MATHSUITE __JBURKARDT inline void * H_PS4PIT2PUL (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PS4PIT2PUL(accessShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDTPIT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT2PITPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT2PITPDTPIT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_3PDT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3PDT3PIT(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2IPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2IPI(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3IT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTFITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTFITPIT(mpfr_get_ui(args[0], MPFR_RNDN), getFID(args[1]), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTFIT2ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTFIT2ITPIT(mpfr_get_ui(args[0], MPFR_RNDN), getFID(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTIT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTIT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITIT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITIT2PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPITDT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPITDT4PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_IT2PITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT2PITDT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_IT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT3PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_IT3PITDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT3PITDTPI(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPIT2ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPIT2ITPIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_IT3PITDTPI2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT3PITDTPI2PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDTPI2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDTPI2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2ITPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2ITPIPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPIPB(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessBoolContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPII (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPII(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DT6PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT6PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT2PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_2IT2PITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2IT2PITPDT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2ITPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_FIT2ITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_FIT2ITDT(getFID(args[0]), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_IPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IPI(mpfr_get_si(args[0], MPFR_RNDN), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_4DT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4DT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT3PIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT3PIPB(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessBoolContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_FDTDT2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_FDTDT2PDT(accessUShortContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), accessUShortContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_IT2DT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT2DT2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITDTPITPS (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITDTPITPS(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITDTPS (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITDTPS(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITPDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITPDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PIT2ITPDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PIT2ITPDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), accessUShortContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITI2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITI2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PITDTPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PITDTPIPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITITPI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2PITITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PITITPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4FITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4FITPIT(mpfr_get_ui(args[0], MPFR_RNDN), getFID(args[1]), getFID(args[2]), getFID(args[3]), getFID(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3FITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3FITPIT(mpfr_get_ui(args[0], MPFR_RNDN), getFID(args[1]), getFID(args[2]), getFID(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITFIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITFIT(mpfr_get_ui(args[0], MPFR_RNDN), getFID(args[1]), getFID(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4IT2FITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4IT2FITPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), getFID(args[5]), getFID(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITDTFIT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITDTFIT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), getFID(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPI3PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPI3PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDT2PITDT2PITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDT2PITDT2PITDT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), mpfr_get_ui(args[8], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_6IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_6IT2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3DT2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DT2IT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_4DT3IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4DT3IT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2PDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PDTPIT(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessFloatContainerPointer(args[2])));
}
__MATHSUITE __JBURKARDT inline void * H_PIT3DTIT2PIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT3DTIT2PIPIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT3DTITPI3PITPIPITI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT3DTITPI3PITPIPITI(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessIntContainerPointer(args[9]), accessFloatContainerPointer(args[10]), mpfr_get_si(args[11], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT3DTIT2PI3PITPIPITI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT3DTIT2PI3PITPIPITI(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessIntContainerPointer(args[10]), accessFloatContainerPointer(args[11]), mpfr_get_si(args[12], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT4DT3PITPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT4DT3PITPIPIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessUShortContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessFloatContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_IT5PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT5PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_3IT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3IT4PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPIT4ITPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPIT4ITPDTPIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessUShortContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT3ITDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT3ITDTPIT(accessFloatContainerPointer(args[0]), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_IT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT4PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPIT2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPIT2IT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3IT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3IT3PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2PIT2ITPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PIT2ITPDTPIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessUShortContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_CHITDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_CHITDT2PIT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPIT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPIT2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2PITITDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PITITDTPI(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessUShortContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT2DTPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT2DTPDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_ITPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPITPI(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPITPDTPI3DTPIT2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPITPDTPI3DTPIT2PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessIntContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN), accessFloatContainerPointer(args[9]), accessUShortContainerPointer(args[10]), accessUShortContainerPointer(args[11])));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPIT6PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPIT6PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessUShortContainerPointer(args[5]), accessUShortContainerPointer(args[6]), accessUShortContainerPointer(args[7]), accessUShortContainerPointer(args[8]), accessUShortContainerPointer(args[9]), accessUShortContainerPointer(args[10])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT5ITDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT5ITDTPIT(accessFloatContainerPointer(args[0]), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_6IT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_6IT3PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_4IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4IT2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_ITDT3ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITDT3ITPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_4DTPIT3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4DTPIT3PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_3PITPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3PITPB(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessBoolContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_4ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4ITPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_6ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_6ITPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_4IT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4IT3PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_6ITPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_6ITPDTPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessUShortContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_8IT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_8IT3PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10])));
}

__MATHSUITE __JBURKARDT inline void * H_6IT5PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_6IT5PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10])));
}

__MATHSUITE __JBURKARDT inline void * H_3PITDT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3PITDT3PIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DT5PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT5PIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_4PITDT2PITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4PITDT2PITPDT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessUShortContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_10ITPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_10ITPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN), mpfr_get_d(args[9], MPFR_RNDN), accessFloatContainerPointer(args[10])));
}

__MATHSUITE __JBURKARDT inline void * H_2PIT4IT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PIT4IT3PIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_4ITPITPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4ITPITPIPIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PI2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PI2PIT(accessIntContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPDTPI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT2DTPDTDTPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT2DTPDTDTPDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessUShortContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPIT2PDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPIT2PDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPITIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPITIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT4PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITIT2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITIT2PDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_IT3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT3PI(mpfr_get_d(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIPDTPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIPDTPB(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessBoolContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PIPB(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessBoolContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PDT(mpfr_get_ui(args[0], MPFR_RNDN), accessUShortContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessUShortContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_3DT2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DT2PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPIT4DTPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPIT4DTPDTPIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), accessUShortContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPIT3DTPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPIT3DTPDTPIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessUShortContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPITDTPITDTPIT2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPITDTPITDTPIT2IT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_IT6PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT6PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PIT3PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PIT3PDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_4ITDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4ITDT2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PIPDTPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIPDTPB(accessIntContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessBoolContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_2PITITCH (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PITITCH(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2PITPCHIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PITPCHIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessBoolContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3PITPBPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3PITPBPIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessBoolContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT3ITPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT3ITPDTPIT(accessFloatContainerPointer(args[0]), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessUShortContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2PITPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PITPB(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessBoolContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4PDT(mpfr_get_ui(args[0], MPFR_RNDN), accessUShortContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_4DT4PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4DT4PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5]), accessUShortContainerPointer(args[6]), accessUShortContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_PDTPI3DTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDTPI3DTPI(accessUShortContainerPointer(args[0]), accessIntContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT7PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT7PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4PITDT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4PITDT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPDT(mpfr_get_ui(args[0], MPFR_RNDN), accessUShortContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPIPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPIPDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2PDT(mpfr_get_ui(args[0], MPFR_RNDN), accessUShortContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessUShortContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPI2PDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPI2PDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessUShortContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessUShortContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT2DTPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT2DTPDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT2DTPITPDTDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT2DTPITPDTDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessUShortContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2IT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT4IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT4IT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPIIPITPI2PITDT4IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPIIPITPI2PITDT4IT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), mpfr_get_ui(args[8], MPFR_RNDN), mpfr_get_d(args[9], MPFR_RNDN), mpfr_get_d(args[10], MPFR_RNDN), mpfr_get_d(args[11], MPFR_RNDN), mpfr_get_d(args[12], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPDT2PITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPDT2PITPDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTITPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTITPITPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_ITI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITI(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_IT4PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IT4PI(mpfr_get_d(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2I2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2I2PI(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_I2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_I2PI(mpfr_get_si(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PPPI3DT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPPI3DT(accessIntContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PPI2DT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPI2DT(accessIntContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPI2DT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPI2DT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPII2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPII2PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPI2I2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPI2I2PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), accessUShortContainerPointer(args[5]), accessUShortContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIPS (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPIPS(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessShortContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTI2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTI2PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), accessUShortContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_3I3PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3I3PDT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_3I2F (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3I2F(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_flt(args[3], MPFR_RNDN), mpfr_get_flt(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3I2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3I2IT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPI2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPI2PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_3DT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DT2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2I2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2I2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2IPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2IPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIPDT(mpfr_get_ui(args[0], MPFR_RNDN), accessUShortContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTIPIIPII (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTIPIIPII(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), mpfr_get_si(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPII2PS (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPII2PS(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN), accessShortContainerPointer(args[3]), accessShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT3PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPI2DTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPI2DTPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDTPI(accessUShortContainerPointer(args[0]), accessIntContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_PDTPII (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDTPII(accessUShortContainerPointer(args[0]), accessIntContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTI2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTI2PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PII3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PII3PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIPDT2PIPS (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIPDT2PIPS(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_3PII (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3PII(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DT4PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2PIPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PIPDT(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessUShortContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIPDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIPDTPI(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessUShortContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPI2I (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPI2I(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIIPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT3PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), mpfr_get_si(args[4], MPFR_RNDN), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIDTPDT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIDTPDT2PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessUShortContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PII (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PII(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_IDTPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IDTPB(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessBoolContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_3DT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIT2DT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPIT2DT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITIT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITIT2PI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_PCDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PCDT(accessBoolContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DTITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTITPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPIPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITI2PIT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITI2PIT2PI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PDT2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDT2IT(accessUShortContainerPointer(args[0]), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DT4PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT4PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPITPDT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessUShortContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2ITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2ITPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT3PITDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT3PITDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_9PLI4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_9PLI4PIT(accessIntContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10]), accessFloatContainerPointer(args[11]), accessFloatContainerPointer(args[12])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITDT2PIDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITDT2PIDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT3PI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PIDTPITDT3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PIDTPITDT3PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDTPITDT6PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDTPITDT6PI(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8]), accessIntContainerPointer(args[9]), accessIntContainerPointer(args[10])));
}

__MATHSUITE __JBURKARDT inline void * H_2PITPI2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PITPI2PIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDT3PITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDT3PITPDT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessUShortContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_5DT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_5DT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_7DT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_7DT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_9DT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_9DT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10]), accessFloatContainerPointer(args[11])));
}

__MATHSUITE __JBURKARDT inline void * H_3IPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3IPI(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIPB(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessBoolContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIPB2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIPB2PDT(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessBoolContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_3I3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3I3PI(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPDT2PIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPDT2PIPB(mpfr_get_ui(args[0], MPFR_RNDN), accessUShortContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessBoolContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4PIPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4PIPDT(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4I2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4I2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPITPI2PDTPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPITPI2PDTPITPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPDT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPDT2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4PIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4PIPB(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessBoolContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PIPITPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PIPITPB(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessBoolContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2PIPITITPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2PIPITITPB(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4]), mpfr_get_d(args[5], MPFR_RNDN), accessBoolContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPI2DTPIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPI2DTPIPB(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessBoolContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT3PI2ITPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT3PI2ITPB(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), accessBoolContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2PI2ITPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2PI2ITPB(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessBoolContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PIPB(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessBoolContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITDTPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITDTPIPIT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2ITI3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2ITI3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2ITI4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2ITI4PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPDT2PITDTPITDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPDT2PITDTPITDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessUShortContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2ITDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2ITDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2ITDTPITPDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2ITDTPITPDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessUShortContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2ITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2ITDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT2DT2PIT2DT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT2DT2PIT2DT4PIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT4DTPIT2DT2PITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT4DTPIT2DT2PITDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), mpfr_get_ui(args[10], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT4DTPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT4DTPDTPIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessUShortContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT4DTPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT4DTPIPIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT4DTPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT4DTPDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PIT4DTPDTPITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT4DTPDTPITDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessUShortContainerPointer(args[5]), accessFloatContainerPointer(args[6]), mpfr_get_ui(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT2DTPI2PITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT2DTPI2PITDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT3DTPITPIPITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT3DTPITPIPITDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6]), mpfr_get_ui(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT3DT7PITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT3DT7PITDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10]), mpfr_get_ui(args[11], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITPDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2PI3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PI3PIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessFloatContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_PIT3DT3PITDTPITDTPITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIT3DT3PITDTPITDTPITDT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), mpfr_get_ui(args[7], MPFR_RNDN), accessFloatContainerPointer(args[8]), mpfr_get_ui(args[9], MPFR_RNDN), accessFloatContainerPointer(args[10]), mpfr_get_ui(args[11], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPCXDTPCX (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return NULL; // unimplemented
	// return nodes->data.function.rfptr(C_2DT2PCX(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDTPITDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDTPITDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTITPCXIT2PCX (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return NULL; // unimplemented
	// return nodes->data.function.rfptr(C_2DTITPCXIT2PCX(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTITPCX (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return NULL; // unimplemented
	// return nodes->data.function.rfptr(C_2DTITPCX(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPCX (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return NULL; // unimplemented
	// return nodes->data.function.rfptr(C_2DTPCX(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPCX (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return NULL; // unimplemented
	// return nodes->data.function.rfptr(C_DTPCX(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2PITIT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PITIT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_flt(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPITPI3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPITPI3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT4PIT2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT4PIT2IT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DT3PPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT3PPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPIPDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessUShortContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4IT2FITPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4IT2FITPDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), getFID(args[5]), getFID(args[6]), accessUShortContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_4ITDT5PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4ITDT5PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDT2PITIT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDT2PITIT3PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), mpfr_get_d(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_ITDT8PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITDT8PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDT8PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDT8PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPITIT2PIT2DT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPITIT2PIT2DT2PIT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITDTITDT2IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITDTITDT2IT2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2ITPPIT2DT2PPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2ITPPIT2DT2PPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_BDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_BDTPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTS2IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTS2IT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_3DT2PI2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DT2PI2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2DTPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2DTPIPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITDTITPITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITDTITPITPDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessUShortContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_3DT6PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DT6PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_DT4IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT4IT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPB(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessBoolContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_PITPIPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PITPIPITPI(accessFloatContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPIPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPIPITPI(mpfr_get_d(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2ITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2ITPI(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_ITIPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITIPI(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_IDTIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IDTIT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3IT2I (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3IT2I(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN), mpfr_get_si(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3ITI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3ITI(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_I2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_I2IT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIIPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_si(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIT2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPIT2PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPIIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPIIPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTIPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPITDTIT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPITDTIT2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DTIPDTDT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTIPDTDT2PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PPIT2DT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PPIT2DT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIDTPITPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIDTPITPIPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPIITPDTPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPIITPDTPITPI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN), accessUShortContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPIPDTPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPIPDTPITPI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_PDTPITPIIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDTPITPIIT(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPIIT3PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPIIT3PDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5]), accessUShortContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITPDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITPDTPI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPI2IT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPI2IT2PI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PDTPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDTPITPI(accessUShortContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2IT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2IT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2PI2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2PI2IT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITDTIT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITDTIT2PI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITITPDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITITPDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), accessUShortContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITITDTPDTPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITITDTPDTPITPI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessUShortContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPB(mpfr_get_ui(args[0], MPFR_RNDN), accessBoolContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PITITDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PITITDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2PITITDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PITITDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_I2DT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_I2DT2PIT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_IDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_IDT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPIPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2ITFIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2ITFIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), getFID(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITDTPIIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITDTPIIPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), mpfr_get_si(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2DTPI3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2DTPI3PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PITDT2PIDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PITDT2PIDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_7DTPDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_7DTPDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), accessUShortContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_ITPIT3DT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_ITPIT3DT(mpfr_get_d(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DT3IT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3IT4PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PITPI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPI2PITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPI2PITPI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_3ITDTPIT2DT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3ITDTPIT2DT2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_4IT2FITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_4IT2FITPI(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), getFID(args[4]), getFID(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_PITDT2IT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PITDT2IT4PIT(accessFloatContainerPointer(args[0]), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIPSPIB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIPSPIB(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessIntContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2PIDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PIDTPI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPITPI2DTPIPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPITPI2DTPIPITPI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessIntContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_ui(args[5], MPFR_RNDN), accessIntContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIDTPIT2PIDTPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIDTPIT2PIDTPIT(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), accessFloatContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_DTI2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTI2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_3ITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3ITPI(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT4IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT4IT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPIDTPI2DT2PITDTIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPIDTPI2DT2PITDTIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), mpfr_get_ui(args[9], MPFR_RNDN), mpfr_get_d(args[10], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2PITDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PITDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT2ITDT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2ITDT2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_d(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2PITPI2DTPI4DT2ITDT2IT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PITPI2DTPI4DT2ITDT2IT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN), mpfr_get_ui(args[8], MPFR_RNDN), mpfr_get_ui(args[9], MPFR_RNDN), mpfr_get_d(args[10], MPFR_RNDN), mpfr_get_d(args[11], MPFR_RNDN), mpfr_get_ui(args[12], MPFR_RNDN), mpfr_get_d(args[13], MPFR_RNDN), mpfr_get_d(args[14], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_2PITPI2DTPIDTPIT3DT2ITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PITPI2DTPIDTPIT3DT2ITDT(accessFloatContainerPointer(args[0]), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN), accessFloatContainerPointer(args[7]), mpfr_get_ui(args[8], MPFR_RNDN), mpfr_get_ui(args[9], MPFR_RNDN), mpfr_get_ui(args[10], MPFR_RNDN), mpfr_get_d(args[11], MPFR_RNDN), mpfr_get_d(args[12], MPFR_RNDN), mpfr_get_ui(args[13], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPIDTPI2DT4ITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPIDTPI2DT4ITDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN), mpfr_get_d(args[9], MPFR_RNDN), mpfr_get_d(args[10], MPFR_RNDN), mpfr_get_ui(args[11], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITPIDTPI2DT2ITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITPIDTPI2DT2ITDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessIntContainerPointer(args[2]), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), mpfr_get_ui(args[6], MPFR_RNDN), mpfr_get_d(args[7], MPFR_RNDN), mpfr_get_d(args[8], MPFR_RNDN), mpfr_get_ui(args[9], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPITDT4ITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPITDT4ITDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), mpfr_get_d(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), mpfr_get_d(args[6], MPFR_RNDN), mpfr_get_ui(args[7], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3DTITPITDT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTITPITDT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_d(args[3], MPFR_RNDN), accessFloatContainerPointer(args[4]), mpfr_get_ui(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_2C3DTITPITDTPITDTITPITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2C3DTITPITDTPITDTITPITDT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), mpfr_get_ui(args[4], MPFR_RNDN), mpfr_get_d(args[5], MPFR_RNDN), accessFloatContainerPointer(args[6]), mpfr_get_ui(args[7], MPFR_RNDN), accessFloatContainerPointer(args[8]), mpfr_get_ui(args[9], MPFR_RNDN), mpfr_get_d(args[10], MPFR_RNDN), accessFloatContainerPointer(args[11]), mpfr_get_ui(args[12], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_3DT6PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DT6PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessIntContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_DTIPIPB (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTIPIPB(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessBoolContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DT5PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT5PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_DT6PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT6PI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_DT2PIPDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT2PIPDTPI(mpfr_get_ui(args[0], MPFR_RNDN), accessIntContainerPointer(args[1]), accessIntContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessIntContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_DT3PITDT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT3PITDT3PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_2IT2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2IT2PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_DT5PITDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DT5PITDT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), mpfr_get_ui(args[6], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2DTPI4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2DTPI4PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2DTPI5PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2DTPI5PIT(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessFloatContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT void * H_2DT2PIDT3PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT2PIDT3PI(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), accessIntContainerPointer(args[3]), mpfr_get_ui(args[4], MPFR_RNDN), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7])));
}

__MATHSUITE __JBURKARDT inline void * H_3I5PIDT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3I5PIDT2PI(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), mpfr_get_si(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessIntContainerPointer(args[6]), accessIntContainerPointer(args[7]), mpfr_get_ui(args[8], MPFR_RNDN), accessIntContainerPointer(args[9]), accessIntContainerPointer(args[10])));
}

__MATHSUITE __JBURKARDT inline void * H_PIPIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PIPIT(accessIntContainerPointer(args[0]), accessFloatContainerPointer(args[1])));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2DT2PIPITPIPITPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2DT2PIPITPIPITPI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessIntContainerPointer(args[9])));
}

__MATHSUITE __JBURKARDT inline void * H_2DTPII (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DTPII(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessIntContainerPointer(args[2]), mpfr_get_si(args[3], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTPIT2DT2PIPITPI3PIT2PI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTPIT2DT2PIPITPI3PIT2PI(mpfr_get_ui(args[0], MPFR_RNDN), accessFloatContainerPointer(args[1]), mpfr_get_ui(args[2], MPFR_RNDN), mpfr_get_ui(args[3], MPFR_RNDN), accessIntContainerPointer(args[4]), accessIntContainerPointer(args[5]), accessFloatContainerPointer(args[6]), accessIntContainerPointer(args[7]), accessFloatContainerPointer(args[8]), accessFloatContainerPointer(args[9]), accessFloatContainerPointer(args[10]), accessIntContainerPointer(args[11]), accessIntContainerPointer(args[12])));
}

__MATHSUITE __JBURKARDT inline void * H_PDTPI3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_PDTPI3PIT(accessUShortContainerPointer(args[0]), accessIntContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_2PDTPU (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PDTPU(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_3PDTPU (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3PDTPU(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessUShortContainerPointer(args[2]), accessIntContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2PDT4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PDT4PIT(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_2PDTPI (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PDTPI(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessIntContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_2PDTPS (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PDTPS(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessShortContainerPointer(args[2])));
}

__MATHSUITE __JBURKARDT inline void * H_2PDTPS4PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2PDTPS4PIT(accessUShortContainerPointer(args[0]), accessUShortContainerPointer(args[1]), accessShortContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4]), accessFloatContainerPointer(args[5]), accessFloatContainerPointer(args[6])));
}

__MATHSUITE __JBURKARDT inline void * H_2DT3PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2DT3PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), accessUShortContainerPointer(args[2]), accessUShortContainerPointer(args[3]), accessUShortContainerPointer(args[4])));
}

__MATHSUITE __JBURKARDT inline void * H_3DTPI2PDT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_3DTPI2PDT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_ui(args[1], MPFR_RNDN), mpfr_get_ui(args[2], MPFR_RNDN), accessIntContainerPointer(args[3]), accessUShortContainerPointer(args[4]), accessUShortContainerPointer(args[5])));
}

__MATHSUITE __JBURKARDT inline void * H_SIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_SIT(mpfr_get_si(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN)));
}

__MATHSUITE __JBURKARDT inline void * H_DTS2PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_DTS2PIT(mpfr_get_ui(args[0], MPFR_RNDN), mpfr_get_si(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3])));
}

__MATHSUITE __JBURKARDT inline void * H_2IT3PIT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[])
{
	return nodes->data.function.rfptr(C_2IT3PIT(mpfr_get_d(args[0], MPFR_RNDN), mpfr_get_d(args[1], MPFR_RNDN), accessFloatContainerPointer(args[2]), accessFloatContainerPointer(args[3]), accessFloatContainerPointer(args[4])));
}

