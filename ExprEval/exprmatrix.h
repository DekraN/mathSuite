#ifndef WRAPPER_EXPRMATRIX_H_INCLUDED
#define WRAPPER_EXPRMATRIX_H_INCLUDED

__MATHSUITE  void * FUNCNAME_MATRIXDET(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE  void * FUNCNAME_MATRIXNORM(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE  void * FUNCNAME_MATRIXTRACE(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE  void * FUNCNAME_MATRIXRANK(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE  void * FUNCNAME_MATRIXILLCHK(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE  void * FUNCNAME_SCALARPROD(EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);

#endif
