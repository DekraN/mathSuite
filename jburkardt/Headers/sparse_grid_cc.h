#ifndef WRAPPER_SPARSE_GRID_CC_H_INCLUDED
#define WRAPPER_SPARSE_GRID_CC_H_INCLUDED

#define i4_mop(a) R_SHRT(FUNCNAME_I4MOP(C_SUSHRT(a)))
#define cc_abscissa(a,b) R_DBL(FUNCNAME_CCABSCISSA(C_PUSHRT2(a,b)))
#define cc_weights(a) FUNCNAME_CCWEIGHTS(C_SUSHRT(a))
#define abscissa_level_closed_nd(a,b,c,d) FUNCNAME_ABSCISSALEVELCLOSEDND(C_3DTPI(a,b,c,d))
#define index_to_level_closed(a,b,c,d) R_INT(FUNCNAME_INDEXTOLEVELCLOSED(C_3DTPI(a,c,d,b)))
#define level_to_order_closed(a,b,c) FUNCNAME_LEVELTOORDERCLOSED(C_DT2PI(a,b,c))
#define multigrid_index0(a,b,c) FUNCNAME_MULTIGRIDINDEX0(C_2DTPI(a,c,b))
#define multigrid_scale_closed(a,b,c,d,e) FUNCNAME_MULTIGRIDSCALECLOSED(C_3DT2PI(a,b,c,d,e))
#define product_weights_cc(a,b,c) FUNCNAME_PRODUCTWEIGHTSCC(C_2DTPI(a,c,b))
#define sparse_grid_cc(a,b,c,d,e) FUNCNAME_SPARSEGRIDCC(C_3DT2PIT(a,b,c,d,e))
#define sparse_grid_cc_index(a,b,c) FUNCNAME_SPARSEGRIDCCINDEX(C_PUSHRT3(a,b,c))
#define sparse_grid_cc_size_old(a,b) FUNCNAME_SPARSEGRIDCCSIZEOLD(C_PUSHRT2(a,b))
#define sparse_grid_cc_weights(a,b,c,d,e) FUNCNAME_SPARSEGRIDCCWEIGHTS(C_3DTPIPIT(a,b,c,d,e))
#define sparse_grid_ccs_size(a,b) FUNCNAME_SPARSEGRIDCCSSIZE(C_PUSHRT2(a,b))
#define vec_colex_next2(a,b,c,d) FUNCNAME_VECCOLEXNEXT2(C_DT2PIPB(a,b,c,d))

__MATHSUITE __JBURKARDT void * FUNCNAME_I4MOP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_INDEXTOLEVELCLOSED(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LEVELTOORDERCLOSED(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MULTIGRIDSCALECLOSED(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPARSEGRIDCC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPARSEGRIDCCWEIGHTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECCOLEXNEXT2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CCWEIGHTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ABSCISSALEVELCLOSEDND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MULTIGRIDINDEX0(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PRODUCTWEIGHTSCC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPARSEGRIDCCINDEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPARSEGRIDCCSIZEOLD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPARSEGRIDCCSSIZE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CCABSCISSA(void *);

#endif