#ifndef WRAPPER_FEM2D_SERENE_H_INCLUDED
#define WRAPPER_FEM2D_SERENE_H_INCLUDED

#define basis_serene(a,b,c,d,e,f,g,h) FUNCNAME_BASISSERENE(C_6IT2PIT(a,b,c,d,e,f,g,h))
#define basis_dx_serene(a,b,c,d,e,f,g,h) FUNCNAME_BASISDXSERENE(C_6IT2PIT(a,b,c,d,e,f,g,h))
#define basis_dy_serene(a,b,c,d,e,f,g,h) FUNCNAME_BASISDYSERENE(C_6IT2PIT(a,b,c,d,e,f,g,h))
#define fem2d_bvp_serene_node_num(a,b) R_DBL(FUNCNAME_FEM2DBVPSERENENODENUM(C_PUSHRT2(a,b)))
#define not1(a,b,c) R_DBL(FUNCNAME_NOT1(C_PDBL3(a,b,c)))
#define not1d(a,b) R_DBL(FUNCNAME_NOT1D(C_PDBL2(a,b)))
#define not2(a,b,c,d,e,f,g,h) R_DBL(FUNCNAME_NOT2(C_PDBL8(a,b,c,d,e,f,g,h)))
#define not2dx(a,b,c,d,e,f) R_DBL(FUNCNAME_NOT2DX(C_PDBL6(a,b,c,d,e,f)))
#define not2dy(a,b,c,d,e,f) R_DBL(FUNCNAME_NOT2DY(C_PDBL6(a,b,c,d,e,f)))


__MATHSUITE __JBURKARDT void * FUNCNAME_BASISSERENE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BASISDXSERENE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BASISDYSERENE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FEM2DBVPSERENENODENUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NOT1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NOT1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NOT2DX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NOT2DY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NOT2(void *);

#endif

