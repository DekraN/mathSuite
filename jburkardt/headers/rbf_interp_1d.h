#ifndef WRAPPER_RBF_INTERP_1D_H_INCLUDED
#define WRAPPER_RBF_INTERP_1D_H_INCLUDED

#define phi1(a,b,c,d) FUNCNAME_PHI1(C_DTIT2PIT(a,c,b,d))
#define phi2(a,b,c,d) FUNCNAME_PHI2(C_DTIT2PIT(a,c,b,d))
#define phi3(a,b,c,d) FUNCNAME_PHI3(C_DTIT2PIT(a,c,b,d))
#define phi4(a,b,c,d) FUNCNAME_PHI4(C_DTIT2PIT(a,c,b,d))
#define r8mat_solve_svd(a,b,c,d) FUNCNAME_R8MATSOLVESVD(C_2DT2PIT(a,b,c,d))

__MATHSUITE __JBURKARDT void * FUNCNAME_PHI1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PHI2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PHI3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PHI4(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATSOLVESVD(void *);

__MATHSUITE __JBURKARDT  ityp   *rbf_interp ( const register dim_typ m, const register dim_typ nd, ityp [static m*nd], const register ityp,
  void phix ( int, ityp [], ityp , ityp [] ), ityp [static nd], const register dim_typ ni, ityp [static m*ni] );
__MATHSUITE __JBURKARDT  ityp   *rbf_weight ( const register dim_typ m, const register dim_typ nd, ityp [static m*nd], const register ityp,void ( int , ityp [], ityp, ityp [] ),ityp [static nd] );

#endif
