#ifndef FEM2D_BVP_LINEAR_H_INCLUDED
#define FEM2D_BVP_LINEAR_H_INCLUDED

__MATHSUITE __JBURKARDT ityp   *fem2d_bvp_linear (const register dim_typ nx, const register dim_typ ny, ityp ( ityp, ityp ),ityp ( ityp, ityp ), ityp ( ityp, ityp ),ityp [static nx], ityp [static ny] );
__MATHSUITE __JBURKARDT ityp   fem2d_h1s_error_linear (const register dim_typ nx, const register dim_typ ny, ityp [static nx], ityp [static ny],
  ityp [static nx*ny], ityp ( ityp, ityp ),ityp ( ityp, ityp ) );
__MATHSUITE __JBURKARDT ityp   fem2d_l1_error (const register dim_typ nx, const register dim_typ ny, ityp [static nx], ityp [static ny], ityp [static nx*ny],ityp ( ityp, ityp ) );
__MATHSUITE __JBURKARDT ityp   fem2d_l2_error_linear (const register dim_typ nx, const register dim_typ ny, ityp [static nx], ityp [static ny],ityp [static nx*ny], ityp ( ityp, ityp ) );

#endif // FEM2D_BVP_LINEAR_H_INCLUDED
