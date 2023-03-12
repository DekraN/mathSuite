#ifndef __DISABLEDEEP_FEM1DBVPQUADRATIC

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _fem1d_bvp_quadratic ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEM1D_BVP_QUADRATIC solves a two point boundary value problem.
  Discussion:
    The finite element method is used, with a mesh of piecewise quadratic
    elements.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of nodes.
    Input, double A ( double X ), evaluates a(x);
    Input, double C ( double X ), evaluates c(x);
    Input, double F ( double X ), evaluates f(x);
    Input, double X[N], the mesh points.
    Output, double FEM1D_BVP_QUADRATIC[N], the finite element coefficients,
    which are also the value of the computed solution at the mesh points.
*/
{
	const dt3fitpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp (* a)(ityp) = s_data->a1;
	ityp (* c)(ityp) = s_data->a2;
	ityp (* f)(ityp) = s_data->a3;
	ityp * x = s_data->a4;
	
    # define QUAD_NUM 3

    ityp abscissa[QUAD_NUM] =
    {
        -0.774596669241483377035853079956,
        0.000000000000000000000000000000,
        0.774596669241483377035853079956
    };
    ityp al;
    ityp am;
    ityp ar;
    ityp *amat;
    ityp axq;
    ityp *b;
    ityp bm;
    ityp cxq;
    dim_typ e;
    dim_typ e_num;
    ityp fxq;
    dim_typ i;
    dim_typ ierror;
    dim_typ j;
    dim_typ l;
    dim_typ m;
    dim_typ q;
    dim_typ quad_num = QUAD_NUM;
    dim_typ r;
    ityp *u;
    ityp weight[QUAD_NUM] =
    {
        0.555555555555555555555555555556,
        0.888888888888888888888888888889,
        0.555555555555555555555555555556
    };
    ityp wq;
    ityp vl;
    ityp vlp;
    ityp vm;
    ityp vmp;
    ityp vr;
    ityp vrp;
    ityp xl;
    ityp xm;
    ityp xq;
    ityp xr;
    /*
    Zero out the matrix and right hand side.
    */
    amat = r8mat_zero_new ( n, n );
    b = r8vec_zero_new ( n );
    /*
    Integrate over element E.
    */
    e_num = ( n - 1 ) / 2;

    for ( e = 0; e < e_num; ++e )
    {
        /*
        Element E uses nodes
        L = 2 * E
        M = 2 * E + 1
        R = 2 * E + 2
        */
        l = e<<1;
        m = (e<<1) + 1;
        r = (e<<1) + 2;

        xl = x[l];
        xm = x[m];
        xr = x[r];

        for ( q = 0; q < quad_num; ++q )
        {

            xq = ( ( 1.00 - abscissa[q] ) * xl+ ( 1.00 + abscissa[q] ) * xr )/   2.00;
            wq = weight[q] * ( xr - xl ) / 2.0;
            axq = a ( xq );
            cxq = c ( xq );
            fxq = f ( xq );
            vl = ( ( xq - xm ) / ( xl - xm ) )
            * ( ( xq - xr ) / ( xl - xr ) );

            vm = ( ( xq - xl ) / ( xm - xl ) )* ( ( xq - xr ) / ( xm - xr ) );

            vr = ( ( xq - xl ) / ( xr - xl ) )* ( ( xq - xm ) / ( xr - xm ) );

            vlp = (         1.00 / ( xl - xm ) )* ( ( xq - xr ) / ( xl - xr ) )+ ( ( xq - xm ) / ( xl - xm ) )* (         1.00 / ( xl - xr ) );

            vmp = (         1.00 / ( xm - xl ) )* ( ( xq - xr ) / ( xm - xr ) )+ ( ( xq - xl ) / ( xm - xl ) )* (         1.00 / ( xm - xr ) );

            vrp = (         1.00 / ( xr - xl ) )* ( ( xq - xm ) / ( xr - xm ) )+ ( ( xq - xl ) / ( xr - xl ) )* (         1.00 / ( xr - xm ) );

            amat[l+l*n] += wq * ( vlp * axq * vlp + vl * cxq * vl );
            amat[l+m*n] += wq * ( vlp * axq * vmp + vl * cxq * vm );
            amat[l+r*n] += wq * ( vlp * axq * vrp + vl * cxq * vr );
            b[l]   += wq * ( vl * fxq );

            amat[m+l*n] += wq * ( vmp * axq * vlp + vm * cxq * vl );
            amat[m+m*n] += wq * ( vmp * axq * vmp + vm * cxq * vm );
            amat[m+r*n] += wq * ( vmp * axq * vrp + vm * cxq * vr );
            b[m] += wq * ( vm * fxq );

            amat[r+l*n] += wq * ( vrp * axq * vlp + vr * cxq * vl );
            amat[r+m*n] += wq * ( vrp * axq * vmp + vr * cxq * vm );
            amat[r+r*n] += wq * ( vrp * axq * vrp + vr * cxq * vr );
            b[r] += wq * ( vr * fxq );
        }
    }
    /*
    Equation 0 is the left boundary condition, U(0.0) = 0.0;
    */
    i = 0;
    for ( j = 0; j < n; ++j )
        amat[i+j*n] = 0.00;
    amat[i+i*n] = 1.00;
    b[i] = 0.00;
    /*
    Equation N-1 is the right boundary condition, U(1.0) = 0.0;
    */
    i = n - 1;
    for ( j = 0; j < n; ++j )
        amat[i+j*n] = 0.00;
    amat[i+i*n] = 1.00;
    b[i] = 0.00;
    /*
    Solve the linear system.
    */
    u = r8mat_solve2 ( n, amat, b, &ierror );

    free ( amat );
    free ( b );
    return u;
    # undef QUAD_NUM
}

#endif
