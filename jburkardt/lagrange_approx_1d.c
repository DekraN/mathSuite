#ifndef __DISABLEDEEP_LAGRANGEAPPROX1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_approx_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_APPROX_1D evaluates the Lagrange approximant of degree M.
  Discussion:
    The Lagrange approximant L(M,ND,XD,YD)(X) is a polynomial of
    degree M which approximates the data (XD(I),YD(I)) for I = 1 to ND.
    We can represent any polynomial of degree M+1 as the sum of the Lagrange
    basis functions at the M+1 Chebyshev points.
      L(M)(X) = sum ( 1 <= I <= M+1 ) C(I) LB(M,XC)(X)
    Given our data, we can seek the M+1 unknown coefficients C which minimize
    the norm of || L(M)(XD(1:ND)) - YD(1:ND) ||.
    Given the coefficients, we can then evaluate the polynomial at the
    points XI.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the polynomial degree.
    Input, int ND, the number of data points.
    ND must be at least 1.
    Input, double XD[ND], the data points.
    Input, double YD[ND], the data values.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double LAGRANGE_APPROX_1D[NI], the interpolated values.
*/
{
	const _3dt3pit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ nd = s_data->a1;
	const register dim_typ ni = s_data->a2;
	ityp * xd = s_data->a3;
	ityp * yd = s_data->a4;
	ityp * xi = s_data->a5; 
	
    ityp a;
    ityp b;
    ityp *ld;
    ityp *li;
    dim_typ nc;
    ityp *xc;
    ityp *yc;
    ityp *yi;

    nc = m + 1;
    /*
    Evaluate the Chebyshev points.
    */
    a = -1.00;
    b = +1.00;
    // WARNING!!
    // xc = r8vec_chebyspace_new ( nc, a, b );
    /*
    Evaluate the Lagrange basis functions for the Chebyshev points
    at the data points.
    */
    ld = lagrange_basis_1d ( nc, xc, nd, xd );
    /*
    The value of the Lagrange approximant at each data point should
    approximate the data value: LD * YC = YD, where YC are the unknown
    coefficients.
    */
    yc = qr_solve ( nd, nc, ld, yd );
    /*
    Now we want to evaluate the Lagrange approximant at the "interpolant
    points": LI * YC = YI
    */
    li = lagrange_basis_1d ( nc, xc, ni, xi );

    yi = r8mat_mv_new ( ni, nc, li, yc );

    free ( ld );
    free ( li );
    free ( xc );
    free ( yc );

    return yi;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_basis_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_BASIS_1D evaluates the Lagrange basis polynomials.
  Discussion:
    Given ND distinct abscissas, XD(1:ND),
    the I-th Lagrange basis polynomial LB(I)(T) is defined as the polynomial of
    degree ND - 1 which is 1 at XD(I) and 0 at the ND - 1
    other abscissas.
    A formal representation is:
      LB(I)(T) = Product ( 1 <= J <= ND, I /= J )
 ( T - T(J) ) / ( T(I) - T(J) )
    This routine accepts a set of NI values at which all the Lagrange
    basis polynomials should be evaluated.
    Given data values YD at each of the abscissas, the value of the
    Lagrange interpolating polynomial at each of the interpolation points
    is then simple to compute by matrix multiplication:
      YI(1:NI) = LB(1:NI,1:ND) * YD(1:ND)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data points.
    ND must be at least 1.
    Input, double XD[ND], the data points.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double LAGRANGE_BASIS[NI*ND], the values
    of the Lagrange basis polynomials at the interpolation points.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ ni = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * xi = s_data->a3;
	
    dim_typ i, j, k;
    ityp *lb;
    /*
    Evaluate the polynomial.
    */
    lb = ( ityp * ) malloc ( ni * nd * sizeof ( ityp ) );

    for ( j = 0; j < nd; ++j )
        for( i = 0; i < ni; ++i )
            lb[i+j*ni] = 1.00;

    for ( i = 0; i < nd; ++i )
        for ( j = 0; j < nd; ++j )
            if ( j != i )
                for ( k = 0; k < ni; ++k )
                    lb[k+i*ni] *= ( xi[k] - xd[j] ) / ( xd[i] - xd[j] );

    return lb;
}

#endif
