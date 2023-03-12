#ifndef __DISABLEDEEP_VANDERMONDEINTERP2D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _vandermonde_interp_2d_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    VANDERMONDE_INTERP_2D_MATRIX computes a Vandermonde 2D interpolation matrix.
  Discussion:
    We assume the approximating function has the form of a polynomial
    in X and Y of total degree M.
      p(x,y) = c00
             + c10 * x                + c01 * y
             + c20 * x^2   + c11 * xy + c02 * y^2
             + ...
             + cm0 * x^(m) + ...      + c0m * y^m.
    If we let T(K) = the K-th triangular number
            = sum ( 1 <= I <= K ) I
    then the number of coefficients in the above polynomial is T(M+1).
    We have n data locations (x(i),y(i)) and values z(i) to approximate:
      p(x(i),y(i)) = z(i)
    and we assume that N = T(M+1).
    This can be cast as an NxN linear system for the polynomial
    coefficients:
      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
      [ ...................... ] [ ... ] = [ ... ]
      [ 1 xn yn  xn^2 ... yn^m ] [ c0n ] = [  zn ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data points.  It is necessary
    that N = T(M+1), where T(K) is the K-th triangular number.
    Input, int M, the degree of the polynomial.
    Input, double X[N], Y[N], the data locations.
    Output, double VANDERMONDE_INTERP_2D_MATRIX[N*N], the Vandermonde matrix for X.
*/
{
	const _2dt2pit * const s_data = data; 
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * x = s_data->a2;
	ityp * y = s_data->a3;
	
    ityp *a;
    dim_typ ex;
    int ey;
    dim_typ i;
    dim_typ j;
    dim_typ s;

    if ( n != triangle_num ( m + 1 ) )
        return NULL;

    a = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
    j = 0;

    for ( s = 0; s <= m; ++s )
        for ( ex = s; 0 <= ex; --ex )
        {
            ey = s - ex;
            for ( i = 0; i < n; ++i )
                a[i+j*n] = pow ( x[i], ex ) * pow ( y[i], ey );
            ++ j;
        }

    return a;
}

#endif
