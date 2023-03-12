#ifndef __DISABLEDEEP_VANDERMONDEAPPROX1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _vandermonde_approx_1d_coef ( void * data)
/******************************************************************************/
/*
  Purpose:
    VANDERMONDE_APPROX_1D_COEF computes a 1D polynomial approximant.
  Discussion:
    We assume the approximating function has the form
      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
    We have n data values (x(i),y(i)) which must be approximated:
      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
    This can be cast as an Nx(M+1) linear system for the polynomial
    coefficients:
      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
      [ .................. ] [ ... ] = [ ... ]
      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
    In the typical case, N is greater than M+1 (we have more data and equations
    than degrees of freedom) and so a least squares solution is appropriate,
    in which case the computed polynomial will be a least squares approximant
    to the data.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data points.
    Input, int M, the degree of the polynomial.
    Input, double X[N], Y[N], the data values.
    Output, double VANDERMONDE_APPROX_1D_COEF[M+1], the coefficients of
    the approximating polynomial.  C(0) is the constant term, and C(M)
    multiplies X^M.
*/
{
	const _2dt2pit * const s_data = data; 
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * x = s_data->a2;
	ityp * y = s_data->a3;
	
    ityp *a = vandermonde_approx_1d_matrix ( n, m, x );
    ityp *c = qr_solve ( n, m + 1, a, y );
    free ( a );
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _vandermonde_approx_1d_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    VANDERMONDE_APPROX_1D_MATRIX computes a Vandermonde 1D approximation matrix.
  Discussion:
    We assume the approximant has the form
      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
    We have n data values (x(i),y(i)) which must be approximated:
      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
    This can be cast as an Nx(M+1) linear system for the polynomial
    coefficients:
      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
      [ .................. ] [ ... ] = [ ... ]
      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data points.
    Input, int M, the degree of the polynomial.
    Input, double X(N), the data values.
    Output, double VANDERMONDE_APPROX_1D_MATRIX[N*(M+1)], the Vandermonde matrix for X.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * x = s_data->a2;
	
    dim_typ i, j;
    ityp *a = ( ityp * ) malloc ( n * ( m + 1 ) * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        a[i] = 1.00;

    for ( j = 1; j <= m; ++j )
        for ( i = 0; i < n; ++i)
            a[i+j*n] = a[i+(j-1)*n] * x[i];

    return a;
}

#endif
