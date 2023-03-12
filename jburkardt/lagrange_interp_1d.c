#ifndef __DISABLEDEEP_LAGRANGEINTERP1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_value_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_VALUE_1D evaluates the Lagrange interpolant.
  Discussion:
    The Lagrange interpolant L(ND,XD,YD)(X) is the unique polynomial of
    degree ND-1 which interpolates the points (XD(I),YD(I)) for I = 1
    to ND.
    The Lagrange interpolant can be constructed from the Lagrange basis
    polynomials.  Given ND distinct abscissas, XD(1:ND), the I-th Lagrange
    basis polynomial LB(ND,XD,I)(X) is defined as the polynomial of degree
    ND - 1 which is 1 at  XD(I) and 0 at the ND - 1 other abscissas.
    Given data values YD at each of the abscissas, the value of the
    Lagrange interpolant may be written as
      L(ND,XD,YD)(X) = sum ( 1 <= I <= ND ) LB(ND,XD,I)(X) * YD(I)
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
    Input, double YD[ND], the data values.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double LAGRANGE_VALUE_1D[NI], the interpolated values.
*/
{
	const _2dt3pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ ni = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	ityp * xi = s_data->a4;
	
    ityp *lb;
    ityp *yi;
    lb = lagrange_basis_1d ( nd, xd, ni, xi );
    yi = r8mat_mv_new ( ni, nd, lb, yd );
    free ( lb );
    return yi;
}

#endif
