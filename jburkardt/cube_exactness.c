#ifndef __DISABLEDEEP_CUBEEXACTNESS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _legendre_3d_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_3D_MONOMIAL_INTEGRAL the Legendre integral of a monomial.
  Discussion:
    The Legendre integral to be evaluated has the form
      I(f) = integral ( z1 <= z <= z2 )
             integral ( y1 <= y <= y2 )
             integral ( x1 <= x <= x2 ) x^i y^j z^k dx dy dz
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A[3], the lower limits of integration.
    Input, double B[3], the upper limits of integration.

    Input, int P[3], the exponents of X and Y.

    Output, double LEGENDRE_3D_MONOMIAL_INTEGRAL, the value of the exact integral.
*/
{
	static ityp result = MAX_VAL;
	
	const pdt2pit * const s_data = data;
	
	dim_typ * p = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
	result = ( pow ( b[0], p[0] + 1 ) - pow ( a[0], p[0] + 1 ) ) / ( ityp ) ( p[0] + 1 )
        * ( pow ( b[1], p[1] + 1 ) - pow ( a[1], p[1] + 1 ) ) / ( ityp ) ( p[1] + 1 )
        * ( pow ( b[2], p[2] + 1 ) - pow ( a[2], p[2] + 1 ) ) / ( ityp ) ( p[2] + 1 );
    return &result;
}

#endif
