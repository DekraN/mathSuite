#ifndef __DISABLEDEEP_LINECVTLLOYD

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_ccvt_lloyd_step ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_CCVT_LLOYD_STEP takes one step of Lloyd"s constrained CVT algorithm.
  Discussion:
    Each step of Lloyd"s algorithm replaces a point by the center of mass
    of the associated region.  For points on a line, with a uniform
    density, the associated region is demarcated by the midways between
    successive points.
    Here, we include the additional constraint that we want the first and last
    points to be fixed at the endpoints of the line, that is, X(1) = A
    and X(2) = B.  In that case, the calculation of the updates for the
    first two and last two points must be handled differently.
    For points away from the boundary, a step of Lloyd"s method can be
    regarded as replacing each point by the average of the left and right
    midways.  The midways, of course, are the average of two points.
    So for point J, we have:
      M(J-1,J) = ( X(J-1) + X(J) ) / 2
      M(J,J+1) = ( X(J) + X(J+1) ) / 2
      X*(J) = ( M(J-1,J) + M(J,J+1) ) / 2 = ( X(J-1) + 2 X(J) + X(J+1) ) / 4
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    1 <= N.
    Input, double A, B, the left and right endpoints.
    Input/output, double X[N], the point locations.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	ityp * x = s_data->a3;
	
    dim_typ j;
    ityp *x_old = r8vec_copy_new ( n, x );

    if ( n == 1 )
        x[0] = ( a + b ) / 2.00;
    else if ( n == 2 )
    {
        x[0] = a;
        x[1] = b;
    }
    else
    {
        x[0] = a;

        for ( j = 1; j < n - 1; ++j )
            x[j] = ( 0.50 * ( x_old[j-1] + x_old[j] ) + 0.50 * ( x_old[j] + x_old[j+1] ) ) / 2.00;
        x[n-1] = b;
    }

    free ( x_old );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_cvt_energy ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_CVT_ENERGY computes the CVT energy for a given set of generators.
  Discussion:
    Given a set of generators G over the line [A,B], then the energy
    is defined as
      E = integral ( a <= x <= b ) ( x - g(x) )^2 dx
    where g(x) is the nearest generator to the point x.
    For the 1D case, this integral can be evaluated exactly as the
    sum of integrals over each subinterval:
      E(i) = integral ( xl <= x <= xr ) ( x - x(i) )^2 dx
           = ( ( x(i) - xl )^3 + ( xr - x(i) )^3 ) / 3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of generators.
    Input, double A, B, the left and right endpoints.
    Input, double X[N], the generator locations.
    Output, double LINE_CVT_ENERGY, the energy of the generator distribution.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	ityp * x = s_data->a3;
	
    dim_typ j;
    ityp xl;
    ityp xr;
    ityp e = 0.00;

    for ( j = 0; j < n; ++j)
    {
        xl = j == 0 ? a : ( x[j-1] + x[j] ) / 2.00;
        xr = j == n-1 ? b : ( x[j] + x[j+1] ) / 2.00;
        e += ( pow ( x[j] - xl, 3 ) + pow ( xr - x[j], 3 )  ) / 3.0;
    }

	result = e;
    return &result;
}

#endif
