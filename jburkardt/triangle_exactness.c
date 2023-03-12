#ifndef __DISABLEDEEP_TRIANGLEEXACTNESS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE01_MONOMIAL_INTEGRAL integrates a monomial over the unit triangle.
  Discussion:
    This routine evaluates a monomial of the form
      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
    where the exponents are nonnegative integers.  Note that
    if the combination 0^0 is encountered, it should be treated
    as 1.
    Integral ( over unit triangle ) x^m y^n dx dy = m! * n! / ( m + n + 2 )!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 July 2007
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int EXPON[DIM_NUM], the exponents.
    Output, double TRIANGLE01_MONOMIAL_INTEGRAL, the value of the integral
    of the monomial.
*/
{	
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	int * expon = s_data->a1;
	
    /*
    The first computation ends with VALUE = 1.0;
    */
    ityp value = 1.00;
    dim_typ i, k = 0;

    for ( i = 1; i <= expon[0]; i++, ++k);

    for ( i = 1; i <= expon[1]; ++i)
    {
        ++ k;
        value *= ( ityp ) ( i ) / ( ityp ) ( k );
    }

    ++ k;
    value /= ( ityp ) ( k );

    ++ k;
    value /= ( ityp ) ( k );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle01_monomial_quadrature ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRANGLE01_MONOMIAL_QUADRATURE applies quadrature to a monomial in a triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 July 2007
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int EXPON[DIM_NUM], the exponents.
    Input, int POINT_NUM, the number of points in the rule.
    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
    Input, double WEIGHT[POINT_NUM], the quadrature weights.
    Output, double TRIANGLE01_MONOMIAL_QUADRATURE, the quadrature error.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitdtpipit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	ityp * expon = s_data->a1;
	const register dim_typ point_num = s_data->a2;
	int * x = s_data->a3;
	ityp * weight = s_data->a4;
	
    ityp area;
    ityp exact;
    ityp quad;
    ityp quad_error;
    ityp scale;
    ityp *value;
    /*
    Get the exact value of the integral of the unscaled monomial.
    */
    scale = triangle01_monomial_integral ( dim_num, (int*)expon );
    /*
    Evaluate the monomial at the quadrature points.
    */
    value = monomial_value ( dim_num, point_num, expon, x );
    /*
    Compute the weighted sum and divide by the exact value.
    */
    area = 0.50;
    quad = area * r8vec_dot_product ( point_num, weight, value ) / scale;
    /*
    Error:
    */
    exact = 1.00;
    quad_error = fabs ( quad - exact );

    free ( value );
	
	result = quad_error;
    return &result;
}

#endif
