#ifndef __DISABLEDEEP_TETRAHEDRONEXACTNESS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tet01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET01_MONOMIAL_INTEGRAL integrates a monomial over the unit tetrahedron.
  Discussion:
    This routine evaluates a monomial of the form
      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
    where the exponents are nonnegative integers.  Note that
    if the combination 0^0 is encountered, it should be treated
    as 1.
    Integral ( over unit tetrahedron ) x^l y^m z^n dx dy
    = l! * m! * n! / ( l + m + n + 3 )!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 July 2007
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int EXPON[DIM_NUM], the exponents.
    Output, double TET01_MONOMIAL_INTEGRAL, the value of the integral of the
    monomial.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	int * expon = s_data->a1;
	
    dim_typ i, k = 0;
    /*
    The first computation ends with VALUE = 1.0;
    */
    ityp value = 1.00;


    for ( i = 1; i <= expon[0]; ++i,++k);

    for ( i = 1; i <= expon[1]; ++i )
    {
        ++ k;
        value *= ( ityp ) ( i ) / ( ityp ) ( k );
    }

    for ( i = 1; i <= expon[2]; ++i )
    {
        ++ k;
        value *= ( ityp ) ( i ) / ( ityp ) ( k );
    }

    ++ k;
    value /= ( ityp ) ( k );

    ++ k;
    value /= ( ityp ) ( k );

    ++ k;
    value /= ( ityp ) ( k );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tet01_monomial_quadrature ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET01_MONOMIAL_QUADRATURE applies quadrature to a monomial in a tetrahedron.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 July 2007
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int EXPON[DIM_NUM], the exponents.
    Input, int POINT_NUM, the number of points in the rule.
    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
    Input, double WEIGHT[POINT_NUM], the quadrature weights.
    Output, double MONOMIAL_QUADRATURE, the quadrature error.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitdtpipit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	ityp * expon = s_data->a1;
	const register dim_typ point_num = s_data->a2;
	int * x = s_data->a3;
	ityp * weight = s_data->a4;
	
    ityp exact;
    ityp quad;
    ityp  quad_error;
    ityp *value;
    /*
    Get the exact value of the integral of the unscaled monomial.
    */
    /*
    Evaluate the monomial at the quadrature points.
    */
    value = monomial_value ( dim_num, point_num, expon, x );
    /*
    Compute the weighted sum and divide by the exact value.
    */
    quad = (1.00 / 6.00) * r8vec_dot_product ( point_num, weight, value ) / tet01_monomial_integral ( dim_num, (int*)expon );
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
