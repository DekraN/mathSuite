#ifndef __DISABLEDEEP_LINEINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line01_length ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE01_LENGTH: length of the unit line in 1D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2014
  Author:
    John Burkardt
  Parameters:
    Output, double LINE01_LENGTH, the length.
*/
{
	static ityp result = 1.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE01_MONOMIAL_INTEGRAL: integrals on the unit line in 1D.
  Discussion:
    The integration region is
      0 <= X <= 1.
    The monomial is F(X) = X^E.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2014
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Academic Press, 1984, page 263.
  Parameters:
    Input, int E, the exponent.  E must be nonnegative.
    Output, double LINE01_MONOMIAL_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ e = *(dim_typ *) data;
	
	result = e == -1 ? MAX_VAL : 1.00 / ( ityp ) ( e + 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _line01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE01_SAMPLE samples points on the unit line in 1D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2014
  Author:
    John Burkardt
  Reference:
    Russell Cheng,
    Random Variate Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998, pages 168.
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.
  Parameters:
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double X[N], the points.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1; 
	
    return r8vec_uniform_01_new ( n, seed );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _monomial_value_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONOMIAL_VALUE_1D evaluates a monomial in 1D.
  Discussion:
    This routine evaluates a monomial of the form
      x^e
    where the exponent is a nonnegative integer.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points at which the
    monomial is to be evaluated.
    Input, int E, the exponent.
    Input, double X[M*N], the point coordinates.
    Output, double MONOMIAL_VALUE_1D[N], the value of the monomial.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ e = s_data->a1;
	ityp * x = s_data->a2;
	
    ityp *v = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ j = 0; j < n; ++j )
        v[j] = pow ( x[j], e );
    return v;
}

#endif
