#ifndef __DISABLEDEEP_TETRAHEDRONINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON01_MONOMIAL_INTEGRAL: integrals in the unit tetrahedron in 3D.
  Discussion:
    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 January 2014
  Author:
    John Burkardt
  Parameters:
    Input, int E[3], the exponents.
    Each exponent must be nonnegative.
    Output, double TETRAHEDRON01_MONOMIAL_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	int * e = data;
	
    ityp integral;
    dim_typ i, j, k;
    const register dim_typ m = 3;

    for ( i = 0; i < m; ++i)
        if ( e[i] < 0 )
        {
        	result = MAX_VAL;
            return &result;
        }

    k = 00;
    integral = 1.00;

    for ( i = 0; i < m; ++i)
        for ( j = 1; j <= e[i]; ++j)
        {
            ++ k;
            integral *= ( ityp ) ( j ) / ( ityp ) ( k );
        }

    for ( i = 0; i < m; ++i )
    {
        ++ k;
        integral /= ( ityp ) ( k );
    }

	result = integral;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tetrahedron01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON01_SAMPLE samples the unit tetrahedron in 3D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 January 2014
  Author:
    John Burkardt
  Reference:
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
    Output, double TETRAHEDRON01_SAMPLE[3*N], the points.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    ityp *e;
    ityp e_sum;
    dim_typ i, j;
    const register dim_typ m = 3;
    ityp *x = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = 0; j < n; j++ )
    {
        e = r8vec_uniform_01_new ( m + 1, seed );

        for ( i = 0; i < m + 1; ++i )
            e[i] = - log ( e[i] );
        e_sum = r8vec_sum ( m + 1, e );

        for ( i = 0; i < m; ++i )
            x[i+j*m] = e[i] / e_sum;
        free ( e );
    }

    return x;
}

#endif
