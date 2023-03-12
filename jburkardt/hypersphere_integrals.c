#ifndef __DISABLEDEEP_HYPERSPHEREINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hypersphere01_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE01_AREA returns the surface area of the unit hypersphere.
  Discussion:
     M   Area
     2    2        * M_PI
     3    4        * M_PI
     4 ( 2 /   1) * M_PI^2
     5 ( 8 /   3) * M_PI^2
     6 ( 1 /   1) * M_PI^3
     7 (16 /  15) * M_PI^3
     8 ( 1 /   3) * M_PI^4
     9 (32 / 105) * M_PI^4
    10 ( 1 /  12) * M_PI^5
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 January 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Output, double HYPERSPHERE01_AREA, the area.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ m = *(dim_typ *) data;
	
    ityp area;
    dim_typ i, m_half;

    if ( ( m % 2 ) == 0 )
    {
        m_half = m / 2;
        area = 2.00 * pow ( M_PI, m_half );
        for ( i = 1; i <= m_half - 1; ++i )
            area /= ( ityp ) ( i );
    }
    else
    {
        m_half = ( m - 1 ) / 2;
        area = pow ( M_PI, m_half ) * pow ( 2.00, m );
        for ( i = m_half + 1; i <= (m_half<<1); ++i )
            area /= ( ityp ) ( i );
    }

	result = area;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hypersphere01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE01_MONOMIAL_INTEGRAL: monomial integrals on the unit hypersphere.
  Discussion:
    The integration region is
      sum ( 1 <= I <= M ) X(I)^2 = 1.
    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 January 2014
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Academic Press, 1984, page 263.
  Parameters:
    Input, int M, the spatial dimension.
    Input, int E[M], the exponents.  Each exponent must be nonnegative.
    Output, double HYPERSPHERE01_MONOMIAL_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	int * e = s_data->a1;
	
    ityp arg;
    dim_typ i;
    ityp integral;

    for ( i = 0; i < m; ++i)
        if ( e[i] < 0 )
        {
        	result = MAX_VAL;
            return &result;
        }

    for ( i = 0; i < m; ++i )
        if ( ( e[i] % 2 ) == 1 )
        {
        	result = 0.00;
            return &result;
        }

    integral = 2.00;

    for ( i = 0; i < m; ++i )
    {
        arg = 0.50 * ( ityp ) ( e[i] + 1 );
        integral *= r8_gamma ( arg );
    }
    
    result = integral / r8_gamma ( 0.50 * ( ityp ) ( i4vec_sum ( m, e ) + m ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hypersphere01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE01_SAMPLE uniformly samples the surface of the unit hypersphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 January 2014
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
    Input, int M, the spatial dimension.
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double X[M*N], the points.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
    dim_typ i, j;
    ityp norm;
    ityp *x;

    x = r8mat_normal_01_new ( m, n, seed );

    for ( j = 0; j < n; ++j )
    {
        /*
        Compute the length of the vector.
        */
        norm = 0.00;
        for ( i = 0; i < m; ++i )
            norm += pow ( x[i+j*m], 2 );
        norm = sqrt ( norm );
        /*
        Normalize the vector.
        */
        for ( i = 0; i < m; ++i)
            x[i+j*m] /= norm;
    }
    return x;
}

#endif
