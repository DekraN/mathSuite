#ifndef __DISABLEDEEP_HYPERBALLINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_normal_01_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_NORMAL_01_NEW returns a unit pseudonormal R8MAT.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    03 October 2005
  Author:
    John Burkardt
  Reference:
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.
    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.
    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.
  Parameters:
    Input, int M, N, the number of rows and columns in the array.
    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.
    Output, double R8MAT_NORMAL_01_NEW[M*N], the array of pseudonormal values.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
	return r8vec_normal_01_new ( m * n, seed );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hyperball01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERBALL01_MONOMIAL_INTEGRAL: integrals in the unit hyperball in M dimensions.
  Discussion:
    The integration region is
      sum ( 1 <= I <= M ) X(I)^2 <= 1.
    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 January 2014
  Author:
    John Burkardt
  Reference:
    Philip Dav__MATHSUITE __JBURKARDT  ityp  is, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Academic Press, 1984, page 263.
  Parameters:
    Input, int M, the spatial dimension.
    Input, int E[M], the exponents.  Each exponent must be nonnegative.
    Output, double HYPERBALL01_MONOMIAL_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	int * e = s_data->a1;
	
    dim_typ arg;
    dim_typ i;
    dim_typ integral;
    const dim_typ r = 1.00;
    dim_typ s;

    for ( i = 0; i < m; ++i )
        if ( e[i] < 0 )
        {
        	result = MAX_VAL;
            return &result;
        }

    for ( i = 0; i < m; ++i )
        if ( ( e[i] % 2 ) == 1 )
        {
        	result = MAX_VAL;
            return &result;
        }

    integral = 2.00;

    for ( i = 0; i < m; ++i )
    {
        arg = 0.50 * ( ityp ) ( e[i] + 1 );
        integral *= r8_gamma ( arg );
    }

    s = i4vec_sum ( m, e ) + m;
    arg = 0.50 * ( dim_typ ) ( s );
    integral /= r8_gamma ( arg );
    /*
    The surface integral is now adjusted to give the volume integral.
    */

	result = integral*pow ( r, s  ) / ( dim_typ ) ( s );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hyperball01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERBALL01_SAMPLE samples points from the unit hyperball in M dimensions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 January 2014
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
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double X[M*N], the points.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
    ityp exponent;
    dim_typ i, j;
    ityp norm;
    ityp r;
    ityp *x;

    exponent = 1.00 / ( ityp ) ( m );

    x = r8mat_normal_01_new ( m, n, seed );

    for ( j = 0; j < n; ++j)
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
        for ( i = 0; i < m; ++i )
            x[i+j*m] /= norm;

        /*
        Now compute a value to map the point ON the sphere INTO the sphere.
        */
        r = r8_uniform_01 ( seed );

        for ( i = 0; i < m; ++i )
            x[i+j*m] = pow ( r, exponent ) * x[i+j*m];
    }
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hyperball01_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERBALL01_VOLUME returns the volume of the unit hyperball in M dimensions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 January 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Output, double HYPERBALL01_VOLUME, the volume.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ m = *(dim_typ *) data;
	
    dim_typ i, m_half;
    ityp volume;

    if ( ( m % 2 ) == 0 )
    {
        m_half = m / 2;
        volume = pow ( M_PI, m_half );
        for ( i = 1; i <= m_half; ++i )
            volume /= ( ityp ) ( i );
    }
    else
    {
        m_half = ( m - 1 ) / 2;
        volume = pow ( M_PI, m_half ) * pow ( 2.00, m );
        for ( i = m_half + 1; i <= (m_half<<1) + 1; ++i)
        volume /= ( ityp ) ( i );
    }

	result = volume;
    return &result;
}

#endif
