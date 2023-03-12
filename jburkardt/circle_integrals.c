#ifndef __DISABLEDEEP_CIRCLEINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle01_length ( void * data )
/******************************************************************************/
/*
  Purpose:
    CIRCLE01_LENGTH: length of the circumference of the unit circle in 2D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2014
  Author:
    John Burkardt
  Parameters:
    Output, double CIRCLE01_LENGTH, the length.
*/
{
	static ityp result = (M_2TPI);
    return &result;
}

#define CIRCLE01MONINT_NEGEXPERROR MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE01_MONOMIAL_INTEGRAL: integrals on circumference of unit circle in 2D.
  Discussion:
    The integration region is
      X^2 + Y^2 = 1.
    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2014
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Academic Press, 1984, page 263.
  Parameters:
    Input, int E[2], the exponents of X and Y in the
    monomial.  Each exponent must be nonnegative.
    Output, double CIRCLE01_MONOMIAL_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * e = data;
	
    dim_typ i;
    ityp integral;

    if ( ( e[0] % 2 ) == 1 || ( e[1] % 2 ) == 1 )
        integral = 0.00;
    else
    {
        integral = 2.00;

        #pragma omp parallel for
        for ( i = 0; i < 2; ++i )
            integral *= r8_gamma ( 0.50 * ( ityp ) ( e[i] + 1 ) );

        integral /= r8_gamma ( 0.50 * ( ityp ) ( e[0] + e[1] + 2 ) );
    }
    
    result = integral;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _circle01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE01_SAMPLE samples points on the circumference of the unit circle in 2D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2014
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
    Output, double X[2*N], the points.
*/
{
	const dtpitpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int * seed = s_data->a2;
	
    const ityp c[2] =
    {
        0.00,
        0.00
    };
    dim_typ j;
    const register ityp r = 1.00;
    ityp *theta;

    theta = r8vec_uniform_01_new ( n, seed );

    for ( j = 0; j < n; ++j )
    {
        x[0+(j<<1)] = c[0] + r * cos ( M_PI * theta[j] );
        x[1+(j<<1)] = c[0] + r * sin ( M_PI * theta[j] );
    }

    free ( theta );
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_uniform_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_UNIFORM_01 returns a unit pseudorandom r8VEC.
  Discussion:
    This routine implements the recursion
      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )
    The integer arithmetic never requires more than 32 bits,
    including a sign bit.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 August 2004
  Author:
    John Burkardt
  Reference:
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.
    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.
    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.
    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.
  Parameters:
    Input, int N, the number of entries in the vector.
    Input/output, int *SEED, a seed for the random number generator.
    Output, ityp R[N], the vector of pseudorandom values.
*/
{
	const dtpitpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * r = s_data->a1;
	int * seed = s_data->a2;
	
    dim_typ i, k;
    
    if ( *seed == 0 )
        return NULL;

    for ( i = 0; i < n; ++i )
    {
        k = *seed / 127773;
        *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
        if ( *seed < 0 )
            *seed += i4_huge;
        r[i] = ( ityp ) ( *seed ) * 4.656612875E-10;
    }

    return NULL;
}

#endif
