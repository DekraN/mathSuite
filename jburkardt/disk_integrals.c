#ifndef __DISABLEDEEP_DISKINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _disk01_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    DISK01_AREA returns the area of the unit disk in 2D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 January 2014
  Author:
    John Burkardt
  Parameters:
    Output, double DISK01_AREA, the area of the unit disk.
*/
{
	static ityp result = M_PI;
    return &result;
}

#define DISK01MONOMIALINTEGRAL_INVALIDRETURNVALUE MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _disk01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    DISK01_MONOMIAL_INTEGRAL returns monomial integrals in the unit disk in 2D.
  Discussion:
    The integration region is
      X^2 + Y^2 <= 1.
    The monomial is F(X,Y) = X^E(1) * Y^E(2).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 January 2014
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
    Output, double DISK01_MONOMIAL_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * e = data;
	
    ityp arg;
    dim_typ i;
    ityp integral;
    const ityp r = 1.00;
    ityp s;

    if ( e[0] < 0 || e[1] < 0 )
    {
    	result = DISK01MONOMIALINTEGRAL_INVALIDRETURNVALUE;
        return &result;
    }

    if( ( e[0] % 2 ) == 1 || ( e[1] % 2 ) == 1 )
        integral = 0.00;
    else
    {

        #pragma omp parallel for num_threads(2)
        for ( i = 0; i < 2; ++i )
        {
            arg = 0.50 * ( ityp ) ( e[i] + 1 );
            integral *= r8_gamma ( arg );
        }
        arg = 0.50 * ( ityp ) ( e[0] + e[1] + 2 );
        integral /= r8_gamma ( arg );
    }
    /*
    Adjust the surface integral to get the volume integral.
    */
    result = integral * pow ( r, s ) / ( ityp ) ( e[0] + e[1] + 2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _disk01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    DISK01_SAMPLE uniformly samples the unit disk in 2D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 January 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double X[2*N], the points.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    dim_typ i, j;
    ityp norm;
    ityp r;
    ityp *v;
    ityp *x = ( ityp * ) malloc ( (n<<1) * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
    {
        v = r8vec_normal_01_new ( 2, seed );
        /*
        Compute the length of the vector.
        */
        norm = sqrt ( pow ( v[0], 2 ) + pow ( v[1], 2 ) );
        /*
        Normalize the vector.
        */
        #pragma omp parallel for num_threads(2)
        for ( i = 0; i < 2; ++i )
            v[i] /= norm;
        /*
        Now compute a value to map the point ON the circle INTO the circle.
        */
        r = r8_uniform_01 ( seed );

        #pragma omp parallel for num_threads(2)
        for ( i = 0; i < 2; ++i )
            x[i+j*2] = sqrt ( r ) * v[i];
        free ( v );
    }
    return x;
}

#endif
