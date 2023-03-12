#ifndef __DISABLEDEEP_SPHEREINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere01_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_MONOMIAL_INTEGRAL: integrals on the surface of the unit sphere in 3D.
  Discussion:
    The integration region is 
      X^2 + Y^2 + Z^2 = 1.
    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    05 January 2014
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Academic Press, 1984, page 263.
  Parameters:
    Input, int E[3], the exponents of X, Y and Z in the 
    monomial.  Each exponent must be nonnegative.
    Output, double SPHERE01_MONOMIAL_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	int * e = data;
	
	ityp arg;
	dim_typ i;
	ityp integral;
	
	for ( i = 0; i < 3; ++i )
		if ( e[i] < 0 )
		{
			result = MAX_VAL;
			return &result;
		}
	
	for ( i = 0; i < 3; ++i )
		if ( ( e[i] % 2 ) == 1 )
		{
			result = 0.00;
			return &result;
		}
	
	integral = 2.00;

	#pragma omp parallel for num_threads(3)
	for ( i = 0; i < 3; ++i )
	{
		arg = 0.50 * ( ityp ) ( e[i] + 1 );
		integral *= r8_gamma ( arg );
	}
	
	result = integral / r8_gamma ( 0.50 * ( ityp ) ( e[0] + e[1] + e[2] + 3 ) );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_SAMPLE uniformly sample from the surface of the unit sphere in 3D.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    25 September 2010
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
    Output, double X[3*N], the points.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
	dim_typ i, j;
	ityp norm;
	ityp *x = r8mat_normal_01_new ( 3, n, seed );
	
	for ( j = 0; j < n; ++j )
	{
		/*
		Compute the length of the vector.
		*/
		norm = sqrt ( pow ( x[0+j*3], 2 ) + pow ( x[1+j*3], 2 )+ pow ( x[2+j*3], 2 ) );
		/*
		Normalize the vector.
		*/
		#pragma omp parallel for num_threads(3)
		for ( i = 0; i < 3; ++i)
			x[i+j*3] /= norm;
	}
	return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere01_monomial_quadrature ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_MONOMIAL_QUADRATURE applies quadrature to a monomial in a sphere.
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
    Input, double XYZ[DIM_NUM*POINT_NUM], the quadrature points.
    Input, double W[POINT_NUM], the quadrature weights.
    Output, double SPHERE01_MONOMIAL_QUADRATURE, the quadrature error.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitdtpipit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	ityp * expon = s_data->a1;
	const register dim_typ point_num = s_data->a2;
	int * xyz = s_data->a3;
	ityp * w = s_data->a4;

    ityp exact;
    ityp quad;
    ityp quad_error;
    ityp *value;
    ityp volume;
    /*
    Get the exact value of the integral.
    */
    exact = sphere01_monomial_integral ( (int*)expon );
    /*
    Evaluate the monomial at the quadrature points.
    */
    value = monomial_value ( 3, point_num, expon, xyz );
    /*
    Compute the weighted sum.
    */
    quad = r8vec_dot_product ( point_num, w, value );
    /*
    Error:
    */
    quad_error = fabs ( quad - exact );
    free ( value );
    
    result = quad_error;
    return &result;
}

#endif
