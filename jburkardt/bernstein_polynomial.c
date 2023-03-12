#ifndef __DISABLEDEEP_BERNSTEINPOLYNOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bernstein_matrix( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNSTEIN_MATRIX returns the Bernstein matrix.
  Discussion:
    The Bernstein matrix of order N is an NxN matrix A which can be used to
    transform a vector of power basis coefficients C representing a polynomial
    P(X) to a corresponding Bernstein basis coefficient vector B:
      B = A * C
    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N
    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
    For N = 5, the matrix has the form:
      1 -4   6  -4  1
      0  4 -12  12 -4
      0  0   6 -12  6
      0  0   0   4 -4
      0  0   0   0  1
    and the numbers in each column represent the coefficients in the power
    series expansion of a Bernstein polynomial, so that
      B(5,4) = - 4 x^4 + 12 x^3 - 12 x^2 + 4 x
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Output, double BERNSTEIN_MATRIX[N*N], the Bernstein matrix.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
	dim_typ i;
	#pragma omp parallel for
	for (dim_typ j = 0; j < n; ++j)
	{
		#pragma omp parallel for
		for (i = 0; i <= j; ++i)
			*(a+i*n+j) = r8_mop(j-i)*r8_choose(n-1-i, j-i)*r8_choose(n-1,i);
		#pragma omp parallel for
		for ( i = j + 1; i < n; ++i)
			*(a+i*n+j) = 0.00;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _bernstein_matrix_inverse( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNSTEIN_MATRIX_INVERSE returns the inverse Bernstein matrix.
  Discussion:
    The inverse Bernstein matrix of order N is an NxN matrix A which can
    be used to transform a vector of Bernstein basis coefficients B
    representing a polynomial P(X) to a corresponding power basis
    coefficient vector C:
      C = A * B
    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N
    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
    For N = 5, the matrix has the form:
      1   1    1    1   1
      0  1/4  1/2  3/4  1
      0   0   1/6  1/2  1
      0   0    0   1/4  1
      0   0    0    0   1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Output, double BERNSTEIN_MATRIX_INVERSE[N*N], the inverse Bernstein matrix.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
	dim_typ i;
	#pragma omp parallel for
	for (dim_typ j = 0; j < n; ++j)
	{
		#pragma omp parallel for
		for (i = 0; i <= j; ++i)
			*(a+i*n+j) = r8_choose(j,i)/r8_choose(n-1,i);
		#pragma omp parallel for
		for (i = j + 1; i < n; ++i)
			*(a+i*n+j) = 0.00;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bernstein_poly_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNSTEIN_POLY_01 evaluates the Bernstein polynomials based in [0,1].
  Discussion:
    The Bernstein polynomials are assumed to be based on [0,1].
    The formula is:
      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
  First values:
    B(0,0)(X) = 1
    B(1,0)(X) =      1-X
    B(1,1)(X) =                X
    B(2,0)(X) =  (1-X)^2
    B(2,1)(X) = 2 * (1-X)    * X
    B(2,2)(X) =                X^2
    B(3,0)(X) =  (1-X)^3
    B(3,1)(X) = 3 * (1-X)^2 * X
    B(3,2)(X) = 3 * (1-X)   * X^2
    B(3,3)(X) =               X^3
    B(4,0)(X) =  (1-X)^4
    B(4,1)(X) = 4 * (1-X)^3 * X
    B(4,2)(X) = 6 * (1-X)^2 * X^2
    B(4,3)(X) = 4 * (1-X)   * X^3
    B(4,4)(X) =               X^4
  Special values:
    B(N,I)(X) has a unique maximum value at X = I/N.
    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
    B(N,I)(1/2) = C(N,K) / 2^N
    For a fixed X and N, the polynomials add up to 1:
      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the degree of the Bernstein polynomials
    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
    each of degree N, which form a basis for polynomials on [0,1].
    Input, double X, the evaluation point.
    Output, double BERNSTEIN_POLY[N+1], the values of the N+1
    Bernstein polynomials at X.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * bern = s_data->a1;
	const register ityp x = s_data->a2;
	
	dim_typ j;
	if(!n)
		bern[0] = 1.00;
	else if (0 < n)
	{
		bern[0] = 1.00 - x;
		bern[1] = x;

		for (dim_typ i=2; i <= n; ++i)
		{
			bern[i] = x * bern[i-1];
			for (j = i - 1; 1 <= j; --j)
				bern[j] = x*bern[j-1]+(1.00-x)*bern[j];
			bern[0] = ( 1.0 - x ) * bern[0];
		}
	}
	return NULL;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bernstein_poly_ab( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNSTEIN_POLY_AB evaluates the Bernstein polynomials based in [A,B].
  Discussion:
    The formula is:
      BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
  First values:
    B(0,0)(X) =   1
    B(1,0)(X) = (      B-X                ) / (B-A)
    B(1,1)(X) = (                 X-A     ) / (B-A)
    B(2,0)(X) = (  (B-X)^2             ) / (B-A)^2
    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
    B(2,2)(X) = (             (X-A)^2  ) / (B-A)^2
    B(3,0)(X) = (  (B-X)^3             ) / (B-A)^3
    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
    B(3,3)(X) = (             (X-A)^3  ) / (B-A)^3
    B(4,0)(X) = (  (B-X)^4             ) / (B-A)^4
    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
    B(4,4)(X) = (             (X-A)^4  ) / (B-A)^4
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the degree of the Bernstein polynomials
    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
    each of degree N, which form a basis for polynomials on [A,B].
    Input, double A, B, the endpoints of the interval on which the
    polynomials are to be based.  A and B should not be equal.
    Input, double X, the point at which the polynomials
    are to be evaluated.
    Output, double BPAB[N+1], the values of the N+1
    Bernstein polynomials at X.
*/
{
	static bool result = 2;
	
	const dt3itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	const register ityp x = s_data->a3;
	ityp * bern = s_data->a4;
	
	if(b == a)
	{
		result = BERNSTEIN_DOMAIN_ERROR;
		return &result;
	}

	if(!n)
		bern[0] = 1.00;
	else if (0 < n)
	{
		dim_typ j;
		bern[0] = (b-x)/(b-a);
		bern[1] = (x-a)/(b-a);

		for (dim_typ i=2; i <= n; ++i)
		{
			bern[i] = (x-a ) * bern[i-1] / (b-a);
			for ( j = i - 1; 1 <= j; --j)
				bern[j] = ((b-x)*bern[j]+(x-a)*bern[j-1])/(b-a);
			bern[0] = (b-x)*bern[0]/(b-a);
		}
	}

	result = BERNSTEIN_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bernstein_poly_01_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNSTEIN_POLY_01_VALUES returns some values of the Bernstein polynomials.
  Discussion:
    The Bernstein polynomials are assumed to be based on [0,1].
    The formula for the Bernstein polynomials is
      B(N,I)(X) = [N!/(I(N-I)!)] * (1-X)^(N-I) * X^I
    In Mathematica, the function can be evaluated by:
      Binomial[n,i] * (1-x)^(n-i) * x^i
  First values:
    B(0,0)(X) = 1
    B(1,0)(X) =      1-X
    B(1,1)(X) =                X
    B(2,0)(X) =  (1-X)^2
    B(2,1)(X) = 2 * (1-X)    * X
    B(2,2)(X) =                X^2
    B(3,0)(X) =  (1-X)^3
    B(3,1)(X) = 3 * (1-X)^2  * X
    B(3,2)(X) = 3 * (1-X)    * X^2
    B(3,3)(X) =                X^3
    B(4,0)(X) =  (1-X)^4
    B(4,1)(X) = 4 * (1-X)^3  * X
    B(4,2)(X) = 6 * (1-X)^2  * X^2
    B(4,3)(X) = 4 * (1-X)    * X^3
    B(4,4)(X) =                X^4
  Special values:
    B(N,I)(X) has a unique maximum value at X = I/N.
    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
    B(N,I)(1/2) = C(N,K) / 2^N
    For a fixed X and N, the polynomials add up to 1:
      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2012
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the degree of the polynomial.
    Output, int *K, the index of the polynomial.
    Output, double *X, the argument of the polynomial.
    Output, double *B, the value of the polynomial B(N,K)(X).
*/
{
	const _3pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	dim_typ * k = s_data->a2;
	ityp * x = s_data->a3;
	ityp * b = s_data->a4;
	
    # define N_MAX 15

    static ityp b_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.7500000000000000E+00,
        0.2500000000000000E+00,
        0.5625000000000000E+00,
        0.3750000000000000E+00,
        0.6250000000000000E-01,
        0.4218750000000000E+00,
        0.4218750000000000E+00,
        0.1406250000000000E+00,
        0.1562500000000000E-01,
        0.3164062500000000E+00,
        0.4218750000000000E+00,
        0.2109375000000000E+00,
        0.4687500000000000E-01,
        0.3906250000000000E-02
    };

    static dim_typ k_vec[N_MAX] =
    {
        0,
        0, 1,
        0, 1, 2,
        0, 1, 2, 3,
        0, 1, 2, 3, 4
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,
        1, 1,
        2, 2, 2,
        3, 3, 3, 3,
        4, 4, 4, 4, 4
    };

    static ityp x_vec[N_MAX] =
    {
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = *k = 0;
        *x = *b = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *k = k_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *b = b_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_choose ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
  Discussion:
    The value is calculated in such a way as to avoid overflow and
    roundoff.  The calculation is done in R8 arithmetic.
    The formula used is:
      C(N,K) = N! / ( K! * (N-K)! )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2008
  Author:
    John Burkardt
  Reference:
    ML Wolfson, HV Wright,
    Algorithm 160:
    Combinatorial of M Things Taken N at a Time,
    Communications of the ACM,
    Volume 6, Number 4, April 1963, page 161.
  Parameters:
    Input, int N, K, the values of N and K.
    Output, double R8_CHOOSE, the number of combinations of N
    things taken K at a time.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ k = a_data[1];
	
	dim_typ mx;
	ityp value;
	const register dim_typ mn = MIN (k, n - k);

	if(mn <= 0)
		value = 0.00 + (!mn);
	else
	{
		mx = MAX ( k, n - k );
		value = (ityp) (mx + 1);

		for (dim_typ i = 2; i <= mn; ++i)
			value *= (ityp)(mx+i)/(ityp)i;
	}

	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_sign( void * data)
/******************************************************************************/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	result = -1 + ((x>=0.00)<<1);
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_linspace_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
  Discussion:
    An R8VEC is a vector of R8's.
    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 March 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, double A, B, the first and last entries.
    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
*/
{
	const _2itdt * const s_data = data;
	
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    dim_typ i;
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
        x[0] = ( a + b ) / 2.00;
    else
        for ( i = 0; i < n; ++i )
            x[i] = ( ( ityp ) ( n - 1 - i ) * a + ( ityp ) (         i ) * b ) / ( ityp ) ( n - 1     );
   return x;
}

#endif
