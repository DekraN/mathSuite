#ifndef __DISABLEDEEP_JBURKARDT

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8poly_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8POLY_VALUE evaluates a double precision polynomial.
  Discussion:
    For sanity's sake, the value of N indicates the NUMBER of
    coefficients, or more precisely, the ORDER of the polynomial,
    rather than the DEGREE of the polynomial.  The two quantities
    differ by 1, but cause a great deal of confusion.
    Given N and A, the form of the polynomial is:
      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 March 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the polynomial.
    Input, double A[N], the coefficients of the polynomial.
    A[0] is the constant term.
    Input, double X, the point at which the polynomial is to be evaluated.
    Output, double R8POLY_VALUE, the value of the polynomial at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp x = s_data->a2;
	
	ityp value = 0.00;
	for (dim_typ i = n - 1; 0 <= i; --i )
		value = value*x + a[i];
		
	result = value; 
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _p_polynomial_prime2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    P_POLYNOMIAL_PRIME2: second derivative of Legendre polynomials P(n,x).
  Discussion:
    P(0,X) = 1
    P(1,X) = X
    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
    P'(0,X) = 0
    P'(1,X) = 1
    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
    P"(0,X) = 0
    P"(1,X) = 0
    P"(N,X) = ( (2*N-1)*(2*P'(N-1,X)+X*P"(N-1,X)-(N-1)*P"(N-2,X) ) / N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.
    Input, double X[M], the evaluation points.
    Output, double P_POLYNOMIAL_PRIME2[M*(N+1)], the second derivative
    of the Legendre polynomials of order 0 through N at the points.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	ityp v[n+1];
	ityp vp[n+1];
	ityp vpp[n+1];

	if(n < 0)
	{
		result = LEGENDRE_INVALIDRETURNVALUE;
		return &result;
	}

	if(!n)
	{
		result = 0.00;
		return &result;
	}

	vpp[0] = 0.00;
	v[0] = 1.00;
	vp[0] = 0.00;

	v[1] = x;
	vp[1] = 1.00;
	vpp[1] = 0.00;

	for (dim_typ i = 2; i <= n; ++i )
	{
		v[i] = ((ityp)((i<<1)-1)*x*v[i-1]-(ityp)(i-1)*v[i-2])/(ityp)(i);
		vp[i] = ((ityp)((i<<1)-1)*(v[i-1]+x*vp[i-1])-(ityp)(i-1)*vp[i-2])/(ityp)(i);
		vpp[i] = ((ityp)((i<<1)-1)*(2.00*vp[i-1]+x*vpp[i-1])-(ityp)(i-1)*vpp[i-2])/(ityp)(i);
	}

	result = vpp[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _p_polynomial_value( void * data)
/******************************************************************************/
/*
  Purpose:
    P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
  Discussion:
    P(n,1) = 1.
    P(n,-1) = (-1)^N.
    | P(n,x) | <= 1 in [-1,1].
    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
    quadrature of the integral of a function F(X) with weight function 1
    over the interval [-1,1].
    The Legendre polynomials are orthogonal under the inner product defined
    as integration from -1 to 1:
      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
        = 0 if I =/= J
        = 2 / ( 2*I+1 ) if I = J.
    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
    A function F(X) defined on [-1,1] may be approximated by the series
      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
    where
      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
    The formula is:
      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
  Differential equation:
 (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
  First terms:
    P( 0,x) =      1
    P( 1,x) =      1 X
    P( 2,x) = (    3 X^2 -       1)/2
    P( 3,x) = (    5 X^3 -     3 X)/2
    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
  Recursion:
    P(0,x) = 1
    P(1,x) = x
    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
    P'(0,x) = 0
    P'(1,x) = 1
    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.
    Input, double X[M], the evaluation points.
    Output, double P_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
    polynomials of order 0 through N.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	ityp v[n+1];

	if(n < 0)
	{
		result = LEGENDRE_INVALIDRETURNVALUE;
		return &result;
	}

	if(n < 1)
	{
		result = 1.00;
		return &result;
	}

	v[0] = 1.00;
	v[1] = x;

	for (dim_typ i=2; i <= n; ++i )
		v[i] = ((ityp)((i<<1)-1)*x*v[i-1]-(ityp)(i-1)*v[i-2])/(ityp)(i);

	result = v[n]; 
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pm_polynomial_value( void * data)
/******************************************************************************/
/*
  Purpose:
    PM_POLYNOMIAL_VALUE evaluates the Legendre polynomials Pm(n,m,x).
  Differential equation:
 (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
  First terms:
    M = 0 ( = Legendre polynomials of first kind P(N,X) )
    Pm(0,0,x) =    1
    Pm(1,0,x) =    1 X
    Pm(2,0,x) = (  3 X^2 -   1)/2
    Pm(3,0,x) = (  5 X^3 -   3 X)/2
    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
    M = 1
    Pm(0,1,x) =   0
    Pm(1,1,x) =   1 * SQRT(1-X^2)
    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
    M = 2
    Pm(0,2,x) =   0
    Pm(1,2,x) =   0
    Pm(2,2,x) =   3 * (1-X^2)
    Pm(3,2,x) =  15 * (1-X^2) * X
    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)
    M = 3
    Pm(0,3,x) =   0
    Pm(1,3,x) =   0
    Pm(2,3,x) =   0
    Pm(3,3,x) =  15 * (1-X^2)^1.5
    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X
    M = 4
    Pm(0,4,x) =   0
    Pm(1,4,x) =   0
    Pm(2,4,x) =   0
    Pm(3,4,x) =   0
    Pm(4,4,x) = 105 * (1-X^2)^2
  Recursion:
    if N < M:
      Pm(N,M,x) = 0
    if N = M:
      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
      all the odd integers less than or equal to N.
    if N = M+1:
      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
    if M+1 < N:
      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input, int MM, the number of evaluation points.
    Input, int N, the maximum first index of the Legendre
    function, which must be at least 0.
    Input, int M, the second index of the Legendre function,
    which must be at least 0, and no greater than N.
    Input, double X[MM], the point at which the function is to be
    evaluated.
    Output, double PM_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register ityp x = s_data->a2;
	
	dim_typ i;
	ityp fact;
	ityp v[n+1];

	#pragma omp parallel for
	for (i=0; i < n+1; ++i )
		v[i] = 0.00;

	/*
	J = M is the first nonzero function.
	*/
	if ( m <= n )
	{
		v[m] = 1.00;
		fact = 1.00;
		#pragma omp parallel for
		for ( i= 0; i < m; ++i )
		{
			v[m] *= -(fact*sqrt (1.00-(x*x)));
			fact += 2;
		}
	}

	/*
	J = M + 1 is the second nonzero function.
	*/
	if (m+1 <= n)
		v[m+1] = x*(ityp)((m<<1)+1)*v[m];

	/*
	Now we use a three term recurrence.
	*/
	for (i=m+2; i<=n; ++i)
		v[i] = ((ityp)((i<<1)-1)*x*v[i-1]+(ityp)(-i-m+1)*v[i-2])/(ityp)(i-m);

	result = v[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pmn_polynomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    PMN_POLYNOMIAL_VALUE: normalized Legendre polynomial Pmn(n,m,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input, int MM, the number of evaluation points.
    Input, int N, the maximum first index of the Legendre
    function, which must be at least 0.
    Input, int M, the second index of the Legendre function,
    which must be at least 0, and no greater than N.
    Input, double X[MM], the evaluation points.
    Output, double PMN_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register ityp x = s_data->a2;
	
	const register ityp v = pm_polynomial_value(n, m, x);
	
	result = (v != LEGENDRE_INVALIDRETURNVALUE ? v*sqrt(((ityp)((n<<1)+1)*r8_factorial(n-m))/(2.00*r8_factorial(n+m))) : v);
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pmns_polynomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    PMNS_POLYNOMIAL_VALUE: sphere-normalized Legendre polynomial Pmn(n,m,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input, int MM, the number of evaluation points.
    Input, int N, the maximum first index of the Legendre
    function, which must be at least 0.
    Input, int M, the second index of the Legendre function,
    which must be at least 0, and no greater than N.
    Input, double X[MM], the evaluation points.
    Output, double PMNS_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register ityp x = s_data->a2;
	
  const register ityp v = pm_polynomial_value(n, m, x );
  
  result = (v != LEGENDRE_INVALIDRETURNVALUE ? v*sqrt(((ityp)((n<<1)+1)*r8_factorial(n-m))/(4.0*M_PI*r8_factorial(n + m))) : v);
  return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pn_polynomial_value( void * data)
/******************************************************************************/
/*
  Purpose:
    PN_POLYNOMIAL_POLYNOMIAL evaluates the normalized Legendre polynomials Pn(n,x).
  Discussion:
    The normalized Legendre polynomials are orthonormal under the inner product
    defined as integration from -1 to 1:
      Integral ( -1 <= x <= +1 ) Pn(i,x) * Pn(j,x) dx = delta(i,j)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.
    Input, double X[M], the evaluation points.
    Output, double PN_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
    polynomials of order 0 through N.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	const register ityp v = p_polynomial_value(n, x);
	
	result = (v != LEGENDRE_INVALIDRETURNVALUE ? (v/sqrt (2/(ityp)((n<<1)+1))) : LEGENDRE_INVALIDRETURNVALUE);
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pn_polynomial_coefficients ( void * data)
/******************************************************************************/
/*
  Purpose:
    PN_POLYNOMIAL_COEFFICIENTS: coefficients of normalized Legendre Pn(n,x).
  Discussion:
    Pn(n,x) = P(n,x) * sqrt ( (2n+1)/2 )
          1       x       x^2     x^3     x^4      x^5    x^6     x^7
    0   0.707
    1   0.000   1.224
    2  -0.790   0.000   2.371
    3   0.000  -2.806   0.000   4.677
    4   0.795   0.000  -7.954   0.000   9.280
    5   0.000   4.397   0.000 -20.520   0.000   18.468
    6  -0.796   0.000  16.731   0.000 -50.193    0.000  36.808
    7   0.000  -5.990   0.000  53.916   0.000 -118.616   0.000  73.429
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 October 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.
    Output, double PN_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients of
    the normalized Legendre polynomials of degree 0 through N.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * c = s_data->a1;
	
	if (n<0)
	{
		result = LEGENDRE_DOMAIN_ERROR;
		return &result;
	}

	dim_typ i, j;

	#pragma omp parallel for
	for (i=0; i<=n; ++i )
		#pragma omp parallel for
		for (j=0; j<=n; ++j )
			*(c+i*(n+1)+j) = 0.00;

	*c = 1.00;

	if (!n)
	{
		result = LEGENDRE_SUCCESS;
		return &result;
	}

	*(c+(n+1)+1) = 1.00;

	for (i=2; i <= n; ++i)
	{
		for ( j = 0; j <= i-2; ++j)
		*(c+i*(n+1)+j) = (ityp)(-i+1)* *(c+(i-2)*(n+1)+j)/(ityp) i;

		for ( j = 1; j <= i; ++j)
			*(c+i*(n+1)+j) = *(c+i*(n+1)+j)+(ityp)(i+i-1)* *(c+(i-1)*(n+1)+j-1)/(ityp)i;
	}

	register ityp t;

	/*
	Normalize them.
	*/
	#pragma omp parallel for
	for (i=0; i <= n; ++i)
	{
		t = sqrt((ityp)((i<<1)+1)/2.00);
		#pragma omp parallel for
		for (j = 0; j <= i; ++j)
			*(c+i*(n+1)+j) *= t;
	}

	result = LEGENDRE_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_normal_01_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
  Discussion:
    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.
    This routine can generate a vector of values on one call.  It
    has the feature that it should provide the same results
    in the same order no matter how we break up the task.
    Before calling this routine, the user may call RANDOM_SEED
    in order to set the seed of the random number generator.
    The Box-Muller method is used, which is efficient, but
    generates an even number of values each time.  On any call
    to this routine, an even number of new values are generated.
    Depending on the situation, one value may be left over.
    In that case, it is saved for the next call.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of values desired.  If N is negative,
    then the code will flush its internal memory; in particular,
    if there is a saved value to be used on the next call, it is
    instead discarded.  This is useful if the user has reset the
    random number seed, for instance.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
  Local parameters:
    Local, int MADE, records the number of values that have
    been computed.  On input with negative N, this value overwrites
    the return value of N, so the user can get an accounting of
    how much work has been done.
    Local, double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.
    Local, int SAVED, is 0 or 1 depending on whether there is a
    single saved value left over from the previous call.
    Local, int X_LO, X_HI, records the range of entries of
    X that we need to compute.  This starts off as 1:N, but is adjusted
    if we have a saved value that can be immediately stored in X(1),
    and so on.
    Local, double Y, the value saved from the previous call, if
    SAVED is 1.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;

	dim_typ i, m;
	static dim_typ  made = 0;
	double *r;
	static dim_typ saved = 0;
	ityp * x = malloc(sizeof(ityp)*n);
	dim_typ x_hi, x_lo;
	static ityp y = 0.00;

	/*
	I'd like to allow the user to reset the internal data.
	But this won't work properly if we have a saved value Y.
	I'm making a crock option that allows the user to signal
	explicitly that any internal memory should be flushed,
	by passing in a negative value for N.
	*/
	if ( n < 0 )
	{
		made = 0;
		saved = 0;
		y = 0.00;
		return false;
	}
	else if (!n)
		return false;
	/*
	Record the range of X we need to fill in.
	*/
	x_lo = 1;
	x_hi = n;
	/*
	Use up the old value, if we have it.
	*/
	if ( saved == 1 )
	{
		x[0] = y;
		saved = 0;
		x_lo = 2;
	}
	/*
	Maybe we don't need any more values.
	*/
	if (!(x_hi - x_lo + 1));

	/*
	If we need just one new value, do that here to avoid null arrays.
	*/
	else if ( x_hi - x_lo + 1 == 1 )
	{
		r = r8vec_uniform_01_new ( 2, seed );

		x[x_hi-1] = sqrt ( - 2.00 * log ( r[0] ) ) * cos ( M_2TPI * r[1] );
		y =         sqrt ( - 2.00 * log ( r[0] ) ) * sin ( M_2TPI * r[1] );

		saved = 1;
		made += 2;
		free ( r );
	}
	/*
	If we require an even number of values, that's easy.
	*/
	else if ( !(( x_hi - x_lo + 1 ) % 2))
	{
		m = ( x_hi - x_lo + 1 ) / 2;
		r = r8vec_uniform_01_new ( m<<1, seed );

		for ( i = 0; i <= 2*m-2; i = i + 2 )
		{
			x[x_lo+i-1] = sqrt ( - 2.00 * log ( r[i] ) ) * cos ( M_2TPI * r[i+1] );
			x[x_lo+i  ] = sqrt ( - 2.00 * log ( r[i] ) ) * sin ( M_2TPI * r[i+1] );
		}
		made += x_hi - x_lo + 1;
		free ( r );
	}
	/*
	If we require an odd number of values, we generate an even number,
	and handle the last pair specially, storing one in X(N), and
	saving the other for later.
	*/
	else
	{
		x_hi = x_hi - 1;
		m = ( x_hi - x_lo + 1 ) / 2 + 1;
		r = r8vec_uniform_01_new ( 2*m, seed );

		for ( i = 0; i <= (m<<1)-4; i += 2 )
		{
			x[x_lo+i-1] = sqrt ( - 2.00 * log ( r[i] ) ) * cos ( M_2TPI * r[i+1] );
			x[x_lo+i  ] = sqrt ( - 2.00 * log ( r[i] ) ) * sin ( M_2TPI * r[i+1] );
		}

		i = (m<<1) - 2;

		x[x_lo+i-1] = sqrt ( - 2.00 * log ( r[i] ) ) * cos ( M_2TPI * r[i+1] );
		y           = sqrt ( - 2.00 * log ( r[i] ) ) * sin ( M_2TPI * r[i+1] );

		saved = 1;

		made += x_hi - x_lo + 2;
		free ( r );
	}
	return x;
}

#define I4VECUNIFORMABNEW_ZEROSEEDERROR NULL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_uniform_ab_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_UNIFORM_AB_NEW returns a scaled pseudorandom I4VEC.
  Discussion:
    The pseudorandom numbers should be uniformly distributed
    between A and B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 January 2014
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
    Input, integer N, the dimension of the vector.
    Input, int A, B, the limits of the interval.
    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.
    Output, int I4VEC_UNIFORM_AB_NEW[N], a vector of random values
    between A and B.
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	register dim_typ a = s_data->a1;
	register dim_typ b = s_data->a2;
	int * seed = s_data->a3; 
	
    dim_typ c;
    dim_typ i;
    dim_typ k;
    ityp r;
    dim_typ value;
    int *x;

    if ( !(*seed) )
        return I4VECUNIFORMABNEW_ZEROSEEDERROR;
    /*
    Guaranteee A <= B.
    */
    if ( b < a )
    {
        c = a;
        a = b;
        b = c;
    }

    x = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
    {
        k = *seed / 127773;
        *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

        if ( *seed < 0 )
        {

            r = ( ityp ) ( *seed ) * 4.656612875E-10;
            /*
            Scale R to lie between A-0.5 and B+0.5.
            */
            r = ( 1.0 - r ) * ( ( ityp ) a - 0.5 ) +         r   * ( ( ityp ) b + 0.5 );
            /*
            Use rounding to convert R to an integer between A and B.
            */
            value = round ( r );
            /*
            Guarantee A <= VALUE <= B.
            */
            value = value<a ?a:b;
        }

        x[i] = value;
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gauss ( void * data)
/******************************************************************************/
/*
  Purpose:
    GAUSS computes a Gauss quadrature rule.
  Discussion:
    Given a weight function W encoded by the first N recurrence coefficients
    ALPHA and BETA for the associated orthogonal polynomials, the call
      call gauss ( n, alpha, beta, x, w )
    generates the nodes and weights of the N-point Gauss quadrature rule
    for the weight function W.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2013
  Author:
    Original MATLAB version by Walter Gautschi.
    C version by John Burkardt.
  Reference:
    Walter Gautschi,
    Orthogonal Polynomials: Computation and Approximation,
    Oxford, 2004,
    ISBN: 0-19-850672-4,
    LC: QA404.5 G3555.
  Parameters:
    Input, int N, the order of the desired quadrature rule.
    Input, double ALPHA[N], BETA[N], the alpha and beta recurrence
    coefficients for the othogonal polynomials associated with the
    weight function.
    Output, double X[N], W[N], the nodes and  weights of the desired
    quadrature rule.  The nodes are listed in increasing order.
*/
{
	const dt4pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * beta = s_data->a2;
	ityp * x = s_data->a3;
	ityp * w = s_data->a4;
	
    dim_typ i, j;
    dim_typ it_max;
    int it_num;
    int rot_num;
    ityp a[n*n];
    ityp v[n*n];
    /*
    Define the tridiagonal Jacobi matrix.
    */

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            if ( i == j )
                a[i+j*n] = alpha[i];
            else if ( i == j - 1 )
                a[i+j*n] = sqrt ( beta[j] );
            else if ( i - 1 == j )
                a[i+j*n] = sqrt ( beta[i] );
            else
                a[i+j*n] = 0.00;
    /*
    Get the eigenvectors and eigenvalues.
    */
    it_max = 100;
    jacobi_eigenvalue ( n, a, it_max, v, x, &it_num, &rot_num );

    for ( j = 0; j < n; ++j )
        w[j] = beta[0] * v[0+j*n] * v[0+j*n];
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r_jacobi ( void * data)
/******************************************************************************/
/*
  Purpose:
    R_JACOBI computes recurrence coefficients for monic Jacobi polynomials.
  Discussion:
    This function generates the first N recurrence coefficients for monic
    Jacobi polynomials with parameters A and B.
    These polynomials are orthogonal on [-1,1] relative to the weight
      w(x) = (1.0-x)^A * (1.0+x)^B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 July 2013
  Author:
    Original MATLAB version by Dirk Laurie, Walter Gautschi.
    C version by John Burkardt.
  Reference:
    Walter Gautschi,
    Orthogonal Polynomials: Computation and Approximation,
    Oxford, 2004,
    ISBN: 0-19-850672-4,
    LC: QA404.5 G3555.
  Parameters:
    Input, int N, the number of coefficients desired.
    Input, double A, B, the parameters for the Jacobi polynomial.
    -1.0 < A, -1.0 < B.
    Output, double ALPHA[N], BETA[N], the first N recurrence
    coefficients.
*/
{
	const dt2it2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	ityp * alpha = s_data->a3;
	ityp * beta = s_data->a4;
	
    int i;
    double i_r8;
    ityp mu;
    double nab;
    double nu;

    if ( a <= -1.00 || b <= -1.00)
        return NULL;

    nu = ( b - a ) / ( a + b + 2.00 );

    mu = pow ( 2.00, a + b + 1.00 )* r8_gamma ( a + 1.00 ) * r8_gamma ( b + 1.00 ) / r8_gamma ( a + b + 2.00 );

    alpha[0] = nu;
    beta[0] = mu;

    if ( n == 1 )
        return NULL;

    for ( i = 1; i < n; ++i)
    {
        i_r8 = ( ityp ) ( i + 1 );
        alpha[i] = ( b - a ) * ( b + a ) / ( 2.00 * ( i_r8 - 1.00 ) + a + b ) / ( 2.00 * i_r8 + a + b );
    }

    beta[1] = 4.00 * ( a + 1.00 ) * ( b + 1.00 ) / ( a + b + 2.00 ) / ( a + b + 2.00 )/ ( a + b + 3.00 );

    for ( i = 2; i < n; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        nab = 2.00 * ( i_r8 - 1.00 ) + a +  b;
        beta[i] = 4.00 * ( i_r8 - 1.00 + a ) * ( i_r8 - 1.00 + b )
        * ( i_r8 - 1.00 ) * ( i_r8 - 1.00 + a + b ) / nab / nab/ ( nab + 1.00 ) / ( nab - 1.00 );
    }

    return NULL;
}

#define R8NORMAL01_INVALIDRETURNVALUE MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _diaphony_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIAPHONY_COMPUTE evaluates the diaphony of a N-dimensional point set.
  Discussion:
    The diaphony is analogous to, and related to, the discrepancy,
    and is a measure of how well spread a set of point is.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 January 2012
  Author:
    John Burkardt
  Reference:
    Peter Heelekalek, Harald Niederreiter,
    The Weighted Spectral Test: Diaphony,
    ACM Transactions on Modeling and Computer Simulation,
    Volume 8, Number 1, January 1998, pages 43-60.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int POINT_NUM, the number of points.
    Input, double X[DIM_NUM*POINT_NUM], the point set, which is
    presumed to lie in the DIM_NUM dimensional unit hypercube.
    Output, double DIAPHONY_COMPUTE, the value of the diaphony.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ point_num = s_data->a1;
	ityp * x = s_data->a2;
	
    ityp bot;
    ityp d;
    dim_typ i, j, k;
    ityp prod;
    ityp z;

    d = 0.00;

    for ( i = 0; i < point_num; ++i )
    {
        for ( j = 0; j < point_num; ++j )
        {
            prod = 1.00;
            for ( k = 0; k < dim_num; ++k )
            {
                z = r8_modp ( x[k+i*dim_num] - x[k+j*dim_num], 1.00 );
                prod = prod * ( 1.00 + M_2TPI * M_PI * ( z * z - z + 1.00 / 6.00 ) );
            }
            d += prod - 1.00;
        }
    }
    
    result = sqrt ( d/((ityp) powi ( point_num, 2 )* ( pow ( 1.00 + M_PI * M_PI / 3.00, dim_num ) - 1.00 )));
    return &result;
}

#define R8MODP_INVALIDRETURNVALUE MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_modp ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_MODP returns the nonnegative remainder of R8 division.
  Formula:
    If
      REM = R8_MODP ( X, Y )
      RMULT = ( X - REM ) / Y
    then
      X = Y * RMULT + REM
    where REM is always nonnegative.
  Discussion:
    The MOD function computes a result with the same sign as the
    quantity being divided.  Thus, suppose you had an angle A,
    and you wanted to ensure that it was between 0 and 360.
    Then mod(A,360.0) would do, if A was positive, but if A
    was negative, your result would be between -360 and 0.
    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
  Example:
        I         J     MOD  R8_MODP   R8_MODP Factorization
      107        50       7       7    107 =  2 *  50 + 7
      107       -50       7       7    107 = -2 * -50 + 7
     -107        50      -7      43   -107 = -3 *  50 + 43
     -107       -50      -7      43   -107 =  3 * -50 + 43
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2006
  Author:
    John Burkardt
  Parameters:
    Input, double X, the number to be divided.
    Input, double Y, the number that divides X.
    Output, double R8_MODP, the nonnegative remainder when X is divided by Y.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp y = a_data[1];
	
    ityp value;

    if ( y == 0.00 )
    {
    	result = R8MODP_INVALIDRETURNVALUE;
        return &result;
    }

    value = x - ( ( ityp ) ( ( int ) ( x / y ) ) ) * y;

    if ( value < 0.00 )
        value +=  fabs ( y );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _inverse_error ( void * data)
/******************************************************************************/
/*
  Purpose:
    INVERSE_ERROR determines the error in an inverse matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the matrix.
    Input, double B[N*N], the inverse.
    Output, double ERROR_FROBENIUS, the Frobenius norm
    of (A*B-I) + (B*A-I).
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    ityp *c;
    dim_typ j;
    ityp value;

    c = r8mat_mm_new ( n, n, n, a, b );

    for ( j = 0; j < n; ++j )
        c[j+j*n] -= 1.0;

    value = r8mat_norm_fro ( n, n, c );
    free ( c );
    c = r8mat_mm_new ( n, n, n, b, a );

    for ( j = 0; j < n; ++j)
    c[j+j*n] -= 1.00;

    value += r8mat_norm_fro ( n, n, c );

    free ( c );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8po_sl2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8PO_SL solves a linear system that has been factored by R8PO_FA.
  Discussion:
    The R8PO storage format is appropriate for a symmetric positive definite
    matrix and its inverse. (The Cholesky factor of a R8PO matrix is an
    upper triangular matrix, so it will be in R8GE storage format.)
    Only the diagonal and upper triangle of the square array are used.
    This same storage format is used when the matrix is factored by
    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
    is set to zero.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 February 2013
  Author:
    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
    C version by John Burkardt.
  Reference:
    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A_LU[N*N], the Cholesky factor from R8PO_FA.
    Input, double B[N], the right hand side.
    Output, double R8PO_SL[N], the solution vector.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a_lu = s_data->a1;
	ityp * b = s_data->a2;
	
    dim_typ i, k;
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );


    for ( k = 0; k < n; ++k )
        x[k] = b[k];
    /*
    Solve R' * y = b.
    */
    for ( k = 0; k < n; ++k )
    {
        for ( i = 0; i < k; ++i)
            x[k]  -= x[i] * a_lu[i+k*n];
        x[k] /= a_lu[k+k*n];
    }
    /*
    Solve R * x = y.
    */
    for ( k = n-1; 0 <= k; --k )
    {
        x[k] /= a_lu[k+k*n];
        for ( i = 0; i < k; ++i )
            x[i] -= a_lu[i+k*n] * x[k];
    }
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _chebyshev1_exactness ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV1_EXACTNESS: monomial exactness for the Chebyshev1 integral.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points in the rule.
    Input, double X[N], the quadrature points.
    Input, double W[N], the quadrature weights.
    Input, int P_MAX, the maximum exponent.
    0 <= P_MAX.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ p_max = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	
    ityp e;
    int i;
    int p;
    ityp q;
    ityp s;
    ityp *v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( p = 0; p <= p_max; ++p )
    {
        s = chebyshev1_integral ( p );
        for ( i = 0; i < n; ++i )
            v[i] = pow ( x[i], p );

        q = r8vec_dot_product ( n, w, v );
        e = s == 0.00 ? fabs(q): fabs ( q - s ) / fabs ( s );
    }

    free ( v );
    return NULL;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _chebyshev1_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV1_INTEGRAL evaluates a monomial Chebyshev type 1 integral.
  Discussion:
    The integral:
      integral ( -1 <= x <= +1 ) x^n / sqrt ( 1 - x^2 ) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2008
  Author:
    John Burkardt
  Parameters:
    Input, int EXPON, the exponent.
    Output, double CHEBYSHEV1_INTEGRAL, the value of the exact integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ expon = *(dim_typ *) data;
	
    ityp bot;
    ityp exact;
    ityp top;
    /*
    Get the exact value of the integral.
    */
    if ( ( expon % 2 ) == 0 )
    {
        top = 1;
        bot = 1;
        for (dim_typ i = 2; i <= expon; i += 2 )
        {
            top *= ( i - 1 );
            bot *=   i;
        }
        exact = M_PI * ( ityp ) ( top ) / ( ityp ) ( bot );
    }
    else
        exact = 0.00;
        
    result = exact;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _chebyshev2_exactness ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV2_EXACTNESS: monomial exactness for the Chebyshev2 integral.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points in the rule.
    Input, double X[N], the quadrature points.
    Input, double W[N], the quadrature weights.
    Input, int P_MAX, the maximum exponent.
    0 <= P_MAX.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ p_max = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	
    ityp e;
    dim_typ i, p;
    ityp q;
    ityp s;
    ityp *v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( p = 0; p <= p_max; p++ )
    {
        s = chebyshev2_integral ( p );

        for ( i = 0; i < n; i++ )
            v[i] = pow ( x[i], p );

        q = r8vec_dot_product ( n, w, v );
        e = s == 0.00 ? fabs(q) : fabs ( q - s ) / fabs ( s );

    }

    free ( v );
    return NULL;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _chebyshev2_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV2_INTEGRAL evaluates a monomial Chebyshev type 2 integral.
  Discussion:
    The integral:
      integral ( -1 <= x <= +1 ) x^n * sqrt ( 1 - x^2 ) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2008
  Author:
    John Burkardt
  Parameters:
    Input, int EXPON, the exponent.
    Output, double CHEBYSHEV2_INTEGRAL, the value of the exact integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ expon = *(dim_typ *) data;
	
  ityp bot;
  ityp exact;
  ityp top;
/*
  Get the exact value of the integral.
*/
    if ( ( expon % 2 ) == 0 )
    {
        top = bot = 1;
        for (dim_typ i = 2; i <= expon; i += 2 )
        {
            top *= ( i - 1 );
            bot *=   i;
        }
        bot *= ( ityp ) ( expon + 2 );
        exact = M_PI * ( ityp ) ( top ) / ( ityp ) ( bot );
    }
    else
        exact = 0.00;
        
    result = exact;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_cosd ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_COSD returns the cosine of an angle given in degrees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, double DEGREES, the angle in degrees.
    Output, double R8_COSD, the cosine of the angle.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp degrees = *(ityp *) data;
	
	result = cos ( M_PI * ( degrees / 180.0 ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cotd ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_COTD returns the cotangent of an angle given in degrees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, double DEGREES, the angle in degrees.
    Output, double R8_COTD, the cotangent of the angle.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp degrees = *(ityp *) data;
	
    const register ityp radians = M_PI * ( degrees / 180.0 );
    
    result = cos ( radians ) / sin ( radians );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_cscd ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CSCD returns the cosecant of an angle given in degrees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, double DEGREES, the angle in degrees.
    Output, double R8_CSCD, the cosecant of the angle.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp degrees = *(ityp *) data;
	
	result = 1.00 / sin ( M_PI * ( M_PI / 180.00 ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_secd ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SECD returns the secant of an angle given in degrees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, double DEGREES, the angle in degrees.
    Output, double R8_SECD, the secant of the angle.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp degrees = *(ityp *) data;
	
	result = 1.00 / cos ( M_PI * ( degrees / 180.00 ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_sind (  void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SIND returns the sine of an angle given in degrees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, double DEGREES, the angle in degrees.
    Output, double R8_SIND, the sine of the angle.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp degrees = *(ityp *) data;
	
	result = sin ( M_PI * ( degrees / 180.00 ));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_tand ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_TAND returns the tangent of an angle given in degrees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, double DEGREES, the angle in degrees.
    Output, double R8_TAND, the tangent of the angle.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp degrees = *(ityp *) data;
	
    const register ityp radians = M_PI * ( degrees / 180.00 );
    
    result = sin ( radians ) / cos ( radians );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _c4_cos ( void * data)
/******************************************************************************/
/*
  Purpose:
    C4_COS evaluates the cosine of a C4 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, complex Z, the argument.
    Output, complex C4_COS, the cosine of Z.
*/
{
	static complex result = I;
	
	const register complex z = *(complex *) data;
	
    ityp x = creal ( z );
    ityp y = cimag ( z );
    
    result = ( complex ) ( cos ( x ) * cosh ( y ), - sin ( x ) * sinh ( y ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _c4_sin ( void * data)
/******************************************************************************/
/*
  Purpose:
    C4_SIN evaluates the sine of a C4 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, complex Z, the argument.
    Output, complex C4_SIN, the sine of Z.
*/
{
	static complex result = I;
	
	const register complex z = *(complex *) data;
	
    ityp x = creal ( z );
    ityp y = cimag ( z );
    
    result = ( complex ) ( sin ( x ) * cosh ( y ), cos ( x ) * sinh ( y ) );
    return &result;
}

/******************************************************************************/
 __MATHSUITE __JBURKARDT void *   _r8_rand ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_RAND is a portable pseudorandom number generator.
  Discussion:
    This pseudo-random number generator is portable amoung a wide
    variety of computers.  It is undoubtedly not as good as many
    readily available installation dependent versions, and so this
    routine is not recommended for widespread usage.  Its redeeming
    feature is that the exact same random numbers (to within final round-
    off error) can be generated from machine to machine.  Thus, programs
    that make use of random numbers can be easily transported to and
    checked in a new environment.
    The random numbers are generated by the linear congruential
    method described by Knuth in seminumerical methods (p.9),
    addison-wesley, 1969.  Given the i-th number of a pseudo-random
    sequence, the i+1 -st number is generated from
      x(i+1) = (a*x(i) + c) mod m,
    where here m = 2^22 = 4194304, c = 1731 and several suitable values
    of the multiplier a are discussed below.  Both the multiplier a and
    random number x are represented in double precision as two 11-bit
    words.  The constants are chosen so that the period is the maximum
    possible, 4194304.
    In order that the same numbers be generated from machine to
    machine, it is necessary that 23-bit ints be reducible modulo
    2^11 exactly, that 23-bit ints be added exactly, and that 11-bit
    ints be multiplied exactly.  Furthermore, if the restart option
    is used (where r is between 0 and 1), then the product r*2^22 =
    r*4194304 must be correct to the nearest int.
    The first four random numbers should be
      0.0004127026,
      0.6750836372,
      0.1614754200,
      0.9086198807.
    The tenth random number is
      0.5527787209.
    The hundredth random number is
      0.3600893021.
    The thousandth number should be
      0.2176990509.
    In order to generate several effectively independent sequences
    with the same generator, it is necessary to know the random number
    for several widely spaced calls.  The I-th random number times 2^22,
    where I=K*P/8 and P is the period of the sequence (P = 2^22), is
    still of the form L*P/8.  In particular we find the I-th random
    number multiplied by 2^22 is given by
      I   =  0  1*p/8  2*p/8  3*p/8  4*p/8  5*p/8  6*p/8  7*p/8  8*p/8
      RAND=  0  5*p/8  2*p/8  7*p/8  4*p/8  1*p/8  6*p/8  3*p/8  0
    thus the 4*P/8 = 2097152 random number is 2097152/2^22.
    Several multipliers have been subjected to the spectral test
 (see Knuth, p. 82).  Four suitable multipliers roughly in order of
    goodness according to the spectral test are
      3146757 = 1536*2048 + 1029 = 2^21 + 2^20 + 2^10 + 5
      2098181 = 1024*2048 + 1029 = 2^21 + 2^10 + 5
      3146245 = 1536*2048 +  517 = 2^21 + 2^20 + 2^9 + 5
      2776669 = 1355*2048 + 1629 = 5^9 + 7^7 + 1
    In the table below log10(NU(I)) gives roughly the number of
    random decimal digits in the random numbers considered I at a time.
    C is the primary measure of goodness.  In both cases bigger is better.
                     log10 nu(i)              c(i)
         a       i=2  i=3  i=4  i=5    i=2  i=3  i=4  i=5
      3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
      2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
      3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
      2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
     best
      possible   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    C version by John Burkardt.
  Parameters:
    Input, real R, determines the action.
    * R = 0.0, the next random number of the sequence is generated.
    * R < 0.0, the last generated number will be returned for
    possible use in a restart procedure.
    * R > 0.0, the sequence of random numbers will start with the
    seed ( R mod 1 ).  This seed is also returned as the value of
    r8_RAND provided the arithmetic is done exactly.
    Output, real r8_RAND, a pseudo-random number between 0.0 and 1.0.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp r = *(ityp *) data;
	
    static dim_typ ia0 = 1029;
    static dim_typ ia1 = 1536;
    static dim_typ ia1ma0 = 507;
    static dim_typ ic = 1731;
    static dim_typ ix0 = 0;
    static dim_typ ix1 = 0;
    dim_typ iy0;
    dim_typ iy1;
    ityp value;

    if ( r == 0.0E+00 )
    {
        iy0 = ia0 * ix0;
        iy1 = ia1 * ix1 + ia1ma0 * ( ix0 - ix1 ) + iy0;
        iy0 = iy0 + ic;
        ix0 = ( iy0 % 2048 );
        iy1 = iy1 + ( iy0 - ix0 ) / 2048;
        ix1 = ( iy1 % 2048 );
    }

    if ( 0.00 < r )
    {
        ix1 = ( dim_typ ) ( r8_mod ( r, 1.0E+00 ) * 4194304.0 + 0.5E+00 );
        ix0 = ( ix1 % 2048 );
        ix1 = ( ix1 - ix0 ) / 2048;
    }

	result = ( ityp ) ( ix1 * 2048 + ix0 ) / 4194304.0E+0;
    return &result;
}

/******************************************************************************/
 __MATHSUITE __JBURKARDT void *   _r8_randgs ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_RANDGS generates a normally distributed random number.
  Discussion:
    This function generate a normally distributed random number, that is,
    it generates random numbers with a Gaussian distribution.  These
    random numbers are not exceptionally good, especially in the tails
    of the distribution, but this implementation is simple and suitable
    for most applications.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    C version by John Burkardt.
  Reference:
    Richard Hamming,
    Numerical Methods for Scientists and Engineers,
    Dover, 1986,
    ISBN: 0486652416,
    LC: QA297.H28.
  Parameters:
    Input, ityp XMEAN, the mean of the Gaussian distribution.
    Input, ityp SD, the standard deviation of the Gaussian function.
    Output, ityp r8_RANDGS, a normally distributed random number.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp xmean = a_data[0];
	const register ityp sd = a_data[1];
	
    ityp value = - 6.0E+00;
    dim_typ i;
    #pragma omp parallel for num_threads(12)
    for ( i = 1; i <= 12; ++i )
        value += r8_rand ( 0.0E+00 );
        
    result = xmean + sd * value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_admp ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_ADMP: modulus and phase of the derivative of the Airy function.
  Description:
    This function must only be called when X <= -1.0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double *AMPL, *PHI, the modulus and phase of the
    derivative of the Airy function.
*/
{
	const it2pit * const s_data = data;
	const register ityp x = s_data->a0;
	ityp * ampl = s_data->a1;
	ityp * phi = s_data->a2;
	
    static ityp an20cs[57] =
    {
        0.0126732217145738027154610751034240,
        -0.0005212847072615621184780942309478,
        -0.0000052672111140370429809074052969,
        -0.0000001628202185026483752632460680,
        -0.0000000090991442687371386325973075,
        -0.0000000007438647126242192890685403,
        -0.0000000000795494751591469486122822,
        -0.0000000000104050944288303742803960,
        -0.0000000000015932425598414551523990,
        -0.0000000000002770648272341913946674,
        -0.0000000000000535342629237606295104,
        -0.0000000000000113061541781728314051,
        -0.0000000000000025772190078943167788,
        -0.0000000000000006278033116032485076,
        -0.0000000000000001621295400189939757,
        -0.0000000000000000440992985240675353,
        -0.0000000000000000125655516553258972,
        -0.0000000000000000037336906988015204,
        -0.0000000000000000011524626926724671,
        -0.0000000000000000003683081499099144,
        -0.0000000000000000001215206965331797,
        -0.0000000000000000000412916177724016,
        -0.0000000000000000000144177364239347,
        -0.0000000000000000000051631842875864,
        -0.0000000000000000000018931242668250,
        -0.0000000000000000000007096054668569,
        -0.0000000000000000000002715406646904,
        -0.0000000000000000000001059486979400,
        -0.0000000000000000000000421030035685,
        -0.0000000000000000000000170233781664,
        -0.0000000000000000000000069966677028,
        -0.0000000000000000000000029206643813,
        -0.0000000000000000000000012373128203,
        -0.0000000000000000000000005315871095,
        -0.0000000000000000000000002314622618,
        -0.0000000000000000000000001020779922,
        -0.0000000000000000000000000455706227,
        -0.0000000000000000000000000205831071,
        -0.0000000000000000000000000094015189,
        -0.0000000000000000000000000043405874,
        -0.0000000000000000000000000020247792,
        -0.0000000000000000000000000009539214,
        -0.0000000000000000000000000004537234,
        -0.0000000000000000000000000002178016,
        -0.0000000000000000000000000001054823,
        -0.0000000000000000000000000000515242,
        -0.0000000000000000000000000000253763,
        -0.0000000000000000000000000000125983,
        -0.0000000000000000000000000000063030,
        -0.0000000000000000000000000000031771,
        -0.0000000000000000000000000000016131,
        -0.0000000000000000000000000000008248,
        -0.0000000000000000000000000000004246,
        -0.0000000000000000000000000000002200,
        -0.0000000000000000000000000000001147,
        -0.0000000000000000000000000000000602,
        -0.0000000000000000000000000000000318
    };
    static ityp an21cs[60] =
    {
        0.0198313155263169394420342483165643,
        -0.0029376249067087533460593745594484,
        -0.0001136260695958195549872611137182,
        -0.0000100554451087156009750981645918,
        -0.0000013048787116563250421785598252,
        -0.0000002123881993150664830666079609,
        -0.0000000402270833384269040347850109,
        -0.0000000084996745953161799142201792,
        -0.0000000019514839426178614099532934,
        -0.0000000004783865343840384282992480,
        -0.0000000001236733992099450501137105,
        -0.0000000000334137486398754232219789,
        -0.0000000000093702823540766329897780,
        -0.0000000000027130128156139564687240,
        -0.0000000000008075953800583479535949,
        -0.0000000000002463214304700125252160,
        -0.0000000000000767655689109321564410,
        -0.0000000000000243882598807354919791,
        -0.0000000000000078831466358760308462,
        -0.0000000000000025882400995585864077,
        -0.0000000000000008619457862945690828,
        -0.0000000000000002907994739663128534,
        -0.0000000000000000992846796122890484,
        -0.0000000000000000342720229187774480,
        -0.0000000000000000119511048205515026,
        -0.0000000000000000042069729043678359,
        -0.0000000000000000014939697762818400,
        -0.0000000000000000005348981161589517,
        -0.0000000000000000001929877577826238,
        -0.0000000000000000000701313701018203,
        -0.0000000000000000000256585738509682,
        -0.0000000000000000000094475894562734,
        -0.0000000000000000000034996401941465,
        -0.0000000000000000000013037622466397,
        -0.0000000000000000000004883334163346,
        -0.0000000000000000000001838477586152,
        -0.0000000000000000000000695527324058,
        -0.0000000000000000000000264351910209,
        -0.0000000000000000000000100918094655,
        -0.0000000000000000000000038688924289,
        -0.0000000000000000000000014892036525,
        -0.0000000000000000000000005754342426,
        -0.0000000000000000000000002231725971,
        -0.0000000000000000000000000868607480,
        -0.0000000000000000000000000339220403,
        -0.0000000000000000000000000132910128,
        -0.0000000000000000000000000052239309,
        -0.0000000000000000000000000020594383,
        -0.0000000000000000000000000008142614,
        -0.0000000000000000000000000003228473,
        -0.0000000000000000000000000001283529,
        -0.0000000000000000000000000000511622,
        -0.0000000000000000000000000000204451,
        -0.0000000000000000000000000000081901,
        -0.0000000000000000000000000000032886,
        -0.0000000000000000000000000000013235,
        -0.0000000000000000000000000000005338,
        -0.0000000000000000000000000000002158,
        -0.0000000000000000000000000000000874,
        -0.0000000000000000000000000000000355
    };
    static ityp an22cs[74] =
    {
        0.0537418629629794329091103360917783,
        -0.0126661435859883193466312085036450,
        -0.0011924334106593006840848916913681,
        -0.0002032327627275654552687155176363,
        -0.0000446468963075163979516164905945,
        -0.0000113359036053123490416997893086,
        -0.0000031641352378546107356671355827,
        -0.0000009446708886148939120888532442,
        -0.0000002966562236471765527900905456,
        -0.0000000969118892024367799908661433,
        -0.0000000326822538653274091533072559,
        -0.0000000113144618963583865900447294,
        -0.0000000040042691001741501738278050,
        -0.0000000014440333683907423778522199,
        -0.0000000005292853746152611585663541,
        -0.0000000001967763373707889528245726,
        -0.0000000000740800095755849858816731,
        -0.0000000000282016314294661982842740,
        -0.0000000000108440066463128331337590,
        -0.0000000000042074800682644236920617,
        -0.0000000000016459149670634819724739,
        -0.0000000000006486826705121018896077,
        -0.0000000000002574095003354105832300,
        -0.0000000000001027889029407822132143,
        -0.0000000000000412845827195222720128,
        -0.0000000000000166711029332862509726,
        -0.0000000000000067656696165608023403,
        -0.0000000000000027585448232693576823,
        -0.0000000000000011296397915297168938,
        -0.0000000000000004644848225457314333,
        -0.0000000000000001917198035033912928,
        -0.0000000000000000794197570111893530,
        -0.0000000000000000330116492300368930,
        -0.0000000000000000137658057726549714,
        -0.0000000000000000057578093720012791,
        -0.0000000000000000024152700858632017,
        -0.0000000000000000010159301700933666,
        -0.0000000000000000004284434955330055,
        -0.0000000000000000001811344052168016,
        -0.0000000000000000000767602045619422,
        -0.0000000000000000000326026346758614,
        -0.0000000000000000000138773806682627,
        -0.0000000000000000000059191627103729,
        -0.0000000000000000000025297256431944,
        -0.0000000000000000000010832077293819,
        -0.0000000000000000000004646674880404,
        -0.0000000000000000000001996797783865,
        -0.0000000000000000000000859524108705,
        -0.0000000000000000000000370584152073,
        -0.0000000000000000000000160027503479,
        -0.0000000000000000000000069208124999,
        -0.0000000000000000000000029974448994,
        -0.0000000000000000000000013000356362,
        -0.0000000000000000000000005646100942,
        -0.0000000000000000000000002455341103,
        -0.0000000000000000000000001069119686,
        -0.0000000000000000000000000466095090,
        -0.0000000000000000000000000203441579,
        -0.0000000000000000000000000088900866,
        -0.0000000000000000000000000038891813,
        -0.0000000000000000000000000017032637,
        -0.0000000000000000000000000007467295,
        -0.0000000000000000000000000003277097,
        -0.0000000000000000000000000001439618,
        -0.0000000000000000000000000000633031,
        -0.0000000000000000000000000000278620,
        -0.0000000000000000000000000000122743,
        -0.0000000000000000000000000000054121,
        -0.0000000000000000000000000000023884,
        -0.0000000000000000000000000000010549,
        -0.0000000000000000000000000000004663,
        -0.0000000000000000000000000000002063,
        -0.0000000000000000000000000000000913,
        -0.0000000000000000000000000000000405
    };
    static ityp aph0cs[53] =
    {
        -0.0855849241130933256920124260179491,
        0.0011214378867065260735786722471124,
        0.0000042721029353664113951573742015,
        0.0000000817607381483243644018062323,
        0.0000000033907645000492724207816418,
        0.0000000002253264422619113939845276,
        0.0000000000206284209229015251256990,
        0.0000000000023858762828130887627258,
        0.0000000000003301618105886705480628,
        0.0000000000000527009648508328581123,
        0.0000000000000094555482203813492868,
        0.0000000000000018709426951344836908,
        0.0000000000000004023980041825392741,
        0.0000000000000000930192879258983167,
        0.0000000000000000229038635402379945,
        0.0000000000000000059634359822083386,
        0.0000000000000000016320279659403399,
        0.0000000000000000004671145658861339,
        0.0000000000000000001392334415363502,
        0.0000000000000000000430642670285155,
        0.0000000000000000000137781416318755,
        0.0000000000000000000045476710480396,
        0.0000000000000000000015448420203026,
        0.0000000000000000000005389770551212,
        0.0000000000000000000001927726737155,
        0.0000000000000000000000705659320166,
        0.0000000000000000000000263985084827,
        0.0000000000000000000000100791301805,
        0.0000000000000000000000039228928481,
        0.0000000000000000000000015547422955,
        0.0000000000000000000000006268306372,
        0.0000000000000000000000002568563962,
        0.0000000000000000000000001068858883,
        0.0000000000000000000000000451347253,
        0.0000000000000000000000000193267262,
        0.0000000000000000000000000083865369,
        0.0000000000000000000000000036857386,
        0.0000000000000000000000000016396202,
        0.0000000000000000000000000007379298,
        0.0000000000000000000000000003358392,
        0.0000000000000000000000000001544891,
        0.0000000000000000000000000000718013,
        0.0000000000000000000000000000337026,
        0.0000000000000000000000000000159710,
        0.0000000000000000000000000000076382,
        0.0000000000000000000000000000036855,
        0.0000000000000000000000000000017935,
        0.0000000000000000000000000000008800,
        0.0000000000000000000000000000004353,
        0.0000000000000000000000000000002170,
        0.0000000000000000000000000000001090,
        0.0000000000000000000000000000000551,
        0.0000000000000000000000000000000281
    };
    static ityp aph1cs[58] =
    {
        -0.1024172908077571694021123321813917,
        0.0071697275146591248047211649144704,
        0.0001209959363122328589813856491397,
        0.0000073361512841219912080297845684,
        0.0000007535382954271607069982903869,
        0.0000001041478171741301926885109155,
        0.0000000174358728518545691858907606,
        0.0000000033399795033346451660184961,
        0.0000000007073075174363527083399508,
        0.0000000001619187515189773266792272,
        0.0000000000394539981881954889879668,
        0.0000000000101192281734227133292631,
        0.0000000000027092778259520332198030,
        0.0000000000007523806418422548885854,
        0.0000000000002156368733008966357328,
        0.0000000000000635282777126068410174,
        0.0000000000000191756972641501729345,
        0.0000000000000059143072446464891558,
        0.0000000000000018597128517275028357,
        0.0000000000000005950444923946103668,
        0.0000000000000001934229956430180252,
        0.0000000000000000637843021489504324,
        0.0000000000000000213127290087312393,
        0.0000000000000000072081380656728500,
        0.0000000000000000024652494144769247,
        0.0000000000000000008519110570266154,
        0.0000000000000000002972384468491170,
        0.0000000000000000001046426648811446,
        0.0000000000000000000371493036347327,
        0.0000000000000000000132923247793472,
        0.0000000000000000000047912837925909,
        0.0000000000000000000017390619859336,
        0.0000000000000000000006353585173501,
        0.0000000000000000000002335643614263,
        0.0000000000000000000000863643881606,
        0.0000000000000000000000321123006944,
        0.0000000000000000000000120031540983,
        0.0000000000000000000000045091488699,
        0.0000000000000000000000017020228580,
        0.0000000000000000000000006453744630,
        0.0000000000000000000000002457788564,
        0.0000000000000000000000000939897684,
        0.0000000000000000000000000360863150,
        0.0000000000000000000000000139077884,
        0.0000000000000000000000000053797184,
        0.0000000000000000000000000020882551,
        0.0000000000000000000000000008133371,
        0.0000000000000000000000000003178080,
        0.0000000000000000000000000001245700,
        0.0000000000000000000000000000489742,
        0.0000000000000000000000000000193099,
        0.0000000000000000000000000000076349,
        0.0000000000000000000000000000030269,
        0.0000000000000000000000000000012032,
        0.0000000000000000000000000000004795,
        0.0000000000000000000000000000001915,
        0.0000000000000000000000000000000767,
        0.0000000000000000000000000000000308
    };
    static ityp aph2cs[72] =
    {
        -0.2057088719781465106973648665602125,
        0.0422196961357771921673114980369460,
        0.0020482560511207275042660577813334,
        0.0002607800735165005631187879922652,
        0.0000474824268004728875381750519293,
        0.0000105102756431611743473630026955,
        0.0000026353534014667945109314041983,
        0.0000007208824863499147299790783731,
        0.0000002103236664473352859749477082,
        0.0000000644975634555295598437362273,
        0.0000000205802377264368507978116888,
        0.0000000067836273920906428963513918,
        0.0000000022974015284009400168343792,
        0.0000000007961306765491187534883226,
        0.0000000002813860609741591719003632,
        0.0000000001011749056931973922841793,
        0.0000000000369306737952476559097060,
        0.0000000000136615066127098031778842,
        0.0000000000051142751416045045119388,
        0.0000000000019351688931706516247975,
        0.0000000000007393606916493224217271,
        0.0000000000002849792219222743597555,
        0.0000000000001107280782459648335733,
        0.0000000000000433412199370134633169,
        0.0000000000000170800825265670367471,
        0.0000000000000067733080195631114673,
        0.0000000000000027016904789262414108,
        0.0000000000000010834720751810782141,
        0.0000000000000004367060312970286167,
        0.0000000000000001768511738053366608,
        0.0000000000000000719359213093645717,
        0.0000000000000000293823610002933154,
        0.0000000000000000120482811525848357,
        0.0000000000000000049586659491091389,
        0.0000000000000000020479438315847217,
        0.0000000000000000008486019944410629,
        0.0000000000000000003527351765384506,
        0.0000000000000000001470563996804903,
        0.0000000000000000000614817826902188,
        0.0000000000000000000257737706565077,
        0.0000000000000000000108323903590042,
        0.0000000000000000000045638898024998,
        0.0000000000000000000019273635403662,
        0.0000000000000000000008157668569775,
        0.0000000000000000000003460202828346,
        0.0000000000000000000001470726482427,
        0.0000000000000000000000626356074088,
        0.0000000000000000000000267261292780,
        0.0000000000000000000000114246948763,
        0.0000000000000000000000048923460516,
        0.0000000000000000000000020985807810,
        0.0000000000000000000000009016618807,
        0.0000000000000000000000003880129464,
        0.0000000000000000000000001672282170,
        0.0000000000000000000000000721790800,
        0.0000000000000000000000000311982573,
        0.0000000000000000000000000135035015,
        0.0000000000000000000000000058524861,
        0.0000000000000000000000000025397686,
        0.0000000000000000000000000011035457,
        0.0000000000000000000000000004800788,
        0.0000000000000000000000000002090956,
        0.0000000000000000000000000000911743,
        0.0000000000000000000000000000397998,
        0.0000000000000000000000000000173923,
        0.0000000000000000000000000000076083,
        0.0000000000000000000000000000033316,
        0.0000000000000000000000000000014604,
        0.0000000000000000000000000000006407,
        0.0000000000000000000000000000002814,
        0.0000000000000000000000000000001237,
        0.0000000000000000000000000000000544
    };
    ityp eta;
    static dim_typ nan20 = 0;
    static dim_typ nan21 = 0;
    static dim_typ nan22 = 0;
    static dim_typ naph0 = 0;
    static dim_typ naph1 = 0;
    static dim_typ naph2 = 0;
    static ityp pi34 = 2.35619449019234492884698253745962716313;
    ityp sqrtx;
    static ityp xsml = 0.00;
    ityp z;

    if ( nan20 == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        nan20 = r8_inits ( an20cs, 57, eta );
        nan21 = r8_inits ( an21cs, 60, eta );
        nan22 = r8_inits ( an22cs, 74, eta );
        naph0 = r8_inits ( aph0cs, 53, eta );
        naph1 = r8_inits ( aph1cs, 58, eta );
        naph2 = r8_inits ( aph2cs, 72, eta );
        xsml = - pow ( 128.00 / r8_mach ( 3 ), 0.3333 );
    }

    if ( x < xsml )
    {
        z = 1.00;
        *ampl = 0.3125 + r8_csevl ( z, an20cs, nan20 );
        *phi = - 0.625 + r8_csevl ( z, aph0cs, naph0 );
    }
    else if ( x < - 4.00 )
    {
        z = 128.0 / x / x / x + 1.00;
        *ampl = 0.3125 + r8_csevl ( z, an20cs, nan20 );
        *phi = - 0.625 + r8_csevl ( z, aph0cs, naph0 );
    }
    else if ( x < - 2.0 )
    {
        z = ( 128. / x / x / x + 9.00 ) / 7.00;
        *ampl = 0.3125 + r8_csevl ( z, an21cs, nan21 );
        *phi = - 0.625 + r8_csevl ( z, aph1cs, naph1 );
    }
    else if ( x <= - 1.00 )
    {
        z = ( 16.0 / x / x / x + 9.00 ) / 7.00;
        *ampl = 0.3125 + r8_csevl ( z, an22cs, nan22 );
        *phi = - 0.625 + r8_csevl ( z, aph2cs, naph2 );
    }
    else
        return NULL;

    sqrtx = sqrt ( - x );
    *ampl = sqrt ( *ampl * sqrtx );
    *phi = pi34 - x * sqrtx * *phi;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_ai ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_AI evaluates the Airy function Ai of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_AI, the Airy function Ai of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp aifcs[13] =
    {
        -0.37971358496669997496197089469414E-01,
        +0.59191888537263638574319728013777E-01,
        +0.98629280577279975365603891044060E-03,
        +0.68488438190765667554854830182412E-05,
        +0.25942025962194713019489279081403E-07,
        +0.61766127740813750329445749697236E-10,
        +0.10092454172466117901429556224601E-12,
        +0.12014792511179938141288033225333E-15,
        +0.10882945588716991878525295466666E-18,
        +0.77513772196684887039238400000000E-22,
        +0.44548112037175638391466666666666E-25,
        +0.21092845231692343466666666666666E-28,
        +0.83701735910741333333333333333333E-32
    };
    static ityp aigcs[13] =
    {
        +0.18152365581161273011556209957864E-01,
        +0.21572563166010755534030638819968E-01,
        +0.25678356987483249659052428090133E-03,
        +0.14265214119792403898829496921721E-05,
        +0.45721149200180426070434097558191E-08,
        +0.95251708435647098607392278840592E-11,
        +0.13925634605771399051150420686190E-13,
        +0.15070999142762379592306991138666E-16,
        +0.12559148312567778822703205333333E-19,
        +0.83063073770821340343829333333333E-23,
        +0.44657538493718567445333333333333E-26,
        +0.19900855034518869333333333333333E-29,
        +0.74702885256533333333333333333333E-33
    };
    static dim_typ naif = 0;
    static dim_typ naig = 0;
    ityp theta;
    ityp value;
    static ityp x3sml = 0.00;
    ityp xm;
    static ityp xmax = 0.00;
    ityp z;

    if ( naif == 0 )
    {
        naif = r8_inits ( aifcs, 13, 0.1 * r8_mach ( 3 ) );
        naig = r8_inits ( aigcs, 13, 0.1 * r8_mach ( 3 ) );
        x3sml = pow ( r8_mach ( 3 ), 0.3334 );
        xmax = pow ( - 1.5 * log ( r8_mach ( 1 ) ), 0.6667 );
        xmax -= xmax * log ( xmax ) /( 4.00 * xmax * sqrt ( xmax ) + 1.00 ) - 0.01;
    }

    if ( x < - 1.00 )
    {
        r8_aimp ( x, &xm, &theta );
        value = xm * cos ( theta );
    }
    else if ( abs ( x ) <= x3sml )
    {
        z = 0.00;
        value = 0.375 + ( r8_csevl ( z, aifcs, naif )- x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
    }
    else if ( x <= 1.0 )
    {
        z = x * x * x;
        value = 0.375 + ( r8_csevl ( z, aifcs, naif )- x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
    }
    else if ( x <= xmax )
    	value = r8_aie ( x ) * exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
    else
    	value = 0.00;

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_aid ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_AID evaluates the derivative of the Airy function Ai of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_AID, the derivative of the Airy function
    Ai of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp aifcs[13] =
    {
        0.105274612265314088088970057325134114,
        0.011836136281529978442889292583980840,
        0.000123281041732256643051689242469164,
        0.000000622612256381399016825658693579,
        0.000000001852988878441452950548140821,
        0.000000000003633288725904357915995625,
        0.000000000000005046217040440664768330,
        0.000000000000000005223816555471480985,
        0.000000000000000000004185745090748989,
        0.000000000000000000000002672887324883,
        0.000000000000000000000000001392128006,
        0.000000000000000000000000000000602653,
        0.000000000000000000000000000000000220
    };
    static ityp aigcs[13] =
    {
        0.0212338781509186668523122276848937,
        0.0863159303352144067524942809461604,
        0.0017975947203832313578033963225230,
        0.0000142654998755506932526620687495,
        0.0000000594379952836832010488787064,
        0.0000000001524033664794478945214786,
        0.0000000000002645876603490435305100,
        0.0000000000000003315624296815020591,
        0.0000000000000000003139789757594792,
        0.0000000000000000000002325767379040,
        0.0000000000000000000000001384384231,
        0.0000000000000000000000000000676629,
        0.0000000000000000000000000000000276
    };
    ityp eta;
    static dim_typ naif = 0;
    static dim_typ naig = 0;
    ityp phi;
    ityp value;
    ityp x2;
    static ityp x2sml = 0.00;
    ityp x3;
    static ityp x3sml = 0.00;
    ityp xn;

    if ( naif == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        naif = r8_inits ( aifcs, 13, eta );
        naig = r8_inits ( aigcs, 13, eta );
        x3sml = pow ( r8_mach ( 3 ), 0.3334 );
        x2sml = sqrt ( r8_mach ( 3 ) );
    }

    if ( x < - 1.00 )
    {
        r8_admp ( x, &xn, &phi );
        value = xn * cos ( phi );
    }
    else if ( abs ( x ) <= x2sml )
    {
        x2 = x3 = 0.00;
        value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) )- r8_csevl ( x3, aigcs, naig ) ) - 0.25;
    }
    else if ( abs ( x ) <= x3sml )
    {
        x2 = x * x;
        x3 = 0.00;
        value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) )- r8_csevl ( x3, aigcs, naig ) ) - 0.25;
    }
    else if ( x <= 1.00 )
    {
        x2 = x * x;
        x3 = x * x * x;
        value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) )- r8_csevl ( x3, aigcs, naig ) ) - 0.25;
    }
    else
        value = r8_aide ( x ) * exp ( - 2.0 * x * sqrt ( x ) / 3.00 );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_aide ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_AIDE: exponentially scaled derivative, Airy function Ai of an R8 argument.
  Discussion:
    if X <= 0,
      R8_AIDE ( X ) = R8_AID ( X )
    else
      R8_AIDE ( X ) = R8_AID ( X ) * exp ( 2/3 * X^(3/2) )
    Thanks to Aleksandra Piper for pointing out a correction involving
    the computation of Z in the last two cases, 02 February 2012.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 February 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_AIDE, the exponentially scaled derivative of
    the Airy function Ai of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp aifcs[13] =
    {
        0.105274612265314088088970057325134114,
        0.011836136281529978442889292583980840,
        0.000123281041732256643051689242469164,
        0.000000622612256381399016825658693579,
        0.000000001852988878441452950548140821,
        0.000000000003633288725904357915995625,
        0.000000000000005046217040440664768330,
        0.000000000000000005223816555471480985,
        0.000000000000000000004185745090748989,
        0.000000000000000000000002672887324883,
        0.000000000000000000000000001392128006,
        0.000000000000000000000000000000602653,
        0.000000000000000000000000000000000220
    };
    ityp aigcs[13] =
    {
        0.0212338781509186668523122276848937,
        0.0863159303352144067524942809461604,
        0.0017975947203832313578033963225230,
        0.0000142654998755506932526620687495,
        0.0000000594379952836832010488787064,
        0.0000000001524033664794478945214786,
        0.0000000000002645876603490435305100,
        0.0000000000000003315624296815020591,
        0.0000000000000000003139789757594792,
        0.0000000000000000000002325767379040,
        0.0000000000000000000000001384384231,
        0.0000000000000000000000000000676629,
        0.0000000000000000000000000000000276
    };
    ityp aip1cs[57] =
    {
        0.0358865097808301537956710489261688,
        0.0114668575627764898572700883121766,
        -0.0007592073583861400301335647601603,
        0.0000869517610893841271948619434021,
        -0.0000128237294298591691789607600486,
        0.0000022062695681038336934376250420,
        -0.0000004222295185920749486945988432,
        0.0000000874686415726348479356130376,
        -0.0000000192773588418365388625693417,
        0.0000000044668460054492719699777137,
        -0.0000000010790108051948168015747466,
        0.0000000002700029446696248083071434,
        -0.0000000000696480108007915257318929,
        0.0000000000184489907003246687076806,
        -0.0000000000050027817358071698301149,
        0.0000000000013852243366012168297298,
        -0.0000000000003908218466657048253473,
        0.0000000000001121536072524563451273,
        -0.0000000000000326861522579502522443,
        0.0000000000000096619179010090805752,
        -0.0000000000000028934767442698434271,
        0.0000000000000008770086661150897069,
        -0.0000000000000002688046261195853754,
        0.0000000000000000832498823872342992,
        -0.0000000000000000260343254786947057,
        0.0000000000000000082159528142686287,
        -0.0000000000000000026150406704984940,
        0.0000000000000000008390563463261051,
        -0.0000000000000000002712685618629660,
        0.0000000000000000000883333375271942,
        -0.0000000000000000000289603206822333,
        0.0000000000000000000095562185928676,
        -0.0000000000000000000031727463569051,
        0.0000000000000000000010595576960768,
        -0.0000000000000000000003558253765402,
        0.0000000000000000000001201334680517,
        -0.0000000000000000000000407666883800,
        0.0000000000000000000000139016944446,
        -0.0000000000000000000000047628165730,
        0.0000000000000000000000016391265551,
        -0.0000000000000000000000005665491354,
        0.0000000000000000000000001966381969,
        -0.0000000000000000000000000685230229,
        0.0000000000000000000000000239706939,
        -0.0000000000000000000000000084166831,
        0.0000000000000000000000000029659364,
        -0.0000000000000000000000000010487947,
        0.0000000000000000000000000003721150,
        -0.0000000000000000000000000001324570,
        0.0000000000000000000000000000472976,
        -0.0000000000000000000000000000169405,
        0.0000000000000000000000000000060855,
        -0.0000000000000000000000000000021924,
        0.0000000000000000000000000000007920,
        -0.0000000000000000000000000000002869,
        0.0000000000000000000000000000001042,
        -0.0000000000000000000000000000000379
    };
    ityp aip2cs[37] =
    {
        0.0065457691989713756794276979067064,
        0.0023833724120774591992772552886923,
        -0.0000430700770220585862775012110584,
        0.0000015629125858629202330785369063,
        -0.0000000815417186162706965112501015,
        0.0000000054103738056935918208008783,
        -0.0000000004284130882614696528766222,
        0.0000000000389497962832286424862198,
        -0.0000000000039623161264979257658071,
        0.0000000000004428184214405989602353,
        -0.0000000000000536296527150689675318,
        0.0000000000000069649872139936028200,
        -0.0000000000000009619636286095319210,
        0.0000000000000001403454967784808032,
        -0.0000000000000000215097136525875715,
        0.0000000000000000034471230632678283,
        -0.0000000000000000005753907621819442,
        0.0000000000000000000997001165824168,
        -0.0000000000000000000178811436021458,
        0.0000000000000000000033110307923551,
        -0.0000000000000000000006315885529506,
        0.0000000000000000000001238666952364,
        -0.0000000000000000000000249324053394,
        0.0000000000000000000000051426030999,
        -0.0000000000000000000000010854236402,
        0.0000000000000000000000002341316852,
        -0.0000000000000000000000000515542099,
        0.0000000000000000000000000115758841,
        -0.0000000000000000000000000026479669,
        0.0000000000000000000000000006165328,
        -0.0000000000000000000000000001459931,
        0.0000000000000000000000000000351331,
        -0.0000000000000000000000000000085863,
        0.0000000000000000000000000000021297,
        -0.0000000000000000000000000000005358,
        0.0000000000000000000000000000001367,
        -0.0000000000000000000000000000000353
    };
    ityp eta;
    static dim_typ naif = 0;
    static dim_typ naig = 0;
    static dim_typ naip1 = 0;
    static dim_typ naip2 = 0;
    ityp phi;
    ityp sqrtx;
    ityp value;
    ityp x2;
    static ityp x2sml = 0.00;
    ityp x3;
    static ityp x32sml = 0.00;
    static ityp x3sml = 0.00;
    static ityp xbig = 0.00;
    ityp xn;
    ityp z;

    if ( naif == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        naif = r8_inits ( aifcs, 13, eta );
        naig = r8_inits ( aigcs, 13, eta );
        naip1 = r8_inits ( aip1cs, 57, eta );
        naip2 = r8_inits ( aip2cs, 37, eta );
        x2sml = sqrt ( eta );
        x3sml = pow ( eta, 0.3333 );
        x32sml = 1.3104 * x3sml * x3sml;
        xbig = pow ( r8_mach ( 2 ), 0.6666 );
    }

    if ( x < - 1.00 )
    {
        r8_admp ( x, &xn, &phi );
        value = xn * cos ( phi );
    }
    else if ( abs ( x ) < x2sml )
    {
        x2 = x3 = 0.00;
        value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) )- r8_csevl ( x3, aigcs, naig ) ) - 0.25;
        if ( x32sml < x )
            value *= exp ( 2.0 * x * sqrt ( x ) / 3.0 );
    }
    else if ( abs ( x ) < x3sml )
    {
        x2 = x * x;
        x3 = 0.00;
        value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) )- r8_csevl ( x3, aigcs, naig ) ) - 0.25;
        if ( x32sml < x )
            value *= exp ( 2.0 * x * sqrt ( x ) / 3.0 );
    }
    else if ( x <= 1.00 )
    {
        x2 = x * x;
        x3 = x * x;
        value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) )- r8_csevl ( x3, aigcs, naig ) ) - 0.25;
        if ( x32sml < x )
            value *= exp ( 2.0 * x * sqrt ( x ) / 3.0 );
    }
    else if ( x <= 4.00 )
    {
        sqrtx = sqrt ( x );
        z = ( 16.0  / ( x * sqrtx ) - 9.0 ) / 7.0;
        value = ( - 0.28125 - r8_csevl ( z, aip1cs, naip1 ) ) * sqrt ( sqrtx );
    }
    else if ( x < xbig )
    {
        sqrtx = sqrt ( x );
        z = 16.00  / ( x * sqrtx ) - 1.00;
        value = ( - 0.28125 - r8_csevl ( z, aip2cs, naip2 ) ) * sqrt ( sqrtx );
    }
    else
    {
        sqrtx = sqrt ( x );
        z = - 1.00;
        value = ( - 0.28125 - r8_csevl ( z, aip2cs, naip2 ) ) * sqrt ( sqrtx );
    }

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_aie ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_AIE evaluates the exponentially scaled Airy function Ai of an R8 argument.
  Discussion:
    if X <= 0,
      R8_AIE ( X ) = R8_AI ( X )
    else
      R8_AIE ( X ) = R8_AI ( X ) * exp ( 2/3 * X^(3/2) )
    Thanks to Aleksandra Piper for pointing out a correction involving a
    missing assignment to SQRTX, 27 January 2012.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_AIE, the exponentially scaled Airy function Ai of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp aifcs[13] =
    {
        -0.37971358496669997496197089469414E-01,
        +0.59191888537263638574319728013777E-01,
        +0.98629280577279975365603891044060E-03,
        +0.68488438190765667554854830182412E-05,
        +0.25942025962194713019489279081403E-07,
        +0.61766127740813750329445749697236E-10,
        +0.10092454172466117901429556224601E-12,
        +0.12014792511179938141288033225333E-15,
        +0.10882945588716991878525295466666E-18,
        +0.77513772196684887039238400000000E-22,
        +0.44548112037175638391466666666666E-25,
        +0.21092845231692343466666666666666E-28,
        +0.83701735910741333333333333333333E-32
    };
    static ityp aigcs[13] =
    {
        +0.18152365581161273011556209957864E-01,
        +0.21572563166010755534030638819968E-01,
        +0.25678356987483249659052428090133E-03,
        +0.14265214119792403898829496921721E-05,
        +0.45721149200180426070434097558191E-08,
        +0.95251708435647098607392278840592E-11,
        +0.13925634605771399051150420686190E-13,
        +0.15070999142762379592306991138666E-16,
        +0.12559148312567778822703205333333E-19,
        +0.83063073770821340343829333333333E-23,
        +0.44657538493718567445333333333333E-26,
        +0.19900855034518869333333333333333E-29,
        +0.74702885256533333333333333333333E-33
    };
    static ityp aip1cs[57] =
    {
        -0.2146951858910538455460863467778E-01,
        -0.7535382535043301166219720865565E-02,
        +0.5971527949026380852035388881994E-03,
        -0.7283251254207610648502368291548E-04,
        +0.1110297130739299666517381821140E-04,
        -0.1950386152284405710346930314033E-05,
        +0.3786973885159515193885319670057E-06,
        -0.7929675297350978279039072879154E-07,
        +0.1762247638674256075568420122202E-07,
        -0.4110767539667195045029896593893E-08,
        +0.9984770057857892247183414107544E-09,
        -0.2510093251387122211349867730034E-09,
        +0.6500501929860695409272038601725E-10,
        -0.1727818405393616515478877107366E-10,
        +0.4699378842824512578362292872307E-11,
        -0.1304675656297743914491241246272E-11,
        +0.3689698478462678810473948382282E-12,
        -0.1061087206646806173650359679035E-12,
        +0.3098414384878187438660210070110E-13,
        -0.9174908079824139307833423547851E-14,
        +0.2752049140347210895693579062271E-14,
        -0.8353750115922046558091393301880E-15,
        +0.2563931129357934947568636168612E-15,
        -0.7950633762598854983273747289822E-16,
        +0.2489283634603069977437281175644E-16,
        -0.7864326933928735569664626221296E-17,
        +0.2505687311439975672324470645019E-17,
        -0.8047420364163909524537958682241E-18,
        +0.2604097118952053964443401104392E-18,
        -0.8486954164056412259482488834184E-19,
        +0.2784706882142337843359429186027E-19,
        -0.9195858953498612913687224151354E-20,
        +0.3055304318374238742247668225583E-20,
        -0.1021035455479477875902177048439E-20,
        +0.3431118190743757844000555680836E-21,
        -0.1159129341797749513376922463109E-21,
        +0.3935772844200255610836268229154E-22,
        -0.1342880980296717611956718989038E-22,
        +0.4603287883520002741659190305314E-23,
        -0.1585043927004064227810772499387E-23,
        +0.5481275667729675908925523755008E-24,
        -0.1903349371855047259064017948945E-24,
        +0.6635682302374008716777612115968E-25,
        -0.2322311650026314307975200986453E-25,
        +0.8157640113429179313142743695359E-26,
        -0.2875824240632900490057489929557E-26,
        +0.1017329450942901435079714319018E-26,
        -0.3610879108742216446575703490559E-27,
        +0.1285788540363993421256640342698E-27,
        -0.4592901037378547425160693022719E-28,
        +0.1645597033820713725812102485333E-28,
        -0.5913421299843501842087920271360E-29,
        +0.2131057006604993303479369509546E-29,
        -0.7701158157787598216982761745066E-30,
        +0.2790533307968930417581783777280E-30,
        -0.1013807715111284006452241367039E-30,
        +0.3692580158719624093658286216533E-31
    };
    static ityp aip2cs[37] =
    {
        -0.174314496929375513390355844011E-02,
        -0.167893854325541671632190613480E-02,
        +0.359653403352166035885983858114E-04,
        -0.138081860273922835457399383100E-05,
        +0.741122807731505298848699095233E-07,
        -0.500238203900133013130422866325E-08,
        +0.400693917417184240675446866355E-09,
        -0.367331242795905044199318496207E-10,
        +0.376034439592373852439592002918E-11,
        -0.422321332718747538026564938968E-12,
        +0.513509454033657070919618754120E-13,
        -0.669095850390477595651681356676E-14,
        +0.926667545641290648239550724382E-15,
        -0.135514382416070576333397356591E-15,
        +0.208115496312830995299006549335E-16,
        -0.334116499159176856871277570256E-17,
        +0.558578584585924316868032946585E-18,
        -0.969219040152365247518658209109E-19,
        +0.174045700128893206465696557738E-19,
        -0.322640979731130400247846333098E-20,
        +0.616074471106625258533259618986E-21,
        -0.120936347982490059076420676266E-21,
        +0.243632763310138108261570095786E-22,
        -0.502914221497457468943403144533E-23,
        +0.106224175543635689495470626133E-23,
        -0.229284284895989241509856324266E-24,
        +0.505181733929503744986884778666E-25,
        -0.113498123714412404979793920000E-25,
        +0.259765565985606980698374144000E-26,
        -0.605124621542939506172231679999E-27,
        +0.143359777966772800720295253333E-27,
        -0.345147757060899986280721066666E-28,
        +0.843875190213646740427025066666E-29,
        -0.209396142298188169434453333333E-29,
        +0.527008873478945503182848000000E-30,
        -0.134457433014553385789030399999E-30,
        +0.347570964526601147340117333333E-31
    };
    ityp eta;
    static dim_typ naif = 0;
    static dim_typ naig =  0;
    static dim_typ naip1 = 0;
    static dim_typ naip2 = 0;
    ityp sqrtx;
    ityp theta;
    ityp value;
    static ityp x32sml = 0.00;
    static ityp x3sml = 0.00;
    static ityp xbig = 0.00;
    ityp xm;
    ityp z;

    if ( naif == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        naif = r8_inits ( aifcs, 13, eta );
        naig = r8_inits ( aigcs, 13, eta );
        naip1 = r8_inits ( aip1cs, 57, eta );
        naip2 = r8_inits ( aip2cs, 37, eta );
        x3sml = pow ( eta, 0.3333 );
        x32sml = 1.3104 * x3sml * x3sml;
        xbig = pow ( r8_mach ( 2 ), 0.6666 );
    }

    if ( x < - 1.00 )
    {
        r8_aimp ( x, &xm, &theta );
        value = xm * cos ( theta );
    }
    else if ( 0.00 <= x && x <= x32sml )
    {
        z = 0.00;
        value = 0.3750 + ( r8_csevl ( z, aifcs, naif )- x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
    }
    else if ( abs ( x ) <= x3sml )
    {
        z = 0.00;
        value = 0.3750 + ( r8_csevl ( z, aifcs, naif )- x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
        value *= exp ( 2.0 * x * sqrt ( x ) / 3.0 );
    }
    else if ( x <= 1.00 )
    {
        z = x * x * x;
        value = 0.3750 + ( r8_csevl ( z, aifcs, naif )- x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
        value *= exp ( 2.0 * x * sqrt ( x ) / 3.0 );
    }
    else if ( x <= 4.00 )
    {
        sqrtx = sqrt ( x );
        z = ( 16.0 / ( x * sqrtx ) - 9.0 ) / 7.0;
        value = ( 0.28125 + r8_csevl ( z, aip1cs, naip1 ) ) / sqrt ( sqrtx );
    }
    else if ( x < xbig )
    {
        sqrtx = sqrt ( x );
        z = 16.00 / ( x * sqrtx ) - 1.00;
        value = ( 0.28125 + r8_csevl ( z, aip2cs, naip2 ) ) / sqrt ( sqrtx );
    }
    else
    {
        sqrtx = sqrt ( x );
        z = - 1.00;
        value = ( 0.28125 + r8_csevl ( z, aip2cs, naip2 ) ) / sqrt ( sqrtx );
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_aimp ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_AIMP evaluates the modulus and phase of the Airy function.
  Description:
    This function must only be called when X <= -1.0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double *AMPL, *PHI, the modulus and phase of the
    Airy function.
*/
{
	const it2pit * const s_data = data;
	const register ityp x = s_data->a0;
	ityp * ampl = s_data->a1;
	ityp * theta = s_data->a2;
	
    static ityp am20cs[57] =
    {
        +0.108716749086561856615730588125E-01,
        +0.369489228982663555091728665146E-03,
        +0.440680100484689563667507001327E-05,
        +0.143686762361911153929183952833E-06,
        +0.824275552390078308670628855353E-08,
        +0.684426758893661606173927278180E-09,
        +0.739566697282739287731004740213E-10,
        +0.974595633696825017638702600847E-11,
        +0.150076885829405775650973119497E-11,
        +0.262147910221527634206252854802E-12,
        +0.508354111376487180357278966914E-13,
        +0.107684753358811440492985997070E-13,
        +0.246091286618433429335914062617E-14,
        +0.600786380358656418436110373550E-15,
        +0.155449156102388071150651388384E-15,
        +0.423535125035576604426382780182E-16,
        +0.120862166289299840154401109189E-16,
        +0.359609651214658240861499706423E-17,
        +0.111134218386395638261774604677E-17,
        +0.355559532432366609893680289225E-18,
        +0.117433021600139309998766947387E-18,
        +0.399397454661077561389162200966E-19,
        +0.139576671528916310425606325640E-19,
        +0.500240055309236041393459280716E-20,
        +0.183552760958132679184834866457E-20,
        +0.688490998179202743197790112404E-21,
        +0.263631035611417012359996885105E-21,
        +0.102924890237338360287153563785E-21,
        +0.409246966671594885489762960571E-22,
        +0.165558573406734651039727903828E-22,
        +0.680797467063033356116599685727E-23,
        +0.284326559934079832419751134476E-23,
        +0.120507398348965255097287818819E-23,
        +0.517961243287505217976613610424E-24,
        +0.225622613427562816303268640887E-24,
        +0.995418801147745168832117078246E-25,
        +0.444551696397342424308280582053E-25,
        +0.200865195461501101425916097338E-25,
        +0.917786344151775165973885645402E-26,
        +0.423872958105589240661672197948E-26,
        +0.197789272007846092370846251490E-26,
        +0.932116351284620665680435253373E-27,
        +0.443482133249918099955611379722E-27,
        +0.212945672365573895594589552837E-27,
        +0.103158569651075977552209344907E-27,
        +0.504023773022591199157904590029E-28,
        +0.248301304570155945304046541005E-28,
        +0.123301783128562196054198238560E-28,
        +0.617033449920521746121976730507E-29,
        +0.311092617415918897233869792213E-29,
        +0.157983085201706173015269071503E-29,
        +0.807931987538283607678121339092E-30,
        +0.415997394138667562722951360052E-30,
        +0.215610934097716900471935862504E-30,
        +0.112468857265869178296752823613E-30,
        +0.590331560632838091123040811797E-31,
        +0.311735667692928562046280505333E-31 };
        static ityp am21cs[60] = {
        +0.592790266721309588375717482814E-02,
        +0.200569405393165186428695217690E-02,
        +0.911081850262275893553072526291E-04,
        +0.849894306372047155633172107475E-05,
        +0.113297908976913076637929215494E-05,
        +0.187517946100666496180950627804E-06,
        +0.359306519018245832699035211192E-07,
        +0.765757714071683864039093517470E-08,
        +0.176999967168039173925953460744E-08,
        +0.436259555654598932720546585535E-09,
        +0.113291641337853230035520085219E-09,
        +0.307257690982419244137868398126E-10,
        +0.864482416482201075541200465766E-11,
        +0.251015250060924402115104562212E-11,
        +0.749102496764440371601802227751E-12,
        +0.228996928487994073089565214432E-12,
        +0.715113658927987694949327491175E-13,
        +0.227607924959566841946395165061E-13,
        +0.736942142760886513969953227782E-14,
        +0.242328675267827490463991742006E-14,
        +0.808153774548239869283406558403E-15,
        +0.273008079804356086659174563386E-15,
        +0.933236070891385318473519474326E-16,
        +0.322508099681084622213867546973E-16,
        +0.112581932346444541217757573416E-16,
        +0.396699463986938821660259459530E-17,
        +0.141006567944319504660865034527E-17,
        +0.505302086537851213375537393032E-18,
        +0.182461523215945141197999102789E-18,
        +0.663584568262130466928029121642E-19,
        +0.242963731631276179741747455826E-19,
        +0.895238915123687802013669922963E-20,
        +0.331845289350050791260229250755E-20,
        +0.123706196188658315384437905922E-20,
        +0.463636677012390840306767734243E-21,
        +0.174653135947764475469758765989E-21,
        +0.661116810234991176307910643111E-22,
        +0.251409918994072486176125666459E-22,
        +0.960274995571732568694034386998E-23,
        +0.368324952289296395686436898078E-23,
        +0.141843138269159136145535939553E-23,
        +0.548342674276935830106345800990E-24,
        +0.212761054623118806650372562616E-24,
        +0.828443700849418591487734760953E-25,
        +0.323670563926127001421028600927E-25,
        +0.126868882963286057355055062493E-25,
        +0.498843818992121626935068934362E-26,
        +0.196734584467649390967119381790E-26,
        +0.778135971020326957713212064836E-27,
        +0.308633941498911152919192968451E-27,
        +0.122744647045453119789338037234E-27,
        +0.489431279134292205885241216204E-28,
        +0.195646879802909821175925099724E-28,
        +0.783988952922426171166311492266E-29,
        +0.314896914002484223748298978099E-29,
        +0.126769763137250681307067842559E-29,
        +0.511470691906900141641632107724E-30,
        +0.206801709795538770250900316706E-30,
        +0.837891344768519001325996867583E-31,
        +0.340168991971489802052339079577E-31
    };
    static ityp am22cs[74] =
    {
        -0.156284448062534112753545828583E-01,
        +0.778336445239681307018943100334E-02,
        +0.867057770477189528406072812110E-03,
        +0.156966273156113719469953482266E-03,
        +0.356396257143286511324100666302E-04,
        +0.924598335425043154495080090994E-05,
        +0.262110161850422389523194982066E-05,
        +0.791882216516012561489469982263E-06,
        +0.251041527921011847803162690862E-06,
        +0.826522320665407734472997712940E-07,
        +0.280571166281305264396384290014E-07,
        +0.976821090484680786674631273890E-08,
        +0.347407923227710343287279035573E-08,
        +0.125828132169836914219092738164E-08,
        +0.462988260641895264497330784625E-09,
        +0.172728258813604072468143128696E-09,
        +0.652319200131154135148574124970E-10,
        +0.249047168520982056019881087112E-10,
        +0.960156820553765948078189890126E-11,
        +0.373448002067726856974776596757E-11,
        +0.146417565032053391722216189678E-11,
        +0.578265471168512825475827881553E-12,
        +0.229915407244706118560254184494E-12,
        +0.919780711231997257150883662365E-13,
        +0.370060068813090065807504045556E-13,
        +0.149675761698672987823326345205E-13,
        +0.608361194938461148720451399443E-14,
        +0.248404087115121397635425326873E-14,
        +0.101862476526769080727914465339E-14,
        +0.419383856352753989429640310957E-15,
        +0.173318901762930756149702493501E-15,
        +0.718821902388508517820445406811E-16,
        +0.299123633598403607712470896113E-16,
        +0.124868990433238627855713110880E-16,
        +0.522829344609483661928651193632E-17,
        +0.219532961724713396595998454359E-17,
        +0.924298325229777281154410024332E-18,
        +0.390157708236091407825543197309E-18,
        +0.165093892693863707213759030367E-18,
        +0.700221815715994367565716554487E-19,
        +0.297651833616786915573214963506E-19,
        +0.126796539086902072571134261229E-19,
        +0.541243400697077628687581725061E-20,
        +0.231487350218155252296382133283E-20,
        +0.991920288386566563462623851167E-21,
        +0.425803015323732357158897608174E-21,
        +0.183101842973024501678402003088E-21,
        +0.788678712311075375564526811022E-22,
        +0.340254607386229874956582997235E-22,
        +0.147020881405712530791860892535E-22,
        +0.636211018324916957733348071767E-23,
        +0.275707050680980721919395987768E-23,
        +0.119645858090104071356261780457E-23,
        +0.519912545729242147981768210567E-24,
        +0.226217674847104475260575286850E-24,
        +0.985526113754431819448565068283E-25,
        +0.429870630332508717223681286187E-25,
        +0.187723641661580639829657670189E-25,
        +0.820721941772842137268801052115E-26,
        +0.359214665604615507812767944463E-26,
        +0.157390594612773315611458940587E-26,
        +0.690329781039333834965319153586E-27,
        +0.303092079078968534607859331415E-27,
        +0.133204934160481219185689121944E-27,
        +0.585978836851523490117937981442E-28,
        +0.258016868489487806338425080457E-28,
        +0.113712433637283667223632182863E-28,
        +0.501592557226068509236430548549E-29,
        +0.221445829395509373322569708484E-29,
        +0.978470283886507289984691416411E-30,
        +0.432695414934180170112000952983E-30,
        +0.191497288193994570612929860440E-30,
        +0.848164622402392354171298331562E-31,
        +0.375947065173955919947455052934E-31
    };
    static ityp ath0cs[53] =
    {
        -0.8172601764161634499840208700543E-01,
        -0.8004012824788273287596481113068E-03,
        -0.3186525268782113203795553628242E-05,
        -0.6688388266477509330741698865033E-07,
        -0.2931759284994564516506822463184E-08,
        -0.2011263760883621669049030307186E-09,
        -0.1877522678055973426074008166652E-10,
        -0.2199637137704601251899002199848E-11,
        -0.3071616682592272449025746605586E-12,
        -0.4936140553673418361025600985389E-13,
        -0.8902833722583660416935236969866E-14,
        -0.1768987764615272613656814199467E-14,
        -0.3817868689032277014678199609600E-15,
        -0.8851159014819947594156286509984E-16,
        -0.2184818181414365953149677679568E-16,
        -0.5700849046986452380599442295119E-17,
        -0.1563121122177875392516031795495E-17,
        -0.4481437996768995067906688776353E-18,
        -0.1337794883736188022044566044098E-18,
        -0.4143340036874114453776852445442E-19,
        -0.1327263385718805025080481164652E-19,
        -0.4385728589128440522215756835955E-20,
        -0.1491360695952818067686201743956E-20,
        -0.5208104738630711377154238188773E-21,
        -0.1864382222390498923872526604979E-21,
        -0.6830263751167969012975435381881E-22,
        -0.2557117058029329629296207591347E-22,
        -0.9770158640254300218246907254046E-23,
        -0.3805161433416679084068428254886E-23,
        -0.1509022750737054063493926482995E-23,
        -0.6087551341242424929005568014525E-24,
        -0.2495879513809711495425982124058E-24,
        -0.1039157654581920948909588084274E-24,
        -0.4390235913976846536974594969051E-25,
        -0.1880790678447990211675826820582E-25,
        -0.8165070764199462948863022205753E-26,
        -0.3589944503749750514266435585041E-26,
        -0.1597658126632132872981291608708E-26,
        -0.7193250175703823969113802835305E-27,
        -0.3274943012727856506209351132721E-27,
        -0.1507042445783690665816975047272E-27,
        -0.7006624198319904717843967949140E-28,
        -0.3289907402983718226528815678356E-28,
        -0.1559518084365146526445322711496E-28,
        -0.7460690508208254582833851119721E-29,
        -0.3600877034824662020563277249431E-29,
        -0.1752851437473772257350402219197E-29,
        -0.8603275775188512909623778628724E-30,
        -0.4256432603226946534668039480105E-30,
        -0.2122161865044262927723650698206E-30,
        -0.1065996156704879052472060798561E-30,
        -0.5393568608816949116410688086892E-31,
        -0.2748174851043954822278496517870E-31
    };
    static ityp ath1cs[58] =
    {
        -0.6972849916208883845888148415037E-01,
        -0.5108722790650044987073448077961E-02,
        -0.8644335996989755094525334749512E-04,
        -0.5604720044235263542188698916125E-05,
        -0.6045735125623897409156376640077E-06,
        -0.8639802632488334393219721138499E-07,
        -0.1480809484309927157147782480780E-07,
        -0.2885809334577236039999449908712E-08,
        -0.6191631975665699609309191231800E-09,
        -0.1431992808860957830931365259879E-09,
        -0.3518141102137214721504616874321E-10,
        -0.9084761919955078290070339808051E-11,
        -0.2446171672688598449343283664767E-11,
        -0.6826083203213446240828996710264E-12,
        -0.1964579931194940171278546257802E-12,
        -0.5808933227139693164009191265856E-13,
        -0.1759042249527441992795400959024E-13,
        -0.5440902932714896613632538945319E-14,
        -0.1715247407486806802622358519451E-14,
        -0.5500929233576991546871101847161E-15,
        -0.1791878287739317259495152638754E-15,
        -0.5920372520086694197778411062231E-16,
        -0.1981713027876483962470972206590E-16,
        -0.6713232347016352262049984343790E-17,
        -0.2299450243658281116122358619832E-17,
        -0.7957300928236376595304637145634E-18,
        -0.2779994027291784157172290233739E-18,
        -0.9798924361326985224406795480814E-19,
        -0.3482717006061574386702645565849E-19,
        -0.1247489122558599057173300058084E-19,
        -0.4501210041478228113487751824452E-20,
        -0.1635346244013352135596114164667E-20,
        -0.5980102897780336268098762265941E-21,
        -0.2200246286286123454028196295475E-21,
        -0.8142463073515085897408205291519E-22,
        -0.3029924773660042537432330709674E-22,
        -0.1133390098574623537722943969689E-22,
        -0.4260766024749295719283049889791E-23,
        -0.1609363396278189718797500634453E-23,
        -0.6106377190825026293045330444287E-24,
        -0.2326954318021694061836577887573E-24,
        -0.8903987877472252604474129558186E-25,
        -0.3420558530005675024117914752341E-25,
        -0.1319026715257272659017212100607E-25,
        -0.5104899493612043091316191177386E-26,
        -0.1982599478474547451242444663466E-26,
        -0.7725702356880830535636111851519E-27,
        -0.3020234733664680100815776863573E-27,
        -0.1184379739074169993712946380800E-27,
        -0.4658430227922308520573252840106E-28,
        -0.1837554188100384647157502006613E-28,
        -0.7268566894427990953321876684800E-29,
        -0.2882863120391468135527089875626E-29,
        -0.1146374629459906350417591664639E-29,
        -0.4570031437748533058179991688533E-30,
        -0.1826276602045346104809934028799E-30,
        -0.7315349993385250469111066350933E-31,
        -0.2936925599971429781637815773866E-31
    };
    static ityp ath2cs[72] =
    {
        +0.4405273458718778997061127057775E-02,
        -0.3042919452318454608483844239873E-01,
        -0.1385653283771793791602692842653E-02,
        -0.1804443908954952302670486910952E-03,
        -0.3380847108327308671057465323618E-04,
        -0.7678183535229023055257676817765E-05,
        -0.1967839443716035324690935417077E-05,
        -0.5483727115877700361586143659281E-06,
        -0.1625461550532612452712696212258E-06,
        -0.5053049981268895015277637842078E-07,
        -0.1631580701124066881183851715617E-07,
        -0.5434204112348517507963436694817E-08,
        -0.1857398556409900325763850109630E-08,
        -0.6489512033326108816213513640676E-09,
        -0.2310594885800944720482995987079E-09,
        -0.8363282183204411682819329546745E-10,
        -0.3071196844890191462660661303891E-10,
        -0.1142367142432716819409514579892E-10,
        -0.4298116066345803065822470108971E-11,
        -0.1633898699596715440601646086632E-11,
        -0.6269328620016619432123443754076E-12,
        -0.2426052694816257357356159203991E-12,
        -0.9461198321624039090742527765052E-13,
        -0.3716060313411504806847798281269E-13,
        -0.1469155684097526763170138810309E-13,
        -0.5843694726140911944556401363094E-14,
        -0.2337502595591951298832675034934E-14,
        -0.9399231371171435401160167358411E-15,
        -0.3798014669372894500076335263715E-15,
        -0.1541731043984972524883443681775E-15,
        -0.6285287079535307162925662365202E-16,
        -0.2572731812811455424755383992774E-16,
        -0.1057098119354017809340974866555E-16,
        -0.4359080267402696966695992699964E-17,
        -0.1803634315959978013953176945540E-17,
        -0.7486838064380536821719431676914E-18,
        -0.3117261367347604656799597209985E-18,
        -0.1301687980927700734792871620696E-18,
        -0.5450527587519522468973883909909E-19,
        -0.2288293490114231872268635931903E-19,
        -0.9631059503829538655655060440088E-20,
        -0.4063281001524614089092195416434E-20,
        -0.1718203980908026763900413858510E-20,
        -0.7281574619892536367415322473328E-21,
        -0.3092352652680643127960680345790E-21,
        -0.1315917855965440490383417023254E-21,
        -0.5610606786087055512664907412668E-22,
        -0.2396621894086355206020304337895E-22,
        -0.1025574332390581200832954423924E-22,
        -0.4396264138143656476403607323663E-23,
        -0.1887652998372577373342508719450E-23,
        -0.8118140359576807603579433230445E-24,
        -0.3496734274366286856375952089214E-24,
        -0.1508402925156873215171751475867E-24,
        -0.6516268284778671059787773834341E-25,
        -0.2818945797529207424505942114583E-25,
        -0.1221127596512262744598094464505E-25,
        -0.5296674341169867168620011705073E-26,
        -0.2300359270773673431358870971744E-26,
        -0.1000279482355367494781220348930E-26,
        -0.4354760404180879394806893162179E-27,
        -0.1898056134741477522515482827030E-27,
        -0.8282111868712974697554009309315E-28,
        -0.3617815493066569006586213484374E-28,
        -0.1582018896178003654858941843636E-28,
        -0.6925068597802270011772820383247E-29,
        -0.3034390239778629128908629727335E-29,
        -0.1330889568166725224761977446509E-29,
        -0.5842848522173090120487606971706E-30,
        -0.2567488423238302631121274357678E-30,
        -0.1129232322268882185791505819151E-30,
        -0.4970947029753336916550570105023E-31
    };
    ityp eta;
    static dim_typ nam20 = 0;
    static dim_typ nam21 = 0;
    static dim_typ nam22 = 0;
    static dim_typ nath0 = 0;
    static dim_typ nath1 = 0;
    static dim_typ nath2 = 0;
    static ityp pi4 = 0.78539816339744830961566084581988;
    ityp sqrtx;
    static ityp xsml = 0.00;
    ityp z;

    if ( nam20 == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        nam20 = r8_inits ( am20cs, 57, eta );
        nath0 = r8_inits ( ath0cs, 53, eta );
        nam21 = r8_inits ( am21cs, 60, eta );
        nath1 = r8_inits ( ath1cs, 58, eta );
        nam22 = r8_inits ( am22cs, 74, eta );
        nath2 = r8_inits ( ath2cs, 72, eta );
        xsml = - pow ( 128.0 / r8_mach ( 3 ), 0.3333 );
    }

    if ( x <= xsml )
    {
        z = 1.00;
        *ampl = 0.3125 + r8_csevl ( z, am20cs, nam20 );
        *theta = - 0.625 + r8_csevl ( z, ath0cs, nath0 );
    }
    else if ( x < - 4.00 )
    {
        z = 128.00 / x / x / x + 1.0;
        *ampl = 0.3125 + r8_csevl ( z, am20cs, nam20 );
        *theta = - 0.625 + r8_csevl ( z, ath0cs, nath0 );
    }
    else if ( x < - 2.00 )
    {
        z = ( 128.00 / x / x / x + 9.00 ) / 7.00;
        *ampl = 0.3125 + r8_csevl ( z, am21cs, nam21 );
        *theta = - 0.625 + r8_csevl ( z, ath1cs, nath1 );
    }
    else if ( x <= - 1.0 )
    {
        z = ( 16.00 / x / x / x + 9.00 ) / 7.00;
        *ampl = 0.3125 + r8_csevl ( z, am22cs, nam22 );
        *theta = - 0.625 + r8_csevl ( z, ath2cs, nath2 );
    }
    else
        return NULL;

    sqrtx = sqrt ( - x );
    *ampl = sqrt ( *ampl / sqrtx );
    *theta = pi4 - x * sqrtx * *theta;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_aint ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_AINT truncates an R8 argument to an integer.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    John Burkardt.
  Parameters:
    Input, double X, the argument.
    Output, double R8_AINT, the truncated version of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	result = x<0.0E+00 ? - ( ityp ) ( ( dim_typ ) ( abs ( x ) ) ) : ( ityp ) ( ( dim_typ ) ( abs ( x ) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_asinh ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_ASINH evaluates the arc-sine of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_ASINH, the arc-hyperbolic sine of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp aln2 = 0.69314718055994530941723212145818;
    static ityp asnhcs[39] =
    {
        -0.12820039911738186343372127359268E+00,
        -0.58811761189951767565211757138362E-01,
        +0.47274654322124815640725249756029E-02,
        -0.49383631626536172101360174790273E-03,
        +0.58506207058557412287494835259321E-04,
        -0.74669983289313681354755069217188E-05,
        +0.10011693583558199265966192015812E-05,
        -0.13903543858708333608616472258886E-06,
        +0.19823169483172793547317360237148E-07,
        -0.28847468417848843612747272800317E-08,
        +0.42672965467159937953457514995907E-09,
        -0.63976084654366357868752632309681E-10,
        +0.96991686089064704147878293131179E-11,
        -0.14844276972043770830246658365696E-11,
        +0.22903737939027447988040184378983E-12,
        -0.35588395132732645159978942651310E-13,
        +0.55639694080056789953374539088554E-14,
        -0.87462509599624678045666593520162E-15,
        +0.13815248844526692155868802298129E-15,
        -0.21916688282900363984955142264149E-16,
        +0.34904658524827565638313923706880E-17,
        -0.55785788400895742439630157032106E-18,
        +0.89445146617134012551050882798933E-19,
        -0.14383426346571317305551845239466E-19,
        +0.23191811872169963036326144682666E-20,
        -0.37487007953314343674570604543999E-21,
        +0.60732109822064279404549242880000E-22,
        -0.98599402764633583177370173440000E-23,
        +0.16039217452788496315232638293333E-23,
        -0.26138847350287686596716134399999E-24,
        +0.42670849606857390833358165333333E-25,
        -0.69770217039185243299730773333333E-26,
        +0.11425088336806858659812693333333E-26,
        -0.18735292078860968933021013333333E-27,
        +0.30763584414464922794065920000000E-28,
        -0.50577364031639824787046399999999E-29,
        +0.83250754712689142224213333333333E-30,
        -0.13718457282501044163925333333333E-30,
        +0.22629868426552784104106666666666E-31
    };
    static dim_typ nterms = 0;
    static ityp sqeps = 0.00;
    ityp value;
    static ityp xmax = 0.00;
    ityp y;

    if ( nterms == 0 )
    {
        nterms = r8_inits ( asnhcs, 39, 0.1 * r8_mach ( 3 ) );
        sqeps = sqrt ( r8_mach ( 3 ) );
        xmax = 1.00 / sqeps;
    }

    y = abs ( x );

    if ( y <= sqeps )
        value = x;
    else if ( y <= 1.00 )
        value = x * ( 1.00 + r8_csevl ( 2.00 * x * x - 1.00, asnhcs, nterms ) );
    else if ( y < xmax )
    {
        value = log ( y + sqrt ( y * y + 1.00 ) );
        if ( x < 0.00 )
        value *= -1;
    }
    else
    {
        value = aln2 + log ( y );
        if(x < 0.00)
            value *= -1;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_b1mp ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_B1MP evaluates the modulus and phase for the Bessel J1 and Y1 functions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double *AMPL, *THETA, the modulus and phase.
*/
{
	const it2pit * const s_data = data;
	const register ityp x = s_data->a0;
	ityp * ampl = s_data->a1;
	ityp * theta = s_data->a2;
	
    static ityp bm12cs[40] =
    {
        +0.9807979156233050027272093546937E-01,
        +0.1150961189504685306175483484602E-02,
        -0.4312482164338205409889358097732E-05,
        +0.5951839610088816307813029801832E-07,
        -0.1704844019826909857400701586478E-08,
        +0.7798265413611109508658173827401E-10,
        -0.4958986126766415809491754951865E-11,
        +0.4038432416421141516838202265144E-12,
        -0.3993046163725175445765483846645E-13,
        +0.4619886183118966494313342432775E-14,
        -0.6089208019095383301345472619333E-15,
        +0.8960930916433876482157048041249E-16,
        -0.1449629423942023122916518918925E-16,
        +0.2546463158537776056165149648068E-17,
        -0.4809472874647836444259263718620E-18,
        +0.9687684668292599049087275839124E-19,
        -0.2067213372277966023245038117551E-19,
        +0.4646651559150384731802767809590E-20,
        -0.1094966128848334138241351328339E-20,
        +0.2693892797288682860905707612785E-21,
        -0.6894992910930374477818970026857E-22,
        +0.1830268262752062909890668554740E-22,
        -0.5025064246351916428156113553224E-23,
        +0.1423545194454806039631693634194E-23,
        -0.4152191203616450388068886769801E-24,
        +0.1244609201503979325882330076547E-24,
        -0.3827336370569304299431918661286E-25,
        +0.1205591357815617535374723981835E-25,
        -0.3884536246376488076431859361124E-26,
        +0.1278689528720409721904895283461E-26,
        -0.4295146689447946272061936915912E-27,
        +0.1470689117829070886456802707983E-27,
        -0.5128315665106073128180374017796E-28,
        +0.1819509585471169385481437373286E-28,
        -0.6563031314841980867618635050373E-29,
        +0.2404898976919960653198914875834E-29,
        -0.8945966744690612473234958242979E-30,
        +0.3376085160657231026637148978240E-30,
        -0.1291791454620656360913099916966E-30,
        +0.5008634462958810520684951501254E-31
    };
    static ityp bm1cs[37] =
    {
        +0.1069845452618063014969985308538,
        +0.3274915039715964900729055143445E-02,
        -0.2987783266831698592030445777938E-04,
        +0.8331237177991974531393222669023E-06,
        -0.4112665690302007304896381725498E-07,
        +0.2855344228789215220719757663161E-08,
        -0.2485408305415623878060026596055E-09,
        +0.2543393338072582442742484397174E-10,
        -0.2941045772822967523489750827909E-11,
        +0.3743392025493903309265056153626E-12,
        -0.5149118293821167218720548243527E-13,
        +0.7552535949865143908034040764199E-14,
        -0.1169409706828846444166290622464E-14,
        +0.1896562449434791571721824605060E-15,
        -0.3201955368693286420664775316394E-16,
        +0.5599548399316204114484169905493E-17,
        -0.1010215894730432443119390444544E-17,
        +0.1873844985727562983302042719573E-18,
        -0.3563537470328580219274301439999E-19,
        +0.6931283819971238330422763519999E-20,
        -0.1376059453406500152251408930133E-20,
        +0.2783430784107080220599779327999E-21,
        -0.5727595364320561689348669439999E-22,
        +0.1197361445918892672535756799999E-22,
        -0.2539928509891871976641440426666E-23,
        +0.5461378289657295973069619199999E-24,
        -0.1189211341773320288986289493333E-24,
        +0.2620150977340081594957824000000E-25,
        -0.5836810774255685901920938666666E-26,
        +0.1313743500080595773423615999999E-26,
        -0.2985814622510380355332778666666E-27,
        +0.6848390471334604937625599999999E-28,
        -0.1584401568222476721192960000000E-28,
        +0.3695641006570938054301013333333E-29,
        -0.8687115921144668243012266666666E-30,
        +0.2057080846158763462929066666666E-30,
        -0.4905225761116225518523733333333E-31
    };
    static ityp bt12cs[39] =
    {
        +0.73823860128742974662620839792764,
        -0.33361113174483906384470147681189E-02,
        +0.61463454888046964698514899420186E-04,
        -0.24024585161602374264977635469568E-05,
        +0.14663555577509746153210591997204E-06,
        -0.11841917305589180567005147504983E-07,
        +0.11574198963919197052125466303055E-08,
        -0.13001161129439187449366007794571E-09,
        +0.16245391141361731937742166273667E-10,
        -0.22089636821403188752155441770128E-11,
        +0.32180304258553177090474358653778E-12,
        -0.49653147932768480785552021135381E-13,
        +0.80438900432847825985558882639317E-14,
        -0.13589121310161291384694712682282E-14,
        +0.23810504397147214869676529605973E-15,
        -0.43081466363849106724471241420799E-16,
        +0.80202544032771002434993512550400E-17,
        -0.15316310642462311864230027468799E-17,
        +0.29928606352715568924073040554666E-18,
        -0.59709964658085443393815636650666E-19,
        +0.12140289669415185024160852650666E-19,
        -0.25115114696612948901006977706666E-20,
        +0.52790567170328744850738380799999E-21,
        -0.11260509227550498324361161386666E-21,
        +0.24348277359576326659663462400000E-22,
        -0.53317261236931800130038442666666E-23,
        +0.11813615059707121039205990399999E-23,
        -0.26465368283353523514856789333333E-24,
        +0.59903394041361503945577813333333E-25,
        -0.13690854630829503109136383999999E-25,
        +0.31576790154380228326413653333333E-26,
        -0.73457915082084356491400533333333E-27,
        +0.17228081480722747930705920000000E-27,
        -0.40716907961286507941068800000000E-28,
        +0.96934745136779622700373333333333E-29,
        -0.23237636337765716765354666666666E-29,
        +0.56074510673522029406890666666666E-30,
        -0.13616465391539005860522666666666E-30,
        +0.33263109233894654388906666666666E-31
    };
    static ityp bth1cs[44] =
    {
        +0.74749957203587276055443483969695,
        -0.12400777144651711252545777541384E-02,
        +0.99252442404424527376641497689592E-05,
        -0.20303690737159711052419375375608E-06,
        +0.75359617705690885712184017583629E-08,
        -0.41661612715343550107630023856228E-09,
        +0.30701618070834890481245102091216E-10,
        -0.28178499637605213992324008883924E-11,
        +0.30790696739040295476028146821647E-12,
        -0.38803300262803434112787347554781E-13,
        +0.55096039608630904934561726208562E-14,
        -0.86590060768383779940103398953994E-15,
        +0.14856049141536749003423689060683E-15,
        -0.27519529815904085805371212125009E-16,
        +0.54550796090481089625036223640923E-17,
        -0.11486534501983642749543631027177E-17,
        +0.25535213377973900223199052533522E-18,
        -0.59621490197413450395768287907849E-19,
        +0.14556622902372718620288302005833E-19,
        -0.37022185422450538201579776019593E-20,
        +0.97763074125345357664168434517924E-21,
        -0.26726821639668488468723775393052E-21,
        +0.75453300384983271794038190655764E-22,
        -0.21947899919802744897892383371647E-22,
        +0.65648394623955262178906999817493E-23,
        -0.20155604298370207570784076869519E-23,
        +0.63417768556776143492144667185670E-24,
        -0.20419277885337895634813769955591E-24,
        +0.67191464220720567486658980018551E-25,
        -0.22569079110207573595709003687336E-25,
        +0.77297719892989706370926959871929E-26,
        -0.26967444512294640913211424080920E-26,
        +0.95749344518502698072295521933627E-27,
        -0.34569168448890113000175680827627E-27,
        +0.12681234817398436504211986238374E-27,
        -0.47232536630722639860464993713445E-28,
        +0.17850008478186376177858619796417E-28,
        -0.68404361004510395406215223566746E-29,
        +0.26566028671720419358293422672212E-29,
        -0.10450402527914452917714161484670E-29,
        +0.41618290825377144306861917197064E-30,
        -0.16771639203643714856501347882887E-30,
        +0.68361997776664389173535928028528E-31,
        -0.28172247861233641166739574622810E-31
    };
    ityp eta;
    static dim_typ nbm1 = 0;
    static dim_typ nbm12 = 0;
    static dim_typ nbt12 = 0;
    static dim_typ nbth1 = 0;
    static ityp pi4 = 0.785398163397448309615660845819876;
    static ityp xmax = 0.0;
    ityp z;

    if ( nbm1 == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        nbm1 = r8_inits ( bm1cs, 37, eta );
        nbt12 = r8_inits ( bt12cs, 39, eta );
        nbm12 = r8_inits ( bm12cs, 40, eta );
        nbth1 = r8_inits ( bth1cs, 44, eta );
        xmax = 1.00 / r8_mach ( 4 );
    }

    if ( x < 4.00 )
        return NULL;
    else if ( x <= 8.0 )
    {
        z = ( 128.00 / x / x - 5.00 ) / 3.00;
        *ampl = ( 0.75 + r8_csevl ( z, bm1cs, nbm1 ) ) / sqrt ( x );
        *theta = x - 3.00 * pi4 + r8_csevl ( z, bt12cs, nbt12 ) / x;
    }
    else
    {
        z = 128.00 / x / x - 1.00;
        *ampl = ( 0.75 + r8_csevl ( z, bm12cs, nbm12 ) ) / sqrt ( x );
        *theta = x - 3.00 * pi4 + r8_csevl ( z, bth1cs, nbth1 ) / x;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_bi ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BI evaluates the Airy function Bi of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BI, the Airy function Bi of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp bifcs[13] =
    {
        -0.16730216471986649483537423928176E-01,
        +0.10252335834249445611426362777757,
        +0.17083092507381516539429650242013E-02,
        +0.11862545467744681179216459210040E-04,
        +0.44932907017792133694531887927242E-07,
        +0.10698207143387889067567767663628E-09,
        +0.17480643399771824706010517628573E-12,
        +0.20810231071761711025881891834399E-15,
        +0.18849814695665416509927971733333E-18,
        +0.13425779173097804625882666666666E-21,
        +0.77159593429658887893333333333333E-25,
        +0.36533879617478566399999999999999E-28,
        +0.14497565927953066666666666666666E-31
    };
    static ityp bif2cs[15] =
    {
        +0.0998457269381604104468284257993,
        +0.47862497786300553772211467318231,
        +0.25155211960433011771324415436675E-01,
        +0.58206938852326456396515697872216E-03,
        +0.74997659644377865943861457378217E-05,
        +0.61346028703493836681403010356474E-07,
        +0.34627538851480632900434268733359E-09,
        +0.14288910080270254287770846748931E-11,
        +0.44962704298334641895056472179200E-14,
        +0.11142323065833011708428300106666E-16,
        +0.22304791066175002081517866666666E-19,
        +0.36815778736393142842922666666666E-22,
        +0.50960868449338261333333333333333E-25,
        +0.60003386926288554666666666666666E-28,
        +0.60827497446570666666666666666666E-31
    };
    static ityp bigcs[13] =
    {
        +0.22466223248574522283468220139024E-01,
        +0.37364775453019545441727561666752E-01,
        +0.44476218957212285696215294326639E-03,
        +0.24708075636329384245494591948882E-05,
        +0.79191353395149635134862426285596E-08,
        +0.16498079851827779880887872402706E-10,
        +0.24119906664835455909247501122841E-13,
        +0.26103736236091436985184781269333E-16,
        +0.21753082977160323853123792000000E-19,
        +0.14386946400390433219483733333333E-22,
        +0.77349125612083468629333333333333E-26,
        +0.34469292033849002666666666666666E-29,
        +0.12938919273216000000000000000000E-32
    };
    static ityp big2cs[15] =
    {
        +0.033305662145514340465176188111647,
        +0.161309215123197067613287532084943,
        +0.631900730961342869121615634921173E-02,
        +0.118790456816251736389780192304567E-03,
        +0.130453458862002656147116485012843E-05,
        +0.937412599553521729546809615508936E-08,
        +0.474580188674725153788510169834595E-10,
        +0.178310726509481399800065667560946E-12,
        +0.516759192784958180374276356640000E-15,
        +0.119004508386827125129496251733333E-17,
        +0.222982880666403517277063466666666E-20,
        +0.346551923027689419722666666666666E-23,
        +0.453926336320504514133333333333333E-26,
        +0.507884996513522346666666666666666E-29,
        +0.491020674696533333333333333333333E-32
    };
    ityp eta;
    static dim_typ nbif = 0;
    static dim_typ nbif2 = 0;
    static dim_typ nbig = 0;
    static dim_typ nbig2 = 0;
    ityp theta;
    ityp value;
    static ityp x3sml = 0.00;
    ityp xm;
    static ityp xmax = 0.00;
    ityp z;

    if ( nbif == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        nbif = r8_inits ( bifcs, 13, eta );
        nbig = r8_inits ( bigcs, 13, eta );
        nbif2 = r8_inits ( bif2cs, 15, eta );
        nbig2 = r8_inits ( big2cs, 15, eta );
        x3sml = pow ( eta, 0.3333 );
        xmax = pow ( 1.50 * log ( r8_mach ( 2 ) ), 0.6666 );
    }

    if ( x < - 1.00 )
    {
        r8_aimp ( x, &xm, &theta );
        value = xm * sin ( theta );
    }
    else if ( abs ( x ) <= x3sml )
    {
        z = 0.00;
        value = 0.625 + r8_csevl ( z, bifcs, nbif )+ x * ( 0.4375 + r8_csevl ( z, bigcs, nbig ) );
    }
    else if ( x <= 1.00 )
    {
        z = x * x * x;
        value = 0.625 + r8_csevl ( z, bifcs, nbif )+ x * ( 0.4375 + r8_csevl ( z, bigcs, nbig ) );
    }
    else if ( x <= 2.00 )
    {
        z = ( 2.0 * x * x * x - 9.00 ) / 7.00;
        value = 1.125 + r8_csevl ( z, bif2cs, nbif2 )+ x * ( 0.625 + r8_csevl ( z, big2cs, nbig2 ) );
    }
    else
        value = r8_bie ( x ) * exp ( 2.0 * x * sqrt ( x ) / 3.00 );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_bid ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BID evaluates the derivative of the Airy function Bi of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BID, the derivative of the Airy function Bi of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp bif2cs[15] =
    {
        0.32349398760352203352119193596266015,
        0.08629787153556355913888835323811100,
        0.00299402555265539742613821050727155,
        0.00005143052836466163720464316950821,
        0.00000052584025003681146026033098613,
        0.00000000356175137395770028102730600,
        0.00000000001714686400714584830518308,
        0.00000000000006166351969232555406693,
        0.00000000000000017191082154315985806,
        0.00000000000000000038236889518803943,
        0.00000000000000000000069424173624884,
        0.00000000000000000000000104833932510,
        0.00000000000000000000000000133721972,
        0.00000000000000000000000000000145986,
        0.00000000000000000000000000000000138
    };
    static ityp bifcs[13] =
    {
        0.115353679082857024267474446284908879,
        0.020500789404919287530357789445940252,
        0.000213529027890287581892679619451158,
        0.000001078396061467683042209155523569,
        0.000000003209470883320666783353670420,
        0.000000000006293040671833540390213316,
        0.000000000000008740304300063083340121,
        0.000000000000000009047915683496049529,
        0.000000000000000000007249923164709251,
        0.000000000000000000000004629576649604,
        0.000000000000000000000000002411236436,
        0.000000000000000000000000000001043825,
        0.000000000000000000000000000000000382
    };
    static ityp big2cs[16] =
    {
        1.606299946362129457759284537862622883,
        0.744908881987608865201476685194753972,
        0.047013873861027737964095177635353019,
        0.001228442206254823907016188785848091,
        0.000017322241225662362670987355613727,
        0.000000152190165236801893711508366563,
        0.000000000911356024911957704145528786,
        0.000000000003954791842356644201722554,
        0.000000000000013001737033862320007309,
        0.000000000000000033493506858269079763,
        0.000000000000000000069419094403694057,
        0.000000000000000000000118248256604581,
        0.000000000000000000000000168462493472,
        0.000000000000000000000000000203684674,
        0.000000000000000000000000000000211619,
        0.000000000000000000000000000000000191
    };
    static ityp bigcs[13] =
    {
        -0.0971964404164435373897790974606802,
        0.1495035768431670665710843445326264,
        0.0031135253871213260419419176839631,
        0.0000247085705798212967777021920569,
        0.0000001029496277313786081987324295,
        0.0000000002639703739869432892676778,
        0.0000000000004582792707803206608181,
        0.0000000000000005742829740893447321,
        0.0000000000000000005438275385238549,
        0.0000000000000000000004028347267083,
        0.0000000000000000000000002397823826,
        0.0000000000000000000000000001171956,
        0.0000000000000000000000000000000479
    };
    ityp eta;
    static dim_typ nbif = 0;
    static dim_typ nbif2 = 0;
    static dim_typ nbig = 0;
    static dim_typ nbig2 = 0;
    ityp phi;
    ityp value;
    ityp x2;
    static ityp x2sml = 0.00;
    ityp x3;
    static ityp x3sml = 0.00;
    static ityp xmax = 0.00;
    ityp xn;
    ityp z;

    if ( nbif == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        nbif = r8_inits ( bifcs, 13, eta );
        nbig = r8_inits ( bigcs, 13, eta );
        nbif2 = r8_inits ( bif2cs, 15, eta );
        nbig2 = r8_inits ( big2cs, 16, eta );
        x2sml = sqrt ( eta );
        x3sml = pow ( eta, 0.3333 );
        xmax = pow ( 1.50 * log ( r8_mach ( 2 ) ), 0.6666 );
    }

    if ( x < - 1.00 )
    {
        r8_admp ( x, &xn, &phi );
        value = xn * sin ( phi );
    }
    else if ( abs ( x ) <= x2sml )
    {
        x2 = x3 = 0.00;
        value = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25 )+ r8_csevl ( x3, bigcs, nbig ) + 0.50;
    }
    else if ( abs ( x ) <= x3sml )
    {
        x2 = x * x;
        x3 = 0.00;
        value = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25 )+ r8_csevl ( x3, bigcs, nbig ) + 0.50;
    }
    else if ( x <= 1.00 )
    {
        x2 = x * x;
        x3 = x * x * x;
        value = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25 )+ r8_csevl ( x3, bigcs, nbig ) + 0.50;
    }
    else if ( x <= 2.00 )
    {
        z = ( 2.00 * x * x * x - 9.00 ) / 7.00;
        value = x * x * ( r8_csevl ( z, bif2cs, nbif2 ) + 0.25 ) +r8_csevl ( z, big2cs, nbig2 ) + 0.50;
    }
    else
        value = r8_bide ( x ) * exp ( 2.00 * x * sqrt ( x ) / 3.00 );
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_bide ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BIDE: exponentially scaled derivative, Airy function Bi of an R8 argument.
  Discussion:
    if X < 0,
      R8_BIDE ( X ) = R8_BID ( X )
    else
      R8_BIDE ( X ) = R8_BID ( X ) * exp ( - 2/3 * X^(3/2) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BIDE, the exponentially scaled derivative of
    the Airy function Bi of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp atr = 8.75069057084843450880771988210148;
    static ityp bif2cs[15] =
    {
        0.32349398760352203352119193596266015,
        0.08629787153556355913888835323811100,
        0.00299402555265539742613821050727155,
        0.00005143052836466163720464316950821,
        0.00000052584025003681146026033098613,
        0.00000000356175137395770028102730600,
        0.00000000001714686400714584830518308,
        0.00000000000006166351969232555406693,
        0.00000000000000017191082154315985806,
        0.00000000000000000038236889518803943,
        0.00000000000000000000069424173624884,
        0.00000000000000000000000104833932510,
        0.00000000000000000000000000133721972,
        0.00000000000000000000000000000145986,
        0.00000000000000000000000000000000138
    };
    static ityp bifcs[13] =
    {
        0.115353679082857024267474446284908879,
        0.020500789404919287530357789445940252,
        0.000213529027890287581892679619451158,
        0.000001078396061467683042209155523569,
        0.000000003209470883320666783353670420,
        0.000000000006293040671833540390213316,
        0.000000000000008740304300063083340121,
        0.000000000000000009047915683496049529,
        0.000000000000000000007249923164709251,
        0.000000000000000000000004629576649604,
        0.000000000000000000000000002411236436,
        0.000000000000000000000000000001043825,
        0.000000000000000000000000000000000382
    };
    static ityp big2cs[16] =
    {
        1.606299946362129457759284537862622883,
        0.744908881987608865201476685194753972,
        0.047013873861027737964095177635353019,
        0.001228442206254823907016188785848091,
        0.000017322241225662362670987355613727,
        0.000000152190165236801893711508366563,
        0.000000000911356024911957704145528786,
        0.000000000003954791842356644201722554,
        0.000000000000013001737033862320007309,
        0.000000000000000033493506858269079763,
        0.000000000000000000069419094403694057,
        0.000000000000000000000118248256604581,
        0.000000000000000000000000168462493472,
        0.000000000000000000000000000203684674,
        0.000000000000000000000000000000211619,
        0.000000000000000000000000000000000191
    };
    static ityp bigcs[13] =
    {
        -0.0971964404164435373897790974606802,
        0.1495035768431670665710843445326264,
        0.0031135253871213260419419176839631,
        0.0000247085705798212967777021920569,
        0.0000001029496277313786081987324295,
        0.0000000002639703739869432892676778,
        0.0000000000004582792707803206608181,
        0.0000000000000005742829740893447321,
        0.0000000000000000005438275385238549,
        0.0000000000000000000004028347267083,
        0.0000000000000000000000002397823826,
        0.0000000000000000000000000001171956,
        0.0000000000000000000000000000000479
    };
    static ityp bip1cs[47] =
    {
        -0.17291873510795537186124679823741003,
        -0.01493584929846943639486231021818675,
        -0.00054711049516785663990658697874460,
        0.00015379662929584083449573727856666,
        0.00001543534761921794131028948022869,
        -0.00000654341138519060129226087106765,
        0.00000037280824078787032232152275240,
        0.00000020720783881887480080810710514,
        -0.00000006581733364696191689495883922,
        0.00000000749267463539288212986048985,
        0.00000000111013368840707147698890101,
        -0.00000000072651405529159512323880794,
        0.00000000017827235598470153962165668,
        -0.00000000002173463524809506269656807,
        -0.00000000000203020349653882594017049,
        0.00000000000193118272294077519319859,
        -0.00000000000060449525048290296023117,
        0.00000000000012094496248933664277802,
        -0.00000000000001251088360074479784619,
        -0.00000000000000199173832424881344036,
        0.00000000000000151540816342864303038,
        -0.00000000000000049768927059816240250,
        0.00000000000000011545959731810501403,
        -0.00000000000000001863286862907983871,
        0.00000000000000000099330392344759104,
        0.00000000000000000068182083667412417,
        -0.00000000000000000034854456479650551,
        0.00000000000000000010860382134235961,
        -0.00000000000000000002599290185240166,
        0.00000000000000000000476895370459000,
        -0.00000000000000000000051946940777177,
        -0.00000000000000000000005925575044912,
        0.00000000000000000000005746008970972,
        -0.00000000000000000000002186119806494,
        0.00000000000000000000000624124294738,
        -0.00000000000000000000000146003421785,
        0.00000000000000000000000027493893904,
        -0.00000000000000000000000003474678018,
        -0.00000000000000000000000000109303694,
        0.00000000000000000000000000261972744,
        -0.00000000000000000000000000112365018,
        0.00000000000000000000000000035152059,
        -0.00000000000000000000000000009167601,
        0.00000000000000000000000000002040203,
        -0.00000000000000000000000000000373038,
        0.00000000000000000000000000000046070,
        0.00000000000000000000000000000001748
    };
    static ityp bip2cs[88] =
    {
        -0.13269705443526630494937031210217135,
        -0.00568443626045977481306046339037428,
        -0.00015643601119611609623698471216660,
        -0.00001136737203679562267336053207940,
        -0.00000143464350991283669643136951338,
        -0.00000018098531185164131868746481700,
        0.00000000926177343610865546229511422,
        0.00000001710005490720592181887296162,
        0.00000000476698163503781708252686849,
        -0.00000000035195022023163141945397159,
        -0.00000000058890614315886871574147635,
        -0.00000000006678499607795537597612089,
        0.00000000006395565101720391190697713,
        0.00000000001554529427064394106403245,
        -0.00000000000792396999744612971684001,
        -0.00000000000258326242689717798947525,
        0.00000000000121655047787849117978773,
        0.00000000000038707207172899985942258,
        -0.00000000000022487045479618229130656,
        -0.00000000000004953476515684046293493,
        0.00000000000004563781601526912756017,
        0.00000000000000332998314345014118494,
        -0.00000000000000921750185832874202719,
        0.00000000000000094156670658958205765,
        0.00000000000000167153952640716157721,
        -0.00000000000000055134268782182410852,
        -0.00000000000000022368651572006617795,
        0.00000000000000017486948976520089209,
        0.00000000000000000206518666352329750,
        -0.00000000000000003973060018130712479,
        0.00000000000000001154836935724892335,
        0.00000000000000000553906053678276421,
        -0.00000000000000000457174427396478267,
        0.00000000000000000026567111858284432,
        0.00000000000000000101599148154167823,
        -0.00000000000000000044821231272196246,
        -0.00000000000000000007959149661617295,
        0.00000000000000000014583615616165794,
        -0.00000000000000000004015127893061405,
        -0.00000000000000000002079152963743616,
        0.00000000000000000001972630449634388,
        -0.00000000000000000000336033404001683,
        -0.00000000000000000000376504832685507,
        0.00000000000000000000269935508825595,
        -0.00000000000000000000026985946069808,
        -0.00000000000000000000061794011788222,
        0.00000000000000000000038782693311711,
        -0.00000000000000000000002420094005071,
        -0.00000000000000000000009844051058925,
        0.00000000000000000000005954353358494,
        -0.00000000000000000000000361274446366,
        -0.00000000000000000000001552634578088,
        0.00000000000000000000000977819380304,
        -0.00000000000000000000000092239447509,
        -0.00000000000000000000000241545903934,
        0.00000000000000000000000169558652255,
        -0.00000000000000000000000026762408641,
        -0.00000000000000000000000036188116265,
        0.00000000000000000000000030372404951,
        -0.00000000000000000000000007422876903,
        -0.00000000000000000000000004930678544,
        0.00000000000000000000000005468790028,
        -0.00000000000000000000000001920315188,
        -0.00000000000000000000000000516335154,
        0.00000000000000000000000000957723167,
        -0.00000000000000000000000000463659079,
        -0.00000000000000000000000000004509226,
        0.00000000000000000000000000155617519,
        -0.00000000000000000000000000104156509,
        0.00000000000000000000000000019565323,
        0.00000000000000000000000000021335380,
        -0.00000000000000000000000000021461958,
        0.00000000000000000000000000007875791,
        0.00000000000000000000000000001713768,
        -0.00000000000000000000000000003917137,
        0.00000000000000000000000000002233559,
        -0.00000000000000000000000000000269383,
        -0.00000000000000000000000000000577764,
        0.00000000000000000000000000000519650,
        -0.00000000000000000000000000000183361,
        -0.00000000000000000000000000000045763,
        0.00000000000000000000000000000099235,
        -0.00000000000000000000000000000058938,
        0.00000000000000000000000000000009568,
        0.00000000000000000000000000000013758,
        -0.00000000000000000000000000000014066,
        0.00000000000000000000000000000005964,
        0.00000000000000000000000000000000437
    };
    static ityp btr = -2.09383632135605431360096498526268;
    ityp eta;
    static dim_typ nbif = 0;
    static dim_typ nbif2 = 0;
    static dim_typ nbig = 0;
    static dim_typ nbig2 = 0;
    static dim_typ nbip1 = 0;
    static dim_typ nbip2 = 0;
    ityp phi;
    ityp sqrtx;
    ityp value;
    ityp x2;
    static ityp x2sml = 0.00;
    ityp x3;
    static ityp x3sml = 0.00;
    static ityp x32sml = 0.00;
    static ityp xbig = 0.00;
    ityp xn;
    ityp z;

    if ( nbif == 0 )
    {
    eta = 0.10 * r8_mach ( 3 );
        nbif = r8_inits ( bifcs, 13, eta );
        nbig = r8_inits ( bigcs, 13, eta );
        nbif2 = r8_inits ( bif2cs, 15, eta );
        nbig2 = r8_inits ( big2cs, 16, eta );
        nbip1 = r8_inits ( bip1cs, 47, eta );
        nbip2 = r8_inits ( bip2cs, 88, eta );
        x2sml = sqrt ( eta );
        x3sml = pow ( eta, 0.3333 );
        x32sml = 1.3104 * x3sml * x3sml;
        xbig = pow ( r8_mach ( 2 ), 0.6666 );
    }

    if ( x < -1.00 )
    {
        r8_admp ( x, &xn, &phi );
        value = xn * sin ( phi );
    }
    else if ( abs ( x ) <= x2sml )
    {
        x2 = x3 = 0.00;
        value = x2 * ( r8_csevl ( x3, bifcs, nbif )+ 0.25 ) + r8_csevl ( x3, bigcs, nbig ) + 0.50;
        if ( x32sml < x )
            value *= exp ( - 2.00 * x * sqrt ( x ) / 3.00 );
    }
    else if ( abs ( x ) <= x3sml )
    {
        x2 = x * x;
        x3 = 0.00;
        value = x2 * ( r8_csevl ( x3, bifcs, nbif )+ 0.25 ) + r8_csevl ( x3, bigcs, nbig ) + 0.50;
        if ( x32sml < x )
            value *= exp ( - 2.00 * x * sqrt ( x ) / 3.00 );
    }
    else if ( x <= 1.00 )
    {
        x2 = x * x;
        x3 = x * x * x;
        value = x2 * ( r8_csevl ( x3, bifcs, nbif )+ 0.25 ) + r8_csevl ( x3, bigcs, nbig ) + 0.5;
        if ( x32sml < x )
            value *= exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
    }
    else if ( x <= 2.00 )
    {
        z = ( 2.0 * x * x * x - 9.00 ) / 7.00;
        value = exp ( - 2.00 * x * sqrt ( x ) / 3.00 )* ( x * x * ( 0.25 + r8_csevl ( z, bif2cs, nbif2 ) )+ 0.50 + r8_csevl ( z, big2cs, nbig2 ) );
    }
    else if ( x <= 4.00 )
    {
        sqrtx = sqrt ( x );
        z = atr / x / sqrtx + btr;
        value = ( 0.625 + r8_csevl ( z, bip1cs, nbip1 ) ) * sqrt ( sqrtx );
    }
    else if ( x <= xbig )
    {
        sqrtx = sqrt ( x );
        z = 16.00 / x / sqrtx - 1.00;
        value = ( 0.625 + r8_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx );
    }
    else
    {
        sqrtx = sqrt ( x );
        z = -1.00;
        value = ( 0.625 + r8_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx );
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_bie ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BIE evaluates the exponentially scaled Airy function Bi of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BIE, the exponentially scaled Airy function Bi of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp atr = 8.75069057084843450880771988210148;
    static ityp bif2cs[15] =
    {
        +0.0998457269381604104468284257993,
        +0.47862497786300553772211467318231,
        +0.25155211960433011771324415436675E-01,
        +0.58206938852326456396515697872216E-03,
        +0.74997659644377865943861457378217E-05,
        +0.61346028703493836681403010356474E-07,
        +0.34627538851480632900434268733359E-09,
        +0.14288910080270254287770846748931E-11,
        +0.44962704298334641895056472179200E-14,
        +0.11142323065833011708428300106666E-16,
        +0.22304791066175002081517866666666E-19,
        +0.36815778736393142842922666666666E-22,
        +0.50960868449338261333333333333333E-25,
        +0.60003386926288554666666666666666E-28,
        +0.60827497446570666666666666666666E-31
    };
    static ityp bifcs[13] =
    {
        -0.16730216471986649483537423928176E-01,
        +0.10252335834249445611426362777757,
        +0.17083092507381516539429650242013E-02,
        +0.11862545467744681179216459210040E-04,
        +0.44932907017792133694531887927242E-07,
        +0.10698207143387889067567767663628E-09,
        +0.17480643399771824706010517628573E-12,
        +0.20810231071761711025881891834399E-15,
        +0.18849814695665416509927971733333E-18,
        +0.13425779173097804625882666666666E-21,
        +0.77159593429658887893333333333333E-25,
        +0.36533879617478566399999999999999E-28,
        +0.14497565927953066666666666666666E-31
    };
    static ityp big2cs[15] =
    {
        +0.033305662145514340465176188111647,
        +0.161309215123197067613287532084943,
        +0.631900730961342869121615634921173E-02,
        +0.118790456816251736389780192304567E-03,
        +0.130453458862002656147116485012843E-05,
        +0.937412599553521729546809615508936E-08,
        +0.474580188674725153788510169834595E-10,
        +0.178310726509481399800065667560946E-12,
        +0.516759192784958180374276356640000E-15,
        +0.119004508386827125129496251733333E-17,
        +0.222982880666403517277063466666666E-20,
        +0.346551923027689419722666666666666E-23,
        +0.453926336320504514133333333333333E-26,
        +0.507884996513522346666666666666666E-29,
        +0.491020674696533333333333333333333E-32
    };
    static ityp bigcs[13] =
    {
        +0.22466223248574522283468220139024E-01,
        +0.37364775453019545441727561666752E-01,
        +0.44476218957212285696215294326639E-03,
        +0.24708075636329384245494591948882E-05,
        +0.79191353395149635134862426285596E-08,
        +0.16498079851827779880887872402706E-10,
        +0.24119906664835455909247501122841E-13,
        +0.26103736236091436985184781269333E-16,
        +0.21753082977160323853123792000000E-19,
        +0.14386946400390433219483733333333E-22,
        +0.77349125612083468629333333333333E-26,
        +0.34469292033849002666666666666666E-29,
        +0.12938919273216000000000000000000E-32
    };
    static ityp bip1cs[47] =
    {
        -0.83220474779434474687471864707973E-01,
        +0.11461189273711742889920226128031E-01,
        +0.42896440718911509494134472566635E-03,
        -0.14906639379950514017847677732954E-03,
        -0.13076597267876290663136340998881E-04,
        +0.63275983961030344754535716032494E-05,
        -0.42226696982681924884778515889433E-06,
        -0.19147186298654689632835494181277E-06,
        +0.64531062845583173611038157880934E-07,
        -0.78448546771397719289748310448628E-08,
        -0.96077216623785085879198533565432E-09,
        +0.70004713316443966339006074402068E-09,
        -0.17731789132814932022083128056698E-09,
        +0.22720894783465236347282126389311E-10,
        +0.16540456313972049847032860681891E-11,
        -0.18517125559292316390755369896693E-11,
        +0.59576312477117290165680715534277E-12,
        -0.12194348147346564781055769498986E-12,
        +0.13347869253513048815386347813597E-13,
        +0.17278311524339746664384792889731E-14,
        -0.14590732013016720735268871713166E-14,
        +0.49010319927115819978994989520104E-15,
        -0.11556545519261548129262972762521E-15,
        +0.19098807367072411430671732441524E-16,
        -0.11768966854492179886913995957862E-17,
        -0.63271925149530064474537459677047E-18,
        +0.33861838880715361614130191322316E-18,
        -0.10725825321758625254992162219622E-18,
        +0.25995709605617169284786933115562E-19,
        -0.48477583571081193660962309494101E-20,
        +0.55298913982121625361505513198933E-21,
        +0.49421660826069471371748197444266E-22,
        -0.55162121924145707458069720814933E-22,
        +0.21437560417632550086631884499626E-22,
        -0.61910313387655605798785061137066E-23,
        +0.14629362707391245659830967336959E-23,
        -0.27918484471059005576177866069333E-24,
        +0.36455703168570246150906795349333E-25,
        +0.58511821906188711839382459733333E-27,
        -0.24946950487566510969745047551999E-26,
        +0.10979323980338380977919579477333E-26,
        -0.34743388345961115015034088106666E-27,
        +0.91373402635349697363171082240000E-28,
        -0.20510352728210629186247720959999E-28,
        +0.37976985698546461748651622399999E-29,
        -0.48479458497755565887848448000000E-30,
        -0.10558306941230714314205866666666E-31
    };
    static ityp bip2cs[88] =
    {
        -0.11359673758598867913797310895527,
        +0.41381473947881595760052081171444E-02,
        +0.13534706221193329857696921727508E-03,
        +0.10427316653015353405887183456780E-04,
        +0.13474954767849907889589911958925E-05,
        +0.16965374054383983356062511163756E-06,
        -0.10096500865641624301366228396373E-07,
        -0.16729119493778475127836973095943E-07,
        -0.45815364485068383217152795613391E-08,
        +0.37366813665655477274064749384284E-09,
        +0.57669303201452448119584643502111E-09,
        +0.62181265087850324095393408792371E-10,
        -0.63294120282743068241589177281354E-10,
        -0.14915047908598767633999091989487E-10,
        +0.78896213942486771938172394294891E-11,
        +0.24960513721857797984888064000127E-11,
        -0.12130075287291659477746664734814E-11,
        -0.37404939108727277887343460402716E-12,
        +0.22377278140321476798783446931091E-12,
        +0.47490296312192466341986077472514E-13,
        -0.45261607991821224810605655831294E-13,
        -0.30172271841986072645112245876020E-14,
        +0.91058603558754058327592683478908E-14,
        -0.98149238033807062926643864207709E-15,
        -0.16429400647889465253601245251589E-14,
        +0.55334834214274215451182114635164E-15,
        +0.21750479864482655984374381998156E-15,
        -0.17379236200220656971287029558087E-15,
        -0.10470023471443714959283909313604E-17,
        +0.39219145986056386925441403311462E-16,
        -0.11621293686345196925824005665910E-16,
        -0.54027474491754245533735411307773E-17,
        +0.45441582123884610882675428553304E-17,
        -0.28775599625221075729427585480086E-18,
        -0.10017340927225341243596162960440E-17,
        +0.44823931215068369856332561906313E-18,
        +0.76135968654908942328948982366775E-19,
        -0.14448324094881347238956060145422E-18,
        +0.40460859449205362251624847392112E-19,
        +0.20321085700338446891325190707277E-19,
        -0.19602795471446798718272758041962E-19,
        +0.34273038443944824263518958211738E-20,
        +0.37023705853905135480024651593154E-20,
        -0.26879595172041591131400332966712E-20,
        +0.28121678463531712209714454683364E-21,
        +0.60933963636177797173271119680329E-21,
        -0.38666621897150844994172977893413E-21,
        +0.25989331253566943450895651927228E-22,
        +0.97194393622938503767281175216084E-22,
        -0.59392817834375098415630478204591E-22,
        +0.38864949977113015409591960439444E-23,
        +0.15334307393617272869721512868769E-22,
        -0.97513555209762624036336521409724E-23,
        +0.96340644440489471424741339383726E-24,
        +0.23841999400208880109946748792454E-23,
        -0.16896986315019706184848044205207E-23,
        +0.27352715888928361222578444801478E-24,
        +0.35660016185409578960111685025730E-24,
        -0.30234026608258827249534280666954E-24,
        +0.75002041605973930653144204823232E-25,
        +0.48403287575851388827455319838748E-25,
        -0.54364137654447888432698010297766E-25,
        +0.19281214470820962653345978809756E-25,
        +0.50116355020532656659611814172172E-26,
        -0.95040744582693253786034620869972E-26,
        +0.46372646157101975948696332245611E-26,
        +0.21177170704466954163768170577046E-28,
        -0.15404850268168594303692204548726E-26,
        +0.10387944293201213662047889194441E-26,
        -0.19890078156915416751316728235153E-27,
        -0.21022173878658495471177044522532E-27,
        +0.21353099724525793150633356670491E-27,
        -0.79040810747961342319023537632627E-28,
        -0.16575359960435585049973741763592E-28,
        +0.38868342850124112587625586496537E-28,
        -0.22309237330896866182621562424717E-28,
        +0.27777244420176260265625977404382E-29,
        +0.57078543472657725368712433782772E-29,
        -0.51743084445303852800173371555280E-29,
        +0.18413280751095837198450927071569E-29,
        +0.44422562390957094598544071068647E-30,
        -0.98504142639629801547464958226943E-30,
        +0.58857201353585104884754198881995E-30,
        -0.97636075440429787961402312628595E-31,
        -0.13581011996074695047063597884122E-30,
        +0.13999743518492413270568048380345E-30,
        -0.59754904545248477620884562981118E-31,
        -0.40391653875428313641045327529856E-32
    };
    static ityp btr = -2.09383632135605431360096498526268;
    ityp eta;
    static dim_typ nbif = 0;
    static dim_typ nbif2 = 0;
    static dim_typ nbig = 0;
    static dim_typ nbig2 = 0;
    static dim_typ nbip1 = 0;
    static dim_typ nbip2 = 0;
    ityp sqrtx;
    ityp theta;
    ityp value;
    static ityp x32sml;
    static ityp x3sml;
    static ityp xbig;
    ityp xm;
    ityp z;

    if ( nbif == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        nbif = r8_inits ( bifcs, 13, eta );
        nbig = r8_inits ( bigcs, 13, eta );
        nbif2 = r8_inits ( bif2cs, 15, eta );
        nbig2 = r8_inits ( big2cs, 15, eta );
        nbip1 = r8_inits ( bip1cs, 47, eta );
        nbip2 = r8_inits ( bip2cs, 88, eta );
        x3sml = pow ( eta, 0.3333 );
        x32sml = 1.3104 * x3sml * x3sml;
        xbig = pow ( r8_mach ( 2 ), 0.6666 );
    }

    if ( x < - 1.00 )
    {
        r8_aimp ( x, &xm, &theta );
        value = xm * sin ( theta );
    }
    else if ( abs ( x ) <= x3sml )
    {
        z = 0.00;
        value = 0.625 + r8_csevl ( z, bifcs, nbif )+ x * ( 0.4375 + r8_csevl ( z, bigcs, nbig ) );
        if (  x32sml <= x )
            value *= exp ( - 2.00 * x * sqrt ( x ) / 3.00 );
    }
    else if ( x <= 1.00 )
    {
        z = x * x * x;
        value = 0.625 + r8_csevl ( z, bifcs, nbif )+ x * ( 0.4375 + r8_csevl ( z, bigcs, nbig ) );
        if (  x32sml <= x )
            value *= exp ( - 2.00 * x * sqrt ( x ) / 3.00 );
    }
    else if ( x <= 2.00 )
    {
        z = ( 2.00 * x * x * x - 9.00 ) / 7.00;
        value = exp ( - 2.00 * x * sqrt ( x ) / 3.00 )* ( 1.125 + r8_csevl ( z, bif2cs, nbif2 )+ x * ( 0.625 + r8_csevl ( z, big2cs, nbig2 ) ) );
    }
    else if ( x <= 4.00 )
    {
        sqrtx = sqrt ( x );
        z = atr / x / sqrtx + btr;
        value = ( 0.625 + r8_csevl ( z, bip1cs, nbip1 ) ) / sqrt ( sqrtx );
    }
    else if ( x < xbig )
    {
        sqrtx = sqrt ( x );
        z = 16.00 / ( x * sqrtx ) - 1.00;
        value = ( 0.625 + r8_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx );
    }
    else
    {
        sqrtx = sqrt ( x );
        z = - 1.00;
        value = ( 0.625 + r8_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx );
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_binom ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BINOM evaluates the binomial coefficient using R8 arithmetic.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, int N, M, the arguments.
    Output, double R8_BINOM, the binomial coefficient.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ m = a_data[1];
	
    static ityp bilnmx = 0.00;
    ityp corr;
    static ityp fintmx = 0.00;
    int i, k;
    static ityp sq2pil = 0.91893853320467274178032973640562;
    ityp value;
    ityp xk;
    ityp xn;
    ityp xnk;

    if ( bilnmx == 0.00 )
    {
        bilnmx = log ( r8_mach ( 2 ) ) - 0.0001;
        fintmx = 0.90 / r8_mach ( 3 );
    }

    if ( n < 0 || m < 0 || n < m )
    {
    	result = MAX_VAL;
        return &result;
    }

    k = MIN ( m, n - m );

    if ( k <= 20 &&( ityp ) ( k ) * log ( ( ityp ) ( MAX ( n, 1 ) ) ) <= bilnmx )
    {
        value = 1.00;
        for ( i = 1; i <= k; ++i )
            value *= ( ityp ) ( n - i + 1 ) / ( ityp ) ( i );
    }
    else
    {
        if ( k < 9 )
        {
        	result = MAX_VAL;
            return &result;
        }

        xn = ( ityp ) ( n + 1 );
        xk = ( ityp ) ( k + 1 );
        xnk = ( ityp ) ( n - k + 1 );

        corr = r8_lgmc ( xn ) - r8_lgmc ( xk ) - r8_lgmc ( xnk );
        value = xk * log ( xnk / xk )- xn * r8_lnrel ( - ( xk - 1.00 ) / xn )- 0.50 * log ( xn * xnk / xk ) + 1.00 - sq2pil + corr;

        if ( bilnmx < value )
        {
        	result = MAX_VAL;
            return &result;
        }
        value = exp ( value );
    }

    if ( value < fintmx )
        value = r8_aint ( value + 0.50 );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_chi ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CHI evaluates the hyperbolic cosine integral of an R8 argument.
  Discussion:
    The hyperbolic cosine integral is defined by
      CHI(X) = gamma + log ( x )
        + integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
    where gamma is Euler's constant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_CHI, the hyperbolic cosine integral
    evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	result = 0.50 * ( r8_ei ( x ) - r8_e1 ( x ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_chu ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CHU evaluates the confluent hypergeometric function of R8 arguments.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, B, the parameters.
    Input, double X, the argument.
    Output, double R8_CHU, the function value.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp b = a_data[1];
	const register ityp x = a_data[2];
	
    ityp a0;
    ityp aintb;
    ityp alnx;
    ityp b0;
    ityp beps;
    ityp c0;
    static ityp eps = 0.00;
    ityp factor;
    ityp gamri1;
    ityp gamrni;
    dim_typ i;
    dim_typ istrt;
    dim_typ m;
    dim_typ n;
    ityp pch1ai;
    ityp pch1i;
    ityp pochai;
    ityp sum;
    ityp t;
    ityp value;
    ityp xeps1;
    ityp xi;
    ityp xi1;
    ityp xn;
    ityp xtoeps;

    if ( eps == 0.00 )
        eps = r8_mach ( 3 );

    if ( x < 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( x == 0.00 )
    {
        if ( 1.00 <= b )
        {
        	result = MAX_VAL;
        	return &result;
        }
        
        result = r8_gamma ( 1.0 - b ) / r8_gamma ( 1.00 + a - b );
        return &result;
    }

    if ( MAX ( abs ( a ), 1.00 )* MAX ( abs ( 1.00 + a - b ), 1.00 ) < 0.99 * abs ( x ) )
    {
    	result = pow ( x, - a ) * r8_chu_scaled ( a, b, x );
        return &result;
    }
    /*
    The ascending series will be used, because the descending rational
    approximation (which is based on the asymptotic series) is unstable.
    */
    aintb = r8_aint(0.00<=b?b+0.50:b-0.50);
    beps = b - aintb;
    n = ( dim_typ ) aintb;
    alnx = log ( x );
    xtoeps = exp ( - beps * alnx );
    /*
    Evaluate the finite sum.

    Consider the case b < 1.0 first.
    */
    if ( n < 1 )
    {
        sum = t = 1.00;
        m = - n;
        for ( i = 1; i <= m; ++i )
        {
            xi1 = ( ityp ) ( i - 1 );
            t *= ( a + xi1 ) * x / ( ( b + xi1 ) * ( xi1 + 1.0 ) );
            sum += t;
        }
        sum = r8_poch ( 1.00 + a - b, - a ) * sum;
    }
    /*
    Now consider the case 1 <= b.
    */
    else
    {
        sum = 0.00;
        m = n - 2;

        if ( 0 <= m )
        {
            t = sum = 1.00;

            for ( i = 1; i <= m; ++i )
            {
                xi = ( ityp ) ( i );
                t *= ( a - b + xi ) * x / ( ( 1.00 - b + xi ) * xi );
                sum +=  t;
            }

            sum = r8_gamma ( b - 1.00 ) * r8_gamr ( a )* pow ( x, ( ityp ) ( 1 - n ) ) * xtoeps * sum;
        }
    }
    /*
    Next evaluate the infinite sum.
    */
    istrt = 0 + (n<1)*(1-n);
    xi = ( ityp ) ( istrt );
    factor = r8_mop ( n ) * r8_gamr ( 1.0 + a - b ) * pow ( x, xi );

    if ( beps != 0.00 )
        factor *= beps * M_PI / sin ( beps * M_PI );

    pochai = r8_poch ( a, xi );
    gamri1 = r8_gamr ( xi + 1.00 );
    gamrni = r8_gamr ( aintb + xi );
    b0 = factor * r8_poch ( a, xi - beps )* gamrni * r8_gamr ( xi + 1.00 - beps );
    /*
    x^(-beps) is close to 1.0, so we must be careful in evaluating the
    differences.
    */
    if ( abs ( xtoeps - 1.00 ) <= 0.50 )
    {
        pch1ai = r8_poch1 ( a + xi, -beps );
        pch1i = r8_poch1 ( xi + 1.00 - beps, beps );
        c0 = factor * pochai * gamrni * gamri1 * (- r8_poch1 ( b + xi,- beps ) + pch1ai- pch1i + beps * pch1ai * pch1i );
        /*
        xeps1 = (1.0 - x^(-beps))/beps = (x^(-beps) - 1.0)/(-beps)
        */
        xeps1 = alnx* r8_exprel ( - beps * alnx );

        value = sum + c0 + xeps1 * b0;
        xn = ( ityp ) ( n );

        for ( i = 1; i <= 1000; ++i )
        {
            xi = ( ityp ) ( istrt + i );
            xi1 = ( ityp ) ( istrt + i - 1 );
            b0 = ( a + xi1 - beps ) * b0 * x/ ( ( xn + xi1 ) * ( xi - beps ) );
            c0 = ( a + xi1 ) * c0 * x / ( ( b + xi1) * xi )- ( ( a - 1.0 ) * ( xn + 2.0 * xi - 1.0 )+ xi * ( xi - beps ) ) * b0/ ( xi * ( b + xi1 ) * ( a + xi1 - beps ) );
            t = c0 + xeps1 * b0;
            value += t;
            if ( abs ( t ) < eps * abs ( value ) )
            {
            	result = value;
            }
                return &result;
        }

		result = MAX_VAL;  
        return &result;
    }
    /*
    x^(-beps) is very different from 1.0, so the straightforward
    formulation is stable.
    */
    a0 = factor * pochai * r8_gamr ( b + xi ) * gamri1 / beps;
    b0 = xtoeps * b0 / beps;

    value = sum + a0 - b0;

    for ( i = 1; i <= 1000; ++i )
    {
        xi = ( ityp ) ( istrt + i );
        xi1 = ( ityp ) ( istrt + i - 1 );
        a0 = ( a + xi1 ) * a0 * x / ( ( b + xi1 ) * xi );
        b0 = ( a + xi1 - beps ) * b0 * x/ ( ( aintb + xi1 ) * ( xi - beps ) );
        t = a0 - b0;
        value += t;
        if ( abs ( t ) < eps * abs ( value ) )
        {
        	result = value;
            return &result;
        }
    }
    
    result = MAX_VAL;  
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_chu_scaled ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CHU_SCALED: scaled confluent hypergeometric function of R8 arguments.
  Discussion:
    Evaluate, for large z, z^a * u(a,b,z)  where U is the logarithmic
    confluent hypergeometric function.  A rational approximation due to
    Y L Luke is used.  When U is not in the asymptotic region, that is, when A
    or B is large compared with Z, considerable significance loss occurs.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, B, the parameters.
    Input, double Z, the argument.
    Output, double R8_CHU_SCALED, the function value.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp b = a_data[1];
	const register ityp z = a_data[2];
	
    ityp aa[4];
    ityp ab;
    ityp anbn;
    ityp bb[4];
    ityp bp;
    ityp c2;
    ityp ct1;
    ityp ct2;
    ityp ct3;
    ityp d1z;
    static ityp eps = 0.00;
    ityp g1;
    ityp g2;
    ityp g3;
    dim_typ i, j;
    ityp sab;
    static ityp sqeps = 0.0;
    ityp value;
    ityp x2i1;

    if ( eps == 0.00 )
    {
        eps = 4.00 * r8_mach ( 4 );
        sqeps = sqrt ( r8_mach ( 4 ) );
    }

    bp = 1.00 + a - b;
    ab = a * bp;
    ct2 = 2.00 * ( z - ab );
    sab = a + bp;

    bb[0] = aa[0] = 1.00;

    ct3 = sab + 1.00 + ab;
    bb[1] = 1.00 + 2.00 * z / ct3;
    aa[1] = 1.00 + ct2 / ct3;

    anbn = ct3 + sab + 3.00;
    ct1 = 1.00 + 2.00 * z / anbn;
    bb[2] = 1.00 + 6.00 * ct1 * z / ct3;
    aa[2] = 1.00 + 6.00 * ab / anbn + 3.00 * ct1 * ct2 / ct3;

    for ( i = 4; i <= 300; ++i )
    {
        x2i1 = ( ityp ) ( (i<<1) - 3 );
        ct1 = x2i1 / ( x2i1 - 2.00 );
        anbn = anbn + x2i1 + sab;
        ct2 = ( x2i1 - 1.00 ) /anbn;
        c2 = x2i1 * ct2 - 1.00;
        d1z = x2i1 * 2.00 * z / anbn;

        ct3 = sab * ct2;
        g1 = d1z + ct1 * ( c2 + ct3 );
        g2 = d1z - c2;
        g3 = ct1 * ( 1.00 - ct3 - 2.0 * ct2 );

        bb[3] = g1 * bb[2] + g2 * bb[1] + g3 * bb[0];
        aa[3] = g1 * aa[2] + g2 * aa[1] + g3 * aa[0];

        value = aa[3] / bb[3];

        if ( abs ( value - aa[0] / bb[0] ) < eps * abs ( value ) )
        {
        	result = value;
            return &result;
        }

        #pragma omp parallel for num_threads(3)
        for ( j = 0; j < 3; ++j )
        {
            aa[j] = aa[j+1];
            bb[j] = bb[j+1];
        }
    }
    
    result = MAX_VAL;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_ci ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CI evaluates the cosine integral Ci of an R8 argument.
  Discussion:
    The cosine integral is defined by
      CI(X) = - integral ( X <= T < Infinity ) ( cos ( T ) ) / T  dT
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_CI, the cosine integral Ci evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp cics[19] =
    {
        -0.34004281856055363156281076633129873,
        -1.03302166401177456807159271040163751,
        0.19388222659917082876715874606081709,
        -0.01918260436019865893946346270175301,
        0.00110789252584784967184098099266118,
        -0.00004157234558247208803840231814601,
        0.00000109278524300228715295578966285,
        -0.00000002123285954183465219601280329,
        0.00000000031733482164348544865129873,
        -0.00000000000376141547987683699381798,
        0.00000000000003622653488483964336956,
        -0.00000000000000028911528493651852433,
        0.00000000000000000194327860676494420,
        -0.00000000000000000001115183182650184,
        0.00000000000000000000005527858887706,
        -0.00000000000000000000000023907013943,
        0.00000000000000000000000000091001612,
        -0.00000000000000000000000000000307233,
        0.00000000000000000000000000000000926
    };
    ityp f;
    ityp g;
    static dim_typ nci = 0;
    ityp sinx;
    ityp value;
    static ityp xsml = 0.00;
    ityp y;

    if ( nci == 0 )
    {
        nci = r8_inits ( cics, 19, 0.10 * r8_mach ( 3 ) );
        xsml = sqrt ( r8_mach ( 3 ) );
    }

    if ( x <= 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }
    else if ( x <= xsml )
    {
        y = - 1.00;
        value = log ( x ) - 0.50 + r8_csevl ( y, cics, nci );
    }
    else if ( x <= 4.00 )
    {
        y = ( x * x - 8.00 ) * 0.125;
        value = log ( x ) - 0.50 + r8_csevl ( y, cics, nci );
    }
    else
    {
        r8_sifg ( x, &f, &g );
        sinx = sin ( x );
        value = f * sinx - g * cos ( x );
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_cin ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CIN evaluates the alternate cosine integral Cin of an R8 argument.
  Discussion:
    CIN(X) = gamma + log(X)
      + integral ( 0 <= T <= X ) ( cos ( T ) - 1 ) / T  dT
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_CIN, the cosine integral Cin evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp absx;
    static ityp cincs[18] =
    {
        0.37074501750909688741654801228564992,
        -0.05893574896364446831956864397363697,
        0.00538189642113569124048745326203340,
        -0.00029860052841962135319594906563410,
        0.00001095572575321620077031054467306,
        -0.00000028405454877346630491727187731,
        0.00000000546973994875384912457861806,
        -0.00000000008124187461318157083277452,
        0.00000000000095868593117706609013181,
        -0.00000000000000920266004392351031377,
        0.00000000000000007325887999017895024,
        -0.00000000000000000049143726675842909,
        0.00000000000000000000281577746753902,
        -0.00000000000000000000001393986788501,
        0.00000000000000000000000006022485646,
        -0.00000000000000000000000000022904717,
        0.00000000000000000000000000000077273,
        -0.00000000000000000000000000000000233
    };
    static ityp eul = 0.57721566490153286060651209008240;
    ityp f;
    ityp g;
    static dim_typ ncin = 0;
    ityp sinx;
    ityp value;
    static ityp xmin = 0.00;

    if ( ncin == 0 )
    {
        ncin = r8_inits ( cincs, 18, 0.10 * r8_mach ( 3 ) );
        xmin = sqrt ( r8_mach ( 1 ) );
    }

    absx = abs ( x );


    if ( absx <= xmin )
        value = 0.00;
    else if ( absx <= 4.00 )
        value = r8_csevl ( ( x * x - 8.00 ) * 0.125, cincs, ncin ) * x * x;
    else
    {
        r8_sifg ( absx, &f, &g );
        sinx = sin ( absx );
        value = - f * sinx + g * cos ( absx ) + log ( absx ) + eul;
    }
    
    result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_cinh ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CINH: alternate hyperbolic cosine integral Cinh of an R8 argument.
  Discussion:
    Cinh ( x ) = Integral ( 0 <= t <= x ) ( cosh ( t ) - 1 ) dt / t
    The original text of this program had a mistake:
      y = x * x / 9.0 - 1.0
    has been corrected to
      y = x * x / 4.5 - 1.0
    JVB, 27 March 2010
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_CINH, the hyperbolic cosine integral Cinh
    evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp absx;
    static ityp cinhcs[16] =
    {
        0.1093291636520734431407425199795917,
        0.0573928847550379676445323429825108,
        0.0028095756978830353416404208940774,
        0.0000828780840721356655731765069792,
        0.0000016278596173914185577726018815,
        0.0000000227809519255856619859083591,
        0.0000000002384484842463059257284002,
        0.0000000000019360829780781957471028,
        0.0000000000000125453698328172559683,
        0.0000000000000000663637449497262300,
        0.0000000000000000002919639263594744,
        0.0000000000000000000010849123956107,
        0.0000000000000000000000034499080805,
        0.0000000000000000000000000094936664,
        0.0000000000000000000000000000228291,
        0.0000000000000000000000000000000484
    };
    static ityp eul = 0.57721566490153286060651209008240;
    static dim_typ ncinh = 0;
    ityp value;
    static ityp xmin = 0.00;
    static ityp xsml = 0.00;
    ityp y;

    if ( ncinh == 0 )
    {
        ncinh = r8_inits ( cinhcs, 16, 0.10 * r8_mach ( 3 ) );
        xsml = sqrt ( r8_mach ( 3 ) );
        xmin = 2.00 * sqrt ( r8_mach ( 1 ) );
    }

    absx = abs ( x );

    if ( x == 0.00  || absx <= xmin)
        value = 0.00;
    else if ( x <= xsml )
    {
        y = - 1.00;
        value = x * x * ( 0.25 + r8_csevl ( y, cinhcs, ncinh ) );
    }
    else if ( x <= 3.00 )
    {
        y = x * x / 4.50 - 1.00;
        value = x * x * ( 0.25 + r8_csevl ( y, cinhcs, ncinh ) );
    }
    else
        value = r8_chi ( absx ) - eul - log ( absx );
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_cos_deg ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_COS_DEG evaluates the cosine of an R8 argument in degrees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument in degrees.
    Output, double R8_COS_DEG, the cosine of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    dim_typ n;
    static ityp raddeg = 0.017453292519943295769236907684886;
    ityp value =cos ( raddeg * x );

    if ( fmod ( x, 90.00 ) == 0.00 )
    {
        n = ( dim_typ ) ( abs ( x ) / 90.00 + 0.50 );
        n = ( n % 2 );

        if ( n == 1 )
            value = 0.00;
        else if ( value < 0.00 )
            value = - 1.00;
        else
            value = + 1.00;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_dawson ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_DAWSON evaluates Dawson's integral of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_DAWSON, the value of Dawson's integral at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp daw2cs[45] =
    {
        -0.56886544105215527114160533733674E-01,
        -0.31811346996168131279322878048822,
        +0.20873845413642236789741580198858,
        -0.12475409913779131214073498314784,
        +0.67869305186676777092847516423676E-01,
        -0.33659144895270939503068230966587E-01,
        +0.15260781271987971743682460381640E-01,
        -0.63483709625962148230586094788535E-02,
        +0.24326740920748520596865966109343E-02,
        -0.86219541491065032038526983549637E-03,
        +0.28376573336321625302857636538295E-03,
        -0.87057549874170423699396581464335E-04,
        +0.24986849985481658331800044137276E-04,
        -0.67319286764160294344603050339520E-05,
        +0.17078578785573543710504524047844E-05,
        -0.40917551226475381271896592490038E-06,
        +0.92828292216755773260751785312273E-07,
        -0.19991403610147617829845096332198E-07,
        +0.40963490644082195241210487868917E-08,
        -0.80032409540993168075706781753561E-09,
        +0.14938503128761465059143225550110E-09,
        -0.26687999885622329284924651063339E-10,
        +0.45712216985159458151405617724103E-11,
        -0.75187305222043565872243727326771E-12,
        +0.11893100052629681879029828987302E-12,
        -0.18116907933852346973490318263084E-13,
        +0.26611733684358969193001612199626E-14,
        -0.37738863052129419795444109905930E-15,
        +0.51727953789087172679680082229329E-16,
        -0.68603684084077500979419564670102E-17,
        +0.88123751354161071806469337321745E-18,
        -0.10974248249996606292106299624652E-18,
        +0.13261199326367178513595545891635E-19,
        -0.15562732768137380785488776571562E-20,
        +0.17751425583655720607833415570773E-21,
        -0.19695006967006578384953608765439E-22,
        +0.21270074896998699661924010120533E-23,
        -0.22375398124627973794182113962666E-24,
        +0.22942768578582348946971383125333E-25,
        -0.22943788846552928693329592319999E-26,
        +0.22391702100592453618342297600000E-27,
        -0.21338230616608897703678225066666E-28,
        +0.19866196585123531518028458666666E-29,
        -0.18079295866694391771955199999999E-30,
        +0.16090686015283030305450666666666E-31
    };
    static ityp dawacs[75] =
    {
        +0.1690485637765703755422637438849E-01,
        +0.8683252278406957990536107850768E-02,
        +0.2424864042417715453277703459889E-03,
        +0.1261182399572690001651949240377E-04,
        +0.1066453314636176955705691125906E-05,
        +0.1358159794790727611348424505728E-06,
        +0.2171042356577298398904312744743E-07,
        +0.2867010501805295270343676804813E-08,
        -0.1901336393035820112282492378024E-09,
        -0.3097780484395201125532065774268E-09,
        -0.1029414876057509247398132286413E-09,
        -0.6260356459459576150417587283121E-11,
        +0.8563132497446451216262303166276E-11,
        +0.3033045148075659292976266276257E-11,
        -0.2523618306809291372630886938826E-12,
        -0.4210604795440664513175461934510E-12,
        -0.4431140826646238312143429452036E-13,
        +0.4911210272841205205940037065117E-13,
        +0.1235856242283903407076477954739E-13,
        -0.5788733199016569246955765071069E-14,
        -0.2282723294807358620978183957030E-14,
        +0.7637149411014126476312362917590E-15,
        +0.3851546883566811728777594002095E-15,
        -0.1199932056928290592803237283045E-15,
        -0.6313439150094572347334270285250E-16,
        +0.2239559965972975375254912790237E-16,
        +0.9987925830076495995132891200749E-17,
        -0.4681068274322495334536246507252E-17,
        -0.1436303644349721337241628751534E-17,
        +0.1020822731410541112977908032130E-17,
        +0.1538908873136092072837389822372E-18,
        -0.2189157877645793888894790926056E-18,
        +0.2156879197938651750392359152517E-20,
        +0.4370219827442449851134792557395E-19,
        -0.8234581460977207241098927905177E-20,
        -0.7498648721256466222903202835420E-20,
        +0.3282536720735671610957612930039E-20,
        +0.8858064309503921116076561515151E-21,
        -0.9185087111727002988094460531485E-21,
        +0.2978962223788748988314166045791E-22,
        +0.1972132136618471883159505468041E-21,
        -0.5974775596362906638089584995117E-22,
        -0.2834410031503850965443825182441E-22,
        +0.2209560791131554514777150489012E-22,
        -0.5439955741897144300079480307711E-25,
        -0.5213549243294848668017136696470E-23,
        +0.1702350556813114199065671499076E-23,
        +0.6917400860836148343022185660197E-24,
        -0.6540941793002752512239445125802E-24,
        +0.6093576580439328960371824654636E-25,
        +0.1408070432905187461501945080272E-24,
        -0.6785886121054846331167674943755E-25,
        -0.9799732036214295711741583102225E-26,
        +0.2121244903099041332598960939160E-25,
        -0.5954455022548790938238802154487E-26,
        -0.3093088861875470177838847232049E-26,
        +0.2854389216344524682400691986104E-26,
        -0.3951289447379305566023477271811E-27,
        -0.5906000648607628478116840894453E-27,
        +0.3670236964668687003647889980609E-27,
        -0.4839958238042276256598303038941E-29,
        -0.9799265984210443869597404017022E-28,
        +0.4684773732612130606158908804300E-28,
        +0.5030877696993461051647667603155E-29,
        -0.1547395051706028239247552068295E-28,
        +0.6112180185086419243976005662714E-29,
        +0.1357913399124811650343602736158E-29,
        -0.2417687752768673088385304299044E-29,
        +0.8369074582074298945292887587291E-30,
        +0.2665413042788979165838319401566E-30,
        -0.3811653692354890336935691003712E-30,
        +0.1230054721884951464371706872585E-30,
        +0.4622506399041493508805536929983E-31,
        -0.6120087296881677722911435593001E-31,
        +0.1966024640193164686956230217896E-31
    };
    static ityp dawcs[21] =
    {
        -0.6351734375145949201065127736293E-02,
        -0.2294071479677386939899824125866,
        +0.2213050093908476441683979161786E-01,
        -0.1549265453892985046743057753375E-02,
        +0.8497327715684917456777542948066E-04,
        -0.3828266270972014924994099521309E-05,
        +0.1462854806250163197757148949539E-06,
        -0.4851982381825991798846715425114E-08,
        +0.1421463577759139790347568183304E-09,
        -0.3728836087920596525335493054088E-11,
        +0.8854942961778203370194565231369E-13,
        -0.1920757131350206355421648417493E-14,
        +0.3834325867246327588241074439253E-16,
        -0.7089154168175881633584099327999E-18,
        +0.1220552135889457674416901120000E-19,
        -0.1966204826605348760299451733333E-21,
        +0.2975845541376597189113173333333E-23,
        -0.4247069514800596951039999999999E-25,
        +0.5734270767391742798506666666666E-27,
        -0.7345836823178450261333333333333E-29,
        +0.8951937667516552533333333333333E-31
    };
    ityp eps;
    static dim_typ ntdaw = 0;
    static dim_typ ntdaw2 = 0;
    static dim_typ ntdawa = 0;
    ityp value;
    static ityp xbig = 0.00;
    static ityp xmax = 0.00;
    static ityp xsml = 0.00;
    ityp y;

    if ( ntdaw == 0 )
    {
        eps = r8_mach ( 3 );
        ntdaw  = r8_inits ( dawcs,  21, 0.10 * eps );
        ntdaw2 = r8_inits ( daw2cs, 45, 0.10 * eps );
        ntdawa = r8_inits ( dawacs, 75, 0.10 * eps );
        xsml = sqrt ( 1.50 * eps );
        xbig = sqrt ( 0.50 / eps );
        xmax = exp ( MIN ( - log ( 2.0 * r8_mach ( 1 ) ),
        log ( r8_mach ( 2 ) ) ) - 0.01 );
    }

    y = abs ( x );

    if ( y <= xsml )
        value = x;
    else if ( y <= 1.00 )
        value = x * ( 0.75 + r8_csevl ( 2.00 * y * y - 1.00, dawcs, ntdaw ) );
    else if ( y <= 4.0 )
        value = x * ( 0.25 + r8_csevl ( 0.125 * y * y - 1.00, daw2cs, ntdaw2 ) );
    else if ( y < xbig )
        value = ( 0.50 + r8_csevl ( 32.00 / y / y - 1.00, dawacs, ntdawa ) ) / x;
    else if ( y <= xmax )
        value = 0.50 / x;
    else
        value = 0.0;

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_e1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_E1 evaluates the exponential integral E1 for an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_E1, the exponential integral E1 evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp ae10cs[50] =
    {
        +0.3284394579616699087873844201881E-01,
        -0.1669920452031362851476184343387E-01,
        +0.2845284724361346807424899853252E-03,
        -0.7563944358516206489487866938533E-05,
        +0.2798971289450859157504843180879E-06,
        -0.1357901828534531069525563926255E-07,
        +0.8343596202040469255856102904906E-09,
        -0.6370971727640248438275242988532E-10,
        +0.6007247608811861235760831561584E-11,
        -0.7022876174679773590750626150088E-12,
        +0.1018302673703687693096652346883E-12,
        -0.1761812903430880040406309966422E-13,
        +0.3250828614235360694244030353877E-14,
        -0.5071770025505818678824872259044E-15,
        +0.1665177387043294298172486084156E-16,
        +0.3166753890797514400677003536555E-16,
        -0.1588403763664141515133118343538E-16,
        +0.4175513256138018833003034618484E-17,
        -0.2892347749707141906710714478852E-18,
        -0.2800625903396608103506340589669E-18,
        +0.1322938639539270903707580023781E-18,
        -0.1804447444177301627283887833557E-19,
        -0.7905384086522616076291644817604E-20,
        +0.4435711366369570103946235838027E-20,
        -0.4264103994978120868865309206555E-21,
        -0.3920101766937117541553713162048E-21,
        +0.1527378051343994266343752326971E-21,
        +0.1024849527049372339310308783117E-22,
        -0.2134907874771433576262711405882E-22,
        +0.3239139475160028267061694700366E-23,
        +0.2142183762299889954762643168296E-23,
        -0.8234609419601018414700348082312E-24,
        -0.1524652829645809479613694401140E-24,
        +0.1378208282460639134668480364325E-24,
        +0.2131311202833947879523224999253E-26,
        -0.2012649651526484121817466763127E-25,
        +0.1995535662263358016106311782673E-26,
        +0.2798995808984003464948686520319E-26,
        -0.5534511845389626637640819277823E-27,
        -0.3884995396159968861682544026146E-27,
        +0.1121304434507359382850680354679E-27,
        +0.5566568152423740948256563833514E-28,
        -0.2045482929810499700448533938176E-28,
        -0.8453813992712336233411457493674E-29,
        +0.3565758433431291562816111116287E-29,
        +0.1383653872125634705539949098871E-29,
        -0.6062167864451372436584533764778E-30,
        -0.2447198043989313267437655119189E-30,
        +0.1006850640933998348011548180480E-30,
        +0.4623685555014869015664341461674E-31
    };
    static ityp ae11cs[60] =
    {
        +0.20263150647078889499401236517381,
        -0.73655140991203130439536898728034E-01,
        +0.63909349118361915862753283840020E-02,
        -0.60797252705247911780653153363999E-03,
        -0.73706498620176629330681411493484E-04,
        +0.48732857449450183453464992488076E-04,
        -0.23837064840448290766588489460235E-05,
        -0.30518612628561521027027332246121E-05,
        +0.17050331572564559009688032992907E-06,
        +0.23834204527487747258601598136403E-06,
        +0.10781772556163166562596872364020E-07,
        -0.17955692847399102653642691446599E-07,
        -0.41284072341950457727912394640436E-08,
        +0.68622148588631968618346844526664E-09,
        +0.53130183120506356147602009675961E-09,
        +0.78796880261490694831305022893515E-10,
        -0.26261762329356522290341675271232E-10,
        -0.15483687636308261963125756294100E-10,
        -0.25818962377261390492802405122591E-11,
        +0.59542879191591072658903529959352E-12,
        +0.46451400387681525833784919321405E-12,
        +0.11557855023255861496288006203731E-12,
        -0.10475236870835799012317547189670E-14,
        -0.11896653502709004368104489260929E-13,
        -0.47749077490261778752643019349950E-14,
        -0.81077649615772777976249734754135E-15,
        +0.13435569250031554199376987998178E-15,
        +0.14134530022913106260248873881287E-15,
        +0.49451592573953173115520663232883E-16,
        +0.79884048480080665648858587399367E-17,
        -0.14008632188089809829248711935393E-17,
        -0.14814246958417372107722804001680E-17,
        -0.55826173646025601904010693937113E-18,
        -0.11442074542191647264783072544598E-18,
        +0.25371823879566853500524018479923E-20,
        +0.13205328154805359813278863389097E-19,
        +0.62930261081586809166287426789485E-20,
        +0.17688270424882713734999261332548E-20,
        +0.23266187985146045209674296887432E-21,
        -0.67803060811125233043773831844113E-22,
        -0.59440876959676373802874150531891E-22,
        -0.23618214531184415968532592503466E-22,
        -0.60214499724601478214168478744576E-23,
        -0.65517906474348299071370444144639E-24,
        +0.29388755297497724587042038699349E-24,
        +0.22601606200642115173215728758510E-24,
        +0.89534369245958628745091206873087E-25,
        +0.24015923471098457555772067457706E-25,
        +0.34118376888907172955666423043413E-26,
        -0.71617071694630342052355013345279E-27,
        -0.75620390659281725157928651980799E-27,
        -0.33774612157467324637952920780800E-27,
        -0.10479325703300941711526430332245E-27,
        -0.21654550252170342240854880201386E-28,
        -0.75297125745288269994689298432000E-30,
        +0.19103179392798935768638084000426E-29,
        +0.11492104966530338547790728833706E-29,
        +0.43896970582661751514410359193600E-30,
        +0.12320883239205686471647157725866E-30,
        +0.22220174457553175317538581162666E-31
    };
    static ityp ae12cs[41] =
    {
        +0.63629589796747038767129887806803,
        -0.13081168675067634385812671121135,
        -0.84367410213053930014487662129752E-02,
        +0.26568491531006685413029428068906E-02,
        +0.32822721781658133778792170142517E-03,
        -0.23783447771430248269579807851050E-04,
        -0.11439804308100055514447076797047E-04,
        -0.14405943433238338455239717699323E-05,
        +0.52415956651148829963772818061664E-08,
        +0.38407306407844323480979203059716E-07,
        +0.85880244860267195879660515759344E-08,
        +0.10219226625855003286339969553911E-08,
        +0.21749132323289724542821339805992E-10,
        -0.22090238142623144809523503811741E-10,
        -0.63457533544928753294383622208801E-11,
        -0.10837746566857661115340539732919E-11,
        -0.11909822872222586730262200440277E-12,
        -0.28438682389265590299508766008661E-14,
        +0.25080327026686769668587195487546E-14,
        +0.78729641528559842431597726421265E-15,
        +0.15475066347785217148484334637329E-15,
        +0.22575322831665075055272608197290E-16,
        +0.22233352867266608760281380836693E-17,
        +0.16967819563544153513464194662399E-19,
        -0.57608316255947682105310087304533E-19,
        -0.17591235774646878055625369408853E-19,
        -0.36286056375103174394755328682666E-20,
        -0.59235569797328991652558143488000E-21,
        -0.76030380926310191114429136895999E-22,
        -0.62547843521711763842641428479999E-23,
        +0.25483360759307648606037606400000E-24,
        +0.25598615731739857020168874666666E-24,
        +0.71376239357899318800207052800000E-25,
        +0.14703759939567568181578956800000E-25,
        +0.25105524765386733555198634666666E-26,
        +0.35886666387790890886583637333333E-27,
        +0.39886035156771301763317759999999E-28,
        +0.21763676947356220478805333333333E-29,
        -0.46146998487618942367607466666666E-30,
        -0.20713517877481987707153066666666E-30,
        -0.51890378563534371596970666666666E-31
    };
    static ityp ae13cs[50] =
    {
        -0.60577324664060345999319382737747,
        -0.11253524348366090030649768852718,
        +0.13432266247902779492487859329414E-01,
        -0.19268451873811457249246838991303E-02,
        +0.30911833772060318335586737475368E-03,
        -0.53564132129618418776393559795147E-04,
        +0.98278128802474923952491882717237E-05,
        -0.18853689849165182826902891938910E-05,
        +0.37494319356894735406964042190531E-06,
        -0.76823455870552639273733465680556E-07,
        +0.16143270567198777552956300060868E-07,
        -0.34668022114907354566309060226027E-08,
        +0.75875420919036277572889747054114E-09,
        -0.16886433329881412573514526636703E-09,
        +0.38145706749552265682804250927272E-10,
        -0.87330266324446292706851718272334E-11,
        +0.20236728645867960961794311064330E-11,
        -0.47413283039555834655210340820160E-12,
        +0.11221172048389864324731799928920E-12,
        -0.26804225434840309912826809093395E-13,
        +0.64578514417716530343580369067212E-14,
        -0.15682760501666478830305702849194E-14,
        +0.38367865399315404861821516441408E-15,
        -0.94517173027579130478871048932556E-16,
        +0.23434812288949573293896666439133E-16,
        -0.58458661580214714576123194419882E-17,
        +0.14666229867947778605873617419195E-17,
        -0.36993923476444472706592538274474E-18,
        +0.93790159936721242136014291817813E-19,
        -0.23893673221937873136308224087381E-19,
        +0.61150624629497608051934223837866E-20,
        -0.15718585327554025507719853288106E-20,
        +0.40572387285585397769519294491306E-21,
        -0.10514026554738034990566367122773E-21,
        +0.27349664930638667785806003131733E-22,
        -0.71401604080205796099355574271999E-23,
        +0.18705552432235079986756924211199E-23,
        -0.49167468166870480520478020949333E-24,
        +0.12964988119684031730916087125333E-24,
        -0.34292515688362864461623940437333E-25,
        +0.90972241643887034329104820906666E-26,
        -0.24202112314316856489934847999999E-26,
        +0.64563612934639510757670475093333E-27,
        -0.17269132735340541122315987626666E-27,
        +0.46308611659151500715194231466666E-28,
        -0.12448703637214131241755170133333E-28,
        +0.33544574090520678532907007999999E-29,
        -0.90598868521070774437543935999999E-30,
        +0.24524147051474238587273216000000E-30,
        -0.66528178733552062817107967999999E-31
    };
    static ityp ae14cs[64] =
    {
        -0.1892918000753016825495679942820,
        -0.8648117855259871489968817056824E-01,
        +0.7224101543746594747021514839184E-02,
        -0.8097559457557386197159655610181E-03,
        +0.1099913443266138867179251157002E-03,
        -0.1717332998937767371495358814487E-04,
        +0.2985627514479283322825342495003E-05,
        -0.5659649145771930056560167267155E-06,
        +0.1152680839714140019226583501663E-06,
        -0.2495030440269338228842128765065E-07,
        +0.5692324201833754367039370368140E-08,
        -0.1359957664805600338490030939176E-08,
        +0.3384662888760884590184512925859E-09,
        -0.8737853904474681952350849316580E-10,
        +0.2331588663222659718612613400470E-10,
        -0.6411481049213785969753165196326E-11,
        +0.1812246980204816433384359484682E-11,
        -0.5253831761558460688819403840466E-12,
        +0.1559218272591925698855028609825E-12,
        -0.4729168297080398718476429369466E-13,
        +0.1463761864393243502076199493808E-13,
        -0.4617388988712924102232173623604E-14,
        +0.1482710348289369323789239660371E-14,
        -0.4841672496239229146973165734417E-15,
        +0.1606215575700290408116571966188E-15,
        -0.5408917538957170947895023784252E-16,
        +0.1847470159346897881370231402310E-16,
        -0.6395830792759094470500610425050E-17,
        +0.2242780721699759457250233276170E-17,
        -0.7961369173983947552744555308646E-18,
        +0.2859308111540197459808619929272E-18,
        -0.1038450244701137145900697137446E-18,
        +0.3812040607097975780866841008319E-19,
        -0.1413795417717200768717562723696E-19,
        +0.5295367865182740958305442594815E-20,
        -0.2002264245026825902137211131439E-20,
        +0.7640262751275196014736848610918E-21,
        -0.2941119006868787883311263523362E-21,
        +0.1141823539078927193037691483586E-21,
        -0.4469308475955298425247020718489E-22,
        +0.1763262410571750770630491408520E-22,
        -0.7009968187925902356351518262340E-23,
        +0.2807573556558378922287757507515E-23,
        -0.1132560944981086432141888891562E-23,
        +0.4600574684375017946156764233727E-24,
        -0.1881448598976133459864609148108E-24,
        +0.7744916111507730845444328478037E-25,
        -0.3208512760585368926702703826261E-25,
        +0.1337445542910839760619930421384E-25,
        -0.5608671881802217048894771735210E-26,
        +0.2365839716528537483710069473279E-26,
        -0.1003656195025305334065834526856E-26,
        +0.4281490878094161131286642556927E-27,
        -0.1836345261815318199691326958250E-27,
        +0.7917798231349540000097468678144E-28,
        -0.3431542358742220361025015775231E-28,
        +0.1494705493897103237475066008917E-28,
        -0.6542620279865705439739042420053E-29,
        +0.2877581395199171114340487353685E-29,
        -0.1271557211796024711027981200042E-29,
        +0.5644615555648722522388044622506E-30,
        -0.2516994994284095106080616830293E-30,
        +0.1127259818927510206370368804181E-30,
        -0.5069814875800460855562584719360E-31
    };
    static ityp e11cs[29] =
    {
        -0.16113461655571494025720663927566180E+02,
        +0.77940727787426802769272245891741497E+01,
        -0.19554058188631419507127283812814491E+01,
        +0.37337293866277945611517190865690209,
        -0.56925031910929019385263892220051166E-01,
        +0.72110777696600918537847724812635813E-02,
        -0.78104901449841593997715184089064148E-03,
        +0.73880933562621681878974881366177858E-04,
        -0.62028618758082045134358133607909712E-05,
        +0.46816002303176735524405823868362657E-06,
        -0.32092888533298649524072553027228719E-07,
        +0.20151997487404533394826262213019548E-08,
        -0.11673686816697793105356271695015419E-09,
        +0.62762706672039943397788748379615573E-11,
        -0.31481541672275441045246781802393600E-12,
        +0.14799041744493474210894472251733333E-13,
        -0.65457091583979673774263401588053333E-15,
        +0.27336872223137291142508012748799999E-16,
        -0.10813524349754406876721727624533333E-17,
        +0.40628328040434303295300348586666666E-19,
        -0.14535539358960455858914372266666666E-20,
        +0.49632746181648636830198442666666666E-22,
        -0.16208612696636044604866560000000000E-23,
        +0.50721448038607422226431999999999999E-25,
        -0.15235811133372207813973333333333333E-26,
        +0.44001511256103618696533333333333333E-28,
        -0.12236141945416231594666666666666666E-29,
        +0.32809216661066001066666666666666666E-31,
        -0.84933452268306432000000000000000000E-33
    };
    static ityp e12cs[25] =
    {
        -0.3739021479220279511668698204827E-01,
        +0.4272398606220957726049179176528E-01,
        -0.130318207984970054415392055219726,
        +0.144191240246988907341095893982137E-01,
        -0.134617078051068022116121527983553E-02,
        +0.107310292530637799976115850970073E-03,
        -0.742999951611943649610283062223163E-05,
        +0.453773256907537139386383211511827E-06,
        -0.247641721139060131846547423802912E-07,
        +0.122076581374590953700228167846102E-08,
        -0.548514148064092393821357398028261E-10,
        +0.226362142130078799293688162377002E-11,
        -0.863589727169800979404172916282240E-13,
        +0.306291553669332997581032894881279E-14,
        -0.101485718855944147557128906734933E-15,
        +0.315482174034069877546855328426666E-17,
        -0.923604240769240954484015923200000E-19,
        +0.255504267970814002440435029333333E-20,
        -0.669912805684566847217882453333333E-22,
        +0.166925405435387319431987199999999E-23,
        -0.396254925184379641856000000000000E-25,
        +0.898135896598511332010666666666666E-27,
        -0.194763366993016433322666666666666E-28,
        +0.404836019024630033066666666666666E-30,
        -0.807981567699845120000000000000000E-32
    };
    ityp eta;
    static dim_typ ntae10 = 0;
    static dim_typ ntae11 = 0;
    static dim_typ ntae12 = 0;
    static dim_typ ntae13 = 0;
    static dim_typ ntae14 = 0;
    static dim_typ nte11 = 0;
    static dim_typ nte12 = 0;
    ityp value;
    static ityp xmax = 0.00;

    if ( ntae10 == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        ntae10 = r8_inits ( ae10cs, 50, eta );
        ntae11 = r8_inits ( ae11cs, 60, eta );
        ntae12 = r8_inits ( ae12cs, 41, eta );
        nte11 = r8_inits ( e11cs, 29, eta );
        nte12 = r8_inits ( e12cs, 25, eta );
        ntae13 = r8_inits ( ae13cs, 50, eta );
        ntae14 = r8_inits ( ae14cs, 64, eta );
        xmax = - log ( r8_mach ( 1 ) );
        xmax -= log ( xmax );
    }

    if ( x <= - 32.00 )
        value = exp ( - x ) / x * ( 1.00+ r8_csevl ( 64.00 / x + 1.00, ae10cs, ntae10 ) );
    else if ( x <= - 8.00 )
        value = exp ( - x ) / x * ( 1.00+ r8_csevl ( ( 64.00 / x + 5.00 ) / 3.00, ae11cs, ntae11 ) );
    else if ( x <= - 4.00 )
        value = exp ( - x ) / x * (1.00+ r8_csevl ( 16.00 / x + 3.00, ae12cs, ntae12 ) );
    else if ( x <= - 1.00 )
        value = - log ( - x )+ r8_csevl ( ( 2.0 * x + 5.0 ) / 3.0, e11cs, nte11 );
    else if ( x == 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }
    else if ( x <= 1.00 )
        value = ( - log ( abs ( x ) ) - 0.6875 + x )+ r8_csevl ( x, e12cs, nte12 );
    else if ( x <= 4.00 )
        value = exp ( - x ) / x * ( 1.00+ r8_csevl ( ( 8.00 / x - 5.00 ) / 3.00, ae13cs, ntae13 ) );
    else if ( x <= xmax )
        value = exp ( - x ) / x * ( 1.00+ r8_csevl ( 8.00 / x - 1.00, ae14cs, ntae14 ) );
    else
    	value = 0.00;

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_ei ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_EI evaluates the exponential integral Ei for an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_EI, the exponential integral Ei evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	result = - r8_e1 ( - x );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_erf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_ERF evaluates the error function of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_ERF, the error function of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp erfcs[21] =
    {
        -0.49046121234691808039984544033376E-01,
        -0.14226120510371364237824741899631,
        +0.10035582187599795575754676712933E-01,
        -0.57687646997674847650827025509167E-03,
        +0.27419931252196061034422160791471E-04,
        -0.11043175507344507604135381295905E-05,
        +0.38488755420345036949961311498174E-07,
        -0.11808582533875466969631751801581E-08,
        +0.32334215826050909646402930953354E-10,
        -0.79910159470045487581607374708595E-12,
        +0.17990725113961455611967245486634E-13,
        -0.37186354878186926382316828209493E-15,
        +0.71035990037142529711689908394666E-17,
        -0.12612455119155225832495424853333E-18,
        +0.20916406941769294369170500266666E-20,
        -0.32539731029314072982364160000000E-22,
        +0.47668672097976748332373333333333E-24,
        -0.65980120782851343155199999999999E-26,
        +0.86550114699637626197333333333333E-28,
        -0.10788925177498064213333333333333E-29,
        +0.12811883993017002666666666666666E-31
    };
    static dim_typ nterf = 0;
    static ityp sqeps = 0.00;
    static ityp sqrtpi = 1.77245385090551602729816748334115;
    ityp value;
    static ityp xbig = 0.00;
    ityp y;

    if ( nterf == 0 )
    {
        nterf = r8_inits ( erfcs, 21, 0.10 * r8_mach ( 3 ) );
        xbig = sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) );
        sqeps = sqrt ( 2.00 * r8_mach ( 3 ) );
    }

    y = abs ( x );

    if ( y <= sqeps )
        value = 2.00 * x / sqrtpi;
    else if ( y <= 1.00 )
        value = x * ( 1.00 + r8_csevl ( 2.00 * x * x - 1.00, erfcs, nterf ) );
    else if ( y <= xbig )
    {
        value = 1.00 - r8_erfc ( y );
        if ( x < 0.00 )
            value *= -1;
    }
    else
    {
        value = 1.00;
        if ( x < 0.0 )
            value *= -1;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_erfc ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_ERFC evaluates the co-error function of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_ERFC, the co-error function of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp erc2cs[49] =
    {
        -0.6960134660230950112739150826197E-01,
        -0.4110133936262089348982212084666E-01,
        +0.3914495866689626881561143705244E-02,
        -0.4906395650548979161280935450774E-03,
        +0.7157479001377036380760894141825E-04,
        -0.1153071634131232833808232847912E-04,
        +0.1994670590201997635052314867709E-05,
        -0.3642666471599222873936118430711E-06,
        +0.6944372610005012589931277214633E-07,
        -0.1371220902104366019534605141210E-07,
        +0.2788389661007137131963860348087E-08,
        -0.5814164724331161551864791050316E-09,
        +0.1238920491752753181180168817950E-09,
        -0.2690639145306743432390424937889E-10,
        +0.5942614350847910982444709683840E-11,
        -0.1332386735758119579287754420570E-11,
        +0.3028046806177132017173697243304E-12,
        -0.6966648814941032588795867588954E-13,
        +0.1620854541053922969812893227628E-13,
        -0.3809934465250491999876913057729E-14,
        +0.9040487815978831149368971012975E-15,
        -0.2164006195089607347809812047003E-15,
        +0.5222102233995854984607980244172E-16,
        -0.1269729602364555336372415527780E-16,
        +0.3109145504276197583836227412951E-17,
        -0.7663762920320385524009566714811E-18,
        +0.1900819251362745202536929733290E-18,
        -0.4742207279069039545225655999965E-19,
        +0.1189649200076528382880683078451E-19,
        -0.3000035590325780256845271313066E-20,
        +0.7602993453043246173019385277098E-21,
        -0.1935909447606872881569811049130E-21,
        +0.4951399124773337881000042386773E-22,
        -0.1271807481336371879608621989888E-22,
        +0.3280049600469513043315841652053E-23,
        -0.8492320176822896568924792422399E-24,
        +0.2206917892807560223519879987199E-24,
        -0.5755617245696528498312819507199E-25,
        +0.1506191533639234250354144051199E-25,
        -0.3954502959018796953104285695999E-26,
        +0.1041529704151500979984645051733E-26,
        -0.2751487795278765079450178901333E-27,
        +0.7290058205497557408997703680000E-28,
        -0.1936939645915947804077501098666E-28,
        +0.5160357112051487298370054826666E-29,
        -0.1378419322193094099389644800000E-29,
        +0.3691326793107069042251093333333E-30,
        -0.9909389590624365420653226666666E-31,
        +0.2666491705195388413323946666666E-31
    };
    static ityp erfccs[59] =
    {
        +0.715179310202924774503697709496E-01,
        -0.265324343376067157558893386681E-01,
        +0.171115397792085588332699194606E-02,
        -0.163751663458517884163746404749E-03,
        +0.198712935005520364995974806758E-04,
        -0.284371241276655508750175183152E-05,
        +0.460616130896313036969379968464E-06,
        -0.822775302587920842057766536366E-07,
        +0.159214187277090112989358340826E-07,
        -0.329507136225284321486631665072E-08,
        +0.722343976040055546581261153890E-09,
        -0.166485581339872959344695966886E-09,
        +0.401039258823766482077671768814E-10,
        -0.100481621442573113272170176283E-10,
        +0.260827591330033380859341009439E-11,
        -0.699111056040402486557697812476E-12,
        +0.192949233326170708624205749803E-12,
        -0.547013118875433106490125085271E-13,
        +0.158966330976269744839084032762E-13,
        -0.472689398019755483920369584290E-14,
        +0.143587337678498478672873997840E-14,
        -0.444951056181735839417250062829E-15,
        +0.140481088476823343737305537466E-15,
        -0.451381838776421089625963281623E-16,
        +0.147452154104513307787018713262E-16,
        -0.489262140694577615436841552532E-17,
        +0.164761214141064673895301522827E-17,
        -0.562681717632940809299928521323E-18,
        +0.194744338223207851429197867821E-18,
        -0.682630564294842072956664144723E-19,
        +0.242198888729864924018301125438E-19,
        -0.869341413350307042563800861857E-20,
        +0.315518034622808557122363401262E-20,
        -0.115737232404960874261239486742E-20,
        +0.428894716160565394623737097442E-21,
        -0.160503074205761685005737770964E-21,
        +0.606329875745380264495069923027E-22,
        -0.231140425169795849098840801367E-22,
        +0.888877854066188552554702955697E-23,
        -0.344726057665137652230718495566E-23,
        +0.134786546020696506827582774181E-23,
        -0.531179407112502173645873201807E-24,
        +0.210934105861978316828954734537E-24,
        -0.843836558792378911598133256738E-25,
        +0.339998252494520890627359576337E-25,
        -0.137945238807324209002238377110E-25,
        +0.563449031183325261513392634811E-26,
        -0.231649043447706544823427752700E-26,
        +0.958446284460181015263158381226E-27,
        -0.399072288033010972624224850193E-27,
        +0.167212922594447736017228709669E-27,
        -0.704599152276601385638803782587E-28,
        +0.297976840286420635412357989444E-28,
        -0.126252246646061929722422632994E-28,
        +0.539543870454248793985299653154E-29,
        -0.238099288253145918675346190062E-29,
        +0.109905283010276157359726683750E-29,
        -0.486771374164496572732518677435E-30,
        +0.152587726411035756763200828211E-30
    };
    static ityp erfcs[21] =
    {
        -0.49046121234691808039984544033376E-01,
        -0.14226120510371364237824741899631,
        +0.10035582187599795575754676712933E-01,
        -0.57687646997674847650827025509167E-03,
        +0.27419931252196061034422160791471E-04,
        -0.11043175507344507604135381295905E-05,
        +0.38488755420345036949961311498174E-07,
        -0.11808582533875466969631751801581E-08,
        +0.32334215826050909646402930953354E-10,
        -0.79910159470045487581607374708595E-12,
        +0.17990725113961455611967245486634E-13,
        -0.37186354878186926382316828209493E-15,
        +0.71035990037142529711689908394666E-17,
        -0.12612455119155225832495424853333E-18,
        +0.20916406941769294369170500266666E-20,
        -0.32539731029314072982364160000000E-22,
        +0.47668672097976748332373333333333E-24,
        -0.65980120782851343155199999999999E-26,
        +0.86550114699637626197333333333333E-28,
        -0.10788925177498064213333333333333E-29,
        +0.12811883993017002666666666666666E-31
    };
    ityp eta;
    static int nterc2 = 0;
    static int nterf = 0;
    static int nterfc = 0;
    static ityp sqeps = 0.00;
    static ityp sqrtpi = 1.77245385090551602729816748334115;
    ityp value;
    static ityp xmax = 0.00;
    static ityp xsml = 0.00;
    ityp y;

    if ( nterf == 0 )
    {
        eta = 0.10 * r8_mach ( 3 );
        nterf = r8_inits ( erfcs, 21, eta );
        nterfc = r8_inits ( erfccs, 59, eta );
        nterc2 = r8_inits ( erc2cs, 49, eta );

        xsml = - sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) );
        xmax = sqrt (- log ( sqrtpi * r8_mach ( 1 ) ) );
        xmax = xmax - 0.50 * log ( xmax ) / xmax - 0.01;
        sqeps = sqrt ( 2.00 * r8_mach ( 3 ) );
    }

    if ( x <= xsml )
    {
    	result = 2.00;
        return &result;
    }

    if ( xmax < x )
    {
    	result = MAX_VAL;
        return &result;
    }

    y = abs ( x );

    if ( y < sqeps )
    {
    	result = 2.00 * x / sqrtpi;
        return &result;
    }
    else if ( y <= 1.0 )
    {
    	result = 1.00 - x * ( 1.00+ r8_csevl ( 2.00 * x * x - 1.00, erfcs, nterf ) );
        return &result;
	}

    y = y * y;
    value = y<=4.00 ? exp ( - y ) / abs ( x ) * ( 0.50+ r8_csevl ( ( 8.00 / y - 5.00 ) / 3.00, erc2cs, nterc2 ) ) : exp ( - y ) / abs ( x ) * ( 0.5+ r8_csevl ( 8.0 / y - 1.0, erfccs, nterfc ) );

    if ( x < 0.00 )
        value = 2.00 - value;

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_exprel ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_EXPREL evaluates the exponential relative error term of an R8 argument.
  Discussion:
    The relative error term is ( exp ( x ) - 1 ) / x.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_EXPREL, the exponential relative error term
    at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp absx;
    ityp alneps;
    dim_typ i;
    static dim_typ nterms = 0;
    ityp value;
    static ityp xbnd = 0.00;
    ityp xln;
    ityp xn;

    if ( nterms == 0 )
    {
        alneps = log ( r8_mach ( 3 ) );
        xn = 3.72 - 0.30 * alneps;
        xln = log ( ( xn + 1.00 ) / 1.36 );
        nterms = ( dim_typ ) ( xn - ( xn * xln + alneps ) / ( xln + 1.36 ) + 1.5 );
        xbnd = r8_mach ( 3 );
    }

    absx = abs ( x );

    if ( absx < xbnd )
    value = 1.00;
    else if ( absx <= 0.50 )
    {
        value = 0.00;
        for ( i = 1; i <= nterms; ++i )
            value = 1.00 + value * x / ( ityp ) ( nterms + 2 - i );
    }
    else
        value = ( exp ( x ) - 1.00 ) / x;

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_fac ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_FAC evaluates the factorial of an I4 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, int N, the argument.
    Output, double R8_FAC, the factorial of N.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
    static ityp facn[31] =
    {
        +0.100000000000000000000000000000000E+01,
        +0.100000000000000000000000000000000E+01,
        +0.200000000000000000000000000000000E+01,
        +0.600000000000000000000000000000000E+01,
        +0.240000000000000000000000000000000E+02,
        +0.120000000000000000000000000000000E+03,
        +0.720000000000000000000000000000000E+03,
        +0.504000000000000000000000000000000E+04,
        +0.403200000000000000000000000000000E+05,
        +0.362880000000000000000000000000000E+06,
        +0.362880000000000000000000000000000E+07,
        +0.399168000000000000000000000000000E+08,
        +0.479001600000000000000000000000000E+09,
        +0.622702080000000000000000000000000E+10,
        +0.871782912000000000000000000000000E+11,
        +0.130767436800000000000000000000000E+13,
        +0.209227898880000000000000000000000E+14,
        +0.355687428096000000000000000000000E+15,
        +0.640237370572800000000000000000000E+16,
        +0.121645100408832000000000000000000E+18,
        +0.243290200817664000000000000000000E+19,
        +0.510909421717094400000000000000000E+20,
        +0.112400072777760768000000000000000E+22,
        +0.258520167388849766400000000000000E+23,
        +0.620448401733239439360000000000000E+24,
        +0.155112100433309859840000000000000E+26,
        +0.403291461126605635584000000000000E+27,
        +0.108888694504183521607680000000000E+29,
        +0.304888344611713860501504000000000E+30,
        +0.884176199373970195454361600000000E+31,
        +0.265252859812191058636308480000000E+33
    };
    static dim_typ nmax = 0;
    static ityp sq2pil = 0.91893853320467274178032973640562;
    ityp value;
    ityp x;
    ityp xmax;
    ityp xmin;

    if ( nmax == 0 )
    {
        r8_gaml ( &xmin, &xmax );
        nmax = ( dim_typ ) ( xmax - 1.00 );
    }

   	if ( n <= 30 )
        value = facn[n];
    else if ( n <= nmax )
    {
        x = ( ityp ) ( n + 1 );
        value = exp ( ( x - 0.50 ) * log ( x ) - x + sq2pil + r8_lgmc ( x ) );
    }
    else
    {
    	result = MAX_VAL;
        return &result;
    }

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gami ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAMI evaluates the incomplete gamma function for an R8 argument.
  Discussion:
    GAMI = Integral ( 0 <= T <= X ) exp ( - t ) * t^( a - 1 )  dt
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, the parameter.
    Input, double X, the argument.
    Output, double R8_GAMI, the value of the incomplete gamma function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp x = a_data[1];
	
    ityp factor;
    ityp value;

    if ( a <= 0.00 || x < 0.00)
    {
    	result = MAX_VAL;
        return &result;
    }
    else if ( x == 0.00 )
    value = 0.00;
    else
    {
        factor = exp ( r8_lngam ( a ) + a * log ( x ) );
        value = factor * r8_gamit ( a, x );
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gamic ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAMIC evaluates the complementary incomplete gamma function.
  Discussion:
    GAMIC = integral ( x <= t < oo ) exp(-t) * t^(a-1) dt
    GAMIC is evaluated for arbitrary real values of A and non-negative
    values X (even though GAMIC is defined for X < 0.0), except that
    for X = 0 and A <= 0.0, GAMIC is undefined.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
    Walter Gautschi,
    A Computational Procedure for Incomplete Gamma Functions,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 4, December 1979, pages 466-481.
  Parameters:
    Input, double A, the parameter.
    Input, double X, the evaluation point.
    Output, double R8_GAMIC, the value of the incomplete
    gamma function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp x = a_data[1];
	
    ityp aeps;
    ityp ainta;
    ityp algap1;
    static ityp alneps = 0.00;
    ityp alngs;
    ityp alx;
    static ityp bot = 0.00;
    ityp e;
    static ityp eps = 0.00;
    ityp gstar;
    ityp h;
    dim_typ izero;
    dim_typ ma;
    ityp sga;
    ityp sgng;
    ityp sgngam;
    ityp sgngs;
    static ityp sqeps = 0.00;
    ityp t;
    ityp value;

    if ( eps == 0.00 )
    {
        eps = 0.50 * r8_mach ( 3 );
        sqeps = sqrt ( r8_mach ( 4 ) );
        alneps = - log ( r8_mach ( 3 ) );
        bot = log ( r8_mach ( 1 ) );
    }

    if ( x < 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( x == 0.00 )
    {
        if ( a <= 0.00 )
        {
        	result = MAX_VAL;
        	return &result;
        }
        
        result = exp ( r8_lngam ( a + 1.00 ) - log ( a ) );
        return &result;
    }

    alx = log ( x );
    sga = 1.00 - ((a<0.00)<<1);
    ainta = r8_aint ( a + 0.50 * sga );
    aeps = a - ainta;

    izero = 0;

    if ( x < 1.00 )
    {
        if ( a <= 0.50 && abs ( aeps ) <= 0.001 )
        {
            e = -ainta<=1.00 ? 2.00 : 2.00 * ( - ainta + 2.00 ) / ( ainta * ainta - 1.00 );
            e -= alx * pow ( x, - 0.001 );
            if ( e * abs ( aeps ) <= eps )
            {
            	result = r8_gmic ( a, x, alx );
                return &result;
            }
        }

        r8_lgams ( a + 1.0, &algap1, &sgngam );
        gstar = r8_gmit ( a, x, algap1, sgngam, alx );

        if ( gstar == 0.00 )
            izero = 1;
        else
        {
            alngs = log ( abs ( gstar ) );
            sgngs = r8_sign ( gstar );
        }
    }
    else
    {
        if ( a < x )
        {
        	result =exp ( r8_lgic ( a, x, alx ) );
        	return &result;
    	}

        sgngam = 1.00;
        algap1 = r8_lngam ( a + 1.0 );
        sgngs = 1.00;
        alngs = r8_lgit ( a, x, algap1 );
    }

    h = 1.00;

    if ( izero != 1 )
    {
        t = a * alx + alngs;

        if ( alneps < t )
        {
            sgng = - sgngs * sga * sgngam;
            result = sgng * exp ( t + algap1 - log ( abs ( a ) ) );
            return &result;
        }

        if ( - alneps < t )
            h = 1.00 - sgngs * exp ( t );
    }
    
    result = r8_sign ( h ) * sga * sgngam * exp ( log ( abs ( h ) ) + algap1 - log ( abs ( a ) ));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gamit ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAMIT evaluates Tricomi's incomplete gamma function for an R8 argument.
  Discussion:
      GAMIT = x^(-a) / gamma(a)
        * Integral ( 0 <= t <= x ) exp(-t) * t^(a-1) dt
    with analytic continuation for a <= 0.0.  Gamma(x) is the complete
    gamma function of X.  GAMIT is evaluated for arbitrary real values of
    A and for non-negative values of X (even though GAMIT is defined for
    X < 0.0).
    A slight deterioration of 2 or 3 digits accuracy will occur when
    gamit is very large or very small in absolute value, because log-
    arithmic variables are used.  Also, if the parameter A is very close
    to a negative integer (but not a negative integer), there is a loss
    of accuracy.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
    Walter Gautschi,
    A Computational Procedure for Incomplete Gamma Functions,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 4, December 1979, pages 466-481.
  Parameters:
    Input, double A, the parameter.
    Input, double X, the argument.
    Output, double R8_GAMIT, the function value.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp x = a_data[1];
	
    ityp aeps;
    ityp ainta;
    ityp algap1;
    static ityp alneps = 0.00;
    ityp alng;
    ityp alx;
    static ityp bot = 0.00;
    ityp h;
    ityp sga;
    ityp sgngam;
    static ityp sqeps = 0.00;
    ityp t;
    ityp value;

    if ( alneps == 0.00 )
    {
        alneps = - log ( r8_mach ( 3 ) );
        sqeps = sqrt ( r8_mach ( 4 ) );
        bot = log ( r8_mach ( 1 ) );
    }

    if ( x < 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    else
        alx = 0.00 + (x!=0)*log(x);

    sga = 1.00 - ((a<0.00)<<1);
    ainta = r8_aint ( a + 0.50 * sga );
    aeps = a - ainta;

    if ( x == 0.00 )
    {
    	result = 0.00 + (0.00 < ainta || aeps != 0.00)*r8_gamr ( a + 1.00 );
        return &result;
    }

    if ( x <= 1.00 )
    {
        if ( - 0.50 <= a || aeps != 0.00 )
            r8_lgams ( a + 1.00, &algap1, &sgngam );
            
        result = r8_gmit ( a, x, algap1, sgngam, alx );
        return &result;
    }

    if ( x <= a )
    {
    	result = exp ( r8_lgit (a, x, r8_lngam ( a + 1.00 ) ) );
        return &result;
    }

    alng = r8_lgic ( a, x, alx );
    /*
    Evaluate in terms of log (r8_gamic (a, x))
    */
    h = 1.00;

    if ( aeps != 0.00 || 0.00 < ainta )
    {
        r8_lgams ( a + 1.00, &algap1, &sgngam );
        t = log ( abs ( a ) ) + alng - algap1;

        if ( alneps < t )
        {
        	result = - sga * sgngam * exp ( t - a * alx );
            return &result;
        }

        if ( - alneps < t )
            h = 1.00 - sga * sgngam * exp ( t );
    }
    t = - a * alx + log ( abs ( h ) );
    
    result =h<0.00 ? -exp(t):exp(t);
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gamr ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAMR evaluates the reciprocal gamma function of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_GAMR, the value of the reciprocal gamma
    function at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp alngx;
    ityp value;
    ityp sgngx;

    if ( x <= 0.00 && r8_aint ( x ) == x )
        value = 0.00;
    else if ( abs ( x ) <= 10.00 )
        value = 1.00 / r8_gamma ( x );
    else
    {
        r8_lgams ( x, &alngx, &sgngx );
        value = sgngx * exp ( - alngx );
    }
    
    result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gmic ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GMIC: complementary incomplete gamma, small X, A near negative int.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, the parameter.
    Input, double X, the argument.
    Input, double ALX, the logarithm of X.
    Output, double R8_GMIC, the complementary incomplete
    gamma function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp x = a_data[1];
	const register ityp alx = a_data[2];
	
    ityp alng;
    static ityp bot = 0.0;
    bool converged;
    static ityp eps = 0.0;
    static ityp euler = 0.57721566490153286060651209008240;
    ityp fk;
    ityp fkp1;
    ityp fm;
    dim_typ k;
    dim_typ m;
    dim_typ ma;
    dim_typ mm1;
    ityp s;
    ityp sgng;
    ityp t;
    ityp te;
    ityp value;

    if ( eps == 0.00 )
    {
        eps = 0.50 * r8_mach ( 3 );
        bot = log ( r8_mach ( 1 ) );
    }

    if ( 0.00 < a || x <= 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    m = - ( dim_typ ) ( a - 0.50 );
    fm = ( ityp ) ( m );

    te = t = 1.00;
    s = t;
    converged = false;

    for ( k = 1; k <= 200; ++k )
    {
        fkp1 = ( ityp ) ( k + 1 );
        te = - x * te / ( fm + fkp1 );
        t = te / fkp1;
        s += t;
        if ( abs ( t ) < eps * s )
        {
            converged = true;
            break;
        }
    }

    if ( !converged )
    {
    	result = MAX_VAL;
        return &result;
    }

    value = - alx - euler + x * s / ( fm + 1.00 );

    if ( m == 0 )
    {
    	result = value;
        return &result;
    }
    else if ( m == 1 )
    {
    	result = - value - 1.00 + 1.00 / x;
        return &result;
    }

    te = fm;
    t = 1.00;
    s = t;
    mm1 = m - 1;

    for ( k = 1; k <= mm1; ++k)
    {
        fk = ( ityp ) ( k );
        te = - x * te / fk;
        t = te / ( fm - fk );
        s += t;
        if ( abs ( t ) < eps * abs ( s ) )
            break;
    }

    for ( k = 1; k <= m; ++k )
        value += 1.00 / ( ityp ) ( k );

    sgng = 1.00 - (((m%2)==1)<<1);
    alng = log ( value ) - r8_lngam ( fm + 1.00 );

    value = 0.00 + (bot<alng)*sgng*exp(alng);

    if ( s != 0.00 )
        value += r8_sign ( s ) * exp ( - fm * alx + log ( abs ( s ) / fm ) );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gmit ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GMIT: Tricomi's incomplete gamma function for small X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, the parameter.
    Input, double X, the argument.
    Input, double ALGAP1, the logarithm of Gamma ( A + 1 ).
    Input, double SGNGAM, the sign of Gamma ( A + 1 ).
    Input, double ALX, the logarithm of X.
    Output, double R8_GMIT, the Tricomi incomplete gamma function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	ityp a = a_data[0];
	ityp x = a_data[1];
	ityp algap1 = a_data[2];
	ityp sgngam = a_data[3];
	ityp alx = a_data[4];
	
    ityp ae;
    ityp aeps;
    ityp alg2;
    ityp algs;
    static ityp bot = 0.00;
    bool converged;
    static ityp eps = 0.00;
    ityp fk;
    dim_typ k;
    dim_typ m;
    dim_typ ma;
    ityp s;
    ityp sgng2;
    ityp t;
    ityp te;
    ityp value;

    if ( eps == 0.00 )
    {
        eps = 0.50 * r8_mach ( 3 );
        bot = log ( r8_mach ( 1 ) );
    }

    if ( x <= 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    ma = a<0.00 ? ( dim_typ ) ( a - 0.50 ) : ( dim_typ ) ( a + 0.50 );
    aeps = a - ( ityp ) ( ma );
    ae = a<-0.5 ? aeps:a;
    t = 1.00;
    te = ae;
    s = t;
    converged = false;

    for ( k = 1; k <= 200; ++k )
    {
        fk = ( ityp ) ( k );
        te = - x * te / fk;
        t = te / ( ae + fk );
        s += t;
        if ( abs ( t ) < eps * abs ( s ) )
        {
            converged = true;
            break;
        }
    }

    if ( !converged )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( - 0.50 <= a )
    {
    	result = exp ( - algap1 + log ( s ) );
        return &result;
    }

    algs = - r8_lngam ( 1.00 + aeps ) + log ( s );
    s = t = 1.00;
    m = - ma - 1;

    for ( k = 1; k <= m; ++k )
    {
        t = x * t / ( aeps - ( ityp ) ( m + 1 - k ) );
        s += t;
        if ( abs ( t ) < eps * abs ( s ) )
            break;
    }

    value = 0.00;
    algs = - ( ityp ) ( ma ) * log ( x ) + algs;

    if ( s == 0.00 || aeps == 0.00 )
    {
    	result = exp ( algs );
        return &result;
    }

    sgng2 = sgngam * r8_sign ( s );
    alg2 = - x - algap1 + log ( abs ( s ) );

    if ( bot < alg2 )
        value = sgng2 * exp ( alg2 );

    if ( bot < algs )
    	value += exp ( algs );

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_int ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_INT returns the integer part of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_INT, the integer part of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    dim_typ i;
    dim_typ ibase;
    dim_typ ipart;
    static dim_typ npart = 0;
    ityp part;
    static ityp scale = 0.00;
    ityp value;
    static ityp xbig = 0.00;
    static ityp xmax = 0.00;
    ityp xscl;

    if ( npart == 0 )
    {
        ibase = i4_mach ( 10 );
        xmax = 1.00 / r8_mach ( 4 );
        xbig = MIN ( ( ityp ) ( i4_mach ( 9 ) ), 1.00 / r8_mach ( 4 ) );
        scale = ( ityp ) i4_pow ( ibase, ( dim_typ ) ( log ( xbig ) / log ( ( ityp ) ( ibase ) ) - 0.50 ) );
        npart = log ( xmax ) / log ( scale ) + 1.00;
    }
    /*
    X may be too small.
    */
    if ( x < - xmax )
        value = x;
    else if ( x < - xbig )
    {
        xscl = - x;

        for ( i = 1; i <= npart; ++i )
            xscl /= scale;

        value = 0.00;
        for ( i = 1; i <= npart; i++ )
        {
            xscl *= scale;
            ipart = ( dim_typ ) ( xscl );
            part = ( ityp ) ( ipart );
            xscl -= part;
            value *= scale + part;
        }
        value *= -1;
    }
    else if ( x <= xbig )
        value = ( dim_typ ) ( x );
    else if ( x <= xmax )
    {
        xscl = x;

        for ( i = 1; i <= npart; ++i )
            xscl /= scale;

        value = 0.00;
        for ( i = 1; i <= npart; ++i )
        {
            xscl *= scale;
            ipart = ( dim_typ ) ( xscl );
            part = ( ityp ) ( ipart );
            xscl -= part;
            value *= scale + part;
        }
    }
    /*
    X may be too large.
    */
    else
        value = x;

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_lbeta ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_LBETA evaluates the logarithm of the beta function of R8 arguments.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, B, the arguments.
    Output, double R8_LBETA, the logarithm of the beta function of A
    and B.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp b = a_data[1];
	
    ityp corr;
    ityp p;
    ityp q;
    static ityp sq2pil = 0.91893853320467274178032973640562;
    ityp value;

    p = MIN ( a, b );
    q = MAX ( a, b );

    if ( p <= 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }
    else if ( p < 10.00 && q <= 10.00 )
        value = log ( r8_gamma ( p ) * ( r8_gamma ( q ) / r8_gamma ( p + q ) ) );
    else if ( p < 10.00 )
    {
        corr = r8_lgmc ( q ) - r8_lgmc ( p + q );
        value = r8_lngam ( p ) + corr + p - p * log ( p + q ) +( q - 0.5 ) * r8_lnrel ( - p / ( p + q ) );
    }
    else
    {
        corr = r8_lgmc ( p ) + r8_lgmc ( q ) - r8_lgmc ( p + q );
        value = - 0.50 * log ( q ) + sq2pil + corr+ ( p - 0.50 ) * log ( p / ( p + q ) )+ q * r8_lnrel ( - p / ( p + q ) );
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_lgams ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_LGAMS evaluates the log of |gamma(x)| and sign, for an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double *ALGAM, the logarithm of the absolute value of
    gamma ( X ).
    Output, double *SGNGAM, the sign (+1 or -1) of gamma ( X ).
*/
{
	const it2pit * const s_data = data;
	const register ityp x  = s_data->a0;
	ityp * algam = s_data->a1;
	ityp * sgngam = s_data->a2;
	
    dim_typ k;

    *algam = r8_lngam ( x );
    *sgngam = 1.00;

    if ( x <= 0.00 )
    {
        k = ( int ) ( fmod ( - r8_aint ( x ), 2.00 ) + 0.10 );

        if ( k == 0 )
            *sgngam = - 1.00;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_lgic ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_LGIC evaluates the log complementary incomplete gamma function for large X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, the parameter.
    Input, double X, the argument.
    Input, double ALX, the logarithm of X.
    Output, double R8_LGIC, the log complementary incomplete
    gamma function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp x = a_data[1];
	const register ityp alx = a_data[2];
	
    static ityp eps = 0.00;
    ityp fk;
    dim_typ k;
    ityp p;
    ityp r;
    ityp s;
    ityp t;
    ityp value;
    ityp xma;
    ityp xpa;

    if ( eps == 0.00 )
        eps = 0.50 * r8_mach ( 3 );

    xpa = x + 1.00 - a;
    xma = x - 1.00 - a;

    r = 0.00;
    p = 1.00;
    s = p;
    for ( k = 1; k <= 300; ++k )
    {
        fk = ( ityp ) ( k );
        t = fk * ( a - fk ) * ( 1.00 + r );
        r = - t / ( ( xma + 2.00 * fk ) * ( xpa + 2.00 * fk ) + t );
        p = r * p;
        s += p;
        if ( abs ( p ) < eps * s )
        {
        	result = a * alx - x + log ( s / xpa );
            return &result;
        }
    }
    
    result = MAX_VAL;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_lgit ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_LGIT evaluates the log of Tricomi's incomplete gamma function.
  Discussion:
    Perron's continued fraction is used for large X and X <= A.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, the parameter.
    Input, double X, the argument.
    Input, double ALGAP1, the logarithm of A+1.
    Output, double R8_LGIT, the log of Tricomi's incomplete
    gamma function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp x = a_data[1];
	const register ityp algap1 = a_data[2];
	
    ityp a1x;
    ityp ax;
    static ityp eps = 0.00;
    ityp fk;
    ityp hstar;
    dim_typ k;
    ityp p;
    ityp r;
    ityp s;
    static ityp sqeps = 0.00;
    ityp t;
    ityp value;

    if ( eps == 0.00 )
    {
        eps = 0.50 * r8_mach ( 3 );
        sqeps = sqrt ( r8_mach ( 4 ) );
    }

    if ( x <= 0.0 || a < x )
    {
    	result = MAX_VAL;
        return &result;
    }

    ax = a + x;
    a1x = ax + 1.00;
    r = 0.00;
    s = p = 1.00;

    for ( k = 1; k <= 200; ++k )
    {
        fk = ( ityp ) ( k );
        t = ( a + fk ) * x * ( 1.00 + r );
        r = t / ( ( ax + fk ) * ( a1x + fk ) - t );
        p = r * p;
        s += p;
        if ( abs ( p ) < eps * s )
        {
        	result = - x - algap1 - log ( 1.00 - x * s / a1x );
            return &result;
        }
    }

    result = MAX_VAL;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_li ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_LI evaluates the logarithmic integral for an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_LI, the logarithmic integral evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp sqeps = 0.00;
    ityp value;
    if ( sqeps == 0.00 )
        sqeps = sqrt ( r8_mach ( 3 ) );
    if ( x < 0.0 || x == 1.00 )
    {
    	result = MAX_VAL;
        return &result;
    }
    if ( x == 0.00 )
    {
    	result = 0.00;
        return &result;
    }
    
    result = r8_ei ( log ( x ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_lngam ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_LNGAM: log of the absolute value of gamma of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of N:umerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_LNGAM, the logarithm of the absolute value of
    the gamma function of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp dxrel = 0.00;
    ityp sinpiy;
    static ityp sq2pil = 0.91893853320467274178032973640562;
    ityp value;
    static ityp xmax = 0.00;
    ityp y;

    if ( xmax == 0.00 )
    {
        xmax = r8_mach ( 2 ) / log ( r8_mach ( 2 ) );
        dxrel = sqrt ( r8_mach ( 4 ) );
    }

    y = abs ( x );

    if ( y <= 10.00 )
    {
    	result = log ( abs ( r8_gamma ( x ) ) );
        return &result;
    }

    if ( xmax < y )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( 0.00 < x )
    {
    	result = sq2pil + ( x - 0.50 ) * log ( x ) - x + r8_lgmc ( y );
        return &result;
    }

    sinpiy = abs ( sin ( M_PI * y ) );

    if ( sinpiy == 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

	result = M_SQRTPI2 + ( x - 0.50 ) * log ( y ) - x - log ( sinpiy ) - r8_lgmc ( y );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_lnrel ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_LNREL evaluates log ( 1 + X ) for an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_LNREL, the value of LOG ( 1 + X ).
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp alnrcs[43] =
    {
        +0.10378693562743769800686267719098E+01,
        -0.13364301504908918098766041553133,
        +0.19408249135520563357926199374750E-01,
        -0.30107551127535777690376537776592E-02,
        +0.48694614797154850090456366509137E-03,
        -0.81054881893175356066809943008622E-04,
        +0.13778847799559524782938251496059E-04,
        -0.23802210894358970251369992914935E-05,
        +0.41640416213865183476391859901989E-06,
        -0.73595828378075994984266837031998E-07,
        +0.13117611876241674949152294345011E-07,
        -0.23546709317742425136696092330175E-08,
        +0.42522773276034997775638052962567E-09,
        -0.77190894134840796826108107493300E-10,
        +0.14075746481359069909215356472191E-10,
        -0.25769072058024680627537078627584E-11,
        +0.47342406666294421849154395005938E-12,
        -0.87249012674742641745301263292675E-13,
        +0.16124614902740551465739833119115E-13,
        -0.29875652015665773006710792416815E-14,
        +0.55480701209082887983041321697279E-15,
        -0.10324619158271569595141333961932E-15,
        +0.19250239203049851177878503244868E-16,
        -0.35955073465265150011189707844266E-17,
        +0.67264542537876857892194574226773E-18,
        -0.12602624168735219252082425637546E-18,
        +0.23644884408606210044916158955519E-19,
        -0.44419377050807936898878389179733E-20,
        +0.83546594464034259016241293994666E-21,
        -0.15731559416479562574899253521066E-21,
        +0.29653128740247422686154369706666E-22,
        -0.55949583481815947292156013226666E-23,
        +0.10566354268835681048187284138666E-23,
        -0.19972483680670204548314999466666E-24,
        +0.37782977818839361421049855999999E-25,
        -0.71531586889081740345038165333333E-26,
        +0.13552488463674213646502024533333E-26,
        -0.25694673048487567430079829333333E-27,
        +0.48747756066216949076459519999999E-28,
        -0.92542112530849715321132373333333E-29,
        +0.17578597841760239233269760000000E-29,
        -0.33410026677731010351377066666666E-30,
        +0.63533936180236187354180266666666E-31
    };
    static dim_typ nlnrel = 0;
    ityp value;
    static ityp xmin = 0.00;

    if ( nlnrel == 0 )
    {
        nlnrel = r8_inits ( alnrcs, 43, 0.10 * r8_mach ( 3 ) );
        xmin = - 1.00 + sqrt ( r8_mach ( 4 ) );
    }

    if ( x <= - 1.0 )
    {
    	result = MAX_VAL; 
        return &result;
    }

	result = abs ( x ) <= 0.375 ? x * ( 1.00 - x * r8_csevl ( x / 0.375, alnrcs, nlnrel ) ) : log ( 1.0 + x );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_mod ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_MOD returns the remainder of R8 division.
  Discussion:
    If
      REM = R8_MOD ( X, Y )
      RMULT = ( X - REM ) / Y
    then
      X = Y * RMULT + REM
    where REM has the same sign as X, and abs ( REM ) < Y.
  Example:
        X         Y     R8_MOD   R8_MOD  Factorization
      107        50       7     107 =  2 *  50 + 7
      107       -50       7     107 = -2 * -50 + 7
     -107        50      -7    -107 = -2 *  50 - 7
     -107       -50      -7    -107 =  2 * -50 - 7
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 June 2007
  Author:
    John Burkardt
  Parameters:
    Input, double X, the number to be divided.
    Input, double Y, the number that divides X.
    Output, double R8_MOD, the remainder when X is divided by Y.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp y = a_data[1];
	
    ityp value;
    if ( y == 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }
    
    value = x - ( ( ityp ) ( ( dim_typ ) ( x / y ) ) ) * y;
    
    result = value + (x < 0.00 && 0.00 < value) ? -abs(y):abs(y);
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_pak ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_PAK packs a base 2 exponent into an R8.
  Discussion:
    This routine is almost the inverse of R8_UPAK.  It is not exactly
    the inverse, because abs(x) need not be between 0.5 and 1.0.
    If both R8_PAK and 2.0^n were known to be in range, we could compute
    R8_PAK = x * 2.0^n .
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    C version by John Burkardt.
  Parameters:
    Input, double Y, the mantissa.
    Input, int N, the exponent.
    Output, double R8_PAK, the packed value.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp y = s_data->a1;
	
    static ityp aln210 = 3.321928094887362347870319429489;
    ityp aln2b;
    static dim_typ nmax = 0;
    static dim_typ nmin = 0;
    dim_typ nsum, ny;
    ityp value;

    if ( nmin == 0 )
    {
        aln2b = 1.00;
        if ( i4_mach ( 10 ) != 2 );
            aln2b = r8_mach ( 5 ) * aln210;
        nmin = aln2b * ( ityp ) ( i4_mach ( 15 ) );
        nmax = aln2b * ( ityp ) ( i4_mach ( 16 ) );
    }

    r8_upak ( y, &value, &ny );

    nsum = n + ny;

    if ( nsum < nmin || nmax < nsum )
    {
    	result = 0.00;
        return &result;
    }

    while ( nsum < 0 )
    {
        value = 0.50 * value;
        ++ nsum;
    }

    while ( 0 < nsum )
    {
        value = 2.00 * value;
        -- nsum;
    }

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_poch ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_POCH evaluates Pochhammer's function of R8 arguments.
  Discussion:
    POCH ( A, X ) = Gamma ( A + X ) / Gamma ( A ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, X, the arguments.
    Output, double R8_POCH, the Pochhammer function of A and X.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp x = a_data[1];
	
    ityp absa;
    ityp absax;
    ityp alnga;
    ityp alngax;
    ityp ax;
    ityp b;
    ityp cospia;
    ityp cospix;
    ityp den;
    static ityp eps = 0.00;
    ityp err;
    ityp errpch;
    dim_typ i, ia, n;
    ityp sgnga;
    ityp sgngax;
    ityp sinpia;
    ityp sinpix;
    static ityp sqeps = 0.00;
    ityp value;

    if ( eps == 0.00 )
    {
        eps = r8_mach ( 4 );
        sqeps = sqrt ( eps );
    }

    ax = a + x;

    if ( ax <= 0.00 && r8_aint ( ax ) == ax )
    {
        if ( 0.00 < a || r8_aint ( a ) != a )
        {
        	result = MAX_VAL;
            return &result;
        }

        /*
        We know here that both A+X and A are non-positive integers.
        */
        if ( x == 0.00 )
            value = 1.00;
        else if ( - 20.00 < MIN ( a + x, a ) )
        {
            n = (dim_typ ) ( x );
            ia = ( dim_typ ) ( a );
            value = r8_mop ( n ) * r8_fac ( - ia ) / r8_fac ( - ia - n );
        }
        else
        {
            n = ( dim_typ ) ( x );
            value = r8_mop ( n ) * exp ( ( a - 0.50 )* r8_lnrel ( x / ( a - 1.00 ) )+ x * log ( - a + 1.00 - x ) - x+ r8_lgmc ( - a + 1.0 )- r8_lgmc ( - a - x + 1.0 ) );
        }
        
        result = value;
        return &result;
    }
    /*
    A + X is not zero or a negative integer.
    */
    if ( a <= 0.00 && r8_aint ( a ) == a )
    {
    	result = 0.00;
        return &result;
    }

    n = abs ( x );
    /*
    X is a small non-positive integer, presummably a common case.
    */
    if ( ( ityp ) ( n ) == x && n <= 20 )
    {
        value = 1.00;
        for ( i = 1; i <= n; ++i)
            value *= ( a + ( ityp ) ( i - 1 ) );
            
        result = value;
        return &result;
    }

    absax = abs ( a + x );
    absa = abs ( a );

    if ( MAX ( absax, absa ) <= 20.00 )
    {
    	result = r8_gamma ( a + x ) * r8_gamr ( a );
        return &result;
    }

    if ( 0.50 * absa < abs ( x ) )
    {
        r8_lgams ( a + x, &alngax, &sgngax );
        r8_lgams ( a, &alnga, &sgnga );
        result = sgngax * sgnga * exp ( alngax - alnga ); 
        return &result;
    }
    /*
    abs(x) is small and both abs(a+x) and abs(a) are large.  thus,
    a+x and a must have the same sign.  for negative a, we use
    gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
    sin(M_PI*a)/sin(M_PI*(a+x))
    */
    b = a<0.00 ? -a-x+1.00 : a;
    value = exp ( ( b - 0.50 ) * r8_lnrel ( x / b )+ x * log ( b + x ) - x + r8_lgmc ( b + x ) - r8_lgmc ( b ) );

    if ( 0.00 <= a || value == 0.00 )
    {
    	result = value;
        return &result;
    }

    cospix = cos ( M_PI * x );
    sinpix = sin ( M_PI * x );
    cospia = cos ( M_PI * a );
    sinpia = sin ( M_PI * a );

    err = ( abs ( x ) * ( abs ( sinpix )+ abs ( cospia * cospix / sinpia ) )+ abs ( a * sinpix ) / sinpia / sinpia ) * M_PI;
    err = abs ( x ) * ( 1.00 + log ( b ) ) + err / abs ( den );
    
    result = value / cospix + cospia * sinpix / sinpia;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _comp_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMP_NEXT computes the compositions of the integer N into K parts.
  Discussion:
    A composition of the integer N into K parts is an ordered sequence
    of K nonnegative integers which sum to N.  The compositions (1,2,1)
    and (1,1,2) are considered to be distinct.
    The routine computes one composition on each call until there are no more.
    For instance, one composition of 6 into 3 parts is
    3+2+1, another would be 6+0+0.
    On the first call to this routine, set MORE = FALSE.  The routine
    will compute the first element in the sequence of compositions, and
    return it, as well as setting MORE = TRUE.  If more compositions
    are desired, call again, and again.  Each time, the routine will
    return with a new composition.
    However, when the LAST composition in the sequence is computed
    and returned, the routine will reset MORE to FALSE, signaling that
    the end of the sequence has been reached.
    This routine originally used a STATICE statement to maintain the
    variables H and T.  I have decided (based on an wasting an
    entire morning trying to track down a problem) that it is safer
    to pass these variables as arguments, even though the user should
    never alter them.  This allows this routine to safely shuffle
    between several ongoing calculations.
    There are 28 compositions of 6 into three parts.  This routine will
    produce those compositions in the following order:
     I         A
     -     ---------
     1     6   0   0
     2     5   1   0
     3     4   2   0
     4     3   3   0
     5     2   4   0
     6     1   5   0
     7     0   6   0
     8     5   0   1
     9     4   1   1
    10     3   2   1
    11     2   3   1
    12     1   4   1
    13     0   5   1
    14     4   0   2
    15     3   1   2
    16     2   2   2
    17     1   3   2
    18     0   4   2
    19     3   0   3
    20     2   1   3
    21     1   2   3
    22     0   3   3
    23     2   0   4
    24     1   1   4
    25     0   2   4
    26     1   0   5
    27     0   1   5
    28     0   0   6
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2008
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the integer whose compositions are desired.
    Input, int K, the number of parts in the composition.
    Input/output, int A[K], the parts of the composition.
    Input/output, int *MORE.
    Set MORE = FALSE on first call.  It will be reset to TRUE on return
    with a new composition.  Each new call returns another composition until
    MORE is set to FALSE when the last composition has been computed
    and returned.
    Input/output, int *H, *T, two internal parameters needed for the
    computation.  The user should allocate space for these in the calling
    program, include them in the calling sequence, but never alter them!
*/
{
	const _2dt4pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	int * a = s_data->a2;
	int * more = s_data->a3;
	int * h = s_data->a4;
	int * t = s_data->a5;

	if ( !( *more ) )
	{
		*t = n;
		*h = 0;
		a[0] = n;
		for (dim_typ i = 1; i < k; ++i )
			a[i] = 0;
	}
	else
	{
		if ( 1 < *t )
			*h = 0;
		*h = *h + 1;
		*t = a[*h-1];
		a[*h-1] = 0;
		a[0] = *t - 1;
		a[*h] = a[*h] + 1;
	}

	*more = ( a[k-1] != n );

	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_poch1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_POCH1 evaluates a quantity related to Pochhammer's symbol.
  Discussion:
    Evaluate a generalization of Pochhammer's symbol for special
    situations that require especially accurate values when x is small in
      poch1(a,x) = (poch(a,x)-1)/x
                 = (gamma(a+x)/gamma(a) - 1.0)/x .
    This specification is particularly suited for stably computing
    expressions such as
   (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
           = poch1(a,x) - poch1(b,x)
    Note that poch1(a,0.0) = psi(a)
    When abs(x) is so small that substantial cancellation will occur if
    the straightforward formula is used, we  use an expansion due
    to fields and discussed by y. l. luke, the special functions and their
    approximations, vol. 1, academic press, 1969, page 34.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double A, the parameter.
    Input, double X, the evaluation point.
    Output, double R8_POCH1, the value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp x = a_data[1];
	
    ityp absa;
    ityp absx;
    static ityp alneps = 0.00;
    ityp alnvar;
    ityp b;
    static ityp bern[20] =
    {
        +0.833333333333333333333333333333333E-01,
        -0.138888888888888888888888888888888E-02,
        +0.330687830687830687830687830687830E-04,
        -0.826719576719576719576719576719576E-06,
        +0.208767569878680989792100903212014E-07,
        -0.528419013868749318484768220217955E-09,
        +0.133825365306846788328269809751291E-10,
        -0.338968029632258286683019539124944E-12,
        +0.858606205627784456413590545042562E-14,
        -0.217486869855806187304151642386591E-15,
        +0.550900282836022951520265260890225E-17,
        -0.139544646858125233407076862640635E-18,
        +0.353470703962946747169322997780379E-20,
        -0.895351742703754685040261131811274E-22,
        +0.226795245233768306031095073886816E-23,
        -0.574472439520264523834847971943400E-24,
        +0.145517247561486490186626486727132E-26,
        -0.368599494066531017818178247990866E-28,
        +0.933673425709504467203255515278562E-30,
        -0.236502241570062993455963519636983E-31
    };
    ityp binv;
    ityp bp;
    ityp gbern[21];
    ityp gbk;
    dim_typ i;
    dim_typ ii;
    dim_typ incr;
    dim_typ j;
    dim_typ k;
    dim_typ ndx;
    dim_typ nterms;
    ityp poly1;
    ityp q;
    ityp rho;
    ityp sinpxx;
    ityp sinpx2;
    static ityp sqtbig = 0.00;
    ityp term;
    ityp trig;
    ityp value;
    ityp var;
    ityp var2;

    if ( sqtbig == 0.00 )
    {
        sqtbig = 1.00 / sqrt ( 24.00 * r8_mach ( 1 ) );
        alneps = log ( r8_mach ( 3 ) );
    }

    if ( x == 0.00 )
    {
    	result = r8_psi ( a );
        return &result;
    }

    absx = abs ( x );
    absa = abs ( a );

    if ( 0.10 * absa < absx || 0.10 < absx * log ( MAX ( absa, 2.00 ) ) )
    {
    	result = ( r8_poch ( a, x ) - 1.00 ) / x;
        return &result;
    }

    bp = a<-0.5 ? 1.00-a-x : a;
    incr = 0 + (bp<10.00)*(r8_aint ( 11.0 - bp ));
    b = bp + ( ityp ) ( incr );
    var = b + 0.50 * ( x - 1.00 );
    alnvar = log ( var );
    q = x * alnvar;
    poly1 = 0.00;

    if ( var < sqtbig )
    {
        var2 = 1.00 / var / var;

        rho = 0.50 * ( x + 1.00 );
        gbern[0] = 1.00;
        gbern[1] = - rho / 12.00;
        term = var2;
        poly1 = gbern[1] * term;

        nterms = ( dim_typ ) ( - 0.50 * alneps / alnvar + 1.0 );

        if ( 20 < nterms )
        {
        	result = MAX_VAL;
            return &result;
        }

        for ( k = 2; k <= nterms; ++k )
        {
            gbk = 0.00;
            for ( j = 1; j <= k; ++j )
            {
                ndx = k - j + 1;
                gbk += bern[ndx-1] * gbern[j-1];
            }
            gbern[k] = - rho * gbk / ( ityp ) ( k );
            term *= ( ( ityp ) ( (k<<1) - 2 ) - x )* ( ( ityp ) ( (k<<1) - 1 ) - x ) * var2;
            poly1 += gbern[k] * term;
        }
    }
    poly1 = ( x - 1.00 ) * poly1;
    value = r8_exprel ( q ) * ( alnvar + q * poly1 ) + poly1;
    /*
    We have r8_poch1(b,x), but bp is small, so we use backwards recursion
    to obtain r8_poch1(bp,x).
    */
    for ( ii = 1; ii <= incr; ++ii )
    {
        i = incr - ii;
        binv = 1.00 / ( bp + ( ityp ) ( i ) );
        value = ( value - binv ) / ( 1.00 + x * binv );
    }

    if ( bp == a )
    {
    	result = value;
        return &result;
    }
    /*
    We have r8_poch1(bp,x), but a is lt -0.5.  we therefore use a reflection
    formula to obtain r8_poch1(a,x).
    */
    sinpxx = sin ( M_PI * x ) / x;
    sinpx2 = sin ( 0.50 * M_PI * x );
    trig = sinpxx * cotan ( M_PI * b ) - 2.00 * sinpx2 * ( sinpx2 / x );
    
    result = trig + ( 1.00 + x * trig ) * value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_ren ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_REN is a simple random number generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    C version by John Burkardt.
  Reference:
    Malcolm Pike, David Hill,
    Algorithm 266:
    Pseudo-Random Numbers,
    Communications of the ACM,
    Volume 8, Number 10, October 1965, page 605.
  Parameters:
    Output, double R8_REN, the random value.
*/
{
	static ityp result = MAX_VAL;
	
    static int iy = 100001;
    ityp value;
    iy = iy * 125;
    
    result = ( ityp ) ( iy - ( iy / 2796203 ) * 2796203 ) / 2796203.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_shi ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SHI evaluates the hyperbolic sine integral Shi of an R8 argument.
  Discussion:
    Shi ( x ) = Integral ( 0 <= t <= x ) sinh ( t ) dt / t
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_SHI, the hyperbolic sine integral
    Shi evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp absx;
    static dim_typ nshi = 0;
    static ityp shics[10] =
    {
        0.0078372685688900950695200984317332E+00,
        0.0039227664934234563972697574427225E+00,
        0.0000041346787887617266746747908275E+00,
        0.0000000024707480372882742135145302E+00,
        0.0000000000009379295590763630457157E+00,
        0.0000000000000002451817019520867353E+00,
        0.0000000000000000000467416155257592E+00,
        0.0000000000000000000000067803072389E+00,
        0.0000000000000000000000000007731289E+00,
        0.0000000000000000000000000000000711E+00
    };
    ityp value;
    static ityp xsml = 0.00;

    if ( nshi == 0 )
    {
        nshi = r8_inits ( shics, 10, 0.1 * r8_mach ( 3 ) );
        xsml = sqrt ( r8_mach ( 3 ) );
    }

    absx = abs ( x );

    if ( absx <= xsml )
        value = x;
    else if ( absx <= 0.375 )
        value = x * ( 1.00 + r8_csevl ( 128.00 * x * x / 9.00 - 1.00, shics, nshi ) );
    else
        value = 0.50 * ( r8_ei ( x ) + r8_e1 ( x ) );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_si ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SI evaluates the sine integral Si of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_SI, the sine integral Si evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp absx;
    ityp cosx;
    ityp f;
    ityp g;
    static int nsi = 0;
    static ityp sics[18] =
    {
        -0.1315646598184841928904275173000457,
        -0.2776578526973601892048287660157299,
        0.0354414054866659179749135464710086,
        -0.0025631631447933977658752788361530,
        0.0001162365390497009281264921482985,
        -0.0000035904327241606042670004347148,
        0.0000000802342123705710162308652976,
        -0.0000000013562997692540250649931846,
        0.0000000000179440721599736775567759,
        -0.0000000000001908387343087145490737,
        0.0000000000000016669989586824330853,
        -0.0000000000000000121730988368503042,
        0.0000000000000000000754181866993865,
        -0.0000000000000000000004014178842446,
        0.0000000000000000000000018553690716,
        -0.0000000000000000000000000075166966,
        0.0000000000000000000000000000269113,
        -0.0000000000000000000000000000000858
    };
    ityp value;
    static ityp xsml = 0.00;

    if ( nsi == 0 )
    {
        nsi = r8_inits ( sics, 18, 0.1 * r8_mach ( 3 ) );
        xsml = sqrt ( r8_mach ( 3 ) );
    }

    absx = abs ( x );

    if ( absx < xsml )
        value = x;
    else if ( absx <= 4.00 )
        value = x * ( 0.75 + r8_csevl ( ( x * x - 8.00 ) * 0.125, sics, nsi ) );
    else
    {
        r8_sifg ( absx, &f, &g );
        cosx = cos ( absx );
        value = M_PI_2 - f * cosx - g * sin ( x );
        if ( x < 0.0 )
            value *= -1;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_sifg ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SIFG is a utility routine.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double *F, *G.
*/
{
	const it2pit * const s_data = data;
	const register ityp x  = s_data->a0;
	ityp * f = s_data->a1;
	ityp * g = s_data->a2;
	
    static ityp f1cs[43] =
    {
        -0.1191081969051363610348201965828918,
        -0.0247823144996236247590074150823133,
        0.0011910281453357821268120363054457,
        -0.0000927027714388561748308600360706,
        0.0000093373141568270996868204582766,
        -0.0000011058287820557143938979426306,
        0.0000001464772071460162169336550799,
        -0.0000000210694496287689532601227548,
        0.0000000032293492366848236382857374,
        -0.0000000005206529617529375828014986,
        0.0000000000874878884570278750268316,
        -0.0000000000152176187056123668294574,
        0.0000000000027257192405419573900583,
        -0.0000000000005007053075968556290255,
        0.0000000000000940240902726068511779,
        -0.0000000000000180014444791803678336,
        0.0000000000000035062621432741785826,
        -0.0000000000000006935282926769149709,
        0.0000000000000001390925136454216568,
        -0.0000000000000000282486885074170585,
        0.0000000000000000058031305693579081,
        -0.0000000000000000012046901573375820,
        0.0000000000000000002525052443655940,
        -0.0000000000000000000533980268805594,
        0.0000000000000000000113855786274122,
        -0.0000000000000000000024462861505259,
        0.0000000000000000000005293659320439,
        -0.0000000000000000000001153184940277,
        0.0000000000000000000000252786568318,
        -0.0000000000000000000000055738645378,
        0.0000000000000000000000012358245621,
        -0.0000000000000000000000002754350842,
        0.0000000000000000000000000616906808,
        -0.0000000000000000000000000138817443,
        0.0000000000000000000000000031375329,
        -0.0000000000000000000000000007121249,
        0.0000000000000000000000000001622778,
        -0.0000000000000000000000000000371206,
        0.0000000000000000000000000000085221,
        -0.0000000000000000000000000000019633,
        0.0000000000000000000000000000004538,
        -0.0000000000000000000000000000001052,
        0.0000000000000000000000000000000245
    };
    static ityp f2cs[99] =
    {
        -0.03484092538970132330836049733745577,
        -0.01668422056779596873246786312278676,
        0.00067529012412377385045207859239727,
        -0.00005350666225447013628785577557429,
        0.00000626934217790075267050759431626,
        -0.00000095266388019916680677790414293,
        0.00000017456292242509880425504427666,
        -0.00000003687954030653093307097646628,
        0.00000000872026777051395264075816938,
        -0.00000000226019703919738748530423167,
        0.00000000063246249765250612520444877,
        -0.00000000018889118884717869240911480,
        0.00000000005967746729997813372620472,
        -0.00000000001980443117372239011196007,
        0.00000000000686413954772103383713264,
        -0.00000000000247310193070199106074890,
        0.00000000000092263594549941404196042,
        -0.00000000000035523634999261784497297,
        0.00000000000014076049625351591461820,
        -0.00000000000005726228499747652794311,
        0.00000000000002386537545413171810106,
        -0.00000000000001017141890764597142232,
        0.00000000000000442594531078364424968,
        -0.00000000000000196344933049189761979,
        0.00000000000000088688748314810461024,
        -0.00000000000000040743345027311546948,
        0.00000000000000019016837215675339859,
        -0.00000000000000009009707297478042442,
        0.00000000000000004329211274095668667,
        -0.00000000000000002108144465322479526,
        0.00000000000000001039637907026452274,
        -0.00000000000000000518891007948931936,
        0.00000000000000000261955324869899371,
        -0.00000000000000000133690399951301570,
        0.00000000000000000068941057702931664,
        -0.00000000000000000035905362610437250,
        0.00000000000000000018878077255791706,
        -0.00000000000000000010016125265594380,
        0.00000000000000000005360725691578228,
        -0.00000000000000000002893198974944827,
        0.00000000000000000001574065100202625,
        -0.00000000000000000000863027106431206,
        0.00000000000000000000476715602862288,
        -0.00000000000000000000265222739998504,
        0.00000000000000000000148582865063866,
        -0.00000000000000000000083797235923135,
        0.00000000000000000000047565916422711,
        -0.00000000000000000000027169073353112,
        0.00000000000000000000015612738881686,
        -0.00000000000000000000009024555078347,
        0.00000000000000000000005246097049119,
        -0.00000000000000000000003066450818697,
        0.00000000000000000000001801996250957,
        -0.00000000000000000000001064443050752,
        0.00000000000000000000000631942158881,
        -0.00000000000000000000000377013812246,
        0.00000000000000000000000225997542918,
        -0.00000000000000000000000136100844814,
        0.00000000000000000000000082333232003,
        -0.00000000000000000000000050025986091,
        0.00000000000000000000000030526245684,
        -0.00000000000000000000000018705164021,
        0.00000000000000000000000011508404393,
        -0.00000000000000000000000007108714611,
        0.00000000000000000000000004408065533,
        -0.00000000000000000000000002743760867,
        0.00000000000000000000000001714144851,
        -0.00000000000000000000000001074768860,
        0.00000000000000000000000000676259777,
        -0.00000000000000000000000000426981348,
        0.00000000000000000000000000270500637,
        -0.00000000000000000000000000171933331,
        0.00000000000000000000000000109636138,
        -0.00000000000000000000000000070132573,
        0.00000000000000000000000000045001784,
        -0.00000000000000000000000000028963835,
        0.00000000000000000000000000018697009,
        -0.00000000000000000000000000012104646,
        0.00000000000000000000000000007859065,
        -0.00000000000000000000000000005116867,
        0.00000000000000000000000000003340627,
        -0.00000000000000000000000000002186851,
        0.00000000000000000000000000001435340,
        -0.00000000000000000000000000000944523,
        0.00000000000000000000000000000623117,
        -0.00000000000000000000000000000412101,
        0.00000000000000000000000000000273208,
        -0.00000000000000000000000000000181558,
        0.00000000000000000000000000000120934,
        -0.00000000000000000000000000000080737,
        0.00000000000000000000000000000054022,
        -0.00000000000000000000000000000036227,
        0.00000000000000000000000000000024348,
        -0.00000000000000000000000000000016401,
        0.00000000000000000000000000000011074,
        -0.00000000000000000000000000000007497,
        0.00000000000000000000000000000005091,
        -0.00000000000000000000000000000003470,
        0.00000000000000000000000000000002377
    };
    static ityp g1cs[44] =
    {
        -0.3040578798253495954499726682091083,
        -0.0566890984597120587731339156118269,
        0.0039046158173275643919984071554082,
        -0.0003746075959202260618619339867489,
        0.0000435431556559843679552220840065,
        -0.0000057417294453025046561970723475,
        0.0000008282552104502629741937616492,
        -0.0000001278245892594642727883913223,
        0.0000000207978352948687884439257529,
        -0.0000000035313205921990798042032682,
        0.0000000006210824236308951068631449,
        -0.0000000001125215474446292649336987,
        0.0000000000209088917684421605267019,
        -0.0000000000039715831737681727689158,
        0.0000000000007690431314272089939005,
        -0.0000000000001514696742731613519826,
        0.0000000000000302892146552359684119,
        -0.0000000000000061399703834708825400,
        0.0000000000000012600605829510933553,
        -0.0000000000000002615029250939483683,
        0.0000000000000000548278844891796821,
        -0.0000000000000000116038182129526571,
        0.0000000000000000024771654107129795,
        -0.0000000000000000005330672753223389,
        0.0000000000000000001155666075598465,
        -0.0000000000000000000252280547744957,
        0.0000000000000000000055429038550786,
        -0.0000000000000000000012252208421297,
        0.0000000000000000000002723664318684,
        -0.0000000000000000000000608707831422,
        0.0000000000000000000000136724874476,
        -0.0000000000000000000000030856626806,
        0.0000000000000000000000006995212319,
        -0.0000000000000000000000001592587569,
        0.0000000000000000000000000364051056,
        -0.0000000000000000000000000083539465,
        0.0000000000000000000000000019240303,
        -0.0000000000000000000000000004446816,
        0.0000000000000000000000000001031182,
        -0.0000000000000000000000000000239887,
        0.0000000000000000000000000000055976,
        -0.0000000000000000000000000000013100,
        0.0000000000000000000000000000003074,
        -0.0000000000000000000000000000000723
    };
    static ityp g2cs[44] =
    {
        -0.1211802894731646263541834046858267,
        -0.0316761386394950286701407923505610,
        0.0013383199778862680163819429492182,
        -0.0000895511011392252425531905069518,
        0.0000079155562961718213115249467924,
        -0.0000008438793322241520181418982080,
        0.0000001029980425677530146647227274,
        -0.0000000139295750605183835795834444,
        0.0000000020422703959875980400677594,
        -0.0000000003196534694206427035434752,
        0.0000000000528147832657267698615312,
        -0.0000000000091339554672671033735289,
        0.0000000000016426251238967760444819,
        -0.0000000000003055897039322660002410,
        0.0000000000000585655825785779717892,
        -0.0000000000000115229197730940120563,
        0.0000000000000023209469119988537310,
        -0.0000000000000004774355834177535025,
        0.0000000000000001000996765800180573,
        -0.0000000000000000213533778082256704,
        0.0000000000000000046277190777367671,
        -0.0000000000000000010175807410227657,
        0.0000000000000000002267657399884672,
        -0.0000000000000000000511630776076426,
        0.0000000000000000000116767014913108,
        -0.0000000000000000000026935427672470,
        0.0000000000000000000006275665841146,
        -0.0000000000000000000001475880557531,
        0.0000000000000000000000350145314739,
        -0.0000000000000000000000083757732152,
        0.0000000000000000000000020191815152,
        -0.0000000000000000000000004903567705,
        0.0000000000000000000000001199123348,
        -0.0000000000000000000000000295170610,
        0.0000000000000000000000000073113112,
        -0.0000000000000000000000000018217843,
        0.0000000000000000000000000004565148,
        -0.0000000000000000000000000001150151,
        0.0000000000000000000000000000291267,
        -0.0000000000000000000000000000074125,
        0.0000000000000000000000000000018953,
        -0.0000000000000000000000000000004868,
        0.0000000000000000000000000000001256,
        -0.0000000000000000000000000000000325
    };
    static ityp g3cs[56] =
    {
        -0.0280574367809472928402815264335299,
        -0.0137271597162236975409100508089556,
        0.0002894032638760296027448941273751,
        -0.0000114129239391197145908743622517,
        0.0000006813965590726242997720207302,
        -0.0000000547952289604652363669058052,
        0.0000000055207429918212529109406521,
        -0.0000000006641464199322920022491428,
        0.0000000000922373663487041108564960,
        -0.0000000000144299088886682862611718,
        0.0000000000024963904892030710248705,
        -0.0000000000004708240675875244722971,
        0.0000000000000957217659216759988140,
        -0.0000000000000207889966095809030537,
        0.0000000000000047875099970877431627,
        -0.0000000000000011619070583377173759,
        0.0000000000000002956508969267836974,
        -0.0000000000000000785294988256492025,
        0.0000000000000000216922264368256612,
        -0.0000000000000000062113515831676342,
        0.0000000000000000018384568838450977,
        -0.0000000000000000005610887482137276,
        0.0000000000000000001761862805280062,
        -0.0000000000000000000568111050541451,
        0.0000000000000000000187786279582313,
        -0.0000000000000000000063531694151124,
        0.0000000000000000000021968802368238,
        -0.0000000000000000000007754666550395,
        0.0000000000000000000002791018356581,
        -0.0000000000000000000001023178525247,
        0.0000000000000000000000381693403919,
        -0.0000000000000000000000144767895606,
        0.0000000000000000000000055779512634,
        -0.0000000000000000000000021817239071,
        0.0000000000000000000000008656646309,
        -0.0000000000000000000000003482157895,
        0.0000000000000000000000001419188130,
        -0.0000000000000000000000000585714314,
        0.0000000000000000000000000244660482,
        -0.0000000000000000000000000103387099,
        0.0000000000000000000000000044177299,
        -0.0000000000000000000000000019080079,
        0.0000000000000000000000000008326038,
        -0.0000000000000000000000000003669553,
        0.0000000000000000000000000001632875,
        -0.0000000000000000000000000000733357,
        0.0000000000000000000000000000332327,
        -0.0000000000000000000000000000151906,
        0.0000000000000000000000000000070020,
        -0.0000000000000000000000000000032539,
        0.0000000000000000000000000000015240,
        -0.0000000000000000000000000000007193,
        0.0000000000000000000000000000003420,
        -0.0000000000000000000000000000001638,
        0.0000000000000000000000000000000790,
        -0.0000000000000000000000000000000383
    };
    static dim_typ nf1 = 0;
    static dim_typ nf2 = 0;
    static dim_typ ng1 = 0;
    static dim_typ ng2 = 0;
    static dim_typ ng3 = 0;
    ityp tol;
    static ityp xbig = 0.00;
    static ityp xbnd = 0.00;
    static ityp xbndg = 0.00;
    static ityp xmaxf = 0.00;
    static ityp xmaxg = 0.00;

    if ( nf1 == 0 )
    {
        tol = 0.10 * r8_mach ( 3 );
        nf1 = r8_inits ( f1cs, 43, tol );
        nf2 = r8_inits ( f2cs, 99, tol );
        ng1 = r8_inits ( g1cs, 44, tol );
        ng2 = r8_inits ( g2cs, 44, tol );
        ng3 = r8_inits ( g3cs, 56, tol );
        xbig = sqrt ( 1.00 / r8_mach ( 3 ) );
        xmaxf = exp ( MIN ( - log ( r8_mach ( 1 ) ),
        log ( r8_mach ( 2 ) ) ) - 0.01 );
        xmaxg = 1.00 / sqrt ( r8_mach ( 1 ) );
        xbnd = sqrt ( 50.00 );
        xbndg = sqrt ( 200.00 );
    }

    if ( x < 4.00 )
        return NULL;
    else if ( x <= xbnd )
    {
        *f = ( 1.00 + r8_csevl ( ( 1.00 / x / x - 0.04125 )/ 0.02125, f1cs, nf1 ) ) / x;
        *g = ( 1.00 + r8_csevl ( ( 1.00 / x / x - 0.04125 )/ 0.02125, g1cs, ng1 ) ) / x / x;
    }
    else if ( x <= xbig )
    {
        *f = ( 1.00 + r8_csevl ( 100.00 / x / x - 1.00, f2cs, nf2 ) ) / x;
        *g = x<=xbndg ? ( 1.00 + r8_csevl ( ( 10000.00 / x / x - 125.00 )/ 75.00, g2cs, ng2 ) ) / x / x : ( 1.00 + r8_csevl ( 400.00 / x / x - 1.00, g3cs, ng3 ) ) / x / x;
    }
    else
    {
        *f = 0.00 + (x<xmaxf)*(1.00/x);
        *g = 0.00 + (x<xmaxg)*(1.00/x/x);
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_sin_deg ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SIN_DEG evaluates the sine of an R8 argument in degrees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument in degrees.
    Output, double R8_SIN_DEG, the sine of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    dim_typ n;
    static ityp raddeg = 0.017453292519943295769236907684886E+00;
    ityp value;

    value = sin ( raddeg * x );

    if ( r8_mod ( x, 90.0E+00 ) == 0.0E+00 )
    {
        n = ( int ) ( abs ( x ) / 90.0E+00 + 0.5E+00 );
        n = ( n % 2 );
        if ( n == 0 )
            value = 0.0E+00;
        else if ( value < 0.0E+00 )
            value = - 1.0E+00;
        else
            value = + 1.0E+00;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_spence ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SPENCE evaluates a form of Spence's function for an R8 argument.
  Discussion:
    This function evaluates a form of Spence's function defined by
      f(x) = Integral ( 0 <= y <= x ) - log ( abs ( 1 - y ) ) / y dy
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions, page 1004,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
    K Mitchell,
    Tables of the function Integral ( 0 < y < x ) - log | 1 - y | dy / y
    with an account of some properties of this and related functions,
    Philosophical Magazine,
    Volume 40, pages 351-368, 1949.
  Parameters:
    Input, double X, the argument.
    Output, double R8_SPENCE, Spence's function evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp aln;
    static dim_typ nspenc = 0;
    static ityp pi26 = +1.644934066848226436472415166646025189219;
    ityp spencs[38] =
    {
        +0.1527365598892405872946684910028,
        +0.8169658058051014403501838185271E-01,
        +0.5814157140778730872977350641182E-02,
        +0.5371619814541527542247889005319E-03,
        +0.5724704675185826233210603054782E-04,
        +0.6674546121649336343607835438589E-05,
        +0.8276467339715676981584391689011E-06,
        +0.1073315673030678951270005873354E-06,
        +0.1440077294303239402334590331513E-07,
        +0.1984442029965906367898877139608E-08,
        +0.2794005822163638720201994821615E-09,
        +0.4003991310883311823072580445908E-10,
        +0.5823462892044638471368135835757E-11,
        +0.8576708692638689278097914771224E-12,
        +0.1276862586280193045989483033433E-12,
        +0.1918826209042517081162380416062E-13,
        +0.2907319206977138177795799719673E-14,
        +0.4437112685276780462557473641745E-15,
        +0.6815727787414599527867359135607E-16,
        +0.1053017386015574429547019416644E-16,
        +0.1635389806752377100051821734570E-17,
        +0.2551852874940463932310901642581E-18,
        +0.3999020621999360112770470379519E-19,
        +0.6291501645216811876514149171199E-20,
        +0.9933827435675677643803887752533E-21,
        +0.1573679570749964816721763805866E-21,
        +0.2500595316849476129369270954666E-22,
        +0.3984740918383811139210663253333E-23,
        +0.6366473210082843892691326293333E-24,
        +0.1019674287239678367077061973333E-24,
        +0.1636881058913518841111074133333E-25,
        +0.2633310439417650117345279999999E-26,
        +0.4244811560123976817224362666666E-27,
        +0.6855411983680052916824746666666E-28,
        +0.1109122433438056434018986666666E-28,
        +0.1797431304999891457365333333333E-29,
        +0.2917505845976095173290666666666E-30,
        +0.4742646808928671061333333333333E-31
    };
    ityp value;
    static ityp xbig = 0.00;

    if ( nspenc == 0 )
    {
        nspenc = r8_inits ( spencs, 38, 0.10 * r8_mach ( 3 ) );
        xbig = 1.00 / r8_mach ( 3 );
    }

    if ( x <= - xbig )
    {
        aln = log ( 1.00 - x );
        value = - pi26  - 0.50 * aln * ( 2.00 * log ( - x ) - aln );
    }
    else if ( x <= - 1.00 )
    {
        aln = log ( 1.00 - x );
        value = - pi26 - 0.50 * aln * ( 2.00* log ( - x ) - aln ) + ( 1.00 + r8_csevl (4.00 / ( 1.00 - x ) - 1.00, spencs, nspenc ) ) / ( 1.0 - x );
    }
    else if ( x <= 0.0 )
        value = - 0.50 * log ( 1.00 - x )* log ( 1.00 - x ) - x * ( 1.00 + r8_csevl (4.00 * x / ( x - 1.00 ) - 1.00, spencs, nspenc ) ) / ( x - 1.00 );
    else if ( x <= 0.5 )
        value = x * ( 1.00 + r8_csevl ( 4.00 * x - 1.00, spencs, nspenc ) );
    else if ( x < 1.0 )
        value = pi26 - log ( x ) * log ( 1.00 - x )- ( 1.00 - x ) * ( 1.00 + r8_csevl ( 4.00* ( 1.00 - x ) - 1.00, spencs, nspenc ) );
    else if ( x == 1.0 )
        value = pi26;
    else if ( x <= 2.0 )
        value = pi26 - 0.50 * log ( x )* log ( ( x - 1.00 ) * ( x - 1.00 ) / x )+ ( x - 1.00 ) * ( 1.00 + r8_csevl ( 4.00* ( x - 1.00 ) / x - 1.00, spencs, nspenc ) ) / x;
    else if ( x < xbig )
        value = 2.00 * pi26 - 0.50 * log ( x ) * log ( x )- ( 1.00 + r8_csevl ( 4.00 / x - 1.00, spencs, nspenc ) ) / x;
    else
        value = 2.00 * pi26 - 0.50 * log ( x ) * log ( x );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_sqrt ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SQRT computes the square root of an R8.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the number whose square root is desired.
    Output, double R8_SQRT, the square root of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    dim_typ irem;
    dim_typ iter;
    dim_typ ixpnt;
    dim_typ n;
    static dim_typ niter = 0;
    static ityp sqrt2[3] =
    {
        0.70710678118654752440084436210485,
        1.00,
        1.41421356237309504880168872420970
    };
    ityp value;
    ityp y;

    if ( niter == 0 )
        niter = 1.443 * log ( - 0.104 * log ( 0.10 * r8_mach ( 3 ) ) ) + 1.00;

    if ( x < 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }
    else if ( x == 0.00 )
        value = 0.00;
    else
    {
        r8_upak ( x, &y, &n );
        ixpnt = n / 2;
        irem = n - (ixpnt<<1) + 2;
        value = 0.261599 + y * ( 1.114292 + y * ( -0.516888 + y * 0.141067 ) );

        for ( iter = 1; iter <= niter; ++iter )
            value += 0.50 * ( y - value * value ) / value;
        value = r8_pak ( sqrt2[irem-1] * value, ixpnt );
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_tan ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_TAN evaluates the tangent of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_TAN, the tangent of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp ainty;
    ityp ainty2;
    dim_typ ifn;
    static dim_typ nterms = 0;
    static ityp pi2rec = 0.011619772367581343075535053490057;
    ityp prodbg;
    static ityp sqeps = 0.00;
    static ityp tancs[19] =
    {
        +0.22627932763129357846578636531752,
        +0.43017913146548961775583410748067E-01,
        +0.68544610682565088756929473623461E-03,
        +0.11045326947597098383578849369696E-04,
        +0.17817477903926312943238512588940E-06,
        +0.28744968582365265947529646832471E-08,
        +0.46374854195902995494137478234363E-10,
        +0.74817609041556138502341633308215E-12,
        +0.12070497002957544801644516947824E-13,
        +0.19473610812823019305513858584533E-15,
        +0.31417224874732446504614586026666E-17,
        +0.50686132555800153941904891733333E-19,
        +0.81773105159836540043979946666666E-21,
        +0.13192643412147384408951466666666E-22,
        +0.21283995497042377309866666666666E-24,
        +0.34337960192345945292800000000000E-26,
        +0.55398222121173811200000000000000E-28,
        +0.89375227794352810666666666666666E-30,
        +0.14419111371369130666666666666666E-31
    };
    ityp value;
    static ityp xmax = 0.00;
    static ityp xsml = 0.00;
    ityp y;
    ityp yrem;

    if ( nterms == 0 )
    {
        nterms = r8_inits ( tancs, 19, 0.10 * r8_mach ( 3 ) );
        xmax = 1.00 / r8_mach ( 4 );
        xsml = sqrt ( 3.00 * r8_mach ( 3 ) );
        sqeps = sqrt ( r8_mach ( 4 ) );
    }

    y = abs ( x );

    if ( xmax < y )
    {
    	result = 0.00;
        return &result;
    }
    /*
    Carefully compute y * (2/M_PI) = (aint(y) + rem(y)) * (.625 + pi2rec)
    = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
    = aint(.625*y) + aint(z) + rem(z)
    */
    ainty = r8_aint ( y );
    yrem = y - ainty;
    prodbg = 0.625 * ainty;
    ainty = r8_aint ( prodbg );
    y = ( prodbg - ainty ) + 0.625 * yrem + pi2rec * y;
    ainty2 = r8_aint ( y );
    ainty += ainty2;
    y -= ainty2;

    ifn = ( dim_typ ) fmod ( ainty, 2.0 );

    if ( ifn == 1 )
    y = 1.00 - y;


    if ( y == 1.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( y <= 0.25 )
    {
        value = y;
        if ( xsml < y )
            value = y * ( 1.5 + r8_csevl ( 32.00 * y * y - 1.00, tancs, nterms ) );
    }
    else if ( y <= 0.5 )
    {
        value = 0.50 * y * ( 1.50 + r8_csevl (8.00 * y * y - 1.00, tancs, nterms ) );
        value = 2.00 * value / ( 1.00 - value * value );
    }
    else
    {
        value = 0.25 * y * ( 1.50 + r8_csevl (2.00 * y * y - 1.00, tancs, nterms ) );
        value = 2.00 * value / ( 1.00 - value * value );
        value = 2.00 * value / ( 1.00 - value * value );
    }

    value = x<0.00 ? -abs(value):abs(value);

    if ( ifn == 1 )
        value *= -1;
	
	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_tanh ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_TANH evaluates the hyperbolic tangent of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_TANH, the hyperbolic tangent of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static dim_typ nterms = 0;
    static ityp sqeps = 0.00;
    static ityp tanhcs[31] =
    {
        -0.25828756643634710438338151450605,
        -0.11836106330053496535383671940204,
        +0.98694426480063988762827307999681E-02,
        -0.83579866234458257836163690398638E-03,
        +0.70904321198943582626778034363413E-04,
        -0.60164243181207040390743479001010E-05,
        +0.51052419080064402965136297723411E-06,
        -0.43320729077584087216545467387192E-07,
        +0.36759990553445306144930076233714E-08,
        -0.31192849612492011117215651480953E-09,
        +0.26468828199718962579377758445381E-10,
        -0.22460239307504140621870997006196E-11,
        +0.19058733768288196054319468396139E-12,
        -0.16172371446432292391330769279701E-13,
        +0.13723136142294289632897761289386E-14,
        -0.11644826870554194634439647293781E-15,
        +0.98812684971669738285540514338133E-17,
        -0.83847933677744865122269229055999E-18,
        +0.71149528869124351310723506176000E-19,
        -0.60374242229442045413288837119999E-20,
        +0.51230825877768084883404663466666E-21,
        -0.43472140157782110106047829333333E-22,
        +0.36888473639031328479423146666666E-23,
        -0.31301874774939399883325439999999E-24,
        +0.26561342006551994468488533333333E-25,
        -0.22538742304145029883494399999999E-26,
        +0.19125347827973995102208000000000E-27,
        -0.16228897096543663117653333333333E-28,
        +0.13771101229854738786986666666666E-29,
        -0.11685527840188950118399999999999E-30,
        +0.99158055384640389120000000000000E-32
    };
    ityp value;
    static ityp xmax = 0.00;
    ityp y;
    ityp yrec;

    if ( nterms == 0 )
    {
        nterms = r8_inits ( tanhcs, 31, 0.10 * r8_mach ( 3 ) );
        sqeps = sqrt ( 3.00 * r8_mach ( 3 ) );
        xmax = - 0.50 * log ( r8_mach ( 3 ) );
    }

    y = abs ( x );

    if ( y <= sqeps )
        value = x;
    else if ( y <= 1.00 )
        value = x * ( 1.00 + r8_csevl ( 2.00 * x * x - 1.00, tanhcs, nterms ) );
    else if ( y <= xmax )
    {
        y = exp ( y );
        yrec = 1.00 / y;
        value = ( y - yrec ) / ( y + yrec );

        if ( x < 0.00 )
            value *= -1;
    }
    else
        value = 1.00 - ((x<0.00)<<1);
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_upak ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_UPAK unpacks an R8 into a mantissa and exponent.
  Discussion:
    This function unpacks a floating point number x so that
      x = y * 2.0^n
    where
      0.5 <= abs ( y ) < 1.0 .
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    C version by John Burkardt.
  Parameters:
    Input, double X, the number to be unpacked.
    Output, double *Y, the mantissa.
    Output, int *N, the exponent.
*/
{
	const itpitpdt * const s_data = data;
	const register ityp x = s_data->a0;
	ityp * y = s_data->a1;
	dim_typ * n = s_data->a2;
	
    ityp absx;

    absx = abs ( x );
    *n = 0;
    *y = 0.00;

    if ( x == 0.00 )
        return NULL;

    while ( absx < 0.50 )
    {
        -- *n;
        absx *= 2.00;
    }

    while ( 1.00 <= absx )
    {
        ++ *n;
        absx *= 0.50;
    }

    *y = absx - ((x<0.00)<<1);
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_uniform_ab_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.
  Discussion:
    Each dimension ranges from A to B.
    This routine implements the recursion
      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )
    The integer arithmetic never requires more than 32 bits,
    including a sign bit.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2005
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
    Input, double A, B, the lower and upper limits of the pseudorandom values.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
*/
{
	const dt2itpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1; 
	const register ityp b = s_data->a2;
	int * seed = s_data->a3;
	
    dim_typ i, k;
    ityp *r;

    if ( *seed == 0 )
        return NULL;

    r = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        k = *seed / 127773;
        *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

        if ( *seed < 0 )
            *seed += i4_huge;

        r[i] = a + ( b - a ) * ( ityp ) ( *seed ) * 4.656612875E-10;
    }
    return r;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sorted_nearest0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SORTED_NEAREST0 returns the nearest element in a sorted R8VEC.
  Discussion:
    An R8VEC is a vector of R8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, double A[N], a sorted vector.
    Input, double VALUE, the value whose nearest vector entry is sought.
    Output, int R8VEC_SORTED_NEAREST0, the index of the nearest
    entry in the vector.
*/
{
	static short result = SHRT_MAX;
	
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp value = s_data->a2;
	
    dim_typ hi;
    dim_typ lo;
    dim_typ mid;

    if ( n < 1 )
    {
    	result = (-1);
        return &result;
    }

    if ( n == 1 )
    {
    	result = 0;
        return &result;
    }

    if ( a[0] < a[n-1] )
    {
        if ( value < a[0] )
        {
        	result = 0;
        	return &result;
        }
        else if ( a[n-1] < value )
        {
        	result = n -1;
        	return &result;
        }
        /*
        Seek an interval containing the value.
        */
        lo = 1;
        hi = n;

        while ( lo < hi - 1 )
        {
            mid = ( lo + hi ) / 2;

            if ( value == a[mid-1] )
            {
            	result = mid -1;
        		return &result;
            }
            else if ( value < a[mid-1] )
                hi = mid;
            else
                lo = mid;
        }
        /*
        Take the nearest.
        */
        if ( abs ( value - a[lo-1] ) < abs ( value - a[hi-1] ) )
        {
        	result = lo -1;
        	return &result;
        }
        else
        {
        	result = hi -1;
            return &result;
        }
    }
    /*
    A descending sorted vector A.
    */
    else
    {
        if ( value < a[n-1] )
        {
        	result = n -1;
            return &result;
        }
        else if ( a[0] < value )
        {
        	result = 0;
            return &result;
        }
        /*
        Seek an interval containing the value.
        */
        lo = n;
        hi = 1;

        while ( lo < hi - 1 )
        {
            mid = ( lo + hi ) / 2;

            if ( value == a[mid-1] )
            {
            	result = mid -1;
            	return &result;
            }
            else if ( value < a[mid-1] )
                hi = mid - 1;
            else
                lo = mid - 1;
        }
        /*
        Take the nearest.
        */
        
        result = abs ( value - a[lo-1] ) < abs ( value - a[hi-1] ) ? lo-1:hi-1;
        return &result;
    }
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _jacobi1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI1 carries out one step of the Jacobi iteration.
  Discussion:
    The linear system A*x=b is to be solved.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N,N], the matrix.
    Input, double B[N], the right hand side.
    Input, double X[N], the current solution estimate.
    Output, double JACOBI1[N], the solution estimate updated by
    one step of the Jacobi iteration.
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x = s_data->a3;
	
    dim_typ i, j;
    ityp * x_new = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        x_new[i] = b[i];
        for ( j = 0; j < n; ++j )
            if ( j != i )
                x_new[i] -= a[i+j*n] * x[j];
        x_new[i] /= a[i+i*n];
    }

    return x_new;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lagrange_basis_function_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_BASIS_FUNCTION_1D evaluates a 1D Lagrange basis function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int MX, the degree of the basis function.
    Input, double XD[MX+1], the interpolation nodes.
    Input, int I, the index of the basis function.
    0 <= I <= MX.
    Input, double XI, the evaluation point.
    Output, double LAGRANGE_BASIS_FUNCTION_1D, the value of the I-th Lagrange 1D
    basis function for the nodes XD, evaluated at XI.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpitit * const s_data = data;
	
	const register dim_typ mx = s_data->a0;
	const register dim_typ i = s_data->a1;	
	ityp * xd = s_data->a2;
	const register ityp xi = s_data->a3;
	
    ityp yi = 1.00;
    if ( xi != xd[i] )
        for (dim_typ j = 0; j < mx + 1; ++j )
            if ( j != i )
                yi *= ( xi - xd[j] ) / ( xd[i] - xd[j] );

	result = yi;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_interp_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_INTERP_2D evaluates the Lagrange interpolant for a product grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int MX, MY, the polynomial degree in X and Y.
    Input, double XD_1D[MX+1], YD_1D[MY+1], the 1D data locations.
    Input, double ZD[(MX+1)*(MY+1)], the 2D array of data values.
    Input, int NI, the number of 2D interpolation points.
    Input, double XI[NI], YI[NI], the 2D interpolation points.
    Output, double LAGRANGE_INTERP_2D[NI], the interpolated values.
*/
{
	const _2dt3pitdt2pit * const s_data = data;
	const register dim_typ mx = s_data->a0;
	const register dim_typ my = s_data->a1;
	ityp * xd_1d = s_data->a2;
	ityp * yd_1d = s_data->a3;
	ityp * zd = s_data->a4;
	const register dim_typ ni = s_data->a5;
	ityp * xi = s_data->a6;
	ityp * yi = s_data->a7;
	
	
    dim_typ i, j, k, l;
    ityp lx;
    ityp ly;
    ityp *zi;

    zi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( k = 0; k < ni; ++k )
    {
        l = 0;
        zi[k] = 0.00;
        for ( j = 0; j < my + 1; ++j)
            for ( i = 0; i < mx + 1; ++i )
            {
                lx = lagrange_basis_function_1d ( mx, xd_1d, i, xi[k] );
                ly = lagrange_basis_function_1d ( my, yd_1d, j, yi[k] );
                zi[k] += zd[l] * lx * ly;
                ++ l;
            }
    }
    return zi;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laguerre_monomial_quadrature ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 February 2008
  Author:
    John Burkardt
  Parameters:
    Input, int EXPON, the exponent.
    Input, int ORDER, the number of points in the rule.
    Input, double W[ORDER], the quadrature weights.
    Input, double X[ORDER], the quadrature points.
    Output, double LAGUERRE_MONOMIAL_QUADRATURE, the quadrature error.
*/
{
	static ityp result = MAX_VAL;
	
	const _3dt2pit * const s_data = data;
	const register dim_typ expon = s_data->a0;
	const register dim_typ order = s_data->a1;
	const register dim_typ option = s_data->a2;
	ityp * w = s_data->a3;
	ityp * x = s_data->a4;
	
    ityp exact;
    dim_typ i;
    ityp quad;
    ityp quad_error;
    /*
    Get the exact value of the integral of the unscaled monomial.
    */
    exact = laguerre_integral ( expon );
    /*
    Evaluate the unweighted monomial at the quadrature points.
    */
    quad = 0.00;

    if ( option == 0 )
        for ( i = 0; i < order; ++i )
            quad +=                  w[i] * pow ( x[i], expon );
    else
        for ( i = 0; i < order; ++i )
            quad +=  exp ( - x[i] ) * w[i] * pow ( x[i], expon );
    /*
    Error:
    */

	result = fabs ( quad - exact ) / exact;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lu_error ( void * data)
/******************************************************************************/
/*
  Purpose:
    LU_ERROR determines the error in an LU factorization.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the matrix.
    Input, double L[N*N], U[N*N], the LU factorization.
    Output, double LU_ERROR, the Frobenius norm
    of the difference matrix A - L * U.
*/
{
	static ityp result = MAX_VAL;
	
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * l = s_data->a2;
	ityp * u = s_data->a3;
	
    ityp *d;
    ityp *lu;
    ityp value;
    lu = r8mat_mm_new ( n, n, n, l, u );
    d = r8mat_sub_new ( n, n, a, lu );
    value = r8mat_norm_fro ( n, n, d );
    free ( d );
    free ( lu );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_linspace2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_LINSPACE2_NEW creates a vector of linearly spaced values.
  Discussion:
    An R8VEC is a vector of R8's.
    4 points evenly spaced between 0 and 12 will yield 2, 4, 6, 8, 10.
    In other words, the interval is divided into N+1 even subintervals,
    and the endpoints of internal intervals are used as the points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, double A, B, the first and last entries.
    Output, double X[N], a vector of linearly spaced data.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp a = s_data->a1;
	ityp b = s_data->a2;
	ityp * x = s_data->a3;
		
    if ( n == 1 )
        x[0] = ( a + b ) / 2.00;
    else
        for (dim_typ i = 0; i < n; ++i)
            x[i] = ( ( ityp ) ( n - i     ) * a+ ( ityp ) (     i + 1 ) * b ) / ( ityp ) ( n     + 1 );

    return NULL;
}

/*********************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_machar ( void * data)
/*********************************************************************/
/*
  Purpose:
    R8_MACHAR computes machine constants for double floating point arithmetic.
  Discussion:
    This routine determines the parameters of the floating-point
    arithmetic system specified below.  The determination of the first
    three uses an extension of an algorithm due to Malcolm,
    incorporating some of the improvements suggested by Gentleman and
    Marovich.
    A FORTRAN version of this routine appeared as ACM algorithm 665.
    This routine is a C translation of the FORTRAN code, and appeared
    as part of ACM algorithm 722.
    An earlier version of this program was published in Cody and Waite.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 April 2006
  Author:
    Original C version by William Cody,
    This C version by John Burkardt.
  Reference:
    William Cody,
    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine
    machine parameters,
    ACM Transactions on Mathematical Software,
    Volume 14, Number 4, pages 303-311, 1988.
    William Cody, W Waite,
    Software Manual for the Elementary Functions,
    Prentice Hall, 1980.
    M Gentleman, S Marovich,
    Communications of the ACM,
    Volume 17, pages 276-277, 1974.
    M. Malcolm,
    Communications of the ACM,
    Volume 15, pages 949-951, 1972.
  Parameters:
    Output, long int* ibeta, the radix for the floating-point representation.
    Output, long int* it, the number of base IBETA digits in the floating-point
    significand.
    Output, long int* irnd:
    0, if floating-point addition chops.
    1, if floating-point addition rounds, but not in the IEEE style.
    2, if floating-point addition rounds in the IEEE style.
    3, if floating-point addition chops, and there is partial underflow.
    4, if floating-point addition rounds, but not in the IEEE style, and
      there is partial underflow.
    5, if floating-point addition rounds in the IEEE style, and there is
      partial underflow.
    Output, long int* ngrd, the number of guard digits for multiplication with
    truncating arithmetic.  It is
    0, if floating-point arithmetic rounds, or if it truncates and only
      IT base IBETA digits participate in the post-normalization shift of the
      floating-point significand in multiplication;
   1, if floating-point arithmetic truncates and more than IT base IBETA
      digits participate in the post-normalization shift of the floating-point
      significand in multiplication.
    Output, long int* MACHEP, the largest negative integer such that
      1.0 + ( double ) IBETA ^ MACHEP != 1.0,
    except that MACHEP is bounded below by - ( IT + 3 ).
    Output, long int* NEGEPS, the largest negative integer such that
      1.0 - ( double ) IBETA ) ^ NEGEPS != 1.0,
    except that NEGEPS is bounded below by - ( IT + 3 ).
    Output, long int* IEXP, the number of bits (decimal places if IBETA = 10)
    reserved for the representation of the exponent (including the bias or
    sign) of a floating-point number.
    Output, long int* MINEXP, the largest in magnitude negative integer such
    that
  ( double ) IBETA ^ MINEXP
    is positive and normalized.
    Output, long int* MAXEXP, the smallest positive power of BETA that overflows.
    Output, double* EPS, the smallest positive floating-point number such
    that
      1.0 + EPS != 1.0.
    in particular, if either IBETA = 2  or IRND = 0,
      EPS = ( double ) IBETA ^ MACHEP.
    Otherwise,
      EPS = ( ( double ) IBETA ^ MACHEP ) / 2.
    Output, double* EPSNEG, a small positive floating-point number such that
      1.0 - EPSNEG != 1.0.
    In particular, if IBETA = 2 or IRND = 0,
      EPSNEG = ( double ) IBETA ^ NEGEPS.
    Otherwise,
      EPSNEG = ( double ) IBETA ^ NEGEPS ) / 2.
    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
    smallest number that can alter 1.0 by subtraction.
    Output, double* XMIN, the smallest non-vanishing normalized floating-point
    power of the radix:
      XMIN = ( double ) IBETA ^ MINEXP
    Output, ityp* XMAX, the largest finite floating-point number.  In
    particular,
      XMAX = ( 1.0 - EPSNEG ) * ( double ) IBETA ^ MAXEXP
    On some machines, the computed value of XMAX will be only the second,
    or perhaps third, largest number, being too small by 1 or 2 units in
    the last digit of the significand.
*/
{
	const _9pli4pit * const s_data = data;
	long int * ibeta = s_data->a0;
	long int * it = s_data->a1;
	long int * irnd = s_data->a2;
	long int * ngrd = s_data->a3;
  	long int * machep = s_data->a4;
	long int * negep = s_data->a5;
	long int * iexp = s_data->a6;
	long int * minexp = s_data->a7;
	long int * maxexp = s_data->a8;
	ityp * eps = s_data->a9;
	ityp * epsneg = s_data->a10;
	ityp * xmin = s_data->a11;
	ityp * xmax = s_data->a12;
	
    ityp a;
    ityp b;
    ityp beta;
    ityp betah;
    ityp betainx;
    int i;
    int itmp;
    int iz;
    int j;
    int k;
    int mx;
    int nxres;
    ityp one;
    ityp t;
    ityp tmp;
    ityp tmp1;
    ityp tmpa;
    ityp two;
    ityp y;
    ityp z;
    ityp zero;

 (*irnd) = 1;
    one = (ityp) (*irnd);
    two = one + one;
    a = two;
    b = a;
    zero = 0.0e0;
    /*
    Determine IBETA and BETA ala Malcolm.
    */
    tmp = ( ( a + one ) - a ) - one;

    while ( tmp == zero )
    {
        a += a;
        tmp = a + one;
        tmp1 = tmp - a;
        tmp = tmp1 - one;
    }

    tmp = a + b;
    itmp = ( int ) ( tmp - a );

    while ( itmp == 0 )
    {
        b += b;
        tmp = a + b;
        itmp = ( int ) ( tmp - a );
    }

    *ibeta = itmp;
    beta = ( ityp ) ( *ibeta );
    /*
    Determine IRND, IT.
    */
 ( *it ) = 0;
    b = one;
    tmp = ( ( b + one ) - b ) - one;

    while ( tmp == zero )
    {
        ++ *it;
        b *= beta;
        tmp = b + one;
        tmp1 = tmp - b;
        tmp = tmp1 - one;
    }

    *irnd = 0;
    betah = beta / two;
    tmp = a + betah;
    tmp1 = tmp - a;

    if ( tmp1 != zero )
        *irnd = 1;

    tmpa = a + beta;
    tmp = tmpa + betah;

    if ( ( *irnd == 0 ) && ( tmp - tmpa != zero ) )
        *irnd = 2;
    /*
    Determine NEGEP, EPSNEG.
    */
 (*negep) = (*it) + 3;
    betainx = one / beta;
    a = one;

    for ( i = 1; i <= (*negep); ++i )
        a *= betainx;

    b = a;
    tmp = ( one - a );
    tmp = tmp - one;

    while ( tmp == zero )
    {
        a *= beta;
        -- *negep;
        tmp1 = one - a;
        tmp = tmp1 - one;
    }

 (*negep) = -(*negep);
 (*epsneg) = a;
    /*
    Determine MACHEP, EPS.
    */
 (*machep) = -(*it) - 3;
    a = b;
    tmp = one + a;

    while ( tmp - one == zero)
    {
        a *= beta;
        ++ *machep;
        tmp = one + a;
    }

    *eps = a;
    /*
    Determine NGRD.
    */
 (*ngrd) = 0;
    tmp = one + *eps;
    tmp = tmp * one;

    if ( ( (*irnd) == 0 ) && ( tmp - one ) != zero )
  (*ngrd) = 1;
    /*
    Determine IEXP, MINEXP and XMIN.

    Loop to determine largest I such that (1/BETA) ** (2**(I))
    does not underflow.  Exit from loop is signaled by an underflow.
    */
    i = nxres = 0;
    k = 1;
    z = betainx;
    t = one + *eps;

    for ( ; ; )
    {
        y = z;
        z = y * y;
        /*
        Check for underflow
        */
        a = z * one;
        tmp = z * t;

        if ( ( a + a == zero ) || ( abs ( z ) > y ) )
            break;

        tmp1 = tmp * betainx;

        if ( tmp1 * beta == z )
            break;

        ++ i;
        k += k;
    }
    /*
    Determine K such that (1/BETA)**K does not underflow.
    First set  K = 2 ** I.
    */
 (*iexp) = i + 1;
    mx = k + k;
    /*
    For decimal machines only.
    */
    if ( *ibeta == 10 )
    {
  (*iexp) = 2;
        iz = *ibeta;
        while ( iz <= k )
        {
            iz *= ( *ibeta );
            ++ (*iexp);
        }
        mx = iz + iz - 1;
    }
    /*
    Loop to determine MINEXP, XMIN.
    Exit from loop is signaled by an underflow.
    */
    for ( ; ; )
    {
  (*xmin) = y;
        y = y * betainx;
        a = y * one;
        tmp = y * t;
        tmp1 = a + a;

        if ( ( tmp1 == zero ) || ( abs ( y ) >= ( *xmin ) ) )
            break;

        ++ k;
        tmp1 = tmp * betainx;
        tmp1 = tmp1 * beta;

        if ( ( tmp1 == y ) && ( tmp != y ) )
        {
            nxres = 3;
            *xmin = y;
            break;
        }
    }

 (*minexp) = -k;
    /*
    Determine MAXEXP, XMAX.
    */
    if ( ( mx <= k + k - 3 ) && ( ( *ibeta ) != 10 ) )
    {
        mx += mx;
        ++ (*iexp);
    }

 (*maxexp) = mx + (*minexp);
    /*
    Adjust IRND to reflect partial underflow.
    */
 (*irnd) += nxres;
    /*
    Adjust for IEEE style machines.
    */
    if ( ( *irnd) >= 2 )
  (*maxexp) = (*maxexp) - 2;
    /*
    Adjust for machines with implicit leading bit in binary
    significand and machines with radix point at extreme
    right of significand.
    */
    i = (*maxexp) + (*minexp);

    if ( ( ( *ibeta ) == 2 ) && ( i == 0 ) )
        -- (*maxexp);

    if ( i > 20 )
        -- (*maxexp);

    if ( a != y )
  (*maxexp) -= 2;

 (*xmax) = one - (*epsneg);
    tmp = (*xmax) * one;

    if ( tmp != (*xmax) )
  (*xmax) = one - beta * (*epsneg);

 (*xmax) /= ( beta * beta * beta * (*xmin) );
    i = (*maxexp) + (*minexp) + 3;

    if ( i > 0 )
    {
        for ( j = 1; j <= i; ++j )
        {
            if ( (*ibeta) == 2 )
          (*xmax) = (*xmax) + (*xmax);
            if ( (*ibeta) != 2 )
          (*xmax) *= beta;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _explode ( void * data)
/******************************************************************************/
/*
  Purpose:
    EXPLODE reports the step when the Mandelbrot iteration at (x,y) "explodes".
  Discussion:
    We assume that the iteration has exploded if an iterate leaves the
    rectangle -2 <= x <= +2, -2 <= y <= +2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 July 2012
  Author:
    John Burkardt
  Parameters:
    Input, double X, Y, the coordinates of the point.

    Input, int COUNT_MAX, the maximum number of steps to consider.
    Output, int EXPLODE, the step on which the iteration 'exploded',
    or 0 if it did not explode in COUNT_MAX steps.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2itdt * const s_data = data;
	const register ityp x = s_data->a0;
	const register ityp y = s_data->a1;
	const register dim_typ count_max = s_data->a2;
	
    dim_typ k;
    int value = 0;
    ityp x1;
    ityp x2;
    ityp y1;
    ityp y2;

    x1 = x;
    y1 = y;

    for ( k = 1; k <= count_max; ++k)
    {
        x2 = x1 * x1 - y1 * y1 + x;
        y2 = 2.00 * x1 * y1 + y;
        if ( x2 < -2.00 || 2.00 < x2 || y2 < -2.00 || 2.0 < y2 )
        {
            value = k;
            break;
        }
        x1 = x2;
        y1 = y2;
    }
    
    result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _grid_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRID_2D returns a regular 2D grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 January 2015
  Author:
    John Burkardt
  Parameters:
    Input, int X_NUM, the number of X values to use.
    Input, double X_LO, X_HI, the range of X values.
    Input, int Y_NUM, the number of Y values to use.
    Input, double Y_LO, Y_HI, the range of Y values.
    Output, double X[X_NUM*Y_NUM], Y[X_NUM*Y_NUM],
    the coordinates of the grid.
*/
{
	const _2dt4it2pit * const s_data = data;
	
	const register dim_typ x_num = s_data->a0;
	const register dim_typ y_num = s_data->a1;
	ityp x_lo = s_data->a2;
	ityp x_hi = s_data->a3;
	ityp y_lo = s_data->a4;
	ityp y_hi = s_data->a5;
	ityp * x = s_data->a6;
	ityp * y = s_data->a7;
	
    dim_typ i, j;
    ityp xi;
    ityp yj;

    if ( x_num == 1 )
        for ( j = 0; j < y_num; ++j )
            for ( i = 0; i < x_num; ++i)
                x[i+j*x_num] = ( x_lo + x_hi ) / 2.00;
    else
        for ( i = 0; i < x_num; ++i )
        {
            xi = ( ( ityp ) ( x_num - i - 1 ) * x_lo+ ( ityp ) (         i     ) * x_hi )/ ( ityp ) ( x_num     - 1 );
            for ( j = 0; j < y_num; ++j )
                x[i+j*x_num] = xi;
        }

    if ( y_num == 1 )
        for ( j = 0; j < y_num; ++j )
            for ( i = 0; i < x_num; ++i )
                y[i+j*x_num] = ( y_lo + y_hi ) / 2.00;
    else
        for ( j = 0; j < y_num; ++j )
        {
            yj = ( ( ityp ) ( y_num - j - 1 ) * y_lo+ ( ityp ) (         j     ) * y_hi )/ ( ityp ) ( y_num     - 1 );
            for ( i = 0; i < x_num; ++i)
                y[i+j*x_num] = yj;
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _polygon_triangulate ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLYGON_TRIANGULATE determines a triangulation of a polygon.
  Discussion:
    There are N-3 triangles in the triangulation.
    For the first N-2 triangles, the first edge listed is always an
    internal diagonal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2014
  Author:
    Original C version by Joseph ORourke.
    This C version by John Burkardt.
  Reference:
    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
  Parameters:
    Input, int N, the number of vertices.
    Input, double X[N], Y[N], the coordinates of each vertex.
    Output, int TRIANGLES[3*(N-2)], the triangles of the
    triangulation.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	
    ityp area;
    int *ear;
    int first;
    int i;
    int i0;
    int i1;
    int i2;
    int i3;
    int i4;
    int *next;
    dim_typ node;
    dim_typ node_m1;
    int *prev;
    int triangle_num;
    int *triangles;
    /*
    We must have at least 3 vertices.
    */
    if ( n < 3 )
        return NULL;
    /*
    Consecutive vertices cannot be equal.
    */
    node_m1 = n - 1;
    for ( node = 0; node < n; ++node )
    {
        if ( x[node_m1] == x[node] && y[node_m1] == y[node] )
            return NULL;
        node_m1 = node;
    }
    /*
    Area must be positive.
    */
    area = 0.00;
    for ( node = 0; node < n - 2; ++node)
        area += 0.50 * (( x[node+1] - x[node] ) * ( y[node+2] - y[node] )- ( x[node+2] - x[node] ) * ( y[node+1] - y[node] ));

    if ( area <= 0.0 )
        return NULL;

    triangles = ( int * ) malloc ( 3 * ( n - 2 ) * sizeof ( int ) );
    /*
    PREV and NEXT point to the previous and next nodes.
    */
    prev = ( int * ) malloc ( n * sizeof ( int ) );
    next = ( int * ) malloc ( n * sizeof ( int ) );

    i = 0;
    prev[i] = n - 1;
    next[i] = i + 1;

    for ( i = 1; i < n - 1; ++i )
    {
        prev[i] = i - 1;
        next[i] = i + 1;
    }

    i = n - 1;
    prev[i] = i - 1;
    next[i] = 0;
    /*
    EAR indicates whether the node and its immediate neighbors form an ear
    that can be sliced off immediately.
    */
    ear = ( int * ) malloc ( n * sizeof ( int ) );
    for ( i = 0; i < n; ++i )
        ear[i] = diagonal ( prev[i], next[i], n, prev, next, x, y );

    triangle_num = 0;

    i2 = 0;

    while ( triangle_num < n - 3 )
    {
        /*
        If I2 is an ear, gather information necessary to carry out
        the slicing operation and subsequent "healing".
        */
        if ( ear[i2] )
        {
            i3 = next[i2];
            i4 = next[i3];
            i1 = prev[i2];
            i0 = prev[i1];
            /*
            Make vertex I2 disappear.
            */
            next[i1] = i3;
            prev[i3] = i1;
            /*
            Update the earity of I1 and I3, because I2 disappeared.
            */
            ear[i1] = diagonal ( i0, i3, n, prev, next, x, y );
            ear[i3] = diagonal ( i1, i4, n, prev, next, x, y );
            /*
            Add the diagonal [I3, I1, I2] to the list.
            */
            triangles[0+triangle_num*3] = i3;
            triangles[1+triangle_num*3] = i1;
            triangles[2+triangle_num*3] = i2;
            triangle_num = triangle_num + 1;
        }
        /*
        Try the next vertex.
        */
        i2 = next[i2];
    }
    /*
    The last triangle is formed from the three remaining vertices.
    */
    i3 = next[i2];
    i1 = prev[i2];

    triangles[0+triangle_num*3] = i3;
    triangles[1+triangle_num*3] = i1;
    triangles[2+triangle_num*3] = i2;
    triangle_num = triangle_num + 1;

    free ( ear );
    free ( next );
    free ( prev );

    return triangles;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _polygon_integral_xx ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLYGON_INTEGRAL_XX integrates the function X*X over a polygon.
  Discussion:
    INTEGRAL = (1/12) * SUM ( I = 1 to N )
  ( X[I]^3 + X[I]^2 * X[I-1] + X[I] * X[I-1]^2 + X[I-1]^3 )
      * ( Y[I] - Y[I-1] )
    where X[N] and Y[N] should be replaced by X[0] and Y[0].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2010
  Author:
    John Burkardt
  Reference:
    SF Bockman,
    Generalizing the Formula for Areas of Polygons to Moments,
    American Mathematical Society Monthly,
    1989, pages 131-132.
  Parameters:
    Input, int N, the number of vertices of the polygon.
    N should be at least 3 for a nonzero result.
    Input, double V[2*N], the coordinates of the vertices.
    Output, double POLYGON_INTEGRAL_XX, the value of the integral.
*/
{
	static ityp _result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * v = s_data->a1;
	
    dim_typ i;
    dim_typ im1;

    ityp result = 0.00;

    if ( n < 3 )
    {
    	_result = MAX_VAL;
        return &_result;
    }

    for ( i = 0; i < n; ++i )
    {
        im1 = i == 0 ? n-1:i-1;
        result += ( v[0+i*2]   * v[0+i*2]   * v[0+i*2] + v[0+i*2]   * v[0+i*2]   * v[0+im1*2] + v[0+i*2]   * v[0+im1*2] * v[0+im1*2] + v[0+im1*2] * v[0+im1*2] * v[0+im1*2] ) * ( v[1+i*2] - v[1+im1*2] );
    }

	_result = result / 12.00;
    return &_result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _polygon_integral_xy ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLYGON_INTEGRAL_XY integrates the function X*Y over a polygon.
  Discussion:
    INTEGRAL = (1/24) * SUM (I=1 to N)
  ( Y[I] *
  ( 3 * X[I]^2 + 2 * X[I] * X[I-1] + X[I-1]^2 )
      + Y[I-1] *
  ( X[I]^2 + 2 * X[I] * X[I-1] + 3 * X[I-1]^2 )
      ) * ( Y[I] - Y[I-1] )
    where X[N] and Y[N] should be replaced by X[0] and Y[0].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2010
  Author:
    John Burkardt
  Reference:
    SF Bockman,
    Generalizing the Formula for Areas of Polygons to Moments,
    American Mathematical Society Monthly,
    1989, pages 131-132.
  Parameters:
    Input, int N, the number of vertices of the polygon.
    N should be at least 3 for a nonzero result.
    Input, double V[2*N], the coordinates of the vertices.
    Output, double POLYGON_INTEGRAL_XY, the value of the integral.
*/
{
	static ityp _result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * v = s_data->a1;
	
    dim_typ i;
    dim_typ im1;
    ityp result = 0.00;

    if ( n < 3 )
    {
    	_result = MAX_VAL;
        return &_result;
    }

    for ( i = 0; i < n; ++i )
    {
        im1 = i == 0 ? n-1:i-1;
        result += (v[1+(i<<1)]   * ( 3.00 *   v[0+(i<<1)]   * v[0+(i<<1)]+ 2.00 *   v[0+(i<<1)]   * v[0+(im1<<1)]+         v[0+(im1<<1)] * v[0+(im1<<1)] )+ v[1+(im1<<1)] * (         v[0+(i<<1)]   * v[0+(i<<1)]+ 2.00 *   v[0+(i<<1)]   * v[0+(im1<<1)]+ 3.00 *   v[0+(im1<<1)] * v[0+(im1<<1)] )) * ( v[1+(i<<1)] - v[1+(im1<<1)] );
    }

	_result = result / 24.00;
    return &_result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _polygon_integral_y ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLYGON_Y integrates the function Y over a polygon.
  Discussion:
    INTEGRAL = (1/6) * SUM ( I = 1 to N )
      - ( Y[I]^2 + Y[I] * Y[I-1] + Y[I-1]^2 ) * ( X[I] - X[I-1] )
    where X[N] and Y[N] should be replaced by X[0] and Y[0].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2010
  Author:
    John
  Reference:
    SF Bockman,
    Generalizing the Formula for Areas of Polygons to Moments,
    American Mathematical Society Monthly,
    1989, pages 131-132.
  Parameters:
    Input, int N, the number of vertices of the polygon.
    N should be at least 3 for a nonzero result.
    Input, double V[2*N], the coordinates of the vertices.
    Output, double POLYGON_INTEGRAL_Y, the value of the integral.
*/
{
	static ityp _result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * v = s_data->a1;
	
    dim_typ i;
    dim_typ im1;
    ityp result = 0.00;

    if ( n < 3 )
    {
    	_result = MAX_VAL;
        return &_result;
    }

    for ( i = 0; i < n; ++i )
    {
        im1 = i == 0 ? n-1:i-1;
        result -= ( v[1+(i<<1)]   * v[1+(i<<1)]+ v[1+(i<<1)]   * v[1+(im1<<1)]+ v[1+(im1<<1)] * v[1+(im1<<1)] )* ( v[0+(i<<1)] - v[0+(im1<<1)] );
    }

	_result = result / 6.00;
    return &_result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _polygon_integral_yy ( void * data) 
/******************************************************************************/
/*
  Purpose:
    POLYGON_INTEGRAL_YY integrates the function Y*Y over a polygon.
  Discussion:
    INTEGRAL = (1/12) * SUM ( I = 1 to N )
      - ( Y[I]^3 + Y[I]^2 * Y[I-1] + Y[I] * Y[I-1]^2 + Y[I-1]^3 )
      * ( X[I] - X[I-1] )
    where X[N] and Y[N] should be replaced by X[0] and Y[0].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2010
  Author:
    John Burkardt
  Reference:
    SF Bockman,
    Generalizing the Formula for Areas of Polygons to Moments,
    American Mathematical Society Monthly,
    1989, pages 131-132.
  Parameters:
    Input, int N, the number of vertices of the polygon.
    N should be at least 3 for a nonzero result.
    Input, double V[2*N], the coordinates of the vertices.
    Output, double POLYGON_INTEGRAL_YY, the value of the integral.

*/
{
	static ityp _result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * v = s_data->a1;
	
    dim_typ i;
    dim_typ im1;
    ityp result = 0.00;

    if ( n < 3 )
    {
    	_result = MAX_VAL;
        return &_result;
    }

    for ( i = 0; i < n; ++i )
    {
        im1 = i == 0 ? n-1:i-1;
        result -= ( v[1+(i<<1)]   * v[1+(i<<1)]   * v[1+(i<<1)]+ v[1+(i<<1)]   * v[1+(i<<1)]   * v[1+(im1<<1)]+ v[1+(i<<1)]   * v[1+(im1<<1)] * v[1+(im1<<1)]+ v[1+(im1<<1)] * v[1+(im1<<1)] * v[1+(im1<<1)] )* ( v[0+(i<<1)] - v[0+(im1<<1)] );
    }

	_result = result / 12.00;
    return &_result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _diaedg ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIAEDG chooses a diagonal edge.
  Discussion:
    The routine determines whether 0--2 or 1--3 is the diagonal edge
    that should be chosen, based on the circumcircle criterion, where
 (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
    quadrilateral in counterclockwise order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 August 2003
  Author:
    Original FORTRAN77 version by Barry Joe.
    C version by John Burkardt.
  Reference:
    Barry Joe,
    GEOMPACK - a software package for the generation of meshes
    using geometric algorithms,
    Advances in Engineering Software,
    Volume 13, pages 325-331, 1991.
  Parameters:
    Input, double X0, Y0, X1, Y1, X2, Y2, X3, Y3, the coordinates of the
    vertices of a quadrilateral, given in counter clockwise order.
    Output, int DIAEDG, chooses a diagonal:
    +1, if diagonal edge 02 is chosen;
    -1, if diagonal edge 13 is chosen;
     0, if the four vertices are cocircular.
*/
{
	static short result = SHRT_MAX;
	
	ityp * const a_data = data;
	ityp x0 = a_data[0];
	ityp y0 = a_data[1];
	ityp x1 = a_data[2];
	ityp y1 = a_data[3];
	ityp x2 = a_data[4];
	ityp y2 = a_data[5];
	ityp x3 = a_data[6];
	ityp y3 = a_data[7];
	
    ityp ca;
    ityp cb;
    ityp dx10;
    ityp dx12;
    ityp dx30;
    ityp dx32;
    ityp dy10;
    ityp dy12;
    ityp dy30;
    ityp dy32;
    ityp s;
    ityp tol;
    ityp tola;
    ityp tolb;
    short value;

    tol = 100.00 * r8_epsilon ( );

    dx10 = x1 - x0;
    dy10 = y1 - y0;
    dx12 = x1 - x2;
    dy12 = y1 - y2;
    dx30 = x3 - x0;
    dy30 = y3 - y0;
    dx32 = x3 - x2;
    dy32 = y3 - y2;

    tola = tol * MAX ( fabs ( dx10 ),
    MAX ( fabs ( dy10 ),
    MAX ( fabs ( dx30 ), fabs ( dy30 ) ) ) );

    tolb = tol * MAX ( fabs ( dx12 ),
    MAX ( fabs ( dy12 ),
    MAX ( fabs ( dx32 ), fabs ( dy32 ) ) ) );

    ca = dx10 * dx30 + dy10 * dy30;
    cb = dx12 * dx32 + dy12 * dy32;

    if ( tola < ca && tolb < cb )
        value = -1;
    else if ( ca < -tola && cb < -tolb )
        value = 1;
    else
    {
        tola = MAX ( tola, tolb );s = ( dx10 * dy30 - dx30 * dy10 ) * cb+ ( dx32 * dy12 - dx12 * dy32 ) * ca;

        if ( tola < s )
            value = -1;
        else if ( s < -tola )
            value = 1;
        else
            value = 0;
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lrline ( void * data)
/******************************************************************************/
/*
  Purpose:
    LRLINE determines where a point lies in relation to a directed line.
  Discussion:
    LRLINE determines whether a point is to the left of, right of,
    or on a directed line parallel to a line through given points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 August 2003
  Author:
    Original FORTRAN77 version by Barry Joe.
    C version by John Burkardt.
  Reference:
    Barry Joe,
    GEOMPACK - a software package for the generation of meshes
    using geometric algorithms,
    Advances in Engineering Software,
    Volume 13, pages 325-331, 1991.
  Parameters:
    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
    directed line is parallel to and at signed distance DV to the left of
    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
    which the position relative to the directed line is to be determined.
    Input, double DV, the signed distance, positive for left.
    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
    to the right of, on, or left of the directed line.  LRLINE is 0 if
    the line degenerates to a point.
*/
{
	static short result = SHRT_MAX;
	
	ityp * const a_data = data;
	ityp xu = a_data[0];
	ityp yu = a_data[1];
	ityp xv1 = a_data[2];
	ityp yv1 = a_data[3];
	ityp xv2 = a_data[4];
	ityp yv2 = a_data[5];
	ityp dv = a_data[6];
	
    ityp dx;
    ityp dxu;
    ityp dy;
    ityp dyu;
    ityp t;
    ityp tol = 0.0000001;
    ityp tolabs;

    dx = xv2 - xv1;
    dy = yv2 - yv1;
    dxu = xu - xv1;
    dyu = yu - yv1;

    tolabs = tol * MAX ( fabs ( dx ),
    MAX ( fabs ( dy ),
    MAX ( fabs ( dxu ),
    MAX ( fabs ( dyu ), fabs ( dv ) ) ) ) );

    t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy );

    if ( tolabs < t )
    {
    	result = 1;	
        return &result;
    }
    else if ( -tolabs <= t )
    {
    	result = 0;
        return &result;
    }

	result = -1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pwl_interp_2d_scattered_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    PWL_INTERP_2D_SCATTERED_VALUE evaluates a 2d interpolant of scattered data
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data points.
    Input, double XYD[2*ND], the data point coordinates.
    Input, double ZD[ND], the data values.
    Input, int T_NUM, the number of triangles.
    Input, int T[3*T_NUM], the triangle information.
    Input, int T_NEIGHBOR[3*T_NUM], the triangle neighbors.
    Input, int NI, the number of interpolation points.
    Input, double XYI[2*NI], the interpolation point coordinates.
    Output, double PWL_INTERP_2D_SCATTERED_VALUE[NI], the interpolated values.
*/
{
	const dt2pitdt2pidtpit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xyd = s_data->a1;
	ityp * zd = s_data->a2;
	const register dim_typ t_num = s_data->a3;
	int * t = s_data->a4;
	int * t_neighbor = s_data->a5;
	const register dim_typ ni = s_data->a6;
	ityp * xyi = s_data->a7;
	
    ityp alpha;
    ityp beta;
    ityp gamma;
    int edge;
    int i, j;
    int step_num;
    ityp *zi;

    zi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( i = 0; i < ni; ++i)
    {
        triangulation_search_delaunay ( nd, xyd, 3, t_num, t, t_neighbor,
        xyi+(i<<1), &j, &alpha, &beta, &gamma, &edge, &step_num );

        if ( j == -1 )
            zi[i] = -1.0;

        zi[i] = alpha * zd[t[0+j*3]]+ beta  * zd[t[1+j*3]]+ gamma * zd[t[2+j*3]];
    }
    return zi;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8tris2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8TRIS2 constructs a Delaunay triangulation of 2D vertices.
  Discussion:
    The routine constructs the Delaunay triangulation of a set of 2D vertices
    using an incremental approach and diagonal edge swaps.  Vertices are
    first sorted in lexicographically increasing (X,Y) order, and
    then are inserted one at a time from outside the convex hull.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 January 2004
  Author:
    Original FORTRAN77 version by Barry Joe.
    C version by John Burkardt.
  Reference:
    Barry Joe,
    GEOMPACK - a software package for the generation of meshes
    using geometric algorithms,
    Advances in Engineering Software,
    Volume 13, pages 325-331, 1991.
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Output, int *TRIANGLE_NUM, the number of triangles in the triangulation;
    TRIANGLE_NUM is equal to 2*node_num - NB - 2, where NB is the number
    of boundary vertices.
    Output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up each
    triangle.  The elements are indices of NODE_XY.  The vertices of the
    triangles are in counterclockwise order.
    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list.
    Positive elements are indices of TIL; negative elements are used for links
    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
    where I, J = triangle, edge index; TRIANGLE_NEIGHBOR[I,J] refers to
    the neighbor along edge from vertex J to J+1 (mod 3).
    Output, int R8TRIS2, is 0 for no error.
*/
{
	static bool result = 2;
	
	const dtpit3pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	ityp * node_xy = s_data->a1;
	int * triangle_num = s_data->a2;
	int * triangle_node = s_data->a3;
	int * triangle_neighbor = s_data->a4;
	
    int base;
    ityp cmax;
    int e;
    int error;
    dim_typ i;
    int *indx;
    dim_typ j;
    dim_typ k;
    dim_typ l;
    int ledg;
    int lr;
    int ltri;
    int m;
    int m1;
    int m2;
    dim_typ n;
    int redg;
    int rtri;
    int *stack;
    int t;
    ityp tol = 100.00 * r8_epsilon ( );
    int top;

    stack = ( int * ) malloc ( node_num * sizeof ( int ) );

    /*
    Sort the vertices by increasing (x,y).
    */
    base = 0;
    indx = r82vec_sort_heap_index_a ( node_num, base, node_xy );

    r82vec_permute ( node_num, indx, base, node_xy );
    /*
    Make sure that the nodes are "reasonably" distinct.
    */
    m1 = 1;

    for ( i = 2; i <= node_num; ++i )
    {
        m = m1;
        m1 = i;

        k = -1;

        for ( j = 0; j <= 1; ++j)
        {
            cmax = MAX ( fabs ( node_xy[((m-1)<<1)+j] ),
            fabs ( node_xy[2*(m1-1)+j] ) );

            if ( tol * ( cmax + 1.0 )< fabs ( node_xy[2*(m-1)+j] - node_xy[2*(m1-1)+j] ) )
            {
                k = j;
                break;
            }
        }

        if ( k == -1 )
        {
        	result = true;
            return &result;
        }

    }
    /*
    Starting from nodes M1 and M2, search for a third point M that
    makes a "healthy" triangle (M1,M2,M)
    */
    m1 = 1;
    m2 = 2;
    j = 3;

    for ( ; ; )
    {
        if ( node_num < j )
        {
        	result = true;
            return &result;
        }

        m = j;

        lr = lrline ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1],
        node_xy[2*(m1-1)+0], node_xy[2*(m1-1)+1],
        node_xy[2*(m2-1)+0], node_xy[2*(m2-1)+1], 0.0 );

        if ( lr != 0 )
            break;

        ++ j;

    }
    /*
    Set up the triangle information for (M1,M2,M), and for any other
    triangles you created because nodes were collinear with M1, M2.
    */
    *triangle_num = j - 2;

    if ( lr == -1 )
    {
        triangle_node[3*0+0] = m1;
        triangle_node[3*0+1] = m2;
        triangle_node[3*0+2] = m;
        triangle_neighbor[3*0+2] = -3;

        for ( i = 2; i <= *triangle_num; ++i )
        {
            m1 = m2;
            m2 = i+1;
            triangle_node[3*(i-1)+0] = m1;
            triangle_node[3*(i-1)+1] = m2;
            triangle_node[3*(i-1)+2] = m;
            triangle_neighbor[3*(i-1)+0] = -3 * i;
            triangle_neighbor[3*(i-1)+1] = i;
            triangle_neighbor[3*(i-1)+2] = i - 1;
        }

        triangle_neighbor[3*(*triangle_num-1)+0] = -3 * (*triangle_num) - 1;
        triangle_neighbor[3*(*triangle_num-1)+1] = -5;
        ledg = 2;
        ltri = *triangle_num;
    }
    else
    {
        triangle_node[3*0+0] = m2;
        triangle_node[3*0+1] = m1;
        triangle_node[3*0+2] = m;
        triangle_neighbor[3*0+0] = -4;

        for ( i = 2; i <= *triangle_num; ++i )
        {
            m1 = m2;
            m2 = i+1;
            triangle_node[3*(i-1)+0] = m2;
            triangle_node[3*(i-1)+1] = m1;
            triangle_node[3*(i-1)+2] = m;
            triangle_neighbor[3*(i-2)+2] = i;
            triangle_neighbor[3*(i-1)+0] = -3 * i - 3;
            triangle_neighbor[3*(i-1)+1] = i - 1;
        }

        triangle_neighbor[3*(*triangle_num-1)+2] = -3 * (*triangle_num);
        triangle_neighbor[3*0+1] = -3 * (*triangle_num) - 2;
        ledg = 2;
        ltri = 1;

    }
    /*
    Insert the vertices one at a time from outside the convex hull,
    determine visible boundary edges, and apply diagonal edge swaps until
    Delaunay triangulation of vertices (so far) is obtained.
    */
    top = 0;

    for ( i = j+1; i <= node_num; ++i )
    {
        m = i;
        m1 = triangle_node[3*(ltri-1)+ledg-1];
        m2 = ledg<=2 ? triangle_node[3*(ltri-1)+ledg] : triangle_node[3*(ltri-1)+0];
        lr = lrline ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1],
        node_xy[2*(m1-1)+0], node_xy[2*(m1-1)+1],
        node_xy[2*(m2-1)+0], node_xy[2*(m2-1)+1], 0.0 );

        if ( 0 < lr )
        {
            rtri = ltri;
            redg = ledg;
            ltri = 0;
        }
        else
        {
            l = -triangle_neighbor[3*(ltri-1)+ledg-1];
            rtri = l / 3;
            redg = (l % 3) + 1;
        }

        vbedg ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1], node_num,
        node_xy, *triangle_num, triangle_node, triangle_neighbor,
        &ltri, &ledg, &rtri, &redg );

        n = *triangle_num + 1;
        l = -triangle_neighbor[3*(ltri-1)+ledg-1];

        for ( ; ; )
        {
            t = l / 3;
            e = ( l % 3 ) + 1;
            l = -triangle_neighbor[3*(t-1)+e-1];
            m2 = triangle_node[3*(t-1)+e-1];
            m1 = e<=2 ? triangle_node[3*(t-1)+e] : triangle_node[3*(t-1)+0];

            ++ *triangle_num;
            triangle_neighbor[3*(t-1)+e-1] = *triangle_num;
            triangle_node[3*(*triangle_num-1)+0] = m1;
            triangle_node[3*(*triangle_num-1)+1] = m2;
            triangle_node[3*(*triangle_num-1)+2] = m;
            triangle_neighbor[3*(*triangle_num-1)+0] = t;
            triangle_neighbor[3*(*triangle_num-1)+1] = *triangle_num - 1;
            triangle_neighbor[3*(*triangle_num-1)+2] = *triangle_num + 1;
            ++ top;

            if ( node_num < top )
            {
            	result = UCHAR_MAX;
                return &result;
            }

            stack[top-1] = *triangle_num;

            if ( t == rtri && e == redg )
                break;
        }

        triangle_neighbor[3*(ltri-1)+ledg-1] = -3 * n - 1;
        triangle_neighbor[3*(n-1)+1] = -3 * (*triangle_num) - 2;
        triangle_neighbor[3*(*triangle_num-1)+2] = -l;
        ltri = n;
        ledg = 2;

        error = swapec ( m, &top, &ltri, &ledg, node_num, node_xy, *triangle_num,
        triangle_node, triangle_neighbor, stack );

        if ( error != 0 )
        {
        	result = true;
            return &result;
        }
    }
    /*
    Undo the sorting.
    */
    for ( i = 0; i < 3; ++i)
        for ( j = 0; j < *triangle_num; ++j )
            triangle_node[i+j*3] = indx [ triangle_node[i+j*3] - 1 ];

    perm_inverse ( node_num, indx );

    r82vec_permute ( node_num, indx, base, node_xy );

    free ( indx );
    free ( stack );

	result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _swapec ( void * data)
/******************************************************************************/
/*
  Purpose:
    SWAPEC swaps diagonal edges until all triangles are Delaunay.
  Discussion:
    The routine swaps diagonal edges in a 2D triangulation, based on
    the empty circumcircle criterion, until all triangles are Delaunay,
    given that I is the index of the new vertex added to the triangulation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 September 2003
  Author:
    Original FORTRAN77 version by Barry Joe.
    C version by John Burkardt.
  Reference:
    Barry Joe,
    GEOMPACK - a software package for the generation of meshes
    using geometric algorithms,
    Advances in Engineering Software,
    Volume 13, pages 325-331, 1991.
  Parameters:
    Input, int I, the index of the new vertex.
    Input/output, int *TOP, the index of the top of the stack.
    On output, TOP is zero.
    Input/output, int *BTRI, *BEDG; on input, if positive, are the
    triangle and edge indices of a boundary edge whose updated indices
    must be recorded.  On output, these may be updated because of swaps.
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input/output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the triangle incidence
    list.  May be updated on output because of swaps.
    Input/output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor
    list; negative values are used for links of the counter-clockwise linked
    list of boundary edges;  May be updated on output because of swaps.
      LINK = -(3*I + J-1) where I, J = triangle, edge index.
    Workspace, int STACK[MAXST]; on input, entries 1 through TOP
    contain the indices of initial triangles (involving vertex I)
    put in stack; the edges opposite I should be in interior;  entries
    TOP+1 through MAXST are used as a stack.
    Output, int SWAPEC, is set to 8 for abnormal return.
*/
{
	static short result = SHRT_MAX;
	
	const dt3pidtpitdt3pi * const s_data = data;
	const register dim_typ i = s_data->a0;
	int * top = s_data->a1;
	int * btri = s_data->a2;
	int * bedg = s_data->a3;
	const register dim_typ node_num = s_data->a4;
	ityp * node_xy = s_data->a5;
	const register dim_typ triangle_num = s_data->a6;
	int * triangle_node = s_data->a7;
	int * triangle_neighbor = s_data->a8;
	int * stack = s_data->a9;
	
    int a;
    int b;
    int c;
    int e;
    int ee;
    int em1;
    int ep1;
    int f;
    int fm1;
    int fp1;
    int l;
    int r;
    int s;
    int swap;
    int t;
    int tt;
    int u;
    ityp x;
    ityp y;
    /*
    Determine whether triangles in stack are Delaunay, and swap
    diagonal edge of convex quadrilateral if not.
    */
    x = node_xy[2*(i-1)+0];
    y = node_xy[2*(i-1)+1];

    for ( ; ; )
    {
        if ( *top <= 0 )
            break;
        t = stack[(*top)-1];
        -- *top;

        if ( triangle_node[3*(t-1)+0] == i )
        {
            e = 2;
            b = triangle_node[3*(t-1)+2];
        }
        else if ( triangle_node[3*(t-1)+1] == i )
        {
            e = 3;
            b = triangle_node[3*(t-1)+0];
        }
        else
        {
            e = 1;
            b = triangle_node[3*(t-1)+1];
        }

        a = triangle_node[3*(t-1)+e-1];
        u = triangle_neighbor[3*(t-1)+e-1];

        if ( triangle_neighbor[3*(u-1)+0] == t )
        {
            f = 1;
            c = triangle_node[3*(u-1)+2];
        }
        else if ( triangle_neighbor[3*(u-1)+1] == t )
        {
            f = 2;
            c = triangle_node[3*(u-1)+0];
        }
        else
        {
            f = 3;
            c = triangle_node[3*(u-1)+1];
        }

        swap = diaedg ( x, y,
        node_xy[2*(a-1)+0], node_xy[2*(a-1)+1],
        node_xy[2*(c-1)+0], node_xy[2*(c-1)+1],
        node_xy[2*(b-1)+0], node_xy[2*(b-1)+1] );

        if ( swap == 1 )
        {
            em1 = i4_wrap ( e - 1, 1, 3 );
            ep1 = i4_wrap ( e + 1, 1, 3 );
            fm1 = i4_wrap ( f - 1, 1, 3 );
            fp1 = i4_wrap ( f + 1, 1, 3 );

            triangle_node[3*(t-1)+ep1-1] = c;
            triangle_node[3*(u-1)+fp1-1] = i;
            r = triangle_neighbor[3*(t-1)+ep1-1];
            s = triangle_neighbor[3*(u-1)+fp1-1];
            triangle_neighbor[3*(t-1)+ep1-1] = u;
            triangle_neighbor[3*(u-1)+fp1-1] = t;
            triangle_neighbor[3*(t-1)+e-1] = s;
            triangle_neighbor[3*(u-1)+f-1] = r;

            if ( 0 < triangle_neighbor[3*(u-1)+fm1-1] )
            {
                *top = *top + 1;
                stack[(*top)-1] = u;
            }

            if ( 0 < s )
            {
                if ( triangle_neighbor[3*(s-1)+0] == u )
                    triangle_neighbor[3*(s-1)+0] = t;
                else if ( triangle_neighbor[3*(s-1)+1] == u )
                    triangle_neighbor[3*(s-1)+1] = t;
                else
                    triangle_neighbor[3*(s-1)+2] = t;

                ++ *top ;

                if ( node_num < *top )
                {
                	result = 8;
                    return &result;
                }
                stack[(*top)-1] = t;
            }
            else
            {
                if ( u == *btri && fp1 == *bedg )
                {
                    *btri = t;
                    *bedg = e;
                }

                l = - ( 3 * t + e - 1 );
                tt = t;
                ee = em1;

                while ( 0 < triangle_neighbor[3*(tt-1)+ee-1] )
                {
                    tt = triangle_neighbor[3*(tt-1)+ee-1];
                    if ( triangle_node[3*(tt-1)+0] == a )
                        ee = 3;
                    else if ( triangle_node[3*(tt-1)+1] == a )
                        ee = 1;
                    else
                        ee = 2;
                }

                triangle_neighbor[3*(tt-1)+ee-1] = l;

            }

            if ( 0 < r )
            {
                if ( triangle_neighbor[3*(r-1)+0] == t )
                    triangle_neighbor[3*(r-1)+0] = u;
                else if ( triangle_neighbor[3*(r-1)+1] == t )
                    triangle_neighbor[3*(r-1)+1] = u;
                else
                    triangle_neighbor[3*(r-1)+2] = u;
            }
            else
            {
                if ( t == *btri && ep1 == *bedg )
                {
                    *btri = u;
                    *bedg = f;
                }

                l = - ( 3 * u + f - 1 );
                tt = u;
                ee = fm1;

                while ( 0 < triangle_neighbor[3*(tt-1)+ee-1] )
                {
                    tt = triangle_neighbor[3*(tt-1)+ee-1];

                    if ( triangle_node[3*(tt-1)+0] == b )
                        ee = 3;
                    else if ( triangle_node[3*(tt-1)+1] == b )
                        ee = 1;
                    else
                        ee = 2;
                }

                triangle_neighbor[3*(tt-1)+ee-1] = l;

            }
        }
    }
	
	result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _vbedg ( void * data)
/******************************************************************************/
/*
  Purpose:
    VBEDG determines which boundary edges are visible to a point.
  Discussion:
    The point (X,Y) is assumed to be outside the convex hull of the
    region covered by the 2D triangulation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 September 2008
  Author:
    Original FORTRAN77 version by Barry Joe.
    C version by John Burkardt.
  Reference:
    Barry Joe,
    GEOMPACK - a software package for the generation of meshes
    using geometric algorithms,
    Advances in Engineering Software,
    Volume 13, pages 325-331, 1991.
  Parameters:
    Input, double X, Y, the coordinates of a point outside the convex hull
    of the current triangulation.
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the triangle incidence list.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list;
    negative values are used for links of a counter clockwise linked list
    of boundary edges;
      LINK = -(3*I + J-1) where I, J = triangle, edge index.
    Input/output, int *LTRI, *LEDG.  If LTRI != 0 then these values are
    assumed to be already computed and are not changed, else they are updated.
    On output, LTRI is the index of boundary triangle to the left of the
    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
    edge of triangle LTRI to the left of the leftmost boundary edge visible
    from (X,Y).  1 <= LEDG <= 3.
    Input/output, int *RTRI.  On input, the index of the boundary triangle
    to begin the search at.  On output, the index of the rightmost boundary
    triangle visible from (X,Y).
    Input/output, int *REDG, the edge of triangle RTRI that is visible
    from (X,Y).  1 <= REDG <= 3.
*/
{
	const _2itdtpitdt6pi * const s_data = data;
	ityp x = s_data->a0;
	ityp y = s_data->a1;
	const register dim_typ node_num = s_data->a2;
	ityp * node_xy = s_data->a3;
	const register dim_typ triangle_num = s_data->a4;
	int * triangle_node = s_data->a5;
	int * triangle_neighbor = s_data->a6;
	int * ltri = s_data->a7;
	int * ledg = s_data->a8;
	int * rtri = s_data->a9;
	int * redg = s_data->a10;
	
	
    int a;
    ityp ax;
    ityp ay;
    int b;
    ityp bx;
    ityp by;
    int done;
    int e;
    int l;
    int lr;
    int t;
    /*
    Find the rightmost visible boundary edge using links, then possibly
    leftmost visible boundary edge using triangle neighbor information.
    */
    if ( *ltri == 0 )
    {
        done = 0;
        *ltri = *rtri;
        *ledg = *redg;
    }
    else
        done = 1;

    for ( ; ; )
    {
        l = -triangle_neighbor[3*((*rtri)-1)+(*redg)-1];
        t = l / 3;
        e = 1 + l % 3;
        a = triangle_node[3*(t-1)+e-1];

        b = triangle_node[e<=2 ? 3*(t-1)+e : 3*(t-1)+0];

        ax = node_xy[((a-1)<<1)+0];
        ay = node_xy[((a-1)<<1)+1];

        bx = node_xy[((b-1)<<1)+0];
        by = node_xy[((b-1)<<1)+1];

        lr = lrline ( x, y, ax, ay, bx, by, 0.00 );

        if ( lr <= 0 )
            break;

        *rtri = t;
        *redg = e;

    }

    if ( done )
        return NULL;

    t = *ltri;
    e = *ledg;

    for ( ; ; )
    {
        b = triangle_node[3*(t-1)+e-1];
        e = i4_wrap ( e-1, 1, 3 );

        while ( 0 < triangle_neighbor[3*(t-1)+e-1] )
        {
            t = triangle_neighbor[3*(t-1)+e-1];

            if ( triangle_node[3*(t-1)+0] == b )
                e = 3;
            else if ( triangle_node[3*(t-1)+1] == b )
                e = 1;
            else
                e = 2;

        }

        a = triangle_node[3*(t-1)+e-1];
        ax = node_xy[2*(a-1)+0];
        ay = node_xy[2*(a-1)+1];

        bx = node_xy[2*(b-1)+0];
        by = node_xy[2*(b-1)+1];

        lr = lrline ( x, y, ax, ay, bx, by, 0.00 );

        if ( lr <= 0 )
            break;

    }

    *ltri = t;
    *ledg = e;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_swtb ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SWTB computes a "slow" backward wavelet transform of an R8VEC.
  Discussion:
    This function inverts the D4 Daubechies wavelet.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 June 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input/output, double S[(N+1)/2], D[(N+1)/2], the transformed data.
    On output, S and D have been overwritten.
    Output, double R8VEC_SWTB[N], the original data sequence.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * s = s_data->a1;
	ityp * d = s_data->a2;
	
    dim_typ i;
    int im1;
    int ip1;
    int n2;
    int nh;
    int np1h;
    ityp *x;
    ityp *y;

    n2 = n + ( n % 2 ) == 1;

    np1h = ( n + 1 ) / 2;
    nh = n / 2;

    for ( i = 0; i < np1h; ++i )
    {
        d[i] /= ( ( sqrt ( 3.00 ) + 1.00 ) / sqrt ( 2.00 ) );
        s[i] /= ( ( sqrt ( 3.00 ) - 1.00 ) / sqrt ( 2.00 ) );
    }

    for ( i = 0; i < np1h; ++i )
    {
        ip1 = i4_wrap ( i + 1, 0, np1h - 1 );
        s[i] += d[ip1];
    }

    y = ( ityp * ) malloc ( n2 * sizeof ( ityp ) );

    for ( i = 0; i < np1h; ++i )
    {
        im1 = i4_wrap ( i - 1, 0, np1h - 1 );
        y[(i<<1)+1] = d[i] + sqrt ( 3.00 ) / 4.00 * s[i]+ ( sqrt ( 3.00 ) - 2.00 ) / 4.00 * s[im1];
    }

    for ( i = 0; i < np1h; ++i )
        y[i<<1] = s[i] - sqrt ( 3.00 ) * y[(i<<1)+1];

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        x[i] = y[i];

    free ( y );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_swtf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SWTF computes a "slow" forward wavelet transform of an R8VEC.
  Discussion:
    This function applies the D4 Daubechies wavelet.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 June 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, double X[N], the data sequence.
    Output, double S[(N+1)/2], D[(N+1)/2], the transformed data.
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * s = s_data->a2;
	ityp * d = s_data->a3;
	
    dim_typ i;
    int im1;
    int ip1;
    int n2;
    int np1h;
    ityp *y;

    n2 = n + ( n % 2 ) == 1;
    y = ( ityp * ) malloc ( n2 * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        y[i] = x[i];

    if ( n < n2 )
        y[n] = 0.00;

    np1h = ( n + 1 ) / 2;

    for ( i = 0; i < np1h; ++i )
        s[i] = y[i<<1] + sqrt ( 3.00 ) * y[(i<<1)+1];

    for ( i = 0; i < np1h; ++i)
    {
        im1 = i4_wrap ( i - 1, 0, np1h - 1 );
        d[i] = y[(i<<1)+1] - sqrt ( 3.00 ) / 4.00 * s[i]- ( sqrt ( 3.00 ) - 2.00 ) / 4.00 * s[im1];
    }

    for ( i = 0; i < np1h; ++i )
    {
        ip1 = i4_wrap ( i + 1, 0, np1h - 1 );
        s[i] -= d[ip1];
    }

    for ( i = 0; i < np1h; ++i )
    {
        s[i] *= ( sqrt ( 3.00 ) - 1.00 ) / sqrt ( 2.00 );
        d[i] *= ( sqrt ( 3.00 ) + 1.00 ) / sqrt ( 2.00 );
    }

    free ( y );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sparse_grid_cfn_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPARSE_GRID_CFN_SIZE sizes a sparse grid using Closed Fully Nested rules.
  Discussion:
    The grid is defined as the sum of the product rules whose LEVEL
    satisfies:
      0 <= LEVEL <= LEVEL_MAX.
    This calculation is much faster than a previous method.  It simply
    computes the number of new points that are added at each level in the
    1D rule, and then counts the new points at a given DIM_NUM dimensional
    level vector as the product of the new points added in each dimension.
    This approach will work for nested families, and may be extensible
    to other families, and to mixed rules.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Output, int SPARSE_GRID_CFN_SIZE, the number of points in the grid.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ dim_num = a_data[0];
	const register dim_typ level_max = a_data[1];
	
    dim_typ dim;
    int h;
    dim_typ j;
    dim_typ l;
    dim_typ level;
    int *level_1d;
    int more;
    int *new_1d;
    dim_typ point_num;
    int t;
    dim_typ v;
    /*
    Special case.
    */

    if ( level_max == 0 )
    {
    	result = 1;
        return &result;
    }
    /*
    Construct the vector that counts the new points in the 1D rule.
    */
    new_1d = ( int * ) malloc ( ( level_max + 1 ) * sizeof ( int ) );

    new_1d[0] = 1;
    new_1d[1] = 2;

    j = 1;
    for ( l = 2; l <= level_max; ++l )
    {
        j = j<<1;
        new_1d[l] = j;
    }
    /*
    Count the number of points by counting the number of new points
    associated with each level vector.
    */
    level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );

    point_num = 0;

    for ( level = 0; level <= level_max; ++level )
    {
        more = h = t = 0;

        for ( ; ;)
        {
            comp_next ( level, dim_num, level_1d, &more, &h, &t );

            v = 1;
            for ( dim = 0; dim < dim_num; ++dim )
                v *= new_1d[level_1d[dim]];

            point_num += v;

            if ( !more )
                break;
        }
    }
    free ( level_1d );
    free ( new_1d );

	result = point_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _square_monomial ( void * data)
/******************************************************************************/
/*
  Purpose:
    SQUARE_MONOMIAL integrates a monomial over a square in 2D.
  Discussion:
    This routine integrates a monomial of the form
      product ( 1 <= dim <= 2 ) x(dim)^expon(dim)
    where the exponents are nonnegative integers.  Note that
    if the combination 0^0 is encountered, it should be treated
    as 1.
    The integration region is:
      A(1) <= X <= B(1)
      A(2) <= Y <= B(2)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A[2], B[2], the lower and upper limits.
    Input, int EXPON[2], the exponents.
    Output, double SQUARE_MONOMIAL, the integral of the monomial.
*/
{
	static ityp result = MAX_VAL;
	
	const pi2pit * const s_data = data;
	
	int * expon = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    dim_typ i;
    ityp value = 1.00;
    for ( i = 0; i < 2; ++i )
        if ( expon[i] == -1 )
        {
        	result = MAX_VAL;
            return &result;
        }

    #pragma omp parallel for num_threads(2)
    for ( i = 0; i < 2; ++i )
        value *= ( pow ( b[i], expon[i] + 1 ) - pow ( a[i], expon[i] + 1 ) ) / ( ityp ) ( expon[i] + 1 );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _square_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    SQUARE_RULE returns a quadrature rule for a square in 2D.
  Discussion:
    The integration region is:
      A(1) <= X <= B(1)
      A(2) <= Y <= B(2)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2014
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Input, double A[2], B[2], the lower and upper limits.
    Input, int ORDER_1D[2], the order of the rule in
    each dimension.  1 <= ORDER_1D(I) <= 5.
    Output, double W[ORDER_1D[0]*ORDER_1D[1]], the weights.
    Output, double XY[2*ORDER_1D[0]*ORDER_1D[1]], the abscissas.
*/
{
	const _2pitpi2pit * const s_data = data; 
	ityp * a = s_data->a0;
	ityp * b = s_data->a1;
	int * order_1d = s_data->a2;
	ityp * w = s_data->a3;
	ityp * xy = s_data->a4;
	
    dim_typ i, j, o;
    dim_typ order;
    ityp *w_1d;
    ityp *x_1d;

    order = order_1d[0] * order_1d[1];

    for ( i = 0; i < 2; ++i )
    {
        o = order_1d[i];

        w_1d = ( ityp * ) malloc ( o * sizeof ( ityp ) );
        x_1d = ( ityp * ) malloc ( o * sizeof ( ityp ) );

        if ( o == 1 )
            line_unit_o01 ( w_1d, x_1d );
        else if ( o == 2 )
            line_unit_o02 ( w_1d, x_1d );
        else if ( o == 3 )
            line_unit_o03 ( w_1d, x_1d );
        else if ( o == 4 )
            line_unit_o04 ( w_1d, x_1d );
        else if ( o == 5 )
            line_unit_o05 ( w_1d, x_1d );
        else
            return NULL;
        /*
        Transform from [-1,+1] to [Ai,Bi]
        */
        for ( j = 0; j < o; ++j )
        {
            w_1d[j] *= ( b[i] - a[i] ) / 2.00;
            x_1d[j] = ( ( 1.00 - x_1d[j] ) * a[i]   + ( 1.00 + x_1d[j] ) * b[i] ) /   2.00;
        }
        /*
        Add this information to the rule.
        */
        r8vec_direct_product ( i, o, x_1d, 2, order, xy );
        r8vec_direct_product2 ( i, o, w_1d, 2, order, w );

        free ( w_1d );
        free ( x_1d );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _square_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    SQUARE_VOLUME: volume of a unit quadrilateral.
  Discussion:
    The integration region is:
      A(1) <= X <= B(1)
      A(2) <= Y <= B(2)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A[2], B[2], the lower and upper limits.
    Output, double SQUARE_VOLUME, the volume.
*/
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	
	result = ( b[0] - a[0] ) * ( b[1] - a[1] );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _angle_shift_deg ( void * data)
/******************************************************************************/
/*
  Purpose:
    ANGLE_SHIFT_DEG shifts angle ALPHA to lie between BETA and BETA+360.
  Discussion:
    The input angle ALPHA is shifted by multiples of 360 to lie
    between BETA and BETA+360.
    The resulting angle GAMMA has all the same trigonometric function
    values as ALPHA.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 June 2007
  Author:
    John Burkardt
  Parameters:
    Input, double ALPHA, the angle to be shifted.
    Input, double BETA, defines the lower endpoint of
    the angle range.
    Output, double ANGLE_SHIFT, the shifted angle.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp alpha = a_data[0];
	const register ityp beta = a_data[1];
	
	result = alpha < beta ? beta - fmod ( beta - alpha, 360.00 ) + 360.00 :  beta + fmod ( alpha - beta, 360.00 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _angle_to_rgb ( void * data)
/******************************************************************************/
/*
  Purpose:
    ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2006
  Author:
    John Burkardt
  Parameters:
    Input, double ANGLE, the angle in the color hexagon.  The sextants are
    defined by the following points:
        0 degrees, 1, 0, 0, red;
       60 degrees, 1, 1, 0, yellow;
      120 degrees, 0, 1, 0, green;
      180 degrees, 0, 1, 1, cyan;
      240 degrees, 0, 0, 1, blue;
      300 degrees, 1, 0, 1, magenta.
    Output, double ANGLE_TO_RGB[3], the RGB specifications for the
    color that lies at the given angle, on the perimeter of the
    color hexagon.  One value will be 1, and one value will be 0.
*/
{
	register ityp angle = *(ityp *) data;
	
    # define DEGREES_TO_RADIANS ( M_PI / 180.0 )

    ityp angle2;
    ityp *rgb;

    rgb = ( ityp * ) malloc ( 3 * sizeof ( ityp ) );
    angle = r8_modp ( angle, 360.00 );

    if ( angle <= 60.00)
    {
        angle2 = DEGREES_TO_RADIANS * 3.00 * angle / 4.00;
        rgb[0] = 1.00;
        rgb[1] = tan ( angle2 );
        rgb[2] = 0.00;
    }
    else if ( angle <= 120.00 )
    {
        angle2 = DEGREES_TO_RADIANS * 3.00 * angle / 4.00;
        rgb[0] = cos ( angle2 ) / sin ( angle2 );
        rgb[1] = 1.00;
        rgb[2] = 0.00;
    }
    else if ( angle <= 180.0 )
    {
        angle2 = DEGREES_TO_RADIANS * 3.00 * ( angle - 120.00 ) / 4.00;
        rgb[0] = 0.00;
        rgb[1] = 1.00;
        rgb[2] = tan ( angle2 );
    }
    else if ( angle <= 240.00 )
    {
        angle2 = DEGREES_TO_RADIANS * 3.00 * ( angle - 120.00 ) / 4.00;
        rgb[0] = 0.00;
        rgb[1] = cos ( angle2 ) / sin ( angle2 );
        rgb[2] = 1.00;
    }
    else if ( angle <= 300.00 )
    {
        angle2 = DEGREES_TO_RADIANS * 3.00 * ( angle - 240.00 ) / 4.00;
        rgb[0] = tan ( angle2 );
        rgb[1] = 0.00;
        rgb[2] = 1.00;
    }
    else if ( angle <= 360.00 )
    {
        angle2 = DEGREES_TO_RADIANS * 3.00 * ( angle - 240.00 ) / 4.00;
        rgb[0] = 1.00;
        rgb[1] = 0.00;
        rgb[2] = cos ( angle2 ) / sin ( angle2 );
    }

    return rgb;
    # undef DEGREES_TO_RADIANS
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _axis_limits ( void * data)
/******************************************************************************/
/*
  Purpose:
    AXIS_LIMITS returns "nice" axis limits for a plot.
  Discussion:
    The routine is given information about the range of a variable, and
    the number of divisions desired.  It returns suggestions for
    labeling a plotting axis for the variable, including the
    starting and ending points, the length of a single division,
    and a suggested tick marking for the axis.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 March 2004
  Author:
    John Burkardt
  Parameters:
    Input, double XMIN, XMAX, the lower and upper values that must be
    included on the axis.  XMIN must be less than XMAX.
    Input, int NDIVS, the number of divisions desired along
    the axis.
    Output, double *PXMIN, *PXMAX, the recommended lower and upper axis
    bounds.  It will be the case that PXMIN <= XMIN < XMAX <= PXMAX.
    Output, double *PXDIV, the recommended size of a single division.
    Output, int *NTICKS, a suggested number of ticks to use,
    if subdividing each of the NDIVS divisions of the axis.
*/
{
	const _2itdt3pitpdt * const s_data = data;
	register ityp xmin = s_data->a0;
 	register ityp xmax = s_data->a1;
 	register dim_typ ndivs = s_data->a2;
 	ityp * pxmin = s_data->a3;
 	ityp * pxmax = s_data->a4;
 	ityp * pxdiv = s_data->a5;
 	dim_typ * nticks = s_data->a6;
	
    # define NSTEPS 5

    ityp best;
    ityp good;
    dim_typ i;
    dim_typ ihi;
    dim_typ ilo;
    int intlog;
    dim_typ iticks[NSTEPS] =
    {
        5, 4, 4, 5, 5
    };
    int ival;
    int j;
    ityp pxmax2;
    ityp pxmin2;
    ityp pxdiv2;
    ityp reldif;
    ityp steps[NSTEPS] =
    {
        1.00,  2.00,  4.0,  5.00, 10.00
    };
    ityp temp;

    if ( xmin == xmax )
    {
        xmin -= 0.50;
        xmax += 0.50;
    }
    else if ( xmax < xmin )
    {
        temp = xmin;
        xmin = xmax;
        xmax = temp;
    }

    if ( ndivs <= 0 )
        ndivs = 5;
    /*
    Set RELDIF, the size of the X interval divided by the largest X.
    */
    reldif = 0.00 + (xmax != xmin)*( xmax - xmin ) / MAX ( fabs ( xmax ), fabs ( xmin ) );
    /*
    If RELDIF tells us that XMIN and XMAX are extremely close,
    do some simple things.
    */
    if ( reldif < 0.00001 )
    {
        if ( xmax == 0.00 )
            *pxdiv = 1.00;
        else
        {
            intlog = ( int ) ( log10 ( xmax ) );

            if ( intlog < 0 )
                -- intlog;

            *pxdiv = pow ( 10.00, intlog );

            if ( 1.00 < *pxdiv )
                *pxdiv = 1.00;
        }

        *nticks = 5;
        *pxmin = xmax - ( ityp ) ( ndivs / 2 ) * (*pxdiv);
        *pxmax = xmax + ( ityp ) ( ndivs - ( ndivs / 2 ) ) * (*pxdiv);
    }
    /*
    But now handle the more general case, when XMIN and XMAX
    are relatively far apart.
    */
    else
    {
        best = -999.0;
        /*
        On second loop, increase INTLOG by 1.
        */
        for ( j = 1; j <= 2; ++j)
        {
            /*
            Compute INTLOG, roughly the logarithm base 10 of the range
            divided by the number of divisions.
            */
            intlog = ( int ) ( log10 ( ( xmax - xmin ) / ( ityp ) ( ndivs ) ) )+ ( j - 1 );

            if ( xmax - xmin  < ( ityp ) ( ndivs ) )
                -- intlog;
            /*
            Now consider taking 1, 2, 4, 5 or 10 steps of size 10**INTLOG:
            */
            for ( i = 1; i <= NSTEPS; ++i )
            {
                /*
                Compute the size of each step.
                */
                pxdiv2 = steps[i-1] * pow ( 10.00, intlog );
                /*
                Make sure NDIVS steps can reach from XMIN to XMAX, at least.
                */
                if ( xmax <= xmin + ndivs * pxdiv2 )
                {
                    /*
                    Now decide where to start the axis.
                    Start the axis at PXMIN2, to the left of XMIN, and
                    representing a whole number of steps of size PXDIV2.
                    */
                    ival = 0.00 <= xmin ? ( dim_typ ) ( xmin / pxdiv2 ) : ( dim_typ ) ( xmin / pxdiv2 ) - 1;
                    pxmin2 = ival * pxdiv2;
                    /*
                    PXMAX2 is, of course, NDIVS steps above PXMIN2.
                    */
                    pxmax2 = pxmin2 + ndivs * pxdiv2;
                    /*
                    Only consider going on if PXMAX2 is at least XMAX.
                    */
                    if ( xmax <= pxmax2 )
                    {
                        /*
                        Now judge this grid by the relative amount of wasted axis length.
                        */
                        good = ( xmax - xmin ) / ( pxmax2 - pxmin2 );

                        if ( best < good )
                        {
                            best = good;
                            *pxmax = pxmax2;
                            *pxmin = pxmin2;
                            *pxdiv = pxdiv2;
                            *nticks = iticks[i-1];
                        }
                    }
                }
            }
        }
    }
    /*
    If necessary, adjust the locations of PXMIN and PXMAX so that the
    interval is more symmetric in containing XMIN through XMAX.
    */
    for ( ; ; )
    {
        ilo = ( dim_typ ) ( ( xmin - *pxmin ) / (*pxdiv) );
        ihi = ( dim_typ ) ( ( *pxmax - xmax ) / (*pxdiv) );

        if ( ihi < ilo + 2 )
            break;

        *pxmin -= (*pxdiv);
        *pxmax -= (*pxdiv);

    }

    return NULL;
    # undef NSTEPS
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bar_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    BAR_CHECK computes the heck digit for a barcode.
  Formula:
    CHECK = SUM ( I = 1, 11, by 2's ) DIGIT(I)
       + 3 * SUM ( I = 2, 10, by 2's ) DIGIT(I)
    CHECK = MOD ( 10 - MOD ( CHECK, 10 ), 10 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2004
  Author:
    John Burkardt
  Parameters:
    Input, int DIGIT[12], entries 1 through 11 of DIGIT contain
    the digits of the bar code.  Each entry must be between 0 and 9.
    The 12th digit should be the check digit.
    Output, int BAR_CHECK, the correct check digit.  If the bar code
    is correct, then DIGIT(12) should equal BAR_CHECK.
*/
{
	static bool result = 2;
	
	dim_typ * digit = data;
	
	result = ( 10 - ( (( digit[0] + digit[2] + digit[4] + digit[6] + digit[8] + digit[10] )+ 3 * ( digit[1] + digit[3] + digit[5] + digit[7] + digit[9] )) % 10 ) ) % 10;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _bar_code ( void * data)
/******************************************************************************/
/*
  Purpose:
    BAR_CODE constructs the 113 character barcode from 11 digits.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2004
  Author:
    John Burkardt
  Parameters:
    Input/output, int DIGIT(12).
    On input, the first 11 entries of DIGIT contain a code to be
    turned into a barcode.
    On output, the 12-th entry of DIGIT is a check digit.
    Output, char BAR_CODE[114], the 113 character bar code corresponding
    to the digit information.
*/
{
	dim_typ * digit = data;
	
    char *bar;
    bool check;
    char *codel;
    char *coder;
    dim_typ i;

    bar = ( char * ) malloc ( 114 * sizeof ( char ) );
    /*
    9 character quiet zone.
    */
    strcpy ( bar, "000000000" );
    /*
    3 character guard pattern.
    */
    strcpy ( bar+9, "101" );
    /*
    7 character product category.
    */
    codel = bar_digit_code_left ( digit[0] );
    strcpy ( bar+12, codel );
    free ( codel );
    /*
    35 characters contain the 5 digit manufacturer code.
    */
    for ( i = 0; i < 5; i++)
    {
        codel = bar_digit_code_left ( digit[i+1] );
        strcpy ( bar+19+i*7, codel );
        free ( codel );
    }
    /*
    Center guard pattern.
    */
    strcpy ( bar+54, "01010" );
    /*
    35 characters contain the 5 digit product code.
    */
    for ( i = 0; i < 5; ++i )
    {
        coder = bar_digit_code_right ( digit[i+6] );
        strcpy ( bar+59+i*7, coder );
        free ( coder );
    }
    /*
    Compute the check digit.
    */
    check = bar_check ( digit );
    digit[11] = check;

    coder = bar_digit_code_right ( check );
    strcpy ( bar+94, coder );
    /*
    Guard pattern.
    */
    strcpy ( bar+101, "101" );
    /*
    Quiet zone.
    */
    strcpy ( bar+104, "000000000" );
    bar[113] = '\0';

    return bar;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _bar_digit_code_left ( void * data)
/******************************************************************************/
/*
  Purpose:
    BAR_DIGIT_CODE_LEFT returns the 7 character left bar code for a digit.
  Example:
    DIGIT = 3
    CODEL = '0111101'
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2004
  Author:
    John Burkardt
  Parameters:
    Input, int DIGIT, the digit, between 0 and 9.
    Output, char BAR_CODE_DIGIT_LEFT[8], the 7 character left code for the digit.
*/
{
	const register dim_typ digit = *(dim_typ *) data;
	
    char *codel = ( char * ) malloc ( 8 * sizeof ( char ) );

    switch(digit)
    {
        case 0:
            strcpy ( codel, "0001101" );
            break;
        case 1:
            strcpy ( codel, "0011001" );
            break;
        case 2:
            strcpy ( codel, "0010011" );
            break;
        case 3:
            strcpy ( codel, "0111101" );
            break;
        case 4:
            strcpy ( codel, "0100011" );
            break;
        case 5:
            strcpy ( codel, "0110001" );
            break;
        case 6:
            strcpy ( codel, "0101111" );
            break;
        case 7:
            strcpy ( codel, "0111011" );
            break;
        case 8:
            strcpy ( codel, "0110111" );
            break;
        case 9:
            strcpy ( codel, "0001011" );
            break;
        default:
            strcpy ( codel, "???????" );
    }

    return codel;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _bar_digit_code_right ( void * data)
/******************************************************************************/
/*
  Purpose:
    BAR_DIGIT_CODE_RIGHT returns the 7 character right bar code for a digit.
  Example:
    DIGIT = 3
    CODER = "1000010"
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2004
  Author:
    John Burkardt
  Parameters:
    Input, int DIGIT, the digit, between 0 and 9.
    Output, char BAR_DIGIT_CODE_RIGHT[8], the 7 character right code.
*/
{
	const register dim_typ digit = *(dim_typ *) data;
	
    char *coder = ( char * ) malloc ( 8 * sizeof ( char ) );

    switch(digit)
    {
        case 0:
            strcpy ( coder, "1110010" );
            break;
        case 1:
            strcpy ( coder, "1100110" );
            break;
        case 2:
            strcpy ( coder, "1101100" );
            break;
        case 3:
            strcpy ( coder, "1000010" );
            break;
        case 4:
            strcpy ( coder, "1011100" );
            break;
        case 5:
            strcpy ( coder, "1001110" );
            break;
        case 6:
            strcpy ( coder, "1010000" );
            break;
        case 7:
            strcpy ( coder, "1000100" );
            break;
        case 8:
            strcpy ( coder, "1001000" );
            break;
        case 9:
            strcpy ( coder, "1110100" );
            break;
        default:
            strcpy ( coder, "???????" );
    }

    return coder;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bmi_english ( void * data)
/******************************************************************************/
/*
  Purpose:
    BMI_ENGLISH computes the body mass index given English measurements.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2004
  Author:
    John Burkardt
  Parameters:
    Input, double W_LB, the body weight in pounds.
    Input, double H_FT, H_IN, the body height in feet and inches
    Output, double BMI_ENGLISH, the body mass index.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp w_lb = a_data[0]; 
	const register ityp h_ft = a_data[1]; 
	const register ityp h_in = a_data[2]; 
	
	result = bmi_metric ( feet_to_meters ( h_ft + ( h_in / 12.0 ) ), feet_to_meters ( h_ft + ( h_in / 12.00 ) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _bmi_metric ( void * data)
/******************************************************************************/
/*
  Purpose:
    BMI_METRIC computes the body mass index given metric measurements.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2004
  Author:
    John Burkardt
  Parameters:
    Input, double W_KG, the body weight in kilograms.
    Input, double H_M, the body height in meters.
    Output, double BMI_METRIC, the body mass index.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp w_kg = a_data[0]; 
	const register ityp h_m = a_data[1]; 
	
	result = w_kg / ( h_m * h_m );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _euler_constant ( void * data)
/******************************************************************************/
/*
  Purpose:
    EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
  Discussion:
    The Euler-Mascheroni constant is often denoted by a lower-case
    Gamma.  Gamma is defined as
      Gamma = limit ( M -> Infinity ) ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2005
  Author:
    John Burkardt
  Parameters:
    Output, double EULER_CONSTANT, the value of the
    Euler-Mascheroni constant.
*/
{
	static ityp result = 0.5772156649015328;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _feet_to_meters ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEET_TO_METERS converts a measurement in feet to meters.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2004
  Author:
    John Burkardt
  Parameters:
    Input, double FT, the length in feet.
    Output, double FEET_TO_METERS, the corresponding length in meters.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp ft = *(ityp *) data; 
	
	result = 0.0254 * 12.0 * ft;
  	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gauss_sum ( void * data)
/******************************************************************************/
/*
  Purpose:
    GAUSS_SUM evaluates a function that is the sum of Gaussians.
  Discussion:
    Gauss_Sum(X) = Sum ( 1 <= J <= Ngauss ) Amplitude(I) * exp ( -Arg )
    where
      Arg = sum ( 1 <= I <= NDIM ) ( ( ( X(I) - Center(I,J) ) / Width(J) )^2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int NDIM, the spatial dimension.
    Input, int N, the number of component Gaussian functions.
    Input, double AMPLITUDE[N], CENTER[NDIM*N], WIDTH[N],
    the amplitude, center and width for the component Gaussian functions.
    Input, double X[NDIM], the point at which the function
    is to be evaluated.
    Output, double GAUSS_SUM, the value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dt4pit * const s_data = data;
	const register dim_typ ndim = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * amplitude = s_data->a2;
	ityp * center = s_data->a3;
	ityp * width = s_data->a4;
	ityp * x = s_data->a5;
	
    ityp arg;
    dim_typ i, j;
    ityp value = 0.00;

    for ( j = 0; j < n; ++j )
    {
        arg = 0.00;
        for ( i = 0; i < ndim; ++i )
            arg += pow ( ( x[i] - center[i+j*ndim] ) / width[j], 2 );
        value += amplitude[j] * exp ( -arg );
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _grid1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRID1 finds grid points between X1 and X2 in N dimensions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the dimension of the points X1 and X2.
    Input, int NSTEP, the number of points to be generated.
    NSTEP must be at least 2.
    Input, double X1[DIM_NUM], X2[DIM_NUM], the first and last
    points, between which the equally spaced points are
    to be computed.
    Output, double X[DIM_NUM*NSTEP], the set of equally spaced
    points.  Each column of X represents one point, with X[*,1] = X1
    and X[*,NSTEP] = X2.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ nstep = s_data->a1;
	ityp * x1 = s_data->a2;
	ityp * x2 = s_data->a3;
	
    dim_typ i, j;
    ityp *x;

    if ( dim_num < 1 || nstep < 2 )
        return NULL;

    x = ( ityp * ) malloc ( dim_num * nstep * sizeof ( ityp ) );

    for ( j = 1; j <= nstep; ++j )
        for ( i = 0; i < dim_num; ++i )
            x[i+(j-1)*dim_num] =( ( ityp ) ( nstep - j     ) * x1[i]+ ( ityp ) (         j - 1 ) * x2[i] )/ ( ityp ) ( nstep     - 1 );
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _grid1n ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRID1N finds the I-th grid point between X1 and X2 in N dimensions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int J, the number of the desired point.
    Normally J would be between 1 and NSTEP, but that is
    not necessary.  Note that J = 1 returns X1 and J = NSTEP
    returns X2.
    Input, int DIM_NUM, the dimension of the points X, X1 and X2.
    Input, int NSTEP, this is the number of equally
    spaced points that are between X1 and X2.  NSTEP must
    be at least 2, because X1 and X2 are always included
    in the set of points.
    Input, double X1[DIM_NUM], X2[DIM_NUM], the first and last
    points, between which the equally spaced points lie.
    Output, double GRID1N[DIM_NUM], the J-th grid point between X1
    and X2.
*/
{
	const _3dt2pit * const s_data = data;
	const register dim_typ j = s_data->a0;
	const register dim_typ dim_num = s_data->a1;
	const register dim_typ nstep = s_data->a2;
	ityp * x1 = s_data->a3;
	ityp * x2 = s_data->a4;
	
    dim_typ i;
    ityp *x;

    if ( dim_num < 1 || nstep < 2 )
        return NULL;

    x = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    for ( i = 0; i < dim_num; ++i)
        x[i] = ( ( ityp ) ( nstep - j     ) * x1[i]+ ( ityp ) (         j - 1 ) * x2[i] )/ ( ityp ) ( nstep     - 1 );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _grid2 ( void * data)
/******************************************************************************/
/*
//  Purpose:
//    GRID2 computes grid points between X1 and X2 in N dimensions.
//  Discussion:
//    GRID2 computes grid points between X1 and X2 in N dimensions.
//    However, X1 need not be the first point computed, nor X2 the last.
//    The user must specify the steps on which X1 and X2 are passed
//    through.  These steps may even be outside the range of 1 through NSTEP.
//    We assume that a set of equally spaced points have
//    been drawn on the line through X1 and X2, and that
//    they have been numbered, with X1 labeled J1 and X2
//    labeled J2.  J1 or J2 may be between 1 and NSTEP,
//    in which case X1 or X2 will actually be returned in the
//    X array, but there is no requirement that J1 or J2
//    satisfy this condition.
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    23 July 2010
//  Author:
//    John Burkardt
//  Parameters:
//    Input, int J1, J2.  J1 specifies the step on which
//    X1 would be computed, and similarly for J2.
//    J1 and J2 must be distinct.
//    Input, int DIM_NUM, the dimension of the points X1 and X2.
//    Input, int NSTEP, this is the number of equally
//    spaced points that are to be generated.
//    NSTEP should be at least 1.
//    Input, double X1[DIM_NUM], X2[DIM_NUM], the points that define
//    the line along which the equally spaced points are generated, and
//    which may or may not be included in the set of computed points.
//    Output, double GRID2[DIM_NUM*NSTEP], the set of equally spaced
//    points.  Each column of X represents one point.
//    If 1 <= J1 <= NSTEP, then X(*,J1) = X1, and similarly for J2.
*/
{
	const _4dt2pit * const s_data = data;
	dim_typ j1 = s_data->a0;
	dim_typ j2 = s_data->a1;
	dim_typ dim_num = s_data->a2;
	dim_typ nstep = s_data->a3;
	ityp * x1 = s_data->a4;
	ityp * x2 = s_data->a5;
	
    dim_typ i, j;
    ityp *x;

    if ( dim_num < 1 || j1 == j2 )
        return NULL;

    x = ( ityp * ) malloc ( nstep * dim_num * sizeof ( ityp ) );

    for ( j = 1; j <= nstep; ++j )
        for ( i = 0; i < dim_num; ++i )
            x[i+(j-1)*dim_num] = ( ( ityp ) ( j2 - j      ) * x1[i]+ ( ityp ) (      j - j1 ) * x2[i] )/ ( ityp ) ( j2     - j1 );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _grid2n ( void * data)
/******************************************************************************/
/*
//  Purpose:
//    GRID2N computes one grid point between X1 and X2 in N dimensions.
//  Discussion:
//    However, X1 need not be the first point computed, nor X2 the last.
//    The user must specify the steps on which X1 and X2 are passed through.
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    23 July 2010
//  Author:
//    John Burkardt
//  Parameters:
//    Input, int J, the J coordinate of the desired point.
//    Note that if J = J1, X will be returned as X1, and if
//    J = J2, X will be returned as X2.
//    Input, int J1, J2.  J1 specifies the step on which
//    X1 would be computed, and similarly for J2.  That is,
//    we assume that a set of equally spaced points have
//    been drawn on the line through X1 and X2, and that
//    they have been numbered, with X1 labeled J1 and X2
//    labeled J2.  J1 and J2 must be distinct.
//    Input, int DIM_NUM, the dimension of the points X1 and X2.
//    Input, double X1[DIM_NUM], X2[DIM_NUM], the points that define
//    the line along which the equally spaced points are
//    generated, and which may or may not be included in the
//    set of computed points.
//    Output, double GRID_2N[DIM_NUM].  X(I) is the J-th point from the
//    set of equally spaced points.
*/
{
	const _4dt2pit * const s_data = data;
	dim_typ j = s_data->a0;
	dim_typ j1 = s_data->a1;
	dim_typ j2 = s_data->a2;
	dim_typ dim_num = s_data->a3;
	ityp * x1 = s_data->a4;
	ityp * x2 = s_data->a5;
	
    dim_typ i;
    ityp *x;

    if ( dim_num < 1 || j1 == j2 )
        return NULL;

    x = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    for ( i = 0; i < dim_num; ++i)
        x[i] = ( ( ityp ) ( j2 - j      ) * x1[i]+ ( ityp ) (      j - j1 ) * x2[i] )/ ( ityp ) ( j2     - j1 );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _grid3 ( void * data)
/******************************************************************************/
/*
//  Purpose:
//
//    GRID3 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
//  Discussion:
//    The line between X1 and X2 will have NSTEP1 points generated along
//    it, and the line between X1 and X3 will have NSTEP2 points generated
//    along it.
//    Fixing the second and third indices of X represents one point, with
//    the following special values:
//      X(*,1,1)      = X1
//      X(*,NSTEP1,1) = X2
//      X(*,1,NSTEP2) = X3.
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    23 July 2010
//  Author:
//    John Burkardt
//  Parameters:
//    Input, int DIM_NUM, the dimension of the points X1, X2 and X3.
//    Input, int NSTEP1, NSTEP2.  These are the number of
//    equally spaced points to generate in the first and second
//    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
//    X3 are always included in the set of points.
//    Input, double X1[DIM_NUM], X2[DIM_NUM], X3[DIM_NUM], the points
//    which define three corners of the parallelogram on
//    which the grid will be generated.
//    Output, double GRID3[DIM_NUM*NSTEP1*NSTEP2], the set of equally
//    spaced points.
*/
{
	const _3dt3pit * const s_data = data;
	dim_typ dim_num = s_data->a0;
	dim_typ nstep1 = s_data->a1;
	dim_typ nstep2 = s_data->a2;
	ityp * x1 = s_data->a3;
	ityp * x2 = s_data->a4;
	ityp * x3 = s_data->a5;
	
    dim_typ i, j, k;
    ityp psi1;
    ityp psi2;
    ityp psi3;
    ityp *x;

    if ( dim_num < 1 || nstep1 < 2 || nstep2 < 2 )
        return NULL;

    x = ( ityp * ) malloc ( nstep1 * nstep2 * dim_num * sizeof ( ityp ) );

    for ( j = 1; j <= nstep1; ++j )
    {
        psi2 = ( ityp ) ( j      - 1 )/ ( ityp ) ( nstep1 - 1 );

        for ( k = 1; k <= nstep2; ++k )
        {
            psi3 = ( ityp ) (      k - 1 )/ ( ityp ) ( nstep2 - 1 );
            psi1 = 1.00 - psi2 - psi3;

            for ( i = 0; i < dim_num; ++i )
                x[i+(j-1)*dim_num+(k-1)*dim_num*nstep1] =psi1 * x1[i]+ psi2 * x2[i]+ psi3 * x3[i];
        }
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _grid3n ( void * data)
/******************************************************************************/
/*
//  Purpose:
//    GRID3N computes a parallelogram grid on 3 points in N dimensions.
//  Discussion:
//    The line between X1 and X2 will have NSTEP1
//    points generated along it, and the line between X1 and
//    X3 will have NSTEP2 points generated along it.
//    The following special values are:
//      J       K         X
//      1       1         X1
//      NSTEP1  1         X2
//      1       NSTEP2    X3
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    23 July 2010
//  Author:
//    John Burkardt
//  Parameters:
//    Input, int J, K, the parallelogram coordinates
//    of the point.  J measures steps from X1 to X2, and
//    K measures steps from X1 to X3.  Normally, J would
//    be between 1 and NSTEP1, K between 1 and NSTEP2,
//    but this is not necessary.
//    Input, int DIM_NUM, the dimension of the points X1, X2 and X3.
//    Input, int NSTEP1, NSTEP2.  These are the number of
//    equally spaced points to generate in the first and second
//    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
//    X3 are always included in the set of points.
//    Input, double X1[DIM_NUM], X2[DIM_NUM], X3[DIM_NUM], the points
//    which define three corners of the parallelogram on
//    which the grid will be generated.
//    Output, double GRID3N[DIM_NUM], the point with coordinates (J,K)
//    from the the set of equally spaced points.
*/
{
	const _5dt3pit * const s_data = data;
	dim_typ j = s_data->a0;
	dim_typ k = s_data->a1;
	dim_typ dim_num = s_data->a2;
	dim_typ nstep1 = s_data->a3;
	dim_typ nstep2 = s_data->a4;
	ityp * x1 = s_data->a5;
	ityp * x2 = s_data->a6;
	ityp * x3 = s_data->a7;
	
	
    dim_typ i;
    ityp psi1;
    ityp psi2;
    ityp psi3;
    ityp *x;

    if ( dim_num < 1 || nstep1 < 2 || nstep2 < 2 )
        return NULL;

    x = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    psi2 = ( ityp ) ( j - 1  ) / ( ityp ) ( nstep1 - 1 );
    psi3 = ( ityp ) ( k - 1  ) / ( ityp ) ( nstep2 - 1 );
    psi1 = 1.00 - psi2 - psi3;

    for ( i = 0; i < dim_num; ++i )
        x[i] = psi1 * x1[i] + psi2 * x2[i] + psi3 * x3[i];

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _grid4 ( void * data)
/******************************************************************************/
/*
//  Purpose:
//    GRID4 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
//  Discussion:
//    Unlike GRID3, GRID4 does not necessarily place X1 at the
//    "origin" of the parallelogram, with X2 and X3 set at the
//    extreme J and K coordinates.  Instead, the user is free
//    to specify the J and K coordinates of the points, although
//    they are required to lie on a subparallelogram of the
//    larger one.
//    The line through X1 and X2 will have NSTEP1
//    points generated along it, and the line through X1 and
//    X3 will have NSTEP2 points generated along it.
//    If we imagine that the
//    main parallelogram is drawn first, with coordinate
//    ranges 1 <= J <= NSTEP1 and 1 <= K <= NSTEP2, then
//    these indices determine the (J,K) coordinates of the
//    three points, namely:
//      X1 : (J1,K1)
//      X2 : (J2,K1)
//      X3 : (J1,K2)
//    Of course, we actually start with the points X1, X2,
//    and X3, and they define a parallelogram and a (J,K)
//    coordinate system over the plane containing them.  We
//    then are free to consider the parallelogram defined
//    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
//    which may or may not contain any of the points X1, X2
//    and X3.
//    Assuming that the indices J1, J2, K1 and K2 are "within
//    bounds", the following special values will be computed:
//      X(*,J1,K1) = X1
//      X(*,J2,K1) = X2
//      X(*,J1,K2) = X3.
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    23 July 2010
//  Author:
//    John Burkardt
//  Parameters:
//    Input, int J1, J2, K1, K2, the indices.
//    Input, int DIM_NUM, the dimension of the points X1, X2 and X3.
//    Input, int NSTEP1, NSTEP2.  These are the number of
//    equally spaced points to generate in the first and second
//    directions.  NSTEP1 and NSTEP2 should be at least 1.
//    Input, double X1[DIM_NUM], X2[DIM_NUM], X3[DIM_NUM], the points
//    which define three corners of the parallelogram on
//    which the grid will be generated.
//    Output, double X[DIM_NUM*NSTEP1*NSTEP2], the set of equally
//    spaced points.  Fixing the second and third indices
//    of X represents one point.
*/
{
	const _7dt3pit * const s_data = data;
	dim_typ j1 = s_data->a0;
	dim_typ j2 = s_data->a1;
	dim_typ k1 = s_data->a2;
	dim_typ k2 = s_data->a3;
	dim_typ dim_num = s_data->a4;
	dim_typ nstep1 = s_data->a5;
	dim_typ nstep2 = s_data->a6;
	ityp * x1 = s_data->a7;
	ityp * x2 = s_data->a8;
	ityp * x3 = s_data->a9;
	
    dim_typ i, j, k;
    ityp psi1;
    ityp psi2;
    ityp psi3;
    ityp *x;

    if ( dim_num < 1 || nstep1 < 2 || nstep2 < 2 || k1 == k2 )
        return NULL;

    x = ( ityp * ) malloc ( nstep1 * nstep2 * dim_num * sizeof ( ityp ) );

    for ( j = 1; j <= nstep1; ++j)
    {
        psi2 = ( ityp ) (  j - j1 )/ ( ityp ) ( j2 - j1 );

        for ( k = 1; k <= nstep2; ++k )
        {
            psi3 = ( ityp ) (  k - k1 )/ ( ityp ) ( k2 - k1 );
            psi1 = 1.00 - psi2 - psi3;

            for ( i = 0; i < dim_num; ++i )
                x[i+(j-1)*dim_num+(k-1)*dim_num*nstep1] =psi1 * x1[i]+ psi2 * x2[i]+ psi3 * x3[i];
        }
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   *_grid4n ( void * data)
/******************************************************************************/
/*
//  Purpose:
//    GRID4N computes a single point on a parallelogram grid in N space.
//  Discussion:
//    The computation is identical to that of GRID4, except that
//    only one point at a time is computed.
//    The line through X1 and X2 will have NSTEP1
//    points generated along it, and the line through X1 and
//    X3 will have NSTEP2 points generated along it.
//    The following special values will be computed:
//      J  K  X
//      J1 K1 X1
//      J2 K2 X2
//      J1 K2 X3
//    If we imagine that the main parallelogram is drawn first, with
//    coordinate ranges 1 <= J <= NSTEP1 and 1 <= K <= NSTEP2, then
//    the indices J and K determine the (J,K) coordinates of the
//    three points X1, X2, and X3, namely:
//      X1 : (J1,K1)
//      X2 : (J2,K1)
//      X3 : (J1,K2)
//    Of course, we actually start with the points X1, X2,
//    and X3, and they define a parallelogram and an (J,K)
//    coordinate system over the plane containing them.  We
//    then are free to consider the parallelogram defined
//    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
//    which may or may not contain any of the points X1, X2
//    and X3.
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    23 July 2010
//  Author:
//    John Burkardt
//  Parameters:
//    Input, int J, the J coordinate of the point X.
//    Input, int J1, J2.  See discussion.
//    Input, int K, the K coordinate of the point X.
//    Input, int K1, K2.  See discussion.
//    Input, int DIM_NUM, the dimension of the points X, X1, X2 and X3.
//    Input, int NSTEP1, NSTEP2.  These are the number of
//    equally spaced points generated in the first and second
//    directions.
//    NSTEP1 and NSTEP2 should be at least 1.
//    Input, double X1[DIM_NUM], X2[DIM_NUM], X3[DIM_NUM], the points
//    which define three corners of the parallelogram on
//    which the grid will be generated.
//    Output, double GRID4N[DIM_NUM], the point whose parallelogram
//    coordinates are (J,K).
*/
{
	const _9dt3pit * const s_data = data;
	dim_typ j = s_data->a0;
	dim_typ j1 = s_data->a1;
	dim_typ j2 = s_data->a2;
	dim_typ k = s_data->a3;
	dim_typ k1 = s_data->a4;
	dim_typ k2 = s_data->a5;
	dim_typ dim_num = s_data->a6;
	dim_typ nstep1 = s_data->a7;
	dim_typ nstep2 = s_data->a8;
	ityp * x1 = s_data->a9;
	ityp * x2 = s_data->a10;
	ityp * x3 = s_data->a11;
	
    dim_typ i;
    ityp psi1;
    ityp psi2;
    ityp psi3;
    ityp *x;

    if ( dim_num < 1 || nstep1 < 2 || nstep2 < 2 || k1 == k2 )
        return NULL;

    psi2 = ( ityp ) ( j  - j1 ) / ( ityp ) ( j2 - j1 );
    psi3 = ( ityp ) ( k  - k1 ) / ( ityp ) ( k2 - k1 );
    psi1 = 1.00 - psi2 - psi3;
    x = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    for ( i = 0; i < dim_num; ++i)
        x[i] = psi1 * x1[i] + psi2 * x2[i] + psi3 * x3[i];

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i2_reverse_bytes ( void * data)
/******************************************************************************/
/*
    Purpose:
    I2_REVERSE_BYTES reverses the two bytes in an I2.
    Licensing:
    This code is distributed under the GNU LGPL license.
    Modified:
    12 May 2007
    Author:
    John Burkardt
    Parameters:
    Input, short int X, a value whose bytes are to be reversed.
    Output, short int I2_REVERSE_BYTES, a value with
    bytes in reverse order.
*/
{
	static short result = SHRT_MAX;
	
	const register short x = *(short *) data;
	
    char c;
    union
    {
        short int yshortint;
        char ychar[2];
    } y;

    y.yshortint = x;

    c = y.ychar[0];
    y.ychar[0] = y.ychar[1];
    y.ychar[1] = c;

	result = ( y.yshortint );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pounds_to_kilograms ( void * data)
/******************************************************************************/
/*
  Purpose:
    POUNDS_TO_KILOGRAMS converts a measurement in pounds to kilograms.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 March 2004
  Author:
    John Burkardt
  Parameters:
    Input, double LB, the weight in pounds.
    Output, double POUNDS_TO_KILOGRAMS, the corresponding weight in kilograms.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp lb = *(ityp *) data;
	
	result = 0.4535924 * lb;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _versine_pulse ( void * data)
/******************************************************************************/
/*
  Purpose:
    VERSINE_PULSE adds a versine pulse to a constant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 2010
  Author:
    John Burkardt
  Parameters:
    Input, double T, the current time.
    Input, double TA, the time at which the pulse begins.
    Input, double TB, the time at which the pulse finishes.
    Input, double V1, the constant value.
    Input, double AMP, the amplitude of the pulse.
    Output, double VERSINE_PULSE, the value of the signal at time T.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	ityp t = a_data[0];
	ityp ta = a_data[1];
	ityp tb = a_data[2];
	ityp v1 = a_data[3];
	ityp amp = a_data[4];
	
	result = v1 + (ta <= t && t <= tb)*( 0.50 * amp * ( 1.00 - cos ( M_2TPI * ( t - ta ) / ( tb - ta ) ) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _asm_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    ASM_ENUM returns the number of alternating sign matrices of a given order.
  Discussion:
    N     ASM_NUM
    0       1
    1       1
    2       2
    3       7
    4      42
    5     429
    6    7436
    7  218348
    A direct formula is
      ASM_NUM ( N ) = product ( 0 <= I <= N-1 ) ( 3 * I + 1 )! / ( N + I )!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 May 2008
  Author:
    John Burkardt
  Parameters:
    Input,int N, the order of the matrices.
    Output, int ASM_ENUM, the number of alternating sign matrices of
    order N.
*/
{
	static int result = INT_MAX;
	
	const register int n = *(int *) data;
	
    int *a;
    int asm_num;
    int *b;
    int *c;
    dim_typ i;
    dim_typ nn;
    int value;

    if ( n <= -1 )
    {
    	result = 0;
        return &result;
    }
    /*
    Row 1
    */
    /*
    Row 2
    */
    if ( n == 0 || n == 1 )
    {
    	result = 1;
        return &result;
    }

    value = 0;

    a = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
    b = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
    c = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    b[0] = 2;
    c[0] = 2;
    a[0] = 1;
    a[1] = 1;
    /*
    Row 3 and on.
    */
    for ( nn = 3; nn <= n; ++nn )
    {
        b[nn-2] = nn;
        for ( i = nn - 2; 2 <= i; --i )
            b[i-1] += b[i-2];
        b[0] = 2;

        c[nn-2] = 2;
        for ( i = nn - 2; 2 <= i; --i)
            c[i-1] += c[i-2];
        c[0] = nn;

        for ( i = 2; i <= nn-1; ++i )
            a[0] += a[i-1];

        for ( i = 2; i <= nn; i++ )
            a[i-1] = a[i-2] * c[i-2] / b[i-2];

    }

    asm_num = 0;
    for ( i = 0; i < n; ++i )
        asm_num += a[i];

    free ( a );
    free ( b );
    free ( c );

	result = asm_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _asm_triangle ( void * data)
/******************************************************************************/
/*
  Purpose:
    ASM_TRIANGLE returns a row of the alternating sign matrix triangle.
  Discussion:
    The first seven rows of the triangle are as follows:
          1      2      3      4      5      6     7
    0     1
    1     1      1
    2     2      3      2
    3     7     14     14      7
    4    42    105    135    105     42
    5   429   1287   2002   2002   1287    429
    6  7436  26026  47320  56784  47320  26026  7436
    For a given N, the value of A(J) represents entry A(I,J) of
    the triangular matrix, and gives the number of alternating sign matrices
    of order N in which the (unique) 1 in row 1 occurs in column J.
    Thus, of alternating sign matrices of order 3, there are
    2 with a leading 1 in column 1:
      1 0 0  1 0 0
      0 1 0  0 0 1
      0 0 1  0 1 0
    3 with a leading 1 in column 2, and
      0 1 0  0 1 0  0 1 0
      1 0 0  0 0 1  1-1 1
      0 0 1  1 0 0  0 1 0
    2 with a leading 1 in column 3:
      0 0 1  0 0 1
      1 0 0  0 1 0
      0 1 0  1 0 0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the desired row.
    Output, int A[N+1], the entries of the row.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int *b;
    int *c;
    dim_typ i;
    dim_typ nn;

    /*
    Row 1
    */
    a[0] = 1;

    if ( n == 0)
        return NULL;
    /*
    Row 2
    */
    a[0] = a[1] = 1;

    if ( n == 1)
        return NULL;
    /*
    Row 3 and on.
    */
    b = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
    c = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    b[0] = c[0] = 2;

    for ( nn = 3; nn <= n + 1; ++nn )
    {

        b[nn-2] = nn;
        for ( i = nn - 2; 2 <= i; --i )
            b[i-1] += b[i-2];

        b[0] = 2;

        c[nn-2] = 2;
        for ( i = nn-2; 2 <= i; --i )
            c[i-1] += c[i-2];
        c[0] = nn;

        for ( i = 2; i <= nn-1; ++i )
            a[0] += a[i-1];

        for ( i = 2; i <= nn; ++i)
            a[i-1] = a[i-2] * c[i-2] / b[i-2];

    }

    free ( b );
    free ( c );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _bell ( void * data)
/******************************************************************************/
/*
  Purpose:
    BELL returns the Bell numbers from 0 to N.
  Discussion:
    The Bell number B(N) is the number of restricted growth functions
    on N.
    Note that the Stirling numbers of the second kind, S^m_n, count the
    number of partitions of N objects into M classes, and so it is
    true that
      B(N) = S^1_N + S^2_N + ... + S^N_N.
  Definition:
    The Bell number B(N) is defined as the number of partitions (of
    any size) of a set of N distinguishable objects.
    A partition of a set is a division of the objects of the set into
    subsets.
  Examples:
    There are 15 partitions of a set of 4 objects:
   (1234), (123)(4), (124)(3), (12)(34), (12)(3)(4),
   (134)(2), (13)(24), (13)(2)(4), (14)(23), (1)(234),
   (1)(23)(4), (14)(2)(3), (1)(24)(3), (1)(2)(34), (1)(2)(3)(4)
    and so B(4) = 15.
  First values:
     N         B(N)
     0           1
     1           1
     2           2
     3           5
     4          15
     5          52
     6         203
     7         877
     8        4140
     9       21147
    10      115975
  Recursion:
    B(I) = sum ( 1 <= J <= I ) Binomial ( I-1, J-1 ) * B(I-J)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 May 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of Bell numbers desired.
    Output, int B[N+1], the Bell numbers from 0 to N.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * b = s_data->a1;
	
    dim_typ i, j;

    b[0] = 1;

    for ( i = 1; i <= n; ++i)
    {
        b[i] = 0;
        for ( j = 1; j <= i; ++j )
            b[i] += b[i-j] * i4_choose ( i-1, j-1 );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _change_greedy ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHANGE_GREEDY makes change for a given total using the biggest coins first.
  Discussion:
    The algorithm is simply to use as many of the largest coin first,
    then the next largest, and so on.
    It is assumed that there is always a coin of value 1.  The
    algorithm will otherwise fail!
  Example:
    Total = 17
    COIN_NUM = 3
    COIN_VALUE = (/ 1, 5, 10 /)
    #  CHANGE              COIN_VALUE(CHANGE)
    4  3 2 1 1             10 5 1 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int TOTAL, the total for which change is to be made.
    Input, int COIN_NUM, the number of types of coins.
    Input, int COIN_VALUE[COIN_NUM], the value of each coin.
    The values should be in ascending order, and if they are not,
    they will be sorted.
    Output, int *CHANGE_NUM, the number of coins given in change.
    Output, int CHANGE[TOTAL], the indices of the coins will be
    in entries 1 through CHANGE_NUM.
*/
{
	const _2dt3pi * const s_data = data;
	register dim_typ total = s_data->a0;
	const register dim_typ coin_num = s_data->a1;
	int * coin_value = s_data->a2;
	int * change_num = s_data->a3;
	int * change = s_data->a4;
	
    dim_typ j;
    *change_num = 0;
    /*
    Find the largest coin smaller than the total.
    */
    j = coin_num - 1;

    while ( 0 <= j )
    {
        if ( coin_value[j] <= total )
            break;
        -- j;
    }

    if ( j < 0 )
        return NULL;
    /*
    Subtract the current coin from the total.
    Once that coin is too big, use the next coin.
    */
    while ( 0 < total )
    {
        if ( coin_value[j] <= total )
        {
            total -= coin_value[j];
            change[*change_num] = j;
            ++ *change_num;
        }
        else
        {
            -- j;
            if ( j < 0 )
                break;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _change_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHANGE_NEXT computes the next set of change for a given sum.
  Example:
    Total = 17
    COIN_NUM = 3
    COIN_VALUE = { 1, 5, 10 }
        #  CHANGE              COIN_VALUE(CHANGE)
    1   4  3 2 1 1             10 5 1 1
    2   8  3 1 1 1 1 1 1 1     10 1 1 1 1 1 1 1
    3   5  2 2 2 1 1            5 5 5 1 1
    4   9  2 2 1 1 1 1 1 1 1    5 5 1 1 1 1 1 1 1
    5  13  2 1 1 1 1 1 1 1 1 1  5 1 1 1 1 1 1 1 1 1
           1 1 1                1 1 1
    6  17  1 1 1 1 1 1 1 1 1 1  1 1 1 1 1 1 1 1 1 1 1
           1 1 1 1 1 1 1        1 1 1 1 1 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 June 2004
  Author:
    John Burkardt
  Parameters:
    Input, int TOTAL, the total for which change is to be made.
    Input, int COIN_NUM, the number of types of coins.
    Input, int COIN_VALUE[COIN_NUM], the value of each coin.
    The values must be in ascending order.
    Input/output, int *CHANGE_NUM, the number of coins given in change
    for this form of the change.
    Input/output, int CHANGE[CHANGE_NUM], the indices of the coins.
    The user must dimension this array to have dimension TOTAL!
    Input/output, int *DONE.  The user sets DONE = .TRUE. on
    first call to tell the routine this is the beginning of a computation.
    The program resets DONE to .FALSE. and it stays that way until
    the last possible change combination is made, at which point the
    program sets DONE to TRUE again.
*/
{
	const _2dt3pipb * const s_data = data;
	const register dim_typ total = s_data->a0;
	const register dim_typ coin_num = s_data->a1;
	int * coin_value = s_data->a2;
	int * change_num = s_data->a3;
	int * change = s_data->a4;
	bool * done = s_data->a5;
	
    int change_num2;
    dim_typ coin_num2;
    dim_typ i;
    int last;
    int total2;

    if ( *done )
    {
        /*
        Make sure the coin values are sorted into ascending order.
        */
        if ( !i4vec_ascends ( coin_num, coin_value ) )
            return NULL;
        /*
        Start with the greedy change.
        */
        change_greedy ( total, coin_num, coin_value, change_num, change );
        /*
        In a few cases, like change for 4 cents, we're done after the first call.
        */
        *done = *change_num == total;
        return NULL;

    }
    /*
    Find the last location in the input change which is NOT a penny.
    */
    last = -1;

    for ( i = *change_num-1; 0 <= i; --i )
        if ( change[i] != 0 )
        {
            last = i;
            break;
        }
    /*
    If that location is still -1, an error was made.
    */
    if ( last == -1 )
    {
        *done = 1;
        return NULL;
    }
    /*
    Sum the entries from that point to the end.
    */
    total2 = 0;
    for ( i = last; i <= *change_num-1; ++i )
        total2 +=  coin_value [ change[i] ];
    /*
    Make greedy change for the partial sum using coins smaller than that one.
    */
    coin_num2 = change[last];

    change_greedy ( total2, coin_num2, coin_value, &change_num2, change+last );

    *change_num = last + change_num2;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _chinese_check ( void * data) 
/******************************************************************************/
/*
  Purpose:
    CHINESE_CHECK checks the Chinese remainder moduluses.
  Discussion:
    For a Chinese remainder representation, the moduluses M(I) must
    be positive and pairwise prime.  Also, in case this is not obvious,
    no more than one of the moduluses may be 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of moduluses.
    Input, int M[N], the moduluses.  These should be positive
    and pairwise prime.
    Output, int CHINESE_CHECK, is TRUE if an error was detected.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * m = s_data->a1;
	
    dim_typ i, j;
    /*
    Do not allow nonpositive entries.
    */
    for ( i = 0; i < n; ++i )
        if ( m[i] <= 0 )
        {
        	result = 1;
        	return &result;
        }
    /*
    Allow one entry to be 1, but not two entries.
    */
    for ( i = 0; i < n; ++i)
        for ( j = i+1; j < n; ++j )
            if ( m[i] == 1 && m[j] == 1 )
            {
            	result = 2;
        		return &result;
            }
    /*
    Now check pairwise primeness.
    */
    if ( !i4vec_pairwise_prime ( n, m ) )
    {
    	result = 3;
        return &result;
    }

    result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _congruence ( void * data)
/******************************************************************************/
/*
  Purpose:
    CONGRUENCE solves a congruence of the form A * X = C ( mod B ).
  Discussion:
    A, B and C are given integers.  The equation is solvable if and only
    if the greatest common divisor of A and B also divides C.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 November 2004
  Author:
    John Burkardt
  Reference:
    Eric Weisstein, editor,
    CRC Concise Encylopedia of Mathematics,
    CRC Press, 1998, page 446.
  Parameters:
    Input, int A, B, C, the coefficients of the Diophantine equation.
    Output, int *ERROR, error flag, is 1 if an error occurred..
    Output, int CONGRUENCE, the solution of the Diophantine equation.
    X will be between 0 and B-1.
*/
{
	static int result = INT_MAX;
	
	const _3ipi * const s_data = data;
	const register int a = s_data->a0;
	const register int b = s_data->a1;
	const register int c = s_data->a2;
	int * error = s_data->a3;
	
    # define N_MAX 100

    int a_copy;
    int a_mag;
    int a_sign;
    int b_copy;
    int b_mag;
    int b_sign;
    int c_copy;
    int g;
    int k;
    int n;
    int q[N_MAX];
    int swap;
    int x;
    int y;
    int z;
    /*
    Defaults for output parameters.
    */
    *error = x = y = 0;
    /*
    Special cases.
    */
    if ( a == 0 && b == 0 && c == 0 || a == 0 && b != 0 && c == 0 || a != 0 && b == 0 && c == 0 || a != 0 && b == 0 && c == 0 || a != 0 && b != 0 && c == 0 )
    {
    	result = 0;
        return &result;
    }
    else if ( a == 0 && b == 0 && c != 0 )
    {
        *error = 1;
        result = 0;
        return &result;
    }
    else if ( a == 0 && b != 0 && c != 0 )
    {
        x = 0;
        if ( ( c % b ) != 0 )
            *error = 2;
        result = x;
        return &result;
    }
    else if ( a != 0 && b == 0 && c != 0 )
    {
        x = c / a;
        if ( ( c % a ) != 0 )
            *error = 3;
        result = x;
        return &result;
    }
    
    /*
    Now handle the "general" case: A, B and C are nonzero.

    Step 1: Compute the GCD of A and B, which must also divide C.
    */
    g = i4_gcd ( a, b );

    if ( ( c % g ) != 0 )
    {
        *error = 4;
        result = x;
        return &result;
    }

    a_copy = a / g;
    b_copy = b / g;
    c_copy = c / g;
    /*
    Step 2: Split A and B into sign and magnitude.
    */
    a_mag = abs ( a_copy );
    a_sign = i4_sign ( a_copy );
    b_mag = abs ( b_copy );
    b_sign = i4_sign ( b_copy );
    /*
    Another special case, A_MAG = 1 or B_MAG = 1.
    */
    if ( a_mag == 1 )
    {
    	result = a_sign * c_copy;
        return &result;
    }
    else if ( b_mag == 1 )
    {
    	result = 0;
        return &result;
    }
    /*
    Step 3: Produce the Euclidean remainder sequence.
    */
    if ( b_mag <= a_mag )
    {
        swap = 0;
        q[0] = a_mag;
        q[1] = b_mag;
    }
    else
    {
        swap = 1;
        q[0] = b_mag;
        q[1] = a_mag;
    }

    n = 3;

    for ( ; ; )
    {
        q[n-1] = ( q[n-3] % q[n-2] );

        if ( q[n-1] == 1 )
            break;
        ++ n;

        if ( N_MAX < n )
        {
        	result = INT_MAX;
            return &result;
        }
    }
    /*
    Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
    */
    y = 0;
    for ( k = n; 2 <= k; --k )
    {
        x = y;
        y = ( 1 - x * q[k-2] ) / q[k-1];
    }
    /*
    Step 5: Undo the swapping.
    */
    if ( swap == 1 )
    {
        z = x;
        x = y;
        y = z;
    }
    /*
    Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
    */
    x *= a_sign;
    /*
    Step 7: Multiply by C, so that X * A + Y * B = C.
    */
    x *=  c_copy;
    /*
    Step 8: Now force 0 <= X < B.
    */
    x %= b;
    /*
    Step 9: Force positivity.
    */
    if ( x < 0 )
        x +=  b;

	result = x;
    return &result;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _count_pose_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    COUNT_POSE_RANDOM poses a problem for the game "The Count is Good"
  Discussion:
    The French television show "The Count is Good" has a game that goes
    as follows:
      A number is chosen at random between 100 and 999.  This is the GOAL.
      Six numbers are randomly chosen from the set 1, 2, 3, 4, 5, 6, 7, 8,
      9, 10, 25, 50, 75, 100.  These numbers are the BLOCKS.
      The player must construct a formula, using some or all of the blocks,
   (but not more than once), and the operations of addition, subtraction,
      multiplication and division.  Parentheses should be used to remove
      all ambiguity.  However, it is forbidden to use subtraction in a
      way that produces a negative result, and all division must come out
      exactly, with no remainder.
    This routine poses a sample problem from the show.  The point is,
    to determine how to write a program that can solve such a problem.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 June 2003
  Author:
    John Burkardt
  Reference:
    Raymond Seroul,
    Programming for Mathematicians,
    Springer Verlag, 2000, page 355-357.
  Parameters:
    Input/output, int *SEED, a seed for the random number generator.
    Output, int BLOCKS[6], the six numbers available for the formula.
    Output, int *GOAL, the goal number.
*/
{
	int ** const a_data = data;
	int * seed = a_data[0];
	int * blocks = a_data[1];
	int * goal = a_data[2];
	
    dim_typ i;
    int ind[6];
    const int stuff[14] =
    {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 75, 100
    };

    *goal = i4_uniform_ab ( 100, 999, seed );
    ksub_random ( 14, 6, seed, ind );

    #pragma omp parallel for num_threads(6)
    for ( i = 0; i < 6; ++i)
        blocks[i] = stuff[ind[i]-1];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _debruijn ( void * data)
/******************************************************************************/
/*
  Purpose:
    DEBRUIJN constructs a de Bruijn sequence.
  Discussion:
    Suppose we have an alphabet of M letters, and we are interested in
    all possible strings of length N.  If M = 2 and N = 3, then we are
    interested in the M**N strings:
      000
      001
      010
      011
      100
      101
      110
      111
    Now, instead of making a list like this, we prefer, if possible, to
    write a string of letters, such that every consecutive sequence of
    N letters is one of the strings, and every string occurs once, if
    we allow wraparound.
    For the above example, a suitable sequence would be the 8 characters:
      00011101(00...
    where we have suggested the wraparound feature by repeating the first
    two characters at the end.
    Such a sequence is called a de Bruijn sequence.  It can easily be
    constructed by considering a directed graph, whose nodes are all
    M^(N-1) strings of length N-1.  A node I has a directed edge to
    node J (labeled with character K) if the string at node J can
    be constructed by beheading the string at node I and adding character K.
    In this setting, a de Bruijn sequence is simply an Eulerian circuit
    of the directed graph, with the edge labels being the entries of the
    sequence.  In general, there are many distinct de Bruijn sequences
    for the same parameter M and N.  This program will only find one
    of them.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of letters in the alphabet.
    Input, int N, the number of letters in a codeword.
    Output, int STRING[M**N], a deBruijn string.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * string = s_data->a2;
	
    dim_typ i;
    int iedge;
    int *inode;
    int *ivec;
    dim_typ j;
    int *jnode;
    int *jvec;
    dim_typ k;
    int *knode;
    dim_typ nedge;
    dim_typ nnode;
    int success;
    int *trail;
    /*
    Construct the adjacency information.
    */
    nnode = powi ( m, n-1 );
    nedge = powi ( m, n );

    inode = ( int * ) malloc ( nedge * sizeof ( int ) );
    ivec = ( int * ) malloc ( ( n - 1 ) * sizeof ( int ) );
    jnode = ( int * ) malloc ( nedge * sizeof ( int ) );
    jvec = ( int * ) malloc ( ( n - 1 ) * sizeof ( int ) );
    knode = ( int * ) malloc ( nedge * sizeof ( int ) );

    iedge = 0;

    for ( i = 1; i <= nnode; ++i )
    {
        index_unrank0 ( n-1, m, i, ivec );

        for ( k = 1; k <= m; ++k )
        {
            /*
            Shift N-2 entries of IVEC down.
            */
            for ( j = 0; j < n-2; ++j )
                jvec[j] = ivec[j+1];
            jvec[n-2] = k;

            j = index_rank0 ( n-1, m, jvec );

            inode[iedge] = i;
            jnode[iedge] = j;
            knode[iedge] = k;
            ++ iedge;
        }
    }

    free ( ivec );
    free ( jvec );
    /*
    Determine a circuit.
    */
    trail = ( int * ) malloc ( nedge * sizeof ( int ) );
    digraph_arc_euler ( nnode, nedge, inode, jnode, &success, trail );
    /*
    The string is constructed from the labels of the edges in the circuit.
    */
    for ( i = 0; i < nedge; ++i )
        string[i] = knode[trail[i]-1];

    free ( inode );
    free ( jnode );
    free ( knode );
    free ( trail );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _derange_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    DERANGE_ENUM returns the number of derangements of N objects.
  Discussion:
    A derangement of N objects is a permutation which leaves no object
    unchanged.
    A derangement of N objects is a permutation with no fixed
    points.  If we symbolize the permutation operation by "P",
    then for a derangment, P(I) is never equal to I.
    The number of derangements of N objects is sometimes called
    the subfactorial function, or the derangement number D(N).
    D(N) is the number of ways of placing N non-attacking rooks on
    an N by N chessboard with one diagonal deleted.
    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
    The number of permutations with exactly K items in the right
    place is COMB(N,K) * D(N-K).
    The formula:
      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
    based on the inclusion/exclusion law.
  Recursion:
      D(0) = 1
      D(1) = 0
      D(2) = 1
      D(N) = (N-1) * ( D(N-1) + D(N-2) )
    or
      D(0) = 1
      D(1) = 0
      D(N) = N * D(N-1) + (-1)**N
  First values:
     N         D(N)
     0           1
     1           0
     2           1
     3           2
     4           9
     5          44
     6         265
     7        1854
     8       14833
     9      133496
    10     1334961
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects to be permuted.
    Output, int DERANGE_ENUM, the number of derangements of N objects.
*/
{
	static int result = INT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    int i;
    int value;
    int value1;
    int value2;

    if ( n < 0 )
        value = 0;
    else if ( n == 0 )
        value = 1;
    else if ( n == 1 )
        value = 0;
    else if ( n == 2 )
        value = 1;
    else
    {
        value1 = 0;
        value = 1;

        for ( i = 3; i <= n; ++i )
        {
            value2 = value1;
            value1 = value;
            value = ( i - 1 ) * ( value1 + value2 );
        }
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _derange_enum2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    DERANGE_ENUM2 returns the number of derangements of 0 through N objects.
  Discussion:
    A derangement of N objects is a permutation which leaves no object
    unchanged.
    A derangement of N objects is a permutation with no fixed
    points.  If we symbolize the permutation operation by "P",
    then for a derangment, P(I) is never equal to I.
    The number of derangements of N objects is sometimes called
    the subfactorial function, or the derangement number D(N).
    D(N) is the number of ways of placing N non-attacking rooks on
    an N by N chessboard with one diagonal deleted.
    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
    The number of permutations with exactly K items in the right
    place is COMB(N,K) * D(N-K).
    The formula is:
      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
    based on the inclusion/exclusion law.
  Recursion:
      D(0) = 1
      D(1) = 0
      D(2) = 1
      D(N) = (N-1) * ( D(N-1) + D(N-2) )
    or
      D(0) = 1
      D(1) = 0
      D(N) = N * D(N-1) + (-1)^N
  Example:
     N         D(N)
     0           1
     1           0
     2           1
     3           2
     4           9
     5          44
     6         265
     7        1854
     8       14833
     9      133496
    10     1334961
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the maximum number of objects to be permuted.
    Output, int D[N+1]; D(I) is the number of derangements of
    I objects.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * d = s_data->a1;
	
    d[0] = 1;
    d[1] = 0;

    for (dim_typ i = 2; i <= n; ++i )
        d[i] = ( i - 1 ) * ( d[i-1] + d[i-2] );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _derange_enum3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    DERANGE_ENUM3 returns the number of derangements of 0 through N objects.
  Discussion:
    A derangement of N objects is a permutation which leaves no object
    unchanged.
    A derangement of N objects is a permutation with no fixed
    points.  If we symbolize the permutation operation by "P",
    then for a derangment, P(I) is never equal to I.
    The number of derangements of N objects is sometimes called
    the subfactorial function, or the derangement number D(N).
    D(N) is the number of ways of placing N non-attacking rooks on
    an N by N chessboard with one diagonal deleted.
    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
    The number of permutations with exactly K items in the right
    place is COMB(N,K) * D(N-K).
    The formula is:
      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
    based on the inclusion/exclusion law.
    D(N) = nint ( N! / E )
  Recursion:
      D(0) = 1
      D(1) = 0
      D(2) = 1
      D(N) = (N-1) * ( D(N-1) + D(N-2) )
    or
      D(0) = 1
      D(1) = 0
      D(N) = N * D(N-1) + (-1)^N
  Example:
     N         D(N)
     0           1
     1           0
     2           1
     3           2
     4           9
     5          44
     6         265
     7        1854
     8       14833
     9      133496
    10     1334961
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the maximum number of objects to be permuted.
    Output, int DERANGE_ENUM3, the number of derangements of N objects.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = n==0 ? 1 : n == 1 ? 0 : ( dim_typ ) ( 0.50 + ( r8_factorial ( n ) / M_E ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _derange0_back_candidate ( void * data)
/******************************************************************************/
/*
  Purpose:
    DERANGE0_BACK_CANDIDATE finds possible K-th entries of a derangement.
  Discussion:
    A derangement of N objects is a permutation of (0,...,N-1) which leaves
    no object unchanged.
    A derangement of N objects is a permutation with no fixed
    points.  If we symbolize the permutation operation by "P",
    then for a derangment, P(I) is never equal to I.
    The number of derangements of N objects is sometimes called
    the subfactorial function, or the derangement number D(N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 June 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the derangement.
    Input, int A[N].  The first K-1 entries of A
    record the currently set values of the derangement.
    Input, int K, the entry of the derangement for which candidates
    are to be found.
    Input/output, int *NSTACK, the length of the stack.
    Input/output, int STACK[(N*(N+1))/2].  On output, we have added
    the candidates for entry K to the end of the stack.

    Input/output, int NCAN[N], the number of candidates for each level.
*/
{
	const _2dt4pi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	int * a = s_data->a2;
	int * nstack = s_data->a3;
	int * stack = s_data->a4;
	int * ncan = s_data->a5;
	
    dim_typ ican;
    int *ifree;
    dim_typ nfree;
    /*
    Consider all the integers from 1 through N that have not been used yet.
    */
    nfree = n - k + 1;
    ifree = ( int * ) malloc ( n * sizeof ( int ) );

    perm0_free ( k - 1, a, nfree, ifree );
    /*
    Everything but K is a legitimate candidate for the K-th entry.
    */
    ncan[k-1] = 0;

    for ( ican = 0; ican < nfree; ++ican )

        if ( ifree[ican] != k - 1 )
        {
            ncan[k-1] = ncan[k-1] + 1;
            stack[*nstack] = ifree[ican];
            ++ *nstack;
        }

    free ( ifree );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _derange0_back_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    DERANGE0_BACK_NEXT returns the next derangement of N items.
  Discussion:
    A derangement of N objects is a permutation of (0,...,N-1) which leaves
    no object unchanged.
    A derangement of N objects is a permutation with no fixed
    points.  If we symbolize the permutation operation by "P",
    then for a derangment, P(I) is never equal to I.
    The number of derangements of N objects is sometimes called
    the subfactorial function, or the derangement number D(N).
    This routine uses backtracking.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items to be deranged.  N should be 2 or more.
    Input/output, int A[N].
    On the first call, the input value of A is not important.
    On return with MORE = TRUE, A contains the next derangement.
    On subsequent input, A should not be changed.
    Input/output, int *MORE.
    On first call, set MORE to FALSE and do not alter it after.
    On return, MORE is TRUE if another derangement is being returned in A,
    and FALSE if no more derangements could be found.
*/
{
	const dtpipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	bool * more = s_data->a2;
	
    dim_typ i;
    static int indx = -1;
    static int k = -1;
    static int *ncan = NULL;
    static int *stack = NULL;
    static int stack_max = -1;
    static int stack_num = -1;

    if ( !( *more ) )
    {
        if ( n < 2 )
        {
            *more = 0;
            return NULL;
        }

        indx = k = stack_num = 0;
        stack_max = ( n * ( n + 1 ) ) >> 1;

        if ( stack )
            free ( stack );

        stack = ( int * ) malloc ( stack_max * sizeof ( int ) );

        for ( i = 0; i < stack_max; ++i )
            stack[i] = 0;

        if ( ncan )
            free ( ncan );

        ncan = ( int * ) malloc ( n * sizeof ( int ) );

        for ( i = 0; i < n; ++i )
            ncan[i] = 0;
        *more = true;
    }

    for ( ; ; )
    {
        i4vec_backtrack ( n, stack_max, a, &indx, &k, &stack_num, stack, ncan );

        if ( indx == 1 )
            break;
        else if ( indx == 2 )
            derange0_back_candidate ( n, a, k, &stack_num, stack, ncan );
        else
        {
            *more = false;
            free ( ncan );
            ncan = NULL;
            free ( stack );
            stack = NULL;
            break;
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _derange0_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    DERANGE0_CHECK is TRUE if a permutation is a derangement of (0,...,N-1).
  Discussion:
    A derangement of N objects is a permutation which leaves no object
    unchanged.
    A derangement of N objects is a permutation with no fixed
    points.  If we symbolize the permutation operation by "P",
    then for a derangment, P(I) is never equal to I.
    The number of derangements of N objects is sometimes called
    the subfactorial function, or the derangement number D(N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 June 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects permuted.
    Input, int A[N], a permutation of (0,...,N-1).
    Output, int DERANGE0_CHECK is TRUE if there was an error.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ check;
    dim_typ i, j;
    /*
    Values must be between 0 and N-1.
    */
    for ( i = 0; i < n; ++i )
        if ( a[i] < 0 || n - 1 < a[i] )
        {
        	result = 0;
            return &result;
        }
    /*
    Every value must be represented.
    */
    for ( j = 0; j < n; ++j )
    {
        check = 0;
        for ( i = 0; i < n; ++i )
            if ( a[i] == j )
            {
                check = 1;
                break;
            }
        if ( ! check )
        {
        	result = check;
            return &result;
        }
    }
    /*
    Values must be deranged.
    */
    for ( i = 0; i < n; ++i )
        if ( a[i] == i )
        {
        	result = 0;
            return &result;
        }

    result = 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _derange0_weed_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    DERANGE0_WEED_NEXT computes derangements of (0,...,N-1).
  Discussion:
    A derangement of N objects is a permutation which leaves no object
    unchanged.
    A derangement of N objects is a permutation with no fixed
    points.  If we symbolize the permutation operation by "P",
    then for a derangment, P(I) is never equal to I.
    The number of derangements of N objects is sometimes called
    the subfactorial function, or the derangement number D(N).
    This routine simply generates all permutations, one at a time,
    and weeds out those that are not derangements.
  Example:
    Here are the derangements when N = 4:
    1032
    1230
    1302
    2031
    2301
    2310
    3012
    3201
    3210
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 June 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects being permuted.
    Input/output, int A[N].
    On first call, the input contents of A are unimportant.  But
    on the second and later calls, the input value of A should be
    the output value returned on the previous call.
    On output, A contains the next derangement.
    Input/output, int *MORE.
    Set MORE = FALSE before the first call.
    MORE will be reset to TRUE and a derangement will be returned.
    Each new call produces a new derangement until MORE is returned FALSE.
    Input/output, int *MAXDER, *NUMDER, two parameters
    used by the program for bookkeeping.  The user should declare these
    variables, and pass the output values from one call to the next,
    but should not alter them.
*/
{
	const dtpipb2pdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	bool * more = s_data->a2;
	dim_typ * maxder = s_data->a3;
	dim_typ * numder = s_data->a4;
	
    bool check;
    /*
    Initialization on call with MORE = FALSE.
    */
    if ( !( *more ) )
    {
        *maxder = derange_enum ( n );
        *numder = 0;
    }
    /*
    Watch out for cases where there are no derangements.
    */
    if ( *maxder == 0 )
    {
        *more = 0;
        return NULL;
    }
    /*
    Get the next permutation.
    */
    for ( ; ; )
    {
        perm0_lex_next ( n, a, more );
        /*
        See if it is a derangment.
        */
        check = derange0_check ( n, a );

        if ( check )
            break;
    }

    ++ *numder;

    if ( *maxder <= *numder)
        *more = false;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _digraph_arc_euler ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIGRAPH_ARC_EULER returns an Euler circuit in a digraph.
  Discussion:
    An Euler circuit of a digraph is a path which starts and ends at
    the same node and uses each directed edge exactly once.  A digraph is
    eulerian if it has an Euler circuit.  The problem is to decide whether
    a given digraph is eulerian and to find an Euler circuit if the
    answer is affirmative.
  Method:
    A digraph has an Euler circuit if and only if the number of incoming
    edges is equal to the number of outgoing edges at each node.
    This characterization gives a straightforward procedure to decide whether
    a digraph is eulerian.  Furthermore, an Euler circuit in an eulerian
    digraph G of NEDGE edges can be determined by the following method:
      STEP 1: Choose any node U as the starting node, and traverse any edge
  ( U, V ) incident to node U, and than traverse any unused edge incident
      to node U.  Repeat this process of traversing unused edges until the
      starting node U is reached.  Let P be the resulting walk consisting of
      all used edges.  If all edges of G are in P, than stop.
      STEP 2: Choose any unused edge ( X,  Y) in G such that X is
      in P and Y is not in P.  Use node X as the starting node and
      find another walk Q using all unused edges as in step 1.
      STEP 3: Walk P and walk Q share a common node X, they can be merged
      to form a walk R by starting at any node S of P and to traverse P
      until node X is reached; than, detour and traverse all edges of Q
      until node X is reached and continue to traverse the edges of P until
      the starting node S is reached.  Set P = R.
      STEP 4: Repeat steps 2 and 3 until all edges are used.
    The running time of the algorithm is O ( NEDGE ).
    The digraph is assumed to be connected.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2009
  Author:
    Original FORTRAN77 version by Hang Tong Lau.
    C version by John Burkardt.
  Reference:
    Hang Tong Lau,
    Algorithms on Graphs,
    Tab Books, 1989.
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int NEDGE, the number of edges.
    Input, int INODE[NEDGE], JNODE(NEDGE); the I-th edge starts at node
    INODE(I) and ends at node JNODE(I).
    Output, int *SUCCESS, is TRUE if an Euler circuit was found.
    Output, int TRAIL[NEDGE].  TRAIL[I] is the edge number of the I-th
    edge in the Euler circuit.
*/
{
	const _2dt4pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	const register dim_typ nedge = s_data->a1;
	int * inode = s_data->a2;
	int * jnode = s_data->a3;
	int * success = s_data->a4;
	int * trail = s_data->a5;
	
    int *candid;
    int *endnod;
    int istak;
    dim_typ i, j, k, l;
    dim_typ len;
    dim_typ lensol;
    dim_typ lenstk;
    dim_typ *stack;
    /*
    Check if the digraph is eulerian.
    */
    for ( i = 0; i < nedge; ++i )
        trail[i] = 0;

    endnod = ( int * ) malloc ( nedge * sizeof ( int ) );

    for ( i = 0; i < nedge; i++ )
        endnod[i] = 0;

    for ( i = 1; i <= nedge; ++i )
    {
        j = inode[i-1];
        ++ trail[j-1];
        j = jnode[i-1];
        ++ endnod[j-1];
    }

    for ( i = 1; i <= nnode; ++i)
        if ( trail[i-1] != endnod[i-1] )
        {
            *success = false;
            free ( endnod );
            return NULL;
        }
    /*
    The digraph is eulerian; find an Euler circuit.
    */
    *success = 1;
    lensol = 1;
    lenstk = 0;

    candid = ( int * ) malloc ( nedge * sizeof ( int ) );
    stack = ( dim_typ * ) malloc (  nedge * sizeof ( dim_typ ) << 1 );
    /*
    Find the next edge.
    */
    for ( ; ; )
    {
        if ( lensol == 1 )
        {
            endnod[0] = inode[0];
            stack[0] = stack[1] = 1;

            lenstk = 2;
        }
        else
        {
            l = lensol - 1;

            if ( lensol != 2 )
            endnod[l-1] = inode[trail[l-1]-1] + jnode[trail[l-1]-1] - endnod[l-2];

            k = endnod[l-1];

            for ( i = 1; i <= nedge; ++i )
                candid[i-1] = ( k == jnode[i-1] );

            for ( i = 1; i <= l; ++i )
                candid[trail[i-1]-1] = 0;

            len = lenstk;

            for ( i = 1; i <= nedge; ++i )
                if ( candid[i-1] )
                {
                    ++ len;
                    stack[len-1] = i;
                }

            stack[len] = len - lenstk;
            lenstk = len + 1;
        }

        for ( ; ; )
        {
            istak = stack[lenstk-1];
            lenstk = lenstk - 1;

            if ( istak != 0 )
                break;

            -- lensol;

            if ( lensol == 0 )
            {
                i4vec_reverse ( nedge, trail );

                free ( candid );
                free ( endnod );
                free ( stack );
                return NULL;
            }
        }

        trail[lensol-1] = stack[lenstk-1];
        stack[lenstk-1] = istak - 1;

        if ( lensol == nedge )
            break;
        ++ lensol;
    }

    i4vec_reverse ( nedge, trail );

    free ( candid );
    free ( endnod );
    free ( stack );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _diophantine ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIOPHANTINE solves a Diophantine equation A * X + B * Y = C.
  Discussion:
    Given integers A, B and C, produce X and Y so that
      A * X + B * Y = C.
    In general, the equation is solvable if and only if the
    greatest common divisor of A and B also divides C.
    A solution (X,Y) of the Diophantine equation also gives the solution
    X to the congruence equation:
      A * X = C mod ( B ).
    Generally, if there is one nontrivial solution, there are an infinite
    number of solutions to a Diophantine problem.
    If (X0,Y0) is a solution, then so is ( X0+T*B/D, Y0-T*A/D ) where
    T is any integer, and D is the greatest common divisor of A and B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2009
  Author:
    John Burkardt
  Reference:
    Eric Weisstein, editor,
    CRC Concise Encylopedia of Mathematics,
    CRC Press, 1998, page 446.
  Parameters:
    Input, int A, B, C, the coefficients of the Diophantine equation.
    Output, int *ERROR, is TRUE if an error occurred.
    Output, int *X, *Y, the solution of the Diophantine equation.
    Note that the algorithm will attempt to return a solution with
    smallest Euclidean norm.
*/
{
	const _3i3pi * const s_data = data;
	const register int a = s_data->a0;
	const register int b = s_data->a1;
	const register int c = s_data->a2;
	int * error = s_data->a3;
	int * x = s_data->a4;
	int * y = s_data->a5;
	
    # define N_MAX 100

    int a_copy;
    int a_mag;
    int a_sign;
    int b_copy;
    int b_mag;
    int b_sign;
    int c_copy;
    dim_typ g;
    dim_typ k;
    dim_typ n;
    int q[N_MAX];
    bool swap;
    /*
    Defaults for output parameters.
    */
    *error = *x = *y = 0;
    /*
    Special cases.
    */
    if ( a == 0 && b == 0 && c == 0 )
    {
        *x = *y = 0;
        return NULL;
    }
    else if ( a == 0 && b == 0 && c != 0 )
    {
        *error = 1;
        *x = *y = 0;
        return NULL;
    }
    else if ( a == 0 && b != 0 && c == 0 )
    {
        *x = *y = 0;
        return NULL;
    }
    else if ( a == 0 && b != 0 && c != 0 )
    {
        *x = 0;
        *y = c / b;
        if ( ( c % b ) != 0 )
            *error = 1;
        return NULL;
    }
    else if ( a != 0 && b == 0 && c == 0 )
    {
        *x = *y = 0;
        return NULL;
    }
    else if ( a != 0 && b == 0 && c != 0 )
    {
        *x = c / a;
        *y = 0;
        if ( ( c % a ) != 0 )
        *error = 1;
        return NULL;
    }
    else if ( a != 0 && b != 0 && c == 0 )
    {
        g = i4_gcd ( a, b );
        *x = b / g;
        *y = - a / g;
        return NULL;
    }
    /*
    Now handle the "general" case: A, B and C are nonzero.

    Step 1: Compute the GCD of A and B, which must also divide C.
    */
    g = i4_gcd ( a, b );

    if ( ( c % g ) != 0 )
    {
        *error = 1;
        return NULL;
    }

    a_copy = a / g;
    b_copy = b / g;
    c_copy = c / g;
    /*
    Step 2: Split A and B into sign and magnitude.
    */
    a_mag = abs ( a_copy );
    a_sign = i4_sign ( a_copy );
    b_mag = abs ( b_copy );
    b_sign = i4_sign ( b_copy );
    /*
    Another special case, A_MAG = 1 or B_MAG = 1.
    */
    if ( a_mag == 1 )
    {
        *x = a_sign * c_copy;
        *y = 0;
        return NULL;
    }
    else if ( b_mag == 1 )
    {
        *x = 0;
        *y = b_sign * c_copy;
        return NULL;
    }
    /*
    Step 3: Produce the Euclidean remainder sequence.
    */
    if ( b_mag <= a_mag )
    {
        swap = false;
        q[0] = a_mag;
        q[1] = b_mag;
    }
    else
    {
        swap = true;
        q[0] = b_mag;
        q[1] = a_mag;
    }

    n = 3;

    for ( ; ; )
    {
        q[n-1] = q[n-3] % q[n-2];

        if ( q[n-1] == 1 )
            break;

        ++ n;

        if ( N_MAX < n )
            return NULL;
    }
    /*
    Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
    */
    *y = 0;
    for ( k = n; 2 <= k; --k )
    {
        *x = *y;
        *y = ( 1 - *x * q[k-2] ) / q[k-1];
    }
    /*
    Step 5: Undo the swapping.
    */
    if ( swap )
        i4_swap ( x, y );
    /*
    Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
    */
    *x *= a_sign;
    *y *= b_sign;
    /*
    Step 7: Multiply by C, so that X * A + Y * B = C.
    */
    *x *= c_copy;
    *y *= c_copy;
    /*
    Step 8: Given a solution (X,Y), try to find the solution of
    minimal magnitude.
    */
    diophantine_solution_minimize ( a_copy, b_copy, x, y );

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _diophantine_solution_minimize ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIOPHANTINE_SOLUTION_MINIMIZE seeks a minimal solution of a Diophantine equation.
  Discussion:
    Given a solution (X,Y) of a Diophantine equation:
      A * X + B * Y = C.
    then there are an infinite family of solutions of the form
  ( X(i), Y(i) ) = ( X + i * B, Y - i * A )
    An integral solution of minimal Euclidean norm can be found by
    tentatively moving along the vectors (B,-A) and (-B,A) one step
    at a time.
    When large integer values are input, the real arithmetic used
    is essential.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2009
  Author:
    John Burkardt
  Reference:
    Eric Weisstein, editor,
    CRC Concise Encylopedia of Mathematics,
    CRC Press, 1998, page 446.
  Parameters:
    Input, int A, B, the coefficients of the Diophantine equation.
    A and B are assumed to be relatively prime.
    Input/output, int *X, *Y, on input, a solution of the Diophantine
    equation.  On output, a solution of minimal Euclidean norm.
*/
{
	const _2i2pi * const s_data = data;
	const register int a = s_data->a0;
	const register int b = s_data->a1;
	int * x = s_data->a2;
	int * y = s_data->a3;
	
    ityp fa;
    ityp fb;
    ityp fx;
    ityp fy;
    ityp norm;
    ityp norm_new;
    ityp t;
    int xnew;
    int ynew;
    /*
    Compute the minimum for T real, and then look nearby.
    */
    fa = ( ityp ) a;
    fb = ( ityp ) b;
    fx = ( ityp ) ( *x );
    fy = ( ityp ) ( *y );

    t = ( - fb * fx + fa * fy ) / ( fa * fa + fb * fb );

    *x += r8_nint ( t ) * b;
    *y -= r8_nint ( t ) * a;
    /*
    Now look nearby.
    */
    norm = ( fx * fx + fy * fy );

    for ( ; ; )
    {
        xnew = *x + b;
        ynew = *y - a;

        fx = ( ityp ) xnew;
        fy = ( ityp ) ynew;

        norm_new = ( fx * fx + fy * fy );

        if ( norm <= norm_new )
            break;

        *x = xnew;
        *y = ynew;
        norm = norm_new;
    }

    for ( ; ; )
    {
        xnew = *x - b;
        ynew = *y + a;

        fx = ( ityp ) xnew;
        fy = ( ityp ) ynew;

        norm_new = ( fx * fx + fy * fy );

        if ( norm <= norm_new )
            break;

        *x = xnew;
        *y = ynew;
        norm = norm_new;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _equiv_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    EQUIV_NEXT computes the partitions of a set one at a time.
  Discussion:
    A partition of a set assigns each element to exactly one subset.
    The number of partitions of a set of size N is the Bell number B(N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, number of elements in the set to be partitioned.
    Output, int *NPART, number of subsets in the partition.
    Output, int JARRAY[N].  JARRAY[I] is the number of elements
    in the I-th subset of the partition.
    Output, int IARRAY[N].  IARRAY(I) is the class to which
    element I belongs.
    Input/output, int *MORE.  Set MORE = FALSE before first call.
    It is reset and held at TRUE as long as
    the partition returned is not the last one.
    When MORE is returned FALSE, all the partitions
    have been computed and returned.
*/
{
	const dtpdt2pipb * const s_data = data;
	const register dim_typ n = s_data->a0; 
	dim_typ * npart = s_data->a1;
	int * jarray = s_data->a2;
	int * iarray = s_data->a3;
	bool * more = s_data->a4;
	
    dim_typ i;
    dim_typ l;
    dim_typ m;

    if ( !( *more ) )
    {
        *npart = 1;
        for ( i = 0; i < n; ++i )
            iarray[i] = 1;
        jarray[0] = n;
    }
    else
    {
        m = n;

        while ( jarray[iarray[m-1]-1] == 1 )
        {
            iarray[m-1] = 1;
            -- m;
        }

        l = iarray[m-1];
        *npart = *npart + m - n;
        jarray[0] += n - m;

        if ( l == *npart )
        {
            ++ *npart;
            jarray[*npart-1] = 0;
        }
        iarray[m-1] = l + 1;
        -- jarray[l-1];
        ++ jarray[l];
    }

    *more = ( *npart != n );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _equiv_next2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    EQUIV_NEXT2 computes, one at a time, the partitions of a set.
  Discussion:
    A partition of a set assigns each element to exactly one subset.
    The number of partitions of a set of size N is the Bell number B(N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    John Burkardt.
  Parameters:
    Input/output, int *DONE.  Before the very first call, the
    user should set DONE to TRUE, which prompts the program
    to initialize its data, and return the first partition.
    Thereafter, the user should call again, for the next
    partition, and so on, until the routine returns with DONE
    equal to TRUE, at which point there are no more partitions
    to compute.
    Input/output, int IARRAY[N], contains the information
    defining the current partition.  The user should not alter
    IARRAY between calls.  Except for the very first
    call, the routine uses the previous output value of IARRAY to compute
    the next value.
    The entries of IARRAY are the partition subset to which each
    element of the original set belongs.  If there are NPART distinct
    parts of the partition, then each entry of IARRAY will be a
    number between 1 and NPART.  Every number from 1 to NPART will
    occur somewhere in the list.  If the entries of IARRAY are
    examined in order, then each time a new partition subset occurs,
    it will be the next unused integer.
    For instance, for N = 4, the program will describe the set
    where each element is in a separate subset as 1, 2, 3, 4,
    even though such a partition might also be described as
    4, 3, 2, 1 or even 1, 5, 8, 19.
    Input, int N, the number of elements in the set.
*/
{
	const dtpipb * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	int * iarray = s_data->a1;
	bool * done = s_data->a2;
	
    dim_typ i;
    dim_typ imax;
    dim_typ j;
    dim_typ jmax;

    if ( *done )
    {
        *done = false;
        for ( i = 0; i < n; ++i )
            iarray[i] = 1;
    }
    else
    {
        /*
        Find the last element J that can be increased by 1.
        This is the element that is not equal to its maximum possible value,
        which is the maximum value of all preceding elements +1.
        */
        jmax = iarray[0];
        imax = 1;

        for ( j = 2; j <= n; ++j )
            jmax = jmax < iarray[j-1] ? iarray[j-1] : j;
        /*
        If no element can be increased by 1, we are done.
        */
        if ( imax == 1 )
        {
            *done = 1;
            return NULL;
        }
        /*
        Increase the value of the IMAX-th element by 1, set its successors to 1.
        */
        *done = false;
        ++ iarray[imax-1];
        for ( i = imax; i < n; ++i )
            iarray[i] = 1;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _equiv_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    EQUIV_RANDOM selects a random partition of a set.
  Discussion:
    The user does not control the number of parts in the partition.
    The equivalence classes are numbered in no particular order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 May 2015
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the number of elements in the set to be partitioned.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int *NPART, the number of classes or parts in the
    partition.  NPART will be between 1 and N.
    Output, int A[N], indicates the class to which each element
    is assigned.
*/
{
	const dtpipdtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	dim_typ * npart = s_data->a2;
	int * a = s_data->a3;
	
    ityp *b;
    dim_typ i;
    dim_typ j;
    dim_typ k;
    dim_typ l;
    dim_typ m;
    ityp sum1;
    dim_typ t;
    ityp z;
    ityp zhi = 1.00;
    ityp zlo = 0.00;

    b = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    b[0] = 1.00;

    for ( l = 1; l <= n-1; ++l)
    {
        sum1 = 1.00 / ( ityp ) l;
        for ( k = 1; k <= l-1; k++ )
            sum1 = ( sum1 + b[k-1] ) / ( ityp  ) ( l - k );
        b[l] = ( sum1 + b[l-1] ) / ( ityp ) ( l + 1 );
    }

    m = n;
    *npart = 0;

    for ( ; ; )
    {
        z = r8_uniform_ab ( zlo, zhi, seed );
        z = ( ityp ) ( m ) * b[m-1] * z;
        k = 0;
        ++ *npart;

        while ( 0.00 <= z )
        {
            a[m-1] = *npart;
            -- m;

            if ( m == 0 )
                break;

            z -= b[m-1];
            ++ k;
            z *= k;
        }

        if ( m == 0 )
            break;
    }
    /*
    Randomly permute the assignments.
    */
    for ( i = 0; i < n - 1; ++i )
    {
        j = i4_uniform_ab ( i, n - 1, seed );
        t    = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    /*
    Free memory.
    */
    free ( b );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _euler_row ( void * data)
/******************************************************************************/
/*
  Purpose:
    EULER_ROW returns the N-th row of Euler's triangle.
  Discussion:
    E(N,K) counts the number of permutations of the N digits that have
    exactly K "ascents", that is, K places where the Ith digit is
    less than the (I+1)th digit.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 June 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the row of Euler's triangle desired.
    Output, int IEULER[N+1], the N-th row of Euler's
    triangle, IEULER[K] contains the value of E(N,K).  Note
    that IEULER[0] should be 1 and IEULER[N] should be 0.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * ieuler = s_data->a1;
	
    dim_typ irow;
    dim_typ k;

    ieuler[0] = 1;

    if ( 0 < n )
    {
        ieuler[1] = 0;

        for ( irow = 2; irow <= n; ++irow )
        {
            ieuler[irow] = 0;

            for ( k = irow-1; 1 <= k; --k )
                ieuler[k] = ( k + 1 ) * ieuler[k] + ( irow - k ) * ieuler[k-1];
            ieuler[0] = 1;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _frobenius_number_order2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FROBENIUS_NUMBER_ORDER2 returns the Frobenius number for order 2.
  Discussion:
    The Frobenius number of order N is the solution of the Frobenius
    coin sum problem for N coin denominations.
    The Frobenius coin sum problem assumes the existence of
    N coin denominations, and asks for the largest value that cannot
    be formed by any combination of coins of these denominations.
    The coin denominations are assumed to be distinct positive integers.
    For general N, this problem is fairly difficult to handle.
    For N = 2, it is known that:
    * if C1 and C2 are not relatively prime, then
      there are infinitely large values that cannot be formed.
    * otherwise, the largest value that cannot be formed is
      C1 * C2 - C1 - C2, and that exactly half the values between
      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
    As a simple example, if C1 = 2 and C2 = 7, then the largest
    unrepresentable value is 5, and there are (5+1)/2 = 3
    unrepresentable values, namely 1, 3, and 5.
    For a general N, and a set of coin denominations C1, C2, ..., CN,
    the Frobenius number F(N, C(1:N) ) is defined as the largest value
    B for which the equation
      C1*X1 + C2*X2 + ... + CN*XN = B
    has no nonnegative integer solution X(1:N).
    In the Mathematica Package "NumberTheory", the Frobenius number
    can be determined by
    <<NumberTheory`Frobenius`
    FrobeniusF[ {C1,...,CN} ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2008
  Author:
    John Burkardt
  Reference:
    James Sylvester,
    Question 7382,
    Mathematical Questions with their Solutions,
    Educational Times,
    Volume 41, page 21, 1884.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input, int C1, C2, the coin denominations. C1 and C2
    should be positive and relatively prime.
    Output, int FROBENIUS_NUMBER_ORDER2, the Frobenius number of (C1,C2).
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int c1 = a_data[0];
	const register int c2 = a_data[1];
	
	result = c1 <= 0 || c2 <= 0 || i4_gcd ( c1, c2 ) != 1 ? i4_huge : c1 * c2 - c1 - c2;
    return  &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gray_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAY_NEXT generates the next Gray code by switching one item at a time.
  Discussion:
    On the first call only, the user must set CHANGE = -N.
    This initializes the routine to the Gray code for N zeroes.
    Each time it is called thereafter, it returns in CHANGE the index
    of the item to be switched in the Gray code.  The sign of CHANGE
    indicates whether the item is to be added or subtracted (or
    whether the corresponding bit should become 1 or 0).  When
    CHANGE is equal to N+1 on output, all the Gray codes have been
    generated.
  Example:
    N  CHANGE         Subset in/out   Binary Number
                      Interpretation  Interpretation
                       1 2 4 8
   --  ---------      --------------  --------------
    4   -4 / 0         0 0 0 0         0
        +1             1 0 0 0         1
          +2           1 1 0 0         3
        -1             0 1 0 0         2
            +3         0 1 1 0         6
        +1             1 1 1 0         7
          -2           1 0 1 0         5
        -1             0 0 1 0         4
              +4       0 0 1 1        12
        +1             1 0 1 1        13
          +2           1 1 1 1        15
        -1             0 1 1 1        14
            -3         0 1 0 1        10
        +1             1 1 0 1        11
          -2           1 0 0 1         9
        -1             0 0 0 1         8
              -4       0 0 0 0         0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 June 2015
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the order of the total set from which
    subsets will be drawn.
    Input/output, int *CHANGE.  This item is used for input only
    on the first call for a particular sequence of Gray codes,
    at which time it must be set to -N.  This corresponds to
    all items being excluded, or all bits being 0, in the Gray code.
    On output, CHANGE indicates which of the N items must be "changed",
    and the sign indicates whether the item is to be added or removed
 (or the bit is to become 1 or 0).  Note that on return from the
    first call, CHANGE is set to 0, indicating that we begin with
    the empty set.
    Input/output, int *K, a bookkeeping variable.
    The user must declare this variable before the first call.
    The output value from one call should be the input value for the next.
    The user should not change this variable.
    Input/output, int A[N], a bookkeeping variable.
    The user must declare this variable before the first call.
    The output value from one call should be the input value for the next.
    The user should not change this variable.
*/
{
	const dtpipdtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * change = s_data->a1;
	dim_typ * k = s_data->a2;
	int * a = s_data->a3;
	
    dim_typ i;

    if ( n <= 0 )
        return NULL;

    if ( *change == -n )
    {
        for ( i = 0; i < n; ++i )
            a[i] = 0;

        *k = 1;
        *change = 0;
        return NULL;
    }
    /*
    First determine WHICH item is to be changed.
    */
    if ( ( *k % 2 ) == 1 )
        *change = 1;
    else
        for ( i = 1; i <= n; ++i )
            if ( a[i-1] == 1 )
            {
                *change = i + 1;
                break;
            }
    /*
    Take care of the terminal case CHANGE = N + 1.
    */
    if ( *change == n + 1 )
        *change = n;
    /*
    Now determine HOW the item is to be changed.
    */
    if ( a[*change-1] == 0 )
        a[*change-1] = 1;
    else if ( a[*change-1] == 1 )
    {
        a[*change-1] = 0;
        *change = -( *change );
    }
    /*
    Update the counter.
    */
    ++ *k;
    /*
    If the output CHANGE = -N, then we're done.
    */
    if ( *change == -n )
        *k = 0;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gray_rank ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAY_RANK ranks a Gray code.
  Discussion:
    Given the number GRAY, its ranking is the order in which it would be
    visited in the Gray code ordering.  The Gray code ordering begins
    Rank  Gray  Gray
       (Dec) (Bin)
       0     0  0000
       1     1  0001
       2     3  0011
       3     2  0010
       4     6  0110
       5     7  0111
       6     5  0101
       7     4  0100
       8    12  0110
       etc
   This routine is given a Gray code, and has to return the rank.
  Example:
    Gray  Gray  Rank
 (Dec) (Bin)
     0       0     0
     1       1     1
     2      10     3
     3      11     2
     4     100     7
     5     101     6
     6     110     4
     7     111     5
     8    1000    15
     9    1001    14
    10    1010    12
    11    1011    13
    12    1100     8
    13    1101     9
    14    1110    11
    15    1111    10
    16   10000    31
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int GRAY, the Gray code to be ranked.
    Output, int GRAY_RANK, the rank of GRAY, and the integer whose Gray
    code is GRAY.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register int gray = *(int *) data;
	
    dim_typ i;
    register dim_typ rank;

    rank = 0;

    if ( i4_btest ( gray, 31 ) )
        rank = i4_bset ( rank, 31 );

    for ( i = 30; 0 <= i; --i )
        if ( i4_btest ( rank, i+1 ) != i4_btest ( gray, i ) )
            rank = i4_bset ( rank, i );
            
    result = rank;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gray_rank2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAY_RANK2 ranks a Gray code.
  Discussion:
    In contrast to GRAY_RANK, this routine is entirely arithmetical,
    and does not require access to bit testing and setting routines.
    Given the number GRAY, its ranking is the order in which it would be
    visited in the Gray code ordering.  The Gray code ordering begins
    Rank  Gray  Gray
       (Dec) (Bin)
       0     0  0000
       1     1  0001
       2     3  0011
       3     2  0010
       4     6  0110
       5     7  0111
       6     5  0101
       7     4  0100
       8    12  0110
       etc
   This routine is given a Gray code, and has to return the rank.
  Example:
    Gray  Gray  Rank
 (Dec) (Bin)
     0       0     0
     1       1     1
     2      10     3
     3      11     2
     4     100     7
     5     101     6
     6     110     4
     7     111     5
     8    1000    15
     9    1001    14
    10    1010    12
    11    1011    13
    12    1100     8
    13    1101     9
    14    1110    11
    15    1111    10
    16   10000    31
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int GRAY, the Gray code to be ranked.
    Output, int GRAY_RANK, the rank of GRAY, and the integer whose Gray
    code is GRAY.
*/
{
	static dim_typ result = USHRT_MAX;
	
	register int gray = *(int *) data;	
	
    dim_typ k;
    dim_typ last;
    dim_typ next;
    dim_typ rank;
    dim_typ two_k;

    if ( gray < 0 )
    {
    	result = USHRT_MAX;
        return &result;
    }

    if ( gray == 0 )
    {
    	result = 0;
        return &result;
    }
    /*
    Find TWO_K, the largest power of 2 less than or equal to GRAY.
    */
    k = 0;
    two_k = 1;
    while ( two_k << 1 <= gray )
    {
        two_k <<= 1;
        ++ k;
    }

    rank = two_k;
    last = 1;
    gray -= two_k;

    while ( 0 < k )
    {
        two_k /= 2;
        -- k;

        next = ( two_k <= gray && gray < (two_k<<1));

        if ( next )
            gray -= two_k;

        if ( next != last )
        {
            rank += two_k;
            last = 1;
        }
        else
            last = 0;
    }
    
    result = rank;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _gray_unrank ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAY_UNRANK unranks a Gray code.
  Discussion:
    The binary values of the Gray codes of successive integers differ in
    just one bit.
    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
    Hamiltonian cycle on a graph of the cube in N dimensions.
  Example:
    Rank  Gray  Gray
       (Dec) (Bin)
     0     0       0
     1     1       1
     2     3      11
     3     2      10
     4     6     110
     5     7     111
     6     5     101
     7     4     100
     8    12    1100
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int RANK, the integer whose Gray code is desired.
    Output, int GRAY_UNRANK, the Gray code of the given rank.
*/
{
	static int result = INT_MAX;
	
	const register dim_typ rank = *(dim_typ *) data;
	
    register int gray = 0;
    dim_typ i;

    if ( i4_btest ( rank, 31 ) )
        gray = i4_bset ( gray, 31 );

    for ( i = 30; 0 <= i; --i )
        if ( i4_btest ( rank, i+1 ) !=  i4_btest ( rank, i ) )
            gray = i4_bset ( gray, i );
            
	result = gray;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gray_unrank2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAY_UNRANK2 unranks a Gray code.
  Discussion:
    In contrast to GRAY_UNRANK, this routine is entirely arithmetical,
    and does not require access to bit testing and setting routines.
    The binary values of the Gray codes of successive integers differ in
    just one bit.
    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
    Hamiltonian cycle on a graph of the cube in N dimensions.
  Example:
    Rank  Gray  Gray
       (Dec) (Bin)
     0     0       0
     1     1       1
     2     3      11
     3     2      10
     4     6     110
     5     7     111
     6     5     101
     7     4     100
     8    12    1100
     9    14    1001
    10    12    1010
    11    13    1011
    12     8    1100
    13     9    1101
    14    11    1110
    15    10    1111
    16    31   10000
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int RANK, the integer whose Gray code is desired.
    Output, int GRAY_UNRANK2, the Gray code of the given rank.
*/
{
	static int result = INT_MAX;
	
	register dim_typ rank = *(dim_typ *) data;
	
    register int gray;
    dim_typ k;
    dim_typ last;
    dim_typ next;
    dim_typ two_k;

    if ( rank == 0)
    {
    	result = 0;
        return &result;
    }

    k = 0;
    two_k = 1;
    while ( (two_k<<1) <= rank )
    {
        two_k <<= 1;
        ++ k;
    }

    gray = two_k;
    rank -= two_k;
    next = 1;

    while ( 0 < k )
    {
        two_k /= 2;
        -- k;

        last = next;
        next = ( two_k <= rank && rank <= (two_k<<1) );

        if ( next != last )
            gray += two_k;

        if ( next )
            rank -= two_k;
    }
    
    result = gray;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4_partition_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_PARTITION_RANDOM selects a random partition of the int N.
  Discussion:
    Note that some elements of the partition may be 0.  The partition is
    returned as (MULT(I),I), with NPART nonzero entries in MULT, and
      N = sum ( 1 <= I <= N ) MULT(I) * I.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the integer to be partitioned.
    Input, int TABLE[N], the number of partitions of each integer
    from 1 to N.  This table may be computed by I4_PARTITION_COUNT2.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int A[N], contains in A(1:NPART) the parts of the partition.
    Output, int MULT[N], contains in MULT(1:NPART) the multiplicity
    of the parts.
    Output, int *NPART, the number of parts in the partition chosen,
    that is, the number of integers I with nonzero multiplicity MULT(I).
*/
{
	const dt4pipdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * table = s_data->a1;
	int * seed = s_data->a2;
	int * a = s_data->a3;
	int * mult = s_data->a4;
	dim_typ * npart = s_data->a5;
	
    dim_typ i;
    dim_typ i1;
    dim_typ id;
    dim_typ j;
    dim_typ m;
    ityp z;

    m = n;
    *npart = 0;
    for ( i = 0; i < n; ++i)
        mult[i] = 0;

    while ( 0 < m )
    {
        z = r8_uniform_01 ( seed );
        z *= m * table[m-1];
        id = 1;
        i1 = m;
        j = 0;

        for ( ; ; )
        {
            ++ j;
            i1 -= id;

            if ( i1 < 0 )
            {
                ++ id;
                i1 = m;
                j = 0;
                continue;
            }

            if ( i1 == 0 )
            {
                z -= id;
                if ( 0.00 < z )
                {
                    ++ id;
                    i1 = m;
                    j = 0;
                    continue;
                }
                else
                    break;
            }

            if ( 0 < i1 )
            {
                z -= id * table[i1-1];
                if ( z <= 0.0 )
                    break;
            }
        }

        mult[id-1] += j;
        *npart += j;
        m = i1;
    }
    /*
    Reformulate the partition in the standard form.
    NPART is the number of distinct parts.
    */
    *npart = 0;

    for ( i = 1; i <= n; ++i )
        if ( mult[i-1] != 0 )
        {
            ++ *npart;
            a[*npart-1] = i;
            mult[*npart-1] = mult[i-1];
        }

    for ( i = *npart+1; i <= n; ++i )
        mult[i-1] = 0;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _josephus ( void * data)
/******************************************************************************/
/*
  Purpose:
    JOSEPHUS returns the position X of the K-th man to be executed.
  Discussion:
    The classic Josephus problem concerns a circle of 41 men.
    Every third man is killed and removed from the circle.  Counting
    and executing continues until all are dead.  Where was the last
    survivor sitting?
    Note that the first person killed was sitting in the third position.
    Moreover, when we get down to 2 people, and we need to count the
    "third" one, we just do the obvious thing, which is to keep counting
    around the circle until our count is completed.
    The process may be regarded as generating a permutation of
    the integers from 1 to N.  The permutation would be the execution
    list, that is, the list of the executed men, by position number.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 March 2009
  Author:
     John Burkardt.
  Reference:
    W W Rouse Ball,
    Mathematical Recreations and Essays,
    Macmillan, 1962, pages 32-36.
    Donald Knuth,
    The Art of Computer Programming,
    Volume 1, Fundamental Algorithms,
    Addison Wesley, 1968, pages 158-159.
    Donald Knuth,
    The Art of Computer Programming,
    Volume 3, Sorting and Searching,
    Addison Wesley, 1968, pages 18-19.
  Parameters:
    Input, int N, the number of men.
    N must be positive.
    Input, int M, the counting index.
    M must not be zero.  Ordinarily, M is positive, and no greater than N.
    Input, int K, the index of the executed man of interest.
    K must be between 1 and N.
    Output, int JOSEPHUS, the position of the K-th man.
    The value will be between 1 and N.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ m = a_data[1];
	const register dim_typ k = a_data[2];
	
    dim_typ m2;
    dim_typ x;

    if ( n == 0 || m == 0 || k == 0 || n < k )
    {
    	result = USHRT_MAX;
        return &result;
    }
    /*
    In case M is bigger than N, or negative, get the
    equivalent positive value between 1 and N.
    You can skip this operation if 1 <= M <= N.
    */
    m2 = i4_modp ( m, n );
    x = k * m2;

    while ( n < x )
        x = ( m2 * ( x - n ) - 1 ) / ( m2 - 1 );

	result = x;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _matrix_product_opt ( void * data)
/******************************************************************************/
/*
  Purpose:
    MATRIX_PRODUCT_OPT determines the optimal cost of a matrix product.
  Discussion
    The cost of multiplying an LxM matrix by an M by N matrix is
    assessed as L*M*N.
    Any particular order of multiplying a set of N matrices is equivalent
    to parenthesizing an expression of N objects.
    The actual number of ways of parenthesizing an expression
    of N objects is C(N), the N-th Catalan number.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    John Burkardt
  Reference:
    Robert Sedgewick,
    Algorithms,
    Addison-Wesley, 1984, pages 486-489.
  Parameters:
    Input, int N, the number of matrices to be multiplied.
    Input, int RANK[N+1], the rank information for the matrices.
    Matrix I has RANK[I] rows and RANK[I+1] columns.
    Output, int *COST, the cost of the multiplication if the optimal
    order is used.
    Output, int ORDER[N-1], indicates the order in which the N-1
    multiplications are to be carried out.  ORDER[0] is the first
    multiplication to do, and so on.
*/
{
	const dtpipdtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * rank = s_data->a1;
	dim_typ * cost = s_data->a2;
	int * order = s_data->a3;
	
    # define STACK_MAX 100

    int *best;
    int *cost2;
    dim_typ cost3;
    dim_typ i;
    dim_typ i1;
    dim_typ i2;
    dim_typ i3;
    dim_typ j;
    dim_typ k;
    int stack[STACK_MAX];
    dim_typ stack_num;
    dim_typ step;
    /*
    Initialize the cost matrix.
    */
    best = ( int * ) malloc ( n * n * sizeof ( int ) );
    cost2 = ( int * ) malloc ( n * n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
    {
        for ( j = 0; j <= i; ++j )
            cost2[i+j*n] = 0;
        for ( j = i+1; j < n; ++j )
            cost2[i+j*n] = i4_huge;
    }
    /*
    Initialize the BEST matrix.
    */
    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j)
        best[i+j*n] = 0;
    /*
    Compute the cost and best matrices.
    */
    for ( j = 1; j <= n-1; ++j )
        for ( i = 1; i <= n-j; ++i )
            for ( k = i+1; k <= i+j; ++k )
            {
                cost3 = cost2[i-1+(k-2)*n] + cost2[k-1+(i+j-1)*n]+ rank[i-1] * rank[k-1] * rank[i+j];
                if ( cost3 < cost2[i-1+(i+j-1)*n] )
                {
                    cost2[i-1+(i+j-1)*n] = cost3;
                    best[i-1+(i+j-1)*n] = k;
                }
            }
    /*
    Pick off the optimal cost.
    */
    *cost = cost2[0+(n-1)*n];
    /*
    Backtrack to determine the optimal order.
    */
    stack_num = 0;

    i1 = 1;
    i2 = n;

    if ( i1+1 < i2 )
    {
        stack[stack_num] = i1;
        ++ stack_num;
        stack[stack_num] = i2;
        ++ stack_num;
    }

    step = n - 1;
    /*
    Take an item off the stack.
    */
    while ( 0 < stack_num )
    {
        -- stack_num;
        i3 = stack[stack_num];
        -- stack_num;
        i1 = stack[stack_num];

        i2 = best[i1-1+(i3-1)*n];

        step = step - 1;
        order[step] = i2 - 1;
        /*
        The left chunk is matrices (I1...I2-1)
        */
        if ( i1 == i2-1 );
        else if ( i1+1 == i2-1 )
        {
            -- step;
            order[step] = i2 - 2;
        }
        else
        {
            stack[stack_num] = i1;
            ++ stack_num;
            stack[stack_num] = i2 - 1;
            ++ stack_num;
        }
    /*
    The right chunk is matrices (I2...I3)
    */
    if ( i2 == i3 );
    else if ( i2+1 == i3 )
    {
        -- step;
        order[step] = i2;
    }
    else
    {
        stack[stack_num] = i2;
        ++ stack_num;
        stack[stack_num] = i3;
        ++ stack_num;
    }

    }
    free ( best );
    free ( cost2 );

    return NULL;
    # undef STACK_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _moebius_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOEBIUS_MATRIX finds the Moebius matrix from a covering relation.
  Discussion:
    This routine can be called with A and MU being the same matrix.
    The routine will correctly compute the Moebius matrix, which
    will, in this case, overwrite the input matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 June 2015
  Author:
    John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, number of elements in the partially ordered set.
    Input, int A[N*N].  A(I,J) = 1 if I is covered by J,
    0 otherwise.
    Output, int MU[N*N], the Moebius matrix as computed by the routine.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * mu = s_data->a2;
	
    dim_typ i, j;
    int *p1;
    int *p2;

    p1 = ( int * ) malloc ( n * sizeof ( n ) );
    /*
    Compute a reordering of the elements of the partially ordered matrix.
    */
    triang ( n, a, p1 );
    /*
    Copy the matrix.
    */
    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j)
            mu[i+j*n] = a[i+j*n];
    /*
    Apply the reordering to MU.
    */
    i4mat_2perm0 ( n, n, mu, p1, p1 );
    /*
    Negate the (strict) upper triangular elements of MU.
    */
    for ( i = 0; i < n-1; ++i )
        for ( j = i+1; j < n; ++j )
            mu[i+j*n] *= -1;
    /*
    Compute the inverse of MU.
    */
    i4mat_u1_inverse ( n, mu );
    /*
    All nonzero elements are reset to 1.
    */
    for ( i = 0; i < n; ++i )
        for ( j = i; j < n; ++j )
            if ( mu[i+j*n] )
                mu[i+j*n] = 1;
    /*
    Invert the matrix again.
    */
    i4mat_u1_inverse ( n, mu);
    /*
    Compute the inverse permutation.
    */
    p2 = perm0_inverse ( n, p1 );
    /*
    Unpermute the rows and columns of MU.
    */
    i4mat_2perm0 ( n, n, mu, p2, p2 );

    free ( p1 );
    free ( p2 );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _nim_sum ( void * data)
/******************************************************************************/
/*
  Purpose:
    NIM_SUM computes the Nim sum of two integers.
  Discussion:
    If K is the Nim sum of I and J, then each bit of K is the exclusive
    OR of the corresponding bits of I and J.
  Example:
     I     J     K       I_2        J_2         K_2
   ----  ----  ----  ----------  ----------  ----------
      0     0     0           0           0           0
      1     0     1           1           0           1
      1     1     0           1           1           0
      2     7     5          10         111         101
     11    28    23        1011       11100       10111
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, unsigned int UI, UJ, the integers to be Nim-summed.
    Output, unsigned int NIM_SUM, the Nim sum of I and J.
*/
{
	static unsigned result = UINT_MAX;
	
	const unsigned * const a_data = data;
	const register unsigned ui = a_data[0];
	const register unsigned uj = a_data[1];
	
    # define NBITS 32

    int bvec1[NBITS];
    int bvec2[NBITS];
    int bvec3[NBITS];


    ui4_to_ubvec ( ui, NBITS, bvec1 );
    ui4_to_ubvec ( uj, NBITS, bvec2 );
    ubvec_xor ( NBITS, bvec1, bvec2, bvec3 );
    
	result = ubvec_to_ui4 ( NBITS, bvec3 );
    return &result;
    # undef NBITS
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _padovan ( void * data)
/******************************************************************************/
/*
  Purpose:
    PADOVAN returns the first N values of the Padovan sequence.
  Discussion:
    The Padovan sequence has the initial values:
      P(0) = 1
      P(1) = 1
      P(2) = 1
    and subsequent entries are generated by the recurrence
      P(I+1) = P(I-1) + P(I-2)
  Example:
    0   1
    1   1
    2   1
    3   2
    4   2
    5   3
    6   4
    7   5
    8   7
    9   9
   10  12
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    John Burkardt
  Reference:
    Ian Stewart,
    "A Neglected Number",
    Scientific American, Volume 274, pages 102-102, June 1996.
    Ian Stewart,
    Math Hysteria,
    Oxford, 2004.
  Parameters:
    Input, int N, the number of terms.
    Output, int P[N], terms 0 though N-1 of the sequence.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
    dim_typ i;

    if ( n == 0)
        return NULL;

    p[0] = 1;

    if ( n < 2 )
        return NULL;

    p[1] = 1;

    if ( n < 3 )
        return NULL;

    p[2] = 1;

    for ( i = 4; i <= n; ++i )
        p[i-1] = p[i-3] + p[i-4];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _pell_basic ( void * data)
/******************************************************************************/
/*
  Purpose:
    PELL_BASIC returns the fundamental solution for Pell's basic equation.
  Discussion:
    Pell's equation has the form:
      X**2 - D * Y**2 = 1
    where D is a given non-square integer, and X and Y may be assumed
    to be positive integers.
  Example:
     D   X0   Y0
     2    3    2
     3    2    1
     5    9    4
     6    5    2
     7    8    3
     8    3    1
    10   19    6
    11   10    3
    12    7    2
    13  649  180
    14   15    4
    15    4    1
    17   33    8
    18   17    4
    19  170   39
    20    9    2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
   John Burkardt
  Reference:
    Mark Herkommer,
    Number Theory, A Programmer's Guide,
    McGraw Hill, 1999, pages 294-307
  Parameters:
    Input, int D, the coefficient in Pell's equation.  D should be
    positive, and not a perfect square.
    Output, int *X0, *Y0, the fundamental or 0'th solution.
    If X0 = Y0 = 0, then the calculation was canceled because of an error.
    Both X0 and Y0 will be nonnegative.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ d = s_data->a0;
	int * x0 = s_data->a1;
	int * y0 = s_data->a2;
	
    # define TERM_MAX 100

    int b[TERM_MAX+1];
    dim_typ i;
    int p;
    int pm1;
    int pm2;
    int q;
    int qm1;
    int qm2;
    int r;
    int term_num;
    /*
    Check D.
    */
    if ( d == 0 )
    {
        *x0 = *y0 = 0;
        return NULL;
    }

    i4_sqrt ( d, &q, &r );

    if ( r == 0 )
    {
        *x0 = *y0 = 0;
        return NULL;
    }
    /*
    Find the continued fraction representation of sqrt ( D ).
    */
    i4_sqrt_cf ( d, TERM_MAX, &term_num, b );
    /*
    If necessary, go for two periods.
    */
    if ( ( term_num % 2 ) == 1 )
    {
        for ( i = term_num+1; i <=  term_num << 1; ++i)
            b[i] = b[i-term_num];
        term_num <<= 1;
    }
    /*
    Evaluate the continued fraction using the forward recursion algorithm.
    */
    pm2 = qm1 = 0;
    pm1 = qm2 = 1;

    for ( i = 0; i < term_num; ++i )
    {
        p = b[i] * pm1 + pm2;
        q = b[i] * qm1 + qm2;
        pm2 = pm1;
        pm1 = p;
        qm2 = qm1;
        qm1 = q;
    }
    /*
    Get the fundamental solution.
    */
    *x0 = p;
    *y0 = q;

    return NULL;
    # undef TERM_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pell_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    PELL_NEXT returns the next solution of Pell's equation.
  Discussion:
    Pell's equation has the form:
      X**2 - D * Y**2 = 1
    where D is a given non-square integer, and X and Y may be assumed
    to be positive integers.
    To compute X0, Y0, call PELL_BASIC.
    To compute X1, Y1, call this routine, with XN and YN set to X0 and Y0.
    To compute further solutions, call again with X0, Y0 and the previous
    solution.
  Example:
    ------INPUT--------  --OUTPUT--
    D  X0  Y0   XN   YN  XNP1  YNP1
    2   3   2    3    2    17    12
    2   3   2   17   12    99    70
    2   3   2   99   70   577   408
    2   3   2  577  408  3363  2378
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
   John Burkardt
  Reference:
    Mark Herkommer,
    Number Theory, A Programmer's Guide,
    McGraw Hill, 1999, pages 294-307
  Parameters:
    Input, int D, the coefficient in Pell's equation.
    Input, int X0, Y0, the fundamental or 0'th solution.
    Input, int XN, YN, the N-th solution.
    Output, int *XNP1, *YNP1, the N+1-th solution.
*/
{
	const dt4i2pi * const s_data = data;
	const register dim_typ d = s_data->a0;
	int x0 = s_data->a1;
	int y0 = s_data->a2;
	int xn = s_data->a3;
	int yn = s_data->a4;
	int * xnp1 = s_data->a5;
	int * ynp1 = s_data->a6;
	
    *xnp1 = x0 * xn + d * y0 * yn;
    *ynp1 = x0 * yn +     y0 * xn;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perrin ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERRIN returns the first N values of the Perrin sequence.
  Discussion:
    The Perrin sequence has the initial values:
      P(0) = 3
      P(1) = 0
      P(2) = 2
    and subsequent entries are generated by the recurrence
      P(I+1) = P(I-1) + P(I-2)
    Note that if N is a prime, then N must evenly divide P(N).
  Example:
    0   3
    1   0
    2   2
    3   3
    4   2
    5   5
    6   5
    7   7
    8  10
    9  12
   10  17
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    John Burkardt
  Reference:
    Ian Stewart,
    "A Neglected Number",
    Scientific American, Volume 274, pages 102-102, June 1996.
    Ian Stewart,
    Math Hysteria,
    Oxford, 2004.
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 1999.
  Parameters:
    Input, integer N, the number of terms.
    Output, integer P(N), the terms 0 through N-1 of the sequence.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
    dim_typ i;

    if ( n == 0 )
        return NULL;

    p[0] = 3;

    if ( n < 2 )
        return NULL;

    p[1] = 0;

    if ( n < 3 )
        return NULL;

    p[2] = 2;

    for ( i = 4; i <= n; ++i )
        p[i-1] = p[i-3] + p[i-4];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pord_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    PORD_CHECK checks a matrix representing a partial ordering.
  Discussion:
    The array A is supposed to represent a partial ordering of
    the elements of a set of N objects.
    For distinct indices I and J, the value of A(I,J) is:
      1, if I << J
      0, otherwise ( I and J may be unrelated, or perhaps J << I).
    Diagonal elements of A are ignored.
    This routine checks that the values of A do represent
    a partial ordering.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements in the set.
    Input, int A[N*N], the partial ordering.  A[I+J*N] is
    1 if I is less than J in the partial ordering,
    0 otherwise.
    Output, int PORD_CHECK, is 1 if an error was detected.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;

    if ( n <= 0 )
    {
    	result = true;
        return &result;
    }

    for ( i = 0; i < n; ++i )
        for ( j = i+1; j < n; ++j )
            if (  0 < a[i+j*n] && 0 < a[j+i*n] )
            {
            	result = true;
        		return &result;
            }

    result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _power_series1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    POWER_SERIES1 computes a power series for a function G(Z) = (1+F(Z))**ALPHA.
  Discussion:
    The power series for F(Z) is given.
    The form of the power series are:
      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the number of terms in the power series.
    Input, double ALPHA, the exponent of 1+F(Z) in the definition of G(Z).
    Input, double A[N], the power series coefficients for F(Z).
    Output, double B[N], the power series coefficients for G(Z).
*/
{
	const dtit2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp alpha = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	
    dim_typ i, j;
    ityp v;

    for ( j = 1; j <= n; ++j )
    {
        v = 0.00;
        for ( i = 1; i <= j-1; ++i )
            v += b[i-1] * a[j-i-1] * ( alpha * ( j - i ) - i );

        b[j-1] = alpha * a[j-1] + v / ( ( ityp ) j );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _power_series2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    POWER_SERIES2 computes the power series for a function G(Z) = EXP(F(Z)) - 1.
  Discussion:
    The power series for F(Z) is given.
    The power series have the form:
      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the number of terms in the power series.
    Input, double A[N], the power series coefficients for F(Z).
    Output, double B[N], the power series coefficients for G(Z).
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    dim_typ i, j;
    ityp v;

    for ( j = 1; j <= n; ++j )
    {
        v = 0.00;

        for ( i = 1; i <= j-1; ++i )
            v += b[i-1] * a[j-i-1] * ( ityp ) ( j - i );

        b[j-1] = a[j-1] + v / ( ityp ) j;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _power_series3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    POWER_SERIES3 computes the power series for a function H(Z) = G(F(Z)).
  Discussion:
    The power series for F and G are given.
    We assume that
      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the number of terms in the power series.
    Input, double A[N], the power series for F
    Input, double B[N], the power series for G.
    Output, double C[N], the power series for H.
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * c = s_data->a3;
	
    dim_typ i;
    dim_typ iq;
    dim_typ j;
    dim_typ m;
    ityp r;
    ityp v;
    ityp *work = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        work[i] = b[0] * a[i];
    /*
    Search for IQ, the index of the first nonzero entry in A.
    */
    iq = 0;

    for ( i = 1; i <= n; ++i )
        if ( a[i-1] != 0.00 )
        {
            iq = i;
            break;
        }

    if ( iq != 0 )
    {
        m = 1;

        for ( ; ; )
        {
            ++ m;

            if ( n < m * iq )
                break;

            if ( b[m-1] == 0.00 )
                continue;

            r = b[m-1] * pow ( a[iq-1], m );
            work[m*iq-1] = work[m*iq-1] + r;

            for ( j = 1; j <= n-m*iq; ++j )
            {
                v = 0.00;
                for ( i = 1; i <= j-1; ++i )
                v += c[i-1] * a[j-i+iq-1] * ( ityp ) ( m * ( j - i ) - i );

                c[j-1] = ( ( ityp ) m * a[j-1] + v / ( ityp ) j ) / a[iq-1];

            }

            for ( i = 1; i <= n-m*iq; ++i )
                work[i+m*iq-1] += c[i-1] * r;
        }
    }

    for ( i = 0; i < n; ++i )
        c[i] = work[i];

    free ( work );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _power_series4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    POWER_SERIES4 computes the power series for a function H(Z) = G ( 1/F(Z) ).
  Discussion:
    POWER_SERIES4 is given the power series for the functions F and G.
    We assume that
      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the number of terms in the power series.
    Input, double A[N], the power series for F.  For this problem, A(1)
    may not be 0.0.
    Input, double B(N), the power series for G.
    Output, double C(N), the power series for H.
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * c = s_data->a3;
	
    dim_typ i;
    dim_typ l;
    dim_typ m;
    ityp s;
    ityp t;
    ityp *work;

    if ( a[0] == 0.00 )
        return NULL;

    work = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    t = 1.00;

    for ( i = 0; i < n; ++i )
    {
        t /= a[0];
        c[i] = b[i] * t;
        work[i] = a[i] * t;
    }

    for ( m = 2; m <= n; ++m )
    {
        s = -work[m-1];
        for ( i = m; i <= n; ++i )
            for ( l = i; l <= n; ++l )
            {
                    c[l-1] = c[l-1] + s * c[l-m];
                    work[l-1] += s * work[l-m];
            }
    }

    free ( work );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_permanent ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_PERMANENT computes the permanent of an R8MAT.
  Discussion:
    The permanent function is similar to the determinant.  Recall that
    the determinant of a matrix may be defined as the sum of all the
    products:
      S * A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
    where I is any permutation of the columns of the matrix, and S is the
    sign of the permutation.  By contrast, the permanent function is
    the (unsigned) sum of all the products
      A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
    where I is any permutation of the columns of the matrix.  The only
    difference is that there is no permutation sign multiplying each summand.
    Symbolically, then, the determinant of a 2 by 2 matrix
      a b
      c d
    is a*d-b*c, whereas the permanent of the same matrix is a*d+b*c.
    The permanent is invariant under row and column permutations.
    If a row or column of the matrix is multiplied by S, then the
      permanent is likewise multiplied by S.
    If the matrix is square, then the permanent is unaffected by
      transposing the matrix.
    Unlike the determinant, however, the permanent does change if
      one row is added to another, and it is not 1 that the
      permanent of the product is the product of the permanents.
    Note that if A is a matrix of all 1's and 0's, then the permanent
    of A counts exactly which permutations hit exactly 1's in the matrix.
    This fact can be exploited for various combinatorial purposes.
    For instance, setting the diagonal of A to 0, and the offdiagonals
    to 1, the permanent of A counts the number of derangements of N
    objects.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, number of rows and columns in matrix.
    Input, double A[N*N], the matrix whose permanent is desired.
    Output, double R8MAT_PERMANENT, the value of the permanent of A.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data; 
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    int iadd;
    int *iwork;
    dim_typ j;
    int more;
    int ncard;
    ityp p;
    ityp perm;
    ityp prod;
    ityp sgn;
    ityp *work;
    ityp z;

    more = false;

    iwork = ( int * ) malloc ( n * sizeof ( int ) );
    work = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 1; i <= n; ++i )
    {
        work[i-1] = a[i-1+(n-1)*n];
        for ( j = 1; j <= n; ++j )
            work[i-1] -= 0.50 * a[i-1+(j-1)*n];
    }

    p = 0.00;
    sgn = -1.00;

    for ( ; ; )
    {
        sgn *= -1;

        subset_gray_next ( n-1, iwork, &more, &ncard, &iadd );

        if ( ncard != 0 )
        {
            z = ( ityp ) ( (iwork[iadd-1]<<1) - 1 );
            for ( i = 1; i <= n; ++i )
                work[i-1] += z * a[i-1+(iadd-1)*n];
        }

        prod = 1.00;
        for ( i = 0; i < n; ++i )
            prod *= work[i];
        p += sgn * prod;

        if ( !more )
            break;
    }

    free ( iwork );
    free ( work );

    perm = p * ( ityp ) ( ( ( n % 2 ) << 2) - 2 );
	
	result = perm;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_backtrack ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_BACKTRACK supervises a backtrack search for a vector.
  Discussion:
    The routine tries to construct a vector one index at a time,
    using possible candidates as supplied by the user.
    At any time, the partially constructed vector may be discovered to be
    unsatisfactory, but the routine records information about where the
    last arbitrary choice was made, so that the search can be
    carried out efficiently, rather than starting out all over again.
    First, call the routine with INDX = 0 so it can initialize itself.
    Now, on each return from the routine, if INDX is:
      1, you've just been handed a complete candidate vector;
         Admire it, analyze it, do what you like.
      2, please determine suitable candidates for position X(K).
         Return the number of candidates in NCAN(K), adding each
         candidate to the end of STACK, and increasing NSTACK.
      3, you're done.  Stop calling the routine;
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the number of positions to be filled in the vector.
    Input, int MAXSTACK, the maximum length of the stack.
    Input/output, double X[N], the partial or complete candidate vector.
    Input/output, int *INDX, a communication flag.
    On input,
      0 to start a search.
    On output:
      1, a complete output vector has been determined and returned in X(1:N);
      2, candidates are needed for position X(K);
      3, no more possible vectors exist.
    Input/output, int *K, if INDX=2, the current vector index being considered.
    Input/output, int *NSTACK, the current length of the stack.
    Input/output, double STACK[MAXSTACK], a list of all current candidates for
    all positions 1 through K.
    Input/output, int NCAN[N], lists the current number of candidates for
    positions 1 through K.
*/
{
	const _2dtpitpi2pdtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ maxstack = s_data->a1;
	ityp * x = s_data->a2;
	int * indx = s_data->a3;
	dim_typ * k = s_data->a4;
	dim_typ * nstack = s_data->a5; 
	ityp * stack = s_data->a6;
	int * ncan = s_data->a7;
	
    /*
    If this is the first call, request a candidate for position 1.
    */
    if ( *indx == 0 )
    {
        *k = 1;
        *nstack = 0;
        *indx = 2;
        return NULL;
    }
    /*
    Examine the stack.
    */
    for ( ; ; )
    {
        /*
        If there are candidates for position K, take the first available
        one off the stack, and increment K.

        This may cause K to reach the desired value of N, in which case
        we need to signal the user that a complete set of candidates
        is being returned.
        */
        if ( 0 < ncan[(*k)-1] )
        {
            x[(*k)-1] = stack[(*nstack)-1];
            -- *nstack;

            -- ncan[(*k)-1];

            if ( *k != n )
            {
                ++ *k ;
                *indx = 2;
            }
            else
                *indx = 1;

                break;
        }
        /*
        If there are no candidates for position K, then decrement K.
        If K is still positive, repeat the examination of the stack.
        */
        else
        {
            -- *k;

            if ( *k <= 0 )
            {
                *indx = 3;
                break;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _rat_farey ( void * data)
/******************************************************************************/
/*
  Purpose:
    RAT_FAREY computes the N-th row of the Farey fraction table.
  Example:
    N = 5
    NUM_FRAC = 11
    A =  0  1  1  1  2  1  3  2  3  4  1
    B =  1  5  4  3  5  2  5  3  4  5  1
  Discussion:
    In this form of the Farey fraction table, fractions in row N lie between
    0 and 1, are in lowest terms, and have a denominator that is no greater
    than N.  Row N is computed directly, and does not require the computation
    of previous rows.
    The data satisfy the relationship:
      A(K+1) * B(K) - A(K) * B(K+1) = 1
    The number of items in the N-th row is roughly N**2 / M_PI**2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 March 2009
  Author:
    John Burkardt
  Reference:
    Donald Knuth,
    The Art of Computer Programming,
    Volume 1, Fundamental Algorithms,
    Addison Wesley, 1968, page 157.
  Parameters:
    Input, int N, the desired row number.  N must be positive.
    Input, int MAX_FRAC, the maximum number of fractions to compute.
    Output, int NUM_FRAC, the number of fractions computed.
    Output, int A[MAX_FRAC], B[MAX_FRAC], contains the NUM_FRAC
    numerators and denominators of the N-th row of the Farey fraction table.
*/
{
	const _2dtpdt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ max_frac = s_data->a1;
	dim_typ * num_frac = s_data->a2;
	int * a = s_data->a3;
	int * b = s_data->a4;
	
    dim_typ c;
    dim_typ k;

    if ( n == 0 )
    {
        *num_frac = 0;
        return NULL;
    }

    k = 0;

    if ( max_frac == 0 )
    {
        *num_frac = k;
        return NULL;
    }

    a[k] = 0;
    b[k] = k = 1;

    if ( max_frac <= 1 )
    {
        *num_frac = k;
        return NULL;
    }

    a[k] = 1;
    b[k] = n;
    k = 2;

    while ( k < max_frac )
    {
        if ( a[k-1] == 1 && b[k-1] == 1 )
            break;

        c = ( b[k-2] + n ) / b[k-1];
        a[k] = c * a[k-1] - a[k-2];
        b[k] = c * b[k-1] - b[k-2];
        ++ k;
    }

    *num_frac = k;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rat_farey2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    RAT_FAREY2 computes the next row of the Farey fraction table.
  Example:
    Input:
      N = 3
      A =  0  1  1  2  1
      B =  1  3  2  3  1
    Output:
      A =  0  1  1  2  1  3  2  3  1
      B =  1  4  3  5  2  5  3  4  1
  Discussion:
    In this form of the Farey fraction table, fractions in row N lie between
    0 and 1, and are in lowest terms.  For every adjacent pair of input
    fractions, A1/B1 and A2/B2, the mediant (A1+A2)/(B1+B2) is computed
    and inserted between them.
    The number of items in the N-th row is 1+2^(N-1).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the input row number.  N must be nonnegative.
    If N is zero, then the input is ignored, and the entries of
    row 1 are computed directly.
    Input/output, int A[1+2^N], B[1+2^N].
    On input, entries 1 through 1+2^(N-1) contain the entries of row N.
    On output, entries 1 through 1+2^N contain the entries of row N+1.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * b = s_data->a2;
	
    dim_typ i;
    dim_typ ihi;

    if ( n == 0 )
    {
        a[0] = 0;
        b[0] = a[1] = b[1];
        return NULL;
    }
    /*
    Shift the current data.
    */
    ihi = powi ( 2, n - 1 );
    for ( i = ihi; 0 <= i; --i)
    {
        a[i<<1] = a[i];
        b[i<<1] = b[i];
    }
    /*
    Compute the mediants.
    */
    ihi = powi ( 2, n ) - 1;

    for ( i = 1; i <= ihi; i +=2 )
    {
        a[i] = a[i-1] + a[i+1];
        b[i] = b[i-1] + b[i+1];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _schroeder ( void * data)
/******************************************************************************/
/*
  Purpose:
    SCHROEDER generates the Schroeder numbers.
  Discussion:
    The Schroeder number S(N) counts the number of ways to insert
    parentheses into an expression of N items, with two or more items within
    a parenthesis.
    Note that the Catalan number C(N) counts the number of ways
    to legally arrange a set of N left and N right parentheses.
    The formula is:
      S(N) = ( P(N)(3.0) - 3 P(N-1)(3.0) ) / ( 4 * ( N - 1 ) )
    where P(N)(X) is the N-th Legendre polynomial.
  Example:
    N = 4
    1234
    12(34)
    1(234)
    1(2(34))
    1(23)4
    1((23)4)
 (123)4
 (12)34
 (12)(34)
 (1(23))4
 ((12)3)4
  First Values:
           1
           1
           3
          11
          45
         197
         903
        4279
       20793
      103049
      518859
     2646723
    13648869
    71039373
  Recursion:
    S(1) = 1
    S(2) = 1
    S(N) = ( ( 6 * N - 9 ) * S(N-1) - ( N - 3 ) * S(N-2) ) / N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    John Burkardt
  Reference:
    R P Stanley,
    Hipparchus, Plutarch, Schroeder, and Hough,
    American Mathematical Monthly,
    Volume 104, Number 4, 1997, pages 344-350.
    Laurent Habsieger, Maxim Kazarian, Sergei Lando,
    On the Second Number of Plutarch,
    American Mathematical Monthly, May 1998, page 446.
  Parameters:
    Input, int N, the number of Schroeder numbers desired.

    Outpt, int S[N], the Schroeder numbers.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * s = s_data->a1;
	
    dim_typ i;

    if ( n == 0)
        return NULL;

    s[0] = 1;

    if ( n <= 1)
        return NULL;

    s[1] = 1;

    if ( n <= 2 )
        return NULL;

    for ( i = 3; i <= n; ++i )
        s[i-1] = ( ( 6 * i - 9 ) * s[i-2] - ( i - 3 ) * s[i-3] ) / i;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subset_gray_rank ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_GRAY_RANK ranks a subset of an N set, using the Gray code ordering.
  Example:
    N = 4
       A       Rank
    -------   -----
    0 0 0 0       1
    0 0 0 1       2
    0 0 1 1       3
    0 0 1 0       4
    0 1 1 0       5
    0 1 1 1       6
    0 1 0 1       7
    0 1 0 0       8
    1 1 0 0       9
    1 1 0 1      10
    1 1 1 1      11
    1 1 1 0      12
    1 0 1 0      13
    1 0 1 1      14
    1 0 0 1      15
    1 0 0 0      16
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the total set from which
    subsets will be drawn.
    Input, int A[N]; A(I) is 1 if element I is in the set,
    and 0 otherwise.
    Output, int SUBSET_GRAY_RANK, the rank of the subset in the
    Gray code ordering.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
	result = gray_rank ( ( int ) ubvec_to_ui4 ( n, a ) +1);
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subset_gray_unrank ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_GRAY_UNRANK produces a subset of an N set of the given Gray code rank.
  Example:
    N = 4
     Rank     A
    -----  -------
        1  0 0 0 0
        2  0 0 0 1
        3  0 0 1 1
        4  0 0 1 0
        5  0 1 1 0
        6  0 1 1 1
        7  0 1 0 1
        8  0 1 0 0
        9  1 1 0 0
       10  1 1 0 1
       11  1 1 1 1
       12  1 1 1 0
       13  1 0 1 0
       14  1 0 1 1
       15  1 0 0 1
       16  1 0 0 0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 December 2014
  Author
    John Burkardt
  Parameters:
    Input, int RANK, the rank of the subset in the Gray code ordering.
    Input, int N, the order of the total set from which
    subsets will be drawn.
    Output, int A[N]; A(I) is 1 if element I is in the set,
    and 0 otherwise.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ rank = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    ui4_to_ubvec ( ( unsigned ) gray_unrank ( rank - 1 ), n, a );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subset_lex_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_LEX_NEXT generates the subsets of a set of N elements, one at a time.
  Discussion:
    The subsets are generated in lexicographical order.
    The routine can also be forced to generate only those subsets whose
    size is no greater than some user-specified maximum.
  Example:
    N = 5, JMP = ( K == 3 )
    1
    1 2
    1 2 3
    1 2 4
    1 2 5
    1 3
    1 3 4
    1 3 5
    1 4
    1 4 5
    1 5
    2
    2 3
    2 3 4
    2 3 5
    2 4
    2 4 5
    2 5
    3
    3 4
    3 4 5
    3 5
    4
    4 5
    5
    empty set.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the order of the main set from which subsets
    are chosen.
    Input, int JMP.  In the simplest case, set JMP = .FALSE. for
    a normal computation.  But to jump over supersets of the input set,
    set JMP = TRUE.  Setting JMP = ( K == 3 ) before every new call
    will, for example, force all the subsets returned
    to have cardinality 3 or less.
    Input, int NDIM, the allowed storage for A.  If NDIM < N,
    JMP must be used to avoid creation of a subset too large to store in A.
    Input/output, int *K.  On first call, the user must set K = 0 as
    a startup signal to the program.  Thereafter, the routine returns
    the size of the computed subset in K.  On the last return,
    the empty set is returned and K is 0, which is a signal to
    the user that the computation is complete.
    Input/output, int A[NDIM].  A(I) is the I-th element of the
    subset, listed in increasing order, with 0's in entries
    beyond entry K.
*/
{
	const _3dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ jmp = s_data->a1;
	const register dim_typ ndim = s_data->a2;
	int * k = s_data->a3;
	int * a = s_data->a4;

    dim_typ is;

    if ( *k <= 0 )
    {
        if ( jmp )
            return NULL;
        is = 0;
        *k = a[0] = 1;
    }
    else if ( a[*k-1] != n )
    {
        is = a[*k-1];

        if ( !jmp )
            ++ *k;

        a[*k-1] = is + 1;
    }
    else
    {
        -- *k;

        if ( *k != 0 )
            ++ a[*k-1];
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subset_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_RANDOM selects a random subset of an N-set.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 May 2003
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the size of the full set.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int A[N].  A vector to hold the information about
    the set chosen.  On return, if A[I] = 1, then
    I is in the random subset, otherwise, A[I] = 0
    and I is not in the random subset.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	int * a = s_data->a2;
	
    for (dim_typ i = 0; i < n; ++i )
        a[i] = i4_uniform_ab ( 0, 1, seed );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _thue_binary_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    THUE_BINARY_NEXT returns the next element in a binary Thue sequence.
  Discussion:
    Thue demonstrated that arbitrarily long sequences of 0's and
    1's could be generated which had the "cubefree" property.  In
    other words, for a given string S, there was no substring W
    such that S contained "WWW".  In fact, a stronger result holds:
    if "a" is the first letter of W, it is never the case that S
    contains the substring "WWa".
    In this example, the digits allowed are binary, that is, just
    "0" and "1".  The replacement rules are:
    "0" -> "01"
    "1" -> "10"
    This routine produces the next binary Thue sequence in a given series.
    However, the input sequence must be a Thue sequence in order for
    us to guarantee that the output sequence will also have the
    cubic nonrepetition property.
    Also, enough space must be set aside in THUE to hold the
    output sequence.  This will always be twice the input
    value of N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N.  On input, the length of the input sequence.
    On output, the length of the output sequence.
    Input, int THUE[*N].  On input, the initial Thue sequence, and on
    output, the result of applying the substitution rules once.
*/
{
	const pdtpi * const s_data = data;
	dim_typ * n = s_data->a0;
	int * thue = s_data->a1;
	
    dim_typ i;
    dim_typ n_out;
    int *thue_out;

    n_out = 0;
    thue_out = ( int * ) malloc ( (*n) * sizeof ( int ) << 1);

    for ( i = 0; i < *n; i++ )
    {
        if ( thue[i] == 0 )
        {
            thue_out[n_out] = 0;
            n_out = n_out + 1;
            thue_out[n_out] = 1;
            n_out = n_out + 1;
        }
        else if ( thue[i] == 1 )
        {
            thue_out[n_out] = 1;
            ++ n_out;
            thue_out[n_out] = 0;
            ++ n_out;
        }
        else
            return NULL;
    }

    *n = n_out;

    for ( i = 0; i < *n; ++i)
        thue[i] = thue_out[i];

    free ( thue_out );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _thue_ternary_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    THUE_TERNARY_NEXT returns the next element in a ternary Thue sequence.
  Discussion:
    Thue was interested in showing that there were arbitrarily long
    sequences of digits which never displayed a pair of contiguous
    repetitions of any length.  That is, there was no occurrence of
    "00" or "1010" or "121121", anywhere in the string.  This makes
    the string "squarefree".
    To do this, he demonstrated a way to start with a single digit,
    and to repeatedly apply a series of transformation rules to each
    digit of the sequence, deriving nonrepeating sequences of ever
    greater length.
    In this example, the digits allowed are ternary, that is, just
    "0", "1" and "2".  The replacement rules are:
    "0" -> "12"
    "1" -> "102"
    "2" -> "0"
    This routine produces the next Thue sequence in a given series.
    However, the input sequence must be a Thue sequence in order for
    us to guarantee that the output sequence will also have the
    nonrepetition property.
    Also, enough space must be set aside in THUE to hold the
    output sequence.  This will never be more than 3 times the input
    value of N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2009
  Author:
    John Burkardt
  Reference:
    Brian Hayes,
    Third Base,
    American Scientist,
    Volume 89, Number 6, pages 490-494, November-December 2001.
  Parameters:
    Input/output, int *N.  On input, the length of the input sequence.
    On output, the length of the output sequence.
    Input, int THUE[*N].  On input, the initial Thue sequence, and on
    output, the result of applying the substitution rules once.
*/
{
	const pdtpi * const s_data = data;
	dim_typ * n = s_data->a0;
	int * thue = s_data->a1;
	
    dim_typ i;
    dim_typ n_out;
    int *thue_out;

    n_out = 0;
    thue_out = ( int * ) malloc ( 3 * (*n) * sizeof ( int ) );

    for ( i = 0; i < *n; ++i )
    {

        if ( thue[i] == 0 )
        {
            thue_out[n_out] = 1;
            ++ n_out;
            thue_out[n_out] = 2;
            ++ n_out;
        }
        else if ( thue[i] == 1 )
        {
            thue_out[n_out] = 1;
            ++ n_out;
            thue_out[n_out] = 0;
            ++ n_out;
            thue_out[n_out] = 2;
           ++ n_out;
        }
        else if ( thue[i] == 2 )
        {
            thue_out[n_out] = 0;
            ++ n_out;
        }
        else
            return NULL;
    }

    *n = n_out;
    for ( i = 0; i < n_out; ++i )
        thue[i] = thue_out[i];

    free ( thue_out );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _vector_constrained_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    VECTOR_CONSTRAINED_NEXT returns the "next" constrained vector.
  Discussion:
    We consider all vectors of dimension N whose components
    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
    We are only interested in the subset of these vectors which
    satisfy the following constraint:
      sum ( 1 <= I <= N ) ( ( X(I) - 1 ) / X_MAX(I) ) <= 1
    We can carry out this check using integer arithmetic if we
    multiply through by P = product ( X_MAX(1:N) ):
      sum ( 1 <= I <= N ) ( ( X(I) - 1 ) * ( P / X_MAX(I) ) ) <= P.
    This routine returns, one at a time, and in right-to-left
    lexicographic order, exactly those vectors which satisfy
    the constraint.
  Example:
    N = 3
    X_MIN:   2   2   1
    X_MAX:   4   5   3
    P = 60
    #  X(1)  X(2)  X(3)  CONSTRAINT
    1    2     2     1       27
    2    3     2     1       42
    3    4     2     1       57
    4    2     3     1       39
    5    3     3     1       54
    6    2     4     1       51
    7    2     2     2       47
    8    2     3     2       59
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of components in the vector.
    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
    values allowed in each component.
    Input/output, int X[N].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Output, int *CONSTRAINT, the constraint value for X.  Valid vectors X
    will have a CONSTRAINT value between product(X_MIN(1:N)) (automatically)
    and product(X_MAX(1:N)) (because we skip over vectors with a
    constraint larger than this value).
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dt4pipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x_min = s_data->a1;
	int * x_max = s_data->a2;
	int * x = s_data->a3;
	int * constraint = s_data->a4;
	bool * more = s_data->a5;
	
    dim_typ i;
    dim_typ  j;
    int x_prod = 1;

    for ( j = 0; j < n; ++j)
        x_prod *= x_max[j];

    if ( !( *more ) )
    {
        for ( j = 0; j < n; ++j )
            x[j] = x_min[j];

        *constraint = 0;
        for ( j = 0; j < n; ++j )
            *constraint +=( ( x[j] - 1 ) * ( x_prod / x_max[j] ) );

        *more = x_prod >= *constraint;
        return NULL;
    }
    else
    {
        i = 0;

        for ( ; ; )
        {
            if ( x[i] < x_max[i] )
            {
                ++ x[i];

                *constraint = 0;
                for ( j = 0; j < n; ++j )
                    *constraint += ( ( x[j] - 1 ) * ( x_prod / x_max[j] ) );

                if ( *constraint <= x_prod )
                    break;
            }

            x[i] = x_min[i];

            ++ i;

            if ( n <= i )
            {
                *more = false;
                break;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _vector_constrained_next2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VECTOR_CONSTRAINED_NEXT2 returns the "next" constrained vector.
  Discussion:
    We consider all vectors of dimension N whose components
    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
    We are only interested in the subset of these vectors which
    satisfy the following constraint:
      sum ( 1 <= I <= N ) ( X(I) / X_MAX(I) ) <= 1
    We can carry out this check using integer arithmetic if we
    multiply through by P = product ( X_MAX(1:N) ):
      sum ( 1 <= I <= N ) ( X(I)  * ( P / X_MAX(I) ) ) <= P.
    This routine returns, one at a time, and in right-to-left
    lexicographic order, exactly those vectors which satisfy
    the constraint.
  Example:
    N = 3
    X_MIN:   1   1   1
    X_MAX:   5   6   4
    P = 120
    #  X(1)  X(2)  X(3)  CONSTRAINT
    1    1     1     1       74
    2    2     1     1       98
    3    1     2     1       94
    4    2     2     1      119
    5    1     3     1      114
    6    1     1     2      104
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of components in the vector.
    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
    values allowed in each component.
    Input/output, int X[N].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Output, int *CONSTRAINT, the constraint value for X.  Valid vectors X
    will have a CONSTRAINT value between product(X_MIN(1:N)) (automatically)
    and product(X_MAX(1:N)) (because we skip over vectors with a
    constraint larger than this value).
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dt4pipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x_min = s_data->a1;
	int * x_max = s_data->a2;
	int * x = s_data->a3;
	int * constraint = s_data->a4;
	bool * more = s_data->a5;
	
    dim_typ i;
    dim_typ j;
    int x_prod = 1;

    for ( j = 0; j < n; ++j)
        x_prod *= x_max[j];

    if ( !(*more) )
    {
        for ( j = 0; j < n; ++j )
            x[j] = x_min[j];

        *constraint = 0;
        for ( j = 0; j < n; ++j )
            *constraint += ( x[j] * ( x_prod / x_max[j] ) );

        *more = x_prod >= *constraint;

        return NULL;
    }
    else
    {
        i = 0;

        for ( ; ; )
        {
            if ( x[i] < x_max[i] )
            {
                ++ x[i];

                *constraint = 0;
                for ( j = 0; j < n; ++j )
                    *constraint += ( x[j] * ( x_prod / x_max[j] ) );

                if ( *constraint <= x_prod )
                    break;
            }

            x[i] = x_min[i];

            ++ i ;

            if ( n <= i )
            {
                *more = false;
                break;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _vector_constrained_next3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VECTOR_CONSTRAINED_NEXT3 returns the "next" constrained vector.
  Discussion:
    This routine addresses the same problem as VECTOR_CONSTRAINED_NEXT2,
    and differs only in that real arithmetic is used, rather than
    integer arithmetic.  Integer arithmetic allows us to do an exact
    calculation, but we run into overflow problems in simple cases
    where N is 10 and the X_MAX entries are of order 10, for instance.
    We consider all vectors of dimension N whose components
    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
    We are only interested in the subset of these vectors which
    satisfy the following constraint:
      sum ( 1 <= I <= N ) ( X(I) / X_MAX(I) ) <= 1
  Example:
    N = 3
    X_MIN:   1   1   1
    X_MAX:   5   6   4
    P = 120
    #  X(1)  X(2)  X(3)  CONSTRAINT
    1    1     1     1       0.62
    2    2     1     1       0.82
    3    1     2     1       0.78
    4    2     2     1       0.98
    5    1     3     1       0.95
    6    1     1     2       0.87
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of components in the vector.
    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
    values allowed in each component.
    Input/output, int X[N].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Output, double *CONSTRAINT, the constraint value for X.  Valid vectors
    X will have a CONSTRAINT value between
      product(X_MIN(1:N)) / product(X_MAX(1:N))
    and 1.0.
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dt3pipitpb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x_min = s_data->a1;
	int * x_max = s_data->a2;
	int * x = s_data->a3;
	ityp * constraint = s_data->a4;
	bool * more = s_data->a5;
	
    dim_typ i, j;

    if ( !( *more ) )
    {
        for ( j = 0; j < n; ++j)
            x[j] = x_min[j];

        *constraint = 0.00;
        for ( j = 0; j < n; ++j )
            *constraint += ( ityp ) x[j] / ( ityp ) x_max[j];

        *more = 1.00 >= *constraint;

        return NULL;
    }
    else
    {
        i = 0;

        for ( ; ; )
        {
            if ( x[i] < x_max[i] )
            {
                ++ x[i];

                *constraint = 0;
                for ( j = 0; j < n; ++j )
                    *constraint += ( ityp ) x[j] / ( ityp ) x_max[j];

                if ( *constraint <= 1.0 )
                    break;
            }

            x[i] = x_min[i];

            ++ i;

            if ( n <= i )
            {
                *more = false;
                break;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _vector_constrained_next4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VECTOR_CONSTRAINED_NEXT4 returns the "next" constrained vector.
  Discussion:
    This routine is similar to VECTOR_CONSTRAINED_NEXT2 and
    VECTOR_CONSTRAINED_NEXT3.
    We consider all vectors X of dimension N whose components
    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
    We are only interested in the subset of these vectors which
    satisfy the following constraint:
      sum ( 1 <= I <= N ) ( ALPHA(I) * X(I) ) <= Q
  Example:
    N = 3
    ALPHA    4.0  3.0  5.0
    Q       20.0
    X_MIN:   1   1   1
    X_MAX:   5   6   4
    P = 120
    #  X(1)  X(2)  X(3)      Total
    1    1     1     1       12.0
    2    2     1     1       20.0
    3    1     2     1       15.0
    4    2     2     1       19.0
    5    1     3     1       18.0
    6    1     1     2       17.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of components in the vector.
    Input, double ALPHA[N], the coefficient vector.
    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
    values allowed in each component.
    Input/output, int X[N].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Input, double Q, the limit on the sum.
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dtpit2pipititpb * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	int * x_min = s_data->a2;
	int * x_max = s_data->a3;
	ityp * x = s_data->a4;
	ityp q = s_data->a5;
	bool * more = s_data->a6;
	
    dim_typ i, j;
    ityp total;

    if ( ! (*more) )
    {
        for ( j = 0; j < n; ++j )
            x[j] = x_min[j];

        total = 0.00;
        for ( j = 0; j < n; ++j )
            total += alpha[j] * ( ityp ) x[j];

        *more = q >= total;
        return NULL;
    }
    else
    {
        i = 0;

        for ( ; ; )
        {
            if ( x[i] < x_max[i] )
            {
                ++ x[i];

                total = 0;
                for ( j = 0; j < n; ++j )
                    total += alpha[j] * ( ityp ) x[j];

                if ( total <= q )
                    break;
            }

            x[i] = x_min[i];

            ++ i;

            if ( n <= i )
            {
                *more = 0;
                break;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _vector_constrained_next5 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VECTOR_CONSTRAINED_NEXT5 returns the "next" constrained vector.
  Discussion:
    We consider all positive integer vectors of dimension N whose
    components satisfy SUM_MIN <= X(1:N) <= SUM_MAX.
    This routine returns, one at a time, and in right-to-left
    lexicographic order, exactly those vectors which satisfy
    the constraint.
  Example:
    N = 3
    SUM_MIN = 5
    SUM_MAX = 6
    #  X(1)  X(2)  X(3)     SUM
    1    3     1     1        5
    2    2     2     1        5
    3    2     1     2        5
    4    1     3     1        5
    5    1     2     2        5
    6    1     1     3        5
    7    4     1     1        6
    8    3     2     1        6
    9    3     1     2        6
   10    2     3     1        6
   11    2     2     2        6
   12    2     1     3        6
   13    1     4     1        6
   14    1     3     2        6
   15    1     2     3        6
   16    1     1     4        6
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of components in the vector.
    Input, int SUM_MIN, SUM_MAX, the minimum and maximum sums..
    Input/output, int X(N).  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Input/output, int *BASE, a value controlled by this function.
    The user must declare this variable in the calling program,
    and should pass the output value from one call as input to the next,
    but should not alter its value.
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dtpi2dtpipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	const register dim_typ sum_min = s_data->a2;
	const register dim_typ sum_max = s_data->a3;
	int * base = s_data->a4;
	bool * more = s_data->a5;
	
    dim_typ i, j;
    /*
    Initialization.
    */
    if ( !( *more ) )
    {

        if ( sum_max < n || sum_max < sum_min )
        {
            *more = 0;
            return NULL;
        }

        *more = true;

        *base = sum_min;
        if ( *base < n )
            *base = n;

        x[0] = *base - n + 1;
        for ( i = 1; i < n; ++i )
            x[i] = 1;
        return NULL;
    }
    /*
    Next element.
    */
    else
    {
        /*
        Search from the right, seeking an index I < N for which 1 < X(I).
        */
        for ( i = n-2; 0 <= i; --i )
        {
            /*
            If you find such an I, decrease X(I) by 1, and add that to X(I+1).
            */
            if ( 1 < x[i] )
            {
                -- x[i];
                ++ x[i+1];
                /*
                Now grab all the "excess" 1's from the entries to the right of X(I+1).
                */
                for ( j = i+2; j < n; ++j )
                    if ( 1 < x[j] )
                    {
                        x[i+1] += x[j] - 1;
                        x[j] = 1;
                    }
                    return NULL;
            }
        }
        /*
        The current vector is (1,1,1,...BASE-N+1).
        If BASE < SUM_MAX, then increase BASE by 1, and start the new series.
        */
        if ( *base < sum_max )
        {
            ++ *base;
            x[0] = *base - n + 1;
            for ( i = 1; i < n; ++i )
                x[i] = 1;
            return NULL;
        }
        /*
        We returned the last legal vector on the previouis call.
        The calculation is done.
        */
        *more = false;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _vector_constrained_next6 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VECTOR_CONSTRAINED_NEXT6 returns the "next" constrained vector.
  Discussion:
    This routine is similar to VECTOR_CONSTRAINED_NEXT2,
    VECTOR_CONSTRAINED_NEXT3, and VECTOR_CONSTRAINED_NEXT4.
    We consider all vectors X of dimension N whose components
    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
    We are only interested in the subset of these vectors which
    satisfy the following constraint:
      Q_MIN <= sum ( 1 <= I <= N ) ALPHA(I) * X(I) <= Q_MAX
    This routine returns, one at a time, and in right-to-left
    lexicographic order, exactly those vectors which satisfy
    the constraint.
  Example:
    N = 3
    ALPHA    4.0  3.0  5.0
    Q_MIN   16.0
    Q_MAX   20.0
    X_MIN:   1   1   1
    X_MAX:   5   6   4
    #  X(1)  X(2)  X(3)     Total
    1    2     1     1       20.0
    2    2     2     1       19.0
    3    1     3     1       18.0
    4    1     1     2       17.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of components in the vector.
    Input, double ALPHA[N], the coefficient vector.
    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
    values allowed in each component.
    Input/output, int X[N].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Input, double Q_MIN, Q_MAX, the lower and upper
    limits on the sum.
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dtpit3pi2itpb * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	int * x_min = s_data->a2;
	int * x_max = s_data->a3;
	int * x = s_data->a4;
	const register ityp q_min = s_data->a5;
	const register ityp q_max = s_data->a6;
	bool * more = s_data->a7;
	
    dim_typ i, j;
    ityp total;

    if ( ! ( *more ) )
    {
        *more = true;
        for ( i = 0; i < n; ++i )
            x[i] = x_min[i];

        total = 0.00;
        for ( i = 0; i < n; ++i )
        total += alpha[i] * ( ityp ) ( x[i] );

        if ( q_min <= total && total <= q_max )
            return NULL;
    }

    for ( ; ; )
    {
        j = n - 1;

        for ( ; ; )
        {
            if ( x[j] < x_max[j] )
                break;

            if ( j <= 0 )
            {
                *more = 0;
                return NULL;
            }
                -- j;
        }

        ++ x[j];
        for ( i = j + 1; i < n; ++i )
            x[i] = x_min[i];

        total = 0.00;
        for ( i = 0; i < n; ++i )
        total += alpha[i] * ( ityp ) ( x[i] );

        if ( q_min <= total && total <= q_max )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _vector_constrained_next7 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VECTOR_CONSTRAINED_NEXT7 returns the "next" constrained vector.
  Discussion:
    We consider vectors X of dimension N satisfying:
      0 <= X(1:N) <= X_MAX(1:N).
    We are only interested in the subset of these vectors which
    satisfy the following constraint:
      Q_MIN < sum ( 1 <= I <= N ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
    This routine returns, one at a time, and in right-to-left
    lexicographic order, exactly those vectors which satisfy
    the constraint.
  Example:
    N = 3
    LEVEL_WEIGHT    4.0  3.0  5.0
    Q_MIN   16.0
    Q_MAX   20.0
    X_MAX:   5   6   4
    #  X(1)  X(2)  X(3)     Total
    1    2     1     1       20.0
    2    2     2     1       19.0
    3    1     3     1       18.0
    4    1     1     2       17.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of components in the vector.
    Input, double LEVEL_WEIGHT[N], the coefficient vector.
    Input, int X_MAX[N], the maximum values allowed in each component.
    Input/output, int X[N].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Input, double Q_MIN, Q_MAX, the lower and upper
    limits on the sum.
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dtpit2pi2itpb * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * level_weight = s_data->a1;
	int * x_max = s_data->a2;
	int * x = s_data->a3;
	const register ityp q_min = s_data->a4;
	const register ityp q_max = s_data->a5;
	bool * more = s_data->a6;
	
    dim_typ i, j;
    ityp total;

    if ( ! ( *more ) )
    {
        *more = true;
        for ( i = 0; i < n; ++i)
        x[i] = 0;


        total = 0.00;
        for ( i = 0; i < n; ++i )
            total += level_weight[i] * ( ityp ) ( x[i] );

        if ( q_min < total && total <= q_max )
            return NULL;
    }

    for ( ; ; )
    {
        j = n - 1;

        for ( ; ; )
        {
            if ( x[j] < x_max[j] )
                break;

            if ( j <= 0 )
            {
                *more = 0;
                return NULL;
            }
                -- j;
        }

        ++ x[j];
        for ( i = j + 1; i < n; ++i )
            x[i] = 0;

        total = 0.00;
        for ( i = 0; i < n; ++i)
            total += level_weight[i] * ( ityp ) ( x[i] );

        if ( q_min < total && total <= q_max )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _vector_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    VECTOR_NEXT returns the "next" integer vector between two ranges.
  Discussion:
    We consider all integer vectors of dimension N satisfying:
      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
    This routine returns, one at a time, and in right-to-left
    lexicographic order, all these vectors.
  Example:
    N = 3
    X_MIN:   2   2   0
    X_MAX:   4   3   1
    #  X(1)  X(2)  X(3)
    1    2     2     0
    2    3     2     0
    3    4     2     0
    4    2     3     0
    5    3     3     0
    6    4     3     0
    7    2     2     1
    8    3     2     1
    9    4     2     1
   10    2     3     1
   11    3     3     1
   12    4     3     1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of components in the vector.
    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
    values allowed in each component.
    Input/output, int X[N].  On first call, with
    MORE = FALSE, the input value of X is not important.  On subsequent calls,
    the input value of X should be the output value from the previous call.
    On output, with MORE = TRUE, the value of X will be the "next"
    vector in the reverse lexicographical list of vectors.  However, on
    output with MORE = FALSE, the vector X is meaningless, because there
    are no more vectors in the list.
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dt3pipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x_min = s_data->a1;
	int * x_max = s_data->a2;
	int * x = s_data->a3;
	bool * more = s_data->a4;
	
    dim_typ i;

    if ( !( *more ) )
    {
        for ( i = 0; i < n; ++i )
            x[i] = x_min[i];
        *more = true;
    }
    else
    {
        i = 0;

        for ( ; ; )
        {
            if ( x[i] < x_max[i] )
            {
                ++ x[i];
                break;
            }

            x[i] = x_min[i];

            ++ i;

            if ( n <= i )
            {
                *more = false;
                break;
            }
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _ytb_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    YTB_ENUM enumerates the Young tableau of size N.
  Discussion:
    If A(N) is the number of Young tableau of size N, then A(1) = 1,
    A(2) = 2, and
    A(N) = A(N-1) + (N-1) * A(N-2).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the integer which is to be partitioned.
    Output, int YTB_ENUM, the number of Young tableau of N.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ a1;
    dim_typ a2;
    dim_typ a3;
    dim_typ i;
    dim_typ num;

    if ( n == 0 )
        num = 0;
    else if ( n == 1 )
        num = 1;
    else if ( n == 2 )
        num = 2;
    else
    {
        a2 = 1;
        a3 = 2;
        for ( i = 3; i <= n; ++i )
        {
            a1 = a2;
            a2 = a3;
            a3 = a2 + ( i - 1 ) * a1;
        }
        num = a3;
    }

	result = num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ytb_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    YTB_NEXT computes the next Young tableau for a given shape.
  Discussion:
    When the routine is called with MORE = .FALSE. (the first time), and
    if LAMBDA on this call has M parts, with M<N, then the user
    must also make sure that LAMBDA(M+1) = 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2008
  Author:
    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the integer which is to be partitioned.
    Output, int LAMBDA[N], contains a partition of N, that is,
    the entries are positive integers that sum to N.
    Output, int A[N].  A[I] is the row containing I
    in the output tableau.
    Input/output, int *MORE.  Set MORE FALSE before first call.
    It is reset to TRUE as the program returns a new tableau
    on each call, until the last tableau is computed, when
    the program also sets MORE = FALSE.
*/
{
	const dt2pipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * lambda = s_data->a1;
	int * a = s_data->a2;
	bool * more = s_data->a3;
	
    dim_typ i;
    dim_typ ir;
    dim_typ it;
    dim_typ j;
    dim_typ k;
    dim_typ isave;

    it = n;

    if ( *more )
    {
        lambda[0] = 1;
        for ( i = 1; i < n; ++i )
            lambda[i] = 0;

        isave = 0;

        for ( i = 2; i <= n; ++i )
        {
            ++ lambda[a[i-1]-1];

            if ( a[i-1] < a[i-2] )
            {
                isave = i;
                break;
            }

        }

        if ( isave == 0 )
        {
            *more = false;
            return NULL;
        }

        it = lambda[a[isave-1]];

        for ( i = n; 1 <= i; --i)
        {
            if ( lambda[i-1] == it )
            {
                a[isave-1] = i;
                -- lambda[i-1];
                it = isave - 1;
                break;
            }
        }
    }

    k = ir = 1;

    for ( ; ; )
    {
        if ( n < ir )
            break;

        if ( lambda[ir-1] != 0 )
        {
            a[k-1] = ir;
            -- lambda[ir-1];
            ++ k;
            ++ ir;
            continue;
        }

        if ( it < k )
            break;

        ir = 1;

    }

    if ( n == 1 )
    {
        *more = false;
        return NULL;
    }

    for ( j = 1; j < n; ++j )
        if ( a[j] < a[j-1] )
        {
            *more = true;
            return NULL;
        }

    *more = false;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ytb_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    YTB_RANDOM selects a random Young tableau of a given shape.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2008
  Author:
    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the integer which has been partitioned.
    Input, int LAMBDA[N], is a partition of N, that is,
    N = sum ( 0 <= I < N ) LAMBDA[I].
    Input/output, int *SEED, a seed for the random number generator.
    Output, int A[N], the vector describing the Young tableau.
*/
{
	const dt3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * lambda = s_data->a1;
	int * seed = s_data->a2;
	int * a = s_data->a3;
	
    dim_typ i;
    dim_typ ih;
    dim_typ j;
    dim_typ k;
    dim_typ m;

    for ( i = 0; i < n; ++i )
        a[i] = 0;

    i = k = 0;

    for ( ; ; )
    {
        ++ i;

        for ( j = 0; j < lambda[i-1]; ++j )
        {
            ++ a[j];
            ++ k;
        }

        if ( n <= k )
            break;
    }

    for ( m = 1; m <= n; ++m )
    {

        for ( ; ; )
        {
            i = i4_uniform_ab ( 1, a[0], seed );
            j = i4_uniform_ab ( 1, lambda[0], seed );

            if ( i <= a[j-1] && j <= lambda[i-1] )
                break;
        }

        for ( ; ; )
        {
            ih = a[j-1] + lambda[i-1] - i - j;

            if ( ih == 0 )
                break;

            k = i4_uniform_ab ( 1, ih, seed );

            if ( k <= lambda[i-1] - j )
                j += k;
            else
                i += k - lambda[i-1] + j;
        }

        -- lambda[i-1];
        -- a[j-1];
        a[n-m] = i;

    }

    for ( i = 0; i < n; ++i )
        ++ lambda[a[i]-1];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _simple_f ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLE_F returns the right hand side of the three body ODE system.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, double T, the value of the independent variable.
    Input, double Y[NEQN], the value of the dependent variable.
    Output, double YP[NEQN], the value of the derivatives.
*/
{
	const it2pit * const s_data = data;
	const register ityp t = s_data->a0;
	ityp * y = s_data->a1;
	ityp * yp = s_data->a2;
	
    ityp m0;
    ityp m1;
    ityp m2;
    ityp n0;
    ityp n1;
    ityp n2;
    ityp x0;
    ityp x1;
    ityp x2;
    ityp y0;
    ityp y1;
    ityp y2;

    m0 = 5.00;
    m1 = 3.00;
    m2 = 4.00;

    x0 = y[0];
    y0 = y[1];

    x1 = y[4];
    y1 = y[5];

    x2 = y[8];
    y2 = y[9];

    n0 = sqrt ( pow ( pow ( x2 - x1, 2 ) + pow ( y2 - y1, 2 ), 3 ) );
    n1 = sqrt ( pow ( pow ( x0 - x2, 2 ) + pow ( y0 - y2, 2 ), 3 ) );
    n2 = sqrt ( pow ( pow ( x1 - x0, 2 ) + pow ( y1 - y0, 2 ), 3 ) );

    yp[0]  =  y[2];
    yp[1]  =  y[3];
    yp[2]  = - m1 * ( x0 - x1 ) / n2 - m2 * ( x0 - x2 ) / n1;
    yp[3]  = - m1 * ( y0 - y1 ) / n2 - m2 * ( y0 - y2 ) / n1;
    yp[4]  =  y[6];
    yp[5]  =  y[7];
    yp[6]  = - m2 * ( x1 - x0 ) / n0 - m0 * ( x1 - x2 ) / n2;
    yp[7]  = - m2 * ( y1 - y0 ) / n0 - m0 * ( y1 - y2 ) / n2;
    yp[8]  = y[10];
    yp[9]  = y[11];
    yp[10] = - m0 * ( x2 - x0 ) / n1 - m1 * ( x2 - x1 ) / n0;
    yp[11] = - m0 * ( y2 - y0 ) / n1 - m1 * ( y2 - y1 ) / n0;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cbt_traverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    CBT_TRAVERSE traverses a complete binary tree of given depth.
  Discussion:
    There will be 2^DEPTH terminal nodes of the complete binary tree.
    This function traverses the tree, and prints out a binary code of 0's
    and 1's each time it encounters a terminal node.  This results in a
    printout of the binary digits from 0 to 2^DEPTH - 1.
    The function is intended as a framework to be used to traverse a binary
    tree.  Thus, in practice, a user would insert some action when a terminal
    node is encountered.
    Another use would occur when a combinatorial search is being made, for
    example in a knapsack problem.  Each binary string then represents which
    objects are to be included in the knapsack.  In that case, the traversal
    could be speeded up by noticing cases where a nonterminal node has been
    reached, but the knapsack is already full, in which case the only solution
    uses none of the succeeding items, or overfull, in which case no solutions
    exist that include this initial path segment.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int DEPTH, the depth of the tree.
*/
{
	const register dim_typ depth = *(dim_typ *) data;
	
    int *b;
    dim_typ direction;
    dim_typ DOWNLEFT = 1;
    dim_typ i;
    dim_typ k;
    dim_typ p;
    dim_typ UP = 3;
    dim_typ UPDOWNRIGHT = 2;

    if ( depth == 0 )
        return NULL;

    b = ( int * ) malloc ( ( depth + 1 ) * sizeof ( int ) );

    for ( i = 0; i <= depth; i++ )
    {
    b[i] = 0;
    }
    p = 0;
    direction = DOWNLEFT;
    k = 0;

    for ( ; ; )
    {
        /*
        Try going in direction DOWNLEFT.
        */
        if ( direction == DOWNLEFT )
        {
            ++ p;
            b[p-1] = 0;
            if ( p < depth );
            else
            {
                ++ k;
                direction = UPDOWNRIGHT;
            }
        }
        /*
        Try going in direction UPDOWNRIGHT.
        */
        if ( direction == UPDOWNRIGHT )
        {
            b[p-1] = + 1;
            if ( p < depth )
                direction = DOWNLEFT;
            else
            {
                ++ k;
                direction = UP;
            }
        }
        /*
        Try going in direction UP.
        */
        if ( direction == UP )
        {
            -- p;
            if ( 1 <= p && b[p-1] == 0 )
                direction = UPDOWNRIGHT;
            else
                break;
        }
    }

    free ( b );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _triangle_order6_physical_to_reference ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE maps a physical point to a reference point.
  Discussion:
    Given the vertices of an order 6 physical triangle and a point
 (X,Y) in the physical triangle, the routine computes the value
    of the corresponding image point (R,S) in reference space.
    The mapping from (R,S) to (X,Y) has the form:
      X(R,S) = A1 * R * R + B1 * R * S + C1 * S * S
             + D1 * R     + E1 * S     + F1

      Y(R,S) = A2 * R * R + B2 * R * S + C2 * S * S
             + D2 * R     + E2 * S     + F2
  Reference Element T3:
    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 December 2006
  Author:
    John Burkardt
  Parameters:
    Input, double T(2,6), the coordinates of the vertices.
    The vertices are assumed to be the images of (0,0), (1,0), (0,1),
 (1/2,0), (1/2,1/2) and (0,1/2), in that order.
    Input, int N, the number of points to transform.
    Input, double PHY(2,N), the coordinates of points in the
    physical space.
    Output, double REF(2,N), the coordinates of the corresponding
    points in the reference space.
*/
{
	const dt3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * phy = s_data->a1;
	ityp * ref = s_data->a2;
	ityp * t = s_data->a3;
	
    ityp a[2];
    ityp b[2];
    ityp c[2];
    ityp d[2];
    ityp det;
    ityp dx[2];
    ityp e[2];
    ityp f[2];
    ityp fun[2];
    ityp fun_norm;
    dim_typ i, j, it;
    ityp jac[4];
    /*
    Set iteration parameters.
    */
    #pragma omp parallel for num_threads(2)
    for ( i = 0; i < 2; ++i )
    {
        a[i] =   2.00 * t[i+0] + 2.00 * t[i+2] - 4.00 * t[i+6];
        b[i] =   4.00 * t[i+0] - 4.00 * t[i+6] + 4.00 * t[i+8] - 4.00 * t[i+10];
        c[i] =   2.00 * t[i+0] + 2.00 * t[i+4] - 4.00 * t[i+10];

        d[i] = - 3.00 * t[i] - t[i+2] + 4.00 * t[i+6];
        e[i] = - 3.00 * t[i] - t[i+4] + 4.00 * t[i+10];

        f[i] =   t[i];
    }
    /*
    Initialize the points by inverting the linear map.
    */
    triangle_order3_physical_to_reference ( t, n, phy, ref );
    /*
    Carry out the Newton iteration.
    */
    for ( j = 0; j < n; ++j )
    {
        for ( it = 0; it < 10; ++it )
        {
			#pragma omp parallel for num_threads(2)
            for ( i = 0; i < 2; ++i )
                fun[i] = a[i] * ref[0+(j<<1)] * ref[0+(j<<1)]+ b[i] * ref[0+(j<<1)] * ref[1+(j<<1)]+ c[i] * ref[1+(j<<1)] * ref[1+(j<<1)]+ d[i] * ref[(j<<1)]+ e[i] * ref[1+(j<<1)]+ f[i]- phy[i+(j<<1)];

            fun_norm = sqrt ( pow ( fun[0], 2 ) + pow ( fun[1], 2 ) );

            if ( fun_norm <= 0.000001 )
                break;

            jac[0+0*2] = 2.00 * a[0] * ref[(j<<1)] + b[0] * ref[1+(j<<1)] + d[0];
            jac[1+0*2] = 2.00 * a[1] * ref[(j<<1)] + b[1] * ref[1+(j<<1)] + d[1];
            jac[0+1*2] = b[0] * ref[(j<<1)] + 2.00 * c[0] * ref[1+(j<<1)] + e[0];
            jac[3] = b[1] * ref[(j<<1)] + 2.00 * c[1] * ref[1+(j<<1)] + e[1];

            det = jac[0] * jac[3] - jac[2] * jac[1];

            if ( det == 0.00 )
                return NULL;

            dx[0] = (  jac[3] * fun[0] - jac[2] * fun[1] ) / det;
            dx[1] = ( -jac[1] * fun[0] + jac[0] * fun[1] ) / det;

            ref[0+(j<<1)] -= dx[0];
            ref[1+(j<<1)] -= dx[1];
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _triangle_order6_reference_to_physical ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL maps reference points to physical points.
  Discussion:
    Given the vertices of an order 6 physical triangle and a point
 (XSI,ETA) in the reference triangle, the routine computes the value
    of the corresponding image point (X,Y) in physical space.
    The mapping from (XSI,ETA) to (X,Y) has the form:
      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
                 + D1 * XSI    + E1 * ETA     + F1
      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
                 + D2 * XSI    + E2 * ETA     + F2
  Reference Element T6:
    |
    1  3
    |  |\
    |  | \
    S  6  5
    |  |   \
    |  |    \
    0  1--4--2
    |
    +--0--R--1-->
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 June 2005
  Author:
    John Burkardt
  Parameters:
    Input, double T[2*6], the coordinates of the vertices.
    The vertices are assumed to be the images of (0,0), (1,0),
 (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
    Input, integer N, the number of points to transform.
    Input, double REF[2*N], points in the reference triangle.
    Output, double PHY[2*N], corresponding points in the
    physical triangle.
*/
{
	const dt3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * ref = s_data->a1;
	ityp * phy = s_data->a2;
	ityp * t = s_data->a3;
	
    ityp a[2];
    ityp b[2];
    ityp c[2];
    ityp d[2];
    ityp e[2];
    ityp f[2];
    dim_typ i, j;

    #pragma omp parallel for num_threads(2)
    for ( i = 0; i < 2; ++i )
    {
        a[i] =   2.00 * t[i] + 2.00 * t[i+2]- 4.00 * t[i+6];
        b[i] =   4.00 * t[i]- 4.00 * t[i+6] + 4.00 * t[i+8] - 4.00 * t[i+10];
        c[i] =   2.00 * t[i]                  + 2.00 * t[i+4]- 4.00 * t[i+10];
        d[i] = - 3.00 * t[i] -       t[i+2]+ 4.0 * t[i+3*2];
        e[i] = - 3.00 * t[i]                  -       t[i+4]+ 4.00 * t[i+10];
        f[i] =         t[i];

    }

    for ( j = 0; j < n; ++j )
        #pragma omp parallel for num_threads(2)
        for ( i = 0; i < 2; ++i )
            phy[i+(j<<1)] = a[i] * ref[(j<<1)] * ref[(j<<1)]+ b[i] * ref[(j<<1)] * ref[1+(j<<1)]+ c[i] * ref[1+(j<<1)] * ref[1+(j<<1)]+ d[i] * ref[(j<<1)]+ e[i] * ref[1+(j<<1)]+ f[i];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _vand1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VAND1 returns the Vandermonde1 matrix A with 1's on the first row.
  Formula:
    A(I,J) = X(J)^(I-1)
  Example:
    N = 5, X = ( 2, 3, 4, 5, 6 )
    1  1   1   1   1
    2  3   4   5   6
    4  9  16  25  36
    8 27  64 125  216
   16 81 256 625 1296
  Properties:
    A is generally not symmetric: A' /= A.
    A is nonsingular if, and only if, the X values are distinct.
    det ( A ) = product ( 1 <= I <= N ) ( 1 <= J .lt. I ) ( X(I) - X(J) ).
             = product ( 1 <= J <= N ) X(J)
             * product ( 1 <= I .lt. J ) ( X(J) - X(I) ).
    A is generally ill-conditioned.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 February 2014
  Author:
    John Burkardt
  Reference:
    Robert Gregory, David Karney,
    A Collection of Matrices for Testing Computational Algorithms,
    Wiley, 1969, page 27,
    LC: QA263.G68.
    Nicholas Higham,
    Stability analysis of algorithms for solving confluent
    Vandermonde-like systems,
    SIAM Journal on Matrix Analysis and Applications,
    Volume 11, 1990, pages 23-41.
  Parameters:
    Input, int N, the order of the matrix desired.
    Input, double X[N], the values that define A.
    Output, double VAND1[N*N], the matrix.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, j;

    ityp *a = (ityp * ) malloc ( n * n * sizeof (ityp) );

    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j )
            a[i+j*n] = (i == 0 && x[j] == 0.00) ? 1.00 : pow ( x[j], i );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _vandermonde_value_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    VANDERMONDE_VALUE_1D evaluates a Vandermonde interpolant.
  Discussion:
    The polynomial
      p(x) = cd0 + cd1 * x + cd2 * x^2 + ... + cd(nd-1) * x^(nd-1)
    is to be evaluated at the vector of NI values XI.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 July 2015
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data values.
    Input, double CD[ND], the polynomial coefficients.
    C[I] is the coefficient of X^I.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double VANDERMONDE_VALUE_1D[NI], the interpolation values.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ ni = s_data->a1;
	ityp * cd = s_data->a2;
	ityp * xi = s_data->a3;
	
    dim_typ i, j;
    ityp *yi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( j = 0; j < ni; ++j )
        yi[j] = cd[nd-1];

    for ( i = nd - 2; 0 <= i; --i )
        for ( j = 0; j < ni; ++j )
            yi[j] *= xi[j] + cd[i];
    return yi;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cg_gb ( void * data)
/******************************************************************************/
/*
  Purpose:
    CG_GB uses the conjugate gradient method for a general banded (GB) matrix.
  Discussion:
    The linear system has the form A*x=b, where A is a positive-definite
    symmetric matrix.
    The method is designed to reach the solution to the linear system
      A * x = b
    after N computational steps.  However, roundoff may introduce
    unacceptably large errors for some problems.  In such a case,
    calling the routine a second time, using the current solution estimate
    as the new starting guess, should result in improved results.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 June 2014
  Author:
    John Burkardt
  Reference:
    Frank Beckman,
    The Solution of Linear Equations by the Conjugate Gradient Method,
    in Mathematical Methods for Digital Computers,
    edited by John Ralston, Herbert Wilf,
    Wiley, 1967,
    ISBN: 0471706892,
    LC: QA76.5.R3.
  Parameters:
    Input, int N, the order of the matrix.
    Input, int ML, MU, the lower and upper bandwidths.
    Input, double A[(2*ML+MU+1)*N], the band matrix.
    Input, double B[N], the right hand side vector.
    Input/output, double X[N].
    On input, an estimate for the solution, which may be 0.
    On output, the approximate solution vector.
*/
{
	const _3dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ ml = s_data->a1;
	const register dim_typ mu = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	ityp * x =  s_data->a5;
	
    ityp alpha;
    ityp *ap;
    ityp beta;
    dim_typ i;
    dim_typ it;
    ityp *p;
    ityp pap;
    ityp pr;
    ityp *r;
    ityp rap;
    /*
    Initialize
    AP = A * x,
    R  = b - A * x,
    P  = b - A * x.
    */
    ap = mv_gb ( n, n, ml, mu, a, x );

    r = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    p = ( ityp *) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
    {
        r[i] = b[i] - ap[i];
        p[i] = b[i] - ap[i];
    }

    /*
    Do the N steps of the conjugate gradient method.
    */
    for ( it = 1; it <= n; ++it )
    {
        /*
        Compute the matrix*vector product AP = A*P.
        */
        free ( ap );
        ap = mv_gb ( n, n, ml, mu, a, p );
        /*
        Compute the dot products
        PAP = P*AP,
        PR  = P*R
        Set
        ALPHA = PR / PAP.
        */
        pap = r8vec_dot_product ( n, p, ap );
        pr = r8vec_dot_product ( n, p, r );

        if ( pap == 0.00 )
            break;
        alpha = pr / pap;
        /*
        Set
        X = X + ALPHA * P
        R = R - ALPHA * AP.
        */
        for ( i = 0; i < n; ++i )
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        /*
        Compute the vector dot product
        RAP = R*AP
        Set
        BETA = - RAP / PAP.
        */
        rap = r8vec_dot_product ( n, r, ap );

        beta = - rap / pap;
        /*
        Update the perturbation vector
        P = R + BETA * P.
        */
        for ( i = 0; i < n; ++i )
            p[i] = r[i] + beta * p[i];
    }

    free ( ap );
    free ( p );
    free ( r );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cg_ge ( void * data)
/******************************************************************************/
/*
  Purpose:
    CG_GE uses the conjugate gradient method for a general storage (GE) matrix.
  Discussion:
    The linear system has the form A*x=b, where A is a positive-definite
    symmetric matrix, stored as a full storage matrix.
    The method is designed to reach the solution to the linear system
      A * x = b
    after N computational steps.  However, roundoff may introduce
    unacceptably large errors for some problems.  In such a case,
    calling the routine a second time, using the current solution estimate
    as the new starting guess, should result in improved results.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 June 2014
  Author:
    John Burkardt
  Reference:
    Frank Beckman,
    The Solution of Linear Equations by the Conjugate Gradient Method,
    in Mathematical Methods for Digital Computers,
    edited by John Ralston, Herbert Wilf,
    Wiley, 1967,
    ISBN: 0471706892,
    LC: QA76.5.R3.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the matrix.
    Input, double B[N], the right hand side vector.
    Input/output, double X[N].
    On input, an estimate for the solution, which may be 0.
    On output,  the approximate solution vector.
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x =  s_data->a3;
	
    ityp alpha;
    ityp *ap;
    ityp beta;
    dim_typ i;
    dim_typ it;
    ityp *p;
    ityp pap;
    ityp pr;
    ityp *r;
    ityp rap;
    /*
    Initialize
    AP = A * x,
    R  = b - A * x,
    P  = b - A * x.
    */
    ap = mv_ge ( n, n, a, x );

    r = ( ityp * ) malloc ( n * sizeof ( ityp) );
    p = ( ityp *) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
    {
        r[i] = b[i] - ap[i];
        p[i] = b[i] - ap[i];
    }

    /*
    Do the N steps of the conjugate gradient method.
    */
    for ( it = 1; it <= n; it++ )
    {
    /*
    Compute the matrix*vector product AP = A*P.
    */
    free ( ap );
    ap = mv_ge ( n, n, a, p );
    /*
    Compute the dot products
    PAP = P*AP,
    PR  = P*R
    Set
    ALPHA = PR / PAP.
    */
    pap = r8vec_dot_product ( n, p, ap );
    pr = r8vec_dot_product ( n, p, r );

    if ( pap == 0.00 )
        break;

    alpha = pr / pap;
    /*
    Set
    X = X + ALPHA * P
    R = R - ALPHA * AP.
    */
    for ( i = 0; i < n; ++i )
    {
        x[i] += alpha * p[i];
        r[i] -= alpha * ap[i];
    }
    /*
    Compute the vector dot product
    RAP = R*AP
    Set
    BETA = - RAP / PAP.
    */
    rap = r8vec_dot_product ( n, r, ap );

    beta = - rap / pap;
    /*
    Update the perturbation vector
    P = R + BETA * P.
    */
    for ( i = 0; i < n; ++i )
        p[i] = r[i] + beta * p[i];
    }

    free ( ap );
    free ( p );
    free ( r );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rnorm ( void * data)
/******************************************************************************/
/*
  Purpose:
    RNORM returns two independent standard random normal deviates.
  Discussion:
    This routine sets U1 and U2 to two independent standardized
    random normal deviates.   This is a version of the
    method given in Knuth.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 April 2014
  Author:
    Original FORTRAN77 version by William Smith, Ronald Hocking.
    This C version by John Burkardt.
  Reference:
    Donald Knuth,
    The Art of Computer Programming,
    Volume 2, Seminumerical Algorithms,
    Third Edition,
    Addison Wesley, 1997,
    ISBN: 0201896842,
    LC: QA76.6.K64.
  Parameters:
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double *U1, *U2, two standard random normal deviates.
*/
{
	static ityp result = MAX_VAL;
	
	const pi2pit * const s_data = data; 
	int * seed = s_data->a0;
	ityp * u1 = s_data->a1;
	ityp * u2 = s_data->a2;
	
	ityp s, x, y;

	for ( ; ; )
	{
		x = r8_uniform_01 ( seed );
		y = r8_uniform_01 ( seed );
		x = 2.00 * x - 1.00;
		y = 2.00 * y - 1.00;
		s = x * x + y * y;

		if ( s <= 1.00 )
		{
			s = sqrt ( - 2.00 * log ( s ) / s );
			*u1 = x * s;
			*u2 = y * s;
			break;
		}
	}

	result = s;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _wshrt ( void * data)
/******************************************************************************/
/*
  Purpose:
    WSHRT returns a random Wishart variate.
  Discussion:
    This routine is a Wishart variate generator.
    On output, SA is an upper-triangular matrix of size NP * NP,
    written in linear form, column ordered, whose elements have a
    Wishart(N, SIGMA) distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 April 2014
  Author:
    Original FORTRAN77 version by William Smith, Ronald Hocking.
    This C version by John Burkardt.
  Reference:
    William Smith, Ronald Hocking,
    Algorithm AS 53, Wishart Variate Generator,
    Applied Statistics,
    Volume 21, Number 3, pages 341-345, 1972.
  Parameters:
    Input, double D[NP*(NP+1)/2], the upper triangular array that
    represents the Cholesky factor of the correlation matrix SIGMA.
    D is stored in column-major form.
    Input, int N, the number of degrees of freedom.
    1 <= N <= NP.
    Input, int NP, the size of variables.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double WSHART[NP*(NP+1)/2], a sample from the
    Wishart distribution.
*/
{
	const dtpitdtpipit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * d = s_data->a1;
	const register dim_typ np = s_data->a2;
	int * seed = s_data->a3;
	ityp * sa = s_data->a4;
	
	ityp c;
	ityp df;
	dim_typ i;
	dim_typ ii;
	dim_typ ip;
	dim_typ j;
	dim_typ k;
	dim_typ nnp;
	dim_typ nq;
	dim_typ nr;
	dim_typ ns;
	ityp rn;
	ityp u1;
	ityp u2;

	k = 0;
	nnp = ( np * ( np + 1 ) ) / 2;
	/*
	Load SB with independent normal (0, 1) variates.
	*/

	ityp sb[nnp];

	while ( k < nnp )
	{
		rnorm ( seed, &u1, &u2 );
		sb[k] = u1;
		++ k;

		if ( k < nnp )
		{
			sb[k] = u2;
			++ k;
		}
	}
	/*
	Load diagonal elements with square root of chi-square variates.
	*/
	ns = 0;

	for ( i = 1; i <= np; ++i )
	{
		df = (ityp) (np - i + 1);
		ns += i;
		u1 = 2.00 / ( 9.00 * df );
		u2 = 1.00 - u1;
		u1 = sqrt ( u1 );
		/*
		Wilson-Hilferty formula for approximating chi-square variates:
		The original code did not take the absolute value!
		*/
		sb[ns-1] = sqrt ( df * fabs ( pow ( u2 + sb[ns-1] * u1, 3 ) ) );
	}

	rn = (ityp)(n);
	nr = 1;

	for ( i = 1; i <= np; ++i )
	{
		nr += i - 1;
		for ( j = i; j <= np; ++j )
		{
			ip = nr;
			nq = ( j * ( j - 1 ) ) / 2 + i - 1;
			c = 0.00;
			for ( k = i; k <= j; k++ )
			{
				ip += k - 1;
				++ nq;
				c += sb[ip-1] * d[nq-1];
			}
			sa[ip-1] = c;
		}
	}

	for ( i = 1; i <= np; ++i )
	{
		ii = np - i + 1;
		nq = nnp - np;
		for ( j = 1; j <= i; ++j )
		{
			ip = ( ii * ( ii - 1 ) ) / 2;
			c = 0.00;
			for ( k = i; k <= np; ++k )
			{
				++ ip;
				++ nq;
				c += sa[ip-1] * sa[nq-1];
			}
			sa[nq-1] = c / rn;
			nq = nq - (np<<1) + i + j - 1;
		}
	}
	return NULL;
}

#endif
