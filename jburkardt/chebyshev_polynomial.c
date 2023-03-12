#ifndef __DISABLEDEEP_CHEBYSHEVPOLYNOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_uniform ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_UNIFORM returns a scaled pseudorandom I4.
  Discussion:
    The pseudorandom number should be uniformly distributed
    between A and B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 November 2006
  Author:
    John Burkardt
  Reference:
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.
    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.
    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.
    Peter Lewis, Allen Goodman, James Miller
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.
  Parameters:
    Input, int A, B, the limits of the interval.
    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.
    Output, int I4_UNIFORM, a number between A and B.
*/
{
	static int result = INT_MAX;
	
	const _2ipi * const s_data = data;
	const register int a = s_data->a0;
	const register int b = s_data->a1;
	int * seed = s_data->a2;
	
    int k;
    ityp r;
    int value;

    if ( *seed == 0 )
    {
    	result = INT_MAX;
        return &result;
    }

    k = *seed / 127773;
    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
        *seed += 2147483647;

    r = ( ityp ) ( *seed ) * 4.656612875E-10;
    /*
    Scale R to lie between A-0.5 and B+0.5.
    */
    r = ( 1.00 - r ) * ( ( ityp ) ( MIN ( a, b ) ) - 0.50 )+         r   * ( ( ityp ) ( MAX ( a, b ) ) + 0.5 );
    /*
    Use rounding to convert R to an integer between A and B.
    */
    value = r8_nint ( r );
    value = MAX ( value, MIN ( a, b ) );
    value = MIN ( value, MAX ( a, b ) );
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_copy_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_COPY_NEW copies one r8MAT to a "new" r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8's, which
    may be stored as a vector in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A1[M*N], the matrix to be copied.
    Output, ityp r8MAT_COPY_NEW[M*N], the copy of A1.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a1 = s_data->a2;
	
    ityp *a2;
    dim_typ i, j;

    a2 = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i)
            a2[i+j*m] = a1[i+j*m];

    return a2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_uniform_01_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_UNIFORM_01_NEW returns a unit pseudorandom r8VEC.
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
    Output, ityp r8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
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
        r[i] = ( ityp ) ( *seed ) * 4.656612875E-10;
    }

    return r;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_double_product_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)/sqrt(1-x^2) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the polynomial indices.
    0 <= I, J.
    Output, double T_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * dim = data;
	result = ((dim[XROW] < 0 || dim[YROW] < 0) ? CHEBYSHEV_INVALIDRETURNVALUE : (dim[XROW] == dim[YROW] ? M_PI/(1+(0<dim[XROW])) : 0.00));
	return &result;
}



/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_INTEGRAL: integral ( -1 <= x <= +1 ) x^e dx / sqrt ( 1 - x^2 ).
  Discussion:
    Set
      x = cos ( theta ),
      dx = - sin ( theta ) d theta = - sqrt ( 1 - x^2 ) d theta
    to transform the integral to
      integral ( 0 <= theta <= M_PI ) - ( cos ( theta ) )^e d theta
    which becomes
      0 if E is odd,
   (1/2^e) * choose ( e, e/2 ) * M_PI if E is even.
  Licensing:
    This code is distributed under the GNU LGPL license.

  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int E, the exponent of X.
    0 <= E.
    Output, double T_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ e = *(dim_typ *) data;
	result = (e%2) == 1 ? 0.00 : r8_choose(e, e/2)*M_PI/pow(2.00,e);
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_polynomial( void * data)
/******************************************************************************/
/*
  Purpose:
    T_POLYNOMIAL evaluates Chebyshev polynomials T(n,x).
  Discussion:
    Chebyshev polynomials are useful as a basis for representing the
    approximation of functions since they are well conditioned, in the sense
    that in the interval [-1,1] they each have maximum absolute value 1.
    Hence an error in the value of a coefficient of the approximation, of
    size epsilon, is exactly reflected in an error of size epsilon between
    the computed approximation and the theoretical approximation.
    Typical usage is as follows, where we assume for the moment
    that the interval of approximation is [-1,1].  The value
    of N is chosen, the highest polynomial to be used in the
    approximation.  Then the function to be approximated is
    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
    Chebyshev polynomial.  Let these values be denoted by F(XJ).
    The coefficients of the approximation are now defined by
      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
    except that C(0) is given a value which is half that assigned
    to it by the above formula,
    and the representation is
    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
    Now note that, again because of the fact that the Chebyshev polynomials
    have maximum absolute value 1, if the higher order terms of the
    coefficients C are small, then we have the option of truncating
    the approximation by dropping these terms, and we will have an
    exact value for maximum perturbation to the approximation that
    this will cause.
    It should be noted that typically the error in approximation
    is dominated by the first neglected basis function (some multiple of
    T(N+1,X) in the example above).  If this term were the exact error,
    then we would have found the minimax polynomial, the approximating
    polynomial of smallest maximum deviation from the original function.
    The minimax polynomial is hard to compute, and another important
    feature of the Chebyshev approximation is that it tends to behave
    like the minimax polynomial while being easy to compute.
    To evaluate a sum like
      sum ( 0 <= J <= N ) C(J) T(J,X),
    Clenshaw's recurrence formula is recommended instead of computing the
    polynomial values, forming the products and summing.
    Assuming that the coefficients C(J) have been computed
    for J = 0 to N, then the coefficients of the representation of the
    indefinite integral of the function may be computed by
      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1,
    with
      C(N+1)=0
      B(0) arbitrary.
    Also, the coefficients of the representation of the derivative of the
    function may be computed by:
      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0,
    with
      D(N+1) = D(N)=0.
    Some of the above may have to adjusted because of the irregularity of C(0).
    The formula is:
      T(N,X) = COS(N*acos(X))
  Differential equation:
 (1-X*X) Y'' - X Y' + N N Y = 0
  First terms:
    T(0,X) =  1
    T(1,X) =  1 X
    T(2,X) =  2 X^2 -   1
    T(3,X) =  4 X^3 -   3 X
    T(4,X) =  8 X^4 -   8 X^2 +  1
    T(5,X) = 16 X^5 -  20 X^3 +  5 X
    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
  Inequality:
    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
  Orthogonality:
    For integration over [-1,1] with weight
      W(X) = 1 / sqrt(1-X*X),
    if we write the inner product of T(I,X) and T(J,X) as
      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
    then the result is:
      0 if I /= J
      M_PI/2 if I == J /= 0
      M_PI if I == J == 0
    A discrete orthogonality relation is also satisfied at each of
    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
                              = 0 if I /= J
                              = N/2 if I == J /= 0
                              = N if I == J == 0
  Recursion:
    T(0,X) = 1,
    T(1,X) = X,
    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
  Special values:
    T(N,1) = 1
    T(N,-1) = (-1)^N
    T(2N,0) = (-1)^N
    T(2N+1,0) = 0
    T(N,X) = (-1)**N * T(N,-X)
  Zeroes:
    M-th zero of T(N,X) is cos((2*M-1)*M_PI/(2*N)), M = 1 to N
  Extrema:
    M-th extremum of T(N,X) is cos(M_PI*M/N), M = 0 to N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest polynomial to compute.
    Input, double X[M], the evaluation points.
    Output, double T_POLYNOMIAL[M*(N+1)], the values of the Chebyshev polynomials.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	if(n<0)
	{
		result = CHEBYSHEV_INVALIDRETURNVALUE;
		return &result;
	}
	
	if(!n)
	{
		result = 1.00;
		return &result;
	}

	ityp v[n+1];

	v[0] = 1.00;
	v[1] = x;

	for (dim_typ i = 2; i <= n; ++i)
		v[i] = 2.00 * x * v[i-1] - v[i-2];

	result = v[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_polynomial_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_POLYNOMIAL_AB: Chebyshev polynomials T(n,x) in [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, the domain of definition.
    Input, int M, the number of evaluation points.
    Input, int N, the highest polynomial to compute.
    Input, double X[M], the evaluation points.
    It must be the case that A <= X(*) <= B.
    Output, double T_POLYNOMIAL_AB[M*(N+1)], the values.
*/
{
	static ityp result = MAX_VAL;
	
	const dt3it * const s_data = data;
	
	const register dim_typ n = s_data->a2;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register ityp x =  s_data->a3;
	
	result = t_polynomial(n, ((b-x)-(x-a))/(b-a));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_polynomial_coefficients( void * data)
/******************************************************************************/
/*
  Purpose:
    T_POLYNOMIAL_COEFFICIENTS: coefficients of the Chebyshev polynomial T(n,x).
  First terms:
    N/K     0     1      2      3       4     5      6    7      8    9   10
     0      1
     1      0     1
     2     -1     0      2
     3      0    -3      0      4
     4      1     0     -8      0       8
     5      0     5      0    -20       0    16
     6     -1     0     18      0     -48     0     32
     7      0    -7      0     56       0  -112      0    64
  Recursion:
    T(0,X) = 1,
    T(1,X) = X,
    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input, int N, the highest order polynomial to compute.
    Note that polynomials 0 through N will be computed.
    Output, double T_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients
    of the Chebyshev T polynomials.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * c = s_data->a1;

	if(n<0)
	{
		result = CHEBYSHEV_DOMAIN_ERROR;
		return &result;
	}

	dim_typ i, j;

	#pragma omp parallel for
	for(i = 0; i <= n; ++i)
		#pragma omp parallel for
		for (j = 0; j <= n; ++j)
			*(c+i*(n+1)+j) = 0.00;

	*c = 1.00;

	if(!n)
	{
		result = CHEBYSHEV_SUCCESS;
		return &result;
	}

	*(c+(n+1)+1) = 1.00;

	for ( i = 2; i <= n; i++ )
	{
		*(c+i*(n+1)) = - *(c+(i-2)*(n+1));
		for ( j = 1; j <= i - 2; j++ )
			*(c+i*(n+1)+j) = 2.00 * *(c+(i-1)*(n+1)+(j-1)) - *(c+(i-2)*(n+1)+j);
		*(c+i*(n+1)+(i-1)) = 2.00 * *(c+(i-1)*(n+1)+(i-2));
		*(c+i*(n+1)+i) = 2.00 * *(c+(i-1)*(n+1)+(i-1));
	}

	result = CHEBYSHEV_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial T(n,x).
  Discussion:
    The I-th zero is cos((2*I-1)*M_PI/(2*N)), I = 1 to N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the polynomial.
    Output, double T_POLYNOMIAL_ZEROS[N], the zeroes.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * z = s_data->a1;
	
	#pragma omp parallel for
	for (dim_typ i = 0; i < n; ++i)
		z[i] = cos((ityp)((i<<1)-1)*M_PI/(ityp)((n<<1)));
		
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_project_coefficients ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_PROJECT_COEFFICIENTS: function projected onto Chebyshev polynomials T(n,x).
  Discussion:
    It is assumed that the interval of definition is -1 <= x <= +1.
    Over this interval, f(x) will be well approximated by
      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,x)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the highest order polynomial to compute.
    Input, double F ( double X ), evaluates the function.
    Output, double T_PROJECT_COEFFICIENTS[N+1], the projection coefficients
    of f(x) onto T(0,x) through T(n,x).
*/
{	
	const dtfitpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp (* f)(ityp) = s_data->a1;
	ityp * c = s_data->a2;
	
	ityp y;
	ityp d[n+1];
	dim_typ k;
	#pragma omp parallel for
	for ( k = 0; k <= n; ++k)
	{
		y = cos ( M_PI * ( ( ityp ) (k)+0.5)/(ityp)(n+1));
		d[k] = f(y);
	}

	ityp total;
	const register ityp fac = 2.00 / (ityp)(n+1);

	dim_typ j;

	#pragma omp parallel for
	for (j = 0; j <= n; ++j)
	{
		total = 0.00;
		#pragma omp parallel for
		for (k = 0; k <= n; ++k)
			total += d[k]*cos((M_PI*(ityp)(j))*(((ityp)(k)+0.50)/(ityp)(n+1)));
		c[j] = fac * total;
	}

	c[0] /= 2.00;
 
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_project_coefficients_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_PROJECT_COEFFICIENTS_AB: function projected onto T(n,x) over [a,b]
  Discussion:
    It is assumed that the interval of definition is a <= x <= b.
    Over this interval, f(x) will be well approximated by
      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,(2x-a-b)/(b-a))
    where x* = ( x - b - a ) / ( b - a )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the highest order polynomial to compute.
    Input, double F ( double X ), evaluates the function.
    Input, double A, B, the interval of definition.
    Output, double T_PROJECT_COEFFICIENTS_AB[N+1], the projection coefficients
    of f(x) onto T(0,x) through T(n,x).
*/
{
	const dtfit2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp (* f)(ityp) = s_data->a1;
	const register ityp a = s_data->a2;
	const register ityp b = s_data->a3;
	ityp * c = s_data->a4;
	
	
	ityp t, y;
	ityp d[n+1];

	dim_typ k;
	#pragma omp parallel for
	for(k = 0; k <= n; ++k )
	{
		t = cos(M_PI*((ityp)(k)+0.50)/(ityp)(n+1));
		y = ((1.00+t)*b+(1.00-t)*a)/2.00;
		d[k] = f(y);
	}

	ityp total;
	const register ityp fac = 2.00/(ityp)(n+1);

	dim_typ j;

	#pragma omp parallel for
	for(j = 0; j <= n; ++j)
	{
		total = 0.00;
		#pragma omp parallel for
		for(k = 0; k <= n; ++k)
			total += d[k]*cos((M_PI*(ityp)(j))*(((ityp)(k)+0.50)/(ityp)(n+1)));
		c[j] = fac * total;
	}

	c[0] /= 2.00;

	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_project_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_PROJECT_VALUE evaluates an expansion in Chebyshev polynomials T(n,x).
  Discussion:
    The projection is assumed to be based on the interval [-1,+1].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to compute.
    Input, double X[M], the evaluation points.
    Input, double C[N+1], the expansion coefficients.
    Output, double T_PROJECT_VALUE[M], the value of the Chebyshev function.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * c = s_data->a1;
	const register ityp x = s_data->a2;
	
	ityp b0, b1, b2;
	b0 = b1 = 0.00;

	for (dim_typ i=n; 0 <= i ; --i)
	{
		b2 = b1;
		b1 = b0;
		b0 = c[i] + (2.00*x*b1) - b2;
	}

	result = (0.50*(c[0]+b0-b2));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_project_value_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_PROJECT_VALUE_AB evaluates an expansion in Chebyshev polynomials T(n,x).
  Discussion:
    The projection is assumed to be based on the interval [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to compute.
    Input, double X[M], the evaluation points.
    Input, double C[N+1], the expansion coefficients.
    Input, double A, B, the interval of definition.
    Output, double T_PROJECT_VALUE_AB[M], the value of the Chebyshev function.
*/
{
	static ityp result = MAX_VAL;
	
	const dt3itpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;	
	const register ityp x = s_data->a1;
	const register ityp a = s_data->a2;
	const register ityp b = s_data->a3;
	ityp * c = s_data->a4;
	
	ityp b0, b1, b2;
	b0 = b1 = 0.00;

	for(dim_typ i=n; 0 <= i; --i)
	{
		b2 = b1;
		b1 = b0;
		b0 = c[i] + (2.00/(b-a)*(2.00*x-a-b)*b1) -b2;
	}

	result = (0.5*(c[0]+b0-b2));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_quadrature_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_QUADRATURE_RULE: quadrature rule for T(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the rule.
    Output, double T[N], W[N], the points and weights of the rule.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * t = s_data->a1;
	ityp * w = s_data->a2;
	
	dim_typ i;
	ityp bj[n];
	bj[0] = sqrt(0.50);

	#pragma omp parallel for
	for (i = 1; i < n; ++i)
	{
		bj[i] = 0.5;
		t[i] = w[i] = 0.00;
	}

	t[0] = 0.00;
	w[0] = sqrt(M_PI);

	const bool status = imtqlx ( n, t, bj, w );

	#pragma omp parallel for
	for(i = 0; i < n; ++i)
		w[i] *= w[i];

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _t_triple_product_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    T_TRIPLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)*T(k,x)/sqrt(1-x^2) dx
  Discussion:
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Reference:
    John Mason, David Handscomb,
    Chebyshev Polynomials,
    CRC Press, 2002,
    ISBN: 0-8493-035509,
    LC: QA404.5.M37.
  Parameters:
    Input, int I, J, K, the polynomial indices.
    0 <= I, J.
    Output, double T_TRIPLE_PRODUCT_INTEGRAL, the integral.
*/
{ 
	static ityp result = MAX_VAL;

	dim_typ * dim = data;
	result = ((dim[XROW] < 0 || dim[YROW] < 0 || dim[ZROW] < 0) ? CHEBYSHEV_INVALIDRETURNVALUE : 0.50*(t_double_product_integral (dim[XROW]+dim[YROW],dim[ZROW])+t_double_product_integral(dim[XROW]-dim[YROW],dim[ZROW])));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _u_double_product_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    U_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) U(i,x)*U(j,x)*sqrt(1-x^2) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the polynomial indices.
    0 <= I, J.
    Output, double U_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * dim = data;
	result = ((dim[XROW] < 0 || dim[YROW] < 0) ? CHEBYSHEV_INVALIDRETURNVALUE : ((dim[XROW]==dim[YROW])?M_PI_2 : 0.00));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _u_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    U_INTEGRAL: integral ( -1 <= x <= +1 ) x^e sqrt ( 1 - x^2 ) dx.
  Discussion:
     E    U_INTEGRAL
    --    --------------
     0         M_PI /    2
     2         M_PI /    8
     4         M_PI /   16
     6     5 * M_PI /  128
     8     7 * M_PI /  256
    10    21 * M_PI / 1024
    12    33 * M_PI / 2048
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int E, the exponent of X.
    0 <= E.
    Output, double U_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ e = *(dim_typ *) data;
	result = ((e%2) == 1 ? 0.00 : 0.50*sqrt (M_PI)*tgamma(0.50*(ityp)(1+e))/tgamma(2.00+0.50*(ityp)(e)));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _u_polynomial ( void * data)
/******************************************************************************/
/*
  Purpose:
    U_POLYNOMIAL evaluates Chebyshev polynomials U(n,x).
  Differential equation:
 (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
  First terms:
    U(0,X) =   1
    U(1,X) =   2 X
    U(2,X) =   4 X^2 -   1
    U(3,X) =   8 X^3 -   4 X
    U(4,X) =  16 X^4 -  12 X^2 +  1
    U(5,X) =  32 X^5 -  32 X^3 +  6 X
    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
  Recursion:
    U(0,X) = 1,
    U(1,X) = 2 * X,
    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
  Norm:
    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = M_PI/2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest polynomial to compute.
    Input, double X[M], the evaluation points.
    Output, double U_POLYNOMIAL[M*(N+1)], the values of the N+1 Chebyshev polynomials.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;  
	
	if(n<0)
	{
		result = CHEBYSHEV_INVALIDRETURNVALUE;
		return &result;
	}

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	ityp v[n+1];
	v[0] = 1.00;
	v[1] = 2.00 * x;


	for(dim_typ i = 2; i <= n; ++i)
		v[i] = 2.00 * x * v[i-1] - v[i-2];

	result = v[n];
	return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _u_polynomial_coefficients( void * data)
/******************************************************************************/
/*
  Purpose:
    U_POLYNOMIAL_COEFFICIENTS evaluates coefficients of Chebyshev polynomials U(n,x).
  First terms:
    N/K     0     1      2      3       4     5      6    7      8    9   10
     0      1
     1      0     2
     2     -1     0      4
     3      0    -4      0      8
     4      1     0    -12      0      16
     5      0     6      0    -32       0    32
     6     -1     0     24      0     -80     0     64
     7      0    -8      0     80       0  -192      0   128
  Recursion:
    U(0,X) = 1,
    U(1,X) = 2*X,
    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input, int N, the highest order polynomial to compute.
    Note that polynomials 0 through N will be computed.
    Output, double U_POLYNOMIAL_COEFFICIENTS[(N+1)*((N+1)], the coefficients
    of the Chebyshev U polynomials.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * c = s_data->a1;
	
	if(n<0)
	{
		result = CHEBYSHEV_DOMAIN_ERROR;
		return &result;
	}

	dim_typ i, j;

	#pragma omp parallel for
	for(i = 0; i <= n; ++i)
		for(j = 0; j <= n; ++j)
			*(c+i*(n+1)+j) = 0.00;

	*c = 1.00;

	if(!n)
	{
		result = CHEBYSHEV_SUCCESS;
		return &result;
	}

	*(c+(n+1)+1) = 2.00;

	for ( i = 2; i <= n; i++ )
	{
		*(c+i*(n+1)) = -*(c+(i-2)*(n+1));
		for ( j = 1; j <= i-2; j++ )
			*(c+i*(n+1)+j) = 2.00 * *(c+(i-1)*(n+1)+(j-1)) - *(c+(i-2)*(n+1)+j);

		*(c+i*(n+1)+(i-1)) = 2.00 * *(c+(i-1)*(n+1)+(i-2));
		*(c+i*(n+1)+i) = 2.00 * *(c+(i-1)*(n+1)+(i-1));
	}

	result = CHEBYSHEV_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _u_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    U_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials U(n,x).
  Discussion:
    The I-th zero is cos((I-1)*M_PI/(N-1)), I = 1 to N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the polynomial.
    Output, double U_POLYNOMIAL_ZEROS[N], the zeroes.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * z = s_data->a1;
	
	#pragma omp parallel for
	for (dim_typ i = 0; i < n; ++i)
		z[i] = cos((ityp)(i)*M_PI/(ityp)(n+1));
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _u_quadrature_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    U_QUADRATURE_RULE: quadrature rule for U(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the rule.
    Output, double T[N], W[N], the points and weights of the rule.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * t = s_data->a1;
	ityp * w = s_data->a2; 
	
	dim_typ i;
	ityp bj[n];
	#pragma omp parallel for
	for ( i = 0; i < n; i++ )
	{
		bj[i] = 0.5;
		t[i] = w[i] = 0.00;
	}

	w[0]=sqrt (M_PI_2);
	const bool status = imtqlx ( n, t, bj, w );

	#pragma omp parallel for
	for (i = 0; i < n; ++i)
		w[i] *= w[i];

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _v_double_product_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    V_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) V(i,x)*V(j,x)*sqrt(1+x)/sqrt(1-x) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the polynomial indices.
    0 <= I, J.
    Output, double V_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * dim = data; 
	result = ((dim[XROW] < 0 || dim[YROW] < 0) ? CHEBYSHEV_INVALIDRETURNVALUE : (dim[XROW]==dim[YROW] ? M_PI : 0.00));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _v_polynomial( void * data)
/******************************************************************************/
/*
  Purpose:
    V_POLYNOMIAL evaluates Chebyshev polynomials V(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest polynomial to compute.
    Input, double X[M], the evaluation points.
    Output, double V_POLYNOMIAL[M*(N+1)], the values of the N+1 Chebyshev polynomials.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	if(n<0)
	{
		result = CHEBYSHEV_INVALIDRETURNVALUE; 
		return &result;
	}

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	ityp v[n+1];
	v[0] = 1.00;
	v[1] = 2.00*x - 1.00;

	#pragma omp parallel for
	for (dim_typ i=2; i <= n; ++i)
		v[i] = 2.00 * x * v[i-1] - v[i-2];

	result = v[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _v_polynomial_zeros ( void * data) 
/******************************************************************************/
/*
  Purpose:
    V_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial V(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the polynomial.
    Output, double V_POLYNOMIAL_ZEROS[N], the zeroes.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * z = s_data->a1;
	
	#pragma omp parallel for
	for (dim_typ i=0; i < n; ++i)
		z[i] = cos ((ityp)((n<<1)-(i<<1)-1)*M_PI/(ityp)((n<<1)+1));
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _w_double_product_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    W_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) W(i,x)*W(j,x)*sqrt(1-x)/sqrt(1+x) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the polynomial indices.
    0 <= I, J.
    Output, double W_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * dim = data; 
	result = ((dim[XROW] < 0 || dim[YROW] < 0) ? CHEBYSHEV_INVALIDRETURNVALUE : (dim[XROW]==dim[YROW] ? M_PI : 0.00));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _w_polynomial ( void * data)
/******************************************************************************/
/*
  Purpose:
    W_POLYNOMIAL evaluates Chebyshev polynomials W(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest polynomial to compute.
    Input, double X[M], the evaluation points.
    Output, double W_POLYNOMIAL[M*(N+1)], the values of the N+1 Chebyshev polynomials.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	if(n<0)
	{
		result = CHEBYSHEV_INVALIDRETURNVALUE;
		return &result;
	}

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	ityp v[n+1];
	v[0] = 1.00;
	v[1] = 2.00*x + 1.00;

	#pragma omp parallel for
	for (dim_typ i=2; i <= n; ++i)
		v[i] = 2.00 * x * v[i-1] - v[i-2];

	result = v[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _w_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    W_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial W(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the polynomial.
    Output, double W_POLYNOMIAL_ZEROS[N], the zeroes.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * z = s_data->a1;
	
	#pragma omp parallel for
	for (dim_typ i = 0; i < n; ++i)
		z[i] = cos((ityp)((n-i)<<1)*M_PI/(ityp)((n<<1)+1));
	return NULL;
}

#endif
