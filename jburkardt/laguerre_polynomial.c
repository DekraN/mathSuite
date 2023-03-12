#ifndef __DISABLEDEEP_LAGUERREPOLYNOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l_polynomial( void * data)
/******************************************************************************/
/*
  Purpose:
    L_POLYNOMIAL evaluates the Laguerre polynomials L(n,x).
  Differential equation:
    X * Y'' + (1-X) * Y' + N * Y = 0
  First terms:
      1
     -X    +  1
 (  X^2 -  4 X     +  2 ) / 2
 ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
 (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
 ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
 (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
 ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
      + 52920 X^2 - 35280 X + 5040 ) / 5040
  Recursion:
    L(0,X) = 1,
    L(1,X) = 1-X,
    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
  Orthogonality:
    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,X) * L(M,X) dX
    = 0 if N /= M
    = 1 if N == M
  Special values:
    L(N,0) = 1.
  Relations:
    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )

  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:

    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to compute.
    Input, double X[M], the evaluation points.
    Output, double L_POLYNOMIAL[M*(N+1)], the function values.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	ityp v[n+1];

	if(n < 0)
	{
		result = LAGUERRE_INVALIDRETURNVALUE;
		return &result;
	}

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	v[0] = 1.00;
	v[1] = 1.00 - x;

	for(dim_typ i = 2; i <= n; ++i)
		v[i] = (((ityp)((i<<1)-1) - x)*v[i-1] + (ityp)(-i+1)*v[i-2])/(ityp)(i);

	result = v[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l_polynomial_coefficients( void * data)
/******************************************************************************/
/*
  Purpose:
    L_POLYNOMIAL_COEFFICIENTS: coeffs for Laguerre polynomial L(n,x).
  First terms:
    0: 1
    1: 1  -1
    2: 1  -2  1/2
    3: 1  -3  3/2  1/6
    4: 1  -4  4   -2/3  1/24
    5: 1  -5  5   -5/3  5/24  -1/120
  Recursion:
    L(0,X) = ( 1,  0, 0, ..., 0 )
    L(1,X) = ( 1, -1, 0, ..., 0 )
    L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X) / N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
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

    Output, double L_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients
    of the Laguerre polynomials of degree 0 through N.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * c = s_data->a1; 
	
	if ( n < 0 )
	{
		result = LAGUERRE_DOMAIN_ERROR;
		return &result;
	}

	dim_typ i, j;

	#pragma omp parallel for
	for(i = 0; i <= n; ++i)
		#pragma omp parallel for
		for(j = 0; j <= n; ++j)
			*(c+i*(n+1)+j) = 0.00;

	#pragma omp parallel for
	for(i=0; i<=n; ++i)
		*(c+i*(n+1)) = 1.00;

	if(!n)
	{
		result = LAGUERRE_SUCCESS;
		return &result;
	}

	*(c+(2+n)*(n+1)) = -1.00;

	for(i = 2; i <= n; ++i )
		for(j = 1; j <= n; ++j )
			*(c+i*(n+1)+j) = ((ityp) ((i<<1) - 1 ) * *(c+(i-1)*(n+1)+j) + (ityp) (-i+1) * *(c+(i-2)*(n+1)+j) - *(c+(i-1)*(n+1)+j-1)/(ityp)i);

	result = LAGUERRE_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lm_polynomial ( void * data)
/******************************************************************************/
/*
  Purpose:
    LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
  First terms:
    M = 0
    Lm(0,0,X) =   1
    Lm(1,0,X) =  -X   +  1
    Lm(2,0,X) =   X^2 -  4 X   +  2
    Lm(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
    Lm(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
    Lm(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
    Lm(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
    M = 1
    Lm(0,1,X) =    0
    Lm(1,1,X) =   -1,
    Lm(2,1,X) =    2 X - 4,
    Lm(3,1,X) =   -3 X^2 + 18 X - 18,
    Lm(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
    M = 2
    Lm(0,2,X) =    0
    Lm(1,2,X) =    0,
    Lm(2,2,X) =    2,
    Lm(3,2,X) =   -6 X + 18,
    Lm(4,2,X) =   12 X^2 - 96 X + 144
    M = 3
    Lm(0,3,X) =    0
    Lm(1,3,X) =    0,
    Lm(2,3,X) =    0,
    Lm(3,3,X) =   -6,
    Lm(4,3,X) =   24 X - 96
    M = 4
    Lm(0,4,X) =    0
    Lm(1,4,X) =    0
    Lm(2,4,X) =    0
    Lm(3,4,X) =    0
    Lm(4,4,X) =   24
  Recursion:
    Lm(0,M,X)   = 1
    Lm(1,M,X)   = (M+1-X)
    if 2 <= N:
      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X)
                   +  (1-M-N)    * Lm(N-2,M,X) ) / N
  Special values:
    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal
    to the Laguerre polynomials L(N,X).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
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
    Input, int N, the highest order polynomial to compute.
    Input, int M, the parameter.  M must be nonnegative.
    Input, double X[MM], the evaluation points.
    Output, double LM_POLYNOMIAL[MM*(N+1)], the function values.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register ityp x = s_data->a2;
	
	ityp v[n+1];

	if(n < 0)
	{
		result = LAGUERRE_INVALIDRETURNVALUE;
		return &result;
	}

	dim_typ i, j;

	#pragma omp parallel for
	for(i=0; i<=n; ++i)
		v[i] = 0.00;

	if(!n)
	{
		result = 1.00; 
		return &result;
	}

	v[0] = 1.00;
	v[1] = (ityp) (m+1) - x;

	for (i=2; i <= n; ++i )
		v[i] = (((ityp)(m+(i<<1)-1)-x)*v[i-1]+(ityp)(-m-i+1)*v[i-2])/(ityp)(i);

	result = v[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lm_polynomial_coefficients( void * data)
/******************************************************************************/
/*
  Purpose:
    LM_POLYNOMIAL_COEFFICIENTS: coefficients of Laguerre polynomial Lm(n,m,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author
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
    Input, int M, the parameter.
    Output, double LM_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients
    of the Laguerre polynomials of degree 0 through N.
*/
{
	static bool result = 2;
	
	const _2dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * c = s_data->a2;
	
	if (n < 0)
	{
		result = LAGUERRE_DOMAIN_ERROR;
		return &result;
	}

	dim_typ i, j;

	#pragma omp parallel for
	for(i = 0; i <= n; ++i)
		#pragma omp parallel for
		for(j = 0; j <= n; ++j)
			*(c+i*(n+1)*j) = 0.00;

	*c = 1.00;

	if(!n)
	{
		result = LAGUERRE_SUCCESS;
		return &result;
	}

	*(c+(n+1)) = (ityp)(m+1);
	*(c+(n+1)+1) = -1.00;

	for ( i = 2; i <= n; ++i )
	{
		for ( j = 0; j <= i; ++j )
			*(c+i*(n+1)+j) = ((ityp)(m+(i<<1)-1) * *(c+(i-1)*(n+1)+j)+(ityp)(-m-i+1) * *(c+(i-2)*(n+1)+j))/(ityp)i;
		for ( j = 1; j <= i; ++j )
			*(c+i*(n+1)+j) -= *(c+(i-1)*(n+1)+j-1) / (ityp) i;
	}
	
	result = LAGUERRE_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lf_function( void * data)
/******************************************************************************/
/*
  Purpose:
    LF_FUNCTION evaluates the Laguerre function Lf(n,alpha,x).
  Recursion:
    Lf(0,ALPHA,X) = 1
    Lf(1,ALPHA,X) = 1+ALPHA-X
    Lf(N,ALPHA,X) = (2*N-1+ALPHA-X)/N * Lf(N-1,ALPHA,X)
                      - (N-1+ALPHA)/N * Lf(N-2,ALPHA,X)
  Restrictions:
    -1 < ALPHA
  Special values:
    Lf(N,0,X) = L(N,X).
    Lf(N,ALPHA,X) = LM(N,ALPHA,X) for ALPHA integral.
  Norm:
    Integral ( 0 <= X < +oo ) exp ( - X ) * Lf(N,ALPHA,X)^2 dX
    = Gamma ( N + ALPHA + 1 ) / N!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
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
    Input, int N, the highest order polynomial to compute.
    Input, double ALPHA, the parameter.  -1.0 < ALPHA.
    Input, double X[MM], the evaluation points.
    Output, double LM_POLYNOMIAL[MM*(N+1)], the function values.
*/
{
	static ityp result = MAX_VAL;
	
	const _2itdt * const s_data = data;
	
	const register ityp alpha = s_data->a0;
	const register ityp x = s_data->a1;
	const register dim_typ n = s_data->a2;
	
	ityp v[n+1];

	if (n < 0)
	{
		result = LAGUERRE_INVALIDRETURNVALUE;
		return &result;
	}

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	v[0] = 1.00;
	v[1] = alpha + 1.00 - x;

	for (dim_typ i = 2; i<= n; ++i )
		v[i] = ((alpha+(ityp)((i<<1)-1)-x)*v[i-1]+(-alpha+(ityp)(-i+1))*v[i-2])/(ityp)(i);
	
	result = v[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l_quadrature_rule( void * data)
/******************************************************************************/
/*
  Purpose:
    L_QUADRATURE_RULE: Gauss-Laguerre quadrature based on L(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int N, the order.
    Output, double X[N], the abscissas.
    Output, double W[N], the weights.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
	dim_typ i;
	/*
	Define the zero-th moment.
	*/
	const register ityp zemu = 1.00;
	/*
	Define the Jacobi matrix.
	*/
	ityp bj[n];

	#pragma omp parallel for
	for (i=0; i<n; ++i)
	{
		bj[i] = (ityp)(i+1);
		x[i] = (ityp) ((i<<1)+1);
		w[i] = 0.00;
	}

	w[0] = sqrt(zemu);
	/*
	Diagonalize the Jacobi matrix.
	*/
	const bool status = imtqlx ( n, x, bj, w );

	#pragma omp parallel for
	for (i=0; i<n; ++i)
		w[i] *= w[i];

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lf_quadrature_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    LF_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lf(n,alpha,x);
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int N, the order.
    Input, double ALPHA, the exponent of the X factor.
    ALPHA must be nonnegative.
    Output, double X[N], the abscissas.
    Output, double W[N], the weights.
*/
{
	static bool result = 2;
	
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp alpha = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	dim_typ i;
	ityp i_r8;
	/*
	Define the zero-th moment.
	*/
	const register ityp zemu = r8_gamma(alpha+1.00);
	/*
	Define the Jacobi matrix.
	*/
	ityp bj[n];

	#pragma omp parallel for
	for(i = 0; i < n; ++i)
	{
		i_r8 = (ityp)(i);
		bj[i] = ( i_r8 + 1.00 )*( i_r8 + 1.00 + alpha );
		x[i] = (ityp)(2.00*i_r8+1.00+alpha);
		w[i] = 0.00;
	}

	w[0] = sqrt(zemu);

	/*
	Diagonalize the Jacobi matrix.
	*/
	const bool status = imtqlx(n, x, bj, w);

	#pragma omp parallel for
	for (i = 0; i < n; ++i)
		w[i] *= w[i];

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lm_quadrature_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    LM_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lm(n,m,x);
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int N, the order.
    Input, int M, the parameter.
    0 <= M.
    Output, double X[N], the abscissas.
    Output, double W[N], the weights.
*/
{
	static bool result = 2;
	
	const _2dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	dim_typ i;
	/*
	Define the zero-th moment.
	*/
	const register ityp zemu = r8_factorial ( m );
	/*
	Define the Jacobi matrix.
	*/
	ityp bj[n];


	#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		bj[i] = (ityp)(i+1)*(i+1+m);
		x[i] = (ityp) ((i<<1)+1+m);
		w[i] = 0.00;
	}

	w[0] = sqrt (zemu);
	/*
	Diagonalize the Jacobi matrix.
	*/
	const bool status = imtqlx ( n, x, bj, w );

	#pragma omp parallel for
	for (i = 0; i < n; ++i)
		w[i] *= w[i];

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l_integral( void * data)
/******************************************************************************/
/*
  Purpose:
    L_INTEGRAL evaluates a monomial integral associated with L(n,x).
  Discussion:
    The integral:
      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the exponent.
    0 <= N.
    Output, double L_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = r8_factorial(n);
  	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lf_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    LF_INTEGRAL evaluates a monomial integral associated with Lf(n,alpha,x).
  Discussion:
    The integral:
      integral ( 0 <= x < +oo ) x^n * x^alpha * exp ( -x ) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the exponent.
    0 <= N.
    Input, double ALPHA, the exponent of X in the weight function.
    Output, double LF_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp alpha = s_data->a1;
	
	result = r8_gamma(alpha+(ityp)(n+1));
  	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lm_integral( void * data)
/******************************************************************************/
/*
  Purpose:
    LM_INTEGRAL evaluates a monomial integral associated with Lm(n,m,x).
  Discussion:
    The integral:
      integral ( 0 <= x < +oo ) x^n * x^m * exp ( -x ) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the exponent.
    0 <= N.
    Input, int M, the parameter.
    0 <= M.
    Output, double LM_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ m = a_data[1];
	
	result = r8_factorial(n+m);
  	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    L_POLYNOMIAL_ZEROS: zeros of the Laguerre polynomial L(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int N, the order of the polynomial.
    Output, double L_POLYNOMIAL_ZEROS[N], the zeros.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1; 
	
	dim_typ i;
	ityp w[n];
	/*
	Define the zero-th moment.
	*/
	ityp zemu = 1.00;
	/*
	Define the Jacobi matrix.
	*/
	ityp bj[n];

	#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		bj[i] = (ityp) (i+1);
		x[i] = (ityp) ((i<<1)+1);
		w[i] = 0.00;
	}

	w[0] = sqrt ( zemu );
	/*
	Diagonalize the Jacobi matrix.
	*/
	
	result =imtqlx ( n, x, bj, w );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lf_function_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    LF_FUNCTION_ZEROS returns the zeros of Lf(n,alpha,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int N, the order.
    Input, double ALPHA, the exponent of the X factor.
    ALPHA must be nonnegative.
    Output, double LF_FUNCTION_ZEROS[N], the zeros.
*/
{
	static bool result = 2;
	
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	const register ityp alpha = s_data->a2;
	
	dim_typ i;
	ityp i_r8;
	ityp w[n];
	/*
	Define the zero-th moment.
	*/
	ityp zemu = r8_gamma(alpha+1.0);
	/*
	Define the Jacobi matrix.
	*/
	ityp bj[n];

	#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		i_r8 = ( ityp ) ( i );
		bj[i] = ( i_r8 + 1.00 ) * ( i_r8 + 1.00 + alpha );
		x[i] = ( ityp ) ( 2.0 * i_r8 + 1.0 + alpha );
		w[i] = 0.00;
	}

	w[0] = sqrt (zemu);
	/*
	Diagonalize the Jacobi matrix.
	*/
	const bool status = imtqlx ( n, x, bj, w );

	#pragma omp parallel for
	for (i = 0; i < n; ++i )
		w[i] *= w[i];

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lm_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    LM_POLYNOMIAL_ZEROS returns the zeros for Lm(n,m,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int N, the order.
    Input, int M, the parameter.
    0 <= M.
    Output, double X[N], the zeros.
*/
{
	static bool result = 2;
	
	const _2dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * x = s_data->a2;
	
	dim_typ i;
	ityp w[n];
	/*
	Define the zero-th moment.
	*/
	const register ityp zemu = r8_factorial (m);
	/*
	Define the Jacobi matrix.
	*/
	ityp bj[n];

	#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		bj[i] = (ityp)(i + 1)*(i+1+ m );
		x[i] = ( ityp ) ((i<<1)+1+m);
		w[i] = 0.00;
	}

	w[0] = sqrt (zemu);

	/*
	Diagonalize the Jacobi matrix.
	*/
	result = imtqlx ( n, x, bj, w );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lm_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LM_POLYNOMIAL_VALUES: some values of the Laguerre polynomial Lm(n,m,x).
  Discussion:
    In Mathematica, the function can be evaluated by:
      LaguerreL[n,m,x]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: Q76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the function.
    Output, int *M, the parameter.
    Output, double *X, the point where the function is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _3pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	dim_typ * m = s_data->a2;
	ityp * x = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 20

    static ityp fx_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1500000000000000E+01,
        0.1625000000000000E+01,
        0.1479166666666667E+01,
        0.1148437500000000E+01,
        0.4586666666666667E+00,
        0.2878666666666667E+01,
        0.8098666666666667E+01,
        0.1711866666666667E+02,
        0.1045328776041667E+02,
        0.1329019368489583E+02,
        0.5622453647189670E+02,
        0.7484729341779436E+02,
        0.3238912982762806E+03,
        0.4426100000097533E+03,
        0.1936876572288250E+04
    };

    static dim_typ m_vec[N_MAX] =
    {
        0, 0, 0, 0,
        0, 1, 1, 1,
        1, 0, 1, 2,
        3, 2, 2, 3,
        3, 4, 4, 5
    };

    static dim_typ n_vec[N_MAX] =
    {
        1,  2,  3,  4,
        5,  1,  2,  3,
        4,  3,  3,  3,
        3,  4,  5,  6,
        7,  8,  9, 10
    };

    static ityp x_vec[N_MAX] =
    {
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.50E+00,
        0.50E+00,
        0.50E+00,
        0.50E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
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
        *n_data = *n = *m = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *m = m_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
