#ifndef __DISABLEDEEP_HERMITEPOLYNOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _h_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    H_INTEGRAL evaluates the integral of H(i,x).
  Discussion:
    H(i,x) is the physicist's Hermite polynomial of degree I.
    The integral computed is:
      integral ( -oo < x < +oo ) H(i,x) exp(-x^2) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 March 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the integral.
    0 <= N.
    Output, double H_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ((n%2) == 1 ? 0.00 : r8_factorial2 (n-1)*sqrt (M_PI)/pow(2.00, n/2));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _h_polynomial_coefficients ( void * data)
/******************************************************************************/
/*
  Purpose:
    H_POLYNOMIAL_COEFFICIENTS: coefficients of H(i,x).
  Discussion:
    H(i,x) is the physicist's Hermite polynomial of degree I.
  First terms:
    N/K     0     1      2      3       4     5      6    7      8    9   10
     0      1
     1      0     2
     2     -2     0      4
     3      0   -12      0      8
     4     12     0    -48      0      16
     5      0   120      0   -160       0    32
     6   -120     0    720      0    -480     0     64
     7      0 -1680      0   3360       0 -1344      0   128
     8   1680     0 -13440      0   13440     0  -3584     0    256
     9      0 30240      0 -80640       0 48384      0 -9216      0 512
    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 March 2012
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
    Output, double HERMITE_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the
    coefficients of the Hermite polynomials of orders 0 through N.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * c = s_data->a1;

	if(n<0)
	{
		result = HERMITE_DOMAIN_ERROR;
		return &result;
	}

	dim_typ i, j;

	#pragma omp parallel for
	for ( i = 0; i <= n; ++i )
		#pragma omp parallel for
		for ( j = 0; j <= n; ++j )
			*(c+i*(n+1)+j) = 0.00;

	*c = 1.00;

	if(!n)
	{
		result = HERMITE_SUCCESS;
		return &result;
	}

	*(c+(n+1)+1) = 2.00;

	for ( i = 2; i <= n; ++i )
	{
		*(c+i*(n+1)) = - 2.00 * (ityp) (i-1) * *(c+(i-2)*(n+1));
		for ( j = 1; j <= i - 2; ++j )
			*(c+i*(n+1)+j) = 2.00 * *(c+(i-1)*(n+1)+(j-1))-2.00*(ityp)(i-1) * *(c+(i-2)*(n+1)+j);
		*(c + i*(n+1) + (i-1)) = 2.00 * *(c + (i-1)*(n+1)+(i-2));
		*(c + i*(n+1) + i) = 2.00 * *(c + (i-1)*(n+1) + (i-1));
	}

	result = HERMITE_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _h_polynomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    H_POLYNOMIAL_VALUE evaluates H(i,x).
  Discussion:
    H(i,x) is the physicist's Hermite polynomial of degree I.
  Differential equation:
    Y'' - 2 X Y' + 2 N Y = 0
  First terms:
      1
      2 X
      4 X^2     -  2
      8 X^3     - 12 X
     16 X^4     - 48 X^2     + 12
     32 X^5    - 160 X^3    + 120 X
     64 X^6    - 480 X^4    + 720 X^2    - 120
    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
  Recursion:
    H(0,X) = 1,
    H(1,X) = 2*X,
    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
  Norm:
    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
    = sqrt ( M_PI ) * 2^N * N!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
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
    Note that polynomials 0 through N will be computed.
    Input, double X[M], the evaluation points.
    Output, double H_POLYNOMIAL_VALUE[M*(N+1)], the values of the first N+1
    Hermite polynomials at the evaluation points.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	if(n<0)
	{
		result = HERMITE_INVALIDRETURNVALUE;
		return &result;
	}

	ityp p[n+1];

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	p[0] = 1.00;
	p[1] = 2.00 * x;

	for (dim_typ i=2; i <= n; ++i)
		p[i] = 2.00*x*p[i-1]-2.00*(ityp)(i-1)*p[i-2];

	result = p[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _h_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    H_POLYNOMIAL_ZEROS: zeros of H(i,x).
  Discussion:
    H(i,x) is the physicist's Hermite polynomial of degree I.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int NT, the degree of the polynomial.
    Output, double H_POLYNOMIAL_ZEROS[NT], the zeros of the polynomial.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * z = s_data->a1;
	
	ityp bj[nt];
	ityp wts[nt];

	#pragma omp parallel for
	for(dim_typ i = 0; i < nt; ++i )
	{
		bj[i] = sqrt((ityp)(i+1)/2.00);
		wts[i] = z[i] = 0.00;
	}

	wts[0] = sqrt(sqrt(M_PI));

	result = imtqlx ( nt, z, bj, wts );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _h_quadrature_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    H_QUADRATURE_RULE: quadrature for H(i,x).
  Discussion:
    H(i,x) is the physicist's Hermite polynomial of degree I.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int NT, the order of the rule.
    Output, double T[NT], WTS[NT], the points and weights
    of the rule.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	ityp * wts = s_data->a2;
	
	dim_typ i;
	ityp bj[nt];

	#pragma omp parallel for
	for(i = 0; i < nt; ++i)
	{
		bj[i] = sqrt ((ityp)(i+1)/2.00);
		wts[i] = t[i] = 0.00;
	}

	wts[0] = sqrt(sqrt(M_PI));

	const bool status = imtqlx ( nt, t, bj, wts );

	#pragma omp parallel for
	for(i = 0; i < nt; ++i)
		wts[i] *= wts[i];

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _he_double_product_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    HE_DOUBLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*e^(-x^2/2).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 March 2012
  Author:
    John Burkardt
  Reference:
    Dongbin Xiu,
    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
    Princeton, 2010,
    ISBN13: 978-0-691-14212-8,
    LC: QA274.23.X58.
  Parameters:
    Input, int I, J, the polynomial indices.
    Output, double HE_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * dim = data;
	
	result = (dim[XROW] == dim[YROW] ? r8_factorial(dim[XROW]) : 0.00);
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _he_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    HE_INTEGRAL evaluates the integral of He(i,x).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 March 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the integral.
    0 <= N.
    Output, double HE_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ((n%2) == 1 ? 0.00 : r8_factorial2(n-1)*sqrt(M_2TPI));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _he_polynomial_coefficients ( void * data)
/******************************************************************************/
/*
  Purpose:
    HE_POLYNOMIAL_COEFFICIENTS: coefficients of He(i,x).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
  First terms:
    N/K     0     1      2      3       4     5      6    7      8    9   10
     0      1
     1      0     1
     2     -1     0      1
     3      0    -3      0      1
     4      3     0     -6      0       1
     5      0    15      0    -10       0     1
     6    -15     0     45      0     -15     0      1
     7      0  -105      0    105       0   -21      0     1
     8    105     0   -420      0     210     0    -28     0      1
     9      0   945      0  -1260       0   378      0   -36      0   1
    10   -945     0   4725      0   -3150     0    630     0    -45   0    1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 March 2012
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
    Output, double HE_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients
    of the Hermite polynomials of orders 0 through N.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * c = s_data->a1;
	
	if(n<0)
	{
		result = HERMITE_DOMAIN_ERROR;
		return &result;
	}

	dim_typ i, j;

	#pragma omp parallel for
	for ( i = 0; i <= n; ++i )
		#pragma omp parallel for
		for ( j = 0; j <= n; ++j )
			*(c+i*(n+1)+j) = 0.00;

	*c = 1.00;

	if(!n)
	{
		result = HERMITE_SUCCESS;
		return &result;
	}

	*(c+(n+1)+1) = 1.00;

	for ( i = 2; i <= n; ++i)
	{
		*(c+i*(n+1)) = - ( ityp ) ( i - 1 ) * *(c+(i-2)*(n+1));

		for ( j = 1; j <= i - 2; ++j)
			*(c+i*(n+1)+j) = *(c+(i-1)*(n+1)+(j-1))-(ityp)(i-1) * *(c+(i-2)*(n+1)+j);

		*(c+i*(n+1)+(i-1)) = *(c+(i-1)*(n+1)+(i-2));
		*(c+i*(n+1)+i) = *(c+(i-1)*(n+1)+(i-1));
	}

	result = HERMITE_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _he_polynomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    HE_POLYNOMIAL_VALUE evaluates He(i,x).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
  Differential equation:
 ( exp ( - 0.5 * x^2 ) * He(n,x)' )' + n * exp ( - 0.5 * x^2 ) * He(n,x) = 0
  First terms:
   1
   X
   X^2  -  1
   X^3  -  3 X
   X^4  -  6 X^2 +   3
   X^5  - 10 X^3 +  15 X
   X^6  - 15 X^4 +  45 X^2 -   15
   X^7  - 21 X^5 + 105 X^3 -  105 X
   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
  Recursion:
    He(0,X) = 1,
    He(1,X) = X,
    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
  Orthogonality:
    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX
    = sqrt ( 2 * M_PI ) * N// * delta ( N, M )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2012
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
    NIST Handbook of Mathematical Functions,
    Cambridge University Press, 2010,
    ISBN: 978-0521192255,
    LC: QA331.N57.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to compute.
    Note that polynomials 0 through N will be computed.
    Input, double X[M], the evaluation points.
    Output, double HE_POLYNOMIAL_VALUE[M*(N+1)], the values of the
    probabilist's Hermite polynomials of index 0 through N.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	if(n<0)
	{
		result = HERMITE_INVALIDRETURNVALUE;
		return &result;
	}

	ityp p[n+1];

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	p[0] = 1.00;
	p[1] = x;

	dim_typ i;

	for (i=2; i <= n; ++i)
		p[i] = x*p[i-1]-(ityp)(i-1)*p[i-2];

	result = p[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _he_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    HE_POLYNOMIAL_ZEROS: zeros of He(i,x).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int NT, the degree of the polynomial.
    Output, double Z[NT], the zeros of the polynomial.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * z = s_data->a1;
	
	dim_typ i;
	ityp bj[nt];
	ityp wts[nt];

	#pragma omp parallel for
	for(dim_typ i = 0; i < nt; ++i)
	{
		bj[i] = sqrt((ityp)(i+1)/2.00);
		wts[i] = z[i] = 0.00;
	}

	wts[0] = sqrt(sqrt(M_PI));

	const bool status = imtqlx ( nt, z, bj, wts );

	#pragma omp parallel for
	for ( i = 0; i < nt; ++i)
		z[i] *= sqrt(2.00);

	result = status; 
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _he_quadrature_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    HE_QUADRATURE_RULE: quadrature for He(i,x).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int NT, the order of the rule.
    Output, double T[NT], WTS[NT], the points and weights
    of the rule.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	ityp * wts = s_data->a2;
	
	dim_typ i;
	ityp bj[nt];

	#pragma omp parallel for
	for (i = 0; i < nt; ++i)
	{
		bj[i] = sqrt ((ityp)(i+1/2.00));
		wts[i] = t[i] = 0.00;
	}

	wts[0] = sqrt(sqrt(M_PI));

	const bool status = imtqlx ( nt, t, bj, wts );

	#pragma omp parallel for
	for ( i = 0; i < nt; ++i )
	{
		t[i] *= sqrt(2.00);
		wts[i] *= wts[i] * sqrt(2.00);
	}

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _he_triple_product_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    HE_TRIPLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 March 2012
  Author:
    John Burkardt
  Reference:
    Dongbin Xiu,
    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
    Princeton, 2010,
    ISBN13: 978-0-691-14212-8,
    LC: QA274.23.X58.
  Parameters:
    Input, int I, J, K, the polynomial indices.
    Output, double HE_TRIPLE_PRODUCT_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * dim = data;
	
	const register dim_typ sm = dim[XROW]+dim[YROW]+dim[ZROW];
	const register dim_typ s = sm/2;
	
	result = ((s < dim[XROW] || s < dim[YROW] || s < dim[ZROW] || sm%2) ? 0.00 : r8_factorial(dim[XROW])/r8_factorial(s-dim[XROW])*r8_factorial(dim[YROW])/r8_factorial(s-dim[YROW])*r8_factorial ( dim[ZROW] ) / r8_factorial ( s - dim[ZROW]));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hen_polynomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    HEN_POLYNOMIAL_VALUE: evaluates Hen(i,x).
  Discussion:
    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
    These polynomials satisfy the orthonormality condition:
      Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * Hen(M,X) Hen(N,X) dX
      = delta ( N, M )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2012
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
    NIST Handbook of Mathematical Functions,
    Cambridge University Press, 2010,
    ISBN: 978-0521192255,
    LC: QA331.N57.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to compute.
    Note that polynomials 0 through N will be computed.
    Input, double X[M], the evaluation points.
    Output, double HEN_POLYNOMIAL_VALUE[M*(N+1)], the values of the
    polynomials of index 0 through N.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;

	if (n < 0)
	{
		result = HERMITE_INVALIDRETURNVALUE;
		return &result;
	}

	ityp p[n+1];

	if(!n)
	{
		result = 1.00;
		return &result;
	}
	
	p[0] = 1.00;
	p[1] = x;

	dim_typ i;

	for (i=2; i <= n; ++i)
		p[i] = x*p[i-1]-(ityp)(i-1)*p[i-2];

	/*
	Normalize.
	*/
	ityp fact = 1.00;
	for ( i = 0; i <= n; ++i)
	{
		p[i] /= sqrt(fact*sqrt(M_2TPI));
		fact *= (ityp)(i+1);
	}
	
	result = p[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hf_function_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    HF_FUNCTION_VALUE evaluates Hf(i,x).
  Discussion:
    Hf(i,x) represents the Hermite function of "degree" I.
    The Hermite function of degree n is related to the physicist's
    Hermite polynomial H(n,x):
      Hf(n,x) = H(n,x) * exp ( - 0.5 * x^2 ) / sqrt ( 2^n n// sqrt ( M_PI ) )
    The Hermite functions are orthonormal:
      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2012
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
    NIST Handbook of Mathematical Functions,
    Cambridge University Press, 2010,
    ISBN: 978-0521192255,
    LC: QA331.N57.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to compute.
    Note that polynomials 0 through N will be computed.
    Input, double X[M], the point at which the polynomials are
    to be evaluated.
    Output, double HF_FUNCTION_VALUE[M*(N+1)], the values of the Hermite
    functions of index 0 through N at the evaluation points.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	ityp f[n+1];

	f[0] = exp(-0.5*x*x)/sqrt(sqrt(M_PI));

	if(!n)
	{
		result = f[0];
		return &result;
	}
	
	f[1] = 2.00*exp(-0.5*x*x)*x/sqrt(2.00*sqrt(M_PI));


	for (dim_typ i = 2; i <= n; ++i)
		f[i] = (sqrt(2.00)*x*f[i-1]-sqrt((ityp)(i-1))*f[i-2])/sqrt((ityp)(i));

	result = f[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hf_quadrature_rule( void * data)
/******************************************************************************/
/*
  Purpose:
    HF_QUADRATURE_RULE: quadrature for Hf(i,x).
  Discussion:
    Hf(i,x) represents the Hermite function of "degree" I.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int NT, the order of the rule.
    Output, double T[NT], WTS[NT], the points and weights
    of the rule.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	ityp * wts = s_data->a2;

	dim_typ i;
	ityp bj[nt];

	#pragma omp parallel for
	for ( i = 0; i < nt; ++i)
	{
		bj[i] = sqrt((ityp)(i+1)/2.00);
		wts[i] = t[i] = 0.00;
	}

	wts[0] = sqrt(sqrt(M_PI));

	const bool status = imtqlx ( nt, t, bj, wts );

	#pragma omp parallel for
	for (i = 0; i < nt; ++i)
		wts[i] *= wts[i]*exp(t[i]*t[i]);

	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hn_polynomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    HN_POLYNOMIAL_VALUE evaluates Hn(i,x).
  Discussion:
    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.
    These polynomials satisfy the orthonormality condition:
      Integral ( -oo < X < +oo )
        exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2012
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
    Note that polynomials 0 through N will be computed.
    Input, double X[M], the evaluation points.
    Output, double HN_POLYNOMIAL_VALUE[M*(N+1)], the values of the first N+1
    Hermite polynomials at the evaluation points.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	dim_typ i;
	if (n < 0)
	{
		result = HERMITE_INVALIDRETURNVALUE;
		return &result;
	}

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	ityp p[n+1];

	p[0] = 1.00;
	p[1] = 2.00 * x;

	for (i=2; i <= n; ++i)
		p[i] = 2.00*x*p[i-1]-2.00*(ityp)(i-1)*p[i-2];

	/*
	Normalize.
	*/
	ityp fact, two_power;
	fact = two_power = 1.00;
	for ( i = 0; i <= n; ++i )
	{
		p[i] /= sqrt (fact*two_power*sqrt(M_PI));
		fact *= (ityp)(i+1);
		two_power *= 2.00;
	}
	
	result = p[n];
	return &result;
}

#endif
