#ifndef __DISABLEDEEP_LEGENDREPOLYNOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _p_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    P_INTEGRAL evaluates a monomial integral associated with P(n,x).
  Discussion:
    The integral:
      integral ( -1 <= x < +1 ) x^n dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the exponent.
    0 <= N.
    Output, double P_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ((n%2) == 1 ? 0.00 : 2.00/(ityp)(n+1));
  	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _p_quadrature_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    P_QUADRATURE_RULE: quadrature for Legendre function P(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
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
		bj[i] = sqrt((ityp)((i+1)*(i+1) )/ (ityp)(((i+1)<<2)*(i+1)-1));
		wts[i] = t[i] = 0.00;
	}

	wts[0] = sqrt (2.00);

	const bool status = imtqlx ( nt, t, bj, wts );

	#pragma omp parallel for
	for (i = 0; i < nt; ++i)
		wts[i] = pow ( wts[i], 2);


	result = status;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _p_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    P_POLYNOMIAL_ZEROS: zeros of Legendre function P(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NT, the order of the rule.
    Output, double P_POLYNOMIAL_ZEROS[NT], the zeros.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	
	dim_typ i;
	ityp wts[nt];
	ityp bj[nt];
	#pragma omp parallel for
	for(i = 0; i < nt; ++i)
	{
		bj[i] = (ityp)((i+1)*(i+1))/(ityp)(((i+1)<<2)*(i+1)-1);
		bj[i] = sqrt(bj[i]);
		wts[i] = t[i] = 0.00;
	}

	wts[0] = sqrt (2.00);
	
	result = imtqlx ( nt, t, bj, wts );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _p_polynomial_prime ( void * data)
/******************************************************************************/
/*
  Purpose:
    P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P(n,x).
  Discussion:
    P(0,X) = 1
    P(1,X) = X
    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
    P'(0,X) = 0
    P'(1,X) = 1
    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
    Thanks to Dimitriy Morozov for pointing out a memory leak caused by
    not deleting the work array V before return, 19 March 2013.
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
    Output, double P_POLYNOMIAL_PRIME[M*(N+1)], the values of the derivatives
    of the Legendre polynomials of order 0 through N at the points.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	ityp v[n+1];
	ityp vp[n+1];

	if ( n < 0 )
	{
		result = LEGENDRE_INVALIDRETURNVALUE;
		return &result;
	}

	if ( n < 1 )
	{
		result = 0.00;
		return &result;
	}

	vp[0] = 0.00;
	v[0] = 1.00;

	v[1] = x;
	vp[1] = 1.00;

	for (dim_typ i=2; i <= n; ++i )
	{
		v[i] = ((ityp)((i<<1)-1)*x*v[i-1]-(ityp)(i-1)*v[i-2])/(ityp)(i);
		vp[i] = ( (ityp) ((i<<1)-1)*(v[i-1]+x*vp[i-1])-(ityp) (i-1)*vp[i-2])/(ityp)(i);
	}

	result = vp[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _p_polynomial_coefficients ( void * data)
/******************************************************************************/
/*
  Purpose:
    P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).
  First terms:
     1
     0     1
    -1/2   0      3/2
     0    -3/2    0     5/2
     3/8   0    -30/8   0     35/8
     0    15/8    0   -70/8    0     63/8
    -5/16  0    105/16  0   -315/16   0    231/16
     0   -35/16   0   315/16   0   -693/16   0    429/16
     1.00000
     0.00000  1.00000
    -0.50000  0.00000  1.50000
     0.00000 -1.50000  0.00000  2.5000
     0.37500  0.00000 -3.75000  0.00000  4.37500
     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
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
    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Output, double P_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients of
    the Legendre polynomials of degree 0 through N.
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
	for (i=0; i<=n; ++i)
		#pragma omp parallel for
		for (j=0; j<=n; ++j)
			*(c+i*(n+1)+j) = 0.00;

	*c = 1.00;

	if (!n)
	{
		result = LEGENDRE_SUCCESS;
		return &result;
	}

	*(c+(n+1)+1) = 1.00;

	for (i = 2; i <= n; ++i )
	{
		for ( j = 0; j <= i-2; ++j )
		*(c+i*(n+1)+j) = (ityp)(-i+1)* *(c+(i-2)*(n+1)+j)/(ityp) i;

		for ( j = 1; j <= i; ++j )
			*(c+i*(n+1)+j) = *(c+i*(n+1)+j)+(ityp)(i+i-1)* *(c+(i-1)*(n+1)+j-1)/(ityp)i;
	}

	result = LEGENDRE_SUCCESS;
	return &result;
}



/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _legendre_function_q_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_FUNCTION_Q_VALUES returns values of the Legendre Q function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 February 2003
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input/output, int *N_DATA.
    On input, if N_DATA is 0, the first test data is returned, and N_DATA
    is set to 1.  On each subsequent call, the input value of N_DATA is
    incremented and that test data item is returned, if available.  When
    there is no more test data, N_DATA is set to 0.
    Output, int *N, the order of the function.
    Output, double *X, the point where the function is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 12

    ityp fx_vec[N_MAX] =
    {
        0.00000000, -1.00000000,  0.00000000,
        0.66666667, -0.40634921,  0.00000000,
        0.54930614, -0.72534693, -0.81866327,
        -0.19865477, -0.11616303,  0.29165814
    };
    dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  9, 10,
        0,  1,  2,
        3,  9, 10
    };
    ityp x_vec[N_MAX] =
    {
        0.00,  0.00,  0.00,
        0.00,  0.00,  0.00,
        0.50,  0.50,  0.50,
        0.50,  0.50,  0.50
    };

    if ( N_MAX <= *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data];
        *x = x_vec[*n_data];
        *fx = fx_vec[*n_data];
        ++ *n_data;
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _p_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    P_POLYNOMIAL_VALUES returns values of the Legendre polynomials P(n,x).
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
    Output, int *N, the order of the function.
    Output, double *X, the point where the function is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 22

    static ityp fx_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2500000000000000E+00,
        -0.4062500000000000E+00,
        -0.3359375000000000E+00,
        0.1577148437500000E+00,
        0.3397216796875000E+00,
        0.2427673339843750E-01,
        -0.2799186706542969E+00,
        -0.1524540185928345E+00,
        0.1768244206905365E+00,
        0.2212002165615559E+00,
        0.0000000000000000E+00,
        -0.1475000000000000E+00,
        -0.2800000000000000E+00,
        -0.3825000000000000E+00,
        -0.4400000000000000E+00,
        -0.4375000000000000E+00,
        -0.3600000000000000E+00,
        -0.1925000000000000E+00,
        0.8000000000000000E-01,
        0.4725000000000000E+00,
        0.1000000000000000E+01
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10,  3,
        3,  3,  3,
        3,  3,  3,
        3,  3,  3,
        3
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
        0.00E+00,
        0.10E+00,
        0.20E+00,
        0.30E+00,
        0.40E+00,
        0.50E+00,
        0.60E+00,
        0.70E+00,
        0.80E+00,
        0.90E+00,
        1.00E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pm_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    PM_POLYNOMIAL_VALUES returns values of Legendre polynomials Pm(n,m,x).
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
    Output, int *N, int *M, double *X,
    the arguments of the function.
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
        0.0000000000000000E+00,
        -0.5000000000000000E+00,
        0.0000000000000000E+00,
        0.3750000000000000E+00,
        0.0000000000000000E+00,
        -0.8660254037844386E+00,
        -0.1299038105676658E+01,
        -0.3247595264191645E+00,
        0.1353164693413185E+01,
        -0.2800000000000000E+00,
        0.1175755076535925E+01,
        0.2880000000000000E+01,
        -0.1410906091843111E+02,
        -0.3955078125000000E+01,
        -0.9997558593750000E+01,
        0.8265311444100484E+02,
        0.2024442836815152E+02,
        -0.4237997531890869E+03,
        0.1638320624828339E+04,
        -0.2025687389227225E+05
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
