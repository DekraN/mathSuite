#ifndef __DISABLEDEEP_ASA226

#include "../dutils.h"

//****************************************************************************80
__MATHSUITE __JBURKARDT  void *   _betanc ( void * data)
//****************************************************************************80
//  Purpose:
//    BETANC computes the tail of the noncentral Beta distribution.
//  Discussion:
//    This routine returns the cumulative probability of X for the non-central
//    Beta distribution with parameters A, B and non-centrality LAMBDA.
//    Note that if LAMBDA = 0, the standard Beta distribution is defined.
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    24 January 2008
//  Author:
//    Original FORTRAN77 version by Russell Lenth.
//    C++ version by John Burkardt.
//  Reference:
//    Russell Lenth,
//    Algorithm AS 226:
//    Computing Noncentral Beta Probabilities,
//    Applied Statistics,
//    Volume 36, Number 2, 1987, pages 241-244.
//    H Frick,
//    Algorithm AS R84:
//    A Remark on Algorithm AS 226:
//    Computing Noncentral Beta Probabilities,
//    Applied Statistics,
//    Volume 39, Number 2, 1990, pages 311-312.
//  Parameters:
//    Input, double X, the value defining the cumulative
//    probability lower tail.  Normally, 0 <= X <= 1, but any value
//    is allowed.
//    Input, double A, B, the parameters of the distribution.
//    0 < A, 0 < B.
//    Input, double LAMBDA, the noncentrality parameter
//    of the distribution.  0 <= LAMBDA.  The program can produce reasonably
//    accurate results for values of LAMBDA up to about 100.
//    Output, int *IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//    Output, double BETANC, the cumulative probability
//    of X.
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	const register ityp lambda = a_data[3];
	
	ityp a0;
	ityp ax;
	ityp beta;
	ityp c;
	ityp errbd;
	ityp gx;
	ityp q;
	ityp sumq;
	ityp temp;
	ityp value;
	ityp x0;
	ityp xj;

	if ( lambda < 0.0 || a <= 0.0 || b <= 0.0 || x <= 0.00 || 1.00 <= x)
	{
		result = BETANC_INVALIDRETURNVALUE;
		return &result;
	}

	c = 0.50 * lambda;
	//  Initialize the series.
	beta = alngam ( a ) + alngam ( b ) - alngam ( a + b );
	temp = betain ( x, a, b, beta );
	gx = exp ( a * log ( x ) + b * log ( 1.00 - x ) - beta - log ( a ) );
	q = exp ( - c );
	xj = 0.00;
	ax = q * temp;
	sumq = 1.00 - q;
	value = ax;

	//  Recur over subsequent terms until convergence is achieved.

	for ( ; ; )
	{
		xj += 1.00;
		temp -= gx;
		gx = x * ( a + b + xj - 1.0 ) * gx / ( a + xj );
		q *= c / xj;
		sumq -= q;
		ax = temp * q;
		value += ax;
		//  Check for convergence and act accordingly.
		errbd = abs ( ( temp - gx ) * sumq );

		if ( errbd <= 1.0E-07 )
			break;

		if (  150 < ( dim_typ ) xj )
		break;
	}
	
	result = value;
	return &result;
}

//****************************************************************************80
__MATHSUITE __JBURKARDT  void * _beta_noncentral_cdf_values ( void * data)
//****************************************************************************80
//  Purpose:
//    BETA_NONCENTRAL_CDF_VALUES returns some values of the noncentral Beta CDF.
//  Discussion:
//    The values presented here are taken from the reference, where they
//    were given to a limited number of decimal places.
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    24 January 2008
//  Author:
//    John Burkardt
//  Reference:
//    R Chattamvelli, R Shanmugam,
//    Algorithm AS 310:
//    Computing the Non-central Beta Distribution Function,
//    Applied Statistics,
//    Volume 46, Number 1, 1997, pages 146-156.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA by 1,
//    and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, double *A, *B, the shape parameters.
//    Output, double *LAMBDA, the noncentrality parameter.
//    Output, double *X, the argument of the function.
//    Output, double *FX, the value of the function.
{
	const _4pitpdtpit * const s_data = data;
	
	ityp * a = s_data->a0;
	ityp * b = s_data->a1;
	ityp * lambda = s_data->a2;
	ityp * x = s_data->a3;
	dim_typ * n_data = s_data->a4;
	ityp * fx = s_data->a5;
	
    # define N_MAX 25

    ityp a_vec[N_MAX] =
    {
        5.00,
        5.00,
        5.00,
        10.00,
        10.00,
        10.00,
        20.00,
        20.00,
        20.00,
        10.00,
        10.00,
        15.00,
        20.00,
        20.00,
        20.00,
        30.00,
        30.0,
        10.00,
        10.00,
        10.00,
        15.00,
        10.00,
        12.00,
        30.00,
        35.00
    };

    ityp b_vec[N_MAX] =
    {
        5.00,
        5.00,
        5.00,
        10.00,
        10.00,
        10.00,
        20.00,
        20.00,
        20.00,
        20.0,
        10.0,
        5.00,
        10.00,
        30.00,
        50.00,
        20.00,
        40.00,
        5.00,
        10.00,
        30.00,
        20.00,
        5.00,
        17.00,
        30.00,
        30.00
    };
    ityp fx_vec[N_MAX] =
    {
        0.4563021,
        0.1041337,
        0.6022353,
        0.9187770,
        0.6008106,
        0.0902850,
        0.9998655,
        0.9925997,
        0.9641112,
        0.9376626573,
        0.7306817858,
        0.1604256918,
        0.1867485313,
        0.6559386874,
        0.9796881486,
        0.1162386423,
        0.9930430054,
        0.0506899273,
        0.1030959706,
        0.9978417832,
        0.2555552369,
        0.0668307064,
        0.0113601067,
        0.7813366615,
        0.8867126477
    };
    ityp lambda_vec[N_MAX] =
    {
        54.00,
        140.00,
        170.00,
        54.00,
        140.00,
        250.00,
        54.00,
        140.00,
        250.00,
        150.00,
        120.00,
        80.00,
        110.00,
        65.00,
        130.00,
        80.00,
        130.00,
        20.00,
        54.00,
        80.00,
        120.00,
        55.00,
        64.00,
        140.00,
        20.00
    };
    ityp x_vec[N_MAX] =
    {
        0.8640,
        0.9000,
        0.9560,
        0.8686,
        0.9000,
        0.9000,
        0.8787,
        0.9000,
        0.9220,
        0.868,
        0.900,
        0.880,
        0.850,
        0.660,
        0.720,
        0.720,
        0.800,
        0.644,
        0.700,
        0.780,
        0.760,
        0.795,
        0.560,
        0.800,
        0.670
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *a = *b = *lambda = *x = *fx = 0.00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *lambda = lambda_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
