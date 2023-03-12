#ifndef __DISABLEDEEP_PINKNOISE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cdelay2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    CDELAY2 is a circular buffer implementation of M-fold delay.
  Example:
    Suppose we call CDELAY2 12 times, always with M = 3, and with
    Q having the input value 3 on the first call.  Q will go through
    the following sequence of values over the 12 calls:
    I   M  Qin  Qout
    1   3   3   2
    2   3   2   1
    3   3   1   0
    4   3   0   3
    5   3   3   2
    6   3   2   1
    7   3   1   0
    8   3   0   3
    9   3   3   2
   10   3   2   1
   11   3   1   0
   12   3   0   3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 May 2010
  Author:
    Original C version by Sophocles Orfanidis.
    This C version by John Burkardt.
  Reference:
    Sophocles Orfanidis,
    Introduction to Signal Processing,
    Prentice-Hall, 1995,
    ISBN: 0-13-209172-0,
    LC: TK5102.5.O246.
  Parameters:
    Input, int M, the maximum value that Q can have.
    Input/output, int *Q, a counter which is decremented on every call.
    However, the value "after" 0 is M.
*/
{
	const ipi * const s_data = data;
	const register int m = s_data->a0;
	int * q = s_data->a1;
	
    /*
    Decrement the offset.
    */
    -- *q;
    /*
    Q = - 1 wraps to Q = M.
    */
    wrap2 ( m, q );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _corr ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORR computes the sample correlation of a signal sample.
  Discussion:
    The sample correlation is defined, for 0 <= i < N, as
      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) X(i+j) * X(j)
    The sample correlation is an estimate of the correlation function.
    It is usually the case that the signal X is assumed to
    have zero mean.  Here, we compute the mean and adjust the
    calculation accordingly:
      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i )
  ( X(i+j) - Xbar ) * ( X(j) - Xbar )
    Experience suggests that only the first 5 or 10 percent of
    the lags are statistically reliable, so that one might choose
    M = N / 20 or M = N / 10, for instance.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 June 2010
  Author:
    John Burkardt
  Reference:
    Sophocles Orfanidis,
    Introduction to Signal Processing,
    Prentice-Hall, 1995,
    ISBN: 0-13-209172-0,
    LC: TK5102.5.O246.
  Parameters:
    Input, int N, the number of equally spaced signal
    samples.
    Input, double X[N], the signal samples.
    Input, int M, the maximum lag to consider.
    0 <= M < N.
    Output, double CORR[M+1], the sample correlations.
*/
{
	const _2dtpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * x = s_data->a2;
	
    dim_typ i, j;
    ityp *r;
    ityp xbar;

    r = ( ityp * ) malloc ( ( m + 1 ) * sizeof ( ityp ) );

    for ( i = 0; i <= m; ++i)
        r[i] = 0.00;

    xbar = 0.00;
    for ( j = 0; j < n; ++j )
        xbar += x[j];

    xbar /= ( ityp ) ( n );

    for ( i = 0; i <= m; ++i )
        for ( j = 0; j < n - i; ++j )
            r[i] += ( x[i+j] - xbar ) * ( x[j] - xbar );

    for ( i = 0; i <= m; ++i )
        r[i] /= ( ityp ) ( n );

    return r;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cross_corr ( void * data)
/******************************************************************************/
/*
  Purpose:
    CROSS_CORR computes the sample cross correlation between two signal samples.
  Discussion:
    The sample cross correlation is defined, for 0 <= i < N, as
      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i ) X(i+j) * Y(j)
    The sample cross correlation is an estimate of the cross
    correlation function.
    It is usually the case that the signals X and Y are assumed to
    have zero mean.  Here, we compute the means and adjust the
    calculation accordingly:
      R(i) = 1/N * sum ( 0 <= j <= N - 1 - i )
  ( X(i+j) - Xbar ) * ( Y(j) - Ybar )
    Experience suggests that only the first 5 or 10 percent of
    the lags are statistically reliable, so that one might choose
    M = N / 20 or M = N / 10, for instance.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 June 2010
  Author:
    John Burkardt
  Reference:
    Sophocles Orfanidis,
    Introduction to Signal Processing,
    Prentice-Hall, 1995,
    ISBN: 0-13-209172-0,
    LC: TK5102.5.O246.
  Parameters:
    Input, int N, the number of equally spaced signal
    samples.
    Input, double X[N], Y[N], the signal samples.
    Input, int M, the maximum lag to consider.
    0 <= M < N.
    Output, double CROSS_CORR[M+1], the sample correlations.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * x = s_data->a2;
	ityp * y = s_data->a3;
	
	
    dim_typ i, j;
    ityp *r;
    ityp xbar;
    ityp ybar;

    r = ( ityp * ) malloc ( ( m + 1 ) * sizeof ( ityp ) );

    for ( i = 0; i <= m; ++i )
        r[i] = 0.00;

    xbar = 0.00;
    for ( j = 0; j < n; ++j )
        xbar += x[j];


    xbar /=( ityp ) ( n );

    ybar = 0.00;
    for ( j = 0; j < n; ++j )
        ybar += x[j];

    ybar /= ( ityp ) ( n );

    for ( i = 0; i <= m; i++ )
        for ( j = 0; j < n - i; ++j )
            r[i] += ( x[i+j] - xbar ) * ( y[j] - ybar );

    for ( i = 0; i <= m; ++i )
        r[i] /= ( ityp ) ( n );

    return r;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ran1f ( void * data)
/******************************************************************************/
/*
  Purpose:
    RAN1F is a 1/F random number generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 May 2010
  Author:
    Original C version by Sophocles Orfanidis.
    This C version by John Burkardt.
  Reference:
    Sophocles Orfanidis,
    Introduction to Signal Processing,
    Prentice-Hall, 1995,
    ISBN: 0-13-209172-0,
    LC: TK5102.5.O246.
  Parameters:
    Input, int B, the number of signals to combine.
    For this algorithm, B cannot be more than 31!
    Input/output, double U[B], the signals to combine.  It is expected
    that each of the initial values of U will be drawn from a distribution
    with zero mean.
    Input/output, int Q[B], a set of counters that determine when each
    entry of U is to be updated.
    Output, double RAN1F, the value.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitpi * const s_data = data;
	const register dim_typ b = s_data->a0;
	ityp * u = s_data->a1;
	int * q = s_data->a2;
	
    dim_typ i, j;
    ityp y;

    if ( 31 < b )
    {
    	result = MAX_VAL;
        return &result;
    }

    y = 0.00;

    j = 1;
    for ( i = 0; i < b; ++i)
    {
        y += ranh ( j, u+i, q+i );
        j *= 2;
    }

	result = 0<b ? y/b:y;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ranh ( void * data)
/******************************************************************************/
/*
  Purpose:
    RANH is a hold random number generator of period D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2010
  Author:
    Original C version by Sophocles Orfanidis.
    This C version by John Burkardt.
  Reference:
    Sophocles Orfanidis,
    Introduction to Signal Processing,
    Prentice-Hall, 1995,
    ISBN: 0-13-209172-0,
    LC: TK5102.5.O246.
  Parameters:
    Input, int D, the hold period.  D must be at least 1.
    Input/output, double *U, a value to be held until Q has decremented
    to 0, when Q will be reset to D, and U will be randomly reset.
    Input/output, int *Q, a counter which is decremented by 1 on each call
    until reaching 0.
    Output, double RANH, the input value of U.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitpi * const s_data = data;
	const register dim_typ d = s_data->a0;
	ityp * u = s_data->a1;
	int * q = s_data->a2;
	
    ityp y;

    if ( d < 1 )
    {
    	result = MAX_VAL;
		return &result;
	}
    /*
    Hold this sample for D calls.
    */
    y = *u;
    /*
    Decrement Q and wrap mod D.
    */
    cdelay2 ( d - 1, q );
    /*
    Every D calls, get a new U with zero mean.
    */
    if ( *q == 0 )
        *u = 2.00 * ( ityp ) rand ( ) / ( ityp ) ( RAND_MAX ) - 1.00;
        
    result = y;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wrap2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    WRAP2 is a circular wrap of the pointer offset Q.
  Discussion:
    Input values of Q between 0 and M are "legal".
    Values of Q below 0 are incremented by M + 1 until they are legal.
    Values of Q above M are decremented by M + 1 until they become legal.
    The legal value is the output value of the function.
  Example:
    M  Qin  Qout
    3  -5   3
    3  -4   0
    3  -3   1
    3  -2   2
    3  -1   3
    3   0   0
    3   1   1
    3   2   2
    3   3   3
    3   4   0
    3   5   1
    3   6   2
    3   7   3
    3   8   0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 May 2010
  Author:
    Original C version by Sophocles Orfanidis.
    This C version by John Burkardt.
  Reference:
    Sophocles Orfanidis,
    Introduction to Signal Processing,
    Prentice-Hall, 1995,
    ISBN: 0-13-209172-0,
    LC: TK5102.5.O246.
  Parameters:
    Input, int M, the maximum acceptable value for outputs.
    M must be at least 0.
    Input/output, int *Q, the value to be wrapped.
*/
{
	const ipi * const s_data = data;
	const register int m = s_data->a0;
	int * q = s_data->a1;
	
    if ( m < 0 )
        return NULL;
    /*
    When Q = M + 1, it wraps to Q = 0.
    */
    while ( m < *q )
        *q -= m - 1;
    /*
    When Q = - 1, it wraps to Q = M.
    */
    while ( *q < 0 )
        *q += m + 1;

    return NULL;
}

#endif
