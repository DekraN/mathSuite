#ifndef __DISABLEDEEP_ASA121

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _trigamma ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIGAMMA calculates trigamma(x) = d^2 log(gamma(x)) / dx^2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 November 2010
  Author:
    Original FORTRAN77 version by BE Schneider.
    C version by John Burkardt.
  Reference:
    BE Schneider,
    Algorithm AS 121:
    Trigamma Function,
    Applied Statistics,
    Volume 27, Number 1, pages 97-99, 1978.
  Parameters:
    Input, double X, the argument of the trigamma function.
    0 < X.
    Output, int *IFAULT, error flag.
    0, no error.
    1, X <= 0.
    Output, double TRIGAMMA, the value of the trigamma function at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	ityp value;
	ityp y;
	ityp z;
	/*
	Check the input.
	*/
	if ( x <= 0.00 )
	{
		result = TRIGAMMA_INVALIDRETURNVALUE;
		return &result;
	}

	z = x;
	/*
	Use small value approximation if X <= A.
	*/
	if ( x <= 0.0001 )
	{
		result = 1.00 / x / x;
		return &result;
	}
	/*
	Increase argument to ( X + I ) >= B.
	*/
	value = 0.00;

	while ( z < 5.0 )
	{
		value += 1.00 / z / z;
		z += 1.00;
	}
	/*
	Apply asymptotic formula if argument is B or greater.
	*/
	y = 1.00 / z / z;
	result = value + 0.5 * y + ( 1.0 	+ y * ( 0.1666666667  + y * ( -0.03333333333  + y * ( 0.02380952381  + y *   -0.03333333333 )))) / z;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _trigamma_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIGAMMA_VALUES returns some values of the TriGamma function.
  Discussion:
    In Mathematica, the function can be evaluated by:
      PolyGamma[1,x]
    TriGamma(X) = d^2 ln ( Gamma ( X ) ) / d X^2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 September 2004
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
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 11

    ityp fx_vec[N_MAX] =
    {
        0.1644934066848226E+01,
        0.1433299150792759E+01,
        0.1267377205423779E+01,
        0.1134253434996619E+01,
        0.1025356590529597E+01,
        0.9348022005446793E+00,
        0.8584318931245799E+00,
        0.7932328301639984E+00,
        0.7369741375017002E+00,
        0.6879720582426356E+00,
        0.6449340668482264E+00
    };

    ityp x_vec[N_MAX] =
    {
        1.0E+00,
        1.1E+00,
        1.2E+00,
        1.3E+00,
        1.4E+00,
        1.5E+00,
        1.6E+00,
        1.7E+00,
        1.8E+00,
        1.9E+00,
        2.0E+00
    };


    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
