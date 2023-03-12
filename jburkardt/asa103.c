#ifndef __DISABLEDEEP_ASA103

#include "../dutils.h"


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _digama ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 June 2013
  Author:
    Original FORTRAN77 version by Jose Bernardo.
    C version by John Burkardt.
  Reference:
    Jose Bernardo,
    Algorithm AS 103:
    Psi ( Digamma ) Function,
    Applied Statistics,
    Volume 25, Number 3, 1976, pages 315-317.
  Parameters:
    Input, double X, the argument of the digamma function.
    0 < X.
    Output, int *IFAULT, error flag.
    0, no error.
    1, X <= 0.
    Output, double DIGAMA, the value of the digamma function at X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	ityp r;
	ityp value;
	ityp x2;
	/*
	Check the input.
	*/
	if ( x <= 0.0 )
	{
		result = 0.00;
		return &result;
	}
	/*
	Initialize.
	*/
	x2 = x;
	value = 0.00;
	/*
	Use approximation for small argument.
	*/
	if ( x2 <= 0.00001 )
	{
		result = - 0.57721566490153286060 - 1.00 / x2;
		return &result;
	}
	/*
	Reduce to DIGAMA(X + N).
	*/
	while ( x2 < 8.50 )
	{
		value -= 1.00/x2;
		x2 += 1.00;
	}
	/*
	Use Stirling's (actually de Moivre's) expansion.
	*/
	r = 1.00 / x2;
	value += log ( x2 ) - 0.50 * r;
	r *= r;
	value -= r * ( 1.00 / 12.00 - r * ( 1.00 / 120.00 - r *   1.00 / 252.00 ) );
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _psi_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    PSI_VALUES returns some values of the Psi or Digamma function.
  Discussion:
    In Mathematica, the function can be evaluated by:
      PolyGamma[x]
    or
      Polygamma[0,x]
    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
    PSI(1) = -Euler's constant.
    PSI(X+1) = PSI(X) + 1 / X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2004
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
        -0.5772156649015329E+00,
        -0.4237549404110768E+00,
        -0.2890398965921883E+00,
        -0.1691908888667997E+00,
        -0.6138454458511615E-01,
        0.3648997397857652E-01,
        0.1260474527734763E+00,
        0.2085478748734940E+00,
        0.2849914332938615E+00,
        0.3561841611640597E+00,
        0.4227843350984671E+00
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
