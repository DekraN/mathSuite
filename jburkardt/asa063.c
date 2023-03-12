#ifndef __DISABLEDEEP_ASA063

#include "../dutils.h"



/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _betain ( void * data)
/******************************************************************************/
/*
  Purpose:
    BETAIN computes the incomplete Beta function ratio.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 November 2010
  Author:
    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
    C version by John Burkardt.
  Reference:
    KL Majumder, GP Bhattacharjee,
    Algorithm AS 63:
    The incomplete Beta Integral,
    Applied Statistics,
    Volume 22, Number 3, 1973, pages 409-411.
  Parameters:
    Input, double X, the argument, between 0 and 1.
    Input, double P, Q, the parameters, which
    must be positive.
    Input, double BETA, the logarithm of the complete
    beta function
    Output, int *IFAULT, error flag.
    0, no error.
    nonzero, an error occurred.
    Output, double BETAIN, the value of the incomplete
    Beta function ratio.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp p = a_data[1];
	const register ityp q = a_data[2];
	const register ityp beta = a_data[3];
	
	ityp ai;
	ityp betainx;
	ityp cx;
	dim_typ indx;
	dim_typ ns;
	ityp pp;
	ityp psq;
	ityp qq;
	ityp rx;
	ityp temp;
	ityp term;
	ityp value;
	ityp xx;

	value = x;
	/*
	Check the input arguments.
	*/
	if (p <= 0.00 || q <= 0.00 || x < 0.00 || 1.00 < x || x == 0.00 || x == 1.00)
	{
		result = BETAIN_INVALIDRETURNVALUE;
		return &result;
	}

	/*
	Change tail if necessary and determine S.
	*/
	psq = p + q;
	cx = 1.00-x;

	if (p < psq * x )
	{
		xx = cx;
		cx = x;
		pp = q;
		qq = p;
		indx = 1;
	}
	else
	{
		xx = x;
		pp = p;
		qq = q;
		indx = 0;
	}

	term = ai = value = 1.00;
	ns = (dim_typ)(qq + cx*psq);
	/*
	Use the Soper reduction formula.
	*/
	rx = xx / cx;
	temp = qq - ai;

	if (!ns)
		rx = xx;

	for ( ; ; )
	{
		term = term * temp * rx / ( pp + ai );
		value = value + term;
		temp = fabs (term);

		if ( temp <= 0.1E-14 && temp <= 0.1E-14 * value )
		{
			value *= exp ( pp * log ( xx ) + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;

			if (indx)
				value = 1.00 - value;
			break;
		}

		ai += 1.00;
		-- ns;

		if ( 0 <= ns )
		{
			temp = qq - ai;
			if (!ns)
			rx = xx;
		}
		else
		{
			temp = psq;
			psq += 1.00;
		}
	}

	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _beta_inc_values ( void *data)
/******************************************************************************/
/*
  Purpose:
    BETA_INC_VALUES returns some values of the incomplete Beta function.
  Discussion:
    The incomplete Beta function may be written
      BETA_INC(A,B,X) = Integral (0 <= T <= X) T^(A-1) * (1-T)^(B-1) dT
                      / Integral (0 <= T <= 1) T^(A-1) * (1-T)^(B-1) dT
    Thus,
      BETA_INC(A,B,0.0) = 0.0;
      BETA_INC(A,B,1.0) = 1.0
    The incomplete Beta function is also sometimes called the
    "modified" Beta function, or the "normalized" Beta function
    or the Beta CDF (cumulative density function.
    In Mathematica, the function can be evaluated by:
      BETA[X,A,B] / BETA[A,B]
    The function can also be evaluated by using the Statistics package:
      Needs["Statistics`ContinuousDistributions`"]
      dist = BetaDistribution [ a, b ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 April 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Karl Pearson,
    Tables of the Incomplete Beta Function,
    Cambridge University Press, 1968.
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
    Output, double *A, B, the parameters of the function.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x = s_data->a3; 
	ityp * fx = s_data->a4;
	
    # define N_MAX 45

    static ityp a_vec[N_MAX] =
    {
        0.5E+00,
        0.5E+00,
        0.5E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        5.5E+00,
        10.0E+00,
        10.0E+00,
        10.0E+00,
        10.0E+00,
        20.0E+00,
        20.0E+00,
        20.0E+00,
        20.0E+00,
        20.0E+00,
        30.0E+00,
        30.0E+00,
        40.0E+00,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.2E+01,
        0.3E+01,
        0.4E+01,
        0.5E+01,
        1.30625,
        1.30625,
        1.30625
    };

    static ityp b_vec[N_MAX] =
    {
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        1.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        5.0E+00,
        0.5E+00,
        5.0E+00,
        5.0E+00,
        10.0E+00,
        5.0E+00,
        10.0E+00,
        10.0E+00,
        20.0E+00,
        20.0E+00,
        10.0E+00,
        10.0E+00,
        20.0E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.2E+01,
        0.3E+01,
        0.4E+01,
        0.5E+01,
        0.2E+01,
        0.2E+01,
        0.2E+01,
        0.2E+01,
        11.7562,
        11.7562,
        11.7562
    };

    static ityp fx_vec[N_MAX] =
    {
        0.6376856085851985E-01,
        0.2048327646991335E+00,
        0.1000000000000000E+01,
        0.0000000000000000E+00,
        0.5012562893380045E-02,
        0.5131670194948620E-01,
        0.2928932188134525E+00,
        0.5000000000000000E+00,
        0.2800000000000000E-01,
        0.1040000000000000E+00,
        0.2160000000000000E+00,
        0.3520000000000000E+00,
        0.5000000000000000E+00,
        0.6480000000000000E+00,
        0.7840000000000000E+00,
        0.8960000000000000E+00,
        0.9720000000000000E+00,
        0.4361908850559777E+00,
        0.1516409096347099E+00,
        0.8978271484375000E-01,
        0.1000000000000000E+01,
        0.5000000000000000E+00,
        0.4598773297575791E+00,
        0.2146816102371739E+00,
        0.9507364826957875E+00,
        0.5000000000000000E+00,
        0.8979413687105918E+00,
        0.2241297491808366E+00,
        0.7586405487192086E+00,
        0.7001783247477069E+00,
        0.5131670194948620E-01,
        0.1055728090000841E+00,
        0.1633399734659245E+00,
        0.2254033307585166E+00,
        0.3600000000000000E+00,
        0.4880000000000000E+00,
        0.5904000000000000E+00,
        0.6723200000000000E+00,
        0.2160000000000000E+00,
        0.8370000000000000E-01,
        0.3078000000000000E-01,
        0.1093500000000000E-01,
        0.918884684620518,
        0.21052977489419,
        0.1824130512500673
    };

    static ityp x_vec[N_MAX] =
    {
        0.01E+00,
        0.10E+00,
        1.00E+00,
        0.00E+00,
        0.01E+00,
        0.10E+00,
        0.50E+00,
        0.50E+00,
        0.10E+00,
        0.20E+00,
        0.30E+00,
        0.40E+00,
        0.50E+00,
        0.60E+00,
        0.70E+00,
        0.80E+00,
        0.90E+00,
        0.50E+00,
        0.90E+00,
        0.50E+00,
        1.00E+00,
        0.50E+00,
        0.80E+00,
        0.60E+00,
        0.80E+00,
        0.50E+00,
        0.60E+00,
        0.70E+00,
        0.80E+00,
        0.70E+00,
        0.10E+00,
        0.20E+00,
        0.30E+00,
        0.40E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.30E+00,
        0.30E+00,
        0.30E+00,
        0.30E+00,
        0.225609,
        0.0335568,
        0.0295222
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *a = *x = *b = *fx = 0.00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
