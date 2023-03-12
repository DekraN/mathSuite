#ifndef __DISABLEDEEP_ASA245

#include "../dutils.h"

//****************************************************************************80
__MATHSUITE __JBURKARDT  void *   _alngam ( void * data)
//****************************************************************************80
//
//  Purpose:
//    ALNGAM computes the logarithm of the gamma function.
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    13 January 2008
//  Author:
//    Original FORTRAN77 version by Allan Macleod.
//    C++ version by John Burkardt.
//  Reference:
//    Allan Macleod,
//    Algorithm AS 245,
//    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
//    Applied Statistics,
//    Volume 38, Number 2, 1989, pages 397-402.
//  Parameters:
//    Input, double XVALUE, the argument of the Gamma function.
//    Output, int IFAULT, error flag.
//    0, no error occurred.
//    1, XVALUE is less than or equal to 0.
//    2, XVALUE is too big.
//    Output, double ALNGAM, the logarithm of the gamma function of X.
{
	static ityp result = MAX_VAL;
	
	const register ityp xvalue = *(ityp*) data;
	
	ityp r1[9] =
	{
		-2.66685511495,
		-24.4387534237,
		-21.9698958928,
		11.1667541262,
		3.13060547623,
		0.607771387771,
		11.9400905721,
		31.4690115749,
		15.2346874070
	};
	ityp r2[9] =
	{
		-78.3359299449,
		-142.046296688,
		137.519416416,
		78.6994924154,
		4.16438922228,
		47.0668766060,
		313.399215894,
		263.505074721,
		43.3400022514
	};
	ityp r3[9] =
	{
		-2.12159572323E+05,
		2.30661510616E+05,
		2.74647644705E+04,
		-4.02621119975E+04,
		-2.29660729780E+03,
		-1.16328495004E+05,
		-1.46025937511E+05,
		-2.42357409629E+04,
		-5.70691009324E+02
	};
	ityp r8[5] =
	{
		0.279195317918525,
		0.4917317610505968,
		0.0692910599291889,
		3.350343815022304,
		6.012459259764103
	};

	ityp value = 0.00;
	ityp x, y, x1, x2;
	x = xvalue;

	//
	//  Check the input.
	//
	if ( 1.0E+30 <= x  || x <= 0.00)
	{
		result = ALNGAM_INVALIDRETURNVALUE;
		return &result;
	}

	//
	//  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
	//
	if ( x < 1.5 )
	{
		if(x < 0.50)
		{
			if (x + 1.00 == 1.00 )
			{
				result = -log(x);
				return &result;
			}
		}
		else
		{
			value = 0.00;
			y = x;
			x = ( x - 0.50 ) - 0.50;
		}

		result = value + x * ((((r1[4] * y + r1[3] ) * y + r1[2] ) * y + r1[1] ) * y + r1[0] ) / ((((y + r1[8] ) * y + r1[7] ) * y + r1[6] ) * y + r1[5] );
		return &result;
	}

	//
	//  Calculation for 1.5 <= X < 4.0.
	//
	if ( x < 4.00 )
	{
		y = ( x - 1.00 ) - 1.00;
		value = y * ((((r2[4] * x + r2[3] ) * x + r2[2] ) * x + r2[1] ) * x + r2[0] ) / ((((x + r2[8] ) * x + r2[7] ) * x + r2[6] ) * x + r2[5] );
	}
	//
	//  Calculation for 4.0 <= X < 12.0.
	//
	else if ( x < 12.0 )
		value = ((((r3[4]*x+ r3[3] ) * x + r3[2] ) * x + r3[1] ) * x + r3[0] ) / ((((x + r3[8] ) * x + r3[7] ) * x + r3[6] ) * x + r3[5] );
	//
	//  Calculation for 12.0 <= X.
	//
	else
	{
		y = log ( x );
		value = x * ( y - 1.0 ) - 0.5 * y + 0.918938533204673;

		if ( x <= 510000.0 )
		{
			x1 = 1.00 / x;
			x2 = x1 * x1;
			value += x1 * ( (r8[2] * x2 + r8[1] ) * x2 + r8[0] ) / ( (x2 + r8[4] ) * x2 + r8[3] );
		}
	}
	
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _gamma_log_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    GAMMA_LOG_VALUES returns some values of the Log Gamma function.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Log[Gamma[x]]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 August 2004
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
	
    # define N_MAX 20

    ityp fx_vec[N_MAX] =
    {
        0.1524063822430784E+01,
        0.7966778177017837E+00,
        0.3982338580692348E+00,
        0.1520596783998375E+00,
        0.0000000000000000E+00,
        -0.4987244125983972E-01,
        -0.8537409000331584E-01,
        -0.1081748095078604E+00,
        -0.1196129141723712E+00,
        -0.1207822376352452E+00,
        -0.1125917656967557E+00,
        -0.9580769740706586E-01,
        -0.7108387291437216E-01,
        -0.3898427592308333E-01,
        0.00000000000000000E+00,
        0.69314718055994530E+00,
        0.17917594692280550E+01,
        0.12801827480081469E+02,
        0.39339884187199494E+02,
        0.71257038967168009E+02
    };

    ityp x_vec[N_MAX] =
    {
        0.20E+00,
        0.40E+00,
        0.60E+00,
        0.80E+00,
        1.00E+00,
        1.10E+00,
        1.20E+00,
        1.30E+00,
        1.40E+00,
        1.50E+00,
        1.60E+00,
        1.70E+00,
        1.80E+00,
        1.90E+00,
        2.00E+00,
        3.00E+00,
        4.00E+00,
        10.00E+00,
        20.00E+00,
        30.00E+00
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
