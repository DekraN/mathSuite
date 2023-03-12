#ifndef __DISABLEDEEP_ASA091

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gammad ( void * data)
/******************************************************************************/
/*
  Purpose:
    GAMMAD computes the Incomplete Gamma Integral
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 November 2010
  Author:
    Original FORTRAN77 version by B Shea.
    C version by John Burkardt.
  Reference:
    B Shea,
    Algorithm AS 239:
    Chi-squared and Incomplete Gamma Integral,
    Applied Statistics,
    Volume 37, Number 3, 1988, pages 466-473.
  Parameters:
    Input, double X, P, the parameters of the incomplete
    gamma ratio.  0 <= X, and 0 < P.
    Output, int IFAULT, error flag.
    0, no error.
    1, X < 0 or P <= 0.
    Output, double GAMMAD, the value of the incomplete
    Gamma integral.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp p = a_data[1];
	
	ityp a;
	ityp an;
	ityp arg;
	ityp b;
	ityp c;
	ityp pn1;
	ityp pn2;
	ityp pn3;
	ityp pn4;
	ityp pn5;
	ityp pn6;
	ityp rn;
	dim_typ upper;
	ityp value = 0.00;

	/*
	Check the input.
	*/
	if ( x <= 0.00 || p <= 0.00)
	{
		result = GAMMAD_INVALIDRETURNVALUE;
		return &result;
	}

	/*
	If P is large, use a normal approximation.
	*/
	if ( 1000.0 < p )
	{
		result = alnorm ( 3.00 * sqrt ( p ) * ( pow ( x / p, 1.00 / 3.00 )  + 1.00 / ( 9.00 * p ) - 1.00 ), 0);
		return &result;
	}
	/*
	If X is large set value = 1.
	*/
	if ( 1.0E+08 < x )
	{
		result = GAMMAD_INVALIDRETURNVALUE;
		return &result; // X huge
	}
	
	/*
	Use Pearson's series expansion.
	(Note that P is not large enough to force overflow in ALOGAM).
	No need to test IFAULT on exit since P > 0.
	*/
	if ( x <= 1.0 || x < p )
	{
		arg = p * log ( x ) - x - lgamma ( p + 1.0 );
		c = 1.0;
		value = 1.0;
		a = p;

		for ( ; ; )
		{
			a += 1.00;
			c *= x / a;
			value += c;

			if ( c <= 1.0E-14 )
				break;
		}

		arg += log ( value );
		value = - 88.0 <= arg ? exp ( arg ) : 0.00;
	}
	/*
	Use a continued fraction expansion.
	*/
	else
	{
		arg = p * log ( x ) - x - lgamma ( p );
		a = 1.00 - p;
		b = a + x + 1.00;
		c = 0.00;
		pn1 = 1.00;
		pn2 = x;
		pn3 = x + 1.00;
		pn4 = x * b;
		value = pn3 / pn4;

		for ( ; ; )
		{
			a += 1.00;
			b +=  2.00;
			c += 1.00;
			an = a * c;
			pn5 = b * pn3 - an * pn1;
			pn6 = b * pn4 - an * pn2;

			if ( pn6 != 0.0 )
			{
				rn = pn5 / pn6;

				if ( fabs ( value - rn ) <= MIN ( 1.0E-14, 1.0E-14 * rn ) )
					break;
				value = rn;
			}

			pn1 = pn3;
			pn2 = pn4;
			pn3 = pn5;
			pn4 = pn6;
			/*
			Re-scale terms in continued fraction if terms are large.
			*/
			if ( 1.0E+37 <= fabs ( pn5 ) )
			{
				pn1 /= 1.0E+37;
				pn2 /= 1.0E+37;
				pn3 /= 1.0E+37;
				pn4 /= 1.0E+37;
			}
		}

		arg += log ( value );
	}

	result = 1.00 - exp(arg)*(- 88.0 <= arg);
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ppchi2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    PPCHI2 evaluates the percentage points of the Chi-squared PDF.
  Discussion
    Incorporates the suggested changes in AS R85 (vol.40(1),
    pages 233-5, 1991) which should eliminate the need for the limited
    range for P, though these limits have not been removed
    from the routine.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2013
  Author:
    Original FORTRAN77 version by Donald Best, DE Roberts.
    C version by John Burkardt.
  Reference:
    Donald Best, DE Roberts,
    Algorithm AS 91:
    The Percentage Points of the Chi-Squared Distribution,
    Applied Statistics,
    Volume 24, Number 3, 1975, pages 385-390.
  Parameters:
    Input, double P,  value of the chi-squared cumulative
    probability density function.
    0.000002 <= P <= 0.999998.
    Input, double V, the parameter of the chi-squared probability
    density function.
    0 < V.
    Input, double G, the value of log ( Gamma ( V / 2 ) ).
    Output, int *IFAULT, is nonzero if an error occurred.
    0, no error.
    1, P is outside the legal range.
    2, V is not positive.
    3, an error occurred in GAMMAD.
    4, the result is probably as accurate as the machine will allow.
    Output, double PPCHI2, the value of the chi-squared random
    deviate with the property that the probability that a chi-squared random
    deviate with parameter V is less than or equal to PPCHI2 is P.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp p = a_data[0];
	const register ityp v = a_data[1];
	const register ityp g = a_data[2];
	
	ityp a;
	ityp b;
	ityp c;
	ityp ch;
	dim_typ i;
	dim_typ if1;
	ityp p1;
	ityp p2;
	ityp q;
	ityp s1;
	ityp s2;
	ityp s3;
	ityp s4;
	ityp s5;
	ityp s6;
	ityp t;
	ityp x;
	ityp xx;
	/*
	Test arguments and initialize.
	*/

	if ( p < 0.000002 || 0.999998 < p || v <= 0.00)
	{
		result = PPCHI2_INVALIDRETURNVALUE;
		return &result;
	}

	xx = 0.5 * v;
	c = xx - 1.0;
	/*
	Starting approximation for small chi-squared
	*/
	if ( v < - 1.24 * log ( p ) )
	{
		ch = pow ( p * xx * exp ( g + xx * 0.6931471806 ), 1.0 / xx );
		if ( ch < 0.5E-06 )
		{
			result = ch;
			return &result;
		}
	}
	/*
	Starting approximation for V less than or equal to 0.32
	*/
	else if ( v <= 0.32 )
	{
		ch = 0.40;
		a = log ( 1.00 - p );
		for ( ; ; )
		{
			q = ch;
			p1 = 1.00 + ch * ( 4.67 + ch );
			p2 = ch * (6.73 + ch * ( 6.66 + ch ) );
			t = - 0.50 + (4.67 + 2.00 * ch ) / p1 - ( 6.73 + ch * ( 13.32 + 3.0 * ch ) ) / p2;
			ch -= ( 1.00 - exp ( a + g + 0.5 * ch + c * 0.6931471806 ) * p2 / p1) / t;
			if ( fabs ( q / ch - 1.00 ) <= 0.01 )
				break;
		}
	}
	else
	{
		/*
		Call to algorithm AS 111 - note that P has been tested above.
		AS 241 could be used as an alternative.
		*/
		x = ppnd ( p );
		/*
		Starting approximation using Wilson and Hilferty estimate
		*/
		p1 = 0.222222 / v;
		ch = v * pow ( x * sqrt ( p1 ) + 1.0 - p1, 3 );
		/*
		Starting approximation for P tending to 1.
		*/
		if ( 2.20 * v + 6.00 < ch )
			ch = - 2.0 * ( log ( 1.0 - p ) - c * log ( 0.5 * ch ) + g );
	}
	/*
	Call to algorithm AS 239 and calculation of seven term
	Taylor series
	*/
	for ( i = 1; i <= 20; ++i )
	{
		q = ch;
		p1 = 0.50 * ch;
		p2 = p - gammad ( p1, xx );

		if (if1)
		{
			result = PPCHI2_INVALIDRETURNVALUE;
			return &result;
		}

		t = p2 * exp ( xx * 0.6931471806 + g + p1 - c * log ( ch ) );
		b = t / ch;
		a = 0.50 * t - b * c;
		s1 = ( 210.0 + a * ( 140.0 + a * ( 105.0 + a * ( 84.0 + a * ( 70.0 +
		60.0 * a ))))) / 420.0;
		s2 = ( 420.0 + a * ( 735.0 + a * ( 966.0 + a * ( 1141.0 + 1278.0 * a )))) / 2520.0;
		s3 = ( 210.0 + a * ( 462.0 + a * ( 707.0 + 932.0 * a ))) / 2520.0;
		s4 = ( 252.0 + a * ( 672.0 + 1182.0 * a) + c * ( 294.0 + a * ( 889.0 + 1740.0 * a ))) / 5040.0;
		s5 = ( 84.0 + 264.0 * a + c * ( 175.0 + 606.0 * a )) / 2520.0;
		s6 = ( 120.0 + c * ( 346.0 + 127.0 * c )) / 5040.0;
		ch = ch + t * ( 1.00 + 0.50 * t * s1 - b * c * ( s1 - b *
		( s2 - b * ( s3 - b * ( s4 - b * ( s5 - b * s6 ))))));

		if ( 0.5E-06 < fabs ( q / ch - 1.0 ) )
		{
			result = ch;
			return &result;
		}
	}

	result = ch;
	return &result;
}

#endif
