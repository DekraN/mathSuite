#ifndef __DISABLEDEEP_ASA109

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _xinbta ( void * data)
/******************************************************************************/
/*
  Purpose:
    XINBTA computes inverse of the incomplete Beta function.
  Discussion:
    The accuracy exponent SAE was loosened from -37 to -30, because
    the code would not otherwise accept the results of an iteration
    with p = 0.3, q = 3.0, alpha = 0.2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2014
  Author:
    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
    C version by John Burkardt.
  Reference:
    GW Cran, KJ Martin, GE Thomas,
    Remark AS R19 and Algorithm AS 109:
    A Remark on Algorithms AS 63: The Incomplete Beta Integral
    and AS 64: Inverse of the Incomplete Beta Integeral,
    Applied Statistics,
    Volume 26, Number 1, 1977, pages 111-114.
  Parameters:
    Input, double P, Q, the parameters of the incomplete
    Beta function.
    Input, double BETA, the logarithm of the value of
    the complete Beta function.
    Input, double ALPHA, the value of the incomplete Beta
    function.  0 <= ALPHA <= 1.
    Output, int *IFAULT, error flag.
    0, no error occurred.
    nonzero, an error occurred.
    Output, double XINBTA, the argument of the incomplete
    Beta function which produces the value ALPHA.
  Local Parameters:
    Local, double SAE, accuracy is requested to about 10^SAE.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp p = a_data[0];
	const register ityp q = a_data[1];
	const register ityp beta = a_data[2];
	const register ityp alpha = a_data[3];
	
	ityp a;
	ityp acu;
	ityp adj;
	ityp fpu;
	ityp g;
	ityp h;
	bool indx;
	ityp pp;
	ityp prev;
	ityp qq;
	ityp r;
	ityp s;
	ityp sq;
	ityp t;
	ityp tx;
	ityp value;
	ityp w;
	ityp xin;
	ityp y;
	ityp yprev;

	fpu = pow ( 10.0, -30.0 );

	value = alpha;
	/*
	Test for admissibility of parameters.
	*/
	if ( p <= 0.00 || q <= 0.00 || alpha <= 0.00 || alpha >= 1.00)
	{
		result = XINBTA_INVALIDRETURNVALUE;
		return &result;
	}

	/*
	Change tail if necessary.
	*/
	if ( 0.5 < alpha )
	{
		a = 1.00 - alpha;
		pp = q;
		qq = p;
		indx = true;
	}
	else
	{
		a = alpha;
		pp = p;
		qq = q;
		indx = false;
	}
	/*
	Calculate the initial approximation.
	*/
	r = sqrt ( - log ( a * a ) );
	y = r - ( 2.30753 + 0.27061 * r ) / ( 1.00 + ( 0.99229 + 0.04481 * r ) * r );

	if ( 1.00 < pp && 1.00 < qq )
	{
		r = ( y * y - 3.00 ) / 6.00;
		s = 1.00 / ( pp + pp - 1.00 );
		t = 1.00 / ( qq + qq - 1.00 );
		h = 2.00 / ( s + t );
		w = y * sqrt ( h + r ) / h - ( t - s )
		* ( r + 5.00 / 6.00 - 2.00 / ( 3.00 * h ) );
		value = pp / ( pp + qq * exp ( w + w ) );
	}
	else
	{
		r = qq + qq;
		t = 1.00 / ( 9.00 * qq );
		t = r * pow ( 1.00 - t + y * sqrt ( t ), 3 );

		if ( t <= 0.00 )
			value = 1.00 - exp ( ( log ( ( 1.00 - a ) * qq ) + beta ) / qq );
		else
		{
			t = ( 4.00 * pp + r - 2.00 ) / t;
			value = t<= 1.00 ? exp ( ( log ( a * pp ) + beta ) / pp ) : 1.0 - 2.0 / ( t + 1.0 );
		}
	}
	/*
	Solve for X by a modified Newton-Raphson method,
	using the function BETAIN.
	*/
	r = 1.00 - pp;
	t = 1.00 - qq;
	yprev = 0.00;
	sq = 1.00;
	prev = 1.00;

	if ( value < 0.0001 )
		value = 0.0001;

	if ( 0.9999 < value )
		value = 0.9999;

	acu = pow ( 10.0, MAX ( - 5.0 / pp / pp - 1.0 / pow ( a, 0.2 ) - 13.0, -30.0 ) );
	/*
	Iteration loop.
	*/
	for ( ; ; )
	{
		y = betain ( value, pp, qq, beta );

		if (y == BETAIN_INVALIDRETURNVALUE)
		{
			result = XINBTA_INVALIDRETURNVALUE;
			return &result;
		}

		xin = value;
		y = ( y - a ) * exp ( beta + r * log ( xin ) + t * log ( 1.00 - xin ) );

		if ( y * yprev <= 0.00 )
			prev = MAX ( sq, fpu );

		g = 1.00;

		for ( ; ; )
		{
		/*
		Choose damping factor.
		*/
			for ( ; ; )
			{
				adj = g * y;
				sq = adj * adj;

				if ( sq < prev )
				{
					tx = value - adj;

					if ( 0.00 <= tx && tx <= 1.00 )
						break;
				}
				g /= 3.00;
			}
			/*
			Check whether current estimate is acceptable.
			The change "VALUE = TX" was suggested by Ivan Ukhov.
			*/
			if ( prev <= acu || y * y <= acu )
			{
				result = indx ? 1.00 - value : tx;
				return &result;
			}

			if ( tx != 0.00 && tx != 1.00 )
				break;

			g /= 3.00;
		}

		if ( tx == value )
			break;

		value = tx;
		yprev = y;
	}
	
	result = indx ? 1.00 - value : value;
	return &result;
}

#endif
