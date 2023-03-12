#ifndef __DISABLEDEEP_ASA111

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ppnd ( void * data)
/******************************************************************************/
/*
  Purpose:
    PPND produces the normal deviate value corresponding to lower tail area = P.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2013
  Author:
    Original FORTRAN77 version by J Beasley, S Springer.
    C version by John Burkardt.
  Reference:
    J Beasley, S Springer,
    Algorithm AS 111:
    The Percentage Points of the Normal Distribution,
    Applied Statistics,
    Volume 26, Number 1, 1977, pages 118-121.
  Parameters:
    Input, double P, the value of the cumulative probability
    densitity function.  0 < P < 1.
    Output, integer *IFAULT, error flag.
    0, no error.
    1, P <= 0 or P >= 1.  PPND is returned as 0.
    Output, double PPND, the normal deviate value with the property that
    the probability of a standard normal deviate being less than or
    equal to PPND is P.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp p = *(ityp *) data;
		
	ityp r;
	ityp value;

	/*
	0.08 < P < 0.92
	*/
	if ( fabs ( p - 0.5 ) <= 0.42 )
	{
		r = ( p - 0.5 ) * ( p - 0.5 );
		value = ( p - 0.5 ) * ( ( ( -25.44106049637   * r + 41.39119773534 ) * r + -18.61500062529 ) * r + 2.50662823884 ) / ( ( ( ( 3.13082909833   * r + -21.06224101826 ) * r + 23.08336743743 ) * r + -8.47351093090 ) * r + 1.0 );
	}
	/*
	P < 0.08 or P > 0.92,
	R = MIN ( P, 1-P )
	*/
	else if ( 0.0 < p && p < 1.0 )
	{
		r = sqrt ( - log ( 0.5 < p ? 1.0 - p : p ) );
		value = ( ( ( 2.32121276858   * r + 4.85014127135 ) * r + -2.29796479134 ) * r + -2.78718931138 ) / ( ( 1.63706781897   * r + 3.54388924762 ) * r + 1.0 )*(1-((p<0.50)<<1));
	}
	/*
	P <= 0.0 or 1.0 <= P
	*/
	else
		value = 0.00;

	result = value;
	return &result;
}

#endif
