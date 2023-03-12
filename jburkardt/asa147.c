#ifndef __DISABLEDEEP_ASA147

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gammds( void * data)
/******************************************************************************/
/*
  Purpose:
    GAMMDS computes the incomplete Gamma integral.
  Discussion:
    The parameters must be positive.  An infinite series is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2010
  Author:
    Original FORTRAN77 version by Chi Leung Lau.
    C version by John Burkardt.
  Reference:
    Chi Leung Lau,
    Algorithm AS 147:
    A Simple Series for the Incomplete Gamma Integral,
    Applied Statistics,
    Volume 29, Number 1, 1980, pages 113-114.
  Parameters:
    Input, double X, P, the arguments of the incomplete
    Gamma integral.  X and P must be greater than 0.
    Output, int *IFAULT, error flag.
    0, no errors.
    1, X <= 0 or P <= 0.
    2, underflow during the computation.
    Output, double GAMMDS, the value of the incomplete
    Gamma integral.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp p = a_data[1];
		
	ityp a;
	ityp arg;
	ityp c;
	ityp f;
	ityp value;
	/*
	Check the input.
	*/
	if ( x <= 0.0 || p <= 0.00 || (arg = p * log ( x ) - lgamma ( p + 1.00 ) - x) < log ( 1.0E-37 ) || !f)
	{
		result = GAMMDS_INVALIDRETURNVALUE;
		return &result;
	}
	
	/*
	LGAMMA is the natural logarithm of the gamma function.
	*/

	f = exp ( arg );

	/*
	Series begins.
	*/
	c = 1.00;
	value = 1.00;
	a = p;

	for ( ; ; )
	{
		a += 1.00;
		c *= x/a;
		value +=c;

		if ( c <= 1.0E-09 * value )
			break;
	}
	
	result = value*f;
	return &result;
}
#endif
