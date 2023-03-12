#ifndef __DISABLEDEEP_ASA243

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tnc ( void * data)
/******************************************************************************/
/*
  Purpose:
    TNC computes the tail of the noncentral T distribution.
  Discussion:
    This routine computes the cumulative probability at T of the
    non-central T-distribution with DF degrees of freedom (which may
    be fractional) and non-centrality parameter DELTA.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 November 2010
  Author:
    Original FORTRAN77 version by Russell Lenth.
    C version by John Burkardt.
  Reference:
    Russell Lenth,
    Algorithm AS 243:
    Cumulative Distribution Function of the Non-Central T Distribution,
    Applied Statistics,
    Volume 38, Number 1, 1989, pages 185-189.
    William Guenther,
    Evaluation of probabilities for the noncentral distributions and
    difference of two T-variables with a desk calculator,
    Journal of Statistical Computation and Simulation,
    Volume 6, Number 3-4, 1978, pages 199-206.
  Parameters:
    Input, double T, the point whose cumulative probability
    is desired.
    Input, double DF, the number of degrees of freedom.
    Input, double DELTA, the noncentrality parameter.
    Output, int *IFAULT, error flag.
    0, no error.
    nonzero, an error occcurred.
    Output, double TNC, the tail of the noncentral
    T distribution.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp t = a_data[0];
	const register ityp df = a_data[1];
	const register ityp delta = a_data[2];
	
	ityp a;
	ityp albeta;
	ityp b;
	ityp del = delta;
	ityp en;
	ityp errbd;
	ityp geven;
	ityp godd;
	ityp half;
	ityp lambda;
	bool negdel = false;
	ityp one;
	ityp p;
	ityp q;
	ityp rxb;
	ityp s;
	ityp tt = t;
	ityp two;
	ityp value = 0.00;
	ityp x;
	ityp xeven;
	ityp xodd;
	ityp zero;

	if ( df <= 0.00)
	{
		result = TNC_INVALIDRETURNVALUE;
		return &result;
	}

	if ( t < 0.00 )
	{
		negdel = true;
		tt *= -1;
		del *= -1;
	}
	/*
	Initialize twin series.
	*/
	en = 1.00;
	x = t * t / ( t * t + df );

	if ( x <= 0.00 )
	{
		value += alnorm(del, 1);

		if (negdel)
		{
			result = 1.00 - value;
			return &result;
		}
	}

	lambda = del * del;
	p = 0.50 * exp ( - 0.50 * lambda );
	q = 0.79788456080286535588 * p * del;
	s = 0.50 - p;
	a = 0.50;
	b = 0.50 * df;
	rxb = pow ( 1.00 - x, b );
	albeta = 0.57236494292470008707 + lgamma ( b ) - lgamma ( a + b );
	xodd = betain ( x, a, b, albeta);
	godd = 2.00 * rxb * exp ( a * log ( x ) - albeta );
	xeven = 1.00 - rxb;
	geven = b * x * rxb;
	value = p * xodd + q * xeven;
	/*
	Repeat until convergence.
	*/
	for ( ; ; )
	{
		++ a;
		xodd = xodd - godd;
		xeven = xeven - geven;
		godd = godd * x * ( a + b - 1.00 ) / a;
		geven = geven * x * ( a + b - 0.50 ) / ( a + 0.50 );
		p *= lambda / ( 2.00 * en );
		q *= lambda / ( 2.00 * en + 1.00 );
		s -= p;
		++ en;
		value =+ p * xodd + q * xeven;
		errbd = 2.00 * s * ( xodd - godd );

		if ( errbd <= 1.0E-10 || 100 < en)
			break;
	}

	value += alnorm ( del, 1 );
	result = negdel ? 1.00 - value : value;
	return &result;
}

#endif
