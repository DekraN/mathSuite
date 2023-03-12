#ifndef __DISABLEDEEP_ASA005

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _alnorm (void *data) 
/******************************************************************************/
/*
  Purpose:
    ALNORM computes the cumulative density of the standard normal distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2010
  Author:
    Original FORTRAN77 version by David Hill.
    C version by John Burkardt.
  Reference:
    David Hill,
    Algorithm AS 66:
    The Normal Integral,
    Applied Statistics,
    Volume 22, Number 3, 1973, pages 424-427.
  Parameters:
    Input, double X, is one endpoint of the semi-infinite interval
    over which the integration takes place.
    Input, int UPPER, determines whether the upper or lower
    interval is to be integrated:
    1  => integrate from X to + Infinity;
    0 => integrate from - Infinity to X.
    Output, double ALNORM, the integral of the standard normal
    distribution over the desired interval.
*/
{
	static ityp result = MAX_VAL;
	
	const itb * const s_data = data;
	register ityp x = s_data->a0;
	bool upper = s_data->a1;
	
	ityp value;
	ityp y;

	if(x < 0.00)
	{
		upper = !upper;
		x *= -1;
	}

	if ( 7.00 < x && ( ( !upper ) || 18.66 < x) )
	{
		result = 0.00 + !upper;
		return &result;
	}

	y = 0.5*x*x;
	value = x <= 1.28 ? 0.5 - x * ( 0.398942280444 - 0.39990348504 * y / ( y + 5.75885480458 + -29.8213557807 / ( y + 2.62433121679 + 48.6959930692 / ( y + 5.92885724438 )))) :
	0.398942280385 * exp ( - y ) / ( x + -0.000000038052 + 1.00000615302 / ( x + 0.000398064794 + 1.98615381364 / ( x + -0.151679116635 + 5.29330324926 / ( x + 4.8385912808 + -15.1508972451 / ( x + 0.742380924027 + 30.789933034 / ( x + 3.99019417011 ))))));

	if (!upper)
		value = 1.00 - value;

	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _prncst (void *data)
/******************************************************************************/
/*
  Purpose:
    PRNCST computes the lower tail of noncentral T distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2010
  Author:
    Original FORTRAN77 version by BE Cooper.
    C version by John Burkardt.
  Reference:
    BE Cooper,
    Algorithm AS 5:
    The Integral of the Non-Central T-Distribution,
    Applied Statistics,
    Volume 17, Number 2, 1968, page 193.
  Parameters:
    Input, double ST, the argument.
    Input, int IDF, the number of degrees of freedom.
    Input, double D, the noncentrality parameter.
    Output, int *IFAULT, error flag.
    0, no error occurred.
    nonzero, an error occurred.
    Output, double PRNCST, the value of the lower tail of
    the noncentral T distribution.
  Local Parameters:
    Local, double G1, 1.0 / sqrt(M_2TPI)
    Local, double G2, 1.0 / (M_2TPI)
    Local, double G3, sqrt(M_2TPI)
*/
{
	static ityp result = MAX_VAL;
	
	const _2itdt * const s_data = data;
	
	const register ityp st = s_data->a0;
	const register ityp d = s_data->a1;
	const register dim_typ idf = s_data->a2;
	
	
	ityp a;
	ityp ak;
	ityp b;
	ityp da;
	ityp drb;
	ityp f;
	ityp fk;
	ityp fkm1;
	ityp fmkm1;
	ityp fmkm2;
	dim_typ ioe;
	dim_typ k;
	ityp rb;
	ityp sum;

	f = (ityp) (idf);
	/*
	For very large IDF, use the normal approximation.
	*/
	if (100 < idf )
	{
		a = sqrt ( 0.50 * f ) * exp ( lgamma ( 0.50 * ( f - 1.00 ) ) - lgamma ( 0.50 * f ) ) * d;
		result = alnorm ( ( st - a ) / sqrt ( f * ( 1.00 + d * d ) / ( f - 2.00 ) - a * a ), 0 );
		return &result;
	}

	ioe = (idf % 2);
	a = st / sqrt (f);
	b = f / (f + st * st);
	rb = sqrt (b);
	da = d * a;
	drb = d * rb;

	if (idf == 1)
	{
		// tfn ( drb, a );
		result = alnorm ( drb, 1 ) + 2.00 * tfn( drb, a);
		return &result;
	}

	fmkm2 = fabs ( drb ) < 12.50 ? a * rb * exp ( - 0.50 * drb * drb ) * alnorm ( a * drb, 0 ) * 0.3989422804 : 0.00;
	fmkm1 = b*da*fmkm2 + (fabs ( d ) < 12.50)*(fmkm1+b*a*0.1591549431*exp(-0.50*d*d));
	sum = ioe ? fmkm1 : fmkm2;

	ak = 1.00;
	fk = 2.00;

	for ( k = 2; k <= idf - 2; k +=2 )
	{
		fkm1 = fk - 1.00;
		fmkm2 = b * ( da * ak * fmkm1 + fmkm2 ) * fkm1 / fk;
		ak = 1.00 / ( ak * fkm1 );
		fmkm1 = b * ( da * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0 );

		sum += (ioe ? fmkm1 : fmkm2);
		ak = 1.00 / ( ak * fk );
		fk += 2.00;
	}

	result = (ioe ? alnorm ( drb, 1 ) + 2.00 * ( sum + tfn ( drb, a ) ) : alnorm ( d, 1 ) + sum * 2.5066282746);
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tfn (void *data)
/******************************************************************************/
/*
  Purpose:
    TFN calculates the T-function of Owen.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 January 2008
  Author:
    Original FORTRAN77 version by JC Young, Christoph Minder.
    C version by John Burkardt.
  Reference:
    MA Porter, DJ Winstanley,
    Remark AS R30:
    A Remark on Algorithm AS76:
    An Integral Useful in Calculating Noncentral T and Bivariate
    Normal Probabilities,
    Applied Statistics,
    Volume 28, Number 1, 1979, page 113.
    JC Young, Christoph Minder,
    Algorithm AS 76:
    An Algorithm Useful in Calculating Non-Central T and
    Bivariate Normal Distributions,
    Applied Statistics,
    Volume 23, Number 3, 1974, pages 455-457.
  Parameters:
    Input, double X, FX, the parameters of the function.
    Output, double TFN, the value of the T-function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp fx = a_data[1];
	
	ityp fxs;
	dim_typ i;
	const ityp r[5] =
	{
		0.1477621,
		0.1346334,
		0.1095432,
		0.0747257,
		0.0333357
	};
	ityp r1;
	ityp r2;
	ityp rt;
	const ityp u[5] =
	{
		0.0744372,
		0.2166977,
		0.3397048,
		0.4325317,
		0.4869533
	};
	ityp x1;
	ityp x2;
	ityp xs;
	/*
	Test for X near zero.
	*/
	if (fabs (x) < 1.0E-35)
	{
		result = 0.159155 * atan ( fx );
		return &result;
	}
	/*
	Test for large values of abs(X).
	*/
	/*
	Test for FX near zero.
	*/
	if (15.00 < fabs (x) || fabs ( fx ) < 1.0E-35)
	{
		result = 0.00;
		return &result;
	}

	/*
	Test whether abs ( FX ) is so large that it must be truncated.
	*/
	xs = - 0.50*x*x;
	x2 = fx;
	fxs = fx*fx;
	/*
	Computation of truncation point by Newton iteration.
	*/
	if ( 15.00 <= log ( 1.00 + fxs ) - xs * fxs )
	{
		x1 = 0.50 * fx;
		fxs = 0.25 * fxs;

		for ( ; ; )
		{
			rt = fxs + 1.0;
			x2 = x1 + ( xs * fxs + 15.00 - log ( rt ) ) / ( 2.00 * x1 * ( 1.00 / rt - xs ) );
			fxs = x2 * x2;

			if ( fabs ( x2 - x1 ) < 1.0E-05)
				break;
			x1 = x2;
		}
	}
	/*
	Gaussian quadrature.
	*/
	rt = 0.00;
	for ( i = 0; i < 5; ++i )
	{
		r1 = 1.00 + fxs * pow ( 0.50 + u[i], 2 );
		r2 = 1.00 + fxs * pow ( 0.50 - u[i], 2 );

		rt += r[i] * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 );
	}

	result = rt * x2 * 0.159155;
	return &result;
}

#endif
