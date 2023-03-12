#ifndef __DISABLEDEEP_ASA066

#include "../dutils.h"


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_normal_01_cdf_inverse( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
  Discussion:
    The result is accurate to about 1 part in 10**7.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 December 2004
  Author:
    Original FORTRAN77 version by Michael Wichura.
    C version by John Burkardt.
  Reference:
    Michael Wichura,
    The Percentage Points of the Normal Distribution,
    Algorithm AS 241,
    Applied Statistics,
    Volume 37, Number 3, pages 477-484, 1988.
  Parameters:
    Input, ityp P, the value of the cumulative probability densitity function.
    0 < P < 1.  If P is outside this range, an "infinite" result is returned.
    Output, ityp r8_NORMAL_01_CDF_INVERSE, the normal deviate value with the
    property that the probability of a standard normal deviate being less than or
    equal to this value is P.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp p = *(ityp*) data;
	
	static ityp a[4] =
	{
		3.3871327179,
		50.434271938,
		159.29113202,
		59.109374720
	};
	static ityp b[4] =
	{
		1.0,
		17.895169469,
		78.757757664,
		67.187563600
	};
	static ityp c[4] =
	{
		1.4234372777,
		2.7568153900,
		1.3067284816,
		0.17023821103
	};
	static ityp d[3] =
	{
		1.0,
		0.73700164250,
		0.12021132975
	};
	static ityp e[4] =
	{
		6.6579051150,
		3.0812263860,
		0.42868294337,
		0.017337203997
	};
	static ityp f[3] =
	{
		1.0,
		0.24197894225,
		0.012258202635
	};
	ityp q;
	ityp r;
	ityp value;

	if ( p <= 0.00 || 1.00 <= p)
	{
		result = R8NORMAL01CDFINVERSE_INVALIDRETURNVALUE;
		return &result;
	}

	q = p - 0.50;

	if ( fabs ( q ) <= 0.425 )
	{
		r = 0.180625 - q * q;
		value = q * r8poly_value ( 4, a, r ) / r8poly_value ( 4, b, r );
	}
	else
	{
		r = q < 0.00 ? p : 1.00 - p;
		
		if ( r <= 0.00 )
		{
			result = R8NORMAL01CDFINVERSE_INVALIDRETURNVALUE;
			return &result;
		}

		r = sqrt ( -log ( r ) );

		if ( r <= 5.00 )
		{
			r = r - 1.60;
			value = r8poly_value ( 4, c, r ) / r8poly_value ( 3, d, r );
		}
		else
		{
			r = r - 5.00;
			value = r8poly_value ( 4, e, r ) / r8poly_value ( 3, f, r );
		}

		if ( q < 0.00 )
			value *= -1;
	}

	result = value; 
	return &result;
}

#endif
