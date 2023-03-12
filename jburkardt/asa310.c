#ifndef __DISABLEDEEP_ASA310

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ncbeta ( void * data)
/******************************************************************************/
/*
  Purpose:
    NCBETA computes the noncentral Beta CDF.
  Discussion:
    Three corrections needed to be made to the text of this routine.
    They are noted in the comments below.
    Two of these corrections were errors in transcription made when
    producing the online copy distributed by APSTAT.
    One error, an error of omission, occurred in the original printed
    copy of the routine, and was carried over into the online copy.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 February 2008
  Author:
    Original FORTRAN77 version by R Chattamvelli, R Shanmugam.
    C version by John Burkardt.
  Reference:
    R Chattamvelli, R Shanmugam,
    Algorithm AS 310:
    Computing the Non-central Beta Distribution Function,
    Applied Statistics,
    Volume 46, Number 1, 1997, pages 146-156.
  Parameters:
    Input, double A, B, the shape parameters.
    0 <= A, 0 <= B.
    Input, double LAMBDA, the noncentrality parameter.
    0 <= LAMBDA.
    Input, double X, the value at which the CDF is desired.
    Input, double ERRMAX, the precision tolerance.
    Output, int *IFAULT, error flag.
    0, no error occurred.
    1, X is 0 or 1.
    2, X < 0 or 1 < X.
    3, A, B or LAMBDA is less than 0.
    Output, double NCBETA, the value of the noncentral Beta CDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp b = a_data[1];
	const register ityp lambda = a_data[2];
	const register ityp x = a_data[3];
	ityp errmax = a_data[4];
	
	ityp beta;
	ityp c;
	ityp ebd;
	ityp errbd;
	ityp ftemp;
	ityp fx;
	ityp gx;
	dim_typ i;
	dim_typ iter1;
	dim_typ iter2;
	dim_typ iterhi;
	dim_typ iterlo;
	dim_typ j;
	dim_typ m;
	ityp mr;
	ityp psum;
	ityp q;
	ityp r;
	ityp s;
	ityp s0;
	ityp s1;
	ityp sum;
	ityp t;
	ityp t0;
	ityp t1;
	ityp temp;
	dim_typ xj;

	/*
	Check parameters.
	*/
	if ( lambda <= 0.0 || a <= 0.00 || b <= 0.00 || x <= 0.00 || 1.00 <= x)
	{
		result = NCBETA_INVALIDRETURNVALUE;
		return &result;
	}

	c = 0.50 * lambda;
	xj = 0.00;
	/*
	AS 226 as it stands is sufficient in this situation.
	*/
	if ( lambda < 54.0 )
	{
		result = betanc ( x, a, b, lambda );
		return &result;
	}
	else
	{
		m = ( dim_typ ) ( c + 0.50 );
		mr = ( ityp ) ( m );
		iterlo = m - ( int ) ( 5.00 * sqrt ( mr ) );
		iterhi = m + ( int ) ( 5.00 * sqrt ( mr ) );
		t = - c + mr * log ( c ) - lgamma ( mr + 1.00 );
		q = exp ( t );
		r = q;
		psum = q;
		beta = lgamma ( a + mr ) + lgamma ( b ) - lgamma ( a + mr + b );

		s1 = ( a + mr ) * log ( x )+ b * log ( 1.0 - x ) - log ( a + mr ) - beta;
		gx = exp ( s1 );
		fx = gx;
		temp = betain ( x, a + mr, b, beta );
		ftemp = temp;
		xj = xj + 1.0;
		/*
		The online copy of AS 310 has "SUM = Q - TEMP" which is incorrect.
		*/
		sum = q * temp;
		iter1 = m;
		/*
		The first set of iterations starts from M and goes downwards
		*/
		for ( ; ; )
		{
			if ( iter1 < iterlo || q < errmax)
				break;
			/*
			The online copy of AS 310 has "Q = Q - ITER1 / C" which is incorrect.
			*/
			q *= iter1 / c;
			xj += 1.00;
			gx = ( a + iter1 ) / ( x * ( a + b + iter1 - 1.0 ) ) * gx;
			-- iter1;
			temp += gx;
			psum +=  q;
			sum += q * temp;
		}

		t0 = lgamma ( a + b ) - lgamma ( a + 1.00 ) - lgamma ( b );
		s0 = a * log ( x ) + b * log ( 1.00 - x );
		/*
		Both the online copy of AS 310 and the text printed in the reference
		did not initialize the variable S to zero, which is incorrect.
		JVB, 12 January 2008.
		*/
		s = 0.00;
		for ( i = 1; i <= iter1; ++i )
		{
			j = i - 1;
			s += exp ( t0 + s0 + j * log ( x ) );
			t1 = log ( a + b + j ) - log ( a + 1.00 + j ) + t0;
			t0 = t1;
		}
		/*
		Compute the first part of error bound.
		*/
		errbd = ( 1.00 - gammad ( c, ( ityp ) ( iter1 ) ) ) * ( temp + s );
		q = r;
		temp = ftemp;
		gx = fx;
		iter2 = m;

		for ( ; ; )
		{
			ebd = errbd + ( 1.0 - psum ) * temp;

			if ( ebd < errmax || iterhi <= iter2 )
				break;
			++ iter2;
			xj += 1.00;
			q *= c / iter2;
			psum += q;
			temp -= gx;
			gx = x * ( a + b + iter2 - 1.0 ) / ( a + iter2 ) * gx;
			sum += q * temp;
		}
	}
	
	result = sum;
	return &result;
}

#endif
