#ifndef __DISABLEDEEP_BRENT

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _glomin ( void * data)
/******************************************************************************/
/*
  Purpose:
    GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
  Discussion:
    This function assumes that F(X) is twice continuously differentiable
    over [A,B] and that F''(X) <= M for all X in [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2008
  Author:
    Original FORTRAN77 version by Richard Brent.
    C version by John Burkardt.
  Reference:
    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.
  Parameters:
    Input, double A, B, the endpoints of the interval.
    It must be the case that A < B.
    Input, double C, an initial guess for the global
    minimizer.  If no good guess is known, C = A or B is acceptable.
    Input, double M, the bound on the second derivative.
    Input, double MACHEP, an estimate for the relative machine
    precision.
    Input, double E, a positive tolerance, a bound for the
    absolute error in the evaluation of F(X) for any X in [A,B].
    Input, double T, a positive error tolerance.
    Input, double F ( double x ), a user-supplied
    function whose global minimum is being sought.
    Output, double *X, the estimated value of the abscissa
    for which F attains its global minimum value in [A,B].
    Output, double GLOMIN, the value F(X).
*/
{
	static ityp result = MAX_VAL;
	
	const _7itfitpit * const s_data = data;
	ityp a = s_data->a0;
	ityp b = s_data->a1;
	ityp c = s_data->a2;
	ityp m = s_data->a3;
	ityp machep = s_data->a4;
	ityp e = s_data->a5;
	ityp t = s_data->a6;
	ityp (* f)(ityp) = s_data->a7;
	ityp * x = s_data->a8;
	
	ityp a0;
	ityp a2;
	ityp a3;
	ityp d0;
	ityp d1;
	ityp d2;
	ityp h;
	dim_typ k;
	ityp m2;
	ityp p;
	ityp q;
	ityp qs;
	ityp r;
	ityp s;
	ityp sc;
	ityp y;
	ityp y0;
	ityp y1;
	ityp y2;
	ityp y3;
	ityp yb;
	ityp z0;
	ityp z1;
	ityp z2;

	a0 = b;
	*x = a0;
	a2 = a;
	y0 = f ( b );
	yb = y0;
	y2 = f ( a );
	y = y2;


	if ( y0 < y )
		y = y0;
	else
		*x = a;

	if ( m <= 0.00 || b <= a )
	{
		result = y;
		return &result;
	}

	m2 = 0.50 * ( 1.00 + 16.00 * machep ) * m;

	sc = (c <= a || b <= c) ? 0.50 * ( a + b ) : c;

	y1 = f ( sc );
	k = 3;
	d0 = a2 - sc;
	h = 9.00 / 11.00;

	if ( y1 < y )
	{
		*x = sc;
		y = y1;
	}

	for ( ; ; )
	{
		d1 = a2 - a0;
		d2 = sc - a0;
		z2 = b - a2;
		z0 = y2 - y1;
		z1 = y2 - y0;
		r = d1 * d1 * z0 - d0 * d0 * z1;
		p = r;
		qs = 2.00 * ( d0 * z1 - d1 * z0 );
		q = qs;

		if ( k < 1000000 || y2 <= y )
		{
			for ( ; ; )
			{
				if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < z2 * m2 * r * ( z2 * q - r ) )
				{
					a3 = a2 + r / q;
					y3 = f ( a3 );
					if ( y3 < y )
					{
						*x = a3;
						y = y3;
					}
				}
				k = ( ( 1611 * k ) % 1048576 );
				q = 1.00;
				r = ( b - a ) * 0.00001 * ( double ) ( k );

				if ( z2 <= r )
					break;
			}
		}
		else
		{
			k = (( 1611*k)%1048576);
			q = 1.00;
			r = ( b - a ) * 0.00001 * ( ityp ) ( k );

			while ( r < z2 )
			{
				if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) < z2 * m2 * r * ( z2 * q - r ) )
				{
					a3 = a2 + r / q;
					y3 = f ( a3 );

					if ( y3 < y )
					{
						*x = a3;
						y = y3;
					}
				}
				k = ( ( 1611 * k ) % 1048576 );
				q = 1.0;
				r = ( b - a ) * 0.00001 * ( ityp ) ( k );
			}
		}

		r = m2 * d0 * d1 * d2;
		s = sqrt ( ( ( y2 - y ) + t ) / m2 );
		h = 0.50 * ( 1.00 + h );
		p = h * ( p + 2.00 * r * s );
		q = q + 0.50 * qs;
		r = - 0.50 * ( d0 + ( z0 + 2.01 * e ) / ( d0 * m2 ) );

		r = a2 + ( r < s || d0 < 0.00 ) ? s:r;
		a3 = (0.0 < p * q) ? a2+p/q : r;

		for ( ; ; )
		{
			a3 = MAX ( a3, r );

			if ( b <= a3 )
			{
				a3 = b;
				y3 = yb;
			}
			else
				y3 = f ( a3 );

			if ( y3 < y )
			{
				*x = a3;
				y = y3;
			}

			d0 = a3 - a2;

			if ( a3 <= r )
				break;

			p = 2.00 * ( y2 - y3 ) / ( m * d0 );

			if ( ( 1.00 + 9.00 * machep ) * d0 <= fabs ( p ) )
				break;

			if ( 0.5 * m2 * ( d0 * d0 + p * p ) <= ( y2 - y ) + ( y3 - y ) + 2.0 * t )
				break;
			a3 = 0.5 * ( a2 + a3 );
			h = 0.9 * h;

		}

		if ( b <= a3 )
			break;

		a0 = sc;
		sc = a2;
		a2 = a3;
		y0 = y1;
		y1 = y2;
		y2 = y3;
	}
	
	result = y;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _local_min ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
  Discussion:
    The method used is a combination of golden section search and
    successive parabolic interpolation.  Convergence is never much slower
    than that for a Fibonacci search.  If F has a continuous second
    derivative which is positive at the minimum (which is not at A or
    B), then convergence is superlinear, and usually of the order of
    about 1.324....
    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
    F is never evaluated at two points closer than TOL.
    If F is a unimodal function and the computed values of F are always
    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
    LOCAL_MIN approximates the abscissa of the global minimum of F on the
    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
    If F is not unimodal, then LOCAL_MIN may approximate a local, but
    perhaps non-global, minimum to the same accuracy.
    Thanks to Jonathan Eggleston for pointing out a correction to the
    golden section step, 01 July 2013.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2013
  Author:
    Orignal FORTRAN77 version by Richard Brent.
    C version by John Burkardt.
  Reference:
    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.
  Parameters:
    Input, double A, B, the endpoints of the interval.
    Input, double EPS, a positive relative error tolerance.
    EPS should be no smaller than twice the relative machine precision,
    and preferably not much less than the square root of the relative
    machine precision.
    Input, double T, a positive absolute error tolerance.
    Input, double F ( double x ), a user-supplied
    function whose local minimum is being sought.
    Output, double *X, the estimated value of an abscissa
    for which F attains a local minimum value in [A,B].
    Output, double LOCAL_MIN, the value F(X).
*/
{
	static ityp result = MAX_VAL;
	
	const _4pitfitpit * const s_data = data;
	ityp a = s_data->a0;
	ityp b = s_data->a1;
	ityp eps = s_data->a2;
	ityp t = s_data->a3;
	ityp (* f)(ityp) = s_data->a4;
	ityp * x = s_data->a5;
	
	ityp c;
	ityp d;
	ityp e;
	ityp fu;
	ityp fv;
	ityp fw;
	ityp fx;
	ityp m;
	ityp p;
	ityp q;
	ityp r;
	ityp sa;
	ityp sb;
	ityp t2;
	ityp tol;
	ityp u;
	ityp v;
	ityp w;
	/*
	C is the square of the inverse of the golden ratio.
	*/
	c = 0.50 * ( 3.00 - sqrt ( 5.00 ) );

	sa = a;
	sb = b;
	*x = sa + c * ( b - a );
	w = *x;
	v = w;
	e = 0.0;
	fx = f ( *x );
	fw = fx;
	fv = fw;

	for ( ; ; )
	{
		m = 0.50 * ( sa + sb ) ;
		tol = eps * fabs ( *x ) + t;
		t2 = 2.00 * tol;
		/*
		Check the stopping criterion.
		*/
		if ( fabs ( *x - m ) <= t2 - 0.50 * ( sb - sa ) )
			break;
		/*
		Fit a parabola.
		*/
		r = 0.00;
		q = r;
		p = q;

		if ( tol < fabs ( e ) )
		{
			r = ( *x - w ) * ( fx - fv );
			q = ( *x - v ) * ( fx - fw );
			p = ( *x - v ) * q - ( *x - w ) * r;
			q = 2.00 * ( q - r );
			if ( 0.00 < q )
				p *= -1;
			q = fabs ( q );
			r = e;
			e = d;
		}

		if ( fabs ( p ) < fabs ( 0.5 * q * r ) && q * ( sa - *x ) < p &&  p < q * ( sb - *x ) )
		{
			/*
			Take the parabolic interpolation step.
			*/
			d = p / q;
			u = *x + d;
			/*
			F must not be evaluated too close to A or B.
			*/
			if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
				d = tol *(-1+((*x < m)<<1));
		}
		/*
		A golden-section step.
		*/
		else
		{
			e = sb - (*x < m ? sb:sa) - *x;
			d = c * e;
		}
		/*
		F must not be evaluated too close to X.
		*/
		u = *x + tol <= fabs ( d ) ? d : tol*(-1+((0.00 < d)<<1));
		fu = f ( u );
		/*
		Update A, B, V, W, and X.
		*/
		if ( fu <= fx )
		{
			if ( u < *x )
				sb = *x;
			else
				sa = *x;
			v = w;
			fv = fw;
			w = *x;
			fw = fx;
			*x = u;
			fx = fu;
		}
		else
		{
			if ( u < *x )
				sa = u;
			else
				sb = u;

			if ( fu <= fw || w == *x )
			{
				v = w;
				fv = fw;
				w = u;
				fw = fu;
			}
			else if ( fu <= fv || v == *x || v == w )
			{
				v = u;
				fv = fu;
			}
		}
	}
	
	result = fx;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT inline void *   _r8_epsilon(void * data)
/******************************************************************************/
/*
  Purpose:
    R8_EPSILON returns the R8 round off unit.
  Discussion:
    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 September 2012
  Author:
    John Burkardt
  Parameters:
    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
	static ityp result = 2.220446049250313E-016;;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _zero ( void * data)
/******************************************************************************/
/*
  Purpose:
    ZERO seeks the root of a function F(X) in an interval [A,B].
  Discussion:
    The interval [A,B] must be a change of sign interval for F.
    That is, F(A) and F(B) must be of opposite signs.  Then
    assuming that F is continuous implies the existence of at least
    one value C between A and B for which F(C) = 0.
    The location of the zero is determined to within an accuracy
    of 6 * MACHEPS * abs ( C ) + 2 * T.
    Thanks to Thomas Secretin for pointing out a transcription error in the
    setting of the value of P, 11 February 2013.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2013
  Author:
    Original FORTRAN77 version by Richard Brent.
    C version by John Burkardt.
  Reference:
    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.
  Parameters:
    Input, double A, B, the endpoints of the change of sign interval.
    Input, double MACHEP, an estimate for the relative machine
    precision.
    Input, double T, a positive error tolerance.
    Input, double F ( double x ), a user-supplied function whose zero
    is being sought.
    Output, double ZERO, the estimated value of a zero of
    the function F.
*/
{
	static ityp result = MAX_VAL;
	
	const _4pitfit * const s_data = data;
	ityp a = s_data->a0;
	ityp b = s_data->a1;
	ityp machep = s_data->a2;
	ityp t = s_data->a3;
	ityp (* f)(ityp) = s_data->a4;
	
	ityp c;
	ityp d;
	ityp e;
	ityp fa;
	ityp fb;
	ityp fc;
	ityp m;
	ityp p;
	ityp q;
	ityp r;
	ityp s;
	ityp sa;
	ityp sb;
	ityp tol;
	/*
	Make local copies of A and B.
	*/
	sa = a;
	sb = b;
	fa = f ( sa );
	fb = f ( sb );

	c = sa;
	fc = fa;
	e = sb - sa;
	d = e;

	for ( ; ; )
	{
		if ( fabs ( fc ) < fabs ( fb ) )
		{
			sa = sb;
			sb = c;
			c = sa;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol = 2.00 * machep * fabs ( sb ) + t;
		m = 0.50 * ( c - sb );

		if ( fabs ( m ) <= tol || fb == 0.0 )
			break;

		if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
		{
			e = m;
			d = e;
		}
		else
		{
			s = fb / fa;

			if ( sa == c )
			{
				p = 2.00 * m * s;
				q = 1.00 - s;
			}
			else
			{
				q = fa / fc;
				r = fb / fc;
				p = s * ( 2.00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.00 ) );
				q = ( q - 1.00 ) * ( r - 1.00 ) * ( s - 1.00 );
			}

			if ( 0.00 < p )
				q *= -1;
			else
				p *= -1;

			s = e;
			e = d;

			if ( 2.00 * p < 3.00 * m * q - fabs ( tol * q ) && p < fabs ( 0.50 * s * q ) )
				d = p / q;
			else
			{
				e = m;
				d = e;
			}
		}
		sa = sb;
		fa = fb;
		sb += tol < fabs ( d ) ? d : tol*(-1+((0.00 < m)<<1));
		fb = f ( sb );

		if ( ( 0.00 < fb && 0.00 < fc ) || ( fb <= 0.00 && fc <= 0.00 ) )
		{
			c = sa;
			fc = fa;
			e = sb - sa;
			d = e;
		}
	}
	
	result = sb;
	return &result;
}

#endif
