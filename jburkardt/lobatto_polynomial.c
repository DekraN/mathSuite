#ifndef __DISABLEDEEP_LOBATTOPOLYNOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lobatto_polynomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOBATTO_POLYNOMIAL_VALUE evaluates the completed Lobatto polynomials Lo(n,x).
  Discussion:
    L(N,X) = ( 1 - X^2 ) * P'(N,X)
           = N * ( P(N-1,X) - X * P(N,X) )
    The Lobatto polynomials are 0 at -1 and +1.
   (1-x^2) * 1
   (1-x^2) * 3X
   (1-x^2) * ( -3 + 15x^2 ) / 2
   (1-x^2) * ( -60x + 140x^3 ) / 8
   (1-x^2) * ( -15 - 210x^2 + 315x^4 ) / 8
   (1-x^2) * ( 210x - 1260x^3 + 1386x^5 ) / 16
   (1-x^2) * ( -35 + 945x^2 - 3465x^4 + 3003x^6 ) / 16
   (1-x^2) * ( -2520x + 27720x^3 - 72072x^5 + 51480x^7 ) / 128
   (1-x^2) * ( 315 - 13860x^2 + 90090x^4 - 180180x^6 + 109395x^8 ) / 128
   (1-x^2) * ( 6930x - 120120x^3 + 540540x^5 - 875160x^7 + 461890x^9 ) / 256
    Mathematica: (replacing "n" by desired index):
      Expand [ ( 1-x^2) * D [ LegendreP[n,x], {x} ] ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 November 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Larry Andrews,
    Special Functions of Mathematics for Engineers,
    Second Edition,
    Oxford University Press, 1998,
    ISBN: 0819426164,
    LC: QA351.A75.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.
    Input, double X[M], the evaluation points.
    Output, double LOBATTO_POLYNOMIAL_VALUE[M*N], the values of the
    completed Lobatto polynomials of order 1 through N at the point X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	ityp l[n];
	ityp p[n+1];

	if(1 <= n)
	{
		l[0] = 1.00 - x*x;
		if(2 <= n)
		{
			p[0] = 1.00;
			p[1] = x;

			dim_typ i;

			for ( i=1; i< n; ++i)
				p[i+1] = ((ityp)((i<<1)+1)*x*p[i]-(ityp)(i)*p[i-1])/(ityp)(i+1);

			#pragma omp parallel for
			for (i=1; i < n; ++i)
				l[i] = (ityp)(i+1)*(p[i]-x*p[i+1]);
		}
	}

	result = l[n-1];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lobatto_polynomial_derivative ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOBATTO_POLYNOMIAL_DERIVATIVE: derivative of completed Lobatto polynomial.
  Discussion:
    L(N,X)  =  N * ( P(N-1,X) - X * P(N,X) )
    L'(N,X) =  N * ( P'(N-1,X) - P(N,X) - X * P'(N,X) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 November 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Larry Andrews,
    Special Functions of Mathematics for Engineers,
    Second Edition,
    Oxford University Press, 1998,
    ISBN: 0819426164,
    LC: QA351.A75.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.
    Input, double X[M], the evaluation points.
    Output, double LOBATTO_POLYNOMIAL_DERIVATIVE[M*N], the derivative of
    the completed Lobatto polynomials of order 1 through N at the point X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	ityp lp[n];
	ityp p[n+1];
	ityp pp[n+1];

	if ( 1 <= n )
	{
		lp[0] = -2.00 * x;

		if(2 <= n)
		{
			p[0] = 1.00;
			p[1]= x;

			dim_typ i;

			for(i = 1; i < n; ++i)
				p[i+1] = ((ityp)((i<<1)+1)*x*p[i]-(ityp)(i)*p[i-1])/(ityp)(i+1);

			pp[0] = 0.00;
			pp[1] = 1.00;

			for(i = 1; i < n; ++i)
				pp[i+1] = ((ityp)((i<<1)+1)*(p[i]+x*pp[i])-(ityp)(i)*pp[i-1])/(ityp)(i+1);

			#pragma omp parallel for
			for(i = 1; i < n; ++i)
				lp[i] = (ityp)(i+1)*(pp[i]-p[i+1]-x*pp[i+1]);
		}
	}

	result = lp[n-1];
  	return &result;
}

#endif
