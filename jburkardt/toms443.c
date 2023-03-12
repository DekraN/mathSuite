#ifndef __DISABLEDEEP_TOMS443

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wew_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEW_A estimates Lambert's W function.
  Discussion:
    For a given X, this routine estimates the solution W of Lambert's
    equation:
      X = W * EXP ( W )
    This routine has higher accuracy than WEW_B.
  Modified:
    11 June 2014
  Reference:
    Fred Fritsch, R Shafer, W Crowley,
    Algorithm 443: Solution of the transcendental equation w e^w = x,
    Communications of the ACM,
    October 1973, Volume 16, Number 2, pages 123-124.
  Parameters:
    Input, double X, the argument of W(X)
    Output, double *EN, the last relative correction to W(X).
    Output, double WEW_A, the estimated value of W(X).
*/
{
	static ityp result = MAX_VAL;
	
	const itpit * const s_data = data;
	const register ityp x = s_data->a0;
	ityp * en = s_data->a1;
	
    const register ityp c1 = 4.00 / 3.00;
    const register ityp c2 = 7.00 / 3.00;
    const register ityp c3 = 5.00 / 6.00;
    const register ityp c4 = 2.00 / 3.00;
    ityp f;
    ityp temp;
    ityp temp2;
    ityp wn;
    ityp y;
    ityp zn;
    /*
    Initial guess.
    */
    f = log ( x );

    if ( x <= 6.46 )
    {
        wn = x * ( 1.00 + c1 * x ) / ( 1.00 + x * ( c2 + c3 * x ) );
        zn = f - wn - log ( wn );
    }
    else
    {
        wn = f;
        zn = - log ( wn );
    }
    /*
    Iteration 1.
    */
    temp = 1.00 + wn;
    y = 2.00 * temp * ( temp + c4 * zn ) - zn;
    wn *= ( 1.00 + zn * y / ( temp * ( y - zn ) ) );
    /*
    Iteration 2.
    */
    zn = f - wn - log ( wn );
    temp = 1.00 + wn;
    temp2 = temp + c4 * zn;
    *en = zn * temp2 / ( temp * temp2 - 0.5 * zn );
    
    result = wn * ( 1.00 + *en );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wew_b ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEW_B estimates Lambert's W function.
  Discussion:
    For a given X, this routine estimates the solution W of Lambert's
    equation:
      X = W * EXP ( W )
    This routine has lower accuracy than WEW_A.
  Modified:
    11 June 20149
  Reference:
    Fred Fritsch, R Shafer, W Crowley,
    Algorithm 443: Solution of the transcendental equation w e^w = x,
    Communications of the ACM,
    October 1973, Volume 16, Number 2, pages 123-124.
  Parameters:
    Input, double X, the argument of W(X)
    Output, double *EN, the last relative correction to W(X).
    Output, double WEW_B, the estimated value of W(X).
*/
{
	static ityp result = MAX_VAL;
	
	const itpit * const s_data = data;
	const register ityp x = s_data->a0;
	ityp * en = s_data->a1;
	
    const register ityp c1 = 4.00 / 3.00;
    const register ityp c2 = 7.00 / 3.00;
    const register ityp c3 = 5.00 / 6.00;
    const register ityp c4 = 2.00 / 3.00;
    ityp f;
    ityp temp;
    ityp wn;
    ityp y;
    ityp zn;
    /*
    Initial guess.
    */
    f = log ( x );

    wn = x <= 0.7385 ? x * ( 1.00 + c1 * x ) / ( 1.00 + x * ( c2 + c3 * x ) ) : f - 24.00 * ( ( f + 2.00 ) * f - 3.00 ) / ( ( 0.70 * f + 58.00 ) * f + 127.00 );
    /*
    Iteration 1.
    */
    zn = f - wn - log ( wn );
    temp = 1.00 + wn;
    y = 2.00 * temp * ( temp + c4 * zn ) - zn;
    *en = zn * y / ( temp * ( y - zn ) );
    wn *= ( 1.00 + *en );

	result = wn;
    return &result;
}

#endif
