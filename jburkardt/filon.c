#ifndef __DISABLEDEEP_FILON

#include "../dutils.h"

#define FILON_INVALIDRETURNVALUE MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   filon_fun_cos (const register dim_typ n, ityp *f ( dim_typ n, ityp x[] ), ityp a, ityp b, ityp t )
/******************************************************************************/
/*
  Purpose:
    FILON_FUN_COS uses Filon's method on integrals with a cosine factor.
  Discussion:
    The integral to be approximated has the form:
      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
    where T is user specified.
    The function is interpolated over each subinterval by
    a parabolic arc.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Chase, Lloyd Fosdick,
    An Algorithm for Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 453-457.
    Stephen Chase, Lloyd Fosdick,
    Algorithm 353:
    Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 457-458.
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, int N, the number of data points.
    N must be odd, and greater than 1.
    Input, double *F ( int n, double x[] ), the function which evaluates the
    integrand.
    Input, double A, B, the limits of integration.
    Input, double T, the multiplier of the X argument of the cosine.
    Output, double FILON_FUN_COS, the approximate value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
    ityp alpha;
    ityp beta;
    ityp c2n;
    ityp c2nm1;
    ityp cost;
    ityp *ftab;
    ityp gamma;
    ityp h;
    dim_typ i;
    ityp sint;
    ityp theta;
    ityp value;
    ityp *x;

    if ( a == b )
    {
    	result = 0.00;
        return &result;
    }

    if ( n <= 1 || n%2)
    {
    	result = FILON_INVALIDRETURNVALUE;
        return &result;
    }

    /*
    Set the X values.
    */
    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        x[i] = ( ( ityp ) ( n - i - 1 ) * a   + ( ityp ) (     i     ) * b ) / ( ityp ) ( n     - 1 );


    h = ( b - a ) / ( ityp ) ( n - 1 );
    theta = t * h;
    sint = sin ( theta );
    cost = cos ( theta );

    if ( 6.00 * fabs ( theta ) <= 1.0000 )
    {
        alpha = 2.00 * pow ( theta, 3 ) /   45.0 - 2.00 * pow ( theta, 5 ) /  315.0 + 2.00 * pow ( theta, 7 ) / 4725.00;
        beta =  2.00                    /     3.00 + 2.00 * pow ( theta, 2 ) /    15.00 - 4.00 * pow ( theta, 4 ) /   105.00 + 2.0 * pow ( theta, 6 ) /   567.00 - 4.000 * pow ( theta, 8 ) / 22275.00;
        gamma = 4.00                    /      3.0 - 2.00 * pow ( theta, 2 ) /     15.00 +       pow ( theta, 4 ) /    210.00 -       pow ( theta, 6 ) /  11340.00;
    }
    else
    {
        alpha = ( pow ( theta, 2 ) + theta * sint * cost - 2.00 * sint * sint ) / pow ( theta, 3 );
        beta = ( 2.00 * theta + 2.00 * theta * cost * cost- 4.00 * sint * cost ) / pow ( theta, 3 );
        gamma = 4.00 * ( sint - theta * cost ) / pow ( theta, 3 );
    }
    /*
    Tabulate the function.
    */
    ftab = f ( n, x );

    c2n = 0.50 * ftab[0] * cos ( t * x[0] );
    for ( i = 2; i < n - 1; i += 2 )
        c2n += ftab[i] * cos ( t * x[i] );
    c2n = c2n + 0.50 * ftab[n-1] * cos ( t * x[n-1] );

    c2nm1 = 0.00;
    for ( i = 1; i <= n - 2; i += 2)
        c2nm1 += ftab[i] * cos ( t * x[i] );

    value = h * ( alpha * ( ftab[n-1] * sin ( t * x[n-1] )  - ftab[0]   * sin ( t * x[0] ) ) + beta * c2n + gamma * c2nm1 );

    free ( ftab );
    free ( x );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _filon_tab_cos ( void * data)
/******************************************************************************/
/*
  Purpose:
    FILON_TAB_COS uses Filon's method on integrals with a cosine factor.
  Discussion:
    The integral to be approximated has the form:
      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
    where T is user specified.
    The function is interpolated over each subinterval by
    a parabolic arc.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Chase, Lloyd Fosdick,
    An Algorithm for Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 453-457.
    Stephen Chase, Lloyd Fosdick,
    Algorithm 353:
    Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 457-458.
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, int N, the number of data points.
    N must be odd, and greater than 1.
    Input, double FTAB[N], contains the value of the function
    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(N-1).
    Input, double A, B, the limits of integration.
    Input, double T, the multiplier of the X argument of the cosine.
    Output, double FILON_TAB_COS, the approximate value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const dt3itpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp a = s_data->a1;
	ityp b = s_data->a2;
	ityp t = s_data->a3;
	ityp * ftab = s_data->a4;
	
    ityp alpha;
    ityp beta;
    ityp c2n;
    ityp c2nm1;
    ityp cost;
    ityp gamma;
    ityp h;
    dim_typ i;
    ityp sint;
    ityp theta;
    ityp value;
    ityp *x;

    if ( a == b )
    {
    	result = 0.00;
        return &result;
    }

    if ( n <= 1 || n%2)
    {
    	result = FILON_INVALIDRETURNVALUE;
        return &result;
    }

    /*
    Set the X values.
    */
    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        x[i] = ( ( ityp ) ( n - i - 1 ) * a+ ( ityp ) (     i     ) * b )/ ( ityp ) ( n     - 1 );

    h = ( b - a ) / ( ityp ) ( n - 1 );
    theta = t * h;
    sint = sin ( theta );
    cost = cos ( theta );

    if ( 6.0 * fabs ( theta ) <= 1.00 )
    {
        alpha = 2.00 * pow ( theta, 3 ) /   45.00- 2.00 * pow ( theta, 5 ) /  315.00+ 2.00 * pow ( theta, 7 ) / 4725.00;
        beta =  2.00                    /     3.00+ 2.00 * pow ( theta, 2 ) /    15.00- 4.00 * pow ( theta, 4 ) /   105.00+ 2.00 * pow ( theta, 6 ) /   567.00- 4.00 * pow ( theta, 8 ) / 22275.00;
        gamma = 4.00                    /      3.00- 2.00 * pow ( theta, 2 ) /     15.00+       pow ( theta, 4 ) /    210.00-       pow ( theta, 6 ) /  11340.00;
    }
    else
    {
        alpha = ( pow ( theta, 2 ) + theta * sint * cost - 2.00 * sint * sint )/ pow ( theta, 3 );
        beta = ( 2.00 * theta + 2.00 * theta * cost * cost- 4.00 * sint * cost ) / pow ( theta, 3 );
        gamma = 4.00 * ( sint - theta * cost ) / pow ( theta, 3 );
    }

    c2n = + 0.50 * ftab[0] * cos ( t * x[0] );
    for ( i = 2; i < n - 1; i += 2 )
        c2n += ftab[i] * cos ( t * x[i] );
    c2n = c2n + 0.50 * ftab[n-1] * cos ( t * x[n-1] );

    c2nm1 = 0.00;
    for ( i = 1; i <= n - 2; i += 2 )
        c2nm1 += ftab[i] * cos ( t * x[i] );

    value = h * (
    alpha * ( ftab[n-1] * sin ( t * x[n-1] )- ftab[0]   * sin ( t * x[0] ) )+ beta * c2n+ gamma * c2nm1 );

    free ( x );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   filon_fun_sin (const register dim_typ n, ityp *f ( dim_typ n, ityp x[] ), ityp a,ityp b, ityp t )
/******************************************************************************/
/*
  Purpose:
    FILON_FUN_SIN uses Filon's method on integrals with a sine factor.
  Discussion:
    The integral to be approximated has the form
      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
    where T is user specified.
    The function is interpolated over each subinterval by
    a parabolic arc.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Chase, Lloyd Fosdick,
    An Algorithm for Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 453-457.
    Stephen Chase, Lloyd Fosdick,
    Algorithm 353:
    Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 457-458.
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, int N, the number of data points,
    including the endpoints.  N must be odd, and greater than 1.
    Input, external F, the subroutine which evaluates the integrand,
    of the form subroutine F ( N, X, FX ).
    Input, double A, B, the limits of integration.
    Input, double T, multiplier of the X argument of the sine.
    Output, double FILON_FUN_SIN, the approximate value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
    ityp alpha;
    ityp beta;
    ityp cost;
    ityp *ftab;
    ityp gamma;
    ityp h;
    dim_typ i;
    ityp s2n;
    ityp s2nm1;
    ityp sint;
    ityp theta;
    ityp value;
    ityp *x;

    if ( a == b )
    {
    	result = 0.00;
        return &result;
    }
    
    if ( n <= 1 || n%2 )
    {
    	result = FILON_INVALIDRETURNVALUE;
        return &result;
    }
    /*
    Set the X values.
    */
    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        x[i] = ( ( ityp ) ( n - i - 1 ) * a+ ( ityp ) (     i     ) * b )/ ( ityp ) ( n     - 1 );

    h = ( b - a ) / ( ityp ) ( n - 1 );
    theta = t * h;
    sint = sin ( theta );
    cost = cos ( theta );

    if ( 6.00 * fabs ( theta ) <= 1.00 )
    {
        alpha = 2.00 * pow ( theta, 3 ) /   45.00- 2.00* pow ( theta, 5 ) /  315.00+ 2.00 * pow ( theta, 7 ) / 4725.00;
        beta =  2.00                    /     3.00+ 2.00 * pow ( theta, 2 ) /    15.00- 4.00 * pow ( theta, 4 ) /   105.00+ 2.00 * pow ( theta, 6 ) /   567.00- 4.00 * pow ( theta, 8 ) / 22275.00;
        gamma = 4.00                    /      3.00- 2.00 * pow ( theta, 2 ) /     15.00+       pow ( theta, 4 ) /    210.00-       pow ( theta, 6 ) /  11340.00;
    }
    else
    {
        alpha = ( pow ( theta, 2 ) + theta * sint * cost- 2.00 * sint * sint ) / pow ( theta, 3 );
        beta = ( 2.00 * theta + 2.00 * theta * cost * cost- 4.00 * sint * cost ) / pow ( theta, 3 );
        gamma = 4.00 * ( sint - theta * cost ) / pow ( theta, 3 );
    }
    /*
    Tabulate the function.
    */
    ftab = f ( n, x );

    s2n = + 0.50 * ftab[0] * sin ( t * x[0] );
    for ( i = 2; i < n - 1; i += 2 )
        s2n += ftab[i] * sin ( t * x[i] );
    s2n += 0.50 * ftab[n-1] * sin ( t * x[n-1] );

    s2nm1 = 0.00;
    for ( i = 1; i <= n - 2; i += 2 )
        s2nm1 += ftab[i] * sin ( t * x[i] );

    value = h * (
    alpha * ( ftab[0]   * cos ( t * x[0] )- ftab[n-1] * cos ( t * x[n-1] ) )+ beta * s2n+ gamma * s2nm1 );

    free ( ftab );
    free ( x );
    
    result = value;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _filon_tab_sin ( void * data)
/******************************************************************************/
/*
  Purpose:
    FILON_TAB_SIN uses Filon's method on integrals with a sine factor.
  Discussion:
    The integral to be approximated has the form
      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
    where T is user specified.
    The function is interpolated over each subinterval by
    a parabolic arc.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Chase, Lloyd Fosdick,
    An Algorithm for Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 453-457.
    Stephen Chase, Lloyd Fosdick,
    Algorithm 353:
    Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 457-458.
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, int N, the number of data points,
    including the endpoints.  N must be odd, and greater than 1.
    Input, double FTAB[N], contains the value of the function
    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(N-1).
    Input, double A, B, the limits of integration.
    Input, double T, multiplier of the X argument of the sine.
    Output, double FILON_TAB_SIN, the approximate value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const dt3itpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp a = s_data->a1;
	ityp b = s_data->a2;
	ityp t = s_data->a3;
	ityp * ftab = s_data->a4;
	
    ityp alpha;
    ityp beta;
    ityp cost;
    ityp gamma;
    ityp h;
    dim_typ i;
    ityp s2n;
    ityp s2nm1;
    ityp sint;
    ityp theta;
    ityp value;
    ityp *x;

    if ( a == b )
    {
    	result = 0.00;
    	return &result;
    }

    if ( n <= 1 || n%2)
    {
    	result = FILON_INVALIDRETURNVALUE; 
    	return &result;
    }
    /*
    Set the X values.
    */
    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    
    if ( a == b )
    {
    	result = 0.00;
    	return &result;
    }

    if ( n <= 1 || n%2)
    {
    	result = FILON_INVALIDRETURNVALUE;
        return &result;
    }
        
    for ( i = 0; i < n; ++i )
		x[i] = ( ( ityp ) ( n - i - 1 ) * a   + ( ityp ) (     i     ) * b ) / ( ityp ) ( n     - 1 );
	
    h = ( b - a ) / ( ityp ) ( n - 1 );
    theta = t * h;
    sint = sin ( theta );
    cost = cos ( theta );

    if ( 6.00 * fabs ( theta ) <= 1.00 )
    {
        alpha = 2.00 * pow ( theta, 3 ) /   45.00- 2.00 * pow ( theta, 5 ) /  315.00+ 2.00 * pow ( theta, 7 ) / 4725.00;
        beta =  2.00                    /     3.00+ 2.00 * pow ( theta, 2 ) /    15.00- 4.00 * pow ( theta, 4 ) /   105.00+ 2.00 * pow ( theta, 6 ) /   567.00- 4.00 * pow ( theta, 8 ) / 22275.00;
        gamma = 4.00                    /      3.00- 2.00 * pow ( theta, 2 ) /     15.0+       pow ( theta, 4 ) /    210.00-       pow ( theta, 6 ) /  11340.00;
    }
    else
    {
        alpha = ( pow ( theta, 2 ) + theta * sint * cost- 2.00 * sint * sint ) / pow ( theta, 3 );
        beta = ( 2.00 * theta + 2.00 * theta * cost * cost- 4.00 * sint * cost ) / pow ( theta, 3 );
        gamma = 4.00 * ( sint - theta * cost ) / pow ( theta, 3 );
    }

    s2n = + 0.50 * ftab[0] * sin ( t * x[0] );
    for ( i = 2; i < n - 1; i += 2)
        s2n += ftab[i] * sin ( t * x[i] );
    s2n += 0.50 * ftab[n-1] * sin ( t * x[n-1] );
    s2nm1 = 0.00;
    for ( i = 1; i <= n - 2; i += 2 )
        s2nm1 += ftab[i] * sin ( t * x[i] );

    value = h * (alpha * ( ftab[0]   * cos ( t * x[0] )- ftab[n-1] * cos ( t * x[n-1] ) )+ beta * s2n+ gamma * s2nm1 );
    free ( x );
    
    result = value;
    return &result;
}

#define FILON_INVALIDRETURNVALUE MAX_VAL

#endif
