#ifndef __DISABLEDEEP_CHEBYSHEVSERIES

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _echebser0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ECHEBSER0 evaluates a Chebyshev series.
  Discussion:
    This function implements a modification and extension of
    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
    gives an example for treating the first derivative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double ECHEBSER0, the value of the Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data; 
	
	const register dim_typ nc = s_data->a0;
	ityp * coef =  s_data->a1;
	const register ityp x = s_data->a2;
	
    ityp b0;
    ityp b1 = 0.0;
    ityp b2 = 0.0;
    dim_typ i;
    ityp x2;

    b0 = coef[nc-1];
    x2 = 2.0 * x;
    for ( i = nc - 2; 0 <= i; --i )
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + x2 * b1;
    }

	result = 0.5 * ( b0 - b2 );
    return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _echebser1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ECHEBSER1 evaluates a Chebyshev series and first derivative.
  Discussion:
    This function implements a modification and extension of
    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
    gives an example for treating the first derivative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double ECHEBSER1, the value of the Chebyshev series at X.
    Output, double *Y1, the value of the 1st derivative of the
    Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit2pit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	const register ityp x = s_data->a1;
	ityp * coef =  s_data->a2;
	ityp * y1 = s_data->a3;
	
    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    ityp c0;
    ityp c1 = 0.00;
    ityp c2 = 0.00;
    dim_typ i;
    ityp value;
    ityp x2;

    b0 = coef[nc-1];
    c0 = coef[nc-1];

    x2 = 2.00 * x;

    for ( i = nc - 2; 0 <= i; --i )
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + x2 * b1;

        if ( 0 < i )
        {
            c2 = c1;
            c1 = c0;
            c0 = b0 - c2 + x2 * c1;
        }
    }
    *y1 = c0 - c2;
    result = 0.50 * ( b0 - b2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _echebser2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ECHEBSER2 evaluates a Chebyshev series and two derivatives.
  Discussion:
    This function implements a modification and extension of
    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
    gives an example for treating the first derivative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double ECHEBSER2, the value of the Chebyshev series at X.
    Output, double *Y1, the value of the 1st derivative of the
    Chebyshev series at X.
    Output, double *Y2, the value of the 2nd derivative of the
    Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit3pit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	const register ityp x = s_data->a1;
	ityp * y1 = s_data->a2;
	ityp * y2 = s_data->a3;
	ityp * coef =  s_data->a4;
	
    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    ityp c0;
    ityp c1 = 0.00;
    ityp c2 = 0.00;
    ityp d0;
    ityp d1 = 0.00;
    ityp d2 = 0.00;
    dim_typ i;
    ityp value;
    ityp x2;

    b0 = coef[nc-1];
    c0 = coef[nc-1];
    d0 = coef[nc-1];

    x2 = 2.00 * x;

    for ( i = nc - 2; 0 <= i; --i)
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + x2 * b1;

        if ( 0 < i )
        {
            c2 = c1;
            c1 = c0;
            c0 = b0 - c2 + x2 * c1;
        }
        if ( 1 < i )
        {
            d2 = d1;
            d1 = d0;
            d0 = c0 - d2 + x2 * d1;
        }
    }

    *y2 = ( d0 - d2 ) * 4.0;
    *y1 = c0 - c2;
    
    result = 0.5 * ( b0 - b2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _echebser3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ECHEBSER3 evaluates a Chebyshev series and three derivatives.
  Discussion:
    This function implements a modification and extension of
    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
    gives an example for treating the first derivative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double ECHEBSER3, the value of the Chebyshev series at X.
    Output, double *Y1, the value of the 1st derivative of the
    Chebyshev series at X.
    Output, double *Y2, the value of the 2nd derivative of the
    Chebyshev series at X.
    Output, double *Y3, the value of the 3rd derivative of the
    Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pitit2pit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	ityp * y1 = s_data->a1;
	ityp * y2 = s_data->a2;
	const register ityp x = s_data->a3;
	ityp * y3 = s_data->a4;
	ityp * coef =  s_data->a5;
	
    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    ityp c0;
    ityp c1 = 0.00;
    ityp c2 = 0.00;
    ityp d0;
    ityp d1 = 0.00;
    ityp d2 = 0.00;
    ityp e0;
    ityp e1 = 0.00;
    ityp e2 = 0.00;
    dim_typ i;
    ityp value;
    ityp x2;

    b0 = coef[nc-1];
    c0 = coef[nc-1];
    d0 = coef[nc-1];
    e0 = coef[nc-1];

    x2 = 2.00 * x;

    for ( i = nc - 2; 0 <= i; --i)
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + x2 * b1;

        if ( 0 < i )
        {
            c2 = c1;
            c1 = c0;
            c0 = b0 - c2 + x2 * c1;
        }
        if ( 1 < i )
        {
            d2 = d1;
            d1 = d0;
            d0 = c0 - d2 + x2 * d1;
        }
        if ( 2 < i )
        {
            e2 = e1;
            e1 = e0;
            e0 = d0 - e2 + x2 * e1;
        }
    }

    *y3 = ( e0 - e2 ) * 24.00;
    *y2 = ( d0 - d2 ) * 4.00;
    *y1 = c0 - c2;
    
    result = 0.50 * ( b0 - b2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _echebser4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ECHEBSER4 evaluates a Chebyshev series and four derivatives.
  Discussion:
    This function implements a modification and extension of
    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
    gives an example for treating the first derivative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double ECHEBSER3, the value of the Chebyshev series at X.
    Output, double *Y1, *Y2, *Y3, *Y4, the value of the first four
    derivatives of the Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const itpitdt4pit * const s_data = data;
	const register ityp x = s_data->a0;
	ityp * coef =  s_data->a1;
	const register dim_typ nc = s_data->a2;
	ityp * y1 = s_data->a3;
	ityp * y2 = s_data->a4;
	ityp * y3 = s_data->a5;
	ityp * y4 = s_data->a6;

    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    ityp c0;
    ityp c1 = 0.00;
    ityp c2 = 0.00;
    ityp d0;
    ityp d1 = 0.00;
    ityp d2 = 0.00;
    ityp e0;
    ityp e1 = 0.00;
    ityp e2 = 0.00;
    ityp f0;
    ityp f1 = 0.00;
    ityp f2 = 0.00;
    dim_typ i;
    ityp y0;

    b0 = coef[nc-1];
    c0 = coef[nc-1];
    d0 = coef[nc-1];
    e0 = coef[nc-1];
    f0 = coef[nc-1];

    for ( i = nc - 2; 0 <= i; --i )
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + 2.00 * x * b1;

        if ( 0 < i )
        {
            c2 = c1;
            c1 = c0;
            c0 = b0 - c2 + 2.00 * x * c1;
        }
        if ( 1 < i )
        {
            d2 = d1;
            d1 = d0;
            d0 = c0 - d2 + 2.00 * x * d1;
        }
        if ( 2 < i )
        {
            e2 = e1;
            e1 = e0;
            e0 = d0 - e2 + 2.00 * x * e1;
        }
        if ( 3 < i )
        {
            f2 = f1;
            f1 = f0;
            f0 = e0 - f2 + 2.00 * x * f1;
        }
    }
    *y1 =   c0 - c2;
    *y2 = ( d0 - d2 ) *  2.00 * 2.00;
    *y3 = ( e0 - e2 ) *  6.00 * 4.00;
    *y4 = ( f0 - f2 ) * 24.00 * 8.00;

	result = ( b0 - b2 )        / 2.00;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _evenchebser0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    EVENCHEBSER0 evaluates an even Chebyshev series.
  Discussion:
    This function implements Clenshaw's modification of his
    algorithm for even series.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double EVENCHEBSER0, the value of the Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	ityp * coef =  s_data->a1;
	const register ityp x = s_data->a2;
	
    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    dim_typ i;
    ityp value;
    ityp x2;

    b0 = coef[nc-1];

    x2 = 4.00 * x * x - 2.00;

    for ( i = nc - 2; 0 <= i; --i )
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + x2 * b1;
    }

	result = 0.50 * ( b0 - b2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _evenchebser1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    EVENCHEBSER1 evaluates an even Chebyshev series and first derivative.
  Discussion:
    This function implements a modification and extension of
    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
    gives an example for treating the first derivative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double EVENCHEBSER1, the value of the Chebyshev series at X.
    Output, double *Y1, the value of the 1st derivative of the
    Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit2pit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	const register ityp x = s_data->a1;
	ityp * coef =  s_data->a2;
	ityp * y1 = s_data->a3;
	
    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    ityp c0;
    ityp c1 = 0.00;
    ityp c2 = 0.00;
    dim_typ i;
    ityp value;
    ityp x2;

    b0 = coef[nc-1];
    c0 = coef[nc-1];

    x2 = 4.00 * x * x - 2.00;

    for ( i = nc - 2; 0 <= i; --i )
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + x2 * b1;
        if ( 0 < i )
        {
            c2 = c1;
            c1 = c0;
            c0 = b0 - c2 + x2 * c1;
        }
    }

    *y1 = ( c0 - c2 ) * 4.00 * x;
    
    result = 0.50 * ( b0 - b2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _evenchebser2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    EVENCHEBSER2 evaluates an even Chebyshev series and first two derivatives
  Discussion:
    This function implements a modification and extension of
    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
    gives an example for treating the first derivative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double EVENCHEBSER2, the value of the Chebyshev series at X.
    Output, double *Y1, the value of the 1st derivative of the
    Chebyshev series at X.
    Output, double *Y2, the value of the 2nd derivative of the
    Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit3pit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	const register ityp x = s_data->a1;
	ityp * y1 = s_data->a2;
	ityp * y2 = s_data->a3;
	ityp * coef =  s_data->a4;

    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    ityp c0;
    ityp c1 = 0.00;
    ityp c2 = 0.00;
    ityp d0;
    ityp d1 = 0.00;
    ityp d2 = 0.00;
    dim_typ i;
    ityp value;
    ityp x2;

    b0 = coef[nc-1];
    c0 = coef[nc-1];
    d0 = coef[nc-1];

    x2 = 4.00 * x * x - 2.00;

    for ( i = nc - 2; 0 <= i; --i)
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + x2 * b1;
        if ( 0 < i )
        {
            c2 = c1;
            c1 = c0;
            c0 = b0 - c2 + x2 * c1;
        }
        if ( 1 < i )
        {
            d2 = d1;
            d1 = d0;
            d0 = c0 - d2 + x2 * d1;
        }
    }

    *y2 = ( d0 - d2 ) * 64.00 * x * x + ( c0 - c2 ) * 4.00;
    *y1 = ( c0 - c2 ) * 4.00 * x;
    
    result = 0.5 * ( b0 - b2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _oddchebser0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ODDCHEBSER0 evaluates an odd Chebyshev series.
  Discussion:
    This function implements Clenshaw's modification of  his algorithm
    for odd series.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double ODDCHEBSER0, the value of the Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	ityp * coef =  s_data->a1;
	const register ityp x = s_data->a2;
	
    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    dim_typ i;
    ityp value;
    ityp x2;

    b0 = coef[nc-1];

    x2 = 4.00* x * x - 2.00;

    for ( i = nc - 2; 0 <= i; --i)
    {
        b2 = b1;
        b1 = b0;
        b0 = coef[i] - b2 + x2 * b1;
    }
    
    result = ( b0 - b1 ) * x;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _oddchebser1 ( void * data) 
/******************************************************************************/
/*
  Purpose:
    ODDCHEBSER1 evaluates an odd Chebyshev series and the first derivative.
  Discussion:
    This function implements a modification and extension of
    Clenshaw's algorithm.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double ODDCHEBSER1, the value of the Chebyshev series at X.
    Output, double *Y1, the value of the 1st derivative of the
    Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit2pit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	const register ityp x = s_data->a1;
	ityp * coef =  s_data->a2;
	ityp * y1 = s_data->a3;

    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    ityp c0;
    ityp c1 = 0.00;
    ityp c2 = 0.00;
    ityp coefi;
    dim_typ i;
    ityp value;
    ityp x2;

    coefi = 2.00 * coef[nc-1];
    b0 = coefi;
    c0 = coefi;

    x2 = 4.00 * x * x - 2.00;

    for ( i = nc - 2; 0 <= i; --i )
    {
        b2 = b1;
        b1 = b0;
        coefi = 2.0 * coef[i] - coefi;
        b0 = coefi - b2 + x2 * b1;
        if ( 0 < i )
        {
            c2 = c1;
            c1 = c0;
            c0 = b0 - c2 + x2 * c1;
        }
    }
    *y1 = ( c0 - c2 ) * 4.00 * x * x + ( b0 - b2 ) * 0.50;
    
    result = ( b0 - b2 ) * 0.50 * x;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _oddchebser2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ODDCHEBSER2 evaluates an odd Chebyshev series and first two derivatives.
  Discussion:
    This function implements a modification and extension of
    Clenshaw's algorithm.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2014
  Author:
    Manfred Zimmer
  Reference:
    Charles Clenshaw,
    Mathematical Tables, Volume 5,
    Chebyshev series for mathematical functions,
    London, 1962.
    Gerhard Maess,
    Vorlesungen ueber Numerische Mathematik II, Analysis,
    Berlin, Akademie_Verlag, 1984-1988,
    ISBN: 978-3764318840,
    LC: QA297.M325.��
  Parameters:
    Input, double X, the evaluation point.
    -1 <= X <= +1.
    Input, double COEF[NC], the Chebyshev series.
    Input, int NC, the number of terms in the Chebyshev series.
    0 < NC.
    Output, double ODDCHEBSER2, the value of the Chebyshev series at X.
    Output, double *Y1, the value of the 1st derivative of the
    Chebyshev series at X.
    Output, double *Y2, the value of the 2nd derivative of the
    Chebyshev series at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit3pit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	const register ityp x = s_data->a1;
	ityp * y1 = s_data->a2;
	ityp * y2 = s_data->a3;
	ityp * coef =  s_data->a4;

    ityp b0;
    ityp b1 = 0.00;
    ityp b2 = 0.00;
    ityp c0;
    ityp c1 = 0.00;
    ityp c2 = 0.00;
    ityp d0;
    ityp d1 = 0.00;
    ityp d2 = 0.00;
    ityp coefi;
    dim_typ i;
    ityp value;
    ityp x2;

    coefi = 2.00 * coef[nc-1];
    b0 = coefi;
    c0 = coefi;
    d0 = coefi;

    x2 = 4.00 * x * x - 2.00;

    for ( i = nc - 2; 0 <= i; --i)
    {
        b2 = b1;
        b1 = b0;
        coefi = 2.0* coef[i] - coefi;
        b0 = coefi - b2 + x2 * b1;
        if ( 0 < i )
        {
            c2 = c1;
            c1 = c0;
            c0 = b0 - c2 + x2 * c1;
        }
        if ( 1 < i )
        {
            d2 = d1;
            d1 = d0;
            d0 = c0 - d2 + x2 * d1;
        }
    }

    x2 = x * x;
    *y1 = ( c0 - c2 ) * 4.00 * x2  + ( b0 - b2 ) * 0.50;
    *y2 = ( ( d0 - d2 ) * 64.00 * x2 + ( c0 - c2 ) * 12.00 ) * x;

	result = ( b0 - b2 ) * 0.50 * x;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _a_to_i4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    A_TO_I4 returns the index of an alphabetic character.
  Example:
    CH  A_TO_I4
    'A'   1
    'B'   2
    ...
    'Z'  26
    'a'  27
    'b'  28
    ...
    'z'  52
    '$'   0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2008
  Author:
    John Burkardt
  Parameters:
    Input, char CH, a character.
    Output, int A_TO_I4, is the alphabetic index of the character,
    between 1 and 26 if the character is a capital letter,
    between 27 and 52 if it is lower case, and 0 otherwise.
*/
{
	static int result = INT_MAX;
	
	char ch = *(char *) data;
	
    if ( 'A' <= ch && ch <= 'Z' )
    {
    	result = ( ( int ) ( ch - 'A' + 1 ) );
    	return &result;
    }
    else if ( 'a' <= ch && ch <= 'z' )
    {
    	result = ( ( int ) ( ch - 'a' + 26 + 1 ) );
        return &result;
    }
    
    result = 0;
    return &result;
}

#endif
