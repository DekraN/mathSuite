#ifndef __DISABLEDEEP_EXACTNESS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hermite_exactness ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_EXACTNESS investigates exactness of Hermite quadrature.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points in the rule.
    Input, double X[N], the quadrature points.
    Input, double W[N], the quadrature weights.
    Input, int P_MAX, the maximum exponent.
    0 <= P_MAX.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ p_max = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	
    ityp e;
    dim_typ i, p;
    ityp q;
    ityp s;
    ityp *v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( p = 0; p <= p_max; p++ )
    {
        s = hermite_integral ( p );

        for ( i = 0; i < n; ++i )
            v[i] = pow ( x[i], p );

        q = r8vec_dot_product ( n, w, v );
        e = s == 0.00 ? fabs(q) : fabs ( q - s ) / fabs ( s );
    }

    free ( v );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hermite_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_INTEGRAL evaluates a monomial Hermite integral.
  Discussion:
    Integral ( -oo < x < oo ) x^p exp(-x^2) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, int P, the exponent of the monomial.
    0 <= P.
    Output, double HERMITE_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ p = *(dim_typ *) data;
	
	result = p%2 ? 0.00 : r8_factorial2 ( p - 1 ) * sqrt (M_PI) / pow ( 2.00, p / 2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laguerre_exactness ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_EXACTNESS investigates exactness of Laguerre quadrature.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points in the rule.
    Input, double X[N], the quadrature points.
    Input, double W[N], the quadrature weights.
    Input, int P_MAX, the maximum exponent.
    0 <= P_MAX.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ p_max = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	
    ityp e;
    dim_typ i, p;
    ityp q;
    ityp s;
    ityp *v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( p = 0; p <= p_max; ++p )
    {
        s = laguerre_integral ( p );

        for ( i = 0; i < n; ++i )
            v[i] = pow ( x[i], p );

        q = r8vec_dot_product ( n, w, v );
        e = s == 0.00 ? fabs(q) : fabs ( q - s ) / fabs ( s );
    }

    free ( v );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laguerre_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.
  Discussion:
    The integral being computed is
      integral ( 0 <= x < +oo ) x^p * exp ( -x ) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, int P, the exponent.
    0 <= P.
    Output, double LAGUERRE_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ p = *(dim_typ *) data;	
	
	result = r8_factorial ( p );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _legendre_exactness ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_EXACTNESS investigates exactness of Legendre quadrature.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points in the rule.
    Input, double X[N], the quadrature points.
    Input, double W[N], the quadrature weights.
    Input, int P_MAX, the maximum exponent.
    0 <= P_MAX.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ p_max = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	
    ityp e;
    dim_typ i, p;
    ityp q;
    ityp s;
    ityp *v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( p = 0; p <= p_max; ++p )
    {
        s = legendre_integral ( p );

        for ( i = 0; i < n; ++i )
            v[i] = pow ( x[i], p );

    q = r8vec_dot_product ( n, w, v );
    e = s == 0.00 ? fabs(q) : fabs ( q - s ) / fabs ( s );
    }

    free ( v );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _legendre_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.
  Discussion:
    Integral ( -1 <= x <= +1 ) x^p dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 February 2008
  Author:
    John Burkardt
  Parameters:
    Input, int P, the exponent.
    0 <= P.
    Outp__MATHSUITE __JBURKARDT  void  ut, double LEGENDRE_INTEGRAL, the value of the exact integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ p = *(dim_typ *) data;	
	
	result = p%2 ? 0.00 : 2.00 / ( ityp ) ( p + 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _legendre_monomial_quadrature ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 February 2008
  Author:
    John Burkardt
  Parameters:
    Input, int EXPON, the exponent.
    Input, int ORDER, the number of points in the rule.
    Input, double W[ORDER], the quadrature weights.
    Input, double X[ORDER], the quadrature points.
    Output, double LEGENDRE_MONOMIAL_QUADRATURE, the quadrature error.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dt2pit * const s_data = data;
	const register dim_typ expon = s_data->a0;
	const register dim_typ order = s_data->a1;
	ityp * w = s_data->a2;
	ityp * x = s_data->a3;
	
    /*
    Get the exact value of the integral.
    */
    ityp exact = legendre_integral ( expon );
    /*
    Evaluate the monomial at the quadrature points
    and compute the weighted sum.
    */
    ityp quad = 0.00;
    for (dim_typ i = 0; i < order; ++i )
        quad += w[i] * pow ( x[i], expon );
    /*
    Error:
    */
    
    result = exact == 0.00 ? fabs ( quad - exact ) : fabs ( ( quad - exact ) / exact );
    return &result;
}

#endif
