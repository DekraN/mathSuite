#ifndef __DISABLEDEEP_LEBESGUE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _chebyshev1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV1 returns the Type 1 Chebyshev points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double CHEBYSHEV1[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp angle;
    dim_typ i;

    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        angle = M_PI * ( ityp ) ( (i<<1) + 1 ) / ( ityp ) ( n<<1 );
        x[i] = cos ( angle );
    }
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _chebyshev2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV2 returns the Type 2 Chebyshev points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double CHEBYSHEV2[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp angle;
    dim_typ i;

    ityp * x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
        x[0] = 0.00;
    else
    {
        for ( i = 0; i < n; ++i)
        {
            angle = M_PI * ( ityp ) ( n - i - 1 ) / ( ityp ) ( n - 1 );
            x[i] = cos ( angle );
        }
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _chebyshev3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV3 returns the Type 3 Chebyshev points.
  Discussion:
    Note that this point set is NOT symmetric in [-1,+1].
    It is sometimes augmented by the value -1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double CHEBYSHEV3[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp angle;
    dim_typ i;
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        angle = M_PI * ( ityp ) ( (n<<1) - (i<<1) - 1 ) / ( ityp ) ( (n<<1)         + 1 );
        x[i] = cos ( angle );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _chebyshev4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV4 returns the Type 4 Chebyshev points.
  Discussion:
    Note that this point set is NOT symmetric in [-1,+1].
    It is sometimes augmented by the value +1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double CHEBYSHEV4[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp angle;
    dim_typ i;
    ityp * x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i)
    {
        angle = M_PI * ( ityp ) ( (n<<1) - (i<<1) )/ ( ityp ) ( (n<<1) + 1 );
        x[i] = cos ( angle );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _equidistant1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    EQUIDISTANT1 returns the Type 1 Equidistant points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double EQUIDISTANT1[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        x[i] = ( ityp ) ( - n + 1 + (i<<1) ) / ( ityp ) ( n + 1 );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _equidistant2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    EQUIDISTANT2 returns the Type 2 Equidistant points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double EQUIDISTANT2[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp *x = ( double * ) malloc ( n * sizeof ( double ) );

    if ( n == 1 )
        x[0] = 0.00;
    else
        for ( i = 0; i < n; ++i )
            x[i] = ( ityp ) ( - n + 1 + (i<<1) ) / ( ityp ) ( n - 1 );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _equidistant3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    EQUIDISTANT3 returns the Type 3 Equidistant points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double EQUIDISTANT3[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        x[i] = ( ityp ) ( - n + 1 + (i<<1) ) / ( ityp ) ( n );
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _fejer1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEJER1 returns the Type 1 Fejer points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double FEJER1[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp theta;
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; i++ )
    {
        theta = M_PI * ( ityp ) ( (n<<1) - 1 - (i<<1) )/ ( ityp ) ( n<<1 );
        x[i] = cos ( theta );
    }
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _fejer2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEJER2 returns the Type 2 Fejer points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2018
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the number of points.
    Input, double FEJER2[N], the points.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp theta;
    ityp * x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i)
    {
        theta = M_PI * ( ityp ) ( n - i )/ ( ityp ) ( n + 1 );
        x[i] = cos ( theta );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lebesgue_constant ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEBESGUE_CONSTANT estimates the Lebesgue constant for a set of points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2014
  Author:
    John Burkardt.
  Parameters:
    Jean-Paul Berrut, Lloyd Trefethen,
    Barycentric Lagrange Interpolation,
    SIAM Review,
    Volume 46, Number 3, September 2004, pages 501-517.
  Parameters:
    Input, int N, the number of interpolation points.
    Input, double X[N], the interpolation points.
    Input, int NFUN, the number of evaluation points.
    Input, double XFUN[NFUN], the evaluation points.
    Output, double LEBESGUE_CONSTANT, an estimate of the Lebesgue constant
    for the points.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ nfun = s_data->a1;
	ityp * x = s_data->a2;
	ityp * xfun = s_data->a3;
	
    ityp *lfun;
    ityp lmax;
    lfun = lebesgue_function ( n, x, nfun, xfun );
    _MAX ( &lmax, nfun, lfun );
    free ( lfun );
    
    result = lmax;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lebesgue_function ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEBESGUE_FUNCTION evaluates the Lebesgue function for a set of points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2014
  Author:
    Original MATLAB version by Greg von Winckel.
    C version by John Burkardt.
  Parameters:
    Jean-Paul Berrut, Lloyd Trefethen,
    Barycentric Lagrange Interpolation,
    SIAM Review,
    Volume 46, Number 3, September 2004, pages 501-517.
  Parameters:
    Input, int N, the number of interpolation points.
    Input, double X[N], the interpolation points.
    Input, int NFUN, the number of evaluation points.
    Input, double XFUN[NFUN], the evaluation points.
    Output, double LEBESGUE_FUNCTION[NFUN], the Lebesgue function values.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ nfun = s_data->a1;
	ityp * x = s_data->a2;
	ityp * xfun = s_data->a3;
	
    dim_typ i, j;
    ityp *lfun;
    ityp *llfun;
    ityp t;

    lfun = ( ityp * ) malloc ( nfun * sizeof ( ityp ) );
    /*
    Handle special case.
    */
    if ( n == 1 )
    {
        for ( j = 0; j < nfun; ++j )
            lfun[j] = 1.00;
        return lfun;
    }

    llfun = lagrange_value ( n, x, nfun, xfun );

    for ( j = 0; j < nfun; ++j )
    {
        t = 0.00;
        for ( i = 0; i < n; ++i )
            t += fabs ( llfun[i+j*n] );
        lfun[j] = t;
    }

    free ( llfun );
    return lfun;
}

#endif
