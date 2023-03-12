#ifndef __DISABLEDEEP_INTERP

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cc_abscissas ( void * data)
/******************************************************************************/
/*
  Purpose:
    CC_ABSCISSAS computes the Clenshaw Curtis abscissas.
  Discussion:
    The interval is [ -1, 1 ].
    The abscissas are the cosines of equally spaced angles between
    180 and 0 degrees, including the endpoints.
      X(I) = cos ( ( ORDER - I ) * M_PI / ( ORDER - 1 ) )
    except for the basic case ORDER = 1, when
      X(1) = 0.
    If the value of ORDER is increased in a sensible way, then
    the new set of abscissas will include the old ones.  One such
    sequence would be ORDER(K) = 2*K+1 for K = 0, 1, 2, ...
    When doing interpolation with Lagrange polynomials, the Clenshaw Curtis
    abscissas can be a better choice than regularly spaced abscissas.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Reference:
    Charles Clenshaw, Alan Curtis,
    A Method for Numerical Integration on an Automatic Computer,
    Numerische Mathematik,
    Volume 2, Number 1, December 1960, pages 197-205.
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  Parameters:
    Input, int N, the order of the rule.
    Output, double CC_ABSCISSAS[N], the abscissas.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp theta;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
    {
        x[0] = 0.00;
        return x;
    }

    for ( i = 0; i < n; ++i )
    {
        theta = M_PI * ( ityp ) ( n - 1 - i )/ ( ityp ) ( n - 1 );
        x[i] = cos ( theta );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cc_abscissas_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    CC_ABSCISSAS_AB computes Clenshaw Curtis abscissas for the interval [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, the endpoints of the interval.
    Input, int N, the order of the rule.
    Output, double CC_ABSCISSAS_AB[N], the abscissas.
*/
{
	const _2itdt * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    dim_typ i;
    ityp theta;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
    {
        x[0] = 0.50 * ( b + a );
        return x;
    }

    for ( i = 0; i < n; ++i )
    {
        theta = M_PI * ( ityp ) ( n - 1 - i )/ ( ityp ) ( n - 1 );
        x[i] = 0.50 * ( ( b + a ) + ( b - a ) * cos ( theta ) );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _f1_abscissas ( void * data)
/******************************************************************************/
/*
  Purpose:
    F1_ABSCISSAS computes Fejer type 1 abscissas.
  Discussion:
    The interval is [ -1, +1 ].
    The abscissas are the cosines of equally spaced angles, which
    are the midpoints of N equal intervals between 0 and M_PI.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  Parameters:
    Input, int N, the order of the rule.
    Output, double F1_ABSCISSAS[N], the abscissas.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp theta;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
    {
        x[0] = 0.00;
        return x;
    }

    for ( i = 0; i < n; ++i )
    {
        theta = ( ityp ) ( (n<<1) - (i<<1) - 1 ) * M_PI / ( ityp ) ( n<<1             );
        x[i] = cos ( theta );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _f1_abscissas_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    F1_ABSCISSAS_AB computes Fejer type 1 abscissas for the interval [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  Parameters:
    Input, double A, B, the endpoints of the interval.
    Input, int N, the order of the rule.
    Output, double F1_ABSCISSAS_AB[N], the abscissas.
*/
{
	const _2itdt * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    dim_typ i;
    ityp theta;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
    {
        x[0] = 0.50 * ( b + a );
        return x;
    }

    for ( i = 0; i < n; ++i )
    {
        theta = ( ityp ) ((n<<1) - (i<<1) - 1 ) * M_PI / ( ityp ) ( n<<1             );
        x[i] = 0.50 * ( ( b + a ) + ( b - a ) * cos ( theta ) );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _f2_abscissas ( void * data)
/******************************************************************************/
/*
  Purpose:
    F2_ABSCISSAS computes Fejer Type 2 abscissas.
  Discussion:
    The interval is [-1,+1].
    The abscissas are the cosines of equally spaced angles.
    The angles are computed as N+2 equally spaced values between 0 and M_PI,
    but with the first and last angle omitted.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  Parameters:
    Input, int N, the order of the rule.
    Output, double F2_ABSCISSAS[N], the abscissas.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp theta;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
    {
        x[0] = 0.00;
        return x;
    }
    else if ( n == 2 )
    {
        x[0] = -0.50;
        x[1] =  0.50;
        return x;
    }

    for ( i = 0; i < n; ++i )
    {
        theta = ( ityp ) ( n - i ) *M_PI / ( ityp ) ( n + 1 );
        x[i] = cos ( theta );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _f2_abscissas_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    F2_ABSCISSAS_AB computes Fejer Type 2 abscissas for the interval [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, the endpoints of the interval.
    Input, int N, the order of the rule.
    Output, double F2_ABSCISSAS_AB[N], the abscissas.
*/
{
	const _2itdt * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    dim_typ i;
    ityp theta;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        theta = ( ityp ) ( n - i ) * M_PI / ( ityp ) ( n + 1 );
        x[i] = 0.50 * ( ( b + a ) + ( b - a ) * cos ( theta ) );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _interp_lagrange ( void * data)
/******************************************************************************/
/*
  Purpose:
    INTERP_LAGRANGE: Lagrange polynomial interpolant to a curve in M dimensions.
  Discussion:
    From a space of M dimensions, we are given a sequence of
    DATA_NUM points, which are presumed to be successive samples
    from a curve of points P.
    We are also given a parameterization of this data, that is,
    an associated sequence of DATA_NUM values of a variable T.
    Thus, we have a sequence of values P(T), where T is a scalar,
    and each value of P is of dimension M.
    We are then given INTERP_NUM values of T, for which values P
    are to be produced, by linear interpolation of the data we are given.
    The user may request extrapolation.  This occurs whenever
    a T_INTERP value is less than the minimum T_DATA or greater than the
    maximum T_DATA.  In that case, extrapolation is used.
    For each spatial component, a polynomial of degree
 ( DATA_NUM - 1 ) is generated for the interpolation.  In most cases,
    such a polynomial interpolant begins to oscillate as DATA_NUM
    increases, even if the original data seems well behaved.  Typically,
    values of DATA_NUM should be no greater than 10
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int DATA_NUM, the number of data points.
    Input, double T_DATA[DATA_NUM], the value of the
    independent variable at the sample points.
    Input, double P_DATA[M*DATA_NUM], the value of the
    dependent variables at the sample points.
    Input, int INTERP_NUM, the number of points
    at which interpolation is to be done.
    Input, double T_INTERP[INTERP_NUM], the value of the
    independent variable at the interpolation points.
    Output, double INTERP_LAGRANGE[M*DATA_NUM], the interpolated
    values of the dependent variables at the interpolation points.
*/
{
	const _3dt3pit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ data_num = s_data->a1;
	const register dim_typ interp_num = s_data->a2;
	ityp * t_data = s_data->a3;
	ityp * p_data = s_data->a4;
	ityp * t_interp = s_data->a5;	

    ityp *l_interp;
    ityp *p_interp;
    /*
    Evaluate the DATA_NUM Lagrange polynomials associated with T_DATA(1:DATA_NUM)
    for the interpolation points T_INTERP(1:INTERP_NUM).
    */
    l_interp = lagrange_value ( data_num, t_data, interp_num, t_interp );
    /*
    Multiply P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
    to get P_INTERP(1:M,1:INTERP_NUM).
    */
    p_interp = r8mat_mm_new ( m, data_num, interp_num, p_data, l_interp );
    free ( l_interp );
    return p_interp;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _interp_linear ( void * data)
/******************************************************************************/
/*
  Purpose:
    INTERP_LINEAR: piecewise linear interpolation to a curve in M dimensions.
  Discussion:
    From a space of M dimensions, we are given a sequence of
    DATA_NUM points, which are presumed to be successive samples
    from a curve of points P.
    We are also given a parameterization of this data, that is,
    an associated sequence of DATA_NUM values of a variable T.
    The values of T are assumed to be strictly increasing.
    Thus, we have a sequence of values P(T), where T is a scalar,
    and each value of P is of dimension M.
    We are then given INTERP_NUM values of T, for which values P
    are to be produced, by linear interpolation of the data we are given.
    Note that the user may request extrapolation.  This occurs whenever
    a T_INTERP value is less than the minimum T_DATA or greater than the
    maximum T_DATA.  In that case, linear extrapolation is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int DATA_NUM, the number of data points.
    Input, double T_DATA[DATA_NUM], the value of the
    independent variable at the sample points.  The values of T_DATA
    must be strictly increasing.
    Input, double P_DATA[M*DATA_NUM], the value of the
    dependent variables at the sample points.
    Input, int INTERP_NUM, the number of points
    at which interpolation is to be done.
    Input, double T_INTERP[INTERP_NUM], the value of the
    independent variable at the interpolation points.
    Output, double INTERP_LINEAR[M*DATA_NUM], the interpolated
    values of the dependent variables at the interpolation points.
*/
{
	const _3dt3pit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ data_num = s_data->a1;
	const register dim_typ interp_num = s_data->a2;
	ityp * t_data = s_data->a3;
	ityp * p_data = s_data->a4;
	ityp * t_interp = s_data->a5;	
	
    dim_typ i;
    int interp;
    int left;
    ityp *p_interp;
    int right;
    ityp t;

    if ( ! r8vec_ascends_strictly ( data_num, t_data ) )
        return NULL;

    p_interp = ( ityp * ) malloc ( m * interp_num * sizeof ( ityp ) );

    for ( interp = 0; interp < interp_num; ++interp )
    {
        t = t_interp[interp];
        /*
        Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
        nearest to, TVAL.
        */
        r8vec_bracket0 ( data_num, t_data, t, &left, &right );

        for ( i = 0; i < m; ++i )
            p_interp[i+interp*m] = ( ( t_data[right] - t                ) * p_data[i+left*m]   + (                 t - t_data[left] ) * p_data[i+right*m] ) / ( t_data[right]     - t_data[left] );
    }

    return p_interp;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _interp_nearest ( void * data)
/******************************************************************************/
/*
  Purpose:
    INTERP_NEAREST: Nearest neighbor interpolation to a curve in M dimensions.
  Discussion:
    From a space of M dimensions, we are given a sequence of
    DATA_NUM points, which are presumed to be successive samples
    from a curve of points P.
    We are also given a parameterization of this data, that is,
    an associated sequence of DATA_NUM values of a variable T.
    Thus, we have a sequence of values P(T), where T is a scalar,
    and each value of P is of dimension M.
    We are then given INTERP_NUM values of T, for which values P
    are to be produced, by nearest neighbor interpolation.
    The user may request extrapolation.  This occurs whenever
    a T_INTERP value is less than the minimum T_DATA or greater than the
    maximum T_DATA.  In that case, extrapolation is used.
    The resulting interpolant is piecewise constant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int DATA_NUM, the number of data points.
    Input, double T_DATA[DATA_NUM], the value of the
    independent variable at the sample points.
    Input, double P_DATA[M*DATA_NUM], the value of the
    dependent variables at the sample points.
    Input, int INTERP_NUM, the number of points
    at which interpolation is to be done.
    Input, double T_INTERP[INTERP_NUM], the value of the
    independent variable at the interpolation points.
    Output, double INTERP_NEAREST[M*DATA_NUM], the interpolated
    values of the dependent variables at the interpolation points.
*/
{
	const _3dt3pit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ data_num = s_data->a1;
	const register dim_typ interp_num = s_data->a2;
	ityp * t_data = s_data->a3;
	ityp * p_data = s_data->a4;
	ityp * t_interp = s_data->a5;	
	
    dim_typ i;
    dim_typ jd;
    dim_typ ji;
    ityp *p_interp;
    /*
    For each interpolation point, find the index of the nearest data point.
    */
    p_interp = ( double * ) malloc ( m * interp_num * sizeof ( double ) );

    for ( ji = 0; ji < interp_num; ++ji )
    {
        jd = r8vec_sorted_nearest0 ( data_num, t_data, t_interp[ji] );
        for ( i = 0; i < m; ++i )
            p_interp[i+ji*m] = p_data[i+jd*m];
    }

    return p_interp;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_VALUE evaluates the Lagrange polynomials.
  Discussion:
    Given DATA_NUM distinct abscissas, T_DATA(1:DATA_NUM),
    the I-th Lagrange polynomial L(I)(T) is defined as the polynomial of
    degree DATA_NUM - 1 which is 1 at T_DATA(I) and 0 at the DATA_NUM - 1
    other abscissas.
    A formal representation is:
      L(I)(T) = Product ( 1 <= J <= DATA_NUM, I /= J )
 ( T - T(J) ) / ( T(I) - T(J) )
    This routine accepts a set of INTERP_NUM values at which all the Lagrange
    polynomials should be evaluated.
    Given data values P_DATA at each of the abscissas, the value of the
    Lagrange interpolating polynomial at each of the interpolation points
    is then simple to compute by matrix multiplication:
      P_INTERP(1:INTERP_NUM) =
        P_DATA(1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
    or, in the case where P is multidimensional:
      P_INTERP(1:M,1:INTERP_NUM) =
        P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int DATA_NUM, the number of data points.
    DATA_NUM must be at least 1.
    Input, double T_DATA[DATA_NUM], the data points.
    Input, int INTERP_NUM, the number of
    interpolation points.
    Input, double T_INTERP[INTERP_NUM], the
    interpolation points.
    Output, double LAGRANGE_VALUE[DATA_NUM*INTERP_NUM], the values
    of the Lagrange polynomials at the interpolation points.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ data_num = s_data->a0;
	const register dim_typ interp_num = s_data->a1;
	ityp * t_data = s_data->a2;
	ityp * t_interp = s_data->a3;	
	
    dim_typ i;
    dim_typ i1;
    dim_typ i2;
    dim_typ j;
    ityp *l_interp;

    l_interp = ( ityp * ) malloc ( data_num * interp_num * sizeof ( ityp ) );
    /*
    Evaluate the polynomial.
    */
    for ( j = 0; j < interp_num; ++j )
        for ( i = 0; i < data_num; ++i )
            l_interp[i+j*data_num] = 1.00;

    for ( i1 = 0; i1 < data_num; ++i1 )
        for ( i2 = 0; i2 < data_num; ++i2  )
            if ( i1 != i2 )
                for ( j = 0; j < interp_num; ++j )
                    l_interp[i1+j*data_num] *= ( t_interp[j] - t_data[i2] ) / ( t_data[i1] - t_data[i2] );

    return l_interp;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ncc_abscissas ( void * data)
/******************************************************************************/
/*
  Purpose:
    NCC_ABSCISSAS computes the Newton Cotes Closed abscissas.
  Discussion:
    The interval is [ -1, 1 ].
    The abscissas are the equally spaced points between -1 and 1,
    including the endpoints.
    If N is 1, however, the single abscissas is X = 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the rule.
    Output, double NCC_ABSCISSAS[N], the abscissas.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
    {
        x[0] = 0.00;
        return x;
    }

    for ( i = 0; i < n; ++i )
        x[i] = ( ( ityp ) ( n - i - 1 ) * ( -1.00 )+ ( ityp ) (     i     ) * ( +1.00 ) )/ ( ityp ) ( n     - 1 );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ncc_abscissas_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    NCC_ABSCISSAS_AB computes the Newton Cotes Closed abscissas for [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, the endpoints of the interval.
    Input, int N, the order of the rule.
    Output, double NCC_ABSCISSAS_AB[N], the abscissas.
*/
{
	const _2itdt * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    dim_typ i;
    ityp *x;

    if ( n < 1 )
        return NULL;
    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
    {
        x[0] = 0.50 * ( b + a );
        return x;
    }

    for ( i = 0; i < n; ++i )
        x[i] = ( ( ityp ) ( n - i - 1 ) * ( a )+ ( ityp ) (     i     ) * ( b ) )/ ( ityp ) ( n     - 1 );

    return x;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _nco_abscissas ( void * data)
/******************************************************************************/
/*
  Purpose:
    NCO_ABSCISSAS computes the Newton Cotes Open abscissas.
  Discussion:
    The interval is [ -1, 1 ].
    The abscissas are the equally spaced points between -1 and 1,
    not including the endpoints.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the rule.
    Output, double NCO_ABSCISSAS[N], the abscissas.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        x[i] = ( ( ityp ) ( n - i     ) * ( -1.00 )+ ( ityp ) (     i + 1 ) * ( +1.00 ) )/ ( ityp ) ( n     + 1 );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _nco_abscissas_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    NCO_ABSCISSAS_AB computes the Newton Cotes Open abscissas for [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, the endpoints of the interval.
    Input, int N, the order of the rule.
    Output, double NCO_ABSCISSAS_AB[N], the abscissas.
*/
{
	const _2itdt * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    dim_typ i;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        x[i] = ( ( ityp ) ( n - i     ) * ( a )+ ( ityp ) (     i + 1 ) * ( b ) )/ ( ityp ) ( n     + 1 );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _parameterize_arc_length ( void * data)
/******************************************************************************/
/*
  Purpose:
    PARAMETERIZE_ARC_LENGTH parameterizes data by pseudo-arclength.
  Discussion:
    A parameterization is required for the interpolation.
    This routine provides a parameterization by computing the
    pseudo-arclength of the data, that is, the Euclidean distance
    between successive points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int DATA_NUM, the number of data points.
    Input, double P_DATA[M*DATA_NUM], the data values.
    Output, double PARAMETERIZE_ARC_LENGTH[DATA_NUM], parameter values
    assigned to the data.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ data_num = s_data->a1;
	ityp * p_data = s_data->a2;
	
    dim_typ i, j;
    ityp t;
    ityp *t_data;

    t_data = ( ityp * ) malloc ( data_num * sizeof ( ityp ) );
    t_data[0] = 0.00;

    for ( j = 1; j < data_num; ++j )
    {
        t = 0.00;
        for ( i = 0; i < m; ++i )
            t += pow ( p_data[i+j*m] - p_data[i+(j-1)*m], 2 );
        t_data[j] = t_data[j-1] + sqrt ( t );
    }

    return t_data;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _parameterize_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    PARAMETERIZE_INDEX parameterizes data by its index.
  Discussion:
    A parameterization is required for the interpolation.
    This routine provides a naive parameterization by vector index.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int DATA_NUM, the number of data points.
    Input, double P_DATA[M*DATA_NUM], the data values.
    Output, double PARAMETERIZE_INDEX[DATA_NUM], parameter values
    assigned to the data.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ data_num = s_data->a1;
	ityp * p_data = s_data->a2;
	
    ityp * t_data = ( ityp * ) malloc ( data_num * sizeof ( ityp ) );
    for (dim_typ j = 0; j < data_num; ++j)
        t_data[j] = ( ityp ) j;
    return t_data;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_expand_linear2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
    In this version of the routine, the expansion is indicated
    by specifying the dimensions of the expanded array.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 August 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in A.
    Input, double A(M,N), a "small" M by N array.
    Input, int M2, N2, the number of rows and columns in A2.
    Output, double R8MAT_EXPAND_LINEAR2[M2*N2], the expanded array,
    which contains an interpolated version of the data in A.
*/
{
	const _2dtpit2dt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	const register dim_typ m2 = s_data->a3;
	const register dim_typ n2 = s_data->a4;
	
    ityp *a2;
    dim_typ i;
    dim_typ i1;
    dim_typ i2;
    dim_typ j;
    dim_typ j1;
    dim_typ j2;
    ityp r;
    ityp r1;
    ityp r2;
    ityp s;
    ityp s1;
    ityp s2;

    a2 = ( ityp * ) malloc ( m2 * n2 * sizeof ( ityp ) );

    for ( i = 1; i <= m2; ++i )
    {
        r = m2 == 1 ? 0.50 : ( ityp ) ( i - 1 ) / ( ityp ) ( m2 - 1 );
        i1 = 1 + ( int ) ( r * ( ityp ) ( m - 1 ) );
        i2 = i1 + 1;

        if ( m < i2 )
        {
            i1 = m - 1;
            i2 = m;
        }

        r1 = ( ityp ) ( i1 - 1 ) / ( ityp ) ( m - 1 );
        r2 = ( ityp ) ( i2 - 1 ) / ( ityp ) ( m - 1 );

        for ( j = 1; j <= n2; ++j )
        {
            s = n2 == 1 ? 0.50 : ( ityp ) ( j - 1 ) / ( ityp ) ( n2 - 1 );
            j1 = 1 + ( int ) ( s * ( ityp ) ( n - 1 ) );
            j2 = j1 + 1;

            if ( n < j2 )
            {
                j1 = n - 1;
                j2 = n;
            }

            s1 = ( ityp ) ( j1 - 1 ) / ( ityp ) ( n - 1 );
            s2 = ( ityp ) ( j2 - 1 ) / ( ityp ) ( n - 1 );

            a2[i-1+(j-1)*m2] =( ( r2 - r ) * ( s2 - s ) * a[i1-1+(j1-1)*m]+ ( r - r1 ) * ( s2 - s ) * a[i2-1+(j1-1)*m]+ ( r2 - r ) * ( s - s1 ) * a[i1-1+(j2-1)*m]+ ( r - r1 ) * ( s - s1 ) * a[i2-1+(j2-1)*m] )/ ( ( r2 - r1 ) * ( s2 - s1 ) );
        }
    }

    return a2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_ascends_strictly ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
  Discussion:
    An R8VEC is a vector of R8's.
    Notice the effect of entry number 6 in the following results:
      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the array.
    Input, double X[N], the array to be examined.
    Output, int R8VEC_ASCENDS_STRICTLY, is TRUE if the
    entries of X strictly ascend.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    for (dim_typ i = 0; i < n - 1; ++i )
        if ( x[i+1] <= x[i] )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_bracket0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_BRACKET searches a sorted array for successive brackets of a value.
  Discussion:
    An R8VEC is a vector of R8's.
    If the values in the vector are thought of as defining intervals
    on the real line, then this routine searches for the interval
    nearest to or containing the given value.
    It is always true that RIGHT = LEFT+1.
    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
      XVAL   < X[0] < X[1];
    If X(1) <= XVAL < X[N-1], 
    
      X[LEFT-1] <= XVAL < X[RIGHT-1];
    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
    This routine computes indices RIGHT and LEFT that are 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, length of input array.
    Input, double X[N], an array that has been sorted into ascending order.
    Input, double XVAL, a value to be bracketed.
    Output, int *LEFT, *RIGHT, the results of the search.
*/
{
	const dtpitit2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	const register ityp xval = s_data->a2;
	int * left = s_data->a3;
	int * right = s_data->a4;
	
    for (dim_typ i = 2; i <= n - 1; ++i )
        if ( xval < x[i-1] )
        {
            *left = i - 2;
            *right = i - 1;
            return NULL;
        }

    *left = n - 2;
    *right = n - 1;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_expand_linear2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_EXPAND_LINEAR2 linearly interpolates new data into an R8VEC.
  Discussion:
    An R8VEC is a vector of R8's.
    This routine starts with a vector of data.
    The intent is to "fatten" the data, that is, to insert more points
    between successive values of the original data.
    There will also be extra points placed BEFORE the first original
    value and AFTER that last original value.
    The "fattened" data is equally spaced between the original points.
    The BEFORE data uses the spacing of the first original interval,
    and the AFTER data uses the spacing of the last original interval.
  Example:
    N = 3
    BEFORE = 3
    FAT = 2
    AFTER = 1
    X    = (/                   0.0,           6.0,             7.0       /)
    XFAT = (/ -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0, 7.66 /)
            3 "BEFORE's"        Old  2 "FATS"  Old    2 "FATS"  Old  1 "AFTER"
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of input data values.
    N must be at least 2.
    Input, double X[N], the original data.
    Input, int BEFORE, the number of "before" values.
    Input, int FAT, the number of data values to interpolate
    between each pair of original data values.
    Input, int AFTER, the number of "after" values.
    Output, double R8VEC_EXPAND_LINEAR2[BEFORE+(N-1)*(FAT+1)+1+AFTER], the
    "fattened" data.
*/
{
	const _2dtpit2dt * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ before = s_data->a1;
	ityp * x = s_data->a2;
	const register dim_typ fat = s_data->a3;
	const register dim_typ after = s_data->a4;
	
    dim_typ i, j, k;
    ityp * xfat = ( ityp * ) malloc ( (before+(n-1)*(fat+1)+1+after) * sizeof ( ityp ) );

    k = 0;
    /*
    Points BEFORE.
    */
    for ( j = 1 - before + fat; j <= fat; ++j )
    {
        xfat[k] = ( ( ityp ) ( fat - j + 1 ) * ( x[0] - ( x[1] - x[0] ) ) + ( ityp ) (       j     ) *   x[0]          ) / ( ityp ) ( fat     + 1 );
        ++ k;
    }
    /*
    Original points and FAT points.
    */
    for ( i = 0; i < n - 1; ++i )
    {
        xfat[k] = x[0];
        ++ k;
        for ( j = 1; j <= fat; ++j )
        {
            xfat[k] = ( ( ityp ) ( fat - j + 1 ) * x[i]+ ( ityp ) (       j     ) * x[i+1] ) / ( ityp ) ( fat     + 1 );
            ++ k;
        }
    }

    xfat[k] = x[n-1];
    ++ k;
    /*
    Points AFTER.
    */
    for ( j = 1; j <= after; ++j )
    {
        xfat[k] = ( ( ityp ) ( fat - j + 1 ) * x[n-1]+ ( ityp ) (       j     ) * ( x[n-1] + ( x[n-1] - x[n-2] ) ) ) / ( ityp ) ( fat     + 1 );
        ++ k;
    }

    return xfat;
}

#endif
