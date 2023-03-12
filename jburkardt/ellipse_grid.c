#ifndef __DISABLEDEEP_ELLIPSEGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ellipse_grid ( void * data)
/******************************************************************************/
/*
  Purpose:
    ELLIPSE_GRID generates grid points inside an ellipse.
  Discussion:
    The ellipse is specified as
  ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
    The user supplies a number N.  There will be N+1 grid points along
    the shorter axis.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of subintervals.
    Input, double R[2], the half axis lengths.
    Input, double C[2], the center of the ellipse.
    Input, int NG, the number of grid points inside the ellipse.
    Output, double ELLIPSE_GRID[2*NG], the grid points.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ ng = s_data->a1;
	ityp * r = s_data->a2;
	ityp * c = s_data->a3;
	
	
    ityp h;
    dim_typ i;
    dim_typ j;
    dim_typ ni;
    dim_typ nj;
    dim_typ p;
    ityp x;
    ityp *xy = ( ityp * ) malloc ( (ng<<1) * sizeof ( ityp ) );
    ityp y;

    if ( r[0] < r[1] )
    {
        h = 2.00 * r[0] / ( ityp ) ( (n<<1) + 1 );
        ni = n;
        nj = i4_ceiling ( r[1] / r[0] ) * ( ityp ) ( n );
    }
    else
    {
        h = 2.00 * r[1] / ( ityp ) ( (n<<1) + 1 );
        nj = n;
        ni = i4_ceiling ( r[0] / r[1] ) * ( ityp ) ( n );
    }

    p = 0;

    for ( j = 0; j <= nj; ++j )
    {
        i = 0;
        x = c[0];
        y = c[1] + ( ityp ) ( j ) * h;

        xy[0+p*2] = x;
        xy[1+p*2] = y;
        p = p + 1;

        if ( 0 < j )
        {
            xy[0+p*2] = x;
            xy[1+p*2] = 2.0 * c[1] - y;
            ++ p;
        }

        for ( ; ; )
        {
            ++ i;
            x = c[0] + ( ityp ) ( i ) * h;

            if ( 1.00 < pow ( ( x - c[0] ) / r[0], 2 ) + pow ( ( y - c[1] ) / r[1], 2 ) )
                break;

            xy[0+p*2] = x;
            xy[1+p*2] = y;
            p = p + 1;
            xy[0+p*2] = 2.00 * c[0] - x;
            xy[1+p*2] = y;
            ++ p;

            if ( 0 < j )
            {
                xy[0+p*2] = x;
                xy[1+p*2] = 2.00 * c[1] - y;
                ++ p;
                xy[0+p*2] = 2.00 * c[0] - x;
                xy[1+p*2] = 2.00 * c[1] - y;
                ++ p;
            }
        }
    }
    return xy;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ellipse_grid_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    ELLIPSE_GRID_COUNT counts the grid points inside an ellipse.
  Discussion:
    The ellipse is specified as
  ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
    The user supplies a number N.  There will be N+1 grid points along
    the shorter axis.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of subintervals.
    Input, double R[2], the half axis lengths.
    Input, double C[2], the center of the ellipse.
    Output, int ELLIPSE_GRID)_COUNT, the number of grid points inside
    the ellipse.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * r = s_data->a1;
	ityp * c = s_data->a2;
	
    ityp h;
    dim_typ i;
    dim_typ j;
    dim_typ ni;
    dim_typ nj;
    dim_typ p;
    ityp x;
    ityp y;

    if ( r[0] < r[1] )
    {
        h = 2.00 * r[0] / ( ityp ) ( (n<<1) + 1 );
        ni = n;
        nj = i4_ceiling ( r[1] / r[0] ) * ( ityp ) ( n );
    }
    else
    {
        h = 2.00 * r[1] / ( ityp ) ( (n<<1) + 1 );
        nj = n;
        ni = i4_ceiling ( r[0] / r[1] ) * ( ityp ) ( n );
    }

    p = 0;

    for ( j = 0; j <= nj; ++j )
    {
        i = 0;
        x = c[0];
        y = c[1] + ( ityp ) ( j ) * h;

        ++ p;

        if ( 0 < j )
            ++ p;

        for ( ; ; )
        {
            ++ i;
            x = c[0] + ( ityp ) ( i ) * h;

            if ( 1.00 < pow ( ( x - c[0] ) / r[0], 2 ) + pow ( ( y - c[1] ) / r[1], 2 ) )
                break;
            ++ p;
            ++ p;

            if ( 0 < j )
            {
                ++ p;
                ++ p;
            }
        }
    }
    
    result = p;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_ceiling ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_CEILING rounds an R8 up to the nearest I4.
  Discussion:
    The "ceiling" of X is the value of X rounded towards plus infinity.
  Example:
    X        I4_CEILING(X)
   -1.1      -1
   -1.0      -1
   -0.9       0
   -0.1       0
    0.0       0
    0.1       1
    0.9       1
    1.0       1
    1.1       2
    2.9       3
    3.0       3
    3.14159   4
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, double X, the number whose ceiling is desired.
    Output, int I4_CEILING, the ceiling of X.
*/
{
	static int result = INT_MAX;
	
	const register ityp x = *(ityp *) data;
	
    int value = ( int ) x;
    if ( value < x )
        ++ value;
        
    result = value;
    return &result;
}

#endif
