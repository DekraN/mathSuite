#ifndef __DISABLEDEEP_NACA

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _naca4_cambered ( void * data)
/******************************************************************************/
/*
  Purpose:
    NACA4_CAMBERED: (xu,yu), (xl,yl) for a NACA cambered 4-digit airfoil.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2014
  Author:
    John Burkardt
  Reference:
    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
    "The characteristics of 78 related airfoil sections from tests in
    the variable-density wind tunnel",
    NACA Report 460, 1933.
  Parameters:
    Input, double M, the maximum camber.
    0.0 < M.
    Input, double P, the location of maximum camber.
    0.0 < P < 1.0
    Input, double T, the maximum relative thickness.
    0.0 < T <= 1.0
    Input, double C, the chord length.
    0.0 < C.
    Input, int N, the number of sample points.
    Input, double XC[N], points along the chord length.
    0.0 <= XC(*) <= C.
    Output, double XU[N], YU[N], XL[N], YL[N], for each value of
    XC, measured along the camber line, the corresponding values (XU,YU)
    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
*/
{
	const _4itdt5pit * const s_data = data;
	ityp m = s_data->a0;
	ityp p = s_data->a1;
	ityp t = s_data->a2;
	ityp c = s_data->a3;
	const register dim_typ n = s_data->a4;
	ityp * xc = s_data->a5;
	ityp * xu = s_data->a6;
	ityp * yu = s_data->a7;
	ityp * xl = s_data->a8;
	ityp * yl = s_data->a9;
	
    ityp divisor;
    ityp dycdx;
    ityp theta;
    ityp yc;
    ityp yt;

    for (dim_typ i = 0; i < n; ++i )
    {
        if ( 0.0 <= xc[i] / c && xc[i] / c <= p )
            divisor = p * p;
        else if ( p <= xc[i] / c && xc[i] / c <= 1.00 )
            divisor = pow ( 1.00 - p, 2 );
        else
            divisor = 1.00;

        dycdx = 2.00 * m * ( p - xc[i] / c ) / divisor;
        theta = atan ( dycdx );
        yt = 5.00 * t * c * ( 0.2969 * sqrt ( xc[i] / c ) + (((( - 0.1015 ) * ( xc[i] / c ) + 0.2843 ) * ( xc[i] / c ) - 0.3516 ) * ( xc[i] / c ) - 0.1260 ) * ( xc[i] / c ) );

        if ( 0.00 <= xc[i] / c && xc[i] / c <= p )
            yc = m * xc[i] * ( 2.00 * p - xc[i] / c ) / p / p;
        else if ( p <= xc[i] / c && xc[i] / c <= 1.00 )
            yc = m * ( xc[i] - c ) * ( 2.00 * p - xc[i] / c - 1.00 )/ ( 1.00 - p ) / ( 1.00 - p );
        else
            yc = 0.00;

        xu[i] = xc[i] - yt * sin ( theta );
        yu[i] = yc + yt * cos ( theta );
        xl[i] = xc[i] + yt * sin ( theta );
        yl[i] = yc - yt * cos ( theta );
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _naca4_symmetric ( void * data)
/******************************************************************************/
/*
  Purpose:
    NACA4_SYMMETRIC evaluates y(x) for a NACA symmetric 4-digit airfoil.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 May 2014
  Author:
    John Burkardt
  Reference:
    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
    "The characteristics of 78 related airfoil sections from tests in
    the variable-density wind tunnel",
    NACA Report 460, 1933.
  Parameters:
    Input, double T, the maximum relative thickness.
    Input, double C, the chord length.
    Input, int N, the number of sample points.
    Input, double X[N], points along the chord length.
    0.0 <= X(*) <= C.
    Output, double NACA4_SYMMETRIC[N], for each value of X, the corresponding
    value of Y so that (X,Y) is on the upper wing surface, and (X,-Y) is on the
    lower wing surface.
*/
{
	const dt2itpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp t = s_data->a1;
	const register ityp c = s_data->a2;
	ityp * x = s_data->a3;
	
    ityp *y = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i)
        y[i] = 5.0 * t * c * (0.2969 * sqrt ( x[i] / c )+ ((((- 0.1015 ) * ( x[i] / c )+ 0.2843 ) * ( x[i] / c )- 0.3516 ) * ( x[i] / c )- 0.1260 ) * ( x[i] / c ) );

    return y;
}

#endif
