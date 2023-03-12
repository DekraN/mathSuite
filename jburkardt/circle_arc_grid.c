#ifndef __DISABLEDEEP_CIRCLEARCGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _circle_arc_grid ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_ARC_GRID computes grid points along a circular arc.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    Input, double C[2], the coordinates of the center.
    Input, double A[2], the angle of the first and last
    points, in DEGREES.
    Input, int N, the number of points.
    Output, double CIRCLE_ARC_GRID[2*N], the grid points.
*/
{
	const dtit3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp r = s_data->a1;
	ityp * c = s_data->a2;
	ityp * a = s_data->a3;
	ityp * xy = s_data->a4;
	
    ityp aj;
    dim_typ j;

    for ( j = 0; j < n; ++j )
    {
        aj = ( ( ityp ) ( n - j - 1 ) * a[0] + ( ityp ) (     j     ) * a[1] ) / ( ityp ) ( n     - 1 );

        xy[0+(j<<1)] = c[0] + r * cos ( aj * M_PI / 180.00 );
        xy[1+(j<<1)] = c[1] + r * sin ( aj * M_PI / 180.00 );
    }

    return xy;
}

#endif
