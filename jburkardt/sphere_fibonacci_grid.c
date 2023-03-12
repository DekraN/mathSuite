#ifndef __DISABLEDEEP_SPHEREFIBONACCIGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_fibonacci_grid_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_FIBONACCI_GRID_POINTS computes sphere points on a Fibonacci spiral.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2015
  Author:
    John Burkardt
  Reference:
    Richard Swinbank, James Purser,
    Fibonacci grids: A novel approach to global modelling,
    Quarterly Journal of the Royal Meteorological Society,
    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
  Parameters:
    Input, int NG, the number of points.
    Output, double SPHERE_FIBONACCI_GRID_POINTS[3*NG], the Fibonacci
  spiral points.
*/
{
	const register dim_typ ng = *(dim_typ *) data;
	
    ityp cphi;
    dim_typ i, j;
    ityp i_r8;
    ityp ng_r8;
    ityp r8_phi;
    ityp sphi;
    ityp theta;
    ityp *xyz;

    xyz = ( ityp * ) malloc ( 3 * ng * sizeof ( ityp ) );

    r8_phi = ( 1.00 + sqrt ( 5.00 ) ) / 2.00;
    ng_r8 = ( ityp ) ( ng );

    for ( j = 0; j < ng; j++ )
    {
        i_r8 = ( ityp ) ( - ng + 1 + (j<<1) );
        theta = M_2TPI * i_r8 / r8_phi;
        sphi = i_r8 / ng_r8;
        cphi = sqrt ( ( ng_r8 + i_r8 ) * ( ng_r8 - i_r8 ) ) / ng_r8;
        xyz[0+j*3] = cphi * sin ( theta );
        xyz[1+j*3] = cphi * cos ( theta );
        xyz[2+j*3] = sphi;
    }

    return xyz;
}

#endif

