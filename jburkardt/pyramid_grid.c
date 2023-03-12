#ifndef __DISABLEDEEP_PYRAMIDGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pyramid_unit_grid ( void * data)
/******************************************************************************/
/*
  Purpose
    PYRAMID_UNIT_GRID computes grid points in the unit pyramid.
  Discussion:
    The unit pyramid has base (-1,-1,0), (+1,1,0), (+1,+1,0), (-1,+1,0)
    and vertex (0,0,1).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of subintervals.
    Input, int NG, the number of nodes to generate,
    as determined by pyramid_grid_count().
    Output, double PYRAMID_UNIT_GRID[3*NG], the grid point coordinates.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ ng = a_data[1];
	
    dim_typ g;
    int hi;
    dim_typ i, j, k;
    int lo;
    ityp *pg = ( ityp * ) malloc ( 3 * ng * sizeof ( ityp ) );
    g = 0;

    for ( k = n; 0 <= k; --k )
    {
        hi = n - k;
        lo = - hi;
        for ( j = lo; j <= hi; j += 2 )
            for ( i = lo; i <= hi; i += 2 )
            {
                pg[0+g*3] = ( ityp ) ( i ) / ( ityp ) ( n );
                pg[1+g*3] = ( ityp ) ( j ) / ( ityp ) ( n );
                pg[2+g*3] = ( ityp ) ( k ) / ( ityp ) ( n );
                ++ g;
            }
    }

    return pg;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *    _pyramid_unit_vertices ( void * data)
/******************************************************************************/
/*
  Purpose:
    PYRAMID_UNIT_VERTICES returns the vertices of the unit pyramid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 August 2014
  Author:
    John Burkardt
  Parameters:
    Output, double V1[3], V2[3], V3[3], V4[3], V5[3], the vertices.
*/
{
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	ityp * v3 = a_data[2];
	ityp * v4 = a_data[3];
	ityp * v5 = a_data[4];
	
    v1[0] = v1[1] = v2[2] = v3[2] = v4[2] = v5[2] = 0.00;
    v1[2] = v3[0] = v4[0] = v4[1] = v5[1] = 1.00;
    v2[0] = v2[1] = v3[1] = v5[0] = -1.00;
    return NULL;
}
#endif
