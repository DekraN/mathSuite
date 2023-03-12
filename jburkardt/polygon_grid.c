#ifndef __DISABLEDEEP_POLYGONGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _polygon_grid_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLYGON_GRID_POINTS computes points on a polygonal grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of subintervals.
    Input, int NV, the number of vertices in the polygon.
    Input, double V[2*NV], the coordinates of the vertices.
    Input, int NG, the number of grid points.
    Output, double POLYGON_GRID_POINTS[2*NG], the coordinates of the
    grid points.
*/
{
	const _3dtpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ nv = s_data->a1;
	const register dim_typ ng = s_data->a2;
	ityp * v = s_data->a3;
	
    dim_typ i, j, k;
    int l;
    int lp1;
    int p;
    ityp vc[2];
    ityp *xg;

    xg = ( ityp * ) malloc ( (ng<<1) * sizeof ( ityp ) );
    p = 0;
    /*
    Determine the centroid.
    */
    vc[0] = vc[1] = 0.00;
    for ( j = 0; j < nv; ++j )
    {
        vc[0] += v[0+(j<<1)];
        vc[1] += v[1+(j<<1)];
    }
    vc[0] /= ( ityp ) ( nv );
    vc[1] /= ( ityp ) ( nv );
    /*
    The centroid is the first point.
    */
    xg[0+(p<<1)] = vc[0];
    xg[1+(p<<1)] = vc[1];
    ++ p;
    /*
    Consider each triangle formed by two consecutive vertices and the centroid,
    but skip the first line of points.
    */
    for ( l = 0; l < nv; ++l )
    {
        lp1 = ( ( l + 1 ) % nv );
        for ( i = 1; i <= n; ++i )
            for ( j = 0; j <= n - i; ++j)
            {
                k = n - i - j;
                xg[0+(p<<1)] = ( ( ityp ) ( i ) * v[0+(l<<1)]
                + ( ityp ) ( j ) * v[0+lp1*2]
                + ( ityp ) ( k ) * vc[0] )
                / ( ityp ) ( n );
                xg[1+(p<<1)] = ( ( ityp ) ( i ) * v[1+(l<<1)]   + ( ityp ) ( j ) * v[1+(lp1<<1)] + ( ityp ) ( k ) * vc[1] )  / ( ityp ) ( n );
                ++ p;
            }
    }

    return xg;
}

#endif
