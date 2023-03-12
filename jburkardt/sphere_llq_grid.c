#ifndef __DISABLEDEEP_SPHERELLQGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_llq_grid_point_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LLQ_GRID_POINT_COUNT counts points for a latitude/longitude grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int LAT_NUM, LONG_NUM, the number of latitude
    and longitude lines to draw.  The latitudes do not include the North and
    South poles, which will be included automatically, so LAT_NUM = 5, for
    instance, will result in points along 7 lines of latitude.
    Output, int SPHERE_LLQ_GRID_POINT_COUNT, the number of grid points.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ lat_num = a_data[0];
	const register dim_typ long_num = a_data[1];
	
	result = 2 + lat_num * long_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_llq_grid_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LLQ_GRID_POINTS produces points on a latitude/longitude grid.
  Discussion:
    A SPHERE LLQ grid imposes a grid of quadrilaterals on a sphere,
    using latitude and longitude lines.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the sphere.
    Input, double PC[3], the coordinates of the center of the sphere.
    Input, int LAT_NUM, LON_NUM, the number of latitude and longitude
    lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so LAT_NUM = 5, for instance,
    will result in points along 7 lines of latitude.
    Input, int POINT_NUM, the number of points.
    Output, double SPHERE_LLQ_GRID_POINTS[3*POINT_NUM], the coordinates
    of the grid points.
*/
{
	const itpit3dt * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * pc = s_data->a1;
	const register dim_typ lat_num = s_data->a2;
	const register dim_typ lon_num = s_data->a3;
	const register dim_typ point_num = s_data->a4;
	
    dim_typ lat;
    dim_typ lon;
    dim_typ n;
    ityp *p;
    ityp phi;
    ityp theta;

    p = ( ityp * ) malloc ( 3 * point_num * sizeof ( ityp ) );
    n = 0;
    /*
    The north pole.
    */
    theta = phi = 0.00;

    p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
    p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
    p[2+n*3] = pc[2] + r * cos ( phi );
    n = n + 1;
    /*
    Do each intermediate ring of latitude.
    */
    for ( lat = 1; lat <= lat_num; ++lat )
    {
        phi = ( ityp ) ( lat ) * M_PI / ( ityp ) ( lat_num + 1 );
        /*
        Along that ring of latitude, compute points at various longitudes.
        */
        for ( lon = 0; lon < lon_num; ++lon )
        {
            theta = ( ityp ) ( lon ) * M_2TPI / ( ityp ) ( lon_num );

            p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
            p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
            p[2+n*3] = pc[2] + r * cos ( phi );
            ++ n;
        }
    }
    /*
    The south pole.
    */
    theta = 0.00;
    phi = M_PI;
    p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
    p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
    p[2+n*3] = pc[2] + r * cos ( phi );
    ++ n;

    return p;
}

#endif
