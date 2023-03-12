#ifndef __DISABLEDEEP_SPHERELLTGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_llt_grid_line_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LLT_GRID_LINE_COUNT counts latitude/longitude triangle grid lines.
  Discussion:
    A SPHERE LLT grid imposes a grid of triangles on a sphere,
    using latitude and longitude lines.
    The number returned is the number of pairs of points to be connected.
  Licensing:
    This code is distributed under the GNU LGPL license.   Modified:
    20 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int LAT_NUM, LONG_NUM, the number of latitude and
    longitude lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so LAT_NUM = 5, for instance,
    will result in points along 7 lines of latitude.
    Output, int SPHERE_LLT_GRID_LINE_COUNT, the number of grid lines.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ lat_num = a_data[0];
	const register dim_typ long_num = a_data[1];
	
	result = long_num * ( lat_num + 1 ) + long_num *   lat_num + long_num * ( lat_num - 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_llt_grid_lines ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LLT_GRID_LINES: latitude/longitude triangle grid lines.
  Discussion:
    A SPHERE LLT grid imposes a grid of triangles on a sphere,
    using latitude and longitude lines.
    The point numbering system is the same used in SPHERE_LLT_POINTS,
    and that routine may be used to compute the coordinates of the points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int NLAT, NLONG, the number of latitude and longitude
    lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so NLAT = 5, for instance,
    will result in points along 7 lines of latitude.
    Input, int LINE_NUM, the number of grid lines.
    Output, int SPHERE_LLT_GRID_LINES[2*LINE_NUM], contains pairs of point
    indices for line segments that make up the grid.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ nlat = a_data[0];
	const register dim_typ nlong = a_data[1];
	const register dim_typ line_num = a_data[2];
	
    dim_typ i;
    dim_typ j;
    dim_typ l;
    int *line;
    dim_typ next;
    dim_typ newcol;
    dim_typ old;

    line = ( int * ) malloc ( line_num * sizeof ( int ) << 1 );

    l = 0;
    /*
    "Vertical" lines.
    */
    for ( j = 0; j <= nlong - 1; ++j )
    {
        old = 0;
        next = j + 1;
        line[0+(l<<1)] = old;
        line[1+(l<<1)] = next;
        ++ l;

        for ( i = 1; i <= nlat-1; ++i )
        {
            old = next;
            next = old + nlong;
            line[0+(l<<1)] = old;
            line[1+(l<<1)] = next;
            ++ l;
        }

        old = next;
        line[0+(l<<1)] = old;
        line[1+(l<<1)] = 1 + nlat * nlong;
        ++ l;
    }
    /*
    "Horizontal" lines.
    */
    for ( i = 1; i <= nlat; ++i )
    {
        next = ( i - 1 ) * nlong + 1;

        for ( j = 0; j <= nlong - 2; ++j )
        {
            old = next;
            next = old + 1;
            line[0+(l<<1)] = old;
            line[1+(l<<1)] = next;
            ++ l;
        }

        old = next;
        next = ( i - 1 ) * nlong + 1;
        line[0+(l<<1)] = old;
        line[1+(l<<1)] = next;
        ++ l;
    }
    /*
    "Diagonal" lines.
    */
    for ( j = 0; j < nlong; ++j )
    {
        old = 0;
        next = j + 1;
        newcol = j;

        for ( i = 1; i < nlat; ++i )
        {
            old = next;
            next = old + nlong + 1;
            ++ newcol;
            if ( nlong - 1 < newcol )
            {
                newcol = 0;
                next -= nlong;
            }

            line[0+(l<<1)] = old;
            line[1+(l<<1)] = next;
            ++ l;
        }
    }

    return line;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_llt_grid_point_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LLT_GRID_POINT_COUNT counts points for a latitude/longitude grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int LAT_NUM, LONG_NUM, the number of latitude
    and longitude lines to draw.  The latitudes do not include the North and
    South poles, which will be included automatically, so LAT_NUM = 5, for
    instance, will result in points along 7 lines of latitude.
    Output, int SPHERE_LLT_GRID_POINT_COUNT, the number of grid points.
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
__MATHSUITE __JBURKARDT  void   * _sphere_llt_grid_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LLT_GRID_POINTS produces points on a latitude/longitude grid.
  Discussion:
    A SPHERE LLT grid imposes a grid of triangles on a sphere,
    using latitude and longitude lines.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2015
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
    Output, double SPHERE_LLT_GRID_POINTS[3*POINT_NUM], the coordinates
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
