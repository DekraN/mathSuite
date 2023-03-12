#ifndef __DISABLEDEEP_SPHEREGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _atan4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ATAN4 computes the inverse tangent of the ratio Y / X.
  Discussion:
    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
    the built in functions ATAN and ATAN2 already do.
    However:
    * ATAN4 always returns a positive angle, between 0 and 2 PI,
      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
      and [-PI,+PI] respectively;
    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
     function by contrast always returns an angle in the first or fourth
     quadrants.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, double Y, X, two quantities which represent the tangent of
    an angle.  If Y is not zero, then the tangent is (Y/X).
    Output, double ATAN4, an angle between 0 and 2 * PI, whose tangent is
 (Y/X), and which lies in the appropriate quadrant so that the signs
    of its cosine and sine match those of X and Y.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp y = a_data[0];
	const register ityp x = a_data[1];
	
	/*
	Special cases:
	*/
	if ( x == 0.00 )
	{
		if ( 0.00 < y )
		{
			result = ( M_PI / 2.00 );
			return &result;
		}
		else if ( y < 0.00 )
		{
			result = ( 3.00 * M_PI / 2.00 );
			return &result;
		}
		else if ( y == 0.00 )
		{
			result = ( 0.00 );
			return &result;
		}
	} 
	else if ( y == 0.00 )
	{
		if ( 0.00 < x )
		{
			result = 0.00;
			return &result;
		}
		else if ( x < 0.00 )
		{
			result = M_PI;
			return &result;
		}
	}
	/*
	We assume that ATAN2 is reliable when both arguments are positive.
	*/
	if  ( 0.00 < x && 0.00 < y )
	{
		result = atan2 (  y,  x );
		return &result;
	}
	else if ( x < 0.00 && 0.00 < y )
	{
		result = (           M_PI - atan2 (  y, -x ) );
		return &result;
	}
	else if ( x < 0.00 && y < 0.00 )
	{
		result = (           M_PI + atan2 ( -y, -x ) );
		return &result;
	}
	else if ( 0.00 < x && y < 0.00 )
	{
		result = ( 2.00 * M_PI - atan2 ( -y,  x ) );
		return &result;
	}
	
	result = 0.00; 
	return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _icos_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    ICOS_NUM gives "sizes" for an icosahedron.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    19 July 2007
  Author:
    John Burkardt
  Parameters:
    Output, int *POINT_NUM, the number of points.
    Output, int *EDGE_NUM, the number of edges.
    Output, int *FACE_NUM, the number of faces.
    Output, int *FACE_ORDER_MAX, the maximum order of any face.
*/
{
	int ** const a_data = data;
	int * point_num = a_data[0];
	int * edge_num = a_data[1];
	int * face_num = a_data[2];
	int * face_order_max = a_data[3];
	
	*point_num = 12;
	*edge_num = 30;
	*face_num = 20;
	*face_order_max = 3;
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_distance_xyz ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_DISTANCE_XYZ computes great circle distances on a sphere.
  Discussion:
    XYZ coordinates are used.
    We assume the points XYZ1 and XYZ2 lie on the same sphere.
    This computation is a special form of the Vincenty formula.
    It should be less sensitive to errors associated with very small 
    or very large angular separations.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    29 August 2010
  Author:
    John Burkardt
  Reference:
    "Great-circle distance",
    Wikipedia.
  Parameters:
    Input, double XYZ1[3], the coordinates of the first point.
    Input, double XYZ2[3], the coordinates of the second point.
    Output, double DIST, the great circle distance between the points.
*/
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * xyz1 = a_data[0];
	ityp * xyz2 = a_data[1];
	
	ityp dist;
	ityp lat1;
	ityp lat2;
	ityp lon1;
	ityp lon2;
	ityp r;
	ityp top;
	
	r = r8vec_norm ( 3, xyz1 );
	
	lat1 = asin ( xyz1[2] );
	lon1 = atan4 ( xyz1[1], xyz1[0] );
	
	lat2 = asin ( xyz2[2] );
	lon2 = atan4 ( xyz2[1], xyz2[0] );
	
	top = pow ( cos ( lat2 ) * sin ( lon1 - lon2 ), 2 )+ pow ( cos ( lat1 ) * sin ( lat2 ) - sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ), 2 );
	top = sqrt ( top );
	
	
	result = r * atan2 ( sqrt ( top ), sin ( lat1 ) * sin ( lat2 ) + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_cubed_ijk_to_xyz ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_CUBED_IJK_TO_XYZ: cubed sphere IJK to XYZ coordinates.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of sections into which each
    face of the cube is to be divided.
    Input, int I, J, K, indices between 0 and N.  Normally,
    at least one of the indices should have the value 0 or N.
    Output, double XYZ[3], coordinates of the point.
*/
{
	const _2dtpit2dt * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	dim_typ i = s_data->a1;
	ityp * xyz = s_data->a2;
	dim_typ j = s_data->a3;
	dim_typ k = s_data->a4;
	
    ityp xc;
    ityp xyzn;
    ityp yc;
    ityp zc;

    xc = i == 0 ? -1.00 : i == n ? 1.00 : tan ( ( ityp ) ( (i<<1) - n ) * 0.25 * M_PI / ( ityp ) ( n ) );
    yc = j == 0 ? -1.00 : j == n ? 1.00 : tan ( ( ityp ) ( (j<<1) - n ) * 0.25 * M_PI / ( ityp ) ( n ) );
    zc = k == 0 ? -1.00 : k == n ? 1.00 : tan ( ( ityp ) ( (k<<1) - n ) * 0.25 * M_PI / ( ityp ) ( n ) );
    xyzn = sqrt ( xc * xc + yc * yc + zc * zc );

    xyz[0] = xc / xyzn;
    xyz[1] = yc / xyzn;
    xyz[2] = zc / xyzn;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_cubed_line_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_CUBED_LINE_NUM counts lines on a cubed sphere grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of sections into which each
    face of the cube is to be divided.
    Output, int LINE_NUM, the number of lines.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ line_num = 0;
    /*
    If N = 1, the corners form 12 lines.
    */
    if ( n == 1 )
    {
    	result = 12;
        return &result;
    }
    /*
    If 1 < N, each of 8 corners connects to three neighboring edges.
    */
    else
        line_num += 24;

    /*
    If 2 < N, then each of the 12 edges includes lines.
    */
    if ( 2 < n )
        line_num += 12 * ( n - 2 );
    /*
    Lines that belong to one of the six faces.
    */
    
    result = line_num + (1<n)*12 * n * ( n - 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_cubed_lines ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_CUBED_LINES computes the lines on a cubed sphere grid.
  Licensing
    This code is distributed under the GNU LGPL license.
  Modified:
    15 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of sections into which
    each face of the cube is to be divided.
    Input, int LINE_NUM, the number of lines.
    Output, double SPHERE_CUBED_LINES[3*2*LINE_NUM], distinct points
    on the unit sphere generated by a cubed sphere grid.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ line_num = a_data[1];
	
    dim_typ i, j, l = 0;
    ityp *line_data = ( ityp * ) malloc ( 6 * line_num );
    /*
    If N = 1, the corners form 12 lines.
    */
    if ( n == 1 )
    {
        sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, n, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, n, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, 0, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, 0, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, n, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, n, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, n, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, 0, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, n, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, n, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, n, n, line_data+0+1*3+l*6 );
        ++ l;
        return line_data;
    }
    /*
    If 1 < N, each of 8 corners connects to three neighboring edges.
    */
    else
    {
        sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 1, 0, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, 1, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, 0, 1, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n,   0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n-1, 0, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, 1, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, 0, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, 0, 1, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n,   n, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n-1, n, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, n,   0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, n-1, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, n, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, n, 1, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 1, n, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n,   0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, n-1, 0, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n, 0, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, n, 1, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 1, 0, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, 0, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, 1, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, 0, n,   line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, 0, n-1, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n,   0, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n-1, 0, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, 0, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, 1, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, 0, n,   line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, 0, n-1, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n,   n, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n-1, n, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, n,   n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, n-1, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, n, n, n,   line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, n, n, n-1, line_data+0+1*3+l*6 );

        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n, n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 1, n, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n,   n, line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, n-1, n, line_data+0+1*3+l*6 );
        ++ l;
        sphere_cubed_ijk_to_xyz ( n, 0, n, n,   line_data+0+0*3+l*6 );
        sphere_cubed_ijk_to_xyz ( n, 0, n, n-1, line_data+0+1*3+l*6 );
        ++ l;
    }
    /*
    If 2 < N, then each of the 12 edges includes lines.
    */
    if ( 2 < n )
    {
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, i,   0, 0, line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, i+1, 0, 0, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, n,   i, 0, line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, n, i+1, 0, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, n-i,   n, 0, line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, n-i-1, n, 0, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, 0, n-i,   0, line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, 0, n-i-1, 0, line_data+0+1*3+l*6 );
            ++ l;
        }

        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, i,   0, n, line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, i+1, 0, n, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, n,   i, n, line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, n, i+1, n, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, n-i,   n, n, line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, n-i-1, n, n, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, 0, n-i,   n, line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, 0, n-i-1, n, line_data+0+1*3+l*6 );
            ++ l;
        }

        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, 0, 0, i,   line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, 0, 0, i+1, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, n, 0, i,   line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, n, 0, i+1, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
            sphere_cubed_ijk_to_xyz ( n, n, n, i,   line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, n, n, i+1, line_data+0+1*3+l*6 );
            ++ l;
        }
        for ( i = 1; i <= n - 2; ++i )
        {
        sphere_cubed_ijk_to_xyz ( n, 0, n, i,   line_data+0+0*3+l*6 );
            sphere_cubed_ijk_to_xyz ( n, 0, n, i+1, line_data+0+1*3+l*6 );
            ++ l;
        }
    }
    /*
    Lines that belong to one of the six faces.
    */
    if ( 1 < n )
    {
        /*
        000, nn0
        */
        for ( i = 1; i <= n - 1; ++i )
        {
            for ( j = 0; j <= n - 1; ++j )
            {
                sphere_cubed_ijk_to_xyz ( n, i, j,   0, line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, i, j+1, 0, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        for ( j = 1; j <= n - 1; ++j )
        {
            for ( i = 0; i <= n - 1; ++i )
            {
                sphere_cubed_ijk_to_xyz ( n, i,   j, 0, line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, i+1, j, 0, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        /*
        00n, nnn
        */
        for ( i = 1; i <= n - 1; ++i )
        {
            for ( j = 0; j <= n - 1; ++j )
            {
                sphere_cubed_ijk_to_xyz ( n, i, j,   n, line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, i, j+1, n, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        for ( j = 1; j <= n - 1; ++j )
        {
            for ( i = 0; i <= n - 1; ++i )
            {
                sphere_cubed_ijk_to_xyz ( n, i,   j, n, line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, i+1, j, n, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
    /*
    000:n0n
    */
        for ( i = 1; i <= n - 1; ++i )
        {
            for ( j = 0; j <= n - 1; ++j )
            {
                sphere_cubed_ijk_to_xyz ( n, i, 0, j,   line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, i, 0, j+1, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        for ( j = 1; j <= n - 1; ++j )
        {
            for ( i = 0; i <= n - 1; ++i )
            {
                sphere_cubed_ijk_to_xyz ( n, i,   0, j, line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, i+1, 0, j, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
    /*
        0n0:nnn
        */
        for ( i = 1; i <= n - 1; ++i )
        {
            for ( j = 0; j <= n - 1; ++j )
            {
                sphere_cubed_ijk_to_xyz ( n, i, n, j,   line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, i, n, j+1, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        for ( j = 1; j <= n - 1; ++j )
        {
            for ( i = 0; i <= n - 1; ++i )
            {
                sphere_cubed_ijk_to_xyz ( n, i,   n, j, line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, i+1, n, j, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        /*
        000:0nn
        */
        for ( i = 1; i <= n - 1; ++i )
        {
            for ( j = 0; j <= n - 1; ++j )
            {
                sphere_cubed_ijk_to_xyz ( n, 0, i, j,   line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, 0, i, j+1, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        for ( j = 1; j <= n - 1; ++j )
        {
            for ( i = 0; i <= n - 1; ++i )
            {
                sphere_cubed_ijk_to_xyz ( n, 0, i,   j, line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, 0, i+1, j, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        /*
        n00:nnn
        */
        for ( i = 1; i <= n - 1; ++i )
        {
            for ( j = 0; j <= n - 1; ++j )
            {
                sphere_cubed_ijk_to_xyz ( n, n, i, j,   line_data+0+0*3+l*6 );
                sphere_cubed_ijk_to_xyz ( n, n, i, j+1, line_data+0+1*3+l*6 );
                ++ l;
            }
        }
        for ( j = 1; j <= n - 1; ++j )
        {
            for ( i = 0; i <= n - 1; ++i )
            {
                    sphere_cubed_ijk_to_xyz ( n, n, i,   j, line_data+0+0*3+l*6 );
                    sphere_cubed_ijk_to_xyz ( n, n, i+1, j, line_data+0+1*3+l*6 );
                    ++ l;
            }
        }

    }

    if ( l != line_num )
        return NULL;

    return line_data;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_cubed_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_CUBED_POINTS computes the points on a cubed sphere grid.
  Discussion:
    For a value of N = 3, for instance, each of the 6 cube faces will
    be divided into 3 sections, so that a single cube face will have
 (3+1)x(3+1) points:
      X---X---X---X
      | 1 | 4 | 7 |
      X---X---X---X
      | 2 | 5 | 8 |
      X---X---X---X
      | 3 | 6 | 9 |
      X---X---X---X
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of sections into which each
    face of the cube is to be divided.
    Input, int NS, the number of points.
    Output, double SPHERE_CUBED_POINTS[3*NS], distinct points on the
    unit sphere generated by a cubed sphere grid.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ ns = a_data[1];
	
    dim_typ ns2 = 0;
    double *xyz =( ityp * ) malloc ( 3 * ns * sizeof ( ityp ) );
    /*
    Bottom full.
    */
    sphere_cubed_points_face ( n, 0, 0, 0, n, n, 0, &ns2, xyz );
    /*
    To avoid repetition, draw the middles as grids of n-2 x n-1 points.
    */
    sphere_cubed_points_face ( n, 0, 0, 1, 0,   n-1, n-1, &ns2, xyz );
    sphere_cubed_points_face ( n, 0, n, 1, n-1, n,   n-1, &ns2, xyz );
    sphere_cubed_points_face ( n, n, 1, 1, n,   n,   n-1, &ns2, xyz );
    sphere_cubed_points_face ( n, 1, 0, 1, n,   0,   n-1, &ns2, xyz );
    /*
    Top full.
    */
    sphere_cubed_points_face ( n, 0, 0, n, n, n, n, &ns2, xyz );
    return ns2!=ns ? NULL:xyz;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_cubed_points_face ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_CUBED_POINTS_FACE: points on one face of a cubed sphere grid.
  Discussion:
    This routine starts with NS = 0, and is called repeatedly to
    add points for another face.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of sections into which each face
    of the cube is to be divided.
    Input, int I1, J1, K1, I2, J2, K2, the logical indices,
    between 0 and N, of two corners of the face grid.  It is guaranteed that
    I1 <= I2, J1 <= J2, and K1 <= K2.
    Input/output, int *NS, the number of points.
    Input/output, double XYZ[3*NS], distinct points on the unit sphere
    generated by a cubed sphere grid.
*/
{
	const _7dtpdtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	dim_typ i1 = s_data->a1;
	dim_typ j1 = s_data->a2;
	dim_typ k1 = s_data->a3;
	dim_typ i2 = s_data->a4;
	dim_typ j2 = s_data->a5;
	dim_typ k2 = s_data->a6;
	dim_typ * ns = s_data->a7;
	ityp * xyz = s_data->a8;
	
    dim_typ i, j, k;
    ityp xyzn;
    ityp xc;
    ityp yc;
    ityp zc;

    for ( i = i1; i <= i2; ++i )
    {
        if ( i1 < i2 )
            xc = tan ( ( ityp ) ( (i<<1) - n ) * 0.25 * M_PI / ( ityp ) ( n ) );
        else if ( i1 == 0 )
            xc = -1.00;
        else if ( i1 == n )
            xc = +1.00;
        else
            xc = 0.00;

        for ( j = j1; j <= j2; ++j )
        {
            if ( j1 < j2 )
                yc = tan ( ( ityp ) ( (j<<1) - n ) * 0.25 * M_PI / ( ityp ) ( n ) );
            else if ( j1 == 0 )
                yc = -1.00;
            else if ( j1 == n )
                yc = +1.00;
            else
                yc = 0.00;

            for ( k = k1; k <= k2; ++k )
            {
                if ( k1 < k2 )
                    zc = tan ( ( ityp ) ( (k<<1) - n ) * 0.25 * M_PI / ( ityp ) ( n ) );
                else if ( k1 == 0 )
                    zc = -1.00;
                else if ( k1 == n )
                    zc = +1.00;
                else
                    zc = 0.00;

                xyzn = sqrt ( xc * xc + yc * yc + zc * zc );

                xyz[0+ *ns *3] = xc / xyzn;
                xyz[1+ *ns *3] = yc / xyzn;
                xyz[2+ *ns *3] = zc / xyzn;
                ++ *ns;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_cubed_point_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_CUBED_POINT_NUM counts the points on a cubed sphere grid.
  Discussion:
    For a value of N = 3, for instance, each of the 6 cube faces will
    be divided into 3 sections, so that a single cube face will have
 (3+1)x(3+1) points:
      X---X---X---X
      | 1 | 4 | 7 |
      X---X---X---X
      | 2 | 5 | 8 |
      X---X---X---X
      | 3 | 6 | 9 |
      X---X---X---X
    The number of points is simply (N+1)^3 - (N-1)^3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of sections into which
    each face of the cube is to be divided.
    Output, int SPHERE_CUBED_POINT_NUM, the number of points.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = powi ( n + 1, 3 ) - powi ( n - 1, 3 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_grid_q4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_GRID_Q4: rectangular grid on a sphere.
  Discussion:
    The point numbering system is the same used in SPHERE_GRIDPOINTS,
    and that routine may be used to compute the coordinates of the points.
    A sphere in 3D satisfies the equation:
      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R * R
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int LAT_NUM, the number of "rows" of rectangles to
    be created.  LAT_NUM must be at least 2.
    Input, int LONG_NUM, the number of "columns" of
    rectangles to be created.
    Output, int SPHERE_GRID_Q4[4*LAT_NUM*LONG_NUM],
    the indices of the nodes that make up each rectangle.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ lat_num = a_data[0];
	const register dim_typ long_num = a_data[1];
	
    dim_typ i;
    dim_typ j;
    dim_typ n;
    dim_typ n_max;
    dim_typ n_min;
    dim_typ ne;
    dim_typ nw;
    dim_typ s;
    dim_typ s_max;
    dim_typ s_min;
    dim_typ se;
    dim_typ sw;
    int *rectangle_node;
    dim_typ rectangle_num;

    rectangle_node = ( int * )
    malloc ( ( lat_num * long_num << 2 ) * sizeof ( int ) );

    rectangle_num = n = 0;
    /*
    The first row.
    */

    sw = 1;
    se = sw + 1;

    s_min = 1;
    s_max = long_num;

    for ( j = 1; j <= long_num; ++j )
    {
        rectangle_node[0+(rectangle_num<<2)] = sw;
        rectangle_node[1+(rectangle_num<<2)] = se;
        rectangle_node[2+(rectangle_num<<2)] = n;
        rectangle_node[3+(rectangle_num<<2)] = n;
        rectangle_num = rectangle_num + 1;

        sw = se = se == s_max ? s_min:se+1;
    }
    /*
    The intermediate rows.
    */
    for ( i = 2; i < lat_num; ++i )
    {
        n_max = s_max;
        n_min = s_min;

        s_max = s_max + long_num;
        s_min = s_min + long_num;

        nw = n_min;
        ne = nw + 1;
        sw = s_min;
        se = sw + 1;

        for ( j = 1; j <= long_num; ++j )
        {

            rectangle_node[0+(rectangle_num<<2)] = sw;
            rectangle_node[1+(rectangle_num<<2)] = se;
            rectangle_node[2+(rectangle_num<<2)] = ne;
            rectangle_node[3+(rectangle_num<<2)] = nw;
            rectangle_num = rectangle_num + 1;

            sw = se;
            nw = ne = ne == n_max ? n_min:ne+1;
            se == se == s_max ? s_min:se+1;
        }
    }
    /*
    The last row.
    */
    n_max = s_max;
    n_min = s_min;

    s = n_max + 1;

    nw = n_min;
    ne = nw + 1;

    for ( j = 1; j <= long_num; ++j )
    {
        rectangle_node[0+rectangle_num*4] = ne;
        rectangle_node[1+rectangle_num*4] = nw;
        rectangle_node[2+rectangle_num*4] = s;
        rectangle_node[3+rectangle_num*4] = s;
        rectangle_num = rectangle_num + 1;

        nw = ne = ne == n_max ? n_min:ne+1;
    }
    return rectangle_node;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_grid_t3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_GRID_T3 produces a triangle grid on a sphere.
  Discussion:
    The point numbering system is the same used in SPHERE_GRIDPOINTS,
    and that routine may be used to compute the coordinates of the points.
    A sphere in 3D satisfies the equation:
      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R * R
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int LAT_NUM, LONG_NUM, the number of latitude
    and longitude lines to draw.  The latitudes do not include the North
    and South poles, which will be included automatically, so LAT_NUM = 5,
    for instance, will result in points along 7 lines of latitude.
    Output, int SPHERE_GRID_T3[3*2*(LAT_NUM+1)*LONG_NUM], the
    triangle vertices.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ lat_num = a_data[0];
	const register dim_typ long_num = a_data[1];
	
    dim_typ i;
    dim_typ j;
    dim_typ n;
    dim_typ n_max;
    dim_typ n_min;
    dim_typ ne;
    dim_typ nw;
    dim_typ s;
    dim_typ s_max;
    dim_typ s_min;
    dim_typ se;
    dim_typ sw;
    int *triangle_node;
    dim_typ triangle_num;

    triangle_node = ( int * )
    malloc ( 6 * ( lat_num + 1 ) * long_num * sizeof ( int ) );

    triangle_num = 0;
    /*
    The first row.
    */
    /*
    Not working!:

    n = 0;
    sw = 1;
    se = sw + 1;
    s_min = 1;
    s_max = long_num;
    */
    n = 1;
    sw = s_min = 2;
    se = sw + 1;
    s_max = long_num + 1;

    for ( j = 0; j < long_num; ++j)
    {
        triangle_node[0+triangle_num*3] = sw - 1;
        triangle_node[1+triangle_num*3] = se - 1;
        triangle_node[2+triangle_num*3] = n - 1;
        ++ triangle_num;

        sw = se = se == s_max?s_min:se+1;
    }
    /*
    The intermediate rows.
    */
    for ( i = 1; i <= lat_num; ++i )
    {
        n_max = s_max;
        n_min = s_min;

        s_max = s_max + long_num;
        s_min = s_min + long_num;

        nw = n_min;
        ne = nw + 1;
        sw = s_min;
        se = sw + 1;

        for ( j = 0; j < long_num; ++j )
        {
            triangle_node[0+triangle_num*3] = sw - 1;
            triangle_node[1+triangle_num*3] = se - 1;
            triangle_node[2+triangle_num*3] = nw - 1;
            triangle_num = triangle_num + 1;

            triangle_node[0+triangle_num*3] = ne - 1;
            triangle_node[1+triangle_num*3] = nw - 1;
            triangle_node[2+triangle_num*3] = se - 1;
            ++ triangle_num;

            sw = se = se == s_max?s_min:ne+1;
            nw = ne = ne == n_max?n_min:ne+1;
        }
    }
    /*
    The last row.
    */
    n_max = s_max;
    n_min = s_min;

    s = n_max + 1;

    nw = n_min;
    ne = nw + 1;

    for ( j = 0; j < long_num; ++j )
    {
        triangle_node[0+triangle_num*3] = ne - 1;
        triangle_node[1+triangle_num*3] = nw - 1;
        triangle_node[2+triangle_num*3] = s - 1;
        ++ triangle_num;

        nw = ne = ne == n_max?n_min:ne+1;
    }
    return triangle_node;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_icos_edge_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_ICOS_EDGE_NUM sizes an icosahedral grid on a sphere.
  Discussion:
    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
    With FACTOR = 2, each triangle of the icosahedron is subdivided into
    2x2 subtriangles, resulting in 80 faces, 120 edges, and
    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
    With FACTOR = 3, each triangle of the icosahedron is subdivided into
    3x3 subtriangles, resulting in 180 faces, 270 edges and
    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
      12
    + 20 * 3          * (FACTOR-1) / 2
    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.
    Output, int SPHERE_ICOS_EDGE_NUM, the number of edges.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ factor = *(dim_typ *) data;
	
	result = 30 * factor * factor;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_icos_face_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_ICOS_FACE_NUM sizes an icosahedral grid on a sphere.
  Discussion:
    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
    With FACTOR = 2, each triangle of the icosahedron is subdivided into
    2x2 subtriangles, resulting in 80 faces, 120 edges, and
    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
    With FACTOR = 3, each triangle of the icosahedron is subdivided into
    3x3 subtriangles, resulting in 180 faces, 270 edges and
    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
      12
    + 20 * 3          * (FACTOR-1) / 2
    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.
    Output, int SPHERE_ICOS_FACE_NUM, the number of triangles.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ factor = *(dim_typ *) data;
	
	result = 20 * factor * factor;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_icos_point_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_ICOS_POINT_NUM sizes an icosahedral grid on a sphere.
  Discussion:
    With FACTOR = 1, the grid has 20 triangular faces, 30 edges, and 12 nodes.
    With FACTOR = 2, each triangle of the icosahedron is subdivided into
    2x2 subtriangles, resulting in 80 faces, 120 edges, and
    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
    With FACTOR = 3, each triangle of the icosahedron is subdivided into
    3x3 subtriangles, resulting in 180 faces, 270 edges and
    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
    resulting in 20 * FACTOR * FACTOR faces, 30 * FACTOR * FACTOR edges, and
      12
    + 20 * 3          * (FACTOR-1) / 2
    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.

    Output, int SPHERE_ICOS_NODE_NUM, the number of nodes.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ factor = *(dim_typ *) data;
	
	result = 12 + 10 * 3              * ( factor - 1 )+ 10 * ( factor - 2 ) * ( factor - 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_icos1_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_ICOS1_POINTS returns icosahedral grid points on a sphere.
  Discussion:
    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
    With FACTOR = 2, each triangle of the icosahedron is subdivided into
    2x2 subtriangles, resulting in 80 faces and
    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
    With FACTOR = 3, each triangle of the icosahedron is subdivided into
    3x3 subtriangles, resulting in 180 faces and
    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
    resulting in 20 * FACTOR * FACTOR faces and
      12
    + 20 * 3          * (FACTOR-1) / 2
    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 July 2007
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.
    Input, int NODE_NUM, the number of nodes, as reported
    by SPHERE_GRID_ICOS_NUM.
    Output, double SPHERE_GRIDPOINTS_ICOS1[3*NODE_NUM], the node coordinates.
  Local Parameters:
    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters
    associated with the icosahedron, and POINT_COORD, EDGE_POINT,
    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
    We need to refer to this data to generate the grid.
    NODE counts the number of nodes we have generated so far.  At the
    end of the routine, it should be equal to NODE_NUM.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ factor = a_data[0];
	const register dim_typ node_num = a_data[1];
	
    int a;
    int b;
    int c;
    dim_typ dim;
    int edge;
    int edge_num;
    int *edge_point;
    int f;
    dim_typ f1;
    dim_typ f2;
    int face;
    int face_num;
    int *face_order;
    int *face_point;
    int face_order_max;
    dim_typ node;
    ityp node_norm;
    ityp *node_xyz;
    dim_typ point;
    ityp *point_coord;
    int point_num;

    node_xyz = ( ityp * ) malloc ( 3 * node_num * sizeof ( ityp ) );
    /*
    Size the icosahedron.
    */
    icos_num ( &point_num, &edge_num, &face_num, &face_order_max );
    /*
    Set the icosahedron.
    */
    point_coord = ( ityp * ) malloc ( 3 * point_num * sizeof ( ityp ) );
    edge_point = ( int * ) malloc ( edge_num * sizeof ( int ) <<1 );
    face_order = ( int * ) malloc ( face_num * sizeof ( int ) );
    face_point = ( int * ) malloc ( face_order_max * face_num * sizeof ( int ) );

    icos_shape ( point_num, edge_num, face_num, face_order_max,
    point_coord, edge_point, face_order, face_point );
    /*
    Generate the point coordinates.

    A.  Points that are the icosahedral vertices.
    */
    node = 0;
    for ( point = 0; point < point_num; ++point )
    {
        #pragma omp parallel for num_threads(3)
        for ( dim = 0; dim < 3; ++dim)
            node_xyz[dim+node*3] = point_coord[dim+point*3];
        ++ node ;
    }
    /*
    B. Points in the icosahedral edges, at
    1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
    */
    for ( edge = 0; edge < edge_num; ++edge )
    {
        a = edge_point[0+edge*2];
        b = edge_point[1+edge*2];

        for ( f = 1; f < factor; ++f )
        {
            #pragma omp parallel for num_threads(3)
            for ( dim = 0; dim < 3; ++dim )
                node_xyz[dim+node*3] =( ( ityp ) ( factor - f ) * point_coord[dim+a*3]+ ( ityp ) (          f ) * point_coord[dim+b*3] )/ ( ityp ) ( factor     );

            node_norm = r8vec_norm ( 3, node_xyz+node*3 );

            #pragma omp parallel for num_threads(3)
            for ( dim = 0; dim < 3; ++dim )
                node_xyz[dim+node*3] /= node_norm;
            ++ node;
        }
    }
    /*
    C.  Points in the icosahedral faces.
    */
    for ( face = 0; face < face_num; ++face )
    {
        a = face_point[0+face*3];
        b = face_point[1+face*3];
        c = face_point[2+face*3];

        for ( f1 = 1; f1 < factor; ++f1 )
        {
            for ( f2 = 1; f2 < factor - f1; ++f2 )
            {
                #pragma omp parallel for num_threads(3)
                for ( dim = 0; dim < 3; ++dim )
                    node_xyz[dim+node*3] =( ( ityp ) ( factor - f1 - f2 ) * point_coord[dim+a*3]+ ( ityp ) (          f1      ) * point_coord[dim+b*3]+ ( ityp ) (               f2 ) * point_coord[dim+c*3] )/ ( ityp ) ( factor           );
                node_norm = r8vec_norm ( 3, node_xyz+node*3 );

                #pragma omp parallel for num_threads(3)
                for ( dim = 0; dim < 3; ++dim )
                node_xyz[dim+node*3] /= node_norm;
                ++ node;
            }
        }
    }
    /*
    Discard allocated memory.
    */
    free ( edge_point );
    free ( face_order );
    free ( face_point );
    free ( point_coord );

    return node_xyz;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_icos2_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_ICOS2_POINTS returns icosahedral grid points on a sphere.
  Discussion:
    With FACTOR = 1, the grid has 20 triangular faces and 12 nodes.
    With FACTOR = 2, each triangle of the icosahedron is subdivided into
    2x2 subtriangles, resulting in 80 faces and
    42 = 12 + 20 * 3 * (1)/2 + 20 * 0 ) nodes.
    With FACTOR = 3, each triangle of the icosahedron is subdivided into
    3x3 subtriangles, resulting in 180 faces and
    92 ( = 12 + 20 * 3 * (2)/2 + 20 * 1 ) nodes.
    In general, each triangle is subdivided into FACTOR*FACTOR subtriangles,
    resulting in 20 * FACTOR * FACTOR faces and
      12
    + 20 * 3          * (FACTOR-1) / 2
    + 20 * (FACTOR-2) * (FACTOR-1) / 2 nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.
    Input, int NODE_NUM, the number of nodes, as reported
    by SPHERE_IMP_GRID_ICOS_NUM.
    Output, double SPHERE_GRIDPOINTS_ICOS2[3*NODE_NUM], the node coordinates.
  Local Parameters:
    POINT_NUM, EDGE_NUM, FACE_NUM and FACE_ORDER_MAX are counters
    associated with the icosahedron, and POINT_COORD, EDGE_POINT,
    FACE_ORDER and FACE_POINT are data associated with the icosahedron.
    We need to refer to this data to generate the grid.
    NODE counts the number of nodes we have generated so far.  At the
    end of the routine, it should be equal to NODE_NUM.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ factor = a_data[0];
	const register dim_typ node_num = a_data[1];
	
    int a;
    ityp angle;
    ityp ab[3];
    ityp ac[3];
    ityp acn[3];
    ityp acn_norm;
    ityp acp[3];
    int b;
    ityp bn[3];
    ityp bn_norm;
    ityp bp[3];
    int c;
    ityp cn[3];
    ityp cn_norm;
    ityp cp[3];
    int edge;
    int edge_num;
    int *edge_point;
    int f;
    int fa;
    int fbc;
    dim_typ face;
    int face_num;
    int *face_order;
    int *face_point;
    int face_order_max;
    dim_typ i;
    dim_typ j;
    dim_typ node;
    ityp *node_xyz;
    ityp p;
    dim_typ point;
    ityp *point_coord;
    int point_num;
    ityp theta;
    ityp theta_ab;
    ityp theta_ac;
    ityp theta_bc;

    node_xyz = ( ityp * ) malloc ( 3 * node_num * sizeof ( ityp ) );
    /*
    Size the icosahedron.
    */
    icos_num ( &point_num, &edge_num, &face_num, &face_order_max );
    /*
    Set the icosahedron.
    */
    point_coord = ( ityp * ) malloc ( 3 * point_num * sizeof ( ityp ) );
    edge_point = ( int * ) malloc ( edge_num * sizeof ( int ) << 1 );
    face_order = ( int * ) malloc ( face_num * sizeof ( int ) );
    face_point = ( int * ) malloc ( face_order_max * face_num * sizeof ( int ) );

    icos_shape ( point_num, edge_num, face_num, face_order_max,
    point_coord, edge_point, face_order, face_point );
    /*
    Generate the point coordinates.

    A.  Points that are the icosahedral vertices.
    */
    node = 0;
    for ( j = 0; j < point_num; ++j )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            node_xyz[i+j*3] = point_coord[i+j*3];
        ++ node;
    }
    /*
    B. Points in the icosahedral edges, at
    1/FACTOR, 2/FACTOR, ..., (FACTOR-1)/FACTOR.
    */
    for ( edge = 0; edge < edge_num; ++edge )
    {
        a = edge_point[0+edge*2];
        b = edge_point[1+edge*2];
        /*
        Determine the "distance" = angle between points A and B.
        */
        theta = sphere_distance_xyz ( point_coord+a*3, point_coord+b*3 );
        /*
        Polarize B into BP + BN and normalize BN.
        */
        r8vec_polarize ( 3, point_coord+b*3, point_coord+a*3, bn, bp );
        bn_norm = r8vec_norm ( 3, bn );
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            bn[i] /= bn_norm;
        /*
        March from A to B, by taking equally spaced angles from 0 to THETA.
        F = 0      => ANGLE = 0     => A
        F = FACTOR => ANGLE = THETA => B
        */
        for ( f = 1; f < factor; ++f )
        {
            angle = ( ( ityp ) ( f ) * theta ) / ( ityp ) ( factor );
            #pragma omp parallel for num_threads(3)
            for ( i = 0; i < 3; ++i )
                node_xyz[i+node*3] = cos ( angle ) * point_coord[i+a*3]+ sin ( angle ) * bn[i];
            ++ node;
        }
    }
    /*
    C.  Points in the icosahedral faces.
    */
    for ( face = 0; face < face_num; ++face )
    {
        a = face_point[0+face*3];
        b = face_point[1+face*3];
        c = face_point[2+face*3];
        /*
        Determine the "distance" = angle between points A and B, A and C.
        */
        theta_ab = sphere_distance_xyz ( point_coord+a*3, point_coord+b*3 );
        theta_ac = sphere_distance_xyz ( point_coord+a*3, point_coord+c*3 );
        /*
        Polarize B = BP + BN and normalize BN, C = CP + CN, and normalize CN.
        */
        r8vec_polarize ( 3, point_coord+b*3, point_coord+a*3, bn, bp );
        bn_norm = r8vec_norm ( 3, bn );
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            bn[i] /= bn_norm;

        r8vec_polarize ( 3, point_coord+c*3, point_coord+a*3, cn, cp );
        cn_norm = r8vec_norm ( 3, cn );
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            cn[i] /= cn_norm;
        /*
        March AB from A to B:
        FA = 0      => ANGLE = 0        => AB = A
        FA = FACTOR => ANGLE = THETA_AB => AB = B

        March AC from A to C:
        FA = 0      => ANGLE = 0        => AC = A
        FA = FACTOR => ANGLE = THETA_AC => AC = C
        */
        for ( fa = 2; fa < factor; ++fa )
        {
            /*
            Determine points AB and AC that use cos ( FA / FACTOR ) of A
            and cos ( ( FACTOR - FA ) / FACTOR ) of B or C.
            */
            angle = ( (ityp ) ( fa ) * theta_ab ) / ( ityp ) ( factor );
            #pragma omp parallel for num_threads(3)
            for ( i = 0; i < 3; ++i )
                ab[i] = cos ( angle ) * point_coord[i+a*3] + sin ( angle ) * bn[i];

            angle = ( ( ityp ) ( fa ) * theta_ac ) / ( ityp ) ( factor );
            #pragma omp parallel for num_threads(3)
            for ( i = 0; i < 3; ++i )
                ac[i] = cos ( angle ) * point_coord[i+a*3] + sin ( angle ) * cn[i];
            /*
            Determine the "distance" = angle between points AB and AC.
            */
            theta_bc = sphere_distance_xyz ( ab, ac );
            /*
            Polarize AC into ACP + ACN and normalize ACN.
            */
            r8vec_polarize ( 3, ac, ab, acn, acp );
            acn_norm = r8vec_norm ( 3, acn );
            #pragma omp parallel for num_threads(3)
            for ( i = 0; i < 3; ++i )
                acn[i] /= acn_norm;
            /*
            The interval between AB and AC is broken into FA intervals.
            Go from 1 to FA - 1.
            */
            for ( fbc = 1; fbc < fa; ++fbc)
            {
                angle = ( ( ityp ) ( fbc ) * theta_bc ) / ( ityp ) ( fa );
                #pragma omp parallel for num_threads(3)
                for ( i = 0; i < 3; ++i )
                    node_xyz[i+node*3] = cos ( angle ) * ab[i] + sin ( angle ) * acn[i];
                ++ node;
            }
        }
    }
    /*
    Discard allocated memory.
    */
    free ( edge_point );
    free ( face_order );
    free ( face_point );
    free ( point_coord );

    return node_xyz;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_line_project ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LINE_PROJECT projects a line onto an implicit sphere.
  Discussion:
    The implicit form of a sphere in 3D is:
        pow ( P[0] - PC[0], 2 )
      + pow ( P[1] - PC[1], 2 )
      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
    The line to be projected is specified as a sequence of points.
    If two successive points subtend a small angle, then the second
    point is essentially dropped.  If two successive points subtend
    a large angle, then intermediate points are inserted, so that
    the projected line stays closer to the sphere.
    Note that if any P coincides with the center of the sphere, then
    its projection is mathematically undefined.  P will
    be returned as PC.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2005
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the sphere.  If R is
    zero, PP will be returned as PC, and if R is
    negative, points will end up diametrically opposite from where
    you would expect them for a positive R.
    Input, double PC[3], the coordinates of the center of the sphere.
    Input, int N, the number of points on the line that is
    to be projected.
    Input, double P[3*N], the coordinates of the points
    on the line that is to be projected.
    Input, int MAXPNT2, the maximum number of points on the projected
    line.  Even if the routine thinks that more points are needed,
    no more than MAXPNT2 will be generated.
    Output, double PP[3*N2], the coordinates of the
    points representing the projected line.  The value N2 is returned
    as the function value of this routine.  These points lie on the
    sphere.  Successive points are separated by at least THETAMIN
    radians, and by no more than THETAMAX radians.
    Input, double THETAMIN, THETAMAX, the minimum and maximum angular
    projections allowed between successive projected points.
    If two successive points on the original line have projections
    separated by more than THETAMAX radians, then intermediate points
    will be inserted, in an attempt to keep the line closer to the
    sphere.  If two successive points are separated by less than
    THETAMIN radians, then the second point is dropped, and the
    line from the first point to the next point is considered.
    Output, int SPHERE_LINE_PROJECT, the number of points on
    the projected line.  This value can be zero, if the line has an
    angular projection of less than THETAMIN radians.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const itpitdtpitdtpit2it * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * pc = s_data->a1;
	dim_typ n = s_data->a2;
	ityp * p = s_data->a3;
	dim_typ maxpnt2 = s_data->a4;
	ityp * pp = s_data->a5;
	const register ityp thetamin = s_data->a6;
	const register ityp thetamax = s_data->a7;
	
    # define DIM_NUM 3

    ityp alpha;
    ityp ang3d;
    ityp dot;
    dim_typ i;
    dim_typ j;
    dim_typ nfill;
    dim_typ n2;
    ityp tnorm;
    ityp p1[DIM_NUM];
    ityp p2[DIM_NUM];
    ityp pi[DIM_NUM];
    /*
    Check the input.
    */
    if ( r == 0.00 )
    {
    	result = 0;
        return &result;
    }
    r8vec_copy ( DIM_NUM, pc, p1 );
    r8vec_copy ( DIM_NUM, pc, p2 );

    n2 = 0;

    for ( i = 0; i < n; ++i)
    {
        if ( r8vec_eq ( DIM_NUM, p, pc ) );
        else
        {
            r8vec_copy ( DIM_NUM, p2, p1 );

            alpha = sqrt ( pow ( p[0+i*3] - pc[0], 2 )+ pow ( p[1+i*3] - pc[1], 2 )+ pow ( p[2+i*3] - pc[2], 2 ) );

            p2[0] = pc[0] + r * ( p[0+i*3] - pc[0] ) / alpha;
            p2[1] = pc[1] + r * ( p[1+i*3] - pc[1] ) / alpha;
            p2[2] = pc[2] + r * ( p[2+i*3] - pc[2] ) / alpha;
            /*
            If we haven't gotten any points yet, take this point as our start.
            /*/
            if ( n2 == 0 )
            {
                pp[0+n2*3] = p2[0];
                pp[1+n2*3] = p2[1];
                pp[2+n2*3] = p2[2];
                ++ n2;
            }
            /*
            Compute the angular projection of P1 to P2.
            */
            else if ( 1 <= n2 )
            {
                dot = ( p1[0] - pc[0] ) * ( p2[0] - pc[0] )+ ( p1[1] - pc[1] ) * ( p2[1] - pc[1] )+ ( p1[2] - pc[2] ) * ( p2[2] - pc[2] );
                ang3d = acos (  dot / ( r * r ) );
                /*
                If the angle is at least THETAMIN, (or it's the last point),
                then we will draw a line segment.
                */
                if ( thetamin < abs ( ang3d ) || i == n )
                {
                    /*
                    Now we check to see if the line segment is too long.
                    */
                    if ( thetamax < abs ( ang3d ) )
                    {
                        nfill = ( dim_typ ) ( abs ( ang3d ) / thetamax );

                        for ( j = 1; j < nfill; ++j )
                        {
                            pi[0] = ( ( ityp ) ( nfill - j ) * ( p1[0] - pc[0] )+ ( ityp ) (         j ) * ( p2[0] - pc[0] ) );
                            pi[1] = ( ( ityp ) ( nfill - j ) * ( p1[1] - pc[1] )+ ( ityp ) (         j ) * ( p2[1] - pc[1] ) );
                            pi[2] = ( ( ityp ) ( nfill - j ) * ( p1[2] - pc[2] )+ ( ityp ) (         j ) * ( p2[2] - pc[2] ) );

                            tnorm = r8vec_norm ( DIM_NUM, pi );

                            if ( tnorm != 0.00 )
                            {
                                pi[0] = pc[0] + r * pi[0] / tnorm;
                                pi[1] = pc[1] + r * pi[1] / tnorm;
                                pi[2] = pc[2] + r * pi[2] / tnorm;
                                pp[0+n2*3] = pi[0];
                                pp[1+n2*3] = pi[1];
                                pp[2+n2*3] = pi[2];
                                ++ n2;
                            }
                        }
                    }
                    /*
                    Now tack on the projection of point 2.
                    */
                    pp[0+n2*3] = p2[0];
                    pp[1+n2*3] = p2[1];
                    pp[2+n2*3] = p2[2];
                    ++ n2;
                }
            }
        }
    }
    
    result = n2;
    return &result;
    # undef DIM_NUM
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_ll_lines ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LL_LINES produces lines on a latitude/longitude grid.
  Discussion:
    The point numbering system is the same used in SPHERE_LL_POINTS,
    and that routine may be used to compute the coordinates of the points.
    The implicit form of a sphere in 3D is:
        pow ( P[0] - PC[0], 2 )
      + pow ( P[1] - PC[1], 2 )
      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int NLAT, NLONG, the number of latitude and longitude
    lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so NLAT = 5, for instance,
    will result in points along 7 lines of latitude.
    Input, int LINE_NUM, the number of grid lines.
    Output, int SPHERE_LL_LINES[2*LINE_NUM], contains pairs of point indices for
    line segments that make up the grid.
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
        old = 1;
        next = j + 2;
        line[0+(l<<1)] = old;
        line[1+(l<<1)] = next;
        l = l + 1;

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
        line[1+(l<<1)] = 1 + nlat * nlong + 1;
        ++ l;
    }
    /*
    "Horizontal" lines.
    */
    for ( i = 1; i <= nlat; ++i )
    {
        next = 1 + ( i - 1 ) * nlong + 1;

        for ( j = 0; j <= nlong-2; ++j )
        {
            old = next;
            next = old + 1;
            line[0+(l<<1)] = old;
            line[1+(l<<1)] = next;
            ++ l;
        }

        old = next;
        next = 1 + ( i - 1 ) * nlong + 1;
        line[0+(l<<1)] = old;
        line[1+(l<<1)] = next;
        ++ l;
    }
    /*
    "Diagonal" lines.
    */
    for ( j = 0; j <= nlong - 1; ++j )
    {
        old = 1;
        next = j + 2;
        newcol = j;

        for ( i = 1; i <= nlat - 1; ++i )
        {
            old = next;
            next = old + nlong + 1;

            ++ newcol;
            if ( nlong - 1 < newcol )
            {
                newcol = 0;
                next = next - nlong;
            }

            line[0+l*2] = old;
            line[1+l*2] = next;
            ++ l;
        }
    }
    return line;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_ll_line_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LL_LINE_NUM counts lines for a latitude/longitude grid.
  Discussion:
    The number returned is the number of pairs of points to be connected.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int LAT_NUM, LONG_NUM, the number of latitude and
    longitude lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so LAT_NUM = 5, for instance,
    will result in points along 7 lines of latitude.
    Output, int SPHERE_LL_LINE_NUM, the number of grid lines.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ lat_num = a_data[0];
	const register dim_typ long_num = a_data[1];
	
	result = long_num * ( lat_num + 1 )+ lat_num * long_num+ long_num * ( lat_num - 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_ll_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LL_POINTS produces points on a latitude/longitude grid.
  Discussion:
    The implicit form of a sphere in 3D is:
        pow ( P[0] - PC[0], 2 )
      + pow ( P[1] - PC[1], 2 )
      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2012
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
    Output, double SPHERE_LL_POINTS[3*POINT_NUM], the coordinates
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
    ++ n;
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

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_ll_point_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LL_POINT_NUM counts points for a latitude/longitude grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int LAT_NUM, LONG_NUM, the number of latitude
    and longitude lines to draw.  The latitudes do not include the North and
    South poles, which will be included automatically, so LAT_NUM = 5, for
    instance, will result in points along 7 lines of latitude.
    Output, int SPHERE_LL_POINT_NUM, the number of grid points.
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
__MATHSUITE __JBURKARDT  void   * _sphere_llq_lines ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LLQ_LINES: latitude/longitude quadrilateral grid lines.
  Discussion:
    The point numbering system is the same used in SPHERE_LL_POINTS,
    and that routine may be used to compute the coordinates of the points.
    The implicit form of a sphere in 3D is:
        pow ( P[0] - PC[0], 2 )
      + pow ( P[1] - PC[1], 2 )
      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int NLAT, NLONG, the number of latitude and longitude
    lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so NLAT = 5, for instance,
    will result in points along 7 lines of latitude.
    Input, int LINE_NUM, the number of grid lines.
    Output, int SPHERE_LLQ_LINES[2*LINE_NUM], contains pairs of point indices for
    line segments that make up the grid.
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
        old = 1;
        next = j + 2;
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
        line[1+(l<<1)] = 1 + nlat * nlong + 1;
        ++ l;
    }
    /*
    "Horizontal" lines.
    */
    for ( i = 1; i <= nlat; ++i )
    {
        next = 1 + ( i - 1 ) * nlong + 1;

        for ( j = 0; j <= nlong-2; ++j )
        {
            old = next;
            next = old + 1;
            line[0+(l<<1)] = old;
            line[1+(l<<1)] = next;
            ++ l ;
        }

        old = next;
        next = 1 + ( i - 1 ) * nlong + 1;
        line[0+(l<<1)] = old;
        line[1+(l<<1)] = next;
        ++ l;
    }

    return line;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere_llq_line_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_LLQ_LINE_NUM counts lines for a latitude/longitude quadrilateral grid.
  Discussion:
    The number returned is the number of pairs of points to be connected.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int LAT_NUM, LONG_NUM, the number of latitude and
    longitude lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so LAT_NUM = 5, for instance,
    will result in points along 7 lines of latitude.
    Output, int SPHERE_LLQ_LINE_NUM, the number of grid lines.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ lat_num = a_data[0];
	const register dim_typ long_num = a_data[1];
	
	result = long_num * ( lat_num + 1 )+ lat_num * long_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_spiralpoints ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_SPIRALPOINTS produces spiral points on an implicit sphere.
  Discussion:
    The points should be arranged on the sphere in a pleasing design.
    The implicit form of a sphere in 3D is:
        pow ( P[0] - PC[0], 2 )
      + pow ( P[1] - PC[1], 2 )
      + pow ( P[2] - PC[2], 2 ) = pow ( R, 2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 July 2005
  Author:
    John Burkardt
  Reference:
    Edward Saff, Arno Kuijlaars,
    Distributing Many Points on a Sphere,
    The Mathematical Intelligencer,
    Volume 19, Number 1, 1997, pages 5-11.
  Parameters:
    Input, double R, the radius of the sphere.
    Input, double PC[3], the coordinates of the center of the sphere.
    Input, int N, the number of points to create.
    Output, double SPHERE_SPIRALPOINTS[3*N], the coordinates of the grid points.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * pc = s_data->a1;
	const register ityp r = s_data->a2;
	
    ityp cosphi = 0.00;
    dim_typ i;
    ityp *p;
    ityp sinphi = 0.0;
    ityp theta = 0.00;

    p = ( ityp * ) malloc ( 3 * n * sizeof ( ityp ) );

    for ( i = 0; i < n; i++ )
    {
        cosphi = ( ( ityp ) ( n - i - 1 ) * ( -1.0 )+ ( ityp ) (     i     ) * (  1.0 ) )/ ( ityp ) ( n     - 1 );

        sinphi = sqrt ( 1.0 - cosphi * cosphi );

        if ( i == 0 || i == n - 1 )
            theta = 0.0;
        else
        {
            theta += 3.60 / ( sinphi * sqrt ( ( ityp ) n ) );
            theta = r8_modp ( theta, M_2TPI );
        }
        p[0+i*3] = pc[0] + r * sinphi * cos ( theta );
        p[1+i*3] = pc[1] + r * sinphi * sin ( theta );
        p[2+i*3] = pc[2] + r * cosphi;
    }

    return p;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_unit_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_UNIT_SAMPLE picks a random point on the unit sphere.
  Discussion:
    The unit sphere in 3D satisfies the equation:
      X * X + Y * Y + Z * Z = 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points to generate.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double SPHERE_UNIT_SAMPLE[3*N], the sample point.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    dim_typ j;
    ityp phi;
    ityp theta;
    ityp vdot;
    ityp *x;

    x = ( ityp * ) malloc ( 3 * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
    {
        /*
        Pick a uniformly random VDOT, which must be between -1 and 1.
        This represents the dot product of the random vector with the Z unit vector.

        this works because the surface area of the sphere between
        Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
        a patch of area uniformly.
        */
        vdot = 2.00 * r8_uniform_01 ( seed ) - 1.0;

        phi = acos ( vdot );
        /*
        Pick a uniformly random rotation between 0 and 2 Pi around the
        axis of the Z vector.
        */
        theta = M_2TPI * r8_uniform_01 ( seed );

        x[0+j*3] = cos ( theta ) * sin ( phi );
        x[1+j*3] = sin ( theta ) * sin ( phi );
        x[2+j*3] = cos ( phi );
    }

    return x;
}

#endif
