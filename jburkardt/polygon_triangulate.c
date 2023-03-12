#ifndef __DISABLEDEEP_POLYGONTRIANGULATE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _between ( void * data)
/******************************************************************************/
/*
  Purpose:
    BETWEEN is TRUE if vertex C is between vertices A and B.
  Discussion:
    The points must be (numerically) collinear.
    Given that condition, we take the greater of XA - XB and YA - YB
    as a "scale" and check where C's value lies.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2014
  Author:
    Original C version by Joseph ORourke.
    This C version by John Burkardt.
  Reference:
    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
  Parameters:
    Input, double XA, YA, XB, YB, XC, YC, the coordinates of
    the vertices.
    Output, int BETWEEN, is TRUE if C is between A and B.
*/
{
	static bool result = 2;
	
	ityp * const a_data = data;
	ityp xa = a_data[0];
	ityp ya = a_data[1];
	ityp xb = a_data[2];
	ityp yb = a_data[3];
	ityp xc = a_data[4];
	ityp yc = a_data[5];
	
    bool value;
    ityp xmax;
    ityp xmin;
    ityp ymax;
    ityp ymin;

    if ( ! collinear ( xa, ya, xb, yb, xc, yc ) )
        value = false;
    else if ( fabs ( ya - yb ) < fabs ( xa - xb ) )
    {
        xmax = MAX ( xa, xb );
        xmin = MIN ( xa, xb );
        value = ( xmin <= xc && xc <= xmax );
    }
    else
    {
        ymax = MAX ( ya, yb );
        ymin = MIN ( ya, yb );
        value = ( ymin <= yc && yc <= ymax );
    }

	result = value;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _collinear ( void * data)
/******************************************************************************/
/*
  Purpose:
    COLLINEAR returns a measure of collinearity for three points.
  Discussion:
    In order to deal with collinear points whose coordinates are not
    numerically exact, we compare the area of the largest square
    that can be created by the line segment between two of the points
    to (twice) the area of the triangle formed by the points.
    If the points are collinear, their triangle has zero area.
    If the points are close to collinear, then the area of this triangle
    will be small relative to the square of the longest segment.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2014
  Author:
    Original C version by Joseph ORourke.
    This C version by John Burkardt.
  Reference:
    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
  Parameters:
    Input, double XA, YA, XB, YB, XC, YC, the coordinates of
    the vertices.
    Output, int COLLINEAR, is TRUE if the points are judged
    to be collinear.
*/
{
	static bool result = 2;
	
	ityp * const a_data = data;
	ityp xa = a_data[0];
	ityp ya = a_data[1];
	ityp xb = a_data[2];
	ityp yb = a_data[3];
	ityp xc = a_data[4];
	ityp yc = a_data[5];
	
    ityp area;
    const register ityp r8_eps = 2.220446049250313E-016;
    ityp side_max_sq;
    bool value;

    area = 0.50 * ( ( xb - xa ) * ( yc - ya ) - ( xc - xa ) * ( yb - ya ) );
    side_max_sq = MAX ( pow ( xa - xb, 2 ) + pow ( ya - yb, 2 ), MAX ( pow ( xb - xc, 2 ) + pow ( yb - yc, 2 ), pow ( xc - xa, 2 ) + pow ( yc - ya, 2 ) ) );
    
	result = side_max_sq <= r8_eps || 2.00 * fabs ( area ) <= r8_eps * side_max_sq;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _diagonal ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIAGONAL: VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2014
  Author:
    Original C version by Joseph ORourke.
    This C version by John Burkardt.
  Reference:
    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
  Parameters:
    Input, int IM1, IP1, the indices of two vertices.
    Input, int N, the number of vertices.
    Input, int PREV[N], the previous neighbor of each vertex.
    Input, int NEXT[N], the next neighbor of each vertex.
    Input, double X[N], Y[N], the coordinates of each vertex.
    Output, int DIAGONAL, the value of the test.
*/
{
	static bool result = 2;
	
	const _3dt2pi2pit * const s_data = data;
	const register dim_typ im1 = s_data->a0;
	const register dim_typ ip1 = s_data->a1;
	const register dim_typ n = s_data->a2;
	int * prev = s_data->a3;
	int * next = s_data->a4;
	ityp * x = s_data->a5;
	ityp * y = s_data->a6;
	
	result = ( in_cone ( im1, ip1, n, prev, next, x, y ) && in_cone ( ip1, im1, n, prev, next, x, y ) && diagonalie ( im1, ip1, n, next, x, y ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _diagonalie ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIAGONALIE is true if VERTEX(IM1):VERTEX(IP1) is a proper diagonal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2014
  Author:
    Original C version by Joseph ORourke.
    This C version by John Burkardt.
  Reference:
    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
  Parameters:
    Input, int IM1, IP1, the indices of two vertices.
    Input, int N, the number of vertices.
    Input, int NEXT[N], the next neighbor of each vertex.
    Input, double X[N], Y[N], the coordinates of each vertex.
    Output, int DIAGONALIE, the value of the test.
*/
{
	static bool result = 2;
	
	const dtpit2dtpipit * const s_data = data;
	
	const register dim_typ im1 = s_data->a0;
	ityp * x = s_data->a1;
	const register dim_typ ip1 = s_data->a2;
	const register dim_typ n = s_data->a3;
	int * next = s_data->a4;
	ityp * y = s_data->a5;
	
    int first;
    dim_typ j;
    int jp1;
    bool value;
    int value2;

    first = im1;
    j = first;
    jp1 = next[first];

    value = true;
    /*
    For each edge VERTEX(J):VERTEX(JP1) of the polygon:
    */
    while ( true )
    {
	    /*
	    Skip any edge that includes vertex IM1 or IP1.
	    */
	    if ( j == im1 || j == ip1 || jp1 == im1 || jp1 == ip1 );
	    else
	    {
	        value2 = intersect ( x[im1], y[im1], x[ip1], y[ip1], x[j], y[j],
	        x[jp1], y[jp1] );
	
	        if ( value2 )
	        {
	            value = false;
	            break;
	        }
	    }
        j = jp1;
        jp1 = next[j];

        if ( j == first )
            break;
    }

	result = value;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _in_cone ( void * data)
/******************************************************************************/
/*
  Purpose:
    IN_CONE is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2014
  Author:
    Original C version by Joseph ORourke.
    This C version by John Burkardt.
  Reference:
    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
  Parameters:
    Input, int IM1, IP1, the indices of two vertices.
    Input, int N, the number of vertices.
    Input, int PREV[N], the previous neighbor of each vertex.
    Input, int NEXT[N], the next neighbor of each vertex.
    Input, double X[N], Y[N], the coordinates of each vertex.
    Output, int IN_CONE, the value of the test.
*/
{
	static bool result = 2;
	
	const _3dt2pi2pit * const s_data = data;
	const register dim_typ im1 = s_data->a0;
	const register dim_typ ip1 = s_data->a1;
	const register dim_typ n = s_data->a2;
	int * prev = s_data->a3;
	int * next = s_data->a4;
	ityp * x = s_data->a5;
	ityp * y = s_data->a6;
	
    dim_typ i;
    int im2;
    ityp t1;
    ityp t2;
    ityp t3;
    ityp t4;
    ityp t5;
    bool value;

    im2 = prev[im1];
    i = next[im1];

    t1 = triangle_area_ptriang ( x[im1], y[im1], x[i], y[i], x[im2], y[im2] );

    if ( 0.00 <= t1 )
    {
        t2 = triangle_area_ptriang ( x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2] );
        t3 = triangle_area_ptriang ( x[ip1], y[ip1], x[im1], y[im1], x[i], y[i] );
        value = ( ( 0.00 < t2 ) && ( 0.00 < t3 ) );
    }
    else
    {
        t4 = triangle_area_ptriang ( x[im1], y[im1], x[ip1], y[ip1], x[i], y[i] );
        t5 = triangle_area_ptriang ( x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2] );
        value = ! ( ( 0.00 <= t4 ) && ( 0.00 <= t5 ) );
    }
    
    result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _intersect ( void * data)
/******************************************************************************/
/*
  Purpose:
    INTERSECT is true if lines VA:VB and VC:VD intersect.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2014
  Author:
    Original C version by Joseph ORourke.
    This C version by John Burkardt.
  Reference:
    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
  Parameters:
    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y
    coordinates of the four vertices.
    Output, int INTERSECT, the value of the test.
*/
{
	static bool result = 2;
	
	ityp * const a_data = data;
	ityp xa = a_data[0];
	ityp ya = a_data[1];
	ityp xb = a_data[2];
	ityp yb = a_data[3];
	ityp xc = a_data[4];
	ityp yc = a_data[5];
	ityp xd = a_data[6];
	ityp yd = a_data[7];
	
	result = intersect_prop ( xa, ya, xb, yb, xc, yc, xc, yd ) ||
         between ( xa, ya, xb, yb, xc, yc ) ||
         between ( xa, ya, xb, yb, xd, yd ) ||
         between ( xc, yc, xd, yd, xa, ya ) ||
         between ( xc, yc, xd, yd, xb, yb );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _intersect_prop ( void * data)
/******************************************************************************/
/*
  Purpose:
    INTERSECT_PROP is TRUE if lines VA:VB and VC:VD have a proper intersection.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 May 2014
  Author:
    Original C version by Joseph ORourke.
    This C version by John Burkardt.
  Reference:
    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
  Parameters:
    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y
    coordinates of the four vertices.
    Output, int INTERSECT_PROP, the result of the test.
*/
{
	static bool result = 2;
	
	ityp * const a_data = data;
	ityp xa = a_data[0];
	ityp ya = a_data[1];
	ityp xb = a_data[2];
	ityp yb = a_data[3];
	ityp xc = a_data[4];
	ityp yc = a_data[5];
	ityp xd = a_data[6];
	ityp yd = a_data[7];
	
    ityp t1;
    ityp t2;
    ityp t3;
    ityp t4;
    bool value;
    bool value1;
    bool value2;
    bool value3;
    bool value4;

    if ( collinear ( xa, ya, xb, yb, xc, yc )  ||
    collinear ( xa, ya, xb, yb, xd, yd ) ||
    collinear ( xa, ya, xb, yb, xd, yd ) ||
    collinear ( xc, yc, xd, yd, xa, ya ) ||
    collinear ( xc, yc, xd, yd, xb, yb ))
    {
    	result = true;
        return &result;
    }
    else
    {
        t1 = triangle_area_ptriang ( xa, ya, xb, yb, xc, yc );
        t2 = triangle_area_ptriang ( xa, ya, xb, yb, xd, yd );
        t3 = triangle_area_ptriang ( xc, yc, xd, yd, xa, ya );
        t4 = triangle_area_ptriang ( xc, yc, xd, yd, xb, yb );

        value1 = ( 0.00 < t1 );
        value2 = ( 0.00 < t2 );
        value3 = ( 0.00 < t3 );
        value4 = ( 0.00 < t4 );

        value = ( l4_xor ( value1, value2 ) ) && ( l4_xor ( value3, value4 ) );
    }

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_area_ptriang ( void * data)
/******************************************************************************/
/*
  Purpose:
    triangle_area_ptriang computes the signed area of a triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, double XA, YA, XB, YB, XC, YC, the coordinates of
    the vertices of the triangle, given in counterclockwise order.
    Output, double TRIANGLE_AREA_PTRIANG, the signed area of the triangle.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	ityp xa = a_data[0];
	ityp ya = a_data[1];
	ityp xb = a_data[2];
	ityp yb = a_data[3];
	ityp xc = a_data[4];
	ityp yc = a_data[5];
	
	result = 0.50 * ( ( xb - xa ) * ( yc - ya ) - ( xc - xa ) * ( yb - ya ) );
	return &result;
}

#endif
