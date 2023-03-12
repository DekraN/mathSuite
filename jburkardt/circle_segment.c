#ifndef __DISABLEDEEP_CIRCLESEGMENT

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_angle_from_chord ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_ANGLE_FROM_CHORD computes the angle of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the center of the circle.
    Input, double P1[2], P2[2], the ends of the chord.
    Output, double CIRCLE_SEGMENT_ANGLE_FROM_CHORD, the value of THETA,
    the angle of the circle segment.  0 <= THETA < 2 * M_PI.
*/
{
	static ityp result = MAX_VAL;
	
	const it3pit * const s_data = data;
	ityp r = s_data->a0;
	ityp * c = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;
	
    ityp theta;
    ityp v1[2];
    ityp v2[2];
    /*
    Compute the radial vectors V1 and V2.
    */
    v1[0] = p1[0] - c[0];
    v1[1] = p1[1] - c[1];
    v2[0] = p2[0] - c[0];
    v2[1] = p2[1] - c[1];
    /*
    The arc cosine will only give us an answer between 0 and M_PI.
    */
    theta = atan2 ( v2[1], v2[0] ) - atan2 ( v1[1], v1[0] );
    /*
    Force 0 <= THETA < 2 * M_PI.
    */
    while ( theta < 0.00 )
        theta  +=  M_2TPI;

    while ( M_2TPI <= theta )
        theta -= M_2TPI;
        
    result = theta;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_angle_from_chord_angles ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_ANGLE_FROM_CHORD_ANGLES computes angle of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double OMEGA1, OMEGA2, the angles of the points P1
    and P2.  OMEGA1 <= OMEGA2.
    Output, double CIRCLE_SEGMENT_ANGLE_FROM_CHORD_ANGLES, the angle THETA
    of the circle segment.  Essentially, THETA = OMEGA2 - OMEGA1.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp omega1 = a_data[0];
	register ityp omega2 = a_data[1];
	
    while ( omega2 < omega1 )
        omega2 += M_2TPI;
    
    result = omega2-omega1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _circle_segment_angle_from_height ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT computes the angle of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double H, the "height" of the circle segment.
    0 <= H <= 2 * R.
    Output, double CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT, the angle THETA
    of the circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp r = a_data[0];
	const register ityp h = a_data[1];
	
    ityp theta;

    if ( h <= 0.00 )
        theta = 0.00;
    else if ( h <= r )
    {
        theta = 2.00 * acos ( ( r - h ) / r );
        theta = 2.00 * asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
    }
    else if ( h <= 2.00 * r )
    {
        theta = 2.00 * acos ( ( r - h ) / r );
        theta = 2.00 * asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
        theta = M_2TPI - theta;
    }
    else
        theta = M_2TPI;
        
    result = theta;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_area_from_angle ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_AREA_FROM_ANGLE computes the area of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burr * r * ( theta - sin ( theta ) ) / 2.0kardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double THETA, the angle of the circle segment.
    Output, double CIRCLE_SEGMENT_AREA_FROM_ANGLE, the area of the
    circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp r = a_data[0];
	const register ityp theta = a_data[1];
	
	result = r * r * ( theta - sin ( theta ) ) / 2.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_area_from_chord ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_AREA_FROM_CHORD computes the area of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the center of the circle.
    Input, double P1[2], P2[2], the ends of the chord.
    Output, double CIRCLE_SEGMENT_AREA_FROM_CHORD, the area of the
    circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	const it3pit * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * c = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;
	
    const register ityp theta = circle_segment_angle_from_chord ( r, c, p1, p2 );
    
    result = r * r * ( theta - sin ( theta ) ) / 2.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_area_from_height ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double H, the height of the circle segment.
    0 <= H <= 2 * R.
    Output, double CIRCLE_SEGMENT_AREA_FROM_HEIGHT, the area of the
    circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp r = a_data[0];
	const register ityp h = a_data[1];
	
    ityp area, theta;

    if ( h <= 0.00 )
        area = 0.00;
    else if ( h <= r )
    {
        theta = 2.00 * asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
        area = r * r * ( theta - sin ( theta ) ) / 2.00;
    }
    else if ( h <= 2.00 * r )
    {
        theta = 2.00 * asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r );
        theta = M_2TPI - theta;
        area = r * r * ( theta - sin ( theta ) ) / 2.00;
    }
    else
        area = M_PI * r * r;

	result = area;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_area_from_sample ( void * data )
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_AREA_FROM_SAMPLE computes the area of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.

  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the center of the circle.
    Input, double P1[2], P2[2], the ends of the chord.
    Input, int N, the number of sample points.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double CIRCLE_SEGMENT_AREA_FROM_SAMPLE, the area of the
    circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	const it3pitdtpi * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * c = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;
	const register dim_typ n = s_data->a4;
	int * seed = s_data->a5;
	
    ityp *angle;
    ityp area;
    int i;
    int m;
    ityp omega1;
    ityp omega2;
    ityp p[2];
    ityp *r2;
    ityp rmh;
    ityp vdotp[n];
    ityp x[n];
    ityp y[n];
    /*
    Determine the angles of the chord endpoints.
    */
    omega1 = atan2 ( p1[1] - c[1], p1[0] - c[0] );
    while ( omega1 < 0.00 )
        omega1 += M_2TPI;

    omega2 = atan2 ( p2[1] - c[1], p2[0] - c[0] );
    while ( omega2 < omega1 )
        omega2 +=  M_2TPI;
    /*
    Get N random points in the circle.
    To simplify angle measurement, take OMEGA1 as your smallest angle.
    That way, the check OMEGA1 <= ANGLE <= OMEGA2 will be legitimate.
    */
    angle = r8vec_uniform_01_new ( n, seed );
    for ( i = 0; i < n; ++i )
        angle[i] = omega1 + M_2TPI * angle[i];
    r2 = r8vec_uniform_01_new ( n, seed );
    for ( i = 0; i < n; ++i)
        r2[i] = sqrt ( r2[i] );
    for ( i = 0; i < n; ++i )
    {
        x[i] = c[0] + r2[i] * cos ( angle[i] );
        y[i] = c[0] + r2[i] * sin ( angle[i] );
    }
    /*
    Determine the vector that touches the circle segment base.
    */
    p[0] = 0.50 * ( p1[0] + p2[0] ) - c[0];
    p[1] = 0.50 * ( p1[1] + p2[1] ) - c[1];

    rmh = sqrt ( p[0] * p[0] + p[1] * p[1] );

    p[0] /= rmh;
    p[1] /= rmh;

    if ( M_PI < omega2 - omega1 )
    {
        p[0] = - p[0];
        p[1] = - p[1];
        rmh *= -1;
    }
    /*
    Compute the projection of each point onto P.
    */
    for ( i = 0; i < n; ++i )
        vdotp[i] = ( x[i] - c[0] ) * p[0] + ( y[i] - c[1] ) * p[1];
    /*
    Points in the segment lie in the sector, and project at least RMH onto P.
    */
    m = 0;
    for ( i = 0; i < n; ++i)
        if ( omega1 < angle[i] && angle[i] < omega2 && rmh < vdotp[i] )
            ++ m;
    /*
    The area of the segment is its relative share of the circle area.
    */
    free ( angle );
    free ( r2 );
    
    result = M_PI * r * r * ( ityp ) ( m ) / ( ityp ) ( n );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT void *   _circle_segment_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_CDF computes a CDF related to a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
    Now, suppose we want to assign a cumulative density function or CDF
    based on a variable H2 which measures the height of the circle segment
    formed by an arbitrary point in the circle segment.  CDF(H2) will
    measure the probability that a point drawn uniformly at random from
    the circle segment defines a (smaller) circle segment of height H2.
    If we can define this CDF, then we will be able to sample uniformly
    from the circle segment, since our "Y" value can be determined from H2,
    and our X value is chosen uniformly over the horizontal chord
    through (0,Y).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double H, the "height" of the circle segment.
    0 <= H <= 2 * R.
    Input, double H2, the "height" of the new circle segment
    defined by a given point in the circle segment.  0 <= H2 <= H.
    Output, double CDF, the cumulative density function for H2,
    the probability that a point chosen at random in the circle segment
    would define a smaller circle segment of height H2 or less.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp r = a_data[0];
	const register ityp h = a_data[1];
	const register ityp h2 = a_data[2];
	
    ityp a, a2, cdf;

    if ( h2 <= 0.00 )
        cdf = 0.00;
    else if ( h <= h2 )
        cdf = 1.00;
    else
    {
        a = circle_segment_area_from_height ( r, h  );
        a2 = circle_segment_area_from_height ( r, h2 );
        cdf = a2 / a;
    }
    
    result = cdf;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _circle_segment_centroid_from_chord ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_CENTROID_FROM_CHORD computes the centroid of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
    For this function, we assume that the center of the circle is at (0,0),
    that the chord is horizontal, and that the circle segment is at the top.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the center of the circle.
    Input, double P1[2], P2[2], the coordinates of the endpoints
    of the chord.
    Output, double CIRCLE_SEGMENT_CENTROID_FROM_CHORD[2], the coordinates
    of the centroid.
*/
{
	const it3pit * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * c = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;
	
	
    ityp *d;
    ityp s;
    ityp theta;
    ityp thetah;
    ityp v1[2];
    /*
    Get the angle subtended by P1:P2.
    */
    theta = circle_segment_angle_from_chord ( r, c, p1, p2 );
    /*
    Construct V1, the vector from C to P1.
    */
    v1[0] = p1[0] - c[0];
    v1[1] = p1[1] - c[1];
    /*
    Rotate V1 through THETA / 2.
    */
    thetah = theta / 2.00;

    d = ( ityp * ) malloc ( sizeof ( ityp ) << 1);
    d[0] = cos ( thetah ) * v1[0] - sin ( thetah ) * v1[1];
    d[1] = sin ( thetah ) * v1[0] + cos ( thetah ) * v1[1];
    /*
    Scale this vector so it represents the distance to the centroid
    relative to R.
    */
    s = 4.00 * pow ( sin ( theta / 2.00 ), 3 ) / 3.00 / ( theta - sin ( theta ) );

    d[0] = s * d[0];
    d[1] = s * d[1];
    /*
    Add the center.
    */
    d[0] += c[0];
    d[1] += c[1];

    return d;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _circle_segment_centroid_from_height ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_CENTROID_FROM_HEIGHT computes centroid of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
    For this function, we assume that the center of the circle is at (0,0).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double H, the "height" of the circle segment.
    0 <= H <= 2 * R.
    Output, double CIRCLE_SEGMENT_CENTROID_FROM_HEIGHT[2], the coordinates
    of the centroid.
*/
{
	ityp * const a_data = data;
	const register ityp r = a_data[0];
	const register ityp h = a_data[1];
	
    ityp *d;
    const register ityp theta = circle_segment_angle_from_height ( r, h );
    d = ( ityp * ) malloc ( sizeof ( ityp ) << 1 );
    d[0] = 0.00;
    d[1] = 4.00 * r * pow ( sin ( theta / 2.0 ), 3 ) / 3.0 / ( theta - sin ( theta ) );
    return d;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _circle_segment_centroid_from_sample ( void * data )
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE estimates a circle segment centroid.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the center of the circle.
    Input, double P1[2], P2[2], the ends of the chord.
    Input, int N, the number of sample points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE[2], the estimated
    centroid of the circle segment.
*/
{
	const it3pitdtpi * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * c = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;
	const register dim_typ n = s_data->a4;
	int * seed = s_data->a5;
	
    ityp *d;
    ityp x[n];
    ityp y[n];
    circle_segment_sample_from_chord ( r, c, p1, p2, n, seed, x, y );

    d = ( ityp * ) malloc ( sizeof ( ityp ) << 1);

    d[0] = r8vec_sum ( n, x ) / ( ityp ) ( n );
    d[1] = r8vec_sum ( n, y ) / ( ityp ) ( n );
    free ( x );
    free ( y );
    return d;
}

#define CIRCLESEGMENTCONTAINSPOINT_NEGERROR false

/******************************************************************************/
__MATHSUITE __JBURKARDT  void  * _circle_segment_contains_point ( void * data) 
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_CONTAINS_POINT reports whether a point is in a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
    In this function, we allow the circle to have an arbitrary center C,
    arbitrary radius R, and we describe the points P1 and P2 by specifying
    their angles OMEGA1 and OMEGA2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the center of the circle.
    Input, double OMEGA1, OMEGA2, the angles of the two points on
    the circumference of the circle that define the circle segment.
    OMEGA1 < OMEGA2 <= OMEGA1 + 2 * M_PI
    Input, double XY[2], a point.
    Output, int CIRCLE_SEGMENT_CONTAINS_POINT, is TRUE if the point is inside
    the circle segment.
*/
{
	static bool result = 2;
	
	const itpit2itpit * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * c = s_data->a1;
	const register ityp omega1 = s_data->a2;
	register ityp omega2 = s_data->a3;
	ityp * xy = s_data->a4;
	
    ityp h;
    ityp omegah;
    ityp theta;
    ityp v[2];
    ityp v_omega;
    ityp v_project;
    ityp v_r;
    bool value;

    if ( r <= 0.00 )
    {
    	result = CIRCLESEGMENTCONTAINSPOINT_NEGERROR;
        return &result;
    }

    while ( omega2 < omega1 )
    omega2 += M_2TPI;
    /*
    Compute the vector V = XY - C:
    */
    v[0] = xy[0] - c[0];
    v[1] = xy[1] - c[1];
    /*
    a: Point must be inside the circle, so ||V|| <= R.
    */
    v_r = sqrt ( v[0] * v[0] + v[1] * v[1] );

    if ( r < v_r )
	{
		result = false;
        return &result;
    }
        
    /*
    b: Angle made by the vector V must be between OMEGA1 and OMEGA2.
    */
    v_omega = atan2 ( v[1], v[0] );

    while ( omega1 <= v_omega + M_2TPI )
        v_omega -= M_2TPI;

    while ( v_omega + M_2TPI <= omega1 )
        v_omega += M_2TPI;

    if ( omega2 < v_omega )
	{
		result = false;
        return &result;
	}
    /*
    c: Projection of V onto unit centerline must be at least R-H.
    */
    omegah = 0.50 * ( omega1 + omega2 );
    v_project = v[0] * cos ( omegah ) + v[1] * sin ( omegah );

    theta = omega2 - omega1;
    h = circle_segment_height_from_angle ( r, theta );

    if ( v_project < r - h )
    {
    	result = false;
        return &result;
    }

	result = true;
    return &result;
}

#define CIRLCESEGMENTHEIGHTFROMANGLE_ERROR MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_height_from_angle ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE: height of a circle segment from angle.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
    This function is given the radius R and angle of the segment, and
    determines the corresponding height.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double ANGLE, the angle of the circle segment.
    0 <= ANGLE <= M_2TPI.
    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE, the height of the
    circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp r = a_data[0];
	const register ityp angle = a_data[1];
	
    if ( angle < 0.00 || angle > M_2TPI || r <= 0 )
    {
    	result = CIRLCESEGMENTHEIGHTFROMANGLE_ERROR;
        return &result;
    }
    if(angle == 0.00)
    {
    	result = 0.00;
        return &result;
    }
    if(angle == M_2TPI)
    {
    	result = 2.00*r;
        return &result;
    }
    
    result = r * angle <= M_PI ? ( 1.00 - cos (              angle   / 2.00 ) ) : ( 1.00 + cos ( ( M_2TPI - angle ) / 2.00 ) );
    return &result;
}

#define CIRCLESEGMENTHEIGHTFROMAREA_ERROR MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_height_from_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_HEIGHT_FROM_AREA: height of a circle segment from area.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
    This function is given the radius R and area of the segment, and
    determines the corresponding height.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double AREA, the area of the circle segment.
    0 <= AREA <= M_2TPI * R^2.
    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_AREA, the height of the
    circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp r = a_data[0];
	const register ityp area = a_data[1];
	
    ityp a;
    ityp a1;
    ityp a2;
    ityp area_circle;
    ityp eps;
    ityp h;
    ityp h1;
    ityp h2;
    dim_typ it;

    if(area < 0.00 || area > area_circle || r <= 0.00)
    {
    	result = CIRCLESEGMENTHEIGHTFROMAREA_ERROR;
        return &result;
    }


    area_circle = M_2TPI * r * r;

    if ( area == 0.0 )
    {
    	result = 0.00;
        return &result;
    }

    if ( area == area_circle )
    {
    	result = 2.00*r;
        return &result;
    }

    h1 = 0.00;
    a1 = circle_segment_area_from_height ( r, h1 );
    h2 = 20.0 * r;
    a2 = circle_segment_area_from_height ( r, h2 );

    it = 0;
    eps = r8_epsilon ( );

    while ( it < 30 )
    {
        h = 0.50 * ( h1 + h2 );
        a = circle_segment_area_from_height ( r, h );
        it = it + 1;

        if ( fabs ( a - area ) < sqrt ( eps ) * area )
            break;

        if ( a < area )
        {
            h1 = h;
            a1 = a;
        }
        else
        {
            h2 = h;
            a2 = a;
        }
    }
    
    result = h;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_height_from_chord ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_HEIGHT_FROM_CHORD: height of a circle segment from chord.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the coordinates of the circle center.
    Input, double P1[2], P2[2], the coordinates of the
    chord endpoints.
    Output, double CIRCLE_SEGMENT_HEIGHT_FROM_CHORD, the height of the circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	const it3pit * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * c = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;
	
	result = circle_segment_height_from_angle ( r, circle_segment_angle_from_chord ( r, c, p1, p2 ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_rotation_from_chord ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_ROTATION_FROM_CHORD computes the rotation of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the center of the circle.
    Input, double P1[2], P2[2], the ends of the chord.
    Warning! If P1 = P2, we can't tell whether the segment is the whole
    circle or none of it!
    Output, double CIRCLE_SEGMENT_ROTATION_FROM_CHORD, the rotation of the
    circle segment.  0 <= ALPHA < 2 * M_PI.
*/
{
	static ityp result = MAX_VAL;

	const it3pit * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * c = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;	
	
    ityp alpha;
    ityp rho1;
    ityp rho2;
    ityp theta;
    ityp v1[2];
    ityp v2[2];
    /*
    Compute the radial vectors V1 and V2.
    */
    v1[0] = p1[0] - c[0];
    v1[1] = p1[1] - c[1];
    v2[0] = p2[0] - c[0];
    v2[1] = p2[1] - c[1];
    /*
    Use R8_ATAN to guarantee that 0 <= RHO1, RHO2 <= 2 * M_PI.
    */
    rho1 = atan2 ( v1[1], v1[0] );
    rho2 = atan2 ( v2[1], v2[0] );
    /*
    Force RHO2 to be bigger than RHO1.
    */
    while ( rho2 <= rho1 )
        rho2 += M_2TPI;
    /*
    Compute THETA.
    */
    theta = rho2 - rho1;
    /*
    ALPHA is RHO1, plus half of the angular distance between P1 and P2.
    */
    alpha = rho1 + 0.50 * theta;

    while ( M_2TPI <= alpha )
        alpha -= M_2TPI;

	result = alpha;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _circle_segment_sample_from_chord ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_SAMPLE_FROM_CHORD samples points from a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double C[2], the center of the circle.
    Input, double P1[2], P2[2], the endpoints of the chord.
    Input, int N, the number of sample points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double X[N], Y[N], the sample points.
*/
{
	const it3pitdtpi2pit * const s_data = data;
	const register ityp r = s_data->a0;
	ityp * c = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;
	const register dim_typ n = s_data->a4;
	int * seed = s_data->a5;
	ityp * x = s_data->a6;
	ityp * y = s_data->a7;
	
    ityp c2[2];
    ityp h;
    dim_typ i;
    ityp t;
    ityp vc[2];
    ityp vr[2];
    ityp xi[n];
    ityp eta[n];
    /*
    Determine unit vectors VR and VC.
    VR points to the center of the chord from the radius.
    VC points along the chord, from P1 to P2.
    */
    vr[0] = 0.50 * ( p1[0] + p2[0] ) - c[0];
    vr[1] = 0.50 * ( p1[1] + p2[1] ) - c[1];

    t = sqrt ( vr[0] * vr[0] + vr[1] * vr[1] );
    vr[0] /=  t;
    vr[1] /= t;

    vc[0] = p2[0] - p1[0];
    vc[1] = p2[1] - p1[1];

    t = sqrt ( vc[0] * vc[0] + vc[1] * vc[1] );
    vc[0] /=  t;
    vc[1] /= t;
    /*
    Get the height of the circle segment.
    */
    c2[0] = 0.00;
    c2[1] = 0.00;
    h = circle_segment_height_from_chord ( r, c2, p1, p2 );
    /*
    Sample (xi,eta) in the reference coordinates, where the chord
    is horizontal.
    */
    circle_segment_sample_from_height ( r, h, n, seed, xi, eta );
    /*
    XI is the left/right coordinate along VC.
    ETA is the distance along VR.
    */
    for ( i = 0; i < n; ++i )
    {
        x[i] = c[0] + eta[i] * vr[0] + xi[i] * vc[0];
        y[i] = c[1] + eta[i] * vr[1] + xi[i] * vc[1];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _circle_segment_sample_from_height ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT samples points from a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double H, the height of the circle segment.
    0 <= H <= 2 * R.
    Input, int N, the number of sample points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double X[N], Y[N], the sample points.
*/
{
	const _2itdtpi2pit * const s_data = data;
	const register ityp r = s_data->a0;
	const register ityp h = s_data->a1;
	const register dim_typ n = s_data->a2;
	int * seed = s_data->a3;
	ityp * x = s_data->a4;
	ityp * y = s_data->a5;
	
    ityp area;
    dim_typ i;
    ityp *u;
    ityp wh[n];
    ityp area2[n];
    ityp h2[n];

    area = circle_segment_area_from_height ( r, h );
    /*
    Pick CDF's randomly.
    */
    u = r8vec_uniform_01_new ( n, seed );
    /*
    Choose points randomly by choosing ordered areas randomly.
    */
    for ( i = 0; i < n; ++i )
        area2[i] = u[i] * area;
    /*
    Each area corresponds to a height H2.  Find it.
    */
    for ( i = 0; i < n; ++i )
        h2[i] = circle_segment_height_from_area ( r, area2[i] );
    /*
    Determine the half-width WH of the segment for each H2.
    */
    for ( i = 0; i < n; ++i )
        wh[i] = sqrt ( h2[i] * ( 2.0 * r - h2[i] ) );
    /*
    Choose an X position randomly in [-WH,+WH].
    */
    free ( u );
    u = r8vec_uniform_01_new ( n, seed );

    for ( i = 0; i < n; ++i )
        x[i] = ( 2.00 * u[i] - 1.00 ) * wh[i];
    /*
    Our circle center is at (0,0).  Our height of H2 is subtracted
    from the height R at the peak of the circle.  Determine the Y
    coordinate using this fact.
    */
    for ( i = 0; i < n; ++i )
        y[i] = r - h2[i];

    free ( u );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _circle_segment_width_from_height ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_SEGMENT_WIDTH_FROM_HEIGHT computes the width of a circle segment.
  Discussion:
    Begin with a circle of radius R.  Choose two points P1 and P2 on the
    circle, and draw the chord P1:P2.  This chord divides the circle
    into two pieces, each of which is called a circle segment.
    Consider one of the pieces.  The "angle" of this segment is the angle
    P1:C:P2, where C is the center of the circle.  Let Q be the point on
    the chord P1:P2 which is closest to C.  The "height" of the segment
    is the distance from Q to the perimeter of the circle.  The "width"
    of the circle segment is the length of P1:P2.
    This function is given the radius R and height H of the segment, and
    determines the corresponding width W.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double R, the radius of the circle.
    0 < R.
    Input, double H, the height of the circle segment.
    0 <= H <= 2 * R.
    Output, double CIRCLE_SEGMENT_WIDTH_FROM_HEIGHT, the width of the
    circle segment.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp r = a_data[0];
	const register ityp h = a_data[1];
	
	result = 2.00 * sqrt ( h * ( 2.00 * r - h ) );
    return &result;
}

#endif
