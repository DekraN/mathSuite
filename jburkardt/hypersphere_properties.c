#ifndef __DISABLEDEEP_HYPERSPHEREPROPERTIES

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cartesian_to_hypersphere ( void * data)
/******************************************************************************/
/*
  Purpose:
    CARTESIAN_TO_HYPERSPHERE: Cartesian to hypersphere coordinate transform.
  Discussion:
    We allow the trivial case M = 1; in that case alone, the value R
    must be assumed to have a sign.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 May 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    1 <= M.
    Input, int N, the number of points to transform.
    Input, double C[M], the center of the hypersphere.
    Input, double X[M*N], the Cartesian coordinates of the points.
    Output, double R[N], the radius of the points on the
    hypersphere.  Except for the trivial case M = 1, R is assumed nonnegative.
    Output, double THETA[(M-1)*N], the coordinate angles of the
    points, measured in radians.
*/
{
	const _2dt4pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * c = s_data->a2;
	ityp * x = s_data->a3;
	ityp * r = s_data->a4;
	ityp * theta = s_data->a5;
	
    dim_typ i, i1, j;
    ityp t;
    ityp *top;
    ityp *x2;
    /*
    Handle special case of M = 1.
    */
    if ( m == 1 )
    {
        for ( j = 0; j < n; ++j )
            r[j] = x[0+j*m] - c[0];
        return NULL;
    }
    /*
    Subtract the center.
    */
    x2 = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i)
            x2[i+j*m] = x[i+j*m] - c[i];
    /*
    Compute R.
    */
    for ( j = 0; j < n; ++j )
    {
        t = 0.00;
        for ( i = 0; i < m; ++i )
            t +=  pow ( x2[i+m*j], 2 );
        r[j] = sqrt ( t );
    }
    /*
    Compute M-2 components of THETA.
    */
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m - 1; ++i )
            theta[i+j*(m-1)] = 0.00;

    for ( j = 0; j < n; ++j )
        for ( i = 1; i < m - 1; ++i)
            for ( i1 = 0; i1 <= i - 1; ++i1 )
                theta[i1+j*(m-1)] += pow ( x2[i+j*m], 2 );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m - 2; ++i )
            theta[i+j*(m-1)] += pow ( x2[m-1+j*m], 2 );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m - 2; ++i )
            theta[i+j*(m-1)] = atan2 ( sqrt ( theta[i+j*(m-1)] ), x2[i+j*m] );
    /*
    Compute last component of THETA.
    */
    top = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
        top[j] = sqrt ( pow ( x2[m-1+j*m], 2 ) + pow ( x2[m-2+j*m], 2 ) ) + x2[m-2+j*m];

    for ( j = 0; j < n; ++j )
        theta[m-2+j*(m-1)] = 2.00 * atan2 ( x2[m-1+j*m], top[j] );

    free ( top );
    free ( x2 );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hypersphere_01_interior_uniform ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_01_INTERIOR_UNIFORM: uniform points inside unit hypersphere.
  Discussion:
    The hypersphere has center 0 and radius 1.
    This routine is valid for any spatial dimension.
    We first generate a point ON the hypersphere, and then distribute it
    IN the hypersphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 December 2013
  Author:
    John Burkardt
  Reference:
    Russell Cheng,
    Random Variate Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998, pages 168.
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.
  Parameters:
    Input, int M, the dimension of the space.
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double HYPERSPHERE_01_INTERIOR_UNIFORM[M*N], the points.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
    ityp exponent;
    dim_typ i, j;
    ityp norm;
    ityp r;
    ityp *v;
    ityp *x;

    x = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    exponent = 1.00 / ( double ) ( m );

    for ( j = 0; j < n; ++j )
    {
    /*
        Fill a vector with normally distributed values.
        */
        v = r8vec_normal_01_new ( m, seed );
        /*
        Compute the length of the vector.
        */
        norm = 0.00;
        for ( i = 0; i < m; ++i )
            norm += pow ( v[i], 2 );
        norm = sqrt ( norm );
        /*
        Normalize the vector.
        */
        for ( i = 0; i < m; ++i)
            x[i+j*m] = v[i] / norm;
        /*
        Now compute a value to map the point ON the hypersphere INTO the hypersphere.
        */
        r = r8_uniform_01 ( seed );

        for ( i = 0; i < m; ++i )
            x[i+j*m] = pow ( r, exponent ) * x[i+j*m];
        free ( v );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hypersphere_01_surface_uniform ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_01_SURFACE_UNIFORM: uniform points on unit hypersphere surface.
  Discussion:
    The hypersphere has center 0 and radius 1.
    This procedure is valid for any spatial dimension.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 November 2010
  Author:
    John Burkardt
  Reference:
    Russell Cheng,
    Random Variate Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998, pages 168.
    George Marsaglia,
    Choosing a point from the surface of a sphere,
    Annals of Mathematical Statistics,
    Volume 43, Number 2, April 1972, pages 645-646.
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.
  Parameters:
    Input, int M, the dimension of the space.
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double HYPERSPHERE_01_UNIFORM_SURFACE[M*N], the points.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
    dim_typ i, j;
    dim_typ norm;
    ityp *x;
    /*
    Fill a matrix with normally distributed values.
    */
    x = r8mat_normal_01_new ( m, n, seed );
    /*
    Normalize each column.
    */
    for ( j = 0; j < n; ++j )
    {
        /*
        Compute the length of the vector.
        */
        norm = 0.00;
        for ( i = 0; i < m; ++i)
            norm += pow ( x[i+j*m], 2 );
            norm = sqrt ( norm );
        /*
        Normalize the vector.
        */
        for ( i = 0; i < m; ++i )
            x[i+j*m] /= norm;
    }
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hypersphere_01_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_01_VOLUME computes the volume of a unit hypersphere.
  Discussion:
    The unit hypersphere satisfies the equation:
      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
     M  Volume
     1    2
     2    1        * M_PI
     3 ( 4 /   3) * M_PI
     4 ( 1 /   2) * M_PI^2
     5 ( 8 /  15) * M_PI^2
     6 ( 1 /   6) * M_PI^3
     7 (16 / 105) * M_PI^3
     8 ( 1 /  24) * M_PI^4
     9 (32 / 945) * M_PI^4
    10 ( 1 / 120) * M_PI^5
    For the unit hypersphere, Volume(M) = 2 * M_PI * Volume(M-2)/ M
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 September 2003
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of the space.
    Output, double HYPERSPHERE_01_VOLUME, the volume.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ m = *(dim_typ *) data;
	
    dim_typ i, m2;
    ityp volume;

    if ( m % 2== 0 )
    {
        m2 = m / 2;
        volume = 1.00;
        for ( i = 1; i <= m2; ++i )
            volume *= M_PI / ( ( ityp ) i );
    }
    else
    {
        m2 = ( m - 1 ) / 2;
        volume = pow ( M_PI, m2 ) * pow ( 2.00, m );
        for ( i = m2 + 1; i <= (m2<<1) + 1; i++ )
            volume /= ( ( ityp ) i );
    }
    
    result = volume;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hypersphere_01_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_01_AREA computes the surface area of a unit hypersphere.
  Discussion:
    The unit hypersphere satisfies the equation:
      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
    M   Area
     2    2        * PI
     3    4        * PI
     4 ( 2 /   1) * PI^2
     5 ( 8 /   3) * PI^2
     6 ( 1 /   1) * PI^3
     7 (16 /  15) * PI^3
     8 ( 1 /   3) * PI^4
     9 (32 / 105) * PI^4
    10 ( 1 /  12) * PI^5
    For the unit hypersphere, Area(M) = M * Volume(M)
    Sphere_Unit_Area ( M ) = 2 * PI^(M/2) / Gamma ( M / 2 )
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    05 September 2003
  Author:
    John Burkardt
  Parame
    Input, int M, the dimension of the space.
    Output, double HYPERSPHERE_01_AREA, the area.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ m = *(dim_typ *) data;
	
	ityp area;
	int i;
	int m2;
	
	if ( ( m % 2 ) == 0 )
	{
		m2 = m / 2;
		area = 2.00 * pow ( M_PI, m2 );
		for ( i = 1; i <= m2 - 1; ++i )
			area /= ( ( ityp ) i );
	}
	else
	{
		m2 = ( m - 1 ) / 2;
		area = pow ( 2.00, m ) * pow ( M_PI, m2 );
		for ( i = m2 + 1; i <= (m2<<1); ++i )
			area /= ( ( ityp ) i );
	}
	
	result = area;
	return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hypersphere_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_AREA computes the surface area of hypersphere.
  Discussion:
    A hypersphere satisfies the equation:
      sum ( ( P(1:M) - C(1:M) )^2 ) = R^2
    M   Area
    2      2       * M_PI   * R
    3      4       * M_PI   * R^2
    4      2       * M_PI^2 * R^3
    5   (8/3)   * M_PI^2 * R^4
    6                M_PI^3 * R^5
    7   (16/15) * M_PI^3 * R^6
    Sphere_Area ( M, R ) = 2 * M_PI^(M/2) * R^(M-1) / Gamma ( M / 2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2003
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of the space.
    Input, double R, the radius.
    Output, double HYPERSPHERE_AREA, the area.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register ityp r = s_data->a1;
	
	result = pow ( r, m - 1  ) * hypersphere_01_area ( m );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hypersphere_stereograph ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_STEREOGRAPH: stereographic mapping of points on a hypersphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 May 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    M must be at least 2.
    Input, int N, the number of points.
    Input, double X[M*N], the points to be mapped.
    Output, double HYPERSPHERE_STEREOGRAPH[M-1)*N], the stereographically
    mapped points.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	
    dim_typ i, j;
    ityp *x2;

    x2 = ( ityp * ) malloc ( ( m - 1 ) * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m - 1; ++i )
            x2[i+j*(m-1)] = x[i+j*m];

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m - 1; ++i )
            x2[i+j*(m-1)] /= ( 1.00 - x[m-1+j*m] );

    return x2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hypersphere_stereograph_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_STEREOGRAPH_INVERSE inverts a stereographic map.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 May 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    M must be at least 2.
    Input, int N, the number of points.
    Input, double X2[(M-1)*N], points in the plane.
    Input, double HYPERSPHERE_STEREOGRAPH_INVERSE[M*N], points mapped
    back to the hypersphere.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x2 = s_data->a2;
	
    ityp *d;
    dim_typ i, j;
    ityp *x;

    x = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m - 1; ++i )
            x[i+j*m] = 2.00 * x2[i+j*(m-1)];

    d = ( ityp *  ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
    {
        d[j] = 0.00;
        for ( i = 0; i < m - 1; ++i )
        d[j] += pow ( x2[i+j*(m-1)], 2 );
    }

    for ( j = 0; j < n; ++j )
        x[m-1+j*m] = d[j] - 1.00;

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < m; ++i )
            x[i+j*m] /= ( d[j] + 1.00 );

    free ( d );
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hypersphere_surface_uniform ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_SURFACE_UNIFORM: uniform hypersphere surface samples
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 May 2013
  Author:
    John Burkardt
  Reference:
    Russell Cheng,
    Random Variate Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998, pages 168.
    George Marsaglia,
    Choosing a point from the surface of a sphere,
    Annals of Mathematical Statistics,
    Volume 43, Number 2, April 1972, pages 645-646.
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Wiley, 1986, page 234.
  Parameters:
    Input, int M, the dimension of the space.
    Input, int N, the number of points.
    Input, real R, the radius.
    Input, real C[M], the center.
    Input/output, int *SEED, a seed for the random number generator.
    Output, real HYPERSPHERE_SURFACE_UNIFORM[M*N], the points.
*/
{
	const _2dtitpitpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register ityp r = s_data->a2;
	ityp * c = s_data->a3;
	int * seed = s_data->a4;
	
    dim_typ i, j;
    ityp *x;

    x = hypersphere_01_surface_uniform ( m, n, seed );
    /*
    Scale by the radius.
    */
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            x[i+j*m] *= r;
    /*
    Shift to the center.
    */
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            x[i+j*m] += c[i];

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hypersphere_to_cartesian ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_TO_CARTESIAN: hypersphere to Cartesian coordinate transform.
  Discussion:
    We allow the trivial case M = 1; in that case alone, the value R
    must be assumed to have a sign.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, integer M, the spatial dimension.
    1 <= M.
    Input, integer N, the number of points to transform.
    Input, real C[M], the center of the hypersphere.
    Input, real R[N], the radius of the points on the hypersphere.
    Except for the trivial case M = 1, R is assumed nonnegative.
    Input, real THETA[(M-1)*N], the coordinate angles of the points,
    measured in radians.
    Output, real HYPERSPHERE_TO_CARTESIAN[M*N], the Cartesian
    coordinates of the points.
*/
{
	const _2dt3pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * c = s_data->a2;
	ityp * r = s_data->a3;
	ityp * theta = s_data->a4;
	
    dim_typ i, i1, i2, j;
    ityp *x;

    x = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    if ( m == 1 )
    {
        for ( j = 0; j < n; ++j )
            x[0+j*m] = r[j];
    }
    else
    {
        for ( j = 0; j < n; ++j )
            for ( i = 0; i < m; ++i )
            x[i+j*m] = r[j];

        for ( j = 0; j < n; ++j )
            for ( i1 = 0; i1 < m - 1; ++i1 )
            {
                x[i1+j*m] = x[i1+j*m] * cos ( theta[i1+j*(m-1)] );
                for ( i2 = i1 + 1; i2 < m; ++i2 )
                    x[i2+j*m] *= sin ( theta[i1+j*(m-1)] );
            }

    }
    /*
    Add the center.
    */
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            x[i+j*m] += c[i];

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hypersphere_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_VOLUME computes the volume of a hypersphere.
  Discussion:
    A hypersphere satisfies the equation:
      sum ( ( X(1:N) - PC(1:N) )^2 ) = R^2
    where R is the radius and PC is the center.
    Results include:
    M     Volume
    -     -----------------------
    2                M_PI   * R^2
    3  (4/3)    * M_PI   * R^3
    4  (1/2)    * M_PI^2 * R^4
    5  (8/15)   * M_PI^2 * R^5
    6  (1/6)    * M_PI^3 * R^6
    7  (16/105) * M_PI^3 * R^7
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2003
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of the space.
    Input, double R, the radius.
    Output, double HYPERSPHERE_VOLUME, the volume.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register ityp r = s_data->a1;
	
	result = pow ( r, m ) * hypersphere_01_volume ( m );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_stereograph ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_STEREOGRAPH computes the stereographic image of points on a sphere.
  Discussion:
    We start with a sphere of radius 1 and center (0,0,0).
    The north pole N = (0,0,1) is the point of tangency to the sphere
    of a plane, and the south pole S = (0,0,-1) is the focus for the
    stereographic projection.
    For any point P on the sphere, the stereographic projection Q of the
    point is defined by drawing the line from S through P, and computing
    Q as the intersection of this line with the plane.
    Actually, we allow the spatial dimension M to be arbitrary.  Values
    of M make sense starting with 2.  The north and south poles are
    selected as the points (0,0,...,+1) and (0,0,...,-1).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2010
  Author:
    John Burkardt
  Reference:
    C F Marcus,
    The stereographic projection in vector notation,
    Mathematics Magazine,
    Volume 39, Number 2, March 1966, pages 100-102.
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of points.
    Input, double P[M*N], a set of points on the unit sphere.
    Output, double SPHERE_STEREOGRAPH[M*N], the coordinates of the
    image points.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * p = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    ityp *q = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
    {
        for ( i = 0; i < m - 1; ++i )
            q[i+j*m] = 2.00 * p[i+j*m] / ( 1.00 + p[m-1+j*m] );
        q[m-1+j*m] = 1.00;
    }

    return q;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_stereograph_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_STEREOGRAPH_INVERSE computes stereographic preimages of points.
  Discussion:
    We start with a sphere of radius 1 and center (0,0,0).
    The north pole N = (0,0,1) is the point of tangency to the sphere
    of a plane, and the south pole S = (0,0,-1) is the focus for the
    stereographic projection.
    For any point Q on the plane, the stereographic inverse projection
    P of the point is defined by drawing the line from S through Q, and
    computing P as the intersection of this line with the sphere.
    Actually, we allow the spatial dimension M to be arbitrary.  Values
    of M make sense starting with 2.  The north and south poles are
    selected as the points (0,0,...,+1) and (0,0,...,-1).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2010
  Author:
    John Burkardt
  Reference:
    C F Marcus,
    The stereographic projection in vector notation,
    Mathematics Magazine,
    Volume 39, Number 2, March 1966, pages 100-102.
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of points.
    Input, double Q[M*N], the points, which are presumed to lie
    on the plane Z = 1.
    Output, double SPHERE_STEREOGRAPH_INVERSE[M*N], the stereographic
    inverse projections of the points.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * q = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    ityp *p = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );
    ityp qn;


    for ( j = 0; j < n; ++j )
    {
        qn = 0.00;
        for ( i = 0; i < m - 1; ++i )
            qn += pow ( q[i+j*m], 2 );
        for ( i = 0; i < m - 1; ++i )
            p[i+j*m] = 4.00 * q[i+j*m] / ( 4.00 + qn );
        p[m-1+j*m] = ( 4.00 - qn ) / ( 4.00 + qn );
    }

    return p;
}

#endif
