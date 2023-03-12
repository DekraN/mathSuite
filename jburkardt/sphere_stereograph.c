#ifndef __DISABLEDEEP_SPHERESTEREOGRAPH

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_stereograph2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_STEREOGRAPH2 computes the stereographic image of points on a sphere.
  Discussion:
    We start with a sphere of center C.
    F is a point on the sphere which is the focus of the mapping,
    and the antipodal point 2*C-F is the point of tangency
    to the sphere of a plane.
    For any point P on the sphere, the stereographic projection Q of the
    point is defined by drawing the line from F through P, and computing
    Q as the intersection of this line with the plane.
    The spatial dimension M is arbitrary, but should be at least 2.
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
    Input, double FOCUS[M], the coordinates of the focus point.
    Input, double CENTER[M], the coordinates of the center of the sphere.
    Output, double SPHERE_STEREOGRAPH2[M*N], the coordinates of the
    image points,
*/
{
	const _2dt3pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * p = s_data->a2;
	ityp * focus = s_data->a3;
	ityp * center = s_data->a4;
	
    ityp cf_dot_pf;
    ityp cf_normsq;
    dim_typ i, j;
    ityp *q=( ityp * ) malloc ( m * n * sizeof ( ityp ) );
    ityp s;

    for ( j = 0; j < n; ++j )
    {
        cf_normsq = cf_dot_pf = 0.00;
        for ( i = 0; i < m; ++i )
        {
            cf_normsq += pow ( center[i] - focus[i], 2 );
            cf_dot_pf += ( center[i] - focus[i] ) * ( p[i+j*m] - focus[i] );
        }
        s = 2.00 * cf_normsq / cf_dot_pf;
        for ( i = 0; i < m; ++i )
            q[i+j*m] = s * p[i+j*m] + ( 1.00 - s ) * focus[i];
    }
    return q;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere_stereograph2_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE_STEREOGRAPH2_INVERSE computes stereographic preimages of points.
  Discussion:
    We start with a sphere of center C.
    F is a point on the sphere which is the focus of the mapping,
    and the antipodal point 2*C-F is the point of tangency
    to the sphere of a plane.
    For any point Q on the plane, the stereographic inverse projection
    P of the point is defined by drawing the line from F through Q, and
    computing P as the intersection of this line with the sphere.
    The spatial dimension M is arbitrary, but should be at least 2.
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
    on the plane.
    Input, double FOCUS[M], the coordinates of the focus point.
    Input, double CENTER[M], the coordinates of the center of the sphere.
    Output, double SPHERE_STEREOGRAPH2_INVERSE[M*N], the stereographic
    inverse projections of the points.
*/
{
	const _2dt3pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * q = s_data->a2;
	ityp * focus = s_data->a3;
	ityp * center = s_data->a4;
	
    ityp cf_dot_qf;
    dim_typ i, j;
    ityp *p = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );
    ityp qf_normsq;
    ityp s;

    for ( j = 0; j < n; ++j )
    {
        cf_dot_qf = qf_normsq = 0.00;
        for ( i = 0; i < m; ++i )
        {
            cf_dot_qf += ( center[i] - focus[i] ) * ( q[i+j*m] - focus[i] );
            qf_normsq += pow ( q[i+j*m] - focus[i], 2 );
        }
        s = 2.0 * cf_dot_qf / qf_normsq;
        for ( i = 0; i < m; ++i )
            p[i+j*m] = s * q[i+j*m] + ( 1.00 - s ) * focus[i];
    }
    return p;
}

#endif
