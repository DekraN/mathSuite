#ifndef __DISABLEDEEP_LAGRANGEINTERPND

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void    * _cc_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    CC_COMPUTE_POINTS computes Clenshaw Curtis quadrature points.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    This rule is defined on [-1,1].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the rule.
    Output, double CC_COMPUTE_POINTS[N], the abscissas.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    ityp *x;

    if ( n < 1 )
        return NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
        x[0] = 0.00;
    else
    {
        for ( i = 1; i <= n; ++i )
            x[i-1] =  cos ( ( ityp ) ( n - i ) * M_PI / ( ityp ) ( n - 1     ) );
        x[0] = -1.00;
        if ( ( n % 2 ) == 1 )
            x[(n-1)/2] = 0.00;
        x[n-1] = +1.00;
    }
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_interp_nd_grid ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N_1D[M], the order of the 1D rule to be used
    in each dimension.
    Input, double A[M], B[M], the upper and lower limits.
    Input, int ND, the number of points in the product grid.
    Output, double LAGRANGE_INTERP_ND_GRID[M*ND], the points at which data
    is to be sampled.
*/
{
	const pit2dtpdtpit * const s_data = data;
	
	ityp * a = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ nd = s_data->a2;
	dim_typ * n_1d = s_data->a3;
	ityp * b = s_data->a4;

	
    dim_typ i, j, n;
    ityp *x_1d;
    ityp *xd;
    /*
    Compute the data points.
    */
    xd = ( ityp * ) malloc ( m * nd * sizeof ( ityp ) );

    for ( j = 0; j < nd; ++j )
        for ( i = 0; i < m; ++i )
            xd[i+j*m] = 0.00;

    for ( i = 0; i < m; ++i )
    {
        n = n_1d[i];
        x_1d = cc_compute_points ( n );
        for ( j = 0; j < n; ++j )
            x_1d[j] = 0.50 * ( ( 1.00 - x_1d[j] ) * a[i]+ ( 1.00 + x_1d[j] ) * b[i] );
        r8vec_direct_product ( i, n, x_1d, m, nd, xd );
        free ( x_1d );
    }

    return xd;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_interp_nd_grid2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int IND[M], the index or level of the 1D rule
    to be used in each dimension.
    Input, double A[M], B[M], the upper and lower limits.
    Input, int ND, the number of points in the product grid.
    Output, double LAGRANGE_INTERP_ND_GRID2[M*ND], the points at which data
    was sampled.
*/
{
	const pit2dtpdtpit * const s_data = data;
	
	ityp * a = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ nd = s_data->a2;
	dim_typ * ind = s_data->a3;
	ityp * b = s_data->a4;
	
	
    dim_typ i, j, n;
    ityp *x_1d;
    ityp *xd;
    /*
    Compute the data points.
    */
    xd = ( ityp * ) malloc ( m * nd * sizeof ( ityp ) );

    for ( j = 0; j < nd; ++j)
        for ( i = 0; i < m; ++i )
            xd[i+j*m] = 0.00;

    for ( i = 0; i < m; ++i )
    {
        n = order_from_level_135 ( ind[i] );
        x_1d = cc_compute_points ( n );
        for ( j = 0; j < n; ++j )
            x_1d[j] = 0.50 * ( ( 1.00 - x_1d[j] ) * a[i]+ ( 1.00 + x_1d[j] ) * b[i] );
        r8vec_direct_product ( i, n, x_1d, m, nd, xd );
        free ( x_1d );
    }

    return xd;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lagrange_interp_nd_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N_1D[M], the order of the 1D rule to be used
    in each dimension.
    Output, int LAGRANGE_INTERP_ND_SIZE, the number of points in the product grid.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	int * n_1d = s_data->a1;
	
    /*
    Determine the number of data points.
    */
    result = i4vec_product ( m, n_1d );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lagrange_interp_nd_size2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int IND[M], the index or level of the 1D rule
    to be used in each dimension.

    Output, int LAGRANGE_INTERP_ND_SIZE2, the number of points in the product grid.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	dim_typ * ind = s_data->a1;
	
    /*
    Determine the number of data points.
    */
    dim_typ nd = 1;
    for (dim_typ i = 0; i < m; ++i )
        nd *= order_from_level_135 ( ind[i] );
        
    result = nd;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_interp_nd_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N_1D[M], the order of the 1D rule to be used
    in each dimension.
    Input, double A[M], B[M], the upper and lower limits.
    Input, int ND, the number of points in the product grid.
    Input, double ZD[ND], the function evaluated at the points XD.
    Input, int NI, the number of points at which the
    interpolant is to be evaluated.
    Input, double XI[M*NI], the points at which the interpolant is
    to be evaluated.
    Output, double LAGRANGE_INTERP_ND_VALUE[NI], the interpolant evaluated
    at the points XI.
*/
{
	const dtpdt2pitdtpitdtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	dim_typ * n_1d = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	const register dim_typ nd = s_data->a4;
	ityp * zd = s_data->a5;
	const register dim_typ ni = s_data->a6;
	ityp * xi = s_data->a7;
	
    dim_typ i, j, k, n;
    ityp *value;
    ityp *w;
    ityp *x_1d;
    ityp *zi;

    w = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    zi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( j = 0; j < ni; ++j)
    {
        for ( i = 0; i < nd; ++i )
            w[i] = 1.00;
        for ( i = 0; i < m; i++ )
        {
            n = n_1d[i];
            x_1d = cc_compute_points ( n );
            for ( k = 0; k < n; ++k )
                x_1d[k] = 0.50 * ( ( 1.00 - x_1d[k] ) * a[i]+ ( 1.00 + x_1d[k] ) * b[i] );
            value = lagrange_basis_1d ( n, x_1d, 1, xi+i+j*m );
            r8vec_direct_product2 ( i, n, value, m, nd, w );
            free ( value );
            free ( x_1d );
        }
        zi[j] = r8vec_dot_product ( nd, w, zd );
    }

    free ( w );
    return zi;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lagrange_interp_nd_value2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int IND[M], the index or level of the 1D rule
    to be used in each dimension.
    Input, double A[M], B[M], the upper and lower limits.
    Input, int ND, the number of points in the product grid.
    Input, double ZD[ND], the function evaluated at the points XD.
    Input, int NI, the number of points at which the
    interpolant is to be evaluated.
    Input, double XI[M*NI], the points at which the interpolant
    is to be evaluated.
    Output, double ZI[NI], the interpolant evaluated at the
    points XI.
*/
{
	const dtpdt2pitdtpitdtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	dim_typ * ind = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	const register dim_typ nd = s_data->a4;
	ityp * zd = s_data->a5;
	const register dim_typ ni = s_data->a6;
	ityp * xi = s_data->a7;
	
    dim_typ i, j, k, n;
    ityp *value;
    ityp *w;
    ityp *x_1d;
    ityp *zi;

    w = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    zi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( j = 0; j < ni; ++j)
    {
        for ( i = 0; i < nd; ++i)
            w[i] = 1.00;

        for ( i = 0; i < m; ++i )
        {
            n = order_from_level_135 ( ind[i] );
            x_1d = cc_compute_points ( n );
            for ( k = 0; k < n; ++k )
                x_1d[k] = 0.50 * ( ( 1.00 - x_1d[k] ) * a[i]+ ( 1.00 + x_1d[k] ) * b[i] );

            value = lagrange_basis_1d ( n, x_1d, 1, xi+i+j*m );
            r8vec_direct_product2 ( i, n, value, m, nd, w );
            free ( value );
            free ( x_1d );
        }
        zi[j] = r8vec_dot_product ( nd, w, zd );
    }

    free ( w );
    return zi;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _order_from_level_135 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
  Discussion:
    Clenshaw Curtis rules, and some others, often use the following
    scheme:
    L: 0  1  2  3   4   5
    N: 1  3  5  9  17  33 ... 2^L+1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int L, the level, which should be 0 or greater.
    Output, int ORDER_FROM_LEVEL_135, the order.
*/
{
	static int result = INT_MAX;
	
	const register int l = *(int *) data;
	
	result = l<0 ? INT_MAX : l == 0 ? 1:powi ( 2, l ) + 1;
    return &result;
}

#endif

