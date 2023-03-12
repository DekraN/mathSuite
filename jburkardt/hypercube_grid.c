#ifndef __DISABLEDEEP_HYPERCUBEGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hypercube_grid ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERCUBE_GRID: grid points over the interior of a hypercube in M dimensions.
  Discussion:
    In M dimensional space, a logically rectangular grid is to be created.
    In the I-th dimension, the grid will use S(I) points.
    The total number of grid points is
      N = product ( 1 <= I <= M ) S(I)
    Over the interval [A(i),B(i)], we have 5 choices for grid centering:
      1: 0,   1/3, 2/3, 1
      2: 1/5, 2/5, 3/5, 4/5
      3: 0,   1/4, 2/4, 3/4
      4: 1/4, 2/4, 3/4, 1
      5: 1/8, 3/8, 5/8, 7/8
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of points.
    N = product ( 1 <= I <= M ) NS(I).
    Input, int NS[M], the number of points along
    each dimension.
    Input, double A[M], B[M], the endpoints for each dimension.
    Input, int C[M], the grid centering for each dimension.
    1 <= C(*) <= 5.
    Output, double HYPERCUBE_GRID[M*N] = X(M*S(1),S(2),...,S(M)), the points.
*/
{
	const _2dtpdt2pitpdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	dim_typ * ns = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	dim_typ * c = s_data->a5;
	
    dim_typ i, j, s;
    ityp *x;
    ityp *xs;

    x = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );
    /*
    Create the 1D grids in each dimension.
    */
    for ( i = 0; i < m; ++i )
    {
        s = ns[i];

        xs = ( ityp * ) malloc ( s * sizeof ( ityp ) );

        for ( j = 0; j < s; j++ )
        {
            if ( c[i] == 1 )
            {
                if ( s == 1 )
                    xs[j] = 0.50 * ( a[i] + b[i] );
                else
                    xs[j] = ( ( ityp ) ( s - j - 1 ) * a[i]   + ( ityp ) (     j     ) * b[i] ) / ( ityp ) ( s     - 1 );
            }
            else if ( c[i] == 2 )
                xs[j] = ( ( ityp ) ( s - j     ) * a[i]   + ( ityp ) (     j + 1 ) * b[i] ) / ( ityp ) ( s     + 1 );
            else if ( c[i] == 3 )
                xs[j] = ( ( ityp ) ( s - j     ) * a[i]   + ( ityp ) (     j - 2 ) * b[i] ) / ( ityp ) ( s         );
            else if ( c[i] == 4 )
                xs[j] = ( ( ityp ) ( s - j - 1 ) * a[i]   + ( ityp ) (     j + 1 ) * b[i] ) / ( ityp ) ( s         );
            else if ( c[i] == 5 )
                xs[j] = ( ( ityp ) ( (s<<1) - (j<<1) - 1 ) * a[i]   + ( ityp ) (      (j<<1) + 1 ) * b[i] ) / ( ityp ) ( s<<1             );
        }

        r8vec_direct_product ( i, s, xs, m, n, x );
        free ( xs );
    }

    return x;
}

#endif
