#ifndef __DISABLEDEEP_SQUAREGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _square_grid ( void * data)
/******************************************************************************/
/*
  Purpose:
    SQUARE_GRID: grid points over the interior of a square in 2D.
  Discussion:
    In 2D, a logically rectangular grid is to be created.
    In the I-th dimension, the grid will use S(I) points.
    The total number of grid points is
      N = product ( 1 <= I <= 2 ) S(I)
    Over the interval [A(i),B(i)], we have 5 choices for grid centering:
      1: 0,   1/3, 2/3, 1
      2: 1/5, 2/5, 3/5, 4/5
      3: 0,   1/4, 2/4, 3/4
      4: 1/4, 2/4, 3/4, 1
      5: 1/8, 3/8, 5/8, 7/8
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N = product ( 1 <= I <= 2 ) NS(I).
    Input, int NS[2], the number of points along
    each dimension.
    Input, double A[2], B[2], the endpoints for each dimension.
    Input, int C[2], the grid centering for each dimension.
    1 <= C(*) <= 5.
    Output, double SQUARE_GRID[2*N] = X(2*S(0)*S(1)), the points.
*/
{
	const dtpi2pitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * ns = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	int * c = s_data->a4;
	
    dim_typ i, j;
    static int m = 2;
    dim_typ s;
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

        for ( j = 0; j < s; ++j )
        {
            if ( c[i] == 1 )
                xs[j] = s == 1 ? 0.50 * ( a[i] + b[i] ) : ( ( ityp ) ( s - j - 1 ) * a[i]+ ( ityp ) (     j     ) * b[i] )/ ( ityp ) ( s     - 1 );
            else if ( c[i] == 2 )
                xs[j] = ( ( ityp ) ( s - j     ) * a[i]+ ( ityp ) (     j + 1 ) * b[i] )/ ( ityp ) ( s     + 1 );
            else if ( c[i] == 3 )
                xs[j] = ( ( ityp ) ( s - j     ) * a[i]+ ( ityp ) (     j - 2 ) * b[i] )/ ( ityp ) ( s         );
            else if ( c[i] == 4 )
                xs[j] = ( ( ityp ) ( s - j - 1 ) * a[i]+ ( ityp ) (     j + 1 ) * b[i] )/ ( ityp ) ( s         );
            else if ( c[i] == 5 )
                xs[j] = ( ( ityp ) ( (s<<1)- (j<<1) - 1 ) * a[i]+ ( ityp ) (      (j<<1) + 1 ) * b[i] )/ ( ityp ) ( (s<<1)             );
        }

        r8vec_direct_product ( i, s, xs, m, n, x );

        free ( xs );
    }

    return x;
}

#endif
