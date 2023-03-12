#ifndef __DISABLEDEEP_TRIANGLEGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangle_grid ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_GRID computes points on a triangular grid.
  Discussion:
    The grid is defined by specifying the coordinates of an enclosing
    triangle T, and the number of subintervals each side of the triangle
    should be divided into.
    Choosing N = 10, for instance, breaks each side into 10 subintervals,
    and produces a grid of ((10+1)*(10+2))/2 = 66 points.
              X
             9 X
            8 9 X
           7 8 9 X
          6 7 8 9 X
         5 6 7 8 9 X
        4 5 6 7 8 9 X
       3 4 5 6 7 8 9 X
      2 3 4 5 6 7 8 9 X
     1 2 3 4 5 6 7 8 9 X
    0 1 2 3 4 5 6 7 8 9 X
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of subintervals.
    Input, double T[2*3], the coordinates of the points
    defining the triangle.
    Output, double *TRIANGLE_GRID[2*((N+1)*(N+2))/2], the coordinates
    of the points in the triangle.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * t = s_data->a1;
	
    dim_typ i;
    ityp ir;
    dim_typ j;
    ityp jr;
    dim_typ k;
    ityp kr;
    dim_typ l;
    dim_typ ng;
    ityp nr;
    dim_typ p;
    ityp *tg;

    ng = ( ( n + 1 ) * ( n + 2 ) ) / 2;
    tg = ( ityp * ) malloc ( ng * sizeof ( ityp ) << 1 );

    p = 0;
    nr = ( ityp ) ( n );

    for ( i = 0; i <= n; ++i )
    {
        ir = ( ityp ) ( i );
        for ( j = 0; j <= n - i; ++j )
        {
            jr = ( ityp ) ( j );
            k = n - i - j;
            kr = ( ityp ) ( k );
            #pragma omp parallel for num_threads(2)
            for ( l = 0; l < 2; ++l )
                tg[l+(p<<1)] = ( ir * t[l] + jr * t[l+2] + kr * t[l+4] ) / nr;
            ++ p;
        }
    }

    return tg;
}

#endif
