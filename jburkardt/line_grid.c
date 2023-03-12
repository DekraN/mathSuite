#ifndef __DISABLEDEEP_LINEGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _line_grid ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_GRID: grid points over the interior of a line segment in 1D.
  Discussion:
    In 1D, a grid is to be created using N points.
    Over the interval [A(i),B(i)], we have 5 choices for grid centering:
      1: 0,   1/3, 2/3, 1
      2: 1/5, 2/5, 3/5, 4/5
      3: 0,   1/4, 2/4, 3/4
      4: 1/4, 2/4, 3/4, 1
      5: 1/8, 3/8, 5/8, 7/8
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, double A, B, the endpoints.
    Input, int C, the grid centering.
    1 <= C <= 5.
    Output, double LINE_GRID[N], the points.
*/
{
	const dt2itdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp a = s_data->a1;
	ityp b = s_data->a2;
	const register dim_typ c = s_data->a3;
	
    ityp *x = ( double * ) malloc ( n * sizeof ( double ) );
    /*
    Create the 1D grids in each dimension.
    */
    for (dim_typ j = 0; j < n; ++j )
    {
        if ( c == 1 )
            x[j] = n == 1 ? 0.50 * ( a + b ) : ( ( ityp ) ( n - j - 1 ) * a   + ( ityp ) (     j     ) * b ) / ( ityp ) ( n    - 1 );
        else if ( c == 2 )
            x[j] = ( ( ityp ) ( n - j     ) * a   + ( ityp ) (     j + 1 ) * b ) / ( ityp ) ( n     + 1 );
        else if ( c == 3 )
            x[j] = ( ( ityp ) ( n - j     ) * a   + ( ityp ) (     j - 2 ) * b ) / ( ityp ) ( n         );
        else if ( c == 4 )
            x[j] = ( ( ityp ) ( n - j - 1 ) * a   + ( ityp ) (     j + 1 ) * b ) / ( ityp ) ( n         );
        else if ( c == 5 )
            x[j] = ( ( ityp ) ( (n<<1) - (j<<1) - 1 ) * a   + ( ityp ) (      (j<<1) + 1 ) * b ) / ( ityp ) ( n<<1            );
    }

    return x;
}

#endif
