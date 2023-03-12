#ifndef __DISABLEDEEP_TETRAHEDRONGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tetrahedron_grid ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_GRID computes points on a tetrahedral grid.
  Discussion:
    The grid is defined by specifying the coordinates of an enclosing
    tetrahedron T, and the number of subintervals each edge of the
    tetrahedron should be divided into.
    Choosing N = 10, for instance, breaks each side into 10 subintervals,
    and produces a grid of ((10+1)*(10+2)*(10+3))/6 = 286 points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of subintervals.
    Input, double T[3*4], the vertices of the tetrahedron.
    Input, int NG, the number of grid points.
    Output, double TETRAHEDRON_GIRD[3*NG], the tetrahedron grid points.
*/
{
	const _2dtpit * const s_data = data; 
	
	const register dim_typ n = s_data->a0;
	const register dim_typ ng = s_data->a1;
	ityp * t = s_data->a2;
	
    dim_typ i;
    dim_typ ii;
    dim_typ j;
    dim_typ k;
    dim_typ l;
    dim_typ p;
    ityp *tg = ( ityp * ) malloc ( 3 * ng * sizeof ( double ) );
    p = 0;

    for ( i = 0; i <= n; ++i )
        for ( j = 0; j <= n - i; ++j )
            for ( k = 0; k <= n - i - j; ++k )
            {
                l = n - i - j - k;
                #pragma omp parallel for num_threads(3)
                for ( ii = 0; ii < 3; ++ii )
                    tg[ii+p*3] = ( ( ityp ) ( i ) * t[ii+0*3]+ ( ityp ) ( j ) * t[ii+1*3]+ ( ityp ) ( k ) * t[ii+2*3]+ ( ityp ) ( l ) * t[ii+3*3] ) / ( ityp ) ( n );
                ++ p;
            }

    return tg;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_grid_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_GRID_COUNT counts the grid points inside a tetrahedron.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of subintervals.
    Output, int TETRAHEDRON_GRID_COUNT, the number of grid points.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ( ( n + 1 ) * ( n + 2 ) * ( n + 3 ) ) / 6;
    return &result;
}

#endif
