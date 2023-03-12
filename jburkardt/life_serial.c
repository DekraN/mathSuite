#ifndef __DISABLEDEEP_LIFESERIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _life_init ( void * data)
/******************************************************************************/
/*
  Purpose:
    LIFE_INIT initializes the life grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, double PROB, the probability that a grid cell
    should be alive.
    Input, int M, N, the number of rows and columns
    of interior grid cells.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, int LIFE_INIT[(1+M+1)*(1+N+1)], the initial grid.
*/
{
	const _2dtitpi * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register ityp prob = s_data->a2;
	int * seed = s_data->a3;
	
    int *grid;
    dim_typ i, j;
    ityp r;

    grid = ( int * ) malloc ( ( m + 2 ) * ( n + 2 ) * sizeof ( int ) );
    for ( j = 0; j <= n + 1; ++j)
        for ( i = 0; i <= m + 1; ++i )
            grid[i+j*(m+2)] = 0;

    for ( j = 1; j <= n; ++j )
        for ( i = 1; i <= m; ++i )
        {
            r = r8_uniform_01 ( seed );
            if ( r <= prob )
                grid[i+j*(m+2)] = 1;
        }

    return grid;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _life_update ( void * data)
/******************************************************************************/
/*
  Purpose:
    LIFE_UPDATE updates a Life grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns
    of interior grid cells.
    Input/output, int GRID[(1+M+1)*(1+N+1)], the data.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * grid = s_data->a2;
	
    dim_typ i, j;
    int *s;

    s = ( int * ) malloc ( m * n * sizeof ( int ) );

    for ( j = 1; j <= n; ++j)
        for ( i = 1; i <= m; ++i)
            s[i-1+(j-1)*m] = grid[i-1+(j-1)*(m+2)] + grid[i-1+j*(m+2)] + grid[i-1+(j+1)*(m+2)]+ grid[i  +(j-1)*(m+2)]                     + grid[i  +(j+1)*(m+2)]+ grid[i+1+(j-1)*(m+2)] + grid[i+1+j*(m+2)] + grid[i+1+(j+1)*(m+2)];
    /*
    Any dead cell with 3 live neighbors becomes alive.
    Any living cell with less than 2 or more than 3 neighbors dies.
    */
    for ( j = 1; j <= n; ++j )
        for ( i = 1; i <= m; ++i )
            if ( grid[i+j*(m+2)] == 0 )
            {
                if ( s[i-1+(j-1)*m] == 3 )
                    grid[i+j*(m+2)] = 1;
            }
            else if ( grid[i+j*(m+2)] == 1 )
            {
                    if ( s[i-1+(j-1)*m] < 2 || 3 < s[i-1+(j-1)*m] )
                        grid[i+j*(m+2)] = 0;
            }

    free ( s );
    return NULL;
}

#endif

