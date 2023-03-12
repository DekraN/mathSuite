#ifndef __DISABLEDEEP_TOMS097

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_shortest_path ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_SHORTEST_PATH computes the shortest distance between all pairs of points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 2014
  Author:
    John Burkardt
  Reference:
    Robert Floyd,
    Algorithm 97, Shortest Path,
    Communications of the ACM,
    Volume 5, Number 6, June 1962, page 345.
  Parameters:
    Input, int N, the number of points.
    Input/output, int M[N*N].
    On input, M(I,J) contains the length of the direct link between
    nodes I and J, or HUGE if there is no direct link.
    On output, M(I,J) contains the distance between nodes I and J,
    that is, the length of the shortest path between them.  If there
    is no such path, then M(I,J) will remain HUGE.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * m = s_data->a1;
	
    dim_typ i, j, k, s;

    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j )
            if ( m[j+i*n] < 2147483647 )
                for ( k = 0; k < n; ++k)
                    if ( m[i+k*n] < 2147483647 )
                    {
                        s = m[j+i*n] + m[i+k*n];
                        if ( s < m[j+k*n] )
                            m[j+k*n] = s;
                    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8mat_shortest_path ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_SHORTEST_PATH computes the shortest distance between all pairs of points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 2014
  Author:
    John Burkardt
  Reference:
    Robert Floyd,
    Algorithm 97, Shortest Path,
    Communications of the ACM,
    Volume 5, Number 6, June 1962, page 345.
  Parameters:
    Input, int N, the number of points.
    Input/output, double M[N*N].
    On input, M(I,J) contains the length of the direct link between
    nodes I and J, or HUGE if there is no direct link.
    On output, M(I,J) contains the distance between nodes I and J,
    that is, the length of the shortest path between them.  If there
    is no such path, then M(I,J) will remain HUGE.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * m = s_data->a1;
	
    dim_typ i, j, k;
    ityp s;

    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j)
            if ( m[j+i*n] < 1.0E+30 )
                for ( k = 0; k < n; ++k )
                    if ( m[i+k*n] < 1.0E+30 )
                    {
                        s = m[j+i*n] + m[i+k*n];
                        if ( s < m[j+k*n] )
                            m[j+k*n] = s;
                    }
    return NULL;
}

#endif
