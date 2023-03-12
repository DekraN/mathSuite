#ifndef __DISABLEDEEP_SIMPLEXGRID

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ksub_random_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    KSUB_RANDOM selects a random subset of size K from a set of size N.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    30 April 2003
  Author:
    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the size of the set from which subsets are drawn.
    Input, int K, number of elements in desired subsets.  K must
    be between 0 and N.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int KSUB_RANDOM_NEW[K].  A(I) is the I-th element of the
    output set.  The elements of A are in order.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	int * seed = s_data->a2;
	
	int *a;
	dim_typ i;
	int ids;
	int ihi;
	dim_typ ip;
	int ir;
	int is;
	dim_typ ix;
	int l;
	int ll;
	dim_typ m;
	int m0;
	
	if ( n < k )
		return NULL;
	
	a = ( int * ) malloc ( k * sizeof ( int ) );
	
	if ( k == 0 )
		return a;
	
	for ( i = 1; i <= k; ++i )
		a[i-1] = ( ( i - 1 ) * n ) / k;
	
	for ( i = 1; i <= k; ++i )
	{
		for ( ; ; )
		{
			ix = i4_uniform_ab ( 1, n, seed );
			
			l = 1 + ( ix * k - 1 ) / n;
			
			if ( a[l-1] < ix )
				break;
		}
	
		++ a[l-1];
	
	}
	
	ip = 0;
	is = k;
	
	for ( i = 1; i <= k; ++i )
	{
		m = a[i-1];
		a[i-1] = 0;
		
		if ( m != ( ( i - 1 ) * n ) / k )
		{
			++ ip;
			a[ip-1] = m;
		}
	
	}
	
	ihi = ip;
	
	for ( i = 1; i <= ihi; ++i )
	{
		ip = ihi + 1 - i;
		l = 1 + ( a[ip-1] * k - 1 ) / n;
		ids = a[ip-1] - ( ( l - 1 ) * n ) / k;
		a[ip-1] = 0;
		a[is-1] = l;
		is = is - ids;
	}
	
	for ( ll = 1; ll <= k; ++ll )
	{
		l = k + 1 - ll;
		
		if ( a[l-1]  )
		{
			ir = l;
			m0 = 1 + ( ( a[l-1] - 1 ) * n ) / k;
			m = ( a[l-1] * n ) / k - m0 + 1;
		}
		
		ix = i4_uniform_ab ( m0, m0 + m - 1, seed );
		
		i = l + 1;
		
		while ( i <= ir )
		{
			if ( ix < a[i-1] )
				break;
	
			++ ix;
			a[i-2] = a[i-1];
			++ i;
		}
		a[i-2] = ix;
		-- m;
	}
	return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _comp_random_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMP_RANDOM selects a random composition of the integer N into K parts.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    30 April 2003
  Author:
    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the integer to be decomposed.
    Input, int K, the number of parts in the composition.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int COMP_RANDOM_NEW[K], the parts of the composition.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	int * seed = s_data->a2;
	
	int *a;
	int l;
	int m;
	
	a = ksub_random_new ( n+k-1, k-1, seed );
	
	a[k-1] = n + k;
	l = 0;
	
	for (dim_typ i = 0; i < k; ++i )
	{
		m = a[i];
		a[i] -= l + 1;
		l = m;
	}
	
	return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _simplex_grid_index_all ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_GRID_INDEX_ALL returns all the simplex grid indices.
  Discussion:
    The number of grid indices can be determined by calling
      ng = simplex_grid_size ( m, n )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of subintervals.
    Input, int NG, the number of values in the grid.
    Output, int SIMPLEX_GRID_INDEX_ALL[(M+1)*NG], the current, and then the next,
    grid index.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	const register dim_typ ng = a_data[2];
	
    dim_typ *g;
    int *grid;
    dim_typ i, k;

    g = ( dim_typ * ) malloc ( ( m + 1 ) * sizeof ( dim_typ ) );

    for ( i = 0; i < m; ++i )
        g[i] = 0;
    g[m] = n;

    grid = ( int * ) malloc ( ( m + 1 ) * ng * sizeof ( int ) );

    k = 0;
    for ( i = 0; i <= m; ++i)
        grid[i+k*(m+1)] = g[i];

    while ( k < ng )
    {
        comp_next_grlex ( m + 1, g );
        ++ k;
        for ( i = 0; i <= m; ++i )
            grid[i+k*(m+1)] = g[i];
    }

    free ( g );

    return grid;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _simplex_grid_index_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_GRID_INDEX_NEXT returns the next simplex grid index.
  Discussion:
    The vector G has dimension M+1.  The first M entries may be regarded
    as grid coordinates.  These coordinates must have a sum between 0 and N.
    The M+1 entry contains the remainder, that is N minus the sum of the
    first M coordinates.
    Each time the function is called, it is given a current grid index, and
    computes the next one.  The very first index is all zero except for a
    final value of N, and the very last index has all zero except for an'
    intial value of N.
    For example, here are the coordinates in order for M = 3, N = 3:
     0  0 0 0 3
     1  0 0 1 2
     2  0 0 2 1
     3  0 0 3 0
     4  0 1 0 2
     5  0 1 1 1
     6  0 1 2 0
     7  0 2 0 1
     8  0 2 1 0
     9  0 3 0 0
    10  1 0 0 2
    11  1 0 1 1
    12  1 0 2 0
    13  1 1 0 1
    14  1 1 1 0
    15  1 2 0 0
    16  2 0 0 1
    17  2 0 1 0
    18  2 1 0 0
    19  3 0 0 0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of subintervals.
    Input/output, int G[M+1], the current, and then the next,
    grid index.
*/
{
	const _2dtpdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	dim_typ * g = s_data->a2;
	
    comp_next_grlex ( m + 1, g );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _simplex_grid_index_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_GRID_INDEX_SAMPLE returns a random simplex grid index.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of subintervals in
    each dimension.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int SIMPLEX_GRID_INDEX_SAMPLE[M+1], a randomly selected index
    in the simplex grid.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
    return comp_random_new ( n, m + 1, seed );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _simplex_grid_index_to_point ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_GRID_INDEX_TO_POINT returns  points corresponding to simplex indices.
  Discussion:
    The M-dimensional simplex is defined by M+1 vertices.
    Given a regular grid that uses N subintervals along the edge between
    each pair of vertices, a simplex grid index G is a set of M+1 values
    each between 0 and N, and summing to N.
    This function determines the coordinates X of the point corresponding
    to the index G.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of subintervals.
    Input, int NG, the number of grid indices to be converted.
    Input, int G[(M+1)*NG], the grid indices of 1
    or more points.
    Input, double V[M*(M+1)], the coordinates of the vertices
    of the simplex.
    Output, double SIMPLEX_GRID_INDEX_TO_POINT[M*NG], the coordinates of one
    or more points.
*/
{
	const _3dtpipit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register dim_typ ng = s_data->a2;
	int * g = s_data->a3;
	ityp * v = s_data->a4;
	
    dim_typ i, j, k;
    double *x;

    x = ( ityp * ) malloc ( m * ng * sizeof ( ityp ) );

    for ( j = 0; j < ng; ++j )
    {
        for ( i = 0; i < m; ++i )
        {
            x[i+j*m] = 0.0;
            for ( k = 0; k < m + 1; ++k )
                x[i+j*m] += v[i+k*m] * ( ityp ) ( g[k+j*(m+1)] );
            x[i+j*m] /= ( ityp ) ( n );
        }
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _simplex_grid_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_GRID_SIZE counts the grid points inside a simplex.
  Discussion:
    The size of a grid with parameters M, N is C(M+N,N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of subintervals.
    Output, int SIMPLEX_GRID_SIZE, the number of grid points.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
    int ng = 1;
    for (dim_typ i = 1; i <= m; ++i )
        ng = ( ng * ( n + i ) ) / i;
        
    result = ng;
    return &result;
}

#endif
