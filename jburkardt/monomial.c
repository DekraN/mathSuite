#ifndef __DISABLEDEEP_MONOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_between_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_BETWEEN_ENUM enumerates monomials in D dimensions of degrees in a range.
  Discussion:
    For D = 3, we have the following table:
     N2 0  1  2  3  4  5  6   7   8
   N1 +----------------------------
    0 | 1  4 10 20 35 56 84 120 165
    1 | 0  3  9 19 34 55 83 119 164
    2 | 0  0  6 16 31 52 80 116 161
    3 | 0  0  0 10 25 46 74 110 155
    4 | 0  0  0  0 15 36 64 100 145
    5 | 0  0  0  0  0 21 49  85 130
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N1, N2, the minimum and maximum degrees.
    0 <= N1 <= N2.
    Output, int MONO_BETWEEN_ENUM, the number of monomials
    in D variables, of total degree between N1 and N2 inclusive.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n1 = a_data[1];
	const register dim_typ n2 = a_data[2];
	
    dim_typ n0;
    dim_typ n1_copy;
    dim_typ value;

    n1_copy = MAX ( n1, 0 );

    if ( n2 < n1_copy )
    {
    	result = 0;
        return &result;
    }

    if ( n1_copy == 0 )

        value = i4_choose ( n2 + m, n2 );
    else if ( n1_copy == n2 )
        value = i4_choose ( n2 + m - 1, n2 );
    else
    {
        n0 = n1_copy - 1;
        value = i4_choose ( n2 + m, n2 ) - i4_choose ( n0 + m, n0 );
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_between_next_grevlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_BETWEEN_NEXT_GREVLEX: grevlex next monomial, total degree between N1 and N2.
  Discussion:
    We consider all monomials in an M-dimensional space, with total
    degree N between N1 and N2, inclusive.
    For example:
    M = 3
    N1 = 2
    N2 = 3
    #  X(1)  X(2)  X(3)  Degree
      +------------------------
    1 |  0     0     2        2
    2 |  0     1     1        2
    3 |  1     0     1        2
    4 |  0     2     0        2
    5 |  1     1     0        2
    6 |  2     0     0        2
      |
    7 |  0     0     3        3
    8 |  0     1     2        3
    9 |  1     0     2        3
   10 |  0     2     1        3
   11 |  1     1     1        3
   12 |  2     0     1        3
   13 |  0     3     0        3
   14 |  1     2     0        3
   15 |  2     1     0        3
   16 |  3     0     0        3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N1, N2, the minimum and maximum degrees.
    0 <= N1 <= N2.
    Input, int X[M], the current monomial.
    To start the sequence, set X = [ 0, 0, ..., 0, N1 ].
    Output, int X[M], the next monomial.
    The last value in the sequence is X = [ N2, 0, ..., 0, 0 ].
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n1 = s_data->a1;
	const register dim_typ n2 = s_data->a2;
	int * x = s_data->a3;
	
    dim_typ i;
    dim_typ j;
    dim_typ t;

    if ( n1 < 0 || n2 < n1 || i4vec_sum ( m, x ) < n1 || n2 < i4vec_sum ( m, x ) || n2 == 0 )
        return NULL;

    if ( x[0] == n2 )
    {
        x[0] = 0;
        x[m-1] = n1;
    }
    else
        mono_next_grevlex ( m, x );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_between_next_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_BETWEEN_NEXT_GRLEX: grlex next monomial, degree between N1 and N2.
  Discussion:
    We consider all monomials in an M-dimensional space, with total
    degree N between N1 and N2, inclusive.
    For example:
    M = 3
    N1 = 2
    N2 = 3
    #  X(1)  X(2)  X(3)  Degree
      +------------------------
    1 |  0     0     2        2
    2 |  0     1     1        2
    3 |  0     2     0        2
    4 |  1     0     1        2
    5 |  1     1     0        2
    6 |  2     0     0        2
      |
    7 |  0     0     3        3
    8 |  0     1     2        3
    9 |  0     2     1        3
   10 |  0     3     0        3
   11 |  1     0     2        3
   12 |  1     1     1        3
   13 |  1     2     0        3
   14 |  2     0     1        3
   15 |  2     1     0        3
   16 |  3     0     0        3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N1, N2, the minimum and maximum degrees.
    0 <= N1 <= N2.
    Input/output, int X[M], the current monomial.
    To start the sequence, set X = [ 0, 0, ..., 0, N1 ].
    The last value in the sequence is X = [ N2, 0, ..., 0, 0 ].
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n1 = s_data->a1;
	const register dim_typ n2 = s_data->a2;
	int * x = s_data->a3;
	
    dim_typ i;
    dim_typ im1;
    dim_typ j;
    dim_typ t;

    if ( n1 < 0 || n2 < n1 || i4vec_sum ( m, x ) < n1 || n2 < i4vec_sum ( m, x ) || n2 == 0  )
        return NULL;

    if ( x[0] == n2 )
    {
        x[0] = 0;
        x[m-1] = n1;
    }
    else
        mono_next_grlex ( m, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _mono_between_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_BETWEEN_RANDOM: random monomial with total degree between N1 and N2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N1, N2, the minimum and maximum degrees.
    0 <= N1 <= N2.
    Input/output, int *SEED, the random number seed.
    Output int *RANK, the rank of the monomial.
    Output int MONO_BETWEEN_RANDOM[M], the random monomial.
*/
{
	const _3dtpipdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n1 = s_data->a1;
	const register dim_typ n2 = s_data->a2;
	int * seed = s_data->a3;
	dim_typ * rank = s_data->a4;
	
    int n1_copy;
    dim_typ rank_max;
    dim_typ rank_min;
    int *x;
    n1_copy = MAX ( n1, 0 );
    rank_min = mono_upto_enum ( m, n1_copy - 1 ) + 1;
    rank_max = mono_upto_enum ( m, n2 );
    *rank = i4_uniform_ab ( rank_min, rank_max, seed );
    x = mono_unrank_grlex ( m, *rank );
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_next_grevlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_NEXT_GREVLEX: grevlex next monomial.
  Discussion:
    Example:
    M = 3
    #  X(1)  X(2)  X(3)  Degree
      +------------------------
    1 |  0     0     0        0
      |
    2 |  0     0     1        1
    3 |  0     1     0        1
    4 |  1     0     0        1
      |
    5 |  0     0     2        2
    6 |  0     1     1        2
    7 |  1     0     1        2
    8 |  0     2     0        2
    9 |  1     1     0        2
   10 |  2     0     0        2
      |
   11 |  0     0     3        3
   12 |  0     1     2        3
   13 |  1     0     2        3
   14 |  0     2     1        3
   15 |  1     1     1        3
   16 |  2     0     1        3
   17 |  0     3     0        3
   18 |  1     2     0        3
   19 |  2     1     0        3
   20 |  3     0     0        3
    Thanks to Stefan Klus for pointing out a discrepancy in a previous
    version of this code, 05 February 2015.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 February 2015
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input/output, int X[M], the current monomial.
    The first element is X = [ 0, 0, ..., 0, 0 ].
*/
{
	const dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	int * x = s_data->a1;
	
    dim_typ i;
    dim_typ j;
    dim_typ t;

    if ( i4vec_sum ( m, x ) < 0 )
        return NULL;
    /*
    Seeking the first index 1 < I for which 0 < X(I).
    */
    j = 0;

    for ( i = 1; i < m; ++i )
        if ( 0 < x[i] )
        {
            j = i;
            break;
        }

    if ( j == 0 )
    {
    t = x[0];
        x[0] = 0;
        x[m-1] = t + 1;
    }
    else if ( j < m - 1 )
    {
        x[j] = x[j] - 1;
        t = x[0] + 1;
        x[0] = 0;
        x[j-1] = x[j-1] + t;
    }
    else if ( j == m - 1 )
    {
        t = x[0];
        x[0] = 0;
        x[j-1] = t + 1;
        x[j] = x[j] - 1;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_next_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_NEXT_GRLEX returns the next monomial in grlex order.
  Discussion:
    Example:
    M = 3
    #  X(1)  X(2)  X(3)  Degree
      +------------------------
    1 |  0     0     0        0
      |
    2 |  0     0     1        1
    3 |  0     1     0        1
    4 |  1     0     0        1
      |
    5 |  0     0     2        2
    6 |  0     1     1        2
    7 |  0     2     0        2
    8 |  1     0     1        2
    9 |  1     1     0        2
   10 |  2     0     0        2
      |
   11 |  0     0     3        3
   12 |  0     1     2        3
   13 |  0     2     1        3
   14 |  0     3     0        3
   15 |  1     0     2        3
   16 |  1     1     1        3
   17 |  1     2     0        3
   18 |  2     0     1        3
   19 |  2     1     0        3
   20 |  3     0     0        3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input/output, int X[M], the current monomial.
    The first element is X = [ 0, 0, ..., 0, 0 ].
*/
{
	const dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	int * x = s_data->a1;
	
    dim_typ i;
    dim_typ im1;
    dim_typ j;
    dim_typ t;
    /*
    Ensure that 1 <= D.
    */
    if ( m < 1 )
        return NULL;
    /*
    Ensure that 0 <= X(I).
    */
    for ( i = 0; i < m; ++i )
        if ( x[i] < 0 )
            return NULL;
    /*
    Find I, the index of the rightmost nonzero entry of X.
    */
    i = 0;
    for ( j = m; 1 <= j; --j )
        if ( 0 < x[j-1] )
        {
            i = j;
            break;
        }
    /*
    set T = X(I)
    set X(I) to zero,
    increase X(I-1) by 1,
    increment X(M) by T-1.
    */
    if ( i == 0 )
    {
        x[m-1] = 1;
        return NULL;
    }
    else if ( i == 1 )
    {
        t = x[0] + 1;
        im1 = m;
    }
    else if ( 1 < i )
    {
        t = x[i-1];
        im1 = i - 1;
    }

    x[i-1] = 0;
    x[im1-1] = x[im1-1] + 1;
    x[m-1] = x[m-1] + t - 1;

    return NULL;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_rank_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_RANK_GRLEX computes the graded lexicographic rank of a monomial.
  Discussion:
    The graded lexicographic ordering is used, over all monomials in
    M dimensions, for total degree = 0, 1, 2, ...
    For example, if M = 3, the ranking begins:
    Rank  Sum    1  2  3
    ----  ---   -- -- --
       1    0    0  0  0
       2    1    0  0  1
       3    1    0  1  0
       4    1    1  0  1
       5    2    0  0  2
       6    2    0  1  1
       7    2    0  2  0
       8    2    1  0  1
       9    2    1  1  0
      10    2    2  0  0
      11    3    0  0  3
      12    3    0  1  2
      13    3    0  2  1
      14    3    0  3  0
      15    3    1  0  2
      16    3    1  1  1
      17    3    1  2  0
      18    3    2  0  1
      19    3    2  1  0
      20    3    3  0  0
      21    4    0  0  4
      ..   ..   .. .. ..
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    1 <= M.
    Input, int X[M], the monomial.
    For each 1 <= I <= M, we have 0 <= X(I).
    Output, int MONO_RANK_GRLEX, the rank.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	int * x = s_data->a1;
	
    dim_typ i;
    dim_typ j;
    dim_typ ks;
    dim_typ n;
    dim_typ nm;
    dim_typ ns;
    dim_typ rank;
    dim_typ tim1;
    int *xs;
    /*
    Ensure that 1 <= M.
    */
    if ( m < 1 )
    {
    	result = USHRT_MAX;
        return &result;
    }
    /*
    Ensure that 0 <= X(I).
    */
    for ( i = 0; i < m; ++i )
        if ( x[i] < 0 )
        {
        	result = USHRT_MAX;
            return &result;
    	}
    /*
    NM = sum ( X )
    */
    nm = i4vec_sum ( m, x );
    /*
    Convert to KSUBSET format.
    */
    ns = nm + m - 1;
    ks = m - 1;
    xs = ( int * ) malloc ( ks * sizeof ( int ) );
    xs[0] = x[0] + 1;
    for ( i = 2; i < m; ++i )
        xs[i-1] += x[i-1] + 1;
    /*
    Compute the rank.
    */
    rank = 1;

    for ( i = 1; i <= ks; +i )
    {
        tim1 = 0 +(i!=1)*xs[i-2];

        if ( tim1 + 1 <= xs[i-1] - 1 )
            for ( j = tim1 + 1; j <= xs[i-1] - 1; ++j )
                rank += i4_choose ( ns - j, ks - i );

    }

    for ( n = 0; n < nm; ++n )
        rank += i4_choose ( n + m - 1, n );

    free ( xs );
    
    result = rank;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_total_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_TOTAL_ENUM enumerates monomials in M dimensions of degree equal to N.
  Discussion:
    For M = 3, we have the following values:
    N  VALUE
    0    1
    1    3
    2    6
    3   10
    4   15
    5   21
    In particular, VALUE(3,3) = 10 because we have the 10 monomials:
      x^3, x^2y, x^2z, xy^2, xyz, xz^3, y^3, y^2z, yz^2, z^3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the maximum degree.
    Output, int MONO_TOTAL_ENUM, the number of monomials in M variables,
    of total degree N.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
	result = i4_choose ( n + m - 1, n );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_total_next_grevlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_TOTAL_NEXT_GREVLEX: grevlex next monomial with total degree equal to N.
  Discussion:
    We consider all monomials in an M dimensional space, with total degree N.
    For example:
    M = 3
    N = 3
    #  X(1)  X(2)  X(3)  Degree
      +------------------------
    1 |  0     0     3        3
    2 |  0     1     2        3
    3 |  1     0     2        3
    4 |  0     2     1        3
    5 |  1     1     1        3
    6 |  2     0     1        3
    7 |  0     3     0        3
    8 |  1     2     0        3
    9 |  2     1     0        3
   10 |  3     0     0        3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the degree.
    0 <= N.
    Input, int X[M], the current monomial.
    To start the sequence, set X = [ 0, 0, ..., 0, N ].
    Output, int X[M], the next monomial.
    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * x = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    dim_typ t;

    if ( n <= 0 || i4vec_sum ( m, x ) != n )
        return NULL;

    if ( x[0] == n )
    {
        x[0] = 0;
        x[m-1] = n;
    }
    else
        mono_next_grevlex ( m, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_total_next_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_TOTAL_NEXT_GRLEX: grlex next monomial with total degree equal to N.
  Discussion:
    We consider all monomials in an M-dimensional space, with total degree N.
    For example:
    M = 3
    N = 3
    #  X(1)  X(2)  X(3)  Degree
      +------------------------
    1 |  0     0     3        3
    2 |  0     1     2        3
    3 |  0     2     1        3
    4 |  0     3     0        3
    5 |  1     0     2        3
    6 |  1     1     1        3
    7 |  1     2     0        3
    8 |  2     0     1        3
    9 |  2     1     0        3
   10 |  3     0     0        3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the degree.
    0 <= N.
    Input/output, int X[M], the current monomial.
    To start the sequence, set X = [ 0, 0, ..., 0, N ].
    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * x = s_data->a2;
	
    dim_typ i;
    dim_typ im1;
    dim_typ j;
    dim_typ t;

    if ( n <= 0 || i4vec_sum ( m, x ) != n )
        return NULL;

    if ( x[0] == n )
    {
        x[0] = 0;
        x[m-1] = n;
    }
    else
        mono_next_grlex ( m, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _mono_total_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_TOTAL_RANDOM: random monomial with total degree equal to N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the degree.
    0 <= N.
    Input/output, int *SEED, the random number seed.
    Output, int *RANK, the rank of the monomial.
    Output, int MONO_TOTAL_RANDOM[M], the random monomial.
*/
{
	const _2dtpipdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	dim_typ * rank = s_data->a3;
	
    dim_typ rank_max;
    dim_typ rank_min;
    rank_min = mono_upto_enum ( m, n - 1 ) + 1;
    rank_max = mono_upto_enum ( m, n );
    *rank = i4_uniform_ab ( rank_min, rank_max, seed );
    return mono_unrank_grlex ( m, *rank );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _mono_unrank_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_UNRANK_GRLEX computes the monomial of given grlex rank.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 January 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    1 <= M.
    Input, int RANK, the rank of the monomial.
    1 <= RANK.
    Output, int MONO_UNRANK_GRLEX[M], the monomial X of the given rank.
    For each I, 0 <= X[I] <= NM, and
    sum ( 1 <= I <= M ) X[I] = NM.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ rank = a_data[1];
	
    dim_typ i;
    dim_typ j;
    dim_typ ks;
    dim_typ nksub;
    dim_typ nm;
    dim_typ ns;
    dim_typ r;
    dim_typ rank1;
    dim_typ rank2;
    int *x;
    int *xs;
    /*
    Ensure that 1 <= M.
    */
    if ( m < 1 || rank < 1)
        return NULL;
    /*
    Special case M == 1.
    */
    if ( m == 1 )
    {
        x = ( int * ) malloc ( m * sizeof ( int ) );
        x[0] = rank - 1;
        return x;
    }
    /*
    Determine the appropriate value of NM.
    Do this by adding up the number of compositions of sum 0, 1, 2,
    ..., without exceeding RANK.  Moreover, RANK - this sum essentially
    gives you the rank of the composition within the set of compositions
    of sum NM.  And that's the number you need in order to do the
    unranking.
    */
    rank1 = 1;
    nm = -1;
    for ( ; ; )
    {
        ++ nm;
        r = i4_choose ( nm + m - 1, nm );
        if ( rank < rank1 + r )
            break;
        rank1 += r;
    }

    rank2 = rank - rank1;
    /*
    Convert to KSUBSET format.
    Apology: an unranking algorithm was available for KSUBSETS,
    but not immediately for compositions.  One day we will come back
    and simplify all this.
    */
    ks = m - 1;
    ns = nm + m - 1;
    xs = ( int * ) malloc ( ks * sizeof ( int ) );

    nksub = i4_choose ( ns, ks );
    j = 1;

    for ( i = 1; i <= ks; ++i )
    {
        r = i4_choose ( ns - j, ks - i );

        while ( r <= rank2 && 0 < r )
        {
            rank2 -= r;
            ++ j;
            r = i4_choose ( ns - j, ks - i );
        }
        xs[i-1] = j;
        j = j + 1;
    }
    /*
    Convert from KSUBSET format to COMP format.
    */
    x = ( int * ) malloc ( m * sizeof ( int ) );
    x[0] = xs[0] - 1;
    for ( i = 2; i < m; ++i)
        x[i-1] -= xs[i-2] - 1;
    x[m-1] = ns - xs[ks-1];

    free ( xs );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_upto_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_UPTO_ENUM enumerates monomials in M dimensions of degree up to N.
  Discussion:
    For M = 2, we have the following values:
    N  VALUE
    0    1
    1    3
    2    6
    3   10
    4   15
    5   21
    In particular, VALUE(2,3) = 10 because we have the 10 monomials:
      1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the maximum degree.
    Output, int MONO_UPTO_ENUM, the number of monomials in
    M variables, of total degree N or less.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
	result = i4_choose ( n + m, n );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_upto_next_grevlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_UPTO_NEXT_GREVLEX: grevlex next monomial with total degree up to N.
  Discussion:
    We consider all monomials in an M-dimensional space, with total
    degree up to N.
    For example:
    M = 3
    N = 3
    #  X(1)  X(2)  X(3)  Degree
      +------------------------
    1 |  0     0     0        0
      |
    2 |  0     0     1        1
    3 |  0     1     0        1
    4 |  1     0     0        1
      |
    5 |  0     0     2        2
    6 |  0     1     1        2
    7 |  1     0     1        2
    8 |  0     2     0        2
    9 |  1     1     0        2
   10 |  2     0     0        2
      |
   11 |  0     0     3        3
   12 |  0     1     2        3
   13 |  1     0     2        3
   14 |  0     2     1        3
   15 |  1     1     1        3
   16 |  2     0     1        3
   17 |  0     3     0        3
   18 |  1     2     0        3
   19 |  2     1     0        3
   20 |  3     0     0        3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the maximum degree.
    0 <= N.
    Input, int X[M], the current monomial.
    To start the sequence, set X = [ 0, 0, ..., 0, 0 ].
    Output, int X[M], the next monomial.
    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * x = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    dim_typ t;

    if ( n <= 0 || i4vec_sum ( m, x ) < 0 || n < i4vec_sum ( m, x ) )
        return NULL;

    if ( x[0] == n )
        x[0] = 0;
    else
        mono_next_grevlex ( m, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mono_upto_next_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_UPTO_NEXT_GRLEX: grlex next monomial with total degree up to N.
  Discussion:
    We consider all monomials in an M-dimensional space, with total
    degree up to N.
    For example:
    M = 3
    N = 3
    #  X(1)  X(2)  X(3)  Degree
      +------------------------
    1 |  0     0     0        0
      |
    2 |  0     0     1        1
    3 |  0     1     0        1
    4 |  1     0     0        1
      |
    5 |  0     0     2        2
    6 |  0     1     1        2
    7 |  0     2     0        2
    8 |  1     0     1        2
    9 |  1     1     0        2
   10 |  2     0     0        2
      |
   11 |  0     0     3        3
   12 |  0     1     2        3
   13 |  0     2     1        3
   14 |  0     3     0        3
   15 |  1     0     2        3
   16 |  1     1     1        3
   17 |  1     2     0        3
   18 |  2     0     1        3
   19 |  2     1     0        3
   20 |  3     0     0        3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the maximum degree.
    0 <= N.
    Input/output, int X[M], the current monomial.
    To start the sequence, set X = [ 0, 0, ..., 0, 0 ].
    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * x = s_data->a2;
	
    dim_typ i;
    dim_typ im1;
    dim_typ j;
    dim_typ t;

    if ( n <= 0 || i4vec_sum ( m, x ) < 0 || n < i4vec_sum ( m, x ))
        return NULL;

    if ( x[0] == n )
        x[0] = 0;
    else
    mono_next_grlex ( m, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _mono_upto_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_UPTO_RANDOM: random monomial with total degree less than or equal to N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the degree.
    0 <= N.
    Input/output, int *SEED, the random number seed.
    Output, int *RANK, the rank of the monomial.
    Output, int MONO_UPTO_RANDOM[M], the random monomial.
*/
{
	const _2dtpipdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	dim_typ * rank = s_data->a3;
	
    dim_typ rank_max;
    dim_typ rank_min;
    rank_min = 1;
    rank_max = mono_upto_enum ( m, n );
    *rank = i4_uniform_ab ( rank_min, rank_max, seed );
    return mono_unrank_grlex ( m, *rank );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _mono_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONO_VALUE evaluates a monomial.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of evaluation points.
    Input, int F[M], the exponents of the monomial.
    Input, double X[M*N], the coordinates of the evaluation points.
    Output, double MONO_VALUE[N], the value of the monomial at X.
*/
{
	const _2dtpipit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * f = s_data->a2;
	ityp * x = s_data->a3;
	
    dim_typ i, j;
    ityp *v = ( double * ) malloc ( n * sizeof ( double ) );

    for ( j = 0; j < n; ++j )
    {
        v[j] = 1.00;
        for ( i = 0; i < m; ++i )
            v[j] *= pow ( x[i+j*m], f[i] );
    }

    return v;
}

#endif
