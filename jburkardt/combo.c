#ifndef __DISABLEDEEP_COMBO

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _bell_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BELL_VALUES returns some values of the Bell numbers.
  Discussion:
    The Bell number B(N) is the number of restricted growth functions on N.
    Note that the Stirling numbers of the second kind, S^m_n, count the
    number of partitions of N objects into M classes, and so it is
    true that
      B(N) = S^1_N + S^2_N + ... + S^N_N.
    The Bell numbers were named for Eric Temple Bell.
    In Mathematica, the function can be evaluated by
      Sum[StirlingS2[n,m],{m,1,n}]
  Definition:
    The Bell number B(N) is defined as the number of partitions (of
    any size) of a set of N distinguishable objects.
    A partition of a set is a division of the objects of the set into
    subsets.
  Examples:
    There are 15 partitions of a set of 4 objects:
   (1234),
   (123) (4),
   (124) (3),
   (12) (34),
   (12) (3) (4),
   (134) (2),
   (13) (24),
   (13) (2) (4),
   (14) (23),
   (1) (234),
   (1) (23) (4),
   (14) (2) (3),
   (1) (24) (3),
   (1) (2) (34),
   (1) (2) (3) (4).
    and so B(4) = 15
  First values:
     N         B(N)
     0           1
     1           1
     2           2
     3           5
     4          15
     5          52
     6         203
     7         877
     8        4140
     9       21147
    10      115975
  Recursion:
    B(I) = sum ( 1 <= J <=I ) Binomial ( I-1, J-1 ) * B(I-J)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 February 2003
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the Bell number.
    Output, int *C, the value of the Bell number.
*/
{
	dim_typ ** const a_data = data;
	dim_typ * n_data = a_data[0];
	dim_typ * n = a_data[1];
	dim_typ * c = a_data[2];
	
    # define N_MAX 11

    static unsigned c_vec[N_MAX] =
    {
        1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975
    };

    static unsigned n_vec[N_MAX] =
    {
        0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_search_binary_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
  Discussion:
    An I4VEC is a vector of I4's.
    Binary search is used.
    Note that zero-based indexing is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 October 2015
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Algorithm 1.9,
    Combinatorial Algorithms,
    CRC Press, 1998, page 26.
  Parameters:
    Input, int N, the number of elements in the vector.
    Input, int A[N], the array to be searched.  A must
    be sorted in ascending order.
    Input, int B, the value to be searched for.
    Output, int I4VEC_SEARCH_BINARY_A, the result of the search.
    -1, B does not occur in A.
    I, A[I] = B.
*/
{
	static int result = INT_MAX;
	
	const dtpii * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int b = s_data->a2;
	
    int high;
    int index;
    int low;
    int mid;
    /*
    Check.
    */
    if ( n <= 0 )
    {
    	result = -2;
        return &result;
    }

    index = -1;
    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ( low + high ) / 2;

        if ( a[mid] == b )
        {
            index = mid;
            break;
        }
        else if ( a[mid] < b )
            low = mid + 1;
        else if ( b < a[mid] )
            high = mid - 1;
    }
    
    result = index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *    _i4vec_search_binary_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
  Discussion:
    An I4VEC is a vector of I4's.
    Binary search is used.
    Note that zero-based indexing is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 October 2015
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Algorithm 1.9,
    Combinatorial Algorithms,
    CRC Press, 1998, page 26.
  Parameters:
    Input, int N, the number of elements in the vector.
    Input, int A[N], the array to be searched.  A must
    be sorted in descending order.
    Input, int B, the value to be searched for.
    Output, int I4VEC_SEARCH_BINARY_D, the result of the search.
    -1, B does not occur in A.
    I, A[I] = B.
*/
{
	static int result = INT_MAX;
	
	const dtpii * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int b = s_data->a2;
	
    int high;
    int index;
    int low;
    int mid;
    /*
    Check.
    */
    if ( n <= 0 )
    {
    	result = -2;
        return &result;
    }

    index = -1;
    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ( low + high ) / 2;

        if ( a[mid] == b )
        {
            index = mid;
            break;
        }
        else if ( b < a[mid] )
            low = mid + 1;
        else if ( a[mid] < b )
            high = mid - 1;
    }
    
    result = index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sort_insert_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Algorithm 1.1,
    Combinatorial Algorithms,
    CRC Press, 1998, page 11.
  Parameters:
    Input, int N, the number of items in the vector.
    N must be positive.
    Input/output, int A[N].
    On input, A contains data to be sorted.
    On output, the entries of A have been sorted in ascending order.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;
    int x;

    for ( i = 1; i < n; ++i )
    {
        x = a[i];

        j = i;

        while ( 1 <= j && x < a[j-1] )
        {
            a[j] = a[j-1];
            -- j;
        }

        a[j] = x;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4vec_sort_insert_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Algorithm 1.1,
    Combinatorial Algorithms,
    CRC Press, 1998, page 11.
  Parameters:
    Input, int N, the number of items in the vector.
    N must be positive.
    Input/output, int A[N].
    On input, A contains data to be sorted.
    On output, the entries of A have been sorted in ascending order.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;
    int x;

    for ( i = 1; i < n; ++i )
    {
        x = a[i];
        j = i;

        while ( 1 <= j && a[j-1] < x )
        {
            a[j] = a[j-1];
            -- j;
        }
        a[j] = x;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_backtrack ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_BACKTRACK supervises a backtrack search for an I4VEC.
  Discussion:
    The routine tries to construct an integer vector one index at a time,
    using possible candidates as supplied by the user.
    At any time, the partially constructed vector may be discovered to be
    unsatisfactory, but the routine records information about where the
    last arbitrary choice was made, so that the search can be
    carried out efficiently, rather than starting out all over again.
    First, call the routine with INDX = 0 so it can initialize itself.
    Now, on each return from the routine, if INDX is:
      1, you've just been handed a complete candidate vector;
         Admire it, analyze it, do what you like.
      2, please determine suitable candidates for position X(K).
         Return the number of candidates in NCAN(K), adding each
         candidate to the end of STACK, and increasing NSTACK.
      3, you're done.  Stop calling the routine;
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 July 2004
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C++ version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the number of positions to be filled in the vector.
    Input, int MAXSTACK, the maximum length of the stack.
    Input, int STACK[MAXSTACK], a list of all current candidates for
    all positions 1 through K.
    Input/output, int X[N], the partial or complete candidate vector.
    Input/output, int *INDX, a communication flag.
    On input,
      0 to start a search.
    On output:
      1, a complete output vector has been determined and returned in X(1:N);
      2, candidates are needed for position X(K);
      3, no more possible vectors exist.
    Input/output, int *K, if INDX=2, the current vector index being considered.
    Input/output, int *NSTACK, the current length of the stack.
    Input/output, int NCAN[N], lists the current number of candidates for
    positions 1 through K.
*/
{
	const _2dt6pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ maxstack = s_data->a1;
	int * stack = s_data->a2;
	int * x = s_data->a3;
	int * indx = s_data->a4;
	int * k = s_data->a5;
	int * nstack = s_data->a6;
	int * ncan = s_data->a7;
	
	/*
	If this is the first call, request a candidate for position 1.
	*/
	if ( *indx == 0 )
	{
		*k = 1;
		*nstack = 0;
		*indx = 2;
		return NULL;
	}
	/*
	Examine the stack.
	*/
	for ( ; ; )
	{
		/*
		If there are candidates for position K, take the first available
		one off the stack, and increment K.
		
		This may cause K to reach the desired value of N, in which case
		we need to signal the user that a complete set of candidates
		is being returned.
		*/
		if ( 0 < ncan[(*k)-1] )
		{
			x[(*k)-1] = stack[(*nstack)-1];
			-- *nstack;
			
			ncan[(*k)-1] = ncan[(*k)-1] - 1;
			
			if ( *k != n )
			{
				++ *k;
				*indx = 2;
			}
			else
				*indx = 1;
			break;
		}
		/*
		If there are no candidates for position K, then decrement K.
		If K is still positive, repeat the examination of the stack.
		*/
		else
		{
			-- *k;
	
			if ( *k <= 0 )
			{
				*indx = 3;
				break;
			}
		}
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sum ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SUM sums the entries of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    Input:
      A = ( 1, 2, 3, 4 )
    Output:
      I4VEC_SUM = 10
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector to be summed.
    Output, int I4VEC_SUM, the sum of the entries of A.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int sum = 0;
    for (dim_typ i = 0; i < n; ++i )
        sum += a[i];
        
    result = sum;
    return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _knapsack_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    KNAPSACK_01 seeks a solution of the 0/1 Knapsack problem.
  Discussion:
    In the 0/1 knapsack problem, a knapsack of capacity C is given,
    as well as N items, with the I-th item of weight W(I).
    A selection is "acceptable" if the total weight is no greater than C.
    It is desired to find an optimal acceptable selection, that is,
    an acceptable selection such that there is no acceptable selection
    of greater weight.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of weights.
    Input, inte W[N], the weights.
    Input, int C, the maximum weight.
    Output, int KNAPSACK_01[N], is a binary vector which defines an
    optimal selection.  It is 1 for the weights to be selected, and
    0 otherwise.
*/
{
	const _2dtpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ c = s_data->a1;
	int * w = s_data->a2;
	
    dim_typ i;
    int iadd;
    int more;
    int ncard;
    int *s;
    int *s_test;
    int t;
    int t_test;

    s = ( int * ) malloc ( n * sizeof ( int ) );
    s_test = ( int * ) malloc ( n * sizeof ( int ) );

    more = ncard = 0;

    for ( i = 0; i < n; ++i )
        s_test[i] = s[i] = 0;

    t_test = t = 0;

    for ( ; ; )
    {
        subset_gray_next ( n, s_test, &more, &ncard, &iadd );
        t_test = 0;
        for ( i = 0; i < n; ++i )
            t_test += s_test[i] * w[i];

        if ( t < t_test && t_test <= c )
        {
            t = t_test;
            for ( i = 0; i < n; ++i )
                s[i] = s_test[i];
        }

        if ( ! more )
            break;
    }

    free ( s_test );
    return s;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_CHECK checks that a vector represents a permutation.
  Discussion:
    The routine verifies tat each of the integers from 1
    to N occurs among the N entries of the permutation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int P[N], the array to check.
    Output, int PERM_CHECK, is TRUE if the permutation is OK.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
    dim_typ found, i, seek;

    for ( seek = 1; seek <= n; ++seek )
    {
        found = 0;

        for ( i = 0; i < n; ++i )
            if ( p[i] == seek )
                break;

        if ( !found )
        {
        	result = false;
            return &result;
        }

    }

    result = true;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm_inv ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_INV inverts a permutation "in place".
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects being permuted.
    Input/output, int P[N], the permutation, in standard index form.
    On output, P describes the inverse permutation
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
    dim_typ i;
    dim_typ i0;
    dim_typ i1;
    dim_typ i2;
    dim_typ is;

    if ( n <= 0 || !perm_check ( n, p ))
        return NULL;

    is = 1;

    for ( i = 1; i <= n; ++i )
    {
        i1 = p[i-1];

        while ( i < i1 )
        {
            i2 = p[i1-1];
            p[i1-1] = -i2;
            i1 = i2;
        }

        is = - i4_sign ( p[i-1] );
        p[i-1] = i4_sign ( is ) * abs ( p[i-1] );
    }

    for ( i = 1; i <= n; ++i )
    {
        i1 = -p[i-1];

        if ( 0 <= i1 )
        {
            i0 = i;

            for ( ; ; )
            {
                i2 = p[i1-1];
                p[i1-1] = i0;

                if ( i2 < 0 )
                    break;
                i0 = i1;
                i1 = i2;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _perm_lex_rank ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_LEX_RANK computes the lexicographic rank of a permutation.
  Discussion:
    The original code altered the input permutation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2011
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.
  Parameters:
    Input, int N, the number of values being permuted.
    N must be positive.
    Input, int P[N], describes the permutation.
    P[I] is the item which is permuted into the I-th place
    by the permutation.
    Output, int PERM_LEX_RANK, the rank of the permutation.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
    dim_typ i;
    dim_typ j;
    int *pcopy;
    dim_typ rank;
    /*
    Check.
    */
    perm_check ( n, p );

    rank = 0;
    pcopy = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i)
        pcopy[i] = p[i];

    for ( j = 0; j < n; ++j)
    {
        rank += ( pcopy[j] - 1 ) * i4_factorial ( n - 1 - j );
        for ( i = j + 1; i < n; ++i)
            if ( pcopy[j] < pcopy[i] )
                pcopy[i] *= -1;
    }
    free ( pcopy );

	result = rank;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _perm_lex_unrank ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_LEX_UNRANK computes the permutation of given lexicographic rank
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 June 2012
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.
  Parameters:
    Input, int N, the number of values being permuted.
    N must be positive.
    Input, int RANK, the rank of the permutation
    Output, int PERM_LEX_UNRANK[N], describes the permutation.
*/
{
	dim_typ * const a_data = data; 
	const register dim_typ n = a_data[0];
	const register dim_typ rank = a_data[1];
	
    int d;
    dim_typ i;
    dim_typ j;
    dim_typ nperm;
    int *p;
    dim_typ rank_copy;
    /*
    Check.
    */
    if ( n < 1 )
        return NULL;

    nperm = perm_enum ( n );

    if ( rank < 0 || nperm < rank )
        return NULL;

    rank_copy = rank;

    p = ( int * ) malloc ( n * sizeof ( int ) );

    p[n-1] = 1;

    for ( j = 1; j <= n - 1; ++j )
    {
        d = ( rank_copy % i4_factorial ( j + 1 ) ) / i4_factorial ( j );
        rank_copy -= d * i4_factorial ( j );
        p[n-j-1] = d + 1;

        for ( i = n - j + 1; i <= n; ++i )
            if ( d < p[i-1] )
                ++ p[i-1];
    }
    return p;
}

#endif

