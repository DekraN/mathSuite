#ifndef __DISABLEDEEP_UNICYCLE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT void * _i4vec_indicator ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDICATOR sets an I4VEC to the indicator vector.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, int A[N], the initialized array.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
for (dim_typ i = 0; i < n; ++i )
	a[i] = i + 1;
  return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_indicator_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDICATOR_NEW sets an I4VEC to the indicator vector.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, int I4VEC_INDICATOR_NEW[N], the array.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
	int *a = ( int * ) malloc ( n * sizeof ( int ) );
	
	for (dim_typ i = 0; i < n; ++i )
		a[i] = i + 1;
	return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _perm_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_INVERSE computes the inverse of a permutation.
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
    P(I) is the item which is permuted into the I-th place
    by the permutation.
    Output, int PERM_INVERSE[N], the inverse permutation.
*/ 
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
	/*
	Check.
	*/ 
	perm_check ( n, p );
	int * pinv = ( int * ) malloc ( n * sizeof ( int ) );
	
	for (dim_typ i = 0; i < n; ++i )
		pinv[p[i]-1] = i + 1;
	
	return pinv;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm_lex_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_LEX_NEXT computes the lexicographic permutation successor.
  Example:
    RANK  Permutation
       0  1 2 3 4
       1  1 2 4 3
       2  1 3 2 4
       3  1 3 4 2
       4  1 4 2 3
       5  1 4 3 2
       6  2 1 3 4
       ...
      23  4 3 2 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 July 2011
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
    Input/output, int P[N], describes the permutation.
    P(I) is the item which is permuted into the I-th place
    by the permutation.
    Input/output, int *RANK, the rank.
    If RANK = -1 on input, then the routine understands that this is
    the first call, and that the user wishes the routine to supply
    the first element in the ordering, which has RANK = 0.
    In general, the input value of RANK is increased by 1 for output,
    unless the very last element of the ordering was input, in which
    case the output value of RANK is 0.
*/
{
	const dtpipdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	dim_typ * rank = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    int temp;
    /*
    Return the first element.
    */
    if ( *rank == -1 )
    {
        i4vec_indicator ( n, p );
        rank = 0;
        return NULL;
    }
    /*
    Check.
    */
    perm_check ( n, p );
    /*
    Seek I, the highest index for which the next element is bigger.
    */
    i = n - 1;

    for ( ; ; )
    {
        if ( p[i-1] <= p[i] || i <= 0 )
            break;
        -- i;
    }
    /*
    If no I could be found, then we have reach the final permutation,
    N, N-1, ..., 2, 1.  Time to start over again.
    */
    if ( i == 0 )
    {
        i4vec_indicator ( n, p );
        *rank = -1;
    }
    else
    {
        /*
        Otherwise, look for the the highest index after I whose element
        is bigger than I''s.  We know that I+1 is one such value, so the
        loop will never fail.
        */
        j = n;
        while ( p[j-1] < p[i-1] )
            -- j;
        /*
        Interchange elements I and J.
        */
        temp = p[i-1];
        p[i-1] = p[j-1];
        p[j-1] = temp;
        /*
        Reverse the elements from I+1 to N.
        */
        i4vec_reverse ( n - i, p+i );

        ++ *rank;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_ENUM enumerates the permutations on N digits.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 July 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of values being permuted.
    N must be nonnegative.
    Output, int PERM_ENUM, the number of distinct elements.
*/ 
{
	static int result = INT_MAX;
	
	const register int n = *(int *) data;
	
	result = i4_factorial ( n );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _unicycle_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_CHECK checks that a vector represents a unicycle.
  Discussion:
    A unicycle is a permutation with a single cycle.  This might be called
    a cyclic permutation, except that that name is used with at least two
    other meanings.  So, to be clear, a unicycle is a permutation of N
    objects in which each object is returned to itself precisely after
    N applications of the permutation.
    This routine verifies that each of the integers from 1
    to N occurs among the N entries of the permutation.
    Any permutation of the integers 1 to N describes a unicycle.
    The permutation ( 3, 4, 2, 1 ) indicates that the unicycle
    sends 3 to 4, 4 to 2, 2 to 1 and 1 to 3.  This is the sequential
    description of a unicycle.
    The standard sequence "rotates" the permutation so that it begins
    with 1.  The above sequence becomes a standard sequence when
    written as ( 1, 3, 4, 2 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 June 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int U[N], the unicycle sequence vector
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * u = s_data->a1;
	
    bool error;
    dim_typ i;
    dim_typ iseek;

    for ( iseek = 1; iseek <= n; ++iseek )
    {
        error = true;

        for ( i = 0; i < n; ++i )
            if ( u[i] == iseek )
            {
                error = 0;
                break;
            }

        if ( error )
            return NULL;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _unicycle_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_ENUM enumerates the unicycles.
  Discussion:
    Each standard sequence corresponds to a unique unicycle.  Since the
    first element of a standard sequence is always 1, the number of standard
    sequences, and hence the number of unicycles, is (n-1)!.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 June 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the unicyle.
    Output, int UNICYCLE_ENUM, the number of unicycles.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(int *) data;
	
	result = i4_factorial ( n - 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _unicycle_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_INDEX returns the index form of a unicycle.
  Example:
    N = 4
    U       = 1 3 4 2
    U_INDEX = 3 1 4 2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 June 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the unicycles.
    Input, int U[N], the unicycle sequence vector.
    Output, int UNICYCLE_INDEX[N], the unicycle index vector.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * u = s_data->a1;
	
    int ip1;
    dim_typ *u_index = ( dim_typ * ) malloc ( n * sizeof ( dim_typ ) );

    for (dim_typ i = 0; i < n; ++i )
    {
        ip1 = i4_wrap ( i + 1, 0, n - 1 );
        u_index[u[i]-1] = u[ip1];
    }

    return u_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _unicycle_index_to_sequence ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_INDEX_TO_SEQUENCE converts a unicycle from index to sequence form.
  Example:
    N = 4
    U_INDEX = 3 1 4 2
    U       = 1 3 4 2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 June 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the unicycles.
    Output, int U_INDEX(N), the unicycle index vector.
    Input, int U(N), the unicycle sequence vector.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * u_index = s_data->a1;
	
    dim_typ i, j;
    dim_typ *u = ( dim_typ * ) malloc ( n * sizeof ( dim_typ ) );
    u[0] = i = 1;

    for ( j = 1; j < n; ++j )
    {
        i = u_index[i-1];
        u[j] = i;

        if ( i == 1 )
            return NULL;
    }

    return u;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _unicycle_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_INVERSE returns the inverse of a unicycle.
  Example:
    N = 4
    U         = 1 3 4 2
    U_INVERSE = 1 2 4 3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 June 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the unicycles.
    Input, int U[N], the unicycle sequence vector.
    Output, int UNICYCLE_INVERSE[N], the inverse unicycle.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * u = s_data->a1;
	
    dim_typ i;
    int *u_inverse;

    u_inverse = ( int * ) malloc ( n * sizeof ( int ) );
    u_inverse[0] = 1;
    for ( i = 1; i < n; ++i)
        u_inverse[i] = u[n-i];
    return u_inverse;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _unicycle_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_NEXT generates unicycles in lexical order, one at a time.
  Example:
    N = 4
    1   1 2 3 4
    2   1 2 4 3
    3   1 3 2 4
    4   1 3 4 2
    5   1 4 2 3
    6   1 4 3 2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 August 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the unicycles.
    Input/output, int U[N]; on first call with MORE = FALSE,
    this value is not used.  Otherwise, the input value is the previous
    unicycle.  The output value is the next unicycle.
    Input/output, int *RANK, the rank.
    If RANK = -1 on input, then the routine understands that this is
    the first call, and that the user wishes the routine to supply
    the first element in the ordering, which has RANK = 0.
    In general, the input value of RANK is increased by 1 for output,
    unless the very last element of the ordering was input, in which
    case the output value of RANK is -1.
*/
{
	const dtpipdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * u = s_data->a1;
	dim_typ * rank = s_data->a2;
	
    dim_typ i;
    int *p = ( int * ) malloc ( ( n - 1 ) * sizeof ( int ) );

    if ( *rank == -1 )
        u[0] = 1;
    else
        for ( i = 0; i < n - 1; ++i)
            p[i] = u[i+1] - 1;

    perm_lex_next ( n - 1, p, rank );

    for ( i = 0; i < n - 1; ++i )
        u[i+1] = p[i] + 1;
    free ( p );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _unicycle_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_RANDOM selects a random unicycle of N objects.
  Discussion:
    The routine assumes the objects are labeled 1, 2, ... N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 June 2012
  Author:
    John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the number of objects to be permuted.
    Input/output, int *SEED, a seed for the random number
    generator.
    Output, int UNICYCLE_RANDOM[N], a unicycle in sequence form.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    dim_typ i;
    int j;
    int *u;
    int t;

    u = i4vec_indicator_new ( n );

    for ( i = 1; i < n; ++i )
    {
        j = i4_uniform_ab ( i, n - 1, seed );
        t = u[i];
        u[i] = u[j];
        u[j] = t;
    }

    return u;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _unicycle_rank ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_RANK computes the rank of a unicycle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 June 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the order of the unicycle.
    Input, int U[N], a unicycle in sequence form.
    Output, int UNICYLE_RANK, the rank of the unicycle.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * u = s_data->a1;
	
    dim_typ i;
    int *p;
    dim_typ rank;

    p = ( int * ) malloc ( ( n - 1 ) * sizeof ( int ) );

    for ( i = 0; i < n - 1; ++i )
        p[i] = u[i+1] - 1;

    rank = perm_lex_rank ( n - 1, p );

    free ( p );

	result = rank;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _unicycle_unrank ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNICYCLE_UNRANK "unranks" a unicycle.
  Discussion:
    That is, given a rank, it computes the corresponding unicycle.
    The value of the rank should be between 0 and (N-1)!-1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 June 2012
  Author:
    John Burkardt.
  Reference:
    Dennis Stanton, Dennis White,
    Constructive Combinatorics,
    Springer, 1986,
    ISBN: 0387963472,
    LC: QA164.S79.
  Parameters:
    Input, int N, the number of elements in the set.
    Input, int RANK, the desired rank of the permutation.
    Output, int UNICYCLE_UNRANK[N], the unicycle.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ rank = a_data[1];
	
    dim_typ i;
    int *p = perm_lex_unrank ( n - 1, rank );
    int *u = ( int * ) malloc ( n * sizeof ( int ) );

    u[0] = 1;
    for ( i = 0; i < n - 1; ++i )
        u[i+1] = p[i] + 1;
    free ( p );

    return u;
}

#endif
