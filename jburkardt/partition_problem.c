#ifndef __DISABLEDEEP_PARTITIONPROBLEM

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _partition_brute ( void * data)
/******************************************************************************/
/*
  Purpose:
    PARTITION_BRUTE approaches the partition problem using brute force.
  Discussion:
    We are given a set of N integers W.
    We seek to partition W into subsets W0 and W1, such that the subsets
    have equal sums.
    The "discrepancy" is the absolute value of the difference between the
    two sums, and will be zero if we have solved the problem.
    For a given set of integers, there may be zero, one, or many solutions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the set.
    Input, int W[N], the integers.
    Output, int C[N], indicates the proposed solution.
    C(I) is 0 for items in set W0 and 1 for items in set W1.
    Output, int *DISCREPANCY, the discrepancy.
*/
{
	const dt3pi * const s_data = data; 
	const register dim_typ n = s_data->a0;
	int * w = s_data->a1;
	int * c = s_data->a2;
	int * discrepancy = s_data->a3;
	
    int *d;
    int d_discrepancy;
    int rank;
    int w_sum;

    w_sum = i4vec_sum ( n, w );
    *discrepancy = w_sum;

    rank = -1;
    d = ( int * ) malloc ( n * sizeof ( int ) );

    while ( 1 )
    {
        partition_subset_next ( n, d, &rank );

        if ( rank == -1 )
            break;

        d_discrepancy = i4_abs ( w_sum - (i4vec_dot_product ( n, d, w ) << 1) );

        if ( d_discrepancy < *discrepancy )
        {
            *discrepancy = d_discrepancy;
            i4vec_copy ( n, d, c );
        }

        if ( *discrepancy == 0 )
            break;
    }
    free ( d );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _partition_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    PARTITION_COUNT counts the solutions to a partition problem.
  Discussion:
    We are given a set of N integers W.
    We seek to partition W into subsets W0 and W1, such that the subsets
    have equal sums.
    The "discrepancy" is the absolute value of the difference between the
    two sums, and will be zero if we have solved the problem.
    For a given set of integers, there may be zero, one, or many solutions.
    In the case where the weights are distinct, the count returned by this
    function may be regarded as twice as big as it should be, since the
    partition (W0,W1) is counted a second time as (W1,W0).  A more serious
    overcount can occur if the set W contains duplicate elements - in the
    extreme case, W might be entirely 1's, in which case there is really
    only one (interesting) solution, but this function will count many.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the set.
    Input, int W[N], the integers.
    Output, int PARTITION_COUNT, the number of solutions.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * w = s_data->a1;
	
    int *c;
    dim_typ count;
    int discrepancy;
    int rank;
    int w_sum;

    w_sum = i4vec_sum ( n, w );

    c = ( int * ) malloc ( n * sizeof ( int ) );
    rank = -1;
    count = 0;

    while ( true )
    {
        partition_subset_next ( n, c, &rank );

        if ( rank == -1 )
            break;

        discrepancy = i4_abs ( w_sum - (i4vec_dot_product ( n, c, w ) << 1) );

        if ( discrepancy == 0 )
            ++ count;
    }

    free ( c );

	result = count;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _partition_subset_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_NEXT computes the subset lexicographic successor.
  Discussion:
    This is a lightly modified version of "subset_lex_successor()" from COMBO.
  Example:
    On initial call, N is 5 and the input value of RANK is -1.
    Then here are the successive outputs from the program:
   Rank   T1   T2   T3   T4   T5
   ----   --   --   --   --   --
      0    0    0    0    0    0
      1    0    0    0    0    1
      2    0    0    0    1    0
      3    0    0    0    1    1
     ..   ..   ..   ..   ..   ..
     30    1    1    1    1    0
     31    1    1    1    1    1
     -1    0    0    0    0    0  <-- Reached end of cycle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2012
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.
  Parameters:
    Input, int N, the number of elements in the master set.
    N must be positive.
    Input/output, int T[N], describes a subset.  T(I) is 0 if
    the I-th element of the master set is not in the subset, and is
    1 if the I-th element is part of the subset.
    On input, T describes a subset.
    On output, T describes the next subset in the ordering.
    Input/output, int *RANK, the rank.
    If RANK = -1 on input, then the routine understands that this is
    the first call, and that the user wishes the routine to supply
    the first element in the ordering, which has RANK = 0.
    In general, the input value of RANK is increased by 1 for output,
    unless the very last element of the ordering was input, in which
    case the output value of RANK is -1.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * t = s_data->a1;
	int * rank = s_data->a2;
	
    dim_typ i;
    /*
    Return the first element.
    */
    if ( *rank == -1 )
    {
        for ( i = 0; i < n; ++i )
            t[i] = 0;
        *rank = 0;
        return NULL;
    }

    for ( i = n - 1; 0 <= i; --i )
        if ( t[i] == 0 )
        {
            t[i] = 1;
            ++ *rank;
            return NULL;
        }
        else
            t[i] = 0;
    *rank = -1;

    return NULL;
}

#endif
