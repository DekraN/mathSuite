#ifndef __DISABLEDEEP_KNAPSACK01

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subset_gray_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_GRAY_NEXT generates all subsets of a set of order N, one at a time.
  Discussion:
    It generates the subsets one at a time, by adding or subtracting
    exactly one element on each step.
    The user should set MORE = .FALSE. and the value of N before
    the first call.  On return, the user may examine A which contains
    the definition of the new subset, and must check .MORE., because
    as soon as it is .FALSE. on return, all the subsets have been
    generated and the user probably should cease calling.
    The first set returned is the empty set.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the order of the total set from which
    subsets will be drawn.
    Input/output, int A[N].  On each return, the Gray code for the newly
    generated subset.  A[I] = 0 if element I is in the subset, 1 otherwise.
    Input/output, int *MORE.  Set this variable FALSE before
    the first call.  Normally, MORE will be returned TRUE but once
    all the subsets have been generated, MORE will be
    reset FALSE on return and you should stop calling the program.
    Input/output, int *NCARD, the cardinality of the set returned,
    which may be any value between 0 (the empty set) and N (the
    whole set).
    Output, int *IADD, the element which was added or removed to the
    previous subset to generate the current one.  Exception:
    the empty set is returned on the first call, and IADD is set to -1.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * more = s_data->a2;
	int * ncard = s_data->a3;
	int * iadd = s_data->a4;
	
    int i;
    /*
    First set returned is the empty set.
    */
    if ( !(*more) )
    {
        for ( i = 0; i < n; ++i )
            a[i] = 0;

        *iadd = *ncard = 0;
        *more = -1;
    }
    else
    {
        *iadd = 1;

        if ( ( *ncard % 2 ) != 0 )
        {
            for ( ; ; )
            {
                ++ *iadd;
                if ( a[*iadd-2] != 0 )
                    break;
            }
        }

        a[*iadd-1] = 1 - a[*iadd-1];
        *ncard += (a[*iadd-1]<<1) - 1;
        /*
        Last set returned is the singleton A(N).
        */
        if ( *ncard == a[n-1] )
            *more = 0;
    }
    return NULL;
}

#endif
