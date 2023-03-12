#ifndef __DISABLEDEEP_COMBINATIONLOCK

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *    _bicycle_lock ( void * data)
/******************************************************************************/
/*
  Purpose:
    BICYCLE_LOCK finds the combination on a typical bicycle lock.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int C, the combination, a value between 0 and 999.
    Output, int BICYCLE_LOCK, the step on which the combination
    was found.  A value of -1 means the combination was not found.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ c = *(dim_typ *) data;
	
    dim_typ a;
    dim_typ step = -1;

    for ( a = 0; a <= 999; ++a )
    {
        if ( a == c )
        {
            step = a + 1;
            break;
        }
    }
    
    result = step;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _combination_lock ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMBINATION_LOCK determines the combination of a lock.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of dials.
    Input, int N, the number of symbols on each dial.
    We assume the symbols are the integers 0 to N-1.
    Input, int C[M], the combination.
    Output, int COMBINATION_LOCK, the step on which the combination
    was found.  A value of -1 means the combination was not found.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * c = s_data->a2;
	
    int a[m];
    dim_typ i;
    bool more;
    dim_typ step;
    /*
    Starting with the guess (0, 0, ... 0),
    generate every possible combination, in order, and try it.
    */
    more = false;
    step = 0;

    while ( true )
    {
        combination_next ( m, n, a, &more );

        if ( !more )
        {
            step = -1;
            break;
        }

        ++ step;

        if ( i4vec_eq ( m, a, c ) )
            break;
    }
    
    result = step;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _combination_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMBINATION_NEXT generates lock combinations in lex order.
  Discussion:
    The vectors are produced in lexical order, starting with
 (0,0,...,0),
 (0,0,...,1),
    ...
 (BASE-1,BASE-1,...,BASE-1).
  Example:
    M = 2,
    BASE = 3
    0   0
    0   1
    0   2
    1   0
    1   1
    1   2
    2   0
    2   1
    2   2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:
    John Burkardt
  Reference:
    Dennis Stanton, Dennis White,
    Constructive Combinatorics,
    Springer, 1986,
    ISBN: 0387963472,
    LC: QA164.S79.
  Parameters:
    Input, int M, the size of the vectors to be used.
    Input, int BASE, the base to be used.  BASE = 2 will
    give vectors of 0's and 1's, for instance.
    Input/output, int A[M].  The input value of A is
    not important on the first call.  Thereafter, it should simply be the
    output value from the previous call.  The output value is the next vector
    in the sequence.

    Inpt/output, bool *MORE.  The input value should be FALSE on the first
    call, and TRUE on subsequent calls.  The output value will be TRUE as long
    as the next vector could be computed, and FALSE once there are no more.
*/
{
	const _2dtpipb * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ base = s_data->a1;
	int * a = s_data->a2;
	bool * more = s_data->a3;
	
    int i;

    if ( !(*more) )
    {
        for ( i = 0; i < m; ++i )
            a[i] = 0;
        *more = true;
    }
    else
    {
        for ( i = m - 1; 0 <= i; --i )
        {
            ++ a[i];

            if ( a[i] < base )
                return NULL;
            a[i] = 0;
        }
        *more = false;
    }
    return NULL;
}

#endif
