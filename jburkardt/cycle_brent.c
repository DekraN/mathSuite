#ifndef __DISABLEDEEP_CYCLEBRENT

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _cycle_brent ( void * data)
/******************************************************************************/
/*
  Purpose:
    CYCLE_BRENT finds a cycle in an iterated mapping using Brent's method.
  Discussion:
    Suppose we a repeatedly apply a function f(), starting with the argument
    x0, then f(x0), f(f(x0)) and so on.  Suppose that the range of f is finite.
    Then eventually the iteration must reach a cycle.  Once the cycle is reached,
    succeeding values stay within that cycle.
    Starting at x0, there is a "nearest element" of the cycle, which is
    reached after MU applications of f.
    Once the cycle is entered, the cycle has a length LAM, which is the number
    of steps required to first return to a given value.
    This function uses Brent's method to determine the values of MU and LAM,
    given F and X0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 June 2012
  Author:
    John Burkardt
  Reference:
    Richard Brent,
    An improved Monte Carlo factorization algorithm,
    BIT,
    Volume 20, Number 2, 1980, pages 176-184.
  Parameters:
    Input, int F ( int i ), the name of the function
    to be analyzed.
    Input, int X0, the starting point.
    Output, int *LAM, the length of the cycle.
    Output, int *MU, the index in the sequence starting
    at X0, of the first appearance of an element of the cycle.
*/
{
	const fdtdt2pdt * const s_data = data;
	dim_typ (* f)(dim_typ) = s_data->a0;
	const register dim_typ x0 = s_data->a1; 
	dim_typ * lam = s_data->a2;
	dim_typ * mu = s_data->a3;
	
    dim_typ hare;
    dim_typ i;
    dim_typ power;
    dim_typ tortoise;

    power = 1;
    *lam = 1;
    tortoise = x0;
    hare = f ( x0 );

    while ( tortoise != hare )
    {
        if ( power == *lam )
        {
            tortoise = hare;
            power <<= 1;
            *lam = 0;
        }
        hare = f ( hare );
        ++ *lam;
    }

    *mu = 0;
    tortoise = x0;
    hare = x0;

    for ( i = 0; i < *lam; ++i )
        hare = f ( hare );

    while ( tortoise != hare )
    {
        tortoise = f ( tortoise );
        hare = f ( hare );
        ++ *mu;
    }

    return NULL;
}

#endif
