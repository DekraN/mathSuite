#ifndef __DISABLEDEEP_LLSQ

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _llsq ( void * data)
/******************************************************************************/
/*
  Purpose:
    LLSQ solves a linear least squares problem matching a line to data.
  Discussion:
    A formula for a line of the form Y = A * X + B is sought, which
    will minimize the root-mean-square error to N data points ( X[I], Y[I] );
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 March 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, double X[N], Y[N], the coordinates of the data points.
    Output, double *A, *B, the slope and Y-intercept of the least-squares
    approximant to the data.
*/
{
	const dt4pit * const s_data = data;
	const register dim_typ n = s_data->a0; 
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	
    ityp bot;
    dim_typ i;
    ityp top;
    ityp xbar;
    ityp ybar;
    /*
    Special case.
    */
    if ( n == 1 )
    {
        *a = 0.00;
        *b = y[0];
        return NULL;
    }
    /*
    Average X and Y.
    */
    ybar = ybar = 0.00;
    for ( i = 0; i < n; ++i )
    {
        xbar += x[i];
        ybar += y[i];
    }
    xbar /= ( ityp ) n;
    ybar /= ( ityp ) n;
    /*
    Compute Beta.
    */
    top = bot = 0.00;
    for ( i = 0; i < n; ++i )
    {
        top += ( x[i] - xbar ) * ( y[i] - ybar );
        bot += ( x[i] - xbar ) * ( x[i] - xbar );
    }
    *a = top / bot;
    *b = ybar - *a * xbar;
    return NULL;
}

#endif
