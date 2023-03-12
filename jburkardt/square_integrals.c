#ifndef __DISABLEDEEP_SQUAREINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _square01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    SQUARE01_SAMPLE samples points in the unit square in 2D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 January 2014
  Author:
    John Burkardt
  Reference:
    Russell Cheng,
    Random Variate Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998, pages 168.
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.
  Parameters:
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double X[2*N], the points.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    return r8mat_uniform_01_new ( 2, n, seed );
}

#endif
