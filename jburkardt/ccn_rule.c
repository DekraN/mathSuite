#ifndef __DISABLEDEEP_CCNRULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ccn_compute_points_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    CCN_COMPUTE_POINTS: compute Clenshaw Curtis Nested points.
  Discussion:
    We want to compute the following sequence:
    1/2,
    0, 1
    1/4, 3/4
    1/8, 3/8, 5/8, 7/8,
    1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.
    But we would prefer that the numbers in each row be regrouped in pairs
    that are symmetric about 1/2, with the number above 1/2 coming first.
    Thus, the last row might become:
 (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).
    Once we have our sequence, we apply the Chebyshev transformation
    which maps [0,1] to [-1,+1].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 March 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements to compute.
    Output, double CCN_COMPUTE_POINTS_NEW[N], the elements of the sequence.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
	dim_typ d, i, k, m;
	dim_typ td, tu;
	/*
	Handle first three entries specially.
	*/
	if ( 1 <= n )
		x[0] = 0.50;

	if ( 2 <= n )
		x[1] = 1.00;

	if ( 3 <= n )
		x[2] = 0.00;

	m = 3;
	d = 2;

	while ( m < n )
	{
		tu = d + 1;
		td = d - 1;

		k = MIN ( d, n - m );

		for ( i = 1; i <= k; i++ )
		{
			if ( ( i % 2 ) == 1 )
			{
				x[m+i-1] = tu / 2.00 / ( ityp ) ( k );
				tu = tu + 2;
			}
			else
			{
				x[m+i-1] = td / 2.00 / ( ityp ) ( k );
				td -= 2;
			}
		}
		m += k;
		d *= 2;
	}
	/*
	Apply the Chebyshev transformation.
	*/
	for ( i = 0; i < n; ++i )
		x[i] = cos ( x[i] * M_PI );
	x[0] = 0.00;

	if ( 2 <= n )
		x[1] = -1.00;

	if ( 3 <= n )
		x[2] = +1.00;

	return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _rescale ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESCALE rescales a quadrature rule from [-1,+1] to [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 October 2012
  Author:
    John Burkardt.
  Parameters:
    Input, double A, B, the endpoints of the new interval.
    Input, int N, the order.
    Input/output, double X[N], on input, the abscissas for [-1,+1].
    On output, the abscissas for [A,B].
    Input/output, double W[N], on input, the weights for [-1,+1].
    On output, the weights for [A,B].
*/
{
	const dt2it2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	ityp * x = s_data->a3;
	ityp * w = s_data->a4;
	
    for (dim_typ i = 0; i < n; ++i)
    {
        x[i] = ( ( a + b ) + ( b - a ) * x[i] ) / 2.00;
        w[i] = ( b - a ) * w[i] / 2.00;
    }
    return NULL;
}

#endif
