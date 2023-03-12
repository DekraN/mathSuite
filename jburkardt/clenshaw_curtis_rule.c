#ifndef __DISABLEDEEP_CLENSHAWCURTISRULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _clenshaw_curtis_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
  Discussion:
    The integration interval is [ -1, 1 ].
    The weight function is w(x) = 1.0.
    The integral to approximate:
      Integral ( -1 <= X <= 1 ) F(X) dX
    The quadrature rule:
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    1 <= ORDER.
    Output, double X[ORDER], the abscissas.
    Output, double W[ORDER], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0;
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
    ityp b;
    dim_typ i, j;
    ityp theta;

    if ( order < 1 )
        return NULL;
    else if ( order == 1 )
    {
        x[0] = 0.00;
        w[0] = 2.00;
    }
    else
    {
        for ( i = 0; i < order; ++i)
            x[i] =  cos ( ( ityp ) ( order - 1 - i ) * M_PI/ ( ityp ) ( order - 1     ) );
        x[0] = -1.00;
        if ( ( order % 2 ) == 1 )
            x[(order-1)/2] = 0.00;
        x[order-1] = +1.00;

        for ( i = 0; i < order; ++i )
        {
            theta = ( ityp ) ( i ) * M_PI / ( ityp ) ( order - 1 );
            w[i] = 1.00;

            for ( j = 1; j <= ( order - 1 ) / 2; ++j )
            {
                b = (j<<1) == ( order - 1 ) ? 1.00:2.00;

                w[i] -= b * cos ( 2.00 * ( ityp ) ( j ) * theta )/ ( ityp ) ( (j<<2) * j - 1 );
            }
        }

        w[0] /= ( ityp ) ( order - 1 );
        for ( i = 1; i < order - 1; ++i )
            w[i] = 2.0 * w[i] / ( ityp ) ( order - 1 );

        w[order-1] /= ( ityp ) ( order - 1 );
    }

    return NULL;
}

#endif
