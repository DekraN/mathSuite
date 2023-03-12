#ifndef __DISABLEDEEP_LINENCORULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_nco_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_NCO_RULE computes a Newton-Cotes Open (NCO) quadrature rule.
  Discussion:
    The integral:
      Integral ( A <= X <= B ) F(X) dx
    The quadrature rule:
      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order.
    Input, double A, B, the endpoints of the interval.
    Input, double X[N], the abscissas.
    Output, double W[N], the weights.
*/
{
	const dt2it2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp a = s_data->a1;
	ityp b = s_data->a2;
	ityp * x = s_data->a3;
	ityp * w = s_data->a4;
	
    ityp *d;
    dim_typ i, j, k;
    ityp y_a;
    ityp y_b;
    /*
    Define the points X.
    */
    r8vec_linspace2 ( n, a, b, x );

    d = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        /*
        Compute the Lagrange basis polynomial which is 1 at XTAB(I),
        and zero at the other nodes.
        */
        for ( j = 0; j < n; ++j)
            d[j] = 0.00;
        d[i] = 1.0;

        for ( j = 2; j <= n; ++j )
            for ( k = j; k <= n; ++k )
                d[n+j-k-1] = ( d[n+j-k-2] - d[n+j-k-1] ) / ( x[n-k] - x[n+j-k-1] );

        for ( j = 1; j <= n - 1; ++j )
            for ( k = 1; k <= n - j; ++k )
                d[n-k-1] -= x[n-k-j] * d[n-k];
        /*
        Evaluate the antiderivative of the polynomial at the endpoints.
        */
        y_a = d[n-1] / ( ityp ) ( n );
        for ( j = n - 2; 0 <= j; --j )
            y_a *= a + d[j] / ( ityp ) ( j + 1 );
        y_a *= a;

        y_b = d[n-1] / ( ityp ) ( n );
        for ( j = n - 2; 0 <= j; --j )
            y_b *= b + d[j] / ( ityp ) ( j + 1 );
        y_b *= b;

        w[i] = y_b - y_a;
    }

    free ( d );

    return NULL;
}

#endif
