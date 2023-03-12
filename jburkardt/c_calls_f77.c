#ifndef __DISABLEDEEP_CCALLSF77

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _kronrod ( void * data)
/******************************************************************************/
/*
  Purpose:
    KRONROD adds N+1 points to an N-point Gaussian rule.
  Discussion:
    This subroutine calculates the abscissas and weights of the 2N+1
    point Gauss Kronrod quadrature formula which is obtained from the
    N point Gauss quadrature formula by the optimal addition of N+1 points.
    The optimally added points are called Kronrod abscissas.  The
    abscissas and weights for both the Gauss and Gauss Kronrod rules
    are calculated for integration over the interval [-1,+1].
    Since the quadrature formula is symmetric with respect to the origin,
    only the nonnegative abscissas are calculated.
    Note that the code published in Mathematics of Computation
    omitted the definition of the variable which is here called COEF2.
  Storage:
    Given N, let M = ( N + 1 ) / 2.
    The Gauss-Kronrod rule will include 2*N+1 points.  However, by symmetry,
    only N + 1 of them need to be listed.
    The arrays X, W1 and W2 contain the nonnegative abscissas in decreasing
    order, and the weights of each abscissa in the Gauss-Kronrod and
    Gauss rules respectively.  This means that about half the entries
    in W2 are zero.
    For instance, if N = 3, the output is:
    I      X               W1              W2
    1    0.960491        0.104656         0.000000
    2    0.774597        0.268488         0.555556
    3    0.434244        0.401397         0.000000
    4    0.000000        0.450917         0.888889
    and if N = 4, (notice that 0 is now a Kronrod abscissa)
    the output is
    I      X               W1              W2
    1    0.976560        0.062977        0.000000
    2    0.861136        0.170054        0.347855
    3    0.640286        0.266798        0.000000
    4    0.339981        0.326949        0.652145
    5    0.000000        0.346443        0.000000
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 August 2010
  Author:
    Original FORTRAN77 version by Robert Piessens, Maria Branders.
    C version by John Burkardt.
  Reference:
    Robert Piessens, Maria Branders,
    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
    of Gauss and Lobatto,
    Mathematics of Computation,
    Volume 28, Number 125, January 1974, pages 135-139.
  Parameters:
    Input, int N, the order of the Gauss rule.
    Input, double EPS, the requested absolute accuracy of the
    abscissas.
    Output, double X[N+1], the abscissas.
    Output, double W1[N+1], the weights for
    the Gauss-Kronrod rule.
    Output, double W2[N+1], the weights for
    the Gauss rule.
*/
{
	const dtit3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp eps = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w1 = s_data->a3;
	ityp * w2 = s_data->a4;
	
    ityp ak;
    ityp an;
    ityp *b;
    ityp bb;
    ityp c;
    ityp coef;
    ityp coef2;
    ityp d;
    bool even;
    dim_typ i;
    dim_typ k;
    dim_typ l;
    dim_typ ll;
    dim_typ m;
    ityp s;
    ityp *tau;
    ityp x1;
    ityp xx;
    ityp y;

    b = ( ityp * ) malloc ( (((n+1)/2)+1) * sizeof ( ityp ) );
    tau = ( ityp * ) malloc ( ( (n+1)/2 ) * sizeof ( ityp ) );

    m = ( n + 1 ) / 2;
    even = ( (m<<1) == n );

    d = 2.00;
    an = 0.00;
    for ( k = 1; k <= n; ++k)
    {
        ++ an;
        d *= an / ( an + 0.5 );
    }
    /*
    Calculation of the Chebyshev coefficients of the orthogonal polynomial.
    */
    tau[0] = ( an + 2.00 ) / ( an + an + 3.00 );
    b[m-1] = tau[0] - 1.00;
    ak = an;

    for ( l = 1; l < m; ++l )
    {
        ak += 2.00;
        tau[l] = ( ( ak - 1.00 ) * ak- an * ( an + 1.00 ) ) * ( ak + 2.00 ) * tau[l-1]/ ( ak * ( ( ak + 3.00 ) * ( ak + 2.00 )- an * ( an + 1.00 ) ) );
        b[m-l-1] = tau[l];

        for ( ll = 1; ll <= l; ++l )
            b[m-l-1] += tau[ll-1] * b[m-l+ll-1];
    }

    b[m] = 1.00;
    /*
    Calculation of approximate values for the abscissas.
    */
    bb = sin ( 1.570796 / ( an + an + 1.00 ) );
    x1 = sqrt ( 1.00 - bb * bb );
    s = 2.00 * bb * x1;
    c = sqrt ( 1.00 - s * s );
    coef = 1.00 - ( 1.00 - 1.00 / an ) / ( 8.00 * an * an );
    xx = coef * x1;
    /*
    Coefficient needed for weights.

    COEF2 = 2^(2*n+1) * n! * n! / (2n+1)!
    */
    coef2 = 2.00 / ( ityp ) ( (n<<1) + 1 );
    for ( i = 1; i <= n; i++ )
        coef2 *= 4.00 * ( ityp ) ( i ) / ( ityp ) ( n + i );
    /*
    Calculation of the K-th abscissa (a Kronrod abscissa) and the
    corresponding weight.
    */
    for ( k = 1; k <= n; k +=2 )
    {
        abwe1 ( n, m, eps, coef2, even, b, &xx, w1+k-1 );
        w2[k-1] = 0.00;

        x[k-1] = xx;
        y = x1;
        x1 = y * c - bb * s;
        bb = y * s + bb * c;

        xx = 0.00 + (k!=n)*coef*x1;
        /*
        Calculation of the K+1 abscissa (a Gaussian abscissa) and the
        corresponding weights.
        */
        abwe2 ( n, m, eps, coef2, even, b, &xx, w1+k, w2+k );

        x[k] = xx;
        y = x1;
        x1 = y * c - bb * s;
        bb = y * s + bb * c;
        xx = coef * x1;
    }
    /*
    If N is even, we have one more Kronrod abscissa to compute,
    namely the origin.
    */
    if ( even )
    {
        xx = w2[n] = 0.00;
        abwe1 ( n, m, eps, coef2, even, b, &xx, w1+n );
        x[n] = xx;
    }

    free ( b );
    free ( tau );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _kronrod_adjust ( void * data)
/******************************************************************************/
/*
  Purpose:
    KRONROD_ADJUST adjusts a Gauss-Kronrod rule from [-1,+1] to [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2015
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, the endpoints of the new interval.
    Input, int N, the order of the rule.
    Input/output, double X[N+1], W1[N+1], W2[N+1], the abscissas
    and weights.
*/
{
	const _2itdt3pit * const s_data = data;
	ityp a = s_data->a0;
	ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * w1 = s_data->a4;
	ityp * w2 = s_data->a5;
	
    for ( dim_typ i = 0; i < n + 1; ++i )
    {
        x[i] = ( ( 1.00 - x[i] ) * a+ ( 1.00 + x[i] ) * b )/ 2.00;
        w1[i] = ( ( b - a ) / 2.00 ) * w1[i];
        w2[i] = ( ( b - a ) / 2.00 ) * w2[i];
    }
    return NULL;
}

#endif
