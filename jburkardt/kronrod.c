#ifndef __DISABLEDEEP_KRONROD

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _abwe1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ABWE1 calculates a Kronrod abscissa and weight.
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
    Input, int M, the value of ( N + 1 ) / 2.
    Input, double EPS, the requested absolute accuracy of the
    abscissas.
    Input, double COEF2, a value needed to compute weights.
    Input, int EVEN, is TRUE if N is even.
    Input, double B[M+1], the Chebyshev coefficients.
    Input/output, double *X; on input, an estimate for
    the abscissa, and on output, the computed abscissa.
    Output, double *W, the weight.
*/
{
	const _2dt2iti3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp eps = s_data->a2;
	ityp coef2 = s_data->a3;
	int even = s_data->a4;
	ityp * b = s_data->a5;
	ityp * x = s_data->a6;
	ityp * w = s_data->a7;
	
    ityp ai;
    ityp b0;
    ityp b1;
    ityp b2;
    ityp d0;
    ityp d1;
    ityp d2;
    ityp delta;
    ityp dif;
    ityp f;
    ityp fd;
    dim_typ i;
    dim_typ iter;
    dim_typ k;
    dim_typ ka;
    ityp yy;

    ka = *x == 0.00;
    /*
    Iterative process for the computation of a Kronrod abscissa.
    */
    for ( iter = 1; iter <= 50; ++iter)
    {
        b1 = d1 = 0.00;
        b2 = b[m];
        yy = 4.00 * (*x) * (*x) - 2.00;

        if ( even )
        {
            ai = m + m + 1;
            d2 = ai * b[m];
            dif = 2.00;
        }
        else
        {
            ai = m + 1;
            d2 = 0.00;
            dif = 1.00;
        }

        for ( k = 1; k <= m; ++k )
        {
            ai -= dif;
            i = m - k + 1;
            b0 = b1;
            b1 = b2;
            d0 = d1;
            d1 = d2;
            b2 = yy * b1 - b0 + b[i-1];
            if ( !even )
                ++ i;
            d2 = yy * d1 - d0 + ai * b[i-1];
        }

        if ( even )
        {
            f = ( *x ) * ( b2 - b1 );
            fd = d2 + d1;
        }
        else
        {
            f = 0.50 * ( b2 - b0 );
            fd = 4.00 * ( *x ) * d2;
        }
        /*
        Newton correction.
        */
        delta = f / fd;
        *x -= delta;

        if ( ka == 1 )
            break;

        if ( abs ( delta ) <= eps )
            ka = 1;
    }
    /*
    Catch non-convergence.
    */
    if ( ka != 1 )
        return NULL;
    /*
    Computation of the weight.
    */
    d0 = 1.00;
    d1 = *x;
    ai = 0.00;
    for ( k = 2; k <= n; ++k )
    {
        ++ ai;
        d2 = ( ( ai + ai + 1.00 ) * ( *x ) * d1 - ai * d0 ) / ( ai + 1.0 );
        d0 = d1;
        d1 = d2;
    }

    *w = coef2 / ( fd * d2 );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _abwe2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ABWE2 calculates a Gaussian abscissa and two weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 April 2013
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
    Input, int M, the value of ( N + 1 ) / 2.
    Input, double EPS, the requested absolute accuracy of the
    abscissas.
    Input, double COEF2, a value needed to compute weights.
    Input, int EVEN, is TRUE if N is even.
    Input, double B[M+1], the Chebyshev coefficients.
    Input/output, double *X; on input, an estimate for
    the abscissa, and on output, the computed abscissa.
    Output, double *W1, the Gauss-Kronrod weight.
    Output, double *W2, the Gauss weight.
*/
{
	const _2dt2iti4pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp eps = s_data->a2;
	ityp coef2 = s_data->a3;
	int even = s_data->a4;
	ityp * b = s_data->a5;
	ityp * x = s_data->a5;
	ityp * w1 = s_data->a6;
	ityp * w2 = s_data->a7;
	
    ityp ai;
    ityp an;
    ityp delta;
    dim_typ i;
    dim_typ iter;
    dim_typ k;
    dim_typ ka;
    ityp p0;
    ityp p1;
    ityp p2;
    ityp pd0;
    ityp pd1;
    ityp pd2;
    ityp yy;

    ka = *x == 0.00;
    /*
    Iterative process for the computation of a Gaussian abscissa.
    */
    for ( iter = 1; iter <= 50; ++iter )
    {
        p0 = pd1 = 1.00;
        p1 = *x;
        pd0 = 0.00;
        /*
        When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
        */
        if ( n <= 1 )
        {
            if ( r8_epsilon ( ) < abs ( *x ) )
            {
                p2 = ( 3.00 * ( *x ) * ( *x ) - 1.00 ) / 2.00;
                pd2 = 3.0 * ( *x );
            }
            else
            {
                p2 = 3.00 * ( *x );
                pd2 = 3.00;
            }
        }

    ai = 0.00;
    for ( k = 2; k <= n; ++k )
    {
        ++ ai;
        p2 = ( ( ai + ai + 1.00 ) * (*x) * p1 - ai * p0 ) / ( ai + 1.00 );
        pd2 = ( ( ai + ai + 1.00 ) * ( p1 + (*x) * pd1 ) - ai * pd0 )/ ( ai + 1.0 );
        p0 = p1;
        p1 = p2;
        pd0 = pd1;
        pd1 = pd2;
    }
    /*
    Newton correction.
    */
    delta = p2 / pd2;
    *x -= delta;

    if ( ka == 1 )
        break;

    if ( abs ( delta ) <= eps )
        ka = 1;
    }
    /*
    Catch non-convergence.
    */
    if ( ka != 1 )
        return NULL;
    /*
    Computation of the weight.
    */
    an = n;

    *w2 = 2.00 / ( an * pd2 * p0 );

    p1 = 0.00;
    p2 = b[m];
    yy = 4.00 * (*x) * (*x) - 2.00;
    for ( k = 1; k <= m; ++k )
    {
        i = m - k + 1;
        p0 = p1;
        p1 = p2;
        p2 = yy * p1 - p0 + b[i-1];
    }

    *w1 = *w2 + even ? coef2 / ( pd2 * (*x) * ( p2 - p1 ) ) : 2.00 * coef2 / ( pd2 * ( p2 - p0 ) );
    return NULL;
}

#endif
