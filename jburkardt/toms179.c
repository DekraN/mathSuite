#ifndef __DISABLEDEEP_TOMS179

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mdbeta ( void * data)
/******************************************************************************/
/*
  Purpose:
    MDBETA evaluates the incomplete beta function.
  Modified:
    30 January 2008
  Author:
    Oliver Ludwig
    Modifications by John Burkardt
  Reference:
    Oliver Ludwig,
    Algorithm 179:
    Incomplete Beta Ratio,
    Communications of the ACM,
    Volume 6, Number 6, June 1963, page 314.
  Parameters:
    Input, double X, the value to which function is to be
    integrated.  X must be in the range [0,1] inclusive.
    Input, double P, the first parameter.  P must be greater
    than 0.0.
    Input, double Q, the second parameter.  Q must be greater
    than 0.0.
    Output, int *IER, error parameter.
    0, normal exit.
    1, X is not in the range [0,1] inclusive.
    2, P or Q is less than or equal to 0.
    Output, double MDBETA.  The probability that a random variable
    from a Beta distribution having parameters P and Q will be less than
    or equal to X.
  Local parameters:
    Local, double ALEPS, the logarithm of EPS1.
    Local, double EPS, the machine precision.
    Local, double EPS1, the smallest representable number.
*/
{
	static ityp result = MAX_VAL;
	
	const _3itpi * const s_data = data; 
	const register ityp x = s_data->a0;
	register ityp p = s_data->a1;
	register ityp q = s_data->a2;
	int * ier = s_data->a3;
	
    ityp aleps = - 179.6016;
    ityp c;
    ityp cnt;
    ityp d4;
    ityp dp;
    ityp dq;
    ityp eps = 2.2E-16;
    ityp eps1 = 1.0E-78;
    ityp finsum;
    dim_typ ib;
    ityp infsum;
    int interval;
    ityp p1;
    ityp pq;
    ityp prob;
    ityp ps;
    ityp px;
    ityp temp;
    ityp wh;
    ityp xb;
    ityp y;
    /*
    Check ranges of the arguments.
    */
    prob = 0.00;
    y = x;

    if ( x < 0.00 || 1.00 < x || p <= 0.00 || q <= 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( x <= 0.50 )
        interval = 0;
    else
    {
        interval = 1;
        temp = p;
        p = q;
        q = temp;
        y = 1.00 - y;
    }

    if ( x == 0.00 || x == 1.00 )
    {
        prob = 0.00;

        if ( interval != 0 )
        {
            prob = 1.00 - prob;
            temp = p;
            p = q;
            q = temp;
        }
        
        result = prob;
        return &result;
    }

    ib = q;
    temp = ib;
    ps = q - ( ityp ) ( ib );

    if ( q == temp )
        ps = 1.00;

    dp = p;
    dq = q;
    px = dp * log ( y );
    pq = alogam ( dp + dq );
    p1 = alogam ( dp );
    c = alogam ( dq );
    d4 = log ( dp );
    xb = px + alogam ( ps + dp ) - alogam ( ps ) - d4 - p1;
    /*
    Scaling
    */
    ib = ( dim_typ ) ( xb / aleps );
    infsum = 0.00;
    /*
    First term of a decreasing series will underflow.
    */
    if ( ib == 0 )
    {
        infsum = exp ( xb );
        cnt = infsum * dp;
        /*
        CNT will equal dexp ( temp ) * ( 1.d0 - ps ) * i * p * y**i / factorial ( i ).
        */
        wh = 0.00;

        for ( ; ; )
        {
            ++ wh;
            cnt *= ( wh - ps ) * y / wh;
            xb = cnt / ( dp + wh );
            infsum += xb;

            if ( xb / eps < infsum )
                break;
        }
    }

    finsum = 0.00;

    if ( dq <= 1.00 )
    {
        prob = finsum + infsum;

        if ( interval != 0 )
        {
            prob = 1.00 - prob;
            temp = p;
            p = q;
            q = temp;
        }

		result = prob;
        return &result;
    }

    xb = px + dq * log ( 1.00 - y ) + pq - p1 - log ( dq ) - c;
    /*
    Scaling.
    */
    ib = ( dim_typ ) ( xb / aleps );

    if ( ib < 0 )
        ib = 0;

    c = 1.00 / ( 1.00 - y );
    cnt = exp ( xb - ( ityp ) ( ib ) * aleps );
    ps = dq;
    wh = dq;

    for ( ; ; )
    {
        -- wh;

        if ( wh <= 0.00 )
        {
            prob = finsum + infsum;

            if ( interval != 0 )
            {
                prob = 1.00 - prob;
                temp = p;
                p = q;
                q = temp;
            }
            break;
    }

    px = ( ps * c ) / ( dp + wh );

    if ( px <= 1.00 )
        if ( cnt / eps <= finsum || cnt <= eps1 / px )
        {
            prob = finsum + infsum;

            if ( interval != 0 )
            {
                prob = 1.00 - prob;
                temp = p;
                p = q;
                q = temp;
            }
            break;
        }
    cnt *= px;
    /*
    Rescale.
    */
    if ( 1.00 < cnt )
    {
        -- ib;
        cnt *= eps1;
    }

    ps = wh;

    if ( ib == 0 )
        finsum += cnt;
    }

	result = prob;
    return &result;
}

#endif
