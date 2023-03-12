#ifndef __DISABLEDEEP_TOMS322

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fisher ( void * data)
/******************************************************************************/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register ityp x = s_data->a2;
	
    dim_typ a;
    dim_typ b;
    ityp d;
    dim_typ i, j;
    ityp p;
    ityp w;
    ityp y;
    ityp z;
    ityp zk;

    a = (( m / 2 )<<1) - m + 2;
    b = (( n / 2 )<<1) - n + 2;
    w = ( x * m ) / n;
    z = 1.00 / ( 1.00 + w );

    if ( a == 1)
    {
        if ( b == 1 )
        {
            p = sqrt ( w );
            y = 0.3183098862;
            d = y * z / p;
            p = 2.00 * y * atan ( p );
        }
        else
        {
            p = sqrt ( w * z );
            d = 0.50 * p * z / w;
        }
    }
    else if ( b == 1 )
    {
        p = sqrt ( z );
        d = 0.50 * z * p;
        p = 1.00 - p;
    }
    else
    {
        d = z * z;
        p = w * z;
    }
    y = 2.00 * w / z;
    if ( a == 1 )
        for(j = b + 2; j <= n; j += 2 )
        {
            d *= ( 1.00 + 1.00 / ( j - 2 ) ) * z;
            p += d * y / ( j - 1 );
        }
    else
    {
        zk = pow ( z, ( ityp ) ( ( n - 1 ) / 2 ) );
        d *= ( zk * n ) / b;
        p = p * zk + w * z * ( zk - 1.00 ) / ( z - 1.00 );
    }
    y = w * z;
    z = 2.00 / z;
    b = n - 2;

    for ( i = a + 2; i <= m; i += 2 )
    {
        j = i + b;
        d *= ( y * j ) / ( i - 2 );
        p -= z * d / j;
    }

    if ( p < 0.00 )
       p = 0.00;
    if ( 1.00 < p )
        p = 1.00;

	result = p;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  inline void *   _student ( void * data)
/******************************************************************************/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ df = s_data->a0;
	const register ityp t = s_data->a1;
	
	result = fisher ( 1, df, t*t );
    return &result;
}

#endif
