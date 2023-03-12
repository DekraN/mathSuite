#ifndef __DISABLEDEEP_MGS

#include "../dutils.h"

__MATHSUITE __JBURKARDT  void *   _mgs ( void * data)
{
	const _2dt3ppit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp ** c = s_data->a2;
	ityp ** r = s_data->a3;
	ityp ** q = s_data->a4;
	
    dim_typ i, j, k;
    ityp *x;
    ityp xn;

    x = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    for ( k = 0; k < n; ++k )
    {
        for ( j = 0; j < m; ++j)
            x[j] = c[j][k];
        xn = 0;
        for ( j = 0; j < m; ++j )
            xn += x[j]*x[j];
        r[k][k] = sqrt(xn);
        if ( 0.00 < r[k][k] )
            for ( j = 0; j < m; ++j )
                q[j][k] = c[j][k]/r[k][k];
        else
            for ( j = 0; j < m; ++j )
                q[j][k] = 0.00;
        for ( j = k + 1; j < n; j++ )
        {
            r[k][j] = 0;
            for ( i = 0; i < m; ++i)
                r[k][j] +=  q[i][k]*c[i][j];
            for ( i = 0; i < m; ++i )
                c[i][j] -= q[i][k]*r[k][j];
        }
    }
    free ( x );
    return NULL;
}

#endif
