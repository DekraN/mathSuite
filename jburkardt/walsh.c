#ifndef __DISABLEDEEP_WALSH

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ffwt ( void * data)
/******************************************************************************/
/*
  Purpose:
    FFWT performs an in-place fast Walsh transform.
  Discussion:
    This routine performs a fast Walsh transform on an input series X
    leaving the transformed results in X.
    X is dimensioned N, which must be a power of 2.
    The results of this Walsh transform are in sequency order.
    The output sequence could be normalized by dividing by N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 March 2011
  Author:
    Ken Beauchamp
  Reference:
    Ken Beauchamp,
    Walsh functions and their applications,
    Academic Press, 1975,
    ISBN: 0-12-084050-2,
    LC: QA404.5.B33.
  Parameters:
    Input, int N, the number of items in X.
    N must be a power of 2.
    Input/output, double X[N], the data to be transformed.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp hold;
    dim_typ i;
    dim_typ ii;
    dim_typ j;
    int j2;
    int js;
    dim_typ k;
    int l;
    int m;
    int mw;
    int mw1;
    int nw;
    int nz;
    int nz2;
    int nzi;
    int nzn;
    int *two_power;
    ityp z;

    m = i4_log_2 ( n );

    two_power = ( int * ) malloc ( m * sizeof ( int ) );

    for ( i = 0; i < m; ++i )
        two_power[i] = powi ( 2, m - 1 - i );

    for ( l = 0; l < m; ++l )
    {
        nz = powi ( 2, l );
        nzi = nz << 1;
        nzn = n / nzi;
        nz2 = nz / 2;
        if ( nz2 == 0 )
            nz2 = 1;

        for ( i = 0; i < nzn; ++i )
        {
            js = i * nzi;
            z = 1.00;
            for ( ii = 0; ii < 2; ++ii )
            {
                for ( j = 0; j < nz2; ++j )
                {
                    ++ js;
                    j2 = js + nz;
                    hold = x[js-1] + z * x[j2-1];
                    z *= -1;
                    x[j2-1] = x[js-1] + z * x[j2-1];
                    x[js-1] = hold;
                    z *= -1;
                }
                if ( l == 0 )
                    break;
                z = - 1.00;
            }
        }
    }
    /*
    Bit reversal section.
    */
    nw = 0;
    for ( k = 0; k < n; ++k )
    {
        /*
        Choose correct index and switch elements if not already switched.
        */
        if ( k < nw )
        {
            hold = x[nw];
            x[nw] = x[k];
            x[k] = hold;
        }
        /*
        Bump up series by 1.
        */
        for ( i = 0; i < m; ++i )
        {
            ii = i;
            if ( nw < two_power[i] )
                break;
            mw = nw / two_power[i];
            mw1 = mw / 2;
            if ( mw <= mw1<<1 )
                break;
            nw -= two_power[i];
        }
        nw += two_power[ii];
    }

    free ( two_power );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _fwt ( void * data)
/******************************************************************************/
/*
  Purpose:
    FWT performs a fast Walsh transform.
  Discussion:
    This routine performs a fast Walsh transform on an input series X
    leaving the transformed results in X.
    X is dimensioned N, which must be a power of 2.
    The results of this Walsh transform are in sequency order.
    The output sequence could be normalized by dividing by N.
    Note that the program text in the reference included the line
      y(jd) = abs ( x(j) - x(j2) )
    which has been corrected to:
      y(jd) = x(j) - x(j2)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 March 2011
  Author:
    Ken Beauchamp
  Reference:
    Ken Beauchamp,
    Walsh functions and their applications,
    Academic Press, 1975,
    ISBN: 0-12-084050-2,
    LC: QA404.5.B33.
  Parameters:
    Input, int N, the number of items in X.
    N must be a power of 2.
    Input/output, double X[N], the data to be transformed.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i;
    dim_typ j;
    int j2;
    dim_typ jd;
    dim_typ js;
    dim_typ l;
    int m;
    int n2;
    int nx;
    int ny;
    int nz;
    int nzi;
    int nzn;
    ityp *y = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    n2 = n / 2;
    m = i4_log_2 ( n );

    for ( l = 1; l <= m; ++l )
    {
        ny = 0;
        nz = powi ( 2, l - 1 );
        nzi = nz<<1;
        nzn = n / nzi;
        for ( i = 1; i <= nzn; ++i)
        {
            nx = ny + 1;
            ny += nz;
            js = ( i - 1 ) * nzi;
            jd = js + nzi + 1;
            for ( j = nx; j <= ny; ++j )
            {
                ++ js;
                j2 = j + n2;
                y[js-1] = x[j-1] + x[j2-1];
                -- jd;
                y[jd-1] = x[j-1] - x[j2-1];
            }
        }
        r8vec_copy ( n, y, x );
    }
    free ( y );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _walsh ( void * data)
/******************************************************************************/
/*
  Purpose:
    WALSH performs a fast Walsh transform.
  Discussion:
    This routine performs a fast Wash transform on an input series X
    leaving the transformed results in X.  The array Y is working space.
    X and Y are dimensioned N, which must be a power of 2.
    The results of this Walsh transform are in sequency order.
    The output sequence could be normalized by dividing by N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 March 2011
  Author:
    Ken Beauchamp
  Reference:
    Ken Beauchamp,
    Walsh functions and their applications,
    Academic Press, 1975,
    ISBN: 0-12-084050-2,
    LC: QA404.5.B33.
  Parameters:
    Input, int N, the number of items in X.
    N must be a power of 2.
    Input/output, double X[N], the data to be transformed.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp a;
    dim_typ i;
    dim_typ i1;
    int is;
    dim_typ j;
    int j1;
    int l;
    int m;
    int n1;
    int n2;
    ityp w;
    ityp *y;
    ityp z;

    n2 = n / 2;
    y = ( ityp * ) malloc ( n2 * sizeof ( ityp ) );
    m = i4_log_2 ( n );
    z = - 1.00;

    for ( j = 1; j <= m; ++j )
    {
        n1 = powi ( 2, m - j + 1 );
        j1 = powi ( 2, j - 1 );
        for ( l = 1; l <= j1; ++l )
        {
            is = ( l - 1 ) * n1 + 1;
            i1 = 0;
            w = z;
            for ( i = is; i <= is + n1 - 1; i += 2 )
            {
                a = x[i-1];
                x[is+i1-1] = a + x[i];
                ++ i1;
                y[i1-1] = ( x[i] - a ) * w;
                w *= z;
            }
            for ( i = 1; i <= n1 / 2; ++i )
                x[n1/2+is+i-2] = y[i-1];
        }
    }

    free ( y );

    return NULL;
}

#endif
