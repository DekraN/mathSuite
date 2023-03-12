#ifndef __DISABLEDEEP_TOEPLITZCHOLESKY

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _toep_cholesky_lower ( void * data)
/******************************************************************************/
/*
  Purpose:
    TOEP_CHOLESKY_LOWER: lower Cholesky factor of a compressed Toeplitz matrix.
  Discussion:
    The Toeplitz matrix A is supplied in a compressed form G.
    Theityp Toeplitz matrix must be positive semi-definite.
    After factorization, A = L * L'.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 November 2012
  Author:
    John Burkardt
  Reference:
    Michael Stewart,
    Cholesky factorization of semi-definite Toeplitz matrices.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double G[2*N], the compressed Toeplitz matrix.
    G(1,1:N) contains the first row.
    G(2,2:N) contains the first column.
    Output, double TOEP_CHOLESKY_LOWER[N*N], the lower Cholesky factor.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * g = s_data->a1;
	
    ityp div;
    ityp g1j;
    ityp g2j;
    dim_typ i, j;
    ityp *l;
    ityp rho;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
        l[i+0*n] = g[0+(i<<1)];

    for ( j = n - 1; 1 <= j; --j )
        g[0+(j<<1)] = g[0+((j-1)<<1)];

    g[0] = 0.00;

    for ( i = 1; i < n; ++i )
    {
        rho = - g[1+(i<<1)] / g[0+(i<<1)];
        div = sqrt ( ( 1.00 - rho ) * ( 1.00 + rho ) );
        for ( j = i; j < n; ++j )
        {
            g1j = g[0+(j<<1)];
            g2j = g[1+(j<<1)];
            g[0+(j<<1)] = (       g1j + rho * g2j ) / div;
        g[1+(j<<1)] = ( rho * g1j +       g2j ) / div;
        }
        for ( j = i; j < n; ++j )
            l[j+i*n] = g[0+(j<<1)];
        for ( j = n - 1; i < j; --j )
            g[0+(j<<1)] = g[0+((j-1)<<1)];
        g[0+(i<<1)] = 0.00;
    }
    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _toep_cholesky_upper ( void * data)
/******************************************************************************/
/*
  Purpose:
    TOEP_CHOLESKY_UPPER: upper Cholesky factor of a compressed Toeplitz matrix.
  Discussion:
    The Toeplitz matrix A is supplied in a compressed form G.
    The Toeplitz matrix must be positive semi-definite.
    After factorization, A = R' * R.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 November 2012
  Author:
    John Burkardt
  Reference:
    Michael Stewart,
    Cholesky factorization of semi-definite Toeplitz matrices.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double G[2*N}, the compressed Toeplitz matrix.
    G(1,1:N) contains the first row.
    G(2,2:N) contains the first column.
    Output, double TOEP_CHOLESKY_UPPER[N*N], the upper Cholesky factor.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * g = s_data->a1;
	
    ityp div;
    ityp g1j;
    ityp g2j;
    dim_typ i, j;
    ityp *r;
    ityp rho;

    r = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < n; ++i )
            r[i+j*n] = 0.00;

    for ( j = 0; j < n; ++j)
        r[0+j*n] = g[0+(j<<1)];
    for ( j = n - 1; 1 <= j; --j )
        g[0+(j<<1)] = g[0+((j-1)<<1)];

    g[0] = 0.00;

    for ( i = 1; i < n; ++i )
    {
        rho = - g[1+(i<<1)] / g[0+(i<<1)];
        div = sqrt ( ( 1.00 - rho ) * ( 1.00 + rho ) );

        for ( j = i; j < n; ++j )
        {
            g1j = g[0+(j<<1)];
            g2j = g[1+(j<<1)];
            g[0+(j<<1)] = (       g1j + rho * g2j ) / div;
            g[1+(j<<1)] = ( rho * g1j +       g2j ) / div;
        }
        for ( j = i; j < n; ++j )
            r[i+j*n] = g[0+(j<<1)];
        for ( j = n - 1; i < j; --j )
            g[0+(j<<1)] = g[0+((j-1)<<1)];
        g[0+(i<<1)] = 0.00;
    }

    return r;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _toeplitz_cholesky_lower ( void * data)
/******************************************************************************/
/*
  Purpose:
    TOEPLITZ_CHOLESKY_LOWER: lower Cholesky factor of a Toeplitz matrix.
  Discussion:
    The Toeplitz matrix must be positive semi-definite.
    After factorization, A = L * L'.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 November 2012
  Author:
    John Burkardt
  Reference:
    Michael Stewart,
    Cholesky factorization of semi-definite Toeplitz matrices.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the Toeplitz matrix.
    Output, double TOEPLITZ_CHOLESKY_LOWER[N*N], the lower Cholesky factor.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp div;
    ityp *g;
    ityp g1j;
    ityp g2j;
    dim_typ i, j;
    ityp *l;
    ityp rho;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    g = ( ityp * ) malloc ( n * sizeof ( ityp ) << 1 );

    for ( j = 0; j < n; ++j )
        g[0+(j<<1)] = a[0+j*n];
    g[1] = 0.00;
    for ( j = 1; j < n; ++j )
        g[1+(j<<1)] = a[j+0*n];

    for ( i = 0; i < n; ++i )
    l[i+0*n] = g[0+(i<<1)];
    for ( j = n - 1; 1 <= j; --j )
        g[0+(j<<1)] = g[0+((j-1)<<1)];

    g[0] = 0.00;

    for ( i = 1; i < n; ++i )
    {
        rho = - g[1+(i<<1)] / g[0+(i<<1)];
        div = sqrt ( ( 1.00 - rho ) * ( 1.00 + rho ) );

        for ( j = i; j < n; ++j)
        {
            g1j = g[0+(j<<1)];
            g2j = g[1+(j<<1)];
            g[0+(j<<1)] = (       g1j + rho * g2j ) / div;
            g[1+(j<<1)] = ( rho * g1j +       g2j ) / div;
        }

        for ( j = i; j < n; ++j )
            l[j+i*n] = g[0+(j<<1)];
        for ( j = n - 1; i < j; --j )
            g[0+(j<<1)] = g[0+((j-1)<<1)];
            g[0+(i<<1)] = 0.00;
    }

    free ( g );

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _toeplitz_cholesky_upper ( void * data)
/******************************************************************************/
/*
  Purpose:
    TOEPLITZ_CHOLESKY_UPPER: upper Cholesky factor of a Toeplitz matrix.
  Discussion:
    The Toeplitz matrix must be positive semi-definite.
    After factorization, A = R' * R.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 November 2012
  Author:
    John Burkardt
  Reference:
    Michael Stewart,
    Cholesky factorization of semi-definite Toeplitz matrices.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the Toeplitz matrix.
    Output, double TOEPLITZ_CHOLESKY_UPPER[N*N], the upper Cholesky factor.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp div;
    ityp *g;
    ityp g1j;
    ityp g2j;
    dim_typ i, j;
    ityp *r;
    ityp rho;

    r = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            r[i+j*n] = 0.00;

    g = ( ityp * ) malloc ( n * sizeof ( ityp ) << 1);

    for ( j = 0; j < n; ++j )
    g[0+(j<<1)] = a[0+j*n];

    g[1] = 0.00;

    for ( j = 1; j < n; ++j )
        g[1+(j<<1)] = a[j+0*n];
    for ( j = 0; j < n; ++j )
        r[0+j*n] = g[0+(j<<1)];
    for ( j = n - 1; 1 <= j; --j )
        g[0+(j<<1)] = g[0+((j-1)<<1)];

    g[0] = 0.00;

    for ( i = 1; i < n; ++i )
    {
        rho = - g[1+(i<<1)] / g[0+(i<<1)];
        div = sqrt ( ( 1.00 - rho ) * ( 1.00 + rho ) );
        for ( j = i; j < n; j++ )
        {
            g1j = g[0+(j<<1)];
            g2j = g[1+(j<<1)];
            g[0+(j<<1)] = (       g1j + rho * g2j ) / div;
            g[1+(j<<1)] = ( rho * g1j +       g2j ) / div;
        }
        for ( j = i; j < n; ++j )
            r[i+j*n] = g[0+(j<<1)];
        for ( j = n - 1; i < j; --j )
            g[0+(j<<1)] = g[0+((j-1)<<1)];
            g[0+(i<<1)] = 0.00;
    }

    free ( g );

    return r;
}

#endif
