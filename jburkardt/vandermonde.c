#ifndef __DISABLEDEEP_VANDERMONDE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _bivand1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    BIVAND1 returns a bidimensional Vandermonde1 matrix.
  Discussion:
    N = 3, ALPHA = ( 1, 2, 3 ), BETA = ( 10, 20, 30 )
 (x,y)   | (1,10) (2,10) (3,10) (1,20) (2,20) (1,30)
    --------+-----------------------------------------------
    1       |     1       1       1       1       1       1
    x       |     1       2       3       1       2       1
       y    |    10      10      10      20      20      30
    x^2     |     1       4       9       1       4       1
    x  y    |    10      20      30      20      40      30
    x^2y^2  |   100     100     100     400     400     900
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 February 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the data vectors.
    Input, double ALPHA[N], BETA[N], the values that define A.
    Output, double BIVAND1[((N+1)*N)/2*((N+1)*N)/2], the matrix.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * beta = s_data->a2;
	
    ityp *a;
    dim_typ e;
    dim_typ e1;
    dim_typ e2;
    int i1;
    int i2;
    dim_typ ii;
    dim_typ j1;
    dim_typ j2;
    dim_typ jj;
    int n2;

    n2 = ( n * ( n + 1 ) ) / 2;
    a = ( ityp * ) malloc ( n2 * n2 * sizeof ( ityp ) );

    e1 = e2 = e = 0;

    for ( ii = 0; ii < n2; ++ii)
    {
        j1 = j2 = 0;
        for ( jj = 0; jj < n2; ++jj )
        {
            a[ii+jj*n2] = ii == 0 ? 1.00 : pow ( alpha[j1], e1 ) * pow ( beta[j2], e2 );

            if ( j1 + j2 < n - 1 )
                ++ j1;
            else
            {
                j1 = 0;
                ++ j2;
            }
        }

        if ( e2 < e )
        {
            -- e1;
            ++ e2;
        }
        else
        {
            ++ e;
            e1 = e;
            e2 = 0;
        }
    }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _bivand2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    BIVAND2 returns a bidimensional Vandermonde1 matrix.
  Discussion:
    N = 3, ALPHA = ( 1, 2, 3 ), BETA = ( 10, 20, 30 )
 (x,y)   | (1,10) (2,10) (3,10) (1,20) (2,20) (3,20) (1,30) (2,30) (3,30)
    --------+---------------------------------------------------------------
    1       |     1      1      1      1      1      1      1      1      1
    x       |     1      2      3      1      2      1      1      2      3
    x^2     |     1      4      9      1      4      1      1      4      9
       y    |    10     10     10     20     20     20     30     30     30
    x  y    |    10     20     30     20     40     60     30     60     90
    x^2y    |    10     40     90     20     80    180     30    120    270
       y^2  |   100    100    100    400    400    400    900    900    900
    x  y^2  |   100    200    300    400    800   1200    900   1800   2700
    x^2y^2  |   100    400    900    400   1600   3600    900   3600   8100
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 May 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the data vectors.
    Input, double ALPHA[N], BETA[N], the values that define A.
    Output, double BIVAND2[(N*N)*(N*N)], the matrix.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * beta = s_data->a2;
	
    ityp *a;
    dim_typ i;
    dim_typ ix;
    dim_typ iy;
    dim_typ j;
    dim_typ jx;
    dim_typ jy;

    a = ( ityp * ) malloc ( n * n * n * n * sizeof ( ityp ) );

    i = 0;
    for ( iy = 0; iy < n; ++iy )
        for ( ix = 0; ix < n; ++ix )
        {
            j = 0;
            for ( jy = 0; jy < n; ++jy )
                for ( jx = 0; jx < n; ++jx )
                {
                    a[i+j*n*n] = pow ( alpha[jx], ix ) * pow ( beta[jy], iy );
                    ++ j;
                }
            ++ i;
        }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _dvand ( void * data)
/******************************************************************************/
/*
  Purpose:
    DVAND solves a Vandermonde system A' * x = b.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 February 2014
  Author:
    John Burkardt
  Reference:
    Ake Bjorck, Victor Pereyra,
    Solution of Vandermonde Systems of Equations,
    Mathematics of Computation,
    Volume 24, Number 112, October 1970, pages 893-903.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double ALPHA[N], the parameters that define the matrix.
    The values should be distinct.
    Input, double B[N], the right hand side of the linear system.
    Output, double DVAND[N], the solution of the linear system.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * b = s_data->a2;
	
    dim_typ j;
    dim_typ k;
    ityp *x;

    x = r8vec_copy_new ( n, b );

    for ( k = 0; k < n - 1; ++k )
        for ( j = n - 1; k < j; --j)
            x[j] = ( x[j] - x[j-1] ) / ( alpha[j] - alpha[j-k-1] );

    for ( k = n - 2; 0 <= k; --k )
        for ( j = k; j < n - 1; ++j )
            x[j] -= alpha[k] * x[j+1];


    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dvandprg ( void * data)
/******************************************************************************/
/*
  Purpose:
    DVANDPRG solves a Vandermonde system A' * x = f progressively.
  Discussion:
    This function receives the solution to the system of equations A' * x = f
    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
    and new values alpha(n) and f(n).  It updates the solution.
    To solve a system of Nbig equations, this function may be called
    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the
    current subsystem is returned.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 April 2014
  Author:
    John Burkardt
  Reference:
    Ake Bjorck, Victor Pereyra,
    Solution of Vandermonde Systems of Equations,
    Mathematics of Computation,
    Volume 24, Number 112, October 1970, pages 893-903.
  Parameters:
    Input, int N, the new order of the matrix, which is 1
    larger than on the previous call.  For the first call, N must be 1.
    Input, double ALPHA[N], the parameters that define the matrix.
    The values should be distinct.  The value ALPHA(N) has just been
    added to the system.
    Input, double B[N], the right hand side of the linear system.
    Input/output, double X[N].  On input, the first N-1 entries
    contain the solution of the N-1xN-1 linear system.  On output, the
    solution to the NxN linear system.
    Input/output, double C[N], M[N].  On input, the first N-1
    entries contain factorization data for the N-1xN-1 linear system.  On
    output, factorization data for the NxN linear system.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x = s_data->a3;
	ityp * c = s_data->a4;
	ityp * m = s_data->a5;
	
    ityp cn;
    dim_typ j;

    c[n-1] = b[n-1];
    for ( j = n - 1; 1 <= j; --j)
        c[j-1] = ( c[j] - c[j-1] ) / ( alpha[n-1] - alpha[j-1] );

    m[n-1] = n == 1;

    cn = c[0];
    x[n-1] = c[0];

    for ( j = n - 1; 1 <= j; --j )
    {
        m[j] -= alpha[n-2] * m[j-1];
        x[n-j-1] += m[j] * cn;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pvand ( void * data)
/******************************************************************************/
/*
  Purpose:
    PVAND solves a Vandermonde system A * x = b.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 February 2014
  Author:
    John Burkardt
  Reference:
    Ake Bjorck, Victor Pereyra,
    Solution of Vandermonde Systems of Equations,
    Mathematics of Computation,
    Volume 24, Number 112, October 1970, pages 893-903.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double ALPHA[N], the parameters that define the matrix.
    The values should be distinct.
    Input, double B[N], the right hand side of the linear system.
    Output, double PVAND[N], the solution of the linear system.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * b = s_data->a2;
	
    dim_typ j;
    dim_typ k;
    ityp *x;

    x = r8vec_copy_new ( n, b );

    for ( k = 0; k < n - 1; ++k )
        for ( j = n - 1; k < j; --j )
            x[j] -= alpha[k] * x[j-1];

    for ( k = n - 2; 0 <= k; --k )
    {
        for ( j = k + 1; j < n; ++j )
            x[j] /= (alpha[j] - alpha[j-k-1] );
        for ( j = k; j < n - 1; ++j )
            x[j] -= x[j+1];
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pvandprg ( void * data)
/******************************************************************************/
/*
  Purpose:
    PVANDPRG solves a Vandermonde system A * x = f progressively.
  Discussion:
    This function receives the solution to the system of equations A * x = f
    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
    and new values alpha(n) and f(n).  It updates the solution.
    To solve a system of Nbig equations, this function may be called
    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the
    current subsystem is returned.
    Note that the reference, which lists an Algol version of this algorithm,
    omits a minus sign, writing
      u[j] := u[j] x delta;
    where
      u[j] := - u[j] x delta;
    is actually necessary.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 April 2014
  Author:
    John Burkardt
  Reference:
    Ake Bjorck, Victor Pereyra,
    Solution of Vandermonde Systems of Equations,
    Mathematics of Computation,
    Volume 24, Number 112, October 1970, pages 893-903.
  Parameters:
    Input, int N, the new order of the matrix, which is 1
    larger than on the previous call.  For the first call, N must be 1.
    Input, double ALPHA[N], the parameters that define the matrix.
    The values should be distinct.  The value ALPHA(N) has just been
    added to the system.
    Input, double B[N], the right hand side of the linear system.
    Input/output, double X[N]; on input, the solution of the
    N-1xN-1 linear system.  On output, the solution of the NxN linear system.
    Input/output, double D[N], U[N]; on input, factorization data
    for the N-1xN-1 linear system.  On output, factorization data for the
    NxN linear system.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x = s_data->a3;
	ityp * d = s_data->a4;
	ityp * u = s_data->a5;
	
    ityp delta;
    ityp dn;
    dim_typ j;

    d[n-1] = b[n-1];
    for ( j = n - 1; 1 <= j; --j )
        d[j-1] = d[j] - alpha[n-j-1] * d[j-1];

    dn = d[0];
    u[n-1] = 1.00;

    for ( j = 1; j <= n - 1; ++j )
    {
        delta = alpha[n-1] - alpha[j-1];
        u[j-1] *= - delta;
        u[n-1] *= delta;
        x[j-1] += dn / u[j-1];
    }

    x[n-1] = dn / u[n-1];

    return NULL;
}

#endif
