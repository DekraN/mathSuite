#ifndef __DISABLEDEEP_LAPLACIAN

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cholesky_upper_error ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHOLESKY_UPPER_ERROR determines the error in an upper Cholesky factor.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the matrix.
    Input, double C[N*N], the upper triangular Cholesky factor.
    Output, double CHOLESKY_UPPER_ERROR, the Frobenius norm
    of the difference matrix A - C' * C.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * c = s_data->a2;
	
    ityp *ctc;
    ityp *d;
    ityp value;
    ctc = r8mat_mtm_new ( n, n, n, c, c );
    d = r8mat_sub_new ( n, n, a, ctc );
    value = r8mat_norm_fro ( n, n, d );
    free ( ctc );
    free ( d );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_mtm_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_MTM_NEW computes C = A' * B.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
    For this routine, the result is returned as the function value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, N3, the order of the matrices.
    Input, double A[N2*N1], double B[N2*N3], the matrices to multiply.
    Output, double R8MAT_MTM_NEW[N1*N3], the product matrix C = A' * B.
*/
{
	const _3dt2pit * const s_data = data;
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	const register dim_typ n3 = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	
	ityp *c;
	dim_typ i, j, k;
	
	c = ( ityp * ) malloc ( n1 * n3 * sizeof ( ityp ) );
	
	for ( i = 0; i < n1; ++i )
		for ( j = 0; j < n3; ++j )
		{
			c[i+j*n1] = 0.00;
			for ( k = 0; k < n2; ++k )
			c[i+j*n1] += a[k+i*n2] * b[k+j*n2];
		}


	return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _eigen_error ( void * data)
/******************************************************************************/
/*
  Purpose:
    EIGEN_ERROR determines the error in a (right) eigensystem.
  Discussion:
    An R8MAT is a matrix of double values.
    This routine computes the Frobenius norm of
      A * X - X * LAMBDA
    where
      A is an N by N matrix,
      X is an N by K matrix (each of K columns is an eigenvector)
      LAMBDA is a K by K diagonal matrix of eigenvalues.
    This routine assumes that A, X and LAMBDA are all real!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, int K, the number of eigenvectors.
    K is usually 1 or N.
    Input, double A[N*N], the matrix.
    Input, double X[N*K]), the K eigenvectors.
    Input, double LAMBDA[K], the K eigenvalues.
    Output, double EIGEN_ERROR, the Frobenius norm
    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	ityp * a = s_data->a2;
	ityp * x = s_data->a3;
	ityp * lambda = s_data->a4;
	
    ityp *c;
    dim_typ i, j;
    ityp value;

    c = r8mat_mm_new ( n, n, k, a, x );

    for ( j = 0; j < k; j++ )
        for ( i = 0; i < n; i++ )
            c[i+n*j] -= lambda[j] * x[i+n*j];

    value = r8mat_norm_fro ( n, k, c );
    free ( c );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1dd_apply ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DD_APPLY applies the 1D DD Laplacian to a vector.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
    at both ends of [0,6] is applied to a vector of 7 values, with a spacing
    of H = 6/(N+1) = 1 at the points X:
      0  1  2  3  4  5  6
    and has the matrix form L:
       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Input, double U[N], the value at each point.
    Output, double L1DD_APPLY[N], the Laplacian evaluated at each point.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * u = s_data->a1;
	const register ityp h = s_data->a2;
	
    dim_typ i;
    ityp *lu;

    if ( n < 3 )
        return NULL;

    lu = ( double * ) malloc ( n * sizeof ( double ) );

    i = 0;
    lu[i] = ( 2.00 * u[i] - u[i+1] ) / h / h;
    for ( i = 1; i < n - 1; ++i )
        lu[i] = ( - u[i-1] + 2.0 * u[i] - u[i+1] ) / h / h;
    i = n - 1;
    lu[i] = ( - u[i-1] + 2.00 * u[i] ) / h / h;
    return lu;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1dd_cholesky ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DD_CHOLESKY computes the Cholesky factor of the 1D DD Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1DD_CHOLESKY[N*N], the Cholesky factor.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    ityp *c;
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    c = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            c[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
        c[i+i*n] = sqrt ( ( ityp ) ( i + 2 ) ) / sqrt ( ( ityp ) ( i + 1 ) );

    for ( i = 0; i < n - 1; ++i )
        c[i+(i+1)*n] = - sqrt ( ( ityp ) ( i + 1 ) )/ sqrt ( ( ityp ) ( i + 2 ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            c[i+j*n] /= h;

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l1dd_eigen ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DD_EIGEN returns eigen information for the 1D DD Laplacian.
  Discussion:
    The grid points are assumed to be evenly spaced by H.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, double H, the spacing between points.
    Output, double V[N*N], the eigenvectors.
    Output, double LAMBDA[N], the eigenvalues.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * v = s_data->a2;
	ityp * lambda = s_data->a3; 
	
    dim_typ i, j;
    ityp i_r8;
    ityp j_r8;
    ityp n_r8;
    ityp s;
    ityp theta;

    n_r8 = ( ityp ) ( n );

    for ( j = 0; j < n; ++j )
    {
        j_r8 = ( ityp ) ( j + 1 );
        theta = 0.50 * M_PI * j_r8 / ( n_r8 + 1.00 );
        lambda[j] = pow ( 2.00 * sin ( theta ) / h, 2 );
        for ( i = 0; i < n; ++i )
        {
                i_r8 = ( ityp ) ( i + 1 );
                theta = M_PI * i_r8 * j_r8 / ( n_r8 + 1.00 );
                v[i+j*n] = sqrt ( 2.00 / ( n_r8 + 1.00 ) ) * sin ( theta );
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1dd ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DD stores the 1D DD Laplacian as a full matrix.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
    at both ends of [0,6] has the matrix form L:
       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1DD[N*N], the Laplacian matrix.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    dim_typ i, j;
    ityp *l;

    if ( n < 3 )
        return NULL;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    i = 0;
    l[i+i*n]     =  2.00 / h / h;
    l[i+(i+1)*n] = -1.00 / h / h;

    for ( i = 1; i < n - 1; ++i )
    {
        l[i+(i-1)*n] = -1.00 / h / h;
        l[i+i*n] =      2.00 / h / h;
        l[i+(i+1)*n] = -1.00 / h / h;
    }

    i = n - 1;
    l[i+(i-1)*n] = -1.00 / h / h;
    l[i+i*n] =      2.00 / h / h;

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1dd_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DD_INVERSE stores the inverse of the 1D DD Laplacian.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
    at both ends of [0,6] has the matrix form L:
       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1DD_INVERSE[N*N], the inverse of the Laplacian matrix.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    dim_typ i, j;
    ityp *l;

    if ( n < 3 )
        return NULL;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = ( ityp ) ( MIN ( i + 1, j + 1 ) * ( n - MAX ( i, j ) ) )* h * h / ( ityp ) ( n + 1 );

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _l1dd_lu ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DD_LU computes the LU factors of the 1D DD Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L[N*N], U[N*N], the LU factors.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * l = s_data->a2;
	ityp * u = s_data->a3; 
	
    ityp i_r8;
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;


    for ( i = 0; i < n; ++i )
        l[i+i*n] = 1.00;

    for ( i = 1; i < n; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        l[i+(i-1)*n] = - ( i_r8 - 1.00 ) / i_r8;
    }

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] /= h;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            u[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        u[i+i*n] = ( i_r8 + 1.00 ) / i_r8;
    }

    for ( i = 0; i < n - 1; ++i )
        u[i+(i+1)*n] = - 1.00;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i)
            u[i+j*n] /= h;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1dn_apply ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DN_APPLY applies the 1D DN Laplacian to a vector.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with left Dirichlet and right
    Neumann condition on [0,6] has the matrix form L:
       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Input, double U[N], the value at each point.
    Output, double L1DN_APPLY[N], the Laplacian evaluated at each point.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * u = s_data->a1;
	const register ityp h = s_data->a2;
	
    dim_typ i;
    ityp *lu;

    if ( n < 3 )
        return NULL;

    lu = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    i = 0;
    lu[i] = ( 2.00 * u[i] - u[i+1] ) / h / h;
    for ( i = 1; i < n - 1; ++i )
        lu[i] = ( - u[i-1] + 2.00 * u[i] - u[i+1] ) / h / h;
    i = n - 1;
    lu[i] = ( - u[i-1] + u[i] ) / h / h;

    return lu;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1dn_cholesky ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DN_CHOLESKY computes the Cholesky factor of the 1D DN Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1DN_CHOLESKY[N*N], the Cholesky factor.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    ityp *c;
    ityp i_r8;
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    c = ( double * ) malloc ( n * n * sizeof ( double ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            c[i+j*n] = 0.00;

    for ( i = 0; i < n - 1; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        c[i+i*n]   =   sqrt ( i_r8 + 1.00 ) / sqrt ( i_r8 );
        c[i+(i+1)*n] = - sqrt ( i_r8 ) / sqrt ( i_r8 + 1.00 );
    }
    c[n-1+(n-1)*n] = 1.00 / sqrt ( ( ityp ) ( n ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            c[i+j*n] /= h;

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _l1dn_eigen ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DN_EIGEN returns eigeninformation for the 1D DN Laplacian.
  Discussion:
    The grid points are assumed to be evenly spaced by H.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, double H, the spacing between points.
    Output, double V[N*N], the eigenvectors.
    Output, double LAMBDA[N], the eigenvalues.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * v = s_data->a2;
	ityp * lambda = s_data->a3; 
	
    dim_typ i, j;
    ityp i_r8;
    ityp j_r8;
    ityp n_r8;
    ityp s;
    ityp theta;

    n_r8 = ( ityp ) ( n );

    for ( j = 0; j < n; ++j )
    {
        j_r8 = ( ityp ) ( j + 1 );
        theta = M_PI * ( j_r8 - 0.50 ) / ( 2.00 * n_r8 + 1.00 );
        lambda[j] = pow ( 2.00 * sin ( theta ) / h, 2 );
        for ( i = 0; i < n; ++i )
        {
            i_r8 = ( ityp ) ( i + 1 );
            theta = M_PI * i_r8 * ( 2.00 * j_r8 - 1.00 ) /( 2.00 * n_r8 + 1.00 );
            v[i+j*n] = sqrt ( 2.00 / ( n_r8 + 0.50 ) ) * sin ( theta );
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1dn ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DN stores the 1D DN Laplacian as a full matrix.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with left Dirichlet and right
    Neumann condition on [0,6] has the matrix form L:
       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1DN[N*N], the Laplacian matrix.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    dim_typ i, j;
    ityp *l;

    if ( n < 3 )
        return NULL;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    i = 0;
    l[i+i*n]     =  2.00 / h / h;
    l[i+(i+1)*n] = -1.00 / h / h;

    for ( i = 1; i < n - 1; ++i )
    {
        l[i+(i-1)*n] = -1.00 / h / h;
        l[i+i*n] =      2.00 / h / h;
        l[i+(i+1)*n] = -1.00 / h / h;
    }

    i = n - 1;
    l[i+(i-1)*n] = -1.00 / h / h;
    l[i+i*n] =      1.00 / h / h;

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1dn_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DN_INVERSE stores the inverse of the 1D DN Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1DN_INVERSE[N*N], the inverse of the Laplacian matrix.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    dim_typ i, j;
    ityp *l;

    if ( n < 3 )
        return NULL;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = ( ityp ) ( MIN ( i + 1, j + 1 ) ) * h * h;

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l1dn_lu ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1DD_LU computes the LU factors of the 1D DN Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L[N*N], U[N*N], the LU factors.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * l = s_data->a2;
	ityp * u = s_data->a3; 
	
    ityp i_r8;
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
        l[i+i*n] = 1.00;

    for ( i = 1; i < n; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        l[i+(i-1)*n] = - ( i_r8 - 1.00 ) / i_r8;
    }

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i)
            l[i+j*n] /= h;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            u[i+j*n] = 0.00;

    for ( i = 0; i < n - 1; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        u[i+i*n] = ( i_r8 + 1.00 ) / i_r8;
    }
    i = n - 1;
    i_r8 = ( ityp ) ( i + 1 );
    u[i+i*n] = 1.00 / i_r8;

    for ( i = 0; i < n - 1; ++i )
        u[i+(i+1)*n] = - 1.00;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            u[i+j*n] /= h;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1nd_apply ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1ND_APPLY applies the 1D ND Laplacian to a vector.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with left Neumann and right Dirichlet
    boundary conditions on [0,6] has the matrix form L:
       1 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Input, double U[N], the value at each point.
    Output, double L1ND_APPLY[N], the Laplacian evaluated at each point.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * u = s_data->a1;
	const register ityp h = s_data->a2;
	
    dim_typ i;
    ityp *lu;

    if ( n < 3 )
        return NULL;

    lu = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    i = 0;
    lu[i] = ( u[i] - u[i+1] ) / h / h;
    for ( i = 1; i < n - 1; ++i)
        lu[i] = ( - u[i-1] + 2.00 * u[i] - u[i+1] ) / h / h;
    i = n - 1;
    lu[i] = ( - u[i-1] + 2.00 * u[i] ) / h / h;
    return lu;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1nd_cholesky ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1ND_CHOLESKY computes the Cholesky factor of the 1D ND Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1ND_CHOLESKY[N*N], the Cholesky factor.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    ityp *c;
    ityp i_r8;
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    c = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            c[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
        c[i+i*n] = 1.00;

    for ( i = 0; i < n - 1; ++i )
        c[i+(i+1)*n] = - 1.00;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i)
            c[i+j*n] /= h;

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l1nd_eigen ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1ND_EIGEN returns eigeninformation for the 1D ND Laplacian.
  Discussion:
    The grid points are assumed to be evenly spaced by H.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, double H, the spacing between points.
    Output, double V[N*N], the eigenvectors.
    Output, double LAMBDA[N], the eigenvalues.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * v = s_data->a2;
	ityp * lambda = s_data->a3; 
	
    dim_typ i, j;
    ityp i_r8;
    ityp j_r8;
    ityp n_r8;
    ityp s;
    ityp theta;

    n_r8 = ( ityp ) ( n );

    for ( j = 0; j < n; ++j )
    {
        j_r8 = ( ityp ) ( j + 1 );
        theta = M_PI * ( j_r8 - 0.50 ) / ( 2.00 * n_r8 + 1.0 );
        lambda[j] = pow ( 2.00 * sin ( theta ) / h, 2 );
        for ( i = 0; i < n; ++i )
        {
            i_r8 = ( ityp ) ( i + 1 );
            theta = M_PI * ( i_r8 - 0.50 ) * ( 2.00 * j_r8 - 1.00 ) /( 2.0 * n_r8 + 1.00 );
            v[i+j*n] = sqrt ( 2.00 / ( n_r8 + 0.50 ) ) * cos ( theta );
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1nd ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1ND stores the 1D ND Laplacian as a full matrix.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with left Neumann and right Dirichlet
    boundary conditions on [0,6] has the matrix form L:
       1 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1ND[N*N], the Laplacian matrix.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    dim_typ i, j;
    ityp *l;

    if ( n < 3 )
        return NULL;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    i = 0;
    l[i+i*n]     =  1.00 / h / h;
    l[i+(i+1)*n] = -1.00 / h / h;

    for ( i = 1; i < n - 1; ++i )
    {
        l[i+(i-1)*n] = -1.00 / h / h;
        l[i+i*n] =      2.00 / h / h;
        l[i+(i+1)*n] = -1.00 / h / h;
    }

    i = n - 1;
    l[i+(i-1)*n] = -1.00 / h / h;
    l[i+i*n] =      2.00 / h / h;

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1nd_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1ND_INVERSE stores the inverse of the 1D ND Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1ND_INVERSE[N*N], the inverse of the Laplacian matrix.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    dim_typ i, j;
    ityp *l;

    if ( n < 3 )
        return NULL;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = ( ityp ) ( n - MAX ( i, j ) ) * h * h;

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l1nd_lu ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1ND_LU computes the LU factors of the 1D ND Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L[N*N], U[N*N], the LU factors.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * l = s_data->a2;
	ityp * u = s_data->a3; 
	
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
        l[i+i*n] = 1.00;

    for ( i = 1; i < n; ++i )
        l[i+(i-1)*n] = - 1.00;

    for ( j = 0; j < n; j++ )
        for ( i = 0; i < n; i++ )
            l[i+j*n] /= h;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            u[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
        u[i+i*n] = 1.00;

    for ( i = 0; i < n - 1; ++i )
        u[i+(i+1)*n] = - 1.00;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            u[i+j*n] /= h;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1nn_apply ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1NN_APPLY applies the 1D NN Laplacian to a vector.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with left Neumann and right Neumann
    boundary conditions on [0,6] has the matrix form L:
       1 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Input, double U[N], the value at each point.
    Output, double L1NN_APPLY[N], the Laplacian evaluated at each point.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * u = s_data->a1;
	const register ityp h = s_data->a2;
	
	
    dim_typ i;
    ityp *lu;

    if ( n < 3 )
        return NULL;

    lu = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    i = 0;
    lu[i] = ( u[i] - u[i+1] ) / h / h;
    for ( i = 1; i < n - 1; ++i )
        lu[i] = ( - u[i-1] + 2.00 * u[i] - u[i+1] ) / h / h;
    i = n - 1;
    lu[i] = ( - u[i-1] +  u[i] ) / h / h;
    return lu;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1nn_cholesky ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1NN_CHOLESKY computes the Cholesky factor of the 1D NN Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1NN_CHOLESKY[N*N], the Cholesky factor.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    ityp *c;
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    c = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            c[i+j*n] = 0.00;

    for ( i = 0; i < n - 1; ++i )
    {
        c[i+i*n]     = + 1.00;
        c[i+(i+1)*n] = - 1.00;
    }

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < n; ++i )
            c[i+j*n] /= h;


    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l1nn_eigen ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1NN_EIGEN returns eigeninformation for the 1D NN Laplacian.
  Discussion:
    The grid points are assumed to be evenly spaced by H.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, double H, the spacing between points.
    Output, double V[N*N], the eigenvectors.
    Output, double LAMBDA[N], the eigenvalues.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * v = s_data->a2;
	ityp * lambda = s_data->a3; 
	
    dim_typ i, j;
    ityp i_r8;
    ityp j_r8;
    ityp n_r8;
    ityp s;
    ityp theta;

    n_r8 = ( ityp ) ( n );

    for ( j = 0; j < n; ++j )
    {
        j_r8 = ( ityp ) ( j + 1 );
        theta = M_PI * ( j_r8 - 1.00 ) / ( 2.00 * n_r8 );
        lambda[j] = pow ( 2.00 * sin ( theta ) / h, 2 );
        if ( j == 0 )
        {
            for ( i = 0; i < n; ++i )
            v[i+j*n] = sqrt ( n_r8 );
        }
        else
        {
            for ( i = 0; i < n; ++i)
            {
                i_r8 = ( ityp ) ( i + 1 );
                theta = M_PI * ( i_r8 - 0.50 ) * ( j_r8 - 1.00 ) / n_r8;
                v[i+j*n] = sqrt ( 2.00 / n_r8 ) * cos ( theta );
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1nn ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1NN stores the 1D NN Laplacian as a full matrix.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with Neumann boundary conditions
    at both ends of [0,6] has the matrix form L:
       1 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1NN[N*N], the Laplacian matrix.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    dim_typ i, j;
    ityp *l;

    if ( n < 3 )
        return NULL;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;


    i = 0;
    l[i+i*n]     =  1.00 / h / h;
    l[i+(i+1)*n] = -1.00 / h / h;

    for ( i = 1; i < n - 1; ++i )
    {
        l[i+(i-1)*n] = -1.00 / h / h;
        l[i+i*n] =      2.00 / h / h;
        l[i+(i+1)*n] = -1.00 / h / h;
    }

    i = n - 1;
    l[i+(i-1)*n] = -1.00 / h / h;
    l[i+i*n] =      1.00 / h / h;

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _l1nn_lu ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1NN_LU computes the LU factors of the 1D NN Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L[N*N], U[N*N], the LU factors.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * l = s_data->a2;
	ityp * u = s_data->a3; 
	
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i)
        l[i+i*n] = 1.00;

    for ( i = 1; i < n; ++i )
        l[i+(i-1)*n] = - 1.00;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] /= h;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            u[i+j*n] = 0.00;

    for ( i = 0; i < n - 1; ++i )
        u[i+i*n] = 1.00;

    i = n - 1;
    u[i+i*n] = 0.00;

    for ( i = 0; i < n - 1; ++i )
        u[i+(i+1)*n] = - 1.00;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            u[i+j*n] /= h;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1pp_apply ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1PP_APPLY applies the 1D PP Laplacian to a vector.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with periodic boundary conditions
    on [0,6] has the matrix form L:
       2 -1  0  0 -1
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
      -1  0  0 -1  2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Input, double U[N], the value at each point.
    Output, double L1PP_APPLY[N], the Laplacian evaluated at each point.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * u = s_data->a1;
	const register ityp h = s_data->a2;
	
    dim_typ i;
    ityp *lu;

    if ( n < 3 )
        return NULL;

    lu = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    i = 0;
    lu[i] = ( - u[n-1] + 2.00 * u[i] - u[i+1] ) / h / h;
    for ( i = 1; i < n - 1; ++i )
        lu[i] = ( - u[i-1] + 2.00 * u[i] - u[i+1] ) / h / h;
    i = n - 1;
    lu[i] = ( - u[i-1] + 2.00 * u[i] - u[0] ) / h / h;

    return lu;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1pp_cholesky ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1PP_CHOLESKY computes the Cholesky factor of the 1D PP Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1PP_CHOLESKY[N*N], the Cholesky factor.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    ityp *c;
    ityp i_r8;
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    c = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < n; ++i )
            c[i+j*n] = 0.00;

    for ( i = 0; i < n - 1; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        c[i+i*n] = sqrt ( i_r8 + 1.00 ) / sqrt ( i_r8 );
    }

    for ( i = 0; i < n - 2; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        c[i+(i+1)*n] = - i_r8 / ( i_r8 + 1.00 ) * sqrt ( i_r8 + 1.00 )/ sqrt ( i_r8 );
    }

    for ( i = 0; i < n - 2; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        c[i+(n-1)*n] = - 1.00 / ( i_r8 + 1.00 ) * sqrt ( i_r8 + 1.00 )/ sqrt ( i_r8 );
    }

    i = n - 2;
    i_r8 = ( ityp ) ( i + 1 );
    c[i+(n-1)*n] = - ( ityp ) ( n ) / ( i_r8 + 1.00 )* sqrt ( i_r8 + 1.00 ) / sqrt ( i_r8 );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            c[i+j*n] /= h;

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l1pp_eigen ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1PP_EIGEN returns eigeninformation for the 1D PP Laplacian.
  Discussion:
    The grid points are assumed to be evenly spaced by H.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, double H, the spacing between points.
    Output, double V[N*N], the eigenvectors.
    Output, double LAMBDA[N], the eigenvalues.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * v = s_data->a2;
	ityp * lambda = s_data->a3; 
	
    dim_typ i, j;
    ityp i_r8;
    ityp j_r8;
    ityp n_r8;
    ityp s;
    ityp theta;

    n_r8 = ( ityp ) ( n );

    for ( j = 0; j < n; ++j)
    {
        j_r8 = ( ityp ) ( j + 1 );
        theta = M_PI * ((j%2) == 0)?( j_r8 - 1.0 ) / ( 2.0 * n_r8 ) : j_r8         / ( 2.0 * n_r8 );
        lambda[j] = pow ( 2.00 * sin ( theta ) / h, 2 );

        if ( ( j % 2 ) == 0 )
        {
            if ( j == 0 )
                for ( i = 0; i < n; ++i )
                    v[i+j*n] = 1.00 / sqrt ( n_r8 );
            else
            {
                for ( i = 0; i < n; i++ )
                {
                    i_r8 = ( double ) ( i + 1 );
                    theta = M_PI * ( i_r8 - 0.5 ) * ( j_r8 - 1.0 ) /  n_r8;
                    v[i+j*n] = sqrt ( 2.0 / n_r8 ) * cos ( theta );
                }
            }
        }
        else
        {
            if ( j == n - 1 )
            {
                s = - 1.0 / sqrt ( n_r8 );
                for ( i = 0; i < n; ++i )
                    v[i+j*n] = s = -s;
            }
            else
            {
                for ( i = 0; i < n; ++i )
                {
                    i_r8 = ( ityp ) ( i + 1 );
                    theta = M_PI * ( i_r8 - 0.50 ) * j_r8 / n_r8;
                    v[i+j*n] = sqrt ( 2.00 / n_r8 ) * sin ( theta );
                }
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _l1pp ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1PP stores the 1D PP Laplacian as a full matrix.
  Discussion:
    The N grid points are assumed to be evenly spaced by H.
    For N = 5, the discrete Laplacian with periodic boundary conditions
    has the matrix form L:
       2 -1  0  0 -1
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
      -1  0  0 -1  2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified
    30 October 2013.
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L1PP[N*N], the Laplacian matrix.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    dim_typ i, j;
    ityp *l;

    if ( n < 3 )
        return NULL;

    l = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            l[i+j*n] = 0.00;

    i = 0;
    l[i+i*n]       =  2.00 / h / h;
    l[i+(i+1)*n]   = -1.00 / h / h;
    l[i+(n-1)*n]   = -1.00 / h / h;

    for ( i = 1; i < n - 1; ++i )
    {
        l[i+(i-1)*n] = -1.00 / h / h;
        l[i+i*n] =      2.00 / h / h;
        l[i+(i+1)*n] = -1.00 / h / h;
    }

    i = n - 1;
    l[i+0*n] =     -1.00 / h / h;
    l[i+(i-1)*n] = -1.00 / h / h;
    l[i+i*n] =      2.00 / h / h;

    return l;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l1pp_lu ( void * data)
/******************************************************************************/
/*
  Discussion:
    L1PP_LU computes the LU factors of the 1D PP Laplacian.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    N must be at least 3.
    Input, double H, the spacing between points.
    Output, double L[N*N], U[N*N], the LU factors.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	ityp * l = s_data->a2;
	ityp * u = s_data->a3; 
	
    ityp i_r8;
    dim_typ i, j;

    if ( n < 3 )
        return NULL;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i)
            l[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
        l[i+i*n] = 1.00;

    for ( i = 1; i < n - 1; ++i)
    {
        i_r8 = ( ityp ) ( i + 1 );
        l[i+(i-1)*n] =   - ( i_r8 - 1.00 ) / i_r8;
        l[n-1+(i-1)*n] =          - 1.00   / i_r8;
    }
    i = n - 1;
    l[i+(i-1)*n] = -1.00;

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < n; ++i )
            l[i+j*n] /= h;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            u[i+j*n] = 0.00;

    for ( i = 0; i < n - 2; ++i )
    {
        i_r8 = ( ityp ) ( i + 1 );
        u[i+i*n] = ( i_r8 + 1.0 ) / i_r8;
        u[i+(i+1)*n] = - 1.00;
        u[i+(n-1)*n] =   - 1.00 / i_r8;
    }

    i = n - 2;
    i_r8 = ( ityp ) ( i + 1 );
    u[i+i*n] = ( i_r8 + 1.00 ) / i_r8;
    u[i+(i+1)*n] = - ( i_r8 + 1.00 ) / i_r8;

    i = n - 1;
    u[i+i*n] = 0.00;

    for ( j = 0; j < n; j++ )
        for ( i = 0; i < n; ++i )
            u[i+j*n] /= h;

    return NULL;
}

#endif
