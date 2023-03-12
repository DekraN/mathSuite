#ifndef __DISABLEDEEP_JACOBI

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dif2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIF2 returns the DIF2 matrix.
  Example:
    N = 5
    2 -1  .  .  .
   -1  2 -1  .  .
    . -1  2 -1  .
    .  . -1  2 -1
    .  .  . -1  2
  Properties:
    A is banded, with bandwidth 3.
    A is tridiagonal.
    Because A is tridiagonal, it has property A (bipartite).
    A is a special case of the TRIS or tridiagonal scalar matrix.
    A is integral, therefore det ( A ) is integral, and
    det ( A ) * inverse ( A ) is integral.
    A is Toeplitz: constant along diagonals.
    A is symmetric: A' = A.
    Because A is symmetric, it is normal.
    Because A is normal, it is diagonalizable.
    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
    A is positive definite.
    A is an M matrix.
    A is weakly diagonally dominant, but not strictly diagonally dominant.
    A has an LU factorization A = L * U, without pivoting.
      The matrix L is lower bidiagonal with subdiagonal elements:
        L(I+1,I) = -I/(I+1)
      The matrix U is upper bidiagonal, with diagonal elements
        U(I,I) = (I+1)/I
      and superdiagonal elements which are all -1.
    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
      L(I,I) =    sqrt ( (I+1) / I )
      L(I,I-1) = -sqrt ( (I-1) / I )
    The eigenvalues are
      LAMBDA(I) = 2 + 2 * COS(I*M_PI/(N+1))
                = 4 SIN^2(I*M_PI/(2*N+2))
    The corresponding eigenvector X(I) has entries
       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*M_PI/(N+1) ).
    Simple linear systems:
      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
    det ( A ) = N + 1.
    The value of the determinant can be seen by induction,
    and expanding the determinant across the first row:
      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
                = 2 * N - (N-1)
                = N + 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 September 2010
  Author:
    John Burkardt
  Reference:
    Robert Gregory, David Karney,
    A Collection of Matrices for Testing Computational Algorithms,
    Wiley, 1969,
    ISBN: 0882756494,
    LC: QA263.68
    Morris Newman, John Todd,
    Example A8,
    The evaluation of matrix inversion programs,
    Journal of the Society for Industrial and Applied Mathematics,
    Volume 6, Number 4, pages 466-476, 1958.
    John Todd,
    Basic Numerical Mathematics,
    Volume 2: Numerical Algebra,
    Birkhauser, 1980,
    ISBN: 0817608117,
    LC: QA297.T58.
    Joan Westlake,
    A Handbook of Numerical Matrix Inversion and Solution of
    Linear Equations,
    John Wiley, 1968,
    ISBN13: 978-0471936756,
    LC: QA263.W47.
  Parameters:
    Input, int M, N, the order of the matrix.
    Output, double DIF2[M*N], the matrix.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
    dim_typ i, j;
    ityp *a = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            if ( j == i - 1 )
                a[i+j*m] = -1.00;
            else if ( j == i )
                a[i+j*m] = 2.00;
            else
                a[i+j*m] = 0.00 - (j == i + 1);

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_residual_norm ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_RESIDUAL_NORM returns the norm of A*x-b.
  Discussion:
    A is an MxN R8MAT, a matrix of R8's.
    X is an N R8VEC, and B is an M R8VEC.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 June 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of the matrix.
    Input, double A[M,N], the M by N matrix.
    Input, double X[N], the vector to be multiplied by A.
    Input, double B[M], the right hand side vector.
    Output, double R8MAT_RESIDUAL_NORM, the norm of A*x-b.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dt3pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * x = s_data->a3;
	ityp * b = s_data->a4;
	
    dim_typ i, j;
    ityp r_norm;
    ityp * r = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
    {
        r[i] = - b[i];
        for ( j = 0; j < n; ++j )
            r[i] += a[i+j*m] * x[j];
    }

    r_norm = 0.00;
    for ( i = 0; i < m; ++i )
        r_norm += r[i] * r[i];
    free ( r );
    
    result = sqrt ( r_norm );
    return &result;
}

#endif
