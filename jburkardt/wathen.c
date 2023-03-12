#ifndef __DISABLEDEEP_WATHEN

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _mv_st ( void * data)
/******************************************************************************/
/*
  Purpose:
    MV_ST multiplies a sparse triple matrix times a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int NZ_NUM, the number of nonzero values.
    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and 
    column indices.
    Input, double A[NZ_NUM], the nonzero values in the matrix.
    Input, double X[N], the vector to be multiplied.
    Output, double MV_ST[M], the product A*X.
*/
{
	const _3dt2pi2pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register dim_typ nz_num = s_data->a2;
	int * row = s_data->a3;
	int * col = s_data->a4;
	ityp * a = s_data->a5;
	ityp * x = s_data->a6;
	
	ityp *b;
	dim_typ i, k; 
	
	b = ( ityp * ) malloc ( m * sizeof ( ityp ) );
	
	for ( i = 0; i < m; ++i )
		b[i] = 0.00;
	
	for ( k = 0; k < nz_num; ++k )
		b[row[k]] += a[k] * x[col[k]];
	
	return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _mv_gb ( void * data)
/******************************************************************************/
/*
  Purpose:
    MV_GB multiplies a banded matrix by an R8VEC.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    06 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of the matrix.
    M must be positive.
    Input, int N, the number of columns of the matrix.
    N must be positive.
    Input, int ML, MU, the lower and upper bandwidths.
    Input, double A[(2*ML+MU+1)*N], the matrix.
    Input, double X[N], the vector to be multiplied by A.
    Output, double MV_GB[M], the product A * x.
*/
{
	const _4dt2pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register dim_typ ml = s_data->a2;
	const register dim_typ mu = s_data->a3;
	ityp * a = s_data->a4;
	ityp * x = s_data->a5;
	
	ityp *b;
	dim_typ i, j;
	int jhi;
	int jlo;
	
	b = ( ityp * ) malloc ( m * sizeof ( ityp ) );
	
	for ( i = 0; i < m; ++i )
		b[i] = 0.00;
	
	for ( i = 0; i < m; ++i )
	{
		jlo = MAX ( 0, i - ml );
		jhi = MIN ( n - 1, i + mu );
		for ( j = jlo; j <= jhi; ++j )
			b[i] += a[i-j+ml+mu+j*((ml<<1)+mu+1)] * x[j];
	}
	
	return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _mv_ge ( void * data)
/******************************************************************************/
/*
  Purpose:
    MV_GE multiplies a GE matrix by an R8VEC.
  Discussion:
    The GE storage format is used for a general M by N matrix.  A storage 
    space is made for each entry.  The two dimensional logical
    array can be thought of as a vector of M*N entries, starting with
    the M entries in the column 1, then the M entries in column 2
    and so on.  Considered as a vector, the entry A(I,J) is then stored
    in vector location I+(J-1)*M.
    GE storage is used by LINPACK and LAPACK.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    06 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of the matrix.
    M must be positive.
    Input, int N, the number of columns of the matrix.
    N must be positive.
    Input, double A(M,N), the matrix.
    Input, double X[N], the vector to be multiplied by A.
    Output, double MV_GE[M], the product A * x.
*/
{
	const _2dt2pit * const s_data = data; 
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * a = s_data->a2;
	ityp * x = s_data->a3;
	
	ityp *b;
	dim_typ i, j;
	
	b = ( ityp * ) malloc ( m * sizeof ( ityp ) );
	
	for ( i = 0; i < m; ++i )
	{
		b[i] = 0.00;
		for ( j = 0; j < n; ++j)
			b[i] +=a[i+j*m] * x[j];
	}

  return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cg_st ( void * data)
/******************************************************************************/
/*
  Purpose:
    CG_ST uses the conjugate gradient method for a sparse triplet (ST) matrix.
  Discussion:
    The linear system has the form A*x=b, where A is a positive-definite
    symmetric matrix, stored as a full storage matrix.
    The method is designed to reach the solution to the linear system
      A * x = b
    after N computational steps.  However, roundoff may introduce
    unacceptably large errors for some problems.  In such a case,
    calling the routine a second time, using the current solution estimate
    as the new starting guess, should result in improved results.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 June 2014
  Author:
    John Burkardt
  Reference:
    Frank Beckman,
    The Solution of Linear Equations by the Conjugate Gradient Method,
    in Mathematical Methods for Digital Computers,
    edited by John Ralston, Herbert Wilf,
    Wiley, 1967,
    ISBN: 0471706892,
    LC: QA76.5.R3.
  Parameters:
    Input, int N, the order of the matrix.
    Input, int NZ_NUM, the number of nonzeros.
    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column
    indices of the nonzero entries.
    Input, double A[NZ_NUM], the nonzero entries.
    Input, double B[N], the right hand side vector.
    Input/output, double X[N].
    On input, an estimate for the solution, which may be 0.
    On output, the approximate solution vector.
*/
{
	const _2dt2pi3pit * const s_data = data;
	const register dim_typ n = s_data->a0; 
	const register dim_typ nz_num = s_data->a1;
	int * row = s_data->a2;
	int * col = s_data->a3;
	ityp * a = s_data->a4;
	ityp * b = s_data->a5;
	ityp * x = s_data->a6;
	
    ityp alpha;
    ityp *ap;
    ityp beta;
    dim_typ i;
    dim_typ it;
    ityp *p;
    ityp pap;
    ityp pr;
    ityp *r;
    ityp rap;
    /*
    Initialize
    AP = A * x,
    R  = b - A * x,
    P  = b - A * x.
    */
    ap = mv_st ( n, n, nz_num, row, col, a, x );

    r = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    p = ( ityp *) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i)
    {
        r[i] = b[i] - ap[i];
        p[i] = b[i] - ap[i];
    }

    /*
    Do the N steps of the conjugate gradient method.
    */
    for ( it = 1; it <= n; ++it )
    {
        /*
        Compute the matrix*vector product AP = A*P.
        */
        free ( ap );
        ap = mv_st ( n, n, nz_num, row, col, a, p );
        /*
        Compute the dot products
        PAP = P*AP,
        PR  = P*R
        Set
        ALPHA = PR / PAP.
        */
        pap = r8vec_dot_product ( n, p, ap );
        pr = r8vec_dot_product ( n, p, r );

        if ( pap == 0.00 )
            break;

        alpha = pr / pap;
        /*
        Set
        X = X + ALPHA * P
        R = R - ALPHA * AP.
        */
        for ( i = 0; i < n; ++i )
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        /*
        Compute the vector dot product
        RAP = R*AP
        Set
        BETA = - RAP / PAP.
        */
        rap = r8vec_dot_product ( n, r, ap );

        beta = - rap / pap;
        /*
        Update the perturbation vector
        P = R + BETA * P.
        */
        for ( i = 0; i < n; ++i)
            p[i] = r[i] + beta * p[i];
    }

    free ( ap );
    free ( p );
    free ( r );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wathen_bandwidth ( void * data)
/******************************************************************************/
/*
  Purpose:
    WATHEN_BANDWIDTH returns the bandwidth of the WATHEN matrix.
  Discussion:
    The bandwidth measures the minimal number of contiguous diagonals,
    including the central diagonal, which contain all the nonzero elements
    of a matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 June 2014
  Author:
    John Burkardt
  Reference:
    Nicholas Higham,
    Algorithm 694: A Collection of Test Matrices in MATLAB,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 3, September 1991, pages 289-305.
    Andrew Wathen,
    Realistic eigenvalue bounds for the Galerkin mass matrix,
    IMA Journal of Numerical Analysis,
    Volume 7, 1987, pages 449-457.
  Parameters:
    Input, int NX, NY, values which determine the size of A.
    Output, int *L, *D, *U, the lower, diagonal, and upper
    bandwidths of the matrix,
*/
{
	const _2dt3pdt * const s_data = data;
	const register dim_typ nx = s_data->a0; 
	const register dim_typ ny = s_data->a1;
	dim_typ * l = s_data->a2;
	dim_typ * d = s_data->a3;
	dim_typ * u = s_data->a4;
	
    *l = 3 * nx + 4;
    *d = 1;
    *u = 3 * nx + 4;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wathen_gb ( void * data)
/******************************************************************************/
/*
  Purpose:
    WATHEN_GB returns the Wathen matrix, using general banded (GB) storage.
  Discussion:
    The Wathen matrix is a finite element matrix which is sparse.
    The entries of the matrix depend in part on a physical quantity
    related to density.  That density is here assigned random values between
    0 and 100.
    The matrix order N is determined by the input quantities NX and NY,
    which would usually be the number of elements in the X and Y directions.
    The value of N is
      N = 3*NX*NY + 2*NX + 2*NY + 1,
    The matrix is the consistent mass matrix for a regular NX by NY grid
    of 8 node serendipity elements.
    The local element numbering is
      3--2--1
      |     |
      4     8
      |     |
      5--6--7
    Here is an illustration for NX = 3, NY = 2:
     23-24-25-26-27-28-29
      |     |     |     |
     19    20    21    22
      |     |     |     |
     12-13-14-15-16-17-18
      |     |     |     |
      8     9    10    11
      |     |     |     |
      1--2--3--4--5--6--7
    For this example, the total number of nodes is, as expected,
      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
    The matrix is symmetric positive definite for any positive values of the
    density RHO(X,Y).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 June 2014
  Author:
    John Burkardt
  Reference:
    Nicholas Higham,
    Algorithm 694: A Collection of Test Matrices in MATLAB,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 3, September 1991, pages 289-305.
    Andrew Wathen,
    Realistic eigenvalue bounds for the Galerkin mass matrix,
    IMA Journal of Numerical Analysis,
    Volume 7, Number 4, October 1987, pages 449-457.
  Parameters:
    Input, int NX, NY, values which determine the size
    of the matrix.
    Input, int N, the number of rows and columns.
    Input/output, int *SEED, the random number seed.
    Output, double WATHEN_GB[(9*NX+13)*N], the matrix.
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ nx = s_data->a0;
	const register dim_typ ny = s_data->a1;
	const register dim_typ n = s_data->a2;
	int * seed = s_data->a3;
	
    ityp *a;
    const ityp em[64] =
    {
        6.00, -6.00,  2.00, -8.00,  3.00, -8.00,  2.00, -6.00,
        -6.00, 32.00, -6.00, 20.00, -8.00, 16.00, -8.00, 20.00,
        2.00, -6.00,  6.00, -6.00,  2.00, -8.00,  3.00, -8.00,
        -8.00, 20.00, -6.00, 32.00, -6.00, 20.00, -8.00, 16.00,
        3.00, -8.00,  2.00, -6.00,  6.00, -6.00,  2.00, -8.00,
        -8.00, 16.00, -8.00, 20.00, -6.00, 32.00, -6.00, 20.00,
        2.00, -8.00,  3.00, -8.00,  2.00, -6.00,  6.00, -6.00,
        -6.00, 20.00, -8.00, 16.00, -8.00, 20.00, -6.00, 32.00
    };
    dim_typ i;
    dim_typ ii;
    dim_typ j;
    dim_typ jj;
    dim_typ kcol;
    dim_typ krow;
    dim_typ lda;
    dim_typ ml;
    dim_typ mu;
    int node[8];
    ityp rho;

    ml = 3 * nx + 4;
    mu = 3 * nx + 4;
    lda = (ml<<1) + mu + 1;
    a = ( ityp * ) malloc ( lda * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < lda; ++i )
            a[i+j*lda] = 0.00;

    for ( j = 0; j < nx; ++j )
    {
        for ( i = 0; i < nx; ++i)
        {
            node[0] = 3 * ( j + 1 ) * nx + (( j + 1 )<<1) + (( i + 1 )<<1);
            node[1] = node[0] - 1;
            node[2] = node[0] - 2;
            node[3] = ( 3 * ( j + 1 ) - 1 ) * nx + (( j + 1 )<<1) + ( i + 1 ) - 2;
            node[4] = ( 3 * ( j + 1 ) - 3 ) * nx + (( j + 1 )<<1) + (( i + 1 )<<1) - 4;
            node[5] = node[4] + 1;
            node[6] = node[4] + 2;
            node[7] = node[3] + 1;

            rho = 100.0 * r8_uniform_01 ( seed );

            #pragma omp parallel for num_threads(8)
            for ( krow = 0; krow < 8; ++krow )
            {
                ii = node[krow];
                #pragma omp parallel for num_threads(8)
                for ( kcol = 0; kcol < 8; ++kcol )
                {
                    jj = node[kcol];
                    a[ii-jj+ml+mu+jj*lda] = a[ii-jj+ml+mu+jj*lda]+ rho * em[krow+kcol*8];
                }
            }
        }
    }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wathen_ge ( void * data)
/******************************************************************************/
/*
  Purpose:
    WATHEN_GE returns the Wathen matrix as a general storage (GE) matrix.
  Discussion:
    The Wathen matrix is a finite element matrix which is sparse.
    The entries of the matrix depend in part on a physical quantity
    related to density.  That density is here assigned random values between
    0 and 100.
    The matrix order N is determined by the input quantities NX and NY,
    which would usually be the number of elements in the X and Y directions.
    The value of N is
      N = 3*NX*NY + 2*NX + 2*NY + 1,
    The matrix is the consistent mass matrix for a regular NX by NY grid
    of 8 node serendipity elements.
    The local element numbering is
      3--2--1
      |     |
      4     8
      |     |
      5--6--7
    Here is an illustration for NX = 3, NY = 2:
     23-24-25-26-27-28-29
      |     |     |     |
     19    20    21    22
      |     |     |     |
     12-13-14-15-16-17-18
      |     |     |     |
      8     9    10    11
      |     |     |     |
      1--2--3--4--5--6--7
    For this example, the total number of nodes is, as expected,
      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
    The matrix is symmetric positive definite for any positive values of the
    density RHO(X,Y).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2014
  Author:
    John Burkardt
  Reference:
    Nicholas Higham,
    Algorithm 694: A Collection of Test Matrices in MATLAB,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 3, September 1991, pages 289-305.
    Andrew Wathen,
    Realistic eigenvalue bounds for the Galerkin mass matrix,
    IMA Journal of Numerical Analysis,
    Volume 7, Number 4, October 1987, pages 449-457.
  Parameters:
    Input, int NX, NY, values which determine the size
    of the matrix.
    Input, int N, the number of rows and columns.
    Input/output, int *SEED, the random number seed.
    Output, double WATHEN_GE[N*N], the matrix.
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ nx = s_data->a0;
	const register dim_typ ny = s_data->a1;
	const register dim_typ n = s_data->a2;
	int * seed = s_data->a3;
	
    ityp *a;
    const ityp em[64] =
    {
        6.00, -6.00,  2.00, -8.00,  3.00, -8.00,  2.00, -6.00,
        -6.00, 32.00, -6.00, 20.00, -8.00, 16.00, -8.00, 20.00,
        2.00, -6.00,  6.00, -6.00,  2.00, -8.00,  3.00, -8.00,
        -8.00, 20.00, -6.00, 32.00, -6.00, 20.00, -8.00, 16.00,
        3.00, -8.00,  2.00, -6.00,  6.00, -6.00,  2.00, -8.00,
        -8.00, 16.00, -8.00, 20.00, -6.00, 32.00, -6.00, 20.00,
        2.00, -8.00,  3.00, -8.00,  2.00, -6.00,  6.00, -6.00,
        -6.00, 20.00, -8.00, 16.00, -8.00, 20.00, -6.00, 32.00
    };
    dim_typ i;
    dim_typ ii;
    dim_typ j;
    dim_typ jj;
    dim_typ kcol;
    dim_typ krow;
    int node[8];
    ityp rho;

    a = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            a[i+j*n] = 0.00;

    for ( j = 0; j < nx; ++j)
    {
        for ( i = 0; i < nx; ++i )
        {
            node[0] = 3 * ( j + 1 ) * nx + (( j + 1 )<<1) + (( i + 1 )<<1);
            node[1] = node[0] - 1;
            node[2] = node[0] - 2;
            node[3] = ( 3 * ( j + 1 ) - 1 ) * nx + (( j + 1 )<<1) + ( i + 1 ) - 2;
            node[4] = ( 3 * ( j + 1 ) - 3 ) * nx + (( j + 1 )<<1) + (( i + 1 )<<1) - 4;
            node[5] = node[4] + 1;
            node[6] = node[4] + 2;
            node[7] = node[3] + 1;

            rho = 100.00 * r8_uniform_01 ( seed );

            #pragma omp parallel for num_threads(8)
            for ( krow = 0; krow < 8; ++krow )
            {
                ii = node[krow];
                #pragma omp parallel for num_threads(8)
                for ( kcol = 0; kcol < 8; ++kcol)
                {
                    jj = node[kcol];
                    a[ii+jj*n] = a[ii+jj*n] + rho * em[krow+kcol*8];
                }
            }
        }
    }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wathen_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    WATHEN_ORDER returns the order of the WATHEN matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 June 2014
  Author:
    John Burkardt
  Reference:
    Nicholas Higham,
    Algorithm 694: A Collection of Test Matrices in MATLAB,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 3, September 1991, pages 289-305.
    Andrew Wathen,
    Realistic eigenvalue bounds for the Galerkin mass matrix,
    IMA Journal of Numerical Analysis,
    Volume 7, 1987, pages 449-457.
  Parameters:
    Input, int NX, NY, values which determine the size of A.
    Output, int WATHEN_ORDER, the order of the matrix,
    as determined by NX and NY.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ nx = a_data[0];
	const register dim_typ ny = a_data[1];
	
	result = 3 * nx * ny + (nx<<1) + (ny<<1) + 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wathen_st ( void * data)
/******************************************************************************/
/*
  Purpose:
    WATHEN_ST: Wathen matrix stored in sparse triplet (ST) format.
  Discussion:
    When dealing with sparse matrices in MATLAB, it can be much more efficient
    to work first with a triple of I, J, and X vectors, and only once
    they are complete, convert to MATLAB's sparse format.
    The Wathen matrix is a finite element matrix which is sparse.
    The entries of the matrix depend in part on a physical quantity
    related to density.  That density is here assigned random values between
    0 and 100.
    The matrix order N is determined by the input quantities NX and NY,
    which would usually be the number of elements in the X and Y directions.
    The value of N is
      N = 3*NX*NY + 2*NX + 2*NY + 1,
    The matrix is the consistent mass matrix for a regular NX by NY grid
    of 8 node serendipity elements.
    The local element numbering is
      3--2--1
      |     |
      4     8
      |     |
      5--6--7
    Here is an illustration for NX = 3, NY = 2:
     23-24-25-26-27-28-29
      |     |     |     |
     19    20    21    22
      |     |     |     |
     12-13-14-15-16-17-18
      |     |     |     |
      8     9    10    11
      |     |     |     |
      1--2--3--4--5--6--7
    For this example, the total number of nodes is, as expected,
      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
    The matrix is symmetric positive definite for any positive values of the
    density RHO(X,Y).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2014
  Author:
    John Burkardt.
  Reference:
    Nicholas Higham,
    Algorithm 694: A Collection of Test Matrices in MATLAB,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 3, September 1991, pages 289-305.
    Andrew Wathen,
    Realistic eigenvalue bounds for the Galerkin mass matrix,
    IMA Journal of Numerical Analysis,
    Volume 7, Number 4, October 1987, pages 449-457.
  Parameters:
    Input, int NX, NY, values which determine the size of
    the matrix.
    Input, int NZ_NUM, the number of values used to
    describe the matrix.
    Input/output, int *SEED, the random number seed.
    Output, int ROW[NZ_NUM], COL[NZ_NUM], the row and
    column indices of the nonzero entries.
    Output, double WATHEN_ST[NZ_NUM], the nonzero entries of the matrix.
*/
{
	const _3dtpi2pdt * const s_data = data;
	const register dim_typ nx = s_data->a0;
	const register dim_typ ny = s_data->a1;
	const register dim_typ nz_num = s_data->a2;
	int * seed = s_data->a3;
	dim_typ * row = s_data->a4;
	dim_typ * col = s_data->a5;
	
    ityp *a;
    const ityp em[64] =
    {
        6.00, -6.00,  2.00, -8.00,  3.00, -8.00,  2.00, -6.00,
        -6.00, 32.00, -6.00, 20.00, -8.00, 16.00, -8.00, 20.00,
        2.00, -6.00,  6.00, -6.00,  2.00, -8.00,  3.00, -8.00,
        -8.00, 20.00, -6.00, 32.00, -6.00, 20.00, -8.00, 16.00,
        3.00, -8.00,  2.00, -6.00,  6.00, -6.00,  2.00, -8.00,
        -8.00, 16.00, -8.00, 20.00, -6.00, 32.00, -6.00, 20.00,
        2.00, -8.00,  3.00, -8.00,  2.00, -6.00,  6.00, -6.00,
        -6.00, 20.00, -8.00, 16.00, -8.00, 20.00, -6.00, 32.00
    };
    dim_typ i;
    dim_typ j;
    dim_typ k;
    dim_typ kcol;
    dim_typ krow;
    int node[8];
    ityp rho;

    a = ( ityp * ) malloc ( nz_num * sizeof ( ityp ) );

    for ( k = 0; k < nz_num; ++k )
    {
        row[k] = col[k] = 0;
        a[k] = 0.00;
    }

    k = 0;

    for ( j = 0; j < nx; ++j )
    {
        for ( i = 0; i < nx; ++i )
        {
            node[0] = 3 * ( j + 1 ) * nx + (( j + 1 )<<1) + (( i + 1 )<<1);
            node[1] = node[0] - 1;
            node[2] = node[0] - 2;
            node[3] = ( 3 * ( j + 1 ) - 1 ) * nx + (( j + 1 )<<1) + ( i + 1 ) - 2;
            node[4] = ( 3 * ( j + 1 ) - 3 ) * nx + (( j + 1 )<<1) + (( i + 1 )<<1) - 4;
            node[5] = node[4] + 1;
            node[6] = node[4] + 2;
            node[7] = node[3] + 1;

            rho = 100.00 * r8_uniform_01 ( seed );

            #pragma omp parallel for num_threads(8)
            for ( krow = 0; krow < 8; ++krow )
                #pragma omp parallel for num_threads(8)
                for ( kcol = 0; kcol < 8; ++kcol )
                {
                    row[k] = node[krow];
                    col[k] = node[kcol];
                    a[k] = rho * em[krow+(kcol<<3)];
                    ++ k;
                }
        }
    }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wathen_st_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    WATHEN_ST_SIZE: Size of Wathen matrix stored in sparse triplet format.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 June 2014
  Author:
    John Burkardt.
  Reference:
    Nicholas Higham,
    Algorithm 694: A Collection of Test Matrices in MATLAB,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 3, September 1991, pages 289-305.
    Andrew Wathen,
    Realistic eigenvalue bounds for the Galerkin mass matrix,
    IMA Journal of Numerical Analysis,
    Volume 7, Number 4, October 1987, pages 449-457.
  Parameters:
    Input, integer NX, NY, values which determine the size of the matrix.
    Output, integer NZ_NUM, the number of items of data used to describe
    the matrix.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ nx = a_data[0];
	const register dim_typ ny = a_data[1];
	
	result = nx * (ny << 6);
    return &result;
}

#endif
