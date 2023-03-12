#ifndef __DISABLEDEEP_ELLIPSOIDMONTECARLO

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ellipsoid_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    ELLIPSOID_SAMPLE samples points uniformly from an ellipsoid.
  Discussion:
    The points X in the ellipsoid are described by a M by M
    positive definite symmetric matrix A, a "center" V, and
    a "radius" R, such that
   (X-V)' * A * (X-V) <= R * R
    The algorithm computes the Cholesky factorization of A:
      A = U' * U.
    A set of uniformly random points Y is generated, satisfying:
      Y' * Y <= R * R.
    The appropriate points in the ellipsoid are found by solving
      U * X = Y
      X = X + V
    Thanks to Dr Karl-Heinz Keil for pointing out that the original
    coding was actually correct only if A was replaced by its inverse.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 August 2014
  Author:
    John Burkardt
  Reference:
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.
  Parameters:
    Input, int M, the dimension of the space.
    Input, int N, the number of points.
    Input, double A[M*M], the matrix that describes
    the ellipsoid.
    Input, double V[M], the "center" of the ellipsoid.
    Input, double R, the "radius" of the ellipsoid.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double ELLIPSE_SAMPLE[M*N], the points.
*/
{
	const _2dt2pititpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * v = s_data->a3;
	const register ityp r = s_data->a4;
	int * seed = s_data->a5;
	
    dim_typ i, j;
    ityp *t;
    ityp *u = r8po_factor ( m, a );
    ityp *x;
    /*
    Get the Cholesky factor U.
    */

    if ( !u )
        return NULL;
    /*
    Get the points Y that satisfy Y' * Y <= 1.
    */
    x = uniform_in_sphere01_map ( m, n, seed );
    /*
    Get the points Y that satisfy Y' * Y <= R * R.
    */
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i)
            x[i+j*m] *= r;
    /*
    Solve U * X = Y.
    */
    for ( j = 0; j < n; ++j)
    {
        t = r8po_sl2 ( m, u, x + j * m );
        for ( i = 0; i < m; ++i )
            x[i+j*m] = t[i];
        free ( t );
    }
    /*
    X = X + V.
    */
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            x[i+j*m] += v[i];

    free ( u );
    return x;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ellipsoid_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    ELLIPSOID_VOLUME returns the volume of an ellipsoid.
  Discussion:
    The points X in the ellipsoid are described by an M by M
    positive definite symmetric matrix A, an M-dimensional point V,
    and a "radius" R, such that
   (X-V)' * A * (X-V) <= R * R
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, double A[M*M], the matrix that describes
    the ellipsoid.  A must be symmetric and positive definite.
    Input, double V[M], the "center" of the ellipse.
    The value of V is not actually needed by this function.
    Input, double R, the "radius" of the ellipse.
    Output, double ELLIPSOID_VOLUME, the volume of the ellipsoid.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit2pit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register ityp r = s_data->a1;
	ityp * a = s_data->a2;
	ityp * v = s_data->a3;
	
	
    ityp sqrt_det;
    ityp *u = r8po_factor ( m, a );
    sqrt_det = 1.00;
    for (dim_typ i = 0; i < m; ++i)
        sqrt_det *= u[i+i*m];
    free ( u );
    
    result = pow ( r, m ) * hypersphere_unit_volume ( m ) / sqrt_det;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hypersphere_unit_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_UNIT_VOLUME: volume of a unit hypersphere in M dimensions.
  Discussion:
    The unit sphere in M dimensions satisfies the equation:
      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
     M  Volume
     1    2
     2    1        * M_PI
     3 ( 4 /   3) * M_PI
     4 ( 1 /   2) * M_PI^2
     5 ( 8 /  15) * M_PI^2
     6 ( 1 /   6) * M_PI^3
     7 (16 / 105) * M_PI^3
     8 ( 1 /  24) * M_PI^4
     9 (32 / 945) * M_PI^4
    10 ( 1 / 120) * M_PI^5
    For the unit sphere, Volume(M) = 2 * M_PI * Volume(M-2)/ M
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of the space.
    Output, double HYPERSPHERE_UNIT_VOLUME, the volume of the sphere.
*/
{
	static ityp result = MAX_VAL;
	
 	const register dim_typ m = *(dim_typ *) data;
 
    dim_typ i, m2;
    ityp volume;
    if ( !(m%2))
    {
        m2 = m / 2;
        volume = 1.00;
        for ( i = 1; i <= m2; ++i )
            volume *= M_PI / ( ( ityp ) i );
    }
    else
    {
        m2 = ( m - 1 ) / 2;
        volume = pow ( M_PI, m2 ) * pow ( 2.00, m );
        for ( i = m2 + 1; i <= (m2<<1) + 1; ++i )
            volume /= ( ( ityp ) i );
    }
    
    result = volume;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8po_factor ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8PO_FA factors a R8PO matrix.
  Discussion:
    The R8PO storage format is appropriate for a symmetric positive definite
    matrix and its inverse. (The Cholesky factor of a R8PO matrix is an
    upper triangular matrix, so it will be in R8GE storage format.)
    Only the diagonal and upper triangle of the square array are used.
    This same storage format is used when the matrix is factored by
    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
    is set to zero.
    The positive definite symmetric matrix A has a Cholesky factorization
    of the form:
      A = R' * R
    where R is an upper triangular matrix with positive elements on
    its diagonal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 February 2013
  Author:
    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
    C version by John Burkardt.
  Reference:
    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the matrix in R8PO storage.
    Output, double R8PO_FA[N*N], the Cholesky factor in SGE
    storage, or NULL if there was an error.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp *b = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
    dim_typ i, j, k;
    ityp s;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            b[i+j*n] = a[i+j*n];

    for ( j = 0; j < n; j++ )
    {
        for ( k = 0; k <= j-1; k++ )
        {
                for ( i = 0; i <= k-1; i++ )
                    b[k+j*n] -= b[i+k*n] * b[i+j*n];
                b[k+j*n] /= b[k+k*n];
        }

        s = b[j+j*n];
        for ( i = 0; i <= j-1; ++i)
            s -= b[i+j*n] * b[i+j*n];

        if ( s <= 0.00)
        {
            free ( b );
            return NULL;
        }

        b[j+j*n] = sqrt ( s );
    }
    /*
    Since the Cholesky factor is in R8GE format, zero out the lower triangle.
    */
    for ( i = 0; i < n; ++i )
        for ( j = 0; j < i; ++j )
            b[i+j*n] = 0.00;

    return b;
}

#endif
