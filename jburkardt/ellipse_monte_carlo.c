#ifndef __DISABLEDEEP_ELLIPSEMONTECARLO

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ellipse_area1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ELLIPSE_AREA1 returns the area of an ellipse defined by a matrix.
  Discussion:
    The points X in the ellipse are described by a 2 by 2
    positive definite symmetric matrix A, and a "radius" R, such that
      X' * A * X <= R * R
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A[2*2], the matrix that describes
    the ellipse.  A must be symmetric and positive definite.
    Input, double R, the "radius" of the ellipse.
    Output, double ELLIPSE_AREA1, the area of the ellipse.
*/
{
	static ityp result = MAX_VAL;
	
	const itpit * const s_data = data;
	
	const register ityp r = s_data->a0;
	ityp * a = s_data->a1;
	
	result = r * r * M_PI / sqrt ( ( a[0+0*2] * a[1+1*2] - a[1+0*2] * a[0+1*2] ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ellipse_area2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    ELLIPSE_AREA2 returns the area of an ellipse defined by an equation.
  Discussion:
    The ellipse is described by the formula
      a x^2 + b xy + c y^2 = d
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, C, coefficients on the left hand side.
    Input, double D, the right hand side.
    Output, double ELLIPSE_AREA2, the area of the ellipse.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	ityp a = a_data[0];
	ityp b = a_data[1];
	ityp c = a_data[2];
	ityp d = a_data[3];
	
	result = d * M_PI / ( 4.00 * a * c - b * b );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ellipse_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    ELLIPSE_SAMPLE samples points in an ellipse.
  Discussion:
    The points X in the ellipsoid are described by a 2 by 2 positive
    definite symmetric matrix A, and a "radius" R, such that
      X' * A * X <= R * R
    The algorithm computes the Cholesky factorization of A:
      A = U' * U.
    A set of uniformly random points Y is generated, satisfying:
      Y' * Y <= R * R.
    The appropriate points in the ellipsoid are found by solving
      U * X = Y
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 May 2005
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
    Input, int N, the number of points.
    Input, double A[2*2], the matrix that describes the ellipse.
    Input, double R, the right hand side of the ellipse equation.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double ELLIPSE_SAMPLE[2*N], the points.
*/
{
	const dtpititpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp r = s_data->a2;
	int * seed = s_data->a3;
	
    dim_typ i;
    dim_typ info;
    dim_typ j;
    dim_typ k;
    static dim_typ m = 2;
    ityp *u = ( ityp * ) malloc ( m * m * sizeof ( dim_typ ) );
    ityp *x;
    /*
    Get the upper triangular Cholesky factor U of A.
    */

    for ( j = 0; j < m; ++j)
        for ( i = 0; i < m; ++i)
            u[i+j*m] = a[i+j*m];

    info = r8po_fa ( m, m, u);

    if ( info != 0 )
        return NULL;
    /*
    Get the points Y that satisfy Y' * Y <= R * R.
    */
    x = uniform_in_sphere01_map ( m, n, seed );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < m; ++i )
            x[i+j*m] *= r;
    /*
    Solve U * X = Y.
    */
    for ( j = 0; j < n; ++j)
        r8po_sl ( m, m, u, x+j*m );
    free ( u );
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8po_fa ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8PO_FA factors a real symmetric positive definite matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 May 2005
  Author:
    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
    C++ version by John Burkardt.
  Reference:
    Dongarra, Moler, Bunch and Stewart,
    LINPACK User's Guide,
    SIAM, (Society for Industrial and Applied Mathematics),
    3600 University City Science Center,
    Philadelphia, PA, 19104-2688.
    ISBN 0-89871-172-X
  Parameters:
    Input/output, double A[LDA*N].  On input, the symmetric matrix
    to be  factored.  Only the diagonal and upper triangle are used.
    On output, an upper triangular matrix R so that A = R'*R
    where R' is the transpose.  The strict lower triangle is unaltered.
    If INFO /= 0, the factorization is not complete.
    Input, int LDA, the leading dimension of the array A.
    Input, int N, the order of the matrix.
    Output, int R8PO_FA, error flag.
    0, for normal return.
    K, signals an error condition.  The leading minor of order K is not
    positive definite.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpit * const s_data = data;
	const register dim_typ lda = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ info;
    dim_typ j, k;
    ityp s;
    ityp t;

    for ( j = 1; j <= n; ++j )
    {
        s = 0.00;

        for ( k = 1; k <= j-1; ++k )
        {
            t = a[k-1+(j-1)*lda] - ddot ( k-1, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
            t /= a[k-1+(k-1)*lda];
            a[k-1+(j-1)*lda] = t;
            s += t * t;
        }

        s = a[j-1+(j-1)*lda] - s;

        if ( s <= 0.0 )
        {
            info = j;
            result = info;
            return &result;
        }

        a[j-1+(j-1)*lda] = sqrt ( s );
    }

    info = 0;
    result = info;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8po_sl ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8PO_SL solves a linear system factored by R8PO_FA.
  Discussion:
    A division by zero will occur if the input factor contains
    a zero on the diagonal.  Technically this indicates
    singularity but it is usually caused by improper subroutine
    arguments.  It will not occur if the subroutines are called
    correctly and INFO == 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 May 2005
  Author:
    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
    C++ version by John Burkardt.
  Reference:
    Dongarra, Moler, Bunch and Stewart,
    LINPACK User's Guide,
    SIAM, (Society for Industrial and Applied Mathematics),
    3600 University City Science Center,
    Philadelphia, PA, 19104-2688.
    ISBN 0-89871-172-X
  Parameters:
    Input, double A[LDA*N], the output from R8PO_FA.
    Input, int LDA, the leading dimension of the array A.
    Input, int N, the order of the matrix.
    Input/output, double B[N].  On input, the right hand side.
    On output, the solution.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ lda = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	
    dim_typ k;
    ityp t;
    /*
    Solve R' * Y = B.
    */
    for ( k = 1; k <= n; ++k )
    {
        t = ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
        b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
    }
    /*
    Solve R * X = Y.
    */
    for ( k = n; 1 <= k; --k)
    {
        b[k-1] /= a[k-1+(k-1)*lda];
        t = -b[k-1];
        daxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _uniform_in_sphere01_map ( void * data)
/******************************************************************************/
/*
  Purpose:
    UNIFORM_IN_SPHERE01_MAP maps uniform points into the unit sphere.
  Discussion:
    The sphere has center 0 and radius 1.
    We first generate a point ON the sphere, and then distribute it
    IN the sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2004
  Author:
    John Burkardt
  Reference:
    Russell Cheng,
    Random Variate Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998, pages 168.
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.
  Parameters:
    Input, int DIM_NUM, the dimension of the space.
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double X[DIM_NUM*N], the points.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
    ityp exponent;
    dim_typ i, j;
    ityp norm;
    ityp r;
    ityp *v = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );
    ityp *x = ( ityp * ) malloc ( dim_num * n * sizeof ( ityp ) );

    exponent = 1.00 / ( ityp ) ( dim_num );

    for ( j = 0; j < n; ++j )
    {
        /*
        Fill a vector with normally distributed values.
        */
        r8vec_normal_01 ( dim_num, seed, v );
        /*
        Compute the length of the vector.
        */
        norm = 0.00;
        for ( i = 0; i < dim_num; ++i)
            norm += pow ( v[i], 2 );
        norm = sqrt ( norm );
        /*
        Normalize the vector.
        */
        for ( i = 0; i < dim_num; ++i )
            v[i] /= norm;
        /*
        Now compute a value to map the point ON the sphere INTO the sphere.
        */
        r = r8_uniform_01 ( seed );
        r = pow ( r, exponent );

        for ( i = 0; i < dim_num; ++i )
            x[i+j*dim_num] = r * v[i];
    }

    free ( v );
    return x;
}

#endif
