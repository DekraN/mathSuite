#ifndef __DISABLEDEEP_JACOBIEIGENVALUE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _jacobi_eigenvalue ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
  Discussion:
    This function computes the eigenvalues and eigenvectors of a
    real symmetric matrix, using Rutishauser's modfications of the classical
    Jacobi rotation method with threshold pivoting.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 September 2013
  Author:
    C version by John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the matrix, which must be square, real,
    and symmetric.
    Input, int IT_MAX, the maximum number of iterations.
    Output, double V[N*N], the matrix of eigenvectors.
    Output, double D[N], the eigenvalues, in descending order.
    Output, int *IT_NUM, the total number of iterations.
    Output, int *ROT_NUM, the total number of rotations.
*/
{
	const dtpiti2pit2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	int it_max = s_data->a2;
	ityp * v = s_data->a3;
	ityp * d = s_data->a4;
	int * it_num = s_data->a5;
	int * rot_num = s_data->a6;
	
    ityp *bw;
    ityp c;
    ityp g;
    ityp gapq;
    ityp h;
    dim_typ i;
    dim_typ j;
    dim_typ k;
    dim_typ l;
    dim_typ m;
    dim_typ p;
    dim_typ q;
    ityp s;
    ityp t;
    ityp tau;
    ityp term;
    ityp termp;
    ityp termq;
    ityp theta;
    ityp thresh;
    ityp w;
    ityp *zw;

    r8mat_identity ( n, v );
    r8mat_diag_get_vector ( n, a, d );

    bw = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    zw = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        bw[i] = d[i];
        zw[i] = 0.0;
    }
    *it_num = *rot_num = 0;

    while ( *it_num < it_max )
    {
        ++ *it_num;
        /*
        The convergence threshold is based on the size of the elements in
        the strict upper triangle of the matrix.
        */
        thresh = 0.00;
        for ( j = 0; j < n; ++j )
            for ( i = 0; i < j; ++i )
                thresh += a[i+j*n] * a[i+j*n];

        thresh = sqrt ( thresh ) / ( ityp ) ( n<<2 );

        if ( thresh == 0.00 )
            break;

        for ( p = 0; p < n; ++p )
        {
            for ( q = p + 1; q < n; ++q )
            {
                gapq = 10.00 * fabs ( a[p+q*n] );
                termp = gapq + fabs ( d[p] );
                termq = gapq + fabs ( d[q] );
                /*
                Annihilate tiny offdiagonal elements.
                */
                if ( 4 < *it_num &&termp == fabs ( d[p] ) && termq == fabs ( d[q] ) )
                    a[p+q*n] = 0.00;
                /*
                Otherwise, apply a rotation.
                */
                else if ( thresh <= fabs ( a[p+q*n] ) )
                {
                    h = d[q] - d[p];
                    term = fabs ( h ) + gapq;

                    if ( term == fabs ( h ) )
                        t = a[p+q*n] / h;
                    else
                    {
                        theta = 0.50 * h / a[p+q*n];
                        t = 1.00 / ( fabs ( theta ) + sqrt ( 1.00 + theta * theta ) );
                        if ( theta < 0.00 )
                            t *= -1;
                    }
                    c = 1.00 / sqrt ( 1.00 + t * t );
                    s = t * c;
                    tau = s / ( 1.00 + c );
                    h = t * a[p+q*n];
                    /*
                    Accumulate corrections to diagonal elements.
                    */
                    zw[p] = zw[p] - h;
                    zw[q] = zw[q] + h;
                    d[p] = d[p] - h;
                    d[q] = d[q] + h;

                    a[p+q*n] = 0.0;
                    /*
                    Rotate, using information from the upper triangle of A only.
                    */
                    for ( j = 0; j < p; ++j )
                    {
                        g = a[j+p*n];
                        h = a[j+q*n];
                        a[j+p*n] = g - s * ( h + g * tau );
                        a[j+q*n] = h + s * ( g - h * tau );
                    }

                    for ( j = p + 1; j < q; ++j )
                    {
                        g = a[p+j*n];
                        h = a[j+q*n];
                        a[p+j*n] = g - s * ( h + g * tau );
                        a[j+q*n] = h + s * ( g - h * tau );
                    }

                    for ( j = q + 1; j < n; ++j )
                    {
                        g = a[p+j*n];
                        h = a[q+j*n];
                        a[p+j*n] = g - s * ( h + g * tau );
                        a[q+j*n] = h + s * ( g - h * tau );
                    }
                    /*
                    Accumulate information in the eigenvector matrix.
                    */
                    for ( j = 0; j < n; ++j)
                    {
                        g = v[j+p*n];
                        h = v[j+q*n];
                        v[j+p*n] = g - s * ( h + g * tau );
                        v[j+q*n] = h + s * ( g - h * tau );
                    }
                    ++ *rot_num;
                }
            }
        }

        for ( i = 0; i < n; ++i )
        {
            bw[i] = bw[i] + zw[i];
            d[i] = bw[i];
            zw[i] = 0.00;
        }
    }
    /*
    Restore upper triangle of input matrix.
    */
    for ( j = 0; j < n; ++j)
        for ( i = 0; i < j; ++i)
            a[i+j*n] = a[j+i*n];
    /*
    Ascending sort the eigenvalues and eigenvectors.
    */
    for ( k = 0; k < n - 1; ++k )
    {
        m = k;
        for ( l = k + 1; l < n; ++l )
            if ( d[l] < d[m] )
                m = l;

        if ( m != k )
        {
            t    = d[m];
            d[m] = d[k];
            d[k] = t;
            for ( i = 0; i < n; ++i )
            {
                w        = v[i+m*n];
                v[i+m*n] = v[i+k*n];
                v[i+k*n] = w;
            }
        }
    }

    free ( bw );
    free ( zw );

    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_diag_get_vector ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of rows and columns of the matrix.
    Input, double A[N*N], the N by N matrix.
    Output, double V[N], the diagonal entries
    of the matrix.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * v = s_data->a2;
	
    for (dim_typ i = 0; i < n; ++i )
        v[i] = a[i+i*n];
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _r8mat_identity ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_IDENTITY sets an R8MAT to the identity matrix.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of A.
    Output, double A[N*N], the N by N identity matrix.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i, j, k = 0;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
        {
            a[k] = i == j;
            ++ k;
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_is_eigen_right ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_IS_EIGEN_RIGHT determines the error in a (right) eigensystem.
  Discussion:
    An R8MAT is a matrix of doubles.
    This routine computes the Frobenius norm of
      A * X - X * LAMBDA
    where
      A is an N by N matrix,
      X is an N by K matrix (each of K columns is an eigenvector)
      LAMBDA is a K by K diagonal matrix of eigenvalues.
    This routine assumes that A, X and LAMBDA are all real.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, int K, the number of eigenvectors.
    K is usually 1 or N.
    Input, double A[N*N], the matrix.
    Input, double X[N*K], the K eigenvectors.
    Input, double LAMBDA[K], the K eigenvalues.
    Output, double R8MAT_IS_EIGEN_RIGHT, the Frobenius norm
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
	
    ityp error_frobenius;
    dim_typ i, j, l;
    ityp * c = ( ityp * ) malloc ( n * k * sizeof ( ityp ) );

    for ( j = 0; j < k; ++j )
        for ( i = 0; i < n; ++i )
        {
            c[i+j*n] = 0.0;
            for ( l = 0; l < n; ++l )
                c[i+j*n] += a[i+l*n] * x[l+j*n];
        }

        for ( j = 0; j < k; ++j )
            for ( i = 0; i < n; ++i )
                c[i+j*n] -= lambda[j] * x[i+j*n];

    error_frobenius = r8mat_norm_fro ( n, k, c );
    free ( c );
    
    result = error_frobenius;
    return &result;
}

#endif

