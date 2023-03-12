#ifndef __DISABLEDEEP_DIFFER

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _differ_backward ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFER_BACKWARD computes backward difference coefficients.
  Discussion:
    We determine coefficients C to approximate the derivative at X0
    of order O and precision P, using equally spaced backward
    differences, so that
      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x-ih) + O(h^(p))
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, double H, the spacing.  0 < H.
    Input, int O, the order of the derivative to be
    approximated.  1 <= O.(const register ityp h, const register dim_typ o, const register dim_typ p, ityp c[static o+p], ityp x[static o+p] )
    Input, int P, the order of the error, as a power of H.
    Output, double C[O+P], the coefficients.
    Output, double X[O+P], the evaluation points.
*/
{
	const it2dt2pit * const s_data = data;
	const register ityp h = s_data->a0;
	const register dim_typ o = s_data->a1;
	const register dim_typ p = s_data->a2;
	ityp * c = s_data->a3;
	ityp * x = s_data->a4;
	
    ityp *b;
    dim_typ i;
    dim_typ info;
    dim_typ job;
    dim_typ n = o+p;
    ityp t;

    for ( i = 0; i < n; ++i )
        x[i] = ( ityp ) ( i + 1 - n ) * h;

	b = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        b[i] = 0.00;

    b[o] = 1.00;

    job = 0;
    r8vm_sl ( n, x, b, job, c, &info );

    if(info)
        return NULL;

    t = r8_factorial ( o );
    for ( i = 0; i < n; ++i )
    c[i] *= t;

    free ( b );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _differ_central ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFER_CENTRAL computes central difference coefficients.
  Discussion:
    We determine coefficients C to approximate the derivative at X0
    of order O and precision P, using equally spaced central
    differences, so that
      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x+(2*i-o-p+1)*h/2)
        + O(h^(p))
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, double H, the spacing.  0 < H.
    Input, int O, the order of the derivative to
    be approximated.  1 <= O.
    Input, int P, the order of the error, as a power of H.
    Output, double C[O+P], the coefficients.
    Output, double X[O+P], the evaluation points.
*/
{
	const it2dt2pit * const s_data = data;
	const register ityp h = s_data->a0;
	const register dim_typ o = s_data->a1;
	const register dim_typ p = s_data->a2;
	ityp * c = s_data->a3;
	ityp * x = s_data->a4;
	
    ityp *b;
    dim_typ i;
    dim_typ info;
    dim_typ job;
    dim_typ n = o+p;
    ityp t;

    for ( i = 0; i < n; ++i )
        x[i] = ( ityp ) ( - n + 1 + (i<<1) ) * h / 2.00;

    b = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        b[i] = 0.00;
    b[o] = 1.00;

    job = 0;
    r8vm_sl ( n, x, b, job, c, &info );

    if(info)
        return NULL;

    t = r8_factorial ( o );
    for ( i = 0; i < n; ++i )
        c[i] *= t;

    free ( b );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _differ_forward ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFER_FORWARD computes forward difference coefficients.
  Discussion:
    We determine coefficients C to approximate the derivative at X0
    of order O and precision P, using equally spaced forward
    differences, so that
      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x+ih) + O(h^(p))
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, real H, the spacing.  0 < H.
    Input, integer O, the order of the derivative to be approximated.
    1 <= O.
    Input, integer P, the order of the error, as a power of H.
    Output, real C[O+P], the coefficients.
    Output, real X[O+P], the evaluation points.
*/
{
	const it2dt2pit * const s_data = data;
	const register ityp h = s_data->a0;
	const register dim_typ o = s_data->a1;
	const register dim_typ p = s_data->a2;
	ityp * c = s_data->a3;
	ityp * x = s_data->a4;
	
    ityp *b; 
    dim_typ i;
    dim_typ info;
    dim_typ job;
    dim_typ n = o+p;
    ityp t;

    for ( i = 0; i < n; ++i )
        x[i] = ( ityp ) ( i ) * h;
        
	b = ( ityp * ) malloc ( n * sizeof ( ityp ) );
	
    for ( i = 0; i < n; ++i )
        b[i] = 0.00;
    b[o] = 1.00;

    job = 0;
    r8vm_sl ( n, x, b, job, c, &info );

    if(info)
        return NULL;

    t = r8_factorial ( o );
    for ( i = 0; i < n; ++i )
        c[i] *= t;

    free ( b );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _differ_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFER_INVERSE returns the inverse of the DIFFER matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, double STENCIL[N], the values that define A.
    Output, double DIFFER_INVERSE[N*N], the matrix.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * stencil = s_data->a1;
	
    ityp *a = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
    dim_typ i;
    dim_typ indx;
    dim_typ j;
    dim_typ k;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i)
            a[i+j*n] = j == 0;

    for ( i = 0; i < n; ++i)
    {
	    indx = 0;
	
	    for ( k = 0; k < n; ++k )
	        if ( k != i )
	        {
	            for ( j = indx + 1; 0 <= j; --j)
	            {
	                a[i+j*n] = - stencil[k] * a[i+j*n] / ( stencil[i] - stencil[k] );
	                if ( 0 < j )
	                a[i+j*n] += a[i+(j-1)*n] / ( stencil[i] - stencil[k] );
	            }
	            ++ indx;
	        }
	
	    for ( j = 0; j < n; ++j )
	        for ( i = 0; i < n; ++i )
	            a[i+j*n] /= stencil[i];
	}

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _differ_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFER_MATRIX computes the stencil matrix from the stencil vector.
  Discussion:
    If N = 4, and STENCIL = ( -3, -2, -1, 1 ), then A will be
    -3  -2  -1  1
     9   4   1  1
   -27  -8  -1  1
    81  16   1  1
    This matrix is a generalized form of a Vandermonde matrix A2:
     1   1   1  1
    -3  -2  -1  1
     9   4   1  1
   -27  -8  -1  1
    and if A * x = b, the A2 * x2 = b, where x2(i) = x(i) * stencil(i)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of stencil points.
    Input, double STENCIL[N], the stencil vector.
    The entries in this vector must be distinct.
    No entry of STENCIL may be 0.
    Output, double DIFFER_MATRIX[N*N], the stencil matrix.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * stencil = s_data->a1;
	
    ityp *a = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
    dim_typ i, j;

    for ( j = 0; j < n; ++j )
    {
        a[0+j*n] = stencil[j];
        for ( i = 1; i < n; ++i )
            a[i+j*n] *= stencil[j];
    }
    return a;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _differ_solve ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFER_SOLVE solves for finite difference coefficients.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of stencil points.
    Input, double STENCIL[N], the stencil vector.
    The entries in this vector must be distinct.
    No entry of STENCIL may be 0.
    Input, int ORDER, the order of the derivative to
    be approximated.  1 <= ORDER <= N.
    Output, double DIFFER_SOLVE[N], the coefficients to be used
    to multiply U(STENCIL(I))-U(0), so that the sum forms an
    approximation to the derivative of order ORDER, with error
    of order H^(N+1-ORDER).
*/
{
	const _2dtpit * const s_data = data; 
	
	const register dim_typ n = s_data->a0;
	const register dim_typ order = s_data->a1;
	ityp * stencil = s_data->a2;
	
    ityp *a = differ_matrix(n,stencil);
    ityp *b = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    ityp *c;
    dim_typ i;

    for ( i = 0; i < n; ++i )
        b[i] = 0.00;
    b[order-1] = 1.00;
    /*
    Solve A * C = B.
    */
    c = r8mat_fs_new ( n, a, b );

    free ( a );
    free ( b );
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _differ_stencil (void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFER_STENCIL computes finite difference coefficients.
  Discussion:
    We determine coefficients C to approximate the derivative at X0
    of order O and precision P, using finite differences, so that
      d^o f(x)/dx^o (x0) = sum ( 0 <= i <= o+p-1 ) c(i) f(x(i))
        + O(h^(p))
    where H is the maximum spacing between X0 and any X(I).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, double X0, the point where the derivative is to
    be approximated.
    Input, int O, the order of the derivative to be
    approximated.  1 <= O.
    Input, int P, the order of the error, as a power of H.
    Input, double X[O+P], the evaluation points.
    Output, double C[O+P], the coefficients.
*/
{
	const it2dt2pit * const s_data = data;
	const register ityp x0 = s_data->a0;
	const register dim_typ o = s_data->a1;
	const register dim_typ p = s_data->a2;
	ityp * x = s_data->a3;
	ityp * c = s_data->a4;
	
    ityp *b;
    ityp *dx;
    dim_typ i;
    dim_typ info;
    dim_typ job;
    dim_typ n = o+p;
    ityp t;

	dx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
	
    for ( i = 0; i < n; ++i )
        dx[i] = x[i] - x0;
        
    b = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        b[i] = 0.00;
    b[o] = 1.00;

    job = 0;
    r8vm_sl ( n, dx, b, job, c, &info );

    if(info)
        return NULL;

    t = r8_factorial ( o );
    for ( i = 0; i < n; ++i )
        c[i] *= t;

    free ( b );
    free ( dx );
    return NULL;
}

// LIGHT FACTORIAL VERSION
// simple version, with memoizers systems
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_factorial ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_FACTORIAL computes the factorial of N.
  Discussion:
    factorial ( N ) = product ( 1 <= I <= N ) I
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 June 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the argument of the factorial function.
    If N is less than 1, the function value is returned as 1.
    Output, double R8_FACTORIAL, the factorial of N.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
	ityp value = 1.00;

	#pragma omp parallel for
	for(dim_typ i=1; i<=n; ++i)
		value *= i;

	result = value;
	return &result;
}

#define R8MAXFSNEW_INVALIDRETURNVALUE NULL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_fs_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_FS_NEW factors and solves a system with one right hand side.
  Discussion:
    This routine uses partial pivoting, but no pivot vector is required.
    This routine differs from R8MAT_FSS_NEW in two ways:
    * only one right hand side is allowed;
    * the input matrix A is not modified.
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 January 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be positive.
    Input, double A[N*N], the coefficient matrix of the linear system.
    Input, double B[N], on input, the right hand side of the
    linear system.
    Output, double R8MAT_FS[N], the solution of the linear system.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    ityp *a2 = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
    dim_typ i;
    dim_typ ipiv;
    dim_typ j;
    dim_typ jcol;
    ityp piv;
    ityp t;
    ityp *x = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < n; ++i )
            a2[i+j*n] = a[i+j*n];

    for ( i = 0; i < n; ++i)
        x[i] = b[i];

    for ( jcol = 1; jcol <= n; ++jcol )
    {
    /*
    Find the maximum element in column I.
    */
        piv = abs ( a2[jcol-1+(jcol-1)*n] );
        ipiv = jcol;
        for ( i = jcol+1; i <= n; ++i)
        {
            if ( piv < abs ( a2[i-1+(jcol-1)*n] ) )
            {
                piv = abs ( a2[i-1+(jcol-1)*n] );
                ipiv = i;
            }
        }

        if ( piv == 0.00 )
            return R8MAXFSNEW_INVALIDRETURNVALUE;
        /*
        Switch rows JCOL and IPIV, and X.
        */
        if ( jcol != ipiv )
        {
            for ( j = 1; j <= n; ++j)
            {
                t                  = a2[jcol-1+(j-1)*n];
                a2[jcol-1+(j-1)*n] = a2[ipiv-1+(j-1)*n];
                a2[ipiv-1+(j-1)*n] = t;
            }
            t         = x[jcol-1];
            x[jcol-1] = x[ipiv-1];
            x[ipiv-1] = t;
        }
        /*
        Scale the pivot row.
        */
        t = a2[jcol-1+(jcol-1)*n];
        a2[jcol-1+(jcol-1)*n] = 1.0;
        for ( j = jcol+1; j <= n; ++j )
            a2[jcol-1+(j-1)*n] = a2[jcol-1+(j-1)*n] / t;
        x[jcol-1] /= t;
        /*
        Use the pivot row to eliminate lower entries in that column.
        */
        for ( i = jcol+1; i <= n; ++i )
        {
            if ( a2[i-1+(jcol-1)*n] != 0.00 )
            {
                t = - a2[i-1+(jcol-1)*n];
                a2[i-1+(jcol-1)*n] = 0.00;
                for ( j = jcol+1; j <= n; ++j)
                    a2[i-1+(j-1)*n] += t * a2[jcol-1+(j-1)*n];
                x[i-1] += t * x[jcol-1];
            }
        }
    }
    /*
    Back solve.
    */
    for ( jcol = n; 2 <= jcol; --jcol )
        for ( i = 1; i < jcol; i++ )
            x[i-1] -= a2[i-1+(jcol-1)*n] * x[jcol-1];

    free ( a2 );
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  inline void   * _r8mat_sub_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_SUB_NEW computes C = A - B.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
    For this routine, the result is returned as the function value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the order of the matrices.
    Input, double A[M*N], double B[M*N], the matrices.
    Output, double R8MAT_SUB_NEW[M*N], the value of A-B.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	
    ityp *c = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );
    dim_typ i, j;
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i)
            c[i+j*m] -= b[i+j*m];
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vm_sl ( void * data) 
/******************************************************************************/
/*
  Purpose:
    R8VM_SL solves a R8VM system.
  Discussion
    The R8VM storage format is used for an M by N Vandermonde matrix.
    An M by N Vandermonde matrix is defined by the values in its second
    row, which will be written here as X(1:N).  The matrix has a first
    row of 1's, a second row equal to X(1:N), a third row whose entries
    are the squares of the X values, up to the M-th row whose entries
    are the (M-1)th powers of the X values.  The matrix can be stored
    compactly by listing just the values X(1:N).
    Vandermonde systems are very close to singularity.  The singularity
    gets worse as N increases, and as any pair of values defining
    the matrix get close.  Even a system as small as N = 10 will
    involve the 9th power of the defining values.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 January 2013
  Author:
    Original FORTRAN77 version by Golub, VanLoan.
    C version by John Burkardt.
  Reference:
    Gene Golub, Charles Van Loan,
    Matrix Computations,
    Third Edition,
    Johns Hopkins, 1996.
  Parameters:
    Input, int N, the number of rows and columns of the matrix.
    Input, double A[N], the R8VM matrix.
    Input, double B[N], the right hand side.
    Input, int JOB, specifies the system to solve.
    0, solve A * x = b.
    nonzero, solve A' * x = b.
    Output, double X[N], the solution of the linear system.
    Output, int *INFO.
    0, no error.
    nonzero, at least two of the values in A are equal.
*/
{
	const dt2pitdtpitps * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	const register dim_typ job = s_data->a3;
	ityp * x = s_data->a4;
	short * info = s_data->a5;
	
    dim_typ i, j;
    /*
    Check for explicit singularity.
    */
    *info = 0;

    for ( j = 0; j < n; ++j )
    {
        for ( i = j+1; i < n; ++i )
        {
            if ( a[i] == a[j] )
            {
                *info = 1;
                return NULL;
            }
        }
    }

    for ( i = 0; i < n; ++i )
        x[i] = b[i];

    if ( !job )
    {
        for ( j = 1; j <= n-1; ++j )
            for ( i = n; j+1 <= i; --i )
                x[i-1] -= a[j-1] * x[i-2];

        for ( j = n-1; 1 <= j; --j )
        {
            for ( i = j+1; i <= n; ++i )
                x[i-1] /= ( a[i-1] - a[i-j-1] );

            for ( i = j; i <= n-1; ++i )
                x[i-1] -= x[i];
        }
    }
    else
    {
        for ( j = 1; j <= n-1; ++j )
            for ( i = n; j+1 <= i; --i )
                x[i-1] = ( x[i-1] - x[i-2] ) / ( a[i-1] - a[i-j-1] );

        for ( j = n-1; 1 <= j; --i )
            for ( i = j; i <= n-1; ++i )
                x[i-1] -= x[i] * a[j-1];
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vm_sl_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VM_SL_NEW solves a R8VM system.
  Discussion:
    The R8VM storage format is used for an M by N Vandermonde matrix.
    An M by N Vandermonde matrix is defined by the values in its second
    row, which will be written here as X(1:N).  The matrix has a first
    row of 1's, a second row equal to X(1:N), a third row whose entries
    are the squares of the X values, up to the M-th row whose entries
    are the (M-1)th powers of the X values.  The matrix can be stored
    compactly by listing just the values X(1:N).
    Vandermonde systems are very close to singularity.  The singularity
    gets worse as N increases, and as any pair of values defining
    the matrix get close.  Even a system as small as N = 10 will
    involve the 9th power of the defining values.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 January 2013
  Author:
    Original FORTRAN77 version by Golub, VanLoan.
    C version by John Burkardt.
  Reference:
    Gene Golub, Charles Van Loan,
    Matrix Computations,
    Third Edition,
    Johns Hopkins, 1996.
  Parameters:
    Input, int N, the number of rows and columns of the matrix.
    Input, double A[N], the R8VM matrix.
    Input, double B[N], the right hand side.
    Input, int JOB, specifies the system to solve.
    0, solve A * x = b.
    nonzero, solve A' * x = b.
    Output, int *INFO.
    0, no error.
    nonzero, at least two of the values in A are equal.
    Output, double R8VM_SL[N], the solution of the linear system.
*/
{
	const dt2pitdtps * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	const register dim_typ job = s_data->a3;
	short * info = s_data->a4;
	
    dim_typ i, j;
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    /*
    Check for explicit singularity.
    */
    *info = 0;

    for ( j = 0; j < n; ++j )
    {
        for ( i = j+1; i < n; ++i )
        {
            if ( a[i] == a[j] )
            {
                *info = 1;
                return NULL;
            }
        }
    }

    for ( i = 0; i < n; ++i )
        x[i] = b[i];

    if(!job)
    {
        for ( j = 1; j <= n-1; ++j )
            for ( i = n; j+1 <= i; --i )
                x[i-1] -= a[j-1] * x[i-2];

        for ( j = n-1; 1 <= j; --j )
        {
            for ( i = j+1; i <= n; ++i )
                x[i-1] /= ( a[i-1] - a[i-j-1] );

            for ( i = j; i <= n-1; ++i )
                x[i-1] -= x[i];
        }
    }
    else
    {
        for ( j = 1; j <= n-1; ++j )
            for ( i = n; j+1 <= i; --i )
                x[i-1] = ( x[i-1] - x[i-2] ) / ( a[i-1] - a[i-j-1] );

	    for ( j = n-1; 1 <= j; --j)
	        for ( i = j; i <= n-1; ++i )
	            x[i-1] -=  x[i] * a[j-1];

    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_mm_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_MM_NEW multiplies two matrices.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values,  stored as a vector
    in column-major order.
    For this routine, the result is returned as the function value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 April 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, N3, the order of the matrices.
    Input, ityp A[N1*N2], ityp B[N2*N3], the matrices to multiply.
    Output, ityp r8MAT_MM[N1*N3], the product matrix C = A * B.
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

    for ( i = 0; i < n1; ++i)
        for ( j = 0; j < n3; ++j )
        {
            c[i+j*n1] = 0.00;
            for ( k = 0; k < n2; ++k )
                c[i+j*n1] += a[i+k*n1] * b[k+j*n2];
        }

    return c;
}

#endif
