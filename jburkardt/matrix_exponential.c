#ifndef __DISABLEDEEP_MATRIXEXPONENTIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_norm_li ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_NORM_LI returns the matrix L-oo norm of an R8MAT.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
    The matrix L-oo norm is defined as:
      R8MAT_NORM_LI =  MAX ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
    The matrix L-oo norm is derived from the vector L-oo norm,
    and satisifies:
      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 December 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, double A[M*N], the matrix whose L-oo
    norm is desired.
    Output, double R8MAT_NORM_LI, the L-oo norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data; 
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
	dim_typ i, j;
	ityp row_sum;
	ityp value = 0.00;

	
	for ( i = 0; i < m; ++i )
	{
		row_sum = 0.00;
		for ( j = 0; j < n; ++j )
			row_sum += fabs ( a[i+j*m] );
		value = MAX ( value, row_sum );
	}
	
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_fss_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_FSS_NEW factors and solves a system with multiple right hand sides.
  Discussion:
    This routine uses partial pivoting, but no pivot vector is required.
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    28 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be positive.
    Input/output, double A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.
    Input, int NB, the number of right hand sides.
    Input, double B[N*NB], the right hand sides of the linear systems.
    Output, double R8MAT_FSS_NEW[N*NB], the solutions of the linear systems.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ nb = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	
	dim_typ i;
	int ipiv;
	dim_typ j;
	dim_typ jcol;
	ityp piv;
	ityp t;
	ityp *x = ( ityp * ) malloc ( n * nb * sizeof ( ityp ) );
	
	for ( j = 0; j < nb; ++j )
		for ( i = 0; i < n; ++i )
			x[i+j*n] = b[i+j*n];
	for ( jcol = 1; jcol <= n; ++jcol )
	{
		/*
		Find the maximum element in column I.
		*/
		piv = fabs ( a[jcol-1+(jcol-1)*n] );
		ipiv = jcol;
		for ( i = jcol+1; i <= n; ++i )
		{
			if ( piv < fabs ( a[i-1+(jcol-1)*n] ) )
			{
				piv = fabs ( a[i-1+(jcol-1)*n] );
				ipiv = i;
			}
		}
	
		if ( piv == 0.00 )
			return NULL;
		/*
		Switch rows JCOL and IPIV, and X.
		*/
		if ( jcol != ipiv )
		{
			for ( j = 1; j <= n; ++j )
			{
				t                 = a[jcol-1+(j-1)*n];
				a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
				a[ipiv-1+(j-1)*n] = t;
			}
			for ( j = 0; j < nb; ++j )
			{
				t            = x[jcol-1+j*n];
				x[jcol-1+j*n] = x[ipiv-1+j*n];
				x[ipiv-1+j*n] = t;
			}
		}
		/*
		Scale the pivot row.
		*/
		t = a[jcol-1+(jcol-1)*n];
		a[jcol-1+(jcol-1)*n] = 1.00;
		for ( j = jcol+1; j <= n; ++j )
			a[jcol-1+(j-1)*n] /= t;
		for ( j = 0; j < nb; ++j )
			x[jcol-1+j*n] /= t;
		/*
		Use the pivot row to eliminate lower entries in that column.
		*/
		for ( i = jcol+1; i <= n; ++i )
		{
			if ( a[i-1+(jcol-1)*n] != 0.00 )
			{
				t = - a[i-1+(jcol-1)*n];
				a[i-1+(jcol-1)*n] = 0.00;
				for ( j = jcol+1; j <= n; ++j )
					a[i-1+(j-1)*n] += t * a[jcol-1+(j-1)*n];
				for ( j = 0; j < nb; ++j )
					x[i-1+j*n] += t * x[jcol-1+j*n];
			}
		}
	}
	/*
	Back solve.
	*/
	for ( jcol = n; 2 <= jcol; --jcol )
		for ( i = 1; i < jcol; ++i )
			for ( j = 0; j < nb; ++j )
				x[i-1+j*n] -= a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
	
	return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _c8mat_copy ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_COPY copies one C8MAT to another.
  Discussion:
    A C8MAT is a doubly dimensioned array of C8's, which
    may be stored as a vector in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    03 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, double complex A[M*N], the matrix to be copied.
    Output, double complex B[M*N], the copy.
*/
{
	const dtpcxdtpcx * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	double complex * a = s_data->a1;
	const register dim_typ n = s_data->a2;
	double complex * b = s_data->a3;
	
	dim_typ i, j;
	
	for ( j = 0; j < n; ++j)
		for ( i = 0; i < m; ++i )
			b[i+j*m] = a[i+j*m];
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _c8mat_fss ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_FSS factors and solves a system with multiple right hand sides.
  Discussion:
    This routine uses partial pivoting, but no pivot vector is required.
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    05 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be positive.
    Input/output, double complex A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.
    Input, int NB, the number of right hand sides.
    Input/output, double complex X[N*NB], on input, the right hand sides of the
    linear systems.  On output, the solutions of the linear systems.
*/
{
	const dtpcxdtpcx * const s_data = data;
	const register dim_typ n = s_data->a0;
	double complex * a = s_data->a1;
	const register dim_typ nb = s_data->a2;
	double complex * x = s_data->a3;
	
	dim_typ i;
	int ipiv;
	dim_typ j;
	dim_typ jcol;
	ityp piv;
	double complex t;
	
	for ( jcol = 1; jcol <= n; ++jcol )
	{
		/*
		Find the maximum element in column I.
		*/
		piv = cabs ( a[jcol-1+(jcol-1)*n] );
		ipiv = jcol;
		for ( i = jcol+1; i <= n; ++i )
		{
			if ( piv < cabs ( a[i-1+(jcol-1)*n] ) )
			{
				piv = cabs ( a[i-1+(jcol-1)*n] );
				ipiv = i;
			}
		}
	
		if ( piv == 0.00 )
			return NULL;
	/*
		Switch rows JCOL and IPIV, and X.
		*/
		if ( jcol != ipiv )
		{
			for ( j = 1; j <= n; ++j )
			{
				t                 = a[jcol-1+(j-1)*n];
				a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
				a[ipiv-1+(j-1)*n] = t;
			}
			for ( j = 0; j < nb; ++j )
			{
				t            = x[jcol-1+j*n];
				x[jcol-1+j*n] = x[ipiv-1+j*n];
				x[ipiv-1+j*n] = t;
			}
		}
		/*
		Scale the pivot row.
		*/
		t = a[jcol-1+(jcol-1)*n];
		a[jcol-1+(jcol-1)*n] = 1.0;
		for ( j = jcol+1; j <= n; ++j )
			a[jcol-1+(j-1)*n] /=t;
		for ( j = 0; j < nb; ++j )
			x[jcol-1+j*n] /= t;
		/*
		Use the pivot row to eliminate lower entries in that column.
		*/
		for ( i = jcol+1; i <= n; i++ )
		{
			if ( a[i-1+(jcol-1)*n] != 0.00 )
			{
				t = - a[i-1+(jcol-1)*n];
				a[i-1+(jcol-1)*n] = 0.00;
				for ( j = jcol+1; j <= n; ++j)
					a[i-1+(j-1)*n] += t * a[jcol-1+(j-1)*n];
				for ( j = 0; j < nb; ++j )
					x[i-1+j*n] += t * x[jcol-1+j*n];
			}
		}
	}
	/*
	Back solve.
	*/
	for ( jcol = n; 2 <= jcol; --jcol)
		for ( i = 1; i < jcol; ++i)
			for ( j = 0; j < nb; ++j )
				x[i-1+j*n] -= a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _c8mat_fss_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_FSS_NEW factors and solves a system with multiple right hand sides.
  Discussion:
    This routine uses partial pivoting, but no pivot vector is required.
    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    02 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be positive.
    Input/output, double complex A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.
    Input, int NB, the number of right hand sides.
    Input, double complex B[N*NB], the right hand sides of the linear systems.
    Output, double complex C8MAT_FSS_NEW[N*NB], the solutions of the 
    linear systems.
*/
{
	const dtpcxdtpcx * const s_data = data;
	const register dim_typ n = s_data->a0;
	double complex * a = s_data->a1;
	const register dim_typ nb = s_data->a2;
	double complex * b = s_data->a3;
	
	dim_typ i;
	int ipiv;
	dim_typ j;
	dim_typ jcol;
	ityp piv;
	double complex t;
	double complex *x = ( double complex * ) malloc ( n * nb * sizeof ( double complex ) );
	
	for ( j = 0; j < nb; ++j)
		for ( i = 0; i < n; ++i )
			x[i+j*n] = b[i+j*n];
	for ( jcol = 1; jcol <= n; ++jcol )
	{
		/*
		Find the maximum element in column I.
		*/
		piv = cabs ( a[jcol-1+(jcol-1)*n] );
		ipiv = jcol;
		for ( i = jcol + 1; i <= n; ++i )
			if ( piv < cabs ( a[i-1+(jcol-1)*n] ) )
			{
				piv = cabs ( a[i-1+(jcol-1)*n] );
				ipiv = i;
			}
	
		if ( piv == 0.00 )
			return NULL;
		/*
		Switch rows JCOL and IPIV, and X.
		*/
		if ( jcol != ipiv )
		{
			for ( j = 1; j <= n; ++j )
			{
				t                 = a[jcol-1+(j-1)*n];
				a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
				a[ipiv-1+(j-1)*n] = t;
			}
			for ( j = 0; j < nb; ++j )
			{
				t            = x[jcol-1+j*n];
				x[jcol-1+j*n] = x[ipiv-1+j*n];
				x[ipiv-1+j*n] = t;
			}
		}
		/*
		Scale the pivot row.
		*/
		t = a[jcol-1+(jcol-1)*n];
		a[jcol-1+(jcol-1)*n] = 1.0;
		for ( j = jcol + 1; j <= n; ++j )
			a[jcol-1+(j-1)*n] /= t;
		for ( j = 0; j < nb; ++j )
			x[jcol-1+j*n] /= t;
		/*
		Use the pivot row to eliminate lower entries in that column.
		*/
		for ( i = jcol+1; i <= n; ++i )
		{
			if ( a[i-1+(jcol-1)*n] != 0.00 )
			{
				t = - a[i-1+(jcol-1)*n];
				a[i-1+(jcol-1)*n] = 0.0;
				for ( j = jcol+1; j <= n; ++j )
					a[i-1+(j-1)*n] += t * a[jcol-1+(j-1)*n];
				for ( j = 0; j < nb; ++j )
				x[i-1+j*n] += t * x[jcol-1+j*n];
			}
		}
	}
	/*
	Back solve.
	*/
	for ( jcol = n; 2 <= jcol; --jcol )
		for ( i = 1; i < jcol; ++i )
			for ( j = 0; j < nb; ++j )
				x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
	
	return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8mat_significant ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_SIGNIFICANT determines if an R8MAT is significant compared to another.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the dimension of the matrices.
    Input, double R[M*N], the vector to be compared against.
    Input, double S[M*N], the vector to be compared.
    Output, int R8MAT_SIGNIFICANT, is TRUE if S is significant
    compared to R.
*/
{
	static bool result = 2;
	
	const _2dt2pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * r = s_data->a2;
	ityp * s = s_data->a3;
	
	dim_typ i, j;
	ityp t;
	ityp tol;
	bool value = false;

	
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < m; ++i )
		{
			t = r[i+j*m] + s[i+j*m];
			tol = r8_epsilon ( ) * fabs ( r[i+j*m] );
			
			if ( tol < fabs ( r[i+j*m] - t ) )
			{
				value = true;
				break;
			}
		}
		
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8mat_minvm ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_MINVM returns C = inverse(A) * B for R8MAT's.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 December 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, the order of the matrices.
    Input, double A[N1*N1], B[N1*N2], the matrices.
    Output, double C[N1*N2], the result, C = inverse(A) * B.
*/
{
	const _2dt3pit * const s_data = data; 
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	ityp * c = s_data->a4;
	
	ityp *alu;
	ityp *d;
	
	alu = r8mat_copy_new ( n1, n1, a );
	d = r8mat_fss_new ( n1, alu, n2, b );
	r8mat_copy ( n1, n2, d, c );
	free ( alu );
	free ( d );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_add ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_ADD computes C = alpha * A + beta * B for R8MAT's.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
    Thanks to Kjartan Halvorsen for pointing out an error, 09 November 2012.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 November 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, double ALPHA, the multiplier for A.
    Input, double A[M*N], the first matrix.
    Input, double BETA, the multiplier for A.
    Input, double B[M*N], the second matrix.
    Output, double C[M*N], the sum of alpha*A+beta*B.
*/
{
	const _2itdtpitdt2pit * const s_data = data;
	
	const register ityp alpha = s_data->a0;
	const register ityp beta = s_data->a1;
	const register dim_typ m = s_data->a2;
	ityp * a = s_data->a3;
	const register dim_typ n = s_data->a4;
	ityp * b = s_data->a5;
	ityp * c = s_data->a6;
	
	dim_typ i, j;
	
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < m; ++i )
			c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_scale ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_SCALE multiplies an R8MAT by a scalar.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 December 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, double S, the scale factor.
    Input/output, double A[M*N], the matrix to be scaled.
*/
{
	const _2dtpitit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	const register ityp s = s_data->a3;
	
	dim_typ i, j;
	
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < m; ++i )
			a[i+j*m] *= s;
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _c8mat_minvm ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_MINVM returns inverse(A) * B for C8MAT's.
  Discussion:
    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, the order of the matrices.
    Input, double complex A[N1*N1], B[N1*N2], the matrices.
    Output, double complex C[N1*N2], the result, C = inverse(A) * B.
*/
{
	const _2dt3pcx * const s_data = data;
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	double complex * a = s_data->a2;
	double complex * b = s_data->a3;
	double complex * c = s_data->a4;
	
	double complex *alu;
	alu = c8mat_copy_new ( n1, n1, a );
	c8mat_copy ( n1, n2, b, c );
	c8mat_fss ( n1, alu, n2, c );
	free ( alu );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _c8mat_mm ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_MM multiplies two C8MAT's.
  Discussion:
    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
    in column-major order.
    For this routine, the result is returned as the function value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, N3, the order of the matrices.
    Input, double complex A[N1*N2], double complex B[N2*N3], 
    the matrices to multiply.
    Output, double complex C[N1*N3], the product matrix C = A * B.
*/
{
	const _3dt3pcx * const s_data = data;
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	const register dim_typ n3 = s_data->a2;
	double complex * a = s_data->a3;
	double complex * b = s_data->a4;
	double complex * c = s_data->a5;
	
	double complex *c1;
	dim_typ i, j, k;
	c1 = ( double complex * ) malloc ( n1 * n3 * sizeof ( double complex ) );
	
	for ( i = 0; i < n1; ++i )
		for ( j = 0; j < n3; ++j )
		{
			c1[i+j*n1] = 0.0;
			for ( k = 0; k < n2; ++k)
				c1[i+j*n1] += a[i+k*n1] * b[k+j*n2];
		}

	c8mat_copy ( n1, n3, c1, c );
	
	free ( c1 );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _c8mat_add_r8 ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_ADD_R8 combines two C8MAT's with real scale factors.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    03 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, double ALPHA, the first scale factor.
    Input, double complex A[M*N], the first matrix.
    Input, double BETA, the second scale factor.
    Input, double complex B[M*N], the second matrix.
    Output, double complex C[M*N], the result.
*/
{
	const _2dtitpcxit2pcx * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register ityp alpha = s_data->a2;
	double complex * a = s_data->a3;
	const register ityp beta = s_data->a4;
	double complex * b = s_data->a5;
	double complex * c = s_data->a6;
	
	dim_typ i, j;
	
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < m; ++i)
			c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _c8mat_identity_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_IDENTITY_NEW sets a C8MAT to the identity.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 July 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Output, double complex C8MAT_IDENTITY_NEW[N*N], the matrix.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
	double complex *a;
	dim_typ i, j;
	a = ( double complex * ) malloc ( n * n * sizeof ( double complex ) );
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < n; ++i )
			a[i+j*n] = i == j;
	return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _c8mat_scale_r8 ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_SCALE_R8 scales a C8MAT by a real scale factor
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    03 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, double ALPHA, the scale factor.
    Input/output, double complex A[M*N], the matrix to be scaled.
*/
{
	const _2dtitpcx * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register ityp alpha = s_data->a2;
	double complex * a = s_data->a3;
	
	dim_typ i, j;
	
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < m; ++i )
			a[i+j*m] *=  alpha;
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _c8mat_norm_li ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_NORM_LI returns the L-infinity norm of a C8MAT.
  Discussion:
    A C8MAT is a doubly dimensioned array of C8 values, stored as a vector
    in column-major order.
    The matrix L-oo norm is defined as:
      C8MAT_NORM_LI =  MAX ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
    The matrix L-oo norm is derived from the vector L-oo norm,
    and satisifies:
      c8vec_norm_li ( A * x ) <= c8mat_norm_li ( A ) * c8vec_norm_li ( x ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the order of the matrix.
    Input, double complex A[M*N], the matrix.
    Output, double C8MAT_NORM_LI, the L-infinity norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpcx * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	double complex * a = s_data->a2;
	
	dim_typ i, j;
	ityp row_sum;
	ityp value = 0.00;
	for ( i = 0; i < m; ++i)
	{
		row_sum = 0.00;
		for ( j = 0; j < n; ++j )
			row_sum = row_sum + cabs ( a[i+j*m] );
		value = MAX ( value, row_sum );
	}
	
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _c8mat_copy_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_COPY_NEW copies one C8MAT to a "new" C8MAT.
  Discussion:
    A C8MAT is a doubly dimensioned array of C8's, which
    may be stored as a vector in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    06 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, double complex A1[M*N], the matrix to be copied.
    Output, double complex C8MAT_COPY_NEW[M*N], the copy of A1.
*/
{
	const _2dtpcx * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	double complex * a1 = s_data->a2;
	
	double complex *a2;
	dim_typ i, j;
	
	a2 = ( double complex * ) malloc ( m * n * sizeof ( double complex ) );
	
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < m; ++i)
			a2[i+j*m] = a1[i+j*m];
	
	return a2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _c8mat_expm1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8MAT_EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 March 2013
  Author:
    Cleve Moler, Charles Van Loan
  Reference:
    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.
  Parameters:
    Input, int N, the dimension of the matrix.
    Input, double complex A[N*N], the matrix.
    Output, double complex C8MAT_EXPM1[N*N], the estimate for exp(A).
*/
{
	const dtpcx * const s_data = data; 
	const register dim_typ n = s_data->a0;
	double complex * a = s_data->a1;
	
    double complex *a2;
    ityp a_norm;
    ityp c;
    double complex *d;
    double complex *e;
    int ee;
    dim_typ k;
    const ityp one = 1.00;
    int p;
    const int q = 6;
    int s;
    ityp t;
    double complex *x;

    a2 = c8mat_copy_new ( n, n, a );
    a_norm = c8mat_norm_li ( n, n, a2 );
    ee = ( int ) ( log2 ( a_norm ) ) + 1;
    s = MAX ( 0, ee + 1 );
    t = 1.00 / pow ( 2.00, s );
    c8mat_scale_r8 ( n, n, t, a2 );
    x = c8mat_copy_new ( n, n, a2 );
    c = 0.50;
    e = c8mat_identity_new ( n );
    c8mat_add_r8 ( n, n, one, e, c, a2, e );
    d = c8mat_identity_new ( n );
    c8mat_add_r8 ( n, n, one, d, -c, a2, d );
    p = 1;

    for ( k = 2; k <= q; ++k )
    {
        c *= ( ityp ) ( q - k + 1 ) / ( ityp ) ( k * ( (q<<1) - k + 1 ) );

        c8mat_mm ( n, n, n, a2, x, x );
        c8mat_add_r8 ( n, n, c, x, one, e, e );

        if ( p )
            c8mat_add_r8 ( n, n, c, x, one, d, d );
        else
            c8mat_add_r8 ( n, n, -c, x, one, d, d );

        p = !p;
    }
    /*
    E -> inverse(D) * E
    */
    c8mat_minvm ( n, n, d, e, e );
    /*
    E -> E^(2*S)
    */
    for ( k = 1; k <= s; ++k )
        c8mat_mm ( n, n, n, e, e, e );

    free ( a2 );
    free ( d );
    free ( x );

    return e;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_expm1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 December 2011
  Author:
    Cleve Moler, Charles Van Loan
  Reference:
    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.
  Parameters:
    Input, int N, the dimension of the matrix.
    Input, double A[N*N], the matrix.
    Output, double R8MAT_EXPM1[N*N], the estimate for exp(A).
*/
{
	const dtpit * const s_data = data; 
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    double *a2;
    ityp a_norm;
    ityp c;
    ityp *d;
    ityp *e;
    int ee;
    dim_typ k;
    const ityp one = 1.00;
    int p;
    const dim_typ q = 6;
    int s;
    ityp t;
    ityp *x;

    a2 = r8mat_copy_new ( n, n, a );
    a_norm = r8mat_norm_li ( n, n, a2 );
    ee = ( int ) ( log2 ( a_norm ) ) + 1;
    s = MAX ( 0, ee + 1 );
    t = 1.00 / pow ( 2.00, s );
    r8mat_scale ( n, n, t, a2 );
    x = r8mat_copy_new ( n, n, a2 );
    c = 0.50;
    e = r8mat_identity_new ( n );
    r8mat_add ( n, n, one, e, c, a2, e );
    d = r8mat_identity_new ( n );
    r8mat_add ( n, n, one, d, -c, a2, d );
    p = 1;

    for ( k = 2; k <= q; ++k )
    {
        c *= ( ityp ) ( q - k + 1 ) / ( ityp ) ( k * ( (q<<1) - k + 1 ) );

        r8mat_mm ( n, n, n, a2, x, x );
        r8mat_add ( n, n, c, x, one, e, e );

        if ( p )
            r8mat_add ( n, n, c, x, one, d, d );
        else
            r8mat_add ( n, n, -c, x, one, d, d );

        p = !p;
    }
    /*
    E -> inverse(D) * E
    */
    r8mat_minvm ( n, n, d, e, e );
    /*
    E -> E^(2*S)
    */
    for ( k = 1; k <= s; ++k )
        r8mat_mm ( n, n, n, e, e, e );

    free ( a2 );
    free ( d );
    free ( x );
    return e;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_expm2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_EXPM2 uses the Taylor series for the matrix exponential.
  Discussion:
    Formally,
      exp ( A ) = I + A + 1/2 A^2 + 1/3! A^3 + ...
    This function sums the series until a tolerance is satisfied.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 December 2011
  Author:
    Cleve Moler, Charles Van Loan
  Reference:
    Cleve Moler, Charles VanLoan,
    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
    Twenty-Five Years Later,
    SIAM Review,
    Volume 45, Number 1, March 2003, pages 3-49.
  Parameters:
    Input, int N, the dimension of the matrix.
    Input, double A[N*N], the matrix.
    Output, double R8MAT_EXPM2[N*N], the estimate for exp(A).
*/
{
	const dtpit * const s_data = data; 
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp *e;
    ityp *f;
    dim_typ k = 1;
    const ityp one = 1.0;
    ityp s;

    e = r8mat_zero_new ( n, n );
    f = r8mat_identity_new ( n );

    while ( r8mat_significant ( n, n, e, f ) )
    {
        r8mat_add ( n, n, one, e, one, f, e );
        r8mat_mm ( n, n, n, a, f, f );
        s = 1.00 / ( ityp ) ( k );
        r8mat_scale ( n, n, s, f );
        ++ k ;
    }

    free ( f );
    return e;
}

#endif
