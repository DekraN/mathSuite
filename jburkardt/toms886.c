#ifndef __DISABLEDEEP_TOMS886

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dgemm ( void * data)
/******************************************************************************/
/*
  Purpose:
    DGEMM computes C = alpha * A * B and related operations.
  Discussion:
    DGEMM performs one of the matrix-matrix operations
     C := alpha * op ( A ) * op ( B ) + beta * C,
    where op ( X ) is one of
      op ( X ) = X   or   op ( X ) = X',
    ALPHA and BETA are scalars, and A, B and C are matrices, with op ( A )
    an M by K matrix, op ( B ) a K by N matrix and C an N by N matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 February 2014
  Author:
    Original FORTRAN77 version by Jack Dongarra.
    C version by John Burkardt.
  Parameters:
    Input, char TRANSA, specifies the form of op( A ) to be used in
    the matrix multiplication as follows:
    'N' or 'n', op ( A ) = A.
    'T' or 't', op ( A ) = A'.
    'C' or 'c', op ( A ) = A'.
    Input, char TRANSB, specifies the form of op ( B ) to be used in
    the matrix multiplication as follows:
    'N' or 'n', op ( B ) = B.
    'T' or 't', op ( B ) = B'.
    'C' or 'c', op ( B ) = B'.
    Input, int M, the number of rows of the  matrix op ( A ) and of the  
    matrix C.  0 <= M.
    Input, int N, the number  of columns of the matrix op ( B ) and the 
    number of columns of the matrix C.  0 <= N.
    Input, int K, the number of columns of the matrix op ( A ) and the 
    number of rows of the matrix op ( B ).  0 <= K.
    Input, double ALPHA, the scalar multiplier 
    for op ( A ) * op ( B ).
    Input, double A(LDA,KA), where:
    if TRANSA is 'N' or 'n', KA is equal to K, and the leading M by K
    part of the array contains A;
    if TRANSA is not 'N' or 'n', then KA is equal to M, and the leading
    K by M part of the array must contain the matrix A.
    Input, int LDA, the first dimension of A as declared in the calling 
    routine.  When TRANSA = 'N' or 'n' then LDA must be at least MAX ( 1, M ), 
    otherwise LDA must be at least MAX ( 1, K ).
    Input, double B(LDB,KB), where:
    if TRANSB is 'N' or 'n', kB is N, and the leading K by N 
    part of the array contains B;
    if TRANSB is not 'N' or 'n', then KB is equal to K, and the leading
    N by K part of the array must contain the matrix B.
    Input, int LDB, the first dimension of B as declared in the calling 
    routine.  When TRANSB = 'N' or 'n' then LDB must be at least MAX ( 1, K ), 
    otherwise LDB must be at least MAX ( 1, N ).
    Input, double BETA, the scalar multiplier for C.
    Input, double C(LDC,N).
    Before entry, the leading M by N part of this array must contain the 
    matrix C, except when BETA is 0.0, in which case C need not be set 
    on entry.
    On exit, the array C is overwritten by the M by N matrix
      alpha * op ( A ) * op ( B ) + beta * C.
    Input, int LDC, the first dimension of C as declared in the calling 
    routine.  MAX ( 1, M ) <= LDC.
*/
{
	const _2c3dtitpitdtpitdtitpitdt * const s_data = data;
	char transa = s_data->a0;
	char transb = s_data->a1;
	const register dim_typ m = s_data->a2;
	const register dim_typ n = s_data->a3;
	const register dim_typ k = s_data->a4;
	ityp alpha = s_data->a5;
	ityp * a = s_data->a6;
	const register dim_typ lda = s_data->a7;
	ityp * b = s_data->a8;
	const register dim_typ ldb = s_data->a9;
	ityp beta = s_data->a10;
	ityp * c = s_data->a11;
	const register dim_typ ldc = s_data->a12;
	
	dim_typ i;
	int info;
	int j;
	int l;
	int ncola;
	int nrowa;
	int nrowb;
	int nota;
	int notb;
	ityp temp;
	/*
	Set NOTA and NOTB as true if A and B respectively are not
	transposed and set NROWA, NCOLA and NROWB as the number of rows
	and columns of A and the number of rows of B respectively.
	*/
	nota = ( ( transa == 'N' ) || ( transa == 'n' ) );
	
	if ( nota )
	{
		nrowa = m;
		ncola = k;
	}
	else
	{
		nrowa = k;
		ncola = m;
	}
	
	notb = ( ( transb == 'N' ) || ( transb == 'n' ) );
	
	nrowb = notb ? k:n;
	/*
	Test the input parameters.
	*/
	info = 0;
	
	if ( ! ( transa == 'N' || transa == 'n' ||
	transa == 'C' || transa == 'c' ||
	transa == 'T' || transa == 't' ) )
		return NULL;
	
	if ( ! ( transb == 'N' || transb == 'n' ||
	transb == 'C' || transb == 'c' ||
	transb == 'T' || transb == 't' ) )
		return NULL;
	
	if ( lda < MAX ( 1, nrowa ) )
		return NULL; 
	
	if ( ldb < MAX ( 1, nrowb ) )
		return NULL;
	
	if ( ldc < MAX ( 1, m ) )
		return NULL;
	/*
	Quick return if possible.
	*/
	if ( m*n == 0)
		return NULL; 
	
	if ( ( alpha == 0.00 || k == 0 ) && ( beta == 1.00 ) )
		return NULL;
	/*
	And if alpha is 0.0.
	*/
	if ( alpha == 0.00 )
	{
		if ( beta == 0.00 )
			for ( j = 0; j < n; ++j )
				for ( i = 0; i < m; ++i)
				c[i+j*ldc] = 0.00;
		else
		{
			for ( j = 0; j < n; ++j )
				for ( i = 0; i < m; ++i )
				c[i+j*ldc] *= beta;
		}
		return NULL;
	}
	/*
	Start the operations.
	*/
	if ( notb )
	{
		/*
		Form  C := alpha*A*B + beta*C.
		*/
		if ( nota )
		{
			for ( j = 0; j < n; j++ )
			{
				if ( beta == 0.0 )
				{
					for ( i = 0; i < m; ++i )
					c[i+j*ldc] = 0.00;
				}
				else if ( beta != 1.0 )
				{
					for ( i = 0; i < m; i++ )
						c[i+j*ldc] *= beta;
				}
				
				for ( l = 0; l < k; ++l )
				{
					if ( b[l+j*ldb] != 0.0 )
					{
						temp = alpha * b[l+j*ldb];
						for ( i = 0; i < m; ++i )
						c[i+j*ldc] +=temp * a[i+l*lda];
					}
				}
	
			}
		}
		/*
		Form  C := alpha*A'*B + beta*C
		*/
		else
		{
			for ( j = 0; j < n; ++j )
			{
				for ( i = 0; i < m; ++i )
				{
					temp = 0.00;
					for ( l = 0; l < k; ++l)
						temp += a[l+i*lda] * b[l+j*ldb];
	
					c[i+j*ldc] = beta == 0.00 ? alpha*temp : alpha * temp + beta * c[i+j*ldc];
				}
			}
		}
	}
	/*
	Form  C := alpha*A*B' + beta*C
	*/
	else
	{
		if ( nota )
		{
			for ( j = 0; j < n; ++j )
			{
				if ( beta == 0.00 )
				{
					for ( i = 0; i < m; ++i )
					{
					c[i+j*ldc] = 0.00;
					}
				}
				else if ( beta != 1.00 )
				{
					for ( i = 0; i < m; i++ )
						c[i+j*ldc] *= beta;
				}
	
				for ( l = 0; l < k; ++l )
				{
					if ( b[j+l*ldb] != 0.0 )
					{
						temp = alpha * b[j+l*ldb];
						for ( i = 0; i < m; ++i )
						{
							c[i+j*ldc] += temp * a[i+l*lda];
						}
					}
				}
			}
		}	
		/*
		Form  C := alpha*A'*B' + beta*C
		*/
		else
		{
			for ( j = 0; j < n; ++j )
			{
				for ( i = 0; i < m; ++i )
				{
					temp = 0.00;
					for ( l = 0; l < k; ++l )
						temp += a[l+i*lda] * b[j+l*ldb];
					c[i+j*ldc] = beta == 0.00 ? alpha*temp : alpha * temp + beta * c[i+j*ldc];
				}
			}
		}
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cheb ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEB computes normalized Chebyshev polynomials.
  Discussion:
    This subroutine computes the array TCHEB of normalized Chebyshev 
    polynomials from degree 0 to DEG:
      T_0(x)=1, 
      T_j(x) = sqrt(2) * cos ( j * acos(x) ) 
    at the point x = PT.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 February 2014
  Author:
    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
    Marco Vianello.
    This C version by John Burkardt.
  Reference:
    Marco Caliari, Stefano de Marchi, Marco Vianello,
    Algorithm 886:
    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
    ACM Transactions on Mathematical Software,
    Volume 35, Number 3, October 2008, Article 21, 11 pages.
  Parameters:
    Input, int DEG, the degree.
    0 <= DEG.
    Input, double PT, the evaluation point.
    Output, double TCHEB[DEG+1], the value of the normalized
    Chebyshev polynomials of degrees 0 through DEG at the point PT.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ deg = s_data->a0;
	ityp * tcheb = s_data->a1;
	const register ityp pt = s_data->a2;
	
	const ityp sqrt2 = 1.4142135623730951;

	if ( deg < 0 )
		return NULL;
	
	tcheb[0] = 1.00;
	
	if ( deg < 1 )
		return NULL; 
	
	tcheb[1] = sqrt2 * pt;
	
	if ( deg < 2 )
		return NULL; 
	
	tcheb[2] = 2.00 * pt * tcheb[1] - sqrt2 * tcheb[0];
	/*
	Chebyshev recurrence.
	*/
	for (dim_typ j = 3; j <= deg; ++j )
		tcheb[j] = 2.00 * pt * tcheb[j-1] - tcheb[j-2];
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _franke ( void * data)
/******************************************************************************/
/*
  Purpose:
    FRANKE returns the value of the Franke function #1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 February 2014
  Author:
    John Burkardt
  Reference:
    Richard Franke,
    Scattered Data Interpolation: Tests of Some Methods,
    Mathematics of Computation,
    Volume 38, Number 157, January 1982, pages 181-200.
  Parameters:
    Input, double X, Y, the evalution points.
    Output, double FRANKE, the function values.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp y = a_data[1];
	
    result = 0.75 * exp (
    - ( pow ( 9.00 * x - 2.00, 2 )
    + pow ( 9.00 * y - 2.00, 2 ) ) / 4.00 )
    + 0.75 * exp (
    - ( pow ( 9.0 * x + 1.00, 2 ) ) / 49.00
    - ( 9.00 * y + 1.00 )      / 10.00 )
    + 0.50  * exp (
    - ( pow ( 9.00 * x - 7.00, 2 )
    + pow ( 9.00 * y - 3.00, 2 ) ) / 4.00 )
    - 0.20  * exp (
    - pow ( 9.00 * x - 4.00, 2 )
    - pow ( 9.00 * y - 7.00, 2 ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _padua2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    PADUA2 computes the Padua interpolation coefficient matrix.
  Discussion:
    This function computes the coefficient matrix C0, in the
    orthonormal Chebyshev basis T_j(x)T_{k-j}(y), 0 <= j <= k <= DEG,
    T_0(x)=1, T_j(x) = sqrt(2) * cos(j * acos(x)), of the
    interpolation polynomial of degree DEG of the function values FPD
    at the set of NPD Padua points (PD1,PD2) in the square [-1,1]^2.
    The interpolant may be evaluated at an arbitrary point by the
    function PD2VAL. PD1, PD2 and WPD are the Padua points and weights
    computed by PDPTS.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 February 2014
  Author:
    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi,
    Marco Vianello.
    This C version by John Burkardt.
  Reference:
    Marco Caliari, Stefano de Marchi, Marco Vianello,
    Algorithm 886:
    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
    ACM Transactions on Mathematical Software,
    Volume 35, Number 3, October 2008, Article 21, 11 pages.
  Parameters:
    Input, int DEG, the degree of approximation.
    Input, int DEGMAX, the maximum degree allowed.
    Input, int NPD, the number of Padua points.
    Input, double WPD[NPD], the weights.
    Input, double FPD[NPD], the value at the Padua points
    of the function to be interpolated.
    Workspace, double RAUX1[(DEGMAX+1)*(DEG+2)].
    Workspace, double RAUX2[(DEGMAX+1)*(DEG+2)].
    Output, double C0[(DEGMAX+2)*(DEG+1)], the coefficient matrix.
    Output, double *ESTERR, the estimated error.
*/
{
	const _3dt6pit * const s_data = data;
	const register dim_typ deg = s_data->a0;
	const register dim_typ degmax = s_data->a1;
	const register dim_typ npd = s_data->a2;
	ityp * wpd = s_data->a3;
	ityp * fpd = s_data->a4;
	ityp * raux1 = s_data->a5;
	ityp * raux2 = s_data->a6;
	ityp * c0 = s_data->a7;
	ityp * esterr = s_data->a8;
	
    ityp angle;
    dim_typ i, j, k;
    ityp pt;
    /*
    Build the matrix P_2 and store it in RAUX2.
    */
    for ( i = 0; i <= deg + 1; ++i )
    {
        angle = ( ityp ) ( i ) * M_PI / ( ityp ) ( deg + 1 );
        pt = - cos ( angle );
        cheb ( deg, pt, raux2 + i * ( degmax + 1 ) );
    }
    /*
    Build the matrix G(f) and store it in C0.
    */
    for ( j = 0; j <= deg + 1; ++j )
        for ( i = 0; i <= degmax + 1; ++i)
            c0[i+j*(degmax+2)] = 0.00;

    k = 0;
    for ( j = 0; j <= deg + 1; ++j )
        for ( i = 0; i <= deg; ++i)
            if ( ( i + j ) % 2 == 0 )
            {
                c0[i+j*(degmax+2)] = fpd[k] * wpd[k];
                ++ k;
            }
            else
                c0[i+j*(degmax+2)] = 0.00;
    /*
    Compute the matrix-matrix product G(f)*P_2' and store it in RAUX1.
    */
    dgemm ( 'n', 't', deg + 1, deg + 1, deg + 2, 1.0,
    c0, degmax + 2, raux2, degmax + 1, 0.0, raux1, degmax + 1 );
    /*
    Build the matrix P_1 and store it in RAUX2.
    */
    for ( i = 0; i <= deg; ++i)
    {
        angle = ( ityp ) ( i ) * M_PI / ( ityp ) ( deg );
        pt = - cos ( angle );
        cheb ( deg, pt, raux2 + i * ( degmax + 1 ) );
    }
    /*
    Compute the matrix-matrix product C(f) = P_1 * ( G(f) * P_2' )
    and store it in C0.
    */
    dgemm ( 'n', 'n', deg + 1, deg + 1, deg + 1, 1.00, raux2, degmax + 1, raux1, degmax + 1, 0.00, c0, degmax + 2 );
    c0[deg+0*(degmax+2)] = c0[deg+0*(degmax+2)] / 2.00;
    /*
    Estimate the error.
    */
    *esterr = 0.00;
    for ( i = 0; i <= deg - j; ++i )
        #pragma omp parallel for num_threads(3)
        for ( j = 0; j <= 2; ++j )
            *esterr += fabs ( c0[i+(deg-i-j)*(degmax+2)] );
    *esterr *= 2.00;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pd2val ( void * data)
/******************************************************************************/
/*
  Purpose:
    PD2VAL evaluates the Padua2 interpolant.
  Discussion:
    This function returns the value of the interpolant at (TG1,TG2).
    C0 is the matrix of the coefficients computed by PADUA2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 February 2014
  Author:
    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi,
    Marco Vianello.
    This C version by John Burkardt.
  Reference:
    Marco Caliari, Stefano de Marchi, Marco Vianello,
    Algorithm 886:
    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
    ACM Transactions on Mathematical Software,
    Volume 35, Number 3, October 2008, Article 21, 11 pages.
  Parameters:
 int DEG, the degree of approximation.
    Input, int DEGMAX, the maximum degree allowed.
    Input, double C0[(0:DEGMAX+1)*(0:DEG)], the coefficient matrix.
    Input, double TG1, TG2, the first and second coordinates of
    the target point.
    Output, double PD2VAL, the value of the interpolant at
    the target point.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2itdtpit * const s_data = data;
	
	const register dim_typ deg = s_data->a0;
	const register ityp tg1 = s_data->a1;
	const register ityp tg2 = s_data->a2;
	const register dim_typ degmax = s_data->a3;
	ityp * c0 = s_data->a4;
	
    dim_typ i, j;
    ityp t;
    ityp *ttg1;
    ityp *ttg2;
    ityp value;
    /*
    Compute the normalized Chebyshev polynomials at the target point.
    */
    ttg1 = ( ityp * ) malloc ( ( deg + 1 ) * sizeof ( ityp ) );
    cheb ( deg, tg1, ttg1 );

    ttg2 = ( ityp * ) malloc ( ( deg + 1 ) * sizeof ( ityp ) );
    cheb ( deg, tg2, ttg2 );
    /*
    Evaluate the interpolant
    */
    value = 0.00;
    for ( i = deg; 0 <= i; --i )
    {
        t = 0.00;
        for ( j = 0; j <= i; ++j )
            t += ttg1[j] * c0[j+(deg-i)*(degmax+2)];
        value += ttg2[deg-i] * t;
    }

    free ( ttg1 );
    free ( ttg2 );

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _pdpts ( void * data)
/******************************************************************************/
/*
  Purpose:
    PDPTS returns the points and weights for Padua interpolation.
  Discussion:
    This subroutine computes the first family of Padua points and
    weights corresponding to degree DEG.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 February 2014
  Author:
    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi,
    Marco Vianello.
    This C version by John Burkardt.
  Reference:
    Marco Caliari, Stefano de Marchi, Marco Vianello,
    Algorithm 886:
    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
    ACM Transactions on Mathematical Software,
    Volume 35, Number 3, October 2008, Article 21, 11 pages.
  Parameters:
    Input, int DEG, the degree of approximation.
    Output, double PD1[NPD], PD2[NPD], the first and second
    coordinates of the Padua points
    Output, double WPD[NPD], the weights.
    Output, int *NPD, the number of Padua points.
    NPD = ( DEG + 1 ) * ( DEG + 2 ) / 2.
*/
{
	const dt3pitpi * const s_data = data;
	const register dim_typ deg = s_data->a0;
	ityp * pd1 = s_data->a1;
	ityp * pd2 = s_data->a2;
	ityp * wpd = s_data->a3;
	int * npd = s_data->a4;
	
    dim_typ itemp0;
    dim_typ j, k;
    ityp rtemp0;
    /*
    Compute the Padua points of the first family at degree DEG.
    */
    if ( deg == 0 )
    {
        pd1[0] = pd2[0] = -1.00;
        wpd[0] = 2.00;
        *npd = 1;
        return NULL;
    }

    *npd = 0;
    itemp0 = deg * ( deg + 1 );
    rtemp0 = M_PI / ( ityp ) ( itemp0 );

    for ( j = 0; j <= deg + 1; ++j )
        for ( k = ( j % 2 ); k <= deg; k += 2 )
        {
            pd1[*npd] = - cos ( ( ityp ) ( ( deg + 1 ) * k ) * rtemp0 );
            pd2[*npd] = - cos ( ( ityp ) ( deg * j ) * rtemp0 );
            wpd[*npd] = 2.0 / ( ityp ) ( itemp0 );

            if ( k == 0 || k == deg || j == 0 || j == deg + 1 )
                wpd[*npd] /= 2.00;
            ++ *npd;
        }

    return NULL;
}

#endif
