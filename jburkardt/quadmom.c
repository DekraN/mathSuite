#ifndef __DISABLEDEEP_QUADMOM

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *    _moment_method ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENT_METHOD computes a quadrature rule by the method of moments.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 September 2013
  Author:
    John Burkardt
  Reference:
    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.
  Parameters:
    Input, int N, the order of the quadrature rule.
    Input, double MOMENT[2*N+1], moments 0 through 2*N.
    Output, double X[N], W[N], the points and weights of the quadrature rule.
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * moment = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
    ityp *alpha;
    ityp *beta;
    bool flag;
    ityp *h;
    dim_typ i;
    dim_typ it_max;
    int it_num;
    dim_typ j;
    ityp *jacobi;
    ityp *r;
    int rot_num;
    ityp *v;


    /*
    Define the N+1 by N+1 Hankel matrix H(I,J) = moment(I+J).
    */
    h = ( ityp * ) malloc ( ( n + 1 ) * ( n + 1 ) * sizeof ( ityp ) );

    for ( i = 0; i <= n; ++i )
        for ( j = 0; j <= n; ++j )
            h[i+j*(n+1)] = moment[i+j];

    /*
    Compute R, the upper triangular Cholesky factor of H.
    */
    r = r8mat_cholesky_factor_upper ( n + 1, h, &flag );

    if ( flag != 0 )
        return NULL;
    /*
    Compute ALPHA and BETA from R, using Golub and Welsch's formula.
    */
    alpha = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    alpha[0] = r[0+1*(n+1)] / r[0+0*(n+1)];
    for ( i = 1; i < n; ++i )
        alpha[i] = r[i+(i+1)*(n+1)] / r[i+i*(n+1)] - r[i-1+i*(n+1)] / r[i-1+(i-1)*(n+1)];
    beta = ( ityp * ) malloc ( ( n - 1 ) * sizeof ( ityp ) );

    for ( i = 0; i < n - 1; ++i )
        beta[i] = r[i+1+(i+1)*(n+1)] / r[i+i*(n+1)];
    /*
    Compute the points and weights from the moments.
    */
    jacobi = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            jacobi[i+j*n] = 0.00;

    for ( i = 0; i < n; ++i )
        jacobi[i+i*n] = alpha[i];

    for ( i = 0; i < n - 1; ++i )
    {
        jacobi[i+(i+1)*n] = beta[i];
        jacobi[i+1+i*n] = beta[i];
    }
    /*
    Get the eigendecomposition of the Jacobi matrix.
    */
    it_max = 100;
    v = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    jacobi_eigenvalue ( n, jacobi, it_max, v, x, &it_num, &rot_num );

    for ( i = 0; i < n; ++i )
        w[i] = moment[0] * pow ( v[0+i*n], 2 );

    free ( alpha );
    free ( beta );
    free ( h );
    free ( jacobi );
    free ( r );
    free ( v );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _moments_laguerre ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENTS_LAGUERRE returns moments of the Laguerre distribution.
  Discussion:
    pdf(x) = exp ( -x )
    mu(k) = integral ( 0 <= x < +oo ) x^k pdf(x) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of moments desired.
    Output, double W(0:M-1), the weighted integrals of X^0
    through X^(M-1).
*/
{
	const register dim_typ m = *(dim_typ *) data;
	
	ityp *w = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    for (dim_typ k = 0; k < m; ++k)
        w[k] = r8_factorial ( k );
    return w;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _moments_legendre ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENTS_LEGENDRE returns moments of the Legendre weight on [A,B].
  Discussion:
    mu(k) = integral ( a <= x <= b ) x^k dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of moments desired.
    Input, double A, B, the left and right endpoints
    of the interval.
    Output, double W(0:M-1), the weighted integrals of X^0
    through X^(M-1).
*/
{
	const _2itdt * const s_data = data;
	
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ m = s_data->a2;
	
    ityp ak;
    ityp bk;

    ityp *w = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    bk = bk = 1.00;
    for (dim_typ k = 0; k < m; ++k )
    {
        bk *= b;
        ak *=  a;
        w[k] = ( bk - ak ) / ( ityp ) ( k + 1 );
    }

    return w;
}



/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _moments_normal ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENTS_NORMAL returns moments of the standard Normal distribution.
  Discussion:
    pdf(x) = exp ( -((x-mu)/sigma)^2/2 ) / sigma / sqrt ( M_PI * 2 )
    mu(k) = integral ( -oo < x < +oo ) x^k pdf(x) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of moments desired.
    Input, double MU, SIGMA, the mean and standard deviation.
    Output, double W(0:M-1), the weighted integrals of X^0
    through X^(M-1).
*/
{
	const _2itdt * const s_data = data;
	
	const register ityp mu = s_data->a0;
	const register ityp sigma = s_data->a1;
	const register dim_typ m = s_data->a2;
	
    int j;
    int j_hi;
    int k;
    double t;
    double *w;

    w = ( double * ) malloc ( m * sizeof ( double ) );

    for ( k = 0; k < m; ++k )
    {
        t = 0.00;
        j_hi = k / 2;
        for ( j = 0; j <= j_hi; ++j )
            t = t + r8_choose ( k, j<<1 ) * r8_factorial2 ( (j<<1) - 1 )* pow ( sigma, (j<<1) ) * pow ( mu, k - (j<<1) );
        w[k] = t;
    }

    return w;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _moments_truncated_normal_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENTS_TRUNCATED_NORMAL_AB: moments of the truncated Normal distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of moments desired.
    Input, double MU, SIGMA, the mean and standard deviation.
    Input, double A, B, the lower and upper truncation limits.
    Output, double W(0:M-1), the weighted integrals of X^0
    through X^(M-1).
*/
{
	const dt4it * const s_data = data;
	const register dim_typ m = s_data->a0;
	ityp mu = s_data->a1;
	ityp sigma = s_data->a2;
	ityp a = s_data->a3;
	ityp b = s_data->a4;
	
    ityp *w = ( double * ) malloc ( m * sizeof ( double ) );
    for (dim_typ order = 0; order < m; ++order )
        w[order] = truncated_normal_ab_moment ( order, mu, sigma, a, b );

    return w;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _moments_truncated_normal_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENTS_TRUNCATED_NORMAL_A: moments of the lower truncated Normal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of moments desired.
    Input, double MU, SIGMA, the mean and standard deviation.
    Input, double A, the lower truncation limit.
    Output, double W(0:M-1), the weighted integrals of X^0
    through X^(M-1).
*/
{
	const dt3it * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register ityp mu = s_data->a1;
	const register ityp sigma = s_data->a2;
	const register ityp a = s_data->a3;
	
    ityp *w = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    for (dim_typ order = 0; order < m; ++order )
        w[order] = truncated_normal_a_moment ( order, mu, sigma, a );
    return w;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _moments_truncated_normal_b ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENTS_TRUNCATED_NORMAL_B: moments of the upper truncated Normal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of moments desired.
    Input, double MU, SIGMA, the mean and standard deviation.
    Input, double B, the upper truncation limit.
    Output, double W(0:M-1), the weighted integrals of X^0
    through X^(M-1).
*/
{
	const dt3it * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register ityp mu = s_data->a1;
	const register ityp sigma = s_data->a2;
	const register ityp b = s_data->a3;
	
    ityp *w = ( double * ) malloc ( m * sizeof ( double ) );
    for (dim_typ order = 0; order < m; ++order )
        w[order] = truncated_normal_b_moment ( order, mu, sigma, b );
    return w;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _normal_01_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    NORMAL_01_CDF evaluates the Normal 01 CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 February 1999
  Author:
    John Burkardt
  Reference:
    A G Adams,
    Areas Under the Normal Curve,
    Algorithm 39,
    Computer j.,
    Volume 12, pages 197-198, 1969.
  Parameters:
    Input, double X, the argument of the CDF.
    Output, double CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	ityp a1 = 0.398942280444;
	ityp a2 = 0.399903438504;
	ityp a3 = 5.75885480458;
	ityp a4 = 29.8213557808;
	ityp a5 = 2.62433121679;
	ityp a6 = 48.6959930692;
	ityp a7 = 5.92885724438;
	ityp b0 = 0.398942280385;
	ityp b1 = 3.8052E-08;
	ityp b2 = 1.00000615302;
	ityp b3 = 3.98064794E-04;
	ityp b4 = 1.98615381364;
	ityp b5 = 0.151679116635;
	ityp b6 = 5.29330324926;
	ityp b7 = 4.8385912808;
	ityp b8 = 15.1508972451;
	ityp b9 = 0.742380924027;
	ityp b10 = 30.789933034;
	ityp b11 = 3.99019417011;
	ityp cdf;
	ityp q;
	ityp y;
	/*
	|X| <= 1.28.
	*/
	if ( fabs ( x ) <= 1.28 )
	{
		y = 0.50 * x * x;
		q = 0.50 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5+ a6 / ( y + a7 ) ) ) );
		/*
		1.28 < |X| <= 12.7
		*/
	}
	else if ( fabs ( x ) <= 12.70 )
	{
		y = 0.50 * x * x;
		q = exp ( - y ) * b0 / ( fabs ( x ) - b1+ b2 / ( fabs ( x ) + b3+ b4 / ( fabs ( x ) - b5+ b6 / ( fabs ( x ) + b7- b8 / ( fabs ( x ) + b9+ b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
		/*
		12.7 < |X|
		*/
	}
	else
		q = 0.00;
	/*
	Take account of negative X.
	*/
	
	result = x<0.00 ? q : 1.00-q;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _normal_01_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    NORMAL_01_PDF evaluates the Normal 01 PDF.
  Discussion:
    The Normal 01 PDF is also called the "Standard Normal" PDF, or
    the Normal PDF with 0 mean and variance 1.
    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 September 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the PDF.
    Output, double PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	result = exp ( -0.50 * x * x ) / sqrt ( M_2TPI );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _truncated_normal_ab_moment ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_AB_MOMENT: moments of the truncated Normal PDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 September 2013
  Author:
    John Burkardt
  Reference:
    Phoebus Dhrymes,
    Moments of Truncated Normal Distributions,
    May 2005.
  Parameters:
    Input, int ORDER, the order of the moment.
    0 <= ORDER.
    Input, double MU, S, the mean and standard deviation of the
    parent Normal distribution.
    0.0 < S.
    Input, double A, B, the lower and upper truncation limits.
    A < B.
    Output, double TRUNCATED_NORMAL_AB_MOMENT, the moment of the PDF.
*/
{
	static ityp result = MAX_VAL;	
	
	const dt4it * const s_data = data;
	const register dim_typ order = s_data->a0;
	ityp mu = s_data->a1;
	ityp s = s_data->a2;
	ityp a = s_data->a3;
	ityp b = s_data->a4;
	
    ityp a_h;
    ityp a_cdf;
    ityp a_pdf;
    ityp b_h;
    ityp b_cdf;
    ityp b_pdf;
    ityp ir;
    ityp irm1;
    ityp irm2;
    ityp moment;
    dim_typ r;

    if ( s <= 0.0 || b <= a )
    {
    	result = MAX_VAL;
        return &result;
    }

    a_h = ( a - mu ) / s;
    a_pdf = normal_01_pdf ( a_h );
    a_cdf = normal_01_cdf ( a_h );

    if ( a_cdf == 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    b_h = ( b - mu ) / s;
    b_pdf = normal_01_pdf ( b_h );
    b_cdf = normal_01_cdf ( b_h );

    if ( b_cdf == 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    moment = irm2 = irm1 = 0.00;

    for ( r = 0; r <= order; ++r )
    {
        if ( r == 0 )
        ir = 1.00;
        else if ( r == 1 )
            ir = - ( b_pdf - a_pdf ) / ( b_cdf - a_cdf );
        else
            ir = ( ityp ) ( r - 1 ) * irm2- ( pow ( b_h, r - 1 ) * b_pdf - pow ( a_h, r - 1 ) * a_pdf )/ ( b_cdf - a_cdf );

        moment = moment + r8_choose ( order, r ) * pow ( mu, order - r )* pow ( s, r ) * ir;

        irm2 = irm1;
        irm1 = ir;
    }
	
	result = moment; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _truncated_normal_a_moment ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_A_MOMENT: moments of the lower truncated Normal PDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 September 2013
  Author:
    John Burkardt
  Reference:
    Phoebus Dhrymes,
    Moments of Truncated Normal Distributions,
    May 2005.
  Parameters:
    Input, int ORDER, the order of the moment.
    0 <= ORDER.
    Input, double MU, S, the mean and standard deviation of the
    parent Normal distribution.
    Input, double A, the lower truncation limit.
    Output, double TRUNCATED_NORMAL_A_MOMENT, the moment of the PDF.
*/
{
	static ityp result = MAX_VAL;	
	
	const dt3it * const s_data = data;
	const register dim_typ order = s_data->a0;
	const register ityp mu = s_data->a1;
	const register ityp s = s_data->a2;
	const register ityp a = s_data->a3;
	
	result = r8_mop ( order )* truncated_normal_b_moment ( order, - mu, s, - a );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _truncated_normal_b_moment ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_B_MOMENT: moments of the upper truncated Normal PDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 September 2013
  Author:
    John Burkardt
  Reference:
    Phoebus Dhrymes,
    Moments of Truncated Normal Distributions,
    May 2005.
  Parameters:
    Input, int ORDER, the order of the moment.
    0 <= ORDER.
    Input, double MU, S, the mean and standard deviation of the
    parent Normal distribution.
    Input, double B, the upper truncation limit.
    Output, double TRUNCATED_NORMAL_B_MOMENT, the moment of the PDF.
*/
{
	static ityp result = MAX_VAL;	
	
	const dt3it * const s_data = data;
	const register dim_typ order = s_data->a0;
	const register ityp mu = s_data->a1;
	const register ityp s = s_data->a2;
	const register ityp b = s_data->a3;
	
    ityp f;
    ityp h;
    ityp h_cdf;
    ityp h_pdf;
    ityp ir;
    ityp irm1;
    ityp irm2;
    ityp moment;
    dim_typ r;

    if ( order < 0 )
    {
    	result = MAX_VAL;
        return &result;
    }

    h = ( b - mu ) / s;
    h_pdf = normal_01_pdf ( h );
    h_cdf = normal_01_cdf ( h );

    if ( h_cdf == 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    f = h_pdf / h_cdf;

    moment = irm2 = irm1 = 0.00;

    for ( r = 0; r <= order; ++r )
    {
        if ( r == 0 )
            ir = 1.00;
        else if ( r == 1 )
            ir = - f;
        else
        ir = - pow ( h, r - 1 ) * f + ( ityp ) ( r - 1 ) * irm2;

        moment = moment + r8_choose ( order, r ) * pow ( mu, order - r )* pow ( s, r ) * ir;

        irm2 = irm1;
        irm1 = ir;
    }

	result = moment;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_cholesky_factor_upper ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_CHOLESKY_FACTOR_UPPER: upper Cholesky factor of a symmetric R8MAT.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
    The matrix must be symmetric and positive semidefinite.
    For a positive semidefinite symmetric matrix A, the Cholesky factorization
    is an upper triangular matrix R such that:
      A = R' * R
    Note that the usual Cholesky factor is a LOWER triangular matrix L
    such that
      A = L * L'
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of rows and columns of the matrix A.
    Input, double A[N*N], the N by N matrix.
    Output, int *FLAG, an error flag.
    0, no error was detected.
    1, the matrix was not positive definite.  A NULL factor was returned.
    Output, double R8MAT_CHOLESKY_FACTOR_UPPER[N*N], the N by N upper triangular
    "Choresky" factor.
*/
{
	const dtpitpb * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	bool * flag = s_data->a2;
	
    ityp *c;
    dim_typ i, j, k;
    ityp sum2;

    *flag = 0;

    c = r8mat_copy_new ( n, n, a );

    for ( j = 0; j < n; ++j )
    {
        for ( i = 0; i < j; ++i)
            c[j+i*n] = 0.0;
        for ( i = j; i < n; i++ )
        {
            sum2 = c[i+j*n];
            for ( k = 0; k < j; ++k )
                sum2 -= c[k+j*n] * c[k+i*n];
            if ( i == j )
            {
                if ( sum2 <= 0.0 )
                {
                    *flag = 1;
                    return NULL;
                }
                c[j+i*n] = sqrt ( sum2 );
            }
            else
            {
                c[j+i*n] = 0.00 + (c[j+j*n] != 0.00)*sum2 / c[j+j*n];
            }
        }
    }

  	return c;
}

#endif
