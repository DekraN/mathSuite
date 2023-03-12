#ifndef __DISABLEDEEP_RBFINTERP1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _phi1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    PHI1 evaluates the multiquadric radial basis function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2012
  Author:
    John Burkardt
  Reference:
    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    Third Edition,
    Cambridge University Press, 2007,
    ISBN13: 978-0-521-88068-8,
    LC: QA297.N866.
  Parameters:
    Input, int N, the number of points.
    Input, double R[N], the radial separation.
    0 < R.
    Input, double R0, a scale factor.
    Output, double V[N], the value of the radial basis function.
*/
{
	const dtit2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp r0 = s_data->a1;
	ityp * r = s_data->a2;
	ityp * v = s_data->a3;
	
    for (dim_typ i = 0; i < n; ++i )
        v[i] = sqrt ( r[i] * r[i] + r0 * r0 );
        
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _phi2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    PHI2 evaluates the inverse multiquadric radial basis function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2012
  Author:
    John Burkardt
  Reference:
    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    Third Edition,
    Cambridge University Press, 2007,
    ISBN13: 978-0-521-88068-8,
    LC: QA297.N866.
  Parameters:
    Input, int N, the number of points.
    Input, double R[N], the radial separation.
    0 < R.
    Input, double R0, a scale factor.
    Output, double V[N], the value of the radial basis function.
*/
{
	const dtit2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp r0 = s_data->a1;
	ityp * r = s_data->a2;
	ityp * v = s_data->a3;
	
    for (dim_typ i = 0; i < n; ++i )
        v[i] = 1.00 / sqrt ( r[i] * r[i] + r0 * r0 );
        
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _phi3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    PHI3 evaluates the thin-plate spline radial basis function.
  Discussion:
    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
    it may be desirable to choose a value of R0 smaller than any possible R.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2012
  Author:
    John Burkardt
  Reference:
    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    Third Edition,
    Cambridge University Press, 2007,
    ISBN13: 978-0-521-88068-8,
    LC: QA297.N866.
  Parameters:
    Input, int N, the number of points.
    Input, double R[N], the radial separation.
    0 < R.
    Input, double R0, a scale factor.
    Output, double V[N], the value of the radial basis function.
*/
{
	const dtit2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp r0 = s_data->a1;
	ityp * r = s_data->a2;
	ityp * v = s_data->a3;
	
    for (dim_typ i = 0; i < n; ++i )
        v[i] = 0.00 + (r[i]>0.00)*r[i] * r[i] * log ( r[i] / r0 );
        
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _phi4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    PHI4 evaluates the gaussian radial basis function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2012
  Author:
    John Burkardt
  Reference:
    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    Third Edition,
    Cambridge University Press, 2007,
    ISBN13: 978-0-521-88068-8,
    LC: QA297.N866.
  Parameters:
    Input, int N, the number of points.
    Input, double R[N], the radial separation.
    0 < R.
    Input, double R0, a scale factor.
    Output, double V[N], the value of the radial basis function.
*/
{
	const dtit2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp r0 = s_data->a1;
	ityp * r = s_data->a2;
	ityp * v = s_data->a3;
	
    for (dim_typ i = 0; i < n; ++i )
        v[i] = exp ( - 0.50 * r[i] * r[i] / r0 / r0 );
        
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_solve_svd ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_SOLVE_SVD solves a linear system A*x=b using the SVD.
  Discussion:
    When the system is determined, the solution is the solution in the
    ordinary sense, and A*x = b.
    When the system is overdetermined, the solution minimizes the
    L2 norm of the residual ||A*x-b||.
    When the system is underdetermined, ||A*x-b|| should be zero, and
    the solution is the solution of minimum L2 norm, ||x||.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns
    in the matrix A.
    Input, double A[M,*N], the matrix.
    Input, double B[M], the right hand side.
    Output, double R8MAT_SOLVE_SVD[N], the solution.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	
    ityp *a_copy;
    ityp *a_pseudo;
    ityp *e;
    dim_typ i, j, k, l;
    int info;
    int lda;
    int ldu;
    int ldv;
    int job;
    int lwork;
    ityp *s;
    ityp *sp;
    ityp *sdiag;
    ityp *u;
    ityp *v;
    ityp *work;
    ityp *x;
    /*
    Compute the SVD decomposition.
    */
    a_copy = r8mat_copy_new ( m, n, a );
    lda = m;
    sdiag = ( ityp * ) malloc ( MAX ( m + 1, n ) * sizeof ( ityp ) );
    e = ( ityp * ) malloc ( MAX ( m + 1, n ) * sizeof ( ityp ) );
    u = ( ityp * ) malloc ( m * m * sizeof ( ityp ) );
    ldu = m;
    v = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
    ldv = n;
    work = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    job = 11;

    info = dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job );

    if ( info != 0 )
        return NULL;

    s = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
    for ( i = 0; i < m; ++i)
        s[i+j*m] = sp[i+j*n] = 0.00;

    for ( i = 0; i < MIN ( m, n ); ++i )
        s[i+i*m] = sdiag[i];

    /*
    Compute the pseudo inverse.
    */
    sp = ( ityp * ) malloc ( n * m * sizeof ( ityp ) );

    for ( i = 0; i < MIN ( m, n ); ++i )
        if ( s[i+i*m] != 0.00 )
            sp[i+i*n] = 1.00 / s[i+i*m];

    a_pseudo = ( ityp * ) malloc ( n * m * sizeof ( ityp ) );

    for ( j = 0; j < m; ++j )
        for ( i = 0; i < n; ++i )
        {
            a_pseudo[i+j*n] = 0.0;
            for ( k = 0; k < n; ++k)
                for ( l = 0; l < m; ++l )
                    a_pseudo[i+j*n] = a_pseudo[i+j*n] + v[i+k*n] * sp[k+l*n] * u[j+l*m];
        }
    /*
    Compute x = A_pseudo * b.
    */
    x = r8mat_mv_new ( n, m, a_pseudo, b );

    free ( a_copy );
    free ( a_pseudo );
    free ( e );
    free ( s );
    free ( sdiag );
    free ( sp );
    free ( u );
    free ( v );
    free ( work );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   *rbf_interp ( const register dim_typ m, const register dim_typ nd, ityp xd[static m*nd], const register ityp r0,
  void phix ( int n, ityp r[], ityp r0, ityp v[] ), ityp w[static nd], const register dim_typ ni, ityp xi[static m*ni] )
/******************************************************************************/
/*
  Purpose:
    RBF_INTERP evaluates a radial basis function interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2012
  Author:
    John Burkardt
  Reference:
    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    Third Edition,
    Cambridge University Press, 2007,
    ISBN13: 978-0-521-88068-8,
    LC: QA297.N866.
  Parameters:
    Input, int M, the spatial dimension.
    Input, int ND, the number of data points.
    Input, double XD[M*ND], the data points.
    Input, double R0, a scale factor.  R0 should be larger than the typical
    separation between points, but smaller than the maximum separation.
    The value of R0 has a significant effect on the resulting interpolant.
    Input, void PHI ( int N, double R[], double R0, double V[] ), a
    function to evaluate the radial basis functions.
    Input, double W[ND], the weights, as computed by RBF_WEIGHTS.
    Input, int NI, the number of interpolation points.
    Input, double XI[M*NI], the interpolation points.
    Output, double RBF_INTERP[NI], the interpolated values.
*/
{
    ityp *fi;
    dim_typ i, j, k;
    ityp *r;
    ityp *v;

    fi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );
    r = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    v = ( ityp * ) malloc ( nd * sizeof ( ityp ) );

    for ( i = 0; i < ni; ++i )
    {
        for ( j = 0; j < nd; ++j )
        {
            r[j] = 0.00;
            for ( k = 0; k < m; ++k )
                r[j] = r[j] + pow ( xi[k+i*m] - xd[k+j*m], 2 );
            r[j] = sqrt ( r[j] );
        }
        phix ( nd, r, r0, v );
        fi[i] = r8vec_dot_product ( nd, v, w );
    }

    free ( r );
    free ( v );

    return fi;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   *rbf_weight ( const register dim_typ m, const register dim_typ nd, ityp xd[static m*nd], const register ityp r0,void phix ( int n, ityp r[], ityp r0, ityp v[] ),ityp fd[static nd] )
/******************************************************************************/
/*
  Purpose:
    RBF_WEIGHT computes weights for radial basis function interpolation.
  Discussion:
    We assume that there are N (nonsingular) equations in N unknowns.
    However, it should be clear that, if we are willing to do some kind
    of least squares calculation, we could allow for singularity,
    inconsistency, or underdetermine systems.  This could be associated
    with data points that are very close or repeated, a smaller number
    of data points than function values, or some other ill-conditioning
    of the system arising from a peculiarity in the point spacing.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2012
  Author:
    John Burkardt
  Reference:
    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    Third Edition,
    Cambridge University Press, 2007,
    ISBN13: 978-0-521-88068-8,
    LC: QA297.N866.
  Parameters:
    Input, int M, the spatial dimension.
    Input, int ND, the number of data points.
    Input, double XD[M*ND], the data points.
    Input, double R0, a scale factor.  R0 should be larger than the typical
    separation between points, but smaller than the maximum separation.
    The value of R0 has a significant effect on the resulting interpolant.
    Input, void PHI ( int N, double R[], double R0, double V[] ), a
    function to evaluate the radial basis functions.
    Input, double FD[ND], the function values at the data points.
    Output, double RBF_WEIGHT[ND], the weights.
*/
{
    ityp *a;
    dim_typ i, j, k;
    ityp *r;
    ityp *v;
    ityp *w;

    a = ( ityp * ) malloc ( nd * nd * sizeof ( ityp ) );
    r = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    v = ( ityp * ) malloc ( nd * sizeof ( ityp ) );

    for ( i = 0; i < nd; ++i )
    {
        for ( j = 0; j < nd; ++j )
        {
            r[j] = 0.0;
            for ( k = 0; k < m; ++k )
                r[j] += pow ( xd[k+i*m] - xd[k+j*m], 2 );
            r[j] = sqrt ( r[j] );
        }
        phix ( nd, r, r0, v );

        for ( j = 0; j < nd; ++j)
            a[i+j*nd] = v[j];
    }
    /*
    Solve for the weights.
    */
    w = r8mat_solve_svd ( nd, nd, a, fd );

    free ( a );
    free ( r );
    free ( v );

    return w;
}

#endif
