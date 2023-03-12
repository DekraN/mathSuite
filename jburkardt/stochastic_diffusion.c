#ifndef __DISABLEDEEP_STOCHASTICDIFFUSION

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _diffusivity_1d_xk ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFUSIVITY_1D_XK evaluates a 1D stochastic diffusivity function.
  Discussion:
    The 1D diffusion equation has the form
      - d/dx ( DC(X) Del U(X) ) = F(X)
    where DC(X) is a function called the diffusivity.
    In the stochastic version of the problem, the diffusivity function
    includes the influence of stochastic parameters:
      - d/dx ( DC(X;OMEGA) d/dx U(X) ) = F(X).
    In this function, the domain is assumed to be the unit interval [0.1].
    For DC0 = 1 and F(X) = 0, with boundary conditions U(0:OMEGA) = 0,
    U(1;OMEGA) = 1, the exact solution is
    If OMEGA ~= 0:
      U(X;OMEGA) = log ( 1 + OMEGA * X ) / log ( 1 + OMEGA )
    If OMEGA = 0:
      U(X;OMEGA) = X
    In the numerical experiments described in the paper, OMEGA was taken
    to be a random variable with a Beta, or Uniform, or Gaussian or
    Poisson or Binomial distribution.
    For the Gaussian and Poisson distributions, the positivity requirement
    could not be guaranteed, and the experiments were simply made with a
    "small" variance of 0.1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 August 2013
  Author:
    John Burkardt
  Reference:
    Dongbin Xiu, George Karniadakis,
    Modeling uncertainty in steady state diffusion problems via
    generalized polynomial chaos,
    Computer Methods in Applied Mechanics and Engineering,
    Volume 191, 2002, pages 4927-4948.
  Parameters:
    Input, double DC0, the constant term in the expansion of the
    diffusion coefficient.
    Input, int M, the number of stochastic parameters.
    Input, double OMEGA[M], the stochastic parameters.
    Input, int N, the number of evaluation points.
    Input, double X[N], the point where the diffusion coefficient
    is to be evaluated.
    Output, double DIFFUSIVITY_1D_XK[N], the value of the diffusion coefficient
    at X.
*/
{
	const it2dt2pit * const s_data = data;
	
	const register ityp dc0 = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * omega = s_data->a3;
	ityp * x = s_data->a4;
	
    ityp *dc;
    dim_typ j, k;
    ityp w;

    k = 0;
    w = 1.00;

    dc = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        dc[j] = 0.00;

    while ( k < m )
    {
        if ( k < m )
        {
            ++ k;
            for ( j = 0; j < n; ++j)
                dc[j] += omega[k-1] * sin ( w * M_PI * x[j] );
        }

        if ( k < m )
        {
            ++ k;
            for ( j = 0; j < n; ++j )
                dc[j] += omega[k-1] * cos ( w * M_PI * x[j] );
        }

        ++ w;

    }

    for ( j = 0; j < n; ++j )
        dc[j] = dc0 + exp ( exp ( - 0.125 ) * dc[j] );

    return dc;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _diffusivity_2d_bnt ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFUSIVITY_2D_BNT evaluates a 2D stochastic diffusivity function.
  Discussion:
    The 2D diffusion equation has the form
      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
    where DC(X,Y) is a function called the diffusivity.
    In the stochastic version of the problem, the diffusivity function
    includes the influence of stochastic parameters:
      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
    In this function, the domain is the rectangle [-1.5,0]x[-0.4,0.8].
    The four stochastic parameters OMEGA(1:4) are assumed to be independent
    identically distributed random variables with mean value zero and
    variance 1.  The distribution is typically taken to be Gaussian or
    uniform.
    A collocation approach to this problem would then use the roots of
    Hermite or Legendre polynomials.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 August 2013
  Author:
    John Burkardt
  Reference:
    Ivo Babuska, Fabio Nobile, Raul Tempone,
    A stochastic collocation method for elliptic partial differential equations
    with random input data,
    SIAM Journal on Numerical Analysis,
    Volume 45, Number 3, 2007, pages 1005-1034.
  Parameters:
    Input, double DC0, the constant term in the expansion of the
    diffusion coefficient.  Take DC0 = 10.
    Input, double OMEGA[4], the stochastic parameters.
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the points where the diffusion
    coefficient is to be evaluated.
    Output, double DIFFUSIVITY_2D_BNT[N], the value of the diffusion coefficient
    at (X,Y).
*/
{
	const dtit3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp dc0 = s_data->a1;
	ityp * x = s_data->a2;
	ityp * y = s_data->a3;
	ityp * omega =  s_data->a4;
	
    ityp *arg;
    ityp *dc;
    dim_typ j;

    arg = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dc = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( j = 0; j < n; ++j)
        dc[j] = dc0 + exp ( exp ( - 0.125 ) * omega[0] * cos ( M_PI * x[j] )+ omega[1] * sin ( M_PI * x[j] )+ omega[2] * cos ( M_PI * y[j] )+ omega[3] * sin ( M_PI * y[j] ) );

    free ( arg );
    return dc;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _diffusivity_2d_elman ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFUSIVITY_2D_ELMAN evaluates a 2D stochastic diffusivity function.
  Discussion:
    The 2D diffusion equation has the form
      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
    where DC(X,Y) is a function called the diffusivity.
    In the stochastic version of the problem, the diffusivity function
    includes the influence of stochastic parameters:
      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
    In this function, the domain is assumed to be the square [-A,+A]x[-A,+A].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 August 2013
  Author:
    John Burkardt
  Reference:
    Howard Elman, Darran Furnaval,
    Solving the stochastic steady-state diffusion problem using multigrid,
    IMA Journal on Numerical Analysis,
    Volume 27, Number 4, 2007, pages 675-688.
    Roger Ghanem, Pol Spanos,
    Stochastic Finite Elements: A Spectral Approach,
    Revised Edition,
    Dover, 2003,
    ISBN: 0486428184,
    LC: TA347.F5.G56.
  Parameters:
    Input, double A, the "radius" of the square region.  The region
    is assumed to be [-A,+A]x[-A,+A].
    0 < A.
    Input, double CL, the correlation length.
    0 < CL.
    Input, double DC0, the constant term in the expansion of the
    diffusion coefficient.  Take DC0 = 10.
    Input, int M_1D, the first and second dimensions of the
    stochastic parameter array.
    Input, double OMEGA[M_1D*M_1D], the stochastic parameters.
    Input, int N1, N2, the dimensions of the X and Y arrays.
    Input, double X[N1*N2], Y[N1*N2], the points where the diffusion
    coefficient is to be evaluated.
    Output, double DIFFUSIVITY_2D_ELMAN[N1*N2], the value of the diffusion
    coefficient at X.
*/
{
	const _3itdtpit2dt2pit * const s_data = data;
	ityp a = s_data->a0;
	ityp cl = s_data->a1;
	ityp dc0 = s_data->a2;
	const register dim_typ m_1d = s_data->a3;
	ityp * omega = s_data->a4;
	const register dim_typ n1 = s_data->a5;
	const register dim_typ n2 = s_data->a6;
	ityp * x = s_data->a7;
	ityp * y = s_data->a8;
	
    ityp *c_1dx;
    ityp *c_1dy;
    ityp *dc;
    dim_typ i;
    dim_typ i1;
    dim_typ i2;
    dim_typ j;
    dim_typ k;
    ityp *lambda_1d;
    dim_typ m;
    ityp *theta_1d;

    m = m_1d * m_1d;
    /*
    Compute THETA.
    */
    theta_1d = theta_solve ( a, cl, m_1d );
    /*
    Compute LAMBDA_1D.
    */
    lambda_1d = ( ityp * ) malloc ( m_1d * sizeof ( ityp ) );

    for ( i = 0; i < m_1d; ++i )
        lambda_1d[i] = 2.00 * cl / ( 1.00 + cl * cl * theta_1d[i] * theta_1d[i] );
    /*
    Compute C_1DX(1:M1D) and C_1DY(1:M1D) at (X,Y).
    */
    c_1dx = ( ityp * ) malloc ( m_1d * n1 * n2 * sizeof ( ityp ) );
    c_1dy = ( ityp * ) malloc ( m_1d * n1 * n2 * sizeof ( ityp ) );

    for ( k = 0; k < n2; ++k )
        for ( j = 0; j < n1; ++j )
            for ( i = 0; i < m_1d; ++i )
            c_1dx[i+j*m_1d+k*m_1d*n1] = c_1dy[i+j*m_1d+k*m_1d*n1] = 0.00;


    i = 0;

    for ( ; ; )
    {
        if ( m_1d <= i )
            break;

        for ( k = 0; k < n2; ++k )
            for ( j = 0; j < n1; ++j )
            {
                c_1dx[i+j*m_1d+k*m_1d*n1] = cos ( theta_1d[i] * a * x[j+k*n1] )/ sqrt ( a + sin ( 2.00 * theta_1d[i] * a )/ ( 2.00 * theta_1d[i] ) );
                c_1dy[i+j*m_1d+k*m_1d*n1] = cos ( theta_1d[i] * a * y[j+k*n1] )/ sqrt ( a + sin ( 2.00 * theta_1d[i] * a )/ ( 2.00 * theta_1d[i] ) );
            }

        ++ i;

        if ( m_1d <= i )
            break;

        for ( k = 0; k < n2; ++k )
            for ( j = 0; j < n1; ++j )
            {
                c_1dx[i+j*m_1d+k*m_1d*n1] = sin ( theta_1d[i] * a * x[j+k*n1] )/ sqrt ( a - sin ( 2.00 * theta_1d[i] * a )/ ( 2.00 * theta_1d[i] ) );
                c_1dy[i+j*m_1d+k*m_1d*n1] = sin ( theta_1d[i] * a * y[j+k*n1] )/ sqrt ( a - sin ( 2.00 * theta_1d[i] * a )/ ( 2.00 * theta_1d[i] ) );
            }

            ++ i;
    }
    /*
    Evaluate the diffusion coefficient DC at (X,Y).
    */
    dc = ( ityp * ) malloc ( n1 * n2 * sizeof ( ityp ) );

    for ( k = 0; k < n2; ++k )
        for ( j = 0; j < n1; ++j )
        {
            dc[j+k*n1] = dc0;
            for ( i2 = 0; i2 < m_1d; ++i2)
                for ( i1 = 0; i1 < m_1d; ++i1 )
                    dc[j+k*n1] = dc[j+k*n1] + sqrt ( lambda_1d[i1] * lambda_1d[i2] )* c_1dx[i1+j*m_1d+k*m_1d*n1] * c_1dy[i2+j*m_1d+k*m_1d*n1] * omega[i1+i2*m_1d];
        }

    free ( c_1dx );
    free ( c_1dy );
    free ( lambda_1d );
    free ( theta_1d );
    return dc;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _diffusivity_2d_ntw ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIFFUSIVITY_2D_NTW evaluates a 2D stochastic diffusivity function.
  Discussion:
    The 2D diffusion equation has the form
      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
    where DC(X,Y) is a function called the diffusivity.
    In the stochastic version of the problem, the diffusivity function
    includes the influence of stochastic parameters:
      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
    In this function, the domain is the rectangle [0,D]x[0,D] where D = 1.
    Note that in this problem the diffusivity has a one-dimensional
    spatial dependence on X, but not on Y
    The random variables OMEGA are independent, have zero mean and unit
    variance, and are uniformly distributed in [-sqrt(3),+sqrt(3)].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 August 2013
  Author:
    John Burkardt
  Reference:
    Xiang Ma, Nicholas Zabaras,
    An adaptive hierarchical sparse grid collocation algorithm for the solution
    of stochastic differential equations,
    Journal of Computational Physics,
    Volume 228, pages 3084-3113, 2009.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, double CL, the desired physical correlation length for
    the coefficient.
    Input, double DC0, the constant term in the expansion of the
    diffusion coefficient.  Take DC0 = 0.5.
    Input, int M, the number of terms in the expansion.
    Input, double OMEGA[M], the stochastic parameters.
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the points where the diffusion
    coefficient is to be evaluated.
    Output, double DIFFUSIVITY_2D_NTW[N], the value of the diffusion coefficient
    at (X,Y).
*/
{
	const _2itdtpitdt2pit * const s_data = data;
	const register ityp cl = s_data->a0;
	const register ityp dc0 = s_data->a1;
	const register dim_typ m = s_data->a2;
	ityp * omega = s_data->a3;
	const register dim_typ n = s_data->a4;
	ityp * x = s_data->a5;
	ityp * y = s_data->a6;
	
    ityp d;
    ityp *dc;
    ityp *dc_arg;
    dim_typ i, j;
    ityp ihalf_r8;
    ityp l;
    ityp lj;
    ityp lp;
    ityp *phi;
    ityp zeta;
    ityp zeta_arg;

    d = 1.00;
    lp = MAX ( d, 2.00 * cl );
    l = cl / lp;

    dc_arg = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        dc_arg[j] = 1.00 + omega[0] * sqrt ( sqrt ( M_PI ) * l / 2.00 );

    dc = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    phi = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 2; i <= m; i++ )
    {
        ihalf_r8 = ( ityp ) ( i / 2 );
        zeta_arg = - pow ( ihalf_r8 * M_PI * l, 2 ) / 8.00;
        zeta = sqrt ( sqrt ( M_PI ) * l ) * exp ( zeta_arg );

        if ( !( i % 2 ) )
            for ( j = 0; j < n; ++j )
                phi[j] = sin ( ihalf_r8 * M_PI * x[j] / lp );
        else
            for ( j = 0; j < n; ++j )
                phi[j] = cos ( ihalf_r8 * M_PI * x[j] / lp );

        for ( j = 0; j < n; ++j )
            dc_arg[j] += zeta * phi[j] * omega[i-1];
    }

    for ( j = 0; j < n; ++j )
        dc[j] = dc0 + exp ( dc_arg[j] );

    free ( dc_arg );
    free ( phi );
    return dc;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_mesh_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_MESH_2D creates a 2D mesh from X and Y vectors.
  Discussion:
    An R8VEC is a vector of R8's.
    NX = 2
    XVEC = ( 1, 2, 3 )
    NY = 3
    YVEC = ( 4, 5 )
    XMAT = (
      1, 2, 3
      1, 2, 3 )
    YMAT = (
      4, 4, 4
      5, 5, 5 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 July 2013
  Parameters:
    Input, int NX, NY, the number of X and Y values.
    Input, double XVEC[NX], YVEC[NY], the X and Y coordinate
    values.
    Output, double XMAT[NX*NY], YMAT[NX*NY], the coordinate
    values of points on an NX by NY mesh.
*/
{
	const _2dt4pit * const s_data = data;
	const register dim_typ nx = s_data->a0;
	const register dim_typ ny = s_data->a1;
	ityp * xvec = s_data->a2;
	ityp * yvec = s_data->a3;
	ityp * xmat = s_data->a4;
	ityp * ymat = s_data->a5;
	
    dim_typ i;
    dim_typ j;

    for ( j = 0; j < ny; ++j )
        for ( i = 0; i < nx; ++i )
        {
            xmat[i+j*nx] = xvec[i];
            ymat[i+j*nx] = yvec[j];
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _theta_solve ( void * data)
/******************************************************************************/
/*
  Purpose:
    THETA_SOLVE solves a pair of transcendental equations.
  Discussion:
    The vector THETA returned by this function is needed in order to define
    the terms in a Karhunen-Loeve expansion of a diffusion coefficient.
    The two equations are:
      1/CL - THETA * TAN ( A * THETA ) = 0
      THETA - 1/CL * TAN ( A * THETA ) = 0
    A and CL are taken to be positive.  Over each open interval
  ( n - 1/2 M_PI, n + 1/2 M_PI ) / A, for N = 0, 1, ...
    the function TAN ( A * THETA ) monotonically rises from -oo to +00;
    therefore, it can be shown that there is one root of each equation
    in every interval of this form.  Moreover, because of the positivity
    of A and CL, we can restrict our search to the interval
      [ n M_PI, n + 1/2 M_PI ) / A, for N = 0, 1, ...
    This function computes K such roots, starting in the first interval,
    finding those two roots, moving to the next interval, and so on, until
    the requested number of roots have been found.  Odd index roots will
    correspond to the first equation, and even index roots to the second.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 August 2013
  Author:
    John Burkardt
  Reference:
    Howard Elman, Darran Furnival,
    Solving the Stochastic Steady-State Diffusion Problem Using Multigrid,
    University of Maryland Department of Computer Science,
    Technical Report TR-4786.
  Parameters:
    Input, double A, the "radius" of the domain, D = (-A,A)x(-A,A).
    0 < A.
    Input, double CL, the correlation length.
    0 < CL.
    Input, int M, the number of values to compute.
    Output, double THETA_SOLVE[M], the values of Theta.
*/
{
	const _2itdt * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp cl = s_data->a1;
	const register dim_typ m = s_data->a2;
	
    ityp bmatol;
    ityp eps;
    ityp fa;
    ityp fb;
    ityp fc;
    ityp ftol;
    dim_typ k;
    ityp *theta;
    ityp xa;
    ityp xa_init;
    ityp xb;
    ityp xb_init;
    ityp xc;

    theta = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    for ( k = 0; k < m; ++k )
        theta[k] = 0.00;
    /*
    [ XA_INIT, XB_INIT] = [ n * M_PI, n+1/2 M_PI ] / a, n = 0, 1, 2, ...
    */
    xa_init = 0.00;
    xb_init = ( M_PI / 2.00 ) / a;
    eps = r8_epsilon ( );

    k = 0;
    for ( ; ; )
    {
        /*
        Seek root of equation 1 in interval.
        */
        if ( m <= k )
            break;
        ++ k;
        xa = xa_init;
        fa = 1.00 / cl - xa * tan ( a * xa );
        ftol = eps * ( fabs ( fa ) + 1.00 );
        xb = xb_init;
        fb = - fa;
        fc = fa;
        bmatol = 100.00 * eps * ( fabs ( xa ) + fabs ( xb ) );

        while ( bmatol < xb - xa )
        {
            xc = ( xa + xb ) / 2.00;
            fc = 1.00 / cl - xc * tan ( a * xc );

            if ( fabs ( fc ) <= ftol )
                break;
            else if ( 0.0 < fc )
                xa = xc;
            else
                xb = xc;
        }

        theta[k-1] = xc;
        /*
        Seek root of equation 2 in interval.
        */
        if ( m <= k )
            break;

        ++ k;
        /*
        In the first interval, we need to skip the zero root of equation 2.
        */
        if ( k == 2 )
        {
            -- k;
        }
        else
        {
            xa = xa_init;
            fa = xa - tan ( a * xa ) / cl;
            ftol = eps * ( fabs ( fa ) + 1.00 );
            xb = xb_init;
            fb = - fa;

            while ( bmatol < xb - xa )
            {
                xc = ( xa + xb ) / 2.0;
                fc = xc - tan ( a * xc ) / cl;

                if ( fabs ( fc ) <= ftol )
                    break;
                else if ( 0.00 < fc )
                    xa = xc;
                else
                    xb = xc;
            }
            theta[k-1] = xc;
        }
        /*
        Advance the interval.
        */
        xa_init += M_PI / a;
        xb_init += M_PI / a;
    }

    return theta;
}

#endif
