#ifndef __DISABLEDEEP_MULTIGRIDPOISSON1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ctof ( void * data)
/******************************************************************************/
/*
  Purpose:
    CTOF transfers data from a coarse to a finer grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 December 2011
  Author:
    John Burkardt
  Reference:
    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.
  Parameters:
    Input, int NC, the number of coarse nodes.
    Input, double UC[NC], the coarse correction data.
    Input, int NF, the number of fine nodes.
    Input/output, double UF[NF], on input, the fine grid data.
    On output, the data has been updated with prolonged coarse
    correction data.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ nc = s_data->a0;
	const register dim_typ nf = s_data->a1;
	ityp * uc = s_data->a2;
	ityp * uf = s_data->a3;
	
    dim_typ ic;
    dim_typ iff;

    for ( ic = 0; ic < nc; ++ic )
    {
        iff = ic<<1;
        uf[iff] += uc[ic];
    }

    for ( ic = 0; ic < nc - 1; ++ic )
    {
        iff = (ic<<1) + 1;
        uf[iff] += 0.5 * ( uc[ic] + uc[ic+1] );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ftoc ( void * data)
/******************************************************************************/
/*
  Purpose:
    FTOC transfers data from a fine grid to a coarser grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 December 2011
  Author:
    John Burkardt
  Reference:
    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.
  Parameters:
    Input, int NF, the number of fine nodes.
    Input, double UF[NF], the fine data.
    Input, double RF[NF], the right hand side for the fine grid.
    Input, int NC, the number of coarse nodes.
    Output, double UC[NC], the coarse grid data, set to zero.
    Output, double RC[NC], the right hand side for the coarse grid.
*/
{
	const _2dt4pit * const s_data = data;
	
	const register dim_typ nf = s_data->a0;
	const register dim_typ nc = s_data->a1;
	ityp * uf = s_data->a2;
	ityp * rf = s_data->a3;
	ityp * uc = s_data->a4;
	ityp * rc = s_data->a5;
	
	
    dim_typ ic;
    dim_typ iff;

    for ( ic = 0; ic < nc; ++ic )
        uc[ic] = 0.00;

    rc[0] = 0.00;
    for ( ic = 1; ic < nc - 1; ++ic )
    {
        iff = (ic<<1);
        rc[ic] = 4.00 * ( rf[iff] + uf[iff-1] - 2.00 * uf[iff] + uf[iff+1] );
    }
    rc[nc-1] = 0.0;

    return NULL;
} 

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gauss_seidel ( void * data)
/******************************************************************************/
/*
  Purpose:
    GAUSS_SEIDEL carries out one step of a Gauss-Seidel iteration.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 December 2011
  Author:
    John Burkardt
  Reference:
    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.
  Parameters:
    Input, int N, the number of unknowns.
    Input, double R[N], the right hand side.
    Input/output, double U[N], the estimated solution.
    Output, double *DIF_L1, the L1 norm of the difference between the
    input and output solution estimates.
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * r = s_data->a1;
	ityp * u = s_data->a2;
	ityp * dif_l1 = s_data->a3;
	
    dim_typ i;
    ityp u_old;

    *dif_l1 = 0.00;

    for ( i = 1; i < n - 1; ++i )
    {
        u_old = u[i];
        u[i] = 0.50 * ( u[i-1] + u[i+1] + r[i] );*dif_l1 = *dif_l1 + fabs ( u[i] - u_old );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _monogrid_poisson_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONOGRID_POISSON_1D solves a 1D PDE, using the Gauss-Seidel method.
  Discussion:
    This routine solves a 1D boundary value problem of the form
      - U''(X) = F(X) for A < X < B,
    with boundary conditions U(A) = UA, U(B) = UB.
    The Gauss-Seidel method is used.
    This routine is provided primarily for comparison with the
    multigrid solver.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 July 2014
  Author:
    John Burkardt
  Reference:
    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.
  Parameters:
    Input, int N, the number of intervals.
    Input, double A, B, the left and right endpoints.
    Input, double UA, UB, the left and right boundary values.
    Input, double FORCE ( double x ), the name of the function
    which evaluates the right hand side.
    Input, double EXACT ( double x ), the name of the function
    which evaluates the exact solution.
    Output, int *IT_NUM, the number of iterations.
    Output, double U[N+1], the computed solution.
*/
{
	const dt4it2fitpdtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp a = s_data->a1;
	ityp b = s_data->a2;
	ityp ua = s_data->a3;
	ityp ub = s_data->a4;
	ityp (* force)(ityp) = s_data->a5;
	ityp (* exact)(ityp) = s_data->a6;
	dim_typ * it_num = s_data->a7;
	ityp * u = s_data->a8;
	
    ityp d1;
    ityp h;
    int i;
    ityp *r;
    ityp tol = 0.0001;
    ityp *x;
    /*
    Initialization.
    */
    /*
    Set the nodes.
    */
    x = r8vec_linspace_new ( n + 1, a, b );
    /*
    Set the right hand side.
    */
    r = ( ityp * ) malloc ( ( n + 1 ) * sizeof ( ityp ) );

    r[0] = ua;
    h = ( b - a ) / ( ityp ) ( n );
    for ( i = 1; i < n; ++i )
        r[i] = h * h * force ( x[i] );
    r[n] = ub;

    for ( i = 0; i <= n; ++i )
        u[i] = 0.00;
    *it_num = 0;
    /*
    Gauss-Seidel iteration.
    */
    for ( ; ; )
    {
        ++ *it_num;

        gauss_seidel ( n + 1, r, u, &d1 );

        if ( d1 <= tol )
            break;
    }
    /*
    Free memory.
    */
    free ( r );
    free ( x );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _multigrid_poisson_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    MULTIGRID_POISSON_1D solves a 1D PDE using the multigrid method.
  Discussion:
    This routine solves a 1D boundary value problem of the form
      - U''(X) = F(X) for A < X < B,
    with boundary conditions U(A) = UA, U(B) = UB.
    The multigrid method is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 July 2014
  Author:
    Original FORTRAN77 version by William Hager.
    C version by John Burkardt.
  Reference:
    William Hager,
    Applied Numerical Linear Algebra,
    Prentice-Hall, 1988,
    ISBN13: 978-0130412942,
    LC: QA184.H33.
  Parameters:
    Input, int N, the number of intervals.
    N must be a power of 2.
    Input, double A, B, the left and right endpoints.
    Input, double UA, UB, the left and right boundary values.
    Input, double FORCE ( double x ), the name of the function
    which evaluates the right hand side.
    Input, double EXACT ( double x ), the name of the function
    which evaluates the exact solution.
    Output, int *IT_NUM, the number of iterations.
    Output, double U[N+1], the computed solution.
*/
{
	const dt4it2fitpdtpit * const s_data = data;
	register dim_typ n = s_data->a0;
	ityp a = s_data->a1;
	ityp b = s_data->a2;
	ityp ua = s_data->a3;
	ityp ub = s_data->a4;
	ityp (* force)(ityp) = s_data->a5;
	ityp (* exact)(ityp) = s_data->a6;
	dim_typ * it_num = s_data->a7;
	ityp * u = s_data->a8;
	
    ityp d0;
    ityp d1;
    ityp h;
    dim_typ i;
    dim_typ it;
    dim_typ j;
    dim_typ k;
    dim_typ l;
    dim_typ ll;
    dim_typ m;
    dim_typ nl;
    ityp *r;
    ityp s;
    ityp tol;
    ityp utol;
    ityp *uu;
    ityp *x;
    /*
    Determine if we have enough storage.
    */
    k = i4_log_2 ( n );

    if ( n != powi ( 2, k ) )
        return NULL;

    nl = n + n + k - 2;
    /*
    Initialization.
    */
    it = 4;
    *it_num = 0;
    tol = 0.0001;
    utol = 0.70;
    m = n;
    /*
    Set the nodes.
    */
    x = r8vec_linspace_new ( n + 1, a, b );
    /*
    Set the right hand side.
    */
    r = ( ityp * ) malloc ( nl * sizeof ( ityp ) );
    r[0] = ua;
    h = ( b - a ) / ( ityp ) ( n );
    for ( i = 1; i < n; ++i)
        r[i] = h * h * force ( x[i] );
    r[n] = ub;

    uu = ( ityp * ) malloc ( nl * sizeof ( ityp ) );

    for ( i = 0; i < nl; ++i )
        uu[i] = 0.00;
    /*
    L points to first entry of solution
    LL points to penultimate entry.
    */
    l = j = 0;
    ll = n - 1;
    /*
    Gauss-Seidel iteration
    */
    d1 = 0.00;

    for ( ; ; )
    {
        d0 = d1;
        ++ j;
        gauss_seidel ( n + 1, r + l, uu + l, &d1 );
        ++ *it_num;
        /*
        Do at least 4 iterations at each level.
        */
        if ( j < it )
            continue;
        /*
        Enough iterations, satisfactory decrease, on finest grid, exit.
        */
        else if ( d1 < tol && n == m )
            break;
        /*
        Enough iterations, satisfactory convergence, go finer.
        */
        else if ( d1 < tol )
        {
            ctof ( n + 1, uu + l, n + n + 1, uu + (l-1-n-n) );
            n += n;
            ll = l - 2;
            l -= 1 - n;
            j = 0;
        }
        /*
        Enough iterations, slow convergence, 2 < N, go coarser.
        */
        else if ( utol * d0 <= d1 && 2 < n )
        {
            ftoc ( n + 1, uu + l, r + l, (n/2)+1, uu+(l+n+1), r+(l+n+1) );

            n /= 2;
            l = ll + 2;
            ll += n + 1;
            j = 0;
        }
    }

    for ( i = 0; i < n + 1; ++i )
        u[i] = uu[i];

    /*
    Free memory.
    */
    free ( r );
    free ( uu );
    free ( x );

    return NULL;
}

#endif
