#ifndef __DISABLEDEEP_LORENZODE

#include "../dutils.h"



/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _lorenz_rhs ( void * data)
/******************************************************************************/
/*
  Purpose:
    LORENZ_RHS evaluates the right hand side of the Lorenz ODE.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, double T, the value of the independent variable.
    Input, int M, the spatial dimension.
    Input, double X[M], the values of the dependent variables
    at time T.
    Output, double DXDT[M], the values of the derivatives
    of the dependent variables at time T.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	ityp * x = s_data->a1;
	const register ityp t = s_data->a2;
	
    ityp beta = 8.00 / 3.00;
    ityp *dxdt;
    ityp rho = 28.00;
    ityp sigma = 10.00;
    dxdt = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = x[0] * ( rho - x[2] ) - x[1];
    dxdt[2] = x[0] * x[1] - beta * x[2];
    return dxdt;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   *rk4vec ( const register ityp t0, int m, ityp u0[static m], const register ityp dt, ityp *f ( ityp t, int n, ityp u[] ) )
/******************************************************************************/
/*
  Purpose:
    RK4VEC takes one Runge-Kutta step for a vector ODE.
  Discussion:
    It is assumed that an initial value problem, of the form
      du/dt = f ( t, u )
      u(t0) = u0
    is being solved.
    If the user can supply current values of t, u, a stepsize dt, and a
    function to evaluate the derivative, this function can compute the
    fourth-order Runge Kutta estimate to the solution at time t+dt.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, double T0, the current time.
    Input, int M, the dimension of the space.
    Input, double U0[M], the solution estimate at the current time.
    Input, double DT, the time step.
    Input, double *F ( double T, int M, double U[] ), a function which evaluates
    the derivative, or right hand side of the problem.
    Output, double RK4VEC[M], the fourth-order Runge-Kutta solution estimate
    at time T0+DT.
*/
{
    ityp *f0;
    ityp *f1;
    ityp *f2;
    ityp *f3;
    int i;
    ityp t1;
    ityp t2;
    ityp t3;
    ityp *u;
    ityp *u1;
    ityp *u2;
    ityp *u3;
    /*
    Get four sample values of the derivative.
    */
    f0 = f ( t0, m, u0 );
    t1 = t0 + dt / 2.00;
    u1 = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    for ( i = 0; i < m; ++i)
        u1[i] = u0[i] + dt * f0[i] / 2.00;
    f1 = f ( t1, m, u1 );
    t2 = t0 + dt / 2.0;
    u2 = ( ityp * ) malloc ( m * sizeof (ityp ) );
    for ( i = 0; i < m; ++i )
        u2[i] = u0[i] + dt * f1[i] / 2.00;
    f2 = f ( t2, m, u2 );

    t3 = t0 + dt;
    u3 = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    for ( i = 0; i < m; ++i)
        u3[i] = u0[i] + dt * f2[i];
    f3 = f ( t3, m, u3 );
    /*
    Combine them to estimate the solution.
    */
    u = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    for ( i = 0; i < m; ++i)
        u[i] = u0[i] + dt * ( f0[i] + 2.00 * f1[i] + 2.00 * f2[i] + f3[i] ) / 6.00;
    /*
    Free memory.
    */
    free ( f0 );
    free ( f1 );
    free ( f2 );
    free ( f3 );
    free ( u1 );
    free ( u2 );
    free ( u3 );
    return u;
}

#endif
