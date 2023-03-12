#ifndef __DISABLEDEEP_RK4

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   rk4 ( const register ityp t0, const register ityp u0, const register ityp dt, ityp f ( ityp t, ityp u ) )
/******************************************************************************/
/*
  Purpose:
    RK4 takes one Runge-Kutta step for a scalar ODE.
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
    Input, double U0, the solution estimate at the current time.
    Input, double DT, the time step.
    Input, double F ( double T, double U ), a function which evaluates
    the derivative, or right hand side of the problem.
    Output, double RK4, the fourth-order Runge-Kutta solution estimate
    at time T0+DT.
*/
{
    ityp f0;
    ityp f1;
    ityp f2;
    ityp f3;
    ityp t1;
    ityp t2;
    ityp t3;

    ityp u1;
    ityp u2;
    ityp u3;
    /*
    Get four sample values of the derivative.
    */
    f0 = f ( t0, u0 );

    t1 = t0 + dt / 2.0;
    u1 = u0 + dt * f0 / 2.0;
    f1 = f ( t1, u1 );

    t2 = t0 + dt / 2.0;
    u2 = u0 + dt * f1 / 2.0;
    f2 = f ( t2, u2 );

    t3 = t0 + dt;
    u3 = u0 + dt * f2;
    f3 = f ( t3, u3 );
    /*
    Combine them to estimate the solution.
    */
    return u0 + dt * ( f0 + 2.00 * f1 + 2.00 * f2 + f3 ) / 6.00;
}
#endif
