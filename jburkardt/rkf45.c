#ifndef __DISABLEDEEP_RKF45

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   r8_fehl ( void f ( ityp t, ityp y[], ityp yp[] ), const register dim_typ neqn,
  ityp y[static neqn], const register ityp t, const register ityp h, ityp yp[static neqn], ityp f1[static neqn], ityp f2[static neqn], ityp f3[static neqn], ityp f4[static neqn], ityp f5[static neqn], ityp s[static neqn] )
/******************************************************************************/
/*
  Purpose:
    R8_FEHL takes one Fehlberg fourth-fifth order step.
  Discussion:
    This version of the routine uses DOUBLE real arithemtic.
    This routine integrates a system of NEQN first order ordinary differential
    equations of the form
      dY(i)/dT = F(T,Y(1:NEQN))
    where the initial values Y and the initial derivatives
    YP are specified at the starting point T.
    The routine advances the solution over the fixed step H and returns
    the fifth order (sixth order accurate locally) solution
    approximation at T+H in array S.
    The formulas have been grouped to control loss of significance.
    The routine should be called with an H not smaller than 13 units of
    roundoff in T so that the various independent arguments can be
    distinguished.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 April 2011
  Author:
    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
    C++ version by John Burkardt.
  Reference:
    Erwin Fehlberg,
    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
    NASA Technical Report R-315, 1969.
    Lawrence Shampine, Herman Watts, S Davenport,
    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
    SIAM Review,
    Volume 18, pages 376-411, 1976.
  Parameters:
    Input, external F, a user-supplied subroutine to evaluate the
    derivatives Y'(T), of the form:
      void f ( double t, double y[], double yp[] )
    Input, int NEQN, the number of equations to be integrated.
    Input, double Y[NEQN], the current value of the dependent variable.
    Input, double T, the current value of the independent variable.
    Input, double H, the step size to take.
    Input, double YP[NEQN], the current value of the derivative of the
    dependent variable.
    Output, double F1[NEQN], F2[NEQN], F3[NEQN], F4[NEQN], F5[NEQN], derivative
    values needed for the computation.
    Output, double S[NEQN], the estimate of the solution at T+H.
*/
{
    ityp ch;
    dim_typ i;

    ch = h / 4.0;

    for ( i = 0; i < neqn; ++i)
        f5[i] = y[i] + ch * yp[i];

    f ( t + ch, f5, f1 );

    ch = 3.00 * h / 32.00;

    for ( i = 0; i < neqn; ++i)
        f5[i] = y[i] + ch * ( yp[i] + 3.00 * f1[i] );

    f ( t + 3.00 * h / 8.00, f5, f2 );
    ch = h / 2197.00;

    for ( i = 0; i < neqn; ++i )
        f5[i] = y[i] + ch * ( 1932.00 * yp[i] + ( 7296.00 * f2[i] - 7200.00 * f1[i] ) );

    f ( t + 12.00 * h / 13.00, f5, f3 );
    ch = h / 4104.00;

    for ( i = 0; i < neqn; ++i )
        f5[i] = y[i] + ch * ( ( 8341.00 * yp[i] - 845.00 * f3[i] ) + ( 29440.00 * f2[i] - 32832.00 * f1[i] ) );

    f ( t + h, f5, f4 );

    ch = h / 20520.00;

    for ( i = 0; i < neqn; ++i )
        f1[i] = y[i] + ch * ( ( -6080.00 * yp[i] + ( 9295.00 * f3[i] - 5643.00 * f4[i] ) ) + ( 41040.00 * f1[i] - 28352.00 * f2[i] ) );

    f ( t + h / 2.0, f1, f5 );
    /*
    Ready to compute the approximate solution at T+H.
    */
    ch = h / 7618050.00;

    for ( i = 0; i < neqn; ++i )
        s[i] = y[i] + ch * ( ( 902880.0 * yp[i] + ( 3855735.0 * f3[i] - 1371249.0 * f4[i] ) ) + ( 3953664.0 * f2[i] + 277020.0 * f5[i] ) );

    return;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  short   r8_rkf45 ( void f ( ityp t, ityp y[], ityp yp[] ), int neqn,
  ityp y[static neqn], ityp yp[static neqn], ityp *t, const register ityp tout, ityp *relerr,const register ityp abserr, int flag )
/******************************************************************************/
/*
  Purpose:
    R8_RKF45 carries out the Runge-Kutta-Fehlberg method.
  Discussion:
    This version of the routine uses DOUBLE real arithmetic.
    This routine is primarily designed to solve non-stiff and mildly stiff
    differential equations when derivative evaluations are inexpensive.
    It should generally not be used when the user is demanding
    high accuracy.
    This routine integrates a system of NEQN first-order ordinary differential
    equations of the form:
      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
    where the Y(1:NEQN) are given at T.
    Typically the subroutine is used to integrate from T to TOUT but it
    can be used as a one-step integrator to advance the solution a
    single step in the direction of TOUT.  On return, the parameters in
    the call list are set for continuing the integration.  The user has
    only to call again (and perhaps define a new value for TOUT).
    Before the first call, the user must
    * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
      and declare F in an EXTERNAL statement;
    * initialize the parameters:
      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
      In particular, T should initially be the starting point for integration,
      Y should be the value of the initial conditions, and FLAG should
      normally be +1.
    Normally, the user only sets the value of FLAG before the first call, and
    thereafter, the program manages the value.  On the first call, FLAG should
    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
    have been reset by the program to the value of 2 (or -2 in single
    step mode), and the user can continue to call the routine with that
    value of FLAG.
 (When the input magnitude of FLAG is 1, this indicates to the program
    that it is necessary to do some initialization work.  An input magnitude
    of 2 lets the program know that that initialization can be skipped,
    and that useful information was computed earlier.)
    The routine returns with all the information needed to continue
    the integration.  If the integration reached TOUT, the user need only
    define a new TOUT and call again.  In the one-step integrator
    mode, returning with FLAG = -2, the user must keep in mind that
    each step taken is in the direction of the current TOUT.  Upon
    reaching TOUT, indicated by the output value of FLAG switching to 2,
    the user must define a new TOUT and reset FLAG to -2 to continue
    in the one-step integrator mode.
    In some cases, an error or difficulty occurs during a call.  In that case,
    the output value of FLAG is used to indicate that there is a problem
    that the user must address.  These values include:
    * 3, integration was not completed because the input value of RELERR, the
      relative error tolerance, was too small.  RELERR has been increased
      appropriately for continuing.  If the user accepts the output value of
      RELERR, then simply reset FLAG to 2 and continue.
    * 4, integration was not completed because more than MAXNFE derivative
      evaluations were needed.  This is approximately (MAXNFE/6) steps.
      The user may continue by simply calling again.  The function counter
      will be reset to 0, and another MAXNFE function evaluations are allowed.
    * 5, integration was not completed because the solution vanished,
      making a pure relative error test impossible.  The user must use
      a non-zero ABSERR to continue.  Using the one-step integration mode
      for one step is a good way to proceed.
    * 6, integration was not completed because the requested accuracy
      could not be achieved, even using the smallest allowable stepsize.
      The user must increase the error tolerances ABSERR or RELERR before
      continuing.  It is also necessary to reset FLAG to 2 (or -2 when
      the one-step integration mode is being used).  The occurrence of
      FLAG = 6 indicates a trouble spot.  The solution is changing
      rapidly, or a singularity may be present.  It often is inadvisable
      to continue.
    * 7, it is likely that this routine is inefficient for solving
      this problem.  Too much output is restricting the natural stepsize
      choice.  The user should use the one-step integration mode with
      the stepsize determined by the code.  If the user insists upon
      continuing the integration, reset FLAG to 2 before calling
      again.  Otherwise, execution will be terminated.
    * 8, invalid input parameters, indicates one of the following:
      NEQN <= 0;
      T = TOUT and |FLAG| /= 1;
      RELERR < 0 or ABSERR < 0;
      FLAG == 0  or FLAG < -2 or 8 < FLAG.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2004
  Author:
    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
    C++ version by John Burkardt.
  Reference:
    Erwin Fehlberg,
    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
    NASA Technical Report R-315, 1969.
    Lawrence Shampine, Herman Watts, S Davenport,
    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
    SIAM Review,
    Volume 18, pages 376-411, 1976.
  Parameters:
    Input, external F, a user-supplied subroutine to evaluate the
    derivatives Y'(T), of the form:
      void f ( double t, double y[], double yp[] )
    Input, int NEQN, the number of equations to be integrated.
    Input/output, double Y[NEQN], the current solution vector at T.
    Input/output, double YP[NEQN], the derivative of the current solution
    vector at T.  The user should not set or alter this information!
    Input/output, double *T, the current value of the independent variable.
    Input, double TOUT, the output point at which solution is desired.
    TOUT = T is allowed on the first call only, in which case the routine
    returns with FLAG = 2 if continuation is possible.
    Input, double *RELERR, ABSERR, the relative and absolute error tolerances
    for the local error test.  At each step the code requires:
      abs ( local error ) <= RELERR * abs ( Y ) + ABSERR
    for each component of the local error and the solution vector Y.
    RELERR cannot be "too small".  If the routine believes RELERR has been
    set too small, it will reset RELERR to an acceptable value and return
    immediately for user action.
    Input, int FLAG, indicator for status of integration. On the first call,
    set FLAG to +1 for normal use, or to -1 for single step mode.  On
    subsequent continuation steps, FLAG should be +2, or -2 for single
    step mode.
    Output, int RKF45_D, indicator for status of integration.  A value of 2
    or -2 indicates normal progress, while any other value indicates a
    problem that should be addressed.
*/
{
    # define MAXNFE 3000

    static ityp abserr_save = -1.00;
    ityp ae;
    ityp dt;
    ityp ee;
    ityp eeoet;
    ityp eps;
    ityp esttol;
    ityp et;
    ityp *f1;
    ityp *f2;
    ityp *f3;
    ityp *f4;
    ityp *f5;
    int flag_return;
    static int flag_save = -1000;
    static ityp h = -1.00;
    int hfaild;
    ityp hmin;
    int i;
    static int init = -1000;
    int k;
    static int kflag = -1000;
    static int kop = -1;
    int mflag;
    static int nfe = -1;
    int output;
    double relerr_min;
    static ityp relerr_save = -1.00;
    static ityp remin = 1.0E-12;
    ityp s;
    ityp scale;
    ityp tol;
    ityp toln;
    ityp ypk;

    flag_return = flag;
    /*
    Check the input parameters.
    */
    eps = r8_epsilon ( );

    if ( neqn < 1 || (*relerr) < 0.00 || abserr < 0.00 || flag_return == 0 || 8 < flag_return  || flag_return < -2  )
        return 8;

    mflag = abs ( flag_return );
    /*
    Is this a continuation call?
    */
    if ( mflag != 1 )
    {
        if ( *t == tout && kflag != 3 )
            return 8;
        /*
        FLAG = -2 or +2:
        */
        if ( mflag == 2 )
        {
            if ( kflag == 3 )
            {
                flag_return = flag_save;
                mflag = abs ( flag_return );
            }
            else if ( init == 0 )
                flag_return = flag_save;
            else if ( kflag == 4 )
                nfe = 0;
            else if ( kflag == 5 && abserr == 0.0 )
                return SHRT_MAX;
            else if ( kflag == 6 && (*relerr) <= relerr_save && abserr <= abserr_save )
                return SHRT_MAX;
        }
        /*
        FLAG = 3, 4, 5, 6, 7 or 8.
        */
        else
        {
            if ( flag_return == 3 )
            {
                flag_return = flag_save;
                if ( kflag == 3 )
                mflag = abs ( flag_return );
            }
            else if ( flag_return == 4 )
            {
                nfe = 0;
                flag_return = flag_save;
                if ( kflag == 3 )
                    mflag = abs ( flag_return );
            }
            else if ( flag_return == 5 && 0.0 < abserr )
            {
                flag_return = flag_save;
                if ( kflag == 3 )
                    mflag = abs ( flag_return );
            }
            /*
            Integration cannot be continued because the user did not respond to
            the instructions pertaining to FLAG = 5, 6, 7 or 8.
            */
            else
                return SHRT_MAX;
        }
    }
    /*
    Save the input value of FLAG.
    Set the continuation flag KFLAG for subsequent input checking.
    */
    flag_save = flag_return;
    kflag = 0;
    /*
    Save RELERR and ABSERR for checking input on subsequent calls.
    */
    relerr_save = (*relerr);
    abserr_save = abserr;
    /*
    Restrict the relative error tolerance to be at least

    2*EPS+REMIN

    to avoid limiting precision difficulties arising from impossible
    accuracy requests.
    */
    relerr_min = 2.00 * r8_epsilon ( ) + remin;
    /*
    Is the relative error tolerance too small?
    */
    if ( (*relerr) < relerr_min )
    {
  (*relerr) = relerr_min;
        kflag = 3;
        return 3;
    }

    dt = tout - *t;
    /*
    Initialization:

    Set the initialization completion indicator, INIT;
    set the indicator for too many output points, KOP;
    evaluate the initial derivatives
    set the counter for function evaluations, NFE;
    estimate the starting stepsize.
    */
    f1 = ( ityp * ) malloc ( neqn * sizeof ( ityp ) );
    f2 = ( ityp * ) malloc ( neqn * sizeof ( ityp ) );
    f3 = ( ityp * ) malloc ( neqn * sizeof ( ityp ) );
    f4 = ( ityp * ) malloc ( neqn * sizeof ( ityp ) );
    f5 = ( ityp * ) malloc ( neqn * sizeof ( ityp ) );

    if ( mflag == 1 )
    {
        init = kop = 0;
        f ( *t, y, yp );
        nfe = 1;

        if ( *t == tout )
            return 2;
    }

    if ( init == 0 )
    {
        init = 1;
        h = abs ( dt );
        toln = 0.00;

        for ( k = 0; k < neqn; ++k )
        {
            tol = (*relerr) * abs ( y[k] ) + abserr;
            if ( 0.00 < tol )
            {
                toln = tol;
                ypk = abs ( yp[k] );
                if ( tol < ypk * pow ( h, 5 ) )
                    h = pow ( ( tol / ypk ), 0.2 );
            }
        }

        if ( toln <= 0.00 )
            h = 0.00;

        h = MAX ( h, 26.0 * eps * MAX ( abs ( *t ), abs ( dt ) ) );

        flag_save = 2 - ((flag_return<0)<<2);
    }
    /*
    Set stepsize for integration in the direction from T to TOUT.
    */
    h = r8_sign ( dt ) * abs ( h );
    /*
    Test to see if too may output points are being requested.
    */
    if ( 2.00 * abs ( dt ) <= abs ( h ) )
        ++ kop;
    /*
    Unnecessary frequency of output.
    */
    if ( kop == 100 )
    {
        kop = 0;
        free ( f1 );
        free ( f2 );
        free ( f3 );
        free ( f4 );
        free ( f5 );
        return 7;
    }
    /*
    If we are too close to the output point, then simply extrapolate and return.
    */
    if ( abs ( dt ) <= 26.00 * eps * abs ( *t ) )
    {
        *t = tout;
        for ( i = 0; i < neqn; ++i )
            y[i] += dt * yp[i];
        f ( *t, y, yp );
        nfe = nfe + 1;

        free ( f1 );
        free ( f2 );
        free ( f3 );
        free ( f4 );
        free ( f5 );
        return 2;
    }
    /*
    Initialize the output point indicator.
    */
    output = 0;
    /*
    To avoid premature underflow in the error tolerance function,
    scale the error tolerances.
    */
    scale = 2.00 / (*relerr);
    ae = scale * abserr;
    /*
    Step by step integration.
    */
    for ( ; ; )
    {
        hfaild = 0;
        /*
        Set the smallest allowable stepsize.
        */
        hmin = 26.00 * eps * abs ( *t );
        /*
        Adjust the stepsize if necessary to hit the output point.

        Look ahead two steps to avoid drastic changes in the stepsize and
        thus lessen the impact of output points on the code.
        */
        dt = tout - *t;

        if ( 2.00 * abs ( h ) <= abs ( dt ) );
        else
        /*
        Will the next successful step complete the integration to the output point?
        */
        {
            if ( abs ( dt ) <= abs ( h ) )
            {
                output = 1;
                h = dt;
            }
            else
                h = 0.50 * dt;
        }
        /*
        Here begins the core integrator for taking a single step.

        The tolerances have been scaled to avoid premature underflow in
        computing the error tolerance function ET.
        To avoid problems with zero crossings, relative error is measured
        using the average of the magnitudes of the solution at the
        beginning and end of a step.
        The error estimate formula has been grouped to control loss of
        significance.

        To distinguish the various arguments, H is not permitted
        to become smaller than 26 units of roundoff in T.
        Practical limits on the change in the stepsize are enforced to
        smooth the stepsize selection process and to avoid excessive
        chattering on problems having discontinuities.
        To prevent unnecessary failures, the code uses 9/10 the stepsize
        it estimates will succeed.

        After a step failure, the stepsize is not allowed to increase for
        the next attempted step.  This makes the code more efficient on
        problems having discontinuities and more effective in general
        since local extrapolation is being used and extra caution seems
        warranted.

        Test the number of derivative function evaluations.
        If okay, try to advance the integration from T to T+H.
        */
        for ( ; ; )
        {
            /*
            Have we done too much work?
            */
            if ( MAXNFE < nfe )
            {
                kflag = 4;
                free ( f1 );
                free ( f2 );
                free ( f3 );
                free ( f4 );
                free ( f5 );
                return 4;
            }
            /*
            Advance an approximate solution over one step of length H.
            */
            r8_fehl ( f, neqn, y, *t, h, yp, f1, f2, f3, f4, f5, f1 );
            nfe = nfe + 5;
            /*
            Compute and test allowable tolerances versus local error estimates
            and remove scaling of tolerances.  The relative error is
            measured with respect to the average of the magnitudes of the
            solution at the beginning and end of the step.
            */
            eeoet = 0.00;

            for ( k = 0; k < neqn; ++k )
            {
                et = abs ( y[k] ) + abs ( f1[k] ) + ae;

                if ( et <= 0.0 )
                {
                    free ( f1 );
                    free ( f2 );
                    free ( f3 );
                    free ( f4 );
                    free ( f5 );
                    return 5;
                }
                ee = abs( ( -2090.0 * yp[k]+ ( 21970.0 * f3[k] - 15048.0 * f4[k] ))+ ( 22528.0 * f2[k] - 27360.0 * f5[k] ));
                eeoet = MAX ( eeoet, ee / et );
            }

            esttol = abs ( h ) * eeoet * scale / 752400.0;

            if ( esttol <= 1.00 )
                break;
            /*
            Unsuccessful step.  Reduce the stepsize, try again.
            The decrease is limited to a factor of 1/10.
            */
            hfaild = 1;
            output = 0;

            s = esttol < 59049.00 ? 0.90 / pow ( esttol, 0.20 ) : 0.10;
            h = s * h;

            if ( abs ( h ) < hmin )
            {
                kflag = 6;
                free ( f1 );
                free ( f2 );
                free ( f3 );
                free ( f4 );
                free ( f5 );
                return 6;
            }
        }
        /*
        We exited the loop because we took a successful step.
        Store the solution for T+H, and evaluate the derivative there.
        */
        *t += h;
        for ( i = 0; i < neqn; ++i )
            y[i] = f1[i];
        f ( *t, y, yp );
        ++ nfe;
        /*
        Choose the next stepsize.  The increase is limited to a factor of 5.
        If the step failed, the next stepsize is not allowed to increase.
        */

        s = 0.0001889568 < esttol ? 0.90 / pow ( esttol, 0.20 ) : 5.00;

        if ( hfaild )
            s = MIN ( s, 1.00 );
        h = r8_sign ( h ) * MAX ( s * abs ( h ), hmin );
        /*
        End of core integrator

        Should we take another step?
        */
        if ( output )
        {
            *t = tout;
            free ( f1 );
            free ( f2 );
            free ( f3 );
            free ( f4 );
            free ( f5 );
            return 2;
        }

        if ( flag_return <= 0 )
        {
            free ( f1 );
            free ( f2 );
            free ( f3 );
            free ( f4 );
            free ( f5 );
            return -2;
        }

    }
    # undef MAXNFE
}
#endif
