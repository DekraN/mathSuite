#ifndef __DISABLEDEEP_SIMPLERKF45

#include "../dutils.h"


// UNIMPLEMENTED
// NOT PROPERLY WORKING
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _simple_rkf45_run ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLE_RKF45_RUN runs the simple three body ODE system.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 October 2012
  Author:
    John Burkardt
*/
{
    ityp abserr;
    dim_typ flag;
    dim_typ i;
    ityp m0;
    ityp m1;
    ityp m2;
    const register dim_typ neqn = 12;
    ityp relerr;
    dim_typ step;
    const register dim_typ step_num = 630;
    ityp t;
    ityp t_out;
    ityp t_start;
    ityp t_stop;
    ityp *ts;
    ityp *y;
    ityp *yp;
    ityp *ys;

    ts = ( ityp * ) malloc ( ( step_num + 1 ) * sizeof ( ityp ) );
    y = ( ityp * ) malloc ( neqn * sizeof ( ityp ) );
    yp = ( ityp * ) malloc ( neqn * sizeof ( ityp ) );
    ys = ( ityp * ) malloc ( neqn * ( step_num + 1 ) * sizeof ( ityp ) );

    m0 = 5.00;
    m1 = 3.00;
    m2 = 4.00;

    abserr = 1.0E-10;
    relerr = 1.0E-10;

    flag = 1;

    t_start = t = t_out = 0;
    t_stop = 63.00;

    y[0] = y[4] = 1.00;
    y[1] = y[9] = -1.00;
    y[2] = y[3] = y[6] = y[7] = y[10] = y[11] = 0.00;
    y[5] =  3.0;
    y[8] = -2.0;

    simple_f ( t, y, yp );

    for ( i = 0; i < neqn; ++i )
   		ys[i+0*neqn] = y[i];
    ts[0] = t;

    for ( step = 1; step <= step_num; ++step )
    {
        t = ( ( ityp ) ( step_num - step + 1 ) * t_start + ( ityp ) (            step - 1 ) * t_stop ) / ( ityp ) ( step_num            );
        t_out = ( ( ityp ) ( step_num - step ) * t_start + ( ityp ) (            step ) * t_stop ) / ( ityp ) ( step_num        );
        // flag = r8_rkf45 ( simple_f, neqn, y, yp, &t, t_out, &relerr, abserr, flag );
	// UNIMPLEMENTED
        if ( flag == 7 )
            flag = 2;

        for ( i = 0; i < neqn; ++i )
            ys[i+step*neqn] = y[i];
        ts[step] = t_out;
    }

    free ( ts );
    free ( y );
    free ( yp );
    free ( ys );

    return NULL;
}
#endif
