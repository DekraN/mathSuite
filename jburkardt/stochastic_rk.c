#ifndef __DISABLEDEEP_STOCHASTICRK

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rk1_ti_step ( void * data)
/******************************************************************************/
/*
  Purpose:
    RK1_TI_STEP takes one step of a stochastic Runge Kutta scheme.
  Discussion:
    The Runge-Kutta scheme is first-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.
    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 July 2010
  Author:
    John Burkardt
  Reference:
    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
  Parameters:
    Input, double X, the value at the current time.
    Input, double T, the current time.
    Input, double H, the time step.
    Input, double Q, the spectral density of the input white noise.
    Input, double FI ( double X ), the name of the deterministic
    right hand side function.
    Input, double GI ( double X ), the name of the stochastic
    right hand side function.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double RK1_TI_STEP, the value at time T+H.
*/
{
	static ityp result = MAX_VAL;
	
	const _4it2fitpi * const s_data = data;
	ityp x = s_data->a0;
	ityp y = s_data->a1;
	ityp h = s_data->a2;
	ityp q = s_data->a3;
	ityp (* fi)(ityp) = s_data->a4;
	ityp (* gi)(ityp) = s_data->a5;
	int * seed = s_data->a6;
	
	result = x + h * fi ( x ) + h * gi ( x ) * r8_normal_01 ( seed ) * sqrt ( q/h );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rk2_ti_step ( void * data)
/******************************************************************************/
/*
  Purpose:
    RK2_TI_STEP takes one step of a stochastic Runge Kutta scheme.
  Discussion:
    The Runge-Kutta scheme is second-order, and suitable for time-invariant
    systems.
    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 July 2010
  Author:
    John Burkardt
  Reference:
    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
  Parameters:
    Input, double X, the value at the current time.
    Input, double T, the current time.
    Input, double H, the time step.
    Input, double Q, the spectral density of the input white noise.
    Input, double FI ( double X ), the name of the deterministic
    right hand side function.
    Input, double GI ( double X ), the name of the stochastic
    right hand side function.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double RK2_TI_STEP, the value at time T+H.
*/
{
	static ityp result = MAX_VAL;
	
	const _4it2fitpi * const s_data = data;
	ityp x = s_data->a0;
	ityp t = s_data->a1;
	ityp h = s_data->a2;
	ityp q = s_data->a3;
	ityp (* fi)(ityp) = s_data->a4;
	ityp (* gi)(ityp) = s_data->a5;
	int * seed = s_data->a6;
	
	ityp  a21;
	ityp a31;
	ityp a32;
	ityp k1;
	ityp k2;
	ityp q1;
	ityp q2;
	ityp t1;
	ityp t2;
	ityp w1;
	ityp w2;
	ityp x1;
	ityp x2;
	
	a21 =   1.00;
	a31 = a32 = 0.50;
	q1 = q2 = 2.00;
	
	t1 = t;
	x1 = x;
	w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
	k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1;
	
	t2 = t1 + a21 * h;
	x2 = x1 + a21 * k1;
	w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
	k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2;
	
	result = x1 + a31 * k1 + a32 * k2;
  	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rk3_ti_step ( void * data)
/******************************************************************************/
/*
  Purpose:
    RK3_TI_STEP takes one step of a stochastic Runge Kutta scheme.
  Discussion:
    The Runge-Kutta scheme is third-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.
    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 July 2010
  Author:
    John Burkardt
  Reference:
    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
  Parameters:
    Input, double X, the value at the current time.
    Input, double T, the current time.
    Input, double H, the time step.
    Input, double Q, the spectral density of the input white noise.
    Input, double FI ( double X ), the name of the deterministic
    right hand side function.
    Input, double GI ( double X ), the name of the stochastic
    right hand side function.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double RK3_TI_STEP, the value at time T+H.
*/
{
	static ityp result = MAX_VAL;
	
	const _4it2fitpi * const s_data = data;
	ityp x = s_data->a0;
	ityp t = s_data->a1;
	ityp h = s_data->a2;
	ityp q = s_data->a3;
	ityp (* fi)(ityp) = s_data->a4;
	ityp (* gi)(ityp) = s_data->a5;
	int * seed = s_data->a6;
	
    ityp a21;
    ityp a31;
    ityp a32;
    ityp a41;
    ityp a42;
    ityp a43;
    ityp k1;
    ityp k2;
    ityp k3;
    ityp q1;
    ityp q2;
    ityp q3;
    ityp t1;
    ityp t2;
    ityp t3;
    ityp w1;
    ityp w2;
    ityp w3;
    ityp x1;
    ityp x2;
    ityp x3;

    a21 =   1.52880952525675;
    a31 =   0.00;
    a32 =   0.51578733443615;
    a41 =   0.53289582961739;
    a42 =   0.25574324768195;
    a43 =   0.21136092270067;

    q1 = 1.87653936176981;
    q2 = 3.91017166264989;
    q3 = 4.73124353935667;

    t1 = t;
    x1 = x;
    w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
    k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1;

    t2 = t1 + a21 * h;
    x2 = x1 + a21 * k1;
    w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
    k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2;

    t3 = t1 + a31 * h  + a32 * h;
    x3 = x1 + a31 * k1 + a32 * k2;
    w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h );
    k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3;
	
	result = x1 + a41 * k1 + a42 * k2 + a43 * k3;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rk4_ti_step ( void * data)
/******************************************************************************/
/*
  Purpose:
    RK4_TI_STEP takes one step of a stochastic Runge Kutta scheme.
  Discussion:
    The Runge-Kutta scheme is fourth-order, and suitable for time-invariant
    systems in which F and G do not depend explicitly on time.
    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 July 2010
  Author:
    John Burkardt
  Reference:
    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
  Parameters:
    Input, double X, the value at the current time.
    Input, double T, the current time.
    Input, double H, the time step.
    Input, double Q, the spectral density of the input white noise.
    Input, double FI ( double X ), the name of the deterministic
    right hand side function.
    Input, double GI ( double X ), the name of the stochastic
    right hand side function.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double RK4_TI_STEP, the value at time T+H.
*/
{
	static ityp result = MAX_VAL;
	
	const _4it2fitpi * const s_data = data;
	ityp x = s_data->a0;
	ityp t = s_data->a1;
	ityp h = s_data->a2;
	ityp q = s_data->a3;
	ityp (* fi)(ityp) = s_data->a4;
	ityp (* gi)(ityp) = s_data->a5;
	int * seed = s_data->a6;
	
    ityp a21;
    ityp a31;
    ityp a32;
    ityp a41;
    ityp a42;
    ityp a43;
    ityp a51;
    ityp a52;
    ityp a53;
    ityp a54;
    ityp k1;
    ityp k2;
    ityp k3;
    ityp k4;
    ityp q1;
    ityp q2;
    ityp q3;
    ityp q4;
    ityp t1;
    ityp t2;
    ityp t3;
    ityp t4;
    ityp w1;
    ityp w2;
    ityp w3;
    ityp w4;
    ityp x1;
    ityp x2;
    ityp x3;
    ityp x4;

    a21 =   2.71644396264860;
    a31 = - 6.95653259006152;
    a32 =   0.78313689457981;
    a41 =   0.0;
    a42 =   0.48257353309214;
    a43 =   0.26171080165848;
    a51 =   0.47012396888046;
    a52 =   0.36597075368373;
    a53 =   0.08906615686702;
    a54 =   0.07483912056879;

    q1 =   2.12709852335625;
    q2 =   2.73245878238737;
    q3 =  11.22760917474960;
    q4 =  13.36199560336697;

    t1 = t;
    x1 = x;
    w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
    k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1;

    t2 = t1 + a21 * h;
    x2 = x1 + a21 * k1;
    w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
    k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2;

    t3 = t1 + a31 * h  + a32 * h;
    x3 = x1 + a31 * k1 + a32 * k2;
    w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h );
    k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3;

    t4 = t1 + a41 * h  + a42 * h + a43 * h;
    x4 = x1 + a41 * k1 + a42 * k2;
    w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h );
    k4 = h * fi ( x4 ) + h * gi ( x4 ) * w4;

	result = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  inline ityp   rk1_tv_step ( ityp x, ityp t, ityp h, ityp q,ityp fv ( ityp t, ityp x ), ityp gv ( ityp t, ityp x ),int *seed )
/******************************************************************************/
/*
  Purpose:
    RK1_TV_STEP takes one step of a stochastic Runge Kutta scheme.
  Discussion:
    The Runge-Kutta scheme is first-order, and suitable for time-varying
    systems.
    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 July 2010
  Author:
    John Burkardt
  Reference:
    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
  Parameters:
    Input, double X, the value at the current time.
    Input, double T, the current time.
    Input, double H, the time step.
    Input, double Q, the spectral density of the input white noise.
    Input, double FV ( double T, double X ), the name of the deterministic
    right hand side function.
    Input, double GV ( double T, double X ), the name of the stochastic
    right hand side function.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double RK1_TV_STEP the value at time T+H.
*/
{
    return x + h * fv ( t, x ) + h * gv ( t, x ) * r8_normal_01 ( seed ) * sqrt ( q / h );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   rk2_tv_step ( ityp x, ityp t, ityp h, ityp q,ityp fv ( ityp t, ityp x ), ityp gv ( ityp t, ityp x ),int *seed )
/******************************************************************************/
/*
  Purpose:
    RK2_TV_STEP takes one step of a stochastic Runge Kutta scheme.
  Discussion:
    The Runge-Kutta scheme is second-order, and suitable for time-varying
    systems.
    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 July 2010
  Author:
    John Burkardt
  Reference:
    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
  Parameters:
    Input, double X, the value at the current time.
    Input, double T, the current time.
    Input, double H, the time step.
    Input, double Q, the spectral density of the input white noise.
    Input, double FV ( double T, double X ), the name of the deterministic
    right hand side function.
    Input, double GV ( double T, double X ), the name of the stochastic
    right hand side function.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double RK2_TV_STEP, the value at time T+H.
*/
{
    ityp a21;
    ityp a31;
    ityp a32;
    ityp k1;
    ityp k2;
    ityp q1;
    ityp q2;
    ityp t1;
    ityp t2;
    ityp w1;
    ityp w2;
    ityp x1;
    ityp x2;

    a21 =   1.00;
    a31 = a32 = 0.50;

    q1 = q2 = 2.00;

    t1 = t;
    x1 = x;
    w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
    k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1;

    t2 = t1 + a21 * h;
    x2 = x1 + a21 * k1;
    w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
    k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2;

    return x1 + a31 * k1 + a32 * k2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   rk4_tv_step ( ityp x, ityp t, ityp h, ityp q,ityp fv ( ityp t, ityp x ), ityp gv ( ityp t, ityp x ),int *seed )
/******************************************************************************/
/*
  Purpose:
    RK4_TV_STEP takes one step of a stochastic Runge Kutta scheme.
  Discussion:
    The Runge-Kutta scheme is fourth-order, and suitable for time-varying
    systems.
    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 July 2010
  Author:
    John Burkardt
  Reference:
    Jeremy Kasdin,
    Runge-Kutta algorithm for the numerical integration of
    stochastic differential equations,
    Journal of Guidance, Control, and Dynamics,
    Volume 18, Number 1, January-February 1995, pages 114-120.
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
  Parameters:
    Input, double X, the value at the current time.
    Input, double T, the current time.
    Input, double H, the time step.
    Input, double Q, the spectral density of the input white noise.
    Input, double FV ( double T, double X ), the name of the deterministic
    right hand side function.
    Input, double GV ( double T, double X ), the name of the stochastic
    right hand side function.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double RK4_TV_STEP, the value at time T+H.
*/
{
    ityp a21;
    ityp a31;
    ityp a32;
    ityp a41;
    ityp a42;
    ityp a43;
    ityp a51;
    ityp a52;
    ityp a53;
    ityp a54;
    ityp k1;
    ityp k2;
    ityp k3;
    ityp k4;
    ityp q1;
    ityp q2;
    ityp q3;
    ityp q4;
    ityp t1;
    ityp t2;
    ityp t3;
    ityp t4;
    ityp w1;
    ityp w2;
    ityp w3;
    ityp w4;
    ityp x1;
    ityp x2;
    ityp x3;
    ityp x4;

    a21 =   0.66667754298442;
    a31 =   0.63493935027993;
    a32 =   0.00342761715422;
    a41 = - 2.32428921184321;
    a42 =   2.69723745129487;
    a43 =   0.29093673271592;
    a51 =   0.25001351164789;
    a52 =   0.67428574806272;
    a53 = - 0.00831795169360;
    a54 =   0.08401868181222;

    q1 = 3.99956364361748;
    q2 = 1.64524970733585;
    q3 = 1.59330355118722;
    q4 = 0.26330006501868;

    t1 = t;
    x1 = x;
    w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h );
    k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1;

    t2 = t1 + a21 * h;
    x2 = x1 + a21 * k1;
    w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h );
    k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2;

    t3 = t1 + a31 * h  + a32 * h;
    x3 = x1 + a31 * k1 + a32 * k2;
    w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h );
    k3 = h * fv ( t3, x3 ) + h * gv ( t3, x3 ) * w3;

    t4 = t1 + a41 * h  + a42 * h  + a43 * h;
    x4 = x1 + a41 * k1 + a42 * k2 + a43 * k3;
    w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h );
    k4 = h * fv ( t4, x4 ) + h * gv ( t4, x4 ) * w4;

    return x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4;
}

#endif
