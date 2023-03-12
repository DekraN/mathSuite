#ifndef __DISABLEDEEP_PCEODEHERMITE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _pce_ode_hermite ( void * data)
/******************************************************************************/
/*
  Purpose:
    PCE_ODE_HERMITE applies the polynomial chaos expansion to a scalar ODE.
  Discussion:
    The deterministic equation is
      du/dt = - alpha * u,
      u(0) = u0
    In the stochastic version, it is assumed that the decay coefficient
    ALPHA is a Gaussian random variable with mean value ALPHA_MU and variance
    ALPHA_SIGMA^2.
    The exact expected value of the stochastic equation will be
      u(t) = u0 * exp ( t^2/2)
    This should be matched by the first component of the polynomial chaos
    expansion.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, double TI, TF, the initial and final times.
    Input, int  NT, the number of output points.
    Input, double UI, the initial condition.
    Input, int NP, the degree of the expansion.  Polynomials
    of degree 0 through NP will be used.
    Input, double ALPHA_MU, ALPHA_SIGMA, the mean and standard
    deviation of the decay coefficient.
    Output, double T[NT+1], U[(NT+1)*(NP+1)], the times and the PCE
    coefficients at the successive time steps.
*/
{
	const _2itdtitdt2it2pit * const s_data = data;
	ityp ti = s_data->a0;
	ityp tf = s_data->a1;
	const register dim_typ nt = s_data->a2;
	ityp ui = s_data->a3;
	const register dim_typ np = s_data->a4;
	ityp alpha_mu = s_data->a5;
	ityp alpha_sigma = s_data->a6;
	ityp * t = s_data->a7;
	ityp * u = s_data->a8;
	
    ityp dp;
    ityp dt;
    dim_typ i;
    dim_typ it;
    dim_typ j;
    dim_typ k;
    ityp t1;
    ityp t2;
    ityp term;
    ityp tp;
    ityp *u1;
    ityp *u2;

    u1 = ( ityp * ) malloc ( ( np + 1 ) * sizeof ( ityp ) );
    u2 = ( ityp * ) malloc ( ( np + 1 ) * sizeof ( ityp ) );

    dt = ( tf - ti ) / ( ityp ) ( nt );
    /*
    Set the PCE coefficients for the initial time.
    */
    t1 = ti;

    u1[0] = ui;
    for ( j = 1; j <= np; ++j )
        u1[j] = 0.0;
    /*
    Copy into the output arrays.
    */
    t[0] = t1;
    for ( j = 0; j <= np; ++j )
        u[0+j*(nt+1)] = u1[j];
    /*
    Time integration.
    */
    for ( it = 1; it <= nt; ++it)
    {
        t2 = ( ( ityp ) ( nt - it ) * ti   + ( ityp ) (      it ) * tf ) / ( ityp ) ( nt      );

        for ( k = 0; k <= np; ++k)
        {
        	dim_typ initializer1[2] =
        	{
        		k,
        		k
        	};
            dp = he_double_product_integral ( initializer1 );
            term = - alpha_mu * u1[k];

            i = 1;
            for ( j = 0; j <= np; ++j)
            {
            	dim_typ initializer2[3] =
            	{
            		i,
            		j,
            		k
            	};
                tp = he_triple_product_integral ( initializer2) ;
                term -= alpha_sigma * u1[j] * tp / dp;
            }
            u2[k] = u1[k] + dt * term;
        }
        /*
        Prepare for next step.
        */
        t1 = t2;
        for ( j = 0; j <= np; ++j )
            u1[j] = u2[j];
        /*
        Copy into the output arrays.
        */
        t[it] = t1;
        for ( j = 0; j <= np; ++j )
            u[it+j*(nt+1)] = u1[j];
    }

    free ( u1 );
    free ( u2 );

    return NULL;
}

#endif
