#ifndef __DISABLEDEEP_BROWNIANMOTIONSIMULATION

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _brownian_motion_simulation ( void * data)
/******************************************************************************/
/*
  Purpose:
    BROWNIAN_MOTION_SIMULATION simulates Brownian motion.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of time steps to take, plus 1.
    Input, double D, the diffusion coefficient.
    Input, double T, the total time.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double BROWNIAN_MOTION_SIMULATION[M*N], the initial position at
    time 0.0, and the N-1 successive locations of the particle.
*/
{
	const _2dt2itpitpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp d = s_data->a2;
	ityp t = s_data->a3;
	ityp * x = s_data->a4;
	int * seed = s_data->a5;
	
	ityp dt;
	ityp *dx;
	dim_typ i, j;
	ityp norm_dx, s;

	/*
	Set the time step.
	*/
	dt = t / ( ityp ) ( n - 1 );
	/*
	Start at the origin.
	*/
	#pragma omp parallel for
	for ( i = 0; i < m; ++i)
		x[i+0*m] = 0.0;
	/*
	Take N - 1 steps.
	*/
	for ( j = 1; j < n; ++j )
	{
		/*
		S is the stepsize.
		*/
		s = sqrt ( d * dt ) * r8_normal_01 ( seed );
		/*
		Direction DX is random, unit norm.
		*/
		dx = r8vec_normal_01_new ( m, seed );
		norm_dx = 0.00;
		for ( i = 00; i < m; ++i )
			norm_dx = norm_dx + pow ( dx[i], 2 );
		norm_dx = sqrt ( norm_dx );
		for ( i = 0; i < m; ++i )
			dx[i] *= s / norm_dx;
	/*
	Add the step to the current position.
	*/
		for ( i = 0; i < m; ++i )
	  		x[i+j*m] = x[i+(j-1)*m] + dx[i];
	}
	return NULL;
}

#endif
