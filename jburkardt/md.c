#ifndef __DISABLEDEEP_MD

#include "../dutils.h"



/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _md_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMPUTE computes the forces and energies.
  Discussion:
    The computation of forces and energies is fully parallel.
    The potential function V(X) is a harmonic well which smoothly
    saturates to a maximum value at M_PI/2:
      v(x) = ( sin ( MIN ( x, M_PI_2 ) ) )**2
    The derivative of the potential is:
      dv(x) = 2.0 * sin ( MIN ( x, M_PI_2 ) ) * cos ( MIN ( x, M_PI_2 ) )
            = sin ( 2.0 * MIN ( x, M_PI_2 ) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2007
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    Input, int NP, the number of particles.
    Input, int ND, the number of spatial dimensions.
    Input, double POS[ND*NP], the position of each particle.
    Input, double VEL[ND*NP], the velocity of each particle.
    Input, double MASS, the mass of each particle.
    Output, double F[ND*NP], the forces.
    Output, double *POT, the total potential energy.
    Output, double *KIN, the total kinetic energy.
*/
{
	const _2dt2pitit3pit * const s_data = data;
	const register dim_typ np = s_data->a0;
	const register dim_typ nd = s_data->a1;
	ityp * pos = s_data->a2;
	ityp * vel = s_data->a3; 
	const register ityp mass = s_data->a4;
	ityp * f = s_data->a5;
	ityp * pot = s_data->a6;
	ityp * kin = s_data->a7;
	
    ityp d;
    ityp d2;
    dim_typ i, j, k;
    ityp ke;
    ityp pe;
    ityp rij[3];

    pe = ke = 0.00;
    
	/*
    # pragma omp parallel \
    shared ( f, nd, np, pos, vel ) \
    private ( i, j, k, rij, d, d2 )


    # pragma omp for reduction ( + : pe, ke )
    */
    for ( k = 0; k < np; ++k )
    {
        /*
        Compute the potential energy and forces.
        */
        for ( i = 0; i < nd; ++i )
            f[i+k*nd] = 0.00;

        for ( j = 0; j < np; ++j )
        {
            if ( k != j )
            {
                d = md_dist ( nd, pos+k*nd, pos+j*nd, rij );
                d2 = d<M_PI_2?d:M_PI_2;
                /*
                Attribute half of the potential energy to particle J.
                */

                pe += 0.50 * pow ( sin ( d2 ), 2 );

                for ( i = 0; i < nd; ++i )
                    f[i+k*nd] -= rij[i] * sin ( 2.00 * d2 ) / d;
            }
        }
        /*
        Compute the kinetic energy.
        */
        for ( i = 0; i < nd; ++i )
            ke += vel[i+k*nd] * vel[i+k*nd];
    }

    ke = ke * 0.50 * mass;
    *pot = pe;
    *kin = ke;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _md_dist ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIST computes the displacement (and its norm) between two particles.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2007
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    Input, int ND, the number of spatial dimensions.
    Input, double R1[ND], R2[ND], the positions of the particles.
    Output, double DR[ND], the displacement vector.
    Output, double D, the Euclidean norm of the displacement.
*/
{
	static ityp result = MAX_VAL;
	
	const dt3pit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * r1 = s_data->a1;
	ityp * r2 = s_data->a2; 
	ityp * dr = s_data->a3; 
	
    ityp d = 0.00;
    for (dim_typ i = 0; i < nd; ++i )
    {
        dr[i] = r1[i] - r2[i];
        d += dr[i] * dr[i];
    }

	result = sqrt ( d );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _md_initialize ( void * data)
/******************************************************************************/
/*
  Purpose:
    INITIALIZE initializes the positions, velocities, and accelerations.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2007
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    Input, int NP, the number of particles.
    Input, int ND, the number of spatial dimensions.
    Input, double BOX[ND], specifies the maximum position
    of particles in each dimension.
    Input, int *SEED, a seed for the random number generator.
    Output, double POS[ND*NP], the position of each particle.
    Output, double VEL[ND*NP], the velocity of each particle.
    Output, double ACC[ND*NP], the acceleration of each particle.
*/
{
	const _2dtpitpi3pit * const s_data = data;
	const register dim_typ np = s_data->a0;
	const register dim_typ nd = s_data->a1;
	ityp * box = s_data->a2;
	int * seed = s_data->a3; 
	ityp * pos = s_data->a4;
	ityp * vel = s_data->a5;
	ityp * acc = s_data->a6;
	
    dim_typ i, j;
    /*
    Give the particles random positions within the box.
    */

    for ( j = 0; j < np; ++j )
        for ( i = 0; i < nd; ++i )
        {
            vel[i+j*nd] = acc[i+j*nd] = 0.00;
            pos[i+j*nd] = box[i] * r8_uniform_01 ( seed );
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _md_update ( void * data)
/******************************************************************************/
/*
  Purpose:
    UPDATE updates positions, velocities and accelerations.
  Discussion:
    The time integration is fully parallel.
    A velocity Verlet algorithm is used for the updating.
    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
    a(t+dt) = f(t) / m
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    Input, int NP, the number of particles.
    Input, int ND, the number of spatial dimensions.
    Input/output, double POS[ND*NP], the position of each particle.
    Input/output, double VEL[ND*NP], the velocity of each particle.
    Input, double F[ND*NP], the force on each particle.
    Input/output, double ACC[ND*NP], the acceleration of each particle.
    Input, double MASS, the mass of each particle.
    Input, double DT, the time step.
*/
{
	const _2dt4pit2it * const s_data = data;
	const register dim_typ np = s_data->a0;
	const register dim_typ nd = s_data->a1;
	ityp * pos = s_data->a2;
	ityp * vel = s_data->a3; 
	ityp * f = s_data->a4;
	ityp * acc = s_data->a5;
	ityp mass = s_data->a6;
	ityp dt = s_data->a7;
	
    dim_typ i, j;
    ityp rmass = 1.00 / mass;

	/*
    # pragma omp parallel \
    shared ( acc, dt, f, nd, np, pos, rmass, vel ) \
    private ( i, j )

    # pragma omp for
    */
    for ( j = 0; j < np; j++ )
        for ( i = 0; i < nd; i++ )
        {
            pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.50 * acc[i+j*nd] * dt * dt;
            vel[i+j*nd] = vel[i+j*nd] + 0.50 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
            acc[i+j*nd] = f[i+j*nd] * rmass;
        }

    return NULL;
}

#endif
