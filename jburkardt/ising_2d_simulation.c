#ifndef __DISABLEDEEP_ISING2DSIMULATION

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ising_2d_agree ( void * data)
/******************************************************************************/
/*
  Purpose:
    ISING_2D_AGREE returns the number of neighbors agreeing with each cell.
  Discussion:
    The count includes the cell itself, so it is between 1 and 5.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 Noveber 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of cells in each
    spatial dimension.
    Input, int C1[M*N], an array of 1's and -1's.
    Output, int C5[M*N], the number of neighbors
    that agree.  1, 2, 3, 4, or 5.
*/
{
	const _2dt2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * c1 = s_data->a2;
	int * c5 = s_data->a3;
	
    dim_typ i;
    int im;
    int ip;
    dim_typ j;
    int jm;
    int jp;

    for ( j = 0; j < n; ++j )
    {
        jp = i4_wrap ( j + 1, 0, n - 1 );
        jm = i4_wrap ( j - 1, 0, n - 1 );
        for ( i = 0; i < m; ++i )
        {
            ip = i4_wrap ( i + 1, 0, m - 1 );
            im = i4_wrap ( i - 1, 0, m - 1 );
            c5[i+j*m] = c1[i+j*m] + c1[ip+j*m] + c1[im+j*m] + c1[i+jm*m] + c1[i+jp*m];
            c5[i+j*m] = 0 < c1[i+j*m] ? ( 5 + c5[i+j*m] ) / 2 : ( 5 - c5[i+j*m] ) / 2;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ising_2d_initialize ( void * data)
/******************************************************************************/
/*
  Purpose:
    ISING_2D_INITIALIZE initializes the Ising array.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, double THRESH, the threshhold.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, in ISING_2D_INITIALIZE[M*N], the initial Ising array.
*/
{
	const _2dtitpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register ityp tresh = s_data->a2;
	int * seed = s_data->a3;
	
    int *c1;
    dim_typ i, j;
    ityp *r = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    r8mat_uniform_01 ( m, n, seed, r );

    c1 = ( int * ) malloc ( m * n * sizeof ( int ) );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < m; ++i )
            c1[i+j*m] = 1 - ((r[i+j*m] <= tresh)<<1);
    free ( r );

    return c1;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_uniform_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom values.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
    in column-major order.
    This routine implements the recursion
      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )
    The integer arithmetic never requires more than 32 bits,
    including a sign bit.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Reference:
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.
    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.
    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0, otherwise the output value of SEED
    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
    been updated.
    Output, double R[M*N], a matrix of pseudorandom values.
*/
{
	const _2dtpipit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	ityp * r = s_data->a3;
	
	dim_typ i, j;
	int k;
	
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < m; ++i )
		{
			k = *seed / 127773;
			
			*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
			
			if ( *seed < 0 )
				*seed += 2147483647;
			
			r[i+j*m] = ( ityp ) ( *seed ) * 4.656612875E-10;
		}
	
	return NULL;
}

#endif
