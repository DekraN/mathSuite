#ifndef __DISABLEDEEP_LATINRANDOM

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_uniform_01_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_UNIFORM_01_NEW fills an R8MAT with pseudorandom values scaled to [0,1].
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
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
    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
	const _2dtpi * const s_data = data; 
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
	dim_typ i, j;
	int k;
	ityp *r = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );
	
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < m; ++i )
		{
			k = *seed / 127773;
			
			*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
			
			if ( *seed < 0 )
				*seed += i4_huge;
			r[i+j*m] = ( ityp ) ( *seed ) * 4.656612875E-10;
		}
	
	return r;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _latin_random_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    LATIN_RANDOM_NEW returns points in a Latin Random square.
  Discussion:
    In each spatial dimension, there will be exactly one
    point whose coordinate value lies between consecutive
    values in the list:
  ( 0, 1, 2, ..., point_num ) / point_num
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 April 2003
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int POINT_NUM, the number of points.
    Input/output, int *SEED, a seed for UNIFORM.
    Output, double LATIN_RANDOM_NEW[DIM_NUM,POINT_NUM], the points.
*/
{
	const _2dtpi * const s_data = data; 
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ point_num = s_data->a1;
	int * seed = s_data->a2;
	
    dim_typ i, j;
    int *perm;
    ityp r;
    ityp * x = r8mat_uniform_01_new ( dim_num, point_num, seed );
    /*
    For spatial dimension I,
    pick a random permutation of 1 to POINT_NUM,
    force the corresponding I-th components of X to lie in the
    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
    */
    for ( i = 0; i < dim_num; ++i )
    {
        perm = perm_uniform_new ( point_num, seed );

        for ( j = 0; j < point_num; ++j )
            x[i+j*dim_num] = ( ( ( dim_typ ) perm[j] ) + x[i+j*dim_num] ) / ( ( dim_typ ) point_num );
        free ( perm );
    }
    return x;
}

#endif
