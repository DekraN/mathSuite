#ifndef __DISABLEDEEP_POISSONOPENMP

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _sweep ( void * data)
/******************************************************************************/
/*
  Purpose:
   SWEEP carries out one step of the Jacobi iteration.
  Discussion:
    Assuming DX = DY, we can approximate
      - ( d/dx d/dx + d/dy d/dy ) U(X,Y)
    by
  ( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) - 4*U(i,j) ) / dx / dy
    The discretization employed below will not be correct in the general
    case where DX and DY are not equal.  It's only a little more complicated
    to allow DX and DY to be different, but we're not going to worry about
    that right now.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 December 2011
  Author:
    John Burkardt
  Parameters:
    Input, int NX, NY, the X and Y grid dimensions.
    Input, double DX, DY, the spacing between grid points.
    Input, double F[NX][NY], the right hand side data.
    Input, int ITOLD, the iteration index on input.
    Input, int ITNEW, the desired iteration index
    on output.
    Input, double U[NX][NY], the solution estimate on
    iteration ITNEW-1.
    Input/output, double UNEW[NX][NY], on input, the solution
    estimate on iteration ITOLD.  On output, the solution estimate on
    iteration ITNEW.
*/
{
	const _2dt2itppit2dt2ppit * const s_data = data;
	dim_typ nx = s_data->a0;
	dim_typ ny = s_data->a1;
	ityp dx = s_data->a2;
	ityp dy = s_data->a3;
	ityp ** f = s_data->a4;
	dim_typ itold = s_data->a5;
	dim_typ itnew = s_data->a6;
	ityp ** u = s_data->a7;
	ityp ** unew = s_data->a8;
	
    dim_typ i;
    dim_typ it;
    dim_typ j;

    # pragma omp parallel shared ( dx, dy, f, itnew, itold, nx, ny, u, unew ) private ( i, it, j )

    for ( it = itold + 1; it <= itnew; ++it )
    {
        /*
        Save the current estimate.
        */
        # pragma omp for
        for ( j = 0; j < ny; ++j )
            for ( i = 0; i < nx; ++i )
                u[i][j] = unew[i][j];
        /*
        Compute a new estimate.
        */
        # pragma omp for
        for ( j = 0; j < ny; ++j )
            for ( i = 0; i < nx; ++i )
                unew[i][j] = i == 0 || j == 0 || i == nx - 1 || j == ny - 1 ? f[i][j] : 0.25 * ( u[i-1][j] + u[i][j+1] + u[i][j-1] + u[i+1][j] + f[i][j] * dx * dy );

    }
    return NULL;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _u_exact ( void * data)
/******************************************************************************/
/*
  Purpose:
    U_EXACT evaluates the exact solution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, double X, Y, the coordinates of a point.
    Output, double U_EXACT, the value of the exact solution
    at (X,Y).
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp y = a_data[1];
	
	result = sin ( M_PI * x * y );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _uxxyy_exact ( void * data)
/******************************************************************************/
/*
  Purpose:
    UXXYY_EXACT evaluates ( d/dx d/dx + d/dy d/dy ) of the exact solution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, double X, Y, the coordinates of a point.
    Output, double UXXYY_EXACT, the value of
 ( d/dx d/dx + d/dy d/dy ) of the exact solution at (X,Y).
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp y = a_data[1];
	
	result = - M_PI * M_PI * ( x * x + y * y ) * sin ( M_PI * x * y );
    return &result;
}
# undef NX
# undef NY

#endif
