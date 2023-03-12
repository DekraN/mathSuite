#ifndef __DISABLEDEEP_FD1DBVP

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _fd1d_bvp ( void * data)
/******************************************************************************/
/*
  Purpose:
    FD1D_BVP solves a two point boundary value problem.
  Discussion:
    The program uses the finite difference method to solve a BVP
 (boundary value problem) in one dimension.
    The problem is defined on the region X[0] <= x <= X[N-1].
    The following differential equation is imposed:
      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
    where a(x), c(x), and f(x) are given functions.  We write out
    the equation in full as
      - a(x) * u''(x) - a'(x) * u'(x) + c(x) * u(x) = f(x)
    At the boundaries, the following conditions are applied:
      u(X[0]) = 0.0
      u(X[N-1]) = 0.0
    We replace the function U(X) by a vector of N values U associated
    with the nodes.
    The first and last values of U are determined by the boundary conditions.
    At each interior node I, we write an equation to help us determine
    U(I).  We do this by approximating the derivatives of U(X) by
    finite differences.  Let us write XL, XM, and XR for X(I-1), X(I) and X(I+1).
    Similarly we have UL, UM, and UR.  Other quantities to be evaluated at
    X(I) = XM will also be labeled with an M:
      - AM * ( UL - 2 UM + UR ) / DX^2 - A'M * ( UL - UR ) / ( 2 * DX ) = FM

    These N-2 linear equations for the unknown coefficients complete the
    linear system and allow us to compute the finite difference approximation
    to the solution of the BVP.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of nodes.
    Input, double A ( double x ), evaluates a(x);
    Input, double APRIME ( double x ), evaluates a'(x);
    Input, double C ( double x ), evaluates c(x);
    Input, double F ( double x ), evaluates f(x);
    Input, double X[N], the mesh points, which may be nonuniformly spaced.
    Output, double FD1D_BVP[N], the value of the finite difference
    approximation to the solution.
*/
{
	const dt4fitpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp (* a)(ityp) = s_data->a1;
	ityp (* aprime)(ityp) = s_data->a2;
	ityp (* c)(ityp) = s_data->a3;
	ityp (* f)(ityp) = s_data->a4;
	ityp * x = s_data->a5;
	
    ityp am;
    ityp apm;
    ityp cm;
    ityp fm;
    dim_typ i;
    ityp *rhs = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    ityp *tri = ( ityp * ) malloc ( 3 * n * sizeof ( ityp ) );
    ityp *u;
    ityp x1;
    ityp x2;
    ityp xm;
    /*
    Equation 1 is the left boundary condition, U(X[0]) = 0.0;
    */

    tri[0+0*3] = 0.00;
    tri[1+0*3] = 1.00;
    tri[2+0*3] = 0.00;
    rhs[0] = 0.00;
    /*
    Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1),
    and so on.
    */
    for ( i = 1; i < n - 1; ++i )
    {
        xm  = x[i];
        am  = a ( xm );
        apm = aprime ( xm );
        cm  = c ( xm );
        fm  = f ( xm );

        tri[0+i*3] = - 2.00 * am / ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] )+ apm / ( x[i+1] - x[i-1] );
        tri[1+i*3] = + 2.00 * am / ( x[i] - x[i-1] ) / ( x[i+1] - x[i] )+ cm;
        tri[2+i*3] = - 2.00 * am / ( x[i+1] - x[i] ) / ( x[i+1] - x[i-1] )- apm / ( x[i+1] - x[i-1] );
        rhs[i]   = fm;
    }
        /*
        Equation N is the right boundary condition, U(X[N-1]) = 0.0;
        */
    tri[0+(n-1)*3] = 0.00;
    tri[1+(n-1)*3] = 1.00;
    tri[2+(n-1)*3] = 0.00;
    rhs[n-1] = 0.00;
    /*
    Solve the linear system.
    */
    u = r83np_fs ( n, tri, rhs );
    free ( rhs );
    free ( tri );
    return u;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r83np_fs ( void * data)
/******************************************************************************/
/*
  Purpose:
    R83NP_FS factors and solves an R83NP system.
  Discussion:
    The R83NP storage format is used for a tridiagonal matrix.
    The subdiagonal   is in entries (0,1:N-1),
    the diagonal      is in entries (1,0:N-1),
    the superdiagonal is in entries (2,0:N-2).
    This algorithm requires that each diagonal entry be nonzero.
    It does not use pivoting, and so can fail on systems that
    are actually nonsingular.
    The "R83NP" format used for this routine is different from the R83 format.
    Here, we insist that the nonzero entries
    for a given row now appear in the corresponding column of the
    packed array.
  Example:
    Here is how a R83 matrix of order 5 would be stored:
       *  A21 A32 A43 A54
      A11 A22 A33 A44 A55
      A12 A23 A34 A45  *
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the linear system.
    Input/output, double A[3*N].
    On input, the nonzero diagonals of the linear system.
    On output, the data in these vectors has been overwritten
    by factorization information.
    Input, double B[N], the right hand side.
    Output, double R83NP_FS[N], the solution of the linear system.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    dim_typ i;
    ityp *x;
    /*
    Check.
    */
    for ( i = 0; i < n; ++i)
        if ( !a[1+i*3])
        return  NULL;

    x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i)
        x[i] = b[i];

    for ( i = 1; i < n; ++i )
    {
        a[1+i*3] -= a[2+(i-1)*3] * a[0+i*3] / a[1+(i-1)*3];
        x[i]     -= x[i-1]       * a[0+i*3] / a[1+(i-1)*3];
    }

    x[n-1] = x[n-1] / a[1+(n-1)*3];
    for ( i = n-2; 0 <= i; --i )
        x[i] = ( x[i] - a[2+i*3] * x[i+1] ) / a[1+i*3];

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_even ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_EVEN returns an r8VEC of values evenly spaced between ALO and AHI.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 February 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of values.
    Input, ityp ALO, AHI, the low and high values.
    Output, ityp A[N], N evenly spaced values.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp alo = s_data->a1;
	const register ityp ahi = s_data->a2;
	ityp * a = s_data->a3;
	
    if ( n == 1 )
        a[0] = 0.50 * ( alo + ahi );
    else
        for (dim_typ i = 1; i <= n; ++i )
            a[i-1] = ( ( ityp ) ( n - i     ) * alo+ ( ityp ) (     i - 1 ) * ahi )/ ( ityp ) ( n     - 1 );

    return NULL;
}

#endif
