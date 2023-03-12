#ifndef __DISABLEDEEP_FEM1DBVPLINEAR

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_zero_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    05 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
	int *a = ( int * ) malloc ( n * sizeof ( int ) );
	for (dim_typ i = 0; i < n; ++i)
	a[i] = 0;
	return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _fem1d_bvp_linear ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEM1D_BVP_LINEAR solves a two point boundary value problem.
  Location:
    http://people.sc.fsu.edu/~jburkardt/c_src/fem1d_bvp_linear/fem1d_bvp_linear.c
  Discussion:
    The program uses the finite element method, with piecewise linear basis
    functions to solve a boundary value problem in one dimension.
    The problem is defined on the region 0 <= x <= 1.
    The following differential equation is imposed between 0 and 1:
      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
    where a(x), c(x), and f(x) are given functions.
    At the boundaries, the following conditions are applied:
      u(0.0) = 0.0
      u(1.0) = 0.0
    A set of N equally spaced nodes is defined on this
    interval, with 0 = X(1) < X(2) < ... < X(N) = 1.0.
    At each node I, we associate a piecewise linear basis function V(I,X),
    which is 0 at all nodes except node I.  This implies that V(I,X) is
    everywhere 0 except that
    for X(I-1) <= X <= X(I):
      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) )
    for X(I) <= X <= X(I+1):
      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )
    We now assume that the solution U(X) can be written as a linear
    sum of these basis functions:
      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)
    where U(X) on the left is the function of X, but on the right,
    is meant to indicate the coefficients of the basis functions.
    To determine the coefficient U(J), we multiply the original
    differential equation by the basis function V(J,X), and use
    integration by parts, to arrive at the I-th finite element equation:
        Integral A(X) * U'(X) * V'(I,X) + C(X) * U(X) * V(I,X) dx
      = Integral F(X) * V(I,X) dx
    We note that the functions U(X) and U'(X) can be replaced by
    the finite element form involving the linear sum of basis functions,
    but we also note that the resulting integrand will only be nonzero
    for terms where J = I - 1, I, or I + 1.
    By writing this equation for basis functions I = 2 through N - 1,
    and using the boundary conditions, we have N linear equations
    for the N unknown coefficients U(1) through U(N), which can
    be easily solved.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of nodes.
    Input, double A ( double X ), evaluates a(x);
    Input, double C ( double X ), evaluates c(x);
    Input, double F ( double X ), evaluates f(x);
    Input, double X[N], the mesh points.
    Output, double FEM1D_BVP_LINEAR[N], the finite element coefficients,
    which are also the value of the computed solution at the mesh points.
*/
{
	const dt3fitpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp (* a)(ityp) = s_data->a1;
	ityp (* c)(ityp) = s_data->a2;
	ityp (* f)(ityp) = s_data->a3;
	ityp * x = s_data->a4;
	
    # define QUAD_NUM 2

    ityp abscissa[QUAD_NUM] =
    {
        -0.577350269189625764509148780502,
        +0.577350269189625764509148780502
    };
    ityp *amat;
    ityp axq;
    ityp *b;
    ityp cxq;
    dim_typ e;
    dim_typ e_num;
    ityp fxq;
    dim_typ i;
    dim_typ ierror;
    dim_typ j;
    dim_typ l;
    dim_typ q;
    dim_typ quad_num = QUAD_NUM;
    dim_typ r;
    ityp *u;
    ityp weight[QUAD_NUM] =
    {
        1.00,
        1.00
    };
    ityp wq;
    ityp vl;
    ityp vlp;
    ityp vr;
    ityp vrp;
    ityp xl;
    ityp xq;
    ityp xr;
    /*
    Zero out the matrix and right hand side.
    */
    amat = r8mat_zero_new ( n, n );
    b = r8vec_zero_new ( n );

    e_num = n - 1;

    for ( e = 0; e < e_num; ++e )
    {
        l = e;
        r = e + 1;

        xl = x[l];
        xr = x[r];

        for ( q = 0; q < quad_num; ++q )
        {
            xq = ( ( 1.00 - abscissa[q] ) * xl   + ( 1.00 + abscissa[q] ) * xr ) /   2.00;

            wq = weight[q] * ( xr - xl ) / 2.00;

            vl = ( xr - xq ) / ( xr - xl );
            vlp =      - 1.00  / ( xr - xl );

            vr = ( xq - xl ) / ( xr - xl );
            vrp =  + 1.00      / ( xr - xl );

            axq = a ( xq );
            cxq = c ( xq );
            fxq = f ( xq );

            amat[l+l*n] += wq * ( vlp * axq * vlp + vl * cxq * vl );
            amat[l+r*n] += wq * ( vlp * axq * vrp + vl * cxq * vr );
            b[l]        += wq * ( vl * fxq );

            amat[r+l*n] += wq * ( vrp * axq * vlp + vr * cxq * vl );
            amat[r+r*n] += wq * ( vrp * axq * vrp + vr * cxq * vr );
            b[r]        += wq * ( vr * fxq );
        }
    }
    /*
    Equation 1 is the left boundary condition, U(0.0) = 0.0;
    */
    for ( j = 0; j < n; ++j )
        amat[0+j*n] = 0.00;
    b[0] = 0.00;
    for ( i = 1; i < n; ++i )
    b[i] -= amat[i+0*n] * b[0];
    for ( i = 0; i < n; ++i )
        amat[i+0*n] = 0.00;
    amat[0+0*n] = 1.00;
    /*
    Equation N is the right boundary condition, U(1.0) = 0.0;
    */
    for ( j = 0; j < n; ++j )
        amat[n-1+j*n] = 0.00;
    b[n-1] = 0.00;
    for ( i = 0; i < n - 1; ++i )
        b[i] -= amat[i+(n-1)*n] * b[n-1];
    for ( i = 0; i < n; ++i )
        amat[i+(n-1)*n] = 0.00;
    amat[n-1+(n-1)*n] = 1.00;
    /*
    Solve the linear system.
    */
    u = r8mat_solve2 ( n, amat, b, &ierror );
    free ( amat );
    free ( b );
    return u;
    # undef QUAD_NUM
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _h1s_error_linear ( void * data)
/******************************************************************************/
/*
  Purpose:
    H1S_ERROR_LINEAR estimates the seminorm error of a finite element solution.
  Discussion:
    We assume the finite element method has been used, over an interval [A,B]
    involving N nodes, with piecewise linear elements used for the basis.
    The coefficients U(1:N) have been computed, and a formula for the
    exact derivative is known.
    This function estimates the seminorm of the error:
      SEMINORM = Integral ( A <= X <= B ) ( dU(X)/dx - EXACT_UX(X) )^2 dX
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of nodes.
    Input, double X(N), the mesh points.
    Input, double U(N), the finite element coefficients.
    Input, function EQ = EXACT_UX ( X ), returns the value of the exact
    derivative at the point X.
    Output, double H1S_ERROR_LINEAR, the estimated seminorm of
    the error.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pitfit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * u = s_data->a2;
	ityp (* exact_ux)(ityp) = s_data->a3;
	
    # define QUAD_NUM 2

    ityp abscissa[QUAD_NUM] =
    {
        -0.577350269189625764509148780502,
        +0.577350269189625764509148780502
    };
    ityp exq;
    ityp h1s;
    dim_typ i;
    dim_typ q;
    dim_typ quad_num = QUAD_NUM;
    ityp ul;
    ityp ur;
    ityp uxq;
    ityp weight[QUAD_NUM] =
    {
        1.0,
        1.0
    };
    ityp wq;
    ityp xl;
    ityp xq;
    ityp xr;

    h1s = 0.00;
    /*
    Integrate over each interval.
    */
    for ( i = 0; i < n - 1; ++i )
    {
        xl = x[i];
        xr = x[i+1];
        ul = u[i];
        ur = u[i+1];

        for ( q = 0; q < quad_num; ++q )
        {
            xq = ( ( 1.00 - abscissa[q] ) * xl+ ( 1.00 + abscissa[q] ) * xr )/   2.00;
            wq = weight[q] * ( xr - xl ) / 2.00;
            /*
            The piecewise linear derivative is a constant in the interval.
            */
            uxq = ( ur - ul ) / ( xr - xl );
            exq = exact_ux ( xq );
            h1s += wq * pow ( uxq - exq, 2);
        }
    }
    
    result = sqrt ( h1s ); 
    return &result;
    # undef QUAD_NUM
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l1_error ( void * data)
/******************************************************************************/
/*
  Purpose:
    L1_ERROR estimates the l1 error norm of a finite element solution.
  Discussion:
    We assume the finite element method has been used, over an interval [A,B]
    involving N nodes.
    The coefficients U(1:N) have been computed, and a formula for the
    exact solution is known.
    This function estimates the little l1 norm of the error:
      L1_NORM = sum ( 1 <= I <= N ) abs ( U(i) - EXACT(X(i)) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of nodes.
    Input, double X[N], the mesh points.
    Input, double U[N], the finite element coefficients.
    Input, function EQ = EXACT ( X ), returns the value of the exact
    solution:
    We assume the finite element method has been used, over an interval [A,B]
    involving N nodes.
    The cion at the point X.
    Output, double L1_ERROR, the little l1 norm of the error.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pitfit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * u = s_data->a2;
	ityp (* exact)(ityp) = s_data->a3;
	
    ityp e1 = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        e1 += fabs ( u[i] - exact ( x[i] ) );
        
    result = e1 / ( ityp ) n; 
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _l2_error_linear ( void * data)
/******************************************************************************/
/*
  Purpose:
    L2_ERROR_LINEAR estimates the L2 error norm of a finite element solution.
  Discussion:
    We assume the finite element method has been used, over an interval [A,B]
    involving N nodes, with piecewise linear elements used for the basis.
    The coefficients U(1:N) have been computed, and a formula for the
    exact solution is known.
    This function estimates the L2 norm of the error:
      L2_NORM = Integral ( A <= X <= B ) ( U(X) - EXACT(X) )^2 dX
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of nodes.
    Input, double X[N], the mesh points.
    Input, double U[N], the finite element coefficients.
    Input, function EQ = EXACT ( X ), returns the value of the exact
    solution at the point X.
    Output, double L2_ERROR_LINEAR, the estimated L2 norm of the error.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pitfit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * u = s_data->a2;
	ityp (* exact)(ityp) = s_data->a3;
	
    # define QUAD_NUM 2

    ityp abscissa[QUAD_NUM] =
    {
        -0.577350269189625764509148780502,
        +0.577350269189625764509148780502
    };
    ityp e2;
    ityp eq;
    dim_typ i;
    dim_typ q;
    dim_typ quad_num = QUAD_NUM;
    ityp ul;
    ityp ur;
    ityp uq;
    ityp weight[QUAD_NUM] =
    {
        1.00,
        1.00
    };
    ityp wq;
    ityp xl;
    ityp xq;
    ityp xr;

    e2 = 0.00;
    /*
    Integrate over each interval.
    */
    for ( i = 0; i < n - 1; ++i )
    {
        xl = x[i];
        xr = x[i+1];
        ul = u[i];
        ur = u[i+1];

        for ( q = 0; q < quad_num; ++q )
        {
            xq = ( ( 1.00 - abscissa[q] ) * xl   + ( 1.00 + abscissa[q] ) * xr ) /   2.00;
            wq = weight[q] * ( xr - xl ) / 2.00;
            /*
            Use the fact that U is a linear combination of piecewise linears.
            */
            uq = ( ( xr - xq      ) * ul + (      xq - xl ) * ur ) / ( xr      - xl );
            eq = exact ( xq );
            e2 += wq * pow ( uq - eq, 2 );
        }
    }
    
    result = sqrt ( e2 );
    return &result;
    # undef QUAD_NUM
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _max_error_linear ( void * data)
/******************************************************************************/
/*
  Purpose:
    MAX_ERROR_LINEAR estimates the MAX error norm of a finite element solution.
  Discussion:
    We assume the finite element method has been used, over an interval [A,B]
    involving N nodes, with piecewise linear elements used for the basis.
    The coefficients U(1:N) have been computed, and a formula for the
    exact solution is known.
    This function estimates the MAX norm of the error:
      MAX_NORM = Integral ( A <= X <= B ) MAX ( abs ( U(X) - EXACT(X) ) ) dX
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 July 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of nodes.
    Input, double X[N], the mesh points.
    Input, double U[N], the finite element coefficients.
    Input, function EQ = EXACT ( X ), returns the value of the exact
    solution at the point X.
    Output, double MAX_ERROR_LINEAR, the estimated MAX norm of the error.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pitfit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * u = s_data->a2;
	ityp (* exact)(ityp) = s_data->a3;
	
    dim_typ e;
    dim_typ e_num;
    ityp eq;
    dim_typ i;
    dim_typ l;
    dim_typ q;
    dim_typ quad_num = 8;
    dim_typ r;
    ityp ul;
    ityp ur;
    ityp uq;
    ityp value;
    ityp wq;
    ityp xl;
    ityp xq;
    ityp xr;

    value = 0.00;
    /*
    Integrate over each interval.
    */
    e_num = n - 1;

    for ( e = 0; e < e_num; ++e )
    {
        l = e;
        xl = x[l];
        ul = u[l];

        r = e + 1;
        xr = x[r];
        ur = u[r];

        for ( q = 0; q < quad_num; q++ )
        {
            xq = ( ( ityp ) ( quad_num - q ) * xl+ ( ityp ) (            q ) * xr )/ ( ityp ) ( quad_num );
            /*
            Use the fact that U is a linear combination of piecewise linears.
            */
            uq = ( ( xr - xq      ) * ul+ (      xq - xl ) * ur )/ ( xr      - xl );
            eq = exact ( xq );
            value = MAX ( value, fabs ( uq - eq ) );
        }
    }
    /*
    For completeness, check last node.
    */
    /*
    Integral approximation requires multiplication by interval length.
    */
    
    result = MAX ( value, fabs ( x[n-1] - exact(u[n-1]) ) ) * ( x[n-1] - x[0] );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_solve2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_SOLVE2 computes the solution of an N by N linear system.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
    The linear system may be represented as
      A*X = B
    If the linear system is singular, but consistent, then the routine will
    still produce a solution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of equations.
    Input/output, double A[N*N].
    On input, A is the coefficient matrix to be inverted.
    On output, A has been overwritten.
    Input/output, double B[N].
    On input, B is the right hand side of the system.
    On output, B has been overwritten.
    Output, double R8MAT_SOLVE2[N], the solution of the linear system.
    Output, int *IERROR.
    0, no error detected.
    1, consistent singularity.
    2, inconsistent singularity.
*/
{
	const dt2pitpdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	dim_typ * ierror = s_data->a3; 
	
    ityp amax;
    dim_typ i;
    dim_typ imax;
    dim_typ j;
    dim_typ k;
    int *piv;
    ityp *x;

    *ierror = 0;
    piv = i4vec_zero_new ( n );
    x = r8vec_zero_new ( n );
    /*
    Process the matrix.
    */
    for ( k = 1; k <= n; ++k )
    {
    /*
    In column K:
    Seek the row IMAX with the properties that:
    IMAX has not already been used as a pivot;
    A(IMAX,K) is larger in magnitude than any other candidate.
    */
        amax = 0.00;
        imax = 0;
        for ( i = 1; i <= n; ++i )
        {
            if ( piv[i-1] == 0 )
            {
                if ( amax < fabs ( a[i-1+(k-1)*n] ) )
                {
                    imax = i;
                    amax = fabs ( a[i-1+(k-1)*n] );
                }
            }
        }
        /*
        If you found a pivot row IMAX, then,
        eliminate the K-th entry in all rows that have not been used for pivoting.
        */
        if ( imax != 0 )
        {
            piv[imax-1] = k;
            for ( j = k+1; j <= n; ++j )
                a[imax-1+(j-1)*n] /= a[imax-1+(k-1)*n];
            b[imax-1] /= a[imax-1+(k-1)*n];
            a[imax-1+(k-1)*n] = 1.00;

            for ( i = 1; i <= n; i++ )
            {
                if ( piv[i-1] == 0 )
                {
                    for ( j = k+1; j <= n; ++j )
                        a[i-1+(j-1)*n] -= a[i-1+(k-1)*n] * a[imax-1+(j-1)*n];
                    b[i-1] -= a[i-1+(k-1)*n] * b[imax-1];
                    a[i-1+(k-1)*n] = 0.00;
                }
            }
        }
    }
    /*
    Now, every row with nonzero PIV begins with a 1, and
    all other rows are all zero.  Begin solution.
    */
    for ( j = n; 1 <= j; --j )
    {
        imax = 0;
        for ( k = 1; k <= n; ++k )
            if ( piv[k-1] == j )
                imax = k;

        if ( !imax )
        {
            x[j-1] = 0.00;
            *ierror = b[j-1] == 0.00 ? 1:2;
        }
        else
        {
            x[j-1] = b[imax-1];

            for ( i = 1; i <= n; ++i )
                if ( i != imax )
                    b[i-1] -= a[i-1+(j-1)*n] * x[j-1];
        }
    }

    free ( piv );
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_zero_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_ZERO_NEW creates and zeroes an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, ityp r8VEC_ZERO_NEW[N], a vector of zeroes.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        a[i] = 0.00;
    return a;
}

#endif
