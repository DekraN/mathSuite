#ifndef __DISABLEDEEP_FEMBASIS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fem_basis_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEM_BASIS_1D evaluates an arbitrary 1D basis function.
  Discussion:
    Given the maximum degree D for the polynomial basis defined
    on a reference interval, we have D + 1 monomials
    of degree at most D.  In each barycentric coordinate, we define
    D+1 points, so that 0 <= I, J <= D and I+J = D, with
 (I,J) corresponding to
    * the basis point X(I,J) = ( I/D );
    * the basis monomial P(I,J)(X) = X^I.
    For example, with D = 2, we have simply:
      A---B---C
    with
       I J    X      P(I,J)(X)
    A (0 2) ( 0.0 )  1
    B (1 1) ( 0.5 )  x
    C (2 0) ( 1.0 )  x^2
    Now instead of the monomials P(I,J)(X), we want a set of
    polynomials L(I,J)(X) which span the same space, but have
    the Lagrange property, namely L(I,J) (X) is 1 if X is
    equal to X(I,J), and 0 if X is equal to any other
    of the basis points.
    This is easily arranged.  Given an index (I,J), we compute
    1) I factors of the form (   X -0/D) * (   X -1/D) * ... * (   X -(I-1)/D);
    2) J factors of the form ((1-X)-0/D) * ((1-X)-1/D) * ... * ((1-X)-(J-1)/D).
    This results in the product of I+J linear factors, in other words,
    a polynomial of degree D.  This polynomial is 0 at all basis points
    except X(I,J).  If we divide this polynomial by its value at
    the basis point, we arrive at the desired Lagrange polynomial
    L(I,J)(X).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the integer barycentric coordinates of
    the basis function, 0 <= I, J.  The polynomial degree D = I + J.
    Input, double X, the evaluation point.
    Output, double FEM_BASIS_1D, the value of the basis function at X.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	register dim_typ i = s_data->a0;
	const register dim_typ j = s_data->a1;
	const register ityp x = s_data->a2;
	
    ityp c;
    dim_typ d=i+j;
    ityp lij;
    dim_typ p;
    ityp w;

    lij = c = 1.00;

    for ( p = 0; p <= i - 1; ++i )
    {
        lij *= ( d * x - p );
        c *= (     i - p );
    }
    w = 1.00 - x;
    for ( p = 0; p <= j - 1; ++p )
    {
        lij*= ( d * w - p );
        c *= (     j - p );
    }
    
    result = lij/c; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fem_basis_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEM_BASIS_2D evaluates an arbitrary triangular basis function.
  Discussion:
    Given the maximum degree D for the polynomial basis defined
    on a reference triangle, we have ( ( D + 1 ) * ( D + 2 ) ) / 2 monomials
    of degree at most D.  In each barycentric coordinate, we define
    D+1 planes, so that 0 <= I, J, K <= D and I+J+K = D, with
 (I,J,K) corresponding to
    * the basis point (X,Y)(I,J,K) = ( I/D, J/D );
    * the basis monomial P(I,J,K)(X,Y) = X^I Y^J.
    For example, with D = 2, we have simply:
    F
    |\
    C-E
    |\|\
    A-B-D
    with
       I J K    X    Y    P(I,J,K)(X,Y)
    A (0 0 2) (0.0, 0.0)  1
    B (1 0 1) (0.5, 0.0)  x
    C (0 1 1) (0.0, 0.5)  y
    D (2 0 0) (1.0, 0.0)  x^2
    E (1 1 0) (0.5, 0.5)  x y
    F (0 2 0) (0.0, 1.0)  y^2
    Now instead of the monomials P(I,J,K)(X,Y), we want a set of
    polynomials L(I,J,K)(X,Y) which span the same space, but have
    the Lagrange property, namely L(I,J,K) (X,Y) is 1 if (X,Y) is
    equal to (X,Y)(I,J,K), and 0 if (X,Y) is equal to any other
    of the basis points.
    This is easily arranged.  Given an index (I,J,K), we compute
    1) I factors of the form (X-0)   * (X-1/D)   * ... * (X-(I-1)/D);
    2) J factors of the form (Y-0)   * (Y-1/D)   * ... * (Y-(J-1)/D);
    3) K factors of the form ((1-X-Y)-0/D) * ((1-X-Y)-1/D) * ...
       * ((1-X-Y)-(K-1)/D).
    This results in the product of I+J+K linear factors, in other words,
    a polynomial of degree D.  This polynomial is 0 at all basis points
    except (X,Y)(I,J,K).  If we divide this polynomial by its value at
    the basis point, we arrive at the desired Lagrange polynomial
    L(I,J,K)(X,Y).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, K, the integer barycentric coordinates of
    the basis function, 0 <= I, J, K.  The polynomial degree D = I + J + K.
    Input, double X, Y, the evaluation point.
    Output, double FEM_BASIS_2D, the value of the basis function at (X,Y).
*/
{
	static ityp result = MAX_VAL;
	
	const _3dt2it * const s_data = data;
	const register dim_typ i = s_data->a0;
	const register dim_typ j = s_data->a1;
	const register dim_typ k = s_data->a2;
	const register ityp x = s_data->a3;
	const register ityp y = s_data->a4;
	
    ityp c;
    dim_typ d = i+j+k;
    ityp lijk;
    dim_typ p;
    ityp w;

    lijk = c = 1.00;

    for ( p = 0; p <= i - 1; ++p )
    {
        lijk *= ( d * x - p );
        c *= (     i - p );
    }
    for ( p = 0; p <= j - 1; ++p )
    {
        lijk *= ( d * y - p );
        c *= (     j - p );
    }
    w = 1.00 - x - y;
    for ( p = 0; p <= k - 1; ++p )
    {
        lijk *= ( d * w - p );
        c *= (     k - p );
    }

	result = lijk/c;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _fem_basis_3d ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEM_BASIS_3D evaluates an arbitrary tetrahedral basis function.
  Discussion:
    Given the maximum degree D for the polynomial basis defined
    on a reference tetrahedron, we have
 ( D + 1 ) * ( D + 2 ) * ( D + 3 ) / 6 monomials
    of degree at most D.  In each barycentric coordinate, we define
    D+1 planes, so that 0 <= I, J, K, L <= D and I+J+K+L = D, with
 (I,J,K,L) corresponding to
    * the basis point (X,Y,Z)(I,J,K,L) = ( I/D, J/D, K/D );
    * the basis monomial P(I,J,K,L)(X,Y,Z) = X^I Y^J Z^K.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, K, L, the integer barycentric
    coordinates of the basis function, 0 <= I, J, K, L.
    The polynomial degree D = I + J + K + L.
    Input, double X, Y, Z, the evaluation point.
    Output, double FEM_BASIS_3D, the value of the basis function at (X,Y,Z).
*/
{
	static ityp result = MAX_VAL;
	
	const _4dt3it * const s_data = data;
	const register dim_typ i = s_data->a0;
	const register dim_typ j = s_data->a1;
	const register dim_typ k = s_data->a2;
	const register dim_typ l = s_data->a3;
	const register ityp x = s_data->a4;
	const register ityp y = s_data->a5;
	const register ityp z = s_data->a6;

    ityp c;
    dim_typ d = i+j+k+l;
    ityp lijkl;
    dim_typ p;
    ityp w;


    lijkl = c = 1.00;
    for ( p = 0; p <= i - 1; ++p)
    {
        lijkl *= ( d * x - p );
        c *= (     i - p );
    }
    for ( p = 0; p <= j - 1; ++p )
    {
        lijkl *= ( d * y - p );
        c *= (     j - p );
    }
    for ( p = 0; p <= k - 1; ++p )
    {
        lijkl *= ( d * z - p );
        c *= (     k - p );
    }
    w = 1.00 - x - y - z;
    for ( p = 0; p <= l - 1; ++p )
    {
        lijkl *= ( d * w - p );
        c *= (     l - p );
    }

	result = lijkl/c;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fem_basis_md ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEM_BASIS_MD evaluates an arbitrary M-dimensional basis function.
  Discussion:
    Given the maximum degree D for the polynomial basis defined
    on a reference tetrahedron, we have
 ( D + 1 ) * ( D + 2 ) * ( D + 3 ) / 6 monomials
    of degree at most D.  In each barycentric coordinate, we define
    D+1 planes, so that 0 <= I, J, K, L <= D and I+J+K+L = D, with
 (I,J,K,L) corresponding to
    * the basis point (X,Y,Z)(I,J,K,L) = ( I/D, J/D, K/D );
    * the basis monomial P(I,J,K,L)(X,Y,Z) = X^I Y^J Z^K.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int I[M+1], the integer barycentric coordinates of the
    basis function.  The polynomial degree D = sum ( I );
    Input, double X[M], the evaluation point.
    Output, double FEM_BASIS_MD, the value of the basis function at (X,Y,Z).
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitpi * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	ityp * x =  s_data->a1;
	int * i = s_data->a2;
	
    ityp c;
    dim_typ d = i4vec_sum ( m + 1, i );
    ityp l;
    dim_typ p, q;
    ityp w;

    for ( q = 0; q < m; ++q )
    {
        for ( p = 0; p < i[q]; ++p )
        {
            l *= ( d * x[q] - p );
            c *= (     i[q] - p );
        }
    }

    w = 1.00 - r8vec_sum ( m, x );

    for ( p = 0; p < i[m]; ++p )
    {
        l *= ( d * w    - p );
        c *= (     i[m] - p );
    }

	result = l/c;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fem_basis_prism_triangle ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEM_BASIS_PRISM_TRIANGLE evaluates a triangular prism basis function.
  Discussion:
    The element is a 3D prism, formed from a triangular base in the
    XY plane that is extended vertically in the Z direction.
    I[*] are the integer barycentric coordinates of a point in the
    triangle.  I[0] + I[1] + I[2] = DI, the degree of the triangular
    basis function BI.  X = I[0] / DI, Y = I[1] / DI.
    The triangle is assumed to be the unit reference
    triangle 0 <= X <= 1, 0 <= Y <= 1, 0 <= X + Y <= 1.
    J[*] are the integer barycentric coordinates of a point in the
    line segment.  J[0] + J[1] = DJ, the degree of the linear basis
    function BJ.  Z = J[0] / DJ.
    The line is assumed to be the unit line 0 <= Z <= 1.
    The degree of the basis function B = BI * BJ is D = DI + DJ.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int I[3], the integer barycentric coordinates of
    the triangular basis function, 0 <= I[*].
    The polynomial degree DI = I[0] + I[1] + I[2].
    Input, int J[2], the integer barycentric coordinates of
    the linear basis function, 0 <= J[*].
    The polynomial degree DJ = J[0] + J[1].
    Input, double XYZ[3], the evaluation point.
    Output, double B, the value of the basis function at XYZ.
*/
{
	static ityp result = MAX_VAL;
	
	const _2pdtpit * const s_data = data;
	dim_typ * i = s_data->a0;
	dim_typ * j = s_data->a1;
	ityp * xyz = s_data->a2;
	
	result = fem_basis_2d ( i[0], i[1], i[2], xyz[0], xyz[1] )*fem_basis_1d ( j[0], j[1], xyz[2] );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_fraction ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_FRACTION uses real arithmetic on an integer ratio.
  Discussion:
    Given integer variables I and J, both FORTRAN and C will evaluate
    an expression such as "I/J" using what is called "integer division",
    with the result being an integer.  It is often convenient to express
    the parts of a fraction as integers but expect the result to be computed
    using real arithmetic.  This function carries out that operation.
  Example:
       I     J   I/J  r8_FRACTION
       1     2     0  0.5
       7     4     1  1.75
       8     4     2  2.00
       9     4     2  2.25
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the arguments.
    Output, ityp r8_FRACTION, the value of the ratio.
*/
{
	static ityp result = MAX_VAL;
	
	int * const a_data = data;
	const register int i = a_data[0];
	const register int j = a_data[1];
	
	result = ( ityp ) ( i ) / ( ityp ) ( j );
    return &result;
}

#endif
