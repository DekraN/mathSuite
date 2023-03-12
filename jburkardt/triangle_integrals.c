#ifndef __DISABLEDEEP_TRIANGLEINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle01_monomial_integral_ ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE01_MONOMIAL_INTEGRAL: monomial integrals in the unit triangle.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    19 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the exponents.  
    Each exponent must be nonnegative.
    Output, double TRIANGLE01_MONOMIAL_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	int * const a_data = data;
	const register int i = a_data[0];
	const register int j = a_data[1];
	
	int k, l;
	ityp q = 1.00;
	
	for ( k = 0, l = 1; l <= i; ++l, ++k )
		q *= ( ityp ) ( l ) / ( ityp ) ( k );
	
	for ( l = 1; l <= j; ++l, ++k )
		q *= ( ityp ) ( l ) / ( ityp ) ( k );
	
	for ( l = 1; l <= 2; ++l, ++k )
		q /= ( ityp ) ( k );
		
	result = q;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle01_poly_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE01_POLY_INTEGRAL: polynomial integral over the unit triangle.
  Discussion:
    The unit triangle is T = ( (0,0), (1,0), (0,1) ).
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    21 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, integer D, the degree of the polynomial.
    Input, double P[M], the polynomial coefficients.
    M = ((D+1)*(D+2))/2.
    Output, double TRIANGLE01_POLY_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ d = s_data->a0;
	ityp * p = s_data->a1;
	
	int i;
	int j;
	int k;
	dim_typ km1;
	int m;
	ityp q;
	
	m = ( ( d + 1 ) * ( d + 2 ) ) / 2;
	
	q = 0.00;
	for ( km1 = 0; km1 < m; ++km1 )
	{
		k = km1 + 1;
		i4_to_pascal ( k, &i, &j );
		q += p[km1] * triangle01_monomial_integral_ ( i, j );
	}
	
	result = q;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _poly_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLY_PRODUCT computes P3(x,y) = P1(x,y) * P2(x,y) for polynomials.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int D1, the degree of factor 1.
    Input, double P1[M1], the factor 1 coefficients.
    M1 = ((D1+1)*(D1+2))/2.
    Input, int D2, the degree of factor 2.
    Input, double P2[M2], the factor2 coefficients.
    M2 = ((D2+1)*(D2+2))/2.
    Output, double POLY_PRODUCT[M3], the result coefficients.
    D3 = D1 + D2;
    M3 = ((D3+1)*(D3+2))/2.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ d1 = s_data->a0;
	const register dim_typ d2 = s_data->a1;
	ityp * p1 = s_data->a2;
	ityp * p2 = s_data->a3;
	
	int d3;
	int i1;
	int i2;
	int i3;
	int j1;
	int j2;
	int j3;
	int k1;
	dim_typ k1m1;
	int k2;
	dim_typ k2m1;
	int k3;
	dim_typ k3m1;
	int m1;
	int m2;
	int m3;
	double *p3;
	
	m1 = ( ( d1 + 1 ) * ( d1 + 2 ) ) / 2;
	m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
	/*
	Consider each entry in P1:
	P1(K1) * X^I1 * Y^J1
	and multiply it by each entry in P2:
	P2(K2) * X^I2 * Y^J2
	getting 
	P3(K3) = P3(K3) + P1(K1) * P2(X2) * X^(I1+I2) * Y(J1+J2)
	*/
	d3 = d1 + d2;
	m3 = ( ( d3 + 1 ) * ( d3 + 2 ) ) / 2;
	p3 = ( ityp * ) malloc ( m3 * sizeof ( ityp ) );
	
	for ( k3m1 = 0; k3m1 < m3; ++k3m1 )
		p3[k3m1] = 0.00;
	
	for ( k1m1 = 0; k1m1 < m1; ++k1m1)
	{
		k1 = k1m1 + 1;
		i4_to_pascal ( k1, &i1, &j1 );
		for ( k2m1 = 0; k2m1 < m2; ++k2m1)
		{
			k2 = k2m1 + 1;
			i4_to_pascal ( k2, &i2, &j2 );
			i3 = i1 + i2;
			j3 = j1 + j2;
			k3 = pascal_to_i4 ( i3, j3 );
			k3m1 = k3 - 1;
			p3[k3m1] += p1[k1m1] * p2[k2m1];
		}
	}
	
	return p3;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _poly_power_linear ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLY_POWER_LINEAR computes the polynomial ( a + b*x + c*y ) ^ n.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int D1, the degree of the linear polynomial,
    which should be 1 (or possibly 0).
    Input, double P1(M1), the coefficients of the linear polynomial.
    M1 = ( (D1+1)*(D1+2) ) / 2, which should be 3.
    Input, int N, the power to which the polynomial is to be 
    raised.  0 <= N.
    Output, double P2(M2), the coefficients of the power polynomial.
    D2 = N * D1;
    M2 = ( (D2+1)*(D2+2) ) / 2.
*/
{
	const _2dtpit * const s_data = data; 
	
	const register dim_typ d1 = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * p1 = s_data->a2;
	
	int d2;
	dim_typ i, j, k;
	int l;
	int m2;
	ityp *p2;
	
	d2 = n * d1;
	m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
	p2 = ( ityp * ) malloc ( m2 * sizeof ( ityp ) );
	
	if ( d1 == 0 )
	{
		p2[0] = pow ( p1[0], n );
		return p2;
	}
	
	if ( n == 0 )
	{
		p2[0] = 1.00;
		return p2;
	}
	/*
	Use the Trinomial formula.
	*/
	for ( i = 0; i <= n; ++i)
		for ( j = 0; j <= n - i; ++j )
			for ( k = 0; k <= n - i - j; ++k )
			{
				/*
				We store X^J Y^K in location L.
				*/
				l = pascal_to_i4 ( j, k );
				p2[l-1] = ( ityp ) ( trinomial ( i, j, k ) ) * pow ( p1[0], i ) * pow ( p1[1], j ) * pow ( p1[2], k );
			}
	
	return p2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rs_to_xy_map ( void * data)
/******************************************************************************/
/*
  Purpose:
    RS_TO_XY_MAP returns the linear map from reference to physical triangle.
  Discussion:
    This function returns the coefficients of the linear map that sends
    the vertices of the reference triangle, (0,0), (1,0) and (0,1), to
    the vertices of a physical triangle T, of the form:
      X = A + B * R + C * S;
      Y = D + E * R + F * S.
  Reference Element:
    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, double T[2,3], the coordinates of the vertices.  The
    vertices are assumed to be the images of (0,0), (1,0) and (0,1) 
    respectively.
    Output, double *A, *B, *C, *D, *E, *F, the mapping coefficients.
*/
{
	ityp ** const a_data = data;
	ityp * t = a_data[0];
	ityp * a = a_data[1];
	ityp * b = a_data[2];
	ityp * c = a_data[3];
	ityp * d = a_data[4];
	ityp * e = a_data[5];
	ityp * f = a_data[6];
	
	*a = t[0+0*2];
	*b = t[0+1*2] - t[0+0*2];
	*c = t[0+2*2] - t[0+0*2];
	
	*d = t[1+0*2];
	*e = t[1+1*2] - t[1+0*2];
	*f = t[1+2*2] - t[1+0*2];
	
	return NULL;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_area_from_vertex ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_AREA returns the area of a triangle.
  Discussion:
    If the vertices are given in counter clockwise order, the area
    will be positive.
  Licensing:
    This code is distributed under the GNU GPL license.
  Modified:
    21 April 2015
  Author:
    John Burkardt.
  Parameters:
    Input, double T[2*3], the vertices of the triangle.
    Output, double TRIANGLE_AREA, the area of the triangle.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * t = data;
	
	result = 0.50 * 	( ( t[0+1*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_monomial_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_MONOMIAL_INTEGRAL integrates a monomial over an arbitrary triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the exponents of X and Y in the monomial.
    0 <= I, J.
    Input, double T[2*3], the vertices of the triangle.
    Output, double TRIANGLE_MONOMIAL_INTEGRAL, the integral 
    of X^I * Y^J over triangle T.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data;
	const register dim_typ i = s_data->a0;
	const register dim_typ j = s_data->a1;
	ityp * t = s_data->a2;
	
	ityp a;
	ityp b;
	ityp c;
	ityp d;
	int d1;
	int d2;
	int d3;
	int d4;
	int d5;
	ityp e;
	ityp f;
	int m1; 
	int m2;
	ityp *p1;
	ityp *p2;
	ityp *p3;
	ityp *p4;
	ityp *p5;
	ityp q;
	/*
	Get map coefficients from reference RS triangle to general XY triangle.
	R = a+b*X+c*Y
	S = d+e*X+f*Y
	*/
	rs_to_xy_map ( t, &a, &b, &c, &d, &e, &f );
	/*
	Set
	P1(R,S) = a+b*R+c*S
	P2(R,S) = d+e*R+f*S
	*/
	d1 = 1;
	m1 = ( ( d1 + 1 ) * ( d1 + 2 ) ) / 2;
	p1 = ( ityp * ) malloc ( m1 * sizeof ( ityp ) );
	p1[0] = a;
	p1[1] = b;
	p1[2] = c;
	
	d2 = 1;
	m2 = ( ( d2 + 1 ) * ( d2 + 2 ) ) / 2;
	p2 = ( ityp * ) malloc ( m2 * sizeof ( ityp ) );
	p2[0] = d;
	p2[1] = e;
	p2[2] = f;
	/*
	Exponentiate:
	P3(R,S) = P1(R,S)^i
	P4(R,S) = P2(R,S)^j
	*/
	d3 = i * d1;
	p3 = poly_power_linear ( d1, p1, i );
	
	d4 = j * d2;
	p4 = poly_power_linear ( d2, p2, j );
	/*
	Compute the product 
	P5(R,S) = P3(R,S) * P4(R,S)
	*/
	d5 = d3 + d4;
	p5 = poly_product ( d3, p3, d4, p4 );
	/*
	Compute the integral of P5(R,S) over the reference triangle.
	*/
	q = triangle01_poly_integral ( d5, p5 );
	/*
	Multiply by the area of the physical triangle T(X,Y) divided by
	the area of the reference triangle.
	*/
	q *= triangle_area_from_vertex ( t ) / 0.50;
	
	free ( p1 );
	free ( p2 );
	free ( p3 );
	free ( p4 );
	free ( p5 );
	
	result = q;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_poly_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_POLY_INTEGRAL: polynomial integral over a triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int D, the degree of the polynomial.
    Input, double P[M], the polynomial coefficients.
    M = ((D+1)*(D+2))/2.
    Input, double T[2*3], the vertices of the triangle.
    Output, double TRIANGLE_POLY_INTEGRAL, the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ d = s_data->a0;
	ityp * p = s_data->a1;
	ityp * t = s_data->a2;
	
    int i, j;
    dim_typ k;
    dim_typ km1;
    dim_typ m;
    ityp q;

    m = ( ( d + 1 ) * ( d + 2 ) ) / 2;
    q = 0.00;
    for ( km1 = 0; km1 <= m; ++km1)
    {
        k = km1 + 1;
        i4_to_pascal ( k, &i, &j );
        q += p[km1] * triangle_monomial_integral ( i, j, t );
    }

	result = q;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_xy_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_XY_INTEGRAL computes the integral of XY over a triangle.
  Discussion:
    This function was written as a special test case for the general
    problem of integrating a monomial x^alpha * y^beta over a general
    triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of the
    triangle vertices.
    Output, double TRIANGLE_XY_INTEGRAL, the integral of X*Y
    over the triangle.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	ityp x1 = a_data[0];
	ityp y1 = a_data[1];
	ityp x2 = a_data[2];
	ityp y2 = a_data[3];
	ityp x3 = a_data[4];
	ityp y3 = a_data[5];
	
    ityp det;
    ityp p00;
    ityp p01;
    ityp p02;
    ityp p10;
    ityp p11;
    ityp p20;
    ityp q;
    /*
    x = x1 * ( 1.0 - xi - eta )
    + x2 *         xi
    + x3 *              eta;

    y = y1 * ( 1.0 - xi - eta )
    + y2 *         xi
    + y3 *              eta;

    Rewrite as linear polynomials in (xi,eta):

    x = x1 + ( x2 - x1 ) * xi + ( x3 - x1 ) * eta
    y = y1 + ( y2 - y1 ) * xi + ( y3 - y1 ) * eta

    Jacobian:

    J = [ ( x2 - x1 ) ( x3 - x1 ) ]
    [ ( y2 - y1 ) ( y3 - y1 ) ]

    det J = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )

    Integrand

    x * y = ( x1 + ( x2 - x1 ) * xi + ( x3 - x1 ) * eta )
    * ( y1 + ( y2 - y1 ) * xi + ( y3 - y1 ) * eta )

    Rewrite as linear combination of monomials:

    x * y = 1      * x1 * y1
    + eta    * ( x1 * ( y3 - y1 ) + ( x3 - x1 ) * y1 )
    + xi     * ( x1 * ( y2 - y1 ) + ( x2 - x1 ) * y1 )
    + eta^2  * ( x3 - x1 ) * ( y3 - y1 )
    + xi*eta * ( ( x2 - x1 ) * ( y3 - y1 ) + ( x3 - x1 ) * ( y2 - y1 ) )
    + xi^2   * ( x2 - x1 ) * ( y2 - y1 )
    */
    det = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 );

    p00 = x1 * y1;

    p01 = x1 * ( y3 - y1 ) + ( x3 - x1 ) * y1;
    p10 = x1 * ( y2 - y1 ) + ( x2 - x1 ) * y1;

    p02 = ( x3 - x1 ) * ( y3 - y1 );
    p11 = ( x2 - x1 ) * ( y3 - y1 ) + ( x3 - x1 ) * ( y2 - y1 );
    p20 = ( x2 - x1 ) * ( y2 - y1 );

    q = 0.00;
    q += p00 * triangle01_monomial_integral_ ( 0, 0 );
    q += p10 * triangle01_monomial_integral_ ( 1, 0 );
    q += p01 * triangle01_monomial_integral_ ( 0, 1 );
    q += p20 * triangle01_monomial_integral_ ( 2, 0 );
    q += p11 * triangle01_monomial_integral_ ( 1, 1 );
    q += p02 * triangle01_monomial_integral_ ( 0, 2 );

	result = q * det;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _xy_to_rs_map ( void * data)
/******************************************************************************/
/*
  Purpose:
    XY_TO_RS_MAP returns the linear map from physical to reference triangle.
  Discussion:
    Given the vertices T of an arbitrary triangle in the (X,Y) coordinate
    system, this function returns the coefficients of the linear map
    that sends the vertices of T to (0,0), (1,0) and (0,1) respectively
    in the reference triangle with coordinates (R,S):
      R = A + B * X + C * Y;
      S = D + E * X + F * Y.
  Reference Element T3:
    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2015
  Author:
    John Burkardt
  Parameters
    Input, double T[2*3], the X and Y coordinates
    of the vertices.  The vertices are assumed to be the images of
 (0,0), (1,0) and (0,1) respectively.

    Output, double *A, *B, *C, *D, *E, *F, the mapping coefficients.
*/
{
	ityp ** const a_data = data;
	ityp * t = a_data[0];
	ityp * a = a_data[1];
	ityp * b = a_data[2];
	ityp * c = a_data[3];
	ityp * d = a_data[4];
	ityp * e = a_data[5];
	ityp * f = a_data[6];
	
    ityp g;

    g = ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2] - t[0+0*2] )
    - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) );

    *a = ( - ( t[1+2*2] - t[1+0*2] ) * t[0+0*2]
    + ( t[0+2*2] - t[0+0*2] ) * t[1+0*2] ) / g;

    *b = ( t[1+2*2] - t[1+0*2] ) / g;

    *c =   - ( t[0+2*2] - t[0+0*2] ) / g;

    *d = ( ( t[1+1*2] - t[1+0*2] ) * t[0+0*2]
    - ( t[0+1*2] - t[0+0*2] ) * t[1+0*2] ) / g;

    *e =   - ( t[1+1*2] - t[1+0*2] ) / g;

    *f = ( t[0+1*2] - t[0+0*2] ) / g;

    return NULL;
}

#endif
