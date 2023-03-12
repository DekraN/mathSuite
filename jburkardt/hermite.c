#ifndef __DISABLEDEEP_HERMITE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _dif_shift_x ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIF_SHIFT_X replaces one abscissa of a divided difference table with a new one.
  Discussion:
    This routine shifts the representation of a divided difference polynomial by
    dropping the last X value in XD, and adding a new X value to the
    beginning of the XD array, suitably modifying the coefficients stored
    in YD.
    The representation of the polynomial is changed, but the polynomial itself
    should be identical.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 June 2011
  Author:
    John Burkardt
  Reference:
    Carl deBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.
  Parameters:
    Input, int ND, the number of divided difference coefficients, and
    the number of entries in XD.
    Input/output, double XD[ND], the X values used in the representation of
    the divided difference polynomial.
    After a call to this routine, the last entry of XD has been dropped, the other
    entries have shifted up one index, and XV has been inserted at the
    beginning of the array.
    Input/output, double YD[ND], the divided difference coefficients
    corresponding to the XD array.  On output, this array has been
    adjusted.
    Input, double XV, a new X value which is to be used in the representation
    of the polynomial.  On output, XD[0] equals XV and the representation
    of the polynomial has been suitably changed.
    Note that XV does not have to be distinct from any of the original XD
    values.
*/
{
	const dtit2pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	ityp xv = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	
    dim_typ i;
    /*
    Recompute the divided difference coefficients.
    */
    for ( i = nd - 2; 0 <= i; --i )
        yd[i] += ( xv - xd[i] ) * yd[i+1];
    /*
    Shift the XD values up one position and insert XV.
    */
    for ( i = nd - 1; 0 < i; --i)
        xd[i] = xd[i-1];

    xd[0] = xv;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _dif_shift_zero ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIF_SHIFT_ZERO shifts a divided difference table so that all abscissas are zero.
  Discussion:
    When the abscissas are changed, the coefficients naturally
    must also be changed.
    The resulting pair (XD, YD) still represents the
    same polynomial, but the entries in YD are now the
    standard polynomial coefficients.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 June 2011
  Author:
    John Burkardt
  Reference:
    Carl deBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.
  Parameters:
    Input, int ND, the length of the XD and YD arrays.
    Input/output, double XD[ND], the X values that correspond to the
    divided difference table.  On output, XD contains only zeroes.
    Input/output, double YD[ND], the divided difference table
    for the polynomial.  On output, YD is also
    the coefficient array for the standard representation
    of the polynomial.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xd = s_data->a1;
	ityp * yd = s_data->a2;	
	
    ityp xv = 0.00;
    for (dim_typ i = 1; i <= nd; ++i )
        dif_shift_x ( nd, xd, yd, xv );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _dif_to_r8poly ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIF_TO_R8POLY converts a divided difference table to a standard polynomial.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 February 2011
  Author:
    John Burkardt
  Reference:
    Carl deBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.
  Parameters:
    Input, int ND, the number of coefficients, and abscissas.
    Input, double XD[ND], the X values used in the divided difference
    representation of the polynomial.
    Input, double YD[ND], the divided difference table.
    Output, double C[ND], the standard form polyomial coefficients.
    C[0] is the constant term, and C[ND-1] is the coefficient
    of X^(ND-1).
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xd = s_data->a1;
	ityp * yd = s_data->a2;	
	ityp * c = s_data->a3;
	
    dim_typ i, j;

    for ( i = 0; i < nd; ++i )
        c[i] = yd[i];
    /*
    Recompute the divided difference coefficients.
    */
    for ( j = 1; j <= nd - 1; ++j )
        for ( i = 1; i <= nd - j; ++i )
            c[nd-i-1] -= xd[nd-i-j] * c[nd-i];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _dif_vals ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIF_VALS evaluates a divided difference polynomial at a set of points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 May 2011
  Author:
    John Burkardt
  Reference:
    Carl deBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.
  Parameters:
    Input, int ND, the order of the difference table.
    Input, double XD[ND], the X values of the difference table.
    Input, double YD[ND], the divided differences.
    Input, int NV, the number of evaluation points.
    Input, double XV[NV], the evaluation points.
    Output, double DIF_VALS[NV], the value of the divided difference
    polynomial at the evaluation points.
*/
{
	const _2dt3pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ nv = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	ityp * xv = s_data->a4;
	
    dim_typ i, j;
    ityp *yv = (ityp * ) malloc ( nv * sizeof ( ityp ) );

    for ( j = 0; j < nv; ++j )
    {
        yv[j] = yd[nd-1];
        for ( i = 2; i <= nd; ++i )
            yv[j] = yd[nd-i] + ( xv[j] - xd[nd-i] ) * yv[j];
    }
    return yv;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hermite_basis_0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_BASIS_0 evaluates a zero-order Hermite interpolation basis function.
  Discussion:
    Given ND points XD, with values YD and derivative values YPD, the
    Hermite interpolant can be written as:
      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
    and H1(I;X) is the I-th first order Hermite interpolation basis function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of abscissas.
    Input, double X[N], the abscissas.
    Input, int I, the index of the first-order basis function.
    Indices are 0-based
    Input, double XV, the evaluation point.
    Output, double HERMITE_BASIS_0, the value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ i = s_data->a1;
	ityp * x = s_data->a2;
	ityp xv = s_data->a3;
	
    ityp *factor;
    dim_typ j;
    ityp li;
    ityp lp;
    ityp lpp;
    ityp value;

    if ( i < 0 || n - 1 < i )
    {
    	result = MAX_VAL;
        return &result;
    }

    factor = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    /*
    L(X) = product ( X - X(1:N) )

    L'(X(I)).
    */
    for ( j = 0; j < n; ++j )
        factor[j] = x[i] - x[j];
    factor[i] = 1.00;

    lp = r8vec_product ( n, factor );
    /*
    LI(X) = L(X) / ( X - X(I) ) / L'(X(I))
    */
    for ( j = 0; j < n; ++j )
        factor[j] = xv - x[j];
    factor[i] = 1.00;

    li = r8vec_product ( n, factor ) / lp;
    /*
    L''(X(I)).
    */
    lpp = 0.00;
    for ( j = 0; j < n; ++j )
        factor[j] = x[i] - x[j];

    factor[i] = 1.00;

    for ( j = 0; j < n; ++j )
        if ( j != i )
        {
            factor[j] = 1.00;
            lpp +=  2.00 * r8vec_product ( n, factor );
            factor[j] = x[i] - x[j];
        }

    free ( factor );

	result = ( 1.00 - ( xv - x[i] ) * lpp / lp ) * li * li;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hermite_basis_1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_BASIS_1 evaluates a first-order Hermite interpolation basis function.
  Discussion:
    Given ND points XD, with values YD and derivative values YPD, the
    Hermite interpolant can be written as:
      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
    and H1(I;X) is the I-th first order Hermite interpolation basis function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of abscissas.
    Input, double X[N], the abscissas.
    Input, int I, the index of the first-order basis function.
    Indices are 0-based
    Input, double XV, the evaluation point.
    Output, double VALUE, the value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ i = s_data->a1;
	ityp * x = s_data->a2;
	ityp xv = s_data->a3;
	
    ityp bot;
    ityp *factor;
    dim_typ j;
    ityp top;
    ityp value;

    if ( i < 0 || n - 1 < i )
    {
    	result = MAX_VAL;
        return &result;
    }

    factor = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        factor[j] = xv - x[j];
    factor[i] = 1.00;
    top = r8vec_product ( n, factor );

    for ( j = 0; j < n; ++j )
        factor[j] = x[i] - x[j];
    factor[i] = 1.00;
    bot = r8vec_product ( n, factor );
    free ( factor );
    
    result = ( xv - x[i] ) * ( top / bot ) * ( top / bot );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _hermite_interpolant ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_INTERPOLANT sets up a divided difference table from Hermite data.
  Discussion:
    The polynomial represented by the divided difference table can be
    evaluated by calling HERMITE_INTERPOLANT_VALUE.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 October 2011.
  Author:
    John Burkardt
  Reference:
    Carl deBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.
  Parameters:
    Input, int N, the number of items of data
 ( X(I), Y(I), YP(I) ).
    Input, double X[N], the abscissas.
    These values must be distinct.
    Input, double Y[N], YP[N], the function and derivative values.
    Output, double XD[2*N], YD[2*N], the divided difference table
    for the interpolant value.
    Output, double XDP[2*N-1], YDP[2*N-1], the divided difference
    table for the interpolant derivative.
*/
{
	const dt7pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * yp = s_data->a3;
	ityp * xd = s_data->a4;
	ityp * yd = s_data->a5;
	ityp * xdp = s_data->a6;
	ityp * ydp = s_data->a7;
	
    dim_typ i, j, nd, ndp;
    /*
    Copy the data.
    */
    nd = n<<1;

    for ( i = 0; i < n; ++i )
    {
        xd[0+(i<<1)] = x[i];
        xd[1+(i<<1)] = x[i];
    }
    /*
    Carry out the first step of differencing.
    */
    yd[0] = y[0];
    for ( i = 1; i < n; ++i)
    {
        yd[0+(i<<1)] = ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
        yd[1+(j<<1)] = yp[i];
    }
    /*
    Carry out the remaining steps in the usual way.
    */
    for ( i = 2; i < nd; ++i )
        for ( j = nd - 1; i <= j; --j)
            yd[j] = ( yd[j] - yd[j-1] ) / ( xd[j] - xd[j-i] );
    /*
    Compute the difference table for the derivative.
    */
    dif_deriv ( nd, xd, yd, &ndp, xdp, ydp );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hermite_interpolant_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_INTERPOLANT_RULE: quadrature rule for a Hermite interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of abscissas.
    Input, double A, B, the integration limits.
    Input, double X[N], the abscissas.
    Output, double HERMITE_INTERPOLANT_RULE[2*N], the quadrature
    coefficients, given as pairs for function and derivative values
    at each abscissa.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
 	const register ityp b = s_data->a2;
 	ityp * x = s_data->a3;
	
    ityp a_value;
    ityp b_value;
    ityp *c;
    dim_typ i, j, k;
    dim_typ nd, ndp;
    ityp *w;
    ityp *xd;
    ityp *xdp;
    ityp *y;
    ityp *yd;
    ityp *ydp;
    ityp *yp;

    nd = n<<1;
    ndp = nd - 1;

    c = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    w = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    xd = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    xdp = ( ityp * ) malloc ( ndp * sizeof ( ityp ) );
    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    yd = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    ydp = ( ityp * ) malloc ( ndp * sizeof ( ityp ) );
    yp = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        y[i] = yp[i] = 0.00;

    k = 0;

    for ( i = 0; i < n; ++i )
    {
        y[i] = 1.00;
        hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
        dif_to_r8poly ( nd, xd, yd, c );
        a_value = r8poly_ant_val ( n, c, a );
        b_value = r8poly_ant_val ( n, c, b );
        w[k] = b_value - a_value;
        y[i] = 0.00;
        ++ k;

        yp[i] = 1.00;
        hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
        dif_to_r8poly ( nd, xd, yd, c );
        a_value = r8poly_ant_val ( n, c, a );
        b_value = r8poly_ant_val ( n, c, b );
        w[k] = b_value - a_value;
        yp[i] = 0.00;
        ++ k;
    }

    free ( c );
    free ( xd );
    free ( xdp );
    free ( y );
    free ( yd );
    free ( ydp );
    free ( yp );
    return w;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hermite_interpolant_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_INTERPOLANT_VALUE evaluates the Hermite interpolant polynomial.
  Discussion:
    In fact, this function will evaluate an arbitrary polynomial that is
    represented by a difference table.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 October 2011
  Author:
    John Burkardt
  Reference:
    Carl deBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.
  Parameters:
    Input, int ND, the order of the difference table.
    Input, double XD[ND], YD[ND], the difference table for the
    interpolant value.
    Input, double XDP[ND-1], YDP[ND-1], the difference table for
    the interpolant derivative.
    Input, int NV, the number of evaluation points.
    Input, double XV[NV], the evaluation points.
    Output, double YV[NV], YVP[NV], the value of the interpolant and
    its derivative at the evaluation points.
*/
{
	const dt4pitdt3pit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xd = s_data->a1;
	ityp * yd = s_data->a2;
	ityp * xdp = s_data->a3; 
	ityp * ydp = s_data->a4; 
	const register dim_typ nv = s_data->a5;
	ityp * xv = s_data->a6;
	ityp * yv = s_data->a7;
	ityp * yvp = s_data->a8;
	
    dim_typ i, j, ndp;

    ndp = nd - 1;

    for ( j = 0; j < nv; ++j )
    {
        yv[j] = yd[nd-1];
        for ( i = nd - 2; 0 <= i; --i )
            yv[j] = yd[i] + ( xv[j] - xd[i] ) * yv[j];

        yvp[j] = ydp[ndp-1];
        for ( i = ndp - 2; 0 <= i; --i )
            yvp[j] = ydp[i] + ( xv[j] - xdp[i] ) * yvp[j];
    }
    return NULL;

}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8poly_ant_val ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8POLY_ANT_VAL evaluates the antiderivative of an R8POLY in standard form.
  Discussion:
    The constant term of the antiderivative is taken to be zero.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 February 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the polynomial.
    Input, double POLY_COF[N], the polynomial coefficients.  POLY_COF[0]
    is the constant term, and POLY_COF[N-1] is the coefficient of X**(N-1).
    Input, double XVAL, the point where the antiderivative is to be
    evaluated.
    Output, double R8POLY_ANT_VAL, the value of the antiderivative of the polynomial
    at XVAL.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * poly_cof = s_data->a1;
	const register ityp xval = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = n-1; 0 <= i; --i)
        value = ( value + poly_cof[i] / ( ityp ) ( i + 1 ) ) * xval;
    
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8poly_degree ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8POLY_DEGREE returns the degree of a polynomial.
  Discussion:
    The degree of a polynomial is the index of the highest power
    of X with a nonzero coefficient.
    The degree of a constant polynomial is 0.  The degree of the
    zero polynomial is debatable, but this routine returns the
    degree as 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int NA, the dimension of A.
    Input, double A[NA+1], the coefficients of the polynomials.
    Output, int R8POLY_DEGREE, the degree of A.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpit * const s_data = data;
	const register dim_typ na = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ degree = na;
    while ( 0 < degree )
    {
        if ( a[degree] )
        {
        	result = degree;
        	return &result;
        	
        }
        -- degree;
    }
    
    result = degree;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_chebyshev_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_CHEBYSHEV_NEW creates a vector of Chebyshev spaced values.
  Discussion:
    An R8VEC is a vector of R8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 June 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, double A_FIRST, A_LAST, the first and last entries.
    Output, double R8VEC_CHEBYSHEV_NEW[N], a vector of Chebyshev spaced data.
*/
{
	const _2itdt * const s_data = data;
	
	const register ityp a_first = s_data->a0;
	const register ityp a_last = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    ityp *a;
    ityp c;
    dim_typ i;
    ityp theta;

    a = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
        a[0] = ( a_first + a_last ) / 2.00;
    else
    {
        for ( i = 0; i < n; ++i )
        {
            theta = ( ityp ) ( n - i - 1 ) * M_PI / ( ityp ) ( n - 1 );
            c = cos ( theta );

            if ( ( n % 2 ) == 1 && (i<<1) + 1 == n )
                c = 0.00;
            a[i] = ( ( 1.00 - c ) * a_first+ ( 1.00 + c ) * a_last )/   2.00;

        }
    }

    return a;
}

#endif
