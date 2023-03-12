#ifndef __DISABLEDEEP_TOMS655

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _scmm ( void * data)
/******************************************************************************/
/*
  Purpose:
    SCMM computes moments of a classical weight function scaled to [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int M, the number of moments.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, double A, B, the interval endpoints.
    Output, double W(M), the scaled moments.
*/
{
	const _2dt4it * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ kind = s_data->a1;
	ityp alpha = s_data->a2;
	ityp beta = s_data->a3;
	ityp a = s_data->a4;
	ityp b = s_data->a5;
	
	ityp al;
	ityp be;
	dim_typ i;
	ityp p;
	ityp q;
	ityp temp;
	ityp tmp;
	ityp *w;
	
	temp = r8_epsilon ( );
	
	if ( kind == 1 )
	{
		al = be = 0.00;
		if ( ( b - a ) <= temp )
			return NULL;
		q = ( b - a ) / 2.00;
		p = pow ( q, al + be + 1.00 );
	}
	else if ( kind == 2 )
	{
		al = be = -0.50;
		be = -0.5;
		if ( ( b - a ) <= temp )
			return NULL;
		q = ( b - a ) / 2.00;
		p = pow ( q, al + be + 1.00 );
	}
	else if ( kind == 3 )
	{
		al = be = alpha;
		if ( ( b - a ) <= temp )
			return NULL;
		q = ( b - a ) / 2.00;
		p = pow ( q, al + be + 1.00 );
	}
	else if ( kind == 4 )
	{
		al = alpha;
		be = beta;
		if ( ( b - a ) <= temp )
			return NULL;
		q = ( b - a ) / 2.00;
		p = pow ( q, al + be + 1.00 );
	}
	else if ( kind == 5 )
	{
		if ( b <= 0.0 )
			return NULL;
		q = 1.00 / b;
		p = pow ( q, alpha + 1.00 );
	}
	else if ( kind == 6 )
	{
		if ( b <= 0.00 )
			return NULL;
		q = 1.00 / sqrt ( b );
		p = pow ( q, alpha + 1.00 );
	}
	else if ( kind == 7 )
	{
		al = alpha;
		be = 0.00;
		if ( ( b - a ) <= temp )
			return NULL;
		q = ( b - a ) / 2.00;
		p = pow ( q, al + be + 1.00 );
	}
	else if ( kind == 8 )
	{
		if ( a + b <= 0.00 )
			return NULL;
		q = a + b;
		p = pow ( q, alpha + beta + 1.00 );
	}
	/*
	Compute the moments in W.
	*/
	w = wm ( m, kind, alpha, beta );
	
	tmp = p;
	
	for ( i = 0; i < m; ++i)
	{
		w[i] *= tmp;
		tmp *= q;
	}
	
	return w;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cawiq ( void * data)
/******************************************************************************/
/*
   Purpose:
     CAWIQ computes quadrature weights for a given set of knots.
   Discussion:
     This routine is given a set of distinct knots, T, their multiplicities MLT,
     the Jacobi matrix associated with the polynomials orthogonal with respect
     to the weight function W(X), and the zero-th moment of W(X).
     It computes the weights of the quadrature formula
       sum ( 1 <= J <= NT ) sum ( 0 <= I <= MLT(J) - 1 ) wts(j) d^i/dx^i f(t(j))
     which is to approximate
       integral ( a < x < b ) f(t) w(t) dt
     The routine makes various checks, as indicated below, sets up
    various vectors and, if necessary, calls for the diagonalization
    of the Jacobi matrix that is associated with the polynomials
    orthogonal with respect to W(X) on the interval A, B.
    Then for each knot, the weights of which are required, it calls the
    routine CWIQD which to compute the weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, int MLT[NT], the multiplicity of the knots.
    Input, int NWTS, the number of weights.
    Input/output, int NDX[NT], associates with each distinct
    knot T(J), an integer NDX(J) which is such that the weight to the I-th
    derivative value of F at the J-th knot, is stored in
      WTS(abs(NDX(J))+I) for J = 1,2,...,NT, and I = 0,1,2,...,MLT(J)-1.
    The sign of NDX includes the following information:
    > 0, weights are wanted for this knot
    < 0, weights not wanted for this knot but it is included in the quadrature
    = 0. means ignore this knot completely.

    Input, int KEY, indicates structure of WTS and NDX.
    KEY is an integer with absolute value between 1 and 4.
    The sign of KEY choosed the form of WTS:
    0 < KEY, WTS in standard form.
    0 > KEY, J]WTS(J) required.
    The absolute value has the following effect:
    1, set up pointers in NDX for all knots in T array (routine CAWIQ does
    this).  the contents of NDX are not tested on input and weights are
    packed sequentially in WTS as indicated above.
    2, set up pointers only for knots which have nonzero NDX on input.  All
    knots which have a non-zero flag are allocated space in WTS.
    3, set up pointers only for knots which have NDX > 0 on input.  Space in
    WTS allocated only for knots with NDX > 0.
    4, NDX assumed to be preset as pointer array on input.
    Input, int NST, the dimension of the Jacobi matrix.
    NST should be between (N+1)/2 and N.  The usual choice will be (N+1)/2.
    Input/output, double AJ[NST], BJ[NST].
    If JDF = 0 then AJ contains the  diagonal of the Jacobi matrix and
    BJ(1:NST-1) contains the subdiagonal.
    If JDF = 1, AJ contains the eigenvalues of the Jacobi matrix and
    BJ contains the squares of the elements of the first row of U, the
    orthogonal matrix which diagonalized the Jacobi matrix as U*D*U'.
    Input/output, int JDF, indicates whether the Jacobi
    matrix needs to be diagonalized.
    0, diagonalization required;
    1, diagonalization not required.
    Input, double ZEMU, the zero-th moment of the weight
    function W(X).
    Output, double CAWIQ[NWTS], the weights.
*/
{
	const dtpitpidtpi2dt2pitdtit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	int * mlt = s_data->a2;
	const register dim_typ nwts = s_data->a3;
	int * ndx = s_data->a4;
	const register dim_typ key = s_data->a5;
	const register dim_typ nst = s_data->a6;
	ityp * aj = s_data->a7;
	ityp * bj = s_data->a8;
	register dim_typ jdf = s_data->a9;
	const register ityp zemu = s_data->a10; 
	
    dim_typ i;
    dim_typ ip;
    dim_typ j;
    dim_typ jj;
    dim_typ jp;
    dim_typ k;
    dim_typ l;
    dim_typ m;
    dim_typ mnm;
    dim_typ mtj;
    dim_typ n;
    ityp p;
    ityp prec;
    ityp *r;
    ityp tmp;
    ityp *xk;
    ityp *wtmp;
    ityp *wts;
    ityp *z;

    prec = r8_epsilon ( );

    if ( nt == 0)
        return NULL;
    /*
    Check for indistinct knots.
    */
    if ( 1 < nt )
    {
        k = nt - 1;
        for ( i = 1; i <= k; ++i )
        {
            tmp = t[i-1];
            l = i + 1;
            for ( j = l; j <= nt; ++j )
                if ( abs ( tmp - t[j-1] ) <= prec )
                    return NULL;
        }
    }
    /*
    Check multiplicities,
    Set up various useful parameters and
    set up or check pointers to WTS array.
    */
    l = abs ( key );

    if ( l < 1 || 4 < l )
        return NULL;

    k = 1;

    if ( l == 1 )
    {
        for ( i = 1; i <= nt; ++i )
        {
            ndx[i-1] = k;
            if ( mlt[i-1] < 1 )
                return NULL;
            k += mlt[i-1];
        }
        n = k - 1;
    }
    else if ( l == 2 || l == 3 )
    {
        n = 0;

        for ( i = 1; i <= nt; ++i )
        {
            if ( ndx[i-1] == 0 )
                continue;

            if ( mlt[i-1] < 1 )
                return NULL;

            n += mlt[i-1];

            if ( ndx[i-1] < 0 && l == 3 )
                continue;

            ndx[i-1] = abs ( k ) * i4_sign ( ndx[i-1] );
            k += mlt[i-1];
        }

        if ( nwts + 1 < k )
            return NULL;
    }
    else if ( l == 4 )
    {
        for ( i = 1; i <= nt; ++i )
        {
            ip = abs ( ndx[i-1] );

            if ( ip == 0 )
                continue;

            if ( nwts < ip + mlt[i-1] )
                return NULL;

            if ( i == nt )
                break;

            l = i + 1;
            for ( j = l; j <= nt; ++j )
            {
                jp = abs ( ndx[j-1] );
                if ( jp != 0 && jp <= ip + mlt[i-1] && ip <= jp + mlt[j-1] )
                    break;
            }
        }
    }
    /*
    Test some parameters.
    */

    if ( nst < ( n + 1 ) / 2 || zemu <= 0.0 )
        return NULL;

    wts = ( ityp * ) malloc ( nwts * sizeof ( ityp ) );
    /*
    Treat a quadrature formula with 1 simple knot first.
    */
    if ( n <= 1 )
        for ( i = 0; i < nt; ++i )
            if ( 0 < ndx[i] )
            {
                wts[ abs ( ndx[i] ) - 1 ] = zemu;
                return wts;
            }
    /*
    Carry out diagonalization if not already done.
    */
    if ( jdf == 0 )
    {
        /*
        Set unit vector in work field to get back first row of Q.
        */
        z = ( ityp * ) malloc ( nst * sizeof ( ityp ) );

        for ( i = 0; i < nst; ++i )
            z[i] = 0.00;
        z[0] = 1.00;
        /*
        Diagonalize the Jacobi matrix.
        */
        imtqlx ( nst, aj, bj, z );
        /*
        Signal Jacobi matrix now diagonalized successfully.
        */
        jdf = 1;
        /*
        Save squares of first row of U in subdiagonal array.
        */
        for ( i = 0; i < nst; ++i )
            bj[i] = z[i] * z[i];
        free ( z );
    }
    /*
    Find all the weights for each knot flagged.
    */
    for ( i = 1; i <= nt; ++i )
    {
        if ( ndx[i-1] <= 0 )
            continue;
        m = mlt[i-1];
        mnm = MAX ( n - m, 1 );
        l = MIN ( m, n - m + 1 );
        /*
        Set up K-hat matrix for CWIQD with knots according to their multiplicities.
        */
        xk = ( ityp * ) malloc ( mnm * sizeof ( ityp ) );

        k = 1;
        for ( j = 1; j <= nt; ++j )
                if (ndx[j-1] && j != i )
                    for ( jj = 1; jj <= mlt[j-1]; ++jj, ++k )
                        xk[k-1] = t[j-1];

        /*
        Set up the right principal vector.
        */
        r = ( ityp * ) malloc ( l * sizeof ( ityp ) );

        r[0] = 1.00 / zemu;
        for ( j = 1; j < l; ++j )
            r[j] = 0.00;
        /*
        Pick up pointer for the location of the weights to be output.
        */
        k = ndx[i-1];
        /*
        Find all the weights for this knot.
        */
        wtmp = cwiqd ( m, mnm, l, t[i-1], xk, nst, aj, bj, r );

        free ( r );
        free ( xk );

        for ( j = 0; j < m; ++j )
            wts[k-1+j] = wtmp[j];
        free ( wtmp );

        if ( key < 0 )
            continue;
        /*
        Divide by factorials for weights in standard form.
        */
        tmp = 1.00;
        for ( j = 1; j < m - 1; ++j )
        {
            p = j;
            tmp *= p;
            wts[k-1+j] = wts[k-1+j] / tmp;
        }
    }
    return wts;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   cegqfs ( const register dim_typ nt, const register dim_typ kind, const register ityp alpha, const register ityp beta,ityp f ( ityp x, dim_typ i ) )
/******************************************************************************/
/*
  Purpose:
    CEGQFS estimates an integral using a standard quadrature formula.
  Discussion:
    The user chooses one of the standard quadrature rules
    with the default values of A and B.  This routine determines
    the corresponding weights and evaluates the quadrature formula
    on a given function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, double F ( double X, int I ), the name of a routine which
    evaluates the function and some of its derivatives.  The routine
    must return in F the value of the I-th derivative of the function
    at X.  The value  I will always be 0.  The value X will always be a knot.
    Output, double CEGQFS, the value of the quadrature formula
    applied to F.
*/
{
    dim_typ lu = 0;
    ityp qfsum;
    ityp *t;
    ityp *wts;

    t = ( ityp * ) malloc ( nt * sizeof ( ityp ) );
    wts = ( ityp * ) malloc ( nt * sizeof ( ityp ) );

    cgqfs ( nt, kind, alpha, beta, lu, t, wts );
    /*
    Evaluate the quadrature sum.
    */
    qfsum = eiqfs ( nt, t, wts, f );

    free ( t );
    free ( wts );

    return qfsum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   ceiqf ( const register dim_typ nt, ityp t[static nt], int mlt[static nt], const register dim_typ kind, const register ityp alpha,ityp beta, ityp a, ityp b, ityp f ( ityp x, dim_typ i ) )
/******************************************************************************/
/*
  Purpose:
    CEIQF constructs and applies a quadrature formula based on user knots.
  Discussion:
    The knots may have multiplicity.  The quadrature interval is over
    any valid A, B.  A classical weight function is selected by the user.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, int MLT[NT], the multiplicity of the knots.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, double A, B, the interval endpoints.
    Input, double F ( double X, int I ), the name of a routine which
    evaluates the function and some of its derivatives.  The routine
    must return in F the value of the I-th derivative of the function
    at X.  The highest value of I will be the maximum value in MLT minus
    one.  The value X will always be a knot.
    Output, double CEIQF, the value of the quadrature formula
    applied to F.
*/
{
    dim_typ i;
    dim_typ key;
    dim_typ lu;
    dim_typ m;
    dim_typ n;
    int *ndx;
    ityp qfsum;
    ityp *wts;

    lu = n = 0;
    for ( i = 0; i < nt; ++i )
        n += mlt[i];

    key = 1;
    ndx = ( int  * ) malloc ( nt * sizeof ( int ) );

    wts = ciqf ( nt, t, mlt, n, ndx, key, kind, alpha, beta, a, b, lu );

    qfsum = eiqf ( nt, t, mlt, wts, n, ndx, key, f );

    free ( ndx );
    free ( wts );

    return qfsum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   ceiqfs ( const register dim_typ nt, ityp t[], int mlt[], const register dim_typ kind, const register ityp alpha,const register dim_typ beta, ityp f ( ityp x, dim_typ i ) )
/******************************************************************************/
/*
  Purpose:
    CEIQFS computes and applies a quadrature formula based on user knots.
  Discussion:
    The knots may have multiplicity.  The quadrature interval is over
    the standard interval A, B for the classical weight function selected
    by the user.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, int MLT[NT], the multiplicity of the knots.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, double F ( double X, int I ), the name of a routine which
    evaluates the function and some of its derivatives.  The routine
    must return in F the value of the I-th derivative of the function
    at X.  The highest value of I will be the maximum value in MLT minus
    one.  The value X will always be a knot.
    Output, double CEIQFS, the value of the quadrature formula
    applied to F.
*/
{
    dim_typ i;
    dim_typ key;
    dim_typ lu;
    dim_typ n;
    int *ndx;
    ityp qfsum;
    ityp *wts;

    lu = n = 0;
    for ( i = 0; i < nt; ++i)
        n += mlt[i];
    ndx = ( int * ) malloc ( nt * sizeof ( int ) );
    key = 1;

    wts = ciqfs ( nt, t, mlt, n, ndx, key, kind, alpha, beta, lu );

    qfsum = eiqf ( nt, t, mlt, wts, n, ndx, key, f );

    free ( ndx );
    free ( wts );

    return qfsum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cgqfs ( void * data)
/******************************************************************************/
/*
  Purpose:
    CGQFS computes knots and weights of a Gauss quadrature formula.
  Discussion:
    This routine computes the knots and weights of a Gauss quadrature
    formula with:
    * a classical weight function with default values for A and B;
    * only simple knots
    * optionally print knots and weights and a check of the moments
    Use routine EIQFS to evaluate a quadrature formula computed by
    this routine.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, int LO, selects the action.
    > 0, compute and print knots and weights.  Print moments check.
    = 0, compute knots and weights.
    < 0, compute and print knots and weights.
    Output, double T[NT], the knots.
    Output, double WTS[NT], the weights.
*/
{
	const _2dt2itdt2pit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	const register dim_typ kind = s_data->a1;
	ityp alpha = s_data->a2;
	ityp beta = s_data->a3;
	const register dim_typ lo = s_data->a4;
	ityp * t = s_data->a5;
	ityp * wts = s_data->a6;
	
    dim_typ i;
    dim_typ key;
    dim_typ m;
    dim_typ mex;
    int *mlt;
    dim_typ mmex;
    dim_typ mop;
    int *ndx;
    ityp *w;
    /*
    Check there is enough workfield and assign workfield
    */
    key = 1;
    mop = nt<<1;
    m = mop + 1;
    mex = m + 2;
    mmex = MAX ( mex, 1 );
    /*
    Compute the Gauss quadrature formula for default values of A and B.
    */
    cdgqf ( nt, kind, alpha, beta, t, wts );
    /*
    Exit if no print required.
    */
    if ( lo != 0 )
    {
        mlt = ( int * ) malloc ( nt * sizeof ( int ) );
        for ( i = 0; i < nt; ++i)
            mlt[i] = 1;
        ndx = ( int * ) malloc ( nt * sizeof ( int ) );
        for ( i = 0; i < nt; ++i)
            ndx[i] = i + 1;
        w = ( ityp * ) malloc ( mmex * sizeof ( ityp ) );

        chkqfs ( t, wts, mlt, nt, nt, ndx, key, w, mop, mmex, kind, alpha,
        beta, lo );

        free ( mlt );
        free ( ndx );
        free ( w );
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _chkqf ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHKQF computes and prints the moments of a quadrature formula.
  Discussion:
    The quadrature formula is based on a clasical weight function with
    any valid A, B.
    No check can be made for non-classical weight functions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, double T[NT], the knots.
    Input, double WTS[NWTS], the weights.
    Input, int MLT[NT], the multiplicity of the knots.
    Input, int NT, the number of knots.
    Input, int NWTS, the number of weights.
    Input, int NDX[NT], used to index the array WTS.
    If KEY = 1, then NDX need not be preset.  For more details see the
    comments in CAWIQ.
    Input, int KEY, indicates the structure of the WTS
    array.  It will normally be set to 1.  This will cause the weights to be
    packed sequentially in array WTS.  For more details see the comments
    in CAWIQ.
    Input, int MOP, the expected order of precision of the
    quadrature formula.
    Input, int MEX, the number of moments required to be
    tested.  Set MEX = 1 and LO < 0 for no moments check.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, int LO, selects the action to carry out.
     > 0, print weights and moment tests.
     = 0, print nothing. compute moment test.
     < 0, print weights only. don't compute moment tests.
    Input, double A, B, the interval endpoints.
*/
{
	const _2pitpi2dtpi4dt2itdt2it * const s_data = data;
	ityp * t = s_data->a0;
	ityp * wts = s_data->a1;
	int * mlt = s_data->a2;
	dim_typ nt = s_data->a3;
	dim_typ nwts = s_data->a4;
	int * ndx = s_data->a5;
	dim_typ key = s_data->a6; 
	dim_typ mop = s_data->a7; 
	dim_typ mex = s_data->a8; 
	dim_typ kind = s_data->a9; 
	ityp alpha = s_data->a10;
	ityp beta = s_data->a11;
	dim_typ lo = s_data->a12; 
	ityp a = s_data->a13;
	ityp b = s_data->a14;
	
    dim_typ i;
    dim_typ izero;
    short neg;
    ityp *t2;
    ityp tmp;
    ityp *w = ( ityp * ) malloc ( mex * sizeof ( ityp ) );

    if ( lo != 0 )
    {
        izero = 0;
        chkqfs ( t, wts, mlt, nt, nwts, ndx, key, w, mop, mex, izero,
        alpha, beta, - abs ( lo ) );
    }

    if ( 0 <= lo )
    {
        /*
        Compute the moments in W.
        */
        w = scmm ( mex, kind, alpha, beta, a, b );
        tmp = ( kind == 1 || kind == 2 || kind == 3 || kind == 4 || kind == 7 ) ? ( b + a ) / 2.00 : a;

        t2 = (ityp * ) malloc ( nt * sizeof ( ityp ) );

        for ( i = 0; i < nt; ++i )
            t2[i] = t[i] - tmp;

        neg = -1;
        /*
        Check moments.
        */
        chkqfs ( t2, wts, mlt, nt, nwts, ndx, key, w, mop, mex, neg, alpha, beta,lo );

        free ( t2 );
    }

    free ( w );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _chkqfs ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHKQFS checks the polynomial accuracy of a quadrature formula.
  Discussion:
    This routine will optionally print weights, and results of a moments test.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, double T[NT], the knots.
    Input, double WTS[NWTS], the weights.
    Input, int MLT[NT], the multiplicity of the knots.
    Input, int NT, the number of knots.
    Input, int NWTS, the number of weights.
    Input, int NDX[NT], used to index the array WTS.
    If KEY = 1, then NDX need not be preset.  For more details see the
    comments in CAWIQ.
    Input, int KEY, indicates the structure of the WTS
    array.  It will normally be set to 1.  This will cause the weights to be
    packed sequentially in array WTS.  For more details see the comments
    in CAWIQ.
    Input/output, double W[MEX], the moments array.
    This is input only if KIND = 0.
    Input, int MOP, the expected order of precision of the
    quadrature formula.
    Input, int MEX, the number of moments to be tested.
    MEX must be at least 1.  Set MEX = 1 and LO < 0 for no moment check.
    Input, int KIND, the rule.
    0, unknown weight function (the user must set the first MEX moments in
       array W in this case.)
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, int LO, selects the action to carry out.
     > 0, print weights and moment tests.
     = 0, print nothing.  Dompute moment test.
     < 0, print weights only.  Do not compute moment tests.
  Local Parameters:
    Local, double E[MEX], ER[MEX], the absolute and relative
    errors of the quadrature formula applied to (X-DEL)^n.
    Local, double QM[MEX], the values of the quadrature formula
    applied to (X-DEL)^N.
*/
{
	const _2pitpi2dtpidtpit3dt2itdt * const s_data = data;
	ityp * t = s_data->a0;
	ityp * wts = s_data->a1;
	int * mlt = s_data->a2;
	dim_typ nt = s_data->a3;
	dim_typ nwts = s_data->a4;
	int * ndx = s_data->a5;
	dim_typ key = s_data->a6;
	ityp * w = s_data->a7;
	dim_typ mop = s_data->a8;
	dim_typ mex = s_data->a9;
	dim_typ kind = s_data->a10;
	ityp alpha = s_data->a11;
	ityp beta = s_data->a12;
	dim_typ lo = s_data->a13; 
	
    ityp *e;
    ityp ek;
    ityp emn;
    ityp emx;
    ityp erest;
    ityp ern;
    ityp erx;
    ityp *er;
    dim_typ i;
    dim_typ j;
    dim_typ jl;
    dim_typ k;
    dim_typ kindp;
    dim_typ kjl;
    dim_typ l;
    dim_typ m;
    dim_typ mx;
    ityp px;
    ityp tmp;
    ityp tmpx;
    ityp prec;
    ityp *qm;
    /*
    KIND may be set to -1 to allow printing of moments only.

    This feature is only used internally, by CHKQF.
    */
    kindp = MAX ( 0, kind );

    if ( lo < 0 )
        return NULL;
    /*
    Compute the moments in W.
    */
    if ( kindp != 0 )
        w = wm ( mex, kindp, alpha, beta );

    e = ( ityp * ) malloc ( mex * sizeof ( ityp ) );
    er = ( ityp * ) malloc ( mex * sizeof ( ityp ) );
    qm = ( ityp * ) malloc ( mex * sizeof ( ityp ) );

    for ( j = 0; j < mex; ++j)
        qm[j] = 0.00;
    erest = 0.00;

    for ( k = 1; k <= nt; ++k )
    {
        tmp = 1.00;
        l = abs ( ndx[k-1] );
        if ( l == 0 )
            continue;

        erest += abs ( wts[l-1] );
        for ( j = 1; j <= mex; ++j )
        {
            qm[j-1] = qm[j-1] + tmp * wts[l-1];
            tmpx = tmp;
            px = 1.00;
            for( jl = 2; jl <= MIN ( mlt[k-1], mex - j + 1 ); ++jl )
            {
                kjl = j + jl - 1;
                tmpx *= ( kjl - 1 );
                qm[kjl-1] += tmpx * wts[l+jl-2] / px;
                if ( key <= 0 )
                    px *= jl;
            }
            tmp *= t[k-1];
        }

    }
    for ( j = 0; j < mex; j++ )
    {
        e[j] = w[j] - qm[j];
        er[j] = e[j] / ( abs ( w[j] ) + 1.00 );
    }
    /*
    For some strange weight functions W(1) may vanish.
    */
    erest /= ( abs ( w[0] ) + 1.00 );

    if ( 0 < lo )
    {
        m = mop + 1;
        mx = MIN ( mop, mex );

        emx = abs ( e[0] );
        emn = emx;
        erx = abs ( er[0] );
        ern = erx;
        for ( k = 1; k < mx; ++k )
        {
            emx = MAX ( abs ( e[k] ), emx );
            emn = MIN ( abs ( e[k] ), emn );
            erx = MAX ( abs ( er[k] ), erx );
            ern = MIN ( abs ( er[k] ), ern );
        }

        if ( m <= mex )
        {
            ek = e[m-1];
            for ( j = 1; j <= mop; ++j )
                k /= ( ityp ) ( j );

        }
    }

    free ( e );
    free ( er );
    free ( qm );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ciqf ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIQF computes weights for a classical weight function and any interval.
  Discussion:
    This routine compute somes or all the weights of a quadrature formula
    for a classical weight function with any valid A, B and a given set of
    knots and multiplicities.
    The weights may be packed into the output array WTS according to a
    user-defined pattern or sequentially.
    The routine will also optionally print knots and weights and a check
    of the moments.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, int MLT[NT], the multiplicity of the knots.
    Input, int NWTS, the number of weights.
    Input/output, int NDX[NT], used to index the output
    array WTS.  If KEY = 1, then NDX need not be preset.  For more
    details see the comments in CAWIQ.
    Input, int KEY, indicates the structure of the WTS
    array.  It will normally be set to 1.  This will cause the weights to be
    packed sequentially in array WTS.  For more details see the comments
    in CAWIQ.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, double A, B, the interval endpoints.
    Input, int LO, selects the actions to perform.
     > 0, compute and print weights.  Print moments check.
     = 0, compute weights.
     < 0, compute and print weights.
    Output, double CIQF[NWTS], the weights.
*/
{
	const dtpitpidtpi2dt4itdt * const s_data = data; 
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	int * mlt = s_data->a2;
	dim_typ nwts = s_data->a3;
	int * ndx = s_data->a4;
	dim_typ key = s_data->a5;
	dim_typ kind = s_data->a6;
	ityp alpha = s_data->a7;
	ityp beta = s_data->a8;
	ityp a = s_data->a9;
	ityp b = s_data->a10;
	dim_typ lo = s_data->a11;
	
    dim_typ j;
    dim_typ k;
    dim_typ l;
    dim_typ lu;
    dim_typ m;
    dim_typ mex;
    dim_typ mop;
    ityp *st;
    ityp *wts;

    m = 1;
    l = abs ( key );

    for ( j = 0; j < nt; ++j )
    if ( l == 1 || abs ( ndx[j] ) != 0 )
        m += mlt[j];

    if ( nwts + 1 < m )
        return NULL;

    mex = 2 + m;
    /*
    Scale the knots to default A, B.
    */
    st = sct ( nt, t, kind, a, b );
    /*
    Compute the weights.
    */
    lu = 0;

    wts = ciqfs ( nt, st, mlt, nwts, ndx, key, kind, alpha, beta, lu );
    /*
    Don't scale user's knots - only scale weights.
    */
    scqf ( nt, st, mlt, nwts, wts, ndx, wts, st, kind, alpha, beta, a, b );

    if ( lo != 0 )
    {
        mop = m - 1;
        chkqf ( t, wts, mlt, nt, nwts, ndx, key, mop, mex, kind,
        alpha, beta, lo, a, b );
    }

    return wts;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _ciqfs ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIQFS computes some weights of a quadrature formula in the default interval.
  Discussion:
    This routine computes some or all the weights of a quadrature formula
    for a classical weight function with default values of A and B,
    and a given set of knots and multiplicities.
    The weights may be packed into the output array WTS according to a
    user-defined pattern or sequentially.
    The routine will also optionally print knots and weights and a check of
    the moments.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, int MLT[NT], the multiplicity of the knots.
    Input, int NWTS, the number of weights.
    Input/output, int NDX[NT],  used to index the output
    array WTS.  If KEY = 1, then NDX need not be preset.  For more
    details see the comments in CAWIQ.
    Input, int KEY, indicates the structure of the WTS
    array.  It will normally be set to 1.  For more details see
    the comments in CAWIQ.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, int LO, selects the actions to perform.
     > 0, compute and print weights.  Print moments check.
     = 0, compute weights.
     < 0, compute and print weights.
    Output, double CIQFS[NWTS], the weights.
*/
{
	const dtpitpidtpi2dt2itdt * const s_data = data; 
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	int * mlt = s_data->a2;
	dim_typ nwts = s_data->a3;
	int * ndx = s_data->a4;
	dim_typ key = s_data->a5;
	dim_typ kind = s_data->a6;
	ityp alpha = s_data->a7;
	ityp beta = s_data->a8;
	dim_typ lo = s_data->a9;
	
    ityp *aj;
    ityp *bj;
    dim_typ j;
    dim_typ jdf;
    dim_typ l;
    dim_typ m;
    dim_typ mex;
    dim_typ mmex;
    dim_typ mop;
    dim_typ n;
    dim_typ nst;
    ityp *w;
    ityp *wts;
    ityp zemu;

    jdf = n = 0;
    l = abs ( key );

    for ( j = 0; j < nt; ++j)
        if ( l == 1 || abs ( ndx[j] ) != 0 )
            n += mlt[j];
    /*
    N knots when counted according to multiplicity.
    */
    if ( nwts < n )
        return NULL;

    m = n + 1;
    mex = 2 + m;
    nst = m / 2;
    /*
    Get the Jacobi matrix.
    */
    aj = ( ityp * ) malloc ( nst * sizeof ( ityp ) );
    bj = ( ityp * ) malloc ( nst * sizeof ( ityp ) );

    zemu = class_matrix ( kind, nst, alpha, beta, aj, bj );
    /*
    Call weights routine.
    */
    wts = cawiq ( nt, t, mlt, n, ndx, key, nst, aj, bj, jdf, zemu );

    free ( aj );
    free ( bj );
    /*
    Call checking routine.
    */
    if ( lo != 0 )
    {
        mop = m - 1;
        w = ( ityp * ) malloc ( mex * sizeof ( ityp ) );
        chkqfs ( t, wts, mlt, nt, n, ndx, key, w, mop, mex, kind,alpha, beta, lo );
        free ( w );
    }
    return wts;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cliqf ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLIQF computes a classical quadrature formula, with optional printing.
  Discussion:
    This routine computes all the weights of an interpolatory
    quadrature formula with
    1. only simple knots and
    2. a classical weight function with any valid A and B, and
    3. optionally prints the knots and weights and a check of the moments.
    To evaluate this quadrature formula for a given function F,
    call routine EIQFS.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, double A, B, the interval endpoints.
    Input, int LO, indicates what is to be done.
    > 0, compute and print weights and moments check.
    = 0, compute weights.
    < 0, compute and print weights.
    Output, double CLIQF[NT], the weights.
*/
{
	const dtpitdt4itdt * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	dim_typ kind = s_data->a2;
	ityp alpha = s_data->a3;
	ityp beta = s_data->a4;
	ityp a = s_data->a5;
	ityp b = s_data->a6;
	dim_typ lo = s_data->a7;
	
    dim_typ i;
    dim_typ key;
    int *mlt;
    int *ndx;
    ityp *wts;

    key = 1;
    mlt = ( int * ) malloc ( nt * sizeof ( int ) );
    for ( i = 0; i < nt; ++i )
        mlt[i] = 1;
    ndx = ( int * ) malloc ( nt * sizeof ( int ) );

    wts = ciqf ( nt, t, mlt, nt, ndx, key, kind, alpha, beta, a, b, lo );

    free ( mlt );
    free ( ndx );

    return wts;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cliqfs ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLIQFS computes the weights of a quadrature formula in the default interval.
  Discussion:
    This routine computes the weights of an interpolatory quadrature formula
    with a classical weight function, in the default interval A, B,
    using only simple knots.
    It can optionally print knots and weights and a check of the moments.
    To evaluate a quadrature computed by CLIQFS, call EIQFS.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Input, int LO, chooses the printing option.
     > 0, compute weights, print them, print the moment check results.
     0, compute weights.
     < 0, compute weights and print them.
    Output, double CLIQFS[NT], the weights.
*/
{
	const dtpitdt4itdt * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	dim_typ kind = s_data->a2;
	ityp alpha = s_data->a3;
	ityp beta = s_data->a4;
	dim_typ lo = s_data->a5;
	
    dim_typ i;
    dim_typ key;
    int *mlt;
    int *ndx;
    ityp *wts;

    key = 1;
    mlt = ( int * ) malloc ( nt * sizeof ( int ) );

    for ( i = 0; i < nt; ++i )
        mlt[i] = 1;
    ndx = ( int * ) malloc ( nt * sizeof ( int ) );

    wts = ciqfs ( nt, t, mlt, nt, ndx, key, kind, alpha, beta, lo );

    free ( mlt );
    free ( ndx );

    return wts;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cwiqd ( void * data)
/******************************************************************************/
/*
  Purpose:
    CWIQD computes all the weights for a given knot.
  Discussion:
    The variable names correspond to the 1982 reference, and explanations of
    some of the terminology may be found there.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
    Jaroslav Kautsky, Sylvan Elhay,
    Calculation of the Weights of Interpolatory Quadratures,
    Numerische Mathematik,
    Volume 40, 1982, pages 407-422.
  Parameters:
    Input, int M, the multiplicity of the knot in question.
    Input, int NM, is equal to MAX ( N - M, 1 ), where N is
    the number of knots used, counted according to multiplicity.
    Input, int L, MIN ( M, N - M + 1), where N is the number
    of knots used, counted according to multiplicity.
    Input, double V,  the knot in question.
    Input, double XK[NM], all but the last M entries in the
    diagonal of K-hat.
    Input, int NSTAR, the dimension of the Jacobi matrix.
    Input, double PHI[NSTAR], the eigenvalues of the Jacobi matrix.
    Input,double A[NSTAR], the square of the first row of the
    orthogonal matrix that diagonalizes the Jacobi matrix.
    Input, double R[L], used to compute the right
    principal vectors.
    Output, double CWIQD[M], the weights.
*/
{
	const _3dtitpitdt3pit * const s_data = data;
	dim_typ m = s_data->a0;
	dim_typ nm = s_data->a1;
	dim_typ l = s_data->a2;
	ityp v = s_data->a3;
	ityp * xk = s_data->a4;
	dim_typ nstar = s_data->a5; 
	ityp * phi = s_data->a6;
	ityp * a = s_data->a7;
	ityp * r = s_data->a8;
	
    ityp *d;
    dim_typ i, j;
    dim_typ jr;
    dim_typ k;
    dim_typ last;
    dim_typ minil;
    ityp sum;
    ityp tmp;
    ityp *wf;
    ityp *y;
    ityp *z;

    d = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    wf = ( ityp * ) malloc ( nstar * sizeof ( ityp ) );
    y = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    z = ( ityp * ) malloc ( m * sizeof ( ityp ) );
    /*
    Compute products required for Y-hat.
    */
    for ( j = 0; j < nstar; ++j )
    {
        wf[j] = a[j];
        for (i = 0; i < nm; ++i )
            wf[j] *= ( phi[j] - xk[i] );
    }
    /*
    Compute Y-hat.
    */
    for ( i = 0; i < m; ++i)
    {
    sum = 0.00;
        for ( j = 0; j  < nstar; ++j )
        {
            sum += wf[j];
            wf[j] *= ( phi[j] - v );
        }
        y[i] = sum;
    }
    /*
    If N = 1 the right principal vector is already in R.
    Otherwise compute the R-principal vector of grade M-1.
    */
    for ( i = 1; i <= nm; ++i )
    {
        tmp = v - xk[i-1];
        last = MIN ( l, i + 1 );
        for ( jr = 2; jr <= last; ++jr )
        {
            j = last - jr + 2;
            r[j-1] = tmp * r[j-1] + r[j-2];
        }
        r[0] = tmp * r[0];
    }
    /*
    Compute left principal vector(s) and weight for highest derivative.
    The following statement contains the only division in this
    routine.  Any test for overflow should be made after it.
    */
    d[m-1] = y[m-1] / r[0];

    if ( m == 1 )
    {
        free ( wf );
        free ( y );
        free ( z );
        return d;
    }
    /*
    Compute left principal vector.
    */
    z[0] = 1.00 / r[0];
    for ( i = 2; i <= m; ++i )
    {
        sum = 0.00;
        minil = MIN ( i, l );
        for ( j = 2; j <= minil; ++j )
        {
            k = i - j + 1;
            sum += r[j-1] * z[k-1];
        }
        z[i-1] = - sum * z[0];
    }
    /*
    Accumulate weights.
    */
    for ( i = 2; i <= m; ++i )
    {
        sum = 0.00;
        for ( j = 1; j <= i; ++j )
        {
            k = m - i + j;
            sum += z[j-1] * y[k-1];
        }
        k = m - i + 1;
        d[k-1] = sum;
    }

    free ( wf );
    free ( y );
    free ( z );

    return d;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   eiqf ( const register dim_typ nt, ityp t[static nt], int mlt[static nt], ityp wts[], const register dim_typ nwts, int ndx[static nt],dim_typ key, ityp f ( ityp x, dim_typ i ) )
/******************************************************************************/
/*
  Purpose:
    EIQF evaluates an interpolatory quadrature formula.
  Discussion:
   The knots, weights and integrand are supplied.
   All knots with nonzero NDX are used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, int MLT[NT], the multiplicity of the knots.
    Input, double WTS[NWTS], the weights.
    Input, int NWTS, the number of weights.
    Input, int NDX[NT], used to index the array WTS.
    If KEY = 1, then NDX need not be preset.  For more details see the
    comments in CAWIQ.
    Input, int KEY, indicates the structure of the WTS
    array.  It will normally be set to 1.  This will cause the weights to be
    packed sequentially in array WTS.  For more details see the comments
    in CAWIQ.
    Input, double F ( double X, int I ), the name of a routine which
    evaluates the function and some of its derivatives.  The routine
    must return in F the value of the I-th derivative of the function
    at X.  The highest value of I will be the maximum value in MLT minus
    one.  The value X will always be a knot.
    Output, double EIQF, the value of the quadrature formula
    applied to F.
*/
{
    dim_typ i, j, l;
    ityp p;
    ityp qfsum;

    l = abs ( key );

    if ( l < 1 || 4 < l )
        return MAX_VAL;

    qfsum = 0.00;
    for ( j = 0; j < nt; ++j )
    {
        l = abs ( ndx[j] );
        if ( l != 0 )
        {
            p = 1.00;
            for ( i = 0; i < mlt[j]; ++i )
            {
                qfsum += wts[l+i-1] * f ( t[j], i ) / p;
                if ( key <= 0 )
                    p *= ( i + 1 );
            }
        }
    }
    return qfsum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   eiqfs ( const register dim_typ nt, ityp t[static nt], ityp wts[static nt], ityp f ( ityp x, dim_typ i ) )
/******************************************************************************/
/*
  Purpose:
    EIQFS evaluates a quadrature formula defined by CLIQF or CLIQFS.
  Discussion:
    This routine evaluates an interpolatory quadrature formula with all knots
    simple and all knots included in the quadrature.  This routine will be used
    typically after CLIQF or CLIQFS has been called.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the knots.
    Input, double WTS[NT], the weights.
    Input, double F ( double X, int I ), the name of a routine which
    evaluates the function and some of its derivatives.  The routine
    must return in F the value of the I-th derivative of the function
    at X.  The value of I will always be 0.  The value X will always be a knot.
    Output, double EIQFS, the value of the quadrature formula
    applied to F.
*/
{
    ityp qfsum = 0.00;
    for (dim_typ j = 0; j < nt; ++j )
        qfsum += wts[j] * f ( t[j], 0 );
    return qfsum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sct ( void * data)
/******************************************************************************/
/*
  Purpose:
    SCT rescales distinct knots to an interval [A,B].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int NT, the number of knots.
    Input, double T[NT], the original knots.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double A, B, the interval endpoints for which the
    knots ST should be scaled.
    Output, double SCT[NT], the scaled knots.
*/
{
	const dt2itdtpit * const s_data = data;
	
	const register dim_typ kind = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	const register dim_typ nt = s_data->a3;
	ityp * t = s_data->a4;
	
    ityp bma;
    dim_typ i;
    ityp shft;
    ityp slp;
    ityp *st;
    ityp tmp;

    if ( kind <= 0 || 8 < kind )
        return NULL;

    if ( kind == 1 || kind == 2 || kind == 3 || kind == 4 || kind == 7 )
    {
        tmp = r8_epsilon ( );
        bma = b - a;

        if ( bma <= tmp )
            return NULL;
        slp = 2.00 / bma;
        shft = - ( a + b ) / bma;
    }
    else if ( kind == 5 )
    {
        if ( b < 0.0 )
            return NULL;
        slp = b;
        shft = - a * b;
    }
    else if ( kind == 6 )
    {
        if ( b < 0.0 )
            return NULL;
        slp = sqrt ( b );
        shft = - a * slp;
    }
    else if ( kind == 8 )
    {
        slp = 1.00 / ( a + b );

        if ( slp <= 0.00 )
            return NULL;
        shft = - a * slp;
    }

    st = ( ityp * ) malloc ( nt * sizeof ( ityp ) );

    for ( i = 0; i < nt; ++i)
        st[i] = shft + slp * t[i];

    return st;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wm ( void * data)
/******************************************************************************/
/*
  Purpose:
    WM evaluates the first M moments of classical weight functions.
  Discussion:
    W(K) = Integral ( A <= X <= B ) X^(K-1) * W(X) dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int M, the number of moments to evaluate.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Output, double WM[M], the first M moments.
*/
{
	const dt2itdt * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register ityp alpha = s_data->a1;
	const register ityp beta = s_data->a2;
	const register dim_typ kind = s_data->a3;
	
    ityp als;
    dim_typ i;
    dim_typ ja;
    dim_typ jb;
    dim_typ k;
    ityp rk;
    ityp sum;
    ityp tmpa;
    ityp tmpb;
    ityp trm;
    ityp *w;

    w = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    for ( k = 2; k <= m; k += 2 )
        w[k-1] = 0.0;

    if ( kind == 1 )
        for ( k = 1; k <= m; k += 2 )
        {
            rk = ( ityp ) ( k );
            w[k-1] = 2.00 / rk;
        }
    else if ( kind == 2 )
    {
        w[0] = M_PI;
        for ( k = 3; k <= m; k += 2 )
        {
            rk = ( ityp ) ( k );
            w[k-1] = w[k-3] * ( rk - 2.00 ) / ( rk - 1.00 );
        }
    }
    else if ( kind == 3 )
    {
        w[0] = sqrt ( M_PI ) * r8_gamma ( alpha + 1.00 )/ r8_gamma ( alpha + 3.00 / 2.00 );

        for ( k = 3; k <= m; k += 2 )
        {
            rk = ( ityp ) ( k );
            w[k-1] = w[k-3] * ( rk - 2.00 ) / ( 2.00 * alpha + rk );
        }
    }
    else if ( kind == 4 )
    {
        als = alpha + beta + 1.00;
        w[0] = pow ( 2.00, als ) * r8_gamma ( alpha + 1.00 )/ r8_gamma ( als + 1.00 ) * r8_gamma ( beta + 1.00 );

        for ( k = 2; k <= m; k++ )
        {
            sum = 0.00;
            trm = 1.00;
            rk = ( ityp ) ( k );

            for ( i = 0; i <= ( k - 2 ) / 2; ++i )
            {
                tmpa = trm;
                for ( ja = 1; ja <= (i<<1); ++ja )
                    tmpa *= ( alpha + ja ) / ( als + ja );

                for ( jb = 1; jb <= k - (i<<1) - 1; ++jb )
                    tmpa *= ( beta + jb ) / ( als + (i<<1) + jb );

                tmpa /= ( (i<<1) + 1.00 ) *( (i<<1) * ( beta + alpha ) + beta - ( rk - 1.00 ) * alpha ) / ( beta + rk - (i<<1) - 1.00 );
                sum += tmpa;
                trm *= ( rk - (i<<1) - 1.00 )/ ( (i<<1) + 1.00 ) * ( rk - (i<<1) - 2.00 ) / ( (i<<1) + 2.00 );
            }

            if ( ( k % 2 ) != 0 )
            {
                tmpb = 1.00;
                for ( i = 1; i <= k - 1; ++i )
                    tmpb *= ( alpha + i ) / ( als + i );
                sum += tmpb;
            }
            w[k-1] = sum * w[0];
        }
    }
    else if ( kind == 5 )
    {
        w[0] = r8_gamma ( alpha + 1.00 );

        for ( k = 2; k <= m; ++k )
        {
            rk = ( ityp ) ( k );
            w[k-1] = ( alpha + rk - 1.00 ) * w[k-2];
        }
    }
    else if ( kind == 6 )
    {
        w[0] = r8_gamma ( ( alpha + 1.00 ) / 2.00 );

        for ( k = 3; k <= m; k += 2 )
        {
            rk = ( ityp ) ( k );
            w[k-1] = w[k-3] * ( alpha + rk - 2.00 ) / 2.00;
        }
    }
    else if ( kind == 7 )
    {
        als = alpha;
        for ( k = 1; k <= m; k += 2 )
        {
            rk = ( ityp ) ( k );
            w[k-1] = 2.00 / ( rk + als );
        }
    }
    else if ( kind == 8 )
    {
        w[0] = r8_gamma ( alpha + 1.00 )* r8_gamma ( - alpha - beta - 1.00 )/ r8_gamma ( - beta );

        for ( k = 2; k <= m; ++k )
        {
            rk = ( ityp ) ( k );
            w[k-1] = - w[k-2] * ( alpha + rk - 1.00 ) / ( alpha + beta + rk );
        }

    }

    return w;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wtfn ( void * data)
/******************************************************************************/
/*
  Purpose:
    WTFN evaluates the classical weight functions at given points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, double T[NT], the points where the weight function
    is to be evaluated.
    Input, int NT, the number of evaluation points.
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Output, double WTFN[NT], the value of the weight function.
*/
{
	const dt2itdtpit * const s_data = data;
	
	const register dim_typ nt = s_data->a0;
	const register ityp alpha = s_data->a1;
	const register ityp beta = s_data->a2;
	const register dim_typ kind = s_data->a3;
	ityp * t = s_data->a4;
	
    dim_typ i;
    ityp *w;

    w = ( ityp * ) malloc ( nt * sizeof ( ityp ) );

    if ( kind == 1 )
        for ( i = 0; i < nt; ++i)
            w[i] = 1.00;
    else if ( kind == 2 )
        for ( i = 0; i < nt; ++i)
            w[i] = 1.00 / sqrt ( ( 1.00 - t[i] ) * ( 1.00 + t[i] ) );
    else if ( kind == 3 )
    {
        if ( alpha == 0.00 )
            for ( i = 0; i < nt; ++i )
                w[i] = 1.00;
        else
            for ( i = 0; i < nt; ++i)
                w[i] = pow ( ( 1.00 - t[i] ) * ( 1.0 + t[i] ), alpha );
    }
    else if ( kind == 4 )
    {
        if ( alpha == 0.00 )
            for ( i = 0; i < nt; ++i )
                w[i] = 1.00;
        else
            for ( i = 0; i < nt; ++i)
                w[i] = pow ( 1.00 - t[i], alpha );
        if ( beta != 0.00 )
            for ( i = 0; i < nt; ++i )
                w[i] *= pow ( 1.00 + t[i], beta );
    }
    else if ( kind == 5 )
    {
        if ( alpha == 0.00 )
            for ( i = 0; i < nt; ++i )
                w[i] = exp ( - t[i] );
        else
            for ( i = 0; i < nt; ++i)
                w[i] = exp ( - t[i] ) * pow ( t[i], alpha );
    }
    else if ( kind == 6 )
    {
        if ( alpha == 0.00 )
            for ( i = 0; i < nt; ++i )
                w[i] = exp ( - t[i] * t[i] );
        else
            for ( i = 0; i < nt; ++i )
                w[i] = exp ( - t[i] * t[i] ) * pow ( abs ( t[i] ), alpha );
    }
    else if ( kind == 7 )
    {
        if ( alpha != 0.00 )
            for ( i = 0; i < nt; ++i )
                w[i] = pow ( abs ( t[i] ), alpha );
        else
            for ( i = 0; i < nt; ++i )
                w[i] = 1.0;
    }
    else if ( kind == 8 )
    {
        if ( alpha == 0.0 )
            for ( i = 0; i < nt; ++i )
                w[i] = 1.00;
        else
            for ( i = 0; i < nt; ++i)
                w[i] = pow ( t[i], alpha );
        if ( beta != 0.00 )
            for ( i = 0; i < nt; ++i )
                w[i] *= pow ( 1.00 + t[i], beta );
    }
    return w;
}

#endif
