#ifndef __DISABLEDEEP_HERMITERULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cdgqf ( void * data)
/******************************************************************************/
/*
  Purpose:
    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
  Discussion:
    This routine computes all the knots and weights of a Gauss quadrature
    formula with a classical weight function with default values for A and B,
    and only simple knots.
    There are no moments checks and no printing is done.
    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
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
    Output, double T[NT], the knots.
    Output, double WTS[NT], the weights.
*/
{
	const _2dt2it2pit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	const register dim_typ kind = s_data->a1;
	const register ityp alpha = s_data->a2;
	const register ityp beta = s_data->a3;
	ityp * t = s_data->a4;
	ityp * wts = s_data->a5;
	
    ityp *aj;
    ityp *bj;
    ityp zemu;
    
    /*
    Get the Jacobi matrix and zero-th moment.
    */
    aj = ( ityp * ) malloc ( nt * sizeof ( ityp ) );
    bj = ( ityp * ) malloc ( nt * sizeof ( ityp ) );

    zemu = class_matrix ( kind, nt, alpha, beta, aj, bj );
    /*
    Compute the knots and weights.
    */
    sgqf ( nt, aj, bj, zemu, t, wts );

    free ( aj );
    free ( bj );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cgqf ( void * data)
/******************************************************************************/
/*
  Purpose:
    CGQF computes knots and weights of a Gauss quadrature formula.
  Discussion:
    The user may specify the interval (A,B).
    Only simple knots are produced.
    The user may request that the routine print the knots and weights,
    and perform a moment check.
    Use routine EIQFS to evaluate this quadrature formula.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 September 2013
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
    Input, double A, B, the interval endpoints.
    Output, double T[NT], the knots.
    Output, double WTS[NT], the weights.
*/
{
	const _2dt4it2pit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	const register dim_typ kind = s_data->a1;
	ityp alpha = s_data->a2;
	ityp beta = s_data->a3;
	ityp a = s_data->a4;
	ityp b = s_data->a5;
	ityp * t = s_data->a6;
	ityp * wts = s_data->a7;
	
    dim_typ i;
    int *mlt;
    int *ndx;
    /*
    Compute the Gauss quadrature formula for default values of A and B.
    */
    cdgqf ( nt, kind, alpha, beta, t, wts );
    /*
    Prepare to scale the quadrature formula to other weight function with
    valid A and B.
    */
    mlt = ( int * ) malloc ( nt * sizeof ( int ) );
    for ( i = 0; i < nt; ++i )
        mlt[i] = 1;
    ndx = ( int * ) malloc ( nt * sizeof ( int ) );
    for ( i = 0; i < nt; ++i )
        ndx[i] = i + 1;
    scqf ( nt, t, mlt, nt, wts, ndx, wts, t, kind, alpha, beta, a, b );

    free ( mlt );
    free ( ndx );

    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _class_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
  Discussion:
    This routine computes the diagonal AJ and sub-diagonal BJ
    elements of the order M tridiagonal symmetric Jacobi matrix
    associated with the polynomials orthogonal with respect to
    the weight function specified by KIND.
    For weight functions 1-7, M elements are defined in BJ even
    though only M-1 are needed.  For weight function 8, BJ(M) is
    set to zero.
    The zero-th moment of the weight function is returned in ZEMU.
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
    Input, int KIND, the rule.
    1, Legendre,          (a,b)       1.0
    2, Chebyshev Type 1,  (a,b) ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,        (a,b) ((b-x)*(x-a))^alpha
    4, Jacobi,            (a,b)    (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)  (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite, (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,       (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,          (a,inf)  (x-a)^alpha*(x+b)^beta
    Input, int M, the order of the Jacobi matrix.
    Input, double ALPHA, the value of Alpha, if needed.
    Input, double BETA, the value of Beta, if needed.
    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
    of the Jacobi matrix.
    Output, double CLASS_MATRIX, the zero-th moment.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dt2it2pit * const s_data = data;
	const register dim_typ kind = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register ityp alpha = s_data->a2;
	const register ityp beta = s_data->a3;
	ityp * aj = s_data->a4;
	ityp * bj = s_data->a5;
	
    ityp a2b2;
    ityp ab;
    ityp aba;
    ityp abi;
    ityp abj;
    ityp abti;
    ityp apone;
    int i;
    ityp temp;
    ityp temp2;
    ityp zemu;

    temp = r8_epsilon ( );
    temp2 = 0.50;

    if ( 500.0 * temp < abs ( pow ( r8_gamma ( temp2 ), 2 ) - M_PI ) )
    {
    	result = MAX_VAL;
        return &result;
    }

    switch(kind)
    {

        case 1:
        ab = 0.00;

        zemu = 2.00 / ( ab + 1.00 );

        for ( i = 0; i < m; ++i)
            aj[i] = 0.00;

        for ( i = 1; i <= m; ++i )
        {
            abi = i + ab * ( i % 2 );
            abj = (i<<1) + ab;
            bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.00 ) );
        }
        break;
        case 2:
            zemu = M_PI;

            for ( i = 0; i < m; ++i )
                aj[i] = 0.00;


            bj[0] =  sqrt ( 0.50 );
            for ( i = 1; i < m; ++i )
                bj[i] = 0.50;
                break;
        case 3:
            ab = alpha * 2.00;
            zemu = pow ( 2.00, ab + 1.00 ) * pow ( r8_gamma ( alpha + 1.00 ), 2 )/ r8_gamma ( ab + 2.00 );

            for ( i = 0; i < m; ++i )
                aj[i] = 0.00;


            bj[0] = sqrt ( 1.00 / ( 2.00 * alpha + 3.00 ) );
            for ( i = 2; i <= m; ++i )
                bj[i-1] = sqrt ( i * ( i + ab ) / ( 4.00 * pow ( i + alpha, 2 ) - 1.00 ) );
            break;
        case 4:
            ab = alpha + beta;
            abi = 2.00 + ab;
            zemu = pow ( 2.00, ab + 1.00 ) * r8_gamma ( alpha + 1.00 ) * r8_gamma ( beta + 1.00 ) / r8_gamma ( abi );
            aj[0] = ( beta - alpha ) / abi;
            bj[0] = sqrt ( 4.00 * ( 1.00 + alpha ) * ( 1.00 + beta ) / ( ( abi + 1.00 ) * abi * abi ) );
            a2b2 = beta * beta - alpha * alpha;

            for ( i = 2; i <= m; ++i )
            {
                abi = 2.00 * i + ab;
                aj[i-1] = a2b2 / ( ( abi - 2.00 ) * abi );
                abi = abi * abi;
                bj[i-1] = sqrt ( 4.00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) / ( ( abi - 1.00 ) * abi ) );
            }
            break;
        case 5:
            zemu = r8_gamma ( alpha + 1.0 );

            for ( i = 1; i <= m; ++i )
            {
                aj[i-1] = 2.00 * i - 1.00 + alpha;
                bj[i-1] = sqrt ( i * ( i + alpha ) );
            }
            break;
        case 6:
            zemu = r8_gamma ( ( alpha + 1.00 ) / 2.00 );

            for ( i = 0; i < m; ++i )
                aj[i] = 0.00;
            for ( i = 1; i <= m; ++i )
                bj[i-1] = sqrt ( ( i + alpha * ( i % 2 ) ) / 2.00 );
            break;
        case 7:
            ab = alpha;
            zemu = 2.00 / ( ab + 1.00 );

            for ( i = 0; i < m; ++i )
                aj[i] = 0.00;

            for ( i = 1; i <= m; ++i)
            {
                abi = i + ab * ( i % 2 );
                abj = (i<<1) + ab;
                bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.00 ) );
            }
            break;
        case 8:

            ab = alpha + beta;
            zemu = r8_gamma ( alpha + 1.00 ) * r8_gamma ( - ( ab + 1.00 ) ) / r8_gamma ( - beta );
            apone = alpha + 1.00;
            aba = ab * apone;
            aj[0] = - apone / ( ab + 2.00 );
            bj[0] = - aj[0] * ( beta + 1.00 ) / ( ab + 2.00 ) / ( ab + 3.00 );
            for ( i = 2; i <= m; ++i )
            {
                abti = ab + 2.00 * i;
                aj[i-1] = aba + 2.00 * ( ab + i ) * ( i - 1 );
                aj[i-1] = - aj[i-1] / abti / ( abti - 2.00 );
            }

            for ( i = 2; i <= m - 1; ++i )
            {
                abti = ab + 2.00 * i;
                bj[i-1] = i * ( alpha + i ) / ( abti - 1.00 ) * ( beta + i ) / ( abti * abti ) * ( ab + i ) / ( abti + 1.00 );
            }
            bj[m-1] = 0.00;
            for ( i = 0; i < m; ++i )
                bj[i] =  sqrt ( bj[i] );
            break;
        default:
        	{
        		result = MAX_VAL;
            	return &result;
            }
    }

	result = zemu;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _scqf ( void * data)
/******************************************************************************/
/*
  Purpose:
    SCQF scales a quadrature formula to a nonstandard interval.
  Discussion:
    The arrays WTS and SWTS may coincide.
    The arrays T and ST may coincide.
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
    Input, int MLT[NT], the multiplicity of the knots.
    Input, double WTS[NWTS], the weights.
    Input, int NWTS, the number of weights.
    Input, int NDX[NT], used to index the array WTS.
    For more details see the comments in CAWIQ.
    Output, double SWTS[NWTS], the scaled weights.
    Output, double ST[NT], the scaled knots.
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
*/
{
	const dtpitpiipitpi2pitdt4it * const s_data = data;
	dim_typ nt = s_data->a0;
	ityp * t = s_data->a1;
	int * mlt = s_data->a2;
	int nwts = s_data->a3;
	ityp * wts = s_data->a4;
	int * ndx = s_data->a5;
	ityp * swts = s_data->a6;
	ityp * st = s_data->a7;
	dim_typ kind = s_data->a8;
	ityp alpha = s_data->a9;
	ityp beta = s_data->a10;
	ityp a = s_data->a11;
	ityp b = s_data->a12;
	
    ityp al;
    ityp be;
    dim_typ i, k, l;
    ityp p;
    ityp shft;
    ityp slp;
    ityp temp;
    ityp tmp;

    temp = r8_epsilon ( );

    switch(kind)
    {
        case 1:
            al = be = 0.00;
            if ( ( b - a ) <= temp )
                return NULL;
            shft = ( a + b ) / 2.00;
            slp = ( b - a ) / 2.00;
            break;
        case 2:
            al = be = -0.50;
            if ( ( b - a ) <= temp )
                return NULL;
            shft = ( a + b ) / 2.00;
            slp = ( b - a ) / 2.00;
            break;
        case 3:
            al = be = alpha;
            if ( ( b - a ) <= temp )
                return NULL;
            shft = ( a + b ) / 2.00;
            slp = ( b - a ) / 2.00;
            break;
        case 4:
            al = alpha;
            be = beta;

            if ( ( b - a ) <= temp )
                return NULL;
            shft = ( a + b ) / 2.00;
            slp = ( b - a ) / 2.00;
            break;
        case 5:
            if ( b <= 0.0 )
                return NULL;
            shft = a;
            slp = 1.00 / b;
            al = alpha;
            be = 0.00;
            break;
        case 6:
            if ( b <= 0.0 )
                return NULL;
            shft = a;
            slp = 1.00 / sqrt ( b );
            al = alpha;
            be = 0.00;
            break;
        case 7:
            al = alpha;
            be = 0.00;
            if ( ( b - a ) <= temp )
                return NULL;
            shft = ( a + b ) / 2.00;
            slp = ( b - a ) / 2.00;
            break;
        case 8:
            if ( a + b <= 0.0 )
                return NULL;
            shft = a;
            slp = a + b;
            al = alpha;
            be = beta;
            break;
    }

    p = pow ( slp, al + be + 1.00 );

    for ( k = 0; k < nt; ++k )
    {
        st[k] = shft + slp * t[k];
        l = abs ( ndx[k] );

        if ( l != 0 )
        {
            tmp = p;
            for ( i = l - 1; i <= l - 1 + mlt[k] - 1; ++i )
            {
                swts[i] = wts[i] * tmp;
                tmp *= slp;
            }
        }
    }
    
    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sgqf ( void * data)
/******************************************************************************/
/*
  Purpose:
    SGQF computes knots and weights of a Gauss Quadrature formula.
  Discussion:
    This routine computes all the knots and weights of a Gauss quadrature
    formula with simple knots from the Jacobi matrix and the zero-th
    moment of the weight function, using the Golub-Welsch technique.
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
    Input, double AJ[NT], the diagonal of the Jacobi matrix.
    Input/output, double BJ[NT], the subdiagonal of the Jacobi
    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
    Input, double ZEMU, the zero-th moment of the weight function.
    Output, double T[NT], the knots.
    Output, double WTS[NT], the weights.
*/
{
	const dt2pitit2pit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * aj = s_data->a1;
	ityp * bj = s_data->a2;
	const register ityp zemu = s_data->a3;
	ityp * t = s_data->a4;
	ityp * wts = s_data->a5;
	
    dim_typ i;
    /*
    Exit if the zero-th moment is not positive.
    */
    if ( zemu <= 0.0 )
        return NULL;
    /*
    Set up vectors for IMTQLX.
    */
    for ( i = 0; i < nt; ++i )
        t[i] = aj[i];
    wts[0] = sqrt ( zemu );
    for ( i = 1; i < nt; ++i )
        wts[i] = 0.00;
    /*
    Diagonalize the Jacobi matrix.
    */
    imtqlx ( nt, t, bj, wts );

    for ( i = 0; i < nt; ++i )
        wts[i] *= wts[i];

    return NULL;
}

#endif

