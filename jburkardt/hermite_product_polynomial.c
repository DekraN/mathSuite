#ifndef __DISABLEDEEP_HERMITEPRODUCTPOLYNOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _comp_next_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMP_NEXT_GRLEX returns the next composition in grlex order.
  Discussion:
    Example:
    KC = 3
    #   XC(1) XC(2) XC(3)  Degree
      +------------------------
    1 |  0     0     0        0
      |
    2 |  0     0     1        1
    3 |  0     1     0        1
    4 |  1     0     0        1
      |
    5 |  0     0     2        2
    6 |  0     1     1        2
    7 |  0     2     0        2
    8 |  1     0     1        2
    9 |  1     1     0        2
   10 |  2     0     0        2
      |
   11 |  0     0     3        3
   12 |  0     1     2        3
   13 |  0     2     1        3
   14 |  0     3     0        3
   15 |  1     0     2        3
   16 |  1     1     1        3
   17 |  1     2     0        3
   18 |  2     0     1        3
   19 |  2     1     0        3
   20 |  3     0     0        3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int KC, the number of parts of the composition.
    1 <= KC.
    Input/output, int XC[KC], the current composition.
    Each entry of XC must be nonnegative.
    On return, XC has been replaced by the next composition in the
    grlex order.
*/
{
	const dtpdt * const s_data = data;
	const register dim_typ kc = s_data->a0;
	dim_typ * xc = s_data->a1;
	
    dim_typ i;
    dim_typ im1;
    dim_typ j;
    dim_typ t;
    /*
    Ensure that 1 <= KC.
    */
    if ( kc < 1 )
        return NULL;
    /*
    Ensure that 0 <= XC(I).
    */
    for ( i = 0; i < kc; ++i )
        if ( xc[i] < 0 )
            return NULL;
    /*
    Find I, the index of the rightmost nonzero entry of X.
    */
    i = 0;
    for ( j = kc; 1 <= j; --j )
    {
        if ( 0 < xc[j-1] )
        {
            i = j;
            break;
        }
    }
    /*
    set T = X(I)
    set XC(I) to zero,
    increase XC(I-1) by 1,
    increment XC(KC) by T-1.
    */
    if ( i == 0 )
    {
        xc[kc-1] = 1;
        return NULL;
    }
    else if ( i == 1 )
    {
        t = xc[0] + 1;
        im1 = kc;
    }
    else if ( 1 < i )
    {
        t = xc[i-1];
        im1 = i - 1;
    }

    xc[i-1] = 0;
    ++ xc[im1-1];
    xc[kc-1] += t - 1;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _comp_random_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMP_RANDOM_GRLEX: random composition with degree less than or equal to NC.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int KC, the number of parts in the composition.
    Input, int RANK1, RANK2, the minimum and maximum ranks.
    1 <= RANK1 <= RANK2.
    Input/output, int *SEED, the random number seed.
    Output, int *RANK, the rank of the composition.
    Output, int COMP_RANDOM_GRLEX[KC], the random composition.
*/
{
	const _3dtpipdt * const s_data = data;
	const register dim_typ kc = s_data->a0;
	const register dim_typ rank1 = s_data->a1;
	const register dim_typ rank2 = s_data->a2;
	int * seed = s_data->a3;
	dim_typ * rank = s_data->a4;
	
    int *xc;
    /*
    Ensure that 1 <= KC.
    */
    if ( kc < 1 || rank1<1 || rank2<rank1 )
        return NULL;
    /*
    Choose RANK between RANK1 and RANK2.
    */
    *rank = i4_uniform_ab ( rank1, rank2, seed );
    /*
    Recover the composition of given RANK.
    */
    return comp_unrank_grlex ( kc, *rank );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _comp_rank_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMP_RANK_GRLEX computes the graded lexicographic rank of a composition.
  Discussion:
    The graded lexicographic ordering is used, over all KC-compositions
    for NC = 0, 1, 2, ...
    For example, if KC = 3, the ranking begins:
    Rank  Sum    1  2  3
    ----  ---   -- -- --
       1    0    0  0  0

       2    1    0  0  1
       3    1    0  1  0
       4    1    1  0  1
       5    2    0  0  2
       6    2    0  1  1
       7    2    0  2  0
       8    2    1  0  1
       9    2    1  1  0
      10    2    2  0  0
      11    3    0  0  3
      12    3    0  1  2
      13    3    0  2  1
      14    3    0  3  0
      15    3    1  0  2
      16    3    1  1  1
      17    3    1  2  0
      18    3    2  0  1
      19    3    2  1  0
      20    3    3  0  0
      21    4    0  0  4
      ..   ..   .. .. ..
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int KC, the number of parts in the composition.
    1 <= KC.
    Input, int XC[KC], the composition.
    For each 1 <= I <= KC, we have 0 <= XC(I).
    Output, int COMP_RANK_GRLEX, the rank of the composition.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpdt * const s_data = data;
	const register dim_typ kc = s_data->a0;
	dim_typ * xc = s_data->a1;
	
    dim_typ i;
    dim_typ j;
    dim_typ ks;
    dim_typ n;
    dim_typ nc;
    dim_typ ns;
    dim_typ rank;
    dim_typ tim1;
    dim_typ *xs;
    /*
    Ensure that 1 <= KC.
    */
    if ( kc < 1 )
    {
    	result = USHRT_MAX;
        return &result;
    }
    /*
    Ensure that 0 <= XC(I).
    */
    for ( i = 0; i < kc; ++i )
        if ( xc[i] < 0 )
        {
        	result = USHRT_MAX;
            return &result;
        }
    /*
    NC = sum ( XC )
    */
    nc = i4vec_sum ( kc, xc );
    /*
    Convert to KSUBSET format.
    */
    ns = nc + kc - 1;
    ks = kc - 1;
    xs = ( dim_typ * ) malloc ( ks * sizeof ( dim_typ ) );
    xs[0] = xc[0] + 1;
    for ( i = 2; i < kc; ++i )
        xs[i-1] = xs[i-2] + xc[i-1] + 1;
    /*
    Compute the rank.
    */
    rank = 1;

    for ( i = 1; i <= ks; ++i)
    {
        tim1 = 0 + (i!=1)*xs[i-2];
        if ( tim1 + 1 <= xs[i-1] - 1 )
            for ( j = tim1 + 1; j <= xs[i-1] - 1; ++j )
                rank += i4_choose ( ns - j, ks - i );
    }
    for ( n = 0; n < nc; ++n )
        rank += i4_choose ( n + kc - 1, n );
    free ( xs );

	result = rank;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _comp_unrank_grlex ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMP_UNRANK_GRLEX computes the composition of given grlex rank.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int KC, the number of parts of the composition.
    1 <= KC.
    Input, int RANK, the rank of the composition.
    1 <= RANK.
    Output, int COMP_UNRANK_GRLEX[KC], the composition XC of the given rank.
    For each I, 0 <= XC[I] <= NC, and
    sum ( 1 <= I <= KC ) XC[I] = NC.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ kc = a_data[0];
	const register dim_typ rank = a_data[1];
	
    dim_typ i;
    dim_typ j;
    dim_typ ks;
    dim_typ nc;
    dim_typ nksub;
    dim_typ ns;
    dim_typ r;
    dim_typ rank1;
    dim_typ rank2;
    dim_typ *xc;
    dim_typ *xs;
    /*
    Ensure that 1 <= KC.
    */
    if ( kc < 1 || rank < 1 )
        return NULL;
    /*
    Ensure that 1 <= RANK.
    */
    /*
    Determine the appropriate value of NC.
    Do this by adding up the number of compositions of sum 0, 1, 2,
    ..., without exceeding RANK.  Moreover, RANK - this sum essentially
    gives you the rank of the composition within the set of compositions
    of sum NC.  And that's the number you need in order to do the
    unranking.
    */
    rank1 = 1;
    nc = -1;
    for ( ; ; )
    {
        ++ nc;
        r = i4_choose ( nc + kc - 1, nc );
        if ( rank < rank1 + r )
            break;
        rank1 += r;
    }

    rank2 = rank - rank1;
    /*
    Convert to KSUBSET format.
    Apology: an unranking algorithm was available for KSUBSETS,
    but not immediately for compositions.  One day we will come back
    and simplify all this.
    */
    ks = kc - 1;
    ns = nc + kc - 1;
    xs = ( dim_typ * ) malloc ( ks * sizeof ( dim_typ ) );

    nksub = i4_choose ( ns, ks );

    j = 1;

    for ( i = 1; i <= ks; ++i )
    {
        r = i4_choose ( ns - j, ks - i );

        while ( r <= rank2 && 0 < r )
        {
            rank2 = rank2 - r;
            ++ j;
            r = i4_choose ( ns - j, ks - i );
        }
        xs[i-1] = j;
        ++ j;
    }
    /*
    Convert from KSUBSET format to COMP format.
    */
    xc = ( dim_typ * ) malloc ( kc * sizeof ( dim_typ ) );
    xc[0] = xs[0] - 1;
    for ( i = 2; i < kc; ++i )
        xc[i-1] = xs[i-1] - xs[i-2] - 1;
    xc[kc-1] = ns - xs[ks-1];

    free ( xs );
    return xc;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _hep_coefficients ( void * data)
/******************************************************************************/
/*
  Purpose:
    HEP_COEFFICIENTS: coefficients of Hermite polynomials He(n,x).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
  First terms:
    N/K     0     1      2      3       4     5      6    7      8    9   10
     0      1
     1      0     1
     2     -1     0      1
     3      0    -3      0      1
     4      3     0     -6      0       1
     5      0    15      0    -10       0     1
     6    -15     0     45      0     -15     0      1
     7      0  -105      0    105       0   -21      0     1
     8    105     0   -420      0     210     0    -28     0      1
     9      0   945      0  -1260       0   378      0   -36      0   1
    10   -945     0   4725      0   -3150     0    630     0    -45   0    1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 October 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.
    Output, int *O, the number of coefficients.
    Output, double C[(N+2)/2], the coefficients of the Hermite
    polynomial of degree N.

    Output, int F[(N+2)/2], the exponents.
*/
{
	const dtpit2pdt * const s_data = data;
	
	const register dim_typ n = s_data->a0; 
	ityp * c = s_data->a1;
	dim_typ * o = s_data->a2;
	dim_typ * f = s_data->a3;
	
    ityp *ct;
    dim_typ i, j, k;

    ct = ( ityp * ) malloc ( (n+1)*(n+1) * sizeof ( ityp ) );

    for ( i = 0; i <= n; ++i )
        for ( j = 0; j <= n; ++j )
            ct[i+j*(n+1)] = 0.00;

    ct[0+0*(n+1)] = 1.00;

    if ( 0 < n )
    {
        ct[1+1*(n+1)] = 1.00;

        for ( i = 2; i <= n; ++i )
        {
            ct[i+0*(n+1)] =                       - ( ityp ) ( i - 1 ) * ct[i-2+0*(n+1)];
            for ( j = 1; j <= i - 2; ++j)
                ct[i+j*(n+1)] = ct[i-1+(j-1)*(n+1)] - ( ityp ) ( i - 1 ) * ct[i-2+j*(n+1)];
            ct[i+(i-1)*(n+1)] =   ct[i-1+(i-2)*(n+1)];
            ct[i+ i   *(n+1)] =   ct[i-1+(i-1)*(n+1)];
        }
    }
    /*
    Extract the nonzero data from the alternating columns of the last row.
    */
    *o = ( n + 2 ) / 2;

    k = *o;
    for ( j = n; 0 <= j; j -= 2 )
    {
        -- k;
        c[k] = ct[n+j*(n+1)];
        f[k] = j;
    }

    free ( ct );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hep_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    HEP_VALUE evaluates the Hermite polynomials He(i,x).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
    1
    X
    X^2  -  1
    X^3  -  3 X
    X^4  -  6 X^2 +   3
    X^5  - 10 X^3 +  15 X
    X^6  - 15 X^4 +  45 X^2 -   15
    X^7  - 21 X^5 + 105 X^3 -  105 X
    X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
    X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
    X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, int N, the number of evaluation points.
    Input, int O, the degree of the polynomial.
    Input, double X[N], the evaluation points.
    Output, double HEP_VALUE[N], the value of the Hermite polynomial
    of degree N at the points X.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
 	const register dim_typ o = s_data->a1;
 	ityp * x = s_data->a2;
	
    dim_typ i, j;
    ityp *v;
    ityp *vtable;

    vtable = ( ityp * ) malloc ( n*(o+1) * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        vtable[i+0*n] = 1.00;

    if ( 1 <= o )
    {
        for ( i = 0; i < n; ++i )
            vtable[i+1*n] = x[i];


        for ( j = 2; j <= o; ++j )
            for ( i = 0; i < n; ++i )
                vtable[i+j*n] = x[i] * vtable[i+(j-1)*n]- ( ityp ) ( j - 1 ) * vtable[i+(j-2)*n];
    }

    v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        v[i] = vtable[i+o*n];

    free ( vtable );

    return v;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hepp_to_polynomial ( void * data)
/******************************************************************************/
/*
  Purpose:
    HEPP_TO_POLYNOMIAL writes a Hermite Product Polynomial as a polynomial.
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
    For example, if
      M = 3,
      L = ( 1, 0, 2 ),
    then
      He(1,0,2)(X,Y,Z)
      = He(1)(X) * He(0)(Y) * He(2)(Z)
      = X * 1 * ( Z^3-3Z)
      = - 3XZ + X Z^3
    so
      O = 2 (2 nonzero terms)
      C = -3.0
           1.0
      E =  8   <-- index in 3-space of exponent (1,0,1)
          23   <-- index in 3-space of exponent (1,0,3)
    The output value of O is no greater than
      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int L[M], the index of each Hermite product
    polynomial factor.  0 <= L(*).
    Input, int O_MAX, an upper limit on the size of the
    output arrays.
      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2.
    Output, int *O, the "order" of the polynomial product.
    Output, double C[O], the coefficients of the polynomial product.
    Output, int E[O], the indices of the exponents of the
    polynomial product.
*/
{
	const _2dtpi2pdtpit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ o_max = s_data->a1;
	int * e = s_data->a2; 
	dim_typ * l = s_data->a3;
	dim_typ * o = s_data->a4;
	ityp * c = s_data->a5;
	
	
    ityp *c1;
    ityp *c2;
    dim_typ *e1;
    dim_typ *e2;
    dim_typ *f2;
    dim_typ i;
    dim_typ i1;
    dim_typ i2;
    dim_typ j1;
    dim_typ j2;
    dim_typ o1;
    dim_typ o2;
    int *p;
    int *pp;

    c1 = ( ityp * ) malloc ( o_max * sizeof ( ityp ) );
    c2 = ( ityp * ) malloc ( o_max * sizeof ( ityp ) );
    e1 = ( dim_typ * ) malloc ( o_max * sizeof ( dim_typ ) );
    e2 = ( dim_typ * ) malloc ( o_max * sizeof ( dim_typ ) );
    f2 = ( dim_typ * ) malloc ( o_max * sizeof ( dim_typ ) );
    pp = ( int * ) malloc ( m * sizeof ( int ) );

    o1 = e1[0] = 1;
    c1[0] = 1.00;
    /*
    Implicate one factor at a time.
    */
    for ( i = 0; i < m; ++i)
    {
        hep_coefficients ( l[i], &o2, c2, f2 );

        *o = 0;

        for ( j2 = 0; j2 < o2; ++j2)
            for ( j1 = 0; j1 < o1; ++j1 )
            {
                c[*o] = c1[j1] * c2[j2];
                if ( 0 < i )
                    p = mono_unrank_grlex ( i, e1[j1] );

                for ( i2 = 0; i2 < i; ++i2)
                    pp[i2] = p[i2];
                pp[i] = f2[j2];
                e[*o] = mono_rank_grlex ( i + 1, pp );
                ++ *o;
                if ( 0 < i )
                    free ( p );
            }

        polynomial_sort ( *o, c, e );
        polynomial_compress ( *o, c, e, o, c, e );

        o1 = *o;
        for ( i1 = 0; i1 < *o; ++i1 )
        {
            c1[i1] = c[i1];
            e1[i1] = e[i1];
        }
    }

    free ( c1 );
    free ( c2 );
    free ( e1 );
    free ( e2 );
    free ( f2 );
    free ( pp );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _hepp_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    HEPP_VALUE evaluates a Hermite Product Polynomial.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of evaluation points.
    Input, int O[M], the degree of the polynomial factors.
    0 <= O(*).
    Input, double X[M*N], the evaluation points.
    Output, double HEPP_VALUE[N], the value of the Hermite Product
    Polynomial of degree O at the points X.
*/
{
	const pit2dtpdt * const s_data = data;
	
	ityp * x = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	dim_typ * o = s_data->a3; 
	
    dim_typ i, j;    
    ityp *v;
    ityp *vi;
    ityp *xi;

    v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        v[j] = 1.00;

    xi = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
    {
        for ( j = 0; j < n; ++j )
            xi[j] = x[i+j*m];
        vi = hep_value ( n, o[i], xi );
        for ( j = 0; j < n; ++j )
            v[j] *= vi[j];
        free ( vi );
    }

    free ( xi );

    return v;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4_uniform_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
  Discussion:
    The pseudorandom number should be uniformly distributed
    between A and B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 May 2012
  Author:
    John Burkardt
  Reference:
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.
    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.
    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.
    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.
  Parameters:
    Input, int A, B, the limits of the interval.
    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.
    Output, int I4_UNIFORM_AB, a number between A and B.
*/
{
	static int result = INT_MAX;
	
	const _2ipi * const s_data = data;
	register int a = s_data->a0;
	register int b = s_data->a1;
	int * seed = s_data->a2;
	
    int c;
    int k;
    ityp r;
    int value;

    if ( *seed == 0 )
    {
    	result = INT_MAX;
        return &result;
    }
    /*
    Guaranteee A <= B.
    */
    if ( b < a )
    {
        c = a;
        a = b;
        b = c;
    }

    k = *seed / 127773;
    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
        *seed += i4_huge;

    r = ( ityp ) ( *seed ) * 4.656612875E-10;
    /*
    Scale R to lie between A-0.5 and B+0.5.
    */
    r = ( 1.00 - r ) * ( ( ityp ) ( a ) - 0.50 )+         r   * ( ( ityp ) ( b ) + 0.50 );
    /*
    Round R to the nearest integer.
    */
    /*
    Guarantee that A <= VALUE <= B.
    */
    
    result = round ( r )<a?a:b;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _polynomial_compress ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLYNOMIAL_COMPRESS compresses a polynomial.
  Discussion:
    The function polynomial_sort ( ) should be called first.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 January 2014
  Author:
    John Burkardt
  Parameters:
    Input, int O1, the "order" of the polynomial.
    Input, double C1[O1], the coefficients of the polynomial.
    Input, int E1[O1], the indices of the exponents of
    the polynomial.
    Output, int *O2, the "order" of the polynomial.
    Output, double C2[*O2], the coefficients of the polynomial.
    Output, int E2[*O2], the indices of the exponents of
    the polynomial.
*/
{
	const dtpitpipdtpitpi * const s_data = data;
	const register dim_typ o1 = s_data->a0;
	ityp * c1 = s_data->a1;
	int * e1 = s_data->a2;
	dim_typ * o2 = s_data->a3;
	ityp * c2 = s_data->a4;
	int * e2 = s_data->a5;
	
    dim_typ get, put;
    const register ityp r8_epsilon_sqrt = 0.1490116119384766E-07;

    get = put = 0;

    while ( get < o1 )
    {
        ++ get;

        if ( fabs ( c1[get-1] ) <= r8_epsilon_sqrt )
            continue;

        if ( 0 == put )
        {
            ++ put;
            c2[put-1] = c1[get-1];
            e2[put-1] = e1[get-1];
        }
        else
            if ( e2[put-1] == e1[get-1] )
                c2[put-1] = c2[put-1] + c1[get-1];
            else
            {
                ++ put;
                c2[put-1] = c1[get-1];
                e2[put-1] = e1[get-1];
            }
    }

    *o2 = put;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _polynomial_sort ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLYNOMIAL_SORT sorts the information in a polynomial.
  Discussion
    The coefficients C and exponents E are rearranged so that
    the elements of E are in ascending order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int O, the "order" of the polynomial.
    Input/output, double C[O], the coefficients of the polynomial.
    Input/output, int E[O], the indices of the exponents of
    the polynomial.
*/
{
	const dtpitpi * const s_data = data;
	const register dim_typ o = s_data->a0;
	ityp * c = s_data->a1;
	int * e = s_data->a2;
	
    int *indx = i4vec_sort_heap_index_a ( o, e );
    i4vec_permute ( o, indx, e );
    r8vec_permute ( o, indx, 0, c );
    free ( indx );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _polynomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    POLYNOMIAL_VALUE evaluates a polynomial.
  Discussion:
    The polynomial is evaluated term by term, and no attempt is made to
    use an approach such as Horner's method to speed up the process.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int O, the "order" of the polynomial.
    Input, double C[O], the coefficients of the polynomial.
    Input, int E[O], the indices of the exponents
    of the polynomial.
    Input, int NX, the number of evaluation points.
    Input, double X[M*NX], the coordinates of the evaluation points.
    Output, double POLYNOMIAL_VALUE[NX], the value of the polynomial at X.
*/
{
	const pit2dtpitpdtdt * const s_data = data;
	
	ityp * c = s_data->a0;
	const register dim_typ m = s_data->a1; 
	const register dim_typ o = s_data->a2;
	ityp * x = s_data->a3;
	dim_typ * e = s_data->a4;
	const register dim_typ nx = s_data->a5;
	
    int *f;
    dim_typ j, k;
    ityp *p;
    ityp *v;

    p = ( ityp * ) malloc ( nx * sizeof ( ityp ) );

    for ( k = 0; k < nx; ++k )
        p[k] = 0.00;

    for ( j = 0; j < o; ++j)
    {
        f = mono_unrank_grlex ( m, e[j] );
        v = mono_value ( m, nx, f, x );
        for ( k = 0; k < nx; ++k )
            p[k] += c[j] * v[k];
        free ( f );
    }
    return p;
}

#endif
