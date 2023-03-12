#ifndef __DISABLEDEEP_NINTLIB

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   box_nd (ityp func ( dim_typ dim_num, ityp x[] ), const register dim_typ dim_num, const register dim_typ order, ityp xtab[static order], ityp weight[static order], int *eval_num )
/******************************************************************************/
/*
  Purpose:
    BOX_ND estimates a multidimensional integral using a product rule.
  Discussion:
    The routine creates a DIM_NUM-dimensional product rule from a 1D rule
    supplied by the user.  The routine is fairly inflexible.  If
    you supply a rule for integration from -1 to 1, then your product
    box must be a product of DIM_NUM copies of the interval [-1,1].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 February 2007
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, double FUNC ( int dim_num, double x[] ), evaluates
    the function to be integrated.
    Input, int DIM_NUM, the spatial dimension.
    Input, int ORDER, the number of points used in the 1D rule.
    Input, double XTAB[ORDER], the abscissas of the 1D rule.
    Input, double WEIGHT[ORDER], the weights of the 1D rule.
    Output, int *EVAL_NUM, the number of function evaluations.
    Output, double BOX_ND, the approximate value of the integral.
*/
{
    dim_typ dim;
    int *indx;
    dim_typ k;
    ityp result;
    ityp w;
    ityp *x;

    *eval_num = 0;

    if ( dim_num < 1 || order < 1 )
        return MAX_VAL;

    k = 0;
    result = 0.00;

    indx = ( int * ) malloc ( dim_num * sizeof ( int ) );
    x = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    for ( ; ; )
    {
        tuple_next ( 1, order, dim_num, &k, indx );

        if ( k == 0  )
            break;

        w = 1.00;
        for ( dim = 0; dim < dim_num; ++dim )
        {
            w *= weight[indx[dim]-1];
            x[dim] = xtab[indx[dim]-1];
        }

        result += w * func ( dim_num, x );
        ++ *eval_num;
    }

    free ( indx );
    free ( x );

    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   monte_carlo_nd ( ityp func ( dim_typ dim_num, ityp x[] ), const register dim_typ dim_num,
  ityp a[static dim_num], ityp b[static dim_num ], const register dim_typ eval_num, int *seed )
/******************************************************************************/
/*
  Purpose:
    MONTE_CARLO_ND estimates a multidimensional integral using Monte Carlo.
  Discussion:
    Unlike the other routines, this routine requires the user to specify
    the number of function evaluations as an INPUT quantity.
    No attempt at error estimation is made.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 February 2007
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, double FUNC ( int dim_num, double x[] ), evaluates
    the function to be integrated.
    Input, int DIM_NUM, the spatial dimension.
    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits
    Input, int EVAL_NUM, the number of function evaluations.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double MONTE_CARLO_ND, the approximate value of the integral.
*/
{
    dim_typ i, dim;
    ityp result;
    ityp volume;
    ityp *x;

    result = 0.00;

    for ( i = 0; i < eval_num; ++i)
    {
        x = r8vec_uniform_01_new ( dim_num, seed );
        result += func ( dim_num, x );
        free ( x );
    }

    volume = 1.00;
    for ( dim = 0; dim < dim_num; ++dim )
        volume *= ( b[dim] - a[dim] );

    return result * volume / ( ityp ) ( eval_num );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   p5_nd ( ityp func ( int dim_num, ityp x[] ), const register dim_typ dim_num,ityp a[static dim_num], ityp b[static dim_num], int *eval_num )
/******************************************************************************/
/*
  Purpose:
    P5_ND estimates a multidimensional integral with a formula of exactness 5.
  Discussion:
    The routine uses a method which is exact for polynomials of total
    degree 5 or less.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 February 2007
  Author:
    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, double FUNC ( int dim_num, double x[] ), evaluates
    the function to be integrated.
    Input, int DIM_NUM, the spatial dimension.
    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
    Output, int *EVAL_NUM, the number of function evaluations.
    Output, double P5_ND, the approximate value of the integral.
*/
{
    ityp a0;
    ityp a1;
    ityp a2;
    ityp a3;
    ityp a4;
    ityp a5;
    ityp en;
    dim_typ dim, i, j;
    ityp result;
    ityp sum1;
    ityp sum2;
    ityp sum3;
    ityp volume;
    ityp *work;

    *eval_num = 0;

    if ( dim_num < 1 )
        return MAX_VAL;

    a2 = 25.00 / 324.00;
    a3 = sqrt ( 0.60 );
    en = ( ityp ) ( dim_num );
    a0 = ( 25.00 * en * en - 115.00 * en + 162.00 ) / 162.00;
    a1 = ( 70.00 - 25.00 * en ) / 162.00;

    volume = 1.00;
    work = ( double * ) malloc ( dim_num * sizeof ( double ) );
    for ( dim = 0; dim < dim_num; ++dim )
    {
        volume *= ( b[dim] - a[dim] );
        work[dim] = 0.50 * ( a[dim] + b[dim] );
    }

    result = 0.00;
    if ( volume == 0.00 )
        return MAX_VAL;


    sum1 = a0 * func ( dim_num, work );
    *eval_num = *eval_num + 1;

    sum2 = sum3 = 0.00;

    for ( i = 0; i < dim_num; ++i )
    {
        work[i] = 0.50 * ( ( a[i] + b[i] ) + a3 * ( b[i] - a[i] ) );
        sum2 = sum2 + func ( dim_num, work );
        ++ *eval_num;

        work[i] = 0.50 * ( ( a[i] + b[i] ) - a3 * ( b[i] - a[i] ) );
        sum2 = sum2 + func ( dim_num, work );
        ++ *eval_num;

        work[i] = 0.50 * ( a[i] + b[i] );
    }

    if ( 1 < dim_num )
    {
        a4 = a3;

        for ( ; ; )
        {
            for ( i = 0; i < dim_num - 1; ++i )
            {
                work[i] = 0.5 * ( ( a[i] + b[i] ) + a4 * ( b[i] - a[i] ) );
                a5 = a3;

                for ( ; ; )
                {
                    for ( j = i + 1; j < dim_num; ++j )
                    {
                        work[j] = 0.50 * ( ( a[j] + b[j] ) + a5 * ( b[j] - a[j] ) );
                        sum3 = sum3 + func ( dim_num, work );
                        *eval_num = *eval_num + 1;
                        work[j] = 0.50 * ( a[j]+ b[j] );
                    }

                    a5 = -a5;

                    if ( 0.0 <= a5 )
                        break;

                }
                work[i] = 0.50 * ( a[i] + b[i] );
            }

            a4 *= -1;

            if ( 0.00 <= a4 )
                break;
        }
    }

    free ( work );
    return volume * ( sum1 + a1 * sum2 + a2 * sum3 );
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   romberg_nd ( ityp func ( dim_typ dim_num, ityp x[] ), ityp a[],
  ityp b[], const register dim_typ dim_num, dim_typ sub_num[static dim_num], const register dim_typ it_max, const register ityp tol, int *ind,dim_typ *eval_num )
/******************************************************************************/
/*
  Purpose:
    ROMBERG_ND estimates a multidimensional integral using Romberg integration.
  Discussion:
    The routine uses a Romberg method based on the midpoint rule.
    In the reference, this routine is called "NDIMRI".
    Thanks to Barak Bringoltz for pointing out problems in a previous
    FORTRAN90 implementation of this routine.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 February 2007
  Author:
    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, double FUNC ( int dim_num, double x[] ), evaluates
    the function to be integrated.
    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
    Input, int DIM_NUM, the spatial dimension.
    Input, int SUB_NUM[DIM_NUM], the number of subintervals into
    which the I-th integration interval (A(I), B(I)) is
    initially subdivided.  SUB_NUM(I) must be greater than 0.
    Input, int IT_MAX, the maximum number of iterations to
    be performed.  The number of function evaluations on
    iteration J is at least J**DIM_NUM, which grows very rapidly.
    IT_MAX should be small!
    Input, double TOL, an error tolerance for the approximation
    of the integral.
    Output, int *IND, error return flag.
    IND = -1 if the error tolerance could not be achieved.
    IND = 1 if the error tolerance was achieved.
    Output, int *EVAL_NUM, the number of function evaluations.
    Output, double ROMBERG_ND, the approximate value of the integral.
  Local Parameters:
    Local, int IWORK[DIM_NUM], a pointer used to generate all the
    points X in the product region.
    Local, int IWORK2[IT_MAX], a counter of the number of points
    used at each step of the Romberg iteration.
    Local, int SUB_NUM2[DIM_NUM], the number of subintervals used
    in each direction, a refinement of the user's input SUB_NUM.
    Local, double TABLE[IT_MAX], the difference table.
    Local, double X[DIM_NUM], an evaluation point.
*/
{
    dim_typ dim;
    ityp en;
    ityp factor;
    dim_typ i;
    dim_typ it;
    int *iwork;
    int *iwork2;
    int kdim;
    int ll;
    ityp result;
    ityp result_old;
    ityp rnderr;
    ityp submid;
    int *sub_num2;
    ityp sum1;
    ityp weight;
    ityp *table;
    ityp *x;

    *eval_num = 0;

    if ( dim_num < 1 || it_max < 1 )
        return MAX_VAL;

    for ( i = 0; i < dim_num; ++i )
        if ( sub_num[i] <= 0 )
            return MAX_VAL;

    iwork = ( int * ) malloc ( dim_num * sizeof ( int ) );
    iwork2 = ( int * ) malloc ( it_max * sizeof ( int ) );
    sub_num2 = ( int * ) malloc ( dim_num * sizeof ( int ) );
    table = ( ityp * ) malloc ( it_max * sizeof ( ityp ) );
    x = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    *ind = 0;
    rnderr = r8_epsilon ( );
    iwork2[0] = 1;

    for ( dim = 0; dim < dim_num; ++dim )
        sub_num2[dim] = sub_num[dim];
    if ( 1 < it_max )
        iwork2[1] = 2;

    it = 1;

    for ( ; ; )
    {
        sum1 = 0.00;

        weight = 1.0;
        for ( dim = 0; dim < dim_num; ++dim )
        {
        weight = weight * ( b[dim] - a[dim] ) / ( double ) sub_num2[dim];
        }
        /*
        Generate every point X in the product region, and evaluate F(X).
        */
        for ( dim = 0; dim < dim_num; ++dim )
            iwork[dim] = 1;

        for ( ; ; )
        {
            for ( dim = 0; dim < dim_num; ++dim )
                x[dim] =( ( ityp ) ( (sub_num2[dim]<<1) - (iwork[dim]<<1) + 1 ) * a[dim]+ ( ityp ) (                   + (iwork[dim]<<1) - 1 ) * b[dim] )/ ( ityp ) ( (sub_num2[dim]<<1)                     );


            sum1 += func ( dim_num, x );
            *eval_num = *eval_num + 1;

            kdim = dim_num;

            while ( 0 < kdim )
            {
                if ( iwork[kdim-1] < sub_num2[kdim-1] )
                {
                    iwork[kdim-1] = iwork[kdim-1] + 1;
                    break;
                }
                iwork[kdim-1] = 1;
                -- kdim;
            }

            if ( kdim == 0 )
                break;
        }
        /*
        Done with summing.
        */
        table[it-1] = weight * sum1;

        if ( it <= 1 )
        {
            result = table[0];
            result_old = result;

            if ( it_max <= it )
            {
                *ind = 1;
                break;
            }

            ++ it;
            for ( dim = 0; dim < dim_num; ++dim )
                sub_num2[dim] = iwork2[it-1] * sub_num2[dim];
            continue;
        }
        /*
        Compute the difference table for Richardson extrapolation.
        */
        for ( ll = 2; ll <= it; ++ll )
        {
            i = it + 1 - ll;
            factor = ( ityp ) ( iwork2[i-1] * iwork2[i-1] )/ ( ityp ) ( iwork2[it-1] * iwork2[it-1] - iwork2[i-1] * iwork2[i-1] );
            table[i] += ( table[i] - table[i-1] ) * factor;
        }

        result = table[0];
        /*
        Terminate successfully if the estimated error is acceptable.
        */
        if ( abs ( result - result_old ) <= abs ( result * ( tol + rnderr ) ) )
        {
            *ind = 1;
            break;
        }
        /*
        Terminate unsuccessfully if the iteration limit has been reached.
        */
        if ( it_max <= it )
        {
            *ind = -1;
            break;
        }
        /*
        Prepare for another step.
        */
        result_old = result;
        ++ it;

        iwork2[it-1] = ( int ) ( 1.50 * ( ityp ) ( iwork2[it-2] ) );

        for ( dim = 0; dim < dim_num; ++dim )
            sub_num2[dim] = ( int ) ( 1.50 * ( ityp ) ( sub_num2[dim] ) );
    }

    free ( iwork );
    free ( iwork2 );
    free ( sub_num2 );
    free ( table );
    free ( x );

    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   sample_nd ( ityp func ( dim_typ dim_num, ityp x[] ), const register dim_typ k1, const register dim_typ k2,
  const register dim_typ dim_num, ityp est1[static k2], ityp err1[static k2], ityp dev1[static k2], ityp est2[static k2],ityp err2[static k2], ityp dev2[static k2], int *eval_num )
/******************************************************************************/
/*
  Purpose:
    SAMPLE_ND estimates a multidimensional integral using sampling.
  Discussion:
    This routine computes two sequences of integral estimates, EST1
    and EST2, for indices K going from K1 to K2.  These estimates are
    produced by the generation of 'random' abscissas in the region.
    The process can become very expensive if high accuracy is needed.
    The total number of function evaluations is
    4*(K1^DIM_NUM+(K1+1)^DIM_NUM+...+(K2-1)^DIM_NUM+K2^DIM_NUM), and K2
    should be chosen so as to make this quantity reasonable.
    In most situations, EST2(K) are much better estimates than EST1(K).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 2007
  Author:
    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, double FUNC ( int dim_num, double x[] ), evaluates
    the function to be integrated.
    Input, int K1, the beginning index for the iteration.
    1 <= K1 <= K2.
    Input, int K2, the final index for the iteration.  K1 <= K2.
    Increasing K2 increases the accuracy of the calculation,
    but vastly increases the work and running time of the code.
    Input, int DIM_NUM, the spatial dimension.  1 <= DIM_NUM <= 10.
    Output, double EST1[K2].  Entries K1 through K2 contain
    successively better estimates of the integral.
    Output, double ERR1[K2].  Entries K1 through K2 contain
    the corresponding estimates of the integration errors.
    Output, double DEV1[K2].  Entries K1 through K2 contain
    estimates of the reliability of the the integration.
    If consecutive values DEV1(K) and DEV1(K+1) do not differ
    by more than 10 percent, then ERR1(K) can be taken as
    a reliable upper bound on the difference between EST1(K)
    and the true value of the integral.
    Output, double EST2[K2].  Entries K2 through K2 contain
    successively better estimates of the integral.
    Output, double ERR2[K2].  Entries K2 through K2 contain
    the corresponding estimates of the integration errors.
    Output, double DEV2[K2].  Entries K2 through K2 contain
    estimates of the reliability of the the integration.
    If consecutive values DEV2(K) and DEV2(K+2) do not differ
    by more than 10 percent, then ERR2(K) can be taken as
    a reliable upper bound on the difference between EST2(K)
    and the true value of the integral.
    Output, int *EVAL_NUM, the number of function evaluations.
*/
{
    # define DIM_MAX 10

    ityp ak;
    ityp ak1;
    ityp akn;
    ityp al[DIM_MAX] =
    {
        0.4142135623730950,
        0.7320508075688773,
        0.2360679774997897,
        0.6457513110645906,
        0.3166247903553998,
        0.6055512754639893,
        0.1231056256176605,
        0.3589989435406736,
        0.7958315233127195,
        0.3851648071345040
    };
    ityp b;
    ityp *be;
    ityp bk;
    ityp d1;
    ityp d2;
    ityp *dex;
    dim_typ dim;
    ityp g;
    ityp *ga;
    dim_typ i;
    dim_typ j;
    dim_typ k;
    int key;
    bool more;
    ityp *p1;
    ityp *p2;
    ityp *p3;
    ityp *p4;
    ityp s1;
    ityp s2;
    ityp t;
    ityp y1;
    ityp y2;
    ityp y3;
    ityp y4;

    *eval_num = 0;
    /*
    Check input.
    */
    if ( dim_num < 1 || DIM_MAX < dim_num || k1 < 1 || k2 < k1 )
        return;

    be = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );
    dex = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );
    ga = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );
    p1 = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );
    p2 = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );
    p3 = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );
    p4 = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    for ( dim = 0; dim < dim_num; ++dim)
    {
        be[dim] = ga[dim] = al[dim];
        dex[dim] = 0.00;
    }

    for ( k = k1; k <= k2; ++k )
    {
        ak = ( ityp ) ( k );
        key = 0;
        ak1 = ak - 1.10;
        s1 = d1 = d2 = 0.00;
        akn = pow ( ak, dim_num );
        t = sqrt ( pow ( ak, dim_num ) ) * ak;
        bk = 1.00 / ak;

        for ( ; ; )
        {
            ++ key;

            if ( key != 1 )
            {
                -- key;
                more = 0;
                for ( j = 0; j < dim_num; ++j )
                {
                    if ( dex[j] <= ak1 )
                    {
                        dex[j] = dex[j] + 1.00;
                        more = 1;
                        break;
                    }

                    dex[j] = 0.0;
                }

                if ( !more )
                    break;
            }

            for ( i = 0; i < dim_num; ++i )
            {
                b = be[i] + al[i];
                if ( 1.00 < b )
                    -- b;

                g = ga[i] + b;
                if ( 1.00 < g )
                    -- g;

                be[i] = b + al[i];
                if ( 1.00 < be[i] )
                    -- be[i];

                ga[i] = be[i] + g;
                if ( 1.00 < ga[i] )
                    -- ga[i];

                p1[i] = ( dex[i] + g ) * bk;
                p2[i] = ( dex[i] + 1.00 - g ) * bk;
                p3[i] = ( dex[i] + ga[i] ) * bk;
                p4[i] = ( dex[i] + 1.00 - ga[i] ) * bk;
            }

            y1 = func ( dim_num, p1 );
            ++ *eval_num ;
            /*
            There may be an error in the next two lines,
            but oddly enough, that is how the original reads
            */
            y3 = func ( dim_num, p2 );
            ++ *eval_num;
            y2 = func ( dim_num, p3 );
            ++ *eval_num;
            y4 = func ( dim_num, p4 );
            ++ *eval_num;

            s1 = s1 + y1 + y2;
            d1 = d1 + ( y1 - y2 ) * ( y1 - y2 );
            s2 = s2 + y3 + y4;
            d2 = d2 + ( y1 + y3 - y2 - y4 ) * ( y1 + y3 - y2 - y4 );
        }

        est1[k-1] = 0.50 * s1 / akn;
        err1[k-1] = 1.50 * sqrt ( d1 ) / akn;
        dev1[k-1] = err1[k-1] * t;
        est2[k-1] = 0.25 * ( s1 + s2 ) / akn;
        err2[k-1] = 0.75 * sqrt ( d2 ) / akn;
        dev2[k-1] = err2[k-1] * t * ak;
    }

    free ( be );
    free ( dex );
    free ( ga );
    free ( p1 );
    free ( p2 );
    free ( p3 );
    free ( p4 );

    return;
    # undef DIM_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sum2_nd ( ityp func ( dim_typ dim_num, ityp x[] ), ityp xtab[],ityp weight[], dim_typ order[], const register dim_typ dim_num, int *eval_num )
/******************************************************************************/
/*
  Purpose:
    SUM2_ND estimates a multidimensional integral using a product rule.
  Discussion:
    The routine uses a product rule supplied by the user.
    The region may be a product of any combination of finite,
    semi-infinite, or infinite intervals.
    For each factor in the region, it is assumed that an integration
    rule is given, and hence, the region is defined implicitly by
    the integration rule chosen.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 February 2007
  Author:
    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, double FUNC ( int dim_num, double x[] ), evaluates
    the function to be integrated.
    Input, double XTAB[DIM_NUM*ORDER_MAX].  XTAB(I,J) is the
    I-th abscissa of the J-th rule.
    Input, double WEIGHT[DIM_NUM*ORDER_MAX].  WEIGHT(I,J) is the
    I-th weight for the J-th rule.
    Input, int ORDER[DIM_NUM].  ORDER(I) is the number of
    abscissas to be used in the J-th rule.  ORDER(I) must be
    greater than 0 and less than or equal to ORDER_MAX.
    Input, int DIM_NUM, the spatial dimension.
    Output, int EVAL_NUM, the number of function evaluations.
    Output, double SUM2_ND, the approximate value of the integral.
*/
{
    dim_typ dim;
    dim_typ i;
    int *iwork;
    dim_typ k;
    dim_typ m1;
    ityp result;
    ityp w1;
    ityp *work;
    /*
    Default values.
    */
    result = 0.00;
    *eval_num = 0;

    if ( dim_num < 1 )
        return MAX_VAL;


    for ( i = 0; i < dim_num; ++i )
        if ( order[i] < 1 )
            return MAX_VAL;

    iwork = ( int * ) malloc ( dim_num * sizeof ( int ) );
    work = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    for ( dim = 0; dim < dim_num; ++dim )
    iwork[dim] = 1;

    for ( ; ; )
    {
        k = 1;
        w1 = 1.00;
        for ( i = 0; i < dim_num; ++i )
        {
            m1 = iwork[i];
            work[i] = xtab[i+(m1-1)*dim_num];
            w1 = w1 * weight[i+(m1-1)*dim_num];
        }

        result += w1 * func ( dim_num, work );
        ++ *eval_num;

        while ( iwork[k-1] == order[k-1] )
        {
            iwork[k-1] = 1;
            k = k + 1;

            if ( dim_num < k )
                return result;
        }
        ++ iwork[k-1];
    }

    free ( iwork );
    free ( work );

    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tuple_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    TUPLE_NEXT computes the next element of a tuple space.
  Discussion:
    The elements are N vectors.  Each entry is constrained to lie
    between M1 and M2.  The elements are produced one at a time.
    The first element is
   (M1,M1,...,M1),
    the second element is
   (M1,M1,...,M1+1),
    and the last element is
   (M2,M2,...,M2)
    Intermediate elements are produced in lexicographic order.
  Example:
    N = 2, M1 = 1, M2 = 3
    INPUT        OUTPUT
    -------      -------
    Rank  X      Rank   X
    ----  ---    -----  ---
    0     * *    1      1 1
    1     1 1    2      1 2
    2     1 2    3      1 3
    3     1 3    4      2 1
    4     2 1    5      2 2
    5     2 2    6      2 3
    6     2 3    7      3 1
    7     3 1    8      3 2
    8     3 2    9      3 3
    9     3 3    0      0 0
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    29 April 2003
  Author:
    John Burkardt
  Parameters:
    Input, int M1, M2, the minimum and maximum entries.
    Input, int N, the number of components.
    Input/output, int *RANK, counts the elements.
    On first call, set RANK to 0.  Thereafter, the output value of RANK
    will indicate the order of the element returned.  When there are no
    more elements, RANK will be returned as 0.
    Input/output, int X[N], on input the previous tuple.
    On output, the next tuple.
*/
{
	const _3dtpipdt * const s_data = data;
	
	const register dim_typ m1 = s_data->a0;
	const register dim_typ m2 = s_data->a1;
	const register dim_typ n = s_data->a2;
	int * x = s_data->a3;
	dim_typ * rank = s_data->a4;

	
	int i;
	int j;
	
	if ( m2 < m1 )
	{
		*rank = 0;
		return NULL;
	}
	
	if ( *rank <= 0 )
	{
		for ( i = 0; i < n; ++i )
			x[i] = m1;
		*rank = 1;
	}
	else
	{
		++ *rank;
		i = n - 1;
	
		for ( ; ; )
		{
	
			if ( x[i] < m2 )
			{
				++ x[i];
				break;
			}
	
			x[i] = m1;
	
			if ( i == 0 )
			{
				*rank = 0;
				for ( j = 0; j < n; ++j)
					x[j] = m1;
				break;
			}
			-- i;
		}
	}
	
	return NULL;
}

#endif
