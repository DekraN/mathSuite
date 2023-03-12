#ifndef __DISABLEDEEP_VALUES

#include "dutils.h"


/******************************************************************************/
__MATHSUITE  void *   _student_noncentral_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = NoncentralStudentTDistribution [ df, lambda ]
      CDF [ dist, x ]
    Mathematica seems to have some difficulty computing this function
    to the desired number of digits.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 September 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *DF, double *LAMBDA, the parameters of the
    function.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdtpi3pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	int * df = s_data->a1;
	ityp * lambda = s_data->a2;
	ityp * x = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 30

    int df_vec[N_MAX] =
    {
        1,  2,  3,
        1,  2,  3,
        1,  2,  3,
        1,  2,  3,
        1,  2,  3,
        15, 20, 25,
        1,  2,  3,
        10, 10, 10,
        10, 10, 10,
        10, 10, 10
    };

    ityp fx_vec[N_MAX] =
    {
        0.8975836176504333E+00,
        0.9522670169E+00,
        0.9711655571887813E+00,
        0.8231218864E+00,
        0.9049021510E+00,
        0.9363471834E+00,
        0.7301025986E+00,
        0.8335594263E+00,
        0.8774010255E+00,
        0.5248571617E+00,
        0.6293856597E+00,
        0.6800271741E+00,
        0.20590131975E+00,
        0.2112148916E+00,
        0.2074730718E+00,
        0.9981130072E+00,
        0.9994873850E+00,
        0.9998391562E+00,
        0.168610566972E+00,
        0.16967950985E+00,
        0.1701041003E+00,
        0.9247683363E+00,
        0.7483139269E+00,
        0.4659802096E+00,
        0.9761872541E+00,
        0.8979689357E+00,
        0.7181904627E+00,
        0.9923658945E+00,
        0.9610341649E+00,
        0.8688007350E+00
    };

    ityp lambda_vec[N_MAX] =
    {
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        4.0E+00,
        4.0E+00,
        4.0E+00,
        7.0E+00,
        7.0E+00,
        7.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        2.0E+00,
        3.0E+00,
        4.0E+00,
        2.0E+00,
        3.0E+00,
        4.0E+00,
        2.0E+00,
        3.0E+00,
        4.0E+00
    };

    ityp x_vec[N_MAX] =
    {
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        3.00E+00,
        15.00E+00,
        15.00E+00,
        15.00E+00,
        0.05E+00,
        0.05E+00,
        0.05E+00,
        4.00E+00,
        4.00E+00,
        4.00E+00,
        5.00E+00,
        5.00E+00,
        5.00E+00,
        6.00E+00,
        6.00E+00,
        6.00E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *df = 0;
        *lambda = *x = *fx = 0.00;
    }
    else
    {
        *df = df_vec[*n_data-1];
        *lambda = lambda_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _normal_01_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = NormalDistribution [ 0, 1 ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 August 2004
  Author
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 17

    static ityp fx_vec[N_MAX] =
    {
        0.5000000000000000E+00,
        0.5398278372770290E+00,
        0.5792597094391030E+00,
        0.6179114221889526E+00,
        0.6554217416103242E+00,
        0.6914624612740131E+00,
        0.7257468822499270E+00,
        0.7580363477769270E+00,
        0.7881446014166033E+00,
        0.8159398746532405E+00,
        0.8413447460685429E+00,
        0.9331927987311419E+00,
        0.9772498680518208E+00,
        0.9937903346742239E+00,
        0.9986501019683699E+00,
        0.9997673709209645E+00,
        0.9999683287581669E+00
    };

    static ityp x_vec[N_MAX] =
    {
        0.0000000000000000E+00,
        0.1000000000000000E+00,
        0.2000000000000000E+00,
        0.3000000000000000E+00,
        0.4000000000000000E+00,
        0.5000000000000000E+00,
        0.6000000000000000E+00,
        0.7000000000000000E+00,
        0.8000000000000000E+00,
        0.9000000000000000E+00,
        0.1000000000000000E+01,
        0.1500000000000000E+01,
        0.2000000000000000E+01,
        0.2500000000000000E+01,
        0.3000000000000000E+01,
        0.3500000000000000E+01,
        0.4000000000000000E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _gamma_values ( void * data)
/******************************************************************************/
//
//  Purpose:
//    GAMMA_VALUES returns some values of the Gamma function.
//  Definition:
//    GAMMA(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) EXP(-T) dT
//  Recursion:
//    GAMMA(X+1) = X*GAMMA(X)
//  Restrictions:
//    0 < X ( a software restriction).
//  Special values:
//    GAMMA(0.5) = sqrt(M_PI)
//    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
//  Modified:
//    31 May 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, double *X, the argument of the function.
//    Output, double *FX, the value of the function.
//
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 18

    ityp fx_vec[N_MAX] =
    {
        4.590845E+00,     2.218160E+00,     1.489192E+00,     1.164230E+00,
        1.0000000000E+00, 0.9513507699E+00, 0.9181687424E+00, 0.8974706963E+00,
        0.8872638175E+00, 0.8862269255E+00, 0.8935153493E+00, 0.9086387329E+00,
        0.9313837710E+00, 0.9617658319E+00, 1.0000000000E+00, 3.6288000E+05,
        1.2164510E+17,    8.8417620E+30
    };
    ityp x_vec[N_MAX] =
    {
        0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00,
        1.0E+00,  1.1E+00,  1.2E+00,  1.3E+00,
        1.4E+00,  1.5E+00,  1.6E+00,  1.7E+00,
        1.8E+00,  1.9E+00,  2.0E+00, 10.0E+00,
        20.0E+00, 30.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.0E+00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }
    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _i4_fall_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_FALL_VALUES returns values of the integer falling factorial function.
  Discussion:
    The definition of the falling factorial function is
   (m)_n = (m)! / (m-n)!
            = ( m ) * ( m - 1 ) * ( m - 2 ) ... * ( m - n + 1 )
            = Gamma ( m + 1 ) / Gamma ( m - n + 1 )
    We assume 0 <= N <= M.
    In Mathematica, the function can be evaluated by:
      FactorialPower[m,n]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 December 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *M, N, the arguments of the function.
    Output, int *FMN, the value of the function.
*/
{
	dim_typ ** const a_data = data;

	dim_typ * n_data = a_data[0];
	dim_typ * m = a_data[1];
	dim_typ * n = a_data[2];
	dim_typ * fmn = a_data[3];
	
    # define N_MAX 15

    static unsigned fmn_vec[N_MAX] =
    {
        1, 5, 20, 60, 120,
        120, 0, 1, 10, 4000,
        90, 4896, 24, 912576, 0
    };

    static dim_typ m_vec[N_MAX] =
    {
        5, 5, 5, 5, 5,
        5, 5, 50, 10, 4000,
        10, 18, 4, 98, 1
    };

    static dim_typ n_vec[N_MAX] =
    {
        0, 1, 2, 3, 4,
        5, 6, 0, 1, 1,
        2, 3, 4, 3, 7
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *m = *n = *fmn = 0;
    else
    {
        *m = m_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fmn = fmn_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _bessel_j0_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BESSEL_J0_VALUES returns some values of the J0 Bessel function.
  Discussion:
    In Mathematica, the function can be evaluated by:
      BesselJ[0,x]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 21

    static ityp fx_vec[N_MAX] =
    {
        -0.1775967713143383E+00,
        -0.3971498098638474E+00,
        -0.2600519549019334E+00,
        0.2238907791412357E+00,
        0.7651976865579666E+00,
        0.1000000000000000E+01,
        0.7651976865579666E+00,
        0.2238907791412357E+00,
        -0.2600519549019334E+00,
        -0.3971498098638474E+00,
        -0.1775967713143383E+00,
        0.1506452572509969E+00,
        0.3000792705195556E+00,
        0.1716508071375539E+00,
        -0.9033361118287613E-01,
        -0.2459357644513483E+00,
        -0.1711903004071961E+00,
        0.4768931079683354E-01,
        0.2069261023770678E+00,
        0.1710734761104587E+00,
        -0.1422447282678077E-01
    };

    static ityp x_vec[N_MAX] =
    {
        -5.0E+00,
        -4.0E+00,
        -3.0E+00,
        -2.0E+00,
        -1.0E+00,
        0.0E+00,
        1.0E+00,
        2.0E+00,
        3.0E+00,
        4.0E+00,
        5.0E+00,
        6.0E+00,
        7.0E+00,
        8.0E+00,
        9.0E+00,
        10.0E+00,
        11.0E+00,
        12.0E+00,
        13.0E+00,
        14.0E+00,
        15.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _erfc_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    ERFC_VALUES returns some values of the ERFC or "complementary error" function.
  Discussion:
    The complementary error function is defined by:
      ERFC(X) = 1 - ( 2 / sqrt ( M_PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
    In Mathematica, the function can be evaluated by:
      Erfc[x]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 May 2007
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 21

    ityp fx_vec[N_MAX] =
    {
        1.000000000000000E+00,
        0.7772974107895215E+00,
        0.5716076449533315E+00,
        0.3961439091520741E+00,
        0.2578990352923395E+00,
        0.1572992070502851E+00,
        0.08968602177036462E+00,
        0.04771488023735119E+00,
        0.02365161665535599E+00,
        0.01090949836426929E+00,
        0.004677734981047266E+00,
        0.001862846297981891E+00,
        0.0006885138966450786E+00,
        0.0002360344165293492E+00,
        0.00007501319466545902E+00,
        0.00002209049699858544E+00,
        6.025761151762095E-06,
        1.521993362862285E-06,
        3.558629930076853E-07,
        7.700392745696413E-08,
        1.541725790028002E-08
    };

    ityp x_vec[N_MAX] =
    {
        0.0E+00,
        0.2E+00,
        0.4E+00,
        0.6E+00,
        0.8E+00,
        1.0E+00,
        1.2E+00,
        1.4E+00,
        1.6E+00,
        1.8E+00,
        2.0E+00,
        2.2E+00,
        2.4E+00,
        2.6E+00,
        2.8E+00,
        3.0E+00,
        3.2E+00,
        3.4E+00,
        3.6E+00,
        3.8E+00,
        4.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _h_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    H_POLYNOMIAL_VALUES: tabulated values of H(i,x).
  Discussion:
    H(i,x) is the physicist's Hermite polynomial of degree I.
    In Mathematica, the function can be evaluated by:
      HermiteH[n,x]
  Differential equation:
    Y'' - 2 X Y' + 2 N Y = 0;
  First terms:
      1
      2 X
      4 X^2     -  2
      8 X^3     - 12 X
     16 X^4     - 48 X^2     + 12
     32 X^5    - 160 X^3    + 120 X
     64 X^6    - 480 X^4    + 720 X^2    - 120
    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
  Recursion:
    H(0,X) = 1,
    H(1,X) = 2*X,
    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
  Norm:
    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
    = sqrt ( M_PI ) * 2^N * N!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 February 2012
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the polynomial.
    Output, double *X, the point where the polynomial is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 18

    static ityp fx_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+02,
        0.9800000000000000E+02,
        0.9400000000000000E+03,
        0.8812000000000000E+04,
        0.8060000000000000E+05,
        0.7178800000000000E+06,
        0.6211600000000000E+07,
        0.5206568000000000E+08,
        0.4212712000000000E+09,
        0.3275529760000000E+10,
        0.2432987360000000E+11,
        0.1712370812800000E+12,
        0.0000000000000000E+00,
        0.4100000000000000E+02,
        -0.8000000000000000E+01,
        0.3816000000000000E+04,
        0.3041200000000000E+07
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10, 11,
        12,  5,  5,
        5,  5,  5
    };

    static ityp x_vec[N_MAX] =
    {
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        0.0E+00,
        0.5E+00,
        1.0E+00,
        3.0E+00,
        1.0E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _he_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    HE_POLYNOMIAL_VALUES: tabulated values of He(i,x).
  Discussion:
    He(i,x) represents the probabilist's Hermite polynomial.
    In Mathematica, the function can be evaluated by:
      He(n,x) = HermiteH[n,x/Sqrt[2]] / Sqrt [ 2^n ]
  First terms:
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
  Recursion:
    He(0,X) = 1,
    He(1,X) = X,
    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
  Norm:
    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX
    = sqrt ( 2 * M_PI ) * N! * delta ( M, N )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 February 2012
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the polynomial.
    Output, double *X, the point where the polynomial is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 18

    static ityp fx_vec[N_MAX] =
    {
        1.000000000000000E+00,
        5.000000000000000E+00,
        24.00000000000000E+00,
        110.0000000000000E+00,
        478.0000000000000E+00,
        1950.000000000000E+00,
        7360.000000000000E+00,
        25100.00000000000E+00,
        73980.00000000000E+00,
        169100.0000000000E+00,
        179680.0000000000E+00,
        -792600.0000000000E+00,
        -5939480.000000000E+00,
        0.000000000000000E+00,
        6.281250000000000E+00,
        6.000000000000000E+00,
        18.00000000000000E+00,
        90150.00000000000E+00
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10, 11,
        12,  5,  5,
        5,  5,  5
    };

    static ityp x_vec[N_MAX] =
    {
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        0.0E+00,
        0.5E+00,
        1.0E+00,
        3.0E+00,
        1.0E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _hf_function_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    HF_FUNCTION_VALUES: tabulated values of Hf(i,x).
   Discussion:
    Hf(i,x) represents the Hermite function of "degree" I.
    In Mathematica, the function can be evaluated by:
      Hf(n,x) = HermiteH[n,x]
        * Exp [ -1/2 * x^2] / Sqrt [ 2^n * n! * Sqrt[Pi] ]
    The Hermite functions are orthonormal:
      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 February 2012
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the polynomial.
    Output, double *X, the point where the polynomial is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 23

    static ityp fx_vec[N_MAX] =
    {
        0.7511255444649425E+00,  0.0000000000000000E+00, -0.5311259660135985E+00,
        0.0000000000000000E+00,  0.4599685791773266E+00,  0.0000000000000000E+00,
        0.4555806720113325E+00,  0.6442883651134752E+00,  0.3221441825567376E+00,
        -0.2630296236233334E+00, -0.4649750762925110E+00, -0.5881521185179581E-01,
        0.3905052515434106E+00,  0.2631861423064045E+00, -0.2336911435996523E+00,
        -0.3582973361472840E+00,  0.6146344487883041E-01,  0.3678312067984882E+00,
        0.9131969309166278E-01,  0.4385750950032321E+00, -0.2624689527931006E-01,
        0.5138426125477819E+00,  0.9355563118061758E-01
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10, 11,
        12,  5,  5,
        5,  5
    };

    static ityp x_vec[N_MAX] =
    {
        0.0E+00, 0.0E+00, 0.0E+00,
        0.0E+00, 0.0E+00, 0.0E+00,
        1.0E+00, 1.0E+00, 1.0E+00,
        1.0E+00, 1.0E+00, 1.0E+00,
        1.0E+00, 1.0E+00, 1.0E+00,
        1.0E+00, 1.0E+00, 1.0E+00,
        1.0E+00, 0.5E+00, 2.0E+00,
        3.0E+00, 4.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _hep_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    HEP_VALUES returns values of the Hermite polynomials He(n,x).
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
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the function.
    Output, double *X, the point where the function is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 18

    static ityp fx_vec[N_MAX] =
    {
        1.000000000000000E+00,
        5.000000000000000E+00,
        24.00000000000000E+00,
        110.0000000000000E+00,
        478.0000000000000E+00,
        1950.000000000000E+00,
        7360.000000000000E+00,
        25100.00000000000E+00,
        73980.00000000000E+00,
        169100.0000000000E+00,
        179680.0000000000E+00,
        -792600.0000000000E+00,
        -5939480.000000000E+00,
        0.000000000000000E+00,
        6.281250000000000E+00,
        6.000000000000000E+00,
        18.00000000000000E+00,
        90150.00000000000E+00
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10, 11,
        12,  5,  5,
        5,  5,  5
    };

    static ityp x_vec[N_MAX] =
    {
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        5.0E+00,
        0.0E+00,
        0.5E+00,
        1.0E+00,
        3.0E+00,
        1.0E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _hypersphere_01_area_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_01_AREA_VALUES returns some areas of the unit hypersphere.
  Discussion:
    The unit hypersphere satisfies the equation:
      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
     M    Area
     2    2        * M_PI
     3 ( 4 /    ) * M_PI
     4 ( 2 /   1) * M_PI^2
     5 ( 8 /   3) * M_PI^2
     6 ( 1 /   1) * M_PI^3
     7 (16 /  15) * M_PI^3
     8 ( 1 /   3) * M_PI^4
     9 (32 / 105) * M_PI^4
    10 ( 1 /  12) * M_PI^5
    For the unit hypersphere, Area(M) = M * Volume(M)
    Sphere_Unit_Area ( M ) = 2 * M_PI^(M/2) / Gamma ( M / 2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 December 2013
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N_DATA.
    On input, if N_DATA is 0, the first test data is returned, and
    N_DATA is set to the index of the test data.  On each subsequent
    call, N_DATA is incremented and that test data is returned.  When
    there is no more test data, N_DATA is set to 0.
    Output, int *M, the spatial dimension.
    Output, double *AREA, the area of the unit hypersphere.
*/
{
	const _2pdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * m = s_data->a1;
	ityp * area = s_data->a2;
	
    # define NMAX 20

    ityp area_vec[NMAX] =
    {
        0.2000000000000000E+01,
        0.6283185307179586E+01,
        0.1256637061435917E+02,
        0.1973920880217872E+02,
        0.2631894506957162E+02,
        0.3100627668029982E+02,
        0.3307336179231981E+02,
        0.3246969701133415E+02,
        0.2968658012464836E+02,
        0.2550164039877345E+02,
        0.2072514267328890E+02,
        0.1602315322625507E+02,
        0.1183817381218268E+02,
        0.8389703410491089E+01,
        0.5721649212349567E+01,
        0.3765290085742291E+01,
        0.2396678817591364E+01,
        0.1478625959000308E+01,
        0.8858104195716824E+00,
        0.5161378278002812E+00
    };
    dim_typ m_vec[NMAX] =
    {
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20
    };

    if ( NMAX <= *n_data )
    {
        *n_data = *m = 0;
        *area = 0.00;
    }
    else
    {
        *m = m_vec[*n_data];
        *area = area_vec[*n_data];
        ++ *n_data;
    }

    return NULL;
    # undef NMAX
}

/******************************************************************************/
__MATHSUITE  void *   _hypersphere_01_volume_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERSPHERE_01_VOLUME_VALUES returns some volumes of the unit hypersphere.
  Discussion:
    The unit hypersphere satisfies the equation:
      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
     M  Volume
     1    1
     2    1        * M_PI
     3 ( 4 /   3) * M_PI
     4 ( 1 /   2) * M_PI^2
     5 ( 8 /  15) * M_PI^2
     6 ( 1 /   6) * M_PI^3
     7 (16 / 105) * M_PI^3
     8 ( 1 /  24) * M_PI^4
     9 (32 / 945) * M_PI^4
    10 ( 1 / 120) * M_PI^5
    For the unit hypersphere, Volume(M) = 2 * M_PI * Volume(M-2) / M
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2003
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N_DATA.
    On input, if N_DATA is 0, the first test data is returned, and
    N_DATA is set to the index of the test data.  On each subsequent
    call, N_DATA is incremented and that test data is returned.  When
    there is no more test data, N_DATA is set to 0.
    Output, int *M, the spatial dimension.
    Output, double *VOLUME, the volume.
*/
{
	const _2pdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * m = s_data->a1;
	ityp * volume = s_data->a2;
	
    # define N_MAX 20

    dim_typ m_vec[N_MAX] =
    {
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20
    };

    ityp volume_vec[N_MAX] =
    {
        0.2000000000000000E+01,
        0.3141592653589793E+01,
        0.4188790204786391E+01,
        0.4934802200544679E+01,
        0.5263789013914325E+01,
        0.5167712780049970E+01,
        0.4724765970331401E+01,
        0.4058712126416768E+01,
        0.3298508902738707E+01,
        0.2550164039877345E+01,
        0.1884103879389900E+01,
        0.1335262768854589E+01,
        0.9106287547832831E+00,
        0.5992645293207921E+00,
        0.3814432808233045E+00,
        0.2353306303588932E+00,
        0.1409811069171390E+00,
        0.8214588661112823E-01,
        0.4662160103008855E-01,
        0.2580689139001406E-01
    };

    if ( N_MAX <= *n_data )
    {
        *n_data = *m = 0;
        *volume = 0.00;
    }
    else
    {
        *m = m_vec[*n_data];
        *volume = volume_vec[*n_data];
        ++ *n_data;
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _i4_factorial2_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_FACTORIAL2_VALUES returns values of the double factorial function.
  Formula:
    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 ) (N even)
                    = Product ( N * (N-2) * (N-4) * ... * 1 ) (N odd)
    In Mathematica, the function can be evaluated by:
      n!!
  Example:
     N    N!!
     0     1
     1     1
     2     2
     3     3
     4     8
     5    15
     6    48
     7   105
     8   384
     9   945
    10  3840
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
    Daniel Zwillinger,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996, page 16.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the function.
    Output, int *FN, the value of the function.
*/
{
	const _2pdtpu * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	unsigned * fn = s_data->a2;
	
    # define N_MAX 16

    static unsigned fn_vec[N_MAX] =
    {
        1,
        1,
        2,
        3,
        8,
        15,
        48,
        105,
        384,
        945,
        3840,
        10395,
        46080,
        135135,
        645120,
        2027025
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,
        1,  2,  3,  4,  5,
        6,  7,  8,  9, 10,
        11, 12, 13, 14, 15
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *fn = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *fn = fn_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _i4_rise_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_RISE_VALUES returns values of the integer rising factorial function.
  Discussion:
    The integer rising factorial function is sometimes symbolized by (m)_n.
    The definition i
   (m)_n = (m-1+n)! / (m-1)!
            = ( m ) * ( m + 1 ) * ( m + 2 ) ... * ( m - 1 + n )
            = Gamma ( m + n ) / Gamma ( m )
    We assume 0 <= N <= M.
    In Mathematica, the function can be evaluated by:
      Pochhammer[m,n]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 December 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *M, N, the arguments of the function.
    Output, int *FMN, the value of the function.
*/
{
	const _3pdtpu * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * m = s_data->a1;
	dim_typ * n = s_data->a2;
	unsigned * fmn = s_data->a3;
	
    # define N_MAX 15

    static unsigned fmn_vec[N_MAX] =
    {
        1, 5, 30, 210, 1680,
        15120, 151200, 1, 10, 4000,
        110, 6840, 840, 970200, 5040
    };

    static dim_typ m_vec[N_MAX] =
    {
        5, 5, 5, 5, 5,
        5, 5, 50, 10, 4000,
        10, 18, 4, 98, 1
    };

    static dim_typ n_vec[N_MAX] =
    {
        0, 1, 2, 3, 4,
        5, 6, 0, 1, 1,
        2, 3, 4, 3, 7
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *m = *n = *fmn = 0;
    else
    {
        *m = m_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fmn = fmn_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _j_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    J_POLYNOMIAL_VALUES returns some values of the Jacobi polynomial.
  Discussion:
    In Mathematica, the function can be evaluated by:
      JacobiP[ n, a, b, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the degree of the polynomial.
    Output, double *A, *B, parameters of the function.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt4pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	ityp * x = s_data->a4;
	ityp * fx = s_data->a5;
	
    # define N_MAX 26

    static ityp a_vec[N_MAX] =
    {
        0.00, 0.00, 0.00, 0,
        0.00, 0.00, 1.00, 2,
        3.00, 4.00, 5.00, 0,
        0.00, 0.00, 0.00, 0,
        0.00, 0.00, 0.00, 0,
        0.00, 0.00, 0.00, 0,
        0.00, 0.00
    };

    static ityp b_vec[N_MAX] =
    {
        1.00, 1.00, 1.00, 1.00,
        1.00, 1.00, 1.00, 1.00,
        1.00, 1.00, 1.00, 2.00,
        3.00, 4.0, 5.00, 1.00,
        1.00, 1.00, 1.00, 1.00,
        1.00, 1.00, 1.00, 1.00,
        1.00, 1.00
    };

    static ityp fx_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2500000000000000E+00,
        -0.3750000000000000E+00,
        -0.4843750000000000E+00,
        -0.1328125000000000E+00,
        0.2753906250000000E+00,
        -0.1640625000000000E+00,
        -0.1174804687500000E+01,
        -0.2361328125000000E+01,
        -0.2616210937500000E+01,
        0.1171875000000000E+00,
        0.4218750000000000E+00,
        0.5048828125000000E+00,
        0.5097656250000000E+00,
        0.4306640625000000E+00,
        -0.6000000000000000E+01,
        0.3862000000000000E-01,
        0.8118400000000000E+00,
        0.3666000000000000E-01,
        -0.4851200000000000E+00,
        -0.3125000000000000E+00,
        0.1891200000000000E+00,
        0.4023400000000000E+00,
        0.1216000000000000E-01,
        -0.4396200000000000E+00,
        0.1000000000000000E+01
    };

    static dim_typ n_vec[N_MAX] =
    {
        0, 1, 2, 3,
        4, 5, 5, 5,
        5, 5, 5, 5,
        5, 5, 5, 5,
        5, 5, 5, 5,
        5, 5, 5, 5,
        5, 5
    };

    static ityp x_vec[N_MAX] =
    {
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        -1.0E+00,
        -0.8E+00,
        -0.6E+00,
        -0.4E+00,
        -0.2E+00,
        0.0E+00,
        0.2E+00,
        0.4E+00,
        0.6E+00,
        0.8E+00,
        1.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *a = *b = *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _l_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    L_POLYNOMIAL_VALUES returns some values of the Laguerre polynomial L(n,x).
  Discussion:
    In Mathematica, the function can be evaluated by:
      LaguerreL[n,x]
  Differential equation:
    X * Y'' + (1-X) * Y' + N * Y = 0;
  First terms:
      1
     -X    +  1
 (  X^2 -  4 X     +  2 ) / 2
 ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
 (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
 ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
 (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
 ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
      + 52920 X^2 - 35280 X + 5040 ) / 5040
  Recursion:
    L(0,X) = 1,
    L(1,X) = 1-X,
    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
  Orthogonality:
    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX
    = 0 if N /= M
    = 1 if N == M
  Special values:
    L(N,0) = 1.
  Relations:
    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the polynomial.
    Output, double *X, the point where the polynomial is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 17

    static ityp fx_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.0000000000000000E+00,
        -0.5000000000000000E+00,
        -0.6666666666666667E+00,
        -0.6250000000000000E+00,
        -0.4666666666666667E+00,
        -0.2569444444444444E+00,
        -0.4047619047619048E-01,
        0.1539930555555556E+00,
        0.3097442680776014E+00,
        0.4189459325396825E+00,
        0.4801341790925124E+00,
        0.4962122235082305E+00,
        -0.4455729166666667E+00,
        0.8500000000000000E+00,
        -0.3166666666666667E+01,
        0.3433333333333333E+02
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10, 11,
        12,  5,  5,
        5,  5
    };

    static ityp x_vec[N_MAX] =
    {
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        0.5E+00,
        3.0E+00,
        5.0E+00,
        1.0E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
# undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _lf_function_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LF_FUNCTION_VALUES: some values of the Laguerre function Lf(n,alpha,x).
  Discussion:
    In Mathematica, the function can be evaluated by:
      LaguerreL[n,a,x]
    The functions satisfy the following differential equation:
      X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0;
    Function values can be generated by the recursion:
      Lf(0,ALPHA,X) = 1
      Lf(1,ALPHA,X) = 1+ALPHA-X
      Lf(N,ALPHA,X) = ( (2*N-1+ALPHA-X) * Lf(N-1,ALPHA,X)
                     - (N-1+ALPHA) * Lf(N-2,ALPHA,X) ) / N
    The parameter ALPHA is required to be greater than -1.
    For ALPHA = 0, the generalized Laguerre function Lf(N,ALPHA,X)
    is equal to the Laguerre polynomial L(N,X).
    For ALPHA integral, the generalized Laguerre function
    Lf(N,ALPHA,X) equals the associated Laguerre polynomial Lm(N,ALPHA,X).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the function.
    Output, double *A, the parameter.
    Output, double *X, the point where the function is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pitpdtpit * const s_data = data; 
	
	dim_typ * n_data = s_data->a0;
	ityp * a = s_data->a1;
	ityp * x = s_data->a2;
	dim_typ * n = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 20

    static ityp a_vec[N_MAX] =
    {
        0.00E+00,
        0.25E+00,
        0.50E+00,
        0.75E+00,
        1.50E+00,
        2.50E+00,
        5.00E+00,
        1.20E+00,
        1.20E+00,
        1.20E+00,
        1.20E+00,
        1.20E+00,
        1.20E+00,
        5.20E+00,
        5.20E+00,
        5.20E+00,
        5.20E+00,
        5.20E+00,
        5.20E+00,
        5.20E+00
    };

    static ityp fx_vec[N_MAX] =
    {
        0.3726399739583333E-01,
        0.3494791666666667E+00,
        0.8710042317708333E+00,
        0.1672395833333333E+01,
        0.6657625325520833E+01,
        0.2395726725260417E+02,
        0.2031344319661458E+03,
        0.1284193996800000E+02,
        0.5359924801587302E+01,
        0.9204589064126984E+00,
        -0.1341585114857143E+01,
        -0.2119726307555556E+01,
        -0.1959193658349206E+01,
        0.1000000000000000E+01,
        0.5450000000000000E+01,
        0.1720125000000000E+02,
        0.4110393750000000E+02,
        0.8239745859375000E+02,
        0.1460179186171875E+03,
        0.2359204608298828E+03
    };

    static dim_typ n_vec[N_MAX] =
    {
        5,
        5,
        5,
        5,
        5,
        5,
        5,
        8,
        8,
        8,
        8,
        8,
        8,
        0,
        1,
        2,
        3,
        4,
        5,
        6
    };

    static ityp x_vec[N_MAX] =
    {
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.00E+00,
        0.20E+00,
        0.40E+00,
        0.60E+00,
        0.80E+00,
        1.00E+00,
        0.75E+00,
        0.75E+00,
        0.75E+00,
        0.75E+00,
        0.75E+00,
        0.75E+00,
        0.75E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
    *n_data = *n = 0;
    *a = *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *a = a_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _pmn_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    PMN_POLYNOMIAL_VALUES: normalized Legendre polynomial Pmn(n,m,x).
  Discussion:
    In Mathematica, the unnormalized function can be evaluated by:
      LegendreP [ n, m, x ]
    The function is normalized by dividing by
      sqrt ( 2 * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
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
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, int *M, double *X,
    the arguments of the function.
    Output, double *FX, the value of the function.
*/
{
	const _3pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	dim_typ * m = s_data->a2;
	ityp * x = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 21

    static ityp fx_vec[N_MAX] =
    {
        0.7071067811865475E+00,
        0.6123724356957945E+00,
        -0.7500000000000000E+00,
        -0.1976423537605237E+00,
        -0.8385254915624211E+00,
        0.7261843774138907E+00,
        -0.8184875533567997E+00,
        -0.1753901900050285E+00,
        0.9606516343087123E+00,
        -0.6792832849776299E+00,
        -0.6131941618102092E+00,
        0.6418623720763665E+00,
        0.4716705890038619E+00,
        -0.1018924927466445E+01,
        0.6239615396237876E+00,
        0.2107022704608181E+00,
        0.8256314721961969E+00,
        -0.3982651281554632E+00,
        -0.7040399320721435E+00,
        0.1034723155272289E+01,
        -0.5667412129155530E+00
    };

    static dim_typ m_vec[N_MAX] =
    {
        0, 0, 1, 0,
        1, 2, 0, 1,
        2, 3, 0, 1,
        2, 3, 4, 0,
        1, 2, 3, 4,
        5
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,  1,  1,  2,
        2,  2,  3,  3,
        3,  3,  4,  4,
        4,  4,  4,  5,
        5,  5,  5,  5,
        5
    };

    static ityp x_vec[N_MAX] =
    {
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = *m = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *m = m_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _pmns_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    PMNS_POLYNOMIAL_VALUES: sphere-normalized Legendre polynomial Pmns(n,m,x).
  Discussion:
    In Mathematica, the unnormalized function can be evaluated by:
      LegendreP [ n, m, x ]
    The function is normalized for the sphere by dividing by
      sqrt ( 4 * M_PI * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
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
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, int *M, double *X,
    the arguments of the function.
    Output, double *FX, the value of the function.
*/
{
	const _3pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	dim_typ * m = s_data->a2;
	ityp * x = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 21

    static ityp fx_vec[N_MAX] =
    {
        0.2820947917738781,
        0.2443012559514600,
        -0.2992067103010745,
        -0.07884789131313000,
        -0.3345232717786446,
        0.2897056515173922,
        -0.3265292910163510,
        -0.06997056236064664,
        0.3832445536624809,
        -0.2709948227475519,
        -0.2446290772414100,
        0.2560660384200185,
        0.1881693403754876,
        -0.4064922341213279,
        0.2489246395003027,
        0.08405804426339821,
        0.3293793022891428,
        -0.1588847984307093,
        -0.2808712959945307,
        0.4127948151484925,
        -0.2260970318780046
    };

    static dim_typ m_vec[N_MAX] =
    {
        0, 0, 1, 0,
        1, 2, 0, 1,
        2, 3, 0, 1,
        2, 3, 4, 0,
        1, 2, 3, 4,
        5
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,  1,  1,  2,
        2,  2,  3,  3,
        3,  3,  4,  4,
        4,  4,  4,  5,
        5,  5,  5,  5,
        5
    };

    static ityp x_vec[N_MAX] =
    {
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50,
        0.50
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = *m = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *m = m_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _lobatto_polynomial_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOBATTO_POLYNOMIAL_VALUES: values of the completed Lobatto polynomials.
  Discussion:
    In Mathematica, the function can be evaluated by:
      n * LegendreP [ n - 1, x ] - n * x * LegendreP [ n, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 May 2013
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0
    before the first call.  On each call, the routine increments N_DATA by 1,
    and returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the function.
    Output, double *X, the point where the function is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 31

    static ityp fx_vec[N_MAX] =
    {
        0.9375000000000000,
        0.7031250000000000,
        -0.9667968750000000,
        -1.501464843750000,
        0.3639221191406250,
        2.001914978027344,
        0.6597948074340820,
        -1.934441328048706,
        -1.769941113889217,
        1.215243665501475,
        0.000000000000000,
        0.8692500000000000,
        1.188000000000000,
        1.109250000000000,
        0.7680000000000000,
        0.2812500000000000,
        -0.2520000000000000,
        -0.7507500000000000,
        -1.152000000000000,
        -1.410750000000000,
        -1.500000000000000,
        -1.410750000000000,
        -1.152000000000000,
        -0.7507500000000000,
        -0.2520000000000000,
        0.2812500000000000,
        0.7680000000000000,
        1.109250000000000,
        1.188000000000000,
        0.8692500000000000,
        0.000000000000000
    };

    static dim_typ n_vec[N_MAX] =
    {
        1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10,  3,
        3,  3,  3,
        3,  3,  3,
        3,  3,  3,
        3,  3,  3,
        3,  3,  3,
        3,  3,  3,
        3,  3
    };

    static ityp x_vec[N_MAX] =
    {
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        -1.00,
        -0.90,
        -0.80,
        -0.70,
        -0.60,
        -0.50,
        -0.40,
        -0.30,
        -0.20,
        -0.10,
        0.00,
        0.10,
        0.20,
        0.30,
        0.40,
        0.50,
        0.60,
        0.70,
        0.80,
        0.90,
        1.00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n  = n_vec[*n_data-1];
        *x  = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _agm_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    AGM_VALUES returns some values of the AGM.
  Discussion:
    The AGM is defined for nonnegative A and B.
    The AGM of numbers A and B is defined by setting
      A(0) = A,
      B(0) = B
      A(N+1) = ( A(N) + B(N) ) / 2
      B(N+1) = sqrt ( A(N) * B(N) )
    The two sequences both converge to AGM(A,B).
    In Mathematica, the AGM can be evaluated by
      ArithmeticGeometricMean [ a, b ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 February 2008
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *A, *B, the argument ofs the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt3pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 15

    ityp a_vec[N_MAX] =
    {
        22.00,
        83.00,
        42.00,
        26.00,
        4.00,
        6.00,
        40.00,
        80.00,
        90.00,
        9.00,
        53.00,
        1.00,
        1.00,
        1.00,
        1.50
    };
    ityp b_vec[N_MAX] =
    {
        96.00,
        56.00,
        7.00,
        11.00,
        63.00,
        45.00,
        75.00,
        0.00,
        35.00,
        1.00,
        53.00,
        2.00,
        4.00,
        8.00,
        8.00
    };
    ityp fx_vec[N_MAX] =
    {
        52.274641198704240049,
        68.836530059858524345,
        20.659301196734009322,
        17.696854873743648823,
        23.867049721753300163,
        20.717015982805991662,
        56.127842255616681863,
        0.000000000000000000,
        59.269565081229636528,
        3.9362355036495554780,
        53.000000000000000000,
        1.4567910310469068692,
        2.2430285802876025701,
        3.6157561775973627487,
        4.0816924080221632670
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *a = *b = *fx = 0.00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _bernoulli_number_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNOULLI_NUMBER_VALUES returns some values of the Bernoulli numbers.
  Discussion:
    The Bernoulli numbers are rational.
    If we define the sum of the M-th powers of the first N integers as:
      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
    and let C(I,J) be the combinatorial coefficient:
      C(I,J) = I! / ( ( I - J )! * J! )
    then the Bernoulli numbers B(J) satisfy:
      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)**(M+1-J)
    In Mathematica, the function can be evaluated by:
      BernoulliB[n]
  First values:
   B0  1                   =         1.00000000000
   B1 -1/2                 =        -0.50000000000
   B2  1/6                 =         1.66666666666
   B3  0                   =         0
   B4 -1/30                =        -0.03333333333
   B5  0                   =         0
   B6  1/42                =         0.02380952380
   B7  0                   =         0
   B8 -1/30                =        -0.03333333333
   B9  0                   =         0
  B10  5/66                =         0.07575757575
  B11  0                   =         0
  B12 -691/2730            =        -0.25311355311
  B13  0                   =         0
  B14  7/6                 =         1.16666666666
  B15  0                   =         0
  B16 -3617/510            =        -7.09215686274
  B17  0                   =         0
  B18  43867/798           =        54.97117794486
  B19  0                   =         0
  B20 -174611/330          =      -529.12424242424
  B21  0                   =         0
  B22  854,513/138         =      6192.123
  B23  0                   =         0
  B24 -236364091/2730      =    -86580.257
  B25  0                   =         0
  B26  8553103/6           =   1425517.16666
  B27  0                   =         0
  B28 -23749461029/870     = -27298231.0678
  B29  0                   =         0
  B30  8615841276005/14322 = 601580873.901
  Recursion:
    With C(N+1,K) denoting the standard binomial coefficient,
    B(0) = 1.0
    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
  Special Values:
    Except for B(1), all Bernoulli numbers of odd index are 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the Bernoulli number.
    Output, double *C, the value of the Bernoulli number.
*/
{
	const _2pdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * c = s_data->a2;
	
    # define N_MAX 10

    ityp c_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        -0.5000000000000000E+00,
        0.1666666666666667E+00,
        0.0000000000000000E+00,
        -0.3333333333333333E-01,
        -0.2380952380952380E-01,
        -0.3333333333333333E-01,
        0.7575757575757575E-01,
        -0.5291242424242424E+03,
        0.6015808739006424E+09
    };

    dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,  3,  4, 6,  8, 10, 20, 30
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *c = 0.0E+00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _bernstein_poly_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNSTEIN_POLY_VALUES returns some values of the Bernstein polynomials.
  Discussion:
    The Bernstein polynomials are assumed to be based on [0,1].
    The formula for the Bernstein polynomials is
      B(N,I,X) = [N!/(I(N-I)!)] * (1-X)**(N-I) * X^I
    In Mathematica, the function can be evaluated by:
      Binomial[n,i] * (1-x)^(n-i) * x^i
  First values:
    B(0,0,X) = 1
    B(1,0,X) =      1-X
    B(1,1,X) =                X
    B(2,0,X) =  (1-X)^2
    B(2,1,X) = 2 * (1-X)   * X
    B(2,2,X) =               X^2
    B(3,0,X) =  (1-X)^3
    B(3,1,X) = 3 * (1-X)^2 * X
    B(3,2,X) = 3 * (1-X)   * X^2
    B(3,3,X) =               X^3
    B(4,0,X) =  (1-X)^4
    B(4,1,X) = 4 * (1-X)^3 * X
    B(4,2,X) = 6 * (1-X)^2 * X^2
    B(4,3,X) = 4 * (1-X)   * X^3
    B(4,4,X) =               X^4
  Special values:
    B(N,I,X) has a unique maximum value at X = I/N.
    B(N,I,X) has an I-fold zero at 0 and and N-I fold zero at 1.
    B(N,I,1/2) = C(N,K) / 2^N
    For a fixed X and N, the polynomials add up to 1:
      Sum ( 0 <= I <= N ) B(N,I,X) = 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 August 2004
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the degree of the polynomial.
    Output, int *K, the index of the polynomial.
    Output, double *X, the argument of the polynomial.
    Output, double *B, the value of the polynomial B(N,K,X).
*/
{
	const _3pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	dim_typ * k = s_data->a2;
	ityp * x = s_data->a3;
	ityp * b = s_data->a4;
	
    # define N_MAX 15

    ityp b_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.7500000000000000E+00,
        0.2500000000000000E+00,
        0.5625000000000000E+00,
        0.3750000000000000E+00,
        0.6250000000000000E-01,
        0.4218750000000000E+00,
        0.4218750000000000E+00,
        0.1406250000000000E+00,
        0.1562500000000000E-01,
        0.3164062500000000E+00,
        0.4218750000000000E+00,
        0.2109375000000000E+00,
        0.4687500000000000E-01,
        0.3906250000000000E-02
    };

    dim_typ k_vec[N_MAX] =
    {
        0,
        0, 1,
        0, 1, 2,
        0, 1, 2, 3,
        0, 1, 2, 3, 4
    };

    dim_typ n_vec[N_MAX] =
    {
        0,
        1, 1,
        2, 2, 2,
        3, 3, 3, 3,
        4, 4, 4, 4, 4
    };

    ityp x_vec[N_MAX] =
    {
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00,
        0.25E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = *k = 0;
        *x = *b = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *k = k_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *b = b_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}


/******************************************************************************/
__MATHSUITE  void *   _beta_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BETA_VALUES returns some values of the Beta function.
  Discussion:
    Beta(X,Y) = ( Gamma(X) * Gamma(Y) ) / Gamma(X+Y)
    Both X and Y must be greater than 0.
    In Mathematica, the function can be evaluated by:
      Beta[X,Y]
  Properties:
    Beta(X,Y) = Beta(Y,X).
    Beta(X,Y) = Integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
    Beta(X,Y) = Gamma(X) * Gamma(Y) / Gamma(X+Y)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, Y, the arguments of the function.
    Output, double *FXY, the value of the function.
*/
{
	const pdt3pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * fxy = s_data->a3;
	
    # define N_MAX 17

    ityp b_vec[N_MAX] =
    {
        0.5000000000000000E+01,
        0.2500000000000000E+01,
        0.1666666666666667E+01,
        0.1250000000000000E+01,
        0.5000000000000000E+01,
        0.2500000000000000E+01,
        0.1000000000000000E+01,
        0.1666666666666667E+00,
        0.3333333333333333E-01,
        0.7142857142857143E-02,
        0.1587301587301587E-02,
        0.2380952380952381E-01,
        0.5952380952380952E-02,
        0.1984126984126984E-02,
        0.7936507936507937E-03,
        0.3607503607503608E-03,
        0.8325008325008325E-04
    };

    ityp x_vec[N_MAX] =
    {
        0.2E+00,
        0.4E+00,
        0.6E+00,
        0.8E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        2.0E+00,
        3.0E+00,
        4.0E+00,
        5.0E+00,
        6.0E+00,
        6.0E+00,
        6.0E+00,
        6.0E+00,
        6.0E+00,
        7.0E+00
    };

    ityp y_vec[N_MAX] =
    {
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        0.2E+00,
        0.4E+00,
        1.0E+00,
        2.0E+00,
        3.0E+00,
        4.0E+00,
        5.0E+00,
        2.0E+00,
        3.0E+00,
        4.0E+00,
        5.0E+00,
        6.0E+00,
        7.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *y = *fxy = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *y = y_vec[*n_data-1];
        *fxy = b_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _catalan_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    CATALAN_VALUES returns some values of the Catalan numbers.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Binomial[2*n,n] / ( n + 1 )
  First values:
     C(0)     1
     C(1)     1
     C(2)     2
     C(3)     5
     C(4)    14
     C(5)    42
     C(6)   132
     C(7)   429
     C(8)  1430
     C(9)  4862
    C(10) 16796
  Formula:
    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
         = 1 / (N+1) * COMB ( 2N, N )
         = 1 / (2N+1) * COMB ( 2N+1, N+1).
  Recursion:
    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
  Discussion:
    The Catalan number C(N) counts:
    1) the number of binary trees on N vertices;
    2) the number of ordered trees on N+1 vertices;
    3) the number of full binary trees on 2N+1 vertices;
    4) the number of well formed sequences of 2N parentheses;
    5) the number of ways 2N ballots can be counted, in order,
       with N positive and N negative, so that the running sum
       is never negative;
    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
    7) the number of monotone functions from [1..N} to [1..N} which
       satisfy f(i) <= i for all i;
    8) the number of ways to triangulate a polygon with N+2 vertices.
  Example:
    N = 3
 ()()()
 ()(())
 (()())
 (())()
 ((()))
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 February 2003
  Author
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the Catalan number.
    Output, int *C, the value of the Catalan number.
*/
{
	dim_typ ** const a_data = data; 
	dim_typ * n_data = a_data[0];
	dim_typ * n = a_data[1];
	dim_typ *c = a_data[2];
	
    # define N_MAX 11

    dim_typ c_vec[N_MAX] =
    {
        1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796
    };

    dim_typ n_vec[N_MAX] =
    {
 		0 ,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10
    };


    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _cheby_t_poly_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBY_T_POLY_VALUES returns values of Chebyshev polynomials T(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the function.
    Output, double *X, the point where the function is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 13

    ityp fx_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.8000000000000000E+00,
        0.2800000000000000E+00,
        -0.3520000000000000E+00,
        -0.8432000000000000E+00,
        -0.9971200000000000E+00,
        -0.7521920000000000E+00,
        -0.2063872000000000E+00,
        0.4219724800000000E+00,
        0.8815431680000000E+00,
        0.9884965888000000E+00,
        0.7000513740800000E+00,
        0.1315856097280000E+00
    };

    dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10, 11,
        12
    };

    ityp x_vec[N_MAX] =
    {
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _cheby_u_poly_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBY_U_POLY_VALUES returns values of Chebyshev polynomials U(n,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 April 2012
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the function.
    Output, double *X, the point where the function is evaluated.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 13

    ityp fx_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1600000000000000E+01,
        0.1560000000000000E+01,
        0.8960000000000000E+00,
        -0.1264000000000000E+00,
        -0.1098240000000000E+01,
        -0.1630784000000000E+01,
        -0.1511014400000000E+01,
        -0.7868390400000000E+00,
        0.2520719360000000E+00,
        0.1190154137600000E+01,
        0.1652174684160000E+01,
        0.1453325357056000E+01
    };

    dim_typ n_vec[N_MAX] =
    {
    0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10, 11,
        12
    };

    ityp x_vec[N_MAX] =
    {
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00,
        0.8E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _collatz_count_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    COLLATZ_COUNT_VALUES returns some values of the Collatz count function.
  Discussion:
    The rules for generation of the Collatz sequence are recursive.
    If T is the current entry of the sequence, (T is
    assumed to be a positive integer), then the next
    entry, U is determined as follows:
      if T is 1 (or less)
        terminate the sequence;
      else if T is even
        U = T/2.
      else (if T is odd and not 1)
        U = 3*T+1;
    The Collatz count is the length of the Collatz sequence for a given
    starting value.  By convention, we include the initial value in the
    count, so the minimum value of the count is 1.
     N  Sequence                                                 Count
     1                                                               1
     2   1                                                           2
     3  10,  5, 16,  8,  4,  2,  1                                   8
     4   2   1                                                       3
     5  16,  8,  4,  2,  1                                           6
     6   3, 10,  5, 16,  8,  4,  2,  1                               9
     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
     8   4,  2,  1                                                   4
     9  28, 14,  7, ...                                             20
    10   5, 16,  8,  4,  2,  1                                       7
    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 March 2006
  Author:
    John Burkardt
  Reference:
    Eric Weisstein,
    "The Collatz Problem",
    CRC Concise Encyclopedia of Mathematics,
    CRC 1998.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the initial value of a Collatz sequence.
    Output, int *COUNT, the length of the Collatz sequence starting
    with N.
*/
{
	dim_typ ** const a_data = data; 
	dim_typ * n_data = a_data[0];
	dim_typ * n = a_data[1];
	dim_typ * count = a_data[2];
	
    # define N_MAX 20

    dim_typ count_vec[N_MAX] =
    {
        1,   2,   8,   3,   6,   9,   17,   4,  20,   7,
        112,  25,  26,  27,  17,  28,  111,  18,  83,  29
    };
    dim_typ n_vec[N_MAX] =
    {
        1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
        27,  50, 100, 200, 300, 400, 500, 600, 700, 800
    };


    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *count = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *count = count_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _cos_power_int_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    COS_POWER_INT_VALUES returns some values of the cosine power integral.
  Discussion:
    The function has the form
      COS_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( cos(T) )^N dt
    In Mathematica, the function can be evaluated by:
      Integrate [ ( Cos[x] )^n, { x, a, b } ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2012
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *A, B, the limits of integration.
    Output, int *N, the power.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pitpdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	dim_typ * n = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 11

    static ityp a_vec[N_MAX] =
    {
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00
    };

    static ityp b_vec[N_MAX] =
    {
        M_PI,
        M_PI,
        M_PI,
        M_PI,
        M_PI,
        M_PI,
        M_PI,
        M_PI,
        M_PI,
        M_PI,
        M_PI
    };

    static ityp fx_vec[N_MAX] =
    {
        M_PI,
        0.0,
        1.570796326794897,
        0.0,
        1.178097245096172,
        0.0,
        0.9817477042468104,
        0.0,
        0.8590292412159591,
        0.0,
        0.7731263170943632
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *a = *b = *fx = 0.00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _euler_number_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    EULER_NUMBER_VALUES returns some values of the Euler numbers.
  Discussion:
    In Mathematica, the function can be evaluated by:
      EulerE[n]
    These numbers rapidly get too big to store in an ordinary integer!
    The terms of odd index are 0.
    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
  First terms:
    E0  = 1
    E1  = 0;
    E2  = -1
    E3  = 0;
    E4  = 5
    E5  = 0;
    E6  = -61
    E7  = 0;
    E8  = 1385
    E9  = 0;
    E10 = -50521
    E11 = 0;
    E12 = 2702765
    E13 = 0;
    E14 = -199360981
    E15 = 0;
    E16 = 19391512145
    E17 = 0;
    E18 = -2404879675441
    E19 = 0;
    E20 = 370371188237525
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 February 2015
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order of the Euler number.
    Output, int *C, the value of the Euler number.
*/
{
	const _2pdtpi * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	int * c = s_data->a2;
	
    # define N_MAX 8

    int c_vec[N_MAX] =
    {
        1, 0, -1, 5, -61, 1385, -50521, 2702765
    };

    dim_typ n_vec[N_MAX] =
    {
        0, 1, 2, 4, 6, 8, 10, 12
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _gegenbauer_poly_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEGENBAUER_POLY_VALUES returns some values of the Gegenbauer polynomials.
  Discussion:
    The Gegenbauer polynomials are also known as the "spherical
    polynomials" or "ultraspherical polynomials".
    In Mathematica, the function can be evaluated by:
      GegenbauerC[n,m,x]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the order parameter of the function.
    Output, double *A, the real parameter of the function.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pitpdtpit * const s_data = data; 
	
	dim_typ * n_data = s_data->a0;
	ityp * a = s_data->a1;
	ityp * x = s_data->a2;
	dim_typ * n = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 38

    ityp a_vec[N_MAX] =
    {
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.5E+00,
        0.0E+00,
        1.0E+00,
        2.0E+00,
        3.0E+00,
        4.0E+00,
        5.0E+00,
        6.0E+00,
        7.0E+00,
        8.0E+00,
        9.0E+00,
        10.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00
    };

    ityp fx_vec[N_MAX] =
    {
        1.0000000000E+00,
        0.2000000000E+00,
        -0.4400000000E+00,
        -0.2800000000E+00,
        0.2320000000E+00,
        0.3075200000E+00,
        -0.0805760000E+00,
        -0.2935168000E+00,
        -0.0395648000E+00,
        0.2459712000E+00,
        0.1290720256E+00,
        0.0000000000E+00,
        -0.3600000000E+00,
        -0.0800000000E+00,
        0.8400000000E+00,
        2.4000000000E+00,
        4.6000000000E+00,
        7.4400000000E+00,
        10.9200000000E+00,
        15.0400000000E+00,
        19.8000000000E+00,
        25.2000000000E+00,
        -9.0000000000E+00,
        -0.1612800000E+00,
        -6.6729600000E+00,
        -8.3750400000E+00,
        -5.5267200000E+00,
        0.0000000000E+00,
        5.5267200000E+00,
        8.3750400000E+00,
        6.6729600000E+00,
        0.1612800000E+00,
        -9.0000000000E+00,
        -15.4252800000E+00,
        -9.6969600000E+00,
        22.4409600000E+00,
        100.8892800000E+00,
        252.0000000000E+00
    };

    dim_typ n_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10,  2,
        2,  2,  2,
        2,  2,  2,
        2,  2,  2,
        2,  5,  5,
        5,  5,  5,
        5,  5,  5,
        5,  5,  5,
        5,  5,  5,
        5,  5
    };

    ityp x_vec[N_MAX] =
    {
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.20E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        0.40E+00,
        -0.50E+00,
        -0.40E+00,
        -0.30E+00,
        -0.20E+00,
        -0.10E+00,
        0.00E+00,
        0.10E+00,
        0.20E+00,
        0.30E+00,
        0.40E+00,
        0.50E+00,
        0.60E+00,
        0.70E+00,
        0.80E+00,
        0.90E+00,
        1.00E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *a = *x = *fx = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *a = a_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _gud_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    GUD_VALUES returns some values of the Gudermannian function.
  Discussion:
    The Gudermannian function relates the hyperbolic and trigonomentric
    functions.  For any argument X, there is a corresponding value
    GD so that
      SINH(X) = TAN(GD).
    This value GD is called the Gudermannian of X and symbolized
    GD(X).  The inverse Gudermannian function is given as input a value
    GD and computes the corresponding value X.
    GD(X) = 2 * arctan ( exp ( X ) ) - M_PI / 2
    In Mathematica, the function can be evaluated by:
      2 * Atan[Exp[x]] - Pi/2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 August 2004
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 13

    ityp fx_vec[N_MAX] =
    {
        -0.1301760336046015E+01,
        -0.8657694832396586E+00,
        0.0000000000000000E+00,
        0.9983374879348662E-01,
        0.1986798470079397E+00,
        0.4803810791337294E+00,
        0.8657694832396586E+00,
        0.1131728345250509E+01,
        0.1301760336046015E+01,
        0.1406993568936154E+01,
        0.1471304341117193E+01,
        0.1510419907545700E+01,
        0.1534169144334733E+01
    };

    ityp x_vec[N_MAX] =
    {
        -2.0E+00,
        -1.0E+00,
        0.0E+00,
        0.1E+00,
        0.2E+00,
        0.5E+00,
        1.0E+00,
        1.5E+00,
        2.0E+00,
        2.5E+00,
        3.0E+00,
        3.5E+00,
        4.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _hyper_2f1_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric function 2F1.
  Discussion:
    In Mathematica, the function can be evaluated by:
      fx = Hypergeometric2F1 [ a, b, c, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 September 2007
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
    Shanjie Zhang, Jianming Jin,
    Computation of Special Functions,
    Wiley, 1996,
    ISBN: 0-471-11963-6,
    LC: QA351.C45
    Daniel Zwillinger,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition, CRC Press, 1996, pages 651-652.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *A, *B, *C, *X, the parameters of the function.
    Output, double *FX, the value of the function.
*/
{
	const _4pitpdtpit * const s_data = data;
	
	ityp * a = s_data->a0;
	ityp * b = s_data->a1;
	ityp * c = s_data->a2;
	ityp * x = s_data->a3;
	dim_typ * n_data = s_data->a4;
	ityp * fx = s_data->a5;
	
    # define N_MAX 24

    ityp a_vec[N_MAX] =
    {
        -2.50,
        -0.50,
        0.50,
        2.50,
        -2.50,
        -0.50,
        0.50,
        2.50,
        -2.50,
        -0.50,
        0.50,
        2.50,
        3.30,
        1.10,
        1.10,
        3.30,
        3.30,
        1.10,
        1.10,
        3.30,
        3.30,
        1.10,
        1.10,
        3.30
    };
    ityp b_vec[N_MAX] =
    {
        3.30,
        1.10,
        1.10,
        3.30,
        3.30,
        1.10,
        1.10,
        3.30,
        3.30,
        1.10,
        1.10,
        3.30,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70
    };
    ityp c_vec[N_MAX] =
    {
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        6.70,
        -5.50,
        -0.50,
        0.50,
        4.50,
        -5.50,
        -0.50,
        0.50,
        4.50,
        -5.50,
        -0.50,
        0.50,
        4.50
    };
    ityp fx_vec[N_MAX] =
    {
        0.72356129348997784913,
        0.97911109345277961340,
        1.0216578140088564160,
        1.4051563200112126405,
        0.46961431639821611095,
        0.95296194977446325454,
        1.0512814213947987916,
        2.3999062904777858999,
        0.29106095928414718320,
        0.92536967910373175753,
        1.0865504094806997287,
        5.7381565526189046578,
        15090.669748704606754,
        -104.31170067364349677,
        21.175050707768812938,
        4.1946915819031922850,
        1.0170777974048815592E+10,
        -24708.635322489155868,
        1372.2304548384989560,
        58.092728706394652211,
        5.8682087615124176162E+18,
        -4.4635010147295996680E+08,
        5.3835057561295731310E+06,
        20396.913776019659426
    };
    ityp x_vec[N_MAX] =
    {
        0.25,
        0.25,
        0.25,
        0.25,
        0.55,
        0.55,
        0.55,
        0.55,
        0.85,
        0.85,
        0.85,
        0.85,
        0.25,
        0.25,
        0.25,
        0.25,
        0.55,
        0.55,
        0.55,
        0.55,
        0.85,
        0.85,
        0.85,
        0.85
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *a = *b = *c = *x = *fx = 0.00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *c = c_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _lerch_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LERCH_VALUES returns some values of the Lerch transcendent function.
  Discussion:
    The Lerch function is defined as
      Phi(z,s,a) = Sum ( 0 <= k < +oo ) z^k / ( a + k )^s
    omitting any terms with ( a + k ) = 0.
    In Mathematica, the function can be evaluated by:
      LerchPhi[z,s,a]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 September 2004
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *Z, the parameters of the function.
    Output, int *S, the parameters of the function.
    Output, double *A, the parameters of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pitpdtpit * const s_data = data; 
	
	dim_typ * n_data = s_data->a0;
	ityp * z = s_data->a1;
	ityp * a = s_data->a2;
	dim_typ * s = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 12

    ityp a_vec[N_MAX] =
    {
        0.0E+00,
        0.0E+00,
        0.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        2.0E+00,
        2.0E+00,
        2.0E+00,
        3.0E+00,
        3.0E+00,
        3.0E+00
    };

    ityp fx_vec[N_MAX] =
    {
        0.1644934066848226E+01,
        0.1202056903159594E+01,
        0.1000994575127818E+01,
        0.1164481052930025E+01,
        0.1074426387216080E+01,
        0.1000492641212014E+01,
        0.2959190697935714E+00,
        0.1394507503935608E+00,
        0.9823175058446061E-03,
        0.1177910993911311E+00,
        0.3868447922298962E-01,
        0.1703149614186634E-04
    };

    dim_typ s_vec[N_MAX] =
    {
        2, 3, 10,
        2, 3, 10,
        2, 3, 10,
        2, 3, 10
    };

    ityp z_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.5000000000000000E+00,
        0.5000000000000000E+00,
        0.5000000000000000E+00,
        0.3333333333333333E+00,
        0.3333333333333333E+00,
        0.3333333333333333E+00,
        0.1000000000000000E+00,
        0.1000000000000000E+00,
        0.1000000000000000E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *s = 0;
        *a = *fx = *z = 0.00;
    }
    else
    {
        *z = z_vec[*n_data-1];
        *s = s_vec[*n_data-1];
        *a = a_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _mertens_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    MERTENS_VALUES returns some values of the Mertens function.
  Discussion:
    The Mertens function M(N) is the sum from 1 to N of the Moebius
    function MU.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 October 2007
  Author:
    John Burkadt
  Reference:
    M Deleglise, J Rivat,
    Computing the Summation of the Moebius Function,
    Experimental Mathematics,
    Volume 5, 1996, pages 291-295.
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 2002,
    Second edition,
    ISBN: 1584883472,
    LC: QA5.W45.
  Parameters:
    Input/output, int *N_DATA.
    On input, if N_DATA is 0, the first test data is returned, and N_DATA
    is set to 1.  On each subsequent call, the input value of N_DATA is
    incremented and that test data item is returned, if available.  When
    there is no more test data, N_DATA is set to 0.
    Output, int *N, the argument of the Mertens function.
    Output, int *C, the value of the Mertens function.
*/
{
	const _2pdtps * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	short * c = s_data->a2;
	
    # define N_MAX 15

    short c_vec[N_MAX] =
    {
        1,   0,  -1,   -1,  -2,  -1,  -2,  -2,   -2,  -1,
        -2,  -2,   1,    2, -23
    };
    dim_typ n_vec[N_MAX] =
    {
        1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
        11,  12,  100, 1000, 10000
    };

    if ( N_MAX <= *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data];
        *c = c_vec[*n_data];
        ++ *n_data;
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _moebius_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOEBIUS_VALUES returns some values of the Moebius function.
  Discussion:
    MU(N) is defined as follows:
      MU(N) = 1 if N = 1;
              0 if N is divisible by the square of a prime;
        (-1)**K, if N is the product of K distinct primes.
    In Mathematica, the function can be evaluated by:
      MoebiusMu[n]
  First values:
     N  MU(N)
     1    1
     2   -1
     3   -1
     4    0
     5   -1
     6    1
     7   -1
     8    0
     9    0
    10    1
    11   -1
    12    0
    13   -1
    14    1
    15    1
    16    0
    17   -1
    18    0
    19   -1
    20    0
  Note:
    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
    if N is a square, cube, etc.
  Formula:
    The Moebius function is related to Euler's totient function:
      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 February 2003
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the Moebius function.
    Output, int *C, the value of the Moebius function.
*/
{
	const _2pdtps * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	short * c = s_data->a2;
	
    # define N_MAX 20

    short c_vec[N_MAX] =
    {
        1,  -1,  -1,   0,  -1,   1,  -1,   0,   0,   1,
        -1,   0,  -1,   1,   1,   0,  -1,   0,  -1,   0
    };

    dim_typ n_vec[N_MAX] =
    {
        1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
        11,  12,  13,  14,  15,  16,  17,  18,  19,  20
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _omega_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    OMEGA_VALUES returns some values of the OMEGA function.
  Discussion:
    In Mathematica, the function can be evaluated by
      Length [ FactorInteger [ n ] ]
  First values:
     N   OMEGA(N)
     1    1
     2    1
     3    1
     4    1
     5    1
     6    2
     7    1
     8    1
     9    1
    10    2
    11    1
    12    2
    13    1
    14    2
    15    2
    16    1
    17    1
    18    2
    19    1
    20    2
  Formula:
    If N = 1, then
      OMEGA(N) = 1
    else if the prime factorization of N is
      N = P1**E1 * P2**E2 * ... * PM**EM,
    then
      OMEGA(N) = M
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2003
  Author:

    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the OMEGA function.
    Output, int *C, the value of the OMEGA function.
*/
{
	const _2pdtpi * const s_data = data;
	
	dim_typ * n_data = s_data->a0;
	dim_typ * c = s_data->a1;
	int * n = s_data->a2;
	
    # define N_MAX 23

    dim_typ c_vec[N_MAX] =
    {
        1,   1,   1,   1,   1,
        2,   1,   1,   1,   2,
        3,   1,   4,   4,   3,
        1,   5,   2,   2,   1,
        6,   7,   8
    };

    int n_vec[N_MAX] =
    {
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        30,
        101,
        210,
        1320,
        1764,
        2003,
        2310,
        2827,
        8717,
        12553,
        30030,
        510510,
        9699690
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _partition_distinct_count_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    PARTITION_DISTINCT_COUNT_VALUES returns some values of Q(N).
  Discussion:
    A partition of an int *N is a representation of the integer
    as the sum of nonzero positive integers.  The order of the summands
    does not matter.  The number of partitions of N is symbolized
    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
    following partitions:
    5 = 5
      = 4 + 1
      = 3 + 2
      = 3 + 1 + 1
      = 2 + 2 + 1
      = 2 + 1 + 1 + 1
      = 1 + 1 + 1 + 1 + 1
    However, if we require that each member of the partition
    be distinct, so that no nonzero summand occurs more than once,
    we are computing something symbolized by Q(N).
    The number 5 has Q(N) = 3, because it has the following partitions
    into distinct parts:
    5 = 5
      = 4 + 1
      = 3 + 2
    In Mathematica, the function can be evaluated by
      PartitionsQ[n]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 February 2003
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the integer.
    Output, int *C, the number of partitions of the integer
    into distinct parts.
*/
{
	dim_typ ** const a_data = data;
	dim_typ * n_data = a_data[0];
	dim_typ * n = a_data[1];
	dim_typ * c = a_data[2];
	
    # define N_MAX 21

    dim_typ c_vec[N_MAX] =
    {
        1,
        1,   1,   2,   2,   3,   4,   5,   6,   8,  10,
        12,  15,  18,  22,  27,  32,  38,  46,  54,  64
    };

    dim_typ n_vec[N_MAX] =
    {
        0,
        1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _phi_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    PHI_VALUES returns some values of the PHI function.
  Discussion:
    PHI(N) is the number of integers between 1 and N which are
    relatively prime to N.  I and J are relatively prime if they
    have no common factors.  The function PHI(N) is known as
    "Euler's totient function".
    By convention, 1 and N are relatively prime.
    In Mathematica, the function can be evaluated by:
      EulerPhi[n]
  First values:
     N  PHI(N)
     1    1
     2    1
     3    2
     4    2
     5    4
     6    2
     7    6
     8    4
     9    6
    10    4
    11   10
    12    4
    13   12
    14    6
    15    8
    16    8
    17   16
    18    6
    19   18
    20    8
  Formula:
    PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
    PHI(P**K) = P**(K-1) * ( P - 1 ) if P is prime.
    PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
    N = Sum ( D divides N ) PHI(D).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 February 2003
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the PHI function.
    Output, int *C, the value of the PHI function.
*/
{
	dim_typ ** const a_data = data;
	dim_typ * n_data = a_data[0];
	dim_typ * n = a_data[1];
	dim_typ * c = a_data[2];
	
    # define N_MAX 20

    dim_typ c_vec[N_MAX] =
    {
        1,   1,   2,   2,   4,   2,   6,   4,   6,   4,
        8,   8,  16,  20,  16,  40, 148, 200, 200, 648
    };

    dim_typ n_vec[N_MAX] =
    {
        1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
        20,  30,  40,  50,  60, 100, 149, 500, 750, 999
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _r8_factorial_log_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_FACTORIAL_LOG_VALUES returns values of log(n!).
  Discussion:
    The function log(n!) can be written as
     log(n!) = sum ( 1 <= i <= n ) log ( i )
    In Mathematica, the function can be evaluated by:
      Log[n!]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2004
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the function.
    Output, double *FN, the value of the function.
*/
{
	const _2pdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * fn = s_data->a2;
	
    # define N_MAX 27

    ityp fn_vec[N_MAX] =
    {
        0.0000000000000000E+00,
        0.0000000000000000E+00,
        0.6931471805599453E+00,
        0.1791759469228055E+01,
        0.3178053830347946E+01,
        0.4787491742782046E+01,
        0.6579251212010101E+01,
        0.8525161361065414E+01,
        0.1060460290274525E+02,
        0.1280182748008147E+02,
        0.1510441257307552E+02,
        0.1750230784587389E+02,
        0.1998721449566189E+02,
        0.2255216385312342E+02,
        0.2519122118273868E+02,
        0.2789927138384089E+02,
        0.3067186010608067E+02,
        0.3350507345013689E+02,
        0.3639544520803305E+02,
        0.3933988418719949E+02,
        0.4233561646075349E+02,
        0.5800360522298052E+02,
        0.1484777669517730E+03,
        0.3637393755555635E+03,
        0.6050201058494237E+03,
        0.2611330458460156E+04,
        0.5912128178488163E+04
    };

    dim_typ n_vec[N_MAX] =
    {
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        25,
        50,
        100,
        150,
        500,
        1000
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *fn = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *fn = fn_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _sigma_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIGMA_VALUES returns some values of the Sigma function.
  Discussion:
    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
    In Mathematica, the function can be evaluated by:
      DivisorSigma[1,n]
  First values:
     N  SIGMA(N)
     1    1
     2    3
     3    4
     4    7
     5    6
     6   12
     7    8
     8   15
     9   13
    10   18
    11   12
    12   28
    13   14
    14   24
    15   24
    16   31
    17   18
    18   39
    19   20
    20   42
  Formula:
    SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
    SIGMA(P**K) = ( P**(K+1) - 1 ) / ( P - 1 ) if P is prime.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2003
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the Sigma function.
    Output, int *C, the value of the Sigma function.
*/
{
	dim_typ ** const a_data = data;
	dim_typ * n_data = a_data[0];
	dim_typ * n = a_data[1];
	dim_typ * c = a_data[2];
	
    # define N_MAX 20

    dim_typ c_vec[N_MAX] =
    {
        1,    3,    4,    7,    6,   12,    8,   15,   13,   18,
        72,  128,  255,  176,  576, 1170,  618,  984, 2232, 2340
    };

    dim_typ n_vec[N_MAX] =
    {
        1,   2,   3,   4,   5,   6,   7,   8,   9,   10,
        30, 127, 128, 129, 210, 360, 617, 815, 816, 1000
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _sin_power_int_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIN_POWER_INT_VALUES returns some values of the sine power integral.
  Discussion:
    The function has the form
      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin(T) )^N dt
    In Mathematica, the function can be evaluated by:
      Integrate [ ( Sin[x] )^n, { x, a, b } ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 September 2004
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *A, B, the limits of integration.
    Output, int *N, the power.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pitpdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	dim_typ * n = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 10

    ityp a_vec[N_MAX] =
    {
        0.10E+02,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.00E+00,
        0.10E+01,
        0.00E+00,
        0.00E+00
    };

    double b_vec[N_MAX] =
    {
        0.20E+02,
        0.10E+01,
        0.10E+01,
        0.10E+01,
        0.10E+01,
        0.10E+01,
        0.20E+01,
        0.20E+01,
        0.10E+01,
        0.10E+01
    };

    double fx_vec[N_MAX] =
    {
        0.10000000000000000000E+02,
        0.45969769413186028260E+00,
        0.27267564329357957615E+00,
        0.17894056254885809051E+00,
        0.12402556531520681830E+00,
        0.88974396451575946519E-01,
        0.90393123848149944133E+00,
        0.81495684202992349481E+00,
        0.21887522421729849008E-01,
        0.17023439374069324596E-01
    };

    dim_typ n_vec[N_MAX] =
    {
        0,
        1,
        2,
        3,
        4,
        5,
        5,
        5,
        10,
        11
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *a = *b = *fx = 0.00;
        *n = 0;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *  _spherical_harmonic_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERICAL_HARMONIC_VALUES returns values of spherical harmonic functions.
  Discussion:
    In Mathematica, the function can be evaluated by
      SphericalHarmonicY [ l, m, theta, phi ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 March 2005
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 1998.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *L, int *M, double THETA, PHI, the arguments
    of the function.
    Output, double *YR, *YI, the real and imaginary parts of
    the function.
*/
{
	const _2pdtps4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	dim_typ * l = s_data->a1;
	short * m = s_data->a2;
	ityp * theta = s_data->a3;
	ityp * phi = s_data->a4;
	ityp * yr = s_data->a5;
	ityp * yi = s_data->a6;
	
    # define N_MAX 20

    dim_typ l_vec[N_MAX] =
    {
        0,  1,  2,
        3,  4,  5,
        5,  5,  5,
        5,  4,  4,
        4,  4,  4,
        3,  3,  3,
        3,  3
    };
    short m_vec[N_MAX] =
    {
        0,  0,  1,
        2,  3,  5,
        4,  3,  2,
        1,  2,  2,
        2,  2,  2,
        -1, -1, -1,
        -1, -1
    };
    ityp phi_vec[N_MAX] =
    {
        0.1047197551196598E+01, 0.1047197551196598E+01, 0.1047197551196598E+01,
        0.1047197551196598E+01, 0.1047197551196598E+01, 0.6283185307179586E+00,
        0.6283185307179586E+00, 0.6283185307179586E+00, 0.6283185307179586E+00,
        0.6283185307179586E+00, 0.7853981633974483E+00, 0.7853981633974483E+00,
        0.7853981633974483E+00, 0.7853981633974483E+00, 0.7853981633974483E+00,
        0.4487989505128276E+00, 0.8975979010256552E+00, 0.1346396851538483E+01,
        0.1795195802051310E+01, 0.2243994752564138E+01
    };
    ityp theta_vec[N_MAX] =
    {
        0.5235987755982989E+00, 0.5235987755982989E+00, 0.5235987755982989E+00,
        0.5235987755982989E+00, 0.5235987755982989E+00, 0.2617993877991494E+00,
        0.2617993877991494E+00, 0.2617993877991494E+00, 0.2617993877991494E+00,
        0.2617993877991494E+00, 0.6283185307179586E+00, 0.1884955592153876E+01,
        0.3141592653589793E+01, 0.4398229715025711E+01, 0.5654866776461628E+01,
        0.3926990816987242E+00, 0.3926990816987242E+00, 0.3926990816987242E+00,
        0.3926990816987242E+00, 0.3926990816987242E+00
    };
    ityp yi_vec[N_MAX] =
    {
        0.0000000000000000E+00,  0.0000000000000000E+00, -0.2897056515173922E+00,
        0.1916222768312404E+00,  0.0000000000000000E+00,  0.0000000000000000E+00,
        0.3739289485283311E-02, -0.4219517552320796E-01,  0.1876264225575173E+00,
        -0.3029973424491321E+00,  0.4139385503112256E+00, -0.1003229830187463E+00,
        0.0000000000000000E+00, -0.1003229830187463E+00,  0.4139385503112256E+00,
        -0.1753512375142586E+00, -0.3159720118970196E+00, -0.3940106541811563E+00,
        -0.3940106541811563E+00, -0.3159720118970196E+00
    };
    ityp yr_vec[N_MAX] =
    {
        0.2820947917738781E+00,  0.4231421876608172E+00, -0.1672616358893223E+00,
        -0.1106331731112457E+00,  0.1354974113737760E+00,  0.5390423109043568E-03,
        -0.5146690442951909E-02,  0.1371004361349490E-01,  0.6096352022265540E-01,
        -0.4170400640977983E+00,  0.0000000000000000E+00,  0.0000000000000000E+00,
        0.0000000000000000E+00,  0.0000000000000000E+00,  0.0000000000000000E+00,
        0.3641205966137958E+00,  0.2519792711195075E+00,  0.8993036065704300E-01,
        -0.8993036065704300E-01, -0.2519792711195075E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *l = *m = 0;
        *theta = *phi = *yr = *yi = 0.00;
    }
    else
    {
        *l = l_vec[*n_data-1];
        *m = m_vec[*n_data-1];
        *theta = theta_vec[*n_data-1];
        *phi = phi_vec[*n_data-1];
        *yr = yr_vec[*n_data-1];
        *yi = yi_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _tau_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    TAU_VALUES returns some values of the Tau function.
  Discussion:
    TAU(N) is the number of divisors of N, including 1 and N.
    In Mathematica, the function can be evaluated by:
      DivisorSigma[1,n]
  First values:
     N   TAU(N)
     1    1
     2    2
     3    2
     4    3
     5    2
     6    4
     7    2
     8    4
     9    3
    10    4
    11    2
    12    6
    13    2
    14    4
    15    4
    16    5
    17    2
    18    6
    19    2
    20    6
  Formula:
    If the prime factorization of N is
      N = P1^E1 * P2^E2 * ... * PM^EM,
    then
      TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 February 2003
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the Tau function.
    Output, int *C, the value of the Tau function.
*/
{
	dim_typ ** const a_data = data;
	dim_typ * n_data = a_data[0];
	dim_typ * n = a_data[1];
	dim_typ * c = a_data[2];
	
    # define N_MAX 20

    dim_typ c_vec[N_MAX] =
    {
        1,  2,  2,  3,  2,  4,  2,  4,  3,  4,
        2, 12, 12,  4, 18, 24,  2,  8, 14, 28
    };

    dim_typ n_vec[N_MAX] =
    {
        1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
        23,  72, 126, 226, 300, 480, 521, 610, 832, 960
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *n = *c = 0;
    else
    {
        *n = n_vec[*n_data-1];
        *c = c_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _zeta_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    ZETA_VALUES returns some values of the Riemann Zeta function.
  Discussion:
    ZETA(N) = sum ( 1 <= I < +oo ) 1 / I**N
    In Mathematica, the function can be evaluated by:
      Zeta[n]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the Zeta function.
    Output, double *ZETA, the value of the Zeta function.
*/
{
	const _2pdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * zeta = s_data->a2;
	
    # define N_MAX 15

    dim_typ n_vec[N_MAX] =
    {
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        16,
        20,
        30,
        40
    };

    ityp zeta_vec[N_MAX] =
    {
        0.164493406684822643647E+01,
        0.120205690315959428540E+01,
        0.108232323371113819152E+01,
        0.103692775514336992633E+01,
        0.101734306198444913971E+01,
        0.100834927738192282684E+01,
        0.100407735619794433939E+01,
        0.100200839292608221442E+01,
        0.100099457512781808534E+01,
        0.100049418860411946456E+01,
        0.100024608655330804830E+01,
        0.100001528225940865187E+01,
        0.100000095396203387280E+01,
        0.100000000093132743242E+01,
        0.100000000000090949478E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *zeta = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *zeta = zeta_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _bessel_i0_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BESSEL_I0_VALUES returns some values of the I0 Bessel function.
  Discussion:
    The modified Bessel functions In(Z) and Kn(Z) are solutions of
    the differential equation
      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
    The modified Bessel function I0(Z) corresponds to N = 0.
    In Mathematica, the function can be evaluated by:
      BesselI[0,x]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 20

    ityp fx_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1010025027795146E+01,
        0.1040401782229341E+01,
        0.1092045364317340E+01,
        0.1166514922869803E+01,
        0.1266065877752008E+01,
        0.1393725584134064E+01,
        0.1553395099731217E+01,
        0.1749980639738909E+01,
        0.1989559356618051E+01,
        0.2279585302336067E+01,
        0.3289839144050123E+01,
        0.4880792585865024E+01,
        0.7378203432225480E+01,
        0.1130192195213633E+02,
        0.1748117185560928E+02,
        0.2723987182360445E+02,
        0.6723440697647798E+02,
        0.4275641157218048E+03,
        0.2815716628466254E+04
    };

    ityp x_vec[N_MAX] =
    {
        0.00,
        0.20,
        0.40,
        0.60,
        0.80,
        0.10E+01,
        0.12E+01,
        0.14E+01,
        0.16E+01,
        0.18E+01,
        0.20E+01,
        0.25E+01,
        0.30E+01,
        0.35E+01,
        0.40E+01,
        0.45E+01,
        0.50E+01,
        0.60E+01,
        0.80E+01,
        0.10E+02
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _bessel_i1_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BESSEL_I1_VALUES returns some values of the I1 Bessel function.
  Discussion:
    The modified Bessel functions In(Z) and Kn(Z) are solutions of
    the differential equation
      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
    In Mathematica, the function can be evaluated by:
      BesselI[1,x]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 20

    ityp fx_vec[N_MAX] =
    {
        0.0000000000000000,
        0.1005008340281251,
        0.2040267557335706,
        0.3137040256049221,
        0.4328648026206398,
        0.5651591039924850,
        0.7146779415526431,
        0.8860919814143274,
        0.1084810635129880E+01,
        0.1317167230391899E+01,
        0.1590636854637329E+01,
        0.2516716245288698E+01,
        0.3953370217402609E+01,
        0.6205834922258365E+01,
        0.9759465153704450E+01,
        0.1538922275373592E+02,
        0.2433564214245053E+02,
        0.6134193677764024E+02,
        0.3998731367825601E+03,
        0.2670988303701255E+04
    };

    ityp x_vec[N_MAX] =
    {
        0.00,
        0.20,
        0.40,
        0.60,
        0.80,
        0.10E+01,
        0.12E+01,
        0.14E+01,
        0.16E+01,
        0.18E+01,
        0.20E+01,
        0.25E+01,
        0.30E+01,
        0.35E+01,
        0.40E+01,
        0.45E+01,
        0.50E+01,
        0.60E+01,
        0.80E+01,
        0.10E+02
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _bessel_ix_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BESSEL_IX_VALUES returns some values of the Ix Bessel function.
  Discussion:
    This set of data considers the less common case in which the
    index of the Bessel function In is actually not an integer.
    We may suggest this case by occasionally replacing the symbol
    "In" by "Ix".
    The modified Bessel functions In(Z) and Kn(Z) are solutions of
    the differential equation
      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
    In Mathematica, the function can be evaluated by:
      BesselI[n,x]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2007
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *NU, the order of the function.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt3pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * nu = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 28

    ityp fx_vec[N_MAX] =
    {
        0.3592084175833614E+00,
        0.9376748882454876E+00,
        2.046236863089055E+00,
        3.053093538196718E+00,
        4.614822903407601E+00,
        26.47754749755907E+00,
        2778.784603874571E+00,
        4.327974627242893E+07,
        0.2935253263474798E+00,
        1.099473188633110E+00,
        21.18444226479414E+00,
        2500.906154942118E+00,
        2.866653715931464E+20,
        0.05709890920304825E+00,
        0.3970270801393905E+00,
        13.76688213868258E+00,
        2028.512757391936E+00,
        2.753157630035402E+20,
        0.4139416015642352E+00,
        1.340196758982897E+00,
        22.85715510364670E+00,
        2593.006763432002E+00,
        2.886630075077766E+20,
        0.03590910483251082E+00,
        0.2931108636266483E+00,
        11.99397010023068E+00,
        1894.575731562383E+00,
        2.716911375760483E+20
    };

    ityp nu_vec[N_MAX] =
    {
        0.50E+00,
        0.50E+00,
        0.50E+00,
        0.50E+00,
        0.50E+00,
        0.50E+00,
        0.50E+00,
        0.50E+00,
        1.50E+00,
        1.50E+00,
        1.50E+00,
        1.50E+00,
        1.50E+00,
        2.50E+00,
        2.50E+00,
        2.50E+00,
        2.50E+00,
        2.50E+00,
        1.25E+00,
        1.25E+00,
        1.25E+00,
        1.25E+00,
        1.25E+00,
        2.75E+00,
        2.75E+00,
        2.75E+00,
        2.75E+00,
        2.75E+00
    };

    ityp x_vec[N_MAX] =
    {
        0.2E+00,
        1.0E+00,
        2.0E+00,
        2.5E+00,
        3.0E+00,
        5.0E+00,
        10.0E+00,
        20.0E+00,
        1.0E+00,
        2.0E+00,
        5.0E+00,
        10.0E+00,
        50.0E+00,
        1.0E+00,
        2.0E+00,
        5.0E+00,
        10.0E+00,
        50.0E+00,
        1.0E+00,
        2.0E+00,
        5.0E+00,
        10.0E+00,
        50.0E+00,
        1.0E+00,
        2.0E+00,
        5.0E+00,
        10.0E+00,
        50.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *nu = *x = *fx = 0.00;
    }
    else
    {
        *nu = nu_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _cauchy_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    CAUCHY_CDF_VALUES returns some values of the Cauchy CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = CauchyDistribution [ mu, sigma ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *MU, the mean of the distribution.
    Output, double *SIGMA, the variance of the distribution.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * sigma = s_data->a2;
	ityp * x = s_data->a3; 
	ityp * fx = s_data->a4;
	
    # define N_MAX 12

    ityp fx_vec[N_MAX] =
        {
        0.5000000000000000,
        0.8524163823495667,
        0.9220208696226307,
        0.9474315432887466,
        0.6475836176504333,
        0.6024163823495667,
        0.5779791303773693,
        0.5628329581890012,
        0.6475836176504333,
        0.5000000000000000,
        0.3524163823495667,
        0.2500000000000000
    };

    ityp mu_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01
    };

    ityp sigma_vec[N_MAX] =
    {
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01
    };

    ityp x_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01 };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *mu = *sigma = *x = *fx = 0.00;
    }
    else
    {
        *mu = mu_vec[*n_data-1];
        *sigma = sigma_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
# undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _exponential_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    EXPONENTIAL_CDF_VALUES returns some values of the Exponential CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = ExponentialDistribution [ lambda ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *LAMBDA, the parameter of the distribution.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt3pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * lambda = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 9

    ityp fx_vec[N_MAX] =
    {
        0.3934693402873666,
        0.6321205588285577,
        0.7768698398515702,
        0.8646647167633873,
        0.8646647167633873,
        0.9816843611112658,
        0.9975212478233336,
        0.9996645373720975,
        0.9999546000702375
    };

    ityp lambda_vec[N_MAX] =
    {
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01
    };

    ityp x_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *lambda = *x = *fx = 0.00;
    }
    else
    {
        *lambda = lambda_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _extreme_values_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    EXTREME_VALUES_CDF_VALUES returns some values of the Extreme Values CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = ExtremeValuesDistribution [ alpha, beta ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *ALPHA, the first parameter of the distribution.
    Output, double *BETA, the second parameter of the distribution.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * beta = s_data->a2;
	ityp * x = s_data->a3; 
	ityp * fx = s_data->a4;
	
    # define N_MAX 12

    ityp alpha_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01
    };

    ityp beta_vec[N_MAX] =
    {
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01
    };

    ityp fx_vec[N_MAX] =
    {
        0.3678794411714423,
        0.8734230184931166,
        0.9818510730616665,
        0.9975243173927525,
        0.5452392118926051,
        0.4884435800065159,
        0.4589560693076638,
        0.4409910259429826,
        0.5452392118926051,
        0.3678794411714423,
        0.1922956455479649,
        0.6598803584531254E-01
    };

    ityp x_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *alpha = *beta = *x = *fx = 0.00;
    }
    else
    {
        *alpha = alpha_vec[*n_data-1];
        *beta = beta_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *    _geometric_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEOMETRIC_CDF_VALUES returns values of the geometric CDF.
  Discussion:
    The geometric or Pascal probability density function gives the
    probability that the first success will happen on the X-th Bernoulli
    trial, given that the probability of a success on a single trial is P.
    The value of CDF ( X, P ) is the probability that the first success
    will happen on or before the X-th trial.
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`DiscreteDistributions`]
      dist = GeometricDistribution [ p ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 August 2004
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
    Daniel Zwillinger and Stephen Kokoska,
    CRC Standard Probability and Statistics Tables and Formulae,
    Chapman and Hall / CRC Press, 2000.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *X, the number of trials.
    Output, double *P, the probability of success
    on one trial.
    Output, double *CDF, the cumulative density function.
*/
{
	const _2pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * x = s_data->a1;
	ityp * p = s_data->a2;
	ityp * cdf = s_data->a3;
	
    # define N_MAX 14

    ityp cdf_vec[N_MAX] =
    {
        0.1900000000000000,
        0.2710000000000000,
        0.3439000000000000,
        0.6861894039100000,
        0.3600000000000000,
        0.4880000000000000,
        0.5904000000000000,
        0.9141006540800000,
        0.7599000000000000,
        0.8704000000000000,
        0.9375000000000000,
        0.9843750000000000,
        0.9995117187500000,
        0.9999000000000000
    };

    ityp p_vec[N_MAX] =
    {
        0.10,
        0.10,
        0.10,
        0.10,
        0.20,
        0.20,
        0.20,
        0.20,
        0.30,
        0.40,
        0.50,
        0.50,
        0.50,
        0.90
    };

    dim_typ x_vec[N_MAX] =
    {
        1,  2,  3, 10, 1,
        2,  3, 10,  3, 3,
        3,  5, 10,  3
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *x = 0;
        *p = *cdf = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *p = p_vec[*n_data-1];
        *cdf = cdf_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _laplace_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAPLACE_CDF_VALUES returns some values of the Laplace CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = LaplaceDistribution [ mu, beta ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *MU, the mean of the distribution.
    Output, double *BETA, the shape parameter.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * beta = s_data->a2;
	ityp * x = s_data->a3; 
	ityp * fx = s_data->a4;
	
    # define N_MAX 12

    ityp beta_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01
    };

    ityp fx_vec[N_MAX] =
    {
        0.5000000000000000,
        0.8160602794142788,
        0.9323323583816937,
        0.9751064658160680,
        0.6967346701436833,
        0.6417343447131054,
        0.6105996084642976,
        0.5906346234610091,
        0.5000000000000000,
        0.3032653298563167,
        0.1839397205857212,
        0.1115650800742149
    };

    ityp mu_vec[N_MAX] =
    {
        0.0000000000000000E+01,
        0.0000000000000000E+01,
        0.0000000000000000E+01,
        0.0000000000000000E+01,
        0.0000000000000000E+01,
        0.0000000000000000E+01,
        0.0000000000000000E+01,
        0.0000000000000000E+01,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01
    };

    ityp x_vec[N_MAX] =
    {
        0.0000000000000000E+01,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *mu = *beta = *x = *fx = 0.00;
    }
    else
    {
        *mu = mu_vec[*n_data-1];
        *beta = beta_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _log_normal_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOG_NORMAL_CDF_VALUES returns some values of the Log Normal CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = LogNormalDistribution [ mu, sigma ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *MU, the mean of the distribution.
    Output, double *SIGMA, the shape parameter of the distribution.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * sigma = s_data->a2;
	ityp * x = s_data->a3; 
	ityp * fx = s_data->a4;
	
    # define N_MAX 12

    ityp fx_vec[N_MAX] =
    {
        0.2275013194817921E-01,
        0.2697049307349095,
        0.5781741008028732,
        0.7801170895122241,
        0.4390310097476894,
        0.4592655190218048,
        0.4694258497695908,
        0.4755320473858733,
        0.3261051056816658,
        0.1708799040927608,
        0.7343256357952060E-01,
        0.2554673736161761E-01
    };

    ityp mu_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01
    };

    ityp sigma_vec[N_MAX] =
    {
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01
    };

    ityp x_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *mu = *sigma = *x = *fx = 0.00;
    }
    else
    {
        *mu = mu_vec[*n_data-1];
        *sigma = sigma_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _log_series_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOG_SERIES_CDF_VALUES returns some values of the log series CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`DiscreteDistributions`]
      dist = LogSeriesDistribution [ t ]
      CDF [ dist, n ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 August 2004
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *T, the parameter of the function.
    Output, int *N, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data;
	
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * t = s_data->a2;
	ityp * fx =  s_data->a3;
	
    # define N_MAX 29

    ityp fx_vec[N_MAX] =
    {
        0.9491221581029903,
        0.9433541128559735,
        0.9361094611773272,
        0.9267370278044118,
        0.9141358246245129,
        0.8962840235449100,
        0.8690148741955517,
        0.8221011541254772,
        0.7213475204444817,
        0.6068261510845583,
        0.5410106403333613,
        0.4970679476476894,
        0.4650921887927060,
        0.4404842934597863,
        0.4207860535926143,
        0.4045507673897055,
        0.3908650337129266,
        0.2149757685421097,
        0.0000000000000000,
        0.2149757685421097,
        0.3213887739704539,
        0.3916213575531612,
        0.4437690508633213,
        0.4850700239649681,
        0.5191433267738267,
        0.5480569580144867,
        0.5731033910767085,
        0.5951442521714636,
        0.6147826594068904
    };

    dim_typ n_vec[N_MAX] =
    {
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 0, 1,
        2, 3, 4, 5, 6,
        7, 8, 9, 10
    };

    ityp t_vec[N_MAX] =
    {
        0.1000000000000000,
        0.1111111111111111,
        0.1250000000000000,
        0.1428571428571429,
        0.1666666666666667,
        0.2000000000000000,
        0.2500000000000000,
        0.3333333333333333,
        0.5000000000000000,
        0.6666666666666667,
        0.7500000000000000,
        0.8000000000000000,
        0.8333333333333333,
        0.8571485714857149,
        0.8750000000000000,
        0.8888888888888889,
        0.9000000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000,
        0.9900000000000000
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *t = *fx = 0.00;
    }
    else
    {
        *t = t_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _logistic_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOGISTIC_CDF_VALUES returns some values of the Logistic CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = LogisticDistribution [ mu, beta ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *MU, the mean of the distribution.
    Output, double *BETA, the shape parameter of the distribution.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * beta = s_data->a2;
	ityp * x = s_data->a3; 
	ityp * fx = s_data->a4;
	
    # define N_MAX 12

    ityp beta_vec[N_MAX] =
    {
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01
    };

    ityp fx_vec[N_MAX] =
    {
        0.5000000000000000,
        0.8807970779778824,
        0.9820137900379084,
        0.9975273768433652,
        0.6224593312018546,
        0.5825702064623147,
        0.5621765008857981,
        0.5498339973124779,
        0.6224593312018546,
        0.5000000000000000,
        0.3775406687981454,
        0.2689414213699951
    };

    ityp mu_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01
    };

    ityp x_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *mu = *beta = *x = *fx = 0.00;
    }
    else
    {
        *mu = mu_vec[*n_data-1];
        *beta = beta_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _rayleigh_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    RAYLEIGH_CDF_VALUES returns some values of the Rayleigh CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = RayleighDistribution [ sigma ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 September 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *SIGMA, the shape parameter of the distribution.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt3pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * sigma = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3; 
	
    # define N_MAX 9

    ityp fx_vec[N_MAX] =
    {
        0.8646647167633873,
        0.9996645373720975,
        0.9999999847700203,
        0.999999999999987,
        0.8646647167633873,
        0.3934693402873666,
        0.1992625970831920,
        0.1175030974154046,
        0.7688365361336422E-01
    };

    ityp sigma_vec[N_MAX] =
    {
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01
    };

    ityp x_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *sigma = *x = *fx = 0.00;
    }
    else
    {
        *sigma = sigma_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _von_mises_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    VON_MISES_CDF_VALUES returns some values of the von Mises CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2004
  Author:
    John Burkardt
  Reference:
    Kanti Mardia, Peter Jupp,
    Directional Statistics,
    Wiley, 2000, QA276.M335
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *A, B, the parameters of the function.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x = s_data->a3; 
	ityp * fx = s_data->a4;
	
    # define N_MAX 23

    ityp a_vec[N_MAX] =
    {
        0.00,
        0.00,
        0.00,
        0.00,
        0.00,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        -0.2E+01,
        -0.1E+01,
        0.0E+01,
        0.1E+01,
        0.2E+01,
        0.3E+01,
        0.00,
        0.00,
        0.00,
        0.00,
        0.00,
        0.00
    };

    ityp b_vec[N_MAX] =
    {
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.1E+01,
        0.2E+01,
        0.2E+01,
        0.2E+01,
        0.2E+01,
        0.2E+01,
        0.2E+01,
        0.3E+01,
        0.3E+01,
        0.3E+01,
        0.3E+01,
        0.3E+01,
        0.3E+01,
        0.00,
        0.1E+01,
        0.2E+01,
        0.3E+01,
        0.4E+01,
        0.5E+01
    };

    ityp fx_vec[N_MAX] =
    {
        0.2535089956281180E-01,
        0.1097539041177346,
        0.5000000000000000,
        0.8043381312498558,
        0.9417460124555197,
        0.5000000000000000,
        0.6018204118446155,
        0.6959356933122230,
        0.7765935901304593,
        0.8410725934916615,
        0.8895777369550366,
        0.9960322705517925,
        0.9404336090170247,
        0.5000000000000000,
        0.5956639098297530E-01,
        0.3967729448207649E-02,
        0.2321953958111930E-03,
        0.6250000000000000,
        0.7438406999109122,
        0.8369224904294019,
        0.8941711407897124,
        0.9291058600568743,
        0.9514289900655436
    };

    ityp x_vec[N_MAX] =
    {
        -0.2617993977991494E+01,
        -0.1570796326794897E+01,
        0.0000000000000000,
        0.1047197551196598E+01,
        0.2094395102393195E+01,
        0.1000000000000000E+01,
        0.1200000000000000E+01,
        0.1400000000000000E+01,
        0.1600000000000000E+01,
        0.1800000000000000E+01,
        0.2000000000000000E+01,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        0.0000000000000000,
        0.7853981633974483,
        0.7853981633974483,
        0.7853981633974483,
        0.7853981633974483,
        0.7853981633974483,
        0.7853981633974483
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *a = *b = *x = *fx = 0.00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _weibull_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEIBULL_CDF_VALUES returns some values of the Weibull CDF.
  Discussion:
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`ContinuousDistributions`"]
      dist = WeibullDistribution [ alpha, beta ]
      CDF [ dist, x ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 August 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *ALPHA, the first parameter of the distribution.
    Output, double *BETA, the second parameter of the distribution.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * alpha = s_data->a1;
	ityp * beta = s_data->a2;
	ityp * x = s_data->a3; 
	ityp * fx = s_data->a4;
	
    # define N_MAX 12

    ityp alpha_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01
    };

    ityp beta_vec[N_MAX] =
    {
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.5000000000000000,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.5000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01
    };

    ityp fx_vec[N_MAX] =
    {
        0.8646647167633873,
        0.9816843611112658,
        0.9975212478233336,
        0.9996645373720975,
        0.6321205588285577,
        0.4865828809674080,
        0.3934693402873666,
        0.3296799539643607,
        0.8946007754381357,
        0.9657818816883340,
        0.9936702845725143,
        0.9994964109502630
    };

    ityp x_vec[N_MAX] =
    {
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.4000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01,
        0.3000000000000000E+01
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *alpha = *beta = *x = *fx = 0.00;
    }
    else
    {
        *alpha = alpha_vec[*n_data-1];
        *beta = beta_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _frobenius_number_order2_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    FROBENIUS_NUMBER_ORDER2_VALUES returns values of the order 2 Frobenius number.
  Discussion:
    The Frobenius number of order N is the solution of the Frobenius
    coin sum problem for N coin denominations.
    The Frobenius coin sum problem assumes the existence of
    N coin denominations, and asks for the largest value that cannot
    be formed by any combination of coins of these denominations.
    The coin denominations are assumed to be distinct positive integers.
    For general N, this problem is fairly difficult to handle.
    For N = 2, it is known that:
    * if C1 and C2 are not relatively prime, then
      there are infinitely large values that cannot be formed.
    * otherwise, the largest value that cannot be formed is
      C1 * C2 - C1 - C2, and that exactly half the values between
      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
    As a simple example, if C1 = 2 and C2 = 7, then the largest
    unrepresentable value is 5, and there are (5+1)/2 = 3
    unrepresentable values, namely 1, 3, and 5.
    For a general N, and a set of coin denominations C1, C2, ..., CN,
    the Frobenius number F(N, C(1:N) ) is defined as the largest value
    B for which the equation
      C1*X1 + C2*X2 + ... + CN*XN = B
    has no nonnegative integer solution X(1:N).
    In Mathematica, the Frobenius number can be determined by
      FrobeniusNumber[ {C1,...,CN} ]
  Modified:
    28 November 2007
  Author:
    John Burkardt
  Reference:
    James Sylvester,
    Question 7382,
    Mathematical Questions with their Solutions,
    Educational Times,
    Volume 41, page 21, 1884.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *C1, *C2, the parameters of the function.
    Output, int *F, the value of the function.
*/
{
	dim_typ ** const a_data = data;
	dim_typ * n_data = a_data[0];
	dim_typ * c1 = a_data[1];
	dim_typ * c2 = a_data[2];
	dim_typ * f = a_data[3];
	
    # define N_MAX 6

    dim_typ c1_vec[N_MAX] =
    {
        2,
        3,
        4,
        5,
        12,
        99
    };
    dim_typ c2_vec[N_MAX] =
    {
        5,
        17,
        19,
        13,
        11,
        100
    };
    dim_typ f_vec[N_MAX] =
    {
        3,
        31,
        53,
        47,
        109,
        9701
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *c1 = *c2 = *f = 0;
    else
    {
        *c1 = c1_vec[*n_data-1];
        *c2 = c2_vec[*n_data-1];
        *f = f_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _r8_fall_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_FALL_VALUES returns some values of the falling factorial function.
  Discussion:
    In Mathematica, the function can be evaluated by:
      FactorialPower[X,N]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 December 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, int *N, the arguments of the function.
    Output, double *F, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data;
	
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * f =  s_data->a3;
	
    # define N_MAX 15

    const ityp f_vec[N_MAX] =
    {
        120.0000000000000,
        163.1601562500000,
        216.5625000000000,
        281.6601562500000,
        360.0000000000000,
        1.000000000000000,
        7.500000000000000,
        48.75000000000000,
        268.1250000000000,
        1206.562500000000,
        4222.968750000000,
        10557.42187500000,
        15836.13281250000,
        7918.066406250000,
        -3959.03320312500
    };

    const dim_typ n_vec[N_MAX] =
    {
        4,
        4,
        4,
        4,
        4,
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9
    };

    const ityp x_vec[N_MAX] =
    {
        5.00,
        5.25,
        5.50,
        5.75,
        6.00,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50
    };


    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *f = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *f = f_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _r8_rise_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_RISE_VALUES returns some values of the rising factorial function.
  Discussion:
    Pochhammer(X,Y) = Gamma(X+Y) / Gamma(X)
    For integer arguments, Pochhammer(M,N) = ( M + N - 1 )! / ( N - 1 )!
    In Mathematica, the function can be evaluated by:
      Pochammer[X,N]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 December 2014
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, int *N, the arguments of the function.
    Output, double *F, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data;
	
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * f =  s_data->a3;
	
    # define N_MAX 15

    const ityp f_vec[N_MAX] =
    {
        1680.000000000000,
        1962.597656250000,
        2279.062500000000,
        2631.972656250000,
        3024.000000000000,
        1.000000000000000,
        7.500000000000000,
        63.75000000000000,
        605.6250000000000,
        6359.062500000000,
        73129.21875000000,
        914115.2343750000,
        1.234055566406250E+07,
        1.789380571289063E+08,
        2.773539885498047E+09
    };

    const dim_typ n_vec[N_MAX] =
    {
        4,
        4,
        4,
        4,
        4,
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9
    };

    const ityp x_vec[N_MAX] =
    {
        5.00,
        5.25,
        5.50,
        5.75,
        6.00,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50,
        7.50
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *x = *f = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *f = f_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _lambert_w_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAMBERT_W_VALUES returns some values of the Lambert W function.
  Discussion:
    The function W(X) is defined implicitly by:
      W(X) * e^W(X) = X
    The function is also known as the "Omega" function.
    In Mathematica, the function can be evaluated by:
      W = ProductLog [ X ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 February 2005
  Author:
    John Burkardt
  Reference:
    Brian Hayes,
    "Why W?",
    The American Scientist,
    Volume 93, March-April 2005, pages 104-108.
    Eric Weisstein,
    "Lambert's W-Function",
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 1998.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 22

    static ityp fx_vec[N_MAX] =
    {
        0.0000000000000000E+00,
        0.3517337112491958E+00,
        0.5671432904097839E+00,
        0.7258613577662263E+00,
        0.8526055020137255E+00,
        0.9585863567287029E+00,
        0.1000000000000000E+01,
        0.1049908894964040E+01,
        0.1130289326974136E+01,
        0.1202167873197043E+01,
        0.1267237814307435E+01,
        0.1326724665242200E+01,
        0.1381545379445041E+01,
        0.1432404775898300E+01,
        0.1479856830173851E+01,
        0.1524345204984144E+01,
        0.1566230953782388E+01,
        0.1605811996320178E+01,
        0.1745528002740699E+01,
        0.3385630140290050E+01,
        0.5249602852401596E+01,
        0.1138335808614005E+02
    };
    static ityp x_vec[N_MAX] =
    {
        0.0000000000000000E+00,
        0.5000000000000000E+00,
        0.1000000000000000E+01,
        0.1500000000000000E+01,
        0.2000000000000000E+01,
        0.2500000000000000E+01,
        0.2718281828459045E+01,
        0.3000000000000000E+01,
        0.3500000000000000E+01,
        0.4000000000000000E+01,
        0.4500000000000000E+01,
        0.5000000000000000E+01,
        0.5500000000000000E+01,
        0.6000000000000000E+01,
        0.6500000000000000E+01,
        0.7000000000000000E+01,
        0.7500000000000000E+01,
        0.8000000000000000E+01,
        0.1000000000000000E+02,
        0.1000000000000000E+03,
        0.1000000000000000E+04,
        0.1000000000000000E+07
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x  = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _r8_factorial2_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_FACTORIAL2_VALUES returns values of the double factorial function.
  Formula:
    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 ) (N even)
                    = Product ( N * (N-2) * (N-4) * ... * 1 ) (N odd)
    In Mathematica, the function can be evaluated by:
      n!!
  Example:
     N    N!!
     0     1
     1     1
     2     2
     3     3
     4     8
     5    15
     6    48
     7   105
     8   384
     9   945
    10  3840
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 February 2015
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
    Daniel Zwillinger,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996, page 16.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *N, the argument of the function.
    Output, double *F, the value of the function.
*/
{
	const _2pdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * n = s_data->a1;
	ityp * f = s_data->a2;
	
    # define N_MAX 16

    static ityp f_vec[N_MAX] =
    {
        1.00,
        1.00,
        2.00,
        3.00,
        8.00,
        15.00,
        48.00,
        105.00,
        384.00,
        945.00,
        3840.00,
        10395.00,
        46080.00,
        135135.00,
        645120.00,
        2027025.00
    };

    static dim_typ n_vec[N_MAX] =
    {
        0,
        1,  2,  3,  4,  5,
        6,  7,  8,  9, 10,
        11, 12, 13, 14, 15
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *n = 0;
        *f = 0.00;
    }
    else
    {
        *n = n_vec[*n_data-1];
        *f = f_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE  void *   _bivariate_normal_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    BIVARIATE_NORMAL_CDF_VALUES returns some values of the bivariate normal CDF.
  Discussion:
    FXY is the probability that two variables A and B, which are
    related by a bivariate normal distribution with correlation R,
    respectively satisfy A <= X and B <= Y.
    Mathematica can evaluate the bivariate normal CDF via the commands:
      <<MultivariateStatistics`
      cdf = CDF[MultinormalDistribution[{0,0}{{1,r},{r,1}}],{x,y}]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 November 2010
  Author:
    John Burkardt
  Reference:
    National Bureau of Standards,
    Tables of the Bivariate Normal Distribution and Related Functions,
    NBS, Applied Mathematics Series, Number 50, 1959.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *X, *Y, the parameters of the function.
    Output, double *R, the correlation value.
    Output, double *FXY, the value of the function.
*/
{
	const pdt4pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * r = s_data->a3; 
	ityp * fxy = s_data->a4;
	
    #  define N_MAX 41

    static ityp fxy_vec[N_MAX] =
    {
        0.02260327218569867E+00,
        0.1548729518584100E+00,
        0.4687428083352184E+00,
        0.7452035868929476E+00,
        0.8318608306874188E+00,
        0.8410314261134202E+00,
        0.1377019384919464E+00,
        0.1621749501739030E+00,
        0.1827411243233119E+00,
        0.2010067421506235E+00,
        0.2177751155265290E+00,
        0.2335088436446962E+00,
        0.2485057781834286E+00,
        0.2629747825154868E+00,
        0.2770729823404738E+00,
        0.2909261168683812E+00,
        0.3046406378726738E+00,
        0.3183113449213638E+00,
        0.3320262544108028E+00,
        0.3458686754647614E+00,
        0.3599150462310668E+00,
        0.3742210899871168E+00,
        0.3887706405282320E+00,
        0.4032765198361344E+00,
        0.4162100291953678E+00,
        0.6508271498838664E+00,
        0.8318608306874188E+00,
        0.0000000000000000,
        0.1666666666539970,
        0.2500000000000000,
        0.3333333333328906,
        0.5000000000000000,
        0.7452035868929476,
        0.1548729518584100,
        0.1548729518584100,
        0.06251409470431653,
        0.7452035868929476,
        0.1548729518584100,
        0.1548729518584100,
        0.06251409470431653,
        0.6337020457912916
    };
    static ityp r_vec[N_MAX] =
    {
        0.500,  0.500,  0.500,  0.500,  0.500,
        0.500, -0.900, -0.800, -0.700, -0.600,
        -0.500, -0.400, -0.300, -0.200, -0.100,
        0.000,  0.100,  0.200,  0.300,  0.400,
        0.500,  0.600,  0.700,  0.800,  0.900,
        0.673,  0.500, -1.000, -0.500,  0.000,
        0.500,  1.000,  0.500,  0.500,  0.500,
        0.500,  0.500,  0.500,  0.500,  0.500,
        0.500
    };
    static ityp x_vec[N_MAX] =
    {
        -2.00, -1.00,  0.00,  1.00,  2.00,
        3.00, -0.20, -0.20, -0.20, -0.20,
        -0.20, -0.20, -0.20, -0.20, -0.20,
        -0.20, -0.20, -0.20, -0.20, -0.20,
        -0.20, -0.20, -0.20, -0.20, -0.20,
        1.00,  2.00,  0.00,  0.00,  0.00,
        0.00,  0.00,  1.00,  1.00, -1.00,
        -1.00,  1.00,  1.00, -1.00, -1.00,
        0.7071067811865475
    };
    static ityp y_vec[N_MAX] =
    {
        1.00,  1.00,  1.00,  1.00,  1.00,
        1.00,  0.50,  0.50,  0.50,  0.50,
        0.50,  0.50,  0.50,  0.50,  0.50,
        0.50,  0.50,  0.50,  0.50,  0.50,
        0.50,  0.50,  0.50,  0.50,  0.50,
        0.50,  1.00,  0.00,  0.00,  0.00,
        0.00,  0.00,  1.00, -1.00,  1.00,
        -1.00,  1.00, -1.00,  1.00, -1.00,
        0.7071067811865475
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *r = *x = *y = *fxy = 0.00;
    }
    else
    {
        *r = r_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *y = y_vec[*n_data-1];
        *fxy = fxy_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
