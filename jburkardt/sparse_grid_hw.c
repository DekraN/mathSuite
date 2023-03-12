#ifndef __DISABLEDEEP_SPARSEGRIDHW

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8cvv_offset ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8CVV_OFFSET determines the row offsets of an R8CVV.
  Discussion:
    An R8CVV is a "vector of vectors" of R8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in the array.
    Input, int NR(M), the row sizes.
    Output, int R8CVV_OFFSET[M+1], the row offsets.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	int * nr = s_data->a1;
	
	int *roff = ( int * ) malloc ( ( m + 1 ) * sizeof ( int ) );
	roff[0] = 0;
	for (dim_typ i = 0; i < m; ++i )
		roff[i+1] = roff[i] + nr[i];
	return roff;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8cvv_rget_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8CVV_RGET_NEW gets row I from an R8CVV.
  Discussion:
    An R8CVV is a "vector of vectors" of R8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2012
  Author:
    John Burkardt
  Parameters:
    Input, int MN, the size of the cell array.
    Input, double A[MN], the cell array.
    Input, int M, the number of rows in the array.
    Input, int ROFF[M+1], the row offsets.
    Input, int I, the row.
    0 <= I < M.
    Output, double R8CVV_RGET_NEW[NR[I]], the value of A(I,*).
*/
{
	const _2dtpiipit * const s_data = data;
	
	const register dim_typ mn = s_data->a0;
	const register dim_typ m = s_data->a1;
	int * roff = s_data->a2;
	int i = s_data->a3;
	ityp * a = s_data->a4;
	
	ityp *ai;
	int k1;
	int k2;
	int nv;
	
	k1 = roff[i];
	k2 = roff[i+1];
	nv = k2 - k1;
	ai = ( ityp * ) malloc ( nv * sizeof ( ityp ) );
	for (dim_typ j = 0; j < nv; ++j )
		ai[j] = a[k1+j];
	
	return ai;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8cvv_rset ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8CVV_RSET sets row I from an R8CVV.
  Discussion:
    An R8CVV is a "vector of vectors" of R8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2012
  Author:
    John Burkardt
  Parameters:
    Input, int MN, the size of the cell array.
    Input/output, double A[MN], the cell array.
    Input, int M, the number of rows in the array.
    Input, int ROFF[M+1], the row offsets.
    Input, int I, the row.
    0 <= I < M.
    Input, double AI[NR[I]], the new value of A(I,*).
*/
{
	const dtpitdtpiipit * const s_data = data;
	const register dim_typ mn = s_data->a0;
	ityp * a = s_data->a1;
	const register dim_typ m = s_data->a2;
	int * roff = s_data->a3;
	int i = s_data->a4;
	ityp * ai = s_data->a5;
	
	int k1;
	int k2;
	int nv;
	
	k1 = roff[i];
	k2 = roff[i+1];
	nv = k2 - k1;
	for (dim_typ j = 0; j < nv; ++j )
		a[k1+j] = ai[j];
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cce_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    CCE_ORDER: order of a Clenshaw-Curtis Exponential rule from the level.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 February 2014
  Author:
    John Burkardt.
  Parameters:
    Input, int L, the level of the rule.
    1 <= L.
    Output, int CCE_ORDER, the order of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ l = *(dim_typ *) data;	
	
	result = l<1?USHRT_MAX : 1 + (l!=1)*powi ( 2, l - 1 ) + 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ccl_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    CCL_ORDER computes the order of a CCL rule from the level.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 February 2014
  Author:
    John Burkardt.
  Parameters:
    Input, int L, the level of the rule.
    1 <= L.
    Output, int CCL_ORDER, the order of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ l = *(dim_typ *) data;	
	
	result = l<1?USHRT_MAX : (l<<1)-1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ccs_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    CCS_ORDER: order of a "slow growth" Clenshaw Curtis quadrature rule.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    The rule is defined on [0,1].
    The integral to approximate:
      Integral ( 0 <= X <= 1 ) F(X) dX
    The quadrature rule:
      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
    The input value L requests a rule of precision at least 2*L-1.
    In order to preserve nestedness, this function returns the order
    of a rule which is the smallest value of the form 1+2^E which
    is greater than or equal to 2*L-1.
     L  2*L-1   N
    --  -----  --
     1      1   1
     2      3   3
     3      5   5
     4      7   9
     5      9   9
     6     11  17
     7     13  17
     8     15  17
     9     17  17
    10     19  33
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt
  Parameters:
    Input, int L, the level of the rule.
    1 <= L.
    Output, int CCS_ORDER, the appropriate order.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ l = *(dim_typ *) data;	
	
    int n;

    if ( l < 1 )
    {
    	result = USHRT_MAX;
        return &result;
    }
    /*
    Find the order N that satisfies the precision requirement.
    */
    if ( l == 1 )
        n = 1;
    else
    {
        n = 3;
        while ( n < (l<<1) - 1 )
            n = (n<<1) - 1;
    }

	result = n;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fn_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    FN_INTEGRAL is the integral of the Hermite test function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int D, the spatial dimension.
    Output, double FN_INTEGRAL, the integral value.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ d = *(dim_typ *) data;	
	
	result = ( ityp ) ( i4_factorial2 ( 5 ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _fn_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    FN_VALUE is a Hermite test function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 May 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int D, the spatial dimension.
    Input, int N, the number of points.
    Input, double X[D*N], the points.
    Output, double FN_VALUE[N], the function values.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ d = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	
    ityp *fx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i)
        fx[i] = pow ( x[0+i*d], 6 );
    return fx;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fu_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    FU_INTEGRAL is the integral of the test function for the [0,1]^D interval.
  Discussion:
    The same function, integrated over [-1,+1]^D, has an integral
    that is 2^D times larger.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int D, the spatial dimension.
    Output, double FU_INTEGRAL, the integral value.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ d = *(dim_typ *) data;	
	
	result = pow ( 0.50 * erf ( 0.50 / sqrt ( 2.00 ) ), d );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _fu_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    FU_VALUE is a sample function for the [0,1]^D interval.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int D, the spatial dimension.
    Input, int N, the number of points.
    Input, double X[D*N], the points.
    Output, double FU_VALUE[N], the function values.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ d = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	
    dim_typ i, j;
    ityp *fx = ( double * ) malloc ( n * sizeof ( double ) );

    for ( j = 0; j < n; ++j )
    {
        fx[j] = 1.00;
        for ( i = 0; i < d; ++i )
            fx[j] *= exp ( - pow ( x[i+j*d] / 2.00, 2 ) / 2.00 )/ 2.00 / sqrt ( M_2TPI );
    }

    return fx;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _get_seq ( void * data)
/******************************************************************************/
/*
  Purpose:
    GET_SEQ generates all positive integer D-vectors that sum to NORM.
  Discussion:
    This function computes a list, in reverse dictionary order, of
    all D-vectors of positive values that sum to NORM.
    For example, call get_seq ( 3, 5, 6, fs ) returns
      3  1  1
      2  2  1
      2  1  2
      1  3  1
      1  2  2
      1  1  3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    Original MATLAB version by Florian Heiss, Viktor Winschel.
    C version by John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
  Parameters:
    Input, int D, the dimension.
    1 <= D.
    Input, int NORM, the value that each row must sum to.
    D <= NORM.
    Input, int SEQ_NUM, the number of rows of FS.
    Output, int GET_SEQ[SEQ_NUM*D].  Each row of FS represents
    one vector with all elements positive and summing to NORM.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ d = a_data[0];
	const register dim_typ norm = a_data[1];
	const register dim_typ seq_num = a_data[2];
	
    dim_typ a;
    dim_typ c;
    int *fs;
    dim_typ i;
    dim_typ row;
    int *seq;

    seq = ( int * ) malloc ( d * sizeof ( int ) );
    fs = ( int * ) malloc ( seq_num * d * sizeof ( int ) );

    if ( norm < d )
        return NULL;

    for ( i = 0; i < d; ++i )
        seq[i] = 0;
    /*
    The algorithm is written to work with vectors whose minimum value is
    allowed to be zero.  So we subtract D from NORM at the beginning and
    then increment the result vectors by 1 at the end!
    */
    a = norm - d;
    seq[0] = a;

    row = 0;
    for ( i = 0; i < d; ++i )
        fs[row+i*seq_num] = seq[i] + 1;
    c = 0;

    while ( seq[d-1] < a )
    {
        if ( c == d - 1 )
        {
            for ( i = c - 1; 0 <= i; --i)
            {
                c = i;
                if ( seq[i] != 0 )
                    break;
            }
        }

        -- seq[c];
        ++ c;
        seq[c] = a;
        for ( i = 0; i < c; ++i )
            seq[c] -= seq[i];

        if ( c < d - 1 )
            for ( i = c + 1; i < d; ++i )
                seq[i] = 0;

        ++ row;
        for ( i = 0; i < d; ++i )
            fs[row+i*seq_num] = seq[i] + 1;
    }

    free ( seq );
    return fs;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gqn ( void * data)
/******************************************************************************/
/*
  Purpose:
    GQN provides data for Gauss quadrature with a normal weight.
  Discussion:
    This data assumes integration over the interval (-oo,+oo) with
    weight function w(x) = exp(-x*x/2)/sqrt(2*M_PI).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 December 2012
  Author:
    Original MATLAB version by Florian Heiss, Viktor Winschel.
    C version by John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
  Parameters:
    Input, int N, the number of points and weights.
    1 <= N <= 25.
    Output, double X[N], the nodes.
    Output, double W[N], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
    ityp x01[1] =
    {
        0.00000000000000000
    };
    ityp w01[1] =
    {
        1.00000000000000000
    };
    ityp x02[2] =
    {
        1.00000000000000000,
        1.00000000000000000
    };
    ityp w02[2] =
    {
        0.50000000000000000,
        0.50000000000000000
    };
    ityp x03[3] =
    {
        -1.73205080756887719,
        0.00000000000000000,
        1.73205080756887719
    };
    ityp w03[3] =
    {
        0.166666666666666741,
         0.66666666666666663,
         0.166666666666666741
    };
    ityp x04[4] =
    {
        -2.33441421833897733,
        -0.741963784302725915,
        0.741963784302725915,
        2.33441421833897733
    };
    ityp w04[4] =
    {
        0.0458758547680684983,
        0.454124145231931453,
        0.454124145231931453,
        0.0458758547680684983
    };
    ityp x05[5] =
    {
        -2.85697001387280558,
        -1.35562617997426593,
        0.00000000000000000,
        1.35562617997426593,
        2.85697001387280558
    };
    ityp w05[5] =
    {
        0.011257411327720691,
        0.22207592200561263,
        0.533333333333333437,
        0.22207592200561263,
        0.011257411327720691
    };
    ityp x06[6] =
    {
        -3.32425743355211933,
        -1.88917587775371087,
        -0.616706590192594217,
        0.616706590192594217,
        1.88917587775371087,
        3.32425743355211933
    };
    ityp w06[6] =
    {
        0.00255578440205624308,
        0.0886157460419145226,
        0.408828469556029195,
        0.408828469556029195,
        0.0886157460419145226,
        0.00255578440205624308
    };
    ityp x07[7] =
    {
        -3.75043971772574247,
        -2.36675941073454155,
        -1.15440539473996817,
        0.00000000000000000,
        1.15440539473996817,
        2.36675941073454155,
        3.75043971772574247
    };
    ityp w07[7] =
    {
        0.000548268855972218754,
        0.0307571239675864909,
        0.240123178605012505,
        0.457142857142857573,
        0.240123178605012505,
        0.0307571239675864909,
        0.000548268855972218754
    };
    ityp x08[8] =
    {
        -4.14454718612589446,
        -2.80248586128754162,
        -1.63651904243510815,
        -0.539079811351375171,
        0.539079811351375171,
        1.63651904243510815,
        2.80248586128754162,
        4.14454718612589446
    };
    ityp w08[8] =
    {
        0.000112614538375367836,
        0.00963522012078826297,
        0.117239907661758971,
        0.373012257679077364,
        0.373012257679077364,
        0.117239907661758971,
        0.00963522012078826297,
        0.000112614538375367836
    };
    ityp x09[9] =
    {
        -4.51274586339978256,
        -3.20542900285647026,
        -2.07684797867783022,
        -1.02325566378913257,
        0.00000000000000000,
        1.02325566378913257,
        2.07684797867783022,
        3.20542900285647026,
        4.51274586339978256
    };
    ityp w09[9] =
    {
        2.23458440077465626E-05,
        0.0027891413212317675,
        0.0499164067652179688,
        0.244097502894939089,
        0.406349206349206848,
        0.244097502894939089,
        0.0499164067652179688,
        0.0027891413212317675,
        2.23458440077465626E-05
    };
    ityp x10[10] =
    {
        -4.85946282833231269,
        -3.58182348355192692,
        -2.48432584163895465,
        -1.46598909439115821,
        -0.484935707515497638,
        0.484935707515497638,
        1.46598909439115821,
        2.48432584163895465,
        3.58182348355192692,
        4.85946282833231269
    };
    ityp w10[10] =
    {
        4.3106526307183106E-06,
        0.000758070934312219725,
        0.0191115805007703171,
        0.135483702980267295,
        0.344642334932019401,
        0.344642334932019401,
        0.135483702980267295,
        0.0191115805007703171,
        0.000758070934312219725,
        4.3106526307183106E-06
    };
    ityp x11[11] =
    {
        -5.18800122437487143,
        -3.93616660712997746,
        -2.8651231606436447,
        -1.87603502015484591,
        -0.928868997381063877,
        0.00000000000000000,
        0.928868997381063877,
        1.87603502015484591,
        2.8651231606436447,
        3.93616660712997746,
        5.18800122437487143
    };
    ityp w11[11] =
    {
        8.12184979021490357E-07,
        0.000195671930271223244,
        0.0067202852355372697,
        0.0661387460710576441,
        0.242240299873970027,
        0.369408369408369575,
        0.242240299873970027,
        0.0661387460710576441,
        0.0067202852355372697,
        0.000195671930271223244,
        8.12184979021490357E-07
    };
    ityp x12[12] =
    {
        -5.50090170446774795,
        -4.27182584793228148,
        -3.22370982877009737,
        -2.25946445100079929,
        -1.34037519715161668,
        -0.444403001944139009,
        0.444403001944139009,
        1.34037519715161668,
        2.25946445100079929,
        3.22370982877009737,
        4.27182584793228148,
        5.50090170446774795
    };
    ityp w12[12] =
    {
        1.49992716763715968E-07, 4.83718492259060763E-05, 0.00220338068753318491,
        0.0291166879123641378,   0.146967048045329951,   0.321664361512830066,
        0.321664361512830066,   0.146967048045329951,  0.0291166879123641378,
        0.00220338068753318491, 4.83718492259060763E-05, 1.49992716763715968E-07
    };
    ityp x13[13] =
    {
        -5.8001672523865011,   -4.59139844893652072,   -3.56344438028163468,
        -2.62068997343221488,    -1.7254183795882394,  -0.856679493519450053,
        0.00000000000000000,   0.856679493519450053,     1.7254183795882394,
        2.62068997343221488,    3.56344438028163468,    4.59139844893652072,
        5.8001672523865011
    };
    ityp w13[13] =
    {
        2.72262764280590389E-08, 1.15265965273338848E-05, 0.000681236350442926191,
        0.0117705605059965426,  0.0791689558604501409,   0.237871522964135884,
        0.340992340992341492,   0.237871522964135884,  0.0791689558604501409,
        0.0117705605059965426, 0.000681236350442926191, 1.15265965273338848E-05,
        2.72262764280590389E-08
    };
    ityp x14[14] =
    {
        -6.08740954690129144,   -4.89693639734556463,   -3.88692457505976963,
        -2.96303657983866753,   -2.08834474570194439,   -1.24268895548546432,
        -0.412590457954601808,   0.412590457954601808,    1.24268895548546432,
        2.08834474570194439,    2.96303657983866753,    3.88692457505976963,
        4.89693639734556463,    6.08740954690129144
    };
    ityp w14[14] =
    {
        4.86816125774838718E-09, 2.66099134406763342E-06, 0.00020033955376074381,
        0.00442891910694740657,  0.0386501088242534319,   0.154083339842513656,
        0.302634626813019447,   0.302634626813019447,   0.154083339842513656,
        0.0386501088242534319, 0.00442891910694740657, 0.00020033955376074381,
        2.66099134406763342E-06, 4.86816125774838718E-09
    };
    ityp x15[15] =
    {
        -6.36394788882983775,   -5.19009359130478209,   -4.19620771126901548,
        -3.28908242439876641,   -2.43243682700975805,   -1.60671006902873015,
        -0.799129068324548109,    0.00000000000000000,   0.799129068324548109,
        1.60671006902873015,    2.43243682700975805,    3.28908242439876641,
        4.19620771126901548,    5.19009359130478209,    6.36394788882983775
    };
    ityp w15[15] =
    {
        8.58964989963318053E-10, 5.97541959792059611E-07, 5.64214640518901565E-05,
        0.00156735750354995707,  0.0173657744921376159,  0.0894177953998444436,
        0.232462293609732223,   0.318259518259518204,   0.232462293609732223,
        0.0894177953998444436,  0.0173657744921376159, 0.00156735750354995707,
        5.64214640518901565E-05, 5.97541959792059611E-07, 8.58964989963318053E-10
    };
    ityp x16[16] =
    {
        -6.63087819839312953,   -5.47222570594934332,   -4.49295530252001196,
        -3.60087362417154866,   -2.76024504763070189,   -1.95198034571633361,
        -1.1638291005549648,  -0.386760604500557381,   0.386760604500557381,
        1.1638291005549648,    1.95198034571633361,    2.76024504763070189,
        3.60087362417154866,    4.49295530252001196,    5.47222570594934332,
        6.63087819839312953
    };
    ityp w16[16] =
    {
        1.49781472316183141E-10, 1.30947321628682029E-07, 1.53000321624872858E-05,
        0.000525984926573909786,  0.0072669376011847411,  0.0472847523540140674,
        0.158338372750949252,   0.286568521238012408,   0.286568521238012408,
        0.158338372750949252,  0.0472847523540140674,  0.0072669376011847411,
        0.000525984926573909786, 1.53000321624872858E-05, 1.30947321628682029E-07,
        1.49781472316183141E-10
    };
    ityp x17[17] =
    {
        -6.88912243989533302,   -5.74446007865940711,   -4.77853158962998403,
        -3.90006571719801043,   -3.07379717532819408,   -2.28101944025298886,
        -1.50988330779674085,  -0.751842600703896302,    0.00000000000000000,
        0.751842600703896302,    1.50988330779674085,    2.28101944025298886,
        3.07379717532819408,    3.90006571719801043,    4.77853158962998403,
        5.74446007865940711,    6.88912243989533302
    };
    ityp w17[17] =
    {
        2.58431491937491514E-11, 2.80801611793057831E-08, 4.0126794479798725E-06,
        0.000168491431551339447, 0.00285894606228464989,   0.023086657025711152,
        0.0974063711627180806,   0.226706308468978768,   0.299538370126607556,
        0.226706308468978768,  0.0974063711627180806,   0.023086657025711152,
        0.00285894606228464989, 0.000168491431551339447, 4.0126794479798725E-06,
        2.80801611793057831E-08, 2.58431491937491514E-11
    };
    ityp x18[18] =
    {
        -7.13946484914647961,   -6.00774591135959746,   -5.05407268544274046,
        -4.1880202316294044,   -3.37473653577809074,   -2.59583368891124033,
        -1.83977992150864567,   -1.09839551809150127,  -0.365245755507697667,
        0.365245755507697667,    1.09839551809150127,    1.83977992150864567,
        2.59583368891124033,    3.37473653577809074,     4.1880202316294044,
        5.05407268544274046,    6.00774591135959746,    7.13946484914647961
    };
    ityp w18[18] =
    {
        4.41658876935870775E-12, 5.90548847883654844E-09, 1.02155239763698159E-06,
        5.17989614411619621E-05, 0.00106548479629164959,  0.0105165177519413525,
        0.0548966324802226541,   0.160685303893512627,   0.272783234654287887,
        0.272783234654287887,   0.160685303893512627,  0.0548966324802226541,
        0.0105165177519413525, 0.00106548479629164959, 5.17989614411619621E-05,
        1.02155239763698159E-06, 5.90548847883654844E-09, 4.41658876935870775E-12
    };
    ityp x19[19] =
    {
        -7.38257902403043165,    -6.2628911565132519,   -5.32053637733603857,
        -4.46587262683103159,   -3.66441654745063827,   -2.89805127651575356,
        -2.15550276131693508,    -1.4288766760783731,  -0.712085044042379933,
        0.00000000000000000,   0.712085044042379933,     1.4288766760783731,
        2.15550276131693508,    2.89805127651575356,    3.66441654745063827,
        4.46587262683103159,    5.32053637733603857,     6.2628911565132519,
        7.38257902403043165
    };
    ityp w19[19] =
    {
        7.4828300540572308E-13, 1.22037084844747862E-09, 2.53222003209286807E-07,
        1.53511459546667444E-05, 0.000378502109414267593, 0.00450723542034203555,
        0.0286666910301184956,   0.103603657276143998,   0.220941712199143658,
        0.283773192751521075,   0.220941712199143658,   0.103603657276143998,
        0.0286666910301184956, 0.00450723542034203555, 0.000378502109414267593,
        1.53511459546667444E-05, 2.53222003209286807E-07, 1.22037084844747862E-09,
        7.4828300540572308E-13
    };
    ityp x20[20] =
    {
        -7.61904854167975909,   -6.51059015701365507,   -5.57873880589320148,
        -4.73458133404605519,    -3.9439673506573163,   -3.18901481655339003,
        -2.45866361117236787,   -1.74524732081412703,   -1.04294534880275092,
        -0.346964157081355917,   0.346964157081355917,    1.04294534880275092,
        1.74524732081412703,    2.45866361117236787,    3.18901481655339003,
        3.9439673506573163,    4.73458133404605519,    5.57873880589320148,
        6.51059015701365507,    7.61904854167975909
    };
    ityp w20[20] =
    {
        1.25780067243793047E-13, 2.4820623623151838E-10, 6.12749025998295973E-08,
        4.40212109023086457E-06, 0.000128826279961928981, 0.00183010313108049175,
        0.0139978374471010428,  0.0615063720639760295,   0.161739333984000255,
        0.26079306344955544,    0.26079306344955544,   0.161739333984000255,
        0.0615063720639760295,  0.0139978374471010428, 0.00183010313108049175,
        0.000128826279961928981, 4.40212109023086457E-06, 6.12749025998295973E-08,
        2.4820623623151838E-10, 1.25780067243793047E-13
    };
    ityp x21[21] =
    {
        -7.84938289511382248,   -6.75144471871746088,   -5.82938200730447065,
        -4.99496394478202532,   -4.21434398168842161,   -3.46984669047537642,
        -2.75059298105237326,    -2.0491024682571628,   -1.35976582321123041,
        -0.678045692440644054,     0.00000000000000000,   0.678045692440644054,
        1.35976582321123041,     2.0491024682571628,    2.75059298105237326,
        3.46984669047537642,    4.21434398168842161,    4.99496394478202532,
        5.82938200730447065,    6.75144471871746088,    7.84938289511382248
    };
    ityp w21[21] =
    {
        2.09899121956566525E-14, 4.97536860412174643E-11, 1.45066128449307397E-08,
        1.22535483614825217E-06, 4.21923474255158655E-05, 0.000708047795481537364,
        0.00643969705140877684,  0.0339527297865428387,    0.10839228562641938,
        0.215333715695059824,   0.270260183572877066,   0.215333715695059824,
        0.10839228562641938,  0.0339527297865428387, 0.00643969705140877684,
        0.000708047795481537364, 4.21923474255158655E-05, 1.22535483614825217E-06,
        1.45066128449307397E-08, 4.97536860412174643E-11, 2.09899121956566525E-14
    };
    ityp x22[22] =
    {
        -8.07402998402171157,   -6.98598042401881525,   -6.07307495112289786,
        -5.2477244337144251,   -4.47636197731086849,   -3.74149635026651772,
        -3.03240422783167629,   -2.34175999628770803,    -1.6641248391179071,
        -0.995162422271215541,  -0.331179315715273814,   0.331179315715273814,
        0.995162422271215541,     1.6641248391179071,    2.34175999628770803,
        3.03240422783167629,    3.74149635026651772,    4.47636197731086849,
        5.2477244337144251,    6.07307495112289786,    6.98598042401881525,
        8.07402998402171157
    };
    ityp w22[22] =
    {
        3.47946064787714279E-15, 9.84137898234601051E-12, 3.36651415945821088E-09,
        3.31985374981400429E-07, 1.33459771268087124E-05, 0.000262283303255964159,
        0.00280876104757721073,  0.0175690728808057736,  0.0671963114288898905,
        0.161906293413675378,   0.250243596586935013,   0.250243596586935013,
        0.161906293413675378,  0.0671963114288898905,  0.0175690728808057736,
        0.00280876104757721073, 0.000262283303255964159, 1.33459771268087124E-05,
        3.31985374981400429E-07, 3.36651415945821088E-09, 9.84137898234601051E-12,
        3.47946064787714279E-15
    };
    ityp x23[23] =
    {
        -8.29338602741735365,   -7.21465943505186225,   -6.31034985444839958,
        -5.49347398647179475,   -4.73072419745147332,   -4.00477532173330442,
        -3.30504002175296518,   -2.62432363405918201,    -1.9573275529334242,
        -1.29987646830397896,  -0.648471153534495803,     0.00000000000000000,
        0.648471153534495803,    1.29987646830397896,     1.9573275529334242,
        2.62432363405918201,    3.30504002175296518,    4.00477532173330442,
        4.73072419745147332,    5.49347398647179475,    6.31034985444839958,
        7.21465943505186225,    8.29338602741735365
    };
    ityp w23[23] =
    {
        5.73238316780208728E-16, 1.92293531156779128E-12, 7.67088886239990765E-10,
        8.77506248386171607E-08, 4.08997724499215494E-06, 9.34081860903129835E-05,
        0.00116762863749786134, 0.00857967839146566401,  0.0388671837034809467,
        0.112073382602620911,   0.209959669577542613,   0.258509740808839039,
        0.209959669577542613,   0.112073382602620911,  0.0388671837034809467,
        0.00857967839146566401, 0.00116762863749786134, 9.34081860903129835E-05,
        4.08997724499215494E-06, 8.77506248386171607E-08, 7.67088886239990765E-10,
        1.92293531156779128E-12, 5.73238316780208728E-16
    };
    ityp x24[24] =
    {
        -8.50780351919525835,   -7.43789066602166304,   -6.54167500509863409,
        -5.73274717525120092,   -4.97804137463912078,   -4.26038360501990532,
        -3.56930676407356096,   -2.89772864322331403,   -2.24046785169175244,
        -1.59348042981642024,  -0.953421922932109256,  -0.317370096629452314,
        0.317370096629452314,   0.953421922932109256,    1.59348042981642024,
        2.24046785169175244,    2.89772864322331403,    3.56930676407356096,
        4.26038360501990532,    4.97804137463912078,    5.73274717525120092,
        6.54167500509863409,    7.43789066602166304,    8.50780351919525835
    };
    ityp w24[24] =
    {
        9.39019368904192022E-17, 3.71497415276241595E-13, 1.71866492796486901E-10,
        2.26746167348046514E-08, 1.21765974544258296E-06, 3.20950056527459886E-05,
        0.000464718718779397633, 0.00397660892918131129,  0.0211263444089670287,
        0.0720693640171784361,   0.161459512867000249,   0.240870115546640562,
        0.240870115546640562,   0.161459512867000249,  0.0720693640171784361,
        0.0211263444089670287, 0.00397660892918131129, 0.000464718718779397633,
        3.20950056527459886E-05, 1.21765974544258296E-06, 2.26746167348046514E-08,
        1.71866492796486901E-10, 3.71497415276241595E-13, 9.39019368904192022E-17
    };
    ityp x25[25] =
    {
        -8.71759767839958855,   -7.65603795539307619,   -6.76746496380971685,
        -5.96601469060670198,   -5.21884809364427937,   -4.50892992296728501,
        -3.82590056997249173,   -3.16277567938819271,   -2.51447330395220581,
        -1.8770583699478387,   -1.24731197561678919,  -0.622462279186076106,
        0.00000000000000000,   0.622462279186076106,    1.24731197561678919,
        1.8770583699478387,    2.51447330395220581,    3.16277567938819271,
        3.82590056997249173,    4.50892992296728501,    5.21884809364427937,
        5.96601469060670198,    6.76746496380971685,    7.65603795539307619,
        8.71759767839958855
    };
    ityp w25[25] =
    {
        1.53003899799868247E-17, 7.10210303700392527E-14, 3.79115000047718706E-11,
        5.7380238688993763E-09, 3.53015256024549785E-07, 1.06721949052025363E-05,
        0.0001777669069265266, 0.00175785040526379608,  0.0108567559914623159,
        0.0433799701676449712,   0.114880924303951637,   0.204851025650340413,
        0.248169351176485475,   0.204851025650340413,   0.114880924303951637,
        0.0433799701676449712,  0.0108567559914623159, 0.00175785040526379608,
        0.0001777669069265266, 1.06721949052025363E-05, 3.53015256024549785E-07,
        5.7380238688993763E-09, 3.79115000047718706E-11, 7.10210303700392527E-14,
        1.53003899799868247E-17
    };

    if ( n == 1 )
    {
        r8vec_copy ( n, x01, x );
        r8vec_copy ( n, w01, w );
    }
    else if ( n == 2 )
    {
        r8vec_copy ( n, x02, x );
        r8vec_copy ( n, w02, w );
    }
    else if ( n == 3 )
    {
        r8vec_copy ( n, x03, x );
        r8vec_copy ( n, w03, w );
    }
    else if ( n == 4 )
    {
        r8vec_copy ( n, x04, x );
        r8vec_copy ( n, w04, w );
    }
    else if ( n == 5 )
    {
        r8vec_copy ( n, x05, x );
        r8vec_copy ( n, w05, w );
    }
    else if ( n == 6 )
    {
        r8vec_copy ( n, x06, x );
        r8vec_copy ( n, w06, w );
    }
    else if ( n == 7 )
    {
        r8vec_copy ( n, x07, x );
        r8vec_copy ( n, w07, w );
    }
    else if ( n == 8 )
    {
        r8vec_copy ( n, x08, x );
        r8vec_copy ( n, w08, w );
    }
    else if ( n == 9 )
    {
        r8vec_copy ( n, x09, x );
        r8vec_copy ( n, w09, w );
    }
    else if ( n == 10 )
    {
        r8vec_copy ( n, x10, x );
        r8vec_copy ( n, w10, w );
    }
    else if ( n == 11 )
    {
        r8vec_copy ( n, x11, x );
        r8vec_copy ( n, w11, w );
    }
    else if ( n == 12 )
    {
        r8vec_copy ( n, x12, x );
        r8vec_copy ( n, w12, w );
    }
    else if ( n == 13 )
    {
        r8vec_copy ( n, x13, x );
        r8vec_copy ( n, w13, w );
    }
    else if ( n == 14 )
    {
        r8vec_copy ( n, x14, x );
        r8vec_copy ( n, w14, w );
    }
    else if ( n == 15 )
    {
        r8vec_copy ( n, x15, x );
        r8vec_copy ( n, w15, w );
    }
    else if ( n == 16 )
    {
        r8vec_copy ( n, x16, x );
        r8vec_copy ( n, w16, w );
    }
    else if ( n == 17 )
    {
        r8vec_copy ( n, x17, x );
        r8vec_copy ( n, w17, w );
    }
    else if ( n == 18 )
    {
        r8vec_copy ( n, x18, x );
        r8vec_copy ( n, w18, w );
    }
    else if ( n == 19 )
    {
        r8vec_copy ( n, x19, x );
        r8vec_copy ( n, w19, w );
    }
    else if ( n == 20 )
    {
        r8vec_copy ( n, x20, x );
        r8vec_copy ( n, w20, w );
    }
    else if ( n == 21 )
    {
        r8vec_copy ( n, x21, x );
        r8vec_copy ( n, w21, w );
    }
    else if ( n == 22 )
    {
        r8vec_copy ( n, x22, x );
        r8vec_copy ( n, w22, w );
    }
    else if ( n == 23 )
    {
        r8vec_copy ( n, x23, x );
        r8vec_copy ( n, w23, w );
    }
    else if ( n == 24 )
    {
        r8vec_copy ( n, x24, x );
        r8vec_copy ( n, w24, w );
    }
    else if ( n == 25 )
    {
        r8vec_copy ( n, x25, x );
        r8vec_copy ( n, w25, w );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gqn_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    GQN_ORDER computes the order of a GQN rule from the level.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int L, the level of the rule.
    1 <= L.
    Output, int GQN_ORDER, the order of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ l = *(dim_typ *) data;
	
	result = l<1 || l> 25 ? USHRT_MAX : l;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gqn2_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    GQN2_ORDER computes the order of a GQN rule from the level.
  Discussion:
    For this version of the order routine, we have
      n = 2 * l - 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 February 2014
  Author:
    John Burkardt.
  Parameters:
    Input, int L, the level of the rule.
    1 <= L.
    Output, int GQN2_ORDER, the order of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ l = *(dim_typ *) data;
	
	result = l<1 || l>13 ? USHRT_MAX : (l<<1)-1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gqu ( void * data)
/******************************************************************************/
/*
  Purpose:
    GQU provides data for Gauss quadrature with a uniform weight.
  Discussion:
    This data assumes integration over the interval [0,1] with
    weight function w(x) = 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 December 2012
  Author:
    John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
  Parameters:
    Input, int N, the number of points and weights.
    1 <= N <= 25.
    Output, double X[N], the nodes.
    Output, double W[N], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
    ityp x01[1] =
    {
        0.500000000000000000
    };
    ityp w01[1] =
    {
        1.000000000000000000
    };
    ityp x02[2] =
    {
        0.211324865405187134,  0.788675134594812866
    };
    ityp w02[2] =
    {
        0.500000000000000000,  0.500000000000000000
    };
    ityp x03[3] =
    {
        0.112701665379258298,  0.500000000000000000,  0.887298334620741702
    };
    ityp w03[3] =
    {
        0.277777777777777124,  0.444444444444445697,  0.277777777777777124
    };
    ityp x04[4] =
    {
        0.0694318442029737692,  0.330009478207571871,  0.669990521792428129,
        0.930568155797026231
    };
    ityp w04[4] =
    {
        0.173927422568724843,  0.326072577431275157,  0.326072577431275157,
        0.173927422568724843
    };
    ityp x05[5] =
    {
        0.0469100770306680737,  0.230765344947158502, 0.500000000000000000,
        0.769234655052841498,  0.953089922969331926
    };
    ityp w05[5] =
    {
        0.118463442528091739,  0.239314335249685012,  0.284444444444446554,
        0.239314335249685012,  0.118463442528091739
    };
    ityp x06[6] =
    {
        0.033765242898423975,  0.169395306766867648,  0.380690406958401506,
        0.619309593041598494,  0.830604693233132352,  0.966234757101576025
    };
    ityp w06[6] =
    {
        0.0856622461895818338,  0.180380786524070719,  0.233956967286347461,
        0.233956967286347461,  0.180380786524070719, 0.0856622461895818338
    };
    ityp x07[7] =
    {
        0.0254460438286208124,   0.12923440720030277,  0.297077424311301463,
        0.500000000000000000,  0.702922575688698537,   0.87076559279969723,
        0.974553956171379188
    };
    ityp w07[7] =
    {
        0.0647424830844317012,  0.139852695744639349,  0.190915025252560905,
        0.2089795918367362,  0.190915025252560905,  0.139852695744639349,
        0.0647424830844317012
    };
    ityp x08[8] =
    {
        0.0198550717512319119,  0.101666761293186525,  0.237233795041835505,
        0.408282678752175054,  0.591717321247824946,  0.762766204958164495,
        0.898333238706813475,  0.980144928248768088
    };
    ityp w08[8] =
    {
        0.0506142681451851803,  0.111190517226687935,  0.156853322938944689,
        0.181341891689182133,  0.181341891689182133,  0.156853322938944689,
        0.111190517226687935, 0.0506142681451851803
    };
    ityp x09[9] =
    {
        0.0159198802461869571, 0.0819844463366821152,  0.193314283649704821,
        0.337873288298095487,  0.500000000000000000,  0.662126711701904513,
        0.806685716350295179,  0.918015553663317885,  0.984080119753813043
    };
    ityp w09[9] =
    {
        0.0406371941807845832, 0.0903240803474292531,  0.130305348201468441,
        0.156173538520002264,  0.165119677500630752,  0.156173538520002264,
        0.130305348201468441, 0.0903240803474292531, 0.0406371941807845832
    };
    ityp x10[10] =
    {
        0.0130467357414141283, 0.0674683166555076763,  0.160295215850487782,
        0.283302302935376393,   0.42556283050918442,   0.57443716949081558,
        0.716697697064623607,  0.839704784149512218,  0.932531683344492324,
        0.986953264258585872
    };
    ityp w10[10] =
    {
        0.0333356721543420012, 0.0747256745752905988,  0.109543181257991576,
        0.13463335965499873,  0.147762112357377129,  0.147762112357377129,
        0.13463335965499873,  0.109543181257991576, 0.0747256745752905988,
        0.0333356721543420012
    };
    ityp x11[11] =
    {
        0.010885670926971569, 0.0564687001159522861,  0.134923997212975322,
        0.240451935396594152,  0.365228422023827548,  0.500000000000000000,
        0.634771577976172452,  0.759548064603405848,  0.865076002787024678,
        0.943531299884047714,  0.989114329073028431
    };
    ityp w11[11] =
        {
        0.0278342835580849164, 0.0627901847324526252, 0.0931451054638675197,
        0.11659688229599563,  0.131402272255123881,  0.136462543388950863,
        0.131402272255123881,   0.11659688229599563, 0.0931451054638675197,
        0.0627901847324526252, 0.0278342835580849164
    };
    ityp x12[12] =
    {
        0.00921968287664043373, 0.0479413718147625456,  0.115048662902847654,
        0.206341022856691314,  0.316084250500909936,  0.437383295744265488,
        0.562616704255734512,  0.683915749499090064,  0.793658977143308686,
        0.884951337097152346,  0.952058628185237454,  0.990780317123359566
    };
    ityp w12[12] =
    {
        0.0235876681932543145, 0.0534696629976592758, 0.0800391642716734436,
        0.101583713361533282,  0.116746268269177805,  0.124573522906701886,
        0.124573522906701886,  0.116746268269177805,  0.101583713361533282,
        0.0800391642716734436, 0.0534696629976592758, 0.0235876681932543145
    };
    ityp x13[13] =
    {
        0.00790847264070593248, 0.0412008003885109275, 0.0992109546333450609,
        0.178825330279829942,  0.275753624481776649,  0.384770842022432613,
        0.500000000000000000,  0.615229157977567387,  0.724246375518223351,
        0.821174669720170058,  0.900789045366654939,  0.958799199611489072,
        0.992091527359294068
    };
    ityp w13[13] =
    {
        0.0202420023826562281, 0.0460607499188643785, 0.0694367551098938746,
        0.0890729903809732021,  0.103908023768444616,  0.113141590131449032,
        0.116275776615437407,  0.113141590131449032,  0.103908023768444616,
        0.0890729903809732021, 0.0694367551098938746, 0.0460607499188643785,
        0.0202420023826562281
    };
    ityp x14[14] =
    {
        0.00685809565159378742, 0.0357825581682131855, 0.0863993424651174902,
        0.156353547594157316,   0.24237568182092295,  0.340443815536055183,
        0.445972525646328166,  0.554027474353671834,  0.659556184463944817,
        0.75762431817907705,  0.843646452405842684,   0.91360065753488251,
        0.964217441831786815,  0.993141904348406213
    };
    ityp w14[14] =
    {
        0.0175597301658745736, 0.0400790435798802913, 0.0607592853439517105,
        0.0786015835790969952, 0.0927691987389691608,  0.102599231860648107,
        0.107631926731579161,  0.107631926731579161,  0.102599231860648107,
        0.0927691987389691608, 0.0786015835790969952, 0.0607592853439517105,
        0.0400790435798802913, 0.0175597301658745736
    };
    ityp x15[15] =
    {
        0.0060037409897573113, 0.0313633037996470243, 0.0758967082947863414,
        0.137791134319914965,  0.214513913695730585,  0.302924326461218252,
        0.399402953001282701,  0.500000000000000000,  0.600597046998717299,
        0.697075673538781748,  0.785486086304269415,  0.862208865680085035,
        0.924103291705213659,  0.968636696200352976,  0.993996259010242689
    };
    ityp w15[15] =
    {
        0.0153766209980574341, 0.0351830237440541593, 0.0535796102335861571,
        0.0697853389630773147,  0.083134602908497196, 0.0930805000077812861,
        0.0992157426635560391,  0.101289120962780907, 0.0992157426635560391,
        0.0930805000077812861,  0.083134602908497196, 0.0697853389630773147,
        0.0535796102335861571, 0.0351830237440541593, 0.0153766209980574341
    };
    ityp x16[16] =
    {
        0.00529953250417503074, 0.0277124884633836999, 0.0671843988060840669,
        0.122297795822498445,  0.191061877798678115,  0.270991611171386371,
        0.359198224610370542,  0.452493745081181231,  0.547506254918818769,
        0.640801775389629458,  0.729008388828613629,  0.808938122201321885,
        0.877702204177501555,  0.932815601193915933,    0.9722875115366163,
        0.994700467495824969
    };
    ityp w16[16] =
    {
        0.0135762297058759553, 0.0311267619693239538, 0.0475792558412465455,
        0.062314485627767105, 0.0747979944082885623, 0.0845782596975014622,
        0.0913017075224620001, 0.0947253052275344315, 0.0947253052275344315,
        0.0913017075224620001, 0.0845782596975014622, 0.0747979944082885623,
        0.062314485627767105, 0.0475792558412465455, 0.0311267619693239538,
        0.0135762297058759553
    };
    ityp x17[17] =
    {
        0.00471226234279131795,  0.024662239115616158, 0.0598804231365070994,
        0.109242998051599316,  0.171164420391654692,  0.243654731456761531,
        0.324384118273061794,  0.410757909252076114,  0.500000000000000000,
        0.589242090747923886,  0.675615881726938206,  0.756345268543238469,
        0.828835579608345308,  0.890757001948400684,  0.940119576863492901,
        0.975337760884383842,  0.995287737657208682
    };
    ityp w17[17] =
    {
        0.01207415143427314, 0.0277297646869936118, 0.0425180741585896443,
        0.0559419235967020534, 0.0675681842342628902, 0.0770228805384053083,
        0.0840020510782251428, 0.0882813526834964474, 0.0897232351781034193,
        0.0882813526834964474, 0.0840020510782251428, 0.0770228805384053083,
        0.0675681842342628902, 0.0559419235967020534, 0.0425180741585896443,
        0.0277297646869936118,   0.01207415143427314
    };
    ityp x18[18] =
    {
        0.00421741578953449547,  0.022088025214301199, 0.0536987667512220934,
        0.0981475205137384288,  0.154156478469823388,   0.22011458446302623,
        0.294124419268578685,  0.374056887154247231,  0.457612493479132354,
        0.542387506520867646,  0.625943112845752769,  0.705875580731421315,
        0.77988541553697377,  0.845843521530176612,  0.901852479486261571,
        0.946301233248777907,  0.977911974785698801,  0.995782584210465505
    };
    ityp w18[18] =
    {
        0.0108080067632407191, 0.0248572744474849679, 0.0382128651274446646,
        0.0504710220531437159,  0.061277603355739306, 0.0703214573353254518,
        0.0773423375631328014, 0.0821382418729165037,  0.084571191481571939,
        0.084571191481571939, 0.0821382418729165037, 0.0773423375631328014,
        0.0703214573353254518,  0.061277603355739306, 0.0504710220531437159,
        0.0382128651274446646, 0.0248572744474849679, 0.0108080067632407191
    };
    ityp x19[19] =
    {
        0.00379657807820787951, 0.0198959239325849913,  0.048422048192590994,
        0.0886426717314285906,  0.139516911332385307,   0.19972734766915945,
        0.267714629312019503,  0.341717950018185057,  0.419820677179887358,
        0.500000000000000000,  0.580179322820112642,  0.658282049981814943,
        0.732285370687980497,   0.80027265233084055,  0.860483088667614693,
        0.911357328268571409,  0.951577951807409006,  0.980104076067415009,
        0.99620342192179212
    };
    ityp w19[19] =
    {
        0.00973089411486243415, 0.0224071133828498206, 0.0345222713688206687,
        0.0457450108112251244, 0.0557833227736671128,  0.064376981269668232,
        0.0713033510868034126, 0.0763830210329299597, 0.0794844216969773365,
        0.0805272249243919463, 0.0794844216969773365, 0.0763830210329299597,
        0.0713033510868034126,  0.064376981269668232, 0.0557833227736671128,
        0.0457450108112251244, 0.0345222713688206687, 0.0224071133828498206,
        0.00973089411486243415
    };
    ityp x20[20] =
    {
        0.00343570040745255767, 0.0180140363610430398, 0.0438827858743370269,
        0.0804415140888905533,  0.126834046769924602,  0.181973159636742432,
        0.244566499024586381,  0.313146955642290226,  0.386107074429177466,
        0.461736739433251331,  0.538263260566748669,  0.613892925570822534,
        0.686853044357709774,  0.755433500975413619,  0.818026840363257568,
        0.873165953230075398,  0.919558485911109447,  0.956117214125662973,
        0.98198596363895696,  0.996564299592547442
    };
    ityp w20[20] =
    {
        0.0088070035695753026, 0.0203007149001935561, 0.0313360241670545686,
        0.0416383707883524329, 0.0509650599086203179, 0.0590972659807592476,
        0.0658443192245883463, 0.0710480546591911871, 0.0745864932363019956,
        0.0763766935653631129, 0.0763766935653631129, 0.0745864932363019956,
        0.0710480546591911871, 0.0658443192245883463, 0.0590972659807592476,
        0.0509650599086203179, 0.0416383707883524329, 0.0313360241670545686,
        0.0203007149001935561, 0.0088070035695753026
    };
    ityp x21[21] =
    {
        0.00312391468980521836, 0.0163865807168468436, 0.0399503329247996586,
        0.0733183177083414073,  0.115780018262161111,  0.166430597901293886,
        0.224190582056390086,  0.287828939896280556,  0.355989341598799469,
        0.427219072919552412,  0.500000000000000000,  0.572780927080447588,
        0.644010658401200531,  0.712171060103719444,  0.775809417943609914,
        0.833569402098706114,  0.884219981737838889,  0.926681682291658593,
        0.960049667075200341,  0.983613419283153156,  0.996876085310194782
    };
    ityp w21[21] =
    {
        0.00800861412888644909,  0.018476894885426285, 0.0285672127134286406,
        0.0380500568141897075, 0.0467222117280169935, 0.0543986495835743558,
        0.06091570802686435, 0.0661344693166688452, 0.0699436973955366581,
        0.0722622019949851341,  0.073040566824845346, 0.0722622019949851341,
        0.0699436973955366581, 0.0661344693166688452,   0.06091570802686435,
        0.0543986495835743558, 0.0467222117280169935, 0.0380500568141897075,
        0.0285672127134286406,  0.018476894885426285,0.00800861412888644909
    };
    ityp x22[22] =
    {
        0.0028527072588003799, 0.0149697510822857094,  0.036521613906413064,
        0.0670937111398499653,  0.106091597010395944,  0.152756368406658627,
        0.206179798246544199,  0.265322081006621469,  0.329032089553957907,
        0.396069786655889322,  0.465130363340138908,  0.534869636659861092,
        0.603930213344110678,  0.670967910446042093,  0.734677918993378531,
        0.793820201753455801,  0.847243631593341373,  0.893908402989604056,
        0.932906288860150035,  0.963478386093586936,  0.985030248917714291,
        0.99714729274119962
    };
    ityp w22[22] =
    {
        0.00731399764913532799, 0.0168874507924071104, 0.0261466675763416916,
        0.0348982342122602998, 0.0429708031085339753, 0.0502070722214406004,
        0.0564661480402697119, 0.0616261884052562506, 0.0655867523935313168,
        0.0682707491730076971, 0.0696259364278161291, 0.0696259364278161291,
        0.0682707491730076971, 0.0655867523935313168, 0.0616261884052562506,
        0.0564661480402697119, 0.0502070722214406004, 0.0429708031085339753,
        0.0348982342122602998, 0.0261466675763416916, 0.0168874507924071104,
        0.00731399764913532799
    };
    ityp x23[23] =
    {
        0.00261533250122392147,  0.013728764390942394, 0.0335144565869919253,
        0.061623820864779244, 0.0975557991905799948,  0.140669318434024859,
        0.190195062118176939,  0.245249261076996294,  0.304849480984854537,
        0.367932159514827495,  0.433371587850766904,  0.500000000000000000,
        0.566628412149233096,  0.632067840485172505,  0.695150519015145463,
        0.754750738923003706,  0.809804937881823061,  0.859330681565975141,
        0.902444200809420005,  0.938376179135220756,  0.966485543413008075,
        0.986271235609057606,  0.997384667498776079
    };
    ityp w23[23] =
    {
        0.0067059297435702412,  0.015494002928489686, 0.0240188358655423692,
        0.0321162107042629943, 0.0396407058883595509, 0.0464578830300175633,
        0.052446045732270824,  0.057498320111205814, 0.0615245421533648154,
        0.06445286109404115, 0.0662310197023484037, 0.0668272860930531759,
        0.0662310197023484037,   0.06445286109404115, 0.0615245421533648154,
        0.057498320111205814,  0.052446045732270824, 0.0464578830300175633,
        0.0396407058883595509, 0.0321162107042629943, 0.0240188358655423692,
        0.015494002928489686, 0.0067059297435702412
    };
    ityp x24[24] =
    {
        0.00240639000148934468, 0.0126357220143452631, 0.0308627239986336566,
        0.0567922364977995198, 0.0899990070130485265,  0.129937904210722821,
        0.175953174031512227,  0.227289264305580163,  0.283103246186977464,
        0.342478660151918302,  0.404440566263191803,  0.467971553568697241,
        0.532028446431302759,  0.595559433736808197,  0.657521339848081698,
        0.716896753813022536,  0.772710735694419837,  0.824046825968487773,
        0.870062095789277179,  0.910000992986951474,   0.94320776350220048,
        0.969137276001366343,  0.987364277985654737,  0.997593609998510655
    };
    ityp w24[24] =
    {
        0.00617061489999283508,  0.014265694314466934, 0.0221387194087098796,
        0.0296492924577183847, 0.0366732407055402054, 0.0430950807659766927,
        0.0488093260520570393, 0.0537221350579829143, 0.0577528340268628829,
        0.0608352364639017928, 0.0629187281734143178, 0.0639690976733762462,
        0.0639690976733762462, 0.0629187281734143178, 0.0608352364639017928,
        0.0577528340268628829, 0.0537221350579829143, 0.0488093260520570393,
        0.0430950807659766927, 0.0366732407055402054, 0.0296492924577183847,
        0.0221387194087098796,  0.014265694314466934,0.00617061489999283508
    };
    ityp x25[25] =
    {
        0.00222151510475088187, 0.0116680392702412927, 0.0285127143855128384,
        0.052504001060862393, 0.0832786856195830705,  0.120370368481321099,
        0.163216815763265854,  0.211168534879388581,  0.263498634277142485,
        0.319413847095306069,  0.378066558139505737,  0.438567653694644788,
        0.500000000000000000,  0.561432346305355212,  0.621933441860494263,
        0.680586152904693931,  0.736501365722857515,  0.788831465120611419,
        0.836783184236734146,  0.879629631518678901,   0.91672131438041693,
        0.947495998939137607,  0.971487285614487162,  0.988331960729758707,
        0.997778484895249118
    };
    ityp w25[25] =
    {
        0.00569689925051255347, 0.0131774933075161083, 0.0204695783506531476,
        0.0274523479879176906, 0.0340191669061785454, 0.0400703501675005319,
        0.0455141309914819034, 0.0502679745335253628, 0.0542598122371318672,
        0.0574291295728558623, 0.0597278817678924615, 0.0611212214951551217,
        0.061588026863357799, 0.0611212214951551217, 0.0597278817678924615,
        0.0574291295728558623, 0.0542598122371318672, 0.0502679745335253628,
        0.0455141309914819034, 0.0400703501675005319, 0.0340191669061785454,
        0.0274523479879176906, 0.0204695783506531476, 0.0131774933075161083,
        0.00569689925051255347
    };

    if ( n == 1 )
    {
        r8vec_copy ( n, x01, x );
        r8vec_copy ( n, w01, w );
    }
    else if ( n == 2 )
    {
        r8vec_copy ( n, x02, x );
        r8vec_copy ( n, w02, w );
    }
    else if ( n == 3 )
    {
        r8vec_copy ( n, x03, x );
        r8vec_copy ( n, w03, w );
    }
    else if ( n == 4 )
    {
        r8vec_copy ( n, x04, x );
        r8vec_copy ( n, w04, w );
    }
    else if ( n == 5 )
    {
        r8vec_copy ( n, x05, x );
        r8vec_copy ( n, w05, w );
    }
    else if ( n == 6 )
    {
        r8vec_copy ( n, x06, x );
        r8vec_copy ( n, w06, w );
    }
    else if ( n == 7 )
    {
        r8vec_copy ( n, x07, x );
        r8vec_copy ( n, w07, w );
    }
    else if ( n == 8 )
    {
        r8vec_copy ( n, x08, x );
        r8vec_copy ( n, w08, w );
    }
    else if ( n == 9 )
    {
        r8vec_copy ( n, x09, x );
        r8vec_copy ( n, w09, w );
    }
    else if ( n == 10 )
    {
        r8vec_copy ( n, x10, x );
        r8vec_copy ( n, w10, w );
    }
    else if ( n == 11 )
    {
        r8vec_copy ( n, x11, x );
        r8vec_copy ( n, w11, w );
    }
    else if ( n == 12 )
    {
        r8vec_copy ( n, x12, x );
        r8vec_copy ( n, w12, w );
    }
    else if ( n == 13 )
    {
        r8vec_copy ( n, x13, x );
        r8vec_copy ( n, w13, w );
    }
    else if ( n == 14 )
    {
        r8vec_copy ( n, x14, x );
        r8vec_copy ( n, w14, w );
    }
    else if ( n == 15 )
    {
        r8vec_copy ( n, x15, x );
        r8vec_copy ( n, w15, w );
    }
    else if ( n == 16 )
    {
        r8vec_copy ( n, x16, x );
        r8vec_copy ( n, w16, w );
    }
    else if ( n == 17 )
    {
        r8vec_copy ( n, x17, x );
        r8vec_copy ( n, w17, w );
    }
    else if ( n == 18 )
    {
        r8vec_copy ( n, x18, x );
        r8vec_copy ( n, w18, w );
    }
    else if ( n == 19 )
    {
        r8vec_copy ( n, x19, x );
        r8vec_copy ( n, w19, w );
    }
    else if ( n == 20 )
    {
        r8vec_copy ( n, x20, x );
        r8vec_copy ( n, w20, w );
    }
    else if ( n == 21 )
    {
        r8vec_copy ( n, x21, x );
        r8vec_copy ( n, w21, w );
    }
    else if ( n == 22 )
    {
        r8vec_copy ( n, x22, x );
        r8vec_copy ( n, w22, w );
    }
    else if ( n == 23 )
    {
        r8vec_copy ( n, x23, x );
        r8vec_copy ( n, w23, w );
    }
    else if ( n == 24 )
    {
        r8vec_copy ( n, x24, x );
        r8vec_copy ( n, w24, w );
    }
    else if ( n == 25 )
    {
        r8vec_copy ( n, x25, x );
        r8vec_copy ( n, w25, w );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gqu_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    GQU_ORDER computes the order of a GQU rule from the level.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int L, the level of the rule.
    1 <= L <= 25.
    Output, int GQU_ORDER, the order of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ l = *(dim_typ *) data;
	
	result = l<1 || l>25 ? USHRT_MAX:l;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _kpn ( void * data)
/******************************************************************************/
/*
  Purpose:
    KPN provides data for Kronrod-Patterson quadrature with a normal weight.
  Discussion:
    This data assumes integration over the interval (-oo,+oo) with
    weight function w(x) = exp(-x*x/2)/sqrt(2*M_PI).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 December 2012
  Author:
    Original MATLAB version by Florian Heiss, Viktor Winschel.
    C version by John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
    Alan Genz, Bradley Keister,
    Fully symmetric interpolatory rules for multiple integrals
    over infinite regions with Gaussian weight,
    Journal of Computational and Applied Mathematics,
    Volume 71, 1996, pages 299-309.
    Thomas Patterson,
    The optimal addition of points to quadrature formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.
  Parameters:
    Input, int N, the order of the rule.
    Output, double X[N], the nodes.
    Output, double W[N], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
    ityp x01[1] =
    {
        0.0000000000000000
    };
    ityp w01[1] =
    {
        1.0000000000000000
    };
    ityp x03[3] =
    {
        -1.73205080756887719, 0.000000000000000000, 1.73205080756887719
    };
    ityp w03[3] =
    {
        0.166666666666666657, 0.66666666666666663, 0.166666666666666657
    };
    ityp x07[7] =
    {
        -4.18495601767273229, -1.73205080756887719, -0.741095349994540853,
        0.00000000000000000,  0.741095349994540853, 1.73205080756887719,
        4.18495601767273229
    };
    ityp w07[7] =
    {
        0.000695684158369139867, 0.138553274729749237, 0.13137860698313561,
        0.458744868257491889,    0.13137860698313561,  0.138553274729749237,
        0.000695684158369139867
    };
    ityp x09[9] =
    {
        -4.18495601767273229, -2.86127957605705818, -1.73205080756887719,
        -0.741095349994540853, 0.00000000000000000,  0.741095349994540853,
        1.73205080756887719,  2.86127957605705818,  4.18495601767273229
    };
    ityp w09[9] =
    {
        9.42694575565174701E-05,   0.00799632547089352934, 0.0948509485094851251,
        0.270074329577937755,  0.253968253968254065,   0.270074329577937755,
        0.0948509485094851251, 0.00799632547089352934, 9.42694575565174701E-05
    };
    ityp x17[17] =
    {
        -6.36339449433636961,   -5.18701603991365623,       -4.18495601767273229,
        -2.86127957605705818,   -2.59608311504920231,       -1.73205080756887719,
        -1.23042363402730603,   -0.741095349994540853,       0.00000000000000000,
        0.741095349994540853,   1.23042363402730603,        1.73205080756887719,
        2.59608311504920231,    2.86127957605705818,        4.18495601767273229,
        5.18701603991365623,    6.36339449433636961
    };
    ityp w17[17] =
    {
        2.11364995054242569E-08,    -8.20492075415092169E-07,     0.000105637836154169414,
        0.00703348023782790748,  0.0019656770938777492,   0.0886810021520280101,
        0.0141926548264493645,   0.254561232041712215,    0.266922230335053023,
        0.254561232041712215,    0.0141926548264493645,   0.0886810021520280101,
        0.0019656770938777492,   0.00703348023782790748,  0.000105637836154169414,
        -8.20492075415092169E-07,     2.11364995054242569E-08
    };
    ityp x19[19] =
    {
        -6.36339449433636961,    -5.18701603991365623,    -4.18495601767273229,
        -3.20533379449919442,    -2.86127957605705818,    -2.59608311504920231,
        -1.73205080756887719,    -1.23042363402730603,    -0.741095349994540853,
        0.0000000000000000,      0.741095349994540853,    1.23042363402730603,
        1.73205080756887719,     2.59608311504920231,     2.86127957605705818,
        3.20533379449919442,     4.18495601767273229,     5.18701603991365623,
        6.36339449433636961
    };
    ityp w19[19] =
    {
        8.62968460222986318E-10,     6.09480873146898402E-07,     6.01233694598479965E-05,
        0.00288488043650675591, -0.00633722479337375712,  0.0180852342547984622,
        0.0640960546868076103,   0.0611517301252477163,   0.208324991649608771,
        0.303467199854206227,    0.208324991649608771,    0.0611517301252477163,
        0.0640960546868076103,   0.0180852342547984622,  -0.00633722479337375712,
        0.00288488043650675591,  6.01233694598479965E-05,     6.09480873146898402E-07,
        8.62968460222986318E-10
    };
    ityp x31[31] =
    {
        -9.0169397898903032,     -7.98077179859056063,    -7.12210670080461661,
        -6.36339449433636961,    -5.18701603991365623,    -4.18495601767273229,
        -3.63531851903727832,    -3.20533379449919442,    -2.86127957605705818,
        -2.59608311504920231,    -2.23362606167694189,    -1.73205080756887719,
        -1.23042363402730603,    -0.741095349994540853,   -0.248992297579960609,
        0.00000000000000000,     0.248992297579960609,    0.741095349994540853,
        1.23042363402730603,     1.73205080756887719,     2.23362606167694189,
        2.59608311504920231,     2.86127957605705818,     3.20533379449919442,
        3.63531851903727832,     4.18495601767273229,     5.18701603991365623,
        6.36339449433636961,     7.12210670080461661,     7.98077179859056063,
        9.0169397898903032
    };
    ityp w31[31] =
    {
        1.26184642808151181E-15,    -1.4840835740298868E-13,      5.11580531055042083E-12,
        7.92982678648693382E-10,     6.14358432326179133E-07,     5.94749611639316215E-05,
        1.50442053909142189E-05,     0.00272984304673340016, -0.00556100630683581572,
        0.0165924926989360101,   0.00176084755813180017,  0.0617185325658671791,
        0.0654173928360925611,   0.199688635117345498,    0.0281281015400331666,
        0.25890005324151566,     0.0281281015400331666,   0.199688635117345498,
        0.0654173928360925611,   0.0617185325658671791,   0.00176084755813180017,
        0.0165924926989360101,  -0.00556100630683581572,  0.00272984304673340016,
        1.50442053909142189E-05,     5.94749611639316215E-05,     6.14358432326179133E-07,
        7.92982678648693382E-10,     5.11580531055042083E-12,    -1.4840835740298868E-13,
        1.26184642808151181E-15
    };
    ityp x33[33] =
    {
        -9.0169397898903032,     -7.98077179859056063,    -7.12210670080461661,
        -6.36339449433636961,    -5.69817776848810986,    -5.18701603991365623,
        -4.18495601767273229,    -3.63531851903727832,    -3.20533379449919442,
        -2.86127957605705818,    -2.59608311504920231,    -2.23362606167694189,
        -1.73205080756887719,    -1.23042363402730603,    -0.741095349994540853,
        -0.248992297579960609,    0.00000000000000000,     0.248992297579960609,
        0.741095349994540853,    1.23042363402730603,     1.73205080756887719,
        2.23362606167694189,     2.59608311504920231,     2.86127957605705818,
        3.20533379449919442,     3.63531851903727832,     4.18495601767273229,
        5.18701603991365623,     5.69817776848810986,     6.36339449433636961,
        7.12210670080461661,     7.98077179859056063,     9.0169397898903032
    };
    ityp w33[33] =
    {
        -9.93139132868224651E-16,     2.66406251662316506E-13,    -1.93413050008809555E-11,
        1.5542195992782658E-09,     -1.34860173485429301E-08,     6.90862611791137378E-07,
        5.56911589810814793E-05,     8.32360452957667447E-05,     0.00212022595595963252,
        -0.00277121890077892431,  0.01152924706539879,     0.00735301102049550764,
        0.0546775561434630422,   0.0774436027462994808,   0.176075987415714591,
        0.103876871255742839,    0.139110222363380387,    0.103876871255742839,
        0.176075987415714591,    0.0774436027462994808,   0.0546775561434630422,
        0.00735301102049550764,  0.01152924706539879,    -0.00277121890077892431,
        0.00212022595595963252,  8.32360452957667447E-05,     5.56911589810814793E-05,
        6.90862611791137378E-07,    -1.34860173485429301E-08,     1.5542195992782658E-09,
        -1.93413050008809555E-11,     2.66406251662316506E-13,    -9.93139132868224651E-16
    };
    ityp x35[35] =
    {
        -9.0169397898903032,     -7.98077179859056063,    -7.12210670080461661,
        -6.36339449433636961,    -5.69817776848810986,    -5.18701603991365623,
        -4.73643308595229673,    -4.18495601767273229,    -3.63531851903727832,
        -3.20533379449919442,    -2.86127957605705818,    -2.59608311504920231,
        -2.23362606167694189,    -1.73205080756887719,    -1.23042363402730603,
        -0.741095349994540853,   -0.248992297579960609,    0.00000000000000000,
        0.248992297579960609,    0.741095349994540853,    1.23042363402730603,
        1.73205080756887719,     2.23362606167694189,     2.59608311504920231,
        2.86127957605705818,     3.20533379449919442,     3.63531851903727832,
        4.18495601767273229,     4.73643308595229673,     5.18701603991365623,
        5.69817776848810986,     6.36339449433636961,     7.12210670080461661,
        7.98077179859056063,     9.0169397898903032
    };
    ityp w35[35] =
    {
        1.05413265823340136E-18,     5.45004126506381281E-15,     3.09722235760629949E-12,
        4.60117603486559168E-10,     2.13941944795610622E-08,     2.46764213457981401E-07,
        2.73422068011878881E-06,     3.57293481989753322E-05,     0.000275242141167851312,
        0.000818953927502267349, 0.00231134524035220713,  0.00315544626918755127,
        0.015673473751851151,    0.0452736854651503914,   0.0923647267169863534,
        0.148070831155215854,    0.191760115888044341,    0.000514894508069213769,
        0.191760115888044341,    0.148070831155215854,    0.0923647267169863534,
        0.0452736854651503914,   0.015673473751851151,    0.00315544626918755127,
        0.00231134524035220713,  0.000818953927502267349, 0.000275242141167851312,
        3.57293481989753322E-05,     2.73422068011878881E-06,     2.46764213457981401E-07,
        2.13941944795610622E-08,     4.60117603486559168E-10,     3.09722235760629949E-12,
        5.45004126506381281E-15,     1.05413265823340136E-18
    };

    if ( n == 1 )
    {
        r8vec_copy ( n, x01, x );
        r8vec_copy ( n, w01, w );
    }
    else if ( n == 3 )
    {
        r8vec_copy ( n, x03, x );
        r8vec_copy ( n, w03, w );
    }
    else if ( n == 7 )
    {
        r8vec_copy ( n, x07, x );
        r8vec_copy ( n, w07, w );
    }
    else if ( n == 9 )
    {
        r8vec_copy ( n, x09, x );
        r8vec_copy ( n, w09, w );
    }
    else if ( n == 17 )
    {
        r8vec_copy ( n, x17, x );
        r8vec_copy ( n, w17, w );
    }
    else if ( n == 19 )
    {
        r8vec_copy ( n, x19, x );
        r8vec_copy ( n, w19, w );
    }
    else if ( n == 31 )
    {
        r8vec_copy ( n, x31, x );
        r8vec_copy ( n, w31, w );
    }
    else if ( n == 33 )
    {
        r8vec_copy ( n, x33, x );
        r8vec_copy ( n, w33, w );
    }
    else if ( n == 35 )
    {
        r8vec_copy ( n, x35, x );
        r8vec_copy ( n, w35, w );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _kpn_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    KPN_ORDER computes the order of a KPN rule from the level.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int L, the level of the rule.
    1 <= L <= 25
    Output, int KPN_ORDER, the order of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ l = *(dim_typ *) data;
	
    dim_typ n;

    if ( l < 1 )
    {
    	result = USHRT_MAX;
        return &result;
    }
    else if ( l == 1 )
        n = 1;
    else if ( l <= 3 )
        n = 3;
    else if ( l == 4 )
        n = 7;
    else if ( l <= 8 )
        n = 9;
    else if ( l == 9 )
        n = 17;
    else if ( l <= 15 )
        n = 19;
    else if ( l == 16 )
        n = 31;
    else if ( l == 17 )
        n = 33;
    else if ( l <= 25 )
        n = 35;

	result = n;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _kpu ( void * data)
/******************************************************************************/
/*
  Purpose:
    KPU provides data for Kronrod-Patterson quadrature with a uniform weight.
  Discussion:
    This data assumes integration over the interval [0,1] with
    weight function w(x) = 1.
    This data was originally supplied with only 7 digit accuracy.
    It has been replaced by higher accuracy data, which is defined over [-1,+1],
    but adjusted to the interval [0,1] before return.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 December 2012
  Author:
    John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
    Alan Genz, Bradley Keister,
    Fully symmetric interpolatory rules for multiple integrals
    over infinite regions with Gaussian weight,
    Journal of Computational and Applied Mathematics,
    Volume 71, 1996, pages 299-309.
    Thomas Patterson,
    The optimal addition of points to quadrature formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.
  Parameters:
    Input, int N, the order of the rule.
    Only 1, 3, 7, 15, 31 and 63 are legal input values for N.
    Output, double X[N], the nodes.
    Output, double W[N], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
    dim_typ i;

    ityp x01[1] =
    {
        0.0000000
    };
    ityp w01[1] =
    {
        2.0000000
    };
    ityp x03[3] =
    {
        -0.77459666924148337704,
        0.00,
        0.77459666924148337704
    };
    ityp w03[3] =
    {
        0.555555555555555555556,
        0.888888888888888888889,
        0.555555555555555555556
    };
    ityp x07[7] =
    {
        -0.96049126870802028342,
        -0.77459666924148337704,
        -0.43424374934680255800,
        0.00,
        0.43424374934680255800,
        0.77459666924148337704,
        0.96049126870802028342
    };
    ityp w07[7] =
    {
        0.104656226026467265194,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194
    };
    ityp x15[15] =
    {
        -0.99383196321275502221,
        -0.96049126870802028342,
        -0.88845923287225699889,
        -0.77459666924148337704,
        -0.62110294673722640294,
        -0.43424374934680255800,
        -0.22338668642896688163,
        0.00,
        0.22338668642896688163,
        0.43424374934680255800,
        0.62110294673722640294,
        0.77459666924148337704,
        0.88845923287225699889,
        0.96049126870802028342,
        0.99383196321275502221
    };
    ityp w15[15] =
    {
        0.0170017196299402603390,
        0.0516032829970797396969,
        0.0929271953151245376859,
        0.134415255243784220360,
        0.171511909136391380787,
        0.200628529376989021034,
        0.219156858401587496404,
        0.225510499798206687386,
        0.219156858401587496404,
        0.200628529376989021034,
        0.171511909136391380787,
        0.134415255243784220360,
        0.0929271953151245376859,
        0.0516032829970797396969,
        0.0170017196299402603390
    };
    ityp x31[31] =
    {
        -0.99909812496766759766,
        -0.99383196321275502221,
        -0.98153114955374010687,
        -0.96049126870802028342,
        -0.92965485742974005667,
        -0.88845923287225699889,
        -0.83672593816886873550,
        -0.77459666924148337704,
        -0.70249620649152707861,
        -0.62110294673722640294,
        -0.53131974364437562397,
        -0.43424374934680255800,
        -0.33113539325797683309,
        -0.22338668642896688163,
        -0.11248894313318662575,
        0.00,
        0.11248894313318662575,
        0.22338668642896688163,
        0.33113539325797683309,
        0.43424374934680255800,
        0.53131974364437562397,
        0.62110294673722640294,
        0.70249620649152707861,
        0.77459666924148337704,
        0.83672593816886873550,
        0.88845923287225699889,
        0.92965485742974005667,
        0.96049126870802028342,
        0.98153114955374010687,
        0.99383196321275502221,
        0.99909812496766759766
    };
    ityp w31[31] =
    {
        0.00254478079156187441540,
        0.00843456573932110624631,
        0.0164460498543878109338,
        0.0258075980961766535646,
        0.0359571033071293220968,
        0.0464628932617579865414,
        0.0569795094941233574122,
        0.0672077542959907035404,
        0.0768796204990035310427,
        0.0857559200499903511542,
        0.0936271099812644736167,
        0.100314278611795578771,
        0.105669893580234809744,
        0.109578421055924638237,
        0.111956873020953456880,
        0.112755256720768691607,
        0.111956873020953456880,
        0.109578421055924638237,
        0.105669893580234809744,
        0.100314278611795578771,
        0.0936271099812644736167,
        0.0857559200499903511542,
        0.0768796204990035310427,
        0.0672077542959907035404,
        0.0569795094941233574122,
        0.0464628932617579865414,
        0.0359571033071293220968,
        0.0258075980961766535646,
        0.0164460498543878109338,
        0.00843456573932110624631,
        0.00254478079156187441540
    };
    ityp x63[63] =
    {
        -0.99987288812035761194,
        -0.99909812496766759766,
        -0.99720625937222195908,
        -0.99383196321275502221,
        -0.98868475754742947994,
        -0.98153114955374010687,
        -0.97218287474858179658,
        -0.96049126870802028342,
        -0.94634285837340290515,
        -0.92965485742974005667,
        -0.91037115695700429250,
        -0.88845923287225699889,
        -0.86390793819369047715,
        -0.83672593816886873550,
        -0.80694053195021761186,
        -0.77459666924148337704,
        -0.73975604435269475868,
        -0.70249620649152707861,
        -0.66290966002478059546,
        -0.62110294673722640294,
        -0.57719571005204581484,
        -0.53131974364437562397,
        -0.48361802694584102756,
        -0.43424374934680255800,
        -0.38335932419873034692,
        -0.33113539325797683309,
        -0.27774982202182431507,
        -0.22338668642896688163,
        -0.16823525155220746498,
        -0.11248894313318662575,
        -0.056344313046592789972,
        0.00,
        0.056344313046592789972,
        0.11248894313318662575,
        0.16823525155220746498,
        0.22338668642896688163,
        0.27774982202182431507,
        0.33113539325797683309,
        0.38335932419873034692,
        0.43424374934680255800,
        0.48361802694584102756,
        0.53131974364437562397,
        0.57719571005204581484,
        0.62110294673722640294,
        0.66290966002478059546,
        0.70249620649152707861,
        0.73975604435269475868,
        0.77459666924148337704,
        0.80694053195021761186,
        0.83672593816886873550,
        0.86390793819369047715,
        0.88845923287225699889,
        0.91037115695700429250,
        0.92965485742974005667,
        0.94634285837340290515,
        0.96049126870802028342,
        0.97218287474858179658,
        0.98153114955374010687,
        0.98868475754742947994,
        0.99383196321275502221,
        0.99720625937222195908,
        0.99909812496766759766,
        0.99987288812035761194
    };
    ityp w63[63] =
    {
        0.000363221481845530659694,
        0.00126515655623006801137,
        0.00257904979468568827243,
        0.00421763044155885483908,
        0.00611550682211724633968,
        0.00822300795723592966926,
        0.0104982469096213218983,
        0.0129038001003512656260,
        0.0154067504665594978021,
        0.0179785515681282703329,
        0.0205942339159127111492,
        0.0232314466399102694433,
        0.0258696793272147469108,
        0.0284897547458335486125,
        0.0310735511116879648799,
        0.0336038771482077305417,
        0.0360644327807825726401,
        0.0384398102494555320386,
        0.0407155101169443189339,
        0.0428779600250077344929,
        0.0449145316536321974143,
        0.0468135549906280124026,
        0.0485643304066731987159,
        0.0501571393058995374137,
        0.0515832539520484587768,
        0.0528349467901165198621,
        0.0539054993352660639269,
        0.0547892105279628650322,
        0.0554814043565593639878,
        0.0559784365104763194076,
        0.0562776998312543012726,
        0.0563776283603847173877,
        0.0562776998312543012726,
        0.0559784365104763194076,
        0.0554814043565593639878,
        0.0547892105279628650322,
        0.0539054993352660639269,
        0.0528349467901165198621,
        0.0515832539520484587768,
        0.0501571393058995374137,
        0.0485643304066731987159,
        0.0468135549906280124026,
        0.0449145316536321974143,
        0.0428779600250077344929,
        0.0407155101169443189339,
        0.0384398102494555320386,
        0.0360644327807825726401,
        0.0336038771482077305417,
        0.0310735511116879648799,
        0.0284897547458335486125,
        0.0258696793272147469108,
        0.0232314466399102694433,
        0.0205942339159127111492,
        0.0179785515681282703329,
        0.0154067504665594978021,
        0.0129038001003512656260,
        0.0104982469096213218983,
        0.00822300795723592966926,
        0.00611550682211724633968,
        0.00421763044155885483908,
        0.00257904979468568827243,
        0.00126515655623006801137,
        0.000363221481845530659694
    };

    if ( n == 1 )
    {
        r8vec_copy ( n, x01, x );
        r8vec_copy ( n, w01, w );
    }
    else if ( n == 3 )
    {
        r8vec_copy ( n, x03, x );
        r8vec_copy ( n, w03, w );
    }
    else if ( n == 7 )
    {
        r8vec_copy ( n, x07, x );
        r8vec_copy ( n, w07, w );
    }
    else if ( n == 15 )
    {
        r8vec_copy ( n, x15, x );
        r8vec_copy ( n, w15, w );
    }
    else if ( n == 31 )
    {
        r8vec_copy ( n, x31, x );
        r8vec_copy ( n, w31, w );
    }
    else if ( n == 63 )
    {
        r8vec_copy ( n, x63, x );
        r8vec_copy ( n, w63, w );
    }

    /*
    The rule as stored is for the interval [-1,+1].
    Adjust it to the interval [0,1].
    */
    rule_adjust ( -1.00, +1.00, 0.00, +1.00, n, x, w );

  return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _kpu_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    KPU_ORDER computes the order of a KPU rule from the level.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int L, the level of the rule.
    1 <= L <= 25
    Output, int KPU_ORDER, the order of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ l = *(dim_typ *) data;
	
    dim_typ n;

    if ( l < 1 )
    {
    	result = USHRT_MAX;
        return &result;
    }
    else if ( l == 1 )
        n = 1;
    else if ( l <= 3 )
        n = 3;
    else if ( l <= 6 )
        n = 7;
    else if ( l <= 12 )
        n = 15;
    else if ( l <= 24 )
        n = 31;
    else if ( l <= 25 )

        n = 63;
	
	result = n;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _num_seq ( void * data)
/******************************************************************************/
/*
  Purpose:
    NUM_SEQ returns the number of compositions of the integer N into K parts.
  Discussion:
    A composition of the integer N into K parts is an ordered sequence
    of K nonnegative integers which sum to N.  The compositions (1,2,1)
    and (1,1,2) are considered to be distinct.
    The 28 compositions of 6 into three parts are:
      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
      0 3 3,  0 2 4,  0 1 5,  0 0 6.
    The formula for the number of compositions of N into K parts is
      Number = ( N + K - 1 )! / ( N! * ( K - 1 )! )
    Describe the composition using N '1's and K-1 dividing lines '|'.
    The number of distinct permutations of these symbols is the number
    of compositions.  This is equal to the number of permutations of
    N+K-1 things, with N identical of one kind and K-1 identical of another.
    Thus, for the above example, we have:
      Number = ( 6 + 3 - 1 )! / ( 6! * (3-1)! ) = 28
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the integer whose compositions are desired.
    Input, int K, the number of parts in the composition.
    Output, int NUM_SEQ, the number of compositions of N
    into K parts.
*/
{
	static dim_typ result = USHRT_MAX;
	
	int * const a_data = data;
	const register int n = a_data[0];
	const register int k = a_data[1];
	
	result = i4_choose ( n + k - 1, n );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   nwspgr ( void rule ( dim_typ n, ityp x[], ityp w[] ),
  dim_typ rule_order ( dim_typ l ), const register dim_typ dim, const register dim_typ k, const register dim_typ r_size, dim_typ *s_size,ityp nodes[static dim*r_size], ityp weights[static r_size] )
/******************************************************************************/
/*
  Purpose:
    NWSPGR generates nodes and weights for sparse grid integration.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 April 2013
  Author:
    Original MATLAB version by Florian Heiss, Viktor Winschel.
    C version by John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
  Parameters:
    Input, void RULE ( int n, double x[], double w[] ), the name of a function
    which is given the order N and returns the points X and weights W of the
    corresponding 1D quadrature rule.
    Input, int RULE_ORDER ( int l ), the name of a function which
    is given the level L and returns the order N of the corresponding 1D rule.
    Input, int DIM, the spatial dimension.
    Input, int K, the level of the sparse rule.
    Input, int R_SIZE, the "size" of the sparse rule.
    Output, int *S_SIZE, the size of the sparse rule, after
    duplicate points have been merged.
    Output, double NODES[DIM*R_SIZE], the nodes of the sparse rule.
    Output, double WEIGHTS[R_SIZE], the weights of the sparse rule.
*/
{
    int bq;
    bool equal;
    dim_typ i;
    int *is;
    dim_typ j;
    dim_typ j2;
    int *keep;
    int level;
    int *lr;
    dim_typ maxq;
    dim_typ minq;
    dim_typ n;
    int *n1d;
    dim_typ n1d_total;
    dim_typ nc;
    dim_typ np;
    int *nr;
    dim_typ q;
    dim_typ r;
    int *roff;
    int *rq;
    dim_typ seq_num;
    ityp t;
    ityp *w;
    ityp *w1d;
    int *w1d_off;
    ityp *wc;
    ityp *wp;
    ityp *wr;
    ityp *x;
    ityp *x1d;
    int *x1d_off;
    ityp *xc;
    ityp *xp;
    ityp *xr;

    for ( j = 0; j < r_size; ++j)
    {
        for ( i = 0; i < dim; ++i )
            nodes[i+j*dim] = 0.00;
        weights[j] = 0.00;
    }
    /*
    Create cell arrays that will contain the points and weights
    for levels 1 through K.
    */
    n1d = ( int * ) malloc ( k * sizeof ( int ) );
    x1d_off = ( int * ) malloc ( ( k + 1 ) * sizeof ( int ) );
    w1d_off = ( int * ) malloc ( ( k + 1 ) * sizeof ( int ) );

    x1d_off[0] = w1d_off[0] = 0;

    for ( level = 1; level <= k; ++level )
    {
        n = rule_order ( level );
        n1d[level-1] = n;
        x1d_off[level] = x1d_off[level-1] + n;
        w1d_off[level] = w1d_off[level-1] + n;
    }

    n1d_total = x1d_off[k];
    /*
    Calculate all the 1D rules needed.
    */
    x1d = ( ityp * ) malloc ( n1d_total * sizeof ( ityp ) );
    w1d = ( ityp * ) malloc ( n1d_total * sizeof ( ityp ) );

    for ( level = 1; level <= k; ++level)
    {
        n = n1d[level-1];

        x = ( ityp * ) malloc ( n * sizeof ( ityp ) );
        w = ( ityp * ) malloc ( n * sizeof ( ityp ) );

        rule ( n, x, w );
        r8cvv_rset ( n1d_total, x1d, k, x1d_off, level - 1, x );
        r8cvv_rset ( n1d_total, w1d, k, w1d_off, level - 1, w );

        free ( x );
        free ( w );
    }
    /*
    Construct the sparse grid.
    */
    minq = MAX ( 0, k - dim );
    maxq = k - 1;
    /*
    Q is the level total.
    */
    lr = ( int * ) malloc ( dim * sizeof ( int ) );
    nr = ( int * ) malloc ( dim * sizeof ( int ) );

    r = 0;

    for ( q = minq; q <= maxq; ++q )
    {
        /*
        BQ is the combinatorial coefficient applied to the component
        product rules which have level Q.
        */
        bq = i4_mop ( maxq - q ) * i4_choose ( dim - 1, dim + q - k );
        /*
        Compute the D-dimensional row vectors that sum to DIM+Q.
        */
        seq_num = num_seq ( q, dim );

        is = get_seq ( dim, q + dim, seq_num );
        /*
        Allocate new rows for nodes and weights.
        */
        rq = ( int * ) malloc ( seq_num * sizeof ( int ) );

        for ( j = 0; j < seq_num; ++j)
        {
            rq[j] = 1;
            for ( i = 0; i < dim; ++i)
            {
                level = is[j+i*seq_num] - 1;
                rq[j] = rq[j] * n1d[level];
            }
        }
        /*
        Generate every product rule whose level total is Q.
        */
        for ( j2 = 0; j2 < seq_num; ++j2 )
        {
            for ( i = 0; i < dim; ++i )
            {
                lr[i] = is[j2+i*seq_num];
                nr[i] = rule_order ( lr[i] );
            }


            roff = r8cvv_offset ( dim, nr );

            nc = i4vec_sum ( dim, nr );
            wc = ( ityp * ) malloc ( nc * sizeof ( ityp ) );
            xc = ( ityp * ) malloc ( nc * sizeof ( ityp ) );
            /*
            Corrected first dimension to N1D_TOTAL in following two calls
            to R8CVV_GET_NEW, 19 April 2013.
            */
            for ( i = 0; i < dim; ++i )
            {
                xr = r8cvv_rget_new ( n1d_total, x1d, k, x1d_off, lr[i] - 1 );
                wr = r8cvv_rget_new ( n1d_total, w1d, k, w1d_off, lr[i] - 1 );
                r8cvv_rset ( nc, xc, dim, roff, i, xr );
                r8cvv_rset ( nc, wc, dim, roff, i, wr );
                free ( wr );
                free ( xr );
            }

            np = rq[j2];
            wp = ( ityp * ) malloc ( np * sizeof ( ityp ) );
            xp = ( ityp * ) malloc ( dim * np * sizeof ( ityp ) );

            tensor_product_cell ( nc, xc, wc, dim, nr, roff, np, xp, wp );
            /*
            Append the new nodes and weights to the arrays.
            */
            for ( j = 0; j < np; ++j )
            {
                for ( i = 0; i < dim; ++i )
                    nodes[i+(r+j)*dim] = xp[i+j*dim];
                weights[r+j] = bq*wp[j];
            }
            /*
            Update the count.
            */
            r = r + rq[j2];

            free ( roff );
            free ( wc );
            free ( wp );
            free ( xc );
            free ( xp );
        }
    free ( is );
    free ( rq );
    }

    free ( lr );
    free ( n1d );
    free ( nr );
    free ( w1d );
    free ( w1d_off );
    free ( x1d );
    free ( x1d_off );
    /*
    Reorder the rule so the points are in ascending lexicographic order.
    */
    rule_sort ( dim, r_size, nodes, weights );
    /*
    Suppress duplicate points and merge weights.
    */
    r = 0;
    for ( j = 1; j < r_size; ++j)
    {
        equal = true;
        for ( i = 0; i < dim; ++i)
        {
            if ( nodes[i+r*dim] != nodes[i+j*dim] )
            {
                equal = false;
                break;
            }
        }
        if ( equal )
            weights[r] += weights[j];
        else
        {
            ++ r;
            weights[r] = weights[j];
            for ( i = 0; i < dim; ++i )
                nodes[i+r*dim] = nodes[i+j*dim];
        }
    }

    ++ r;
    *s_size = r;
    /*
    Zero out unneeded entries.
    */
    for ( j = r; j < r_size; ++j )
    {
        for ( i = 0; i < dim; ++i)
            nodes[i+j*dim] = 0.00;
        weights[j] = 0.00;
    }

    /*
    Normalize the weights to sum to 1.
    */
    t = r8vec_sum ( r, weights );

    for ( j = 0; j < r; ++j )
        weights[j] /= t;

    return;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  dim_typ   nwspgr_size ( dim_typ rule_order ( dim_typ l ), const register dim_typ dim, const register dim_typ k )
/******************************************************************************/
/*
  Purpose:
    NWSPGR_SIZE determines the size of a sparse grid rule.
  Discussion:
    This routine does a "raw" count, that is, it does not notice that many
    points may be repeated, in which case, the size of the rule could be
    reduced by merging repeated points and combining the corresponding weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 December 2012
  Author:
    John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
  Parameters:
    Input, int RULE_ORDER ( int l ), the name of a function which
    is given the level L and returns the order N of the corresponding rule.
    Input, int DIM, the dimension of the integration problem.
    Input, int K, the level.  When using the built in 1D
    rules, the resulting sparse grid will be exact for polynomials up to total
    order 2*K-1.  When using the built-in quadrature rules, the maximum value
    of K that is available is 25.
    Output, int NWSPGR_SIZE, the "size" of the rule, that is,
    the number of weights and multidimensional quadrature points that will
    be needed.  The size of the rule will be reduced when duplicate points
    are merged.
*/
{
    dim_typ i;
    int *is;
    dim_typ j;
    dim_typ level;
    dim_typ maxq;
    dim_typ minq;
    dim_typ n;
    int *n1d;
    dim_typ n1d_total;
    dim_typ q;
    dim_typ r_size;
    int *rq;
    dim_typ seq_num;
    /*
    Determine the size of each 1D rule.
    */
    n1d = ( int * ) malloc ( k * sizeof ( int ) );

    for ( level = 1; level <= k; ++level )
    {
        n = rule_order ( level );
        n1d[level-1] = n;
    }

    n1d_total = i4vec_sum ( k, n1d );
    /*
    Go through the motions of generating the rules.
    */
    minq = MAX ( 0, k - dim );
    maxq = k - 1;
    r_size = 0;

    for ( q = minq; q <= maxq; ++q)
    {
        /*
        Compute the D-dimensional vectors that sum to Q+DIM.
        */
        seq_num = num_seq ( q, dim );

        is = get_seq ( dim, q + dim, seq_num );
        /*
        Determine the size of each rule.
        */
        rq = ( int * ) malloc ( seq_num * sizeof ( int ) );

        for ( j = 0; j < seq_num; ++j )
        {
            rq[j] = 1;
            for ( i = 0; i < dim; ++i )
                rq[j] *= n1d[is[j+i*seq_num]-1];
        }
        /*
        Add the sizes to the total.
        */
        r_size = r_size + i4vec_sum ( seq_num, rq );

        free ( is );
        free ( rq );
    }

    free ( n1d );

    return r_size;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rule_adjust ( void * data)
/******************************************************************************/
/*
  Purpose:
    RULE_ADJUST adjusts a 1D quadrature rule from [A,B] to [C,D].
  Discussion:
    This function is only appropriate for cases involving finite intervals
    and a uniform weighting function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 February 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, the left and right endpoints of the
    original interval.
    Input, double C, D, the left and right endpoints of the
    new interval.
    Input, int N, the order of the rule.
    1 <= N.
    Input/output, double X[N], the abscissas.
    Input/output, double W[N], the weights.
*/
{
	const _4itdt2pit * const s_data = data;
	ityp a = s_data->a0;
	ityp b = s_data->a1;
	ityp c = s_data->a2;
	ityp d = s_data->a3;
	const register dim_typ n = s_data->a4;
	ityp * x = s_data->a5;
	ityp * w = s_data->a6;
	
    dim_typ i;
    ityp s;

    for ( i = 0; i < n; ++i )
        x[i] = ( ( b - x[i]     ) * c+ (     x[i] - a ) * d )/ ( b        - a );

    s = ( d - c ) / ( b - a );

    for ( i = 0; i < n; ++i )
        w[i] *= s;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _rule_sort ( void * data)
/******************************************************************************/
/*
  Purpose:
    RULE_SORT sorts a multidimensional quadrature rule.
  Discussion:
    This routine reorders the items in the rule so that the points
    occur in increasing lexicographic order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, double X[M*N], W[N].
    The points and weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
    int i, j1, j2;
    int indx;
    short isgn;
    ityp t;
    ityp ww;

    if ( m == 0 || n == 1 )
        return NULL;

    /*
    Initialize.
    */
    indx = isgn = j1 = j2 = 0;
    /*
    Call the external heap sorter.
    */
    while ( true )
    {
        sort_heap_external ( n, &indx, &j1, &j2, isgn );
        /*
        Interchange columns J1 and J2.
        */
        if ( 0 < indx )
        {
            for ( i = 0; i < m; ++i )
            {
                t             = x[i+(j1-1)*m];
                x[i+(j1-1)*m] = x[i+(j2-1)*m];
                x[i+(j2-1)*m] = t;
            }

            ww      = w[j1-1];
            w[j1-1] = w[j2-1];
            w[j2-1] = ww;
        }
        /*
        Compare columns J1 and J2.
        */
        else if ( indx < 0 )
        {
            isgn = 0;
            for ( i = 0; i < m; ++i )
            {
                if ( x[i+(j1-1)*m] < x[i+(j2-1)*m] )
                {
                    isgn = -1;
                    break;
                }
                else if ( x[i+(j2-1)*m] < x[i+(j1-1)*m] )
                {
                    isgn = +1;
                    break;
                }
            }
        }
        /*
        The columns are sorted.
        */
        else if ( indx == 0 )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _symmetric_sparse_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    SYMMETRIC_SPARSE_SIZE sizes a symmetric sparse rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2012
  Author:
    John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
  Parameters:
    Input, int DIM, the dimension.
    Input, int NR, the dimension of the rule in the
    positive orthant.
    Input, double NODES[NR*DIM], the nodes for the positive orthant.
    Input, double X0, the point of symmetry for the 1D rule,
    typically 0.
    Output, int SYMMETRIC_SPARSE_SIZE, the dimension of the rule
    when "unfolded" to the full space.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpitit * const s_data = data;
	const register dim_typ nr = s_data->a0;
	const register dim_typ dim = s_data->a1;
	ityp * nodes = s_data->a2;
	const register ityp x0 = s_data->a3;
	
    dim_typ count;
    dim_typ j;
    dim_typ nr2 = 0;
    dim_typ r;
    /*
    Count the size of the full rule.
    */

    for ( r = 0; r < nr; ++r)
    {
        count = 1;
        for ( j = 0; j < dim; ++j )
            if ( nodes[r+j*nr] != x0 )
                count = count<<1;
        nr2 += count;
    }
	
	result = nr2;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _tensor_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    TENSOR_PRODUCT generates a tensor product quadrature rule.
  Discussion:
    The Kronecker product of an K by L matrix A and an M by N matrix B
    is the K*M by L*N matrix formed by
      a(1,1) * B,  a(1,2) * B,  ..., a(1,l) * B
      a(2,1) * B,  a(2,2) * B,  ..., a(2,l) * B
      ..........   ..........   .... ..........
      a(k,1) * B,  a(k,2) * B,  ..., a(k,l) * B
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 December 2012
  Author:
    Original MATLAB version by Florian Heiss, Viktor Winschel.
    C version by John Burkardt.
  Reference:
    Florian Heiss, Viktor Winschel,
    Likelihood approximation by numerical integration on sparse grids,
    Journal of Econometrics,
    Volume 144, 2008, pages 62-80.
  Parameters:
    Input, int D, the spatial dimension.
    Input, int ORDER1D[D], the order of each 1D rule.
    Input, int N1D, the number of 1D items.
    Input, double X1D[N1D], the 1D nodes.
    Input, double W1D[N1D], the 1D weights.
    Input, int N, the number of N-dimensional items.
    Output, double XND[D*N], the nodes.
    Output, double WND[N], the weights.
*/
{
	const dtpit2dtpi3pit * const s_data = data;
	
	const register dim_typ d = s_data->a0;
	ityp * x1d = s_data->a1;
	const register dim_typ n1d = s_data->a2;
	const register dim_typ n = s_data->a3;	
	int * order1d = s_data->a4;
	ityp * w1d = s_data->a5;
	ityp * xnd = s_data->a6;
	ityp * wnd = s_data->a7;
	
    dim_typ i;
    dim_typ i1;
    short i2;
    /*
    Compute the weights.
    */
    i2 = -1;
    for ( i = 0; i < d; ++i )
    {
        i1 = i2 + 1;
        i2 += order1d[i];
        r8vec_direct_product2 ( i, order1d[i], w1d+i1, d, n, wnd );
    }
    /*
    Compute the points.
    */
    i2 = -1;
    for ( i = 0; i < d; ++i )
    {
        i1 = i2 + 1;
        i2 += order1d[i];
        r8vec_direct_product ( i, order1d[i], x1d+i1, d, n, xnd );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tensor_product_cell ( void * data)
/******************************************************************************/
/*
  Purpose:
    TENSOR_PRODUCT_CELL generates a tensor product quadrature rule.
  Discussion:
    The Kronecker product of an K by L matrix A and an M by N matrix B
    is the K*M by L*N matrix formed by
      a(1,1) * B,  a(1,2) * B,  ..., a(1,l) * B
      a(2,1) * B,  a(2,2) * B,  ..., a(2,l) * B
      ..........   ..........   .... ..........
      a(k,1) * B,  a(k,2) * B,  ..., a(k,l) * B
    The 1D factors are stored in a kind of cell array structure.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 December 2012
  Author:
    John Burkardt.
  Parameters:
    Input, int NC, the number of items in the cell arrays.
    Input, double XC[NC], a cell array containing points for
    1D rules.
    Input, double WC[NC], a cell array containing weights for
    1D rules.
    Input, int DIM, the spatial dimension.
    Input, int NR[DIM], the length of each row of the
    cell array.
    Input, int ROFF[DIM+1], offsets for the cell arrays.
    Input, int NP, the number of points in the product rule.
    Output, double XP[DIM*NP], the nodes.
    Output, double WP[NP], the weights.
*/
{
	const dt2pitdt2pidt2pit * const s_data = data;
	const register dim_typ nc = s_data->a0;
	ityp * xc = s_data->a1;
	ityp * wc = s_data->a2;
	const register dim_typ dim = s_data->a3;
	int * nr = s_data->a4;
	int * roff = s_data->a5;
	const register dim_typ np = s_data->a6;
	ityp * xp = s_data->a7;
	ityp * wp = s_data->a8;
	
	
    dim_typ i;
    dim_typ n1d;
    ityp *w1d;
    ityp *x1d;
    /*
    Compute the weights.
    */
    for ( i = 0; i < dim; ++i )
    {
        n1d = nr[i];
        w1d = r8cvv_rget_new ( nc, wc, dim, roff, i );
        r8vec_direct_product2 ( i, n1d, w1d, dim, np, wp );
        free ( w1d );
    }
    /*
    Compute the points.
    */
    for ( i = 0; i < dim; ++i )
    {
        n1d = nr[i];
        x1d = r8cvv_rget_new ( nc, xc, dim, roff, i );
        r8vec_direct_product ( i, n1d, x1d, dim, np, xp );
        free ( x1d );
    }

    return NULL;
}

#endif
