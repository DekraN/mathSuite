#ifndef __DISABLEDEEP_WAVELET

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_convolution ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_CONVOLUTION returns the convolution of two R8VEC's.
  Discussion:
    An R8VEC is a vector of R8's.
    The I-th entry of the convolution can be formed by summing the products 
    that lie along the I-th diagonal of the following table:
    Y3 | 3   4   5   6   7
    Y2 | 2   3   4   5   6
    Y1 | 1   2   3   4   5
       +------------------
        X1  X2  X3  X4  X5
    which will result in:
    Z = ( X1 * Y1,
          X1 * Y2 + X2 * Y1,
          X1 * Y3 + X2 * Y2 + X3 * Y1,
                    X2 * Y3 + X3 * Y2 + X4 * Y1,
                              X3 * Y3 + X4 * Y2 + X5 * Y1,
                                        X4 * Y3 + X5 * Y2,
                                                  X5 * Y3 )     
  Example:
    Input:
      X = (/ 1, 2, 3, 4 /)
      Y = (/ -1, 5, 3 /)
    Output:
      Z = (/ -1, 3, 10, 17, 29, 12 /)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:

    05 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of X.
    Input, double X[M], the first vector to be convolved.
    Input, int N, the dimension of Y.
    Input, double Y[N], the second vector to be convolved.
    Output, double R8VEC_CONVOLUTION[M+N-1], the convolution of X and Y.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * y = s_data->a3;
	
	dim_typ i, j;
	ityp *z = ( ityp * ) malloc ( ( m + n - 1 ) * sizeof ( ityp ) );
	
	for ( i = 0; i < m + n - 1; ++i )
		z[i] = 0.00;
	
	for ( j = 0; j < n; ++j)
		for ( i = 0; i < m; ++i )
			z[j+i] += x[i] * y[j];

	return z;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cascade ( void * data)
/******************************************************************************/
/*
  Purpose:
    CASCADE carries out the cascade algorithm.
  Discussion:
    The value of T3 computed by
      call cascade ( 3, t_length, t0, c_length, c, t3 )
    will be the same if computed in three steps by:
      call cascade ( 1, t_length, t0, c_length, c, t1 );
      call cascade ( 1, t_length, t1, c_length, c, t2 );
      call cascade ( 1, t_length, t2, c_length, c, t3 );
    If C represents a vector of Daubechies filter coefficients, then
      call cascade ( 5, c_length, c, c_length, c, c5 );
    computes an approximation to the corresponding scaling function, and
      call r8vec_conjugate ( c_length, c, w )
      call cascade ( 5, c_length, w, c_length, c, w5 );
    computes an approximation to the corresponding wavelet.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of iterations to carry out.
    0 <= N.
    Input, int T_LENGTH, the length of T.
    Input, double T[T_LENGTH], the initial value of the quantity,
    or the value of the quantity at the integers 0 through P-1.
    Input, int C_LENGTH, the number of transform coefficients.
    Input, double C[C_LENGTH], the transform coefficients.
    Output, double CASCADE[2^N * T_LENGTH + (2^N-1)*C_LENGTH - 2*(2^N-1)],
    the values of  the function.
*/
{
	const _3dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ t_length = s_data->a1;
	const register dim_typ c_length = s_data->a2;
	ityp * t = s_data->a3;
	ityp * c = s_data->a4;
	
    dim_typ i, j;
    ityp *s;
    dim_typ s_length;
    dim_typ x_length;
    ityp *x;

    s_length = t_length;

    s = r8vec_copy_new ( t_length, t );

    for ( i = 0; i < n; ++i )
    {
        x_length = (s_length <<1) - 1;

        x = ( ityp * ) malloc ( x_length * sizeof ( ityp ) );

        j = 0;
        for ( i = 0; i < x_length; i += 2 )
        {
            x[i] = s[j];
            ++ j;
        }

        for ( i = 1; i < x_length - 1; i += 2 )
            x[i] = 0.0;

        free ( s );
        s = r8vec_convolution ( x_length, x, c_length, c );
        free ( x );
        s_length = x_length + c_length - 1;
    }

    return s;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub_coefficients ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB_COEFFICIENTS returns a set of Daubechies coefficients.
  Discussion:
    Often, the uses to which these coefficients are applied require that they
    be rescaled, by being multiplied by sqrt ( 2 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the coefficient set.
    2 <= N <= 20, and N must be even.
    Output, double DAUB_COEFFICIENTS[N], the coefficients.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *c;
    static ityp c02[2] =
    {
        7.071067811865475E-01,
        7.071067811865475E-01
    };
    static ityp c04[4] =
    {
        0.4829629131445341E+00,
        0.8365163037378079E+00,
        0.2241438680420133E+00,
        - 0.1294095225512603E+00
    };
    static ityp c06[6] =
    {
        0.3326705529500826E+00,
        0.8068915093110925E+00,
        0.4598775021184915E+00,
        - 0.1350110200102545E+00,
        - 0.8544127388202666E-01,
        0.3522629188570953E-01
    };
    static ityp c08[8] =
    {
        0.2303778133088965E+00,
        0.7148465705529156E+00,
        0.6308807679298589E+00,
        -0.2798376941685985E-01,
        -0.1870348117190930E+00,
        0.3084138183556076E-01,
        0.3288301166688519E-01,
        -0.1059740178506903E-01
    };
    static ityp c10[10] =
    {
        0.1601023979741929E+00,
        0.6038292697971896E+00,
        0.7243085284377729E+00,
        0.1384281459013207E+00,
        -0.2422948870663820E+00,
        -0.3224486958463837E-01,
        0.7757149384004571E-01,
        -0.6241490212798274E-02,
        -0.1258075199908199E-01,
        0.3335725285473771E-02
    };
    static ityp c12[12] =
    {
        0.1115407433501094E+00,
        0.4946238903984530E+00,
        0.7511339080210953E+00,
        0.3152503517091976E+00,
        -0.2262646939654398E+00,
        -0.1297668675672619E+00,
        0.9750160558732304E-01,
        0.2752286553030572E-01,
        -0.3158203931748602E-01,
        0.5538422011614961E-03,
        0.4777257510945510E-02,
        -0.1077301085308479E-02
    };
    static ityp c14[14] =
    {
        7.785205408500917E-02,
        3.965393194819173E-01,
        7.291320908462351E-01,
        4.697822874051931E-01,
        -1.439060039285649E-01,
        -2.240361849938749E-01,
        7.130921926683026E-02,
        8.061260915108307E-02,
        -3.802993693501441E-02,
        -1.657454163066688E-02,
        1.255099855609984E-02,
        4.295779729213665E-04,
        -1.801640704047490E-03,
        3.537137999745202E-04
    };
    static ityp c16[16] =
    {
        5.441584224310400E-02,
        3.128715909142999E-01,
        6.756307362972898E-01,
        5.853546836542067E-01,
        -1.582910525634930E-02,
        -2.840155429615469E-01,
        4.724845739132827E-04,
        1.287474266204784E-01,
        -1.736930100180754E-02,
        -4.408825393079475E-02,
        1.398102791739828E-02,
        8.746094047405776E-03,
        -4.870352993451574E-03,
        -3.917403733769470E-04,
        6.754494064505693E-04,
        -1.174767841247695E-04
    };
    static ityp c18[18] =
    {
        3.807794736387834E-02,
        2.438346746125903E-01,
        6.048231236901111E-01,
        6.572880780513005E-01,
        1.331973858250075E-01,
        -2.932737832791749E-01,
        -9.684078322297646E-02,
        1.485407493381063E-01,
        3.072568147933337E-02,
        -6.763282906132997E-02,
        2.509471148314519E-04,
        2.236166212367909E-02,
        -4.723204757751397E-03,
        -4.281503682463429E-03,
        1.847646883056226E-03,
        2.303857635231959E-04,
        -2.519631889427101E-04,
        3.934732031627159E-05
    };
    static ityp c20[20] =
    {
        2.667005790055555E-02,
        1.881768000776914E-01,
        5.272011889317255E-01,
        6.884590394536035E-01,
        2.811723436605774E-01,
        -2.498464243273153E-01,
        -1.959462743773770E-01,
        1.273693403357932E-01,
        9.305736460357235E-02,
        -7.139414716639708E-02,
        -2.945753682187581E-02,
        3.321267405934100E-02,
        3.606553566956169E-03,
        -1.073317548333057E-02,
        1.395351747052901E-03,
        1.992405295185056E-03,
        -6.858566949597116E-04,
        -1.164668551292854E-04,
        9.358867032006959E-05,
        -1.326420289452124E-05
    };

    switch(n)
    {
        case 2:
            c = r8vec_copy_new ( n, c02 );
            break;
        case 4:
            c = r8vec_copy_new ( n, c04 );
            break;
        case 6:
            c = r8vec_copy_new ( n, c06 );
            break;
        case 8:
            c = r8vec_copy_new ( n, c08 );
            break;
        case 10:
            c = r8vec_copy_new ( n, c10 );
            break;
        case 12:
            c = r8vec_copy_new ( n, c12 );
            break;
        case 14:
            c = r8vec_copy_new ( n, c14 );
            break;
        case 16:
            c = r8vec_copy_new ( n, c16 );
            break;
        case 18:
            c = r8vec_copy_new ( n, c18 );
            break;
        case 20:
            c = r8vec_copy_new ( n, c20 );
            break;
        default:
            return NULL;
    }

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub2_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB2_MATRIX returns the DAUB2 matrix.
  Discussion:
    The DAUB2 matrix is the Daubechies wavelet transformation matrix
    with 2 coefficients.
    The DAUB2 matrix is also known as the Haar matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be at least 2 and a multiple of 2.
    Output, double DAUB2_MATRIX[N*N], the matrix.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a;
    ityp *c;
    dim_typ i, m;

    if ( n < 2 || ( n % 2 ) != 0 )
        return NULL;

    a = r8mat_zero_new ( n, n );
    c = daub_coefficients ( 2 );

    for ( i = 0; i < n - 1; i += 2 )
    {
        a[i+i*n]       =   c[0];
        a[i+(i+1)*n]   =   c[1];

        a[i+1+i*n]     =   c[1];
        a[i+1+(i+1)*n] = - c[0];
    }

    free ( c );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _daub2_scale ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB2_SCALE recursively evaluates the DAUB2 scaling function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the recursion level.
    Input, double X, the point at which the function is to
    be evaluated.
    Output, double DAUB2_SCALE, the estimated value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const sit * const s_data = data;
	const register short n = s_data->a0;
	const register ityp x = s_data->a1;
	
    ityp c[2] =
    {
        7.071067811865475E-01,
        7.071067811865475E-01
    };
    
    result = n>0 ? sqrt ( 2.00 ) *( c[0] * daub2_scale ( n - 1, 2.00 * x       )+ c[1] * daub2_scale ( n - 1, 2.00 * x - 1.00 ) ) : 0.00 <= x && x < 1.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub2_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB2_TRANSFORM computes the DAUB2 transform of a vector.
  Discussion:
    DAUB2 is better known as the Haar transform.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB2_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[2] =
    {
        7.071067811865475E-01,
        7.071067811865475E-01
    };
    dim_typ i, m;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );

    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        z[i] = 0.00;

    m = n;

    while ( 2 <= m )
    {
        m /= 2;

        for ( i = 0; i < m; ++i )
        {
            z[i]   = c[0] * ( y[(i<<1)] + y[(i<<1)+1] );
            z[i+m] = c[1] * ( y[(i<<1)] - y[(i<<1)+1] );
        }

        for ( i = 0; i < (m<<1); ++i )
            y[i] = z[i];
    }
    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub2_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB2_TRANSFORM_INVERSE inverts the DAUB2 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2.
    Input, double Y[N], the transformed vector.
    Output, double DAUB2_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[2] =
    {
        7.071067811865475E-01,
        7.071067811865475E-01
    };
    dim_typ i, m;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        z[i] = 0.00;

    m = 1;
    
	
	while ( (m<<1) <= n )
	{
		for ( i = 0; i < m; ++i )
		{
			z[i<<1]   = c[0] * ( x[i] + x[i+m] );
			z[(i<<1)+1] = c[1] * ( x[i] - x[i+m] );
		}
		
		for ( i = 0; i < (m<<1); ++i )
			x[i] = z[i];
		m <<= 1;
	}
	
	free ( z );
	
	return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub4_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB4_MATRIX returns the DAUB4 matrix.
  Discussion:
    The DAUB4 matrix is the Daubechies wavelet transformation matrix
    with 4 coefficients.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be at least 4 and a multiple of 2.
    Output, double DAUB4_MATRIX[N*N], the matrix.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a;
    ityp *c;
    dim_typ i, j;

    if ( n < 4 || ( n % 2 ) != 0 )
        return NULL;

    a = r8mat_zero_new ( n, n );
    c = daub_coefficients ( 4 );

    for ( i = 0; i < n - 1; i += 2 )
    {
        j = i;
        a[i+j*n] = c[0];
        j = i + 1;
        a[i+j*n] = c[1];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+j*n] = c[2];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+j*n] = c[3];

        j = i;
        a[i+1+j*n] =   c[3];
        j = i + 1;
        a[i+1+j*n] = - c[2];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+1+j*n] =   c[1];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+1+j*n] = - c[0];
    }

    free ( c );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _daub4_scale ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB4_SCALE recursively evaluates the DAUB4 scaling function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:

    John Burkardt
  Parameters:
    Input, int N, the recursion level.
    Input, double X, the point at which the function is to
    be evaluated.
    Output, double DAUB4_SCALE, the estimated value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const sit * const s_data = data;
	const register short n = s_data->a0;
	const register ityp x = s_data->a1;
	
    ityp c[4] =
    {
        0.4829629131445341E+00,
        0.8365163037378079E+00,
        0.2241438680420133E+00,
        -0.1294095225512603E+00
    };
    

    result = 0<n ? sqrt ( 2.0 ) *( c[0] * daub4_scale ( n - 1, 2.00 * x       )+ c[1] * daub4_scale ( n - 1, 2.00 * x - 1.00 )+ c[2] * daub4_scale ( n - 1, 2.00 * x - 2.00 )+ c[3] * daub4_scale ( n - 1, 2.00 * x - 3.00 ) ) :
        0.00 <= x && x < 1.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub4_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB4_TRANSFORM computes the DAUB4 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB4_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[4] =
    {
        0.4829629131445341,
        0.8365163037378079,
        0.2241438680420133,
        -0.1294095225512603
    };
    dim_typ i, j;
    int j0;
    int j1;
    int j2;
    int j3;
    int m;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        z[i] = 0.0;
    m = n;

    while ( 4 <= m )
    {
        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            j0 = i4_wrap ( j,     0, m - 1 );
            j1 = i4_wrap ( j + 1, 0, m - 1 );
            j2 = i4_wrap ( j + 2, 0, m - 1 );
            j3 = i4_wrap ( j + 3, 0, m - 1 );

            z[i]     = c[0] * y[j0] + c[1] * y[j1]+ c[2] * y[j2] + c[3] * y[j3];
            z[i+m/2] = c[3] * y[j0] - c[2] * y[j1]+ c[1] * y[j2] - c[0] * y[j3];

            ++ i;
        }

        for ( i = 0; i < m; ++i )
            y[i] = z[i];

        m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub4_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB4_TRANSFORM_INVERSE inverts the DAUB4 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB4_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[4] =
    {
        0.4829629131445341,
        0.8365163037378079,
        0.2241438680420133,
        -0.1294095225512603
    };
    dim_typ i;
    int i0;
    int i1;
    int i2;
    int i3;
    dim_typ j, m;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        z[i] = 0.00;

    m = 4;

    while ( m <= n )
    {
        j = 0;

        for ( i = 0; i < m / 2; ++i )
        {
            i0 = i4_wrap ( i - 1,          0,     m / 2 - 1 );
            i2 = i4_wrap ( i,              0,     m / 2 - 1 );

            i1 = i4_wrap ( i + m / 2 - 1,  m / 2, m - 1 );
            i3 = i4_wrap ( i + m / 2,      m / 2, m - 1 );

            z[j]   = c[2] * x[i0] + c[1] * x[i1]+ c[0] * x[i2] + c[3] * x[i3];

            z[j+1] = c[3] * x[i0] - c[0] * x[i1]+ c[1] * x[i2] - c[2] * x[i3];

            j += 2;
        }

        for ( i = 0; i < m; ++i )
            x[i] = z[i];
        m <<= 1;
    }

    free ( z );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub6_matrix ( void * data)
/*****************************************************************************/
/*
  Purpose:
    DAUB6_MATRIX returns the DAUB6 matrix.
  Discussion:
    The DAUB6 matrix is the Daubechies wavelet transformation matrix
    with 6 coefficients.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be at least 6 and a multiple of 2.
    Output, double DAUB6_MATRIX[N*N], the matrix.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a;
    ityp *c;
    dim_typ i, j;

    if ( n < 6 || ( n % 2 ) != 0 )
        return NULL;

    a = r8mat_zero_new ( n, n );
    c = daub_coefficients ( 6 );

    for ( i = 0; i < n - 1; i += 2 )
    {
        j = i;
        a[i+j*n] = c[0];
        j = i + 1;
        a[i+j*n] = c[1];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+j*n] = c[2];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+j*n] = c[3];
        j = i4_wrap ( i + 4, 0, n - 1 );
        a[i+j*n] = c[4];
        j = i4_wrap ( i + 5, 0, n - 1 );
        a[i+j*n] = c[5];

        j = i;
        a[i+1+j*n] =   c[5];
        j = i + 1;
        a[i+1+j*n] = - c[4];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+1+j*n] =   c[3];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+1+j*n] = - c[2];
        j = i4_wrap ( i + 4, 0, n - 1 );
        a[i+1+j*n] =   c[1];
        j = i4_wrap ( i + 5, 0, n - 1 );
        a[i+1+j*n] = - c[0];
    }

    free ( c );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _daub6_scale ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB6_SCALE recursively evaluates the DAUB6 scaling function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the recursion level.
    Input, double X, the point at which the function is to
    be evaluated.
    Output, double DAUB6_SCALE, the estimated value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const sit * const s_data = data;
	const register short n = s_data->a0;
	const register ityp x = s_data->a1;
	
    ityp c[6] =
    {
        0.3326705529500826E+00,
        0.8068915093110925E+00,
        0.4598775021184915E+00,
        - 0.1350110200102545E+00,
        - 0.08544127388202666E+00,
        0.03522629188570953E+00
    };
    
    result = 0<n ?  sqrt ( 2.00 ) *( c[0] * daub6_scale ( n - 1, 2.00 * x       )+ c[1] * daub6_scale ( n - 1, 2.00 * x - 1.00 )+ c[2] * daub6_scale ( n - 1, 2.00 * x - 2.00 )+ c[3] * daub6_scale ( n - 1, 2.00 * x - 3.00 )+ c[4] * daub6_scale ( n - 1, 2.00 * x - 4.00 )+ c[5] * daub6_scale ( n - 1, 2.00 * x - 5.00 ) ) :0.00 <= x && x < 1.00;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub6_transform ( void * data)
/*******************************************************************************/
/*
  Purpose:
    DAUB6_TRANSFORM computes the DAUB6 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB6_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[6] =
    {
        0.3326705529500826,
        0.8068915093110925,
        0.4598775021184915,
        - 0.1350110200102545,
        - 0.08544127388202666,
        0.03522629188570953
    };
    dim_typ i, j;
    int j0;
    int j1;
    int k;
    dim_typ m;
    dim_typ p = 5;
    dim_typ q;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = n;
    q = ( p - 1 ) / 2;

    while ( 4 <= m )
    {
        for ( i = 0; i < m; ++i )
            z[i] = 0.00;

        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            for ( k = 0; k < p; k += 2 )
            {
                j0 = i4_wrap ( j + k,     0, m - 1 );
                j1 = i4_wrap ( j + k + 1, 0, m - 1 );
                z[i]     = z[i]     + c[  k] * y[j0] + c[  k+1] * y[j1];
                z[i+m/2] = z[i+m/2] + c[p-k] * y[j0] - c[p-k-1] * y[j1];
            }
            ++ i;
        }

        for ( i = 0; i < m; ++i )
            y[i] = z[i];
        m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub6_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB6_TRANSFORM_INVERSE inverts the DAUB6 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB6_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[6] =
    {
        0.3326705529500826,
        0.8068915093110925,
        0.4598775021184915,
        - 0.1350110200102545,
        - 0.08544127388202666,
        0.03522629188570953
    };
    dim_typ i;
    int i0;
    int i1;
    dim_typ j;
    dim_typ k;
    dim_typ m;
    dim_typ p = 5;
    dim_typ q;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = 4;
    q = ( p - 1 ) / 2;

    while ( m <= n )
    {
        for ( i = 0; i < n; ++i )
            z[i] = 0.0;
        j = 0;

        for ( i = - q; i < m / 2 - q; ++i )
        {
            for ( k = 0; k < p; k += 2 )
            {
                i0 = i4_wrap ( i         + k / 2,     0,     m / 2 - 1 );
                i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2, m     - 1 );
                z[j]   = z[j]   + c[p-k-1] * x[i0] + c[k+1] * x[i1];
                z[j+1] = z[j+1] + c[p-k]   * x[i0] - c[k]   * x[i1];
            }
            j += 2;
        }

        for ( i = 0; i < m; ++i )
            x[i] = z[i];
        m <<= 1;
    }

    free ( z );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub8_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB8_MATRIX returns the DAUB8 matrix.
  Discussion:
    The DAUB8 matrix is the Daubechies wavelet transformation matrix
    with 8 coefficients.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be at least 8 and a multiple of 2.
    Output, double DAUB8_MATRIX[N*N], the matrix.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a;
    ityp *c;
    dim_typ i, j;

    if ( n < 8 || ( n % 2 ) != 0 )
        return NULL;

    a = r8mat_zero_new ( n, n );
    c = daub_coefficients ( 8 );

    for ( i = 0; i < n - 1; i += 2 )
    {
        j = i;
        a[i+j*n] = c[0];
        j = i + 1;
        a[i+j*n] = c[1];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+j*n] = c[2];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+j*n] = c[3];
        j = i4_wrap ( i + 4, 0, n - 1 );
        a[i+j*n] = c[4];
        j = i4_wrap ( i + 5, 0, n - 1 );
        a[i+j*n] = c[5];
        j = i4_wrap ( i + 6, 0, n - 1 );
        a[i+j*n] = c[6];
        j = i4_wrap ( i + 7, 0, n - 1 );
        a[i+j*n] = c[7];

        j = i;
        a[i+1+j*n] =   c[7];
        j = i + 1;
        a[i+1+j*n] = - c[6];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+1+j*n] =   c[5];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+1+j*n] = - c[4];
        j = i4_wrap ( i + 4, 0, n - 1 );
        a[i+1+j*n] =   c[3];
        j = i4_wrap ( i + 5, 0, n - 1 );
        a[i+1+j*n] = - c[2];
        j = i4_wrap ( i + 6, 0, n - 1 );
        a[i+1+j*n] =   c[1];
        j = i4_wrap ( i + 7, 0, n - 1 );
        a[i+1+j*n] = - c[0];
    }

    free ( c );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _daub8_scale ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB8_SCALE recursively evaluates the DAUB8 scaling function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the recursion level.
    Input, double X, the point at which the function is to
    be evaluated.
    Output, double DAUB8_SCALE, the estimated value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const sit * const s_data = data;
	const register short n = s_data->a0;
	const register ityp x = s_data->a1;
	
    ityp c[8] =
    {
        0.2303778133088964E+00,
        0.7148465705529154E+00,
        0.6308807679298587E+00,
        -0.0279837694168599E+00,
        -0.1870348117190931E+00,
        0.0308413818355607E+00,
        0.0328830116668852E+00,
        -0.0105974017850690E+00
    };
    result = 0<n? sqrt ( 2.00 ) *
  ( c[0] * daub8_scale ( n - 1, 2.00 * x       )
        + c[1] * daub8_scale ( n - 1, 2.00 * x - 1.00 )
        + c[2] * daub8_scale ( n - 1, 2.00 * x - 2.00 )
        + c[3] * daub8_scale ( n - 1, 2.00 * x - 3.00 )
        + c[4] * daub8_scale ( n - 1, 2.00 * x - 4.00 )
        + c[5] * daub8_scale ( n - 1, 2.00 * x - 5.00 )
        + c[6] * daub8_scale ( n - 1, 2.00 * x - 6.00 )
        + c[7] * daub8_scale ( n - 1, 2.00 * x - 7.00 ) ) :
            0.00 <= x && x < 1.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub8_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB8_TRANSFORM computes the DAUB8 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB8_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[8] =
    {
        0.2303778133088964,
        0.7148465705529154,
        0.6308807679298587,
        -0.02798376941685985,
        -0.1870348117190931,
        0.03084138183556076,
        0.03288301166688519,
        -0.01059740178506903
    };
    dim_typ i, j;
    int j0;
    int j1;
    dim_typ k;
    dim_typ m;
    dim_typ p = 7;
    dim_typ q;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = n;
    q = ( p - 1 ) / 2;

    while ( 4 <= m )
    {
        for ( i = 0; i < m; ++i)
            z[i] = 0.00;
        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            for ( k = 0; k < p; k += 2 )
            {
                j0 = i4_wrap ( j + k,     0, m - 1 );
                j1 = i4_wrap ( j + k + 1, 0, m - 1 );
                z[i]     = z[i]     + c[  k] * y[j0] + c[  k+1] * y[j1];
                z[i+m/2] = z[i+m/2] + c[p-k] * y[j0] - c[p-k-1] * y[j1];
            }
                ++ i;
        }

        for ( i = 0; i < m; ++i )
            y[i] = z[i];
        m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub8_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB8_TRANSFORM_INVERSE inverts the DAUB8 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 April 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB8_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[8] =
    {
        0.2303778133088964,
        0.7148465705529154,
        0.6308807679298587,
        -0.02798376941685985,
        -0.1870348117190931,
        0.03084138183556076,
        0.03288301166688519,
        -0.01059740178506903
    };
    dim_typ i;
    int i0;
    int i1;
    dim_typ j;
    dim_typ k;
    dim_typ m;
    dim_typ p = 7;
    dim_typ q;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = 4;
    q = ( p - 1 ) / 2;

    while ( m <= n )
    {
        for ( i = 0; i < n; ++i )
            z[i] = 0.00;

        j = 0;

        for ( i = - q; i < m / 2 - q; ++i )
        {
            for ( k = 0; k < p; k += 2 )
            {
                i0 = i4_wrap ( i         + k / 2,     0,     m / 2 - 1 );
                i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2, m     - 1 );
                z[j]   = z[j]   + c[p-k-1] * x[i0] + c[k+1] * x[i1];
                z[j+1] = z[j+1] + c[p-k]   * x[i0] - c[k]   * x[i1];
            }
            j += 2;
        }

        for ( i = 0; i < m; ++i )
            x[i] = z[i];
        m <<= 1;
    }

    free ( z );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub10_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB10_MATRIX returns the DAUB10 matrix.
  Discussion:
    The DAUB10 matrix is the Daubechies wavelet transformation matrix
    with 10 coefficients.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be at least 10 and a multiple of 2.
    Output, double DAUB10_MATRIX[N*N], the matrix.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a;
    ityp *c;
    dim_typ i, j;

    if ( n < 10 || ( n % 2 ) != 0 )
        return NULL;

    a = r8mat_zero_new ( n, n );

    c = daub_coefficients ( 10 );

    for ( i = 0; i < n - 1; i += 2 )
    {
        j = i;
        a[i+j*n] = c[0];
        j = i + 1;
        a[i+j*n] = c[1];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+j*n] = c[2];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+j*n] = c[3];
        j = i4_wrap ( i + 4, 0, n - 1 );
        a[i+j*n] = c[4];
        j = i4_wrap ( i + 5, 0, n - 1 );
        a[i+j*n] = c[5];
        j = i4_wrap ( i + 6, 0, n - 1 );
        a[i+j*n] = c[6];
        j = i4_wrap ( i + 7, 0, n - 1 );
        a[i+j*n] = c[7];
        j = i4_wrap ( i + 8, 0, n - 1 );
        a[i+j*n] = c[8];
        j = i4_wrap ( i + 9, 0, n - 1 );
        a[i+j*n] = c[9];

        j = i;
        a[i+1+j*n] =   c[9];
        j = i + 1;
        a[i+1+j*n] = - c[8];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+1+j*n] =   c[7];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+1+j*n] = - c[6];
        j = i4_wrap ( i + 4, 0, n - 1 );
        a[i+1+j*n] =   c[5];
        j = i4_wrap ( i + 5, 0, n - 1 );
        a[i+1+j*n] = - c[4];
        j = i4_wrap ( i + 6, 0, n - 1 );
        a[i+1+j*n] =   c[3];
        j = i4_wrap ( i + 7, 0, n - 1 );
        a[i+1+j*n] = - c[2];
        j = i4_wrap ( i + 8, 0, n - 1 );
        a[i+1+j*n] =   c[1];
        j = i4_wrap ( i + 9, 0, n - 1 );
        a[i+1+j*n] = - c[0];
    }

    free ( c );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _daub10_scale ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB10_SCALE recursively evaluates the DAUB10 scaling function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the recursion level.
    Input, double X, the point at which the function is to
    be evaluated.
    Output, double DAUB10_SCALE, the estimated value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const sit * const s_data = data;
	const register short n = s_data->a0;
	const register ityp x = s_data->a1;
	
    ityp c[10] =
    {
        0.1601023979741929E+00,
        0.6038292697971895E+00,
        0.7243085284377726E+00,
        0.1384281459013203E+00,
        -0.2422948870663823E+00,
        -0.0322448695846381E+00,
        0.0775714938400459E+00,
        -0.0062414902127983E+00,
        -0.0125807519990820E+00,
        0.0033357252854738E+00
    };
    result = 0<n ?  sqrt ( 2.00 ) *
 ( c[0] * daub10_scale ( n - 1, 2.00 * x       )
    + c[1] * daub10_scale ( n - 1, 2.00 * x - 1.00 )
    + c[2] * daub10_scale ( n - 1, 2.00 * x - 2.00 )
    + c[3] * daub10_scale ( n - 1, 2.00 * x - 3.00 )
    + c[4] * daub10_scale ( n - 1, 2.00 * x - 4.00 )
    + c[5] * daub10_scale ( n - 1, 2.00 * x - 5.00 )
    + c[6] * daub10_scale ( n - 1, 2.00 * x - 6.00 )
    + c[7] * daub10_scale ( n - 1, 2.00 * x - 7.00 )
    + c[8] * daub10_scale ( n - 1, 2.00 * x - 8.00 )
    + c[9] * daub10_scale ( n - 1, 2.00 * x - 9.00 ) ) :
        0.00 <= x && x < 1.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub10_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB10_TRANSFORM computes the DAUB10 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB10_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[10] =
    {
        1.601023979741929E-01,
        6.038292697971896E-01,
        7.243085284377729E-01,
        1.384281459013207E-01,
        -2.422948870663820E-01,
        -3.224486958463837E-02,
        7.757149384004571E-02,
        -6.241490212798274E-03,
        -1.258075199908199E-02,
        3.335725285473771E-03
    };
    dim_typ i, j;
    int j0;
    int j1;
    dim_typ k;
    dim_typ m;
    dim_typ p = 9;
    dim_typ q;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = n;
    q = ( p - 1 ) / 2;

    while ( 4 <= m )
    {
        for ( i = 0; i < m; ++i )
            z[i] = 0.0;
        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            for ( k = 0; k < p; k += 2 )
            {
                j0 = i4_wrap ( j + k,     0, m - 1 );
                j1 = i4_wrap ( j + k + 1, 0, m - 1 );
                z[i]     = z[i]     + c[  k] * y[j0] + c[  k+1] * y[j1];
                z[i+m/2] = z[i+m/2] + c[p-k] * y[j0] - c[p-k-1] * y[j1];
            }
            ++ i;
        }

        for ( i = 0; i < m; ++i)
            y[i] = z[i];

        m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub10_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB10_TRANSFORM_INVERSE inverts the DAUB10 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB10_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[10] =
    {
        1.601023979741929E-01,
        6.038292697971896E-01,
        7.243085284377729E-01,
        1.384281459013207E-01,
        -2.422948870663820E-01,
        -3.224486958463837E-02,
        7.757149384004571E-02,
        -6.241490212798274E-03,
        -1.258075199908199E-02,
        3.335725285473771E-03
    };
    dim_typ i;
    int i0;
    int i1;
    dim_typ j;
    dim_typ k;
    dim_typ m;
    dim_typ p = 9;
    ityp q;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = 4;
    q = ( p - 1 ) / 2;

    while ( m <= n )
    {
        for ( i = 0; i < n; ++i)
            z[i] = 0.00;

        j = 0;

        for ( i = - q; i < m / 2 - q; ++i )
        {
            for ( k = 0; k < p; k += 2 )
            {
                i0 = i4_wrap ( i         + k / 2,     0,     m / 2 - 1 );
                i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2, m     - 1 );
                z[j]   = z[j]   + c[p-k-1] * x[i0] + c[k+1] * x[i1];
                z[j+1] = z[j+1] + c[p-k]   * x[i0] - c[k]   * x[i1];
            }
            j += 2;
        }

        for ( i = 0; i < m; ++i)
            x[i] = z[i];
        m <<= 1;
    }

    free ( z );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub12_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB12_MATRIX returns the DAUB12 matrix.
  Discussion:
    The DAUB12 matrix is the Daubechies wavelet transformation matrix
    with 12 coefficients.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    N must be at least 12 and a multiple of 2.
    Output, double DAUB12_MATRIX[N*N], the matrix.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a;
    ityp *c;
    dim_typ i, j;

    if ( n < 12 || ( n % 2 ) != 0 )
        return NULL;

    a = r8mat_zero_new ( n, n );
    c = daub_coefficients ( 12 );

    for ( i = 0; i < n - 1; i += 2 )
    {
        j = i;
        a[i+j*n] = c[0];
        j = i + 1;
        a[i+j*n] = c[1];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+j*n] = c[2];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+j*n] = c[3];
        j = i4_wrap ( i + 4, 0, n - 1 );
        a[i+j*n] = c[4];
        j = i4_wrap ( i + 5, 0, n - 1 );
        a[i+j*n] = c[5];
        j = i4_wrap ( i + 6, 0, n - 1 );
        a[i+j*n] = c[6];
        j = i4_wrap ( i + 7, 0, n - 1 );
        a[i+j*n] = c[7];
        j = i4_wrap ( i + 8, 0, n - 1 );
        a[i+j*n] = c[8];
        j = i4_wrap ( i + 9, 0, n - 1 );
        a[i+j*n] = c[9];
        j = i4_wrap ( i + 10, 0, n - 1 );
        a[i+j*n] = c[10];
        j = i4_wrap ( i + 11, 0, n - 1 );
        a[i+j*n] = c[11];

        j = i;
        a[i+1+j*n] =   c[11];
        j = i + 1;
        a[i+1+j*n] = - c[10];
        j = i4_wrap ( i + 2, 0, n - 1 );
        a[i+1+j*n] =   c[9];
        j = i4_wrap ( i + 3, 0, n - 1 );
        a[i+1+j*n] = - c[8];
        j = i4_wrap ( i + 4, 0, n - 1 );
        a[i+1+j*n] =   c[7];
        j = i4_wrap ( i + 5, 0, n - 1 );
        a[i+1+j*n] = - c[6];
        j = i4_wrap ( i + 6, 0, n - 1 );
        a[i+1+j*n] =   c[5];
        j = i4_wrap ( i + 7, 0, n - 1 );
        a[i+1+j*n] = - c[4];
        j = i4_wrap ( i + 8, 0, n - 1 );
        a[i+1+j*n] =   c[3];
        j = i4_wrap ( i + 9, 0, n - 1 );
        a[i+1+j*n] = - c[2];
        j = i4_wrap ( i + 10, 0, n - 1 );
        a[i+1+j*n] =   c[1];
        j = i4_wrap ( i + 11, 0, n - 1 );
        a[i+1+j*n] = - c[0];
    }

    free ( c );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub12_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB12_TRANSFORM computes the DAUB12 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB12_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[12] =
    {
        0.1115407433501095E+00,
        0.4946238903984533E+00,
        0.7511339080210959E+00,
        0.3152503517091982E+00,
        -0.2262646939654400E+00,
        -0.1297668675672625E+00,
        0.0975016055873225E+00,
        0.0275228655303053E+00,
        -0.0315820393174862E+00,
        0.0005538422011614E+00,
        0.0047772575109455E+00,
        -0.0010773010853085E+00
    };
    dim_typ i, j;
    int j0;
    int j1;
    dim_typ k;
    dim_typ m;
    dim_typ p = 11;
    dim_typ q;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = n;
    q = ( p - 1 ) / 2;

    while ( 4 <= m )
    {
        for ( i = 0; i < m; ++i )
            z[i] = 0.00;
        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            for ( k = 0; k < p; k += 2 )
            {
                j0 = i4_wrap ( j + k,     0, m - 1 );
                j1 = i4_wrap ( j + k + 1, 0, m - 1 );
                z[i]     = z[i]     + c[  k] * y[j0] + c[  k+1] * y[j1];
                z[i+m/2] = z[i+m/2] + c[p-k] * y[j0] - c[p-k-1] * y[j1];
            }
            ++ i;
        }

        for ( i = 0; i < m; ++i)
            y[i] = z[i];

 	   m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub12_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB12_TRANSFORM_INVERSE inverts the DAUB12 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB12_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[12] =
    {
        0.1115407433501095E+00,
        0.4946238903984533E+00,
        0.7511339080210959E+00,
        0.3152503517091982E+00,
        -0.2262646939654400E+00,
        -0.1297668675672625E+00,
        0.0975016055873225E+00,
        0.0275228655303053E+00,
        -0.0315820393174862E+00,
        0.0005538422011614E+00,
        0.0047772575109455E+00,
        -0.0010773010853085E+00
    };
    dim_typ i;
    int i0;
    int i1;
    dim_typ j;
    dim_typ k;
    dim_typ m;
    dim_typ p = 11;
    dim_typ q;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = 4;
    q = ( p - 1 ) / 2;

    while ( m <= n )
    {
        for ( i = 0; i < n; ++i )
            z[i] = 0.00;

        j = 0;

        for ( i = - q; i < m / 2 - q; ++i )
        {
            for ( k = 0; k < p; k += 2 )
            {
                i0 = i4_wrap ( i         + k / 2,     0,     m / 2 - 1 );
                i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2, m     - 1 );
                z[j]   = z[j]   + c[p-k-1] * x[i0] + c[k+1] * x[i1];
                z[j+1] = z[j+1] + c[p-k]   * x[i0] - c[k]   * x[i1];
            }
            j += 2;
        }

        for ( i = 0; i < m; ++i )
            x[i] = z[i];
        m <<= 1;
    }

    free ( z );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub14_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB14_TRANSFORM computes the DAUB14 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB14_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[14] =
    {
        7.785205408500917E-02,
        3.965393194819173E-01,
        7.291320908462351E-01,
        4.697822874051931E-01,
        -1.439060039285649E-01,
        -2.240361849938749E-01,
        7.130921926683026E-02,
        8.061260915108307E-02,
        -3.802993693501441E-02,
        -1.657454163066688E-02,
        1.255099855609984E-02,
        4.295779729213665E-04,
        -1.801640704047490E-03,
        3.537137999745202E-04
    };
    dim_typ i, j;
    int j0;
    int j1;
    dim_typ k;
    dim_typ m;
    dim_typ p = 13;
    dim_typ q;
    double *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = n;
    q = ( p - 1 ) / 2;

    while ( 4 <= m )
    {
        for ( i = 0; i < m; ++i)
            z[i] = 0.00;
        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            for ( k = 0; k < p; k += 2 )
            {
                j0 = i4_wrap ( j + k,     0, m - 1 );
                j1 = i4_wrap ( j + k + 1, 0, m - 1 );
                z[i]     = z[i]     + c[  k] * y[j0] + c[  k+1] * y[j1];
                z[i+m/2] = z[i+m/2] + c[p-k] * y[j0] - c[p-k-1] * y[j1];
            }
            ++ i;
        }

        for ( i = 0; i < m; ++i )
            y[i] = z[i];

        m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub14_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB14_TRANSFORM_INVERSE inverts the DAUB14 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB14_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[14] =
    {
        7.785205408500917E-02,
        3.965393194819173E-01,
        7.291320908462351E-01,
        4.697822874051931E-01,
        -1.439060039285649E-01,
        -2.240361849938749E-01,
        7.130921926683026E-02,
        8.061260915108307E-02,
        -3.802993693501441E-02,
        -1.657454163066688E-02,
        1.255099855609984E-02,
        4.295779729213665E-04,
        -1.801640704047490E-03,
        3.537137999745202E-04
    };
    dim_typ i;
    int i0;
    int i1;
    dim_typ j;
    dim_typ k;
    dim_typ m;
    dim_typ p = 13;
    dim_typ q;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = 4;
    q = ( p - 1 ) / 2;

    while ( m <= n )
    {
        for ( i = 0; i < n; ++i )
            z[i] = 0.00;

        j = 0;

        for ( i = - q; i < m / 2 - q; ++i )
        {
            for ( k = 0; k < p; k += 2 )
            {
                i0 = i4_wrap ( i         + k / 2,     0,     m / 2 - 1 );
                i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2, m     - 1 );
                z[j]   = z[j]   + c[p-k-1] * x[i0] + c[k+1] * x[i1];
                z[j+1] = z[j+1] + c[p-k]   * x[i0] - c[k]   * x[i1];
            }
            j += 2;
        }

        for ( i = 0; i < m; ++i)
            x[i] = z[i];
        m <<= 1;
    }

    free ( z );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub16_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB16_TRANSFORM computes the DAUB16 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB16_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[16] =
    {
        5.441584224310400E-02,
        3.128715909142999E-01,
        6.756307362972898E-01,
        5.853546836542067E-01,
        -1.582910525634930E-02,
        -2.840155429615469E-01,
        4.724845739132827E-04,
        1.287474266204784E-01,
        -1.736930100180754E-02,
        -4.408825393079475E-02,
        1.398102791739828E-02,
        8.746094047405776E-03,
        -4.870352993451574E-03,
        -3.917403733769470E-04,
        6.754494064505693E-04,
        -1.174767841247695E-04
    };
    dim_typ i, j;
    int j0;
    int j1;
    dim_typ k;
    dim_typ m;
    dim_typ p = 15;
    dim_typ q;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = n;
    q = ( p - 1 ) / 2;

    while ( 4 <= m )
    {
        for ( i = 0; i < m; ++i )
            z[i] = 0.0;
        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            for ( k = 0; k < p; k += 2 )
            {
                j0 = i4_wrap ( j + k,     0, m - 1 );
                j1 = i4_wrap ( j + k + 1, 0, m - 1 );
                z[i]     = z[i]     + c[  k] * y[j0] + c[  k+1] * y[j1];
                z[i+m/2] = z[i+m/2] + c[p-k] * y[j0] - c[p-k-1] * y[j1];
            }
            ++ i;
        }

        for ( i = 0; i < m; ++i )
            y[i] = z[i];

        m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub16_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB16_TRANSFORM_INVERSE inverts the DAUB16 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB16_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[16] =
    {
        5.441584224310400E-02,
        3.128715909142999E-01,
        6.756307362972898E-01,
        5.853546836542067E-01,
        -1.582910525634930E-02,
        -2.840155429615469E-01,
        4.724845739132827E-04,
        1.287474266204784E-01,
        -1.736930100180754E-02,
        -4.408825393079475E-02,
        1.398102791739828E-02,
        8.746094047405776E-03,
        -4.870352993451574E-03,
        -3.917403733769470E-04,
        6.754494064505693E-04,
        -1.174767841247695E-04
    };
    dim_typ i;
    int i0;
    int i1;
    dim_typ j;
    dim_typ k;
    dim_typ m;
    dim_typ p = 15;
    dim_typ q;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = 4;
    q = ( p - 1 ) / 2;

    while ( m <= n )
    {
        for ( i = 0; i < n; ++i)
            z[i] = 0.0;

        j = 0;

        for ( i = - q; i < m / 2 - q; ++i )
        {
            for ( k = 0; k < p; k += 2 )
            {
                i0 = i4_wrap ( i         + k / 2,     0,     m / 2 - 1 );
                i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2, m     - 1 );
                z[j]   = z[j]   + c[p-k-1] * x[i0] + c[k+1] * x[i1];
                z[j+1] = z[j+1] + c[p-k]   * x[i0] - c[k]   * x[i1];
            }
            j += 2;
        }

        for ( i = 0; i < m; ++i)
            x[i] = z[i];
        m <<= 1;
    }

    free ( z );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub18_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB18_TRANSFORM computes the DAUB18 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB18_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[18] =
    {
        3.807794736387834E-02,
        2.438346746125903E-01,
        6.048231236901111E-01,
        6.572880780513005E-01,
        1.331973858250075E-01,
        -2.932737832791749E-01,
        -9.684078322297646E-02,
        1.485407493381063E-01,
        3.072568147933337E-02,
        -6.763282906132997E-02,
        2.509471148314519E-04,
        2.236166212367909E-02,
        -4.723204757751397E-03,
        -4.281503682463429E-03,
        1.847646883056226E-03,
        2.303857635231959E-04,
        -2.519631889427101E-04,
        3.934732031627159E-05
    };
    dim_typ i, j;
    int j0;
    int j1;
    dim_typ k;
    dim_typ m;
    dim_typ p = 17;
    dim_typ q;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = n;
    q = ( p - 1 ) / 2;

    while ( 4 <= m )
    {
        for ( i = 0; i < m; ++i)
            z[i] = 0.00;
        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            for ( k = 0; k < p; k += 2 )
            {
                j0 = i4_wrap ( j + k,     0, m - 1 );
                j1 = i4_wrap ( j + k + 1, 0, m - 1 );
                z[i]     = z[i]     + c[  k] * y[j0] + c[  k+1] * y[j1];
                z[i+m/2] = z[i+m/2] + c[p-k] * y[j0] - c[p-k-1] * y[j1];
            }
            ++ i;
        }

        for ( i = 0; i < m; ++i )
            y[i] = z[i];

        m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub18_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB18_TRANSFORM_INVERSE inverts the DAUB18 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB18_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[18] =
    {
        3.807794736387834E-02,
        2.438346746125903E-01,
        6.048231236901111E-01,
        6.572880780513005E-01,
        1.331973858250075E-01,
        -2.932737832791749E-01,
        -9.684078322297646E-02,
        1.485407493381063E-01,
        3.072568147933337E-02,
        -6.763282906132997E-02,
        2.509471148314519E-04,
        2.236166212367909E-02,
        -4.723204757751397E-03,
        -4.281503682463429E-03,
        1.847646883056226E-03,
        2.303857635231959E-04,
        -2.519631889427101E-04,
        3.934732031627159E-05
    };
    dim_typ i;
    int i0;
    int i1;
    dim_typ j;
    dim_typ k;
    dim_typ m;
    dim_typ p = 17;
    dim_typ q;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = 4;
    q = ( p - 1 ) / 2;

    while ( m <= n )
    {
        for ( i = 0; i < n; ++i )
            z[i] = 0.00;

        j = 0;

        for ( i = - q; i < m / 2 - q; ++i )
        {
            for ( k = 0; k < p; k += 2 )
            {
                i0 = i4_wrap ( i         + k / 2,     0,     m / 2 - 1 );
                i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2, m     - 1 );
                z[j]   = z[j]   + c[p-k-1] * x[i0] + c[k+1] * x[i1];
                z[j+1] = z[j+1] + c[p-k]   * x[i0] - c[k]   * x[i1];
            }
            j += 2;
        }

        for ( i = 0; i < m; ++i )
            x[i] = z[i];

        m <<= 1;
    }

    free ( z );

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub20_transform ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB20_TRANSFORM computes the DAUB20 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double X[N], the vector to be transformed.
    Output, double DAUB20_TRANSFORM[N], the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp c[20] =
    {
        2.667005790055555E-02,
        1.881768000776914E-01,
        5.272011889317255E-01,
        6.884590394536035E-01,
        2.811723436605774E-01,
        -2.498464243273153E-01,
        -1.959462743773770E-01,
        1.273693403357932E-01,
        9.305736460357235E-02,
        -7.139414716639708E-02,
        -2.945753682187581E-02,
        3.321267405934100E-02,
        3.606553566956169E-03,
        -1.073317548333057E-02,
        1.395351747052901E-03,
        1.992405295185056E-03,
        -6.858566949597116E-04,
        -1.164668551292854E-04,
        9.358867032006959E-05,
        -1.326420289452124E-05
    };
    dim_typ i, j;
    int j0;
    int j1;
    dim_typ k;
    dim_typ m;
    dim_typ p = 19;
    dim_typ q;
    ityp *y;
    ityp *z;

    y = r8vec_copy_new ( n, x );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = n;
    q = ( p - 1 ) / 2;

    while ( 4 <= m )
    {
        for ( i = 0; i < m; ++i )
            z[i] = 0.0;
        i = 0;

        for ( j = 0; j < m - 1; j += 2 )
        {
            for ( k = 0; k < p; k += 2 )
            {
                j0 = i4_wrap ( j + k,     0, m - 1 );
                j1 = i4_wrap ( j + k + 1, 0, m - 1 );
                z[i]     = z[i]     + c[  k] * y[j0] + c[  k+1] * y[j1];
                z[i+m/2] = z[i+m/2] + c[p-k] * y[j0] - c[p-k-1] * y[j1];
            }
            ++ i;
        }

        for ( i = 0; i < m; ++i)
            y[i] = z[i];

        m /= 2;
    }

    free ( z );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _daub20_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAUB20_TRANSFORM_INVERSE inverts the DAUB20 transform of a vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    N must be a power of 2 and at least 4.
    Input, double Y[N], the transformed vector.
    Output, double DAUB20_TRANSFORM_INVERSE[N], the original vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * y = s_data->a1;
	
    ityp c[20] =
    {
        2.667005790055555E-02,
        1.881768000776914E-01,
        5.272011889317255E-01,
        6.884590394536035E-01,
        2.811723436605774E-01,
        -2.498464243273153E-01,
        -1.959462743773770E-01,
        1.273693403357932E-01,
        9.305736460357235E-02,
        -7.139414716639708E-02,
        -2.945753682187581E-02,
        3.321267405934100E-02,
        3.606553566956169E-03,
        -1.073317548333057E-02,
        1.395351747052901E-03,
        1.992405295185056E-03,
        -6.858566949597116E-04,
        -1.164668551292854E-04,
        9.358867032006959E-05,
        -1.326420289452124E-05
    };
    dim_typ i;
    int i0;
    int i1;
    dim_typ j;
    dim_typ k;
    dim_typ m;
    dim_typ p = 19;
    dim_typ q;
    ityp *x;
    ityp *z;

    x = r8vec_copy_new ( n, y );
    z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    m = 4;
    q = ( p - 1 ) / 2;

    while ( m <= n )
    {
        for ( i = 0; i < n; ++i )
            z[i] = 0.00;
        j = 0;

        for ( i = - q; i < m / 2 - q; ++i)
        {
            for ( k = 0; k < p; k += 2 )
            {
                i0 = i4_wrap ( i         + k / 2,     0,     m / 2 - 1 );
                i1 = i4_wrap ( i + m / 2 + k / 2,     m / 2, m     - 1 );
                z[j]   = z[j]   + c[p-k-1] * x[i0] + c[k+1] * x[i1];
                z[j+1] = z[j+1] + c[p-k]   * x[i0] - c[k]   * x[i1];
            }
            j += 2;
        }

        for ( i = 0; i < m; ++i )
            x[i] = z[i];
        m <<= 1;
    }

    free ( z );

    return x;
}

#endif
