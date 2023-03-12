#ifndef __DISABLEDEEP_HAAR

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _haar_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    HAAR_1D computes the Haar transform of a vector.
  Discussion:
    For the classical Haar transform, N should be a power of 2.
    However, this is not required here.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    Input/output, double X[N], on input, the vector to be transformed.
    On output, the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, k;
    ityp s;
    ityp *y;

    s = sqrt ( 2.00 );
    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    /*
    Initialize.
    */
    for ( i = 0; i < n; ++i )
        y[i] = 0.00;
    /*
    Determine K, the largest power of 2 such that K <= N.
    */
    k = 1;
    while ( (k<<1) <= n )
        k *= 2;

    while ( 1 < k )
    {
        k /= 2;
        for ( i = 0; i < k; ++i)
        {
            y[i]   = ( x[i<<1] + x[(i<<1)+1] ) / s;
            y[i+k] = ( x[i<<1] - x[(i<<1)+1] ) / s;
        }
        for ( i = 0; i < (k<<1); ++i )
            x[i] = y[i];
    }

    free ( y );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _haar_1d_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    HAAR_1D_INVERSE computes the inverse Haar transform of a vector.
  Discussion:
    For the classical Haar transform, N should be a power of 2.
    However, this is not required here.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    Input/output, double X[N], on input, the vector to be transformed.
    On output, the transformed vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, k;
    ityp s;
    ityp *y;

    s = sqrt ( 2.00 );
    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    /*
    Initialize.
    */
    for ( i = 0; i < n; ++i )
        y[i] = 0.00;

    k = 1;
    while ( (k<<1) <= n )
    {
        for ( i = 0; i < k; ++i)
        {
            y[i<<1]   = ( x[i] + x[i+k] ) / s;
            y[(i<<1)+1] = ( x[i] - x[i+k] ) / s;
        }
        for ( i = 0; i < (k<<1); ++i )
            x[i] = y[i];
        k *= 2;
    }

    free ( y );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _haar_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    HAAR_2D computes the Haar transform of an array.
  Discussion:
    For the classical Haar transform, M and N should be a power of 2.
    However, this is not required here.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the dimensions of the array.
    Input/output, double U[M*N], the array to be transformed.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1; 
	ityp * u = s_data->a2;
	
    dim_typ i, j, k;
    ityp s;
    ityp *v;

    s = sqrt ( 2.00 );
    v = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            v[i+j*m] = u[i+j*m];
    /*
    Determine K, the largest power of 2 such that K <= M.
    */
    k = 1;
    while ( (k<<1) <= m )
        k *= 2;
    /*
    Transform all columns.
    */
    while ( 1 < k )
    {
        k /= 2;

        for ( j = 0; j < n; ++j )
        {
            for ( i = 0; i < k; ++i )
            {
                v[i  +j*m] = ( u[(i<<1)+j*m] + u[(i<<1)+1+j*m] ) / s;
                v[k+i+j*m] = ( u[(i<<1)+j*m] - u[(i<<1)+1+j*m] ) / s;
            }
    }
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < (i<<1); ++i )
            u[i+j*m] = v[i+j*m];
    }
    /*
    Determine K, the largest power of 2 such that K <= N.
    */
    k = 1;
    while ( (k<<1) <= n )
        k *= 2;
    /*
    Transform all rows.
    */
    while ( 1 < k )
    {
        k /= 2;

        for ( j = 0; j < k; ++j)
            for ( i = 0; i < m; ++i )
            {
                v[i+(  j)*m] = ( u[i+(j<<1)*m] + u[i+((j<<1)+1)*m] ) / s;
                v[i+(k+j)*m] = ( u[i+(j<<1)*m] - u[i+((j<<1)+1)*m] ) / s;
            }

        for ( j = 0; j < (k<<1); ++j )
            for ( i = 0; i < m; ++i )
                u[i+j*m] = v[i+j*m];
    }
    free ( v );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _haar_2d_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    HAAR_2D_INVERSE inverts the Haar transform of an array.
  Discussion:
    For the classical Haar transform, M and N should be a power of 2.
    However, this is not required here.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the dimensions of the array.
    Input/output, double U[M*N], the array to be transformed.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1; 
	ityp * u = s_data->a2;
	
    dim_typ i, j, k;
    ityp s;
    ityp *v;

    s = sqrt ( 2.00 );
    v = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            v[i+j*m] = u[i+j*m];
    /*
    Inverse transform of all rows.
    */
    k = 1;

    while ( (k<<1) <= n )
    {
        for ( j = 0; j < k; ++j )
            for ( i = 0; i < m; ++i )
            {
                v[i+(j<<1  )*m] = ( u[i+j*m] + u[i+(k+j)*m] ) / s;
                v[i+((j<<1)+1)*m] = ( u[i+j*m] - u[i+(k+j)*m] ) / s;
            }

        for ( j = 0; j < (k<<1); ++j )
            for ( i = 0; i < m; ++i )
                u[i+j*m] = v[i+j*m];
        k *= 2;
    }
    /*
    Inverse transform of all columns.
    */
    k = 1;

    while ( (k<<1) <= m )
    {
        for ( j = 0; j < n; ++j )
            for ( i = 0; i < k; ++i )
            {
                v[(i<<1)  +j*m] = ( u[i+j*m] + u[k+i+j*m] ) / s;
                v[(i<<1)+1+j*m] = ( u[i+j*m] - u[k+i+j*m] ) / s;
            }

        for ( j = 0; j < n; ++j )
            for ( i = 0; i < (k<<1); ++i )
                u[i+j*m] = v[i+j*m];
        k <<= 1;
    }
    free ( v );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _dif_deriv ( void * data)
/******************************************************************************/
/*
  Purpose:
    DIF_DERIV computes the derivative of a polynomial in divided difference form.
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
    Input, int ND, the size of the input table.
    Input, double XD[ND], the abscissas for the divided
    difference table.
    Input, double YD[ND], the divided difference table.
    Output, int *NDP, the size of the output table, which is ND-1.
    Input, double XDP[NDP], the abscissas for the divided
    difference table for the derivative.
    Output, double YDP[NDP], the divided difference
    table for the derivative.
*/
{
	const dt2pitpdt2pit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xd = s_data->a1;
	ityp * yd = s_data->a2;
	dim_typ * ndp = s_data->a3;
	ityp * xdp = s_data->a4;
	ityp * ydp = s_data->a5;
	
    dim_typ i;
    /*
    Using a temporary copy of the difference table, shift the
    abscissas to zero.
    */
    ityp *xd_temp = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    ityp *yd_temp = ( ityp * ) malloc ( nd * sizeof ( ityp ) );

    for ( i = 0; i < nd; ++i )
    {
        xd_temp[i] = xd[i];
        yd_temp[i] = yd[i];
    }

    dif_shift_zero ( nd, xd_temp, yd_temp );
    /*
    Construct the derivative.
    */
    *ndp = nd - 1;

    for ( i = 0; i < *ndp; ++i )
    {
        xdp[i] = 0.00;
        ydp[i] = ( double ) ( i + 1 ) * yd_temp[i+1];
    }
    free ( xd_temp );
    free ( yd_temp );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_copy_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_COPY_NEW copies an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, ityp A1[N], the vector to be copied.
    Output, ityp r8VEC_COPY_NEW[N], the copy of A1.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	
    ityp *a2 = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for (dim_typ i = 0; i < n; ++i )
        a2[i] = a1[i];
    return a2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_ones_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_ONES_NEW creates a vector of 1's.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 March 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, ityp r8VEC_ONES_NEW[N], a vector of 1's.
*/
{
	const register dim_typ n = *(dim_typ *) data;

    ityp *a = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        a[i] = 1.00;
    return a;
}

#endif
