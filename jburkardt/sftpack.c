#ifndef __DISABLEDEEP_SFTPACK

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT void   * _c8vec_sftb ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8VEC_SFTB computes a "slow" backward Fourier transform of complex data.
  Discussion:
    SFTF and SFTB are inverses of each other.  If we begin with data
    X and apply SFTF to get Y, and then apply SFTB to Y,
    we should get back the original X.
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines that use fast Fourier techniques.
    For 0 <= I <= N - 1
      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 M_PI i I J / N )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, double complex Y[N], the Fourier coefficients.
    Output, double complex C8VEC_SFTB, the data.
*/
{
	const dtpcx * const s_data = data;
	const register dim_typ n = s_data->a0;
	double complex * y = s_data->a1;
	
    dim_typ i, j;
    ityp theta;
    double complex *x = ( double complex * ) malloc ( n * sizeof ( double complex ) );

    for ( i = 0; i < n; ++i )
    {
        __real__ x[i] = __imag__ x[i]= 0.00;
        for ( j = 0; j < n; ++j )
        {
            theta = - M_2TPI * ( ityp ) ( i * j ) / ( ityp ) ( n );
            __real__ x[i] += creal(y[j]) * cos ( theta )- cimag(y[j]) * sin ( theta );
            __imag__ x[i] += cimag(y[j]) * cos ( theta )+ creal(y[j]) * sin ( theta );
        }
        __real__ x[i] /= ( ityp ) ( n );
        __imag__ x[i] /= ( ityp ) ( n );
    }
    return x;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT void   * _c8vec_sftf ( void * data)
/******************************************************************************/
/*
  Purpose:
    C8VEC_SFTF computes a "slow" forward Fourier transform of complex data.
  Discussion:
    SFTF and SFTB are inverses of each other.  If we begin with data
    X and apply SFTF to get Y, and then apply SFTB to Y,
    we should get back the original X.
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines that use fast Fourier techniques.
    For 0 <= I <= N - 1
      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 M_PI i I J / N )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 June 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, double complex X[N], the data to be transformed.
    Output, double complex C8VEC_SFTF[N], the Fourier coefficients.
*/
{
	const dtpcx * const s_data = data;
	const register dim_typ n = s_data->a0;
	double complex * x = s_data->a1;
	
    dim_typ i, j;
    ityp theta;
    double complex *y = ( double complex * ) malloc ( n * sizeof ( double complex ) );

    for ( i = 0; i < n; ++i )
    {
        __real__ y[i] = __imag__ y[i] = 0.00;
        for ( j = 0; j < n; ++j )
        {
    		theta = - M_2TPI * ( ityp ) ( i * j ) / ( ityp ) ( n );
            __real__ y[i] += creal(x[j]) * cos ( theta )
                            + cimag(x[j]) * sin ( theta );
            __imag__ y[i] += cimag(x[j]) * cos ( theta )
                            - creal(x[j]) * sin ( theta );
        }
    }
    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sct ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SCT computes a "slow" cosine transform on real data.
  Discussion:
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines that use fast Fourier techniques.
      Y(1) = Sum ( 1 <= J <= N ) X(J)
      For 2 <= I <= N-1:
        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J)
          * cos ( M_PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
    Applying the routine twice in succession should yield the original data,
    multiplied by 2 * ( N + 1 ).  This is a good check for correctness
    and accuracy.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, double X[N], the data sequence.
    Output, double SCT[N], the transformed data.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp angle;
    dim_typ i, j;
    ityp *y = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        y[i] = x[0] / 2.00;

        for ( j = 1; j < n - 1; ++j)
        {
            angle = M_PI * ( ityp ) ( ( i * j ) % ( ( n - 1 ) <<1 ) )/ ( ityp ) ( n - 1 );
            y[i] = y[i] + x[j] * cos ( angle );
        }

        j = n - 1;

        angle = M_PI * ( ityp ) ( ( i * j ) % ( ( n - 1 ) <<1 ) )/ ( ityp ) ( n - 1 );

        y[i] += x[j] * cos ( angle ) / 2.00;
    }

    for ( i = 0; i < n; ++i )
        y[i] = 2.00 * y[i] * sqrt ( ( ityp ) ( n ) / ( ityp ) ( n - 1 ) );


    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sftb ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SFTB computes a "slow" backward Fourier transform of real data.
  Discussion:
    SFTB and SFTF are inverses of each other.  If we begin with data
    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
    resulting R vector, we should get back the original AZERO, A and B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, double AZERO, the constant Fourier coefficient.
    Input, double A[N/2], B[N/2], the Fourier coefficients.
    Output, double SFTB[N], the reconstructed data sequence.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp azero = s_data->a1;
	ityp * a = s_data->a2;
	ityp * b = s_data->a3;
	
    dim_typ i, k;
    ityp *r;
    ityp theta;
    r = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        r[i] = azero;
        for ( k = 0; k < ( n / 2 ); ++k )
        {
            theta = ( ityp ) ( ( k + 1 ) * i <<1 ) * M_PI / ( ityp ) ( n );
            r[i] += a[k] * cos ( theta ) + b[k] * sin ( theta );
        }
    }

    return r;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sftf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SFTF computes a "slow" forward Fourier transform of real data.
  Discussion:
    SFTF and SFTB are inverses of each other.  If we begin with data
    R and apply SFTB to it, and then apply SFTB to the resulting AZERO,
    A, and B, we should get back the original R.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, double R[N], the data to be transformed.
    Output, double *AZERO, = sum ( 1 <= I <= N ) R(I) / N.
    Output, double A[N/2], B[N/2], the Fourier coefficients.
*/
{
	const dt4pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * r = s_data->a1;
	ityp * azero = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	
    dim_typ i, j;
    ityp theta;

    *azero = 0.00;
    for ( i = 0; i < n; ++i )
        *azero += r[i];
    *azero /= ( ityp ) ( n );

    for ( i = 0; i < ( n / 2 ); ++i )
    {
        a[i] = b[i] = 0.00;

        for ( j = 0; j < n; ++j )
        {
            theta = ( ityp ) ( ( i + 1 ) * j << 1 ) * M_PI / ( ityp ) ( n );
            a[i] += r[j] * cos ( theta );
            b[i] +=  r[j] * sin ( theta );
        }

        a[i] /= ( ityp ) ( n );
        b[i] /= ( ityp ) ( n );

        if ( i != ( n / 2 - 1 ) )
        {
            a[i] *= 2.00;
            b[i] *= 2.00;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sht ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SHT computes a "slow" Hartley transform of real data.
  Discussion:
    The discrete Hartley transform B of a set of data A is
      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*M_PI*I*J/N)
    Here, the data and coefficients are indexed from 0 to N-1.
    With the above normalization factor of 1/sqrt(N), the Hartley
    transform is its own inverse.
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Reference:
    Ralph Hartley,
    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
    Proceedings of the Institute of Radio Engineers,
    Volume 30, pages 144-150, 1942.
  Parameters:
    Input, int N, the number of data values.
    Input, double A[N], the data to be transformed.
    Output, double SHT[N], the transformed data.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    double *b;
    dim_typ i, j;
    ityp theta;

    b = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        b[i] = 0.00;

    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j )
        {
            theta = M_2TPI * ( ityp ) ( ( i * j ) % n ) / ( ityp ) ( n );
            b[i] += a[j] * ( cos ( theta ) + sin ( theta ) );
        }

    for ( i = 0; i < n; ++i )
        b[i] /= sqrt ( ( ityp ) ( n  ) );

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sqctb ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SQCTB computes a "slow" quarter cosine transform backward on real data.
  Discussion:
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines that use fast Fourier techniques.
    For 0 <= I <= N-1,
      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( M_PI * J * (I+1/2) / N )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Reference:
    William Briggs, Van Emden Henson,
    The Discrete Fourier Transform,
    SIAM, 1995,
    LC: QA403.5 B75
  Parameters:
    Input, int N, the number of data values.
    Input, double X[N], the data sequence.
    Output, double SQCTB[N], the transformed data.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, j;
    ityp theta;
    ityp *y;

    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        y[i] = x[0];

    for ( i = 0; i < n; ++i )
        for ( j = 1; j < n; ++j )
        {
            theta = 0.50 * M_PI * ( ityp ) ( j * ( (i<<1) + 1 ) ) / ( ityp ) ( n );
            y[i] += 2.00 * x[j] * cos ( theta  );
        }

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sqctf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SQCTF computes a "slow" quarter cosine transform forward on real data.
  Discussion:
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines that use fast Fourier techniques.
    For 0 <= I <= N-1,
      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( M_PI * I * (J+1/2) / N )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Reference:
    William Briggs, Van Emden Henson,
    The Discrete Fourier Transform,
    SIAM, 1995,
    QA403.5 B75
  Parameters:
    Input, int N, the number of data values.
    Input, double X[N], the data sequence.
    Output, double SQCTF[N], the transformed data.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, j;
    ityp theta;
    ityp *y;

    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        y[i] = 0.00;

    for ( i = 0; i < n; ++i)
        for ( j = 0; j < n; ++j )
        {
            theta = 0.50 * M_PI * ( ityp ) ( i * ( (j<<1) + 1 ) ) / ( ityp ) ( n );
            y[i] += x[j] * cos ( theta  );
        }

    for ( i = 0; i < n; ++i )
        y[i] /= ( ityp ) ( n );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sqstb ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SQSTB computes a "slow" quarter sine transform backward on real data.
  Discussion:
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines that use fast Fourier techniques.
    For 0 <= I <= N-1,
      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( M_PI * J * (I+1/2) / N )
             - X(N) * cos ( M_PI * I )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Reference:
    William Briggs, Van Emden Henson,
    The Discrete Fourier Transform,
    SIAM, 1995,
    QA403.5 B75
  Parameters:
    Input, int N, the number of data values.
    Input, double X[N], the data sequence.
    Output, double SQSTB[N], the transformed data.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, j;
    ityp theta;
    ityp *y;

    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i)
        y[i] = 0.00;

    for ( i = 0; i < n; ++i )
    {
        for ( j = 0; j < n - 1; ++j )
        {
            theta = 0.50 * M_PI * ( ityp ) ( ( j + 1 ) * ( (i<<1) + 1 ) )/ ( ityp ) ( n );
            y[i] -= 2.00 * x[j] * sin ( theta  );
        }
        theta = M_PI * ( ityp ) ( i );
        y[i] -= x[n-1] * cos ( theta );
    }
    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void    * _r8vec_sqstf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SQSTF computes a "slow" quarter sine transform forward on real data.
  Discussion:
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines that use fast Fourier techniques.
    For 1 <= I <= N,
      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( M_PI * I * (J+1/2) / N )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Reference:
    William Briggs, Van Emden Henson,
    The Discrete Fourier Transform,
    SIAM, 1995,
    QA403.5 B75
  Parameters:
    Input, int N, the number of data values.
    Input, double X[N], the data sequence.

    Outut, double SQSTF{N], the transformed data.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, j;
    ityp theta;
    ityp *y = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        y[i] = 0.00;

    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j )
        {
            theta = 0.50 * M_PI * ( ityp ) ( ( i + 1 ) * ( (j<<1) + 1 ) )/ ( ityp ) ( n );
            y[i] += x[j] * sin ( theta );
        }

    for ( i = 0; i < n; ++i )
        y[i] /= -( ityp ) ( n );

    return y;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sst ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_SST computes a "slow" sine transform on real data.
  Discussion:
    This routine is provided for illustration and testing.  It is inefficient
    relative to optimized routines that use fast Fourier techniques.
    For 1 <= I <= N,
      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( M_PI * I * J / ( N + 1 ) )
    Applying the routine twice in succession should yield the original data,
    multiplied by N / 2.  This is a good check for correctness and accuracy.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, double X[N], the data sequence.
    Output, double SST[N], the transformed data.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, j;
    ityp *theta;
    ityp *y;

    theta = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        theta[i] = M_PI * ( ityp ) ( i + 1 ) / ( ityp ) ( n + 1 );


    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for ( i = 0; i < n; ++i )
        y[i] = 0.00;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            y[i] += 2.00 * x[j] * sin ( ( ityp ) ( j + 1 ) * theta[i] );

    free ( theta );

    return y;
}

#endif
