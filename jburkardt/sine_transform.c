#ifndef __DISABLEDEEP_SINETRANSFORM

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sine_transform_data ( void * data)
/******************************************************************************/
/*
  Purpose:
    SINE_TRANSFORM_DATA does a sine transform on a vector of data.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, integer N, the number of data points.
    Input, double D[N], the vector of data.
    Output, double SINE_TRANSFORM_DATA[N], the sine transform coefficients.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * d = s_data->a1;
	
    ityp angle;
    dim_typ i, j;
    ityp *s;

    s = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        s[i] = 0.00;
        for ( j = 0; j < n; ++j )
        {
            angle = M_PI * ( ityp ) ( ( i + 1 ) * ( j + 1 ) ) / ( ityp ) ( n + 1 );
            s[i] += sin ( angle ) * d[j];
        }
        s[i] *= sqrt ( 2.00 / ( ityp ) ( n + 1 ) );
    }
    return s;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sine_transform_function ( void * data)
/******************************************************************************/
/*
  Purpose:
    SINE_TRANSFORM_FUNCTION does a sine transform on functional data.
  Discussion:
    The interval [A,B] is divided into N+1 intervals using N+2 points,
    which are indexed by 0 through N+1.
    The original function F(X) is regarded as the sum of a linear function
    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
    which is 0 at A and B.
    The sine transform coefficients for F2 are then computed.
    To recover the interpolant of F(X), it is necessary to combine the
    linear part F1 with the sine transform interpolant:
      Interp(F)(X) = F1(X) + F2(X)
    This can be done by calling SINE_TRANSFORM_INTERPOLANT().
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data points.
    Input, double A, B, the interval endpoints.
    Input, double F ( double X ), a pointer to the function.
    Output, SINE_TRANSFORM_FUNCTION[N], the sine transform coefficients.
*/
{
	const dt2itfit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	ityp (* f)(ityp) = s_data->a3;
	
    ityp angle;
    ityp *f2;
    ityp fa;
    ityp fb;
    dim_typ i, j;
    ityp *s;
    ityp x;

    fa = f ( a );
    fb = f ( b );

    f2 = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        x = ( ( ityp ) ( n - i     ) * a+ ( ityp ) (     i + 1 ) * b )/ ( ityp ) ( n     + 1 );
        f2[i] = f ( x )                - ( ( b - x     ) * fa   + (     x - a ) * fb ) / ( b     - a );
    }

    s = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        s[i] = 0.00;

        for ( j = 0; j < n; +j )
        {
            angle = M_PI * ( ityp ) ( ( i + 1 ) * ( j + 1 ) ) / ( ityp ) ( n + 1 );
            s[i] = s[i] + sin ( angle ) * f2[j];
        }
        s[i] *= sqrt ( 2.00 / ( ityp ) ( n + 1 ) );
    }

    free ( f2 );

    return s;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sine_transform_interpolant ( void * data)
/******************************************************************************/
/*
  Purpose:
    SINE_TRANSFORM_INTERPOLANT evaluates the sine transform interpolant.
  Discussion:
    The interval [A,B] is divided into N+1 intervals using N+2 points,
    which are indexed by 0 through N+1.
    The original function F(X) is regarded as the sum of a linear function
    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
    which is 0 at A and B.
    The function F2 has been approximated using the sine transform,
    and the interpolant is then evaluated as:
      Interp(F)(X) = F1(X) + F2(X)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of terms in the approximation.
    Input, double A, B, the interval over which the approximant
    was defined.
    Input, double FA, FB, the function values at A and B.
    Input, double S[N], the approximant coefficients.
    Input, int NX, the number of evaluation points.
    Input, double X[NX], the evaluation points.
    Output, double SINE_TRANSFORM_INTERPOLANT[NX], the value of the interpolant.
*/
{
	const _2dt4it2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ nx = s_data->a1;
	ityp a = s_data->a2;
	ityp b = s_data->a3;
	ityp fa = s_data->a4;
	ityp fb = s_data->a5;
	ityp * s = s_data->a6;
	ityp * x = s_data->a7;
	
    ityp angle;
    ityp f1;
    ityp f2;
    dim_typ i, j;
    ityp *value;

    value = ( ityp * ) malloc ( nx * sizeof ( ityp ) );

    for ( i = 0; i < nx; ++i)
    {
        f1 = ( ( b - x[i]     ) * fa   + (     x[i] - a ) * fb ) / ( b           - a );
        f2 = 0.00;
        for ( j = 0; j < n; ++j)
        {
            angle = ( ityp ) ( j + 1 ) * ( x[i] - a ) * M_PI / ( b - a );
            f2 = f2 + s[j] * sin ( angle );
        }
        f2 *= sqrt ( 2.00 / ( ityp ) ( n + 1 ) );
        value[i] = f1 + f2;
    }
    return value;
}

#endif
