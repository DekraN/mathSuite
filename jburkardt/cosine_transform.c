#ifndef __DISABLEDEEP_COSINETRANSFORM

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cosine_transform_data ( void * data)
/******************************************************************************/
/*
  Purpose:
    COSINE_TRANSFORM_DATA does a cosine transform on a vector of data.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2015
  Author:
    John Burkardt
  Parameters:
    Input, integer N, the number of data points.
    Input, double D[N], the vector of data.
    Output, double COSINE_TRANSFORM_DATA[N], the transform coefficients.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;	
	ityp * d = s_data->a1;
	
	
    ityp angle;
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dim_typ i, j;
    for ( i = 0; i < n; ++i )
    {
        c[i] = 0.00;
        for ( j = 0; j < n; ++j )
        {
            angle = M_PI * ( ityp ) ( i * ( (j<<1) + 1 ) ) / ( ityp ) ( n<<1 );
            c[i] += cos ( angle ) * d[j];
        }
        c[i] *= sqrt ( 2.00 / ( ityp ) ( n ) );
    }
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cosine_transform_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    COSINE_TRANSFORM_INVERSE does an inverse cosine transform.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 August 2015
  Author:
    John Burkardt
  Parameters:
    Input, integer N, the number of data points.
    Input, double C[N], the vector of transform coefficients
    Output, double COSINE_TRANSFORM_INVERSE[N], the original data.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;	
	ityp * c = s_data->a1;
	
    ityp angle;
    ityp *d = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dim_typ i, j;

    for ( i = 0; i < n; ++i )
    {
    d[i] = c[0] / 2.00;
        for ( j = 1; j < n; ++j )
        {
            angle = M_PI * ( ityp ) ( ( (i<<1) + 1 ) * j ) / ( ityp ) ( (n<<1) );
            d[i] += cos ( angle ) * c[j];
        }
        d[i] *= sqrt ( 2.00 / ( ityp ) ( n ) );
    }
    return d;
}

#endif
