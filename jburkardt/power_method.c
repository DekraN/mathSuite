#ifndef __DISABLEDEEP_POWERMETHOD

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _power_method ( void * data)
/******************************************************************************/
/*
  Purpose:
    POWER_METHOD applies several steps of the power method.
  Discussion:
    For a given NxN matrix A and an N vector Y, the power method produces
    a series of estimates for LAMBDA, the largest eigenvalue, and Y,
    the eigenvector corresponding to LAMBDA.
    The iteration repeats the following steps
      AY     = A * Y
      LAMBDA = || AY ||
      Y      = AY / LAMBDA
    If the matrix A has a single real eigenvalue of maximum modulus,
    then this iteration will generally produce a good estimate for that
    eigenvalue and its corresponding eigenvector.
    If there are multiple distinct eigenvalues of the same modulus,
    perhaps two values of opposite sign, or complex eigenvalues, then
    the situation is more complicated.
    Separate issues:
    * when estimating the value of LAMBDA, we use the Rayleigh quotient,
    LAMBDA = ( y' * A * y ) / ( y' * y ).  Since we normalize Y, the
    bottom of the fraction is 1.  Using this estimate allows us to
    easily capture the sign of LAMDBA.  Using the eucldean norm
    instead, for instance, would always give a positive value.
    * If the dominant eigenvalue is negative, then the iteration
    as given will produce eigenvector iterates that alternate in sign.

    * It is worth knowing whether the successive eigenvector estimates
    are tending to some value.  Since an eigenvector is really a direction,
    we need to normalize the vectors, and we need to somehow treat both
    a vector and its negative as holding the same information.  This
    means that the proper way of measuring the difference between two
    eigenvector estimates is to normalize them both, and then compute
    the cosine between them as y1'y2, followed by the sine, which is
    sqrt ( 1 - ( y1'y2)^2 ).  If this sine is small, the vectors y1 and y2
    are "close" in the sense of direction.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the matrix.
    Input/output, double Y[N], the estimate for the eigenvector.
    Input, int IT_MAX, the maximum number of iterations to take.
    1 <= IT_MAX.
    Input, double TOL, an error tolerance.
    Output, double *LAMBDA, the estimate for the eigenvalue.
    Output, int *IT_NUM, the number of iterations taken.
*/
{
	const dt2pitdtitpitpdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * y = s_data->a2;
	const register dim_typ it_max = s_data->a3;
	const register ityp tol = s_data->a4;
	ityp * lambda = s_data->a5;
	dim_typ * it_num = s_data->a6;
	
    ityp *ay;
    ityp cos_y1y2;
    dim_typ i;
    dim_typ it;
    ityp lambda_old;
    ityp norm;
    ityp sin_y1y2;
    ityp val_dif;
    ityp *y_old;

    ay = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    y_old = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    /*
    Force Y to be a vector of unit norm.
    */
    norm = r8vec_norm_l2 ( n, y );

    for ( i = 0; i < n; ++i)
        y[i] /= norm;
    it = 0;

    for ( i = 0; i < n; ++i )
        y_old[i] = y[i];

    r8mat_mv ( n, n, a, y, ay );
    *lambda = r8vec_dot_product ( n, y, ay );
    norm = r8vec_norm_l2 ( n, ay );
    for ( i = 0; i < n; ++i )
        y[i] = ay[i] / norm;

    if ( *lambda < 0.00 )
        for ( i = 0; i < n; ++i)
            y[i] *= -1;


    val_dif = 0.00;
    cos_y1y2 = r8vec_dot_product ( n, y, y_old );
    sin_y1y2 = sqrt ( ( 1.00 - cos_y1y2 ) * ( 1.00 + cos_y1y2 ) );

    for ( it = 1; it <= it_max; it++ )
    {
        lambda_old = *lambda;
        for ( i = 0; i < n; ++i )
            y_old[i] = y[i];

        r8mat_mv ( n, n, a, y, ay );
        *lambda = r8vec_dot_product ( n, y, ay );
        norm = r8vec_norm_l2 ( n, ay );
        for ( i = 0; i < n; ++i )
            y[i] = ay[i] / norm;
        if ( *lambda < 0.00 )
            for ( i = 0; i < n;++i )
            y[i] *= -1;

        val_dif = abs ( *lambda - lambda_old );
        cos_y1y2 = r8vec_dot_product ( n, y, y_old );
        sin_y1y2 = sqrt ( ( 1.00 - cos_y1y2 ) * ( 1.00 + cos_y1y2 ) );

        if ( val_dif <= tol )
            break;
    }

    for ( i = 0; i < n; ++i )
        y[i] = ay[i] / *lambda;

    *it_num = it;

    free ( ay );
    free ( y_old );

    return NULL;
}

#endif
