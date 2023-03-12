#ifndef __DISABLEDEEP_PWLAPPROX1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pwl_approx_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    PWL_APPROX_1D determines the control values for a PWL approximant.
  Discussion:
    The piecewise linear approximant is defined by NC control pairs
 (XC(I),YC(I)) and approximates ND data pairs (XD(I),YD(I)).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data points.
    ND must be at least 1.
    Input, double XD[ND], the data points.
    Input, double YD[ND], the data values.
    Input, int NC, the number of control points.
    NC must be at least 1.
    Input, double XC[NC], the control points.  Set these with a
    command like
      xc = r8vec_linspace_new ( nc, xmin, xmax );
    Output, double PWL_APPROX_1D[NC], the control values.
*/
{
	const _2dt3pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ nc = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	ityp * xc = s_data->a4;
	
    double *a;
    double *yc;
    /*
    Define the NDxNC linear system that determines the control values.
    */
    a = pwl_approx_1d_matrix ( nd, xd, yd, nc, xc );
    /*
    Solve the system.
    */
    yc = qr_solve ( nd, nc, a, yd );

    free ( a );

    return yc;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pwl_approx_1d_matrix ( void * data)
/******************************************************************************/
/*
  Purpose:
    PWL_APPROX_1D_MATRIX returns the matrix for the PWL approximant controls.
  Discussion:
    The value of the piecewise linear approximant, using control points XC
    and control values YC, evaluated at the point XD, can be represented by
      YD = A * YC
    where A is a matrix whose values depend on XC and XD.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data points.
    ND must be at least 1.
    Input, double XD[ND], the data points.
    Input, double YD[ND], the data values.
    Input, int NC, the number of control points.
    NC must be at least 1.
    Input, double XC[NC], the control points.
    Output, double PWL_APPROX_1D_MATRIX[ND*NC], the matrix.
*/
{
	const _2dt3pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ nc = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	ityp * xc = s_data->a4;
	
    ityp *a;
    dim_typ i, j, k;
    ityp t;

    a = ( ityp * ) malloc ( nd * nc * sizeof ( ityp ) );

    for ( j = 0; j < nc; ++j )
        for ( i = 0; i < nd; ++i)
            a[i+j*nd] = 0.00;

    for ( i = 0; i < nd; ++i )
    {
        k = nc - 2;
        for ( j = 1; j < nc - 1; ++j )
        {
            if ( xd[i] < xc[j] )
            {
                k = j - 1;
                break;
            }
        }
        t = ( xd[i] - xc[k] ) / ( xc[k+1] - xc[k] );
        a[i+k*nd]     = 1.00 - t;
        a[i+(k+1)*nd] =       t;
    }

    return a;
}

#endif
