#ifndef __DISABLEDEEP_PWLINTERP1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT void   * _pwl_interp_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    PWL_INTERP_1D evaluates the piecewise linear interpolant.
  Discussion:
    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
    linear function which interpolates the data (XD(I),YD(I)) for I = 1
    to ND.
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
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double PWL_INTERP_1D[NI], the interpolated values.
*/
{
	const _2dt3pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ ni = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	ityp * xi = s_data->a4;
	
    dim_typ i, k;
    ityp t;
    ityp *yi;

    yi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    if ( nd == 1 )
    {
        for ( i = 0; i < ni; ++i )
            yi[i] = yd[0];
        return yi;
    }

    for ( i = 0; i < ni; ++i )
    {
        if ( xi[i] <= xd[0] )
        {
            t = ( xi[i] - xd[0] ) / ( xd[1] - xd[0] );
            yi[i] = ( 1.00 - t ) * yd[0] + t * yd[1];
        }
        else if ( xd[nd-1] <= xi[i] )
        {
            t = ( xi[i] - xd[nd-2] ) / ( xd[nd-1] - xd[nd-2] );
            yi[i] = ( 1.00 - t ) * yd[nd-2] + t * yd[nd-1];
        }
        else
            for ( k = 1; k < nd; ++k )
                if ( xd[k-1] <= xi[i] && xi[i] <= xd[k] )
                    t = ( xi[i] - xd[k-1] ) / ( xd[k] - xd[k-1] );
                    yi[i] = ( 1.00 - t ) * yd[k-1] + t * yd[k];
                    break;
    }

    return yi;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pwl_basis_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    PWL_BASIS_1D evaluates a 1D piecewise linear basis function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2015
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data points.
    Input, double XD[ND], the data points.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double PW_BASIS_1D[NI*ND], the basis function at the
    interpolation points.
*/
{
	const _2dt2pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ ni = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * xi = s_data->a3;
	
    ityp *bk;
    dim_typ i, j;
    ityp t;

    bk = ( ityp * ) malloc ( ni * nd * sizeof ( ityp ) );

    for ( j = 0; j < nd; ++j )
        for ( i = 0; i < ni; ++i )
            bk[i+j*ni] = 0.00;

    if ( nd == 1 )
    {
        for ( j = 0; j < nd; ++j )
            for ( i = 0; i < ni; ++i )
            bk[i+j*ni] = 1.00;
        return bk;
    }

    for ( i = 0; i < ni; ++i )
    {
        for ( j = 0; j < nd; ++j )
        {
            if ( j == 0 && xi[i] <= xd[j] )
            {
                t = ( xi[i] - xd[j] ) / ( xd[j+1] - xd[j] );
                bk[i+j*ni] = 1.00 - t;
            }
            else if ( j == nd - 1 && xd[j] <= xi[i] )
            {
                t = ( xi[i] - xd[j-1] ) / ( xd[j] - xd[j-1] );
                bk[i+j*ni] = t;
            }
            else if ( xd[j-1] < xi[i] && xi[i] <= xd[j] )
            {
                t = ( xi[i] - xd[j-1] ) / ( xd[j] - xd[j-1] );
                bk[i+j*ni] = t;
            }
            else if ( xd[j] <= xi[i] && xi[i] < xd[j+1] )
            {
                t = ( xi[i] - xd[j] ) / ( xd[j+1] - xd[j] );
                bk[i+j*ni] = 1.00 - t;
            }
        }
    }

    return bk;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pwl_value_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    PWL_VALUE_1D evaluates the piecewise linear interpolant.
  Discussion:
    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
    linear function which interpolates the data (XD(I),YD(I)) for I = 1
    to ND.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data points.
    ND must be at least 1.
    Input, double XD[ND], the data points.
    Input, double YD[ND], the data values.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double PWL_VALUE_1D[NI], the interpolated values.
*/
{
	const _2dt3pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ ni = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	ityp * xi = s_data->a4;
	
    dim_typ i, k;
    ityp t;
    ityp *yi;

    yi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( i = 0; i < ni; ++i )
        yi[i] = 0.0;

    if ( nd == 1 )
    {
        for ( i = 0; i < ni; ++i)
            yi[i] = yd[0];
        return yi;
    }

    for ( i = 0; i < ni; ++i )
    {
        if ( xi[i] <= xd[0] )
        {
            t = ( xi[i] - xd[0] ) / ( xd[1] - xd[0] );
            yi[i] = ( 1.00 - t ) * yd[0] + t * yd[1];
        }
        else if ( xd[nd-1] <= xi[i] )
        {
            t = ( xi[i] - xd[nd-2] ) / ( xd[nd-1] - xd[nd-2] );
            yi[i] = ( 1.00 - t ) * yd[nd-2] + t * yd[nd-1];
        }
        else
            for ( k = 1; k < nd; ++k )
                if ( xd[k-1] <= xi[i] && xi[i] <= xd[k] )
                    t = ( xi[i] - xd[k-1] ) / ( xd[k] - xd[k-1] );
                    yi[i] = ( 1.00 - t ) * yd[k-1] + t * yd[k];
                    break;
    }

    return yi;
}

#endif
