#ifndef __DISABLEDEEP_NEARESTINTERP1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _nearest_interp_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    NEAREST_INTERP_1D evaluates the nearest neighbor interpolant.
  Discussion:
    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
    constant function which interpolates the data (XD(I),YD(I)) for I = 1
    to ND.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data points.
    ND must be at least 1.
    Input, double XD[ND], the data points.
    Input, double YD[ND], the data values.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double NEAREST_INTERP_1D[NI], the interpolated values.
*/
{
	const _2dt3pit * const s_data = data;
	
	const register dim_typ nd = s_data->a0;
	const register dim_typ ni = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	ityp * xi = s_data->a4;
	
    ityp d;
    ityp d2;
    dim_typ i, j, k;
    ityp *yi;

    yi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( i = 0; i < ni; ++i )
    {
        k = 0;
        d = abs ( xi[i] - xd[k] );
        for ( j = 1; j < nd; ++j )
        {
            d2 = abs ( xi[i] - xd[j] );
            if ( d2 < d )
            {
                k = j;
                d = d2;
            }
        }
        yi[i] = yd[k];
    }

    return yi;
}

#endif
