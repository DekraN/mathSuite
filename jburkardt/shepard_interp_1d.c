#ifndef __DISABLEDEEP_SHEPARDINTERP1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _shepard_basis_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    SHEPARD_BASIS_1D evaluates a 1D Shepard basis function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2015
  Author:
    John Burkardt
  Reference:
    Donald Shepard,
    A two-dimensional interpolation function for irregularly spaced data,
    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
    ACM, pages 517-524, 1969.
  Parameters:
    Input, int ND, the number of data points.
    Input, double XD[ND], the data points.
    Input, double P, the power.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double SHEPARD_BASIS_1D[NI*ND], the basis function at the interpolation
    points.
*/
{
	const it2dt2pit * const s_data = data;
	
	const register ityp p = s_data->a0;
	const register dim_typ ni = s_data->a1;
	const register dim_typ nd = s_data->a2;
	ityp * xi = s_data->a3;
	ityp * xd = s_data->a4;
	
    ityp *bk;
    dim_typ i, j;
    ityp s;
    ityp *w;
    int z;

    w = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    bk = ( ityp * ) malloc ( ni * nd * sizeof (ityp ) );

    for ( i = 0; i < ni; ++i )
    {
        if ( p == 0.00 )
            for ( j = 0; j < nd; ++j)
                w[j] = 1.00 / ( ityp ) ( nd );
        else
        {
            z = -1;
            for ( j = 0; j < nd; ++j )
            {
                w[j] = abs ( xi[i] - xd[j] );
                if ( w[j] == 0.00 )
                {
                    z = j;
                    break;
                }
            }

            if ( z != -1 )
            {
                for ( j = 0; j < nd; ++j)
                    w[j] = 0.00;
                w[z] = 1.00;
            }
            else
            {
                for ( j = 0; j < nd; ++j )
                    w[j] = 1.00 / pow ( w[j], p );
                s = r8vec_sum ( nd, w );
                for ( j = 0; j < nd; ++j)
                    w[j] /= s;
            }
        }
        for ( j = 0; j < nd; ++j )
            bk[i+j*ni] = w[j];
    }

    free ( w );

    return bk;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _shepard_value_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    SHEPARD_VALUE_1D evaluates a 1D Shepard interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 October 2012
  Author:
    John Burkardt
  Reference:
    Donald Shepard,
    A two-dimensional interpolation function for irregularly spaced data,
    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
    ACM, pages 517-524, 1969.
  Parameters:
    Input, int ND, the number of data points.
    Input, double XD[ND], the data points.
    Input, double YD[ND], the data values.
    Input, double P, the power.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double SHEPARD_VALUE_1D[NI], the interpolated values.
*/
{
	const dt2pititdtpit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xd = s_data->a1;
	ityp * yd = s_data->a2;
	const register ityp p = s_data->a3;
	const register dim_typ ni = s_data->a4;
	ityp * xi = s_data->a5;
	
    dim_typ i, j, k;
    ityp s;
    ityp *w;
    ityp *yi;
    int z;

    w = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    yi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( i = 0; i < ni; ++i )
    {
        if ( p == 0.00 )
            for ( j = 0; j < nd; ++j )
            w[j] = 1.00 / ( ityp ) ( nd );
        else
        {
            z = -1;
            for ( j = 0; j < nd; ++j)
            {
                w[j] = abs ( xi[i] - xd[j] );
                if ( w[j] == 0.0 )
                {
                    z = j;
                    break;
                }
            }

            if ( z != -1 )
            {
                for ( j = 0; j < nd; ++j )
                    w[j] = 0.00;
                w[z] = 1.00;
            }
            else
            {
                for ( j = 0; j < nd; ++j )
                    w[j] = 1.00 / pow ( w[j], p );
                s = r8vec_sum ( nd, w );
                for ( j = 0; j < nd; ++j)
                    w[j] /= s;
            }
        }
        yi[i] = r8vec_dot_product ( nd, w, yd );
    }
    free ( w );

    return yi;
}

#endif
