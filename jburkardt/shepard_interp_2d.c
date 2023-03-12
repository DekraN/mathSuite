#ifndef __DISABLEDEEP_SHEPARDINTERP2D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _shepard_interp_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    SHEPARD_INTERP_2D evaluates a 2D Shepard interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2012
  Author:
    John Burkardt
  Reference:
    Donald Shepard,
    A two-dimensional interpolation function for irregularly spaced data,
    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
    ACM, pages 517-524, 1969.
  Parameters:
    Input, int ND, the number of data points.
    Input, double XD[ND], YD[ND], the data points.
    Input, double ZD[ND], the data values.
    Input, double P, the power.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], YI[NI], the interpolation points.
    Output, double SHEPARD_INTERP_2D[NI], the interpolated values.
*/
{
	const dt3pititdt2pit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xd = s_data->a1;
	ityp * yd = s_data->a2;
	ityp * zd = s_data->a3;
	const register ityp p = s_data->a4;
	const register dim_typ ni = s_data->a5;
	ityp * xi = s_data->a6;
	ityp * yi = s_data->a7;
	
    dim_typ i, j, k;
    ityp s;
    ityp *w;
    int z;
    ityp *zi;

    w = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    zi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( i = 0; i < ni; ++i )
    {
        if ( p == 0.00 )
            for ( j = 0; j < nd; ++j )
                w[j] = 1.00 / ( ityp ) ( nd );
        else
        {
            z = -1;
            for ( j = 0; j < nd; ++j )
            {
                w[j] = sqrt ( pow ( xi[i] - xd[j], 2 )+ pow ( yi[i] - yd[j], 2 ) );
                if ( w[j] == 0.00 )
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
                for ( j = 0; j < nd; ++j )
                    w[j] /= s;
            }
        }
        zi[i] = r8vec_dot_product ( nd, w, zd );
    }
    free ( w );

    return zi;
}

#endif
