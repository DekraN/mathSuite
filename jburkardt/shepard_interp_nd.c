#ifndef __DISABLEDEEP_SHEPARDINTERPND

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _shepard_interp_nd ( void * data)
/******************************************************************************/
/*
  Purpose:
    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
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
    Input, int M, the spatial dimension.
    Input, int ND, the number of data points.
    Input, double XD[M*ND], the data points.
    Input, double ZD[ND], the data values.
    Input, double P, the power.
    Input, int NI, the number of interpolation points.
    Input, double XI[M*NI], the interpolation points.
    Output, double SHEPARD_INTERP_ND[NI], the interpolated values.
*/
{ 
	const _2dt2pititdtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ nd = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * zd = s_data->a3;
	const register ityp p = s_data->a4;
	const register dim_typ ni = s_data->a5;
	ityp * xi = s_data->a6;
	
	dim_typ i;
	dim_typ i2;
	dim_typ j;
	dim_typ k;
	ityp s;
	ityp t;
	ityp *w;
	dim_typ z;
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
				t = 0.00;
				for ( i2 = 0; i2 < m; ++i2 )
					t += pow ( xi[i2+i*m] - xd[i2+j*m], 2 );
				w[j] = sqrt ( t );
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
				s = 0.00;
				for ( j = 0; j < nd; ++j )
					s += w[j];
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
