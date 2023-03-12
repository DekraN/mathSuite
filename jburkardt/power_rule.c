#ifndef __DISABLEDEEP_POWERRULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _power_rule_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    POWER_RULE_SET sets up a power rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 February 2014
  Author:
    John Burkardt
  Parameters:
    Input, int POINT_NUM_1D, the order of the 1D rule.
    Input, double X_1D[POINT_NUM_1D], the points of the 1D rule.
    Input, double W_1D[POINT_NUM_1D], the weights of the
    1D rule.
    Input, double R_1D[2], the extreme points that define
    the range of the 1D region.
    Input, int DIM_NUM, the spatial dimension.
    Input, int POINT_NUM, the number of points in the rule.
    Output, double X[DIM_NUM*POINT_NUM], the points of the rule.
    Output, double W[POINT_NUM], the weights of the rule.
    Output, double R[DIM_NUM*2], the extreme points
    that define the range of the product rule region.
*/
{
	const _3dt6pit * const s_data = data;
	
	const register dim_typ point_num_1d = s_data->a0;
	const register dim_typ dim_num = s_data->a1;
	const register dim_typ point_num = s_data->a2;
	ityp * x_1d = s_data->a3;
	ityp * w_1d = s_data->a4;
	ityp * r_1d = s_data->a5;
	ityp * x = s_data->a6;
	ityp * w = s_data->a7;
	ityp * r = s_data->a8;
	
    dim_typ dim;
    dim_typ k;

    int *indx = ( int * ) malloc ( dim_num * sizeof ( int ) );
    k = 0;

    for ( ; ; )
    {
        tuple_next ( 0, point_num_1d-1, dim_num, &k, indx );

        if ( k == 0 )
            break;

        w[k-1] = 1.00;

        for ( dim = 0; dim < dim_num; ++dim)
        {
            w[k-1] *= w_1d[indx[dim]];
            x[dim+(k-1)*dim_num] = x_1d[indx[dim]];
        }
    }

    free ( indx );

    for ( dim = 0; dim < dim_num; ++dim )
    {
        r[dim+0*dim_num] = r_1d[0];
        r[dim+1*dim_num] = r_1d[1];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _power_rule_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    POWER_RULE_SIZE returns the size of a power rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 Febraury 2014
  Author:
    John Burkardt
  Parameters:
    Input, int POINT_NUM_1D, the number of points in the 1D rule.
    Input, int DIM_NUM, the spatial dimension.
    Output, int POWER_RULE_SIZE, the number of points in the rule.
*/
{
	static int result = INT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ point_num_1d = a_data[0];
	const register dim_typ dim_num = a_data[1];
	
	result = powi ( point_num_1d, dim_num );
    return &result;
}

#endif
