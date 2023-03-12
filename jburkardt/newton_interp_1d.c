#ifndef __DISABLEDEEP_NEWTONINTERP1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _newton_coef_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    NEWTON_COEF_1D computes coefficients of a Newton 1D interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 July 2015
  Author:
    John Burkardt
  Reference:
    Carl deBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.
  Parameters:
    Input, int ND, the number of data points.
    Input, double XD[ND], the X values at which data was taken.
    These values must be distinct.
    Input, double YD[ND], the corresponding Y values.
    Output, double NEWTON_COEF_1D[ND], the divided difference coefficients.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xd = s_data->a1;
	ityp * yd = s_data->a2;
	
    dim_typ i, j;
    /*
    Copy the data values.
    */
    ityp *cd = ( double * ) malloc ( nd * sizeof ( double ) );

    for ( i = 0; i < nd; ++i )
        cd[i] = yd[i];
    /*
    Make sure the abscissas are distinct.
    */
    for ( i = 0; i < nd; ++i )
        for ( j = i + 1; j < nd; ++j )
            if ( xd[i] - xd[j] == 0.00 )
                return NULL;
    /*
    Compute the divided differences.
    */
    for ( i = 1; i <= nd - 1; ++i )
        for ( j = nd - 1; i <= j; --j)
                cd[j] = ( cd[j] - cd[j-1] ) / ( xd[j] - xd[j-i] );
    return cd;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _newton_value_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    NEWTON_VALUE_1D evaluates a Newton 1D interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 July 2015
  Author:
    John Burkardt
  Reference:
    Carl deBoor,
    A Practical Guide to Splines,
    Springer, 2001,
    ISBN: 0387953663,
    LC: QA1.A647.v27.
  Parameters:
    Input, int ND, the order of the difference table.
    Input, double XD[ND], the X values of the difference table.
    Input, double CD[ND], the divided differences.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], the interpolation points.
    Output, double NEWTON_VALUE_1D[NI], the interpolation values.
*/
{
	const _2dt3pit * const s_data = data;
	
	const register dim_typ ni = s_data->a0;
	const register dim_typ nd = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * cd = s_data->a3;
	ityp * xi = s_data->a4;
	
    dim_typ i, j;
    ityp *yi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( j = 0; j < ni; ++j )
    {
        yi[j] = cd[nd-1];
        for ( i = 2; i <= nd; ++i )
            yi[j] = cd[nd-i] + ( xi[j] - xd[nd-i] ) * yi[j];
    }
    return yi;
}

#endif
