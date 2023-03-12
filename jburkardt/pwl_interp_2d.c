#ifndef __DISABLEDEEP_PWLINTERP2D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_bracket5 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
  Discussion:
    We assume XD is sorted.
    If XI is contained in the interval [XD(1),XD(N)], then the returned 
    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
    This code implements a version of binary search which is perhaps more
    understandable than the usual ones.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ND, the number of data values.
    Input, double XD[N], the sorted data.
    Input, double XD, the query value.
    Output, int R8VEC_BRACKET5, the bracket information.
*/
{
	static int result = INT_MAX;
	
	const dtpitit * const s_data = data;
	const register dim_typ nd = s_data->a0;
	ityp * xd = s_data->a1;
	const register ityp xi = s_data->a2;
	
	int b;
	int l;
	int m;
	int r;
	
	if ( xi < xd[0] || xd[nd-1] < xi )
		b = -1;
	else
	{
		l = 0;
		r = nd - 1;
		
		while ( l + 1 < r )
		{
			m = ( l + r ) / 2;
			if ( xi < xd[m] )
				r = m;
			else
				l = m;
		}
		b = l;
	}
	
	result = b;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pwl_interp_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 October 2012
  Author:
    John Burkardt
  Parameters:
    Input, int NXD, NYD, the number of X and Y data values.
    Input, double XD[NXD], YD[NYD], the sorted X and Y data.
    Input, double ZD[NXD*NYD}, the Z data.
    Input, int NI, the number of interpolation points.
    Input, double XI[NI], YI[NI], the coordinates of the
    interpolation points.
    Output, double PWL_INTERP_2D[NI], the value of the interpolant.
*/
{
	const _2dt3pitdt2pit * const s_data = data;
	const register dim_typ nxd = s_data->a0;
	const register dim_typ nyd = s_data->a1;
	ityp * xd = s_data->a2;
	ityp * yd = s_data->a3;
	ityp * zd = s_data->a4;
	const register dim_typ ni = s_data->a5;
	ityp * xi = s_data->a6;
	ityp * yi = s_data->a7;
	
    ityp alpha;
    ityp beta;
    ityp det;
    ityp dxa;
    ityp dxb;
    ityp dxi;
    ityp dya;
    ityp dyb;
    ityp dyi;
    ityp gamma;
    dim_typ i, j, k;
    ityp *zi;

    zi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );

    for ( k = 0; k < ni; ++k )
    {
        i = r8vec_bracket5 ( nxd, xd, xi[k] );
        if ( i == -1 )
        {
            zi[k] = r8_huge;
            continue;
        }

        j = r8vec_bracket5 ( nyd, yd, yi[k] );
        if ( j == -1 )
        {
            zi[k] = r8_huge;
            continue;
        }

        if ( yi[k] < yd[j+1] + ( yd[j] - yd[j+1] ) * ( xi[i] - xd[i] ) / ( xd[i+1] - xd[i] ) )
        {
            dxa = xd[i+1] - xd[i];
            dya = yd[j]   - yd[j];

            dxb = xd[i]   - xd[i];
            dyb = yd[j+1] - yd[j];

            dxi = xi[k]   - xd[i];
            dyi = yi[k]   - yd[j];

            det = dxa * dyb - dya * dxb;

            alpha = ( dxi * dyb - dyi * dxb ) / det;
            beta = ( dxa * dyi - dya * dxi ) / det;
            gamma = 1.00 - alpha - beta;

            zi[k] = alpha * zd[i+1+j*nxd] + beta * zd[i+(j+1)*nxd] + gamma * zd[i+j*nxd];
        }
        else
        {
            dxa = xd[i]   - xd[i+1];
            dya = yd[j+1] - yd[j+1];

            dxb = xd[i+1] - xd[i+1];
            dyb = yd[j]   - yd[j+1];

            dxi = xi[k]   - xd[i+1];
            dyi = yi[k]   - yd[j+1];

            det = dxa * dyb - dya * dxb;

            alpha = ( dxi * dyb - dyi * dxb ) / det;
            beta = ( dxa * dyi - dya * dxi ) / det;
            gamma = 1.00 - alpha - beta;

            zi[k] = alpha * zd[i+(j+1)*nxd] + beta * zd[i+1+j*nxd] + gamma * zd[i+1+(j+1)*nxd];
        }
    }

    return zi;
}

#endif
