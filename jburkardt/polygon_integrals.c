#ifndef __DISABLEDEEP_POLYGONINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _moment ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENT computes an unnormalized moment of a polygon.
  Discussion:
    Nu(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 October 2012
  Author:
    John Burkardt
  Reference:
    Carsten Steger,
    On the calculation of arbitrary moments of polygons,
    Technical Report FGBV-96-05,
    Forschungsgruppe Bildverstehen, Informatik IX,
    Technische Universitaet Muenchen, October 1996.
  Parameters:
    Input, int N, the number of vertices of the polygon.
    Input, double X[N], Y[N], the vertex coordinates.
    Input, int P, Q, the indices of the moment.
    Output, double MOMENT, the unnormalized moment Nu(P,Q).
*/
{
	static ityp result = MAX_VAL;
	
	const _3dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ p = s_data->a1;
	const register dim_typ q = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	
	
    dim_typ i, k, l;
    ityp nu_pq;
    ityp s_pq;
    ityp xi;
    ityp xj;
    ityp yi;
    ityp yj;

    nu_pq = 0.00;

    xj = x[n-1];
    yj = y[n-1];

    for ( i = 0; i < n; ++i )
    {
        xi = x[i];
        yi = y[i];

        s_pq = 0.00;
        for ( k = 0; k <= p; ++k )
            for ( l = 0; l <= q; ++l )
                s_pq = s_pq+ r8_choose ( k + l, l ) * r8_choose ( p + q - k - l, q - l )* pow ( xi, k ) * pow ( xj, p - k )* pow ( yi, l ) * pow ( yj, q - l );

        nu_pq += ( xj * yi - xi * yj ) * s_pq;

        xj = xi;
        yj = yi;
    }

	result = nu_pq / ( ityp ) ( p + q + 2 ) / ( ityp ) ( p + q + 1 )/ r8_choose ( p + q, p );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _moment_central ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENT_CENTRAL computes central moments of a polygon.
  Discussion:
    The central moment Mu(P,Q) is defined by
      Mu(P,Q) = Integral ( polygon ) (x-Alpha(1,0))^p (y-Alpha(0,1))^q dx dy
              / Area ( polygon )
    where
      Alpha(1,0) = Integral ( polygon ) x dx dy / Area ( polygon )
      Alpha(0,1) = Integral ( polygon ) y dx dy / Area ( polygon )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 October 2012
  Author:
    John Burkardt
  Reference:
    Carsten Steger,
    On the calculation of arbitrary moments of polygons,
    Technical Report FGBV-96-05,
    Forschungsgruppe Bildverstehen, Informatik IX,
    Technische Universitaet Muenchen, October 1996.
  Parameters:
    Input, int N, the number of vertices of the polygon.
    Input, double X[N], Y[N], the vertex coordinates.
    Input, int P, Q, the indices of the moment.
    Output, double MOMENT_CENTRAL, the unnormalized moment Mu(P,Q).
*/
{
	static ityp result = MAX_VAL;
	
	const _3dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ p = s_data->a1;
	const register dim_typ q = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	
	
    ityp alpha_01;
    ityp alpha_10;
    ityp alpha_ij;
    dim_typ i, j;
    ityp mu_pq;

    alpha_10 = moment_normalized ( n, x, y, 1, 0 );
    alpha_01 = moment_normalized ( n, x, y, 0, 1 );

    mu_pq = 0.00;

    for ( i = 0; i <= p; ++i )
        for ( j = 0; j <= q; ++j )
        {
            alpha_ij = moment_normalized ( n, x, y, i, j );
            mu_pq = mu_pq + r8_mop ( p + q - i - j )* r8_choose ( p, i ) * r8_choose ( q, j )* pow ( alpha_10, p - i ) * pow ( alpha_01, q - j ) * alpha_ij;
        }

	result = mu_pq;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _moment_normalized ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOMENT_NORMALIZED computes a normalized moment of a polygon.
  Discussion:
    Alpha(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy / Area ( polygon )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 October 2012
  Author:
    John Burkardt
  Reference:
    Carsten Steger,
    On the calculation of arbitrary moments of polygons,
    Technical Report FGBV-96-05,
    Forschungsgruppe Bildverstehen, Informatik IX,
    Technische Universitaet Muenchen, October 1996.
  Parameters:
    Input, int N, the number of vertices of the polygon.
    Input, double X[N], Y[N], the vertex coordinates.
    Input, int P, Q, the indices of the moment.
    Output, double MOMENT_NORMALIZED, the normalized moment Alpha(P,Q).
*/
{
	static ityp result = MAX_VAL;
	
	const _3dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ p = s_data->a1;
	const register dim_typ q = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	
	
	result = moment ( n, x, y, p, q ) / moment ( n, x, y, 0, 0 );
    return &result;
}

#endif
