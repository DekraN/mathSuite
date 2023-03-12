#ifndef __DISABLEDEEP_STOKES2DEXACT

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _resid_stokes1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_STOKES1 returns residuals of the exact Stokes solution #1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 January 2015
  Author:
    John Burkardt
  Reference:
    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.
  Parameters:
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Output, double UR[N], VR[N], PR[N], the residuals in the U,
    V and P equations.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * ur = s_data->a3;
	ityp * vr = s_data->a4;
	ityp * pr = s_data->a5;
	
    ityp *f;
    ityp *g;
    ityp *h;
    dim_typ i;
    ityp u;
    ityp ux;
    ityp uxx;
    ityp uy;
    ityp uyy;
    ityp v;
    ityp vx;
    ityp vxx;
    ityp vy;
    ityp vyy;
    /*
    Get the right hand sides.
    */
    f = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    g = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    h = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    rhs_stokes1 ( n, x, y, f, g, h );
    /*
    Form the functions and derivatives.
    */
    for ( i = 0; i < n; ++i )
    {
        u = - 2.00 * pow ( x[i], 2 ) * pow ( x[i] - 1.00, 2 ) * y[i] * ( y[i] - 1.00 ) * ( 2.00 * y[i] - 1.00 );
        ux = - 2.00 * ( 4.00 * pow ( x[i], 3 ) - 6.00 * pow ( x[i], 2 ) + 2.00 * x[i] ) * y[i] * ( y[i] - 1.00 ) * ( 2.00 * y[i] - 1.00 );
        uxx = - 2.00 * ( 12.00 * pow ( x[i], 2 ) - 12.00 * x[i] + 2.00 ) * ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        uy = - 2.00  * pow ( x[i], 2 ) * pow ( x[i] - 1.00, 2 )  * ( 6.00 * pow ( y[i], 2 ) - 3.00 * y[i] + 1.00 );
        uyy = - 2.00 * ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) ) * ( 12.00 * y[i] - 6.00 );
        v =   2.00 * x[i] * ( x[i] - 1.00 ) * ( 2.00 * x[i] - 1.0 ) * pow ( y[i], 2 ) * pow ( y[i] - 1.00, 2 );
        vx =   2.00 * ( 6.00 * pow ( x[i], 2 ) - 6.00 * x[i] + 1.00 ) * pow ( y[i], 2 ) * pow ( y[i] - 1.00, 2 );
        vxx =   2.00 * ( 12.00 * x[i] - 6.00 ) * pow ( y[i], 2 ) * pow ( y[i] - 1.00, 2 );
        vy =   2.00 * x[i] * ( x[i] - 1.0 ) * ( 2.00 * x[i] - 1.00 ) * ( 4.00 * pow ( y[i], 3 ) - 6.00 * pow ( y[i], 2 )  + 2.00 * y[i] );

        vyy =   2.00 * x[i] * ( x[i] - 1.0 ) * ( 2.00 * x[i] - 1.00 ) * ( 12.00 * pow ( y[i], 2 ) - 12.00 * y[i] + 2.00 );

        ur[i] = - ( uxx + uyy ) - f[i];
        vr[i] = - ( vxx + vyy ) - g[i];
        pr[i] = ux + vy - h[i];
    }
    /*
    Deallocate memory.
    */
    free ( f );
    free ( g );
    free ( h );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _resid_stokes2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_STOKES2 returns residuals of the exact Stokes solution #2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 January 2015
  Author:
    John Burkardt
  Reference:
    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.
  Parameters:
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Output, double UR[N], VR[N], PR[N], the residuals in the U,
    V and P equations.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * ur = s_data->a3;
	ityp * vr = s_data->a4;
	ityp * pr = s_data->a5;
	
    ityp *f;
    ityp *g;
    ityp *h;
    dim_typ i;
    ityp u;
    ityp ux;
    ityp uxx;
    ityp uy;
    ityp uyy;
    ityp v;
    ityp vx;
    ityp vxx;
    ityp vy;
    ityp vyy;
    /*
    Get the right hand sides.
    */
    f = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    g = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    h = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    rhs_stokes2 ( n, x, y, f, g, h );
    /*
    Form the functions and derivatives.
    */
    for ( i = 0; i < n; ++i )
    {
        u =   2.00* sin ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        ux =   4.00 * M_PI* cos ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        uxx = - 8.00 * pow ( M_PI, 2 )* sin ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        uy = - 4.00 * M_PI* sin ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );
        uyy = - 8.00 * pow ( M_PI, 2 )* sin ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        v = - 2.00* cos ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );
        vx =   4.00 * M_PI* sin ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );
        vxx =   8.00 * pow ( M_PI, 2 )* cos ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );
        vy = - 4.00 * M_PI* cos ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        vyy =   8.00 * pow ( M_PI, 2 )* cos ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );

        ur[i] = 2.00 * x[i] - ( uxx + uyy ) - f[i];
        vr[i] = 2.00 * y[i] - ( vxx + vyy ) - g[i];
        pr[i] = ux + vy - h[i];
    }
    /*
    Deallocate memory.
    */
    free ( f );
    free ( g );
    free ( h );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _resid_stokes3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_STOKES3 returns residuals of the exact Stokes solution #3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2015
  Author:
    John Burkardt
  Reference:
    Howard Elman, Alison Ramage, David Silvester,
    Finite Elements and Fast Iterative Solvers with
    Applications in Incompressible Fluid Dynamics,
    Oxford, 2005,
    ISBN: 978-0198528678,
    LC: QA911.E39.
  Parameters:
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Output, double UR[N], VR[N], PR[N], the residuals in the U,
    V and P equations.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * ur = s_data->a3;
	ityp * vr = s_data->a4;
	ityp * pr = s_data->a5;
	
    ityp *f;
    ityp *g;
    ityp *h;
    dim_typ i;
    ityp u;
    ityp ux;
    ityp uxx;
    ityp uy;
    ityp uyy;
    ityp v;
    ityp vx;
    ityp vxx;
    ityp vy;
    ityp vyy;
    /*
    Get the right hand sides.
    */
    f = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    g = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    h = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    rhs_stokes3 ( n, x, y, f, g, h );
    /*
    Form the functions and derivatives.
    */
    for ( i = 0; i < n; ++i )
    {
        u =   20.00 * x[i] * pow ( y[i], 3 );
        ux = 20.00 * pow ( y[i], 3 );
        uxx = 0.00;
        uy = 60.00 * x[i] * pow ( y[i], 2 );
        uyy = 120.00 * x[i] * y[i];

        v = 5.00 * ( pow ( x[i], 4 )  - pow ( y[i], 4 ) );
        vx = 20.00 * pow ( x[i], 3 );
        vxx = 60.00 * pow ( x[i], 2 );
        vy = - 20.00 * pow ( y[i], 3 );
        vyy = - 60.00 * pow ( y[i], 2 );

        ur[i] = 120.00 * x[i] * y[i] - ( uxx + uyy ) - f[i];
        vr[i] = 60.00 * pow ( x[i], 2 ) - 60.00 * pow ( y[i], 2 ) - ( vxx + vyy ) - g[i];
        pr[i] = ux + vy - h[i];
    }
    /*
    Deallocate memory.
    */
    free ( f );
    free ( g );
    free ( h );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _rhs_stokes1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_STOKES1 returns the right hand sides of the exact Stokes solution #1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 January 2015
  Author:
    John Burkardt
  Reference:
    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.
  Parameters:
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Output, double F[N], G[N], H[N], the right hand sides in the U,
    V and P equations.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * f = s_data->a3;
	ityp * g = s_data->a4;
	ityp * h = s_data->a5;
	
    for (dim_typ i = 0; i < n; ++i )
    {
        f[i] = + 2.00* ( 12.00 * pow ( x[i], 2 ) - 12.00 * x[i] + 2.00 )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] )+ 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 12.00 * y[i] - 6.00 );
        g[i] = - 2.00* ( 12.00 * x[i] - 6.00 )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) )- 2.00* ( 2.0 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( 12.00 * pow ( y[i], 2 ) - 12.00 * y[i] + 2.00 );
        h[i] = 0.00;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _rhs_stokes2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_STOKES2 returns the right hand sides of the exact Stokes solution #2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 January 2015
  Author:
    John Burkardt
  Reference:
    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.
  Parameters:
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Output, double F[N], G[N], H[N], the right hand sides in the U,
    V and P equations.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * f = s_data->a3;
	ityp * g = s_data->a4;
	ityp * h = s_data->a5;
	
    dim_typ i;
    ityp u;
    ityp ux;
    ityp uxx;
    ityp uy;
    ityp uyy;
    ityp v;
    ityp vx;
    ityp vxx;
    ityp vy;
    ityp vyy;

    for ( i = 0; i < n; ++i )
    {
        u =   2.00* sin ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );

        ux =   4.00 * M_PI* cos ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        uxx = - 8.00 * pow ( M_PI, 2 )* sin ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        uy = - 4.00 * M_PI* sin ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );
        uyy = - 8.00 * pow ( M_PI, 2 )* sin ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        v = - 2.00* cos ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );
        vx =   4.00 * M_PI* sin ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );
        vxx =   8.00 * pow ( M_PI, 2 )* cos ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );
        vy = - 4.00 * M_PI* cos ( M_2TPI * x[i] )* cos ( M_2TPI * y[i] );
        vyy =   8.00 * pow ( M_PI, 2 )* cos ( M_2TPI * x[i] )* sin ( M_2TPI * y[i] );

        f[i] = 2.00 * x[i] - ( uxx + uyy );
        g[i] = 2.00 * y[i] - ( vxx + vyy );
        h[i] = ux + vy;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rhs_stokes3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_STOKES3 returns the right hand sides of the exact Stokes solution #3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2015
  Author:
    John Burkardt
  Reference:
    Howard Elman, Alison Ramage, David Silvester,
    Finite Elements and Fast Iterative Solvers with
    Applications in Incompressible Fluid Dynamics,
    Oxford, 2005,
    ISBN: 978-0198528678,
    LC: QA911.E39.
  Parameters:
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Output, double F[N], G[N], H[N], the right hand sides in the U,
    V and P equations.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	ityp * f = s_data->a3;
	ityp * g = s_data->a4;
	ityp * h = s_data->a5;
	
    for (dim_typ i = 0; i < n; ++i )
        f[i] = g[i] = h[i] = 0.00;
    return NULL;
}

#endif
