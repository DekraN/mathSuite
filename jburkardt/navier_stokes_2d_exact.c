#ifndef __DISABLEDEEP_NAVIERSTOKES2DEXACT

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _resid_lucas ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_LUCAS returns Lucas Bystricky residuals.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 March 2015
  Author:
    John Burkardt
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the density.
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Input, double T, the time coordinate or coordinates.
    Output, double UR[N], VR[N], PR[N], the residuals in the U,
    V and P equations.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * ur = s_data->a6;
	ityp * vr = s_data->a7;
	ityp * pr = s_data->a8;
	
    ityp dpdx;
    ityp dpdy;
    ityp dudt;
    ityp dudx;
    ityp dudxx;
    ityp dudy;
    ityp dudyy;
    ityp dvdt;
    ityp dvdx;
    ityp dvdxx;
    ityp dvdy;
    ityp dvdyy;
    ityp *f;
    ityp *g;
    ityp *h;
    dim_typ i;
    ityp p;
    ityp u;
    ityp v;
    /*
    Get the right hand sides.
    */
    f = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    g = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    h = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    rhs_lucas ( nu, rho, n, x, y, t, f, g, h );
    /*
    Form the functions and derivatives of the left hand side.
    */
    for ( i = 0; i < n; ++i )
    {
        u = - cos ( M_PI * x[i] ) / M_PI;
        dudx = sin ( M_PI * x[i] );
        dudxx = M_PI * cos ( M_PI * x[i] );
        dudt = dudy = dudyy = 0.00;
        v = - y[i] * sin ( M_PI * x[i] );
        dvdx = - M_PI * y[i] * cos ( M_PI * x[i] );
        dvdxx = + M_PI * M_PI * y[i] * sin ( M_PI * x[i] );
        dvdy = - sin ( M_PI * x[i] );

        dvdt = dvdyy = p = dpdx = dpdy = 0.00;
        /*
        Evaluate the residuals.
        */
        ur[i] = dudt + u * dudx + v * dudy+ ( 1.00 / rho ) * dpdx - nu * ( dudxx + dudyy ) - f[i];
        vr[i] = dvdt + u * dvdx + v * dvdy+ ( 1.00 / rho ) * dpdy - nu * ( dvdxx + dvdyy ) - g[i];
        pr[i] = dudx + dvdy - h[i];
    }
    /*
    Free memory.
    */
    free ( f );
    free ( g );
    free ( h );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _resid_poiseuille ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_POISEUILLE returns Poiseuille flow residuals.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 July 2015
  Author:
    John Burkardt
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the density.
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Input, double T, the time coordinate or coordinates.
    Output, double UR[N], VR[N], PR[N], the residuals in the U,
    V and P equations.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * ur = s_data->a6;
	ityp * vr = s_data->a7;
	ityp * pr = s_data->a8;
	
    ityp dpdx;
    ityp dpdy;
    ityp dudt;
    ityp dudx;
    ityp dudxx;
    ityp dudy;
    ityp dudyy;
    ityp dvdt;
    ityp dvdx;
    ityp dvdxx;
    ityp dvdy;
    ityp dvdyy;
    ityp *f;
    ityp *g;
    ityp *h;
    dim_typ i;
    ityp p;
    ityp u;
    ityp v;
    /*
    Get the right hand sides.
    */
    f = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    g = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    h = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    rhs_poiseuille ( nu, rho, n, x, y, t, f, g, h );
    /*
    Form the functions and derivatives of the left hand side.
    */
    for ( i = 0; i < n; ++i )
    {
        u = 1.00 - y[i] * y[i];
        dudxx = 0.00;
        dudy = - 2.00 * y[i];
        dudyy = - 2.00;
        dudt = dudx = 0.00;
        v = dvdt = dvdx = dvdxx = dvdy = dvdyy = dpdy = 0.00;

        p = - 2.00 * rho * nu * x[i];
        dpdx = - 2.00 * rho * nu;
        /*
        Evaluate the residuals.
        */
        ur[i] = dudt + u * dudx + v * dudy+ ( 1.00 / rho ) * dpdx - nu * ( dudxx + dudyy ) - f[i];
        vr[i] = dvdt + u * dvdx + v * dvdy+ ( 1.00 / rho ) * dpdy - nu * ( dvdxx + dvdyy ) - g[i];
        pr[i] = dudx + dvdy - h[i];
    }
    /*
    Free memory.
    */
    free ( f );
    free ( g );
    free ( h );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _resid_spiral ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_SPIRAL returns Spiral residuals.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2015
  Author:
    John Burkardt
  Reference:
    Maxim Olshanskii, Leo Rebholz,
    Application of barycenter refined meshes in linear elasticity
    and incompressible fluid dynamics,
    ETNA: Electronic Transactions in Numerical Analysis,
    Volume 38, pages 258-274, 2011.
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the fluid density.
    Input, int N, the number of nodes.
    Input, double X[N], Y[N], the coordinates of nodes.
    Input, double T, the current time.
    Output, double UR[N], VR[N], PR[N], the residuals sides.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * ur = s_data->a6;
	ityp * vr = s_data->a7;
	ityp * pr = s_data->a8;
	
    ityp *dpdx;
    ityp *dpdy;
    ityp *dudt;
    ityp *dudx;
    ityp *dudxx;
    ityp *dudy;
    ityp *dudyy;
    ityp *dvdt;
    ityp *dvdx;
    ityp *dvdxx;
    ityp *dvdy;
    ityp *dvdyy;
    ityp *f;
    ityp *g;
    ityp *h;
    dim_typ i;
    ityp *p;
    ityp *u;
    ityp *v;

    dpdx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dpdy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudt = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudxx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudyy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdt = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdxx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdyy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    f = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    g = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    h = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    p = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    u = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    v = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    /*
    Get the right hand side functions.
    */
    rhs_spiral ( nu, rho, n, x, y, t, f, g, h );
    /*
    Form the functions and derivatives for the left hand side.
    */
    for ( i = 0; i < n; ++i )
    {
        u[i] = ( 1.00 + nu * t ) * 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        dudt[i] = nu * 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        dudx[i] = ( 1.00 + nu * t ) * 2.00* ( 4.00 * pow ( x[i], 3 ) - 6.00 * pow ( x[i], 2 ) + 2.00 * x[i] )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        dudxx[i] = ( 1.00 + nu * t ) * 2.00* ( 12.00 * pow ( x[i], 2 ) - 12.00 * x[i] + 2.00 )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        dudy[i] = ( 1.00 + nu * t ) * 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 6.00 * pow ( y[i], 2 ) - 6.00 * y[i] + 1.00 );
        dudyy[i] = ( 1.00 + nu * t ) * 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 12.00 * y[i] - 6.00 );
        v[i] = - ( 1.00 + nu * t ) * 2.00* ( 2.00 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );
        dvdt[i] = - nu * 2.00* ( 2.00 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );
        dvdx[i] = - ( 1.0 + nu * t ) * 2.00* ( 6.00 * pow ( x[i], 2 ) - 6.00 * x[i] + 1.00 )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );
        dvdxx[i] = - ( 1.0 + nu * t ) * 2.00* ( 12.00 * x[i] - 6.00 )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );
        dvdy[i] = - ( 1.00 + nu * t ) * 2.00* ( 2.00 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( 4.00 * pow ( y[i], 3 ) - 6.00 * pow ( y[i], 2 ) + 2.00 * y[i] );
        dvdyy[i] = - ( 1.00 + nu * t ) * 2.00* ( 2.00 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( 12.00 * pow ( y[i], 2 ) - 12.00 * y[i] + 2.00 );

        p[i] = rho * y[i];
        dpdx[i] = 0.00;
        dpdy[i] = rho;
        /*
        Evaluate the residuals.
        */
        ur[i] = dudt[i] - nu * ( dudxx[i] + dudyy[i] )+ u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho - f[i];
        vr[i] = dvdt[i] - nu * ( dvdxx[i] + dvdyy[i] )+ u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho - g[i];
        pr[i] = dudx[i] + dvdy[i] - h[i];
    }

    free ( dpdx );
    free ( dpdy );
    free ( dudt );
    free ( dudx );
    free ( dudxx );
    free ( dudy );
    free ( dudyy );
    free ( dvdt );
    free ( dvdx );
    free ( dvdxx );
    free ( dvdy );
    free ( dvdyy );
    free ( f );
    free ( g );
    free ( h );
    free ( p );
    free ( u );
    free ( v );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _resid_taylor ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_TAYLOR returns Taylor residuals.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2015
  Author:
    John Burkardt
  Reference:
    Geoffrey Taylor,
    On the decay of vortices in a viscous fluid,
    Philosophical Magazine,
    Volume 46, 1923, pages 671-674.

    Geoffrey Taylor, A E Green,
    Mechanism for the production of small eddies from large ones,
    Proceedings of the Royal Society of London,
    Series A, Volume 158, 1937, pages 499-521.
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the density.
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Input, double T, the time coordinate or coordinates.
    Output, double UR[N], VR[N], PR[N], the residuals in the U,
    V and P equations.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * ur = s_data->a6;
	ityp * vr = s_data->a7;
	ityp * pr = s_data->a8;
	
    ityp dpdx;
    ityp dpdy;
    ityp dudt;
    ityp dudx;
    ityp dudxx;
    ityp dudy;
    ityp dudyy;
    ityp dvdt;
    ityp dvdx;
    ityp dvdxx;
    ityp dvdy;
    ityp dvdyy;
    ityp *f;
    ityp *g;
    ityp *h;
    int i;
    ityp p;
    ityp u;
    ityp v;

    /*
    Get the right hand sides.
    */
    f = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    g = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    h = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    rhs_taylor ( nu, rho, n, x, y, t, f, g, h );
    /*
    Form the functions and derivatives of the left hand side.
    */
    for ( i = 0; i < n; ++i )
    {
        u  =  -
        cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dudx =                       M_PI* sin ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dudxx =              M_PI * M_PI* cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dudy =  -                    M_PI* cos ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dudyy =              M_PI * M_PI* cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dudt =  + 2.00 * nu * M_PI * M_PI* cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );

        v  =
        sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dvdx =                       M_PI* cos ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dvdxx = -            M_PI * M_PI* sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dvdy =  -                    M_PI* sin ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dvdyy = -            M_PI * M_PI* sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dvdt =  - 2.00 * nu * M_PI * M_PI* sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );

        p =   - 0.25 * rho *( cos ( M_2TPI * x[i] ) + cos ( M_2TPI * y[i] ) );
        dpdx =  + 0.50  * rho * M_PI * sin ( M_2TPI * x[i] );
        dpdy =  + 0.50  * rho * M_PI * sin ( M_2TPI * y[i] );
        /*
        Time scaling.
        */
        u     = u     * exp ( - M_2TPI * M_PI * nu * t );
        dudx  = dudx  * exp ( - M_2TPI * M_PI * nu * t );
        dudxx = dudxx * exp ( - M_2TPI * M_PI * nu * t );
        dudy  = dudy  * exp ( - M_2TPI * M_PI * nu * t );
        dudyy = dudyy * exp ( - M_2TPI * M_PI * nu * t );
        dudt  = dudt  * exp ( - M_2TPI * M_PI * nu * t );

        v     = v     * exp ( - M_2TPI * M_PI * nu * t );
        dvdx  = dvdx  * exp ( - M_2TPI * M_PI * nu * t );
        dvdxx = dvdxx * exp ( - M_2TPI * M_PI * nu * t );
        dvdy  = dvdy  * exp ( - M_2TPI * M_PI * nu * t );
        dvdyy = dvdyy * exp ( - M_2TPI * M_PI * nu * t );
        dvdt  = dvdt  * exp ( - M_2TPI * M_PI * nu * t );

        p =     p     * exp ( - 4.00 * M_PI * M_PI * nu * t );
        dpdx =  dpdx  * exp ( - 4.00 * M_PI * M_PI * nu * t );
        dpdy =  dpdy  * exp ( - 4.00 * M_PI * M_PI * nu * t );
        /*
        Evaluate the residuals.
        */
        ur[i] = dudt + u * dudx + v * dudy+ ( 1.00 / rho ) * dpdx - nu * ( dudxx + dudyy ) - f[i];
        vr[i] = dvdt + u * dvdx + v * dvdy+ ( 1.00 / rho ) * dpdy - nu * ( dvdxx + dvdyy ) - g[i];
        pr[i] = dudx + dvdy - h[i];
    }
    /*
    Free memory.
    */
    free ( f );
    free ( g );
    free ( h );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _resid_vortex ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_VORTEX returns Vortex residuals.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2015
  Author:
    John Burkardt
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the density.
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Input, double T, the time coordinate or coordinates.
    Output, double UR[N], VR[N], PR[N], the residuals in the U,
    V and P equations.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * ur = s_data->a6;
	ityp * vr = s_data->a7;
	ityp * pr = s_data->a8;
	
    ityp dpdx;
    ityp dpdy;
    ityp dudt;
    ityp dudx;
    ityp dudxx;
    ityp dudy;
    ityp dudyy;
    ityp dvdt;
    ityp dvdx;
    ityp dvdxx;
    ityp dvdy;
    ityp dvdyy;
    ityp *f;
    ityp *g;
    ityp *h;
    dim_typ i;
    ityp p;
    ityp u;
    ityp v;

    /*
    Get the right hand sides.
    */
    f = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    g = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    h = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    rhs_taylor ( nu, rho, n, x, y, t, f, g, h );
    /*
    Form the functions and derivatives of the left hand side.
    */
    for ( i = 0; i < n; ++i )
    {
        u  =  -
        cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dudx =                       M_PI* sin ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dudxx =              M_PI * M_PI* cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dudy =  -                    M_PI* cos ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dudyy =              M_PI * M_PI* cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dudt =  + 2.00 * nu * M_PI * M_PI* cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );

        v  =
        sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dvdx =                       M_PI* cos ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dvdxx = -            M_PI * M_PI* sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dvdy =  -                    M_PI* sin ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        dvdyy = -            M_PI * M_PI* sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        dvdt =  - 2.00 * nu * M_PI * M_PI* sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );

        p =   - 0.25 * rho *( cos ( M_2TPI * x[i] ) + cos ( M_2TPI * y[i] ) );
        dpdx =  + 0.50  * rho * M_PI * sin ( M_2TPI * x[i] );
        dpdy =  + 0.50  * rho * M_PI * sin ( M_2TPI * y[i] );
        /*
        Evaluate the residuals.
        */
        ur[i] = dudt + u * dudx + v * dudy+ ( 1.00 / rho ) * dpdx - nu * ( dudxx + dudyy ) - f[i];
        vr[i] = dvdt + u * dvdx + v * dvdy+ ( 1.00 / rho ) * dpdy - nu * ( dvdxx + dvdyy ) - g[i];
        pr[i] = dudx + dvdy - h[i];
    }
    /*
    Free memory.
    */
    free ( f );
    free ( g );
    free ( h );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rhs_lucas ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_LUCAS evaluates Lucas Bystricky's right hand sides.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 March 2015
  Author:
    John Burkardt
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the fluid density.
    Input, int N, the number of nodes.
    Input, double X[N], Y[N], the coordinates of nodes.
    Input, double T, the current time.
    Output, double F[N], G[N], H[N], the right hand sides.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * f = s_data->a6;
	ityp * g = s_data->a7;
	ityp * h = s_data->a8;
	
    ityp *dpdx;
    ityp *dpdy;
    ityp *dudt;
    ityp *dudx;
    ityp *dudxx;
    ityp *dudy;
    ityp *dudyy;
    ityp *dvdt;
    ityp *dvdx;
    ityp *dvdxx;
    ityp *dvdy;
    ityp *dvdyy;
    dim_typ i;
    ityp *p;
    ityp *u;
    ityp *v;

    dpdx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dpdy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudt = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudxx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudyy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdt = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdxx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdyy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    p = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    u = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        u[i] = - cos ( M_PI * x[i] ) / M_PI;
        dudx[i] = sin ( M_PI * x[i] );
        dudxx[i] = M_PI * cos ( M_PI * x[i] );
        dudt[i] = dudx[i] = dudyy[i] = 0.00;
        v[i] = - y[i] * sin ( M_PI * x[i] );
        dvdx[i] = - M_PI * y[i] * cos ( M_PI * x[i] );
        dvdxx[i] = + M_PI * M_PI * y[i] * sin ( M_PI * x[i] );
        dvdy[i] = - sin ( M_PI * x[i] );

        p[i] = dvdt[i] = dvdyy[i] = dpdx[i] = dpdy[i] = 0.00;

        f[i] = dudt[i] - nu * ( dudxx[i] + dudyy[i] )+ u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho;
        g[i] = dvdt[i] - nu * ( dvdxx[i] + dvdyy[i] )+ u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho;
        h[i] = dudx[i] + dvdy[i];
    }

    free ( dpdx );
    free ( dpdy );
    free ( dudt );
    free ( dudx );
    free ( dudxx );
    free ( dudy );
    free ( dudyy );
    free ( dvdt );
    free ( dvdx );
    free ( dvdxx );
    free ( dvdy );
    free ( dvdyy );
    free ( p );
    free ( u );
    free ( v );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rhs_poiseuille ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_POISEUILLE evaluates Poiseuille right hand sides.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 July 2015
  Author:
    John Burkardt
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the fluid density.
    Input, int N, the number of nodes.
    Input, double X[N], Y[N], the coordinates of nodes.
    Input, double T, the current time.
    Output, double F[N], G[N], H[N], the right hand sides.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * f = s_data->a6;
	ityp * g = s_data->a7;
	ityp * h = s_data->a8;
	
    for (dim_typ i = 0; i < n; ++i )
        f[i] = g[i] = h[i] = 0.00;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rhs_spiral ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_SPIRAL returns Spiral right hand sides.
  Discussion:
    The right hand side is artificially determined by the requirement
    that the specified values of U, V and P satisfy the discretized
    Navier Stokes and continuity equations.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2015
  Author:
    John Burkardt
  Reference:
    Maxim Olshanskii, Leo Rebholz,
    Application of barycenter refined meshes in linear elasticity
    and incompressible fluid dynamics,
    ETNA: Electronic Transactions in Numerical Analysis,
    Volume 38, pages 258-274, 2011.
  Parameters:
    Input,double NU, the kinematic viscosity.
    Input, double RHO, the fluid density.
    Input, int N, the number of nodes.
    Input, double X[N], Y[N], the coordinates of nodes.
    Input, double T, the current time.
    Output, double F[N], G[N], H[N], the right hand sides.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * f = s_data->a6;
	ityp * g = s_data->a7;
	ityp * h = s_data->a8;
	
    ityp *dpdx;
    ityp *dpdy;
    ityp *dudt;
    ityp *dudx;
    ityp *dudxx;
    ityp *dudy;
    ityp *dudyy;
    ityp *dvdt;
    ityp *dvdx;
    ityp *dvdxx;
    ityp *dvdy;
    ityp *dvdyy;
    dim_typ i;
    ityp *p;
    ityp *u;
    ityp *v;

    dpdx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dpdy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudt = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudxx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dudyy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdt = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdxx = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    dvdyy = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    p = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    u = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    v = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i)
    {
        u[i] = ( 1.00 + nu * t ) * 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        dudt[i] = nu * 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        dudx[i] = ( 1.00 + nu * t ) * 2.00* ( 4.00 * pow ( x[i], 3 ) - 6.00 * pow ( x[i], 2 ) + 2.00 * x[i] )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        dudxx[i] = ( 1.00 + nu * t ) * 2.00* ( 12.00 * pow ( x[i], 2 ) - 12.00 * x[i] + 2.00 )* ( 2.00 * pow ( y[i], 3 ) - 3.00 * pow ( y[i], 2 ) + y[i] );
        dudy[i] = ( 1.00 + nu * t ) * 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 6.00 * pow ( y[i], 2 ) - 6.00 * y[i] + 1.00 );
        dudyy[i] = ( 1.00 + nu * t ) * 2.00* ( pow ( x[i], 4 ) - 2.00 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )* ( 12.00 * y[i] - 6.00 );
        v[i] = - ( 1.00 + nu * t ) * 2.00* ( 2.00 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );
        dvdt[i] = - nu * 2.00* ( 2.00 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );
        dvdx[i] = - ( 1.00 + nu * t ) * 2.00* ( 6.00 * pow ( x[i], 2 ) - 6.00 * x[i] + 1.00 )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );
        dvdxx[i] = - ( 1.00 + nu * t ) * 2.00* ( 12.00 * x[i] - 6.00 )* ( pow ( y[i], 4 ) - 2.00 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );
        dvdy[i] = - ( 1.00 + nu * t ) * 2.0* ( 2.00 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( 4.00 * pow ( y[i], 3 ) - 6.00 * pow ( y[i], 2 ) + 2.00 * y[i] );
        dvdyy[i] = - ( 1.00 + nu * t ) * 2.00* ( 2.00 * pow ( x[i], 3 ) - 3.00 * pow ( x[i], 2 ) + x[i] )* ( 12.00 * pow ( y[i], 2 ) - 12.00 * y[i] + 2.00 );
        p[i] = rho * y[i];
        dpdx[i] = 0.00;
        dpdy[i] = rho;
        f[i] = dudt[i] - nu * ( dudxx[i] + dudyy[i] )+ u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho;
        g[i] = dvdt[i] - nu * ( dvdxx[i] + dvdyy[i] )+ u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho;
        h[i] = dudx[i] + dvdy[i];
    }

    free ( dpdx );
    free ( dpdy );
    free ( dudt );
    free ( dudx );
    free ( dudxx );
    free ( dudy );
    free ( dudyy );
    free ( dvdt );
    free ( dvdx );
    free ( dvdxx );
    free ( dvdy );
    free ( dvdyy );
    free ( p );
    free ( u );
    free ( v );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rhs_taylor ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_TAYLOR returns Taylor right hand sides.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2015
  Author:
    John Burkardt
  Reference:
    Geoffrey Taylor,
    On the decay of vortices in a viscous fluid,
    Philosophical Magazine,
    Volume 46, 1923, pages 671-674.
    Geoffrey Taylor, A E Green,
    Mechanism for the production of small eddies from large ones,
    Proceedings of the Royal Society of London,
    Series A, Volume 158, 1937, pages 499-521.
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the density.
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Input, double T, the time coordinate or coordinates.
    Output, double F[N], G[N], H[N], the residuals in the U,
    V and P equations.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * f = s_data->a6;
	ityp * g = s_data->a7;
	ityp * h = s_data->a8;
	
    for (dim_typ i = 0; i < n; ++i)
        f[i] = g[i] = h[i] = 0.00;
    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _rhs_vortex ( void * data)
/******************************************************************************/
/*
  Purpose:
    RHS_VORTEX returns Vortex right hand sides.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2015
  Author:
    John Burkardt
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, double RHO, the density.
    Input, int N, the number of evaluation points.
    Input, double X[N], Y[N], the coordinates of the points.
    Input, double T, the time coordinate or coordinates.
    Output, double F[N], G[N], H[N], the residuals in the U,
    V and P equations.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	ityp nu = s_data->a0;
	ityp rho = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp t = s_data->a5;
	ityp * f = s_data->a6;
	ityp * g = s_data->a7;
	ityp * h = s_data->a8;

    for (dim_typ i = 0; i < n; ++i )
    {
        f[i] = - 2.00 * nu * M_PI * M_PI* cos ( M_PI * x[i] ) * sin ( M_PI * y[i] );
        g[i] =   2.00 * nu * M_PI * M_PI* sin ( M_PI * x[i] ) * cos ( M_PI * y[i] );
        h[i] = 0.00;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _resid_burgers ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_BURGERS evaluates the Burgers residual.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 July 2015
  Author:
    John Burkardt
  Reference:
    Martin Bazant, Henry Moffatt,
    Exact solutions of the Navier-Stokes equations having steady vortex structures,
    Journal of Fluid Mechanics,
    Volume 541, pages 55-64, 2005.
    Johannes Burgers,
    A mathematical model illustrating the theory of turbulence,
    Advances in Applied Mechanics,
    Volume 1, pages 171-199, 1948.
  Parameters:
    Input, double NU, the kinematic viscosity.
    Input, int N, the number of points at which the solution
    is to be evaluated.
    Input, double X[N], Y[N], Z[N], the coordinates of the points.
    Input, double T[N], the time coordinates.
    Output, double UR[N], VR[N], WR[N], PR[N], the residuals.
*/
{
	const itdt8pit * const s_data = data;
	const register ityp nu = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	ityp * y = s_data->a3;
	ityp * z = s_data->a4;
	ityp * t = s_data->a5;
	ityp * ur = s_data->a6;
	ityp * vr = s_data->a7;
	ityp * wr = s_data->a8;
	ityp * pr = s_data->a9;
	
    dim_typ i;
    ityp p;
    ityp px;
    ityp py;
    ityp pz;
    ityp u;
    ityp ut;
    ityp ux;
    ityp uxx;
    ityp uy;
    ityp uyy;
    ityp uz;
    ityp uzz;
    ityp v;
    ityp vt;
    ityp vx;
    ityp vxx;
    ityp vy;
    ityp vyy;
    ityp vz;
    ityp vzz;
    ityp w;
    ityp wt;
    ityp wx;
    ityp wxx;
    ityp wy;
    ityp wyy;
    ityp wz;
    ityp wzz;

    for ( i = 0; i < n; ++i )
    {
        /*
        Form the functions and derivatives.
        */
        u =   2.0 * x[i];
        ux =  2.0;

        uxx = uy = uyy = uz = uzz = ut = 0.00;

        v =   - 2.0 * y[i];
        vx = vxx = vyy = vz = vzz = vt = 0.00;
        vy =  - 2.0;

        w =   r8_erf ( y[i] / sqrt ( nu ) );
        wy =    2.00 * sqrt ( 1.00 / nu / M_PI )        * exp ( - y[i] * y[i] / nu );
        wyy = - 4.00 * sqrt ( 1.00 / nu / M_PI ) * y[i] * exp ( - y[i] * y[i] / nu ) / nu;
        wz =  wzz = wt = wx = wxx = 0.00;

        p = - 2.00 * ( x[i] * x[i] + y[i] * y[i] );
        px = - 4.00 * x[i];
        py = - 4.00 * y[i];
        pz = 0.00;
        /*
        Evaluate the residuals.
        */
        ur[i] = ut + u * ux + v * uy + w * uz + px - nu * ( uxx + uyy + uzz );
        vr[i] = vt + u * vx + v * vy + w * vz + py - nu * ( vxx + vyy + vzz );
        wr[i] = wt+ u * wx + v * wy + w * wz + pz - nu * ( wxx + wyy + wzz );
        pr[i] = ux + vy + wz;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _resid_ethier ( void * data)
/******************************************************************************/
/*
  Purpose:
    RESID_ETHIER evaluates the Ethier residual.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2015
  Author:
    John Burkardt
  Reference:
    C Ross Ethier, David Steinman,
    Exact fully 3D Navier-Stokes solutions for benchmarking,
    International Journal for Numerical Methods in Fluids,
    Volume 19, Number 5, March 1994, pages 369-375.
  Parameters:
    Input, double A, D, the parameters.  Sample values are A = M_PI/4
    and D = M_PI/2.
    Input, int N, the number of points at which the solution
    is to be evaluated.
    Input, double X[N], Y[N], Z[N], the coordinates of the points.
    Input, double T[N], the time coordinates.
    Output, double UR[N], VR[N], WR[N], PR[N], the residuals.
*/
{
	const _2itdt8pit * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp d = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp * z = s_data->a5;
	ityp * t = s_data->a6;
	ityp * ur = s_data->a7;
	ityp * vr = s_data->a8;
	ityp * wr = s_data->a9;
	ityp * pr = s_data->a10;
	
    ityp cxy;
    ityp cyz;
    ityp czx;
    ityp e2t;
    ityp e2x;
    ityp e2y;
    ityp e2z;
    ityp e4t;
    ityp ex;
    ityp exy;
    ityp ey;
    ityp eyz;
    ityp ez;
    ityp ezx;
    dim_typ i;
    ityp p;
    ityp px;
    ityp py;
    ityp pz;
    ityp sxy;
    ityp syz;
    ityp szx;
    ityp u;
    ityp ut;
    ityp ux;
    ityp uxx;
    ityp uy;
    ityp uyy;
    ityp uz;
    ityp uzz;
    ityp v;
    ityp vt;
    ityp vx;
    ityp vxx;
    ityp vy;
    ityp vyy;
    ityp vz;
    ityp vzz;
    ityp w;
    ityp wt;
    ityp wx;
    ityp wxx;
    ityp wy;
    ityp wyy;
    ityp wz;
    ityp wzz;
    /*
    Make some temporaries.
    */
    for ( i = 0; i < n; ++i )
    {
        ex = exp ( a * x[i] );
        ey = exp ( a * y[i] );
        ez = exp ( a * z[i] );

        e2x = exp ( 2.00 * a * x[i] );
        e2y = exp ( 2.00 * a * y[i] );
        e2z = exp ( 2.00 * a * z[i] );

        e2t = exp ( -       d * d * t[i] );
        e4t = exp ( - 2.00 * d * d * t[i] );

        exy = exp ( a * ( x[i] + y[i] ) );
        eyz = exp ( a * ( y[i] + z[i] ) );
        ezx = exp ( a * ( z[i] + x[i] ) );

        sxy = sin ( a * x[i] + d * y[i] );
        syz = sin ( a * y[i] + d * z[i] );
        szx = sin ( a * z[i] + d * x[i] );

        cxy = cos ( a * x[i] + d * y[i] );
        cyz = cos ( a * y[i] + d * z[i] );
        czx = cos ( a * z[i] + d * x[i] );
        /*
        Form the functions and derivatives.
        */
        u =   -         a * (           ex * syz+         ez * cxy ) * e2t;
        ux =  -         a * (       a * ex * syz-     a * ez * sxy ) * e2t;
        uxx = -         a * (   a * a * ex * syz- a * a * ez * cxy ) * e2t;
        uy =  -         a * (       a * ex * cyz-     d * ez * sxy ) * e2t;
        uyy = -         a * ( - a * a * ex * syz- d * d * ez * cxy ) * e2t;
        uz =  -         a * (       d * ex * cyz+     a * ez * cxy ) * e2t;
        uzz =  -        a * ( - d * d * ex * syz+ a * a * ez * cxy ) * e2t;
        ut =  + d * d * a * (           ex * syz+         ez * cxy ) * e2t;

        v =   -         a * (           ey * szx+         ex * cyz ) * e2t;
        vx =  -         a * (       d * ey * czx+     a * ex * cyz ) * e2t;
        vxx = -         a * ( - d * d * ey * szx+ a * a * ex * cyz ) * e2t;
        vy =  -         a * (       a * ey * szx-     a * ex * syz ) * e2t;
        vyy = -         a * (   a * a * ey * szx- a * a * ex * cyz ) * e2t;
        vz =  -         a * (       a * ey * czx-     d * ex * syz ) * e2t;
        vzz =  -        a * ( - a * a * ey * szx- d * d * ex * cyz ) * e2t;
        vt =  + d * d * a * (           ey * szx+         ex * cyz ) * e2t;

        w =   -         a * (           ez * sxy+         ey * czx ) * e2t;
        wx =  -         a * (       a * ez * cxy-     d * ey * szx ) * e2t;
        wxx = -         a * ( - a * a * ez * sxy- d * d * ey * czx ) * e2t;
        wy =  -         a * (       d * ez * cxy+     a * ey * czx ) * e2t;
        wyy = -         a * ( - d * d * ez * sxy+ a * a * ey * czx ) * e2t;
        wz =  -         a * (       a * ez * sxy-     a * ey * szx ) * e2t;
        wzz = -         a * (   a * a * ez * sxy- a * a * ey * czx ) * e2t;
        wt =  + d * d * a * (           ez * sxy+         ey * czx ) * e2t;

        p = - 0.50 * a * a * e4t * (+ e2x + 2.00 * sxy * czx * eyz+ e2y + 2.00 * syz * cxy * ezx+ e2z + 2.00 * szx * cyz * exy );
        px = - 0.50 * a * a * e4t * (+ 2.00 * a * e2x+ 2.00 * a * cxy * czx * eyz- 2.00 * d * sxy * szx * eyz- 2.00 * a * syz * sxy * ezx+ 2.00 * a * syz * cxy * ezx+ 2.00 * d * czx * cyz * exy+ 2.00 * a * szx * cyz * exy );
        py = - 0.50 * a * a * e4t * (+ 2.00 * d * cxy * czx * eyz+ 2.00 * a * sxy * czx * eyz+ 2.00 * a * e2y+ 2.00 * a * cyz * cxy * ezx- 2.00 * d * syz * sxy * ezx- 2.00 * a * szx * syz * exy+ 2.00 * a * szx * cyz * exy );
        pz = - 0.50 * a * a * e4t * (- 2.00 * a * sxy * szx * eyz+ 2.00 * a * sxy * czx * eyz+ 2.00 * d * cyz * cxy * ezx+ 2.00 * a * syz * cxy * ezx+ 2.00 * a * e2z+ 2.00 * a * czx * cyz * exy- 2.00 * d * szx * syz * exy );
        /*
        Evaluate the residuals.
        */
        ur[i] = ut+ u * ux + v * uy + w * uz + px- ( uxx + uyy + uzz );
        vr[i] = vt+ u * vx + v * vy + w * vz + py- ( vxx + vyy + vzz );
        wr[i] = wt+ u * wx + v * wy + w * wz + pz- ( wxx + wyy + wzz );
        pr[i] = ux + vy + wz;
    }

    return NULL;
}

#endif
