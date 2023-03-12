#ifndef __DISABLEDEEP_FEM2DBVPQUADRATIC

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   *fem2d_bvp_quadratic (const register dim_typ nx, const register dim_typ ny, ityp a ( ityp x, ityp y ),
  ityp c ( ityp x, ityp y ), ityp f ( ityp x, ityp y ), ityp x[static nx], ityp y[static ny] )
/******************************************************************************/
/*
  Purpose:
    FEM2D_BVP_QUADRATIC solves a boundary value problem on a rectangle.
  Discussion:
    The procedure uses the finite element method, with piecewise quadratic
    basis functions to solve a 2D boundary value problem over a rectangle
    The following differential equation is imposed inside the region:
      - d/dx a(x,y) du/dx - d/dy a(x,y) du/dy + c(x,y) * u(x,y) = f(x,y)
    where a(x,y), c(x,y), and f(x,y) are given functions.
    On the boundary, the solution is constrained to have the value 0.
    The finite element method will use a regular grid of NX nodes in X, and
    NY nodes in Y.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int NX, NY, the number of X and Y grid values.
    NX and NY must be odd and at least 3.
    Input, double A ( double X, double Y ), evaluates a(x,y);
    Input, double C ( double X, double Y ), evaluates c(x,y);
    Input, double F ( double X, double Y ), evaluates f(x,y);
    Input, double X[NX], Y[NY], the grid coordinates.
    Output, double FEM1D_BVP_QUADRATIC[NX*NY], the finite element coefficients,
    which are also the value of the computed solution at the mesh points.
*/
{
    # define QUAD_NUM 3

    ityp abscissa[QUAD_NUM] =
	{
        -0.774596669241483377035853079956,
        0.000000000000000000000000000000,
        0.774596669241483377035853079956
    };
    ityp *amat;
    ityp aq;
    ityp *b;
    dim_typ cc;
    ityp cq;
    dim_typ e;
    dim_typ ex;
    dim_typ ex_num;
    dim_typ ey;
    dim_typ ey_num;
    ityp fq;
    dim_typ i;
    dim_typ ierror;
    dim_typ ii;
    dim_typ il;
    dim_typ il2;
    dim_typ il3;
    dim_typ j;
    dim_typ jj;
    dim_typ jl;
    dim_typ jl2;
    dim_typ jl3;
    dim_typ k;
    dim_typ mm;
    dim_typ mn = nx*ny;
    dim_typ n;
    dim_typ node[9];
    dim_typ quad_num = QUAD_NUM;
    dim_typ qx;
    dim_typ qy;
    dim_typ s;
    ityp t;
    ityp *u;
    ityp v[9];
    ityp vx[9];
    ityp vy[9];
    dim_typ w;
    ityp weight[QUAD_NUM] =
    {
        0.555555555555555555555555555556,
        0.888888888888888888888888888889,
        0.555555555555555555555555555556
    };
    ityp wq;
    ityp xq;
    ityp xx[3];
    ityp yq;
    ityp yy[3];

    amat = r8mat_zero_new ( mn, mn );
    b = r8vec_zero_new ( mn );

    ex_num = ( nx - 1 ) / 2;
    ey_num = ( ny - 1 ) / 2;

    for ( ex = 0; ex < ex_num; ++ex )
    {
        w = ex<<1;
        cc = (ex<<1) + 1;
        e = (ex<<1) + 2;

        xx[0] = x[w];
        xx[1] = x[cc];
        xx[2] = x[e];

        for ( ey = 0; ey < ey_num; ++ey )
        {
            s = ey<<1;
            mm = (ey<<1) + 1;
            n = (ey<<1)+ 2;

            yy[0] = y[s];
            yy[1] = y[mm];
            yy[2] = y[n];
            /*
            Node indices

            7  8  9   wn cn en
            4  5  6   wm cm em
            1  2  3   ws cs es
            */
            node[0] = ( (ey<<1)     ) * nx + (ex<<1) + 0;
            node[1] = ( (ey<<1)     ) * nx + (ex<<1) + 1;
            node[2] = ( (ey<<1)     ) * nx + (ex<<1) + 2;
            node[3] = ( (ey<<1) + 1 ) * nx + (ex<<1) + 0;
            node[4] = ( (ey<<1) + 1 ) * nx + (ex<<1) + 1;
            node[5] = ( (ey<<1) + 1 ) * nx + (ex<<1) + 2;
            node[6] = ( (ey<<1) + 2 ) * nx + (ex<<1) + 0;
            node[7] = ( (ey<<1) + 2 ) * nx + (ex<<1) + 1;
            node[8] = ( (ey<<1) + 2 ) * nx + (ex<<1) + 2;

            for ( qx = 0; qx < quad_num; ++qx )
            {
                xq = ( ( 1.00 - abscissa[qx] ) * xx[0]+ ( 1.00 + abscissa[qx] ) * xx[2] ) / 2.00;

                for ( qy = 0; qy < quad_num; ++qy )
                {
                    yq = ( ( 1.00 - abscissa[qy] ) * yy[0]   + ( 1.00 + abscissa[qy] ) * yy[2] ) / 2.00;
                    wq = weight[qx] * ( xx[2] - xx[0] ) / 2.00 * weight[qy] * ( yy[2] - yy[0] ) / 2.00;
                    k = 0;

                    for ( jl = 0; jl < 3; ++jl )
                    {
                        for ( il = 0; il < 3; ++il )
                        {
                            v[k] = 1.00;
                            vx[k] = vy[k] = 0.00;
                            for ( il2 = 0; il2 < 3; ++il2)
                            {
                                if ( il2 != il )
                                {
                                    v[k] *= ( xq - xx[il2] ) / ( xx[il] - xx[il2] );
                                    t = 1.00 / ( xx[il] - xx[il2] );
                                    for ( il3 = 0; il3 < 3; ++il3 )
                                        if ( il3 != il && il3 != il2 )
                                            t *= ( xq - xx[il3] ) / ( xx[il] - xx[il3] );
                                    for ( jl2 = 0; jl2 < 3; ++jl2  )
                                        if ( jl2 != jl )
                                            t *= ( yq - yy[jl2] ) / ( yy[jl] - yy[jl2] );
                                    vx[k] += t;
                                }
                            }

                            for ( jl2 = 0; jl2 < 3; ++jl2 )
                            {
                                if ( jl2 != jl )
                                {
                                    v[k] = v[k] * ( yq - yy[jl2] ) / ( yy[jl] - yy[jl2] );
                                    t = 1.00 / ( yy[jl] - yy[jl2] );
                                    for ( il2 = 0; il2 < 3; ++il2 )
                                        if ( il2 != il )
                                            t *= ( xq - xx[il2] ) / ( xx[il] - xx[il2] );
                                    for ( jl3 = 0; jl3 < 3; jl3++ )
                                        if ( jl3 != jl && jl3 != jl2 )
                                            t = t * ( yq - yy[jl3] ) / ( yy[jl] - yy[jl3] );
                                    vy[k] +=t;
                                }
                            }
                            ++ k;
                        }
                    }

                    aq = a ( xq, yq );
                    cq = c ( xq, yq );
                    fq = f ( xq, yq );

                    for ( i = 0; i < 9; ++i )
                    {
                        ii = node[i];
                        for ( j = 0; j < 9; ++j )
                        {
                            jj = node[j];
                            amat[ii+jj*mn] = amat[ii+jj*mn] + wq * (
                            vx[i] * aq * vx[j] + vy[i] * aq * vy[j] + v[i]  * cq * v[j] );
                        }
                        b[ii] += wq * ( v[i] * fq );
                    }
                }
            }
        }
    }
    /*
    Where a node is on the boundary,
    replace the finite element equation by a boundary condition.
    */
    k = 0;
    for ( j = 0; j < ny; ++j )
    {
        for ( i = 0; i < nx; ++i )
        {
        if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
        {
            for ( jj = 0; jj < mn; ++j )
                amat[k+jj*mn] = 0.00;
            for ( ii = 0; ii < mn; ++ii )
                amat[ii+k*mn] = 0.00;
            amat[k+k*mn] = 1.00;
            b[k] = 0.00;
            }
            ++ k;
        }
    }
    /*
    Solve the linear system.
    */
    u = r8mat_solve2 ( mn, amat, b, &ierror );
    free ( amat );
    free ( b );
    return u;
    # undef QUAD_NUM
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   fem2d_h1s_error_quadratic (const register dim_typ nx, const register dim_typ ny, ityp x[static nx], ityp y[static ny],
  ityp u[static nx*ny], ityp exact_ux ( ityp x, ityp y ),ityp exact_uy ( ityp x, ityp y ) )
/******************************************************************************/
/*
  Purpose:
    FEM2D_H1S_ERROR_QUADRATIC: seminorm error of a finite element solution.
  Discussion:
    The finite element method has been used, over a rectangle,
    involving a grid of NX*NY nodes, with piecewise quadratic elements used
    for the basis.
    The finite element solution U(x,y) has been computed, and formulas for the
    exact derivatives Vx(x,y) and Vy(x,y) are known.
    This function estimates the H1 seminorm of the error:
      H1S = sqrt ( integral ( x, y ) ( Ux(x,y) - Vx(x,y) )^2
                                     + ( Uy(x,y) - Vy(x,y) )^2 dx dy )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int NX, NY, the number of X and Y grid values.
    Input, double X[NX], Y[NY], the grid coordinates.
    Input, double U[NX*NY], the finite element coefficients.
    Input, function EXACT_UX(X,Y), EXACT_UY(X,Y) return the
    value of the derivatives of the exact solution with respect to
    X and Y, respectively, at the point (X,Y).
    Output, double FEM2D_H1S_ERROR_QUADRATIC, the estimated seminorm of
    the error.
*/
{
    # define QUAD_NUM 3

    ityp abscissa[QUAD_NUM] =
    {
        -0.774596669241483377035853079956,
        0.000000000000000000000000000000,
        0.774596669241483377035853079956
    };
    ityp aq;
    dim_typ cc;
    ityp cq;
    dim_typ e;
    dim_typ ex;
    dim_typ ex_num;
    ityp exq;
    ityp eyq;
    dim_typ ey;
    dim_typ ey_num;
    ityp fq;
    ityp h1s;
    dim_typ i;
    dim_typ ierror;
    dim_typ ii;
    dim_typ il;
    dim_typ il2;
    dim_typ il3;
    dim_typ j;
    dim_typ jj;
    dim_typ jl;
    dim_typ jl2;
    dim_typ jl3;
    dim_typ k;
    dim_typ mm;
    dim_typ mn = nx*ny;
    dim_typ n;
    dim_typ node[9];
    dim_typ quad_num = QUAD_NUM;
    dim_typ qx;
    dim_typ qy;
    dim_typ s;
    ityp t;
    ityp uxq;
    ityp uyq;
    ityp vx[9];
    ityp vy[9];
    dim_typ w;
    ityp weight[QUAD_NUM] =
    {
        0.555555555555555555555555555556,
        0.888888888888888888888888888889,
        0.555555555555555555555555555556
    };
    ityp wq;
    ityp xq;
    ityp xx[3];
    ityp yq;
    ityp yy[3];

    ex_num = ( nx - 1 ) / 2;
    ey_num = ( ny - 1 ) / 2;

    for ( ex = 0; ex < ex_num; ++ex )
    {
        w = ex<<1;
        cc = w + 1;
        e = w + 2;

        xx[0] = x[w];
        xx[1] = x[cc];
        xx[2] = x[e];

        for ( ey = 0; ey < ey_num; ey++ )
        {
            s = ey<<1;
            mm = s + 1;
            n = s + 2;

            yy[0] = y[s];
            yy[1] = y[mm];
            yy[2] = y[n];
            /*
            Node indices

            7  8  9   wn cn en
            4  5  6   wm cm em
            1  2  3   ws cs es
            */
            node[0] = ( s     ) * nx + w + 0;
            node[1] = ( s     ) * nx + w + 1;
            node[2] = ( s     ) * nx + w + 2;
            node[3] = ( s + 1 ) * nx + w + 0;
            node[4] = ( s + 1 ) * nx + w + 1;
            node[5] = ( s + 1 ) * nx + w + 2;
            node[6] = ( s + 2 ) * nx + w + 0;
            node[7] = ( s + 2 ) * nx + w + 1;
            node[8] = ( s + 2 ) * nx + w + 2;

            for ( qx = 0; qx < quad_num; ++qx )
            {
                xq = ( ( 1.00 - abscissa[qx] ) * xx[0]+ ( 1.00 + abscissa[qx] ) * xx[2] )/ 2.0;

                for ( qy = 0; qy < quad_num; ++qy )
                {
                    yq = ( ( 1.00 - abscissa[qy] ) * yy[0]+ ( 1.00 + abscissa[qy] ) * yy[2] )/ 2.00;
                    wq = weight[qx] * ( xx[2] - xx[0] ) / 2.00* weight[qy] * ( yy[2] - yy[0] ) / 2.00;
                    uxq = uyq = 0.00;

                    k = 0;

                    for ( jl = 0; jl < 3; ++jl )
                    {
                        for ( il = 0; il < 3; ++il )
                        {
                            vx[k] = vy[k] = 0.00;
                            for ( il2 = 0; il2 < 3; ++il2 )
                            {
                                if ( il2 != il )
                                {
                                    t = 1.00 / ( xx[il] - xx[il2] );
                                    for ( il3 = 0; il3 < 3; ++il3 )
                                        if ( il3 != il && il3 != il2 )
                                            t *= ( xq - xx[il3] ) / ( xx[il] - xx[il3] );
                                    for ( jl2 = 0; jl2 < 3; jl2++ )
                                        if ( jl2 != jl )
                                            t *= ( yq - yy[jl2] ) / ( yy[jl] - yy[jl2] );
                                    vx[k] += t;
                                }
                            }

                            for ( jl2 = 0; jl2 < 3; ++jl2 )
                            {
                                if ( jl2 != jl )
                                {
                                    t = 1.00 / ( yy[jl] - yy[jl2] );
                                    for ( il2 = 0; il2 < 3; ++il2 )
                                        if ( il2 != il )
                                            t *= ( xq - xx[il2] ) / ( xx[il] - xx[il2] );
                                    for ( jl3 = 0; jl3 < 3; ++jl3 )
                                        if ( jl3 != jl && jl3 != jl2 )
                                            t *= ( yq - yy[jl3] ) / ( yy[jl] - yy[jl3] );

                                    vy[k] += t;
                                }
                            }
                            uxq += u[node[k]] * vx[k];
                            uyq += u[node[k]] * vy[k];
                            ++ k;
                        }
                    }

                    exq = exact_ux ( xq, yq );
                    eyq = exact_uy ( xq, yq );

                    h1s += wq * ( pow ( uxq - exq, 2 ) + pow ( uyq - eyq, 2 ) );
                }
            }
        }
    }
    return sqrt ( h1s );
    # undef QUAD_NUM
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   fem2d_l2_error_quadratic (const register dim_typ nx, const register dim_typ ny, ityp x[static nx], ityp y[static ny],
  ityp u[static nx*ny], ityp exact ( ityp x, ityp y ) )
/******************************************************************************/
/*
  Purpose:
    FEM2D_L2_ERROR_QUADRATIC: L2 error norm of a finite element solution.
  Discussion:
    The finite element method has been used, over a rectangle,
    involving a grid of NX*NY nodes, with piecewise quadratic elements used
    for the basis.
    The finite element coefficients have been computed, and a formula for the
    exact solution is known.
    This function estimates E2, the L2 norm of the error:
      E2 = Integral ( X, Y ) ( U(X,Y) - EXACT(X,Y) )^2 dX dY
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 June 2014
  Author:
    John Burkardt
  Parameters:
    Input, int NX, NY, the number of X and Y grid values.
    Input, double X[NX], Y[NY], the grid coordinates.
    Input, double U[NX*NY], the finite element coefficients.
    Input, function EQ = EXACT(X,Y), returns the value of the exact
    solution at the point (X,Y).
    Output, double FEM2D_L2_ERROR_QUADRATIC, the estimated L2 norm of the error.
*/
{
    # define QUAD_NUM 3

    ityp abscissa[QUAD_NUM] =
    {
        -0.774596669241483377035853079956,
        0.000000000000000000000000000000,
        0.774596669241483377035853079956
    };
    ityp aq;
    dim_typ cc;
    ityp cq;
    dim_typ e;
    ityp e2;
    ityp eq;
    dim_typ ex;
    dim_typ ex_num;
    dim_typ ey;
    dim_typ ey_num;
    ityp fq;
    dim_typ i;
    dim_typ ierror;
    dim_typ ii;
    dim_typ il;
    dim_typ il2;
    dim_typ il3;
    dim_typ j;
    dim_typ jj;
    dim_typ jl;
    dim_typ jl2;
    dim_typ jl3;
    dim_typ k;
    dim_typ mm;
    dim_typ mn = nx*ny;
    dim_typ n;
    dim_typ node[9];
    dim_typ quad_num = QUAD_NUM;
    dim_typ qx;
    dim_typ qy;
    dim_typ s;
    ityp t;
    ityp uq;
    ityp v[9];
    dim_typ w;
    ityp weight[QUAD_NUM] =
    {
        0.555555555555555555555555555556,
        0.888888888888888888888888888889,
        0.555555555555555555555555555556
    };
    ityp wq;
    ityp xq;
    ityp xx[3];
    ityp yq;
    ityp yy[3];

    ex_num = ( nx - 1 ) / 2;
    ey_num = ( ny - 1 ) / 2;

    for ( ex = 0; ex < ex_num; ++e )
    {
        w = ex<<1;
        cc = w + 1;
        e = w + 2;

        xx[0] = x[w];
        xx[1] = x[cc];
        xx[2] = x[e];

        for ( ey = 0; ey < ey_num; ++ey )
        {
            s = ey<<1;
            mm = s + 1;
            n = s + 2;

            yy[0] = y[s];
            yy[1] = y[mm];
            yy[2] = y[n];
            /*
            Node indices

            7  8  9   wn cn en
            4  5  6   wm cm em
            1  2  3   ws cs es
            */
            node[0] = ( s     ) * nx + w + 0;
            node[1] = ( s     ) * nx + w + 1;
            node[2] = ( s     ) * nx + w + 2;
            node[3] = ( s + 1 ) * nx + w + 0;
            node[4] = ( s + 1 ) * nx + w + 1;
            node[5] = ( s + 1 ) * nx + w + 2;
            node[6] = ( s + 2 ) * nx + w;
            node[7] = ( s + 2 ) * nx + w + 1;
            node[8] = ( s + 2 ) * nx + w + 2;

            for ( qx = 0; qx < quad_num; ++qx )
            {
                xq = ( ( 1.00 - abscissa[qx] ) * xx[0]+ ( 1.00 + abscissa[qx] ) * xx[2] )/ 2.00;

                for ( qy = 0; qy < quad_num; qy++ )
                {
                    yq = ( ( 1.00 - abscissa[qy] ) * yy[0]+ ( 1.00 + abscissa[qy] ) * yy[2] )/ 2.00;
                    wq = weight[qx] * ( xx[2] - xx[0] ) / 2.00* weight[qy] * ( yy[2] - yy[0] ) / 2.00;
                    uq = 0.00;
                    k = 0;

                    for ( jl = 0; jl < 3; ++jl )
                    {
                        for ( il = 0; il < 3; il++ )
                        {
                            v[k] = 1.00;
                            for ( il2 = 0; il2 < 3; ++il2 )
                                if ( il2 != il )
                                    v[k] *= ( xq - xx[il2] ) / ( xx[il] - xx[il2] );

                            for ( jl2 = 0; jl2 < 3; ++jl2 )
                                if ( jl2 != jl )
                                    v[k] *= ( yq - yy[jl2] ) / ( yy[jl] - yy[jl2] );
                            uq += u[node[k]] * v[k];
                            ++ k;
                        }
                    }

                    eq = exact ( xq, yq );

                    e2 += wq * pow ( uq - eq, 2 );
                }
            }
        }
    }
    return sqrt ( e2 );
    # undef QUAD_NUM
}
#endif
