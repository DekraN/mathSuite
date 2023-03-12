#ifndef __DISABLEDEEP_PIECEWISELINEARPRODUCTINTEGRAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _piecewise_linear_product_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    PIECEWISE_LINEAR_PRODUCT_INTEGRAL: piecewise linear product integral.
  Discussion:
    We are given two piecewise linear functions F(X) and G(X) and we wish
    to compute the exact value of the integral
      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
    The functions F(X) and G(X) are defined as tables of coordinates X and
    values V.  A piecewise linear function is evaluated at a point X by
    evaluating the interpolant to the data at the endpoints of the interval
    containing X.
    It must be the case that A <= B.
    It must be the case that the node coordinates F_X(*) and G_X(*) are
    given in ascending order.
    It must be the case that:
      F_X(1) <= A and B <= F_X(F_NUM)
      G_X(1) <= A and B <= G_X(G_NUM)

  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, the limits of integration.
    Input, int F_NUM, the number of nodes for F.
    Input, double F_X[F_NUM], the node coordinates for F.
    Input, double F_V[F_NUM], the nodal values for F.
    Input, int G_NUM, the number of nodes for G.
    Input, double G_X[G_NUM], the node coordinates for G.
    Input, double G_V[G_NUM], the nodal values for G.
    Output, double INTEGRAL, the integral of F(X) * G(X)
    from A to B.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dt4pit2it * const s_data = data;
	
	const register dim_typ f_num = s_data->a0; 
	const register dim_typ g_num = s_data->a1; 
	ityp * f_x = s_data->a2;
	ityp * f_v = s_data->a3;
	ityp * g_x = s_data->a4;
	ityp * g_v = s_data->a5;
	ityp a = s_data->a6;
	ityp b = s_data->a7;
	
    ityp bit;
    int f_left;
    ityp f0;
    ityp f1;
    ityp fl;
    ityp fr;
    int g_left;
    ityp g0;
    ityp g1;
    ityp gl;
    ityp gr;
    ityp h0;
    ityp h1;
    ityp h2;
    dim_typ i;
    ityp integral;
    ityp xl;
    ityp xr;
    ityp xr_max;

    if ( f_x[f_num-1] <= a || g_x[g_num-1] <= a || f_num < 2 || g_num < 2 )
    {
    	result = 0.00; 
		return &result;
	}

    integral = 0.00;
    xr = a;

    f_left = 0;
    r8vec_bracket3 ( f_num, f_x, xr, &f_left );
    fr = f_v[f_left] + ( xr - f_x[f_left] ) * ( f_v[f_left+1] - f_v[f_left] ) / ( f_x[f_left+1] - f_x[f_left] );

    g_left = 0;
    r8vec_bracket3 ( g_num, g_x, xr, &g_left );
    gr = g_v[g_left] + ( xr - g_x[g_left] ) * ( g_v[g_left+1] - g_v[g_left] ) / ( g_x[g_left+1] - g_x[g_left] );

    xr_max = b;
    xr_max = MIN ( xr_max, f_x[f_num-1] );
    xr_max = MIN ( xr_max, g_x[g_num-1] );

    while ( xr < xr_max )
    {
        /*
        Shift right values to left.
        */
        xl = xr;
        fl = fr;
        gl = gr;
        /*
        Determine the new right values.
        The hard part is figuring out how to advance XR some, but not too much.
        */
        xr = xr_max;

        for ( i = 1; i <= 2; ++i )
            if ( f_left + i <= f_num - 1 )
                if ( xl < f_x[f_left+i] && f_x[f_left+i] < xr )
                {
                    xr = f_x[f_left+i];
                    break;
                }

        for ( i = 1; i <= 2; ++i)
            if ( g_left + i <= g_num - 1 )
                if ( xl < g_x[g_left+i] && g_x[g_left+i] < xr )
                {
                    xr = g_x[g_left+i];
                    break;
                }

        r8vec_bracket3 ( f_num, f_x, xr, &f_left );
        fr = f_v[f_left] + ( xr - f_x[f_left] ) * ( f_v[f_left+1] - f_v[f_left] ) / ( f_x[f_left+1] - f_x[f_left] );

        r8vec_bracket3 ( g_num, g_x, xr, &g_left );
        gr = g_v[g_left] + ( xr - g_x[g_left] ) * ( g_v[g_left+1] - g_v[g_left] ) / ( g_x[g_left+1] - g_x[g_left] );
        /*
        Form the linear polynomials for F(X) and G(X) over [XL,XR],
        then the product H(X), integrate H(X) and add to the running total.
        */
        if ( r8_epsilon ( ) <= abs ( xr - xl ) )
        {
            f1 = fl - fr;
            f0 = fr * xl - fl * xr;

            g1 = gl - gr;
            g0 = gr * xl - gl * xr;

            h2 = f1 * g1;
            h1 = f1 * g0 + f0 * g1;
            h0 = f0 * g0;

            h2 /= 3.00;
            h1 /= 2.00;

            bit = ( ( h2 * xr + h1 ) * xr + h0 ) * xr - ( ( h2 * xl + h1 ) * xl + h0 ) * xl;

            integral += bit / ( xr - xl ) / ( xr - xl );
        }
    }

	result = integral;
    return &result;
}

#endif
