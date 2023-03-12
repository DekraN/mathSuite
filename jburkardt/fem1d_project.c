#ifndef __DISABLEDEEP_FEM1DPROJECT

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _piecewise_linear_product_quad ( void * data)
/******************************************************************************/
/*
  Purpose:
    PIECEWISE_LINEAR_PRODUCT_QUAD: estimate piecewise linear product integral.
  Discussion:
    We are given two piecewise linear functions F(X) and G(X) and we wish
    to estimate the value of the integral
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
    Input, int QUAD_NUM, the number of quadrature points.
    Output, double PIECEWISE_LINEAR_PRODUCT_QUAD, an estimate for the integral
    of F(X) * G(X) from A to B.
*/
{
	static ityp result = MAX_VAL;
	
	const _2itdt2pitdt2pitdt * const s_data = data;
	ityp a = s_data->a0;
	ityp b = s_data->a1;
	const register dim_typ f_num = s_data->a2;
	ityp * f_x = s_data->a3;
	ityp * f_v = s_data->a4;
	const register dim_typ g_num = s_data->a5;
	ityp * g_x = s_data->a6;
	ityp * g_v = s_data->a7;
	const register dim_typ quad_num = s_data->a8;
	
    ityp a2;
    ityp b2;
    int f_left;
    ityp fq;
    int g_left;
    ityp gq;
    dim_typ i;
    ityp quad;
    ityp xq;

    quad = 0.00;
    f_left = g_left = 0;

    a2 = a;
    a2 = MAX ( a2, f_x[0] );
    a2 = MAX ( a2, g_x[0] );

    b2 = b;
    b2 = MIN ( b2, f_x[f_num-1] );
    b2 = MIN ( b2, g_x[g_num-1] );

    for ( i = 1; i <= quad_num; ++i )
    {
        xq = ( ( ityp ) (             (i<<1) - 1 ) * b2+ ( ityp ) ( (quad_num<<1) - (i<<1) + 1 ) * a2 )/ ( ityp ) ( (quad_num<<1)             );
        r8vec_bracket3 ( f_num, f_x, xq, &f_left );
        fq = f_v[f_left] + ( xq - f_x[f_left] ) * ( f_v[f_left+1] - f_v[f_left] )/ ( f_x[f_left+1] - f_x[f_left] );
        r8vec_bracket3 ( g_num, g_x, xq, &g_left );
        gq = g_v[g_left] + ( xq - g_x[g_left] ) * ( g_v[g_left+1] - g_v[g_left] )/ ( g_x[g_left+1] - g_x[g_left] );
        quad += fq * gq;
    }

	result = quad * ( b - a ) / ( ityp ) ( quad_num );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_bracket3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_BRACKET3 finds the interval containing or nearest a given value.
  Discussion:
    An R8VEC is a vector of R8's.
    The routine always returns the index LEFT of the sorted array
    T with the property that either
    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
    *  T < T[LEFT] = T[0], or
    *  T > T[LEFT+1] = T[N-1].
    The routine is useful for interpolation problems, where
    the abscissa must be located within an interval of data
    abscissas for interpolation, or the "nearest" interval
    to the (extreme) abscissa must be found so that extrapolation
    can be carried out.
    This version of the function has been revised so that the value of
    LEFT that is returned uses the 0-based indexing natural to C++.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, length of the input array.
    Input, double T[N], an array that has been sorted into ascending order.
    Input, double TVAL, a value to be bracketed by entries of T.
    Input/output, int *LEFT.
    On input, if 0 <= LEFT <= N-2, LEFT is taken as a suggestion for the
    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
    is searched first, followed by the appropriate interval to the left
    or right.  After that, a binary search is used.
    On output, LEFT is set so that the interval [ T[LEFT], T[LEFT+1] ]
    is the closest to TVAL; it either contains TVAL, or else TVAL
    lies outside the interval [ T[0], T[N-1] ].
*/
{
	const dtpititpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * t = s_data->a1; 
	const register ityp tval = s_data->a2;
	int * left = s_data->a3;
	
    int high;
    int low;
    int mid;
    /*
    Check the input data.
    */
    if ( n < 2 )
        return NULL;
    /*
    If *LEFT is not between 0 and N-2, set it to the middle value.
    */
    if ( *left < 0 || n - 2 < *left )
        *left = ( n - 1 ) / 2;
    /*
    CASE 1: TVAL < T[*LEFT]:
    Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
    */
    if ( tval < t[*left] )
    {
        if ( *left == 0 )
            return NULL;
        else if ( *left == 1 )
        {
            *left = 0;
            return NULL;
        }
        else if ( t[*left-1] <= tval )
        {
            -- *left;
            return NULL;
        }
        else if ( tval <= t[1] )
        {
            *left = 0;
            return NULL;
        }
        /*
        ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
        */
        low = 1;
        high = *left - 2;

        for ( ; ; )
        {
            if ( low == high )
            {
                *left = low;
                return NULL;
            }

            mid = ( low + high + 1 ) / 2;

            if ( t[mid] <= tval )
                low = mid;
            else
                high = mid - 1;

        }
    }
    /*
    CASE 2: T[*LEFT+1] < TVAL:
    Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
    */
    else if ( t[*left+1] < tval )
    {
        if ( *left == n - 2 )
            return NULL;
        else if ( *left == n - 3 )
        {
            ++ *left;
            return NULL;
        }
        else if ( tval <= t[*left+2] )
        {
            ++ *left;
            return NULL;
        }
        else if ( t[n-2] <= tval )
        {
            *left = n - 2;
            return NULL;
        }
        /*
        ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
        */
        low = *left + 2;
        high = n - 3;

        for ( ; ; )
        {

            if ( low == high )
            {
                *left = low;
                return NULL;
            }

            mid = ( low + high + 1 ) / 2;

            if ( t[mid] <= tval )
                low = mid;
            else
                high = mid - 1;
        }
    }
    /*
    CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
    T is just where the user said it might be.
    */
    else;

    return NULL;
}

#endif
