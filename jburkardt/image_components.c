#ifndef __DISABLEDEEP_IMAGECOMPONENTS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4block_components ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4BLOCK_COMPONENTS assigns contiguous nonzero pixels to a common component.
  Discussion:
    On input, the A array contains values of 0 or 1.
    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
    into connected components.
    The pixel A(I,J,K) is "connected" to the pixels:
      A(I-1,J,  K  ),  A(I+1,J,  K  ),
      A(I,  J-1,K  ),  A(I,  J+1,K  ),
      A(I,  J,  K-1),  A(I,  J,  K+1),
    so most pixels have 6 neighbors.
    On output, COMPONENT_NUM reports the number of components of nonzero
    data, and the array C contains the component assignment for
    each nonzero pixel, and is 0 for zero pixels.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int L, M, N, the order of the array.
    Input, int A[L*M*N], the pixel array.
    Output, int C[L*M*N], the component array.
    Output, int I4BLOCK_COMPONENTS, the number of components
    of nonzero data.
*/
{
	static int result = INT_MAX;
	
	const _3dt2pi * const s_data = data;
	const register dim_typ l = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	int * a = s_data->a3;
	int * c = s_data->a4; 
	
    dim_typ b;
    dim_typ c1;
    dim_typ c2;
    int component;
    int component_num;
    dim_typ i, j, k;
    dim_typ north;
    int *p;
    int *q;
    dim_typ up;
    dim_typ west;
    /*
    Initialization.
    */
    for ( k = 0; k < n; ++k )
        for ( j = 0; j < m; ++j )
            for ( i = 0; i < l; ++i )
                c[i+j*l+k*l*m] = 0;
    component_num = 0;
    /*
    P is simply used to store the component labels.  The dimension used
    here is, of course, usually an absurd overestimate.
    */
    p = ( int * ) malloc ( ( l * m * n + 1 ) * sizeof ( int ) );
    for ( i = 0; i <= l * m * n; ++i )
        p[i] = i;
    /*
    "Read" the array one pixel at a time.  If a (nonzero) pixel has a north or
    west neighbor with a label, the current pixel inherits it.
    In case the labels disagree, we need to adjust the P array so we can
    later deal with the fact that the two labels need to be merged.
    */
    for ( i = 0; i < l; i++ )
    {
        for ( j = 0; j < m; j++ )
        {
            for ( k = 0; k < n; k++ )
            {
                north = 0 + (i!=0)*c[i-1+j*l+k*l*m];
                west = 0 + (j!=0)*c[i+(j-1)*l+k*l*m];
                up = 0 + (k!=0)*c[i+j*l+(k-1)*l*m];
                if ( a[i+j*l+k*l*m] != 0 )
                {
                /*
                New component?
                */
                    if ( north == 0 && west == 0 && up == 0 )
                    {
                        ++ component_num;
                        c[i+j*l+k*l*m] = component_num;
                    }
                    /*
                    One predecessor is labeled.
                    */
                    else if ( north != 0 && west == 0 && up == 0 )
                        c[i+j*l+k*l*m] = north;
                    else if ( north == 0 && west != 0 && up == 0 )
                        c[i+j*l+k*l*m] = west;
                    else if ( north == 0 && west == 0 && up != 0 )
                        c[i+j*l+k*l*m] = up;
                    /*
                    Two predecessors are labeled.
                    */
                    else if ( north == 0 && west != 0 && up != 0 )
                    {
                        c[i+j*l+k*l*m] = MIN ( west, up );
                        c1 = MIN ( p[west], p[up] );
                        p[west] = c1;
                        p[up] = c1;
                    }
                    else if ( north != 0 && west == 0 && up != 0 )
                    {
                        c[i+j*l+k*l*m] = MIN ( north, up );
                        c1 = MIN ( p[north], p[up] );
                        p[north] = c1;
                        p[up] = c1;
                    }
                    else if ( north != 0 && west != 0 && up == 0 )
                    {
                        c[i+j*l+k*l*m] = MIN ( north, west );
                        c1 = MIN ( p[north], p[west] );
                        p[north] = c1;
                        p[west] = c1;
                    }
                    /*
                    Three predecessors are labeled.
                    */
                    else if ( north != 0 && west != 0 && up != 0 )
                    {
                        c[i+j*l+k*l*m] = MIN ( north, MIN ( west, up ) );
                        c1 = MIN ( p[north], MIN ( p[west], p[up] ) );
                        p[north] = c1;
                        p[west] = c1;
                        p[up] = c1;
                    }
                }
            }
        }
    }
    /*
    When a component has multiple labels, have the higher labels
    point to the lowest one.
    */
    for ( component = component_num; 1 <= component; --component )
    {
        b = component;
        while ( p[b] != b )
            b = p[b];
        p[component] = b;
    }
    /*
    Locate the minimum label for each component.
    Assign these mininum labels new consecutive indices.
    */
    q = ( int * ) malloc ( ( component_num + 1 ) * sizeof ( int ) );

    for ( j = 0; j <= component_num; ++j )
        q[j] = 0;
    i = 0;
    for ( component = 1; component <= component_num; ++component)
        if ( p[component] == component )
        {
            ++ i;
            q[component] = i;
        }
    component_num = i;
    /*
    Replace the labels by consecutive labels.
    */
    for ( i = 0; i < l; ++i )
        for ( j = 0; j < m; ++j )
            for ( k = 0; k < n; ++k )
                c[i+j*l+k*l*m] = q [ p [ c[i+j*l+k*l*m] ] ];

    free ( p );
    free ( q );

	result = component_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_components ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_COMPONENTS assigns contiguous nonzero pixels to a common component.
  Discussion:
    On input, the A array contains values of 0 or 1.
    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
    into connected components.
    The pixel A(I,J) is "connected" to the pixels A(I-1,J), A(I+1,J),
    A(I,J-1) and A(I,J+1), so most pixels have 4 neighbors.
 (Another choice would be to assume that a pixel was connected
    to the other 8 pixels in the 3x3 block containing it.)
    On output, COMPONENT_NUM reports the number of components of nonzero
    data, and the array C contains the component assignment for
    each nonzero pixel, and is 0 for zero pixels.
  Picture:
    Input A:
      0  2  0  0 17  0  3
      0  0  3  0  1  0  4
      1  0  4  8  8  0  7
      3  0  6 45  0  0  0
      3 17  0  5  9  2  5
    Output:
      COMPONENT_NUM = 4
      C:
      0  1  0  0  2  0  3
      0  0  2  0  2  0  3
      4  0  2  2  2  0  3
      4  0  2  2  0  0  0
      4  4  0  2  2  2  2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the order of the array.
    Input, int A[M*N], the pixel array.
    Output, int C[M*N], the component array.
    Output, int I4MAT_COMPONENTS, the number of components
    of nonzero data.
*/
{
	static int result = INT_MAX;
	
	const _2dt2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	int * c = s_data->a3;
	
    dim_typ b;
    int component;
    int component_num;
    dim_typ i, j;
    int north;
    int *p;
    int *q;
    int west;
    /*
    Initialization.
    */
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            c[i+j*m] = 0;
    component_num = 0;
    /*
    P is simply used to store the component labels.  The dimension used
    here is, of course, usually an absurd overestimate.
    */
    p = ( int * ) malloc ( ( m * n + 1 ) * sizeof ( int ) );

    for ( i = 0; i <= m * n; ++i )
        p[i] = i;
    /*
    "Read" the array one pixel at a time.  If a (nonzero) pixel has a north or
    west neighbor with a label, the current pixel inherits it.
    In case the labels disagree, we need to adjust the P array so we can
    later deal with the fact that the two labels need to be merged.
    */
    for ( i = 0; i < m; ++i )
    {
        for ( j = 0; j < n; ++j )
        {
            north = 0+(i!=0)*c[i-1+j*m];
            west = 0+(j!=0)*c[i+(j-1)*m];
            if ( a[i+j*m] != 0 )
            {
                if ( north == 0 )
                {
                    if ( west == 0 )
                    {
                        ++ component_num;
                        c[i+j*m] = component_num;
                    }
                    else
                        c[i+j*m] = west;
                }
                else if ( north != 0 )
                {
                    if ( west == 0 || west == north )
                    {
                        c[i+j*m] = north;
                    }
                    else
                    {
                        c[i+j*m] = MIN ( north, west );
                        if ( north < west )
                            p[west] = north;
                        else
                            p[north] = west;
                    }
                }
            }
        }
    }
    /*
    When a component has multiple labels, have the higher labels
    point to the lowest one.
    */
    for ( component = component_num; 1 <= component; --component )
    {
        b = component;
        while ( p[b] != b )
            b = p[b];
        p[component] = b;
    }
    /*
    Locate the minimum label for each component.
    Assign these mininum labels new consecutive indices.
    */
    q = ( int * ) malloc ( ( component_num + 1 ) * sizeof ( int ) );

    for ( j = 0; j <= component_num; ++j )
        q[j] = 0;

    i = 0;
    for ( component = 1; component <= component_num; ++component )
        if ( p[component] == component )
        {
            ++ i;
            q[component] = i;
        }

    component_num = i;
    /*
    Replace the labels by consecutive labels.
    */
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            c[i+j*m] = q [ p [ c[i+j*m] ] ];

    free ( p );
    free ( q );

	result = component_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_components ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_COMPONENTS assigns contiguous nonzero pixels to a common component.
  Discussion:
    This calculation is trivial compared to the 2D problem, and is included
    primarily for comparison.
    On input, the A array contains values of 0 or 1.
    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
    into connected components.
    The pixel A(I) is "connected" to the pixels A(I-1) and A(I+1).
    On output, COMPONENT_NUM reports the number of components of nonzero
    data, and the array C contains the component assignment for
    each nonzero pixel, and is 0 for zero pixels.
  Picture:
    Input A:
      0 0 1 2 4 0 0 4 0 0 0 8 9 9 1 2 3 0 0 5 0 1 6 0 0 0 4 0
    Output:
      COMPONENT_NUM = 6
      C:
      0 0 1 1 1 0 0 2 0 0 0 3 3 3 3 3 3 0 0 4 0 5 5 0 0 0 6 0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the vector.
    Input, int A(N), the pixel array.
    Output, int C[N], the component array.
    Output, int I4VEC_COMPONENTS, the number of components
    of nonzero data.
*/
{
	static int result = INT_MAX;
	
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * c = s_data->a2;
	
    int component_num;
    dim_typ j;
    int west;
    /*
    Initialization.
    */
    for ( j = 0; j < n; ++j )
        c[j] = 0;
    component_num = 0;
    /*
    "Read" the array one pixel at a time.  If a (nonzero) pixel has a west
    neighbor with a label, the current pixel inherits it.  Otherwise, we have
    begun a new component.
    */
    west = 0;

    for ( j = 0; j < n; ++j )
    {
        if ( a[j] != 0 )
        {
            if ( west == 0 )
                ++ component_num;
            c[j] = component_num;
        }
        west = c[j];
    }
	
	result = component_num;
    return &result;
}

#endif
