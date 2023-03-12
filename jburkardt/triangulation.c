#ifndef __DISABLEDEEP_TRIANGULATION

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_order3_reference_to_physical ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL maps reference points to physical points.
  Discussion:
    Given the vertices of an order 3 physical triangle and a point
 (XSI,ETA) in the reference triangle, the routine computes the value
    of the corresponding image point (X,Y) in physical space.
    Note that this routine may also be appropriate for an order 6
    triangle, if the mapping between reference and physical space
    is linear.  This implies, in particular, that the sides of the
    image triangle are straight and that the "midside" nodes in the
    physical triangle are halfway along the sides of
    the physical triangle.
  Reference Element T3:
    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 June 2005
  Author:
    John Burkardt
  Parameters:
    Input, double T[2*3], the coordinates of the vertices.
    The vertices are assumed to be the images of (0,0), (1,0) and
 (0,1) respectively.
    Input, int N, the number of points to transform.
    Input, double REF[2*N], points in the reference triangle.
    Output, double PHY[2*N], corresponding points in the
    physical triangle.
*/
{
	const dt3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * ref = s_data->a1;
	ityp * phy = s_data->a2; 
	ityp * t = s_data->a3;
	
	dim_typ i, j;

	#pragma omp parallel for num_threads(2)
	for ( i = 0; i < 2; ++i )
		for ( j = 0; j < n; ++j )
		{
			phy[i+(j<<1)] = t[i] * ( 1.00 - ref[j<<1] - ref[1+(j<<1)] )
			+ t[i+2] *       + ref[j<<1]
			+ t[i+4] *                    + ref[1+(j<<1)];
		}
	
	return NULL;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_order3_physical_to_reference ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE maps physical points to reference points.
  Discussion:
    Given the vertices of an order 3 physical triangle and a point
 (X,Y) in the physical triangle, the routine computes the value
    of the corresponding image point (XSI,ETA) in reference space.
    Note that this routine may also be appropriate for an order 6
    triangle, if the mapping between reference and physical space
    is linear.  This implies, in particular, that the sides of the
    image triangle are straight and that the "midside" nodes in the
    physical triangle are halfway along the sides of
    the physical triangle.
  Reference Element T3:
    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 June 2005
  Author:
    John Burkardt
  Parameters:
    Input, double T[2*3], the X and Y coordinates
    of the vertices.  The vertices are assumed to be the images of
 (0,0), (1,0) and (0,1) respectively.
    Input, int N, the number of points to transform.
    Input, double PHY[2*N], the coordinates of physical points
    to be transformed.
    Output, double REF[2*N], the coordinates of the corresponding
    points in the reference space.
*/
{
	const dt3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * phy = s_data->a1;
	ityp * ref = s_data->a2; 
	ityp * t = s_data->a3;
	
	for (dim_typ j = 0; j < n; ++j )
	{
	
		ref[j<<1] = ( ( t[5] - t[1] ) * ( phy[j<<1] - t[0] )- ( t[4] - t[0] ) * ( phy[1+(j<<1)] - t[1] ) )/ ( ( t[5] - t[1] ) * ( t[2]   - t[0] )- ( t[5] - t[0] ) * ( t[3]   - t[1] ) );
		ref[1+(j<<1)] = ( ( t[2] - t[0] ) * ( phy[1+(j<<1)] - t[1] )- ( t[3] - t[1] ) * ( phy[j<<1] - t[0] ) )/ ( ( t[5] - t[1] ) * ( t[2]   - t[0] )- ( t[4] - t[0] ) * ( t[3]   - t[1] ) );
	}
  return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _alpha_measure ( void * data)
/******************************************************************************/
/*
  Purpose:
    ALPHA_MEASURE determines the triangulated pointset quality measure ALPHA.
  Discusion:
    The ALPHA measure evaluates the uniformity of the shapes of the triangles
    defined by a triangulated pointset.
    We compute the minimum angle among all the triangles in the triangulated
    dataset and divide by the maximum possible value (which, in degrees,
    is 60).  The best possible value is 1, and the worst 0.  A good
    triangulation should have an ALPHA score close to 1.
    The code has been modified to 'allow' 6-node triangulations.
    However, no effort is made to actually process the midside nodes.
    Only information from the vertices is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, real ( kind = 8 ) Z(2,N), the points.
    Input, int TRIANGLE_ORDER, the order of the triangles.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
    the triangulation.
    Output, double *ALPHA_MIN, the minimum value of ALPHA over all
    triangles.
    Output, double *ALPHA_AVE, the value of ALPHA averaged over
    all triangles.
    Output, double *ALPHA_AREA, the value of ALPHA averaged over
    all triangles and weighted by area.
*/
{
	const dtpit2dtpi3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * z = s_data->a1;
	const register dim_typ triangle_order = s_data->a2;
	const register dim_typ triangle_num = s_data->a3;
	int * triangle_node = s_data->a4;
	ityp * alpha_min = s_data->a5;
	ityp * alpha_ave = s_data->a6;
	ityp * alpha_area = s_data->a7;
	
    ityp a_angle;
    dim_typ a_index;
    ityp a_x;
    ityp a_y;
    ityp ab_len;
    ityp alpha;
    ityp area;
    ityp area_total;
    ityp b_angle;
    dim_typ b_index;
    ityp b_x;
    ityp b_y;
    ityp bc_len;
    ityp c_angle;
    dim_typ c_index;
    ityp c_x;
    ityp c_y;
    ityp ca_len;
    dim_typ triangle;
    ityp value;

    *alpha_min = r8_huge;
    *alpha_ave = *alpha_area = area_total = 0.00;

    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        a_index = triangle_node[0+triangle*triangle_order];
        b_index = triangle_node[1+triangle*triangle_order];
        c_index = triangle_node[2+triangle*triangle_order];

        a_x = z[0+((a_index-1)<<1)];
        a_y = z[1+((a_index-1)<<1)];
        b_x = z[0+((b_index-1)<<1)];
        b_y = z[1+((b_index-1)<<1)];
        c_x = z[0+((c_index-1)<<1)];
        c_y = z[1+((c_index-1)<<1)];

        area = 0.50 * abs ( a_x * ( b_y - c_y )
        + b_x * ( c_y - a_y )
        + c_x * ( a_y - b_y ) );

        ab_len = sqrt ( pow ( a_x - b_x, 2 ) + pow ( a_y - b_y, 2 ) );
        bc_len = sqrt ( pow ( b_x - c_x, 2 ) + pow ( b_y - c_y, 2 ) );
        ca_len = sqrt ( pow ( c_x - a_x, 2 ) + pow ( c_y - a_y, 2 ) );
        /*
        Take care of a ridiculous special case.
        */
        if ( ab_len == 0.00 && bc_len == 0.00 && ca_len == 0.00 )
        {
            a_angle = M_2TPI / 3.00;
            b_angle = M_2TPI / 3.00;
            c_angle = M_2TPI / 3.00;
        }
        else
        {
            a_angle = ca_len == 0.00|| ab_len == 0.00 ? M_PI : acos (( ca_len * ca_len + ab_len * ab_len - bc_len * bc_len )/ ( 2.00 * ca_len * ab_len ) );
            b_angle = ab_len == 0.00 || bc_len == 0.00 ? M_PI : acos (( ab_len * ab_len + bc_len * bc_len - ca_len * ca_len )/ ( 2.00 * ab_len * bc_len ) );
            c_angle = bc_len == 0.00 || ca_len == 0.00 ? M_PI : acos (( bc_len * bc_len + ca_len * ca_len - ab_len * ab_len )/ ( 2.00 * bc_len * ca_len ) );
        }
        *alpha_min = MIN ( *alpha_min, a_angle );
        *alpha_min = MIN ( *alpha_min, b_angle );
        *alpha_min = MIN ( *alpha_min, c_angle );

        *alpha_ave += *alpha_min;

        *alpha_area += area * *alpha_min;

        area_total += area;
    }
    *alpha_ave /= ( ityp ) ( triangle_num );
    *alpha_area /= area_total;
    /*
    Normalize angles from [0,M_PI/3] radians into qualities in [0,1].
    */
    *alpha_min *= 3.00 / M_PI;
    *alpha_ave *= 3.00 / M_PI;
    *alpha_area *= 3.00 / M_PI;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _area_measure ( void * data)
/******************************************************************************/
/*
  Purpose:
    AREA_MEASURE determines the area ratio quality measure.
  Discussion:
    This measure computes the area of every triangle, and returns
    the ratio of the minimum to the maximum triangle.  A value of
    1 is "perfect", indicating that all triangles have the same area.
    A value of 0 is the worst possible result.
    The code has been modified to 'allow' 6-node triangulations.
    However, no effort is made to actually process the midside nodes.
    Only information from the vertices is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, double Z[2*N], the points.
    Input, int TRIANGLE_ORDER, the order of the triangles.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
    the triangulation.
    Output, double *AREA_MIN, *AREA_MAX, the minimum and maximum
    areas.
    Output, double *AREA_RATIO, the ratio of the minimum to the
    maximum area.
    Output, double *AREA_AVE, the average area.
    Output, double *AREA_STD, the standard deviation of the areas.
*/
{
	const dtpit2dtpi5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * z = s_data->a1;
	const register dim_typ triangle_order = s_data->a2;
	const register dim_typ triangle_num = s_data->a3;
	int * triangle_node = s_data->a4;
	ityp * area_min = s_data->a5;
	ityp * area_max = s_data->a6;
	ityp * area_ratio = s_data->a7;
	ityp * area_ave = s_data->a8;
	ityp * area_std = s_data->a9;
	
    ityp area;
    dim_typ triangle;
    ityp value;
    ityp x1;
    ityp x2;
    ityp x3;
    ityp y1;
    ityp y2;
    ityp y3;

    *area_max = *area_ave = 0.00;
    *area_min = r8_huge;

    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        x1 = z[0+((triangle_node[0+triangle*triangle_order]-1)<<1)];
        y1 = z[1+((triangle_node[0+triangle*triangle_order]-1)<<1)];
        x2 = z[0+((triangle_node[1+triangle*triangle_order]-1)<<1)];
        y2 = z[1+((triangle_node[1+triangle*triangle_order]-1)<<1)];
        x3 = z[0+((triangle_node[2+triangle*triangle_order]-1)<<1)];
        y3 = z[1+((triangle_node[2+triangle*triangle_order]-1)<<1)];

        area = 0.50 * abs ( x1 * ( y2 - y3 )+ x2 * ( y3 - y1 )+ x3 * ( y1 - y2 ) );

        *area_min = MIN ( *area_min, area );
        *area_max = MAX ( *area_max, area );
        *area_ave += area;
    }

    *area_ave /= ( ityp ) ( triangle_num );
    *area_std = 0.00;
    for ( triangle = 0; triangle < triangle_num; ++triangle)
    {
        x1 = z[0+((triangle_node[0+triangle*triangle_order]-1)<<1)];
        y1 = z[1+((triangle_node[0+triangle*triangle_order]-1)<<1)];
        x2 = z[0+((triangle_node[1+triangle*triangle_order]-1)<<1)];
        y2 = z[1+((triangle_node[1+triangle*triangle_order]-1)<<1)];
        x3 = z[0+((triangle_node[2+triangle*triangle_order]-1)<<1)];
        y3 = z[1+((triangle_node[2+triangle*triangle_order]-1)<<1)];

        area = 0.50 * abs ( x1 * ( y2 - y3 )+ x2 * ( y3 - y1 )+ x3 * ( y1 - y2 ) );

        *area_std += pow ( area - *area_ave, 2 );
    }
    *area_std = sqrt ( *area_std / ( ityp ) ( triangle_num ) );


    *area_ratio = 0.00 + (0.00 < *area_max)*(*area_min / *area_max);
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bandwidth ( void * data)
/******************************************************************************/
/*
  Purpose:
    BANDWIDTH determines the bandwidth associated with the finite element mesh.
  Discussion:
    The quantity computed here is the "geometric" bandwidth determined
    by the finite element mesh alone.
    If a single finite element variable is associated with each node
    of the mesh, and if the nodes and variables are numbered in the
    same way, then the geometric bandwidth is the same as the bandwidth
    of a typical finite element matrix.
    The bandwidth M is defined in terms of the lower and upper bandwidths:
      M = ML + 1 + MU
    where
      ML = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but earlier column,
      MU = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but later column.
    Because the finite element node adjacency relationship is symmetric,
    we are guaranteed that ML = MU.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2006
  Author:
    John Burkardt
  Parameters:
    Input, int ELEMENT_ORDER, the order of the elements.
    Input, int ELEMENT_NUM, the number of elements.
    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.
    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.
    Output, int *M, the bandwidth of the matrix.
*/
{
	const _2dtpi3pdt * const s_data = data;
	const register dim_typ element_order = s_data->a0;
	const register dim_typ element_num = s_data->a1;
	int * element_node = s_data->a2; 
	dim_typ * ml = s_data->a3;
	dim_typ * mu = s_data->a4;
	dim_typ * m = s_data->a5;
	
    dim_typ element;
    dim_typ global_i;
    dim_typ global_j;
    dim_typ local_i;
    dim_typ local_j;

    *ml = *mu = 0;

    for ( element = 0; element < element_num; ++element )
        for ( local_i = 0; local_i < element_order; ++local_i )
        {
            global_i = element_node[local_i+element*element_order];
            for ( local_j = 0; local_j < element_order; ++local_j )
            {
                global_j = element_node[local_j+element*element_order];

                *mu = MAX ( *mu, global_j - global_i );
                *ml = MAX ( *ml, global_i - global_j );
            }
        }

    *m = *ml + 1 + *mu;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _delaunay_swap_test ( void * data)
/******************************************************************************/
/*
  Purpose:
    DELAUNAY_SWAP_TEST performs the Delaunay swap test.
  Discussion:
    The current triangles are formed by nodes [0+2,3) and [0+3,4).
    if a swap is recommended, the new triangles will be [0+2,4) and [1+3,4).
      4     ?     4
     / \         /|\
    1---3  ==>  1 | 3
     \ /         \|/
      2           2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 June 2009
  Author:
    John Burkardt
  Reference:
    Graham Carey,
    Computational Grids:
    Generation, Adaptation and Solution Strategies,
    Taylor and Francis, 1997,
    ISBN13: 978-1560326359,
    LC: QA377.C32.
  Parameters:
    Input, double XY[2*4], the coordinates of four points.
    Output, int SWAP, is TRUE if the triangles [0+2,4) and [1+3,4)
    are to replace triangles [0+2,3) and [0+3,4).
*/
{
	static bool result = 2;
	
	ityp * xy = data;
	
    ityp a;
    ityp b;
    ityp c;
    ityp d;
    ityp x13;
    ityp x14;
    ityp x23;
    ityp x24;
    ityp y13;
    ityp y14;
    ityp y23;
    ityp y24;

    x13 = xy[0] - xy[4];
    x14 = xy[0] - xy[6];
    x23 = xy[2] - xy[4];
    x24 = xy[2] - xy[6];

    y13 = xy[1] - xy[5];
    y14 = xy[1] - xy[7];
    y23 = xy[3] - xy[5];
    y24 = xy[3] - xy[7];

    a = x13 * x23 + y13 * y23;
    b = x24 * y14 - x14 * y24;
    c = x23 * y13 - x13 * y23;
    d = x24 * x14 + y14 * y24;
    /*
    The reference gives two initial tests before the
    main one.  However, there seems to be an error
    in at least one of these tests.  Since they are
    intended to avoid error in borderline cases, but
    instead cause real error in common cases, they are
    omitted for now.

    if ( 0.0 <= a && 0.0 <= d )
    {
    swap = 1;
    }
    else if ( a < d && d < 0.0 )
    {
    swap = 1;
    }
    else if ...
    */

	result = a * b < c * d;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ns_adj_col_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    NS_AJ_COL_SET sets the COL array in a Navier Stokes triangulation.
  Discussion:
    This routine also counts the the value and returns the value of
    ADJ_NUM, the number of Navier-Stokes variable adjacencies, which
    should be identical to the value that would have been computed
    by calling NS_ADJ_COUNT.
    This routine is called to set up the ADJ_COL array, which indicates
    the number of entries needed to store each column in the sparse
    compressed column storage to be used for the adjancency matrix.
    The triangulation is assumed to involve 6-node triangles.
    Variables for the horizontal and vertical velocities are associated
    with every node.  Variables for the pressure are associated only with
    the vertex nodes.
    We are interested in determining the number of nonzero entries in the
    stiffness matrix of the Stokes equations, or the jacobian matrix of
    the Navier Stokes equations.  To this end, we will say, somewhat
    too broadly, that two variables are "adjacent" if their associated
    nodes both occur in some common element.  This adjacency of variables
    I and J is taken to be equivalent to the possible nonzeroness of
    matrix entries A(I,J) and A(J,I).
    A sparse compressed column format is used to store the counts for
    the nonzeroes.  In other words, while the value ADJ_NUM reports the
    number of adjacencies, the vector ADJ_COL is sufficient to allow us
    to properly set up a sparse compressed matrix for the actual storage
    of the sparse matrix, if we desire to proceed.
  Local Node Numbering:
       3
    s  |\
    i  | \
    d  |  \
    e  6   5  side 2
       |    \
    3  |     \
       |      \
       1---4---2
         side 1
  Variable Diagram:
      UVP
       |\
       | \
       |  \
      UV   UV
       |    \
       |     \
       |      \
      UVP--UV--UVP
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 September 2006
  Author:
    John Burkardt
  Parameters
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int VARIABLE_NUM, the number of variables.
    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
    make up each triangle.  The first three nodes are the vertices,
    in counterclockwise order.  The fourth value is the midside
    node between nodes 1 and 2; the fifth and sixth values are
    the other midside nodes in the logical order.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
    a triangle, lists the neighboring triangle, or -1 if there is
    no neighbor.
    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
    vertical velocity and pressure variables associated with a node,
    or -1 if no such variable is associated with the node.
    Output, int ADJ_COL[VARIABLE_NUM+1].  Information about variable J
    is stored in entries ADJ_COL[J] through ADJ_COL[J+1]-1 of ADJ.
    Output, int NS_ADJ_COL_SET, the number of Navier Stokes variable
    adjacencies.
*/
{
	static int result = INT_MAX;
	
	const _3dt6pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	const register dim_typ variable_num = s_data->a2;
	int * triangle_node = s_data->a3;
	int * triangle_neighbor = s_data->a4;
	int * node_u_variable = s_data->a5;
	int * node_v_variable = s_data->a6;
	int * node_p_variable = s_data->a7;
	int * adj_col = s_data->a8;
	
    dim_typ adj_num;
    int n1;
    int n2;
    int n3;
    int n4;
    int n5;
    int n6;
    dim_typ node;
    int p1;
    int p2;
    int p3;
    dim_typ triangle;
    dim_typ triangle_order = 6;
    dim_typ triangle2;
    int u1;
    int u2;
    int u3;
    int u4;
    int u5;
    int u6;
    int v1;
    int v2;
    int v3;
    int v4;
    int v5;
    int v6;
    int variable;

    adj_num = 0;
    /*
    Set every variable to be adjacent to itself.
    */
    for ( variable = 0; variable < variable_num; ++variable )
        adj_col[variable] = 1;
    /*
    Set every variable to be adjacent to the other variables associated with
    that node.

    U <=> V
    U <=> P (if there is a P variable)
    V <=> P (if there is a P variable)
    */
    for ( node = 0; node < node_num; ++node )
    {
        u1 = node_u_variable[node] - 1;
        v1 = node_v_variable[node] - 1;
        p1 = node_p_variable[node] - 1;

        ++ adj_col[u1];
        ++ adj_col[v1];

        if ( 0 <= p1 )
        {
            ++ adj_col[u1];
            adj_col[v1] += 1;
            adj_col[p1] += 2;
        }
    }
    /*
    Examine each triangle.
    */
    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        n1 = triangle_node[0+triangle*triangle_order] - 1;
        n2 = triangle_node[1+triangle*triangle_order] - 1;
        n3 = triangle_node[2+triangle*triangle_order] - 1;
        n4 = triangle_node[3+triangle*triangle_order] - 1;
        n5 = triangle_node[4+triangle*triangle_order] - 1;
        n6 = triangle_node[5+triangle*triangle_order] - 1;

        u1 = node_u_variable[n1] - 1;
        v1 = node_v_variable[n1] - 1;
        p1 = node_p_variable[n1] - 1;

        u2 = node_u_variable[n2] - 1;
        v2 = node_v_variable[n2] - 1;
        p2 = node_p_variable[n2] - 1;

        u3 = node_u_variable[n3] - 1;
        v3 = node_v_variable[n3] - 1;
        p3 = node_p_variable[n3] - 1;

        u4 = node_u_variable[n4] - 1;
        v4 = node_v_variable[n4] - 1;

        u5 = node_u_variable[n5] - 1;
        v5 = node_v_variable[n5] - 1;

        u6 = node_u_variable[n6] - 1;
        v6 = node_v_variable[n6] - 1;
        /*
        For sure, we add the new adjacencies:

        U5 V5 <=> U1 V1 P1
        U6 V6 <=> U2 V2 P2
        U4 V4 <=> U3 V3 P3
        U5 V5 <=> U4 V4
        U6 V6 <=> U4 V4
        U6 V6 <=> U5 V5
        */
        adj_col[u1] += 2;
        adj_col[v1] += 2;
        adj_col[p1] += 2;

        adj_col[u2] += 2;
        adj_col[v2] += 2;
        adj_col[p2] += 2;

        adj_col[u3] += 2;
        adj_col[v3] += 2;
        adj_col[p3] += 2;

        adj_col[u4] += 7;
        adj_col[v4] += 7;

        adj_col[u5] += 7;
        adj_col[v5] += 7;

        adj_col[u6] += 7;
        adj_col[v6] += 7;
        /*
        Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
        that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
        or if this triangle is the first of the pair in which the edge
        occurs (TRIANGLE < TRIANGLE2).

        Maybe add

        U1 V1 P1 <=> U2 V2 P2
        U1 V1 P1 <=> U4 V4
        U2 V2 P2 <=> U4 V4
        */
        triangle2 = triangle_neighbor[0+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj_col[u1] += 5;
            adj_col[v1] += 5;
            adj_col[p1] += 5;

            adj_col[u2] += 5;
            adj_col[v2] += 5;
            adj_col[p2] += 5;

            adj_col[u4] += 6;
            adj_col[v4] += 6;
        }
        /*
        Maybe add

        U2 V2 P2 <=> U3 V3 P3
        U2 V2 P2 <=> U5 V5
        U3 V3 P3 <=> U5 V5
        */
        triangle2 = triangle_neighbor[1+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj_col[u2] += 5;
            adj_col[v2] += 5;
            adj_col[p2] += 5;

            adj_col[u3] += 5;
            adj_col[v3] += 5;
            adj_col[p3] += 5;

            adj_col[u5] += 6;
            adj_col[v5] += 6;
        }
        /*
        Maybe add

        U1 V1 P1 <=> U3 V3 P3
        U1 V1 P1 <=> U6 V6
        U3 V3 P3 <=> U6 V6
        */
        triangle2 = triangle_neighbor[2+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj_col[u1] += 5;
            adj_col[v1] += 5;
            adj_col[p1] += 5;

            adj_col[u3] += + 5;
            adj_col[v3] += + 5;
            adj_col[p3] += + 5;

            adj_col[u6] += + 6;
            adj_col[v6] += 6;
        }

    }
    /*
    We used ADJ_COL to count the number of entries in each column.
    Convert it to pointers into the ADJ array.
    */
    for ( variable = variable_num; 0 < variable; --variable)
        adj_col[variable] = adj_col[variable-1];

    adj_col[0] = 1;
    for ( variable = 1; variable <= variable_num; ++variable )
        adj_col[variable] += adj_col[variable-1];

    adj_num = adj_col[variable_num] - 1;

	result = adj_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ns_adj_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    NS_ADJ_COUNT counts adjacencies in a Navier Stokes triangulation.
  Discussion:
    This routine is called to count the adjacencies, so that the
    appropriate amount of memory can be set aside for storage when
    the adjacency structure is created.
    The value of ADJ_NUM computed and returned by this routine should
    be identical to the value computed by NS_ADJ_COL_SET.
    The triangulation is assumed to involve 6-node triangles.
    Variables for the horizontal and vertical velocities are associated
    with every node.  Variables for the pressure are associated only with
    the vertex nodes.
    We are interested in determining the number of nonzero entries in the
    stiffness matrix of the Stokes equations, or the jacobian matrix of
    the Navier Stokes equations.  To this end, we will say, somewhat
    too broadly, that two variables are "adjacent" if their associated
    nodes both occur in some common element.  This adjacency of variables
    I and J is taken to be equivalent to the possible nonzeroness of
    matrix entries A(I,J) and A(J,I).
    A sparse compressed column format is used to store the counts for
    the nonzeroes.  In other words, while the value ADJ_NUM reports the
    number of adjacencies, the vector ADJ_COL is sufficient to allow us
    to properly set up a sparse compressed matrix for the actual storage
    of the sparse matrix, if we desire to proceed.
  Local Node Numbering:
       3
    s  |\
    i  | \
    d  |  \
    e  6   5  side 2
       |    \
    3  |     \
       |      \
       1---4---2
         side 1
  Variable Diagram:
      UVP
       |\
       | \
       |  \
      UV   UV
       |    \
       |     \
       |      \
      UVP--UV--UVP
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 September 2006
  Author:
    John Burkardt
  Parameters
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int VARIABLE_NUM, the number of variables.
    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
    make up each triangle.  The first three nodes are the vertices,
    in counterclockwise order.  The fourth value is the midside
    node between nodes 1 and 2; the fifth and sixth values are
    the other midside nodes in the logical order.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
    a triangle, lists the neighboring triangle, or -1 if there is
    no neighbor.
    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
    vertical velocity and pressure variables associated with a node,
    or -1 if no such variable is associated with the node.
    Output, int NS_ADJ_COUNT, the value of ADJ_NUM, the number of
    Navier Stokes variable adjacencies.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _3dt6pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	const register dim_typ variable_num = s_data->a2;
	int * triangle_node = s_data->a3;
	int * triangle_neighbor = s_data->a4;
	int * node_u_variable = s_data->a5;
	int * node_v_variable = s_data->a6;
	int * node_p_variable = s_data->a7;
	int * adj_col = s_data->a8;
	
    dim_typ adj_num = 0;
    dim_typ node;
    int p1;
    dim_typ triangle;
    dim_typ triangle_order = 6;
    dim_typ triangle2;
    dim_typ variable;

    /*
    Set every variable to be adjacent to itself.
    */
    adj_num = variable_num;
    /*
    Set every variable to be adjacent to the other variables associated with
    that node.

    U <=> V
    U <=> P (if there is a P variable)
    V <=> P (if there is a P variable)
    */
    for ( node = 0; node < node_num; ++node )
    {
        adj_num += 2;
        p1 = node_p_variable[node] - 1;

        if ( 0 <= p1 )
            adj_num += 4;
    }
    /*
    Examine each triangle.
    */
    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        /*
        For sure, we add the new adjacencies:

        U5 V5 <=> U1 V1 P1
        U6 V6 <=> U2 V2 P2
        U4 V4 <=> U3 V3 P3
        U5 V5 <=> U4 V4
        U6 V6 <=> U4 V4
        U6 V6 <=> U5 V5
        */
        adj_num += 60;
        /*
        Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
        that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
        or if this triangle is the first of the pair in which the edge
        occurs (TRIANGLE < TRIANGLE2).

        Maybe add

        U1 V1 P1 <=> U2 V2 P2
        U1 V1 P1 <=> U4 V4
        U2 V2 P2 <=> U4 V4
        */
        triangle2 = triangle_neighbor[0+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
            adj_num += 42;
        /*
        Maybe add

        U2 V2 P2 <=> U3 V3 P3
        U2 V2 P2 <=> U5 V5
        U3 V3 P3 <=> U5 V5
        */
        triangle2 = triangle_neighbor[1+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
            adj_num += 42;
        /*
        Maybe add

        U1 V1 P1 <=> U3 V3 P3
        U1 V1 P1 <=> U6 V6
        U3 V3 P3 <=> U6 V6
        */
        triangle2 = triangle_neighbor[2+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
            adj_num += 42;

    }

	result = adj_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _ns_adj_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    NS_ADJ_INSERT inserts an adjacency into a compressed column adjacency matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2006
  Author:
    John Burkardt
  Parameters:
    Input, int V1, V2, the indices of two items which are adjacent.
    Input, int VARIABLE_NUM, the number of items.
    Input, int ADJ_NUM, the number of entries available in ADJ_ROW.
    Input/output, int ADJ_COL_FREE[VARIABLE_NUM], contains the next free
    location in which an entry for a given column can be stored.  On output,
    two pointers have been updated.
    Input/output, int ADJ_ROW[ADJ_NUM], the row indices of the Navier Stokes
    variable adjacency matrix.  On output, two new entries have been added.
*/
{
	const _2dtpi2dtpi * const s_data = data;
	
	const register dim_typ v1 = s_data->a0;
	const register dim_typ v2 = s_data->a1;
	int * adj_col_free = s_data->a2;
	const register dim_typ variable_num = s_data->a3;
	const register dim_typ adj_num = s_data->a4;
	int * adj_row = s_data->a5;
	
    adj_row[adj_col_free[v1-1]-1] = v2;
    ++ adj_col_free[v1-1];

    if ( v1 == v2 )
        return NULL;

    adj_row[adj_col_free[v2-1]-1] = v1;
    ++ adj_col_free[v2-1];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ns_adj_row_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    NS_ADJ_ROW_SET sets the Navier Stokes sparse compressed column row indices.
  Discussion:
    After NS_ADJ_COUNT has been called to count ADJ_NUM, the number of
    variable adjacencies and to set up ADJ_COL, the compressed column pointer,
    this routine can be called to assign values to ADJ_ROW, the row
    indices of the sparse compressed column adjacency matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 September 2006
  Author:
    John Burkardt
  Parameters
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int VARIABLE_NUM, the number of variables.
    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
    make up each triangle.  The first three nodes are the vertices,
    in counterclockwise order.  The fourth value is the midside
    node between nodes 1 and 2; the fifth and sixth values are
    the other midside nodes in the logical order.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
    a triangle, lists the neighboring triangle, or -1 if there is
    no neighbor.
    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
    vertical velocity and pressure variables associated with a node,
    or -1 if no such variable is associated with the node.
    Input, int ADJ_NUM, the number of Navier Stokes variable adjacencies.
    Input, int ADJ_COL[VARIABLE_NUM+1].  Information about variable J
    is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
    Output, int ADJ_ROW[ADJ_NUM], the row indices of the Navier Stokes
    variable adjacency matrix.
  Local Parameters:
    Local, int ADJ_COL_FREE[VARIABLE_NUM], for each column,
    the location in ADJ_ROW which can store the next row index.
*/
{
	const _3i5pidt2pi * const s_data = data;
	int node_num = s_data->a0;
	int triangle_num = s_data->a1;
	int variable_num = s_data->a2;
	int * triangle_node = s_data->a3;
	int * triangle_neighbor = s_data->a4;
	int * node_u_variable = s_data->a5;
	int * node_v_variable = s_data->a6;
	int * node_p_variable = s_data->a7;
	const register dim_typ adj_num = s_data->a8;
	int * adj_col = s_data->a9;
	int * adj_row = s_data->a10;
	
    int *adj_col_free;
    int k1;
    int k2;
    int n1;
    int n2;
    int n3;
    int n4;
    int n5;
    int n6;
    dim_typ node;
    int p1;
    int p2;
    int p3;
    dim_typ triangle;
    dim_typ triangle_order = 6;
    dim_typ triangle2;
    int u1;
    int u2;
    int u3;
    int u4;
    int u5;
    int u6;
    dim_typ v;
    int v1;
    int v2;
    int v3;
    int v4;
    int v5;
    int v6;

    for ( v = 0; v < adj_num; ++v )
        adj_row[v] = -1;

    adj_col_free = ( int * ) malloc ( variable_num * sizeof ( int ) );

    for ( v = 0; v < variable_num; ++v )
        adj_col_free[v] = adj_col[v];
    /*
    Set every variable to be adjacent to itself.
    Here, we have to be careful to start at index 1.
    */
    for ( v = 1; v <= variable_num; ++v )
        ns_adj_insert ( v, v, variable_num, adj_num, adj_col_free, adj_row );
    /*
    Set every variable to be adjacent to the other variables associated with
    that node.

    U <=> V
    U <=> P (if there is a P variable)
    V <=> P (if there is a P variable)
    */
    for ( node = 0; node < node_num; ++node )
    {
        u1 = node_u_variable[node];
        v1 = node_v_variable[node];
        p1 = node_p_variable[node];

        ns_adj_insert ( u1, v1, variable_num, adj_num, adj_col_free, adj_row );

        if ( 0 < p1 )
        {
            ns_adj_insert ( u1, p1, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, p1, variable_num, adj_num, adj_col_free, adj_row );
        }
    }
    /*
    Examine each triangle.
    */
    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        n1 = triangle_node[0+triangle*6];
        n2 = triangle_node[1+triangle*6];
        n3 = triangle_node[2+triangle*6];
        n4 = triangle_node[3+triangle*6];
        n5 = triangle_node[4+triangle*6];
        n6 = triangle_node[5+triangle*6];

        u1 = node_u_variable[n1-1];
        v1 = node_v_variable[n1-1];
        p1 = node_p_variable[n1-1];

        u2 = node_u_variable[n2-1];
        v2 = node_v_variable[n2-1];
        p2 = node_p_variable[n2-1];

        u3 = node_u_variable[n3-1];
        v3 = node_v_variable[n3-1];
        p3 = node_p_variable[n3-1];

        u4 = node_u_variable[n4-1];
        v4 = node_v_variable[n4-1];

        u5 = node_u_variable[n5-1];
        v5 = node_v_variable[n5-1];

        u6 = node_u_variable[n6-1];
        v6 = node_v_variable[n6-1];
        /*
        For sure, we add the new adjacencies:

        U5 V5 <=> U1 V1 P1
        U6 V6 <=> U2 V2 P2
        U4 V4 <=> U3 V3 P3
        U5 V5 <=> U4 V4
        U6 V6 <=> U4 V4
        U6 V6 <=> U5 V5
        */
        ns_adj_insert ( u5, u1, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u5, v1, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u5, p1, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v5, u1, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v5, v1, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v5, p1, variable_num, adj_num, adj_col_free, adj_row );

        ns_adj_insert ( u6, u2, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u6, v2, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u6, p2, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v6, u2, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v6, v2, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v6, p2, variable_num, adj_num, adj_col_free, adj_row );

        ns_adj_insert ( u4, u3, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u4, v3, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u4, p3, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v4, u3, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v4, v3, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v4, p3, variable_num, adj_num, adj_col_free, adj_row );

        ns_adj_insert ( u5, u4, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u5, v4, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v5, u4, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v5, v4, variable_num, adj_num, adj_col_free, adj_row );

        ns_adj_insert ( u6, u4, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u6, v4, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v6, u4, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v6, v4, variable_num, adj_num, adj_col_free, adj_row );

        ns_adj_insert ( u6, u5, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( u6, v5, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v6, u5, variable_num, adj_num, adj_col_free, adj_row );
        ns_adj_insert ( v6, v5, variable_num, adj_num, adj_col_free, adj_row );
        /*
        Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
        that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
        or if this triangle is the first of the pair in which the edge
        occurs (TRIANGLE < TRIANGLE2).

        Maybe add

        U1 V1 P1 <=> U2 V2 P2
        U1 V1 P1 <=> U4 V4
        U2 V2 P2 <=> U4 V4
        */
        triangle2 = triangle_neighbor[0+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ns_adj_insert ( u1, u2, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u1, v2, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u1, p2, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, u2, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, v2, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, p2, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, u2, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, v2, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, p2, variable_num, adj_num, adj_col_free, adj_row );

            ns_adj_insert ( u1, u4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u1, v4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, u4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, v4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, u4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, v4, variable_num, adj_num, adj_col_free, adj_row );

            ns_adj_insert ( u2, u4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u2, v4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v2, u4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v2, v4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p2, u4, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p2, v4, variable_num, adj_num, adj_col_free, adj_row );
        }
        /*
        Maybe add

        U2 V2 P2 <=> U3 V3 P3
        U2 V2 P2 <=> U5 V5
        U3 V3 P3 <=> U5 V5
        */
        triangle2 = triangle_neighbor[1+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ns_adj_insert ( u2, u3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u2, v3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u2, p3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v2, u3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v2, v3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v2, p3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p2, u3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p2, v3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p2, p3, variable_num, adj_num, adj_col_free, adj_row );

            ns_adj_insert ( u2, u5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u2, v5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v2, u5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v2, v5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p2, u5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p2, v5, variable_num, adj_num, adj_col_free, adj_row );

            ns_adj_insert ( u3, u5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u3, v5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v3, u5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v3, v5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p3, u5, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p3, v5, variable_num, adj_num, adj_col_free, adj_row );
        }
        /*
        Maybe add

        U1 V1 P1 <=> U3 V3 P3
        U1 V1 P1 <=> U6 V6
        U3 V3 P3 <=> U6 V6
        */
        triangle2 = triangle_neighbor[2+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ns_adj_insert ( u1, u3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u1, v3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u1, p3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, u3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, v3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, p3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, u3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, v3, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, p3, variable_num, adj_num, adj_col_free, adj_row );

            ns_adj_insert ( u1, u6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u1, v6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, u6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v1, v6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, u6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p1, v6, variable_num, adj_num, adj_col_free, adj_row );

            ns_adj_insert ( u3, u6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( u3, v6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v3, u6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( v3, v6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p3, u6, variable_num, adj_num, adj_col_free, adj_row );
            ns_adj_insert ( p3, v6, variable_num, adj_num, adj_col_free, adj_row );
        }
    }
    /*
    Ascending sort the entries for each variable.
    */
    for ( v = 0; v < variable_num; ++v )
    {
        k1 = adj_col[v];
        k2 = adj_col[v+1]-1;

        i4vec_sort_heap_a ( k2+1-k1, adj_row+k1-1 );
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _q_measure ( void * data)
/******************************************************************************/
/*
  Purpose:
    Q_MEASURE determines the triangulated pointset quality measure Q.
  Discussion:
    The Q measure evaluates the uniformity of the shapes of the triangles
    defined by a triangulated pointset.
    For a single triangle T, the value of Q(T) is defined as follows:
      TAU_IN = radius of the inscribed circle,
      TAU_OUT = radius of the circumscribed circle,
      Q(T) = 2 * TAU_IN / TAU_OUT
        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
    where A, B and C are the lengths of the sides of the triangle T.
    The Q measure computes the value of Q(T) for every triangle T in the
    triangulation, and then computes the minimum of this
    set of values:
      Q_MEASURE = MIN ( all T in triangulation ) Q(T)
    In an ideally regular mesh, all triangles would have the same
    equilateral shape, for which Q = 1.  A good mesh would have
    0.5 < Q.
    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
    triangles.  Generally, a maximal triangulation is expected, namely,
    a triangulation whose image is a planar graph, but for which the
    addition of any new triangle would mean the graph was no longer planar.
    A Delaunay triangulation is a maximal triangulation which maximizes
    the minimum angle that occurs in any triangle.
    The code has been modified to 'allow' 6-node triangulations.
    However, no effort is made to actually process the midside nodes.
    Only information from the vertices is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 June 2009
  Author:
    John Burkardt
  Reference:
    Max Gunzburger and John Burkardt,
    Uniformity Measures for Point Samples in Hypercubes.
    Per-Olof Persson and Gilbert Strang,
    A Simple Mesh Generator in MATLAB,
    SIAM Review,
    Volume 46, Number 2, pages 329-345, June 2004.
  Parameters:
    Input, int N, the number of points.
    Input, double Z[2*N], the points.
    Input, int TRIANGLE_ORDER, the order of the triangles.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
    the triangulation.
    Output, double *Q_MIN, *Q_MAX, the minimum and maximum values
    of Q over all triangles.
    Output, double *Q_AVE, the average value of Q.
    Output, double *Q_AREA, the average value of Q, weighted by
    the area of each triangle.
*/
{
	const dtpit2dtpi4pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * z = s_data->a1;
	const register dim_typ triangle_order = s_data->a2;
	const register dim_typ triangle_num = s_data->a3;
	int * triangle_node = s_data->a4;
	ityp * q_min = s_data->a5;
	ityp * q_max = s_data->a6;
	ityp * q_ave = s_data->a7;
	ityp * q_area = s_data->a8;
	
    dim_typ a_index;
    ityp ab_length;
    ityp area;
    ityp area_total;
    dim_typ b_index;
    ityp bc_length;
    dim_typ c_index;
    ityp ca_length;
    ityp q;
    dim_typ triangle;
    ityp x1;
    ityp x2;
    ityp x3;
    ityp y1;
    ityp y2;
    ityp y3;

    *q_min = r8_huge;
    *q_max = - r8_huge;
    *q_ave = *q_area = area_total = 0.00;

    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        a_index = triangle_node[0+triangle*triangle_order];
        b_index = triangle_node[1+triangle*triangle_order];
        c_index = triangle_node[2+triangle*triangle_order];

        ab_length = sqrt (pow ( z[0+((a_index-1)<<1)] - z[0+((b_index-1)<<1)], 2 )+ pow ( z[1+((a_index-1)<<1)] - z[1+((b_index-1)<<1)], 2 ) );
        bc_length = sqrt (pow ( z[0+((b_index-1)<<1)] - z[0+((c_index-1)<<1)], 2 )+ pow ( z[1+((b_index-1)<<1)] - z[1+((c_index-1)<<1)], 2 ) );
        ca_length = sqrt (pow ( z[0+((c_index-1)<<1)] - z[0+((a_index-1)<<1)], 2 )+ pow ( z[1+((c_index-1)<<1)] - z[1+((a_index-1)<<1)], 2 ) );

        q = ( bc_length + ca_length - ab_length )* ( ca_length + ab_length - bc_length )* ( ab_length + bc_length - ca_length )/ ( ab_length * bc_length * ca_length );

        x1 = z[0+((triangle_node[0+triangle*triangle_order]-1)<<1)];
        y1 = z[1+((triangle_node[0+triangle*triangle_order]-1)<<1)];
        x2 = z[0+((triangle_node[1+triangle*triangle_order]-1)<<1)];
        y2 = z[1+((triangle_node[1+triangle*triangle_order]-1)<<1)];
        x3 = z[0+((triangle_node[2+triangle*triangle_order]-1)<<1)];
        y3 = z[1+((triangle_node[2+triangle*triangle_order]-1)<<1)];

        area = 0.50 * abs ( x1 * ( y2 - y3 )+ x2 * ( y3 - y1 )+ x3 * ( y1 - y2 ) );

        *q_min = MIN ( *q_min, q );
        *q_max = MAX ( *q_max, q );
        *q_ave += q;
        *q_area += q * area;

        area_total += area;
    }

    *q_ave /= ( ityp ) ( triangle_num );

    *q_area = 0.00 + ( 0.0 < area_total )*(*q_area / area_total);

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _quad_convex_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    QUAD_CONVEX_RANDOM returns a random convex quadrilateral.
  Description:
    The quadrilateral is constrained in that the vertices must all lie
    with the unit square.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 June 2009
  Author:
    John Burkardt
  Parameters:
    Input/output, int *SEED, a seed for the random number
    generator.
    Output, double XY[2*NODE_NUM], the coordinates of the
    nodes of the quadrilateral, given in counterclockwise order.
*/
{
	const pipit * const s_data = data;
	int * seed = s_data->a0;
	ityp * xy = s_data->a1;
	
    int hull[4];
    dim_typ hull_num;
    dim_typ i, j;
    ityp xy_random[8];

    for ( ; ; )
    {
        /*
        Generate 4 random points.
        */
        r8mat_uniform_01 ( 2, 4, seed, xy_random );
        /*
        Determine the convex hull.
        */
        points_hull_2d ( 4, xy_random, &hull_num, hull );
        /*
        If HULL_NUM < 4, then our convex hull is a triangle.
        Try again.
        */
        if ( hull_num == 4 )
            break;
    }
    /*
    Make an ordered copy of the random points.
    */
    #pragma omp parallel for num_threads(2)
    for ( j = 0; j < 4; ++j )
        #pragma omp parallel for num_threads(2)
        for ( i = 0; i < 2; ++i )
            xy[i+(j<<1)] = xy_random[i+((hull[j]-1)<<1)];
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _points_delaunay_naive_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    POINTS_DELAUNAY_NAIVE_2D computes the Delaunay triangulation in 2D.
  Discussion:
    A naive and inefficient (but extremely simple) method is used.
    This routine is only suitable as a demonstration code for small
    problems.  Its running time is of order NODE_NUM^4.  Much faster
    algorithms are available.
    Given a set of nodes in the plane, a triangulation is a set of
    triples of distinct nodes, forming triangles, so that every
    point with the convex hull of the set of  nodes is either one
    of the nodes, or lies on an edge of one or more triangles,
    or lies within exactly one triangle.
    The number of nodes must be at least 3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 June 2005
  Author:
    John Burkardt
  Reference:
    Joseph ORourke,
    Computational Geometry,
    Cambridge University Press,
    Second Edition, 1998, page 187.
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Output, int *TRIANGLE_NUM, the number of triangles.
    Output, int POINTS_DELAUNAY_NAIVE_2D[3*TRIANGLE_NUM], the indices of the
    nodes making each triangle.
*/
{
	const dtpitpdt * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	ityp * node_xy = s_data->a1;
	dim_typ * triangle_num = s_data->a2;
	
    dim_typ count = 0;
    bool flag;
    dim_typ i, j, k, m;
    dim_typ pass;
    int *tri;
    ityp xn;
    ityp yn;
    ityp zn;
    ityp *z = ( ityp * ) malloc ( node_num * sizeof ( ityp ) );

    for ( i = 0; i < node_num; ++i )
        z[i] = node_xy[0+(i<<1)] * node_xy[0+(i<<1)] + node_xy[1+(i<<1)] * node_xy[1+(i<<1)];
    /*
    First pass counts triangles,
    Second pass allocates triangles and sets them.
    */
    for ( pass = 1; pass <= 2; ++pass )
    {
        if ( pass == 2 )
            tri = ( int * ) malloc ( 3 * count * sizeof ( int ) );
        count = 0;
        /*
        For each triple (I,J,K):
        */
        for ( i = 0; i < node_num - 2; ++i )
        {
            for ( j = i+1; j < node_num; ++j )
            {
                for ( k = i+1; k < node_num; ++k )
                {
                    if ( j != k )
                    {
                        xn = ( node_xy[1+(j<<1)] - node_xy[1+(i<<1)] ) * ( z[k] - z[i] )- ( node_xy[1+(k<<1)] - node_xy[1+(i<<1)] ) * ( z[j] - z[i] );
                        yn = ( node_xy[0+(k<<1)] - node_xy[0+(i<<1)] ) * ( z[j] - z[i] )- ( node_xy[0+(j<<1)] - node_xy[0+(i<<1)] ) * ( z[k] - z[i] );
                        zn = ( node_xy[0+(j<<1)] - node_xy[0+(i<<1)] )* ( node_xy[1+(k<<1)] - node_xy[1+(i<<1)] )- ( node_xy[0+(k<<1)] - node_xy[0+(i<<1)] )* ( node_xy[1+(j<<1)] - node_xy[1+(i<<1)] );

                        flag = ( zn < 0 );

                        if ( flag )
                            for ( m = 0; m < node_num; ++m )
                                flag = flag && ( ( node_xy[0+(m<<1)] - node_xy[0+(i<<1)] ) * xn+ ( node_xy[1+(m<<1)] - node_xy[1+(i<<1)] ) * yn+ ( z[m] - z[i] ) * zn <= 0 );

                        if ( flag )
                        {
                            if ( pass == 2 )
                            {
                                tri[0+count*3] = i + 1;
                                tri[1+count*3] = j + 1;
                                tri[2+count*3] = k + 1;
                            }
                            ++ count;
                        }
                    }
                }
            }
        }
    }

    *triangle_num = count;
    free ( z );

    return tri;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_AREA computes the area of a triangulation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 December 2011
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Input, int ELEMENT_ORDER, the order of the triangles.
    Input, int ELEMENT_NUM, the number of triangles.
    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM],
    the nodes making up each triangle.
    Output, double TRIANGULATION_AREA, the area.
*/
{
	static ityp result = MAX_VAL;
	
	const _3dtpipit * const s_data = data;
	
	const register dim_typ node_num = s_data->a0;
	const register dim_typ element_order = s_data->a1;
	const register dim_typ element_num = s_data->a2;
	int * element_node = s_data->a3;
	ityp * node_xy = s_data->a4;
 
    dim_typ element;
    ityp element_area;
    ityp element_xy[6];
    dim_typ j, nj;
    ityp value = 0.00;

    for ( element = 0; element < element_num; ++ element )
    {
        #pragma omp parallel for num_threads(3)
        for ( j = 0; j < 3; ++j )
        {
            nj = element_node[j+element*element_order];
            element_xy[(j<<1)] = node_xy[nj<<1];
            element_xy[1+(j<<1)] = node_xy[1+(nj<<1)];
        }

        element_area = 0.50 * (element_xy[0] * ( element_xy[3] - element_xy[5] ) +element_xy[2] * ( element_xy[5] - element_xy[1] ) +element_xy[4] * ( element_xy[1] - element_xy[3] ) );
        value += element_area;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_areas ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_AREAS computes triangle and triangulation areas.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2009
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes in the
    triangulation.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Input, int TRIANGLE_ORDER, the order of triangles in
    the triangulation.
    Input, int TRIANGLE_NUM, the number of triangles in
    the triangulation.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
    the nodes making up each triangle.
    Output, double TRIANGLE_AREA[TRIANGLE_NUM], the area of
    the triangles.
    Output, double TRIANGULATION_AREAS, the area of the triangulation.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit2dtpipit * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	ityp * node_xy = s_data->a1;
	const register dim_typ triangle_order = s_data->a2;
	const register dim_typ triangle_num = s_data->a3;
	int * triangle_node = s_data->a4;
	ityp * triangle_area = s_data->a5;
	
    dim_typ j, nj;
    ityp t_xy[2*3];
    dim_typ triangle;
    ityp triangulation_area = 0.00;

    for ( triangle = 0; triangle < triangle_num; ++ triangle )
    {
        #pragma omp parallel for num_threads(3)
        for ( j = 0; j < 3; ++j )
        {
                nj = triangle_node[j+triangle*triangle_order];
                t_xy[0+(j<<1)] = node_xy[0+(nj<<1)];
                t_xy[1+(j<<1)] = node_xy[1+(nj<<1)];
        }

        triangle_area[triangle] = 0.50 * (t_xy[0] * ( t_xy[3] - t_xy[5] ) +t_xy[2] * ( t_xy[5] - t_xy[1] ) +t_xy[4] * ( t_xy[1] - t_xy[3] ) );
        triangulation_area += triangle_area[triangle];
    }

	result = triangulation_area;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_delaunay_discrepancy_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_DELAUNAY_DISCREPANCY_COMPUTE reports if a triangulation is Delaunay.
  Discussion:
    A (maximal) triangulation is Delaunay if and only if it is locally
    Delaunay.
    A triangulation is Delaunay if the minimum angle over all triangles
    in the triangulation is maximized.  That is, there is no other
    triangulation of the points which has a larger minimum angle.
    A triangulation is locally Delaunay if, for every pair of triangles
    that share an edge, the minimum angle in the two triangles is larger
    than the minimum angle in the two triangles formed by removing the
    common edge and joining the two opposing vertices.
    This function examines the question of whether a given triangulation
    is locally Delaunay.  It does this by looking at every pair of
    neighboring triangles and comparing the minimum angle attained
    for the current triangle pair and the alternative triangle pair.
    Let A(i,j) be the minimum angle formed by triangles T(i) and T(j),
    which are two triangles in the triangulation which share a common edge.
    Let B(I,J) be the minimum angle formed by triangles S(i) and S(j),
    where S(i) and S(j) are formed by removing the common edge of T(i)
    and T(j), and joining the opposing vertices.
    Then the triangulation is Delaunay if B(i,j) <= A(i,j) for every
    pair of neighbors T(i) and T(j).
    If A(i,j) < B(i,j) for at least one pair of neighbors, the triangulation
    is not a Delaunay triangulation.
    This program returns VALUE = MIN ( A(i,j) - B(i,j) ) over all
    triangle neighbors.  VALUE is scaled to be in degrees, for
    comprehensibility.  If VALUE is negative, then at least one pair
    of triangles violates the Delaunay condition, and so the entire
    triangulation is not a Delaunay triangulation.  If VALUE is nonnegative,
    then the triangulation is a Delaunay triangulation.
    It is useful to return VALUE, rather than a simple True/False value,
    because there can be cases where the Delaunay condition is only
    "slightly" violated.  A simple example is a triangulation formed
    by starting with a mesh of squares and dividing each square into
    two triangles by choosing one of the diagonals of the square.
    The Delaunay discrepancy for this mesh, if computed exactly, is 0,
    but roundoff could easily result in discrepancies that were very
    slightly negative.
    If VALUE is positive, and not very small in magnitude, then every
    pair of triangles in the triangulation satisfies the local Delaunay
    condition, and so the triangulation is a Delaunay triangulation.
    If VALUE is negative, and not very small in magnitude, then at least
    one pair of triangles violates the Delaunay condition, and to a
    significant degree.  The triangulation is not a Delaunay triangulation.
    If the magnitude of VALUE is very close to zero, then the triangulation
    is numerically ambiguous.  At least one pair of triangles violates
    or almost violates the condition, but no triangle violates the
    condition to a great extent.  The user must judge whether the
    violation is significant or not.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 September 2009
  Author:
    John Burkardt.
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Input, int TRIANGLE_ORDER, the order of the triangles.
    Input, int TRIANGLE_NUM, the number of triangles in
    the triangulation.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
    the nodes that make up each triangle.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the
    triangle neighbor list.
    Output, double *ANGLE_MIN, the minimum angle that occurred in
    the triangulation.
    Output, int *ANGLE_MIN_TRIANGLE, the triangle in which
    the minimum angle occurred.
    Output, double *ANGLE_MAX, the maximum angle that occurred in
    the triangulation.
    Output, int *ANGLE_MAX_TRIANGLE, the triangle in which
    the maximum angle occurred.
    Output, double TRIANGULATION_DELAUNAY_DISCREPANCY,
    the minimum value of ( A(i,j) - B(i,j) ).
    POSITIVE indicates the triangulation is Delaunay.
    VERY NEAR ZERO is a numerically ambiguous case.
    NEGATIVE indicates the triangulation is not Delaunay.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit2dt2pipitpipitpi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	ityp * node_xy = s_data->a1;
	const register dim_typ triangle_order = s_data->a2;
	const register dim_typ triangle_num = s_data->a3;
	int * triangle_node = s_data->a4;
	int * triangle_neighbor = s_data->a5;
	ityp * angle_min = s_data->a6;
	int * angle_min_triangle = s_data->a7;
	ityp * angle_max = s_data->a8;
	int * angle_max_triangle = s_data->a9;
	
    ityp angle_min1;
    ityp angle_min2;
    ityp *angles1;
    ityp *angles2;
    dim_typ i;
    int i1;
    int i2;
    int i3;
    int i4;
    int n1;
    int n2;
    int n3;
    int n4;
    int neighbor;
    ityp t[6];
    dim_typ triangle_index;
    dim_typ triangle1;
    dim_typ triangle2;
    double value;

    *angle_max = value = 0.00;
    *angle_max_triangle = - 1;
    *angle_min = M_PI;
    *angle_min_triangle = -1;
    /*
    Consider triangle TRIANGLE1
    */
    for ( triangle1 = 0; triangle1 < triangle_num; ++triangle1 )
    {
        /*
        Consider the side opposite vertex NEIGHBOR.
        */
        for ( neighbor = 0; neighbor < 3; ++neighbor )
        {
            triangle2 = triangle_neighbor[neighbor+triangle1*3];
            /*
            There might be no neighbor on side NEIGHBOR.
            */
            /*
            We only need to check a pair of triangles once.
            */
            if ( triangle2 < 0 || triangle2 < triangle1 )
                continue;
            /*
            List the vertices of the quadrilateral in such a way
            that the nodes of triangle 1 come first.

            We rely on a property of the TRIANGLE_NEIGHBOR array, namely, that
            neighbor #1 is on the side opposite to vertex #1, and so on.
            */
            i1 = i4_wrap ( neighbor + 2, 0, 2 );
            i2 = i4_wrap ( neighbor,     0, 2 );
            i3 = i4_wrap ( neighbor + 1, 0, 2 );

            n1 = triangle_node[i1+triangle1*triangle_order];
            n2 = triangle_node[i2+triangle1*triangle_order];
            n3 = triangle_node[i3+triangle1*triangle_order];
            /*
            The "odd" or "opposing" node of the neighboring triangle
            is the one which follows common node I3.
            */
            n4 = -1;
            for ( i = 0; i < 3; ++i )
                if ( triangle_node[i+triangle2*triangle_order] == n3 )
                {
                    i4 = i + 1;
                    i4 = i4_wrap ( i4, 0, 2 );
                    n4 = triangle_node[i4+triangle2*triangle_order];
                    break;
                }

            if ( n4 == -1 )
            {
            	result = MAX_VAL;
                return &result;
            }
            /*
            Compute the minimum angle for (I1,I2,I3) and (I1,I3,I4).
            */
            t[0] = node_xy[(n1<<1)];
            t[1] = node_xy[1+(n1<<1)];
            t[2] = node_xy[(n2<<1)];
            t[3] = node_xy[1+(n2<<1)];
            t[4] = node_xy[(n3<<1)];
            t[5] = node_xy[1+(n3<<1)];
            angles1 = triangle_angles_2d_new ( t );

            t[0] = node_xy[(n1<<1)];
            t[1] = node_xy[1+(n1<<1)];
            t[2] = node_xy[(n3<<1)];
            t[3] = node_xy[1+(n3<<1)];
            t[4] = node_xy[(n4<<1)];
            t[5] = node_xy[1+(n4<<1)];
            angles2 = triangle_angles_2d_new ( t );
            
            ityp min_angles1, min_angles2, max_angles1, max_angles2;
            
        	_MIN(&min_angles1, 3, angles1);
        	_MAX(&min_angles2, 3,  angles2);
            angle_min1 = MIN(min_angles1, min_angles2 );
            _MAX(&max_angles1, 3, angles1);

            if ( *angle_max < max_angles1 )
            {
                *angle_max = max_angles1;
                *angle_max_triangle = triangle1;
            }
            
            _MAX(&max_angles2, 3, angles2);

            if ( *angle_max < max_angles2 )
            {
                *angle_max = max_angles2;
                *angle_max_triangle = triangle2;
            }

            if ( min_angles1 < *angle_min )
            {
                *angle_min = min_angles1;
                *angle_min_triangle = triangle1;
            }

            if ( min_angles2 < *angle_min )
            {
                *angle_min = min_angles2;
                *angle_min_triangle = triangle2;
            }

            free ( angles1 );
            free ( angles2 );
            /*
            Compute the minimum angle for (I1,I2,I4) and (I2,I3,I4).
            */
            t[0] = node_xy[(n1<<1)];
            t[1] = node_xy[1+(n1<<1)];
            t[2] = node_xy[(n2<<1)];
            t[3] = node_xy[(n2<<1)];
            t[4] = node_xy[(n4<<1)];
            t[5] = node_xy[1+(n4<<1)];
            angles1 = triangle_angles_2d_new ( t );

            t[0+0*2] = node_xy[0+n3*2];
            t[1+0*2] = node_xy[1+n3*2];
            t[0+1*2] = node_xy[0+n3*2];
            t[1+1*2] = node_xy[1+n3*2];
            t[0+2*2] = node_xy[0+n4*2];
            t[1+2*2] = node_xy[1+n4*2];
            angles2 = triangle_angles_2d_new ( t );

			_MIN ( &min_angles1, 3, angles1 );
			_MIN ( &min_angles2, 3, angles2 );
            angle_min2 = MIN ( min_angles1, min_angles2 );

            free ( angles1 );
            free ( angles2 );
            /*
            Compare this value to the current minimum.
            */
            value = MIN ( value, angle_min1 - angle_min2 );
        }
    }
    /*
    Scale the results to degrees.
    */
    value *= 180.00 / M_PI;
    *angle_max *= 180.0 / M_PI;
    *angle_min *= 180.00 / M_PI;

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_neighbor_elements ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_NEIGHBOR_ELEMENTS determines element neighbors.
  Discussion:
    A triangulation of a set of nodes can be completely described by
    the coordinates of the nodes, and the list of nodes that make up
    each triangle.  However, in some cases, it is necessary to know
    triangle adjacency information, that is, which triangle, if any,
    is adjacent to a given triangle on a particular side.
    This routine creates a data structure recording this information.
    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
    data items.
    This routine was modified to work with columns rather than rows.
  Example:
    The input information from TRIANGLE_NODE:
    Triangle   Nodes
    --------   ---------------
     1         3      4      1
     2         3      1      2
     3         3      2      8
     4         2      1      5
     5         8      2     13
     6         8     13      9
     7         3      8      9
     8        13      2      5
     9         9     13      7
    10         7     13      5
    11         6      7      5
    12         9      7      6
    13        10      9      6
    14         6      5     12
    15        11      6     12
    16        10      6     11
    The output information in TRIANGLE_NEIGHBOR:
    Triangle  Neighboring Triangles
    --------  ---------------------
     1        -1     -1      2
     2         1      4      3
     3         2      5      7
     4         2     -1      8
     5         3      8      6
     6         5      9      7
     7         3      6     -1
     8         5      4     10
     9         6     10     12
    10         9      8     11
    11        12     10     14
    12         9     11     13
    13        -1     12     16
    14        11     -1     15
    15        16     14     -1
    16        13     15     -1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 September 2009
  Author:
    John Burkardt
  Parameters:
    Input, int TRIANGLE_ORDER, the order of the triangles.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], the nodes that
    make up each triangle.
    Output, int TRIANGLE_ORDER3_NEIGHBOR_TRIANGLES[3*TRIANGLE_NUM],
    the three triangles
    that are direct neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I)
    is the index of the triangle which touches side 1, defined by nodes 2
    and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative if there is no
    neighbor on that side.  In this case, that side of the triangle lies
    on the boundary of the triangulation.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ triangle_order = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	
    int *col;
    dim_typ i;
    dim_typ icol;
    dim_typ j;
    dim_typ k;
    int side1;
    int side2;
    dim_typ tri;
    dim_typ tri1;
    dim_typ tri2;
    int *triangle_neighbor;

    triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
    col = ( int * ) malloc ( 12 * triangle_num * sizeof ( int ) );
    /*
    Step 1.
    From the list of nodes for triangle T, of the form: (I,J,K)
    construct the three neighbor relations:

 (I,J,3,T) or (J,I,3,T),
 (J,K,1,T) or (K,J,1,T),
 (K,I,2,T) or (I,K,2,T)

    where we choose (I,J,3,T) if I < J, or else (J,I,3,T)
    */
    for ( tri = 0; tri < triangle_num; ++tri)
    {
        i = triangle_node[0+tri*triangle_order];
        j = triangle_node[1+tri*triangle_order];
        k = triangle_node[2+tri*triangle_order];

        if ( i < j )
        {
            col[((3*tri+0)<<2)] = i;
            col[1+((3*tri+0)<<2)] = j;
            col[2+((3*tri+0)<<2)] = 3;
            col[3+((3*tri+0)<<2)] = tri + 1;
        }
        else
        {
            col[((3*tri+0)<<2)] = j;
            col[1+((3*tri+0)<<2)] = i;
            col[2+((3*tri+0)<<2)] = 3;
            col[3+((3*tri+0)<<2)] = tri + 1;
        }

        if ( j < k )
        {
            col[((3*tri+1)<<2)] = j;
            col[1+((3*tri+1)<<2)] = k;
            col[2+((3*tri+1)<<2)] = 1;
            col[3+((3*tri+1)<<2)] = tri + 1;
        }
        else
        {
            col[((3*tri+1)<<2)] = k;
            col[1+((3*tri+1)<<2)] = j;
            col[2+((3*tri+1)<<2)] = 1;
            col[3+((3*tri+1)<<2)] = tri + 1;
        }

        if ( k < i )
        {
            col[((3*tri+2)<<2)] = k;
            col[1+((3*tri+2)<<2)] = i;
            col[2+((3*tri+2)<<2)] = 2;
            col[3+((3*tri+2)<<2)] = tri + 1;
        }
        else
        {
            col[((3*tri+2)<<2)] = i;
            col[1+((3*tri+2)<<2)] = k;
            col[2+((3*tri+2)<<2)] = 2;
            col[3+((3*tri+2)<<2)] = tri + 1;
        }
    }
    /*
    Step 2. Perform an ascending dictionary sort on the neighbor relations.
    We only intend to sort on rows 1 and 2; the routine we call here
    sorts on rows 1 through 4 but that won't hurt us.

    What we need is to find cases where two triangles share an edge.
    Say they share an edge defined by the nodes I and J.  Then there are
    two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
    we make sure that these two columns occur consecutively.  That will
    make it easy to notice that the triangles are neighbors.
    */
    i4col_sort_a ( 4, 3*triangle_num, col );
    /*
    Step 3. Neighboring triangles show up as consecutive columns with
    identical first two entries.  Whenever you spot this happening,
    make the appropriate entries in TRIANGLE_NEIGHBOR.
    */
    for ( j = 0; j < triangle_num; ++j )
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            triangle_neighbor[i+j*3] = -1;

    icol = 1;

    for ( ; ; )
    {
        if ( 3 * triangle_num <= icol )
            break;

        if ( col[((icol-1)<<2)] != col[icol<<2] ||
        col[1+((icol-1)<<2)] != col[1+(icol<<2)] )
        {
            ++ icol;
            continue;
        }

        side1 = col[2+((icol-1)<<2)];
        tri1 =  col[3+((icol-1)<<2)];
        side2 = col[2+ (icol   <<2)];
        tri2 =  col[3+ (icol   <<2)];

        triangle_neighbor[side1-1+(tri1-1)*3] = tri2;
        triangle_neighbor[side2-1+(tri2-1)*3] = tri1;

        icol += 2;
    }

    free ( col );

    return triangle_neighbor;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_node_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_NODE_ORDER determines the order of nodes in a triangulation.
  Discussion:
    The order of a node is the number of triangles that use that node
    as a vertex.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 August 2005
  Author:
    John Burkardt
  Parameters:
    Input, integer TRIANGLE_ORDER, the order of the triangulation.
    Input, integer TRIANGLE_NUM, the number of triangles.
    Input, integer TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], the nodes
    that make up the triangles.
    Input, integer NODE_NUM, the number of nodes.
    Output, integer TRIANGULATION_NODE_ORDER[NODE_NUM], the order of
    each node.
*/
{
	const _2dtpii * const s_data = data;
	const register dim_typ triangle_order = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	int node_num = s_data->a3;
	
    dim_typ i;
    dim_typ node;
    int *node_order;
    dim_typ triangle;

    node_order = ( int * ) malloc ( node_num * sizeof ( int ) );

    for ( node = 0; node < node_num; ++node )
        node_order[node] = 0;

    for ( triangle = 0; triangle < triangle_num; ++triangle)
        for ( i = 0; i < triangle_order; ++i )
        {
            node = triangle_node[i+triangle*triangle_order];

            if ( node < 1 || node_num < node )
                return NULL;
            else
                ++ node_order[node-1];
        }

    return node_order;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order3_adj_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
  Discussion:
    This routine is called to count the adjacencies, so that the
    appropriate amount of memory can be set aside for storage when
    the adjacency structure is created.
    The triangulation is assumed to involve 3-node triangles.
    Two nodes are "adjacent" if they are both nodes in some triangle.
    Also, a node is considered to be adjacent to itself.
  Diagram:
       3
    s  |\
    i  | \
    d  |  \
    e  |   \  side 2
       |    \
    3  |     \
       |      \
       1-------2
         side 1
    The local node numbering
   21-22-23-24-25
    |\ |\ |\ |\ |
    | \| \| \| \|
   16-17-18-19-20
    |\ |\ |\ |\ |
    | \| \| \| \|
   11-12-13-14-15
    |\ |\ |\ |\ |
    | \| \| \| \|
    6--7--8--9-10
    |\ |\ |\ |\ |
    | \| \| \| \|
    1--2--3--4--5
    A sample grid.
    Below, we have a chart that summarizes the adjacency relationships
    in the sample grid.  On the left, we list the node, and its neighbors,
    with an asterisk to indicate the adjacency of the node to itself
 (in some cases, you want to count this self adjacency and in some
    you don't).  On the right, we list the number of adjancencies to
    lower-indexed nodes, to the node itself, to higher-indexed nodes,
    the total number of adjacencies for this node, and the location
    of the first and last entries required to list this set of adjacencies
    in a single list of all the adjacencies.
    N   Adjacencies                Below  Self   Above   Total First  Last
   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
    1:  *  2  6                        0     1       2       3     1     3
    2:  1  *  3  6  7                  1     1       3       5     4     8
    3:  2  *  4  7  8                  1     1       3       5     9    13
    4:  3  *  5  8  9                  1     1       3       5    14    18
    5:  4  *  9 10                     1     1       2       4    19    22
    6:  1  2  *  7 11                  2     1       2       5    23    27
    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
   10:  5  9  * 14 15                  2     1       2       5    49    53
   11:  6  7  * 12 16                  2     1       2       5    54    58
   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
   15: 10 14  * 19 20                  2     1       2       5    80    84
   16: 11 12  * 17 21                  2     1       2       5    85    89
   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
   20: 15 19  * 24 25                  2     1       2       5   111   115
   21: 16 17  * 22                     2     1       1       4   116   119
   22: 17 18 21  * 23                  3     1       1       5   120   124
   23: 18 19 22  * 24                  3     1       1       5   125   129
   24: 19 20 23  * 25                  3     1       1       5   130   134
   25: 20 24  *                        2     1       0       3   135   137
   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 August 2006
  Author:
    John Burkardt
  Parameters
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
    make up each triangle, in counterclockwise order.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
    a triangle, lists the neighboring triangle, or -1 if there is
    no neighbor.
    Output, TRIANGULATION_ORDER3_ADJ_COUNT, the number of adjacencies.
    Output, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dt3pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	int * triangle_neighbor = s_data->a3;
	int * adj_col = s_data->a4;
	
    dim_typ adj_num = 0;
    dim_typ i;
    int n1;
    int n2;
    int n3;
    dim_typ node;
    dim_typ triangle;
    dim_typ triangle_order = 3;
    dim_typ triangle2;

    /*
    Set every node to be adjacent to itself.
    */
    for ( node = 0; node < node_num; ++node )
        adj_col[node] = 1;
    /*
    Examine each triangle.
    */
    for ( triangle = 0; triangle < triangle_num; ++triangle)
    {
        n1 = triangle_node[0+triangle*triangle_order];
        n2 = triangle_node[1+triangle*triangle_order];
        n3 = triangle_node[2+triangle*triangle_order];
        /*
        Add edge (1,2) if this is the first occurrence,
        that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
        or if this triangle is the first of the pair in which the edge
        occurs (TRIANGLE < TRIANGLE2).
        */
        triangle2 = triangle_neighbor[0+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ++ adj_col[n1-1];
            ++ adj_col[n2-1];
        }
        /*
        Add edge (2,3).
        */
        triangle2 = triangle_neighbor[1+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ++ adj_col[n2-1];
            ++ adj_col[n3-1];
        }
        /*
        Add edge (3,1).
        */
        triangle2 = triangle_neighbor[2+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ++ adj_col[n1-1];
            ++ adj_col[n3-1];
        }
    }
    /*
    We used ADJ_COL to count the number of entries in each column.
    Convert it to pointers into the ADJ array.
    */
    for ( node = node_num; 1 <= node; --node )
        adj_col[node] = adj_col[node-1];

    adj_col[0] = 1;
    for ( i = 1; i <= node_num; ++i )
        adj_col[i] += adj_col[i-1];

    adj_num = adj_col[node_num] -1;

	result = adj_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_order3_adj_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
  Discussion:
    This routine is called to set the adjacencies, after the
    appropriate amount of memory has been set aside for storage.
    The triangulation is assumed to involve 3-node triangles.
    Two nodes are "adjacent" if they are both nodes in some triangle.
    Also, a node is considered to be adjacent to itself.
    This routine can be used to create the compressed column storage
    for a linear triangle finite element discretization of
    Poisson's equation in two dimensions.
  Diagram:
       3
    s  |\
    i  | \
    d  |  \
    e  |   \  side 2
       |    \
    3  |     \
       |      \
       1-------2
         side 1
    The local node numbering
   21-22-23-24-25
    |\ |\ |\ |\ |
    | \| \| \| \|
   16-17-18-19-20
    |\ |\ |\ |\ |
    | \| \| \| \|
   11-12-13-14-15
    |\ |\ |\ |\ |
    | \| \| \| \|
    6--7--8--9-10
    |\ |\ |\ |\ |
    | \| \| \| \|
    1--2--3--4--5
    A sample grid
    Below, we have a chart that summarizes the adjacency relationships
    in the sample grid.  On the left, we list the node, and its neighbors,
    with an asterisk to indicate the adjacency of the node to itself
 (in some cases, you want to count this self adjacency and in some
    you don't).  On the right, we list the number of adjancencies to
    lower-indexed nodes, to the node itself, to higher-indexed nodes,
    the total number of adjacencies for this node, and the location
    of the first and last entries required to list this set of adjacencies
    in a single list of all the adjacencies.
    N   Adjacencies                Below  Self    Above  Total First  Last
   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
    1:  *  2  6                        0     1       2       3     1     3
    2:  1  *  3  6  7                  1     1       3       5     4     8
    3:  2  *  4  7  8                  1     1       3       5     9    13
    4:  3  *  5  8  9                  1     1       3       5    14    18
    5:  4  *  9 10                     1     1       2       4    19    22
    6:  1  2  *  7 11                  2     1       2       5    23    27
    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
   10:  5  9  * 14 15                  2     1       2       5    49    53
   11:  6  7  * 12 16                  2     1       2       5    54    58
   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
   15: 10 14  * 19 20                  2     1       2       5    80    84
   16: 11 12  * 17 21                  2     1       2       5    85    89
   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
   20: 15 19  * 24 25                  2     1       2       5   111   115
   21: 16 17  * 22                     2     1       1       4   116   119
   22: 17 18 21  * 23                  3     1       1       5   120   124
   23: 18 19 22  * 24                  3     1       1       5   125   129
   24: 19 20 23  * 25                  3     1       1       5   130   134
   25: 20 24  *                        2     1       0       3   135   137
   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 August 2006
  Author:
    John Burkardt
  Parameters
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
    make up each triangle in counterclockwise order.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
    a triangle, lists the neighboring triangle, or -1 if there is
    no neighbor.
    Input, int ADJ_NUM, the number of adjacencies.
    Input, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
    Output, int TRIANGULATION_ORDER3_ADJ_SET[ADJ_NUM], the adjacency
    information.
*/
{
	const _2dt2pidtpi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	int * triangle_neighbor = s_data->a3;
	const register dim_typ adj_num = s_data->a4;
	int * adj_col = s_data->a5;
	
    int *adj;
    int *adj_copy;
    dim_typ k;
    int k1;
    int k2;
    int n1;
    int n2;
    int n3;
    dim_typ node;
    dim_typ triangle;
    dim_typ triangle2;
    dim_typ triangle_order = 3;

    adj = ( int * ) malloc ( adj_num * sizeof ( int ) );
    for ( k = 0; k < adj_num; ++k )
        adj[k] = -1;


    adj_copy = ( int * ) malloc ( node_num * sizeof ( int ) );
    for ( node = 0; node < node_num; ++node )
        adj_copy[node] = adj_col[node];
    /*
    Set every node to be adjacent to itself.
    */
    for ( node = 1; node <= node_num; ++node )
    {
        adj[adj_copy[node-1]-1] = node;
        ++ adj_copy[node-1];
    }
    /*
    Examine each triangle.
    */
    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        n1 = triangle_node[0+triangle*triangle_order];
        n2 = triangle_node[1+triangle*triangle_order];
        n3 = triangle_node[2+triangle*triangle_order];
        /*
        Add edge (1,2) if this is the first occurrence,
        that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
        or if this triangle is the first of the pair in which the edge
        occurs (TRIANGLE < TRIANGLE2).
        */
        triangle2 = triangle_neighbor[0+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj[adj_copy[n1-1]-1] = n2;
            ++ adj_copy[n1-1];
            adj[adj_copy[n2-1]-1] = n1;
            ++ adj_copy[n2-1];
        }
        /*
        Add edge (2,3).
        */
        triangle2 = triangle_neighbor[1+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj[adj_copy[n2-1]-1] = n3;
            ++ adj_copy[n2-1];
            adj[adj_copy[n3-1]-1] = n2;
            ++ adj_copy[n3-1];
        }
        /*
        Add edge (3,1).
        */
        triangle2 = triangle_neighbor[2+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj[adj_copy[n1-1]-1] = n3;
            ++ adj_copy[n1-1];
            adj[adj_copy[n3-1]-1] = n1;
            ++ adj_copy[n3-1];
        }
    }
    /*
    Ascending sort the entries for each node.
    */
    for ( node = 1; node <= node_num; ++node )
    {
        k1 = adj_col[node-1];
        k2 = adj_col[node]-1;
        i4vec_sort_heap_a ( k2+1-k1, adj+k1-1 );
    }

    free ( adj_copy );

    return adj;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order3_adj_set2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_ADJ_SET2 sets adjacencies in a triangulation.
  Discussion:
    This routine is called to set up the arrays IA and JA that
    record which nodes are adjacent in a triangulation.
    The triangulation is assumed to involve 3-node triangles.
    Two nodes are "adjacent" if they are both nodes in some triangle.
    Also, a node is considered to be adjacent to itself.
    This routine can be used to create the compressed column storage
    for a linear triangle finite element discretization of
    Poisson's equation in two dimensions.
  Diagram:
       3
    s  |\
    i  | \
    d  |  \
    e  |   \  side 2
       |    \
    3  |     \
       |      \
       1-------2
         side 1
    The local node numbering
   21-22-23-24-25
    |\ |\ |\ |\ |
    | \| \| \| \|
   16-17-18-19-20
    |\ |\ |\ |\ |
    | \| \| \| \|
   11-12-13-14-15
    |\ |\ |\ |\ |
    | \| \| \| \|
    6--7--8--9-10
    |\ |\ |\ |\ |
    | \| \| \| \|
    1--2--3--4--5
    A sample grid
    Below, we have a chart that summarizes the adjacency relationships
    in the sample grid.  On the left, we list the node, and its neighbors,
    with an asterisk to indicate the adjacency of the node to itself
 (in some cases, you want to count this self adjacency and in some
    you don't).  On the right, we list the number of adjancencies to
    lower-indexed nodes, to the node itself, to higher-indexed nodes,
    the total number of adjacencies for this node, and the location
    of the first and last entries required to list this set of adjacencies
    in a single list of all the adjacencies.
    N   Adjacencies                Below  Self    Above  Total First  Last
   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
    1:  *  2  6                        0     1       2       3     1     3
    2:  1  *  3  6  7                  1     1       3       5     4     8
    3:  2  *  4  7  8                  1     1       3       5     9    13
    4:  3  *  5  8  9                  1     1       3       5    14    18
    5:  4  *  9 10                     1     1       2       4    19    22
    6:  1  2  *  7 11                  2     1       2       5    23    27
    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
   10:  5  9  * 14 15                  2     1       2       5    49    53
   11:  6  7  * 12 16                  2     1       2       5    54    58
   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
   15: 10 14  * 19 20                  2     1       2       5    80    84
   16: 11 12  * 17 21                  2     1       2       5    85    89
   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
   20: 15 19  * 24 25                  2     1       2       5   111   115
   21: 16 17  * 22                     2     1       1       4   116   119
   22: 17 18 21  * 23                  3     1       1       5   120   124
   23: 18 19 22  * 24                  3     1       1       5   125   129
   24: 19 20 23  * 25                  3     1       1       5   130   134
   25: 20 24  *                        2     1       0       3   135   137
   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
    For this example, the initial portion of the IA and JA arrays will be:
   (1,1), (1,2), (1,6),
   (2,1), (2,2), (2,3), (2,6), (2,7),
   (3,2), (3,3), (3,4), (3,7), (3,8),
     ...
   (25,20), (25,24), (25,25)
    for a total of 137 pairs of values.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 July 2007
  Author:
    John Burkardt
  Parameters:
    Input int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
    make up each triangle in counterclockwise order.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
    a triangle, lists the neighboring triangle, or -1 if there is
    no neighbor.
    Input, int ADJ_NUM, the number of adjacencies.
    Input, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
    Output, int IA[ADJ_NUM], JA[ADJ_NUM], the adjacency information.
*/
{
	const _2dt2pidt3pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	int * triangle_neighbor = s_data->a3;
	const register dim_typ adj_num = s_data->a4;
	int * adj_col = s_data->a5;
	int * ia = s_data->a6;
	int * ja = s_data->a7;
	
    dim_typ adj;
    int *adj_copy;
    dim_typ k;
    int k1;
    int k2;
    int n1;
    int n2;
    int n3;
    dim_typ node;
    dim_typ triangle;
    dim_typ triangle2;
    dim_typ triangle_order = 3;

    for ( adj = 0; adj < adj_num; ++adj)
        ia[adj] = ja[adj] = -1;


    adj_copy = ( int * ) malloc ( node_num * sizeof ( int ) );
    for ( node = 0; node < node_num; ++node)
        adj_copy[node] = adj_col[node];
    /*
    Set every node to be adjacent to itself.
    */
    for ( node = 1; node <= node_num; ++node )
    {
        ia[adj_copy[node-1]-1] = node;
        ja[adj_copy[node-1]-1] = node;
        ++ adj_copy[node-1];
    }
    /*
    Examine each triangle.
    */
    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        n1 = triangle_node[0+triangle*triangle_order];
        n2 = triangle_node[1+triangle*triangle_order];
        n3 = triangle_node[2+triangle*triangle_order];
        /*
        Add edge (1,2) if this is the first occurrence,
        that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
        or if this triangle is the first of the pair in which the edge
        occurs (TRIANGLE < TRIANGLE2).
        */
        triangle2 = triangle_neighbor[0+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ia[adj_copy[n1-1]-1] = n1;
            ja[adj_copy[n1-1]-1] = n2;
            ++ adj_copy[n1-1];

            ia[adj_copy[n2-1]-1] = n2;
            ja[adj_copy[n2-1]-1] = n1;
            ++ adj_copy[n2-1];
        }
        /*
        Add edge (2,3).
        */
        triangle2 = triangle_neighbor[1+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ia[adj_copy[n2-1]-1] = n2;
            ja[adj_copy[n2-1]-1] = n3;
            ++ adj_copy[n2-1];

            ia[adj_copy[n3-1]-1] = n3;
            ja[adj_copy[n3-1]-1] = n2;
            ++ adj_copy[n3-1];
        }
        /*
        Add edge (3,1).
        */
        triangle2 = triangle_neighbor[2+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ia[adj_copy[n1-1]-1] = n1;
            ja[adj_copy[n1-1]-1] = n3;
            ++ adj_copy[n1-1];

            ia[adj_copy[n3-1]-1] = n3;
            ja[adj_copy[n3-1]-1] = n1;
            ++ adj_copy[n3-1];
        }
    }
    /*
    Lexically sort the IA, JA values.
    */
    i4vec2_sort_a ( adj_num, ia, ja );

    free ( adj_copy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_order3_adjacency ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_ADJACENCY computes the full adjacency matrix
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 2014
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes in the
    triangulation.
    Input, int ELEMENT_NUM, the number of triangles in
    the triangulation.
    Input, int ELEMENT_NODE[3*ELEMENT_NUM],
    the nodes making up each triangle.
    Output, int TRIANGULATION_ORDER3_ADJACENCY[NODE_NUM*NODE_NUM], the adjacency
    matrix.  ADJ(I,J) is 1 if nodes I and J are adjacent, that is,
    they are immediate neighbors on an edge of the triangulation.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ element_num = s_data->a1;
	int * element_node = s_data->a2;
	
    int *adj;
    dim_typ element;
    dim_typ i, j, k;

    adj = ( int * ) malloc ( node_num * node_num * sizeof ( int ) );

    for ( j = 0; j < node_num; ++j )
        for ( i = 0; i < node_num; ++i )
            adj[i+j*node_num] = 0;

    for ( element = 0; element < element_num; ++element )
    {
        i = element_node[0+element*3];
        j = element_node[1+element*3];
        k = element_node[2+element*3];

        adj[i+j*node_num] = adj[i+k*node_num] = adj[j+i*node_num] = adj[j+k*node_num] =  adj[k+i*node_num] = adj[k+j*node_num] = 1;
    }

    return adj;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order3_boundary_edge_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the boundary edges.
  Discussion:
    This routine is given a triangulation, an abstract list of triples
    of nodes.  It is assumed that the nodes in each triangle are listed
    in a counterclockwise order, although the routine should work
    if the nodes are consistently listed in a clockwise order as well.
    It is assumed that each edge of the triangulation is either
    * an INTERIOR edge, which is listed twice, once with positive
      orientation and once with negative orientation, or;
    * a BOUNDARY edge, which will occur only once.
    This routine should work even if the region has holes - as long
    as the boundary of the hole comprises more than 3 edges!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 June 2005
  Author:
    John Burkardt
  Parameters:
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
    triangles.  These should be listed in counterclockwise order.
    Output, integer TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT, the number
    of boundary edges.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ triangle_num = s_data->a0;
	int * triangle_node = s_data->a1;
	
    dim_typ boundary_edge_num;
    int e1;
    int e2;
    int *edge;
    dim_typ i;
    dim_typ interior_edge_num;
    dim_typ j;
    dim_typ m;
    dim_typ n;
    dim_typ unique_num;

    m = 2;
    n = 3 * triangle_num;
    /*
    Set up the edge array.
    */
    edge = ( int * ) malloc ( m * n * sizeof ( int ) );

    for ( j = 0; j < triangle_num; ++j )
    {
        edge[0+(j               )*m] = triangle_node[0+j*3];
        edge[1+(j               )*m] = triangle_node[1+j*3];
        edge[0+(j+  triangle_num)*m] = triangle_node[1+j*3];
        edge[1+(j+  triangle_num)*m] = triangle_node[2+j*3];
        edge[0+(j+2*triangle_num)*m] = triangle_node[2+j*3];
        edge[1+(j+2*triangle_num)*m] = triangle_node[0+j*3];
    }
    /*
    In each column, force the smaller entry to appear first.
    */
    for ( j = 0; j < n; ++j )
    {
        e1 = MIN ( edge[0+j*m], edge[1+j*m] );
        e2 = MAX ( edge[0+j*m], edge[1+j*m] );
        edge[0+j*m] = e1;
        edge[1+j*m] = e2;
    }
    /*
    Ascending sort the column array.
    */
    i4col_sort_a ( m, n, edge );
    /*
    Get the number of unique columns in EDGE.
    */
    unique_num = i4col_sorted_unique_count ( m, n, edge );
    interior_edge_num = 3 * triangle_num - unique_num;
    boundary_edge_num = 3 * triangle_num - (interior_edge_num<<1);
    free ( edge );

	result = boundary_edge_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order3_boundary_edge_count_euler ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
  Discussion:
    We assume we are given information about a triangulation
    of a set of nodes in the plane.
    Given the number of nodes and triangles, we are going to apply
    Euler's formula to determine the number of edges that lie on the
    boundary of the set of nodes.
    The number of faces, including the infinite face and internal holes,
    is TRIANGLE_NUM + HOLE_NUM + 1.
    Let BOUNDARY_NUM denote the number of edges on the boundary.
    Each of the TRIANGLE_NUM triangles uses three edges.  Every edge
    occurs in two different faces, so the number of edges must be
 ( 3 * TRIANGLE_NUM + BOUNDARY_NUM ) / 2.
    The number of nodes used in the triangulation is NODE_NUM.
    Euler's formula asserts that, for a simple connected figure in the
    plane with no edge crossings, NODE_NUM nodes, EDGE_NUM edges and
    FACE_NUM faces:
      NODE_NUM - EDGE_NUM + FACE_NUM = 2
    In our context, this becomes
      NODE_NUM - ( 3 * TRIANGLE_NUM + BOUNDARY_NUM ) / 2
      + TRIANGLE_NUM + HOLE_NUM + 1 = 2
    or
      BOUNDARY_NUM = 2 * NODE_NUM + 2 * HOLE_NUM - TRIANGLE_NUM - 2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 June 2005
  Author:
    John Burkardt
  Reference:
    Marc deBerg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
    Computational Geometry,
    Springer, 2000,
    ISBN: 3-540-65620-0.
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int HOLE_NUM, the number of holes.
    Output, int TRIANGULATION_BOUNDARY_COUNT, the number of edges that
    lie on the convex hull of the triangulation.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data; 
	const register dim_typ node_num = a_data[0];
	const register dim_typ triangle_num = a_data[1];
	const register dim_typ hole_num = a_data[2];
	
	result = ( (node_num<<1) + (hole_num<<1) - triangle_num - 2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_order3_boundary_node ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_BOUNDARY_NODE indicates nodes on the boundary.
  Discussion:
    This routine is given a triangulation, an abstract list of triples
    of nodes.  It is assumed that the nodes in each triangle are listed
    in a counterclockwise order, although the routine should work
    if the nodes are consistently listed in a clockwise order as well.
    It is assumed that each edge of the triangulation is either
    * an INTERIOR edge, which is listed twice, once with positive
      orientation and once with negative orientation, or;
    * a BOUNDARY edge, which will occur only once.
    This routine should work even if the region has holes - as long
    as the boundary of the hole comprises more than 3 edges!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 January 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
    triangles.  These should be listed in counterclockwise order.
    Output, int TRIANGULATION_ORDER3_BOUNDARY_NODE[NODE_NUM],
    is TRUE if the node is on a boundary edge.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	
    int e1;
    int e2;
    int *edge;
    bool equal;
    dim_typ i;
    dim_typ j;
    dim_typ m;
    dim_typ n;
    int *node_boundary;

    m = 2;
    n = 3 * triangle_num;
    /*
    Set up the edge array.
    */
    edge = ( int * ) malloc ( m * n * sizeof ( int ) );

    for ( j = 0; j < triangle_num; ++j)
    {
        edge[0+(j               )*m] = triangle_node[0+j*3];
        edge[1+(j               )*m] = triangle_node[1+j*3];
        edge[0+(j+  triangle_num)*m] = triangle_node[1+j*3];
        edge[1+(j+  triangle_num)*m] = triangle_node[2+j*3];
        edge[0+(j+2*triangle_num)*m] = triangle_node[2+j*3];
        edge[1+(j+2*triangle_num)*m] = triangle_node[0+j*3];
    }
    /*
    In each column, force the smaller entry to appear first.
    */
    for ( j = 0; j < n; ++j )
    {
        e1 = MIN ( edge[0+j*m], edge[1+j*m] );
        e2 = MAX ( edge[0+j*m], edge[1+j*m] );
        edge[0+j*m] = e1;
        edge[1+j*m] = e2;
    }
    /*
    Ascending sort the column array.
    */
    i4col_sort_a ( m, n, edge );
    /*
    Records which appear twice are internal edges and can be ignored.
    */
    node_boundary = ( int * ) malloc ( node_num * sizeof ( int ) );

    for ( i = 0; i < node_num; ++i)
        node_boundary[i] = 0;

    j = 0;

    while ( j < 3 * triangle_num )
    {
        ++ j;

        if ( j == 3 * triangle_num )
        {
            for ( i = 0; i < m; ++i )
                node_boundary[edge[i+(j-1)*m]-1] = 1;
            break;
        }

        equal = true;

        for ( i = 0; i < m; ++i )
            if ( edge[i+(j-1)*m] != edge[i+j*m] )
                equal = false;


        if ( equal )
            ++ j;
        else
            for ( i = 0; i < m; ++i )
                node_boundary[edge[i+(j-1)*m]-1] = 1;

    }

    free ( edge );

    return node_boundary;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order3_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_CHECK makes some simple checks on a triangulation.
  Discussion:
    Because this routine does not receive the physical coordinates of
    the nodes, it cannot ensure that the triangulation is maximal,
    that is, that no more triangles can be created.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 June 2005
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
    triangles.  These should be listed in counterclockwise order.
    Output, int TRIANGULATION_CHECK, error flag.
    0, no error occurred.
    nonzero, an error occurred, the triangulation is not valid.
*/
{
	static bool result = UCHAR_MAX;
	
	const _2dtpi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	
    short boundary_num;
    int error;
    int euler;
    dim_typ i;
    dim_typ j;
    int *used;
    /*
    Checks 1 and 2:
    node_num must be at least 3.
    TRIANGLE_NUM must be at least 1.
    */
    if ( node_num < 3 )
    {
    	result = 1;
        return &result;
    }

    if ( triangle_num == 0 )
    {
    	result = 2;
        return &result;
    }
    /*
    Checks 3 and 4:
    Verify that all node values are greater than or equal to 1
    and less than or equal to node_num.
    */
    for ( j = 0; j < triangle_num; ++j)
        for ( i = 0; i < 3; ++i)
            if ( triangle_node[i+j*3] < 1 )
            {
            	result = 3;
        		return &result;
            }

    for ( j = 0; j < triangle_num; ++j )
        for ( i = 0; i < 3; ++i)
            if ( node_num < triangle_node[i+j*3] )
            {
            	result = 4;
        		return &result;
            }

    /*
    Check 5:
    Verify that every node is used at least once.
    */
    used = ( int * ) malloc ( node_num * sizeof ( int ) );

    for ( i = 0; i < node_num; ++i )
        used[i] = 0;

    for ( j = 0; j < triangle_num; ++j )
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            ++ used[triangle_node[i+j*3]-1];


    for ( i = 0; i < node_num; ++i )
        if ( used[i] == 0 )
        {
        	result = 5;
        	return &result;
        }
    free ( used );
    /*
    Check 6:
    Verify that no node is repeated in a triangle.
    */
    for ( j = 0; j < triangle_num; ++j )
        if ( triangle_node[0+j*3] == triangle_node[1+j*3] ||triangle_node[1+j*3] == triangle_node[2+j*3] ||triangle_node[2+j*3] == triangle_node[0+j*3] )
        {
        	result = 6;
        	return &result;
        }

    /*
    Check 7:
    Verify that no edge is repeated, and that repeated edges occur in
    negated pairs.
    */
    boundary_num = triangulation_order3_edge_check ( triangle_num,triangle_node );

    if ( boundary_num < 0 )
    {
    	result = 7;
        return &result;
    }
    /*
    Check 8:
    Does the triangulation satisfy Euler's criterion?
    If not, then the triangulation is not proper. (For instance, there
    might be a hole in the interior.)
    */
    euler = boundary_num + triangle_num + 2 - (node_num<<1);

    if ( euler != 0 )
    {
    	result = 8;
        return &result;
    }

    result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order3_edge_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_EDGE_CHECK checks the edges of a triangulation.
  Discussion:
    Converted from a row-based to a column-based calculation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 February 2007
  Author:
    John Burkardt
  Parameters:
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
    each triangle.
    Output, int TRIANGULATION_EDGE_CHECK is negative if an error was
    detected; otherwise, it is the number of edges that lie on the boundary.
*/
{
	static short result = SHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ triangle_num = s_data->a0;
	int * triangle_node = s_data->a1;
	
    short boundary_num;
    dim_typ i, j, k;
    int *col;
    dim_typ tri;
    dim_typ triangle_order = 3;
    /*
    Step 1.
    From the list of nodes for triangle T, of the form: (I,J,K)
    construct the three neighbor relations:

 (I,J,+1) or (J,I,-1),
 (J,K,+1) or (K,J,-1),
 (K,I,+1) or (I,K,-1)

    where we choose (I,J,+1) if I < J, or else (J,I,-1) and so on.
    */
    col = ( int * ) malloc ( 9 * triangle_num * sizeof ( int ) );

    for ( tri = 0; tri < triangle_num; ++tri )
    {
        i = triangle_node[0+tri*triangle_order];
        j = triangle_node[1+tri*triangle_order];
        k = triangle_node[2+tri*triangle_order];

        if ( i < j )
        {
            col[0+(3*tri+0)*3] =  i;
            col[1+(3*tri+0)*3] =  j;
            col[2+(3*tri+0)*3] = +1;
        }
        else
        {
            col[0+(3*tri+0)*3] =  j;
            col[1+(3*tri+0)*3] =  i;
            col[2+(3*tri+0)*3] = -1;
        }

        if ( j < k )
        {
            col[0+(3*tri+1)*3] =  j;
            col[1+(3*tri+1)*3] =  k;
            col[2+(3*tri+1)*3] = +1;
        }
        else
        {
            col[0+(3*tri+1)*3] =  k;
            col[1+(3*tri+1)*3] =  j;
            col[2+(3*tri+1)*3] = -1;
        }

        if ( k < i )
        {
            col[0+(3*tri+2)*3] =  k;
            col[1+(3*tri+2)*3] =  i;
            col[2+(3*tri+2)*3] = +1;
        }
        else
        {
            col[0+(3*tri+2)*3] =  i;
            col[1+(3*tri+2)*3] =  k;
            col[2+(3*tri+2)*3] = -1;
        }
    }
    /*
    Step 2. Perform an ascending dictionary sort on the neighbor relations.
    */
    i4col_sort_a ( 3, 3*triangle_num, col );
    /*
    Step 3.

    If any record occurs twice, we have an error.
    Unpaired records lie on the convex hull.
    */
    i = boundary_num = 0;

    while ( i < 3 * triangle_num )
    {
        ++ i;

        if ( i == 3 * triangle_num )
            ++ boundary_num;
        else
        {
            if ( col[0+(i-1)*3] == col[0+i*3] &&col[1+(i-1)*3] == col[1+i*3] )
            {
                if ( col[2+(i-1)*3] == col[2+i*3] )
                {
                	result = -1;
                    return &result;
                }
                else
                    ++ i;
            }
            else
                ++ boundary_num;
        }
    }

    free ( col );

	result = boundary_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order3_neighbor ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_NEIGHBOR determines a neighbor of a given triangle.
  Discussion:
    A set of nodes is given.  A triangulation of the nodes has been
    defined and recorded in TRIANGLE_NODE.  The TRIANGLE_NODE data structure
    records triangles as sets of three nodes, N1, N2, N3, that implicitly
    define three sides, being the line segments N1-N2, N2-N3, and N3-N1.
    The nodes of the triangle are listed in counterclockwise order.
    This means that if two triangles share a side, then the nodes
    defining that side occur in the order (N1,N2) for one triangle,
    and (N2,N1) for the other.
    The routine is given a triangle and a side, and asked to find
    another triangle (if any) that shares that side.  The routine
    simply searches the TRIANGLE_NODE structure for an occurrence of the
    nodes in the opposite order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 October 2003
  Author:
    John Burkardt
  Parameters:
    Input, int TRIANGLE_NUM, the number of triangles.
    Input/output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that define
    each triangle.
    Input, int T1, the index of the triangle.
    Input, int S1, the index of the triangle side.
    Output, int *T2, the index of the triangle which is the neighbor
    to T1 on side S1, or -1 if there is no such neighbor.
    Output, int *S2, the index of the side of triangle T2 which
    is shared with triangle T1, or -1 if there is no such neighbor.
*/
{
	const _2dt2pidtpi * const s_data = data;
	
	dim_typ t1 = s_data->a0;
	dim_typ s1 = s_data->a1;
	int * t2 = s_data->a2;
	int * s2 = s_data->a3;
	const register dim_typ triangle_num = s_data->a4;
	int * triangle_node = s_data->a5;
	
    int n1;
    int n2;
    dim_typ s;
    int ss;
    dim_typ t;

    n1 = triangle_node[s1-1+(t1-1)*3];
    ss = i4_wrap ( s1+1, 1, 3 );
    n2 = triangle_node[ss-1+(t1-1)*3];

    for ( t = 0; t < triangle_num; ++t )
        for ( s = 0; s < 3; ++s)
            if ( triangle_node[s+t*3] == n1 )
            {
                ss = i4_wrap ( s-1, 0, 2 );
                if ( triangle_node[ss+t*3] == n2 )
                {
                    *t2 = t + 1;
                    *s2 = ss + 1;
                    return NULL;
                }
            }

    *t2 = *s2 = -1;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   triangulation_order3_quad ( const register dim_typ node_num, ityp node_xy[static node_num<<1],
  const register dim_typ triangle_order, const register dim_typ triangle_num, int triangle_node[static triangle_num*triangle_order],void quad_fun ( int n, ityp xy_vec[], ityp f_vec[] ), const register dim_typ quad_num,ityp quad_xy[static quad_num<<1], ityp quad_w[static quad_num], ityp *quad_value, ityp *region_area )
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_QUAD approximates an integral over a triangulation.
  Discussion:
    The routine will accept triangulations of order higher than 3.
    However, only the first three nodes (the vertices) of each
    triangle will be used.  This will still produce correct results
    for higher order triangulations, as long as the sides of the
    triangle are straight.
    We assume that the vertices of each triangle are listed first
    in the description of higher order triangles, and we assume that
    the vertices are listed in counterclockwise order.
    The approximation of the integral is made using a quadrature rule
    defined on the unit triangle, and supplied by the user.
    The user also supplies the name of a subroutine, here called "QUAD_FUN",
    which evaluates the integrand at a set of points.  The form is:
      void quad_fun ( int n, double xy_vec[], double f_vec[] )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes in the triangulation.
    Input, double NODE_XY(2,NODE_NUM), the coordinates of the nodes.
    Input, int TRIANGLE_ORDER, the order of triangles in the triangulation.
    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
    the nodes making up each triangle.
    Input, void QUAD_FUN ( int N, double XY_VEC[], double F_VEC[] ),
    the name of the function that evaluates the integrand.
    Input, int QUAD_NUM, the order of the quadrature rule.
    Input, double QUAD_XY(2,QUAD_NUM), the abscissas of the
    quadrature rule, in the unit triangle.
    Input, double QUAD_W(QUAD_NUM), the weights of the
    quadrature rule.
    Output, double *QUAD_VALUE, the estimate of the integral
    of F(X,Y) over the region covered by the triangulation.
    Output, double *REGION_AREA, the area of the region.
*/
{
    dim_typ i, j;
    dim_typ quad;
    ityp *quad_f;
    ityp *quad2_xy;
    ityp temp;
    dim_typ triangle;
    ityp triangle_area;
    ityp triangle_xy[6];

    quad_f = ( ityp * ) malloc ( quad_num * sizeof ( ityp ) );
    quad2_xy = ( ityp * ) malloc ( quad_num * sizeof ( ityp ) << 1 );

    *quad_value = *region_area = 0.00;

    for ( triangle = 0; triangle < triangle_num; ++triangle)
    {
        #pragma omp parallel for num_threads(3)
        for ( j = 0; j < 3; ++j )
            #pragma omp parallel for num_threads(2)
            for ( i = 0; i < 2; ++i )
                triangle_xy[i+(j<<1)] = node_xy[i+((triangle_node[j+triangle*3]-1)<<1)];

        triangle_area = triangle_area_2d ( triangle_xy );

        triangle_order3_reference_to_physical ( triangle_xy,quad_num, quad_xy, quad2_xy );

        quad_fun ( quad_num, quad2_xy, quad_f );

        temp = 0.00;
        for ( quad = 0; quad < quad_num; ++quad )
            temp += quad_w[quad] * quad_f[quad];

        *quad_value += triangle_area * temp;

        *region_area += triangle_area;
    }

    free ( quad_f );
    free ( quad2_xy );

    return;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _triangulation_order3_refine_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_REFINE_COMPUTE computes a refined order 3 triangulation.
  Discussion:
    Given a triangle defined by nodes 1, 2, 3, we need to generate
    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
    and T4.
    The task is more complicated by the fact that we are working with
    a mesh of triangles, so that we want to create a node only once,
    even though it may be shared by other triangles.
          3
         / \
        /T3 \
      13----23
      / \T4 / \
     /T1 \ /T2 \
    1----12-----2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM1, the number of nodes.
    Input, int TRIANGLE_NUM1, the number of triangles.
    Input, double NODE_XY1[2*NODE_NUM1], the nodes.
    Input, int TRIANGLE_NODE1[3*TRIANGLE_NUM1], the nodes that make up the
    triangles.  These should be listed in counterclockwise order.
    Input, int NODE_NUM2, the number of nodes in the refined mesh.
    Input, int TRIANGLE_NUM2, the number of triangles in the refined mesh.
    Input, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge information computed
    by TRIANGULATION_ORDER3_REFINE_SIZE.
    Output, double NODE_XY2[2*NODE_NUM2], the refined nodes.
    Output, int TRIANGLE_NODE2[3*TRIANGLE_NUM2], the nodes that make up the
    triangles in the refined mesh.
*/
{
	const _2dtpitpi2dtpipitpi * const s_data = data;
	const register dim_typ node_num1 = s_data->a0;
	const register dim_typ triangle_num1 = s_data->a1;
	ityp * node_xy1 = s_data->a2;
	int * triangle_node1 = s_data->a3;
	const register dim_typ node_num2 = s_data->a4;
	const register dim_typ triangle_num2 = s_data->a5;
	int * edge_data = s_data->a6;
	ityp * node_xy2 = s_data->a7;
	int * triangle_node2 = s_data->a8;
	
    dim_typ edge;
    dim_typ i;
    dim_typ j;
    int n1;
    int n1_old;
    int n2;
    int n2_old;
    dim_typ node;
    dim_typ triangle1;
    int v1;
    int v2;
    /*
    Copy the old nodes.
    */
    for ( j = 0; j < node_num1; ++j )
        #pragma omp parallel for num_threads(2)
        for ( i = 0; i < 2; ++i )
            node_xy2[i+j*2] = node_xy1[i+j*2];
    for ( j = 0; j < triangle_num2; ++j )
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i)
            triangle_node2[i+j*3] = -1;
    /*
    We can assign the existing nodes to the new triangles.
    */
    for ( triangle1 = 0; triangle1 < triangle_num1; ++triangle1 )
    {
        triangle_node2[0+((triangle1<<2)+0)*3] = triangle_node1[0+triangle1*3];
        triangle_node2[1+((triangle1<<2)+1)*3] = triangle_node1[1+triangle1*3];
        triangle_node2[2+((triangle1<<2)+2)*3] = triangle_node1[2+triangle1*3];
    }

    node = node_num1;

    n1_old = -1;
    n2_old = -1;

    for ( edge = 0; edge < 3 * triangle_num1; ++edge )
    {
        n1 = edge_data[0+edge*5] - 1;
        n2 = edge_data[1+edge*5] - 1;
        /*
        If this edge is new, create the coordinates and index for this node.
        */
        if ( n1 != n1_old || n2 != n2_old )
        {

            if ( node_num2 < node )
                return NULL;

            #pragma omp parallel for num_threads(2)
            for ( i = 0; i < 2; ++i )
                node_xy2[i+(node<<1)] = ( node_xy2[i+(n1<<1)] + node_xy2[i+(n2<<1)] ) / 2.00;

            ++ node;

            n1_old = n1;
            n2_old = n2;
        }
        /*
        Assign the node to triangles.
        */
        v1 = edge_data[2+edge*5];
        v2 = edge_data[3+edge*5];
        triangle1 = edge_data[4+edge*5];

        if ( v1 == 1 && v2 == 2 )
        {
            triangle_node2[0+((triangle1<<2)+1)*3] = node;
            triangle_node2[1+((triangle1<<2)+0)*3] = node;
            triangle_node2[2+((triangle1<<2)+3)*3] = node;
        }
        else if ( v1 == 1 && v2 == 3 )
        {
            triangle_node2[0+((triangle1<<2)+2)*3] = node;
            triangle_node2[1+((triangle1<<2)+3)*3] = node;
            triangle_node2[2+((triangle1<<2)+0)*3] = node;
        }
        else if ( v1 == 2 && v2 == 3 )
        {
            triangle_node2[0+((triangle1<<2)+3)*3] = node;
            triangle_node2[1+((triangle1<<2)+2)*3] = node;
            triangle_node2[2+((triangle1<<2)+1)*3] = node;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order3_refine_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_REFINE_SIZE sizes a refined order 3 triangulation.
  Discussion:
    Given a triangle defined by nodes 1, 2, 3, we need to generate
    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
    and T4.
    The task is more complicated by the fact that we are working with
    a mesh of triangles, so that we want to create a node only once,
    even though it may be shared by other triangles.
          3
         / \
        /T3 \
      13----23
      / \T4 / \
     /T1 \ /T2 \
    1----12-----2
    This routine simply determines the sizes of the resulting node
    and triangle arrays.
    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
    data items, one item for every edge of every triangle.  Each
    data item records, for a given edge, the global indices
    of the two endpoints, the local indices of the two endpoints,
    and the index of the triangle.
    Through careful sorting, it is possible to arrange this data in
    a way that allows the proper generation of the interpolated nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM1, the number of nodes in the original mesh.
    Input, int  TRIANGLE_NUM1, the number of triangles in the
    original mesh.
    Input, int TRIANGLE_NODE1[3*TRIANGLE_NUM1], the indices of the nodes
    that form the triangles in the input mesh.
    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
    Output, int *TRIANGLE_NUM2, the number of triangles in the
    refined mesh.
    Output, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge data that will
    be needed by TRIANGULATION_ORDER3_REFINE_COMPUTE.
*/
{
	const _2dt4pi * const s_data = data;
	register dim_typ node_num1 = s_data->a0;
	const register dim_typ triangle_num1 = s_data->a1;
	int * triangle_node1 = s_data->a2;
	int * node_num2 = s_data->a3;
	int * triangle_num2 = s_data->a4;
	int * edge_data = s_data->a5;
	
    int a;
    int b;
    dim_typ edge;
    dim_typ i;
    dim_typ j;
    dim_typ k;
    int n1;
    int n1_old;
    int n2;
    int n2_old;
    dim_typ triangle;
    /*
    Step 1.
    From the list of nodes for triangle T, of the form: (I,J,K)
    construct the edge relations:

 (I,J,1,2,T)
 (I,K,1,3,T)
 (J,K,2,3,T)

    In order to make matching easier, we reorder each pair of nodes
    into ascending order.
    */
    for ( triangle = 0; triangle < triangle_num1; ++triangle )
    {
        i = triangle_node1[0+triangle*3];
        j = triangle_node1[1+triangle*3];
        k = triangle_node1[2+triangle*3];

        a = MIN ( i, j );
        b = MAX ( i, j );

        edge_data[0+5*(3*triangle+0)] = a;
        edge_data[1+5*(3*triangle+0)] = b;
        edge_data[2+5*(3*triangle+0)] = 1;
        edge_data[3+5*(3*triangle+0)] = 2;
        edge_data[4+5*(3*triangle+0)] = triangle;

        a = MIN ( i, k );
        b = MAX ( i, k );

        edge_data[0+5*(3*triangle+1)] = a;
        edge_data[1+5*(3*triangle+1)] = b;
        edge_data[2+5*(3*triangle+1)] = 1;
        edge_data[3+5*(3*triangle+1)] = 3;
        edge_data[4+5*(3*triangle+1)] = triangle;

        a = MIN ( j, k );
        b = MAX ( j, k );

        edge_data[0+5*(3*triangle+2)] = a;
        edge_data[1+5*(3*triangle+2)] = b;
        edge_data[2+5*(3*triangle+2)] = 2;
        edge_data[3+5*(3*triangle+2)] = 3;
        edge_data[4+5*(3*triangle+2)] = triangle;
    }
    /*
    Step 2. Perform an ascending dictionary sort on the neighbor relations.
    We only intend to sort on rows 1:2; the routine we call here
    sorts on the full column but that won't hurt us.

    What we need is to find all cases where triangles share an edge.
    By sorting the columns of the EDGE_DATA array, we will put shared edges
    next to each other.
    */
    i4col_sort_a ( 5, 3*triangle_num1, edge_data );
    /*
    Step 3. All the triangles which share an edge show up as consecutive
    columns with identical first two entries.  Figure out how many new
    nodes there are, and allocate space for their coordinates.
    */
    *node_num2 = node_num1;

    n1_old = n2_old = -1;

    for ( edge = 0; edge < 3 * triangle_num1; ++edge )
    {
        n1 = edge_data[0+edge*5];
        n2 = edge_data[1+edge*5];
        if ( n1 != n1_old || n2 != n2_old )
        {
            ++ *node_num2;
            n1_old = n1;
            n2_old = n2;
        }
    }

    *triangle_num2 = triangle_num1 << 2;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _triangulation_order6_adj_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
  Discussion:
    This routine is called to count the adjacencies, so that the
    appropriate amount of memory can be set aside for storage when
    the adjacency structure is created.
    The triangulation is assumed to involve 6-node triangles.
    Two nodes are "adjacent" if they are both nodes in some triangle.
    Also, a node is considered to be adjacent to itself.
  Diagram:
       3
    s  |\
    i  | \
    d  |  \
    e  6   5  side 2
       |    \
    3  |     \
       |      \
       1---4---2
         side 1
    The local node numbering
   21-22-23-24-25
    |\    |\    |
    | \   | \   |
   16 17 18 19 20
    |   \ |   \ |
    |    \|    \|
   11-12-13-14-15
    |\    |\    |
    | \   | \   |
    6  7  8  9 10
    |   \ |   \ |
    |    \|    \|
    1--2--3--4--5
    A sample grid.
    Below, we have a chart that lists the nodes adjacent to each node, with
    an asterisk to indicate the adjacency of the node to itself
 (in some cases, you want to count this self adjacency and in some
    you don't).
    N   Adjacencies
    1:  *  2  3  6  7 11
    2:  1  *  3  6  7 11
    3:  1  2  *  4  5  6  7  8  9 11 12 13
    4:  3  *  5  8  9 13
    5:  3  4  *  8  9 10 13 14 15
    6:  1  2  3  *  7 11
    7:  1  2  3  6  *  8 11 12 13
    8:  3  4  5  7  *  9 11 12 13
    9:  3  4  5  8  * 10 13 14 15
   10:  5  9  * 13 14 15
   11:  1  2  3  6  7  8  * 12 13 16 17 21
   12:  3  7  8 11  * 13 16 17 21
   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
   14:  5  9 10 13  * 15 18 19 23
   15:  5  9 10 13 14  * 18 19 20 23 24 25
   16: 11 12 13  * 17 21
   17: 11 12 13 16  * 18 21 22 23
   18: 13 14 15 17  * 19 21 22 23
   19: 13 14 15 18  * 20 23 24 25
   20: 15 19  * 23 24 25
   21: 11 12 13 16 17 18  * 22 23
   22: 13 17 18 21  * 23
   23: 13 14 15 17 18 19 20 21 22  * 24 25
   24: 15 19 20 23  * 25
   25: 15 19 20 23 24  *
    Below, we list the number of adjancencies to lower-indexed nodes, to
    the node itself, to higher-indexed nodes, the total number of
    adjacencies for this node, and the location of the first and last
    entries required to list this set of adjacencies in a single list
    of all the adjacencies.
    N   Below  Self   Above   Total First  Last
   --      --    --      --      --   ---     0
    1:      0     1       5       6     1     6
    2:      1     1       4       6     7    12
    3:      2     1       9      12    13    24
    4:      1     1       4       6    25    30
    5:      2     1       6       9    31    39
    6:      3     1       2       6    40    45
    7:      4     1       4       9    46    54
    8:      4     1       4       9    55    63
    9:      4     1       4       9    62    72
   10:      2     1       3       6    73    78
   11:      6     1       5      12    79    90
   12:      4     1       4       9    91    99
   13:      9     1       9      19   100   118
   14:      4     1       4       9   119   127
   15:      5     1       6      12   128   139
   16:      3     1       2       6   140   145
   17:      4     1       4       9   146   154
   18:      4     1       4       9   155   163
   19:      4     1       4       9   164   172
   20:      2     1       3       6   173   178
   21:      6     1       2       9   179   187
   22:      4     1       1       6   188   193
   23:      9     1       2      12   194   205
   24:      4     1       1       6   206   211
   25:      5     1       0       6   212   217
   --      --    --      --      --   218   ---
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 August 2006
  Author:
    John Burkardt
  Parameters
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
    make up each triangle.  The first three nodes are the vertices,
    in counterclockwise order.  The fourth value is the midside
    node between nodes 1 and 2; the fifth and sixth values are
    the other midside nodes in the logical order.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
    a triangle, lists the neighboring triangle, or -1 if there is
    no neighbor.
    Output, int TRIANGULATION_ORDER6_ADJ_COUNT, the number of adjacencies.
    Output, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dt3pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	int * triangle_neighbor = s_data->a3;
	int * adj_col = s_data->a4;
	
    dim_typ adj_num = 0;
    dim_typ i;
    int n1;
    int n2;
    int n3;
    int n4;
    int n5;
    int n6;
    dim_typ node;
    dim_typ triangle;
    dim_typ triangle_order = 6;
    dim_typ triangle2;

    /*
    Set every node to be adjacent to itself.
    */
    for ( node = 0; node < node_num; ++node )
        adj_col[node] = 1;
    /*
    Examine each triangle.
    */
    for ( triangle = 0; triangle < triangle_num; ++triangle )
    {
        n1 = triangle_node[0+triangle*triangle_order];
        n2 = triangle_node[1+triangle*triangle_order];
        n3 = triangle_node[2+triangle*triangle_order];
        n4 = triangle_node[3+triangle*triangle_order];
        n5 = triangle_node[4+triangle*triangle_order];
        n6 = triangle_node[5+triangle*triangle_order];
        /*
        For sure, we add the adjacencies:
        43 / (34)
        51 / (15)
        54 / (45)
        62 / (26)
        64 / (46)
        65 / (56)
        */
        ++ adj_col[n3-1];
        ++ adj_col[n4-1];
        ++ adj_col[n1-1];
        ++ adj_col[n5-1];
        ++ adj_col[n4-1];
        ++ adj_col[n5-1];
        ++ adj_col[n2-1];
        ++ adj_col[n6-1];
        ++ adj_col[n4-1];
        ++ adj_col[n6-1];
        ++ adj_col[n5-1];
        ++ adj_col[n6-1];
        /*
        Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
        that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
        or if this triangle is the first of the pair in which the edge
        occurs (TRIANGLE < TRIANGLE2).

        Maybe add
        21 / 12
        41 / 14
        42 / 24
        */
        triangle2 = triangle_neighbor[0+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ++ adj_col[n1-1];
            ++ adj_col[n2-1];
            ++ adj_col[n1-1];
            ++ adj_col[n4-1];
            ++ adj_col[n2-1];
            ++ adj_col[n4-1];
        }
        /*
        Maybe add
        32 / 23
        52 / 25
        53 / 35
        */
        triangle2 = triangle_neighbor[1+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ++ adj_col[n2-1];
            ++ adj_col[n3-1];
            ++ adj_col[n2-1];
            ++ adj_col[n5-1];
            ++ adj_col[n3-1];
            ++ adj_col[n5-1];
        }
        /*
        Maybe add
        31 / 13
        61 / 16
        63 / 36
        */
        triangle2 = triangle_neighbor[2+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            ++ adj_col[n1-1];
            ++ adj_col[n3-1];
            ++ adj_col[n1-1];
            ++ adj_col[n6-1];
            ++ adj_col[n3-1];
            ++ adj_col[n6-1];
        }
    }
    /*
    We used ADJ_COL to count the number of entries in each column.
    Convert it to pointers into the ADJ array.
    */
    for ( node = node_num; 1 <= node; --node )
        adj_col[node] = adj_col[node-1];
    adj_col[0] = 1;
    for ( i = 1; i <= node_num; ++i )
        adj_col[i] += adj_col[i-1];

	result = adj_col[node_num] - 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_order6_adj_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
  Discussion:
    This routine is called to count the adjacencies, so that the
    appropriate amount of memory can be set aside for storage when
    the adjacency structure is created.
    The triangulation is assumed to involve 6-node triangles.
    Two nodes are "adjacent" if they are both nodes in some triangle.
    Also, a node is considered to be adjacent to itself.
    This routine can be used to create the compressed column storage
    for a quadratic triangle finite element discretization of
    Poisson's equation in two dimensions.
  Diagram:
       3
    s  |\
    i  | \
    d  |  \
    e  6   5  side 2
       |    \
    3  |     \
       |      \
       1---4---2
         side 1
    The local node numbering
   21-22-23-24-25
    |\    |\    |
    | \   | \   |
   16 17 18 19 20
    |   \ |   \ |
    |    \|    \|
   11-12-13-14-15
    |\    |\    |
    | \   | \   |
    6  7  8  9 10
    |   \ |   \ |
    |    \|    \|
    1--2--3--4--5
    A sample grid.
    Below, we have a chart that lists the nodes adjacent to each node, with
    an asterisk to indicate the adjacency of the node to itself
 (in some cases, you want to count this self adjacency and in some
    you don't).
    N   Adjacencies
    1:  *  2  3  6  7 11
    2:  1  *  3  6  7 11
    3:  1  2  *  4  5  6  7  8  9 11 12 13
    4:  3  *  5  8  9 13
    5:  3  4  *  8  9 10 13 14 15
    6:  1  2  3  *  7 11
    7:  1  2  3  6  *  8 11 12 13
    8:  3  4  5  7  *  9 11 12 13
    9:  3  4  5  8  * 10 13 14 15
   10:  5  9  * 13 14 15
   11:  1  2  3  6  7  8  * 12 13 16 17 21
   12:  3  7  8 11  * 13 16 17 21
   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
   14:  5  9 10 13  * 15 18 19 23
   15:  5  9 10 13 14  * 18 19 20 23 24 25
   16: 11 12 13  * 17 21
   17: 11 12 13 16  * 18 21 22 23
   18: 13 14 15 17  * 19 21 22 23
   19: 13 14 15 18  * 20 23 24 25
   20: 15 19  * 23 24 25
   21: 11 12 13 16 17 18  * 22 23
   22: 13 17 18 21  * 23
   23: 13 14 15 17 18 19 20 21 22  * 24 25
   24: 15 19 20 23  * 25
   25: 15 19 20 23 24  *
    Below, we list the number of adjancencies to lower-indexed nodes, to
    the node itself, to higher-indexed nodes, the total number of
    adjacencies for this node, and the location of the first and last
    entries required to list this set of adjacencies in a single list
    of all the adjacencies.
    N   Below  Self   Above   Total First  Last
   --      --    --      --      --   ---     0
    1:      0     1       5       6     1     6
    2:      1     1       4       6     7    12
    3:      2     1       9      12    13    24
    4:      1     1       4       6    25    30
    5:      2     1       6       9    31    39
    6:      3     1       2       6    40    45
    7:      4     1       4       9    46    54
    8:      4     1       4       9    55    63
    9:      4     1       4       9    62    72
   10:      2     1       3       6    73    78
   11:      6     1       5      12    79    90
   12:      4     1       4       9    91    99
   13:      9     1       9      19   100   118
   14:      4     1       4       9   119   127
   15:      5     1       6      12   128   139
   16:      3     1       2       6   140   145
   17:      4     1       4       9   146   154
   18:      4     1       4       9   155   163
   19:      4     1       4       9   164   172
   20:      2     1       3       6   173   178
   21:      6     1       2       9   179   187
   22:      4     1       1       6   188   193
   23:      9     1       2      12   194   205
   24:      4     1       1       6   206   211
   25:      5     1       0       6   212   217
   --      --    --      --      --   218   ---
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 August 2006
  Author:
    John Burkardt
  Parameters
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
    make up each triangle.  The first three nodes are the vertices,
    in counterclockwise order.  The fourth value is the midside
    node between nodes 1 and 2; the fifth and sixth values are
    the other midside nodes in the logical order.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
    a triangle, lists the neighboring triangle, or -1 if there is
    no neighbor.
    Input, int ADJ_NUM, the number of adjacencies.
    Input, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
    Output, int TRIANGULATION_ORDER6_ADJ_SET[ADJ_NUM], the adjacency
    information.
*/
{
	const _2dt2pidtpi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	int * triangle_neighbor = s_data->a3;
	const register dim_typ adj_num = s_data->a4;
	int * adj_col = s_data->a5;
	
    int *adj;
    int *adj_copy;
    dim_typ k;
    int k1;
    int k2;
    int n1;
    int n2;
    int n3;
    int n4;
    int n5;
    int n6;
    dim_typ node;
    dim_typ triangle;
    dim_typ triangle2;
    dim_typ triangle_order = 6;

    adj = ( int * ) malloc ( adj_num * sizeof ( int ) );
    for ( k = 0; k < adj_num; ++k )
        adj[k] = -1;


    adj_copy = ( int * ) malloc ( node_num * sizeof ( int ) );
    for ( node = 0; node < node_num; ++node )

        adj_copy[node] = adj_col[node];
    /*
    Set every node to be adjacent to itself.
    */
    for ( node = 1; node <= node_num; ++node )
    {
        adj[adj_copy[node-1]-1] = node;
        ++ adj_copy[node-1];
    }
    /*
    Examine each triangle.
    */
    for ( triangle = 0; triangle < triangle_num; ++ triangle )
    {
        n1 = triangle_node[0+triangle*triangle_order];
        n2 = triangle_node[1+triangle*triangle_order];
        n3 = triangle_node[2+triangle*triangle_order];
        n4 = triangle_node[3+triangle*triangle_order];
        n5 = triangle_node[4+triangle*triangle_order];
        n6 = triangle_node[5+triangle*triangle_order];
        /*
        For sure, we add the adjacencies:
        43 / (34)
        51 / (15)
        54 / (45)
        62 / (26)
        64 / (46)
        65 / (56)
        */
        adj[adj_copy[n3-1]-1] = n4;
        ++ adj_copy[n3-1];
        adj[adj_copy[n4-1]-1] = n3;
        ++ adj_copy[n4-1];

        adj[adj_copy[n1-1]-1] = n5;
        ++ adj_copy[n1-1];
        adj[adj_copy[n5-1]-1] = n1;
        ++ adj_copy[n5-1];

        adj[adj_copy[n4-1]-1] = n5;
        ++ adj_copy[n4-1];
        adj[adj_copy[n5-1]-1] = n4;
        ++ adj_copy[n5-1];

        adj[adj_copy[n2-1]-1] = n6;
        ++ adj_copy[n2-1];
        adj[adj_copy[n6-1]-1] = n2;
        ++ adj_copy[n6-1];

        adj[adj_copy[n4-1]-1] = n6;
        ++ adj_copy[n4-1];
        adj[adj_copy[n6-1]-1] = n4;
        ++ adj_copy[n6-1];

        adj[adj_copy[n5-1]-1] = n6;
        ++adj_copy[n5-1];
        adj[adj_copy[n6-1]-1] = n5;
        ++ adj_copy[n6-1];
        /*
        Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
        that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
        or if this triangle is the first of the pair in which the edge
        occurs (TRIANGLE < TRIANGLE2).

        Maybe add
        21 / 12
        41 / 14
        42 / 24
        */
        triangle2 = triangle_neighbor[0+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj[adj_copy[n1-1]-1] = n2;
            ++ adj_copy[n1-1];
            adj[adj_copy[n2-1]-1] = n1;
            ++ adj_copy[n2-1];
            adj[adj_copy[n1-1]-1] = n4;
            ++ adj_copy[n1-1];
            adj[adj_copy[n4-1]-1] = n1;
            ++ adj_copy[n4-1];
            adj[adj_copy[n2-1]-1] = n4;
            ++ adj_copy[n2-1];
            adj[adj_copy[n4-1]-1] = n2;
            ++ adj_copy[n4-1];
        }
        /*
        Maybe add
        32 / 23
        52 / 25
        53 / 35
        */
        triangle2 = triangle_neighbor[1+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj[adj_copy[n2-1]-1] = n3;
            ++ adj_copy[n2-1];
            adj[adj_copy[n3-1]-1] = n2;
            ++ adj_copy[n3-1];
            adj[adj_copy[n2-1]-1] = n5;
            ++ adj_copy[n2-1];
            adj[adj_copy[n5-1]-1] = n2;
            ++ adj_copy[n5-1];
            adj[adj_copy[n3-1]-1] = n5;
            ++ adj_copy[n3-1];
            adj[adj_copy[n5-1]-1] = n3;
            ++ adj_copy[n5-1];
        }
        /*
        Maybe add
        31 / 13
        61 / 16
        63 / 36
        */
        triangle2 = triangle_neighbor[2+triangle*3];

        if ( triangle2 < 0 || triangle < triangle2 )
        {
            adj[adj_copy[n1-1]-1] = n3;
            ++ adj_copy[n1-1];
            adj[adj_copy[n3-1]-1] = n1;
            ++ adj_copy[n3-1];
            adj[adj_copy[n1-1]-1] = n6;
            ++ adj_copy[n1-1];
            adj[adj_copy[n6-1]-1] = n1;
            ++ adj_copy[n6-1];
            adj[adj_copy[n3-1]-1] = n6;
            ++ adj_copy[n3-1];
            adj[adj_copy[n6-1]-1] = n3;
            ++ adj_copy[n6-1];
        }
    }
    /*
    Ascending sort the entries for each node.
    */
    for ( node = 1; node <= node_num; ++node )
    {
        k1 = adj_col[node-1];
        k2 = adj_col[node]-1;
        i4vec_sort_heap_a ( k2+1-k1, adj+k1-1 );
    }

    free ( adj_copy );

    return adj;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order6_boundary_edge_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the boundary edges.
  Discussion:
    This routine is given a triangulation, a set of 6-node triangles.
    It is assumed that, in each list of 6 nodes, the vertices are listed
    first, in counterclockwise order, followed by the three midside nodes,
    in counterclockwise order, starting with the node between vertices
    1 and 2.
    It is assumed that each edge of the triangulation is either
    * an INTERIOR edge, which is listed twice, once with positive
      orientation and once with negative orientation, or;
    * a BOUNDARY edge, which will occur only once.
    This routine should work even if the region has holes - as long
    as the boundary of the hole comprises more than 3 edges!
    Except for the dimension of TRIANGLE, this routine is identical
    to the routine for the order 3 case.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 June 2005
  Author:
    John Burkardt
  Parameters:
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up the
    triangles.  These should be listed in counterclockwise order.
    Output, integer TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT, the number
    of boundary edges.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ triangle_num = s_data->a0;
	int * triangle_node = s_data->a1;
	
    dim_typ boundary_edge_num;
    int e1;
    int e2;
    int *edge;
    dim_typ i;
    dim_typ interior_edge_num;
    dim_typ j;
    dim_typ m;
    dim_typ n;
    dim_typ unique_num;

    m = 2;
    n = 3 * triangle_num;
    /*
    Set up the edge array.
    */
    edge = ( int * ) malloc ( m * n * sizeof ( int ) );

    for ( j = 0; j < triangle_num; ++j )
    {
        edge[0+(j               )*m] = triangle_node[0+j*6];
        edge[1+(j               )*m] = triangle_node[1+j*6];
        edge[0+(j+  triangle_num)*m] = triangle_node[1+j*6];
        edge[1+(j+  triangle_num)*m] = triangle_node[2+j*6];
        edge[0+(j+(triangle_num<<1))*m] = triangle_node[2+j*6];
        edge[1+(j+(triangle_num<<1))*m] = triangle_node[0+j*6];
    }
    /*
    In each column, force the smaller entry to appear first.
    */
    for ( j = 0; j < n; ++j )
    {
        e1 = MIN ( edge[0+j*m], edge[1+j*m] );
        e2 = MAX ( edge[0+j*m], edge[1+j*m] );
        edge[0+j*m] = e1;
        edge[1+j*m] = e2;
    }
    /*
    Ascending sort the column array.
    */
    i4col_sort_a ( m, n, edge );
    /*
    Get the number of unique columns in EDGE.
    */
    unique_num = i4col_sorted_unique_count ( m, n, edge );

    interior_edge_num = 3 * triangle_num - unique_num;
    boundary_edge_num = 3 * triangle_num - (interior_edge_num<<1);

    free ( edge );

	result = boundary_edge_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _triangulation_order6_boundary_edge_count_euler ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
  Discussion:
    We assume we are given information about an order 6 triangulation
    of a set of nodes in the plane.
    By ignoring the midside nodes, we can determine the corresponding
    information for an order 3 triangulation, and apply
    Euler's formula to determine the number of edges that lie on the
    boundary of the set of nodes.
    Thus, if we have TRIANGLE_NUM triangles, and NODE_NUM nodes, we
    imagine that each triangle is replaced by 4 triangles, created
    by adding the edges created by joining the midside nodes.
    Thus, for 4 * TRIANGLE_NUM triangles, we can apply Euler's formula
    to compute the number of boundary edges.
    Now, to adjust the data to our order 6 triangles, we divide the
    number of boundary edges by 2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 June 2005
  Author:
    John Burkardt
  Reference:
    Marc deBerg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
    Computational Geometry,
    Springer, 2000,
    ISBN: 3-540-65620-0.
  Parameters:
    Input, integer NODE_NUM, the number of nodes.
    Input, integer TRIANGLE_NUM, the number of triangles.
    Input, integer HOLE_NUM, the number of internal nodes.
    Output, int TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT, the number of
    edges that lie on the boundary of the triangulation.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data; 
	const register dim_typ node_num = a_data[0];
	const register dim_typ triangle_num = a_data[1];
	const register dim_typ hole_num = a_data[2];
	
    result = ( (node_num<<1) + (hole_num<<1) - (triangle_num<<2) - 2 ) >> 1;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_order6_boundary_node ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_BOUNDARY_NODE indicates nodes on the boundary.
  Discussion:
    This routine is given an order 6 triangulation, an abstract list of
    sets of six nodes.  The vertices are listed clockwise, then the
    midside nodes.
    It is assumed that each edge of the triangulation is either
    * an INTERIOR edge, which is listed twice, once with positive
      orientation and once with negative orientation, or;
    * a BOUNDARY edge, which will occur only once.
    This routine should work even if the region has holes - as long
    as the boundary of the hole comprises more than 3 edges!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 January 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up the
    triangles.
    Output, int TRIANGULATION_ORDER6_BOUNDARY_NODE[NODE_NUM],
    is TRUE if the node is on a boundary edge.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	
    int e1;
    int e2;
    int *edge;
    bool equal;
    dim_typ i;
    dim_typ j;
    dim_typ m;
    dim_typ n;
    int *node_boundary;

    m = 3;
    n = 3 * triangle_num;
    /*
    Set up the edge array.
    */
    edge = ( int * ) malloc ( m * n * sizeof ( int ) );

    for ( j = 0; j < triangle_num; ++j )
    {
        edge[0+(j               )*m] = triangle_node[0+j*6];
        edge[1+(j               )*m] = triangle_node[3+j*6];
        edge[2+(j               )*m] = triangle_node[1+j*6];

        edge[0+(j+  triangle_num)*m] = triangle_node[1+j*6];
        edge[1+(j+  triangle_num)*m] = triangle_node[4+j*6];
        edge[2+(j+  triangle_num)*m] = triangle_node[2+j*6];

        edge[0+(j+(triangle_num<<1))*m] = triangle_node[2+j*6];
        edge[1+(j+(triangle_num<<1))*m] = triangle_node[5+j*6];
        edge[2+(j+(triangle_num<<1))*m] = triangle_node[0+j*6];
    }
    /*
    In each column, force the smaller entry to appear first.
    */
    for ( j = 0; j < n; j++ )
    {
        e1 = MIN ( edge[0+j*m], edge[2+j*m] );
        e2 = MAX ( edge[0+j*m], edge[2+j*m] );
        edge[0+j*m] = e1;
        edge[2+j*m] = e2;
    }
    /*
    Ascending sort the column array.
    */
    i4col_sort_a ( m, n, edge );
    /*
    Records which appear twice are internal edges and can be ignored.
    */
    node_boundary = ( int * ) malloc ( node_num * sizeof ( int ) );

    for ( i = 0; i < node_num; ++i )
        node_boundary[i] = 0;


    j = 0;

    while ( j < 3 * triangle_num )
    {
        ++ j;

        if ( j == 3 * triangle_num )
        {
            for ( i = 0; i < m; ++i )
                node_boundary[edge[i+(j-1)*m]-1] = 1;
            break;
        }

        equal = true;

        for ( i = 0; i < m; ++i )
            if ( edge[i+(j-1)*m] != edge[i+j*m] )
                equal = false;


        if ( equal )
            ++ j;
        else
            for ( i = 0; i < m; ++i )
                node_boundary[edge[i+(j-1)*m]-1] = 1;
    }

    free ( edge );

    return node_boundary;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _triangulation_order6_refine_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_REFINE_COMPUTE computes a refined order 6 triangulation.
  Discussion:
    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
    quadratic subtriangles T1, T2, T3 and T4.
    The task is more complicated by the fact that we are working with
    a mesh of triangles, so that we want to create a node only once,
    even though it may be shared by other triangles. (In fact, only
    the new nodes on the edges can be shared, and then only by at most
    one other triangle.)
            3
           / \
          36 35
         / T3  \
        6--56---5
       / \ T4  / \
      16 46  45  25
     / T1  \ / T2  \
    1--14---4--24---2
    This routine is given sorted information defining the edges, and uses
    it to build the new node and triangle arrays.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM1, the number of nodes.
    Input, int TRIANGLE_NUM1, the number of triangles.
    Input, double NODE_XY1[2*NODE_NUM1], the nodes.
    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the nodes that make up the
    triangles.
    Input, int NODE_NUM2, the number of nodes in the refined mesh.
    Input, int TRIANGLE_NUM2, the number of triangles in the refined mesh.
    Input, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge information computed
    by TRIANGULATION_ORDER6_REFINE_SIZE.
    Output, double NODE_XY2[2*NODE_NUM2], the refined nodes.
    Output, int TRIANGLE_NODE2[6*TRIANGLE_NUM2], the nodes that make up the
    triangles in the refined mesh.
*/
{
	const _2dtpitpi2dtpipitpi * const s_data = data;
	const register dim_typ node_num1 = s_data->a0;
	const register dim_typ triangle_num1 = s_data->a1;
	ityp * node_xy1 = s_data->a2;
	int * triangle_node1 = s_data->a3;
	const register dim_typ node_num2 = s_data->a4;
	const register dim_typ triangle_num2 = s_data->a5;
	int * edge_data = s_data->a6;
	ityp * node_xy2 = s_data->a7;
	int * triangle_node2 = s_data->a8;
	
    dim_typ edge;
    int i;
    int j;
    int l1;
    int l2;
    int l3;
    int n1;
    int n1_old;
    int n2;
    int n2_old;
    dim_typ node;
    int t1;
    int t2;
    int t3;
    int t4;
    dim_typ triangle1;
    int v1;
    int v2;
    int v3;
    int v4;
    int v5;
    int v6;
    /*
    Step 1:
    Copy the old nodes.
    */
    for ( j = 0; j < node_num1; ++j )
        #pragma omp parallel for num_threads(2)
        for ( i = 0; i < 2; ++i )
            node_xy2[i+(j<<1)] = node_xy1[i+(j<<1)];
    for ( j = 0; j < triangle_num2; ++j )
        #pragma omp parallel for num_threads(6)
        for ( i = 0; i < 6; ++i )
            triangle_node2[i+j*6] = -1;
    /*
    We can assign the existing nodes to the new triangles.
    */
    for ( triangle1 = 0; triangle1 < triangle_num1; ++triangle1 )
    {
        t1 = (triangle1<<2) + 0;
        t2 = (triangle1<<2) + 1;
        t3 = (triangle1<<2) + 2;
        t4 = (triangle1<<2) + 3;

        triangle_node2[0+t1*6] = triangle_node1[0+triangle1*6];
        triangle_node2[1+t1*6] = triangle_node1[3+triangle1*6];
        triangle_node2[2+t1*6] = triangle_node1[5+triangle1*6];

        triangle_node2[0+t2*6] = triangle_node1[3+triangle1*6];
        triangle_node2[1+t2*6] = triangle_node1[1+triangle1*6];
        triangle_node2[2+t2*6] = triangle_node1[4+triangle1*6];

        triangle_node2[0+t3*6] = triangle_node1[5+triangle1*6];
        triangle_node2[1+t3*6] = triangle_node1[4+triangle1*6];
        triangle_node2[2+t3*6] = triangle_node1[2+triangle1*6];

        triangle_node2[0+t4*6] = triangle_node1[4+triangle1*6];
        triangle_node2[1+t4*6] = triangle_node1[5+triangle1*6];
        triangle_node2[2+t4*6] = triangle_node1[3+triangle1*6];
    }
    /*
    Step 2.
    Examine sorted edge information.  The first time an edge is encountered,
    generate two new nodes, then assign them (usually) to the four subtriangles
    of the two triangles that share that edge.
    */
    node = node_num1;

    n1_old = n2_old = -1;

    for ( edge = 0; edge < 3 * triangle_num1; ++edge )
    {
        n1 = edge_data[0+edge*5] - 1;
        n2 = edge_data[1+edge*5] - 1;

        l1 = edge_data[2+edge*5];
        l3 = edge_data[3+edge*5];

        if ( l1 == 1 && l3 == 2 )
            l2 = 4;
        else if ( l1 == 1 && l3 == 3 )
            l2 = 6;
        else if ( l1 == 2 && l3 == 3 )
            l2 = 5;
        triangle1 = edge_data[4+edge*5];
        /*
        If this is the first time we've encountered this edge,
        create the new nodes.
        */
        if ( n1 != n1_old || n2 != n2_old )
        {
            n1_old = n1;
            n2_old = n2;

            v1 = triangle_node1[l1-1+triangle1*6];
            v2 = triangle_node1[l2-1+triangle1*6];
            v3 = triangle_node1[l3-1+triangle1*6];

            #pragma omp parallel for num_threads(2)
            for ( i = 0; i < 2; ++i )
                node_xy2[i+(node<<1)] = ( node_xy2[i+((v1-1)<<1)]+ node_xy2[i+((v2-1)<<1)] ) / 2.00;
            ++ node;
            v4 = node;

            #pragma omp parallel for num_threads(2)
            for ( i = 0; i < 2; ++i )
                node_xy2[i+(node<<1)] = ( node_xy2[i+((v2-1)<<1)]+ node_xy2[i+((v3-1)<<1)] ) / 2.00;
            ++ node;
            v5 = node;
        }
        t1 = (triangle1 << 2) + 0;
        t2 = (triangle1 << 2) + 1;
        t3 = (triangle1 << 2) + 2;

        if ( l1 == 1 && l3 == 2 )
        {
            if ( triangle_node1[0+triangle1*6] == v1 + 1 )
            {
                triangle_node2[3+t1*6] = v4;
                triangle_node2[3+t2*6] = v5;
            }
            else
            {
                triangle_node2[3+t1*6] = v5;
                triangle_node2[3+t2*6] = v4;
            }
        }
        else if ( l1 == 1 && l3 == 3 )
        {
            if ( triangle_node1[0+triangle1*6] == v1 + 1 )
            {
                triangle_node2[5+t1*6] = v4;
                triangle_node2[5+t3*6] = v5;
            }
            else
            {
                triangle_node2[5+t1*6] = v5;
                triangle_node2[5+t3*6] = v4;
            }
        }
        else if ( l1 == 2 && l3 == 3 )
        {
            if ( triangle_node1[1+triangle1*6] == v1 + 1 )
            {
                triangle_node2[4+t3*6] = v4;
                triangle_node2[4+t2*6] = v5;
            }
            else
            {
                triangle_node2[4+t3*6] = v5;
                triangle_node2[4+t2*6] = v4;
            }
        }
    }
    /*
    Step 3.
    Each old triangle has a single central subtriangle, for which we now
    need to generate three new "interior" nodes.
    */
    for ( triangle1 = 0; triangle1 < triangle_num1; ++triangle1)
    {
        v4 = triangle_node1[3+triangle1*6];
        v5 = triangle_node1[4+triangle1*6];
        v6 = triangle_node1[5+triangle1*6];

        t1 = (triangle1 << 2) + 0;
        t2 = (triangle1 << 2) + 1;
        t3 = (triangle1 << 2) + 2;
        t4 = (triangle1 << 2) + 3;

        node_xy2[0+(node<<1)] = 0.50 * ( node_xy1[0+((v5-1)<<1)] + node_xy1[0+((v6-1)<<1)] );
        node_xy2[1+(node<<1)] = 0.50 * ( node_xy1[1+((v5-1)<<1)] + node_xy1[1+((v6-1)<<1)] );
        ++ node ;
        triangle_node2[3+t4*6] = node;
        triangle_node2[3+t3*6] = node;

        node_xy2[0+(node<<1)] = 0.5 * ( node_xy1[0+((v6-1)<<1)] + node_xy1[0+((v4-1)<<1)] );
        node_xy2[1+(node<<1)] = 0.5 * ( node_xy1[1+((v6-1)<<1)] + node_xy1[1+((v4-1)<<1)] );
        node = node + 1;
        triangle_node2[4+t4*6] = node;
        triangle_node2[4+t1*6] = node;

        node_xy2[0+(node<<1)] = 0.5 * ( node_xy1[0+((v4-1)<<1)] + node_xy1[0+((v5-1)<<1)] );
        node_xy2[1+(node<<1)] = 0.5 * ( node_xy1[1+((v4-1)<<1)] + node_xy1[1+((v5-1)<<1)] );
        ++ node;
        triangle_node2[5+t4*6] = node;
        triangle_node2[5+t2*6] = node;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order6_refine_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_REFINE_SIZE sizes a refined order 6 triangulation.
  Discussion:
    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
    quadratic subtriangles T1, T2, T3 and T4.
    The task is more complicated by the fact that we are working with
    a mesh of triangles, so that we want to create a node only once,
    even though it may be shared by other triangles. (In fact, only
    the new nodes on the edges can be shared, and then only by at most
    one other triangle.)
            3
           / \
          36 35
         / T3  \
        6--56---5
       / \ T4  / \
      16 46  45  25
     / T1  \ / T2  \
    1--14---4--24---2
    This routine determines the sizes of the resulting node and
    triangles, and constructs an edge array that can be used to
    properly number the new nodes.
    The primary work occurs in sorting a list related to the edges.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2007
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM1, the number of nodes in the original mesh.
    Input, int  TRIANGLE_NUM1, the number of triangles in the
    original mesh.
    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the indices of the nodes
    that form the triangles in the input mesh.
    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
    Output, int *TRIANGLE_NUM2, the number of triangles in the
    refined mesh.
    Output, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge data that will
    be needed by TRIANGULATION_ORDER6_REFINE_COMPUTE.
*/
{
	const _2dt4pi * const s_data = data;
	register dim_typ node_num1 = s_data->a0;
	const register dim_typ triangle_num1 = s_data->a1;
	int * triangle_node1 = s_data->a2;
	int * node_num2 = s_data->a3;
	int * triangle_num2 = s_data->a4;
	int * edge_data = s_data->a5;
	
    int a;
    int b;
    dim_typ edge;
    dim_typ i, j, k;
    int n1;
    int n1_old;
    int n2;
    int n2_old;
    dim_typ triangle1;
    /*
    Step 1.
    From the list of nodes for triangle T, of the form: (I,J,K)
    construct the edge relations:

 (I,J,1,2,T)
 (I,K,1,3,T)
 (J,K,2,3,T)

    In order to make matching easier, we reorder each pair of nodes
    into ascending order.
    */
    for ( triangle1 = 0; triangle1 < triangle_num1; ++triangle1 )
    {
        i = triangle_node1[0+triangle1*6];
        j = triangle_node1[1+triangle1*6];
        k = triangle_node1[2+triangle1*6];

        a = MIN ( i, j );
        b = MAX ( i, j );

        edge_data[0+5*(3*triangle1+0)] = a;
        edge_data[1+5*(3*triangle1+0)] = b;
        edge_data[2+5*(3*triangle1+0)] = 1;
        edge_data[3+5*(3*triangle1+0)] = 2;
        edge_data[4+5*(3*triangle1+0)] = triangle1;

        a = MIN ( i, k );
        b = MAX ( i, k );

        edge_data[0+5*(3*triangle1+1)] = a;
        edge_data[1+5*(3*triangle1+1)] = b;
        edge_data[2+5*(3*triangle1+1)] = 1;
        edge_data[3+5*(3*triangle1+1)] = 3;
        edge_data[4+5*(3*triangle1+1)] = triangle1;

        a = MIN ( j, k );
        b = MAX ( j, k );

        edge_data[0+5*(3*triangle1+2)] = a;
        edge_data[1+5*(3*triangle1+2)] = b;
        edge_data[2+5*(3*triangle1+2)] = 2;
        edge_data[3+5*(3*triangle1+2)] = 3;
        edge_data[4+5*(3*triangle1+2)] = triangle1;
    }
    /*
    Step 2. Perform an ascending dictionary sort on the neighbor relations.
    We only intend to sort on rows 1:2; the routine we call here
    sorts on the full column but that won't hurt us.

    What we need is to find all cases where triangles share an edge.
    By sorting the columns of the EDGE_DATA array, we will put shared edges
    next to each other.
    */
    i4col_sort_a ( 5, 3*triangle_num1, edge_data );
    /*
    Step 3. All the triangles which share an edge show up as consecutive
    columns with identical first two entries.  Figure out how many new
    nodes there are, and allocate space for their coordinates.
    */
    *node_num2 = node_num1;

    n1_old = n2_old = -1;

    for ( edge = 0; edge < 3 * triangle_num1; ++edge)
    {
        n1 = edge_data[0+edge*5];
        n2 = edge_data[1+edge*5];
        if ( n1 != n1_old || n2 != n2_old )
        {
            *node_num2 += 2;
            n1_old = n1;
            n2_old = n2;
        }
    }

    *node_num2 += 3 * triangle_num1;

    *triangle_num2 = triangle_num1 << 2;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_order6_to_order3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_TO_ORDER3 linearizes a quadratic triangulation.
  Discussion:
    A quadratic triangulation is assumed to consist of 6-node triangles,
    as in the following:
    11-12-13-14-15
     |\    |\    |
     | \   | \   |
     6  7  8  9 10
     |   \ |   \ |
     |    \|    \|
     1--2--3--4--5
   This routine rearranges information so as to define the 3-node
   triangulation:
    11-12-13-14-15
     |\ |\ |\ |\ |
     | \| \| \| \|
     6--7--8--9-10
     |\ |\ |\ |\ |
     | \| \| \| \|
     1--2--3--4--5
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 March 2005
  Author:
    John Burkardt
  Parameters:
    Input, int TRIANGLE_NUM1, the number of triangles in the quadratic
    triangulation.
    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the quadratic triangulation.
    Output, int TRIANGULATION_ORDER6_TO_ORDER3[3*TRIANGLE_NUM2], the linear
    triangulation.  Here, TRIANGLE_NUM2 = 4 * TRIANGLE_NUM1.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ triangle_num1 = s_data->a0;
	int * triangle_node1 = s_data->a1;
	
    int n1;
    int n2;
    int n3;
    int n4;
    int n5;
    int n6;
    dim_typ triangle_num2;
    dim_typ tri1;
    dim_typ tri2;
    int *triangle_node2;

    triangle_num2 = triangle_num1<<2;
    triangle_node2 = ( int * ) malloc ( 3 * triangle_num2 * sizeof ( int ) );

    tri2 = 0;

    for ( tri1 = 0; tri1 < triangle_num1; ++tri1 )
    {
        n1 = triangle_node1[0+tri1*6];
        n2 = triangle_node1[1+tri1*6];
        n3 = triangle_node1[2+tri1*6];
        n4 = triangle_node1[3+tri1*6];
        n5 = triangle_node1[4+tri1*6];
        n6 = triangle_node1[5+tri1*6];

        triangle_node2[0+tri2*3] = n1;
        triangle_node2[1+tri2*3] = n4;
        triangle_node2[2+tri2*3] = n6;
        ++ tri2;

        triangle_node2[0+tri2*3] = n2;
        triangle_node2[1+tri2*3] = n5;
        triangle_node2[2+tri2*3] = n4;
        ++ tri2;

        triangle_node2[0+tri2*3] = n3;
        triangle_node2[1+tri2*3] = n6;
        triangle_node2[2+tri2*3] = n5;
        ++ tri2;

        triangle_node2[0+tri2*3] = n4;
        triangle_node2[1+tri2*3] = n5;
        triangle_node2[2+tri2*3] = n6;
        ++ tri2;
    }

    return triangle_node2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_order6_vertex_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER6_VERTEX_COUNT counts vertex nodes in a triangulation.
  Discussion:
    In a triangulation of order 6, some nodes are midside nodes and some
    nodes are vertex nodes.
    Especially when an order 6 triangulation is used to handle the
    Navier Stokes equations, it is useful to know the number of
    vertex and midside nodes.
    Note that the number of midside nodes is simple NODE_NUM - VERTEX_NUM.
  Diagram:
       3
    s  |\
    i  | \
    d  |  \
    e  6   5  side 2
       |    \
    3  |     \
       |      \
       1---4---2
         side 1
    The local node numbering.  Local nodes 1, 2 and 3 are "vertex nodes",
    while nodes 4, 5 and 6 are "midside nodes".
   21-22-23-24-25
    |\    |\    |
    | \   | \   |
   16 17 18 19 20
    |   \ |   \ |
    |    \|    \|
   11-12-13-14-15
    |\    |\    |
    | \   | \   |
    6  7  8  9 10
    |   \ |   \ |
    |    \|    \|
    1--2--3--4--5
    A sample grid, which contains 9 vertex nodes and 16 midside nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 August 2006
  Author:
    John Burkardt
  Parameters
    Input, int TRI_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[6*TRI_NUM], lists the nodes that
    make up each triangle.  The first three nodes are the vertices,
    in counterclockwise order.  The fourth value is the midside
    node between nodes 1 and 2; the fifth and sixth values are
    the other midside nodes in the logical order.
    Output, int TRIANGULATION_ORDER6_VERTEX_COUNT, the number of nodes
    which are vertices.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ tri_num = s_data->a0;
	int * triangle_node = s_data->a1;
	
    dim_typ j;
    dim_typ vertex_num;
    int *vertices;

    vertices = ( int * ) malloc ( 3 * tri_num * sizeof ( int ) );

    for ( j = 0; j < tri_num; ++j )
    {
        vertices[j] = triangle_node[0+j*6];
        vertices[tri_num+j] = triangle_node[1+j*6];
        vertices[(tri_num<<1)+j] = triangle_node[2+j*6];
    }

    i4vec_sort_heap_a ( 3*tri_num, vertices );
    vertex_num = i4vec_sorted_unique ( 3*tri_num, vertices );
    free ( vertices );

	result = vertex_num; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangulation_search_naive ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_SEARCH_NAIVE naively searches a triangulation for a point.
  Discussion:
    The algorithm simply checks each triangle to see if point P is
    contained in it.  Surprisingly, this is not the fastest way to
    do the check, at least if the triangulation is Delaunay.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Input, int TRIANGLE_ORDER, the order of the triangles.
    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
    the nodes of each triangle.
    Input, double P[2], the coordinates of a point.
    Output, int TRIANGULATION_SEARCH_NAIVE, the index of the triangle
    containing the point, or -1 if no triangle was found containing
    the point.
*/
{
	static short result = SHRT_MAX;
	
	const dtpit2dtpipit * const s_data = data;
	const register dim_typ node_num = s_data->a0; 
	ityp * node_xy = s_data->a1;
	const register dim_typ triangle_order = s_data->a2;
	const register dim_typ triangle_num = s_data->a3;
	int * triangle_node = s_data->a4; 
	ityp * p = s_data->a5;
	
    int a;
    ityp alpha;
    int b;
    ityp beta;
    int c;
    ityp det;
    ityp dxp;
    ityp dxa;
    ityp dxb;
    ityp dyp;
    ityp dya;
    ityp dyb;
    ityp gamma;
    dim_typ triangle;
    short triangle_index = -1;

    for ( triangle = 0; triangle < triangle_num; ++ triangle )
    {
        /*
        Get the vertices of triangle TRIANGLE.
        */
        a = triangle_node[0+triangle*triangle_order];
        b = triangle_node[1+triangle*triangle_order];
        c = triangle_node[2+triangle*triangle_order];
        /*
        Using vertex C as a base, compute the distances to vertices A and B,
        and the point (X,Y).
        */
        dxa = node_xy[0+a*2] - node_xy[0+c*2];
        dya = node_xy[1+a*2] - node_xy[1+c*2];

        dxb = node_xy[0+b*2] - node_xy[0+c*2];
        dyb = node_xy[1+b*2] - node_xy[1+c*2];

        dxp = p[0]           - node_xy[0+c*2];
        dyp = p[1]           - node_xy[1+c*2];

        det = dxa * dyb - dya * dxb;
        /*
        Compute the barycentric coordinates of the point (X,Y) with respect
        to this triangle.
        */
        alpha = ( dxp * dyb - dyp * dxb ) / det;
        beta = ( dxa * dyp - dya * dxp ) / det;
        gamma = 1.00 - alpha - beta;
        /*
        If the barycentric coordinates are all positive, then the point
        is inside the triangle and we're done.
        */
        if ( 0.00 <= alpha &&0.00 <= beta  &&0.00 <= gamma )
        {
            triangle_index = triangle + 1;
            break;
        }
    }

	result = triangle_index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _voronoi_polygon_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    VORONOI_POLYGON_AREA computes the area of a Voronoi polygon.
  Formula:
    It is assumed that the Voronoi polygon is finite!  Every Voronoi
    diagram includes some regions which are infinite, and for those,
    this formula is not appropriate.
    The routine is given the indices of the nodes that are
    Voronoi "neighbors" of a given node.  These are also the nodes
    that are paired to form edges of Delaunay triangles.
    The assumption that the polygon is a Voronoi polygon is
    used to determine the location of the boundaries of the polygon,
    which are the perpendicular bisectors of the lines connecting
    the center point to each of its neighbors.
    The finiteness assumption is employed in part in the
    assumption that the polygon is bounded by the finite
    line segments from point 1 to 2, 2 to 3, ...,
    M-1 to M, and M to 1, where M is the number of neighbors.
    It is assumed that this routine is being called by a
    process which has computed the Voronoi diagram of a large
    set of nodes, so the arrays X and Y are dimensioned by
    NODE_NUM, which may be much greater than the number of neighbor
    nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 February 2005
  Author:
    John Burkardt
  Reference:
    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
    Second Edition,
    Wiley, 2000, page 485.
  Parameters:
    Input, int NODE, the index of the node whose Voronoi
    polygon is to be measured. 0 <= NODE < NODE_NUM.
    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
    the given node.
    Input, int NEIGHBOR_INDEX[NEIGHBOR_NUM], the indices
    of the neighbor nodes (used to access X and Y).  The neighbor
    nodes should be listed in the (counter-clockwise) order in
    which they occur as one circles the center node.
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Output, double VORONOI_POLYGON_AREA, the area of the Voronoi polygon.
*/
{
	static ityp result = MAX_VAL;
	
	const _3dtpipit * const s_data = data; 
	
	const register dim_typ node = s_data->a0;
	const register dim_typ neighbor_num = s_data->a1;
	const register dim_typ node_num = s_data->a2;
	int * neighbor_index = s_data->a3;
	ityp * node_xy = s_data->a4;
	
    ityp a;
    ityp area;
    ityp b;
    ityp c;
    dim_typ i;
    dim_typ ip1;
    ityp ui;
    ityp uip1;
    ityp vi;
    ityp vip1;
    ityp xc;
    ityp xi;
    ityp xip1;
    ityp yc;
    ityp yi;
    ityp yip1;

    if ( node_num <= node )
    {
    	result = MAX_VAL;
        return &result;
    }

    area = 0.00;

    xc = node_xy[0+(node<<1)];
    yc = node_xy[1+(node<<1)];

    i = neighbor_num - 1;
    i = neighbor_index[i];

    xi = node_xy[i<<1];
    yi = node_xy[1+(i<<1)];

    ip1 = 0;
    ip1 = neighbor_index[ip1];

    xip1 = node_xy[0+(ip1<<1)];
    yip1 = node_xy[1+(ip1<<1)];
    a = ( xi   * xi   + yi   * yi   - xc * xc - yc * yc );
    b = ( xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc );
    c = 2.00 * ( ( xi - xc ) * ( yip1 - yc ) - ( xip1 - xc ) * ( yi - yc ) );
    uip1 = ( a * ( yip1 - yc ) - b * ( yi - yc )  ) / c;
    vip1 = ( a * ( xip1 - xc ) - b * ( xi - xc )  ) / c;

    for ( i = 0; i < neighbor_num; ++i)
    {
        xi = xip1;
        yi = yip1;
        ui = uip1;
        vi = vip1;

        ip1 = i4_wrap ( i+1, 0, neighbor_num-1 );
        ip1 = neighbor_index[ip1];

        xip1 = node_xy[0+(ip1<<1)];
        yip1 = node_xy[1+(ip1<<1)];
        a = ( xi   * xi   + yi   * yi   - xc * xc - yc * yc );
        b = ( xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc );
        c = 2.00 * ( ( xi - xc ) * ( yip1 - yc ) - ( xip1 - xc ) * ( yi - yc ) );
        uip1 = ( a * ( yip1 - yc ) - b * ( yi - yc )  ) / c;
        vip1 = ( a * ( xip1 - xc ) - b * ( xi - xc )  ) / c;

        area += uip1 * vi - ui * vip1;
    }

	result = 0.50 * area;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _voronoi_polygon_centroid ( void * data)
/******************************************************************************/
/*
  Purpose:
    VORONOI_POLYGON_CENTROID_2D computes the centroid of a Voronoi polygon.
  Formula:
    It is assumed that the Voronoi polygon is finite!  Every Voronoi
    diagram includes some regions which are infinite, and for those,
    this formula is not appropriate.
    The routine is given the indices of the nodes that are
    Voronoi "neighbors" of a given node.  These are also the nodes
    that are paired to form edges of Delaunay triangles.
    The assumption that the polygon is a Voronoi polygon is
    used to determine the location of the boundaries of the polygon,
    which are the perpendicular bisectors of the lines connecting
    the center point to each of its neighbors.

    The finiteness assumption is employed in part in the
    assumption that the polygon is bounded by the finite
    line segments from point 1 to 2, 2 to 3, ...,
    M-1 to M, and M to 1, where M is the number of neighbors.

    It is assumed that this routine is being called by a
    process which has computed the Voronoi diagram of a large
    set of nodes, so the arrays X and Y are dimensioned by
    NODE_NUM, which may be much greater than the number of neighbor
    nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 February 2005
  Author:
    John Burkardt
  Reference:
    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
    Second Edition,
    Wiley, 2000, page 490.
  Parameters:
    Input, int NODE, the index of the node whose Voronoi
    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
    the given node.
    Input, int NEIGHBOR_INDEX[NEIGHBOR_NUM], the indices
    of the neighbor nodes.  These indices are used to access the
    X and Y arrays.  The neighbor nodes should be listed in the
 (counter-clockwise) order in which they occur as one circles
    the center node.
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Output, double *VORONOI_POLYGON_CENTROID_2D, a pointer to a 2D array
    containing the coordinates of the centroid of the Voronoi polygon
    of node NODE.
*/
{
	const _3dtpipit * const s_data = data; 
	
	const register dim_typ node = s_data->a0;
	const register dim_typ neighbor_num = s_data->a1;
	const register dim_typ node_num = s_data->a2;
	int * neighbor_index = s_data->a3;
	ityp * node_xy = s_data->a4;
	
    double a;
    ityp area;
    ityp b;
    ityp c;
    ityp *centroid;
    dim_typ i;
    dim_typ ip1;
    ityp ui;
    ityp uip1;
    ityp vi;
    ityp vip1;
    ityp xc;
    ityp xi;
    ityp xip1;
    ityp yc;
    ityp yi;
    ityp yip1;

    centroid = ( ityp * ) malloc ( sizeof ( ityp ) << 1 );

    centroid[0] = centroid[1] = 0.00;

    if (node_num <= node )
        return NULL;

    xc = node_xy[0+(node<<1)];
    yc = node_xy[1+(node<<1)];

    // possible error
    // i = neighbor_num - 1;
    i = neighbor_index[i];

    xi = node_xy[0+(i<<1)];
    yi = node_xy[1+(i<<1)];

    ip1 = 0;
    ip1 = neighbor_index[ip1];

    xip1 = node_xy[0+(ip1<<1)];
    yip1 = node_xy[1+(ip1<<1)];
    a = ( xi   * xi   + yi   * yi   - xc * xc - yc * yc );
    b = ( xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc );
    c = 2.00 * ( ( xi - xc ) * ( yip1 - yc ) - ( xip1 - xc ) * ( yi - yc ) );
    uip1 = ( a * ( yip1 - yc ) - b * ( yi - yc )  ) / c;
    vip1 = ( a * ( xip1 - xc ) - b * ( xi - xc )  ) / c;

    for ( i = 0; i < neighbor_num; ++i )
    {
        xi = xip1;
        yi = yip1;
        ui = uip1;
        vi = vip1;

        ip1 = i4_wrap ( i+1, 0, neighbor_num-1 );
        ip1 = neighbor_index[ip1];

        xip1 = node_xy[0+(ip1<<1)];
        yip1 = node_xy[1+(ip1<<1)];
        a = ( xi   * xi   + yi   * yi   - xc * xc - yc * yc );
        b = ( xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc );
        c = 2.00 * ( ( xi - xc ) * ( yip1 - yc ) - ( xip1 - xc ) * ( yi - yc ) );
        uip1 = ( a * ( yip1 - yc ) - b * ( yi - yc )  ) / c;
        vip1 = ( a * ( xip1 - xc ) - b * ( xi - xc )  ) / c;

        centroid[0] += ( vi - vip1 )* ( ( uip1 + ui ) * ( uip1 + ui ) - uip1 * ui );
        centroid[1] += ( ui - uip1 )* ( ( vip1 + vi ) * ( vip1 + vi ) - vip1 * vi );
    }

    area = voronoi_polygon_area ( node, neighbor_num, neighbor_index,node_num, node_xy );

    centroid[0] /= ( 6.00 * area );
    centroid[1] /= ( 6.00 * area );

    return centroid;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _voronoi_polygon_vertices ( void * data)
/******************************************************************************/
/*
  Purpose:
    VORONOI_POLYGON_VERTICES_2D computes the vertices of a Voronoi polygon.
  Formula:
    This routine is only appropriate for Voronoi polygons that are finite.
    The routine is given the indices of the nodes that are neighbors of a
    given "center" node.  A node is a neighbor of the center node if the
    Voronoi polygons of the two nodes share an edge.  The triangles of the
    Delaunay triangulation are formed from successive pairs of these neighbor
    nodes along with the center node.
    Given only the neighbor node information, it is possible to determine
    the location of the vertices of the polygonal Voronoi region by computing
    the circumcenters of the Delaunay triangles.
    It is assumed that this routine is being called by a process which has
    computed the Voronoi diagram of a large set of nodes, so the arrays X and
    Y are dimensioned by NODE_NUM, which may be much greater than the number
    of neighbor nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 February 2005
  Author:
    John Burkardt
  Reference:
    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
    Second Edition,
    Wiley, 2000.
  Parameters:
    Input, int NODE, the index of the node whose Voronoi
    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
    the given node.
    Input, int NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
    of the neighbor nodes.  These indices are used to access the
    X and Y arrays.  The neighbor nodes should be listed in the
 (counter-clockwise) order in which they occur as one circles
    the center node.
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Output, double V[2*NEIGHBOR_NUM], the vertices of the Voronoi polygon
    around node NODE.
*/
{
	const dtpit2dtpipit * const s_data = data; 
	
	const register dim_typ node = s_data->a0;
	ityp * node_xy = s_data->a1;
	const register dim_typ neighbor_num = s_data->a2;
	const register dim_typ node_num = s_data->a3;
	int * neighbor_index = s_data->a4;
	ityp * v = s_data->a5;
	
    # define DIM_NUM 2

    ityp *center;
    dim_typ i;
    dim_typ ip1;
    ityp t[DIM_NUM*3];

    if (node_num <= node )
        return NULL;

    t[0] = node_xy[0+(node<<1)];
    t[1] = node_xy[1+(node<<1)];

    ip1 = neighbor_index[0];
    t[4] = node_xy[0+(ip1<<1)];
    t[5] = node_xy[1+(ip1<<1)];

    for ( i = 0; i < neighbor_num; ++i )
    {
        t[2] = t[4];
        t[3] = t[5];

        ip1 = i4_wrap ( i+1, 0, neighbor_num-1 );
        ip1 = neighbor_index[ip1];

        t[0+2*2] = node_xy[0+ip1*2];
        t[1+2*2] = node_xy[1+ip1*2];

        center = triangle_circumcenter_2d ( t );

        v[0+(i<<1)] = center[0];
        v[1+(i<<1)] = center[1];

        free ( center );
    }

    return NULL;
    # undef DIM_NUM
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangulation_neighbor_triangles ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines triangle neighbors.
  Discussion:
    A triangulation of a set of nodes can be completely described by
    the coordinates of the nodes, and the list of nodes that make up
    each triangle.  However, in some cases, it is necessary to know
    triangle adjacency information, that is, which triangle, if any,
    is adjacent to a given triangle on a particular side.
    This routine creates a data structure recording this information.
    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
    data items.
    This routine was modified to work with columns rather than rows.
  Example:
    The input information from TRIANGLE_NODE:
    Triangle   Nodes
    --------   ---------------
     1         3      4      1
     2         3      1      2
     3         3      2      8
     4         2      1      5
     5         8      2     13
     6         8     13      9
     7         3      8      9
     8        13      2      5
     9         9     13      7
    10         7     13      5
    11         6      7      5
    12         9      7      6
    13        10      9      6
    14         6      5     12
    15        11      6     12
    16        10      6     11
    The output information in TRIANGLE_NEIGHBOR:
    Triangle  Neighboring Triangles
    --------  ---------------------
     1        -1     -1      2
     2         1      4      3
     3         2      5      7
     4         2     -1      8
     5         3      8      6
     6         5      9      7
     7         3      6     -1
     8         5      4     10
     9         6     10     12
    10         9      8     11
    11        12     10     14
    12         9     11     13
    13        -1     12     16
    14        11     -1     15
    15        16     14     -1
    16        13     15     -1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int TRIANGLE_ORDER, the order of the triangles.
    Input, int TRIANGLE_NUM, the number of triangles.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], the nodes that
    make up each triangle.
    Output, int TRIANGLE_ORDER3_NEIGHBOR_TRIANGLES[3*TRIANGLE_NUM],
    the three triangles
    that are direct neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I)
    is the index of the triangle which touches side 1, defined by nodes 2
    and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative if there is no
    neighbor on that side.  In this case, that side of the triangle lies
    on the boundary of the triangulation.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ triangle_order = s_data->a0;
	const register dim_typ triangle_num = s_data->a1;
	int * triangle_node = s_data->a2;
	
    int *col;
    dim_typ i;
    dim_typ icol;
    dim_typ j;
    dim_typ k;
    int side1;
    int side2;
    dim_typ tri;
    dim_typ tri1;
    dim_typ tri2;
    int *triangle_neighbor;

    triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
    col = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) << 2 );
    /*
    Step 1.
    From the list of nodes for triangle T, of the form: (I,J,K)
    construct the three neighbor relations:

 (I,J,3,T) or (J,I,3,T),
 (J,K,1,T) or (K,J,1,T),
 (K,I,2,T) or (I,K,2,T)

    where we choose (I,J,3,T) if I < J, or else (J,I,3,T)
    */
    for ( tri = 0; tri < triangle_num; ++tri )
    {
        i = triangle_node[0+tri*triangle_order];
        j = triangle_node[1+tri*triangle_order];
        k = triangle_node[2+tri*triangle_order];

        if ( i < j )
        {
            col[0+((3*tri+0)<<2)] = i;
            col[1+((3*tri+0)<<2)] = j;
            col[2+((3*tri+0)<<2)] = 3;
            col[3+((3*tri+0)<<2)] = tri + 1;
        }
        else
        {
            col[0+((3*tri+0)<<2)] = j;
            col[1+((3*tri+0)<<2)] = i;
            col[2+((3*tri+0)<<2)] = 3;
            col[3+((3*tri+0)<<2)] = tri + 1;
        }

        if ( j < k )
        {
            col[0+((3*tri+1)<<2)] = j;
            col[1+((3*tri+1)<<2)] = k;
            col[2+((3*tri+1)<<2)] = 1;
            col[3+((3*tri+1)<<2)] = tri + 1;
        }
        else
        {
            col[0+((3*tri+1)<<2)] = k;
            col[1+((3*tri+1)<<2)] = j;
            col[2+((3*tri+1)<<2)] = 1;
            col[3+((3*tri+1)<<2)] = tri + 1;
        }

        if ( k < i )
        {
            col[0+((3*tri+2)<<2)] = k;
            col[1+((3*tri+2)<<2)] = i;
            col[2+((3*tri+2)<<2)] = 2;
            col[3+((3*tri+2)<<2)] = tri + 1;
        }
        else
        {
            col[0+((3*tri+2)<<2)] = i;
            col[1+((3*tri+2)<<2)] = k;
            col[2+((3*tri+2)<<2)] = 2;
            col[3+((3*tri+2)<<2)] = tri + 1;
        }
    }
    /*
    Step 2. Perform an ascending dictionary sort on the neighbor relations.
    We only intend to sort on rows 1 and 2; the routine we call here
    sorts on rows 1 through 4 but that won't hurt us.

    What we need is to find cases where two triangles share an edge.
    Say they share an edge defined by the nodes I and J.  Then there are
    two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
    we make sure that these two columns occur consecutively.  That will
    make it easy to notice that the triangles are neighbors.
    */
    i4col_sort_a ( 4, 3*triangle_num, col );
    /*
    Step 3. Neighboring triangles show up as consecutive columns with
    identical first two entries.  Whenever you spot this happening,
    make the appropriate entries in TRIANGLE_NEIGHBOR.
    */
    for ( j = 0; j < triangle_num; ++j )
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; i++ )
            triangle_neighbor[i+j*3] = -1;

    icol = 1;

    for ( ; ; )
    {
        if ( 3 * triangle_num <= icol )
            break;

        if ( col[0+((icol-1)<<2)] != col[0+(icol<<2)] ||col[1+((icol-1)<<2)] != col[1+(icol<<2)] )
        {
            ++ icol;
            continue;
        }

        side1 = col[2+((icol-1)<<2)];
        tri1 =  col[3+((icol-1)<<2)];
        side2 = col[2+ (icol   <<2)];
        tri2 =  col[3+ (icol   <<2)];

        triangle_neighbor[side1-1+(tri1-1)*3] = tri2;
        triangle_neighbor[side2-1+(tri2-1)*3] = tri1;

        icol += 2;
    }

    free ( col );

    return triangle_neighbor;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _triangulation_search_delaunay ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGULATION_SEARCH_DELAUNAY searches a triangulation for a point.
  Discussion:
    The algorithm "walks" from one triangle to its neighboring triangle,
    and so on, until a triangle is found containing point P, or P is found
    to be outside the convex hull.
    The algorithm computes the barycentric coordinates of the point with
    respect to the current triangle.  If all three quantities are positive,
    the point is contained in the triangle.  If the I-th coordinate is
    negative, then (X,Y) lies on the far side of edge I, which is opposite
    from vertex I.  This gives a hint as to where to search next.
    For a Delaunay triangulation, the search is guaranteed to terminate.
    For other triangulations, a cycle may occur.
    Note the surprising fact that, even for a Delaunay triangulation of
    a set of nodes, the nearest point to (X,Y) need not be one of the
    vertices of the triangle containing (X,Y)
    The code can be called for triangulations of any order, but only
    the first three nodes in each triangle are considered.  Thus, if
    higher order triangles are used, and the extra nodes are intended
    to give the triangle a polygonal shape, these will have no effect,
    and the results obtained here might be misleading.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2009
  Author:
    Original FORTRAN77 version by Barry Joe.
    C version by John Burkardt.
  Reference:
    Barry Joe,
    GEOMPACK - a software package for the generation of meshes
    using geometric algorithms,
    Advances in Engineering Software,
    Volume 13, pages 325-331, 1991.
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
    Input, int TRIANGLE_ORDER, the order of the triangles.
    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
    the nodes of each triangle.
    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list.
    Input, double P[2], the coordinates of a point.
    Output, int *TRIANGLE_INDEX, the index of the triangle where the search ended.
    If a cycle occurred, then TRIANGLE_INDEX = -1.
    Output, double *ALPHA, *BETA, *GAMMA, the barycentric
    coordinates of the point relative to triangle TRIANGLE_INDEX.
    Output, int *EDGE, indicates the position of the point (X,Y) in
    triangle TRIANGLE:
    0, the interior or boundary of the triangle;
    -1, outside the convex hull of the triangulation, past edge 1;
    -2, outside the convex hull of the triangulation, past edge 2;
    -3, outside the convex hull of the triangulation, past edge 3.
    Output, int *STEP_NUM, the number of steps.
*/
{
	const dtpit2dt2pipitpi3pit2pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	ityp * node_xy = s_data->a1;
	const register dim_typ triangle_order = s_data->a2;
	const register dim_typ triangle_num = s_data->a3;
	int * triangle_node = s_data->a4;
	int * triangle_neighbor = s_data->a5;
	ityp * p = s_data->a6;
	int * triangle_index = s_data->a7;
	ityp * alpha = s_data->a8;
	ityp * beta = s_data->a9;
	ityp * gamma = s_data->a10;
	int * edge = s_data->a11;
	int * step_num = s_data->a12;
	
    int a;
    int b;
    int c;
    ityp det;
    ityp dxp;
    ityp dxa;
    ityp dxb;
    ityp dyp;
    ityp dya;
    ityp dyb;
    static int triangle_index_save = -1;

    *step_num = - 1;
    *edge = 0;

    *triangle_index = triangle_index_save < 0 || triangle_num <= triangle_index_save ? ( triangle_num + 1 ) / 2 - 1 : triangle_index_save;

    for ( ; ; )
    {
        ++ *step_num;

        if ( triangle_num < *step_num )
            return NULL;
        /*
        Get the vertices of triangle TRIANGLE.
        */
        a = triangle_node[0+(*triangle_index)*triangle_order];
        b = triangle_node[1+(*triangle_index)*triangle_order];
        c = triangle_node[2+(*triangle_index)*triangle_order];
        /*
        Using vertex C as a base, compute the distances to vertices A and B,
        and the point (X,Y).
        */
        dxa = node_xy[0+a*2] - node_xy[0+c*2];
        dya = node_xy[1+a*2] - node_xy[1+c*2];

        dxb = node_xy[0+b*2] - node_xy[0+c*2];
        dyb = node_xy[1+b*2] - node_xy[1+c*2];

        dxp = p[0]           - node_xy[0+c*2];
        dyp = p[1]           - node_xy[1+c*2];

        det = dxa * dyb - dya * dxb;
        /*
        Compute the barycentric coordinates of the point (X,Y) with respect
        to this triangle.
        */
        *alpha = ( dxp * dyb - dyp * dxb ) / det;
        *beta = ( dxa * dyp - dya * dxp ) / det;
        *gamma = 1.00 - *alpha - *beta;
        /*
        If the barycentric coordinates are all positive, then the point
        is inside the triangle and we're done.
        */
        if ( 0.00 <= *alpha &&0.00 <= *gamma )
            break;
        /*
        At least one barycentric coordinate is negative.

        If there is a negative barycentric coordinate for which there exists
        an opposing triangle neighbor closer to the point, move to that triangle.

     (Two coordinates could be negative, in which case we could go for the
        most negative one, or the most negative one normalized by the actual
        distance it represents).
        */
        if ( *alpha < 0.00 && 0 <= triangle_neighbor[1+(*triangle_index)*3] )
        {
            *triangle_index = triangle_neighbor[1+(*triangle_index)*3];
            continue;
        }
        else if ( *beta < 0.0 && 0 <= triangle_neighbor[2+(*triangle_index)*3] )
        {
            *triangle_index = triangle_neighbor[2+(*triangle_index)*3];
            continue;
        }
        else if ( *gamma < 0.00 && 0 <= triangle_neighbor[0+(*triangle_index)*3] )
        {
            *triangle_index = triangle_neighbor[0+(*triangle_index)*3];
            continue;
        }
        /*
        All negative barycentric coordinates correspond to vertices opposite
        sides on the convex hull.

        Note the edge and exit.
        */
        if ( *alpha < 0.00 )
        {
            *edge = -2;
            break;
        }
        else if ( *beta < 0.00 )
        {
            *edge = -3;
            break;
        }
        else if ( *gamma < 0.00 )
        {
            *edge = -1;
            break;
        }
        else
            return NULL;
    }
    triangle_index_save = *triangle_index;

    return NULL;
}

#endif
