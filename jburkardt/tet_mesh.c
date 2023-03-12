#ifndef __DISABLEDEEP_TETMESH

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tet_mesh_neighbor_tets ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_NEIGHBOR_TETS determines tetrahedron neighbors.
  Discussion:
    A tet mesh of a set of nodes can be completely described by
    the coordinates of the nodes, and the list of nodes that make up
    each tetrahedron.  In the most common case, four nodes are used.
    There is also a 10 node case, where nodes are also placed on
    the midsides of the tetrahedral edges.
    This routine can handle 4 or 10-node tetrahedral meshes.  The
    10-node case is handled simply by ignoring the six midside nodes,
    which are presumed to be listed after the vertices.
    The tetrahedron adjacency information records which tetrahedron
    is adjacent to a given tetrahedron on a particular face.
    This routine creates a data structure recording this information.
    The primary amount of work occurs in sorting a list of 4 * TETRA_NUM
    data items.
    The neighbor tetrahedrons are indexed by the face they share with
    the tetrahedron.
    Each face of the tetrahedron is indexed by the node which is NOT
    part of the face.  That is:
    * Neighbor 1 shares face 1 defined by nodes 2, 3, 4.
    * Neighbor 2 shares face 2 defined by nodes 1, 3, 4;
    * Neighbor 3 shares face 3 defined by nodes 1, 2, 4;
    * Neighbor 4 shares face 4 defined by nodes 1, 2, 3.
    For instance, if the (transposed) TETRA_NODE array was:
    Row       1      2      3      4
    Col
      1       4      3      5      1
      2       4      2      5      1
      3       4      7      3      5
      4       4      7      8      5
      5       4      6      2      5
      6       4      6      8      5
    then the (transposed) TETRA_NEIGHBOR array should be:
    Row       1      2      3      4
    Col
      1      -1      2     -1      3
      2      -1      1     -1      5
      3      -1      1      4     -1
      4      -1      6      3     -1
      5      -1      2      6     -1
      6      -1      4      5     -1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 December 2006
  Author:
    John Burkardt
  Parameters:
    Input, int TETRA_ORDER, the order of the tetrahedrons.
    Input, int TETRA_NUM, the number of tetrahedrons.
    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
    Output, int TET_MESH_NEIGHBORS[4*TETRA_NUM], the four tetrahedrons that
    are direct neighbors of a given tetrahedron.  If there is no neighbor
    sharing a given face, the index is set to -1.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ tetra_order = s_data->a0;
	const register dim_typ tetra_num = s_data->a1;
	int * tetra_node = s_data->a2;
	
    dim_typ a;
    dim_typ b;
    dim_typ c;
    dim_typ face;
    dim_typ face1;
    dim_typ face2;
    int *faces;
    int i, j, k;
    dim_typ l;
    dim_typ tetra;
    int *tetra_neighbor;
    dim_typ tetra1;
    dim_typ tetra2;

    faces = ( int * ) malloc ( 5 * tetra_num * sizeof ( int ) << 2 );
    tetra_neighbor = ( int * ) malloc ( tetra_num * sizeof ( int ) << 2);
    /*
    Step 1.
    From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
    construct the four face relations:

 (J,K,L,1,T)
 (I,K,L,2,T)
 (I,J,L,3,T)
 (I,J,K,4,T)

    In order to make matching easier, we reorder each triple of nodes
    into ascending order.
    */
    for ( tetra = 0; tetra < tetra_num; ++tetra )
    {
        i = tetra_node[0+tetra*tetra_order];
        j = tetra_node[1+tetra*tetra_order];
        k = tetra_node[2+tetra*tetra_order];
        l = tetra_node[3+tetra*tetra_order];

        i4i4i4_sort_a ( j, k, l, &a, &b, &c );

        faces[tetra*20] = a;
        faces[1+tetra*20] = b;
        faces[2+tetra*20] = c;
        faces[3+tetra*20] = 0;
        faces[4+tetra*20] = tetra;

        i4i4i4_sort_a ( i, k, l, &a, &b, &c );

        faces[5+tetra*20] = a;
        faces[6+tetra*20] = b;
        faces[7+tetra*20] = c;
        faces[8+tetra*20] = 1;
        faces[9+tetra*20] = tetra;

        i4i4i4_sort_a ( i, j, l, &a, &b, &c );

        faces[10+tetra*20] = a;
        faces[11+tetra*20] = b;
        faces[2+10+tetra*20] = c;
        faces[3+10+tetra*20] = 2;
        faces[4+10+tetra*20] = tetra;

        i4i4i4_sort_a ( i, j, k, &a, &b, &c );

        faces[15+tetra*20] = a;
        faces[16+tetra*20] = b;
        faces[17+tetra*20] = c;
        faces[18+tetra*20] = 3;
        faces[19+tetra*20] = tetra;
    }
    /*
    Step 2. Perform an ascending dictionary sort on the neighbor relations.
    We only intend to sort on rows 1:3; the routine we call here
    sorts on rows 1 through 5 but that won't hurt us.

    What we need is to find cases where two tetrahedrons share a face.
    By sorting the columns of the FACES array, we will put shared faces
    next to each other.
    */
    i4col_sort_a ( 5, tetra_num<<2, faces );
    /*
    Step 3. Neighboring tetrahedrons show up as consecutive columns with
    identical first three entries.  Whenever you spot this happening,
    make the appropriate entries in TETRA_NEIGHBOR.
    */
    for ( j = 0; j < tetra_num; ++j )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
            tetra_neighbor[i+j*4] = -1;

    face = 0;

    for ( ; ; )
    {
        if ( (tetra_num<<2) - 1 <= face )
            break;

        if ( faces[0+face*5] == faces[0+(face+1)*5] &&faces[1+face*5] == faces[1+(face+1)*5] &&faces[2+face*5] == faces[2+(face+1)*5] )
        {
            face1 = faces[3+face*5];
            tetra1 = faces[4+face*5];
            face2 = faces[3+(face+1)*5];
            tetra2 = faces[4+(face+1)*5];
            tetra_neighbor[face1+(tetra1<<2)] = tetra2;
            tetra_neighbor[face2+(tetra2<<2)] = tetra1;
            face += 2;
        }
        else
            ++ face1;
    }

    free ( faces );

    return tetra_neighbor;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tet_mesh_node_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_NODE_ORDER: determines the order of nodes.
  Discussion:
    The order of a node is the number of tetrahedrons that use that node
    as a vertex.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 October 2005
  Author:
    John Burkardt
  Parameters:
    Input, int TETRA_ORDER, the order of the tetrahedrons.
    Input, int TETRA_NUM, the number of tetrahedrons.
    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
    Input, int NODE_NUM, the number of nodes.
    Output, int TET_MESH_NODE_ORDER[NODE_NUM], the order of each node.
*/
{
	const _3dtpi * const s_data = data;
	
	const register dim_typ tetra_order = s_data->a0;
	const register dim_typ tetra_num = s_data->a1;
	const register dim_typ node_num = s_data->a2;
	int * tetra_node = s_data->a3;	
	
    dim_typ i;
    int node;
    int *node_order;
    dim_typ tetra;

    node_order = ( int * ) malloc ( node_num * sizeof ( int ) );
    i4vec_zero ( node_num, node_order );

    for ( tetra = 0; tetra < tetra_num; ++tetra)
        for ( i = 0; i < tetra_order; ++i )
        {
            node = tetra_node[i+tetra*tetra_order];
            if ( node < 0 || node_num <= node )
                return NULL;
            else
                ++ node_order[node];
        }

    return node_order;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tet_mesh_order4_adj_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDEr4_ADJ_COUNT counts the number of nodal adjacencies.
  Discussion:
    Assuming that the tet mesh is to be used in a finite element
    computation, we declare that two distinct nodes are "adjacent" if and
    only if they are both included in some tetrahedron.
    It is the purpose of this routine to determine the number of
    such adjacency relationships.
    The initial count gets only the (I,J) relationships, for which
    node I is strictly less than node J.  This value is doubled
    to account for symmetry.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int TETRA_NUM, the number of tetrahedrons.
    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
    Output, int *ADJ_NUM, the total number of adjacency relationships,
    Output, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
*/
{
	const _2dt3pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ tetra_num = s_data->a1;
	int * tetra_node = s_data->a2;
	int * adj_num = s_data->a3;
	int * adj_row = s_data->a4;
	
    int i, j, k;
    int node;
    int *pair;
    dim_typ pair_num;
    dim_typ pair_unique_num;
    dim_typ tetra;
    /*
    Each order 4 tetrahedron defines 6 adjacency pairs.
    */
    pair = ( int * ) malloc ( 6 * tetra_num * sizeof ( int ) << 1 );

    for ( tetra = 0; tetra < tetra_num; ++tetra )
    {
        pair[0+             tetra *2] = tetra_node[0+tetra*4];
        pair[1+             tetra *2] = tetra_node[1+tetra*4];

        pair[0+(  tetra_num+tetra)*2] = tetra_node[0+tetra*4];
        pair[1+(  tetra_num+tetra)*2] = tetra_node[2+tetra*4];

        pair[0+(2*tetra_num+tetra)*2] = tetra_node[0+tetra*4];
        pair[1+(2*tetra_num+tetra)*2] = tetra_node[3+tetra*4];

        pair[0+(3*tetra_num+tetra)*2] = tetra_node[1+tetra*4];
        pair[1+(3*tetra_num+tetra)*2] = tetra_node[2+tetra*4];

        pair[0+(4*tetra_num+tetra)*2] = tetra_node[1+tetra*4];
        pair[1+(4*tetra_num+tetra)*2] = tetra_node[3+tetra*4];

        pair[0+(5*tetra_num+tetra)*2] = tetra_node[2+tetra*4];
        pair[1+(5*tetra_num+tetra)*2] = tetra_node[3+tetra*4];
    }
    pair_num = 6 * tetra_num;
    /*
    Force the nodes of each pair to be listed in ascending order.
    */
    i4col_sort2_a ( 2, pair_num, pair );
    /*
    Rearrange the columns in ascending order.
    */
    i4col_sort_a ( 2, pair_num, pair );
    /*
    Get the number of unique columns.
    */
    pair_unique_num = i4col_sorted_unique_count ( 2, pair_num, pair );
    /*
    The number of adjacencies is TWICE this value, plus the number of nodes.
    */
    *adj_num = pair_unique_num<<1;
    /*
    Now set up the ADJ_ROW counts.
    */
    for ( node = 0; node < node_num; ++node)
        adj_row[node] = 0;

    for ( k = 0; k < pair_num; ++k )
    {
        if (0<k && pair[0+((k-1)<<1)] == pair[0+(k<<1)] &&pair[1+((k-1)<<1)] == pair[1+(k<<1)] )
            continue;
        i = pair[0+k*2];
        j = pair[1+k*2];

        ++ adj_row[i-1];
        ++ adj_row[j-1];
    }
    /*
    We used ADJ_ROW to count the number of entries in each row.
    Convert it to pointers into the ADJ array.
    */
    for ( node = node_num-1; 0 <= node; --node )
        adj_row[node] = adj_row[node+1];

    adj_row[0] = 1;
    for ( node = 1; node <= node_num; ++node )
        adj_row[node] += adj_row[i];

    free ( pair );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tet_mesh_order4_adj_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDEr4_ADJ_SET sets the nodal adjacency matrix.
  Discussion:
    A compressed format is used for the nodal adjacency matrix.
    It is assumed that we know ADJ_NUM, the number of adjacency entries
    and the ADJ_ROW array, which keeps track of the list of slots
    in ADJ where we can store adjacency information for each row.
    We essentially repeat the work of TET_MESH_ORDEr8_ADJ_COUNT, but
    now we have a place to store the adjacency information.
    A copy of the ADJ_ROW array is useful, as we can use it to keep track
    of the next available entry in ADJ for adjacencies associated with
    a given row.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int TETRA_NUM, the number of tetrahedrons.
    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
    Input, int ADJ_NUM, the total number of adjacency relationships,
    Input, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
    Output, int TET_MESH_ORDEr8_ADJ_SET[ADJ_NUM],
    the adjacency information.
*/
{
	const _3dt2pi * const s_data = data;
	
	const register dim_typ node_num = s_data->a0;
	const register dim_typ element_num = s_data->a1;
	const register dim_typ adj_num = s_data->a2;
	int * element_node = s_data->a3;
	int * adj_row = s_data->a4;
	
    int *adj;
    int *adj_row_copy;
    int i, j, k;
    int node;
    int *pair;
    dim_typ pair_num;
    dim_typ tetra;
    /*
    Each order 4 tetrahedron defines 6 adjacency pairs.
    */
    pair = ( int * ) malloc ( 6 * element_num * sizeof ( int ) << 1);

    for ( tetra = 0; tetra < element_num; ++tetra)
    {
        pair[0+          (tetra <<1)] = element_node[(tetra<<2)];
        pair[1+          (tetra <<1)] = element_node[1+(tetra<<2)];

        pair[0+((  element_num+tetra)<<1)] = element_node[(tetra<<2)];
        pair[1+((  element_num+tetra)<<1)] = element_node[2+(tetra<<2)];

        pair[0+(((element_num<<1)+tetra)<<1)] = element_node[0+(tetra<<2)];
        pair[1+(((element_num<<1)+tetra)<<1)] = element_node[3+(tetra<<2)];

        pair[0+((3*element_num+tetra)<<1)] = element_node[1+(tetra<<2)];
        pair[1+((3*element_num+tetra)<<1)] = element_node[2+(tetra<<2)];

        pair[0+(((element_num<<2)+tetra)<<1)] = element_node[1+(tetra<<2)];
        pair[1+(((element_num<<2)+tetra)<<1)] = element_node[3+(tetra<<2)];

        pair[0+((5*element_num+tetra)<<1)] = element_node[2+(tetra<<2)];
        pair[1+((5*element_num+tetra)<<1)] = element_node[3+(tetra<<2)];
    }
    pair_num = 6 * element_num;
    /*
    Force the nodes of each pair to be listed in ascending order.
    */
    i4col_sort2_a ( 2, pair_num, pair );
    /*
    Rearrange the columns in ascending order.
    */
    i4col_sort_a ( 2, pair_num, pair );
    /*
    Mark all entries of ADJ so we will know later if we missed one.
    */
    adj = ( int * ) malloc ( adj_num * sizeof ( int ) );

    for ( i = 0; i < adj_num; ++i )
        adj[i] = -1;
    /*
    Copy the ADJ_ROW array and use it to keep track of the next
    free entry for each row.
    */
    adj_row_copy = ( int * ) malloc ( node_num * sizeof ( int ) );

    for ( node = 0; node < node_num; ++node )
        adj_row_copy[node] = adj_row[node];
    /*
    Now set up the ADJ_ROW counts.
    */
    for ( k = 0; k < pair_num; ++k)
    {
        if (0<k && pair[0+((k-1)<<1)] == pair[0+(k<<1)] &&pair[1+((k-1)<<1)] == pair[1+(k<<1)] )
            continue;
        i = pair[0+(k<<1)];
        j = pair[1+(k<<1)];

        adj[adj_row_copy[i]] = j;
        ++ adj_row_copy[i];
        adj[adj_row_copy[j]] = i;
        ++ adj_row_copy[j];
    }
    free ( adj_row_copy );
    free ( pair );

    return adj;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tet_mesh_order4_boundary_face_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDEr4_BOUNDARY_FACE_COUNT counts the number of boundary faces.
  Discussion:
    This routine is given a tet mesh, an abstract list of
    quadruples of nodes.  It is assumed that the nodes forming each
    face of each tetrahedron are listed in a counterclockwise order,
    although the routine should work if the nodes are consistently
    listed in a clockwise order as well.
    It is assumed that each face of the tet mesh is either
    * an INTERIOR face, which is listed twice, once with positive
      orientation and once with negative orientation, or;
    * a BOUNDARY face, which will occur only once.
    This routine should work even if the region has holes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int TETRA_NUM, the number of tetrahedrons.
    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
    Output, int TET_MESH_ORDEr8_BOUNDARY_FACE_COUNT, the number of
    boundary faces.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ tetra_num = s_data->a0;
	int * tetra_node = s_data->a1;
	
    dim_typ boundary_face_num;
    int *face;
    dim_typ face_num;
    dim_typ interior_face_num;
    dim_typ m;
    dim_typ tet;
    dim_typ unique_face_num;

    face = ( int * ) malloc ( 3 * tetra_num * sizeof ( int ) << 2 );

    m = 3;
    face_num = tetra_num<<2;
    /*
    Set up the face array:
 (Omit node 1)
 (Omit node 2)
 (Omit node 3)
 (Omit node 4)
    */
    for ( tet = 0; tet < tetra_num; ++tet )
    {
        face[0+(            tet)*3] = tetra_node[1+(tet<<2)];
        face[1+(            tet)*3] = tetra_node[2+(tet<<2)];
        face[2+(            tet)*3] = tetra_node[3+(tet<<2)];

        face[0+(  tetra_num+tet)*3] = tetra_node[0+(tet<<2)];
        face[1+(  tetra_num+tet)*3] = tetra_node[2+(tet<<2)];
        face[2+(  tetra_num+tet)*3] = tetra_node[3+(tet<<2)];

        face[0+((tetra_num<<1)+tet)*3] = tetra_node[0+(tet<<2)];
        face[1+((tetra_num<<1)+tet)*3] = tetra_node[1+(tet<<2)];
        face[2+((tetra_num<<1)+tet)*3] = tetra_node[3+(tet<<2)];

        face[0+(3*tetra_num+tet)*3] = tetra_node[0+(tet<<2)];
        face[1+(3*tetra_num+tet)*3] = tetra_node[1+(tet<<2)];
        face[2+(3*tetra_num+tet)*3] = tetra_node[2+(tet<<2)];
    }
    /*
    Force the nodes of each face to be listed in ascending order.
    */
    i4col_sort2_a ( m, face_num, face );
    /*
    Ascending sort the columns.
    */
    i4col_sort_a ( m, face_num, face );
    /*
    Get the number of unique columns.
    */
    /*
    Determine the number of interior and boundary faces.
    */
    boundary_face_num = (tetra_num<<1) - ((tetra_num<<2) - i4col_sorted_unique_count ( m, face_num, face )<<1);

    free ( face );
    
    result = boundary_face_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tet_mesh_order4_edge_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDEr4_EDGE_COUNT counts the number of edges.
  Discussion:
    This routine is given a tet mesh, an abstract list of
    quadruples of nodes.  Each tetrahedron defines 6 edges; however,
    assuming that tetrahedrons are touching each other, most edges
    will be used more than once.  This routine determines the actual
    number of "geometric" edges associated with the tet mesh.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int TETRA_NUM, the number of tetrahedrons.
    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
    Output, int TET_MESH_ORDEr8_EDGE_COUNT, the number of edges.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ tetra_num = s_data->a0;
	int * tetra_node = s_data->a1;
	
    int *edge;
    dim_typ edge_num;
    dim_typ edge_num_raw;
    dim_typ m;
    dim_typ tet;

    edge = ( int * ) malloc ( 6 * tetra_num * sizeof ( int ) << 1 );

    m = 3;
    edge_num_raw = 6 * tetra_num;
    /*
    Set up the raw edge array:
    */
    for ( tet = 0; tet < tetra_num; tet++ )
    {
        edge[0+          (tet <<1)] = tetra_node[0+(tet<<2)];
        edge[1+          (tet <<1)] = tetra_node[1+(tet<<2)];

        edge[0+((  tetra_num+tet)<<1)] = tetra_node[0+(tet<<2)];
        edge[1+((  tetra_num+tet)<<1)] = tetra_node[2+(tet<<2)];

        edge[0+(((tetra_num<<1)+tet)<<1)] = tetra_node[0+(tet<<2)];
        edge[1+(((tetra_num<<1)+tet)<<1)] = tetra_node[3+(tet<<2)];

        edge[0+((3*tetra_num+tet)<<1)] = tetra_node[1+(tet<<2)];
        edge[1+((3*tetra_num+tet)<<1)] = tetra_node[2+(tet<<2)];

        edge[0+(((tetra_num<<2)+tet)<<1)] = tetra_node[1+(tet<<2)];
        edge[1+(((tetra_num<<2)+tet)<<1)] = tetra_node[3+(tet<<2)];

        edge[0+((5*tetra_num+tet)<<1)] = tetra_node[2+(tet<<2)];
        edge[1+((5*tetra_num+tet)<<1)] = tetra_node[3+(tet<<2)];
    }
    /*
    Force the nodes of each face to be listed in ascending order.
    */
    i4col_sort2_a ( m, edge_num_raw, edge );
    /*
    Ascending sort the columns.
    */
    i4col_sort_a ( m, edge_num_raw, edge );
    /*
    Get the number of unique columns.
    */
    edge_num = i4col_sorted_unique_count ( m, edge_num_raw, edge );
    free ( edge );
    
    result = edge_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tet_mesh_order4_refine_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDEr4_REFINE_COMPUTE computes a refined order 4 tet mesh
  Discussion:
    A refined 4-node tet mesh can be derived from a given
    4-node tet mesh by interpolating nodes at the midpoint of
    every edge of the mesh.
    The mesh is described indirectly, as the sum of individual
    tetrahedrons.  A single physical edge may be a logical edge of
    any number of tetrahedrons.  It is important, however, that a
    new node be created exactly once for each edge, assigned an index,
    and associated with every tetrahedron that shares this edge.
    This routine handles that problem.
    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
    data items, one item for every edge of every tetrahedron.  Each
    data item records, for a given tetrahedron edge, the global indices
    of the two endpoints, the local indices of the two endpoints,
    and the index of the tetrahedron.
    Through careful sorting, it is possible to arrange this data in
    a way that allows the proper generation of the interpolated nodes.
    Let us add the new nodes and temporarily assign them local indices
    5 through X, based on the following ordering:
      1, 2, 3, 4, (1+2), (1+3), (1+4), (2+3), (2+4), (3+4).
    Then let us assign these nodes to eight subtetrahedrons as follows:
      1, 5, 6, 7
      2, 5, 8, 9
      3, 6, 8, 9
      4, 7, 9, X
      5, 6, 7, 9
      5, 6, 8, 9
      6, 7, 9, X
      6, 8, 9, X
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 January 2007
  Author:
    John Burkardt
  Reference:
    Anwei Liu, Barry Joe,
    Quality Local Refinement of Tetrahedral Meshes Based
    on 8-Subtetrahedron Subdivision,
    Mathematics of Computation,
    Volume 65, Number 215, July 1996, pages 1183-1200.
  Parameters:
    Input, int NODE_NUM1, the number of nodes in the input mesh.
    Input, int TETRA_NUM1, the number of tetrahedrons in the
    input mesh.
    Input, double NODE_XYZ1[3*NODE_NUM1], the coordinates of
    the nodes that make up the input mesh.
    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
    in the input mesh.
    Input, int NODE_NUM2, the number of nodes in the refined mesh.
    Input, int TETRA_NUM2, the number of tetrahedrons in the
    refined mesh.s
    Input, int EDGE_DATA[5*(6*TETRA_NUM1)], edge data.
    Output, double NODE_XYZ2[3*NODE_NUM2], the coordinates of
    the nodes that make up the refined mesh.
    Output, int TETRA_NODE2[4*TETRA_NUM2], the indices of the nodes
    in the refined mesh.
*/
{
	const _2dtpitpi2dtpipitpi * const s_data = data;
	const register dim_typ node_num1 = s_data->a0;
	const register dim_typ tetra_num1 = s_data->a1;
	ityp * node_xyz1 = s_data->a2;
	int * tetra_node1 = s_data->a3;
	const register dim_typ node_num2 = s_data->a4;
	const register dim_typ tetra_num2 = s_data->a5;
	int * edge_data = s_data->a6;
	ityp * node_xyz2 = s_data->a7;
	int * tetra_node2 = s_data->a8;
	
	
    const register dim_typ dim_num = 3;
    dim_typ edge;
    int i, j;
    int n1;
    int n1_old;
    int n2;
    int n2_old;
    int node;
    const register dim_typ tetra_order = 4;
    dim_typ tetra1;
    dim_typ tetra2;
    dim_typ v;
    dim_typ v1;
    dim_typ v2;
    /*
    Generate the index and coordinates of the new midside nodes,
    and update the tetradehron-node data.
    */
    for ( j = 0; j < node_num1; ++j )
        for ( i = 0; i < dim_num; ++i )
            node_xyz2[i+j*dim_num] = node_xyz1[i+j*dim_num];
    for ( j = 0; j < tetra_num2; ++j )
        for ( i = 0; i < tetra_order; ++i )
            tetra_node2[i+j*tetra_order] = -1;
    /*
    The vertices of the input tetrahedron can be assigned now.
    */
    for ( tetra1 = 0; tetra1 < tetra_num1; ++tetra1 )
    {
        tetra_node2[0+(tetra1*8+0)*tetra_order] = tetra_node1[0+tetra1*tetra_order];
        tetra_node2[0+(tetra1*8+1)*tetra_order] = tetra_node1[1+tetra1*tetra_order];
        tetra_node2[0+(tetra1*8+2)*tetra_order] = tetra_node1[2+tetra1*tetra_order];
        tetra_node2[0+(tetra1*8+3)*tetra_order] = tetra_node1[3+tetra1*tetra_order];
    }
    node = node_num1;

    n1_old = -1;
    n2_old = -1;

    for ( edge = 0; edge < 6 * tetra_num1; ++edge)
    {
        /*
        Read the data defining the edge.
        */
        n1 = edge_data[0+edge*5];
        n2 = edge_data[1+edge*5];
        /*
        If this edge is new, create the coordinates and index.
        */
        if ( n1 != n1_old || n2 != n2_old )
        {
            if ( node_num2 <= node )
                return NULL;

            for ( i = 0; i < dim_num; ++i )
                node_xyz2[i+node*dim_num] =( node_xyz2[i+(n1-1)*dim_num] + node_xyz2[i+(n2-1)*dim_num] ) / 2.00;
            ++ node ;
            n1_old = n1;
            n2_old = n2;
        }
        /*
        Assign the node to the tetrahedron.
        */
        v1 = edge_data[2+edge*5];
        v2 = edge_data[3+edge*5];
        tetra1 = edge_data[4+edge*5];
        /*
        We know the two vertices that bracket this new node.
        This tells us whether it is new node number 5, 6, 7, 8, 9 or 10.
        This tells us which of the new subtetrahedrons it belongs to,
        and what position it occupies.
        */
        if ( v1 == 1 && v2 == 2 )
        {
            tetra_node2[1+((tetra1<<3)+0)*tetra_order] = node;
            tetra_node2[1+((tetra1<<3)+1)*tetra_order] = node;
            tetra_node2[0+((tetra1<<3)+4)*tetra_order] = node;
            tetra_node2[0+((tetra1<<3)+5)*tetra_order] = node;
        }
        else if ( v1 == 1 && v2 == 3 )
        {
            tetra_node2[2+((tetra1<<3)+0)*tetra_order] = node;
            tetra_node2[1+((tetra1<<3)+2)*tetra_order] = node;
            tetra_node2[1+((tetra1<<3)+4)*tetra_order] = node;
            tetra_node2[1+((tetra1<<3)+5)*tetra_order] = node;
            tetra_node2[0+((tetra1<<3)+6)*tetra_order] = node;
            tetra_node2[0+((tetra1<<3)+7)*tetra_order] = node;
        }
        else if ( v1 == 1 && v2 == 4 )
        {
            tetra_node2[3+((tetra1<<3)+0)*tetra_order] = node;
            tetra_node2[1+((tetra1<<3)+3)*tetra_order] = node;
            tetra_node2[2+((tetra1<<3)+4)*tetra_order] = node;
            tetra_node2[1+((tetra1<<3)+6)*tetra_order] = node;
        }
        else if ( v1 == 2 && v2 == 3 )
        {
            tetra_node2[2+((tetra1<<3)+1)*tetra_order] = node;
            tetra_node2[2+((tetra1<<3)+2)*tetra_order] = node;
            tetra_node2[2+((tetra1<<3)+5)*tetra_order] = node;
            tetra_node2[1+((tetra1<<3)+7)*tetra_order] = node;
        }
        else if ( v1 == 2 && v2 == 4 )
        {
            tetra_node2[3+((tetra1<<3)+1)*tetra_order] = node;
            tetra_node2[3+((tetra1<<3)+2)*tetra_order] = node;
            tetra_node2[2+((tetra1<<3)+3)*tetra_order] = node;
            tetra_node2[3+((tetra1<<3)+4)*tetra_order] = node;
            tetra_node2[3+((tetra1<<3)+5)*tetra_order] = node;
            tetra_node2[2+((tetra1<<3)+6)*tetra_order] = node;
            tetra_node2[2+((tetra1<<3)+7)*tetra_order] = node;
        }
        else if ( v1 == 3 && v2 == 4 )
        {
            tetra_node2[3+((tetra1<<3)+3)*tetra_order] = node;
            tetra_node2[3+((tetra1<<3)+6)*tetra_order] = node;
            tetra_node2[3+((tetra1<<3)+7)*tetra_order] = node;
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tet_mesh_order4_refine_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDEr4_REFINE_SIZE sizes a refined order 4 tet mesh.
  Discussion:
    A refined tet mesh can be derived from an existing one by interpolating
    nodes at the midpoint of every edge of the mesh.
    The mesh is described indirectly, as the sum of individual
    tetrahedrons.  A single physical edge may be a logical edge of
    any number of tetrahedrons.  It is important, however, that a
    new node be created exactly once for each edge, assigned an index,
    and associated with every tetrahedron that shares this edge.
    This routine handles that problem.
    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
    data items, one item for every edge of every tetrahedron.  Each
    data item records, for a given tetrahedron edge, the global indices
    of the two endpoints, the local indices of the two endpoints,
    and the index of the tetrahedron.
    Through careful sorting, it is possible to arrange this data in
    a way that allows the proper generation of the interpolated nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM1, the number of nodes in the original mesh.
    Input, int TETRA_NUM1, the number of tetrahedrons in the
    original mesh.
    Input, int TETRA_NODE1[4*TETRA_NUM1], the indices of the nodes
    in the original mesh.
    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
    Output, int *TETRA_NUM2, the number of tetrahedrons in the refined mesh.
    Output, int EDGE_DATA[5*(6*TETRA_NUM1)], edge data.
*/
{
	const _2dt4pi * const s_data = data;
	register dim_typ node_num1 = s_data->a0;
	const register dim_typ tetra_num1 = s_data->a1;
	int * tetra_node1 = s_data->a2;
	int * node_num2 = s_data->a3;
	int * tetra_num2 = s_data->a4;
	int * edge_data = s_data->a5;
	
    dim_typ a;
    dim_typ b;
    dim_typ edge;
    int i, j, k;
    dim_typ l;
    int n1;
    int n1_old;
    int n2;
    int n2_old;
    dim_typ tetra;
    const register dim_typ tetra_order = 4;
    /*
    Step 1.
    From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
    construct the six edge relations:

 (I,J,1,2,T)
 (I,K,1,3,T)
 (I,L,1,4,T)
 (J,K,2,3,T)
 (J,L,2,4,T)
 (K,L,3,4,T)

    In order to make matching easier, we reorder each pair of nodes
    into ascending order.
    */
    for ( tetra = 0; tetra < tetra_num1; ++tetra )
    {
        i = tetra_node1[0+tetra*tetra_order];
        j = tetra_node1[1+tetra*tetra_order];
        k = tetra_node1[2+tetra*tetra_order];
        l = tetra_node1[3+tetra*tetra_order];

        i4i4_sort_a ( i, j, &a, &b );

        edge_data[0+(6*tetra)*5] = a;
        edge_data[1+(6*tetra)*5] = b;
        edge_data[2+(6*tetra)*5] = 1;
        edge_data[3+(6*tetra)*5] = 2;
        edge_data[4+(6*tetra)*5] = tetra;

        i4i4_sort_a ( i, k, &a, &b );

        edge_data[0+(6*tetra+1)*5] = a;
        edge_data[1+(6*tetra+1)*5] = b;
        edge_data[2+(6*tetra+1)*5] = 1;
        edge_data[3+(6*tetra+1)*5] = 3;
        edge_data[4+(6*tetra+1)*5] = tetra;

        i4i4_sort_a ( i, l, &a, &b );

        edge_data[0+(6*tetra+2)*5] = a;
        edge_data[1+(6*tetra+2)*5] = b;
        edge_data[2+(6*tetra+2)*5] = 1;
        edge_data[3+(6*tetra+2)*5] = 4;
        edge_data[4+(6*tetra+2)*5] = tetra;

        i4i4_sort_a ( j, k, &a, &b );

        edge_data[0+(6*tetra+3)*5] = a;
        edge_data[1+(6*tetra+3)*5] = b;
        edge_data[2+(6*tetra+3)*5] = 2;
        edge_data[3+(6*tetra+3)*5] = 3;
        edge_data[4+(6*tetra+3)*5] = tetra;

        i4i4_sort_a ( j, l, &a, &b );

        edge_data[0+(6*tetra+4)*5] = a;
        edge_data[1+(6*tetra+4)*5] = b;
        edge_data[2+(6*tetra+4)*5] = 2;
        edge_data[3+(6*tetra+4)*5] = 4;
        edge_data[4+(6*tetra+4)*5] = tetra;

        i4i4_sort_a ( k, l, &a, &b );

        edge_data[0+(6*tetra+5)*5] = a;
        edge_data[1+(6*tetra+5)*5] = b;
        edge_data[2+(6*tetra+5)*5] = 3;
        edge_data[3+(6*tetra+5)*5] = 4;
        edge_data[4+(6*tetra+5)*5] = tetra;
    }
    /*
    Step 2. Perform an ascending dictionary sort on the neighbor relations.
    We only intend to sort on rows 1:2; the routine we call here
    sorts on the full column but that won't hurt us.

    What we need is to find all cases where tetrahedrons share an edge.
    By sorting the columns of the EDGE_DATA array, we will put shared edges
    next to each other.
    */
    i4col_sort_a ( 5, 6*tetra_num1, edge_data );
    /*
    Step 3. All the tetrahedrons which share an edge show up as consecutive
    columns with identical first two entries.  Figure out how many new
    nodes there are, and allocate space for their coordinates.
    */
    *node_num2 = node_num1;

    n1_old = n2_old=  -1;

    for ( edge = 0; edge < 6 * tetra_num1; ++edge )
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

    *tetra_num2 = tetra_num1 << 3;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tet_mesh_order4_to_order10_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDEr4_TO_ORDER10_COMPUTE computes a quadratic tet mesh from a linear one.
  Discussion:
    A quadratic (10 node) tet mesh can be derived from a linear
 (4 node) tet mesh by interpolating nodes at the midpoint of
    every edge of the mesh.
    The mesh is described indirectly, as the sum of individual
    tetrahedrons.  A single physical edge may be a logical edge of
    any number of tetrahedrons.  It is important, however, that a
    new node be created exactly once for each edge, assigned an index,
    and associated with every tetrahedron that shares this edge.
    This routine handles that problem.
    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
    data items, one item for every edge of every tetrahedron.  Each
    data item records, for a given tetrahedron edge, the global indices
    of the two endpoints, the local indices of the two endpoints,
    and the index of the tetrahedron.
    Through careful sorting, it is possible to arrange this data in
    a way that allows the proper generation of the interpolated nodes.
    The node ordering for the quadratic tetrahedron is somewhat
    arbitrary.  In the current scheme, the vertices are listed
    first, followed by the 6 midside nodes.  Each midside node
    may be identified by the two vertices that bracket it.  Thus,
    the node ordering may be suggested by:
      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int TETRA_NUM, the number of tetrahedrons in the
    linear mesh.
    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
    in the linear mesh.
    Input, int NODE_NUM1, the number of nodes for the linear mesh.
    Input, double NODE_XYZ1[3*NODE_NUM1], the coordinates of
    the nodes that make up the linear mesh.
    Input, int EDGE_DATA[5*(6*TETRA_NUM)], edge data.
    Output, int TETRA_NODE2[10*TETRA_NUM], the indices of the nodes
    in the quadratic mesh.
    Input, int NODE_NUM2, the number of nodes for the quadratic mesh.
    Output, double NODE_XYZ2[3*NODE_NUM2], the coordinates of
    the nodes that make up the quadratic mesh.
*/
{
	const dtpidtpit2pidtpit * const s_data = data;
	const register dim_typ tetra_num = s_data->a0;
	int * tetra_node1 = s_data->a1;
	const register dim_typ node_num1 = s_data->a2;
	ityp * node_xyz1 = s_data->a3;
	int * edge_data = s_data->a4;
	int * tetra_node2 = s_data->a5;
	const register dim_typ node_num2 = s_data->a6;
	ityp * node_xyz2 = s_data->a7;
	
    const register dim_typ dim_num = 3;
    dim_typ edge;
    int i, j;
    int n1;
    int n1_old;
    int n2;
    int n2_old;
    int node;
    dim_typ tetra;
    const register dim_typ tetra_order1 = 4;
    const register dim_typ tetra_order2 = 10;
    dim_typ v;
    dim_typ v1;
    dim_typ v2;
    /*
    Generate the index and coordinates of the new midside nodes,
    and update the tetradehron-node data.
    */
    for ( j = 0; j < node_num1; ++j )
        for ( i = 0; i < dim_num; ++i)
            node_xyz2[i+j*dim_num] = node_xyz1[i+j*dim_num];
    for ( j = 0; j < tetra_num; ++j)
        for ( i = 0; i < tetra_order1; ++i )
            tetra_node2[i+j*tetra_order2] = tetra_node1[i+j*tetra_order1];
    node = node_num1;

    n1_old = n2_old = -1;

    for ( edge = 0; edge < 6 * tetra_num; ++edge )
    {
        /*
        Read the data defining the edge.
        */
        n1 = edge_data[0+edge*5];
        n2 = edge_data[1+edge*5];
        /*
        If this edge is new, create the coordinates and index.
        */
        if ( n1 != n1_old || n2 != n2_old )
        {
            if ( node_num2 <= node )
                return NULL;

            for ( i = 0; i < dim_num; ++i )
                node_xyz2[i+node*dim_num] = ( node_xyz2[i+(n1-1)*dim_num] + node_xyz2[i+(n2-1)*dim_num] ) / 2.00;
            ++ node;
            n1_old = n1;
            n2_old = n2;
        }
        /*
        Assign the node to the tetrahedron.
        */
        v1 = edge_data[2+edge*5];
        v2 = edge_data[3+edge*5];
        /*
        Here is where the local ordering of the nodes is effected:
        */
        if ( v1 == 1 && v2 == 2 )
            v = 5;
        else if ( v1 == 1 && v2 == 3 )
            v = 6;
        else if ( v1 == 1 && v2 == 4 )
            v = 7;
        else if ( v1 == 2 && v2 == 3 )
            v = 8;
        else if ( v1 == 2 && v2 == 4 )
            v = 9;
        else if ( v1 == 3 && v2 == 4 )
            v = 10;

        tetra = edge_data[4+edge*5];
        tetra_node2[v-1+tetra*tetra_order2] = node;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tet_mesh_order4_to_order10_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDEr4_TO_ORDER10_SIZE sizes a quadratic tet mesh from a linear one.
  Discussion:
    A quadratic (10 node) tet mesh can be derived from a linear
 (4 node) tet mesh by interpolating nodes at the midpoint of
    every edge of the mesh.
    The mesh is described indirectly, as the sum of individual
    tetrahedrons.  A single physical edge may be a logical edge of
    any number of tetrahedrons.  It is important, however, that a
    new node be created exactly once for each edge, assigned an index,
    and associated with every tetrahedron that shares this edge.
    This routine handles that problem.
    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
    data items, one item for every edge of every tetrahedron.  Each
    data item records, for a given tetrahedron edge, the global indices
    of the two endpoints, the local indices of the two endpoints,
    and the index of the tetrahedron.
    Through careful sorting, it is possible to arrange this data in
    a way that allows the proper generation of the interpolated nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 December 2006
  Author:
    John Burkardt
  Parameters:
    Input, int TETRA_NUM, the number of tetrahedrons in the
    linear mesh.
    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
    in the linear mesh.
    Input, int NODE_NUM1, the number of nodes for the linear mesh.
    Output, int EDGE_DATA[5*(6*TETRA_NUM)], edge data.
    Output, int *NODE_NUM2, the number of nodes for the quadratic mesh.
*/
{
	const _2dt3pi * const s_data = data;
	
	const register dim_typ tetra_num = s_data->a0;
	const register dim_typ node_num1 = s_data->a1;
	int * tetra_node1 = s_data->a2;
	int * edge_data = s_data->a3;
	int * node_num2 = s_data->a4;
	
    dim_typ a;
    dim_typ b;
    dim_typ edge;
    int i, j, k;
    dim_typ l;
    int n1;
    int n1_old;
    int n2;
    int n2_old;
    dim_typ tetra;
    const register dim_typ tetra_order1 = 4;
    /*
    Step 1.
    From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
    construct the six edge relations:

 (I,J,1,2,T)
 (I,K,1,3,T)
 (I,L,1,4,T)
 (J,K,2,3,T)
 (J,L,2,4,T)
 (K,L,3,4,T)

    In order to make matching easier, we reorder each pair of nodes
    into ascending order.
    */
    for ( tetra = 0; tetra < tetra_num; ++tetra)
    {
        i = tetra_node1[0+tetra*tetra_order1];
        j = tetra_node1[1+tetra*tetra_order1];
        k = tetra_node1[2+tetra*tetra_order1];
        l = tetra_node1[3+tetra*tetra_order1];

        i4i4_sort_a ( i, j, &a, &b );

        edge_data[0+(6*tetra)*5] = a;
        edge_data[1+(6*tetra)*5] = b;
        edge_data[2+(6*tetra)*5] = 1;
        edge_data[3+(6*tetra)*5] = 2;
        edge_data[4+(6*tetra)*5] = tetra;

        i4i4_sort_a ( i, k, &a, &b );

        edge_data[0+(6*tetra+1)*5] = a;
        edge_data[1+(6*tetra+1)*5] = b;
        edge_data[2+(6*tetra+1)*5] = 1;
        edge_data[3+(6*tetra+1)*5] = 3;
        edge_data[4+(6*tetra+1)*5] = tetra;

        i4i4_sort_a ( i, l, &a, &b );

        edge_data[0+(6*tetra+2)*5] = a;
        edge_data[1+(6*tetra+2)*5] = b;
        edge_data[2+(6*tetra+2)*5] = 1;
        edge_data[3+(6*tetra+2)*5] = 4;
        edge_data[4+(6*tetra+2)*5] = tetra;

        i4i4_sort_a ( j, k, &a, &b );

        edge_data[0+(6*tetra+3)*5] = a;
        edge_data[1+(6*tetra+3)*5] = b;
        edge_data[2+(6*tetra+3)*5] = 2;
        edge_data[3+(6*tetra+3)*5] = 3;
        edge_data[4+(6*tetra+3)*5] = tetra;

        i4i4_sort_a ( j, l, &a, &b );

        edge_data[0+(6*tetra+4)*5] = a;
        edge_data[1+(6*tetra+4)*5] = b;
        edge_data[2+(6*tetra+4)*5] = 2;
        edge_data[3+(6*tetra+4)*5] = 4;
        edge_data[4+(6*tetra+4)*5] = tetra;

        i4i4_sort_a ( k, l, &a, &b );

        edge_data[0+(6*tetra+5)*5] = a;
        edge_data[1+(6*tetra+5)*5] = b;
        edge_data[2+(6*tetra+5)*5] = 3;
        edge_data[3+(6*tetra+5)*5] = 4;
        edge_data[4+(6*tetra+5)*5] = tetra;
    }
    /*
    Step 2. Perform an ascending dictionary sort on the neighbor relations.
    We only intend to sort on rows 1:2; the routine we call here
    sorts on the full column but that won't hurt us.

    What we need is to find all cases where tetrahedrons share an edge.
    By sorting the columns of the EDGE_DATA array, we will put shared edges
    next to each other.
    */
    i4col_sort_a ( 5, 6*tetra_num, edge_data );
    /*
    Step 3. All the tetrahedrons which share an edge show up as consecutive
    columns with identical first two entries.  Figure out how many new
    nodes there are, and allocate space for their coordinates.
    */
    *node_num2 = node_num1;

    n1_old = n2_old = -1;

    for ( edge = 0; edge < 6 * tetra_num; ++edge )
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

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tet_mesh_order10_adj_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDER10_ADJ_COUNT counts the number of nodal adjacencies.
  Discussion:
    Assuming that the tet mesh is to be used in a finite element
    computation, we declare that two distinct nodes are "adjacent" if and
    only if they are both included in some tetrahedron.
    It is the purpose of this routine to determine the number of
    such adjacency relationships.
    The initial count gets only the (I,J) relationships, for which
    node I is strictly less than node J.  This value is doubled
    to account for symmetry.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int TET_NUM, the number of tetrahedrons.
    Input, int TET_NODE[10*TET_NUM], the indices of the nodes.
    Output, int *ADJ_NUM, the total number of adjacency relationships,
    Output, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
*/
{
	const _2dt3pi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ tet_num = s_data->a1;
	int * tet_node = s_data->a2;
	int * adj_num = s_data->a3;
	int * adj_row = s_data->a4;
	
    dim_typ i, j, k;
    dim_typ l;
    int node;
    int *pair;
    dim_typ pair_num;
    dim_typ pair_unique_num;
    /*
    Each order 10 tetrahedron defines 45 adjacency pairs.
    */
    pair = ( int * ) malloc ( 45 * tet_num * sizeof ( int ) << 1 );

    k = 0;
    for ( i = 0; i < 9; ++i )
        for ( j = i + 1; j < 10; ++j )
        {
            #pragma omp parallel for num_threads(tet_num)
            for ( l = 0; l < tet_num; ++l )
            {
                pair[0+(k*tet_num+l)*2] = tet_node[i+l*10];
                pair[1+(k*tet_num+l)*2] = tet_node[j+l*10];
            }
            ++ k;
        }
    /*
    Force the nodes of each pair to be listed in ascending order.
    */
    pair_num = 45 * tet_num;

    i4col_sort2_a ( 2, pair_num, pair );
    /*
    Rearrange the columns in ascending order.
    */
    i4col_sort_a ( 2, pair_num, pair );
    /*
    Get the number of unique columns.
    */
    pair_unique_num = i4col_sorted_unique_count ( 2, pair_num, pair );
    /*
    The number of adjacencies is TWICE this value, plus the number of nodes.
    */
    *adj_num = pair_unique_num<<1;
    /*
    Now set up the ADJ_ROW counts.
    */
    for ( node = 0; node < node_num; ++node)
        adj_row[node] = 0;

    for ( k = 0; k < pair_num; ++k )
    {
        if (0<k &&  pair[0+(k-1)*2] == pair[0+k*2] &&pair[1+(k-1)*2] == pair[1+k*2] )
            continue;
        i = pair[0+(k<<1)];
        j = pair[1+(k<<1)];

        ++ adj_row[i-1];
        ++ adj_row[j-1];
    }
    /*
    We used ADJ_ROW to count the number of entries in each row.
    Convert it to pointers into the ADJ array.
    */
    for ( node = node_num-1; 0 <= node; --node )
        adj_row[node] = adj_row[node+1];

    adj_row[0] = 1;
    for ( node = 1; node <= node_num; ++node)
        adj_row[node] = adj_row[node-1] + adj_row[i];

    free ( pair );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tet_mesh_order10_adj_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    TET_MESH_ORDER10_ADJ_SET sets the nodal adjacency matrix.
  Discussion:
    A compressed format is used for the nodal adjacency matrix.
    It is assumed that we know ADJ_NUM, the number of adjacency entries
    and the ADJ_ROW array, which keeps track of the list of slots
    in ADJ where we can store adjacency information for each row.
    We essentially repeat the work of TET_MESH_ORDEr8_ADJ_COUNT, but
    now we have a place to store the adjacency information.
    A copy of the ADJ_ROW array is useful, as we can use it to keep track
    of the next available entry in ADJ for adjacencies associated with
    a given row.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 March 2013
  Author:
    John Brkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int TET_NUM, the number of tetrahedrons.
    Input, int TET_NODE[10*TET_NUM], the indices of the nodes.
    Input, int ADJ_NUM, the total number of adjacency relationships,
    Input, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
    Output, int TET_MESH_ORDEr8_ADJ_SET[ADJ_NUM],
    the adjacency information.
*/
{
	const _2dtpiipi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ tet_num = s_data->a1;
	int * tet_node = s_data->a2;
	int adj_num = s_data->a3;
	int * adj_row = s_data->a4;
	
    int *adj;
    int *adj_row_copy;
    dim_typ i, j, k;
    dim_typ l;
    int node;
    int *pair;
    dim_typ pair_num;
    /*
    Each order 10 tetrahedron defines 45 adjacency pairs.
    */
    pair = ( int * ) malloc ( 90 * tet_num * sizeof ( int ) );

    k = 0;
    for ( i = 0; i < 9; ++i)
        for ( j = i + 1; j < 10; ++j )
        {
            #pragma omp parallel for num_threads(tet_num)
            for ( l = 0; l < tet_num; ++l )
            {
                pair[0+(k*tet_num+l)*2] = tet_node[i+l*10];
                pair[1+(k*tet_num+l)*2] = tet_node[j+l*10];
            }
            ++ k;
        }
    /*
    Force the nodes of each pair to be listed in ascending order.
    */
    pair_num = 45 * tet_num;
    i4col_sort2_a ( 2, pair_num, pair );
    /*
    Rearrange the columns in ascending order.
    */
    i4col_sort_a ( 2, pair_num, pair );
    /*
    Mark all entries of ADJ so we will know later if we missed one.
    */
    adj = ( int * ) malloc ( adj_num * sizeof ( int ) );

    for ( i = 0; i < adj_num; ++i)
        adj[i] = -1;
    /*
    Copy the ADJ_ROW array and use it to keep track of the next
    free entry for each row.
    */
    adj_row_copy = ( int * ) malloc ( node_num * sizeof ( int ) );

    for ( node = 0; node < node_num; ++node)
        adj_row_copy[node] = adj_row[node];
    /*
    Now set up the ADJ_ROW counts.
    */
    for ( k = 0; k < pair_num; ++k )
    {
        if (0<k && pair[0+(k-1)*2] == pair[0+k*2] &&pair[1+(k-1)*2] == pair[1+k*2] )
            continue;
        i = pair[0+(k<<1)];
        j = pair[1+(k<<1)];

        adj[adj_row_copy[i]] = j;
        ++ adj_row_copy[i];
        adj[adj_row_copy[j]] = i;
        ++ adj_row_copy[j];
    }
    free ( adj_row_copy );
    free ( pair );

    return adj;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tetrahedron_order4_physical_to_reference ( void * data)
/*****************************************************************************80
/*
  Purpose:
    TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE maps physical points to reference points.
  Discussion:
    Given the vertices of an order 4 physical tetrahedron and a point
 (X,Y,Z) in the physical tetrahedron, the routine computes the value
    of the corresponding image point (R,S,T) in reference space.
    This routine may be appropriate for an order 10 tetrahedron,
    if the mapping between reference and physical space is linear.
    This implies, in particular, that the edges of the image tetrahedron
    are straight, the faces are flat, and the "midside" nodes in the
    physical tetrahedron are halfway along the sides of the physical
    tetrahedron.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 December 2006
  Author:
    John Burkardt
  Parameters:
    Input, double TETRA[3*4], the coordinates of the vertices.
    The vertices are assumed to be the images of
 (0,0,0), (1,0,0), (0,1,0) and (0,0,1) respectively.
    Input, int N, the number of points to transform.
    Input, double PHY[3*N], the coordinates of physical points
    to be transformed.
    Output, double REF[3*N], the coordinates of the corresponding
    points in the reference space.
*/
{
	const dt3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * phy = s_data->a1;
	ityp * ref = s_data->a2; 
	ityp * tetra = s_data->a3;
	
    ityp a[9];
    ityp det;
    dim_typ i, j;
    /*
    Set up the matrix.
    */
    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i)
        a[i+0*3] = tetra[i+1*3] - tetra[i+0*3];a[i+1*3] = tetra[i+2*3] - tetra[i+0*3];a[i+2*3] = tetra[i+3*3] - tetra[i+0*3];
    /*
    Compute the determinant.
    */
    det =  a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )+ a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )+ a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );
    /*
    If the determinant is zero, bail out.
    */
    if ( det == 0.00 )
    {
        for ( j = 0; j < n; ++j)
            #pragma omp parallel for num_threads(3)
            for ( i = 0; i < 3; ++i)
                ref[i+j*3] = 0.0;
        return NULL;
    }
    /*
    Compute the solution.
    */
    for ( j = 0; j < n; ++j )
    {
        ref[0+j*3] = ( ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )* ( phy[0+j*3] - tetra[0+0*3] )- ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] )* ( phy[1+j*3] - tetra[1+0*3] )+ ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] )* ( phy[2+j*3] - tetra[2+0*3] )) / det;
        ref[1+j*3] = ( - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] )* ( phy[0+j*3] - tetra[0+0*3] )+ ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] )* ( phy[1+j*3] - tetra[1+0*3] )- ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] )* ( phy[2+j*3] - tetra[2+0*3] )) / det;
        ref[2+j*3] = ( ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] )* ( phy[0+j*3] - tetra[0+0*3] )- ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] )* ( phy[1+j*3] - tetra[1+0*3] )+ ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] )* ( phy[2+j*3] - tetra[2+0*3] )) / det;
    }
    return NULL;
}

#endif
