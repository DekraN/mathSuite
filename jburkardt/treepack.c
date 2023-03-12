#ifndef __DISABLEDEEP_TREEPACK

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _vec_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    VEC_NEXT generates all N-vectors of integers modulo a given base.
  Discussion:
    The items are produced one at a time.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    14 July 2013
  Parameters:
    Input, int N, the size of the vectors to be used.
    Input, int IBASE, the base to be used.  IBASE = 2 will
    give vectors of 0's and 1's, for instance.
    Input/output, int IARRAY[N].  On each return,
    IARRAY will contain entries in the range 0 to IBASE-1.
    Input/output, int *MORE.  Set this variable 0 before
    the first call.  Normally, MORE will be returned 1 but
    once all the vectors have been generated, MORE will be
    reset 0 and you should stop calling the program.
*/
{
	const dtipipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int ibase = s_data->a1;
	int * iarray = s_data->a2;
	bool * more = s_data->a3;
	
	dim_typ i;
	static dim_typ kount = 0;
	static dim_typ last = 0;
	dim_typ nn;
	
	if ( ! ( *more ) )
	{
		kount = 1;
		last = powi ( ibase, n );
		*more = true; 
		for ( i = 0; i < n; ++i )
			iarray[i] = 0;
	}
	else
	{
		++ kount;
	
		if ( kount == last )
			*more = false; 
	
		++ iarray[n-1];
	
		for ( i = 1; i <= n; ++i )
		{
			nn = n - i;
			
			if ( iarray[nn] < ibase )
				return NULL;
		
			iarray[nn] = 0;
	
			if ( nn != 0 )
				++ iarray[nn-1];
		}
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _vec_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    VEC_RANDOM selects a random N-vector of integers modulo a given base.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    21 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the vector to be generated.
    Input, int BASE, the base to be used.
    Input/output, int *SEED, a random number seed.
    Output, int A[N], a list of N random values between
    0 and BASE-1.
*/
{
	const dti2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int base = s_data->a1;
	int * seed = s_data->a2;
	int * a = s_data->a3;
	
	for (dim_typ i = 0; i < n; ++i )
		a[i] = i4_uniform_ab ( 0, base-1, seed );
		
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _graph_adj_edge_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAPH_ADJ_EDGE_COUNT counts the number of edges in a graph.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int ADJ[NNODE*NNODE], the adjacency information.
    ADJ(I,J) is 1 if there is an edge from node I to node J.
    Input, int NNODE, the number of nodes.
    Output, int GRAPH_ADJ_EDGE_COUNT, the number of edges in the graph.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	
	const register dim_typ nnode = s_data->a0;
	int * adj = s_data->a1;
	
    dim_typ i, j;
    dim_typ nedge = 0;

    for ( i = 0; i < nnode; ++i )
        for ( j = 0; j < nnode; ++j )
            nedge += adj[i+j*nnode]<<(i == j);
	
	result = nedge >> 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _graph_adj_is_node_connected ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAPH_ADJ_IS_NODE_CONNECTED determines if a graph is nodewise connected.
  Discussion:
    A graph is nodewise connected if, from every node, there is a path
    to any other node.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int ADJ[NNODE*NNODE], the adjacency matrix for the
    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
    Input, int NNODE, the number of nodes.
    Output, int GRAPH_ADJ_IS_NODE_CONNECTED.
    0, the graph is not nodewise connected.
    1, the graph is nodewise connected.
*/
{
	static bool _result = 2;
	
	const dtpi * const s_data = data;
	
	const register dim_typ nnode = s_data->a0;
	int * adj = s_data->a1;
	
    bool *found;
    dim_typ i;
    dim_typ ihi;
    dim_typ ii;
    dim_typ ilo;
    dim_typ j;
    dim_typ jhi;
    dim_typ jlo;
    int *list;
    dim_typ result;
    /*
    FOUND(I) is 1 if node I has been reached.
    LIST(I) contains a list of the nodes as they are reached.
    */
    list = ( int * ) malloc ( nnode * sizeof ( int ) );
    found = ( bool * ) malloc ( nnode * sizeof ( bool ) );

    for ( i = 0; i < nnode; ++i )
        list[i] = found[i] = 0;
    /*
    Start at node 1.
    */
    found[0] = list[0] = ilo = ihi = 1;
    /*
    From the batch of nodes found last time, LIST(ILO:IHI),
    look for unfound neighbors, and store their indices in LIST(JLO:JHI).
    */
    for ( ; ; )
    {
        jlo = ihi + 1;
        jhi = ihi;

        for ( ii = ilo; ii <= ihi; ++ii )
        {
            i = list[ii-1];

            for ( j = 1; j <= nnode; ++j )
                if ( adj[i-1+(j-1)*nnode] != 0 || adj[j-1+(i-1)*nnode] != 0 )
                    if ( found[j-1] == false )
                    {
                        ++ jhi;
                        list[jhi-1] = j;
                        found[j-1] = true;
                    }
        }
        /*
        If no neighbors were found, exit.
        */
        if ( jhi < jlo )
            break;
        /*
        If neighbors were found, then go back and find THEIR neighbors.
        */
        ilo = jlo;
        ihi = jhi;
    }
    /*
    No more neighbors were found.  Have we reached all nodes?
    */
    result = ihi == nnode;

    free ( found );
    free ( list );

	_result = result;
    return &_result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _graph_adj_is_tree ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAPH_ADJ_IS_TREE determines whether a graph is a tree.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int ADJ[NNODE*NNODE], the adjacency matrix for the
    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
    Input, int NNODE, the number of nodes.
    Output, int GRAPH_ADJ_IS_TREE.
    0, the graph is not a tree.
    1, the graph is a tree.
*/
{
	static bool _result = 2;
	
	const dtpi * const s_data = data;
	
	const register dim_typ nnode = s_data->a0;
	int * adj = s_data->a1;
	
    dim_typ nedge;
    bool result;

    if ( nnode <= 1 )
    {
    	_result = 1;
        return &_result;
    }
    /*
    Every node must be connected to every other node.
    */
    result = graph_adj_is_node_connected ( adj, nnode );
    /*
    There must be exactly NNODE-1 edges.
    */

	_result = result == 0 ? result : graph_adj_edge_count ( adj, nnode ) == nnode-1;
    return &_result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _graph_arc_degree ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAPH_ARC_DEGREE determines the degree of the nodes of a graph.
  Discussion:
    The degree of a node is the number of edges that include the node.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int NEDGE, the number of edges.
    Input, int INODE[NEDGE], JNODE[NEDGE], the pairs of nodes
    that form the edges.
    Output, int GRAPH_ARC_DEGREE[NNODE], the degree of each node,
    that is, the number of edges that include the node.
*/
{
	const _2dt2pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	const register dim_typ nedge = s_data->a1;
	int * inode = s_data->a2;
	int * jnode = s_data->a3;
	
    int *degree;
    dim_typ i, n;

    degree = ( int * ) malloc ( nnode * sizeof ( int ) );

    for ( i = 0; i < nnode; ++i )
        degree[i] = 0;

    for ( i = 0; i < nedge; ++i )
    {
        n = inode[i];

        if ( 1 <= n && n <= nnode )
            ++ degree[n-1];
        else
            return NULL;

        n = jnode[i];
        if ( 1 <= n && n <= nnode )
            ++ degree[n-1];
        else
            return NULL;
    }
    return degree;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _graph_arc_is_tree ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAPH_ARC_IS_TREE determines whether a graph is a tree.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and
    JNODE(I) are the start and end nodes of the I-th edge of the graph G.
    Input, int NEDGE, the number of edges in the graph G.
    Input, int NNODE, the number of nodes.
    Output, int GRAPH_ARC_IS_TREE.
    0, the graph is not a tree.
    1, the graph is a tree.
*/
{
	static bool _result = 2;
	
	const _2dt2pi * const s_data = data;
	
	const register dim_typ nedge = s_data->a0;
	const register dim_typ nnode = s_data->a1;
	int * inode = s_data->a2;
	int * jnode = s_data->a3;
	
    int *adj;
    bool result;
    adj = graph_arc_to_graph_adj ( nedge, inode, jnode );
    result = graph_adj_is_tree ( adj, nnode );
    free ( adj );
    
    _result = result;
    return &_result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _graph_arc_node_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAPH_ARC_NODE_COUNT counts the number of nodes in a graph.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NEDGE, the number of edges.
    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and
    JNODE(I) are the start and end nodes of the I-th edge.
    Output, int GRAPH_ARC_NODE_COUNT, the number of distinct nodes.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dt2pi * const s_data = data;
	const register dim_typ nedge = s_data->a0;
	int * inode = s_data->a1;
	int * jnode = s_data->a2;
	
    dim_typ i;
    int *knode;
    dim_typ nnode;
    /*
    Copy all the node labels into KNODE,
    sort KNODE,
    count the unique entries.

    That's NNODE.
    */
    knode = ( int * ) malloc ( nedge * sizeof ( int ) << 1 );

    for ( i = 0; i < nedge; ++i )
    {
        knode[i] = inode[i];
        knode[i+nedge] = jnode[i];
    }

    i4vec_sort_heap_a ( 2*nedge, knode );
    nnode = i4vec_sorted_unique_count ( 2*nedge, knode );
    free ( knode );
    
    result = nnode;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _graph_arc_node_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAPH_ARC_NODE_MAX determines the maximum node label.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NEDGE, the number of edges.
    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and
    JNODE(I) are the start and end nodes of the I-th edge.
    Output, int GRAPH_ARC_NODE_MAX, the maximum node index.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dt2pi * const s_data = data;
	const register dim_typ nedge = s_data->a0;
	int * inode = s_data->a1;
	int * jnode = s_data->a2;
	
    dim_typ node_max = 0;
    for (dim_typ i = 0; i < nedge; ++i)
        node_max = node_max < inode[i] ? inode[i]:jnode[i];
        
    result = node_max;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _graph_arc_to_graph_adj ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAPH_ARC_TO_GRAPH_ADJ converts an arc list graph to an adjacency graph.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NEDGE, the number of edges.
    Input, int INODE[NEDGE], JNODE[NEDGE], the edge array for
    an undirected graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
    Output, int GRAPH_ARC_TO_GRAPH_ADJ[NNODE*NNODE], the adjacency information.

*/
{
	const dt2pi * const s_data = data;
	const register dim_typ nedge = s_data->a0;
	int * inode = s_data->a1;
	int * jnode = s_data->a2;
	
    int *adj;
    dim_typ i, j, k;
    dim_typ nnode;
    /*
    Determine the number of nodes.
    */
    nnode = graph_arc_node_count ( nedge, inode, jnode );
    adj = ( int * ) malloc ( nnode * nnode * sizeof ( int ) );

    for ( j = 0; j < nnode; ++j )
        adj[i+j*nnode] = 0;
            for ( i = 0; i < nnode; ++i )

    for ( k = 0; k < nedge; ++k)
    {
        i = inode[k] - 1;
        j = jnode[k] - 1;
        adj[i+j*nnode] = adj[j+i*nnode] = 1;
    }

    return adj;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _pruefer_to_tree_arc ( void * data)
/******************************************************************************/
/*
  Purpose:
    PRUEFER_TO_TREE_ARC is given a Pruefer code, and computes the tree.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 July 2013
  Reference:
    Dennis Stanton, Dennis White,
    Constructive Combinatorics,
    Springer Verlag, New York, 1986.
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int IARRAY[NNODE-2], the Pruefer code of the tree.
    Output, int INODE[NNODE-1], JNODE[NNODE-1], the edge
    array of the tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
*/
{
	const dt3pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * iarray = s_data->a1;
	int * inode = s_data->a2;
	int * jnode = s_data->a3;
	
    dim_typ i;
    short ii;
    int *iwork;
    dim_typ j;
    /*
    Initialize IWORK(I) to count the number of neighbors of node I.
    The Pruefer code uses each node one less time than its total
    number of neighbors.
    */
    iwork = ( int * ) malloc ( nnode * sizeof ( int ) );

    for ( i = 0; i < nnode; ++i )
        iwork[i] = 1;
    for ( i = 0; i < nnode - 2; ++i )
        ++ iwork[iarray[i]-1];

    for ( i = 0; i < nnode - 1; ++i )
        inode[i] = jnode[i] = -1;
    /*
    Now process each entry in the Pruefer code.
    */
    for ( i = 0; i < nnode - 2; ++i )
    {
        ii = -1;
        for ( j = 0; j < nnode; ++j )
            if ( iwork[j] == 1 )
                ii = j;
        inode[i] = ii + 1;
        jnode[i] = iarray[i];
        iwork[ii] = 0;
        -- iwork[iarray[i]-1];
    }

    inode[nnode-2] = iarray[nnode-3];
    jnode[nnode-2] = 1 + iarray[nnode-3] == 1;

    free ( iwork );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _pruefer_to_tree_2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    PRUEFER_TO_TREE_2 produces the edge list of a tree from its Pruefer code.
  Discussion:
    One can thus exhibit all trees on N nodes, produce
    one at random, find the M-th one on the list, etc, by
    manipulating the Pruefer codes.
    For every labeled tree on N nodes, there is a unique N-2 tuple
    of integers A1 through AN-2, with each A between 1 and N.  There
    are N^(N-2) such sequences, and each one is associated with exactly
    one tree.
    This routine apparently assumes that the Pruefer code is
    generated by taking the LOWEST labeled terminal node each time.
    This is not consistent with PRUEFER_TO_TREE and TREE_TO_PRUEFER.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 June 2013
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis. Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int NNODE, number of nodes in desired tree.
    Input, int IARRAY[NNODE].  IARRAY(I), I = 1, NNODE-2
    is the Pruefer code for the tree.
    Output, int ITREE[NNODE]; the I-th edge of the tree
    joins nodes I and ITREE(I).
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * iarray = s_data->a1;
	int * itree = s_data->a2;
	
    dim_typ i;
    dim_typ ir;
    dim_typ j;
    dim_typ k;
    dim_typ kp;
    dim_typ l;

    for ( i = 0; i < nnode; ++i )
        itree[i] = 0;

    for ( i = nnode - 2; 1 <= i; --i )
    {
        l = iarray[i-1];

        if ( itree[l-1] == 0 )
        {
            iarray[i-1] = - l;
            itree[l-1] = - 1;
        }
    }
    iarray[nnode-2] = nnode;
    /*
    Find next index K so that ITREE(K) is 0.
    */
    k = 1;


    while ( itree[k-1] != 0 )
        ++ k;

    j = 0;
    kp = k;

    for ( ; ; )
    {
        ++ j;
        ir = abs ( iarray[j-1] );
        itree[kp-1] = ir;

        if ( j == nnode - 1 )
            break;

        if ( 0 < iarray[j-1] )
        {
            while ( itree[k-1] != 0 )
                ++ k;
            kp = k;
            continue;
        }

        if ( k < ir )
        {
            itree[ir-1] = 0;
            while ( itree[k-1] != 0 )
                ++ k;
            kp = k;
            continue;
        }
        kp = ir;
    }
    /*
    Restore the signs of IARRAY.
    */
    for ( i = 0; i < nnode - 2; ++i )
        iarray[i] = abs ( iarray[i] );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _pruefer_to_tree_2_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    PRUEFER_TO_TREE_2_NEW produces the edge list of a tree from its Pruefer code.
  Discussion:
    One can thus exhibit all trees on N nodes, produce
    one at random, find the M-th one on the list, etc, by
    manipulating the Pruefer codes.
    For every labeled tree on N nodes, there is a unique N-2 tuple
    of integers A1 through AN-2, with each A between 1 and N.  There
    are N^(N-2) such sequences, and each one is associated with exactly
    one tree.
    This routine apparently assumes that the Pruefer code is
    generated by taking the LOWEST labeled terminal node each time.
    This is not consistent with PRUEFER_TO_TREE and TREE_TO_PRUEFER.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 June 2013
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis. Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int NNODE, number of nodes in desired tree.
    Input, int IARRAY[NNODE].  IARRAY(I), I = 1, NNODE-2
    is the Pruefer code for the tree.
    Output, int PRUEFER_TO_TREE_2_NEW[NNODE]; the I-th edge of the tree
    joins nodes I and ITREE(I).
*/
{
	const dtpi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * iarray = s_data->a1;
	
    int *itree;
    itree = ( int * ) malloc ( nnode * sizeof ( int ) );
    pruefer_to_tree_2 ( nnode, iarray, itree );
    return itree;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tree_arc_center ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ARC_CENTER computes the center, eccentricity, and parity of a tree.
  Discussion:
    A tree is an undirected graph of N nodes, which uses N-1 edges,
    and is connected.
    A graph with N-1 edges is not guaranteed to be a tree, and so this
    routine must first check that condition before proceeding.
    The edge distance between two nodes I and J is the minimum number of
    edges that must be traversed in a path from I and J.
    The eccentricity of a node I is the maximum edge distance between
    node I and the other nodes J in the graph.
    The radius of a graph is the minimum eccentricity over all nodes
    in the graph.
    The diameter of a graph is the maximum eccentricity over all nodes
    in the graph.
    The center of a graph is the set of nodes whose eccentricity is
    equal to the radius, that is, the set of nodes of minimum eccentricity.
    For a tree, the center is either a single node, or a pair of
    neighbor nodes.
    The parity of the tree is 1 if the center is a single node, or 2 if
    the center is 2 nodes.
    The center of a tree can be found by removing all "leaves", that is,
    nodes of degree 1.  This step is repeated until only 1 or 2 nodes
    are left.
    Thanks to Alexander Sax for pointing out that a previous version of the
    code was failing when the tree had an odd parity, that is, a single
    center node, 15 April 2013.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int INODE[NNODE-1], JNODE[NNODE-1], the edges of
    the tree.  Edge I connects nodes INODE(I) and JNODE(I).
    Output, int CENTER[2].  CENTER(1) is the index of the
    first node in the center.  CENTER(2) is 0 if there is only one node
    in the center, or else the index of the second node.
    Output, int *ECCENT, the eccentricity of the nodes in
    the center, and the radius of the the tree.
    Output, int *PARITY, the parity of the tree, which is
    normally 1 or 2.
*/
{
	const dt5pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * inode = s_data->a1;
	int * jnode = s_data->a2;
	int * center = s_data->a3;
	int * eccent = s_data->a4;
	int * parity = s_data->a5;
	
    int *degree;
    dim_typ i;
    dim_typ iedge;
    dim_typ ileaf;
    dim_typ j;
    int *list;
    dim_typ nedge;
    dim_typ nleaf;
    dim_typ nnode2;
    dim_typ result;

    *eccent = center[0] = center[1] = *parity = 0;

    if ( nnode <= 0 )
        return NULL;
    else if ( nnode == 1 )
    {
        *eccent = center[1] = 0;
        center[0] = *parity = 1;
        return NULL;
    }
    else if ( nnode == 2 )
    {
        *eccent = center[0] = 1;
        center[1] = *parity = 2;
        return NULL;
    }
    /*
    Is this graph really a tree?
    */
    nedge = nnode - 1;
    result = graph_arc_is_tree ( nedge, inode, jnode, nnode );

    if ( result == 0 )
        return NULL;
    /*
    Compute the degrees.
    */
    degree = graph_arc_degree ( nnode, nedge, inode, jnode );
    /*
    Defoliate the tree.
    */
    nnode2 = nnode;
    list = ( int * ) malloc ( nnode * sizeof ( int ) );

    for ( ; ; )
    {
        ++ *eccent;
        /*
        Find and mark the leaves.
        */
        nleaf = 0;

        for ( i = 1; i <= nnode; ++i )
            if ( degree[i-1] == 1 )
            {
                ++ nleaf;
                list[nleaf-1] = i;
            }
        /*
        Delete the leaves.
        */
        for ( ileaf = 1; ileaf <= nleaf; ++ileaf )
        {
            i = list[ileaf-1];

            iedge = j = 0;

            for ( ; ; )
            {
                ++ iedge;

                if ( nedge < iedge )
                    return NULL;

                if ( inode[iedge-1] == i )
                {
                    j = jnode[iedge-1];
                    inode[iedge-1] *= -1;
                    jnode[iedge-1] *= -1;
                }
                else if ( jnode[iedge-1] == i )
                {
                    j = inode[iedge-1];
                    inode[iedge-1] *= -1;
                    jnode[iedge-1] *= -1;
                }

                if ( j != 0 )
                    break;
            }

            degree[i-1] = -1;
            -- nnode2;
            -- degree[j-1];
            /*
            If the other node has degree 0, we must have just finished
            stripping all leaves from the tree, leaving a single node.
            Don't kill it here.  It is our odd center.

            if ( degree(j) == 0 )
            {
            nnode2 = nnode2 - 1
            }
            */
        }
        /*
        Find the remaining nodes.
        */
        nnode2 = 0;

        for ( i = 1; i <= nnode; ++i )
            if ( 0 <= degree[i-1] )
            {
                ++ nnode2;
                list[nnode2-1] = i;
            }
        /*
        If at least 3, more pruning is required.
        */
        if ( nnode2 < 3 )
            break;
    }
    /*
    If only one or two nodes left, we are done.
    */
    *parity = nnode2;

    for ( i = 0; i < nnode2; ++i )
        center[i] = list[i];
    for ( i = 0; i < nedge; ++i)
    {
        inode[i] = abs ( inode[i] );
        jnode[i] = abs ( jnode[i] );
    }

    free ( list );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tree_arc_diam ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ARC_DIAM computes the "diameter" of a tree.
  Discussion:
    A tree is an undirected graph of N nodes, which uses N-1 edges,
    and is connected.
    A graph with N-1 edges is not guaranteed to be a tree, and so this
    routine must first check that condition before proceeding.
    The diameter of a graph is the length of the longest possible
    path that never repeats an edge.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int INODE[NNODE-1], JNODE[NNODE-1], the edges
    of the tree.  Edge I connects nodes INODE(I) and JNODE(I).
    Output, int *DIAM, the length of the longest path
    in the tree.
    Output, int LABEL[NNODE], marks the path between
    nodes N1 and N2.  Node I is in this path if LABEL(I) is 1.
    Output, int *N1, *N2, the indices of two nodes in the
    tree which are separated by DIAM edges.
*/
{
	const dt6pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * inode = s_data->a1;
	int * jnode = s_data->a2;
	int * diam = s_data->a3;
	int * label = s_data->a4;
	int * n1 = s_data->a5;
	int * n2 = s_data->a6;
	
    int *degree;
    dim_typ i;
    dim_typ invals;
    dim_typ j;
    dim_typ k;
    dim_typ kstep;
    dim_typ nabe;
    dim_typ nedge;
    dim_typ result;

    if ( nnode <= 0 )
    {
        *diam = 0;
        return NULL;
    }

    if ( nnode == 1 )
    {
        *diam = 0;
        *n1 = *n2 = 1;
        return NULL;
    }
    /*
    Is this graph really a tree?
    */
    nedge = nnode - 1;
    result = graph_arc_is_tree ( nedge, inode, jnode, nnode );

    if ( result == 0 )
        return NULL;

    label = i4vec_indicator_new ( nnode );
    /*
    On step K:
    Identify the terminal and interior nodes.
    If there are no interior nodes left,
    then there are just two nodes left at all.  The diameter is 2*K-1,
    and a maximal path extends between the nodes whose labels are
    contained in the two remaining terminal nodes.
    Else
    The label of each terminal node is passed to its interior neighbor.
    If more than one label arrives, take any one.
    The terminal nodes are removed.
    */
    kstep = 0;
    degree = ( int * ) malloc ( nnode * sizeof ( int ) );

    for ( ; ; )
    {
        ++ kstep;
        /*
        Compute the degree of each node.
        */
        for ( j = 1; j <= nnode; ++j )
            degree[j-1] = 0;
        for ( j = 1; j <= nedge; ++j )
        {
            k = inode[j-1];
            if ( 0 < k )
                ++ degree[k-1];

            k = jnode[j-1];
            if ( 0 < k )
                ++ degree[k-1];
        }
        /*
        Count the number of interior nodes.
        */
        invals = 0;
        for ( i = 1; i <= nnode; ++i )
                if ( 1 < degree[i-1] )
                    ++ invals;
        /*
        If there are 1 or 0 interior nodes, it's time to stop.
        */
        if ( invals == 1 )
        {
            *diam = kstep<<1;
            break;
        }
        else if ( invals == 0 )
        {
            *diam = (kstep<<1) - 1;
            break;
        }
        /*
        If there are at least two interior nodes, then chop off the
        terminal nodes and pass their labels inward.
        */
        for ( k = 1; k <= nnode; ++k )
        {
            if ( degree[k-1] == 1 )
            {
                for ( j = 1; j <= nedge; ++j )
                {
                    if ( inode[j-1] == k )
                    {
                        nabe = jnode[j-1];
                        label[nabe-1] = label[k-1];
                        inode[j-1] *= -1;
                        jnode[j-1] *= -1;
                    }
                    else if ( jnode[j-1] == k )
                    {
                        nabe = inode[j-1];
                        label[nabe-1] = label[k-1];
                        inode[j-1] *= -1;
                        jnode[j-1] *= -1;
                    }
                }
            }
        }
    }
    /*
    Now get the labels from two of the remaining terminal nodes.
    The nodes represented by these labels will be a diameter apart.
    */
    *n1 = *n2 = 0;

    for ( i = 1; i <= nnode; ++i );
        if ( degree[i-1] == 1)
            if ( *n1 == 0 )
                *n1 = label[i-1];
            else if ( *n2 == 0 )
                *n2 = label[i-1];
    /*
    Set the labels of the interior node (if any) and nodes marked
    N1 and N2 to 1, and all others to 0.  This will label the nodes on the path.
    */
    if ( invals == 1 )
        for ( i = 1; i <= nnode; ++i )
            if ( 1 < degree[i-1] )
            label[i-1] = 1;

    for ( i = 1; i <= nnode; ++i )
        label[i-1] = label[i-1] == *n1 || label[i-1] == *n2;
    /*
    Clean up the arrays.
    */
    for ( j = 1; j <= nedge; ++j )
    {
        inode[j-1] = abs ( inode[j-1] );
        k = inode[j-1];
        jnode[j-1] = abs ( jnode[j-1] );
        k = jnode[j-1];
    }

    free ( degree );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tree_arc_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ARC_RANDOM selects a random labeled tree and its Pruefer code.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes.
    Input/output, int *SEED, a seed for the random number
    generator.
    Output, int CODE[NNODE-2], the Pruefer code for the
    labeled tree.
    Output, int INODE[NNODE-1], JNODE[NNODE-1], the edge
    array for the tree.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * seed = s_data->a1;
	int * code = s_data->a2;
	int * inode = s_data->a3;
	int * jnode = s_data->a4;
	
    if ( nnode <= 2 )
        return NULL;
    vec_random ( nnode-2, nnode, seed, code );
    for (dim_typ i = 0; i < nnode - 2; i++, ++code[i]);
    pruefer_to_tree_arc ( nnode, code, inode, jnode );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tree_arc_to_pruefer ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ARC_TO_PRUEFER is given a labeled tree, and computes its Pruefer code.
  Discussion:
    A tree is an undirected graph of N nodes, which uses N-1 edges,
    and is connected.
    A graph with N-1 edges is not guaranteed to be a tree, and so this
    routine must first check that condition before proceeding.
    The Pruefer code is a correspondence between all labeled trees of
    N nodes, and all list of N-2 integers between 1 and N (with repetition
    allowed).  The number of labeled trees on N nodes is therefore N^(N-2).
    The Pruefer code is constructed from the tree as follows:
    A terminal node on the tree is defined as a node with only one neighbor.
    Consider the set of all terminal nodes on the tree.  Take the one
    with the highest label, I.  Record the label of its neighbor, J.
    Delete node I and the edge between node I and J.
    J is the first entry in the Pruefer code for the tree.   Repeat
    the operation a total of N-2 times to get the complete code.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2013
  Reference:
    Dennis Stanton, Dennis White,
    Constructive Combinatorics,
    Springer Verlage, New York, 1986.
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int INODE[NNODE-1], JNODE[NNODE-1], the edge array
    of the tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
    Output, int TREE_ARC_TO_PRUEFER[NNODE-2], the Pruefer code of the tree.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * inode = s_data->a1;
	int * jnode = s_data->a2;
	
    int *code;
    int *degree;
    dim_typ i;
    dim_typ i2;
    dim_typ iterm;
    dim_typ j;
    dim_typ jsave;
    dim_typ nedge;
    dim_typ result;
    /*
    Is this graph really a tree?
    */
    nedge = nnode - 1;
    result = graph_arc_is_tree ( nedge, inode, jnode, nnode );

    if ( result == 0 )
        return NULL; 

    code = ( int * ) malloc ( ( nnode - 2 ) * sizeof ( int ) );
    /*
    Compute the degree of each node.
    */
    nedge = nnode - 1;
    degree =  graph_arc_degree ( nnode, nedge, inode, jnode );
    /*
    Compute the next term of the Pruefer code.
    */
    for ( i = 1; i <= nnode - 2; ++i )
    {
        /*
        Find the terminal node with the highest label.
        */
        iterm = 0;

        for ( j = 1; j <= nnode; ++j )
            if ( degree[j-1] == 1 )
                iterm = j;
        /*
        Find the edge that includes this node, and note the
        index of the other node.
        */
        for ( j = 1; j < nnode - 1; ++j )
        {
            jsave = j;

            if ( inode[j-1] == iterm )
            {
                i2 = 2;
                break;
            }
            else if ( jnode[j-1] == iterm )
            {
                i2 = 1;
                break;
            }
        }
        /*
        Delete the edge from the tree.
        */
        -- degree[inode[jsave-1]-1];
        -- degree[jnode[jsave-1]-1];
        /*
        Add the neighbor of the node to the Pruefer code.
        */
        code[i-1] = i2 == 1 ? inode[jsave-1] : jnode[jsave-1];
        /*
        Negate the nodes in the edge list to mark them as used.
        */
        inode[jsave-1] *= -1;
        jnode[jsave-1] *= -1;
    }
    /*
    Before returning, restore the original form of the edge list.
    */
    for ( i = 1; i <= nnode - 1; ++i )
    {
        inode[i-1] = abs ( inode[i-1] );
        jnode[i-1] = abs ( jnode[i-1] );
    }

    free ( degree );

    return code;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tree_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ENUM enumerates the labeled trees on NNODE nodes.
  Discussion:
    The formula is due to Cauchy.
  Example:
    NNODE      NTREE
    0              1
    1              1
    2              1
    3              3
    4             16
    5            125
    6           1296
    7          16807
    8         262144
    9        4782969
   10      100000000
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes in each tree.
    NNODE must normally be at least 3, but for this routine,
    any value of NNODE is allowed.  Values of NNODE greater than 10
    will probably overflow.
    Output, int TREE_ENUM, the number of distinct labeled trees.
*/
{
	static unsigned result = UINT_MAX;
	
	const register unsigned nnode = *(unsigned *) nnode;
	
	result = nnode<0 ? 0 : nnode == 0 || nnode == 1 || nnode == 2 ? 1 : powi ( nnode, nnode - 2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tree_parent_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_PARENT_NEXT generates, one at a time, all labeled trees.
  Discussion:
    The routine also returns the corresponding Pruefer codes.
  Formula:
    There are N^(N-2) labeled trees on N nodes (Cayley's formula).
    The number of trees in which node I has degree D(I) is the
    multinomial coefficient: ( N-2; D(1)-1, D(2)-1, ..., D(N)-1 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2013
  Parameters:
    Input, int NNODE, the number of nodes to be used in
    the trees.
    Input/output, int CODE[NNODE].  The first NNODE-2 entries
    of CODE contain the Pruefer code for the given labeled tree.
    Output, int ITREE[NNODE].  The first NNODE-1 entries
    of ITREE describe the edges that go between the nodes.  Each pair
 (I, ITREE(I)) represents an edge.  Thus if ITREE(5) = 3,
    there is an edge from node 3 to node 5.
    Input/output, int *MORE.  On the first call only, the
    user is required to set MORE = .FALSE.  Then call the routine, and
    it will return information about the first tree
    as well as setting MORE to the value .TRUE.
    Keep calling to get another tree until MORE is .FALSE.
    on return, at which point there are no more trees.
*/
{
	const dt2pipb * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * code = s_data->a1;
	int * itree = s_data->a2;
	bool * more = s_data->a3;
	
    dim_typ i;
    if ( *more )
        for ( i = 0; i < nnode - 2; i++, --code[i] );
    vec_next ( nnode-2, nnode, code, more );
    for ( i = 0; i < nnode - 2; i++, ++code[i] );
    pruefer_to_tree_2 ( nnode, code, itree );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tree_parent_to_arc ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_PARENT_TO_ARC converts a tree from parent to arc representation.
  Discussion:
    Parent representation lists the parent node of each node.  For a
    tree of N nodes, one node has a parent of 0, representing a null link.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes in the tree.
    Input, int PARENT[NNODE], the parent node representation
    of the tree.
    Output, int *NEDGE, the number of edges, normally NNODE-1.
    Output, int INODE[NEDGE], JNODE[NEDGE], pairs of nodes
    that define the links.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * parent = s_data->a1;
	int * nedge = s_data->a2;
	int * inode = s_data->a3;
	int * jnode = s_data->a4;

    *nedge = 0;

    for (dim_typ i = 1; i <= nnode; ++i )
        if ( parent[i-1] != 0 )
        {
            ++ *nedge;
            inode[*nedge-1] = i;
            jnode[*nedge-1] = parent[i-1];
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tree_rb_lex_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_RB_LEX_NEXT generates rooted binary trees in lexicographic order.
  Discussion:
    The information definining the tree of N nodes is stored in a vector
    of 0's and 1's, in preorder traversal form.  Essentially, the
    shape of the tree is traced out with a pencil that starts at the root,
    and stops at the very last null leaf.  The first time that a (non-null)
    node is encountered, a 1 is added to the vector, and the left
    descendant of the node is visited next.  When the path returns from
    the first descendant, the second descendant is visited.  When then path
    returns again, the path passes back up from the node to its parent.
    A null leaf is encountered only once, and causes a zero to be added to
    the vector, and the path goes back up to the parent node.
    The lexicographic order is used to order the vectors of 1's and 0's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 August 2013
  Reference:
    Frank Ruskey,
    Combinatorial Generation,
    To appear.
  Parameters:
    Input, int N, the number of nodes in the rooted binary
    tree.  N should be odd.
    Input/output, int A[N], the preorder traversal form for
    the previous/next rooted binary tree.
    Output, logical *MORE, is TRUE if the next rooted binary tree was
    returned on this call, or FALSE if there are no more rooted binary
    trees, and the output of the previous call was the last one.
*/
{
	const dtpipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	bool * more = s_data->a2;
	
    dim_typ i;
    dim_typ k;
    dim_typ p;
    dim_typ q;

    if ( ! *more )
    {
        for ( i = 1; i <= n - 2; i += 2 )
            a[i-1] = 1;
        for ( i = 2; i <= n - 1; i += 2 )
            a[i-1] = 0;
        a[n-1] = 0;
        *more = true;
        return NULL;
    }
    /*
    Find the last 1 in A.
    */
    k = n;
    while ( a[k-1] == 0 )
        -- k;
    q = n - k - 1;
    /*
    Find the last 0 preceding the last 1 in A.
    If there is none, then we are done, because 11...1100..00
    is the final element.
    */
    for ( ; ; )
    {
        if ( k == 1 )
        {
            *more = false;
            return NULL;
        }

        if ( a[k-1] == 0 )
            break;
        -- k;
    }

    p = n - k - q - 1;
    a[k-1] = 1;
    for ( i = k + 1; i <= n - (p<<1) + 1; ++i )
        a[i-1] = 0;
    for ( i = n - (p<<1) + 2; i <= n - 2; i += 2 )
        a[i-1] = 1;
    for ( i = n - (p<<1) + 3; i <= n - 1; i += 2 )
        a[i-1] = 0;
    a[n-1] = 0;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tree_rb_to_parent ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_RB_TO_PARENT converts rooted binary tree to parent node representation.
  Discussion:
    Parent node representation of a tree assigns to each node a "parent" node,
    which represents the first link of the path between the node and the
    root node.  The root node itself is assigned a parent of 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of nodes in the tree.
    Input, int A[N], the preorder traversal form for the
    rooted binary tree.
    Output, int TREE_RB_TO_PARENT[N], the parent node representation
    of the tree.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ dad;
    dim_typ k;
    dim_typ node;
    dim_typ node_num;
    int *parent;
    int *use;

    parent = ( int * ) malloc ( n * sizeof ( int ) );
    use = ( int * ) malloc ( n * sizeof ( int ) );

    node = node_num = 0;

    for ( k = 1; k <= n; ++k )
    {
        dad = node;
        node_num = node_num + 1;
        node = node_num;
        parent[node-1] = dad;

        if ( a[k-1] == 1 )

            use[node-1] = 0;
        else
        {
            use[node-1] = 2;

            while ( use[node-1] == 2 )
            {
                node = dad;
                if ( node == 0 )
                    break;
                ++ use[node-1];
                dad = parent[node-1];
            }
        }
    }
    free ( use );

    return parent;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tree_rb_yule ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_RB_YULE adds two nodes to a rooted binary tree using the Yule model.
  Discussion:
    The Yule model is a simulation of how an evolutionary family tree
    develops.  We start with a root node.  The internal nodes of the tree
    are inactive and never change.  Each pendant or leaf node of the
    tree represents a biological family that can spontaneously "fission",
    developing two new distinct sub families.  In graphical terms, the node
    becomes internal, with two new leaf nodes depending from it.
    The tree is stored in inorder traversal form.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 August 2013
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N, the number of nodes in the input
    tree.  On output, this number has been increased, usually by 2.
    Input/output, int *SEED, a seed for the random number
    generator.
    Input/output, int A[*], the preorder traversal form
    for the rooted binary tree.  The number of entries in A is N.
*/
{
	const _2pipdt * const s_data = data;
	
	int * seed = s_data->a0;
	int * a = s_data->a1;
	dim_typ * n = s_data->a2;
	
    dim_typ i;
    dim_typ ileaf;
    dim_typ j;
    dim_typ jleaf;
    dim_typ nleaf;

    if ( *n <= 0 )
    {
        *n = 1;
        a[0] = 0;
        return NULL;
    }
    /*
    Count the expected number of leaves, which are the 0 values.
    */
    nleaf = ( *n + 1 ) / 2;
    /*
    Choose a random number between 1 and NLEAF.
    */
    ileaf = i4_uniform_ab ( 1, nleaf, seed );
    /*
    Locate leaf number ILEAF.
    */
    j = jleaf = 0;
    for ( i = 1; i <= *n; ++i)
    {
        if ( a[i-1] == 0 )
            ++ jleaf;
        if ( jleaf == ileaf )
        {
            j = i;
            break;
        }
    }
    /*
    Replace '0' by '100'
    */
    for ( i = *n; j <= i; --i )
        a[i+1] = a[i-1];
    a[j-1] = 1;
    a[j] = 0;

    *n += 2;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tree_rooted_code ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ROOTED_CODE returns the code of a rooted tree.
  Discussion:
    This code for a rooted tree depends on the node ordering, so it's actually
    the code for a labeled rooted tree.  To eliminate the effects of node
    labeling, one could choose as the code for a tree the maximum of all
    the codes associated with the different possible labelings of the tree.
    There are more effective ways of arriving at this code than simply
    generating all possible codes and comparing them.
    For a tree with NNODES, the code is a list of 2*NNODE 0's and 1's,
    describing a traversal of the tree starting at an imaginary node 0,
    moving "down" to the root (a code entry of 1), and then moving
    "down" (1) or "up" (0) as the tree is traversed in a depth first
    manner.  The final move must be from the root up to the imaginary
    node 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int PARENT[NNODE], is the parent node of each node.
    The node with parent 0 is the root.
    Output, int TREE_ROOTED_CODE[2*NNODE], the code for the tree.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * parent = s_data->a1;
	
    int *code;
    dim_typ father;
    dim_typ i, k;
    dim_typ son;

    code = ( int * ) malloc ( nnode * sizeof ( int ) << 1 );
    /*
    Find the root.
    */
    father = 0;
    for ( i = 1; i <= nnode; ++i )
        if ( parent[i-1] == 0 )
        {
            k = 1;
            code[0] = 1;
            father = i;
            break;
        }

    if ( father == 0 )
        return NULL;

    while ( father != 0 )
    {
        ++ k;
        code[k-1] = 0;
        for ( son = 1; son <= nnode; ++son )
            if ( parent[son-1] == father )
            {
                code[k-1] = 1;
                father = son;
                break;
            }

        if ( code[k-1] == 0 )
        {
            parent[father-1] *= -1;
            father = - parent[father-1];
        }
    }

    for ( i = 0; i < nnode; ++i )
        parent[i] *= -1;

    return code;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tree_rooted_code_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ROOTED_CODE_COMPARE compares a portion of the code for two rooted trees.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int NPART, the number of nodes for which the code
    has been determined.  This determines the portion of the codes to be
    compared.  We expect 0 <= NPART <= NNODE.
    Input, int CODE1[2*NNODE], CODE2[2*NNODE], the two
    rooted tree codes to be compared.
    Output, int TREE_ROOTED_CODE_COMPARE, the result of the comparison.
    -1, CODE1 < CODE2,
     0, CODE1 = CODE2,
    +1, CODE1 > CODE2.
*/
{
	
	static short result = SHRT_MAX;
	
	const _2dt2pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	const register dim_typ npart = s_data->a1;
	int * code1 = s_data->a2;
	int * code2 = s_data->a3;
	
    dim_typ i;
    dim_typ ihi;

    if ( npart <= 0 )
    {
    	result = 0;
        return &result;
    }

    ihi = nnode<<1;
    if ( npart < nnode )
        ihi = npart<<1;

    for ( i = 0; i < ihi; ++i )
    {
        if ( code1[i] < code2[i] )
        {
        	result = -1;
            return &result;
        }
        else if ( code2[i] < code1[i] )
        {
        	result = 1;
            return &result;
        }
    }
	
	result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tree_rooted_depth ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ROOTED_DEPTH returns the depth of a rooted tree.
  Discussion:
    The depth of any node of a rooted tree is the number of edges in
    the shortest path from the root to the node.
    The depth of the rooted tree is the maximum of the depths
    of all the nodes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 July 2013
  Author:
    John Burkardt
  Parameters:
    Input, int NNODE, the number of nodes.
    Input, int PARENT[NNODE], is the parent node of each node.
    The node with parent 0 is the root.
    Output, int *DEPTH, the depth of the tree.
    Output, int DEPTH_NODE[NNODE], the depth of each node.
*/
{
	const dt3pi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * parent = s_data->a1;
	int * depth = s_data->a2;
	int * depth_node = s_data->a3;
	
    dim_typ i, j;
    short root;
    /*
    Find the root.
    */
    root = -1;
    for ( i = 1; i <= nnode; ++i)
        if ( parent[i-1] == 0 )
        {
            root = i;
            break;
        }

    if ( root == -1 )
        return NULL;
    /*
    Determine the depth of each node by moving towards the node.
    If you reach a node whose depth is already known, stop early.
    */
    for ( i = 0; i < nnode; ++i )
        depth_node[i] = 0;

    for ( i = 1; i <= nnode; ++i )
    {
        j = i;

        while ( j != root )
        {
            ++ depth_node[i-1];
            j = parent[j-1];

            if ( 0 < depth_node[j-1] )
            {
                depth_node[i-1] += depth_node[j-1];
                break;
            }
        }
    }
    /*
    Determine the maximum depth.
    */
    *depth = i4vec_max ( nnode, depth_node );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tree_rooted_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ROOTED_ENUM counts the number of unlabeled rooted trees.
  Example:
    Input    Output
      1         1
      2         1
      3         2
      4         4
      5         9
      6        20
      7        48
      8       115
      9       286
     10       719
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 July 2013
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int NNODE, the number of nodes.
    Output, int TREE_ROOTED_ENUM[NNODE].  NTREE(I) is the number of
    rooted, unlabeled trees on I nodes, for I = 1, 2, ... NNODE.
*/
{
	const register dim_typ nnode = *(dim_typ *) data;
	
    dim_typ i;
    dim_typ id;
    dim_typ isum;
    dim_typ itd;
    dim_typ j;
    dim_typ nlast;
    int *ntree;

    ntree = ( int * ) malloc ( nnode * sizeof ( int ) );

    ntree[0] = 1;

    for ( nlast = 2; nlast <= nnode; ++nlast )
    {
        isum = 0;

        for ( id = 1; id <= nlast - 1; ++id )
        {
            i = nlast;
            itd = ntree[id-1] * id;

            for ( j = 1; j <= nlast - 1; ++j )
            {
                i -= id;

                if ( i <= 0 )
                    break;
                isum += ntree[i-1] * itd;
            }
        }
        ntree[nlast-1] = isum / ( nlast - 1 );
    }
    return ntree;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tree_rooted_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    TREE_ROOTED_RANDOM selects a random unlabeled rooted tree.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 July 2013
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int NNODE, the number of nodes.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, int TREE_ROOTED_RANDOM[NNODE]. (I,ITREE(I)) is the I-th edge
    of the output tree for I = 2,NNODE.  ITREE(1)=0.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ nnode = s_data->a0;
	int * seed = s_data->a1;
	
    dim_typ i;
    dim_typ id;
    dim_typ is1;
    dim_typ is2;
    dim_typ itd;
    int *itree;
    dim_typ iz;
    dim_typ j;
    dim_typ l;
    dim_typ ll;
    dim_typ ls;
    dim_typ m;
    dim_typ *ntree;
    dim_typ nval;
    ityp r;
    int *stack;

    if ( nnode <= 0  )
        return NULL;

    itree = ( int * ) malloc ( nnode * sizeof ( int ) );
    stack = ( int * ) malloc ( nnode * sizeof ( int ) << 1 );
    /*
    Compute a table of the number of such trees for a given number of nodes.
    */
    ntree = tree_rooted_enum ( nnode );
    /*
    Now select one such tree at random.
    */

    nval = nnode;
    is1 = is2 = l = 0;

    for ( ; ; )
    {
        while ( 2 < nval )
        {
            r = r8_uniform_01 ( seed );

            iz = ( int ) ( ( nval - 1 ) * ntree[nval-1] * r );

            id = j = 0;

            ++ id;
            itd = id * ntree[id-1];
            m = nval;

            for ( ; ; )
            {
                ++ j;
                m -= id;

                if ( m < 1 )
                {
                    ++ id;
                    itd = id * ntree[id-1];
                    m = nval;
                    j = 0;
                    continue;
                }

                iz -= ntree[m-1] * itd;
                if ( iz < 0 )
                    break;
            }
            is1 = is1 + 1;
            stack[0+((is1-1)<<1)] = j;
            stack[1+((is1-1)<<1)] = id;
            nval = m;
        }

        itree[is2] = l;
        l = is2 + 1;
        is2 += nval;

        if ( 1 < nval )
            itree[is2-1] = is2 - 1;

        for ( ; ; )
        {
            nval = stack[1+((is1-1)<<1)];

            if ( nval != 0 )
            {
                stack[1+((is1-1)<<1)] = 0;
                break;
            }

            j = stack[0+((is1-1)<<1)];
            is1 = is1 - 1;
            m = is2 - l + 1;
            ll = itree[l-1];
            ls = l + ( j - 1 ) * m - 1;

            if ( j != 1 )
                for ( i = l; i <= ls; ++i )
                {
                    itree[i+m-1] = itree[i-1] + m;
                    if ( ( ( i - l ) % m ) == 0 )
                        itree[i+m-1] = ll;
                }

            is2 = ls + m;

            if ( is2 == nnode )
            {
                free ( ntree );
                free ( stack );
                return itree;
            }
            l = ll;
        }
    }

    free ( ntree );
    free ( stack );

    return itree;
}

#endif
