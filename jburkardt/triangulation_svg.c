#ifndef __DISABLEDEEP_TRIANGULATIONSVG

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _mesh_base_zero ( void * data)
/******************************************************************************/
/*
  Purpose:
    MESH_BASE_ZERO ensures that the element definition is zero-based.
  Discussion:
    The ELEMENT_NODE array contains nodes indices that form elements.
    The convention for node indexing might start at 0 or at 1.
    Since a C program will naturally assume a 0-based indexing, it is
    necessary to check a given element definition and, if it is actually
    1-based, to convert it.
    This function attempts to detect 1-based node indexing and correct it.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 October 2014
  Author:
    John Burkardt
  Parameters:
    Input, int NODE_NUM, the number of nodes.
    Input, int ELEMENT_ORDER, the order of the elements.
    Input, int ELEMENT_NUM, the number of elements.
    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
    definitions.
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ node_num = s_data->a0;
	const register dim_typ element_order = s_data->a1;
	const register dim_typ element_num = s_data->a2;
	int * element_node = s_data->a3;
	
    dim_typ element;
    int node;
    int node_max;
    int node_min;
    dim_typ order;

    node_min = + i4_huge;
    node_max = - i4_huge;
    for ( element = 0; element < element_num; ++element )
        for ( order = 0; order < element_order; ++order )
        {
            node = element_node[order+element*element_order];
            if ( node < node_min )
                node_min = node;
            if ( node_max < node )
                node_max = node;
        }

    if ( node_min == 1 && node_max == node_num )
        for ( element = 0; element < element_num; ++element )
            for ( order = 0; order < element_order; ++order)
            -- element_node[order+element*element_order];

    return NULL;
}

#endif
