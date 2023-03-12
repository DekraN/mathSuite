#ifndef __DISABLEDEEP_FEM1DPACK

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bandwidth_mesh ( void * data)
/******************************************************************************/
/*
  Purpose:
    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.
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
    11 September 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ELEMENT_ORDER, the order of the elements.
    Input, int ELEMENT_NUM, the number of elements.
    Input,  ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
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
    int global_i;
    int global_j;
    dim_typ local_i;
    dim_typ local_j;

    *ml = *mu = 0;

    for ( element = 0; element < element_num; ++element )
        for ( local_i = 0; local_i < element_order; ++local_i )
        {
            global_i = element_node[local_i+element*element_order];
            for ( local_j = 0; local_j < element_order; ++local_j)
            {
                global_j = element_node[local_j+element*element_order];
                *mu = MAX ( *mu, global_j - global_i );
                *ml = MAX ( *ml, global_i - global_j );
            }
        }

    *m = *ml + 1 + *mu;

    return NULL;
}

#endif
