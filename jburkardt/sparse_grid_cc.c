#ifndef __DISABLEDEEP_SPARSEGRIDCC

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_mop ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_MOP returns the I-th power of -1 as an I4 value.
  Discussion:
    An I4 is an int value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I, the power of -1.
    Output, int I4_MOP, the I-th power of -1.
*/
{
	static short result = SHRT_MAX;
	
	const register dim_typ i = *(dim_typ *) data;
	
	result = 1 - (((i%2)!= 0)<<1);
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cc_abscissa ( void * data)
/******************************************************************************/
/*
  Purpose:
    CC_ABSCISSA returns the I-th abscissa of the Clenshaw Curtis rule.
  Discussion:
    Our convention is that the abscissas are numbered from left to
    right.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int I, the index of the desired abscissa.  1 <= I <= ORDER.
    Output, double CC_ABSCISSA, the value of the I-th 
    abscissa in the rule of order ORDER.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ * const a_data = data;
	const register dim_typ order = a_data[0];
	const register dim_typ i = a_data[1];
	
	if ( order == 0 || i < 1 || order < i )
	{
		result = MAX_VAL;
		return &result;
	}
	
	if ( order == 1 )
	{
		result = 0.00;
		return &result;
	}
		
	result = (i<<1) - 1 == order ? 0.00 : cos ( ( ityp ) ( order - i ) * M_PI / ( ityp ) ( order - 1 ) );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cc_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    CC_WEIGHTS computes Clenshaw Curtis weights.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Charles Clenshaw, Alan Curtis,
    A Method for Numerical Integration on an Automatic Computer,
    Numerische Mathematik,
    Volume 2, Number 1, December 1960, pages 197-205.
  Parameters:
    Input, int N, the order of the rule.
    Output, double CC_WEIGHTS[N], the weights of the rule.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
	ityp b;
	dim_typ i, j;
	ityp theta;
	ityp *w = ( ityp * ) malloc ( n * sizeof ( ityp ) );
	
	if ( n == 1 )
	{
		w[0] = 2.0;
		return w;
	}
	
	for ( i = 1; i <= n; ++i )
	{
		theta = ( ityp ) ( i - 1 ) * M_PI / ( ityp ) ( n - 1 );
		
		w[i-1] = 1.00;
		
		for ( j = 1; j <= ( n - 1 ) / 2; ++j )
		{
			b = 1.00 + ((j<<1) != n-1); 
			w[i-1] -= b * cos ( 2.00 * ( ityp ) ( j ) * theta ) / ( ityp ) ( (j<<2) * j - 1 );
		}
	}
	
	w[0] /= ( ityp ) ( n - 1 );
	for ( i = 1; i < n-1; ++i )
		w[i] = 2.00 * w[i] / ( ityp ) ( n - 1 );
	w[n-1] /= ( ityp ) ( n - 1 );
	
	return w;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _abscissa_level_closed_nd ( void * data)
/******************************************************************************/
/*
  Purpose:
    ABSCISSA_LEVEL_CLOSED_ND: first level at which an abscissa is generated.
  Discussion:
    We need this routine because the sparse grid is generated as a sum of
    product grids, and many points in the sparse grid will belong to several
    of these product grids, and we need to do something special the very
    first time we encounter such a point - namely, count it.  So this routine
    determines, for any point in the full product grid, the first level
    at which that point would be included.
    We assume an underlying product grid.  In each dimension, this product
    grid has order 2^LEVEL_MAX + 1.
    We will say a sparse grid has total level LEVEL if each point in the
    grid has a total level of LEVEL or less.
    The "level" of a point is determined as the sum of the levels of the
    point in each spatial dimension.
    The level of a point in a single spatial dimension I is determined as
    the level, between 0 and LEVEL_MAX, at which the point's I'th index
    would have been generated.
    This description is terse and perhaps unenlightening.  Keep in mind
    that the product grid is the product of 1D grids,
    that the 1D grids are built up by levels, having
    orders (total number of points ) 1, 3, 5, 9, 17, 33 and so on,
    and that these 1D grids are nested, so that each point in a 1D grid
    has a first level at which it appears.
    Our procedure for generating the points of a sparse grid, then, is
    to choose a value LEVEL_MAX, to generate the full product grid,
    but then only to keep those points on the full product grid whose
    LEVEL is less than or equal to LEVEL_MAX.
    Note that this routine is really just testing out the idea of
    determining the level.  Our true desire is to be able to start
    with a value LEVEL, and determine, in a straightforward manner,
    all the points that are generated exactly at that level, or
    all the points that are generated up to and including that level.
    This allows us to generate the new points to be added to one sparse
    grid to get the next, or to generate a particular sparse grid at once.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int LEVEL_MAX, controls the size of the final sparse grid.
    Input, int DIM_NUM, the spatial dimension.
    Input, int TEST_NUM, the number of points to be tested.
    Input, int TEST_VAL[DIM_NUM*TEST_NUM], the indices of the points
    to be tested.  Normally, each index would be between 0 and 2^LEVEL_MAX.
    Output, int ABSCISSA_LEVEL_ND[TEST_NUM], the value of LEVEL at which the
    point would first be generated, assuming that a standard sequence of
    nested grids is used.
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ level_max = s_data->a0;
	const register dim_typ dim_num = s_data->a1;
	const register dim_typ test_num = s_data->a2;
	int * test_val = s_data->a3;
	
    dim_typ dim;
    dim_typ j;
    dim_typ order;
    dim_typ t;
    int *test_level;

    test_level = ( int * ) malloc ( test_num * sizeof ( int ) );

    if ( level_max == 0 )
    {
        for ( j = 0; j < test_num; ++j)
            test_level[j] = 0;
        return test_level;
    }

    order = powi ( 2, level_max ) + 1;

    for ( j = 0; j < test_num; ++j )
        test_level[j] = index_to_level_closed ( dim_num, test_val+j*dim_num, order, level_max );

    return test_level;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _index_to_level_closed ( void * data)
/******************************************************************************/
/*
  Purpose:
    INDEX_TO_LEVEL_CLOSED determines the level of a point given its index.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int T[DIM_NUM], the grid indices of a point in a 1D closed rule.
    0 <= T[I] <= ORDER.
    Input, int ORDER, the order of the rule.
    Input, int LEVEL_MAX, the level with respect to which the
    index applies.
    Output, int INDEX_TO_LEVEL_CLOSED, the first level on which
    the point associated with the given index will appear.
*/
{
	static int result = INT_MAX;
	
	const _3dtpi * const s_data = data;
	
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ order = s_data->a1;
	const register dim_typ level_max = s_data->a2;
	int * t = s_data->a3;
	
    dim_typ dim;
    dim_typ level;
    dim_typ s;
    dim_typ value = 0;

    for ( dim = 0; dim < dim_num; ++dim )
    {
        s = t[dim];
        s = i4_modp ( s, order );

        if ( s == 0 )
      		level = 0;
        else
        {
            level = level_max;

            while ( ( s % 2 ) == 0 )
            {
                s /= 2;
                --  level;
            }
        }
        value += (level = level == 0);
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _level_to_order_closed ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEVEL_TO_ORDER_CLOSED converts a level to an order for closed rules.
  Discussion:
    Sparse grids can naturally be nested.  A natural scheme is to use
    a series of one-dimensional rules arranged in a series of "levels"
    whose order roughly doubles with each step.
    The arrangement described here works naturally for the Clenshaw Curtis
    and Newton Cotes closed rules.
    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single
    point at the center, and for all values afterwards, we use the
    relationship
      ORDER = 2^LEVEL + 1
    The following table shows how the growth will occur:
    Level    Order
    0          1
    1          3 =  2 + 1
    2          5 =  4 + 1
    3          9 =  8 + 1
    4         17 = 16 + 1
    5         33 = 32 + 1
    For the Clenshaw Curtis and Newton Cotes Closed rules, the point growth
    is nested.  If we have ORDER points on a particular LEVEL, the next
    level includes all these old points, plus ORDER-1 new points, formed
    in the gaps between successive pairs of old points.
    Level    Order = New + Old
    0          1   =  1  +  0
    1          3   =  2  +  1
    2          5   =  2  +  3
    3          9   =  4  +  5
    4         17   =  8  +  9
    5         33   = 16  + 17
    In this routine, we assume that a vector of levels is given,
    and the corresponding orders are desired.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int LEVEL[DIM_NUM], the nesting level.
    Output, int ORDER[DIM_NUM], the order (number of points)
    of the rule.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	int * level = s_data->a1;
	int * order = s_data->a2;
	
    dim_typ dim;

    for ( dim = 0; dim < dim_num; ++dim)
    {
        if ( level[dim] < 0 )
            order[dim] = -1;
        else if ( level[dim] == 0 )
            order[dim] = 1;
        else
            order[dim] = powi ( 2, level[dim] ) + 1 ;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _multigrid_index0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    MULTIGRID_INDEX0 returns an indexed multidimensional grid.
  Discussion:
    For dimension DIM, the second index of INDX may vary from
    0 to ORDER_1D[DIM]-1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int ORDER_1D[DIM_NUM], the order of the
    rule in each dimension.
    Input, int ORDER_ND, the product of the entries of ORDER_1D.
    Output, int INDX[DIM_NUM*ORDER_ND], the indices of the points in
    the grid.  The second dimension of this array is equal to the
    product of the entries of ORDER_1D.
*/
{
	const _2dtpi * const s_data = data;
	
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ order_nd = s_data->a1;
	int * order_1d = s_data->a2;
	
    int *a;
    dim_typ dim;
    bool more;
    int p;
    int *indx;

    indx = ( int * ) malloc ( dim_num * order_nd * sizeof ( int ) );
    a = ( int * ) malloc ( dim_num * sizeof ( int ) );
    more = p = 0;

    for ( ; ; )
    {
        vec_colex_next2 ( dim_num, order_1d, a, &more );

        if ( !more )
            break;

        for ( dim = 0; dim < dim_num; ++dim)
            indx[dim+p*dim_num] = a[dim];
        ++ p;
    }

    free ( a );

    return indx;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _multigrid_scale_closed ( void * data)
/******************************************************************************/
/*
  Purpose:
    MULTIGRID_SCALE_CLOSED renumbers a grid as a subgrid on a higher level.
  Discussion:
    This routine takes a grid associated with a given value of
    LEVEL, and multiplies all the indices by a power of 2, so that
    the indices reflect the position of the same points, but in
    a grid of level LEVEL_MAX.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int ORDER_ND, the number of points in the grid.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Input, int LEVEL_1D[DIM_NUM], the level in each dimension.
    Input/output, int GRID_INDEX[DIM_NUM*POINT_NUM], the index
    values for each grid point.  On input, these indices are based in
    the level for which the grid was generated; on output, the
    indices are appropriate for the grid as a subgrid of a grid
    of level LEVEL_MAX.
*/
{
	const _3dt2pi * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ order_nd = s_data->a1;
	const register dim_typ level_max = s_data->a2;
	int * level_1d = s_data->a3;
	int * grid_index = s_data->a4;
	
    dim_typ dim;
    dim_typ factor;
    dim_typ order;
    dim_typ order_max;

    for ( dim = 0; dim < dim_num; ++dim )
        if ( level_1d[dim] == 0 )
        {
            order_max = 0 == level_max ? 1 : powi ( 2, level_max ) + 1;
            for ( order = 0; order < order_nd; order++ )
                grid_index[dim+order*dim_num] = ( order_max - 1 ) / 2;
        }
        else
        {
            factor = powi ( 2, level_max - level_1d[dim] );
            for ( order = 0; order < order_nd; ++order)
                grid_index[dim+order*dim_num] = grid_index[dim+order*dim_num] * factor;
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _product_weights_cc ( void * data)
/******************************************************************************/
/*
  Purpose:
    PRODUCT_WEIGHTS_CC computes weights for a Clenshaw Curtis product rule.
  Discussion:
    This routine computes the weights for a quadrature rule which is
    a product of closed rules of varying order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
    Input, int ORDER_ND, the order of the product rule.
    Output, double PRODUCT_WEIGHTS_CC[DIM_NUM*ORDER_ND],
    the product rule weights.
*/
{
	const _2dtpi * const s_data = data;
	
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ order_nd = s_data->a1;
	int * order_1d = s_data->a2;
	
    dim_typ dim;
    dim_typ order;
    ityp *w_1d;
    ityp *w_nd = ( ityp * ) malloc ( order_nd * sizeof ( ityp ) );

    for ( order = 0; order < order_nd; ++order )
        w_nd[order] = 1.00;

    for ( dim = 0; dim < dim_num; ++dim )
    {
        w_1d = cc_weights ( order_1d[dim] );
        r8vec_direct_product2 ( dim, order_1d[dim], w_1d, dim_num,
        order_nd, w_nd );
        free ( w_1d );
    }

    return w_nd;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sparse_grid_cc ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPARSE_GRID_CC computes a sparse grid of Clenshaw Curtis points.
  Discussion:
    This program computes a quadrature rule and writes it to a file.
    The quadrature rule is associated with a sparse grid derived from
    a Smolyak construction using a closed 1D quadrature rule.
    The user specifies:
    * the spatial dimension of the quadrature region,
    * the level that defines the Smolyak grid.
    * the closed 1D quadrature rule (Clenshaw-Curtis or Newton-Cotes Closed).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int LEVEL_MAX, controls the size of the final sparse grid.
    Input, int POINT_NUM, the number of points in the grid, as determined
    by SPARSE_GRID_CC_SIZE.
    Output, double GRID_WEIGHTS[POINT_NUM], the weights.
    Output, double GRID_POINTS[DIM_NUM*POINT_NUM], the points.
*/
{
	const _3dt2pit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ level_max = s_data->a1;
	const register dim_typ point_num = s_data->a2;
	ityp * grid_weight = s_data->a3;
	ityp * grid_point = s_data->a4;
	
    dim_typ dim;
    int *grid_index;
    dim_typ order_max;
    dim_typ point;
    /*
    Determine the index vector, relative to the full product grid,
    that identifies the points in the sparse grid.
    */
    grid_index = sparse_grid_cc_index ( dim_num, level_max, point_num );
    /*
    Compute the physical coordinates of the abscissas.
    */
    order_max = 0 == level_max ? 1 : powi ( 2, level_max ) + 1;

    for ( point = 0; point < point_num; ++point )
        for ( dim = 0; dim < dim_num; ++dim )
        {
            grid_point[dim+point*dim_num] =
            cc_abscissa ( order_max, grid_index[dim+point*dim_num] + 1 );
        }
    /*
    Gather the weights.
    */
    sparse_grid_cc_weights ( dim_num, level_max, point_num, grid_index, grid_weight );

    free ( grid_index );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sparse_grid_cc_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPARSE_GRID_CC_INDEX indexes the points forming a sparse grid.
  Discussion:
    The points forming the sparse grid are guaranteed to be a subset
    of a certain product grid.  The product grid is formed by DIM_NUM
    copies of a 1D rule of fixed order.  The orders of the 1D rule,
 (called ORDER_1D) and the order of the product grid, (called ORDER)
    are determined from the value LEVEL_MAX.
    Thus, any point in the product grid can be identified by its grid index,
    a set of DIM_NUM indices, each between 1 and ORDER_1D.
    This routine creates the GRID_INDEX array, listing (uniquely) the
    points of the sparse grid.
    An assumption has been made that the 1D rule is closed (includes
    the interval endpoints) and nested (points that are part of a rule
    of a given level will be part of every rule of higher level).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Input, int POINT_NUM, the total number of points in the grids.
    Output, int SPARSE_GRID_CC_INDEX[DIM_NUM*POINT_NUM], a list of point
    indices, representing a subset of the product grid of level LEVEL_MAX,
    representing (exactly once) each point that will show up in a
    sparse grid of level LEVEL_MAX.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ dim_num = a_data[0];
	const register dim_typ level_max = a_data[1];
	const register dim_typ point_num = a_data[2];
	
    dim_typ dim;
    dim_typ factor;
    int *grid_index;
    int *grid_index2;
    int *grid_level;
    int h;
    dim_typ j;
    dim_typ level;
    int *level_1d;
    int more;
    int *order_1d;
    int order_nd;
    dim_typ point;
    dim_typ point_num2;
    int t;

    grid_index = ( int * ) malloc ( dim_num * point_num * sizeof ( int ) );
    /*
    The outer loop generates LEVELs from 0 to LEVEL_MAX.
    */
    point_num2 = 0;

    level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );

    for ( level = 0; level <= level_max; ++level)
    {
        /*
        The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
        that adds up to LEVEL.
        */
        more = h = t = 0;

        for ( ; ; )
        {
            comp_next ( level, dim_num, level_1d, &more, &h, &t );
            /*
            Transform each 1D level to a corresponding 1D order.
            */
            level_to_order_closed ( dim_num, level_1d, order_1d );
            /*
            The product of the 1D orders gives us the number of points in this grid.
            */
            order_nd = i4vec_product ( dim_num, order_1d );
            /*
            The inner (hidden) loop generates all points corresponding to given grid.
            */
            grid_index2 = multigrid_index0 ( dim_num, order_1d, order_nd );
            /*
            Adjust these grid indices to reflect LEVEL_MAX.
            */
            multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d,grid_index2 );
            /*
            Determine the first level of appearance of each of the points.
            */
            grid_level = abscissa_level_closed_nd ( level_max, dim_num, order_nd,grid_index2 );
            /*
            Only keep those points which first appear on this level.
            */
            for ( point = 0; point < order_nd; ++point )
                if ( grid_level[point] == level )
                {
                    for ( dim = 0; dim < dim_num; ++dim )
                        grid_index[dim+point_num2*dim_num] =grid_index2[dim+point*dim_num];
                    ++ point_num2;
                }

            free ( grid_index2 );
            free ( grid_level );

            if ( !more )
                break;
        }
    }
    free ( level_1d );
    free ( order_1d );

    return grid_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sparse_grid_cc_size_old ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPARSE_GRID_CC_SIZE_OLD sizes a sparse grid of Clenshaw Curtis points.
  Discussion:
    This function has been replaced by a new version which is much faster.
    This version is retained for historical interest.
    The grid is defined as the sum of the product rules whose LEVEL
    satisfies:
      0 <= LEVEL <= LEVEL_MAX.
    This routine works on an abstract set of nested grids.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Output, int SPARSE_GRID_CC_SIZE, the number of points in the grid.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ dim_num = a_data[0];
	const register dim_typ level_max = a_data[1];
	
    dim_typ dim;
    dim_typ factor;
    int *grid_index;
    int *grid_level;
    int h;
    dim_typ j;
    dim_typ level;
    int *level_1d;
    int more;
    int *order_1d;
    dim_typ order_max;
    dim_typ order_nd;
    dim_typ point;
    dim_typ point_num;
    int t;
    /*
    Special case.
    */
    if ( level_max == 0 )
    {
    	result = 1;
        return &result;
    }
    /*
    The outer loop generates LEVELs from 0 to LEVEL_MAX.
    */
    point_num = 0;

    level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );

    for ( level = 0; level <= level_max; ++level)
    {
        /*
        The middle loop generates the next partition that adds up to LEVEL.
        */
        more = h = t = 0;

        for ( ; ; )
        {
            comp_next ( level, dim_num, level_1d, &more, &h, &t );
            /*
            Transform each 1D level to a corresponding 1D order.
            */
            level_to_order_closed ( dim_num, level_1d, order_1d );
            /*
            The product of the 1D orders gives us the number of points in this grid.
            */
            order_nd = i4vec_product ( dim_num, order_1d );
            /*
            The inner (hidden) loop generates all points corresponding to given grid.
            */
            grid_index = multigrid_index0 ( dim_num, order_1d, order_nd );
            /*
            Adjust these grid indices to reflect LEVEL_MAX.
            */
            multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d,grid_index );
            /*
            Determine the first level of appearance of each of the points.
            */
            grid_level = abscissa_level_closed_nd ( level_max, dim_num, order_nd,grid_index );
            /*
            Only keep those points which first appear on this level.
            */
            for ( point = 0; point < order_nd; ++point )
                if ( grid_level[point] == level )
                    ++ point_num;

            free ( grid_index );
            free ( grid_level );

            if ( !more )
                break;
        }
    }

    free ( level_1d );
    free ( order_1d );

	result = point_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sparse_grid_cc_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPARSE_GRID_CC_WEIGHTS gathers the weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Input, int POINT_NUM, the total number of points in the grids.
    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of point indices,
    representing a subset of the product grid of level LEVEL_MAX,
    representing (exactly once) each point that will show up in a
    sparse grid of level LEVEL_MAX.
    Output, double GRID_WEIGHT[POINT_NUM], the weights
    associated with the sparse grid points.
*/
{
	const _3dtpipit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ level_max = s_data->a1;
	const register dim_typ point_num = s_data->a2;
	int * grid_index = s_data->a3;
	ityp * grid_weight = s_data->a4;
	
    bool all_equal;
    dim_typ coeff;
    dim_typ dim;
    dim_typ found;
    int *grid_index2;
    ityp *grid_weight2;
    int h;
    dim_typ level;
    int *level_1d;
    dim_typ level_min;
    int more;
    dim_typ order_nd;
    int *order_1d;
    dim_typ point;
    dim_typ point2;
    int t;

    if ( level_max == 0 )
    {
        for ( point = 0; point < point_num; ++point )
            grid_weight[point] = pow ( 2.00, dim_num );
        return NULL;
    }

    level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );

    for ( point = 0; point < point_num; ++point )
        grid_weight[point] = 0.00;

    level_min = MAX ( 0, level_max + 1 - dim_num );

    for ( level = level_min; level <= level_max; ++level )
    {
        /*
        The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
        that adds up to LEVEL.
        */
        more = h = t = 0;

        for ( ; ; )
        {
            comp_next ( level, dim_num, level_1d, &more, &h, &t );
            /*
            Transform each 1D level to a corresponding 1D order.
            */
            level_to_order_closed ( dim_num, level_1d, order_1d );
            /*
            The product of the 1D orders gives us the number of points in this grid.
            */
            order_nd = i4vec_product ( dim_num, order_1d );
            /*
            Generate the indices of the points corresponding to the grid.
            */
            grid_index2 = multigrid_index0 ( dim_num, order_1d, order_nd );
            /*
            Compute the weights for this grid.
            */
            grid_weight2 = product_weights_cc ( dim_num, order_1d, order_nd );
            /*
            Adjust the grid indices to reflect LEVEL_MAX.
            */
            multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d,grid_index2 );
            /*
            Now determine the coefficient.
            */
            coeff = i4_mop ( level_max - level )* i4_choose ( dim_num - 1, level_max - level );

            for ( point2 = 0; point2 < order_nd; ++point2 )
            {
                found = false;

                for ( point = 0; point < point_num; ++point )
                {
                    all_equal = true;
                    for ( dim = 0; dim < dim_num; ++dim)
                        if ( grid_index2[dim+point2*dim_num] !=grid_index[dim+point*dim_num] )
                        {
                            all_equal = false;
                            break;
                        }
                    if ( all_equal )
                    {
                        grid_weight[point] = grid_weight[point]+ ( ityp ) ( coeff ) * grid_weight2[point2];
                        found = true;
                        break;
                    }
                }
                if ( !found )
                    return NULL;
            }

            free ( grid_index2 );
            free ( grid_weight2 );

            if ( !more )
                break;
        }
    }
    free ( level_1d );
    free ( order_1d );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sparse_grid_ccs_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPARSE_GRID_CCS_SIZE sizes a sparse grid using Clenshaw Curtis Slow rules.
  Discussion:
    The grid is defined as the sum of the product rules whose LEVEL
    satisfies:
      0 <= LEVEL <= LEVEL_MAX.
    This calculation is much faster than a previous method.  It simply
    computes the number of new points that are added at each level in the
    1D rule, and then counts the new points at a given DIM_NUM dimensional
    level vector as the product of the new points added in each dimension.
    This approach will work for nested families, and may be extensible
    to other families, and to mixed rules.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 December 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Output, int SPARSE_GRID_CC_SIZE, the number of points in the grid.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ dim_num = a_data[0];
	const register dim_typ level_max = a_data[1];
	
    dim_typ dim;
    int h;
    dim_typ l;
    dim_typ level;
    int *level_1d;
    int more;
    int *new_1d;
    dim_typ o;
    dim_typ p;
    dim_typ point_num;
    int t;
    dim_typ v;
    /*
    Special case.
    */

    if ( level_max == 0 )
    {
    	result = 1;
        return &result;
    }
    /*
    Construct the vector that counts the new points in the 1D rule.
    */
    new_1d = ( int * ) malloc ( ( level_max + 1 ) * sizeof ( int ) );

    new_1d[0] = 1;
    new_1d[1] = 2;

    p = o = 3;

    for ( l = 2; l <= level_max; ++l )
    {
        p = (l<<1) + 1;
        if ( o < p )
        {
            new_1d[l] = o - 1;
            o = (o<<1) - 1;
        }
        else
            new_1d[l] = 0;
    }
    /*
    Count the number of points by counting the number of new points
    associated with each level vector.
    */
    level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );

    point_num = 0;

    for ( level = 0; level <= level_max; ++level )
    {
        more = h = t = 0;

        for ( ; ;)
        {
            comp_next ( level, dim_num, level_1d, &more, &h, &t );

            v = 1;
            for ( dim = 0; dim < dim_num; ++dim)
                v *= new_1d[level_1d[dim]];

            point_num += v;

            if ( !more )
                break;
        }
    }
    free ( level_1d );
    free ( new_1d );

	result = point_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _vec_colex_next2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VEC_COLEX_NEXT2 generates vectors in colex order.
  Discussion:
    The vectors are produced in colexical order, starting with
 (0,        0,        ...,0),
 (1,        0,        ...,0),
     ...
 (BASE(1)-1,0,        ...,0)
 (0,        1,        ...,0)
 (1,        1,        ...,0)
    ...
 (BASE(1)-1,1,        ...,0)
 (0,        2,        ...,0)
 (1,        2,        ...,0)
    ...
 (BASE(1)-1,BASE(2)-1,...,BASE(DIM_NUM)-1).
  Examples:
    DIM_NUM = 2,
    BASE = { 3, 3 }
    0   0
    1   0
    2   0
    0   1
    1   1
    2   1
    0   2
    1   2
    2   2
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    12 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
    In dimension I, entries will range from 0 to BASE[I]-1.
    Output, int A[DIM_NUM], the next vector.
    Input/output, int *MORE.  Set this variable false before
    the first call.  On return, MORE is TRUE if another vector has
    been computed.  If MORE is returned FALSE, ignore the output 
    vector and stop calling the routine.
*/
{
	const dt2pipb * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	int * base = s_data->a1;
	int * a = s_data->a2;
	bool * more = s_data->a3;
	
	dim_typ i;
	
	if ( !( *more ) )
	{
		for ( i = 0; i < dim_num; ++i)
			a[i] = 0;
		*more = true; 
	}
	else
	{
		for ( i = 0; i < dim_num; ++i )
		{
			++ a[i];
		
			if ( a[i] < base[i] )
				return NULL;
			a[i] = 0;
		}
		*more = false;
	}
	
	return NULL;
}

#endif
