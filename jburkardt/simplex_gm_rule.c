#ifndef __DISABLEDEEP_SIMPLEXGMRULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _simplex_unit_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_UNIT_VOLUME computes the volume of the unit simplex.
  Discussion:
    The formula is simple: volume = 1/DIM_NUM!.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
   04 September 2003
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the dimension of the space.
    Output, double SIMPLEX_UNIT_VOLUME, the volume of the cone.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ dim_num = *(dim_typ *) data;
	
	ityp volume = 1.00;
	for (dim_typ i = 1; i <= dim_num; ++i )
		volume /= ( ( ityp ) i );
		
	result = volume;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _monomial_value ( void * data)
/******************************************************************************/
/*
  Purpose:
    MONOMIAL_VALUE evaluates a monomial.
  Discussion:
    This routine evaluates a monomial of the form
      product ( 1 <= dim <= m ) x(dim)^expon(dim)
    where the exponents are nonnegative integers.  Note that
    if the combination 0^0 is encountered, it should be treated
    as 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of points at which the
    monomial is to be evaluated.
    Input, double X[M*N], the point coordinates.
    Input, int EXPON[M], the exponents.
    Output, double MONOMIAL_VALUE[N], the value of the monomial.
*/
{
	const _2dtpipit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * expon = s_data->a2;
	ityp * x = s_data->a3;
	
	dim_typ i, j;
	ityp *value = ( ityp * ) malloc ( n * sizeof ( ityp ) );
	
	for ( j = 0; j < n; ++j )
		value[j] = 1.00;
	
	for ( i = 0; i < m; ++i )
		if ( 0 != expon[i] )
			for ( j = 0; j < n; ++j )
				value[j] *= pow ( x[i+j*m], expon[i] );
	
	return value;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gm_rule_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    GM_RULE_SET sets a Grundmann-Moeller rule.
  Discussion:
    This is a revised version of the calculation which seeks to compute
    the value of the weight in a cautious way that avoids intermediate
    overflow.  Thanks to John Peterson for pointing out the problem on
    26 June 2008.
    This rule returns weights and abscissas of a Grundmann-Moeller
    quadrature rule for the DIM_NUM-dimensional unit simplex.
    The dimension POINT_NUM can be determined by calling GM_RULE_SIZE.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2012
  Author:
    John Burkardt
  Reference:
    Axel Grundmann, Michael Moeller,
    Invariant Integration Formulas for the N-Simplex
    by Combinatorial Methods,
    SIAM Journal on Numerical Analysis,
    Volume 15, Number 2, April 1978, pages 282-290.
  Parameters:
    Input, int RULE, the index of the rule.
    0 <= RULE.
    Input, int DIM_NUM, the spatial dimension.
    1 <= DIM_NUM.
    Input, int POINT_NUM, the number of points in the rule.
    Output, double W[POINT_NUM], the weights.
    Output, double X[DIM_NUM*POINT_NUM], the abscissas.
*/
{
	const i2dt2pit * const s_data = data;
	int rule = s_data->a0;
	const register dim_typ dim_num = s_data->a1;
	const register dim_typ point_num = s_data->a2;
	ityp * w = s_data->a3;
	ityp * x = s_data->a4;
	
    int *beta;
    int beta_sum;
    int d;
    int dim;
    int h;
    dim_typ i;
    dim_typ j;
    int j_hi;
    dim_typ k;
    int more;
    int n;
    ityp one_pm;
    int s;
    int t;
    ityp weight;

    s = rule;
    d = (s<<1) + 1;
    k = 0;
    n = dim_num;
    one_pm = 1.00;

    beta = ( int * ) malloc ( ( dim_num + 1 ) * sizeof ( int ) );

    for ( i = 0; i <= s; ++i )
    {
        weight = ( ityp ) one_pm;

        j_hi = MAX ( n, MAX ( d, d + n - i ) );

        for ( j = 1; j <= j_hi; ++j )
        {
            if ( j <= n )
                weight *= ( ityp ) ( j );
            if ( j <= d )
                weight *= ( ityp ) ( d + n - (i<<1) );
            if ( j <= (s<<1) )
                weight /= 2.00;
            if ( j <= i )
                weight /= ( ityp ) ( j );
            if ( j <= d + n - i )
                weight /= ( ityp ) ( j );
        }

        one_pm = - one_pm;

        beta_sum = s - i;
        more = h = t = 0;

        for ( ; ; )
        {
            comp_next ( beta_sum, dim_num + 1, beta, &more, &h, &t );

            w[k] = weight;
            for ( dim = 0; dim < dim_num; ++dim )
                x[dim+k*dim_num] = ( ityp ) ( (beta[dim+1]<<1) + 1 )/ ( ityp ) ( d + n - (i<<1) );
            ++ k;

            if ( !more )
                break;
        }
    }

    free ( beta );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gm_rule_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    GM_RULE_SIZE determines the size of a Grundmann-Moeller rule.
  Discussion:
    This rule returns the value of POINT_NUM, the number of points associated
    with a GM rule of given index.
    After calling this rule, the user can use the value of POINT_NUM to
    allocate space for the weight vector as W(POINT_NUM) and the abscissa
    vector as X(DIM_NUM,POINT_NUM), and then call GM_RULE_SET.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2012
  Author:
    John Burkardt
  Reference:
    Axel Grundmann, Michael Moeller,
    Invariant Integration Formulas for the N-Simplex
    by Combinatorial Methods,
    SIAM Journal on Numerical Analysis,
    Volume 15, Number 2, April 1978, pages 282-290.
  Parameters:
    Input, int RULE, the index of the rule.
    0 <= RULE.
    Input, int DIM_NUM, the spatial dimension.
    1 <= DIM_NUM.
    Output, int GM_RULE_SIZE, the number of points in the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const idt * const s_data = data;
	int rule = s_data->a0;
	const register dim_typ dim_num = s_data->a1;
	
	result = i4_choose ( dim_num + rule + 1, rule );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _simplex_unit_monomial_int ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_UNIT_MONOMIAL_INT integrates a monomial over a simplex.
  Discussion:
    This routine evaluates a monomial of the form
      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
    where the exponents are nonnegative integers.  Note that
    if the combination 0^0 is encountered, it should be treated
    as 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2012
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int EXPON[DIM_NUM], the exponents.
    Output, double SIMPLEX_UNIT_MONOMIAL_INT, the value of the integral
    of the monomial.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	int * expon = s_data->a1;
	
    dim_typ dim, i, k;
    ityp value;

    value = 1.00;
    k = 0;

    for ( dim = 0; dim < dim_num; ++dim )
        for ( i = 1; i <= expon[dim]; ++i )
        {
            ++ k;
            value *= ( ityp ) ( i ) / ( ityp ) ( k );
        }

    for ( dim = 0; dim < dim_num; ++dim)
    {
        ++ k;
        value /=( ityp ) ( k );
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _simplex_unit_monomial_quadrature ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_UNIT_MONOMIAL_QUADRATURE: quadrature of monomials in a unit simplex.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2012
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int EXPON[DIM_NUM], the exponents.
    Input, int POINT_NUM, the number of points in the rule.
    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
    Input, double W[POINT_NUM], the quadrature weights.
    Output, double SIMPLEX_UNIT_MONOMIAL_QUADRATURE, the quadrature error.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitdtpipit * const s_data = data;
	
	const register dim_typ dim_num = s_data->a0;
	ityp * x = s_data->a1;
	const register dim_typ point_num = s_data->a2;
	int * expon = s_data->a3;
	ityp * w = s_data->a4;
	
    ityp exact = 1.00;
    ityp quad;
    ityp scale;
    ityp *value;
    ityp volume;
    /*
    Get the exact value of the integral of the unscaled monomial.
    */
    scale = simplex_unit_monomial_int ( dim_num, expon );
    /*
    Evaluate the monomial at the quadrature points.
    */
    value = monomial_value ( dim_num, point_num, x, expon );
    /*
    Compute the weighted sum and divide by the exact value.
    */
    volume = simplex_unit_volume ( dim_num );
    quad = volume * r8vec_dot_product ( point_num, w, value ) / scale;

    free ( value );
    /*
    Error:
    */
    
    result = abs ( quad - exact );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _simplex_unit_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_UNIT_SAMPLE returns uniformly random points from a general simplex.
  Discussion:
    The interior of the unit DIM_NUM dimensional simplex is the set of
    points X(1:DIM_NUM) such that each X(I) is nonnegative, and
    sum(X(1:DIM_NUM)) <= 1.
    This routine is valid for any spatial dimension DIM_NUM.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2012
  Author:
    John Burkardt
  Reference:
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.
  Parameters:
    Input, int DIM_NUM, the dimension of the space.
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double UNIFORM_IN_SIMPLEX01_MAP[DIM_NUM*N], the points.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * seed = s_data->a2;
	
    ityp *e;
    dim_typ i, j;
    ityp total;
    ityp *x;
    /*
    The construction begins by sampling DIM_NUM+1 points from the
    exponential distribution with parameter 1.
    */
    x = ( ityp * ) malloc ( dim_num * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
    {
        e = r8vec_uniform_01_new ( dim_num+1, seed );
        total = 0.00;
        for ( i = 0; i <= dim_num; ++i)
        {
            e[i] = -log ( e[i] );
            total += e[i];
        }

        for ( i = 0; i < dim_num; ++i)
            x[i+dim_num*j] = e[i] / total;
        free ( e );
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _simplex_unit_to_general ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_UNIT_TO_GENERAL maps the unit simplex to a general simplex.
  Discussion:
    Given that the unit simplex has been mapped to a general simplex
    with vertices T, compute the images in T, under the same linear
    mapping, of points whose coordinates in the unit simplex are REF.
    The vertices of the unit simplex are listed as suggested in the
    following:
   (0,0,0,...,0)
   (1,0,0,...,0)
   (0,1,0,...,0)
   (0,0,1,...,0)
  (...........)
   (0,0,0,...,1)
    Thanks to Andrei ("spiritualworlds") for pointing out a mistake in the
    previous implementation of this routine, 02 March 2008.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2012
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int POINT_NUM, the number of points to transform.
    Input, double T[DIM_NUM*(DIM_NUM+1)], the vertices of the
    general simplex.
    Input, double REF[DIM_NUM*POINT_NUM], points in the
    reference triangle.
    Output, double SIMPLEX_UNIT_TO_GENERAL[DIM_NUM*POINT_NUM],
    corresponding points in the physical triangle.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ point_num = s_data->a1;
	ityp * t = s_data->a2;
	ityp * ref = s_data->a3;
	
    dim_typ dim;
    ityp *phy;
    dim_typ point;
    dim_typ vertex;

    phy = ( ityp * ) malloc ( dim_num * point_num * sizeof ( ityp ) );
    //
    //  The image of each point is initially the image of the origin.
    //
    //  Insofar as the pre-image differs from the origin in a given vertex
    //  direction, add that proportion of the difference between the images
    //  of the origin and the vertex.
    //
    for ( point = 0; point < point_num; ++point )
        for ( dim = 0; dim < dim_num; ++dim )
        {
            phy[dim+point*dim_num] = t[dim+0*dim_num];

            for ( vertex = 1; vertex < dim_num + 1; ++vertex )
                phy[dim+point*dim_num] = phy[dim+point*dim_num]+ ( t[dim+vertex*dim_num] - t[dim+0*dim_num] ) * ref[vertex-1+point*dim_num];
        }

    return phy;
}

#endif
