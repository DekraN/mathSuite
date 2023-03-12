#ifndef __DISABLEDEEP_CUBEFELIPPARULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cube_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    CUBE_VOLUME: volume of a cube in 3D.
  Discussion:
    The integration region is:
      A(1) <= X <= B(1)
      A(2) <= Y <= B(2)
      A(3) <= Z <= B(3)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, double A[3], B[3], the lower and upper limits.
    Output, double CUBE_VOLUME, the volume.
*/
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	
    dim_typ i;
    ityp value = 1.00;
    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i)
        value *= ( b[i] - a[i] );
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_PRODUCT multiplies the entries of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    Input:
      A = ( 1, 2, 3, 4 )
    Output:
      I4VEC_PRODUCT = 24
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector
    Output, int I4VEC_PRODUCT, the product of the entries of A.
*/
{
	static int result = INT_MAX;
	
	const ipi * const s_data = data;
	const register int n = s_data->a0;
	int * a = s_data->a1;
	
    int product = 1;
    for (int i = 0; i < n; ++i)
        product *= a[i];

	result = product;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _line_unit_o01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_UNIT_O01 returns a 1 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[1], the weights.
    Output, double X[1], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[1] =
    {
        2.00
    };
    ityp x_save[1] =
    {
        0.00
    };

    r8vec_copy ( 1, w_save, w );
    r8vec_copy ( 1, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_unit_o02 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_UNIT_O02 returns a 2 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[2], the weights.
    Output, double X[2], the abscissas.
*/
{;
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];

    ityp w_save[2] =
    {
        1.0000000000000000000,
        1.0000000000000000000
    };
    ityp x_save[2] =
    {
        -0.57735026918962576451,
        0.57735026918962576451
    };

    r8vec_copy ( 2, w_save, w );
    r8vec_copy ( 2, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_unit_o03 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_UNIT_O03 returns a 3 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[3], the weights.
    Output, double X[3], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[3] =
    {
        0.55555555555555555556,
        0.88888888888888888889,
        0.55555555555555555556
    };
    ityp x_save[3] =
    {
        -0.77459666924148337704,
        0.00000000000000000000,
        0.77459666924148337704
    };

    r8vec_copy ( 3, w_save, w );
    r8vec_copy ( 3, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_unit_o04 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_UNIT_O04 returns a 4 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[4], the weights.
    Output, double X[4], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[4] =
    {
        0.34785484513745385737,
        0.65214515486254614263,
        0.65214515486254614263,
        0.34785484513745385737
    };
    ityp x_save[4] =
    {
        -0.86113631159405257522,
        -0.33998104358485626480,
        0.33998104358485626480,
        0.86113631159405257522
    };

    r8vec_copy ( 4, w_save, w );
    r8vec_copy ( 4, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _line_unit_o05 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_UNIT_O05 returns a 5 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[5], the weights.
    Output, double X[5], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[5] =
    {
        0.23692688505618908751,
        0.47862867049936646804,
        0.56888888888888888889,
        0.47862867049936646804,
        0.23692688505618908751
    };
    ityp x_save[5] =
    {
        -0.90617984593866399280,
        -0.53846931010568309104,
        0.00000000000000000000,
        0.53846931010568309104,
        0.90617984593866399280
    };

    r8vec_copy ( 5, w_save, w );
    r8vec_copy ( 5, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_direct_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
  Discussion:
    An R8VEC is a vector of R8's.
    To explain what is going on here, suppose we had to construct
    a multidimensional quadrature rule as the product of K rules
    for 1D quadrature.
    The product rule will be represented as a list of points and weights.
    The J-th item in the product rule will be associated with
      item J1 of 1D rule 1,
      item J2 of 1D rule 2,
      ...,
      item JK of 1D rule K.
    In particular,
      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
    and
      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
    So we can construct the quadrature rule if we can properly
    distribute the information in the 1D quadrature rules.
    This routine carries out that task.
    Another way to do this would be to compute, one by one, the
    set of all possible indices (J1,J2,...,JK), and then index
    the appropriate information.  An advantage of the method shown
    here is that you can process the K-th set of information and
    then discard it.
  Example:
    Rule 1:
      Order = 4
      X(1:4) = ( 1, 2, 3, 4 )
    Rule 2:
      Order = 3
      X(1:3) = ( 10, 20, 30 )
    Rule 3:
      Order = 2
      X(1:2) = ( 100, 200 )
    Product Rule:
      Order = 24
      X(1:24) =
  ( 1, 10, 100 )
  ( 2, 10, 100 )
  ( 3, 10, 100 )
  ( 4, 10, 100 )
  ( 1, 20, 100 )
  ( 2, 20, 100 )
  ( 3, 20, 100 )
  ( 4, 20, 100 )
  ( 1, 30, 100 )
  ( 2, 30, 100 )
  ( 3, 30, 100 )
  ( 4, 30, 100 )
  ( 1, 10, 200 )
  ( 2, 10, 200 )
  ( 3, 10, 200 )
  ( 4, 10, 200 )
  ( 1, 20, 200 )
  ( 2, 20, 200 )
  ( 3, 20, 200 )
  ( 4, 20, 200 )
  ( 1, 30, 200 )
  ( 2, 30, 200 )
  ( 3, 30, 200 )
  ( 4, 30, 200 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR_INDEX, the index of the factor being processed.
    The first factor processed must be factor 0.
    Input, int FACTOR_ORDER, the order of the factor.
    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values
    for factor FACTOR_INDEX.
    Input, int FACTOR_NUM, the number of factors.
    Input, int POINT_NUM, the number of elements in the direct product.
    Input/output, double X[FACTOR_NUM*POINT_NUM], the elements of the
    direct product, which are built up gradually.
  Local Parameters:
    Local, int START, the first location of a block of values to set.
    Local, int CONTIG, the number of consecutive values to set.
    Local, int SKIP, the distance from the current value of START
    to the next location of a block of values to set.
    Local, int REP, the number of blocks of values to set.
*/
{
	const _4dt2pit * const s_data = data;
	
	dim_typ factor_index = s_data->a0;
	dim_typ factor_order = s_data->a1;
	dim_typ factor_num = s_data->a2;
	dim_typ point_num = s_data->a3;
	ityp * factor_value = s_data->a4;
	ityp * x = s_data->a5;
	
    static dim_typ contig = 0;
    dim_typ i, j, k, start;
    static dim_typ rep = 0;
    static dim_typ skip = 0;

    if ( factor_index == 0 )
    {
        contig = 1;
        skip = 1;
        rep = point_num;
        for ( j = 0; j < point_num; ++j )
            for ( i = 0; i < factor_num; ++i )
                x[i+j*factor_num] = 0.00;
    }

    rep /= factor_order;
    skip = skip * factor_order;

    for ( i = 0; i < factor_order; ++i )
    {
        start = 0 + i * contig;
        for ( k = 1; k <= rep; ++k)
        {
            for ( j = start; j < start + contig; ++j )
                x[factor_index+j*factor_num] = factor_value[i];
            start += skip;
        }
    }
    contig *= factor_order;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_direct_product2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DIRECT_PRODUCT2 creates a direct product of r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    To explain what is going on here, suppose we had to construct
    a multidimensional quadrature rule as the product of K rules
    for 1D quadrature.
    The product rule will be represented as a list of points and weights.
    The J-th item in the product rule will be associated with
      item J1 of 1D rule 1,
      item J2 of 1D rule 2,
      ...,
      item JK of 1D rule K.
    In particular,
      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
    and
      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
    So we can construct the quadrature rule if we can properly
    distribute the information in the 1D quadrature rules.
    This routine carries out that task for the weights W.
    Another way to do this would be to compute, one by one, the
    set of all possible indices (J1,J2,...,JK), and then index
    the appropriate information.  An advantage of the method shown
    here is that you can process the K-th set of information and
    then discard it.
  Example:
    Rule 1:
      Order = 4
      W(1:4) = ( 2, 3, 5, 7 )
    Rule 2:
      Order = 3
      W(1:3) = ( 11, 13, 17 )
    Rule 3:
      Order = 2
      W(1:2) = ( 19, 23 )
    Product Rule:
      Order = 24
      W(1:24) =
  ( 2 * 11 * 19 )
  ( 3 * 11 * 19 )
  ( 4 * 11 * 19 )
  ( 7 * 11 * 19 )
  ( 2 * 13 * 19 )
  ( 3 * 13 * 19 )
  ( 5 * 13 * 19 )
  ( 7 * 13 * 19 )
  ( 2 * 17 * 19 )
  ( 3 * 17 * 19 )
  ( 5 * 17 * 19 )
  ( 7 * 17 * 19 )
  ( 2 * 11 * 23 )
  ( 3 * 11 * 23 )
  ( 5 * 11 * 23 )
  ( 7 * 11 * 23 )
  ( 2 * 13 * 23 )
  ( 3 * 13 * 23 )
  ( 5 * 13 * 23 )
  ( 7 * 13 * 23 )
  ( 2 * 17 * 23 )
  ( 3 * 17 * 23 )
  ( 5 * 17 * 23 )
  ( 7 * 17 * 23 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR_INDEX, the index of the factor being processed.
    The first factor processed must be factor 0.
    Input, int FACTOR_ORDER, the order of the factor.
    Input, ityp FACTOR_VALUE[FACTOR_ORDER], the factor values for
    factor FACTOR_INDEX.
    Input, int FACTOR_NUM, the number of factors.
    Input, int POINT_NUM, the number of elements in the direct product.
    Input/output, ityp W[POINT_NUM], the elements of the
    direct product, which are built up gradually.
  Local Parameters:
    Local, integer START, the first location of a block of values to set.
    Local, integer CONTIG, the number of consecutive values to set.
    Local, integer SKIP, the distance from the current value of START
    to the next location of a block of values to set.
    Local, integer REP, the number of blocks of values to set.
*/
{
	const _4dt2pit * const s_data = data;
	
	const register dim_typ factor_index = s_data->a0;
	const register dim_typ factor_order = s_data->a1;
	dim_typ factor_num = s_data->a2;
	dim_typ point_num = s_data->a3;
	ityp * factor_value = s_data->a4;
	ityp * w = s_data->a5;
	
    static int contig = 0;
    dim_typ i, j, k;
    static int rep = 0;
    static int skip = 0;
    int start;

    if ( factor_index == 0 )
    {
        contig = skip = 1;
        rep = point_num;
        for ( i = 0; i < point_num; ++i )
            w[i] = 1.00;
    }

    rep /= factor_order;
    skip *= factor_order;

    for ( j = 0; j < factor_order; ++j )
    {
        start = j * contig;

        for ( k = 1; k <= rep; ++k )
        {
            for ( i = start; i < start + contig; ++i )
                w[i] *= factor_value[j];
            start += skip;
        }
    }

    contig *= factor_order;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subcomp_next ( void * data )
/******************************************************************************/
/*
  Purpose:
    SUBCOMP_NEXT computes the next subcomposition of N into K parts.
  Discussion:
    A composition of the integer N into K parts is an ordered sequence
    of K nonnegative integers which sum to a value of N.
    A subcomposition of the integer N into K parts is a composition
    of M into K parts, where 0 <= M <= N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2009
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the integer whose subcompositions are desired.
    Input, int K, the number of parts in the subcomposition.
    Input/output, int A[K], the parts of the subcomposition.
    Input/output, int *MORE, set by the user to start the computation,
    and by the routine to terminate it.
    Input/output, int *H, *T, two internal parameters needed for the
    computation.  The user should allocate space for these in the calling
    program, include them in the calling sequence, but never alter them!
*/
{
	const _2dt3pipb * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1; 
	int * a = s_data->a2;
	int * h = s_data->a3;
	int * t = s_data->a4;
	bool * more = s_data->a5;
	
    dim_typ i;
    static int more2 = 0;
    static int n2 = 0;
    /*
    The first computation.
    */
    if ( !( *more ) )
    {
        n2 = 0;

        for ( i = 0; i < k; ++i )
           a[i] = 0;
        more2 = *h = *t = 0;
        *more = true;
    }
    /*
    Do the next element at the current value of N.
    */
    else if ( more2 )
        comp_next ( n2, k, a, &more2, h, t );
    else
    {
        more2 = 00;
        ++ n2;

        comp_next ( n2, k, a, &more2, h, t );
    }
    /*
    Termination occurs if MORE2 = FALSE and N2 = N.
    */
    if ( !more2 && n2 == n )
        *more = false;

    return NULL;
}

#endif
