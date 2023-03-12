#ifndef __DISABLEDEEP_SGMGA

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_copy ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_COPY copies an R8VEC.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    03 July 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, double A1[N], the vector to be copied.
    Input, double A2[N], the copy of A1.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;
	
	for (dim_typ i = 0; i < n; ++i )
		a2[i] = a1[i];
		
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _binary_vector_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    BINARY_VECTOR_NEXT generates the next binary vector.
  Discussion:
    A binary vector is a vector whose entries are 0 or 1.
    The user inputs an initial zero vector to start.  The program returns
    the "next" vector.
    The vectors are produced in the order:
 ( 0, 0, 0, ..., 0 )
 ( 1, 0, 0, ..., 0 ) 
 ( 0, 1, 0, ..., 0 )
 ( 1, 1, 0, ..., 0 )
 ( 0, 0, 1, ..., 0 )
 ( 1, 0, 1, ..., 0 )
               ...
 ( 1, 1, 1, ..., 1)
    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
    we allow wrap around.
  Example:
    N = 3
    Input      Output
    -----      ------
    0 0 0  =>  1 0 0
    1 0 0  =>  0 1 0
    0 1 0  =>  1 1 0
    1 1 0  =>  0 0 1
    0 0 1  =>  1 0 1
    1 0 1  =>  0 1 1
    0 1 1  =>  1 1 1
	1 1 1  =>  0 0 0
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    04 September 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vectors.
    Input/output, int BVEC[N], on output, the successor 
    to the input vector.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * bvec = s_data->a1;
	
	for (dim_typ i = 0; i < n; ++i )
	{  
		if ( bvec[i] == 1 )
			bvec[i] = false;
		else 
		{
			bvec[i] = true;
			break;
		}
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _vec_colex_next3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    VEC_COLEX_NEXT3 generates vectors in colex order.
  Discussion:
    The vectors are produced in colexical order, starting with
 (1,        1,        ...,1),
 (2,        1,        ...,1),
     ...
 (BASE(1),  1,        ...,1)
 (1,        2,        ...,1)
 (2,        2,        ...,1)
    ...
 (BASE(1),  2,        ...,1)
 (1,        3,        ...,1)
 (2,        3,        ...,1)
    ...
 (BASE(1),  BASE(2), ...,BASE(DIM_NUM)).
  Example:
    DIM_NUM = 2,
    BASE = { 3, 3 }
    1   1
    2   1
    3   1
    1   2
    2   2
    3   2
    1   3
    2   3
    3   3
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    22 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
    In dimension I, entries will range from 1 to BASE[I].
    Output, int A[DIM_NUM], the next vector.
    Input/output, int *MORE.  Set this variable 0 before
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
		for ( i = 0; i < dim_num; ++i )
			a[i] = 1;
		*more = true;
	}
	else
	{
		for ( i = 0; i < dim_num; ++i )
		{
			++ a[i];
			
			if ( a[i] <= base[i] )
				return NULL;
			a[i] = 1;
		}
		*more = false; 
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _point_unique_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    POINT_UNIQUE_INDEX indexes unique points.
  Discussion:
    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.
    The goal of this routine is to determine a vector UNDX,
    which points to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.
    This is all done with index vectors, so that the elements of
    A are never moved.
    The first step of the algorithm requires the indexed sorting
    of A, which creates arrays INDX and XDNI. (If all the entries
    of A are unique, then these arrays are the same as UNDX and XDNU.)
    We then use INDX to examine the entries of A in sorted order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.
    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula
      XU(*) = A(UNDX(*)).
    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:
      X(I) = XU(XDNU(I)).
    We could then replace A by the combination of XU and XDNU.
    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.
    Here is an example of a vector A, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.
      I     A  Indx  Xdni       XU  Undx  Xdnu
    ----+-----+-----+-----+--------+-----+-----+
      0 | 11.     0     0 |    11.     0     0
      1 | 22.     2     4 |    22.     1     1
      2 | 11.     5     1 |    33.     3     0
      3 | 33.     8     7 |    55.     4     2
      4 | 55.     1     8 |                  3
      5 | 11.     6     2 |                  0
      6 | 22.     7     5 |                  1
      7 | 22.     3     6 |                  1
      8 | 11.     4     3 |                  0
    INDX(2) = 3 means that sorted item(2) is A(3).
    XDNI(2) = 5 means that A(2) is sorted item(5).
    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
    XU(XDNU(I))) = A(I).
    XU(I)        = A(UNDX(I)).
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    19 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of the data values.
    Input, int N, the number of data values,
    Input, double A[M*N], the data values.
    Input, int UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.
    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
    Output, int XDNU[N], the XDNU vector.
*/
{
	const _3dtpit2pi * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register dim_typ unique_num = s_data->a2;
	ityp * a = s_data->a3;
	int * undx = s_data->a4;
	int * xdnu = s_data->a5;
	
	int base = 0;
	ityp diff;
	dim_typ i;
	int *indx;
	dim_typ j, k;
	/*
	Implicitly sort the array.
	*/
	indx = r8col_sort_heap_index_a ( m, n, base, a );
	/*
	Walk through the implicitly sorted array X.
	*/
	i = j = 0;
	undx[j] = indx[i];
	xdnu[indx[i]] = j;
	
	for ( i = 1; i < n; ++i )
	{
		diff = 0.0;
		for ( k = 0; k < m; ++k )
			diff = MAX ( diff, abs ( a[k+indx[i]*m] - a[k+undx[j]*m] ) );
		if ( 0.00 < diff )
		{
			++ j;
			undx[j] = indx[i];
		}
		xdnu[indx[i]] = j;
	}
	free ( indx );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _clenshaw_curtis_compute_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLENSHAW_CURTIS_COMPUTE_WEIGHTS computes Clenshaw Curtis quadrature weights.
  Discussion:
    The user must preallocate space for the output array W.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Reference:
    Charles Clenshaw, Alan Curtis,
    A Method for Numerical Integration on an Automatic Computer,
    Numerische Mathematik,
    Volume 2, Number 1, December 1960, pages 197-205.
  Parameters:
    Input, int ORDER, the order of the rule.
    Output, double W[ORDER], the weights of the rule.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * w = s_data->a1;
	
	ityp b;
	dim_typ i, j;
	ityp theta;
	
	if ( order == 0)
		return NULL;
	else if ( order == 1 )
	{
		w[0] = 2.00;
		return NULL;
	}
	
	for ( i = 1; i <= order; ++i )
	{
		theta = ( ityp ) ( i - 1 ) * M_PI / ( ityp ) ( order - 1 );
		
		w[i-1] = 1.00;
		
		for ( j = 1; j <= ( order - 1 ) / 2; ++j )
		{
			b = 1.00+(j<<1) != order-1;
			w[i-1] -= b *  cos ( 2.00 * ( ityp ) ( j ) * theta ) / ( ityp ) ( (j<<2) * j - 1 );
		}
	}
	
	w[0] /= ( ityp ) ( order - 1 );
	for ( i = 1; i < order - 1; ++i )
		w[i] = 2.00 * w[i] / ( ityp ) ( order - 1 );
	w[order-1] /= ( ityp ) ( order - 1 );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fejer2_compute_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEJER2_COMPUTE_WEIGHTS computes Fejer type 2 quadrature weights.
  Discussion:
    The user must preallocate space for the output array W.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  Parameters:
    Input, int ORDER, the order.
    Output, double W[ORDER], the weights.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * w = s_data->a1;
	
	dim_typ i, j;
	ityp p;
	ityp theta;
	
	if ( order == 0)
		return NULL;
	else if ( order == 1 )
		w[0] = 2.00;
	else if ( order == 2 )
		w[0] = w[1] = 1.00;
	else
	{
		for ( i = 1; i <= order; ++i )
		{
			theta = ( ityp ) ( order + 1 - i ) * M_PI / ( ityp ) ( order + 1 );
			w[i-1] = 1.00;
		
			for ( j = 1; j <= ( ( order - 1 ) / 2 ); ++j )
				w[i-1] -= 2.00 *  cos ( 2.00 * ( ityp ) ( j ) * theta ) / ( ityp ) ( (j<<2) * j - 1 );
			p = 2.00 * ( ityp ) ( ( ( order + 1 ) / 2 ) ) - 1.00;
			w[i-1] -= cos ( ( p + 1.0 ) * theta ) / p;
		}
		for ( i = 0; i < order; ++i )
			w[i] *= 2.00 / ( ityp ) ( order + 1 );
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fejer2_compute_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEJER2_COMPUTE_WEIGHTS_NP computes Fejer type 2 quadrature weights.
  Discussion:
    The user must preallocate space for the output array W.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  Parameters:
    Input, int ORDER, the order.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double W[ORDER], the weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	fejer2_compute_weights ( order, w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gen_hermite_compute_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_HERMITE_COMPUTE_WEIGHTS computes Generalized Hermite quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, double ALPHA, the exponent of the X factor.
    -1.0 < ALPHA.
    Output, double W[ORDER], the weights.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ order = s_data->a0; 
	ityp * w = s_data->a1;
	const register ityp alpha = s_data->a2; 
	
	
	ityp *x = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	gen_hermite_compute ( order, alpha, x, w );
	free ( x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gen_hermite_compute_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_HERMITE_COMPUTE_WEIGHTS_NP: Generalized Hermite quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor. -1.0 < ALPHA.
    Output, double W[ORDER], the weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	gen_hermite_compute_weights ( order, p[0], w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gen_laguerre_compute_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_LAGUERRE_COMPUTE_WEIGHTS: Generalized Laguerre quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, double ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.
    Output, double W[ORDER], the weights.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ order = s_data->a0; 
	ityp * w = s_data->a1;
	const register ityp alpha = s_data->a2; 
	
	
	ityp *x = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	gen_laguerre_compute ( order, alpha, x, w );
	free ( x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gen_laguerre_compute_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_LAGUERRE_COMPUTE_WEIGHTS_NP: Generalized Laguerre quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.
    Output, double W[ORDER], the weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	gen_laguerre_compute_weights ( order, p[0], w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hermite_compute_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_COMPUTE_WEIGHTS computes Hermite quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Output, double W[ORDER], the weights.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * w = s_data->a1;
	
	ityp *x = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	hermite_compute ( order, x, w );
	free( x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hermite_compute_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_COMPUTE_WEIGHTS_NP computes Hermite quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double W[ORDER], the weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	hermite_compute_weights ( order, w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _jacobi_compute_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_COMPUTE_WEIGHTS computes Jacobi quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, double ALPHA, BETA, the exponents of the (1-X) and (1+X) factors.
    Output, double W[ORDER], the weights.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register ityp alpha = s_data->a1; 
	const register ityp beta = s_data->a2; 
	ityp * w = s_data->a3;
	
	ityp *x = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	jacobi_compute ( order, alpha, beta, x, w );
	free ( x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _jacobi_compute_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_COMPUTE_WEIGHTS_NP computes Jacobi quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameter values.
    P[0] = ALPHA, the exponent of (1-X)
    P[1] = BETA,  the exponent of (1+X).
    -1.0 < ALPHA and -1.0 < BETA are required.
    Output, double W[ORDER], the weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	jacobi_compute_weights ( order, p[0], p[1], w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laguerre_compute_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_COMPUTE_WEIGHTS computes Laguerre quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Output, double W[ORDER], the weights.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * w = s_data->a1;
	
	ityp *x = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	laguerre_compute ( order, x, w );
	free ( x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laguerre_compute_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_COMPUTE_WEIGHTS_NP computes Laguerre quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double W[ORDER], the weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	laguerre_compute_weights ( order, w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _legendre_compute_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_COMPUTE_WEIGHTS computes Legendre quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Output, double W[ORDER], the weights.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * w = s_data->a1;
	
	ityp *x = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	legendre_compute ( order, x, w );
	free ( x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _legendre_compute_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_COMPUTE_WEIGHTS_NP computes Legendre quadrature weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double W[ORDER], the weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	legendre_compute_weights ( order, w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _clenshaw_curtis_compute_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLENSHAW_CURTIS_COMPUTE_WEIGHTS_NP: Clenshaw Curtis quadrature weights.
  Discussion:
    The user must preallocate space for the output array W.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Reference:
    Charles Clenshaw, Alan Curtis,
    A Method for Numerical Integration on an Automatic Computer,
    Numerische Mathematik,
    Volume 2, Number 1, December 1960, pages 197-205.
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double W[ORDER], the weights of the rule.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	clenshaw_curtis_compute_weights ( order, w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fejer2_compute ( void * data)
/******************************************************************************/
/*
  Purpose:

    FEJER2_COMPUTE computes a Fejer type 2 rule.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    The rule is defined on [-1,1].
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the Fejer type 2 rule.
    1 <= ORDER.
    Output, double X[ORDER], the abscissas.
    Output, double W[ORDER], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
	dim_typ i, j;
	ityp p;
	ityp theta;
	
	if ( order == 0 || order == 1)
		return NULL;
	
	for ( i = 0; i < order; ++i )
		x[i] =  cos ( ( ityp ) ( order - i ) * M_PI / ( ityp ) ( order + 1 ) );
	if ( ( order % 2 ) == 1 )
		x[(order-1)/2] = 0.00;
	
	if ( order == 2 )
		w[0] = w[1] = 1.00;
	else
	{
		for ( i = 0; i < order; ++i )
		{
			theta = ( ityp ) ( order - i ) * M_PI / ( ityp ) ( order + 1 );
			
			w[i] = 1.00;
			
			for ( j = 1; j <= ( ( order - 1 ) / 2 ); ++j )
				w[i] -= 2.00 *  cos ( 2.00 * ( ityp ) ( j ) * theta ) / ( ityp ) ( (j<<2) * j - 1 );
			p = 2.00 * ( ityp ) ( ( ( order + 1 ) / 2 ) ) - 1.00;
			w[i] -= cos ( ( p + 1.00 ) * theta ) / p;
		}
		for ( i = 0; i < order; ++i)
			w[i] = 2.00 * w[i] / ( ityp ) ( order + 1 );
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fejer2_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEJER2_COMPUTE_POINTS computes Fejer type 2 quadrature points.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    The rule is defined on [-1,1].
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the Fejer type 2 rule.
    1 <= ORDER.
    Output, double X[ORDER], the abscissas.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	
	if ( order == 0)
		return NULL;
	else if ( order == 1 )
		x[0] = 0.00;
	else
	{
		for (dim_typ index = 1; index <= order; ++index )
			x[index-1] =  cos ( ( ityp ) ( order + 1 - index ) * M_PI / ( ityp ) ( order + 1 ) );
		if ( ( order % 2 ) == 1 )
			x[(order-1)/2] = 0.00;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _fejer2_compute_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEJER2_COMPUTE_POINTS_NP computes Fejer type 2 quadrature points.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    The rule is defined on [-1,1].
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the Fejer type 2 rule.
    1 <= ORDER.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	fejer2_compute_points ( order, x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _patterson_lookup_weights ( void * data)
/******************************************************************************/
/*
  Purpose:
    PATTERSON_LOOKUP_WEIGHTS looks up Patterson quadrature weights.
  Discussion:
    The allowed orders are 1, 3, 7, 15, 31, 63, 127 and 255.
    The weights are positive, symmetric and should sum to 2.
    The user must preallocate space for the output array W.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    17 December 2009
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  Parameters:
    Input, int ORDER, the order of the rule.
    ORDER must be 1, 3, 7, 15, 31, 63, 127 and 255.
    Output, double W[ORDER], the weights.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * w = s_data->a1;
	
	static ityp w_001[1] =
	{
		2.00
	};
	static ityp w_003[3] = 
	{
		0.555555555555555555556,
		0.888888888888888888889,
		0.555555555555555555556
	};
	static ityp w_007[7] =
	{
		0.104656226026467265194,
		0.268488089868333440729,
		0.401397414775962222905,
		0.450916538658474142345,
		0.401397414775962222905,
		0.268488089868333440729,
		0.104656226026467265194
	};
	static ityp w_015[15] =
	{
		0.0170017196299402603390,
		0.0516032829970797396969,
		0.0929271953151245376859,
		0.134415255243784220360,
		0.171511909136391380787,
		0.200628529376989021034,
		0.219156858401587496404,
		0.225510499798206687386,
		0.219156858401587496404,
		0.200628529376989021034,
		0.171511909136391380787,
		0.134415255243784220360,
		0.0929271953151245376859,
		0.0516032829970797396969,
		0.0170017196299402603390
	};
	static ityp w_031[31] =
	{
		0.00254478079156187441540,
		0.00843456573932110624631,
		0.0164460498543878109338,
		0.0258075980961766535646,
		0.0359571033071293220968,
		0.0464628932617579865414,
		0.0569795094941233574122,
		0.0672077542959907035404,
		0.0768796204990035310427,
		0.0857559200499903511542,
		0.0936271099812644736167,
		0.100314278611795578771,
		0.105669893580234809744,
		0.109578421055924638237,
		0.111956873020953456880,
		0.112755256720768691607,
		0.111956873020953456880,
		0.109578421055924638237,
		0.105669893580234809744,
		0.100314278611795578771,
		0.0936271099812644736167,
		0.0857559200499903511542,
		0.0768796204990035310427,
		0.0672077542959907035404,
		0.0569795094941233574122,
		0.0464628932617579865414,
		0.0359571033071293220968,
		0.0258075980961766535646,
		0.0164460498543878109338,
		0.00843456573932110624631,
		0.00254478079156187441540
	};
	static ityp w_063[63] =
	{
		0.000363221481845530659694,
		0.00126515655623006801137,
		0.00257904979468568827243,
		0.00421763044155885483908,
		0.00611550682211724633968,
		0.00822300795723592966926,
		0.0104982469096213218983,
		0.0129038001003512656260,
		0.0154067504665594978021,
		0.0179785515681282703329,
		0.0205942339159127111492,
		0.0232314466399102694433,
		0.0258696793272147469108,
		0.0284897547458335486125,
		0.0310735511116879648799,
		0.0336038771482077305417,
		0.0360644327807825726401,
		0.0384398102494555320386,
		0.0407155101169443189339,
		0.0428779600250077344929,
		0.0449145316536321974143,
		0.0468135549906280124026,
		0.0485643304066731987159,
		0.0501571393058995374137,
		0.0515832539520484587768,
		0.0528349467901165198621,
		0.0539054993352660639269,
		0.0547892105279628650322,
		0.0554814043565593639878,
		0.0559784365104763194076,
		0.0562776998312543012726,
		0.0563776283603847173877,
		0.0562776998312543012726,
		0.0559784365104763194076,
		0.0554814043565593639878,
		0.0547892105279628650322,
		0.0539054993352660639269,
		0.0528349467901165198621,
		0.0515832539520484587768,
		0.0501571393058995374137,
		0.0485643304066731987159,
		0.0468135549906280124026,
		0.0449145316536321974143,
		0.0428779600250077344929,
		0.0407155101169443189339,
		0.0384398102494555320386,
		0.0360644327807825726401,
		0.0336038771482077305417,
		0.0310735511116879648799,
		0.0284897547458335486125,
		0.0258696793272147469108,
		0.0232314466399102694433,
		0.0205942339159127111492,
		0.0179785515681282703329,
		0.0154067504665594978021,
		0.0129038001003512656260,
		0.0104982469096213218983,
		0.00822300795723592966926,
		0.00611550682211724633968,
		0.00421763044155885483908,
		0.00257904979468568827243,
		0.00126515655623006801137,
		0.000363221481845530659694
	};
	static ityp w_127[127] =
	{
		0.0000505360952078625176247,
		0.000180739564445388357820,
		0.000377746646326984660274,
		0.000632607319362633544219,
		0.000938369848542381500794,
		0.00128952408261041739210,
		0.00168114286542146990631,
		0.00210881524572663287933,
		0.00256876494379402037313,
		0.00305775341017553113613,
		0.00357289278351729964938,
		0.00411150397865469304717,
		0.00467105037211432174741,
		0.00524912345480885912513,
		0.00584344987583563950756,
		0.00645190005017573692280,
		0.00707248999543355546805,
		0.00770337523327974184817,
		0.00834283875396815770558,
		0.00898927578406413572328,
		0.00964117772970253669530,
		0.0102971169579563555237,
		0.0109557333878379016480,
		0.0116157233199551347270,
		0.0122758305600827700870,
		0.0129348396636073734547,
		0.0135915710097655467896,
		0.0142448773729167743063,
		0.0148936416648151820348,
		0.0155367755558439824399,
		0.0161732187295777199419,
		0.0168019385741038652709,
		0.0174219301594641737472,
		0.0180322163903912863201,
		0.0186318482561387901863,
		0.0192199051247277660193,
		0.0197954950480974994880,
		0.0203577550584721594669,
		0.0209058514458120238522,
		0.0214389800125038672465,
		0.0219563663053178249393,
		0.0224572658268160987071,
		0.0229409642293877487608,
		0.0234067774953140062013,
		0.0238540521060385400804,
		0.0242821652033365993580,
		0.0246905247444876769091,
		0.0250785696529497687068,
		0.0254457699654647658126,
		0.0257916269760242293884,
		0.0261156733767060976805,
		0.0264174733950582599310,
		0.0266966229274503599062,
		0.0269527496676330319634,
		0.0271855132296247918192,
		0.0273946052639814325161,
		0.0275797495664818730349,
		0.0277407021782796819939,
		0.0278772514766137016085,
		0.0279892182552381597038,
		0.0280764557938172466068,
		0.0281388499156271506363,
		0.0281763190330166021307,
		0.0281888141801923586938,
		0.0281763190330166021307,
		0.0281388499156271506363,
		0.0280764557938172466068,
		0.0279892182552381597038,
		0.0278772514766137016085,
		0.0277407021782796819939,
		0.0275797495664818730349,
		0.0273946052639814325161,
		0.0271855132296247918192,
		0.0269527496676330319634,
		0.0266966229274503599062,
		0.0264174733950582599310,
		0.0261156733767060976805,
		0.0257916269760242293884,
		0.0254457699654647658126,
		0.0250785696529497687068,
		0.0246905247444876769091,
		0.0242821652033365993580,
		0.0238540521060385400804,
		0.0234067774953140062013,
		0.0229409642293877487608,
		0.0224572658268160987071,
		0.0219563663053178249393,
		0.0214389800125038672465,
		0.0209058514458120238522,
		0.0203577550584721594669,
		0.0197954950480974994880,
		0.0192199051247277660193,
		0.0186318482561387901863,
		0.0180322163903912863201,
		0.0174219301594641737472,
		0.0168019385741038652709,
		0.0161732187295777199419,
		0.0155367755558439824399,
		0.0148936416648151820348,
		0.0142448773729167743063,
		0.0135915710097655467896,
		0.0129348396636073734547,
		0.0122758305600827700870,
		0.0116157233199551347270,
		0.0109557333878379016480,
		0.0102971169579563555237,
		0.00964117772970253669530,
		0.00898927578406413572328,
		0.00834283875396815770558,
		0.00770337523327974184817,
		0.00707248999543355546805,
		0.00645190005017573692280,
		0.00584344987583563950756,
		0.00524912345480885912513,
		0.00467105037211432174741,
		0.00411150397865469304717,
		0.00357289278351729964938,
		0.00305775341017553113613,
		0.00256876494379402037313,
		0.00210881524572663287933,
		0.00168114286542146990631,
		0.00128952408261041739210,
		0.000938369848542381500794,
		0.000632607319362633544219,
		0.000377746646326984660274,
		0.000180739564445388357820,
		0.0000505360952078625176247
	};
	static ityp w_255[255] =
	{
		0.69379364324108267170E-05,
		0.25157870384280661489E-04,
		0.53275293669780613125E-04,
		0.90372734658751149261E-04,
		0.13575491094922871973E-03,
		0.18887326450650491366E-03,
		0.24921240048299729402E-03,
		0.31630366082226447689E-03,
		0.38974528447328229322E-03,
		0.46918492424785040975E-03,
		0.55429531493037471492E-03,
		0.64476204130572477933E-03,
		0.74028280424450333046E-03,
		0.84057143271072246365E-03,
		0.94536151685852538246E-03,
		0.10544076228633167722E-02,
		0.11674841174299594077E-02,
		0.12843824718970101768E-02,
		0.14049079956551446427E-02,
		0.15288767050877655684E-02,
		0.16561127281544526052E-02,
		0.17864463917586498247E-02,
		0.19197129710138724125E-02,
		0.20557519893273465236E-02,
		0.21944069253638388388E-02,
		0.23355251860571608737E-02,
		0.24789582266575679307E-02,
		0.26245617274044295626E-02,
		0.27721957645934509940E-02,
		0.29217249379178197538E-02,
		0.30730184347025783234E-02,
		0.32259500250878684614E-02,
		0.33803979910869203823E-02,
		0.35362449977167777340E-02,
		0.36933779170256508183E-02,
		0.38516876166398709241E-02,
		0.40110687240750233989E-02,
		0.41714193769840788528E-02,
		0.43326409680929828545E-02,
		0.44946378920320678616E-02,
		0.46573172997568547773E-02,
		0.48205888648512683476E-02,
		0.49843645647655386012E-02,
		0.51485584789781777618E-02,
		0.53130866051870565663E-02,
		0.54778666939189508240E-02,
		0.56428181013844441585E-02,
		0.58078616599775673635E-02,
		0.59729195655081658049E-02,
		0.61379152800413850435E-02,
		0.63027734490857587172E-02,
		0.64674198318036867274E-02,
		0.66317812429018878941E-02,
		0.67957855048827733948E-02,
		0.69593614093904229394E-02,
		0.71224386864583871532E-02,
		0.72849479805538070639E-02,
		0.74468208324075910174E-02,
		0.76079896657190565832E-02,
		0.77683877779219912200E-02,
		0.79279493342948491103E-02,
		0.80866093647888599710E-02,
		0.82443037630328680306E-02,
		0.84009692870519326354E-02,
		0.85565435613076896192E-02,
		0.87109650797320868736E-02,
		0.88641732094824942641E-02,
		0.90161081951956431600E-02,
		0.91667111635607884067E-02,
		0.93159241280693950932E-02,
		0.94636899938300652943E-02,
		0.96099525623638830097E-02,
		0.97546565363174114611E-02,
		0.98977475240487497440E-02,
		0.10039172044056840798E-01,
		0.10178877529236079733E-01,
		0.10316812330947621682E-01,
		0.10452925722906011926E-01,
		0.10587167904885197931E-01,
		0.10719490006251933623E-01,
		0.10849844089337314099E-01,
		0.10978183152658912470E-01,
		0.11104461134006926537E-01,
		0.11228632913408049354E-01,
		0.11350654315980596602E-01,
		0.11470482114693874380E-01,
		0.11588074033043952568E-01,
		0.11703388747657003101E-01,
		0.11816385890830235763E-01,
		0.11927026053019270040E-01,
		0.12035270785279562630E-01,
		0.12141082601668299679E-01,
		0.12244424981611985899E-01,
		0.12345262372243838455E-01,
		0.12443560190714035263E-01,
		0.12539284826474884353E-01,
		0.12632403643542078765E-01,
		0.12722884982732382906E-01,
		0.12810698163877361967E-01,
		0.12895813488012114694E-01,
		0.12978202239537399286E-01,
		0.13057836688353048840E-01,
		0.13134690091960152836E-01,
		0.13208736697529129966E-01,
		0.13279951743930530650E-01,
		0.13348311463725179953E-01,
		0.13413793085110098513E-01,
		0.13476374833816515982E-01,
		0.13536035934956213614E-01,
		0.13592756614812395910E-01,
		0.13646518102571291428E-01,
		0.13697302631990716258E-01,
		0.13745093443001896632E-01,
		0.13789874783240936517E-01,
		0.13831631909506428676E-01,
		0.13870351089139840997E-01,
		0.13906019601325461264E-01,
		0.13938625738306850804E-01,
		0.13968158806516938516E-01,
		0.13994609127619079852E-01,
		0.14017968039456608810E-01,
		0.14038227896908623303E-01,
		0.14055382072649964277E-01,
		0.14069424957813575318E-01,
		0.14080351962553661325E-01,
		0.14088159516508301065E-01,
		0.14092845069160408355E-01,
		0.14094407090096179347E-01,
		0.14092845069160408355E-01,
		0.14088159516508301065E-01,
		0.14080351962553661325E-01,
		0.14069424957813575318E-01,
		0.14055382072649964277E-01,
		0.14038227896908623303E-01,
		0.14017968039456608810E-01,
		0.13994609127619079852E-01,
		0.13968158806516938516E-01,
		0.13938625738306850804E-01,
		0.13906019601325461264E-01,
		0.13870351089139840997E-01,
		0.13831631909506428676E-01,
		0.13789874783240936517E-01,
		0.13745093443001896632E-01,
		0.13697302631990716258E-01,
		0.13646518102571291428E-01,
		0.13592756614812395910E-01,
		0.13536035934956213614E-01,
		0.13476374833816515982E-01,
		0.13413793085110098513E-01,
		0.13348311463725179953E-01,
		0.13279951743930530650E-01,
		0.13208736697529129966E-01,
		0.13134690091960152836E-01,
		0.13057836688353048840E-01,
		0.12978202239537399286E-01,
		0.12895813488012114694E-01,
		0.12810698163877361967E-01,
		0.12722884982732382906E-01,
		0.12632403643542078765E-01,
		0.12539284826474884353E-01,
		0.12443560190714035263E-01,
		0.12345262372243838455E-01,
		0.12244424981611985899E-01,
		0.12141082601668299679E-01,
		0.12035270785279562630E-01,
		0.11927026053019270040E-01,
		0.11816385890830235763E-01,
		0.11703388747657003101E-01,
		0.11588074033043952568E-01,
		0.11470482114693874380E-01,
		0.11350654315980596602E-01,
		0.11228632913408049354E-01,
		0.11104461134006926537E-01,
		0.10978183152658912470E-01,
		0.10849844089337314099E-01,
		0.10719490006251933623E-01,
		0.10587167904885197931E-01,
		0.10452925722906011926E-01,
		0.10316812330947621682E-01,
		0.10178877529236079733E-01,
		0.10039172044056840798E-01,
		0.98977475240487497440E-02,
		0.97546565363174114611E-02,
		0.96099525623638830097E-02,
		0.94636899938300652943E-02,
		0.93159241280693950932E-02,
		0.91667111635607884067E-02,
		0.90161081951956431600E-02,
		0.88641732094824942641E-02,
		0.87109650797320868736E-02,
		0.85565435613076896192E-02,
		0.84009692870519326354E-02,
		0.82443037630328680306E-02,
		0.80866093647888599710E-02,
		0.79279493342948491103E-02,
		0.77683877779219912200E-02,
		0.76079896657190565832E-02,
		0.74468208324075910174E-02,
		0.72849479805538070639E-02,
		0.71224386864583871532E-02,
		0.69593614093904229394E-02,
		0.67957855048827733948E-02,
		0.66317812429018878941E-02,
		0.64674198318036867274E-02,
		0.63027734490857587172E-02,
		0.61379152800413850435E-02,
		0.59729195655081658049E-02,
		0.58078616599775673635E-02,
		0.56428181013844441585E-02,
		0.54778666939189508240E-02,
		0.53130866051870565663E-02,
		0.51485584789781777618E-02,
		0.49843645647655386012E-02,
		0.48205888648512683476E-02,
		0.46573172997568547773E-02,
		0.44946378920320678616E-02,
		0.43326409680929828545E-02,
		0.41714193769840788528E-02,
		0.40110687240750233989E-02,
		0.38516876166398709241E-02,
		0.36933779170256508183E-02,
		0.35362449977167777340E-02,
		0.33803979910869203823E-02,
		0.32259500250878684614E-02,
		0.30730184347025783234E-02,
		0.29217249379178197538E-02,
		0.27721957645934509940E-02,
		0.26245617274044295626E-02,
		0.24789582266575679307E-02,
		0.23355251860571608737E-02,
		0.21944069253638388388E-02,
		0.20557519893273465236E-02,
		0.19197129710138724125E-02,
		0.17864463917586498247E-02,
		0.16561127281544526052E-02,
		0.15288767050877655684E-02,
		0.14049079956551446427E-02,
		0.12843824718970101768E-02,
		0.11674841174299594077E-02,
		0.10544076228633167722E-02,
		0.94536151685852538246E-03,
		0.84057143271072246365E-03,
		0.74028280424450333046E-03,
		0.64476204130572477933E-03,
		0.55429531493037471492E-03,
		0.46918492424785040975E-03,
		0.38974528447328229322E-03,
		0.31630366082226447689E-03,
		0.24921240048299729402E-03,
		0.18887326450650491366E-03,
		0.13575491094922871973E-03,
		0.90372734658751149261E-04,
		0.53275293669780613125E-04,
		0.25157870384280661489E-04,
		0.69379364324108267170E-05
	};
	
	switch(order)
	{
		case 1: 
			r8vec_copy ( order, w_001, w );
			break;
		case 3: 
			r8vec_copy ( order, w_003, w );
			break;
		case 7:
			r8vec_copy ( order, w_007, w );
			break;
		case 15:
			r8vec_copy ( order, w_015, w );
			break;
		case 31:
			r8vec_copy ( order, w_031, w );
			break;
		case 63:
			r8vec_copy ( order, w_063, w );
			break;
		case 127:
			r8vec_copy ( order, w_127, w );
			break;
		case 255:
			r8vec_copy ( order, w_255, w );
			break;
		default:
			return NULL;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _patterson_lookup_weights_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    PATTERSON_LOOKUP_WEIGHTS_NP looks up Patterson quadrature weights.
  Discussion:
    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.
    The weights are positive, symmetric and should sum to 2.
    The user must preallocate space for the output array W.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  Parameters:
    Input, int ORDER, the order of the rule.
    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double W[ORDER], the weights.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * w = s_data->a3;
	
	patterson_lookup_weights ( order, w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _patterson_lookup ( void * data)
/******************************************************************************/
/*
  Purpose:
    PATTERSON_LOOKUP looks up Patterson quadrature points and weights.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    The rule is defined on [-1,1],
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    19 July 2010
  Author:
    John Burkardt
  Reference:
    Prem Kythe, Michael Schaeferkotter,
    Handbook of Computational Methods for Integration,
    Chapman and Hall, 2004,
    ISBN: 1-58488-428-2,
    LC: QA299.3.K98.
    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.
  Parameters:
    Input, int N, the order.
    Legal values are 1, 3, 7, 15, 31, 63, 127 and 255.
    Output, double X[N], the abscissas.
    Output, double W[N], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0; 
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
	patterson_lookup_points ( n, x );
	patterson_lookup_weights ( n, w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _patterson_lookup_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    PATTERSON_LOOKUP_POINTS looks up Patterson quadrature points.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    The rule is defined on [-1,1],
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    17 December 2009
  Author:
    John Burkardt
  Reference:
    Prem Kythe, Michael Schaeferkotter,
    Handbook of Computational Methods for Integration,
    Chapman and Hall, 2004,
    ISBN: 1-58488-428-2,
    LC: QA299.3.K98.
    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.
  Parameters:
    Input, int ORDER, the order of the rule.
    Legal values are 1, 3, 7, 15, 31, 63, 127 and 255.
    Output, double X[ORDER], the abscissas.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	
	static ityp x_001[1] = 
	{
		0.00 
	};
	static ityp x_003[3] =
	{
		-0.77459666924148337704, 
		0.00, 
		0.77459666924148337704
	};
	static ityp x_007[7] =
	{
		-0.96049126870802028342, 
		-0.77459666924148337704, 
		-0.43424374934680255800, 
		0.00, 
		0.43424374934680255800, 
		0.77459666924148337704, 
		0.96049126870802028342 
	};
	static ityp x_015[15] =
	{
		-0.99383196321275502221, 
		-0.96049126870802028342, 
		-0.88845923287225699889, 
		-0.77459666924148337704, 
		-0.62110294673722640294, 
		-0.43424374934680255800, 
		-0.22338668642896688163, 
		0.00, 
		0.22338668642896688163, 
		0.43424374934680255800, 
		0.62110294673722640294, 
		0.77459666924148337704, 
		0.88845923287225699889, 
		0.96049126870802028342, 
		0.99383196321275502221
	};
	static ityp x_031[31] =
	{
		-0.99909812496766759766, 
		-0.99383196321275502221, 
		-0.98153114955374010687, 
		-0.96049126870802028342, 
		-0.92965485742974005667, 
		-0.88845923287225699889, 
		-0.83672593816886873550, 
		-0.77459666924148337704, 
		-0.70249620649152707861, 
		-0.62110294673722640294, 
		-0.53131974364437562397, 
		-0.43424374934680255800, 
		-0.33113539325797683309, 
		-0.22338668642896688163, 
		-0.11248894313318662575, 
		0.00, 
		0.11248894313318662575, 
		0.22338668642896688163, 
		0.33113539325797683309, 
		0.43424374934680255800, 
		0.53131974364437562397, 
		0.62110294673722640294, 
		0.70249620649152707861, 
		0.77459666924148337704, 
		0.83672593816886873550, 
		0.88845923287225699889, 
		0.92965485742974005667, 
		0.96049126870802028342, 
		0.98153114955374010687, 
		0.99383196321275502221, 
		0.99909812496766759766
	};
	static ityp x_063[63] =
	{
		-0.99987288812035761194, 
		-0.99909812496766759766, 
		-0.99720625937222195908, 
		-0.99383196321275502221, 
		-0.98868475754742947994, 
		-0.98153114955374010687, 
		-0.97218287474858179658, 
		-0.96049126870802028342, 
		-0.94634285837340290515, 
		-0.92965485742974005667, 
		-0.91037115695700429250, 
		-0.88845923287225699889, 
		-0.86390793819369047715, 
		-0.83672593816886873550, 
		-0.80694053195021761186, 
		-0.77459666924148337704, 
		-0.73975604435269475868, 
		-0.70249620649152707861, 
		-0.66290966002478059546, 
		-0.62110294673722640294, 
		-0.57719571005204581484, 
		-0.53131974364437562397, 
		-0.48361802694584102756, 
		-0.43424374934680255800, 
		-0.38335932419873034692, 
		-0.33113539325797683309, 
		-0.27774982202182431507, 
		-0.22338668642896688163, 
		-0.16823525155220746498, 
		-0.11248894313318662575, 
		-0.056344313046592789972, 
		0.00, 
		0.056344313046592789972, 
		0.11248894313318662575, 
		0.16823525155220746498, 
		0.22338668642896688163, 
		0.27774982202182431507, 
		0.33113539325797683309, 
		0.38335932419873034692, 
		0.43424374934680255800, 
		0.48361802694584102756, 
		0.53131974364437562397, 
		0.57719571005204581484, 
		0.62110294673722640294, 
		0.66290966002478059546, 
		0.70249620649152707861, 
		0.73975604435269475868, 
		0.77459666924148337704, 
		0.80694053195021761186, 
		0.83672593816886873550, 
		0.86390793819369047715, 
		0.88845923287225699889, 
		0.91037115695700429250, 
		0.92965485742974005667, 
		0.94634285837340290515, 
		0.96049126870802028342, 
		0.97218287474858179658, 
		0.98153114955374010687, 
		0.98868475754742947994, 
		0.99383196321275502221, 
		0.99720625937222195908, 
		0.99909812496766759766, 
		0.99987288812035761194
	};
	static ityp x_127[127] =
	{
		-0.99998243035489159858, 
		-0.99987288812035761194, 
		-0.99959879967191068325, 
		-0.99909812496766759766, 
		-0.99831663531840739253, 
		-0.99720625937222195908, 
		-0.99572410469840718851, 
		-0.99383196321275502221, 
		-0.99149572117810613240, 
		-0.98868475754742947994, 
		-0.98537149959852037111, 
		-0.98153114955374010687, 
		-0.97714151463970571416, 
		-0.97218287474858179658, 
		-0.96663785155841656709, 
		-0.96049126870802028342, 
		-0.95373000642576113641, 
		-0.94634285837340290515, 
		-0.93832039777959288365, 
		-0.92965485742974005667, 
		-0.92034002547001242073, 
		-0.91037115695700429250, 
		-0.89974489977694003664, 
		-0.88845923287225699889, 
		-0.87651341448470526974, 
		-0.86390793819369047715, 
		-0.85064449476835027976, 
		-0.83672593816886873550, 
		-0.82215625436498040737, 
		-0.80694053195021761186, 
		-0.79108493379984836143, 
		-0.77459666924148337704, 
		-0.75748396638051363793, 
		-0.73975604435269475868, 
		-0.72142308537009891548, 
		-0.70249620649152707861, 
		-0.68298743109107922809, 
		-0.66290966002478059546, 
		-0.64227664250975951377, 
		-0.62110294673722640294, 
		-0.59940393024224289297, 
		-0.57719571005204581484, 
		-0.55449513263193254887, 
		-0.53131974364437562397, 
		-0.50768775753371660215, 
		-0.48361802694584102756, 
		-0.45913001198983233287, 
		-0.43424374934680255800, 
		-0.40897982122988867241, 
		-0.38335932419873034692, 
		-0.35740383783153215238, 
		-0.33113539325797683309, 
		-0.30457644155671404334, 
		-0.27774982202182431507, 
		-0.25067873030348317661, 
		-0.22338668642896688163, 
		-0.19589750271110015392, 
		-0.16823525155220746498, 
		-0.14042423315256017459, 
		-0.11248894313318662575, 
		-0.084454040083710883710, 
		-0.056344313046592789972, 
		-0.028184648949745694339, 
		0.00, 
		0.028184648949745694339, 
		0.056344313046592789972, 
		0.084454040083710883710, 
		0.11248894313318662575, 
		0.14042423315256017459, 
		0.16823525155220746498, 
		0.19589750271110015392, 
		0.22338668642896688163, 
		0.25067873030348317661, 
		0.27774982202182431507, 
		0.30457644155671404334, 
		0.33113539325797683309, 
		0.35740383783153215238, 
		0.38335932419873034692, 
		0.40897982122988867241, 
		0.43424374934680255800, 
		0.45913001198983233287, 
		0.48361802694584102756, 
		0.50768775753371660215, 
		0.53131974364437562397, 
		0.55449513263193254887, 
		0.57719571005204581484, 
		0.59940393024224289297, 
		0.62110294673722640294, 
		0.64227664250975951377, 
		0.66290966002478059546, 
		0.68298743109107922809, 
		0.70249620649152707861, 
		0.72142308537009891548, 
		0.73975604435269475868, 
		0.75748396638051363793, 
		0.77459666924148337704, 
		0.79108493379984836143, 
		0.80694053195021761186, 
		0.82215625436498040737, 
		0.83672593816886873550, 
		0.85064449476835027976, 
		0.86390793819369047715, 
		0.87651341448470526974, 
		0.88845923287225699889, 
		0.89974489977694003664, 
		0.91037115695700429250, 
		0.92034002547001242073, 
		0.92965485742974005667, 
		0.93832039777959288365, 
		0.94634285837340290515, 
		0.95373000642576113641, 
		0.96049126870802028342, 
		0.96663785155841656709, 
		0.97218287474858179658, 
		0.97714151463970571416, 
		0.98153114955374010687, 
		0.98537149959852037111, 
		0.98868475754742947994, 
		0.99149572117810613240, 
		0.99383196321275502221, 
		0.99572410469840718851, 
		0.99720625937222195908, 
		0.99831663531840739253, 
		0.99909812496766759766, 
		0.99959879967191068325, 
		0.99987288812035761194, 
		0.99998243035489159858
	};
	static ityp x_255[255] =
	{
		-0.99999759637974846462, 
		-0.99998243035489159858, 
		-0.99994399620705437576, 
		-0.99987288812035761194, 
		-0.99976049092443204733, 
		-0.99959879967191068325, 
		-0.99938033802502358193, 
		-0.99909812496766759766, 
		-0.99874561446809511470, 
		-0.99831663531840739253, 
		-0.99780535449595727456, 
		-0.99720625937222195908, 
		-0.99651414591489027385, 
		-0.99572410469840718851, 
		-0.99483150280062100052, 
		-0.99383196321275502221, 
		-0.99272134428278861533, 
		-0.99149572117810613240, 
		-0.99015137040077015918, 
		-0.98868475754742947994, 
		-0.98709252795403406719, 
		-0.98537149959852037111, 
		-0.98351865757863272876, 
		-0.98153114955374010687, 
		-0.97940628167086268381, 
		-0.97714151463970571416, 
		-0.97473445975240266776, 
		-0.97218287474858179658, 
		-0.96948465950245923177, 
		-0.96663785155841656709, 
		-0.96364062156981213252, 
		-0.96049126870802028342, 
		-0.95718821610986096274, 
		-0.95373000642576113641, 
		-0.95011529752129487656, 
		-0.94634285837340290515, 
		-0.94241156519108305981, 
		-0.93832039777959288365, 
		-0.93406843615772578800, 
		-0.92965485742974005667, 
		-0.92507893290707565236, 
		-0.92034002547001242073, 
		-0.91543758715576504064, 
		-0.91037115695700429250, 
		-0.90514035881326159519, 
		-0.89974489977694003664, 
		-0.89418456833555902286, 
		-0.88845923287225699889, 
		-0.88256884024734190684, 
		-0.87651341448470526974, 
		-0.87029305554811390585, 
		-0.86390793819369047715, 
		-0.85735831088623215653, 
		-0.85064449476835027976, 
		-0.84376688267270860104, 
		-0.83672593816886873550, 
		-0.82952219463740140018, 
		-0.82215625436498040737, 
		-0.81462878765513741344, 
		-0.80694053195021761186, 
		-0.79909229096084140180, 
		-0.79108493379984836143, 
		-0.78291939411828301639, 
		-0.77459666924148337704, 
		-0.76611781930376009072, 
		-0.75748396638051363793, 
		-0.74869629361693660282, 
		-0.73975604435269475868, 
		-0.73066452124218126133, 
		-0.72142308537009891548, 
		-0.71203315536225203459, 
		-0.70249620649152707861, 
		-0.69281376977911470289, 
		-0.68298743109107922809, 
		-0.67301883023041847920, 
		-0.66290966002478059546, 
		-0.65266166541001749610, 
		-0.64227664250975951377, 
		-0.63175643771119423041, 
		-0.62110294673722640294, 
		-0.61031811371518640016, 
		-0.59940393024224289297, 
		-0.58836243444766254143, 
		-0.57719571005204581484, 
		-0.56590588542365442262, 
		-0.55449513263193254887, 
		-0.54296566649831149049, 
		-0.53131974364437562397, 
		-0.51955966153745702199, 
		-0.50768775753371660215, 
		-0.49570640791876146017, 
		-0.48361802694584102756, 
		-0.47142506587165887693, 
		-0.45913001198983233287, 
		-0.44673538766202847374, 
		-0.43424374934680255800, 
		-0.42165768662616330006, 
		-0.40897982122988867241, 
		-0.39621280605761593918, 
		-0.38335932419873034692, 
		-0.37042208795007823014, 
		-0.35740383783153215238, 
		-0.34430734159943802278, 
		-0.33113539325797683309, 
		-0.31789081206847668318, 
		-0.30457644155671404334, 
		-0.29119514851824668196, 
		-0.27774982202182431507, 
		-0.26424337241092676194, 
		-0.25067873030348317661, 
		-0.23705884558982972721, 
		-0.22338668642896688163, 
		-0.20966523824318119477, 
		-0.19589750271110015392, 
		-0.18208649675925219825, 
		-0.16823525155220746498, 
		-0.15434681148137810869, 
		-0.14042423315256017459, 
		-0.12647058437230196685, 
		-0.11248894313318662575, 
		-0.098482396598119202090, 
		-0.084454040083710883710, 
		-0.070406976042855179063, 
		-0.056344313046592789972, 
		-0.042269164765363603212, 
		-0.028184648949745694339, 
		-0.014093886410782462614, 
		0.00, 
		0.014093886410782462614, 
		0.028184648949745694339, 
		0.042269164765363603212, 
		0.056344313046592789972, 
		0.070406976042855179063, 
		0.084454040083710883710, 
		0.098482396598119202090, 
		0.11248894313318662575, 
		0.12647058437230196685, 
		0.14042423315256017459, 
		0.15434681148137810869, 
		0.16823525155220746498, 
		0.18208649675925219825, 
		0.19589750271110015392, 
		0.20966523824318119477, 
		0.22338668642896688163, 
		0.23705884558982972721, 
		0.25067873030348317661, 
		0.26424337241092676194, 
		0.27774982202182431507, 
		0.29119514851824668196, 
		0.30457644155671404334, 
		0.31789081206847668318, 
		0.33113539325797683309, 
		0.34430734159943802278, 
		0.35740383783153215238, 
		0.37042208795007823014, 
		0.38335932419873034692, 
		0.39621280605761593918, 
		0.40897982122988867241, 
		0.42165768662616330006, 
		0.43424374934680255800, 
		0.44673538766202847374, 
		0.45913001198983233287, 
		0.47142506587165887693, 
		0.48361802694584102756, 
		0.49570640791876146017, 
		0.50768775753371660215, 
		0.51955966153745702199, 
		0.53131974364437562397, 
		0.54296566649831149049, 
		0.55449513263193254887, 
		0.56590588542365442262, 
		0.57719571005204581484, 
		0.58836243444766254143, 
		0.59940393024224289297, 
		0.61031811371518640016, 
		0.62110294673722640294, 
		0.63175643771119423041, 
		0.64227664250975951377, 
		0.65266166541001749610, 
		0.66290966002478059546, 
		0.67301883023041847920, 
		0.68298743109107922809, 
		0.69281376977911470289, 
		0.70249620649152707861, 
		0.71203315536225203459, 
		0.72142308537009891548, 
		0.73066452124218126133, 
		0.73975604435269475868, 
		0.74869629361693660282, 
		0.75748396638051363793, 
		0.76611781930376009072, 
		0.77459666924148337704, 
		0.78291939411828301639, 
		0.79108493379984836143, 
		0.79909229096084140180, 
		0.80694053195021761186, 
		0.81462878765513741344, 
		0.82215625436498040737, 
		0.82952219463740140018, 
		0.83672593816886873550, 
		0.84376688267270860104, 
		0.85064449476835027976, 
		0.85735831088623215653, 
		0.86390793819369047715, 
		0.87029305554811390585, 
		0.87651341448470526974, 
		0.88256884024734190684, 
		0.88845923287225699889, 
		0.89418456833555902286, 
		0.89974489977694003664, 
		0.90514035881326159519, 
		0.91037115695700429250, 
		0.91543758715576504064, 
		0.92034002547001242073, 
		0.92507893290707565236, 
		0.92965485742974005667, 
		0.93406843615772578800, 
		0.93832039777959288365, 
		0.94241156519108305981, 
		0.94634285837340290515, 
		0.95011529752129487656, 
		0.95373000642576113641, 
		0.95718821610986096274, 
		0.96049126870802028342, 
		0.96364062156981213252, 
		0.96663785155841656709, 
		0.96948465950245923177, 
		0.97218287474858179658, 
		0.97473445975240266776, 
		0.97714151463970571416, 
		0.97940628167086268381, 
		0.98153114955374010687, 
		0.98351865757863272876, 
		0.98537149959852037111, 
		0.98709252795403406719, 
		0.98868475754742947994, 
		0.99015137040077015918, 
		0.99149572117810613240, 
		0.99272134428278861533, 
		0.99383196321275502221, 
		0.99483150280062100052, 
		0.99572410469840718851, 
		0.99651414591489027385, 
		0.99720625937222195908, 
		0.99780535449595727456, 
		0.99831663531840739253, 
		0.99874561446809511470, 
		0.99909812496766759766, 
		0.99938033802502358193, 
		0.99959879967191068325, 
		0.99976049092443204733, 
		0.99987288812035761194, 
		0.99994399620705437576, 
		0.99998243035489159858, 
		0.99999759637974846462 
	};
	
	switch(order)
	{
		case 1: 
			r8vec_copy ( order, x_001, x );
			break;
		case 3: 
			r8vec_copy ( order, x_003, x );
			break;
		case 7:
			r8vec_copy ( order, x_007, x );
			break;
		case 15:
			r8vec_copy ( order, x_015, x );
			break;
		case 31:
			r8vec_copy ( order, x_031, x );
			break;
		case 63:
			r8vec_copy ( order, x_063, x );
			break;
		case 127:
			r8vec_copy ( order, x_127, x );
			break;
		case 255:
			r8vec_copy ( order, x_255, x );
			break;
		default:
			return NULL;
	}
	return NULL; 
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _patterson_lookup_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    PATTERSON_LOOKUP_POINTS_NP looks up Patterson quadrature points.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    The rule is defined on [-1,1],
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    17 December 2009
  Author:
    John Burkardt
  Reference:
    Prem Kythe, Michael Schaeferkotter,
    Handbook of Computational Methods for Integration,
    Chapman and Hall, 2004,
    ISBN: 1-58488-428-2,
    LC: QA299.3.K98.
    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.
  Parameters:
    Input, int ORDER, the order of the rule.
    Legal values are 1, 3, 7, 15, 31, 63, 127 and 255.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	patterson_lookup_points ( order, x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _jacobi_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_COMPUTE computes a Jacobi quadrature rule.
  Discussion:
    The integration interval is [ -1, 1 ].
    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.
    The integral to approximate:
      Integral ( -1 <= X <= 1 ) (1-X)^ALPHA * (1+X)^BETA * F(X) dX
    The quadrature rule:
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
    Thanks to Xu Xiang of Fudan University for pointing out that
    an earlier implementation of this routine was incorrect!
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  Reference:
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  Parameters:
    Input, int ORDER, the order of the rule.
    1 <= ORDER.
    Input, double ALPHA, BETA, the exponents of (1-X) and
 (1+X) in the quadrature rule.  For simple Legendre quadrature,
    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
    Output, double X[ORDER], the abscissas.
    Output, double W[ORDER], the weights.
*/
{
	const dt2it2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register ityp alpha = s_data->a1; 
	const register ityp beta = s_data->a2;
	ityp * x = s_data->a3;
	ityp * w = s_data->a4;
	
	ityp an;
	ityp *b;
	ityp bn;
	ityp *c;
	ityp cc;
	ityp delta;
	ityp dp2;
	dim_typ i;
	ityp p1;
	ityp prod;
	ityp r1;
	ityp r2;
	ityp r3;
	ityp temp;
	ityp x0;
	
	if ( order == 0)
		return NULL;
	
	b = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	c = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	/*
	Check ALPHA and BETA.
	*/
	if ( alpha <= -1.00 || beta <= -1.00 )
		return NULL;
	
	/*
	Set the recursion coefficients.
	*/
	for ( i = 1; i <= order; i++ )
	{
		b[i-1] = alpha + beta == 0.0 || beta - alpha == 0.0 ? 0.00 : ( alpha + beta ) * ( beta - alpha ) / 	( ( alpha + beta + ( ityp ) ( i<<1 ) ) * ( alpha + beta + ( ityp ) ( (i<<1) - 2 ) ) );
		c[i-1] = 0.00 + (i!=1)*4.00 * ( ityp ) ( i - 1 ) * ( alpha + ( ityp ) ( i - 1 ) ) * ( beta + ( ityp ) ( i - 1 ) ) * ( alpha + beta + ( ityp ) ( i - 1 ) ) / ( ( alpha + beta + ( ityp ) ( (i<<1) - 1 ) ) * pow ( alpha + beta + ( ityp ) ( (i<<1) - 2 ), 2 ) * ( alpha + beta + ( ityp ) ( (i<<1) - 3 ) ) );
	}
	
	delta = r8_gamma ( alpha        + 1.00 ) * r8_gamma (         beta + 1.00 ) / r8_gamma ( alpha + beta + 2.00 );
	
	prod = 1.00;
	for ( i = 2; i <= order; ++i )
		prod *= c[i-1];
	cc = delta * pow ( 2.00, alpha + beta + 1.00 ) * prod;
	
	for ( i = 1; i <= order; ++i )
	{
		if ( i == 1 )
		{
			an = alpha / ( ityp ) ( order );
			bn = beta / ( ityp ) ( order );
			
			r1 = ( 1.00 + alpha ) * ( 2.78 / ( 4.00 + ( ityp ) ( order * order ) ) + 0.768 * an / ( ityp ) ( order ) );
			r2 = 1.00 + 1.48 * an + 0.96 * bn + 0.452 * an * an + 0.83 * an * bn;
			x0 = ( r2 - r1 ) / r2;
		}
		else if ( i == 2 )
		{
			r1 = ( 4.10 + alpha ) / 	( ( 1.00 + alpha ) * ( 1.00 + 0.156 * alpha ) );
			r2 = 1.00 + 0.06 * ( ( ityp ) ( order ) - 8.00 ) * ( 1.00 + 0.12 * alpha ) / ( ityp ) ( order );
			r3 = 1.00 + 0.012 * beta * ( 1.00 + 0.25 * abs ( alpha ) ) / ( ityp ) ( order );
			x0 = x0 - r1 * r2 * r3 * ( 1.00 - x0 );
		}
		else if ( i == 3 )
		{
			r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );
			r2 = 1.00 + 0.22 * ( ( ityp ) ( order ) - 8.00 ) / ( ityp ) ( order );
			r3 = 1.00 + 8.00 * beta / ( ( 6.28 + beta ) * ( ityp ) ( order * order ) );
			x0 = x0 - r1 * r2 * r3 * ( x[0] - x0 );
		}
		else if ( i < order - 1 )
			x0 = 3.00 * x[i-2] - 3.00 * x[i-3] + x[i-4];
		else if ( i == order - 1 )
		{
			r1 = ( 1.0 + 0.235 * beta ) / ( 0.766 + 0.119 * beta );
			r2 = 1.00 / ( 1.00 + 0.639 * ( ( ityp ) ( order ) - 4.00 ) / ( 1.00 + 0.71 * ( ( ityp ) ( order ) - 4.00 ) ) );
			r3 = 1.00 / ( 1.00 + 20.00 * alpha / ( ( 7.50 + alpha ) * ( ityp ) ( order * order ) ) );
			x0 = x0 + r1 * r2 * r3 * ( x0 - x[i-3] );
		}
		else if ( i == order )
		{
			r1 = ( 1.00 + 0.37 * beta ) / ( 1.67 + 0.28 * beta );
			r2 = 1.00 / 
			( 1.00 + 0.22 * ( ( ityp ) ( order ) - 8.00 ) / ( ityp ) ( order ) );
			
			r3 = 1.0 / ( 1.00 + 8.00 * alpha / ( ( 6.28 + alpha ) * ( ityp ) ( order * order ) ) );
			
			x0 = x0 + r1 * r2 * r3 * ( x0 - x[i-3] );
		}
		
		jacobi_root ( &x0, order, alpha, beta, &dp2, &p1, b, c );
	
		x[i-1] = x0;
		w[i-1] = cc / ( dp2 * p1 );
	}
	/*
	Reverse the order of the values.
	*/
	for ( i = 1; i <= order/2; ++i )
	{
		temp       = x[i-1];
		x[i-1]     = x[order-i];
		x[order-i] = temp;
	}
	
	for ( i = 1; i <=order/2; ++i )
	{
		temp       = w[i-1];
		w[i-1]     = w[order-i];
		w[order-i] = temp;
	}
	
	free ( b );
	free ( c );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _jacobi_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_COMPUTE_POINTS computes Jacobi quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, double ALPHA, BETA, the exponents of the (1-X) and (1+X) factors.
    Output, double X[ORDER], the abscissas.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register ityp alpha = s_data->a1; 
	const register ityp beta = s_data->a2; 
	ityp * x = s_data->a3;
	
	ityp *w = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	jacobi_compute ( order, alpha, beta, x, w );
	free ( w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _jacobi_compute_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_COMPUTE_POINTS_NP computes Jacobi quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameter values.
    P[0] = ALPHA, the exponent of (1-X)
    P[1] = BETA,  the exponent of (1+X).
    -1.0 < ALPHA and -1.0 < BETA are required.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	jacobi_compute_points ( order, p[0], p[1], x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gen_laguerre_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_LAGUERRE_COMPUTE computes a Generalized Laguerre quadrature rule.
  Discussion:
    In the simplest case, ALPHA is 0, and we are approximating the
    integral from 0 to +oo of exp(-X) * F(X).  When this is so,
    it is easy to modify the rule to approximate the integral from
    A to +oo as well.
    If ALPHA is nonzero, then there is no simple way to extend the
    rule to approximate the integral from A to +oo.  The simplest
    procedures would be to approximate the integral from 0 to A.
    The integration interval is [ A, +oo ) or [ 0, +oo ).
    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x^alpha.
    If the integral to approximate is:
        Integral ( A <= X < +oo ) exp ( - X ) * F(X) dX
      or
        Integral ( 0 <= X < +oo ) exp ( - X ) * X^ALPHA * F(X) dX
    then the quadrature rule is:
      exp ( - A ) * Sum ( 1 <= I <= ORDER ) W(I) * F ( A+X(I) )
    or
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
    If the integral to approximate is:
        Integral ( A <= X < +oo ) F(X) dX
      or
        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX
    then the quadrature rule is:
      exp ( - A ) * Sum ( 1 <= I <= ORDER ) 
        W(I) * exp(A+X(I)) * F ( A+X(I) )
    or
      Sum ( 1 <= I <= ORDER ) W(I) * exp(X(I)) * F ( X(I) )
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  Reference:
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  Parameters:
    Input, int ORDER, the order of the quadrature rule to be computed.
    1 <= ORDER.
    Input, double ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.
    Output, double X[ORDER], the abscissas.
    Output, double W[ORDER], the weights.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register ityp alpha = s_data->a1; 
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	ityp *b;
	ityp *c;
	ityp cc;
	ityp dp2;
	dim_typ i;
	ityp p1;
	ityp prod;
	ityp r1;
	ityp r2;
	ityp ratio;
	ityp x0;
	
	if ( order == 0)
		return NULL;
	
	b = ( ityp * ) malloc ( order * sizeof ( double ) );
	c = ( ityp * ) malloc ( order * sizeof ( double ) );
	/*
	Set the recursion coefficients.
	*/
	for ( i = 0; i < order; ++i )
		b[i] = ( alpha + ( ityp ) ( (i<<1) + 1 ) );
	
	for ( i = 0; i < order; ++i )
		c[i] = ( ityp ) ( i ) * ( alpha + ( ityp ) ( i ) );
	prod = 1.00;
	for ( i = 1; i < order; ++i )
		prod *= c[i];
	cc = r8_gamma ( alpha + 1.00 ) * prod;
	
	for ( i = 0; i < order; ++i )
	{
		/*
		Compute an estimate for the root.
		*/
		if ( i == 0 )
			x0 = ( 1.00 + alpha ) * ( 3.00+ 0.92 * alpha ) / ( 1.00 + 2.40 * ( ityp ) ( order ) + 1.80 * alpha );
		else if ( i == 1 )
			x0 = x0 + ( 15.0 + 6.25 * alpha ) / 	( 1.00 + 0.90 * alpha + 2.50 * ( ityp ) ( order ) );
		else
		{
			r1 = ( 1.00 + 2.55 * ( ityp ) ( i - 1 ) ) 	/ ( 1.90 * ( ityp ) ( i - 1 ) );
			r2 = 1.26 * ( ityp ) ( i - 1 ) * alpha / ( 1.00 + 3.50 * ( ityp ) ( i - 1 ) );
			ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );
			x0 += ratio * ( x0 - x[i-2] );
		}
		/*
		Use iteration to find the root.
		*/
		gen_laguerre_root ( &x0, order, alpha, &dp2, &p1, b, c );
		/*
		Set the abscissa and weight.
		*/
		x[i] = x0;
		w[i] = ( cc / dp2 ) / p1;
	}
	
	free ( b );
	free ( c );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gen_laguerre_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_LAGUERRE_COMPUTE_POINTS: Generalized Laguerre quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, double ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.
    Output, double X[ORDER], the abscissas.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	const register ityp alpha = s_data->a2; 
	
	
  ityp *w = ( ityp * ) malloc ( order * sizeof ( ityp ) );
  gen_laguerre_compute ( order, alpha, x, w );
  free ( w );
  return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _gen_laguerre_compute_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_LAGUERRE_COMPUTE_POINTS_NP: Generalized Laguerre quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	gen_laguerre_compute_points ( order, p[0], x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _laguerre_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_COMPUTE computes a Laguerre quadrature rule.
  Discussion:
    The integration interval is [ 0, +oo ).
    The weight function is w(x) = exp ( -x );.
    If the integral to approximate is:
        Integral ( 0 <= X < +oo ) exp ( - X ) * F(X) dX
    then the quadrature rule is:
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
    If the integral to approximate is:
        Integral ( A <= X < +oo ) F(X) dX
    then the quadrature rule is:
      Sum ( 1 <= I <= ORDER ) W(I) * exp ( X(I) ) * F ( X(I) )
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  Reference:
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  Parameters:
    Input, int ORDER, the order of the rule.
    1 <= ORDER.
    Output, double X[ORDER], the abscissas.
    Output, double W[ORDER], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
	ityp *b;
	ityp *c;
	ityp cc;
	ityp dp2;
	dim_typ i;
	ityp p1;
	ityp prod;
	ityp r1;
	ityp x0;
	
	if ( order == 0)
		return NULL;
		
	b = ( ityp * ) malloc ( order * sizeof ( double ) );
	c = ( ityp * ) malloc ( order * sizeof ( double ) );
	/*
	Set the recursion coefficients.
	*/
	for ( i = 0; i < order; ++i )
	{
		b[i] = ( ityp ) ( (i<<1) + 1 );
		c[i] = ( ityp ) ( i * i );
	}

	prod = 1.00;
	for ( i = 1; i < order; ++i )
		prod *= c[i];
	cc = prod;
	
	for ( i = 0; i < order; ++i )
	{
		/*
		Compute an estimate for the root.
		*/
		if ( i == 0 )
			x0 =  3.00 / ( 1.00 + 2.40 * ( ityp ) ( order ) );
		else if ( i == 1 )
			x0 += 15.00 / ( 1.00 + 2.50 * ( ityp ) ( order ) );
		else
		{
			r1 = ( 1.00 + 2.55 * ( ityp ) ( i - 1 ) ) / ( 1.90 * ( ityp ) ( i - 1 ) );
			x0 +=  r1 * ( x0 - x[i-2] );
		}
		/*
		Use iteration to find the root.
		*/
		laguerre_root ( &x0, order, &dp2, &p1, b, c );
		/*
		Set the abscissa and weight.
		*/
		x[i] = x0;
		w[i] = ( cc / dp2 ) / p1;
	}
	
	free ( b );
	free ( c );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _laguerre_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_COMPUTE_POINTS computes Laguerre quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Output, double X[ORDER], the abscissas.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	
	ityp *w = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	laguerre_compute ( order, x, w );
	free ( w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _laguerre_compute_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_COMPUTE_POINTS_NP computes Laguerre quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	laguerre_compute_points ( order, x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hermite_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_COMPUTE computes a Hermite quadrature rule.
  Discussion:
    The abscissas are the zeros of the N-th order Hermite polynomial.
    The integration interval is ( -oo, +oo ).
    The weight function is w(x) = exp ( - x * x ).
    The integral to approximate:
      Integral ( -oo < X < +oo ) exp ( - X * X ) * F(X) dX
    The quadrature rule:
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  Reference:
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  Parameters:
    Input, int ORDER, the order of the rule.
    1 <= ORDER.
    Output, double X[ORDER], the abscissas.
    Output, double W[ORDER], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
	ityp cc;
	ityp dp2;
	dim_typ i;
	ityp p1;
	ityp s;
	ityp temp;
	ityp x0;
	
	if ( order == 0)
		return NULL;
	
	cc = 1.7724538509 * r8_gamma ( ( ityp ) ( order ) ) / pow ( 2.00, order - 1 );
	
	s = pow ( 2.00 * ( ityp ) ( order ) + 1.00, 1.00 / 6.00 );
	
	for ( i = 0; i < ( order + 1 ) / 2; ++i )
	{
		if ( i == 0 )
			x0 = s * s * s - 1.85575 / s;
		else if ( i == 1 )
		x0 = x0 - 1.14 * pow ( ( double ) ( order ), 0.426 ) / x0;
		else if ( i == 2 )
			x0 = 1.86 * x0 - 0.86 * x[0];
		else if ( i == 3 )
			x0 = 1.91 * x0 - 0.91 * x[1];
		else
			x0 = 2.00 * x0 - x[i-2];

		hermite_root ( &x0, order, &dp2, &p1 );
		
		x[i] = x0;
		w[i] = ( cc / dp2 ) / p1;
		
		x[order-i-1] = -x0;
		w[order-i-1] = w[i];
	}
	/*
	Reverse the order of the abscissas.
	*/
	for ( i = 1; i <= order/2; ++i )
	{
		temp       = x[i-1];
		x[i-1]     = x[order-i];
		x[order-i] = temp;
	}
	
	if ( ( order % 2 ) == 1 )
		x[(order-1)/2] = 0.00;
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hermite_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_COMPUTE_POINTS computes Hermite quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Output, double X[ORDER], the abscissas.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	
  ityp *w = ( double * ) malloc ( order * sizeof ( ityp ) );
  hermite_compute ( order, x, w );
  free ( w );
  return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hermite_compute_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_COMPUTE_POINTS_NP computes Hermite quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	hermite_compute_points ( order, x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gen_hermite_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_HERMITE_COMPUTE computes a Generalized Hermite quadrature rule.
  Discussion:
    The integral to be approximated has the form:
      Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, int ORDER, the order of the rule.
    1 <= ORDER.
    Input, double ALPHA, the exponent of the X factor.
    -1.0 < ALPHA.
    Output, double X[ORDER], W[ORDER], the abscissas and weights of the rule.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register ityp alpha = s_data->a1; 
	ityp * x = s_data->a2;
	ityp * w = s_data->a3;
	
	ityp alpha_laguerre;
	ityp arg;
	dim_typ i;
	dim_typ order_laguerre;
	ityp *w_laguerre;
	ityp *x_laguerre;
	
	if ( order == 0)
		return NULL;
	
	if ( order == 1 )
	{
		arg = ( alpha + 1.00 ) / 2.00;
		x[0] = 0.00;
		w[0] = r8_gamma ( arg );
		return NULL;
	}
	
	if ( ( order % 2 ) == 0 ) 
	{
		order_laguerre = order / 2;
		alpha_laguerre = ( alpha - 1.00 ) / 2.00;
	}
	else
	{
		order_laguerre = ( order - 1 ) / 2;
		alpha_laguerre = ( alpha + 1.00 ) / 2.00;
	}
	
	w_laguerre = ( ityp * ) malloc ( order_laguerre * sizeof ( ityp ) );
	x_laguerre = ( ityp * ) malloc ( order_laguerre * sizeof ( ityp ) );
	
	gen_laguerre_compute ( order_laguerre, alpha_laguerre, x_laguerre, w_laguerre );
	
	if ( ( order % 2 ) == 0 )
	{
		for ( i = 0; i < order_laguerre; ++i )
			x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
		for ( i = 0; i < order_laguerre; ++i )
			x[order_laguerre+i] = sqrt ( x_laguerre[i] );
		for ( i = 0; i < order_laguerre; ++i )
			w[i] = 0.50 * w_laguerre[order_laguerre-1-i];
		for ( i = 0; i < order_laguerre; ++i )
			w[order_laguerre+i] = 0.50 * w_laguerre[i];
	}
	else if ( ( order % 2 ) == 1 )
	{
		for ( i = 0; i < order_laguerre; ++i)
			x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
		x[order_laguerre] = 0.00;
		for ( i = 0; i < order_laguerre; ++i )
			x[order_laguerre+1+i] = sqrt ( x_laguerre[i] );
		for ( i = 0; i < order_laguerre; ++i )
			w[i] = 0.50 * w_laguerre[order_laguerre-1-i] / x_laguerre[order_laguerre-1-i];
	
		arg = ( alpha + 1.00 ) / 2.00;
		w[order_laguerre] = r8_gamma ( arg );
		for ( i = 0; i < order_laguerre; ++i )
			w[order_laguerre] = w[order_laguerre] - w_laguerre[i] / x_laguerre[i];

	
		for ( i = 0; i < order_laguerre; ++i )
			w[order_laguerre+1+i] = 0.50 * w_laguerre[i] / x_laguerre[i];
	}
	free ( w_laguerre );
	free ( x_laguerre );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gen_hermite_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_HERMITE_COMPUTE_POINTS computes Generalized Hermite quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, double ALPHA, the exponent of the X factor.
    -1.0 < ALPHA.
    Output, double X[ORDER], the abscissas.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	const register ityp alpha = s_data->a2;
	
  ityp *w = ( ityp * ) malloc ( order * sizeof ( ityp ) );
  gen_hermite_compute ( order, alpha, x, w );
  free ( w );
  return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  inline  void *  _gen_hermite_compute_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_HERMITE_COMPUTE_POINTS_NP: Generalized Hermite quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], contains parameters.
    P[0] = ALPHA, the exponent of the X factor. -1.0 < ALPHA.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	gen_hermite_compute_points ( order, p[0], x );	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _legendre_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_COMPUTE computes a Legendre quadrature rule.
  Discussion:
    The integration interval is [ -1, 1 ].
    The weight function is w(x) = 1.0.
    The integral to approximate:
      Integral ( -1 <= X <= 1 ) F(X) dX
    The quadrature rule:
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt.
  Reference:
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  Parameters:
    Input, int ORDER, the order of the rule.
    1 <= ORDER.
    Output, double X[ORDER], the abscissas of the rule.
    Output, double W[ORDER], the weights of the rule.
    The weights are positive, symmetric, and should sum to 2.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
	ityp d1;
	ityp d2pn;
	ityp d3pn;
	ityp d4pn;
	ityp dp;
	ityp dpn;
	ityp e1;
	ityp fx;
	ityp h;
	dim_typ i;
	int iback;
	int k;
	int m;
	int mp1mi;
	int ncopy;
	int nmove;
	ityp p;
	ityp pk;
	ityp pkm1;
	ityp pkp1;
	ityp t;
	ityp u;
	ityp v;
	ityp x0;
	ityp xtemp;
	
	if ( order == 0)
		return NULL;
	e1 = ( ityp ) ( order * ( order + 1 ) );
	
	m = ( order + 1 ) / 2;
	
	for ( i = 1; i <= m; ++i )
	{
		mp1mi = m + 1 - i;
		t = ( ityp ) ( (i<<2) - 1 ) * M_PI / ( ityp ) ( (order<<2) + 2 );
		
		x0 =  cos ( t ) * ( 1.00 - ( 1.00 - 1.00 / ( ityp ) ( order ) ) / ( ityp ) ( (order<<3) * order ) );
		
		pkm1 = 1.00;
		pk = x0;
		
		for ( k = 2; k <= order; ++k )
		{
			pkp1 = 2.00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( ityp ) ( k );
			pkm1 = pk;
			pk = pkp1;
		}
	
		d1 = ( ityp ) ( order ) * ( pkm1 - x0 * pk );
		dpn = d1 / ( 1.00 - x0 * x0 );
		d2pn = ( 2.00 * x0 * dpn - e1 * pk ) / ( 1.00 - x0 * x0 );
		d3pn = ( 4.00 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );
		d4pn = ( 6.00 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );
		u = pk / dpn;
		v = d2pn / dpn;
		/*
		Initial approximation H:
		*/
		h = -u * ( 1.00 + 0.50 * u * ( v + u * ( v * v - d3pn / ( 3.00 * dpn ) ) ) );
		/*
		Refine H using one step of Newton's method:
		*/
		p = pk + h * ( dpn + 0.50 * h * ( d2pn + h / 3.00 
		* ( d3pn + 0.25 * h * d4pn ) ) );
		dp = dpn + h * ( d2pn + 0.50 * h * ( d3pn + h * d4pn / 3.00 ) );
		
		h += - p / dp;
		xtemp = x0 + h;
		x[mp1mi-1] = xtemp;
	
		fx = d1 - h * e1 * ( pk + 0.50 * h * ( dpn + h / 3.00 * ( d2pn + 0.25 * h * ( d3pn + 0.20 * h * d4pn ) ) ) );

		w[mp1mi-1] = 2.00 * ( 1.00 - xtemp * xtemp ) / ( fx * fx );
	}
	
	if ( ( order % 2 ) == 1 )
	x[0] = 0.00;
	/*
	Shift the data up.
	*/
	nmove = ( order + 1 ) / 2;
	ncopy = order - nmove;
	
	for ( i = 1; i <= nmove; ++i )
	{
		iback = order + 1 - i;
		x[iback-1] = x[iback-ncopy-1];
		w[iback-1] = w[iback-ncopy-1];
	}
	/*
	Reflect values for the negative abscissas.
	*/
	for ( i = 1; i <= order - nmove; ++i)
	{
		x[i-1] = - x[order-i];
		w[i-1] = w[order-i];
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _legendre_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_COMPUTE_POINTS computes Legendre quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Output, double X[ORDER], the abscissas.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	
	ityp *w = ( ityp * ) malloc ( order * sizeof ( ityp ) );
	legendre_compute ( order, x, w );
	free ( w );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _legendre_compute_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_COMPUTE_POINTS_NP computes Legendre quadrature points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	legendre_compute_points ( order, x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _clenshaw_curtis_compute_points ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLENSHAW_CURTIS_COMPUTE_POINTS computes Clenshaw Curtis quadrature points.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    This rule is defined on [-1,1].
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Output, double X[ORDER], the abscissas.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	ityp * x = s_data->a1;
	
	dim_typ index;
	
	if ( order == 0 )
		return NULL;
	else if ( order == 1 )
		x[0] = 0.00;
	else
	{
		for ( index = 1; index <= order; ++index )
			x[index-1] =  cos ( ( ityp ) ( order - index ) * M_PI / ( ityp ) ( order - 1     ) );
		x[0] = -1.00;
		if ( ( order % 2 ) == 1 )
			x[(order-1)/2] = 0.00;
		x[order-1] = +1.00;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _clenshaw_curtis_compute_points_np ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLENSHAW_CURTIS_COMPUTE_POINTS_NP: Clenshaw Curtis quadrature points.
  Discussion:
    Our convention is that the abscissas are numbered from left to right.
    This rule is defined on [-1,1].
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int ORDER, the order of the rule.
    Input, int NP, the number of parameters.
    Input, double P[NP], parameters which are not needed by this function.
    Output, double X[ORDER], the abscissas.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ order = s_data->a0; 
	const register dim_typ np = s_data->a1; 
	ityp * p = s_data->a2;
	ityp * x = s_data->a3;
	
	clenshaw_curtis_compute_points ( order, x );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sgmga_aniso_normalize ( void * data)
/******************************************************************************/
/*
  Purpose:
    SGMGA_ANISO_NORMALIZE normalizes the SGMGA anisotropic weight vector.
  Discussion:
    It is convenient for the user to initialize the anisotropic weight
    vector with any set of positive values.  These values are to be used
    as coefficients of the 1D levels, to evaluate an expression which
    determines which 1D levels will be included in a given rule.
    This means that a relatively LARGE coefficient forces the corresponding
    level to be relatively SMALL.  This is perhaps the opposite of what
    a user might expect.  If a user wishes to use an importance vector,
    so that a relatively large importance should correspond to more levels,
    and hence more points, in that dimension, then the function
    SGMGA_IMPORTANCE_TO_ANISO should be called first!
    Since the weights only represent the relative importance of the
    components, they may be multiplied by any (positive) scale factor.
    Nonetheless, it may be convenient ot choose a particular normalization
    for the weights.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 November 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int OPTION, the normalization option.
    0, no scaling is applied.
    1, the weights are scaled so that the minimum nonzero entry is 1.
    2, the weights are scaled so that they sum to DIM_NUM.
    Input, int DIM_NUM, the spatial dimension.
    Input/output, double LEVEL_WEIGHT[DIM_NUM], the anisotropic
    weights.  The input values must be strictly positive.
    On output, these have been normalized.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ option = s_data->a0; 
	const register dim_typ dim_num = s_data->a1; 
	ityp * level_weight = s_data->a2;
	
    dim_typ dim;
    bool found;
    int level_weight_min;
    int level_weight_sum;
    /*
    Option 0, no normalization.
    */
    if ( option == 0 );
    /*
    Option 1, the minimum nonzero entry is 1.
    */
    else if ( option == 1 )
    {
        level_weight_min = i4_huge;
        found = 0;
        for ( dim = 0; dim < dim_num; ++dim )
            if ( 0.00 < level_weight[dim] )
                if ( level_weight[dim] < level_weight_min )
                {
                    level_weight_min = level_weight[dim];
                    ++ found;
                }

        if ( found == 0 )
            return NULL;

        for ( dim = 0; dim < dim_num; ++dim )
            level_weight[dim] /= level_weight_min;
    }
    /*
    Option 2, rescale so sum of weights is DIM_NUM.
    */
    else if ( option == 2 )
    {
        level_weight_sum = r8vec_sum ( dim_num, level_weight );

        if ( level_weight_sum <= 0.00 )
            return NULL;
        for ( dim = 0; dim < dim_num; ++dim )
            level_weight[dim] = ( ( ityp ) ( dim_num ) * level_weight[dim] )/ level_weight_sum;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sgmga_importance_to_aniso ( void * data)
/******************************************************************************/
/*
  Purpose:
    SGMGA_IMPORTANCE_TO_ANISO: importance vector to anisotropic weight vector.
  Discussion:
    To specify the anisotropy of a multidimensional problem, the user is
    allowed to specify an "importance vector".  This vector can contain
    any set of positive values.  These values represent the relative
    importance of each dimension.  These values, with a suitable normalization,
    will be used to evaluate a constraint of the following form:
      QMIN < Level(1) / Importance(1) + Level(2) / Importance(2) + ...
             Level(N) / Importance(N) <= QMAX
    and a set of levels that satisfies this constraint will then be included
    in a given anistotropic sparse grid rule.  Thus, increasing the
    importance value of a particular dimension allows larger level values
    in that dimension to satisfy the constraint.
    The program actually works with coefficients LEVEL_WEIGHT that are
    the inverse of the importance vector entries, with a suitable
    normalization.  This function is supplied to convert between the
    more natural "importance vector" and the internally useful
    "level_weight" vector.
    This function converts the importance vector to an unnormalized
    anisotropy weight vector.
    Note that some (but not all) of the IMPORTANCE vector entries may be zero.
    This indicates that the corresponding dimension is of "zero" or
    rather "minimal" importance.  In such a case, only a one-point quadrature
    rule will be applied for that dimension, no matter what sparse grid
    level is requested for the overall problem.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 November 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, double IMPORTANCE[DIM_NUM], the importance vector.
    All entries must be nonnegative, and at least one must be positive.
    Output, double LEVEL_WEIGHT[DIM_NUM], the anisotropic
    weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ dim_num = s_data->a0; 
	ityp * importance = s_data->a1;
	ityp * level_weight = s_data->a2;
	
    dim_typ dim;
    bool found;
    ityp level_weight_norm;

    for ( dim = 0; dim < dim_num; ++dim)
        if ( importance[dim] < 0.00 )
            return NULL;

    found = 0;

    for ( dim = 0; dim < dim_num; ++dim)
    {
        if ( 0.0 < importance[dim] )
        {
            level_weight[dim] = 1.00 / importance[dim];
            ++ found;
        }
        else
            level_weight[dim] = 0.0;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   sgmga_index ( const register dim_typ dim_num, ityp level_weight[], const register dim_typ level_max,
  int rule[static dim_num], const register dim_typ point_num, const register dim_typ point_total_num, int sparse_unique_index[static point_total_num],
  void level_to_order ( dim_typ dim_num, int level[], int rule[], int order[] ),int sparse_order[static dim_num*point_num], int sparse_index[static dim_num*point_num] )
/******************************************************************************/
/*
  Purpose:
    SGMGA_INDEX indexes an SGMGA grid.
  Discussion:
    For each "unique" point in the sparse grid, we return its INDEX and ORDER.
    That is, for the I-th unique point P, we determine the product grid which
    first generated this point, and  and we return in SPARSE_ORDER the orders
    of the 1D rules in that grid, and  and in SPARSE_INDEX the component
    indexes in those rules that generated this specific point.
    For instance, say P was first generated by a rule which was a 3D product
    of a 9th order CC rule and  and a 15th order GL rule, and  and that to
    generate P, we used the 7-th point of the CC rule and  and the 3rh point
    of the GL rule.  Then the SPARSE_ORDER information would be (9,15) and
    the SPARSE_INDEX information would be (7,3).  This, combined with the
    information in RULE, is enough to regenerate the value of P.
    The user must preallocate space for the output arrays SPARSE_ORDER and
    SPARSE_INDEX.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 November 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
    Input, int POINT_NUM, the number of unique points
    in the grid.
    Input, int POINT_TOTAL_NUM, the total number of points in the grid.
    Input, int SPARSE_UNIQUE_INDEX[POINT_TOTAL_NUM], associates each
    point in the grid with its unique representative.
    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[],
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or
    "level_to_order_linear".
    Output, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists,
    for each point, the order of the 1D rules used in the grid that
    generated it.
    Output, int SPARSE_INDEX[DIM_NUM*POINT_NUM)] lists, for
    each point, its index in each of the 1D rules in the grid that generated
    it.  The indices are 1-based.
*/
{
    ityp coef;
    dim_typ dim;
    int *level_1d;
    int *level_1d_max;
    ityp level_weight_min_pos;
    bool more_grids;
    bool more_points;
    int *order_1d;
    dim_typ point;
    dim_typ point_count;
    int *point_index;
    dim_typ point_unique;
    ityp q_max;
    ityp q_min;
    /*
    Special cases.
    */

    if ( level_max == 0 )
    {
        point = 0;
        for ( dim = 0; dim < dim_num; ++dim )
            sparse_order[dim+point*dim_num] = sparse_index[dim+point*dim_num] = 1;
        return;
    }
    /*
    Initialize the INDEX and ORDER arrays to -1 to help catch errors.
    */
    for ( point = 0; point < point_num; ++point)
        for ( dim = 0; dim < dim_num; ++dim )
            sparse_order[dim+point*dim_num] = sparse_index[dim+point*dim_num] = -1;

    point_count = 0;

    level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
    order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    point_index = ( int * ) malloc ( dim_num * sizeof ( int ) );
    /*
    Initialization for SGMGA_VCN_ORDERED.
    */
    level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
    q_min = ( ityp ) ( level_max ) * level_weight_min_pos- r8vec_sum ( dim_num, level_weight );
    q_max = ( ityp ) ( level_max ) * level_weight_min_pos;
    for ( dim = 0; dim < dim_num; ++dim )
        level_1d_max[dim] = (0.00 < level_weight[dim])*floor ( q_max / level_weight[dim] ) + 1;
    more_grids = false;
    /*
    Repeatedly call SGMGA_VCN_ORDERED, seeking all vectors LEVEL_1D
    which satisfy the constraint:

    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT )
    < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
    <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
    */
    for ( ; ; )
    {
        sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d, q_min, q_max, &more_grids );

        if ( !more_grids )
            break;
        /*
        Compute the combinatorial coefficient.
        */
        coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, level_1d,
        q_min, q_max );

        if ( coef == 0.00 )
            continue;
        /*
        Transform each 1D level to a corresponding 1D order.
        */
        level_to_order ( dim_num, level_1d, rule, order_1d );
        /*
        The inner loop generates a POINT of the GRID of the LEVEL.
        */
        more_points = false;

        for ( ; ; )
        {
            vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

            if ( !more_points )
                break;
            point_unique = sparse_unique_index[point_count];
            for ( dim = 0; dim < dim_num; ++dim )
            {
                sparse_order[dim+point_unique*dim_num] = order_1d[dim];
                sparse_index[dim+point_unique*dim_num] = point_index[dim];
            }
            ++ point_count;
        }
    }

    free ( level_1d );
    free ( level_1d_max );
    free ( order_1d );
    free ( point_index );

    return;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   sgmga_point ( const register dim_typ dim_num, ityp level_weight[static dim_num], const register dim_typ level_max,
  int rule[static dim_num], int np[static dim_num], ityp p[],
  void ( *gw_compute_points[] ) ( int order, int np, ityp p[], ityp x[] ),
  int point_num, int sparse_order[static dim_num*point_num], int sparse_index[static dim_num*point_num],void level_to_order ( int dim_num, int level[], int rule[], int order[] ),ityp sparse_point[static dim_num*point_num] )
/******************************************************************************/
/*
  Purpose:
    SGMGA_POINT computes the points of an SGMGA rule.
  Discussion:
    The sparse grid is the logical sum of low degree product rules.
    Each product rule is the product of 1D factor rules.
    The user specifies:
    * the spatial dimension of the quadrature region,
    * the level that defines the Smolyak grid.
    * the quadrature rules.
    * the number of points.
    The user must preallocate space for the output array SPARSE_POINT.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 December 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
    Input, int LEVEL_MAX, controls the size of the final
    sparse grid.
    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
    Input, int NP[DIM_NUM], the number of parameters used by each rule.
    Input, double P[sum(NP[*])], the parameters needed by each rule.
    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
    an array of pointers to functions which return the 1D quadrature points
    associated with each spatial dimension for which a Golub Welsch rule
    is used.
    Input, int POINT_NUM, the number of points in the grid,
    as determined by SGMGA_SIZE.
    Input, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists, for each point,
    the order of the 1D rules used in the grid that generated it.
    Input, int SPARSE_INDEX[DIM_NUM*POINT_NUM], lists, for each point,
    its index in each of the 1D rules in the grid that generated it.
    The indices are 1-based.
    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[],
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or
    "level_to_order_linear".
    Output, double SPARSE_POINT[DIM_NUM*POINT_NUM], the points.
*/
{
    dim_typ dim;
    int level;
    int *level_1d_max;
    ityp level_weight_min_pos;
    int order;
    int p_index;
    dim_typ point;
    ityp *points;
    ityp q_max;

    for ( point = 0; point < point_num; ++point )
        for ( dim = 0; dim < dim_num; ++dim)
            sparse_point[dim+point*dim_num] = - r8_huge;
    /*
    Compute the point coordinates.
    */
    level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
    level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
    q_max = ( ityp ) ( level_max ) * level_weight_min_pos;

    p_index = 0;

    for ( dim = 0; dim < dim_num; ++dim )
    {
        level_1d_max[dim] = 0 + (0.00 < level_weight[dim])*floor ( q_max / level_weight[dim] ) + 1;

        for ( level = 0; level <= level_1d_max[dim]; ++level )
        {
            level_to_order ( 1, &level, rule+dim, &order );
            points = ( ityp * ) malloc ( order * sizeof ( ityp ) );

            if ( rule[dim] == 1 )
                clenshaw_curtis_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 2 )
                fejer2_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 3 )
                patterson_lookup_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 4 )
                legendre_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 5 )
                hermite_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 6 )
                gen_hermite_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 7 )
                laguerre_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 8 )
                gen_laguerre_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 9 )
                jacobi_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 10 )
                gw_compute_points[dim] (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 11 )
                clenshaw_curtis_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 12 )
                fejer2_compute_points_np (order, np[dim], p+p_index, points );
            else if ( rule[dim] == 13 )
                patterson_lookup_points_np (order, np[dim], p+p_index, points );
            else
                return;

            for ( point = 0; point < point_num; ++point )
            {
                if ( sparse_order[dim+point*dim_num] == order )
                {
                    sparse_point[dim+point*dim_num] =
                    points[sparse_index[dim+point*dim_num]-1];
                }
            }
            free ( points );
        }
        p_index += np[dim];
    }
    /*
    Check to see if we missed any points.
    */
    for ( point = 0; point < point_num; ++point)
        for ( dim = 0; dim < dim_num; ++dim )
            if ( sparse_point[dim+point*dim_num] == -  r8_huge )
                return;

    free ( level_1d_max );

    return;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   sgmga_product_weight ( const register dim_typ dim_num, int order_1d[static dim_num], const register dim_typ order_nd,
  int rule[], int np[static dim_num], ityp p[],void ( *gw_compute_weights[] ) ( dim_typ order, dim_typ np, ityp p[], ityp w[] ), ityp weight_nd[static order_nd] )
/******************************************************************************/
/*
  Purpose:
    SGMGA_PRODUCT_WEIGHT computes the weights of a mixed product rule.
  Discussion:
    This routine computes the weights for a quadrature rule which is
    a product of 1D rules of varying order and kind.
    The user must preallocate space for the output array WEIGHT_ND.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 December 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
    Input, int ORDER_ND, the order of the product rule.
    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
    Input, int NP[DIM_NUM], the number of parameters used by each rule.
    Input, double P[sum(NP[*])], the parameters needed by each rule.
    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
    an array of pointers to functions which return the 1D quadrature weights
    associated with each spatial dimension for which a Golub Welsch rule
    is used.
    Output, double WEIGHT_ND[ORDER_ND], the product rule weights.
*/
{
    dim_typ dim;
    dim_typ i;
    dim_typ p_index;
    ityp *weight_1d;

    for ( i = 0; i < order_nd; ++i )
        weight_nd[i] = 1.00;

    p_index = 0;

    for ( dim = 0; dim < dim_num; ++dim)
    {
        weight_1d = ( ityp * ) malloc ( order_1d[dim] * sizeof ( ityp ) );

        if ( rule[dim] == 1 )
            clenshaw_curtis_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 2 )
            fejer2_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 3 )
            patterson_lookup_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 4 )
            legendre_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 5 )
            hermite_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 6 )
            gen_hermite_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 7 )
            laguerre_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 8 )
            gen_laguerre_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 9 )
            jacobi_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 10 )
            gw_compute_weights[dim] (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 11 )
            clenshaw_curtis_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 12 )
            fejer2_compute_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else if ( rule[dim] == 13 )
            patterson_lookup_weights_np (order_1d[dim], np[dim], p+p_index, weight_1d );
        else
            return;

        p_index += np[dim];

        r8vec_direct_product2 ( dim, order_1d[dim], weight_1d,
        dim_num, order_nd, weight_nd );

        free ( weight_1d );
    }
    return;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  dim_typ   sgmga_size_total ( const register dim_typ dim_num, ityp level_weight[static dim_num], const register dim_typ level_max,
  int rule[static dim_num],void level_to_order ( dim_typ dim_num, int level[], int rule[], int order[] ) )
/******************************************************************************/
/*
  Purpose:
    SGMGA_SIZE_TOTAL sizes an SGMGA grid, counting duplicates.
  Discussion:
    This routine returns the total point count for an SGMGA
 ( Sparse Grid of Mixed type with Growth rule and Anisotropic weights).
    The sparse grid is the logical sum of product grids.
    The sparse grid has an associated integer index LEVEL_MAX, whose lowest
    value is 0.  LEVEL_MAX = 0 indicates the sparse grid made up of one product
    grid, which in turn is the product of 1D factor grids of the lowest level.
    This usually means the sparse grid with LEVEL_MAX equal to 0 is a
    one point grid.
    We can assign a level to each factor grid, and hence a LEVEL vector
    to the corresponding product grid, and a weighted index
    LEVEL_GRID (which will in general be a real number):
      LEVEL_GRID = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL(I)
    The product grid will participate in the formation of the sparse grid
    if it satisfies the following weighted constraint:
      LEVEL_MAX - DIM_NUM < LEVEL_GRID <= LEVEL_MAX
    This routine determines the total number of abscissas in all the
    product rules used to form the SGMGA associated with the index LEVEL_MAX.
    The count disregards duplication.  If the same multidimensional abcsissa
    occurs in two different product rules that are part of the SGMGA, then
    that single abcissa is counted twice.
    This computation is useful in cases where the entire set of abscissas
    is going to be generated, preparatory to compression to finding, indexing
    and merging the duplicate abcissass.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 November 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[],
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or
    "level_to_order_linear".
    Output, int SGMGA_SIZE_TOTAL, the number of points
    including repetitions.
*/
{
    ityp coef;
    dim_typ dim;
    int *level_1d;
    int *level_1d_max;
    ityp level_weight_min_pos;
    bool more_grids;
    int *order_1d;
    dim_typ point_total_num;
    ityp q_max;
    ityp q_min;
    /*
    Special case.
    */
    if ( level_max == 0 )
        return 1;

    point_total_num = 0;

    level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
    order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    /*
    Initialization for SGMGA_VCN_ORDERED.
    */
    level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
    q_min = ( ityp ) ( level_max ) * level_weight_min_pos- r8vec_sum ( dim_num, level_weight );
    q_max = ( ityp ) ( level_max ) * level_weight_min_pos;
    for ( dim = 0; dim < dim_num; ++dim )
        level_1d_max[dim] = 0 + ( 0.00 < level_weight[dim] )*floor ( q_max / level_weight[dim] ) + 1;

    more_grids = false;
    /*
    Repeatedly call SGMGA_VCN_ORDERED, seeking all vectors LEVEL_1D
    which satisfy the constraint:

    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT )
    < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
    <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
    */
    for ( ; ; )
    {
        sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d,
        q_min, q_max, &more_grids );

        if ( !more_grids )
            break;
        /*
        Compute the combinatorlal coefficient.
        */
        coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, level_1d,
        q_min, q_max );

        if ( coef == 0.00 )
            continue;
        /*
        Transform each 1D level to a corresponding 1D order.
        */
        level_to_order ( dim_num, level_1d, rule, order_1d );
        point_total_num = point_total_num + i4vec_product ( dim_num,order_1d );
    }
    free ( level_1d );
    free ( level_1d_max );
    free ( order_1d );

    return point_total_num;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sgmga_vcn ( void * data)
/******************************************************************************/
/*
  Purpose:
    SGMGA_VCN returns the next constrained vector.
  Discussion:
    We consider vectors X of dimension DIM_NUM satisfying:
      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
    and define
      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
    and seek X satisfying the constraint:
      Q_MIN < Q <= Q_MAX
    For sparse grid applications, we compute
      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
    and assume there is an underlying LEVEL used to index the sets of
    constrained vectors, and that
      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
    This routine returns, one at a time exactly those X which satisfy
    the constraint.  No attempt is made to return the X values in
    any particular order as far as Q goes.
  Example:
    LEVEL_WEIGHT:          1.000000        1.000000
    Q_MIN:        0.000000
    Q_MAX:        2.000000
    X_MAX:                         2         2
         1        1.000000         1         0
         2        2.000000         2         0
         3        1.000000         0         1
         4        2.000000         1         1
         5        2.000000         0         2
    LEVEL_WEIGHT:          1.000000        2.000000
    Q_MIN:       -1.000000
    Q_MAX:        2.000000
    X_MAX:                         2         1
         1        0.000000         0         0
         2        1.000000         1         0
         3        2.000000         2         0
         4        2.000000         0         1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 November 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the number of components in the vector.
    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Input, double Q_MIN, Q_MAX, the lower and upper
    limits on the sum.
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dtpit2pi2itpb * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	ityp * level_weight = s_data->a1;
	int * x_max = s_data->a2;
	int * x = s_data->a3;
	const register ityp q_min = s_data->a4;
	const register ityp q_max = s_data->a5;
	bool * more = s_data->a6;
	
    dim_typ i, j;
    ityp q;

    if ( ! ( *more ) )
    {
        *more = true;
        for ( i = 0; i < dim_num; ++i )
            x[i] = 0;

        q = 0.00;
        for ( i = 0; i < dim_num; ++i )
            q += level_weight[i] * ( ityp ) ( x[i] );

        if ( q_min < q && q <= q_max )
            return NULL;
    }

    for ( ; ; )
    {
        j = 0;

        for ( ; ; )
        {
            if ( x[j] < x_max[j] )
                break;

            if ( dim_num - 1 <= j )
            {
                *more = false;
                return NULL;
            }
            ++ j;
        }

        ++ x[j];
        for ( i = 0; i < j; ++i )
            x[i] = 0;

        q = 0.00;
        for ( i = 0; i < dim_num; ++i )
            q += level_weight[i] * ( ityp ) ( x[i] );

        if ( q_min < q && q <= q_max )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sgmga_vcn_coef ( void * data)
/******************************************************************************/
/*
  Purpose:
    SGMGA_VCN_COEF returns the "next" constrained vector's coefficient.
  Discussion:
    We are given a vector X of dimension DIM_NUM which satisfies:
      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
    and the following constraint:
      Q_MIN < sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
    This routine computes the appropriate coefficient for X in the
    anisotropic sparse grid scheme.
    The coefficient is calculated as follows:
      Let B be a binary vector of length DIM_NUM, and let ||B|| represent
      the sum of the entries of B.
      Coef = sum ( all B such that X+B satisfies constraints ) (-1)^||B||
    Since X+0 satisfies the constraint, there is always at least one
    summand.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 November 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the number of components in the vector.
    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
    Input, int X_MAX[DIM_NUM], the maximum
    values allowed in each component.
    Input, int X[DIM_NUM], a point which satisifies the constraints.
    Input, double Q_MIN, Q_MAX, the lower and upper
    limits on the sum.
    Output, double SGMGA_VCN_COEF, the combinatorial coefficient.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit2pi2it * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	ityp * level_weight = s_data->a1;
	int * x_max = s_data->a2;
	int * x = s_data->a3;
	const register ityp q_min = s_data->a4;
	const register ityp q_max = s_data->a5;
	
    int *b;
    int b_sum;
    ityp coef;
    dim_typ i;
    bool legal;
    ityp q;
    int *x2;

    b = ( int * ) malloc ( dim_num * sizeof ( int ) );
    x2 = ( int * ) malloc ( dim_num * sizeof ( int ) );

    for ( i = 0; i < dim_num; ++i )
        b[i] = 0;
    coef = 1.00;

    for ( ; ; )
    {
        /*
        Generate the next binary perturbation.
        */
        binary_vector_next ( dim_num, b );
        b_sum = i4vec_sum ( dim_num, b );
        /*
        We're done if we've got back to 0.
        */
        if ( b_sum == 0 )
            break;
        /*
        Perturb the vector.
        */
        for ( i = 0; i < dim_num; ++i )
            x2[i] = x[i] + b[i];
        /*
        Does it satisfy the XMAX constraint?
        */
        legal = true;
        for ( i = 0; i < dim_num; ++i )
        {
            if ( x_max[i] < x2[i] )
            {
                legal = false;
                break;
            }
        }
        if ( !legal )
            continue;
        /*
        Does it satisfy the Q_MIN, Q_MAX constraint?
        */
        q = 0.00;
        for ( i = 0; i < dim_num; ++i )
            q += level_weight[i] * ( ityp ) ( x2[i] );

        if ( q_min < q && q <= q_max )
            coef += r8_mop ( b_sum );
    }

    free ( b );
    free ( x2 );

	result = coef;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _sgmga_vcn_ordered ( void * data)
/******************************************************************************/
/*
  Purpose:
    SGMGA_VCN_ORDERED returns the "next" constrained vector, with ordering.
  Discussion:
    We consider vectors X of dimension DIM_NUM satisfying:
      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
    and define
      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
    and seek X's satisfying the constraint:
      Q_MIN < Q <= Q_MAX
    For sparse grid applications, we compute
      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
    and assume there is an underlying LEVEL used to index the sets of
    constrained vectors, and that
      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
    This function returns, one at a time exactly those X which satisfy
    the constraint.
    A weak ordering is imposed on the solution vectors.  This function
    subdivides the range Q_MIN through Q_MAX into subintervals of width 1, so
    that the X vectors returned are roughly sorted (or at least binned)
    by Q value.
  Example:
    If the weights are also integral, then the X vectors are in fact SORTED
    by Q value:
    LEVEL_WEIGHT:          1.000000        1.000000
    Q_MIN:        0.000000
    Q_MAX:        2.000000
    X_MAX:                         2         2
         1        1.000000         1         0
         2        1.000000         0         1
         3        2.000000         2         0
         4        2.000000         1         1
         5        2.000000         0         2
    When the weights are not integral, then the X values are only BINNED
    by Q value, that is, we first get all X's with Q values between Q_MIN
    and Q_MIN+1, then Q_MIN+1 to Q_MIN+2 and so on, as demonstrated here:
    LEVEL_WEIGHT:             1.5               1
    Q_MIN:  0.5
    Q_MAX:  3
    X_MAX:                           2         3
           1             1.5         1         0
           2               1         0         1
           3             2.5         1         1
           4               2         0         2
           5               3         2         0
           6               3         0         3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 November 2009
  Author:
    John Burkardt
  Reference:
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the number of components in the vector.
    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.
    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.
    Input, double Q_MIN, Q_MAX, the lower and upper
    limits on the sum.
    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
	const dtpit2pi2itpb * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	ityp * level_weight = s_data->a1;
	int * x_max = s_data->a2;
	int * x = s_data->a3;
	const register ityp q_min = s_data->a4;
	const register ityp q_max = s_data->a5;
	bool * more = s_data->a6;
	
    ityp q;
    static ityp q_max2;
    static ityp q_min2;
    /*
    On first call, initialize the subrange.
    */
    if ( !(*more) )
    {
        q_min2 = q_min;
        q_max2 = MIN ( q_min + 1.00, q_max );
    }
    /*
    Call a lower level function to search the subrange.
    */
    for ( ; ; )
    {
        sgmga_vcn ( dim_num, level_weight, x_max, x, q_min2, q_max2, more );
        /*
        If another solution was found, return it.
        */
        if ( *more )
            return NULL;
        /*
        If the current subrange is exhausted, try to move to the next one.
        */
        if ( q_max2 < q_max )
        {
            q_min2 = q_max2;
            q_max2 = MIN ( q_max2 + 1.00, q_max );
        }
        /*
        If there are no more subranges, we're done.
        */
        else
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   sgmga_weight ( const register dim_typ dim_num, ityp level_weight[static dim_num], const register dim_typ level_max,
  int rule[static dim_num], int np[static dim_num], ityp p[],void ( *gw_compute_weights[] ) ( dim_typ order, dim_typ np, ityp p[], ityp w[] ),const register dim_typ point_num, const register dim_typ point_total_num, int sparse_unique_index[static point_total_num],void level_to_order ( dim_typ dim_num, int level[], int rule[], int order[] ),ityp sparse_weight[static point_num] )
/******************************************************************************/
/*
  Purpose:
    SGMGA_WEIGHT computes weights for an SGMGA grid.
  Discussion:
    The user must preallocate space for the output array SPARSE_WEIGHT.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 November 2009
  Author:
    John Burkardt
  Reference
    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
    Input, int LEVEL_MAX, the maximum value of LEVEL.
    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.
    Input, int NP[DIM_NUM], the number of parameters used by each rule.
    Input, double P[sum(NP[*])], the parameters needed by each rule.
    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
    an array of pointers to functions which return the 1D quadrature weights
    associated with each spatial dimension for which a Golub Welsch rule
    is used.
    Input, int POINT_NUM, the number of unique points
    in the grid.
    Input, int POINT_TOTAL_NUM, the total number of points
    in the grid.
    Input, int SPARSE UNIQUE_INDEX[POINT_TOTAL_NUM], lists,
    for each (nonunique) point, the corresponding index of the same point in
    the unique listing.
    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[],
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or
    "level_to_order_linear".
    Output, double SPARSE_WEIGHT[POINT_NUM], the weights
    associated with the sparse grid points.
*/
{
    ityp coef;
    dim_typ dim;
    ityp *grid_weight;
    dim_typ level;
    int *level_1d;
    int *level_1d_max;
    ityp level_weight_min_pos;
    bool more_grids;
    dim_typ order;
    int *order_1d;
    dim_typ order_nd;
    dim_typ point;
    dim_typ point_total;
    int point_unique;
    ityp q_max;
    ityp q_min;

    for ( point = 0; point < point_num; ++point )
        sparse_weight[point] = 0.00;

    point_total = 0;

    level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
    level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
    /*
    Initialization for SGMGA_VCN_ORDERED.
    */
    level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
    q_min = ( ityp ) ( level_max ) * level_weight_min_pos- r8vec_sum ( dim_num, level_weight );
    q_max = ( ityp ) ( level_max ) * level_weight_min_pos;
    for ( dim = 0; dim < dim_num; ++dim )
        level_1d_max[dim] = 0 + ( 0.00 < level_weight[dim] )*floor ( q_max / level_weight[dim] ) + 1;
    more_grids = false;
    /*
    Repeatedly call SGMGA_VCN_ORDERED, seeking all vectors LEVEL_1D
    which satisfy the constraint:

    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT )
    < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
    <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
    */
    for ( ; ; )
    {
        sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d,q_min, q_max, &more_grids );

        if ( !more_grids )
            break;
        /*
        Compute the combinatorial coefficient.
        */
        coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, level_1d,q_min, q_max );

        if ( coef == 0.00 )
            continue;
        /*
        Transform each 1D level to a corresponding 1D order.
        */
        level_to_order ( dim_num, level_1d, rule, order_1d );
        /*
        The product of the 1D orders gives us the number of points in this grid.
        */
        order_nd = i4vec_product ( dim_num, order_1d );
        /*
        Compute the weights for this grid.

        The correct transfer of data from the product grid to the sparse grid
        depends on the fact that the product rule weights are stored under colex
        order of the points, and this is the same ordering implicitly used in
        generating the SPARSE_UNIQUE_INDEX array.
        */
        grid_weight = ( ityp * ) malloc ( order_nd * sizeof ( ityp ) );
        sgmga_product_weight ( dim_num, order_1d, order_nd, rule,np, p, gw_compute_weights, grid_weight );
        /*
        Add these weights to the rule.
        */
        for ( order = 0; order < order_nd; ++order)
        {
            point_unique = sparse_unique_index[point_total];
            ++ point_total;
            sparse_weight[point_unique] = sparse_weight[point_unique]+ coef * grid_weight[order];
        }

        free ( grid_weight );
    }

    free ( level_1d );
    free ( level_1d_max );
    free ( order_1d );

    return;
}

#endif
