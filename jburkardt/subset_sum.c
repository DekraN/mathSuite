#ifndef __DISABLEDEEP_SUBSETSUM

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _gen_laguerre_recur ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_LAGUERRE_RECUR evaluates a Generalized Laguerre polynomial.
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
    Output, double *P2, the value of L(ORDER)(X).
    Output, double *DP2, the value of L'(ORDER)(X).
    Output, double *P1, the value of L(ORDER-1)(X).
    Input, double X, the point at which polynomials are evaluated.
    Input, int ORDER, the order of the polynomial to be computed.
    Input, double ALPHA, the exponent of the X factor in the
    integrand.
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
	const pitdt2it4pit * const s_data = data;
	
	ityp * p2 = s_data->a0;
	const register dim_typ order = s_data->a1;
	const register ityp x = s_data->a2;
	const register ityp alpha = s_data->a3;
	ityp * dp2 = s_data->a4;
	ityp * p1 = s_data->a5;
	ityp * b = s_data->a6;
	ityp * c = s_data->a7;
	
	ityp dp0;
	ityp dp1;
	dim_typ i;
	ityp p0;
	
	*p1 = *dp2 = 1.00;
	dp1 = 0.00;
	*p2 = x - alpha - 1.00;
	
	for ( i = 1; i < order; ++i )
	{
		p0 = *p1;
		dp0 = dp1;
		
		*p1 = *p2;
		dp1 = *dp2;
		
		*p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
		*dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hermite_recur ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
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
    Output, double *P2, the value of H(ORDER)(X).
    Output, double *DP2, the value of H'(ORDER)(X).
    Output, double *P1, the value of H(ORDER-1)(X).
    Input, double X, the point at which polynomials are evaluated.
    Input, int ORDER, the order of the polynomial to be computed.
*/
{
	const dtit3pit * const s_data = data;
	
	const register dim_typ order = s_data->a0;
	const register ityp x = s_data->a1;
	ityp * p2 = s_data->a2;
	ityp * dp2 = s_data->a3;
	ityp * p1 = s_data->a4;
	
	
	dim_typ i;
	ityp dq0;
	ityp dq1;
	ityp dq2;
	ityp q0;
	ityp q1;
	ityp q2;
	
	q1 = dq2 = 1.00;
	dq1 = 0.00;
	q2 = x;
	
	for ( i = 2; i <= order; ++i )
	{
		q0 = q1;
		dq0 = dq1;
		
		q1 = q2;
		dq1 = dq2;
		
		q2  = x * q1 - 0.50 * ( ( ityp ) ( i ) - 1.00 ) * q0;
		dq2 = x * dq1 + q1 - 0.50 * ( ( ityp ) ( i ) - 1.00 ) * dq0;
	}
	
	*p2 = q2;
	*dp2 = dq2;
	*p1 = q1;
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _jacobi_recur ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_RECUR evaluates a Jacobi polynomial.
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
    Output, double *P2, the value of J(ORDER)(X).
    Output, double *DP2, the value of J'(ORDER)(X).
    Output, double *P1, the value of J(ORDER-1)(X).
    Input, double X, the point at which polynomials are evaluated.
    Input, int ORDER, the order of the polynomial to be computed.
    Input, double ALPHA, BETA, the exponents of (1-X) and
 (1+X) in the quadrature rule.
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
	const _2itdt2pitit3pit * const s_data = data;
	
	const register ityp alpha = s_data->a0;
	const register ityp beta = s_data->a1;
	const register dim_typ order = s_data->a2;
	ityp * p2 = s_data->a3;
	ityp * dp2 = s_data->a4;
	const register ityp x = s_data->a5;
	ityp * p1 = s_data->a6;
	ityp * b = s_data->a7;
	ityp * c = s_data->a8;
	
	ityp dp0;
	ityp dp1;
	dim_typ i;
	ityp p0;
	
	*p1 = *dp2 = 1.00;
	dp1 = 0.00;
	
	*p2 = x + ( alpha - beta ) / ( alpha + beta + 2.00 );
	
	for ( i = 2; i <= order; ++i)
	{
		p0 = *p1;
		dp0 = dp1;
		
		*p1 = *p2;
		dp1 = *dp2;
		
		*p2 = ( x - b[i-1] ) * ( *p1 ) - c[i-1] * p0;
		*dp2 = ( x - b[i-1] ) * dp1 + ( *p1 ) - c[i-1] * dp0;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laguerre_recur ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_RECUR evaluates a Laguerre polynomial.
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
    Output, double *P2, the value of L(ORDER)(X).
    Output, double *DP2, the value of L'(ORDER)(X).
    Output, double *P1, the value of L(ORDER-1)(X).
    Input, double X, the point at which polynomials are evaluated.
    Input, int ORDER, the order of the polynomial to be computed.
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
	const itpitdt4pit * const s_data = data;
	
	const register ityp x = s_data->a0;
	ityp * p2 = s_data->a1;
	const register dim_typ order = s_data->a2;
	ityp * b = s_data->a3;
	ityp * c = s_data->a4;
	ityp * dp2 = s_data->a5;
	ityp * p1 = s_data->a6;
	
	ityp dp0;
	ityp dp1;
	int i;
	ityp p0;
	
	*p1 = *dp2 = 1.00;
	dp1 = 0.00;
	*p2 = x - 1.00;
	
	for ( i = 1; i < order; ++i )
	{
		p0 = *p1;
		dp0 = dp1;
		
		*p1 = *p2;
		dp1 = *dp2;
		
		*p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
		*dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _gen_laguerre_root ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEN_LAGUERRE_ROOT improves a root of a Generalized Laguerre polynomial.
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
    Input/output, double *X, the approximate root, which
    should be improved on output.
    Input, int ORDER, the order of the polynomial to be computed.
    Input, double ALPHA, the exponent of the X factor.
    Output, double *DP2, the value of L'(ORDER)(X).
    Output, double *P1, the value of L(ORDER-1)(X).
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
	const itpitdt4pit * const s_data = data;
	
	const register ityp alpha = s_data->a0;
	ityp * x = s_data->a1;
	const register dim_typ order = s_data->a2;
	ityp * dp2 = s_data->a3;
	ityp * p1 = s_data->a4;
	ityp * b = s_data->a5;
	ityp * c = s_data->a6;
	
	ityp d;
	ityp eps;
	ityp p2;
	dim_typ step;
	const dim_typ step_max = 10;
	
	eps = r8_epsilon ( );
	
	for ( step = 1; step <= step_max; ++step )
	{
		gen_laguerre_recur ( &p2, dp2, p1, *x, order, alpha, b, c );
		
		d = p2 / ( *dp2 );
		*x -= d;
		
		if ( abs ( d ) <= eps * ( abs ( *x ) + 1.0 ) )
			break;
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hermite_root ( void * data)
/******************************************************************************/
/*
  Purpose:
    HERMITE_ROOT improves an approximate root of a Hermite polynomial.
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
    Input/output, double *X, the approximate root, which
    should be improved on output.
    Input, int ORDER, the order of the Hermite polynomial.
    Output, double *DP2, the value of H'(ORDER)(X).
    Output, double *P1, the value of H(ORDER-1)(X).
*/
{
	const dt3pit * const s_data = data;
	
	const register dim_typ order = s_data->a0;
	ityp * dp2 = s_data->a1;
	ityp * p1 = s_data->a2; 
	ityp * x = s_data->a3;
	
	ityp d;
	ityp eps;
	ityp p2;
	dim_typ step;
	const dim_typ step_max = 10;
	
	eps = r8_epsilon ( );
	
	for ( step = 1; step <= step_max; ++step )
	{
		hermite_recur ( &p2, dp2, p1, *x, order );
		
		d = p2 / ( *dp2 );
		*x -= d;
		
		if ( abs ( d ) <= eps * ( abs ( *x ) + 1.0 ) )
			return NULL;
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _jacobi_root ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_ROOT improves an approximate root of a Jacobi polynomial.
  Licesing:
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
    Input/output, double *X, the approximate root, which
    should be improved on output.
    Input, int ORDER, the order of the polynomial to be computed.
    Input, double ALPHA, BETA, the exponents of (1-X) and
 (1+X) in the quadrature rule.
    Output, double *DP2, the value of J'(ORDER)(X).
    Output, double *P1, the value of J(ORDER-1)(X).
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
	const pitdt2it4pit * const s_data = data;
	ityp * x = s_data->a0;
	const register dim_typ order = s_data->a1;
	const register ityp alpha = s_data->a2;
	const register ityp beta = s_data->a3;
	ityp * dp2 = s_data->a4;
	ityp * p1 = s_data->a5;
	ityp * b = s_data->a6;
	ityp * c = s_data->a7;
	
	ityp d;
	ityp eps;
	ityp p2;
	dim_typ step;
	const dim_typ step_max = 10;
	
	eps = r8_epsilon ( );
	
	for ( step = 1; step <= step_max; ++step )
	{
		jacobi_recur ( &p2, dp2, p1, *x, order, alpha, beta, b, c );
		
		d = p2 / ( *dp2 );
		*x -= d;
		
		if ( abs ( d ) <= eps * ( abs ( *x ) + 1.0 ) )
			return NULL;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laguerre_root ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAGUERRE_ROOT improves a root of a Laguerre polynomial.
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
    Input/output, double *X, the approximate root, which
    should be improved on output.
    Input, int ORDER, the order of the polynomial to be computed.
    Output, double *DP2, the value of L'(ORDER)(X).
    Output, double *P1, the value of L(ORDER-1)(X).
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
	const dt5pit  * const s_data = data; 
	
	const register dim_typ order = s_data->a0;
	ityp * x = s_data->a1;
	ityp * dp2 = s_data->a2;
	ityp * p1 = s_data->a3;
	ityp * b = s_data->a4;
	ityp * c = s_data->a5;
	
	ityp d;
	ityp eps;
	ityp p2;
	dim_typ step;
	const dim_typ step_max = 10;
	
	eps = r8_epsilon ( );
	
	for ( step = 1; step <= step_max; ++step )
	{
		laguerre_recur ( &p2, dp2, p1, *x, order, b, c );
		
		d = p2 / ( *dp2 );
		*x -= d;
	
		if ( abs ( d ) <= eps * ( abs ( *x ) + 1.0 ) )
			break;
	}
	
	return NULL;
}



/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm0_cycle ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM0_CYCLE analyzes a permutation of (0,...,N-1).
  Discussion:
    The routine will count cycles, find the sign of a permutation,
    and tag a permutation.
  Example:
    Input:
      N = 9
      IOPT = 1
      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
    Output:
      NCYCLE = 3
      ISGN = +1
      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    28 March 2009
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
    Input, int N, the number of objects being permuted.
    Input/output, int P[N].  On input, P describes a
    permutation, in the sense that entry I is to be moved to P[I].
    If IOPT = 0, then P will not be changed by this routine.
    If IOPT = 1, then on output, P will be "tagged".  That is,
    one element of every cycle in P will be negated.  In this way,
    a user can traverse a cycle by starting at any entry I1 of P
    which is negative, moving to I2 = ABS(P[I1]), then to
    P[I2], and so on, until returning to I1.
    Output, int *ISGN, the "sign" of the permutation, which is
    +1 if the permutation is even, -1 if odd.  Every permutation
    may be produced by a certain number of pairwise switches.
    If the number of switches is even, the permutation itself is
    called even.
    Output, int *NCYCLE, the number of cycles in the permutation.
    Input, int IOPT, requests tagging.
    0, the permutation will not be tagged.
    1, the permutation will be tagged.
*/
{
	const dtpipspib * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	short * isgn = s_data->a2;
	int * ncycle = s_data->a3; 
	const bool iopt = s_data->a4;
	
	int error;
	int i;
	int i1;
	int i2;
	int ierror;
	int is;
	ierror = perm0_check ( n, p );
	
	if ( ierror != 0 )
		return NULL;
	/*
	Increment.
	*/
	for ( i = 0; i < n; ++i, ++p[i] );
	
	is = 1;
	*ncycle = n;
	
	for ( i = 1; i <= n; ++i )
	{
		i1 = p[i-1];
		
		while ( i < i1 )
		{
			-- *ncycle;
			i2 = p[i1-1];
			p[i1-1] = -i2;
			i1 = i2;
		}
		
		if ( iopt != 0 )
			is = - i4_sign ( p[i-1] );
		p[i-1] = abs ( p[i-1] ) * i4_sign ( is );
	}
	
	*isgn = 1 - (( ( n - *ncycle ) % 2 )<<1);
	/*
	Decrement.
	*/
	for ( i = 0; i < n; ++i,--p[i] );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_sqrt_cf ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_SQRT_CF finds the continued fraction representation of a square root of an integer.
  Discussion:
    The continued fraction representation of the square root of an integer
    has the form
      [ B0, (B1, B2, B3, ..., BM), ... ]
    where
      B0 = int ( sqrt ( real ( N ) ) )
      BM = 2 * B0
      the sequence ( B1, B2, B3, ..., BM ) repeats in the representation.
      the value M is termed the period of the representation.
  Example:
     N  Period  Continued Fraction
     2       1  [ 1, 2, 2, 2, ... ]
     3       2  [ 1, 1, 2, 1, 2, 1, 2... ]
     4       0  [ 2 ]
     5       1  [ 2, 4, 4, 4, ... ]
     6       2  [ 2, 2, 4, 2, 4, 2, 4, ... ]
     7       4  [ 2, 1, 1, 1, 4, 1, 1, 4, 1, 1, 4... ]
     8       2  [ 2, 1, 4, 1, 4, 1, 4, 1, 4, ... ]
     9       0  [ 3 ]
    10       1  [ 3, 6, 6, 6, ... ]
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    25 March 2009
  Author:
   John Burkardt
  Reference:
    Mark Herkommer,
    Number Theory, A Programmer's Guide,
    McGraw Hill, 1999, pages 294-307.
  Parameters:
    Input, int N, the number whose continued fraction square root
    is desired.
    Input, int MAX_TERM, the maximum number of terms that may
    be computed.
    Output, int *N_TERM, the number of terms computed beyond the
    0 term.  The routine should stop if it detects that the period
    has been reached.
    Output, int B[MAX_TERM+1], contains the continued fraction
    coefficients for indices 0 through N_TERM.
*/
{
	const _2dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ max_term = s_data->a1;
	int * n_term = s_data->a2;
	int * b = s_data->a3;
	
	int p;
	int q;
	int r;
	int s;
	
	*n_term = 0;
	
	i4_sqrt ( n, &s, &r );
	b[0] = s;
	
	if ( 0 < r )
	{
		p = 0;
		q = 1;
		
		for ( ; ; )
		{
			p = b[*n_term] * q - p;
			q = ( n - p * p ) / q;
			
			if ( max_term <= *n_term )
				return NULL;
			
			++ *n_term;
			b[*n_term] = ( p + s ) / q;
			
			if ( q == 1 )
				break;
		}
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _ubvec_to_ui4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    UBVEC_TO_UI4 makes an unsigned integer from an unsigned binary vector.
  Discussion:
    A UBVEC is a vector of binary digits representing an unsigned integer.  
    UBVEC[N-1] contains the units digit, UBVEC[N-2]
    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
    so that printing the digits in order gives the binary form of the number.
  Example:
    N = 4
         BVEC   binary  I
    ----------  -----  --
    1  2  3  4
    ----------
    0  0  0  1       1  1
    0  0  1  0      10  2
    0  0  1  1      11  3
    0  1  0  0     100  4
    1  0  0  1    1001  9
    1  1  1  1    1111 15
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    21 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    Input, int BVEC[N], the binary representation.
    Input, unsigned int UBVEC_TO_UI4, the integer value.
*/
{
	static unsigned int result = UINT_MAX;
	
	const dtpi * const s_data = data; 
	const register dim_typ n = s_data->a0;
	int * bvec = s_data->a1;
	
	unsigned int ui4 = 0;
	for ( dim_typ i=0; i < n; ++i)
		ui4 = (ui4<<1) + ( unsigned int ) bvec[i];
	
	result = ui4;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _ubvec_xor ( void * data)
/******************************************************************************/
/*
  Purpose:
    UBVEC_XOR computes the exclusive OR of two UBVECs.
  Discussion:
    A UBVEC is a vector of binary digits representing an unsigned integer.  
    UBVEC[N-1] contains the units digit, UBVEC[N-2]
    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
    so that printing the digits in order gives the binary form of the number.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    01 December 2006
  Author:
    John Burkardt
  Parameters:
    Input, int N, the length of the vectors.
    Input, int BVEC1[N], BVEC2[N], the vectors to be XOR'ed.
    Input, int BVEC3[N], the exclusive OR of the two vectors.
*/
{
	const dt3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * bvec1 = s_data->a1;
	int * bvec2 = s_data->a2;
	int * bvec3 = s_data->a3;
	
	for (dim_typ i = 0; i < n; ++i )
		bvec3[i] = ( bvec1[i] + bvec2[i] ) % 2;
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ui4_to_ubvec ( void * data)
/******************************************************************************/
/*
  Purpose:
    UI4_TO_UBVEC makes a unsigned binary vector from an integer.
  Discussion:
    A UBVEC is a vector of binary digits representing an unsigned integer.  
    UBVEC[N-1] contains the units digit, UBVEC[N-2]
    the coefficient of 2, UBVEC[N-3] the coefficient of 4 and so on,
    so that printing the digits in order gives the binary form of the number.
    To guarantee that there will be enough space for any
    value of I, it would be necessary to set N = 32.
  Example:
    UI4      BVEC         binary
        0  1  2  3  4  5
        1  2  4  8 16 32
    --  ----------------  ------
     1  1  0  0  0  0  1       1
     2  0  1  0  0  1  0      10
     3  1  1  0  0  1  1      11
     4  0  0  0  1  0  0     100
     9  0  0  1  0  0  1    1001
    57  1  1  1  0  1  1  110111
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    21 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, unsigned int UI4, an integer to be represented.
    Input, int N, the dimension of the vector.
    Output, int BVEC[N], the binary representation.
*/
{
	const dtpii * const s_data = data;

	const register dim_typ n = s_data->a0;
	int * bvec = s_data->a1;
	register unsigned int ui4 = s_data->a2;	
	
	int i;
	
	for (dim_typ i = n - 1; 0 <= i; --i )
	{
		bvec[i] = ui4 % 2;
		ui4 >>= 1;
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_sqrt ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_SQRT finds the integer square root of N by solving N = Q^2 + R.
  Discussion:
    The integer square root of N is an integer Q such that
    Q**2 <= N but N < (Q+1)^2.
    A simpler calculation would be something like
      Q = INT ( SQRT ( REAL ( N ) ) )
    but this calculation has the virtue of using only integer arithmetic.
    To avoid the tedium of worrying about negative arguments, the routine
    automatically considers the absolute value of the argument.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    22 March 2009
  Author:
   John Burkardt
  Reference:
    Mark Herkommer,
    Number Theory, A Programmer's Guide,
    McGraw Hill, 1999, pages 294-307.
  Parameters:
    Input, int N, the number whose integer square root is desired.
    Actually, only the absolute value of N is considered.
    Output, int *Q, *R, the integer square root, and positive remainder,
    of N.
*/
{
	const dt2pi * const s_data = data;
	register dim_typ n = s_data->a0;
	int * q = s_data->a1;
	int * r = s_data->a2;
	
	n = abs ( n );
	*q = n;
	
	while ( ( n / *q ) < *q )
		*q = ( *q + ( n / *q ) ) / 2;
	
	*r = n - (*q) * (*q);
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triang ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANG renumbers elements in accordance with a partial ordering.
  Discussion:
    TRIANG is given a partially ordered set.  The partial ordering
    is defined by a matrix ZETA, where element I is partially less than
    or equal to element J if and only if ZETA(I,J) = 1.
    TRIANG renumbers the elements with a permutation P so that if
    element I is partially less than element J in the partial ordering,
    then P(I) < P(J) in the usual, numerical ordering.
    In other words, the elements are relabeled so that their labels
    reflect their ordering.  This is equivalent to relabeling the
    matrix so that, on unscrambling it, the matrix would be upper
    triangular.
    Calling I4MAT_2PERM0 or R8MAT_2PERM0 with P used for both the row
    and column permutations applied to matrix ZETA will result in
    an upper triangular matrix.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    11 June 2015
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
    Input, int N, the number of elements in the set.
    Input, int ZETA[N*N], describes the partial ordering.
    ZETA[I,J] =:
      0, for diagonal elements (I = J), or
         for unrelated elements, or
         if J << I.
      1, if I << J.
    Output, int P[N], a permutation of the elements that reflects
    their partial ordering.  P[I] is the new label of element I, with
    the property that if ZETA[I,J] = 1, that is, I << J,
    then P[I] < P[J] (in the usual ordering).
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * zeta = s_data->a1;
	int * p = s_data->a2;
	
	dim_typ i;
	int error;
	int iq;
	int ir;
	int it;
	int l;
	int m;
	/*
	Make sure ZETA represents a partially ordered set.  In other words,
	if ZETA(I,J) = 1, then ZETA(J,I) must NOT be 1.
	*/
	error = pord_check ( n, zeta );
	
	if ( error )
		return NULL;
	
	m = 1;
	l = 0;
	for ( i = 0; i < n; ++i)
		p[i] = 0;

	it = m + 1;
	ir = m + 1;
	
	for ( ; ; )
	{
		if ( ir <= n )
		{
			if ( p[ir-1] == 0 && zeta[(ir-1)+(m-1)*n] != 0 )
			{
				p[ir-1] = m;
				m = ir;
				ir = it;
			}
			else
				++ ir;
		}
		else
		{
			++ l;
			iq = p[m-1];
			p[m-1] = l;
			
			if ( iq != 0 )
			{
				ir = m + 1;
				m = iq;
			}
			else if ( m == n )
				break;
			else
			{
				for ( ; ; )
				{
					++ m;
					
					if ( p[m-1] == 0 )
						break;
		
					if ( m == n )
					{
						for ( i = 0; i < n; ++i, -- p[i] );
						return NULL;
					}
				}
				it = m + 1;
				ir = m + 1;
			}
		}
	}
	/*
	Decrement the elements of the permutation.
	*/
	for ( i = 0; i < n; ++i,--p[i] );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _perm0_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM0_INVERSE inverts a permutation of (0,...,N-1).
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    08 June 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects being permuted.
    Input, int P1[N], the permutation.
    Output, int PERM0_INVERSE[N], the inverse permutation.
*/
{
	const dtpi * const s_data = data; 
	const register dim_typ n = s_data->a0;
	int * p1 = s_data->a1;
	
	dim_typ i;
	int i0;
	int i1;
	int i2;
	int ierror;
	int is;
	int *p2;
	
	if ( n == 0 )
		return NULL;
	
	ierror = perm0_check ( n, p1 );
	
	if ( ierror != 0 )
		return NULL;
	/*
	Temporary increment.
	*/
	p2 = ( int * ) malloc ( n * sizeof ( int ) );
	
	for ( i = 0; i < n; ++i )
		p2[i] = p1[i] + 1;
	
	is = 1;
	
	for ( i = 1; i <= n; ++i )
	{
		i1 = p2[i-1];
		
		while ( i < i1 )
		{
			i2 = p2[i1-1];
			p2[i1-1] = - i2;
			i1 = i2;
		}
		
		is = - i4_sign ( p2[i-1] );
		p2[i-1] = abs ( p2[i-1] ) * i4_sign ( is );

	}

	for ( i = 1; i <= n; ++i )
	{
		i1 = - p2[i-1];
		
		if ( 0 <= i1 )
		{
			i0 = i;
			
			for ( ; ; )
			{
				i2 = p2[i1-1];
				p2[i1-1] = i0;
				
				if ( i2 < 0 )
					break;
				
				i0 = i1;
				i1 = i2;
			}
		}
	}
	/*
	Decrement.
	*/
	for ( i = 0; i < n; ++i, --p2[i]);
	
	return p2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_2perm0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_2PERM0 permutes the rows and columns of a rectangular I4MAT.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    25 March 2009
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
    Input, int M, number of rows in the matrix.
    Input, int N, number of columns in the matrix.
    Input/output, int A[M*N].
    On input, the matrix to be permuted.
    On output, the permuted matrix.
    Input, int P[M], the row permutation.  P(I) is the new number of row I.
    Input, int Q[N].  The column permutation.  Q(I) is the new number
    of column I.  
*/
{
	const _2dt3pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	int * p = s_data->a3;
	int * q = s_data->a4;
	
	dim_typ i;
	int i1;
	short is;
	int it;
	int j;
	int j1;
	int j2;
	int k;
	int lc;
	int nc;
	int *p1;
	int *q1;
	int temp;
	/*
	Wretched maneuvers to deal with necessity of 1-based values,
	and to handle case where P and Q are same vector.
	*/
	p1 = i4vec_copy_new ( m, p );
	perm0_cycle ( m, p1, &is, &nc, 1 );
	for ( i = 0; i < m; ++i, ++p1[i] );
	
	q1 = i4vec_copy_new ( n, q );
	perm0_cycle ( n, q1, &is, &nc, 1 );
	for ( j = 0; j < n; ++j, ++q1[j]);
	
	for ( i = 1; i <= m; i++ )
	{
		i1 = - p1[i-1];
		
		if ( 0 < i1 )
		{
			lc = 0;
			
			for ( ; ; )
			{
				i1 = p1[i1-1];
				++ lc;
				
				if ( i1 <= 0 )
					break;
			}
		
			i1 = i;
			
			for ( j = 1; j <= n; ++j )
			{
				if ( q1[j-1] <= 0 )
				{
					j2 = j;
					k = lc;
					
					for ( ; ; )
					{
						j1 = j2;
						it = a[i1-1+(j1-1)*m];
						
						for ( ; ; )
						{
							i1 = abs ( p1[i1-1] );
							j1 = abs ( q1[j1-1] );
							
							temp = it;
							it = a[i1-1+(j1-1)*m];
							a[i1-1+(j1-1)*m] = temp;
							
							if ( j1 != j2 )
								continue;
		
							-- k;
		
							if ( i1 == i )
								break; 
						}
							
						j2 = abs ( q1[j2-1] );
						
						if ( k == 0 )
							break;
					}
				}
			}
		}
	}
	/*
	Discard the 1-based permutations.
	*/
	free ( p1 );
	free ( q1 );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_bset ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_BSET returns a copy of an I4 in which the POS-th bit is set to 1.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    23 January 2007
  Author:
    John Burkardt
  Reference:
    Military Standard 1753,
    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
    9 November 1978.
  Parameters:
    Input, int I4, the integer to be tested.
    Input, int POS, the bit position, between 0 and 31.
    Output, int I4_BSET, a copy of I4, but with the POS-th bit
    set to 1.
*/
{
	static int result = INT_MAX;
	
	const idt * const s_data = data; 
	int i4 = s_data->a0;
	const register dim_typ pos = s_data->a1;
	
	int add;
	int j;
	dim_typ k;
	int value;
	
	value = i4;
	
	if ( pos < 0 );
	else if ( pos < 31 )
	{
		add = 1;
		j = 0<=i4 ? i4 : ( i4_huge + i4 ) + 1;
	
		for ( k = 1; k <= pos; ++k )
		{
			j /= 2;
			add <<= 1;
		}
		
		if ( ( j % 2 ) == 0 )
			value = i4 + add;
	}
	else if ( pos == 31 && 0 < i4 )
			value = - ( i4_huge - i4 ) - 1;
	else if ( 31 < pos )
		value = i4;
		
	result = value; 
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_btest ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_BTEST returns TRUE if the POS-th bit of an I4 is 1.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    23 January 2007
  Author:
    John Burkardt
  Reference:
    Military Standard 1753,
    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
    9 November 1978.
  Parameters:
    Input, int I4, the integer to be tested.
    Input, int POS, the bit position, between 0 and 31.
    Output, int I4_BTEST, is TRUE if the POS-th bit of I4 is 1.
*/
{
	static int result = INT_MAX;
	
	const idt * const s_data = data; 
	int i4 = s_data->a0;
	const register dim_typ pos = s_data->a1;
	
	int j;
	dim_typ k;
	int value;
	
	if ( pos == 0 || 31 < pos )
	{
		result = INT_MAX;
		return &result;
	}
	else if ( pos < 31 )
	{
		j =  0 <= i4 ? i4 : ( i4_huge + i4 ) + 1;
		for ( k = 1; k <= pos; ++k, j /=2);
		value = j%2 != 0;
	}
	else if ( pos == 31 )
		value = i4<0;
		
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _index_rank0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    INDEX_RANK0 ranks an index vector within given upper limits.
  Example:
    N = 3,
    HI = 3
    A = ( 3, 1, 2 )
    RANK = 12
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    24 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, int HI, the upper limit for the array indices.
    The lower limit is implicitly 1, and HI should be at least 1.
    Input, int A[N], the index vector to be ranked.
    Output, int INDEX_RANK0, the rank of the index vector, or -1 if A
    is not a legal index.
*/
{
	static short result = SHRT_MAX;
	
	const _2dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ hi = s_data->a1;
	int * a = s_data->a2;
	
	dim_typ i;
	int range;
	int rank;
	
	rank = -1;
	for ( i = 0; i < n; ++i )
		if ( a[i] < 1 || hi < a[i] )
		{
			result = rank;
			return &result;
		}
	
	rank = 0;
	for ( i = n-1; 0 <= i; --i )
		rank = hi * rank + a[i];
	
	rank = range = 1;
	for ( i = 0; i < n; ++i )
	{
		rank += ( a[i] - 1 ) * range;
		range *= hi;
	}
	
	result = rank;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _index_unrank0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    INDEX_UNRANK0 unranks an index vector within given upper limits.
  Example:
    N = 3,
    HI = 3
    RANK = 12
    A = ( 3, 1, 2 )
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    24 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, int HI, the upper limit for the array indices.
    The lower limit is implicitly 1, and HI should be at least 1.
    Input, int RANK, the rank of the desired index vector.
    Output, int A[N], the index vector of the given rank.
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ hi = s_data->a1;
	const register dim_typ rank = s_data->a2;
	int * a = s_data->a3;
	
	dim_typ i, j;
	int k;
	int range;
	
	for ( i = 0; i < n; ++i )
		a[i] = 0;
	/*
	The rank might be too small.
	*/
	if ( rank == 0 )
		return NULL;
	
	range = powi ( hi, n );
	/*
	The rank might be too large.
	*/
	if ( range < rank )
		return NULL;
	
	k = rank - 1;
	
	for ( i = n-1; 0 <= i; --i )
	{
		range /= hi;
		j = k / range;
		a[i] = j + 1;
		k -= j * range;
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _perm0_lex_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM0_LEX_NEXT generates permutations of (0,...,N-1) in lexical order.
  Example:
    N = 3
    1   0 1 2
    2   0 2 1
    3   1 0 2
    4   1 2 0
    5   2 0 1
    6   2 1 0
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    28 March 2009
  Author:
    Original FORTRAN77 version by Mok-Kong Shen.
    C version by John Burkardt.
  Reference:
    Mok-Kong Shen,
    Algorithm 202: Generation of Permutations in Lexicographical Order,
    Communications of the ACM,
    Volume 6, September 1963, page 517.
  Parameters:
    Input, int N, the number of elements being permuted.
    Input/output, int P[N], the permutation, in standard index form.
    The user should not alter the elements of Pbetween successive
    calls.  The routine uses the input value of P to determine
    the output value.
    Input/output, int *MORE.
    On the first call, the user should set MORE =FALSE which signals
    the routine to do initialization.
    On return, if MORE is TRUE then another permutation has been
    computed and returned, while if MORE is FALSE there are no more
    permutations.
*/
{
	const dtpipb * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	bool * more = s_data->a2;
	
	dim_typ j, k;
	int temp;
	int u;
	int w;
	/*
	Initialization.
	*/
	if ( !( *more ) )
	{
		i4vec_indicator0 ( n, p );
		*more = true;
	}
	else
	{
		if ( n <= 1 )
		{
			*more = false;
			return NULL;
		}
	
		w = n;
	
		while ( p[w-1] < p[w-2] )
		{
			if ( w == 2 )
			{
				*more = false;
				return NULL;
			}
	
			-- w;
		}
	
		u = p[w-2];
		
		for ( j = n; w <= j; --j )
		{
			if ( u < p[j-1] )
			{
				p[w-2] = p[j-1];
				p[j-1] = u;
				
				for ( k = 0; k <= ( n - w - 1 ) / 2; ++k )
				{
					temp = p[n-k-1];
					p[n-k-1] = p[w+k-1];
					p[w+k-1] = temp;
				}
				return NULL;
			}
		}
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm0_free ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM0_FREE reports the unused items in a partial permutation of (0,...,N-1).
  Discussion:
    It is assumed that the N objects being permuted are the integers
    from 0 to N-1, and that IPART contains a "partial" permutation, that
    is, the NPART entries of IPART represent the beginning of a
    permutation of all N items.
    The routine returns in IFREE the items that have not been used yet.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    28 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int NPART, the number of entries in IPART.  NPART may be 0.
    Input, int IPART[NPART], the partial permutation, which should
    contain, at most once, some of the integers between 0 and
    NPART+NFREE-1.
    Input, int NFREE, the number of integers that have not been
    used in IPART.  This is simply N - NPART.  NFREE may be zero.
    Output, int IFREE[NFREE], the integers between 1 and NPART+NFREE
    that were not used in IPART.
*/
{
	const _2dt2pi * const s_data = data;
	
	const register dim_typ npart = s_data->a0;
	const register dim_typ nfree = s_data->a1;
	int * ipart = s_data->a2;
	int * ifree = s_data->a3;
	
	dim_typ i, j, k;
	int match;
	int n;
	
	n = npart + nfree;
	
	if ( npart == 0 )
		i4vec_indicator0 ( n, ifree );
	else if ( nfree == 0 )
		return NULL;
	else
	{
		k = 0;
		
		for ( i = 0; i < n; ++i )
		{
			match = -1;
			
			for ( j = 0; j < npart; ++j )
				if ( ipart[j] == i )
				{
					match = j;
					break;
				}
			
			if ( match == -1 )
			{
				++ k;
			
				if ( nfree < k )
					return NULL;
				ifree[k-1] = i;
			}
		}
	}

	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _ksub_random ( void * data)
/******************************************************************************/
/*
  Purpose:
    KSUB_RANDOM selects a random subset of size K from a set of size N.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    30 April 2003
  Author:
    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the size of the set from which subsets are drawn.
    Input, int K, number of elements in desired subsets.  K must
    be between 0 and N.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int A[K].  A(I) is the I-th element of the
    output set.  The elements of A are in order.
*/
{
	const _2dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	int * seed = s_data->a2;
	int * a = s_data->a3;
	
	dim_typ i;
	int ids;
	int ihi;
	int ip;
	int ir;
	int is;
	int ix;
	int l;
	int ll;
	int m;
	int m0;
	
	if ( n < k || k == 0)
		return NULL;
	
	for ( i = 1; i <= k; ++i )
		a[i-1] = ( ( i - 1 ) * n ) / k;
	
	for ( i = 1; i <= k; i++ )
	{
		for ( ; ; )
		{
			ix = i4_uniform_ab ( 1, n, seed );
			l = 1 + ( ix * k - 1 ) / n;
			
			if ( a[l-1] < ix )
				break;
		}
	
		++ a[l-1];
	
	}
	
	ip = 0;
	is = k;
	
	for ( i = 1; i <= k; ++i )
	{
		m = a[i-1];
		a[i-1] = 0;
		
		if ( m != ( ( i - 1 ) * n ) / k )
		{
			++ ip;
			a[ip-1] = m;
		}
	
	}
	
	ihi = ip;
	
	for ( i = 1; i <= ihi; ++i )
	{
		ip = ihi + 1 - i;
		l = 1 + ( a[ip-1] * k - 1 ) / n;
		ids = a[ip-1] - ( ( l - 1 ) * n ) / k;
		a[ip-1] = 0;
		a[is-1] = l;
		is -= ids;
	}
	
	for ( ll = 1; ll <= k; ++ll )
	{
		l = k + 1 - ll;
		
		if ( a[l-1] != 0 )
		{
			ir = l;
			m0 = 1 + ( ( a[l-1] - 1 ) * n ) / k;
			m = ( a[l-1] * n ) / k - m0 + 1;
		}
		/*
		There is something wrong with this algorithm!
		If A[L-1] is zero, then the values of IR, M0, and M are not defined
		on this loop iteration, and hence are either STALE values from the
		previous iteration, or UNDEFINED if this is the first pass.
		JVB, 21 December 2014.
		*/
		ix = i4_uniform_ab ( m0, m0 + m - 1, seed );
		
		i = l + 1;
		
		while ( i <= ir )
		{
			if ( ix < a[i-1] )
				break;
			
			++ ix;
			a[i-2] = a[i-1];
			++ i;
		}
		a[i-2] = ix;
		-- m;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subset_sum_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_SUM_COUNT counts solutions to the subset sum problem in a range.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the subset.
    Input, int W[N], a set of weights.  The length of this
    array must be no more than 31.
    Input, int T, the target value.
    Input, int IND_MIN, IND_MAX, the lower and upper
    limits to be searched.  0 <= IND_MIN <= IND_MAX <= (2^N)-1.
    Output, int SUBSET_SUM_COUNT, the number of distinct
    solutions of the subset sum problem found within the given range.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpi2dt * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ t = s_data->a1;
	int * w = s_data->a2;
	const register dim_typ ind_min = s_data->a3;
	const register dim_typ ind_max = s_data->a4;
	
    int *c;
    dim_typ ind;
    dim_typ ind_max2;
    dim_typ ind_min2;
    dim_typ solution_num;
    /*
    Check the data.
    */

    if ( 31 < n )
    {
    	result = USHRT_MAX;
        return &result;
    }

    ind_min2 = MAX ( ind_min, 0 );
    ind_max2 = MIN ( ind_max, powi ( 2, n ) - 1 );
    /*
    Run through the range.
    */

    solution_num = 0;

    for ( ind = ind_min2; ind <= ind_max2; ++ind )
    {
        /*
        Convert INDEX into vector of indices in W.
        */
        c = i4_to_digits_binary ( ind, n );
        /*
        If the sum of those weights matches the target, return combination.
        */
        if ( i4vec_dot_product ( n, c, w ) == t )
            ++ solution_num ;

        free ( c );
    }

	result = solution_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _subset_sum_find ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_SUM seeks a subset of a set that has a given sum.
  Discussion:
    This function tries to compute a target value as the sum of
    a selected subset of a given set of weights.
    This function works by brute force, that is, it tries every
    possible subset to see if it sums to the desired value.
    Given N weights, every possible selection can be described by
    one of the N-digit binary numbers from 0 to 2^N-1.
    This function includes a range, which allows the user to
    control which subsets are to be checked.  Thus, if there are
    N weights, specifying a range of [ 0, 2^N-1] indicates that
    all subsets should be checked.  On the other hand, this full
    range could be broken down into smaller subranges, each of
    which could be checked independently.
    It is possible that, in the given range, there may be multiple
    solutions of the problem.  This function will only return
    one such solution, if found.  However, the function may be called
    again, with an appropriate restriction of the range, to continue
    the search for other solutions.
  Example:
    w = [ 1, 2, 4, 8, 16, 32 ];
    t = 22;
    r = [ 0, 2^6 - 1 ];
    call subset_sum ( w, t, r, c, ind )
    c = [ 2, 3, 5 ]
    index = 22
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the subset.
    Input, int W[N], a set of weights.  The length of this
    array must be no more than 31.
    Input, int T, the target value.
    Input, int IND_MIN, IND_MAX, the lower and upper
    limits to be searched.  0 <= IND_MIN <= IND_MAX <= (2^N)-1.
    Output, int *IND, the index of the solution.
    If IND is -1, no solution was found in the range.
    Output, int SUBSET_SUM_FIND[N], indicates the solution, assuming
    that IND is not -1.  In that case, the sum T is made by selecting
    those weights W(I) for which C(I) is 1.  In fact,
    T = sum ( 1 <= I <= N ) C(I) * W(I).
*/
{
	const _2dtpi2dtpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ t = s_data->a1;
	int * w = s_data->a2;
	const register dim_typ ind_min = s_data->a3;
	const register dim_typ ind_max = s_data->a4;
	int * ind = s_data->a5;
	
    int *c;
    dim_typ ind_max2;
    dim_typ ind_min2;
    /*
    Check the data.
    */

    if ( 31 < n )
        return NULL;

    ind_min2 = MAX ( ind_min, 0 );
    ind_max2 = MIN ( ind_max, powi ( 2, n ) - 1 );
    /*
    Run through the range.
    */

    for ( *ind = ind_min2; *ind <= ind_max2; ++(*ind) )
    {
        /*
        Convert INDEX into vector of indices in W.
        */
        c = i4_to_digits_binary ( *ind, n );
        /*
        If the sum of those weights matches the target, return combination.
        */
        if ( i4vec_dot_product ( n, c, w ) == t )
            return c;

        free ( c );
    }

    *ind = - 1;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subset_sum_table_to_list_length ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_SUM_TABLE_TO_LIST converts a subset sum table to a list.
  Discussion:
    The subset sum problem seeks to construct the value T by summing a
    subset of the values W.
    This function takes a table computed by subset_sum_table() and converts
    it to the corresponding list of values that form the sum.
  Example:
    t = 22
    n = 6
    w = [ 1, 2, 4, 8, 16, 32 ]
    call subset_sum ( t, n, w, table )
    table = [ 1, 2, 2, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8,
      16, 16, 16, 16, 16, 16, 16 ]
    call subset_sum_table_to_list ( t, table, m, list )
    m = 3
    list = [ 2, 4, 16 ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2015
  Author:
    John Burkardt
  Parameters:
    Input, integer T, the target value.
    Input, integer TABLE[T], the subset sum table.
    Output, int SUBSET_SUM_TABLE_TO_LIST_LENGTH, the number of items in the list.
    If M == 0, then no solution was found and the list is empty.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ t = s_data->a0;
	int * table = s_data->a1;
	
    dim_typ i = t;
    dim_typ m = 0;
    /*
    Determine length of list.
    */
    while ( 0 < i )
    {
        ++ m ;
        i -= table[i];
    }
    
    result = m;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _subset_sum_table ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_SUM_TABLE sets a subset sum table.
  Discussion:
    The subset sum problem seeks to construct the value T by summing a
    subset of the values W.
    This function seeks a solution by constructing a table TABLE of length T,
    so that TABLE(I) = J means that the sum I can be constructed, and that
    the last member of the sum is an entry of W equal to J.
  Example:
    w = [ 1, 2, 4, 8, 16, 32 ];
    t = 22;
    table = subset_sum ( w, t, r )
    table = [ 1, 2, 2, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8,
      16, 16, 16, 16, 16, 16, 16 ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2015
  Author:
    John Burkardt
  Parameters:
    Input, int T, the target value.
    Input, int N, the number of weights.
    Input, int W[N], the weights.
    Output, int SUBSET_SUM_TABLE[T], the subset sum table.  TABLE(I) is
    0 if the target value I cannot be formed.  It is J if the value I can
    be formed, with the last term in the sum being the value J.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ t = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * w = s_data->a2;
	
    dim_typ i, j;
    int *table =table = ( int * ) malloc ( ( t + 1 ) * sizeof ( int ) );

    for ( i = 0; i <= t; ++i )
        table[i] = 0;

    for ( i = 0; i < n; ++i)
        for ( j = t - w[i]; 0 <= j; --j)
            if ( j == 0 && table[w[i]] == 0 )
                table[w[i]] = w[i];
            else if ( table[j] != 0 && table[j+w[i]] == 0 )
                table[j+w[i]] = w[i];

    return table;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _subset_sum_table_to_list ( void * data)
/******************************************************************************/
/*
  Purpose:
    SUBSET_SUM_TABLE_TO_LIST converts a subset sum table to a list.
  Discussion:
    The subset sum problem seeks to construct the value T by summing a
    subset of the values W.
    This function takes a table computed by subset_sum_table() and converts
    it to the corresponding list of values that form the sum.
  Example:
    t = 22
    n = 6
    w = [ 1, 2, 4, 8, 16, 32 ]
    call subset_sum ( t, n, w, table )
    table = [ 1, 2, 2, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8,
      16, 16, 16, 16, 16, 16, 16 ]
    call subset_sum_table_to_list ( t, table, m, list )
    m = 3
    list = [ 2, 4, 16 ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2015
  Author:
    John Burkardt
  Parameters:
    Input, int T, the target value.
    Input, int TABLE(T), the subset sum table.
    Input, int M, the number of items in the list.
    Output, int LIST[M], the weights that sum to T.
*/
{
	const _2dtpi * const s_data = data;
	
	const register dim_typ t = s_data->a0;
	register dim_typ m = s_data->a1;
	int * table = s_data->a2;
	
    dim_typ i = t;
    int *list;
    /*
    Create list.
    */
    list = ( int * ) malloc ( m * sizeof ( int ) );
    /*
    Populate list.
    */
    while ( 0 < i )
    {
        list[m] = table[i];
        ++ m;
        i -= table[i];
    }
    return list;
}

#endif
