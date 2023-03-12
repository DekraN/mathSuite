#ifndef __DISABLEDEEP_WISHART

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_exponential_01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_EXPONENTIAL_01_SAMPLE samples the standard exponential PDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 April 2013
  Author:
    John Burkardt
  Parameters:
    Output, double R8_EXPONENTIAL_01_SAMPLE, a sample of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	result = - log ( r8_uniform_01_sample ( ) );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gamma_01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAMMA_01_SAMPLE samples the standard Gamma distribution.
  Discussion:
    This procedure corresponds to algorithm GD in the reference.
    pdf ( a; x ) = 1/gamma(a) * x^(a-1) * exp ( - x )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 April 2013
  Author:
    Original FORTRAN77 version by Barry Brown, James Lovato.
    C version by John Burkardt.
  Reference:
    Joachim Ahrens, Ulrich Dieter,
    Generating Gamma Variates by a Modified Rejection Technique,
    Communications of the ACM,
    Volume 25, Number 1, January 1982, pages 47-54.
  Parameters:
    Input, double A, the shape parameter of the standard gamma
    distribution.  0.0 < A.
    Output, double R8_GAMMA_01_SAMPLE, a random deviate from the distribution.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp a = *(ityp *) data;
	
	ityp a1 =  0.3333333;
	ityp a2 = -0.2500030;
	ityp a3 =  0.2000062;
	ityp a4 = -0.1662921;
	ityp a5 =  0.1423657;
	ityp a6 = -0.1367177;
	ityp a7 =  0.1233795;
	ityp b;
	ityp c;
	ityp d;
	ityp e;
	ityp e1 = 1.00;
	ityp e2 = 0.4999897;
	ityp e3 = 0.1668290;
	ityp e4 = 0.0407753;
	ityp e5 = 0.0102930;
	ityp p;
	ityp q;
	ityp q0;
	ityp q1 =  0.04166669;
	ityp q2 =  0.02083148;
	ityp q3 =  0.00801191;
	ityp q4 =  0.00144121;
	ityp q5 = -0.00007388;
	ityp q6 =  0.00024511;
	ityp q7 =  0.00024240;
	ityp r;
	ityp s;
	ityp s2;
	ityp si;
	ityp sqrt32 = 5.6568542494923801952;
	ityp t;
	ityp u;
	ityp v;
	ityp value;
	ityp w;
	ityp x;
	
	if ( 1.00 <= a )
	{
		s2 = a - 0.50;
		s = sqrt ( s2 );
		d = sqrt32 - 12.00 * s;
		/*
		Immediate acceptance.
		*/
		t = r8_normal_01_sample ( );
		x = s + 0.5 * t;
		value = x * x;
		
		if ( 0.00 <= t )
		{
			result = value;
			return &result;
		}
		/*
		Squeeze acceptance.
		*/
		u = r8_uniform_01_sample ( );
		if ( d * u <= t * t * t )
		{
			result = value;
			return &result;
		}
		
		r = 1.00 / a;q0 = (((((( q7   * r + q6 ) * r + q5 ) * r + q4 ) * r+ q3 ) * r+ q2 ) * r+ q1 ) * r;
		/*
		Approximation depending on size of parameter A.
		*/
		if ( 13.022 < a )
		{
			b = 1.77;
			si = 0.75;
			c = 0.1515 / s;
		}
		else if ( 3.686 < a )
		{
			b = 1.654 + 0.0076 * s2;
			si = 1.68 / s + 0.275;
			c = 0.062 / s + 0.024;
		}
		else
		{
			b = 0.463 + s + 0.178 * s2;
			si = 1.235;
			c = 0.195 / s - 0.079 + 0.16 * s;
		}
		/*
		Quotient test.
		*/
		if ( 0.00 < x )
		{
			v = 0.50 * t / s;
			q = 0.25 < fabs ( v ) ? q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v ) : q0 + 0.50 * t * t * (((((( a7   * v + a6 ) * v + a5 ) * v + a4 ) * v + a3 ) * v + a2 ) * v + a1 ) * v;
	
			if ( log ( 1.00 - u ) <= q )
			{
				result = value;
				return &result;
			}
		}
	
		for ( ; ; )
		{
			e = r8_exponential_01_sample ( );
			u = 2.00 * r8_uniform_01_sample ( ) - 1.00;
	
			t = 0.00 <= u ? b + fabs ( si * e ) : b - fabs ( si * e );
			/*
			Possible rejection.
			*/
			if ( t < -0.7187449 )
				continue;
			/*
			Calculate V and quotient Q.
			*/
			v = 0.50 * t / s;
			q = 0.25 < fabs(v) ? q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v ) : q0 + 0.50 * t * t * (((((( a7   * v + a6 ) * v + a5 ) * v + a4 ) * v + a3 ) * v + a2 ) * v + a1 ) * v;
			/*
			Hat acceptance.
			*/
			if ( q <= 0.00 )
				continue;
		
			w = 0.50<q ? exp(q) - 1.00 : (((( e5   * q 	+ e4 ) * q + e3 ) * q + e2 ) * q + e1 ) * q;
			/*
			May have to sample again.
			*/
			if ( c * fabs ( u ) <= w * exp ( e - 0.5 * t * t ) )
				break;
		}
	
		x = s + 0.50 * t;
		value = x * x;
	}
	/*
	Method for A < 1.
	*/
	else if ( a < 1.00 )
	{
		b = 1.00 + 0.3678794 * a;
		
		for ( ; ; )
		{
			p = b * r8_uniform_01_sample ( );
			
			if ( p < 1.00 )
			{
				value = exp ( log ( p ) / a );
				if ( value <= r8_exponential_01_sample ( ) )
					break;
			}
			else
			{
				value = - log ( ( b - p ) / a );
				if ( ( 1.00 - a ) * log ( value ) <= r8_exponential_01_sample ( ) )
					break;
			}
		}
	}
	
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gamma_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAMMA_SAMPLE generates a Gamma random deviate.
  Discussion:
    This procedure generates random deviates from the gamma distribution whose
    density is (A^R)/Gamma(R) * X^(R-1) * Exp(-A*X)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 April 2013
  Author:
    Original FORTRAN77 version by Barry Brown, James Lovato.
    C version by John Burkardt.
  Reference:
    Joachim Ahrens, Ulrich Dieter,
    Generating Gamma Variates by a Modified Rejection Technique,
    Communications of the ACM,
    Volume 25, Number 1, January 1982, pages 47-54.
    Joachim Ahrens, Ulrich Dieter,
    Computer Methods for Sampling from Gamma, Beta, Poisson and
    Binomial Distributions,
    Computing,
    Volume 12, Number 3, September 1974, pages 223-246.
  Parameters:
    Input, double A, the location parameter.
    A nonzero.
    Input, double R, the shape parameter.
    0.0 < R.
    Output, double R8_GAMMA_SAMPLE, a random deviate from the distribution.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp r = a_data[1];
	
	result = r8_gamma_01_sample ( r ) / a;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_mmt_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_MMT_NEW computes C = A * B'.
  Discussion:
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
    For this routine, the result is returned as the function value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 November 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, N3, the order of the matrices.
    Input, double A[N1*N2], double B[N3*N2], the matrices to multiply.
    Output, double R8MAT_MMT[N1*N3], the product matrix C = A * B'.
*/
{
	const _3dt2pit * const s_data = data;
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	const register dim_typ n3 = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	
	ityp *c;
	dim_typ i, j, k;
	
	c = ( ityp * ) malloc ( n1 * n3 * sizeof ( ityp ) );
	
	for ( i = 0; i < n1; ++i)
		for ( j = 0; j < n3; ++j )
		{
			c[i+j*n1] = 0.00;
			for ( k = 0; k < n2; ++k )
				c[i+j*n1] += a[i+k*n1] * b[j+k*n3];
		}
	
	return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8ut_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8UT_INVERSE computes the inverse of a R8UT matrix.
  Discussion:
    The R8UT storage format is used for an M by N upper triangular matrix,
    and allocates space even for the zero entries.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    18 February 2013
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6
  Parameters:
    Input, int N, the order of the matrix.
    Input, double A[N*N], the R8UT matrix.
    Output, double R8UT_INVERSE[N*N], the inverse of the upper 
    triangular matrix.
*/
{
	const dtpit * const s_data = data;
	register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
	ityp *b;
	dim_typ i, j, k;
	/*
	Check.
	*/
	for ( i = 0; i < n; ++i)
		if ( a[i+i*n] == 0.00 )
			return NULL;
		
	b = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
	
	for ( j = n-1; 0 <= j; --j)
	{
		for ( i = n-1; 0 <= i; --i )
		{
			if ( j < i )
				b[i+j*n] = 0.00;
			else if ( i == j )
				b[i+j*n] = 1.00 / a[i+j*n];
			else if ( i < j )
			{
				b[i+j*n] = 0.00;
				
				for ( k = i+1; k <= j; ++k)
					b[i+j*n] -= a[i+k*n] * b[k+j*n];
				b[i+j*n] /= a[i+i*n];
			}
		}
	}
	
	return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lg_get ( void * data)
/******************************************************************************/
/*
  Purpose:
    LG_GET queries the LG values for a given generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int G, the index of the generator.
    0 <= G <= 31.
    Output, int *LG1, *LG2, the LG values for generator G.
*/
{
	const i2pi * const s_data = data;
	const register int g = s_data->a0;
	int * lg1 = s_data->a1;
	int * lg2 = s_data->a2;
	
	lg_memory ( -1, g, lg1, lg2 );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lg_memory ( void * data)
/******************************************************************************/
/*
  Purpose:
    LG_MEMORY stores the LG values for all generators.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, the desired action.
    -1, get a value.
    0, initialize all values.
    1, set a value.
    Input, int G, for I = -1 or +1, the index of 
    the generator, with 0 <= G <= 31.
    Input/output, int *LG1, *LG2.  For I = -1, 
    these are output, for I = +1, these are input, for I = 0,
    these arguments are ignored.  When used, the arguments are
    old or new values of the LG parameter for generator G.
*/
{
	const _2i2pi * const s_data = data;
	int i = s_data->a0;
	int g = s_data->a1;
	int * lg1 = s_data->a2;
	int * lg2 = s_data->a3;
	
	# define G_MAX 32
	
	const int g_max = 32;
	dim_typ j;
	static int lg1_save[G_MAX];
	static int lg2_save[G_MAX];
	
	if ( g < 0 || g_max <= g )
		return NULL;
	
	if ( i < 0 )
	{
		*lg1 = lg1_save[g];
		*lg2 = lg2_save[g];
	}
	else if ( i == 0 )
		for ( j = 0; j < g_max; ++j )
			lg1_save[j] = lg2_save[j] = 0;
	else if ( 0 < i )
	{
		lg1_save[g] = *lg1;
		lg2_save[g] = *lg2;
	}
	
	return NULL;
	# undef G_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lg_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    LG_SET sets the LG values for a given generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int G, the index of the generator.
    0 <= G <= 31.
    Input, int LG1, LG2, the LG values for generator G.
*/
{
	int * const a_data = data;
	int g = a_data[0]; 
	int lg1 = a_data[1];
	int lg2 = a_data[2];
	
	lg_memory ( +1, g, &lg1, &lg2 );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ig_get ( void * data)
/******************************************************************************/
/*
  Purpose:
    IG_GET queries the IG values for a given generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int G, the index of the generator.
    0 <= G <= 31.
    Output, int *IG1, *IG2, the IG values for generator G.
*/
{
	const i2pi * const s_data = data;
	const register int g = s_data->a0;
	int * ig1 = s_data->a1;
	int * ig2 = s_data->a2;
	
	ig_memory ( -1, g, ig1, ig2 );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ig_memory ( void * data)
/******************************************************************************/
/*
  Purpose:
    IG_MEMORY stores the IG values for all generators.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, the desired action.
    -1, get a value.
    0, initialize all values.
    1, set a value.
    Input, int G, for I = -1 or +1, the index of 
    the generator, with 0 <= G <= 31.
    Input/output, int *IG1, *IG2.  For I = -1, 
    these are output, for I = +1, these are input, for I = 0,
    these arguments are ignored.  When used, the arguments are
    old or new values of the IG parameter for generator G.
*/
{
	const _2i2pi * const s_data = data;
	int i = s_data->a0;
	int g = s_data->a1;
	int * ig1 = s_data->a2;
	int * ig2 = s_data->a3;
	
	# define G_MAX 32
	
	const int g_max = 32;
	static int ig1_save[G_MAX];
	static int ig2_save[G_MAX];
	int j;
	
	if ( g < 0 || g_max <= g )
		return NULL;
	
	if ( i < 0 )
	{
		*ig1 = ig1_save[g];
		*ig2 = ig2_save[g];
	}
	else if ( i == 0 )
		for ( j = 0; j < g_max; ++j )
			ig1_save[j] = ig2_save[j] = 0;
	else if ( 0 < i )
	{
		ig1_save[g] = *ig1;
		ig2_save[g] = *ig2;
	}
	
	return NULL;
	# undef G_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ig_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    IG_SET sets the IG values for a given generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int G, the index of the generator.
    0 <= G <= 31.
    Input, int IG1, IG2, the IG values for generator G.
*/
{
	int * const a_data = data;
	int g = a_data[0]; 
	int ig1 = a_data[1];
	int ig2 = a_data[2];
	
  ig_memory ( +1, g, &ig1, &ig2 );
  return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _init_generator ( void * data)
/******************************************************************************/
/*
  Purpose:
    INIT_GENERATOR sets the state of generator G to initial, last or new seed.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2013
  Author:
    Original Pascal version by Pierre L'Ecuyer, Serge Cote.
    C version by John Burkardt.
  Reference:
    Pierre LEcuyer, Serge Cote,
    Implementing a Random Number Package with Splitting Facilities,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 1, March 1991, pages 98-111.
  Parameters:
    Input, int T, the seed type:
    0, use the seed chosen at initialization time.
    1, use the last seed.
    2, use a new seed set 2^30 values away.
*/
{
	int t = *(int *) data;
	
	const int a1_w = 1033780774;
	const int a2_w = 1494757890;
	int cg1;
	int cg2;
	int g;
	int ig1;
	int ig2;
	int lg1;
	int lg2;
	const int m1 = 2147483563;
	const int m2 = 2147483399;
	/*
	Check whether the package must be initialized.
	*/
	if ( ! initialized_get ( ) )
		initialize ( );
	/*
	Get the current generator index.
	*/
	g = cgn_get ( );
	/*
	0: restore the initial seed.
	*/
	if ( t == 0 )
	{
		ig_get ( g, &ig1, &ig2 );
		lg1 = ig1;
		lg2 = ig2;
		lg_set ( g, lg1, lg2 );
	}
	/*
	1: restore the last seed.
	*/
	else if ( t == 1 )
		lg_get ( g, &lg1, &lg2 );
	/*
	2: advance to a new seed.
	*/
	else if ( t == 2 )
	{
		lg_get ( g, &lg1, &lg2 );
		lg1 = multmod ( a1_w, lg1, m1 );
		lg2 = multmod ( a2_w, lg2, m2 );
		lg_set ( g, lg1, lg2 );
	}
	else
		return NULL;
	/*
	Store the new seed.
	*/
	cg1 = lg1;
	cg2 = lg2;
	cg_set ( g, cg1, cg2 );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _initialize ( void * data)
/******************************************************************************/
/*
  Purpose:
    INITIALIZE initializes the random number generator library.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2013
  Author:
    Original Pascal version by Pierre L'Ecuyer, Serge Cote.
    C version by John Burkardt.
  Reference:
    Pierre LEcuyer, Serge Cote,
    Implementing a Random Number Package with Splitting Facilities,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 1, March 1991, pages 98-111.
  Parameters:
    None
*/
{
	dim_typ g;
	const int g_max = 32;
	int ig1;
	int ig2;
	int value;
	/*
	Remember that we have called INITIALIZE().
	*/
	initialized_set ( );
	/*
	Initialize all generators to have FALSE antithetic value.
	*/
	value = 0;
	for ( g = 0; g < g_max; ++g )
	{
		cgn_set ( g );
		antithetic_set ( value );
	}
	/*
	Set the initial seeds.
	*/
	ig1 = 1234567890;
	ig2 = 123456789;
	set_initial_seed ( ig1, ig2 );
	/*
	Initialize the current generator index to 0.
	*/
	cgn_set ( 0 );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _multmod ( void * data)
/******************************************************************************/
/*
  Purpose:
    MULTMOD carries out modular multiplication.
  Discussion: 
    This procedure returns 
  ( A * S ) mod M
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2013
  Author:
    Original Pascal version by Pierre L'Ecuyer, Serge Cote.
    C version by John Burkardt.
  Reference:
    Pierre LEcuyer, Serge Cote,
    Implementing a Random Number Package with Splitting Facilities,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 1, March 1991, pages 98-111.
  Parameters:
    Input, int A, S, M, the arguments.
    Output, int MULTMOD, the value of the product of A and S, 
    modulo M.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	int a = a_data[0]; 
	int s = a_data[1];
	int m = a_data[2];
	
	int a0;
	int a1;
	const int h = 32768;
	int k;
	int p;
	int q;
	int qh;
	int rh;
	
	if ( a <= 0 || m <= a || s <= 0 || m <= s)
	{
		result = INT_MAX;
		return &result;
	}
	
	if ( a < h )
	{
		a0 = a;
		p = 0;
	}
	else
	{
		a1 = a / h;
		a0 = a - h * a1;
		qh = m / h;
		rh = m - h * qh;
		
		if ( h <= a1 )
		{
			a1 = a1 - h;
			k = s / qh;
			p = h * ( s - k * qh ) - k * rh;
			
			while ( p < 0 )
				p += m;
		}
		else
			p = 0;
	
		if ( a1 != 0 )
		{
			q = m / a1;
			k = s / q;
			p = p - k * ( m - a1 * q );
			
			if ( 0 < p )
				p -= m;
	
			p += a1 * ( s - k * q );
	
			while ( p < 0 )
				p += m;
		}
		
		k = p / qh;
		p = h * ( p - k * qh ) - k * rh;
	
		while ( p < 0 )
			p += m;
	}
	
	if ( a0 != 0 )
	{
		q = m / a0;
		k = s / q;
		p -= k * ( m - a0 * q );
		
		if ( 0 < p )
			p -= m;
	
		p += a0 * ( s - k * q );
	
		while ( p < 0 )
			p += m;
	}
	
	result = p;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _set_initial_seed ( void * data)
/******************************************************************************/
/*
  Purpose:
    SET_INITIAL_SEED resets the initial seed and state for all generators.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 March 2013
  Author:
    Original Pascal version by Pierre L'Ecuyer, Serge Cote.
    C version by John Burkardt.
  Reference:
    Pierre LEcuyer, Serge Cote,
    Implementing a Random Number Package with Splitting Facilities,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 1, March 1991, pages 98-111.
  Parameters:
    Input, int IG1, IG2, the initial seed values 
    for the first generator.
    1 <= IG1 < 2147483563
    1 <= IG2 < 2147483399
*/
{
	int * const a_data = data;
	int ig1 = a_data[0];
	int ig2 = a_data[1];
	
	const int a1_vw = 2082007225;
	const int a2_vw = 784306273;
	dim_typ g;
	const int g_max = 32;
	int i;
	const int m1 = 2147483563;
	const int m2 = 2147483399;
	int t;
	
	if ( ig1 < 1 || m1 <= ig1 || ig2 < 1 || m2 <= ig2 )
		return NULL;
	/*
	Because INITIALIZE calls SET_INITIAL_SEED, it's not easy to correct
	the error that arises if SET_INITIAL_SEED is called before INITIALIZE.
	So don't bother trying.
	*/
	if ( ! initialized_get ( ) )
		return NULL;
	/*
	Set the initial seed, then initialize the first generator.
	*/
	g = 0;
	cgn_set ( g );
	ig_set ( g, ig1, ig2 );
	t = 0;
	init_generator ( t );
	/*
	Now do similar operations for the other generators.
	*/
	for ( g = 1; g < g_max; ++g )
	{
		cgn_set ( g );
		ig1 = multmod ( a1_vw, ig1, m1 );
		ig2 = multmod ( a2_vw, ig2, m2 );
		ig_set ( g, ig1, ig2 );
		init_generator ( t );
	}
	/*
	Now choose the first generator.
	*/
	g = 0;
	cgn_set ( g );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _antithetic_get ( void * data)
/******************************************************************************/
/*
  Purpose:
    ANTITHETIC_GET queries the antithetic value for a given generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2013
  Author:
    John Burkardt
  Parameters:
    Output, int ANTITHETIC_GET, is TRUE (1) if generator G is antithetic.
*/
{
	static bool result = 2;
	
	int value;
	antithetic_memory ( -1, &value );
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _antithetic_memory ( void * data)
/******************************************************************************/
/*
  Purpose:
    ANTITHETIC_MEMORY stores the antithetic value for all generators.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, the desired action.
    -1, get a value.
    0, initialize all values.
    1, set a value.
    Input/output, int UE.  For I = -1, VALUE is an output
    quantity.  If I = +1, then VALUE is an input quantity.
*/
{
	const ipi * const s_data = data;
	int i = s_data->a0;
	int * value = s_data->a1;	
	
	# define G_MAX 32
	
	static int a_save[G_MAX];
	int g;
	const int g_max = 32;
	dim_typ j;
	
	if ( i < 0 )
	{
		g = cgn_get ( );
		*value = a_save[g];
	}
	else if ( i == 0 )
		for ( j = 0; j < g_max; ++j)
		a_save[j] = 0;
	else if ( 0 < i )
	{
		g = cgn_get ( );
		a_save[g] = *value;
	}
	
	return NULL;
	# undef G_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _antithetic_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    ANTITHETIC_SET sets the antithetic value for a given generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2013
  Author:
    John Burkardt
  Parameters:
    Input, int VALUE, is TRUE (1) if generator G is to be antithetic.
*/
{
	int value = *(int *) data;
		
	antithetic_memory ( +1, &value );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _cgn_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    CGN_SET sets the current generator index.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int G, the current generator index.
    0 <= G <= 31.
*/
{
	int g = *(int *) data;
	
	cgn_memory ( +1, &g );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _initialized_get ( void * data)
/******************************************************************************/
/*
  Purpose:
    INITIALIZED_GET queries the INITIALIZED value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 March 2013
  Author:
    John Burkardt
  Parameters:
    Output, int INITIALIZED_GET, is TRUE (1) if the package has been initialized.
*/
{
	static bool result = 2;
	
	int value;
	initialized_memory ( -1, &value );
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _initialized_memory ( void * data)
/******************************************************************************/
/*
  Purpose:
    INITIALIZED_MEMORY stores the INITIALIZED value for the package.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, the desired action.
    -1, get the value.
    0, initialize the value.
    1, set the value.
    Input/output, int *INITIALIZED.  For I = -1, this is an output
    quantity.  If I = +1, this is an input quantity.  If I = 0, 
    this is ignored.
*/
{
	const ipi * const s_data = data;
	int i = s_data->a0;
	int * initialized = s_data->a1;	
	
	static int initialized_save = 0;
	
	if ( i < 0 )
		*initialized = initialized_save;
	else if ( i == 0 )
		initialized_save = 0;
	else if ( 0 < i )
		initialized_save = *initialized;
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _initialized_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    INITIALIZED_SET sets the INITIALIZED value true.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 March 2013
  Author:
    John Burkardt
  Parameters:
    None
*/
{
	int initialized = 1;
	initialized_memory ( +1, &initialized );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cg_get ( void * data)
/******************************************************************************/
/*
  Purpose:
    CG_GET queries the CG values for a given generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int G, the index of the generator.
    0 <= G <= 31.
    Output, int *CG1, *CG2, the CG values for generator G.
*/
{
	const i2pi * const s_data = data;
	int g = s_data->a0;
	int * cg1 = s_data->a1;
	int * cg2 = s_data->a2;
	
	cg_memory ( -1, g, cg1, cg2 );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cg_memory ( void * data)
/******************************************************************************/
/*
  Purpose:
    CG_MEMORY stores the CG values for all generators.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, the desired action.
    -1, get a value.
    0, initialize all values.
    1, set a value.
    Input, int G, for I = -1 or +1, the index of 
    the generator, with 0 <= G <= 31.
    Input/output, int *CG1, *CG2.  For I = -1, 
    these are output, for I = +1, these are input, for I = 0,
    these arguments are ignored.  When used, the arguments are
    old or new values of the CG parameter for generator G.
*/
{
	const _2i2pi * const s_data = data;
	int i = s_data->a0;
	int g = s_data->a1;
	int * cg1 = s_data->a2;
	int * cg2 = s_data->a3;
	
	# define G_MAX 32
	
	static int cg1_save[G_MAX];
	static int cg2_save[G_MAX];
	const int g_max = 32;
	dim_typ j;
	
	if ( g < 0 || g_max <= g )
		return NULL;
	
	if ( i < 0 )
	{
		*cg1 = cg1_save[g];
		*cg2 = cg2_save[g];
	}
	else if ( i == 0 )
		for ( j = 0; j < g_max; ++j )
			cg1_save[j] = cg2_save[j] = 0;
	else if ( 0 < i )
	{
		cg1_save[g] = *cg1;
		cg2_save[g] = *cg2;
	}
	
	return NULL;
	# undef G_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cg_set ( void * data)
/******************************************************************************/
/*
  Purpose:
    CG_SET sets the CG values for a given generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int G, the index of the generator.
    0 <= G <= 31.
    Input, int CG1, CG2, the CG values for generator G.
*/
{
	int * const a_data = data;
	int g = a_data[0]; 
	int cg1 = a_data[1];
	int cg2 = a_data[2];
	
	cg_memory ( +1, g, &cg1, &cg2 );
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cgn_get ( void * data)
/******************************************************************************/
/*
  Purpose:
    CGN_GET gets the current generator index
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2013
  Author:
    John Burkardt
  Parameters:
    Output, int CGN_GET, the current generator index.
*/
{
	static int result = INT_MAX;
	
	int g;
	cgn_memory ( -1, &g );
	result = g;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cgn_memory ( void * data)
/******************************************************************************/
/*
  Purpose:
    CGN_MEMORY stores the current generator index.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, the desired action.
    -1, get the value.
    0, initialize the value.
    1, set the value.
    Input/output, int *G.  For I = -1 or 0, this is output.
    For I = 1, this is input.
*/
{
	const ipi * const s_data = data;
	int i = s_data->a0;
	int * g = s_data->a1;	
	
	# define G_MAX 32
	
	static int g_save = 0;
	const int g_max = 32;
	int j;
	
	if ( i < 0 )
		*g = g_save;
	else if ( i == 0 )
	{
		g_save = 0;
		*g = g_save;
	}
	else if ( 0 < i )
	{
		if ( *g < 0 || g_max <= *g )
			return NULL;
		g_save = *g;
	}
	
	return NULL;
	# undef G_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_uni ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_UNI generates a random positive integer.
  Discussion:
    This procedure returns a random integer following a uniform distribution 
    over (1, 2147483562) using the current generator.
    The original name of this function was "random()", but this conflicts
    with a standard library function name in C.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 August 2013
  Author:
    Original Pascal version by Pierre L'Ecuyer, Serge Cote.
    C version by John Burkardt.
  Reference:
    Pierre LEcuyer, Serge Cote,
    Implementing a Random Number Package with Splitting Facilities,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 1, March 1991, pages 98-111.
  Parameters:
    Output, int I4_UNI, the random integer.
*/
{
	static int result = INT_MAX;
	
	const int a1 = 40014;
	const int a2 = 40692;
	int cg1;
	int cg2;
	int g;
	int k;
	const int m1 = 2147483563;
	const int m2 = 2147483399;
	int value;
	int z;
	/*
	Check whether the package must be initialized.
	*/
	if ( ! initialized_get ( ) )
		initialize ( );
	/*
	Get the current generator index.
	*/
	g = cgn_get ( );
	/*
	Retrieve the current seeds.
	*/
	cg_get ( g, &cg1, &cg2 );
	/*
	Update the seeds.
	*/
	k = cg1 / 53668;
	cg1 = a1 * ( cg1 - k * 53668 ) - k * 12211;
	
	if ( cg1 < 0 )
		cg1 += m1;
	
	k = cg2 / 52774;
	cg2 = a2 * ( cg2 - k * 52774 ) - k * 3791;
	
	if ( cg2 < 0 )
		cg2 += m2;
	/*
	Store the updated seeds.
	*/
	cg_set ( g, cg1, cg2 );
	/*
	Form the random integer.
	*/
	z = cg1 - cg2;
	
	if ( z < 1 )
		z += m1 - 1;
	/*
	If the generator is antithetic, reflect the value.
	*/
	
	if ( antithetic_get ( ) )
		z = m1 - z;
	
	result = value ? m1-z:z;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_uni_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_UNI_01 returns a uniform random double in [0,1].
  Discussion:
    This procedure returns a random floating point number from a uniform 
    distribution over (0,1), not including the endpoint values, using the
    current random number generator.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 August 2013
  Author:
    Original Pascal version by Pierre L'Ecuyer, Serge Cote.
    C version by John Burkardt.
  Reference:
    Pierre LEcuyer, Serge Cote,
    Implementing a Random Number Package with Splitting Facilities,
    ACM Transactions on Mathematical Software,
    Volume 17, Number 1, March 1991, pages 98-111.
  Parameters:
    Output, double R8_UNI_01, a uniform random value in [0,1].
*/
{
	static ityp result = MAX_VAL;
	
	/*
	Check whether the package must be initialized.
	*/
	if ( ! initialized_get ( ) )
		initialize ( );
	/*
	Get a random integer.
	*/
	/*
	Scale it to [0,1].
	*/
	
	result = ( ityp ) ( i4_uni ( ) ) * 4.656613057E-10;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_uniform_01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_UNIFORM_01_SAMPLE generates a uniform random deviate from [0,1].
  Discussion:
    This function should be the only way that the package accesses random
    numbers.
    Setting OPTION to 0 accesses the R8_UNI_01() function in RNGLIB,
    for which there are versions in various languages, which should result
    in the same values being returned.  This should be the only place in
    this library that accesses a function in RNGLIB.
    Setting OPTION to 1 in the C version calls the system random number
    generator function "rand()".
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 August 2013
  Author:
    John Burkardt.
  Parameters:
    Output, double R8_UNIFORM_01_SAMPLE, a random deviate.
*/
{
	static ityp result = MAX_VAL;
	
	result = r8_uni_01();
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_normal_01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_NORMAL_01_SAMPLE samples the standard normal probability distribution.
  Discussion:
    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.
    The Box-Muller method is used, which is efficient, but
    generates two values at a time.
    Typically, we would use one value and save the other for the next call.
    However, the fact that this function has saved memory makes it difficult
    to correctly handle cases where we want to re-initialize the code,
    or to run in parallel.  Therefore, we will instead use the first value
    and DISCARD the second.
    EFFICIENCY must defer to SIMPLICITY.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 July 2013
  Author:
    John Burkardt
  Parameters:
    Output, double R8_NORMAL_01_SAMPLE, a normally distributed random value.
*/
{
	static ityp result = MAX_VAL;
	
	result = sqrt ( -2.00 * log ( r8_uniform_01_sample ( ) ) ) * cos ( M_2TPI * r8_uniform_01_sample ( ) );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *    _r8_chi_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CHI_SAMPLE generates a Chi-Square random deviate.
  Discussion:
    This procedure generates a random deviate from the chi square distribution
    with DF degrees of freedom random variable.
    The algorithm exploits the relation between chisquare and gamma.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2013
  Author:
    Original FORTRAN77 version by Barry Brown, James Lovato.
    C version by John Burkardt.
  Parameters:
    Input, double DF, the degrees of freedom.
    0.0 < DF.
    Output, double R8_CHI_SAMPLE, a random deviate from the distribution.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp df = *(ityp *) data;
	
	result = df <= 0.00 ? MAX_VAL : 2.00 * r8_gamma_sample ( 1.00, df / 2.00 );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _bartlett_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    BARTLETT_SAMPLE samples the Bartlett distribution.
  Discussion:
    If the matrix T is sampled from the Bartlett distribution, then
    the matrix W = T' * T is a sample from the Wishart distribution.
    This function requires functions from the PDFLIB and RNGLIB libraries.
    The "initialize()" function from RNGLIB must be called before using
    this function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 July 2013
  Author:
    John Burkardt
  Reference:
    Patrick Odell, Alan Feiveson,
    A numerical procedure to generate a sample covariance matrix,
    Journal of the American Statistical Association,
    Volume 61, Number 313, March 1966, pages 199-203.
    Stanley Sawyer,
    Wishart Distributions and Inverse-Wishart Sampling,
    Washington University,
    30 April 2007, 12 pages.
  Parameters:
    Input, int M, the order of the matrix.
    Input, int DF, the number of degrees of freedom.
    M <= DF.
    Input, double SIGMA[M*M], the covariance matrix, which should be
    a symmetric positive definite matrix.
    Output, double BARTLETT_SAMPLE[M*M], the sample matrix from
    the Bartlett distribution.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ df = s_data->a1;
	ityp * sigma = s_data->a2;
	
    bool flag;
    ityp *r;
    ityp *t;
    ityp *tu;

    if ( df < m )
        return NULL;
    /*
    Get the upper triangular Cholesky factor of SIGMA.
    */
    r = r8mat_cholesky_factor_upper ( m, sigma, &flag );

    if ( flag )
        return NULL;
    /*
    Sample the unit Bartlett distribution.
    */
    tu = bartlett_unit_sample ( m, df );
    /*
    Construct the matrix T = TU * R.
    */
    t = r8mat_mm_new ( m, m, m, tu, r );
    /*
    Free memory.
    */
    free ( r );
    free ( tu );

    return t;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _bartlett_unit_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    BARTLETT_UNIT_SAMPLE samples the unit Bartlett distribution.
  Discussion:
    If the matrix T is sampled from the unit Bartlett distribution, then
    the matrix W = T' * T is a sample from the unit Wishart distribution.
    This function requires functions from the PDFLIB and RNGLIB libraries.
    The "initialize()" function from RNGLIB must be called before using
    this function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 October 2013
  Author:
    John Burkardt
  Reference:
    Patrick Odell, Alan Feiveson,
    A numerical procedure to generate a sample covariance matrix,
    Journal of the American Statistical Association,
    Volume 61, Number 313, March 1966, pages 199-203.
    Stanley Sawyer,
    Wishart Distributions and Inverse-Wishart Sampling,
    Washington University,
    30 April 2007, 12 pages.
  Parameters:
    Input, int M, the order of the matrix.
    Input, int DF, the number of degrees of freedom.
    M <= DF.
    Output, double BARTLETT_UNIT_SAMPLE[M*M], the sample matrix from the
    unit Bartlett distribution.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ df = a_data[1];
	
    ityp df_chi;
    dim_typ i, j;
    ityp *t;

    if ( df < m )
        return NULL;

    t = ( ityp * ) malloc ( m * m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i)
    {
        for ( j = 0; j < i; ++j)
            t[i+j*m] = 0.00;
        df_chi = ( ityp ) ( df - i );
        t[i+i*m] = sqrt ( r8_chi_sample ( df_chi ) );
        for ( j = i + 1; j < m; ++j )
            t[i+j*m] = r8_normal_01_sample ( );
    }

    return t;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wishart_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    WISHART_SAMPLE samples the Wishart distribution.
  Discussion:
    This function requires functions from the PDFLIB and RNGLIB libraries.
    The "initialize()" function from RNGLIB must be called before using
    this function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 July 2013
  Author:
    John Burkardt
  Reference:
    Patrick Odell, Alan Feiveson,
    A numerical procedure to generate a sample covariance matrix,
    Journal of the American Statistical Association,
    Volume 61, Number 313, March 1966, pages 199-203.
    Stanley Sawyer,
    Wishart Distributions and Inverse-Wishart Sampling,
    Washington University,
    30 April 2007, 12 pages.
  Parameters:
    Input, int M, the order of the matrix.
    Input, int DF, the number of degrees of freedom.
    M <= DF.
    Input, double SIGMA[M*M], the covariance matrix, which should be
    a symmetric positive definite matrix.
    Output, double WISHART_SAMPLE[M*M], the sample matrix from
    the Wishart distribution.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ df = s_data->a1;
	ityp * sigma = s_data->a2;
	
    ityp *a;
    ityp *au;
    ityp *aur;
    bool flag;
    ityp *r;

    if ( df < m )
        return NULL;
    /*
    Get R, the upper triangular Cholesky factor of SIGMA.
    */
    r = r8mat_cholesky_factor_upper ( m, sigma, &flag );

    if ( flag )
        return NULL;
    /*
    Get AU, a sample from the unit Wishart distribution.
    */
    au = wishart_unit_sample ( m, df );
    /*
    Construct the matrix A = R' * AU * R.
    */
    aur = r8mat_mm_new ( m, m, m, au, r );
    a = r8mat_mtm_new ( m, m, m, r, aur );
    /*
    Free memory.
    */
    free ( au );
    free ( aur );
    free ( r );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wishart_sample_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    WISHART_SAMPLE_INVERSE returns the inverse of a sample Wishart matrix.
  Discussion:
    This function requires functions from the PDFLIB and RNGLIB libraries.
    The "initialize()" function from RNGLIB must be called before using
    this function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 October 2013
  Author:
    John Burkardt
  Reference:
    Patrick Odell, Alan Feiveson,
    A numerical procedure to generate a sample covariance matrix,
    Journal of the American Statistical Association,
    Volume 61, Number 313, March 1966, pages 199-203.
    Stanley Sawyer,
    Wishart Distributions and Inverse-Wishart Sampling,
    Washington University,
    30 April 2007, 12 pages.
  Parameters:
    Input, int M, the order of the matrix.
    Input, int DF, the number of degrees of freedom.
    M <= DF.
    Input, double SIGMA[M*M], the covariance matrix, which should be
    a symmetric positive definite matrix.
    Output, double WISHART_SAMPLE[M*M], the inverse of a sample matrix from
    the Wishart distribution.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ df = s_data->a1;
	ityp * sigma = s_data->a2;
	
    ityp *a;
    bool flag;
    ityp *r;
    ityp *s;
    ityp *ua;
    ityp *uas;

    if ( df < m )
        return NULL;
    /*
    Get R, the upper triangular Cholesky factor of SIGMA.
    */
    r = r8mat_cholesky_factor_upper ( m, sigma, &flag );

    if ( flag )
        return NULL;
    /*
    Get S, the inverse of R.
    */
    s = r8ut_inverse ( m, r );
    /*
    Get UA, the inverse of a sample from the unit Wishart distribution.
    */
    ua = wishart_unit_sample_inverse ( m, df );
    /*
    Construct the matrix A = S * UA * S'.
    */
    uas = r8mat_mmt_new ( m, m, m, ua, s );
    a = r8mat_mm_new ( m, m, m, s, uas );
    /*
    Free memory.
    */
    free ( r );
    free ( s );
    free ( ua );
    free ( uas );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wishart_unit_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    WISHART_UNIT_SAMPLE samples the unit Wishart distribution.
  Discussion:
    This function requires functions from the PDFLIB and RNGLIB libraries.
    The "initialize()" function from RNGLIB must be called before using
    this function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 October 2013
  Author:
    John Burkardt
  Reference:
    Patrick Odell, Alan Feiveson,
    A numerical procedure to generate a sample covariance matrix,
    Journal of the American Statistical Association,
    Volume 61, Number 313, March 1966, pages 199-203.
    Stanley Sawyer,
    Wishart Distributions and Inverse-Wishart Sampling,
    Washington University,
    30 April 2007, 12 pages.
  Parameters:
    Input, int M, the order of the matrix.
    Input, int DF, the number of degrees of freedom.
    M <= DF.
    Output, double WISHART_UNIT_SAMPLE[M*M], the sample matrix from the
    unit Wishart distribution.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ df = a_data[1];
	
    ityp *a;
    ityp *c;
    ityp df_chi;
    dim_typ i, j;

    if ( df < m )
        return NULL;

    c = ( ityp * ) malloc ( m * m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
    {
        for ( j = 0; j < i; ++j )
            c[i+j*m] = 0.00;
        df_chi = ( ityp ) ( df - i );
        c[i+i*m] = sqrt ( r8_chi_sample ( df_chi ) );
        for ( j = i + 1; j < m; ++j )
            c[i+j*m] = r8_normal_01_sample ( );
    }

    a = r8mat_mtm_new ( m, m, m, c, c );
    /*
    Free memory.
    */
    free ( c );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wishart_unit_sample_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    WISHART_UNIT_SAMPLE_INVERSE inverts a unit Wishart sample matrix.
  Discussion:
    This function requires functions from the PDFLIB and RNGLIB libraries.
    The "initialize()" function from RNGLIB must be called before using
    this function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 October 2013
  Author:
    John Burkardt
  Reference:
    Patrick Odell, Alan Feiveson,
    A numerical procedure to generate a sample covariance matrix,
    Journal of the American Statistical Association,
    Volume 61, Number 313, March 1966, pages 199-203.
    Stanley Sawyer,
    Wishart Distributions and Inverse-Wishart Sampling,
    Washington University,
    30 April 2007, 12 pages.
  Parameters:
    Input, int M, the order of the matrix.
    Input, int DF, the number of degrees of freedom.
    M <= DF.
    Output, double WISHART_UNIT_SAMPLE_INVERSE[M*M], the inverse of a
    sample matrix from the unit Wishart distribution.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ df = a_data[1];
	
    ityp *a;
    ityp *b;
    ityp *c;
    ityp df_chi;
    dim_typ i, j;

    if ( df < m )
        return NULL;

    c = ( ityp * ) malloc ( m * m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
    {
        for ( j = 0; j < i; ++j )
            c[i+j*m] = 0.00;
        df_chi = ( ityp ) ( df - i );
        c[i+i*m] = sqrt ( r8_chi_sample ( df_chi ) );
        for ( j = i + 1; j < m; ++j )
            c[i+j*m] = r8_normal_01_sample ( );
    }
    /*
    Compute B, the inverse of C.
    */
    b = r8ut_inverse ( m, c );
    /*
    The inverse of the Wishart sample matrix C'*C is inv(C) * C'.
    */
    a = r8mat_mmt_new ( m, m, m, b, b );
    /*
    Free memory.
    */
    free ( b );
    free ( c );

    return a;
}

#endif
