#ifndef __DISABLEDEEP_POLPAK

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fibonacci_floor ( void * data)
/******************************************************************************/
/*
  Purpose:
    FIBONACCI_FLOOR returns the largest Fibonacci number less or equal to N.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the positive integer whose Fibonacci "floor" is desired.
    Output, int *F, the largest Fibonacci number less than or equal to N.
    Output, int *I, the index of the F.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * f = s_data->a1;
	int * i = s_data->a2;
	
	if ( n == 0 )
		*i = *f = 0; 
	else
	{
		*i = ( int ) ( log ( 0.50 * ( ityp ) ( (n<<1) + 1 ) * sqrt ( 5.00 ) )/ log ( 0.50 * ( 1.00 + sqrt ( 5.00 ) ) ) );
		*f = fibonacci_direct ( *i );
		
		if ( n < *f )
		{
			-- *i;
			*f = fibonacci_direct ( *i );
		}
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _legendre_associated_normalized ( void * data)
/******************************************************************************/
/*
  Purpose:
    LEGENDRE_ASSOCIATED_NORMALIZED: normalized associated Legendre functions.
  Discussion:
    The unnormalized associated Legendre functions P_N^M(X) have
    the property that
      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX
      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
    By dividing the function by the square root of this term,
    the normalized associated Legendre functions have norm 1.
    However, we plan to use these functions to build spherical
    harmonics, so we use a slightly different normalization factor of
      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) )
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    07 March 2005
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input, int N, the maximum first index of the Legendre
    function, which must be at least 0.
    Input, int M, the second index of the Legendre function,
    which must be at least 0, and no greater than N.
    Input, double X, the point at which the function is to be
    evaluated.  X must satisfy -1 <= X <= 1.
    Output, double CX[N+1], the values of the first N+1 function.
*/
{
	const _2dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * cx = s_data->a2;
	const register ityp x = s_data->a3;
	
	ityp factor;
	dim_typ i;
	int mm;
	ityp somx2;
	
	if ( m ==0 || n<m || x <-1.00 || 1.00 < x)
		return NULL;
	
	for ( i = 0; i <= m - 1; ++i )
		cx[i] = 0.00;
	cx[m] = 1.00;
	
	somx2 = sqrt ( 1.00 - x * x );
	
	factor = 1.00;
	for ( i = 1; i <= m; ++i )
	{
		cx[m] = -cx[m] * factor * somx2;
		factor += 2.00;
	}
	
	if ( m + 1 <= n )
		cx[m+1] = x * ( ityp ) ( (m<<1) + 1 ) * cx[m];
	
	for ( i = m + 2; i <= n; ++i )
		cx[i] = ( ( ityp ) ( (i<<1)     - 1 ) * x * cx[i-1] + ( ityp ) (   - i - m + 1 )     * cx[i-2] ) / ( ityp ) (     i - m     );
	/*
	Normalization.
	*/
	for ( mm = m; mm <= n; ++mm )
	{
		factor = sqrt ( ( ( ityp ) ( (mm<<1) + 1 ) * r8_factorial ( mm - m ) ) / ( 4.00 * M_PI * r8_factorial ( mm + m ) ) );
		cx[mm] = cx[mm] * factor;
	}
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_is_triangular ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_IS_TRIANGULAR determines whether an I4 is triangular.
  Discussion:
    The N-th triangular number is equal to the sum of the first
    N integers.
  First Values:
    Index  Value
     0      0
     1      1
     2      3
     3      6
     4     10
     5     15
     6     21
     7     28
     8     36
     9     45
    10     55
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    19 February 2003
  Author:
    John Burkardt
  Parameters:
    Input, int I, the integer to be checked.
    Output, int I4_IS_TRIANGULAR, is TRUE if I is triangular.
*/
{
	static bool result = 2;
	
	const register dim_typ i = *(dim_typ *) data;
	
	int j, k;
	
	if(i)
	{
		i4_to_triangle ( i, &j, &k );
		result = j == k;
		return &result;
	}
	
	result = i;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cheby_t_poly ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBY_T_POLY evaluates Chebyshev polynomials T(n,x).
  Discussion:
    Chebyshev polynomials are useful as a basis for representing the
    approximation of functions since they are well conditioned, in the sense
    that in the interval [-1,1] they each have maximum absolute value 1.
    Hence an error in the value of a coefficient of the approximation, of
    size epsilon, is exactly reflected in an error of size epsilon between
    the computed approximation and the theoretical approximation.
    Typical usage is as follows, where we assume for the moment
    that the interval of approximation is [-1,1].  The value
    of N is chosen, the highest polynomial to be used in the
    approximation.  Then the function to be approximated is
    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
    Chebyshev polynomial.  Let these values be denoted by F(XJ).
    The coefficients of the approximation are now defined by
      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
    except that C(0) is given a value which is half that assigned
    to it by the above formula,
    and the representation is
    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
    Now note that, again because of the fact that the Chebyshev polynomials
    have maximum absolute value 1, if the higher order terms of the
    coefficients C are small, then we have the option of truncating
    the approximation by dropping these terms, and we will have an
    exact value for maximum perturbation to the approximation that
    this will cause.
    It should be noted that typically the error in approximation
    is dominated by the first neglected basis function (some multiple of
    T(N+1,X) in the example above).  If this term were the exact error,
    then we would have found the minimax polynomial, the approximating
    polynomial of smallest maximum deviation from the original function.
    The minimax polynomial is hard to compute, and another important
    feature of the Chebyshev approximation is that it tends to behave
    like the minimax polynomial while being easy to compute.
    To evaluate a sum like
      sum ( 0 <= J <= N ) C(J) T(J,X),
    Clenshaw's recurrence formula is recommended instead of computing the
    polynomial values, forming the products and summing.
    Assuming that the coefficients C(J) have been computed
    for J = 0 to N, then the coefficients of the representation of the
    indefinite integral of the function may be computed by
      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1,
    with
      C(N+1)=0
      B(0) arbitrary.
    Also, the coefficients of the representation of the derivative of the
    function may be computed by:
      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0,
    with
      D(N+1) = D(N)=0.
    Some of the above may have to adjusted because of the irregularity of C(0).
    The formula is:
      T(N,X) = COS(N*acos(X))
  Differential equation:
 (1-X*X) Y'' - X Y' + N N Y = 0
  First terms:
    T(0,X) =  1
    T(1,X) =  1 X
    T(2,X) =  2 X^2 -   1
    T(3,X) =  4 X^3 -   3 X
    T(4,X) =  8 X^4 -   8 X^2 +  1
    T(5,X) = 16 X^5 -  20 X^3 +  5 X
    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
  Inequality:
    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
  Orthogonality:
    For integration over [-1,1] with weight
      W(X) = 1 / sqrt(1-X*X),
    if we write the inner product of T(I,X) and T(J,X) as
      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
    then the result is:
      0 if I /= J
      PI/2 if I == J /= 0
      PI if I == J == 0
    A discrete orthogonality relation is also satisfied at each of
    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
                              = 0 if I /= J
                              = N/2 if I == J /= 0
                              = N if I == J == 0
  Recursion:
    T(0,X) = 1,
    T(1,X) = X,
    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
  Special values:
    T(N,1) = 1
    T(N,-1) = (-1)^N
    T(2N,0) = (-1)^N
    T(2N+1,0) = 0
    T(N,X) = (-1)^N * T(N,-X)
  Zeroes:
    M-th zero of T(N,X) is cos((2*M-1)*PI/(2*N)), M = 1 to N
  Extrema:
    M-th extremum of T(N,X) is cos(PI*M/N), M = 0 to N
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    28 March 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest polynomial to compute.
    Input, double X[M], the evaluation points.
    Output, double CHEBY_T_POLY[M*(N+1)], the values of the Chebyshev polynomials.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	
	dim_typ i, j;
	ityp *v;
	
	if ( n < 0 )
		return NULL;
	
	v = ( ityp * ) malloc ( m * ( n + 1 ) * sizeof ( ityp ) );
	
	for ( i = 0; i < m; ++i)
		v[i+0*m] = 1.00;
	if ( n < 1 )
		return v;
	
	for ( i = 0; i < m; ++i )
		v[i+1*m] = x[i];
	for ( j = 2; j <= n; ++j )
		for ( i = 0; i < m; ++i )
		v[i+j*m] = 2.00 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
	
	return v;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _agud ( void * data)
/******************************************************************************/
/*
  Purpose:
    AGUD evaluates the inverse Gudermannian function.
  Definition:
    The Gudermannian function relates the hyperbolic and trigonomentric
    functions.  For any argument X, there is a corresponding value
    G so that
      SINH(X) = TAN(G).
    This value G(X) is called the Gudermannian of X.  The inverse
    Gudermannian function is given as input a value G and computes
    the corresponding value X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, double G, the value of the Gudermannian.
    Output, double AGUD, the argument of the Gudermannian.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp g = *(ityp *) data;
	
	result = log ( tan ( 0.25 * M_PI + 0.50 * g ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _align_enum ( void * data)
/******************************************************************************/
/*
  Purpose:
    ALIGN_ENUM counts the alignments of two sequences of M and N elements.
  Discussion:
    We assume that we have sequences A and B of M and N characters each.
    An alignment of the two sequences is a rule matching corresponding
    elements of one sequence to another.  Some elements of either sequence
    can be matched to a null element.  If A(I1) and A(I2) are matched
    to B(J1) and B(J2), and I1 < I2, then it must be the case that J1 < J2.
    The 5 alignments of a sequence of 1 to a sequence of 2 are:
          _1_   _2_   __3__   __4__   __5__
      A:  1 -   - 1   - 1 -   - - 1   1 - -
      B:  1 2   1 2   1 - 2   1 2 -   - 1 2
    The formula is:
      F(0,0) = 1
      F(1,0) = 1
      F(0,1) = 1
      F(M,N) = F(M-1,N) + F(M-1,N-1) + F(M,N-1)
    To compute F(M,N), it is not necessary to keep an M+1 by N+1
    array in memory.  A vector of length N will do.
    F(N,N) is approximately ( 1 + sqrt(2) )**(2*N+1) / sqrt ( N )
  Example:
    The initial portion of the table is:
  M/N   0    1    2    3    4       5       6       7       8       9      10
   0    1    1    1    1    1       1       1       1       1       1       1
   1    1    3    5    7    9      11      13      15      17      19      21
   2    1    5   13   25   41      61      85     113     145     181     221
   3    1    7   25   63  129     231     377     575     833    1159    1561
   4    1    9   41  129  321     681    1289    2241    3649    5641    8361
   5    1   11   61  231  681    1683    3653    7183   13073   22363   36365
   6    1   13   85  377 1289    3653    8989   19825   40081   75517  134245
   7    1   15  113  575 2241    7183   19825   48639  108545  224143  433905
   8    1   17  145  833 3649   13073   40081  108545  265729  598417 1256465
   9    1   19  181 1159 5641   22363   75517  224143  598417 1462563 3317445
  10    1   21  221 1561 8361   36365  134245  433905 1256465 3317445 8097453
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2003
  Author:
    John Burkardt
  Reference:
    Michael Waterman,
    Introduction to Computational Biology,
    Chapman and Hall, 1995, pages 186-190.
  Parameters:
    Input, int M, N, the number of elements of the two sequences.
    Output, int ALIGN_ENUM, the number of possible alignments of the
    sequences.
*/
{
	static dim_typ result = USHRT_MAX;
	
	int * const a_data = data;
	const register int m = a_data[0];
	const register int n = a_data[1];
	
    int *fi;
    int fim1j;
    int fim1jm1;
    int i;
    int j;
    int value;

    if ( m < 0 || n < 0 )
    {
    	result = 0;
        return &result;
    }
    else if ( m == 0 || n == 0 )
    {
    	result = 1;
        return &result;
    }

    fi = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
    for ( i = 0; i <= n; ++i )
        fi[i] = 1;

    for ( i = 1; i <= m; ++i )
    {
        fim1jm1 = 1;

        for ( j = 1; j <= n; ++j )
        {

        fim1j = fi[j];
        fi[j] += fi[j-1] + fim1jm1;
        fim1jm1 = fim1j;

        }
    }

    value = fi[n];
    free ( fi );

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _benford ( void * data)
/******************************************************************************/
/*
  Purpose:
    BENFORD returns the Benford probability of one or more significant digits.
  Discussion:
    Benford's law is an empirical formula explaining the observed
    distribution of initial digits in lists culled from newspapers,
    tax forms, stock market prices, and so on.  It predicts the observed
    high frequency of the initial digit 1, for instance.
    Note that the probabilities of digits 1 through 9 are guaranteed
    to add up to 1, since
      LOG10 ( 2/1 ) + LOG10 ( 3/2) + LOG10 ( 4/3 ) + ... + LOG10 ( 10/9 )
      = LOG10 ( 2/1 * 3/2 * 4/3 * ... * 10/9 ) = LOG10 ( 10 ) = 1.
    The formula is:
      Prob ( First significant digits are IVAL ) =
        LOG10 ( ( IVAL + 1 ) / IVAL ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Frank Benford,
    The Law of Anomalous Numbers,
    Proceedings of the American Philosophical Society,
    Volume 78, pages 551-572, 1938.
    T P Hill,
    The First Digit Phenomenon,
    American Scientist,
    Volume 86, July/August 1998, pages 358 - 363.
    R Raimi,
    The Peculiar Distribution of First Digits,
    Scientific American,
    December 1969, pages 109-119.
  Parameters:
    Input, int IVAL, the string of significant digits to be checked.
    If IVAL is 1, then we are asking for the Benford probability that
    a value will have first digit 1.  If IVAL is 123, we are asking for
    the probability that the first three digits will be 123, and so on.
    Note that IVAL must not be 0 or negative.
    Output, double BENFORD, the Benford probability that an item taken
    from a real world distribution will have the initial digits IVAL.
*/
{
	static ityp result = MAX_VAL;
	
	const register int ival = *(int *) data;
	
	result = ival<=0 ? MAX_VAL:log10 ( ( ityp ) ( ival + 1 ) / ( ityp ) ival  );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bernoulli_number ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNOULLI_NUMBER computes the value of the Bernoulli numbers B(0) through B(N).
  Discussion:
    The Bernoulli numbers are rational.
    If we define the sum of the M-th powers of the first N integers as:
      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
    and let C(I,J) be the combinatorial coefficient:
      C(I,J) = I! / ( ( I - J )! * J! )
    then the Bernoulli numbers B(J) satisfy:
      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)**(M+1-J)
  First values:
   B0  1                   =         1.00000000000
   B1 -1/2                 =        -0.50000000000
   B2  1/6                 =         1.66666666666
   B3  0                   =         0
   B4 -1/30                =        -0.03333333333
   B5  0                   =         0
   B6  1/42                =         0.02380952380
   B7  0                   =         0
   B8 -1/30                =        -0.03333333333
   B9  0                   =         0
  B10  5/66                =         0.07575757575
  B11  0                   =         0
  B12 -691/2730            =        -0.25311355311
  B13  0                   =         0
  B14  7/6                 =         1.16666666666
  B15  0                   =         0
  B16 -3617/510            =        -7.09215686274
  B17  0                   =         0
  B18  43867/798           =        54.97117794486
  B19  0                   =         0
  B20 -174611/330          =      -529.12424242424
  B21  0                   =         0
  B22  854,513/138         =      6192.123
  B23  0                   =         0
  B24 -236364091/2730      =    -86580.257
  B25  0                   =         0
  B26  8553103/6           =   1425517.16666
  B27  0                   =         0
  B28 -23749461029/870     = -27298231.0678
  B29  0                   =         0
  B30  8615841276005/14322 = 601580873.901
  Recursion:
    With C(N+1,K) denoting the standard binomial coefficient,
    B(0) = 1.0
    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
  Warning:
    This recursion, which is used in this routine, rapidly results
    in significant errors.
  Special Values:
    Except for B(1), all Bernoulli numbers of odd index are 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the highest Bernoulli number to compute.
    Output, double B[N+1], B(I) contains the I-th Bernoulli number.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * b = s_data->a1;
	
    ityp b_sum;
    int *c;
    dim_typ i, j;

    if ( n < 1)
        return NULL;

    b[0] = 1.00;
    b[1] = -0.50;

    c = ( int * ) malloc ( ( n + 2 ) * sizeof ( int ) );
    c[0] = c[2] = 1;
    c[1] = 2;

    for ( i = 2; i <= n; ++i )
    {
        comb_row_next ( i + 1, c );

        if ( ( i % 2 ) == 1 )
            b[i] = 0.00;
        else
        {
            b_sum = 0.00;
            for ( j = 0; j <= i-1; ++j )
                b_sum += b[j] * ( ityp ) ( c[j] );
            b[i] = - b_sum / ( ityp ) ( c[i] );
        }
    }

    free ( c );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _bernoulli_number2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNOULLI_NUMBER2 evaluates the Bernoulli numbers.
  Discussion:
    The Bernoulli numbers are rational.
    If we define the sum of the M-th powers of the first N integers as:
      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
    and let C(I,J) be the combinatorial coefficient:
      C(I,J) = I! / ( ( I - J )! * J! )
    then the Bernoulli numbers B(J) satisfy:
      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)**(M+1-J)
    Note that the Bernoulli numbers grow rapidly.  Bernoulli number
    62 is probably the last that can be computed on the VAX without
    overflow.
    A different method than that used in BERN is employed.
  First values:
   B0  1                   =         1.00000000000
   B1 -1/2                 =        -0.50000000000
   B2  1/6                 =         1.66666666666
   B3  0                   =         0
   B4 -1/30                =        -0.03333333333
   B5  0                   =         0
   B6  1/42                =         0.02380952380
   B7  0                   =         0
   B8 -1/30                =        -0.03333333333
   B9  0                   =         0
  B10  5/66                =         0.07575757575
  B11  0                   =         0
  B12 -691/2730            =        -0.25311355311
  B13  0                   =         0
  B14  7/6                 =         1.16666666666
  B15  0                   =         0
  B16 -3617/510            =        -7.09215686274
  B17  0                   =         0
  B18  43867/798           =        54.97117794486
  B19  0                   =         0
  B20 -174611/330          =      -529.12424242424
  B21  0                   =         0
  B22  854,513/138         =      6192.123
  B23  0                   =         0
  B24 -236364091/2730      =    -86580.257
  B25  0                   =         0
  B26  8553103/6           =   1425517.16666
  B27  0                   =         0
  B28 -23749461029/870     = -27298231.0678
  B29  0                   =         0
  B30  8615841276005/14322 = 601580873.901
  Recursion:
    With C(N+1,K) denoting the standard binomial coefficient,
    B(0) = 1.0
    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
  Special Values:
    Except for B(1), all Bernoulli numbers of odd index are 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the highest order Bernoulli number to compute.
    Output, double B[N+1], the requested Bernoulli numbers.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * b = s_data->a1;
	
    ityp altpi;
    dim_typ i;
    dim_typ k;
    dim_typ kmax = 400;
    ityp sgn;
    ityp sum2;
    ityp t;
    ityp term;
    ityp tol = 1.0E-06;

    if ( n < 2)
        return NULL;

    b[0] = 1.00;
    b[1] = -0.50;

    altpi = log ( M_2TPI );
    /*
    Initial estimates for B(I), I = 2 to N
    */
    b[2] = log ( 2.00 );

    for ( i = 3; i <= n; ++i )
        b[i]= 0.00 + (( i % 2 ) !=1)*log ( ( ityp ) ( i * ( i - 1 ) ) ) + b[i-2];

    b[2] = 1.00 / 6.00;

    if ( n <= 3 )
        return NULL;

    b[4] = -1.00 / 30.00;
    sgn = -1.00;

    for ( i = 6; i <= n; i += 2 )
    {
        sgn *= -1;
        t = 2.00 * sgn * exp ( b[i] - ( ityp ) ( i ) * altpi );
        sum2 = 1.00;

        for ( k = 2; k <= kmax; ++k )
        {
            term = pow ( ( ityp ) k, ( ityp ) -i );
            sum2 += term;

            if ( term <= tol * sum2 )
                break;

        }

        b[i] = t * sum2;

    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bernoulli_number3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNOULLI_NUMBER3 computes the value of the Bernoulli number B(N).
  Discussion:
    The Bernoulli numbers are rational.
    If we define the sum of the M-th powers of the first N integers as:
      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
    and let C(I,J) be the combinatorial coefficient:
      C(I,J) = I! / ( ( I - J )! * J! )
    then the Bernoulli numbers B(J) satisfy:
      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M )
        C(M+1,J) B(J) * (N+1)**(M+1-J)
  First values:
     B0  1                   =         1.00000000000
     B1 -1/2                 =        -0.50000000000
     B2  1/6                 =         1.66666666666
     B3  0                   =         0
     B4 -1/30                =        -0.03333333333
     B5  0                   =         0
     B6  1/42                =         0.02380952380
     B7  0                   =         0
     B8 -1/30                =        -0.03333333333
     B9  0                   =         0
    B10  5/66                =         0.07575757575
    B11  0                   =         0
    B12 -691/2730            =        -0.25311355311
    B13  0                   =         0
    B14  7/6                 =         1.16666666666
    B15  0                   =         0
    B16 -3617/510            =        -7.09215686274
    B17  0                   =         0
    B18  43867/798           =        54.97117794486
    B19  0                   =         0
    B20 -174611/330          =      -529.12424242424
    B21  0                   =         0
    B22  854513/138          =      6192.123
    B23  0                   =         0
    B24 -236364091/2730      =    -86580.257
    B25  0                   =         0
    B26  8553103/6           =   1425517.16666
    B27  0                   =         0
    B28 -23749461029/870     = -27298231.0678
    B29  0                   =         0
    B30  8615841276005/14322 = 601580873.901
  Recursion:
    With C(N+1,K) denoting the standard binomial coefficient,
    B(0) = 1.0
    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
  Special Values:
    Except for B(1), all Bernoulli numbers of odd index are 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 February 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the Bernoulli number to compute.
    Output, double BERNOULLI_NUMBER3, the desired Bernoulli number.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    dim_typ itmax = 1000;
    ityp sum2;
    ityp term;
    ityp tol = 5.0E-07;
    ityp value;

    if ( n == 0 )
        value = 1.00;
    else if ( n == 1 )
        value = -0.50;
    else if ( n == 2 )
        value = 1.00 / 6.00;
    else if ( ( n % 2 ) == 1 )
        value = 0.00;
    else
    {
        sum2 = 0.00;

        for ( i = 1; i <= itmax; ++i )
        {
            term = 1.00 / pow ( ( ityp ) i, n );
            sum2 += term;

            if ( fabs ( term ) < tol || fabs ( term ) < tol * fabs ( sum2 ) )
                break;

        }

        value = 2.00 * sum2 * r8_factorial ( n )/ pow ( ( M_2TPI ), n );

        if ( ( n % 4 ) == 0 )
            value *= -1;

    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bernoulli_poly ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNOULLI_POLY evaluates the Bernoulli polynomial of order N at X.
  Discussion:
    Thanks to Bart Vandewoestyne for pointing out an error in the previous
    documentation, 31 January 2008.
    Special values of the Bernoulli polynomial include:
      B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
      B'(N,X) = N * B(N-1,X)
      B(N,X+1) - B(N,X) = N * X^(N-1)
      B(N,X) = (-1)**N * B(N,1-X)
    A formula for the Bernoulli polynomial in terms of the Bernoulli
    numbers is:
      B(N,X) = sum ( 0 <= K <= N ) B(K) * C(N,K) * X^(N-K)
    The first few polynomials include:
      B(0,X) = 1
      B(1,X) = X    - 1/2
      B(2,X) = X^2 -   X      +  1/6
      B(3,X) = X^3 - 3/2*X^2 +  1/2*X
      B(4,X) = X^4 - 2*X^3   +      X^2 - 1/30
      B(5,X) = X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
      B(6,X) = X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
      B(7,X) = X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
      B(8,X) = X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the Bernoulli polynomial to
    be evaluated.  N must be 0 or greater.
    Input, double X, the value of X at which the polynomial is to
    be evaluated.
    Output, double BERNOULLI_POLY, the value of B(N,X).
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
    int *c;
    dim_typ i;
    ityp value;
    ityp *work;

    work = ( ityp * ) malloc ( ( n + 1 ) * sizeof ( ityp ) );
    bernoulli_number ( n, work );
    /*
    Get row N of Pascal's triangle.
    */
    c = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
    for ( i = 0; i <= n; ++i )
        comb_row_next ( n, c );

    value = 1.00;
    for ( i = 1; i <= n; ++i )
        value = value * x + work[i] * ( ityp ) c[i];

    free ( c );
    free ( work );

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bernoulli_poly2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    BERNOULLI_POLY2 evaluates the N-th Bernoulli polynomial at X.
  Discussion:
    Thanks to Bart Vandewoestyne for pointing out an error in the previous
    documentation, 31 January 2008.
    Special values of the Bernoulli polynomial include:
      B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
      B'(N,X) = N * B(N-1,X)
      B(N,X+1) - B(N,X) = N * X^(N-1)
      B(N,X) = (-1)**N * B(N,1-X)
    A formula for the Bernoulli polynomial in terms of the Bernoulli
    numbers is:
      B(N,X) = sum ( 0 <= K <= N ) B(K) * C(N,K) * X^(N-K)
    The first few polynomials include:
      B(0,X) = 1
      B(1,X) = X    - 1/2
      B(2,X) = X^2 -   X      +  1/6
      B(3,X) = X^3 - 3/2*X^2 +  1/2*X
      B(4,X) = X^4 - 2*X^3   +      X^2 - 1/30
      B(5,X) = X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
      B(6,X) = X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
      B(7,X) = X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
      B(8,X) = X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the Bernoulli polynomial to
    be evaluated.  N must be 0 or greater.
    Input, double X, the value at which the polynomial is to
    be evaluated.
    Output, double BERNOULLI_POLY2, the value of B(N,X).
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
    ityp fact = 1.00;
    ityp value = bernoulli_number3 ( 0 );

    for (dim_typ i = 1; i <= n; ++i )
    {
        fact *= ( ityp ) ( n + 1 - i ) / ( ityp ) i;
        value = value * x + fact * bernoulli_number3 ( i );
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cardan_poly ( void * data)
/******************************************************************************/
/*
  Purpose:
    CARDAN_POLY evaluates the Cardan polynomials.
  First terms:
     N  C(N,S,X)
     0  2
     1  X
     2  X^2  -  2 S
     3  X^3  -  3 S X
     4  X^4  -  4 S X^2 +  2 S^2
     5  X^5  -  5 S X^3 +  5 S^2 X
     6  X^6  -  6 S X^4 +  9 S^2 X^2 -  2 S^3
     7  X^7  -  7 S X^5 + 14 S^2 X^3 -  7 S^3 X
     8  X^8  -  8 S X^6 + 20 S^2 X^4 - 16 S^3 X^2 +  2 S^4
     9  X^9  -  9 S X^7 + 27 S^2 X^5 - 30 S^3 X^3 +  9 S^4 X
    10  X^10 - 10 S X^8 + 35 S^2 X^6 - 50 S^3 X^4 + 25 S^4 X^2 -  2 S^5
    11  X^11 - 11 S X^9 + 44 S^2 X^7 - 77 S^3 X^5 + 55 S^4 X^3 - 11 S^5 X
  Recursion:
    Writing the N-th polynomial in terms of its coefficients:
      C(N,S,X) = sum ( 0 <= I <= N ) D(N,I) * S^(N-I)/2 * X^I
    then
    D(0,0) = 1
    D(1,1) = 1
    D(1,0) = 0
    D(N,N) = 1
    D(N,K) = D(N-1,K-1) - D(N-2,K)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Thomas Osler,
    Cardan Polynomials and the Reduction of Radicals,
    Mathematics Magazine,
    Volume 74, Number 1, February 2001, pages 26-32.
  Parameters:
    Input, int N, the highest polynomial to compute.
    Input, double X, the point at which the polynomials are to be computed.
    Input, double S, the value of the parameter, which must be positive.
    Output, double CARDAN_POLY[N+1], the values of the Cardan polynomials at X.
*/
{
	const _2itdt * const s_data = data;
	
	const register ityp x = s_data->a0;
	const register ityp s = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    ityp fact;
    ityp s2;
    ityp *v;
    ityp x2[1];

    s2 = sqrt ( s );
    x2[0] = 0.50 * x / s2;
    v = cheby_t_poly ( 1, n, x2 );
    fact = 1.00;

    for (dim_typ i = 0; i <= n; ++i )
    {
        v[i] = 2.00 * fact * v[i];
        fact *= s2;
    }

    return v;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cardan_poly_coef ( void * data)
/******************************************************************************/
/*
  Purpose:
    CARDAN_POLY_COEF computes the coefficients of the N-th Cardan polynomial.
  First terms:
    2
    0       1
   -2 S     0       1
    0      -3 S     0       1
    2 S^2   0      -4 S     0       1
    0       5 S^2   0      -5 S     0       1
   -2 S^3   0       9 S^2   0      -6 S     0       1
    0       7 S^3   0      14 S^2   0      -7 S     0       1
    2 S^4   0     -16 S^3   0      20 S^2   0      -8 S     0        1
    0       9 S^4   0     -30 S^3   0      27 S^2   0      -9 S      0     1
   -2 S^5   0      25 S^4   0     -50 S^3   0      35 S^2   0      -10 S   0   1
    0     -11 S^5   0      55 S^4   0     -77 S^3   0     +44 S^2    0   -11 S 0 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Thomas Osler,
    Cardan Polynomials and the Reduction of Radicals,
    Mathematics Magazine,
    Volume 74, Number 1, February 2001, pages 26-32.
  Parameters:
    Input, int N, the order of the polynomial
    Input, double S, the value of the parameter, which must be positive.
    Output, double C[N+1], the coefficients.  C(0) is the constant term,
    and C(N) is the coefficient of X^N.
*/
{
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * c = s_data->a1;
	const register ityp s = s_data->a2;
	
	
    ityp *cm1;
    ityp *cm2;
    dim_typ i, j;

    c[0] = 2.00;
    for ( i = 1; i <= n; ++i )
        c[i] = 0.00;

    if ( n == 0 )
        return NULL;

    cm1 = ( ityp * ) malloc ( ( n + 1 ) * sizeof ( ityp ) );
    cm2 = ( ityp * ) malloc ( ( n + 1 ) * sizeof ( ityp ) );

    for ( i = 0; i <= n; ++i )
        cm1[i] = c[i];

    c[0] = 0.00;
    c[1] = 1.00;
    for ( i = 2; i <= n; ++i)
        c[i] = 0.0;

    for ( i = 2; i <= n; ++i )
    {

        for ( j = 0; j <= i-2; ++j )
            cm2[j] = cm1[j];

        for ( j = 0; j <= i-1; ++j)
            cm1[j] = c[j];

        c[0] = 0.00;
        for ( j = 1; j <= i; ++j )
            c[j] = cm1[j-1];

        for ( j = 0; j <= i-2; ++j )
            c[j] -= s * cm2[j];

    }

    free ( cm1 );
    free ( cm2 );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cardinal_cos ( void * data)
/******************************************************************************/
/*
  Purpose:
    CARDINAL_COS evaluates the J-th cardinal cosine basis function.
  Discussion:
    The base points are T(I) = M_PI * I / ( M + 1 ), 0 <= I <= M + 1.
    Basis function J is 1 at T(J), and 0 at T(I) for I /= J
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2014
  Author:
    John Burkardt
  Reference:
    John Boyd,
    Exponentially convergent Fourier-Chebyshev quadrature schemes on
    bounded and infinite intervals,
    Journal of Scientific Computing,
    Volume 2, Number 2, 1987, pages 99-109.
  Parameters:
    Input, int J, the index of the basis function.
    0 <= J <= M + 1.
    Input, int M, indicates the size of the basis set.
    Input, int N, the number of sample points.
    Input, double T[N], one or more points in [0,M_PI] where the
    basis function is to be evaluated.
    Output, double CARDINAL_COS[N], the value of the function at T.
*/
{
	const _3dtpit * const s_data = data;
	const register dim_typ j = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * t = s_data->a3;
	 
    ityp *c;
    ityp cj;
    ityp tj;

    c = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    cj = 1.00 + ( j % ( m + 1 ) ) == 0;
    tj = M_PI * ( ityp ) ( j ) / ( ityp ) ( m + 1 );

    for (dim_typ i = 0; i < n; ++i )
        c[i] = fabs ( t[i] - tj ) <= 2.220446049250313E-016 ? 1.00:r8_mop ( j + 1 )* sin ( t[i] )* sin ( ( ityp ) ( m + 1 ) * t[i] )/ cj/ ( ityp ) ( m + 1 )/ ( cos ( t[i] ) - cos ( tj ) );

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _cardinal_sin ( void * data)
/******************************************************************************/
/*
  Purpose:
    CARDINAL_SIN evaluates the J-th cardinal sine basis function.
  Discussion:
    The base points are T(I) = M_PI * I / ( M + 1 ), 0 <= I <= M + 1.
    Basis function J is 1 at T(J), and 0 at T(I) for I /= J
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2014
  Author:
    John Burkardt
  Reference:
    John Boyd,
    Exponentially convergent Fourier-Chebyshev quadrature schemes on
    bounded and infinite intervals,
    Journal of Scientific Computing,
    Volume 2, Number 2, 1987, pages 99-109.
  Parameters:
    Input, int J, the index of the basis function.
    0 <= J <= M + 1.
    Input, int M, indicates the size of the basis set.
    Input, int N, the number of sample points.
    Input, double T[N], one or more points in [0,M_PI] where the
    basis function is to be evaluated.
    Output, double CARDINAL_SIN[N], the value of the function at T.
*/
{
	const _3dtpit * const s_data = data;
	const register dim_typ j = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	ityp * t = s_data->a3;
	
    ityp *s;
    ityp tj;

    s = ( ityp * ) malloc ( n * sizeof ( s ) );
    tj = M_PI * ( ityp ) ( j ) / ( ityp ) ( m + 1 );

    for (dim_typ i = 0; i < n; ++i)
        s[i] = fabs ( t[i] - tj ) <= M_PI ? 1.00 : r8_mop ( j + 1 )* sin ( tj )* sin ( ( ityp ) ( m + 1 ) * t[i] )/ ( ityp ) ( m + 1 )/ ( cos ( t[i] ) - cos ( tj ) );

    return s;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _catalan ( void * data)
/******************************************************************************/
/*
  Purpose:
    CATALAN computes the Catalan numbers, from C(0) to C(N).
  Discussion:
    The Catalan number C(N) counts:
    1) the number of binary trees on N vertices;
    2) the number of ordered trees on N+1 vertices;
    3) the number of full binary trees on 2N+1 vertices;
    4) the number of well formed sequences of 2N parentheses;
    5) the number of ways 2N ballots can be counted, in order,
       with N positive and N negative, so that the running sum
       is never negative;
    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
    7) the number of monotone functions from [1..N} to [1..N} which
       satisfy f(i) <= i for all i;
    8) the number of ways to triangulate a polygon with N+2 vertices.
    The formula is:
      C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
           = 1 / (N+1) * COMB ( 2N, N )
           = 1 / (2N+1) * COMB ( 2N+1, N+1).
  First values:
     C(0)     1
     C(1)     1
     C(2)     2
     C(3)     5
     C(4)    14
     C(5)    42
     C(6)   132
     C(7)   429
     C(8)  1430
     C(9)  4862
    C(10) 16796
  Recursion:
    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
  Example:
    N = 3
 ()()()
 ()(())
 (()())
 (())()
 ((()))
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 May 2003
  Author:
    John Burkardt
  Reference:
    Dennis Stanton and Dennis White,
    Constructive Combinatorics,
    Springer Verlag, New York, 1986.
  Parameters:
    Input, int N, the number of Catalan numbers desired.
    Output, int C[N+1], the Catalan numbers from C(0) to C(N).
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * c = s_data->a1;
	
    int i;
    c[0] = 1;
    /*
    The extra parentheses ensure that the integer division is
    done AFTER the integer multiplication.
    */
    for ( i = 1; i <= n; ++i )
        c[i] = ( c[i-1] * (( (i<<1) - 1 ) << 1) ) / ( i + 1 );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _catalan_row_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    CATALAN_ROW_NEXT computes row N of Catalan's triangle.
  Example:
    I\J 0   1   2   3   4   5   6
    0   1
    1   1   1
    2   1   2   2
    3   1   3   5   5
    4   1   4   9  14  14
    5   1   5  14  28  42  42
    6   1   6  20  48  90 132 132
  Recursion:
    C(0,0) = 1
    C(I,0) = 1
    C(I,J) = 0 for I < J
    C(I,J) = C(I,J-1) + C(I-1,J)
    C(I,I) is the I-th Catalan number.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int NEXT, indicates whether this is a call for
    the 'next' row of the triangle.
    NEXT = 0, this is a startup call.  Row N is desired, but
    presumably this is a first call, or row N-1 was not computed
    on the previous call.
    NEXT = 1, this is not the first call, and row N-1 was computed
    on the previous call.  In this case, much work can be saved
    by using the information from the previous values of IROW
    to build the next values.
    Input, int N, the index of the row of the triangle desired.
    Input/output, int IROW(0:N), the row of coefficients.
    If NEXT = FALSE, then IROW is not required to be set on input.
    If NEXT = TRUE, then IROW must be set on input to the value of
    row N-1.
*/
{
	const bdtpi * const s_data = data;
	const bool next = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * irow = s_data->a2;
	
    dim_typ i, j;

    if ( !next )
    {
        irow[0] = 1;
        for ( i = 1; i <= n; ++i)
            irow[i] = 0;

        for ( i = 1; i <= n; ++i )
        {
            irow[0] = 1;

            for ( j = 1; j <= i-1; ++j )
                irow[j] += irow[j-1];
            irow[i] = irow[i-1];

        }
    }
    else
    {
        irow[0] = 1;
        for ( j = 1; j <= n-1; ++j)
            irow[j] += irow[j-1];

        if ( 1 <= n )
            irow[n] = irow[n-1];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _charlier ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHARLIER evaluates Charlier polynomials at a point.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 March 2009
  Author:
    John Burkardt
  Reference:
    J Simoes Pereira,
    Algorithm 234: Poisson-Charliers Polynomials,
    Communications of the ACM,
    Volume 7, Number 7, page 420, July 1964.
    Walter Gautschi,
    Orthogonal Polynomials: Computation and Approximation,
    Oxford, 2004,
    ISBN: 0-19-850672-4,
    LC: QA404.5 G3555.
    Gabor Szego,
    Orthogonal Polynomials,
    American Mathematical Society, 1975,
    ISBN: 0821810235,
    LC: QA3.A5.v23.
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 2002,
    Second edition,
    ISBN: 1584883472,
    LC: QA5.W45.
  Parameters:
    Input, int N, the maximum order of the polynomial.
    N must be at least 0.
    Input, double A, the parameter.  A must not be 0.
    Input, double X, the evaluation point.
    Output, double VALUE[0:N], the value of the polynomials at X.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp x = s_data->a2;
	ityp * value = s_data->a3;
	
    if ( a == 0.00 )
        return NULL;

    value[0] = 1.00;

    if ( n == 0 )
        return NULL;

    value[1] = - x / a;

    if ( n == 1 )
        return NULL;

    for (dim_typ i = 1; i < n; ++i )
        value[i+1] = ( ( i + a - x ) * value[i] - i * value[i-1] ) / a;
        
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _chebyshev_discrete ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials at a point.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 March 2009
  Author:
    John Burkardt
  Reference:
    Walter Gautschi,
    Orthogonal Polynomials: Computation and Approximation,
    Oxford, 2004,
    ISBN: 0-19-850672-4,
    LC: QA404.5 G3555.
  Parameters:
    Input, int N, the highest order of the polynomials to
    be evaluated.  0 <= N <= M.
    Input, int M, the maximum order of the polynomials.
    0 <= M.
    Input, double X, the evaluation point.
    Output, double V[N+1], the value of the polynomials at X.
*/
{
	const _2dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * v = s_data->a2;
	const register ityp x = s_data->a3;
	
    if ( m < n )
        return NULL;

    v[0] = 1.0;

    if ( n == 0 )
        return NULL;

    v[1] = 2.00 * x + ( ityp ) ( 1 - m );

    if ( n == 1 )
        return NULL;

    for (dim_typ i = 1; i < n; ++i)
        v[i+1] = (( ityp ) ( (i<<1) + 1 )* ( 2.00 * x + ( ityp ) ( 1 - m ) ) * v[i]- ( ityp ) ( i * ( m + i ) * ( m - i ) ) * v[i-1]) / ( ityp ) ( i + 1 );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _collatz_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    COLLATZ_COUNT counts the number of terms in a Collatz sequence.
  Discussion:
    The rules for generation of the Collatz sequence are recursive.
    If T is the current entry of the sequence, (T is
    assumed to be a positive integer), then the next
    entry, U is determined as follows:
      if T is 1 (or less)
        terminate the sequence;
      else if T is even
        U = T/2.
      else (if T is odd and not 1)
        U = 3*T+1;
     N  Sequence                                                Length
     1                                                               1
     2   1                                                           2
     3  10,  5, 16,  8,  4,  2,  1                                   8
     4   2   1                                                       3
     5  16,  8,  4,  2,  1                                           6
     6   3, 10,  5, 16,  8,  4,  2,  1                               9
     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
     8   4,  2,  1                                                   4
     9  28, 14,  7, ...                                             20
    10   5, 16,  8,  4,  2,  1                                       7
    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 March 2006
  Author:
    John Burkardt
  Reference:
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 2002,
    Second edition,
    ISBN: 1584883472,
    LC: QA5.W45.
  Parameters:
    Input, int N, the first element of the sequence.
    Output, int COLLATZ_COUNT, the number of elements in
    the Collatz sequence that begins with N.
*/
{
	static dim_typ result = USHRT_MAX;
	
	register dim_typ n = *(dim_typ *) data;
	
    dim_typ count = 1;

    for ( ; ; )
    {
        if ( n <= 1 )
            break;
        else if ( ( n % 2 ) == 0 )
            n >>= 1;
        else
            n = 3 * n + 1;
        ++ count;
    }

	result = count;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _collatz_count_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    COLLATZ_COUNT_MAX seeks the maximum Collatz count for 1 through N.
  Discussion:
    For each integer I, we compute a sequence of values that
    terminate when we reach 1.  The number of steps required to
    reach 1 is the "rank" of I, and we are searching the numbers
    from 1 to N for the number with maximum rank.
    For a given I, the sequence is produced by:
    1) J = 1, X(J) = I;
    2) If X(J) = 1, stop.
    3) J = J + 1;
       if X(J-1) was even, X(J) = X(J-1)/2;
       else                X(J) = 3 * X(J-1) + 1;
    4) Go to 3
  Example:
            N      I_MAX J_MAX
           10          9    20
          100         97   119
        1,000        871   179
       10,000      6,171   262
      100,000     77,031   351
    1,000,000    837,799   525
   10,000,000  8,400,511   686
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 April 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the maximum integer to check.
    Output, int *I_MAX, *J_MAX, an integer I with the maximum rank,
    and the value of the maximum rank.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * i_max = s_data->a1;
	int * j_max = s_data->a2;
	
    dim_typ i, j, x;
    *i_max = *j_max = -1;

    for ( i = 1; i <= n; ++i )
    {
        j = 1;
        x = i;

        while ( x != 1 )
        {
            ++ j;
            x = ( x % 2 ) == 0 ? x/2 : 3 * x + 1;
        }

        if ( *j_max < j )
        {
            *i_max = i;
            *j_max = j;
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _comb_row_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMB_ROW_NEXT computes the next row of Pascal's triangle.
  Discussion:
    Row N contains the N+1 combinatorial coefficients
      C(N,0), C(N,1), C(N,2), ... C(N,N)
  Discussion:
    The sum of the elements of row N is equal to 2**N.
    The formula is:
      C(N,K) = N! / ( K! * (N-K)! )
  First terms:
     N K:0  1   2   3   4   5   6   7  8  9 10
     0   1
     1   1  1
     2   1  2   1
     3   1  3   3   1
     4   1  4   6   4   1
     5   1  5  10  10   5   1
     6   1  6  15  20  15   6   1
     7   1  7  21  35  35  21   7   1
     8   1  8  28  56  70  56  28   8  1
     9   1  9  36  84 126 126  84  36  9  1
    10   1 10  45 120 210 252 210 120 45 10  1
  Recursion:
    C(N,K) = C(N-1,K-1)+C(N-1,K)
  Special values:
    C(N,0) = C(N,N) = 1
    C(N,1) = C(N,N-1) = N
    C(N,N-2) = sum ( 1 <= I <= N ) N
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 December 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, indicates the desired row.
    Input/output, int ROW[N+1].  On input, row N-1 is
    contained in entries 0 through N-1.  On output, row N is contained
    in entries 0 through N.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * row = s_data->a1;
	
    row[n] = 1;
    for (dim_typ i = n - 1; 1 <= i; --i )
        row[i] += row[i-1];
    row[0] = 1;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _commul ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMMUL computes a multinomial combinatorial coefficient.
  Discussion:
    The multinomial coefficient is a generalization of the binomial
    coefficient.  It may be interpreted as the number of combinations of
    N objects, where IARRAY(1) objects are indistinguishable of type 1,
    ... and IARRAY(K) are indistinguishable of type NFACT.
    The formula is:
      COMMUL = N! / ( IARRAY(1)! IARRAY(2)! ... IARRAY(NFACT)! )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, determines the numerator.
    Input, int NFACT, the number of factors in the numerator.
    Input, int IARRAY(NFACT).
    IARRAY contains the NFACT values used in the denominator.
    Note that the sum of these entries should be N,
    and that all entries should be nonnegative.
    Output, int COMMUL, the value of the multinomial coefficient.
*/
{
	static int result = INT_MAX;
	
	const dtpii * const s_data = data;
	
	const register dim_typ nfact = s_data->a0;
	int * iarray = s_data->a1;
	const register int n = s_data->a2;
	
    ityp arg;
    ityp fack;
    ityp facn;
    dim_typ i;
    int isum;

    for ( i = 0; i < nfact; ++i )
        if ( iarray[i] < 0 )
        {
        	result = INT_MAX;
            return &result;
        }

    isum = 0;
    for ( i = 0; i < nfact; ++i )
        isum += iarray[i];

    if ( isum != n )
    {
    	result = INT_MAX;
        return &result;
    }

    arg = ( ityp ) ( n + 1 );
    facn = lgamma ( arg );

    for ( i = 0; i < nfact; ++i )
    {
        arg = ( ityp ) ( iarray[i] + 1 );
        fack = lgamma ( arg );
        facn = facn - fack;
    }

	result = r8_nint ( exp ( facn ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _complete_symmetric_poly ( void * data)
/******************************************************************************/
/*
  Purpose:
    COMPLETE_SYMMETRIC_POLY evaluates a complete symmetric polynomial.
  Discussion:
    N\R  0   1         2               3
      +--------------------------------------------------------
    0 |  1   0         0               0
    1 |  1   X1        X1^2            X1^3
    2 |  1   X1+X2     X1^2+X1X2+X2^2  X1^3+X1^2X2+X1X2^2+X2^3
    3 |  1   X1+X2+X3  ...
    If X = ( 1, 2, 3, 4, 5, ... ) then
    N\R  0     1     2     3     4 ...
      +--------------------------------------------------------
    0 |  1     0     0     0     0
    1 |  1     1     1     1     1
    2 |  1     3     7    15    31
    3 |  1     6    25    90   301
    4 |  1    10    65   350  1701
    5 |  1    15   140  1050  6951
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of variables.
    0 <= N.
    Input, int R, the degree of the polynomial.
    0 <= R.
    Input, double X[N], the value of the variables.
    Output, double COMPLETE_SYMMETRIC_POLY, the value of TAU(N,R)(X).
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ r = s_data->a1;
	ityp * x = s_data->a2;
	
    dim_typ i;
    dim_typ nn;
    dim_typ rr;
    ityp *tau;
    ityp value;

    if ( n < 0 || r < 0 )
    {
    	result = MAX_VAL;
        return &result;
	}
	
    tau = ( ityp * ) malloc ( ( 1 + MAX ( n, r ) ) * sizeof ( ityp ) );

    for ( i = 0; i <= MAX ( n, r ); ++i )
        tau[i] = 0.00;

    tau[0] = 1.00;
    for ( nn = 1; nn <= n; ++nn )
        for ( rr = 1; rr <= r; ++rr )
            tau[rr] += x[nn-1] * tau[rr-1];

    value = tau[r];
    free ( tau );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cos_power_int ( void * data)
/******************************************************************************/
/*
  Purpose:
    COS_POWER_INT evaluates the cosine power integral.
  Discussion:
    The function is defined by
      COS_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( cos ( t ))^n dt
    The algorithm uses the following fact:
      Integral cos^n ( t ) = -(1/n) * (
        cos^(n-1)(t) * sin(t) + ( n-1 ) * Integral cos^(n-2) ( t ) dt )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 March 2012
  Author:
    John Burkardt
  Parameters
    Input, double A, B, the limits of integration.
    Input, integer N, the power of the sine function.
    Output, double COS_POWER_INT, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const _2itdt * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    ityp ca;
    ityp cb;
    dim_typ m;
    dim_typ mlo;
    ityp sa;
    ityp sb;
    ityp value;

    if ( n < 0 )
    {
    	result = MAX_VAL;
        return &result;
    }

    sa = sin ( a );
    sb = sin ( b );
    ca = cos ( a );
    cb = cos ( b );

    if ( ( n % 2 ) == 0 )
    {
        value = b - a;
        mlo = 2;
    }
    else
    {
        value = sb - sa;
        mlo = 3;
    }

    for ( m = mlo; m <= n; m += 2 )
        value = ( ( ityp ) ( m - 1 ) * value- pow ( ca, (m-1) ) * sa + pow ( cb, (m-1) ) * sb )/ ( ityp ) ( m );

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _euler_number ( void * data)
/******************************************************************************/
/*
  Purpose:
    EULER_NUMBER computes the Euler numbers.
  Discussion:
    The Euler numbers can be evaluated in Mathematica with the call
      EulerE[n]
    These numbers rapidly get too big to store in an ordinary integer!
    The terms of odd index are 0.
    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
  First terms:
    E0  = 1
    E1  = 0
    E2  = -1
    E3  = 0
    E4  = 5
    E5  = 0
    E6  = -61
    E7  = 0
    E8  = 1385
    E9  = 0
    E10 = -50521
    E11 = 0
    E12 = 2702765
    E13 = 0
    E14 = -199360981
    E15 = 0
    E16 = 19391512145
    E17 = 0
    E18 = -2404879675441
    E19 = 0
    E20 = 370371188237525
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 February 2015
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input, int N, the index of the last Euler number to compute.
    Output, int E[N+1], the Euler numbers from index 0 to N.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * e = s_data->a1;
	
    dim_typ i, j;

    e[0] = 1;

    if ( n == 0 )
        return NULL;

    e[1] = 0;

    if ( n == 1 )
        return NULL;

    e[2] = -1;

    for ( i = 3; i <= n; ++i )
    {
        e[i] = 0;
        if ( ( i % 2 ) == 0 )
            for ( j = 2; j <= i; j += 2 )
            e[i] -= i4_choose ( i, j ) * e[i-j];
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _euler_number2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    EULER_NUMBER2 computes the Euler numbers.
  Discussion:
    The Euler numbers can be evaluated in Mathematica with the call
      EulerE[n]
  First terms:
    E0  = 1
    E1  = 0
    E2  = -1
    E3  = 0
    E4  = 5
    E5  = 0
    E6  = -61
    E7  = 0
    E8  = 1385
    E9  = 0
    E10 = -50521
    E11 = 0
    E12 = 2702765
    E13 = 0
    E14 = -199360981
    E15 = 0
    E16 = 19391512145
    E17 = 0
    E18 = -2404879675441
    E19 = 0
    E20 = 370371188237525
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input, int N, the index of the Euler number to compute.
    Output, double EULER_NUMBER2, the value of E(N).
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
    ityp e[7] =
    {
        1.0, -1.0, 5.0, -61.0, 1385.0,
        -50521.0, 2702765.0
    };
    dim_typ i;
    dim_typ itmax = 1000;
    ityp sum1;
    ityp term;
    ityp value;

    if ( n == 0 )
    {
    	result = e[0];
        return &result;
    }

    if ( ( n % 2 ) == 1 )
    {
    	result = 0.00;
        return &result;
    }

    if ( n <= 12 )
    {
    	result = e[n>>1];
        return &result;
    }

    sum1 = 0.00;

    for ( i = 1; i <= itmax; ++i )
    {
        term = 1.00 / pow ( ( ityp ) ( (i<<1) - 1 ), n + 1 );
        sum1 += ( i % 2 ) == 1 ? term:-term;

        if ( fabs ( term ) < 1.0E-10 || fabs ( term ) < 1.0E-08 * fabs ( sum1 ) )
            break;
    }

    value = pow ( 2.0, n + 2 ) * sum1 * r8_factorial ( n )/ pow ( M_PI, n + 1 );
    
    result = ( n % 4 ) != 0 ? -value:value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _euler_poly ( void * data)
/******************************************************************************/
/*
  Purpose:
    EULER_POLY evaluates the N-th Euler polynomial at X.
  First values:
    E(0,X) = 1
    E(1,X) = X - 1/2
    E(2,X) = X^2 - X
    E(3,X) = X^3 - 3/2 X^2 + 1/4
    E(4,X) = X^4 - 2*X^3 + X
    E(5,X) = X^5 - 5/2 X^4 + 5/2 X^2 - 1/2
    E(6,X) = X^6 - 3 X^5 + 5 X^3 - 3 X
    E(7,X) = X^7 - 7/2 X^6 + 35/4 X^4 - 21/2 X^2 + 17/8
    E(8,X) = X^8 - 4 X^7 + 14 X^5 - 28 X^3 + 17 X
  Special values:
    E'(N,X) = N * E(N-1,X)
    E(N,1/2) = E(N) / 2**N, where E(N) is the N-th Euler number.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 February 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the Euler polynomial to
    be evaluated.  N must be 0 or greater.
    Input, double X, the value at which the polynomial is to
    be evaluated.
    Output, double EULER_POLY, the value of E(N,X).
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp x = s_data->a1;
	
	result = 2.00 * ( bernoulli_poly2 ( n+1, x ) - bernoulli_poly2 ( n+1, 0.5 * x ) * pow ( ( ityp ) 2, n+1 ) )/ ( ityp ) ( n + 1 );
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _eulerian ( void * data)
/******************************************************************************/
/*
  Purpose:
    EULERIAN computes the Eulerian number E(N,K).
  Definition:
    A run in a permutation is a sequence of consecutive ascending values.
    E(N,K) is the number of permutations of N objects which contain
    exactly K runs.
  Examples:
     N = 7
     1     0     0     0     0     0     0
     1     1     0     0     0     0     0
     1     4     1     0     0     0     0
     1    11    11     1     0     0     0
     1    26    66    26     1     0     0
     1    57   302   302    57     1     0
     1   120  1191  2416  1191   120     1
  Recursion:
    E(N,K) = K * E(N-1,K) + (N-K+1) * E(N-1,K-1).
  Properties:
    E(N,1) = E(N,N) = 1.
    E(N,K) = 0 if K <= 0 or N < K.
    sum ( 1 <= K <= N ) E(N,K) = N!.
    X^N = sum ( 0 <= K <= N ) COMB(X+K-1, N ) E(N,K)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Dennis Stanton and Dennis White,
    Constructive Combinatorics,
    Springer Verlag, 1986
  Parameters:
    Input, int N, the number of rows desired.
    Output, int E[N*N], the first N rows of Eulerian numbers.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * e = s_data->a1;
	
    dim_typ i, j;

    if ( n < 1 )
        return NULL;
    /*
    Construct rows 1, 2, ..., N of the Eulerian triangle.
    */
    e[1-1+(1-1)*n] = 1;
    for ( j = 2; j <= n; ++j)
        e[1-1+(j-1)*n] = 0;

    for ( i = 2; i <= n; ++i )
    {
        e[i-1+(1-1)*n] = 1;
        for ( j = 2; j <= n; ++j )
            e[i-1+(j-1)*n] = j * e[i-2+(j-1)*n] + ( i - j + 1 ) * e[i-2+(j-2)*n];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _f_hofstadter ( void * data)
/******************************************************************************/
/*
  Purpose:
    F_HOFSTADTER computes the Hofstadter F sequence.
  Discussion:
    F(N) = 0                if N = 0
         = N - F ( N - 1 ), otherwise.
    F(N) is defined for all nonnegative integers, and turns out
    to be equal to int ( ( N + 1 ) / 2 ).
  Table:
     N  F(N)
    --  ----
     0     0
     1     1
     2     1
     3     2
     4     2
     5     3
     6     3
     7     4
     8     4
     9     5
    10     5
    11     6
    12     6
    13     7
    14     7
    15     8
    16     8
    17     9
    18     9
    19    10
    20    10
    21    11
    22    11
    23    12
    24    12
    25    13
    26    13
    27    14
    28    14
    29    15
    30    15
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Douglas Hofstadter,
    Goedel, Escher, Bach,
    Basic Books, 1979.
  Parameters:
    Input, int N, the argument of the function.
    Output, int F_HOFSTADTER, the value of the function.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = 0 + (n>0)*( n - f_hofstadter ( n-1 ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _g_hofstadter ( void * data) 
/******************************************************************************/
/*
  Purpose:
    G_HOFSTADTER computes the Hofstadter G sequence.
  Discussion:
    G(N) = 0                      if N = 0
         = N - G ( G ( N - 1 ) ), otherwise.
    G(N) is defined for all nonnegative integers.
    The value of G(N) turns out to be related to the Zeckendorf
    representation of N as a sum of non-consecutive Fibonacci numbers.
    To compute G(N), determine the Zeckendorf representation:
      N = sum ( 1 <= I <= M ) F(I)
    and reduce the index of each Fibonacci number by 1:
      G(N) = sum ( 1 <= I <= M ) F(I-1)
    However, this is NOT how the computation is done in this routine.
    Instead, a straightforward recursive function call is defined
    to correspond to the definition of the mathematical function.
  Table:
     N  G(N)  Zeckendorf   Decremented
    --  ----  ----------   -----------
     1   1    1            1
     2   1    2            1
     3   2    3            2
     4   3    3 + 1        2 + 1
     5   3    5            3
     6   4    5 + 1        3 + 1
     7   4    5 + 2        3 + 1
     8   5    8            5
     9   6    8 + 1        5 + 1
    10   6    8 + 2        5 + 1
    11   7    8 + 3        5 + 2
    12   8    8 + 3 + 1    5 + 2 + 1
    13   8    13           8
    14   9    13 + 1       8 + 1
    15   9    13 + 2       8 + 1
    16  10    13 + 3       8 + 2
    17  11    13 + 3 + 1   8 + 2 + 1
    18  11    13 + 5       8 + 3
    19  12    13 + 5 + 1   8 + 3 + 1
    20  12    13 + 5 + 2   8 + 3 + 1
    21  13    21           13
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Douglas Hofstadter,
    Goedel, Escher, Bach,
    Basic Books, 1979.
  Parameters:
    Input, int N, the argument of the function.
    Output, int G_HOFSTADTER, the value of the function.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = 0 + (n>0)*( n - g_hofstadter ( g_hofstadter ( n-1 ) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gegenbauer_poly ( void * data)
/******************************************************************************/
/*
  Purpose:
    GEGENBAUER_POLY computes the Gegenbauer polynomials C(0:N,ALPHA,X).
  Discussion:
    The Gegenbauer polynomial can be evaluated in Mathematica with
    the command
      GegenbauerC[n,m,x]
  Differential equation:
 (1-X*X) Y'' - (2 ALPHA + 1) X Y' + N (N + 2 ALPHA) Y = 0
  Recursion:
    C(0,ALPHA,X) = 1,
    C(1,ALPHA,X) = 2*ALPHA*X
    C(N,ALPHA,X) = ( 2*(N-1+ALPHA)*X * C(N-1,ALPHA,X) - (N-2+2*ALPHA) * C(N-2,ALPHA,X) )/N
  Restrictions:
    ALPHA must be greater than -0.5.
  Special values:
    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
    polynomials of the second kind.
  Norm:
    Integral ( -1 <= X <= 1 )
  ( 1 - X^2 )^( ALPHA - 0.5 ) * C(N,ALPHA,X)^2 dX
    = M_PI * 2^( 1 - 2 * ALPHA ) * Gamma ( N + 2 * ALPHA )
      / ( N! * ( N + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 August 2004
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input, int N, the highest order polynomial to compute.
    Note that polynomials 0 through N will be computed.
    Input, double ALPHA, a parameter which is part of the definition of
    the Gegenbauer polynomials.  It must be greater than -0.5.
    Input, double X, the point at which the polynomials are to be evaluated.
    Output, double CX[N+1], the values of the first N+1 Gegenbauer
    polynomials at the point X.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp alpha = s_data->a1;
	const register ityp x = s_data->a2;
	ityp * cx = s_data->a3;
	
    if ( alpha <= -0.5 )
        return NULL;

    cx[0] = 1.00;

    if ( n == 0 )
        return NULL;

    cx[1] = 2.00 * alpha * x;

    for (dim_typ i = 2; i <= n; ++i )
        cx[i] = ( ( ( ityp ) ( (i<<1) - 2 ) + 2.00 * alpha ) * x * cx[i-1]+ ( ( ityp ) (   - i + 2 ) - 2.00 * alpha )     * cx[i-2] )/ ( ityp )       i;
    
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gud ( void * data)
/******************************************************************************/
/*
  Purpose:
    GUD evaluates the Gudermannian function.
  Definition:
    The Gudermannian function relates the hyperbolic and trigonometric
    functions.  For any argument X, there is a corresponding value
    GAMMA so that
      sinh(x) = tan(gamma).
    The value GAMMA is called the Gudermannian of X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the Gudermannian.
    Output, double GUD, the value of the Gudermannian.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	result = 2.00 * atan ( tanh ( 0.50 * x ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _h_hofstadter ( void * data)
/******************************************************************************/
/*
  Purpose:
    H_HOFSTADTER computes the Hofstadter H sequence.
  Discussion:
    H(N) = 0                          if N = 0
         = N - H ( H ( H ( N - 1 ) ), otherwise.
    H(N) is defined for all nonnegative integers.
  Table:
     N  H(N)
    --  ----
     0     0
     1     1
     2     1
     3     2
     4     3
     5     4
     6     4
     7     5
     8     5
     9     6
    10     7
    11     7
    12     8
    13     9
    14    10
    15    10
    16    11
    17    12
    18    13
    19    13
    20    14
    21    14
    22    15
    23    16
    24    17
    25    17
    26    18
    27    18
    28    19
    29    20
    30    20
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Douglas Hofstadter,
    Goedel, Escher, Bach,
    Basic Books, 1979.
  Parameters:
    Input, int N, the argument of the function.
    Output, int H_HOFSTADTER, the value of the function.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = 0 + (n>0)*( n - h_hofstadter ( h_hofstadter ( h_hofstadter ( n-1 ) ) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _hail ( void * data)
/******************************************************************************/
/*
  Purpose:
    HAIL computes the hail function.
  Discussion:
    Starting with a positive integer N, we divide it by 2 if it is
    even, or triple and add 1 if odd, and repeat this process until
    reaching the value 1.  The number of times the process is carried
    out is the value of the hail function for the given starting value.
    Actually, HAIL is not well defined, since it is not known if
    the above process actually terminates, let alone terminates at 1,
    for every starting value N.
  Example:
     N  Sequence                                                  Hail
     1                                                               0
     2   1                                                           1
     3  10,  5, 16,  8,  4,  2,  1                                   7
     4   2   1                                                       2
     5  16,  8,  4,  2,  1                                           5
     6   3, 10,  5, 16,  8,  4,  2,  1                               8
     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   16
     8   4,  2,  1                                                   3
     9  28, 14,  7, ...                                             19
    10   5, 16,  8,  4,  2,  1                                       6
    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          14
    12   6,  3, 10,  5, 16,  8,  4,  2,  1                           9
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the starting value for the hail sequence.
    Output, int HAIL, the number of steps before the hail sequence
    reached 1.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ k=0;
    dim_typ m=n;

    if ( 0 < n )
        while ( m != 1 )
        {
            ++ k;
            m = ( m % 2 ) == 0 ? m/2 : 3*m+1;
        }

	result = k;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_partition_distinct_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_PARTITION_DISTINCT_COUNT returns any value of Q(N).
  Discussion:
    A partition of an integer N is a representation of the integer
    as the sum of nonzero positive integers.  The order of the summands
    does not matter.  The number of partitions of N is symbolized
    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
    following partitions:
    5 = 5
      = 4 + 1
      = 3 + 2
      = 3 + 1 + 1
      = 2 + 2 + 1
      = 2 + 1 + 1 + 1
      = 1 + 1 + 1 + 1 + 1
    However, if we require that each member of the partition
    be distinct, we are computing something symbolized by Q(N).
    The number 5 has Q(N) = 3, because it has the following partitions
    into distinct parts:
    5 = 5
      = 4 + 1
      = 3 + 2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 February 2003
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input, int N, the integer to be partitioned.
    Output, int I4_PARTITION_DISTINCT_COUNT, the number of partitions
    of the integer into distinct parts.
*/
{
	static int result = INT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    int *c;
    dim_typ i;
    dim_typ k;
    dim_typ k2;
    int k_sign;
    int value;

    c = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    c[0] = 1;

    for ( i = 1; i <= n; ++i)
    {
        c[i] = i4_is_triangular ( i );
        k = 0;
        k_sign = -1;

        for ( ; ; )
        {
            ++ k;
            k_sign = -k_sign;
            k2 = k * ( 3 * k + 1 );

            if ( i < k2 )
                break;

            c[i] += k_sign * c[i-k2];
        }

        k = 0;
        k_sign = -1;

        for ( ; ; )
        {
            ++ k;
            k_sign = -k_sign;
            k2 = k * ( 3 * k - 1 );

            if ( i < k2 )
                break;
            c[i] += k_sign * c[i-k2];
        }

    }

    value = c[n];
    free ( c );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _krawtchouk ( void * data)
/******************************************************************************/
/*
//  Purpose:
//    KRAWTCHOUK evaluates the Krawtchouk polynomials at X.
//  Discussion:
//    The polynomial has a parameter P, which must be striclty between
//    0 and 1, and a parameter M which must be a nonnegative integer.
//    The Krawtchouk polynomial of order N, with parameters P and M,
//    evaluated at X, may be written K(N,P,X,M).
//    The first two terms are:
//
//      K(0,P,X,M) = 1
//      K(1,P,X,M) = X - P * M
//    and the recursion, for fixed P and M is
//                       ( N + 1 ) * K(N+1,P,X,M) =
//  ( X - ( N + P * ( M - 2 * N))) * K(N,  P,X,M)
//       - ( M - N + 1 ) * P * ( 1 - P ) * K(N-1,P,X,M)
//  Licensing:
//    This code is distributed under the GNU LGPL license.
//  Modified:
//    18 March 2009
//  Author:
//    John Burkardt
//  Reference:
//    Walter Gautschi,
//    Orthogonal Polynomials: Computation and Approximation,
//    Oxford, 2004,
//    ISBN: 0-19-850672-4,
//    LC: QA404.5 G3555.
//  Parameters:
//    Input, int N, the highest order polynomial to evaluate.
//    0 <= N.
//    Input, double P, the parameter.  0 < P < 1.
//    Input, double X, the evaluation parameter.
//    Input, int M, the parameter.  0 <= M.
//    Output, double V[N+1], the values of the Krawtchouk polynomials
//    of orders 0 through N at X.
*/
{
	const dt2itdtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp p = s_data->a1;
	const register ityp x = s_data->a2;
	const register dim_typ m = s_data->a3;
	ityp * v = s_data->a4;
	
    if ( p <= 0.00 || 1.0 <= p )
        return NULL;

    v[0] = 1.00;

    if ( 1 <= n )
        v[1] = x - p * ( ityp ) ( m );


    for (dim_typ i = 1; i < n; ++i )
        v[i+1] = (( x - ( ( ityp ) ( i ) + p * ( ityp ) ( m - (i<<1)) ) ) * v[i]- ( ityp ) ( m - i + 1 ) * p * ( 1.00 - p ) * v[i-1]) / ( ityp ) ( i + 1 );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _lerch ( void * data)
/******************************************************************************/
/*
  Purpose:
    LERCH estimates the Lerch transcendent function.
  Definition:
    The Lerch transcendent function is defined as:
      LERCH ( Z, S, A ) = Sum ( 0 <= K < +oo ) Z**K / ( A + K )**S
    excluding any term with ( A + K ) = 0.
    In Mathematica, the function can be evaluated by:
      LerchPhi[z,s,a]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2005
  Author:
    John Burkardt
  Reference:
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 2002,
    Second edition,
    ISBN: 1584883472,
    LC: QA5.W45.
  Thanks:
    Oscar van Vlijmen
  Parameters:
    Input, double Z, int S, double A, the parameters of the function.
    Output, double LERCH, an approximation to the Lerch
    transcendent function.
*/
{
	static ityp result = MAX_VAL;
	
	const _2itdt * const s_data = data;
	
	const register ityp z = s_data->a0;
	const register ityp a = s_data->a1;
	const register dim_typ s = s_data->a2;
	
	
    ityp eps = 1.0E-10;
    dim_typ k;
    ityp total;
    ityp term;
    ityp z_k;

    if ( z <= 0.0 )
    {
    	result = 0.00; 
        return &result;
    }

    total = 0.00;
    k = 0;
    z_k = 1.00;

    for ( ; ; )
    {
        if ( a + ( ityp ) ( k ) != 0.00 )
        {
            term = z_k / pow ( a + ( ityp ) k, s );

            total += term;

            if ( fabs ( term ) <= eps * ( 1.00 + fabs ( total ) ) )
                break;
        }

        ++ k ;
        z_k *= z;

    }

	result = total;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _meixner ( void * data)
/******************************************************************************/
/*
  Purpose:
    MEIXNER evaluates Meixner polynomials at a point.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 March 2009
  Author:
    John Burkardt
  Reference:
    Walter Gautschi,
    Orthogonal Polynomials: Computation and Approximation,
    Oxford, 2004,
    ISBN: 0-19-850672-4,
    LC: QA404.5 G3555.
  Parameters:
    Input, int N, the maximum order of the polynomial.
    N must be at least 0.
    Input, double BETA, the Beta parameter.  0 < BETA.
    Input, double C, the C parameter.  0 < C < 1.
    Input, double X, the evaluation point.
    Output, double V[N+1], the value of the polynomials at X.
*/
{
	const dt3itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp beta = s_data->a1;
	ityp c = s_data->a2;
	ityp x = s_data->a3;
	ityp * v = s_data->a4;
	
    if (  beta <= 0.00 || c <= 0.00 || 1.00 <= c )
        return NULL;


    v[0] = 1.00;

    if ( n == 0 )
        return NULL;

    v[1] = ( c - 1.00 ) * x / beta / c + 1.00;

    if ( n == 1 )
        return NULL;

    for (dim_typ i = 1; i < n; ++i )
        v[i+1] = (( ( c - 1.00 ) * x + ( 1.00 + c )* ( ityp ) ( i ) + beta * c ) * v[i]- ( ityp ) ( i ) * v[i-1]) / ( ( ityp ) ( i ) + beta );


    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _mertens ( void * data)
/******************************************************************************/
/*
  Purpose:
    MERTENS evaluates the Mertens function.
  Discussion:
    The Mertens function M(N) is the sum from 1 to N of the Moebius
    function MU.  That is,
    M(N) = sum ( 1 <= I <= N ) MU(I)
        N   M(N)
        --  ----
         1     1
         2     0
         3    -1
         4    -1
         5    -2
         6    -1
         7    -2
         8    -2
         9    -2
        10    -1
        11    -2
        12    -2
       100     1
      1000     2
     10000   -23
    100000   -48
    The determinant of the Redheffer matrix of order N is equal
    to the Mertens function M(N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2007
  Author:
    John Burkardt
  Reference:
    M Deleglise, J Rivat,
    Computing the Summation of the Moebius Function,
    Experimental Mathematics,
    Volume 5, 1996, pages 291-295.
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 2002,
    Second edition,
    ISBN: 1584883472,
    LC: QA5.W45
  Parameters:
    Input, int N, the argument.
    Output, int MERTENS, the value.
*/
{
	static short result = SHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ value = 0;
    for (dim_typ i = 1; i <= n; ++i )
        value += moebius ( i );
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _moebius ( void * data)
/******************************************************************************/
/*
  Purpose
    MOEBIUS returns the value of MU(N), the Moebius function of N.
  Discussion:
    MU(N) is defined as follows:
      MU(N) = 1 if N = 1;
              0 if N is divisible by the square of a prime;
        (-1)**K, if N is the product of K distinct primes.
    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
    if N is a square, cube, etc.
    The Moebius function is related to Euler's totient function:
      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
  First values:
     N  MU(N)
     1    1
     2   -1
     3   -1
     4    0
     5   -1
     6    1
     7   -1
     8    0
     9    0
    10    1
    11   -1
    12    0
    13   -1
    14    1
    15    1
    16    0
    17   -1
    18    0
    19   -1
    20    0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 1999
  Author:
    John Burkardt
  Parameters:
    Input, int N, the value to be analyzed.
    Output, int MOEBIUS, the value of MU(N).
    If N is less than or equal to 0, or there was not enough internal
    space for factoring, MOEBIUS is returned as -1.
*/
{
	static short result = SHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    # define FACTOR_MAX 20

    int factor[FACTOR_MAX];
    dim_typ i;
    int nfactor;
    int nleft;
    int power[FACTOR_MAX];
    int value;

    if ( n == 0 )
    {
    	result = (-1);
        return &result;
    }

    if ( n == 1 )
    {
    	result = 1;
        return &result;
    }
    /*
    Factor N.
    */
    i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

    if ( nleft != 1 )
    {
    	result = SHRT_MAX;
        return &result;
    }

    value = 1;

    for ( i = 0; i < nfactor; ++i )
    {
        value *= -1;
        if ( 1 < power[i] )
        {
        	result = 0;
            return &result;
        }
    }

	result = value;
    return &result;
    # undef FACTOR_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _motzkin ( void * data)
/******************************************************************************/
/*
  Purpose:
    MOTZKIN returns the Motzkin numbers up to order N.
  Discussion:
    The Motzkin number A(N) counts the number of distinct paths
    from (0,0) to (0,N) in which the only steps used are
 (1,1), (1,-1) and (1,0), and the path is never allowed to
    go below the X axis.
  First values:
     N  A(N)
     0    1
     1    1
     2    2
     3    4
     4    9
     5   21
     6   51
     7  127
     8  323
     9  835
    10 2188
  Recursion:
    A(N) = A(N-1) + sum ( 0 <= K <= N-2 ) A(K) * A(N-2-K)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 August 2004
  Author:
    John Burkardt
  Reference:
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 2002,
    Second edition,
    ISBN: 1584883472,
    LC: QA5.W45.
  Parameters:
    Input, int N, the highest order Motzkin number to compute.
    Output, int A[N+1], the Motzkin numbers.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;

    a[0] = 1;

    for ( i = 1; i <= n; ++i )
    {
        a[i] = a[i-1];
        for ( j = 0; j <= i-2; ++j )
            a[i] += a[j] * a[i-2-j];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _omega ( void * data)
/******************************************************************************/
/*
  Purpose:
    OMEGA returns OMEGA(N), the number of distinct prime divisors of N.
  Discussion:
    If N = 1, then
      OMEGA(N) = 1
    else if the prime factorization of N is
      N = P1**E1 * P2**E2 * ... * PM**EM,
    then
      OMEGA(N) = M
  First values:
     N   OMEGA(N)
     1    1
     2    1
     3    1
     4    1
     5    1
     6    2
     7    1
     8    1
     9    1
    10    2
    11    1
    12    2
    13    1
    14    2
    15    2
    16    1
    17    1
    18    2
    19    1
    20    2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 November 2000
  Author:
    John Burkardt
  Parameters:
    Input, int N, the value to be analyzed.  N must be 1 or
    greater.
    Output, int OMEGA, the value of OMEGA(N).  But if N is 0 or
    less, OMEGA is returned as 0, a nonsense value.  If there is
    not enough room for factoring, OMEGA is returned as -1.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    # define FACTOR_MAX 20

    int factor[FACTOR_MAX];
    int nfactor;
    int nleft;
    int power[FACTOR_MAX];

    if ( n == 0)
    {
    	result = 0;
        return &result;
    }

    if ( n == 1 )
    {
    	result = 1;
        return &result;
    }
    /*
    Factor N.
    */
    i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );
    
    result = nleft != 1 ? MAX_VAL:nfactor;
    return &result;

    # undef FACTOR_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pentagon_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    PENTAGON_NUM computes the N-th pentagonal number.
  Discussion:
    The pentagonal number P(N) counts the number of dots in a figure of
    N nested pentagons.  The pentagonal numbers are defined for both
    positive and negative N.
    The formula is:
      P(N) = ( N * ( 3 * N - 1 ) ) / 2
  First values:
    N   P
   -5   40
   -4   26
   -3   15
   -2    7
   -1    2
    0    0
    1    1
    2    5
    3   12
    4   22
    5   35
    6   51
    7   70
    8   92
    9  117
   10  145
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the index of the pentagonal number desired.
    Output, int PENTAGON_NUM, the value of the N-th pentagonal number.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ( n * ( 3 * n - 1 ) ) >> 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _phi ( void * data)
/******************************************************************************/
/*
  Purpose:
    PHI computes the number of relatively prime predecessors of an integer.
  Discussion:
    PHI(N) is the number of integers between 1 and N which are
    relatively prime to N.  I and J are relatively prime if they
    have no common factors.  The function PHI(N) is known as
    "Euler's totient function".
    By convention, 1 and N are relatively prime.
    The formula is:
      PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
      PHI(P^K) = P^(K-1) * ( P - 1 ) if P is prime.
      PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
      N = Sum ( D divides N ) PHI(D).
  First values:
     N  PHI(N)
     1    1
     2    1
     3    2
     4    2
     5    4
     6    2
     7    6
     8    4
     9    6
    10    4
    11   10
    12    4
    13   12
    14    6
    15    8
    16    8
    17   16
    18    6
    19   18
    20    8
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the value to be analyzed.
    Output, int PHI, the value of PHI(N).  If N is less than
    or equal to 0, PHI will be returned as 0.  If there is not enough
    room for full factoring of N, PHI will be returned as -1.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    # define FACTOR_MAX 20

    int factor[FACTOR_MAX];
    dim_typ i;
    int nfactor;
    int nleft;
    int power[FACTOR_MAX];
    int value;

    if ( n == 0)
    {
    	result = 0;
        return &result;
    }

    if ( n == 1 )
	{
		result = 1;
        return &result;
	}
    /*
    Factor N.
    */
    i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

    if ( nleft != 1 )
    {
    	result = USHRT_MAX;
        return &result;
    }

    value = 1;
    for ( i = 0; i < nfactor; ++i)
        value *= ( dim_typ ) pow ( ( ityp ) factor[i], power[i]-1 )* ( factor[i] - 1 );
    
    result = value;
    return &result;
    # undef FACTOR_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _plane_partition_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    PLANE_PARTITION_NUM returns the number of plane partitions of the integer N.
  Discussion:
    A plane partition of a positive integer N is a partition of N in which
    the parts have been arranged in a 2D array that is nonincreasing across
    rows and columns.  There are six plane partitions of 3:
      3   2 1   2   1 1 1   1 1   1
                1           1     1
                                  1
  First Values:
     N PP(N)
     0    1
     1    1
     2    3
     3    6
     4   13
     5   24
     6   48
     7   86
     8  160
     9  282
    10  500
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 April 2014
  Author:
    John Burkardt
  Reference:
    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
    NIST Handbook of Mathematical Functions,
    Cambridge University Press, 2010,
    ISBN: 978-0521140638,
    LC: QA331.N57.
  Parameters:
    Input, int N, the number, which must be at least 0.
    Output, int PLANE_PARTITION_NUM, the number of
    plane partitions of N.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ j;
    dim_typ k;
    dim_typ nn;
    int *pp;
    int s2;
    int value;

    if ( n < 0 )
    {
    	result = USHRT_MAX;
        return &result;
    }

    pp = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
    nn = 0;
    pp[nn] = nn = 1;

    nn = 1;
    if ( nn <= n )
        pp[nn] = 1;

    for ( nn = 2; nn <= n; ++nn)
    {
        pp[nn] = 0;
        for ( j = 1; j <= nn; ++j )
        {
            s2 = 0;
            for ( k = 1; k <= j; ++k )
                if ( ( j % k ) == 0 )
                s2 += k * k;
            pp[nn] = pp[nn] + pp[nn-j] * s2;
        }
        pp[nn] /= nn;
    }

    value = pp[n];
    free ( pp );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pyramid_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    PYRAMID_NUM returns the N-th pyramidal number.
  Discussion:
    The N-th pyramidal number P(N) is formed by the sum of the first
    N triangular numbers T(J):
      T(J) = sum ( 1 <= J <= N ) J
      P(N) = sum ( 1 <= I <= N ) T(I)
    By convention, T(0) = 0.
    The formula is:
      P(N) = ( (N+1)^3 - (N+1) ) / 6
  First Values:
      0
      1
      4
     10
     20
     35
     56
     84
    120
    165
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the index of the desired number, which must be
    at least 0.
    Output, int PYRAMID_NUM, the N-th pyramidal number.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ( ( n + 1 ) * ( n + 1 ) * ( n + 1 ) - ( n + 1 ) ) / 6;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pyramid_square_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    PYRAMID_SQUARE_NUM returns the N-th pyramidal square number.
  Discussion:
    The N-th pyramidal square number PS(N) is formed by the sum of the first
    N squares S:
      S(I) = I^2
      PS(N) = sum ( 1 <= I <= N ) S(I)
    By convention, PS(0) = 0.
    The formula is:
      PS(N) = ( N * ( N + 1 ) * ( 2*N+1 ) ) / 6
    Note that geometrically, this pyramid will have a square base.
  Example:
     0    0
     1    1
     2    5
     3   14
     4   30
     5   55
     6   91
     7  140
     8  204
     9  285
    10  385
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 August 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the index.
    0 <= N.
    Output, int PYRAMID_SQUARE_NUM, the N-th pyramidal square number.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ( n * ( n + 1 ) * ( 2 * n + 1 ) ) / 6;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_agm ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_AGM computes the arithmetic-geometric mean (AGM) of A and B.
  Discussion:
    The AGM is defined for nonnegative A and B.
    The AGM of numbers A and B is defined by setting
      A(0) = A,
      B(0) = B
      A(N+1) = ( A(N) + B(N) ) / 2
      B(N+1) = sqrt ( A(N) * B(N) )
    The two sequences both converge to AGM(A,B).
    In Mathematica, the AGM can be evaluated by
      ArithmeticGeometricMean [ a, b ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 February 2010
  Author:
    John Burkardt
  Reference:
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input, double A, B, the arguments whose AGM is to be computed.
    Output, double R8_AGM, the arithmetic-geometric mean of A and B.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	register ityp a = a_data[0];
	register ityp b = a_data[1];
	
    ityp c;
    ityp d;
    dim_typ it;
    const dim_typ it_max = 1000;
    ityp tol;

    if ( a < 0.00 || b < 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( a == 0.00 || b == 0.00 )
    {
    	result = 0.00;
        return &result;
    }
    
    it = 0;
    tol = 100.00 * r8_epsilon ( );

    for ( ; ; )
    {
        ++ it;

        c = ( a + b ) / 2.00;
        d = sqrt ( a * b );

        if ( fabs ( c - d ) <= tol * ( c + d ) ||  it_max < it )
            break;

        a = c;
        b = d;
    }

	result = c;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_factorial_log ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_FACTORIAL_LOG computes the natural logarithm of the factorial function.
  Discussion:
    LOG ( FACTORIAL ( N ) )
      = LOG ( product ( 1 <= I <= N ) I )
      = sum ( ( 1 <= I <= N ) LOG ( I ) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the argument of the factorial function.
    If N is less than 1, R8_FACTORIAL_LOG is returned as 0.
    Output, double R8_FACTORIAL_LOG, the logarithm of the factorial of N.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
    ityp value =  0.00;
    for (dim_typ i = 1; i <= n; ++i )
        value += log ( ( ityp ) i );
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_hyper_2f1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_HYPER_2F1 evaluates the hypergeometric function 2F1(A,B,C,X).
  Discussion:
    A minor bug was corrected.  The HW variable, used in several places as
    the "old" value of a quantity being iteratively improved, was not
    being initialized.  JVB, 11 February 2008.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 July 2009
  Author:
    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
    C version by John Burkardt.
    The original FORTRAN77 version of this routine is copyrighted by
    Shanjie Zhang and Jianming Jin.  However, they give permission to
    incorporate this routine into a user program provided that the copyright
    is acknowledged.
  Reference:
    Shanjie Zhang, Jianming Jin,
    Computation of Special Functions,
    Wiley, 1996,
    ISBN: 0-471-11963-6,
    LC: QA351.C45
  Parameters:
    Input, double A, B, C, X, the arguments of the function.
    C must not be equal to a nonpositive integer.
    X < 1.
    Output, double R8_HYPER_2F1, the value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	register ityp a = a_data[0];
	register ityp b = a_data[1];
	register ityp c = a_data[2];
	register ityp x = a_data[3];
	
    ityp a0;
    ityp aa;
    ityp bb;
    ityp c0;
    ityp c1;
    ityp el = 0.5772156649015329;
    ityp eps;
    ityp f0;
    ityp f1;
    ityp g0;
    ityp g1;
    ityp g2;
    ityp g3;
    ityp ga;
    ityp gabc;
    ityp gam;
    ityp gb;
    ityp gbm;
    ityp gc;
    ityp gca;
    ityp gcab;
    ityp gcb;
    ityp gm;
    ityp hf;
    ityp hw;
    int j;
    int k;
    int l0;
    int l1;
    int l2;
    int l3;
    int l4;
    int l5;
    int m;
    int nm;
    ityp pa;
    ityp pb;
    ityp r;
    ityp r0;
    ityp r1;
    ityp rm;
    ityp rp;
    ityp sm;
    ityp sp;
    ityp sp0;
    ityp x1;

    l0 = ( c == ( int ) ( c ) ) && ( c < 0.00 );
    l1 = ( 1.00 - x < 1.0E-15 ) && ( c - a - b <= 0.00 );
    l2 = ( a == ( int ) ( a ) ) && ( a < 0.00 );
    l3 = ( b == ( int ) ( b ) ) && ( b < 0.00 );
    l4 = ( c - a == ( int ) ( c - a ) ) && ( c - a <= 0.00 );
    l5 = ( c - b == ( int ) ( c - b ) ) && ( c - b <= 0.00 );

    if ( l0 || l1 )
    {
    	result = MAX_VAL;
        return &result;
    }

    eps = 0.95 < x ? 1.0E-08 : 1.0E-15;

    if ( x == 0.00 || a == 0.00 || b == 0.00 )
    {
    	result = 1.00;
        return &result;
    }
    else if ( 1.00 - x == eps && 0.00 < c - a - b )
    {
    	result = tgamma ( c ) * tgamma ( c - a - b ) / ( tgamma ( c - a ) * tgamma ( c - b ) );
        return &result;
    }
    else if ( 1.00 + x <= eps && fabs ( c - a + b - 1.00 ) <= eps )
    {
    	result = sqrt ( M_PI ) * pow ( 2.00, - a ) * tgamma ( c ) / ( tgamma ( 1.00 + a / 2.00 - b ) * tgamma ( 0.50 + 0.50 * a ) );
        return &result;
    }
    else if ( l2 || l3 )
    {
        if ( l2 )
            nm = ( int ) ( fabs ( a ) );

        if ( l3 )
            nm = ( int ) ( fabs ( b ) );


        hf = r = 1.00;

        for ( k = 1; k <= nm; ++k )
        {
            r *= ( a + k - 1.00 ) * ( b + k - 1.00 )/ ( k * ( c + k - 1.00 ) ) * x;
            hf += r;
        }

		result = hf;
        return &result;
    }
    else if ( l4 || l5 )
    {
        if ( l4 )
            nm = ( int ) ( fabs ( c - a ) );

        if ( l5 )
        nm = ( int ) ( fabs ( c - b ) );

        hf = r = 1.00;
        for ( k = 1; k <= nm; ++k )
        {
            r *= ( c - a + k - 1.0 ) * ( c - b + k - 1.0 )/ ( k * ( c + k - 1.0 ) ) * x;
            hf += r;
        }
        
        result = pow ( 1.00 - x, c - a - b ) * hf;
        return &result;
    }

    aa = a;
    bb = b;
    x1 = x;

    if ( x < 0.00 )
    {
        x /= ( x - 1.00 );
        if ( a < c && b < a && 0.0 < b )
        {
            a = bb;
            b = aa;
        }
        b = c - b;
    }

    if ( 0.75 <= x )
    {
        gm = 0.0;

        if ( fabs ( c - a - b - ( int ) ( c - a - b ) ) < 1.0E-15 )
        {
            m = ( int ) ( c - a - b );
            ga = tgamma ( a );
            gb = tgamma ( b );
            gc = tgamma ( c );
            gam = tgamma ( a + m );
            gbm = tgamma ( b + m );

            pa = r8_psi ( a );
            pb = r8_psi ( b );

            if ( m != 0 )
            gm = 1.00;

            for ( j = 1; j <= abs ( m ) - 1; ++j )
                gm *= j;

            rm = 1.00;
            for ( j = 1; j <= abs ( m ); ++j )
                rm *= j;

            f0 = r0 = r1 = 1.00;
            sp0 = sp = 0.00;

            if ( 0 <= m )
            {
                c0 = gm * gc / ( gam * gbm );
                c1 = - gc * pow ( x - 1.0, m ) / ( ga * gb * rm );

                for ( k = 1; k <= m - 1; ++k)
                {
                    r0 *= ( a + k - 1.0 ) * ( b + k - 1.0 )/ ( k * ( k - m ) ) * ( 1.0 - x );
                    f0 += r0;
                }

                for ( k = 1; k <= m; ++k )
                    sp0 += 1.00 / ( a + k - 1.00 ) + 1.00 / ( b + k - 1.00 )- 1.0 / ( ityp ) ( k );

                f1 = pa + pb + sp0 + 2.00 * el + log ( 1.00 - x );
                hw = f1;

                for ( k = 1; k <= 250; ++k )
                {
                    sp += ( 1.00 - a ) / ( k * ( a + k - 1.00 ) )+ ( 1.00 - b ) / ( k * ( b + k - 1.00 ) );
                    sm = 0.00;
                    for ( j = 1; j <= m; ++j )
                        sm += ( 1.00 - a )/ ( ( j + k ) * ( a + j + k - 1.00 ) )+ 1.00 / ( b + j + k - 1.00 );

                    rp = pa + pb + 2.0 * el + sp + sm + log ( 1.0 - x );

                    r1 *= ( a + m + k - 1.0 ) * ( b + m + k - 1.0 )/ ( k * ( m + k ) ) * ( 1.0 - x );
                    f1 += r1 * rp;

                    if ( fabs ( f1 - hw ) < fabs ( f1 ) * eps )
                        break;
                    hw = f1;
                }
                hf = f0 * c0 + f1 * c1;
            }
            else if ( m < 0 )
            {
                m *= -1;
                c0 = gm * gc / ( ga * gb * pow ( 1.00 - x, m ) );
                c1 = - pow ( - 1.00, m ) * gc / ( gam * gbm * rm );

                for ( k = 1; k <= m - 1; ++k )
                {
                    r0 *= ( a - m + k - 1.0 ) * ( b - m + k - 1.00 )/ ( k * ( k - m ) ) * ( 1.0 - x );
                    f0 += r0;
                }

                for ( k = 1; k <= m; ++k)
                    sp0 += 1.00 / ( ityp ) ( k );

                f1 = pa + pb - sp0 + 2.0 * el + log ( 1.0 - x );
                hw = f1;

                for ( k = 1; k <= 250; ++k )
                {
                    sp += ( 1.00 - a )/ ( k * ( a + k - 1.00 ) )+ ( 1.00 - b ) / ( k * ( b + k - 1.00 ) );
                    sm = 0.00;
                    for ( j = 1; j <= m; ++j )
                        sm += 1.00 / ( ityp ) ( j + k );

                    rp = pa + pb + 2.00 * el + sp - sm + log ( 1.00 - x );
                    r1 *= ( a + k - 1.00 ) * ( b + k - 1.00 )/ ( k * ( m + k ) ) * ( 1.00 - x );

                    f1 += r1 * rp;

                    if ( fabs ( f1 - hw ) < fabs ( f1 ) * eps )
                        break;
                    hw = f1;
                }
                hf = f0 * c0 + f1 * c1;
            }
        }
        else
        {
            ga = tgamma ( a );
            gb = tgamma ( b );
            gc = tgamma ( c );
            gca = tgamma ( c - a );
            gcb = tgamma ( c - b );
            gcab = tgamma ( c - a - b );
            gabc = tgamma ( a + b - c );
            c0 = gc * gcab / ( gca * gcb );
            c1 = gc * gabc / ( ga * gb ) * pow ( 1.00 - x, c - a - b );
            hf = 0.00;
            hw = hf;
            r0 = c0;
            r1 = c1;

            for ( k = 1; k <= 250; ++k )
            {
                r0 *= ( a + k - 1.00 ) * ( b + k - 1.00 )/ ( k * ( a + b - c + k ) ) * ( 1.00 - x );
                r1 *= ( c - a + k - 1.00 ) * ( c - b + k - 1.00 )/ ( k * ( c - a - b + k ) ) * ( 1.00 - x );
                hf += r0 + r1;

                if ( fabs ( hf - hw ) < fabs ( hf ) * eps )
                    break;
                hw = hf;
            }
            hf += c0 + c1;
        }
    }
    else
    {
        a0 = 1.00;

        if ( a < c && c < 2.00 * a && b < c && c < 2.00 * b )
        {
            a0 = pow ( 1.00 - x, c - a - b );
            a = c - a;
            b = c - b;
        }

        hf = r = 1.00;
        hw = hf;

        for ( k = 1; k <= 250; ++k )
        {
            r *= ( a + k - 1.00 ) * ( b + k - 1.00 )/ ( k * ( c + k - 1.00 ) ) * x;
            hf += r;

            if ( fabs ( hf - hw ) <= fabs ( hf ) * eps )
                break;

            hw = hf;
        }
        hf = a0 * hf;
    }

    if ( x1 < 0.00 )
    {
        x = x1;
        c0 = 1.00 / pow ( 1.00 - x, aa );
        hf = c0 * hf;
    }

    a = aa;
    b = bb;
    
    result = 120 < k ? MAX_VAL:hf;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_psi ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_PSI evaluates the function Psi(X).
  Discussion:
    This routine evaluates the logarithmic derivative of the
    Gamma function,
      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
             = d/dX LN ( GAMMA(X) )
    for real X, where either
      - XMAX1 < X < - XMIN, and X is not a negative integer,
    or
      XMIN < X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 February 2008
  Author:
    Original FORTRAN77 version by William Cody.
    C version by John Burkardt.
  Reference:
    William Cody, Anthony Strecok, Henry Thacher,
    Chebyshev Approximations for the Psi Function,
    Mathematics of Computation,
    Volume 27, Number 121, January 1973, pages 123-127.
  Parameters:
    Input, double XX, the argument of the function.
    Output, double R8_PSI, the value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp xx = *(ityp *) data;
	
    ityp aug;
    ityp den;
    dim_typ i;
    dim_typ n;
    dim_typ nq;
    ityp p1[9] =
    {
        4.5104681245762934160E-03,
        5.4932855833000385356,
        3.7646693175929276856E+02,
        7.9525490849151998065E+03,
        7.1451595818951933210E+04,
        3.0655976301987365674E+05,
        6.3606997788964458797E+05,
        5.8041312783537569993E+05,
        1.6585695029761022321E+05
    };
    ityp p2[7] =
    {
        -2.7103228277757834192,
        -1.5166271776896121383E+01,
        -1.9784554148719218667E+01,
        -8.8100958828312219821,
        -1.4479614616899842986,
        -7.3689600332394549911E-02,
        -6.5135387732718171306E-21
    };
    ityp piov4 = 0.78539816339744830962;
    ityp q1[8] =
    {
        9.6141654774222358525E+01,
        2.6287715790581193330E+03,
        2.9862497022250277920E+04,
        1.6206566091533671639E+05,
        4.3487880712768329037E+05,
        5.4256384537269993733E+05,
        2.4242185002017985252E+05,
        6.4155223783576225996E-08
    };
    ityp q2[6] =
    {
        4.4992760373789365846E+01,
        2.0240955312679931159E+02,
        2.4736979003315290057E+02,
        1.0742543875702278326E+02,
        1.7463965060678569906E+01,
        8.8427520398873480342E-01
    };
    ityp sgn;
    ityp upper;
    ityp value;
    ityp w;
    ityp x;
    ityp x01 = 187.00;
    ityp x01d = 128.00;
    ityp x02 = 6.9464496836234126266E-04;
    ityp xinf = 1.70E+38;
    ityp xlarge = 2.04E+15;
    ityp xmax1 = 3.60E+16;
    ityp xmin1 = 5.89E-39;
    ityp xsmall = 2.05E-09;
    ityp z;

    x = xx;
    w = fabs ( x );
    aug = 0.00;
    /*
    Check for valid arguments, then branch to appropriate algorithm.
    */
    if ( xmax1 <= - x || w < xmin1 )
    {
    	result = 0.00<x?-xinf:xinf;
        return &result;
    }

    if ( x < 0.5 )
    {
        /*
        X < 0.5, use reflection formula: psi(1-x) = psi(x) + M_PI * cot(M_PI*x)
        Use 1/X for M_PI*COTAN(M_PI*X)  when  XMIN1 < |X| <= XSMALL.
        */
        if ( w <= xsmall )
        aug = - 1.00 / x;
        /*
        Argument reduction for cotangent.
        */
        else
        {
            sgn = x<0.00 ? piov4:-piov4;
            w = w - ( ityp ) ( ( int ) ( w ) );
            nq = ( int ) ( w * 4.00 );
            w = 4.0 * ( w - ( ityp ) ( nq ) * 0.25 );
            /*
            W is now related to the fractional part of 4.0 * X.
            Adjust argument to correspond to values in the first
            quadrant and determine the sign.
            */
            n = nq >> 1;

            if ( n + n != nq )
                w = 1.00 - w;

            z = piov4 * w;

            if ( ( n % 2 ) != 0 )
                sgn *= -1;
            /*
            Determine the final value for  -M_PI * cotan(M_PI*x).
            */
            n = ( nq + 1 ) / 2;
            if ( ( n % 2 ) == 0 )
            {
                /*
                Check for singularity.
                */
                if ( z == 0.0 )
                {
                	result = 0.00<x?-xinf:xinf;
                    return &result;
                }
                
                aug = sgn * ( 4.00 / tan ( z ) );
            }
            else
                aug = sgn * ( 4.00 * tan ( z ) );
        }
        x = 1.00 - x;
    }
    /*
    0.5 <= X <= 3.0.
    */
    if ( x <= 3.00 )
    {
        den = x;
        upper = p1[0] * x;
        for ( i = 1; i <= 7; ++i )
        {
            den = ( den + q1[i-1] ) * x;
            upper = ( upper + p1[i]) * x;
        }
        
        result = ( upper + p1[8] ) / ( den + q1[7] ) * ( x - x01 / x01d ) - x02 + aug;
        return &result;
    }
    /*
    3.0 < X.
    */
    if ( x < xlarge )
    {
        w = 1.00 / ( x * x );
        den = w;
        upper = p2[0] * w;
        for ( i = 1; i <= 5; ++i )
        {
            den = ( den + q2[i-1] ) * w;
            upper = ( upper + p2[i] ) * w;
        }
        aug = ( upper + p2[6] ) / ( den + q2[5] ) - 0.50 / x + aug;
    }

	result = aug + log ( x ); 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8poly_value_horner ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
  Discussion:
    The polynomial
      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
    is to be evaluated at the value X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 January 2015
  Author:
    John Burkardt
  Parameters:
    Input, int M, the degree of the polynomial.
    Input, double C[M+1], the coefficients of the polynomial.
    C[0] is the constant term.
    Input, double X, the point at which the polynomial is to be evaluated.
    Output, double R8POLY_VALUE_HORNER, the value of the polynomial at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data;
	const register dim_typ m = s_data->a0;
	ityp * c = s_data->a1;
	const register ityp x = s_data->a2;
	
    ityp value = c[m];
    for (dim_typ i = m - 1; 0 <= i; --i )
        value *= x + c[i];

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sigma ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIGMA returns the value of SIGMA(N), the divisor sum.
  Discussion:
    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
    The formula is:
      SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
      SIGMA(P^K) = ( P^(K+1) - 1 ) / ( P - 1 ) if P is prime.
  First values:
     N  SIGMA(N)
     1    1
     2    3
     3    4
     4    7
     5    6
     6   12
     7    8
     8   15
     9   13
    10   18
    11   12
    12   28
    13   14
    14   24
    15   24
    16   31
    17   18
    18   39
    19   20
    20   42
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the value to be analyzed.
    Output, int SIGMA, the value of SIGMA(N).  If N is less than
    or equal to 0, SIGMA will be returned as 0.  If there is not
    enough room for factoring N, SIGMA is returned as -1.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
    # define FACTOR_MAX 20

    int factor[FACTOR_MAX];
    dim_typ i;
    int nfactor;
    int nleft;
    int power[FACTOR_MAX];
    int value;

    if ( n == 0)
    {
    	result = 0;
        return &result;
    }

    if ( n == 1 )
    {
    	result = 1;
        return &result;
    }
    /*
    Factor N.
    */
    i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

    if ( nleft != 1 )
    {
    	result = MAX_VAL;
        return &result;
    }

    value = 1;
    for ( i = 0; i < nfactor; ++i )
        value = ( value *( ( int ) pow ( ( ityp ) factor[i], power[i] + 1 ) - 1 ) )/ ( factor[i] - 1 );
        
    result = value;
    return &result;
    # undef FACTOR_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _simplex_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_NUM evaluates the N-th Simplex number in M dimensions.
  Discussion:
     N\M: 1    2    3    4    5
    --   --   --   --   --   --
     0    0    0    0    0    0
     1    1    1    1    1    1
     2    2    3    4    5    6
     3    3    6   10   15   21
     4    4   10   20   35   56
     5    5   15   35   70  126
     6    6   21   56  126  252
     7    7   28   84  210  462
     8    8   36  120  330  792
     9    9   45  165  495 1287
    10   10   55  220  715 2002
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 February 2015
  Author:
    John Burkardt
  Parameters
    Input, int M, the spatial dimension.
    Input, int N, the index of the number.
    Output, int SIMPLEX_NUM, the desired value.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	dim_typ m = a_data[0];
	dim_typ n = a_data[1];
	
    dim_typ value = 1;
    for (dim_typ i = 1; i <= m; ++i)
    	value = ( value * ( n + i - 1 ) ) / i;

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sin_power_int ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIN_POWER_INT evaluates the sine power integral.
  Discussion:
    The function is defined by
      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
    The algorithm uses the following fact:
      Integral sin^n ( t ) = (1/n) * (
        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 September 2004
  Author:
    John Burkardt
  Parameters
    Input, double A, B, the limits of integration.
    Input, integer N, the power of the sine function.
    Output, double SIN_POWER_INT, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const _2itdt * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    ityp ca;
    ityp cb;
    dim_typ m;
    dim_typ mlo;
    ityp sa;
    ityp sb;
    ityp value;

    sa = sin ( a );
    sb = sin ( b );
    ca = cos ( a );
    cb = cos ( b );

    if ( ( n % 2 ) == 0 )
    {
        value = b - a;
        mlo = 2;
    }
    else
    {
        value = ca - cb;
        mlo = 3;
    }

    for ( m = mlo; m <= n; m += 2 )
        value = ( ( ityp ) ( m - 1 ) * value+ pow ( sa, (m-1) ) * ca - pow ( sb, (m-1) ) * cb )/ ( ityp ) ( m );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _slice ( void * data)
/******************************************************************************/
/*
  Purpose:
    SLICE: maximum number of pieces created by a given number of slices.
  Discussion:
    If we imagine slicing a pizza, each slice produce more pieces.
    The position of the slice affects the number of pieces created, but there
    is a maximum.
    This function determines the maximum number of pieces created by a given
    number of slices, applied to a space of a given dimension.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 August 2011
  Author:
    John Burkardt
  Reference:
    Robert Banks,
    Slicing Pizzas, Racing Turtles, and Further Adventures in
    Applied Mathematics,
    Princeton, 1999,
    ISBN13: 9780691059471,
    LC: QA93.B358.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, int SLICE_NUM, the number of slices.
    Input, int SLICE, the maximum number of pieces that can
    be created by the given number of slices applied in the given dimension.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ dim_num = a_data[0];
	const register dim_typ slice_num = a_data[1];
	
    dim_typ piece_num = 0;
    for (dim_typ j = 0; j <= MIN ( dim_num, slice_num ); ++j )
        piece_num += i4_choose ( slice_num, j );
        
    result = piece_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _spherical_harmonic ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERICAL_HARMONIC evaluates spherical harmonic functions.
  Discussion:
    The spherical harmonic function Y(L,M,THETA,PHI,X) is the
    angular part of the solution to Laplace's equation in spherical
    coordinates.
    Y(L,M,THETA,PHI,X) is related to the associated Legendre
    function as follows:
      Y(L,M,THETA,PHI,X) = FACTOR * P(L,M,cos(THETA)) * exp ( i * M * PHI )
    Here, FACTOR is a normalization factor:
      FACTOR = sqrt ( ( 2 * L + 1 ) * ( L - M )! / ( 4 * M_PI * ( L + M )! ) )
    In Mathematica, a spherical harmonic function can be evaluated by
      SphericalHarmonicY [ l, m, theta, phi ]
    Note that notational tradition in physics requires that THETA
    and PHI represent the reverse of what they would normally mean
    in mathematical notation; that is, THETA goes up and down, and
    PHI goes around.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 March 2005
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 2002,
    Second edition,
    ISBN: 1584883472,
    LC: QA5.W45.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input, int L, the first index of the spherical harmonic function.
    Normally, 0 <= L.
    Input, int M, the second index of the spherical harmonic function.
    Normally, -L <= M <= L.
    Input, double THETA, the polar angle, for which
    0 <= THETA <= M_PI.
    Input, double PHI, the longitudinal angle, for which
    0 <= PHI <= 2*M_PI.
    Output, double C[L+1], S[L+1], the real and imaginary
    parts of the functions Y(L,0:L,THETA,PHI).
*/
{
	const dts2it2pit * const s_data = data;
	const register dim_typ l = s_data->a0;
	const register short m = s_data->a1;
	const register ityp theta = s_data->a2;
	ityp phi = s_data->a3;
	ityp * c = s_data->a4;
	ityp * s = s_data->a5;
	
    ityp angle;
    dim_typ i;
    int m_abs;
    ityp *plm;

    m_abs = abs ( m );

    plm = ( ityp * ) malloc ( ( l + 1 ) * sizeof ( ityp ) );
    legendre_associated_normalized ( l, m_abs, cos ( theta ), plm );
    angle = ( ityp ) ( m ) * phi;

    if ( 0 <= m )
        for ( i = 0; i <= l; ++i )
        {
            c[i] = plm[i] * cos ( angle );
            s[i] = plm[i] * sin ( angle );
        }
    else
        for ( i = 0; i <= l; ++i )
        {
            c[i] = -plm[i] * cos ( angle );
            s[i] = -plm[i] * sin ( angle );
        }

    free ( plm );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tau ( void * data)
/******************************************************************************/
/*
  Purpose:
    TAU returns the value of TAU(N), the number of distinct divisors of N.
  Discussion:
    TAU(N) is the number of divisors of N, including 1 and N.
    The formula is:
      If the prime factorization of N is
        N = P1^E1 * P2^E2 * ... * PM^EM,
      then
        TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
  First values:
     N   TAU(N)
     1    1
     2    2
     3    2
     4    3
     5    2
     6    4
     7    2
     8    4
     9    3
    10    4
    11    2
    12    6
    13    2
    14    4
    15    4
    16    5
    17    2
    18    6
    19    2
    20    6
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 December 1998
  Author:
    John Burkardt
  Parameters:
    Input, int N, the value to be analyzed.  N must be 1 or
    greater.
    Output, int TAU, the value of TAU(N).  But if N is 0 or
    less, TAU is returned as 0, a nonsense value.  If there is
    not enough room for factoring, TAU is returned as -1.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    # define FACTOR_MAX 20

    int factor[FACTOR_MAX];
    dim_typ i;
    int nfactor;
    int nleft;
    int power[FACTOR_MAX];
    int value;

    if ( n == 0)
    {
    	result = 0;
    }
        return &result;

    if ( n == 1 )
    {
    	result = 1;
        return &result;
    }
    /*
    Factor N.
    */
    i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

    if ( nleft != 1 )
    {
    	result = USHRT_MAX;
        return &result;
    }

	value = 1;
    for ( i = 0; i < nfactor; ++i )
   		value *= ( power[i] + 1 );

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NUM returns the N-th tetrahedral number.
  Discussion:
    The N-th tetrahedral number T3(N) is formed by the sum of the first
    N triangular numbers:
      T3(N) = sum ( 1 <= I <= N ) T2(I)
            = sum ( 1 <= I <= N ) sum ( 1 <= J < I ) J
    By convention, T3(0) = 0.
    The formula is:
      T3(N) = ( N * ( N + 1 ) * ( N + 2 ) ) / 6
  First Values:
     0
     1
     4
    10
    20
    35
    56
    84
   120
   165
   220
   275
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the index of the desired number, which must be
    at least 0.
    Output, int TETRAHEDRON_NUM, the N-th tetrahedron number.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ( n * ( n + 1 ) * ( n + 2 ) ) / 6;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NUM returns the N-th triangular number.
  Discussion:
    The N-th triangular number T(N) is formed by the sum of the first
    N integers:
      T(N) = sum ( 1 <= I <= N ) I
    By convention, T(0) = 0.
    The formula is:
      T(N) = ( N * ( N + 1 ) ) / 2
  First Values:
     0
     1
     3
     6
    10
    15
    21
    28
    36
    45
    55
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the index of the desired number, which must be
    at least 0.
    Output, int TRIANGLE_NUM, the N-th triangular number.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ( n * ( n + 1 ) ) >> 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _trinomial ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRINOMIAL computes a trinomial coefficient.
  Discussion:
    The trinomial coefficient is a generalization of the binomial
    coefficient.  It may be interpreted as the number of combinations of
    N objects, where I objects are of type 1, J of type 2, and K of type 3.
    and N = I + J + K.
    T(I,J,K) = N! / ( I! J! K! )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, K, the factors.
    All should be nonnegative.
    Output, int TRINOMIAL, the trinomial coefficient.
*/
{
	static int result = INT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ i = a_data[0];
	const register dim_typ j = a_data[1];
	const register dim_typ k = a_data[2];
	
    dim_typ l;
    int t;
    int value;
    /*
    Each factor must be nonnegative.
    */

    value = t = 1;

    for ( l = 1; l <= i; l++, ++t);

    for ( l = 1; l <= j; l++ )
    {
        value *= t / l;
        ++ t;
    }

    for ( l = 1; l <= k; l++ )
    {
        value *= t / l;
        ++ t;
    }

	result = value;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _v_hofstadter ( void * data)
/******************************************************************************/
/*
  Purpose:
    V_HOFSTADTER computes the Hofstadter V sequence.
  Discussion:
    V(N) = 0                          if N = 0
         = 1                          if 1 <= N <= 4
         = V(N-V(N-1)) + V(N-V(N-4)), otherwise.
    V(N) is defined for all nonnegative integers.
  Table:
     N  V(N)
    --  ----
     0     0
     1     1
     2     1
     3     1
     4     1
     5     2
     6     3
     7     4
     8     5
     9     5
    10     6
    11     6
    12     7
    13     8
    14     8
    15     9
    16     9
    17    10
    18    11
    19    11
    20    11
    21    12
    22    12
    23    13
    24    14
    25    14
    26    15
    27    15
    28    16
    29    17
    30    17
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the argument of the function.
    Output, int V_HOFSTADTER, the value of the function.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = n==0 ? 0: n<=4 ? 1 : (  v_hofstadter ( n - v_hofstadter(n-1) )+ v_hofstadter ( n - v_hofstadter(n-4) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _vibonacci ( void * data)
/******************************************************************************/
/*
  Purpose:
    VIBONACCI computes the first N Vibonacci numbers.
  Discussion:
    The "Vibonacci numbers" are a generalization of the Fibonacci numbers:
      V(N+1) = +/- V(N) +/- V(N-1)
    where the signs are chosen randomly.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Reference:
    Brian Hayes,
    The Vibonacci Numbers,
    American Scientist,
    July-August 1999, Volume 87, Number 4.
    Divakar Viswanath,
    Random Fibonacci sequences and the number 1.13198824,
    Mathematics of Computation, 1998.
  Parameters:
    Input, int N, the highest number to compute.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int V(N), the first N Vibonacci numbers.  By convention,
    V(1) and V(2) are taken to be 1.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	int * v = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    short s1;
    short s2;

    if ( n == 0 )
        return NULL;

    v[0] = 1;

    if ( n <= 1 )
        return NULL;

    v[1] = 1;

    for ( i = 2; i < n; ++i)
    {
        j = i4_uniform_ab ( 0, 1, seed );
        s1 = 1 - ((j==0)<<1);
        j = i4_uniform_ab ( 0, 1, seed );
        s2 = 1 - ((j==0)<<1);
        v[i] = s1 * v[i-1] + s2 * v[i-2];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _zeckendorf ( void * data)
/******************************************************************************/
/*
  Purpose:
    ZECKENDORF produces the Zeckendorf decomposition of a positive integer.
  Discussion:
    Zeckendorf proved that every positive integer can be represented
    uniquely as the sum of non-consecutive Fibonacci numbers.
    N = sum ( 1 <= I <= M ) F_LIST(I)
  Example:
     N    Decomposition
    50    34 + 13 + 3
    51    34 + 13 + 3 + 1
    52    34 + 13 + 5
    53    34 + 13 + 5 + 1
    54    34 + 13 + 5 + 2
    55    55
    56    55 + 1
    57    55 + 2
    58    55 + 3
    59    55 + 3 + 1
    60    55 + 5
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 July 2000
  Author:
    John Burkardt
  Parameters:
    Input, int N, the positive integer to be decomposed.
    Input, int M_MAX, the maximum number of parts in the decomposition.
    Set M_MAX = 100 to be safe.
    Output, int M, the number of parts in the decomposition.
    Output, int I_LIST[M], the index of the Fibonacci numbers
    in the decomposition.
    Output, int F_LIST[M], the value of the Fibonacci numbers
    in the decomposition.
*/
{
	const _2dt3pi * const s_data = data;
	register dim_typ n = s_data->a0;
	const register dim_typ m_max = s_data->a1;
	int * m = s_data->a2;
	int * i_list = s_data->a3;
	int * f_list = s_data->a4;
	
    int f;
    int i, j;

    *m = 0;
    /*
    Extract a sequence of Fibonacci numbers.
    */
    while ( 0 < n && *m < m_max )
    {
        fibonacci_floor ( n, &f, &i );

        i_list[*m] = i;
        ++ *m;
        n -= f;
    }
    /*
    Replace any pair of consecutive indices ( I, I-1 ) by I+1.
    */
    for ( i = *m; 2 <= i; --i)
    {
        if ( i_list[i-2] == i_list[i-1] + 1 )
        {
            ++ i_list[i-2];
            for ( j = i; j <= *m-1; j++ )
                i_list[j-1] = i_list[j];
            i_list[*m-1] = 0;
            -- *m;
        }

    }
    /*
    Fill in the actual values of the Fibonacci numbers.
    */
    for ( i = 0; i < *m; ++i )
        f_list[i] = fibonacci_direct ( i_list[i] );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _zernike_poly ( void * data)
/******************************************************************************/
/*
  Purpose:
    ZERNIKE_POLY evaluates a Zernike polynomial at RHO.
  Discussion:
    This routine uses the facts that:
    *) R^M_N = 0 if M < 0, or N < 0, or N < M.
    *) R^M_M = RHO^M
    *) R^M_N = 0 if mod ( N - M, 2 ) = 1.
    and the recursion:
    R^M_(N+2) = A * [ ( B * RHO^2 - C ) * R^M_N - D * R^M_(N-2) ]
    where
    A = ( N + 2 ) / ( ( N + 2 )^2 - M^2 )
    B = 4 * ( N + 1 )
    C = ( N + M )^2 / N + ( N - M + 2 )^2 / ( N + 2 )
    D = ( N^2 - M^2 ) / N
    I wish I could clean up the recursion in the code, but for
    now, I have to treat the case M = 0 specially.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2005
  Author:
    John Burkardt
  Reference:
    Eric Weisstein,
    CRC Concise Encyclopedia of Mathematics,
    CRC Press, 2002,
    Second edition,
    ISBN: 1584883472,
    LC: QA5.W45.
  Parameters:
    Input, int M, the upper index.
    Input, int N, the lower index.
    Input, double RHO, the radial coordinate.
    Output, double ZERNIKE_POLY, the value of the Zernike
    polynomial R^M_N at the point RHO.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register ityp rho = s_data->a2;
	
    ityp a;
    ityp b;
    ityp c;
    ityp d;
    dim_typ nn;
    ityp z;
    ityp zm2;
    ityp zp2;

    /*
    Do checks.
    */
    if ( ( n - m ) % 2 == 1 || n < m)
    {
    	result = 0.00;
        return &result;
    }

    zm2 = 0.0;
    z = pow ( rho, m );

    if ( m == 0 )
    {
    if ( n == 0 )
    {
    	result = z;
        return &result;
	}
	
    zm2 = z;
    z = 2.00 * rho * rho - 1.00;

    for ( nn = m+2; nn <= n-2; nn += 2 )
    {
    a = ( ityp ) ( nn + 2 )
    / ( ityp ) ( ( nn + 2 ) * ( nn + 2 ) - m * m );

    b = ( ityp ) ( ( nn + 1 ) << 2 );

    c = ( ityp ) ( ( nn + m ) * ( nn + m ) ) / ( ityp ) ( nn )+ ( ityp ) ( ( nn - m + 2 ) * ( nn - m + 2 ) )/ ( ityp ) ( nn + 2 );

    d = ( ityp ) ( nn * nn - m * m ) / ( ityp ) ( nn );

    zp2 = a * ( ( b * rho * rho - c ) * z - d * zm2 );
    zm2 = z;
    z = zp2;
    }
    }
    else
    {
        for ( nn = m; nn <= n-2; nn = nn + 2 )
        {
            a = ( ityp ) ( nn + 2 )/ ( ityp ) ( ( nn + 2 ) * ( nn + 2 ) - m * m );
            b = ( ityp ) ( ( nn + 1 ) << 2 );
            c = ( ityp ) ( ( nn + m ) * ( nn + m ) ) / ( ityp ) ( nn )+ ( ityp ) ( ( nn - m + 2 ) * ( nn - m + 2 ) )/ ( ityp ) ( nn + 2 );
            d = ( ityp ) ( nn * nn - m * m ) / ( ityp ) ( nn );
            zp2 = a * ( ( b * rho * rho - c ) * z - d * zm2 );
            zm2 = z;
            z = zp2;
        }
    }

	result = z;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _zeta ( void * data)
/******************************************************************************/
/*
  Purpose:
    ZETA estimates the Riemann Zeta function.
  Definition:
    For 1 < P, the Riemann Zeta function is defined as:
      ZETA ( P ) = Sum ( 1 <= N < +oo ) 1 / N**P
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 August 2004
  Author:
    John Burkardt
  Reference:
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.
  Parameters:
    Input, double P, the power to which the integers are raised.
    P must be greater than 1.
    Output, double ZETA, an approximation to the Riemann
    Zeta function.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp p = *(ityp *) data;
	
    dim_typ n;
    ityp value;
    ityp value_old;

    if ( p <= 1.00)
    {
    	result =  MAX_VAL;
        return &result;
    }

    value = 0.00;
    n = 0;

    for ( ; ; )
    {
        ++ n ;
        value_old = value;
        value += 1.00 / pow ( ( ityp ) n, p );

        if ( value <= value_old || 50000 <= n )
            break;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4_factor ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_FACTOR factors an integer into prime factors.
  Discussion:
    N = NLEFT * Product ( I = 1 to NFACTOR ) FACTOR(I)**POWER(I).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the integer to be factored.  N may be positive,
    negative, or 0.
    Input, int MAXFACTOR, the maximum number of prime factors for
    which storage has been allocated.
    Output, int *NFACTOR, the number of prime factors of N discovered
    by the routine.
    Output, int FACTOR[MAXFACTOR], the prime factors of N.
    Output, int POWER[MAXFACTOR].  POWER(I) is the power of
    the FACTOR(I) in the representation of N.
    Output, int *NLEFT, the factor of N that the routine could not
    divide out.  If NLEFT is 1, then N has been completely factored.
    Otherwise, NLEFT represents factors of N involving large primes.
*/
{
	const _2dt4pi * const s_data = data;
	register dim_typ n = s_data->a0;
	const register dim_typ maxfactor = s_data->a1;
	int * nfactor = s_data->a2;
	int * factor = s_data->a3;
	int * power = s_data->a4;
	int * nleft = s_data->a5;
	
    dim_typ i;
    int maxprime;
    int p;

    *nfactor = 0;

    for ( i = 0; i < maxfactor; ++i )
        factor[i] = power[i] = 0;

    *nleft = n;

    if ( n == 0 )
        return NULL;

    if ( abs ( n ) == 1 )
    {
        *nfactor = factor[0] = power[0] = 1;
        return NULL;
    }
    /*
    Find out how many primes we stored.
    */
    maxprime = prime ( -1 );
    /*
    Try dividing the remainder by each prime.
    */
    for ( i = 1; i <= maxprime; ++i )
    {
        p = prime ( i );

        if ( abs ( *nleft ) % p == 0 )
        {
            if ( *nfactor < maxfactor )
            {
                ++ *nfactor;
                factor[*nfactor-1] = p;

                for ( ; ; )
                {
                    ++ power[*nfactor-1];
                    *nleft /= p;

                    if ( abs ( *nleft ) % p != 0 )
                        break;

                }
                if ( abs ( *nleft ) == 1 )
                    break;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fibonacci_direct ( void * data)
/******************************************************************************/
/*
  Purpose:
    FIBONACCI_DIRECT computes the N-th Fibonacci number directly.
  Discussion:
    The formula is:
      F(N) = ( PHIP**N - PHIM**N ) / sqrt(5)
    where
      PHIP = ( 1 + sqrt(5) ) / 2,
      PHIM = ( 1 - sqrt(5) ) / 2.
  Example:
     N   F
    --  --
     0   0
     1   1
     2   1
     3   2
     4   3
     5   5
     6   8
     7  13
     8  21
     9  34
    10  55
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the index of the Fibonacci number to compute.
    N should be nonnegative.
    Output, int FIBONACCI_DIRECT, the value of the N-th Fibonacci number.
*/
{
	static int result = INT_MAX;
	
	const register int n = *(int *) data;
	
	result = 0 + (n>=0)*r8_nint ( ( pow ( ( 1.00 + 2.236068 ) / 2.00, n ) - pow ( ( 1.00 - 2.236068 ) / 2.00, n ) ) / sqrt ( 5.00 ) );
    return &result;
}

#endif
