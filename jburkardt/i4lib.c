#ifndef __DISABLEDEEP_I4LIB

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_abs ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_ABS returns the absolute value of an I4.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    John Burkardt
  Parameters:
    Input, int I, an integer.
    Output, int I4_ABS, the absolute value of the integer.
*/
{
	static int result = INT_MAX;
	
	const register int i = *(int *) data;
	
	result = 0<=i ? i:-i;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_mach ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_MACH returns integer machine constants.
  Discussion:
    Input/output unit numbers.
      I4_MACH(1) = the standard input unit.
      I4_MACH(2) = the standard output unit.
      I4_MACH(3) = the standard punch unit.
      I4_MACH(4) = the standard error message unit.
    Words.
      I4_MACH(5) = the number of bits per integer storage unit.
      I4_MACH(6) = the number of characters per integer storage unit.
    Integers.
    Assume integers are represented in the S digit base A form:
      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
    where 0 <= X(1:S-1) < A.
      I4_MACH(7) = A, the base.
      I4_MACH(8) = S, the number of base A digits.
      I4_MACH(9) = A^S-1, the largest integer.
    Floating point numbers
    Assume floating point numbers are represented in the T digit
    base B form:
      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
      I4_MACH(10) = B, the base.
    Single precision
      I4_MACH(11) = T, the number of base B digits.
      I4_MACH(12) = EMIN, the smallest exponent E.
      I4_MACH(13) = EMAX, the largest exponent E.
    Double precision
      I4_MACH(14) = T, the number of base B digits.
      I4_MACH(15) = EMIN, the smallest exponent E.
      I4_MACH(16) = EMAX, the largest exponent E.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 April 2007
  Author:
    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
    C version by John Burkardt.
  Reference:
    Phyllis Fox, Andrew Hall, Norman Schryer,
    Algorithm 528,
    Framework for a Portable Library,
    ACM Transactions on Mathematical Software,
    Volume 4, Number 2, June 1978, page 176-188.
  Parameters:
    Input, int I, chooses the parameter to be returned.
    1 <= I <= 16.
    Output, int I4_MACH, the value of the chosen parameter.
*/
{
	static int result = INT_MAX;
	
	const register int i = *(int *) data;
	
    int value;

    if ( i < 1 || i > 16 )
    {
    	result = INT_MAX;
        return &result;
	}

    switch(i)
    {
        case 1:
            value = 5;
            break;
        case 2:
		case 4:
            value = 6;
            break;
        case 3:
            value = 7;
            break;
        case 5:
            value = 32;
            break;
        case 6:
            value = 4;
            break;
        case 7:
        case 10:
            value = 2;
            break;
        case 8:
            value = 31;
            break;
        case 9:
            value = 2147483647;
            break;
        case 11:
            value = 24;
            break;
        case 12:
            value = -125;
            break;
        case 13:
            value = 128;
            break;
        case 14:
            value = 53;
            break;
        case 15:
            value = -1021;
            break;
        case 16:
            value = 1024;
            break;
    }

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_pow ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_POW returns the value of I^J.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2004
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the base and the power.  J should be nonnegative.
    Output, int I4_POW, the value of I^J.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int i = a_data[0];
	const register int j = a_data[1];
	
    int k;
    int value;

    if ( j < 0 )
    {
        if ( i == 1 )
            value = 1;
        else if ( i == 0 )
        {
        	result = INT_MAX;
            return &result;
        }
        else
            value = 0;
    }
    else if ( j == 0 )
    {
        if(!i)
        {
        	result = INT_MAX;
            return &result;
        }
        else
            value = 1;
    }
    else if ( j == 1 )
        value = i;
    else
    {
        value = 1;
        for ( k = 1; k <= j; ++k )
            value *= i;
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_bit_hi1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_BIT_HI1 returns the position of the high 1 bit base 2 in an I4.
  Example:
       N    Binary    Hi 1
    ----    --------  ----
       0           0     0
       1           1     1
       2          10     2
       3          11     2
       4         100     3
       5         101     3
       6         110     3
       7         111     3
       8        1000     4
       9        1001     4
      10        1010     4
      11        1011     4
      12        1100     4
      13        1101     4
      14        1110     4
      15        1111     4
      16       10000     5
      17       10001     5
    1023  1111111111    10
    1024 10000000000    11
    1025 10000000001    11
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 January 2007
  Author:
    John Burkardt.
  Parameters:
    Input, int N, the integer to be measured.
    N should be nonnegative.  If N is nonpositive, I4_BIT_HI1
    will always be 0.
    Output, int I4_BIT_HI1, the location of the high order bit.
*/
{
	static int result = INT_MAX;
	
	register int n = *(int *) data;
	
    int value = 0;
    while ( 0 < n )
    {
        ++ value;
        n /= 2;
    }
    
    result = value;
    return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_bit_lo0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_BIT_LO0 returns the position of the low 0 bit base 2 in an I4.
  Example:
       N    Binary    Lo 0
    ----    --------  ----
       0           0     1
       1           1     2
       2          10     1
       3          11     3
       4         100     1
       5         101     2
       6         110     1
       7         111     4
       8        1000     1
       9        1001     2
      10        1010     1
      11        1011     3
      12        1100     1
      13        1101     2
      14        1110     1
      15        1111     5
      16       10000     1
      17       10001     2
    1023  1111111111     1
    1024 10000000000     1
    1025 10000000001     1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int N, the integer to be measured.
    N should be nonnegative.
    Output, int I4_BIT_LO0, the position of the low 1 bit.
*/
{
	static int result = INT_MAX;
	
	register int n = *(int *) data;
	
    int bit = 0;
    int n2;
    while ( 1 )
    {
        ++ bit;
        n2 = n / 2;

        if ( n == (n2<<1))
            break;

        n = n2;
    }
    
    result = bit;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_bit_lo1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_BIT_LO1 returns the position of the low 1 bit base 2 in an I4.
  Example:
       N    Binary    Lo 1
    ----    --------  ----
       0           0     0
       1           1     1
       2          10     2
       3          11     1
       4         100     3
       5         101     1
       6         110     2
       7         111     1
       8        1000     4
       9        1001     1
      10        1010     2
      11        1011     1
      12        1100     3
      13        1101     1
      14        1110     2
      15        1111     1
      16       10000     5
      17       10001     1
    1023  1111111111     1
    1024 10000000000    11
    1025 10000000001     1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int N, the integer to be measured.
    N should be nonnegative.
    Output, int I4_BIT_LO1, the position of the low 1 bit.
*/
{
	static int result = INT_MAX;
	
	const register int n = *(int *) data;
	
    int bit = 0;
    int i = n;
    int i2;

    for ( ; ; )
    {
        ++ bit;
        i2 = i / 2;
        if ( i != (i2<<1) )
        break;
        i = i2;
    }

	result = bit;
    return &result;
}

/*******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_bit_reverse ( void * data)
/*******************************************************************************/
/*
  Purpose:
    I4_BIT_REVERSE reverses the bits in an I4.
  Discussion:
    An I4 is an int value.
  Example:
       I      N  2^N     I4_BIT_REVERSE ( I, N )
    ----    --------  -----------------------
       0      0    1     0
       1      0    1     1
       0      3    8     0
       1      3    8     4
       2      3    8     2
       3      3    8     6
       4      3    8     1
       5      3    8     5
       6      3    8     3
       7      3    8     7
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 March 2008
  Author:
    John Burkardt
  Parameters:
    Input, int I, the integer to be bit reversed.
    I should be nonnegative.  Normally I < 2^N.
    Input, int N, indicates the number of bits to
    be reverse (N+1) or the base with respect to which the integer is to
    be reversed (2^N).  N should be nonnegative.
    Output, int I4_BIT_REVERSE, the bit reversed value.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int i = a_data[0];
	const register int n = a_data[1];
	
    int b;
    int j;
    int value;

    value = -1 + ((i<0)<<1);

    if(i >= 0 || n != 0)
    {
        b = powi ( 2, n );
        j = ( i % b );

        value = 0;

        for ( ; ; )
        {
            if ( b == 1 )
            {
                value = value + j;
                j = 0;
                break;
            }
            else
            {
                if ( ( j % 2 ) == 1 )
                {
                    value += b / 2;
                    -- j;
           	 	}
            }

            j >>= 1;
            b >>= 1;
        }
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_characteristic ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_CHARACTERISTIC gives the characteristic for an I4.
  Discussion:
    For any positive integer Q, the characteristic is:
    Q, if Q is a prime;
    P, if Q = P^N for some prime P and some integer N;
    0, otherwise, that is, if Q is negative, 0, 1, or the product
       of more than one distinct prime.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    C version by John Burkardt
  Reference:
    Pa__MATHSUITE __JBURKARDT  int  ul Bratley, Bennett Fox, Harald Niederreiter,
    Algorithm 738:
    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
    ACM Transactions on Mathematical Software,
    Volume 20, Number 4, pages 494-495, 1994.
  Parameters:
    Input, int Q, the value to be tested.
    Output, int I4_CHARACTERISTIC, the characteristic of Q.
*/
{
	static int result = INT_MAX;
	
	register int q = *(int *) data;
	
    int i;
    int i_maximum;

    if ( q <= 1 )
    {
    	result = 0;
        return &result;
    }
    /*
    If Q is not prime, then there is at least one prime factor
    of Q no greater than SQRT(Q)+1.

    A faster code would only consider prime values of I,
    but that entails storing a table of primes and limiting the
    size of Q.  Simplicity and flexibility for now!
    */
    i_maximum = ( int ) ( sqrt ( ( ityp ) ( q ) ) ) + 1;

    for ( i = 2; i <= i_maximum; ++i)
    {
        if ( ( q % i ) == 0 )
        {
            while ( ( q % i ) == 0 )
                q /= i;

			result = 0 + (q==1)*i;
			return &result;
        }
    }
    /*
    If no factor was found, then Q is prime.
    */
    
    result = q;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_choose ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_CHOOSE computes the binomial coefficient C(N,K).
  Discussion:
    The value is calculated in such a way as to avoid overflow and
    roundoff.  The calculation is done in integer arithmetic.
    The formula used is:
      C(N,K) = N! / ( K! * (N-K)! )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 December 2013
  Author:
    John Burkardt
  Reference:
    ML Wolfson, HV Wright,
    Algorithm 160:
    Combinatorial of M Things Taken N at a Time,
    Communications of the ACM,
    Volume 6, Number 4, April 1963, page 161.
  Parameters:
    Input, int N, K, are the values of N and K.
    Output, int I4_CHOOSE, the number of combinations of N
    things taken K at a time.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int n = a_data[0];
	const register int k = a_data[1];
	
    dim_typ i;
    int mn;
    int mx;
    int value;

    if ( k < n - k )
    {
        mn = k;
        mx = n - k;
    }
    else
    {
        mn -= k;
        mx = k;
    }

    if(mn < 0)
    	value = 0;
    else if (mn == 0)
    	value = 1;
	else
    {
        value = mx + 1;
        for ( i = 2; i <= mn; ++i )
            value = ( value * ( mx + i ) ) / i;
    }

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_div_rounded ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_DIV_ROUNDED computes the rounded result of I4 division.
  Discussion:
    This routine computes C = A / B, where A, B and C are integers
    and C is the closest integer value to the exact real result.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, int A, B, the number to be divided,
    and the divisor.
    Output, int I4_DIV_ROUNDED, the rounded result
    of the division.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int a = a_data[0];
	const register int b = a_data[1];
	
    int a_abs;
    int b_abs;
    int value;

    if ( a == 0 && b == 0 )
        value = i4_huge;
    else if ( a == 0 )
        value = 0;
    else if ( b == 0 )
        value = a<0 ? -i4_huge:i4_huge;
    else
    {
        a_abs = abs ( a );
        b_abs = abs ( b );

        value = a_abs / b_abs;
        /*
        Round the value.
        */
        if ( ( (value<<1) + 1 ) * b_abs < (a_abs<<1) )
            ++ value;
    /*
    Set the sign.
    */
        if ( ( a < 0 && 0 < b ) || ( 0 < a && b < 0 ) )
            value *= -1;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_division ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_DIVISION returns the result of integer division.
  Discussion:
    This routine computes C = A / B, where the result is rounded to the
    integer value nearest 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 March 2013
  Author:
    John Burkardt
  Parameters:
    Input, int A, the number to be divided.
    Input, int B, the divisor.
    Output, int I4_DIVISION, the result.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int a = a_data[0];
	const register int b = a_data[1];
	
    int a_abs;
    int b_abs;
    int s;
    int value;
    s = 1 - (((a*b)<0)<<1);
    a_abs = abs ( a );
    b_abs = abs ( b );
    value = s * ( a_abs / b_abs );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_divp ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_DIVP returns the smallest multiple of J greater than or equal to an I4.
  Example:
    I  J  I4_DIVP(I,J)
    0  4    0
    1  4    1
    2  4    1
    3  4    1
    4  4    1
    5  4    2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I, the number to be analyzed.
    Input, int J, the number, multiples of which will
    be compared against I.  J may not be zero.
    Output, int I4_DIVP, the smallest multiple of J that
    is greater than or equal to I.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int i = a_data[0];
	const register int j = a_data[1];
	
	result = j == 0 ? INT_MAX : 1 + ( i - 1 ) / j;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_fall ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_FALL computes the falling factorial function [X]_N.
  Discussion:
    Note that the number of "injections" or 1-to-1 mappings from
    a set of N elements to a set of M elements is [M]_N.
    The number of permutations of N objects out of M is [M]_N.
    Moreover, the Stirling numbers of the first kind can be used
    to convert a falling factorial into a polynomial, as follows:
      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
  Formula:
    [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int X, the argument of the falling factorial function.
    Input, int N, the order of the falling factorial function.
    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
    negative, a "rising" factorial will be computed.
    Output, int I4_FALL, the value of the falling factorial function.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	register int x = a_data[0];
	const register int n = a_data[1];
	
    int i;
    int value;

    value = 1;

    if ( 0 < n )
    {
        for ( i = 1; i <= n; ++i )
        {
            value *= x;
            -- x;
        }
    }
    else if ( n < 0 )
    {
        for ( i = -1; n <= i; --i )
        {
            value *= x;
            -- x;
        }
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_floor ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_FLOOR rounds an R8 down to the nearest I4.
  Example:
    X        I4_FLOOR(X)
   -1.1      -2
   -1.0      -1
   -0.9      -1
   -0.1      -1
    0.0       0
    0.1       0
    0.9       0
    1.0       1
    1.1       1
    2.9       2
    3.0       3
    3.14159   3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, double X, the number whose floor is desired.
    Output, int I4_FLOOR, the floor of X.
*/
{
	static int result = INT_MAX;
	
	const register ityp x = *(ityp *) data; 
	
    int value = ( int ) x;
    if ( x < value )
        -- value;
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4_fraction ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_FRACTION computes a ratio and returns an integer result.
  Discussion:
    Given integer variables I and J, FORTRAN will evaluate the expression
    "I/J" using integer arithmetic.  This routine, which carries out the
    same operation, is thus not needed in FORTRAN.  It is provided simply
    to match the corresponding function in MATLAB, where the default
    result of "I/J" is a real number.
  Example:
       I     J     Real     K = I4_FRACTION ( I, J)
       1     2     0.5      0
       8     4     2.00     2
       9     4     2.25     2
       7     4     1.75     1
      -7     4    -1.75    -1
       7    -4    -1.75    -1
      -7    -4     1.75     1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the arguments.
    Output, int K, the value of the ratio.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int i = a_data[0];
	const register int j = a_data[1];
	
	result = i / j;
    return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_gcdb ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_GCDB finds the greatest common divisor of the form K^N of two I4's.
  Discussion:
    Note that if J is negative, I4_GCDB will also be negative.
    This is because it is likely that the caller is forming
    the fraction I/J, and so any minus sign should be
    factored out of J.
    If I and J are both zero, I4_GCDB is returned as 1.
    If I is zero and J is not, I4_GCDB is returned as J,
    and vice versa.
    If I and J are nonzero, and have no common divisor of the
    form K^N, I4_GCDB is returned as 1.
    Otherwise, I4_GCDB is returned as the largest common divisor
    of the form K^N shared by I and J.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, two numbers whose greatest common divisor K^N
    is desired.
    Input, int K, the possible divisor of I and J.
    Output, int I4_GCDB, the greatest common divisor of
    the form K^N shared by I and J.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	register int i = a_data[0];
	register int j = a_data[1];
	register int k = a_data[2];
	
    int value = 1;
    /*
    If both I and J are zero, I4_GCDB is 1.
    */
    if ( i == 0 && j == 0 )
    {
    	result = 1;
        return &result;
    }
        /*
    If just one of I and J is zero, I4_GCDB is the other one.
    */
    if ( i == 0 )
    {
    	result = j;
        return &result;
    }
    else if ( j == 0 )
    {
    	result = i;
        return &result;
    }
	/*
    Divide out K as long as you can.
    */

    value = -1 + ((0<j)<<1);

    for ( ; ; )
    {
        if ( ( i % k ) != 0 || ( j % k ) != 0 )
            break;

        value *= k;
        i /= k;
        j /= k;
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_huge_normalizer ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_HUGE_NORMALIZER returns the "normalizer" for I4_HUGE.
  Discussion:
    The value returned is 1 / ( I4_HUGE + 1 ).
    For any I4, it should be the case that
     -1 < I4 * I4_HUGE_NORMALIZER < 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 January 2007
  Author:
    John Burkardt
  Parameters:
    Output, double I4_HUGE_NORMALIZER, the "normalizer"
    for I4_HUGE.
*/
{
	static ityp result = 4.656612873077392578125E-10;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_is_even ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_IS_EVEN returns TRUE if an I4 is even.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I, the integer to be tested.
    Output, int I4_IS_EVEN, is TRUE if I is even.
*/
{
	static bool result = 2;
	
	const register int i = *(int *) data; 
	
	result = ( ( i % 2 ) == 0 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_is_odd ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_IS_ODD returns TRUE if an I4 is odd.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I, the integer to be tested.
    Output, int I4_IS_ODD, is TRUE if I is odd.
*/
{
	static bool result = 2;
	
	const register int i = *(int *) data; 
	
	result = ( ( i % 2 ) != 0 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_is_power_of_2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_IS_POWER_OF_2 reports whether an I4 is a power of 2.
  Discussion:
    The powers of 2 are 1, 2, 4, 8, 16, and so on.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int N, the integer to be tested.
    Output, int I4_IS_POWER_OF_2, is TRUE if N is a power of 2.
*/
{
	static bool result = 2;
	
	register int n = *(int *) data; 
	
    int value;

    if ( n <= 0 )
    {
    	result = false;
        return &result;
    }

    while ( n != 1 )
    {
        if ( ( n % 2 ) == 1 )
        {
        	result = false;
        	return &result;
        }
        n >>= 1;
    }

    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_is_prime ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_IS_PRIME reports whether an I4 is prime.
  Discussion:
    A simple, unoptimized sieve of Erasthosthenes is used to
    check whether N can be divided by any integer between 2
    and SQRT(N).
    Note that negative numbers, 0 and 1 are not considered prime.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int N, the integer to be tested.
    Output, int I4_IS_PRIME, is TRUE (1) if N is prime, and FALSE (0)
    otherwise.
*/
{
	static bool result = 2;
	
	const register int n = *(int *) data; 
	
    int i;
    int nhi;

    if ( n <= 0 || n == 1 )
    {
    	result = false;
        return &result;
    }

    if ( n <= 3 )
    {
    	result = true; 
        return &result;
    }

    nhi = ( int ) ( sqrt ( ( ityp ) ( n ) ) );

    for ( i = 2; i <= nhi; i++ )
    {
        if ( ( n % i ) == 0 )
    	{
    		result = false; 
        	return &result;
    	}
    }


    result = true; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_lcm ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_LCM computes the least common multiple of two I4's.
  Discussion:
    The least common multiple may be defined as
      LCM(I,J) = ABS( I * J ) / GCF(I,J)
    where GCF(I,J) is the greatest common factor of I and J.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input,
    Output, int I4_LCM, the least common multiple of I and J.
    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int i = a_data[0];
	const register int j = a_data[1];
	
	result = abs ( i * ( j / i4_gcd ( i, j ) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_log_10 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
  Example:
        I  I4_LOG_10
    -----  --------
        0    0
        1    0
        2    0
        9    0
       10    1
       11    1
       99    1
      100    2
      101    2
      999    2
     1000    3
     1001    3
     9999    3
    10000    4
  Discussion:
    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I, the number whose logarithm base 10 is desired.
    Output, int I4_LOG_10, the integer part of the logarithm base 10 of
    the absolute value of X.
*/
{
	static int result = INT_MAX;
	
	const register int i = *(int *) data; 
	
    int i_abs;
    int ten_pow;
    int value;

    if ( i == 0 )
        value = 0;
    else
    {
        value = 0;
        ten_pow = 10;

        i_abs = abs ( i );

        while ( ten_pow <= i_abs )
        {
            ++ value;
            ten_pow *= 10;
        }
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_log_2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
  Example:
        I  I4_LOG_10
    -----  --------
        0    0
        1    0
        2    1
        3    1
        4    2
        5    2
        7    2
        8    3
        9    3
     1000    9
     1024   10
  Discussion:
    I4_LOG_2 ( I ) + 1 is the number of binary digits in I.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I, the number whose logarithm base 2 is desired.
    Output, int I4_LOG_2, the integer part of the logarithm base 2 of
    the absolute value of X.
*/
{
	static int result = INT_MAX;
	
	const register int i = *(int *) data; 
	
    int i_abs;
    int two_pow;
    int value;

    if ( i == 0 )
        value = 0;
    else
    {
        value = 0;
        two_pow = 2;

        i_abs = abs ( i );

        while ( two_pow <= i_abs )
        {
            ++ value;
            two_pow <<= 1;
        }
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_log_i4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_LOG_I4 returns the logarithm of an I4 to an I4 base.
  Discussion:
    Only the integer part of the logarithm is returned.
    If
      K4 = I4_LOG_J4 ( I4, J4 ),
    then we ordinarily have
      J4^(K4-1) < I4 <= J4^K4.
    The base J4 should be positive, and at least 2.  If J4 is negative,
    a computation is made using the absolute value of J4.  If J4 is
    -1, 0, or 1, the logarithm is returned as 0.
    The number I4 should be positive and at least 2.  If I4 is negative,
    a computation is made using the absolute value of I4.  If I4 is
    -1, 0, or 1, then the logarithm is returned as 0.
    An I4 is an integer ( kind = 4 ) value.
  Example:
    I4  J4  K4
     0   3   0
     1   3   0
     2   3   0
     3   3   1
     4   3   1
     8   3   1
     9   3   2
    10   3   2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I4, the number whose logarithm is desired.
    Input, int J4, the base of the logarithms.
    Output, int I4_LOG_I4, the integer part of the logarithm
    base abs(J4) of abs(I4).
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int i4 = a_data[0];
	const register int j4 = a_data[1];
	
    int i4_abs = abs ( i4 );
    int j4_abs;
    int value = 0;

    if ( 2 <= i4_abs )
    {
        j4_abs = abs ( j4 );

        if ( 2 <= j4_abs )
            while ( j4_abs <= i4_abs )
            {
                i4_abs /= j4_abs;
                ++ value;
            }
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_log_r8 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_LOG_R8 returns the integer part of the logarithm base B of an I4.
  Discussion:
    The base B is normally positive, but in any case only the absolute value
    of B is considered.
    The number X is normally positive, but in any case only the absolute
    value of X is considered.
  Example:
    If B is greater than 1, and X is positive:
    if 1/B^2  <  X <= 1/B   I4_LOG_R8(X) = -1,
    if 1/B    <  X <= 1     I4_LOG_R8(X) = 0,
    if 1      <= X <  B,    I4_LOG_R8(X) = 0,
    if B      <= X <  B^2   I4_LOG_R8(X) = 1,
    if B^2    <= X <  B^3   I4_LOG_R8(X) = 2.
    For positive I4_LOG_R8(X), it should be true that
      ABS(B)^I4_LOG_R8(X) <= ABS(X) < ABS(B)^(I4_LOG_R8(X)+1).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int X, the number whose logarithm base B is desired.
    If X is 0, then I4_LOG_B is returned as -HUGE().
    Input, double B, the absolute value of the base of the
    logarithms.  B must not be -1, 0, or 1.
    Output, int I4_LOG_R8, the integer part of the logarithm
    base abs(B) of abs(X).
*/
{
	static int result = INT_MAX;
	
	const iti * const s_data = data;
	
	const register ityp b = s_data->a0;
	const register int x = s_data->a1;
	
    ityp b_abs;
    int value;
    int value_sign;
    ityp x_abs;

    if ( x == 0 )
    {
    	result = - i4_huge;
        return &result;
    }

    b_abs = fabs ( b );
    value = 0;

    if ( b_abs == 1.0 || b == 0.00 )
    {
    	result = value;
        return &result;
    }


    x_abs = fabs ( ( ityp ) ( x ) );

    if ( b_abs < 1.00 )
    {
        value_sign = -1;
        b_abs = 1.00 / b_abs;
    }
    else
        value_sign = +1;


    if ( 1.00 <= x_abs && x_abs < b_abs )
    {
    	result = value_sign * value;
        return &result;
    }

    while ( b_abs < x_abs )
    {
        x_abs /= b_abs;
        ++ value;
    }

    while ( x_abs * b_abs <= 1.00 )
    {
        x_abs *= b_abs;
        -- value;
    }
    /*
    If the absolute value of the base was less than 1, we inverted
    earlier.  Now negate the logarithm to account for that.
    */
    value *= value_sign;
    
	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_mant ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_MANT computes the "mantissa" of an r8.
  Discussion:
    I4_MANT computes the "mantissa" or "fraction part" of a real
    number X, which it stores as a pair of integers, (J/K).
    It also computes the sign, and the integer part of the logarithm
 (base 2) of X.
    On return:
      X = S * (J/K) * 2^L
    where
      S is +1 or -1,
      K is a power of 2,
      1 <= (J/K) < 2,
      L is an integer.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, double X, the real number to be decomposed.
    Output, int *S, the "sign" of the number.
    S will be -1 if X is less than 0, and +1 if X is greater
    than or equal to zero.
    Output, int J, the top part of the mantissa fraction.
    Output, int K, the bottom part of the mantissa
    fraction.  K is a power of 2.
    Output, int L, the integer part of the logarithm (base 2) of X.
*/
{
	const it4pi * const s_data = data;
	register ityp x = s_data->a0;
	int * s = s_data->a1;
	int * j = s_data->a2;
	int * k = s_data->a3;
	int * l = s_data->a4;
	
    /*
    1: Handle the special case of 0.
    */
    if ( x == 0.00 )
    {
        *s = *k = 1;
        *j = *l = 0;
        return NULL;
    }
    /*
    2: Determine the sign IS.
    */
    if ( 0.00 < x )
        *s = 1;
    else
    {
        *s = -1;
        x = -x;
    }
    /*
    3: Force X to lie between 1 and 2, and compute the logarithm L.
    */
    *l = 0;

    while ( 2.00 <= x )
    {
        x /=  2.00;
        ++ *l;
    }

    while ( x < 1.00 )
    {
        x *= 2.00;
        -- *l;
    }
    /*
    4: Now strip out the mantissa as J/K.
    */
    *j = 0;
    *k = 1;

    for ( ; ; )
    {
        *j = (*j ) << 1;

        if ( 1.0 <= x )
        {
            ++ *j;
            -- x;
        }

        if ( x == 0.0 )
            break;

        *k = ( *k ) << 1;
        x *= 2.00;
    }/*
    1: Handle the special case of 0.
    */
    if ( x == 0.00 )
    {
        *s = *k = 1;
        *j = *l = 0;
        return NULL;
    }
    /*
    2: Determine the sign IS.
    */
    if ( 0.00 < x )
        *s = 1;
    else
    {
        *s = -1;
        x *= -1;
    }
    /*
    3: Force X to lie between 1 and 2, and compute the logarithm L.
    */
    *l = 0;

    while ( 2.00 <= x )
    {
        x /= 2.00;
        ++ *l;
    }

    while ( x < 1.00 )
    {
        x *= 2.00;
        -- *l;
    }
    /*
    4: Now strip out the mantissa as J/K.
    */
    *j = 0;
    *k = 1;

    for ( ; ; )
    {
        *j = ( *j ) << 1;

        for ( ; ; )
        {
            *j = ( *j ) << 1;

	        if ( 1.0 <= x )
	        {
	            ++ *j;
	            -- x;
	        }
	
	        if ( x == 0.0 )
	            break;
	
	        *k = ( *k ) << 1;
	         x *= 2.00;
		}
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_mod_inv ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_MOD_INV calculates the inverse of B mod N.
  Discussion:
    This function uses the extended Euclidean algorithm.
    Unless the algorithm fails, the output value Y will satisfy
  ( B * Y ) mod N = 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2011
  Author:
    Original MATLAB version by Wade Trappe, Lawrence Washington.
    C version by John Burkardt.
  Reference:
    Wade Trappe, Lawrence Washington,
    Introduction to Cryptography with Coding Theory,
    Prentice Hall, 2005,
    ISBN13: 978-0131862395,
    LC: QA268.T73.
  Parameters:
    Input, int B, the value whose inverse is desired.
    B must not be 0, or a multiple of N.  However, B can be negative.
    Input, int N, the value with respect to which the inverse
    is desired.  N must be 2 or greater.
    Output, int I4_MOD_INV, the inverse of B mod N.  However, if the
    inverse does not exist, Y is returned as 0.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int b = a_data[0];
	const register int n = a_data[1];
	
    int b0;
    int n0;
    int q;
    int r;
    int t;
    int t0;
    int temp;
    int y;

    n0 = n;
    b0 = abs ( b );
    t0 = 0;
    t = 1;

    q = ( n0 / b0 );
    r = n0 - q * b0;

    while ( 0 < r )
    {
        temp = t0 - q * t;
        temp = 0<=temp ? temp%n : n - ( ( - temp ) % n );
        n0 = b0;
        b0 = r;
        t0 = t;
        t = temp;
        q = ( n0 / b0 );
        r = n0 - q * b0;
    }

    if ( b0 != 1 )
        y = 0;
    else
    {
        y = t % n;
        if ( b < 0 )
            y *= -1;
    }

	result = y;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_moddiv ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_MODDIV breaks an I4 into a multiple of a divisor and remainder.
  Discussion:
    N = M * D + R
    0 <= || R || < || D ||
    R has the sign of N.
  Example:
    N         D       M      R
   107       50      2      7
   107      -50     -2      7
  -107       50     -2     -7
  -107      -50      2     -7

  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number to be decomposed.
    Input, int D, the divisor.  D may not be zero.
    Output, int *M, the number of times N
    is evenly divided by D.
    Output, int *R, a remainder.
*/
{
	const _2i2pi * const s_data = data;
	const register int n = s_data->a0;
	const register int d = s_data->a1;
	int * m = s_data->a2;
	int * r = s_data->a3;
	
    if ( d == 0 )
        return NULL;
        
    *m = n / d;
    *r = n - d * ( *m );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_reverse_bytes ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_REVERSE_BYTES reverses the bytes in an I4.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int X, a value whose bytes are to be reversed.
    Output, int I4_REVERSE_BYTES, a value whose bytes are
    in reverse order from those in X.
*/
{
	static int result = INT_MAX;
	
	const register int x = *(int *) data;
	
    char c;
    union
    {
        int yint;
        char ychar[4];
    } y;

    y.yint = x;

    c = y.ychar[0];
    y.ychar[0] = y.ychar[3];
    y.ychar[3] = c;

    c = y.ychar[1];
    y.ychar[1] = y.ychar[2];
    y.ychar[2] = c;

	result = ( y.yint );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_rise ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_RISE computes the rising factorial function [X]^N.
  Discussion:
    [X}^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
    Note that the number of ways of arranging N objects in M ordered
    boxes is [M}^N. (Here, the ordering in each box matters).  Thus,
    2 objects in 2 boxes have the following 6 possible arrangements:
      -/12, 1/2, 12/-, -/21, 2/1, 21/-.
    Moreover, the number of non-decreasing maps from a set of
    N to a set of M ordered elements is [M]^N / N!.  Thus the set of
    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int X, the argument of the rising factorial function.
    Input, int N, the order of the rising factorial function.
    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
    negative, a "falling" factorial will be computed.
    Output, int I4_RISE, the value of the rising factorial function.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	register int x = a_data[0];
	const register int n = a_data[1];
	
    int i;
    int value = 1;

    if ( 0 < n )
    {
        for ( i = 1; i <= n; ++i )
        {
            value *= x;
            ++ x;
        }
    }
    else if ( n < 0 )
    {
        for ( i = -1; n <= i; --i)
        {
            value *= x;
            -- x;
        }
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_sign3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_SIGN3 returns the three-way sign of an I4.
  Discussion:
    This is the "three-way" sign, that is:

      I < 0        =>  sign(I) = -1
      I = 0        =>  sign(I) =  0
          0 <  I   =>  sign(I) = +1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int I, the integer whose three-way sign is desired.
    Output, int I4_SIGN3, the sign of I.
*/
{
	static char result = CHAR_MAX;
	
	const register int i = *(int *) data;
	
    if ( i < 0 )
    {
    	result = -1;
        return &result;
    }
    else if ( i == 0 )
    {
    	result = 0;
        return &result;
    }
        
    result = 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4_swap3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_SWAP3 swaps three I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input/output, int *I, *J, *K.  On output, the values of I, J, and K
    have been interchanged.
*/
{
	int ** const a_data = data;
	int * i = a_data[0]; 
	int * j = a_data[1];
	int * k = a_data[2];
	
    int l;
    l = *i;
    *i = *j;
    *j = *k;
    *k =  l;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_to_angle ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_ANGLE maps I4's to points on a circle.
  Discussion:
    The angles are intended to be used to select colors on a color
    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
    magenta.
  Example:
     I   X      ANGLE
     0   0/3      0
     1   1/3    120
     2   2/3    240
     3   1/6     60
     4   3/6    180
     5   5/6    300
     6   1/12    30
     7   3/12    90
     8   5/12   150
     9   7/12   210
    10   9/12   270
    11  11/12   330
    12   1/24    15
    13   3/24    45
    14   5/24    75
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2007
  Author:
    John Burkardt
  Parameters:
    Input, int I, the index of the desired color.
    Output, double I4_TO_ANGLE, an angle, measured in degrees,
    between 0 and 360.
*/
{
	static ityp result = MAX_VAL;
	
	const register int i = *(int *) data;
	
    int i1;
    int i2;
    int i3;
    int i4;
    int value;

    if ( 0 <= abs ( i ) && abs ( i ) <= 2 )
        value = 120.00 * ( ityp ) ( abs ( i ) );
    else
    {
        i1 = i4_log_2 ( abs ( i ) / 3 );
        i2 = abs ( i ) + 1 - 3 * powi ( 2, i1 );
        i3 = (( i2 - 1 ) <<1) + 1;
        i4 = 3 * powi ( 2, ( i1 + 1 ) );
        value = 360.0 * ( ityp ) ( i3 ) / ( ityp ) ( i4 );
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4_to_digits_binary ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_DIGITS_BINARY produces the binary digits of an I4.
  Example:
     I    N     C               Binary
    --  ---   ---         ------------
     0    1   0                      0
     0    2   0, 0                  00
     1    3   1, 0, 0              100
     2    3   0, 1, 0              010
     3    3   1, 1, 0              011
     4    3   0, 0, 1              100
     8    3   0, 0, 0        (1)000
     8    5   0, 0, 0, 1, 0      01000
    -8    5   0, 0, 0, 1, 0 (-) 0100
     0    3   0, 0, 0
     1    3   1, 0, 0
     2    3   0, 1, 0
     3    3   1, 1, 0
     4    3   0, 0, 1
     5    3   1, 0, 1
     6    3   0, 1, 1
     7    3   1, 1, 1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 December 2011
  Author:
    John Burkardt
  Parameters:
    Input, int I, the integer to be analyzed.
    Input, int N, the number of digits to determine.
    Output, int I4_TO_DIGITS_BINARY[N], the first N binary digits of I.
    Entry 0 is the units digit.
*/
{
	int * const a_data = data;
	register int i = a_data[0];
	const register int n = a_data[1];
	
    int *c;
    int j;

    c = ( int * ) malloc ( n * sizeof ( int ) );
    i = abs ( i );

    for ( j = 0; j < n; ++j )
    {
        c[j] = i % 2;
        i = ( i - c[j] ) / 2;
    }

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4_to_digits_decimal ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_DIGITS_DECIMAL determines the last N decimal digits of an I4.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2008
  Author:
    John Burkardt
  Parameters:
    Input, int I, the integer to be analyzed.
    Input, int N, the number of digits to determine.
    Output, int I4_TO_DIGITS_DECIMAL[N], the last N decimal digits of I.
    DIGIT[I-1] is the "coefficient" of 10^(I-1).
*/
{
	int * const a_data = data;
	register int i = a_data[0];
	const register int n = a_data[1];
	
    int *digit;
    int j;

    digit = ( int * ) malloc ( n * sizeof ( int ) );

    i = abs ( i );

    for ( j = 1; j <= n; ++j )
    {
        digit[j-1] = i % 10;
        i = ( i - digit[j-1] ) / 10;
    }

    return digit;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4_to_fac ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_FAC converts an I4 into a product of prime factors.
  Discussion:
    This routine will fail if the input integer is not positive,
    or if PRIME_NUM is too small to account for the factors of the integer.
    The formula is:
      I = Product ( 1 <= J <= PRIME_NUM ) PRIME(J)^NPOWER(J).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int I, the integer to be factored.
    Input, int PRIME_NUM, the number of prime factors for
    which storage has been allocated.
    Output, int I4_TO_FAC[PRIME_NUM], the powers of the primes.
*/
{
	int * const a_data = data;
	register int i = a_data[0];
	const register int prime_num = a_data[1];
	
    int j;
    int *npower;
    int p;

    if ( i <= 0 )
        return NULL;

    npower = ( int * ) malloc ( prime_num * sizeof ( int ) );
    /*
    Try dividing the remainder by each prime.
    */
    for ( j = 1; j <= prime_num; ++j )
    {
        npower[j-1] = 0;

        p = prime ( j );

        while ( ( i % p ) == 0 )
        {
            ++ npower[j-1];
            i /= p;
        }
    }
    return npower;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4_to_halton ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_HALTON computes one element of a leaped Halton subsequence.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Reference:
    John Halton,
    On the efficiency of certain quasi-random sequences of points
    in evaluating multi-dimensional integrals,
    Numerische Mathematik,
    Volume 2, 1960, pages 84-90.
    John Halton, GB Smith,
    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
    Communications of the ACM,
    Volume 7, 1964, pages 701-702.
    Ladislav Kocis, William Whiten,
    Computational Investigations of Low-Discrepancy Sequences,
    ACM Transactions on Mathematical Software,
    Volume 23, Number 2, 1997, pages 266-294.
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    1 <= DIM_NUM is required.
    Input, int STEP, the index of the subsequence element.
    0 <= STEP is required.
    Input, int SEED[DIM_NUM], the Halton sequence index corresponding
    to STEP = 0.
    0 <= SEED(1:DIM_NUM) is required.
    Input, int LEAP[DIM_NUM], the successive jumps in the Halton sequence.
    1 <= LEAP(1:DIM_NUM) is required.
    Input, int BASE[DIM_NUM], the Halton bases.
    1 < BASE(1:DIM_NUM) is required.
    Output, double R[DIM_NUM], the STEP-th element of the leaped
    Halton subsequence.
*/
{
	const _2dtpi2pdtpit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	const register dim_typ step = s_data->a1;
	int * seed = s_data->a2;
	dim_typ * leap = s_data->a3;
	dim_typ * base = s_data->a4;
	ityp * r = s_data->a5;
	
    ityp base_inv;
    ityp digit;
    dim_typ i;
    int seed2;
    /*
    Calculate the data.
    */
    for ( i = 0; i < dim_num; ++i )
    {
        seed2 += step * leap[i];
        r[i] = 0.0;
        base_inv = 1.00 / ( ( ityp ) base[i] );

        while ( seed2 != 0 )
        {
            digit = seed2 % base[i];
            r[i] += ( ( ityp ) digit ) * base_inv;
            base_inv /= ( ( ityp ) base[i] );
            seed2 /= base[i];
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_to_isbn ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_ISBN converts an I4 to an ISBN digit.
  Discussion:
    Only the integers 0 through 10 can be input.  The representation
    of 10 is 'X'.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Reference:
    Book Industry Study Group,
    The Evolution in Product Identification:
    Sunrise 2005 and the ISBN-13,
    http:www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
  Parameters:
    Input, int I, an integer between 0 and 10.
    Output, char I4_TO_ISBN, the ISBN character code of the integer.
    If I is illegal, then I4_TO_ISBN is set to '?'.
*/
{
	static char result = CHAR_MAX;
	
	const register int i = *(int *) data;
	char value = '?';
	
    switch(i)
    {
        case 0:
            value = '0';
            break;
        case 1:
            value = '1';
            break;
        case 2:
            value = '2';
            break;
        case 3:
            value = '3';
            break;
        case 4:
            value = 4;
            break;
        case 5:
            value = '5';
            break;
        case 6:
            value = '6';
            break;
        case 7:
            value = '7';
            break;
        case 8:
            value = '8';
            break;
        case 9:
            value = '9';
            break;
        case 10:
            value = 'X';
            break;
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_to_l4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_L4 converts an I4 to an L4.
  Discussion:
    0 is FALSE, and anything else if TRUE.
    An I4 is an integer value.
    An L4 is a logical value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 January 2012
  Author:
    John Burkardt
  Parameters:
    Input, int I4, an integer.
    Output, int I4_TO_L4, the logical value of I4.
*/
{
	static bool result = 2;
	
	const register int i4 = *(int *) data;
	
	result = ( i4 != 0 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_to_pascal ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_PASCAL converts a linear index to Pascal triangle coordinates.
  Discussion:
    We describe the grid points in Pascal's triangle in two ways:
    As a linear index K:
                     1
                   2   3
                 4   5   6
               7   8   9   10
    As elements (I,J) of Pascal's triangle:
                     0,0
                  1,0   0,1
               2,0   1,1    0,2
            3,0   2,1   1,2    0,3
  Example:
     K  I  J
     1  0  0
     2  1  0
     3  0  1
     4  2  0
     5  1  1
     6  0  2
     7  3  0
     8  2  1
     9  1  2
    10  0  3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int K, the linear index of the (I,J) element.
    1 <= K.
    Output, int *I, *J, the Pascal indices.
*/
{
	const i2pi * const s_data = data;
	const register int k = s_data->a0;
	int * i = s_data->a1;
	int * j = s_data->a2;
	
    int d;
    if ( k <= 0 )
        return NULL;
    d = i4_to_pascal_degree ( k );
    *j = k - ( d * ( d + 1 ) ) / 2 - 1;
    *i = d - *j;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_to_pascal_degree ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_PASCAL_DEGREE converts a linear index to a Pascal triangle degree.
  Discussion:
    We describe the grid points in Pascal's triangle in two ways:
    As a linear index K:
                     1
                   2   3
                 4   5   6
               7   8   9   10
    As elements (I,J) of Pascal's triangle:
                     0,0
                  1,0   0,1
               2,0   1,1    0,2
            3,0   2,1   1,2    0,3
    The quantity D represents the "degree" of the corresponding monomial,
    that is, D = I + J.
    We can compute D directly from K using the quadratic formula.
  Example:
     K  I  J  D
     1  0  0  0
     2  1  0  1
     3  0  1  1
     4  2  0  2
     5  1  1  2
     6  0  2  2
     7  3  0  3
     8  2  1  3
     9  1  2  3
    10  0  3  3
    11  4  0  4
    12  3  1  4
    13  2  2  4
    14  1  3  4
    15  0  4  4
    16  5  0  5
    17  4  1  5
    18  3  2  5
    19  2  3  5
    20  1  4  5
    21  0  5  5
    22  6  0  6
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int K, the linear index of the (I,J) element.
    1 <= K.
    Output, int I4_TO_PASCAL_DEGREE, the degree (sum) of the corresponding
     Pascal indices.
*/
{
	static int result = INT_MAX;
	
	const register int k = *(int *) data;
	
	result = k <= 0 ? INT_MAX : ( int ) ( 0.50 * ( -1.00 + sqrt ( ( ityp ) ( 1 + (( k - 1 )<<3) ) ) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4_to_triangle ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_TRIANGLE converts an integer to triangular coordinates.
  Discussion:
    Triangular coordinates are handy when storing a naturally triangular
    array (such as the lower half of a matrix) in a linear array.
    Thus, for example, we might consider storing
 (0,0)
 (1,0) (1,1)
 (2,0) (2,1) (2,2)
 (3,0) (3,1) (3,2) (3,3)
    as the linear array
 (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) (3,0) (3,1) (3,2) (3,3)
    Here, the quantities in parenthesis represent the natural row and
    column indices of a single number when stored in a rectangular array.
    In this routine, we are given the location K of an item in the
    linear array, and wish to determine the row I and column J
    of the item when stored in the triangular array.
  Example:
    K  I  J
    0  0  0
    1  1  0
    2  1  1
    3  2  0
    4  2  1
    5  2  2
    6  3  0
    7  3  1
    8  3  2
    9  3  3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 January 2009
  Author:
    John Burkardt
  Parameters:
    Input, int K, the linear index of the (I,J) element, which
    must be nonnegative.
    Output, int *I, *J, the row and column indices.
*/
{
	const i2pi * const s_data = data;
	const register int k = s_data->a0;
	int * i = s_data->a1;
	int * j = s_data->a2;
	
    int c, r;

    if ( k < 0 )
        return NULL;
    else if ( k == 0 )
    {
        *i = *j = 0;
        return NULL;
    }
    /*
 ( N - 1 )^2 + ( N - 1 ) < 2 * K <= N^2 + N
    */
    r = ( int ) ( sqrt ( ( ityp ) ( ( k + 1 ) << 1 ) ) );

    if ( r * r + r <( k + 1 ) << 1 )
        ++ r;

    -- r;

    c = k - ( r * ( r + 1 ) ) / 2;

    *i = r;
    *j = c;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4_unswap3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_UNSWAP3 unswaps three I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 April 2004
  Author:
    John Burkardt
  Parameters:
    Input/output, int *I, *J, *K.  On output, the values of I, J, and K
    have been interchanged.
*/
{
	int ** const a_data = data;
	int * i = a_data[0];
	int * j = a_data[1];
	int * k = a_data[2];
	
    int l;
    l = *k;
    *k = *j;
    *j = *i;
    *i =  l;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_walsh_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_WALSH_1D evaluates the Walsh function of a real scalar argument.
  Discussion:
    Consider the binary representation of X, and number the digits
    in descending order, from leading to lowest, with the units digit
    being numbered 0.
    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the Walsh function.
    Input, int DIGIT, the index of the Walsh function.
    Output, int I4_WALSH_1D, the value of the Walsh function.
*/
{
	static int result = INT_MAX;
	
	const iti * const s_data = data;
	register ityp x = s_data->a0;
	const register int digit = s_data->a1;
	
    int n;
    int value;
    /*
    Hide the effect of the sign of X.
    */
    x = abs ( x );
    /*
    If DIGIT is positive, divide by 2 DIGIT times.
    If DIGIT is negative, multiply by 2 (-DIGIT) times.
    */
    x = 0 <= digit ? x / ( ityp ) powi ( 2, digit) : x * ( ityp ) powi ( 2, -digit );
    /*
    Make it an integer.
    Because it's positive, and we're using INT, we don't change the
    units digit.
    */
    n = ( int ) x ;
    /*
    Is the units digit odd or even?
    */

	result = 0 + ((n%2)!=0);
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_width ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_WIDTH returns the "width" of an I4.
  Example:
        I  I4_WIDTH
    -----  -------
    -1234    5
     -123    4
      -12    3
       -1    2
        0    1
        1    1
       12    2
      123    3
     1234    4
    12345    5
  Discussion:
    The width of an integer is the number of characters necessary to print it.
    The width of an integer can be useful when setting the appropriate output
    format for a vector or array of values.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int I, the number whose width is desired.
    Output, int I4_WIDTH, the number of characters necessary to represent
    the integer in base 10, including a negative sign if necessary.
*/
{
	static int result = INT_MAX;
	
	const register int i = *(int *) data;
	
	result = i4_log_10(i) + 1 + 0>i;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_wrap ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_WRAP forces an I4 to lie between given limits by wrapping.
  Example:
    ILO = 4, IHI = 8
    I
    -2     8
    -1     4
     0     5
     1     6
     2     7
     3     8
     4     4
     5     5
     6     6
     7     7
     8     8
     9     4
    10     5
    11     6
    12     7
    13     8
    14     4
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 December 2012
  Author:
    John Burkardt
  Parameters:
    Input, int IVAL, an integer value.
    Input, int ILO, IHI, the desired bounds for the integer value.
    Output, int I4_WRAP, a "wrapped" version of IVAL.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int ival = a_data[0];
	const register  int ilo = a_data[1];
	const register  int ihi = a_data[2];
	
    int jhi;
    int jlo;
    int value;
    int wide;

    if ( ilo < ihi )
    {
        jlo = ilo;
        jhi = ihi;
    }
    else
    {
        jlo = ihi;
        jhi = ilo;
    }

    wide = jhi + 1 - jlo;
    
    result = jlo + (wide!=1)*i4_modp ( ival - jlo, wide );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4_xor ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_XOR calculates the exclusive OR of two I4's.
  Discussion:
    C provides the operator "^" which produces the same result, faster.
    This code is simply provided for illustration.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 January 2010
  Author:
   John Burkardt
  Parameters:
    Input, unsigned int I, J, two values whose exclusive OR is needed.
    Output, unsigned int I4_XOR, the exclusive OR of I and J.
*/
{
	static unsigned int result = UINT_MAX;
	
	const unsigned int * const a_data = data;
	register unsigned int i = a_data[0];
	register unsigned int j = a_data[1];
	
    unsigned int i2;
    unsigned int j2;
    unsigned int k;
    unsigned int l;

    k = 0;
    l = 1;

    while ( i != 0 || j != 0 )
    {
        i2 = i / 2;
        j2 = j / 2;

        if (( ( i == (i2<<1) ) && ( j != (j2<<1) ) ) ||( ( i != (i2<<1) ) && ( j == (j2<<1) ) ) )
            k += l;

        i = i2;
        j = j2;
        l = l<<1;
    }

	result = k;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i43mat_flip_cols ( void * data)
/******************************************************************************/
/*
  Purpose:
    I43MAT_FLIP_COLS swaps the columns of an I43MAT.
  Discussion:
    An I43MAT is a matrix, each of whose entries is an I43, a triple
    of integers.
    An I43MAT can be stored as a 3 x M x N array, where M counts the "columns"
    and N counts the "rows".
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, int A[3*M*N], the matrix whose columns are to be flipped.
*/
{
	const _2dtpdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	dim_typ * a = s_data->a2;
	
    dim_typ b, i, j, k;

    for ( k = 0; k < ( n / 2 ); ++k )
        for ( j = 0; j < m; ++j )
            #pragma omp parallel for num_threads(3)
            for ( i = 0; i < 3; ++i )
            {
                b                    = a[i+j*3+     k *m*3];
                a[i+j*3+     k *m*3] = a[i+j*3+(n-1-k)*m*3];
                a[i+j*3+(n-1-k)*m*3] = b;
            }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i43mat_flip_rows ( void * data)
/******************************************************************************/
/*
  Purpose:
    I43MAT_FLIP_ROWS swaps the rows of an I43MAT.
  Discussion:
    An I43MAT is a matrix, each of whose entries is an I43, a triple of integers.
    An I43MAT can be stored as a 3 x M x N array, where M counts the "columns"
    and N counts the "rows".
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
   21 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, int A[3*M*N], the matrix whose rows are to be flipped.
*/
{
	const _2dtpdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	dim_typ * a = s_data->a2;
	
    dim_typ b, i, j, k;

    for ( k = 0; k < n; ++k )
        for ( j = 0; j < ( m / 2 ); ++j )
            #pragma omp parallel for num_threads(3)
            for ( i = 0; i < 3; ++i )
            {
                b                    = a[i+     j *3+k*m*3];
                a[i+     j *3+k*m*3] = a[i+(m-1-j)*3+k*m*3];
                a[i+(m-1-j)*3+k*m*3] = b;
            }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4block_delete ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4BLOCK_DELETE frees memory associated with an I4BLOCK.
  Discussion:
    This function releases the memory associated with an array that was
    created by a command like
      int ***a;
      a = i4block_new ( l, m, n );
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2012
  Author:
    John Burkardt
  Parameters:
    Input, int ***A, the pointer to the data.
    Input, int L, M, N, the number of rows, columns, and layers in the array.
*/
{
	const pppi3dt * const s_data = data;
	int ***a = s_data->a0;
	const register dim_typ l = s_data->a1;
	const register dim_typ m = s_data->a2;
	const register dim_typ n = s_data->a3;
	
    dim_typ i, j;

    for ( i = 0; i < l; ++i )
        for ( j = 0; j < m; ++j )
            free ( a[i][j] );

    for ( i = 0; i < l; ++i )
        free ( a[i] );

    free ( a );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4block_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4BLOCK_NEW allocates a new I4BLOCK.
  Discussion:
    A declaration of the form

      int ***a;
    is necesary.  Then an assignment of the form:
      a = i4block_new ( l, m, n );
    allows the user to assign entries to the matrix using typical
    3D array notation:
      a[2][3][4] = 17;
      y = a[1][0][3];
    and so on.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 March 2012
  Author:
    John Burkardt
  Parameters:
    Input, int L, M, N, the number of rows, columns and layers.
    Output, int I4BLOCK_NEW[L][M][N], a new block.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ l = a_data[0];
	const register dim_typ m = a_data[1];
	const register dim_typ n = a_data[2];
	
    int ***a;
    dim_typ i, j;

    a = ( int *** ) malloc ( l * sizeof ( int ** ) );

    if ( a == NULL )
        return NULL;

    for ( i = 0; i < l; ++i)
    {
        a[i] = ( int ** ) malloc ( m * sizeof ( int * ) );
        if ( a[i] == NULL )
            return NULL;
    }

    for ( i = 0; i < l; ++i )
        for ( j = 0; j < m; ++j)
        {
            a[i][j] = ( int * ) malloc ( n * sizeof ( int ) );
            if ( a[i][j] == NULL )
                return NULL;
        }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  inline void   * _i4block_zeros_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4BLOCK_ZEROS_NEW returns a new zeroed I4BLOCK.
  Discussion:
    An I4BLOCK is a triple dimensioned array of I4 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 April 2013
  Author:
    John Burkardt
  Parameters:
    Input, int L, M, N, the number of rows, columns, and levels.
    Output, int I4BLOCK_ZEROS_NEW[L*M*N], the new zeroed matrix.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ l = a_data[0];
	const register dim_typ m = a_data[1];
	const register dim_typ n = a_data[2];
	
    dim_typ i, j, k;
    int * a = ( int * ) malloc ( l * m * n * sizeof ( int ) );
    for ( k = 0; k < n; k++ )
        for ( j = 0; j < m; j++ )
            for ( i = 0; i < l; i++ )
            a[i+j*l+k*l*m] = 0;
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4cmat_delete ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4CMAT_DELETE frees the memory set aside by I4CMAT_NEW.
  Discussion:
    This function releases the memory associated with an array that was
    created by a command like
      int **a;
      a = i4cmat_new ( m, n );
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int **A, the array.
    Input, int M, N, the number of rows and columns.
*/
{
	const ppi2dt * const s_data = data;
	int **a = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    for (dim_typ j = 0; j < n; ++j )
        free ( a[j] );
    free ( a );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4cmat_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4CMAT_NEW sets up an I4CMAT of the desired dimensions.
  Discussion:
    An I4CMAT is a column-major dynamically dimensioned 2D integer array.
    A declaration of the form
      int **a;
    is necesary.  Then an assignment of the form:
      a = i4cmat_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17;
      y = a[1][0];
    and so on.
    However, for a column major matrix, the first index is the column, the
    second the row.  Thus a[i][j] references row j, column i!
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Output, int **I4CMAT_NEW, the array.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
    dim_typ j;
    int ** a = ( int ** ) malloc ( m * n * sizeof ( int * ) );
    for ( j = 0; j < n; ++j )
        a[j] = ( int * ) malloc ( m * sizeof ( int ) );
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4col_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_COMPARE compares columns I and J of an I4COL.
  Example:
    Input:
      M = 3, N = 4, I = 2, J = 4
      A = (
        1  2  3  4
        5  6  7  8
        9 10 11 12 )
    Output:
      I4COL_COMPARE = -1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A[M*N], an array of N columns of vectors of length M.
    Input, int I, J, the columns to be compared.
    I and J must be between 1 and N.
    Output, int I4COL_COMPARE, the results of the comparison:
    -1, column I < column J,
     0, column I = column J,
    +1, column J < column I.
*/
{
	static unsigned short result = USHRT_MAX;
	
	const _2dtpi2dt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	const register dim_typ i = s_data->a3;
	const register dim_typ j = s_data->a4;
	
    int k;
    /*
    Check.
    */
    if ( i < 1 || n < i || j < 1 || n < j)
    {
    	result = 2;
        return &result;
    }

    if ( i == j )
    {
    	result = 0;
        return &result;
    }

    k = 1;

    while ( k <= m )
    {
        if ( a[k-1+(i-1)*m] < a[k-1+(j-1)*m] )
        {
        	result = 2;
        	return &result;
        }
        else if ( a[k-1+(j-1)*m] < a[k-1+(i-1)*m] )
    	{
    		result = 1;
        	return &result;
    	}
        ++ k;
    }
	
	result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4col_find ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_FIND seeks a table column equal to an I4COL.
  Example:
    M = 3, N = 4,
    A = (
      1  2  3  4
      5  6  7  8
      9 10 11 12 )
    IVEC = ( 3, 7, 11 )
    ICOL = 3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in
    the table.  M is also the length of IVEC.
    Input, int A[M*N], an array of N columns of vectors
    of length M.
    Input, int IVEC[M], a vector to be matched with the data
    in the array.
    Output, int I4COL_FIND, the index of the first column of the table
    which exactly matches every entry of IVEC, or -1 if no match
    could be found.
*/
{
	static short result = SHRT_MAX;
	
	const _2dt2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	int * ivec = s_data->a3;
	
    dim_typ i;
    short col;
    dim_typ j;

    if ( m <= 0 )
    {
    	result = -1;
        return &result;
    }
    for ( j = 1; j <= n; ++j)
    {
        i = 1;
        while ( ivec[i-1] == a[i-1+(j-1)*m] )
        {
            if ( i == m )
            {
            	result = j;
        		return &result;
            }

            ++ i;
        }
    }

    result = -1;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4col_find_item ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_FIND_ITEM searches an I4COL for a given value.
  Discussion:
    The two dimensional information in TABLE is stored as a one dimensional
    array, by columns.
    The values ROW and COL will be one-based indices.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2005
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns
    in the table.
    Input, int TABLE[M*N], the table to search.
    Input, int ITEM, the value to search for.
    Output, int *ROW, *COL, the row and column indices
    of the first occurrence of the value ITEM.  The search
    is conducted by rows.  If the item is not found, then
    ROW = COL = -1.
*/
{
	const _2dtpii2pdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * table = s_data->a2;
	const register int item = s_data->a3;
	dim_typ * row = s_data->a4;
	dim_typ * col = s_data->a5;
	
    dim_typ i, j;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            if ( table[i+j*m] == item )
            {
                *row = i+1;
                *col = j+1;
                return NULL;
            }

    *row = *col = -1;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4col_find_pair_wrap ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_FIND_PAIR_WRAP wrap searches an I4COL for a pair of items.
  Discussion:
    The items (ITEM1, ITEM2) must occur consecutively.
    However, wrapping is allowed, that is, if ITEM1 occurs
    in the last row, and ITEM2 "follows" it in the first row
    of the same column, a match is declared.
    If the pair of items is not found, then ROW = COL = -1.
    The values ROW and COL will be one-based indices.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns
    in the table.
    Input, int A[M*N], the table to search.
    Input, int ITEM1, ITEM2, the values to search for.
    Output, int *ROW, *COL, the row and column indices
    of the first occurrence of the value ITEM1 followed immediately
    by ITEM2.
*/
{
	const _2dtpi2i2pdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	const register int item1 = s_data->a3;
	const register int item2 = s_data->a4;
	dim_typ * row = s_data->a5;
	dim_typ * col = s_data->a6;
	
    dim_typ i, i2, j;

    for ( j = 1; j <= n; ++j )
        for ( i = 1; i <= m; ++i )
            if ( a[i-1+(j-1)*m] == item1 )
            {
                i2 = i + 1;

                if ( m < i2 )
                    i2 = 1;

                if ( a[i2-1+(j-1)*m] == item2 )
                {
                    *row = i;
                    *col = j;
                    return NULL;
                }
            }

    *row = *col = -1;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4col_first_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_FIRST_INDEX indexes the first occurrence of values in an I4COL.
  Discussion:
    An I4COL is an M by N array of I4 values.
    It is regarded as an array of N columns of length M.
    For element A(1:M,J) of the matrix, FIRST_INDEX(J) is the index in A of
    the first column whose entries are equal to A(1:M,J).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of A.
    The length of an "element" of A, and the number of "elements".
    Input, int A[M*N], the array.
    Output, int I4COL_FIRST_INDEX[N], the first occurrence index.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ diff;
    dim_typ *first_index;
    dim_typ i, j1, j2;
    first_index = ( dim_typ * ) malloc ( n * sizeof ( dim_typ ) );

    for ( j1 = 0; j1 < n; j1++ )
        first_index[j1] = -1;

    for ( j1 = 0; j1 < n; j1++ )
    {
        if ( first_index[j1] == -1 )
        {
            first_index[j1] = j1;

            for ( j2 = j1 + 1; j2 < n; j2++ )
            {
                diff = 0;
                for ( i = 0; i < m; i++ )
                    diff += abs ( a[i+j1*m] - a[i+j2*m] );

                if ( diff == 0 )
                    first_index[j2] = j1;
            }
        }
    }

    return first_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4col_sort_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_SORT_A ascending sorts the columns of an I4COL.
  Discussion:
    In lexicographic order, the statement "X < Y", applied to two
    vectors X and Y of length M, means that there is some index I, with
    1 <= I <= M, with the property that
      X(J) = Y(J) for J < I,
    and
      X(I) < Y(I).
    In other words, X is less than Y if, at the first index where they
    differ, the X value is less than the Y value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of A.
    Input, int N, the number of columns of A.
    Input/output, int A[M*N].
    On input, the array of N columns of M vectors;
    On output, the columns of A have been sorted in ascending
    lexicographic order.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    /*
    Initialize.
    */
    i = indx = isgn = j = 0;
    /*
    Call the external heap sorter.
    */
    for ( ; ; )
    {
        sort_heap_external ( n, &indx, &i, &j, isgn );
        /*
        Interchange the I and J objects.
        */
        if ( 0 < indx )
            i4col_swap ( m, n, a, i, j );
        /*
        Compare the I and J objects.
        */

        else if ( indx < 0 )
            isgn = i4col_compare ( m, n, a, i, j );
        else if ( indx == 0 )
            break;

    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4col_sort_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_SORT_D descending sorts the columns of an I4COL.
  Discussion:
    In lexicographic order, the statement "X < Y", applied to two
    vectors X and Y of length M, means that there is some index I, with
    1 <= I <= M, with the property that
      X(J) = Y(J) for J < I,
    and
      X(I) < Y(I).
    In other words, X is less than Y if, at the first index where they
    differ, the X value is less than the Y value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of A.
    Input, int N, the number of columns of A.
    Input/output, int A[M*N].
    On input, the array of N columns of M vectors;
    On output, the columns of A have been sorted in descending
    lexicographic order.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    /*
    Initialize.
    */
    i = indx = isgn = j = 0;
    /*
    Call the external heap sorter.
    */
    for ( ; ; )
    {
        sort_heap_external ( n, &indx, &i, &j, isgn );
        /*
        Interchange the I and J objects.
        */
        if ( 0 < indx )
            i4col_swap ( m, n, a, i, j );
        /*
        Compare the I and J objects.
        */
        else if ( indx < 0 )
            isgn = -i4col_compare ( m, n, a, i, j );
        else if ( indx == 0 )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4col_sort2_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of A.
    Input, int N, the number of columns of A, and the length
    of a vector of data.
    Input/output, int A[M*N].
    On input, the array of N columns of M vectors.
    On output, the elements of each column of A have been sorted in ascending
    order.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ col;
    int i, j;
    int indx;
    short isgn;
    dim_typ temp;

    if ( m <= 1 || n <= 0 )
        return NULL;

    /*
    Initialize.
    */
    for ( col = 0; col < n; ++col)
    {
        i = indx = isgn = j = 0;
        /*
        Call the external heap sorter.
        */
        for ( ; ; )
        {
            sort_heap_external ( m, &indx, &i, &j, isgn );
            /*
            Interchange the I and J objects.
            */
            if ( 0 < indx )
            {
                temp       = a[i-1+col*m];
                a[i-1+col*m] = a[j-1+col*m];
                a[j-1+col*m] = temp;
            }
            /*
            Compare the I and J objects.
            */
            else if ( indx < 0 )
                isgn = -1 + ((a[j-1+col*m] < a[i-1+col*m])<<1);
            else if ( indx == 0 )
                break;
            }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4col_sort2_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_SORT2_D descending sorts the elements of each column of an I4COL.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of A.
    Input, int N, the number of columns of A, and the length
    of a vector of data.
    Input/output, int A[M*N].
    On input, the array of N columns of M vectors.
    On output, the elements of each column of A have been sorted in descending
    order.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ col;
    int i, j;
    int indx;
    short isgn;
    dim_typ temp;

    if ( m <= 1 || n <= 0 )
        return NULL;
    /*
    Initialize.
    */
    for ( col = 0; col < n; ++col )
    {
        i = indx = isgn = j = 0;
        /*
        Call the external heap sorter.
        */
        for ( ; ; )
        {
            sort_heap_external ( m, &indx, &i, &j, isgn );
            /*
            Interchange the I and J objects.
            */
            if ( 0 < indx )
            {
                temp       = a[i-1+col*m];
                a[i-1+col*m] = a[j-1+col*m];
                a[j-1+col*m] = temp;
            }
            /*
            Compare the I and J objects.
            */
            else if ( indx < 0 )
                isgn = -1 + ((a[i-1+col*m] < a[j-1+col*m]) <<1);
            else if ( indx == 0 )
                break;
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4col_sorted_singleton_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_SORTED_SINGLETON_COUNT counts singletons in an I4COL.
  Discussion:
    The columns of the array may be ascending or descending sorted.
    A "singleton" is an item that occurs exactly once.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 October 2005
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A[M*N], a sorted array, containing
    N columns of data.
    Output, int I4COL_SORTED_SINGLETON_COUNT, the number of singletons.
*/
{
	static int result = INT_MAX;
	
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ differ_from_next;
    dim_typ differ_from_previous;
    dim_typ i, j;
    dim_typ singleton_num;

    singleton_num = 0;

    if ( n <= 0 )
    {
    	result = singleton_num;
        return &result;
    }

    differ_from_next = 1;

    for ( j = 0; j < n; ++j )
    {
        differ_from_previous = differ_from_next;

        if ( j < n )
        {
            differ_from_next = 0;
            for ( i = 0; i < m; ++i )
                if ( a[i+j*m] != a[i+(j+1)*m] )
                {
                    differ_from_next = 1;
                    break;
                }
        }
        else
            differ_from_next = 1;

        if ( differ_from_previous && differ_from_next )
            ++ singleton_num;
    }

	
    result = singleton_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4col_sorted_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_SORTED_UNIQUE keeps unique elements in a sorted I4COL.
  Discussion:
    The array can be sorted into ascending or descending order.
    The important point is that identical elements must be stored
    in adjacent positions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of A, and the length of
    a vector of data.
    Input, int N, the number of columns of A.
    Input/output, int A[M*N].
    On input, the sorted array of N columns of M-vectors.
    On output, a sorted array of columns of M-vectors.
    Output, int *UNIQUE_NUM, the number of unique columns of A.
*/
{
	const _2dtpips * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	short * unique_num = s_data->a3;
	
    dim_typ i;
    dim_typ j1;
    dim_typ j2;
    dim_typ same;

    if ( n <= 0 )
    {
        *unique_num = 0;
        return NULL;
    }

    j1 = 1;

    for ( j2 = 2; j2 <= n; ++j2 )
    {
        same = 1;
        for ( i = 1; i <= m; ++i )
            if ( a[i-1+(j1-1)*m] != a[i-1+(j2-1)] )
            {
                same = 0;
                break;
            }
        if ( !same )
        {
            ++ j1;
            for ( i = 1; i <= m; i++ )
                a[i-1+(j1-1)*m] = a[i-1+(j2-1)*m];
        }
    }

    *unique_num = j1;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4col_sorted_unique_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
  Discussion:
    The columns of the array may be ascending or descending sorted.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A[M*N], a sorted array, containing
    N columns of data.
    Output, int I4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
*/
{
 	static int result = INT_MAX;
 	
 	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
 	
    dim_typ i;
    dim_typ j1;
    dim_typ j2;
    dim_typ unique_num;

    if ( n <= 0 )
    {
    	result = 0;
        return &result;
    }

    unique_num = 1;
    j1 = 0;

    for ( j2 = 1; j2 < n; ++j2 )
        for ( i = 0; i < m; ++i)
            if ( a[i+j1*m] != a[i+j2*m] )
            {
                ++ unique_num;
                j1 = j2;
                break;
            }

	result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4col_swap ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_SWAP swaps two columns of an I4COL.
  Discussion:
    The two dimensional information is stored as a one dimensional
    array, by columns.
    The row indices are 1 based, NOT 0 based.  However, a preprocessor
    variable, called OFFSET, can be reset from 1 to 0 if you wish to
    use 0-based indices.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, int A[M*N], an array of data.
    Input, int ICOL1, ICOL2, the two columns to swap.
    These indices should be between 1 and N.
*/
{
	const _2dtpi2dt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	const register dim_typ icol1 = s_data->a3;
	const register dim_typ icol2 = s_data->a4;
	
    # define OFFSET 1

    dim_typ i, t;
    /*
    Check.
    */
    if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET || icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET || icol1 == icol2 )
        return NULL;

    for ( i = 0; i < m; ++i)
    {
        t                     = a[i+(icol1-OFFSET)*m];
        a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
        a[i+(icol2-OFFSET)*m] = t;
    }

    return NULL;
    # undef OFFSET
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4col_unique_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4COL_UNIQUE_INDEX indexes the unique occurrence of values in an I4COL.
  Discussion:
    An I4COL is an M by N array of I4 values.
    It is regarded as an array of N columns of length M.
    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
    gathered in order, then
      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of A.
    The length of an "element" of A, and the number of "elements".
    Input, int A[M*N], the array.
    Output, int I4COL_UNIQUE_INDEX[N], the unique index.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ diff;
    dim_typ i;
    dim_typ j1;
    dim_typ j2;
    int *unique_index = ( int * ) malloc ( n * sizeof ( int ) );
    dim_typ unique_num;

    for ( j1 = 0; j1 < n; ++j1 )
        unique_index[j1] = -1;

    unique_num = 0;

    for ( j1 = 0; j1 < n; j1++ )
    {
        if ( unique_index[j1] == -1 )
        {
            unique_index[j1] = unique_num;

            for ( j2 = j1 + 1; j2 < n; j2++ )
            {
                diff = 0;
                for ( i = 0; i < m; ++i)
                    diff += abs ( a[i+j1*m] - a[i+j2*m] );
                if ( diff == 0 )
                    unique_index[j2] = unique_num;

            }
            ++ unique_num;
        }
    }
    return unique_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4i4_sort_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4I4_SORT_A ascending sorts a pair of I4's.
  Discussion:
    The program allows the reasonable call:
      i4i4_sort_a ( i1, i2, &i1, &i2 );
    and this will return the reasonable result.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int I1, I2, the values to sort.
    Output, int J1, J2, the sorted values.
*/
{
	const dti2pdt * const s_data = data;
	const register dim_typ i1 = s_data->a0;
	const register int i2 = s_data->a1;
	dim_typ * j1 = s_data->a2;
	dim_typ * j2 = s_data->a3;
	
    dim_typ k1 = i1;
    dim_typ k2 = i2;
    *j1 = MIN ( k1, k2 );
    *j2 = MAX ( k1, k2 );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4i4i4_sort_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4I4I4_SORT_A ascending sorts a triple of I4's.
  Discussion:
    The program allows the reasonable call:
      i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
    and this will return the reasonable result.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int I1, I2, I3, the values to sort.
    Output, int *J1, *J2, *J3, the sorted values.
*/
{
	const _3i3pdt * const s_data = data;
	const register int i1 = s_data->a0;
	const register int i2 = s_data->a1;
	const register int i3 = s_data->a2;
	dim_typ * j1 = s_data->a3;
	dim_typ * j2 = s_data->a4;
	dim_typ * j3 = s_data->a5;
	
    dim_typ  k1 = i1;
    dim_typ k2 = i2;
    dim_typ k3 = i3;
    *j1 = MIN ( MIN ( k1, k2 ), MIN ( k2, k3 ) );
    *j2 = MIN ( MAX ( k1, k2 ),
    MIN ( MAX ( k2, k3 ), MAX ( k3, k1 ) ) );
    *j3 = MAX ( MAX ( k1, k2 ), MAX ( k2, k3 ) );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4int_to_r4int ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4INT_TO_R8INT maps an I4 interval to an R4 interval.
  Discussion:
    The formula is
      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int IMIN, IMAX, the range.
    Input, int I, the integer to be converted.
    Input, float RMIN, RMAX, the range.
    Output, float R, the corresponding value in [RMIN,RMAX].
*/
{
	static float result = FLT_MAX;
	
	const _3i2f * const s_data = data;
	const register int imin = s_data->a0;
	const register int imax = s_data->a1;
	const register int i = s_data->a2;
	float rmin = s_data->a3;
	float rmax = s_data->a4;
	
	result = imax == imin ? 0.50 * ( rmin + rmax ) : ( ( float ) ( imax - i        ) * rmin+ ( float ) (        i - imin ) * rmax )/ ( float ) ( imax     - imin );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4int_to_r8int ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4INT_TO_R8INT maps an I4 interval to an R8 interval.
  Discussion:
    The formula is
      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2005
  Author:
    John Burkardt
  Parameters:
    Input, int IMIN, IMAX, the range.
    Input, int I, the integer to be converted.
    Input, double RMIN, RMAX, the range.
    Output, double R, the corresponding value in [RMIN,RMAX].
*/
{
	static ityp result = MAX_VAL;
	
	const _3i2it * const s_data = data;
	const register int imin = s_data->a0;
	const register int imax = s_data->a1;
	const register int i = s_data->a2;
	ityp rmin = s_data->a3;
	ityp rmax = s_data->a4;
	
	result = imax == imin ? 0.50 * ( rmin + rmax ) : ( ( ityp ) ( imax - i        ) * rmin+ ( ityp ) (        i - imin ) * rmax )/ ( ityp ) ( imax     - imin );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4mat_border_add ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_BORDER_ADD adds a "border" to an I4MAT.
  Discussion:
    We suppose the input data gives values of a quantity on nodes
    in the interior of a 2D grid, and we wish to create a new table
    with additional positions for the nodes that would be on the
    border of the 2D grid.
                  0 0 0 0 0 0
      * * * *     0 * * * * 0
      * * * * --> 0 * * * * 0
      * * * *     0 * * * * 0
                  0 0 0 0 0 0
    The illustration suggests the situation in which a 3 by 4 array
    is input, and a 5 by 6 array is to be output.
    The old data is shifted to its correct positions in the new array.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 May 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of points.
    Input, int TABLE[M*N], the table data.
    Output, int TABLE2[(M+2)*(N+2)], the augmented table data.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * table = s_data->a2;
	
    dim_typ i, j;
    int * table2 = ( int * ) malloc ( ( m + 2 ) * ( n + 2 ) * sizeof ( int ) );

    for ( j = 0; j < n+2; ++j)
        for ( i = 0; i < m+2; ++i )
            table2[i+j*(m+2)] = (!(i == 0 || i == m+1 || j == 0 || j == n+1))*table[(i-1)+(j-1)*m];

    return table2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4mat_border_cut ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_BORDER_CUT cuts the "border" of an I4MAT.
  Discussion:
    We suppose the input data gives values of a quantity on nodes
    on a 2D grid, and we wish to create a new table corresponding only
    to those nodes in the interior of the 2D grid.
      0 0 0 0 0 0
      0 * * * * 0    * * * *
      0 * * * * 0 -> * * * *
      0 * * * * 0    * * * *
      0 0 0 0 0 0
    The illustration suggests the situation in which a 5 by 6 array
    is input, and a 3 by 4 array is to be output.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 May 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of points.
    Input, int TABLE[M*N], the table data.
    Output, int TABLE2[(M-2)*(N-2)], the "interior" table data.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * table = s_data->a2;
	
    dim_typ i, j;
    int *table2;

    if ( m <= 2 || n <= 2 )
        return NULL;

    table2 = ( int * ) malloc ( ( m - 2 ) * ( n - 2 ) * sizeof ( int ) );

    for ( j = 0; j < n-2; ++j )
        for ( i = 0; i < m-2; ++i )
            table2[i+j*(m-2)] = table[(i+1)+(j+1)*m];

    return table2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_copy ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_COPY copies one I4MAT to another.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A1[M*N], the matrix to be copied.
    Output, int A2[M*N], the copy of A1.
*/
{
	const _2dt2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a1 = s_data->a2;
	int * a2 = s_data->a3;
	
    dim_typ i, j;
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            a2[i+j*m] = a1[i+j*m];
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  inline void   * _i4mat_copy_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_COPY_NEW copies an I4MAT to a "new" I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A1[M*N], the matrix to be copied.
    Output, int I4MAT_COPY_NEW[M*N], the copy of A1.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a1 = s_data->a2;
	
    int *a2 = ( int * ) malloc ( m * n * sizeof ( int ) );
    dim_typ i, j;
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i)
            a2[i+j*m] = a1[i+j*m];
    return a2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_elim ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_ELIM carries out exact Gauss elimination on an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input/output, int A[M*N].  On input, the M by N matrix to
    be Gauss eliminated.  On output, the Gauss-eliminated matrix.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ amax;
    dim_typ i;
    dim_typ *icol;
    dim_typ ifact;
    dim_typ imax;
    dim_typ imult;
    dim_typ *irow;
    dim_typ iswap;
    dim_typ j;
    dim_typ jcol;
    dim_typ jmult;

    icol = ( dim_typ * ) malloc ( n * sizeof ( dim_typ ) );
    irow = ( dim_typ * ) malloc ( m * sizeof ( dim_typ ) );
    /*
    Initialize the swap parity counter.
    */
    iswap = 1;
    /*
    For each column JCOL...
    */
    for ( jcol = 0; jcol < MIN ( m, n ); ++jcol )
    {
        /*
        Find the maximum element in rows JCOL through M.
        */
        amax = abs ( a[jcol+jcol*m] );
        imax = jcol;

        for ( i = jcol+1; i < m; ++i )
            if ( amax < abs ( a[i+jcol*m] ) )
            {
                amax = abs ( a[i+jcol*m] );
                imax = i;
            }
        /*
        If the maximum entry is nonzero, then...
        */
        if ( amax != 0 )
        {
            /*
            If the maximum entry does not occur in row JCOL, then swap rows.
            */
            if ( imax != jcol )
            {
                iswap = -iswap;
                i4vec_swap ( n, a+jcol*m, a+imax*m );
            }
            /*
            Eliminate all nonzero entries in column JCOL, below the diagonal entry.
            */
            for ( i = jcol+1; i < m; ++i )
            {
                if ( a[i+jcol*m] != 0 )
                {
                    jmult = a[i+jcol*m];
                    imult = a[jcol+jcol*m];
                    ifact = i4_gcd ( imult, jmult );
                    imult = imult / ifact;
                    jmult = jmult / ifact;

                    for ( j = jcol; j < n; ++j)
                        a[i+j*m] = jmult * a[jcol+j*m] - imult * a[i+j*m];
                }
            }
        /*
        Remove any row or column factors.
        */
        i4mat_red ( m, n, a, irow, icol );
        }
    }

    free ( icol );
    free ( irow );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4mat_flip_cols ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_FLIP_COLS swaps the columns of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    To "flip" the columns of an I4MAT is to start with something like
      11 12 13 14 15
      21 22 23 24 25
      31 32 33 34 35
      41 42 43 44 45
      51 52 53 54 55
    and return
      15 14 13 12 11
      25 24 23 22 21
      35 34 33 32 31
      45 44 43 42 41
      55 54 53 52 51
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, int A[M*N], the matrix whose columns are to be flipped.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ b, i, j;

    for ( i = 0; i < m; ++i )
        for ( j = 0; j < ( n / 2 ); ++j )
        {
            b              = a[i+     j *m];
            a[i+     j *m] = a[i+(n-1-j)*m];
            a[i+(n-1-j)*m] = b;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4mat_flip_rows ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_FLIP_ROWS swaps the rows of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    To "flip" the rows of an I4MAT is to start with something like
      11 12 13 14 15
      21 22 23 24 25
      31 32 33 34 35
      41 42 43 44 45
      51 52 53 54 55
    and return
      51 52 53 54 55
      41 42 43 44 45
      31 32 33 34 35
      21 22 23 24 25
      11 12 13 14 15
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, int A[M*N], the matrix whose rows are to be flipped.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ b, i, j;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < ( m / 2 ); ++i )
        {
            b            = a[    i+j*m];
            a[    i+j*m] = a[m-1-i+j*m];
            a[m-1-i+j*m] = b;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4mat_histogram ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_HISTOGRAM computes a histogram of the elements of an I4MAT.
  Discussion:
    An I4MAT is an array of I4's.
    It is assumed that the entries in the vector A are nonnegative.
    Only values between 0 and HISTO_NUM will be histogrammed.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the order of A.
    Input, int A[M*N], the array to examine.
    Input, int HISTO_NUM, the maximum value for which a
    histogram entry will be computed.
    Output, int I4MAT_HISTOGRAM[HISTO_NUM+1], contains the number of
    entries of A with the values of 0 through HISTO_NUM.
*/
{
	const _3dtpi * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register dim_typ histo_num = s_data->a2;
	int * a = s_data->a3;	
	
	dim_typ i, j;
    int *histo_gram = ( int * ) malloc ( ( histo_num + 1 ) * sizeof ( int ) );

    for (i = 0; i <= histo_num; ++i )
        histo_gram[i] = 0;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            if ( 0 <= a[i+j*m] && a[i+j*m] <= histo_num )
                histo_gram[a[i+j*m]] = histo_gram[a[i+j*m]] + 1;

    return histo_gram;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4mat_indicator_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_INDICATOR_NEW sets up an "indicator" I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    The value of each entry suggests its location, as in:
      11  12  13  14
      21  22  23  24
      31  32  33  34
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 May 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of the matrix.
    M must be positive.
    Input, int N, the number of columns of the matrix.
    N must be positive.
    Output, int I4MAT_INDICATOR_NEW[M*N], the table.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
    int *a;
    dim_typ fac;
    dim_typ i, j;

    a = ( int * ) malloc ( m * n * sizeof ( int ) );
    fac = powi ( 10, i4_log_10 ( n ) + 1 );

    for ( i = 1; i <= m; ++i )
        for ( j = 1; j <= n; ++j )
            a[i-1+(j-1)*m] = fac * i + j;

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4mat_l1_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_L1_INVERSE inverts a unit lower triangular I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    A unit lower triangular matrix is a matrix with only 1's on the main
    diagonal, and only 0's above the main diagonal.
    The inverse of an integer unit lower triangular matrix is also
    an integer unit lower triangular matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, number of rows and columns in the matrix.
    Input, int A[N*N], the unit lower triangular matrix.
    Output, int I4MAT_L1_INVERSE[N*N], the inverse matrix.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int *b;
    dim_typ i, j, k;
    b = ( int * ) malloc ( n * n * sizeof ( int ) );

    for ( i = 1; i <= n; ++i )
    {
        for ( j = 1; j < i; ++j )
        {
            b[i-1+(j-1)*n] = 0;
            for ( k = 1; k < i; k++ )
                b[i-1+(j-1)*n] += a[(i-1)+(k-1)*n] * b[(k-1)+(j-1)*n];
        }
        b[i-1+(i-1)*n] = 1;
        for ( j = i+1; j <= n; j++ )
            b[i-1+(j-1)*n] = 0;
    }

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_MAX returns the maximum of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, int A[M*N], the M by N matrix.
    Output, int I4MAT_MAX, the maximum entry of A.
*/
{
	static int result = INT_MAX;
	
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ i, j;
    int value;

    value = - i4_huge;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            if ( value < a[i+j*m] )
                value = a[i+j*m];
    
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_max_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_MAX_INDEX returns the location of the maximum of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, int A[M*N], the M by N matrix.
    Output, int *I_MAX, *J_MAX, the indices of the maximum entry of A.
*/
{
	const _2dtpi2pdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	dim_typ * i_max = s_data->a3;
	dim_typ * j_max = s_data->a4;
	
    dim_typ i, j;

    *i_max = *j_max = -1;

    for ( j = 1; j <= n; ++j )
        for ( i = 1; i <= m; ++i )
            if ( i == 1 && j == 1 )
            {
                *i_max = i;
                *j_max = j;
            }
            else if ( a[*i_max-1+(*j_max-1)*m] < a[i-1+(j-1)*m] )
            {
                *i_max = i;
                *j_max = j;
            }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_min ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_MIN returns the minimum of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, int A[M*N], the M by N matrix.
    Output, int I4MAT_MIN, the minimum entry of A.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    int value;

    value = i4_huge;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            if ( a[i+j*m] < value )
                value = a[i+j*m];
                
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_min_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_MIN_INDEX returns the location of the minimum of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, int A[M*N], the M by N matrix.
    Output, int *I_MIN, *J_MIN, the indices of the minimum entry of A.
*/
{
	const _2dtpi2pdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	dim_typ * i_min = s_data->a3;
	dim_typ * j_min = s_data->a4;
	
    dim_typ i, j;

    *i_min = -1;
    *j_min = -1;

    for ( j = 1; j <= n; ++j )
        for ( i = 1; i <= m; ++i )
            if ( i == 1 && j == 1 )
            {
                *i_min = i;
                *j_min = j;
            }
            else if ( a[i-1+(j-1)*m] < a[*i_min-1+(*j_min-1)*m] )
            {
                *i_min = i;
                *j_min = j;
            }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4mat_mm ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_MM multiplies two I4MAT's.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, N3, the order of the matrices.
    Input, int A[N1*N2], double B[N2*N#], the matrices to multiply.
    Output, int I4MAT_MM[N1*N3], the product matrix C = A * B.
*/
{
	const _3dt2pi * const s_data = data;
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	const register dim_typ n3 = s_data->a2;
	int * a = s_data->a3;
	int * b = s_data->a4;
	
    int *c;
    dim_typ i, j, k;
    c = ( int * ) malloc ( n1 * n3 * sizeof ( int ) );

    for ( i = 0; i < n1; ++i )
        for ( j = 0; j < n3; ++j )
        {
            c[i+j*n1] = 0;
            for ( k = 0; k < n2; ++k )
                c[i+j*n1] += a[i+k*n1] * b[k+j*n2];
        }

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_perm_uniform ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_PERM_UNIFORM selects a random permutation of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    The matrix is assumed to be square.  A single permutation is
    applied to both rows and columns.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 September 2008
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of rows and columns in the array.
    Input/output, int A[N*N], the N by N array to be permuted.
    Input/output, int *SEED, a seed for the random number generator.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * seed = s_data->a2;
	
    dim_typ i, j, k1, k2;
    int temp;
    /*
    Permute the rows and columns together.
    */
    for ( k1 = 0; k1 < n - 1; ++k1 )
    {
        k2 = i4_uniform_ab ( k1, n - 1, seed );
        for ( j = 0; j < n; ++j )
        {
            temp      = a[k1+j*n];
            a[k1+j*n] = a[k2+j*n];
            a[k2+j*n] = temp;
        }
        for ( i = 0; i < n; ++i)
        {
            temp      = a[i+k1*n];
            a[i+k1*n] = a[i+k2*n];
            a[i+k2*n] = temp;
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_perm2_uniform ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_PERM2_UNIFORM selects a random permutation of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    The matrix may be rectangular.  Separate permutations are
    applied to the rows and columns.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, int A[M*N], the M by N array to be permuted.
    Input/output, int *SEED, a seed for the random number generator.
*/
{
	const _2dt2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	int * seed = s_data->a3;
	
    dim_typ i, j, k1, k2;
    int temp;
    /*
    Permute the rows.
    */
    for ( k1 = 0; k1 < m; ++k1 )
    {
        k2 = i4_uniform_ab ( k1, m - 1, seed );
        for ( j = 0; j < n; j++ )
        {
            temp      = a[k1+j*m];
            a[k1+j*m] = a[k2+j*m];
            a[k2+j*m] = temp;
        }
    }
    /*
    Permute the columns.
    */
    for ( k1 = 0; k1 < n; ++k1 )
    {
    k2 = i4_uniform_ab ( k1, n - 1, seed );
        for ( i = 0; i < m; ++i)
        {
            temp      = a[i+k1*m];
            a[i+k1*m] = a[i+k2*m];
            a[i+k2*m] = temp;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4mat_red ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_RED divides out common factors in a row or column of an I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2005
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in the matrix.
    Input, int N, the number of columns in the matrix.
    Input/output, int A[M*N], on input, the M by N matrix to be reduced.
    On output, A has been reduced.  The greatest common factor in any
    row or column is 1.
    Output, int ROW[M], the row factors that were divided out.
    Output, int COL[N], the column factors that were divided out.
*/
{
	const _2dtpi2pdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	dim_typ * row = s_data->a3;
	dim_typ * col = s_data->a4;
	
    dim_typ i, j;
    int *temp;

    if ( m <= 0 || n <= 0 )
        return NULL;
    /*
    Remove factors common to a column.
    */
    for ( j = 0; j < n; ++j )
        col[j] = i4vec_red ( m, a+j*m );
    /*
    Remove factors common to a row.
    */
    temp = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < m; ++i)
    {
        for ( j = 0; j < n; ++j )
            temp[j] = a[i+j*m];
        row[i] = i4vec_red ( n, temp );
        for ( j = 0; j < n; j++ )
            a[i+j*m] = temp[j];
    }
    free ( temp );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4mat_u1_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_U1_INVERSE inverts a unit upper triangular I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    A unit upper triangular matrix is a matrix with only 1's on the main
    diagonal, and only 0's below the main diagonal.
    The inverse of an integer unit upper triangular matrix is also
    an integer unit upper triangular matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2008
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, number of rows and columns in the matrix.
    Input, int A[N*N], the unit upper triangular matrix.
    Output, int I4MAT_U1_INVERSE[N*N], the inverse matrix.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int *b;
    dim_typ i, j, k;

    b = ( int * ) malloc ( n * n * sizeof ( int ) );

    for ( j = n-1; 0 <= j; --j )
    {
        for ( i = j+1; i < n; i++ )
            b[i+j*n] = 0;
        b[j+j*n] = 1;

        for ( i = j-1; 0 <= i; --i )
        {
            b[i+j*n] = 0;
            for ( k = i+1; k <= j; ++k )
                b[i+j*n] -= a[i+k*n] * b[k+j*n];
        }
    }

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4mat_uniform_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_UNIFORM_AB returns a scaled pseudorandom I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    This routine implements the recursion
      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )
    The integer arithmetic never requires more than 32 bits,
    including a sign bit.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 May 2012
  Author:
    John Burkardt
  Reference:
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.
    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.
    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.
    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A, B, the limits of the pseudorandom values.
    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has
    been updated.
    Output, int X[M*N], a matrix of pseudorandom values.
*/
{
	const _2dt2i2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int a = s_data->a2;
	int b = s_data->a3;
	int * seed = s_data->a4;
	int * x = s_data->a5;
	
    dim_typ c, i, k, j;
    ityp r;
    int value;

    if ( *seed == 0 )
        return NULL;
    /*
    Guaranteee A <= B.
    */
    if ( b < a )
    {
        c = a;
        a = b;
        b = c;
    }

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
        {
            k = *seed / 127773;

            *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

            if ( *seed < 0 )
                *seed  += i4_huge;

            r = ( ityp ) ( *seed ) * 4.656612875E-10;
            /*
            Scale R to lie between A-0.5 and B+0.5.
            */
            r = ( 1.00 - r ) * ( ( ityp ) a - 0.50 )+         r   * ( ( ityp ) b + 0.50 );
            /*
            Use rounding to convert R to an integer between A and B.
            */
            value = round ( r );
            /*
            Guarantee A <= VALUE <= B.
            */
            value = value<a?a:b;

            x[i+j*m] = value;
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4mat_uniform_ab_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4MAT_UNIFORM_AB_NEW returns a scaled pseudorandom I4MAT.
  Discussion:
    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
    This routine implements the recursion
      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )
    The integer arithmetic never requires more than 32 bits,
    including a sign bit.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 May 2012
  Author:
    John Burkardt
  Reference:
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.
    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.
    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.
    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A, B, the limits of the pseudorandom values.
    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has
    been updated.
    Output, int I4MAT_UNIFORM_AB_NEW[M*N], a matrix of pseudorandom values.
*/
{
	const _2dt2ipi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int a = s_data->a2;
	int b = s_data->a3;
	int * seed = s_data->a4;
	
    dim_typ c, i;
    dim_typ j, k;
    ityp r;
    int value;
    int *x;

    if ( *seed == 0 )
        return NULL;
    /*
    Guaranteee A <= B.
    */
    if ( b < a )
    {
        c = a;
        a = b;
        b = c;
    }

    x = ( int * ) malloc ( m * n * sizeof ( int ) );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < m; ++i )
        {
            k = *seed / 127773;

            *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
            if ( *seed < 0 )
                *seed += i4_huge;


            r = ( ityp ) ( *seed ) * 4.656612875E-10;
            /*
            Scale R to lie between A-0.5 and B+0.5.
            */
            r = ( 1.00 - r ) * ( ( ityp ) a - 0.50 )+         r   * ( ( ityp ) b + 0.50 );
            /*
            Use rounding to convert R to an integer between A and B.
            */
            value = round ( r );
            /*
            Guarantee A <= VALUE <= B.
            */
            value = value<a?a:b;

            x[i+j*m] = value;
        }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void  *  _i4rmat_delete ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4RMAT_DELETE frees the memory set aside by I4RMAT_NEW.
  Discussion:
    This function releases the memory associated with an array that was
    created by a command like
      int **a;
      a = i4rmat_new ( m, n );
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int **A, the array.
    Input, int M, N, the number of rows and columns.
*/
{
	const ppi2dt * const s_data = data; 
	int ** a = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    for (dim_typ i = 0; i < m; ++i)
        free ( a[i] );
    free ( a );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4rmat_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4RMAT_NEW sets up an I4RMAT of the desired dimensions.
  Discussion:
    An I4RMAT is a row-major dynamically dimensioned 2D integer array.
    A declaration of the form
      int **a;
    is necesary.  Then an assignment of the form:
      a = i4rmat_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17;
      y = a[1][0];
    and so on.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 September 2013
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Output, int **I4RMAT_NEW, the array.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
    int **a;
    a = ( int ** ) malloc ( m * n * sizeof ( int * ) );
    for (dim_typ i = 0; i < m; ++i )
        a[i] = ( int * ) malloc ( n * sizeof ( int ) );
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4row_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_COMPARE compares two rows of a integer array.
  Discussion:
    The two dimensional information is stored in a one dimensional array,
    by columns.  The entry A(I,J) is stored in A[I+J*M].
    The input arguments I and J are row indices.  They DO NOT use the
    C convention of starting at 0, but rather start at 1.
  Example:
    Input:
      M = 3, N = 4, I = 2, J = 3
      A = (
        1  2  3  4
        5  6  7  8
        9 10 11 12 )
    Output:
      I4ROW_COMPARE = -1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int  A[M*N], the array of data.
    Input, int I, J, the rows to be compared.
    I and J must be between 1 and M.
    Output, int I4ROW_COMPARE, the results of the comparison:
    -1, row I < row J,
     0, row I = row J,
    +1, row J < row I.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpi2dt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	const register dim_typ i = s_data->a3;
	const register dim_typ j = s_data->a4;
	
    dim_typ k;
    /*
    Check that I and J are legal.
    */
    if ( i < 1 || m < i || j < 1 || m < j )
    {
    	result = 2;
        return &result;
    }

    if ( i == j )
    {
    	result = 0;
        return &result;
    }


    for ( k = 0; k < n; ++k)
        if ( a[(i-1)+k*m] < a[(j-1)+k*m] )
        {
        	result = -1;
        	return &result;
        }
        else if ( a[(j-1)+k*m] < a[(i-1)+k*m] )
        {
        	result = +1;
        	return &result;
        }

    result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4row_find_item ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_FIND_ITEM searches the rows of an I4ROW for a given value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A[M*N], the table to search.
    Input, int ITEM, the value to search for.
    Output, int ROW, COL, the row and column indices
    of the first occurrence of the value ITEM.  The search
    is conducted by rows.  If the item is not found, then
    ROW = COL = -1.
*/
{
	const _2dtpii2pdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	const register int item = s_data->a3;
	dim_typ * row = s_data->a4;
	dim_typ * col = s_data->a5;
	
    dim_typ i, j;
    *row = *col = -1;

    for ( i = 0; i < m; ++i )
        for ( j = 0; j < n; ++j)
            if ( a[i+j*m] == item )
            {
                *row = i + 1;
                *col = j + 1;
                return NULL;
            }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4row_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_MAX returns the maximums of an I4ROW.
  Example:
    A =
      1  2  3
      2  6  7
    MAX =
      3
      7
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in the array.
    Input, int A[M*N], the array to be examined.
    Output, int IMAX[M]; IMAX[I] is the column of A in which
    the maximum for row I occurs.
    Output, int I4ROW_MAX[M], the maximums of the rows.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ i, j;
    int *amax;

    amax = ( int * ) malloc ( m * sizeof ( int ) );

    for ( i = 0; i < m; ++i)
    {
        amax[i] = a[i+0*m];
        for ( j = 1; j < n; ++j )
            if ( amax[i] < a[i+j*m] )
                amax[i] = a[i+j*m];
    }

    return amax;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4row_mean ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_MEAN returns the means of an I4ROW.
  Example:
    A =
      1  2  3
      2  6  7
    MEAN =
      2
      5
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in the array.
    Input, int A[M*N], the array to be examined.
    Output, double I4ROW_MEAN[M], the means, or averages, of the rows.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ i, j;
    ityp *mean = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
    {
        mean[i] = 0.00;
        for ( j = 0; j < n; ++j )
            mean[i] += ( ityp ) a[i+j*m];
        mean[i] /= ( ityp ) ( n );
    }

    return mean;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4row_min ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_MIN returns the minimums of an I4ROW.
  Example:
    A =
      1  2  3
      2  6  7
    MIN =
      1
      2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in the array.
    Input, int A[M*N], the array to be examined.
    Output, int I4ROW_MIN[M], the minimums of the rows.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ i, j;
    int *amin = ( int * ) malloc ( m * sizeof ( int ) );

    for ( i = 0; i < m; ++i )
    {
        amin[i] = a[i+0*m];
        for ( j = 1; j < n; ++j)
            if ( a[i+j*m] < amin[i] )
                amin[i] = a[i+j*m];
    }

    return amin;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4row_sort_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_SORT_A ascending sorts the rows of an I4ROW.
  Discussion:
    In lexicographic order, the statement "X < Y", applied to two
    vectors X and Y of length M, means that there is some index I, with
    1 <= I <= M, with the property that
      X(J) = Y(J) for J < I,
    and
      X(I) < Y(I).
    In other words, X is less than Y if, at the first index where they
    differ, the X value is less than the Y value.
  Example:
    Input:
      M = 5, N = 3
      A =
        3  2  1
        2  4  3
        3  1  8
        2  4  2
        1  9  9
    Output:
      A =
        1  9  9
        2  4  2
        2  4  3
        3  1  8
        3  2  1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of A.
    Input, int N, the number of columns of A.
    Input/output, int A[M*N].
    On input, the array of M rows of N-vectors.
    On output, the rows of A have been sorted in ascending
    lexicographic order.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    /*
    Initialize.
    */
    i = indx = isgn = j = 0;
    /*
    Call the external heap sorter.
    */
    for ( ; ; )
    {
        sort_heap_external ( m, &indx, &i, &j, isgn );
        /*
        Interchange the I and J objects.
        */
        if ( 0 < indx )
            i4row_swap ( m, n, a, i, j );
        /*
        Compare the I and J objects.
        */
        else if ( indx < 0 )
            isgn = i4row_compare ( m, n, a, i, j );
        else if ( indx == 0 )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4row_sort_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_SORT_D descending sorts the rows of an I4ROW.
  Discussion:
    In lexicographic order, the statement "X < Y", applied to two real
    vectors X and Y of length M, means that there is some index I, with
    1 <= I <= M, with the property that
      X(J) = Y(J) for J < I,
    and
      X(I) < Y(I).
    In other words, the first time they differ, X is smaller.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows and columns of A.
    Input/output, int A[M*N].
    On input, the array of M rows of N-vectors.
    On output, the rows of A have been sorted in descending
    lexicographic order.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    /*
    Initialize.
    */
    i = indx = isgn = j = 0;
    /*
    Call the external heap sorter.
    */
    for ( ; ; )
    {
        sort_heap_external ( m, &indx, &i, &j, isgn );
        /*
        Interchange the I and J objects.
        */
        if ( 0 < indx )
            i4row_swap ( m, n, a, i, j );
        /*
        Compare the I and J objects.
        */
        else if ( indx < 0 )
            isgn = -i4row_compare ( m, n, a, i, j );
        else if ( indx == 0 )
            break;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4row_sort2_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_SORT2_D descending sorts the elements of each row of an I4ROW.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of A.
    Input, int N, the number of columns of A, and the length
    of a vector of data.
    Input/output, int A[M*N].
    On input, the array of M rows of N-vectors.
    On output, the elements of each row of A have been sorted in descending
    order.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    dim_typ row;
    int temp;

    if ( m <= 1 || n <= 0 )
        return NULL;
    /*
    Initialize.
    */
    for ( row = 0; row < m; ++row )
    {
        i = indx = isgn = 0;
        /*
        Call the external heap sorter.
        */
        for ( ; ; )
        {
            sort_heap_external ( n, &indx, &i, &j, isgn );
            /*
            Interchange the I and J objects.
            */
            if ( 0 < indx )
            {
                temp           = a[row+(i-1)*m];
                a[row+(i-1)*m] = a[row+(j-1)*m];
                a[row+(j-1)*m] = temp;
            }
            /*
            Compare the I and J objects.
            */
            else if ( indx < 0 )
                isgn = -1 + ((a[row+(i-1)*m] < a[row+(j-1)*m])<<1);
            else if ( indx == 0 )
                break;
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4row_sum ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_SUM returns the sums of the rows of an I4ROW.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, int A[M*N], the M by N array.
    Output, int I4ROW_SUM[M], the sum of the entries of
    each row.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ i, j;
    int *rowsum = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < m; ++i )
    {
        rowsum[i] = 0;
        for ( j = 0; j < n; ++j )
            rowsum[i] += a[i+j*m];
    }

    return rowsum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4row_swap ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_SWAP swaps two rows of an I4ROW.
  Discussion:
    The two dimensional information is stored as a one dimensional
    array, by columns.
    The row indices are 1 based, NOT 0 based  However, a preprocessor
    variable, called OFFSET, can be reset from 1 to 0 if you wish to
    use 0-based indices.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, int A[M*N], an array of data.
    Input, int IROW1, IROW2, the two rows to swap.
    These indices should be between 1 and M.
*/
{
	const _2dtpi2dt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	const register dim_typ irow1 = s_data->a3;
	const register dim_typ irow2 = s_data->a4;
	
    # define OFFSET 1

    dim_typ j, t;
    /*
    Check.
    */
    if ( irow1 < 0+OFFSET || m-1+OFFSET < irow1 || irow2 < 0+OFFSET || m-1+OFFSET < irow2 || irow1 == irow2 )
        return NULL;

    for ( j = 0; j < n; ++j )
    {
        t                   = a[irow1-OFFSET+j*m];
        a[irow1-OFFSET+j*m] = a[irow2-OFFSET+j*m];
        a[irow2-OFFSET+j*m] = t;
    }

    return NULL;
    # undef OFFSET
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4row_variance ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4ROW_VARIANCE returns the variances of an I4ROW.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in the array.
    Input, int A[M*N], the array whose variances are desired.
    Output, double I4ROW_VARIANCE[M], the variances of the rows.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * a = s_data->a2;
	
    dim_typ i, j;
    ityp mean;
    ityp *variance = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
    {
        mean = 0.00;
        for ( j = 0; j < n; ++j )
            mean = mean + ( ityp ) a[i+j*m];
        mean /= ( ityp ) ( n );

        variance[i] = 0.00;
        for ( j = 0; j < n; ++j )
        variance[i] += pow ( ( ( ityp ) a[i+j*m] - mean ), 2 );
        variance[i] = 0 + (1<n)*variance[i] / ( double ) ( n - 1 );
    }

    return variance;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_add ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ADD computes C = A + B for I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int A[N], the first vector.
    Input, int B[N], the second vector.
    Output, int C[N], the sum of the vectors.
*/
{
	const dt3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * b = s_data->a2;
	int * c = s_data->a3;
	
    for (dim_typ i = 0; i < n; ++i )
        c[i] = a[i] + b[i];
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_add_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ADD_NEW computes C = A + B for I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 April 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int A[N], the first vector.
    Input, int B[N], the second vector.
    Output, int I4VEC_ADD_NEW[N], the sum of the vectors.
*/
{
	const dtpipdt * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	int * b = s_data->a1;
	dim_typ * a = s_data->a2;
	
    int *c = ( int * ) malloc ( n * sizeof ( int ) );
    for (dim_typ i = 0; i < n; ++i )
        c[i] = a[i] + b[i];
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_all_nonpositive ( void *data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ALL_NONPOSITIVE: ( all ( A <= 0 ) ) for I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int A[N], the vector to check.
    Output, int I4VEC_ALL_NONPOSITIVE is 1 if all entries
    of A are less than or equal to zero.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( 0 < a[i] )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_amax ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_AMAX returns the largest magnitude in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector to be searched.
    Output, int I4VEC_AMAX, the value of the entry of largest magnitude.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int amax = a[0];

    for (dim_typ i = 1; i < n; ++i )
        if ( abs ( amax ) < abs ( a[i] ) )
            amax = a[i];
            
    result = amax;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_amax_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_AMAX_INDEX returns the index of the maximum absolute value in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int A[N], the array.
    Output, int I4VEC_AMAX_INDEX, the index of the entry of largest magnitude.
*/
{
	static short result = SHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int amax;
    short amax_index = 1;
    dim_typ i;

    if ( n <= 0 )
        amax_index = -1;
    else
    {
        amax_index = 1;
        amax = abs ( a[0] );

        for ( i = 2; i <= n; ++i)
            if ( amax < abs ( a[i-1] ) )
            {
                amax_index = i;
                amax = abs ( a[i-1] );
            }
    }
    
    result = amax_index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_amin ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_AMIN returns the smallest magnitude in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector to be checked.
    Output, int I4VEC_AMIN, the value of the smallest magnitude.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int amin;
    dim_typ i;
    amin = a[0];

    for ( i = 1; i < n; ++i     )
        if ( abs ( a[i] ) < abs ( amin ) )
            amin = a[i];

	result = amin;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_amin_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_AMIN_INDEX returns the index of the minimum absolute value in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, double A[N], the array.
    Output, int I4VEC_AMIN_INDEX, the index of the entry of smallest magnitude.
*/
{
	static short result = SHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int amin;
    short amin_index;
    dim_typ i;

    if ( n <= 0 )
        amin_index = -1;
    else
    {
        amin_index = 1;
        amin = abs ( a[0] );

        for ( i = 2; i <= n; ++i )
            if ( abs ( a[i-1] ) < amin )
            {
                amin_index = i;
                amin = abs ( a[i-1] );
            }
    }
    
    result = amin_index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_aminz ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_AMINZ returns the smallest nonzero magnitude in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector to be checked.
    Output, int I4VEC_AMINZ, the value of the smallest nonzero magnitude.
    If all entries are zero, I4VEC_AMINZ is 0.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int aminz = 0;
    dim_typ i;

    for ( i = 0; i < n; ++i)
        if ( a[i] )
        if ( aminz == 0 || abs ( a[i] ) < abs ( aminz ) )
            aminz = a[i];


	result = aminz;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_aminz_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_AMINZ_INDEX returns the smallest nonzero magnitude in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 September 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries to be checked.
    Input, int A[N], the vector to be checked.
    Output, int I4VEC_AMINZ_INDEX, the entry of the smallest nonzero magnitude.
    If all entries are zero, AMINZ_INDEX is 0.
*/
{
	static short result = SHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int aminz = 0;
    dim_typ i;
    short aminz_index;

    for ( i = 1; i <= n; ++i)
        if ( a[i-1] != 0 )
            if ( aminz_index == 0 || abs ( a[i-1] ) < aminz )
                aminz = abs ( a[i-1] );
                aminz_index = i;

	result = aminz_index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_any_lt ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ANY_LT: ( any ( A < B ) ) for I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 April 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int A[N], the first vector.
    Input, int B[N], the second vector.
    Output, int I4VEC_ANY_LT is 1 if any entry
    of A is less than the corresponding entry of B.
*/
{
	static bool result = 2;
	
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * b = s_data->a2;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( a[i] < b[i] )
        {
        	result = true;
            return &result;
        }
        
    result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_any_negative ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ANY_NEGATIVE: ( any A < 0 ) for I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int A[N], the vector to check.
    Output, int I4VEC_ANY_NEGATIVE is 1 if any entry is negative.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( a[i] < 0 )
        {
        	result = true;
            return &result;
        }
        
    result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_any_nonzero ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ANY_NONZERO: ( any A nonzero ) for I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int A[N], the vector to check.
    Output, int I4VEC_ANY_NONZERO is 1 if any entry is nonzero.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( a[i] != 0 )
        {
        	result = true;
            return &result;
        }
        
    result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_ascend_sub ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ASCEND_SUB computes the longest ascending subsequence of vector.
  Discussion:
    An I4VEC is a vector of I4's.
    The subsequence is required to be strictly increasing.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the length of the vector.
    Input,int A[N], the vector to be examined.
    Output, int *LENGTH, the length of the longest increasing subsequence.
    Output, int I4VEC_ASCEND_SUB[LENGTH], a longest increasing
    subsequence of A.
*/
{
	const dtpipdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	dim_typ * length = s_data->a2;
	
    dim_typ i;
    short j;
    dim_typ k;
    int *sub;
    int *top;
    int *top_prev;

    if ( n <= 0 )
        return NULL;

    top =top_prev = ( int * ) malloc ( n * sizeof ( int ) ); ( int * ) malloc ( n * sizeof ( int ) );
    top_prev = ( int * ) malloc ( n * sizeof ( int ) );
    for ( i = 0; i < n; ++i )
    {
        top[i] = 0;
        top_prev[i] = 0;
    }

    *length = 0;

    for ( i = 0; i < n; ++i)
    {
        k = -1;
        for ( j = 0; j < *length; ++j)
        {
            if ( a[i] <= a[top[j]] )
            {
                k = j;
                break;
            }
        }
        if ( k == -1 )
            k = *length ++;

        top[k] = i;
        top_prev[i] = 0<k ? top[k-1]:-1;
    }
    /*
    Extract the subsequence.
    */
    sub = ( int * ) malloc ( ( *length ) * sizeof ( int ) );

    j = top[*length-1];
    sub[*length-1] = a[j];

    for ( i = *length-2; 0 <= i; i-- )
    {
        j = top_prev[j];
        sub[i] = a[j];
    }

    free ( top );
    free ( top_prev );

    return sub;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_ascends ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
  Example:
    X = ( -8, 1, 2, 3, 7, 7, 9 )
    I4VEC_ASCENDS = TRUE
    The sequence is not required to be strictly ascending.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 April 2004
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the array.
    Input, int X[N], the array to be examined.
    Output, int I4VEC_ASCENDS, is TRUE if the entries of X ascend.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	
    for (dim_typ i = 1; i <= n-1; ++i )
        if ( x[i] < x[i-1] )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_axpy ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_AXPY adds IA times the vector IX to the vector IY.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 May 1999
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries of IX and IY.
    Input, int IA, the scalar value by which each entry
    of IX is multiplied before being added to IY.
    Input, int IX[*], the vector, a multiple of which is to be
    added to IY.
    Input, int INCX, the increment between successive entries of IX.
    Input/output, int IY[*].
    On output, each entry of IY has been increased by
    IA times the corresponding entry of IX.
    Input, int INCY, the increment between successive entries of IY.
*/
{
	const dtipiipii * const s_data = data;
	const register dim_typ n = s_data->a0;
	int ia = s_data->a1;
	int * ix = s_data->a2;
	int incx = s_data->a3;
	int * iy = s_data->a4;
	int incy = s_data->a5;
	
    dim_typ i;
    dim_typ indx;
    dim_typ indy;

    indx = indy = 0;

    for ( i = 0; i < n; ++i )
    {
        iy[indy] += ia * ix[indx];
        indx += incx;
        indy += incy;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_bracket ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_BRACKET searches a sorted I4VEC for successive brackets of a value.
  Discussion:
    An I4VEC is a vector of I4's.
    If the values in the vector are thought of as defining intervals
    on the number line, then this routine searches for the interval
    containing the given value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 May 1999
  Author:
    John Burkardt
  Parameters:
    Input, int N, length of input array.
    Input, int X[N], an array that has been sorted into ascending order.
    Input, int XVAL, a value to be bracketed.
    Output, int *LEFT, *RIGHT, the results of the search.
    In the most common case, 0 <= LEFT < LEFT + 1 = RIGHT <= N-1,
    and X[LEFT] <= XVAL <= X[RIGHT].
    Special cases:
      Value is less than all data values:
        LEFT = -1, RIGHT = 0, and XVAL < X[RIGHT].
      Value is greater than all data values:
        LEFT = N-1, RIGHT = -1, and X[LEFT] < XVAL.
      Value is equal to a data value:
        LEFT = RIGHT, and X[LEFT] = X[RIGHT] = XVAL.
      Algorithm failure:
        LEFT = RIGHT = -1.
*/
{
	const dtpii2ps * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	int xval = s_data->a2;
	short * left = s_data->a3;
	short * right = s_data->a4;
	
    dim_typ high;
    dim_typ low;
    dim_typ mid;
    /*
    XVAL < X[0].
    */
    if ( xval < x[0] )
    {
        *left = -1;
        *right = 0;
    }
    /*
    X[N] < XVAL.
    */
    else if ( x[n-1] < xval )
    {
        *left = n-1;
        *right = -1;
    }
    /*
    N = 1
    */
    else if ( n == 1 )
    {
        *left = 0;
        *right = 0;
    }
    /*
    X[0] <= XVAL <= X[N-1].
    */
    else
    {
        low = 0;
        high = n - 2;

        for ( ; ; )
        {
            mid = ( low + high ) >> 1;

            if ( high < low )
                return NULL;

            if ( x[mid] == xval )
            {
                *left = mid;
                *right = mid;
                return NULL;
            }
            else if ( x[mid+1] == xval )
            {
                *left = mid + 1;
                *right = mid + 1;
                return NULL;
            }
            else if ( x[mid] < xval && xval < x[mid+1] )
            {
                *left = mid;
                *right = mid + 1;
                return NULL;
            }
            else if ( x[mid+1] < xval )
                low = mid + 1;
            else if ( xval < x[mid] )
                high = mid - 1;
        }
    }
    
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_COMPARE compares two I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
    The lexicographic ordering is used.
  Example:
    Input:
      A = ( 2, 6, 2 )
      B = ( 2, 8, 12 )
    Output:
      ISGN = -1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, int A[N], B[N], the vectors to be compared.
    Output, int I4VEC_COMPARE, the results of the comparison:
    -1, A is lexicographically less than B,
     0, A is equal to B,
    +1, A is lexicographically greater than B.
*/
{
	static short result = SHRT_MAX;
	
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * b = s_data->a2;
	
    short isgn = 0;
    dim_typ k;
    for ( k = 0; k < n; ++k )
    {
        if ( a[k] < b[k] )
        {
        	result = -1;
            return &result;
        }
        else if ( b[k] < a[k] )
        {
        	result = 1;
            return &result;
        }
    }

	result = isgn;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_concatenate ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_CONCATENATE concatenates two I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 December 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N1, the number of entries in the first vector.
    Input, int A[N1], the first vector.
    Input, int N2, the number of entries in the second vector.
    Input, int B[N2], the second vector.
    Output, int C[N1+N2], the concatenated vector.
*/
{
	const _2dt3pi * const s_data = data;
	
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	int * a = s_data->a2;
	int * b = s_data->a3;
	int * c = s_data->a4;
	
    dim_typ i;

    for ( i = 0; i < n1; ++i )
        c[i] = a[i];
    for ( i = 0; i < n2; ++i )
        c[n1+i] = b[i];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_concatenate_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_CONCATENATE_NEW concatenates two I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 November 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N1, the number of entries in the first vector.
    Input, int A[N1], the first vector.
    Input, int N2, the number of entries in the second vector.
    Input, int B[N2], the second vector.
    Output, int I4VEC_CONCATENATE_NEW[N1+N2], the concatenated vector.
*/
{
	const _2dt2pi * const s_data = data;
	
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	int * a = s_data->a2;
	int * b = s_data->a3;
	
    dim_typ i;
    int *c = ( int * ) malloc ( ( n1 + n2 ) * sizeof ( int ) );

    for ( i = 0; i < n1; ++i )
        c[i] = a[i];
    for ( i = 0; i < n2; ++i )
        c[n1+i] = b[i];

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_copy ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_COPY copies an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 April 2007
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, int A1[N], the vector to be copied.
    Input, int A2[N], the copy of A1.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a1 = s_data->a1;
	int * a2 = s_data->a2;
	
    for (dim_typ i = 0; i < n; ++i)
        a2[i] = a1[i];
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_copy_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_COPY_NEW copies an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 July 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, int A1[N], the vector to be copied.
    Output, int I4VEC_COPY_NEW[N], the copy of A1.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a1 = s_data->a1;
	
    int *a2 = ( int * ) malloc ( n * sizeof ( int ) );
    for (dim_typ i = 0; i < n; ++i )
        a2[i] = a1[i];
    return a2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_cum_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_CUM_NEW computes the cumulative sum of the entries of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    Input:
      A = (/ 1, 2, 3, 4 /)
    Output:
      A_CUM = (/ 0, 1, 3, 6, 10 /)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 December 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector to be summed.
    Output, int I4VEC_CUM_NEW[N], the cumulative sum of the entries of A.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int *a_cum = ( int * ) malloc ( n * sizeof ( int ) );
    a_cum[0] = a[0];
    for (dim_typ i = 1; i <= n - 1; ++i )
        a_cum[i] = a_cum[i-1] + a[i];
    return a_cum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_cum0_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_CUM0_NEW computes the cumulative sum of the entries of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    This routine returns a vector of length N+1, with the first value
    being 0.
  Example:
    Input:
      A = (/ 1, 2, 3, 4 /)
    Output:
      A_CUM = (/ 0, 1, 3, 6, 10 /)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 December 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector to be summed.
    Output, int I4VEC_CUM0_NEW[N+1], the cumulative sum of the entries of A.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int *a_cum =( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
    a_cum[0] = 0;
    for (dim_typ i = 1; i <= n; ++i )
        a_cum[i] = a_cum[i-1] + a[i-1];
    return a_cum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_decrement ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_DECREMENT decrements an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input/output, int A[N], the vector to be decremented.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        -- a[i];
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_descends ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_DESCENDS determines if an I4VEC is (weakly) descending.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    X = ( 9, 7, 7, 3, 2, 1, -8 )
    I4VEC_DESCENDS = TRUE
    The sequence is not required to be strictly descending.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the array.
    Input, int A[N], the array to be examined.
    Output, int I4VEC_DESCENDS, is TRUE if the entries of X descend.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 1; i <= n-1; ++i )
        if ( a[i-1] < a[i] )
        {
        	result = false;
            return &result;
        }

    result = true;
    return &result;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_direct_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_DIRECT_PRODUCT creates a direct product of I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
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
    This routine carries out that task for the points X.
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
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR_INDEX, the index of the factor being processed.
    The first factor processed must be factor 0.
    Input, int FACTOR_ORDER, the order of the factor.
    Input, int FACTOR_VALUE[FACTOR_ORDER], the factor values for
    factor FACTOR_INDEX.
    Input, int FACTOR_NUM, the number of factors.
    Input, int POINT_NUM, the number of elements in the direct product.
    Input/output, int X[FACTOR_NUM*POINT_NUM], the elements of the
    direct product, which are built up gradually.  Before the first call,
    X might be set to 0.  After each factor has been input, X should
    have the correct value.
  Local Parameters:
    Local, integer START, the first location of a block of values to set.
    Local, integer CONTIG, the number of consecutive values to set.
    Local, integer SKIP, the distance from the current value of START
    to the next location of a block of values to set.
    Local, integer REP, the number of blocks of values to set.
*/
{
	const _2dtpi2dtpi * const s_data = data;
	const register dim_typ factor_index = s_data->a0;
	const register dim_typ factor_order = s_data->a1;
	int * factor_value = s_data->a2;
	const register dim_typ factor_num = s_data->a3;
	const register dim_typ point_num = s_data->a4;
	int * x = s_data->a5;
	
    static int contig = 0;
    dim_typ i, j, k;
    static int rep = 0;
    static int skip = 0;
    int start;

    if ( factor_index == 0 )
    {
        contig = 1;
        skip = 1;
        rep = point_num;
    }

    rep /= factor_order;
    skip *= factor_order;

    for ( j = 0; j < factor_order; ++j )
    {
        start = 0 + j * contig;

        for ( k = 1; k <= rep; ++k )
        {
            for ( i = start; i < start + contig; ++i )
                x[factor_index+i*factor_num] = factor_value[j];
            start += skip;
        }
    }

    contig *= factor_order;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_direct_product2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_DIRECT_PRODUCT2 creates a direct product of I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
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
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR_INDEX, the index of the factor being processed.
    The first factor processed must be factor 0.
    Input, int FACTOR_ORDER, the order of the factor.
    Input, int FACTOR_VALUE[FACTOR_ORDER], the factor values for
    factor FACTOR_INDEX.
    Input, int FACTOR_NUM, the number of factors.
    Input, int POINT_NUM, the number of elements in the direct product.
    Input/output, int W[POINT_NUM], the elements of the
    direct product, which are built up gradually.  Before the first call,
    W must be set to 1.
  Local Parameters:
    Local, integer START, the first location of a block of values to set.
    Local, integer CONTIG, the number of consecutive values to set.
    Local, integer SKIP, the distance from the current value of START
    to the next location of a block of values to set.
    Local, integer REP, the number of blocks of values to set.
*/
{
	const _2dtpi2dtpi * const s_data = data;
	const register dim_typ factor_index = s_data->a0;
	const register dim_typ factor_order = s_data->a1;
	int * factor_value = s_data->a2;
	const register dim_typ factor_num = s_data->a3;
	const register dim_typ point_num = s_data->a4;
	int * w = s_data->a5;
	
    static int contig = 0;
    dim_typ i, j, k;
    static int rep = 0;
    static int skip = 0;
    int start;

    if ( factor_index == 0 )
    {
        contig = skip = 1;
        rep = point_num;
    }

    rep /= factor_order;
    skip *= factor_order;

    for ( j = 0; j < factor_order; ++j)
    {
        start = 0 + j * contig;

        for ( k = 1; k <= rep; ++k )
        {
            for ( i = start; i < start + contig; ++i )
                w[i]*=factor_value[j];
            start += skip;
        }
    }

    contig *= factor_order;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_dot_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_DOT_PRODUCT computes the dot product of two I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 December 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the array.
    Input, int X[N], Y[N], the arrays.
    Output, int I4VEC_DOT_PRODUCT, the dot product of X and Y.
*/
{
	static int result = INT_MAX;
	
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	int * y = s_data->a2;
	
    int value = 0;
    for (dim_typ i = 0; i < n; ++i )
        value += x[i] * y[i];
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_eq ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_EQ is true if two I4VEC's are equal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 May 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, int A1[N], A2[N], two vectors to compare.
    Output, int I4VEC_EQ, is TRUE if every pair of elements A1(I) and A2(I)
    are equal, and FALSE otherwise.
*/
{
	static bool result = 2;
	
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a1 = s_data->a1;
	int * a2 = s_data->a2;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( a1[i] != a2[i] )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_even_all ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_EVEN_ALL is TRUE if all entries of an I4VEC are even.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector.
    Output, int I4VEC_EVEN_ALL, TRUE if all entries are even.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( ( a[i] % 2 ) == 1 )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_even_any ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_EVEN_ANY is TRUE if any entry of an I4VEC is even.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector.
    Output, int I4VEC_EVEN_ANY, TRUE if any entry is even.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( ( a[i] % 2 ) == 0 )
        {
        	result = true;
            return &result;
        }
        
    result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_find ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_FIND finds the first occurrence of a value in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 August 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, int A[N], the array.
    Input, int VALUE, the value being sought.
    Output, int I4VEC_FIND, the first location in A where
    VALUE occurs, or -1 if VALUE never occurs.
*/
{
	static short result = SHRT_MAX;
	
	const dtpii * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int value = s_data->a2; 
	
    short location = -1;
    for (short i = 0; i < n; ++i)

        if ( a[i] == value )
        {
        	result = i;
            return &result;
        }
        
	result = -1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_first_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_FIRST_INDEX indexes the first occurrence of values in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    For element A(I) of the vector, FIRST_INDEX(I) is the index in A of
    the first occurrence of the value A(I).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, int A[N], the array.
    Output, int I4VEC_FIRST_INDEX[N], the first occurrence index.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int *first_index;
    dim_typ i, j;
    first_index = ( int * ) malloc ( n * sizeof ( int ) );

    for (short i = 0; i < n; ++i )
        first_index[i] = -1;
    for (short i = 0; i < n; ++i )
        if ( first_index[i] == -1 )
        {
            first_index[i] = i;
            for ( j = i + 1; j < n; ++j)
                if ( a[i] == a[j] )
                    first_index[j] = i;
        }

    return first_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_frac ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_FRAC searches for the K-th smallest entry in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    Hoare's algorithm is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 September 2005
  Parameters:
    Input, int N, the number of elements of A.
    Input/output, int A[N].
    On input, A is the array to search.
    On output, the elements of A have been somewhat rearranged.
    Input, int K, the fractile to be sought.  If K = 1, the minimum
    entry is sought.  If K = N, the maximum is sought.  Other values
    of K search for the entry which is K-th in size.  K must be at
    least 1, and no greater than N.
    Output, double I4VEC_FRAC, the value of the K-th fractile of A.
*/
{
	static int result = INT_MAX;
	
	const _2dtpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	dim_typ k = s_data->a1; 
	int * a = s_data->a2;
	
    int frac;
    dim_typ i;
    dim_typ iryt;
    dim_typ j;
    dim_typ left;
    int temp;
    dim_typ x;

    if ( n <= 0 || k <= 0 || n < k)
    {
    	result = INT_MAX;
        return &result;
    }

    left = 1;
    iryt = n;

    for ( ; ; )
    {
        if ( iryt <= left )
        {
            frac = a[k-1];
            break;
        }

        x = a[k-1];
        i = left;
        j = iryt;

        for ( ; ; )
        {
            if ( j < i )
            {
                if ( j < k )
                    left = i;
                if ( k < i )
                    iryt = j;

                break;
            }
            /*
            Find I so that X <= A(I).
            */
            while ( a[i-1] < x )
                ++ i;
            /*
            Find J so that A(J) <= X.
            */
            while ( x < a[j-1] )
                -- j;

            if ( i <= j )
            {
                temp   = a[i-1];
                a[i-1] = a[j-1];
                a[j-1] = temp;
                ++ i;
                -- j;
            }
        }
    }

	result = frac;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_gcd ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_GCD returns the greatest common divisor of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The value GCD returned has the property that it is the greatest integer
    which evenly divides every entry of V.
    The entries in V may be negative.
    Any zero entries in V are ignored.  If all entries of V are zero,
    GCD is returned as 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of V.
    Input, int V[N], the vector.
    Output, int I4VEC_GCD, the greatest common divisor of V.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * v = s_data->a1;
	
    int gcd = 0;
    dim_typ i;

    for ( i = 0; i < n; ++i)
        if ( v[i] != 0 )
            gcd = gcd ? i4_gcd ( gcd, v[i] ) : abs(v[i]);
    /*
    If GCD is 0, that can only happen because all entries of V are zero.
    */
    
    result = gcd == 0 ? 1:gcd;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_heap_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_HEAP_A reorders an I4VEC into a ascending heap.
  Discussion:
    An I4VEC is a vector of I4's.
    An ascending heap is an array A with the property that, for every index J,
    A[J] <= A[2*J+1] and A[J] <= A[2*J+2], (as long as the indices
    2*J+1 and 2*J+2 are legal).
  Diagram:
                  A(0)
                /      \
            A(1)         A(2)
          /     \        /  \
      A(3)       A(4)  A(5) A(6)
      /  \       /   \
    A(7) A(8)  A(9) A(10)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2008
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the size of the input array.
    Input/Output, int A[N].
    On input, an unsorted array.
    On output, the array has been reordered into a heap.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    dim_typ ifree;
    int key;
    int m;
    /*
    Only nodes (N/2)-1 down to 0 can be "parent" nodes.
    */
    for ( i = ( n / 2 ) - 1; 0 <= i; --i )
    {
        /*
        Copy the value out of the parent node.
        Position IFREE is now "open".
        */
        key = a[i];
        ifree = i;

        for ( ; ; )
        {
            /*
            Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
            IFREE. (One or both may not exist because they equal or exceed N.)
            */
            m = (ifree<<1) + 1;
    /*
            Does the first position exist?
            */
            if ( n <= m )
                break;
            else
            {
                /*
                Does the second position exist?
                */
                if ( m + 1 < n )
                {
                    /*
                    If both positions exist, take the larger of the two values,
                    and update M if necessary.
                    */
                    if ( a[m+1] < a[m] )
                        ++ m;
                }
                /*
                If the large descendant is larger than KEY, move it up,
                and update IFREE, the location of the free position, and
                consider the descendants of THIS position.
                */
                if ( a[m] <= key )
                    break;

                a[ifree] = a[m];
                ifree = m;
            }
        }
        /*
        When you have stopped shifting items up, return the item you
        pulled out back to the heap.
        */
        a[ifree] = key;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4vec_heap_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
  Discussion:
    An I4VEC is a vector of I4's.
    A heap is an array A with the property that, for every index J,
    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
    2*J+1 and 2*J+2 are legal).
  Diagram:
                  A(0)
                /      \
            A(1)         A(2)
          /     \        /  \
      A(3)       A(4)  A(5) A(6)
      /  \       /   \
    A(7) A(8)  A(9) A(10)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 April 1999
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the size of the input array.
    Input/output, int A[N].
    On input, an unsorted array.
    On output, the array has been reordered into a heap.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    dim_typ ifree;
    int key;
    int m;
    /*
    Only nodes (N/2)-1 down to 0 can be "parent" nodes.
    */
    for ( i = ( n / 2 ) - 1; 0 <= i; --i)
    {
        /*
        Copy the value out of the parent node.
        Position IFREE is now "open".
        */
        key = a[i];
        ifree = i;

        for ( ; ; )
        {
            /*
            Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
            IFREE. (One or both may not exist because they equal or exceed N.)
            */
            m = (ifree<<1) + 1;
            /*
            Does the first position exist?
            */
            if ( n <= m )
                break;
            else
            {
                /*
                Does the second position exist?
                */
                if ( m + 1 < n )
                {
                    /*
                    If both positions exist, take the larger of the two values,
                    and update M if necessary.
                    */
                    if ( a[m] < a[m+1] )
                        ++ m;
                }
                /*
                If the large descendant is larger than KEY, move it up,
                and update IFREE, the location of the free position, and
                consider the descendants of THIS position.
                */
                if ( key < a[m] )
                {
                    a[ifree] = a[m];
                    ifree = m;
                }
                else
                    break;
            }
        }
        /*
        When you have stopped shifting items up, return the item you
        pulled out back to the heap.
        */
        a[ifree] = key;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_heap_d_extract ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_HEAP_D_EXTRACT extracts the maximum value from a descending heap.
  Discussion:
    An I4VEC is a vector of I4's.
    In other words, the routine finds the maximum value in the
    heap, returns that value to the user, deletes that value from
    the heap, and restores the heap to its proper form.
    This is one of three functions needed to model a priority queue.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2005
  Author:
    John Burkardt
  Reference:
    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, page 150.
  Parameters:
    Input/output, int *N, the number of items in the heap.
    Input/output, int A[N], the heap.
    Output, int VALUE, the item of maximum value, which has been
    removed from the heap.
*/
{
	static int result = INT_MAX;
	
	const pdtpi * const s_data = data;
	dim_typ * n = s_data->a0;
	int * a = s_data->a1;
	
    int value;

    if ( *n < 1 )
    {
    	result = INT_MAX;
        return &result;
    }
    /*
    Get the maximum value.
    */
    value = a[0];

    if ( *n == 1 )
    {
        *n = 0;
        result = value; 
        return &result;
    }
    /*
    Shift the last value down.
    */
    a[0] = a[*n-1];
    /*
    Restore the heap structure.
    */
    -- *n;
    i4vec_sort_heap_d ( *n, a );
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4vec_heap_d_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_HEAP_D_INSERT inserts a new value into a descending heap.
  Discussion:
    An I4VEC is a vector of I4's.
    This is one of three functions needed to model a priority queue.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Reference:
    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, page 150.
  Parameters:
    Input/output, int *N, the number of items in the heap.
    Input/output, int A[N], the heap.
    Input, int VALUE, the value to be inserted.
*/
{
	const pdtpii * const s_data = data;
	dim_typ * n = s_data->a0;
	int * a = s_data->a1;
	int value = s_data->a2;
	
    dim_typ i;
    int parent;

    ++ *n;
    i = *n;

    while ( 1 < i )
    {
        parent = i / 2;

        if ( value <= a[parent-1] )
            break;
        a[i-1] = a[parent-1];
        i = parent;
    }
    a[i-1] = value;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_heap_d_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_HEAP_D_MAX returns the maximum value in a descending heap of integers.
  Discussion:
    An I4VEC is a vector of I4's.
    This is one of three functions needed to model a priority queue.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Reference:
    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, page 150.
  Parameters:
    Input, int N, the number of items in the heap.
    Input, int A[N], the heap.
    Output, int I4VEC_HEAP_D_MAX, the maximum value in the heap.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
		
	result = a[0];
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_histogram ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_HISTOGRAM computes a histogram of the elements of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    It is assumed that the entries in the vector A are nonnegative.
    Only values between 0 and HISTO_NUM will be histogrammed.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, int A[N], the array to examine.
    Input, int HISTO_NUM, the maximum value for which a
    histogram entry will be computed.
    Output, int I4VEC_HISTOGRAM[HISTO_NUM+1], contains the number of
    entries of A with the values of 0 through HISTO_NUM.
*/
{
	const dtpii * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int histo_num = s_data->a2; 
	
    int *histo_gram;
    dim_typ i;

    histo_gram = ( int * ) malloc ( ( histo_num + 1 ) * sizeof ( int ) );

    for ( i = 0; i <= histo_num; ++i )
        histo_gram[i] = 0;

    for ( i = 0; i < n; ++i )
        if ( 0 <= a[i] && a[i] <= histo_num )
            ++ histo_gram[a[i]];
    return histo_gram;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_histogram_masked ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_HISTOGRAM_MASKED computes a histogram of the elements of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    A histogram will be made for all vector entries with MASK = 1,
    whose values are between 0 and HISTO_NUM - 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, int A[N], the array to examine.
    Input, int HISTO_NUM, the maximum value for which a
    histogram entry will be computed.
    Input, int MASK[N], the mask value for each entry.
    Output, int I4VEC_HISTOGRAM[HISTO_NUM], counts the number
    of entries with MASK = 1, with values of 0 through HISTO_NUM - 1.
*/
{
	const dti2pi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	int histo_num = s_data->a1; 
	int * a = s_data->a2;
	int * mask = s_data->a3;
	
    int *histo_gram;
    dim_typ i;

    histo_gram = ( int * ) malloc ( histo_num * sizeof ( int ) );

    for ( i = 0; i < histo_num; ++i )
    histo_gram[i] = 0;

    for ( i = 0; i < n; ++i )
        if ( mask[i] == 1 )
            if ( 0 <= a[i] && a[i] < histo_num )
                ++ histo_gram[a[i]];
    return histo_gram;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4vec_increment ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INCREMENT increments an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input/output, int A[N], the vector to be incremented.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++ a[i ++] );
    
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX returns the location of the first occurrence of a given value.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 July 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector to be searched.
    Input, int AVAL, the value to be indexed.
    Output, int I4VEC_INDEX, the first location in A which has the
    value AVAL, or -1 if no such index exists.
*/
{
	static short result = SHRT_MAX;
	
	const dtpii * const s_data = data; 
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int aval = s_data->a2; 	
	
    for (dim_typ i = 0; i < n; ++i )
        if ( a[i] == aval )
        {
        	result = i;
            return &result;
        }
        
    result = -1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_index_delete_all ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX_DELETE_ALL deletes all occurrences of a value from an indexed sorted list.
  Discussion:
    An I4VEC is a vector of I4's.
    Note that the value of N is adjusted because of the deletions
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input, int X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Input, int XVAL, the value to be sought.
    Output, int *N2, the size of the current list.
    Output, int X2[N2], the list.
    Output, int INDX2[N2], the sort index of the list.
*/
{
	const dt2pii3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	int * indx = s_data->a2;
	int xval = s_data->a3;
	int * n2 = s_data->a4;
	int * x2 = s_data->a5;
	int * indx2 = s_data->a6;
	
    int equal;
    int equal1;
    int equal2;
    int get;
    dim_typ i;
    int less;
    int more;
    int put;

    if ( n < 1 )
    {
        *n2 = 0;
        return NULL;
    }

    i4vec_copy ( n, indx, indx2 );
    i4vec_copy ( n, x, x2 );
    *n2 = n;

    i4vec_index_search ( *n2, x2, indx2, xval, &less, &equal, &more );

    if ( equal == 0 )
        return NULL;

    equal1 = equal;

    for ( ; ; )
    {
        if ( equal1 <= 1 )
            break;

        if ( x2[indx2[equal1-2]-1] != xval )
            break;
        -- equal1;
    }

    equal2 = equal;

    for ( ; ; )
    {
        if ( *n2 <= equal2 )
            break;

        if ( x2[indx2[equal2]-1] != xval )
            break;

        ++ equal2;
    }
    /*
    Discard certain X values.
    */
    put = 0;

    for ( get = 1; get <= *n2; ++get)
        if ( x2[get-1] != xval )
        {
            ++ put;
            x2[put-1] = x2[get-1];
        }
    /*
    Adjust the INDX values.
    */
    for ( equal = equal1; equal <= equal2; ++equal )
        for ( i = 1; i <= *n2; ++i )
            if ( indx2[equal-1] < indx2[i-1] )
                -- indx2[i-1];
        /*
        Discard certain INDX values.
        */
        for ( i = 0; i <= *n2 - equal2 - 1; ++i)
            indx2[equal1+i-1] = indx2[equal2+i];
        for ( i = *n2 + equal1 - equal2; i <= *n2; ++i)
            indx2[i-1] = 0;
    /*
    Adjust N.
    */
    *n2 = put;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_index_delete_dupes ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The output quantities N2, X2, and INDX2 are computed from the
    input quantities by sorting, and eliminating duplicates.
    The output arrays should be dimensioned of size N, unless the user
    knows in advance what the value of N2 will be.
    The output arrays may be identified with the input arrays.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the input list.
    Input, int X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Output, int *N2, the number of unique entries in X.
    Output, int X2[N2], a copy of the list which has
    been sorted, and made unique.
    Output, int INDX2[N2], the sort index of the new list.
*/
{
	const dtpipdt2pips * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	dim_typ * indx = s_data->a2;
	int * n2 = s_data->a3;
	int * x2 = s_data->a4;
	short * indx2 = s_data->a5; 
	
    dim_typ i;
    int n3;
    int *x3;

    i = n3 = 0;
    x3 = ( int * ) malloc ( n * sizeof ( int ) );

    for ( ; ; )
    {
        ++ i;

        if ( n < i )
            break;

        if ( 1 < i )
            if ( x[indx[i-1]-1] == x3[n3-1] )
                continue;

        ++ n3;
        x3[n3-1] = x[indx[i-1]-1];
    }
    /*
    Set the output data.
    */
    *n2 = n3;
    i4vec_copy ( n3, x3, x2 );
    for ( i = 0; i < n3; ++i )
        indx2[i] = i + 1;

    free ( x3 );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_index_delete_one ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX_DELETE_ONE deletes one copy of an I4 from an indexed sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    If the value occurs in the list more than once, only one copy is deleted.
    Note that the value of N is adjusted because of the deletions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input, int X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Input, int XVAL, the value to be sought.
    Output, int *N2, the size of the current list.
    Output, int X2[N2], the list.
    Output, int INDX2[N2], the sort index of the list.
*/
{
	const dt2pii3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	int * indx = s_data->a2;
	int xval = s_data->a3;
	int * n2 = s_data->a4;
	int * x2 = s_data->a5;
	int * indx2 = s_data->a6;
	
    int equal;
    dim_typ i, j;
    int less;
    int more;

    if ( n < 1 )
    {
        *n2 = 0;
        return NULL;
    }

    *n2 = n;
    i4vec_copy ( *n2, indx, indx2 );
    i4vec_copy ( *n2, x, x2 );

    i4vec_index_search ( *n2, x2, indx2, xval, &less, &equal, &more );

    if ( equal != 0 )
    {
        j = indx2[equal-1];
        for ( i = j; i <= *n2-1; ++i )
            x2[i-1] = x[i];
        for ( i = equal; i <= *n2-1; ++i )
            indx2[i-1] = indx2[i];
        for ( i = 1; i <= *n2 - 1; i++ )
            if ( j < indx2[i-1] )
                -- indx2[i-1];
        -- *n2;
    }

    return NULL;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_index_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX_INSERT inserts an I4 into an indexed sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N, the size of the current list.
    Input, int X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Input, int XVAL, the value to be sought.
*/
{
	const _3pii * const s_data = data;
	int * n = s_data->a0;
	int * x = s_data->a1;
	int * indx = s_data->a2;
	int xval = s_data->a3;
	
    int equal;
    dim_typ i;
    int less;
    int more;

    if ( *n <= 0 )
    {
        *n = indx[0] = 1;
        x[0] = xval;
        return NULL;
    }

    i4vec_index_search ( *n, x, indx, xval, &less, &equal, &more );

    x[*n] = xval;
    for ( i = *n; more <= i; --i )
        indx[i] = indx[i-1];
    indx[more-1] = *n + 1;
    ++ *n;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_index_insert_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX_INSERT_UNIQUE inserts a unique I4 in an indexed sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N, the size of the current list.
    If the input value XVAL does not already occur in X, then N is increased.
    Input/output, int X[N], the list.
    If the input value XVAL does not already occur in X, then it is added
    to X.
    Input/output, int INDX[N], the sort index of the list.
    If the input value XVAL does not already occur in X, then INDX is updated.
    Input, int XVAL, the value which will be inserted into the X
    vector if it is not there already.
*/
{
	const _3pii * const s_data = data;
	int * n = s_data->a0;
	int * x = s_data->a1;
	int * indx = s_data->a2;
	int xval = s_data->a3;
	
    int equal;
    dim_typ i;
    int less;
    int more;

    if ( *n <= 0 )
    {
        *n = indx[0] = 1;
        x[0] = xval;
        return NULL;
    }
    /*
    Does XVAL already occur in X?
    */
    i4vec_index_search ( *n, x, indx, xval, &less, &equal, &more );

    if ( equal == 0 )
    {
        x[*n] = xval;
        for ( i = *n; more <= i; --i)
            indx[i] = indx[i-1];
        indx[more-1] = *n + 1;
        ++ *n;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_index_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX_ORDER sorts an I4VEC using an index vector.
  Discussion:
    An I4VEC is a vector of I4's.
    The index vector itself is not modified.  Therefore, the pair
 (X,INDX) no longer represents an index sorted vector.  If this
    relationship is to be preserved, then simply set INDX(1:N)=(1:N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input/output, int X[N], the list.  On output, the list
    has been sorted.
    Input, int INDX[N], the sort index of the list.
*/
{
	const dtpipdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	dim_typ * indx = s_data->a2;
	
    dim_typ i;
    int *y = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
        y[i] = x[indx[i]-1];
    for ( i = 0; i < n; ++i )
        x[i] = y[i];
    free ( y );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_index_search ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX_SEARCH searches for an I4 in an indexed sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input, int X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Input, int XVAL, the value to be sought.
    Output, int *LESS, *EQUAL, *MORE, the indexes in INDX of the
    entries of X that are just less than, equal to, and just greater
    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
    is the greatest entry of X, then MORE is N+1.
*/
{
	const dt2pii3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	int * indx = s_data->a2;
	int xval = s_data->a3;
	int * less = s_data->a4;
	int * equal = s_data->a5;
	int * more = s_data->a6;
	
    int hi;
    int lo;
    int mid;
    int xhi;
    int xlo;
    int xmid;

    if ( n <= 0 )
    {
        *less = *equal = *more = 0;
        return NULL;
    }

    lo = 1;
    hi = n;
    xlo = x[indx[lo-1]-1];
    xhi = x[indx[hi-1]-1];

    if ( xval < xlo )
    {
        *less = *equal = 0;
        *more = 1;
        return NULL;
    }
    else if ( xval == xlo )
    {
        *less = 0;
        *equal = 1;
        *more = 2;
        return NULL;
    }

    if ( xhi < xval )
    {
        *less = n;
        *equal = 0;
        *more = n + 1;
        return NULL;
    }
    else if ( xval == xhi )
    {
        *less = n - 1;
        *equal = n;
        *more = n + 1;
        return NULL;
    }

    for ( ; ; )
    {
        if ( lo + 1 == hi )
        {
            *less = lo;
            *equal = 0;
            *more = hi;
            return NULL;
        }

        mid = ( lo + hi ) / 2;
        xmid = x[indx[mid-1]-1];

        if ( xval == xmid )
        {
            *equal = mid;
            *less = mid - 1;
            *more = mid + 1;
            return NULL;
        }
        else if ( xval < xmid )
            hi = mid;
        else if ( xmid < xval )
            lo = mid;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_index_sort_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEX_SORT_UNIQUE creates a sort index for an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input, int X[N], the list.
    Output, int *N2, the number of unique elements in X.
    Output, int X2[N2], a list of the unique elements of X.
    Output, int INDX2[N2], the sort index of the list.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	int * n2 = s_data->a2;
	int * x2 = s_data->a3;
	int * indx2 = s_data->a4;
	
    dim_typ i;
    *n2 = 0;

    for ( i = 0; i < n; ++i )
        i4vec_index_insert_unique ( n2, x2, indx2, x[i] );
    for ( i = *n2; i < n; ++i)
        x2[i] = -1;
    for ( i = *n2; i < n; ++i )
        indx2[i] = -1;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_indexed_heap_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEXED_HEAP_D creates a descending heap from an indexed I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
    each referencing an entry of the data vector.
    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
    we have:
      A[INDX[2*J+1]]   <= A[INDX[J]]
    and
      A[INDX[2*J+2]] <= A[INDX[J]]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2010
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.
  Parameters:
    Input, int N, the size of the index array.
    Input, int A[*], the data vector.
    Input/output, int INDX[N], the index array.
    Each entry of INDX must be a valid index for the array A.
    On output, the indices have been reordered into a descending heap.
*/
{
	const dtpipdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	dim_typ * indx = s_data->a2;
	
    dim_typ i;
    int ifree;
    int key;
    int m;
    /*
    Only nodes N/2 - 1 down to 0 can be "parent" nodes.
    */
    for ( i = ( n / 2 ) - 1; 0 <= i; --i )
    {
        /*
        Copy the value out of the parent node.
        Position IFREE is now "open".
        */
        key = indx[i];
        ifree = i;

        for ( ; ; )
        {
        /*
        Positions 2*IFREE+1 and 2*IFREE+2 are the descendants of position
        IFREE. (One or both may not exist because they exceed N-1.)
        */
            m = (ifree<<1)+ 1;
            /*
            Does the first position exist?
            */
            if ( n - 1 < m )
                break;
            /*
            Does the second position exist?
            */
            if ( m + 1 <= n - 1 )
            {
                /*
                If both positions exist, take the larger of the two values,
                and update M if necessary.
                */
                if ( a[indx[m]] < a[indx[m+1]] )
                    ++ m;
            }
            /*
            If the large descendant is larger than KEY, move it up,
            and update IFREE, the location of the free position, and
            consider the descendants of THIS position.
            */
            if ( a[indx[m]] <= a[key] )
                break;

            indx[ifree] = indx[m];
            ifree = m;
        }
        /*
        Once there is no more shifting to do, KEY moves into the free spot IFREE.
        */
        indx[ifree] = key;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_indexed_heap_d_extract ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
    each referencing an entry of the data vector.
    The routine finds the maximum value in the heap, returns that value to the
    user, deletes that value from the heap, and restores the heap to its
    proper form.
    Note that the argument N must be a variable, which will be decremented
    before return, and that INDX will hold one less value on output than it
    held on input.
    This is one of three functions needed to model a priority queue.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2010
  Author:
    John Burkardt
  Reference:
    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, 2001,
    ISBN: 0262032937,
    LC: QA76.C662.
  Parameters:
    Input/output, int *N, the number of items in the index vector.
    Input, int A[*], the data vector.
    Input/output, int INDX[N], the index vector.
    Output, int I4VEC_INDEXED_HEAP_D_EXTRACT, the index in A of the item of
    maximum value, which has now been removed from the heap.
*/
{
	static int result = INT_MAX;
	
	const _2pipdt * const s_data = data;
	int * n = s_data->a0;
	int * a = s_data->a1;
	dim_typ * indx = s_data->a2;
	
    int indx_extract;

    if ( *n < 1 )
    {
    	result = INT_MAX;
        return &result;
    }
    /*
    Get the index of the maximum value.
    */
    indx_extract = indx[0];

    if ( *n == 1 )
    {
        *n = 0;
        result = indx_extract;
        return &result;
    }
    /*
    Shift the last index down.
    */
    indx[0] = indx[*n-1];
    /*
    Restore the heap structure.
    */
    -- *n;
    i4vec_indexed_heap_d ( *n, a, indx );
    
    result = indx_extract;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_indexed_heap_d_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
    each referencing an entry of the data vector.
    Note that the argument N must be a variable, and will be incremented before
    return, and that INDX must be able to hold one more entry on output than
    it held on input.
    This is one of three functions needed to model a priority queue.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2010
  Author:
    John Burkardt
  Reference:
    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, 2001,
    ISBN: 0262032937,
    LC: QA76.C662.
  Parameters:
    Input/output, int *N, the number of items in the index vector.
    Input, int A[*], the data vector.
    Input/output, int INDX[N], the index vector.
    Input, int INDX_INSERT, the index in A of the value
    to be inserted into the heap.
*/
{
	const dtpipdtpi * const s_data = data;
	
	dim_typ indx_insert = s_data->a0;
	int * n = s_data->a1;
	dim_typ * indx = s_data->a2;
	int * a = s_data->a3;
	
    dim_typ i;
    int parent;

    ++ *n;
    i = *n - 1;

    while ( 0 < i )
    {
        parent = ( i - 1 ) / 2;

        if ( a[indx_insert] <= a[indx[parent]] )
            break;

        indx[i] = indx[parent];
        i = parent;
    }

    indx[i] = indx_insert;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_indexed_heap_d_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    An indexed I4VEC is an I4VEC of data values, and an I4VEC of N indices,
    each referencing an entry of the data vector.
    This is one of three functions needed to model a priority queue.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2010
  Author:
    John Burkardt
  Reference:
    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, 2001,
    ISBN: 0262032937,
    LC: QA76.C662.
  Parameters:
    Input, int N, the number of items in the index vector.
    Input, int A[*], the data vector.
    Input, int INDX[N], the index vector.
    Output, int I4VEC_INDEXED_HEAP_D_MAX, the index in A of the maximum value
    in the heap.
*/
{
	static int result = INT_MAX;
	
	const _2pipdt * const s_data = data;
	int * n = s_data->a0;
	int * a = s_data->a1;
	dim_typ * indx = s_data->a2;
	
	result = indx[0];
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_indicator0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDICATOR0 sets an I4VEC to the indicator vector (0,1,2,...).
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, int A[N], the array.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        a[i] = i;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_indicator0_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDICATOR0_NEW sets an I4VEC to the indicator vector (0,1,2,...).
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, int I4VEC_INDICATOR0_NEW[N], the array.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    int * a = ( int * ) malloc ( n * sizeof ( int ) );
    for (dim_typ i = 0; i < n; ++i )
        a[i] = i;
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_indicator1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDICATOR1 sets an I4VEC to the indicator vector (1,2,3,...).
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, int A[N], the array.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        a[i] = i + 1;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_indicator1_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INDICATOR1_NEW sets an I4VEC to the indicator vector (1,2,3,...).
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, int I4VEC_INDICATOR1_NEW[N], the array.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    int * a = ( int * ) malloc ( n * sizeof ( int ) );
    for (dim_typ i = 0; i < n; ++i )
        a[i] = i + 1;
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_INSERT inserts a value into an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 November 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the array on input.
    Input/output, int A[N+1], the array.  On input, A is assumed
    to contain N entries.  On output, A actually contains N+1 entries.
    Input, int POS, the position to be assigned the new entry.
    0 <= POS <= N.
    Input, int VALUE, the value to be inserted at the given position.
*/
{
	const dtpi2i * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int pos = s_data->a2;
	int value = s_data->a3;
	
    int i;

    if ( pos < 0 || n < pos )
        return NULL;
    else
    {
        for ( i = n; pos < i; --i )
            a[i] = a[i-1];
        a[pos] = value;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_lcm ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_LCM returns the least common multiple of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The value LCM returned has the property that it is the smallest integer
    which is evenly divisible by every element of V.
    The entries in V may be negative.
    If any entry of V is 0, then LCM is 0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 July 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of V.
    Input, int V[N], the vector.
    Output, int I4VEC_LCM, the least common multiple of V.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * v = s_data->a1;
	
    int lcm = 1;
    for (dim_typ i = 0; i < n; ++i )
    {
        if ( v[i] == 0 )
        {
            lcm = 0;
            break;
        }
        lcm = i4_lcm ( lcm, v[i] );
    }
    
    result = lcm;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MAX returns the value of the maximum element in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int A[N], the array to be checked.
    Output, int IVEC_MAX, the value of the maximum element.  This
    is set to 0 if N <= 0.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    int value;

    if ( n <= 0 )
    {
    	result = 0;
        return &result;
	}
	
    value = a[0];

    for ( i = 1; i < n; ++i )
        if ( value < a[i] )
            value = a[i];

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_max_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MAX_INDEX returns the index of the maximum value in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int A[N], the array.
    Output, int I4VEC_MAX_INDEX, the index of the largest entry.
*/
{
	static short result = SHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    short max_index;

    if ( n <= 0 )
        max_index = -1;
    else
    {
        max_index = 0;

        for ( i = 1; i < n; ++i )
            if ( a[max_index] < a[i] )
                max_index = i;
    }
    
    result = max_index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_max_index_last ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MAX_INDEX_LAST: index of the last maximum value in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int A[N], the array.
    Output, int I4VEC_MAX_INDEX_LAST, the index of the last largest entry.
*/
{
	static short result = SHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    short max_index_last;

    if ( n <= 0 )
        max_index_last = -1;
    else
    {
        max_index_last = 0;
        for ( i = n-1; 0 <= i; --i )
            if ( a[max_index_last] < a[i] )
                max_index_last = i;
    }

	result = max_index_last;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_mean ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MEAN returns the mean of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 May 1999
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int X[N], the vector whose mean is desired.
    Output, double I4VEC_MEAN, the mean, or average, of the vector entries.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	
    ityp mean = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        mean += ( ityp ) x[i];
        
    result = mean /( ityp ) n;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_median ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MEDIAN returns the median of an unsorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    Hoare's algorithm is used.  The values of the vector are
    rearranged by this routine.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified
    18 September 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input/output, int A[N], the array to search.  On output,
    the order of the elements of A has been somewhat changed.
    Output, int I4VEC_MEDIAN, the value of the median of A.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
	result = i4vec_frac ( n, a, ( n + 1 ) / 2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_merge_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MERGE_A merges two ascending sorted I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
    The elements of A and B should be sorted in ascending order.
    The elements in the output array C will also be in ascending order,
    and unique.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int NA, the dimension of A.
    Input, int A[NA], the first sorted array.
    Input, int NB, the dimension of B.
    Input, int B[NB], the second sorted array.
    Output, int C[NA+NB], the merged unique sorted array.
*/
{
	const _2dt2pi * const s_data = data;
	
	const register dim_typ na = s_data->a0;
	const register dim_typ nb = s_data->a1;
	int * a = s_data->a2;
	int * b = s_data->a3;
	
    int *c;
    dim_typ j;
    dim_typ ja;
    dim_typ jb;
    dim_typ nc;
    dim_typ order;

    ja = jb = nc = 0;

    order = i4vec_order_type ( na, a );

    if ( order < 0 || 2 < order )
        return NULL;

    order = i4vec_order_type ( nb, b );

    if ( order < 0 || 2 < order )
        return NULL;

    c = ( int * ) malloc ( ( na + nb ) * sizeof ( int ) );

    for ( ; ; )
    {
        /*
        If we've used up all the entries of A, stick the rest of B on the end.
        */
        if ( na <= ja )
        {
            for ( j = 1; j <= nb - jb; j++ )
            {
                if ( nc == 0 || c[nc] < b[jb] )
                {
                    c[nc] = b[jb];
                    ++ nc;
                }
                ++ jb;
            }
            break;
        }
        /*
        If we've used up all the entries of B, stick the rest of A on the end.
        */
        else if ( nb <= jb )
        {
            for ( j = 1; j <= na - ja; j++ )
            {
                if ( nc == 0 || c[nc] < a[ja] )
                {
                    c[nc] = a[ja];
                    ++ nc;
                }
                ++ ja;
            }
            break;
        }
        /*
        Otherwise, if the next entry of A is smaller, that's our candidate.
        */
        else if ( a[ja] <= b[jb] )
        {
            if ( nc == 0 || c[nc] < a[ja] )
            {
                c[nc] = a[ja];
                ++ nc;
            }
            ++ ja;
        }
        /*
        ...or if the next entry of B is the smaller, consider that.
        */
        else
        {
            if ( nc == 0 || c[nc] < b[jb] )
            {
                c[nc] = b[jb];
                ++ nc;
            }
            ++ jb;
        }
    }
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_min ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MIN returns the minimum element in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 May 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int A[N], the array to be checked.
    Output, int I4VEC_MIN, the value of the minimum element.  This
    is set to 0 if N <= 0.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int value;

    if ( n <= 0 )
    {
    	result = 0;
    	return &result;
    }

    value = a[0];

    for (dim_typ i = 1; i < n; ++i )
        if ( a[i] < value )
            value = a[i];
            
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_min_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MIN_INDEX returns the index of the minimum value in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int A[N], the array.
    Output, int I4VEC_MIN_INDEX, the index of the smallest entry.
*/
{
	static short result = SHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    short min_index;

    if ( n <= 0 )
        min_index = -1;
    else
    {
        min_index = 0;
        for ( i = 1; i < n; ++i )
            if ( a[i] < a[min_index] )
                min_index = i;
    }

	result = min_index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_min_mv ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_MIN_MV determines U(1:N) /\ V for vectors U and a single vector V.
  Discussion:
    For two vectors U and V, each of length M, we define
  ( U /\ V ) (I) = MIN ( U(I), V(I) ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 January 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of the vectors.
    Input, int N, the number of vectors in U.
    Input, int U[M*N], N vectors, each of length M.
    Input, int V[M], a vector of length M.
    Output, int W[M*N], the value of U /\ W.
*/
{
	const _2dt3pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * u = s_data->a2;
	int * v = s_data->a3;
	int * w = s_data->a4;
	
    dim_typ i, j;
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            w[i+j*m] = MIN ( u[i+j*m], v[i] );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_negone ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_NEGONE sets an I4VEC to -1.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, int A[N], a vector of -1's.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i)
        a[i] = -1;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_negone_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_NEGONE_NEW creates an I4VEC and sets it to -1.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 October 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, int I4VEC_NEGONE_NEW[N], a vector of -1's.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    int *a = ( int * ) malloc ( n * sizeof ( int ) );
    for (dim_typ i = 0; i < n; ++i )
        a[i] = -1;
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_nonzero_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_NONZERO_COUNT counts the nonzero entries in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the input array.
    Input, int A[N], an array.
    Output, int I4VEC_NONZERO_COUNT, the number of nonzero entries.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ nonzero_count = 0;
    for (dim_typ i = 0; i < n; ++i )
        if ( a[i] != 0 )
            ++ nonzero_count;
            
    result = nonzero_count;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_nonzero_first ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_NONZERO_FIRST left-shifts all nonzeros in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The routine preserves the ordering of the nonzero entries.  It counts
    the nonzeros, and returns an index vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input/output, int X[N], the vector to be shifted.
    Output, int *NZ, the number of nonzero entries in the vector.
    Output, int INDX[N], contains the original location of each entry.
*/
{
	const dtpipdtpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	dim_typ * indx = s_data->a2;
	int * nz = s_data->a3;
	
    dim_typ j, k;
    *nz = 0;

    for ( j = 1; j <= n; ++j )
        indx[j-1] = j;

    j = 0;

    while ( j < n )
    {
        ++ j;

        if ( x[j-1] != 0 )
        {
            ++ *nz;
            if ( *nz != j )
            {
                x[*nz-1] = x[j-1];
                x[j-1] = 0;

                k = indx[*nz-1];
                indx[*nz-1] = j;
                indx[j-1] = k;
            }
        }
    }
    return NULL; 
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_norm_l0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_NORM_L0 returns the l0 "norm" of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The l0 "norm" simply counts the number of nonzero entries in the vector.
    It is not a true norm, but has some similarities to one.  It is useful
    in the study of compressive sensing.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 June 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector.
    Output, int I4VEC_NORM_L0, the value of the norm.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ value = 0;
    for (dim_typ i = 0; i < n; ++i )
        if ( a[i] != 0 )
        ++ value;
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_odd_all ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ODD_ALL is TRUE if all entries of an I4VEC are odd.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector.
    Output, int I4VEC_ODD_ALL, TRUE if all entries are odd.
*/
{
	static bool result = 2;

	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( ( a[i] % 2 ) == 0 )
        {
        	result = false;
            return &result;
        }
        
	result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_odd_any ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ODD_ANY is TRUE if any entry of an I4VEC is odd.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector.
    Output, int I4VEC_ODD_ANY, TRUE if any entry is odd.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( ( a[i] % 2 ) == 1 )
        {
        	result = true;
            return &result;
        }
    
    result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_one_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ONE_NEW creates an I4VEC whose entries are 1.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 October 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, int I4VEC_ONE_NEW[N], a vector of zeroes.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    int * a = ( int * ) malloc ( n * sizeof ( int ) );
    for (dim_typ i = 0; i < n; i++ )
        a[i] = 1;
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_order_type ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ORDER_TYPE: is an I4VEC is (non)strictly ascending/descending?
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries of the array.
    Input, int X[N], the array to be checked.
    Output, int I4VEC_ORDER_TYPE, order indicator:
    -1, no discernable order;
    0, all entries are equal;
    1, ascending order;
    2, strictly ascending order;
    3, descending order;
    4, strictly descending order.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	
    dim_typ i;
    dim_typ order;
    /*
    Search for the first value not equal to X(0).
    */
    i = 0;

    for ( ; ; )
    {
        ++ i;
        if ( n-1 < i )
        {
        	result = 0;
            return &result;
        }

        if ( x[0] < x[i] )
        {
            if ( i == 1 )
            {
                order = 2;
                break;
            }
            else
            {
                order = 1;
                break;
            }
        }
        else if ( x[i] < x[0] )
        {
            if ( i == 1 )
            {
                order = 4;
                break;
            }
            else
            {
                order = 3;
                break;
            }
        }
    }
    /*
    Now we have a "direction".  Examine subsequent entries.
    */
    for ( ; ; )
    {
        ++ i;
        if ( n - 1 < i )
            break;

        if ( order == 1 )
        {
            if ( x[i] < x[i-1] )
            {
                order = -1;
                break;
            }
        }
        else if ( order == 2 )
        {
            if ( x[i] < x[i-1] )
            {
                order = -1;
                break;
            }
            else if ( x[i] == x[i-1] )
                order = 1;
        }
        else if ( order == 3 )
        {
            if ( x[i-1] < x[i] )
            {
                order = -1;
                break;
            }
        }
        else if ( order == 4 )
        {
            if ( x[i-1] < x[i] )
            {
                order = -1;
                break;
            }
            else if ( x[i] == x[i-1] )
                order = 3;
        }
    }
    
    result = order;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_pairwise_prime ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_PAIRWISE_PRIME checks whether an I4VEC is pairwise prime.
  Discussion:
    An I4VEC is a vector of I4's.
    Two positive integers I and J are pairwise prime if they have no common
    factor greater than 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of values to check.
    Input, int A[N], the vector of integers.
    Output, int I4VEC_PAIRWISE_PRIME, is TRUE if the vector of integers
    is pairwise prime.
*/
{
	static bool result = 2;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;
    for ( i = 0; i < n; ++i )
        for ( j = i+1; j < n; ++j )
            if ( i4_gcd ( a[i], a[j] ) != 1 )
            {
            	result = false;
                return &result;
            }
            
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_part ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_PART partitions an int NVAL into N nearly equal parts.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    Input:
      N = 5, NVAL = 17
    Output:
      X = ( 4, 4, 3, 3, 3 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int NVAL, the integer to be partitioned.
    NVAL may be positive, zero, or negative.
    Output, int X[N], the partition of NVAL.  The entries of
    X add up to NVAL.  The entries of X are either all equal, or
    differ by at most 1.  The entries of X all have the same sign
    as NVAL, and the "largest" entries occur first.
*/
{
	const dtpii * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	int nval = s_data->a2;
	
    dim_typ i, j;

    for ( i = 0; i < n; ++i )
        x[i] = 0;

    if ( 0 < nval )
    {
        j = 0;
        for ( i = 0; i < nval; ++i )
        {
            ++ x[j];
            ++ j;
            if ( n <= j )
                j = 0;
        }
    }
    else if ( nval < 0 )
    {
        j = 0;
        for ( i = nval; i < 0; ++i)
        {
            -- x[j];
            ++ j;
            if ( n <= j )
                j = 0;
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_part_quick_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_PART_QUICK_A reorders an I4VEC as part of a quick sort.
  Discussion:
    An I4VEC is a vector of I4's.
    I4VEC_PART_QUICK_A reorders the entries of A.  Using A[0] as a
    key, all entries of A that are less than or equal to A[0] will
    precede A[0] which precedes all entries that are greater than A[0].
  Example:
    Input:
      N = 8
      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
    Output:
      L = 3, R = 6
      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
            -------        -------
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries of A.
    Input/output, int A[N].  On input, the array to be checked.
    On output, A has been reordered as described above.
    Output, int L, R, the indices of A that define the three segments.
    Let KEY = the input value of A[0].  Then
    I <= L             A(I) < KEY;
     L < I < R         A(I) = KEY;
             R <= I    A(I) > KEY.
*/
{
	const dt3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * l = s_data->a2;
	int * r = s_data->a3;
	
    dim_typ i;
    int key;
    int m;
    int temp;

    if ( n < 1 )
        return NULL;
    else if ( n == 1 )
    {
        *l = 0;
        *r = 2;
        return NULL;
    }

    key = a[0];
    m = 1;
    /*
    The elements of unknown size have indices between L+1 and R-1.
    */
    *l = 1;
    *r = n + 1;

    for ( i = 2; i <= n; ++i )
    {

        if ( key < a[*l] )
        {
            -- *r;
            temp = a[*r-1];
            a[*r-1] = a[*l];
            a[*l] = temp;
        }
        else if ( a[*l] == key )
        {
            ++ m;
            temp = a[m-1];
            a[m-1] = a[*l];
            a[*l] = temp;
            ++ *l;
        }
        else if ( a[*l] < key )
            ++ *l;

    }
    /*
    Now shift small elements to the left, and KEY elements to center.
    */
    for ( i = 1; i <= *l -m; ++i )
        a[i-1] = a[i+m-1];

    *l -= m;

    for ( i = *l+1; i <= *l+m; ++i )
        a[i-1] = key;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_permute ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_PERMUTE permutes an I4VEC in place.
  Discussion:
    An I4VEC is a vector of I4's.
    This routine permutes an array of integer "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.
  Example:
    Input:
      N = 5
      P = (   1,   3,   4,   0,   2 )
      A = (   1,   2,   3,   4,   5 )
    Output:
      A    = (   2,   4,   5,   1,   3 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input, int P[N], the permutation.  P(I) = J means
    that the I-th element of the output array should be the J-th
    element of the input array.
    Input/output, int A[N], the array to be permuted.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	int * a = s_data->a2;
	
    int a_temp;
    dim_typ i;
    int ierror;
    int iget;
    int iput;
    int istart;

    ierror = perm0_check ( n, p );

    if ( ierror != 0 )
        return NULL;
    /*
    In order for the sign negation trick to work, we need to assume that the
    entries of P are strictly positive.  Presumably, the lowest number is 0.
    So temporarily add 1 to each entry to force positivity.
    */
    for ( i = 0; i < n; ++i )
        ++ p[i];
    /*
    Search for the next element of the permutation that has not been used.
    */
    for ( istart = 1; istart <= n; ++istart )
    {
        if ( p[istart-1] < 0 )
            continue;
        else if ( p[istart-1] == istart )
        {
            p[istart-1] *= -1;
            continue;
        }
        else
        {
            a_temp = a[istart-1];
            iget = istart;
            /*
            Copy the new value into the vacated entry.
            */
            for ( ; ; )
            {
                iput = iget;
                iget = p[iget-1];

                p[iput-1] = - p[iput-1];

                if ( iget < 1 || n < iget )
                    return NULL;
                if ( iget == istart )
                {
                    a[iput-1] = a_temp;
                    break;
                }
                a[iput-1] = a[iget-1];
            }
        }
    }
    /*
    Restore the signs of the entries.
    */
    for ( i = 0; i < n; ++i )
        p[i] = -p[i]-1;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_permute_uniform ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_PERMUTE_UNIFORM randomly permutes an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input/output, int A[N], the array to be permuted.
    Input/output, int *SEED, a seed for the random number generator.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int * seed = s_data->a2;
	
    int *p = p = perm0_uniform_new ( n, seed );
    i4vec_permute ( n, p, a );
    free ( p );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_red ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_RED divides out common factors in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input/output, int A[N], the vector to be reduced.
    On output, the entries have no common factor
    greater than 1.
    Output, int I4VEC_RED, the common factor that was divided out.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int factor;
    dim_typ i;
    /*
    Find the smallest nonzero value.
    */
    factor = 0;

    for ( i = 0; i < n; ++i )
        if ( a[i])
            factor = abs(a[i]);

    if ( factor == 0 )
    {
    	result = factor;
        return &result;
    }
    /*
    Find the greatest common factor of the entire vector.
    */
    for ( i = 0; i < n; ++i )
        factor = i4_gcd ( a[i], factor );

    if ( factor == 1 )
    {
    	result = factor;
        return &result;
    }
    /*
    Divide out the common factor.
    */
    for ( i = 0; i < n; ++i )
        a[i] /= factor;

    result = factor;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_reverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_REVERSE reverses the elements of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    Input:
      N = 5,
      A = ( 11, 12, 13, 14, 15 ).
    Output:
      A = ( 15, 14, 13, 12, 11 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input/output, int A[N], the array to be reversed.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int j;

    for (dim_typ i = 0; i < n / 2; ++i )
    {
        j        = a[i];
        a[i]     = a[n-1-i];
        a[n-1-i] = j;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_rotate ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ROTATE rotates an I4VEC in place.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    Input:
      N = 5, M = 2
      X    = ( 1, 2, 3, 4, 5 )
    Output:
      X    = ( 4, 5, 1, 2, 3 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input, int M, the number of positions to the right that
    each element should be moved.  Elements that shift pass position
    N "wrap around" to the beginning of the array.
    Input/output, int X[N], the array to be rotated.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * x = s_data->a2;
	
    int iget;
    int iput;
    int istart;
    int mcopy;
    int nset;
    int temp;
    /*
    Force M to be positive, between 0 and N-1.
    */
    mcopy = i4_modp ( m, n );

    if ( mcopy == 0 )
        return NULL;

    istart = nset = 0;

    for ( ; ; )
    {
        if ( n <= istart )
            return NULL;

        temp = x[istart];
        iget = istart;
        /*
        Copy the new value into the vacated entry.
        */
        for ( ; ; )
        {
            iput = iget;

            iget = iget - mcopy;

            if ( iget < 0 )
                iget += n;

            if ( iget == istart )
                break;
            x[iput] = x[iget];
            ++ nset;
        }

        x[iput] = temp;
        ++ nset;

        if ( n <= nset )
            break;


        ++ istart;

    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_run_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_RUN_COUNT counts runs of equal values in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    A run is a sequence of equal values.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector to be examined.
    Output, int I4VEC_RUN_COUNT, the number of runs.
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    int run_count;
    int test;

    run_count = 0;

    if ( n < 1 )
    {
    	result = 0;
       return &result;
   }

    test = 0;

    for ( i = 0; i < n; ++i)
        if ( i == 0 || a[i] != test )
        {
            ++ run_count;
            test = a[i];
        }

	result = run_count;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sort_bubble_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_BUBBLE_A ascending sorts an I4VEC using bubble sort.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, length of input array.
    Input/output, int A[N].
    On input, an unsorted array of ints.
    On output, A has been sorted.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;
    int temp;

    for ( i = 0; i < n - 1; ++i )
        for ( j = i + 1; j < n; ++j )
            if ( a[j] < a[i] )
            {
                temp = a[i];
                a[i] = a[j];
                a[j] = temp;
            }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_sort_bubble_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_BUBBLE_D descending sorts an I4VEC using bubble sort.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, length of input array.
    Input/output, int A[N].
    On input, an unsorted array of ints.
    On output, A has been sorted.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;
    int temp;

    for ( i = 0; i < n - 1; ++i )
        for ( j = i + 1; j < n; ++j )
            if ( a[i] < a[j] )
            {
                temp = a[i];
                a[i] = a[j];
                a[j] = temp;
            }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sort_heap_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 April 1999
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of entries in the array.
    Input/output, int A[N].
    On input, the array to be sorted;
    On output, the array has been sorted.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ n1;
    int temp;

    if ( n <= 1 )
        return NULL;
    /*
    1: Put A into descending heap form.
    */
    i4vec_heap_d ( n, a );
    /*
    2: Sort A.

    The largest object in the heap is in A[0].
    Move it to position A[N-1].
    */
    temp = a[0];
    a[0] = a[n-1];
    a[n-1] = temp;
    /*
    Consider the diminished heap of size N1.
    */
    for ( n1 = n - 1; 2 <= n1; --n1)
    {
        /*
        Restore the heap structure of the initial N1 entries of A.
        */
        i4vec_heap_d ( n1, a );
        /*
        Take the largest object from A[0] and move it to A[N1-1].
        */
        temp = a[0];
        a[0] = a[n1-1];
        a[n1-1] = temp;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sort_heap_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_HEAP_D descending sorts an I4VEC using heap sort.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 September 2008
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of entries in the array.
    Input/output, int A[N].
    On input, the array to be sorted;
    On output, the array has been sorted.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ n1;
    int temp;

    if ( n <= 1 )
        return NULL;
    /*
    1: Put A into ascending heap form.
    */
    i4vec_heap_a ( n, a );
    /*
    2: Sort A.

    The smallest object in the heap is in A[0].
    Move it to position A[N-1].
    */
    temp = a[0];
    a[0] = a[n-1];
    a[n-1] = temp;
    /*
    Consider the diminished heap of size N1.
    */
    for ( n1 = n - 1; 2 <= n1; --n1 )
    {
        /*
        Restore the heap structure of the initial N1 entries of A.
        */
        i4vec_heap_a ( n1, a );
        /*
        Take the smallest object from A[0] and move it to A[N1-1].
        */
        temp = a[0];
        a[0] = a[n1-1];
        a[n1-1] = temp;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_sort_heap_index_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.
    Once the index array is computed, the sorting can be carried out
    "implicitly:
      a(indx(*))
    or explicitly, by the call
      i4vec_permute ( n, indx, a )
    after which a(*) is sorted.
    Note that the index vector is 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int A[N], an array to be index-sorted.
    Output, int I4VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int aval;
    dim_typ i;
    int *indx;
    int indxt;
    int ir;
    dim_typ j, l;

    if ( n < 1 )
        return NULL;

    indx = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
        indx[i] = i;

    if ( n == 1 )
        return indx;

    l = n / 2 + 1;
    ir = n;

    for ( ; ; )
    {

        if ( 1 < l )
        {
            -- l;
            indxt = indx[l-1];
            aval = a[indxt];
        }
        else
        {
            indxt = indx[ir-1];
            aval = a[indxt];
            indx[ir-1] = indx[0];
            -- ir ;
            if ( ir == 1 )
            {
                indx[0] = indxt;
                break;
            }
        }

        i = l;
        j = l + l;

        while ( j <= ir )
        {
            if ( j < ir && a[indx[j-1]] < a[indx[j]] )
                ++ j;

            if ( aval < a[indx[j-1]] )
            {
                indx[i-1] = indx[j-1];
                i = j;
                j += j;
            }
            else
                j = ir + 1;
        }
        indx[i-1] = indxt;
    }

    return indx;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void  * _i4vec_sort_heap_index_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.
    Once the index array is computed, the sorting can be carried out
    "implicitly:
      a(indx(*))
    or explicitly, by the call
      i4vec_permute ( n, indx, a )
    after which a(*) is sorted.
    Note that the index vector is 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int A[N], an array to be index-sorted.
    Output, int I4VEC_SORT_HEAP_INDEX_D[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int aval;
    dim_typ i;
    int *indx;
    int indxt;
    int ir;
    dim_typ j, l;

    if ( n < 1 )
        return NULL;

    indx = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
        indx[i] = i;

    if ( n == 1 )
        return indx;

    l = n / 2 + 1;
    ir = n;

    for ( ; ; )
    {
        if ( 1 < l )
        {
            -- l;
            indxt = indx[l-1];
            aval = a[indxt];
        }
        else
        {
            indxt = indx[ir-1];
            aval = a[indxt];
            indx[ir-1] = indx[0];
            -- ir;

            if ( ir == 1 )
            {
                indx[0] = indxt;
                break;
            }
        }

        i = l;
        j = l + l;

        while ( j <= ir )
        {
            if ( j < ir && a[indx[j]] < a[indx[j-1]] )
                ++ j;


            if ( a[indx[j-1]] < aval )
            {
                indx[i-1] = indx[j-1];
                i = j;
                j +=j;
            }
            else
                    j = ir + 1;
        }
        indx[i-1] = indxt;
    }
    return indx;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sort_quick_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_QUICK_A ascending sorts an I4VEC using quick sort.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    Input:
      N = 7
      A = ( 6, 7, 3, 2, 9, 1, 8 )
    Output:
      A = ( 1, 2, 3, 6, 7, 8, 9 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries of A.
    Input/output, int A[N].  On input, the array to be sorted.
    On output, A has been reordered into ascending order.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    # define LEVEL_MAX 30

    int base;
    int l_segment;
    int level;
    dim_typ n_segment;
    int rsave[LEVEL_MAX];
    int r_segment;

    if ( n <= 1 )
        return NULL;

    level = 1;
    rsave[0] = n + 1;
    base = 1;
    n_segment = n;

    while ( 0 < n_segment )
    {
        /*
        Partition the segment.
        */
        i4vec_part_quick_a ( n_segment, a+base-1, &l_segment, &r_segment );
        /*
        If the left segment has more than one element, we need to partition it.
        */
        if ( 1 < l_segment )
        {
            if ( LEVEL_MAX < level )
                return NULL;
            ++ level;
            n_segment = l_segment;
            rsave[level-1] = r_segment + base - 1;
        }
        /*
        The left segment and the middle segment are sorted.
        Must the right segment be partitioned?
        */
        else if ( r_segment < n_segment )
        {
            n_segment += 1 - r_segment;
            base += r_segment - 1;
        }
        /*
        Otherwise, we back up a level if there is an earlier one.
        */
        else
        {
            for ( ; ; )
            {
                if ( 1 < level )
                {
                    base = rsave[level-1];
                    n_segment = rsave[level-2] - rsave[level-1];
                    -- level;
                    if ( 0 < n_segment )
                        break;
                }
                else
                {
                    n_segment = 0;
                    break;
                }
            }
        }
    }

    return NULL;
    # undef LEVEL_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sort_shell_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORT_SHELL_A ascending sorts an I4VEC using Shell's sort.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input/output, int A[N].
    On input, an array to be sorted.
    On output, the sorted array.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int asave;
    int base;
    dim_typ i;
    int ifree;
    int inc;
    int ipow;
    dim_typ j, k;
    int maxpow;

    if ( n <= 1 )
        return NULL;
    /*
    Determine the smallest MAXPOW so that
    N <= ( 3^MAXPOW - 1 ) / 2
    */
    maxpow = 1;
    base = 3;

    while ( base < (n<<1) + 1 )
    {
        ++ maxpow;
        base *= 3;
    }

    if ( 1 < maxpow )
    {
        -- maxpow;
        base /= 3;
    }
    /*
    Now sort groups of size ( 3^IPOW - 1 ) / 2.
    */
    for ( ipow = maxpow; 1 <= ipow; --ipow )
    {
        inc = ( base - 1 ) / 2;
        base = base / 3;
        /*
        Sort the values with indices equal to K mod INC.
        */
        for ( k = 1; k <= inc; ++k )
        {
            /*
            Insertion sort of the items with index
            INC+K, 2*INC+K, 3*INC+K, ...
            */
            for ( i = inc + k; i <= n; i = i + inc )
            {
                asave = a[i-1];
                ifree = i;
                j = i - inc;

                for ( ; ; )
                {
                    if ( j < 1 ||  a[j-1] <= asave )
                        break;

                    ifree = j;
                    a[j+inc-1] = a[j-1];
                    j -=  inc;
                }
                a[ifree-1] = asave;
            }
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sorted_undex ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORTED_UNDEX returns unique sorted indexes for a sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The goal of this routine is to determine a vector UNDX,
    which points, to the unique elements of X, in sorted order,
    and a vector XDNU, which identifies, for each entry of X, the index of
    the unique sorted element of X.
    This is all done with index vectors, so that the elements of
    X are never moved.
    Assuming X is already sorted, we examine the entries of X in order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.
    Once this process has been completed, the vector X could be
    replaced by a compressed vector XU, containing the unique entries
    of X in sorted order, using the formula
      XU(I) = X(UNDX(I)).
    We could then, if we wished, reconstruct the entire vector X, or
    any element of it, by index, as follows:
      X(I) = XU(XDNU(I)).
    We could then replace X by the combination of XU and XDNU.
    Later, when we need the I-th entry of X, we can locate it as
    the XDNU(I)-th entry of XU.
    Here is an example of a vector X, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.
      I    X    XU  Undx  Xdnu
    ----+----+----+-----+-----+
      0 | 11 |  11    0     0
      1 | 11 |  22    4     0
      2 | 11 |  33    7     0
      3 | 11 |  55    8     0
      4 | 22 |              1
      5 | 22 |              1
      6 | 22 |              1
      7 | 33 |              2
      8 | 55 |              3
    INDX(2) = 3 means that sorted item(2) is X(3).
    XDNI(2) = 5 means that X(2) is sorted item(5).
    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
    XU(XDNU(I))) = X(I).
    XU(I)        = X(UNDX(I)).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int X_NUM, the number of data values.
    Input, int X_VAL[X_NUM], the data values.
    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.
    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.
    Output, int XDNU[X_NUM], the XDNU vector.
*/
{
	const _2dt3pi * const s_data = data;
	
	const register dim_typ x_num = s_data->a0;
	const register dim_typ x_unique_num = s_data->a1;
	int * x_val = s_data->a2;
	int * undx = s_data->a3;
	int * xdnu = s_data->a4;
	
    dim_typ i, j;
    /*
    Walk through the sorted array.
    */
    i = j = 0;
    undx[j] = i;
    xdnu[i] = j;

    for ( i = 1; i < x_num; ++i)
    {
        if ( x_val[i] != x_val[undx[j]] )
        {
            ++ j;
            undx[j] = i;
        }
        xdnu[i] = j;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sorted_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements in A.
    Input/output, int A[N].  On input, the sorted
    integer array.  On output, the unique elements in A.
    Output, int I4VEC_SORTED_UNIQUE, the number of unique elements in A.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    dim_typ unique_num;

    if ( n <= 0 )
    {
    	result = 0;
        return &result;
    }

    unique_num = 1;

    for ( i = 1; i < n; ++i )
        if ( a[i] != a[unique_num-1] )
        {
            ++ unique_num;
            a[unique_num-1] = a[i];
        }

	result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sorted_unique_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    Because the array is sorted, this algorithm is O(N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, int A[N], the sorted array to examine.
    Output, int I4VEC_SORTED_UNIQUE_COUNT, the number of unique elements of A.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i;
    dim_typ unique_num;

    if ( n < 1 )
    {
    	result = 0;
        return &result;
    }

    unique_num = 1;

    for ( i = 1; i < n; ++i )
        if ( a[i-1] != a[i] )
            ++ unique_num;

    result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_sorted_unique_hist ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SORTED_UNIQUE_HIST histograms the unique elements of a sorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, int A[N], the array to examine, which must have been
    sorted.
    Input, int MAXUNIQ, the maximum number of unique elements
    that can be handled.  If there are more than MAXUNIQ unique
    elements in A, the excess will be ignored.
    Output, int *UNIQUE_NUM, the number of unique elements of A.
    Output, int AUNIQ[UNIQUE_NUM], the unique elements of A.
    Output, int ACOUNT[UNIQUE_NUM], the number of times each element
    of AUNIQ occurs in A.
*/
{
	const dtpidtpdt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1; 
	const register dim_typ maxuniq = s_data->a2;
	dim_typ * unique_num = s_data->a3;
	int * auniq = s_data->a4;
	int * acount = s_data->a5; 
	
    dim_typ i;
    int index;
    /*
    Start taking statistics.
    */
    index = -1;

    for ( i = 0; i < n; ++i )
    {
        if ( i == 0 )
        {
            index = 0;
            auniq[index] = a[0];
            acount[index] = 1;
        }
        else if ( a[i] == auniq[index] )
            ++ acount[index];
        else if ( index + 1 < maxuniq )
        {
            ++ index ;
            auniq[index] = a[i];
            acount[index] = 1;
        }
    }

    ++ *unique_num;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_split ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SPLIT "splits" an unsorted I4VEC based on a splitting value.
  Discussion:
    An I4VEC is a vector of I4's.
    If the vector is already sorted, it is simpler to do a binary search
    on the data than to call this routine.
    The vector is not assumed to be sorted before input, and is not
    sorted during processing.  If sorting is not needed, then it is
    more efficient to use this routine.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input/output, int A[N], the array to split.  On output,
    all the entries of A that are less than or equal to SPLIT
    are in A(1:I4VEC_SPLIT).
    Input, int SPLIT, the value used to split the vector.
    It is not necessary that any value of A actually equal SPLIT.
    Output, int I4VEC_SPLIT, indicates the position of the last
    entry of the split vector that is less than or equal to SPLIT.
*/
{
	static int result = INT_MAX;
	
	const dtpii * const s_data = data; 
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int split = s_data->a2; 	
	
    dim_typ i;
    dim_typ i1;
    dim_typ i2;
    dim_typ i3;
    dim_typ j1;
    dim_typ j2;
    dim_typ j3;
    int temp;
    /*
    Partition the vector into A1, A2, A3, where
    A1 = A(I1:J1) holds values <= SPLIT,
    A2 = A(I2:J2) holds untested values,
    A3 = A(I3:J3) holds values > SPLIT.
    */
    i1 = i2 = 1;
    j1 = 0;
    j2 = j3 = n;
    i3 = n + 1;
    /*
    Pick the next item from A2, and move it into A1 or A3.
    Adjust indices appropriately.
    */
    for ( i = 0; i < n; ++i)
        if ( a[i2-1] <= split )
        {
            ++ i2;
            ++ j1;
        }
        else
        {
            temp    = a[i2-1];
            a[i2-1] = a[i3-2];
            a[i3-2] = temp;
            -- i3;
            -- j2;
        }
        
    result = j1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_std ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_STD returns the standard deviation of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 August 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int X[N], the vector whose variance is desired.
    Output, double I4VEC_STD, the standard deviation of the vector entries.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	
    dim_typ i;
    ityp mean;
    ityp std;

    if ( n < 2 )
        std = 0.00;
    else
    {
        mean = i4vec_mean ( n, x );
        std = 0.00;
        for ( i = 0; i < n; ++i )
            std += pow ( ( ityp ) x[i] - mean, 2 );

        std = sqrt ( std / ( ityp ) ( n - 1 ) );
    }

	result = std;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_swap ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_SWAP swaps two I4VEC's.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the arrays.
    Input/output, int A1[N], A2[N], the
    two arrays whose entries are to be swapped.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a1 = s_data->a1;
	int * a2 = s_data->a2;
	
    dim_typ i;
    int j;

    for ( i = 0; i < n; ++i )
    {
        j     = a1[i];
        a1[i] = a2[i];
        a2[i] = j;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_undex ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_UNDEX returns unique sorted indexes for an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The goal of this routine is to determine a vector UNDX,
    which points, to the unique elements of X, in sorted order,
    and a vector XDNU, which identifies, for each entry of X, the index of
    the unique sorted element of X.
    This is all done with index vectors, so that the elements of
    X are never moved.
    The first step of the algorithm requires the indexed sorting
    of X, which creates arrays INDX and XDNI. (If all the entries
    of X are unique, then these arrays are the same as UNDX and XDNU.)
    We then use INDX to examine the entries of X in sorted order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.
    Once this process has been completed, the vector X could be
    replaced by a compressed vector XU, containing the unique entries
    of X in sorted order, using the formula
      XU(*) = X(UNDX(*)).
    We could then, if we wished, reconstruct the entire vector X, or
    any element of it, by index, as follows:
      X(I) = XU(XDNU(I)).
    We could then replace X by the combination of XU and XDNU.
    Later, when we need the I-th entry of X, we can locate it as
    the XDNU(I)-th entry of XU.
    Here is an example of a vector X, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.
      I    X  Indx  Xdni      XU  Undx  Xdnu
    ----+----+-----+-----+-------+-----+-----+
      0 | 11     0     0 |    11     0     0
      1 | 22     2     4 |    22     1     1
      2 | 11     5     1 |    33     3     0
      3 | 33     8     7 |    55     4     2
      4 | 55     1     8 |                 3
      5 | 11     6     2 |                 0
      6 | 22     7     5 |                 1
      7 | 22     3     6 |                 1
      8 | 11     4     3 |                 0
    INDX(2) = 3 means that sorted item(2) is X(3).
    XDNI(2) = 5 means that X(2) is sorted item(5).
    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
    XU(XDNU(I))) = X(I).
    XU(I)        = X(UNDX(I)).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int X_NUM, the number of data values.
    Input, int X_VAL[X_NUM], the data values.
    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.
    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.
    Output, int XDNU[X_NUM], the XDNU vector.
*/
{
	const _2dt3pi * const s_data = data;
	
	const register dim_typ x_num = s_data->a0;
	const register dim_typ x_unique_num = s_data->a1;
	int * x_val = s_data->a2;
	int * undx = s_data->a3;
	int * xdnu = s_data->a4;
	
    dim_typ i;
    int *indx;
    dim_typ j;
    /*
    Implicitly sort the array.
    */
    indx = i4vec_sort_heap_index_a ( x_num, x_val );
    /*
    Walk through the implicitly sorted array.
    */
    i = j = 0;
    undx[j] = indx[i];
    xdnu[indx[i]] = j;

    for ( i = 1; i < x_num; ++i )
    {
        if ( x_val[indx[i]] != x_val[undx[j]] )
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
__MATHSUITE __JBURKARDT  void *   _i4vec_unique_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_UNIQUE_COUNT counts the unique elements in an unsorted I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    Because the array is unsorted, this algorithm is O(N^2).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, int A[N], the array to examine, which does NOT have to
    be sorted.
    Output, int I4VEC_UNIQUE_COUNT, the number of unique elements of A.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;
    dim_typ unique_num = 0;

    for ( i = 0; i < n; ++i )
    {
        ++ unique_num;

        for ( j = 0; j < i; ++j )
            if ( a[i] == a[j] )
            {
                -- unique_num;
                break;
            }
    }

	result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_unique_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_UNIQUE_INDEX indexes the unique occurrence of values in an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    For element A(I) of the vector, UNIQUE_INDEX(I) is the uniqueness index
    of A(I).  That is, if A_UNIQUE contains the unique elements of A,
    gathered in order, then
      A_UNIQUE ( UNIQUE_INDEX(I) ) = A(I)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, int A[N], the array.
    Output, int I4VEC_UNIQUE_INDEX[N], the unique index.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    dim_typ i, j;
    int *unique_index;
    int unique_num;

    unique_index = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i)
        unique_index[i] = -1;
    unique_num = 0;

    for ( i = 0; i < n; ++i )
        if ( unique_index[i] == -1 )
        {
            unique_index[i] = unique_num;
            for ( j = i + 1; j < n; ++j )
                if ( a[i] == a[j] )
                    unique_index[j] = unique_num;
            ++ unique_num;
        }
    return unique_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_value_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_VALUE_INDEX indexes I4VEC entries equal to a given value.
  Discussion:
    An I4VEC is a vector of I4's.
  Example:
    Input:
      N = 10
      A = (  2, 3, 1, 3, 2, 4, 2, 3, 5, 3 )
      X_VALUE = 3
    Output:
      N_INDEX = 4
      VALUE_INDEX = ( 2, 4, 8, 10 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input, int A[N], the array to be indexed.
    Input, int VALUE, a value to be searched for.
    Input, int MAX_INDEX, the maximum number of indices to find.
    Output, int *N_INDEX, the number of entries equal to VALUE.
    Output, int I4VEC_VALUE_INDEX[MAX_INDEX], the indices of entries
    equal to VALUE.
*/
{
	const _2dtpiipi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ max_index = s_data->a1;
	int * a = s_data->a2;
	int value = s_data->a3;
	int * n_index = s_data->a4;
	
    dim_typ i;
    int * value_index = ( int * ) malloc ( max_index * sizeof ( int ) );
    *n_index = 0;

    for ( i = 1; i <= n; ++i )
        if ( a[i-1] == value )
        {
            if ( max_index <= *n_index )
                break;
            value_index[*n_index] = i;
            ++ *n_index;
        }

    return value_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_value_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_VALUE_NUM counts I4VEC entries equal to a given value.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input, int A[N], the array to be indexed.
    Input, int VALUE, a value to be searched for.
    Input, int I4VEC_VALUE_NUM, the number of times the value occurs.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpii * const s_data = data; 
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	int value = s_data->a2; 
	
	dim_typ value_num;
    for (dim_typ i = value_num = 0; i < n; ++i)
        if ( a[i] == value )
            ++ value_num;

	result = value_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_variance ( void * data)
/*****************************************************************************/
/*
  Purpose:
    I4VEC_VARIANCE returns the variance of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 May 1999
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int X[N], the vector whose variance is desired.
    Output, double I4VEC_VARIANCE, the variance of the vector entries.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * x = s_data->a1;
	
    dim_typ i;
    ityp mean;
    ityp variance;

    if ( n < 2 )
        variance = 0.00;
    else
    {
        mean = i4vec_mean ( n, x );
        variance = 0.0;
        for ( i = 0; i < n; ++i )
            variance += pow ( ( ityp ) x[i] - mean, 2 );
        variance /= ( ityp ) ( n - 1 );
    }

	result = variance;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec_width ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_WIDTH returns the "width" of an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
    The width of an I4VEC is simply the maximum of the widths of
    its entries.
    The width of a single integer is the number of characters
    necessary to print it.
    The width of an I4VEC can be useful when the vector is
    to be printed.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 September 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, int A[N], the vector.
    Output, int I4VEC_WIDTH, the width of the vector.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    int width = -1;
    for (dim_typ i = 0; i < n; ++i )
        width = MAX ( width, i4_width ( a[i] ) );
        
    result = width;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _i4vec_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ZEROS zeroes an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2008
  Author:
    John Burkardt
    
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, int A[N], a vector of zeroes.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        a[i] = 0;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _i4vec_zeros_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ZEROS_NEW creates and zeroes an I4VEC.
  Discussion:
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, int I4VEC_ZEROS_NEW[N], a vector of zeroes.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    int * a = ( int * ) malloc ( n * sizeof ( int ) );
    for (dim_typ i = 0; i < n; ++i )
        a[i] = 0;
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec2_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC2_COMPARE compares pairs of integers stored in two vectors.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data items.
    Input, int A1[N], A2[N], contain the two components of each item.
    Input, int I, J, the items to be compared.  These values will be
    1-based indices for the arrays A1 and A2.
    Output, int I4VEC2_COMPARE, the results of the comparison:
    -1, item I < item J,
     0, item I = item J,
    +1, item J < item I.
*/
{
	static short result = SHRT_MAX;
	
	const dt2i2pi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	int i = s_data->a1;
	int j = s_data->a2;
	int * a1 = s_data->a3;
	int * a2 = s_data->a4;
	
	
    short isgn = 0;

    if ( a1[i-1] < a1[j-1] )
    isgn = -1;
    else if ( a1[i-1] == a1[j-1] )
    {
        if ( a2[i-1] < a2[j-1] )
            isgn = -1;
        else if ( a2[i-1] < a2[j-1] )
            isgn = 0;
        else if ( a2[j-1] < a2[i-1] )
            isgn = +1;
    }
    else if ( a1[j-1] < a1[i-1] )
        isgn = +1;
        
    result = isgn;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec2_sort_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC2_SORT_A ascending sorts an I4VEC2.
  Discussion:
    Each item to be sorted is a pair of integers (I,J), with the I
    and J values stored in separate vectors A1 and A2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items of data.
    Input/output, int A1[N], A2[N], the data to be sorted..
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a1 = s_data->a1;
	int * a2 = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    int temp;
    /*
    Initialize.
    */
    i = indx = isgn = j = 0;
    /*
    Call the external heap sorter.
    */
    for ( ; ; )
    {
        sort_heap_external ( n, &indx, &i, &j, isgn );
        /*
        Interchange the I and J objects.
        */
        if ( 0 < indx )
        {
            temp    = a1[i-1];
            a1[i-1] = a1[j-1];
            a1[j-1] = temp;

            temp    = a2[i-1];
            a2[i-1] = a2[j-1];
            a2[j-1] = temp;
        }
        /*
        Compare the I and J objects.
        */
        else if ( indx < 0 )
            isgn = i4vec2_compare ( n, a1, a2, i, j );
        else if ( indx == 0 )
            break;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec2_sort_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC2_SORT_D descending sorts an I4VEC2.
  Discussion:
    Each item to be sorted is a pair of integers (I,J), with the I
    and J values stored in separate vectors A1 and A2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items of data.
    Input/output, int A1[N], A2[N], the data to be sorted..
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a1 = s_data->a1;
	int * a2 = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    int temp;
    /*
    Initialize.
    */
    i = indx = isgn = j = 0;
    /*
    Call the external heap sorter.
    */
    for ( ; ; )
    {
        sort_heap_external ( n, &indx, &i, &j, isgn );
        /*
        Interchange the I and J objects.
        */
        if ( 0 < indx )
        {
            temp    = a1[i-1];
            a1[i-1] = a1[j-1];
            a1[j-1] = temp;

            temp    = a2[i-1];
            a2[i-1] = a2[j-1];
            a2[j-1] = temp;
        }
        /*
        Compare the I and J objects.
        */
        else if ( indx < 0 )
            isgn = -i4vec2_compare ( n, a1, a2, i, j );
        else if ( indx == 0 )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec2_sorted_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC2_SORTED_UNIQUE keeps the unique elements in an I4VEC2.
  Discussion:
    Item I is stored as the pair A1(I), A2(I).
    The items must have been sorted, or at least it must be the
    case that equal items are stored in adjacent vector locations.
    If the items were not sorted, then this routine will only
    replace a string of equal values by a single representative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items.
    Input/output, int A1[N], A2[N].
    On input, the array of N items.
    On output, an array of UNIQUE_NUM unique items.
    Output, int *UNIQUE_NUM, the number of unique items.
*/
{
	const dt3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a1 = s_data->a1;
	int * a2 = s_data->a2;
	int * unique_num = s_data->a3;
	
    int itest;
    *unique_num = 0;

    if ( n <= 0 )
        return NULL;

    *unique_num = 1;

    for ( itest = 1; itest < n; ++itest )
        if ( a1[itest] != a1[*unique_num-1] ||a2[itest] != a2[*unique_num-1] )
        {
            a1[*unique_num] = a1[itest];
            a2[*unique_num] = a2[itest];
            ++ *unique_num;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _i4vec2_sorted_unique_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC2_SORTED_UNIQUE_COUNT counts unique elements in an I4VEC2.
  Discussion:
    Item I is stored as the pair A1(I), A2(I).
    The items must have been sorted, or at least it must be the
    case that equal items are stored in adjacent vector locations.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 July 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items.
    Input, int A1[N], A2[N], the array of N items.
    Output, int I4VEC2_SORTED_UNIQUE_NUM, the number of unique items.
*/
{
	static int result = INT_MAX;
	
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * a1 = s_data->a1;
	int * a2 = s_data->a2;
	
    dim_typ i;
    dim_typ iu;
    int unique_num;

    if ( n <= 0 )
    {
    	result = 0;
        return &result;
    }

    iu = 0;
    unique_num = 1;

    for ( i = 1; i < n; ++i )
        if ( a1[i] != a1[iu] ||a2[i] != a2[iu] )
        {
            iu = i;
            ++ unique_num;
        }

	result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l4_to_i4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    L4_TO_I4 converts an L4 to an I4.
  Discussion:
    0 is FALSE, and anything else if TRUE.
    An I4 is an integer value.
    An L4 is a logical value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 January 2012
  Author:
    John Burkardt
  Parameters:
    Input, int L, a logical value.
    Output, int L4_TO_I4, the integer value of L.
*/
{
	static bool result = 2;
	
	const register int l = *(int *) data;
	
	result = l != 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _l4_xor ( void * data)
/*****************************************************************************/
/*
  Purpose:
    L4_XOR returns the exclusive OR of two L4's.
  Discussion:
    An L4 is a logical value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 May 2014
  Author:
   John Burkardt
  Parameters:
    Input, int L1, L2, two values whose exclusive OR is needed.
    Output, int L4_XOR, the exclusive OR of L1 and L2.
*/
{
	static bool result = 2;
	
	int * const a_data = data;
	const register int l1 = a_data[0];
	const register int l2 = a_data[1];
	
	result = ( (     l1   && ( ! l2 ) ) || ( ! l1 ) &&     l2    );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pascal_to_i4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    PASCAL_TO_I4 converts Pacal triangle coordinates to a linear index.
  Discussion:
    We describe the grid points in a Pascal triangle in two ways:
    As a linear index K:
                     1
                   2   3
                 4   5   6
               7   8   9   10
    As elements (I,J) of Pascal's triangle:
                     0,0
                  1,0   0,1
               2,0   1,1    0,2
            3,0   2,1   1,2    0,3
  Example:
     K  I  J
     1  0  0
     2  1  0
     3  0  1
     4  2  0
     5  1  1
     6  0  2
     7  3  0
     8  2  1
     9  1  2
    10  0  3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 April 2015
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the row and column indices.  I and J
    must be nonnegative.
    Output, int PASCAL_TO_I4, the linear index of the (I,J) element.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register int i = a_data[0];
	const register int j = a_data[1];
	
    const register int d = i + j;
    
    result = (i < 0 || j < 0) ? INT_MAX : ( d * ( d + 1 ) ) / 2 + j + 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm_cycle ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_CYCLE analyzes a permutation.
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
	const dt3pii * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	int * isgn = s_data->a2;
	int * ncycle = s_data->a3;
	int iopt = s_data->a4;
	
    dim_typ i;
    int i1;
    int i2;
    int ierror;
    int is;

    ierror = perm0_check ( n, p );

    if ( ierror != 0 )
        return NULL;

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
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm0_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM0_CHECK checks a permutation of (0,...,N-1)
  Discussion:
    The routine verifies that each of the integers from 0 to
    to N-1 occurs among the N entries of the permutation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int P[N], the array to check.
    Output, int PERM0_CHECK:
    0, P is a legal permutation of (0,...,N-1).
    1, P is not a legal permutation of (0,...,N-1);
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
    int ierror = 0;
    int location;
    int value;

    for ( value = 0; value < n; ++value )
    {
        ierror = 1;

        for ( location = 0; location < n; ++location)
            if ( p[location] == value )
            {
                ierror = 0;
                break;
            }

        if ( ierror != 0 )
        {
        	result = INT_MAX;
            return &result;
        }

    }

	result = ierror;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _perm0_uniform_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM0_UNIFORM_NEW selects a random permutation of 0,...,N-1.
  Discussion;
    The algorithm is known as the Fisher-Yates or Knuth shuffle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 May 2015
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of objects to be permuted.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int PERM0_UNIFORM_NEW[N], a permutation of ( 0, 1, ..., N-1).
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    dim_typ i;
    int j, k;
    int * p = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
        p[i] = i;

    for ( i = 0; i < n - 1; ++i )
    {
        j = i4_uniform_ab ( i, n - 1, seed );
        k    = p[i];
        p[i] = p[j];
        p[j] = k;
    }

    return p;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _perm1_check ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM1_CHECK checks a permutation of (1,...,N).
  Discussion:
    The routine verifies that each of the integers from 1 to
    to N occurs among the N entries of the permutation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 May 2015
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, int P[N], the array to check.
    Output, int PERM1_CHECK:
    0, P is a legal permutation of (1,...,N).
    1, P is not a legal permutation of (1,...,N);
*/
{
	static int result = INT_MAX;
	
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
    int ierror = 0;
    int location;
    int value;

    for ( value = 1; value <= n; ++value)
    {
        ierror = 1;
        for ( location = 0; location < n; ++location )
            if ( p[location] == value )
            {
                ierror = 0;
                break;
            }

        if ( ierror != 0 )
            break;
    }

	result = ierror;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _perm1_uniform_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM1_UNIFORM_NEW selects a random permutation of 1,...,N.
  Discussion;
    The algorithm is known as the Fisher-Yates or Knuth shuffle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 May 2015
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of objects to be permuted.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int PERM1_UNIFORM_NEW[N], a permutation of ( 1, ..., N).
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    dim_typ i;
    int j, k;
    int * p = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i)
        p[i] = i + 1;

    for ( i = 0; i < n - 1; ++i )
    {
        j = i4_uniform_ab ( i, n - 1, seed );
        k    = p[i];
        p[i] = p[j];
        p[j] = k;
    }

    return p;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _prime ( void * data)
/******************************************************************************/
/*
  Purpose:
    PRIME returns any of the first PRIME_MAX prime numbers.
  Discussion:
    PRIME_MAX is 1600, and the largest prime stored is 13499.
    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 March 2009
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    US Department of Commerce, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Daniel Zwillinger,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996, pages 95-98.
  Parameters:
    Input, int N, the index of the desired prime number.
    In general, is should be true that 0 <= N <= PRIME_MAX.
    N = -1 returns PRIME_MAX, the index of the largest prime available.
    N = 0 is legal, returning PRIME = 1.
    Output, int PRIME, the N-th prime.  If N is out of range, PRIME
    is returned as -1.
*/
{
	static int result = INT_MAX;
	
	const register dim_typ n = *(dim_typ *) data;
	
    # define PRIME_MAX 1600

    int npvec[PRIME_MAX] =
    {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
        31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
        73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
        127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
        179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
        233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
        283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
        353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
        419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
        467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
        547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
        607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
        661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
        739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
        811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
        877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
        947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
        1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
        1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
        1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
        1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
        1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
        1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
        1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
        1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
        1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
        1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
        1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
        1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
        1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
        1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
        2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
        2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
        2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
        2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
        2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
        2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
        2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
        2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
        2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
        2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
        2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
        2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
        3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
        3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
        3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
        3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
        3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
        3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511,
        3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
        3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
        3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
        3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
        3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
        3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
        4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
        4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
        4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
        4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
        4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
        4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493,
        4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
        4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
        4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
        4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
        4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
        4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
        5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
        5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
        5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
        5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
        5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
        5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
        5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
        5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
        5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
        5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
        5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
        5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
        6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
        6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
        6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
        6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
        6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
        6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571,
        6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
        6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
        6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
        6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
        6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
        7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
        7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
        7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
        7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
        7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
        7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
        7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
        7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
        7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
        7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
        7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
        8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
        8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
        8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
        8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
        8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
        8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
        8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
        8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741,
        8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
        8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
        8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
        9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
        9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
        9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
        9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
        9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
        9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
        9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
        9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
        9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
        9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
        9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
        10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
        10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
        10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
        10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
        10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
        10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
        10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
        10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
        10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
        10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
        10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
        11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
        11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
        11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
        11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
        11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
        11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
        11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
        11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
        11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
        11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
        12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
        12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
        12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
        12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
        12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
        12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
        12569,12577,12583,12589,12601,12611,12613,12619,12637,12641,
        12647,12653,12659,12671,12689,12697,12703,12713,12721,12739,
        12743,12757,12763,12781,12791,12799,12809,12821,12823,12829,
        12841,12853,12889,12893,12899,12907,12911,12917,12919,12923,
        12941,12953,12959,12967,12973,12979,12983,13001,13003,13007,
        13009,13033,13037,13043,13049,13063,13093,13099,13103,13109,
        13121,13127,13147,13151,13159,13163,13171,13177,13183,13187,
        13217,13219,13229,13241,13249,13259,13267,13291,13297,13309,
        13313,13327,13331,13337,13339,13367,13381,13397,13399,13411,
        13417,13421,13441,13451,13457,13463,13469,13477,13487,13499
    };

    if ( n == -1 )
    {
    	result = PRIME_MAX;
        return &result;
    }
    else if ( n == 0 )
    {
    	result = 1;
        return &result;
    }
    else if ( n <= PRIME_MAX )
    {
    	result = npvec[n-1];
        return &result;
    }

	result = 0;
    return &result;
    # undef PRIME_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _sort_heap_external ( void * data)
/******************************************************************************/
/*
  Purpose:
    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
  Discussion:
    The actual list is not passed to the routine.  Hence it may
    consist of integers, reals, numbers, names, etc.  The user,
    after each return from the routine, will be asked to compare or
    interchange two items.
    The current version of this code mimics the FORTRAN version,
    so the values of I and J, in particular, are FORTRAN indices.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 February 2004
  Author:
    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the length of the input list.
    Input/output, int *INDX.
    The user must set INDX to 0 before the first call.
    On return,
      if INDX is greater than 0, the user must interchange
      items I and J and recall the routine.
      If INDX is less than 0, the user is to compare items I
      and J and return in ISGN a negative value if I is to
      precede J, and a positive value otherwise.
      If INDX is 0, the sorting is done.
    Output, int *I, *J.  On return with INDX positive,
    elements I and J of the user's list should be
    interchanged.  On return with INDX negative, elements I
    and J are to be compared by the user.
    Input, int ISGN. On return with INDX negative, the
    user should compare elements I and J of the list.  If
    item I is to precede item J, set ISGN negative,
    otherwise set ISGN positive.
*/
{
	const dt3pii * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * indx = s_data->a1;
	int * i = s_data->a2;
	int * j = s_data->a3;
	int isgn = s_data->a4;
	
    static int i_save = 0;
    static int j_save = 0;
    static int k = 0;
    static int k1 = 0;
    static int n1 = 0;
    /*
    INDX = 0: This is the first call.
    */
    if ( *indx == 0 )
    {

        i_save = j_save = 0;
        k = n / 2;
        k1 = k;
        n1 = n;
    }
    /*
    INDX < 0: The user is returning the results of a comparison.
    */
    else if ( *indx < 0 )
    {
        if ( *indx == -2 )
        {
            if ( isgn < 0 )
                ++ i_save;
            j_save = k1;
            k1 = i_save;
            *indx = -1;
            *i = i_save;
            *j = j_save;
            return NULL;
        }

        if ( 0 < isgn )
        {
            *indx = 2;
            *i = i_save;
            *j = j_save;
            return NULL;
        }

        if ( k <= 1 )
        {
            if ( n1 == 1 )
                i_save = j_save = *indx = 0;
            else
            {
                i_save = n1;
                j_save = 1;
                -- n1;
                *indx = 1;
            }
            *i = i_save;
            *j = j_save;
            return NULL;
        }
        -- k;
        k1 = k;
    }
    /*
    0 < INDX: the user was asked to make an interchange.
    */
    else if ( *indx == 1 )
        k1 = k;

    for ( ; ; )
    {

        i_save = k1<<1;

        if ( i_save == n1 )
        {
            j_save = k1;
            k1 = i_save;
            *indx = -1;
            *i = i_save;
            *j = j_save;
            return NULL;
        }
        else if ( i_save <= n1 )
        {
            j_save = i_save + 1;
            *indx = -2;
            *i = i_save;
            *j = j_save;
            return NULL;
        }

        if ( k <= 1 )
            break;

        -- k;
        k1 = k;
    }

    if ( n1 == 1 )
    {
        i_save = j_save = *indx = 0;
        *i = i_save;
        *j = j_save;
    }
    else
    {
        i_save = n1;
        j_save = *indx = 1;
        -- n1;
        *i = i_save;
        *j = j_save;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_to_i4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_TO_I4 converts a triangular coordinate to an integer.
  Discussion:
    Triangular coordinates are handy when storing a naturally triangular
    array (such as the lower half of a matrix) in a linear array.
    Thus, for example, we might consider storing
 (0,0)
 (1,0) (1,1)
 (2,0) (2,1) (2,2)
 (3,0) (3,1) (3,2) (3,3)
    as the linear array
 (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) (3,0) (3,1) (3,2) (3,3)
    Here, the quantities in parenthesis represent the natural row and
    column indices of a single number when stored in a rectangular array.
    Thus, our goal is, given the row I and column J of the data,
    to produce the value K which indicates its position in the linear
    array.
    The triangular numbers are the indices associated with the
    diagonal elements of the original array, T(0,0), T(1,1), T(2,2),
    T(3,3) and so on.
    The formula is:
      K = J + ( ( I * ( I + 1 ) ) / 2
  Example:
    I  J  K
    0  0  0
    1  0  1
    1  1  2
    2  0  3
    2  1  4
    2  2  5
    3  0  6
    3  1  7
    3  2  8
    3  3  9
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 January 2009
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the row and column indices.  I and J must
    be nonnegative, and J must not be greater than I.
    Output, int TRIANGLE_TO_I4, the linear index of the (I,J) element.
*/
{
	static int result = INT_MAX;
	
	int * const a_data = data;
	const register dim_typ i = a_data[0];
	const register dim_typ j = a_data[1];
	
	result = i < j ? INT_MAX : j + ((  i * ( i + 1 ) ) >> 1);
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4_to_bvec ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4_TO_BVEC converts an integer into a binary vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int I4, the integer.
    Input, int N, the dimension of the vector.
    Output, int BVEC[N], the vector of binary remainders.
*/
{
	const idtpb * const s_data = data;
	int i4 = s_data->a0;
	const register dim_typ n = s_data->a1;
	bool * bvec = s_data->a2; 
	
    for (dim_typ i = n - 1; 0 <= i; --i )
    {
        bvec[i] = i4 % 2;
        i4 >>= 1;
    }
    return NULL;
}

#endif
