#ifndef __DISABLEDEEP_R8LIB

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_factorial2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R4_FACTORIAL2 computes the float factorial function.
  Discussion:
    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 ) (N even)
                    = Product ( N * (N-2) * (N-4) * ... * 1 ) (N odd)
  Example:
     N    Factorial2(N)
     0     1
     1     1
     2     2
     3     3
     4     8
     5    15
     6    48
     7   105
     8   384
     9   945
    10  3840
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the argument of the float factorial
    function.  If N is less than 1, R4_FACTORIAL2 is returned as 1.0.
    Output, float R4_FACTORIAL2, the value of Factorial2(N).
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
	uint64_t n_copy;
	ityp value;
	
	if ( n == 0 )
	{
		result = 1.00;
		return &result;
	}
	
	n_copy = n;
	
	while ( 1 < n_copy )
	{
		value *= ( ityp ) n_copy;
		n_copy -= 2;
	}
	
	result = value; 
	return &result;
}

/*****************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_fractional ( void * data)
/*****************************************************************************/
/*
  Purpose:
    r8_FRACTIONAL returns the fractional part of an r8.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp X, the argument.
    Output, ityp r8_FRACTIONAL, the fraction part of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	result = abs ( x ) - ( ityp ) ( ( int ) abs ( x ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_in_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_IN_01 is 1 if an r8 is in the range [0,1].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp A, the value.
    Output, int r8_IN_01, is 1 if A is between 0 and 1.
*/
{
	static bool result = 2;
	
	const register ityp a = *(ityp *) data;
	
	result = a>= 0.00 && a <= 1.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_is_int ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_IS_INT is 1 if an r8 represents an integer value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp R, the number to be checked.
    Output, int r8_IS_INT, is 1 if R is an integer value.
*/
{
	static bool result = 2;
	
	const register ityp r = *(ityp *) data;
	
    if ( ( ityp ) ( 2147483647 ) < r || r < - ( ityp ) ( 2147483647  ))
    {
    	result = false; 
        return &result;
    }
    else if ( r == ( ityp ) ( ( int ) ( r ) ) )
    {
    	result = true;
        return &result;
    }
        
    result = false; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_mant ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_MANT computes the "mantissa" or "fraction part" of an r8.
  Formula:
    X = S * R * 2**L
    S is +1 or -1,
    R is a real between 1.0 and 2.0,
    L is an integer.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp X, the real number to be decomposed.
    Output, int *S, the "sign" of the number.
    S will be -1 if X is less than 0, and +1 if X is greater
    than or equal to zero.
    Output, ityp *R, the mantissa of X.  R will be greater
    than or equal to 1, and strictly less than 2.  The one
    exception occurs if X is zero, in which case R will also
    be zero.
    Output, int *L, the integer part of the logarithm (base 2) of X.
*/
{
	const itpipitpi * const s_data = data;
	ityp x = s_data->a0;
	int * s = s_data->a1;
	ityp * r = s_data->a2;
	int * l = s_data->a3;
	
    /*
    Determine the sign.
    */
    *s = -1 + ((x>=0.00)<<1);
    /*
    Set R to the absolute value of X, and L to zero.
    Then force R to lie between 1 and 2.
    */
    *r = x<0.00 ? -x:x;

    *l = 0;
    /*
    Time to bail out if X is zero.
    */
    if ( x == 0.00 )
        return NULL;

    while ( 2.00 <= *r )
    {
        *r /= 2.00;
        ++ *l;
    }

    while ( *r < 1.00 )
    {
        *r *= 2.00;
        -- *l;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_mop ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_MOP returns the I-th power of -1 as an r8 value.
  Discussion:
    An r8 is a ityp value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 July 2008
  Author:
    John Burkardt
  Parameters:
    Input, int I, the power of -1.
    Output, ityp r8_MOP, the I-th power of -1.
*/
{
	static ityp result = MAX_VAL;
	
	const register int i = *(int *) data;
	
	result = -1.00 + ((( i % 2 ) == 0)<<1);
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_nint ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_NINT returns the nearest integer to an r8.
  Example:
        X         r8_NINT
      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 May 2006
  Author:
    John Burkardt
  Parameters:
    Input, ityp X, the value.
    Output, int r8_NINT, the nearest integer to X.
*/
{
	static int result = INT_MAX;
	
	const register ityp x = *(ityp *) data;
	
	result = (1 - ((x<0.00)<<1)) * ( int ) ( fabs ( x ) + 0.50 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_normal_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_NORMAL_01 returns a unit pseudonormal r8.
  Discussion:
    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.
    The Box-Muller method is used, which is efficient, but
    generates two values at a time.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 June 2010
  Author:
    John Burkardt
  Parameters:
    Input/output, int *SEED, a seed for the random number generator.
    Output, ityp r8_NORMAL_01, a normally distributed random value.
*/
{
	static ityp result = MAX_VAL;
	
	int * seed = data;
	
    ityp r1;
    ityp r2;
    static int seed1 = 0;
    static int seed2 = 0;
    static int seed3 = 0;
    static int used = -1;
    ityp x;
    static ityp y = 0.00;

    if ( used == -1 )
        used = 0;
    /*
    If we've used an even number of values so far, generate two more, return one,
    and save one.
    */
    if ( ( used % 2 )== 0 )
    {
        seed1 = *seed;
        r1 = r8_uniform_01 ( seed );

        if ( r1 == 0.00 )
        {
        	result = MAX_VAL;
            return &result;
        }

        seed2 = *seed;
        r2 = r8_uniform_01 ( seed );
        seed3 = *seed;
        *seed = seed2;

        x = sqrt ( - 2.00 * log ( r1 ) ) * cos ( M_2TPI * r2 );
        y = sqrt ( - 2.00 * log ( r1 ) ) * sin ( M_2TPI * r2 );
    }
    /*
    Otherwise, return the second, saved, value and the corresponding
    value of SEED.
    */
    else
    {
        x = y;
        *seed = seed3;
    }

    ++ used;

	result = x;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_normal_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_NORMAL_AB returns a scaled pseudonormal r8.
  Discussion:
    The normal probability distribution function (PDF) is sampled,
    with mean A and standard deviation B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, ityp A, the mean of the PDF.
    Input, ityp B, the standard deviation of the PDF.
    Input/output, int *SEED, a seed for the random number generator.
    Output, ityp r8_NORMAL_AB, a sample of the normal PDF.
*/
{
	static ityp result = MAX_VAL;
	
	const _2itpi * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	int * seed = s_data->a2;
	
	result = a + b * r8_normal_01 ( seed );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_power_fast ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_POWER_FAST computes the P-th power of R, for real R and integer P.
  Discussion:
    Obviously, R^P can be computed using P-1 multiplications.
    However, R^P can also be computed using at most 2*LOG2(P) multiplications.
    To do the calculation this way, let N = LOG2(P).
    Compute A, A^2, A^4, ..., A^N by N-1 successive squarings.
    Start the value of R^P at A, and each time that there is a 1 in
    the binary expansion of P, multiply by the current result of the squarings.
    This algorithm is not optimal.  For small exponents, and for special
    cases, the result can be computed even more quickly.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 January 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp R, the base.
    Input, int P, the power, which may be negative.
    Output, int *MULTS, the number of multiplications and divisions.
    Output, ityp r8_POWER_FAST, the value of R^P.
*/
{
	static ityp result = MAX_VAL;
	
	const itipi * const s_data = data;
	const register ityp r = s_data->a0;
	const register int p = s_data->a1;
	int * mults = s_data->a2;
	
    int p_mag;
    int p_sign;
    ityp r2;
    ityp value;

    *mults = 0;
    /*
    Special cases.
    */
    if ( r == 1.0 )
    {
    	result = 1.00;
        return &result;
    }

    if ( r == -1.0 )
    {
    	result = 1.00 - ((( p % 2 ) == 1)<<1);
        return &result;
    }

    if ( r == 0.00 )
    {
    	result = p <= 0 ? MAX_VAL : 0.00;
        return &result;
    }
    /*
    Special powers.
    */
    if ( p == -1 )
    {
        ++ *mults;
        result = 1.00 / r;
        return &result;
    }
    else if ( p == 0 )
    {
    	result = 1.00;
        return &result;
    }
    else if ( p == 1 )
    {
    	result = r;
        return &result;
    }
    /*
    Some work to do.
    */
    p_mag = abs ( p );
    p_sign = i4_sign ( p );

    value = 1.00;
    r2 = r;

    while ( 0 < p_mag )
    {
        if ( ( p_mag % 2 ) == 1 )
        {
            value *= r2;
            ++ *mults;
        }

        p_mag >>= 1;
        r2 *= r2;
        ++ *mults;
    }

    if ( p_sign == -1 )
    {
        value = 1.00 / value;
        ++ *mults;
    }

	result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *    _r8_pythag ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_PYTHAG computes sqrt ( A^2 + B^2 ), avoiding overflow and underflow.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp A, B, the values for which sqrt ( A^2 + B^2 ) is desired.
    Output, ityp r8_PYTHAG, the value of sqrt ( A^2 + B^2 ).
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp b = a_data[1];
	
    ityp a_abs;
    ityp b_abs;
    ityp value;

    a_abs = abs ( a );
    b_abs = abs ( b );

    if ( b_abs < a_abs )
        value = a_abs * sqrt ( 1.00 + pow ( b_abs / a_abs, 2 ) );
    else if ( b_abs == 0.00 )
        value = 0.00;
    else if ( a_abs <= b_abs )
        value = b_abs * sqrt ( 1.00 + pow ( a_abs / b_abs, 2 ) );

    result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_reverse_bytes ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_REVERSE_BYTES reverses the bytes in an r8.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp X, a value whose bytes are to be reversed.
    Output, r8_REVERSE_BYTES, a value with bytes in reverse order;
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    char c;
    union
    {
        ityp yfloat;
        char ychar[4];
    } y;

    y.yfloat = x;

    c = y.ychar[0];
    y.ychar[0] = y.ychar[3];
    y.ychar[3] = c;

    c = y.ychar[1];
    y.ychar[1] = y.ychar[2];
    y.ychar[2] = c;

	result = ( y.yfloat ); 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_round_i4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_ROUND_I4 rounds an r8, returning an I4.
  Example:
        X         Value
      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 April 2013
  Author:
    John Burkardt
  Parameters:
    Input, ityp X, the value.
    Output, int r8_ROUND_I4, the rounded value.
*/
{
	static int result = INT_MAX;
	
	const register ityp x = *(ityp *) data;
	
	result = floor (   x + 0.5 ) *(1-((x<0.00)<<1));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_round2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_ROUND2 rounds an r8 in base 2.
  Discussion:
    Assume that the input quantity X has the form
      X = S * J * 2^L
    where S is plus or minus 1, L is an integer, and J is a binary
    mantissa which is either exactly zero, or greater than or equal
    to 0.5 and less than 1.0.
    Then on return, XROUND = r8_ROUND2 ( NPLACE, X ) will satisfy
      XROUND = S * K * 2^L
    where S and L are unchanged, and K is a binary mantissa which
    agrees with J in the first NPLACE binary digits and is zero
    thereafter.
    If NPLACE is 0, XROUND will always be zero.
    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.
    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
    or 0.75.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int NPLACE, the number of binary digits to
    preserve.  NPLACE should be 0 or positive.
    Input, ityp X, the real number to be decomposed.
    Output, ityp r8_ROUND2, the rounded value of X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ nplace = s_data->a0;
	const register ityp x = s_data->a1;
	
    dim_typ iplace;
    int l, s;
    const ityp two = 2.00;
    ityp xmant;
    ityp xtemp;
    ityp value;

    /*
    1: Handle the special case of 0.
    */

    if ( x == 0.00 || nplace == 0 )
    {
    	result = 0.00;
        return &result;
    }

    value = 0.00;
    /*
    2: Determine the sign S.
    */
    if ( 0.00 < x )
    {
        s = 1;
        xtemp = x;
    }
    else
    {
        s = -1;
        xtemp = -x;
    }
    /*
    3: Force XTEMP to lie between 1 and 2, and compute the logarithm L.
    */
    l = 0;

    while ( 2.00 <= xtemp )
    {
        xtemp /= 2.00;
        ++ l;
    }

    while ( xtemp < 1.00 )
    {
        xtemp *= 2.00;
        -- l;
    }
    /*
    4: Strip out the digits of the mantissa as XMANT, and decrease L.
    */
    xmant = 0.0;
    iplace = 0;

    for ( ; ; )
    {
        xmant *= 2.00;

        if ( 1.00 <= xtemp )
        {
            ++ xmant;
            -- xtemp;
        }

        ++ iplace;

        if ( xtemp == 0.00 || nplace <= iplace )
        {
            value = s * xmant * pow ( two, l );
            break;
        }

        -- l;
        xtemp *= 2.00;
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_roundb ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_ROUNDB rounds an r8 in a given base.
  Discussion:
    The code does not seem to do a good job of rounding when
    the base is negative
    Assume that the input quantity X has the form
      X = S * J * BASE^L
    where S is plus or minus 1, L is an integer, and J is a
    mantissa base BASE which is either exactly zero, or greater
    than or equal to (1/BASE) and less than 1.0.
    Then on return, XROUND will satisfy
      XROUND = S * K * BASE^L
    where S and L are unchanged, and K is a mantissa base BASE
    which agrees with J in the first NPLACE digits and is zero
    thereafter.
    Note that because of rounding, for most bases, most numbers
    with a fractional quantities cannot be stored exactly in the
    computer, and hence will have trailing "bogus" digits.
    If NPLACE is 0, XROUND will always be zero.
    If NPLACE is 1, the mantissa of XROUND will be 0,
    1/BASE, 2/BASE, ..., (BASE-1)/BASE.
    If NPLACE is 2, the mantissa of XROUND will be 0,
    BASE/BASE^2, (BASE+1)/BASE^2, ...,
    BASE^2-2/BASE^2, BASE^2-1/BASE^2.
  Licensing:
    This code is distributed under the GNU LGPL license.  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int BASE, the base of the arithmetic.
    BASE must not be zero.  Theoretically, BASE may be negative.
    Input, int NPLACE, the number of digits base BASE to
    preserve.  NPLACE should be 0 or positive.
    Input, ityp X, the number to be decomposed.
    Output, ityp r8_ROUNDB, the rounded value of X.
*/
{
	static ityp result = MAX_VAL;
	
	const idtit * const s_data = data;
	int base = s_data->a0;
	const register dim_typ nplace = s_data->a1;
	const register ityp x = s_data->a2;
	
    int iplace;
    int is;
    int js;
    int l;
    ityp r8_base;
    ityp value;
    ityp xmant;
    ityp xtemp;

    r8_base = ( ityp ) base;
    /*
    0: Error checks.
    */
    if ( base == 0 )
    {
    	result = MAX_VAL;
        return &result;
    }
    /*
    1: Handle the special case of 0.
    */

    if ( x == 0.00 || nplace == 0 )
    {
    	result = 0.00;
        return &result;
    }
    value = 0.00;
    /*
    2: Determine the sign IS.
    */
    if ( 0.00 < x )
    {
        is = 1;
        xtemp = x;
    }
    else
    {
        is = -1;
        xtemp = -x;
    }
    /*
    3: Force XTEMP to lie between 1 and ABS(BASE), and compute the
    logarithm L.
    */
    l = 0;

    while ( abs ( r8_base ) <= abs ( xtemp ) )
    {
        xtemp /= r8_base;

        if ( xtemp < 0.00 )
        {
            is *= -1;
            xtemp *= -1;
        }
        ++ l;
    }

    while ( abs ( xtemp ) < 1.00 )
    {
        xtemp *= r8_base;

        if ( xtemp < 0.00 )
        {
            is *= -1;
            xtemp *= -1;
        }

        -- l;
    }
    /*
    4: Now strip out the digits of the mantissa as XMANT, and
    decrease L.
    */
    xmant = 0.00;
    iplace = 0;
    js = is;

    for ( ; ; )
    {
        xmant = r8_base * xmant;

        if ( xmant < 0.00 )
        {
            js *= -1;
            xmant *= -1;
        }

        if ( 1.00 <= xtemp )
        {
            xmant += ( int ) ( xtemp );
            xtemp -= ( int ) ( xtemp );
        }

        ++ iplace;

        if ( xtemp == 0.00 || nplace <= iplace )
        {
            value = ( ityp ) js * xmant * pow ( r8_base, l );
            break;
        }

        -- l;
        xtemp *= r8_base;

        if ( xtemp < 0.00 )
        {
            is *= -1;
            xtemp *= -1;
        }
    }
    
    result = value;
  	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_roundx ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_ROUNDX rounds an r8 in base 10.
  Discussion:
    Assume that the input quantity X has the form
      X = S * J * 10^L
    where S is plus or minus 1, L is an integer, and J is a decimal
    mantissa which is either exactly zero, or greater than or equal
    to 0.1 and less than 1.0.
    Then on return, XROUND will satisfy
      XROUND = S * K * 10^L
    where S and L are unchanged, and K is a decimal mantissa which
    agrees with J in the first NPLACE decimal digits and is zero
    thereafter.
    Note that because of rounding, most decimal fraction quantities
    cannot be stored exactly in the computer, and hence will have
    trailing "bogus" digits.
    If NPLACE is 0, XROUND will always be zero.
    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
    0.2, ..., or 0.9.
    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
    0.03, ..., 0.98, 0.99.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int NPLACE, the number of decimal digits to
    preserve.  NPLACE should be 0 or positive.
    Input, ityp X, the number to be decomposed.
    Output, ityp r8_ROUNDX, the rounded value of X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ nplace = s_data->a0;
	const register ityp x = s_data->a1;
	
    int iplace;
    int is;
    int l;
    const ityp ten = 10.0;
    ityp xmant;
    ityp xround;
    ityp xtemp;

    /*
    1: Handle the special case of 0.
    */

    if ( nplace <= 0 || nplace == 0.00)
    {
    	result = 0.00;
        return &result;
    }
    xround = 0.00;
    /*
    2: Determine the sign IS.
    */
    if ( 0.00 < x )
    {
        is = 1;
        xtemp = x;
    }
    else
    {
        is = -1;
        xtemp = -x;
    }
    /*
    3: Force XTEMP to lie between 1 and 10, and compute the logarithm L.
    */
    l = 0;

    while ( 10.00 <= x )
    {
        xtemp /= 10.00;
        ++ l;
    }

    while ( xtemp < 1.00 )
    {
        xtemp *= 10.00;
        -- l;
    }
    /*
    4: Now strip out the digits of the mantissa as XMANT, and
    decrease L.
    */
    xmant = 0.00;
    iplace = 0;

    for ( ; ; )
    {
        xmant *= 10.00;

        if ( 1.00 <= xtemp )
        {
            xmant += ( int ) xtemp;
            xtemp -= ( int ) xtemp;
        }

        iplace = iplace + 1;

        if ( xtemp == 0.00 || nplace <= iplace )
        {
            xround = is * xmant * pow ( ten, l );
            break;
        }

        -- l;
        xtemp *= 10.00;
    }

	result = xround;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_sign_opposite_strict ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_SIGN_OPPOSITE_STRICT is TRUE if two r8's are strictly of opposite sign.
  Discussion:
    This test could be coded numerically as
      if ( r1 * r2 < 0.0 ) ...
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp R1, R2, the values to check.
    Output, int r8_SIGN_OPPOSITE_STRICT, is TRUE if ( R1 < 0 and 0 < R2 )
    or ( R2 < 0 and 0 < R1 ).
*/
{
	static bool result = 2;
	
	ityp * const a_data = data;
	const register ityp r1 = a_data[0];
	const register ityp r2 = a_data[1];
	
	result = ( r1 < 0.00 && 0.00 < r2 ) || ( r2 < 0.00 && 0.00 < r1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_tiny ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_TINY returns a "tiny" r8.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Output, ityp r8_TINY, a "tiny" r8 value.
*/
{
	static ityp result = 0.1175494350822E-37;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_to_dhms ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_TO_DHMS converts an r8 day value into days, hours, minutes, seconds.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp R, a real number representing a time period measured in days.
    Output, int D, H, M, S, the equivalent number of days, hours,
    minutes and seconds.
*/
{
	const it4pi * const s_data = data;
	ityp r = s_data->a0;
	int * d = s_data->a1;
	int * h = s_data->a2;
	int * m = s_data->a3;
	int * s = s_data->a4;
	
    int sign;

    if ( 0.00 <= r )
        sign = 1;
    else if ( r < 0.00 )
    {
        sign = -1;
        r *= -1;
    }

    *d = ( int ) r;

    r -= ( ityp ) *d;
    r = 24.00 * r;
    *h = ( int ) r;

    r -= ( ityp ) *h;
    r = 60.00 * r;
    *m = ( int ) r;

    r -= ( ityp ) *m;
    r = 60.00 * r;
    *s = ( int ) r;

    if ( sign == -1 )
    {
        *d = -(*d);
        *h = -(*h);
        *m = -(*m);
        *s = -(*s);
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_to_i4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
  Discussion:
    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
    IX := MIN ( IX, MAX ( IXMIN, IXMAX ) )
    IX := MAX ( IX, MIN ( IXMIN, IXMAX ) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp X, the real number to be converted.
    Input, ityp XMIN, XMAX, the real range.  XMAX and XMIN must not be
    equal.  It is not necessary that XMIN be less than XMAX.
    Input, int IXMIN, IXMAX, the allowed range of the output
    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
    It is not necessary that IXMIN be less than IXMAX.
    Output, int r8_TO_I4, the value in the range [IXMIN,IXMAX] that
    corresponds to X.
*/
{
	static int result = INT_MAX;
	
	const _3it2i * const s_data = data;
	ityp x = s_data->a0;
	ityp xmin = s_data->a1;
	ityp xmax = s_data->a2;
	int ixmin = s_data->a3;
	int ixmax = s_data->a4;
	
    ityp temp;

    if ( xmax == xmin )
    {
    	result = INT_MAX;
        return &result;
    }

    temp =( ( xmax - x        ) * ( ityp ) ixmin+ (        x - xmin ) * ( ityp ) ixmax )/ ( xmax     - xmin );
    temp += 0.50 * (1 - ((0.00>temp)<<1));

    result = ( int ) temp;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_to_r8_discrete ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_TO_r8_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
  Discussion:
    if ( R < RMIN ) then
      RD = RMIN
    else if ( RMAX < R ) then
      RD = RMAX
    else
      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )
    In the special case where NR = 1, when
      XD = 0.5 * ( RMAX + RMIN )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp R, the number to be converted.
    Input, ityp RMAX, RMIN, the maximum and minimum
    values for RD.
    Input, int NR, the number of allowed values for XD.
    NR should be at least 1.
    Output, ityp RD, the corresponding discrete value.
*/
{
	static ityp result = MAX_VAL;
	
	const _3iti * const s_data = data;
	ityp r = s_data->a0;
	ityp rmin = s_data->a1;
	ityp rmax = s_data->a2;
	int nr = s_data->a3;
	
    int f;
    ityp rd;
    /*
    Check for errors.
    */
    if ( nr < 1 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( nr == 1 )
    {
    	result = 0.50 * ( rmin + rmax );
        return &result;
    }

    if ( rmax == rmin )
    {
    	result = rmax;
        return &result;
    }

    f = r8_nint ( ( ityp ) ( nr ) * ( rmax - r ) / ( rmax - rmin ) );
    f = MAX ( f, 0 );
    f = MIN ( f, nr );

	result = ( ( ityp ) (      f ) * rmin+ ( ityp ) ( nr - f ) * rmax )/ ( ityp ) ( nr     );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_uniform_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_UNIFORM_AB returns a scaled pseudorandom r8.
  Discussion:
    The pseudorandom number should be uniformly distributed
    between A and B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 April 2011
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
    Input, ityp A, B, the limits of the interval.
    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.
    Output, ityp r8_UNIFORM_AB, a number strictly between A and B.
*/
{
	static ityp result = MAX_VAL;
	
	const _2itpi * const s_data = data;
	const register ityp a = s_data->a0;
	const register ityp b = s_data->a1;
	int * seed = s_data->a2;
	
    int k;

    if ( *seed == 0 )
    {
    	result = MAX_VAL;
        return &result;
    }

    k = *seed / 127773;
    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
        *seed += i4_huge;

	result = a + ( b - a ) * ( ityp ) ( *seed ) * 4.656612875E-10;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _r8_unswap3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_UNSWAP3 unswaps three real items.
  Example:
    Input:
      X = 2, Y = 3, Z = 1
    Output:
      X = 1, Y = 2, Z = 3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input/output, ityp *X, *Y, *Z, three values to be swapped.
*/
{
	ityp ** const a_data = data;
	ityp * x = a_data[0];
	ityp * y = a_data[1];
	ityp * z = a_data[2];
	
    ityp w;

    w = *z;
    *z = *y;
    *y = *x;
    *x =  w;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_walsh_1d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_WALSH_1D evaluates the Walsh function of a real scalar argument.
  Discussion:
    Consider the binary representation of X, and number the digits
    in descending order, from leading to lowest, with the units digit
    being numbered 0.
    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp X, the argument of the Walsh function.
    Input, int DIGIT, the index of the Walsh function.
    Output, ityp r8_WALSH_1D, the value of the Walsh function.
*/
{
	static ityp result = MAX_VAL;
	
	const iti * const s_data = data;
	ityp x = s_data->a0;
	int digit = s_data->a1;
	
    int n;
    const ityp two = 2.00;
    /*
    Hide the effect of the sign of X.
    */
    x = fabs ( x );
    /*
    If DIGIT is positive, divide by 2 DIGIT times.
    If DIGIT is negative, multiply by 2 (-DIGIT) times.
    */
    x /= pow ( two, digit );
    /*
    Make it an integer.
    Because it's positive, and we're using INT, we don't change the
    units digit.
    */
    n = ( int ) x;
    /*
    Is the units digit odd or even?
    */

	result = ( n % 2 ) != 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r82_cheby ( void * data)
/******************************************************************************/
/*
  Purpose:
    r82_CHEBY sets up the Chebyshev abscissas in an r8 interval.
  Discussion:
    The routine sets up a vector of X values spaced between the values
    XLO and XHI in a similar way to the spacing of the Chebyshev
    points of the same order in the interval [-1,1].
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points to compute.
    Input, ityp ALO, AHI, the range.
    Output, ityp r82_CHEBY[N], the computed X values.
*/
{
	const i2it * const s_data = data;
	int n = s_data->a0;
	ityp alo = s_data->a1;
	ityp ahi = s_data->a2;
	
    ityp *a;
    ityp arg;

    a = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
        a[0] = 0.50 * ( alo + ahi );
    else if ( 1 < n )
        for (dim_typ i = 0; i < n; ++i )
        {
            arg = ( ityp ) ( (i<<1) + 1 ) * M_PI / ( ityp ) ( n<<1 );
            a[i] = 0.50 * ( ( 1.00 + cos ( arg ) ) * alo+ ( 1.00 - cos ( arg ) ) * ahi );

        }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r82vec_order_type ( void * data)
/******************************************************************************/
/*
  Purpose:
    r82VEC_ORDER_TYPE finds if an r82VEC is (non)strictly ascending/descending.
  Discussion:
    An r82VEC is a vector whose entries are r82's.
    An r82 is a vector of type ityp precision with two entries.
    An r82VEC may be stored as a 2 by N array.
    The dictionary or lexicographic ordering is used.
 (X1,Y1) < (X2,Y2)  <=>  X1 < X2 or ( X1 = X2 and Y1 < Y2).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries of the array.
    Input, ityp A[2*N], the array to be checked.
    Output, int r82VEC_ORDER_TYPE, order indicator:
    -1, no discernable order;
    0, all entries are equal;
    1, ascending order;
    2, strictly ascending order;
    3, descending order;
    4, strictly descending order.
*/
{
	static short result = SHRT_MAX;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a =  s_data->a1;
	
    int i;
    short order;
    /*
    Search for the first value not equal to A(1,1).
    */
    i = 0;

    for ( ; ; )
    {
        ++ i;

        if ( n <= i )
        {
        	result = 0;
            return &result;
        }

        if ( a[0] < a[0+(i<<1)] || ( a[0] == a[0+(i<<1)] && a[1] < a[1+(i<<1)] ) )
        {
            if ( i == 2 )
                order = 2;
            else
                order = 1;
                break;
        }
        else if ( a[0+(i<<1)] < a[0] || ( a[0+(i<<1)] == a[0] && a[1+(i<<1)] < a[1] ) )
        {
            if ( i == 2 )
                order = 4;
            else
                order = 3;
            break;
        }
    }
    /*
    Now we have a "direction".  Examine subsequent entries.
    */
    for ( ; ; )
    {
        ++ i;
        if ( n <= i )
            break;

        if ( order == 1 )
        {
            if ( a[0+(i<<1)] < a[0+((i-1)<<1)] ||( a[0+(i<<1)] == a[0+((i-1)<<1)] && a[1+(i<<1)] < a[1+((i-1)<<1)] ) )
            {
                order = -1;
                break;
            }
        }
        else if ( order == 2 )
        {
            if ( a[0+(i<<1)] < a[0+((i-1)<<1)] ||( a[0+(i<<1)] == a[0+((i-1)<<1)] && a[1+(i<<1)] < a[1+((i-1)<<1)] ) )
            {
                order = -1;
                break;
            }
            else if ( a[0+(i<<1)] == a[0+((i-1)<<1)] && a[1+(i<<1)] == a[1+((i-1)<<1)] )
                order = 1;
        }
        else if ( order == 3 )
        {
            if ( a[0+((i-1)<<1)] < a[0+(i<<1)] ||( a[0+((i-1)<<1)] == a[0+(i<<1)] && a[1+((i-1)<<1)] < a[1+(i<<1)] ) )
            {
                order = -1;
                break;
            }
        }
        else if ( order == 4 )
            if ( a[0+((i-1)<<1)] < a[0+(i<<1)] ||( a[0+((i-1)<<1)] == a[0+(i<<1)] && a[1+((i-1)<<1)] < a[1+(i<<1)] ) )
            {
                order = -1;
                break;
            }
            else if ( a[0+(i<<1)] == a[0+((i-1)<<1)] && a[1+(i<<1)] == a[1+((i-1)<<1)] )
                order = 3;
    }
    
    result = order;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r82vec_permute ( void * data)
/******************************************************************************/
/*
  Purpose:
    r82VEC_PERMUTE permutes an r82VEC in place.
  Discussion:
    An r82VEC is a vector whose entries are r82's.
    An r82 is a vector of r8's with two entries.
    An r82VEC may be stored as a 2 by N array.
    This routine permutes an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.
  Example:
    Input:
      N = 5
      P = (   2,    4,    5,    1,    3 )
      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
       (11.0, 22.0, 33.0, 44.0, 55.0 )
    Output:
      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
       ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input, int P[N], the permutation.  P(I) = J means
    that the I-th element of the output array should be the J-th
    element of the input array.
    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.
    Input/output, ityp A[2*N], the array to be permuted.
*/
{
	const dtpiipit * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	int base = s_data->a2;
	ityp * a = s_data->a3;
	
    ityp a_temp[2];
    dim_typ i;
    int iget;
    int iput;
    int istart;

    if ( !perm_check ( n, p ) )
        return NULL;
    /*
    In order for the sign negation trick to work, we need to assume that the
    entries of P are strictly positive.  Presumably, the lowest number is BASE.
    So temporarily add 1-BASE to each entry to force positivity.
    */
    for ( i = 0; i < n; ++i )
        p[i] += 1 - base;

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
            a_temp[0] = a[0+((istart-1)<<1)];
            a_temp[1] = a[1+((istart-1)<<1)];
            iget = istart;
            /*
            Copy the new value into the vacated entry.
            */
            for ( ; ; )
            {
                iput = iget;
                iget = p[iget-1];

                p[iput-1] *= -1;

                if ( iget < 1 || n < iget )
                    return NULL;

                if ( iget == istart )
                {
                    a[0+((iput-1)<<1)] = a_temp[0];
                    a[1+((iput-1)<<1)] = a_temp[1];
                    break;
                }
                a[0+((iput-1)<<1)] = a[0+((iget-1)<<1)];
                a[1+((iput-1)<<1)] = a[1+((iget-1)<<1)];
            }
        }
    }
    /*
    Restore the signs of the entries.
    */
    /*
    Restore the base of the entries.
    */
    for ( i = 0; i < n; ++i )
    {
        p[i] *= -1;
        p[i] += - 1 + base;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r82vec_sort_heap_index_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an r82VEC.
  Discussion:
    An r82VEC is a vector whose entries are r82's.
    An r82 is a vector of r8's with two entries.
    An r82VEC may be stored as a 2 by N array.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:
      a(*,indx(*))
    or explicitly, by the call
      r82vec_permute ( n, indx, base, a )
    after which a(*,*) is sorted.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, int BASE, the desired indexing for the sort index:
    0 for 0-based indexing,
    1 for 1-based indexing.
    Input, ityp A[2*N], an array to be index-sorted.
    Output, int r82VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
    I-th element of the sorted array is A(0:1,r82VEC_SORT_HEAP_INDEX_A(I)).
*/
{
	const dtpiti * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	int base = s_data->a2;
	
    ityp aval[2];
    dim_typ i;
    int *indx;
    int indxt;
    int ir;
    dim_typ j;
    int l;

    if ( n < 1 )
        return NULL;

    indx = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
        indx[i] = i;

    if ( n == 1 )
    {
        indx[0] += base;
        return indx;
    }

    l = n / 2 + 1;
    ir = n;

    for ( ; ; )
    {
        if ( 1 < l )
        {
            -- l;
            indxt = indx[l-1];
            aval[0] = a[0+(indxt<<1)];
            aval[1] = a[1+(indxt<<1)];
        }
        else
        {
            indxt = indx[ir-1];
            aval[0] = a[0+(indxt<<1)];
            aval[1] = a[1+(indxt<<1)];
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
            if ( j<ir &&( a[0+(indx[j-1]<<1)] <  a[0+(indx[j]<<1)] ||( a[0+(indx[j-1]<<1)] == a[0+(indx[j]<<1)] &&a[1+(indx[j-1]<<1)] <  a[1+(indx[j]<<1)] ) ))
                ++ j;

            if (   aval[0] <  a[0+(indx[j-1]<<1)] ||( aval[0] == a[0+(indx[j-1]<<1)] &&aval[1] <  a[1+(indx[j-1]<<1)] ) )
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
    /*
    Take care of the base.
    */
    for ( i = 0; i < n; ++i )
        indx[i] += base;

    return indx;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8block_zero_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8BLOCK_ZERO_NEW returns a new zeroed r8BLOCK.
  Discussion:
    An r8BLOCK is a triple dimensioned array of r8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 April 2013
  Author:
    John Burkardt
  Parameters:
    Input, int L, M, N, the number of rows, columns, and levels.
    Output, ityp r8BLOCK_ZERO_NEW[L*M*N], the new zeroed matrix.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ l = a_data[0];
	const register dim_typ m = a_data[1];
	const register dim_typ n = a_data[2];
	
    ityp *a;
    dim_typ i, j, k;
    a = ( ityp * ) malloc ( l * m * n * sizeof ( ityp ) );

    for ( k = 0; k < n; ++k )
        for ( j = 0; j < m; ++j )
            for ( i = 0; i < l; ++i )
                a[i+j*l+k*l*m] = 0.00;
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_compare ( void * data)
/******************************************************************************/
/*
    Purpose:
    r8COL_COMPARE compares two columns in an r8COL.
    Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    Example:
    Input:
    M = 3, N = 4, I = 2, J = 4
    A = (
    1.  2.  3.  4.
    5.  6.  7.  8.
    9. 10. 11. 12. )
    Output:
    r8COL_COMPARE = -1
    Licensing:
    This code is distributed under the GNU LGPL license.
    Modified:
    16 November 2009
    Author:
    John Burkardt
    Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the M by N array.
    Input, int I, J, the columns to be compared.
    I and J must be between 1 and N.
    Output, int r8COL_COMPARE, the results of the comparison:
    -1, column I < column J,
    0, column I = column J,
    +1, column J < column I.
    */
{
	static short result = SHRT_MAX;
	
	const _2dtpit2dt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	const register dim_typ i = s_data->a3;
	const register dim_typ j = s_data->a4;
	
    int k;
    /*
    Check.
    */

    if ( j*i == 0 || n < j || n < i)
    {
    	result = SHRT_MAX;
        return &result;
    }

    if ( i == j )
    {
    	result = 0;
        return &result;
    }

    k = 0;

    while ( k < m )
    {
        if ( a[k+(i-1)*m] < a[k+(j-1)*m] )
        {
        	result = -1;
            return &result;
        }
        else if ( a[k+(j-1)*m] < a[k+(i-1)*m] )
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
__MATHSUITE __JBURKARDT  void   * _r8col_duplicates ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_DUPLICATES generates an r8COL with some duplicate columns.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    This routine generates a random r8COL with a specified number of
    duplicate columns.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in each column of A.
    Input, int N, the number of columns in A.
    Input, int N_UNIQUE, the number of unique columns in A.
    1 <= N_UNIQUE <= N.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, ityp r8COL_DUPLICATES[M*N], the array.
*/
{
	const _3dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	dim_typ n_unique = s_data->a2;
	int * seed = s_data->a3;
	
    ityp *a;
    dim_typ i;
    dim_typ j1, j2;
    ityp temp;

    if ( n_unique < 1 || n < n_unique )
        return NULL;

    a = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    r8mat_uniform_01 ( m, n_unique, seed, a );
    /*
    Randomly copy unique columns.
    */
    for ( j1 = n_unique; j1 < n; j1++ )
    {
        j2 = i4_uniform ( 0, n_unique - 1, seed );
        for ( i = 0; i < m; ++i)
            a[i+j1*m] = a[i+j2*m];
    }
    /*
    Permute the columns.
    */
    for ( j1 = 0; j1 < n; ++ j1)
    {
        j2 = i4_uniform ( j1, n - 1, seed );
        for ( i = 0; i < m; ++i )
        {
            temp      = a[i+j1*m];
            a[i+j1*m] = a[i+j2*m];
            a[i+j2*m] = temp;
        }
    }
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_find ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_FIND seeks a column value in an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Example:
    Input:
      M = 3,
      N = 4,
      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )
      x = ( 3.,
            7.,
           11. )
    Output:
      r8COL_FIND = 3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input int M, N, the number of rows and columns.
    Input, ityp A[M*N], a table of numbers, regarded as
    N columns of vectors of length M.
    Input, ityp X[M], a vector to be matched with a column of A.
    Output, int r8COL_FIND, the (one-based) index of the first column of A
    which exactly matches every entry of X, or -1 if no match
    could be found.
*/
{
	static short result = SHRT_MAX;
	
	const _2dt2pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * x = s_data->a3;
	
    short col = -1;
    dim_typ i, j;

    for ( j = 1; j <= n; j++ )
    {
        col = j;

        for ( i = 1; i <= m; ++i )
        {
            if ( x[i-1] != a[i-1+(j-1)*m] )
            {
                col = -1;
                break;
            }
        }
        if ( col != -1 )
        {
        	result = col;
            return &result;
        }
    }
    
    result = col;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_first_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_FIRST_INDEX indexes the first occurrence of values in an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    For element A(1:M,J) of the matrix, FIRST_INDEX(J) is the index in A of
    the first column whose entries are equal to A(1:M,J).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of A.
    The length of an "element" of A, and the number of "elements".
    Input, ityp A[M*N], the array.
    Input, ityp TOL, a tolerance for equality.
    Output, int r8COL_FIRST_INDEX[N], the first occurrence index.
*/
{
	const _2dtpitit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp tol = s_data->a3;
	
    ityp diff;
    int *first_index;
    dim_typ i, j1, j2;

    first_index = ( int * ) malloc ( n * sizeof ( int ) );

    for ( j1 = 0; j1 < n; ++j1 )
        first_index[j1] = -1;
    for ( j1 = 0; j1 < n; ++j1 )
    {
        if ( first_index[j1] == -1 )
        {
            first_index[j1] = j1;

            for ( j2 = j1 + 1; j2 < n; ++j2 )
            {
                diff = 0.00;
                for ( i = 0; i < m; ++i )
                    diff = MAX ( diff, abs ( a[i+j1*m] - a[i+j2*m] ) );
                if ( diff <= tol )
                    first_index[j2] = j1;
            }
        }
    }
    return first_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_INSERT inserts a column into an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Example:
    Input:
      N_MAX = 10,
      M = 3,
      N = 4,
      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )
      X = ( 3., 4., 18. )
    Output:
      N = 5,
      A = (
        1.  2.  3.  3.  4.
        5.  6.  4.  7.  8.
        9. 10. 18. 11. 12. )
      r8COL_INSERT = 3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N_MAX, the maximum number of columns in A.
    Input, int M, the number of rows.
    Input/output, int N, the number of columns.
    If the new column is inserted into the table, then the output
    value of N will be increased by 1.
    Input/output, ityp A[M*N_MAX], a table of numbers, regarded
    as an array of columns.  The columns must have been sorted
    lexicographically.
    Input, ityp X[M], a vector of data which will be inserted
    into the table if it does not already occur.
    Output, int r8COL_INSERT.
    I, X was inserted into column I.
    -I, column I was already equal to X.
    0, N = N_MAX.
*/
{
	static short result = SHRT_MAX;
	
	const _3dt2pit * const s_data = data;
	const register dim_typ n_max = s_data->a0;
	const register dim_typ m = s_data->a1;
	register dim_typ n = s_data->a2;
	ityp * a = s_data->a3;
	ityp * x = s_data->a4;
	
	
    short col;
    int high;
    dim_typ i;
    int isgn;
    dim_typ j;
    int low;
    int mid;
    /*
    Refuse to work if N_MAX <= N.
    */
    if ( n_max <= n )
    {
    	result = 0;
        return &result;
    }
    /*
    Stick X temporarily in column N+1, just so it's easy to use r8COL_COMPARE.
    */
    for ( i = 0; i < m; ++i )
        a[i+n*m] = x[i];
    /*
    Do a binary search.
    */
    low = 1;
    high = n;

    for ( ; ; )
    {
        if ( high < low )
        {
            col = low;
            break;
        }

        mid = ( low + high ) / 2;

        isgn = r8col_compare ( m, n+1, a, mid, n+1 );

        if ( isgn == 0 )
        {
        	result = -mid;
            return &result;
        }
        else if ( isgn == -1 )
            low = mid + 1;
        else if ( isgn == +1 )
            high = mid - 1;
    }
    /*
    Shift part of the table up to make room.
    */
    for ( j = n-1; col-1 <= j; --j )
        for ( i = 0; i < m; ++i )
            a[i+(j+1)*m] = a[i+j*m];
    /*
    Insert the new column.
    */
    for ( i = 0; i < m; ++i )
        a[i+(col-1)*m] = x[i];

    ++ n;

	result = col;
    return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_MAX returns the column maximums of an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the array to be examined.
    Output, ityp r8COL_MAX[N], the maximums of the columns.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    ityp *amax;
    dim_typ i, j;

    amax = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
    {
        amax[j] = a[0+j*m];
        for ( i = 0; i < m; ++i )
            amax[j] = MAX ( amax[j], a[i+j*m] );
    }

    return amax;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_max_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_MAX_INDEX returns the indices of column maximums in an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the array to be examined.
    Output, int r8COL_MAX_INDEX[N]; entry I is the row of A in which
    the maximum for column I occurs.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    ityp amax;
    dim_typ i;
    dim_typ *imax;
    dim_typ j;

    imax = ( dim_typ * ) malloc ( n * sizeof ( dim_typ ) );

    for ( j = 0; j < n; ++j )
    {
        imax[j] = 1;
        amax = a[0+j*m];

        for ( i = 1; i < m; ++i )
            if ( amax < a[i+j*m] )
            {
                imax[j] = i+1;
                amax = a[i+j*m];
            }
    }

    return imax;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8col_max_one ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_MAX_ONE rescales an r8COL so each column maximum is 1.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A[M*N], the array to be rescaled.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i;
    dim_typ i_big;
    dim_typ j;
    ityp temp;

    for ( j = 0; j < n; ++j )
    {
        i_big = 0;
        for ( i = 1; i < m; ++i )
            if ( abs ( a[i_big+j*m] ) < abs ( a[i+j*m] ) )
                i_big = i;
        temp = a[i_big+j*m];

        if ( temp )
            for ( i = 0; i < m; ++i )
                a[i+j*m] /= temp;
    }
    
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_min ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_MIN returns the column minimums of an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the array to be examined.
    Output, ityp r8COL_MIN[N], the minimums of the columns.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    ityp *amin;
    dim_typ i, j;

    amin = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
    {
        amin[j] = a[0+j*m];
        for ( i = 0; i < m; ++i )
            amin[j] = MIN ( amin[j], a[i+j*m] );
    }

    return amin;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_min_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_MIN_INDEX returns the indices of column minimums in an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the array to be examined.
    Output, int r8COL_MIN_INDEX[N]; entry I is the row of A in which
    the minimum for column I occurs.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    ityp amin;
    dim_typ i;
    dim_typ*imin;
    dim_typ j;

    imin = ( dim_typ * ) malloc ( n * sizeof ( dim_typ ) );

    for ( j = 0; j < n; ++j )
    {
        imin[j] = 1;
        amin = a[0+j*m];

        for ( i = 1; i < m; ++i )
            if ( a[i+j*m] < amin )
            {
                imin[j] = i+1;
                amin = a[i+j*m];
            }
    }

    return imin;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_part_quick_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_PART_QUICK_A reorders the columns of an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    The routine reorders the columns of A.  Using A(1:M,1) as a
    key, all entries of A that are less than or equal to the key will
    precede the key, which precedes all entries that are greater than the key.
  Example:
    Input:
      M = 2, N = 8
      A = ( 2  8  6  0 10 10  0  5
            4  8  2  2  6  0  6  8 )
    Output:
      L = 2, R = 4
      A = (  0  0  2  8  6 10 10  4
             2  6  4  8  2  6  0  8 )
             ----     -------------
             LEFT KEY     RIGHT
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, the row dimension of A, and the length of a column.
    Input, int N, the column dimension of A.
    Input/output, ityp A[M*N].  On input, the array to be checked.
    On output, A has been reordered as described above.
    Output, int *L, *R, the indices of A that define the three segments.
    Let KEY = the input value of A(1:M,1).  Then
    I <= L                 A(1:M,I) < KEY;
         L < I < R         A(1:M,I) = KEY;
                 R <= I    KEY < A(1:M,I).
*/
{
	const _2dtpit2pdt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	dim_typ * l = s_data->a3;
	dim_typ * r = s_data->a4;
	
    dim_typ i, j, k;
    ityp *key;

    if ( n == 0)
        return NULL;

    if ( n == 1 )
    {
        *l = 0;
        *r = 2;
        return NULL;
    }

    key = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
        key[i] = a[i];

    k = 1;
    /*
    The elements of unknown size have indices between L+1 and R-1.
    */
    *l = 1;
    *r = n + 1;

    for ( j = 1; j < n; ++j )
    {
        if ( r8vec_gt ( m, a+(*l)*m, key ) )
        {
            -- *r;
            r8vec_swap ( m, a+(*r-1)*m, a+(*l)*m );
        }
        else if ( r8vec_eq ( m, a+(*l)*m, key ) )
        {
            ++ k;
            r8vec_swap ( m, a+(k-1)*m, a+(*l)*m );
            ++ *l;
        }
        else if ( r8vec_lt ( m, a+(*l)*m, key ) )
            ++ *l;
    }
    /*
    Shift small elements to the left.
    */
    for ( j = 0; j < *l - k; ++j)
        for ( i = 0; i < m; ++i )
            a[i+j*m] = a[i+(j+k)*m];
    /*
    Shift KEY elements to center.
    */
    for ( j = *l-k; j < *l; ++j )
        for ( i = 0; i < m; ++i)
                a[i+j*m] = key[i];
    /*
    Update L.
    */
    *l -= k;

    free ( key );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_permute ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_PERMUTE permutes an r8COL in place.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    This routine permutes an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.
  Example:
    Input:
      M = 2
      N = 5
      P = (   2,    4,    5,    1,    3 )
      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
       (11.0, 22.0, 33.0, 44.0, 55.0 )
      BASE = 1
    Output:
      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
       ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, the length of objects.
    Input, int N, the number of objects.
    Input, int P[N], the permutation.  P(I) = J means
    that the I-th element of the output array should be the J-th
    element of the input array.
    Input, int BASE, is 0 for a 0-based permutation and 1 for a
    1-based permutation.
    Input/output, ityp A[M*N], the array to be permuted.
*/
{
	const _2dtpiipit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * p = s_data->a2;
	int base = s_data->a3;
	ityp * a = s_data->a4;
	
    ityp *a_temp;
    dim_typ i;
    int iget;
    int iput;
    dim_typ istart;
    dim_typ j;

    if ( !perm_check ( n, p ) )
        return NULL;
    /*
    In order for the sign negation trick to work, we need to assume that the
    entries of P are strictly positive.  Presumably, the lowest number is BASE.
    So temporarily add 1-BASE to each entry to force positivity.
    */
    for ( i = 0; i < n; ++i )
        p[i] += 1 - base;

    a_temp = ( ityp * ) malloc ( m * sizeof ( ityp ) );
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
            for ( i = 0; i < m; ++i )
                a_temp[i] = a[i+(istart-1)*m];
            iget = istart;
            /*
            Copy the new value into the vacated entry.
            */
            for ( ; ; )
            {
                iput = iget;
                iget = p[iget-1];

                p[iput-1] *= -1;

                if ( iget < 1 || n < iget )
                    return NULL;

                if ( iget == istart )
                {
                    for ( i = 0; i < m; ++i)
                        a[i+(iput-1)*m] = a_temp[i];
                    break;
                }
                for ( i = 0; i < m; ++i)
                    a[i+(iput-1)*m] = a[i+(iget-1)*m];
            }
        }
    }
    /*
    Restore the signs of the entries.
    */
    /*
    Restore the base of the entries.
    */
    for ( i = 0; i < n; ++i )
    {
        p[i] *= -1;
        p[i] += -1 +  base;
    }

    free ( a_temp );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sort_heap_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORT_HEAP_A ascending heapsorts an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
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
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A[M*N].
    On input, the array of N columns of M-vectors.
    On output, the columns of A have been sorted in lexicographic order.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;

    if ( n == 0 || n == 1 || m == 0 || m == 1)
        return NULL;
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
            r8col_swap ( m, n, a, i, j );
        /*
        Compare the I and J objects.
        */
        else if ( indx < 0 )
            isgn = r8col_compare ( m, n, a, i, j );
        else if ( indx == 0 )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_sort_heap_index_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.
    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2) is negative.
    Once the index array is computed, the sorting can be carried out
    "implicitly:
      A(*,INDX(*)) is sorted,
    Note that the index vector is 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in each column of A.
    Input, int N, the number of columns in A.
    Input, int BASE, the desired indexing for the sort index:
    0 for 0-based indexing,
    1 for 1-based indexing.
    Input, ityp A[M*N], the array.
    Output, int r8COL_SORT_HEAP_INDEX_A[N], contains the sort index.  The
    I-th column of the sorted array is A(*,INDX(I)).
*/
{
	const _2dtipit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int base = s_data->a2;
	ityp * a = s_data->a3;
	
    ityp *column;
    dim_typ i;
    int *indx;
    int indxt;
    int ir;
    short isgn;
    dim_typ j;
    dim_typ k;
    dim_typ l;

    if ( n == 0)
        return NULL;

    indx = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
        indx[i] = i;


    if ( n == 1 )
    {
        indx[0] = indx[0] + base;
        return indx;
    }

    column = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    l = n / 2 + 1;
    ir = n;

    for ( ; ; )
    {
        if ( 1 < l )
        {
            -- l;
            indxt = indx[l-1];
            for ( k = 0; k < m; ++k )
                column[k] = a[k+indxt*m];
        }
        else
        {
            indxt = indx[ir-1];
            for ( k = 0; k < m; ++k )
                column[k] = a[k+indxt*m];
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
            if ( j < ir )
            {
                isgn = r8vec_compare ( m, a+indx[j-1]*m, a+indx[j]*m );

                if ( isgn < 0 )
                    ++ j;
            }

            isgn = r8vec_compare ( m, column, a+indx[j-1]*m );

            if ( isgn < 0 )
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
    free ( column );
    /*
    Take care of the base.
    */
    for ( i = 0; i < n; ++i )
        indx[i] += base;

    return indx;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sort_quick_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORT_QUICK_A ascending quick sorts an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, the row order of A, and the length of a column.
    Input, int N, the number of columns of A.
    Input/output, ityp A[M*N].
    On input, the array to be sorted.
    On output, the array has been sorted.
*/
{
	const _2dtipit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int base = s_data->a2;
	ityp * a = s_data->a3;
	
    # define LEVEL_MAX 30

    dim_typ l_segment;
    dim_typ level;
    int n_segment;
    int rsave[LEVEL_MAX];
    dim_typ r_segment;

    if ( n == 1 || n == 0 || m == 0)
        return NULL;

    level = base = 1;
    rsave[level-1] = n + 1;
    n_segment = n;

    for ( ; ; )
    {
        /*
        Partition the segment.
        */
        r8col_part_quick_a ( m, n_segment, a+(base-1)*m, &l_segment, &r_segment );
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
                if ( level <= 1 )
                    return NULL;

                base = rsave[level-1];
                n_segment = rsave[level-2] - rsave[level-1];
                -- level;

                if ( 0 < n_segment )
                    break;
            }
        }
    }
    return NULL;
    # undef LEVEL_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sorted_tol_undex ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORTED_TOL_UNDEX returns tolerably unique indexes for a sorted r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    The goal of this routine is to determine a vector UNDX,
    which points, to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.
    This is all done with index vectors, so that the elements of
    A are never moved.
    Assuming A is already sorted, we examine the entries of A in order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.
    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula
      XU(*) = A(UNDX(*)).
    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:
      A(I) = XU(XDNU(I)).
    We could then replace A by the combination of XU and XDNU.
    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.
    Here is an example of a vector A, the unique sort and inverse unique
    sort vectors and the compressed unique sorted vector.
      I      A      XU  Undx  Xdnu
    ----+------+------+-----+-----+
      0 | 11.0 |  11.0    0     0
      1 | 11.0 |  22.0    4     0
      2 | 11.0 |  33.0    7     0
      3 | 11.0 |  55.0    8     0
      4 | 22.0 |                1
      5 | 22.0 |                1
      6 | 22.0 |                1
      7 | 33.0 |                2
      8 | 55.0 |                3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of the data values.
    Input, int N, the number of data values,
    Input, ityp A[M*N], the data values.
    Input, int UNIQUE_NUM, the number of unique values in A.
    This value is only required for languages in which the size of
    UNDX must be known in advance.
    Input, ityp TOL, a tolerance for equality.
    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
    Output, int XDNU[N], the XDNU vector.
*/
{
	const _2dtpitdtit2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	dim_typ unique_num = s_data->a3;
	ityp tol = s_data->a4;
	int * undx = s_data->a5;
	int * xdnu = s_data->a6;
	
    ityp diff;
    dim_typ i;
    dim_typ i2;
    dim_typ i3;
    dim_typ j;
    dim_typ k;
    bool unique;
    /*
    Consider entry I = 0.
    It is unique, so set the number of unique items to K.
    Set the K-th unique item to I.
    Set the representative of item I to the K-th unique item.
    */
    i = k = 0;
    undx[k] = i;
    xdnu[i] = k;
    /*
    Consider entry I.

    If it is unique, increase the unique count K, set the
    K-th unique item to I, and set the representative of I to K.

    If it is not unique, set the representative of item I to a
    previously determined unique item that is close to it.
    */
    for ( i = 1; i < n; ++i )
    {
        unique = 1;

        for ( j = 0; j <= k; ++j )
        {
            i2 = undx[j];
            diff = 0.00;
            for ( i3 = 0; i3 < m; ++i3 )
                diff = MAX ( diff, abs ( a[i3+i*m] - a[i3+i2*m] ) );
            if ( diff <= tol )
            {
                unique = 0;
                xdnu[i] = j;
                break;
            }
        }
        if ( unique )
        {
            ++ k;
            undx[k] = i;
            xdnu[i] = k;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sorted_tol_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORTED_TOL_UNIQUE keeps tolerably unique elements in a sorted r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    The columns of the array can be ascending or descending sorted.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A(M,N).
    On input, the sorted array of N columns of M-vectors.
    On output, a sorted array of columns of M-vectors.
    Input, ityp TOL, a tolerance for equality.
    Output, int r8COL_SORTED_TOL_UNIQUE, the number of unique columns.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpitit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp tol = s_data->a3;
	
    ityp diff;
    dim_typ i, j, k;
    bool unique;
    dim_typ unique_num;

    if ( n <= 0 )
    {
    	result = 0;
        return &result;
    }

    unique_num = 1;

    for ( i = 1; i < n; ++i )
    {
        unique = true;
        for ( j = 0; j < unique_num; ++j )
        {
            diff = 0.00;
            for ( k = 0; k < m; ++k )
                diff = MAX ( diff, abs ( a[k+i*m] - a[k+j*m] ) );
            if ( diff < tol )
            {
                unique = false;
                break;
            }
        }
        if ( unique )
        {
            for ( k = 0; k < m; ++k )
                a[k+unique_num*m] = a[k+i*m];
            ++ unique_num;
        }
    }
    
    result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sorted_tol_unique_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORTED_TOL_UNIQUE_COUNT counts tolerably unique elements in a sorted r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    The columns of the array may be ascending or descending sorted.
    If the tolerance is large enough, then the concept of uniqueness
    can become ambiguous.  If we have a tolerance of 1.5, then in the
    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
    one unique entry?  That would be because 1 may be regarded as unique,
    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
    be unique and so on.
    This seems wrongheaded.  So I prefer the idea that an item is not
    unique under a tolerance only if it is close to something that IS unique.
    Thus, the unique items are guaranteed to cover the space if we include
    a disk of radius TOL around each one.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], a sorted array, containing
    N columns of data.
    Input, ityp TOL, a tolerance for equality.
    Output, int r8COL_SORTED_TOL_UNIQUE_COUNT, the number of unique columns.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpitit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp tol = s_data->a3;
	
    ityp diff;
    dim_typ i;
    dim_typ i2;
    dim_typ i3;
    dim_typ j;
    dim_typ k;
    int *undx;
    bool unique;

    undx = ( int * ) malloc ( n * sizeof ( int ) );
    /*
    Consider entry I = 0.
    It is unique, so set the number of unique items to K.
    Set the K-th unique item to I.
    Set the representative of item I to the K-th unique item.
    */
    i = k = 0;
    undx[k] = i;
    /*
    Consider entry I.

    If it is unique, increase the unique count K, set the
    K-th unique item to I, and set the representative of I to K.

    If it is not unique, set the representative of item I to a
    previously determined unique item that is close to it.
    */
    for ( i = 1; i < n; ++i )
    {
        unique = true;

        for ( j = 0; j <= k; ++j )
        {
            i2 = undx[j];
            diff = 0.00;
            for ( i3 = 0; i3 < m; ++i3 )
                diff = MAX ( diff, abs ( a[i3+i*m] - a[i3+i2*m] ) );
            if ( diff <= tol )
            {
                unique = false;
                break;
            }
        }
        if ( unique )
        {
            ++ k;
            undx[k] = i;
        }
    }
    free ( undx );

	result = k + 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sorted_undex ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORTED_UNDEX returns unique indexes for a sorted r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    The goal of this routine is to determine a vector UNDX,
    which points, to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.
    This is all done with index vectors, so that the elements of
    A are never moved.
    Assuming A is already sorted, we examine the entries of A in order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.
    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula
      XU(*) = A(UNDX(*)).
    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:
      A(I) = XU(XDNU(I)).
    We could then replace A by the combination of XU and XDNU.
    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.
    Here is an example of a vector A, the unique sort and inverse unique
    sort vectors and the compressed unique sorted vector.
      I      A      XU  Undx  Xdnu
    ----+------+------+-----+-----+
      0 | 11.0 |  11.0    0     0
      1 | 11.0 |  22.0    4     0
      2 | 11.0 |  33.0    7     0
      3 | 11.0 |  55.0    8     0
      4 | 22.0 |                1
      5 | 22.0 |                1
      6 | 22.0 |                1
      7 | 33.0 |                2
      8 | 55.0 |                3
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the dimension of the data values.
    Input, int N, the number of data values,
    Input, ityp A[M*N], the data values.
    Input, int UNIQUE_NUM, the number of unique values in A.
    This value is only required for languages in which the size of
    UNDX must be known in advance.
    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
    Output, int XDNU[N], the XDNU vector.
*/
{
	const dtipdtdt2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	int n = s_data->a1;
	dim_typ * a = s_data->a2; 
	const register dim_typ unique_num = s_data->a3;
	int * undx = s_data->a4;
	int * xdnu = s_data->a5;
	
    ityp diff;
    dim_typ i, j, k;
    /*
    Walk through the sorted array.
    */
    i = j = 0;
    undx[j] = i;
    xdnu[i] = j;

    for ( i = 1; i < n; ++i )
    {
        diff = 0.00;
        for ( k = 0; k < m; ++k )
            diff = MAX ( diff, abs ( a[k+i*m] - a[k+undx[j]*m] ) );
        if ( 0.00 < diff )
        {
            ++ j;
            undx[j] = i;
        }
        xdnu[i] = j;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sorted_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORTED_UNIQUE keeps unique elements in a sorted r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    The columns of the array can be ascending or descending sorted.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A(M,N).
    On input, the sorted array of N columns of M-vectors.
    On output, a sorted array of columns of M-vectors.
    Output, int UNIQUE_NUM, the number of unique columns.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    bool equal;
    dim_typ i;
    dim_typ j1;
    dim_typ j2;

    if ( n <= 0 )
    {
    	result = 0;
        return &result;
    }

    j1 = 0;

    for ( j2 = 1; j2 < n; ++j2 )
    {
        equal = true;
        for ( i = 0; i < m; ++i)
            if ( a[i+j1*m] != a[i+j2*m] )
            {
                equal = false;
                break;
            }
        if ( !equal )
        {
            ++ j1;
            for ( i = 0; i < m; ++i )
                a[i+j1*m] = a[i+j2*m];
        }
    }

	result = j1 + 1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sorted_unique_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORTED_UNIQUE_COUNT counts unique elements in a sorted r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    The columns of the array may be ascending or descending sorted.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], a sorted array, containing
    N columns of data.
    Output, int r8COL_SORTED_UNIQUE_COUNT, the number of unique columns.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    bool equal;
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
    {
        equal = true;
        for ( i = 0; i < m; ++i )
            if ( a[i+j1*m] != a[i+j2*m] )
            {
                equal = false;
                break;
            }
        if ( !equal )
        {
            ++ unique_num;
            j1 = j2;
        }
    }

	result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8col_sortr_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SORTR_A ascending sorts one column of an r8COL, adjusting all entries.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A[M*N].
    On input, an unsorted M by N array.
    On output, rows of the array have been shifted in such
    a way that column KEY of the array is in nondecreasing order.
    Input, int KEY, the column in which the "key" value
    is stored.  On output, column KEY of the array will be
    in nondecreasing order.
*/
{
	const _2dtipit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int key = s_data->a2;
	ityp * a = s_data->a3;
	
    int i, j;
    int indx;
    short isgn;

    if ( m == 0 || key < 1 || n < key )
        return NULL;
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
        r8col_swap ( m, n, a, i, j );
    /*
    Compare the I and J objects.
    */
    else if ( indx < 0 )
        isgn = 1 - ((a[i-1+(key-1)*m] < a[j-1+(key-1)*m])<<1);
    else if ( indx == 0 )
        break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_sum ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SUM sums the columns of an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the array to be examined.
    Output, ityp r8COL_SUM[N], the sums of the columns.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    ityp *colsum;
    dim_typ i, j;

    colsum = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
    {
        colsum[j] = 0.00;
        for ( i = 0; i < m; ++i )
            colsum[j] += a[i+j*m];
    }
    return colsum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8col_swap ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_SWAP swaps columns J1 and J2 of an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
  Example:
    Input:
      M = 3, N = 4, J1 = 2, J2 = 4
      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )
    Output:
      A = (
        1.  4.  3.  2.
        5.  8.  7.  6.
        9. 12. 11. 10. )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A[M*N], the M by N array.
    Input, int J1, J2, the columns to be swapped.
    These columns are 1-based.
*/
{
	const _2dtpit2dt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	const register dim_typ j1 = s_data->a3;
	const register dim_typ j2 = s_data->a4;
	
    ityp temp;

    if ( j1 == 0 || n < j1 || j2 < 1 || n < j2 || j1 == j2 )
        return NULL;

    for (dim_typ i = 0; i < m; ++i)
    {
        temp          = a[i+(j1-1)*m];
        a[i+(j1-1)*m] = a[i+(j2-1)*m];
        a[i+(j2-1)*m] = temp;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_to_r8vec ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_TO_r8VEC converts an r8COL to an r8VEC.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
    each of length M.
    This routine is not really useful in our C++ implementation, since
    we actually store an M by N matrix exactly as a vector already.
  Example:
    M = 3, N = 4
    A =
      11 12 13 14
      21 22 23 24
      31 32 33 34
    r8COL_TO_r8VEC = ( 11, 21, 31, 12, 22, 32, 13, 23, 33, 14, 24, 34 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the M by N array.
    Output, ityp X[M*N], a vector containing the N columns of A.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j, k;
    ityp *x;

    x = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    k = 0;
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
        {
            x[k] = a[i+j*m];
            ++ k;
        }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8col_tol_undex ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_TOL_UNDEX indexes tolerably unique entries in an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded as an array of N columns,
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
      A(I) = XU(XDNU(I)).
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
    Input, ityp A[M*N], the data values.
    Input, int UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.
    Input, ityp TOL, a tolerance for equality.
    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
    Output, int XDNU[N], the XDNU vector.
*/
{
	const _2dtpitdtit2pi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	dim_typ unique_num = s_data->a3;
	ityp tol = s_data->a4;
	int * undx = s_data->a5;
	int * xdnu = s_data->a6;
	
    int base = 0;
    ityp diff;
    dim_typ i;
    dim_typ i2;
    int *indx;
    dim_typ j;
    dim_typ k;
    bool unique;
    /*
    Implicitly sort the array.
    */
    indx = r8col_sort_heap_index_a ( m, n, base, a );
    /*
    Consider entry I = 0.
    It is unique, so set the number of unique items to K.
    Set the K-th unique item to I.
    Set the representative of item I to the K-th unique item.
    */
    i = k = 0;
    undx[k] = indx[i];
    xdnu[indx[i]] = k;
    /*
    Consider entry I.

    If it is unique, increase the unique count K, set the
    K-th unique item to I, and set the representative of I to K.

    If it is not unique, set the representative of item I to a
    previously determined unique item that is close to it.
    */
    for ( i = 1; i < n; ++i )
    {
        unique = true;
        for ( j = 0; j <= k; ++j )
        {
            diff = 0.00;
            for ( i2 = 0; i2 < m; ++i2 )
                diff = MAX ( diff, abs ( a[i2+indx[i]*m] - a[i2+undx[j]*m] ) );
            if ( diff <= tol )
            {
                unique = false;
                xdnu[indx[i]] = j;
                break;
            }
        }
        if ( unique )
        {
            ++ k;
            undx[k] = indx[i];
            xdnu[indx[i]] = k;
        }
    }
    free ( indx );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_copy ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_COPY copies one r8MAT to another.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8's, which
    may be stored as a vector in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input ityp A1[M*N], the matrix to be copied.
    Output, ityp A2[M*N], the copy of A1.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a1 = s_data->a2;
	ityp * a2 = s_data->a3;
	
    dim_typ i, j;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            a2[i+j*m] = a1[i+j*m];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_delete ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_DELETE frees the memory set aside by r8MAT_NEW.
  Discussion:
    This function releases the memory associated with an array that was
    created by a command like
      ityp **a;
      a = r8mat_new ( m, n );
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 September 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp **A, the array.
    Input, int M, N, the number of rows and columns.
*/
{
	const ppit2dt * const s_data = data;
	ityp ** a = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    for (dim_typ i = 0; i < m; ++i )
        free ( a[i] );

    free ( a );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_det ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_DET computes the determinant of an r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    Original FORTRAN77 version by Helmut Spaeth.
    C version by John Burkardt.
  Reference:
    Helmut Spaeth,
    Cluster Analysis Algorithms
    for Data Reduction and Classification of Objects,
    Ellis Horwood, 1980, page 125-127.
  Parameters:
    Input, int N, the order of the matrix.
    Input, ityp A[N*N], the matrix whose determinant is desired.
    Output, ityp r8MAT_DET, the determinant of the matrix.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp *b;
    ityp det;
    int i, j, k;
    dim_typ kk;
    dim_typ m;
    ityp temp;

    b = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i)
            b[i+j*n] = a[i+j*n];

    det = 1.00;

    for ( k = 1; k <= n; ++k )
    {
        m = k;
        for ( kk = k+1; kk <= n; kk++ )
            if ( abs ( b[m-1+(k-1)*n] ) < abs ( b[kk-1+(k-1)*n] ) )
                m = kk;

        if ( m != k )
        {
            det *= -1;

            temp = b[m-1+(k-1)*n];
            b[m-1+(k-1)*n] = b[k-1+(k-1)*n];
            b[k-1+(k-1)*n] = temp;
        }

        det = det * b[k-1+(k-1)*n];

        if ( b[k-1+(k-1)*n] != 0.00 )
        {
            for ( i = k+1; i <= n; ++i )
                b[i-1+(k-1)*n] *= -1/b[k-1+(k-1)*n];

            for ( j = k+1; j <= n; ++j )
            {
                if ( m != k )
                {
                    temp = b[m-1+(j-1)*n];
                    b[m-1+(j-1)*n] = b[k-1+(j-1)*n];
                    b[k-1+(j-1)*n] = temp;
                }
                for ( i = k+1; i <= n; ++i )
                    b[i-1+(j-1)*n] += b[i-1+(k-1)*n] * b[k-1+(j-1)*n];
            }
        }
    }

    free ( b );

	result = det;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_det_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_DET_2D computes the determinant of a 2 by 2 r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
  Discussion:
    The determinant of a 2 by 2 matrix is
      a11 * a22 - a12 * a21.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp A[2*2], the matrix whose determinant is desired.
    Output, ityp r8MAT_DET_2D, the determinant of the matrix.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * a = data;

	result = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_det_3d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_DET_3D computes the determinant of a 3 by 3 r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
    The determinant of a 3 by 3 matrix is
        a11 * a22 * a33 - a11 * a23 * a32
      + a12 * a23 * a31 - a12 * a21 * a33
      + a13 * a21 * a32 - a13 * a22 * a31
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp A[3*3], the matrix whose determinant is desired.
    Output, ityp r8MAT_DET_3D, the determinant of the matrix.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * a = data;
	
	result = a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
    + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
    + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_det_4d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_DET_4D computes the determinant of a 4 by 4 r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp A[4*4], the matrix whose determinant is desired.
    Output, ityp r8MAT_DET_4D, the determinant of the matrix.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * a = data;
	
	result = a[0+0*4] * (
        a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
        - a[0+1*4] * (
        a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
        + a[0+2*4] * (
        a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
        - a[0+3*4] * (
        a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_det_5d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_DET_5D computes the determinant of a 5 by 5 r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp A[5*5], the matrix whose determinant is desired.
    Output, ityp r8MAT_DET_5D, the determinant of the matrix.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * a = data;
	
    ityp b[4*4];
    ityp det = 0.00;
    dim_typ i, j, k;
    char sign = 1; 
    
    /*
    Expand the determinant into the sum of the determinants of the
    five 4 by 4 matrices created by dropping row 1, and column k.
    */

    for ( k = 0; k < 5; ++k)
    {
        for ( i = 0; i < 4; ++i )
            for ( j = 0; j < 4; ++j )
                b[i+(j<<2)] = a[i+1+(j+(j>=k))*5];

        det += sign * a[0+k*5] * r8mat_det_4d ( b );
        sign *= -1;
    }

	result = det;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_flip_cols ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_FLIP_COLS swaps the columns of an r8MAT.
  Discussion:
    An r8MAT is an MxN array of r8's, stored by (I,J) -> [I+J*M].
    To "flip" the columns of an r8MAT is to start with something like
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
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A[M*N], the matrix whose columns are to be flipped.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp t;

    for ( i = 0; i < m; ++i )
        for ( j = 0; j < ( n / 2 ); ++j)
        {
            t              = a[i+     j *m];
            a[i+     j *m] = a[i+(n-1-j)*m];
            a[i+(n-1-j)*m] = t;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_flip_rows ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_FLIP_ROWS swaps the rows of an r8MAT.
  Discussion:
    An r8MAT is an MxN array of r8's, stored by (I,J) -> [I+J*M].
    To "flip" the rows of an r8MAT is to start with something like
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
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A[M*N], the matrix whose rows are to be flipped.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp t;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < ( m / 2 ); ++i)
        {
            t            = a[    i+j*m];
            a[    i+j*m] = a[m-1-i+j*m];
            a[m-1-i+j*m] = t;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_identity_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_IDENTITY sets the square matrix A to the identity.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8's, which
    may be stored as a vector in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of A.
    Output, ityp A[N*N], the N by N identity matrix.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a;
    dim_typ i, j, k;

    a = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = k = 0; j < n; ++j)
        for ( i = 0; i < n; ++i )
        {
            a[k] = i == j;
            ++ k;
        }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_in_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_IN_01 is TRUE if the entries of an r8MAT are in the range [0,1].
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, ityp A[M*N], the matrix.
    Output, int r8MAT_IN_01, is TRUE if every entry of A is
    between 0 and 1.
*/
{
	static bool result = 2;
	
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            if ( a[i+j*m] < 0.00 || 1.00 < a[i+j*m] )
            {
            	result = false;
                return &result;
            }

	result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_indicator_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_INDICATOR_NEW sets up an "indicator" r8MAT.
  Discussion:
    An R8MAT is an array of R8's.
    The value of each entry suggests its location, as in:
      11  12  13  14
      21  22  23  24
      31  32  33  34
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 April 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows of the matrix.
    M must be positive.
    Input, int N, the number of columns of the matrix.
    N must be positive.
    Output, ityp r8MAT_INDICATOR_NEW[M*N], the table.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
    ityp *a;
    int fac;
    dim_typ i, j;

    a = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    fac = powi ( 10, i4_log_10 ( n ) + 1 );

    for ( i = 1; i <= m; ++i )
        for ( j = 1; j <= n; ++j )
            a[i-1+(j-1)*m] = ( ityp ) ( fac * i + j );
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_inverse_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_INVERSE_2D inverts a 2 by 2 r8MAT using Cramer's rule.
  Discussion:
    The two dimensional array is stored as a one dimensional vector,
    by COLUMNS.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp A[2*2], the matrix to be inverted.
    Output, ityp r8MAT_INVERSE_2D[2*2], the inverse of the matrix A.
*/
{
	ityp * a = data;
	
    ityp *b;
    const register ityp det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];
    /*
    Compute the determinant of A.
    */
    /*
    If the determinant is zero, bail out.
    */
    if ( det == 0.00 )
        return NULL;
    /*
    Compute the entries of the inverse matrix using an explicit formula.
    */
    b = ( ityp * ) malloc ( sizeof ( ityp ) << 2 );

    b[0+0*2] = + a[1+1*2] / det;
    b[0+1*2] = - a[0+1*2] / det;
    b[1+0*2] = - a[1+0*2] / det;
    b[1+1*2] = + a[0+0*2] / det;

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_inverse_3d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_INVERSE_3D inverts a 3 by 3 r8MAT using Cramer's rule.
  Discussion:
    The two dimensional array is stored as a one dimensional vector,
    by COLUMNS.
    If the determinant is zero, A is singular, and does not have an
    inverse.  In that case, the output is set to NULL.
    If the determinant is nonzero, its value is an estimate
    of how nonsingular the matrix A is.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp A[3*3], the matrix to be inverted.
    Output, ityp r8MAT3_INVERSE[3*3], the inverse of the matrix A.
*/
{
	ityp * a = data;
	
    ityp *b;
    /*
    Compute the determinant of A.
    */
    const register ityp det =
    a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
    + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
    + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );

    if ( det == 0.00 )
        return NULL;

    b = ( ityp * ) malloc ( 9 * sizeof ( ityp ) );

    b[0+0*3] = ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) / det;
    b[0+1*3] = - ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] ) / det;
    b[0+2*3] = ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] ) / det;

    b[1+0*3] = - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] ) / det;
    b[1+1*3] = ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] ) / det;
    b[1+2*3] = - ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] ) / det;

    b[2+0*3] = ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] ) / det;
    b[2+1*3] = - ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] ) / det;
    b[2+2*3] = ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] ) / det;

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_inverse_4d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_INVERSE_4D inverts a 4 by 4 matrix using Cramer's rule.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 September 2003
  Author:
    John Burkardt
  Parameters:
    Input, ityp A[4][4], the matrix to be inverted.
    Output, ityp r8MAT_INVERSE_4D[4][4], the inverse of the matrix A.
*/
{
	ityp * a = data;
	
    ityp *b;
    ityp det;
    /*
    Compute the determinant of A.
    */
    det = r8mat_det_4d ( a );
    /*
    If the determinant is zero, bail out.
    */
    if ( det == 0.00 )
        return NULL;
    /*
    Compute the entries of the inverse matrix using an explicit formula.
    */
    b = ( ityp * ) malloc ( sizeof ( ityp ) << 2 );

    b[0+0*4] =
        +(
        + a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        + a[1+2*4] * ( a[2+3*4] * a[3+1*4] - a[2+1*4] * a[3+3*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        ) / det;

    b[1+0*4] =
        -(
        + a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        + a[1+2*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        ) / det;

    b[2+0*4] =
        +(
        + a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+1*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
        ) / det;

    b[3+0*4] =
        -(
        + a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        + a[1+1*4] * ( a[2+2*4] * a[3+0*4] - a[2+0*4] * a[3+2*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
        ) / det;

    b[0+1*4] =
        -(
        + a[0+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        + a[0+2*4] * ( a[2+3*4] * a[3+1*4] - a[2+1*4] * a[3+3*4] )
        + a[0+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        ) / det;

    b[1+1*4] =
        +(
        + a[0+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        + a[0+2*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
        + a[0+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        ) / det;

    b[2+1*4] =
        -(
        + a[0+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[0+1*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
        + a[0+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
        ) / det;

    b[3+1*4] =
        +(
        + a[0+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        + a[0+1*4] * ( a[2+2*4] * a[3+0*4] - a[2+0*4] * a[3+2*4] )
        + a[0+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
        ) / det;

    b[0+2*4] =
        +(
        + a[0+1*4] * ( a[1+2*4] * a[3+3*4] - a[1+3*4] * a[3+2*4] )
        + a[0+2*4] * ( a[1+3*4] * a[3+1*4] - a[1+1*4] * a[3+3*4] )
        + a[0+3*4] * ( a[1+1*4] * a[3+2*4] - a[1+2*4] * a[3+1*4] )
        ) / det;

    b[1+2*4] =
        -(
        + a[0+0*4] * ( a[1+2*4] * a[3+3*4] - a[1+3*4] * a[3+2*4] )
        + a[0+2*4] * ( a[1+3*4] * a[3+0*4] - a[1+0*4] * a[3+3*4] )
        + a[0+3*4] * ( a[1+0*4] * a[3+2*4] - a[1+2*4] * a[3+0*4] )
        ) / det;

    b[2+2*4] =
        +(
        + a[0+0*4] * ( a[1+1*4] * a[3+3*4] - a[1+3*4] * a[3+1*4] )
        + a[0+1*4] * ( a[1+3*4] * a[3+0*4] - a[1+0*4] * a[3+3*4] )
        + a[0+3*4] * ( a[1+0*4] * a[3+1*4] - a[1+1*4] * a[3+0*4] )
        ) / det;

    b[3+2*4] =
        -(
        + a[0+0*4] * ( a[1+1*4] * a[3+2*4] - a[1+2*4] * a[3+1*4] )
        + a[0+1*4] * ( a[1+2*4] * a[3+0*4] - a[1+0*4] * a[3+2*4] )
        + a[0+2*4] * ( a[1+0*4] * a[3+1*4] - a[1+1*4] * a[3+0*4] )
        ) / det;

    b[0+3*4] =
        -(
        + a[0+1*4] * ( a[1+2*4] * a[2+3*4] - a[1+3*4] * a[2+2*4] )
        + a[0+2*4] * ( a[1+3*4] * a[2+1*4] - a[1+1*4] * a[2+3*4] )
        + a[0+3*4] * ( a[1+1*4] * a[2+2*4] - a[1+2*4] * a[2+1*4] )
        ) / det;

    b[1+3*4] =
            +(
            + a[0+0*4] * ( a[1+2*4] * a[2+3*4] - a[1+3*4] * a[2+2*4] )
            + a[0+2*4] * ( a[1+3*4] * a[2+0*4] - a[1+0*4] * a[2+3*4] )
            + a[0+3*4] * ( a[1+0*4] * a[2+2*4] - a[1+2*4] * a[2+0*4] )
            ) / det;

    b[2+3*4] =
        -(
        + a[0+0*4] * ( a[1+1*4] * a[2+3*4] - a[1+3*4] * a[2+1*4] )
        + a[0+1*4] * ( a[1+3*4] * a[2+0*4] - a[1+0*4] * a[2+3*4] )
        + a[0+3*4] * ( a[1+0*4] * a[2+1*4] - a[1+1*4] * a[2+0*4] )
        ) / det;

    b[3+3*4] =
        +(
        + a[0+0*4] * ( a[1+1*4] * a[2+2*4] - a[1+2*4] * a[2+1*4] )
        + a[0+1*4] * ( a[1+2*4] * a[2+0*4] - a[1+0*4] * a[2+2*4] )
        + a[0+2*4] * ( a[1+0*4] * a[2+1*4] - a[1+1*4] * a[2+0*4] )
        ) / det;

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_l_inverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_L_INVERSE inverts a lower triangular r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
    A lower triangular matrix is a matrix whose only nonzero entries
    occur on or below the diagonal.
    The inverse of a lower triangular matrix is a lower triangular matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 May 2010
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, number of rows and columns in the matrix.
    Input, ityp A[N*N], the lower triangular matrix.
    Output, ityp r8MAT_L_INVERSE[N*N], the inverse matrix.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp *b;
    dim_typ i, j, k;
    ityp temp;

    b = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i)
            if ( i < j )
                b[i+j*n] = 0.00;
            else if ( j == i )
                b[i+j*n] = 1.00 / a[i+j*n];
            else
            {
                temp = 0.00;
                for ( k = 0; k < i; ++k )
                    temp += a[i+k*n] * b[k+j*n];
                    b[i+j*n] = -temp / a[i+i*n];
            }

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_l_solve ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_L_SOLVE solves a lower triangular linear system.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 June 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of rows and columns of
    the matrix A.
    Input, ityp A[N*N], the N by N lower triangular matrix.
    Input, ityp B[N], the right hand side of the linear system.
    Output, ityp r8MAT_L_SOLVE[N], the solution of the linear system.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    ityp dot;
    dim_typ i, j;
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    /*
    Solve L * x = b.
    */
    for ( i = 0; i < n; ++i )
    {
        dot = 0.00;
        for ( j = 0; j < i; ++j )
            dot += a[i+j*n] * x[j];
        x[i] = ( b[i] - dot ) / a[i+i*n];
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_lt_solve ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_LT_SOLVE solves a transposed lower triangular linear system.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values,  stored as a vector
    in column-major order.
    Given the lower triangular matrix A, the linear system to be solved is:
      A' * x = b
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 April 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of rows and columns of the matrix A.
    Input, ityp A[N*N], the N by N lower triangular matrix.
    Input, ityp B[N], the right hand side of the linear system.
    Output, ityp r8MAT_LT_SOLVE[N], the solution of the linear system.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    dim_typ i, j;
    ityp *x = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = n-1; 0 <= j; --j )
    {
        x[j] = b[j];
        for ( i = j+1; i < n; ++i )
            x[j] -= x[i] * a[i+j*n];
        x[j] /= a[j+j*n];
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_MAX returns the maximum entry of an r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, ityp A[M*N], the M by N matrix.
    Output, ityp r8MAT_MAX, the maximum entry of A.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp value = a[0+0*m];

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            if ( value < a[i+j*m] )
                value = a[i+j*m];
                
	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_min ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_MIN returns the minimum entry of an r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 May 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, ityp A[M*N], the M by N matrix.
    Output, ityp r8MAT_MIN, the minimum entry of A.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp value = a[0+0*m];

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i)
            if ( a[i+j*m] < value )
                value = a[i+j*m];
                
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_mm ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_MM multiplies two matrices.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values,  stored as a vector
    in column-major order.
    For this routine, the result is returned as the function value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 April 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, N3, the order of the matrices.
    Input, ityp A[N1*N2], ityp B[N2*N3], the matrices to multiply.
    Output, ityp C[N1*N3], the product matrix C = A * B.
*/
{
	const _3dt3pit * const s_data = data;
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	const register dim_typ n3 = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	ityp * c = s_data->a5;
	
    dim_typ i, j, k;

    for ( i = 0; i < n1; ++i)
        for ( j = 0; j < n3; ++j )
        {
            c[i+j*n1] = 0.0;
            for ( k = 0; k < n2; ++k )
                c[i+j*n1] += a[i+k*n1] * b[k+j*n2];
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_mv ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_MV multiplies a matrix times a vector.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8's, which
    may be stored as a vector in column-major order.
    For this routine, the result is returned as the function value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 April 2007
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of the matrix.
    Input, ityp A[M,N], the M by N matrix.
    Input, ityp X[N], the vector to be multiplied by A.
    Output, ityp r8MAT_MV[M], the product A*X.
*/
{
	const _2dt3pit * const s_data = data; 
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	
	dim_typ i, j;
	
    for ( i = 0; i < m; ++i )
    {
        y[i] = 0.00;
        for ( j = 0; j < n; ++j )
            y[i] += a[i+j*m] * x[j];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _r8mat_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_NEW sets up an r8MAT of the desired dimensions.
  Discussion:
    A declaration of the form
      ityp **a;
    is necesary.  Then an assignment of the form:
      a = r8mat_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17.0;
      y = a[1][0];
    and so on.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 September 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Output, ityp **r8MAT_NEW, the array.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
    ityp **a;
    dim_typ i;

    a = ( ityp ** ) malloc ( m * n * sizeof ( ityp * ) );

    for ( i = 0; i < m; ++i)
        a[i] = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_norm_eis ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_NORM_EIS returns the EISPACK norm of an r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
    The EISPACK norm is defined as:
      r8MAT_NORM_EIS =
        sum ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, ityp A[M*N], the matrix whose EISPACK norm is desired.
    Output, ityp r8MAT_NORM_EIS, the EISPACK norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp value = 0.00;
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
            value += abs ( a[i+j*m] );

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_norm_fro ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_NORM_FRO returns the Frobenius norm of an r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8's, which
    may be stored as a vector in column-major order.
    The Frobenius norm is defined as
      r8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
    The matrix Frobenius norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:
      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 July 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, the number of rows in A.
    Input, int N, the number of columns in A.
    Input, ityp A[M*N], the matrix whose Frobenius
    norm is desired.
    Output, ityp r8MAT_NORM_FRO, the Frobenius norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp value = 0.00;
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i)
            value += pow ( a[i+j*m], 2 );

	result = sqrt ( value );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8mat_mxm ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_MXM multiplies two matrices.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values,  stored as a vector
    in column-major order.
    For this routine, the result is returned as an argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 October 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N1, N2, N3, the order of the matrices.
    Input, ityp A[N1*N2], ityp B[N2*N3], the matrices to multiply.
    Output, ityp C[N1*N3], the product matrix C = A * B.
*/
{
	const _3dt3pit * const s_data = data;
	const register dim_typ n1 = s_data->a0;
	const register dim_typ n2 = s_data->a1;
	const register dim_typ n3 = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	ityp * c = s_data->a5;
	
    dim_typ i, j, k;

    for ( i = 0; i < n1; ++i )
        for ( j = 0; j < n3; ++j )
        {
            c[i+j*n1] = 0.0;
            for ( k = 0; k < n2; ++k)
                c[i+j*n1] += a[i+k*n1] * b[k+j*n2];
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_nullspace ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_NULLSPACE computes the nullspace of a matrix.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values,  stored as a vector
    in column-major order.
    Let A be an MxN matrix.
    If X is an N-vector, and A*X = 0, then X is a null vector of A.
    The set of all null vectors of A is called the nullspace of A.
    The 0 vector is always in the null space.
    If the 0 vector is the only vector in the nullspace of A, then A
    is said to have maximum column rank. (Because A*X=0 can be regarded
    as a linear combination of the columns of A).  In particular, if A
    is square, and has maximum column rank, it is nonsingular.
    The dimension of the nullspace is the number of linearly independent
    vectors that span the nullspace.  If A has maximum column rank,
    its nullspace has dimension 0.
    This routine uses the reduced row echelon form of A to determine
    a set of NULLSPACE_SIZE independent null vectors.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 October 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of
    the matrix A.
    Input, ityp A[M*N], the matrix to be analyzed.
    Input, int NULLSPACE_SIZE, the size of the nullspace.
    Output, ityp r8MAT_NULLSPACE[N*NULLSPACE_SIZE], vectors that
    span the nullspace.
*/
{
	const _3dtpit * const s_data = data;
	
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	const register dim_typ nullspace_size = s_data->a2;
	ityp * a = s_data->a3;
	
    int *col;
    dim_typ i;
    dim_typ i2;
    dim_typ j;
    dim_typ j2;
    ityp *nullspace;
    int *row;
    ityp *rref;
    /*
    Make a copy of A.
    */
    rref = r8mat_copy_new ( m, n, a );
    /*
    Get the reduced row echelon form of A.
    */
    r8mat_rref ( m, n, rref );
    /*
    Note in ROW the columns of the leading nonzeros.
    COL(J) = +J if there is a leading 1 in that column, and -J otherwise.
    */
    row = ( int * ) malloc ( m * sizeof ( int ) );

    for ( i = 0; i < m; ++i )
        row[i] = 0;

    col = ( int * ) malloc ( n * sizeof ( int ) );

    for ( j = 0; j < n; ++j )
        col[j] = - ( j + 1 );

    for ( i = 0; i < m; ++i )
        for ( j = 0; j < n; ++j )
            if ( rref[i+j*m] == 1.00 )
            {
                row[i] = ( j + 1 );
                col[j] = ( j + 1 );
                break;
            }

    nullspace = r8mat_zero_new ( n, nullspace_size );
    j2 = 0;
    /*
    If column J does not contain a leading 1, then it contains
    information about a null vector.
    */
    for ( j = 0; j < n; ++j)
    {
        if ( col[j] < 0 )
        {
            for ( i = 0; i < m; ++i )
                if ( rref[i+j*m] != 0.00 )
                {
                    i2 = row[i] - 1;
                    nullspace[i2+j2*n] = - rref[i+j*m];
                }
            nullspace[j+j2*n] = 1.00;
            ++ j2;
        }
    }
    free ( col );
    free ( row );
    free ( rref );

    return nullspace;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_nullspace_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values,  stored as a vector
    in column-major order.
    Let A be an MxN matrix.
    If X is an N-vector, and A*X = 0, then X is a null vector of A.
    The set of all null vectors of A is called the nullspace of A.
    The 0 vector is always in the null space.
    If the 0 vector is the only vector in the nullspace of A, then A
    is said to have maximum column rank. (Because A*X=0 can be regarded
    as a linear combination of the columns of A).  In particular, if A
    is square, and has maximum column rank, it is nonsingular.
    The dimension of the nullspace is the number of linearly independent
    vectors that span the nullspace.  If A has maximum column rank,
    its nullspace has dimension 0.
    This routine ESTIMATES the dimension of the nullspace.  Cases of
    singularity that depend on exact arithmetic will probably be missed.
    The nullspace will be estimated by counting the leading 1's in the
    reduced row echelon form of A, and subtracting this from N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 October 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of
    the matrix A.
    Input, ityp A[M*N], the matrix to be analyzed.
    Output, int r8MAT_NULLSPACE_SIZE, the estimated size
    of the nullspace.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    dim_typ leading;
    ityp *rref;
    /*
    Make a copy of A.
    */
    rref = r8mat_copy_new ( m, n, a );
    /*
    Get the reduced row echelon form of A.
    */
    r8mat_rref ( m, n, rref );
    /*
    Count the leading 1's in A.
    */
    leading = 0;
    for ( i = 0; i < m; ++i )
        for ( j = 0; j < n; ++j )
            if ( rref[i+j*m] == 1.00 )
            {
                ++ leading;
                break;
            }

	result = n - leading;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_ref ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_REF computes the row echelon form of a matrix.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values,  stored as a vector
    in column-major order.
    A matrix is in row echelon form if:
    * The first nonzero entry in each row is 1.
    * The leading 1 in a given row occurs in a column to
      the right of the leading 1 in the previous row.
    * Rows which are entirely zero must occur last.
  Example:
    Input matrix:
     1.0  3.0  0.0  2.0  6.0  3.0  1.0
    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
     3.0  9.0  0.0  0.0  6.0  6.0  2.0
    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
    Output matrix:
     1.0  3.0  0.0  2.0  6.0  3.0  1.0
     0.0  0.0  0.0  1.0  2.0  4.5  1.5
     0.0  0.0  0.0  0.0  0.0  1.0  0.3
     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 October 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of
    the matrix A.
    Input/output, ityp A[M*N].  On input, the matrix to be
    analyzed.  On output, the REF form of the matrix.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    dim_typ lead;
    dim_typ r;
    ityp temp;

    lead = 0;

    for ( r = 0; r < m; ++r)
    {
        if ( n - 1 < lead )
            break;

        i = r;

        while ( a[i+lead*m] == 0.0 )
        {
            ++ i;

            if ( m - 1 < i )
            {
                i = r;
                ++ lead;
                if ( n - 1 < lead )
                {
                    lead = -1;
                    break;
                }
            }
        }

        if ( lead < 0 )
            break;

        for ( j = 0; j < n; ++j )
        {
            temp     = a[i+j*m];
            a[i+j*m] = a[r+j*m];
            a[r+j*m] = temp;
        }

        temp = a[r+lead*m];

        for ( j = 0; j < n; ++j )
            a[r+j*m] /= temp;

        for ( i = r + 1; i < m; ++i)
        {
            temp = a[i+lead*m];
            for ( j = 0; j < n; ++j)
                a[i+j*m] = a[i+j*m] - temp * a[r+j*m];
        }
        ++ lead;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8mat_rref ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_RREF computes the reduced row echelon form of a matrix.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values,  stored as a vector
    in column-major order.
    A matrix is in row echelon form if:
    * The first nonzero entry in each row is 1.
    * The leading 1 in a given row occurs in a column to
      the right of the leading 1 in the previous row.
    * Rows which are entirely zero must occur last.
    The matrix is in reduced row echelon form if, in addition to
    the first three conditions, it also satisfies:
    * Each column containing a leading 1 has no other nonzero entries.
  Example:
    Input matrix:
     1.0  3.0  0.0  2.0  6.0  3.0  1.0
    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
     3.0  9.0  0.0  0.0  6.0  6.0  2.0
    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
    Output matrix:
     1.0  3.0  0.0  0.0  2.0  0.0  0.0
     0.0  0.0  0.0  1.0  2.0  0.0  0.0
     0.0  0.0  0.0  0.0  0.0  1.0  0.3
     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 October 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of
    the matrix A.
    Input/output, ityp A[M*N].  On input, the matrix to be
    analyzed.  On output, the RREF form of the matrix.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i;
    dim_typ j;
    dim_typ lead;
    dim_typ r;
    ityp temp;

    lead = 0;

    for ( r = 0; r < m; ++r )
    {
        if ( n - 1 < lead )
            break;

        i = r;

        while ( a[i+lead*m] == 0.00 )
        {
            ++ i;

            if ( m - 1 < i )
            {
                i = r;
                ++ lead;
                if ( n - 1 < lead )
                {
                    lead = -1;
                    break;
                }
            }
        }

        if ( lead < 0 )
            break;

        for ( j = 0; j < n; ++j )
        {
            temp     = a[i+j*m];
            a[i+j*m] = a[r+j*m];
            a[r+j*m] = temp;
        }

        temp = a[r+lead*m];

        for ( j = 0; j < n; ++j )
            a[r+j*m] /= temp;

        for ( i = 0; i < m; ++i )
            if ( i != r )
	        {
	            temp = a[i+lead*m];
	            for ( j = 0; j < n; ++j )
	                a[i+j*m] -= temp * a[r+j*m];
            }
        ++ lead;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_symm_eigen ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_SYMM_EIGEN returns a symmetric matrix with given eigensystem.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
    The user must supply the desired eigenvalue vector, and the desired
    eigenvector matrix.  The eigenvector matrix must be orthogonal.  A
    suitable random orthogonal matrix can be generated by
    r8MAT_ORTH_UNIFORM_NEW.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of A.
    Input, ityp X[N], the desired eigenvalues for the matrix.
    Input, ityp Q[N*N], the eigenvector matrix of A.
    Output, ityp r8MAT_SYMM_EIGEN[N*N], a symmetric N by N matrix with
    eigenvalues X and eigenvectors the columns of Q.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * q = s_data->a2;
	
    ityp *a;
    dim_typ i, j, k;
    /*
    Set A = Q * Lambda * Q'.
    */
    a = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j )
        {
            a[i+j*n] = 0.00;
            for ( k = 0; k < n; ++k )
                a[i+j*n] += q[i+k*n] * x[k] * q[j+k*n];
        }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_symm_jacobi ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_SYMM_JACOBI applies Jacobi eigenvalue iteration to a symmetric matrix.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
    This code was modified so that it treats as zero the off-diagonal
    elements that are sufficiently close to, but not exactly, zero.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 October 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of A.
    Input/output, ityp A[N*N], a symmetric N by N matrix.
    On output, the matrix has been overwritten by an approximately
    diagonal matrix, with the eigenvalues on the diagonal.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp c;
    const register ityp eps = 0.00001;
    dim_typ i;
    dim_typ it;
    dim_typ j;
    dim_typ k;
    ityp s;
    ityp sum2;
    ityp t;
    ityp t1;
    ityp t2;
    ityp u;

    const register ityp norm_fro = r8mat_norm_fro ( n, n, a );

    it = 0;

    for ( ; ; )
    {
        ++ it;

        for ( i = 0; i < n; ++i )
            for ( j = 0; j < i; ++j )

        if ( eps * norm_fro < abs ( a[i+j*n] ) + abs ( a[j+i*n] ) )
        {
            u = ( a[j+j*n] - a[i+i*n] ) / ( a[i+j*n] + a[j+i*n] );

            t = r8_sign ( u ) / ( abs ( u ) + sqrt ( u * u + 1.00 ) );
            c = 1.00 / sqrt ( t * t + 1.00 );
            s = t * c;
            /*
            A -> A * Q.
            */
            for ( k = 0; k < n; ++k )
            {
                t1 = a[i+k*n];
                t2 = a[j+k*n];
                a[i+k*n] = t1 * c - t2 * s;
                a[j+k*n] = t1 * s + t2 * c;
            }
            /*
            A -> QT * A
            */
            for ( k = 0; k < n; ++k )
            {
                t1 = a[k+i*n];
                t2 = a[k+j*n];
                a[k+i*n] = c * t1 - s * t2;
                a[k+j*n] = s * t1 + c * t2;
            }
        }
        /*
        Test the size of the off-diagonal elements.
        */
        sum2 = 0.00;
        for ( i = 0; i < n; ++i )
            for ( j = 0; j < i; ++j )
                sum2 += abs ( a[i+j*n] );

		if ( sum2 <= eps * ( norm_fro + 1.00 ) || 100 <= it )
		    break;

    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_to_r8plu ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_TO_r8PLU factors a general matrix.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
    This routine is a simplified version of the LINPACK routine DGEFA.
    Fortran conventions are used to index doubly-dimensioned arrays.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2010
  Author:
    John Burkardt
  Reference:
    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979
  Parameters:
    Input, int N, the order of the matrix.
    N must be positive.
    Input, ityp A[N*N], the matrix to be factored.
    Output, int PIVOT[N], a vector of pivot indices.
    Output, ityp LU[N*N], an upper triangular matrix U and the multipliers
    L which were used to obtain it.  The factorization can be written
    A = L * U, where L is a product of permutation and unit lower
    triangular matrices and U is upper triangular.
    Output, int r8MAT_TO_r8PLU, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the r8MAT_TO_r8PLU-th step.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const pitdtpipit * const s_data = data;
	
	ityp * a = s_data->a0; 
	const register dim_typ n = s_data->a1;
	int * pivot = s_data->a2;
	ityp * lu = s_data->a3;
	
    dim_typ i;
    dim_typ info;
    dim_typ j;
    dim_typ k;
    dim_typ l;
    ityp temp;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            lu[i+j*n] = a[i+j*n];
    info = 0;

    for ( k = 1; k <= n-1; ++k )
    {
        /*
        Find L, the index of the pivot row.
        */
        l = k;
        for ( i = k+1; i <= n; ++i )
            if ( abs ( lu[l-1+(k-1)*n] ) < abs ( lu[i-1+(k-1)*n] ) )
                l = i;

        pivot[k-1] = l;
        /*
        If the pivot index is zero, the algorithm has failed.
        */
        if ( lu[l-1+(k-1)*n] == 0.00 )
        {
        	result = k;
            return &result;
        }
        /*
        Interchange rows L and K if necessary.
        */
        if ( l != k )
        {
            temp            = lu[l-1+(k-1)*n];
            lu[l-1+(k-1)*n] = lu[k-1+(k-1)*n];
            lu[k-1+(k-1)*n] = temp;
        }
        /*
        Normalize the values that lie below the pivot entry A(K,K).
        */
        for ( i = k+1; i <= n; ++i )
            lu[i-1+(k-1)*n] /= -lu[k-1+(k-1)*n];
        /*
        Row elimination with column indexing.
        */
        for ( j = k+1; j <= n; ++j )
        {
            if ( l != k )
            {
                temp            = lu[l-1+(j-1)*n];
                lu[l-1+(j-1)*n] = lu[k-1+(j-1)*n];
                lu[k-1+(j-1)*n] = temp;
            }

            for ( i = k+1; i <= n; ++i )
                lu[i-1+(j-1)*n] = lu[i-1+(j-1)*n] + lu[i-1+(k-1)*n] * lu[k-1+(j-1)*n];

        }
    }

    pivot[n-1] = n;

    if ( lu[n-1+(n-1)*n] == 0.00 )
        info = n;

	result = info;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_trace ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_TRACE computes the trace of an r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8 values, stored as a vector
    in column-major order.
    The trace of a square matrix is the sum of the diagonal elements.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix A.
    Input, ityp A[N*N], the matrix whose trace is desired.
    Output, ityp r8MAT_TRACE, the trace of the matrix.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i)
        value += a[i+i*n];
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8mat_transpose ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_TRANSPOSE returns the transpose of a matrix.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8's, which
    may be stored as a vector in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 June 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of the matrix A.
    Input, ityp A[M*N], the matrix whose transpose is desired.
    Output, ityp r8MAT_TRANSPOSE[N*M], the transposed matrix.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    ityp *b;
    dim_typ i, j;

    b = ( ityp * ) malloc ( n * m * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
        for ( i = 0; i < m; ++i )
            b[j+i*n] = a[i+j*m];
    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_transpose_in_place ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_TRANSPOSE_IN_PLACE transposes a square matrix in place.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8's, which
    may be stored as a vector in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 June 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of rows and columns of the matrix A.
    Input/output, ityp A[N*N], the matrix to be transposed.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i, j;
    ityp t;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < j; ++i )
        {
            t        = a[i+j*n];
            a[i+j*n] = a[j+i*n];
            a[j+i*n] = t;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_zero ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8MAT_ZERO zeroes an r8MAT.
  Discussion:
    An r8MAT is a doubly dimensioned array of r8's, which
    may be stored as a vector in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 September 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Output, ityp A[M*N], a matrix of zeroes.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;

    for ( j = 0; j < n; ++j )
        for ( i = 0; i < m; ++i)
            a[i+j*m] = 0.00;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8plu_det ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8PLU_DET computes the determinant of a real PLU matrix.
  Discussion:
    The matrix should have been factored by r8MAT_TO_r8PLU.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2010
  Author:
    John Burkardt
  Reference:
    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979
  Parameters:
    Input, int N, the order of the matrix.
    N must be positive.
    Input, int PIVOT[N], the pivot vector computed by r8MAT_TO_r8PLU.
    Input, ityp LU[N*N], the LU factors computed by r8MAT_TO_r8PLU.
    Output, ityp r8PLU_DET, the determinant of the matrix.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * lu = s_data->a1;
	int * pivot = s_data->a2;
	
    ityp det = 1.00;

    for (dim_typ i = 0; i < n; ++i )
    {
        det *= lu[i+i*n];
        if ( pivot[i] != i+1 )
            det *= -1;
    }

	result = det;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8pp_delete ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8PP_DELETE frees the memory set aside by r8PP_NEW.
  Discussion:
    An r8PP is a pointer to pointers to r8's, and is a sort of
    variably-dimensioned matrix.
    This function releases the memory associated with an r8PP that was
    created by a command like:
      ityp **a;
      a = r8pp_new ( m, n );
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, ityp **A, the pointer to the pointers.
    Input, int M, N, the number of rows and columns.
*/
{
	const ppit2dt * const s_data = data;
	ityp ** a = s_data->a0;
	const register dim_typ m = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    for (dim_typ i = 0; i < m; ++i)
        free ( a[i] );
    free ( a );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _r8pp_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8PP_NEW sets up an r8PP.
  Discussion:
    An r8PP is a pointer to pointers to r8's, and is a sort of
    variably-dimensioned matrix.
    A declaration of the form
      ityp **a;
    is necesary.  Then an assignment of the form:
      a = r8pp_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17;
      y = a[1][0];
    and so on.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    07 November 2011
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Output, ityp **r8PP_NEW, a pointer to the pointers to the M by N array.
*/
{
	dim_typ * const a_data = data; 
	const register dim_typ m = a_data[0];
	const register dim_typ n = a_data[1];
	
    ityp **a = ( ityp ** ) malloc ( m * n * sizeof ( ityp * ) );

    for (dim_typ i = 0; i < m; ++i )
        a[i] = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8r8_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8r8_COMPARE compares two r8r8's.
  Discussion:
    An r8r8 is simply a pair of r8 values, stored separately.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp X1, Y1, the first vector.
    Input, ityp X2, Y2, the second vector.
    Output, int r8r8_COMPARE:
    -1, (X1,Y1) < (X2,Y2);
     0, (X1,Y1) = (X2,Y2);
    +1, (X1,Y1) > (X2,Y2).
*/
{
	static short result = SHRT_MAX;
	
	ityp * const a_data = data; 
	ityp x1 = a_data[0];
	ityp y1 = a_data[1];
	ityp x2 = a_data[2];
	ityp y2 = a_data[3];
	
    if ( x2 < x1 || y2 < y1 )
    {
    	result = +1;
        return &result;
    }
    else if ( x1 < x2 || y1 < y2 )
    {
    	result = -1;
        return &result;
    }

    result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8r8r8_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8r8r8_COMPARE compares two r8r8r8's.
  Discussion:
    An r8r8r8 is simply 3 r8 values, stored as scalars.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp X1, Y1, Z1, the first vector.
    Input, ityp X2, Y2, Z2, the second vector.
    Output, int r8r8r8_COMPARE:
    -1, (X1,Y1,Z1) < (X2,Y2,Z2);
     0, (X1,Y1,Z1) = (X2,Y2,Z2);
    +1, (X1,Y1,Z1) > (X2,Y2,Z2).
*/
{
	static short result = SHRT_MAX;
	
	ityp * const a_data = data; 
	ityp x1 = a_data[0];
	ityp y1 = a_data[1];
	ityp z1 = a_data[2];
	ityp x2 = a_data[3];
	ityp y2 = a_data[4];
	ityp z2 = a_data[5];
	
    if ( x1 < x2 || y1 < y2 || z1 < z2 )
    {
    	result = -1;
        return &result;
    }
    else if ( x2 < x1 || y2 < y1 || z2 < z1 )
    {
    	result = 1;
        return &result;
    }
    
    result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8row_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8ROW_MAX returns the row maximums of an r8ROW.
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
    15 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in the array.
    Input, ityp A[M*N], the array to be examined.
    Output, ityp r8ROW_MAX[M], the maximums of the rows.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp *amax = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
    {
        amax[i] = a[i+0*m];

        for ( j = 1; j < n; ++j)
            if ( amax[i] < a[i+j*m] )
                amax[i] = a[i+j*m];
    }

    return amax;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8row_mean ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8ROW_MEAN returns the row means of an r8ROW.
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
    15 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in the array.
    Input, ityp A[M*N], the array to be examined.
    Output, ityp r8ROW_MEAN[M], the means, or averages, of the rows.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp *mean;

    mean = ( ityp * ) malloc ( m * sizeof ( m ) );

    for ( i = 0; i < m; ++i )
    {
        mean[i] = 0.00;
        for ( j = 0; j < n; ++j )
            mean[i] += a[i+j*m];
        mean[i] /= ( ityp ) ( n );
    }

    return mean;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8row_min ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8ROW_MIN returns the row minimums of an r8ROW.
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
    15 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns in the array.
    Input, ityp A[M*N], the array to be examined.
    Output, ityp r8ROW_MIN[M], the minimums of the rows.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp *amin = ( ityp * ) malloc ( m * sizeof ( ityp ) );

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
__MATHSUITE __JBURKARDT  void   * _r8row_sum ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8ROW_SUM returns the sums of the rows of an r8ROW.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the M by N array.
    Output, ityp ROWSUM[M], the sum of the entries of
    each row.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp *rowsum  = ( ityp * ) malloc ( m * sizeof ( ityp ) );

    for ( i = 0; i < m; ++i )
    {
        rowsum[i] = 0.00;
        for ( j = 0; j < n; ++j )
            rowsum[i] += a[i+j*m];
    }

    return rowsum;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8row_swap ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8ROW_SWAP swaps two rows of an r8ROW.
  Discussion:
    The two dimensional information is stored as a one dimensional
    array, by columns.
    The row indices are 1 based, NOT 0 based.  However, a preprocessor
    variable, called OFFSET, can be reset from 1 to 0 if you wish to
    use 0-based indices.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 September 2003
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input/output, ityp A[M*N], an array of data.
    Input, int IROW1, IROW2, the two rows to swap.
    These indices should be between 1 and M.
*/
{
	const _2dtpit2dt * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	const register dim_typ irow1 = s_data->a3;
	const register dim_typ irow2 = s_data->a4;
	
    # define OFFSET 1
    dim_typ j;
    ityp t;
    /*
    Check.
    */

    if ( irow1*irow2 == 0 || m < irow1 || m < irow2 || irow1 == irow2 )
        return NULL;

    for ( j = 0; j < n; ++j )
    {
        t              = a[irow1-1+j*m];
        a[irow1-1+j*m] = a[irow2-1+j*m];
        a[irow2-1+j*m] = t;
    }

    return NULL;
    # undef OFFSET
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8row_to_r8vec ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8ROW_TO_r8VEC converts an r8ROW into an r8VEC.
  Example:
    M = 3, N = 4
    A =
      11 12 13 14
      21 22 23 24
      31 32 33 34
    r8ROW_TO_r8VEC = ( 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the M by N array.
    Output, ityp r8ROW_TO_r8VEC[M*N], a vector containing the M rows of A.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j, k;
    ityp *x = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );

    for ( j = k = 0; j < n; ++j )
        for ( i = 0; i < m; ++i )
        {
            x[k] = a[i+j*m];
            ++ k;
        }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_01_to_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_01_TO_AB shifts and rescales data to lie within given bounds.
  Discussion:
    An r8VEC is a vector of r8's.
    On input, A contains the original data, which is presumed to lie
    between 0 and 1.  However, it is not necessary that this be so.
    On output, A has been shifted and rescaled so that all entries which
    on input lay in [0,1] now lie between AMIN and AMAX.  Other entries will
    be mapped in a corresponding way.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input/output, ityp A[N], the vector to be rescaled.
    Input, ityp AMAX, AMIN, the maximum and minimum values
    allowed for A.
*/
{
	const dt2itpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp amax = s_data->a1;
	const register ityp amin = s_data->a2;
	ityp * a = s_data->a3;
	
	
    ityp amax2;
    ityp amax3;
    ityp amin2;
    ityp amin3;
    dim_typ i;

    if ( amax == amin )
    {
        for ( i = 0; i < n; ++i )
            a[i] = amin;
        return NULL;
    }

    amax2 = MAX ( amax, amin );
    amin2 = MIN ( amax, amin );

    _MIN ( &amin3, n, a );
    _MAX ( &amax3, n, a );

    if ( amax3 != amin3 )
        for ( i = 0; i < n; ++i )
            a[i] = ( ( amax3 - a[i]         ) * amin2+ (         a[i] - amin3 ) * amax2 )/ ( amax3          - amin3 );
    else
        for ( i = 0; i < n; ++i )
            a[i] = 0.50 * ( amax2 + amin2 );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_ab_to_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_AB_TO_01 shifts and rescales data to lie within [0,1].
  Discussion:
    An r8VEC is a vector of r8's.
    On input, A contains the original data.  On output, A has been shifted
    and scaled so that all entries lie between 0 and 1.
    The formula is:
      A(I) := ( A(I) - AMIN ) / ( AMAX - AMIN )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input/output, ityp A[N], the data to be rescaled.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp amax;
    ityp amin;
    dim_typ i;

    _MAX ( &amax, n, a );
    _MIN ( &amin, n, a );

    if ( amin == amax )
        for ( i = 0; i < n; ++i )
            a[i] = 0.50;
    else
        for ( i = 0; i < n; ++i)
            a[i] = ( a[i] - amin ) / ( amax - amin );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_ab_to_cd ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_AB_TO_CD shifts and rescales data to lie within a given pair of bounds.
  Discussion:
    An r8VEC is a vector of r8's.
    The mininum entry of A is mapped to BMIN, the maximum entry
    to BMAX, and values in between are mapped linearly.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data values.
    Input, ityp A[N], the data to be remapped.
    Input, ityp BMIN, BMAX, the values to which MIN(A) and MAX(A)
    are to be assigned.
    Output, ityp r8VEC_AB_TO_CD[N], the remapped data.
*/
{
	const dt2itpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register ityp bmax = s_data->a1;
	const register ityp bmin = s_data->a2;
	ityp * a = s_data->a3;
	
    ityp amax;
    ityp amin;
    ityp *b;
    dim_typ i;

    b = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( bmax == bmin )
    {
        for ( i = 0; i < n; ++i )
            b[i] = bmin;
        return b;
    }

    _MAX ( &amax, n, a );
    _MIN ( &amin, n, a );

    if ( amin == amax )
        for ( i = 0; i < n; ++i)
            b[i] = 0.50 * ( bmax + bmin );
    else
        for ( i = 0; i < n; ++i )
            b[i] = ( ( amax - a[i]        ) * bmin+ (        a[i] - amin ) * bmax )/ ( amax        - amin );

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_amax ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_AMAX returns the maximum absolute value in an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], the array.
    Output, ityp AMAX, the value of the entry
    of largest magnitude.
*/
{
	static ityp result = MAX_VAL; 
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp amax = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        if ( amax < abs ( a[i] ) )
            amax = abs ( a[i] );
            
    result = amax;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_amax_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_AMAX_INDEX returns the index of the maximum absolute value in an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], the array.
    Output, int r8VEC_AMAX_INDEX, the index of the entry of largest magnitude.
*/
{
	static dim_typ result = USHRT_MAX; 
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp amax;
    dim_typ amax_index;

    if ( n <= 0 )
        amax_index = USHRT_MAX;
    else
    {
        amax_index = 1;
        amax = abs ( a[0] );

        for (dim_typ i = 2; i <= n; ++i)
        {
            if ( amax < abs ( a[i-1] ) )
            {
                amax_index = i;
                amax = abs ( a[i-1] );
            }
        }
    }

	result = amax_index;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_amin ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_AMIN returns the minimum absolute value in an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], the array.
    Output, ityp r8VEC_AMIN, the value of the entry
    of smallest magnitude.
*/
{
	static ityp result = MAX_VAL; 
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp amin = r8_huge;
    for (dim_typ i = 0; i < n; ++i )
        if ( abs ( a[i] ) < amin )
            amin = abs ( a[i] );

	result = amin;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_amin_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_AMIN_INDEX returns the index of the minimum absolute value in an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], the array.
    Output, int r8VEC_AMIN_INDEX, the index of the entry of smallest magnitude.
*/
{
	static dim_typ result = USHRT_MAX; 
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp amin;
    dim_typ amin_index;

    if ( n <= 0 )
        amin_index = USHRT_MAX;
    else
    {
        amin_index = 1;
        amin = abs ( a[0] );

        for (dim_typ i = 2; i <= n; ++i )
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
__MATHSUITE __JBURKARDT  void   * _r8vec_any_normal ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_ANY_NORMAL returns some normal vector to V1.
  Discussion:
    An r8VEC is a vector of r8's.
    If DIM_NUM < 2, then no normal vector can be returned.
    If V1 is the zero vector, then any unit vector will do.
    No doubt, there are better, more robust algorithms.  But I will take
    just about ANY reasonable unit vector that is normal to V1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, ityp V1[DIM_NUM], the vector.
    Output, ityp r8VEC_ANY_NORMAL[DIM_NUM], a vector that is
    normal to V2, and has unit Euclidean length.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	ityp * v1 = s_data->a1;
	
    dim_typ i, j, k;
    ityp *v2;
    ityp vj;
    ityp vk;

    if ( dim_num < 2 )
        return NULL;

    v2 = ( ityp * ) malloc ( dim_num * sizeof ( ityp ) );

    if ( r8vec_norm ( dim_num, v1 ) == 0.0 )
    {
        r8vec_zero ( dim_num, v2 );
        v2[0] = 1.0;
        return v2;
    }
    /*
    Seek the largest entry in V1, VJ = V1(J), and the
    second largest, VK = V1(K).

    Since V1 does not have zero norm, we are guaranteed that
    VJ, at least, is not zero.
    */
    j = k = USHRT_MAX;
    vj = vk = 0.00;

    for ( i = 0; i < dim_num; ++i )
        if ( abs ( vk ) < abs ( v1[i] ) || k == -1 )
        {
            if ( abs ( vj ) < abs ( v1[i] ) || j == -1 )
            {
                k = j;
                vk = vj;
                j = i;
                vj = v1[i];
            }
            else
            {
                k = i;
                vk = v1[i];
            }
        }
    /*
    Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
    will just about do the trick.
    */
    r8vec_zero ( dim_num, v2 );

    v2[j] = -vk / sqrt ( vk * vk + vj * vj );
    v2[k] =  vj / sqrt ( vk * vk + vj * vj );

    return v2;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_circular_variance ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_CIRCULAR_VARIANCE returns the circular variance of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, ityp X(N), the vector whose variance is desired.
    Output, ityp r8VEC_CIRCULAR VARIANCE, the circular variance
    of the vector entries.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;	
	
    dim_typ i;
    ityp mean;
    ityp sum_c;
    ityp sum_s;

    mean = mbase_mean ( n, x );

    sum_c = 0.00;
    for ( i = 0; i < n; ++i )
        sum_c += cos ( x[i] - mean );

    sum_s = 0.00;
    for ( i = 0; i < n; ++i )
        sum_s += sin ( x[i] - mean );

	result = 1.00 - sqrt ( sum_c * sum_c + sum_s * sum_s ) / ( ityp ) n;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_COMPARE compares two r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The lexicographic ordering is used.
  Example:
    Input:
      A1 = ( 2.0, 6.0, 2.0 )
      A2 = ( 2.0, 8.0, 12.0 )
    Output
      ISGN = -1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, ityp A[N], B[N], the vectors to be compared.
    Output, int r8VEC_COMPARE, the results of the comparison:
    -1, A is lexicographically less than B,
     0, A is equal to B,
    +1, A is lexicographically greater than B.
*/
{
	static short result = SHRT_MAX;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;

    for ( dim_typ k = 0; k < n; ++k)
        if ( a[k] < b[k] )
        {
        	result = -1;
            return &result;
        }
        else if ( b[k] < a[k] )
        {
        	result = +1;
            return &result;
        }
        
    result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_convolve_circ ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_CONVOLVE_CIRC returns the discrete circular convolution of two r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    z(1+m) = xCCy(m) = sum ( 0 <= k <= n-1 ) x(1+k) * y(1+m-k)
    Here, if the index of Y becomes nonpositive, it is "wrapped around"
    by having N added to it.
    The circular convolution is equivalent to multiplication of Y by a
    circulant matrix formed from the vector X.
  Example:
    Input:
      X = (/ 1, 2, 3, 4 /)
      Y = (/ 1, 2, 4, 8 /)
    Output:
      Circulant form:
      Z = ( 1 4 3 2 ) ( 1 )
    ( 2 1 4 3 ) ( 2 )
    ( 3 2 1 4 ) * ( 4 )
    ( 4 3 2 1 ) ( 8 )
      The formula:
      Z = (/ 1*1 + 2*8 + 3*4 + 4*2,
             1*2 + 2*1 + 3*8 + 4*4,
             1*4 + 2*2 + 3*1 + 4*8,
             1*8 + 2*4 + 3*2 + 4*1 /)
      Result:
      Z = (/ 37, 44, 43, 26 /)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vectors.
    Input, ityp X[N], Y[N], the vectors to be convolved.
    Output, ityp r8VEC_CONVOLVE_CIRC[N], the circular convolution of X and Y.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	
    dim_typ i, m;
    ityp *z = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( m = 1; m <= n; ++m )
    {
        z[m-1] = 0.00;
        for ( i = 1; i <= m; ++i )
            z[m-1] += x[i-1] * y[m-i];
        for ( i = m+1; i <= n; ++i )
            z[m-1] += x[i-1] * y[n+m-i];
    }

    return z;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_correlation ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_CORRELATION returns the correlation of two r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    If X and Y are two nonzero vectors of length N, then
      correlation = (x/||x||)' (y/||y||)
    It is the cosine of the angle between the two vectors.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vectors.
    Input, ityp X[N], Y[N], the vectors to be convolved.
    Output, ityp r8VEC_CORRELATION, the correlation of X and Y.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	
    const register ityp x_norm = r8vec_norm ( n, x );
    const register ityp y_norm = r8vec_norm ( n, y );
    
    result = x_norm == 0.00 || y_norm == 0.00 ? 0.00 : r8vec_dot_product ( n, x, y ) / x_norm / y_norm;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_covar ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_COVAR computes the covariance of two vectors.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 April 2013
  Author:
    John Burkardt.
  Parameters:
    Input, ityp X[N], Y[N], the two vectors.
    Input, int N, the dimension of the two vectors.
    Output, ityp r8VEC_COVAR, the covariance of the two vectors.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	
    dim_typ i;
    ityp value;
    ityp x_average = 0.00;
    ityp y_average = 0.00;

    for ( i = 0; i < n; ++i )
        x_average += x[i];
    x_average /= ( ityp ) ( n );
    for ( i = 0; i < n; ++i )
        y_average += y[i];
    y_average /= ( ityp ) ( n );

    value = 0.00;
    for ( i = 0; i < n; ++i )
    value += ( x[i] - x_average ) * ( y[i] - y_average );

	result = value / ( ityp ) ( n - 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_cross_product_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_CROSS_PRODUCT_2D finds the cross product of a pair of r8VEC's in 2D.
  Discussion:
    Strictly speaking, the vectors lie in the (X,Y) plane, and
    the cross product here is a vector in the Z direction.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp V1[2], V2[2], the vectors.
    Output, ityp r8VEC_CROSS_PRODUCT_2D, the Z component of the cross product
    of V1 and V2.
*/
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	
	result = v1[0] * v2[1] - v1[1] * v2[0];
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_cross_product_affine_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_CROSS_PRODUCT_AFFINE_2D finds the affine cross product in 2D.
  Discussion:
    Strictly speaking, the vectors lie in the (X,Y) plane, and
    the cross product here is a vector in the Z direction.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp V0[2], the base vector.
    Input, ityp V1[2], V2[2], the vectors.
    Output, ityp r8VEC_CROSS_PRODUCT_AFFINE_2D, the Z component of the
    cross product of V1 and V2.
*/
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * v0 = a_data[0];
	ityp * v1 = a_data[1];
	ityp * v2 = a_data[2];
	
	result = ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
        - ( v2[0] - v0[0] ) * ( v1[1] - v0[1] );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_cross_product_3d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_CROSS_PRODUCT_3D computes the cross product of two r8VEC's in 3D.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp V1[3], V2[3], the coordinates of the vectors.
    Output, ityp r8VEC_CROSS_PRODUCT_3D[3], the cross product vector.
*/
{
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	
  ityp *v3 = ( ityp * ) malloc ( 3 * sizeof ( ityp ) );
  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return v3;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_cross_product_affine_3d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_CROSS_PRODUCT_AFFINE_3D computes the affine cross product in 3D.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp V0[3], the base vector.
    Input, ityp V1[3], V2[3], the coordinates of the vectors.
    Output, ityp r8VEC_CROSS_PRODUCT_AFFINE_3D[3], the cross product vector.
*/
{
	ityp ** const a_data = data;
	ityp * v0 = a_data[0];
	ityp * v1 = a_data[1];
	ityp * v2 = a_data[2];
	
    ityp *v3 = ( ityp * ) malloc ( 3 * sizeof ( ityp ) );

    v3[0] =( v1[1] - v0[1] ) * ( v2[2] - v0[2] )- ( v2[1] - v0[1] ) * ( v1[2] - v0[2] );
    v3[1] =( v1[2] - v0[2] ) * ( v2[0] - v0[0] )- ( v2[2] - v0[2] ) * ( v1[0] - v0[0] );
    v3[2] =( v1[0] - v0[0] ) * ( v2[1] - v0[1] )- ( v2[0] - v0[0] ) * ( v1[1] - v0[1] );

    return v3;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_dif ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DIF computes coefficients for estimating the N-th derivative.
  Discussion:
    An r8VEC is a vector of r8's.
    The routine computes the N+1 coefficients for a centered finite difference
    estimate of the N-th derivative of a function.
    The estimate has the form
      FDIF(N,X) = Sum (I = 0 to N) COF(I) * F ( X(I) )
    To understand the computation of the coefficients, it is enough
    to realize that the first difference approximation is
      FDIF(1,X) = F(X+DX) - F(X-DX) ) / (2*DX)
    and that the second difference approximation can be regarded as
    the first difference approximation repeated:
      FDIF(2,X) = FDIF(1,X+DX) - FDIF(1,X-DX) / (2*DX)
         = F(X+2*DX) - 2 F(X) + F(X-2*DX) / (4*DX)
    and so on for higher order differences.
    Thus, the next thing to consider is the integer coefficients of
    the sampled values of F, which are clearly the Pascal coefficients,
    but with an alternating negative sign.  In particular, if we
    consider row I of Pascal's triangle to have entries j = 0 through I,
    then P(I,J) = P(I-1,J-1) - P(I-1,J), where P(*,-1) is taken to be 0,
    and P(0,0) = 1.
       1
      -1  1
       1 -2   1
      -1  3  -3   1
       1 -4   6  -4   1
      -1  5 -10  10  -5  1
       1 -6  15 -20  15 -6 1
    Next, note that the denominator of the approximation for the
    N-th derivative will be (2*DX)**N.
    And finally, consider the location of the N+1 sampling
    points for F:
      X-N*DX, X-(N-2)*DX, X-(N-4)*DX, ..., X+(N-4)*DX, X+(N-2*DX), X+N*DX.
    Thus, a formula for evaluating FDIF(N,X) is
      fdif = 0.0
      do i = 0, n
        xi = x + (2*i-n) * h
        fdif = fdif + cof(i) * f(xi)
      end do
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 May 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the derivative to be approximated.
    N must be 0 or greater.
    Input, ityp H, the half spacing between points.
    H must be positive.
    Output, ityp r8VEC_DIF[N+1], the coefficients needed to approximate
    the N-th derivative of a function F.
*/
{
	const dtit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp h = s_data->a1;
	
    ityp *cof;
    dim_typ i, j;

    if ( n == 0 || h <= 0.00 )
        return NULL;

    cof = ( ityp * ) malloc ( ( n + 1 ) * sizeof ( ityp ) );

    for ( i = 0; i <= n; ++i )
    {
        cof[i] = 1.00;

        for ( j = i - 1; 1 <= j; --j )
            cof[j] = - cof[j] + cof[j-1];

        if ( 0 < i )
            cof[0] *= -1;
    }

    for ( i = 0; i <= n; ++i)
        cof[i] /= pow ( 2.00 * h, n );

    return cof;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_diff_norm ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DIFF_NORM returns the L2 norm of the difference of r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The vector L2 norm is defined as:
      r8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 June 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], B[N], the vectors.
    Output, ityp r8VEC_DIFF_NORM, the L2 norm of A - B.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i)
        value += ( a[i] - b[i] ) * ( a[i] - b[i] );
        
    result = sqrt ( value );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_diff_norm_l1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DIFF_NORM_L1 returns the L1 norm of the difference of r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The vector L1 norm is defined as:
      r8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 April 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], B[N], the vectors.
    Output, ityp r8VEC_DIFF_NORM_L1, the L1 norm of A - B.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        value += abs ( a[i] - b[i] );
        
    result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_diff_norm_l2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DIFF_NORM_L2 returns the L2 norm of the difference of r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The vector L2 norm is defined as:
      r8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 June 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], B[N], the vectors.
    Output, ityp r8VEC_DIFF_NORM_L2, the L2 norm of A - B.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        value += ( a[i] - b[i] ) * ( a[i] - b[i] );

	result = sqrt ( value );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_diff_norm_li ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DIFF_NORM_LI returns the L-oo norm of the difference of r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The vector L-oo norm is defined as:
      r8VEC_NORM_LI = MAX ( 1 <= I <= N ) abs ( A(I) ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 April 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], B[N], the vectors.
    Output, ityp r8VEC_DIFF_NORM_LI, the L-oo norm of A - B.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i)
        value = MAX ( value, abs ( a[i] - b[i] ) );
        
    result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_diff_norm_squared ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DIFF_NORM_SQUARED returns the square of the L2 norm of the difference of r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The vector L2 norm is defined as:
      r8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 June 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], B[N], the vectors.
    Output, ityp r8VEC_DIFF_NORM_L2, the L2 norm of A - B.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        value += ( a[i] - b[i] ) * ( a[i] - b[i] );

    result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_distance ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DISTANCE returns the Euclidean distance between two r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int DIM_NUM, the spatial dimension.
    Input, ityp V1[DIM_NUM], V2[DIM_NUM], the vectors.
    Output, ityp r8VEC_DISTANCE, the Euclidean distance
    between the vectors.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ dim_num = s_data->a0;
	ityp * v1 = s_data->a1;
	ityp * v2 = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < dim_num; ++i )
        value += pow ( v1[i] - v2[i], 2 );

	result = sqrt ( value );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_distinct ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DISTINCT is true if the entries in an r8VEC are distinct.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, ityp X[N], the vector to be checked.
    Output, int r8VEC_DISTINCT is true if all N elements of X
    are distinct.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    dim_typ i, j;

    for ( i = 1; i <= n-1; ++i )
        for ( j = 1; j <= i - 1; ++j )
            if ( x[i] == x[j] )
            {
            	result = false;
                return &result;
            }
            
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_divide ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DIVIDE divides an r8VEC by a nonzero scalar.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input/output, ityp A[N].  On input, the vector to be scaled.
    On output, each entry has been divided by S.
    Input, ityp S, the divisor.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp s = s_data->a2;
	
    for (dim_typ i = 0; i < n; ++i )
        a[i] /= s;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_dot_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DOT_PRODUCT computes the dot product of a pair of r8VEC's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 July 2007
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, ityp A1[N], A2[N], the two vectors to be considered.
    Output, ityp r8VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        value += a1[i] * a2[i];
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_dot_product_affine ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_DOT_PRODUCT_AFFINE computes the affine dot product.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, ityp V0[N], the base vector.
    Input, ityp V1[N], V2[N], the two vectors to be considered.
    Output, ityp r8VEC_DOT_PRODUCT_AFFINE, the dot product of the vectors.
*/
{
	static ityp result = MAX_VAL;
	
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * v0 = s_data->a1;
	ityp * v1 = s_data->a2;
	ityp * v2 = s_data->a3;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        value += ( v1[i] - v0[i] ) * ( v2[i] - v0[i] );
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_eq ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_EQ is true if every pair of entries in two r8VEC's is equal.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 August 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, ityp A1[N], A2[N], two vectors to compare.
    Output, int r8VEC_EQ, is TRUE if every pair of elements A1(I) and A2(I) are equal,
    and FALSE otherwise.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;
	
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
__MATHSUITE __JBURKARDT  void *   _r8vec_even_select ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
  Discussion:
    An r8VEC is a vector of r8's.
    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / ( N - 1 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of values.
    Input, ityp XLO, XHI, the low and high values.
    Input, int IVAL, the index of the desired point.
    IVAL is normally between 1 and N, but may be any integer value.
    Output, ityp r8VEC_EVEN_SELECT, the IVAL-th of N evenly spaced values
    between XLO and XHI.
    Unless N = 1, X(1) = XLO and X(N) = XHI.
    If N = 1, then X(1) = 0.5*(XLO+XHI).
*/
{
	static ityp result = MAX_VAL;
	
	const dt2itdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp xlo = s_data->a1;
	const register ityp xhi = s_data->a2;
	const register dim_typ ival = s_data->a3;
	
    result = n == 1 ? 0.50 * ( xlo + xhi ) : ( ( ityp ) ( n - ival     ) * xlo+ ( ityp ) (     ival - 1 ) * xhi )/ ( ityp ) ( n        - 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_even2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_EVEN2 linearly interpolates new numbers into an r8VECa.
  Discussion:
    An r8VEC is a vector of r8's.
    The number of values created between two old values can vary from
    one pair of values to the next
    The interpolated values are evenly spaced.
    This routine is a generalization of r8VEC_EVEN.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int MAXVAL, the size of the XVAL array, as declared by the
    user.  MAXVAL must be large enough to hold the NVAL values computed by
    this routine.  In other words, MAXVAL must be at least equal to
    NOLD + SUM (1 <= I <= NOLD-1) NFILL(I).
    Input, int NFILL[NOLD-1], the number of values
    to be interpolated between XOLD(I) and XOLD(I+1).
    NFILL(I) does not count the endpoints.  Thus, if
    NFILL(I) is 1, there will be one new point generated
    between XOLD(I) and XOLD(I+1).
    NFILL(I) must be nonnegative.
    Input, int NOLD, the number of values XOLD,
    between which extra values are to be interpolated.
    Input, ityp XOLD[NOLD], the original vector of numbers
    between which new values are to be interpolated.
    Output, int *NVAL, the number of values computed
    in the XVAL array.
    NVAL = NOLD + SUM ( 1 <= I <= NOLD-1 ) NFILL(I)
    Output, ityp XVAL[MAXVAL].  On output, XVAL contains the
    NOLD values of XOLD, as well as the interpolated
    values, making a total of NVAL values.
*/
{
	const dtpidtpitpipit * const s_data = data;
	const register dim_typ maxval = s_data->a0;
	int * nfill = s_data->a1;
	const register dim_typ nold = s_data->a2;
	ityp * xold = s_data->a3;
	int * nval = s_data->a4;
	ityp * xval = s_data->a5;
	
    dim_typ i, j;
    int nadd;

    *nval = 1;

    for ( i = 1; i <= nold - 1; ++i )
    {
        if ( nfill[i-1] < 0 || maxval < *nval + nfill[i-1] + 1 )
            return NULL;

        nadd = nfill[i-1] + 2;

        for ( j = 1; j <= nadd; ++j )
            xval[*nval+j-2] = ( ( ityp ) ( nadd - j     ) * xold[i-1]+ ( ityp ) (        j - 1 ) * xold[i] )/ ( ityp ) ( nadd     - 1 );
        *nval +=  nfill[i-1] + 1;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_even3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_EVEN3 evenly interpolates new data into an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    This routine accepts a short vector of numbers, and returns a longer
    vector of numbers, created by interpolating new values between
    the given values.
    Between any two original values, new values are evenly interpolated.
    Over the whole vector, the new numbers are interpolated in
    such a way as to try to minimize the largest distance interval size.
    The algorithm employed is not "perfect".
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int NOLD, the number of values XOLD, between which extra
    values are to be interpolated.
    Input, int NVAL, the number of values to be computed
    in the XVAL array.  NVAL should be at least NOLD.
    Input, ityp XOLD[NOLD], the original vector of numbers
    between which new values are to be interpolated.
    Output, ityp XVAL[NVAL].  On output, XVAL contains the
    NOLD values of XOLD, as well as interpolated
    values, making a total of NVAL values.
*/
{
	const _2dt2pit * const s_data = data;
	const register dim_typ nold = s_data->a0;
	const register dim_typ nval = s_data->a1;
	ityp * xold = s_data->a2;
	ityp * xval = s_data->a3;
	
    ityp density;
    dim_typ i;
    int ival;
    dim_typ j;
    dim_typ nmaybe;
    dim_typ npts;
    dim_typ ntemp;
    dim_typ ntot;
    ityp xleni;
    ityp xlentot;
    ityp xlen = 0.00;
    for ( i = 1; i <= nold - 1; ++i )
        xlen += abs ( xold[i] - xold[i-1] );

    ntemp = nval - nold;
    density = ( ityp ) ( ntemp ) / xlen;

    ival = 1;
    ntot = 0;
    xlentot = 0.00;

    for ( i = 1; i <= nold - 1; ++i )
    {
        xleni = abs ( xold[i] - xold[i-1] );
        npts = ( int ) ( density * xleni );
        ntot += npts;
        /*
        Determine if we have enough left-over density that it should
        be changed into a point.  A better algorithm would agonize
        more over where that point should go.
        */
        xlentot = xlentot + xleni;
        nmaybe = r8_nint ( xlentot * density );

        if ( ntot < nmaybe )
        {
            npts += nmaybe - ntot;
            ntot = nmaybe;
        }
        for ( j = 1; j <= npts + 2; ++j )
            xval[ival+j-2] = ( ( ityp ) ( npts+2 - j     ) * xold[i-1]+ ( ityp ) (          j - 1 ) * xold[i] )/ ( ityp ) ( npts+2     - 1 );
        ival += npts + 1;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_expand_linear ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_EXPAND_LINEAR linearly interpolates new data into an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of input data values.
    Input, ityp X[N], the original data.
    Input, int FAT, the number of data values to interpolate
    between each pair of original data values.
    Output, ityp r8VEC_EXPAND_LINEAR[(N-1)*(FAT+1)+1], the "fattened" data.
*/
{
	const _2dtpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ fat = s_data->a1;
	ityp * x = s_data->a2;
	
    dim_typ i, j, k;
    ityp *xfat;

    xfat = ( ityp * ) malloc ( ( (n-1) * (fat+1) + 1 ) * sizeof ( ityp ) );

    k = 0;

    for ( i = 0; i < n-1; ++i )
    {
        xfat[k] = x[i];
        ++ k;

        for ( j = 1; j <= fat; ++j )
        {
            xfat[k] = ( ( ityp ) ( fat - j + 1 ) * x[i]+ ( ityp ) (       j     ) * x[i+1] )/ ( ityp ) ( fat     + 1 );
            ++ k;
        }
    }

    xfat[k] = x[n-1];
    ++ k;

    return xfat;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_first_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_FIRST_INDEX indexes the first occurrence of values in an r8VEC.
  Discussion:
    For element A(I) of the vector, FIRST_INDEX(I) is the index in A of
    the first occurrence of the value A(I).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, ityp A[N], the unsorted array to examine.
    Input, ityp TOL, a tolerance for equality.
    Output, int r8VEC_FIRST_INDEX[N], the first occurrence index.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp tol = s_data->a2;
	
    int *first_index;
    dim_typ i, j;

    first_index = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; ++i )
        first_index[i] = -1;
    for ( i = 0; i < n; ++i )
        if ( first_index[i] == -1 )
        {
            first_index[i] = i;
            for ( j = i + 1; j < n; ++j )
                if ( abs ( a[i] - a[j] ) <= tol )
                    first_index[j] = i;
        }
    return first_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_frac ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_FRAC searches for the K-th smallest entry in an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    Hoare's algorithm is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2010
  Parameters:
    Input, int N, the number of elements of A.
    Input/output, ityp A[N].
    On input, A is the array to search.
    On output, the elements of A have been somewhat rearranged.
    Input, int K, the fractile to be sought.  If K = 1, the minimum
    entry is sought.  If K = N, the maximum is sought.  Other values
    of K search for the entry which is K-th in size.  K must be at
    least 1, and no greater than N.
    Output, ityp r8VEC_FRAC, the value of the K-th fractile of A.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	ityp * a = s_data->a2;
	
    ityp frac;
    dim_typ i;
    int iryt;
    dim_typ j;
    int left;
    ityp temp;
    ityp x;

    if (n*k == 0 || n < k )
    {
    	result = MAX_VAL;
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
__MATHSUITE __JBURKARDT  void   * _r8vec_fraction ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_FRACTION returns the fraction parts of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    If we regard a real number as
      r8 = SIGN * ( WHOLE + FRACTION )
    where
      SIGN is +1 or -1,
      WHOLE is a nonnegative integer
      FRACTION is a nonnegative real number strictly less than 1,
    then this routine returns the value of FRACTION.
  Example:
     r8    r8_FRACTION
    0.00      0.00
    1.01      0.01
    2.02      0.02
   19.73      0.73
   -4.34      0.34
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of arguments.
    Input, ityp X[N], the arguments.
    Output, ityp r8_FRACTION[N], the fraction parts.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp *fraction = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        fraction[i] = abs ( x[i] ) - ( ityp ) ( ( int ) ( abs ( x[i] ) ) );
    return fraction;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_gt ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_GT == ( A1 > A2 ) for two r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The comparison is lexicographic.
    A1 > A2  <=>                              A1(1) > A2(1) or
           ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
                 ...
           ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 January 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vectors.
    Input, ityp A1[N], A2[N], the vectors to be compared.
    Output, int r8VEC_GT, is TRUE if and only if A1 > A2.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;

    for (dim_typ i = 0; i < n; ++i )
    {
        if ( a2[i] < a1[i] )
        {
        	result = true;
            return &result;
        }
        else if ( a1[i] < a2[i] )
        {
        	result = false;
            return &result;
        }
    }

    result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_heap_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_HEAP_A reorders an r8VEC into a ascending heap.
  Discussion:
    An r8VEC is a vector of r8's.
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
    31 May 2009
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the size of the input array.
    Input/output, ityp A[N].
    On input, an unsorted array.
    On output, the array has been reordered into a heap.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i;
    dim_typ ifree;
    ityp key;
    dim_typ m;
    /*
    Only nodes (N/2)-1 down to 0 can be "parent" nodes.
    */
    for ( i = (n/2)-1; 0 <= i; --i )
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
__MATHSUITE __JBURKARDT  void *   _r8vec_heap_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_HEAP_D reorders an r8VEC into a descending heap.
  Discussion:
    An r8VEC is a vector of r8's.
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
    31 May 2009
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the size of the input array.
    Input/output, ityp A[N].
    On input, an unsorted array.
    On output, the array has been reordered into a heap.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i;
    dim_typ ifree;
    ityp key;
    dim_typ m;
    /*
    Only nodes (N/2)-1 down to 0 can be "parent" nodes.
    */
    for ( i = (n/2)-1; 0 <= i; --i )
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
__MATHSUITE __JBURKARDT  void *   _i4vec_zero ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_ZERO zeroes an I4VEC.
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
__MATHSUITE __JBURKARDT  void   * _r8vec_histogram ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_HISTOGRAM histograms an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    Values between A_LO and A_HI will be histogrammed into the bins
    1 through HISTO_NUM.  Values below A_LO are counted in bin 0,
    and values greater than A_HI are counted in bin HISTO_NUM+1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, ityp A[N], the array to examine.
    Input, ityp A_LO, A_HI, the lowest and highest
    values to be histogrammed.  These values will also define the bins.
    Input, int HISTO_NUM, the number of bins to use.
    Output, int HISTO_GRAM[HISTO_NUM+2], contains the number of
    entries of A in each bin.
*/
{
	const dt2itdtpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp a_lo = s_data->a1;
	ityp a_hi = s_data->a2;
	const register dim_typ histo_num = s_data->a3;
	ityp * a = s_data->a4;
	
    ityp delta;
    int *histo_gram;
    dim_typ i, j;

    histo_gram = ( int * ) malloc ( ( histo_num + 2 ) * sizeof ( int ) );
    i4vec_zero ( histo_num+2, histo_gram );
    delta = ( a_hi - a_lo ) / ( ityp ) ( histo_num<<1 );

    for ( i = 0; i < n; ++i )
    {
        if ( a[i] < a_lo )
            ++ histo_gram[0];
        else if ( a[i] <= a_hi )
        {
            j = r8_nint (( ( a_hi -       delta - a[i]        ) * ( ityp ) ( 1         )+ (      -       delta + a[i] - a_lo ) * ( ityp ) ( histo_num ) )/ ( a_hi - 2.00 * delta        - a_lo ) );
            ++ histo_gram[j];
        }
        else if ( a_hi < a[i] )
            ++ histo_gram[histo_num+1];
    }

    return histo_gram;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_house_column ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
  Discussion:
    An r8VEC is a vector of r8's.
    The routine returns a vector V that defines a Householder
    premultiplier matrix H(V) that zeros out the subdiagonal entries of
    column K of the matrix A.
       H(V) = I - 2 * v * v'
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the order of the matrix A.
    Input, ityp A[N], column K of the matrix A.
    Input, int K, the column of the matrix to be modified.
    Output, ityp r8VEC_HOUSE_COLUMN[N], a vector of unit L2 norm which
    defines an orthogonal Householder premultiplier matrix H with the property
    that the K-th column of H*A is zero below the diagonal.
*/
{
	const _2dtpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i;
    ityp s;
    ityp *v;

    v = r8vec_zero_new ( n );

    if ( k < 1 || n <= k )
        return v;

    s = r8vec_norm_l2 ( n+1-k, a+k-1 );

    if ( s == 0.00 )
        return v;

    v[k-1] = a[k-1] + abs ( s ) * r8_sign ( a[k-1] );
    r8vec_copy ( n-k, a+k, v+k );
    s = r8vec_norm_l2 ( n-k+1, v+k-1 );

    for ( i = k-1; i < n; ++i )
        v[i] /= s;

    return v;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_i4vec_dot_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_I4VEC_DOT_PRODUCT computes the dot product of an r8VEC and an I4VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    An I4VEC is a vector of I4's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, ityp r8VEC[N], the first vector.
    Input, int I4VEC[N], the second vector.
    Output, ityp r8VEC_I4VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * r8vec = s_data->a1;
	int * i4vec = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        value += r8vec[i] * ( ityp ) ( i4vec[i] );
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_in_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_IN_01 is TRUE if the entries of an r8VEC are in the range [0,1].
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], the vector
    Output, int r8VEC_IN_01, is TRUE if every entry of A is
    between 0 and 1.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;

    for (dim_typ i = 0; i < n; ++i )
        if ( a[i] < 0.00 || 1.00 < a[i] )
        {
        	result = false;
            return &result;
        }

	result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_index_delete_all ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_DELETE_ALL deletes all occurrences of a value from an indexed sorted list.
  Discussion:
    An r8VEC is a vector of r8's.
    Note that the value of N is adjusted because of the deletions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 October 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input, ityp X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Input, ityp XVAL, the value to be sought.
    Output, int *N2, the size of the current list.
    Output, ityp X2[N2], the list.
    Output, int INDX2[N2], the sort index of the list.
*/
{
	const dtpitpiitpdtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int * indx = s_data->a2;
	const register ityp xval = s_data->a3;
	dim_typ * n2 = s_data->a4;
	ityp * x2 = s_data->a5;
	int * indx2 = s_data->a6;
	
    dim_typ equal;
    dim_typ equal1;
    dim_typ equal2;
    dim_typ get;
    dim_typ i;
    dim_typ less;
    dim_typ more;
    dim_typ put;

    if ( n < 1 )
    {
        *n2 = 0;
        return NULL;
    }

    i4vec_copy ( n, indx, indx2 );
    r8vec_copy ( n, x, x2 );
    *n2 = n;

    r8vec_index_search ( *n2, x2, indx2, xval, &less, &equal, &more );

    if ( equal == 0 )
        return NULL;

    equal1 = equal;

    for ( ; ; )
    {
        if (equal1 <= 1 ||  x2[indx2[equal1-2]-1] != xval )
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

    for ( get = 1; get <= *n2; ++get )
        if ( x2[get-1] != xval )
        {
            ++ put;
            x2[put-1] = x2[get-1];
        }
    /*
    Adjust the INDX values.
    */
    for ( equal = equal1; equal <= equal2; ++equal)
        for ( i = 1; i <= *n2; ++i )
            if ( indx2[equal-1] < indx2[i-1] )
                -- indx2[i-1];
    /*
    Discard certain INDX values.
    */
    for ( i = 0; i <= *n2 - equal2 - 1; ++i )
        indx2[equal1+i-1] = indx2[equal2+i];
    for ( i = *n2 + equal1 - equal2; i <= *n2; ++i )
        indx2[i-1] = 0;
    /*
    Adjust N.
    */
    *n2 = put;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_index_delete_dupes ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted list.
  Discussion:
    An r8VEC is a vector of r8's.
    The output quantities N2, X2, and INDX2 are computed from the
    input quantities by sorting, and eliminating duplicates.
    The output arrays should be dimensioned of size N, unless the user
    knows in advance what the value of N2 will be.
    The output arrays may be identified with the input arrays.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 October 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the input list.
    Input, ityp X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Output, int *N2, the number of unique entries in X.
    Output, ityp X2[N2], a copy of the list which has
    been sorted, and made unique.
    Output, int INDX2[N2], the sort index of the new list.
*/
{
	const dtpitpipdtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int * indx = s_data->a2;
	dim_typ * n2 = s_data->a3;
	ityp * x2 = s_data->a4;
	int * indx2 = s_data->a5;
	
    dim_typ i;
    dim_typ n3;
    ityp *x3;

    i = n3 = 0;
    x3 = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( ; ; )
    {
        ++ i;

        if ( n < i )
            break;

        if (1 < i && x[indx[i-1]-1] == x3[n3-1] )
            continue;
        ++ n3;
        x3[n3-1] = x[indx[i-1]-1];
    }
    /*
    Set the output data.
    */
    *n2 = n3;
    r8vec_copy ( n3, x3, x2 );
    for ( i = 0; i < n3; ++i )
        indx2[i] = i + 1;

    free ( x3 );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_index_delete_one ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_DELETE_ONE deletes one copy of a value from an indexed sorted list.
  Discussion:
    An r8VEC is a vector of r8's.
    If the value occurs in the list more than once, only one copy is deleted.
    Note that the value of N is adjusted because of the deletions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 October 2000
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input, ityp X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Input, ityp XVAL, the value to be sought.
    Output, int *N2, the size of the current list.
    Output, ityp X2[N2], the list.
    Output, int INDX2[N2], the sort index of the list.
*/
{
	const dtpitpiitpdtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int * indx = s_data->a2;
	const register ityp xval = s_data->a3;
	dim_typ * n2 = s_data->a4;
	ityp * x2 = s_data->a5;
	int * indx2 = s_data->a6;
	
    dim_typ equal;
    dim_typ i;
    dim_typ j;
    dim_typ less;
    dim_typ more;

    if ( n < 1 )
    {
        *n2 = 0;
        return NULL;
    }

    *n2 = n;
    i4vec_copy ( *n2, indx, indx2 );
    r8vec_copy ( *n2, x, x2 );

    r8vec_index_search ( *n2, x2, indx2, xval, &less, &equal, &more );

    if ( equal != 0 )
    {
        j = indx2[equal-1];
        for ( i = j; i <= *n2-1; ++i )
            x2[i-1] = x[i];
        for ( i = equal; i <= *n2-1; ++i )
            indx2[i-1] = indx2[i];
        for ( i = 1; i <= *n2 - 1; ++i )
            if ( j < indx2[i-1] )
                -- indx2[i-1];
        -- *n2;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_index_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_INSERT inserts a value in an indexed sorted list.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2005
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N, the size of the current list.
    Input, ityp X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Input, ityp XVAL, the value to be sought.
*/
{
	const pdtpitpiit * const s_data = data;
	dim_typ * n = s_data->a0;
	ityp * x = s_data->a1;
	int * indx = s_data->a2;
	const register ityp xval = s_data->a3;
	
    dim_typ equal;
    dim_typ i;
    dim_typ less;
    dim_typ more;

    if ( *n <= 0 )
    {
        *n = indx[0] = 1;
        x[0] = xval;
        return NULL;
    }

    r8vec_index_search ( *n, x, indx, xval, &less, &equal, &more );

    x[*n] = xval;
    for ( i = *n; more <= i; --i )
        indx[i] = indx[i-1];
    indx[more-1] = *n + 1;
    ++ *n;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_index_insert_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_INSERT_UNIQUE inserts a unique value in an indexed sorted list.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2005
  Author:
    John Burkardt
  Parameters:
    Input/output, int *N, the size of the current list.
    If the input value XVAL does not already occur in X, then N is increased.
    Input/output, ityp X[N], the list.
    If the input value XVAL does not already occur in X, then it is added
    to X.
    Input/output, int INDX[N], the sort index of the list.
    If the input value XVAL does not already occur in X, then INDX is updated.
    Input, ityp XVAL, the value which will be inserted into the X
    vector if it is not there already.
*/
{
	const pdtpitpiit * const s_data = data;
	dim_typ * n = s_data->a0;
	ityp * x = s_data->a1;
	int * indx = s_data->a2;
	const register ityp xval = s_data->a3;
	
    dim_typ equal;
    dim_typ i;
    dim_typ less;
    dim_typ more;

    if ( *n <= 0 )
    {
        *n = indx[0] = 1;
        x[0] = xval;
        return NULL;
    }
    /*  
    Does XVAL already occur in X?
    */
    r8vec_index_search ( *n, x, indx, xval, &less, &equal, &more );

    if ( equal == 0 )
    {
        x[*n] = xval;
        for ( i = *n; more <= i; --i )
            indx[i] = indx[i-1];
        indx[more-1] = *n + 1;
        ++ *n;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_index_order ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_ORDER sorts an r8VEC using an index vector.
  Discussion:
    An r8VEC is a vector of r8's.
    The index vector itself is not modified.  Therefore, the pair
 (X,INDX) no longer represents an index sorted vector.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input/output, ityp X[N], the list.  On output, the list
    has been sorted.
    Input, int INDX[N], the sort index of the list.
*/
{
	const dtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int * indx = s_data->a2;
	
    dim_typ i;
    ityp *y;

    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        y[i] = x[indx[i]-1];
        x[i] = y[i];
    }
    free ( y );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_index_search ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_SEARCH searches for a value in an indexed sorted r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input, ityp X[N], the list.
    Input, int INDX[N], the sort index of the list.
    Input, ityp XVAL, the value to be sought.
    Output, int *LESS, *EQUAL, *MORE, the indexes in INDX of the
    entries of X that are just less than, equal to, and just greater
    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
    is the greatest entry of X, then MORE is N+1.
*/
{
	const dtpitpiit3pdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int * indx = s_data->a2;
	const register ityp xval = s_data->a3;
	dim_typ * less = s_data->a4;
	dim_typ * equal = s_data->a5;
	dim_typ * more = s_data->a6;
	
    dim_typ hi;
    dim_typ lo;
    dim_typ mid;
    ityp xhi;
    ityp xlo;
    ityp xmid;

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
__MATHSUITE __JBURKARDT  void *   _r8vec_index_sort_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_SORT_UNIQUE creates a sort index for an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2005
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the current list.
    Input, ityp X[N], the list.
    Output, int *N2, the number of unique elements in X.
    Output, ityp X2[N2], a list of the unique elements of X.
    Output, int INDX2[N2], the sort index of the list.
*/
{
	const dt2pitpdtpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * x2 = s_data->a2;
	dim_typ * n2 = s_data->a3;
	int * indx2 = s_data->a4;
	
    dim_typ i;

    *n2 = 0;

    for ( i = 0; i < n; ++i )
        r8vec_index_insert_unique ( n2, x2, indx2, x[i] );

    for ( i = *n2; i < n; ++i )
        x2[i] = indx2[i] = -1;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_index_sorted_range ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEX_SORTED_RANGE: search index sorted vector for elements in a range.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items in the vector.
    Input, ityp R[N], the index sorted vector.
    Input, int INDX[N], the vector used to sort R.
    The vector R[INDX[*]] is sorted.
    Input, ityp R_LO, R_HI, the limits of the range.
    Output, int *I_LO, *I_HI, the range of indices
    so that I_LO <= I <= I_HI => R_LO <= R[INDX[I]] <= R_HI.  If no
    values in R lie in the range, then I_HI < I_LO will be returned.
*/
{
	const dtpitpi2it2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * r = s_data->a1;
	int * indx = s_data->a2;
	ityp r_lo = s_data->a3;
	ityp r_hi = s_data->a4;
	int * i_lo = s_data->a5;
	int * i_hi = s_data->a6;
	
    dim_typ i1;
    dim_typ i2;
    dim_typ j1;
    dim_typ j2;
    /*
    Cases we can handle immediately.
    */
    if ( r[indx[n-1]] < r_lo )
    {
        *i_lo = n;
        *i_hi = n - 1;
        return NULL;
    }

    if ( r_hi < r[indx[0]] )
    {
        *i_lo = 0;
        *i_hi = -1;
        return NULL;
    }
    /*
    Are there are least two intervals?
    */
    if ( n == 1 )
    {
        if ( r_lo <= r[indx[0]] && r[indx[0]] <= r_hi )
            *i_lo = *i_hi = 1;
        else
        {
            *i_lo = 0;
            *i_hi = -1;
        }
        return NULL;
    }
    /*
    Bracket R_LO.
    */
    if ( r_lo <= r[indx[0]] )
        i_lo = 0;
    else
    {
        /*
        R_LO is in one of the intervals spanned by R(INDX(J1)) to R(INDX(J2)).
        Examine the intermediate interval [R(INDX(I1)), R(INDX(I1+1))].
        Does R_LO lie here, or below or above?
        */
        j1 = 0;
        j2 = n - 1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;

        for ( ; ; )
        {
            if ( r_lo < r[indx[i1]] )
            {
                j2 = i1;
                i1 = ( j1 + j2 - 1 ) / 2;
                i2 = i1 + 1;
            }
            else if ( r[indx[i2]] < r_lo )
            {
                j1 = i2;
                i1 = ( j1 + j2 - 1 ) / 2;
                i2 = i1 + 1;
            }
            else
            {
                *i_lo = i1;
                break;
            }
        }
    }
    /*
    Bracket R_HI.
    */
    if ( r[indx[n-1]] <= r_hi )
        *i_hi = n - 1;
    else
    {
        j1 = *i_lo;
        j2 = n - 1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;

        for ( ; ; )
        {
            if ( r_hi < r[indx[i1]] )
            {
                j2 = i1;
                i1 = ( j1 + j2 - 1 ) / 2;
                i2 = i1 + 1;
            }
            else if ( r[indx[i2]] < r_hi )
            {
                j1 = i2;
                i1 = ( j1 + j2 - 1 ) / 2;
                i2 = i1 + 1;
            }
            else
            {
                *i_hi = i2;
                break;
            }
        }
    }
    /*
    We expect to have computed the largest I_LO and smallest I_HI such that
    R(INDX(I_LO)) <= R_LO <= R_HI <= R(INDX(I_HI))
    but what we want is actually
    R_LO <= R(INDX(I_LO)) <= R(INDX(I_HI)) <= R_HI
    which we can usually get simply by incrementing I_LO and decrementing I_HI.
    */
    if ( r[indx[*i_lo]] < r_lo )
    {
        ++ *i_lo;
        if ( n - 1 < *i_lo )
            *i_hi = *i_lo - 1;
    }

    if ( r_hi < r[indx[*i_hi]] )
    {
        -- *i_hi;
        if ( i_hi < 0 )
            *i_lo = *i_hi + 1;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_indexed_heap_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEXED_HEAP_D creates a descending heap from an indexed r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    An indexed r8VEC is an r8VEC of data values, and an r8VEC of N indices,
    each referencing an entry of the data vector.
    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
    we have:
      A[INDX[2*J+1]]   <= A[INDX[J]]
    and
      A[INDX[2*J+2]] <= A[INDX[J]]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2010
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
    Input, ityp A[*], the data vector.
    Input/output, int INDX[N], the index array.
    Each entry of INDX must be a valid index for the array A.
    On output, the indices have been reordered into a descending heap.
*/
{
	const dtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	int * indx = s_data->a2;
	
    dim_typ i;
    dim_typ ifree;
    int key;
    dim_typ m;
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
            m = (ifree<<1) + 1;
            /*
            Does the first position exist?
            */
            if ( n - 1 < m )
                break;
            /*
            Does the second position exist?
            */
            /*
            If both positions exist, take the larger of the two values,
            and update M if necessary.
            */
            if ( m + 1 <= n - 1 && a[indx[m]] < a[indx[m+1]] )
                ++ m;
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
__MATHSUITE __JBURKARDT  void *   _r8vec_indexed_heap_d_extract ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    An indexed r8VEC is an r8VEC of data values, and an r8VEC of N indices,
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
    18 August 2010
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
    Input, ityp A[*], the data vector.
    Input/output, int INDX[N], the index vector.
    Output, int r8VEC_INDEXED_HEAP_D_EXTRACT, the index in A of the item of
    maximum value, which has now been removed from the heap.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const pdtpitpi * const s_data = data;
	dim_typ * n = s_data->a0;
	ityp * a = s_data->a1;
	int * indx = s_data->a2;
	
    dim_typ indx_extract;

    if ( *n < 1 )
    {
    	result = USHRT_MAX;
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
    r8vec_indexed_heap_d ( *n, a, indx );
    result = indx_extract;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_indexed_heap_d_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    An indexed r8VEC is an r8VEC of data values, and an r8VEC of N indices,
    each referencing an entry of the data vector.
    Note that the argument N must be a variable, and will be incremented before
    return, and that INDX must be able to hold one more entry on output than
    it held on input.
    This is one of three functions needed to model a priority queue.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2010
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
    Input, ityp A[*], the data vector.
    Input/output, int INDX[N], the index vector.
    Input, int INDX_INSERT, the index in A of the value
    to be inserted into the heap.
*/
{
	const dtpitpdtpi * const s_data = data;
	
	const register dim_typ indx_insert = s_data->a0;
	ityp * a = s_data->a1;	
	dim_typ * n = s_data->a2;
	int * indx = s_data->a3;
	
    dim_typ i;
    dim_typ parent;

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
__MATHSUITE __JBURKARDT  void *   _r8vec_indexed_heap_d_max ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    An indexed r8VEC is an r8VEC of data values, and an r8VEC of N indices,
    each referencing an entry of the data vector.
    This is one of three functions needed to model a priority queue.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 August 2010
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
    Input, ityp A[*], the data vector.
    Input, int INDX[N], the index vector.
    Output, int r8VEC_INDEXED_HEAP_D_MAX, the index in A of the maximum value
    in the heap.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	int * indx = s_data->a2;
	
	result = indx[0];
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_indicator0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDICATOR0 sets an r8VEC to the indicator vector {0,1,2...}.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, ityp A[N], the array.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        a[i] = ( ityp ) ( i );
        
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_indicator0_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDICATOR0_NEW sets an r8VEC to the indicator vector {0,1,2...}.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, ityp r8VEC_INDICATOR0_NEW[N], the array.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for (dim_typ i = 0; i < n; ++i )
        a[i] = ( ityp ) ( i );
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_indicator1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDICATOR1 sets an r8VEC to the indicator vector {1,2,3...}.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, ityp A[N], the array.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        a[i] = ( ityp ) ( i + 1 );
        
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_indicator1_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INDICATOR1_NEW sets an r8VEC to the indicator vector {1,2,3...}.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 September 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Output, ityp r8VEC_INDICATOR1_NEW[N], the array.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp *a = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        a[i] = ( ityp ) ( i + 1 );
    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_insert ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_INSERT inserts a value into an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the array on input.
    Input/output, ityp A[N+1], the array.  On input, A is
    assumed to contain only N entries, while on output, A actually
    contains N+1 entries.
    Input, int POS, the position to be assigned the new entry.
    1 <= POS <= N+1.
    Input, ityp VALUE, the value to be inserted.
*/
{
	const _2dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ pos = s_data->a1;	
	ityp * a = s_data->a2;
	ityp value = s_data->a3;

    if ( pos < 1 || n + 1 < pos )
        return NULL;
    else
    {
        for (dim_typ i = n+1; pos+1 <= i; --i )
            a[i-1] = a[i-2];

        a[pos-1] = value;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_is_int ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_IS_INT is TRUE if an r8VEC is integral.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], the vector
    Output, int r8VEC_IS_INT, is TRUE if every entry of A is an integer.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( a[i] != ( ityp ) ( int ) a[i] )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_is_nonnegative ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_IS_NONNEGATIVE is true if all entries in an r8VEC are nonnegative.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, ityp X[N], the vector to be checked.
    Output, int r8VEC_IS_NONNEGATIVE is true if all elements of X
    are nonnegative.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( x[i] < 0.0 )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_is_zero ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_IS_ZERO is true if the entries in an r8VEC are all zero.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, ityp X[N], the vector to be checked.
    Output, int r8VEC_IS_ZERO is true if all N elements of X
    are zero.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        if ( x[i] != 0.00 )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_lt ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_LT == ( A1 < A2 ) for two r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The comparison is lexicographic.
    A1 < A2  <=>                              A1(1) < A2(1) or
           ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
                 ...
           ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 January 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vectors.
    Input, ityp A1[N], A2[N], the vectors to be compared.
    Output, int r8VEC_LT, is TRUE if and only if A1 < A2.
*/
{
	static bool result = 2;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;

    for (dim_typ i = 0; i < n; ++i)
    {
        if ( a1[i] < a2[i] )
        {
        	result = true;
            return &result;
        }
        else if ( a2[i] < a1[i] )
        {
        	result = false;
            return &result;
        }
    }

	result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_max_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_MAX_INDEX returns the index of the maximum value in an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], the array.
    Output, int r8VEC_MAX_INDEX, the index of the largest entry.
*/
{
	static short result = SHRT_MAX;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
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
__MATHSUITE __JBURKARDT  void *   _r8vec_min_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_MIN_INDEX returns the index of the minimum value in an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], the array.
    Output, int r8VEC_MIN_INDEX, the index of the smallest entry.
*/
{
	static short result = SHRT_MAX;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    short value;

    if ( n <= 0 )
        value = - 1;
    else
    {
        value = 0;
        for (dim_typ i = 1; i < n; ++i )
            if ( a[i] < a[value] )
                value = i;
    }
    
    result = value; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_min_pos ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_MIN_POS returns the minimum positive value of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 November 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries.
    Input, ityp A[N], the array.
    Output, ityp r8VEC_MIN_POS, the smallest positive entry.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp value = r8_huge;

    for (dim_typ i = 0; i < n; ++i )
        if ( 0.00 < a[i] )
            if ( a[i] < value )
                value = a[i];
                
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_mirror_next ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_MIRROR_NEXT steps through all sign variations of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    In normal use, the user would set every element of A to be positive.
    The routine will take the input value of A, and output a copy in
    which the signs of one or more entries have been changed.  Repeatedly
    calling the routine with the output from the previous call will generate
    every distinct "variation" of A; that is, all possible sign variations.
    When the output variable DONE is TRUE (or equal to 1), then the
    output value of A_NEW is the last in the series.
    Note that A may have some zero values.  The routine will essentially
    ignore such entries; more exactly, it will not stupidly assume that -0
    is a proper "variation" of 0.
    Also, it is possible to call this routine with the signs of A set
    in any way you like.  The routine will operate properly, but it
    will nonethess terminate when it reaches the value of A in which
    every nonzero entry has negative sign.
    More efficient algorithms using the Gray code seem to require internal
    memory in the routine, which is not one of MATLAB's strong points,
    or the passing back and forth of a "memory array", or the use of
    global variables, or unnatural demands on the user.  This form of
    the routine is about as clean as I can make it.
  Example:
      Input         Output
    ---------    --------------
    A            A         DONE
    ---------    --------  ----
     1  2  3     -1  2  3  false
    -1  2  3      1 -2  3  false
     1 -2  3     -1 -2  3  false
    -1 -2  3      1  2 -3  false
     1  2 -3     -1  2 -3  false
    -1  2 -3      1 -2 -3  false
     1 -2 -3     -1 -2 -3  false
    -1 -2 -3      1  2  3  true
     1  0  3     -1  0  3  false
    -1  0  3      1  0 -3  false
     1  0 -3     -1  0 -3  false
    -1  0 -3      1  0  3  true
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of entries in the vector.
    Input/output, ityp A[N], a vector of real numbers.  On
    output, some signs have been changed.
    Output, int r8VEC_MIRROR_NEXT, is TRUE if the input vector A was
    the last element
    in the series (every entry was nonpositive); the output vector is reset
    so that all entries are nonnegative, but presumably the ride is over.
*/
{
	static int result = INT_MAX;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i;
    /*
    Seek the first strictly positive entry of A.
    */
    short positive = -1;
    for ( i = 0; i < n; ++i)
        if ( 0.00 < a[i] )
        {
            positive = i;
            break;
        }
    /*
    If there is no strictly positive entry of A, there is no successor.
    */
    if ( positive == -1 )
    {
        for ( i = 0; i < n; ++i )
            a[i] *= -1;
        result = true;
        return &result;
    }
    /*
    Otherwise, negate A up to the positive entry.
    */
    for ( i = 0; i <= positive; ++i)
        a[i] *= -1;

	result = false;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_negative_strict ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NEGATIVE_STRICT: all entries of r8VEC are strictly negative.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    Input, ityp A[N], the vector.
    Output, int r8VEC_NEGATIVE_STRICT, is TRUE if every entry of
    A is strictly negative.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i)
        if ( 0 <= a[i] )
        {
        	result = false;
            return &result;
        }
        
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_nint ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NINT rounds the entries of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], the vector to be rounded.
    Output, ityp B[N], the rounded values.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp *b;
    dim_typ i;
    int s;

    b = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        s = 1 - ((a[i]<0.00)<<1);
        b[i] = ( ityp ) ( s * ( int ) ( abs ( a[i] ) + 0.50 ) );
    }

    return b;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_norm ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORM returns the L2 norm of an r8VEC.
  Discussion:
    The vector L2 norm is defined as:
      r8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], the vector whose L2 norm is desired.
    Output, ityp r8VEC_NORM, the L2 norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp v = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        v += a[i] * a[i];

	result = sqrt ( v );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_norm_affine ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORM_AFFINE returns the affine L2 norm of an r8VEC.
  Discussion:
    The affine vector L2 norm is defined as:
      r8VEC_NORM_AFFINE(V0,V1)
        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vectors.
    Input, ityp V0[N], the base vector.
    Input, ityp V1[N], the vector whose affine L2 norm is desired.
    Output, ityp r8VEC_NORM_AFFINE, the affine L2 norm of V1.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * v0 = s_data->a1;
	ityp * v1 = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i)
        value += ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
        
    result = sqrt ( value );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_norm_l1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORM_L1 returns the L1 norm of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    The vector L1 norm is defined as:
      r8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], the vector whose L1 norm is desired.
    Output, ityp r8VEC_NORM_L1, the L1 norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp v = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        v += abs ( a[i] );
        
    result = v;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_norm_l2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORM_L2 returns the L2 norm of an r8VEC.
  Discussion:
    The vector L2 norm is defined as:
      r8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 March 2003
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], the vector whose L2 norm is desired.
    Output, ityp r8VEC_NORM_L2, the L2 norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp v = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        v += a[i] * a[i];

	result = sqrt ( v );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_norm_li ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORM_LI returns the L-oo norm of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    The vector L-oo norm is defined as:
      r8VEC_NORM_LI = MAX ( 1 <= I <= N ) abs ( A(I) ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], the vector whose L-oo norm is desired.
    Output, ityp r8VEC_NORM_LI, the L-oo norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp v2;
    ityp v1 = 0.0;

    for (dim_typ i = 0; i < n; ++i)
    {
        v2 = abs ( a[i] );
        if ( v1 < v2 )
            v1 = v2;
    }

	result = v1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_norm_lp ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORM_LP returns the LP norm of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    The vector LP norm is defined as:
      r8VEC_NORM_LP = ( sum ( 1 <= I <= N ) ( abs ( A(I) ) )^P )^(1/P).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in A.
    Input, ityp A[N], the vector whose LP norm is desired.
    Input, ityp P, the index of the norm.
    Output, ityp r8VEC_NORML_LP, the LP norm of A.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp p = s_data->a2;
	
    dim_typ i;
    ityp v = 0.00;

    if ( p == 1.0 )
        for ( i = 0; i < n; ++i )
            v += abs ( a[i] );
    else if ( p == 2.0 )
    {
        for ( i = 0; i < n; ++i )
            v += a[i] * a[i];
        v = sqrt ( v );
    }
    else
    {
        for ( i = 0; i < n; ++i )
            v += pow ( abs ( a[i] ), p );
        v = pow ( ( ityp ) v, 1.0 / p );
    }

	result = v;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_normalize ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORMALIZE normalizes an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input/output, ityp A[N], the vector to be normalized.
    On output, A should have unit Euclidean norm.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i;
    ityp norm = 0.00;
    for ( i = 0; i < n; ++i)
    {
        norm += a[i] * a[i];
    }
    norm = sqrt ( norm );

    if ( norm == 0.00 )
        return NULL;

    for ( i = 0; i < n; ++i )
        a[i] /= norm;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_normalize_l1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORMALIZE_L1 normalizes an r8VEC to have unit sum.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input/output, ityp A[N], the vector to be normalized.
    On output, the entries of A should have unit sum.  However, if
    the input vector has zero sum, the routine halts.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i;
    ityp a_sum = 0.00;
    for ( i = 0; i < n; ++i )
        a_sum += a[i];

    if ( a_sum == 0.0 )
        return NULL;

    for ( i = 0; i < n; ++i )
        a[i] /= a_sum;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_normsq ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORMSQ returns the squared L2 norm of an r8VEC.
  Discussion:
    The squared vector L2 norm is defined as:
      r8VEC_NORMSQ =  sum ( 1 <= I <= N ) A(I)^2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the vector dimension.
    Input, ityp A[N], the vector.
    Output, ityp r8VEC_NORMSQ, the squared L2 norm.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp v = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        v += a[i] * a[i];
        
    result = v;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_normsq_affine ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORMSQ_AFFINE returns the sqaured affine L2 norm of an r8VEC.
  Discussion:
    The sqaured affine vector L2 norm is defined as:
      r8VEC_NORMSQ_AFFINE(V0,V1)
        = sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vectors.
    Input, ityp V0[N], the base vector.
    Input, ityp V1[N], the vector.
    Output, ityp r8VEC_NORMSQ_AFFINE, the squared affine L2 norm.
*/
{
	static ityp result = MAX_VAL;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * v0 = s_data->a1;
	ityp * v1 = s_data->a2;
	
    ityp value = 0.00;
    for (dim_typ i = 0; i < n; ++i )
        value += ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_order_type ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_ORDER_TYPE determines if an r8VEC is (non)strictly ascending/descending.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries of the array.
    Input, ityp X[N], the array to be checked.
    Output, int r8VEC_ORDER_TYPE, order indicator:
    -1, no discernable order;
    0, all entries are equal;
    1, ascending order;
    2, strictly ascending order;
    3, descending order;
    4, strictly descending order.
*/
{
	static short result = SHRT_MAX;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    short order;
    /*
    Search for the first value not equal to X(0).
    */
    dim_typ i = 0;

    for ( ; ; )
    {
        ++ i;
        if ( n-1 < i )
        {
            order = 0;
            result = order;
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
__MATHSUITE __JBURKARDT  void *   _r8vec_part_quick_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_PART_QUICK_A reorders an r8VEC as part of a quick sort.
  Discussion:
    An r8VEC is a vector of r8's.
    The routine reorders the entries of A.  Using A[0] as a
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
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries of A.
    Input/output, ityp A[N].  On input, the array to be checked.
    On output, A has been reordered as described above.
    Output, int L, R, the indices of A that define the three segments.
    Let KEY = the input value of A[0].  Then
    I <= L             A(I) < KEY;
     L < I < R         A(I) = KEY;
             R <= I    A(I) > KEY.
*/
{
	const dtpit2pdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	dim_typ * l = s_data->a2;
	dim_typ * r = s_data->a3;
	
    dim_typ i;
    ityp key;
    dim_typ m;
    ityp temp;

    if ( n < 1 )
        return NULL;
    else if ( n == 1 )
    {
        *l = 0;
        *r = 2;
        return NULL;
    }

    key = a[0];
    m = *l = 1;
    /*
    The elements of unknown size have indices between L+1 and R-1.
    */
    *r = n + 1;

    for ( i = 2; i <= n; i++ )
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
            ++ m ;
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
__MATHSUITE __JBURKARDT  void *  _r8vec_permute ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_PERMUTE permutes an r8VEC in place.
  Discussion:
    An r8VEC is a vector of r8's.
    This routine permutes an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.
  Example:
    Input:
      N = 5
      P = (   2,   4,   5,   1,   3 )
      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
      BASE = 1
    Output:
      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input, int P[N], the permutation.
    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.
    Input/output, ityp A[N], the array to be permuted.
*/
{
	const dtpiipit * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	int base = s_data->a2;
	ityp * a = s_data->a3;
	
    ityp a_temp;
    dim_typ i;
    dim_typ iget;
    dim_typ iput;
    dim_typ istart;

    if ( !perm_check ( n, p ) )
        return NULL;
    /*
    In order for the sign negation trick to work, we need to assume that the
    entries of P are strictly positive.  Presumably, the lowest number is BASE.
    So temporarily add 1-BASE to each entry to force positivity.
    */
    for ( i = 0; i < n; ++i )
        p[i] += 1 - base;
    /*
    Search for the next element of the permutation that has not been used.
    */
    for ( istart = 1; istart <= n; ++istart )
    {
        if ( p[istart-1] < 0 )
            continue;
        else if ( p[istart-1] == istart )
        {
            p[istart-1] = - p[istart-1];
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
    {
        p[i] *= -1;
        p[i] += - 1 +  base;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_permute_cyclic ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_PERMUTE_CYCLIC performs a cyclic permutation of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    For 0 <= K < N, this function cyclically permutes the input vector
    to have the form
 ( A[K], A[K+1], ..., A[N-1], A[0], ..., A[K-1] )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input, int K, the increment used.
    Input/output, ityp A[N], the array to be permuted.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register dim_typ k = s_data->a1;
	ityp * a = s_data->a2;
	
    ityp *b;
    dim_typ i;
    int ipk;

    b = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        ipk = i4_wrap ( i + k, 0, n - 1 );
        b[i] = a[ipk];
    }

    for ( i = 0; i < n; ++i)
        a[i] = b[i];

    free ( b );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_permute_uniform ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_PERMUTE_UNIFORM randomly permutes an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input/output, ityp A[N], the array to be permuted.
    Input/output, int *SEED, a seed for the random number generator.
*/
{
	const dtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	int * seed = s_data->a2;
	
    int base = 0;
    int *p;
    p = perm_uniform_new ( n, seed );
    r8vec_permute ( n, p, base, a );
    free ( p );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_polarize ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_POLARIZE decomposes an r8VEC into normal and parallel components.
  Discussion:
    An r8VEC is a vector of r8's.
    The (nonzero) vector P defines a direction.
    The vector A can be written as the sum
      A = A_normal + A_parallel
    where A_parallel is a linear multiple of P, and A_normal
    is perpendicular to P.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], the vector to be polarized
    Input, ityp P[N], the polarizing direction.
    Output, ityp A_NORMAL[N], A_PARALLEL[N], the normal
    and parallel components of A.
*/
{
	const dt4pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * p = s_data->a2;
	ityp * a_normal = s_data->a3;
	ityp * a_parallel = s_data->a4;
	
    ityp a_dot_p;
    dim_typ i;
    ityp p_norm = 0.00;

    for ( i = 0; i < n; ++i )
        p_norm += pow ( p[i], 2 );
    p_norm = sqrt ( p_norm );

    if ( p_norm == 0.00 )
    {
        for ( i = 0; i < n; i++ )
        {
            a_normal[i] = a[i];
            a_parallel[i] = 0.00;
        }
        return NULL;
    }
    a_dot_p = 0.00;
    for ( i = 0; i < n; ++i )
        a_dot_p += a[i] * p[i];
    a_dot_p /= p_norm;

    for ( i = 0; i < n; ++i)
    {
        a_parallel[i] = a_dot_p * p[i] / p_norm;
        a_normal[i] = a[i] - a_parallel[i];
    }


    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_positive_strict ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_POSITIVE_STRICT: all entries of r8VEC are strictly positive.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 June 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the dimension of the vector.
    Input, ityp A[N], the vector.
    Output, int r8VEC_POSITIVE_STRICT, is TRUE if every entry of
    A is strictly positive.
*/
{
	static bool result = 2;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    for(dim_typ i = 0; i < n; ++i )
        if ( a[i] <= 0.00 )
        {
        	result = false;
            return &result;
        }
    
    result = true;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _r8vec_range ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_RANGE finds the range of Y's within a restricted X range.
  Discussion:
    An r8VEC is a vector of r8's.
    The routine is given a set of pairs of points (X,Y), and a range
    XMIN to XMAX of valid X values.  Over this range, it seeks
    YMIN and YMAX, the minimum and maximum values of Y for
    valid X's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp X[N], the X array.
    Input, ityp XMIN, XMAX, the range of X values to check.
    Input, ityp Y[N], the Y array.
    Output, ityp *YMIN, *YMAX, the range of Y values whose
    X value is within the X range.
*/
{
	const dtpit2it3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp xmin = s_data->a2;
	ityp xmax = s_data->a3;
	ityp * y = s_data->a4;
	ityp * ymin = s_data->a5;
	ityp * ymax = s_data->a6;

    *ymin =   r8_huge;
    *ymax = - r8_huge;

    for (dim_typ i = 0; i < n; ++i )
        if ( xmin <= x[i] && x[i] <= xmax )
        {
            *ymin = MIN ( *ymin, y[i] );
            *ymax = MAX ( *ymax, y[i] );
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_range_2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_RANGE_2 updates a range to include a new r8VEC
  Discussion:
    An r8VEC is a vector of r8's.
    Given a range AMIN to AMAX, and an array A, the routine will
    decrease AMIN if necessary, or increase AMAX if necessary, so that
    every entry of A is between AMIN and AMAX.
    However, AMIN will not be increased, nor AMAX decreased.
    This routine may be used to compute the maximum and minimum of a
    collection of arrays one at a time.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], the array.
    Input/output, ityp *AMIN, *AMAX.  On input, the
    current legal range of values for A.  On output, AMIN and AMAX
    are either unchanged, or else "widened" so that all entries
    of A are within the range.
*/
{
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * amin = s_data->a2;
	ityp * amax = s_data->a3;
	
	ityp max_a, min_a;
	_MAX ( &max_a, n, a );
	_MIN ( &min_a, n, a);
    MAX ( *amax, max_a );
    MIN ( *amin, min_a );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_reverse ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_REVERSE reverses the elements of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Example:
    Input:
      N = 5, A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
    Output:
      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input/output, ityp A[N], the array to be reversed.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp temp;

    for (dim_typ i = 1; i <= n/2; ++i )
    {
        temp   = a[i-1];
        a[i-1] = a[n-i];
        a[n-i] = temp;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_rotate ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_ROTATE "rotates" the entries of an r8VEC in place.
  Discussion:
    An r8VEC is a vector of r8's.
    This routine rotates an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.
  Example:
    Input:
      N = 5, M = 2
      A    = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
    Output:
      A    = ( 4.0, 5.0, 1.0, 2.0, 3.0 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of objects.
    Input, int M, the number of positions to the right that
    each element should be moved.  Elements that shift pass position
    N "wrap around" to the beginning of the array.
    Input/output, ityp A[N], the array to be rotated.
*/
{
	const _2dtpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ iget;
    dim_typ iput;
    dim_typ istart;
    int mcopy;
    dim_typ nset;
    ityp temp;
    /*
    Force M to be positive, between 0 and N-1.
    */
    mcopy = i4_modp ( m, n );

    if ( mcopy == 0 )
        return NULL;

    istart = nset = 0;

    for ( ; ; )
    {
        ++ istart;

        if ( n < istart )
            break;

        temp = a[istart-1];
        iget = istart;
        /*
        Copy the new value into the vacated entry.
        */
        for ( ; ; )
        {
            iput = iget;

            iget -= mcopy;
            if ( iget < 1 )
                iget += n;

            if ( iget == istart )
                break;

            a[iput-1] = a[iget-1];
            ++ nset;
        }

        a[iput-1] = temp;
        ++ nset;

        if ( n <= nset )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_scalar_triple_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SCALAR_TRIPLE_PRODUCT computes the scalar triple product.
  Discussion:
    STRIPLE = V1 dot ( V2 x V3 ).
    STRIPLE is the volume of the parallelogram whose sides are
    formed by V1, V2 and V3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp V1[3], V2[3], V3[3], the three vectors.
    Output, ityp r8VEC_SCALAR_TRIPLE_PRODUCT, the scalar
    triple product.
*/
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	ityp * v3 = a_data[2];
	
	
    result = v1[0] * ( v2[1] * v3[2] - v2[2] * v3[1] )+ v1[1] * ( v2[2] * v3[0] - v2[0] * v3[2] )+ v1[2] * ( v2[0] * v3[1] - v2[1] * v3[0] );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_search_binary_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SEARCH_BINARY_A searches an ascending sorted r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    Binary search is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Algorithm 1.9,
    Combinatorial Algorithms,
    CRC Press, 1998, page 26.
  Parameters:
    Input, int N, the number of elements in the array.
    Input, ityp A[N], the array to be searched.  The array must
    be sorted in ascending order.
    Input, ityp AVAL, the value to be searched for.
    Output, int r8VEC_SEARCH_BINARY_A, the result of the search.
    -1, AVAL does not occur in the array.
    I, A(I) = AVAL.
*/
{
	static int result = INT_MAX;
	
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp aval = s_data->a2;
	
    dim_typ high;
    int indx;
    dim_typ low;
    dim_typ mid;

    indx = -1;
    low = 1;
    high = n;

    while ( low <= high )
    {
        mid = ( low + high ) / 2;

        if ( a[mid-1] == aval )
        {
            indx = mid;
            break;
        }
        else if ( a[mid-1] < aval )
            low = mid + 1;
        else if ( aval < a[mid-1] )
            high = mid - 1;
    }

	result = indx;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_shift ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SHIFT performs a shift on an r8VEC.
  Discussion:
    An r8VEC is a vector of r8 values.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2011
  Author:
    John Burkardt
  Parameters:
    Input, int SHIFT, the amount by which each entry is to
    be shifted.
    Input, int N, the length of the vector.
    Input/output, ityp X[N], the vector to be shifted.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ shift = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	
    dim_typ i;
    dim_typ ihi;
    dim_typ ilo;
    ityp *y;

    y = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        y[i] = x[i];
        x[i] = 0.00;
    }

    ilo = MAX ( 0, shift );
    ihi = MIN ( n, n + shift );

    for ( i = ilo; i < ihi; ++i )
        x[i] = y[i-shift];
    free ( y );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_shift_circular ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SHIFT_CIRCULAR performs a circular shift on an r8VEC.
  Discussion:
    An r8VEC is a vector of r8 values.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 March 2011
  Author:
    John Burkardt
  Parameters:
    Input, int SHIFT, the amount by which each entry is to
    be shifted.
    Input, int N, the length of the vector.
    Input/output, ityp X[N], the vector to be shifted.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ shift = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * x = s_data->a2;
	
    dim_typ i, j;
    ityp *y = ( ityp  * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        y[i] = x[i];
        j = i4_wrap ( i - shift, 0, n - 1 );
        x[i] = y[j];
    }

    free ( y );

    return NULL; 
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sort_bubble_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_BUBBLE_A ascending sorts an r8VEC using bubble sort.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, length of input array.
    Input/output, ityp A[N].
    On input, an unsorted array of floats.
    On output, A has been sorted.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i, j;
    ityp temp;

    for ( i = 0; i < n-1; ++i)
        for ( j = i+1; j < n; ++j )
        if ( a[j] < a[i] )
        {
            temp = a[i];
            a[i] = a[j];
            a[j] = temp;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sort_bubble_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_BUBBLE_D descending sorts an r8VEC using bubble sort.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, length of input array.
    Input/output, ityp A[N].
    On input, an unsorted array of floats.
    On output, A has been sorted.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i, j;
    ityp temp;

    for ( i = 0; i < n-1; ++i )
        for ( j = i+1; j < n; ++j )
            if ( a[i] < a[j] )
            {
                temp = a[i];
                a[i] = a[j];
                a[j] = temp;
            }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sort_heap_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_HEAP_A ascending sorts an r8VEC using heap sort.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of entries in the array.
    Input/output, ityp A[N].
    On input, the array to be sorted;
    On output, the array has been sorted.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ n1;
    ityp temp;

    if ( n <= 1 )
        return NULL;
    /*
    1: Put A into descending heap form.
    */
    r8vec_heap_d ( n, a );
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
    for ( n1 = n - 1; 2 <= n1; --n1 )
    {
        /*
        Restore the heap structure of the initial N1 entries of A.
        */
        r8vec_heap_d ( n1, a );
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
__MATHSUITE __JBURKARDT  void *   _r8vec_sort_heap_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_HEAP_D descending sorts an r8VEC using heap sort.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 May 2009
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of entries in the array.
    Input/output, ityp A[N].
    On input, the array to be sorted;
    On output, the array has been sorted.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ n1;
    ityp temp;

    if ( n <= 1 )
        return NULL;
    /*
    1: Put A into ascending heap form.
    */
    r8vec_heap_a ( n, a );
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
        r8vec_heap_a ( n1, a );
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
__MATHSUITE __JBURKARDT  void *   _r8vec_sort_heap_index_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an r8VEC
  Discussion:
    An r8VEC is a vector of r8's.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.
    Once the index array is computed, the sorting can be carried out
    "implicitly:
      a(indx(*))
    or explicitly, by the call
      r8vec_permute ( n, indx, 0, a )
    after which a(*) is sorted.
    Note that the index vector is 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], an array to be index-sorted.
    Output, int INDX[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
	const dtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	int * indx = s_data->a2;
	
    ityp aval;
    dim_typ i;
    int indxt;
    dim_typ ir;
    dim_typ j;
    dim_typ l;

    if ( n < 1 )
        return NULL;

    for ( i = 0; i < n; ++i )
        indx[i] = i;

    if ( n == 1 )
        return NULL;

    l = n / 2 + 1;
    ir = n;

    for ( ; ; )
    {
        if ( 1 < l )
        {
            l = l - 1;
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
            if (j < ir &&  a[indx[j-1]] < a[indx[j]] )
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

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sort_heap_index_a_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_HEAP_INDEX_A_NEW does an indexed heap ascending sort of an r8VEC
  Discussion:
    An r8VEC is a vector of r8's.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.
    Once the index array is computed, the sorting can be carried out
    "implicitly:
      a(indx(*))
    or explicitly, by the call
      r8vec_permute ( n, indx, 0, a )
    after which a(*) is sorted.
    Note that the index vector is 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], an array to be index-sorted.
    Output, int r8VEC_SORT_HEAP_INDEX_A_NEW[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp aval;
    dim_typ i;
    int *indx;
    int indxt;
    dim_typ ir;
    dim_typ j;
    dim_typ l;

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
__MATHSUITE __JBURKARDT  void *   _r8vec_sort_heap_index_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.
    Once the index array is computed, the sorting can be carried out
    "implicitly:
      a(indx(*))
    or explicitly, by the call
      r8vec_permute ( n, indx, 0, a )
    after which a(*) is sorted.
    Note that the index vector is 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], an array to be index-sorted.
    Output, int INDX[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
	const dtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	int * indx = s_data->a2;
	
    ityp aval;
    dim_typ i;
    int indxt;
    dim_typ ir;
    dim_typ j;
    dim_typ l;

    if ( n < 1 )
        return NULL;

    for ( i = 0; i < n; ++i )
        indx[i] = i;

    if ( n == 1 )
        return NULL;

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
                j += j;
            }
            else
                j = ir + 1;
        }

        indx[i-1] = indxt;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sort_heap_index_d_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_HEAP_INDEX_D_NEW does an indexed heap descending sort of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.
    Once the index array is computed, the sorting can be carried out
    "implicitly:
      a(indx(*))
    or explicitly, by the call
      r8vec_permute ( n, indx, 0, a )
    after which a(*) is sorted.
    Note that the index vector is 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], an array to be index-sorted.
    Output, int r8VEC_SORT_HEAP_INDEX_D_NEW[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp aval;
    dim_typ i;
    int *indx;
    int indxt;
    dim_typ ir;
    dim_typ j;
    dim_typ l;

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
__MATHSUITE __JBURKARDT  void   * _r8vec_sort_heap_mask_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_HEAP_MASK_A: indexed heap ascending sort of a masked r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    An array A is given.  An array MASK of indices into A is given.
    The routine produces a vector INDX, which is a permutation of the
    entries of MASK, so that:
      A(MASK(INDX(I)) <= A(MASK(INDX(J))
    whenever
      I <= J
    In other words, only the elements of A that are indexed by MASK
    are to be considered, and the only thing that happens is that
    a rearrangment of the indices in MASK is returned that orders the
    masked elements.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], an array to be index-sorted.
    Input, int MASK_NUM, the number of mask elements.
    Input, int MASK[MASK_NUM], the mask array.  This is
    simply a list of indices of A.  The entries of MASK should
    be unique, and each one should be between 1 and N.
    Output, int INDX[MASK_NUM], the sort index.  There are MASK_NUM
    elements of A selected by MASK.  If we want to list those elements
    in order, then the I-th element is A(MASK(INDX(I))).
*/
{
	const _2dtpipit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ mask_num = s_data->a1;
	int * mask = s_data->a2;
	ityp * a = s_data->a3;
	
    ityp aval;
    dim_typ i;
    int *indx;
    int indxt;
    dim_typ ir;
    dim_typ j;
    dim_typ l;

    if ( mask_num*n == 0 )
        return NULL;

    if ( mask_num == 1 )
    {
        indx = ( int * ) malloc ( 1 * sizeof ( int ) );
        indx[0] = 1;
        return indx;
    }

    indx = i4vec_indicator1_new ( mask_num );

    l = mask_num / 2 + 1;
    ir = mask_num;

    for ( ; ; )
    {
        if ( 1 < l )
        {
            -- l;
            indxt = indx[l-1];
            aval = a[mask[indxt-1]-1];
        }
        else
        {
            indxt = indx[ir-1];
            aval = a[mask[indxt-1]-1];
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

            if ( j < ir && a[mask[indx[j-1]-1]-1] < a[mask[indx[j]-1]-1] )
                ++ j;

            if ( aval < a[mask[indx[j-1]-1]-1] )
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
__MATHSUITE __JBURKARDT  void *  _r8vec_sort_insert_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_INSERT_A ascending sorts an r8VEC using an insertion sort.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Algorithm 1.1,
    Combinatorial Algorithms,
    CRC Press, 1998, page 11.
  Parameters:
    Input, int N, the number of items in the vector.
    N must be positive.
    Input/output, ityp A[N].
    On input, A contains data to be sorted.
    On output, the entries of A have been sorted in ascending order.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i, j;
    ityp x;

    for ( i = 1; i < n; ++i)
    {
        x = a[i];
        j = i;

        while ( 1 <= j && x < a[j-1] )
        {
            a[j] = a[j-1];
            -- j;
        }
        a[j] = x;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sort_insert_index_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_INSERT_INDEX_A ascending index sorts an r8VEC using insertion.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Reference:
    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998, page 11.
  Parameters:
    Input, int N, the number of items in the vector.
    N must be positive.
    Input, ityp A[N], the array to be sorted.
    Output, int r8VEC_SORT_INSET_INDEX_A[N], the sorted indices.  The array
    is sorted when listed from A(INDX(1)) through A(INDX(N)).
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    dim_typ i;
    int *indx;
    dim_typ j;
    ityp x;

    if ( n < 1 )
        return NULL;

    indx = i4vec_indicator1_new ( n );

    for ( i = 2; i <= n; ++i )
    {
        x = a[i-1];
        j = i - 1;

        while ( 1 <= j )
        {
            if ( a[indx[j-1]-1] <= x )
                break;

            indx[j] = indx[j-1];
            -- j;
        }
        indx[j] = i;
    }

    return indx;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sort_quick_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_QUICK_A ascending sorts an r8VEC using quick sort.
  Discussion:
    An r8VEC is a vector of r8's.
  Example:
    Input:
      N = 7
      A = ( 6, 7, 3, 2, 9, 1, 8 )
    Output:
      A = ( 1, 2, 3, 6, 7, 8, 9 )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries of A.
    Input/output, ityp A[N].  On input, the array to be sorted.
    On output, A has been reordered into ascending order.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    # define LEVEL_MAX 30

    dim_typ base;
    dim_typ l_segment;
    dim_typ level;
    dim_typ n_segment;
    int rsave[LEVEL_MAX];
    dim_typ r_segment;

    if ( n == 1 || n == 0 )
        return NULL;

    level = base = 1;
    rsave[0] = n + 1;
    n_segment = n;

    while ( 0 < n_segment )
    {
        /*
        Partition the segment.
        */
        r8vec_part_quick_a ( n_segment, a+base-1, &l_segment, &r_segment );
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
__MATHSUITE __JBURKARDT  void *   _r8vec_sort_shell_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORT_SHELL_A ascending sorts an r8VEC using Shell's sort.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input/output, ityp A[N].
    On input, an array to be sorted.
    On output, the sorted array.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    ityp asave;
    dim_typ i;
    dim_typ ifree;
    dim_typ inc;
    dim_typ ipow;
    dim_typ j;
    dim_typ k;
    dim_typ maxpow;
    dim_typ test;

    if ( n <= 1 )
        return NULL;
    /*
    Determine the smallest MAXPOW so that
    N <= ( 3**MAXPOW - 1 ) / 2
    */
    maxpow = 1;
    test = 3;

    while ( test < (n<<1) + 1 )
    {
        ++ maxpow;
        test *= 3;
    }

    if ( 1 < maxpow )
    {
        -- maxpow;
        test /= 3;
    }
    /*
    Now sort groups of size ( 3**IPOW - 1 ) / 2.
    */
    for ( ipow = maxpow; 1 <= ipow; --ipow )
    {
        inc = ( test - 1 ) / 2;
        test /= 3;
        /*
        Sort the values with indices equal to K mod INC.
        */
        for ( k = 1; k <= inc; ++k )
        {
            /*
            Insertion sort of the items with index
            INC+K, 2*INC+K, 3*INC+K, ...
            */
            for ( i = inc+k; i <= n; i = i + inc )
            {
                asave = a[i-1];
                ifree = i;
                j = i - inc;

                for ( ; ; )
                {
                    if ( j<1 || a[j-1] <= asave )
                        break;

                    ifree = j;
                    a[j+inc-1] = a[j-1];
                    j = j - inc;
                }
                a[ifree-1] = asave;
            }
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sorted_merge_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORTED_MERGE_A merges two ascending sorted r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
    The elements of A and B should be sorted in ascending order.
    The elements in the output array C will also be in ascending order,
    and unique.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int NA, the dimension of A.
    Input, ityp A[NA], the first sorted array.
    Input, int NB, the dimension of B.
    Input, ityp B[NB], the second sorted array.
    Output, int *NC, the number of entries in the merged vector.
    Output, ityp r8VEC_SORTED_MERGE_A[NC], the merged unique sorted array.
*/
{
	const pit2dtpdtpit * const s_data = data;
	
	ityp * a = s_data->a0;
	const register dim_typ nb = s_data->a1;
	const register dim_typ na = s_data->a2;
	dim_typ * nc = s_data->a3;
	ityp * b = s_data->a4;
	
    ityp *c;
    ityp *d;
    dim_typ j;
    dim_typ ja;
    dim_typ jb;
    dim_typ na2;
    dim_typ nb2;
    dim_typ nd;
    dim_typ order;

    na2 = na;
    nb2 = nb;

    ja = 0;
    jb = 0;
    *nc = 0;
    nd = 0;
    d = ( ityp * ) malloc ( ( na + nb ) * sizeof ( ityp ) );

    order = r8vec_order_type ( na2, a );

    if ( order < 0 || 2 < order )
        return NULL;

    order = r8vec_order_type ( nb2, b );

    if ( order < 0 || 2 < order )
        return NULL;

    for ( ; ; )
    {
        /*
        If we've used up all the entries of A, stick the rest of B on the end.
        */
        if ( na2 <= ja )
        {
            for ( j = 1; j <= nb2 - jb; ++j )
            {
                ++ jb;
                if ( nd == 0 )
                {
                    ++ nd;
                    d[nd-1] = b[jb-1];
                }
                else if ( d[nd-1] < b[jb-1] )
                {
                    ++ nd;
                    d[nd-1] = b[jb-1];
                }
            }
            break;
        }
        /*
        If we've used up all the entries of B, stick the rest of A on the end.
        */
        else if ( nb2 <= jb )
        {
            for ( j = 1; j <= na2 - ja; ++j )
            {
                ++ ja;
                if ( nd == 0 )
                {
                    ++ nd;
                    d[nd-1] = a[ja-1];
                }
                else if ( d[nd-1] < a[ja-1] )
                {
                    ++ nd;
                    d[nd-1] = a[ja-1];
                }
            }
            break;
        }
        /*
        Otherwise, if the next entry of A is smaller, that's our candidate.
        */
        else if ( a[ja] <= b[jb] )
        {
            ++ ja;
            if ( nd == 0 )
            {
                ++ nd;
                d[nd-1] = a[ja-1];
            }
            else if ( d[nd-1] < a[ja-1] )
            {
                ++ nd;
                d[nd-1] = a[ja-1];
            }
        }
        /*
        ...or if the next entry of B is the smaller, consider that.
        */
        else
        {
            jb = jb + 1;
            if ( nd == 0 )
            {
                ++ nd;
                d[nd-1] = b[jb-1];
            }
            else if ( d[nd-1] < b[jb-1] )
            {
                ++ nd;
                d[nd-1] = b[jb-1];
            }
        }
    }

    *nc = nd;

    c = r8vec_copy_new ( nd, d );

    free ( d );

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sorted_nearest ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORTED_NEAREST returns the nearest element in a sorted r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, ityp A[N], a sorted vector.
    Input, ityp VALUE, the value whose nearest vector entry is sought.
    Output, int r8VEC_SORTED_NEAREST, the index of the nearest
    entry in the vector.
*/
{
	static short result = SHRT_MAX;
	
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp value = s_data->a2;
	
    dim_typ hi;
    dim_typ lo;
    dim_typ mid;

    if ( n < 1 )
    {
    	result = (-1);
        return &result;
	}	

    if ( n == 1 )
    {
    	result = 1;
        return &result;
    }

    if ( a[0] < a[n-1] )
    {
    if ( value < a[0] )
    {
    	result = 1;
        return &result;
    }
    else if ( a[n-1] < value )
    {
    	result = n;
        return &result;
    }
    /*
    Seek an interval containing the value.
    */
    lo = 1;
    hi = n;

    while ( lo < hi - 1 )
    {
        mid = ( lo + hi ) / 2;

        if ( value == a[mid-1] )
        {
        	result = mid;
            return &result;
        }
        else if ( value < a[mid-1] )
            hi = mid;
        else
            lo = mid;
    }
    /*
    Take the nearest.
    */
    if ( abs ( value - a[lo-1] ) < abs ( value - a[hi-1] ) )
    {
    	result = lo;
        return &result;
    }
    else
    {
    	result = hi;
        return &result;
    }
    }
    /*
    A descending sorted vector A.
    */
    else
    {
        if ( value < a[n-1] )
        {
        	result = n;
            return &result;
        }
        else if ( a[0] < value )
        {
        	result = 1;
            return &result;
        }
        /*
        Seek an interval containing the value.
        */
        lo = n;
        hi = 1;

        while ( lo < hi - 1 )
        {
            mid = ( lo + hi ) / 2;

            if ( value == a[mid-1] )
            {
            	result = mid;
                return &result;
            }
            else if ( value < a[mid-1] )
                hi = mid;
            else
                lo = mid;
        }
        /*
        Take the nearest.
        */
        if ( abs ( value - a[lo-1] ) < abs ( value - a[hi-1] ) )
        {
        	result = lo;
            return &result;
        }
        else
        {
        	result = hi;
            return &result;
        }
    }
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sorted_range ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORTED_RANGE searches a sorted vector for elements in a range.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items in the vector.
    Input, ityp R[N], the sorted vector.
    Input, ityp R_LO, R_HI, the limits of the range.
    Output, int *I_LO, *I_HI, the range of indices
    so that I_LO <= I <= I_HI => R_LO <= R(I) <= R_HI.  If no
    values in R lie in the range, then I_HI < I_LO will be returned.
*/
{
	const dtpit2pi2it * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * r = s_data->a1;
	int * i_lo = s_data->a2;
	int * i_hi = s_data->a3;
	ityp r_lo = s_data->a4;
	ityp r_hi = s_data->a5;
	
    dim_typ i1;
    dim_typ i2;
    dim_typ j1;
    dim_typ j2;
    /*
    Cases we can handle immediately.
    */
    if ( r[n-1] < r_lo )
    {
        *i_lo = - 1;
        *i_hi = - 2;
        return NULL;
    }

    if ( r_hi < r[0] )
    {
        *i_lo = - 1;
        *i_hi = - 2;
        return NULL;
    }
    /*
    Are there are least two intervals?
    */
    if ( n == 1 )
    {
        if ( r_lo <= r[0] && r[0] <= r_hi )
            *i_lo = *i_hi = 0;
        else
        {
            *i_lo = - 1;
            *i_hi = - 2;
        }
        return NULL;
    }
    /*
    Bracket R_LO.
    */
    if ( r_lo <= r[0] )
        *i_lo = 0;
    else
    {
        /*
        R_LO is in one of the intervals spanned by R(J1) to R(J2).
        Examine the intermediate interval [R(I1), R(I1+1)].
        Does R_LO lie here, or below or above?
        */
        j1 = 0;
        j2 = n - 1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;

        for ( ; ; )
        {
            if ( r_lo < r[i1] )
            {
                j2 = i1;
                i1 = ( j1 + j2 - 1 ) / 2;
                i2 = i1 + 1;
            }
            else if ( r[i2] < r_lo )
            {
                j1 = i2;
                i1 = ( j1 + j2 - 1 ) / 2;
                i2 = i1 + 1;
            }
            else
            {
                *i_lo = i1;
                break;
            }
        }
    }
    /*
    Bracket R_HI
    */
    if ( r[n-1] <= r_hi )
        *i_hi = n - 1;
    else
    {
        j1 = *i_lo;
        j2 = n - 1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;

        for ( ; ; )
        {
            if ( r_hi < r[i1] )
            {
                j2 = i1;
                i1 = ( j1 + j2 - 1 ) / 2;
                i2 = i1 + 1;
            }
            else if ( r[i2] < r_hi )
            {
                j1 = i2;
                i1 = ( j1 + j2 - 1 ) / 2;
                i2 = i1 + 1;
            }
            else
            {
                *i_hi = i2;
                break;
            }
        }
    }
    /*
    We expect to have computed the largest I_LO and smallest I_HI such that
    R(I_LO) <= R_LO <= R_HI <= R(I_HI)
    but what we want is actually
    R_LO <= R(I_LO) <= R(I_HI) <= R_HI
    which we can usually get simply by incrementing I_LO and decrementing I_HI.
    */
    if ( r[*i_lo] < r_lo )
    {
        ++ *i_lo;
        if ( n - 1 < *i_lo )
            -- *i_hi;
    }

    if ( r_hi < r[*i_hi] )
    {
        -- *i_hi;
        if ( *i_hi < 0 )
            ++ *i_lo;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sorted_split ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORTED_SPLIT "splits" a sorted r8VEC, given a splitting value.
  Discussion:
    An r8VEC is a vector of r8's.
    Given a splitting value SPLIT, the routine seeks indices
    I_LT and I_GT so that
      A(I_LT) < SPLIT < A(I_GT),
    and if there are intermediate index values between I_LT and
    I_GT, then those entries of A are exactly equal to SPLIT.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters
    Input, int N, the number of entries in A.
    Input, ityp A[N], a sorted array.
    Input, ityp SPLIT, a value to which the entries in A are
    to be compared.
    Output, int *I_LT:
    0 if no entries are less than SPLIT;
    N if all entries are less than SPLIT;
    otherwise, the index of the last entry in A less than SPLIT.
    Output, int *I_GT:
    1 if all entries are greater than SPLIT;
    N+1 if no entries are greater than SPLIT;
    otherwise the index of the first entry in A greater than SPLIT.
*/
{
	const dtpitit2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp split = s_data->a2;
	int * i_lt = s_data->a3;
	int * i_gt = s_data->a4;
	
    dim_typ hi;
    dim_typ i;
    dim_typ lo;
    dim_typ mid;

    if ( n < 1 )
    {
        *i_lt = *i_gt = -1;
        return NULL;
    }

    if ( split < a[0] )
    {
        *i_lt = 0;
        *i_gt = 1;
        return NULL;
    }

    if ( a[n-1] < split )
    {
        *i_lt = n;
        *i_gt = n + 1;
        return NULL;
    }

    lo = 1;
    hi = n;

    for ( ; ; )
    {
        if ( lo + 1 == hi )
        {
            *i_lt = lo;
            break;
        }

        mid = ( lo + hi ) / 2;

        if ( split <= a[mid-1] )
            hi = mid;
        else
            lo = mid;
    }

    for ( i = *i_lt + 1; i <= n; ++i )
        if ( split < a[i-1] )
        {
            *i_gt = i;
            return NULL;
        }

    *i_gt = n + 1;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sorted_undex ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORTED_UNDEX returns unique sorted indexes for a sorted r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
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
      I      X      XU  Undx  Xdnu
    ----+------+------+-----+-----+
      0 | 11.0 |  11.0    0     0
      1 | 11.0 |  22.0    4     0
      2 | 11.0 |  33.0    7     0
      3 | 11.0 |  55.0    8     0
      4 | 22.0 |                1
      5 | 22.0 |                1
      6 | 22.0 |                1
      7 | 33.0 |                2
      8 | 55.0 |                3
    INDX(2) = 3 means that sorted item(2) is X(3).
    XDNI(2) = 5 means that X(2) is sorted item(5).
    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
    XU(XDNU(I))) = X(I).
    XU(I)        = X(UNDX(I)).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int X_NUM, the number of data values.
    Input, ityp X_VAL[X_NUM], the data values.
    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.
    Input, ityp TOL, a tolerance for equality.
    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.
    Output, int XDNU[X_NUM], the XDNU vector.
*/
{
	const dtpitdtit2pi * const s_data = data;
	const register dim_typ x_num = s_data->a0;
	ityp * x_val = s_data->a1;
	const register dim_typ x_unique_num = s_data->a2;
	ityp tol = s_data->a3;
	int * undx = s_data->a4;
	int * xdnu = s_data->a5;
	
    dim_typ i, j;
    /*
    Walk through the sorted array.
    */
    i = j = 0;
    undx[j] = i;
    xdnu[i] = j;

    for ( i = 1; i < x_num; ++i)
    {
        if ( tol < abs ( x_val[i] - x_val[undx[j]] ) )
        {
            ++ j;
            undx[j] = i;
        }
        xdnu[i] = j;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_sorted_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORTED_UNIQUE finds the unique elements in a sorted r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    If the data is not sorted, the results of the routine will
    be garbage.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, ityp A[N], the sorted array of N elements;
    Input, ityp TOL, a tolerance for checking equality.
    Output, int *UNIQUE_NUM, the number of unique elements of A.
    Output, ityp r8VEC_SORTED_UNIQUE[UNIQUE_NUM], the unique elements of A.
*/
{
	const dtpititpdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp tol = s_data->a2;
	dim_typ * unique_num = s_data->a3;
	
    ityp *a_unique;
    dim_typ i;
    dim_typ iuniq;

    *unique_num = 0;

    if ( n <= 0 )
        return NULL;
    /*
    Determine the number of unique elements.
    */
    iuniq = 0;
    *unique_num = 1;

    for ( i = 1; i < n; ++i )
        if ( tol < abs ( a[i] - a[iuniq] ) )
        {
            iuniq = i;
            ++ *unique_num;
        }
    /*
    Set aside space for the unique elements.
    */
    a_unique = ( ityp * ) malloc ( *unique_num * sizeof ( ityp ) );
    /*
    Repeat the search, but now store the unique elements.
    */
    *unique_num = 0;

    a_unique[*unique_num] = a[0];
    *unique_num = 1;

    for ( i = 1; i < n; ++i )
        if ( tol < abs ( a[i] - a_unique[*unique_num-1] ) )
        {
            a_unique[*unique_num] = a[i];
            ++ *unique_num;
        }

    return a_unique;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sorted_unique_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    Because the array is sorted, this algorithm is O(N).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, ityp A[N], the sorted array to examine.
    Input, ityp TOL, a tolerance for checking equality.
    Output, int r8VEC_SORTED_UNIQUE_COUNT, the number of unique elements of A.
*/
{
	static dim_typ result = USHRT_MAX;

	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp tol = s_data->a2;
	
    dim_typ i;
    dim_typ unique_num;

    if ( n < 1 )
    {
    	result = 0;
        return &result;
    }

    unique_num = 1;

    for ( i = 1; i < n; ++i)
        if ( tol < abs ( a[i-1] - a[i] ) )
            ++ unique_num;

	result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_sorted_unique_hist ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SORTED_UNIQUE_HIST histograms unique elements of a sorted r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, ityp A[N], the array to examine, which must have been
    sorted.
    Input, ityp TOL, a tolerance for checking equality.
    Input, int MAXUNIQ, the maximum number of unique elements
    that can be handled.  If there are more than MAXUNIQ unique
    elements in A, the excess will be ignored.
    Output, int *UNIQUE_NUM, the number of unique elements of A.
    Output, ityp AUNIQ[UNIQUE_NUM], the unique elements of A.
    Output, int ACOUNT[UNIQUE_NUM], the number of times each element
    of AUNIQ occurs in A.
*/
{
	const dtpititdtpdtpitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp tol = s_data->a2;
	const register dim_typ maxuniq = s_data->a3;
	dim_typ * unique_num = s_data->a4;
	ityp * auniq = s_data->a5;
	int * acount = s_data->a6;
	
    dim_typ i;
    dim_typ index;
    /*
    Start taking statistics.
    */

    for ( i = 0; i < n; ++i )
    {

        if ( i == 0 )
        {
            index = 0;
            auniq[index] = a[0];
            acount[index] = 1;
        }
        else if ( abs ( a[i] - auniq[index] ) <= tol )
            ++ acount[index];
        else if ( index + 1 < maxuniq )
        {
            ++ index;
            auniq[index] = a[i];
            acount[index] = 1;
        }
    }

    *unique_num = index + 1;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT void *   _r8vec_split ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SPLIT "splits" an unsorted r8VEC based on a splitting value.
  Discussion:
    An r8VEC is a vector of r8's.
    If the vector is already sorted, it is simpler to do a binary search
    on the data than to call this routine.
    The vector is not assumed to be sorted before input, and is not
    sorted during processing.  If sorting is not needed, then it is
    more efficient to use this routine.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 August 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input/output, ityp A[N], the array to split.  On output,
    all the entries of A that are less than or equal to SPLIT
    are in A(1:ISPLIT).
    Input, ityp SPLIT, the value used to split the vector.
    It is not necessary that any value of A actually equal SPLIT.
    Output, int r8VEC_SPLIT, indicates the position of the last
    entry of the split vector that is less than or equal to SPLIT.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp split = s_data->a2;
	
    dim_typ i;
    dim_typ i1;
    dim_typ i2;
    dim_typ i3;
    dim_typ isplit;
    dim_typ j1;
    dim_typ j2;
    dim_typ j3;
    ityp temp;
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
    for ( i = 1; i <= n; ++i )
        if ( a[i2-1] <= split )
        {
            ++ i2;
            ++ j1;
        }
        else
        {
            temp = a[i2-1];
            a[i2-1] = a[i3-2];
            a[i3-2] = temp;
            -- i3;
            -- j2;
        }

	result = j1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_stutter ( void * data)
/******************************************************************************/
/*

  Purpose:
    r8VEC_STUTTER makes a "stuttering" copy of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 March 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the input vector.
    Input, ityp A[N], the vector.
    Input, int M, the "stuttering factor".
    Output, ityp AM[M*N], the stuttering vector.
*/
{
	const _2dt2pit * const s_data = data; 
	
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * a = s_data->a2;
	ityp * am = s_data->a3;
	
    dim_typ i, j, k;

    for ( i = k = 0; i < n; ++i )
        for ( j = 0; j < m; ++j )
        {
            am[k] = a[i];
            ++ k;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_stutter_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_STUTTER_NEW makes a "stuttering" copy of an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).
    This code is distributed under the GNU LGPL license.
  Modified:
    28 March 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the size of the input vector.
    Input, ityp A[N], the vector.
    Input, int M, the "stuttering factor".
    Output, ityp r8VEC_STUTTER_NEW[M*N], the stuttering vector.
*/
{
	const _2dtpit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ m = s_data->a1;
	ityp * a = s_data->a2;
	
	
    dim_typ i, j, k;
    ityp *am = ( ityp * ) malloc ( m * n * sizeof ( ityp ) );


    for ( i = k = 0; i < n; ++i)
        for ( j = 0; j < m; ++j )
        {
            am[k] = a[i];
            ++ k;
        }
    return am;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_swap ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_SWAP swaps the entries of two r8VEC's.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 January 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the arrays.
    Input/output, ityp A1[N], A2[N], the vectors to swap.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;
	
    ityp temp;
    for (dim_typ i = 0; i < n; ++i )
    {
        temp  = a1[i];
        a1[i] = a2[i];
        a2[i] = temp;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_undex ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_UNDEX returns unique sorted indexes for an r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
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
      I     X  Indx  Xdni       XU  Undx  Xdnu
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
    INDX(2) = 3 means that sorted item(2) is X(3).
    XDNI(2) = 5 means that X(2) is sorted item(5).
    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
    XU(XDNU(I))) = X(I).
    XU(I)        = X(UNDX(I)).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, int X_NUM, the number of data values.
    Input, ityp X_VAL[X_NUM], the data values.
    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.
    Input, ityp TOL, a tolerance for equality.
    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.
    Output, int XDNU[X_NUM], the XDNU vector.
*/
{
	const dtpitdtit2pi * const s_data = data;
	const register dim_typ x_num = s_data->a0;
	ityp * x_val = s_data->a1;
	const register dim_typ x_unique_num = s_data->a2;
	ityp tol = s_data->a3;
	int * undx = s_data->a4;
	int * xdnu = s_data->a5;
	
    dim_typ i;
    int *indx;
    dim_typ j;
    /*
    Implicitly sort the array.
    */
    indx = r8vec_sort_heap_index_a_new ( x_num, x_val );
    /*
    Walk through the implicitly sorted array X.
    */
    i = j = 0;
    undx[j] = indx[i];
    xdnu[indx[i]] = j;

    for ( i = 1; i < x_num; ++i )
    {
        if ( tol < abs ( x_val[indx[i]] - x_val[undx[j]] ) )
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
__MATHSUITE __JBURKARDT  void *   _r8vec_uniform_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_UNIFORM_AB returns a scaled pseudorandom r8VEC.
  Discussion:
    This routine implements the recursion
      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )
    The integer arithmetic never requires more than 32 bits,
    including a sign bit.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2005
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
    Input, int N, the number of entries in the vector.
    Input, ityp B, C, the lower and upper limits of the pseudorandom values.
    Input/output, int *SEED, a seed for the random number generator.
    Output, ityp R[N], the vector of pseudorandom values.
*/
{
	const dt2itpipit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp b = s_data->a1;
	const register ityp c = s_data->a2;
	int * seed = s_data->a3;
	ityp * r = s_data->a4;
	
    dim_typ i, k;

    if ( *seed == 0 )
        return NULL;

    for ( i = 0; i < n; ++i )
    {
        k = *seed / 127773;
        *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
        if ( *seed < 0 )
            *seed += i4_huge;
        r[i] = b + ( c - b ) * ( ityp ) ( *seed ) * 4.656612875E-10;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_unique_count ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_UNIQUE_COUNT counts the unique elements in an unsorted r8VEC.
  Discussion:
    An r8VEC is a vector of r8's.
    Because the array is unsorted, this algorithm is O(N^2).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 April 2004
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, ityp A[N], the array to examine, which does NOT have to
    be sorted.
    Input, ityp TOL, a tolerance for checking equality.
    Output, int r8VEC_UNIQUE_COUNT, the number of unique elements of A.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp tol = s_data->a2;
	
    dim_typ i, j;
    dim_typ unique_num = 0;

    for ( i = 0; i < n; ++i )
    {
        ++ unique_num;

        for ( j = 0; j < i; ++j )
            if ( abs ( a[i] - a[j] ) <= tol )
            {
                -- unique_num;
                break;
            }
    }
    
    result = unique_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_unique_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_UNIQUE_INDEX indexes the unique occurrence of values in an r8VEC.
  Discussion:
    For element A(I) of the vector, UNIQUE_INDEX(I) is the uniqueness index
    of A(I).  That is, if A_UNIQUE contains the unique elements of A,
    gathered in order, then
      A_UNIQUE ( UNIQUE_INDEX(I) ) = A(I)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of elements of A.
    Input, ityp A[N], the unsorted array to examine.
    Input, ityp TOL, a tolerance for equality.
    Output, int r8VEC_UNIQUE_INDEX[N], the unique index.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp tol = s_data->a2;
	
    dim_typ i, j;
    int *unique_index;
    dim_typ unique_num;

    unique_index = ( int * ) malloc ( n * sizeof ( int ) );

    for ( i = 0; i < n; i++ )
        unique_index[i] = -1;
    unique_num = 0;

    for ( i = 0; i < n; ++i )
        if ( unique_index[i] == -1 )
        {
            unique_index[i] = unique_num;
            for ( j = i + 1; j < n; ++j)
                if ( abs ( a[i] - a[j] ) <= tol )
                    unique_index[j] = unique_num;
            ++ unique_num;
        }
    return unique_index;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_vector_triple_product ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_VECTOR_TRIPLE_PRODUCT computes the vector triple product.
  Discussion:
    VTRIPLE = V1 x (V2 x V3)
    VTRIPLE is a vector perpendicular to V1, lying in the plane
    spanned by V2 and V3.  The norm of VTRIPLE is the product
    of the norms of V1, V2 and V3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 October 2010
  Author:
    John Burkardt
  Parameters:
    Input, ityp V1[3], V2[3], V3[3], the three vectors.
    Output, ityp r8VEC_VECTOR_TRIPLE_PRODUCT[3], the vector triple product.
*/
{
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	ityp * v3 = a_data[2];
	
    ityp *v23 = r8vec_cross_product_3d ( v2, v3 );
    ityp *v123 = r8vec_cross_product_3d ( v1, v23 );
    free ( v23 );
    return v123;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec_zero ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_ZERO zeroes an r8VEC.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 August 2008
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Output, ityp A[N], a vector of zeroes.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	
    for (dim_typ i = 0; i < n; ++i )
        a[i] = 0.00;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec2_compare ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC2_COMPARE compares two elements of an r8VEC2.
  Discussion:
    An r8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of data items.
    Input, ityp A1[N], A2[N], contain the two components of each item.
    Input, int I, J, the items to be compared.  These values will be
    1-based indices for the arrays A1 and A2.
    Output, int r8VEC2_COMPARE, the results of the comparison:
    -1, item I < item J,
     0, item I = item J,
    +1, item J < item I.
*/
{
	static short result = SHRT_MAX;
	
	const _3dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ i = s_data->a1;
	const register dim_typ j = s_data->a2;
	ityp * a1 = s_data->a3;
	ityp * a2 = s_data->a4;
	
	
    if ( a1[i-1] < a1[j-1] || a2[i-1] < a2[j-1])
    {
    	result = -1;
        return &result;
    }
    else if ( a1[i-1] == a1[j-1] || a2[j-1] < a2[i-1] || a1[j-1] < a1[i-1] )
    {
    	result = +1;
        return &result;
    }
    
    result = 0;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec2_sort_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC2_SORT_A ascending sorts an r8VEC2.
  Discussion:
    An r8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.
    Each item to be sorted is a pair of reals (X,Y), with the X
    and Y values stored in separate vectors A1 and A2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items of data.
    Input/output, ityp A1[N], A2[N], the data to be sorted.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    ityp temp;
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
        isgn = r8vec2_compare ( n, a1, a2, i, j );
    else if ( indx == 0 )
        break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec2_sort_d ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC2_SORT_D descending sorts an r8VEC2.
  Discussion:
    An r8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.
    Each item to be sorted is a pair of reals (X,Y), with the X
    and Y values stored in separate vectors A1 and A2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items of data.
    Input/output, ityp A1[N], A2[N], the data to be sorted.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;
	
    int i, j;
    int indx;
    short isgn;
    ityp temp;
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
            isgn = - r8vec2_compare ( n, a1, a2, i, j );
        else if ( indx == 0 )
            break;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec2_sort_heap_index_a ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an r8VEC2.
  Discussion:
    An r8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.
    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.
 ( X(I), Y(I) ) < ( X(J), Y(J) ) if:
    * X(I) < X(J), or
    * X(I) = X(J), and Y(I) < Y(J).
    Once the index array is computed, the sorting can be carried out
    implicitly:
  ( x(indx(*)), y(indx(*) )
    or explicitly, by the calls
      r8vec_permute ( n, indx, 0, x )
      r8vec_permute ( n, indx, 0, y )
    after which ( x(*), y(*) ), is sorted.
    Note that the index vector is 0-based.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp X[N], Y[N], pairs of X, Y coordinates of points.
    Output, int INDX[N], the sort index.  The
    I-th element of the sorted array has coordinates
 ( X(INDX(I)), Y(INDX(I) ).
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * y = s_data->a2;
	
    dim_typ i;
    int *indx;
    dim_typ indxt;
    dim_typ ir;
    dim_typ j;
    dim_typ l;
    ityp xval;
    ityp yval;

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
            xval = x[indxt];
            yval = y[indxt];
        }
        else
        {
            indxt = indx[ir-1];
            xval = x[indxt];
            yval = y[indxt];
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
            if ( j < ir && x[indx[j-1]] < x[indx[j]] ||( x[indx[j-1]] == x[indx[j]] && y[indx[j-1]] < y[indx[j]] ) )
                ++ j;

            if ( xval < x[indx[j-1]] ||( xval == x[indx[j-1]] && yval < y[indx[j-1]] ) )
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
__MATHSUITE __JBURKARDT  void *   _r8vec2_sorted_unique ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC2_SORTED_UNIQUE keeps the unique elements in an r8VEC2.
  Discussion:
    An r8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.
    Item I is stored as the pair A1(I), A2(I).
    The items must have been sorted, or at least it must be the
    case that equal items are stored in adjacent vector locations.
    If the items were not sorted, then this routine will only
    replace a string of equal values by a single representative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items.
    Input/output, ityp A1[N], A2[N].
    On input, the array of N items.
    On output, an array of UNIQUE_NUM unique items.
    Output, int *UNIQUE_NUM, the number of unique items.
*/
{
	const dt2pitpdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;
	dim_typ * unique_num = s_data->a3;
	
    dim_typ itest;
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
__MATHSUITE __JBURKARDT  void *   _r8vec2_sorted_unique_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC2_SORTED_UNIQUE_INDEX indexes unique elements in a sorted r8VEC2.
  Discussion:
    An r8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.
    Item I is stored as the pair A1(I), A2(I).
    The items must have been sorted, or at least it should be the
    case that equal items are stored in adjacent vector locations.
    If the items are not sorted, then this routine will only
    replace a string of equal values by a single representative.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 March 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of items.
    Input/output, ityp A1[N], A2[N].
    On input, the array of N items.
    On output, an array of unique items.
    Output, int *UNIQUE_NUM, the number of unique items.
    Output, int INDX[N], contains in entries 1 through UNIQUE_NUM an index
    array of the unique items.  To build new arrays with no repeated elements:
      B1(*) = A1(INDX(*))
*/
{
	const dt2pitpdtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a1 = s_data->a1;
	ityp * a2 = s_data->a2;
	dim_typ * unique_num = s_data->a3;
	int * indx = s_data->a4;
	
    dim_typ itest;

    if ( n <= 0 )
    {
        *unique_num = 0;
        return NULL;
    }

    i4vec_zero ( n, indx );
    *unique_num = indx[0] = 1;

    for ( itest = 2; itest <= n; ++itest )
        if ( a1[itest-2] != a1[itest-1] || a2[itest-2] != a2[itest-1] )
        {
            ++ *unique_num;
            indx[*unique_num-1] = itest;
        }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8vec2_sum_max_index ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC2_SUM_MAX_INDEX returns the index of the maximum sum of two r8VEC's.
  Discussion:
    An r8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 Mach 2009
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the array.
    Input, ityp A[N], B[N], two arrays whose sum
    is to be examined.
    Output, int r8VEC2_SUM_MAX_INDEX, the index of the largest entry in A+B.
*/
{
	static short result = SHRT_MAX;
	
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	
    dim_typ i;
    ityp sum_max;
    short sum_max_index;

    if ( n <= 0 )
        sum_max_index = -1;
    else
    {
        sum_max_index = 1;
        sum_max = a[0] + b[0];

        for ( i = 2; i <= n; ++i )
            if ( sum_max < a[i-1] + b[i-1] )
            {
                sum_max = a[i-1] + b[i-1];
                sum_max_index = i;
            }
    }
    
    result = sum_max_index;
    return &result;
}

#endif
