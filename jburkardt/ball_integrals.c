#ifndef __DISABLEDEEP_BALLINTEGRALS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _i4vec_uniform_ab ( void * data)
/******************************************************************************/
/*
  Purpose:
    I4VEC_UNIFORM_AB returns a scaled pseudorandom I4VEC.
  Discussion:
    The pseudorandom numbers should be uniformly distributed
    between A and B.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 January 2014
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
    Input, integer N, the dimension of the vector.
    Input, int A, B, the limits of the interval.
    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.
    Output, int X[N], a vector of random values between A and B.
*/
{
	const dt2i2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int a = s_data->a1;
	int b = s_data->a2;
	int * seed = s_data->a3;
	int * x = s_data->a4;
	
    dim_typ c;
    dim_typ i;
    dim_typ k;
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

    for ( i = 0; i < n; ++i )
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
        x[i] = value;
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gamma ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAMMA evaluates Gamma(X) for a real argument.
  Discussion:
    The C math library includes the GAMMA ( X ) function which should generally
    be used instead of this function.
    This routine calculates the gamma function for a real argument X.
    Computation is based on an algorithm outlined in reference 1.
    The program uses rational functions that approximate the gamma
    function to at least 20 significant decimal digits.  Coefficients
    for the approximation over the interval (1,2) are unpublished.
    Those for the approximation for 12 <= X are from reference 2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 January 2010
  Author:
    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.
  Reference:
    William Cody,
    An Overview of Software Development for Special Functions,
    in Numerical Analysis Dundee, 1975,
    edited by GA Watson,
    Lecture Notes in Mathematics 506,
    Springer, 1976.
    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.
  Parameters:
    Input, double X, the argument of the function.
    Output, double R8_GAMMA, the value of the function.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	ityp c[7] =
	{
		-1.910444077728E-03,
		8.4171387781295E-04,
		-5.952379913043012E-04,
		7.93650793500350248E-04,
		-2.777777777777681622553E-03,
		8.333333333333333331554247E-02,
		5.7083835261E-03
	};

	ityp eps = 2.22E-16;
	ityp fact = 1.00;
	dim_typ i, n = 0;
	ityp p[8] =
	{
		-1.71618513886549492533811E+00,
		2.47656508055759199108314E+01,
		-3.79804256470945635097577E+02,
		6.29331155312818442661052E+02,
		8.66966202790413211295064E+02,
		-3.14512729688483675254357E+04,
		-3.61444134186911729807069E+04,
		6.64561438202405440627855E+04
	};

	bool parity = false;

	ityp q[8] =
	{
		-3.08402300119738975254353E+01,
		3.15350626979604161529144E+02,
		-1.01515636749021914166146E+03,
		-3.10777167157231109440444E+03,
		2.25381184209801510330112E+04,
		4.75584627752788110767815E+03,
		-1.34659959864969306392456E+05,
		-1.15132259675553483497211E+05
	};
	ityp res;

	const ityp sqrtpi = 0.9189385332046727417803297;
	ityp sum;
	ityp value;
	ityp xbig = 171.624;
	ityp xden;
	ityp xinf = 1.79E+308;
	ityp xminin = 2.23E-308;
	ityp xnum;
	ityp y = x;
	ityp y1;
	ityp ysq;
	ityp z;
	/*
	Argument is negative.
	*/
	if ( y <= 0.00)
	{
		y = - x;
		y1 = (ityp) (dim_typ) ( y );
		res = y - y1;

		if (res != 0.00)
		{
			if (y1 != (ityp) (dim_typ) ( y1 * 0.5 ) * 2.00 )
				parity = 1;

			fact = - M_PI / sin (M_PI*res);
			y = y + 1.00;
		}
		else
		{
			res = xinf;
			value = res;
			result = value;
			return &result;
		}
	}
	/*
	Argument is positive.
	*/
	if ( y < eps )
	{
		/*
		Argument < EPS.
		*/
		if(xminin <= y)
			res = 1.00/y;
		else
		{
			res = xinf;
			value = res;
			result = value;
			return &result;
		}
	}
	else if(y < 12.00)
	{
		y1 = y;
		/*
		0.0 < argument < 1.0.
		*/
		if (y < 1.00)
		{
			z = y;
			y = y + 1.00;
		}
		/*
		1.0 < argument < 12.0.
		Reduce argument if necessary.
		*/
		else
		{
			n = (dim_typ)(y)-1;
			y = y - (ityp) ( n );
			z = y - 1.00;
		}
		/*
		Evaluate approximation for 1.0 < argument < 2.0.
		*/
		xnum = 0.00;
		xden = 1.00;
		for ( i = 00; i < 8; ++i )
		{
			xnum = ( xnum + p[i] ) * z;
			xden = xden * z + q[i];
		}
		res = xnum / xden + 1.00;
		/*
		Adjust result for case  0.0 < argument < 1.0.
		*/
		if ( y1 < y )
			res /= y1;
		/*
		Adjust result for case 2.0 < argument < 12.0.
		*/
		else if ( y < y1 )
		{
			for ( i = 1; i <= n; ++i )
			{
				res = res * y;
				y = y + 1.00;
			}
		}
	}
	else
	{
		/*
		Evaluate for 12.0 <= argument.
		*/
		if ( y <= xbig )
		{
			ysq = y * y;
			sum = c[6];
			for (i = 0; i < 6; ++i)
				sum = sum / ysq + c[i];

			sum = sum / y - y + sqrtpi;
			sum = sum + ( y - 0.5 ) * log ( y );
			res = exp ( sum );
		}
		else
		{
			res = xinf;
			value = res;
			result = value;
			return &result;
		}
	}
	/*
	Final adjustments and return.
	*/
	if (parity)
		res *= -1;

	result = fact != 1.00 ? (fact/res) : res;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_uniform_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_UNIFORM_01 returns a real pseudorandom r8.
  Discussion:
    This routine implements the recursion
      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )
    The integer arithmetic never requires more than 32 bits,
    including a sign bit.
    If the initial seed is 12345, then the first three computations are
      Input     Output      r8_UNIFORM_01
      SEED      SEED
         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 November 2004
  Author:
    John Burkardt
  Reference:
    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.
    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.
    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.
    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.
  Parameters:
    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.
    Output, ityp r8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
	static ityp result = MAX_VAL;
	
	int * seed = (int *) data;
	
    int k;

    k = *seed / 127773;
    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
        *seed += 2147483647;
    /*
    Although SEED can be represented exactly as a 32 bit integer,
    it generally cannot be represented exactly as a 32 bit real number!
    */

	result = ( ityp ) ( *seed ) * 4.656612875E-10;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_normal_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_NORMAL_01 returns a unit pseudonormal r8VEC.
  Discussion:
    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.
    This routine can generate a vector of values on one call.  It
    has the feature that it should provide the same results
    in the same order no matter how we break up the task.
    Before calling this routine, the user may call RANDOM_SEED
    in order to set the seed of the random number generator.
    The Box-Muller method is used, which is efficient, but
    generates an even number of values each time.  On any call
    to this routine, an even number of new values are generated.
    Depending on the situation, one value may be left over.
    In that case, it is saved for the next call.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of values desired.  If N is negative,
    then the code will flush its internal memory; in particular,
    if there is a saved value to be used on the next call, it is
    instead discarded.  This is useful if the user has reset the
    random number seed, for instance.
    Input/output, int *SEED, a seed for the random number generator.
    Output, ityp X[N], a sample of the standard normal PDF.
  Local parameters:
    Local, int MADE, records the number of values that have
    been computed.  On input with negative N, this value overwrites
    the return value of N, so the user can get an accounting of
    how much work has been done.
    Local, ityp R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.
    Local, int SAVED, is 0 or 1 depending on whether there is a
    single saved value left over from the previous call.
    Local, int X_LO, X_HI, records the range of entries of
    X that we need to compute.  This starts off as 1:N, but is adjusted
    if we have a saved value that can be immediately stored in X(1),
    and so on.
    Local, ityp Y, the value saved from the previous call, if
    SAVED is 1.
*/
{
	const dtpitpi * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int * seed = s_data->a2;
	
    dim_typ i, m;
    static int made = 0;
    ityp *r;
    static int saved = false;
    int x_hi;
    int x_lo;
    static ityp y = 0.00;
    /*
    I'd like to allow the user to reset the internal data.
    But this won't work properly if we have a saved value Y.
    I'm making a crock option that allows the user to signal
    explicitly that any internal memory should be flushed,
    by passing in a negative value for N.
    */
    if ( n < 0 )
    {
        made = saved = 0;
        y = 0.00;
        return NULL;
    }
    else if ( n == 0 )
        return NULL;
    /*
    Record the range of X we need to fill in.
    */
    x_lo = 1;
    x_hi = n;
    /*
    Use up the old value, if we have it.
    */
    if ( saved == 1 )
    {
        x[0] = y;
        saved = 0;
        x_lo = 2;
    }
    /*
    Maybe we don't need any more values.
    */
    if ( x_hi - x_lo + 1 == 0 );
    /*
    If we need just one new value, do that here to avoid null arrays.
    */
    else if ( x_hi - x_lo + 1 == 1 )
    {
        r = r8vec_uniform_01_new ( 2, seed );

        x[x_hi-1] = sqrt ( - 2.00 * log ( r[0] ) ) * cos ( M_2TPI * r[1] );
        y =         sqrt ( - 2.00 * log ( r[0] ) ) * sin ( M_2TPI * r[1] );

        saved = true;
        made += 2;
        free ( r );
    }
    /*
    If we require an even number of values, that's easy.
    */
    else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
    {
        m = ( x_hi - x_lo + 1 ) / 2;

        r = r8vec_uniform_01_new ( m<<1, seed );

        for ( i = 0; i <= (m<<1)-2; i += 2 )
        {
            x[x_lo+i-1] = sqrt ( - 2.00 * log ( r[i] ) ) * cos ( M_2TPI * r[i+1] );
            x[x_lo+i  ] = sqrt ( - 2.00 * log ( r[i] ) ) * sin ( M_2TPI * r[i+1] );
        }
        made += x_hi - x_lo + 1;
        free ( r );
    }
    /*
    If we require an odd number of values, we generate an even number,
    and handle the last pair specially, storing one in X(N), and
    saving the other for later.
    */
    else
    {
        -- x_hi;
        m = ( x_hi - x_lo + 1 ) / 2 + 1;
        r = r8vec_uniform_01_new ( m<<1, seed );

        for ( i = 0; i <= (m<<1)-4; i += 2 )
        {
            x[x_lo+i-1] = sqrt ( - 2.00 * log ( r[i] ) ) * cos ( M_2TPI * r[i+1] );
            x[x_lo+i  ] = sqrt ( - 2.00 * log ( r[i] ) ) * sin ( M_2TPI * r[i+1] );
        }

        i = (m<<1) - 2;

        x[x_lo+i-1] = sqrt ( - 2.00 * log ( r[i] ) ) * cos ( M_2TPI * r[i+1] );
        y           = sqrt ( - 2.00 * log ( r[i] ) ) * sin ( M_2TPI * r[i+1] );

        saved = true;
        made += x_hi - x_lo + 2;
        free ( r );
    }

    return NULL;
}

#endif
