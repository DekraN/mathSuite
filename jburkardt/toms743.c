#ifndef __DISABLEDEEP_TOMS743

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bisect ( void * data)
/******************************************************************************/
/*
  Purpose:
    BISECT approximates the W function using bisection.
  Discussion:
    The parameter TOL, which determines the accuracy of the bisection
    method, is calculated using NBITS (assuming the final bit is lost
    due to rounding error).
    N0 is the maximum number of iterations used in the bisection
    method.
    For XX close to 0 for Wp, the exponential approximation is used.
    The approximation is exact to O(XX^8) so, depending on the value
    of NBITS, the range of application of this formula varies. Outside
    this range, the usual bisection method is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 June 2014
  Author:
    Original FORTRAN77 version by Andrew Barry, S. J. Barry,
    Patricia Culligan-Hensley.
    This C version by John Burkardt.
  Reference:
    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
    Algorithm 743: WAPR - A Fortran routine for calculating real
    values of the W-function,
    ACM Transactions on Mathematical Software,
    Volume 21, Number 2, June 1995, pages 172-181.
  Parameters:
    Input, double XX, the argument.
    Input, int NB, indicates the branch of the W function.
    0, the upper branch;
    nonzero, the lower branch.
    Output, int *NER, the error flag.
    0, success;
    1, the routine did not converge.  Perhaps reduce NBITS and try again.
    Input, int L, the offset indicator.
    1, XX represents the offset of the argument from -exp(-1).
    not 1, XX is the actual argument.
    Output, double BISECT, the value of W(X), as determined
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtitpi * const s_data = data;
	
	const register dim_typ nb = s_data->a0;
	const register dim_typ l = s_data->a1;
	const register ityp xx = s_data->a2;
	int * ner = s_data->a3; 
	
	
    ityp d;
    ityp f;
    ityp fd;
    dim_typ i;
    const register dim_typ n0 = 500;
    static dim_typ nbits = 0;
    ityp r;
    ityp test;
    ityp tol;
    ityp u;
    ityp value;
    ityp x;

    value = 0.00;
    *ner = 0;

    if ( nbits == 0 )
        nbits = nbits_compute ( );

    x = xx + (l==1)*(- exp ( -1.00 ));

    if ( nb == 0 )
    {
        test = 1.0 / pow ( pow ( 2.0, nbits ), ( 1.0 / 7.0 ) );

        if ( fabs ( x ) < test )
        {
        	result = x * exp ( - x * exp ( - x * exp ( - x * exp ( - x * exp ( - x * exp ( - x ))))));
            return &result;
        }
        else
        {
            u = crude ( x, nb ) + 1.0E-03;
            tol = fabs ( u ) / pow ( 2.0, nbits );
            d = MAX ( u - 2.0E-03, -1.0 );

            for ( i = 1; i <= n0; ++i )
            {
                r = 0.50 * ( u - d );
                value = d + r;
                /*
                Find root using w*exp(w)-x to avoid ln(0) error.
                */
                if ( x < exp ( 1.00 ) )
                {
                    f = value * exp ( value ) - x;
                    fd = d * exp ( d ) - x;
                }
                /*
                Find root using ln(w/x)+w to avoid overflow error.
                */
                else
                {
                    f = log ( value / x ) + value;
                    fd = log ( d / x ) + d;
                }

                if ( f == 0.00 || fabs ( r ) <= tol )
                {
                	result = value;
                    return &result;
                }

                if ( 0.0 < fd * f )
                    d = value;
                else
                    u = value;
            }
        }
    }
    else
    {
        d = crude ( x, nb ) - 1.0E-03;
        u = MIN ( d + 2.0E-03, -1.0 );
        tol = fabs ( u ) / pow ( 2.00, nbits );

        for ( i = 1; i <= n0; ++i )
        {
            r = 0.50 * ( u - d );
            value = d + r;
            f = value * exp ( value ) - x;

            if ( f == 0.00 || fabs ( r ) <= tol )
            {
            	result = value;
                return &result;
            }

            fd = d * exp ( d ) - x;

            if ( 0.00 < fd * f )
                d = value;
            else
                u = value;
        }
    }
    /*
    The iteration did not converge.
    */
    *ner = 1;

    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _crude ( void * data)
/******************************************************************************/
/*
  Purpose:
    CRUDE returns a crude approximation for the W function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 June 2014
  Author:
    Original FORTRAN77 version by Andrew Barry, S. J. Barry,
    Patricia Culligan-Hensley.
    This C version by John Burkardt.
  Reference:
    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
    Algorithm 743: WAPR - A Fortran routine for calculating real
    values of the W-function,
    ACM Transactions on Mathematical Software,
    Volume 21, Number 2, June 1995, pages 172-181.
  Parameters:
    Input, double XX, the argument.
    Input, int NB, indicates the desired branch.
    * 0, the upper branch;
    * nonzero, the lower branch.
    Output, double CRUDE, the crude approximation to W at XX.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ nb = s_data->a0;
	const register ityp xx = s_data->a1;
	
    ityp an2;
    static ityp c13;
    ityp crude;
    static ityp em;
    static ityp em2;
    static ityp em9;
    ityp eta;
    static dim_typ init = 0;
    ityp reta;
    static ityp s2;
    static ityp s21;
    static ityp s22;
    static ityp s23;
    ityp t;
    ityp ts;
    ityp value;
    ityp zl;

    value = 0.00;
    /*
    Various mathematical constants.
    */
    if ( init == 0 )
    {
        init = 1;
        em = - exp ( -1.00 );
        em9 = - exp ( -9.00 );
        c13 = 1.00 / 3.00;
        em2 = 2.00 / em;
        s2 = sqrt ( 2.00 );
        s21 = 2.00 * s2 - 3.00;
        s22 = 4.00 - 3.00 * s2;
        s23 = s2 - 2.00;
    }
    /*
    Crude Wp.
    */
    if ( nb == 0 )
    {
        if ( xx <= 20.0 )
        {
            reta = s2 * sqrt ( 1.00 - xx / em );
            an2 = 4.612634277343749 * sqrt ( sqrt ( reta +1.09556884765625 ) );
            value = reta / ( 1.0 + reta / ( 3.00+ (  s21 * an2 + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.00;
        }
        else
        {
            zl = log ( xx );
            value = log ( xx / log ( xx/ pow ( zl, exp ( -1.124491989777808 /( 0.4225028202459761 + zl )) ) ));
        }
    }
    /*
    Crude Wm.
    */
    else
    {
        if ( xx <= em9 )
        {
            zl = log ( -xx );
            t = -1.00 - zl;
            ts = sqrt ( t );
            value = zl - ( 2.00 * ts ) / ( s2 + ( c13 - t/ ( 270.00 + ts * 127.0471381349219 ) ) * ts );
        }
        else
        {
            zl = log ( -xx );
            eta = 2.00 - em2 * xx;
            value = log ( xx / log ( - xx / ( ( 1.0- 0.5043921323068457 * ( zl + 1.00 ) )* ( sqrt ( eta ) + eta / 3.00 ) + 1.00 ) ) );
        }
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _nbits_compute ( void * data)
/******************************************************************************/
/*
  Purpose:
    NBITS_COMPUTE computes the mantissa length minus one.
  Discussion:
    NBITS is the number of bits (less 1) in the mantissa of the
    floating point number number representation of your machine.
    It is used to determine the level of accuracy to which the W
    function should be calculated.
    Most machines use a 24-bit matissa for single precision and
    53-56 bits for double. The IEEE standard is 53
    bits. The Fujitsu VP2200 uses 56 bits. Long word length
    machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
    single precision.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 June 2014
  Author:
    Original FORTRAN77 version by Andrew Barry, S. J. Barry,
    Patricia Culligan-Hensley.
    This C version by John Burkardt.
  Reference:
    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
    Algorithm 743: WAPR - A Fortran routine for calculating real
    values of the W-function,
    ACM Transactions on Mathematical Software,
    Volume 21, Number 2, June 1995, pages 172-181.
  Parameters:
    Output, int NBITS_COMPUTE, the mantissa length, in bits,
    minus one.
*/
{
	static dim_typ result = USHRT_MAX;
	
    ityp b = 1.00;
    dim_typ i;
    dim_typ nbits = 0;
    ityp v;

    for ( ; ; )
    {
        b /= 2.00;
        v = b + 1.00;

        if ( v == 1.0 )
            break;
            ++ nbits;
    }

	result = nbits;
    return &result;
}

#endif
