#ifndef __DISABLEDEEP_PROB

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _zipf_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    ZIPF_CDF evaluates the Zipf CDF.
  Discussion:
    Simple summation is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, int X, the argument of the PDF.
    1 <= N
    Input, double A, the parameter of the PDF.
    1.0 < A.
    Output, double CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ x = s_data->a0;
	const register ityp a = s_data->a1;
	
    ityp c;
    ityp cdf;
    ityp pdf;
    dim_typ y;

    if ( x < 1 )
        cdf = 0.0;
    else
    {
        c = zeta ( a );

        cdf = 0.00;
        for ( y = 1; y <= x; ++y)
        {
            pdf = 1.00 / pow ( y, a ) / c;
            cdf += pdf;
        }
    }
	
	result = cdf;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _zipf_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    ZIPF_PDF evaluates the Zipf PDF.
  Discussion:
    PDF(A;X) = ( 1 / X^A ) / C
    where the normalizing constant is chosen so that
    C = Sum ( 1 <= I < Infinity ) 1 / I^A.
    From observation, the frequency of different words in long
    sequences of text seems to follow the Zipf PDF, with
    parameter A slightly greater than 1.  The Zipf PDF is sometimes
    known as the "discrete Pareto" PDF.
    Lotka's law is a version of the Zipf PDF in which A is 2 or approximately
    2.  Lotka's law describes the frequency of publications by authors in a
    given field, and estimates that the number of authors with X papers is
    about 1/X^A of the number of authors with 1 paper.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 October 2004
  Author:
    John Burkardt
  Reference:
    Alfred Lotka,
    The frequency distribution of scientific productivity,
    Journal of the Washington Academy of Sciences,
    Volume 16, Number 12, 1926, pages 317-324.
  Parameters:
    Input, int X, the argument of the PDF.
    1 <= N
    Input, double A, the parameter of the PDF.
    1.0 < A.
    Output, double PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ x = s_data->a0;
	const register ityp a = s_data->a1;
	
	result = 0.00 + (x!=0)*(1.0 / pow ( x, a ) / zeta ( a ));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _buffon_laplace_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    BUFFON_LAPLACE_PDF evaluates the Buffon-Laplace PDF.
  Discussion:
    In the Buffon-Laplace needle experiment, we suppose that the plane has
    been tiled into a grid of rectangles of width A and height B, and that a
    needle of length L is dropped "at random" onto this grid.
    We may assume that one end, the "eye" of the needle falls at the point
 (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
    ANGLE, the angle that the needle makes is taken to be uniformly random.
    The point of the needle, (X2,Y2), therefore lies at
   (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
    The needle will have crossed at least one grid line if any of the
    following are true:
      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
    If L is larger than sqrt ( A*A + B*B ), then the needle will
    cross every time, and the computation is uninteresting.  However, if
    L is smaller than this limit, then the probability of a crossing on
    a single trial is
      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( M_PI * A * B )
    and therefore, a record of the number of hits for a given number of
    trials can be used as a very roundabout way of estimating M_PI.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 February 2007
  Author:
    John Burkardt
  Reference:
    Sudarshan Raghunathan,
    Making a Supercomputer Do What You Want: High Level Tools for
    Parallel Programming,
    Computing in Science and Engineering,
    Volume 8, Number 5, September/October 2006, pages 70-80.
  Parameters:
    Input, double A, B, the horizontal and vertical dimensions
    of each cell of the grid.  0 <= A, 0 <= B.
    Input, double L, the length of the needle.
    0 <= L <= MIN ( A, B ).
    Output, double BUFFON_LAPLACE_PDF, the Buffon-Laplace PDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp b = a_data[1];
	const register ityp l = a_data[2];
	
    if ( a < 0.0 )
    {
    	result = MAX_VAL;
        return &result;
    }
    else if ( a == 0.0 )
    {
    	result = 1.00;
        return &result;
    }
	
	ityp value = ( 2.00 * l * ( a + b ) - l * l ) / ( M_PI * a * b );
    if ( a == 0.00 || b == 0.00 )
        value = 1.00;
    else if ( l == 0.00 )
        value = 0.00;
    else if (a == 0.00 || l<0.00 || MIN ( a, b ) < l )
        value = MAX_VAL;
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laplace_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAPLACE_CDF evaluates the Laplace CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 February 1999
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the CDF.
    Input, double A, B, the parameters of the PDF.
    0.0 < B.
    Output, double LAPLACE_CDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	
	const register ityp y = (x-a)/b;
	
	result = x <= a ? 0.50 * exp ( y ) : 1.00 - 0.50 * exp ( - y );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laplace_cdf_inv ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAPLACE_CDF_INV inverts the Laplace CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 February 1999
  Author:
    John Burkardt
  Parameters:
    Input, double CDF, the value of the CDF.
    0.0 <= CDF <= 1.0.
    Input, double A, B, the parameters of the PDF.
    0.0 < B.
    Output, double LAPLACE_CDF_INV, the corresponding argument.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp cdf = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	
	result = cdf < 0.00 || 1.00 < cdf ? MAX_VAL : cdf <= 0.50 ? a + b * log ( 2.00 * cdf ) : a - b * log ( 2.00 * ( 1.00 - cdf ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _laplace_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    LAPLACE_PDF evaluates the Laplace PDF.
  Discussion:
    PDF(A,B;X) = exp ( - abs ( X - A ) / B ) / ( 2 * B )
  Discussion:
    The Laplace PDF is also known as the Double Exponential PDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 February 1999
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the PDF.
    Input, double A, B, the parameters of the PDF.
    0.0 < B.
    Output, double LAPLACE_PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	
	result = exp ( - abs ( x - a ) / b ) / ( 2.00 * b );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _log_normal_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOG_NORMAL_CDF evaluates the Lognormal CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 September 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the PDF.
    0.0 < X.
    Input, double A, B, the parameters of the PDF.
    0.0 < B.
    Output, double CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	
	result = 0.00 + (x>0.00)*normal_cdf ( log(x), a, b );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _normal_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    NORMAL_CDF evaluates the Normal CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 September 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the CDF.
    Input, double A, B, the parameters of the PDF.
    0.0 < B.
    Output, double CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	
	result = normal_01_cdf ( ( x - a ) / b );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _log_normal_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOG_NORMAL_PDF evaluates the Lognormal PDF.
  Discussion:
    PDF(A,B;X)
      = exp ( - 0.5 * ( ( log ( X ) - A ) / B )^2 )
        / ( B * X * sqrt ( 2 * M_PI ) )
    The Lognormal PDF is also known as the Cobb-Douglas PDF,
    and as the Antilog_normal PDF.
    The Lognormal PDF describes a variable X whose logarithm
    is normally distributed.
    The special case A = 0, B = 1 is known as Gilbrat's PDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 September 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the PDF.
    0.0 < X
    Input, double A, B, the parameters of the PDF.
    0.0 < B.
    Output, double PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	
	const register ityp y = (log(x)-a)/b;
	
	result = 0.00 + (x>0.00)*exp ( -0.50 * (y * y ) / ( b * x * sqrt ( M_2TPI ) )); 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _log_series_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOG_SERIES_PDF evaluates the Logarithmic Series PDF.
  Discussion:
    PDF(A;X) = - A**X / ( X * log ( 1 - A ) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, int X, the argument of the PDF.
    0 < X
    Input, double A, the parameter of the PDF.
    0.0 < A < 1.0.
    Output, double LOG_SERIES_PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ x = s_data->a0;
	const register ityp a = s_data->a1;
	
	result = 0.00 + (x>0.00)*(- pow ( a, x ) / ( ( ityp ) ( x ) * log ( 1.00 - a ) ));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _genlogistic_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    GENLOGISTIC_CDF evaluates the Generalized Logistic CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the CDF.
    Input, double A, B, C, the parameters of the PDF.
    0.0 < B,
    0.0 < C.
    Output, double CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	const register ityp c = a_data[3];
	
	result = 1.00 / pow ( ( 1.0 + exp ( - ( x - a ) / b ) ), c );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _logistic_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    LOGISTIC_CDF evaluates the Logistic CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the CDF.
    Input, double A, B, the parameters of the PDF.
    0.0 < B.
    Output, double LOGISTIC_CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	
	result = 1.00 / ( 1.00 + exp ( ( a - x ) / b ) ); 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _binomial_coef ( void * data)
/******************************************************************************/
/*
  Purpose:
    BINOMIAL_COEF computes the Binomial coefficient C(N,K).
  Discussion:
    The value is calculated in such a way as to avoid overflow and
    roundoff.  The calculation is done in integer arithmetic.
    CNK = C(N,K) = N! / ( K! * (N-K)! )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 October 2004
  Author:
    John Burkardt
  Reference:
    M L Wolfson, H V Wright,
    Combinatorial of M Things Taken N at a Time,
    ACM Algorithm 160,
    Communications of the ACM,
    April, 1963.
  Parameters:
    Input, int N, K, are the values of N and K.
    Output, int CNK, the number of combinations of N
    things taken K at a time.
*/
{
	static dim_typ result = USHRT_MAX;
	
	dim_typ * const a_data = data;
	const register dim_typ n = a_data[0];
	const register dim_typ k = a_data[1];
	
	int cnk;
	int i;
	int mn;
	int mx;
	
	mn = MIN ( k, n-k );
	
	if ( mn == 0 )
		cnk = 1;
	else
	{
		mx = MAX ( k, n-k );
		cnk = mx + 1;
		
		for ( i = 2; i <= mn; ++i )
			cnk = ( cnk * ( mx + i ) ) / i;
	}
	
	result = cnk; 
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _negative_binomial_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    NEGATIVE_BINOMIAL_CDF evaluates the Negative Binomial CDF.
  Discussion:
    A simple summing approach is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, int X, the argument of the CDF.
    Input, int A, double B, parameters of the PDF.
    0 <= A,
    0 < B <= 1.
    Output, double NEGATIVE_BINOMIAL_CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	const register dim_typ x = s_data->a0;
	const register dim_typ a = s_data->a1;
	const register ityp b = s_data->a2;
	
    ityp cdf = 0.00;
    for (dim_typ y = a; y <= x; ++y )
        cdf += ( ityp ) ( binomial_coef ( y-1, a-1 ) ) * pow ( b, a ) * pow ( 1.00 - b, y - a );
    
	result = cdf;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _negative_binomial_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    NEGATIVE_BINOMIAL_PDF evaluates the Negative Binomial PDF.
  Discussion:
    PDF(A,B;X) = C(X-1,A-1) * B^A * ( 1 - B )^(X-A)
    PDF(A,B;X) is the probability that the A-th success will
    occur on the X-th trial, given that the probability
    of a success on a single trial is B.
    The Negative Binomial PDF is also known as the Pascal PDF or
    the "Polya" PDF.
    NEGATIVE_BINOMIAL_PDF(1,B;X) = GEOMETRIC_PDF(B;X)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, int X, the number of trials.
    A <= X.
    Input, int A, the number of successes required.
    0 <= A <= X, normally.
    Input, double B, the probability of a success on a single trial.
    0.0 < B <= 1.0.
    Output, double NEGATIVE_BINOMIAL_PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	const _2dtit * const s_data = data;
	const register dim_typ x = s_data->a0;
	const register dim_typ a = s_data->a1;
	const register ityp b = s_data->a2;
	
    result = 0.00 + (x>=a)*(( ityp ) ( binomial_coef ( x-1, a-1 ) ) * pow ( b, a ) * pow ( 1.00 - b, x - a ));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _poisson_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    POISSON_PDF evaluates the Poisson PDF.
  Discussion:
    PDF(K,A) is the probability that the number of events observed
    in a unit time period will be K, given the expected number
    of events in a unit time.
    The parameter A is the expected number of events per unit time.
    The Poisson PDF is a discrete version of the exponential PDF.
    The time interval between two Poisson events is a random
    variable with the exponential PDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    01 February 1999
  Author:
    John Burkardt
  Parameters:
    Input, int K, the argument of the PDF.
    Input, double A, the parameter of the PDF.
    0 < A.
    Output, double POISSON_PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	const dtit * const s_data = data;
	const register dim_typ k = s_data->a0;
	const register ityp a = s_data->a1;
	
	result = a<=0.00 ? MAX_VAL : exp ( -a ) * pow ( a, ( ityp  ) k ) / i4_factorial ( k ) ;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rayleigh_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    RAYLEIGH_CDF evaluates the Rayleigh CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the CDF.
    0.0 <= X.
    Input, double A, the parameter of the PDF.
    0.0 < A.
    Output, double RAYLEIGH_CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	
	result = 0.00 + (x>=0.00)*(1.00 - exp ( - x * x / ( 2.00 * a * a ) ));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rayleigh_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    RAYLEIGH_PDF evaluates the Rayleigh PDF.
  Discussion:
    PDF(A;X) = ( X / A^2 ) * EXP ( - X^2 / ( 2 * A^2 ) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 October 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the PDF.
    0.0 <= X
    Input, double A, the parameter of the PDF.
    0 < A.
    Output, double RAYLEIGH_PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	
	result = 0.00 + (x>=0.00)*(( x / ( a * a ) ) * exp ( - x * x / ( 2.00 * a * a ) ));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _von_mises_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    VON_MISES_PDF evaluates the von Mises PDF.
  Discussion:
    PDF(A,B;X) = EXP ( B * COS ( X - A ) ) / ( 2 * M_PI * I0(B) )
    where:
      I0(*) is the modified Bessel function of the first
      kind of order 0.
    The von Mises distribution for points on the unit circle is
    analogous to the normal distribution of points on a line.
    The variable X is interpreted as a deviation from the angle A,
    with B controlling the amount of dispersion.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 October 2004
  Author:
    John Burkardt
  Reference:
    Jerry Banks, editor,
    Handbook of Simulation,
    Engineering and Management Press Books, 1998, page 160.
    D J Best, N I Fisher,
    Efficient Simulation of the von Mises Distribution,
    Applied Statistics,
    Volume 28, Number 2, pages 152-157.
    Kanti Mardia, Peter Jupp,
    Directional Statistics,
    Wiley, 2000, QA276.M335
  Parameters:
    Input, double X, the argument of the PDF.
    A - M_PI <= X <= A + M_PI.
    Input, double A, B, the parameters of the PDF.
    -M_PI <= A <= M_PI,
    0.0 < B.
    Output, double VON_MISES_PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	
    if ( x < a - M_PI )
    {
    	result = 0.00;
        return &result;
    }
    else if ( x <= a + M_PI )
    {
    	result = exp ( b * cos ( x - a ) ) / ( M_2TPI * bessel_i0 ( b ) );
        return &result;
    }
    
    result = 0.00;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bessel_i0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    BESSEL_I0 evaluates the modified Bessel function I0.
  Discussion:
    The main computation evaluates slightly modified forms of
    minimax approximations generated by Blair and Edwards, Chalk
    River (Atomic Energy of Canada Limited) Report AECL-4928,
    October, 1974.  This transportable program is patterned after
    the machine dependent FUNPACK packet NATSI0, but cannot match
    that version for efficiency or accuracy.  This version uses
    rational functions that theoretically approximate I-SUB-0(X)
    to at least 18 significant decimal digits.
  Machine dependent constants:
    beta   = Radix for the floating-point system
    maxexp = Smallest power of beta that overflows
    XSMALL = Positive argument such that 1.0 - X = 1.0 to
             machine precision for all ABS(X) .LE. XSMALL.
    XMAX =   Largest argument acceptable to BESI0;  Solution to
             equation:
               W(X) * (1+1/(8*X)+9/(128*X^2) = beta**maxexp
             where  W(X) = EXP(X)/sqrt(2*PI*X)
    Approximate values for some important machines are:
                             beta       maxexp       XSMALL
    CRAY-1     (S.P.)       2         8191       3.55D-15
    Cyber 180/855
      under NOS  (S.P.)       2         1070       3.55D-15
    IEEE (IBM/XT,
      SUN, etc.) (S.P.)       2          128       2.98D-8
    IEEE (IBM/XT,
      SUN, etc.) (D.P.)       2         1024       5.55D-17
    IBM 3033   (D.P.)      16           63       6.95D-18
    VAX        (S.P.)       2          127       2.98D-8
    VAX D-Format (D.P.)       2          127       6.95D-18
    VAX G-Format (D.P.)       2         1023       5.55D-17
                                  XMAX
    CRAY-1     (S.P.)       5682.810
    Cyber 180/855
      under NOS  (S.P.)       745.893
    IEEE (IBM/XT,
      SUN, etc.) (S.P.)        91.900
    IEEE (IBM/XT,
      SUN, etc.) (D.P.)       713.986
    IBM 3033   (D.P.)       178.182
    VAX        (S.P.)        91.203
    VAX D-Format (D.P.)        91.203
    VAX G-Format (D.P.)       713.293
  Author:
    Original FORTRAN77 version by W. J. Cody and L. Stoltz.
    C version by John Burkardt.
  Parameters:
    Input, double ARG, the argument.
    Output, double BESSEL_I0, the value of the modified Bessel function
    of the first kind.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp arg = *(ityp *) data;
	
	ityp a;
	ityp b;
	ityp exp40 = 2.353852668370199854E+17;
	dim_typ i;
	ityp p[15] = 
	{
		-5.2487866627945699800E-18,
		-1.5982226675653184646E-14,
		-2.6843448573468483278E-11,
		-3.0517226450451067446E-08,
		-2.5172644670688975051E-05,
		-1.5453977791786851041E-02,
		-7.0935347449210549190,
		-2.4125195876041896775E+03,
		-5.9545626019847898221E+05,
		-1.0313066708737980747E+08,
		-1.1912746104985237192E+10,
		-8.4925101247114157499E+11,
		-3.2940087627407749166E+13,
		-5.5050369673018427753E+14,
		-2.2335582639474375249E+15 
	};
	ityp pp[8] = 
	{
		-3.9843750000000000000E-01,
		2.9205384596336793945,
		-2.4708469169133954315,
		4.7914889422856814203E-01,
		-3.7384991926068969150E-03,
		-2.6801520353328635310E-03,
		9.9168777670983678974E-05,
		-2.1877128189032726730E-06 
	};
	ityp q[5] = 
	{
		-3.7277560179962773046E+03,
		6.5158506418655165707E+06,
		-6.5626560740833869295E+09,
		3.7604188704092954661E+12,
		-9.7087946179594019126E+14 
	};
	ityp qq[7] = 
	{
		-3.1446690275135491500E+01,
		8.5539563258012929600E+01,
		-6.0228002066743340583E+01,
		1.3982595353892851542E+01,
		-1.1151759188741312645,
		3.2547697594819615062E-02,
		-5.5194330231005480228E-04 
	};
	const ityp rec15 = 6.6666666666666666666E-02;
	ityp sump;
	ityp sumq;
	ityp value;
	ityp x;
	const ityp xmax = 91.9;
	const ityp xsmall = 2.98E-08;
	ityp xx;
	
	x = abs ( arg );
	
	if ( x < xsmall )
		value = 1.00;
	else if ( x < 15.00 )
	{
		/*
		XSMALL <= ABS(ARG) < 15.0
		*/
		xx = x * x;
		sump = p[0];
		#pragma omp parallel for num_threads(14)
		for ( i = 1; i < 15; ++i )
			sump = sump * xx + p[i];
	
		xx -= 225.00;
		sumq = ((((xx + q[0] )* xx + q[1] )* xx + q[2] )* xx + q[3] )* xx + q[4];
		
		value = sump / sumq;
	}
	else if ( 15.00 <= x  )
	{
		if(xmax < x)
			value = r8_huge;
		else
		{
			/*
			15.0 <= ABS(ARG)
			*/
			xx = 1.00 / x - rec15;
			
			sump = ((((((pp[0]* xx + pp[1] )* xx + pp[2] )* xx + pp[3] )* xx + pp[4] )* xx + pp[5] )* xx + pp[6] )* xx + pp[7];
			sumq = ((((((xx + qq[0] )* xx + qq[1] )* xx + qq[2] )* xx + qq[3] )* xx + qq[4] )* xx + qq[5] )* xx + qq[6];
			
			value = sump / sumq;
			/*
			Calculation reformulated to avoid premature overflow.
			*/
			if ( x <= xmax - 15.00 )
			{
				a = exp ( x );
				b = 1.00;
			}
			else
			{
				a = exp ( x - 40.00 );
				b = exp40;
			}
		
			value = ( ( value * a - pp[0] * a ) / sqrt ( x ) ) * b;
		}
	}
	
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _weibull_cdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEIBULL_CDF evaluates the Weibull CDF.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    18 September 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the CDF.
    A <= X.
    Input, double A, B, C, the parameters of the PDF.
    0.0 < B,
    0.0 < C.
    Output, double CDF, the value of the CDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	const register ityp c = a_data[3];
	
	result = 0.00 + (x>=a)*(1.00 - 1.00 / exp ( pow ( ( x - a ) / b, c ) ));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _weibull_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEIBULL_PDF evaluates the Weibull PDF.
  Discussion:
    PDF(A,B,C;X) = ( C / B ) * ( ( X - A ) / B )^( C - 1 )
     * EXP ( - ( ( X - A ) / B )^C ).
    The Weibull PDF is also known as the Frechet PDF.
    WEIBULL_PDF(A,B,1;X) is the Exponential PDF.
    WEIBULL_PDF(0,1,2;X) is the Rayleigh PDF.
  Modified:
    18 September 2004
  Author:
    John Burkardt
  Parameters:
    Input, double X, the argument of the PDF.
    A <= X
    Input, double A, B, C, the parameters of the PDF.
    0.0 < B,
    0.0 < C.
    Output, double PDF, the value of the PDF.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp a = a_data[1];
	const register ityp b = a_data[2];
	const register ityp c = a_data[3];
	
    ityp y;
    result = 0.00 + (x>=a)*(( c / b ) * pow ( y = ( x - a ) / b, ( c - 1.0 ) )  / exp ( pow ( y, c ) )); 
    return &result;
}

#endif
