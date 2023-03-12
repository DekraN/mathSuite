#ifndef __DISABLEDEEP_ASA266

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _alogam ( void * data)
/******************************************************************************/
/*
  Purpose:
    ALOGAM computes the logarithm of the Gamma function.
  Modified:
    22 January 2008
  Author:
    Malcolm Pike, David Hill
    FORTRAN90 version by John Burkardt
  Reference:
    Malcolm Pike, David Hill,
    Algorithm 291:
    Logarithm of Gamma Function,
    Communications of the ACM,
    Volume 9, Number 9, September 1966, page 684.
  Parameters:
    Input, double X, the argument of the Gamma function.
    X should be greater than 0.
    Output, int *IFAULT, error flag.
    0, no error.
    1, X <= 0.
    Output, double ALOGAM, the logarithm of the Gamma
    function of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp f;
    ityp value;
    ityp y;
    ityp z;

    if ( x <= 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    y = x;

    if ( x < 7.00 )
    {
        f = 1.00;
        z = y;

        while ( z < 7.00 )
        {
            f *= z;
            ++ z;
        }
        y = z;
        f = - log ( f );
    }
    else
        f = 0.0;

    z = 1.00 / y / y;

	result = f + ( y - 0.50 ) * log ( y ) - y + 0.918938533204673 + ((( - 0.000595238095238   * z + 0.000793650793651 ) * z - 0.002777777777778 ) * z + 0.083333333333333 ) / y;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gamain ( void * data)
/******************************************************************************/
/*
  Purpose:
    GAMAIN computes the incomplete gamma ratio.
  Discussion:
    A series expansion is used if P > X or X <= 1.  Otherwise, a
    continued fraction approximation is used.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 June 2014
  Author:
    Original FORTRAN77 version by G Bhattacharjee.
    C version by John Burkardt.
  Reference:
    G Bhattacharjee,
    Algorithm AS 32:
    The Incomplete Gamma Integral,
    Applied Statistics,
    Volume 19, Number 3, 1970, pages 285-287.
  Parameters:
    Input, double X, P, the parameters of the incomplete
    gamma ratio.  0 <= X, and 0 < P.
    Output, int *IFAULT, error flag.
    0, no errors.
    1, P <= 0.
    2, X < 0.
    3, underflow.
    4, error return from the Log Gamma routine.
    Output, double GAMAIN, the value of the incomplete gamma ratio.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp x = a_data[0];
	const register ityp p = a_data[1];
	
	ityp a;
	ityp an;
	ityp arg;
	ityp b;
	ityp dif;
	ityp factor;
	ityp g;
	ityp gin;
	dim_typ i;
	ityp pn[6];
	ityp rn;
	ityp term;

	/*
	Check the input.
	*/
	if (p <= 0.00 || x <= 0.00)
	{
		result = GAMAIN_INVALIDRETURNVALUE;
		return &result;
	}


	g = lgamma(p);
	arg = p * log (x) - x - g;

	if (arg < log(1.0E-37))
	{
		result = GAMAIN_INVALIDRETURNVALUE;
		return &result;
	}

	factor = exp (arg);
	/*
	Calculation by series expansion.
	*/
	if (x <= 1.00 || x < p)
	{
		gin = 1.00;
		term = 1.00;
		rn = p;

		for ( ; ; )
		{
			rn = rn + 1.00;
			term = term * x / rn;
			gin = gin + term;

			if ( term <= 1.0E-08 )
				break;
		}

		result = gin * factor / p;
		return &result;
	}

	ityp value;

	/*
	Calculation by continued fraction.
	*/
	a = 1.00 - p;
	b = a + x + 1.00;
	term = 0.00;

	pn[0] = 1.00;
	pn[1] = x;
	pn[2] = x + 1.00;
	pn[3] = x * b;

	gin = pn[2] / pn[3];

	for ( ; ; )
	{
		a += 1.00;
		b += 2.00;
		term += 1.00;
		an = a * term;

		for ( i = 0; i <= 1; ++i )
			pn[i+4] = b * pn[i+2] - an * pn[i];

		if (pn[5])
		{
			rn = pn[4] / pn[5];
			dif = fabs ( gin - rn );
			/*
			Absolute error tolerance satisfied?
			*/
			if ( dif <= 1.0E-08 )
			{
				/*
				Relative error tolerance satisfied?
				*/
				if ( dif <= 1.0E-08 * rn )
					value = 1.0 - factor * gin;
			}
			gin = rn;
		}

		for (i = 0; i < 4; ++i)
			pn[i] = pn[i+2];

		if(1.0E+37 <= fabs (pn[4]))
			#pragma omp parallel for num_threads(4)
			for ( i = 0; i < 4; ++i )
				pn[i] /= 1.0E+37;
	}

	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8col_mean ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8COL_MEAN returns the column means of an r8COL.
  Discussion:
    An r8COL is an M by N array of r8's, regarded
    as an array of N columns of length M.
  Example:
    A =
      1  2  3
      2  6  7
    r8COL_MEAN =
      1.5  4.0  5.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 October 2008
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns.
    Input, ityp A[M*N], the array to be examined.
    Output, ityp r8COL_MEAN[N], the means, or averages, of the columns.
*/
{
	const _2dtpit * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	
    dim_typ i, j;
    ityp *mean;

    mean = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
    {
        mean[j] = 0.00;
        for ( i = 0; i < m; ++i )
            mean[j] += a[i+j*m];
        mean[j] /= ( ityp ) ( m );
    }

    return mean;
}

#endif
