#ifndef __DISABLEDEEP_ASA152

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _alnfac ( void * data)
/******************************************************************************/
/*
  Purpose:
    ALNFAC computes the logarithm of the factorial of N.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the argument of the factorial.
    Output, double ALNFAC, the logarithm of the factorial of N.
*/
{
	static ityp result = MAX_VAL;
	
	result = lgamma((ityp)((*(dim_typ*)data)+1));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _chyper ( void * data)
/******************************************************************************/
/*
  Purpose:
    CHYPER computes point or cumulative hypergeometric probabilities.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 November 2010
  Author:
    Original FORTRAN77 version by Richard Lund.
    C version by John Burkardt.
  Reference:
    PR Freeman,
    Algorithm AS 59:
    Hypergeometric Probabilities,
    Applied Statistics,
    Volume 22, Number 1, 1973, pages 130-133.
    Richard Lund,
    Algorithm AS 152:
    Cumulative hypergeometric probabilities,
    Applied Statistics,
    Volume 29, Number 2, 1980, pages 221-223.
    BL Shea,
    Remark AS R77:
    A Remark on Algorithm AS 152: Cumulative hypergeometric probabilities,
    Applied Statistics,
    Volume 38, Number 1, 1989, pages 199-204.
  Parameters:
    Input, int POINT, is TRUE if the point probability is desired,
    and FALSE if the cumulative probability is desired.
    Input, int KK, the sample size.
    0 <= KK <= MM.
    Input, int LL, the number of successes in the sample.
    0 <= LL <= KK.
    Input, int MM, the population size that was sampled.
    0 <= MM.
    Input, int NN, the number of "successes" in the population.
    0 <= NN <= MM.
    Output, int *IFAULT, error flag.
    0, no error occurred.
    nonzero, an error occurred.
    Output, double CHYPER, the PDF (point probability) of
    exactly LL successes out of KK samples, or the CDF (cumulative
    probability) of up to LL successes out of KK samples.
*/
{
	static ityp result = MAX_VAL;
	
	const b3dt * const s_data = data;
	const bool point = s_data->a0;
	dim_typ kk = s_data->a1;
	dim_typ ll = s_data->a2;
	dim_typ mm = s_data->a3;
	dim_typ nn = s_data->a4; 
	
	ityp arg;
	bool dir = true;
	dim_typ i;
	dim_typ j;
	dim_typ k = kk+1;
	dim_typ kl;
	dim_typ l = ll+1;
	dim_typ m = mm+1;
	ityp mean;
	dim_typ mnkl;
	dim_typ n = nn+1;
	dim_typ nl;
	ityp p;
	ityp pt;
	ityp sig;
	ityp value = 0.00;

	/*
	Check arguments are within permitted limits.
	*/

	if ( n < 1 || m < n || k < 1 || m < k || l < 1 || m-n < k-l)
	{
		result = CHYPER_INVALIDRETURNVALUE;
		return &result;
	}

	if ( !point )
		value = 1.00;

	if(n < l || k < l)
	{
		result = value;
		return &result;
	}

	if ( k == 1 || k == m || n == 1 || n == m )
	{
		result = 1.00;
		return &result;
	}


	if ( !point && ll == MIN ( kk, nn ) )
	{
		result = 1.00;
		return &result;
	}

	value = 1.00;

	p = ( ityp ) ( nn ) / ( ityp ) ( mm - nn );

	if ( 16.00 * MAX ( p, 1.00 / p )< ( ityp ) ( MIN ( kk, mm - kk ) ) && 1000 < mm && - 100.00 < - 88.00 )
	{
		/*
		Use a normal approximation.
		*/
		mean = ( ityp ) ( kk * nn ) / ( ityp ) ( mm );
		sig = sqrt ( mean * ( ( ityp ) ( mm - nn ) / ( ityp ) ( mm ) ) * ( ( ityp ) ( mm - kk ) / ( ( ityp ) ( mm - 1 ) ) ) );
		value = point ? (- 88.00 <= ((arg = - 0.50 * ( pow ( ( ( ityp ) ( ll ) - mean ) / sig, 2 ) ))) ? exp ( arg ) / ( sig * 2.506628274631001 ) : 0.00) :
			alnorm ( ( ( ityp ) ( ll ) + 0.50 - mean ) / sig, 0 );
	}
	else
	{
		/*
		Calculate exact hypergeometric probabilities.
		Interchange K and N if this saves calculations.
		*/
		if ( MIN ( n - 1, m - n ) < MIN ( k - 1, m - k ) )
		{
			i = k;
			k = n;
			n = i;
		}

		if ( m - k < k - 1 )
		{
			dir = !dir;
			l = n - l + 1;
			k = m - k + 1;
		}

		if ( 600 < mm )
		{
			/*
			Take logarithms of factorials.
			*/
			value = - 88.0 <= ((p = alnfac ( nn ) - alnfac ( mm ) + alnfac ( mm - kk ) + alnfac ( kk ) + alnfac ( mm - nn ) - alnfac ( ll ) - alnfac ( nn - ll ) - alnfac ( kk - ll )- alnfac ( mm - nn - kk + ll ))) ? exp(p) : 0.00;
		}
		else
		{
			/*
			Use Freeman/Lund algorithm.
			*/
			for ( i = 1; i <= l - 1; ++i )
				value *= ( ityp ) ( ( k - i ) * ( n - i ) ) / ( ityp ) ( ( l - i ) * ( m - i ) );
			if ( l != k )
			{
				j = m - n + l;
				for ( i = l; i <= k - 1; i++ )
					value *= ( ityp ) ( j - i ) / ( ityp ) ( m - i );
			}
		}

		if ( point )
		{
			result = value;
			return &result;
		}
	
		if ( value == 0.00 )
		{
		/*
		We must recompute the point probability since it has underflowed.
		*/
			if ( mm <= 600 )
				p = alnfac ( nn ) - alnfac ( mm ) + alnfac ( kk ) + alnfac ( mm - nn ) - alnfac ( ll ) - alnfac ( nn - ll ) - alnfac ( kk - ll ) - alnfac ( mm - nn - kk + ll ) + alnfac ( mm - kk );

			p += log ( 1.0E+35 );

			if ( p < - 88.0 )
			{
				result = ( ( ityp ) ( nn * kk + nn + kk + 1 ) / ( ityp ) ( mm + 2 ) < ( ityp ) ( ll ) ) ? 1.00 : value;
				return &result;
			}
			else
				p = exp ( p );
		}
		else
		/*
		Scale up at this point.
		*/
			p = value * 1.0E+35;

		pt = 0.00;
		nl = n - l;
		kl = k - l;
		mnkl = m - n - kl + 1;

		if ( l <= kl )
		{
			for ( i = 1; i <= l - 1; ++i )
			{
				p = p * ( ityp ) ( ( l - i ) * ( mnkl - i ) ) / ( ityp ) ( ( nl + i ) * ( kl + i ) );
				pt += p;
			}
		}
		else
		{
			dir = !dir;
			for ( j = 0; j <= kl - 1; ++j )
			{
				p *= ( ityp ) ( ( nl - j ) * ( kl - j ) )
				/ ( ityp ) ( ( l + j ) * ( mnkl + j ) );
				pt += p;
			}
		}

		value = dir ? value + ( pt / 1.0E+35 ) : 1.00 - ( pt / 1.0E+35 );
	}

	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hypergeometric_cdf_values ( void *data)
/******************************************************************************/
/*
  Purpose:
    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
  Discussion:
    CDF(X)(A,B) is the probability of at most X successes in A trials,
    given that the probability of success on a single trial is B.
    In Mathematica, the function can be evaluated by:
      Needs["Statistics`DiscreteDistributions`]
      dist = HypergeometricDistribution [ sam, suc, pop ]
      CDF [ dist, n ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 September 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
    Daniel Zwillinger,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition, CRC Press, 1996, pages 651-652.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *SAM, int *SUC, int *POP, the sample size,
    success size, and population parameters of the function.
    Output, int *N, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const _5pdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * sam = s_data->a1;
	dim_typ * suc = s_data->a2;
	dim_typ * pop = s_data->a3;
	dim_typ * n = s_data->a4;
	ityp * fx = s_data->a5;
	
    # define N_MAX 16

    ityp fx_vec[N_MAX] =
    {
        0.6001858177500578E-01,
        0.2615284665839845E+00,
        0.6695237889132748E+00,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.5332595856827856E+00,
        0.1819495964117640E+00,
        0.4448047017527730E-01,
        0.9999991751316731E+00,
        0.9926860896560750E+00,
        0.8410799901444538E+00,
        0.3459800113391901E+00,
        0.0000000000000000E+00,
        0.2088888139634505E-02,
        0.3876752992448843E+00,
        0.9135215248834896E+00
    };

    dim_typ n_vec[N_MAX] =
    {
        7,  8,  9, 10,
        6,  6,  6,  6,
        6,  6,  6,  6,
        0,  0,  0,  0
    };

    dim_typ pop_vec[N_MAX] =
    {
        100, 100, 100, 100,
        100, 100, 100, 100,
        100, 100, 100, 100,
        90,  200, 1000, 10000
    };

    dim_typ sam_vec[N_MAX] =
    {
        10, 10, 10, 10,
        6,  7,  8,  9,
        10, 10, 10, 10,
        10, 10, 10, 10
    };

    dim_typ suc_vec[N_MAX] =
    {
        90, 90, 90, 90,
        90, 90, 90, 90,
        10, 30, 50, 70,
        90, 90, 90, 90
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *sam = *suc = *pop = *n = 0;
        *fx = 0.00;
    }
    else
    {
        *sam = sam_vec[*n_data-1];
        *suc = suc_vec[*n_data-1];
        *pop = pop_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _hypergeometric_pdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    HYPERGEOMETRIC_PDF_VALUES returns some values of the hypergeometric PDF.
  Discussion:
    CDF(X)(A,B) is the probability of X successes in A trials,
    given that the probability of success on a single trial is B.
    In Mathematica, the function can be evaluated by:
      dist = HypergeometricDistribution [ sam, suc, pop ]
      PDF [ dist, n ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 January 2008
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
    Daniel Zwillinger,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition, CRC Press, 1996, pages 651-652.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *SAM, int *SUC, int *POP, the sample size,
    success size, and population parameters of the function.
    Output, int *N, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const _5pdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * sam = s_data->a1;
	dim_typ * suc = s_data->a2;
	dim_typ * pop = s_data->a3;
	dim_typ * n = s_data->a4;
	ityp * fx = s_data->a5;
	
    # define N_MAX 16

    ityp fx_vec[N_MAX] =
    {
        0.05179370533242827E+00,
        0.2015098848089788E+00,
        0.4079953223292903E+00,
        0.3304762110867252E+00,
        0.5223047493549780E+00,
        0.3889503452643453E+00,
        0.1505614239732950E+00,
        0.03927689321042477E+00,
        0.00003099828465518108E+00,
        0.03145116093938197E+00,
        0.2114132170316862E+00,
        0.2075776621999210E+00,
        0.0000000000000000E+00,
        0.002088888139634505E+00,
        0.3876752992448843E+00,
        0.9135215248834896E+00
    };

    dim_typ n_vec[N_MAX] =
    {
        7,  8,  9, 10,
        6,  6,  6,  6,
        6,  6,  6,  6,
        0,  0,  0,  0
    };

    dim_typ pop_vec[N_MAX] =
    {
        100, 100, 100, 100,
        100, 100, 100, 100,
        100, 100, 100, 100,
        90,  200, 1000, 10000
    };

    int sam_vec[N_MAX] =
    {
        10, 10, 10, 10,
        6,  7,  8,  9,
        10, 10, 10, 10,
        10, 10, 10, 10
    };

    int suc_vec[N_MAX] =
    {
        90, 90, 90, 90,
        90, 90, 90, 90,
        10, 30, 50, 70,
        90, 90, 90, 90
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *sam = *suc = *pop = *n = 0;
        *fx = 0.00;
    }
    else
    {
        *sam = sam_vec[*n_data-1];
        *suc = suc_vec[*n_data-1];
        *pop = pop_vec[*n_data-1];
        *n = n_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
