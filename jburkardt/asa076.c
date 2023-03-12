#ifndef __DISABLEDEEP_ASA076

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _owen_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    OWEN_VALUES returns some values of Owen's T function.
  Discussion:
    Owen's T function is useful for computation of the bivariate normal
    distribution and the distribution of a skewed normal distribution.
    Although it was originally formulated in terms of the bivariate
    normal function, the function can be defined more directly as
      T(H,A) = 1 / ( 2 * M_PI ) *
        Integral ( 0 <= X <= A ) e^(H^2*(1+X^2)/2) / (1+X^2) dX
    In Mathematica, the function can be evaluated by:
      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    15 December 2011
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
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *H, a parameter.
    Output, double *A, the upper limit of the integral.
    Output, double *T, the value of the function.
*/
{
	const pdt3pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * h = s_data->a1;
	ityp * a = s_data->a2;
	ityp * t = s_data->a3;
	
    # define N_MAX 28

    static ityp a_vec[N_MAX] =
    {
        0.2500000000000000E+00,
        0.4375000000000000E+00,
        0.9687500000000000E+00,
        0.0625000000000000E+00,
        0.5000000000000000E+00,
        0.9999975000000000E+00,
        0.5000000000000000E+00,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.5000000000000000E+00,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.5000000000000000E+00,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.5000000000000000E+00,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.5000000000000000E+00,
        0.1000000000000000E+01,
        0.2000000000000000E+01,
        0.3000000000000000E+01,
        0.1000000000000000E+02,
        0.1000000000000000E+03
    };

    static ityp h_vec[N_MAX] =
    {
        0.0625000000000000E+00,
        6.5000000000000000E+00,
        7.0000000000000000E+00,
        4.7812500000000000E+00,
        2.0000000000000000E+00,
        1.0000000000000000E+00,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.1000000000000000E+01,
        0.5000000000000000E+00,
        0.5000000000000000E+00,
        0.5000000000000000E+00,
        0.5000000000000000E+00,
        0.2500000000000000E+00,
        0.2500000000000000E+00,
        0.2500000000000000E+00,
        0.2500000000000000E+00,
        0.1250000000000000E+00,
        0.1250000000000000E+00,
        0.1250000000000000E+00,
        0.1250000000000000E+00,
        0.7812500000000000E-02,
        0.7812500000000000E-02,
        0.7812500000000000E-02,
        0.7812500000000000E-02,
        0.7812500000000000E-02,
        0.7812500000000000E-02
    };

    static ityp t_vec[N_MAX] =
    {
        3.8911930234701366E-02,
        2.0005773048508315E-11,
        6.3990627193898685E-13,
        1.0632974804687463E-07,
        8.6250779855215071E-03,
        6.6741808978228592E-02,
        0.4306469112078537E-01,
        0.6674188216570097E-01,
        0.7846818699308410E-01,
        0.7929950474887259E-01,
        0.6448860284750376E-01,
        0.1066710629614485E+00,
        0.1415806036539784E+00,
        0.1510840430760184E+00,
        0.7134663382271778E-01,
        0.1201285306350883E+00,
        0.1666128410939293E+00,
        0.1847501847929859E+00,
        0.7317273327500385E-01,
        0.1237630544953746E+00,
        0.1737438887583106E+00,
        0.1951190307092811E+00,
        0.7378938035365546E-01,
        0.1249951430754052E+00,
        0.1761984774738108E+00,
        0.1987772386442824E+00,
        0.2340886964802671E+00,
        0.2479460829231492E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *h = *a = *t = 0.00;
    }
    else
    {
        *h = h_vec[*n_data-1];
        *a = a_vec[*n_data-1];
        *t = t_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tha ( void * data)
/******************************************************************************/
/*
  Purpose:
    THA computes Owen's T function.
  Discussion:
    This function computes T(H1/H2, A1/A2) for any real numbers H1, H2,
    A1 and A2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 November 2010
  Author:
    Original FORTRAN77 version by JC Young, Christoph Minder.
    C version by John Burkardt.
  Reference:
    Richard Boys,
    Remark AS R80:
    A Remark on Algorithm AS76:
    An Integral Useful in Calculating Noncentral T and Bivariate
    Normal Probabilities,
    Applied Statistics,
    Volume 38, Number 3, 1989, pages 580-582.
    Youn-Min Chou,
    Remark AS R55:
    A Remark on Algorithm AS76:
    An Integral Useful in Calculating Noncentral T and Bivariate
    Normal Probabilities,
    Applied Statistics,
    Volume 34, Number 1, 1985, pages 100-101.
    PW Goedhart, MJW Jansen,
    Remark AS R89:
    A Remark on Algorithm AS76:
    An Integral Useful in Calculating Noncentral T and Bivariate
    Normal Probabilities,
    Applied Statistics,
    Volume 41, Number 2, 1992, pages 496-497.
    JC Young, Christoph Minder,
    Algorithm AS 76:
    An Algorithm Useful in Calculating Noncentral T and
    Bivariate Normal Distributions,
    Applied Statistics,
    Volume 23, Number 3, 1974, pages 455-457.
  Parameters:
    Input, double H1, H2, A1, A2, define the arguments
    of the T function.
    Output, double THA, the value of Owen's T function.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp h1 = a_data[0];
	const register ityp h2 = a_data[1];
	const register ityp a1 = a_data[2];
	const register ityp a2 = a_data[3];
	
	ityp a;
	ityp absa;
	ityp ah;
	ityp c1;
	ityp c2;
	ityp ex;
	ityp g;
	ityp gah;
	ityp gh;
	ityp h;
	ityp lam;
	ityp value;

	if ( h2 == 0.00 )
	{
		result = 0.00;
		return &result;
	}

	h = h1 / h2;

	if ( a2 == 0.00 )
	{
		result = (h<0.00 ? 2.00 : ( 1.00 - alnorm ( h, 0 ) ) / 2.00)*(1-((a1<0.00)<<1));
		return &result;
	}

	a = a1 / a2;

	if ( fabs ( h ) < 0.30 && 7.00 < fabs ( a ) )
	{
		lam = fabs ( a * h );
		ex = exp ( - lam * lam / 2.00 );
		g = alnorm ( lam, 0 );
		c1 = ( ex / lam + sqrt ( 6.2831853071795864769 ) * ( g - 0.5 ) ) / 6.2831853071795864769;
		c2 = ( ( lam * lam + 2.0 ) * ex / lam / lam / lam + sqrt ( 6.2831853071795864769 ) * ( g - 0.5 ) ) / ( 6.0 * 6.2831853071795864769 );
		ah = fabs ( h );
		value = 0.25 - c1 * ah + c2 * ah * ah * ah;
		value = fabs(value)*(1-((a<0.00)<<1));
	}
	else
	{
		absa = fabs ( a );

		if ( absa <= 1.0 )
		{
			result = tfn ( h, a );
			return &result;
		}

		ah = absa * h;
		gh = alnorm ( h, 0 );
		gah = alnorm ( ah, 0 );
		value = (0.5 * ( gh + gah ) - gh * gah - tfn ( ah, 1.0 / absa ))*(1-((a<0.00)<<1));
	}

	result = value;
	return &result;
}

#endif
