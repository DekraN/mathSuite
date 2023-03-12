#ifndef __DISABLEDEEP_PDFLIB

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_gamma_log ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAMMA_LOG evaluates the logarithm of the gamma function.
  Discussion:
    This routine calculates the LOG(GAMMA) function for a positive real
    argument X.  Computation is based on an algorithm outlined in
    references 1 and 2.  The program uses rational functions that
    theoretically approximate LOG(GAMMA) to at least 18 significant
    decimal digits.  The approximation for X > 12 is from reference
    3, while approximations for X < 12.0 are similar to those in
    reference 1, but are unpublished.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 April 2013
  Author:
    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.
  Reference:
    William Cody, Kenneth Hillstrom,
    Chebyshev Approximations for the Natural Logarithm of the
    Gamma Function,
    Mathematics of Computation,
    Volume 21, Number 98, April 1967, pages 198-203.
    Kenneth Hillstrom,
    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
    May 1969.
    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.
  Parameters:
    Input, double X, the argument of the function.
    0.0 < X
    Output, double R8_GAMMA_LOG, the value of the function.
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
	ityp corr;
	const ityp d1 = -5.772156649015328605195174E-01;
	const ityp d2 = 4.227843350984671393993777E-01;
	const ityp d4 = 1.791759469228055000094023;
	const ityp frtbig = 2.25E+76;
	dim_typ i;
	ityp p1[8] = 
	{
		4.945235359296727046734888, 
		2.018112620856775083915565E+02, 
		2.290838373831346393026739E+03, 
		1.131967205903380828685045E+04, 
		2.855724635671635335736389E+04, 
		3.848496228443793359990269E+04, 
		2.637748787624195437963534E+04, 
		7.225813979700288197698961E+03 
	};
	ityp p2[8] = 
	{ 
		4.974607845568932035012064, 
		5.424138599891070494101986E+02, 
		1.550693864978364947665077E+04, 
		1.847932904445632425417223E+05, 
		1.088204769468828767498470E+06, 
		3.338152967987029735917223E+06, 
		5.106661678927352456275255E+06, 
		3.074109054850539556250927E+06
	};
	ityp p4[8] = 
	{
		1.474502166059939948905062E+04, 
		2.426813369486704502836312E+06, 
		1.214755574045093227939592E+08, 
		2.663432449630976949898078E+09, 
		2.940378956634553899906876E+10, 
		1.702665737765398868392998E+11, 
		4.926125793377430887588120E+11, 
		5.606251856223951465078242E+11 
	};
	ityp q1[8] = 
	{ 
		6.748212550303777196073036E+01, 
		1.113332393857199323513008E+03, 
		7.738757056935398733233834E+03, 
		2.763987074403340708898585E+04, 
		5.499310206226157329794414E+04, 
		6.161122180066002127833352E+04, 
		3.635127591501940507276287E+04, 
		8.785536302431013170870835E+03 
	};
	ityp q2[8] = 
	{ 
		1.830328399370592604055942E+02, 
		7.765049321445005871323047E+03, 
		1.331903827966074194402448E+05, 
		1.136705821321969608938755E+06, 
		5.267964117437946917577538E+06, 
		1.346701454311101692290052E+07, 
		1.782736530353274213975932E+07, 
		9.533095591844353613395747E+06 
	};
	ityp q4[8] = 
	{ 
		2.690530175870899333379843E+03, 
		6.393885654300092398984238E+05, 
		4.135599930241388052042842E+07, 
		1.120872109616147941376570E+09, 
		1.488613728678813811542398E+10, 
		1.016803586272438228077304E+11, 
		3.417476345507377132798597E+11, 
		4.463158187419713286462081E+11 
	};
	ityp res;
	const register ityp sqrtpi = 0.9189385332046727417803297;
	const ityp xbig = 2.55E+305;
	ityp xden;
	const ityp xinf = 1.79E+308;
	ityp xm1;
	ityp xm2;
	ityp xm4;
	ityp xnum;
	ityp y;
	ityp ysq;
	
	y = x;
	
	if ( 0.0 < y && y <= xbig )
	{
		if ( y <= r8_epsilon ( ) )
			res = - log ( y );

		/*
		EPS < X <= 1.5.
		*/
		else if ( y <= 1.5 )
		{
			if ( y < 0.6796875 )
			{
				corr = -log ( y );
				xm1 = y;
			}
			else
			{
				corr = 0.00;
				xm1 = ( y - 0.50 ) - 0.50;
			}
	
			if ( y <= 0.50 || 0.6796875 <= y )
			{
				xden = 1.00;
				xnum = 0.00;
				#pragma omp parallel for num_threads(8)
				for ( i = 0; i < 8; ++i )
				{
					xnum = xnum * xm1 + p1[i];
					xden = xden * xm1 + q1[i];
				}
				res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
				}
			else
			{
				xm2 = ( y - 0.50 ) - 0.50;
				xden = 1.00;
				xnum = 0.00;
				#pragma omp parallel for num_threads(8)
				for ( i = 0; i < 8; ++i )
				{
					xnum = xnum * xm2 + p2[i];
					xden = xden * xm2 + q2[i];
				}
				res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );
			}
		}
		/*
		1.5 < X <= 4.0.
		*/
		else if ( y <= 4.00 )
		{
			xm2 = y - 2.00;
			xden = 1.00;
			xnum = 0.00;
			#pragma omp parallel for num_threads(8)
			for ( i = 0; i < 8; ++i )
			{
				xnum = xnum * xm2 + p2[i];
				xden = xden * xm2 + q2[i];
			}
			res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
		}
		/*
		4.0 < X <= 12.0.
		*/
		else if ( y <= 12.00 )
		{
			xm4 = y - 4.00;
			xden = -1.00;
			xnum = 0.00;
			#pragma omp parallel for num_threads(8)
			for ( i = 0; i < 8; ++i )
			{
				xnum = xnum * xm4 + p4[i];
				xden = xden * xm4 + q4[i];
			}
			res = d4 + xm4 * ( xnum / xden );
		}
		/*
		Evaluate for 12 <= argument.
		*/
		else
		{
			res = 0.00;
	
			if ( y <= frtbig )
			{
				res = c[6];
				ysq = y * y;
				#pragma omp parallel for num_threads(6)
				for ( i = 0; i < 6; ++i )
					res = res / ysq + c[i];
			}
			res /= y;
			corr = log ( y );
			res += sqrtpi - 0.50 * corr;
			res += y * ( corr - 1.00 );
		}
	}
	/*
	Return for bad arguments.
	*/
	else
		res = xinf;
	/*
	Final adjustments and return.
	*/
	
	result = res;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_chi_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CHI_PDF evaluates the PDF of a chi-squared distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2013
  Author:
    Original FORTRAN90 version by Guannan Zhang.
    C version by John Burkardt.
  Parameters:
    Input, double DF, the degrees of freedom.
    0.0 < DF.
    Input, double RVAL, the point where the PDF is evaluated.
    Output, double R8_CHI_PDF, the value of the PDF at RVAL.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp df = a_data[0];
	const register ityp rval = a_data[1];
	
    ityp temp2;
    ityp value;

    if ( df <= 0.0 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( rval <= 0.00 )
        value = 0.00;
    else
    {
        temp2 = df * 0.50;
        value = exp ( temp2 - 1.00 ) * log ( rval ) - 0.50 * rval - temp2 * log ( 2.00 ) - r8_gamma_log ( temp2 );
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_exponential_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_EXPONENTIAL_PDF evaluates the PDF of an exponential distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2013
  Author:
    Original FORTRAN90 version by Guannan Zhang.
    C version by John Burkardt.
  Parameters:
    Input, double BETA, the scale value.
    0.0 < BETA.
    Input, double RVAL, the point where the PDF is evaluated.
    Output, double R8_EXPONENTIAL_PDF, the value of the PDF at RVAL.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp beta = a_data[0];
	const register ityp rval = a_data[1];
	
	result = beta <= 0.00 ? MAX_VAL : 0.00 + (rval>=0)*exp ( - rval / beta ) / beta;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_exponential_01_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_EXPONENTIAL_01_PDF: PDF of the standard exponential distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 June 2013
  Author:
    John Burkardt.
  Parameters:
    Input, double RVAL, the point where the PDF is evaluated.
    Output, double R8_EXPONENTIAL_01_PDF, the value of the PDF at RVAL.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp rval = *(ityp *) data;
	
	result = 0.00 + (rval>=0)*exp ( - rval );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_invchi_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_INVCHI_PDF evaluates the PDF of an inverse chi-squared distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2013
  Author:
    Original FORTRAN90 version by Guannan Zhang.
    C version by John Burkardt.
  Parameters:
    Input, double DF, the degrees of freedom.
    0.0 < DF.
    Input, double RVAL, the point where the PDF is evaluated.
    Output, double R8_INVCHI_PDF, the value of the PDF at RVAL.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp df = a_data[0];
	const register ityp rval = a_data[1];
	
    ityp temp2;
    ityp value;

    if ( df <= 0.00 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( rval <= 0.0 )
        value = 0.0;
    else
    {
        temp2 = df * 0.50;
        value = exp ( - temp2 * log ( 2.00 ) - ( temp2 + 1.00 ) * log ( rval ) - 0.50 / rval - r8_gamma_log ( temp2 ) );
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_normal_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_NORMAL_PDF evaluates the PDF of a normal distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2013
  Author:
    Original FORTRAN90 version by Guannan Zhang.
    C version by John Burkardt.
  Parameters:
    Input, double AV, the mean value.
    Input, double SD, the standard deviation.
    0.0 < SD.
    Input, double RVAL, the point where the PDF is evaluated.
    Output, double R8_NORMAL_PDF, the value of the PDF at RVAL.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp av = a_data[0];
	const register ityp sd = a_data[1];
	const register ityp rval = a_data[2];
	
	result = sd <= 0.00 ? MAX_VAL : exp ( - ( rval - av ) * ( rval - av ) * 0.50 / ( sd * sd ) ) / sd / sqrt ( M_2TPI );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_scinvchi_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_SCINVCHI_PDF: PDF for a scaled inverse chi-squared distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2013
  Author:
    Original FORTRAN90 version by Guannan Zhang.
    C version by John Burkardt.
  Parameters:
    Input, double DF, the degrees of freedom.
    0.0 < DF.
    Input, double S, the scale factor.
    0.0 < S.
    Input, double RVAL, the point where the PDF is evaluated.
    Output, double R8_SCINVCHI_PDF, the value of the PDF at RVAL.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp df = a_data[0];
	const register ityp s = a_data[1];
	const register ityp rval = a_data[2];
	
    ityp temp2;
    ityp value;

    if ( df <= 0.0 || s <= 0.00)
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( rval <= 0.00 )
        value = 0.00;
    else
    {
        temp2 = df * 0.50;
        value = exp ( temp2 * log ( temp2 ) + temp2 * log ( s ) - ( temp2 * s / rval ) - ( temp2 + 1.0 ) * log ( rval ) - r8_gamma_log ( temp2 ) );
    }

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_uniform_pdf ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_UNIFORM_PDF evaluates the PDF of a uniform distribution.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 April 2013
  Author:
    Original FORTRAN90 version by Guannan Zhang.
    C version by John Burkardt.
  Parameters:
    Input, double LOWER, UPPER, the lower and upper range limits.
    LOWER < UPPER.
    Input, double RVAL, the point where the PDF is evaluated.
    Output, double R8_UNIFORM_PDF, the value of the PDF at RVAL.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp lower = a_data[0];
	const register ityp upper = a_data[1];
	const register ityp rval = a_data[2];
	
    ityp value;

    if ( upper <= lower )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( rval < lower )
        value = 0.00;
    else if ( rval <= upper )
        value = 1.00 / ( upper - lower );
    else
        value = 0.00;

	result = value;
    return &result;
}

#endif
