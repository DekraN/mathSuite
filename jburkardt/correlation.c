#ifndef __DISABLEDEEP_CORRELATION

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_besselj ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_BESSELJ evaluates the Bessel J correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    ityp rhohat;

    for (dim_typ i = 0; i < n; ++i )
    {
        rhohat = abs ( rho[i] ) / rho0;
        c[i] = r8_besj0 ( rhohat );
    }

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_besselk ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_BESSELK evaluates the Bessel K correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    ityp rhohat;

    for (dim_typ i = 0; i < n; ++i )
    {
        if ( rho[i] == 0.00 )
            c[i] = 1.00;
        else
        {
            rhohat = abs ( rho[i] ) / rho0;
            c[i] = rhohat * r8_besk1 ( rhohat );
        }
    }

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_constant ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_CONSTANT evaluates the constant correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho1 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        c[i] = 1.00;
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_cubic ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_CUBIC evaluates the cubic correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    double rhohat;

    for (dim_typ i = 0; i < n; ++i )
    {
        rhohat = MIN ( abs ( rho[i] ) / rho0, 1.00 );
        c[i] = 1.00 - 7.00  * pow ( rhohat, 2 ) + 8.75 * pow ( rhohat, 3 )- 3.50  * pow ( rhohat, 5 )+ 0.750 * pow ( rhohat, 7 );
    }

    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_damped_cosine ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_DAMPED_COSINE evaluates the damped cosine correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        c[i] = exp ( - abs ( rho[i] ) / rho0 ) * cos ( abs ( rho[i] ) / rho0 );
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_damped_sine ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_DAMPED_SINE evaluates the damped sine correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    double rhohat;

    for (dim_typ i = 0; i < n; ++i )
    {
        if ( rho[i] == 0.00 )
            c[i] = 1.00;
        else
        {
            rhohat = abs ( rho[i] ) / rho0;
            c[i] = sin ( rhohat ) / rhohat;
        }
    }
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_exponential ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_EXPONENTIAL evaluates the exponential correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i)
        c[i] = exp ( - abs ( rho[i] ) / rho0 );
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_gaussian ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_GAUSSIAN evaluates the Gaussian correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c;

    for (dim_typ i = 0; i < n; ++i )
        c[i] = exp ( - pow ( rho[i] / rho0, 2 ) );
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_hole ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_HOLE evaluates the hole correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i)
        c[i] = ( 1.00 - abs ( rho[i] ) / rho0 ) * exp ( - abs ( rho[i] ) / rho0 );
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_linear ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_LINEAR evaluates the linear correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        c[i] = rho0 < abs ( rho[i]) ? 0.0 : ( rho0 - abs ( rho[i] ) ) / rho0;
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_matern ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_MATERN evaluates the Matern correlation function.
  Discussion:
    In order to call this routine under a dummy name, I had to drop NU from
    the parameter list.
    The Matern correlation is
      rho1 = 2 * sqrt ( nu ) * rho / rho0
      c(rho) = ( rho1 )^nu * BesselK ( nu, rho1 )
               / gamma ( nu ) / 2 ^ ( nu - 1 )
    The Matern covariance has the form:
      K(rho) = sigma^2 * c(rho)
    A Gaussian process with Matern covariance has sample paths that are
    differentiable (nu - 1) times.
    When nu = 0.5, the Matern covariance is the exponential covariance.
    As nu goes to +oo, the correlation converges to exp ( - (rho/rho0)^2 ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 November 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    0.0 <= RHO.
    Input, double RHO0, the correlation length.
    0.0 < RHO0.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    ityp nu = 2.50;
	ityp rho1;
	
    for (dim_typ i = 0; i < n; ++i )
    {
        rho1 = 2.00 * sqrt ( nu ) * abs ( rho[i] ) / rho0;
        c[i] = rho1 == 0.00 ? 1.00 : pow ( rho1, nu ) * r8_besk ( nu, rho1 ) / r8_gamma ( nu ) / pow ( 2.00, nu - 1.00 );
    }
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_pentaspherical ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_PENTASPHERICAL evaluates the pentaspherical correlation function.
  Discussion:
    This correlation is based on the volume of overlap of two spheres
    of radius RHO0 and separation RHO.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    ityp rhohat;
    for (dim_typ i = 0; i < n; ++i )
    {
        rhohat = MIN ( abs ( rho[i] ) / rho0, 1.0 );
        c[i] = 1.0 - 1.875 * rhohat + 1.25 * pow ( rhohat, 3 )- 0.375 * pow ( rhohat, 5 );
    }
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_power ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_POWER evaluates the power correlation function.
  Discussion:
    In order to be able to call this routine under a dummy name, I had
    to drop E from the argument list.
    The power correlation is
      C(rho) = ( 1 - |rho| )^e  if 0 <= |rho| <= 1
             = 0                otherwise
      The constraint on the exponent is 2 <= e.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    0.0 <= RHO.
    Input, double RHO0, the correlation length.
    0.0 < RHO0.
    Input, double E, the exponent.
    E has a default value of 2.0;
    2.0 <= E.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    ityp e = 2.00;
    ityp rhohat;
    for (dim_typ i = 0; i < n; ++i )
    {
        rhohat = abs ( rho[i] ) / rho0;
        c[i] = 0.00 + (rhohat<=1.00)*pow ( 1.0 - rhohat, e );
    }
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_rational_quadratic ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_RATIONAL_QUADRATIC: rational quadratic correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        c[i] = 1.00 / ( 1.00 + pow ( rho[i] / rho0, 2 ) );
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_spherical ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_SPHERICAL evaluates the spherical correlation function.
  Discussion:
    This correlation is based on the volume of overlap of two spheres
    of radius RHO0 and separation RHO.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho0 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    ityp rhohat;
    for (dim_typ i = 0; i < n; ++i)
    {
        rhohat = MIN ( abs ( rho[i] ) / rho0, 1.00 );
        c[i] = 1.00 - 1.50 * rhohat + 0.50 * pow ( rhohat, 3 );
    }
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _correlation_white_noise ( void * data)
/******************************************************************************/
/*
  Purpose:
    CORRELATION_WHITE_NOISE evaluates the white noise correlation function.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 November 2012
  Author:
    John Burkardt
  Reference:
    Petter Abrahamsen,
    A Review of Gaussian Random Fields and Correlation Functions,
    Norwegian Computing Center, 1997.
  Parameters:
    Input, int N, the number of arguments.
    Input, double RHO[N], the arguments.
    Input, double RHO0, the correlation length.
    Output, double C[N], the correlations.
*/
{
	const dtpitit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * rho = s_data->a1;
	const register ityp rho1 = s_data->a2;
	
    ityp *c = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    for (dim_typ i = 0; i < n; ++i )
        c[i] = rho[i] = 0.00;
    return c;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_b0mp ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_B0MP evaluates the modulus and phase for the Bessel J0 and Y0 functions.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double *AMPL, *THETA, the modulus and phase.
*/
{
	const it2pit * const s_data = data;
	const register ityp x = s_data->a0;
	ityp * ampl = s_data->a1;
	ityp * theta = s_data->a2;
	
    static ityp bm0cs[37] =
    {
        +0.9211656246827742712573767730182E-01,
        -0.1050590997271905102480716371755E-02,
        +0.1470159840768759754056392850952E-04,
        -0.5058557606038554223347929327702E-06,
        +0.2787254538632444176630356137881E-07,
        -0.2062363611780914802618841018973E-08,
        +0.1870214313138879675138172596261E-09,
        -0.1969330971135636200241730777825E-10,
        +0.2325973793999275444012508818052E-11,
        -0.3009520344938250272851224734482E-12,
        +0.4194521333850669181471206768646E-13,
        -0.6219449312188445825973267429564E-14,
        +0.9718260411336068469601765885269E-15,
        -0.1588478585701075207366635966937E-15,
        +0.2700072193671308890086217324458E-16,
        -0.4750092365234008992477504786773E-17,
        +0.8615128162604370873191703746560E-18,
        -0.1605608686956144815745602703359E-18,
        +0.3066513987314482975188539801599E-19,
        -0.5987764223193956430696505617066E-20,
        +0.1192971253748248306489069841066E-20,
        -0.2420969142044805489484682581333E-21,
        +0.4996751760510616453371002879999E-22,
        -0.1047493639351158510095040511999E-22,
        +0.2227786843797468101048183466666E-23,
        -0.4801813239398162862370542933333E-24,
        +0.1047962723470959956476996266666E-24,
        -0.2313858165678615325101260800000E-25,
        +0.5164823088462674211635199999999E-26,
        -0.1164691191850065389525401599999E-26,
        +0.2651788486043319282958336000000E-27,
        -0.6092559503825728497691306666666E-28,
        +0.1411804686144259308038826666666E-28,
        -0.3298094961231737245750613333333E-29,
        +0.7763931143074065031714133333333E-30,
        -0.1841031343661458478421333333333E-30,
        +0.4395880138594310737100799999999E-31
    };
    static ityp bm02cs[40] =
    {
        +0.9500415145228381369330861335560E-01,
        -0.3801864682365670991748081566851E-03,
        +0.2258339301031481192951829927224E-05,
        -0.3895725802372228764730621412605E-07,
        +0.1246886416512081697930990529725E-08,
        -0.6065949022102503779803835058387E-10,
        +0.4008461651421746991015275971045E-11,
        -0.3350998183398094218467298794574E-12,
        +0.3377119716517417367063264341996E-13,
        -0.3964585901635012700569356295823E-14,
        +0.5286111503883857217387939744735E-15,
        -0.7852519083450852313654640243493E-16,
        +0.1280300573386682201011634073449E-16,
        -0.2263996296391429776287099244884E-17,
        +0.4300496929656790388646410290477E-18,
        -0.8705749805132587079747535451455E-19,
        +0.1865862713962095141181442772050E-19,
        -0.4210482486093065457345086972301E-20,
        +0.9956676964228400991581627417842E-21,
        -0.2457357442805313359605921478547E-21,
        +0.6307692160762031568087353707059E-22,
        -0.1678773691440740142693331172388E-22,
        +0.4620259064673904433770878136087E-23,
        -0.1311782266860308732237693402496E-23,
        +0.3834087564116302827747922440276E-24,
        -0.1151459324077741271072613293576E-24,
        +0.3547210007523338523076971345213E-25,
        -0.1119218385815004646264355942176E-25,
        +0.3611879427629837831698404994257E-26,
        -0.1190687765913333150092641762463E-26,
        +0.4005094059403968131802476449536E-27,
        -0.1373169422452212390595193916017E-27,
        +0.4794199088742531585996491526437E-28,
        -0.1702965627624109584006994476452E-28,
        +0.6149512428936330071503575161324E-29,
        -0.2255766896581828349944300237242E-29,
        +0.8399707509294299486061658353200E-30,
        -0.3172997595562602355567423936152E-30,
        +0.1215205298881298554583333026514E-30,
        -0.4715852749754438693013210568045E-31
    };
    static ityp bt02cs[39] =
    {
        -0.24548295213424597462050467249324,
        +0.12544121039084615780785331778299E-02,
        -0.31253950414871522854973446709571E-04,
        +0.14709778249940831164453426969314E-05,
        -0.99543488937950033643468850351158E-07,
        +0.85493166733203041247578711397751E-08,
        -0.86989759526554334557985512179192E-09,
        +0.10052099533559791084540101082153E-09,
        -0.12828230601708892903483623685544E-10,
        +0.17731700781805131705655750451023E-11,
        -0.26174574569485577488636284180925E-12,
        +0.40828351389972059621966481221103E-13,
        -0.66751668239742720054606749554261E-14,
        +0.11365761393071629448392469549951E-14,
        -0.20051189620647160250559266412117E-15,
        +0.36497978794766269635720591464106E-16,
        -0.68309637564582303169355843788800E-17,
        +0.13107583145670756620057104267946E-17,
        -0.25723363101850607778757130649599E-18,
        +0.51521657441863959925267780949333E-19,
        -0.10513017563758802637940741461333E-19,
        +0.21820381991194813847301084501333E-20,
        -0.46004701210362160577225905493333E-21,
        +0.98407006925466818520953651199999E-22,
        -0.21334038035728375844735986346666E-22,
        +0.46831036423973365296066286933333E-23,
        -0.10400213691985747236513382399999E-23,
        +0.23349105677301510051777740800000E-24,
        -0.52956825323318615788049749333333E-25,
        +0.12126341952959756829196287999999E-25,
        -0.28018897082289428760275626666666E-26,
        +0.65292678987012873342593706666666E-27,
        -0.15337980061873346427835733333333E-27,
        +0.36305884306364536682359466666666E-28,
        -0.86560755713629122479172266666666E-29,
        +0.20779909972536284571238399999999E-29,
        -0.50211170221417221674325333333333E-30,
        +0.12208360279441714184191999999999E-30,
        -0.29860056267039913454250666666666E-31
    };
    static ityp bth0cs[44] =
    {
        -0.24901780862128936717709793789967,
        +0.48550299609623749241048615535485E-03,
        -0.54511837345017204950656273563505E-05,
        +0.13558673059405964054377445929903E-06,
        -0.55691398902227626227583218414920E-08,
        +0.32609031824994335304004205719468E-09,
        -0.24918807862461341125237903877993E-10,
        +0.23449377420882520554352413564891E-11,
        -0.26096534444310387762177574766136E-12,
        +0.33353140420097395105869955014923E-13,
        -0.47890000440572684646750770557409E-14,
        +0.75956178436192215972642568545248E-15,
        -0.13131556016891440382773397487633E-15,
        +0.24483618345240857495426820738355E-16,
        -0.48805729810618777683256761918331E-17,
        +0.10327285029786316149223756361204E-17,
        -0.23057633815057217157004744527025E-18,
        +0.54044443001892693993017108483765E-19,
        -0.13240695194366572724155032882385E-19,
        +0.33780795621371970203424792124722E-20,
        -0.89457629157111779003026926292299E-21,
        +0.24519906889219317090899908651405E-21,
        -0.69388422876866318680139933157657E-22,
        +0.20228278714890138392946303337791E-22,
        -0.60628500002335483105794195371764E-23,
        +0.18649748964037635381823788396270E-23,
        -0.58783732384849894560245036530867E-24,
        +0.18958591447999563485531179503513E-24,
        -0.62481979372258858959291620728565E-25,
        +0.21017901684551024686638633529074E-25,
        -0.72084300935209253690813933992446E-26,
        +0.25181363892474240867156405976746E-26,
        -0.89518042258785778806143945953643E-27,
        +0.32357237479762298533256235868587E-27,
        -0.11883010519855353657047144113796E-27,
        +0.44306286907358104820579231941731E-28,
        -0.16761009648834829495792010135681E-28,
        +0.64292946921207466972532393966088E-29,
        -0.24992261166978652421207213682763E-29,
        +0.98399794299521955672828260355318E-30,
        -0.39220375242408016397989131626158E-30,
        +0.15818107030056522138590618845692E-30,
        -0.64525506144890715944344098365426E-31,
        +0.26611111369199356137177018346367E-31
    };
    ityp eta;
    static dim_typ nbm0 = 0;
    static dim_typ nbm02 = 0;
    static dim_typ nbt02 = 0;
    static dim_typ nbth0 = 0;
    static ityp pi4 = 0.785398163397448309615660845819876;
    static ityp xmax = 0.0;
    ityp z;

    if(!nbm0)
    {
        eta = 0.10 * r8_mach ( 3 );
        nbm0 = r8_inits ( bm0cs, 37, eta );
        nbt02 = r8_inits ( bt02cs, 39, eta );
        nbm02 = r8_inits ( bm02cs, 40, eta );
        nbth0 = r8_inits ( bth0cs, 44, eta );
        xmax = 1.00 / r8_mach ( 4 );
    }

    if ( x < 4.00 )
        return NULL;
    else if ( x <= 8.00 )
    {
        z = ( 128.00 / x / x - 5.00 ) / 3.00;
        *ampl = ( 0.75 + r8_csevl ( z, bm0cs, nbm0 ) ) / sqrt ( x );
        *theta = x - pi4 + r8_csevl ( z, bt02cs, nbt02 ) / x;
    }
    else
    {
        z = 128.00 / x / x - 1.00;
        *ampl = ( 0.75 + r8_csevl ( z, bm02cs, nbm02) ) / sqrt ( x );
        *theta = x - pi4 + r8_csevl ( z, bth0cs, nbth0 ) / x;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_besi1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BESI1 evaluates the Bessel function I of order 1 of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BESI1, the Bessel function I of order 1 of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
	static ityp bi1cs[17] = 
	{
		-0.19717132610998597316138503218149E-02,
		+0.40734887667546480608155393652014,
		+0.34838994299959455866245037783787E-01,
		+0.15453945563001236038598401058489E-02,
		+0.41888521098377784129458832004120E-04,
		+0.76490267648362114741959703966069E-06,
		+0.10042493924741178689179808037238E-07,
		+0.99322077919238106481371298054863E-10,
		+0.76638017918447637275200171681349E-12,
		+0.47414189238167394980388091948160E-14,
		+0.24041144040745181799863172032000E-16,
		+0.10171505007093713649121100799999E-18,
		+0.36450935657866949458491733333333E-21,
		+0.11205749502562039344810666666666E-23,
		+0.29875441934468088832000000000000E-26,
		+0.69732310939194709333333333333333E-29,
		+0.14367948220620800000000000000000E-31 
	};
	static int nti1 = 0;
	ityp  value;
	static ityp xmax = 0.00;
	static ityp xmin = 0.00;
	static ityp xsml = 0.00;
	ityp y;
	
	if ( nti1 == 0 )
	{
		nti1 = r8_inits ( bi1cs, 17, 0.10 * r8_mach ( 3 ) );
		xmin = 2.00 * r8_mach ( 1 );
		xsml = sqrt ( 8.00 * r8_mach ( 3 ) );
		xmax = log ( r8_mach ( 2 ) );
	}
	
	y = abs ( x );
	
	if ( y <= xmin )
		value = 0.00;
	else if ( y <= xsml )
		value = 0.50 * x;
	else if ( y <= 3.00 )
		value = x * ( 0.875 + r8_csevl ( y * y / 4.50 - 1.00, bi1cs, nti1 ) );
	else if ( y <= xmax )
		value = exp ( y ) * r8_besi1e ( x );
	else
	{
		result = MAX_VAL;
		return &result;
	}
		
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_besi1e ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BESI1E evaluates the exponentially scaled Bessel function I1(X).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BESI1E, the exponentially scaled Bessel
    function I1(X).
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp ai12cs[69] =
    {
        +0.2857623501828012047449845948469E-01,
        -0.9761097491361468407765164457302E-02,
        -0.1105889387626237162912569212775E-03,
        -0.3882564808877690393456544776274E-05,
        -0.2512236237870208925294520022121E-06,
        -0.2631468846889519506837052365232E-07,
        -0.3835380385964237022045006787968E-08,
        -0.5589743462196583806868112522229E-09,
        -0.1897495812350541234498925033238E-10,
        +0.3252603583015488238555080679949E-10,
        +0.1412580743661378133163366332846E-10,
        +0.2035628544147089507224526136840E-11,
        -0.7198551776245908512092589890446E-12,
        -0.4083551111092197318228499639691E-12,
        -0.2101541842772664313019845727462E-13,
        +0.4272440016711951354297788336997E-13,
        +0.1042027698412880276417414499948E-13,
        -0.3814403072437007804767072535396E-14,
        -0.1880354775510782448512734533963E-14,
        +0.3308202310920928282731903352405E-15,
        +0.2962628997645950139068546542052E-15,
        -0.3209525921993423958778373532887E-16,
        -0.4650305368489358325571282818979E-16,
        +0.4414348323071707949946113759641E-17,
        +0.7517296310842104805425458080295E-17,
        -0.9314178867326883375684847845157E-18,
        -0.1242193275194890956116784488697E-17,
        +0.2414276719454848469005153902176E-18,
        +0.2026944384053285178971922860692E-18,
        -0.6394267188269097787043919886811E-19,
        -0.3049812452373095896084884503571E-19,
        +0.1612841851651480225134622307691E-19,
        +0.3560913964309925054510270904620E-20,
        -0.3752017947936439079666828003246E-20,
        -0.5787037427074799345951982310741E-22,
        +0.7759997511648161961982369632092E-21,
        -0.1452790897202233394064459874085E-21,
        -0.1318225286739036702121922753374E-21,
        +0.6116654862903070701879991331717E-22,
        +0.1376279762427126427730243383634E-22,
        -0.1690837689959347884919839382306E-22,
        +0.1430596088595433153987201085385E-23,
        +0.3409557828090594020405367729902E-23,
        -0.1309457666270760227845738726424E-23,
        -0.3940706411240257436093521417557E-24,
        +0.4277137426980876580806166797352E-24,
        -0.4424634830982606881900283123029E-25,
        -0.8734113196230714972115309788747E-25,
        +0.4045401335683533392143404142428E-25,
        +0.7067100658094689465651607717806E-26,
        -0.1249463344565105223002864518605E-25,
        +0.2867392244403437032979483391426E-26,
        +0.2044292892504292670281779574210E-26,
        -0.1518636633820462568371346802911E-26,
        +0.8110181098187575886132279107037E-28,
        +0.3580379354773586091127173703270E-27,
        -0.1692929018927902509593057175448E-27,
        -0.2222902499702427639067758527774E-28,
        +0.5424535127145969655048600401128E-28,
        -0.1787068401578018688764912993304E-28,
        -0.6565479068722814938823929437880E-29,
        +0.7807013165061145280922067706839E-29,
        -0.1816595260668979717379333152221E-29,
        -0.1287704952660084820376875598959E-29,
        +0.1114548172988164547413709273694E-29,
        -0.1808343145039336939159368876687E-30,
        -0.2231677718203771952232448228939E-30,
        +0.1619029596080341510617909803614E-30,
        -0.1834079908804941413901308439210E-31
    };
    static ityp ai1cs[46] =
    {
        -0.2846744181881478674100372468307E-01,
        -0.1922953231443220651044448774979E-01,
        -0.6115185857943788982256249917785E-03,
        -0.2069971253350227708882823777979E-04,
        +0.8585619145810725565536944673138E-05,
        +0.1049498246711590862517453997860E-05,
        -0.2918338918447902202093432326697E-06,
        -0.1559378146631739000160680969077E-07,
        +0.1318012367144944705525302873909E-07,
        -0.1448423418183078317639134467815E-08,
        -0.2908512243993142094825040993010E-09,
        +0.1266388917875382387311159690403E-09,
        -0.1664947772919220670624178398580E-10,
        -0.1666653644609432976095937154999E-11,
        +0.1242602414290768265232168472017E-11,
        -0.2731549379672432397251461428633E-12,
        +0.2023947881645803780700262688981E-13,
        +0.7307950018116883636198698126123E-14,
        -0.3332905634404674943813778617133E-14,
        +0.7175346558512953743542254665670E-15,
        -0.6982530324796256355850629223656E-16,
        -0.1299944201562760760060446080587E-16,
        +0.8120942864242798892054678342860E-17,
        -0.2194016207410736898156266643783E-17,
        +0.3630516170029654848279860932334E-18,
        -0.1695139772439104166306866790399E-19,
        -0.1288184829897907807116882538222E-19,
        +0.5694428604967052780109991073109E-20,
        -0.1459597009090480056545509900287E-20,
        +0.2514546010675717314084691334485E-21,
        -0.1844758883139124818160400029013E-22,
        -0.6339760596227948641928609791999E-23,
        +0.3461441102031011111108146626560E-23,
        -0.1017062335371393547596541023573E-23,
        +0.2149877147090431445962500778666E-24,
        -0.3045252425238676401746206173866E-25,
        +0.5238082144721285982177634986666E-27,
        +0.1443583107089382446416789503999E-26,
        -0.6121302074890042733200670719999E-27,
        +0.1700011117467818418349189802666E-27,
        -0.3596589107984244158535215786666E-28,
        +0.5448178578948418576650513066666E-29,
        -0.2731831789689084989162564266666E-30,
        -0.1858905021708600715771903999999E-30,
        +0.9212682974513933441127765333333E-31,
        -0.2813835155653561106370833066666E-31
    };
    static ityp bi1cs[17] =
    {
        -0.19717132610998597316138503218149E-02,
        +0.40734887667546480608155393652014,
        +0.34838994299959455866245037783787E-01,
        +0.15453945563001236038598401058489E-02,
        +0.41888521098377784129458832004120E-04,
        +0.76490267648362114741959703966069E-06,
        +0.10042493924741178689179808037238E-07,
        +0.99322077919238106481371298054863E-10,
        +0.76638017918447637275200171681349E-12,
        +0.47414189238167394980388091948160E-14,
        +0.24041144040745181799863172032000E-16,
        +0.10171505007093713649121100799999E-18,
        +0.36450935657866949458491733333333E-21,
        +0.11205749502562039344810666666666E-23,
        +0.29875441934468088832000000000000E-26,
        +0.69732310939194709333333333333333E-29,
        +0.14367948220620800000000000000000E-31
    };
    ityp eta;
    static dim_typ ntai1 = 0;
    static dim_typ ntai12 = 0;
    static dim_typ nti1 = 0;
    ityp value;
    static ityp xmin = 0.0;
    static ityp xsml = 0.0;
    ityp y;

    if(!nti1)
    {
        eta = 0.10 * r8_mach ( 3 );
        nti1 = r8_inits ( bi1cs, 17, eta );
        ntai1 = r8_inits ( ai1cs, 46, eta );
        ntai12 = r8_inits ( ai12cs, 69, eta );
        xmin = 2.00 * r8_mach ( 1 );
        xsml = sqrt ( 8.00 * r8_mach ( 3 ) );
    }

    y = abs ( x );

    if ( y <= xmin )
        value = 0.00;
    else if ( y <= xsml )
        value = 0.50 * x * exp ( - y );
    else if ( y <= 3.00 )
        value = x * ( 0.875 + r8_csevl ( y * y / 4.50 - 1.00, bi1cs, nti1 ) )* exp ( - y );
    else if ( y <= 8.00 )
    {
        value = ( 0.375 + r8_csevl ( ( 48.00 / y - 11.00) / 5.00, ai1cs, ntai1 ) ) / sqrt ( y );
        if ( x < 0.00 )
            value *= -1;
    }
    else
    {
        value = ( 0.375 + r8_csevl ( 16.00 / y - 1.00, ai12cs, ntai12 ) ) / sqrt ( y );
        if ( x < 0.00 )
            value *= -1;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_besj0 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BESJ0 evaluates the Bessel function J of order 0 of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BESJ0, the Bessel function J of order 0 of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    ityp ampl;
    static ityp bj0cs[19] =
    {
        +0.10025416196893913701073127264074,
        -0.66522300776440513177678757831124,
        +0.24898370349828131370460468726680,
        -0.33252723170035769653884341503854E-01,
        +0.23114179304694015462904924117729E-02,
        -0.99112774199508092339048519336549E-04,
        +0.28916708643998808884733903747078E-05,
        -0.61210858663032635057818407481516E-07,
        +0.98386507938567841324768748636415E-09,
        -0.12423551597301765145515897006836E-10,
        +0.12654336302559045797915827210363E-12,
        -0.10619456495287244546914817512959E-14,
        +0.74706210758024567437098915584000E-17,
        -0.44697032274412780547627007999999E-19,
        +0.23024281584337436200523093333333E-21,
        -0.10319144794166698148522666666666E-23,
        +0.40608178274873322700800000000000E-26,
        -0.14143836005240913919999999999999E-28,
        +0.43910905496698880000000000000000E-31
    };
    static dim_typ ntj0 = 0;
    ityp theta;
    ityp value;
    static ityp xsml = 0.0;
    ityp y;

    if (!ntj0)
    {
        ntj0 = r8_inits ( bj0cs, 19, 0.10 * r8_mach ( 3 ) );
        xsml = sqrt ( 4.00 * r8_mach ( 3 ) );
    }

    y = abs ( x );

    if ( y <= xsml )
        value = 1.00;
    else if ( y <= 4.00 )
        value = r8_csevl ( 0.125 * y * y - 1.00, bj0cs, ntj0 );
    else
    {
        r8_b0mp ( y, &ampl, &theta );
        value = ampl * cos ( theta );
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_besk ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BESK evaluates the Bessel function K of order NU of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 November 2012
  Author:
    John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double NU, the order.
    Input, double X, the argument.
    Output, double R8_BESK, the Bessel function K of order NU at X.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp nu = a_data[0];
	const register ityp x = a_data[1];
	
    ityp *bke;
    dim_typ nin;
    ityp value;
    ityp xnu;
    xnu = nu - ( dim_typ ) ( nu );
    nin = ( dim_typ ) ( nu ) + 1;
    bke = r8_besks ( xnu, x, nin );
    value = bke[nin-1];
    free ( bke );
    
    result = value;
    return &result;
}

#define R8BESK1_INVALIDRETURNVALUE MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_besk1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BESK1 evaluates the Bessel function K of order 1 of an R8 argument.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BESK1, the Bessel function K of order 1 of X.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp bk1cs[16] =
    {
        +0.25300227338947770532531120868533E-01,
        -0.35315596077654487566723831691801,
        -0.12261118082265714823479067930042,
        -0.69757238596398643501812920296083E-02,
        -0.17302889575130520630176507368979E-03,
        -0.24334061415659682349600735030164E-05,
        -0.22133876307347258558315252545126E-07,
        -0.14114883926335277610958330212608E-09,
        -0.66669016941993290060853751264373E-12,
        -0.24274498505193659339263196864853E-14,
        -0.70238634793862875971783797120000E-17,
        -0.16543275155100994675491029333333E-19,
        -0.32338347459944491991893333333333E-22,
        -0.53312750529265274999466666666666E-25,
        -0.75130407162157226666666666666666E-28,
        -0.91550857176541866666666666666666E-31
    };
    static dim_typ ntk1 = 0;
    ityp value;
    static ityp xmax = 0.00;
    static ityp xmin = 0.00;
    static ityp xsml = 0.00;
    ityp y;

    if (!ntk1 )
    {
        ntk1 = r8_inits ( bk1cs, 16, 0.1 * r8_mach ( 3 ) );
        xmin = exp ( MAX ( log ( r8_mach ( 1 ) ), - log ( r8_mach ( 2 ) ) ) + 0.01 );
        xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
        xmax = - log ( r8_mach ( 1 ) );
        xmax -= 0.5 * xmax * log ( xmax ) / ( xmax + 0.50 ) - 0.01;
    }

    if ( x <= 0.00 )
    {
    	result = R8BESK1_INVALIDRETURNVALUE;
        return &result;
    }
    else if ( x <= xsml )
    {
        y = 0.00;
        value = log ( 0.50 * x ) * r8_besi1 ( x ) + ( 0.75 + r8_csevl ( 0.50 * y - 1.00, bk1cs, ntk1 ) ) / x;
    }
    else if ( x <= 2.00 )
    {
        y = x * x;
        value = log ( 0.50 * x ) * r8_besi1 ( x ) + ( 0.75 + r8_csevl ( 0.50 * y - 1.00, bk1cs, ntk1 ) ) / x;
    }
    else if ( x <= xmax )
        value = exp ( - x ) * r8_besk1e ( x );
    else
        value = 0.00;

	result = value;
    return &result;
}

#define BESK1E_INVALIDRETURNVALUE MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_besk1e ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BESK1E evaluates the exponentially scaled Bessel function K1(X).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_BESK1E, the exponentially scaled Bessel
    function K1(X).
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp ak12cs[33] =
    {
        +0.6379308343739001036600488534102E-01,
        +0.2832887813049720935835030284708E-01,
        -0.2475370673905250345414545566732E-03,
        +0.5771972451607248820470976625763E-05,
        -0.2068939219536548302745533196552E-06,
        +0.9739983441381804180309213097887E-08,
        -0.5585336140380624984688895511129E-09,
        +0.3732996634046185240221212854731E-10,
        -0.2825051961023225445135065754928E-11,
        +0.2372019002484144173643496955486E-12,
        -0.2176677387991753979268301667938E-13,
        +0.2157914161616032453939562689706E-14,
        -0.2290196930718269275991551338154E-15,
        +0.2582885729823274961919939565226E-16,
        -0.3076752641268463187621098173440E-17,
        +0.3851487721280491597094896844799E-18,
        -0.5044794897641528977117282508800E-19,
        +0.6888673850418544237018292223999E-20,
        -0.9775041541950118303002132480000E-21,
        +0.1437416218523836461001659733333E-21,
        -0.2185059497344347373499733333333E-22,
        +0.3426245621809220631645388800000E-23,
        -0.5531064394246408232501248000000E-24,
        +0.9176601505685995403782826666666E-25,
        -0.1562287203618024911448746666666E-25,
        +0.2725419375484333132349439999999E-26,
        -0.4865674910074827992378026666666E-27,
        +0.8879388552723502587357866666666E-28,
        -0.1654585918039257548936533333333E-28,
        +0.3145111321357848674303999999999E-29,
        -0.6092998312193127612416000000000E-30,
        +0.1202021939369815834623999999999E-30,
        -0.2412930801459408841386666666666E-31
    };
    static ityp ak1cs[38] =
    {
        +0.27443134069738829695257666227266,
        +0.75719899531993678170892378149290E-01,
        -0.14410515564754061229853116175625E-02,
        +0.66501169551257479394251385477036E-04,
        -0.43699847095201407660580845089167E-05,
        +0.35402774997630526799417139008534E-06,
        -0.33111637792932920208982688245704E-07,
        +0.34459775819010534532311499770992E-08,
        -0.38989323474754271048981937492758E-09,
        +0.47208197504658356400947449339005E-10,
        -0.60478356628753562345373591562890E-11,
        +0.81284948748658747888193837985663E-12,
        -0.11386945747147891428923915951042E-12,
        +0.16540358408462282325972948205090E-13,
        -0.24809025677068848221516010440533E-14,
        +0.38292378907024096948429227299157E-15,
        -0.60647341040012418187768210377386E-16,
        +0.98324256232648616038194004650666E-17,
        -0.16284168738284380035666620115626E-17,
        +0.27501536496752623718284120337066E-18,
        -0.47289666463953250924281069568000E-19,
        +0.82681500028109932722392050346666E-20,
        -0.14681405136624956337193964885333E-20,
        +0.26447639269208245978085894826666E-21,
        -0.48290157564856387897969868800000E-22,
        +0.89293020743610130180656332799999E-23,
        -0.16708397168972517176997751466666E-23,
        +0.31616456034040694931368618666666E-24,
        -0.60462055312274989106506410666666E-25,
        +0.11678798942042732700718421333333E-25,
        -0.22773741582653996232867840000000E-26,
        +0.44811097300773675795305813333333E-27,
        -0.88932884769020194062336000000000E-28,
        +0.17794680018850275131392000000000E-28,
        -0.35884555967329095821994666666666E-29,
        +0.72906290492694257991679999999999E-30,
        -0.14918449845546227073024000000000E-30,
        +0.30736573872934276300799999999999E-31
    };
    static ityp bk1cs[16] =
    {
        +0.25300227338947770532531120868533E-01,
        -0.35315596077654487566723831691801,
        -0.12261118082265714823479067930042,
        -0.69757238596398643501812920296083E-02,
        -0.17302889575130520630176507368979E-03,
        -0.24334061415659682349600735030164E-05,
        -0.22133876307347258558315252545126E-07,
        -0.14114883926335277610958330212608E-09,
        -0.66669016941993290060853751264373E-12,
        -0.24274498505193659339263196864853E-14,
        -0.70238634793862875971783797120000E-17,
        -0.16543275155100994675491029333333E-19,
        -0.32338347459944491991893333333333E-22,
        -0.53312750529265274999466666666666E-25,
        -0.75130407162157226666666666666666E-28,
        -0.91550857176541866666666666666666E-31
    };
    ityp eta;
    static dim_typ ntak1 = 0;
    static dim_typ ntak12 = 0;
    static dim_typ ntk1 = 0;
    ityp value;
    static ityp xmin = 0.0;
    static ityp xsml = 0.0;
    ityp y;

    if (!ntk1 )
    {
        eta = 0.10 * r8_mach ( 3 );
        ntk1 = r8_inits ( bk1cs, 16, eta );
        ntak1 = r8_inits ( ak1cs, 38, eta );
        ntak12 = r8_inits ( ak12cs, 33, eta );
        xmin = exp ( MAX ( log ( r8_mach ( 1 ) ), - log ( r8_mach ( 2 ) ) ) + 0.01 );
        xsml = sqrt ( 4.00 * r8_mach ( 3 ) );
    }

    if ( x <= 0.00 )
    {
    	result = BESK1E_INVALIDRETURNVALUE;
        return &result;
    }
    else if ( x <= xsml )
    {
        y = 0.00;
        value = exp ( x ) * ( log ( 0.50 * x ) * r8_besi1 ( x )+ ( 0.75 + r8_csevl ( 0.50 * y - 1.00, bk1cs, ntk1 ) ) / x );
    }
    else if ( x <= 2.00 )
    {
        y = x * x;
        value = exp ( x ) * ( log ( 0.50 * x ) * r8_besi1 ( x )+ ( 0.75 + r8_csevl ( 0.50 * y - 1.00, bk1cs, ntk1 ) ) / x );
    }
    else if ( x <= 8.00 )
        value = ( 1.25 + r8_csevl ( ( 16.00 / x - 5.00 ) / 3.00, ak1cs, ntak1 ) ) / sqrt ( x );
    else
    {
        value = ( 1.25 +
        r8_csevl ( 16.00 / x - 1.00, ak12cs, ntak12 ) ) / sqrt ( x );
    }
    
    result = value;
    return &result;
}

#define BESK1E_INVALIDRETURNVALUE MAX_VAL
#define BESKES_INVALIDRETURNVALUE NULL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8_beskes ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BESKES: a sequence of exponentially scaled K Bessel functions at X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 November 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double XNU, ?
    |XNU| < 1.
    Input, double X, the argument.
    Input, int NIN, indicates the number of terms to compute.
    Output, double R8_BESKES(abs(NIN)), the exponentially scaled
    K Bessel functions.
*/
{
	ityp * const a_data = data;
	const register ityp xnu = a_data[0];
	const register ityp x = a_data[1];
	const register ityp nin = a_data[2];
	
    ityp *bke = ( ityp * ) malloc ( abs ( nin ) * sizeof ( ityp ) );
    ityp bknu1;
    ityp direct;
    dim_typ i;
    dim_typ iswtch, n;
    ityp v;
    ityp vend;
    ityp vincr;

    v = abs ( xnu );
    n = i4_abs ( nin );

    if ( 1.00 <= v || x <= 0.00 || !n)
        return BESKES_INVALIDRETURNVALUE;

    r8_knus ( v, x, &bke[0], &bknu1, &iswtch );

    if ( n == 1 )
        return bke;

    vincr = +1.00 - ((nin<0)<<1);
    direct = vincr*(1-((xnu<0.00)<<1));
    bke[1] = bknu1;

    if ( direct < 0.00 )
        r8_knus ( abs ( xnu + vincr ), x, &bke[1], &bknu1, &iswtch );

    if ( n == 2 )
        return bke;

    vend = abs ( xnu + ( ityp ) ( nin ) ) - 1.00;

    v = xnu;
    for ( i = 3; i <= n; ++i )
    {
        v += vincr;
        bke[i-1] = 2.00 * v * bke[i-2] / x + bke[i-3];
    }
    return bke;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8_besks ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_BESKS evaluates a sequence of K Bessel functions at X.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 November 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double XNU, ?
    |XNU| < 1.
    Input, double X, the argument.
    Input, int NIN, indicates the number of terms to compute.
    Output, double R8_BESKS(abs(NIN)), the K Bessel functions.
*/
{
	ityp * const a_data = data;
	const register ityp xnu = a_data[0];
	const register ityp x = a_data[1];
	const register ityp nin = a_data[2];
	
    ityp *bk;
    ityp expxi;
    dim_typ i, n;
    static ityp xmax = 0.00;

    if ( xmax == 0.00 )
    {
        xmax = - log ( r8_mach ( 1 ) );
        xmax += 0.50 * log ( 3.14 * 0.50 / xmax );
    }

    bk = r8_beskes ( xnu, x, nin );
    expxi = exp ( - x );
    n = i4_abs ( nin );

    for ( i = 0; i < n; ++i)
        bk[i] = expxi * bk[i];
    return bk;
}

#define R8CSEVL_INVALIDRETURNVALUE MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_csevl ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_CSEVL evaluates a Chebyshev series.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    C version by John Burkardt.
  Reference:
    Roger Broucke,
    Algorithm 446:
    Ten Subroutines for the Manipulation of Chebyshev Series,
    Communications of the ACM,
    Volume 16, Number 4, April 1973, pages 254-256.
  Parameters:
    Input, double X, the evaluation point.
    Input, double CS[N], the Chebyshev coefficients.
    Input, int N, the number of Chebyshev coefficients.
    Output, double R8_CSEVL, the Chebyshev series evaluated at X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpitit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * a = s_data->a1;
	const register ityp x = s_data->a2;
	
    ityp b0;
    ityp b1;
    ityp b2;
    int i;
    ityp twox;
    ityp value;

    if ( n < 1 || 1000<n)
    {
    	result = R8CSEVL_INVALIDRETURNVALUE;
        return &result;
    }

    twox = 2.00 * x;
    b1 = b0 = 0.00;

    for (dim_typ i = n - 1; 0 <= i; --i )
    {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[i];
    }

	result = 0.50 * ( b0 - b2 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _r8_gaml ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_GAML evaluates bounds for an R8 argument of the gamma function.
  Discussion:
    This function calculates the minimum and maximum legal bounds
    for X in the evaluation of GAMMA ( X ).
    XMIN and XMAX are not the only bounds, but they are the only
    non-trivial ones to calculate.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Output, double *XMIN, *XMAX, the bounds.
*/
{
	ityp ** const a_data = data;
	ityp * xmin = a_data[0];
	ityp * xmax = a_data[1];
	
    ityp alnbig;
    ityp alnsml;
    dim_typ i, j;
    ityp xln;
    ityp xold;

    alnsml = log ( r8_mach ( 1 ) );
    *xmin = - alnsml;

    for ( i = 1; i <= 10; ++i )
    {
        xold = *xmin;
        xln = log ( *xmin );
        *xmin -=*xmin * ( ( *xmin + 0.50 ) * xln - *xmin - 0.2258 + alnsml ) / ( *xmin * xln + 0.50 );

        if ( abs ( *xmin - xold ) < 0.005 )
        {
            *xmin = - *xmin + 0.01;

            alnbig = log ( r8_mach ( 2 ) );
            *xmax = alnbig;

            for ( j = 1; j <= 10; ++j )
            {
                xold = *xmax;
                xln = log ( *xmax );
                *xmax -= *xmax * ( ( *xmax - 0.5 ) * xln - *xmax + 0.9189 - alnbig ) / ( *xmax * xln - 0.5 );

                if ( abs ( *xmax - xold ) < 0.005 )
                {
                    *xmax -= 0.01;
                    *xmin = MAX ( *xmin, - *xmax + 1.00 );
                    return NULL;
                }
            }
        }
    }
    return NULL;
}

#define R8INITS_INVALIDRETURNVALUE INT_MAX

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_inits ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_INITS initializes a Chebyshev series.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    C version by John Burkardt.
  Reference:
    Roger Broucke,
    Algorithm 446:
    Ten Subroutines for the Manipulation of Chebyshev Series,
    Communications of the ACM,
    Volume 16, Number 4, April 1973, pages 254-256.
  Parameters:
    Input, double DOS[NOS], the Chebyshev coefficients.
    Input, int NOS, the number of coefficients.
    Input, double ETA, the desired accuracy.
    Output, int R8_INITS, the number of terms of the series needed
    to ensure the requested accuracy.
*/
{
	static int result = INT_MAX;
	
	const dtpitit * const s_data = data;
	
	const register dim_typ nos = s_data->a0;
	ityp * dos = s_data->a1;
	const register ityp eta = s_data->a2;
	
    ityp err;
    int i;
    int value;

    if ( nos < 1 )
    {
    	result = R8INITS_INVALIDRETURNVALUE;
        return &result;
    }
    err = 0.00;

    for ( i = nos - 1; 0 <= i; --i )
    {
        err += abs ( dos[i] );
        if ( eta < err )
        {
            value = i + 1;
            result = value;
            return &result;
        }
    }
    
    result = i;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_knus ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_KNUS computes a sequence of K Bessel functions.
  Discussion:
    This routine computes Bessel functions
      exp(x) * k-sub-xnu (x)
    and
      exp(x) * k-sub-xnu+1 (x)
    for 0.0 <= xnu < 1.0.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double XNU, the order parameter.
    Input, double X, the argument.
    Output, double *BKNU, *BKNU1, the two K Bessel functions.
    Output, int *ISWTCH, ?
*/
{
	const _2it2pitpdt * const s_data = data;
	const register ityp xnu = s_data->a0;
	const register ityp x = s_data->a1;
	ityp * bknu = s_data->a2;
	ityp * bknu1 = s_data->a3;
	dim_typ * iswtch = s_data->a4;
	
    ityp a[32];
    ityp a0;
    static ityp aln2 = 0.69314718055994530941723212145818;
    static ityp alnbig = 0;
    static ityp alneps = 0;
    static ityp alnsml = 0;
    ityp alnz;
    ityp alpha[32];
    ityp an;
    ityp b0;
    ityp beta[32];
    ityp bknu0;
    ityp bknud;
    ityp bn;
    ityp c0;
    static ityp c0kcs[29] =
    {
        +0.60183057242626108387577445180329E-01,
        -0.15364871433017286092959755943124,
        -0.11751176008210492040068229226213E-01,
        -0.85248788891979509827048401550987E-03,
        -0.61329838767496791874098176922111E-04,
        -0.44052281245510444562679889548505E-05,
        -0.31631246728384488192915445892199E-06,
        -0.22710719382899588330673771793396E-07,
        -0.16305644608077609552274620515360E-08,
        -0.11706939299414776568756044043130E-09,
        -0.84052063786464437174546593413792E-11,
        -0.60346670118979991487096050737198E-12,
        -0.43326960335681371952045997366903E-13,
        -0.31107358030203546214634697772237E-14,
        -0.22334078226736982254486133409840E-15,
        -0.16035146716864226300635791528610E-16,
        -0.11512717363666556196035697705305E-17,
        -0.82657591746836959105169479089258E-19,
        -0.59345480806383948172333436695984E-20,
        -0.42608138196467143926499613023976E-21,
        -0.30591266864812876299263698370542E-22,
        -0.21963541426734575224975501815516E-23,
        -0.15769113261495836071105750684760E-24,
        -0.11321713935950320948757731048056E-25,
        -0.81286248834598404082792349714433E-27,
        -0.58360900893453226552829349315949E-28,
        -0.41901241623610922519452337780905E-29,
        -0.30083737960206435069530504212862E-30,
        -0.21599152067808647728342168089832E-31
    };
    ityp eta;
    static ityp euler = 0.57721566490153286060651209008240;
    ityp expx;
    dim_typ i;
    dim_typ ii;
    dim_typ inu;
    dim_typ n;
    static dim_typ ntc0k = 0;
    dim_typ nterms;
    static dim_typ ntznu1 = 0;
    ityp p1;
    ityp p2;
    ityp p3;
    ityp qq;
    ityp result;
    ityp sqrtx;
    ityp v;
    ityp vlnz;
    ityp x2n;
    ityp x2tov;
    ityp xi;
    ityp xmu;
    static ityp xnusml = 0.00;
    static ityp xsml = 0.00;
    ityp z;
    static ityp znu1cs[20] =
    {
        +0.203306756994191729674444001216911,
        +0.140077933413219771062943670790563,
        +0.791679696100161352840972241972320E-02,
        +0.339801182532104045352930092205750E-03,
        +0.117419756889893366664507228352690E-04,
        +0.339357570612261680333825865475121E-06,
        +0.842594176976219910194629891264803E-08,
        +0.183336677024850089184748150900090E-09,
        +0.354969844704416310863007064469557E-11,
        +0.619032496469887332205244342078407E-13,
        +0.981964535680439424960346115456527E-15,
        +0.142851314396490474211473563005985E-16,
        +0.191894921887825298966162467488436E-18,
        +0.239430979739498914162313140597128E-20,
        +0.278890246815347354835870465474995E-22,
        +0.304606650633033442582845214092865E-24,
        +0.313173237042191815771564260932089E-26,
        +0.304133098987854951645174908005034E-28,
        +0.279840384636833084343185097659733E-30,
        +0.244637186274497596485238794922666E-32
    };
    ityp ztov;

    if (!ntc0k )
    {
        eta = 0.10 * r8_mach ( 3 );
        ntc0k = r8_inits ( c0kcs, 29, eta );
        ntznu1 = r8_inits ( znu1cs, 20, eta );
        xnusml = sqrt ( r8_mach ( 3 ) / 8.00 );
        xsml = 0.10 * r8_mach ( 3 );
        alnsml = log ( r8_mach ( 1 ) );
        alnbig = log ( r8_mach ( 2 ) );
        alneps = log ( 0.10 * r8_mach ( 3 ) );
    }

    if ( xnu < 0.00 || 1.00 <= xnu || x <= 0.00)
        return NULL;

    *iswtch = 0;
    /*
    X is small.  Compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
    then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-0.5,+0.5)
    then to (0., .5), because k of negative order (-nu) = k of positive
    order (+nu).
    */
    if ( x <= 2.00 )
    {
        if ( xnu <= 0.50 )
            v = xnu;
        else
            v = 1.00 - xnu;
        /*
        Carefully find (x/2)^xnu and z^xnu where z = x*x/4.
        */
        alnz = 2.00 * ( log ( x ) - aln2 );
        if ( x <= xnu && alnbig < - 0.5 * xnu * alnz - aln2 - log ( xnu ) )
            return NULL;

        vlnz = v * alnz;
        x2tov = exp ( 0.50 * vlnz );

        ztov = vlnz <= alnsml ? 0.00: x2tov*x2tov;

        a0 = 0.50 * r8_gamma ( 1.00 + v );
        b0 = 0.50 * r8_gamma ( 1.00 - v );
        c0 = - euler;

        if ( 0.50 <= ztov && xnusml < v )
            c0 = - 0.75 + r8_csevl ( ( 8.00 * v ) * v - 1.00, c0kcs, ntc0k );

        if ( ztov <= 0.50 )
            alpha[0] = ( a0 - ztov * b0 ) / v;
        else
        {
            alpha[0] = c0 - alnz * ( 0.75 +
            r8_csevl ( vlnz / 0.35 + 1.00, znu1cs, ntznu1 ) ) * b0;
        }

        beta[0] = - 0.50 * ( a0 + ztov * b0 );

        z = x<=xsml ? 0.00 : 0.25*x*x;
        nterms = MAX ( 2, ( dim_typ ) ( 11.00 + ( 8.00 * alnz - 25.19 - alneps ) / ( 4.28 - alnz ) ) );

        for ( i = 2; i <= nterms; ++i )
        {
            xi = ( ityp ) ( i - 1 );
            a0 /= ( xi * ( xi - v ) );
            b0 /= ( xi * ( xi + v ) );
            alpha[i-1] = ( alpha[i-2] + 2.00 * xi * a0 ) / ( xi * ( xi + v ) );
            beta[i-1] = ( xi - 0.50 * v ) * alpha[i-1] - ztov * b0;
        }

        *bknu = alpha[nterms-1];
        bknud = beta[nterms-1];
        for ( ii = 2; ii <= nterms; ++ii)
        {
            i = nterms + 1 - ii;
            *bknu = alpha[i-1] + *bknu * z;
            bknud = beta[i-1] + bknud * z;
        }

        expx = exp ( x );
        *bknu *= expx/x2tov;

        if ( alnbig < - 0.50 * ( xnu + 1.00 ) * alnz - 2.00 * aln2 )
        {
            *iswtch = 1;
            return NULL;
        }

        bknud *= expx * 2.00 / ( x2tov * x );

        if ( xnu <= 0.50 )
        {
            *bknu1 *= v / x - bknud;
            return NULL;
        }
        bknu0 = *bknu;
        *bknu = - v * *bknu / x - bknud;
        *bknu1 *= 2.00 * xnu / x + bknu0;
    }
    /*
    X is large.  find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke-s
    rational expansion.
    */
    else
    {
        sqrtx = sqrt ( x );

        if ( 1.00 / xsml < x )
        {
            *bknu = M_SQRTPI2 / sqrtx;
            *bknu1 = *bknu;
            return NULL;
        }

        an = - 0.60 - 1.02 / x;
        bn = - 0.27 - 0.53 / x;
        nterms = MIN ( 32, MAX ( 3, ( dim_typ ) ( an + bn * alneps ) ) );

        #pragma omp parallel for num_threads(2)
        for ( inu = 1; inu <= 2; ++inu )
            xmu = inu == 1 ? xnu<=xnusml ? 0.00 : (4.00*xnu)*xnu : 4.00 * ( abs ( xnu ) + 1.00 ) * ( abs ( xnu ) + 1.00 );

        a[0] = 1.00 - xmu;
        a[1] = 9.00 - xmu;
        a[2] = 25.00 - xmu;

        if ( a[1] == 0.00 )
        result = M_SQRTPI2 * ( 16.00 * x + xmu + 7.00 ) / ( 16.00 * x * sqrtx );
        else
        {
            alpha[0] = 1.00;
            alpha[1] = ( 16.00 * x + a[1] ) / a[1];
            alpha[2] = ( ( 768.00 * x + 48.00 * a[2] ) * x + a[1] * a[2] ) / ( a[1] * a[2] );

            beta[0] = 1.0;
            beta[1] = ( 16.0 * x + ( xmu + 7.0 ) ) / a[1];
            beta[2] = ( ( 768.0 * x + 48.0 * ( xmu + 23.0 ) ) * x +( ( xmu + 62.0 ) * xmu + 129.0 ) ) / ( a[1] * a[2] );

	
            for ( i = 4; i <= nterms; ++i )
            {
                n = i - 1;
                x2n = ( ityp ) ( (n<<1) - 1 );
                a[i-1] = ( x2n + 2.00 ) * ( x2n + 2.00 ) - xmu;
                qq = 16.00 * x2n / a[i-1];
                p1 = - x2n * ( ( ityp ) ( 12 * n * n - 20 * n ) - a[0] ) / ( ( x2n - 2.00 ) * a[i-1] ) - qq * x;
                p2 = ( ( ityp ) ( 12 * n * n - 28 * n + 8 ) - a[0] ) / a[i-1] - qq * x;
                p3 = - x2n * a[i-4] / ( ( x2n - 2.00 ) * a[i-1] );
                alpha[i-1] = - p1 * alpha[i-2]- p2 * alpha[i-3] - p3 * alpha[i-4];
                beta[i-1] =  - p1 * beta[i-2]- p2 * beta[i-3] - p3 * beta[i-4];

            }
            result = M_SQRTPI2 * beta[nterms-1] / ( sqrtx * alpha[nterms-1] );
        }

        if ( inu == 1 )
            *bknu = result;
        else
            *bknu1 = result;
    }
    return NULL;
}

#define R8LGMC_INVALIDRETURNVALUE MAX_VAL

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_lgmc ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8_LGMC evaluates the log gamma correction factor for an R8 argument.
  Discussion:
    For 10 <= X, compute the log gamma correction factor so that
      log ( gamma ( x ) ) = log ( sqrt ( 2 * M_PI ) )
                          + ( x - 0.5 ) * log ( x ) - x
                          + r8_lgmc ( x )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 January 2012
  Author:
    Original FORTRAN77 version by Wayne Fullerton.
    C version by John Burkardt.
  Reference:
    Wayne Fullerton,
    Portable Special Function Routines,
    in Portability of Numerical Software,
    edited by Wayne Cowell,
    Lecture Notes in Computer Science, Volume 57,
    Springer 1977,
    ISBN: 978-3-540-08446-4,
    LC: QA297.W65.
  Parameters:
    Input, double X, the argument.
    Output, double R8_LGMC, the correction factor.
*/
{
	static ityp result = MAX_VAL;
	
	const register ityp x = *(ityp *) data;
	
    static ityp algmcs[15] =
    {
        +0.1666389480451863247205729650822,
        -0.1384948176067563840732986059135E-04,
        +0.9810825646924729426157171547487E-08,
        -0.1809129475572494194263306266719E-10,
        +0.6221098041892605227126015543416E-13,
        -0.3399615005417721944303330599666E-15,
        +0.2683181998482698748957538846666E-17,
        -0.2868042435334643284144622399999E-19,
        +0.3962837061046434803679306666666E-21,
        -0.6831888753985766870111999999999E-23,
        +0.1429227355942498147573333333333E-24,
        -0.3547598158101070547199999999999E-26,
        +0.1025680058010470912000000000000E-27,
        -0.3401102254316748799999999999999E-29,
        +0.1276642195630062933333333333333E-30
    };
    static dim_typ nalgm = 0;
    ityp value;
    static ityp xbig = 0.00;
    static ityp xmax = 0.00;

    if ( nalgm == 0 )
    {
        nalgm = r8_inits ( algmcs, 15, r8_mach ( 3 ) );
        xbig = 1.00 / sqrt ( r8_mach ( 3 ) );
        xmax = exp ( MIN ( log ( r8_mach ( 2 ) / 12.00 ), - log ( 12.00 * r8_mach ( 1 ) ) ) );
    }

    if ( x < 10.00 )
    {
    	result = R8LGMC_INVALIDRETURNVALUE;
        return &result;
    }
    else if ( x < xbig )
        value = r8_csevl ( 2.00 * ( 10.00 / x ) * ( 10.00 / x ) - 1.00, algmcs, nalgm ) / x;
    else if ( x < xmax )
        value = 1.00 / ( 12.00 * x );
    else
        value = 0.00;

	result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_mach ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8_MACH returns single precision real machine constants.
  Discussion:
    Assume that single precision real numbers are stored with a mantissa
    of T digits in base B, with an exponent whose value must lie
    between EMIN and EMAX.  Then for values of I between 1 and 5,
    r8_MACH will return the following values:
      r8_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
      r8_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
      r8_MACH(3) = B^(-T), the smallest relative spacing.
      r8_MACH(4) = B^(1-T), the largest relative spacing.
      r8_MACH(5) = log10(B)
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
    1 <= I <= 5.
    Output, ityp r8_MACH, the value of the chosen parameter.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ i = *(dim_typ *) data; 
	
    ityp value;

    if ( i < 1 || 5 < i )
    {
    	result = MAX_VAL;
        return &result;
    }

    switch(i)
    {
        case 1:
            value = 1.1754944E-38;
            break;
        case 2:
            value = 3.4028235E+38;
            break;
        case 3:
            value = 5.9604645E-08;
            break;
        case 4:
            value = 1.1920929E-07;
            break;
        case 5:
            value = 0.3010300E+00;
            break;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8vec_linspace ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8VEC_LINSPACE creates a vector of linearly spaced values.
  Discussion:
    An R8VEC is a vector of R8's.
    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 March 2011
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, double A, B, the first and last entries.
    Output, double X[N], a vector of linearly spaced data.
*/
{
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	ityp * x = s_data->a3;
	
	
    if ( n == 1 )
        x[0] = ( a + b ) / 2.00;
    else
        for (dim_typ i = 0; i < n; ++i )
            x[i] = ( ( ityp ) ( n - 1 - i ) * a + ( ityp ) (         i ) * b ) / ( ityp ) ( n - 1     );
            
    return NULL;
}

#endif
