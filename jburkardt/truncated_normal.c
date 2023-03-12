#ifndef __DISABLEDEEP_TRUNCATEDNORMAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _truncated_normal_ab_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_AB_CDF_VALUES: values of the Truncated Normal CDF.
  Discussion:
    The Normal distribution, with mean Mu and standard deviation Sigma,
    is truncated to the interval [A,B].
  Note:
    a = 50.00
    b = 150.00
    mu = 100.00
    sigma = 25.00
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 September 2013
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
    Output, double *MU, the mean of the distribution.
    Output, double *SIGMA, the standard deviation of the distribution.
    Output, double *A, *B, the lower and upper truncation limits.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt6pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * sigma = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	ityp * x = s_data->a5;
	ityp * fx = s_data->a6;
	
    # define N_MAX 11

    static ityp fx_vec[N_MAX] =
    {
        0.3371694242213513,
        0.3685009225506048,
        0.4006444233448185,
        0.4334107066903040,
        0.4665988676496338,
        0.5000000000000000,
        0.5334011323503662,
        0.5665892933096960,
        0.5993555766551815,
        0.6314990774493952,
        0.6628305757786487
    };

    static ityp x_vec[N_MAX] =
    {
        90.00,
        92.00,
        94.00,
        96.00,
        98.00,
        100.00,
        102.00,
        104.00,
        106.00,
        108.00,
        110.00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _truncated_normal_ab_pdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_AB_PDF_VALUES: values of the Truncated Normal PDF.
  Discussion:
    The Normal distribution, with mean Mu and standard deviation Sigma,
    is truncated to the interval [A,B].
  Note:
    a = 50.00
    b = 150.00
    mu = 100.00
    sigma = 25.00
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 September 2013
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
    Output, double *MU, the mean of the distribution.
    Output, double *SIGMA, the standard deviation of the distribution.
    Output, double *A, *B, the lower and upper truncation limits.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt6pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * sigma = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	ityp * x = s_data->a5;
	ityp * fx = s_data->a6;
	
    # define N_MAX 11

    static ityp fx_vec[N_MAX] =
    {
        0.01543301171801836,
        0.01588394472270638,
        0.01624375997031919,
        0.01650575046469259,
        0.01666496869385951,
        0.01671838200940538,
        0.01666496869385951,
        0.01650575046469259,
        0.01624375997031919,
        0.01588394472270638,
        0.01543301171801836
    };

    static ityp x_vec[N_MAX] =
    {
        90.00,
        92.00,
        94.00,
        96.00,
        98.00,
        100.00,
        102.00,
        104.00,
        106.00,
        108.00,
        110.00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _truncated_normal_a_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_A_CDF_VALUES: values of the Lower Truncated Normal CDF.
  Discussion:
    The Normal distribution, with mean Mu and standard deviation Sigma,
    is truncated to the interval [A,+oo).
  Note:
    a = 50.00
    mu = 100.00
    sigma = 25.00
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 September 2013
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
    Output, double *MU, the mean of the distribution.
    Output, double *SIGMA, the standard deviation of the distribution.
    Output, double *A, the lower truncation limit.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt6pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * sigma = s_data->a2;
	ityp * a = s_data->a3;
	ityp * x = s_data->a5;
	ityp * fx = s_data->a6;
	
    # define N_MAX 11

    static ityp fx_vec[N_MAX] =
    {
        0.3293202045481688,
        0.3599223134505957,
        0.3913175216041539,
        0.4233210140873113,
        0.4557365629792204,
        0.4883601253415709,
        0.5209836877039214,
        0.5533992365958304,
        0.5854027290789878,
        0.6167979372325460,
        0.6474000461349729
    };

    static ityp x_vec[N_MAX] =
    {
        90.00,
        92.00,
        94.00,
        96.00,
        98.00,
        100.00,
        102.00,
        104.00,
        106.00,
        108.00,
        110.00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/*****************************************************************************/
__MATHSUITE __JBURKARDT  void *  _truncated_normal_a_pdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_A_PDF_VALUES: values of the Lower Truncated Normal PDF.
  Discussion:
    The Normal distribution, with mean Mu and standard deviation Sigma,
    is truncated to the interval [A,+oo).
  Note:
    a = 50.00
    mu = 100.00
    sigma = 25.00
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 September 2013
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
    Output, double *MU, the mean of the distribution.
    Output, double *SIGMA, the standard deviation of the distribution.
    Output, double *A, the lower truncation limit.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt6pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * sigma = s_data->a2;
	ityp * a = s_data->a3;
	ityp * x = s_data->a5;
	ityp * fx = s_data->a6;
	
    # define N_MAX 11

    static ityp fx_vec[N_MAX] =
    {
        0.01507373507401876,
        0.01551417047139894,
        0.01586560931024694,
        0.01612150073158793,
        0.01627701240029317,
        0.01632918226724295,
        0.01627701240029317,
        0.01612150073158793,
        0.01586560931024694,
        0.01551417047139894,
        0.01507373507401876
    };

    static ityp x_vec[N_MAX] =
    {
        90.00,
        92.00,
        94.00,
        96.00,
        98.00,
        100.00,
        102.00,
        104.00,
        106.00,
        108.00,
        110.00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _truncated_normal_b_cdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_B_CDF_VALUES: values of the Upper Truncated Normal CDF.
  Discussion:
    The Normal distribution, with mean Mu and standard deviation Sigma,
    is truncated to the interval [A,B].
  Note:
    b = 150.00
    mu = 100.00
    sigma = 25.00
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 September 2013
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
    Output, double *MU, the mean of the distribution.
    Output, double *SIGMA, the standard deviation of the distribution.
    Output, double *B, the upper truncation limit.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt6pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * sigma = s_data->a2;
	ityp * b = s_data->a4;
	ityp * x = s_data->a5;
	ityp * fx = s_data->a6;
	
    # define N_MAX 11

    static ityp fx_vec[N_MAX] =
    {
        0.3525999538650271,
        0.3832020627674540,
        0.4145972709210122,
        0.4466007634041696,
        0.4790163122960786,
        0.5116398746584291,
        0.5442634370207796,
        0.5766789859126887,
        0.6086824783958461,
        0.6400776865494043,
        0.6706797954518312
    };

    static ityp x_vec[N_MAX] =
    {
        90.00,
        92.00,
        94.00,
        96.00,
        98.00,
        100.00,
        102.00,
        104.00,
        106.00,
        108.00,
        110.00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _truncated_normal_b_pdf_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRUNCATED_NORMAL_B_PDF_VALUES: values of the Upper Truncated Normal PDF.
  Discussion:
    The Normal distribution, with mean Mu and standard deviation Sigma,
    is truncated to the interval (-oo,B].
  Note:
    b = 150.00
    mu = 100.00
    sigma = 25.00
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 September 2013
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
    Output, double *MU, the mean of the distribution.
    Output, double *SIGMA, the standard deviation of the distribution.
    Output, double *B, the upper truncation limit.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const pdt6pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * mu = s_data->a1;
	ityp * sigma = s_data->a2;
	ityp * b = s_data->a4;
	ityp * x = s_data->a5;
	ityp * fx = s_data->a6;
	
    # define N_MAX 11

    static ityp fx_vec[N_MAX] =
    {
        0.01507373507401876,
        0.01551417047139894,
        0.01586560931024694,
        0.01612150073158793,
        0.01627701240029317,
        0.01632918226724295,
        0.01627701240029317,
        0.01612150073158793,
        0.01586560931024694,
        0.01551417047139894,
        0.01507373507401876
    };

    static ityp x_vec[N_MAX] =
    {
        90.00,
        92.00,
        94.00,
        96.00,
        98.00,
        100.00,
        102.00,
        104.00,
        106.00,
        108.00,
        110.00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
