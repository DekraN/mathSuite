#ifndef __DISABLEDEEP_COLOREDNOISE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *    _f_alpha ( void * data)
/******************************************************************************/
/*
  Purpose:
    F_ALPHA generates a 1/F^ALPHA noise sequence.
  Discussion:
    Thanks to Miro Stoyanov for pointing out that the second half of
    the data returned by the inverse Fourier transform should be
    discarded, 24 August 2010.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 August 2010
  Author:
    Original C version by Todd Walter.
    This C version by John Burkardt.
  Reference:
    Jeremy Kasdin,
    Discrete Simulation of Colored Noise and Stochastic Processes
    and 1/f^a Power Law Noise Generation,
    Proceedings of the IEEE,
    Volume 83, Number 5, 1995, pages 802-827.
  Parameters:
    Input, int N, the number of samples of the noise to generate.
    Input, double Q_D, the variance of the noise.
    Input, double ALPHA, the exponent for the noise.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double F_ALPHA[N], a sequence sampled with the given
    power law.
*/
{
	const dt2itpipit * const s_data = data;
	const register dim_typ n = s_data->a0;
	register ityp q_d = s_data->a1;
	const register ityp alpha = s_data->a2;
	int * seed = s_data->a3;
	ityp * x = s_data->a4;
	
    ityp h_azero;
    ityp h_a[n];
    ityp h_b[n];
    dim_typ i;
    ityp w_a[n];
    ityp w_b[n];
    ityp w_azero;
    ityp wfa[n<<1];
    ityp hfa[n<<1];
    ityp wi;
    ityp wr;
    ityp *x2;
    const register dim_typ n2 = n<<1;
    /*
    Set the deviation of the noise.
    */
    q_d = sqrt ( q_d );
    /*
    Generate the coefficients Hk.
    */
    hfa[0] = 1.00;
    for ( i = 1; i < n; ++i )
        hfa[i] = hfa[i-1] * ( 0.50 * alpha + ( ityp ) ( i - 1 ) ) / ( ( ityp ) ( i ) );
    for ( i = n; i < n2; ++i )
        hfa[i] = 0.00;
    /*
    Fill Wk with white noise.
    */

    for ( i = 0; i < n; ++i )
        wfa[i] = q_d * r8_normal_01 ( seed );
    for ( i = n; i < n2; ++i)
        wfa[i] = 0.00;
    /*
    Perform the discrete Fourier transforms of Hk and Wk.
    */

    r8vec_sftf ( n2, hfa, &h_azero, h_a, h_b );
    r8vec_sftf (  n2, wfa, &w_azero, w_a, w_b );
    /*
    Multiply the two complex vectors.
    */
    w_azero *= h_azero;

    for ( i = 0; i < n; ++i )
    {
        wr = w_a[i];
        wi = w_b[i];
        w_a[i] = wr * h_a[i] - wi * h_b[i];
        w_b[i] = wi * h_a[i] + wr * h_b[i];
    }
    /*
    This scaling is introduced only to match the behavior
    of the Numerical Recipes code...
    */
    w_azero *= n2;
    for ( i = 0; i < n - 1; ++i )
    {
        w_a[i] *= n;
        w_b[i] *= n;
    }
    i = n - 1;
    w_a[i] *= n2;
    w_b[i] *=  n2;
    /*
    Take the inverse Fourier transform of the result.
    */
    x2 = r8vec_sftb ( n2, w_azero, w_a, w_b );
    /*
    Only return the first N inverse Fourier transform values.
    */
    for ( i = 0; i < n; ++i )
        x[i] = x2[i];

    free ( x2 );
    return NULL;
}

#endif
