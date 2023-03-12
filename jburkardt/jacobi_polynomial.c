#ifndef __DISABLEDEEP_JACOBIPOLYNOMIAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _j_double_product_integral( void * data)
/******************************************************************************/
/*
  Purpose:
    J_DOUBLE_PRODUCT_INTEGRAL: integral of J(i,x)*J(j,x)*(1-x)^a*(1+x)^b.
  Discussion:
    VALUE = integral ( -1 <= x <= +1 ) J(i,x)*J(j,x)*(1-x)^a*(1+x)^b dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int I, J, the polynomial indices.
    Input, double A, B, the parameters.
    -1 < A, B.
    Output, double VALUE, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const pdt2it * const s_data = data;
	dim_typ * dim = s_data->a0;
	const register ityp a = s_data->a1;
	const register ityp b = s_data->a2;
	
	ityp value;
	if (dim[XROW] != dim[YROW])
		value= 0.00;
	else
	{
		const register ityp i_r8 = (ityp)(dim[XROW]);
		value=pow ( 2, dim[XROW] + dim[YROW] + 1.00 ) / (2.00*i_r8+dim[XROW]+dim[YROW]+1.00)*tgamma(i_r8+dim[XROW]+1.00)*tgamma(i_r8+dim[YROW]+1.00)/r8_factorial(dim[XROW])/tgamma(i_r8+dim[XROW]+dim[YROW]+1.00);
	}
	
	result = value;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _j_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    J_INTEGRAL evaluates a monomial integral associated with J(n,a,b,x).
  Discussion:
    The integral:
      integral ( -1 <= x < +1 ) x^n dx
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Parameters:
    Input, int N, the exponent.
    0 <= N.
    Output, double J_INTEGRAL, the value of the integral.
*/
{
	static ityp result = MAX_VAL;
	
	const register dim_typ n = *(dim_typ *) data;
	
	result = ((n%2) == 1 ? 0.00 : 2.00/(ityp)(n+1));
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _j_polynomial ( void * data)
/******************************************************************************/
/*
  Purpose:
    JACOBI_POLY evaluates the Jacobi polynomial J(n,a,b,x).
  Differential equation:
 (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
  Recursion:
    P(0,ALPHA,BETA,X) = 1,
    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
    P(N,ALPHA,BETA,X)  =
  (
     (2*N+ALPHA+BETA-1)
        * ((ALPHA^2-BETA**2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X)
        * P(N-1,ALPHA,BETA,X)
        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
  Restrictions:
    -1 < ALPHA
    -1 < BETA
  Norm:
    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA
      * P(N,ALPHA,BETA,X)^2 dX
    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
  ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
  Special values:
    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  Parameters:
    Input, int M, the number of evaluation points.
    Input, int N, the highest order polynomial to compute.  Note
    that polynomials 0 through N will be computed.
    Input, double ALPHA, one of the parameters defining the Jacobi
    polynomials, ALPHA must be greater than -1.
    Input, double BETA, the second parameter defining the Jacobi
    polynomials, BETA must be greater than -1.
    Input, double X[M], the evaluation points.
    Output, double J_POLYNOMIAL[M*(N+1)], the values.
*/
{
	static ityp result = MAX_VAL;
	
	const dt3it * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp alpha = s_data->a1;
	const register ityp beta = s_data->a2;
	const register ityp x = s_data->a3;
	
	ityp c1;
	ityp c2;
	ityp c3;
	ityp c4;

	if(alpha <= -1.00 || beta <= -1.00 || n < 0 )
	{
		result = JACOBI_INVALIDRETURNVALUE;
		return &result;
	}

	ityp v[n+1];

	v[0] = 1.00;

	if(!n)
	{
		result = 1.00;
		return &result;
	}

	v[1] = (1.00+0.50*(alpha+beta ) )*x+0.50*(alpha-beta);

	for (dim_typ i=2; i <= n; ++i)
	{
		c1 = 2.00*(ityp)(i)*((ityp)(i)+alpha+beta)*((ityp)((i<<1)-2)+alpha+beta);
		c2 = ((ityp)((i<<1)-1)+alpha+beta)*((ityp)((i<<1))+alpha+beta)*((ityp)((i<<1)-2)+alpha+beta);
		c3 = ( (ityp)((i<<1)-1)+alpha+beta)*(alpha+beta)*(alpha-beta);
		c4 = -(ityp)(2)*((ityp)(i-1)+alpha)*((ityp)(i-1)+beta)*((ityp)((i<<1))+alpha+beta);
		v[i] = ((c3+c2*x)*v[i-1]+c4*v[i-2])/c1;
	}

	result = v[n];
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _j_polynomial_zeros ( void * data)
/******************************************************************************/
/*
  Purpose:
    J_POLYNOMIAL_ZEROS: zeros of Jacobi polynomial J(n,a,b,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int, N, the order.
    Input, double, ALPHA, BETA, the parameters.
    -1 < ALPHA, BETA.
    Output, double J_POLYNOMIAL_ZEROS[N], the zeros.
*/
{
	static bool result = 2;
	
	const dt2itpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp alpha = s_data->a1;
	const register ityp beta = s_data->a2;
	ityp * x = s_data->a3;
	
	ityp a2b2;
	ityp ab = alpha+beta;
	ityp abi = 2.00 + ab;
	int i;
	ityp i_r8;

	/*
	Define the zero-th moment.
	*/
	const register ityp zemu = pow(2.00,ab+1.00)*tgamma(alpha+1.00)*tgamma(beta+1.00)/tgamma(abi);

	/*
	Define the Jacobi matrix.
	*/

	x[0] = ( beta - alpha ) / abi;

	ityp bj[n];
	ityp w[n];

	bj[0] = 4.00*(1.00+alpha*(1.00+beta)/((abi+1.00)*abi*abi));
	a2b2 = beta*beta - alpha*alpha;



	#pragma omp parallel for
	for (i = 1; i < n; ++i)
	{
		i_r8 = (ityp)(i+1);
		abi = 2.00*i_r8 + ab;
		x[i] = a2b2/((abi-2.00)*abi);
		abi = abi*abi;
		w[i] = x[i] = 0.00;
		bj[i] = sqrt(4.00*i_r8*(i_r8+alpha)*(i_r8+beta)*(i_r8+ab)/((abi-1.00)*abi));
	}

	bj[0] = sqrt(bj[0]);
	w[0] = sqrt (zemu);
	/*
	Diagonalize the Jacobi matrix.
	*/
	
	result = imtqlx ( n, x, bj, w );
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _j_quadrature_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    J_QUADRATURE_RULE: Gauss-Jacobi quadrature based on J(n,a,b,x).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 August 2013
  Author:
    John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int, N, the order.
    Input, double, ALPHA, BETA, the parameters.
    -1 < ALPHA, BETA.
    Output, double X[N], the abscissas.
    Output, double W[N], the weights.
*/
{
	static bool result = 2;
	
	const dt2it2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp alpha = s_data->a1;
	const register ityp beta = s_data->a2;
	ityp * x = s_data->a3;
	ityp * w= s_data->a4;
	
	ityp a2b2;
	ityp ab = alpha+beta;
	ityp abi = 2.00 + ab;
	dim_typ i;
	ityp i_r8;

	/*
	Define the zero-th moment.
	*/
	const register ityp zemu = pow(2.00,ab+1.00)*tgamma(alpha+1.00)*tgamma(beta+1.00)/tgamma(abi);
	ityp bj[n];

	/*
	Define the Jacobi matrix.
	*/
	x[0] = (beta-alpha)/abi;
	bj[0] = 4.00*(1.00+alpha)*(1.00+beta)/((abi+1.00)*abi*abi);

	a2b2 = beta*beta - alpha*alpha;

	#pragma omp parallel for
	for (i = 1; i < n; ++i )
	{
		i_r8 = (ityp) (i+1);
		abi = 2.00*i_r8 + ab;
		x[i] = a2b2 / ((abi-2.0)*abi);
		abi *= abi;
		w[i] = x[i] = 0.00;
		bj[i] = sqrt(4.00*i_r8*(i_r8+alpha)*(i_r8+beta)*(i_r8+ab)/((abi-1.00)*abi));
	}

	w[0] = sqrt (zemu);
	/*
	Diagonalize the Jacobi matrix.
	*/
	const bool status = imtqlx ( n, x, bj, w );

	#pragma omp parallel for
	for ( i = 0; i < n; ++i)
		w[i] *= w[i];

	result = status;
	return &result;
}

#endif
