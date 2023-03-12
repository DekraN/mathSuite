#ifndef __DISABLEDEEP_BURGERSSOLUTION

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _burgers_viscous_time_exact1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    BURGERS_VISCOUS_TIME_EXACT1 evaluates solution #1 to the Burgers equation.
  Discussion:
    The form of the Burgers equation considered here is
      du       du        d^2 u
      -- + u * -- = nu * -----
      dt       dx        dx^2
    for -1.0 < x < +1.0, and 0 < t.
    Initial conditions are u(x,0) = - sin(M_PI*x).  Boundary conditions
    are u(-1,t) = u(+1,t) = 0.  The viscosity parameter nu is taken
    to be 0.01 / M_PI, although this is not essential.
    The authors note an integral representation for the solution u(x,t),
    and present a better version of the formula that is amenable to
    approximation using Hermite quadrature.
    This program library does little more than evaluate the exact solution
    at a user-specified set of points, using the quadrature rule.
    Internally, the order of this quadrature rule is set to 8, but the
    user can easily modify this value if greater accuracy is desired.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 November 2011
  Author:
    John Burkardt.
  Reference:
    Claude Basdevant, Michel Deville, Pierre Haldenwang, J Lacroix,
    J Ouazzani, Roger Peyret, Paolo Orlandi, Anthony Patera,
    Spectral and finite difference solutions of the Burgers equation,
    Computers and Fluids,
    Volume 14, Number 1, 1986, pages 23-41.
  Parameters:
    Input, double NU, the viscosity.
    Input, int VXN, the number of spatial grid points.
    Input, double VX[VXN], the spatial grid points.
    Input, int VTN, the number of time grid points.
    Input, double VT[VTN], the time grid points.
    Output, double BURGERS_VISCOUS_TIME_EXACT1[VXN*VTN], the solution
    of the Burgers equation at each space and time grid point.
*/
{
	const dt2pititdtpit * const s_data = data;
	
	const register dim_typ vxn = s_data->a0;
	ityp * vt = s_data->a1;
	ityp * vu = s_data->a2;
	ityp nu = s_data->a3;
	const register dim_typ vtn = s_data->a4;
	ityp * vx = s_data->a5;
	
	ityp bot, c;
	dim_typ qi;
	dim_typ qn = 8;
	ityp qw[qn];
	ityp qx[qn];
	dim_typ vti, vxi;

	ityp top;
	/*
	Compute the rule.
	*/
	hermite_ek_compute ( qn, qx, qw );
	/*
	Evaluate U(X,T) for later times.
	*/

	for ( vti = 0; vti < vtn; ++vti )
	{
		if ( ISZERO(vt[vti]))
			for ( vxi = 0; vxi < vxn; ++vxi )
				vu[vxi+vti*vxn] = - sin ( M_PI * vx[vxi] );
		else
		{
			for ( vxi = 0; vxi < vxn; ++vxi )
			{
				top = bot = 0.00;
				for ( qi = 0; qi < qn; ++qi )
				{
					c = 2.00 * sqrt ( nu * vt[vti] );
					top = top - qw[qi] * c * sin ( M_PI * ( vx[vxi] - c * qx[qi] ) ) * exp ( - cos ( M_PI * ( vx[vxi] - c * qx[qi]  ) ) / ( M_2TPI * nu ) );
					bot += qw[qi] * c * exp ( - cos ( M_PI * ( vx[vxi] - c * qx[qi]  ) ) / ( M_2TPI * nu ) );
					vu[vxi+vti*vxn] = top / bot;
				}
			}
		}
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   burgers_viscous_time_exact2 ( ityp nu, const register dim_typ xn, ityp x[static xn], const register dim_typ tn,
  ityp t[static tn], ityp u[static tn])
/******************************************************************************/
/*
  Purpose:
    BURGERS_VISCOUS_TIME_EXACT1 evaluates solution #2 to the Burgers equation.
  Discussion:
    The form of the Burgers equation considered here is
      du       du        d^2 u
      -- + u * -- = nu * -----
      dt       dx        dx^2
    for 0.0 < x < 2 Pi and 0 < t.
    The initial condition is
      u(x,0) = 4 - 2 * nu * dphi(x,0)/dx / phi(x,0)
    where
      phi(x,t) = exp ( - ( x-4*t      ) / ( 4*nu*(t+1) ) )
				+ exp ( - ( x-4*t-2*M_PI ) / ( 4*nu*(t+1) ) )
    The boundary conditions are periodic:
      u(0,t) = u(2 Pi,t)
    The viscosity parameter nu may be taken to be 0.01, but other values
    may be chosen.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2015
  Author:
    John Burkardt.
  Parameters:
    Input, double NU, the viscosity.
    Input, int XN, the number of spatial grid points.
    Input, double X[XN], the spatial grid points.
    Input, int TN, the number of time grid points.
    Input, double T[TN], the time grid points.
    Output, double BURGERS_VISCOUS_TIME_EXACT2[XN*TN], the solution
    of the Burgers equation at each space and time grid point.
*/
{
	dim_typ i, j;
	ityp a, b, c, dphi, phi;

	for ( j = 0; j < tn; ++j )
		for ( i = 0; i < xn; ++i )
		{
			a = ( x[i] - 4.00 * t[j] );
			b = ( x[i] - 4.00 * t[j] - M_2TPI );
			c = 4.00 * nu * ( t[j] + 1.00 );
			phi = exp ( - a * a / c ) + exp ( - b * b / c );
			dphi = - 2.00 * a * exp ( - a * a / c ) / c - 2.00 * b * exp ( - b * b / c ) / c;
			u[i+j*xn] = 4.00 - 2.00 * nu * dphi / phi;
		}
	return;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _hermite_ek_compute ( void * data) 
/******************************************************************************/
/*
  Purpose:
    HERMITE_EK_COMPUTE computes a Gauss-Hermite quadrature rule.
  Discussion:
    The code uses an algorithm by Elhay and Kautsky.
    The abscissas are the zeros of the N-th order Hermite polynomial.
    The integral:
      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
    The quadrature rule:
      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 April 2011
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
  Parameters:
    Input, int N, the number of abscissas.
    Output, double X[N], the abscissas.
    Output, double W[N], the weights.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * w = s_data->a2;
	
	ityp arg;
	ityp bj[n];
	dim_typ i;
	ityp zemu;
	/*
	Define the zero-th moment.
	*/
	arg = 0.50;
	zemu = tgamma ( arg );
	/*
	Define the Jacobi matrix.
	*/

	for ( i = 0; i < n; ++i )
		bj[i] = sqrt ( ( ityp ) ( i + 1 ) / 2.00 );

	for ( i = 0; i < n; ++i )
		x[i] = 0.00;

	w[0] = sqrt ( zemu );
	for ( i = 1; i < n; ++i )
		w[i] = 0.00;
	/*
	Diagonalize the Jacobi matrix.
	*/
	imtqlx ( n, x, bj, w );

	for ( i = 0; i < n; ++i )
		w[i] *= w[i];

	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _imtqlx ( void * data)
/******************************************************************************/
/*
  Purpose:
    IMTQLX diagonalizes a symmetric tridiagonal matrix.
  Discussion:
    This routine is a slightly modified version of the EISPACK routine to
    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
    The authors thank the authors of EISPACK for permission to use this
    routine.
    It has been modified to produce the product Q' * Z, where Z is an input
    vector and Q is the orthogonal matrix diagonalizing the input matrix.
    The changes consist (essentialy) of applying the orthogonal transformations
    directly to Z as they are generated.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 January 2010
  Author:
    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.
  Reference:
    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.
    Roger Martin, James Wilkinson,
    The Implicit QL Algorithm,
    Numerische Mathematik,
    Volume 12, Number 5, December 1968, pages 377-383.
  Parameters:
    Input, int N, the order of the matrix.
    Input/output, double D(N), the diagonal entries of the matrix.
    On output, the information in D has been overwritten.
    Input/output, double E(N), the subdiagonal entries of the
    matrix, in entries E(1) through E(N-1).  On output, the information in
    E has been overwritten.
    Input/output, double Z(N).  On input, a vector.  On output,
    the value of Q' * Z, where Q is the matrix that diagonalizes the
    input symmetric tridiagonal matrix.
*/
{
	static bool result = 2;
	
	const dt3pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * d = s_data->a1;
	ityp * e = s_data->a2;
	ityp * z = s_data->a3;
	
	ityp b;
	ityp c;
	ityp f;
	ityp g;
	int i;
	int ii;
	int itn = 30;
	int j;
	int k;
	int l;
	int m;
	int mml;
	ityp p;
	const register ityp prec = r8_epsilon();
	ityp r;
	ityp s;

	if(n == 1)
	{
		result = IMTQLX_SUCCESS;
		return &result;
	}

	e[n-1] = 0.00;

	for (l = 1; l <= n; ++l)
	{
		j = 0;
		for ( ; ; )
		{
			for (m = l; m <= n; ++m )
			{
				if (m == n)
					break;

				if ( abs ( e[m-1] ) <= prec * ( abs ( d[m-1] ) + abs ( d[m] ) ) )
					break;
			}
			p = d[l-1];
			if (m == l)
				break;
			if ( itn <= j )
			{
				printErr (ERROR_INTERRUPTEDFUNCCALL, "\nIMTQLX - Fatal error!\nIteration limit exceeded\n" );
				result = IMTQLX_ITERATIONS_ERROR;
				return &result;
			}
			++ j;
			g = ( d[l] - p ) / ( 2.00 * e[l-1] );
			r =  sqrt ( g * g + 1.00 );
			g = d[m-1] - p + e[l-1] / ( g + abs ( r ) * r8_sign ( g ) );
			s = 1.00;
			c = 1.00;
			p = 0.00;
			mml = m - l;

			for (ii = 1; ii <= mml; ++ii)
			{
				i = m - ii;
				f = s * e[i-1];
				b = c * e[i-1];

				if (abs (g) <= abs (f))
				{
					c = g / f;
					r =  sqrt (c*c + 1.00);
					e[i] = f * r;
					s = 1.00 / r;
					c *= s;
				}
				else
				{
					s = f / g;
					r =  sqrt ( s*s + 1.00 );
					e[i] = g * r;
					c = 1.00 / r;
					s *= c;
				}
				g = d[i] - p;
				r = ( d[i-1] - g ) * s + 2.00 * c * b;
				p = s * r;
				d[i] = g + p;
				g = c * r - b;
				f = z[i];
				z[i] = s * z[i-1] + c * f;
				z[i-1] = c * z[i-1] - s * f;
			}
			d[l-1] = d[l-1] - p;
			e[l-1] = g;
			e[m-1] = 0.00;
		}
	}
	/*
	Sorting.
	*/
	for (ii = 2; ii <= m; ++ii)
	{
		i = ii - 1;
		k = i;
		p = d[i-1];

		for (j = ii; j <= n; ++j)
			if ( d[j-1] < p )
			{
				k = j;
				p = d[j-1];
			}

		if (k != i)
		{
			d[k-1] = d[i-1];
			d[i-1] = p;
			p = z[i-1];
			z[i-1] = z[k-1];
			z[k-1] = p;
		}
	}

	result = IMTQLX_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _r8vec_even_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    r8VEC_EVEN_NEW returns an r8VEC of values evenly spaced between ALO and AHI.
  Discussion:
    An r8VEC is a vector of r8's.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 February 2004
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of values.
    Input, ityp ALO, AHI, the low and high values.
    Output, ityp r8VEC_EVEN_NEW[N], N evenly spaced values.
*/
{
	const _2itdt * s_data = data;
	
	const register ityp alo = s_data->a0;
	const register ityp ahi = s_data->a1;
	const register dim_typ n = s_data->a2;
	
    ityp *a = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    if ( n == 1 )
        a[0] = 0.50 * ( alo + ahi );
    else
        for (dim_typ i = 1; i <= n; ++i )
            a[i-1] = ( ( ityp ) ( n - i     ) * alo+ ( ityp ) (     i - 1 ) * ahi )/ ( ityp ) ( n     - 1 );

    return a;
}

#endif
