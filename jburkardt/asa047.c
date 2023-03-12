#ifndef __DISABLEDEEP_ASA047

#include "../dutils.h"

// Modified Simplex-Method Based for Scalar functions (f:R^n-->R)
/******************************************************************************/
__MATHSUITE __JBURKARDT  void * nelmin ( ityp fn ( ityp * ), const register dim_typ n, ityp * start, ityp xmin[static n],
ityp *ynewlo, const register ityp reqmin, ityp step[static n], dim_typ konvge, dim_typ kcount, dim_typ *icount, dim_typ *numres)
/******************************************************************************/
/*
  Purpose:
    NELMIN minimizes a function using the Nelder-Mead algorithm.
  Discussion:
    This routine seeks the minimum value of a user-specified function.
    Simplex function minimisation procedure due to Nelder+Mead(1965),
    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
    25, 97) and Hill(1978, 27, 380-2)
    The function to be minimized must be defined by a function of
    the form
      function fn ( x, f )
      double fn
      double x(*)
    and the name of this subroutine must be declared EXTERNAL in the
    calling routine and passed as the argument FN.
    This routine does not include a termination test using the
    fitting of a quadratic surface.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    28 October 2010
  Author:
    Original FORTRAN77 version by R ONeill.
    C version by John Burkardt.
  Reference:
    John Nelder, Roger Mead,
    A simplex method for function minimization,
    Computer Journal,
    Volume 7, 1965, pages 308-313.
    R ONeill,
    Algorithm AS 47:
    Function Minimization Using a Simplex Procedure,
    Applied Statistics,
    Volume 20, Number 3, 1971, pages 338-345.
  Parameters:
    Input, double FN ( double x[] ), the name of the routine which evaluates
    the function to be minimized.
    Input, int N, the number of variables.
    Input/output, double START[N].  On input, a starting point
    for the iteration.  On output, this data may have been overwritten.
    Output, double XMIN[N], the coordinates of the point which
    is estimated to minimize the function.
    Output, double YNEWLO, the minimum value of the function.
    Input, double REQMIN, the terminating limit for the variance
    of function values.
    Input, double STEP[N], determines the size and shape of the
    initial simplex.  The relative magnitudes of its elements should reflect
    the units of the variables.
    Input, int KONVGE, the convergence check is carried out
    every KONVGE iterations.
    Input, int KCOUNT, the maximum number of function
    evaluations.
    Output, int *ICOUNT, the number of function evaluations
    used.
    Output, int *NUMRES, the number of restarts.
    Output, int *IFAULT, error indicator.
    0, no errors detected.
    1, REQMIN, N, or KONVGE has an illegal value.
    2, iteration terminated because KCOUNT was exceeded without convergence.
*/
{
	static sel_typ result = UCHAR_MAX;
	
	ityp del;
	ityp dn;
	ityp dnn;
	dim_typ i;
	dim_typ ihi;
	dim_typ ilo;
	dim_typ j;
	dim_typ jcount;
	dim_typ l;
	dim_typ nn;
	ityp rq;
	ityp x;
	ityp y2star;
	ityp ylo;
	ityp ystar;
	ityp z;
	/*
	Check the input parameters.
	*/
	if ( reqmin <= 0.0 )
	{
		result = NELMIN_INVALIDREQMIN;
		return &result;
	}

	if ( n < 1 )
	{
		result = NELMIN_INVALIDDIMENSION;
		return &result;
	}

	if ( konvge < 1 )
	{
		result = NELMIN_INVALIDKONVGE;
		return &result;
	}

	ityp pstar[n];
	ityp p2star[n];
	ityp pbar[n];
	ityp y[n+1];

	bool iflag = false;

	ityp * const p = malloc ( n * ( n + 1 ) * sizeof ( ityp ) );
	errMemEx(p, NELMIN_ALLOC_ERROR);

	*icount = 0;
	*numres = 0;

	jcount = konvge;
	dn = (ityp)(n);
	nn += 1;
	dnn = (ityp)(nn);
	del = 1.00;
	rq = reqmin * dn;
	/*
	Initial or restarted loop.
	*/
	for ( ; ; )
	{
		#pragma omp parallel for
		for ( i = 0; i < n; ++i )
			p[i+n*n] = start[i];
		y[n] = fn (start);
		++ (*icount);

		for ( j = 0; j < n; j++ )
		{
			x = start[j];
			start[j] += step[j] * del;
			#pragma omp parallel for
			for ( i = 0; i < n; ++i)
				p[i+j*n] = start[i];
			y[j] = fn ( start );
			++ (*icount);
			start[j] = x;
		}
		/*
		The simplex construction is complete.

		Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
		the vertex of the simplex to be replaced.
		*/
		ylo = y[0];
		ilo = 0;

		for ( i = 1; i < nn; i++ )
		{
			if ( y[i] < ylo )
			{
				ylo = y[i];
				ilo = i;
			}
		}
		/*
		Inner loop.
		*/
		for ( ; ; )
		{
			if ( kcount <= *icount )
				break;
			*ynewlo = y[0];
			ihi = 0;

			for ( i = 1; i < nn; i++ )
			{
				if ( *ynewlo < y[i] )
				{
					*ynewlo = y[i];
					ihi = i;
				}
			}
			/*
			Calculate PBAR, the centroid of the simplex vertices
			excepting the vertex with Y value YNEWLO.
			*/
			for ( i = 0; i < n; ++i )
			{
				z = 0.00;
				for ( j = 0; j < nn; ++j )
					z += p[i+j*n];
				z -= p[i+ihi*n];
				pbar[i] = z / dn;
			}
			/*
			Reflection through the centroid.
			*/
			#pragma omp parallel for
			for ( i = 0; i < n; ++i )
				pstar[i] = pbar[i] + 1.00 * ( pbar[i] - p[i+ihi*n] );
			ystar = fn ( pstar );
			++ (*icount);
			/*
			Successful reflection, so extension.
			*/
			if ( ystar < ylo )
			{
				#pragma omp parallel for
				for ( i = 0; i < n; ++i )
					p2star[i] = pbar[i] + 2.00 * ( pstar[i] - pbar[i] );
				y2star = fn ( p2star );
				++ (*icount);
				/*
				Check extension.
				*/
				if ( ystar < y2star )
				{
					#pragma omp parallel for
					for ( i = 0; i < n; ++i )
						p[i+ihi*n] = pstar[i];

					y[ihi] = ystar;
				}
				/*
				Retain extension or contraction.
				*/
				else
				{
					#pragma omp parallel for
					for ( i = 0; i < n; ++i )
						p[i+ihi*n] = p2star[i];
					y[ihi] = y2star;
				}
			}
			/*
			No extension.
			*/
			else
			{
				l = 0;
				#pragma omp parallel for
				for ( i = 0; i < nn; ++ i)
					if ( ystar < y[i] )
						++ l;

				if ( 1 < l )
				{
					#pragma omp parallel for
					for ( i = 0; i < n; ++i )
						p[i+ihi*n] = pstar[i];
					y[ihi] = ystar;
				}
				/*
				Contraction on the Y(IHI) side of the centroid.
				*/
				else if (!l)
				{
					#pragma omp parallel for
					for ( i = 0; i < n; ++i )
						p2star[i] = pbar[i] + 0.50 * ( p[i+ihi*n] - pbar[i] );
					y2star = fn ( p2star );
					++ (*icount);
					/*
					Contract the whole simplex.
					*/
					if ( y[ihi] < y2star )
					{
						for ( j = 0; j < nn; ++j )
						{
							#pragma omp parallel for
							for ( i = 0; i < n; ++i )
							{
								p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
								xmin[i] = p[i+j*n];
							}
							y[j] = fn ( xmin );
							++ (*icount);
						}
						ylo = y[0];
						ilo = 0;

						for ( i = 1; i < nn; ++i )
						{
							if ( y[i] < ylo )
							{
								ylo = y[i];
								ilo = i;
							}
						}
						continue;
					}
					/*
					Retain contraction.
					*/
					else
					{
						#pragma omp parallel for
						for ( i = 0; i < n; ++i)
							p[i+ihi*n] = p2star[i];
						y[ihi] = y2star;
					}
				}
				/*
				Contraction on the reflection side of the centroid.
				*/
				else if ( l == 1 )
				{
					#pragma omp parallel for
					for ( i = 0; i < n; ++i )
						p2star[i] = pbar[i] + 0.50 * ( pstar[i] - pbar[i] );
					y2star = fn ( p2star );
					++ (*icount);
					/*
					Retain reflection?
					*/
					if ( y2star <= ystar )
					{
						#pragma omp parallel for
						for ( i = 0; i < n; ++i )
							p[i+ihi*n] = p2star[i];
						y[ihi] = y2star;
					}
					else
					{
						#pragma omp parallel for
						for ( i = 0; i < n; ++i )
							p[i+ihi*n] = pstar[i];
						y[ihi] = ystar;
					}
				}
			}
			/*
			Check if YLO improved.
			*/
			if ( y[ihi] < ylo )
			{
				ylo = y[ihi];
				ilo = ihi;
			}
			-- jcount;

			if ( 0 < jcount )
				continue;
			/*
			Check to see if minimum reached.
			*/
			if ( *icount <= kcount )
			{
				jcount = konvge;
				z = 0.00;
				for ( i = 0; i < nn; ++i)
					z += y[i];
				x = z / dnn;

				z = 0.0;
				for ( i = 0; i < nn; ++i )
					z += pow ( y[i] - x, 2 );

				if ( z <= rq )
					break;
			}
		}
		/*
		Factorial tests to check that YNEWLO is a local minimum.
		*/
		#pragma omp parallel for
		for ( i = 0; i < n; ++i)
			xmin[i] = p[i+ilo*n];
		*ynewlo = y[ilo];

		if ( kcount < *icount )
		{
			iflag = true;
			break;
		}

		for ( i = 0; i < n; i++ )
		{
			del = step[i] * 0.001;
			xmin[i] = xmin[i] + del;
			z = fn ( xmin );
			++ (*icount);
			if ( z < *ynewlo )
			{
				iflag = true;
				break;
			}
			xmin[i] -= del + del;
			z = fn ( xmin );
			++ (*icount);
			if ( z < *ynewlo )
			{
				iflag = true;
				break;
			}
			xmin[i] += del;
		}

		if (!iflag)
			break;

		/*
		Restart the procedure.
		*/
		#pragma omp parallel for
		for ( i = 0; i < n; ++i )
			start[i] = xmin[i];
		del = 0.001;
		++ (*numres);
	}

	free (p);

	result = NELMIN_SUCCESS+iflag;
	return &result;
}
#endif
