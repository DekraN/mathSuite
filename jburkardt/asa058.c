#ifndef __DISABLEDEEP_ASA058

#include "../dutils.h"

 /******************************************************************************/
__MATHSUITE __JBURKARDT  __MATHSUITE __JBURKARDT  void *    _clustr ( void * data)
/******************************************************************************/
/*
  Purpose:
    CLUSTR uses the K-means algorithm to cluster data.
  Discussion:
    Given a matrix of dim[F_OBSERVATIONS] observations on J variables, the
    observations are allocated to N clusters in such a way that the
    within-cluster sum of squares is minimised.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    29 October 2010
  Author:
    Original FORTRAN77 version by David Sparks.
    C version by John Burkardt.
  Reference:
    David Sparks,
    Algorithm AS 58:
    Euclidean Cluster Analysis,
    Applied Statistics,
    Volume 22, Number 1, 1973, pages 126-130.
  Parameters:
    Input, double X[I*J], the observed data.
    Input/output, double D[K*J], the cluster centers.
    On input, the user has chosen these.  On output, they have been
    updated.
    Output, double DEV[K], the sums of squared deviations
    of observations from their cluster centers.
    Output, int B[I], indicates the cluster to which
    each observation has been assigned.
    Workspace, double F[I].
    Output, int E[K], the number of observations assigned
    to each cluster.
    Input, int I, the number of observations.
    Input, int J, the number of variables.
    Input, int N, the number of clusters.
    Input, int NZ, the minimum number of observations
    which any cluster is allowed to have.
    Input, int K, the maximum number of clusters.
*/
{
	const pdt5pitpdt2dt * const s_data = data;
	const dim_typ * dim = s_data->a0;
	ityp * x = s_data->a1;
	ityp * d = s_data->a2;
	ityp * dev = s_data->a3;
	ityp * b = s_data->a4;
	ityp * f = s_data->a5;
	dim_typ * e = s_data->a6;
	const register dim_typ n = s_data->a7;
	const register dim_typ nz = s_data->a8;
	
	ityp da;
	ityp db;
	ityp dc;
	ityp de;
	ityp fl;
	ityp fm;
	ityp fq;
	dim_typ ia;
	dim_typ ic;
	dim_typ id;
	dim_typ ie;
	dim_typ ig;
	dim_typ ih;
	dim_typ ii;
	dim_typ ij;
	dim_typ ik;
	dim_typ il;
	dim_typ in;
	dim_typ ip;
	dim_typ ir;
	dim_typ is;
	dim_typ it;
	dim_typ iu;
	dim_typ iw;
	dim_typ ix;
	dim_typ iy;

	#pragma omp parallel for
	for(ia = 0; ia < n; ++ia )
		e[ia] = 0;
	/*
	For each observation, calculate the distance from each cluster
	center, and assign to the nearest.
	*/
	for ( ic = 1; ic <= dim[F_OBSERVATIONS]; ++ic )
	{
		f[(ic-1)*dim[F_OBSERVATIONS]] = 0.00;
		da = 1.0E+10;

		for ( id = 1; id <= n; id++ )
		{
			db = 0.00;
			for ( ie = 1; ie <= dim[F_VARIABLES]; ie++ )
			{
				dc = x[(ic-1)*dim[F_VARIABLES]+(ie-1)] - d[(id-1)*dim[F_VARIABLES]+(ie-1)];
				db += dc * dc;
			}

			if ( db < da )
			{
				da = db;
				b[ic-1] = id;
			}
		}
		ig = b[(ic-1)*dim[F_OBSERVATIONS]];
		++ e[(ig-1)*dim[F_MAXCLUSTERS]];
	}
	/*
	Calculate the mean and sum of squares for each cluster.
	*/
	for ( ix = 1; ix <= n; ++ix )
	{
		dev[(ix-1)*dim[F_MAXCLUSTERS]] = 0.00;
		#pragma omp parallel for
		for ( iy = 1; iy <= dim[F_VARIABLES]; ++iy )
			d[(ix-1)*dim[F_VARIABLES]+(iy-1)] = 0.00;
	}

	for ( ic = 1; ic <= dim[F_OBSERVATIONS]; ++ic )
	{
		ig = b[ic-1];
		#pragma omp parallel for
		for ( ih = 1; ih <= dim[F_VARIABLES]; ++ih )
			d[(ig-1)*dim[F_VARIABLES]+(ih-1)] += x[(ic-1)*dim[F_VARIABLES]+(ih-1)];
	}

	#pragma omp parallel for
	for ( ij = 1; ij <= dim[F_VARIABLES]; ++ij )
		#pragma omp parallel for
		for ( ii = 1; ii <= n; ++ii )
			d[(ii-1)*dim[F_VARIABLES]+(ij-1)] /= ( ityp ) e[ii-1];

	for ( ij = 1; ij <= dim[F_VARIABLES]; ++ij )
		for ( ik = 1; ik <= dim[F_OBSERVATIONS]; ++ik )
		{
			il = b[(ik-1)*dim[F_OBSERVATIONS]];
			da = x[(ik-1)*dim[F_VARIABLES]+(ij-1)] - d[(il-1)*dim[F_VARIABLES]+(ij-1)];
			db = da * da;
			f[(ik-1)*dim[F_OBSERVATIONS]] = (dev[(il-1)*dim[F_MAXCLUSTERS]] += db);
		}

	for ( ik = 1; ik <= dim[F_OBSERVATIONS]; ++ik )
	{
		il = b[(ik-1)*dim[F_OBSERVATIONS]];
		fl = e[(il-1)*dim[F_MAXCLUSTERS]];
		if ( 2.00 <= fl )
			f[(ik-1)*dim[F_OBSERVATIONS]] *= fl / ( fl - 1.00 );
	}
	/*
	Examine each observation in turn to see if it should be
	reassigned to a different cluster.
	*/
	for ( ; ; )
	{
		iw = 0;

		for ( ik = 1; ik <= dim[F_OBSERVATIONS]; ++ik )
		{
			il = b[(ik-1)*dim[F_OBSERVATIONS]];
			ir = il;
			/*
			If the number of cluster points is less than or equal to the
			specified minimum, NZ, then bypass this iteration.
			*/
			if ( nz < e[(il-1)*dim[F_MAXCLUSTERS]] )
			{
				fl = e[(il-1)*dim[F_MAXCLUSTERS]];
				dc = f[(ik-1)*dim[F_OBSERVATIONS]];

				for ( in = 1; in <= n; ++in )
				{
					if ( in != il )
					{
						fm = e[(in-1)*dim[F_MAXCLUSTERS]];
						fm /= ( fm + 1.00 );

						de = 0.00;
						for ( ip = 1; ip <= dim[F_VARIABLES]; ++ip )
						{
							da = x[(ik-1)*dim[F_VARIABLES]+(ip-1)] - d[(in-1)*dim[F_VARIABLES]+(ip-1)];
							de += da * da * fm;
						}

						if ( de < dc )
						{
							dc = de;
							ir = in;
						}
					}
				}
				/*
				Reassignment is made here if necessary.
				*/
				if ( ir != il )
				{
					fq = e[(ir-1)*dim[F_MAXCLUSTERS]];
					dev[(il-1)*dim[F_MAXCLUSTERS]] -= f[(ik-1)*dim[F_OBSERVATIONS]];
					dev[(ir-1)*dim[F_MAXCLUSTERS]] += dc;
					++ e[(ir-1)*dim[F_MAXCLUSTERS]];
					-- e[(il-1)*dim[F_MAXCLUSTERS]];

					for ( is = 1; is <= dim[F_VARIABLES]; ++is )
					{
						d[(il-1)*dim[F_VARIABLES]+(is-1)] = ( d[(il-1)*dim[F_VARIABLES]+(is-1)]
						* fl - x[(ik-1)*dim[F_VARIABLES]+(is-1)] ) / ( fl - 1.00 );
						d[(ir-1)*dim[F_VARIABLES]+(is-1)] = ( d[(ir-1)*dim[F_VARIABLES]+(is-1)]
						* fq + x[(ik-1)*dim[F_VARIABLES]+(is-1)] ) / ( fq + 1.00 );
					}

					b[(ik-1)*dim[F_OBSERVATIONS]] = ir;

					for ( it = 1; it <= dim[F_OBSERVATIONS]; ++it )
					{
						ij = b[(it-1)*dim[F_OBSERVATIONS]];

						if ( ij == il || ij == ir )
						{
							f[it-1] = 0.00;
							for ( iu = 1; iu <= dim[F_VARIABLES]; ++iu )
							{
								da = x[(it-1)*dim[F_VARIABLES]+(iu-1)] - d[(ij-1)*dim[F_VARIABLES]+(iu-1)];
								f[it-1] += da * da;
							}
							fl = e[(ij-1)*dim[F_MAXCLUSTERS]];
							f[(it-1)*dim[F_OBSERVATIONS]] *= fl / ( fl - 1.0 );
							}
					}
					++ iw;
				}
			}
		}
		/*
		If any reassignments were made on this pass, then do another pass.
		*/
		if (!iw)
			break;
	}
	return NULL;
}

#endif
