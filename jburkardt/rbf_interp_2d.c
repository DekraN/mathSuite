#ifndef __DISABLEDEEP_RBDINTERP2D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   *rbf_interp_2d ( const register dim_typ m, const register dim_typ nd, ityp xd[static m*nd], const register ityp r0, void phix ( int n, ityp r[], ityp r0, ityp v[] ), ityp w[static nd], const register dim_typ ni, ityp xi[static m*ni] )
/******************************************************************************/
/*
  Purpose:
    RBF_INTERP evaluates a radial basis function interpolant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 June 2012
  Author:
    John Burkardt
  Reference:
    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
    Third Edition,
    Cambridge University Press, 2007,
    ISBN13: 978-0-521-88068-8,
    LC: QA297.N866.
  Parameters:
    Input, int M, the spatial dimension.
    Input, int ND, the number of data points.
    Input, double XD[M*ND], the data points.
    Input, double R0, a scale factor.  R0 should be larger than the typical
    separation between points, but smaller than the maximum separation.
    The value of R0 has a significant effect on the resulting interpolant.
    Input, void PHI ( int N, double R[], double R0, double V[] ), a
    function to evaluate the radial basis functions.
    Input, double W[ND], the weights, as computed by RBF_WEIGHTS.
    Input, int NI, the number of interpolation points.
    Input, double XI[M*NI], the interpolation points.
    Output, double RBF_INTERP[NI], the interpolated values.
*/
{
    ityp *fi;
    dim_typ i, j, k;
    ityp *r;
    ityp *v;

    fi = ( ityp * ) malloc ( ni * sizeof ( ityp ) );
    r = ( ityp * ) malloc ( nd * sizeof ( ityp ) );
    v = ( ityp * ) malloc ( nd * sizeof ( ityp ) );

    for ( i = 0; i < ni; ++i )
    {
        for ( j = 0; j < nd; ++j )
        {
            r[j] = 0.00;
            for ( k = 0; k < m; ++k )
                r[j] += pow ( xi[k+i*m] - xd[k+j*m], 2 );
            r[j] = sqrt ( r[j] );
        }
        phix ( nd, r, r0, v );
        fi[i] = r8vec_dot_product ( nd, v, w );
    }

    free ( r );
    free ( v );

    return fi;
}
#endif
