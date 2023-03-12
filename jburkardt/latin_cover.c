#ifndef __DISABLEDEEP_LATINCOVER

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _latin_cover ( void * data)
/******************************************************************************/
/*
  Purpose:
    LATIN_COVER returns a 2D Latin Square Covering.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 August 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, int P[N], a permutation which describes the
    first Latin square.
    Output, int LATIN_COVER[N*N], the Latin cover.  A(I,J) = K
    means that (I,J) is one element of the K-th Latin square.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p = s_data->a1;
	
    int *a;
    dim_typ i;
    int ik;
    int j;
    dim_typ k;

    a = ( int * ) malloc ( n * n * sizeof ( int ) );

    perm_check ( n, p );

    for ( i = 0; i < n; ++i )
        for ( k = 0; k < n; ++k)
        {
            ik = i4_wrap ( i + k, 0, n - 1 );
            j = p[ik] - 1;
            a[i+j*n] = k + 1;
        }

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _latin_cover_2d ( void * data)
/******************************************************************************/
/*
  Purpose:
    LATIN_COVER_2D returns a 2D Latin Square Covering.
  Discussion:
    This procedure has a chance of being extended to M dimensions.
    A basic solution is computed, and the user is permitted to permute
    both the I and J coordinates.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 August 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, int P1[N], P2[N], permutations to be applied
    to the spatial dimensions.
    Output, int LATIN_COVER_2D[N*N], the Latin cover.  A(I,J) = K
    means that (I,J) is one element of the K-th Latin square.
*/
{
	const dt2pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p1 = s_data->a1;
	int * p2 = s_data->a2;
	
    int *a;
    int *b;
    dim_typ i;
    int i1;
    dim_typ j;
    int j1;

    perm_check ( n, p1 );
    perm_check ( n, p2 );

    a = ( int * ) malloc ( n * n * sizeof ( int ) );
    b = ( int * ) malloc ( n * n * sizeof ( int ) );
    /*
    Set up the basic solution.
    */
    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j )
            a[i+j*n] = i4_wrap ( i - j + 1, 1, n );
    /*
    Apply permutation to dimension I.
    */
    for ( i = 0; i < n; ++i)
    {
        i1 = p1[i] - 1;
        for ( j = 0; j < n; ++j)
        b[i1+j*n] = a[i+j*n];
    }
    /*
    Apply permutation to dimension J.
    */
    for ( j = 0; j < n; ++j )
    {
        j1 = p2[j] - 1;
        for ( i = 0; i < n; ++i )
            a[i+j1*n] = b[i+j*n];
    }

    free ( b );

    return a;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _latin_cover_3d ( void * data)
/******************************************************************************/
/*
  Purpose:
    LATIN_COVER_3D returns a 3D Latin Square Covering.
  Discussion:
    A basic solution is computed, and the user is permitted to permute
    I, J and K coordinates.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 June 2012
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of points.
    Input, int P1[N], P2[N], P3[N], permutations to be applied
    to the spatial dimensions.
    Output, int LATIN_COVER_3D[N*N*N], the Latin cover.  A(I,J,K) = L
    means that (I,J,K) is one element of the L-th Latin square.
*/
{
	const dt3pi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * p1 = s_data->a1;
	int * p2 = s_data->a2;
	int * p3 = s_data->a3;
	
    int *a;
    int *b;
    dim_typ i;
    int i1;
    int ik;
    dim_typ j;
    int j1;
    int jk;
    dim_typ k;
    int k1;

    perm_check ( n, p1 );
    perm_check ( n, p2 );
    perm_check ( n, p3 );

    a = ( int * ) malloc ( n * n * n * sizeof ( int ) );
    b = ( int * ) malloc ( n * n * n * sizeof ( int ) );
    /*
    Set up the basic solution.
    */
    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j )
            for ( k = 0; k < n; ++k )
            {
                ik = i4_wrap ( i + 1 - k, 1, n );
                jk = i4_wrap ( j + 1 - k, 1, n );
                b[i+j*n+k*n*n] = ik + ( jk - 1 ) * n;
            }
    /*
    Apply permutation to dimension I.
    */
    for ( i = 0; i < n; ++i )
    {
        i1 = p1[i] - 1;
        for ( j = 0; j < n; ++j )
            for ( k = 0; k < n; ++k)
                a[i1+j*n+k*n*n] = b[i+j*n+k*n*n];
    }
    /*
    Apply permutation to dimension J.
    */
    for ( i = 0; i < n; ++i)
        for ( j = 0; j < n; ++j )
        {
            j1 = p2[j] - 1;
            for ( k = 0; k < n; ++k )
                b[i+j1*n+k*n*n] = a[i+j*n+k*n*n];
        }
    /*
    Apply permutation to dimension K.
    */
    for ( i = 0; i < n; ++i )
        for ( j = 0; j < n; ++j )
            for ( k = 0; k < n; ++k )
            {
                k1 = p3[k] - 1;
                a[i+j*n+k1*n*n] = b[i+j*n+k*n*n];
            }

    free ( b );

    return a;
}

#endif
