#ifndef __DISABLEDEEP_SIMPLEXCOORDINATES

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _simplex_coordinates1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_COORDINATES1 computes the Cartesian coordinates of simplex vertices.
  Discussion:
    The simplex will have its centroid at 0;
    The sum of the vertices will be zero.
    The distance of each vertex from the origin will be 1.
    The length of each edge will be constant.
    The dot product of the vectors defining any two vertices will be - 1 / N.
    This also means the angle subtended by the vectors from the origin
    to any two distinct vertices will be acos ( - 1 / N ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the spatial dimension.

    Output, double SIMPLEX_COORDINATES1[N*(N+1)], the coordinates of the vertices
    of a simplex in N dimensions.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i;
    dim_typ ii;
    dim_typ j;
    ityp s;
    ityp *x;

    x = r8mat_zero_new ( n, n + 1 );

    for ( i = 0; i < n; ++i )
    {
        /*
        Set X(I,I) so that sum ( X(1:I,I)**2 ) = 1.
        */
        s = 0.00;
        for ( ii = 0; ii < i; ++ii )
            s += x[ii+i*n] * x[ii+i*n];
        x[i+i*n] = sqrt ( 1.00 - s );
        /*
        Set X(I,J) for J = I+1 to N+1 by using the fact that XI dot XJ = - 1 / N
        */
        for ( j = i + 1; j < n + 1; ++j )
        {
            s = 0.00;
            for ( ii = 0; ii < i; ++ii)
                s += x[ii+i*n] * x[ii+j*n];
            x[i+j*n] = ( - 1.00 / ( ityp ) ( n ) - s ) / x[i+i*n];
        }
    }

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _simplex_coordinates2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_COORDINATES2 computes the Cartesian coordinates of simplex vertices.
  Discussion:
    This routine uses a simple approach to determining the coordinates of
    the vertices of a regular simplex in n dimensions.
    We want the vertices of the simplex to satisfy the following conditions:
    1) The centroid, or average of the vertices, is 0.
    2) The distance of each vertex from the centroid is 1.
       By 1), this is equivalent to requiring that the sum of the squares
       of the coordinates of any vertex be 1.
    3) The distance between any pair of vertices is equal (and is not zero.)
    4) The dot product of any two coordinate vectors for distinct vertices
       is -1/N; equivalently, the angle subtended by two distinct vertices
       from the centroid is acos ( -1/N).
    Note that if we choose the first N vertices to be the columns of the
    NxN identity matrix, we are almost there.  By symmetry, the last column
    must have all entries equal to some value A.  Because the square of the
    distance between the last column and any other column must be 2 (because
    that's the distance between any pair of columns), we deduce that
 (A-1)^2 + (N-1)*A^2 = 2, hence A = (1-sqrt(1+N))/N.  Now compute the
    centroid C of the vertices, and subtract that, to center the simplex
    around the origin.  Finally, compute the norm of one column, and rescale
    the matrix of coordinates so each vertex has unit distance from the origin.
    This approach devised by John Burkardt, 19 September 2010.  What,
    I'm not the first?
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the spatial dimension.
    Output, double SIMPLEX_COORDINATES2[N*(N+1)], the coordinates of the vertices
    of a simplex in N dimensions.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    ityp a;
    ityp c;
    dim_typ i, j;
    ityp s;
    ityp *x;

    x = r8mat_zero_new ( n, n + 1 );

    for ( i = 0; i < n; ++i )
        x[i+i*n] = 1.00;

    a = ( 1.00 - sqrt ( 1.0 + ( ityp ) ( n ) ) ) / ( ityp ) ( n );

    for ( i = 0; i < n; ++i )
        x[i+n*n] = a;
    /*
    Now adjust coordinates so the centroid is at zero.
    */
    for ( i = 0; i < n; ++i )
    {
        c = 0.00;
        for ( j = 0; j < n + 1; ++j)
            c += x[i+j*n];
        c /= ( ityp ) ( n + 1 );
        for ( j = 0; j < n + 1; ++j )
        x[i+j*n] -= c;
    }
    /*
    Now scale so each column has norm 1.
    */
    s = 0.00;
    for ( i = 0; i < n; ++i)
        s += x[i+0*n] * x[i+0*n];
    s = sqrt ( s );

    for ( j = 0; j < n + 1; ++j )
        for ( i = 0; i < n; ++i )
            x[i+j*n] /= s;

    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _simplex_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMPLEX_VOLUME computes the volume of a simplex.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int N, the spatial dimension.
    Input, double X[N*(N+1)], the coordinates of the vertices
    of a simplex in N dimensions.
    Output, double SIMPLEX_VOLUME, the volume of the simplex.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	
    ityp *a;
    ityp det;
    dim_typ i, j;
    ityp volume;

    a = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
    for ( j = 0; j < n; j++ )
        for ( i = 0; i < n; ++i )
            a[i+j*n] = x[i+j*n];
    for ( j = 0; j < n; ++j )
        for ( i = 0; i < n; ++i )
            a[i+j*n] -= x[i+n*n];

    det = r8mat_det ( n, a );
    volume = abs ( det );
    for ( i = 1; i <= n; ++i )
        volume /= ( ityp ) ( i );
    free ( a );

	result = volume;
    return &result;
}

#endif
