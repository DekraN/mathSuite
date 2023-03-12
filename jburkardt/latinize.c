#ifndef __DISABLEDEEP_LATINIZE

#include "../dutils.h"

/*******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_latinize ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_LATINIZE "latinizes" a real table.
  Discussion:
    On output, each row of the table will have the properties that:
    1) the minimum and maximum row values are the same as on input;
    2) the row contains N evenly spaced values between the
       minimum and maximum;
    3) in each row, the elements retain their ordering.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    04 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, the spatial dimension.
    Input, int N, the number of columns.
    Input/output, double TABLE[M*N].  On input, the dataset to
    be "latinized".  On output, the latinized dataset.
*/
{
	const _2dtpit * const s_data = data; 
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * table = s_data->a2;
	
    dim_typ i;
    int *indx;
    dim_typ j;
    ityp *v;
    ityp v_max;
    ityp v_min;

    v = ( double * ) malloc ( n * sizeof ( double ) );

    for ( i = 0; i < m; ++i )
    {
        for ( j = 0; j < n; ++j )
            v[j] = table[i+j*m];
        _MIN ( &v_min, n, v );
        _MAX ( &v_max, n, v );
        indx = r8vec_sort_heap_index_a_new ( n, v );

        for ( j = 0; j < n; ++j )
            table[i+indx[j]*m] = ( ( ityp ) ( n - j - 1 ) * v_min   + ( ityp ) (     j     ) * v_max ) / ( ityp ) ( n     - 1 );
        free ( indx );
    }
    free ( v );
    return NULL;
}

#endif
