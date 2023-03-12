#ifndef __DISABLEDEEP_IMAGEDENOISE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _gray_median_news ( void * data)
/******************************************************************************/
/*
  Purpose:
    GRAY_MEDIAN_NEWS uses a median NEWS filter on a gray scale image to remove noise.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 February 2012
  Author:
    John Burkardt
  Parameters:
    Input, int M, N, the number of rows and columns of pixels.
    Input, int GRAY[M*N], the noisy grayscale data.
    Output, int GRAY_MEDIAN_NEWS[M*N], the grayscale data for the filtered image.
*/
{
	const _2dtpi * const s_data = data;
	const register dim_typ m = s_data->a0;
	const register dim_typ n = s_data->a1;
	int * gray = s_data->a2;
	
    int *gray2;
    dim_typ i, j;
    int p[5];

    gray2 = ( int * ) malloc ( m * n * sizeof ( int ) );
    /*
    Process the main part of the image:
    */
    for ( i = 1; i < m - 1; ++i)
    {
        for ( j = 1; j < n - 1; ++j )
        {
            p[0] = gray[i-1+ j   *m];
            p[1] = gray[i+1+ j   *m];
            p[2] = gray[i  +(j+1)*m];
            p[3] = gray[i  +(j-1)*m];
            p[4] = gray[i  + j   *m];

            gray2[i+j*m] = i4vec_median ( 5, p );
        }
    }
    /*
    Process the four borders.
    Get an odd number of data points,
    */
    for ( i = 1; i < m - 1; ++i )
    {
        j = 0;
        p[0] = gray[i-1+ j   *m];
        p[1] = gray[i+1+ j   *m];
        p[2] = gray[i  + j   *m];
        p[3] = gray[i  +(j+1)*m];
        p[4] = gray[i  +(j+2)*m];
        gray2[i+j*m] = i4vec_median ( 5, p );

        j = n - 1;
        p[0] = gray[i-1+ j   *m];
        p[1] = gray[i+1+ j   *m];
        p[2] = gray[i  +(j-2)*m];
        p[3] = gray[i  +(j-1)*m];
        p[4] = gray[i  + j   *m];
        gray2[i+j*m] = i4vec_median ( 5, p );
    }

    for ( j = 1; j < n - 1; ++j )
    {
        i = 0;
        p[0] = gray[i  + j   *m];
        p[1] = gray[i+1+ j   *m];
        p[2] = gray[i+2+ j   *m];
        p[3] = gray[i  +(j-1)*m];
        p[4] = gray[i  +(j+1)*m];
        gray2[i+j*m] = i4vec_median ( 5, p );

        i = m - 1;
        p[0] = gray[i-2+ j   *m];
        p[1] = gray[i-1+ j   *m];
        p[2] = gray[i  + j   *m];
        p[3] = gray[i  +(j-1)*m];
        p[4] = gray[i  +(j+1)*m];
        gray2[i+j*m] = i4vec_median ( 5, p );
    }
    /*
    Process the four corners.
    */
    i = j = 0;
    p[0] = gray[i+1+ j   *m];
    p[1] = gray[i  + j   *m];
    p[2] = gray[i  +(j+1)*m];
    gray2[i+j*m] = i4vec_median ( 3, p );

    i = 0;
    j = n - 1;
    p[0] = gray[i+1+ j   *m];
    p[1] = gray[i  + j   *m];
    p[2] = gray[i  +(j-1)*m];
    gray2[i+j*m] = i4vec_median ( 3, p );

    i = m - 1;
    j = 0;
    p[0] = gray[i-1+ j   *m];
    p[1] = gray[i  + j   *m];
    p[2] = gray[i  +(j+1)*m];
    gray2[i+j*m] = i4vec_median ( 3, p );

    i = m - 1;
    j = n - 1;
    p[0] = gray[i-1+ j   *m];
    p[1] = gray[i  + j   *m];
    p[2] = gray[i  +(j-1)*m];
    gray2[i+j*m] = i4vec_median ( 3, p );
    return gray2;
}

#endif
