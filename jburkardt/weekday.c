#ifndef __DISABLEDEEP_WEEKDAY

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _weekday_values ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEEKDAY_VALUES returns the day of the week for various dates.
  Discussion:
    The CE or Common Era calendar is used, under the
    hybrid Julian/Gregorian Calendar, with a transition from Julian
    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
    years BC/BCE are indicated by a negative year value.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 May 2012
  Author:
    John Burkardt
  Reference:
    Edward Reingold, Nachum Dershowitz,
    Calendrical Calculations: The Millennium Edition,
    Cambridge University Press, 2001,
    ISBN: 0 521 77752 6
    LC: CE12.r85.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0
    before the first call.  On each call, the routine increments N_DATA by 1,
    and returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, int *Y, *M, *D, the Common Era date.
    Output, int *W, the day of the week.  Sunday = 1.
*/
{
	dim_typ ** const a_data = data;
	dim_typ * n_data = a_data[0];
	dim_typ * y = a_data[1];
	dim_typ * m = a_data[2];
	dim_typ * d = a_data[3];
	dim_typ * w = a_data[4];
	
    # define N_MAX 34

    static dim_typ d_vec[N_MAX] =
    {
        30,
        8,
        26,
        3,
        7,
        18,
        7,
        19,
        14,
        18,
        16,
        3,
        26,
        20,
        4,
        25,
        31,
        9,
        24,
        10,
        30,
        24,
        19,
        2,
        27,
        19,
        25,
        29,
        19,
        7,
        17,
        25,
        10,
        18
    };
    static dim_typ m_vec[N_MAX] =
    {
        7,
        12,
        9,
        10,
        1,
        5,
        11,
        4,
        10,
        5,
        3,
        3,
        3,
        4,
        6,
        1,
        3,
        9,
        2,
        6,
        6,
        7,
        6,
        8,
        3,
        4,
        8,
        9,
        4,
        10,
        3,
        2,
        11,
        7
    };
    static dim_typ w_vec[N_MAX] =
    {
        1,
        4,
        4,
        1,
        4,
        2,
        7,
        1,
        7,
        1,
        6,
        7,
        6,
        1,
        1,
        4,
        7,
        7,
        7,
        4,
        1,
        6,
        1,
        2,
        4,
        1,
        1,
        2,
        2,
        5,
        3,
        1,
        4,
        1
    };
    static dim_typ y_vec[N_MAX] =
    {
        - 587,
        - 169,
        70,
        135,
        470,
        576,
        694,
        1013,
        1066,
        1096,
        1190,
        1240,
        1288,
        1298,
        1391,
        1436,
        1492,
        1553,
        1560,
        1648,
        1680,
        1716,
        1768,
        1819,
        1839,
        1903,
        1929,
        1941,
        1943,
        1943,
        1992,
        1996,
        2038,
        2094
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
        *n_data = *y = *m = *d = *w = 0;
    else
    {
        *y = y_vec[*n_data-1];
        *m = m_vec[*n_data-1];
        *d = d_vec[*n_data-1];
        *w = w_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
