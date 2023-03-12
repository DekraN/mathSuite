#ifndef __DISABLEDEEP_TRIANGLENCORULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_degree ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_DEGREE returns the degree of an NCO rule for the triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int TRIANGLE_NCO_DEGREE, the polynomial degree of exactness of
    the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
	result = 1 <= rule && rule <= 9 ? rule-1 : USHRT_MAX;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_order_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_ORDER_NUM returns the order of an NCO rule for the triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int TRIANGLE_NCO_ORDER_NUM, the order (number of points)
    of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
    dim_typ order;
    dim_typ order_num;
    int *suborder;
    dim_typ suborder_num;

    suborder_num = triangle_nco_suborder_num ( rule );

    suborder = triangle_nco_suborder ( rule, suborder_num );

    for ( order = order_num = 0; order < suborder_num; ++order )
        order_num += suborder[order];
    free ( suborder );

	result = order_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_RULE returns the points and weights of an NCO rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Input, int ORDER_NUM, the order (number of points) of the rule.
    Output, double XY[2*ORDER_NUM], the points of the rule.
    Output, double W[ORDER_NUM], the weights of the rule.
*/
{
	const _2dt2pit * const s_data = data; 
	const register dim_typ rule = s_data->a0;
	const register dim_typ order_num = s_data->a1;
	ityp * xy = s_data->a2;
	ityp * w = s_data->a3;
	
    dim_typ k;
    dim_typ o;
    dim_typ s;
    int *suborder;
    dim_typ suborder_num;
    ityp *suborder_w;
    ityp *suborder_xyz;
    /*
    Get the suborder information.
    */
    suborder_num = triangle_nco_suborder_num ( rule );

    suborder_xyz = ( ityp * ) malloc ( 3 * suborder_num * sizeof ( ityp ) );
    suborder_w = ( ityp * ) malloc ( suborder_num * sizeof ( ityp ) );

    suborder = triangle_nco_suborder ( rule, suborder_num );

    triangle_nco_subrule ( rule, suborder_num, suborder_xyz, suborder_w );
    /*
    Expand the suborder information to a full order rule.
    */
    o = 0;


    for ( s = 0; s < suborder_num; ++s )
    {
        if ( suborder[s] == 1 )
        {
            xy[0+(o<<1)] = suborder_xyz[0+s*3];
            xy[1+(o<<1)] = suborder_xyz[1+s*3];
            w[o] = suborder_w[s];
            ++ o;
        }
        else if ( suborder[s] == 3 )
        {
            #pragma omp parallel for num_threads(3)
            for ( k = 0; k < 3; ++k )
            {
                xy[0+(o<<1)] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
                xy[1+(o<<1)] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
                w[o] = suborder_w[s];
                ++ o;
            }
        }
        else if ( suborder[s] == 6 )
        {
            #pragma omp parallel for num_threads(3)
            for ( k = 0; k < 3; ++k )
            {
                xy[0+(o<<1)] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
                xy[1+(o<<1)] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
                w[o] = suborder_w[s];
                ++ o;
            }

            #pragma omp parallel for num_threads(3)
            for ( k = 0; k < 3; ++k )
            {
                xy[0+(o<<1)] = suborder_xyz [ i4_wrap(k+1,0,2) + s*3 ];
                xy[1+(o<<1)] = suborder_xyz [ i4_wrap(k,  0,2) + s*3 ];
                w[o] = suborder_w[s];
                ++ o;
            }
        }
        else
            return NULL;
    }

    free ( suborder );
    free ( suborder_xyz );
    free ( suborder_w );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_rule_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_RULE_NUM returns the number of NCO rules available.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Output, int TRIANGLE_NCO_RULE_NUM, the number of rules available.
*/
{
	static dim_typ result = 9;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _triangle_nco_suborder ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBORDER returns the suborders for an NCO rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int TRIANGLE_NCO_SUBORDER[SUBORDER_NUM],
    the suborders of the rule.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ rule = a_data[0];
	const register dim_typ suborder_num = a_data[1];
	
    int *suborder;

    suborder = ( int * ) malloc ( suborder_num * sizeof ( int ) );

    switch(rule)
    {
        case 1:
            suborder[0] = 1;
            break;
        case 2:
            suborder[0] = 3;
            break;
        case 3:
            suborder[0] = suborder[1] = 3;
            break;
        case 4:
            suborder[0] = 3;
            suborder[1] = 6;
            suborder[2] = 1;
            break;
        case 5:
            suborder[0] = suborder[2] = suborder[3] = 3;
            suborder[1] = 6;
            break;
        case 6:
            suborder[0] = suborder[3] = suborder[4] = 3;
            suborder[1] = suborder[2] = 6;
            break;
        case 7:
            suborder[0] = suborder[3] = suborder[4] = 3;
            suborder[1] = suborder[2] = suborder[5] = 6;
            suborder[6] = 1;
            break;
        case 8:
            suborder[0] = suborder[3] = suborder[6] = suborder[7] = 3;
            suborder[1] = suborder[2] = suborder[4] = suborder[5] = 6;
            break;
        case 9:
            suborder[0] = suborder[3] = suborder[6] = suborder[8] = suborder[9] = 3;
            suborder[1] = suborder[2] = suborder[4] = suborder[5] = suborder[7] =6;
            break;
        default:
            return NULL;
    }

    return suborder;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_suborder_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBORDER_NUM returns the number of suborders for an NCO rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int TRIANGLE_NCO_SUBORDER_NUM, the number of suborders
    of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
    dim_typ suborder_num;

    switch(rule)
    {
        case 1:
            suborder_num = 1;
            break;
        case 2:
            suborder_num = 1;
            break;
        case 3:
            suborder_num = 2;
            break;
        case 4:
            suborder_num = 3;
            break;
        case 5:
            suborder_num = 4;
            break;
        case 6:
            suborder_num = 5;
            break;
        case 7:
            suborder_num = 7;
            break;
        case 8:
            suborder_num = 8;
            break;
        case 9:
            suborder_num = 10;
            break;
        default:
            suborder_num = USHRT_MAX;
    }

	result = suborder_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _triangle_nco_subrule ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE returns a compressed NCO rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const _2dt2pit * const s_data = data; 
	const register dim_typ rule = s_data->a0;
	const register dim_typ suborder_num = s_data->a1;
	ityp * suborder_xyz = s_data->a2;
	ityp * suborder_w = s_data->a3;
	
    int i;
    int s;
    int suborder_w_d;
    int *suborder_w_n;
    int suborder_xyz_d;
    int *suborder_xyz_n;

    suborder_xyz_n = ( int * ) malloc ( 3 * suborder_num * sizeof ( int ) );
    suborder_w_n = ( int * ) malloc ( suborder_num * sizeof ( int ) );

    switch(rule)
    {
        case 1:
            triangle_nco_subrule_01 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 2:
            triangle_nco_subrule_02 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 3:
            triangle_nco_subrule_03 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 4:
            triangle_nco_subrule_04 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 5:
            triangle_nco_subrule_05 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 6:
            triangle_nco_subrule_06 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 7:
            triangle_nco_subrule_07 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 8:
            triangle_nco_subrule_08 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 9:
            triangle_nco_subrule_09 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        default:
            return NULL;
    }

    for ( s = 0; s < suborder_num; ++s)
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            suborder_xyz[i+s*3] =( ityp ) ( 1 + suborder_xyz_n[i+s*3] )/ ( ityp ) ( 3 + suborder_xyz_d );
        suborder_w[s] = ( ityp ) suborder_w_n[s] / ( ityp ) suborder_w_d;
    }
    for ( s = 0; s < suborder_num; ++s)
        suborder_w[s] = ( ityp ) suborder_w_n[s] / ( ityp ) suborder_w_d;

    free ( suborder_w_n );
    free ( suborder_xyz_n );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_01 returns a compressed NCO rule 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_01[3] =
    {
        0, 0, 0
    };
    const register dim_typ suborder_xyz_d_01 = 0;
    int suborder_w_n_01[1] =
    {
         1
        };
    const register dim_typ suborder_w_d_01 = 1;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i)
            suborder_xyz_n[i+s*3] = suborder_xyz_n_01[i+s*3];
        suborder_w_n[s] = suborder_w_n_01[s];
    }

    *suborder_xyz_d = suborder_xyz_d_01;
    *suborder_w_d = suborder_w_d_01;

    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_02 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_02 returns a compressed NCO rule 2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_02[3] =
    {
        1, 0, 0
    };
    const register dim_typ suborder_xyz_d_02 = 1;
    int suborder_w_n_02[1] =
    {
         1
    };
    const register dim_typ suborder_w_d_02 = 3;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            suborder_xyz_n[i+s*3] = suborder_xyz_n_02[i+s*3];
        suborder_w_n[s] = suborder_w_n_02[s];
    }

    *suborder_xyz_d = suborder_xyz_d_02;
    *suborder_w_d = suborder_w_d_02;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_03 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_03 returns a compressed NCO rule 3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{	
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_03[36] =
    {
        2, 0, 0,
        1, 1, 0
    };
    const register dim_typ suborder_xyz_d_03 = 2;
    int suborder_w_n_03[2] =
    {
         7, -3
        };
    const register dim_typ suborder_w_d_03 = 12;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            suborder_xyz_n[i+s*3] = suborder_xyz_n_03[i+s*3];
        suborder_w_n[s] = suborder_w_n_03[s];
    }

    *suborder_xyz_d = suborder_xyz_d_03;
    *suborder_w_d = suborder_w_d_03;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_04 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_04 returns a compressed NCO rule 4.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_04[9] =
    {
        3, 0, 0,
        2, 1, 0,
        1, 1, 1
    };
    const register dim_typ suborder_xyz_d_04 = 3;
    int suborder_w_n_04[3] =
    {
        8, 3, -12
    };
    const register dim_typ suborder_w_d_04 = 30;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            suborder_xyz_n[i+s*3] = suborder_xyz_n_04[i+s*3];
        suborder_w_n[s] = suborder_w_n_04[s];
    }

    *suborder_xyz_d = suborder_xyz_d_04;
    *suborder_w_d = suborder_w_d_04;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_05 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_05 returns a compressed NCO rule 5.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_05[12] =
    {
        4, 0, 0,
        3, 1, 0,
        2, 2, 0,
        2, 1, 1
    };
    const register dim_typ suborder_xyz_d_05 = 4;
    int suborder_w_n_05[4] =
    {
        307, -316, 629, -64
    };
    const register dim_typ suborder_w_d_05 = 720;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i)
            suborder_xyz_n[i+s*3] = suborder_xyz_n_05[i+s*3];
        suborder_w_n[s] = suborder_w_n_05[s];
    }

    *suborder_xyz_d = suborder_xyz_d_05;
    *suborder_w_d = suborder_w_d_05;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_06 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_06 returns a compressed NCO rule 6.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_06[15] =
    {
        5, 0, 0,
        4, 1, 0,
        3, 2, 0,
        3, 1, 1,
        2, 2, 1
    };
    const register dim_typ suborder_xyz_d_06 = 5;
    int suborder_w_n_06[5] =
    {
        71, -13, 57, -167, 113
    };
    const register dim_typ suborder_w_d_06 = 315;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            suborder_xyz_n[i+s*3] = suborder_xyz_n_06[i+s*3];
        suborder_w_n[s] = suborder_w_n_06[s];
    }

    *suborder_xyz_d = suborder_xyz_d_06;
    *suborder_w_d = suborder_w_d_06;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_07 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_07 returns a compressed NCO rule 7.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_07[21] =
    {
        6, 0, 0,
        5, 1, 0,
        4, 2, 0,
        4, 1, 1,
        3, 3, 0,
        3, 2, 1,
        2, 2, 2
    };
    const register dim_typ suborder_xyz_d_07 = 6;
    int suborder_w_n_07[7] =
    {
        767, -1257, 2901, 387, -3035, -915, 3509
    };
    const register dim_typ suborder_w_d_07 = 2240;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            suborder_xyz_n[i+s*3] = suborder_xyz_n_07[i+s*3];
        suborder_w_n[s] = suborder_w_n_07[s];
    }

    *suborder_xyz_d = suborder_xyz_d_07;
    *suborder_w_d = suborder_w_d_07;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_08 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_08 returns a compressed NCO rule 8.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_08[24] =
    {
        7, 0, 0,
        6, 1, 0,
        5, 2, 0,
        5, 1, 1,
        4, 3, 0,
        4, 2, 1,
        3, 3, 1,
        3, 2, 2
    };
    const register dim_typ suborder_xyz_d_08 = 7;
    int suborder_w_n_08[8] =
    {
        898, -662, 1573, -2522, -191, 2989, -5726, 1444
    };
    const register dim_typ suborder_w_d_08 = 4536;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            suborder_xyz_n[i+s*3] = suborder_xyz_n_08[i+s*3];
        suborder_w_n[s] = suborder_w_n_08[s];
    }

    *suborder_xyz_d = suborder_xyz_d_08;
    *suborder_w_d = suborder_w_d_08;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_nco_subrule_09 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_NCO_SUBRULE_09 returns a compressed NCO rule 9.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[3*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.
    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.
    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.
    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
	const dt4pi * const s_data = data;
	const register dim_typ suborder_num = s_data->a0;
	int * suborder_xyz_n = s_data->a1;
	int * suborder_xyz_d = s_data->a2;
	int * suborder_w_n = s_data->a3;
	int * suborder_w_d = s_data->a4;
	
    dim_typ i;
    dim_typ s;
    int suborder_xyz_n_09[30] =
    {
        8, 0, 0,
        7, 1, 0,
        6, 2, 0,
        6, 1, 1,
        5, 3, 0,
        5, 2, 1,
        4, 4, 0,
        4, 3, 1,
        4, 2, 2,
        3, 3, 2
    };
    const register dim_typ suborder_xyz_d_09 = 8;
    int suborder_w_n_09[10] =
    {
        1051445, -2366706, 6493915, 1818134, -9986439,-3757007, 12368047,
        478257, 10685542, -6437608
    };
    const register unsigned suborder_w_d_09 = 3628800;

    for ( s = 0; s < suborder_num; ++s )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            suborder_xyz_n[i+s*3] = suborder_xyz_n_09[i+s*3];
        suborder_w_n[s] = suborder_w_n_09[s];
    }

    *suborder_xyz_d = suborder_xyz_d_09;
    *suborder_w_d = suborder_w_d_09;

    return NULL;
}

#endif
