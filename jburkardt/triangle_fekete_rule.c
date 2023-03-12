#ifndef __DISABLEDEEP_TRIANGLEFEKETERULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_degree ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_DEGREE returns the degree of a Fekete rule for the triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int FEKETE_DEGREE, the polynomial degree of exactness of
    the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
    dim_typ degree = USHRT_MAX;

    if ( rule == 1 )
        degree = 3;
    else if ( rule == 2 )
        degree = 6;
    else if ( rule == 3 )
        degree = 9;
    else if ( rule == 4 )
        degree = 12;
    else if ( rule == 5 )
        degree = 12;
    else if ( rule == 6 )
        degree = 15;
    else if ( rule == 7 )
        degree = 18;

	result = degree;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_order_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_ORDER_NUM returns the order of a Fekete rule for the triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int FEKETE_ORDER_NUM, the order (number of points) of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
    dim_typ order;
    dim_typ order_num;
    int *suborder;
    dim_typ suborder_num;

    suborder_num = fekete_suborder_num ( rule );
    suborder = fekete_suborder ( rule, suborder_num );

    for ( order = order_num = 0; order < suborder_num; ++order )
        order_num += suborder[order];

    free ( suborder );

	result = order_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _fekete_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_RULE returns the points and weights of a Fekete rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
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
	ityp *xy = s_data->a2;
	ityp *w = s_data->a3;
	
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
    suborder_num = fekete_suborder_num ( rule );

    suborder_xyz = ( ityp * ) malloc ( 3 * suborder_num * sizeof ( ityp ) );
    suborder_w = ( ityp * ) malloc ( suborder_num * sizeof ( ityp ) );

    suborder = fekete_suborder ( rule, suborder_num );

    fekete_subrule ( rule, suborder_num, suborder_xyz, suborder_w );
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
            for ( k = 0; k < 3; ++k)
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
            for ( k = 0; k < 3; ++k)
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
__MATHSUITE __JBURKARDT  void *   _fekete_rule_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_RULE_NUM returns the number of Fekete rules available.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
    John Burkardt
  Reference:
  Parameters:
    Output, int FEKETE_RULE_NUM, the number of rules available.
*/
{
	static dim_typ result = 7;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _fekete_suborder ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBORDER returns the suborders for a Fekete rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int RULE, the index of the rule.
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int FEKETE_SUBORDER[SUBORDER_NUM], the suborders of the rule.
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
            suborder[1] = 3;
            suborder[2] = 6;
            break;
        case 2:
            suborder[0] = 1;
            suborder[1] = 3;
            suborder[2] = 3;
            suborder[3] = 3;
            suborder[4] = 6;
            suborder[5] = 6;
            suborder[6] = 6;
            break;
        case 3:
            suborder[0] = 1;
            suborder[1] = 3;
            suborder[2] = 3;
            suborder[3] = 3;
            suborder[4] = 3;
            suborder[5] = 6;
            suborder[6] = 6;
            suborder[7] = 6;
            suborder[8] = 6;
            suborder[9] = 6;
            suborder[10] = 6;
            suborder[11] = 6;
            break;
        case 4:
            suborder[0] = 1;
            suborder[1] = 3;
            suborder[2] = 3;
            suborder[3] = 3;
            suborder[4] = 3;
            suborder[5] = 3;
            suborder[6] = 3;
            suborder[7] = 6;
            suborder[8] = 6;
            suborder[9] = 6;
            suborder[10] = 6;
            suborder[11] = 6;
            suborder[12] = 6;
            suborder[13] = 6;
            suborder[14] = 6;
            suborder[15] = 6;
            suborder[16] = 6;
            suborder[17] = 6;
            suborder[18] = 6;
            break;
        case 5:
            suborder[0] = 1;
            suborder[1] = 3;
            suborder[2] = 3;
            suborder[3] = 3;
            suborder[4] = 3;
            suborder[5] = 3;
            suborder[6] = 3;
            suborder[7] = 3;
            suborder[8] = 3;
            suborder[9] = 3;
            suborder[10] = 3;
            suborder[11] = 6;
            suborder[12] = 6;
            suborder[13] = 6;
            suborder[14] = 6;
            suborder[15] = 6;
            suborder[16] = 6;
            suborder[17] = 6;
            suborder[18] = 6;
            suborder[19] = 6;
            suborder[20] = 6;
            break;
        case 6:
            suborder[0] = 1;
            suborder[1] = 3;
            suborder[2] = 3;
            suborder[3] = 3;
            suborder[4] = 3;
            suborder[5] = 3;
            suborder[6] = 3;
            suborder[7] = 3;
            suborder[8] = 3;
            suborder[9] = 3;
            suborder[10] = 6;
            suborder[11] = 6;
            suborder[12] = 6;
            suborder[13] = 6;
            suborder[14] = 6;
            suborder[15] = 6;
            suborder[16] = 6;
            suborder[17] = 6;
            suborder[18] = 6;
            suborder[19] = 6;
            suborder[20] = 6;
            suborder[21] = 6;
            suborder[22] = 6;
            suborder[23] = 6;
            suborder[24] = 6;
            suborder[25] = 6;
            suborder[26] = 6;
            suborder[27] = 6;
            break;
        case 7:
            suborder[0] = 1;
            suborder[1] = 3;
            suborder[2] = 3;
            suborder[3] = 3;
            suborder[4] = 3;
            suborder[5] = 3;
            suborder[6] = 3;
            suborder[7] = 3;
            suborder[8] = 3;
            suborder[9] = 3;
            suborder[10] = 3;
            suborder[11] = 3;
            suborder[12] = 6;
            suborder[13] = 6;
            suborder[14] = 6;
            suborder[15] = 6;
            suborder[16] = 6;
            suborder[17] = 6;
            suborder[18] = 6;
            suborder[19] = 6;
            suborder[20] = 6;
            suborder[21] = 6;
            suborder[22] = 6;
            suborder[23] = 6;
            suborder[24] = 6;
            suborder[25] = 6;
            suborder[26] = 6;
            suborder[27] = 6;
            suborder[28] = 6;
            suborder[29] = 6;
            suborder[30] = 6;
            suborder[31] = 6;
            suborder[32] = 6;
            suborder[33] = 6;
            suborder[34] = 6;
            suborder[35] = 6;
            suborder[36] = 6;
            suborder[37] = 6;
            break;
        default:
            return NULL;
    }

    return suborder;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_suborder_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBORDER_NUM returns the number of suborders for a Fekete rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 October 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int FEKETE_SUBORDER_NUM, the number of suborders of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
    dim_typ suborder_num;

    switch(rule)
    {
        case 1:
            suborder_num = 3;
            break;
        case 2:
            suborder_num = 7;
            break;
        case 3:
            suborder_num = 12;
            break;
        case 4:
            suborder_num = 19;
            break;
        case 5:
            suborder_num = 21;
            break;
        case 6:
            suborder_num = 28;
            break;
        case 7:
            suborder_num = 38;
            break;
        default:
            suborder_num = USHRT_MAX;
    }

	result = suborder_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_subrule ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBRULE returns a compressed Fekete rule.
  Discussion:
    The listed weights are twice what we want...since we want them
    to sum to 1/2, reflecting the area of a unit triangle.  So we
    simple halve the values before exiting this routine.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 January 2007
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
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
	
    int s;

    switch(rule)
    {
        case 1:
            fekete_subrule_1 ( suborder_num, suborder_xyz, suborder_w );
            break;
        case 2:
            fekete_subrule_2 ( suborder_num, suborder_xyz, suborder_w );
            break;
        case 3:
            fekete_subrule_3 ( suborder_num, suborder_xyz, suborder_w );
            break;
        case 4:
            fekete_subrule_4 ( suborder_num, suborder_xyz, suborder_w );
            break;
        case 5:
            fekete_subrule_5 ( suborder_num, suborder_xyz, suborder_w );
            break;
        case 6:
            fekete_subrule_6 ( suborder_num, suborder_xyz, suborder_w );
            break;
        case 7:
            fekete_subrule_7 ( suborder_num, suborder_xyz, suborder_w );
            break;
        default:
            return NULL;
    }

    for ( s = 0; s < suborder_num; ++s )
        suborder_w[s] = 0.50 * suborder_w[s];

    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_subrule_1 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBRULE_1 returns a compressed Fekete rule 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 October 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const dt2pit * const s_data = data; 
	const register dim_typ suborder_num = s_data->a0;
	ityp * suborder_xyz = s_data->a1;
	ityp * suborder_w = s_data->a2;
	
    static ityp suborder_xy_rule_1[9] =
    {
        0.3333333333, 0.3333333333, 0.3333333334,
        1.0000000000, 0.0000000000, 0.0000000000,
        0.0000000000, 0.2763932023, 0.7236067977
    };
    static ityp suborder_w_rule_1[3] =
    {
        0.9000000000,
        0.0333333333,
        0.1666666667
    };

    r8mat_copy ( 3, suborder_num, suborder_xy_rule_1, suborder_xyz );
    r8vec_copy ( suborder_num, suborder_w_rule_1, suborder_w );

    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_subrule_2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBRULE_2 returns a compressed Fekete rule 2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 December 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const dt2pit * const s_data = data; 
	const register dim_typ suborder_num = s_data->a0;
	ityp * suborder_xyz = s_data->a1;
	ityp * suborder_w = s_data->a2;
	
    static ityp suborder_xy_rule_2[21] =
    {
        0.3333333333,  0.3333333333, 0.3333333334,
        0.1063354684,  0.1063354684, 0.7873290632,
        0.5000000000,  0.5000000000, 0.0000000000,
        1.0000000000,  0.0000000000, 0.0000000000,
        0.1171809171,  0.3162697959, 0.5665492870,
        0.0000000000,  0.2655651402, 0.7344348598,
        0.0000000000,  0.0848854223, 0.9151145777
    };
    static ityp suborder_w_rule_2[7] =
    {
        0.2178563571,
        0.1104193374,
        0.0358939762,
        0.0004021278,
        0.1771348660,
        0.0272344079,
        0.0192969460
    };

    r8mat_copy ( 3, suborder_num, suborder_xy_rule_2, suborder_xyz );
    r8vec_copy ( suborder_num, suborder_w_rule_2, suborder_w );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_subrule_3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBRULE_3 returns a compressed Fekete rule 3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 December 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const dt2pit * const s_data = data; 
	const register dim_typ suborder_num = s_data->a0;
	ityp * suborder_xyz = s_data->a1;
	ityp * suborder_w = s_data->a2;
	
    static ityp suborder_xy_rule_3[36] =
    {
        0.3333333333,  0.3333333333, 0.3333333334,
        0.1704318201,  0.1704318201, 0.6591363598,
        0.0600824712,  0.4699587644, 0.4699587644,
        0.0489345696,  0.0489345696, 0.9021308608,
        0.0000000000,  0.0000000000, 1.0000000000,
        0.1784337588,  0.3252434900, 0.4963227512,
        0.0588564879,  0.3010242110, 0.6401193011,
        0.0551758079,  0.1543901944, 0.7904339977,
        0.0000000000,  0.4173602935, 0.5826397065,
        0.0000000000,  0.2610371960, 0.7389628040,
        0.0000000000,  0.1306129092, 0.8693870908,
        0.0000000000,  0.0402330070, 0.9597669930
    };
    static ityp suborder_w_rule_3[12] =
    {
        0.1096011288,
        0.0767491008,
        0.0646677819,
        0.0276211659,
        0.0013925011,
        0.0933486453,
        0.0619010169,
        0.0437466450,
        0.0114553907,
        0.0093115568,
        0.0078421987,
        0.0022457501
    };

    for ( dim_typ s = 0; s < suborder_num; ++s )
    {
        suborder_xyz[0+s*3] = suborder_xy_rule_3[0+s*3];
        suborder_xyz[1+s*3] = suborder_xy_rule_3[1+s*3];
        suborder_xyz[2+s*3] = suborder_xy_rule_3[2+s*3];
        suborder_w[s] = suborder_w_rule_3[s];
    }

    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_subrule_4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBRULE_4 returns a compressed Fekete rule 4.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 December 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const dt2pit * const s_data = data; 
	const register dim_typ suborder_num = s_data->a0;
	ityp * suborder_xyz = s_data->a1;
	ityp * suborder_w = s_data->a2;
	
    static ityp suborder_xy_rule_4[57] =
    {
        0.3333333333,  0.3333333333, 0.3333333334,
        0.1988883477,  0.4005558262, 0.4005558261,
        0.2618405201,  0.2618405201, 0.4763189598,
        0.0807386775,  0.0807386775, 0.8385226450,
        0.0336975736,  0.0336975736, 0.9326048528,
        0.0000000000,  0.5000000000, 0.5000000000,
        0.0000000000,  0.0000000000, 1.0000000000,
        0.1089969290,  0.3837518758, 0.5072511952,
        0.1590834479,  0.2454317980, 0.5954847541,
        0.0887037176,  0.1697134458, 0.7415828366,
        0.0302317829,  0.4071849276, 0.5625832895,
        0.0748751152,  0.2874821712, 0.6376427136,
        0.0250122615,  0.2489279690, 0.7260597695,
        0.0262645218,  0.1206826354, 0.8530528428,
        0.0000000000,  0.3753565349, 0.6246434651,
        0.0000000000,  0.2585450895, 0.7414549105,
        0.0000000000,  0.1569057655, 0.8430942345,
        0.0000000000,  0.0768262177, 0.9231737823,
        0.0000000000,  0.0233450767, 0.9766549233
    };
    static double suborder_w_rule_4[19] =
    {
        0.0626245179,
        0.0571359417,
        0.0545982307,
        0.0172630326,
        0.0142519606,
        0.0030868485,
        0.0004270742,
        0.0455876390,
        0.0496701966,
        0.0387998322,
        0.0335323983,
        0.0268431561,
        0.0237377452,
        0.0177255972,
        0.0043097313,
        0.0028258057,
        0.0030994935,
        0.0023829062,
        0.0009998683
    };

    for (dim_typ s = 0; s < suborder_num; ++s )
    {
        suborder_xyz[0+s*3] = suborder_xy_rule_4[0+s*3];
        suborder_xyz[1+s*3] = suborder_xy_rule_4[1+s*3];
        suborder_xyz[2+s*3] = suborder_xy_rule_4[2+s*3];
        suborder_w[s] = suborder_w_rule_4[s];
    }

    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_subrule_5 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBRULE_5 returns a compressed Fekete rule 5.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 December 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const dt2pit * const s_data = data; 
	const register dim_typ suborder_num = s_data->a0;
	ityp * suborder_xyz = s_data->a1;
	ityp * suborder_w = s_data->a2;
	
    static ityp suborder_xy_rule_5[63] =
    {
        0.3333333333,  0.3333333333, 0.3333333334,
        0.2201371125,  0.3169406831, 0.4629222044,
        0.2201371125,  0.4629222044, 0.3169406831,
        0.1877171129,  0.1877171129, 0.6245657742,
        0.1403402144,  0.4298298928, 0.4298298928,
        0.0833252778,  0.0833252778, 0.8333494444,
        0.0664674598,  0.0252297247, 0.9083028155,
        0.0218884020,  0.4890557990, 0.4890557990,
        0.0252297247,  0.0664674598, 0.9083028155,
        0.0000000000,  0.5000000000, 0.5000000000,
        0.0000000000,  0.0000000000, 1.0000000000,
        0.1157463404,  0.2842319093, 0.6000217503,
        0.0672850606,  0.3971764400, 0.5355384994,
        0.0909839531,  0.1779000668, 0.7311159801,
        0.0318311633,  0.3025963402, 0.6655724965,
        0.0273518579,  0.1733665506, 0.7992815915,
        0.0000000000,  0.3753565349, 0.6246434651,
        0.0000000000,  0.2585450895, 0.7414549105,
        0.0000000000,  0.1569057655, 0.8430942345,
        0.0000000000,  0.0768262177, 0.9231737823,
        0.0000000000,  0.0233450767, 0.9766549233
    };
    static ityp suborder_w_rule_5[21] =
    {
        0.0485965670,
        0.0602711576,
        0.0602711576,
        0.0476929767,
        0.0453940802,
        0.0258019417,
        0.0122004614,
        0.0230003812,
        0.0122004614,
        0.0018106475,
        -0.0006601747,
        0.0455413513,
        0.0334182802,
        0.0324896773,
        0.0299402736,
        0.0233477738,
        0.0065962854,
        0.0021485117,
        0.0034785755,
        0.0013990566,
        0.0028825748
    };

    for ( dim_typ s = 0; s < suborder_num; ++s )
    {
        suborder_xyz[0+s*3] = suborder_xy_rule_5[0+s*3];
        suborder_xyz[1+s*3] = suborder_xy_rule_5[1+s*3];
        suborder_xyz[2+s*3] = suborder_xy_rule_5[2+s*3];
        suborder_w[s] = suborder_w_rule_5[s];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_subrule_6 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBRULE_6 returns a compressed Fekete rule 6.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 December 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const dt2pit * const s_data = data; 
	const register dim_typ suborder_num = s_data->a0;
	ityp * suborder_xyz = s_data->a1;
	ityp * suborder_w = s_data->a2;
	
    static ityp suborder_xy_rule_6[84] =
    {
        0.3333333333,  0.3333333333, 0.3333333334,
        0.2379370518,  0.3270403780, 0.4350225702,
        0.3270403780,  0.2379370518, 0.4350225702,
        0.1586078048,  0.4206960976, 0.4206960976,
        0.2260541354,  0.2260541354, 0.5478917292,
        0.1186657611,  0.1186657611, 0.7626684778,
        0.0477095725,  0.4761452137, 0.4761452138,
        0.0531173538,  0.0531173538, 0.8937652924,
        0.0219495841,  0.0219495841, 0.9561008318,
        0.0000000000,  0.0000000000, 1.0000000000,
        0.1585345951,  0.3013819154, 0.5400834895,
        0.0972525649,  0.3853507643, 0.5173966708,
        0.0875150140,  0.2749910734, 0.6374939126,
        0.1339547708,  0.1975591066, 0.6684861226,
        0.0475622627,  0.3524012205, 0.6000365168,
        0.0596194677,  0.1978887556, 0.7424917767,
        0.0534939782,  0.1162464503, 0.8302595715,
        0.0157189888,  0.4176001732, 0.5666808380,
        0.0196887324,  0.2844332752, 0.6958779924,
        0.0180698489,  0.1759511193, 0.8059790318,
        0.0171941515,  0.0816639421, 0.9011419064,
        0.0000000000,  0.4493368632, 0.5506631368,
        0.0000000000,  0.3500847655, 0.6499152345,
        0.0000000000,  0.2569702891, 0.7430297109,
        0.0000000000,  0.1738056486, 0.8261943514,
        0.0000000000,  0.1039958541, 0.8960041459,
        0.0000000000,  0.0503997335, 0.9496002665,
        0.0000000000,  0.0152159769, 0.9847840231
    };
    static ityp suborder_w_rule_6[28] =
    {
        0.0459710878,
        0.0346650571,
        0.0346650571,
        0.0384470625,
        0.0386013566,
        0.0224308157,
        0.0243531004,
        0.0094392654,
        0.0061105652,
        0.0001283162,
        0.0305412307,
        0.0262101254,
        0.0265367617,
        0.0269859772,
        0.0172635676,
        0.0188795851,
        0.0158224870,
        0.0127170850,
        0.0164489660,
        0.0120018620,
        0.0072268907,
        0.0023599161,
        0.0017624674,
        0.0018648017,
        0.0012975716,
        0.0018506035,
        0.0009919379,
        0.0004893506
    };

    for ( dim_typ s = 0; s < suborder_num; ++s)
    {
        suborder_xyz[0+s*3] = suborder_xy_rule_6[0+s*3];
        suborder_xyz[1+s*3] = suborder_xy_rule_6[1+s*3];
        suborder_xyz[2+s*3] = suborder_xy_rule_6[2+s*3];
        suborder_w[s] = suborder_w_rule_6[s];
    }

    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fekete_subrule_7 ( void * data)
/******************************************************************************/
/*
  Purpose:
    FEKETE_SUBRULE_7 returns a compressed Fekete rule 7.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    05 October 2006
  Author:
    John Burkardt
  Reference:
    Mark Taylor, Beth Wingate, Rachel Vincent,
    An Algorithm for Computing Fekete Points in the Triangle,
    SIAM Journal on Numerical Analysis,
    Volume 38, Number 5, 2000, pages 1707-1720.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM], contains
    Output, double SUBORDER_XYZ[3*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const dt2pit * const s_data = data; 
	const register dim_typ suborder_num = s_data->a0;
	ityp * suborder_xyz = s_data->a1;
	ityp * suborder_w = s_data->a2;
	
    static ityp suborder_xy_rule_7[114] =
    {
        0.3333333333,  0.3333333333, 0.3333333334,
        0.2515553103,  0.3292984162, 0.4191462735,
        0.3292984162,  0.2515553103, 0.4191462735,
        0.1801930996,  0.4099034502, 0.4099034502,
        0.2438647767,  0.2438647767, 0.5122704466,
        0.1512564554,  0.1512564554, 0.6974870892,
        0.0810689493,  0.4594655253, 0.4594655254,
        0.0832757649,  0.0832757649, 0.8334484702,
        0.0369065587,  0.0369065587, 0.9261868826,
        0.0149574850,  0.0149574850, 0.9700850300,
        0.0000000000,  0.5000000000, 0.5000000000,
        0.0000000000,  0.0000000000, 1.0000000000,
        0.1821465920,  0.3095465041, 0.5083069039,
        0.1246901255,  0.3789288931, 0.4963809814,
        0.1179441386,  0.2868915642, 0.5951642972,
        0.1639418454,  0.2204868669, 0.6155712877,
        0.0742549663,  0.3532533654, 0.5724916683,
        0.0937816771,  0.2191980979, 0.6870202250,
        0.0890951387,  0.1446273457, 0.7662775156,
        0.0409065243,  0.4360543636, 0.5230391121,
        0.0488675890,  0.2795984854, 0.6715339256,
        0.0460342127,  0.2034211147, 0.7505446726,
        0.0420687187,  0.1359040280, 0.8220272533,
        0.0116377940,  0.4336892286, 0.5546729774,
        0.0299062187,  0.3585587824, 0.6115349989,
        0.0132313129,  0.2968103667, 0.6899583204,
        0.0136098469,  0.2050279257, 0.7813622274,
        0.0124869684,  0.1232146223, 0.8642984093,
        0.0365197797,  0.0805854893, 0.8828947310,
        0.0118637765,  0.0554881302, 0.9326480933,
        0.0000000000,  0.4154069883, 0.5845930117,
        0.0000000000,  0.3332475761, 0.6667524239,
        0.0000000000,  0.2558853572, 0.7441146428,
        0.0000000000,  0.1855459314, 0.8144540686,
        0.0000000000,  0.1242528987, 0.8757471013,
        0.0000000000,  0.0737697111, 0.9262302889,
        0.0000000000,  0.0355492359, 0.9644507641,
        0.0000000000,  0.0106941169, 0.9893058831
    };
    static double suborder_w_rule_7[38] =
    {
        0.0326079297,
        0.0255331366,
        0.0255331366,
        0.0288093886,
        0.0279490452,
        0.0174438045,
        0.0203594338,
        0.0113349170,
        0.0046614185,
        0.0030346239,
        0.0012508731,
        0.0000782945,
        0.0235716330,
        0.0206304700,
        0.0204028340,
        0.0215105697,
        0.0183482070,
        0.0174161032,
        0.0155972434,
        0.0119269616,
        0.0147074804,
        0.0116182830,
        0.0087639138,
        0.0098563528,
        0.0096342355,
        0.0086477936,
        0.0083868302,
        0.0062576643,
        0.0077839825,
        0.0031415239,
        0.0006513246,
        0.0021137942,
        0.0004393452,
        0.0013662119,
        0.0003331251,
        0.0011613225,
        0.0004342867,
        0.0002031499
    };

    for ( dim_typ s = 0; s < suborder_num; ++s )
    {
        suborder_xyz[0+s*3] = suborder_xy_rule_7[0+s*3];
        suborder_xyz[1+s*3] = suborder_xy_rule_7[1+s*3];
        suborder_xyz[2+s*3] = suborder_xy_rule_7[2+s*3];
        suborder_w[s] = suborder_w_rule_7[s];
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _reference_to_physical_t3 ( void * data)
/******************************************************************************/
/*
  Purpose:
    REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
  Discussion:
    Given the vertices of an order 3 physical triangle and a point
 (XSI,ETA) in the reference triangle, the routine computes the value
    of the corresponding image point (X,Y) in physical space.
    Note that this routine may also be appropriate for an order 6
    triangle, if the mapping between reference and physical space
    is linear.  This implies, in particular, that the sides of the
    image triangle are straight and that the "midside" nodes in the
    physical triangle are literally halfway along the sides of
    the physical triangle.
  Reference Element T3:
    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 June 2005
  Author:
    John Burkardt
  Parameters:
    Input, double T[2*3], the coordinates of the vertices.
    The vertices are assumed to be the images of (0,0), (1,0) and
 (0,1) respectively.
    Input, int N, the number of objects to transform.
    Input, double REF[2*N], points in the reference triangle.
    Output, double PHY[2*N], corresponding points in the
    physical triangle.
*/
{
	
	const dt3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * ref = s_data->a1;
	ityp * phy = s_data->a2; 
	ityp * t = s_data->a3;
	
    dim_typ i, j;

    for ( j = 0; j < n; ++j )
        #pragma omp parallel for num_threads(2)
        for ( i = 0; i < 2; ++i )
            phy[i+(j<<1)] = t[i] * ( 1.0 - ref[(j<<1)] - ref[1+(j<<1)] )+ t[i+2] *       + ref[(j<<1)]+ t[i+4] *                    + ref[1+(j<<1)];

    return NULL;
}

#endif
