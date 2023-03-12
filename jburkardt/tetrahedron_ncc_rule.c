#ifndef __DISABLEDEEP_TETRAHEDRONNCCRULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _reference_to_physical_t4 ( void * data)
/******************************************************************************/
/*
  Purpose:
    REFERENCE_TO_PHYSICAL_T4 maps T4 reference points to physical points.
  Discussion:
    Given the vertices of an order 4 physical tetrahedron and a point
 (R,S,T) in the reference tetrahedron, the routine computes the value
    of the corresponding image point (X,Y,Z) in physical space.
    This routine will also be correct for an order 10 tetrahedron,
    if the mapping between reference and physical space
    is linear.  This implies, in particular, that the sides of the
    image tetrahedron are straight, the faces are flat, and
    the "midside" nodes in the physical tetrahedron are
    halfway along the edges of the physical tetrahedron.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Parameters:
    Input, double T[3*4], the coordinates of the vertices.
    The vertices are assumed to be the images of (0,0,0), (1,0,0),
 (0,1,0) and (0,0,1) respectively.
    Input, int N, the number of objects to transform.
    Input, double REF[3*N], points in the reference tetrahedron.
    Output, double PHY[3*N], corresponding points in the
    physical tetrahedron.
*/
{
	const dt3pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	ityp * ref = s_data->a1;
	ityp * phy = s_data->a2; 
	ityp * t = s_data->a3;
	
    dim_typ i, j;

    for ( j = 0; j < n; ++j )
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            phy[i+j*3] = t[i+0*3] * ( 1.0 - ref[0+j*3] - ref[1+j*3] - ref[2+j*3] )+ t[i+1*3] *       + ref[0+j*3]+ t[i+2*3] *                    + ref[1+j*3]+ t[i+3*3] *                                 + ref[2+j*3];

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_degree ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_DEGREE: degree of an NCC rule for the tetrahedron.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int TETRAHEDRON_NCC_DEGREE, the polynomial degree of exactness of
    the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
	result = 1 <= rule && rule <= 7 ? rule-1 : USHRT_MAX;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_order_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_ORDER_NUM: order of an NCC rule for the tetrahedron.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int TETRAHEDRON_NCC_ORDER_NUM, the order (number of points)
    of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
    dim_typ order;
    dim_typ order_num=0;
    int *suborder;
    dim_typ suborder_num = tetrahedron_ncc_suborder_num ( rule );
    suborder = tetrahedron_ncc_suborder ( rule, suborder_num );
    for ( order = order_num=0; order < suborder_num; ++order )
        order_num += suborder[order];

    free ( suborder );
    
    result = order_num;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tetrahedron_ncc_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_RULE returns the points and weights of an NCC rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
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
    Output, double XYZ[3*ORDER_NUM], the points of the rule.
    Output, double W[ORDER_NUM], the weights of the rule.
*/
{
	const dtit2pit * const s_data = data;
	const register dim_typ rule = s_data->a0;
	int order_num = s_data->a1;
	ityp * xyz = s_data->a2;
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
    suborder_num = tetrahedron_ncc_suborder_num ( rule );

    suborder_xyz = ( ityp * ) malloc ( suborder_num * sizeof ( ityp )<<2);
    suborder_w = ( ityp * ) malloc ( suborder_num * sizeof ( ityp ) );

    suborder = tetrahedron_ncc_suborder ( rule, suborder_num );

    tetrahedron_ncc_subrule ( rule, suborder_num, suborder_xyz, suborder_w );
    /*
    Expand the suborder information to a full order rule.
    */
    o = 0;

    for ( s = 0; s < suborder_num; ++s )
    {
        if ( suborder[s] == 1 )
        {
            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;
        }
        /*
        Fourfold symmetry on (A,A,A,B)

        123 AAA
        124 AAB
        142 ABA
        412 BAA
        */
        else if ( suborder[s] == 4 )
        {
            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;
        }
        /*
        Sixfold symmetry on (A,A,B,B):

        123 (A,A,B)
        132 (A,B,A),
        134 (A,B,B)
        312 (B,A,A)
        314 (B,A,B)
        341 (B,B,A)
        */
        else if ( suborder[s] == 6 )
        {
            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[0+s*4];
            w[o] = suborder_w[s];
            ++ o;
        }
        /*
        Twelvefold symmetry on (A,A,B,C):

        123 (A,A,B)
        124 (A,A,C)
        132 (A,B,A)
        134 (A,B,C)
        142 (A,C,A)
        143 (A,C,B)
        312 (B,A,A)
        314 (B,A,C)
        341 (B,C,A)
        412 (C,A,A)
        413 (C,A,B)
        431 (C,B,A)
        */
        else if ( suborder[s] == 12 )
        {
            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[0+s*4];
            w[o] = suborder_w[s];
            ++ o;
        }
        /*
        24 fold symmetry on (A,B,C,D):

        123 (A,B,C)
        124 (A,B,D)
        132 (A,C,B)
        134 (A,C,D)
        142 (A,D,B)
        143 (A,D,C)
        213 (B,A,C)
        214 (B,A,D)
        231 (B,C,A)
        234 (B,C,D)
        241 (B,D,A)
        243 (B,D,C)
        312 (C,A,B)
        314 (C,A,D)
        321 (C,B,A)
        324 (C,B,D)
        341 (C,D,A)
        342 (C,D,B)
        412 (D,A,B)
        413 (D,A,C)
        421 (D,B,A)
        423 (D,B,C)
        431 (D,C,A)
        432 (D,C,B)
        */
        else if ( suborder[s] == 24 )
        {
            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[0+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[1+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[1+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[4+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[1+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[0+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[1+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[1+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[0+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[1+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[0+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[3+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[0+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[2+s*4];
            xyz[1+o*3] = suborder_xyz[3+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[0+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[0+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[1+s*4];
            xyz[2+o*3] = suborder_xyz[2+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[0+s*4];
            w[o] = suborder_w[s];
            ++ o;

            xyz[0+o*3] = suborder_xyz[3+s*4];
            xyz[1+o*3] = suborder_xyz[2+s*4];
            xyz[2+o*3] = suborder_xyz[1+s*4];
            w[o] = suborder_w[s];
            ++ o;
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
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_rule_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_RULE_NUM returns the number of NCC rules available.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Output, int TETRAHEDRON_NCC_RULE_NUM, the number of rules available.
*/
{
	static dim_typ result = 7;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _tetrahedron_ncc_suborder ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBORDER returns the suborders for an NCC rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
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
    Output, int TETRAHEDRON_NCC_SUBORDER[SUBORDER_NUM],
    the suborders of the rule.
*/
{
	dim_typ * const a_data = data;
	const register dim_typ rule = a_data[0];
	const register dim_typ suborder_num = a_data[1];
	
    int *suborder = ( int * ) malloc ( suborder_num * sizeof ( int ) );

    if ( rule == 1 )
        suborder[0] = 1;
    else if ( rule == 2 )
        suborder[0] = 4;
    else if ( rule == 3 )
    {
        suborder[0] = 4;
        suborder[1] = 6;
    }
    else if ( rule == 4 )
    {
        suborder[0] =  4;
        suborder[1] = 12;
        suborder[2] =  4;
    }
    else if ( rule == 5 )
    {
        suborder[0] =  4;
        suborder[1] = 12;
        suborder[2] =  6;
        suborder[3] = 12;
        suborder[4] =  1;
    }
    else if ( rule == 6 )
    {
        suborder[0] =  4;
        suborder[1] = 12;
        suborder[2] = 12;
        suborder[3] = 12;
        suborder[4] = 12;
        suborder[5] =  4;
    }
    else if ( rule == 7 )
    {
        suborder[0] =  4;
        suborder[1] = 12;
        suborder[2] = 12;
        suborder[3] = 12;
        suborder[4] =  6;
        suborder[5] = 24;
        suborder[6] =  4;
        suborder[7] =  4;
        suborder[8] =  6;
    }
    else
        return NULL;

    return suborder;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_suborder_num ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBORDER_NUM returns the number of suborders for an NCC rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int RULE, the index of the rule.
    Output, int TETRAHEDRON_NCC_SUBORDER_NUM, the number of suborders
    of the rule.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register dim_typ rule = *(dim_typ *) data;
	
	dim_typ value = USHRT_MAX;
    switch(rule)
    {
        case 1:
        case 2:
            value = 1;
            break;
        case 3:
            value = 2;
            break;
        case 4:
            value = 3;
            break;
        case 5:
            value = 5;
            break;
        case 6:
            value = 6;
            break;
        case 7:
            value = 9;
            break;
    }
    
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_subrule ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBRULE returns a compressed NCC rule.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
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
    Output, double SUBORDER_XYZ[4*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.
    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
	const _2dt2pit * const s_data = data; 
	const register dim_typ rule = s_data->a0;
	const register dim_typ suborder_num = s_data->a1;
	ityp * suborder_xyz = s_data->a2;
	ityp * suborder_w = s_data->a3;
	
    dim_typ i, s;
    int suborder_w_d;
    int *suborder_w_n;
    int suborder_xyz_d;
    int *suborder_xyz_n;

    suborder_xyz_n = ( int * ) malloc ( suborder_num * sizeof ( int ) << 2 );
    suborder_w_n = ( int * ) malloc ( suborder_num * sizeof ( int ) );

    switch(rule)
    {

        case 1:
            tetrahedron_ncc_subrule_01 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 2:
            tetrahedron_ncc_subrule_02 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 3:
            tetrahedron_ncc_subrule_03 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 4:
            tetrahedron_ncc_subrule_04 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 5:
            tetrahedron_ncc_subrule_05 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 6:
            tetrahedron_ncc_subrule_06 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        case 7:
            tetrahedron_ncc_subrule_07 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,suborder_w_n, &suborder_w_d );
            break;
        default:
            return NULL;
    }

    for ( s = 0; s < suborder_num; ++s )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
            suborder_xyz[i+(s<<2)] =( ityp ) ( suborder_xyz_n[i+(s<<2)] )/ ( ityp ) ( suborder_xyz_d );

    for ( s = 0; s < suborder_num; ++s)
        suborder_w[s] = ( ityp  ) suborder_w_n[s] / ( ityp ) suborder_w_d;

    free ( suborder_w_n );
    free ( suborder_xyz_n );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_subrule_01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBRULE_01 returns a compressed NCC rule 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
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
	
    dim_typ i, s;
    int suborder_xyz_n_01[4] =
    {
        1, 1, 1, 1
    };
    const register dim_typ suborder_xyz_d_01 = 4;
    int suborder_w_n_01[1] = { 1 };
    dim_typ suborder_w_d_01 = 1;

    for ( s = 0; s < suborder_num; ++s )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
            suborder_xyz_n[i+(s<<2)] = suborder_xyz_n_01[i+(s<<2)];
    *suborder_xyz_d = suborder_xyz_d_01;

    for ( s = 0; s < suborder_num; ++s )
        suborder_w_n[s] = suborder_w_n_01[s];
    *suborder_w_d = suborder_w_d_01;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_subrule_02 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBRULE_02 returns a compressed NCC rule 2.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
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
	
    dim_typ i, s;
    int suborder_xyz_n_02[4] =
    {
        0, 0, 0, 1
    };
    const register dim_typ suborder_xyz_d_02 = 1;
    int suborder_w_n_02[1] =
    {
        1
    };
    const register dim_typ suborder_w_d_02 = 4;

    for ( s = 0; s < suborder_num; ++s )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
            suborder_xyz_n[i+(s<<2)] = suborder_xyz_n_02[i+(s<<2)];
    *suborder_xyz_d = suborder_xyz_d_02;

    for ( s = 0; s < suborder_num; ++s )
        suborder_w_n[s] = suborder_w_n_02[s];
    *suborder_w_d = suborder_w_d_02;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_subrule_03 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBRULE_03 returns a compressed NCC rule 3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
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
	
    dim_typ i, s;
    int suborder_xyz_n_03[8] =
    {
        0, 0, 0, 2,
        1, 1, 0, 0
    };
    const register dim_typ suborder_xyz_d_03 = 2;
    int suborder_w_n_03[2] =
    {
        -1, 4
    };
    const register dim_typ suborder_w_d_03 = 20;

    for ( s = 0; s < suborder_num; ++s )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i)
            suborder_xyz_n[i+(s<<2)] = suborder_xyz_n_03[i+(s<<2)];
    *suborder_xyz_d = suborder_xyz_d_03;

    for ( s = 0; s < suborder_num; ++s )
        suborder_w_n[s] = suborder_w_n_03[s];
    *suborder_w_d = suborder_w_d_03;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_subrule_04 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBRULE_04 returns a compressed NCC rule 4.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
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
	
    dim_typ i, s;
    int suborder_xyz_n_04[12] =
    {
        0, 0, 0, 3,
        0, 0, 1, 2,
        1, 1, 1, 0
    };
    const register dim_typ suborder_xyz_d_04 = 3;
    int suborder_w_n_04[3] =
    {
        1, 0, 9
    };
    const register dim_typ suborder_w_d_04 = 40;

    for ( s = 0; s < suborder_num; ++s )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
                suborder_xyz_n[i+(s<<2)] = suborder_xyz_n_04[i+(s<<2)];
    *suborder_xyz_d = suborder_xyz_d_04;

    for ( s = 0; s < suborder_num; ++s )
    suborder_w_n[s] = suborder_w_n_04[s];
    *suborder_w_d = suborder_w_d_04;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_subrule_05 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBRULE_05 returns a compressed NCC rule 5.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    , int SUBORDER_XYZ_N[4*SUBORDER_NUM],
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
	
    dim_typ i, s;
    int suborder_xyz_n_05[20] =
    {
        0, 0, 0, 4,
        0, 0, 3, 1,
        2, 2, 0, 0,
        1, 1, 0, 2,
        1, 1, 1, 1
    };
    const register dim_typ suborder_xyz_d_05 = 4;
    int suborder_w_n_05[5] =
    {
         -5, 16, -12, 16, 128
    };
    const register dim_typ suborder_w_d_05 = 420;

    for ( s = 0; s < suborder_num; ++s )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
            suborder_xyz_n[i+(s<<2)] = suborder_xyz_n_05[i+(s<<2)];
    *suborder_xyz_d = suborder_xyz_d_05;

    for ( s = 0; s < suborder_num; ++s )
        suborder_w_n[s] = suborder_w_n_05[s];
    *suborder_w_d = suborder_w_d_05;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_ncc_subrule_06 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBRULE_06 returns a compressed NCC rule 6.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
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
	
    dim_typ i, s;
    int suborder_xyz_n_06[24] =
    {
        0, 0, 0, 5,
        0, 0, 4, 1,
        0, 0, 3, 2,
        1, 1, 0, 3,
        2, 2, 1, 0,
        1, 1, 1, 2
    };
    const register dim_typ suborder_xyz_d_06 = 5;
    int suborder_w_n_06[6] =
    {
        33, -35, 35, 275, -75, 375
    };
    const register dim_typ suborder_w_d_06 = 4032;

    for ( s = 0; s < suborder_num; ++s )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
            suborder_xyz_n[i+(s<<2)] = suborder_xyz_n_06[i+(s<<2)];
    *suborder_xyz_d = suborder_xyz_d_06;

    for ( s = 0; s < suborder_num; ++s)
        suborder_w_n[s] = suborder_w_n_06[s];
    *suborder_w_d = suborder_w_d_06;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _tetrahedron_ncc_subrule_07 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_NCC_SUBRULE_07 returns a compressed NCC rule 7.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    31 January 2007
  Author:
    John Burkardt
  Reference:
    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.
  Parameters:
    Input, int SUBORDER_NUM, the number of suborders of the rule.
    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
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
	
    dim_typ i, s;
    int suborder_xyz_n_07[36] =
    {
        0, 0, 0, 6,
        0, 0, 5, 1,
        0, 0, 4, 2,
        1, 1, 0, 4,
        3, 3, 0, 0,
        3, 2, 1, 0,
        1, 1, 1, 3,
        2, 2, 2, 0,
        2, 2, 1, 1
    };
    const register dim_typ suborder_xyz_d_07 = 6;
    int suborder_w_n_07[9] =
    {
        -7, 24, -30, 0, 40, 30, 180, -45, 0
    };
    const register dim_typ suborder_w_d_07 = 1400;

    for ( s = 0; s < suborder_num; ++s )
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
            suborder_xyz_n[i+(s<<2)] = suborder_xyz_n_07[i+(s<<2)];
    *suborder_xyz_d = suborder_xyz_d_07;

    for ( s = 0; s < suborder_num; ++s  )
    suborder_w_n[s] = suborder_w_n_07[s];
    *suborder_w_d = suborder_w_d_07;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _tetrahedron_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    TETRAHEDRON_VOLUME  computes the volume of a tetrahedron.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 August 2005
  Author:
    John Burkardt
  Parameters:
    Input, double TETRA[3*4], the vertices of the tetrahedron.
    Output, double TETRAHEDRON_VOLUME, the volume of the tetrahedron.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * tetra = data;
	
    ityp a[16];
    dim_typ i, j;
    ityp volume;

    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i )
        #pragma omp parallel for num_threads(4)
        for ( j = 0; j < 4; ++j )
            a[i+j*4] = tetra[i+j*3];

    i = 3;
    #pragma omp parallel for num_threads(4)
    for ( j = 0; j < 4; ++j )
    a[i+(j<<2)] = 1.00;

	result = fabs ( r8mat_det_4d ( a ) ) / 6.00;
    return &result;
}

#endif
