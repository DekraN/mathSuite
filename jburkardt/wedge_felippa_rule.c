#ifndef __DISABLEDEEP_WEDGEFELIPPARULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_o01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_O01 returns a 1 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[1], the weights.
    Output, double X[1], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[1] =
    {
        1.00
    };
    ityp x_save[1] =
    {
         0.00
    };

    r8vec_copy ( 1, w_save, w );
    r8vec_copy ( 1, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_o02 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_O02 returns a 2 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[2], the weights.
    Output, double X[2], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[2] =
    {
        0.50,
        0.50
    };
    ityp x_save[2] =
    {
        -0.57735026918962576451,
        0.57735026918962576451
    };

    r8vec_copy ( 2, w_save, w );
    r8vec_copy ( 2, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_o03 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_O03 returns a 3 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[3], the weights.
    Output, double X[3], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[3] =
    {
        0.27777777777777777777,
        0.44444444444444444444,
        0.27777777777777777777
    };
    ityp x_save[3] =
    {
        -0.77459666924148337704,
        0.00000000000000000000,
        0.77459666924148337704
    };

    r8vec_copy ( 3, w_save, w );
    r8vec_copy ( 3, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_o04 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_O04 returns a 4 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[4], the weights.
    Output, double X[4], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[4] =
    {
        0.173927422568727,
        0.326072577431273,
        0.326072577431273,
        0.173927422568727
    };
    ityp x_save[4] =
    {
        -0.86113631159405257522,
        -0.33998104358485626480,
        0.33998104358485626480,
        0.86113631159405257522
    };

    r8vec_copy ( 4, w_save, w);
    r8vec_copy ( 4, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _line_o05 ( void * data)
/******************************************************************************/
/*
  Purpose:
    LINE_O05 returns a 5 point quadrature rule for the unit line.
  Discussion:
    The integration region is:
    - 1.0 <= X <= 1.0
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Output, double W[5], the weights.
    Output, double X[5], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * x = a_data[1];
	
    ityp w_save[5] =
    {
        0.118463442528095,
        0.239314335249683,
        0.284444444444444,
        0.239314335249683,
        0.118463442528095
    };
    ityp x_save[5] =
    {
        -0.90617984593866399280,
        -0.53846931010568309104,
        0.00000000000000000000,
        0.53846931010568309104,
        0.90617984593866399280
    };

    r8vec_copy ( 5, w_save, w );
    r8vec_copy ( 5, x_save, x );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_o01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_O01 returns a 1 point quadrature rule for the unit triangle.
  Discussion:
    This rule is precise for monomials through degree 1.
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 April 2009
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Output, double W[1], the weights.
    Output, double XY[2*1], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * xy = a_data[1];

    ityp w_save[1] =
    {
        1.00
    };
    ityp xy_save[2] =
    {
        0.33333333333333333333,  0.33333333333333333333
    };

    r8vec_copy ( 1, w_save, w );
    r8vec_copy ( 2, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_o03 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_O03 returns a 3 point quadrature rule for the unit triangle.
  Discussion:
    This rule is precise for monomials through degree 2.
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 April 2009
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Output, double W[3], the weights.
    Output, double XY[2*3], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * xy = a_data[1];

    ityp w_save[3] =
    {
        0.33333333333333333333,
        0.33333333333333333333,
        0.33333333333333333333
    };
    ityp xy_save[6] =
    {
        0.66666666666666666667,  0.16666666666666666667,
        0.16666666666666666667,  0.66666666666666666667,
        0.16666666666666666667,  0.16666666666666666667
    };

    r8vec_copy ( 3, w_save, w );
    r8vec_copy ( 6, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_o03b ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_O03B returns a 3 point quadrature rule for the unit triangle.
  Discussion:
    This rule is precise for monomials through degree 2.
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 April 2009
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Output, double W[3], the weights.
    Output, double XY[2*3], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * xy = a_data[1];

    ityp w_save[3] =
    {
        0.33333333333333333333,
        0.33333333333333333333,
        0.33333333333333333333
    };
    ityp xy_save[6] =
    {
        0.00,  0.50,
        0.50,  0.00,
        0.50,  0.50
    };

    r8vec_copy ( 3, w_save, w );
    r8vec_copy ( 6, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_o06 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_O06 returns a 6 point quadrature rule for the unit triangle.
  Discussion:
    This rule is precise for monomials through degree 4.
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 April 2009
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Output, double W[6], the weights.
    Output, double XY[2*6], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * xy = a_data[1];
	
    ityp w_save[6] =
    {
        0.22338158967801146570,
        0.22338158967801146570,
        0.22338158967801146570,
        0.10995174365532186764,
        0.10995174365532186764,
        0.10995174365532186764
    };
    ityp xy_save[12] =
    {
        0.10810301816807022736,  0.44594849091596488632,
        0.44594849091596488632,  0.10810301816807022736,
        0.44594849091596488632,  0.44594849091596488632,
        0.81684757298045851308,  0.091576213509770743460,
        0.091576213509770743460,  0.81684757298045851308,
        0.091576213509770743460,  0.091576213509770743460
    };

    r8vec_copy ( 6, w_save, w );
    r8vec_copy ( 12, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_o06b ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_O06B returns a 6 point quadrature rule for the unit triangle.
  Discussion:
    This rule is precise for monomials through degree 3.
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 April 2009
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Output, double W[6], the weights.
    Output, double XY[2*6], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * xy = a_data[1];
	
    ityp w_save[6] =
    {
        0.30000000000000000000,
        0.30000000000000000000,
        0.30000000000000000000,
        0.033333333333333333333,
        0.033333333333333333333,
        0.033333333333333333333
    };
    ityp xy_save[12] =
    {
        0.66666666666666666667,  0.16666666666666666667,
        0.16666666666666666667,  0.66666666666666666667,
        0.16666666666666666667,  0.16666666666666666667,
        0.00,  0.50,
        0.50,  0.00,
        0.50,  0.50
    };

    r8vec_copy ( 6, w_save, w );
    r8vec_copy ( 12, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_o07 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_O07 returns a 7 point quadrature rule for the unit triangle.
  Discussion:
    This rule is precise for monomials through degree 5.
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    16 April 2009
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Output, double W[7], the weights.

    Output, double XY[2*7], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * xy = a_data[1];

    ityp w_save[7] =
    {
        0.12593918054482715260,
        0.12593918054482715260,
        0.12593918054482715260,
        0.13239415278850618074,
        0.13239415278850618074,
        0.13239415278850618074,
        0.22500000000000000000 }
    ;
    ityp xy_save[14] =
    {
        0.79742698535308732240,  0.10128650732345633880,
        0.10128650732345633880,  0.79742698535308732240,
        0.10128650732345633880,  0.10128650732345633880,
        0.059715871789769820459,  0.47014206410511508977,
        0.47014206410511508977,  0.059715871789769820459,
        0.47014206410511508977,  0.47014206410511508977,
        0.33333333333333333333,  0.33333333333333333333
    };

    r8vec_copy ( 7, w_save, w );
    r8vec_copy (14, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_o12 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_O12 returns a 12 point quadrature rule for the unit triangle.
  Discussion:
    This rule is precise for monomials through degree 6.
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    19 April 2009
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Output, double W[12], the weights.
    Output, double XY[2*12], the abscissas.
*/
{
	ityp ** const a_data = data;
	ityp * w = a_data[0];
	ityp * xy = a_data[1];

    ityp w_save[12] =
    {
        0.050844906370206816921,
        0.050844906370206816921,
        0.050844906370206816921,
        0.11678627572637936603,
        0.11678627572637936603,
        0.11678627572637936603,
        0.082851075618373575194,
        0.082851075618373575194,
        0.082851075618373575194,
        0.082851075618373575194,
        0.082851075618373575194,
        0.082851075618373575194
    };
    ityp xy_save[24] =
    {
        0.87382197101699554332,  0.063089014491502228340,
        0.063089014491502228340,  0.87382197101699554332,
        0.063089014491502228340,  0.063089014491502228340,
        0.50142650965817915742,  0.24928674517091042129,
        0.24928674517091042129,  0.50142650965817915742,
        0.24928674517091042129,  0.24928674517091042129,
        0.053145049844816947353,  0.31035245103378440542,
        0.31035245103378440542,  0.053145049844816947353,
        0.053145049844816947353,  0.63650249912139864723,
        0.31035245103378440542,  0.63650249912139864723,
        0.63650249912139864723,  0.053145049844816947353,
        0.63650249912139864723,  0.31035245103378440542
    };

    r8vec_copy ( 12, w_save, w );
    r8vec_copy ( 24, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wedge_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEDGE_RULE returns a quadrature rule for the unit wedge.
  Discussion:
    It is usually sensible to take LINE_ORDER and TRIG_ORDER so that
    the line and triangle rules are roughly the same precision.  For that
    criterion, we recommend the following combinations:
      TRIANGLE_ORDER  LINE_ORDER  Precision
      --------------  ----------  ---------
          1               1       1 x 1
          3               2       2 x 3
         -3               2       2 x 3
          6               3       4 x 5
         -6               2       3 x 3
          7               3       5 x 5
         12               4       6 x 7
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1
      -1 <= Z <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 April 2009
  Author:
    John Burkardt
  Reference:
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.
  Parameters:
    Input, int LINE_ORDER, the index of the line rule.
    The index of the rule is equal to the order of the rule.
    1 <= LINE_ORDER <= 5.
    Input, int TRIANGLE_ORDER, the indes of the triangle rule.
    The index of the rule is 1, 3, -3, 6, -6, 7 or 12.
    Output, double W[LINE_ORDER*abs(TRIANGLE_ORDER)], the weights.
    Output, double XYZ[3*LINE_ORDER*abs(TRIANGLE_ORDER)], the abscissas.
*/
{
	const dts2pit * const s_data = data;
	const register dim_typ line_order = s_data->a0;
	const register short triangle_order = s_data->a1;
	ityp * w = s_data->a2;
	ityp * xyz = s_data->a3;
	
    dim_typ i, j, k;
    ityp *line_w;
    ityp *line_x;
    ityp *triangle_w;
    ityp *triangle_xy;

    line_w = ( ityp * ) malloc ( line_order * sizeof ( ityp ) );
    line_x = ( ityp * ) malloc ( line_order * sizeof ( ityp ) );

    switch(line_order)
    {
        case 1:
            line_o01 ( line_w, line_x );
            break;
        case 2:
            line_o02 ( line_w, line_x );
            break;
        case 3:
            line_o03 ( line_w, line_x );
            break;
        case 4:
            line_o04 ( line_w, line_x );
            break;
        case 5:
            line_o05 ( line_w, line_x );
            break;
        default:
            return NULL;
    }

    triangle_w = ( ityp * ) malloc ( abs(triangle_order) * sizeof ( ityp ) );
    triangle_xy = ( ityp * ) malloc (  abs(triangle_order) * sizeof ( ityp ) << 1 );

    switch(triangle_order)
    {
        case 1:
            triangle_o01 ( triangle_w, triangle_xy );
            break;
        case 3:
            triangle_o03 ( triangle_w, triangle_xy );
            break;
        case -3:
            triangle_o03b ( triangle_w, triangle_xy );
            break;
        case 6:
            triangle_o06 ( triangle_w, triangle_xy );
            break;
        case -6:
            triangle_o06b ( triangle_w, triangle_xy );
            break;
        case 7:
            triangle_o07 ( triangle_w, triangle_xy );
            break;
        case 12:
            triangle_o12 ( triangle_w, triangle_xy );
            break;
        default:
            return NULL;
    }

    k = 0;
    for ( i = 0; i < line_order; ++i )
        for ( j = 0; j < abs ( triangle_order ); ++j )
        {
            w[k] = line_w[i] * triangle_w[j];
            xyz[0+k*3] = triangle_xy[0+(j<<1)];
            xyz[1+k*3] = triangle_xy[1+(j<<1)];
            xyz[2+k*3] = line_x[i];
                ++ k;
        }

    free ( line_w );
    free ( line_x );
    free ( triangle_w );
    free ( triangle_xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _wedge01_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
   WEDGE01_SAMPLE samples points uniformly from the unit wedge in 3D.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2014
  Author:
    John Burkardt
  Reference:
    Reuven Rubinstein,
    Monte Carlo Optimization, Simulation, and Sensitivity
    of Queueing Networks,
    Krieger, 1992,
    ISBN: 0894647644,
    LC: QA298.R79.
  Parameters:
    Input, int N, the number of points.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double WEDGE01_SAMPLE[3*N], the points.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    ityp *e;
    ityp e_sum;
    dim_typ i, j;
    ityp *x = ( ityp * ) malloc ( 3 * n * sizeof ( ityp ) );

    for ( j = 0; j < n; j++ )
    {
        e = r8vec_uniform_01_new ( 4, seed );
        e_sum = 0.00;
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
        {
            e[i] = - log ( e[i] );
            e_sum += e[i];
        }

        x[0+j*3] = e[0] / e_sum;
        x[1+j*3] = e[1] / e_sum;
        x[2+j*3] = 2.00 * e[3] - 1.00;

        free ( e );
    }

    return x;
}

#endif
