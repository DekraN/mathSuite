#ifndef __DISABLEDEEP_TRIANGLEFELIPPARULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_unit_o01 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_UNIT_O01 returns a 1 point quadrature rule for the unit triangle.
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
        0.33333333333333333333,
        0.33333333333333333333
    };

    r8vec_copy ( 1, w_save, w );
    r8vec_copy ( 2, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_unit_o03 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_UNIT_O03 returns a 3 point quadrature rule for the unit triangle.
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
__MATHSUITE __JBURKARDT  void *  _triangle_unit_o03b ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_UNIT_O03B returns a 3 point quadrature rule for the unit triangle.
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
__MATHSUITE __JBURKARDT  void *   _triangle_unit_o06 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_UNIT_O06 returns a 6 point quadrature rule for the unit triangle.
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
__MATHSUITE __JBURKARDT  void *   _triangle_unit_o06b ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_UNIT_O06B returns a 6 point quadrature rule for the unit triangle.
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
    ityp xy_save[2*6] =
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
__MATHSUITE __JBURKARDT  void * _triangle_unit_o07 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_UNIT_O07 returns a 7 point quadrature rule for the unit triangle.
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
        0.22500000000000000000
    };
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
    r8vec_copy ( 14, xy_save, xy );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _triangle_unit_o12 ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_UNIT_O12 returns a 12 point quadrature rule for the unit triangle.
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
__MATHSUITE __JBURKARDT  void *   _triangle_unit_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    TRIANGLE_UNIT_VOLUME returns the "volume" of the unit triangle in 2D.
  Discussion:
    The "volume" of a triangle is usually called its area.
    The integration region is:
      0 <= X,
      0 <= Y,
      X + Y <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 March 2008
  Author:
    John Burkardt
  Parameters:
    Output, double TRIANGLE_UNIT_VOLUME, the volume.
*/
{
	static ityp result = 0.50;
    return &result;
}

#endif
