#ifndef __DISABLEDEEP_DUELSIMULATION

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _duel_result ( void * data)
/******************************************************************************/
/*
  Purpose:
    DUEL_RESULT returns the outcome of a single duel.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 September 2012
  Author:
    John Burkardt
  Reference:
    Martin Shubik,
    “Does the Fittest Necessarily Survive?”,
    in Readings in Game Theory and Political Behavior,
    edited by Martin Shubik,
    Doubleday, 1954,
    LC: H61.S53.
  Parameters:
    Input, double A_ACCURACY, B_ACCURACY, the probabilities that player A
    and B will hit their opponent in a single shot.
    Output, int DUEL_RESULT, the survivor of the duel.
*/
{
	static dim_typ result = USHRT_MAX;
	
	ityp * const a_data = data;
	const register ityp a_accuracy = a_data[0];
	const register ityp b_accuracy = a_data[1];
	
    ityp r;
    dim_typ winner;

    while (true)
    {
        r = random_double ( );
        if ( r <= a_accuracy )
        {
            winner = 1;
            break;
        }

        r = random_double ( );

        if ( r <= b_accuracy )
            {
            winner = 2;
            break;
        }
    }
    
    result = winner;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _random_double ( void * data)
/******************************************************************************/
/*
  Purpose:
    RANDOM_DOUBLE returns a random number between 0 and 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 November 2009
  Author:
    John Burkardt
  Parameters:
    Output, double RANDOM_DOUBLE, the random value.
*/
{
	static ityp result = MAX_VAL;
	result = ( ityp ) rand ( ) / ( ityp ) RAND_MAX;
	
    return &result;
}

#endif
