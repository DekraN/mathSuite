#ifndef __DISABLEDEEP_HIGHCARDSIMULATION

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _high_card_probability ( void * data)
/******************************************************************************/
/*
  Purpose:
    HIGH_CARD_PROBABILITY: winning probabilities for the high card game.
  Discussion:
    The high card game presents the player with a deck of cards, each
    having an unknown value.  The player is allowed to go throught the
    deck once, looking at the cards one at a time.  At any time, the player
    may decide to take a particular card, winning that amount and stopping
    the game.  If the player continues to the end, by default the last card
    indicates the amount won.
    An optimal strategy for selecting the highest card is as follows:
    * look at, but do not select, the first k-1 cards;
    * stop at the first card, from k to n, that is higher than the first
      k-1 cards.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 February 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of cards.
    Output, double P[N].  P[K] is the probability that a strategy
    that skips K cards will win, given that the deck has N cards.
*/
{
	const register dim_typ n = *(dim_typ *) data;
	
    dim_typ i, j;
    ityp *p;
    ityp t;

    p = ( ityp * ) malloc ( n * sizeof ( ityp ) );

    for ( i = 0; i < n; ++i )
    {
        t = 0.00;
        for ( j = i + 1; j < n; ++j )
            t += 1.00 / ( ityp ) ( j );
        p[i] = ( 1.00 + ( ityp ) ( i ) * t ) / ( ityp ) ( n );
    }

    return p;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _high_card_shuffle ( void * data)
/******************************************************************************/
/*
  Purpose:
    HIGH_CARD_SHUFFLE generates a sequence of numeric "cards" for a game.
  Discussion:
    In this game, you know that the deck contains N cards.  You win by
    choosing the highest card in the deck.  You don't know what this card
    is, and you must choose your card by saying "stop" as, one by one,
    the cards of the deck are exposed.
    A random guesser would get the high card with probability 1/N.
    An intelligent guesser can do much better.
    It is the goal of this program so "shuffle" a deck of cards suitable
    for this game.  The problem is that we know the highest card in an
    ordinary deck.  Let's replace the cards by integers.  Then if we know
    in advance the range of the cards (say, they must lie between 1 and
    1,000), it may be true that we can guess the card that is the maximum.
    However, this program produces a sequence of integer card values for
    which no information can be gained from the values.  It does this
    by regarding the card values as binary integers between 1 and 2^N - 1.
    We can make a perfectly information-free sequence as follows:
      Card 1 sets bit N-1 to 1.
      Card 2 sets bit N-2 to 1, bit  N-1 randomly.
      ...
      Card I sets bit N-I to 1, bits N-1 down to N-I+1 randomly.
      ...
      Card N sets bit N-N to 1, bits N-1 down to 1 randomly.

    The I-th card has equal probability to fall in any of the I intervals
    defined by the I-1 previous cards.  So, knowing the previous cards tells
    you absolutely nothing about where the next card will fall, and each
    card is, at the moment you see it, as good a guess for the maximum as
    the unseen cards.
    For example, the command "high_card_shuffle(7)" might yield
      64    96    80     8     4    82    29
    or
      64    32    16    24    12    58    73
    or
      64    96    48     8   100    26    93
    or
      64    96    16    56    76   122    71
    in which the highest card is #2, #7, #5, or #6 respectively.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 February 2014
  Author:
    John Burkardt
  Parameters:
    Input, int N, the number of cards.  N probably needs to
    be less than 32.
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, int SEQUENCE[N], a set of N integer values
    that can be used as the cards in the high card guessing game.
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    dim_typ c;
    dim_typ i;
    dim_typ j;
    dim_typ k;
    dim_typ *sequence;

    if ( 32 <= n )
        return NULL;

    sequence = ( dim_typ * ) malloc ( n * sizeof ( dim_typ ) );

    for ( i = 0; i < n; ++i )
    {
        c = powi ( 2, n - i - 1 );
        for ( j = 0; j < i; ++j )
        {
            k = i4_uniform_ab ( 0, 1, seed );
            c += k * powi ( 2, n - i + j );
        }
        sequence[i] = c;
    }

    return sequence;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _high_card_simulation ( void * data)
/******************************************************************************/
/*
  Purpose:
    HIGH_CARD_SIMULATION simulates a game of choosing the highest card in a deck.
  Discussion:
    You are given a deck of DECK_SIZE cards.
    Your goal is to select the high card.  For convenience, we can assume
    the cards are a permutation of the integers from 1 to DECK_SIZE, but in
    fact the user mustn't see such values or else it's obvious which is the
    largest card.
    However, your choice is made under the following rules:  You may turn over
    one card at a time.  When a card is turned over, you may declare that to be
    your choice, or else turn over another card.  If you have not chosen a card
    by the end, then your choice is the final card.
    If you have no idea what to do, and simply decide in advance to pick
    a card "at random", that is, for example, you decide to pick the 15th card
    before having seen any cards, then your probability of winning is
    1/DECK_SIZE.
    The question is, can you do better than that?
    Your strategy is as follows: always look at the first SKIP_NUM cards
    without choosing them.  Then choose the very next card you encounter
    that is larger than the cards you skipped.
    Using this program, you can easily see that skipping 5 cards is much better
    than picking one at random, skipping 10 is even better, and so on.
    Of course, you can't skip too many cards, and in fact, the results seem
    to be best for somewhere around 30 to 35 cards skipped.  For problems
    like this, the optimal value is somewhere around 1 / e, where E is the
    base of the natural logarithm system.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 February 2014
  Author:
    John Burkardt
  Parameters:
    Input, int DECK_SIZE, the number of cards in the deck.
    2 <= DECK_SIZE.  Default value is 52;
    Input, int SKIP_NUM, the number of initial cards you plan
    to examine but will NOT select.  If SKIP_NUM is 0, you don't look at any
    cards first.  0 <= SKIP_NUM < DECK_SIZE.
    Input, int TRIAL_NUM, the number of times we will
    simulate this process.
    Input/output, int SEED, a seed for the random
    number generator.
    Output, double HIGH_CARD_SIMULATION, the estimated probability that
    your strategy of skipping SKIP_NUM cards and then selecting the next
    card that is bigger, will result in choosing the highest card.
*/
{
	static ityp result = MAX_VAL;
	
	const _3dtpi * const s_data = data;
	const register dim_typ deck_size = s_data->a0;
	register dim_typ skip_num = s_data->a1;
	const register dim_typ trial_num = s_data->a2;
	int * seed = s_data->a3;
	
    dim_typ card;
    int *cards;
    dim_typ choice;
    dim_typ correct;
    ityp p;
    dim_typ skip_max;
    dim_typ trial;
    dim_typ true_max;
    /*
    Check values.
    */
    if ( deck_size < 2 )
    {
    	result = MAX_VAL;
        return &result;
    }

    if ( skip_num < 0 )
        skip_num = 0;

    if ( deck_size <= skip_num || trial_num < 1 )
    {
    	result = MAX_VAL;
        return &result;
    }

    correct = 0;

    for ( trial = 1; trial <= trial_num; ++trial )
    {
        cards = perm_uniform_new ( deck_size, seed );
        skip_max = 1 <= skip_num ? i4vec_max ( skip_num, cards ) : -i4_huge;
        true_max = i4vec_max ( deck_size, cards );
        /*
        In case you don't encounter a card larger than SKIP_MAX,
        we'll assume you pick the last card in the deck, even though
        you know it's a loser.
        */
        choice = cards[deck_size-1];
        /*
        Turn over the remaining cards in the deck, but stop
        immediately when you find one bigger than SKIP_MAX.
        */
        for ( card = skip_num; card < deck_size; ++card )
        {
        if ( skip_max < cards[card] )
            {
                choice = cards[card];
                break;
            }
        }
        /*
        Record successful choices.
        */
        if ( choice == true_max )
            ++ correct;
        free ( cards );
    }
    /*
    Estimate the probability.
    */
    
    result = ( ityp ) ( correct ) / ( ityp ) ( trial_num );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _perm_uniform_new ( void * data)
/******************************************************************************/
/*
  Purpose:
    PERM_UNIFORM_NEW selects a random permutation of N objects.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 February 2014
  Author:
    John Burkardt
  Reference:
    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.
  Parameters:
    Input, int N, the number of objects to be permuted.
    Input/output, int *SEED, a seed for the random number generator.
    Output, int PERM_UNIFORM_NEW[N], a permutation of
 (BASE, BASE+1, ..., BASE+N-1).
*/
{
	const dtpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	int * seed = s_data->a1;
	
    dim_typ i;
    dim_typ j;
    dim_typ k;
    int * p = ( int * ) malloc ( n * sizeof ( int) );

    for ( i = 0; i < n; ++i )
    p[i] = i;


    for ( i = 0; i < n - 1; ++i)
    {
        j = i4_uniform_ab ( i, n - 1, seed );
        k    = p[i];
        p[i] = p[j];
        p[j] = k;
    }

    return p;
}

#endif
