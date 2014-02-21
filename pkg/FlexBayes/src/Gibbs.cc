#include "R.h"
#include "Rmath.h"

#include "Const.h"
#include "Gibbs.h"


int randomStart( int n );


Gibbs::Gibbs()
{
  burn_in_time = 0;
  simulations_to_keep = 0;
  sample_frequency = 0;
  simulations_performed = 0;
  simulations_kept = 0;
  vars_stats = NULL;
  save_stats = false;
}//end


Gibbs::~Gibbs()
{
  if ( vars_stats != NULL )
  {
    delete vars_stats;
  }

}//end

  

/* constructor.
   initialize loop control variables.
   burn_in_time:  is the number of initial samples to discard.
   simulations_to_keep: is the number of samples to save as output.
   sample_frequency: is the subsample frequency among simulations.
*/
Gibbs::Gibbs( int burn_in, int n_simul, int freq )
{
  burn_in_time = burn_in;
  simulations_to_keep = n_simul;
  sample_frequency = freq;

  simulations_performed = 0;
  simulations_kept = 0;
  vars_stats = NULL;
  save_stats = false;

  //printf("Gibbs: burn in = %d,  n simul = %d, freq = %d\n", burn_in_time,simulations_to_keep,sample_frequency);

}//end


void Gibbs::saveStatistics( bool val )
{
  save_stats = val;
}//end


/* perform Gibbs sampler loop.
   actual sampling is done in gibbsStep() below.
*/
void Gibbs::doGibbs( BayesianModel * model )
{
  //printf( "Gibbs::doGibbs: Started Gibbs sampler for this model [%s].\n", model->name() );
  //fflush(stdout);

  if ( simulations_to_keep > 0 )
  {
    //create vector of statistics on drawings

    //printf("in Gibbs: number of vars is %d\n", model->numberOfVariables() ); fflush(stdout);

    if ( save_stats )
    {
      vars_stats = new CVector( model->numberOfVariables() );
      vars_stats->setToZero();
    }

    //get memory for output
    //printf("create output\n"); fflush(stdout);

    model->createOutput( simulations_to_keep );

    //temporary storage
    //printf("init temp\n"); fflush( stdout );

    model->initializeTemporaryStructures();

    //if data augmentation is used some draws may be necessary
    model->dataAugmentationInitialDraws();

    //printf("starting gibbs\n"); fflush( stdout );

    simulations_kept = 0;
    while ( !IsSamplerDone() )
    {
      gibbsStep( model );

      //if ( simulations_performed % 1 == 0 )
      //{
      //  printf("[ %d ] ", simulations_performed); fflush(stdout);
      //}

      if ( IsASampleToKeep() )
      {      
        model->keepSimulation( simulations_kept );
        simulations_kept++;
        //printf("[%d, %d]  ", simulations_kept, simulations_performed );
      }
    }//end loop

    //printf( "*END Gibbs*\n\n" );

    //printDrawingStats();
  }//end if

}//end


/* returns: true is this gibbs sampler simulations are to be kept as output,
            false, otherwise
*/
bool Gibbs::IsASampleToKeep()
{
  bool keep_draw;
  simulations_performed++;

  keep_draw = ( simulations_performed > burn_in_time 
                && simulations_performed % sample_frequency == 0 )? true : false;

  return keep_draw;
}//end


/* returns: true, if total number of samples to keep has been reached.
 */
bool Gibbs::IsSamplerDone()
{
  bool keep_gibbs;

  keep_gibbs = ( simulations_kept < simulations_to_keep )? false : true;
  return keep_gibbs;
}//end


/* Does one full iteration of Gibbs sampler.
   Drawing from the full conditionals is done in sequential
   order starting from one variable chosen at random.
   This ensures reversibility of the chain, and hence
   better convergence properties.
*/
void Gibbs::gibbsStep( BayesianModel * model )
{
  int iVar, start_index, index;

  //get starting variable
  start_index = randomStart( model->numberOfVariables() );
  //keep statistics of drawings
  if ( save_stats )
  {  
    vars_stats->Val( start_index ) += 1;
  }

  for ( iVar = 0; iVar < model->numberOfVariables(); iVar++ )
  {
    index = ( start_index + iVar ) %  model->numberOfVariables();

    //update index-th distribution 
    //printf("Gibbs: update\n"); fflush(stdout);

    model->fullConditionalUpdateVariable( index );

    //printf("Gibbs: will draw\n"); fflush(stdout); fflush(stdout);

    //draw from it
    model->drawVariable( index );
  }

}//end


void Gibbs::printDrawingStats()
{
  int i;
  double total;

  if ( vars_stats != NULL )
  {
    Rprintf( "\nGibbs: Drawing Statistics: \n" );
    total = vars_stats->Len() * vars_stats->Mean();
    for ( i = 0; i < vars_stats->Len(); i++ )
    {
      Rprintf( "variable[ %d ] = %f (%f out of %f).\n", i, vars_stats->Val(i) / total, vars_stats->Val(i), total );
    }

    Rprintf( "\n" );

  }

}//end


void Gibbs::printDrawingStats( double * stats )
{
  int i;

  if ( vars_stats != NULL )
  {
    for ( i = 0; i < vars_stats->Len(); i++ )
    {
      stats[i] = vars_stats->Val(i);
    }
  }
}//end


int randomStart( int n )
{
  double random_index;
  int index, start_index;

  random_index = runif( 0.0, 1.0 );
  random_index *= (double) n;
  index = (int) random_index;
  if ( ( (double) index ) - random_index == 0 )
  {
    start_index = index;
    if ( start_index == 0 )
    {
      start_index = 1;
    }
  }
  else //index < random_index
  {
    start_index = index + 1;
  }
  start_index--;//it starts at 0

  return start_index;
}//end
