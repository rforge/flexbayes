#include "R.h"
#include "Rmath.h"

#include "Const.h"
#include "MetropolisHastings.h"


extern int randomStart( int n );


MetropolisHastings::MetropolisHastings()
{
  burn_in_time = 0;
  simulations_to_keep = 0;
  sample_frequency = 0;
  simulations_performed = 0;
  simulations_kept = 0;
  vars_stats = NULL;
  accept_count = NULL;
  save_stats = false;
}//end


MetropolisHastings::~MetropolisHastings()
{
  if ( vars_stats != NULL )
  {
    delete vars_stats;
  }

  if ( accept_count != NULL )
  {
    delete accept_count;
  }

}//end

  

/* constructor.
   initialize loop control variables.
   burn_in_time:  is the number of initial samples to discard.
   simulations_to_keep: is the number of samples to save as output.
   sample_frequency: is the subsample frequency among simulations.
*/
MetropolisHastings::MetropolisHastings( int burn_in, int n_simul, int freq )
{
  burn_in_time = burn_in;
  simulations_to_keep = n_simul;
  sample_frequency = freq;

  simulations_performed = 0;
  simulations_kept = 0;
  vars_stats = NULL;
  accept_count = NULL;
  save_stats = false;

  //printf("MetropolisHastings: burn in = %d,  n simul = %d, freq = %d\n", burn_in_time,simulations_to_keep,sample_frequency);

}//end


void MetropolisHastings::saveStatistics( bool val )
{
  save_stats = val;
}//end


/* perform MetropolisHastings sampler loop.
   actual sampling is done in metropolisHastingsStep() below.
*/
void MetropolisHastings::doMetropolisHastings( BayesianModel * model )
{
  //printf( "MetropolisHastings::doMetropolisHastings: Started MetropolisHastings sampler for this model [%s].\n", model->name() ); fflush(stdout);

#ifdef DEBUG1
  printf("MH: Number of variables = %d\n", model->numberOfVariables());  fflush(stdout);
#endif
  if ( simulations_to_keep > 0 )
  {
    //create vector of statistics on drawings
    if ( save_stats )
    {
      vars_stats = new CVector( model->numberOfVariables() );
      vars_stats->setToZero();
      accept_count = new CVector( model->numberOfVariables() );
      accept_count->setToZero();
    }

    //get memory for output
    //printf("create output\n"); fflush(stdout);

    model->createOutput( simulations_to_keep );

    // Names the parameters. 
    //Initializes the linear part of the link function.
    model->initializeTemporaryStructures();

    //if data augmentation is used some draws may be necessary
    model->dataAugmentationInitialDraws();

    //printf("starting metropolisHastings\n"); fflush( stdout );

    simulations_kept = 0;

    while ( !IsSamplerDone() )
    {
      metropolisHastingsStep( model );

      
      //if ( simulations_performed % 1 == 0 )
      //{
      //  printf("[ %d ] ", simulations_performed); fflush(stdout);
      //}
      
      if ( IsASampleToKeep() )
      {      
        //printf("MH: will keep sample %d\n", simulations_kept ); fflush(stdout);

        model->keepSimulation( simulations_kept );
        simulations_kept++;
      }


    }//end loop



    //printf( "*END MetropolisHastings*\n\n" );

    //printDrawingStats();
  }//end if

}//end


/* returns: true is this MetropolisHastings sampler simulations are to be kept as output,
            false, otherwise
*/
bool MetropolisHastings::IsASampleToKeep()
{
  bool keep_draw;
  simulations_performed++;

  keep_draw = ( simulations_performed > burn_in_time 
                && simulations_performed % sample_frequency == 0 )? true : false;

  return keep_draw;
}//end


bool MetropolisHastings::discardSample()
{
  bool discard_draw;

  discard_draw = ( (simulations_performed + 1) <= burn_in_time 
                   || (simulations_performed + 1) % sample_frequency != 0 )? true : false;

  return discard_draw;
}//end


/* returns: true, if total number of samples to keep has been reached.
 */
bool MetropolisHastings::IsSamplerDone()
{
  bool keep_mh;

  keep_mh = ( simulations_kept < simulations_to_keep )? false : true;
  return keep_mh;
}//end


/* Does one full iteration of MetropolisHastings sampler.
   Drawing from the full conditionals is done in sequential
   order starting from one variable chosen at random.
   This ensures reversibility of the chain, and hence
   better convergence properties.
*/
void MetropolisHastings::metropolisHastingsStep( BayesianModel * model )
{
  int iVar, start_index, index;
  double U, log_ratio_proposal, log_ratio_target, alpha;

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
    //printf("MH: will update var from proposal [%d]\n", index); fflush(stdout);

    model->updateVariableForProposal( index );

    //draw from proposal
    //printf("MH: will draw var from proposal [%d]\n", index); fflush(stdout);

    model->drawVariableFromProposal( index );

    //compute acceptance ratio
    //printf("MH: will log Ratio target density [%d]\n", index); fflush(stdout);

    log_ratio_target = model->logRatioTargetDensity( index );

    //printf("MH: will log ratio proposal [%d]\n", index); fflush(stdout);

    log_ratio_proposal = model->logRatioProposal( index );
    alpha = log_ratio_target - log_ratio_proposal;

    //draw a uniform number
    U = runif( 0.0, 1.0 );


    //decide
    if ( log( U ) <= alpha )
    {

      //printf(" MH: will update current  variable [%d]\n", index ); fflush(stdout);

      model->setCurrentVariableFromProposal( index );

      //printf(" MH: updated current  variable [%d]\n", index ); fflush(stdout);

      if ( !discardSample() && save_stats )
      {
        accept_count->Val( index ) = accept_count->Val( index ) + 1;
      }

    }
    else  //reject proposal
    {

      //printf(" MH: will keep variable [%d]\n", index ); fflush(stdout);

      model->keepCurrentVariable( index );

      //printf(" MH: kept variable [%d]\n", index ); fflush(stdout);

    }

    //printf( "MH: end of inner loop[ %d ]\n", index ); fflush(stdout);
  }//end for

  //printf( "MetropolisHastings: finished\n" ); fflush(stdout);

}//end


void MetropolisHastings::printDrawingStats()
{
  int i;
  double total;

  /*
  if ( vars_stats != NULL )
  {
    printf( "\nMetropolisHastings: Drawing Statistics: \n" );
    total = vars_stats->Len() * vars_stats->Mean();
    for ( i = 0; i < vars_stats->Len(); i++ )
    {
      printf( "variable[ %d ] = %f (%f out of %f).\n", i, vars_stats->Val(i) / total, vars_stats->Val(i), total );
    }

    printf( "\n" );

  }
  */

  if ( accept_count != NULL )
  {
    printf( "\nMetropolisHastings: Acceptance Rate Statistics: \n" );
    for ( i = 0; i < accept_count->Len(); i++ )
    {
      total = simulations_to_keep;
      printf( "variable[ %d ] = %f (%f out of %f).\n", i, accept_count->Val(i) / total, accept_count->Val(i), total );
    }

    printf( "\n" );

  }

}//end


void MetropolisHastings::printDrawingStats( double * stats )
{
  int i;

  /*
  if ( vars_stats != NULL )
  {
    for ( i = 0; i < vars_stats->Len(); i++ )
    {
      stats[i] = vars_stats->Val(i);
    }
  }
  */

  if ( accept_count != NULL )
  {
    for ( i = 0; i < accept_count->Len(); i++ )
    {
      stats[i] = accept_count->Val(i);
    }
  }


}//end
