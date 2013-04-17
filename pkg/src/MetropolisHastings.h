
#ifndef BAYES_METROPOLIS_HASTINGS_H
#define BAYES_METROPOLIS_HASTINGS_H

#include "BayesianModel.h"
#include "Vector.h"

class MetropolisHastings
{
  private:
    int burn_in_time;
    int simulations_to_keep;
    int sample_frequency;
    int simulations_performed;
    int simulations_kept;
    CVector * vars_stats;
    CVector * accept_count;
    bool save_stats;
  public:
    MetropolisHastings();
    MetropolisHastings( int burn_in, int n_simul, int freq );
    virtual ~MetropolisHastings();
    void saveStatistics( bool val );
    void printDrawingStats();
    void printDrawingStats( double * stats );
    inline int burnIn() { return burn_in_time; }
    inline int simulationsToKeep() { return simulations_to_keep; }
    inline int simulationsPerformed() { return simulations_performed; }
    inline int simulationsKept() { return simulations_kept; }
    inline int sampleFrequency() { return sample_frequency; }
    void doMetropolisHastings( BayesianModel * model );
    void metropolisHastingsStep( BayesianModel * model );
    bool IsASampleToKeep();
    bool discardSample();
    bool IsSamplerDone();  
}; //end MetropolisHastings class

#endif /* BAYES_METROPOLIS_HASTINGS_H */
