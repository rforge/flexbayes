
#ifndef BAYES_GIBBS_H
#define BAYES_GIBBS_H

#include "BayesianModel.h"
#include "Vector.h"

class Gibbs
{
  private:
    int burn_in_time;
    int simulations_to_keep;
    int sample_frequency;
    int simulations_performed;
    int simulations_kept;
    CVector * vars_stats;
    bool save_stats;
  public:
    Gibbs();
    Gibbs( int burn_in, int n_simul, int freq );
    virtual ~Gibbs();
    void saveStatistics( bool val );
    void printDrawingStats();
    void printDrawingStats( double * stats );
    inline int burnIn() { return burn_in_time; }
    inline int simulationsToKeep() { return simulations_to_keep; }
    inline int simulationsPerformed() { return simulations_performed; }
    inline int simulationsKept() { return simulations_kept; }
    inline int sampleFrequency() { return sample_frequency; }
    void doGibbs( BayesianModel * model );
    void gibbsStep( BayesianModel * model );
    bool IsASampleToKeep();
    bool IsSamplerDone();  
}; //end Gibbs class

#endif /* BAYES_GIBBS_H */
