
#ifndef BAYES_BAYESIAN_MODEL_H
#define BAYES_BAYESIAN_MODEL_H

#include "DistributionParameter.h"

/* abstract declaration of class BayesianModel
 */
class BayesianModel
{
  public:
    virtual ~BayesianModel() { };
    virtual char * name() = 0;
    virtual int numberOfVariables() = 0;

    //necessary methods for Gibbs Sampler
    virtual void drawVariable( int var_index ) = 0;
    virtual void fullConditionalUpdateVariable( int var_index ) = 0;
  
    //necessary method for Metropolis Hastings
    virtual void drawVariableFromProposal( int var_index ) = 0;
    virtual void updateVariableForProposal( int var_index ) = 0;
    virtual double logRatioTargetDensity( int var_index ) = 0;
    virtual double logRatioProposal( int var_index ) = 0;
    virtual void setCurrentVariableFromProposal( int var_index ) = 0;
    virtual void keepCurrentVariable( int var_index ) = 0;
 
    virtual void createOutput( int sim_to_keep ) = 0;
    virtual void initializeTemporaryStructures() = 0;
    virtual void dataAugmentationInitialDraws() = 0;
    virtual void keepSimulation( int simul_number ) = 0;
}; //end BayesianModel class

#endif /* BAYES_BAYESIAN_MODEL_H */
    
    
