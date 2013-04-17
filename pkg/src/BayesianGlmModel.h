
#ifndef BAYES_BAYESIAN_GLM_MODEL_H
#define BAYES_BAYESIAN_GLM_MODEL_H

#include "DistributionParameter.h"

/* abstract declaration of class BayesianModel
 */
class BayesianGlmModel
{
  public:
    virtual ~BayesianGlmModel() { };
    virtual char * name() = 0;
    virtual int numberOfVariables() = 0;

    virtual bool gibbsDrawingsForCoefficients() = 0;
    //necessary methods for Gibbs Sampler
    virtual void drawVariable( int var_index ) = 0;
    virtual void fullConditionalUpdateVariable( int var_index ) = 0;
  
    //necessary method for Metropolis Hastings
    virtual void drawVariableFromProposal( int var_index ) = 0;
    virtual void updateVariableForProposal( int var_index ) = 0;
    virtual double logRatioTargetDensity( int var_index ) = 0;
    virtual double logRatioTargetDensity( int n_pars, DistributionParameter ** pars_predictors, DistributionParameter ** pars_vector ) = 0;
    virtual double logRatioProposal( int var_index ) = 0;
    virtual void setCurrentVariableFromProposal( int var_index ) = 0;
    virtual void keepCurrentVariable( int var_index ) = 0;
    virtual void setUpdateHessians( bool value ) = 0;
 
    virtual void createOutput( int sim_to_keep ) = 0;
    virtual void initializeTemporaryStructures() = 0;
    virtual void dataAugmentationInitialDraws() = 0;
    virtual void keepSimulation( int simul_number ) = 0;

    virtual void simulationsToArray( double * simul_output, int simulations_to_keep, int start_index ) = 0; 
    virtual void simulationsToArray( double * simul_output, int start_index, CMatrix * simulated_mu ) = 0;

    virtual void initializeCountStorage() = 0;
    virtual void considerCounts() = 0;
    virtual void updateCounts( int k ) = 0;
    virtual void setCoefficientDimension( int p_random, int p_fixed ) = 0;
    //required for link
    virtual void metropolisHastingsUpdateModel( DistributionParameter * par_vals ) = 0;
    virtual CMatrix metropolisHastingsHessian( int n_pars, DistributionParameter ** par_vals ) = 0;
    virtual CMatrix metropolisHastingsProposalHessian( int n_pars, DistributionParameter ** par_vals, DistributionParameter ** pars_vector ) = 0;
    virtual CVector metropolisHastingsMean( int n_pars, DistributionParameter ** par_vals, DistributionParameter ** par_vector ) = 0;

    virtual bool needsToUpdateVariable( char * var_name, DistributionParameter ** par_vals ) = 0;

}; //end BayesianGlmModel class

#endif /* BAYES_BAYESIAN_GLM_MODEL_H */
    
    
