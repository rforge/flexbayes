#ifndef BAYESIAN_LINEAR_MODEL_H
#define BAYESIAN_LINEAR_MODEL_H

#include "BayesianModel.h"
#include "DistributionExtended.h"
#include "DistributionMixture.h"

class BayesianLinearModel: public BayesianModel
{
  private:
  // data
    CVector * response;
    CMatrix * predictors;
    
  //distributions
    InvChisqDistribution * sigma2;
    //mixture
    DistributionMixture * sigma2_mix;
    InvChisqDistribution ** sigma2s;

    NormalDistribution * beta;
    //mixture
    bool mixture_prior;
    DistributionMixture * beta_mix;
    NormalDistribution ** betas;
    StudentTDistribution ** betas_t;
    //augmented data 
    InvChisqDistribution * tau2;
    InvChisqDistribution ** tau2_error;   

    int number_of_variables;

    bool t_lkhd;
    bool t_beta;
    bool t_lkhd_beta;
    bool keep_tau2_error;
    bool conjugate_model;
    bool non_informative;
    bool simulations_performed;

  //simulation samples
    CMatrix * simulated_beta;
    CVector * simulated_sigma2;
    CMatrix * simulated_tau2_errors;

  //working matrices for Gibbs sampler
    CMatrix * xTx;
    CVector * xTy;
    CVector * tau2_error_weights;
    CVector * residuals;

    //array of distributions map (maps distributions names to actual objects)
    char ** distr_map;

  public:
    //object declaration
    BayesianLinearModel();
    //the basic model normal-likelihood normal-beta-prior case (??)
    void initialize( CVector * y, CMatrix * x, 
                     CVector * p_beta, CMatrix * p_betaCov,
                     double p_sigma2, double p_nuSigma2 ); 
    //initialization for beta prior consisting of mixture of normals or t's
    void initialize( CVector * y, CMatrix * x, 
                     int number_mixtures, CVector ** p_beta, CMatrix ** p_betaCov,
                     double * gamma, double * d_freedom, char * mixture_type,
                     double p_sigma2, double p_nuSigma, int save_stats );
    //setting up a non-informative beta prior. 
    //Assumes p(beta, sigma2) proportional to sigma2^{-1} (i.e. InvChisq( df = 0 ) )
    //However it accepts any InvChisq prior for sigma2
    void nonInformative( CVector * y, CMatrix * x,
                         double p_sigma2, double p_nuSigma2 ); 
    virtual ~BayesianLinearModel();
    void emptyModel();
    //the t-likelihood normal-beta-prior case
    void t_likelihood( double p_nuError );
    //the normal-likelihood t-beta-prior case
    void t_beta_prior( double p_nuBeta );
    //the t-likelihood t-beta-prior case
    void t_likelihood_and_beta_prior( double p_nuError, double p_nuBeta );
    void conjugatePrior();
    inline bool IsPriorConjugate() { return conjugate_model; }
    inline bool IsPriorNonInformative() { return non_informative; }
    inline bool tLikelihood() { return ( t_lkhd || t_lkhd_beta ); }


    //mixture case
    void updateMixtureProportion( int i );
    void updateSigma2MixtureProportion( int i, double s_beta, double add_df, bool include_Cov );

    //initialize sampler
    void samplerDefaultInitialPoint();
    void samplerBetaInitialPoint( CVector & init_beta );
    void samplerSigma2InitialPoint( double init_sigma2 );
  
    void updateRegressionWeight( int index );
    void gibbsUpdateWorkingMatrices();
    void gibbsUpdateBeta();
    void gibbsUpdateSigma2();
    void gibbsUpdateTau2();
    void gibbsUpdateTau2Error( int index );
    
    //exact sampling methods
    void exactSamplerPosteriorSigma2();
    void exactSamplerStep();
    void exactSampler( int simul_to_keep );

    //output
    double * simulationsToArray( int simulations_to_keep ); 
    void simulationsToArray( double * simul_ouput, int simulations_to_keep ); 
    void saveMixtureDrawingStatistics( double * stats );

    //required methods
    //inline char * name() { return "Bayes Linear Model"; }
    inline int numberOfVariables() { return number_of_variables; }
    void drawVariable( int var_index );
    void fullConditionalUpdateVariable( int var_index );
    void createOutput( int sim_to_keep );
    void initializeTemporaryStructures();
    void dataAugmentationInitialDraws();
    void keepSimulation( int simul_number );

    //required but not used methods
    void drawVariableFromProposal( int var_index ){ };
    void updateVariableForProposal( int var_index ){ };
    double logRatioTargetDensity( int var_index ) { return 0; }
    double logRatioProposal( int var_index ) { return 0; }
    void setCurrentVariableFromProposal( int var_index ) { };
    void keepCurrentVariable( int var_index ) { };    

}; //end BayesianLinearModel


#endif /* BAYESIAN_LINEAR_MODEL_H */
