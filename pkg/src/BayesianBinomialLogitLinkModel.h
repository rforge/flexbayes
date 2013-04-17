#ifndef BAYESIAN_BINOMIAL_LOGIT_LINK_MODEL_H
#define BAYESIAN_BINOMIAL_LOGIT_LINK_MODEL_H

#include "BayesianGlmModel.h"
#include "DistributionExtended.h"

class BayesianBinomialLogitLinkModel: public BayesianGlmModel
{
  private:
    CVector ** trials;
    CVector ** counts;
    CVector ** mu_linear;
    DistributionParameter ** xi_first_draw;

    //Poisson model parameters
    NormalDistribution ** log_xi;
    double * xi_z0;
    int xi_type;  // 0: multiple xi parameters.  
                  // 1: common xi parameter
                  // 2: no overdispersion
    
    int beta_dim;
    int gamma_dim;

    int number_of_observations;
    int number_of_variables;
    int number_of_groups;

    //working structures
    CVector * max_trials;
    CVector ** n_equal;
    CVector ** n_larger_equal;

    CVector * log_xi_hessian;
    CVector ** b1;
    CVector ** b2;
    CVector ** d;

    //simulation samples
    int number_of_simulations;
    CMatrix * simulated_xi;
    CMatrix * simulated_theta;

    CVector * computed_starting_hessian_logXi;
    CVector * computed_starting_hessian_beta;
    int computed_starting_hessian_gamma;

    CMatrix ** hessian_beta;
    CMatrix * hessian_gamma;
    bool update_hessians;

    bool gibbs_for_coefs;

    CVector * computed_starting_logXi_variance;
    CVector * logXi_variance;

  public:
    BayesianBinomialLogitLinkModel();

    void initialize( CVector ** count_v, CVector ** n_trials, int n_groups );
    void xiType( int val ) { if( (val >=0) && (val < 3) )  xi_type = val; }   
    	                                          // 0: multiple xi parameters.  
                                                // 1: common xi parameter
                                                // 2: no overdispersion
    // xiType should be called before logXiPrior, to specify the number of xi parameters
    // If xi_type == 2, then the following prior specification functions do not do 
    // anything since there are no overdispersion parameters
    void logXiPrior( double p_xi );
    void logXiPrior( double * p_xi );
    
    virtual ~BayesianBinomialLogitLinkModel();
    void emptyModel();

    // if xi_type == 2, then the following initial point functions do not do 
    // anything since there are no overdispersion parameters
    void samplerDefaultInitialPoint();
    void samplerXiInitialPoint( double init_xi );
    void samplerXiInitialPoint( CVector & init_xi );

    //require methods for link with mixed effects model
    void initializeCountStorage();
    void considerCounts();
    void updateCounts( int k ) { };
    void setCoefficientDimension( int p_random, int p_fixed );
    void metropolisHastingsUpdateModel( DistributionParameter * par_vals );
    CMatrix metropolisHastingsHessian( int n_pars, DistributionParameter ** par_vals );
    CMatrix metropolisHastingsProposalHessian( int n_pars, DistributionParameter ** par_vals, DistributionParameter ** pars_vector );
    //not used
    CVector metropolisHastingsMean( int n_pars, DistributionParameter ** par_vals, DistributionParameter ** par_vector );

    CVector computeWeights( int var_index, DistributionParameter ** par_vals, CVector & b_star, CVector & b_0 );

    inline bool gibbsDrawingsForCoefficients() { return gibbs_for_coefs; }
    void metropolisHastingsUpdateModel( int i, CVector & mu );
    void metropolisHastingsUpdateModel( int i );
    double logXiGradient( int j, double xi );
    void computeHessianForLogXi( int j, double xi );
    double hessianForlogXi( int i );
    void metropolisHastingsUpdateLogXi( int i );    
    void setUpdateHessians( bool value ) { update_hessians = value; }
    void generateThetas( CMatrix * simulated_mu );

    //output (required)
    void simulationsToArray( double * simul_output, int simulations_to_keep, int start_index ) { }; 
    void simulationsToArray( double * simul_output, int start_index, CMatrix * simulated_mu );

    //required methods
    inline char * name() { return "Bayes Hierarchical Binomial-Beta Conjugate Model"; }
    inline int numberOfVariables() { return number_of_variables; }

    //required but not used
    void drawVariable( int var_index ) { };
    void fullConditionalUpdateVariable( int var_index ) { };


    void createOutput( int sim_to_keep );
    void initializeTemporaryStructures();
    void dataAugmentationInitialDraws();
    void keepSimulation( int simul_number );

    // The following functions perform calculations for the Metropolis Hastings sampler:
    void drawVariableFromProposal( int var_index );
    void updateVariableForProposal( int var_index );
    // Calculates the log ratio target density for a proposed change in xi:
    double logRatioTargetDensity( int var_index );
    // Calculates the log likelihood ratio for a proposed change in beta (the random effects)
    // or gamma (the fixed effects), for the case where there is overdispersion:
    double logRatioTargetDensity( DistributionParameter ** pars_predictors, CVector & b_star, CVector & b_0 );
    // Calculates the log likelihood ratio for a proposed change in beta (the random effects)
    // or gamma (the fixed effects), for the case where there is no overdispersion:
    double logRatioTargetNoOverdisperse( DistributionParameter ** pars_predictors, CVector & b_star, CVector & b_0 );
    // A wrapper for the previous two functions:
    double logRatioTargetDensity( int n_pars, DistributionParameter ** pars_pred, 
      DistributionParameter ** pars_vector );
    double logRatioProposal( int var_index );
    void setCurrentVariableFromProposal( int var_index );
    void keepCurrentVariable( int var_index );



    bool needsToUpdateVariable( char * var_name, DistributionParameter ** pars_test ) { return true; };

}; //end BayesianBinomialLogitLinkModel


#endif /* BAYESIAN_BINOMIAL_LOGIT_LINK_MODEL_H */
