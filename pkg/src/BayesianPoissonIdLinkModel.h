#ifndef BAYESIAN_POISSON_ID_LINK_MODEL_H
#define BAYESIAN_POISSON_ID_LINK_MODEL_H

#include "BayesianGlmModel.h"
#include "DistributionExtended.h"
#include "DistributionMixture.h"

class BayesianPoissonIdLinkModel: public BayesianGlmModel
{
  private:
    CVector ** exposures;
    CVector ** counts;
    CVector ** mu_linear;
    CVector ** residuals;
    DistributionParameter *** log_lambda_previous;
    DistributionParameter ** sigma2_first_draw;

    //Poisson model parameters
    NormalDistribution *** log_lambda;

    DistributionMixture * sigma2;    
    InvChisqDistribution ** sigma2ICS;
    ProperNonInfoPosteriorHLM ** sigma2PNIP;
    double * sigma2NIP;

    bool invChisq_sigma2;
    bool duMouchel_sigma2;
    bool uniformShrinkage_sigma2;
    bool properNonInfoPrior_sigma2;
    bool nonInformativePower_sigma2;
    bool known_sigma2;
    bool common_sigma;


    int beta_dim;
    int gamma_dim;

    int number_of_observations;
    int number_of_variables;
    int number_of_groups;

    //simulation samples
    int number_of_simulations;
    CVector *** simulated_lambda;
    CMatrix * simulated_sigma2;

    bool update_hessians;
    bool gibbs_for_coefs;
    //array of distributions map (maps distributions names to actual objects)
    char ** distr_map;
    char * model_type;
  public:
    BayesianPoissonIdLinkModel();

    void initialize( CVector ** count_v, CVector ** expos, int n_groups );
    void initialize( CVector ** count_v, int n_groups );
    void setModelType( char * type );

    virtual ~BayesianPoissonIdLinkModel();
    void emptyModel();

    void setLambdaPrior();
    void sigma2CommonPrior( bool val ) { common_sigma = val; }
    void sigma2PriorInvChisq( double p_nuSigma2, double p_sigma2 );
    void sigma2PriorInvChisq( double * p_nuSigma2, double * p_sigma2 );
    void sigma2PriorDuMouchel( double p_sigma2 );
    void sigma2PriorDuMouchel( double * p_sigma2 );
    void sigma2PriorUniformShrinkage( double p_sigma2 );
    void sigma2PriorUniformShrinkage( double * p_sigma2 );
    void sigma2PriorNonInformative( double p_power ) throw( rtErr );
    void sigma2PriorNonInformative( double * p_power ) throw( rtErr );
    void sigma2Known( double * p_sigma2 );
    void sigma2Known( double p_sigma2 );

    void samplerDefaultInitialPoint();
    void samplerSigma2InitialPoint( double init_sigma2 );
    void samplerSigma2InitialPoint( CVector & init_sigma2 );

    void updateResiduals( int j, int i );
    void updateResiduals( int j );
    void updateResiduals();

    //require methods for link with mixed effects model
    void initializeCountStorage() { };
    void considerCounts() { };
    void updateCounts( int k ) { };
    void setCoefficientDimension( int p_random, int p_fixed );
    void metropolisHastingsUpdateModel( DistributionParameter * par_vals );
    CMatrix metropolisHastingsHessian( int n_pars, DistributionParameter ** par_vals );

    //not used
    CMatrix metropolisHastingsProposalHessian( int n_pars, DistributionParameter ** par_vals, DistributionParameter ** pars_vector );

    void setUpdateHessians( bool value ) { update_hessians = value; }
    CVector metropolisHastingsMean( int n_pars, DistributionParameter ** par_vals, DistributionParameter ** par_vector );


    inline bool gibbsDrawingsForCoefficients() { return gibbs_for_coefs; }
    void gibbsUpdateSigma2() throw( rtErr );
    void gibbsUpdateSigma2( int index ) throw( rtErr );
    void metropolisHastingsUpdateModel( int i, CVector & mu );
    void metropolisHastingsUpdateModel( int i );
    void metropolisHastingsUpdateModel( int j, int i );
    void metropolisHastingsUpdateLogLambda( int g, int obs );

    //output (required)
    void simulationsToArray( double * simul_output, int simulations_to_keep, int start_index ); 
    void simulationsToArray( double * simul_output, int start_index, CMatrix * simulated_mu ) { };

    //required methods
    inline char * name() { return "Bayes Hierarchical Poisson Log-Normal Model"; }
    inline int numberOfVariables() { return number_of_variables; }

    //required but not used
    void drawVariable( int var_index ) { };
    void fullConditionalUpdateVariable( int var_index ) { };


    void createOutput( int sim_to_keep );
    void initializeTemporaryStructures();
    void dataAugmentationInitialDraws();
    void keepSimulation( int simul_number );

    //required (for Metropolis Hastings sampler)
    void drawVariableFromProposal( int var_index );
    void updateVariableForProposal( int var_index );
    double logRatioTargetDensity( int var_index );
    inline double logRatioTargetDensity( int n_pars, DistributionParameter ** pars_pred, DistributionParameter ** pars_vector ) {return 0; }
    double logRatioProposal( int var_index );
    void setCurrentVariableFromProposal( int var_index );
    void keepCurrentVariable( int var_index );

    bool needsToUpdateVariable( char * var_name, DistributionParameter ** pars_test ) { return true; };

}; //end BayesianPoissonIdLinkModel


#endif /* BAYESIAN_POISSON_ID_LINK_MODEL_H */
