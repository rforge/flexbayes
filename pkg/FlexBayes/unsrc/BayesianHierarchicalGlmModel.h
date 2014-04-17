#ifndef BAYESIAN_HIERARCHICAL_GLM_MODEL_H
#define BAYESIAN_HIERARCHICAL_GLM_MODEL_H

#include "rtErr.h"
#include "BayesianModel.h"
#include "BayesianGlmModel.h"
#include "DistributionExtended.h"
#include "DistributionMixture.h"

class BayesianHierarchicalGlmModel: public BayesianModel
{
  private:
  // data
    DistributionParameter ** random_predictors;
    DistributionParameter ** fixed_predictors;
    CMatrix ** beta_predictors;

  //distributions
    //**** link function ****
    //first stage
    //random effects
    NormalDistribution ** beta;
    //t prior
    InvChisqDistribution ** tau2_betas;

    //fixed effects
    NormalDistribution * gamma;
    InvChisqDistribution * tau2_gamma;

    //second stage
    NormalDistribution * alpha;
    DistributionMixture ** tau2;    
    InvChisqDistribution ** tau2ICS;
    ProperNonInfoPosteriorHLM ** tau2PNIP;
    double *tau2NIP;
    InvWishartDistribution * tau2InvWishart;
    InvChisqDistribution * tau2_alpha;

    DistributionParameter ** tau2_first_draw;

    int number_of_variables;
    int number_of_groups;
    int number_of_observations;
    int dim_beta;

    bool t_gamma;
    bool t_beta;
    bool t_alpha;
    bool non_informative_alpha;
    bool non_informative_beta;
    bool non_informative_gamma;
    bool invWishart_betaCov;

    bool keep_tau2_error;

    bool random_effects;
    bool fixed_effects;
    bool second_stage_effects;

    bool invChisq_tau2;
    bool duMouchel_tau2;
    bool uniformShrinkage_tau2;
    bool properNonInfoPrior_tau2;
    bool nonInformativePower_tau2;

    //connect to glm model
    BayesianGlmModel ** glm_model;

    //flag if proposal is made with update or constant Covariances
    bool update_for_hessians;

  //simulation samples
    int number_of_simulations;
    CMatrix ** simulated_beta;
    CMatrix * simulated_gamma;
    CMatrix * simulated_alpha;
    CMatrix ** simulated_tau2;
    CVector * simulated_tau2_gamma;
    CMatrix * simulated_tau2_beta;
    CVector * simulated_tau2_alpha;

  //working matrices for Metropolis Hastings sampler
    //second stage
    CMatrix ** zTz;
    CVector ** zTb;
    CMatrix * zTvIz;
    CVector * zTvIb;

    CVector ** second_residuals;

    //array of distributions map (maps distributions names to actual objects)
    char ** distr_map;

  public:
    //object declaration
    BayesianHierarchicalGlmModel();

    //the basic model normal-likelihood normal-beta-prior normal-alpha-prior invChisq-variances case 
    void initialize( CVector ** y, int n_groups );
    void randomEffects( CMatrix ** x );
    void fixedEffects( CMatrix ** m );
    void secondStageRandomEffects( CMatrix ** z );
    void glmModel( BayesianGlmModel * glm_vars );
    void setUpdateProposalCovariance( bool value );

    void betaPrior( CMatrix * p_betaCov );
    void betaPrior( CVector * p_beta, CMatrix * p_betaCov );
    void betaPrior( CMatrix * p_beta, CMatrix * p_betaCov );
    void betaPriorNonInformative( int dim ) throw( rtErr );

    void gammaPrior( CVector * p_gamma, CMatrix * p_gammaCov );
    void gammaPriorNonInformative( int dim );

    void alphaPrior( CVector * p_alpha, CMatrix * p_alphaCov ) throw( rtErr );
    void alphaPriorNonInformative( int dim ) throw( rtErr );

    void tau2PriorInvChisq( double *p_nuTau2, double *p_tau2 ) throw( rtErr );
    void tau2PriorDuMouchel( double *p_tau2 ) throw( rtErr );
    void tau2PriorUniformShrinkage( double *p_tau2 ) throw( rtErr );
    void tau2PriorNonInformative( double *p_power ) throw( rtErr );
    void betaCovPriorInvWishart( double p_nuV, CMatrix * p_VCov ) throw( rtErr );

    virtual ~BayesianHierarchicalGlmModel();
    void emptyModel();

    //the t-gamma-prior case
    void gammaTPrior( double p_nuGamma ) throw( rtErr );
    //the t-beta-prior case
    void betaTPrior( double p_nuBeta ) throw( rtErr );
    //the t-alpha-prior case
    void alphaTPrior( double p_nuAlpha ) throw( rtErr );

    inline bool IsAlphaPriorNonInformative() { return non_informative_alpha; }
    inline bool IsGammaT() { return ( t_gamma ); }
    inline bool IsBetaT() { return ( t_beta ); }
    inline bool IsAlphaT() { return ( t_alpha ); }
    inline bool IsBetaCovInvWishart() { return ( invWishart_betaCov ); }

    //initialize sampler
    void samplerDefaultInitialPoint();
    void samplerBetaInitialPoint( CVector & init_beta );
    void samplerBetaInitialPoint( CMatrix & init_beta );
    void samplerGammaInitialPoint( CVector & init_gamma );
    void samplerAlphaInitialPoint( CVector & init_alpha );
    void samplerTau2InitialPoint( double * init_tau2 );
    void samplerTau2InitialPoint( CMatrix & init_tau2 );

    void updateLinearResponse( int index );
    void updateLinearResponse();
    CMatrix betaInverseCovariance( int index );

    void gibbsUpdateSecondStageWorkingMatrices();
    void gibbsUpdateSecondStageResiduals();
    void gibbsUpdateSecondStageResiduals( int index );

    void gibbsUpdateBeta( int index );
    void gibbsUpdateGamma();

    void metropolisHastingsUpdateBeta( int index );
    void metropolisHastingsUpdateGamma();
    void gibbsUpdateAlpha();
    void gibbsUpdateTau2();
    void gibbsUpdateTau2Gamma();
    void gibbsUpdateTau2Beta( int index );
    void gibbsUpdateTau2Alpha();



    //output
    void simulationsToArray( double * simul_output, int simulations_to_keep ); 
    void simulationsToArrayWithMeans( double * simul_output, int simulations_to_keep );

    //required methods
    inline char * name() { return "Bayes Hierarchical Generalized Linear Model"; }
    inline int numberOfVariables() { return number_of_variables; }

    //required but not used
    void drawVariable( int var_index ) { };
    void fullConditionalUpdateVariable( int var_index ) { };

    //required methods
    void createOutput( int sim_to_keep );
    void initializeTemporaryStructures();
    void dataAugmentationInitialDraws();
    void keepSimulation( int simul_number );

    //required (for Metropolis Hastings sampler)
    void drawVariableFromProposal( int var_index );
    void updateVariableForProposal( int var_index );
    double logRatioTargetDensity( int var_index );
    double logRatioProposal( int var_index );
    void setCurrentVariableFromProposal( int var_index );
    void keepCurrentVariable( int var_index );

}; //end BayesianHierarchicalGlmModel


#endif /* BAYESIAN_HIERARCHICAL_GLM_MODEL_H */
