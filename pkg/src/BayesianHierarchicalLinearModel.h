#ifndef BAYESIAN_HIERARCHICAL_LINEAR_MODEL_H
#define BAYESIAN_HIERARCHICAL_LINEAR_MODEL_H

#include "rtErr.h"
#include "BayesianModel.h"
#include "DistributionExtended.h"
#include "DistributionMixture.h"

class BayesianHierarchicalLinearModel: public BayesianModel
{
  private:
  // data
    CVector ** response;
    CMatrix ** random_predictors;
    CMatrix ** fixed_predictors;
    CMatrix ** beta_predictors;
    
  //full conditional distributions
    //first stage
    //random effects
    NormalDistribution ** beta;
    //t prior only:
    InvChisqDistribution ** tau2_betas;

    //fixed effects
    NormalDistribution * gamma;
    //t prior only:
    InvChisqDistribution * tau2_gamma;

    //outcome variance.  This is stored as a "mixture"
    // distribution with only a single mixture component.
    // This could be extended later to accomodate a mixture 
    // prior distribution.  It also makes the syntax easier
    // when one desires to obtain, e.g., the current value of
    // sigma2 without conditioning on the type of prior.
    DistributionMixture * sigma2;    
    InvChisqDistribution ** sigma2ICS;
    ProperNonInfoPosteriorHLM ** sigma2PNIP;
    double * sigma2NIP;

    //second stage
    NormalDistribution * alpha;
    DistributionMixture ** tau2;    
    InvChisqDistribution ** tau2ICS;
    ProperNonInfoPosteriorHLM ** tau2PNIP;
    double *tau2NIP;
    InvWishartDistribution * tau2InvWishart;
    InvChisqDistribution * tau2_alpha;

    DistributionParameter ** sigma2_first_draw;
    DistributionParameter ** tau2_first_draw;

    //t-likelihood
    InvChisqDistribution ** tau2_groups;
    InvChisqDistribution ** tau2_errors;


    bool keep_missing_data;
    CVector * number_missing_response;
    CVector *** missing_response_components;
    CVector * number_missing_random_predictors;
    CVector *** missing_random_predictors_components;
    CVector * number_missing_fixed_predictors;
    CVector *** missing_fixed_predictors_components;


    int number_of_variables;
    int number_of_groups;
    int number_of_observations;
    int dim_beta;

    bool t_lkhd;
    bool t_gamma;
    bool t_beta;
    bool t_alpha;
    bool t_group_lkhd;
    bool non_informative_beta;
    bool non_informative_alpha;
    bool non_informative_gamma;
    bool invWishart_betaCov;

    bool keep_tau2_error;

    bool random_effects;
    bool fixed_effects;
    bool second_stage_effects;

    bool common_sigma;
    bool invChisq_sigma2;
    bool duMouchel_sigma2;
    bool uniformShrinkage_sigma2;
    bool properNonInfoPrior_sigma2;
    bool nonInformativePower_sigma2;
    bool known_sigma2;

    bool invChisq_tau2;
    bool duMouchel_tau2;
    bool uniformShrinkage_tau2;
    bool properNonInfoPrior_tau2;
    bool nonInformativePower_tau2;

  //simulation samples
    int number_of_simulations;
    CMatrix ** simulated_beta;
    CMatrix * simulated_gamma;
    CMatrix * simulated_alpha;
    CMatrix * simulated_sigma2;
    CMatrix ** simulated_tau2;
    CVector * simulated_tau2_gamma;
    CMatrix * simulated_tau2_beta;
    CVector * simulated_tau2_alpha;
    CMatrix * simulated_tau2_errors;
    CMatrix * simulated_tau2_groups;
    CMatrix ** simulated_missingR;
    CMatrix ** simulated_missingRP;
    CMatrix ** simulated_missingFP;

  //working matrices for Gibbs sampler
    //random effects
    CMatrix ** xTx;
    CVector ** xTy;
    //fixed effects
    CMatrix * sMTm;
    CVector * sMTy;
    CMatrix ** mTm;
    CVector ** mTy;
    //second stage
    CMatrix ** zTz;
    CVector ** zTb;
    CMatrix * zTvIz;
    CVector * zTvIb;

    CVector * final_group_index;

    CVector ** tau2_error_weights;
    CVector * tau2_group_weights;
    CVector ** residuals;
    CVector ** second_residuals;


    //array of distributions map (maps distributions names to actual objects)
    char ** distr_map;


    //missing data priors
    CVector ** randomEffects_means;
    CVector ** randomEffects_vars;
    CVector ** fixedEffects_means;
    CVector ** fixedEffects_vars;

    // The constructor calls emptyModel,
    // which sets some of the internal 
    // constants to default values.
    void emptyModel();
    
  public:
    //constructor
    BayesianHierarchicalLinearModel();
    // Destructor
    virtual ~BayesianHierarchicalLinearModel();

    //Initialize the model by specifying the outcome data
    void initialize( CVector ** y, int n_groups );

    /* Specify the covariate data; some of this is 
    for the random effects, while some of it is for the
    fixed effects or for the second-level effects. */
    void randomEffects( CMatrix ** x );
    void fixedEffects( CMatrix ** m );
    void secondStageRandomEffects( CMatrix ** z );

    // Handle missing data
    void initializeResponseMissingData( bool keep_them );
    void responseMissingData( int index, CVector & m_comp );
    void initializeRandomPredictorsMissingData( bool keep_them );
    void randomPredictorsMissingData( int index, CVector & m_comp );
    void initializeFixedPredictorsMissingData( bool keep_them );
    void fixedPredictorsMissingData( int index, CVector & m_comp );

    /* Prior specification functions.  Call these functions only
    after specifying all of the data using the initialize,
    randomEffects, fixedEffects, and secondStageRandomEffects
    functions as appropriate */
    void betaPrior( CMatrix * p_betaCov );
    void betaPrior( CVector * p_beta, CMatrix * p_betaCov );
    void betaPrior( CMatrix * p_beta, CMatrix * p_betaCov );
    void betaPriorNonInformative( int dim );

    void gammaPrior( CVector * p_gamma, CMatrix * p_gammaCov );
    void gammaPriorNonInformative( int dim );

    void alphaPrior( CVector * p_alpha, CMatrix * p_alphaCov ) throw( rtErr );
    void alphaPriorNonInformative( int dim ) throw( rtErr );

    void sigma2CommonPrior( bool val );

    void sigma2PriorInvChisq( double p_nuSigma2, double p_sigma2 );
    void sigma2PriorDuMouchel( double p_sigma2 );
    void sigma2PriorUniformShrinkage( double p_sigma2 );
    void sigma2PriorNonInformative( double p_power ) throw( rtErr );

    void sigma2PriorInvChisq( double * p_nuSigma2, double * p_sigma2 );
    void sigma2PriorDuMouchel( double * p_sigma2 );
    void sigma2PriorUniformShrinkage( double * p_sigma2 );
    void sigma2PriorNonInformative( double * p_power ) throw( rtErr );
    void sigma2Known( double * p_sigma2 );
    void sigma2Known( double p_sigma2 );

    // Prior specification for tau2, the random effect variance parameters.
    // The prior hyperparameters for the Inverse Chisquared, DuMouchel, 
    // UniformShrinkage, and NonInformative priors must be vectors of length
    // equal to the number of random effect elements, i.e. the number of 
    // columns of x in the call to the randomEffects function.
    void tau2PriorInvChisq( double *p_nuTau2, double *p_tau2 ) throw( rtErr );
    void tau2PriorDuMouchel( double *p_tau2 ) throw( rtErr );
    void tau2PriorUniformShrinkage( double *p_tau2 ) throw( rtErr );
    void tau2PriorNonInformative( double *p_power ) throw( rtErr );
    void betaCovPriorInvWishart( double p_nuV, CMatrix * p_VCov ) throw( rtErr );


    //the t-likelihood 
    void tLikelihood( double p_nuError );
    void groupTLikelihood( double p_nuGroup );
    //the t-gamma-prior case
    void gammaTPrior( double p_nuGamma ) throw( rtErr );
    //the t-beta-prior case
    void betaTPrior( double p_nuBeta ) throw( rtErr );
    //the t-alpha-prior case
    void alphaTPrior( double p_nuAlpha ) throw( rtErr );

    inline bool IsAlphaPriorNonInformative() { return non_informative_alpha; }
    inline bool IsTLikelihood() { return ( t_lkhd ); }
    inline bool IsGroupTLikelihood() { return ( t_group_lkhd ); }
    inline bool IsGammaT() { return ( t_gamma ); }
    inline bool IsBetaT() { return ( t_beta ); }
    inline bool IsAlphaT() { return ( t_alpha ); }
    inline bool IsBetaCovInvWishart() { return ( invWishart_betaCov ); }

    //initialize sampler
    void samplerDefaultInitialPoint();
    void samplerSigma2InitialPoint( double init_sigma2 );
    void samplerSigma2InitialPoint( CVector & init_sigma2 );
    void samplerBetaInitialPoint( CVector & init_beta );
    void samplerBetaInitialPoint( CMatrix & init_beta );
    void samplerGammaInitialPoint( CVector & init_gamma );
    void samplerAlphaInitialPoint( CVector & init_alpha );
    void samplerTau2InitialPoint( double *init_tau2 );
    void samplerTau2InitialPoint( CMatrix & init_tau2 );

    void samplerMissingVariablesInitialPoint();
    void samplerMissingResponseInitialPoint();
    void samplerMissingRandomPredictorsInitialPoint();
    void samplerMissingFixedPredictorsInitialPoint();
  
    //retrieve group associated to observation index
    DistributionParameter groupIndex( int index );

    void updateRegressionWeight( int index ) throw( rtErr );
    void updateGroupRegressionWeight( int grp );
    void gibbsUpdateWorkingMatrices();
    void gibbsUpdateWorkingMatrices( int index, char * type );
    void gibbsUpdateSecondStageWorkingMatrices();
    void gibbsUpdateResiduals();
    void gibbsUpdateResiduals( int index );
    void gibbsUpdateSecondStageResiduals();
    void gibbsUpdateSecondStageResiduals( int index );
    void gibbsUpdateBeta( int index );
    void gibbsUpdateGamma();
    void gibbsUpdateAlpha();
    void gibbsUpdateSigma2() throw( rtErr );
    void gibbsUpdateSigma2( int index ) throw( rtErr );
    void gibbsUpdateTau2();
    void gibbsUpdateTau2Gamma();
    void gibbsUpdateTau2Beta( int index );
    void gibbsUpdateTau2Alpha();
    void gibbsUpdateTau2Group( int index );
    void gibbsUpdateTau2Error( int index );

    //output
    void simulationsToArray( double * simul_output, int simulations_to_keep ); 

    //required methods
    inline char * name() { return "Bayes Hierarchical Linear Model"; }
    inline int numberOfVariables() { return number_of_variables; }
    void drawVariable( int var_index );
    void drawMissingResponse( int index );
    void drawMissingRandomPredictors( int index );
    void drawMissingFixedPredictors( int index );
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

}; //end BayesianHierarchicalLinearModel


#endif /* BAYESIAN_HIERARCHICAL_LINEAR_MODEL_H */
