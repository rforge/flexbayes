#ifndef BAYESIAN_MISSING_DATA_HLM_H
#define BAYESIAN_MISSING_DATA_HLM_H

#include "rtErr.h"
#include "BayesianModel.h"
#include "DistributionExtended.h"
#include "DistributionMixture.h"
#include "RejectionSamplingDistributionExtended.h"

class BayesianMissingDataHLM: public BayesianModel
{
  private:
  // data
    CMatrix ** response;
    
  //distributions
    //first stage
    //random effects
    NormalDistribution ** beta;
    //t prior
    InvChisqDistribution ** tau2_betas;

    //fixed effects
    NormalDistribution * gamma;
    InvChisqDistribution * tau2_gamma;

    DistributionMixture * sigma2;    
    InvChisqDistribution ** sigma2ICS;
    ProperNonInfoPosteriorForBHLM ** sigma2PNIP;
    double * sigma2NIP;
    InvWishartDistribution ** sigma2InvWishart;

    //second stage
    NormalDistribution * alpha;
    DistributionMixture * tau2;    
    InvChisqDistribution * tau2ICS;
    ProperNonInfoPosteriorForBHLM * tau2PNIP;
    double tau2NIP;
    InvWishartDistribution * tau2InvWishart;
    InvChisqDistribution * tau2_alpha;

    DistributionParameter ** sigma2_first_draw;
    DistributionParameter * tau2_first_draw;

    //t-likelihood
    InvChisqDistribution ** tau2_groups;
    InvChisqDistribution ** tau2_errors;

    bool keep_missing_data;
    CVector * number_missing_response;
    CVector *** missing_response_components;

    int number_of_variables;
    int number_of_groups;
    int number_of_observations;

    bool t_lkhd;
    bool t_gamma;
    bool t_beta;
    bool t_alpha;
    bool t_group_lkhd;
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
    bool invWishart_sigma2;
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
    CMatrix *** simulated_sigma2;
    CMatrix ** simulated_tau2;
    CVector * simulated_tau2_gamma;
    CMatrix * simulated_tau2_beta;
    CVector * simulated_tau2_alpha;
    CMatrix * simulated_tau2_errors;
    CMatrix * simulated_tau2_groups;
    CMatrix ** simulated_missingR;

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
    CMatrix ** residuals;
    CVector ** second_residuals;


    //array of distributions map (maps distributions gs to actual objects)
    char ** distr_map;

  public:
    //object declaration
    BayesianMissingDataHLM();

    //the basic model normal-likelihood normal-beta-prior normal-alpha-prior invChisq-variances case 
    void initialize( CMatrix ** y, int n_groups );

    void randomEffects();
    void fixedEffects();
    void secondStageRandomEffects();

    void initializeResponseMissingData( bool keep_them );
    void responseMissingData( int index, CVector & m_comp );

    void betaPrior( CMatrix * p_betaCov );
    void betaPrior( CVector * p_beta, CMatrix * p_betaCov );
    void betaPrior( CMatrix * p_beta, CMatrix * p_betaCov );

    void gammaPrior( CVector * p_gamma, CMatrix * p_gammaCov );
    void gammaPriorNonInformative( int dim );

    void alphaPrior( CVector * p_alpha, CMatrix * p_alphaCov ) throw( rtErr );
    void alphaPriorNonInformative( int dim ) throw( rtErr );

    void sigma2CommonPrior( bool val );

    void sigma2PriorInvChisq( double p_nuSigma2, double p_sigma2 );
    void sigma2PriorDuMouchel( double p_sigma2 );
    void sigma2PriorUniformShrinkage( double p_sigma2 );
    void sigma2PriorNonInformative( double p_power ) throw( rtErr );
    void sigma2PriorInvWishart( double p_nu, CMatrix * p_Cov );
    void sigma2Known( CMatrix * p_Cov );

    void sigma2PriorInvChisq( double * p_nuSigma2, double * p_sigma2 );
    void sigma2PriorDuMouchel( double * p_sigma2 );
    void sigma2PriorUniformShrinkage( double * p_sigma2 );
    void sigma2PriorNonInformative( double * p_power ) throw( rtErr );
    void sigma2PriorInvWishart( double * p_nu, CMatrix ** p_Cov );
    void sigma2Known( CMatrix ** p_Cov );

    void tau2PriorInvChisq( double p_nuTau2, double p_tau2 ) throw( rtErr );
    void tau2PriorDuMouchel( double p_tau2 ) throw( rtErr );
    void tau2PriorUniformShrinkage( double p_tau2 ) throw( rtErr );
    void tau2PriorNonInformative( double p_power ) throw( rtErr );
    void betaCovPriorInvWishart( double p_nuV, CMatrix * p_VCov ) throw( rtErr );

    virtual ~BayesianMissingDataHLM();
    void emptyModel();

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
    void samplerSigma2InitialPoint( CMatrix & init_sigma2 );
    void samplerSigma2InitialPoint( CMatrix ** init_sigma2 );
    void samplerBetaInitialPoint( CVector & init_beta );
    void samplerBetaInitialPoint( CMatrix & init_beta );
    void samplerGammaInitialPoint( CVector & init_gamma );
    void samplerAlphaInitialPoint( CVector & init_alpha );
    void samplerTau2InitialPoint( double init_tau2 );
    void samplerTau2InitialPoint( CMatrix & init_tau2 );
  
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
    //inline char * name() { return "Bayes Hierarchical Linear Model"; }
    inline int numberOfVariables() { return number_of_variables; }
    void drawVariable( int var_index );
    void drawMissingResponse( int index );
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

}; //end BayesianMissingDataHLM


#endif /* BAYESIAN_MISSING_DATA_HLM_H */
