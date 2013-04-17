
#ifndef BAYES_BLM_H
#define BAYES_BLM_H

//#include "Matrix.h"
#include "distributions.h"

class BayesLinearModel
{
  private:
  // data
    CVector * response;
    CMatrix * predictors;

  //distributions
    InvChisqObject * sigma2;
    NormalObject * beta;
    //augmented data 
    InvChisqObject * tau2;
    InvChisqObject ** tau2_error;

    bool t_lkhd;
    bool t_beta;
    bool t_lkhd_beta;
    bool conjugate_model;
    bool non_informative;

  //simulation samples
    CMatrix * simulated_beta;
    CVector * simulated_sigma2;

  //Gibbs sampler parameters
    int burn_in_time;
    int simulations_to_keep;
    int sample_frequency;
    int simulations_performed;
  //working matrices for Gibbs sampler
    CMatrix * xTx;
    CVector * xTy;
    CVector * tau2_error_weights;
    CVector * residuals;
  public:
    //object declaration
    BayesLinearModel();
    //the basic model normal-likelihood normal-beta-prior case (??)
    void initialize( CVector * y, CMatrix * x, 
                     CVector * p_beta, CMatrix * p_betaCov,
                     double p_sigma2, double p_nuSigma2 ); 
    //setting up a non-informative beta prior. 
    //Assumes p(beta, sigma2) proportional to sigma2^{-1} (i.e. InvChisq( df = 0 ) )
    //However it accepts any InvChisq prior for sigma2
    void nonInformative( CVector * y, CMatrix * x,
                         double p_sigma2, double p_nuSigma2 ); 
    virtual ~BayesLinearModel();
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

    //initialize sampler
    void samplerDefaultInitialPoint();
    void samplerBetaInitialPoint( CVector & init_beta );
    void samplerSigma2InitialPoint( double init_sigma2 );

    //exact sampler
    void exactSampler( int simul_to_keep );
    void exactSamplerPosteriorSigma2();
    void exactSamplerStep();

    //Gibbs sampler
    void gibbsSamplerInitialize( int burn_in, int simul_to_keep, int sample_freq );
    void doGibbsSampler();
    void gibbsStep();
    bool gibbsSampleToKeep();
    bool gibbsSamplerIsDone( int count );  
    void gibbsUpdateWorkingMatrices();
    void gibbsUpdateBeta();
    void gibbsUpdateSigma2();
    void gibbsUpdateTau2();
    void gibbsUpdateTau2Error( int i );
    //output
    double * simulationsToArray(); 
    void  simulationsToArray( double * simul_ouput ); 

};//end class BayesLinearModel


#endif /* BAYES_BLM_H */
