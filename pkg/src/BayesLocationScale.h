

class BayesianLinearModel: public BayesianModel
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
    CVector * simulated_beta;
    CVector * simulated_sigma2;

  //working matrices for Gibbs sampler
    CMatrix * xTx;
    CVector * xTy;
    CVector * tau2_error_weights;
    CVector * residuals;
