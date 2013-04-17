
#ifndef BAYES_REJECTION_SAMPLING_DISTRIBUTION_EXTENDED_H
#define BAYES_REJECTION_SAMPLING_DISTRIBUTION_EXTENDED_H

#include "Distribution.h"
#include "DistributionParameter.h"
#include "DistributionExtended.h"
#include "RejectionSamplingDistribution.h"

class CInvChisqByRejSampling : public RejectionSamplingDistribution
{
  private:
    CVector * interval;
    bool left_bounded_only;
    bool right_bounded_only;
    InvChisqDistribution * inv_chisq;
  public:
    CInvChisqByRejSampling();
    CInvChisqByRejSampling( double nu0 );
    virtual ~CInvChisqByRejSampling();
    CInvChisqByRejSampling( const CInvChisqByRejSampling & distr );
    CInvChisqByRejSampling & operator=( const CInvChisqByRejSampling & distr );
    void clean();
    inline InvChisqDistribution * boundDensity() { return inv_chisq; }
    void setInterval( double lower, double upper );
    void setInterval( double bound, char * type );
    double upperBound();
    double lowerBound();

    //required methods
    double ratioOfDensities( DistributionParameter & val );
    DistributionParameter draw();
}; //end class



/* 
  implementing Uniform Shrinkage for Bayesian Hierarchical Linear Model
*/

class ProperNonInfoPosteriorForBHLM : public Distribution, public RejectionSamplingDistribution
{
  private:
    double scale_var;
    double dfreedom_a;
    double dfreedom_b;
    double tau02;
    double tau0;
    double const_tau0_ratio_gammas;
    double ratio_gammas_b_over_a;
    double ratio_probs_b_over_a;
    double ratio_a_over_b;
    CVector * min_ratio;
    ConstrainedInvChisqDistribution * inv_chisq_a;
    ConstrainedInvChisqDistribution * inv_chisq_b;
    char * distr_type;
    double last_item_drawn;
  public:
    ProperNonInfoPosteriorForBHLM();
    ProperNonInfoPosteriorForBHLM( double df_a, double df_b, double tau_0_2 );
    virtual ~ProperNonInfoPosteriorForBHLM(); 
    ProperNonInfoPosteriorForBHLM( const ProperNonInfoPosteriorForBHLM & distr );
    ProperNonInfoPosteriorForBHLM & operator=( const ProperNonInfoPosteriorForBHLM & distr );
    Distribution * clone();
    void clean();

    inline double degreesOfFreedomBoundA() { return dfreedom_a; }
    inline double degreesOfFreedomBoundB() { return dfreedom_b; }
    inline double scale() { return scale_var; }
    inline double initialScaleBound() { return tau02; }
    inline char * priorDistribution() { return distr_type; }

    void setScale( double val );
 
    void initialize( double df_a, double df_b, double tau_0_2 );
    void setPriorDistribution( char * type );
    double ratioOfDensities( DistributionParameter & par );
    DistributionParameter draw();

    double logDensity( DistributionParameter & value );
    void update( double add_df, double add_scale );
    void update( DistributionParameter * par_list, int list_size );
    inline double lastItemDrawn() { return last_item_drawn; }
    DistributionParameter lastDraw();
    void setLastDraw( DistributionParameter & par );
    void setLastDraw( double par );
    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance();
}; //end class


#endif /* BAYES_REJECTION_SAMPLING_DISTRIBUTION_EXTENDED_H */
