
#ifndef BAYES_REJECTION_SAMPLING_H
#define BAYES_REJECTION_SAMPLING_H

#include "RejectionSamplingDistribution.h"
#include "DistributionParameter.h"

class RejectionSampling
{
  private:
    int number_of_samples;
    int number_drawings;
    bool save_stats;
    int acceptance_stats;
  public:
    RejectionSampling();
    RejectionSampling( int number_samples );
    RejectionSampling( const RejectionSampling & rej );
    RejectionSampling & operator=( const RejectionSampling & rej );
    virtual ~RejectionSampling();
    DistributionParameter ** drawSamples( RejectionSamplingDistribution * distr );
    DistributionParameter ** drawSamples( int n_samples, RejectionSamplingDistribution * distr );
    DistributionParameter drawASample( RejectionSamplingDistribution * distr );
    inline void saveAcceptanceRateStatistics( bool val ) { acceptance_stats = val; }
    inline int numberOfSamples() { return number_of_samples; }
    DistributionParameter AcceptanceRate();
    
}; //end RejectionSampling class


#endif /* BAYES_REJECTION_SAMPLING_H */
