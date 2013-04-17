#ifndef BAYES_REJECTION_SAMPLING_DISTRIBUTION_H
#define BAYES_REJECTION_SAMPLING_DISTRIBUTION_H

#include "DistributionParameter.h"

/*
  abstract class
  to be used with RejectionSampling class
*/
class RejectionSamplingDistribution
{
  public:
    virtual ~RejectionSamplingDistribution() { };
    virtual double ratioOfDensities( DistributionParameter & par ) = 0;
    virtual DistributionParameter draw() = 0;
}; //end RejectionSamplingDensities class

#endif /* BAYES_REJECTION_SAMPLING_DISTRIBUTION_H */
