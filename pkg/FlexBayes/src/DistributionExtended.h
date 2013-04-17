
#ifndef BAYES_DISTRIBUTION_EXTENDED_H
#define BAYES_DISTRIBUTION_EXTENDED_H

#include "rtErr.h"
#include "Matrix.h"
#include "Distribution.h"


class InvChisqDistribution : public Distribution
{
  private:
    double dfreedom_init;
    double scale_init;
    double dfreedom_updated;
    double scale_updated;
    double last_item_drawn;
  public:
    InvChisqDistribution();
    InvChisqDistribution( double nu0 );
    ~InvChisqDistribution();
    
    inline double degreesOfFreedom() { return dfreedom_updated; }
    inline double scale() { return scale_updated; }
    inline double initialDegreesOfFreedom() { return dfreedom_init; }
    inline double initialScale() { return scale_init; }

    void setDegreesOfFreedom( double nu0 );
    void setScale( double sigma02 ); //initial value of scale = sigma0 * sigma0

    void update( double additional_df, double additional_scale );
    void update( DistributionParameter * par_list, int list_size );

    InvChisqDistribution( const InvChisqDistribution & distr );
    Distribution * clone();
    InvChisqDistribution & operator =(const InvChisqDistribution & distr);

    double drawOneItem();
    DistributionParameter draw();

    inline double lastItemDrawn() { return last_item_drawn; }
    DistributionParameter lastDraw();

    void setLastDraw( double last_draw );
    void setLastDraw( DistributionParameter & last_draw );

    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance();
    double logDensity( DistributionParameter & value );
};//end class


class NormalDistribution : public Distribution
{
  private:
    int dim;
    bool non_informative;
    double inv_determinant_init;
    CVector * very_first_mean;
    CMatrix * very_first_cov;
    CMatrix * very_first_inv_cov;

    CVector * mean_init;
    CMatrix * cov_init;
    CMatrix * inv_cov_init;
    CVector * inv_cov_beta_init;
    CVector * mean_updated;
    CVector * cholesky_invCov_updated;
    CVector * last_item_drawn;
  public:
    NormalDistribution( int dim_vector );
    ~NormalDistribution();
    void clean();
    inline int dimension() { return dim; }

    inline CVector initialMean() { return *mean_init; }
    inline CMatrix initialCovariance() { return *cov_init; }
    inline CMatrix initialInverseCovariance() { return *inv_cov_init; }
    inline CVector initialInvCovBeta() { return *inv_cov_beta_init; }
    inline CVector choleskyInvCov() { return *cholesky_invCov_updated; }
    inline CVector meanVec() { return *mean_updated; }
    inline double initialInvDeterminant() { return inv_determinant_init; }

    void setNonInformative();
    void setMean( CVector * mean );
    void setCovariance( CMatrix * cov ) throw( rtErr );
    void setCovariance( CMatrix * cov, CVector * cholesky_invCov );
    void scaleLastDraw( double scale );
    void updateScale( double scale );

    void updateInitialMean( CVector * mean );
    void scaleInitialMean( double scaleMult );
    void updateInitialCovariance( CMatrix * cov ) throw( rtErr );
    void updateInitialCovariance( CMatrix * cov, CVector * cholesky_invCov );
    void updateVeryFirstCovariance( CMatrix & cov );
    inline CMatrix veryFirstCovariance() { return *very_first_cov; }
    inline CMatrix veryFirstInverseCovariance() { return *very_first_inv_cov; }
    inline CVector veryFirstMean() { return *very_first_mean; }

    NormalDistribution( const NormalDistribution & distr );
    Distribution * clone();
    NormalDistribution & operator =(const NormalDistribution & rhs);

    DistributionParameter mean();
    CMatrix covariance();
    double invCovDeterminant();
    DistributionParameter variance();
    void update( CVector & mean_shift, CMatrix & cov_shift );
    void update( double precision_tau2, CVector & mean_shift, CMatrix & cov_shift );
    void update( CVector & precision_tau2, CVector & mean_shift, CMatrix & cov_shift );
    void update( DistributionParameter * par_list, int list_size );
    CVector drawOneItem();
    CVector drawItemFromPrior();
    DistributionParameter draw();
    void setLastDraw( const CVector& last_beta );
    void setLastDraw( DistributionParameter & last_beta );
    inline CVector lastItemDrawn() { return (*last_item_drawn); }
    DistributionParameter lastDraw();
    virtual DistributionParameter mode();
    virtual double logDensity( DistributionParameter & value );
};//end class



class StudentTDistribution : public Distribution
{
  private:
    InvChisqDistribution * inv_chisq_dist;
    NormalDistribution * normal_dist;
  public:
    StudentTDistribution();
    StudentTDistribution( int dim_vector, double dfreedom );
    ~StudentTDistribution();
    void clean();

    inline int dimension() { return normal_dist->dimension(); }
    inline void setLocation( CVector * mean_v ) { normal_dist->setMean( mean_v ); }
    inline void setScaleMatrix( CMatrix * cov ) { normal_dist->setCovariance( cov ); }
    inline CMatrix scaleMatrix() { return normal_dist->covariance(); }
    inline void setDegreesOfFreedom( double dfreedom ) 
                { inv_chisq_dist->setDegreesOfFreedom( dfreedom ); }
    inline double degreesOfFreedom() { return inv_chisq_dist->degreesOfFreedom(); }
    inline InvChisqDistribution * tau2() { return inv_chisq_dist; }
    inline NormalDistribution * mu() { return normal_dist; }

    StudentTDistribution( const StudentTDistribution & distr );
    Distribution * clone();
    StudentTDistribution & operator =(const StudentTDistribution & distr );

    double logDensity( DistributionParameter & value );
    CVector drawOneItem();
    DistributionParameter draw();
    void update( DistributionParameter * par_list, int list_size );
    void update( CVector & mean_shift, CMatrix & cov_shift );
    CVector lastItemDrawn();
    DistributionParameter lastDraw();
    void setLastDraw( DistributionParameter & par );
    void setLastDraw( CVector & par );
    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance();   
};//end class


class WishartDistribution : public Distribution
{
  private:
    int dim;
    double initial_dfreedom;
    double scale_determinant;

    CMatrix * initial_w;
    CMatrix * inv_initial_w;
    CVector * cholesky_initial_w;

    double dfreedom;
    CVector * cholesky_updated_inv_w;
    CMatrix * inv_w;
    CMatrix * last_item_drawn;
  public:
    WishartDistribution();
    WishartDistribution( double df );
    WishartDistribution( double df, CMatrix * v );
    ~WishartDistribution();
    void clean();
    void cleanVectorsAndMatrices();

    inline int dimension() { return dim; }
    inline double degreesOfFreedom() { return dfreedom; }
    CMatrix scaleMatrix();
    inline double initialDegreesOfFreedom() { return initial_dfreedom; }
    inline CMatrix initialScaleMatrix() { return (*initial_w); }
    inline CMatrix inverseInitialScaleMatrix() { return (*inv_initial_w); }
    inline CMatrix inverseScaleMatrix() { return (*inv_w); }

    void setDegreesOfFreedom( double df );
    void setScaleMatrix( CMatrix * v );
    double scaleDeterminant();

    WishartDistribution( const WishartDistribution & distr );
    Distribution * clone();
    WishartDistribution & operator =(const WishartDistribution & distr );

    double logDensity( DistributionParameter & value );
    CMatrix drawOneItem();
    DistributionParameter draw();
    void update( DistributionParameter * par_list, int list_size );
    void update( double add_df, CMatrix & scale_shift );
    CMatrix lastItemDrawn();
    DistributionParameter lastDraw();
    void setLastDraw( DistributionParameter & par );
    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance(); 
}; //end class



class InvWishartDistribution : public Distribution
{
  private:
    int dim;
    double initial_dfreedom;
    double scale_determinant;

    CMatrix * initial_w;
    CMatrix * inv_initial_w;
    CVector * cholesky_initial_w;

    double dfreedom;
    CVector * cholesky_updated_inv_w;
    CMatrix * inv_w;
    CMatrix * last_item_drawn;
    CMatrix * inv_lastDraw;
  public:
    InvWishartDistribution();
    InvWishartDistribution( double df );
    InvWishartDistribution( double df, CMatrix * v );
    ~InvWishartDistribution();
    void clean();
    void cleanVectorsAndMatrices();

    inline int dimension() { return dim; }
    inline double degreesOfFreedom() { return dfreedom; }
    CMatrix scaleMatrix();
    inline double initialDegreesOfFreedom() { return initial_dfreedom; }
    inline CMatrix initialScaleMatrix() { return (*initial_w); }
    inline CMatrix inverseInitialScaleMatrix() { return (*inv_initial_w); }
    inline CMatrix inverseScaleMatrix() { return (*inv_w); }

    void setDegreesOfFreedom( double df );
    void setScaleMatrix( CMatrix * v );
    double scaleDeterminant();

    InvWishartDistribution( const InvWishartDistribution & distr );
    Distribution * clone();
    InvWishartDistribution & operator =(const InvWishartDistribution & distr );

    double logDensity( DistributionParameter & value );
    CMatrix drawOneItem() throw( rtErr );
    DistributionParameter draw();
    void update( DistributionParameter * par_list, int list_size );
    void update( double add_df, CMatrix & scale_shift );
    CMatrix lastItemDrawn();
    CMatrix inverseLastItemDrawn();
    DistributionParameter lastDraw();
    void setLastDraw( DistributionParameter & par );
    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance(); 
}; //end class




/*
  Constrained Inv Chi-2 based on CDF and inverse CDF on Chi-2 
*/

class ConstrainedInvChisqDistribution: public Distribution 
{
  private:
    double last_item_drawn;
    double initial_dfreedom;
    double initial_scale;
    double dfreedom;
    double scale;
    bool left_bounded_only;
    bool right_bounded_only;
    CVector * interval;
  public:
    ConstrainedInvChisqDistribution();
    ConstrainedInvChisqDistribution( double lower, double upper );
    ConstrainedInvChisqDistribution( double bound, char * type );
    virtual ~ConstrainedInvChisqDistribution();
    void clean();

    void setInterval( double lower, double upper );
    void setInterval( double bound, char * type );
    void setValidLastItemDrawn( double sigma02 );
    double upperBound();
    double lowerBound();
    inline double degreesOfFreedom() { return dfreedom; }
    inline double scaleParameter() { return scale; }
    inline double initialDegreesOfFreedom() { return initial_dfreedom; }
    inline double initialScale() { return initial_scale; }

    void setDegreesOfFreedom( double nu0 );
    void setScale( double sigma02 ); //initial value of scale = sigma0 * sigma0

    double CDFInvChisq( double val );
    double inverseCDFInvChisq( double val ) throw( rtErr );

    void update( double additional_df, double additional_scale );
    void update( DistributionParameter * par_list, int list_size );

    ConstrainedInvChisqDistribution & operator=(const ConstrainedInvChisqDistribution & distr );
    ConstrainedInvChisqDistribution( const ConstrainedInvChisqDistribution & distr );
    Distribution * clone();

    double drawOneItem();
    DistributionParameter draw();

    inline double lastItemDrawn() { return last_item_drawn; }
    DistributionParameter lastDraw();

    void setLastDraw( double last_draw );
    void setLastDraw( DistributionParameter & last_draw );

    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance();
    double logDensity( DistributionParameter & value );

};//end class


//Du Mouchel Distribution
class DuMouchelDistribution : public Distribution
{
  private:
    double tau0;
    double tau02;
    double last_item_drawn;
  public:
    DuMouchelDistribution();
    DuMouchelDistribution( double tau );
    ~DuMouchelDistribution();
    void setParameter( double tau );
    inline double paramater() { return tau0; }

    DuMouchelDistribution & operator=(const DuMouchelDistribution & distr );
    DuMouchelDistribution( const DuMouchelDistribution & distr );
    Distribution * clone();

    void update( DistributionParameter * par_list, int list_size );
    double drawOneItem();
    DistributionParameter draw();

    inline double lastItemDrawn() { return last_item_drawn; }
    DistributionParameter lastDraw();

    void setLastDraw( double last_draw );
    void setLastDraw( DistributionParameter & last_draw );

    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance();
    double logDensity( DistributionParameter & value );
}; //end du mouchel class


//Uniform Shrinkage Distribution
class UniformShrinkageDistribution : public Distribution
{
  private:
    double tau02;
    double last_item_drawn;
  public:
    UniformShrinkageDistribution();
    UniformShrinkageDistribution( double tau2 );
    ~UniformShrinkageDistribution();
    void setParameter( double tau2 );
    inline double parameter() { return tau02; }
    
    UniformShrinkageDistribution & operator=(const UniformShrinkageDistribution & distr );
    UniformShrinkageDistribution( const UniformShrinkageDistribution & distr );
    Distribution * clone();

    void update( DistributionParameter * par_list, int list_size );
    double drawOneItem();
    DistributionParameter draw();

    inline double lastItemDrawn() { return last_item_drawn; }
    DistributionParameter lastDraw();

    void setLastDraw( double last_draw );
    void setLastDraw( DistributionParameter & last_draw );

    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance();
    double logDensity( DistributionParameter & value );
}; //end Uniform Shrinkage class



/* Proper non Informative posteriors for Hierarchical linear model */

class ProperNonInfoPosteriorHLM: public Distribution
{
  private:
    double tau02;
    double tau0;
    double dfreedom;
    double scale;
    double last_item_drawn;
    double mode_x;
    double mean_x;
    double var_x;
    double delta_step;
    char * distr_type;
    CVector * prob;
  public:
    ProperNonInfoPosteriorHLM();
    ProperNonInfoPosteriorHLM( double tau2, double df ) throw( rtErr );
    ~ProperNonInfoPosteriorHLM();
    void clean();
    void computeCDF();
    void initialize( double tau2, double df ) throw( rtErr );
    void setScale( double scl );
    void setPriorDistribution( char * type ) throw( rtErr );
    
    ProperNonInfoPosteriorHLM & operator=( const ProperNonInfoPosteriorHLM & distr ); 
    ProperNonInfoPosteriorHLM( const ProperNonInfoPosteriorHLM & distr );
    Distribution * clone();

    double inverseCDF( double p ) throw( rtErr );
    
    DistributionParameter draw();
    void update( DistributionParameter * par_list, int list_size );
    DistributionParameter lastDraw();

    inline double lastItemDrawn() { return last_item_drawn; }

    void setLastDraw( double par );
    void setLastDraw( DistributionParameter & par );
    DistributionParameter mode();
    DistributionParameter mean();
    DistributionParameter variance();
    double logDensity( DistributionParameter & value );
    double logDensity( double val );

};//end class ProperNonInfoPosteriorHLM



#endif /* BAYES_DISTRIBUTION_EXTENDED_H */
