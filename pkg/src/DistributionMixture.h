
#ifndef BAYES_DISTRIBUTION_MIXTURE_H
#define BAYES_DISTRIBUTION_MIXTURE_H

#include "Distribution.h"

class DistributionMixture: public Distribution
{
  private:
    int n_mixtures;
    int last_component_drawn;
    CVector * initial_props;
    CVector * props;
    CVector * cumu_props;
    CVector * draw_stats;
    bool save_stats;
    Distribution ** distr_list;
  public:
    DistributionMixture();
    DistributionMixture( int n_mixs );
    virtual ~DistributionMixture();
    inline bool statisticsAreSaved() { return save_stats; }
    void saveStatistics( bool val );
    void printDrawingStats();
    void printDrawingStats( double * stats );
    void set( int index, Distribution * distr, double prop );
    DistributionParameter * stratifiedDrawFromMixture( int number_draws );
    inline int numberOfMixtures() { return n_mixtures; }
    inline int lastComponentDrawn() { return last_component_drawn; }
    inline double initialProportion( int i ) { return initial_props->Val(i); }
    inline double proportion( int i ) { return props->Val(i); }
    void updateProportion( int i, double val );
    void normalizeProportions();
    void updateLogProportion( int i, double val );
    void normalizeLogProportions();
    void viewCurrentProportions();
    void viewInitialProportions();

    //required methods
    DistributionMixture( const DistributionMixture & distr );
    Distribution * clone();
    void setLastDraw( DistributionParameter & par );
    void setLastDrawForComponent( int comp_index, DistributionParameter & par );
    DistributionParameter lastDraw();
    DistributionParameter lastDrawFromComponent( int coml_index );
    DistributionParameter draw();
    DistributionParameter drawFromComponent( int comp_index );
    DistributionParameter mean();
    DistributionParameter mode();
    DistributionParameter variance();

    void update( DistributionParameter * par_list, int list_size );
    double logDensity( DistributionParameter & value );

};//end class


#endif /* BAYES_DISTRIBUTION_MIXTURE_H */
