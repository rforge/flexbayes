
#ifndef BAYES_DISTRIBUTION_PARAMETER_H
#define BAYES_DISTRIBUTION_PARAMETER_H

#include "Matrix.h"

class DistributionParameter
{
  protected:
    CMatrix * pars;
  public:
    DistributionParameter();
    DistributionParameter( int iRow, int iCol );
    DistributionParameter( double scalar );
    DistributionParameter( CVector & vector_pars );
    DistributionParameter( CMatrix & matrix_pars );
    DistributionParameter( const DistributionParameter & dist_par );
    virtual ~DistributionParameter();
    DistributionParameter & operator =(const DistributionParameter & rhs);
    void setParameter( double scalar );
    void setParameter( CVector & vector_pars );
    void setParameter( CMatrix & matrix_pars );
    void updateParameter( double scalar );
    void updateParameter( CVector & vector_pars );
    void updateParameter( CMatrix & matrix_pars );
    bool isScalar();
    bool isVector();
    bool isMatrix();
    inline double getScalar() { return pars->Val( 0, 0 ); }
    inline CVector getVector() { return pars->getColumn( 0 ); }
    inline CMatrix getMatrix() { return (*pars); }
    inline CMatrix parameter() { return (*pars); }
    void Print();
    void toArray( double * pvals, int start_index );
    int length();
};//end

#endif /* BAYES_DISTRIBUTION_PARAMETER_H */
