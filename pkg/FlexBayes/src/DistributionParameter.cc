#include "DistributionParameter.h"


DistributionParameter::DistributionParameter()
{
  pars = NULL;
}//end


DistributionParameter::DistributionParameter( int iRow, int iCol )
{
  pars = new CMatrix( iRow, iCol );
}//end


DistributionParameter::DistributionParameter( double scalar )
{
  pars = NULL ;
  setParameter( scalar );
}//end


DistributionParameter::DistributionParameter( CVector & vector_pars )
{
  pars = NULL ;
  setParameter( vector_pars );
}//end


DistributionParameter::DistributionParameter( CMatrix & matrix_pars )
{
  pars = NULL ;
  setParameter( matrix_pars );
}//end


DistributionParameter::DistributionParameter( const DistributionParameter & dist_par )
{
  pars = new CMatrix( (*dist_par.pars) );
}//end


DistributionParameter::~DistributionParameter()
{
  if ( pars != NULL )
  {
    delete pars;
    pars = NULL;
  }
}//end


void DistributionParameter::setParameter( double scalar )
{
  if ( pars != NULL) {
    delete pars ;
    pars = NULL ;
  }
  pars = new CMatrix( 1, 1);
  pars->Val(0,0) = scalar;
}//end


void DistributionParameter::setParameter( CVector & vector_pars )
{
  if ( pars != NULL) {
    delete pars ;
  }
  pars = new CMatrix( vector_pars.Len(), 1 );
  pars->setColumn( 0, vector_pars );
}//end


void DistributionParameter::setParameter( CMatrix & matrix_pars )
{
  if ( pars != NULL) {
    delete pars ;
  }
  pars = new CMatrix( matrix_pars );
}//end



void DistributionParameter::updateParameter( double scalar )
{
  if ( pars == NULL )
  {
    pars = new CMatrix( 1, 1);
  }

  pars->Val(0,0) = scalar;
}//end


void DistributionParameter::updateParameter( CVector & vector_pars )
{
  if ( pars != NULL )
  {
    delete pars;
  }

  pars = new CMatrix( vector_pars.Len(), 1 );
  pars->setColumn( 0, vector_pars );
}//end


void DistributionParameter::updateParameter( CMatrix & matrix_pars )
{
  if ( pars != NULL )
  {
    delete pars;
  }
   
  pars = new CMatrix( matrix_pars );
}//end



bool DistributionParameter::isScalar()
{
  bool is_scalar;
  is_scalar =  ( pars->Row() == 1 && pars->Col() == 1 )? true : false;

  return is_scalar;
}//end


bool DistributionParameter::isVector()
{
  bool is_vector;
  is_vector =  ( pars->Row() > 1 && pars->Col() == 1 )? true : false;

  return is_vector;
}//end


bool DistributionParameter::isMatrix()
{
  bool is_matrix;
  is_matrix =  ( pars->Row() > 1 && pars->Col() > 1 )? true : false;

  return is_matrix;
}//end


void DistributionParameter::Print()
{
  //printf( "\nDistributionParameter::Print\n"); fflush(stdout);
  pars->Print();
}//end


DistributionParameter & DistributionParameter::operator =(const DistributionParameter & rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    delete pars;

    int iCol = rhs.pars->Col();
    int iRow = rhs.pars->Row();
    pars = new CMatrix( iRow, iCol );
    (*pars) = (*(rhs.pars));

    return *this;
}//end


void DistributionParameter::toArray( double * pvals, int start_index )
{
  bool by_row;
  
  by_row = false;
  pars->toArray( pvals, by_row, start_index );    
}//end


int DistributionParameter::length()
{

  return pars->Row() * pars->Col();
}//end
