// Vector.cpp: implementation of the CVector class.
//
//////////////////////////////////////////////////////////////////////
#include "R.h"

#include "Vector.h"
#include "Algebra.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CVector::CVector()
:   m_iLen(0),
    m_pVal(NULL)
{

}

CVector::CVector(int iLen)
:   m_iLen(iLen)

{
    m_pVal = new double[ iLen ];
    // Change by Dawn: the following was:
    //memset((void *)m_pVal, 0, m_iLen * sizeof(double));
    for( int i=0; i<iLen; i++ )
      m_pVal[i] = 0.0;
}

CVector::CVector(const CVector & A)
{
  int i;

    m_iLen = A.Len();
    m_pVal = new double [ m_iLen ];

    for (i = 0; i < m_iLen; i++)
    {
        m_pVal[i] = A.Val(i);
    }
}


CVector::CVector(double *vp, int iLen)
  : m_iLen(iLen)
{
  int i;

  m_pVal = new double[ iLen ];
  for (i = 0; i < m_iLen; i ++)
    {
        m_pVal[i] = vp[i];
    }  
}

CVector::~CVector()
{
    //printf( " in Vector delete: m_iLen=%d\n", m_iLen ); fflush(stdout);
  

    delete [] m_pVal;

    //printf(" in Vector: deleted\n" ); fflush(stdout);

    m_pVal = NULL;
}

double & CVector::Val(int i)
{
  #ifdef BOUNDS_CHECK
    if( (i < 0) || (i >= m_iLen) )
      MESSAGE "Index out of bounds in vector access" ERROR;
  #endif 
  return m_pVal[i];
}

double CVector::Val(int i) const
{
  #ifdef BOUNDS_CHECK
    if( (i < 0) || (i >= m_iLen) )
      MESSAGE "Index out of bounds in vector access" ERROR;
  #endif 
  return m_pVal[i];
}

CVector & CVector::operator =(const CVector & rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    delete [] m_pVal;

    m_iLen = rhs.Len();

    m_pVal = new double[ m_iLen ];

    for (int i = 0; i < m_iLen; i ++)
    {
        m_pVal[i] = rhs.m_pVal[i];
    }

    return *this;
}


bool CVector::includes( double value )
{
  int i;
  bool found;

  found = false;
  i = 0;
  while ( !found && i < Len() )
  {
    if ( Val( i ) == value )
    {
      found = true;
    }
    else
    {
      i++;
    }
  }//end loop

  return found;

}//end


bool CVector::isZero( double eps )
{
  int i;
  bool found;

  found = false;
  i = 0;
  while ( !found && i < Len() )
  {
    if ( fabs( Val( i ) ) > eps )
    {
      found = true;
    }
    else
    {
      i++;
    }
  }//end loop

  return (!found);

}//end



void CVector::setToZero()
{
  int i;
  for ( i = 0; i < m_iLen; i++ )
  {
    m_pVal[i] = 0.0;
  }
}


void CVector::setTo( double value )
{
  int i;
  for ( i = 0; i < m_iLen; i++ )
  {
    m_pVal[i] = value;
  }
}


void CVector::add( CVector & a_vector )
{
  int i;

  if ( a_vector.Len() == m_iLen )
  {
    for ( i = 0; i < m_iLen; i++ )
    {
      m_pVal[i] += a_vector.Val( i );
    }
  }
  else
  {
    Rprintf( "CVector::add: Adding vector dimension [%d] differs from this vector dimension [%d]\n", a_vector.Len(), m_iLen );
    MESSAGE "" ERROR;
  }
}


void CVector::multiplyByScalar( double weight )
{
  int i;
  for ( i = 0; i < m_iLen; i++ )
  {
    m_pVal[i] *= weight;
  }
}


void CVector::setToWeighted( CVector & weights )
{
  int i;
 
  for ( i = 0; i < m_iLen; i++ )
  {
    m_pVal[i] *= weights.Val(i);
  }

}


CVector CVector::weighted( CVector & weights )
{
  int i;
 
  CVector w( m_iLen );

  for ( i = 0; i < m_iLen; i++ )
  {
    w.Val(i) = m_pVal[i] * weights.Val(i);
  }

  return w;
}


CVector CVector::subVector( CVector & index )
{
  int i, j;

  if ( index.Len() > Len() )
  {
    Rprintf( "CVector::subVector: vector of indexes is too long.\n" );
    // Bug fix by Dawn: the following was: exit(1);
    MESSAGE "CVector::subVector: vector of indexes is too long.\n" ERROR;
  }

  CVector sub_v( index.Len() );
  for ( i = 0; i < index.Len(); i++ )
  {
    if ( ((int) index.Val(i)) >= Len() )
    {
      Rprintf( "CVector::subVector: vector of indexes is too long.\n" );
      // Bug fix by Dawn: the following was: exit(1);
      MESSAGE "CVector::subVector: vector of indexes is too long.\n" ERROR;
    }

#ifdef FIX1
    j =  static_cast< int >( index.Val( i ) );
#else
    j =  index.Val( i );
#endif
//    j = index.Val( i );
    sub_v.Val( i ) = Val( j );
  }

  return sub_v;
}//end


CVector CVector::subVector( int from_index, int to_index )
{
  int i;

  if ( to_index - from_index + 1 > Len() )
  {
    Rprintf( "CVector::subVector: list of indexes is too long.\n" );
    // Bug fix by Dawn: the following was:  exit(1);  which causes problems within S-PLUS
    MESSAGE "CVector::subVector: list of indexes is too long.\n" ERROR;
  }

  if ( from_index < 0 )
  {
    Rprintf( "CVector::subVector: negative index.\n" );
    // Bug fix by Dawn: the following was:  exit(1);  which causes problems within S-PLUS
    MESSAGE "CVector::subVector: negative index.\n" ERROR;
  }

  if ( to_index - from_index < 0 )
  {
    Rprintf( "CVector::subVector: negative or zero span for indexes.\n" );
    // Bug fix by Dawn: the following was:  exit(1);  which causes problems within S-PLUS
    MESSAGE "CVector::subVector: negative or zero span for indexes.\n" ERROR;
  }

  CVector sub_v( to_index - from_index + 1 );
  for ( i = from_index; i <= to_index; i++ )
  {
    sub_v.Val( i - from_index ) = Val( i );
  }

  return sub_v;
}//end





CVector CVector::subVectorComplement( CVector & index )
{
  bool found;
  int i, j, k;

  if ( index.Len() >= Len() )
  {
    Rprintf( "CVector::subVector: vector of indexes is too long.\n" );
    // Bug fix by Dawn: the following was: exit(1); which causes problems within SPLUS
    MESSAGE "CVector::subVector: vector of indexes is too long.\n" ERROR;
  }

  CVector sub_v( Len() - index.Len() );
  k = 0;
  for ( i = 0; i < Len(); i++ )
  {
    found = false;
    j = 0;
    while ( !found && j < index.Len() )
    {
      if ( i == ( (int) index.Val( j ) ) )
      {
        found = true;
      }
      else
      {
        j++;
      }
    }//end loop

    if ( !found )
    {
      sub_v.Val( k ) = Val( i );
      k++;
    }
  }//end for i

  return sub_v;
}//end



void CVector::setSubVector( CVector & index, CVector & values )
{
  int i, j;

  if ( index.Len() > Len() )
  {
    Rprintf( "CVector::setSubVector: vector of indexes is too long.\n" );
    // Bug fix by Dawn: the following was:   exit(1);  which causes problems within S-PLUS
    MESSAGE "CVector::setSubVector: vector of indexes is too long.\n" ERROR;
  }

  if ( index.Len() > values.Len() )
  {
    Rprintf( "CVector::setSubVector: vector of indexes is longer than vector of values.\n" );
    // Bug fix by Dawn: the following was:   exit(1);  which causes problems within S-PLUS
    MESSAGE "CVector::setSubVector: vector of indexes is longer than vector of values.\n" ERROR;
  }

  for ( i = 0; i < index.Len(); i++ )
  {
    if ( ((int) index.Val(i)) >= Len() )
    {
      Rprintf( "CVector::setSubVector: vector of indexes is too long.\n" );
      // Bug fix by Dawn: the following was:   exit(1);  which causes problems within S-PLUS
      MESSAGE "CVector::setSubVector: vector of indexes is too long.\n" ERROR;
    }

    if ( index.Val(i) < 0 )
    {
      Rprintf( "CVector::setSubVector: negative index.\n" );
      // Bug fix by Dawn: the following was:   exit(1);  which causes problems within S-PLUS
      MESSAGE "CVector::setSubVector: negative index.\n" ERROR;
    }

#ifdef FIX1
    j =  static_cast< int >( index.Val( i ) );
#else
    j =  index.Val( i );
#endif
//    j =  index.Val( i );
    Val( j ) = values.Val( i );
  }

}//end




/* returns the corresponding (i,j) element of Lower Triangular matrix 
   Storage:
   row i starts at element (i*(i+1))/2 and its followed by the values
   of all columns from 0 until i;
*/
double CVector::lowerTriangularVal( int i, int j )
{
  int index;
  double ltval;

  if ( i >= j )
  {
    index = ( i * (i + 1 ) )/2 + j;
    ltval = Val( index );
  }
  else
  {
    ltval = 0;
  }

  return ltval;
}

void CVector::setLowerTriangularVal( int i, int j, double ltval )
{
  int index;

  if ( i >= j )
  {
    index = ( i * (i + 1 ) )/2 + j;
    Val( index ) = ltval;
  }
  else
  {
    Rprintf( "Indexes (%d, %d) are not valid for lower Triangular matrix\n", i, j );
  }
}


/* assumes vector is  a Triangular matrix L stored as Lower Triangular 
   it returns x = L^t * b
*/
CVector CVector::asUpperTriangularMultiplyRight( CVector & b )
{
  int i, j;

  CVector x( b.Len() );

  for ( i = 0; i < b.Len(); i++ )
  {
    x.Val(i) = 0;
    for ( j = i; j < b.Len(); j++ )
    {
      x.Val(i) += lowerTriangularVal(j, i) * b.Val(j);
    }
  }

  return x;
}      


CVector CVector::asUpperTriangularSolve( CVector & b )
{
  int i, j;
  double sum;
  CVector x( b.Len() );
 
  for ( i = b.Len() - 1; i >= 0; i-- )
  {
    sum = b.Val( i );
    for ( j = i + 1; j < b.Len(); j++ )
    {
      sum -= lowerTriangularVal( j, i ) * x.Val( j );
    }
    x.Val( i ) = sum / lowerTriangularVal( i, i );
  }

  return x;
}


CVector CVector::asLowerTriangularSolve( CVector & b )
{
  int i, j;
  double sum;
  CVector x( b.Len() );
 
  if ( ( b.Len() * (b.Len() + 1) )/2 != Len() )
  {
    Rprintf( "CVector::asLowerTriangularSolve: lengths differ: b.Len = %d, Len = %d\n", b.Len(), Len() );
  }
  else
  {
    for ( i = 0; i < b.Len(); i++ )
    {
      sum = b.Val( i );
      for ( j = i - 1; j >= 0; j-- )
      {
        sum -= lowerTriangularVal( i, j ) * x.Val( j );
      }
      x.Val( i ) = sum / lowerTriangularVal( i, i );
    }
  }

  return x;
}


CVector CVector::asCholeskyDecompositionSolve( CVector & b )
{
  int i, j;
  double * working_matrix;
  double * working_diagonal;
  double * vb;
  double * vx;

  working_matrix = new double [ b.Len() * b.Len() ];
  working_diagonal = new double [ b.Len() ];
  vb = new double [ b.Len() ];
  vx = new double [ b.Len() ];

  for ( i = 0; i < b.Len(); i++ )
  {
    vb[ i ] = b.Val(i);
    working_diagonal[ i ] = lowerTriangularVal( i, i );
    for ( j = i - 1; j >= 0; j-- )
    {
      working_matrix[ i * b.Len() + j ] = lowerTriangularVal( i, j );
      //upper triangle is not needed
      //working_matrix[ j * b.Len() + i ] = working_matrix[ i * b.Len() + j ];
    }
  }

  //solve A vx = vb
  //Algebra::cholsl( working_matrix, (long) b.Len(), working_diagonal, vb, vx );
  Algebra::cholsl( working_matrix, b.Len(), working_diagonal, vb, vx );
  CVector x( vx, b.Len() );

  delete [] working_matrix;
  delete [] working_diagonal;
  delete [] vb;
  delete [] vx;

  return x;  
}


double CVector::asCholeskyDecompositionDeterminant()
{
  int rows, i;
  double det, root_row;

  root_row = ( sqrt( 8 * ((double) m_iLen) + 1 )  - 1 ) / 2;
  rows = (int) root_row;

  det = 1.0;
  for ( i = 0; i < rows; i++ )
  {
    det *= lowerTriangularVal( i, i );
  }
  det = det * det;

  return det;
}


void CVector::viewAsCholeskyDecomposition()
{
  int rows, i, j;
  double root_row;

  root_row = ( sqrt( 8 * ((double) m_iLen) + 1 )  - 1 ) / 2;
  rows = (int) root_row;

  Rprintf( "CVector:viewAsCholeskyDecomposition: root = %f, rows = %d\n", root_row, rows );

  for ( i = 0; i < rows; i++ )
  {
    for ( j = 0; j < rows; j++ )
    {
      Rprintf(  "%f ", lowerTriangularVal( i, j ) );
    }
    Rprintf( "\n" );
  }

}


void CVector::Print()
{
  int i;

    Rprintf("Len = %d\n", m_iLen);

    for ( i = 0; i < m_iLen; i ++)
    {
        if (m_pVal[i] >= 0.0)
        {
            Rprintf(" ");
        }
        Rprintf("%f  ", m_pVal[i]);
    }
    Rprintf("\n\n");
}


double CVector::max()
{
  int i;
  double max_val;

  max_val = m_pVal[ 0 ];
  if ( m_iLen > 1 )
  {
    for ( i = 1; i < m_iLen; i++ )
    {
      if ( m_pVal[ i ] > max_val )
      {  
        max_val = m_pVal[ i ];
      }
    }
  }

  return max_val;
}//end


int CVector::howManyLarger( double a )
{
  int i, count;

  count = 0;
  for ( i = 0; i < m_iLen; i++ )
  {
    if ( m_pVal[ i ] > a )
    {
      count++;
    }
  }

  return count;    
}//end


int CVector::howManyEqual( double a )
{
  int i, count;

  count = 0;
  for ( i = 0; i < m_iLen; i++ )
  {
    if ( m_pVal[ i ] == a )
    {
      count++;
    }
  }

  return count;    
}//end



CVector CVector::setEqual( double a )
{
  int i, count;

  CVector equalSet( m_iLen );

  count = 0;
  for ( i = 0; i < m_iLen; i++ )
  {
    if ( m_pVal[ i ] == a )
    {
      equalSet.Val( count ) = i;
      count++;
    }
  }

  if ( count == 0 )
  {
    equalSet.Val(0) = -1;
    count = 1;
  }

  return equalSet.subVector( 0, count - 1 );    

}//end



double CVector::sum()
{
    double a = 0.0;
    for (int i = 0; i < m_iLen; i ++)
    {
        a += m_pVal[i];
    }

    return a;
}




double CVector::Mean()
{
    double a = 0.0;
    for (int i = 0; i < m_iLen; i ++)
    {
        a += m_pVal[i];
    }

    a = a / m_iLen;
    return a;
}


double CVector::GMean()
{
    double a = 0.0;
    for (int i = 0; i < m_iLen; i ++)
    {
        a += log(m_pVal[i]);
    }

    a = exp( a / m_iLen);
    return a;
}

double CVector::HMean()
{
    double a = 0.0;
    for (int i = 0; i < m_iLen; i ++)
    {
        a += 1.0 /m_pVal[i];
    }

    a = m_iLen / a;
    return a;
}

double CVector::Var()
{
  /*
    double a = 0.0;

    for (int i = 0; i < m_iLen; i ++)
    {
        a += m_pVal[i] * m_pVal[i];
    }

    double dMean = Mean();

    a -= dMean * dMean * m_iLen;

    */

  int i;
  double vector_mean, var;

  vector_mean = Mean();
  var = 0;
  for ( i = 0; i <  m_iLen; i++ )
  {
    var += ( m_pVal[i] - vector_mean ) * ( m_pVal[i] - vector_mean );
  }

  if ( m_iLen > 1 )
  {
    var /= ( (double) (m_iLen - 1) );
  }

  return var;
}


double CVector::Cov( CVector & b )
{
  int i;
  double mean_a, mean_b, covariance;

  mean_a = Mean();
  mean_b = b.Mean();
  
  covariance = 0;
  for ( i = 0; i < m_iLen; i++ )
  {
    covariance += ( Val(i) - mean_a ) * ( b.Val(i) - mean_b );
  }
  if ( m_iLen > 1 )
  {
    covariance /= ( ( double) m_iLen - 1.0 );
  }

  return covariance;
}


CVector CVector::square()
{
  int i;
  CVector vector_square( m_iLen );

  for ( i = 0; i <  m_iLen; i++ )
  {
    vector_square.Val(i) =  m_pVal[i] * m_pVal[i];
  }

  return vector_square;
}


/////////////////////////////////////////////////////////////////////////
//
// Friendly Operators
//
CVector operator +(const CVector &A, const double  dVal)
{
    int Len = A.Len();

    CVector    C(A);

    for (int i = 0; i < Len; i ++)
    {
        C.Val(i) = dVal + A.Val(i);
    }
    return C;
}

CVector operator +(const CVector & A, const CVector & B)
{
    int Len = A.Len();

    CVector    C(Len);

    for (int i = 0; i < Len; i ++)
    {
        C.Val(i) = A.Val(i) + B.Val(i);
    }
    return C;
}

CVector operator -(const CVector &A, const double dVal)
{
    int Len = A.Len();

    CVector    C(A);

    for (int i = 0; i < Len; i ++)
    {
        C.Val(i) =  A.Val(i) - dVal;
    }
    return C;
}

CVector operator -(const CVector & A, const CVector & B)
{
    int Len = A.Len();

    CVector    C(Len);

    for (int i = 0; i < Len; i ++)
    {
        C.Val(i) = A.Val(i) - B.Val(i);
    }
    return C;
}

CVector operator *(const CVector & A,   const double dVal)
{
    int Len = A.Len();

    CVector    C(Len);

    for (int i = 0; i < Len; i ++)
    {
        C.Val(i) = A.Val(i) *dVal;
    }
    return C;
}

double operator *(const CVector & A, const CVector & B)
{
    int Len = A.Len();

    double a = 0.0;
    for (int i = 0; i < Len; i ++)
    {
        a += A.Val(i) * B.Val(i);
    }
    return a;
}



