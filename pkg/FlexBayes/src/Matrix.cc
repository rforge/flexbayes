// Matrix.cpp: implementation of the CMatrix class.
//
//////////////////////////////////////////////////////////////////////
#include "R.h"

#include "Matrix.h"
#include "Const.h"
#include "Algebra.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMatrix::CMatrix()
:   m_iRow(0),
    m_iCol(0),
    m_pVal(NULL)
{
}


CMatrix::CMatrix(int iRow, int iCol)
:   m_iRow(iRow),
    m_iCol(iCol)
{
    int i, j, size = iRow * iCol;
    m_pVal = new double[ size ];

    // reset the data area
    // Change by Dawn: the following was:
    //memset((void *) m_pVal, 0, size * sizeof(double));
    for( i=0; i<iRow; i++ )
      for( j=0; j<iCol; j++ )
        Val( i, j ) = 0.0;
}


CMatrix::CMatrix(const CMatrix & A)
:   m_iRow(A.Row()),
    m_iCol(A.Col())
{
    int size = m_iRow * m_iCol;
    m_pVal = new double [ size ];

    for (int i = 0; i < size; i ++)
    {
        m_pVal[i] = A.m_pVal[i];
    }
}


CMatrix::CMatrix(double *pVal, int iRow, int iCol)
{
    int size = iRow * iCol;
    m_pVal = new double[size];

    memcpy(m_pVal, pVal, size * sizeof(double));
    m_iRow = iRow;
    m_iCol = iCol;
}


CMatrix::CMatrix( double *pVal, int iRow, int iCol, bool by_row )
{
  int i, j, k;
    int size = iRow * iCol;
    m_pVal = new double [ size ];

    m_iRow = iRow;
    m_iCol = iCol;

    if ( by_row )
    {
      memcpy( m_pVal, pVal, size * sizeof( double ) );
    }
    else
    {
      k = 0;
      for ( j = 0; j < m_iCol; j++ )
      {
        for ( i = 0; i < m_iRow; i++ )
        {
          Val(i, j) = pVal[ k ];
          k++;
        }
      }
    }//end by column
}


CMatrix::CMatrix(CVector & pVal, int iRow, int iCol)
{
  int i, j, k;
  int size = iRow * iCol;
  m_pVal = new double[ size ];
  m_iRow = iRow;
  m_iCol = iCol;

  //assume by column
  k = 0;
  for ( j = 0; j < m_iCol; j++ )
  {
    for ( i = 0; i < m_iRow; i++ )
    {
      Val(i, j) = pVal.Val( k );
      k++;
    }
  }

}

double & CMatrix::Val(int iRow, int iCol)
{
  #ifdef BOUNDS_CHECK
    if( (iRow < 0) || (iRow >= m_iRow) || (iCol < 0) || (iCol >= m_iCol) )
    MESSAGE "Index out of bounds in matrix access" ERROR;
  #endif 
  return m_pVal[iRow * m_iCol + iCol];
}

double CMatrix::Val(int iRow, int iCol) const
{
  #ifdef BOUNDS_CHECK
    if( (iRow < 0) || (iRow >= m_iRow) || (iCol < 0) || (iCol >= m_iCol) )
    MESSAGE "Index out of bounds in matrix access" ERROR;
  #endif 
  return m_pVal[iRow * m_iCol + iCol];
}


CVector CMatrix::getRow( int iRow )  throw( rtErr )
{
  int j;

  if ( iRow < Row() && iRow >= 0 )
  {
    CVector row_vector( m_iCol );

    for ( j = 0; j < m_iCol; j++ )
    {
      row_vector.Val( j ) = Val( iRow, j );
    }

    return row_vector;
  }
  else
  {
    Rprintf( "CMatrix::getRow: Row index [%d] is out of range. Max row is %d.\n\n", iRow, Row() );

    char the_error[] = "CMatrix::getRow: Row index is out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}


CVector CMatrix::getColumn( int iCol )  throw( rtErr )
{
  int i;

  if ( iCol < Col() && iCol >= 0 )
  {
    CVector col_vector( m_iRow );

    for ( i = 0; i < m_iRow; i++ )
    {
      col_vector.Val( i ) = Val( i, iCol );
    }

    return col_vector;
  }
  else
  {
    Rprintf( "CMatrix::getColumn: Column index [%d] is out of range. Max column is %d.\n\n", iCol, Col() );

    char the_error[] = "CMatrix::getColumn: Column index is out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}

CMatrix CMatrix::subMatrix( CVector & index )  throw( rtErr )
{
  int i, j, i_index, j_index;
  CMatrix sub_m( index.Len(), index.Len() );

  if ( Col() == Row() && index.Len() <= Col() )
  {
    for ( i = 0; i < index.Len(); i++ )
    {
      for ( j = 0; j < index.Len(); j++ )
      {
        i_index = (int) index.Val( i );
        j_index = (int) index.Val( j );
        sub_m.Val( i, j ) = Val( i_index, j_index );
      }
    }
    return sub_m;
  }
  else if ( Col() != Row() )
  {
    Rprintf( "CMatrix::subMatrix: need a square matrix.\n" );
    char the_error[] = "CMatrix::subMatrix: need a square matrix.";
    rtErr runtime_error( the_error );
    throw runtime_error;    
  }
  else if ( index.Len() > Col() )
  {
    Rprintf( "CMatrix::subMatrix: vector of indexes is too long.\n" );
    char the_error[] = "CMatrix::subMatrix: vector of indexes is too long.";
    rtErr runtime_error( the_error );
    throw runtime_error; 
  }
  else
  {
    Rprintf( "Reached past end of control in CMatrix CMatrix::subMatrix( CVector & index ).\n" );
    return sub_m;
  }

}//end



CMatrix CMatrix::subMatrixComplement( CVector & index )  throw( rtErr )
{
  int i;

  if ( Col() == Row() && index.Len() < Col() )
  {
    CVector all_indexes( Row() );
    for ( i = 0; i < Row(); i++ )
    {
      all_indexes.Val(i) = i;
    }

#ifdef FIX1
    CVector tmpvec = all_indexes.subVectorComplement( index );
    CMatrix sub_m( subMatrix( tmpvec ) );
#else
    CMatrix sub_m( subMatrix( all_indexes.subVectorComplement( index ) ) );
#endif
//    CMatrix sub_m( subMatrix( all_indexes.subVectorComplement( index ) ) );
 
    return sub_m;    
  }
  else if ( Col() != Row() )
  {
    Rprintf( "CMatrix::subMatrix: need a square matrix.\n" );
    char the_error[] = "CMatrix::subMatrix: need a square matrix.";
    rtErr runtime_error( the_error );
    throw runtime_error; 
  }
  else if ( index.Len() >= Col() )
  {
    Rprintf( "CMatrix::subMatrix: vector of indexes is too long.\n" );
    char the_error[] = "CMatrix::subMatrix: vector of indexes is too long.";
    rtErr runtime_error( the_error );
    throw runtime_error; 
  }
  else
  {
    Rprintf( "Reached past end of control in CMatrix CMatrix::subMatrixComplement( CVector & index ).\n" );
  }

  CMatrix sub_m( 0, 0 );
  return sub_m;
}//end


CMatrix CMatrix::subMatrixCross( CVector & index )  throw( rtErr )
{
  int i, j, i_index, j_index;

  if ( Col() == Row() && index.Len() < Col() )
  {
    CVector all_indexes( Row() );
    for ( i = 0; i < Row(); i++ )
    {
      all_indexes.Val(i) = i;
    }

    CVector compl_index( all_indexes.subVectorComplement( index ) );
    CMatrix sub_m( index.Len(), compl_index.Len() );

    for ( i = 0; i < index.Len(); i++ )
    {
      for ( j = 0; j < compl_index.Len(); j++ )
      {
        i_index = (int) index.Val( i );
        j_index = (int) compl_index.Val( j );
        sub_m.Val( i, j ) = Val( i_index, j_index );
      }
    }

    return sub_m;
  }
  else if ( Col() != Row() )
  {
    Rprintf( "CMatrix::subMatrix: need a square matrix.\n" );
    char the_error[] = "CMatrix::subMatrix: need a square matrix.\n";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
  else if ( index.Len() >= Col() )
  {
    Rprintf( "CMatrix::subMatrix: vector of indexes is too long.\n" );
    char the_error[] = "CMatrix::subMatrix: vector of indexes is too long.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
  else
  {
    Rprintf( "Reached past end of control in CMatrix CMatrix::subMatrixCross( CVector & index ).\n" );
  }

  CMatrix sub_m( 0, 0 );
  return sub_m;
}//end






void CMatrix::setRow( int iRow, CVector& row_vector )  throw( rtErr )
{
  int j;

  if ( iRow < Row() && iRow >= 0 )
  {
    for ( j = 0; j < m_iCol; j++ )
    {
      Val( iRow, j ) = row_vector.Val( j );
    }
  }
  else
  {
    Rprintf( "CMatrix::setRow: Row index [%d] is out of range. Max row is %d.\n\n", iRow, Row() );

    char the_error[] = "CMatrix::setRow: Row index is out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}


void CMatrix::setColumn( int iCol, CVector& col_vector ) throw( rtErr )
{
  int i;

  if ( iCol < Col() && iCol >= 0 )
  {
    for ( i = 0; i < m_iRow; i++ )
    {
      Val( i, iCol ) = col_vector.Val( i );
    }
  }
  else
  {
    Rprintf( "CMatrix::setColumn: Column index [%d] is out of range. Max column is %d.\n\n", iCol, Col() );

    char the_error[] = "CMatrix::setColumn: Column index is out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  } 
}



void CMatrix::setSubColumn( int iCol, CVector & index, CVector& col_vector )  throw( rtErr )
{
  int i, idx;

  if ( index.Len() <= Row() )
  {
    for ( i = 0; i < index.Len(); i++ )
    {
      idx = (int) index.Val( i );
      Val( idx, iCol ) = col_vector.Val( i );
    }
  }
  else
  {
    Rprintf( "CMatrix::setSubColumn: wrong index dimensions.\n" );
    char the_error[] = "CMatrix::setSubColumn: wrong index dimensions.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  } 
}//end 


void CMatrix::setSubRow( int iRow, int from_index, int to_index, CVector & row_vector )  throw( rtErr )
{
  int i;

  if ( to_index > Col() || from_index < 0 || to_index < from_index )
  {
    Rprintf( "CMatrix::setSubRow: wrong index dimensions.\n" );
    char the_error[] = "CMatrix::setSubRow: wrong index dimensions.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
  else
  {
    for ( i = from_index; i < to_index; i++ )
    {
      Val( iRow, i ) = row_vector.Val( i - from_index );
    }
  }
}//end 



CMatrix CMatrix::getRows( int init_row, int end_row )  throw( rtErr )
{

  if ( init_row >= 0 && end_row < Row() )
  {
    int i, total_rows;

    total_rows = end_row - init_row + 1;
    CMatrix sub_matrix( total_rows, Col() );

    for ( i = 0; i < total_rows; i++ )
    {
#ifdef FIX1
      CVector tmpvec = getRow( init_row + i );
      sub_matrix.setRow( i , tmpvec );
#else
      sub_matrix.setRow( i, getRow( init_row + i ) );
#endif
//      sub_matrix.setRow( i, getRow( init_row + i ) );
    }
 
    return sub_matrix;
  }
  else
  {
    Rprintf( "CMatrix::getRows: indeces out of range.\n" );
    //return CMatrix();
    char the_error[] = "CMatrix::getRows: indeces out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}



CMatrix CMatrix::getRows( CVector & index )  throw( rtErr )
{
  int i, i_index, j;

  if ( index.Len() <= Row() )
  {
    CMatrix sub_m( index.Len(), Col() );
    for ( i = 0; i < index.Len(); i++ )
    {
#ifdef FIX1
      i_index = static_cast< int >( index.Val( i ) );
#else
      i_index = index.Val( i );
#endif
//      i_index = index.Val( i );
      for ( j = 0; j < Col(); j++ )
      {
        sub_m.Val( i, j ) = Val( i_index, j );
      }
    }

    return sub_m;
  }
  else
  {
    Rprintf( "CMatrix::getRows: indeces out of range.\n" );
    //return CMatrix();
    char the_error[] = "CMatrix::getRows: indeces out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
}//end





void CMatrix::setToZero()
{
  int i, j;

  for ( i = 0; i < m_iRow; i++ )
  {
    for ( j = 0; j < m_iCol; j++ )
    {
      Val( i, j ) = 0.0;
    }
  }

}


void CMatrix::setDiagonal( double value )  throw( rtErr )
{
  int i;

  if ( m_iRow != m_iCol )
  {
    Rprintf( "CMatrix:setDiagonal: matrix need to be a squared matrix.\n" );
    char the_error[] = "CMatrix:setDiagonal: matrix need to be a squared matrix.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }
  else
  {
    for ( i = 0; i < m_iRow; i++ )
    {
      Val( i, i ) = value;
    }
  }

}//end


void CMatrix::weighted( CVector& weights )
{
  int j;

  for ( j = 0; j < Col(); j++ )
  {
    setToWeightedColumn( j, weights );
  }
}


void CMatrix::setToWeightedColumn( int iCol, CVector& weights )  throw( rtErr )
{
  int i;

  if ( iCol < Col() )
  {
    for ( i = 0; i < Row(); i++ )
    {
      Val(i, iCol) *= weights.Val(i);
    }
  }
  else
  {
    Rprintf( "CMatrix::setToWeightedColumn: Column index [%d] is out of range. Max column is %d.\n\n", iCol, Col() );
    char the_error[] = "CMatrix::setToWeightedColumn: Column index is out of range.";
    rtErr runtime_error( the_error );
    throw runtime_error;
  }   
}


void CMatrix::multiplyByScalar( double weight )
{
  int i, j;

  for ( i = 0; i < Row(); i++ )
  {
    for ( j = 0; j < Col(); j++ )
    {
      Val(i, j) *= weight;
    }
  }

}



void CMatrix::add( CMatrix & a )  throw( rtErr )
{
  int i, j;

  if ( a.Row() == Row() && a.Col() == Col() )
  {
    for ( i = 0; i < Row(); i++ )
    {
      for ( j = 0; j < Col(); j++ )
      {
        Val(i, j) += a.Val(i, j);
      }
    }
  }
  else
  {
    Rprintf( "CMatrix::add: Matrix to add has wrong dimensions [%d, %d]. Current dimensions are [%d, %d]\n", a.Row(), a.Col(), Row(), Col() );
    char the_error[] = "CMatrix::add: Matrix to add has wrong dimensions.";
    rtErr runtime_error( the_error );
    throw runtime_error;    
  }

}


CMatrix::~CMatrix()
{
    delete [] m_pVal;
    m_pVal = NULL;
}



CMatrix CMatrix::xTransposedX()
{
  
  int i, j, k;

  CMatrix cov( Col(), Col() );

  for ( i = 0; i < Col(); i++ )
  {
    for ( j = 0; j < Col(); j++ )
    {
      cov.Val(i, j) = 0;
      for ( k = 0; k < Row(); k++ )
      {
        cov.Val(i, j) += Val(k, i) * Val(k, j);
      }
    }
  }

  return cov;
}


CMatrix CMatrix::xTransposedX( CMatrix & weights )  throw( rtErr )
{
  CMatrix cov( Col(), Col() );

  if ( weights.Row() != weights.Col() || weights.Col() != Row() )
  {
    Rprintf( "Matrix::xTransposedX( CMatrix weights ): matrix of weights has wrong dimensions.\n" );
    char the_error[] = "Matrix::xTransposedX( CMatrix weights ): matrix of weights has wrong dimensions.";
    rtErr runtime_error( the_error );
    throw runtime_error;     
  }
  else
  {
    cov = ( T() * weights ) * (*this);

    return cov;
  }
}


CMatrix CMatrix::xTransposedX( CVector & weights )  throw( rtErr )
{
  int i, j, k;

  CMatrix cov( Col(), Col() );

  if ( weights.Len() == Row() )
  {
    for ( i = 0; i < Col(); i++ )
    {
      for ( j = 0; j < Col(); j++ )
      {
        cov.Val(i, j) = 0;
        for ( k = 0; k < Row(); k ++ )
        {
          cov.Val(i, j) += Val(k, i) * Val(k, j) * weights.Val(k);
        }
      }
    }

    return cov;
  }
  else
  {
    Rprintf( "CMatrix::xTransposedX( Cvector weights ): weights vector has wrong dimension [%d]. Row dimension is %d\n\n", weights.Len(), Row() );
    char the_error[] = "CMatrix::xTransposedX( Cvector weights ): weights vector has wrong dimension.";
    rtErr runtime_error( the_error );
    throw runtime_error;  
  }
}



void CMatrix::xxT( CVector & weights )  throw( rtErr )
{
  int i, j;

  if ( weights.Len() == Col() )
  {
    for ( i = 0; i < Col(); i++ )
    {
      for ( j = 0; j < Col(); j++ )
      {
        Val(i, j) += weights.Val( i ) * weights.Val( j );
      }
    }
  }
  else
  {
    Rprintf( "CMatrix::xxT( Cvector weights ): weights vector has wrong dimension [%d]. Column dimension is %d\n\n", weights.Len(), Col() );
    char the_error[] = "CMatrix::xxT( Cvector weights ): weights vector has wrong dimension.";
    rtErr runtime_error( the_error );
    throw runtime_error;  
  }
}//end 



double CMatrix::trace()  throw( rtErr )
{
  int i;
  double sum;

  sum = 0;
  if ( Col() == Row() )
  {
    for ( i = 0; i < Col(); i++ )
    {
      sum += Val( i, i );
    }

    return sum;
  }
  else
  {
    Rprintf( "CMatrix::trace: Invalid call to trace. Matrix is not a square matrix.\n" );
    char the_error[] = "CMatrix::trace: Invalid call to trace. Matrix is not a square matrix.";
    rtErr runtime_error( the_error );
    throw runtime_error; 
  }
}//end


bool CMatrix::isZero( double eps )
{
  int i, j;
  bool found;

  found = false;
  i = 0;
  j = 0;
  while ( !found && i < Row() && j < Col() )
  {
    if ( fabs( Val(i,j) ) > eps )
    {
      found = true;
    }
    else
    {
      if ( j < Col() - 1 )
      {
        j++;
      }
      else
      {
        j = 0;
        i++;
      }
    }
  }//end loop

  return (!found);

}//end



void CMatrix::toArray( double * pvals, bool by_row )
{
  int i, j, k, size;

  if ( by_row )
  {
    size = m_iRow * m_iCol;
    memcpy( pvals, m_pVal, size * sizeof( double ) );
  }
  else //by column
  {
    k = 0;
    for ( j = 0; j < m_iCol; j++ )
    {
      for ( i = 0; i < m_iRow; i++ )
      {
        pvals[ k ] = Val(i, j);
        k++;
      }
    }
  }//end by column   

}//end


void CMatrix::toArray( double * pvals, bool by_row, int start_index )
{
  int i, j, k, size;

  if ( by_row )
  {
    size = m_iRow * m_iCol;
    memcpy( &(pvals[start_index]), m_pVal, size * sizeof( double ) );
  }
  else //by column
  {
    k = 0;
    for ( j = 0; j < m_iCol; j++ )
    {
      for ( i = 0; i < m_iRow; i++ )
      {
        pvals[ start_index + k ] = Val(i, j);
        k++;
      }
    }
  }//end by column   

}//end



CVector CMatrix::toVector()
{
  int i, j, k;

  CVector pvals( m_iCol * m_iRow );
  k = 0;
  for ( j = 0; j < m_iCol; j++ )
  {
    for ( i = 0; i < m_iRow; i++ )
    {
      pvals.Val( k ) = Val(i, j);
      k++;
    }
  }

  return pvals;
}//end



void CMatrix::Print()
{
    Rprintf("iRow = %d, iCol = %d\n", m_iRow, m_iCol);

    if ( m_iRow <= 0 || m_iCol <= 0 )
    {
      Rprintf( "wrong dimensions is matrix.\n" );
    }
    else
    {    
    for (int i = 0; i < m_iRow; i ++)
    {
        for (int j = 0; j < m_iCol; j ++)
        {
            if (m_pVal[i * m_iCol + j] >= 0)
            {
                Rprintf(" ");
            }
            Rprintf(" %f  ", m_pVal[i * m_iCol + j]);
            
        }
        Rprintf("\n");
    }
    Rprintf("\n");
    }
}

/*
void CMatrix::Init(double x)
{
    for (int i = 0; i < m_iRow; i ++)
    {
        for (int j = 0; j < m_iCol; j ++)
        {
            m_pVal[i * m_iCol + j] = x;
            x += 1.0;
        }
    }
}
*/


CMatrix & CMatrix::operator =(const CMatrix & rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    delete [] m_pVal;

    m_iRow = rhs.Row();
    m_iCol = rhs.Col();

    int size = m_iRow * m_iCol;
    m_pVal = new double[ size];

    for (int i = 0; i < size; i ++)
    {
        m_pVal[i] = rhs.m_pVal[i];
    }

    return *this;
}


CMatrix  CMatrix::T()
{
    CMatrix  B(m_iCol, m_iRow);

    for (int i = 0; i < m_iCol; i ++)
    {
        for (int j = 0; j < m_iRow; j ++)
        {
            B.Val(i,j) = Val(j,i);
        }
    }

    return B;
}


////////////////////////////////////////////////////////////////////////////////
//
//        Lower Cholesky
//
//  Abstract:
//
//  Arguments:
//
//      A=LL', A: positive definite; L: lower triangular.
//      A -- (n x n) matrix.
//      n -- the dimensions of A.
//
//  Return Value:
//      0 - successful. A -- L.

CMatrix CMatrix::Cholesky()
{
    if(m_iRow != m_iCol)
    {
      Rprintf( "***CMatrix::Cholesky: needs square matrix. ***");
     
      return *this;
    }

    int i, j, k;
    
    CMatrix B(*this);
    double  s = 0.0;
    for (k = 0; k < m_iRow; k ++)
    {
        s = 0.0;

        for (j = 0; j < k; j ++)
        {
            s += B.Val(k,j)* B.Val(k,j);
        }

        B.Val(k, k) = sqrt(B.Val(k, k) - s);

        if (B.Val(k, k) < DELTA)
        {
	    Rprintf("***CMatrix::Cholesky: nonpositive matrix! ***");
           
            return B;
        }

        for (i = k + 1; i < m_iRow; i ++)
        {
            s = 0.0;
            for ( j = 0; j < k; j ++)
            {
                s += B.Val(i, j) * B.Val(k, j);
            }

            B.Val(i, k) = (B.Val(i, k) - s) / B.Val(k, k);
        }
    }

    for (i = 0; i < m_iRow - 1; i ++)
    {
        for (j = i + 1; j < m_iRow; j ++)
        {
            B.Val(i, j) = 0.0;
        }
    }
 
    return B;
}


////////////////////////////////////////////////////////////////////////////////
//
//        Inverse the lower triangular matrix
//
//  Abstract:
//
//  Arguments:
//
//      B -- pxp lower triangluar matrix.
//
//  Return Value:
//      C -- B^{-1}.
//
CMatrix CMatrix:: Inverse_LT()
{
   
  int i,j;
  if(m_iRow != m_iCol)
    {
      Rprintf("*** Inverse_LT: Need squared matrix. ***");
    }

    CMatrix C(*this);

    for (i = 0; i < m_iCol ; i ++)
    {
        if ( fabs( Val(i,i) ) <= DELTA )
        {
	        Rprintf("*** Inverset_LT: Error: Diagonal element is too small. ***");

            return C;
        }
    }

    double x = 0.0;
    for (j = 0; j < m_iCol; j ++)
    {
        C.Val(j, j) = 1.0 / Val(j, j);

        for (i = j + 1; i < m_iCol ; i ++)
        {
            x = 0.0;

            for (int k = j; k < i; k++)
            {
                x += Val(i, k) * C.Val(k, j);
            }

            C.Val(i,j) = -x / Val(i, i);
         }
    }

    return C;
} // Inverse_LT




CVector CMatrix::choleskyDecomposition() throw( rtErr )
{
  int i, j, starting_index;
  double * working_matrix;
  double * working_diagonal;
  CVector cholesky_dec( ( m_iRow * ( m_iRow + 1 ) ) / 2 );

  working_matrix = new double[ m_iRow * m_iRow ];
  working_diagonal = new double[ m_iRow ];
  for ( i = 0; i < m_iRow; i++ )
  {
    for ( j = 0; j < m_iRow; j++ )
    {
      working_matrix[ i * m_iRow + j ] = Val(i, j);
    }
  }

  try
  {
    Algebra::choldc( working_matrix, m_iRow, working_diagonal, 0.0 );
    //Algebra::choldc( working_matrix, (long) m_iRow, working_diagonal, 0.0 );
  }
  catch( rtErr choldcError )
  {
    Rprintf( "CMatrix::choleskyDecomposition: Matrix does not seem to be symmetric positive definite.\n" );
    char the_error[] = "CMatrix::choleskyDecomposition: Matrix does not seem to be symmetric positive definite.";
    //rtErr runtime_error( the_error );
    //throw runtime_error;
    MESSAGE "" ERROR;

    //PROBLEM choldcError.what WARN;
    //PROBLEM "\nCMatrix::choleskyDecomposition: Matrix is not symmetric positive definite." WARN;
    //throw;
  }
  
  for ( i = 0; i < m_iRow; i++ )
  {
    starting_index = ( i * ( i - 1 ) ) / 2;
    cholesky_dec.setLowerTriangularVal( i, i, working_diagonal[ i ] );
    for ( j = i - 1; j >= 0; j-- )
    {
      cholesky_dec.setLowerTriangularVal( i, j, working_matrix[ i * m_iRow + j ] );
    }
  }

  delete [] working_matrix;
  delete [] working_diagonal;

  return cholesky_dec;

}


bool CMatrix::exploreCholeskyDecomposition(  CVector * cholesky_dec )
{
  bool is_chol_OK;
  int i, j, starting_index;
  double * working_matrix;
  double * working_diagonal;

  //CVector cholesky_dec( ( m_iRow * ( m_iRow + 1 ) ) / 2 );

  working_matrix = new double[ m_iRow * m_iRow ];
  working_diagonal = new double[ m_iRow ];
  for ( i = 0; i < m_iRow; i++ )
  {
    for ( j = 0; j < m_iRow; j++ )
    {
      working_matrix[ i * m_iRow + j ] = Val(i, j);
    }
  }

  is_chol_OK = Algebra::isCholdcOK( working_matrix, m_iRow, working_diagonal, 0.0 );
  //is_chol_OK = Algebra::isCholdcOK( working_matrix, (long) m_iRow, working_diagonal, 0.0 );

  if ( is_chol_OK )
  {  
    for ( i = 0; i < m_iRow; i++ )
    {
      starting_index = ( i * ( i - 1 ) ) / 2;
      cholesky_dec->setLowerTriangularVal( i, i, working_diagonal[ i ] );
      for ( j = i - 1; j >= 0; j-- )
      {
        cholesky_dec->setLowerTriangularVal( i, j, working_matrix[ i * m_iRow + j ] );
      }
    }
  }

  delete [] working_matrix;
  delete [] working_diagonal;

  return ( is_chol_OK );

}//end



void CMatrix::assignLowerTriangular( CVector & lTriangular ) throw( rtErr )
{
  int i, j, lTLen;
  double root_lTLen;

  if( m_iRow != m_iCol )
  {
    Rprintf( "CMatrix::assignLowerTriangular: Invalid call to assignLowerTriangular. Matrix is not a square matrix.\n" );
    char the_error[] = "CMatrix::assignLowerTriangular: Invalid call to assignLowerTriangular. Matrix is not a square matrix.";
    rtErr runtime_error( the_error );
    throw runtime_error; 
  }

  root_lTLen = ( sqrt( 8 * ((double) lTriangular.Len() ) + 1 )  - 1 ) / 2;
  lTLen = (int) root_lTLen;

  if ( lTLen == m_iRow && m_iRow == m_iCol )
  {
    for ( i = 0; i < lTLen; i++ )
    {
      for ( j = 0; j < lTLen; j++ )
      {
        Val( i, j ) = lTriangular.lowerTriangularVal( i, j );
      }
    }
  }
  else
  {
    Rprintf( "CMatrix::assignLowerTriangular: wrong dimensions.\n" );
    char the_error[] = "CMatrix::assignLowerTriangular: wrong dimensions.";
    rtErr runtime_error( the_error );
    throw runtime_error; 
  }

}//end




void CMatrix::assignInverseOfLowerTriangular( CVector & lTriangular ) throw( rtErr )
{
  int i, j, k, lTLen;
  double root_lTLen, sum;

  if( m_iRow != m_iCol )
  {
    Rprintf( "CMatrix::assignInverseOfLowerTriangular: Need squared matrix.\n" );
    char the_error[] = "CMatrix::assignInverseOfLowerTriangular: Need squared matrix.";
    rtErr runtime_error( the_error );
    throw runtime_error; 
  }

  root_lTLen = ( sqrt( 8 * ((double) lTriangular.Len() ) + 1 )  - 1 ) / 2;
  lTLen = (int) root_lTLen;

  if ( lTLen == m_iRow && m_iRow == m_iCol )
  {
    for ( i = 0; i < m_iCol ; i ++ )
    {
        if ( fabs( lTriangular.lowerTriangularVal(i, i) ) <= DELTA )
        {
	  Rprintf( "CMatrix::assignInverseOfLowerTriangular: Error: Diagonal element [%f] is too small.\n", lTriangular.lowerTriangularVal(i, i) );
          char the_error[] = "CMatrix::assignInverseOfLowerTriangular: Error: Diagonal element is too small.";
          rtErr runtime_error( the_error );
          throw runtime_error; 
        }
    }

    sum = 0.0;
    for ( j = 0; j < m_iCol; j++ )
    {
        Val(j, j) = 1.0 / lTriangular.lowerTriangularVal(j, j);

        for ( i = j + 1; i < m_iCol ; i++ )
        {
            sum = 0.0;
            for ( k = j; k < i; k++ )
            {
                sum += lTriangular.lowerTriangularVal(i, k) * Val(k, j);
            }
            Val(i, j) = -sum / lTriangular.lowerTriangularVal(i, i);
	}//end for i
    }//end for j
  }//end if
  else
  {
    Rprintf( "CMatrix::assignInverseOfLowerTriangular: wrong dimensions.\n" );
    char the_error[] = "CMatrix::assignInverseOfLowerTriangular: wrong dimensions.";
    rtErr runtime_error( the_error );
    throw runtime_error; 
  }

  /* another version 
  root_lTLen = ( sqrt( 8 * ((double) lTriangular.Len() ) + 1 )  - 1 ) / 2;
  lTLen = (int) root_lTLen;

  if ( lTLen == m_iRow && m_iRow == m_iCol )
  {
    CVector b( lTLen );

    for ( ibase = 0; ibase < lTLen; ibase++ )
    {
      //get ibase-th canonical base vector
      for ( i = 0; i < lTLen; i++ )
      {
        b.Val( i ) = 0;
      }
      b.Val( ibase ) = 1;

      CVector wv( lTriangular.asLowerTriangularSolve( b ) );
      setColumn( ibase, wv );
    }
  }
  else
  {
    Rprintf( "CMatrix::assignInverseOfLowerTriangular: wrong dimensions.\n" );
  }
  */

}
  

CMatrix CMatrix::inverse() throw( rtErr )
{
  CMatrix inverse_matrix( m_iRow, m_iRow );

  if ( m_iRow > 1 )
  {
    CVector cholLT( ( m_iRow * ( m_iRow + 1) ) / 2 );
    CMatrix inverse_cholLT( m_iRow, m_iRow );

    cholLT = choleskyDecomposition();
    inverse_cholLT.assignInverseOfLowerTriangular( cholLT );

    inverse_matrix = inverse_cholLT.T() * inverse_cholLT;
  }
  else
  {
    if ( m_pVal[0] >= DELTA || m_pVal[0] < - DELTA )
    {
      inverse_matrix.Val( 0, 0 ) = 1.0 / m_pVal[0];
    }
    else
    {
      Rprintf( "CMatrix::inverse: Matrix element [%f] is too small.  Matrix:\n", m_pVal[0] );
      this->Print();
      char the_error[] = "CMatrix::inverse:Matrix element is too small.\n";
      rtErr runtime_error( the_error );
      throw runtime_error;       
    }
  }


  return inverse_matrix;
}


CMatrix CMatrix::Inverse() throw( rtErr )
{
    CMatrix C(*this);

    if ( m_iRow > 1 )
    {
      CMatrix D(m_iRow, m_iCol);
    
      D = (C.Cholesky()).Inverse_LT();
      C = D.T() * D;
    }
    else
    {
      if ( m_pVal[0] >= DELTA || m_pVal[0] < - DELTA )
      {
        C.Val( 0, 0 ) = 1.0 / m_pVal[0];
      }
      else
      {
        Rprintf( "CMatrix::Inverse: Matrix element [%f] is too small.\n", m_pVal[0] );
        char the_error[] = "CMatrix::inverse:Matrix element is too small.\n";
        rtErr runtime_error( the_error );
        throw runtime_error;       
      }
    }
   
    return C;
} // Inverset

double  CMatrix::Det()
{
    CMatrix D(m_iRow, m_iCol);    
    D = Cholesky();
    double a = 1.0;
    for (int i = 0; i < m_iRow; i ++)
      {
	a *= D.Val(i,i);
      }
    a = a * a;
    return a;
} // Determination


////////////////////////////////////////////////////////////////////////////////
//
// Friend operators
//
CMatrix operator +(const CMatrix & A, const CMatrix & B)
{
    CMatrix     C(A.Row(), B.Col());

    for (int i = 0; i < A.Row(); i ++)
    {
        for (int j = 0; j < B.Col(); j ++)
        {
            C.Val(i, j) = A.Val(i, j) + B.Val(i, j);
        }
    }

    return C;
}


CMatrix operator +(const CMatrix & A, const CVector & B)
{
    CMatrix     C(A.Row(), A.Col());

    if ( A.Col() != B.Len() )
    {
      Rprintf( "CMatrix: operator+: Dimensions of Matrix and Vector do not match.\n" );
      // Bug fix by Dawn: the following was:  exit(1);    which causes problems within S-PLUS
      MESSAGE "CMatrix: operator+: Dimensions of Matrix and Vector do not match.\n" ERROR;
    }

    for (int i = 0; i < A.Row(); i ++)
    {
        for (int j = 0; j < B.Len(); j ++)
        {
            C.Val(i, j) = A.Val(i, j) + B.Val(j);
        }
    }

    return C;
}


CMatrix operator -(const CMatrix & A, const CMatrix & B)
{
    CMatrix     C(A.Row(), B.Col());

    for (int i = 0; i < A.Row(); i ++)
    {
        for (int j = 0; j < B.Col(); j ++)
        {
            C.Val(i, j) = A.Val(i, j) - B.Val(i, j);
        }
    }

    return C;
}



CMatrix operator -(const CMatrix & A, const CVector & B)
{
    CMatrix     C(A.Row(), A.Col());

    if ( A.Col() != B.Len() )
    {
      Rprintf( "CMatrix: operator+: Dimensions of Matrix and Vector do not match.\n" );
      // Bug fix by Dawn: the following was:   exit(1);   which causes problems within S-PLUS
      MESSAGE "CMatrix: operator+: Dimensions of Matrix and Vector do not match.\n" ERROR;
    }

    for (int i = 0; i < A.Row(); i ++)
    {
        for (int j = 0; j < B.Len(); j ++)
        {
            C.Val(i, j) = A.Val(i, j) - B.Val(j);
        }
    }

    return C;
}

CMatrix operator *(const double dA, const CMatrix & B)
{
    CMatrix     C(B);

    int iRow = B.Row();
    int iCol = B.Col();

    for (int i = 0; i < iRow; i ++)
    {
        for (int j = 0; j < iCol; j ++)
        {
            C.Val(i, j) = dA * B.Val(i, j);
        }
    }

    return C;
}


CMatrix operator *(const CMatrix & A, const CMatrix & B)
{
    int aCol = A.Col();
    int bCol = B.Col();
    int aRow = A.Row();

    CMatrix     C(aRow, bCol);

    for (int i = 0; i < aRow; i ++)
    {
        for (int j = 0; j < bCol; j ++)
        {
            C.Val(i, j) = 0;
            // B.Col() == n
            for (int k = 0; k < aCol; k ++)
            {
                C.Val(i, j) += A.Val(i, k) * B.Val(k, j);
            }
        }
    }

    return C;
}

CVector operator *(const CMatrix & A, const CVector & B)
{
    int Row = A.Row();
    int Col = B.Len();

    CVector    C(Row);
    if (A.Col() != Col){
    	Rprintf( "operator *: Attempt to multiplying matrix with [%d] columns by a vector of length [%d]\n", A.Col(), Col );
      MESSAGE "operator *: Invalid matrix multiplication" ERROR;
    }
    for (int i = 0; i < Row; i ++)
    {
        // should be inited C.Val(i) = 0.0;
        if(0.0 != C.Val(i)){
        	MESSAGE "operator *: Incorrect matrix initialization" ERROR;
        }
        	

        for (int j = 0; j < Col; j ++)
        {
	        C.Val(i) += A.Val(i, j) * B.Val(j);
        }
    }

    return C;
}

CMatrix TT(const CVector & A, const CVector & B)
{
    int LenA = A.Len();
    int LenB = B.Len();

    CMatrix    C(LenA, LenB);

    for (int i = 0; i < LenA; i ++)
    {
        for (int j = 0; j < LenB; j ++)
        {
	        C.Val(i,j) = A.Val(i) * B.Val(j);
        }
    }

    return C;
}

double M_D(const CMatrix  & A, const CVector & B)
{
  int LenB = B.Len();
  CVector    C(LenB);
  
  C = A * B;
  
  return ( B * C);
}
