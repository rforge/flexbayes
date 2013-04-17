// Matrix.h: interface for the CMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIX_H__531212B9_1ED4_4AC8_937D_F15F79C8C97B__INCLUDED_)
#define AFX_MATRIX_H__531212B9_1ED4_4AC8_937D_F15F79C8C97B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Vector.h"
#include <stdio.h>
#include "rtErr.h"

class CMatrix
{
public:
    CMatrix();
    CMatrix(int iRow, int iCol);
    CMatrix(const CMatrix & A);
    CMatrix(double *pVal, int iRow, int iCol);
    CMatrix(double *pVal, int iRow, int iCol, bool by_row );
    CMatrix(CVector & pVal, int iRow, int iCol);
    virtual ~CMatrix();

    //
    // functions on Matrix
    //
    CMatrix     T();
    CMatrix     Cholesky();
    CMatrix     Inverse_LT();
    CMatrix     Inverse() throw( rtErr );
    double      Det();  
    //
    // utility functions
    //
    CMatrix & operator =(const CMatrix & rhs);

    int Col() { return m_iCol;}
    int Row() { return m_iRow;}
    int Col() const { return m_iCol;}
    int Row() const { return m_iRow;}

    double & Val(int iRow, int iCol);
    double Val(int iRow, int iCol) const;
	    
    CVector getColumn( int iCol ) throw( rtErr );
    CVector getRow( int iRow ) throw( rtErr );
    void setColumn( int iCol, CVector& col_vector ) throw( rtErr );
    void setSubColumn( int iCol, CVector & index, CVector& col_vector ) throw( rtErr );
    void setRow( int iRow, CVector& row_vector ) throw( rtErr );
    void setSubRow( int iRow, int from_index, int to_index, CVector & row_vector ) throw( rtErr );
    CMatrix getRows( int init_row, int end_row ) throw( rtErr );
    CMatrix getRows( CVector & index ) throw( rtErr );
    CMatrix subMatrix( CVector & index ) throw( rtErr );
    CMatrix subMatrixComplement( CVector & index ) throw( rtErr );
    CMatrix subMatrixCross( CVector & index ) throw( rtErr );
    void setToZero();
    void setDiagonal( double value ) throw( rtErr );
    void setToWeightedColumn( int iCol, CVector& weights ) throw( rtErr );
    void weighted( CVector& weights );
    void multiplyByScalar( double weight );
    void add( CMatrix & a_matrix ) throw( rtErr );
    CVector choleskyDecomposition() throw( rtErr );
    bool exploreCholeskyDecomposition(  CVector * cholesky_dec );
    void assignInverseOfLowerTriangular( CVector & lTriangular ) throw( rtErr );
    void assignLowerTriangular( CVector & lTriangular ) throw( rtErr );
    CMatrix inverse() throw( rtErr );
    CMatrix xTransposedX();
    CMatrix xTransposedX( CVector & weights ) throw( rtErr );
    CMatrix xTransposedX( CMatrix & weights ) throw( rtErr );
    void xxT( CVector & weights ) throw( rtErr );
    double trace() throw( rtErr );
    bool isZero( double eps );
    void toArray( double * pvals, bool by_row );
    void toArray( double * pvals, bool by_row, int start_index );
    CVector toVector();

    void Print();
    // void Init(double x);

    //
    // friend operators
    //
    friend CMatrix      operator +(const CMatrix & A, const CMatrix & B);
    friend CMatrix      operator +(const CMatrix & A, const CVector & B);
    friend CMatrix      operator -(const CMatrix & A, const CMatrix & B);
    friend CMatrix      operator -(const CMatrix & A, const CVector & B);
    friend CMatrix      operator *(const double dA,   const CMatrix & B);
    friend CMatrix      operator *(const CMatrix & A, const CMatrix & B);
    friend CVector      operator *(const CMatrix & A,   const CVector & B);
    friend CMatrix      TT(const CVector & A, const CVector & B);
    friend double       M_D(const CMatrix &A, const CVector &B);
protected:
    int         m_iRow;
    int         m_iCol;
    double  *   m_pVal;
};

#endif // !defined(AFX_MATRIX_H__531212B9_1ED4_4AC8_937D_F15F79C8C97B__INCLUDED_)

