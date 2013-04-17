// Vector.h: interface for the CVector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__C13E7228_138A_4B59_B880_39B200A3D203__INCLUDED_)
#define AFX_VECTOR_H__C13E7228_138A_4B59_B880_39B200A3D203__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CVector
{
public:
    CVector();
    CVector(int iLen);
    CVector(const CVector & A);
    CVector(double *vp, int iLen);
    virtual ~CVector();

    //
    // class utility functions
    //
    int Len() { return m_iLen;}
    int Len() const { return m_iLen;}
    
    double & Val(int i);
    double  Val(int i) const;

    CVector & operator =(const CVector & rhs);

    bool includes( double value );
    void setToZero();
    void setTo( double value );
    CVector subVector( CVector & index );
    CVector subVector( int from_index, int to_index );
    CVector subVectorComplement( CVector & index );
    void setSubVector( CVector & index, CVector & values );
    void add( CVector & a_vector );
    void multiplyByScalar( double weight );
    void setToWeighted( CVector & weights );
    CVector weighted( CVector & weights );
    double lowerTriangularVal( int i, int j );
    void setLowerTriangularVal( int i, int j, double ltval );
    CVector asUpperTriangularMultiplyRight( CVector & b );
    CVector asUpperTriangularSolve( CVector & b );
    CVector asLowerTriangularSolve( CVector & b );
    CVector asCholeskyDecompositionSolve( CVector & b );
    double asCholeskyDecompositionDeterminant();
    void viewAsCholeskyDecomposition();
    void    Print();
    bool isZero( double eps );

    //
    // Vector Functions
    //
    double  Mean();
    double  GMean();
    double  HMean();
    double  Var();
    double Cov( CVector & b );
    CVector square();

    double sum();
    double max();
    int howManyLarger( double v );
    int howManyEqual( double v );
    CVector setEqual( double v );

    //
    // friendly operators
    //
    friend CVector      operator -(const CVector & A,   const double dVal);
    friend CVector      operator -(const CVector & A,   const CVector & B);
    friend CVector      operator +(const CVector & A,   const double dVal);
    friend CVector      operator +(const CVector &A,    const CVector &B);
    friend CVector      operator *(const CVector & A,   const double dVal);
    friend double       operator *(const CVector & A,   const CVector & B);

private:
    int         m_iLen;
    double *    m_pVal;
};

#endif // !defined(AFX_VECTOR_H__C13E7228_138A_4B59_B880_39B200A3D203__INCLUDED_)

