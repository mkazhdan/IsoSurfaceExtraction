/* -*- C++ -*-
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/
#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED
#include <cmath>
#include <cassert>
#include <complex>
#include <vector>
#include "Algebra.h"

#define PAN_FIX 1
#define FAST_POINT 1

template< class Real > Real Random( void );

template< class Real , int Dim >
class Point : public InnerProductSpace< Real , Point< Real , Dim > >
{
public:
	/////////////////////////////////
	// Inner product space methods //
	void Add            ( const Point& p );
	void Scale          ( Real s );
	Real InnerProduct   ( const Point< Real , Dim >& p ) const;
	/////////////////////////////////

	Real coords[Dim];
	Point ( void ) { memset( coords , 0 , sizeof(Real)*Dim ); }
	template<class Real2>
	operator Point< Real2, Dim > ( void ) const
	{
		Point< Real2, Dim > p;
		for( int d=0 ; d<Dim ; d++ ) p.coords[d] = Real2( coords[d] ); 
		return p;
	}
	Real& operator [] (int idx) { return coords[idx]; }
	const Real& operator [] (int idx) const { return coords[idx]; }

	////////////////////////////////////////////////////////////////////////////
	/*! Writes an ASCII representation of the point to an output stream.
	//  @param[in]  os   output stream
	//  @param[in]  p    the vector to be output
	//  @return     reference to the output stream for output operator chaining
	*///////////////////////////////////////////////////////////////////////////
	friend std::ostream & operator<<(std::ostream &os, const Point &p)
	{
		os << "[";
		for (int i = 0; i < Dim; ++i)   {
			if (i)
				os << ", ";
			os << p[i];
		}

		os << "]";

		return os;
	}

};

////////////////////////////////////////////////////////////////////////////////
/*! Convenience function to create Point<Real, 2>s
//  @param[in]  x   the first component
//  @param[in]  x   the second component
//  @return     Point<Real, 2> initialized with the two components
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
Point<Real, 2> MakePoint2(Real x, Real y)
{
	Point<Real, 2> result;
	result[0] = x;
	result[1] = y;
	return result;
}

////////////////////////////////////////////////////////////////////////////////
/*! Complex-real vector dot product
//  @param[in]  a   complex vector lhs
//  @param[in]  b   real vector rhs
//  @return     (complex result) a dot b
*///////////////////////////////////////////////////////////////////////////////
template <class Real, int Dim>
inline std::complex<Real> operator*(const Point<std::complex<Real>, Dim> &a
	, const Point<Real              , Dim> &b)
{
	std::complex<Real> dot = 0;
	for (int i = 0; i < Dim; ++i)   {
		dot += a[i] * b[i];
	}
	return dot;
}

template<class Real,int Cols,int Rows>
class Matrix : public InnerProductSpace<Real,Matrix<Real,Cols,Rows> >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add            ( const Matrix& m );
	void Scale          ( Real s );
	Real InnerProduct   ( const Matrix& m ) const;
	//////////////////////////

	Real coords[Cols][Rows];
	Matrix ( void ) { memset( coords , 0 , sizeof( Real ) * Cols * Rows ); }
	template<class Real2>
	operator Matrix< Real2 , Cols , Rows > ( void ) const
	{
		Matrix< Real2, Cols , Rows > m;
		for( int c=0 ; c<Cols ; c++ ) for ( int r=0 ; r<Rows ; r++ ) m.coords[c][r] = Real2( coords[c][r] ); 
		return m;
	}
	template<int C,int R>
	Matrix(const Matrix<Real,C,R>& m)
	{
		for(int i=0;i<Cols && i<C;i++)
			for(int j=0;j<Rows && j<R;j++)
				coords[i][j]=m.coords[i][j];
	}
	Real& operator () (int c,int r) { return coords[c][r]; }
	const Real& operator () (int c,int r) const { return coords[c][r]; }

	template<int Cols1>
	Matrix<Real,Cols1,Rows> operator * ( const Matrix< Real , Cols1 , Cols >& m ) const;

	Matrix<Real,Rows,Cols> transpose( void )    const;

	template<class Real2>
	Point<Real2,Rows> operator * ( const Point< Real2 , Cols >& v ) const;
	template<class Real2>
	Point<Real2,Rows> operator () ( const Point< Real2 , Cols >& v ) const;
};

template< class Real , int Dim >
class SquareMatrix : public Algebra< Real , SquareMatrix< Real , Dim > > , public Matrix< Real , Dim , Dim >
{
public:
	using Matrix< Real , Dim , Dim >::coords;
	////////////////////////////////
	// Additional algebra methods //
	void Multiply (const SquareMatrix& m);
	void SetIdentity(void);
	////////////////////////////////

	SquareMatrix( const Matrix< Real , Dim , Dim >& m ){ memcpy( coords , m.coords , sizeof(Real)*Dim*Dim ); }
	SquareMatrix( void )                               { memset( coords , 0 ,        sizeof(Real)*Dim*Dim ); }
	static SquareMatrix Identity( void ){ SquareMatrix m ; for( int i=0 ; i<Dim ; i++ ) m(i,i) = (Real)1 ; return m; }
	Real subDeterminant(int c,int r) const;
	Real determinant( void ) const;
	Real trace( void ) const;
	SquareMatrix inverse( bool& success ) const;
	SquareMatrix inverse( void ) const;

	using Matrix<Real,Dim,Dim>::operator *;
	using Matrix<Real,Dim,Dim>::operator ();
	template<class Real2>
	Point<Real2,Dim-1> operator () (const Point<Real2,Dim-1>& v) const;
};
template< class V , int Dim , class _R = typename V::R >
class Gradient : public VectorSpace< _R , Gradient< V , Dim , _R > >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add            ( const Gradient& g ) { for( int c=0  ; c<Dim ; c++ ) gradients[c] += g.gradients[c]; }
	void Scale          ( _R s ) { for( int c=0 ; c<Dim ; c++ ) gradients[c] *= s; }
	//                      //
	//////////////////////////

	V gradients[Dim];
	Gradient( void ) { for( int d=0 ; d<Dim ;  d++ ) gradients[d] *= 0; }
	V& operator[] ( int idx ) { return gradients[idx]; }
	const V& operator[] ( int idx ) const { return gradients[idx]; }

	template< class V2 , class _R2>
	operator Gradient< V2, Dim , _R2 > ( void ) const
	{
		Gradient< V2 , Dim , _R2 > g;
		for( int d=0 ; d<Dim ; d++ ) g.gradients[d] = V2( gradients[d] ); 
		return g;
	}

	template< class Real >
	Gradient Project( const Point< Real , Dim >& dir ) const
	{
		V dot;
		Gradient g;
		g *= 0;
		dot *= 0;
		Real len = Real( sqrt( Point< Real , Dim >::SquareNorm( dir ) ) );
		if( !len ) return g;
		Point< Real , Dim > _dir = dir / len;
		for( int d=0 ; d<Dim ; d++ ) dot += gradients[d] * _dir[d];
		for( int d=0 ; d<Dim ; d++ ) g.gradients[d] = dot * _dir[d];
		return g;
	}
};

template< class V , int Dim , class _R = typename V::R >
class ConstantFunction : public VectorSpace< _R , ConstantFunction< V , Dim , _R > >
{
public:
	V value;
	Gradient< V , Dim , _R > gradients;
	ConstantFunction( void ) { value *= 0 , gradients *= 0;}

	template< class Real > V operator( ) ( const Point< Real , Dim >& p ) const { return value; }
	template< class Real > Gradient< V , Dim , _R > gradient( const Point< Real , Dim >& p ) const { return gradients; }

	//////////////////////////
	// Vector space methods //
	void Add            ( const ConstantFunction& cf ) { value += cf.value; }
	void Scale          ( _R s ) { value *= s , this->offset *= s; }
	//////////////////////////
};

template< class V , int Dim , class _R = typename V::R >
class LinearFunction : public VectorSpace< _R , LinearFunction< V , Dim , _R > >
{
public:
	Gradient< V , Dim , _R > gradients;
	V offset;
	LinearFunction( void ) { offset *= 0 ; }
	template< class Real >
	V operator( ) ( const Point< Real , Dim >& p ) const
	{
		V v;
		v *= 0;
		for( int d=0 ; d<Dim ; d++ ) v += gradients[d] * p[d];
		v -= offset;
		return v;
	}
	template< class Real >
	LinearFunction fitToHyperplane( const Point< Real , Dim >& p , const Point< Real , Dim >& n ) const
	{
		LinearFunction f;
		Real len = Point< Real , Dim >::SquareNorm( n );
		if( !len )
		{
			f.gradients *= 0;
			f.offset = -(*this)( p );
		}
		else
		{
			Point< Real , Dim > normal = n / Real( sqrt( double( len ) ) );
			V dot;
			dot *= 0;
			for( int d=0 ; d<Dim ; d++ ) dot += gradients[d] * normal[d];
			for( int d=0 ; d<Dim ; d++ ) f.gradients[d] = gradients[d] - dot * normal[d];
			f.offset *= 0;
			f.offset = -(*this)( p ) + f( p );
		}
		return f;
	}
	template< class V2 , class _R2 >
	operator LinearFunction< V2 , Dim , _R2 > ( void ) const
	{
		LinearFunction< V2 , Dim , _R2 > lf;
		lf.offset = V2 ( offset );
		lf.gradients = Gradient< V2 , Dim , _R2 >( gradients );
		return lf;
	}
	template< class Real >
	Gradient< V , Dim , _R > gradient( const Point< Real , Dim >& p ) const { return gradients; }

	// Warning, this function requires the taking of an inverse, which may fail...
	template< class Real >
	static LinearFunction BestFit( const Point< Real , Dim >* points , const V* values , int count )
	{
		LinearFunction lf;
		V constraint[Dim];
		SquareMatrix< Real , Dim > M , Minv;
		M *= 0;
		for( int d=0 ; d<Dim ; d++ ) constraint[d] *= 0;
		for( int i=0 ; i<count ; i++ )
		{
			for( int k=0 ; k<Dim ; k++ ) for( int l=0 ; l<Dim ; l++ ) M( k , l ) += points[i][k] * points[i][l];
			for( int j=0 ; j<count ; j++ ) for( int k=0 ; k<Dim ; k++ ) for( int l=0 ; l<Dim ; l++ ) M( k , l ) -= points[i][k] * points[j][l] / Real( count ); 

			for( int d=0 ; d<Dim ; d++ ) constraint[d] += values[i] * points[i][d];
			for( int j=0 ; j<count ; j++ ) for( int d=0 ; d<Dim ; d++ ) constraint[d] -= values[j] * points[i][d] / Real( count );
		}
		Minv = M.inverse();

		lf *= 0;
		for( int c=0 ; c<Dim ; c++ ) for( int r=0 ; r<Dim ; r++ ) lf.gradients[r] += constraint[c] * Minv( c , r );
		for( int i=0 ; i<count ; i++ )
		{
			for( int d=0 ; d<Dim ; d++ ) lf.offset += lf.gradients[d] * points[i][d];
			lf.offset -= values[i];
		}
		lf.offset /= Real( count );
		return lf;
	}


	//////////////////////////
	// Vector space methods //
	void Add            ( const LinearFunction& lf ) { this->gradient += lf.gradient , offset += lf.offset; }
	void Scale          ( _R s ) { gradients *= s , offset *= s; }
	//////////////////////////
};

template< class Real , int Dim >
struct OrientedPoint
{
	Point< Real , Dim > position , normal;
	template< class Real2 > operator Point< Real2, Dim > ( void ) const { return Point< Real2 , Dim >( position ); }
};

#if FAST_POINT
template< class Real >
struct Point2D
{
	Point2D( void ){ coords[0] = coords[1] = (Real)0.; }
	Point2D( Real x , Real y ){ coords[0] = x , coords[1] = y; }
	template< class Real2 > Point2D( const Point2D< Real2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1]; }
	template< class Real2 > Point2D( const Point< Real2 , 2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1]; }
	template< class Real2 > Point2D& operator = ( const Point2D< Real2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] ; return *this; }
	template< class Real2 > Point2D& operator = ( const Point< Real2 , 2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] ; return *this; }
	Real& operator[] ( int i ) { return coords[i]; }
	const Real& operator[] ( int i ) const { return coords[i]; }

	Point2D& operator += ( const Point2D& p ){ coords[0] += p.coords[0] , coords[1] += p.coords[1] ; return *this; }
	Point2D& operator -= ( const Point2D& p ){ coords[0] -= p.coords[0] , coords[1] -= p.coords[1] ; return *this; }
	Point2D& operator *= ( Real s ){ coords[0] *= s , coords[1] *= s ; return *this; }
	Point2D& operator /= ( Real s ){ coords[0] /= s , coords[1] /= s ; return *this; }

	Point2D operator - ( void ) const { return Point2D( -coords[0] , -coords[1] ); }

	Point2D operator + ( const Point2D& p ) const { return Point2D( coords[0] + p.coords[0] , coords[1] + p.coords[1] ); }
	Point2D operator - ( const Point2D& p ) const { return Point2D( coords[0] - p.coords[0] , coords[1] - p.coords[1] ); }
	Point2D operator * ( Real s ) const { return Point2D( coords[0] * s , coords[1] * s ); }
	Point2D operator / ( Real s ) const { return Point2D( coords[0] / s , coords[1] / s ); }

	Real squareNorm( void ) const { return coords[0]*coords[0] + coords[1]*coords[1]; }
	static Real Dot( const Point2D& p1 , const Point2D& p2 ){ return p1.coords[0]*p2.coords[0] + p1.coords[1]*p2.coords[1]; }
	static Real SquareNorm( const Point2D& p ){ return p.coords[0]*p.coords[0] + p.coords[1]*p.coords[1]; }
	static Real Length( const Point2D& p ){ return (Real)sqrt( SquareNorm(p) ); }
	Real coords[2];
};
template< class Real > Point2D< Real > operator * ( const SquareMatrix< Real , 2 >& M , const Point2D< Real >& v )
{
	return Point2D< Real >( M.coords[0][0] * v[0] + M.coords[1][0] * v[1] , M.coords[0][1] * v[0] + M.coords[1][1] * v[1] );
}
template< class Real >
struct Point3D
{
	Point3D( void ){ coords[0] = coords[1] = coords[2] = (Real)0.; }
	Point3D( Real x , Real y , Real z ){ coords[0] = x , coords[1] = y , coords[2] = z; }
	template< class Real2 > Point3D( const Point3D< Real2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] , coords[2] = (Real)p.coords[2]; }
	template< class Real2 > Point3D( const Point< Real2 , 3 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] , coords[2] = (Real)p.coords[2]; }
	template< class Real2 > Point3D& operator = ( const Point3D< Real2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] , coords[2] = (Real)p.coords[2] ; return *this; }
	template< class Real2 > Point3D& operator = ( const Point< Real2 , 3 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] , coords[2] = (Real)p.coords[2] ; return *this; }
	Real& operator[] ( int i ) { return coords[i]; }
	const Real& operator[] ( int i ) const { return coords[i]; }

	Point3D& operator += ( const Point3D& p ){ coords[0] += p.coords[0] , coords[1] += p.coords[1] , coords[2] += p.coords[2] ; return *this; }
	Point3D& operator -= ( const Point3D& p ){ coords[0] -= p.coords[0] , coords[1] -= p.coords[1] , coords[2] -= p.coords[2] ; return *this; }
	Point3D& operator *= ( Real s ){ coords[0] *= s , coords[1] *= s , coords[2] *= s ; return *this; }
	Point3D& operator /= ( Real s ){ coords[0] /= s , coords[1] /= s , coords[2] /= s ; return *this; }

	Point3D operator - ( void ) const { return Point3D( -coords[0] , -coords[1] , -coords[2] ); }

	Point3D operator + ( const Point3D& p ) const { return Point3D( coords[0] + p.coords[0] , coords[1] + p.coords[1] , coords[2] + p.coords[2] ); }
	Point3D operator - ( const Point3D& p ) const { return Point3D( coords[0] - p.coords[0] , coords[1] - p.coords[1] , coords[2] - p.coords[2] ); }
	Point3D operator * ( Real s ) const { return Point3D( coords[0] * s , coords[1] * s , coords[2] * s ); }
	Point3D operator / ( Real s ) const { return Point3D( coords[0] / s , coords[1] / s , coords[2] / s ); }

	Real squareNorm( void ) const { return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]; }
	static Real Dot( const Point3D& p1 , const Point3D& p2 ){ return p1.coords[0]*p2.coords[0] + p1.coords[1]*p2.coords[1] + p1.coords[2]*p2.coords[2]; }
	static Real SquareNorm( const Point3D& p ){ return p.coords[0]*p.coords[0] + p.coords[1]*p.coords[1] + p.coords[2]*p.coords[2]; }
	static Real Length( const Point3D& p ){ return (Real)sqrt( SquareNorm(p) ); }
	static Point3D CrossProduct( const Point3D& p1 , const Point3D& p2 )
	{
		return Point3D( p1.coords[1]*p2.coords[2] - p1.coords[2]*p2.coords[1], - p1.coords[0]*p2.coords[2] + p1.coords[2]*p2.coords[0] , p1.coords[0]*p2.coords[1] - p1.coords[1]*p2.coords[0] );
	}
	Real coords[3];
};
template< class Real > Point3D< Real > operator * ( const SquareMatrix< Real , 3 >& M , const Point3D< Real >& v )
{
	return Point3D< Real >( M.coords[0][0] * v[0] + M.coords[1][0] * v[1] + M.coords[2][0] * v[2] , M.coords[0][1] * v[0] + M.coords[1][1] * v[1] + M.coords[2][1] * v[2] , M.coords[0][2] * v[0] + M.coords[1][2] * v[1] + M.coords[2][2] * v[2] );
}
#else // !FAST_POINT
template< class Real >
class Point2D : public Point<Real,2>
{
public:
	Point2D(void) : Point<Real,2>(){;}
	Point2D(const Point<Real,2>& p) {   memcpy(this->coords,p.coords,sizeof(Real)*2);   }
	Point2D( Real v1 , Real v2 ) { coords[0] = v1 , coords[1] = v2; }
};

template< class Real >
class Point3D : public Point<Real,3>
{
public:
	using Point<Real, 3>::coords;
	Point3D(void) : Point<Real,3>(){;}
	Point3D(const Point<Real,3>& p) { *this = p; }
	Point3D( Real v1 , Real v2 , Real v3 ) { coords[0] = v1 , coords[1] = v2 , coords[2] = v3; }
	Point3D(const Real *p)  { *this = p; }

	static Point3D CrossProduct( const Point3D& p1 , const Point3D & p2 );

	////////////////////////////////////////////////////////////////////////////
	/*! Allow assignment from another point of lower or equal dimension. If the
	//  RHS dimension is lower, copy only the number of elements in RHS.
	//  @param[in]  p   vector from which to copy
	//  @return     reference to the LHS (this) to allow chaining
	*///////////////////////////////////////////////////////////////////////////
	template<int RHSDim>
	Point3D &operator= (const Point<Real, RHSDim> &p)
	{
		assert(RHSDim <= 3);
		for (int d = 0; d < RHSDim; d++)
			coords[d] = p.coords[d];

		return *this;
	}

	////////////////////////////////////////////////////////////////////////////
	/*! Allow assignment from a flat array
	//  @param[in]  p   the value array to copy
	//  @return     reference to the LHS (this) to allow chaining
	*///////////////////////////////////////////////////////////////////////////
	Point3D &operator=(const Real *p)
	{
		memcpy(coords, p, sizeof(Real) * 3);
		return *this;
	}

	////////////////////////////////////////////////////////////////////////////
	/*! Allow assignment from a point of a different real type
	//  @return     reference to the LHS (this) to allow chaining
	*///////////////////////////////////////////////////////////////////////////
	template<typename Real2>
	Point3D &operator=(const Point3D<Real2> &p)
	{
		Point3D<Real2> p2;
		for (int d = 0; d < 3; d++)
			coords[d] = (Real) p.coords[d];

		return *this;
	}

	////////////////////////////////////////////////////////////////////////////
	/*! Allow addition with a flat array
	//  @param[in]  p   the value array to accumulate
	//  @return     copy of result
	*///////////////////////////////////////////////////////////////////////////
	Point3D operator+(const Real *p)
	{
		Point3D result(*this);
		for (int i = 0; i < 3; ++i)
			result.coords[i] += p[i];
		return result;
	}
};
#endif // FAST_POINT

////////////////////////////////////////////////////////////////////////////////
/*! Creates a complex 3d point: p1 + i p2
//  @param[in]  p1      real part
//  @param[in]  p1      imaginary part
//  @return     complex 3d point
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
inline Point3D<std::complex<Real> > ComplexPoint3D(const Point3D<Real> &p1
	, const Point3D<Real> &p2)
{
	return Point3D<std::complex<Real> >(std::complex<Real>(p1[0], p2[0])
		, std::complex<Real>(p1[1], p2[1])
		, std::complex<Real>(p1[2], p2[2]));
}

////////////////////////////////////////////////////////////////////////////////
/*! Extracts the real vector part from a complex vector
//  @param[in]  vector to process
//  @return     real part
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
inline Point3D<Real> real(const Point3D<std::complex<Real> > &p)
{
	return Point3D<Real>(real(p[0]), real(p[1]), real(p[2]));
}

////////////////////////////////////////////////////////////////////////////////
/*! Extracts the imaginary vector part from a complex vector
//  @param[in]  vector to process
//  @return     imaginary part
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
inline Point3D<Real> imag(const Point3D<std::complex<Real> > &p)
{
	return Point3D<Real>(imag(p[0]), imag(p[1]), imag(p[2]));
}

////////////////////////////////////////////////////////////////////////////////
/*! Separately normalizes both the real and the imaginary part of a complex 3D
//  vector
//  @param[inout]   p   the vector to normalize
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
inline void NormalizeParts(Point3D<std::complex<Real> > &p)
{
	Point3D<Real> realPart = real(p);
	Point3D<Real> imagPart = imag(p);
	realPart /= sqrt(SquareNorm(realPart));
	imagPart /= sqrt(SquareNorm(imagPart));
	p = ComplexPoint3D(realPart, imagPart);
}

////////////////////////////////////////////////////////////////////////////////
/*! Saves on a sqrt() operation when the vector is unit.
//  @param[inout]   p   the vector to normalize
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
void safe_normalize(Point3D<Real> &p)
{
	Real normSq = SquareNorm(p);
	if ((normSq != Real(0)) && (std::abs(1 - normSq) >= 1.0e-16))
		p *= (1 / sqrt(normSq));
}

////////////////////////////////////////////////////////////////////////////////
/*! Linearly interpolates between two points
//  @param[in]  alpha       blending factor (alpha = 0 ==> p1)
//  @param[in]  p1, p2      points to interpolate
//  @return     interpolation result
*///////////////////////////////////////////////////////////////////////////////
template<typename Num>
Point3D<Num> LinInterp(Num alpha, const Point3D<Num> &p1,
	const Point3D<Num> &p2)
{
	return (1 - alpha) * p1 + alpha * p2;
}

template< class V , class _R = typename V::R >
class Gradient3D : public Gradient< V , 3 , _R >
{
public:
	Gradient3D( void ) : Gradient< V , 3 , _R >( ){ ; }
	Gradient3D( const Gradient< V , 3 , _R >& lf ) { for( int d=0 ; d<3 ; d++ ) this->gradients[d] = lf.gradients[d]; }
};

////////////////////////////////////////////////////////////////////////////////
/*! Class to represent a linear function mapping R^3 to V
//  (as I understand it--JP)
//  @tparam  V  the range type
//  @tparam _R  the floating point representation
*///////////////////////////////////////////////////////////////////////////////
template< class V , class _R = typename V::R >
class LinearFunction3D : public LinearFunction< V , 3 , _R >
{
public:
	LinearFunction3D( void ) : LinearFunction< V , 3 , _R >( ){ ; }
	LinearFunction3D( const LinearFunction< V , 3 , _R >& lf )
	{
		for( int d=0 ; d<3 ; d++ ) this->gradients[d] = lf.gradients[d];
		this->offset = lf.offset;
	}
	template< class Real >
	static LinearFunction3D GetInterpolant( const Point3D< Real >* vertices , const V* values , Point3D< Real > n )
	{
		Point3D< Real > p1 = Point3D< Real >( vertices[1] - vertices[0] );
		Point3D< Real > p2 = Point3D< Real >( vertices[2] - vertices[0] );
		Point3D< Real > v1 = Point3D< Real >::CrossProduct( p2 , n  );
		Point3D< Real > v2 = Point3D< Real >::CrossProduct( p1 , n  );

		Real d1 = Point3D< Real >::Dot( v1 , p1 );
		Real d2 = Point3D< Real >::Dot( v2 , p2 );
		if( !d1 || !d2 )
		{
			LinearFunction3D< V , _R > lf;
			lf.gradients *= 0;
#if PAN_FIX
			lf.offset = -( values[0] + values[1] + values[2] ) / 3;
#else // !PAN_FIX
			lf.offset =  ( values[0] + values[1] + values[2] ) / 3;
#endif // PAN_FIX
			return lf;
		}
		else
		{
			v1 /= d1;
			v2 /= d2;
		}

		LinearFunction3D< V , _R > lf;
		lf.offset =
			- ( values[0] * ( 1.0 + Point3D< Real >::Dot( vertices[0] , v1 ) + Point3D< Real >::Dot( vertices[0] , v2 ) )
			- values[1] * Point3D< Real >::Dot( vertices[0] , v1 )
			- values[2] * Point3D< Real >::Dot( vertices[0] , v2 ) );
		for( int d=0 ; d<3 ; d++ ) lf.gradients[d] = - values[0]*v1[d] - values[0]*v2[d] + values[1]*v1[d] + values[2]*v2[d];
		return lf;
	}
	template< class Real >
	static LinearFunction3D GetInterpolant( const Point3D< Real >& v1 ,
		const Point3D< Real >& v2 ,
		const Point3D< Real >& v3 ,
		const V& s1 , const V& s2 , const V s3 ,
		Point3D< Real > n )
	{
		Point3D< Real > p1 = Point3D< Real >( v2 - v1 );
		Point3D< Real > p2 = Point3D< Real >( v3 - v1 );
		Point3D< Real > _v1 = Point3D< Real >::CrossProduct( p2 , n  );
		Point3D< Real > _v2 = Point3D< Real >::CrossProduct( p1 , n  );

		Real d1 = Point3D< Real >::Dot( _v1 , p1 );
		Real d2 = Point3D< Real >::Dot( _v2 , p2 );
		if( !d1 || !d2 )
		{
			LinearFunction3D< V , _R > lf;
			lf.gradients *= 0;
#if PAN_FIX
			lf.offset = -( s1 + s2 + s3 ) / _R(3);
#else // !PAN_FIX
			lf.offset = ( s1 + s2 + s3 ) / _R(3);
#endif // PAN_FIX
			return lf;
		}
		else
		{
			_v1 /= d1;
			_v2 /= d2;
		}

		LinearFunction3D< V , _R > lf;
		lf.offset = - ( s1 * ( Real(1.0)
			+ Point3D< Real >::Dot( v1 , _v1 )
			+ Point3D< Real >::Dot( v1 , _v2 ))
			- s2 * Point3D< Real >::Dot( v1 , _v1 )
			- s3 * Point3D< Real >::Dot( v1 , _v2 ) );
		for( int d=0 ; d<3 ; d++ ) lf.gradients[d] = - s1*_v1[d] - s1*_v2[d] + s2*_v1[d] + s3*_v2[d];
		return lf;
	}
	template< class Real >
	static LinearFunction3D GetInterpolant( const Point3D< Real >* vertices , const V* values )
	{
		Point3D< Real > p1 = Point3D< Real >( vertices[1] - vertices[0] );
		Point3D< Real > p2 = Point3D< Real >( vertices[2] - vertices[0] );
		Point3D< Real > n  = Point3D< Real >::CrossProduct( p1 , p2 );
		return GetInterpolant( vertices , values , n );
	}
	template< class Real >
	static LinearFunction3D GetInterpolant( const Point3D< Real >& v1 , const Point3D< Real >& v2 , const Point3D< Real >& v3 , const V& s1 , const V& s2 , const V& s3 )
	{
		Point3D< Real > p1 = Point3D< Real >( v2 - v1 );
		Point3D< Real > p2 = Point3D< Real >( v3 - v1 );
		Point3D< Real > n  = Point3D< Real >::CrossProduct( p1 , p2 );
		return GetInterpolant( v1 ,v2 , v3 , s1 , s2 , s3 , n );
	}
};

template<class Real>
class OrientedPoint2D : public OrientedPoint<Real,2>{;};
template<class Real>
class OrientedPoint3D : public OrientedPoint<Real,3>{;};

template<class Real>
class XForm4x4 : public SquareMatrix<Real,4>
{
public:
	XForm4x4(void) : SquareMatrix<Real,4>(){;}
	XForm4x4(const SquareMatrix<Real,4>& xForm) {   memcpy(this->coords,xForm.coords,sizeof(Real)*4*4); };
};
template<class Real>
class XForm3x3 : public SquareMatrix<Real,3>
{
public:
	XForm3x3(void) : SquareMatrix<Real,3>(){;}
	XForm3x3(const SquareMatrix<Real,3>& xForm) {   memcpy(this->coords,xForm.coords,sizeof(Real)*3*3); };
};
template<class Real>
class XForm2x2 : public SquareMatrix<Real,2>
{
public:
	XForm2x2(void) : SquareMatrix<Real,2>(){;}
	XForm2x2(const SquareMatrix<Real,2>& xForm) {   memcpy(this->coords,xForm.coords,sizeof(Real)*2*2); };
};

template< class Real > Point2D< Real > RandomDiskPoint( void );
template< class Real > Point3D< Real > RandomBallPoint( void );

template< class Real > Point2D< Real > RandomCirclePoint( void );
template< class Real > Point3D< Real > RandomSpherePoint( void );

template<class Real>
XForm3x3<Real> RotationMatrix( Real a , Real b , Real c , Real d );

template<class Real>
XForm3x3<Real> RotationMatrix( const Point3D<Real>& axis , const Real& angle );

template<class Real>
XForm3x3<Real> RandomRotationMatrix( void );

template<class Real>
double Length(const Point3D<Real>& p);

template<class Real>
double SquareLength(const Point3D<Real>& p);

template<class Real>
double Distance(const Point3D<Real>& p1,const Point3D<Real>& p2);

template <class Real>
void CrossProduct(const Point3D<Real>& p1,const Point3D<Real>& p2,Point3D<Real>& p);

template <class Real>
Point3D<Real> CrossProduct(const Point3D<Real>& p1,const Point3D<Real>& p2)
{
	// Avoid Point3D<Real>:: in front of every call
	return Point3D<Real>::CrossProduct(p1, p2);
}

template <class Real>
Point3D<Real> CrossProduct(const Point<Real, 3> &p1, const Point<Real, 3>& p2)
{
	// Avoid Point3D<Real>:: in front of every call
	return Point3D<Real>::CrossProduct(p1, p2);
}

template <class Real>
Real Dot(const Point3D<Real>& p1,const Point3D<Real>& p2)
{
	// Avoid Point3D<Real>:: in front of every call
	return Point3D<Real>::Dot(p1, p2);
}

template <class Real>
Real SquareNorm(const Point3D<Real>& p)
{
	// Avoid Point3D<Real>:: in front of every call
	return Point3D<Real>::SquareNorm(p);
}

template<class Real, int Dim>
Real SquareNorm(const Point<Real, Dim> &p)
{
	// Avoid Point<Real, Dim>:: in front of every call
	return Point<Real, Dim>::SquareNorm(p);
}

template< class Real >
double SquareNorm( std::complex< Real > c )
{
	return c.real()*c.real() + c.imag()*c.imag();
}

template<class Real>
void Transform(const XForm4x4<Real>& xForm,const Point3D<Real>& p,Point3D<Real>& q);

template<class Real>
void TransformNoTranslate(const XForm4x4<Real>& xForm,const Point3D<Real>& p,Point3D<Real>& q);


template< class Real >
void BarycentricCoordinates( const Point3D< Real >& p , const Point3D< Real >& v1 , const Point3D< Real >& v2, const Point3D< Real >& v3 , Real& a0 , Real& a1 , Real& a2 );

////////////////////////////////////////////////////////////////////////////////
/*! Uses a barycentric coordinate vector to interpolate three data values
//  @param[in]  coords      bary centric coordinates
//  @param[in]  d0, d1, d2  data values to interpolate
//  @return     interpolated data value
*///////////////////////////////////////////////////////////////////////////////
template<typename BaryCoords, typename DataType>
inline DataType BarycentricInterpolate(const BaryCoords &coords
	, const DataType &d0, const DataType &d1, const DataType &d2);

////////////////////////////////////////////////////////////////////////////////
/*! Computes a triangle's inscribed circle
//  http://en.wikipedia.org/wiki/Incircle
//  @param[in]  p0, p1, p2      triangle vertex positions
//  @param[out] center          incircle center
//  @param[out] radius          incircle radius
*///////////////////////////////////////////////////////////////////////////////
template<typename PointType, typename Real>
inline void Incircle(const PointType &p0, const PointType &p1,
	const PointType &p2, PointType &center, Real &radius);

////////////////////////////////////////////////////////////////////////////////
/*! Computes a triangle's circumscribed circle
//  http://en.wikipedia.org/wiki/Circumscribed_circle
//  @param[in]  p0, p1, p2      triangle vertex positions
//  @param[in]  tri             triangle to process
//  @param[out] center          incircle center
*///////////////////////////////////////////////////////////////////////////////
template<typename PointType, typename Real>
inline void Circumcircle(const PointType &p0, const PointType &p1
	, const PointType &p2, PointType &center, Real &radius);

class Edge
{
public:
	double p[2][2];
	double Length(void) const{
		double d[2];
		d[0]=p[0][0]-p[1][0];
		d[1]=p[0][1]-p[1][1];

		return sqrt(d[0]*d[0]+d[1]*d[1]);
	}
};
class Triangle
{
public:
	double p[3][3];
	double Area(void) const{
		double v1[3],v2[3],v[3];
		for(int d=0;d<3;d++){
			v1[d]=p[1][d]-p[0][d];
			v2[d]=p[2][d]-p[0][d];
		}
		v[0]= v1[1]*v2[2]-v1[2]*v2[1];
		v[1]=-v1[0]*v2[2]+v1[2]*v2[0];
		v[2]= v1[0]*v2[1]-v1[1]*v2[0];
		return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2;
	}
	double AspectRatio(void) const{
		double d=0;
		int i,j;
		for(i=0;i<3;i++){
			for(i=0;i<3;i++)
				for(j=0;j<3;j++){d+=(p[(i+1)%3][j]-p[i][j])*(p[(i+1)%3][j]-p[i][j]);}
		}
		return Area()/d;
	}
};
class EdgeIndex
{
public:
	int v[2];
	int& operator[] ( int idx ) { return v[idx]; }
	const int& operator[] ( int idx ) const { return v[idx]; }
};

class SeamEdge : public EdgeIndex
{
public:
	// NOTE: edge is from index 1 to index 0
	/** The rotation across the seam edge */
	unsigned char rot;
	/** The face with which this seam was associated in the obj file
	*  (used to make writing the seams back to an obj file easy) */
	unsigned int face_idx;
};

class ConeVertex
{
public:
	/** Index of vertex in the vertex list */
	int v;
	/** Vertex index quadrupled so it becomes integer (4 for non-singular
	*  vertices) */
	char two_cone_angle_over_pi;

	ConeVertex()
		: v(-1), two_cone_angle_over_pi(0) { }
	ConeVertex(int i, int doubleConeAngleOverPi)
		: v(i), two_cone_angle_over_pi(doubleConeAngleOverPi)
	{
	}
};

class TriangleIndex
{
protected:
	unsigned int v[3];
public:
	TriangleIndex( void ) { v[0] = v[1] = v[2] = 0; }
	TriangleIndex( unsigned int v0 , unsigned int v1 , unsigned int v2 ){ v[0] = v0; v[1] = v1; v[2] = v2; }
	unsigned int &operator[]( unsigned int idx ) { return v[idx]; }
	unsigned int  operator[]( unsigned int idx ) const { return v[idx]; }
};

template< class Data >
class TriangleIndexWithData : public TriangleIndex
{
public:
	Data data[3];
	TriangleIndexWithData( void ) { v[0] = v[1] = v[2] = 0 ;}
	TriangleIndexWithData( unsigned int v0 , unsigned int v1 , unsigned int v2 ){ v[0] = v0 , v[1] = v1 , v[2] = v2; }
	TriangleIndexWithData( unsigned int v0 , unsigned int v1 , unsigned int v2 , Data d0 , Data d1 , Data d2 ){ v[0] = v0 , v[1] = v1 , v[2] = v2 , data[0] = d0 , data[1] = d1 , data[2] = d2; }
};

////////////////////////////////////////////////////////////////////////////////
/*! @class ParameterizedTriangleIndex
//  Holds holds both the vertex and uv indices of each triangle corner.
*///////////////////////////////////////////////////////////////////////////////
class ParameterizedTriangleIndex : public TriangleIndex
{
	using TriangleIndex::v;
	int m_uv[3];
public: 
	ParameterizedTriangleIndex()
		: TriangleIndex() { m_uv[0] = m_uv[1] = m_uv[2] = 0; }

	ParameterizedTriangleIndex(int v0, int v1, int v2)
		: TriangleIndex(v0, v1, v2) { m_uv[0] = m_uv[1] = m_uv[2] = 0; }
	ParameterizedTriangleIndex(const TriangleIndex &t)
		: TriangleIndex(t) { m_uv[0] = m_uv[1] = m_uv[2] = 0; }
	int &uv(int idx)        { return m_uv[idx]; }
	int  uv(int idx) const  { return m_uv[idx]; }
	////////////////////////////////////////////////////////////////////////////
	/*! Typecast to TriangleIndex strips uv indices
	*///////////////////////////////////////////////////////////////////////////
	operator TriangleIndex() const
	{
		TriangleIndex out(v[0], v[1], v[2]);
		return out;
	}
	operator TriangleIndex&()
	{
		return *((TriangleIndex *) this);
	}
};

template <class Real>
class MinimalAreaTriangulation
{
	double* bestTriangulation;
	int* midPoint;
	double GetArea( const int& i , const int& j , const std::vector<Point3D<Real> >& vertices );
	void GetTriangulation( const int& i , const int& j , const std::vector<Point3D<Real> >& vertices,std::vector<TriangleIndex>& triangles , int& idx);
public:
	MinimalAreaTriangulation(void);
	~MinimalAreaTriangulation(void);
	double GetArea(const std::vector<Point3D<Real> >& vertices);
	void GetTriangulation( const std::vector<Point3D<Real> >& vertices , std::vector<TriangleIndex>& triangles );
};

template<class Vertex>
class Mesh
{
public:
	std::vector<Vertex> vertices;
	std::vector<std::vector<int> > polygons;
};

template< class Real >
void SplitTriangle( std::vector< Point3D< Real > >& vertices , TriangleIndex triangle , Point3D< Real > pNormal , Real pOffset ,
	std::vector< TriangleIndex >& backTriangles , std::vector< TriangleIndex >& frontTriangles );
template< class Real >
void SplitTriangle( std::vector< Point3D< Real > >& vertices , TriangleIndex triangle , int direction , Real offset ,
	std::vector< TriangleIndex >& backTriangles , std::vector< TriangleIndex >& frontTriangles );
#include "Geometry.inl"
#endif // GEOMETRY_INCLUDED
