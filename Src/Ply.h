/*

   Header for PLY polygon files.

   - Greg Turk, March 1994

   A PLY file contains a single polygonal _object_.

   An object is composed of lists of _elements_.  Typical elements are
   vertices, faces, edges and materials.

   Each type of element for a given object has one or more _properties_
   associated with the element type.  For instance, a vertex element may
   have as properties three floating-point values x,y,z and three unsigned
   chars for red, green and blue.

   ---------------------------------------------------------------

   Copyright (c) 1994 The Board of Trustees of The Leland Stanford
   Junior University.  All rights reserved.   

   Permission to use, copy, modify and distribute this software and its   
   documentation for any purpose is hereby granted without fee, provided   
   that the above copyright notice and this permission notice appear in   
   all copies of this software and that you do not sell the software.   

   THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
   EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
   WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.   

*/

// JP: Disable the "deprecated conversion from string constant to char *"
#ifndef WIN32
#pragma GCC diagnostic ignored "-Wwrite-strings"
#endif // WIN32

#ifndef PLY_INCLUDED
#define PLY_INCLUDED
#include <vector>
#include "Geometry.h"
#include "PlyFile.h"
#include "Grid.h"

template< class Real >
class PlyVertex : public VectorSpace< Real , PlyVertex< Real > >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add	(const PlyVertex& p)	{	point+=p.point;	}
	void Scale	(float s)				{	point*=s;		}
	//////////////////////////

	const static int ReadComponents=3;
	const static int WriteComponents=3;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];

	Point3D< Real > point;

	operator Point3D< Real >& ()					{ return point; }
	operator const Point3D< Real >& () const		{ return point; }
	template< class Real2 > operator Point3D< Real2 > ( ) const { return Point3D< Real2 >( point ); }
	PlyVertex( void )								{ point.coords[0]=point.coords[1]=point.coords[2]=0; }
	PlyVertex( const Point3D< Real >& p )			{ point=p; }

	template<class Real2>
	PlyVertex xForm( const XForm4x4<Real2>& xForm ) const
	{
		PlyVertex v;
		v.point  = xForm ( point );
		return v;
	}
};
template<>
PlyProperty PlyVertex< float >::ReadProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex,point.coords[2])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyVertex< float >::WriteProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex,point.coords[2])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyVertex< double >::ReadProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex,point.coords[2])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyVertex< double >::WriteProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex,point.coords[2])), 0, 0, 0, 0}
};
template< class Real >
class PlyVertex2D : public VectorSpace< Real , PlyVertex2D< Real > >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add	(const PlyVertex2D& p)	{	point+=p.point;	}
	void Scale	(float s)				{	point*=s;		}
	//////////////////////////

	const static int ReadComponents=2;
	const static int WriteComponents=2;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];

	Point2D< Real > point;

	operator Point2D< Real >& ()					{ return point; }
	operator const Point2D< Real >& () const		{ return point; }
	template< class Real2 > operator Point2D< Real2 > ( ) const { return Point2D< Real2 >( point ); }
	PlyVertex2D( void )								{ point.coords[0]=point.coords[1]=0; }
	PlyVertex2D( const Point2D< Real >& p )			{ point=p; }
};
template<>
PlyProperty PlyVertex2D< float >::ReadProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex2D,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex2D,point.coords[1])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyVertex2D< float >::WriteProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex2D,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertex2D,point.coords[1])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyVertex2D< double >::ReadProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex2D,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex2D,point.coords[1])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyVertex2D< double >::WriteProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex2D,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyVertex2D,point.coords[1])), 0, 0, 0, 0}
};

template< class Real >
class PlyTextureVertex : public VectorSpace< Real , PlyTextureVertex< Real > >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add	(const PlyTextureVertex& p)	{	point+=p.point;	texture+=p.texture;	}
	void Scale	(Real s)					{	point*=s;		texture*=s;	    	}
	//////////////////////////
	const static int ReadComponents=7;
	const static int WriteComponents=5;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];

	Point3D< Real > point;
	Point2D< Real > texture;

	operator Point3D<Real>& ()					{return point;}
	operator const Point3D<Real>& () const		{return point;}
	template< class Real2 > operator Point3D< Real2 > ( ) const { return Point3D< Real2 >( point ); }
	PlyTextureVertex( void )                   { point[0] = point[1] = point[2] = texture[0] = texture[1] = 0; }
	PlyTextureVertex( const Point3D<Real>& p ) { point = p; }
	PlyTextureVertex( const Point3D<Real>& p , const Point2D< Real >& t ) { point = p , texture = t; }
};
template<>
PlyProperty PlyTextureVertex< float >::ReadProperties[]=
{
	{ "x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,point.coords[0])), 0, 0, 0, 0},
	{ "y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,point.coords[1])), 0, 0, 0, 0},
	{ "z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,point.coords[2])), 0, 0, 0, 0},
	{ "u", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,texture.coords[0])),0, 0, 0, 0},
	{ "v", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,texture.coords[1])),0, 0, 0, 0},
	{ "texture_u", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,texture.coords[0])),0, 0, 0, 0},
	{ "texture_v", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,texture.coords[1])),0, 0, 0, 0},
};
template<>
PlyProperty PlyTextureVertex< float >::WriteProperties[]=
{
	{ "x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,point.coords[0])), 0, 0, 0, 0},
	{ "y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,point.coords[1])), 0, 0, 0, 0},
	{ "z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,point.coords[2])), 0, 0, 0, 0},
	{ "u", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,texture.coords[0])),0, 0, 0, 0},
	{ "v", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyTextureVertex,texture.coords[1])),0, 0, 0, 0},
};
template<>
PlyProperty PlyTextureVertex< double >::ReadProperties[]=
{
	{ "x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,point.coords[0])), 0, 0, 0, 0},
	{ "y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,point.coords[1])), 0, 0, 0, 0},
	{ "z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,point.coords[2])), 0, 0, 0, 0},
	{ "u", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,texture.coords[0])),0, 0, 0, 0},
	{ "v", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,texture.coords[1])),0, 0, 0, 0},
	{ "texture_u", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,texture.coords[0])),0, 0, 0, 0},
	{ "texture_v", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,texture.coords[1])),0, 0, 0, 0},
};
template<>
PlyProperty PlyTextureVertex< double >::WriteProperties[]=
{
	{ "x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,point.coords[0])), 0, 0, 0, 0},
	{ "y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,point.coords[1])), 0, 0, 0, 0},
	{ "z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,point.coords[2])), 0, 0, 0, 0},
	{ "u", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,texture.coords[0])),0, 0, 0, 0},
	{ "v", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyTextureVertex,texture.coords[1])),0, 0, 0, 0},
};

template< class Real >
class PlyValueVertex : public VectorSpace< Real , PlyValueVertex< Real > >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add	(const PlyValueVertex& p)	{	point+=p.point;	value+=p.value;	}
	void Scale	(float s)					{	point*=s;		value*=s;		}
	//////////////////////////
	const static int ReadComponents=4;
	const static int WriteComponents=4;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];

	Point3D< Real > point;
	Real  value;

	operator Point3D< Real >& ()				{return point;}
	operator const Point3D< Real >& () const	{return point;}
	template< class Real2 > operator Point3D< Real2 > ( ) const { return Point3D< Real2 >( point ); }
	PlyValueVertex(void)					{point.coords[0]=point.coords[1]=point.coords[2]=0;value=0;}
	PlyValueVertex(const Point3D<Real>& p)	{point=p;}

	template<class Real2>
	PlyValueVertex xForm(const XForm4x4<Real2>& xForm) const
	{
		PlyValueVertex v;
		v.value=value;
		v.point  = xForm ( point );
		return v;
	}
};
template<>
PlyProperty PlyValueVertex< float >::ReadProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyValueVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyValueVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyValueVertex,point.coords[2])), 0, 0, 0, 0},
	{"value", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyValueVertex,value)), 0, 0, 0, 0}
};
template<>
PlyProperty PlyValueVertex< float >::WriteProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyValueVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyValueVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyValueVertex,point.coords[2])), 0, 0, 0, 0},
	{"value", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyValueVertex,value)), 0, 0, 0, 0}
};
template<>
PlyProperty PlyValueVertex< double >::ReadProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyValueVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyValueVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyValueVertex,point.coords[2])), 0, 0, 0, 0},
	{"value", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyValueVertex,value)), 0, 0, 0, 0}
};
template<>
PlyProperty PlyValueVertex< double >::WriteProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyValueVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyValueVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyValueVertex,point.coords[2])), 0, 0, 0, 0},
	{"value", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyValueVertex,value)), 0, 0, 0, 0}
};

template< class Real >
class PlyOrientedVertex : public VectorSpace< Real , PlyOrientedVertex< Real > >
{
public:
    typedef Real real_type;

	//////////////////////////
	// Vector space methods //
	void Add	(const PlyOrientedVertex& p)	{	point+=p.point;	normal+=p.normal;	}
	void Scale	( Real s)						{	point*=s;		normal*=s;			}
	//////////////////////////
	const static int ReadComponents=6;
	const static int WriteComponents=6;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];

	Point3D< Real > point,normal;

	operator Point3D<Real>& ()					{return point;}
	operator const Point3D<Real>& () const		{return point;}
	template< class Real2 > operator Point3D< Real2 > ( ) const { return Point3D< Real2 >( point ); }
	PlyOrientedVertex(void)						{point.coords[0]=point.coords[1]=point.coords[2]=normal.coords[0]=normal.coords[1]=normal.coords[2]=0;}
	PlyOrientedVertex(const Point3D<Real>& p)	{point=p;}

    /*  BROKEN
     *  fix: xForm2 undeclared
	template<class Real2>
	PlyOrientedVertex xForm(const XForm4x4<Real2>& xForm) const
	{
		XForm3x3<Real2> temp=Matrix<Real,3,3>(xForm2);
		PlyOrientedVertex v;
		v.point  = xForm ( point );
		v.normal = temp * normal;
		return v;
	}
    */
};
template<>
PlyProperty PlyOrientedVertex< float >::ReadProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,normal.coords[2])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyOrientedVertex< float >::WriteProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedVertex,normal.coords[2])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyOrientedVertex< double >::ReadProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,normal.coords[2])), 0, 0, 0, 0}
};
template<>
PlyProperty PlyOrientedVertex< double >::WriteProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedVertex,normal.coords[2])), 0, 0, 0, 0}
};

template< class Real >
class PlyColorVertex : public VectorSpace< Real , PlyColorVertex< Real > >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add	(const PlyColorVertex& p)	{	point+=p.point;	color+=p.color;	}
	void Scale	(Real s)					{	point*=s;		color*=s;		}
	//////////////////////////
	const static int ReadComponents=9;
	const static int WriteComponents=6;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];

	Point3D<Real> point,color;

	operator Point3D<Real>& ()					{return point;}
	operator const Point3D<Real>& () const		{return point;}
	template< class Real2 > operator Point3D< Real2 > ( ) const { return Point3D< Real2 >( point ); }
	PlyColorVertex(void)						{point.coords[0]=point.coords[1]=point.coords[2]=0,color.coords[0]=color.coords[1]=color.coords[2]=0;}
	PlyColorVertex(const Point3D<Real>& p)		{point=p;}

	template<class Real2>
	PlyColorVertex xForm(const XForm4x4<Real2>& xForm) const
	{
		PlyColorVertex< Real > v;
		v.color.coords[0]=color.coords[0];
		v.color.coords[1]=color.coords[1];
		v.color.coords[2]=color.coords[2];
		v.point  = xForm ( point );
		return v;
	}
};
template<>
PlyProperty PlyColorVertex< float >::ReadProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex,point.coords[2])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[0])),0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[1])),0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[2])),0, 0, 0, 0},
	{"diffuse_red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[0])),0, 0, 0, 0},
	{"diffuse_green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[1])),0, 0, 0, 0},
	{"diffuse_blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[2])),0, 0, 0, 0},
};
template<>
PlyProperty PlyColorVertex< float >::WriteProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex,point.coords[2])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[0])),0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[1])),0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex,color.coords[2])),0, 0, 0, 0}
};
template<>
PlyProperty PlyColorVertex< double >::ReadProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex,point.coords[2])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[0])),0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[1])),0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[2])),0, 0, 0, 0},
	{"diffuse_red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[0])),0, 0, 0, 0},
	{"diffuse_green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[1])),0, 0, 0, 0},
	{"diffuse_blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[2])),0, 0, 0, 0},
};
template<>
PlyProperty PlyColorVertex< double >::WriteProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex,point.coords[2])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[0])),0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[1])),0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex,color.coords[2])),0, 0, 0, 0}
};

template< class Real >
class PlyColorVertex2D : public VectorSpace< Real , PlyColorVertex2D< Real > >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add	(const PlyColorVertex2D& p)	{	point+=p.point;	color+=p.color;	}
	void Scale	(Real s)					{	point*=s;		color*=s;		}
	//////////////////////////
	const static int ReadComponents=8;
	const static int WriteComponents=5;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];

	Point2D< Real > point;
	Point3D< Real > color;

	operator Point2D<Real>& ()					{return point;}
	operator const Point2D<Real>& () const		{return point;}
	template< class Real2 > operator Point2D< Real2 > ( ) const { return Point2D< Real2 >( point ); }
	PlyColorVertex2D(void)						{point.coords[0]=point.coords[1]=0,color.coords[0]=color.coords[1]=color.coords[2]=0;}
	PlyColorVertex2D(const Point2D<Real>& p)	{point=p;}
};
template<>
PlyProperty PlyColorVertex2D< float >::ReadProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex2D,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex2D,point.coords[1])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[0])),0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[1])),0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[2])),0, 0, 0, 0},
	{"diffuse_red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[0])),0, 0, 0, 0},
	{"diffuse_green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[1])),0, 0, 0, 0},
	{"diffuse_blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[2])),0, 0, 0, 0},
};
template<>
PlyProperty PlyColorVertex2D< float >::WriteProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex2D,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyColorVertex2D,point.coords[1])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[0])),0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[1])),0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyColorVertex2D,color.coords[2])),0, 0, 0, 0}
};
template<>
PlyProperty PlyColorVertex2D< double >::ReadProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,point.coords[1])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[0])),0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[1])),0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[2])),0, 0, 0, 0},
	{"diffuse_red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[0])),0, 0, 0, 0},
	{"diffuse_green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[1])),0, 0, 0, 0},
	{"diffuse_blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[2])),0, 0, 0, 0},
};
template<>
PlyProperty PlyColorVertex2D< double >::WriteProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,point.coords[1])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[0])),0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[1])),0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyColorVertex2D,color.coords[2])),0, 0, 0, 0}
};

template< class Real >
class PlyOrientedColorVertex : public VectorSpace< Real , PlyOrientedColorVertex< Real > >
{
public:
	//////////////////////////
	// Vector space methods //
	void Add	( const PlyOrientedColorVertex& p )	{	point+=p.point;	normal+=p.normal;	color+=p.color;	}
	void Scale	( Real s )							{	point*=s;		normal*=s;			color*=s;		}
	//////////////////////////
	const static int ReadComponents=12;
	const static int WriteComponents=9;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];

	Point3D< Real > point , normal , color;

	operator Point3D< Real >& ()					{return point;}
	operator const Point3D< Real >& () const		{return point;}
	template< class Real2 > operator Point3D< Real2 > ( ) const { return Point3D< Real2 >( point ); }
	PlyOrientedColorVertex( void )						{ point.coords[0]=point.coords[1]=point.coords[2]=normal.coords[0]=normal.coords[1]=normal.coords[2]=0,color.coords[0]=color.coords[1]=color.coords[2]=0; }
	PlyOrientedColorVertex( const Point3D< Real>& p )	{ point=p; }

	template< class Real2 >
	PlyOrientedColorVertex xForm(const XForm4x4< Real2 >& xForm) const
	{
		XForm3x3< Real2 > temp=XForm3x3< Real2 >( Matrix< Real2 , 3 , 3 >(xForm) );
		PlyOrientedColorVertex v;
		v.color.coords[0]=color.coords[0];
		v.color.coords[1]=color.coords[1];
		v.color.coords[2]=color.coords[2];
		v.point  = xForm ( point );
		v.normal = temp *normal;
		return v;
	}
};
template<>
PlyProperty PlyOrientedColorVertex< float >::ReadProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,normal.coords[2])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[0])),	0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[1])),	0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[2])),	0, 0, 0, 0},
	{"diffuse_red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[0])),	0, 0, 0, 0},
	{"diffuse_green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[1])),	0, 0, 0, 0},
	{"diffuse_blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[2])),	0, 0, 0, 0}
};
template<>
PlyProperty PlyOrientedColorVertex< float >::WriteProperties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,normal.coords[2])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[0])),	0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[1])),	0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_FLOAT, int(offsetof(PlyOrientedColorVertex,color.coords[2])),	0, 0, 0, 0},
};

template<>
PlyProperty PlyOrientedColorVertex< double >::ReadProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,normal.coords[2])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[0])),	0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[1])),	0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[2])),	0, 0, 0, 0},
	{"diffuse_red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[0])),	0, 0, 0, 0},
	{"diffuse_green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[1])),	0, 0, 0, 0},
	{"diffuse_blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[2])),	0, 0, 0, 0}
};
template<>
PlyProperty PlyOrientedColorVertex< double >::WriteProperties[]=
{
	{"x", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_DOUBLE, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,normal.coords[2])), 0, 0, 0, 0},
	{"red",		PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[0])),	0, 0, 0, 0},
	{"green",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[1])),	0, 0, 0, 0},
	{"blue",	PLY_UCHAR, PLY_DOUBLE, int(offsetof(PlyOrientedColorVertex,color.coords[2])),	0, 0, 0, 0},
};

template< class Real >
struct PlyVFFace
{
	unsigned int nr_vertices;
	int *vertices;
	Point3D< Real > v;
	PlyVFFace( void ){ vertices=NULL , nr_vertices=0 ;}
	~PlyVFFace( void ){ resize(0); }
	void resize( unsigned int count )
	{
		if( vertices ) free( vertices );
		vertices = NULL , nr_vertices = 0;
		if( count ) vertices = (int*)malloc( sizeof(int)*count ) , nr_vertices = count;
	}
	int& operator[] ( int idx ){ return vertices[idx]; }
	const int& operator[] ( int idx ) const { return vertices[idx]; }
	int size( void ) const { return nr_vertices; }

	const static int ReadComponents=4;
	const static int WriteComponents=4;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];
};
template<>
PlyProperty PlyVFFace< float >::ReadProperties[] =
{
	{ "vertex_indices" , PLY_INT , PLY_INT , offsetof( PlyVFFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , nr_vertices ) } ,
	{ "vx" , PLY_FLOAT , PLY_FLOAT , (int)offsetof( PlyVFFace , v.coords[0] ) , 0 , 0 , 0 , 0 } ,
	{ "vy" , PLY_FLOAT , PLY_FLOAT , (int)offsetof( PlyVFFace , v.coords[1] ) , 0 , 0 , 0 , 0 } ,
	{ "vz" , PLY_FLOAT , PLY_FLOAT , (int)offsetof( PlyVFFace , v.coords[2] ) , 0 , 0 , 0 , 0 } ,
};
template<>
PlyProperty PlyVFFace< double >::ReadProperties[] =
{
	{ "vertex_indices" , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , nr_vertices ) } ,
	{ "vx" , PLY_DOUBLE , PLY_DOUBLE , (int)offsetof( PlyVFFace , v.coords[0] ) , 0 , 0 , 0 , 0 } ,
	{ "vy" , PLY_DOUBLE , PLY_DOUBLE , (int)offsetof( PlyVFFace , v.coords[1] ) , 0 , 0 , 0 , 0 } ,
	{ "vz" , PLY_DOUBLE , PLY_DOUBLE , (int)offsetof( PlyVFFace , v.coords[2] ) , 0 , 0 , 0 , 0 } ,
};
template<>
PlyProperty PlyVFFace< float >::WriteProperties[] =
{
	{ "vertex_indices" , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , nr_vertices ) } ,
	{ "vx" , PLY_FLOAT , PLY_FLOAT , (int)offsetof( PlyVFFace , v.coords[0] ) , 0 , 0 , 0 , 0 } ,
	{ "vy" , PLY_FLOAT , PLY_FLOAT , (int)offsetof( PlyVFFace , v.coords[1] ) , 0 , 0 , 0 , 0 } ,
	{ "vz" , PLY_FLOAT , PLY_FLOAT , (int)offsetof( PlyVFFace , v.coords[2] ) , 0 , 0 , 0 , 0 } ,
};
template<>
PlyProperty PlyVFFace< double >::WriteProperties[] =
{
	{ "vertex_indices" , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyVFFace , nr_vertices ) } ,
	{ "vx" , PLY_DOUBLE , PLY_DOUBLE , (int)offsetof( PlyVFFace , v.coords[0] ) , 0 , 0 , 0 , 0 } ,
	{ "vy" , PLY_DOUBLE , PLY_DOUBLE , (int)offsetof( PlyVFFace , v.coords[1] ) , 0 , 0 , 0 , 0 } ,
	{ "vz" , PLY_DOUBLE , PLY_DOUBLE , (int)offsetof( PlyVFFace , v.coords[2] ) , 0 , 0 , 0 , 0 } ,
};

template< class Real >
struct PlyTexturedFace
{
	unsigned int nr_vertices , nr_uv_coordinates;
	int *vertices;
	Real *uv_coordinates;
	PlyTexturedFace( void ){ vertices = NULL , uv_coordinates = NULL , nr_vertices = nr_uv_coordinates = 0; }
	~PlyTexturedFace( void ){ resize(0); }
	PlyTexturedFace( const PlyTexturedFace& face )
	{
		vertices = NULL , uv_coordinates = NULL;
		(*this) = face;
	}
	PlyTexturedFace& operator = ( const PlyTexturedFace& face )
	{
		if( vertices ) free( vertices ) , vertices = NULL;
		if( uv_coordinates ) free( uv_coordinates ) , uv_coordinates = NULL;
		nr_vertices = face.nr_vertices , nr_uv_coordinates = face.nr_uv_coordinates;
		if( nr_vertices ) vertices = (int*)malloc( sizeof(int)*nr_vertices );
		else              vertices = NULL;
		if( nr_uv_coordinates ) uv_coordinates = (Real*)malloc( sizeof(Real)*nr_uv_coordinates );
		else                    uv_coordinates = NULL;
		memcpy( vertices , face.vertices , sizeof(int)*nr_vertices );
		memcpy( uv_coordinates , face.uv_coordinates , sizeof(Real)*nr_uv_coordinates );
		return *this;
	}
	void resize( unsigned int count )
	{
		if( vertices ) free( vertices ) , vertices = NULL;
		if( uv_coordinates ) free( uv_coordinates ) , uv_coordinates = NULL;
		nr_vertices = nr_uv_coordinates = 0;
		if( count )
		{
			vertices = (int*)malloc( sizeof(int)*count ) , nr_vertices = count;
			uv_coordinates = (Real*)malloc( sizeof(Real)*2*count ) , nr_uv_coordinates = count*2;
		}
	}
	int& operator[] ( int idx ){ return vertices[idx]; }
	const int& operator[] ( int idx ) const { return vertices[idx]; }
#if 1
	Point2D< Real >  texture( int idx ) const { return Point2D< Real >( uv_coordinates[2*idx] , uv_coordinates[2*idx+1] ); }
	Point2D< Real >& texture( int idx )       { return *( (Point2D< Real >*)(uv_coordinates+2*idx) ); }
#else
	Point2D< Real > texture( int idx ) const { return Point2D< Real >( uv_coordinates[2*idx] , uv_coordinates[2*idx+1] ); }
#endif
	int size( void ) const { return nr_vertices; }

	const static int ReadComponents=2;
	const static int WriteComponents=2;
	static PlyProperty ReadProperties[];
	static PlyProperty WriteProperties[];
};
template<>
PlyProperty PlyTexturedFace< float >::ReadProperties[] =
{
	{ "vertex_indices" , PLY_INT , PLY_INT , offsetof( PlyTexturedFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_vertices ) } ,
	{ "texcoord" , PLY_FLOAT , PLY_FLOAT , (int)offsetof( PlyTexturedFace , uv_coordinates ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_uv_coordinates ) } ,
};
template<>
PlyProperty PlyTexturedFace< double >::ReadProperties[] =
{
	{ "vertex_indices" , PLY_INT , PLY_INT , offsetof( PlyTexturedFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_vertices ) } ,
	{ "texcoord" , PLY_DOUBLE , PLY_DOUBLE , (int)offsetof( PlyTexturedFace , uv_coordinates ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_uv_coordinates ) } ,
};
template<>
PlyProperty PlyTexturedFace< float >::WriteProperties[] =
{
	{ "vertex_indices" , PLY_INT , PLY_INT , offsetof( PlyTexturedFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_vertices ) } ,
	{ "texcoord" , PLY_FLOAT , PLY_FLOAT , (int)offsetof( PlyTexturedFace , uv_coordinates ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_uv_coordinates ) } ,
};
template<>
PlyProperty PlyTexturedFace< double >::WriteProperties[] =
{
	{ "vertex_indices" , PLY_INT , PLY_INT , offsetof( PlyTexturedFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_vertices ) } ,
	{ "texcoord" , PLY_DOUBLE , PLY_DOUBLE , (int)offsetof( PlyTexturedFace , uv_coordinates ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_uv_coordinates ) } ,
};



int PlyDefaultFileType( void ){ return PLY_ASCII; }

template< class Vertex >
int PlyRead
	(
	const char* fileName ,
	std::vector< Vertex >& vertices ,
	std::vector< std::pair< int , int > >* edges ,
	std::vector< std::vector< int > >* polygons ,
	PlyProperty* properties , bool* propertyFlags , int propertyNum ,
	int& file_type ,
	char*** comments=NULL , int* commentNum=NULL
	);
template<class Vertex>
int PlyWrite
	(
	const char* fileName ,
	const std::vector< Vertex >& vertices , 
	const std::vector< std::pair< int , int > >* edges ,
	const std::vector< std::vector< int > >* polygons ,
	PlyProperty* properties , int propertyNum,
	int file_type ,
	char** comments=NULL , const int& commentNum=0
	);

template<class Vertex>
int PlyReadPoints( const char* fileName,
				  std::vector<Vertex>& vertices,
				  PlyProperty* properties , bool* propertyFlags , int propertyNum,
				  int& file_type,
				  char*** comments=NULL,int* commentNum=NULL);
template<class Vertex>
int PlyWritePoints( const char* fileName,
				   const std::vector<Vertex>& vertices,
				   PlyProperty* properties,int propertyNum,
				   int file_type,
				   char** comments=NULL,const int& commentNum=0);

template<class Vertex>
int PlyReadPolygons( const char* fileName,
					std::vector< Vertex >& vertices,std::vector<std::vector<int> >& polygons,
					PlyProperty* properties , bool* propertyFlags , int propertyNum,
					int& file_type,
					char*** comments=NULL,int* commentNum=NULL);
template< class Vertex , class Polygon >
int PlyReadPolygons( const char* fileName,
					std::vector< Vertex >& vertices , std::vector< Polygon >& polygons ,
					PlyProperty*  vertexProperties , bool*  vertexPropertyFlags , int  vertexPropertyNum ,
					PlyProperty* polygonProperties , bool* polygonPropertyFlags , int polygonPropertyNum ,
					int& file_type ,
					char*** comments=NULL , int* commentNum=NULL );

template< class Vertex >
int PlyReadTriangles( const char* fileName ,
					  std::vector<Vertex>& vertices , std::vector< TriangleIndex >& triangles ,
					  PlyProperty* properties , bool* propertyFlags , int propertyNum ,
					  int& file_type ,
					  char*** comments=NULL , int* commentNum=NULL );
template<class Vertex>
int PlyWritePolygons( const char* fileName,
					 const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,
					 PlyProperty* properties,int propertyNum,
					 int file_type,
					 char** comments=NULL,const int& commentNum=0);
template< class Vertex , class Polygon >
int PlyWritePolygons( const char* fileName,
					 const std::vector< Vertex >& vertices,const std::vector< Polygon >& polygons,
					 PlyProperty*  vertexProperties , int  vertexPropertyNum ,
					 PlyProperty* polygonProperties , int polygonPropertyNum ,
					 int file_type,
					 char** comments=NULL,const int& commentNum=0);
template< class Vertex >
int PlyWriteTriangles( const char* fileName ,
					   const std::vector<Vertex>& vertices , const std::vector< TriangleIndex >& triangles ,
					   PlyProperty* properties , int propertyNum ,
					   int file_type ,
					   char** comments=NULL , const int& commentNum=0);

template< class Vertex >
int PlyWriteColorTriangles( const char* fileName ,
					   const std::vector< Vertex >& vertices , const std::vector< std::pair< TriangleIndex , Point3D< float > > >& triangles ,
					   PlyProperty* properties , int propertyNum ,
					   int file_type ,
					   char** comments=NULL , const int& commentNum=0);

template<class Vertex>
int PlyReadGrid( const char* fileName,
				std::vector<Vertex>& vertices,Grid<int>& grid,
				PlyProperty* properties , bool* propertyFlags , int propertyNum,
				int& file_type,
				char*** comments=NULL,int* commentNum=NULL);


template< class Vertex > int  MReadTriangles( const char* fileName ,       std::vector< Vertex >& vertices ,       std::vector< TriangleIndex >& triangles );
template< class Vertex > int MWriteTriangles( const char* fileName , const std::vector< Vertex >& vertices , const std::vector< TriangleIndex >& triangles );
#include "Ply.inl"

#endif // PLY_INCLUDED
