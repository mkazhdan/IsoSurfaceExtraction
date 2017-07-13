/*
Copyright (c) 2015, Michael Kazhdan
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

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unordered_map>
#include <algorithm>
#include <sys/timeb.h>
#ifndef WIN32
#include <sys/time.h>
#endif // WIN32

#include <Src/CmdLineParser.h>
#include <Src/Geometry.h>
#include <Src/Ply.h>
#include <Src/MarchingCubes.h>
#include <Src/Array.h>


const float DEFAULT_DIMENSIONS[] = { 1.f , 1.f , 1.f };
cmdLineParameter< char* > In( "in" ) , Out( "out" );
cmdLineParameterArray< int , 3 > Resolution( "res" );
cmdLineParameter< int > SmoothIterations( "sIters" , 0 );
cmdLineParameter< float > IsoValue( "iso" , 0.f );
cmdLineParameterArray< float , 3 > Dimensions( "dim" , DEFAULT_DIMENSIONS );
cmdLineReadable Float( "float" ) , FullCaseTable( "full" ) , FlipOrientation( "flip" ) , QuadraticFit( "quadratic" ) , Polygons( "polygons" ) , NonManifold( "nonManifold" ) , FlipBytes( "flipBytes" );

cmdLineReadable* params[] = { &In , &Out , &Resolution , &IsoValue , &FullCaseTable , &FlipOrientation , &QuadraticFit , &Polygons , &SmoothIterations , &Float , &Dimensions , &NonManifold , &FlipBytes , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input voxel grid>\n" , In.name );
	printf( "\t --%s <input resolution>\n" , Resolution.name );
	printf( "\t[--%s <output iso-surface>]\n" , Out.name );
	printf( "\t[--%s <iso-value>=%f]\n" , IsoValue.name , IsoValue.value );
	printf( "\t[--%s <smoothing iterations>=%d]\n" , SmoothIterations.name , SmoothIterations.value );
	printf( "\t[--%s <dimensions of a voxel>=%f %f %f]\n" , Dimensions.name , Dimensions.values[0] , Dimensions.values[1] , Dimensions.values[2] );
	printf( "\t[--%s]\n" , FullCaseTable.name );
	printf( "\t[--%s]\n" , FlipOrientation.name );
	printf( "\t[--%s]\n" , QuadraticFit.name );
	printf( "\t[--%s]\n" , Polygons.name );
	printf( "\t[--%s]\n" , NonManifold.name );
	printf( "\t[--%s]\n" , Float.name );
	printf( "\t[--%s]\n" , FlipBytes.name );
}

inline double Time( void )
{
#ifdef WIN32
	struct _timeb t;
	_ftime( &t );
	return double( t.time ) + double( t.millitm ) / 1000.0;
#else // WIN32
	struct timeval t;
	gettimeofday( &t , NULL );
	return t.tv_sec + double( t.tv_usec ) / 1000000;
#endif // WIN32
}


float    LinearInterpolant( float x1 , float x2 , float isoValue ){ return ( isoValue-x1 ) / ( x2-x1 ); }
float QuadraticInterpolant( float x0 , float x1 , float x2 , float x3 , float isoValue )
{
	// Adjust so that we are looking for a zero-crossing
	x0 -= isoValue , x1 -= isoValue , x2 -= isoValue , x3 -= isoValue;
	// Estimate the derivatives at x1 and x2
	float dx1 = (x2-x0) / 2.f , dx2 = (x3-x1) / 2.f;
	// Solve for the quadratic polynomial:
	//		P(x) = a x^2 + b x + c 
	// such that:
	//		P(0) = x1 , P(1) = x2 , and minimizing || P'(0) - dx1 ||^2 + || P'(1) - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || 2*a + b - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || 2*x2 - 2*x1 - b - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || b - ( 2*x2 - 2*x1 - dx2 ) ||^2
	//	=>  c = x1 , b = ( 2*x2 - 2*x1 - dx2 + dx1 ) / 2 , a = x2 - x1 - b
	//	=>  c = x1 , b = ( x2 - x1 ) - ( dx2 - dx1 ) / 2 , a = ( dx2 - dx1 ) / 2

	double a = (dx2-dx1)/2.f , b = (dx1-dx2)/2.f + x2 - x1 , c = x1;
	if( !a )
	{
		// Solve b * x + c = 0
		return (float)( -c / b );
	}
	else
	{
		// Solve a x^2 + b x + c = 0
		b /= a , c /= a;
		double disc = b*b - 4.*c;
		if( disc<0 ) fprintf( stderr , "[ERROR] Negative discriminant: %g\n" , disc ) , exit( 0 );
		disc = sqrt( disc );
		double r1 = ( - b - disc ) / 2. , r2 = ( - b + disc ) / 2.;
		if( r2<0 || r1>1 ) fprintf( stderr , "[ERROR] Roots out of bounds: %g %g\n" , r1 , r2 ) , exit( 0 );
		if( r2>1 ) return (float)r1;
		else       return (float)r2;
	}
}

struct IsoVertex
{
	int dir , idx[3];
	Point3D< float > p;
	IsoVertex( Point3D< float > p , int dir , int x , int y , int z ){ this->p = p , this->dir = dir , idx[0] = x , idx[1] = y , idx[2] = z; }
#define _ABS_( a ) ( (a)<0 ? -(a) : (a) )
	static bool CoFacial( const IsoVertex& t1 , const IsoVertex& t2 )
	{
		int d[] = { _ABS_( t1.idx[0] - t2.idx[0] ) , _ABS_( t1.idx[1] - t2.idx[1] ) , _ABS_( t1.idx[2] - t2.idx[2] ) };
		if( t1.dir==t2.dir ) return d[t1.dir]==0 && ( ( d[(t1.dir+1)%3]==0 && d[(t1.dir+2)%3]<=1 ) || ( d[(t1.dir+2)%3]==0 && d[(t1.dir+1)%3]<=1 ) );
		else                 return d[ 3 - t1.dir - t2.dir ]==0 && d[t1.dir]<=1 && d[t2.dir]<=1;
	}
#undef _ABS_
};
void SetFlags( int resX , int resY , ConstPointer( float ) values , float isoValue , Pointer( unsigned char ) flags )
{
#pragma omp parallel for
	for( int i=0 ; i<resX*resY ; i++ ) flags[i] = MarchingCubes::ValueLabel( values[i] , isoValue );
}
void SetZVertices
(
	int resX , int resY , int z , 
	ConstPointer( float ) values0 , ConstPointer( float ) values1 , ConstPointer( float ) values2 , ConstPointer( float ) values3 , 
	ConstPointer( unsigned char ) flags1 , ConstPointer( unsigned char ) flags2 , 
	float isoValue , bool quadratic ,
	std::unordered_map< long long , int >& isoVertexMap ,
	std::vector< IsoVertex >& vertices
)
{
#define INDEX( x , y ) ( x + (y)*resX )
#pragma omp parallel for
	for( int i=0 ; i<resX ; i++ ) for( int j=0 ; j<resY ; j++ )
	{
		int idx = INDEX( i , j );
		if( flags1[idx]!=flags2[idx] )
		{
			float iso;
			if( quadratic ) iso = QuadraticInterpolant( values0 ? values0[idx] : values1[idx] , values1[idx] , values2[idx] , values3 ? values3[idx] : values2[idx] , isoValue );
			else iso = LinearInterpolant( values1[idx] , values2[idx] , isoValue );
			Point3D< float > p = Point3D< float >( (float)i , (float)j , (float)z + iso );
			long long key = i + j*(resX);
#pragma omp critical
			{
				isoVertexMap[key] = (int)vertices.size();
				vertices.push_back( IsoVertex( p , 2 , i , j , z ) );
			}
		}
	}
#undef INDEX
}
void SetXYVertices
(
	int resX , int resY , int z , 
	ConstPointer( float ) values ,
	ConstPointer( unsigned char ) flags ,
	float isoValue , bool quadratic ,
	std::unordered_map< long long , int >& xIsoVertexMap , std::unordered_map< long long , int >& yIsoVertexMap ,
	std::vector< IsoVertex >& vertices
)
{
#define INDEX( x , y ) ( x + (y)*resX )
#pragma omp parallel for
	for( int i=0 ; i<resX-1 ; i++ ) for( int j=0 ; j<resY ; j++ )
	{
		int idx1 = INDEX( i , j ) , idx2 = INDEX( i+1 , j );
		if( flags[idx1]!=flags[idx2] )
		{
			float iso;
			if( quadratic ) iso = QuadraticInterpolant( i>0 ? values[ INDEX(i-1,j) ] : values[idx1] , values[idx1] , values[idx2] , i+1<resX-1 ? values[ INDEX(i+2,j) ] : values[idx2] , isoValue );
			else iso = LinearInterpolant( values[idx1] , values[idx2] , isoValue );
			Point3D< float > p = Point3D< float >( (float)i + iso , (float)j , (float)z );
			long long key = i + j*(resX);
#pragma omp critical
			{
				xIsoVertexMap[key] = (int)vertices.size();
				vertices.push_back( IsoVertex( p , 0 , i , j , z ) );
			}
		}
	}
#pragma omp parallel for
	for( int i=0 ; i<resX ; i++ ) for( int j=0 ; j<resY-1 ; j++ )
	{
		int idx1 = INDEX( i , j ) , idx2 = INDEX( i , j+1 );
		if( flags[idx1]!=flags[idx2] )
		{
			float iso;
			if( quadratic ) iso = QuadraticInterpolant( j>0 ? values[ INDEX(i,j-1) ] : values[idx1] , values[idx1] , values[idx2] , j+1<resY-1 ? values[ INDEX(i,j+2) ] : values[idx2] , isoValue );
			else iso = LinearInterpolant( values[idx1] , values[idx2] , isoValue );
			Point3D< float > p = Point3D< float >( (float)i , (float)j + iso , (float)z );
			long long key = i + j*(resX);
#pragma omp critical
			{
				yIsoVertexMap[key] = (int)vertices.size();
				vertices.push_back( IsoVertex( p , 1 , i , j , z ) );
			}
		}
	}
#undef INDEX
}
void SetPolygons
(
	int resX , int resY , int z , 
	ConstPointer( float ) values1 , ConstPointer( float ) values2 ,
	float isoValue , bool fullCaseTable , bool flip ,
	const std::unordered_map< long long , int >& xIsoVertexMap1 , const std::unordered_map< long long , int >& xIsoVertexMap2 ,
	const std::unordered_map< long long , int >& yIsoVertexMap1 , const std::unordered_map< long long , int >& yIsoVertexMap2 ,
	const std::unordered_map< long long , int >& zIsoVertexMap ,
	const std::vector< IsoVertex >& vertices , std::vector< std::vector< int > >& polygons
)
{
#define INDEX( x , y ) ( x + (y)*resX )
#pragma omp parallel for
	for( int i=0 ; i<resX-1 ; i++ ) for( int j=0 ; j<resY-1 ; j++ )
	{
		float _values[Cube::CORNERS];
		for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ )
		{
			_values[ Cube::CornerIndex(cx,cy,0) ] = values1[ INDEX(i+cx,j+cy) ];
			_values[ Cube::CornerIndex(cx,cy,1) ] = values2[ INDEX(i+cx,j+cy) ];
		}
		int mcIndex = fullCaseTable ? MarchingCubes::GetFullIndex( _values , isoValue ) : MarchingCubes::GetIndex( _values , isoValue );
		const std::vector< std::vector< int > >& isoPolygons = MarchingCubes::caseTable( mcIndex , fullCaseTable );
		for( int p=0 ; p<isoPolygons.size() ; p++ )
		{
			const std::vector< int >& isoPolygon = isoPolygons[p];
			std::vector< int > polygon( isoPolygon.size() );
			for( int v=0 ; v<isoPolygon.size() ; v++ )
			{
				int orientation , i1 , i2;
				Cube::FactorEdgeIndex( isoPolygon[v] , orientation , i1 , i2 );
				long long key;
				std::unordered_map< long long , int >::const_iterator iter;
				bool success;
				switch( orientation )
				{
				case 0:
					key = (i   ) + (j+i1)*resX;
					if( i2==0 ){ iter = xIsoVertexMap1.find( key ) ; success = iter!=xIsoVertexMap1.end(); }
					else       { iter = xIsoVertexMap2.find( key ) ; success = iter!=xIsoVertexMap2.end(); }
					break;
				case 1:
					key = (i+i1) + (j   )*resX;
					if( i2==0 ){ iter = yIsoVertexMap1.find( key ) ; success = iter!=yIsoVertexMap1.end(); }
					else       { iter = yIsoVertexMap2.find( key ) ; success = iter!=yIsoVertexMap2.end(); }
					break;
				case 2:
					key = (i+i1) + (j+i2)*resX;
					iter = zIsoVertexMap.find( key ) ; success = iter!=zIsoVertexMap.end();
					break;
				}

				if( !success )
				{
					fprintf( stderr , "[ERROR] Couldn't find iso-vertex in map:\n" );
					printf( "\t%d: " , orientation );
					switch( orientation )
					{
					case 0: printf( "%d %d %d\n" , i , j+i1 , z+i2 ) ; break;
					case 1: printf( "%d %d %d\n" , i+i1 , j , z+i2 ) ; break;
					case 2: printf( "%d %d %d\n" , i+i1 , j+i2 , z ) ; break;
					}
					exit( 0 );
				}
				if( flip ) polygon[polygon.size()-1-v] = iter->second;
				else       polygon[v] = iter->second;
			}
#pragma omp critical
			polygons.push_back( polygon );
		}
	}
#undef INDEX
}

void ExtractIsoSurface( int resX , int resY , int resZ , ConstPointer( float ) values , float isoValue , std::vector< IsoVertex >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , bool quadratic , bool flip )
{
	std::unordered_map< long long , int > xIsoVertexMap[2] , yIsoVertexMap[2] , zIsoVertexMap;
	Pointer( unsigned char ) flags[2];
	flags[0] = NewPointer< unsigned char >( resX*resY );
	flags[1] = NewPointer< unsigned char >( resX*resY );

	if( fullCaseTable ) MarchingCubes::SetFullCaseTable();
	else                MarchingCubes::SetCaseTable();

	SetFlags     ( resX , resY ,     values , isoValue , flags[0] );
	SetXYVertices( resX , resY , 0 , values ,            flags[0] , isoValue , quadratic , xIsoVertexMap[0] , yIsoVertexMap[0] , vertices );
	for( int z=0 ; z<resZ-1 ; z++ )
	{
		int z0 = z&1 , z1 = (z+1)&1;
		xIsoVertexMap[z1].clear() , yIsoVertexMap[z1].clear() , zIsoVertexMap.clear();
		SetFlags     ( resX , resY ,       values + (z+1)*resX*resY , isoValue , flags[z1] );
		SetXYVertices( resX , resY , z+1 , values + (z+1)*resX*resY ,            flags[z1] , isoValue , quadratic , xIsoVertexMap[z1] , yIsoVertexMap[z1] , vertices );
		SetZVertices ( resX , resY , z , z>0 ? values + (z-1)*resX*resY : NullPointer< float >() , values + z*resX*resY , values + (z+1)*resX*resY , z+1<resZ-1 ? values + (z+2)*resX*resY : NullPointer< float >() , flags[z0] , flags[z1] , isoValue , quadratic , zIsoVertexMap , vertices );
		SetPolygons  ( resX , resY , z , values + z*resX*resY , values+(z+1)*resX*resY , isoValue , fullCaseTable , flip , xIsoVertexMap[z0] , xIsoVertexMap[z1] , yIsoVertexMap[z0] , yIsoVertexMap[z1] , zIsoVertexMap ,  vertices , polygons );
	}
	DeletePointer( flags[0] );
	DeletePointer( flags[1] );
}
int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }

	FILE* fp = fopen( In.value , "rb" );
	if( !fp )
	{
		fprintf( stderr , "[ERROR] Failed to open file for reading: %s\n" , In.value );
		return EXIT_FAILURE;
	}
	if( !Resolution.set )
	{
		int res;
		if( fread( &res , sizeof(int) , 1 ,fp )!=1 )
		{
			fprintf( stderr , "[ERROR] Failed to read voxel grid resolution from file.\n" );
			return EXIT_FAILURE;
		}
		Resolution.values[0] = Resolution.values[1] = Resolution.values[2] = res;
	}

	Pointer( float ) voxelValues = NewPointer< float >( Resolution.values[0] * Resolution.values[1] * Resolution.values[2] );
	if( !voxelValues )
	{
		fprintf( stderr , "[ERROR] Failed to allocte voxel grid: %d x %d x %d\n" , Resolution.values[0] , Resolution.values[1] , Resolution.values[2] );
		fclose( fp );
		return EXIT_FAILURE;
	}

	if( Float.set )
	{
		if( fread( voxelValues , sizeof( float ) , Resolution.values[0] * Resolution.values[1] * Resolution.values[2] , fp )!=Resolution.values[0]*Resolution.values[1]*Resolution.values[2] )
		{
			fprintf( stderr , "[ERROR] Failed to read voxel grid from file.\n" );
			return EXIT_FAILURE;
		}
		if( FlipBytes.set )
		{
			Pointer( unsigned char ) _voxelValues = ( Pointer( unsigned char ) )voxelValues;
			for( int i=0 ; i<Resolution.values[0] * Resolution.values[1] * Resolution.values[2] ; i++ ) for( int j=0 ; j<sizeof(float)/2 ; j++ )
			{
				unsigned char foo = _voxelValues[ i*sizeof(float) + j ];
				_voxelValues[ i*sizeof(float) + j ] = _voxelValues[ i*sizeof(float) + sizeof(float) - 1 - j ];
				_voxelValues[ i*sizeof(float) + sizeof(float) - 1 - j ] = foo;
			}
		}
	}
	else
	{
		Pointer( unsigned char ) _voxelValues = NewPointer< unsigned char >( Resolution.values[0] * Resolution.values[1] * Resolution.values[2] );
		if( fread( _voxelValues , sizeof( unsigned char ) , Resolution.values[0] * Resolution.values[1] * Resolution.values[2] , fp )!=Resolution.values[0]*Resolution.values[1]*Resolution.values[2] )
		{
			fprintf( stderr , "[ERROR] Failed to read voxel grid from file.\n" );
			return EXIT_FAILURE;
		}
#pragma omp parallel for
		for( int i=0 ; i<Resolution.values[0]*Resolution.values[1]*Resolution.values[2] ; i++ ) voxelValues[i] = (float)_voxelValues[i];
		DeletePointer( _voxelValues );
	}
	fclose( fp );


#define INDEX( x , y , z ) ( (x) + (y)*Resolution.values[0] + (z)*Resolution.values[0]*Resolution.values[1] )

	float min , max;
	min = max = voxelValues[0];
	for( int x=0 ; x<Resolution.values[0] ; x++ ) for( int y=0 ; y<Resolution.values[1] ; y++ ) for( int z=0 ; z<Resolution.values[2] ; z++ )
		min = std::min< float >( min , voxelValues[ INDEX(x,y,z) ] ) , max = std::max< float >( max , voxelValues[ INDEX(x,y,z) ] );
	printf( "Value Range: [%f,%f]\n" , min , max );

	if( SmoothIterations.value>0 )
	{
		float stencil[] = { 0.5f , 1.f , 0.5f };
		Pointer( float ) _voxelValues = NewPointer< float >( Resolution.values[0] * Resolution.values[1] * Resolution.values[2] );
		for( int i=0 ; i<SmoothIterations.value ; i++ )
		{
#pragma omp parallel for
			for( int x=0 ; x<Resolution.values[0] ; x++ ) for( int y=0 ; y<Resolution.values[1] ; y++ ) for( int z=0 ; z<Resolution.values[2] ; z++ )
			{
				float weightSum = 0.f;
				_voxelValues[ INDEX(x,y,z) ] = 0.f;
				for( int xx=-1 ; xx<=1 ; xx++ ) for( int yy=-1 ; yy<=1 ; yy++ ) for( int zz=-1; zz<=1 ; zz++ )
					if( x+xx>=0 && x+xx<Resolution.values[0] && y+yy>=0 && y+yy<Resolution.values[1] && z+zz>=0 && z+zz<Resolution.values[2] )
					{
						_voxelValues[ INDEX(x,y,z) ] += voxelValues[ INDEX(x+xx,y+yy,z+zz) ] * stencil[xx+1] * stencil[yy+1] * stencil[zz+1];
						weightSum += stencil[xx+1] * stencil[yy+1] * stencil[zz+1];
					}
				_voxelValues[ INDEX(x,y,z) ] /= weightSum;
			}
			memcpy( voxelValues , _voxelValues , sizeof(float) * Resolution.values[0] * Resolution.values[1] * Resolution.values[2] );
		}
		DeletePointer( _voxelValues );
	}
#undef INDEX

	std::vector< IsoVertex > vertices;
	std::vector< std::vector< int > > polygons;
	double t = Time();
	ExtractIsoSurface( Resolution.values[0] , Resolution.values[1] , Resolution.values[2] , voxelValues , IsoValue.value , vertices , polygons , FullCaseTable.set , QuadraticFit.set , FlipOrientation.set );
	printf( "Got iso-surface: %.2f(s)\n" , Time()-t );
	DeletePointer( voxelValues );

	if( Out.set )
	{
		std::vector< PlyVertex< float > > _vertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) for( int d=0 ; d<3 ; d++ ) _vertices[i].point[d] = vertices[i].p[d] * Dimensions.values[d];
		if( Polygons.set )
		{
			PlyWritePolygons( Out.value , _vertices , polygons , PlyVertex< float >::WriteProperties , PlyVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
			printf( "Vertices / Polygons: %d / %d\n" , (int)vertices.size() , (int)polygons.size() );
		}
		else
		{
			std::vector< TriangleIndex > triangles;
			MinimalAreaTriangulation< float > mat;

			for( int i=0 ; i<polygons.size() ; i++ )
			{
				// To ensure that we have no more than two triangles adjacent on an edge,
				// we avoid creating a minimial area triangulation when it could introduce a new
				// edge that is on a face of a cube
				bool isCofacial = false;
				if( !NonManifold.set )
					for( int j=0 ; j<(int)polygons[i].size() ; j++ ) for( int k=0 ; k<j ; k++ )
						if( (j+1)%polygons[i].size()!=k && (k+1)%polygons[i].size()!=j )
							if( IsoVertex::CoFacial( vertices[ polygons[i][j] ] , vertices[ polygons[i][k] ] ) ) isCofacial = true;
				if( isCofacial )
				{
					TriangleIndex triangle;
					PlyVertex< float > plyVertex;
					for( int j=0 ; j<(int)polygons[i].size() ; j++ ) plyVertex.point += vertices[ polygons[i][j] ].p;
					plyVertex.point /= (float)polygons[i].size();
					int cIdx = (int)_vertices.size();
					for( int d=0 ; d<3 ; d++ ) plyVertex.point[d] *= Dimensions.values[d];
					_vertices.push_back( plyVertex );
					for( int j=0 ; j<(int)polygons[i].size() ; j++ )
					{
						triangle[0] = polygons[i][j];
						triangle[1] = polygons[i][(j+1)%polygons[i].size()];
						triangle[2] = cIdx;
						triangles.push_back( triangle );
					}
				}
				else
				{
					std::vector< Point3D< float > > _polygon( polygons[i].size() );
					std::vector< TriangleIndex > _triangles;
					for( int j=0 ; j<polygons[i].size() ; j++ ) _polygon[j] = vertices[ polygons[i][j] ].p;
					mat.GetTriangulation( _polygon , _triangles );
					for( int j=0 ; j<_triangles.size() ; j++ )
					{
						TriangleIndex tri;
						for( int k=0 ; k<3 ; k++ ) tri[k] = polygons[i][ _triangles[j][k] ];
						triangles.push_back( tri );
					}
				}
			}

			PlyWriteTriangles( Out.value , _vertices , triangles , PlyVertex< float >::WriteProperties , PlyVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
			printf( "Vertices / Triangles: %d / %d\n" , (int)_vertices.size() , (int)triangles.size() );
		}
	}

	return EXIT_SUCCESS;
}