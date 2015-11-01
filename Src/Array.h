/*
Copyright (c) 2011, Michael Kazhdan and Ming Chuang
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

#ifndef ARRAY_INCLUDED
#define ARRAY_INCLUDED

#include <vector>

#ifndef ARRAY_DEBUG
	#define ARRAY_DEBUG 0
#endif // ARRAY_DEBUG

#ifdef _WIN64
#define ASSERT( x ) { if( !( x ) ) __debugbreak(); }
#else // !_WIN64
#ifdef _WIN32
#define ASSERT( x ) { if( !( x ) ) _asm{ int 0x03 } }
#else // !_WIN32
#define ASSERT( x ) { if( !( x ) ) exit(0); }
#endif // _WIN32
#endif // _WIN64

// Code from http://stackoverflow.com
void* aligned_malloc( size_t size , size_t align )
{
	// Align enough for the data, the alignment padding, and room to store a pointer to the actual start of the memory
	void*  mem = malloc( size + align + sizeof( void* ) );
	// The position at which we could potentially start addressing
	char* amem = ( (char*)mem ) + sizeof( void* );
	// Add align-1 to the start of the address and then zero out at most of the first align-1 bits.
	amem = ( char* )( ( (size_t)( ( (char*)amem ) + (align-1) ) ) & ~( align-1 ) );
	// Pre-write the actual address
	( ( void** ) amem )[-1] = mem;
	return amem;
}
void aligned_free( void* mem ) { free( ( ( void** )mem )[-1] ); }

#if ARRAY_DEBUG
#pragma message ( "[WARNING] Array debugging is enabled" )
#include "Array.inl"
#define      Pointer( ... )      Array< __VA_ARGS__ >
#define ConstPointer( ... ) ConstArray< __VA_ARGS__ >
template< class C > void        FreePointer( Array< C >& a ){ a.Free( ); }
template< class C > void AlignedFreePointer( Array< C >& a ){ a.Free( ); }
template< class C > void       VFreePointer( Array< C >& a ){ a.Free( ); }
template< class C > void      DeletePointer( Array< C >& a ){ a.Delete( ); }

template< class C > Array< C >          NewPointer(                  size_t size ,                    const char* name=NULL ){ return Array< C >::New         (     size ,                     name ); }
template< class C > Array< C >        AllocPointer(                  size_t size ,                    const char* name=NULL ){ return Array< C >::Alloc       (     size ,             false , name ); }
template< class C > Array< C > AlignedAllocPointer(                  size_t size , size_t alignment , const char* name=NULL ){ return Array< C >::AlignedAlloc(     size , alignment , false , name ); }
template< class C > Array< C >      ReAllocPointer( Array< C >&  a , size_t size ,                    const char* name=NULL ){ return Array< C >::ReAlloc     ( a , size ,             false , name ); }

template< class C > Array< C > NullPointer( void ){ return Array< C >( ); }

template< class C >      Array< C > GetPointer(       C& c ) { return      Array< C >::FromPointer( &c , 1 ); }
template< class C > ConstArray< C > GetPointer( const C& c ) { return ConstArray< C >::FromPointer( &c , 1 ); }
template< class C >      Array< C > GetPointer(       std::vector< C >& v ){ return      Array< C >::FromPointer( &v[0] , v.size() ); }
template< class C > ConstArray< C > GetPointer( const std::vector< C >& v ){ return ConstArray< C >::FromPointer( &v[0] , v.size() ); }
//template< class C >      Array< C > GetPointer(      Array< C >& a ) { return a; }
//template< class C > ConstArray< C > GetPointer( ConstArray< C >& a ) { return a; }

template< class C >       C* GetAddress(       Array< C >&  a ){ return a.pointer(); }
template< class C > const C* GetAddress( const Array< C >&  a ){ return a.pointer(); }
template< class C > const C* GetAddress(  ConstArray< C >&  a ){ return a.pointer(); }

#else // !ARRAY_DEBUG
#define      Pointer( ... )       __VA_ARGS__*
#define ConstPointer( ... ) const __VA_ARGS__*

#define        FreePointer( ... ) { if( __VA_ARGS__ )         free( __VA_ARGS__ ) ,                   __VA_ARGS__ = NULL; }
#define AlignedFreePointer( ... ) { if( __VA_ARGS__ ) aligned_free( __VA_ARGS__ ) ,                   __VA_ARGS__ = NULL; }
#define      DeletePointer( ... ) { if( __VA_ARGS__ )      delete[] __VA_ARGS__ ,                     __VA_ARGS__ = NULL; }

template< class C > C*          NewPointer(        size_t size ,                    const char* name=NULL ){ return new C[size]; }
template< class C > C*        AllocPointer(        size_t size ,                    const char* name=NULL ){ return (C*)        malloc(        sizeof(C) * size             ); }
template< class C > C* AlignedAllocPointer(        size_t size , size_t alignment , const char* name=NULL ){ return (C*)aligned_malloc(        sizeof(C) * size , alignment ); }
template< class C > C*      ReAllocPointer( C* c , size_t size ,                    const char* name=NULL ){ return (C*)       realloc( c    , sizeof(C) * size             ); }

template< class C > C* NullPointer( void ){ return NULL; }

template< class C >       C* GetPointer(       C& c ){ return &c; }
template< class C > const C* GetPointer( const C& c ){ return &c; }
template< class C >       C* GetPointer(       std::vector< C >& v ){ return &v[0]; }
template< class C > const C* GetPointer( const std::vector< C >& v ){ return &v[0]; }

template< class C >       C* GetAddress(       C* c ) { return c; }
template< class C > const C* GetAddress( const C* c ) { return c; }
#endif // ARRAY_DEBUG

// [WARNING] The out-of-core array uses a buffer and is not thread-safe
template< class Data , size_t BufferSize=64 >
class OOCArray
{
	FILE* _fp;
	char _fileName[512];
	size_t _currentBufferIndex;
	Data _buffer[BufferSize];

	void _init( void )
	{
#ifdef _WIN32
		const char FileSeparator = '\\';
#else // !_WIN32
		const char FileSeparator = '/';
#endif // _WIN32

		const char prefix[] = "scratch";
		const char current[] = ".";
		{
			const char* dir = NULL;
			if( !dir || !strlen(dir) ) dir = getenv( "TMP" );	// Try to get the scratch directory from the environment variable
			if( !dir || !strlen(dir) ) dir = current;			// Set it to the current directory

			if( dir[strlen(dir)-1]!=FileSeparator ) sprintf( _fileName , "%s%c%s_XXXXXX" , dir , FileSeparator , prefix );
			else                                    sprintf( _fileName ,   "%s%s_XXXXXX" , dir ,                 prefix );
		}
#ifdef _WIN32
		_mktemp( _fileName );
		fopen_s( &_fp , _fileName , "w+b" );
#else // !_WIN32
		mktemp ( _fileName );
		_fp = fopen( _fileName , "w+b" );
#endif // _WIN32
		if( !_fp ) fprintf( stderr , "[ERROR] Failed to open file for writing: %s\n" , _fileName ) , exit( 0 );
		_currentBufferIndex = 0;
	}
public:
	OOCArray( void ){ _init(); }
	OOCArray( size_t sz ) { _init() , resize( sz ); }

	~OOCArray( void )
	{
		fclose( _fp );
		remove( _fileName );
		_fp = NULL;
	}
	void resize( size_t sz )
	{
#ifdef _WIN32
		_fseeki64( _fp , sizeof( Data ) * sz , SEEK_SET );
#else // !_WIN32
		fseek    ( _fp , sizeof( Data ) * sz , SEEK_SET );
#endif // _WIN32
	}
	Data& operator[]( size_t i )
	{
		size_t bufferIndex = i / BufferSize;
		if( bufferIndex!=_currentBufferIndex )
		{
			// Write out the previous buffer
#ifdef _WIN32
			_fseeki64( _fp , sizeof( Data ) * _currentBufferIndex * BufferSize , SEEK_SET );
#else // !_WIN32
			fseek    ( _fp , sizeof( Data ) * _currentBufferIndex * BufferSize , SEEK_SET );
#endif // _WIN32
			fwrite( &_buffer , sizeof( Data ) , BufferSize , _fp );

			// Read in the current buffer
			_currentBufferIndex = bufferIndex;
#ifdef _WIN32
			_fseeki64( _fp , sizeof( Data ) * _currentBufferIndex * BufferSize , SEEK_SET );
#else // !_WIN32
			fseek    ( _fp , sizeof( Data ) * _currentBufferIndex * BufferSize , SEEK_SET );
#endif // _WIN32
			fread( &_buffer , sizeof( Data ) , BufferSize , _fp );
		}
		 return _buffer[ i % BufferSize ];
	}
};

#endif // ARRAY_INCLUDED
