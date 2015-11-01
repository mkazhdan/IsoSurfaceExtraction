#ifndef GRID_INCLUDED
#define GRID_INCLUDED

#include <stdio.h>

template<class Data>
class Grid
{
protected:
	Data* values;
	int resX,resY;
public:
	Grid(void);
	~Grid(void);

	// Returns the dimension of the array
	void resolution(int& resX,int& resY) const;
	int width(void) const;
	int height(void) const;
	// Allocates memory for the array
	int resize(const int& resX,const int& resY);

	// Clears the values of the array to 0
	void clear(void);

	// Returns a reference to the indexed  array element
	Data& operator() (const int& x,const int& y);
	const Data& operator() (const int& x,const int& y) const;
	Data* operator[] (const int& x);
	// Returns the linear interpolation of the value at the spedified index
	Data operator() (const double& x,const double& y);
	Data operator() (const double& x,const double& y) const;

	// Returns the square of the L2-norm of the array elements
	template<class Real>
	Real squareNorm(void) const;

	// Reads in an array from the specified file
	int read(const char* fileName);
	int read(FILE* fp);

	// Writes out the array to the specified file
	int write(const char* fileName) const;
	int write(FILE* fp) const;

	// Returns the square of the L2-difference between two arrays
	template<class Real>
	static Real SquareDifference(const Grid& g1,const Grid& g2);

	// Returns the dot-product of two arrays
	template<class Real>
	static Real Dot(const Grid& g1,const Grid& g2);
};

#include "Grid.inl"
#endif // GRID_INCLUDED
