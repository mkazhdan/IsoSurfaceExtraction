/* -*- C++ -*- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


//////////
// Grid //
//////////
template<class Data>
Grid<Data>::Grid(void)
{
	resX=resY=0;
	values=NULL;
}
template<class Data>
Grid<Data>::~Grid(void){if(values){resize(0,0);}}
template<class Data>
int Grid<Data>::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Data>
int Grid<Data>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Data>
int Grid<Data>::read(FILE* fp){
	int io,resX,resY;
	io=int(fread(&resX,sizeof(int),1,fp));
	if(!io){return 0;}
	io=int(fread(&resY,sizeof(int),1,fp));
	if(!io){return 0;}
	resize(resX,resY);
	io=int(fread(values,sizeof(Data),resX*resY,fp));
	if(io==resX*resY){return 1;}
	else{return 0;}
}
template<class Data>
int Grid<Data>::write(FILE* fp) const {
	int io;
	io=int(fwrite(&resX,sizeof(int),1,fp));
	if(!io){return 0;}
	io=int(fwrite(&resY,sizeof(int),1,fp));
	if(!io){return 0;}
	io=int(fwrite(values,sizeof(Data),resX*resY,fp));
	if(io==resX*resY){return 1;}
	else{return 0;}
}
template<class Data>
void Grid<Data>::resolution(int& rX,int& rY) const
{
	rX=resX;
	rY=resY;
}
template<class Data>
int Grid<Data>::width(void) const {return resX;}
template<class Data>
int Grid<Data>::height(void) const {return resY;}


template<class Data>
int Grid<Data>::resize(const int& rX,const int& rY)
{
	if(rX<0 || rY<0){return 0;}
	else{
		if(values) delete[] values;
		values=NULL;
		resX=resY=0;
		if(rX && rY)
		{
			values=new Data[rX*rY];
			if(!values){return 0;}
			else{resX=rX;resY=rY;}
			clear();
		}
		return 1;
	}
}
template<class Data>
void Grid<Data>::clear(void){if(resX && resY){memset(values,0,sizeof(Data)*resX*resY);}}

template<class Data>
Data* Grid<Data>::operator[] (const int& i){return &values[i*resY];}
template<class Data>
Data& Grid<Data>::operator() (const int& i,const int& j){
	int x=i,y=j;
	if(x<0){x=resX-((-x)%resY);}
	x%=resX;
	if(y<0){y=resX-((-y)%resY);}
	y%=resY;
	return values[x*resY+y];
}
template<class Data>
const Data &Grid<Data>::operator() (const int& i,const int& j) const {
	int x=i,y=j;
	if(x<0){x=resX-((-x)%resX);}
	x%=resX;
	if(y<0){y=resY-((-y)%resY);}
	y%=resY;
	return values[x*resY+y];
}
template<class Data>
Data Grid<Data>::operator() (const double& x,const double& y){
	int x1,y1;
	double dx,dy;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}

	dx=x-x1;
	dy=y-y1;
	return (*this)(x1,y1)*(1.0-dx)*(1.0-dy)+(*this)(x1+1,y1)*dx*(1.0-dy)+(*this)(x1,y1+1)*(1.0-dx)*dy+(*this)(x1+1,y1+1)*dx*dy;
}
template<class Data>
Data Grid<Data>::operator() (const double& x,const double& y) const {
	int x1,y1;
	double dx,dy;

	if(x<0){x1=-int(-x)-1;}
	else{x1=int(x);}
	if(y<0){y1=-int(-y)-1;}
	else{y1=int(y);}

	dx=x-x1;
	dy=y-y1;
	return (*this)(x1,y1)*(Data(1.0)-dx)*(Data(1.0)-dy)+(*this)(x1+1,y1)*dx*(Data(1.0)-dy)+(*this)(x1,y1+1)*(Data(1.0)-dx)*dy+(*this)(x1+1,y1+1)*dx*dy;
}
template<class Data>
template<class Real>
Real Grid<Data>::squareNorm(void) const{return Dot<Real>(*this,*this);}
template<class Data>
template<class Real>
Real Grid<Data>::SquareDifference(const Grid& g1,const Grid& g2){return g1.squareNorm<Real>()+g2.squareNorm<Real>()-2*Dot<Real>(g1,g2);}
template<class Data>
template<class Real>
Real Grid<Data>::Dot(const Grid& g1,const Grid& g2){
	Real d=0;
	if(g1.resX != g2.resX || g1.resY != g2.resY)
	{
		fprintf(stderr,"Could not compare arrays of different sizes: (%d,%d) != (%d,%d)\n",g1.resX,g1.resY,g2.resX,g2.resY);
		exit(0);
	}
	for(int i=0;i<g1.resX*g1.resY;i++){d+=g1.values[i]*g2.values[i];}
	return Real(d/(g1.resX*g1.resY));
}
