#ifndef COMMON
#define COMMON

#include <memory.h>
#include <stdio.h>
#include "math.h"
//#include <direct.h>
#include "stdlib.h"
#include "mkl.h"
#include <assert.h>
#include <complex>
#include <time.h>

//const double SIGMIN=-1.0e+10; //minimum constant
const double SIGMAX=1.0e+7;     //maximum constant
const double EPS=1.0e-7;        //accuracy control constant
const double PI=3.1415926535897932384626433832795;  //PI
const double Z0=120*PI;                             //wave impedance of free space
const double V0=3.0*pow(10.0,8.0);                  //speed of light

//allocate memory for 1D array
template<class Type> Type* Memory1D(Type* a, const int m)
{
	assert(m>0);
	a=new Type[m];
	return a;
}

//release memory for 1D array
template<class Type> Type* RMemory1D(Type* a)
{
	if (a!=NULL)
	{
		delete[] a;
		a=NULL;
	}

	return a;
}


//allocate memory for 2D array
template<class Type> Type** Memory2D(Type** a, const int m, const int n)
{
	assert (m>0);
	assert (n>0);

	a=new Type*[m];
	a[0]=new Type[m*n];  //continuous address
	for (int i=0; i<m; i++)
	{
		a[i]=a[0]+i*n;
	}
	return a;
}

//release memory for 2D array
template<class Type> Type** RMemory2D(Type** a)
{
	assert (a!=NULL);
	delete[] a[0];  //continuous address
	delete[] a;
	a=NULL;
	return a;
}

//allocate memory for 3D array
template<class Type> Type*** Memory3D(Type*** a, const int m, const int n, const int l)
{
    assert (m>0);
	assert (n>0);
	assert (l>0);

	int i,j;
	a=new Type**[m];      //layer 1
	a[0]=new Type*[m*n];  //layer 2
    //layer 2
	for ( i=0; i<m; i++)
	{
		a[i]=a[0]+n*i;
	}

	//layer 3
	a[0][0]=new Type[m*n*l];  //continuous address
    //layer 3
	for ( i=0; i<m; i++)
	{
		for ( j=0; j<n; j++)
		{
			a[i][j]=a[0][0]+(i*n+j)*l;
		}
	}
	return a;
}

//release memory for 3D array
template<class Type> Type*** RMemory3D(Type*** a)
{
	assert (a!=NULL);
	delete[] a[0][0];  //data
	delete[] a[0];     //layer 2
	delete[] a;        //layer 1

	a=NULL;
	return a;
}


#endif
