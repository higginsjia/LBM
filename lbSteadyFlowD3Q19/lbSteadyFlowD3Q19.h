#ifndef _LBSTEADYFLOWD3Q19_H_
#define _LBSTEADYFLOWD3Q19_H_

#include "assert.h"
#include<iostream>
#include<math.h>
#include<fstream>
#include<string>
#include<iomanip>
#include<sstream>
#include <iostream>
#include "omp.h"
using namespace std;

//not show warning C4996, as we use fopen() function in writeData() 
#pragma warning(disable:4996) 

/*
1. preReadSgn: get Nx, Ny, Nz
2. setSgnMem: new mem for sgn
3. readSgn: get info to sgn
4. 

*/


class lbSteadyFlowD3Q19
{
public:
	lbSteadyFlowD3Q19();
	~lbSteadyFlowD3Q19();

	//pre-processing 
	void preReadSgn(char* filename);
	void setSgnMem();
	void readSgn(char* filename);
	void setMem();

	//init
	void initSimuParam();
	void initSimuVarValue();
	void initDfsValue();

	//solver
	void streamDfs();

	void calMarcoVar();

	void collideDfs();

	void treatDfsBoundary();

	void getVelocityField(double* velocityField);
	int getSize();

	double calFeq(double density, double velocity[3],int i);
	double computeLocalRelaxationTime(double tau_, double stressTensorNorm, double smagConstant);

	void resetMem();

	

	FILE *fp;
	//post-processing
	void writeData(const char* filename);

	void Writechar(char t)
	{
		fwrite(&t, sizeof(char), 1, fp);
	};

	void Writeint(int t)
	{
		fwrite(&t, sizeof(int), 1, fp);
	};

	void Writeshort(short t)
	{
		fwrite(&t, sizeof(short), 1, fp);
	};

	void Writefloat(float t)
	{
		fwrite(&t, sizeof(float), 1, fp);
	};

	void Writedouble(double t)
	{
		fwrite(&t, sizeof(double), 1, fp);
	};

	void WriteNchar(char* chars, int length)
	{
		fwrite(chars, sizeof(char), length, fp);
	};

	void WriteString(char* str)
	{
		int t = 0;
		do
		{
			t = (*str);
			Writeint(t);
			str++;
		} while (t != 0);
	};

	
private:

	double gridLength;
	double timeStep;
	double velocityScale;

	int Nx;
	int Ny;
	int Nz;

	double fluidViscosity;
	double phyInletVelocity[3];//unit: m/s
	double latticeInletVelocity[3];//unit: lb

	double relaxationTime;
	double smagorinskyConstant;
	
	int FLUID;
	int INLET;
	int OUTLET;
	int NO_SLIP;

	double* collideField;
	double* streamField;

	int* flag;
	double* latticeUZ;
	double* latticeUY;
	double* latticeUX;
	double* latticeDensity;

	short* sgn;

	const int latticeDim = 19;
	static const double latticeWeight[19];
	static const int latticeVelocity[19][3];

	//assist function for index

	int cellIndex(int x, int y, int z);
	int funcIndex(int x, int y, int z, int i);
};


#endif