

#include "stdafx.h"
#include<iostream>
#include<math.h>
#include<fstream>
#include<string>
#include<iomanip>
#include<sstream>
#include <algorithm>


#include "lbSteadyFlow.h"
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	lbSteadyFlow2d A;
	
	int Nx = 50;
	int Ny = 100;

	int in_x1 = 10;
	int in_x2 = 20;

	int out_x1 = 35;
	int out_x2 = 47;

	double dx = 0.01;
	double Uin[2] = { 0.0, -5.0 };
	double visc = 1.0e-6;
	double CS = 0.13;

	A.setPhysicalParam(Nx, Ny, in_x1, in_x2, out_x1, out_x2, dx, visc, Uin, CS);

	double* Ucur = (double*)malloc(sizeof(double)*Nx*Ny);
	double* Unex = (double*)malloc(sizeof(double)*Nx*Ny);

	A.setMem();

	A.init();

	int MaxIter = 100000;
	int logIter = 1000;
	double velocityAccuracy = 1.0e-6;

	for (int iter = 1; iter <= MaxIter; ++iter)
	{

		A.getVelocityField(Ucur);

		A.streamDFs();

		A.calMarcoVar();

		A.collideDFs();

		A.treatDFsBoundary();

		A.getVelocityField(Unex);

		double u_err = 0.0;
		double u_sum = 0.0;

		for (int i = 0; i < Nx*Ny; ++i)
		{
			u_err += fabs(Ucur[i] - Unex[i]);
			u_sum += Ucur[i];
		}

		u_err /= u_sum;

		if (iter % logIter == 0)
		{

			printf("iter=%6d err=%e\n", iter,u_err);
		}

		if (u_err<velocityAccuracy&&iter>logIter)
		{
			printf("iter=%6d err=%e\n", iter, u_err);
			break;
		}

	}

	A.Out_file(666);
	
	A.resetMem();


	return 0;
}

