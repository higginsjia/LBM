#ifndef _LBSTEADY_FLOW_H_
#define _LBSTEADY_FLOW_H_

/*
DESCRIPTION:
LBM simulate steady flow with LES model
constant velocity inlet and pressure outlet
results output as ACSII Tecplot format

ATTENTION:
pressure outlet may need to be revised
LES coupled LBM model may not accurate for pressure field
and the coutour line of pressure or velocity may not accurate
FUTURE WORK:
use openGL as real time animation and refine the post-processing
*/

#include "assert.h"
#include<iostream>
#include<math.h>
#include<fstream>
#include<string>
#include<iomanip>
#include<sstream>
#include <iostream>

using namespace std;

class lbSteadyFlow2d
{
public:
	lbSteadyFlow2d();
	~lbSteadyFlow2d();

	void setMem();
	void resetMem();
	int cellIndex(int x, int y);
	int funcIndex(int x, int y, int i);
	
	void init();
	void setPhysicalParam(int nx, int ny, int in_x1, int in_x2, int out_x1, int out_x2,
		double dx, double fluidVisc, double Uin[2], double smagCs);
	double calFeq(double density, double velocity[2], int i);
	double computeLocalRelaxationTime(double tau_, double stressTensorNorm, double smagConstant);
	void streamDFs();
	void calMarcoVar();
	void collideDFs();
	void treatDFsBoundary();
	void Out_file(int m);
	void getVelocityField(double* curU);

private:
	double dX;
	double dT;
	double K_scale;

	double fluidVelocity[2];
	double fluidViscosity;

	int Nx;
	int Ny;

	int totalCellNum;

	double velocityIn[2];
	double tau;
	double omega;
	double CS;

	int IN_X1;
	int IN_X2;

	int OUT_X1;
	int OUT_X2;


	double* collideField;
	double* streamField;


	int* flag;
	double* Ux;
	double* Uy;
	double* Density;


	double weight[9];
	double cx[9];
	double cy[9];
	int inv[9];

	int FLUID;
	int INLET;
	int OUTLET;
	int NO_SLIP;
};


lbSteadyFlow2d::lbSteadyFlow2d()
{
	for (int i = 0; i < 9; ++i)
	{
		if (i == 0)
			weight[i] = 4.0 / 9.0;
		else if (i >= 1 && i <= 4)
			weight[i] = 1.0 / 9.0;
		else
			weight[i] = 1.0 / 36.0;
	}

	//cx[9] = { 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0 };
	cx[0] = 0.0;	cx[1] = 1.0;	cx[2] = 0.0;
	cx[3] = -1.0;	cx[4] = 0.0;	cx[5] = 1.0;
	cx[6] = -1.0;	cx[7] = -1.0;	cx[8] = 1.0;

	//cy[9] = { 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0 };
	cy[0] = 0.0; cy[1] = 0.0; cy[2] = 1.0;
	cy[3] = 0.0; cy[4] = -1.0; cy[5] = 1.0;
	cy[6] = 1.0; cy[7] = -1.0; cy[8] = -1.0;

	//inv[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
	inv[0] = 0; inv[1] = 3; inv[2] = 4;
	inv[3] = 1; inv[4] = 2; inv[5] = 7;
	inv[6] = 8; inv[7] = 5; inv[8] = 6;



	FLUID = 0;
	INLET = 1;
	OUTLET = 2;
	NO_SLIP = 3;


	Nx = 40;
	Ny = 40;

	IN_X1 = 8;
	IN_X2 = 18;

	OUT_X1 = 25;
	OUT_X2 = 36;

	tau = 0.501;
	omega = 1.0 / tau;

	velocityIn[0] = 0.0;
	velocityIn[1] = -0.05;

	CS = 0.12;

	K_scale = 100.0;
}

int lbSteadyFlow2d::cellIndex(int x, int y)
{
	assert(x >= 0 && x < Nx);
	assert(y >= 0 && y < Ny);

	return y + x*Ny;
}

int lbSteadyFlow2d::funcIndex(int x, int y, int i)
{
	assert(x >= 0 && x < Nx);
	assert(y >= 0 && y < Ny);
	assert(i >= 0 && i < 9);

	return (y+x*Ny)*9+i;
}

void lbSteadyFlow2d::setPhysicalParam(int nx, int ny,int in_x1,int in_x2,int out_x1,int out_x2,
	double dx, double fluidVisc, double Uin[2], double smagCs)
{
	Nx = nx;
	Ny = ny;

	IN_X1 = in_x1;
	IN_X2 = in_x2;

	OUT_X1 = out_x1;
	OUT_X2 = out_x2;


	dX = dx;
	fluidViscosity = fluidVisc;

	fluidVelocity[0] = Uin[0];
	fluidVelocity[1] = Uin[1];

	CS = smagCs;

	dT = dx / K_scale;
	tau = 3.0*fluidViscosity*dT / dX / dX + 0.5;


	printf("NxxNy:%dx%d\ninlet:%d-%d\noutlet:%d-%d\n",
		Nx, Ny, IN_X1, IN_X2, OUT_X1, OUT_X2);

	printf("tau=%e K=%6.f Uin=(%4.2f,%4.2f)\n",
		tau, K_scale, fluidVelocity[0], fluidVelocity[1]);
}

void lbSteadyFlow2d::getVelocityField(double* curU)
{
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	{
		const int Xcell = cellIndex(x, y);

		curU[Xcell] = sqrt(Ux[Xcell] * Ux[Xcell] + Uy[Xcell] * Uy[Xcell]);
	}

}

void lbSteadyFlow2d::resetMem()
{
	delete[] collideField;
	delete[] streamField;
	delete[] flag;
	delete[] Ux;
	delete[] Uy;
	delete[] Density;

}

void lbSteadyFlow2d::setMem()
{
	const int DfNum = Nx*Ny * 9;
	const int CellNum = Nx*Ny;

	collideField = new double[DfNum];
	streamField = new double[DfNum];

	flag = new int[CellNum];
	Ux = new double[CellNum];
	Uy = new double[CellNum];
	Density = new double[CellNum];


	memset(collideField, 0, sizeof(double)*DfNum);
	memset(streamField, 0, sizeof(double)*DfNum);

	memset(flag, 0, sizeof(int)*CellNum);
	memset(Ux, 0, sizeof(double)*CellNum);
	memset(Uy, 0, sizeof(double)*CellNum);
	memset(Density, 0, sizeof(double)*CellNum);

}

void lbSteadyFlow2d::init()
{

	//init flag
	printf("初始化节点标识开始\n");

	for (int x = 0; x < Nx;++x)
	for (int y = 0; y < Ny; ++y)
	{
		int idx = cellIndex(x, y);
		if (x == 0 || x == Nx - 1 || y == 0 || y == Ny - 1)
			flag[idx] = NO_SLIP;
		else if (x >= IN_X1&&x <= IN_X2&&y == Ny - 2)
			flag[idx] = INLET;
		else if (x >= OUT_X1&&x <= OUT_X2&&y == 1)
			flag[idx] = OUTLET;
		else
			flag[idx] = FLUID;
	}

	printf("初始化节点标识结束\n");

	printf("初始化分布函数和宏观变量开始\n");
#pragma omp parallel for
	for (int x = 0; x < Nx;++x)
	for (int y = 0; y < Ny; ++y)
	{
		const int nodeIdx = cellIndex(x, y);
		Ux[nodeIdx] = 0.0;
		Uy[nodeIdx] = 0.0;
		Density[nodeIdx] = 1.0;

		for (int k = 0; k < 9; ++k)
		{
			int idx = funcIndex(x, y, k);
			collideField[idx] = weight[k];
			streamField[idx] = weight[k];
		}
	}
	printf("初始化分布函数和宏观变量结束\n");
}

void lbSteadyFlow2d::streamDFs()
{
#pragma omp parallel for
	for (int x = 0; x < Nx;++x)
	for (int y = 0; y < Ny; ++y)
	{
		const int curCell = cellIndex(x, y);
		if (flag[curCell] == FLUID)//cal fluid only
		{
			for (int i = 0; i < 9; ++i)
			{
				int inv_ = inv[i];
				int cellX = funcIndex(x, y, i);
				int cellXsubCi = funcIndex(
					x + (int)cx[inv_], 
					y + (int)cy[inv_], 
					i);
				streamField[cellX] = collideField[cellXsubCi];
				
			}//i=0-8
		}//if fluid
	}//out for

	//copy Dfs
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	for (int i = 0; i < 9; ++i)
	{
		int idx = funcIndex(x, y, i);
		collideField[idx] = streamField[idx];
	}
}

void lbSteadyFlow2d::calMarcoVar()
{
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	{
		const int Xcell = cellIndex(x, y);

		Ux[Xcell] = 0.0; Uy[Xcell] = 0.0; Density[Xcell] = 0.0;
		for (int i = 0; i < 9; ++i)
		{
			const int DfIdx = funcIndex(x, y, i);
			double fi = collideField[DfIdx];
			Density[Xcell] += fi;
			Ux[Xcell] += cx[i] * fi;
			Uy[Xcell] += cy[i] * fi;
		}

		//set boundary value for marco var
		if (flag[Xcell] == INLET)
		{
			Density[Xcell] = Density[cellIndex(x, y - 1)];
			Ux[Xcell] = velocityIn[0];
			Uy[Xcell] = velocityIn[1];
		}
		else if (flag[Xcell] == OUTLET)
		{
			Density[Xcell];// = 1.0;//outlet pressure
			Ux[Xcell] = Ux[cellIndex(x, y + 1)];
			Uy[Xcell] = Uy[cellIndex(x, y + 1)];
		}
		else if (flag[Xcell] == NO_SLIP)
		{
			Ux[Xcell] = 0.0;
			Uy[Xcell] = 0.0;
		}

	}
}

double lbSteadyFlow2d::calFeq(double density, double velocity[2], int i)
{
	double cu = cx[i] * velocity[0] + cy[i] * velocity[1];
	double uu = velocity[0] * velocity[0] + velocity[1] * velocity[1];

	double feq = weight[i] * (density + 3.0*cu + 4.5*cu*cu - 1.5*uu);

	return feq;
}

double lbSteadyFlow2d::computeLocalRelaxationTime(double tau_, double stressTensorNorm, double smagConstant)
{
	const double viscosity = (tau_ - 0.5) / 3.0;
	const double smagSqr = smagConstant * smagConstant;
	const double stress =
		(std::sqrt(viscosity * viscosity + 18.0 * smagSqr * stressTensorNorm) - viscosity) /
		(6.0 * smagSqr);

	return 3.0 * (viscosity + smagSqr * stress) + 0.5;
}

void lbSteadyFlow2d::collideDFs()
{
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	{
		const int Xcell = cellIndex(x, y);

		if (flag[Xcell] == FLUID)
		{

			double PI_XX, PI_XY, PI_YY;
			PI_XX = 0.0;
			PI_XY = 0.0;
			PI_YY = 0.0;

			double cellVelo[2] = { Ux[Xcell], Uy[Xcell] };

			for (int i = 0; i < 9; ++i)
			{
				double fi = collideField[funcIndex(x, y, i)];

				double feq = calFeq(Density[Xcell], cellVelo, i);

				PI_XX += cx[i] * cx[i] * (fi - feq);
				PI_XY += cx[i] * cy[i] * (fi - feq);
				PI_YY += cy[i] * cy[i] * (fi - feq);
			}

			double stress = sqrt(PI_XX*PI_XX + 2.0*PI_XY*PI_XY + PI_YY*PI_YY);

			double tau_LES = computeLocalRelaxationTime(tau, stress, CS);

			double omega_LES = 1.0 / tau_LES;

			for (int i = 0; i < 9; ++i)
			{
				const int idx = funcIndex(x, y, i);
				collideField[idx] = (1.0 - omega_LES)*collideField[idx] +
					omega_LES*calFeq(Density[Xcell], cellVelo, i);
			}

		}
	}
}

void lbSteadyFlow2d::treatDFsBoundary()
{
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	{
		int curIdx = cellIndex(x, y);

		if (flag[curIdx] == INLET)
		{
			for (int i = 0; i < 9; ++i)
			{
				int Xci = x + (int)cx[i];
				int Yci = y + (int)cy[i];
				int nearIdx = cellIndex(Xci, Yci);

				double curU[2] = { Ux[curIdx], Uy[curIdx] };
				double nearU[2] = { Ux[nearIdx], Uy[nearIdx] };

				if (flag[nearIdx]==FLUID)
				{
					collideField[funcIndex(x, y, i)] =
						calFeq(Density[curIdx], curU, i) + collideField[funcIndex(Xci,Yci,i)] - calFeq(Density[nearIdx], nearU, i);
				}
			}
		}

		if (flag[curIdx] == OUTLET)
		{
			Density[curIdx] = 1.0;

			for (int i = 0; i < 9; ++i)
			{
				int Xci = x + (int)cx[i];
				int Yci = y + (int)cy[i];
				int nearIdx = cellIndex(Xci, Yci);

				double curU[2] = { Ux[curIdx], Uy[curIdx] };
				double nearU[2] = { Ux[nearIdx], Uy[nearIdx] };

				if (flag[nearIdx] == FLUID)
				{
					collideField[funcIndex(x, y, i)] =
						calFeq(Density[curIdx], curU, i) + collideField[funcIndex(Xci, Yci, i)] - calFeq(Density[nearIdx], nearU, i);
				}
			}
		}


		if (flag[curIdx] == NO_SLIP)
		{
			for (int i = 0; i < 9; ++i)
			{
				int Xci = x + (int)cx[i];
				int Yci = y + (int)cy[i];

				if (Xci >= 0 && Xci < Nx&&Yci >= 0 && Yci < Ny)
				{
					int nearIdx = cellIndex(Xci, Yci);

					double curU[2] = { 0.0, 0.0 };
					double nearU[2] = { Ux[nearIdx], Uy[nearIdx] };

					if (flag[nearIdx] == FLUID)
					{
						collideField[funcIndex(x, y, i)] =
							calFeq(Density[nearIdx], curU, i) + collideField[funcIndex(Xci, Yci, i)] - calFeq(Density[nearIdx], nearU, i);
					}
				}
			}
		}

	}
	
}

void lbSteadyFlow2d::Out_file(int m)
{
	ostringstream name;
	name << "D:\\output\\cavityC_" << setfill('0') << setw(8) << m << ".plt";
	ofstream out(name.str().c_str());
	out << "Title=\"LBM Lid Driven Flow\"\n"
		<< "VARIABLES=X,Y,U,V,UV,Bound,Rho\n"
		<< "ZONE T= \"BOX\", I= "
		<< Nx << ", J=" << Ny << ", F=POINT" << endl;
	int i, j;
	for (j = 0; j <Ny; j++)
	for (i = 0; i <Nx; i++)
	{
		int idx = cellIndex(i, j);
		out << double(i) << " " << double(j) << " "
			<< K_scale*Ux[idx] << " " << K_scale*Uy[idx] << " " 
			<< K_scale*sqrt(Ux[idx] * Ux[idx] + Uy[idx] * Uy[idx])
			<< " " << flag[idx] << " " << Density[idx] << endl;
	}
	printf("write file end.\n");
}

lbSteadyFlow2d::~lbSteadyFlow2d()
{
}



#endif