/*
修改边界条件或平衡分布函数等 
*/
#include "stdafx.h"
#include<iostream>
#include<math.h>
#include<fstream>
#include<string>
#include<iomanip>
#include<sstream>
#include <algorithm>
#include "omp.h"
using namespace std;



/*
#pragma omp parallel for

*/
#define Nx 40
#define Ny 40

#define IN_X1 8
#define IN_X2 18
#define OUT_X1 25
#define OUT_X2 36

#define MaxIter 100000
#define LogStep 1000
const double w[9] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
const double cx[9] = { 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0 };
const double cy[9] = { 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0 };
const int inv[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

const double tau = 0.501;
const double omega = 1.0 / tau;
const double Uin = 0.05;

const int boundary = 1;
const int inlet = 3;
const int fluid = 0;
const int outlet = 4;

double f_old[Nx][Ny][9], f_new[Nx][Ny][9];
double u[Nx][Ny], v[Nx][Ny], rho[Nx][Ny];
int flag[Nx][Ny];

void init();
double feq(double rho_, double u_, double v_, int k);
void evolute();
void Out_file(int m);
double computeLocalRelaxationTime(double tau_, double stressTensorNorm, double smagConstant)
{
	const double viscosity = (tau_ - 0.5) / 3.0;
	const double smagSqr = smagConstant * smagConstant;
	const double stress =
		(std::sqrt(viscosity * viscosity + 18.0 * smagSqr * stressTensorNorm) - viscosity) /
		(6.0 * smagSqr);

	return 3.0 * (viscosity + smagSqr * stress) + 0.5;
}
void calMeanVelocity();
int _tmain(int argc, _TCHAR* argv[])
{
	init();

	for (int i = 0; i<=MaxIter; i++)
	{
		evolute();
		if (i%LogStep == 0)
		{
			printf("iter=%6d  ", i);
			calMeanVelocity();
		}
		
	}

	Out_file(999);

	return 0;
}


void init()
{
	int i, j, k;
	//set flag

	for (i = 0; i<Nx; i++)
	for (j = 0; j<Ny; j++)
	{
		if ((i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1))
			flag[i][j] = boundary;
		else if ((i >= IN_X1 && (i <= IN_X2)) && j == Ny - 2)
			flag[i][j] = inlet;
		else if ((i >= OUT_X1) && (i <= OUT_X2) && (j == 1))
			flag[i][j] = outlet;
		else
			flag[i][j] = fluid;
	}//for
	//
	//set initial value
	printf("init flag\n");

#pragma omp parallel for
	for (i = 0; i<Nx; i++)
	for (j = 0; j<Ny; j++)
	{
		u[i][j] = 0.0;
		v[i][j] = 0.0;
		rho[i][j] = 1.0;
		for (k = 0; k <= 8; k++)
		{
			f_old[i][j][k] = w[k];
			f_new[i][j][k] = f_old[i][j][k];
		}
	}
	printf("init Dfs\n");

}


void calMeanVelocity()
{
	int i, j;
	double sumU = 0.0;

	for (i = 1; i <Nx - 1; i++)
	for (j = 1; j < Ny - 1; j++)
	{
		if (flag[i][j] != boundary)
		{
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			rho[i][j] = 0.0;
			for (int k = 0; k <= 8; k++)
			{
				rho[i][j] += f_new[i][j][k];
				u[i][j] += cx[k] * f_new[i][j][k];
				v[i][j] += cy[k] * f_new[i][j][k];
			}

			sumU += sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]);


		}
	}

	sumU /= (double)(Nx*Ny);
	printf("MeanU=%e\n", sumU);

}

double feq(double rho_, double u_, double v_, int k)
{
	double t1 = u_* u_ + v_ * v_;
	double t2 = u_ * cx[k] + v_ * cy[k];
	return w[k] * (rho_ + 3.0*t2 + 4.50*t2*t2 - 1.50*t1);
}

void evolute()
{
	int i, j, k;
	//stream in fluid zone
//#pragma omp parallel for
	for (i = 0; i <Nx; i++)
	for (j = 0; j <Ny; j++)
	{
		if (flag[i][j] == fluid)
		{
			for (k = 0; k <= 8; k++)
			{
				int inv_ = inv[k];
					f_new[i][j][k] = f_old[i + (int)cx[inv_]][j + (int)cy[inv_]][k];
			}
		}
	}

	//calculate macroscopic var
#pragma omp parallel for
	for (i = 0; i <Nx ; i++)
	for (j = 0; j <Ny ; j++)
	{
			u[i][j] = 0.0;v[i][j] = 0.0;rho[i][j] = 0.0;

			for (k = 0; k <= 8; k++){
				rho[i][j] += f_new[i][j][k];
				u[i][j] += cx[k] * f_new[i][j][k];
				v[i][j] += cy[k] * f_new[i][j][k];
			}

			if (flag[i][j] == inlet)
			{
				rho[i][j] = rho[i][j - 1];
				u[i][j] = 0.0;
				v[i][j] = -Uin;
			}
			if (flag[i][j] == outlet)
			{
				u[i][j] = u[i][j + 1];
				v[i][j] = v[i][j + 1];
				rho[i][j] = 1.0;// rho[i][j + 1];
				
			}
			if (flag[i][j] == boundary)
			{
				u[i][j] = 0.0;
				v[i][j] = 0.0;
			}

			double pi_xx, pi_xy, pi_yy;
			pi_xx = 0.0;
			pi_xy = 0.0;
			pi_yy = 0.0;
			for (int it = 0; it < 9; ++it)
			{
				double fi_ = f_new[i][j][it];
				double f_eq = feq(rho[i][j], u[i][j], v[i][j], it);
				pi_xx += cx[it] * cx[it] * (fi_ - f_eq);
				pi_xy += cx[it] * cy[it] * (fi_ - f_eq);
				pi_yy += cy[it] * cy[it] * (fi_ - f_eq);

			}
			double stress = pi_xx*pi_xx + pi_xy*pi_xy*2.0 + pi_yy*pi_yy;
			double ss_stress = std::sqrt(stress);
			double ome = computeLocalRelaxationTime(tau, ss_stress, 0.12);
			double omegaC = 1.0 /ome;

			double omega_les = omegaC;

			for (k = 0; k <= 8; k++)
				f_old[i][j][k] = (1.0 - omega_les)*f_new[i][j][k]
					+ omega_les*feq(rho[i][j], u[i][j], v[i][j], k);

			//treat boundary for Dfs
			if (flag[i][j] == inlet)
			{

				for (int k = 0; k < 9; ++k)
				{
					if (flag[i + (int)cx[k]][j + (int)cy[k]] == fluid)
						f_old[i][j][k] = feq(rho[i][j], u[i][j], v[i][j], k) +
						f_old[i + (int)cx[k]][j + (int)cy[k]][k] -
						feq(rho[i + (int)cx[k]][j + (int)cy[k]],
						u[i + (int)cx[k]][j + (int)cy[k]],
						v[i + (int)cx[k]][j + (int)cy[k]],
						k);
				}
			}

			if (flag[i][j] == outlet)
			{
				rho[i][j] = 1.0;
				for (int k = 0; k < 9; ++k)
				{
					if (flag[i + (int)cx[k]][j + (int)cy[k]] == fluid)
						f_old[i][j][k] = feq(rho[i][j], u[i][j], v[i][j], k) +
						f_old[i + (int)cx[k]][j + (int)cy[k]][k] -
						feq(rho[i + (int)cx[k]][j + (int)cy[k]],
						u[i + (int)cx[k]][j + (int)cy[k]],
						v[i + (int)cx[k]][j + (int)cy[k]],
						k);
				}
			}

			if (flag[i][j] == boundary)
			{
				//rho[i][j] = 1.0;
				int ip = i + (int)cx[i];
				int jp = j + (int)cy[i];

				if (i >= 0 && i < Nx&&j >= 0 && j < Ny)
				{
					if (flag[i + (int)cx[k]][j + (int)cy[k]] == fluid)
						f_old[i][j][k] = feq(rho[i + (int)cx[k]][j + (int)cy[k]], 0.0, 0.0, k)+
						f_old[i + (int)cx[k]][j + (int)cy[k]][k] -
						feq(rho[i + (int)cx[k]][j + (int)cy[k]],
						u[i + (int)cx[k]][j + (int)cy[k]],
						v[i + (int)cx[k]][j + (int)cy[k]],
						k);
						
				}
			}


	}


	

	//


}

void Out_file(int m)
{
	ostringstream name;
	name << "D:\\output\\cavityA_" << setfill('0') << setw(8) << m << ".plt";
	ofstream out(name.str().c_str());
	out << "Title=\"LBM Lid Driven Flow\"\n"
		<< "VARIABLES=X,Y,U,V,UV,Bound,Rho\n"
		<< "ZONE T= \"BOX\", I= "
		<< Nx  << ", J=" << Ny  << ", F=POINT" << endl;
	int i, j;
	for (j = 0; j <Ny ; j++)
	for (i = 0; i <Nx ; i++)
	{

		out << double(i) << " " << double(j) << " "
			<< u[i][j] << " " << v[i][j] << " " << sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]) << " " << flag[i][j] << " " << rho[i][j] << endl;
	}
	printf("write file end.\n");
}
