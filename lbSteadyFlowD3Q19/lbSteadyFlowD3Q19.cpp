#include "lbSteadyFlowD3Q19.h"


int lbSteadyFlowD3Q19::cellIndex(int x, int y, int z)
{
	assert(x >= 0 && x < Nx);
	assert(y >= 0 && y < Ny);
	assert(z >= 0 && z < Nz);

	//check index according to the max index equal N-1
	return z + y*Nz + x*Ny*Nz;
}

int lbSteadyFlowD3Q19::funcIndex(int x, int y, int z, int i)
{
	assert(x >= 0 && x < Nx);
	assert(y >= 0 && y < Ny);
	assert(z >= 0 && z < Nz);
	assert(i >= 0 && i < latticeDim);

	//check index carefully
	return (z + y*Nz + x*Ny*Nz)*latticeDim + i;
}

void lbSteadyFlowD3Q19::preReadSgn(char* filename)
{
	short nx, ny, nz;

	FILE* fp;
	fp = fopen(filename, "rb");

	if (fp == NULL) { std::cout << "Can not open SGN file" << std::endl; }
	fseek(fp, 12, SEEK_SET);
	fread(&nx, 2, 1, fp);
	fread(&ny, 2, 1, fp);
	fread(&nz, 2, 1, fp);

	Nx = nx;
	Ny = ny;
	Nz = nz;

	fclose(fp);
}

void lbSteadyFlowD3Q19::setSgnMem()
{
	const int cellNum = Nx*Ny*Nz;
	sgn = new short[cellNum];
	memset(sgn, 0, sizeof(short)*cellNum);
}

void lbSteadyFlowD3Q19::readSgn(char* filename)
{
	short nx, ny, nz;

	FILE* fp;
	fp = fopen(filename, "rb");
	if (fp == NULL) { std::cout << "Can not open SGN file" << std::endl; }
	fseek(fp, 12, SEEK_SET);
	fread(&nx, 2, 1, fp);
	fread(&ny, 2, 1, fp);
	fread(&nz, 2, 1, fp);

	int L = nx;
	int M = ny;
	int N = nz;
	int MN = M*N;
	int LMN = L*MN;
	fseek(fp, 238, SEEK_CUR);

	if (sgn == NULL) { std::cout << "ReadSGN: Not enough memory! " << std::endl; }
	std::cout << "Read sgn file" << std::endl;
	for (int i = 0; i<LMN; i++)
	{
		short m;
		fread(&m, 2, 1, fp);
		sgn[i] = m;
	}
	float dxxx = 0;
	fread(&dxxx, 4, 1, fp);
	float dx = dxxx / 1000;
	gridLength = (double)dx;//set gridLength
	std::cout << "dx = " << dx << std::endl;

	printf("L=%d M=%d N=%d\n", L, M, N);
	printf("sgn node=%d\n", L*M*N);
	fclose(fp);
}

void lbSteadyFlowD3Q19::setMem()
{
	const int cellNum = Nx*Ny*Nz;
	const int DfNum = cellNum*latticeDim;

	collideField = new double[DfNum];
	streamField = new double[DfNum];

	flag = new int[cellNum];
	latticeUX = new double[cellNum];
	latticeUY = new double[cellNum];
	latticeUZ = new double[cellNum];
	latticeDensity = new double[cellNum];

	memset(collideField, 0, sizeof(double)*DfNum);
	memset(streamField, 0, sizeof(double)*DfNum);

	memset(flag, 0, sizeof(int)*cellNum);
	memset(latticeUX, 0, sizeof(double)*cellNum);
	memset(latticeUY, 0, sizeof(double)*cellNum);
	memset(latticeUZ, 0, sizeof(double)*cellNum);
	memset(latticeDensity, 0, sizeof(double)*cellNum);

}

void lbSteadyFlowD3Q19::resetMem()
{
	delete[] collideField;
	delete[] streamField;

	delete[] flag;
	delete[] latticeUX;
	delete[] latticeUY;
	delete[] latticeUZ;
	delete[] latticeDensity;

	delete[] sgn;
}

void lbSteadyFlowD3Q19::initSimuParam()
{
	velocityScale = 100.0;
	timeStep = gridLength / velocityScale;


	//we should set phyInletU
	for (int i = 0; i < 3; ++i)
		latticeInletVelocity[i] = phyInletVelocity[i] / velocityScale;

	relaxationTime = 3.0*fluidViscosity*timeStep / gridLength / gridLength + 0.5;


	printf("inlet velocity physical=(%6.4f,%6.4f,%6.4f)m/s\n",
		phyInletVelocity[0], phyInletVelocity[1], phyInletVelocity[2]);

	printf("inlet velocity LB=(%6.4f,%6.4f,%6.4f)\n\n",
		latticeInletVelocity[0], latticeInletVelocity[1], latticeInletVelocity[2]);

	printf("grid size=%6.4f m\ntime step=%e s\n\n",
		gridLength, timeStep);

	printf("fluid viscosity=%e m2/s\ntau=%e\nLES constant=%4.3f\n\n",
		fluidViscosity, relaxationTime, smagorinskyConstant);




}

void lbSteadyFlowD3Q19::initSimuVarValue()
{
	//first we set cell type
	const int cellNum = Nx*Ny*Nz;
	int inletNumFirst = 0;//log inlet node num, sometimes no inlet cell is set in sgn grid
	for (int i = 0; i < cellNum; ++i)
	{
		const int curSgn = sgn[i];
		if (curSgn % 100 == 0 && curSgn != 100)
			flag[i] = FLUID;
		else if (curSgn == 100)
		{
			flag[i] = INLET;
			inletNumFirst++;
		}
		else
			flag[i] = NO_SLIP;
	}

	if (inletNumFirst == 0)
	{
		printf("no inlet cell is set!\n");

		//now we find inlet as top fluid cell
		int findInlet = false;
		int inletPosZ;
		for (int z = Nz - 1; z >= 0 && !findInlet; z--)
		for (int x = 0; x < Nx&&!findInlet; ++x)
		for (int y = 0; y < Ny&&!findInlet; ++y)
		{
			const int idx = cellIndex(x, y, z);
			if (sgn[idx] % 100 == 0)
			{
				inletPosZ = z;
				findInlet = true;
			}
		}
		//now we set inlet cell on top fluid cell
		for (int x = 0; x < Nx; ++x)
		for (int y = 0; y < Ny; ++y)
		{
			const int z = inletPosZ;
			const int idx = cellIndex(x, y, z);
			if (sgn[idx] % 100 == 0)
				flag[idx] = INLET;
		}
	}

	//now the riser is set as fluid, we set outlet according to sgn

	int findTopRiser = false;
	int topRiserZ;
	for (int z = Nz - 1; z >= 0 && !findTopRiser; z--)
	for (int x = 0; x < Nx&&!findTopRiser; ++x)
	for (int y = 0; y < Ny&&!findTopRiser; ++y)
	{
		const int idx = cellIndex(x, y, z);

		if (sgn[idx] == 4000)
		{
			topRiserZ = z;
			findTopRiser = true;
		}
	}

	printf("riser top z index=%d\n", topRiserZ);

	int outletNum = 0;
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	{
		int z = topRiserZ;
		const int idx = cellIndex(x, y, z);
		if (sgn[idx] == 4000)
		{
			flag[idx] = OUTLET;
			outletNum++;
		}
	}

	printf("outlet node num is %d\n", outletNum);

	//now we cal inlet z index and count total num of diff type cell
	int inletIndexZ;
	int inletNumber = 0;
	int outletNumber = 0;
	int fluidNumber = 0;
	int noSlipNumber = 0;

	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	for (int z = 0; z < Nz; ++z)
	{
		const int idx = cellIndex(x, y, z);
		const int curFlag = flag[idx];
		if (curFlag == FLUID)
			fluidNumber++;
		else if (curFlag == INLET)
		{
			inletNumber++;
			inletIndexZ = z;
		}
		else if (curFlag == OUTLET)
			outletNumber++;
		else if (curFlag == NO_SLIP)
			noSlipNumber++;
	}
	int totalNum = inletNumber + outletNumber + noSlipNumber + fluidNumber;
	printf("cell type count:\n");
	printf("inlet=%d\noutlet=%d\nnoSlip=%d\nfluid=%d\ncell sum=%d\n",
		inletNumber, outletNumber, noSlipNumber, fluidNumber, totalNum);
	printf("inlet Z index=%d\n", inletIndexZ);

	if (totalNum == Nx*Ny*Nz)
		printf("cell type sum right.\n");
	else
		printf("cell type sum WRONG!\n");
}

void lbSteadyFlowD3Q19::writeData(const char* filename)
{
	int varnum = 10;
	fp = fopen(filename, "wb");
	assert(fp != NULL);
	WriteNchar("#!TDV102", 8);
	Writeint(1);
	WriteString("");//The TITLE.

	Writeint(varnum);//Number of variables (NumVar) in the datafile.

	WriteString("I");
	WriteString("J");
	WriteString("K");

	WriteString("U");
	WriteString("V");
	WriteString("W");

	WriteString("Velovity");
	WriteString("VOF");
	WriteString("Pressure");
	WriteString("CellType");

	//WriteString("mat");

	Writefloat(299.0);
	WriteString("ZONE 001");//Zone name.
	Writeint(-1);
	Writeint(0);//Zone type: ordered
	Writedouble(0.0);
	Writefloat(0.0);

	Writeint(Nx);
	Writeint(Ny);
	Writeint(Nz);

	Writeint(0);
	Writefloat(357.0);
	Writefloat(299.0);
	for (int it = 0; it<varnum; it++)
	{
		Writeint(1);//Variable data format, 1=Float
	}
	Writeint(0);
	Writeint(-1);


	//write i,j,k
	//X index
	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		Writefloat((float)x);
	}

	//Y index
	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		Writefloat((float)y);
	}

	//Z index
	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		Writefloat((float)z);

	}

	//Ux
	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		const int idx = cellIndex(x, y, z);
		Writefloat((float)(velocityScale *latticeUX[idx]));
	}

	//Uy
	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		const int idx = cellIndex(x, y, z);
		Writefloat((float)(velocityScale* latticeUY[idx]));
	}

	//Uz
	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		const int idx = cellIndex(x, y, z);
		Writefloat((float)(velocityScale* latticeUZ[idx]));
	}

	//UU
	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		const int idx = cellIndex(x, y, z);
		double uMag = sqrt(
			latticeUX[idx] * latticeUX[idx] +
			latticeUY[idx] * latticeUY[idx] +
			latticeUZ[idx] * latticeUZ[idx]);

		Writefloat((float)(velocityScale*uMag));
	}

	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		const int idx = cellIndex(x, y, z);
		Writefloat((float)1.0f);
	}

	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		const int idx = cellIndex(x, y, z);
		Writefloat((float)latticeDensity[idx]);
	}

	for (int z = 0; z < Nz; ++z)
	for (int y = 0; y < Ny; ++y)
	for (int x = 0; x < Nx; ++x)
	{
		const int idx = cellIndex(x, y, z);
		Writefloat((float)flag[idx]);
	}

	printf("write data finished.\n");

	fclose(fp);
}

void lbSteadyFlowD3Q19::initDfsValue()
{
	printf("init DFs begin.\n");

	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	for (int z = 0; z < Nz; ++z)
	{
		const int cellIdx = cellIndex(x, y, z);

		latticeUX[cellIdx] = 0.0;
		latticeUY[cellIdx] = 0.0;
		latticeUZ[cellIdx] = 0.0;

		latticeDensity[cellIdx] = 1.0;//set 1.0 check here

		for (int i = 0; i < latticeDim; ++i)
		{
			const int idx = funcIndex(x, y, z, i);

			streamField[idx] = latticeWeight[i];
			collideField[idx] = latticeWeight[i];
		}
	}

	printf("init DFs finished.\n");
}

void lbSteadyFlowD3Q19::streamDfs()
{
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	for (int z = 0; z < Nz; ++z)
	{
		if (flag[cellIndex(x, y, z)] == FLUID)
		{
			for (int i = 0; i < latticeDim; ++i)
			{
				const int inv = latticeDim - 1 - i;
				const int cellX = funcIndex(x, y, z, i);
				const int cellXsubCi = funcIndex(
					x + latticeVelocity[inv][0],
					y + latticeVelocity[inv][1],
					z + latticeVelocity[inv][2],
					i);

				streamField[cellX] = collideField[cellXsubCi];
			}
		}
	}
#define USE_COPY_
#ifdef USE_COPY_
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	for (int z = 0; z < Nz; ++z)
	for (int i = 0; i < latticeDim; ++i)
	{
		const int dfIndex = funcIndex(x, y, z, i);
		collideField[dfIndex] = streamField[dfIndex];
	}
#else

	double* SW;

	SW = collideField;
	collideField = streamField;
	streamField = SW;

#endif
}

void lbSteadyFlowD3Q19::calMarcoVar()
{
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	for (int z = 0; z < Nz; ++z)
	{
		const int curCell = cellIndex(x, y, z);

		latticeUX[curCell] = 0.0;
		latticeUY[curCell] = 0.0;
		latticeUZ[curCell] = 0.0;
		latticeDensity[curCell] = 0.0;

		for (int i = 0; i < latticeDim; ++i)
		{
			const double fi = collideField[funcIndex(x, y, z, i)];

			latticeDensity[curCell] += fi;
			latticeUX[curCell] += (double)latticeVelocity[i][0] * fi;//check here
			latticeUY[curCell] += (double)latticeVelocity[i][1] * fi;
			latticeUZ[curCell] += (double)latticeVelocity[i][2] * fi;
		}

		if (flag[curCell] == INLET)
		{
			//here we define Z negative inlet velocity,  thus the z-1 near cell is fluid
			//in the future version, we will find the nearby fluid cell 
			//and set density value as near fluid cell
			latticeDensity[curCell] = latticeDensity[cellIndex(x, y, z - 1)];
			latticeUX[curCell] = latticeInletVelocity[0];
			latticeUY[curCell] = latticeInletVelocity[1];
			latticeUZ[curCell] = latticeInletVelocity[2];

		}
		else if (flag[curCell] == OUTLET)
		{
			//here we have the same set direction for outlet 
			//the outlet is at the top of the fluid cell
			//as for density we will set in DFs treatment function 
			latticeUX[curCell] = latticeUX[cellIndex(x, y, z - 1)];
			latticeUY[curCell] = latticeUY[cellIndex(x, y, z - 1)];
			latticeUZ[curCell] = latticeUZ[cellIndex(x, y, z - 1)];
		}
		else if (flag[curCell] == NO_SLIP)
		{
			latticeUX[curCell] = 0.0;
			latticeUY[curCell] = 0.0;
			latticeUZ[curCell] = 0.0;
		}

	}


}

void lbSteadyFlowD3Q19::collideDfs()
{
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	for (int z = 0; z < Nz; ++z)
	{
		const int curCell = cellIndex(x, y, z);
		if (flag[curCell] == FLUID)
		{
			//now we compuate norm stress 
			double curVelocity[3] = { latticeUX[curCell], latticeUY[curCell], latticeUZ[curCell] };

			double stress = 0.0;
			for (int alpha = 0; alpha < 3; ++alpha)
			for (int beta = 0; beta < 3; ++beta)
			{
				double elem = 0.0;
				for (int i = 0; i < latticeDim; ++i)
				{
					double feq = calFeq(latticeDensity[curCell], curVelocity, i);
					double fi = collideField[funcIndex(x, y, z, i)];

					elem += (double)latticeVelocity[i][alpha] * (double)latticeVelocity[i][beta] * (fi - feq);
				}

				stress += elem*elem;
			}

			stress = sqrt(stress);

			//now we compute relaxation time
			double tauLES =
				computeLocalRelaxationTime(relaxationTime, stress, smagorinskyConstant);

			double omegaLES = 1.0 / tauLES;

			//now we compute collide Dfs
			for (int i = 0; i < latticeDim; ++i)
			{
				const int idx = funcIndex(x, y, z, i);

				collideField[idx] = (1.0 - omegaLES)*collideField[idx] +
					omegaLES*calFeq(latticeDensity[curCell], curVelocity, i);
			}
		}
	}

}

void lbSteadyFlowD3Q19::treatDfsBoundary()
{
#pragma omp parallel for
	for (int x = 0; x < Nx; ++x)
	for (int y = 0; y < Ny; ++y)
	for (int z = 0; z < Nz; ++z)
	{
		const int curCell = cellIndex(x, y, z);

		if (flag[curCell] == INLET)
		{
			for (int i = 0; i < latticeDim; ++i)
			{
				const int Xci = x + latticeVelocity[i][0];
				const int Yci = y + latticeVelocity[i][1];
				const int Zci = z + latticeVelocity[i][2];

				const int nearIdx = cellIndex(Xci, Yci, Zci);
				double curU[3] = { latticeUX[curCell], latticeUY[curCell], latticeUZ[curCell] };
				double nearU[3] = { latticeUX[nearIdx], latticeUY[nearIdx], latticeUZ[nearIdx] };

				if (flag[nearIdx] == FLUID)
				{
					collideField[funcIndex(x, y, z, i)] =
						calFeq(latticeDensity[curCell], curU, i) +
						collideField[funcIndex(Xci, Yci, Zci, i)] - calFeq(latticeDensity[nearIdx], nearU, i);
				}

			}//for i=0-Q
		}//if inlet

		if (flag[curCell] == OUTLET)
		{
			latticeDensity[curCell] = 1.0;//here we set density equal 1.0 as atm pressure

			for (int i = 0; i < latticeDim; ++i)
			{
				const int Xci = x + latticeVelocity[i][0];
				const int Yci = y + latticeVelocity[i][1];
				const int Zci = z + latticeVelocity[i][2];

				const int nearIdx = cellIndex(Xci, Yci, Zci);
				double curU[3] = { latticeUX[curCell], latticeUY[curCell], latticeUZ[curCell] };
				double nearU[3] = { latticeUX[nearIdx], latticeUY[nearIdx], latticeUZ[nearIdx] };

				if (flag[nearIdx] == FLUID)
				{
					collideField[funcIndex(x, y, z, i)] =
						calFeq(latticeDensity[curCell], curU, i) +
						collideField[funcIndex(Xci, Yci, Zci, i)] - calFeq(latticeDensity[nearIdx], nearU, i);
				}

			}//for i=0-Q
		}


		if (flag[curCell] == NO_SLIP)
		{
			for (int i = 0; i < latticeDim; ++i)
			{
				const int Xci = x + latticeVelocity[i][0];
				const int Yci = y + latticeVelocity[i][1];
				const int Zci = z + latticeVelocity[i][2];

				if (Xci >= 0 && Xci < Nx&&Yci >= 0 && Yci < Ny&&Zci >= 0 && Zci < Nz)
				{
					const int nearIdx = cellIndex(Xci, Yci, Zci);
					double curU[3] = { 0.0, 0.0, 0.0 };
					double nearU[3] = { latticeUX[nearIdx], latticeUY[nearIdx], latticeUZ[nearIdx] };

					if (flag[nearIdx] == FLUID)
					{
						collideField[funcIndex(x, y, z, i)] =
							calFeq(latticeDensity[nearIdx], curU, i) +
							collideField[funcIndex(Xci, Yci, Zci, i)] - calFeq(latticeDensity[nearIdx], nearU, i);
					}
				}
			}//for i=0-Q
		}

	}

}

double lbSteadyFlowD3Q19::calFeq(double density, double velocity[3], int i)
{
	double cu =
		latticeVelocity[i][0] * velocity[0] +
		latticeVelocity[i][1] * velocity[1] +
		latticeVelocity[i][2] * velocity[2];

	double uu =
		velocity[0] * velocity[0] +
		velocity[1] * velocity[1] +
		velocity[2] * velocity[2];

	double feq = latticeWeight[i] * (density + 3.0*cu + 4.5*cu*cu - 1.5*uu);

	return feq;
}

double lbSteadyFlowD3Q19::computeLocalRelaxationTime(double tau_, double stressTensorNorm, double smagConstant)
{
	const double viscosity = (tau_ - 0.5) / 3.0;
	const double smagSqr = smagConstant * smagConstant;
	const double stress =
		(std::sqrt(viscosity * viscosity + 18.0 * smagSqr * stressTensorNorm) - viscosity) /
		(6.0 * smagSqr);

	return 3.0 * (viscosity + smagSqr * stress) + 0.5;
}

void lbSteadyFlowD3Q19::getVelocityField(double* velocityField)
{
#pragma omp parallel for
	for (int i = 0; i < Nx*Ny*Nz; ++i)
	{
		velocityField[i] = sqrt(
			latticeUX[i] * latticeUX[i] +
			latticeUY[i] * latticeUY[i] +
			latticeUZ[i] * latticeUZ[i]);
	}

}

int lbSteadyFlowD3Q19::getSize()
{
	return Nx*Ny*Nz;
}

lbSteadyFlowD3Q19::lbSteadyFlowD3Q19()
{
	phyInletVelocity[0] = 0.0;
	phyInletVelocity[1] = 0.0;
	phyInletVelocity[2] = -5.0;

	smagorinskyConstant = 0.12;
	fluidViscosity = 1.0e-6;
	gridLength = 0.01;

	FLUID = 0;
	INLET = 1;
	OUTLET = 2;
	NO_SLIP = 3;

	//fp = NULL;
}

lbSteadyFlowD3Q19::~lbSteadyFlowD3Q19()
{
}

const double lbSteadyFlowD3Q19::latticeWeight[19] =
{ (1. / 36), (1. / 36), (2. / 36), (1. / 36), (1. / 36),
(1. / 36), (2. / 36), (1. / 36), (2. / 36), (12. / 36),
(2. / 36), (1. / 36), (2. / 36), (1. / 36), (1. / 36),
(1. / 36), (2. / 36), (1. / 36), (1. / 36) };

const int lbSteadyFlowD3Q19::latticeVelocity[19][3] = {
	{ 0, -1, -1 }, { -1, 0, -1 }, { 0, 0, -1 }, { 1, 0, -1 }, { 0, 1, -1 },
	{ -1, -1, 0 }, { 0, -1, 0 }, { 1, -1, 0 }, { -1, 0, 0 }, { 0, 0, 0 },
	{ 1, 0, 0 }, { -1, 1, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 0, -1, 1 },
	{ -1, 0, 1 }, { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 } };
