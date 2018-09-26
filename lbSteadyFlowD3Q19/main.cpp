#include "stdafx.h"
#include "lbSteadyFlowD3Q19.h"

int _tmain(int argc, _TCHAR* argv[])
{
	lbSteadyFlowD3Q19 A;

	//char* sgnFileName = "castingWithRiser.sgn";
	char* sgnFileName = "d3mm.sgn";

	char* outputFileName = "D:\\output\\results.plt";

	A.preReadSgn(sgnFileName);

	A.setSgnMem();

	A.readSgn(sgnFileName);

	A.setMem();

	A.initSimuParam();

	A.initSimuVarValue();

	A.initDfsValue();

	const int gridSize = A.getSize();
	double* curU = new double[gridSize];
	double* nexU = new double[gridSize];

	int MaxIter=100;
	const int logIter = 5;
	const int writeIter = 50;

	cout << "input max iter num:\n";
	cin >> MaxIter;

	double velocityAccuracy = 1.0e-4;
	cout << "input velocity accuracy:\n";
	cin >> velocityAccuracy;
	
	double t_start = omp_get_wtime();

	for (int iter = 0; iter <= MaxIter; ++iter)
	{
		A.getVelocityField(curU);

		A.streamDfs();

		A.calMarcoVar();

		A.collideDfs();

		A.treatDfsBoundary();

		A.getVelocityField(nexU);

		double u_err = 0.0;
		double u_sum = 0.0;
		for (int i = 0; i < gridSize; ++i)
		{
			u_err += fabs(curU[i] - nexU[i]);
			u_sum += curU[i];
		}

		u_err /= u_sum;
		if (u_sum == 0.0)
			u_err= 1.0;

		if (iter%logIter == 0)
		{
			printf("iter=%7d U err=%e\n", iter,u_err);
		}

		if (iter % writeIter == 0)
		{
			char chNameU[600] = "D:\\output\\steadyFlow";
			char chTailerU[255] = ".plt";
			char str_tmpU[255];
			sprintf(str_tmpU, "%d", 10000+iter/writeIter);
			strcat(chNameU, str_tmpU);
			strcat(chNameU, chTailerU);

			A.writeData(chNameU);
		}
		
		if (u_err < velocityAccuracy)
		{
			printf("calculation finished:\niter=%7d U err=%e\n", iter, u_err);
			break;
		}
	}
	
	double t_end = omp_get_wtime();

	printf("solver time used: %6.2f s\n", t_end - t_start);

	A.writeData(outputFileName);

	A.resetMem();

	return 0;
}