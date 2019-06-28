/*
filename: TSolver.cuh
description: used in Gauss-Newton method
date: 4/25/2019
*/


//system requirements
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

//cv requirements
#include <opencv2/opencv.hpp>

//algebra requirements
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <limits>

//personal requirements
#include "TDefine.hpp"
#include "TInput.h"

//cuda
#include <cusolverSp.h>
#include <cuda_runtime_api.h>

#define MAX 1145141919810

#ifndef TSOLVER_CUH
#define TSOLVER_CUH

struct TSparseMatrixCOO
{
	int m_NonZeros{ 0 };
	int* m_Place_x{ NULL }; //index
	int* m_Place_y{ NULL }; //index
	double* m_Data{ NULL }; //data
	~TSparseMatrixCOO()
	{
		if (m_Place_x != NULL) delete(m_Place_x);
		if (m_Place_y != NULL) delete(m_Place_y);
		if (m_Data != NULL) delete(m_Data);
	}
	void m_Renew()
	{
		if (m_Place_x != NULL) delete(m_Place_x);
		if (m_Place_y != NULL) delete(m_Place_y);
		if (m_Data != NULL) delete(m_Data);
	}
};



class TSolver
{
public:
	const int m_Size;
	const int m_Lambda{ 100 };    //a constant used in calculating the answer
	const double m_MaxTime{ 5 };  //max iteration time
	const double m_HermitDiag{ 5 }; // the value needed to add to the diag of the Hermit
	const double m_MaxDifference; //min reduction of value per time 
	const double m_MaxInSolve; //if r < this, Ax=b is solved

	//hermit of x
	TSparseMatrixCOO m_Hermit_R;
	double* m_B_R;
	double* m_X_R;
	double* m_FX_R;
	double* m_DX_R;
	double* m_Jacobi_R;
	double m_Answer_R{ MAX };
	//double m_JacobiNorm;
	vector<int> m_NonZeroJacobi_R;
	//todo:G,B

public:
	TSolver() :m_Size(g_GroupList.size() * g_PointList.size()), m_MaxDifference(m_Size / 200), m_MaxInSolve(m_Size / 1000)
	{
		m_B_R = new double[m_Size];
		m_X_R = new double[m_Size];
		for (int i = 0; i < m_Size; i++)
		{
			m_X_R[i] = 0;
		}
		m_FX_R = new double[m_Size];
		m_DX_R = new double[m_Size];
		m_Jacobi_R = new double[m_Size];
		
	}
	~TSolver()
	{
		delete(m_X_R);
		delete(m_FX_R);
		delete(m_DX_R);
		delete(m_Jacobi_R);
		m_NonZeroJacobi_R.clear();

	}

	TSolver(int i) :m_Size(2), m_MaxDifference(1), m_MaxInSolve(0)
	{
		m_B_R = new double[m_Size];
		m_X_R = new double[m_Size];
		for (int i = 0; i < m_Size; i++)
		{
			m_X_R[i] = 0;
		}
		m_FX_R = new double[m_Size];
		m_DX_R = new double[m_Size];
		m_Jacobi_R = new double[m_Size];
		m_B_R[0] = 1; m_B_R[1] = 0;
		m_Hermit_R.m_NonZeros = 4;
		m_Hermit_R.m_Place_x = new int[4];
		m_Hermit_R.m_Place_x[0] = 0;m_Hermit_R.m_Place_x[1] = 0;
		m_Hermit_R.m_Place_x[2] = 1;m_Hermit_R.m_Place_x[3] = 1;
		m_Hermit_R.m_Place_y = new int[4];
		m_Hermit_R.m_Place_y[0] = 0;m_Hermit_R.m_Place_y[1] = 1;
		m_Hermit_R.m_Place_y[2] = 0;m_Hermit_R.m_Place_y[3] = 1;
		m_Hermit_R.m_Data = new double[4];
		m_Hermit_R.m_Data[0] = 2; m_Hermit_R.m_Data[1] = -1;
		m_Hermit_R.m_Data[2] = -1; m_Hermit_R.m_Data[3] = 2;
		m_Solve_R();
	}

	//main turn
	void m_GaussNewton();

	//renew the colors after solving
	void m_RenewColor();

	//get F(X)
	void m_GetFX_R();

	//get J(X)
	void m_GetJacobi_R();

	//get total value
	void m_GetAnswer_R();

	//prepare to solve the equation, get H(X) and b
	void m_PrepareSolve_R();

	//solve H(X)dX = b, return the max difference of dx
	double m_Solve_R();

	void m_Test();
};

#endif