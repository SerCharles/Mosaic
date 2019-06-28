/*
filename: TSolver.h
description: used in Gauss-Newton method
date: 2/21/2019
*/

#ifndef TSOLVER_H
#define TSOLVER_H





//system requirements
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

//cv requirements
#include <opencv2/opencv.hpp>

//algebra requirements
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <limits>

//personal requirements
#include "TDefine.hpp"
#include "TInput.h"

#define MAX 1145141919810

class TSolver
{
public:
	const int m_Size;
	const int m_Lambda{ 100 };
	const double m_MaxTime{ 5 };
	const double m_MaxStep{ 5 };
	const double m_MaxDifference;

	Eigen::SparseMatrix<double> m_Hermit_R;
	Eigen::SparseMatrix<double> m_X_R;
	Eigen::SparseMatrix<double> m_FX_R;
	Eigen::SparseMatrix<double> m_DX_R;
	Eigen::SparseMatrix<double> m_Jacobi_R;
	double m_Answer_R{ MAX };

	Eigen::SparseMatrix<double> m_Hermit_G;
	Eigen::SparseMatrix<double> m_X_G;
	Eigen::SparseMatrix<double> m_FX_G;
	Eigen::SparseMatrix<double> m_DX_G;
	Eigen::SparseMatrix<double> m_Jacobi_G;
	double m_Answer_G{ MAX };

	Eigen::SparseMatrix<double> m_Hermit_B;
	Eigen::SparseMatrix<double> m_X_B;
	Eigen::SparseMatrix<double> m_FX_B;
	Eigen::SparseMatrix<double> m_DX_B;
	Eigen::SparseMatrix<double> m_Jacobi_B;
	double m_Answer_B{ MAX };
public:
	TSolver() :m_Size(g_GroupList.size() * g_PointList.size()), m_MaxDifference(m_Size / 1000);
	{
		m_Hermit_R.resize(m_Size, m_Size);
		m_Hermit_G.resize(m_Size, m_Size);
		m_Hermit_B.resize(m_Size, m_Size);
		m_X_R.resize(m_Size,1);
		m_X_G.resize(m_Size, 1);
		m_X_B.resize(m_Size, 1);
		m_DX_R.resize(m_Size, 1);
		m_DX_G.resize(m_Size, 1);
		m_DX_B.resize(m_Size, 1);
		m_FX_R.resize(m_Size, 1);
		m_FX_G.resize(m_Size, 1);
		m_FX_B.resize(m_Size, 1);
		m_Jacobi_R.resize(m_Size, 1);
		m_Jacobi_G.resize(m_Size, 1);
		m_Jacobi_B.resize(m_Size, 1);
		m_GetFX();
	}
	void m_GaussNewton();
	void m_GetAnswer_R();
	void m_GetAnswer_G();
	void m_GetAnswer_B();

	void m_GetFX();
	void m_GetJacobi_R();
	void m_GetJacobi_G();
	void m_GetJacobi_B();


	//return: min value in the answer
	double m_Solve_R();
	double m_Solve_G();
	double m_Solve_B();

	void m_RenewColor();
};
#endif
