/*
filename: TSolver.h
description: used in Gauss-Newton method
date: 2/21/2019
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
#include<Eigen/SparseCholesky>
#include <limits>

//personal requirements
#include "TDefine.hpp"
#include "TInput.h"
#include "TSolver.h"


void TSolver::m_GaussNewton()
{
	int total_time = 0;
	while (1)
	{
		double old_answer = m_Answer_R;
		m_GetAnswer_R();
		if (m_Answer_R >= old_answer) break;
		m_GetJacobi_R();
		double max = m_Solve_R();
		total_time++;
		if (max <= m_MaxDifference || total_time >= m_MaxTime)
		{
			break;
		}
		else
		{
			m_X_R += m_DX_R;
			m_DX_R.setZero();
		}
	}
	total_time = 0;
	while (1)
	{
		double old_answer = m_Answer_G;
		m_GetAnswer_G();
		if (m_Answer_G >= old_answer) break;
		m_GetJacobi_G();
		m_Solve_G();
		double max = m_Solve_G();
		total_time++;
		if (max <= m_MaxDifference || total_time >= m_MaxTime) 
		{
			break;
		}
		else
		{
			m_X_G += m_DX_G;
			m_DX_G.setZero();
		}
	}
	total_time = 0;

	while (1)
	{
		double old_answer = m_Answer_B;
		m_GetAnswer_B();
		if (m_Answer_B >= old_answer) break;
		m_GetJacobi_B();
		m_Solve_B();
		double max = m_Solve_B();
		total_time++;

		if (max <= m_MaxDifference || total_time >= m_MaxTime)
		{
			break;
		}
		else
		{
			m_X_B += m_DX_B;
			m_DX_B.setZero();
		}
	}
}

void TSolver::m_RenewColor()
{
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int group = g_PointList[i].m_GroupList[j];
			int new_color_r = m_X_R.coeff(i * g_GroupList.size() + group, 0);
			int new_color_g = m_X_G.coeff(i * g_GroupList.size() + group, 0);
			int new_color_b = m_X_B.coeff(i * g_GroupList.size() + group, 0);
			Eigen::Vector3i new_color(new_color_r, new_color_g, new_color_b);
			g_PointList[i].m_ColorList[j] += new_color;
		}
	}
}

double TSolver::m_Solve_R()
{
	Eigen::VectorXd vector_r = -1 * m_Answer_R * m_Jacobi_R;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(m_Hermit_R);
	Eigen::VectorXd temp_ans = chol.solve(vector_r);
	vector<Eigen::Triplet<double>> matrix_r;
	double max = 0;
	for (int i = 0; i < m_Size; i++)
	{
		if (temp_ans(i) != 0)
		{
			if (abs(temp_ans(i)) > max)
			{
				max = abs(temp_ans(i));
			}
			Eigen::Triplet<double> element_r(i, 0, temp_ans(i));
			matrix_r.push_back(element_r);
		}
	}
	m_DX_R.setFromTriplets(matrix_r.begin(),matrix_r.end());
	return max;
}

double TSolver::m_Solve_G()
{
	Eigen::VectorXd vector_r = -1 * m_Answer_G * m_Jacobi_G;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(m_Hermit_G);
	Eigen::VectorXd temp_ans = chol.solve(vector_r);
	vector<Eigen::Triplet<double>> matrix_r;
	double max = 0;
	for (int i = 0; i < m_Size; i++)
	{
		if (temp_ans(i) != 0)
		{
			if (abs(temp_ans(i)) > max)
			{
				max = abs(temp_ans(i));
			}
			Eigen::Triplet<double> element_r(i, 0, temp_ans(i));
			matrix_r.push_back(element_r);
		}
	}
	m_DX_G.setFromTriplets(matrix_r.begin(), matrix_r.end());
	return max;
}

double TSolver::m_Solve_B()
{
	Eigen::VectorXd vector_r = -1 * m_Answer_B * m_Jacobi_B;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(m_Hermit_B);
	Eigen::VectorXd temp_ans = chol.solve(vector_r);
	vector<Eigen::Triplet<double>> matrix_r;
	double max = 0;
	for (int i = 0; i < m_Size; i++)
	{
		if (temp_ans(i) != 0)
		{
			if (abs(temp_ans(i)) > max)
			{
				max = abs(temp_ans(i));
			}
			Eigen::Triplet<double> element_r(i, 0, temp_ans(i));
			matrix_r.push_back(element_r);
		}
	}
	m_DX_B.setFromTriplets(matrix_r.begin(), matrix_r.end());
	return max;
}


void TSolver::m_GetFX()
{
	vector<Eigen::Triplet<double>> matrix_r;
	vector<Eigen::Triplet<double>> matrix_g;
	vector<Eigen::Triplet<double>> matrix_b;

	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int num = i * g_GroupList.size() + g_PointList[i].m_GroupList[j];
			double r = g_PointList[i].m_ColorList[j](0);
			double g = g_PointList[i].m_ColorList[j](1);;
			double b = g_PointList[i].m_ColorList[j](2);;
			Eigen::Triplet<double> element_r(num, 0, r);
			Eigen::Triplet<double> element_g(num, 0, g);
			Eigen::Triplet<double> element_b(num, 0, b);
			matrix_r.push_back(element_r);
			matrix_g.push_back(element_g);
			matrix_b.push_back(element_b);

		}
	}
	m_FX_R.setFromTriplets(matrix_r.begin(), matrix_r.end());
	m_FX_G.setFromTriplets(matrix_g.begin(), matrix_g.end());
	m_FX_B.setFromTriplets(matrix_b.begin(), matrix_b.end());

}


bool JudgeAdjacent(int i1, int i2)
{
	for (int i = 0; i < g_PointList[i1].m_PointLink.size(); i++)
	{
		if (i2 == g_PointList[i1].m_PointLink[i]) return 1;
	}
	return 0;
}

void TSolver::m_GetJacobi_R()
{
	vector<Eigen::Triplet<double>> matrix_r;	
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int group = g_PointList[i].m_GroupList[j];
			double ans_r_1 = 0;
			double ans_r_2 = 0;
			double num_r = m_X_R.coeff(i * g_GroupList.size() + group, 0);
			double f_r = m_FX_R.coeff(i * g_GroupList.size() + group, 0);
			for (int k = 0; k < g_GroupList[group].m_PointList.size(); k++)
			{
				int point = g_GroupList[group].m_PointList[k];
				if (JudgeAdjacent(i, point))
				{
					ans_r_1 += num_r;
					ans_r_1 -= m_X_R.coeff(point * g_GroupList.size() + group, 0);
				}
			}
			for (int k = 0; k < g_PointList[i].m_GroupList.size(); k++)
			{
				int group_nova = g_PointList[i].m_GroupList[k];
				ans_r_2 += num_r; ans_r_2 += f_r;

				ans_r_2 -= m_X_R.coeff(i * g_GroupList.size() + group_nova, 0);
				ans_r_2 -= m_FX_R.coeff(i * g_GroupList.size() + group_nova, 0);
			}
			double ans_r = (ans_r_1 + m_Lambda * ans_r_2) * 2;
			if (ans_r != 0)
			{
				Eigen::Triplet<double> new_r(i * g_GroupList.size() + group,  0, ans_r / m_Answer_R);
				matrix_r.push_back(new_r);
			}
		}
	}
	m_Jacobi_R.setFromTriplets(matrix_r.begin(), matrix_r.end());
	m_Hermit_R = (m_Jacobi_R * m_Jacobi_R.transpose()).pruned();

	vector<Eigen::Triplet<double>> matrix_appendix;
	for (int i = 0; i < m_Size; i++)
	{
		Eigen::Triplet<double> diag(i, i, m_MaxStep);
		matrix_appendix.push_back(diag);
	}
	Eigen::SparseMatrix<double> diag_matrix(m_Size, m_Size);
	diag_matrix.setFromTriplets(matrix_appendix.begin(),matrix_appendix.end());
	m_Hermit_R += diag_matrix; //添加对角项限制步长



}

void TSolver::m_GetJacobi_G()
{
	vector<Eigen::Triplet<double>> matrix_r;
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int group = g_PointList[i].m_GroupList[j];
			double ans_r_1 = 0;
			double ans_r_2 = 0;
			double num_r = m_X_G.coeff(i * g_GroupList.size() + group, 0);
			double f_r = m_FX_G.coeff(i * g_GroupList.size() + group, 0);
			for (int k = 0; k < g_GroupList[group].m_PointList.size(); k++)
			{
				int point = g_GroupList[group].m_PointList[k];
				if (JudgeAdjacent(i, point))
				{
					ans_r_1 += num_r;
					ans_r_1 -= m_X_G.coeff(point * g_GroupList.size() + group, 0);
				}
			}
			for (int k = 0; k < g_PointList[i].m_GroupList.size(); k++)
			{
				int group_nova = g_PointList[i].m_GroupList[k];
				ans_r_2 += num_r; ans_r_2 += f_r;

				ans_r_2 -= m_X_G.coeff(i * g_GroupList.size() + group_nova, 0);
				ans_r_2 -= m_FX_G.coeff(i * g_GroupList.size() + group_nova, 0);
			}
			double ans_r = (ans_r_1 + m_Lambda * ans_r_2) * 2 ;
			if (ans_r != 0)
			{
				Eigen::Triplet<double> new_r(i * g_GroupList.size() + group, 0, ans_r / m_Answer_G);
				matrix_r.push_back(new_r);
			}
		}
	}
	m_Jacobi_G.setFromTriplets(matrix_r.begin(), matrix_r.end());
	m_Hermit_G = (m_Jacobi_G * m_Jacobi_G.transpose()).pruned();

	vector<Eigen::Triplet<double>> matrix_appendix;
	for (int i = 0; i < m_Size; i++)
	{
		Eigen::Triplet<double> diag(i, i, m_MaxStep);
		matrix_appendix.push_back(diag);
	}
	Eigen::SparseMatrix<double> diag_matrix(m_Size, m_Size);
	diag_matrix.setFromTriplets(matrix_appendix.begin(), matrix_appendix.end());
	m_Hermit_G += diag_matrix; //添加对角项限制步长
}

void TSolver::m_GetJacobi_B()
{
	vector<Eigen::Triplet<double>> matrix_r;
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int group = g_PointList[i].m_GroupList[j];
			double ans_r_1 = 0;
			double ans_r_2 = 0;
			double num_r = m_X_B.coeff(i * g_GroupList.size() + group, 0);
			double f_r = m_FX_B.coeff(i * g_GroupList.size() + group, 0);
			for (int k = 0; k < g_GroupList[group].m_PointList.size(); k++)
			{
				int point = g_GroupList[group].m_PointList[k];
				if (JudgeAdjacent(i, point))
				{
					ans_r_1 += num_r;
					ans_r_1 -= m_X_B.coeff(point * g_GroupList.size() + group, 0);
				}
			}
			for (int k = 0; k < g_PointList[i].m_GroupList.size(); k++)
			{
				int group_nova = g_PointList[i].m_GroupList[k];
				ans_r_2 += num_r; ans_r_2 += f_r;

				ans_r_2 -= m_X_B.coeff(i * g_GroupList.size() + group_nova, 0);
				ans_r_2 -= m_FX_B.coeff(i * g_GroupList.size() + group_nova, 0);
			}
			double ans_r = (ans_r_1 + m_Lambda * ans_r_2) * 2 ;
			if (ans_r != 0)
			{
				Eigen::Triplet<double> new_r(i * g_GroupList.size() + group, 0, ans_r / m_Answer_B);
				matrix_r.push_back(new_r);
			}
		}
	}
	m_Jacobi_B.setFromTriplets(matrix_r.begin(), matrix_r.end());
	m_Hermit_B = (m_Jacobi_B * m_Jacobi_B.transpose()).pruned();

	vector<Eigen::Triplet<double>> matrix_appendix;
	for (int i = 0; i < m_Size; i++)
	{
		Eigen::Triplet<double> diag(i, i, m_MaxStep);
		matrix_appendix.push_back(diag);
	}
	Eigen::SparseMatrix<double> diag_matrix(m_Size, m_Size);
	diag_matrix.setFromTriplets(matrix_appendix.begin(), matrix_appendix.end());
	m_Hermit_B += diag_matrix; //添加对角项限制步长
}




void TSolver::m_GetAnswer_R()
{
	double ans = 0;
	for (int j = 0; j < g_GroupList.size(); j++)
	{
		for (int i1 = 0; i1 < g_GroupList[j].m_PointList.size() - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < g_GroupList[j].m_PointList.size(); i2++)
			{
				int point_1 = g_GroupList[j].m_PointList[i1];
				int point_2 = g_GroupList[j].m_PointList[i2];
				if (JudgeAdjacent(point_1, point_2))
				{
					double temp_1 = m_X_R.coeff(point_1 * g_GroupList.size() + j, 0);
					double temp_2 = m_X_R.coeff(point_2 * g_GroupList.size() + j, 0);
					ans += (temp_1 - temp_2)*(temp_1 - temp_2);
				}
			}
		}
	}
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j1 = 0; j1 < g_PointList[i].m_GroupList.size() - 1; j1++)
		{
			for (int j2 = j1 + 1; j2 < g_PointList[i].m_GroupList.size(); j2++)
			{
				int group_1 = g_PointList[i].m_GroupList[j1];
				int group_2 = g_PointList[i].m_GroupList[j2];
				double temp_g_1 = m_X_R.coeff(i * g_GroupList.size() + group_1, 0);
				double temp_f_1 = m_FX_R.coeff(i * g_GroupList.size() + group_1, 0);
				double temp_g_2 = m_X_R.coeff(i * g_GroupList.size() + group_2, 0);
				double temp_f_2 = m_FX_R.coeff(i * g_GroupList.size() + group_2, 0);
				ans += m_Lambda * 
					(temp_g_1 + temp_f_1 - temp_f_2 - temp_g_2) * (temp_g_1 + temp_f_1 - temp_f_2 - temp_g_2);
			}
		}
	}
	m_Answer_R = sqrt(ans);
}

void TSolver::m_GetAnswer_G()
{
	double ans = 0;
	for (int j = 0; j < g_GroupList.size(); j++)
	{
		for (int i1 = 0; i1 < g_GroupList[j].m_PointList.size() - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < g_GroupList[j].m_PointList.size(); i2++)
			{
				int point_1 = g_GroupList[j].m_PointList[i1];
				int point_2 = g_GroupList[j].m_PointList[i2];
				if (JudgeAdjacent(point_1, point_2))
				{
					double temp_1 = m_X_G.coeff(point_1 * g_GroupList.size() + j, 0);
					double temp_2 = m_X_G.coeff(point_2 * g_GroupList.size() + j, 0);
					ans += (temp_1 - temp_2)*(temp_1 - temp_2);
				}
			}
		}
	}
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j1 = 0; j1 < g_PointList[i].m_GroupList.size() - 1; j1++)
		{
			for (int j2 = j1 + 1; j2 < g_PointList[i].m_GroupList.size(); j2++)
			{
				int group_1 = g_PointList[i].m_GroupList[j1];
				int group_2 = g_PointList[i].m_GroupList[j2];
				double temp_g_1 = m_X_G.coeff(i * g_GroupList.size() + group_1, 0);
				double temp_f_1 = m_FX_G.coeff(i * g_GroupList.size() + group_1, 0);
				double temp_g_2 = m_X_G.coeff(i * g_GroupList.size() + group_2, 0);
				double temp_f_2 = m_FX_G.coeff(i * g_GroupList.size() + group_2, 0);
				ans += m_Lambda *
					(temp_g_1 + temp_f_1 - temp_f_2 - temp_g_2) * (temp_g_1 + temp_f_1 - temp_f_2 - temp_g_2);
			}
		}
	}
	m_Answer_G = sqrt(ans);
}

void TSolver::m_GetAnswer_B()
{
	double ans = 0;
	for (int j = 0; j < g_GroupList.size(); j++)
	{
		for (int i1 = 0; i1 < g_GroupList[j].m_PointList.size() - 1; i1++)
		{
			for (int i2 = i1 + 1; i2 < g_GroupList[j].m_PointList.size(); i2++)
			{
				int point_1 = g_GroupList[j].m_PointList[i1];
				int point_2 = g_GroupList[j].m_PointList[i2];
				if (JudgeAdjacent(point_1, point_2))
				{
					double temp_1 = m_X_B.coeff(point_1 * g_GroupList.size() + j, 0);
					double temp_2 = m_X_B.coeff(point_2 * g_GroupList.size() + j, 0);
					ans += (temp_1 - temp_2)*(temp_1 - temp_2);
				}
			}
		}
	}
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j1 = 0; j1 < g_PointList[i].m_GroupList.size() - 1; j1++)
		{
			for (int j2 = j1 + 1; j2 < g_PointList[i].m_GroupList.size(); j2++)
			{
				int group_1 = g_PointList[i].m_GroupList[j1];
				int group_2 = g_PointList[i].m_GroupList[j2];
				double temp_g_1 = m_X_B.coeff(i * g_GroupList.size() + group_1, 0);
				double temp_f_1 = m_FX_B.coeff(i * g_GroupList.size() + group_1, 0);
				double temp_g_2 = m_X_B.coeff(i * g_GroupList.size() + group_2, 0);
				double temp_f_2 = m_FX_B.coeff(i * g_GroupList.size() + group_2, 0);
				ans += m_Lambda *
					(temp_g_1 + temp_f_1 - temp_f_2 - temp_g_2) * (temp_g_1 + temp_f_1 - temp_f_2 - temp_g_2);
			}
		}
	}
	m_Answer_B = sqrt(ans);
}