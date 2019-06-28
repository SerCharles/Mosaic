/*
filename: TSolver.cu
description: used in Gauss-Newton method
date: 4/25/2019
*/


//system requirements
#include <time.h>
#include <iostream>
#include <fstream>
#include <windows.h>
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
#include "TSolver.cuh"

//cuda
#include <cusolverSp.h>
#include <cusparse_v2.h>
#include <cublas_v2.h>
#include <cuda_runtime_api.h>

void TSolver::m_GaussNewton()
{
	int total_time = 0;
	m_GetFX_R();
	while (1)
	{
		double old_answer = m_Answer_R;
		m_GetAnswer_R();
		if (m_Answer_R >= old_answer) break;
		m_GetJacobi_R();
		m_PrepareSolve_R();
		//m_Test();
		double max = m_Solve_R();
		total_time++;
		if (max <= m_MaxDifference || total_time >= m_MaxTime)
		{
			break;
		}
		else
		{
			for (int i = 0; i < m_Size; i++)
			{
				m_X_R[i] += m_DX_R[i];
			}
		}
	}
	total_time = 0;
}

void TSolver::m_RenewColor()
{
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int group = g_PointList[i].m_GroupList[j];
			int new_color_r = m_X_R[i * g_GroupList.size() + group];

			///todo: g and b
			//int new_color_g = m_X_G.coeff(i * g_GroupList.size() + group, 0);
			//int new_color_b = m_X_B.coeff(i * g_GroupList.size() + group, 0);
			Eigen::Vector3i new_color(new_color_r, 0, 0);
			g_PointList[i].m_ColorList[j] += new_color;
		}
	}
}

bool JudgeAdjacent(int i1, int i2)
{
	for (int i = 0; i < g_PointList[i1].m_PointLink.size(); i++)
	{
		if (i2 == g_PointList[i1].m_PointLink[i]) return 1;
	}
	return 0;
}

void TSolver::m_GetFX_R()
{
	memset(m_FX_R, 0, m_Size * sizeof(int));

	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int num = i * g_GroupList.size() + g_PointList[i].m_GroupList[j];
			double r = g_PointList[i].m_ColorList[j](0);
			m_FX_R[num] = r;
		}
	}
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
					double temp_1 = m_X_R[point_1 * g_GroupList.size() + j];
					double temp_2 = m_X_R[point_2 * g_GroupList.size() + j];
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
				double temp_g_1 = m_X_R[i * g_GroupList.size() + group_1];
				double temp_f_1 = m_FX_R[i * g_GroupList.size() + group_1];
				double temp_g_2 = m_X_R[i * g_GroupList.size() + group_2];
				double temp_f_2 = m_FX_R[i * g_GroupList.size() + group_2];
				ans += m_Lambda *
					(temp_g_1 + temp_f_1 - temp_f_2 - temp_g_2) * (temp_g_1 + temp_f_1 - temp_f_2 - temp_g_2);
			}
		}
	}
	m_Answer_R = sqrt(ans);
}

void TSolver::m_GetJacobi_R()
{
	m_NonZeroJacobi_R.clear();
	//m_JacobiNorm = 0;
	memset(m_Jacobi_R, 0, m_Size * sizeof(int));
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int group = g_PointList[i].m_GroupList[j];
			double ans_r_1 = 0;
			double ans_r_2 = 0;
			double num_r = m_X_R[i * g_GroupList.size() + group];
			double f_r = m_FX_R[i * g_GroupList.size() + group];
			for (int k = 0; k < g_GroupList[group].m_PointList.size(); k++)
			{
				int point = g_GroupList[group].m_PointList[k];
				if (JudgeAdjacent(i, point))
				{
					ans_r_1 += num_r;
					ans_r_1 -= m_X_R[point * g_GroupList.size() + group];
				}
			}
			for (int k = 0; k < g_PointList[i].m_GroupList.size(); k++)
			{
				int group_nova = g_PointList[i].m_GroupList[k];
				ans_r_2 += num_r; ans_r_2 += f_r;

				ans_r_2 -= m_X_R[i * g_GroupList.size() + group_nova];
				ans_r_2 -= m_FX_R[i * g_GroupList.size() + group_nova];
			}
			double ans_r = (ans_r_1 + m_Lambda * ans_r_2) * 2;

			m_Jacobi_R[i * g_GroupList.size() + group] = ans_r / m_Answer_R;
		}
	}
	for (int i = 0; i < m_Size; i++)
	{
		if (m_Jacobi_R[i] != 0)
		{
			m_NonZeroJacobi_R.push_back(i);
		}
	}



}

void TSolver::m_PrepareSolve_R()
{
	//get b
	for (int i = 0; i < m_Size; i++)
	{
		m_B_R[i] = - m_Jacobi_R[i] * m_Answer_R;
	}



	//get H(x)
	//Actually it is Q(x)=H(x)*J(x).norm2
	m_Hermit_R.m_Renew();
	m_Hermit_R.m_NonZeros = m_NonZeroJacobi_R.size() *  m_NonZeroJacobi_R.size() + 
		m_Size - m_NonZeroJacobi_R.size();

	m_Hermit_R.m_Place_x = new int[m_Hermit_R.m_NonZeros];
	m_Hermit_R.m_Place_y = new int[m_Hermit_R.m_NonZeros];
	m_Hermit_R.m_Data = new double[m_Hermit_R.m_NonZeros];
	int flag = 0;
	for (int i = 0; i < m_NonZeroJacobi_R.size(); i++)
	{
		int start = 0, end = 0;
		if (i == 0)
		{
			start = 0;
			end = m_NonZeroJacobi_R[0];
		}
		else
		{
			start = m_NonZeroJacobi_R[i - 1] + 1;
			end = m_NonZeroJacobi_R[i];
		}
		for (int j = start; j < end; j++)
		{
			m_Hermit_R.m_Place_x[flag] = j;
			m_Hermit_R.m_Place_y[flag] = j;
			m_Hermit_R.m_Data[flag] = m_HermitDiag;
			flag++;
		}
		for (int jj = 0; jj < m_NonZeroJacobi_R.size(); jj++)
		{
			m_Hermit_R.m_Place_x[flag] = m_NonZeroJacobi_R[i];
			m_Hermit_R.m_Place_y[flag] = m_NonZeroJacobi_R[jj];
			double data = m_Jacobi_R[m_NonZeroJacobi_R[i]] * m_Jacobi_R[m_NonZeroJacobi_R[jj]];
			if (i == jj)
			{
				//添加对角项限制步长，保证正定
				data += m_HermitDiag;
			}
			m_Hermit_R.m_Data[flag] = data;
			flag++;
		}
	}
	for (int j = m_NonZeroJacobi_R[m_NonZeroJacobi_R.size() - 1] + 1; j < m_Size; j++)
	{
		m_Hermit_R.m_Place_x[flag] = j;
		m_Hermit_R.m_Place_y[flag] = j;
		m_Hermit_R.m_Data[flag] = m_HermitDiag;
		flag++;
	}


}

void TSolver::m_Test()
{
	ofstream file_1("hermit.txt");
	file_1 << m_Hermit_R.m_NonZeros << endl;
	for (int i = 0; i <= m_Size; i++) file_1 << m_Hermit_R.m_Place_x[i]<<" ";
	file_1 << endl;
	for (int i = 0; i < m_Hermit_R.m_NonZeros; i++) file_1 << m_Hermit_R.m_Place_y[i]<<" ";
	file_1 << endl;
	for(int i = 0; i < m_Hermit_R.m_NonZeros;i++) file_1 << m_Hermit_R.m_Data[i] << " ";
	file_1 << endl;
	ofstream file_2("b.txt");
	for (int i = 0; i < m_Size; i++) file_2 << m_B_R[i] << " ";
	ofstream file_3("jacobi.txt");
	for (int i = 0; i < m_Size; i++) file_3 << m_Jacobi_R[i] << " ";

}

double TSolver::m_Solve_R()
{
	memset(m_DX_R, 0, m_Size * sizeof(int));

	//define status and variables
	cusparseStatus_t status;
	cublasStatus_t status_cub;
	cusparseHandle_t handle = 0;
	cublasHandle_t cublasH = NULL;
	cudaStream_t stream = NULL;
	cusparseMatDescr_t descr = 0;


	cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
	cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t cuda_state1 = cudaSuccess;
	cudaError_t cuda_state2 = cudaSuccess;
	cudaError_t cuda_state3 = cudaSuccess;
	cudaError_t cuda_state4 = cudaSuccess;
	cudaError_t cuda_state5 = cudaSuccess;
	int *csr_gpu_x = NULL;
	int *csr_gpu_offset = NULL;
	int *csr_gpu_y = NULL;
	double *csr_gpu_value = NULL;
	double *csr_gpu_b = NULL; 
	double *csr_gpu_r = NULL;
	double *csr_gpu_d = NULL;
	double *csr_gpu_dx = NULL; 


	//allocate gpu space
	cuda_state1 = cudaMalloc((void**)&csr_gpu_value, sizeof(double) * m_Hermit_R.m_NonZeros);
	cuda_state2 = cudaMalloc((void**)&csr_gpu_y, sizeof(int) * m_Hermit_R.m_NonZeros);
	cuda_state3 = cudaMalloc((void**)&csr_gpu_x, sizeof(int) * m_Hermit_R.m_NonZeros);
	cuda_state4 = cudaMalloc((void**)&csr_gpu_b, sizeof(double) * m_Size);
	cuda_state5 = cudaMalloc((void**)&csr_gpu_dx, sizeof(double) * m_Size);
	assert(cuda_state1 == cudaSuccess);
	assert(cuda_state2 == cudaSuccess);
	assert(cuda_state3 == cudaSuccess);
	assert(cuda_state4 == cudaSuccess);
	assert(cuda_state5 == cudaSuccess);

	cuda_state1 = cudaMalloc((void**)&csr_gpu_offset, sizeof(int) * (m_Size + 1));
	assert(cuda_state1 == cudaSuccess);


	cuda_state1 = cudaMalloc((void**)&csr_gpu_d, sizeof(double) * m_Size);
	cuda_state2 = cudaMalloc((void**)&csr_gpu_r, sizeof(double) * m_Size);
	assert(cuda_state1 == cudaSuccess);
	assert(cuda_state2 == cudaSuccess);

	//copy cpu data to gpu
	cuda_state1 = cudaMemcpy(csr_gpu_value, m_Hermit_R.m_Data, sizeof(double) * m_Hermit_R.m_NonZeros, cudaMemcpyHostToDevice);
	cuda_state2 = cudaMemcpy(csr_gpu_y, m_Hermit_R.m_Place_y, sizeof(int) * m_Hermit_R.m_NonZeros, cudaMemcpyHostToDevice);
	cuda_state3 = cudaMemcpy(csr_gpu_x, m_Hermit_R.m_Place_x, sizeof(int) * m_Hermit_R.m_NonZeros, cudaMemcpyHostToDevice);
	cuda_state4 = cudaMemcpy(csr_gpu_b, m_B_R, sizeof(double) * m_Size, cudaMemcpyHostToDevice);
	assert(cuda_state1 == cudaSuccess);
	assert(cuda_state2 == cudaSuccess);
	assert(cuda_state3 == cudaSuccess);
	assert(cuda_state4 == cudaSuccess);



	//init handles
	/* initialize cusparse library */
	status = cusparseCreate(&handle);
	assert(status == CUSPARSE_STATUS_SUCCESS);

	/* create and setup matrix descriptor */
	status = cusparseCreateMatDescr(&descr);
	assert(status == CUSPARSE_STATUS_SUCCESS);

	/*cublas*/
	status_cub = cublasCreate(&cublasH);
	assert(CUBLAS_STATUS_SUCCESS == status_cub);

	status_cub = cublasSetStream(cublasH, stream);
	assert(CUBLAS_STATUS_SUCCESS == status_cub);

	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	//init csr version H
	status = cusparseXcoo2csr(handle, csr_gpu_x,
		m_Hermit_R.m_NonZeros, m_Size, csr_gpu_offset, CUSPARSE_INDEX_BASE_ZERO);
	assert(status == CUSPARSE_STATUS_SUCCESS);


	//init dx=0,d=b,r=b
	cuda_state1 = cudaMemset(csr_gpu_dx, 0, m_Size * sizeof(double));
	assert(cuda_state1 == cudaSuccess);

	cuda_state1 = cudaMemcpy(csr_gpu_d, csr_gpu_b, sizeof(double) * m_Size, cudaMemcpyDeviceToDevice);
	cuda_state2 = cudaMemcpy(csr_gpu_r, csr_gpu_b, sizeof(double) * m_Size, cudaMemcpyDeviceToDevice);
	assert(cuda_state1 == cudaSuccess);
	assert(cuda_state2 == cudaSuccess);

	//iteration
	for (int i = 1; i <= m_Size; i++)
	{

		//calculate ak
		double alpha = 0;
		//get alpha k
		double rk_norm = 0;//分子
		double dk_norm = 0;//分母
		status_cub = cublasDdot(cublasH, m_Size,
			csr_gpu_r, 1,
			csr_gpu_r, 1,
			&rk_norm);
			assert(CUBLAS_STATUS_SUCCESS == status_cub);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);


		double* temp = NULL; //temp = Qd
		cuda_state1 = cudaMalloc((void**)&temp, sizeof(double) * m_Size);;
		assert(cuda_state1 == cudaSuccess);

		cuda_state1 = cudaMemset(temp, 0, sizeof(double) * m_Size);
		assert(cuda_state1 == cudaSuccess);
		double temp_alpha = 1, temp_beta = 0;


		status = cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
			m_Size, m_Size, m_Hermit_R.m_NonZeros, &temp_alpha,
			descr,
			csr_gpu_value,
			csr_gpu_offset, csr_gpu_y,
			csr_gpu_d, &temp_beta,
			temp);
		assert(status == CUSPARSE_STATUS_SUCCESS);

		//test
		/*double temp_norm = 0;
		status_cub = cublasDnrm2(cublasH, m_Size, temp, 1, &temp_norm);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);*/



		status_cub = cublasDdot(cublasH, m_Size,
			csr_gpu_d, 1,
			temp, 1,
			&dk_norm);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);


		alpha = rk_norm / dk_norm;


		//calculate dxk
		status_cub = cublasDaxpy(cublasH, m_Size,
			&alpha,
			csr_gpu_d, 1,
			csr_gpu_dx, 1);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);
		

		//test
		/*double dx_norm = 0;
		status_cub = cublasDnrm2(cublasH, m_Size, csr_gpu_dx, 1, &dx_norm);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);*/


		//test
		/*double* old_r = NULL;
		cudaMalloc((void**)&old_r, m_Size * sizeof(double));
		cudaMemcpy(old_r, csr_gpu_r, m_Size * sizeof(double), cudaMemcpyDeviceToDevice);
		*/

		//calculate rk
		double minus_alpha = 0 - alpha;
		status_cub = cublasDaxpy(cublasH, m_Size,
			&minus_alpha,
			temp, 1,
			csr_gpu_r, 1);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);
		
		//test
		/*double kk = 0;
		status_cub = cublasDdot(cublasH, m_Size,
			csr_gpu_r, 1,
			old_r, 1,
			&kk);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);
		*/

		double rk_plus1_norm = 0;
		status_cub = cublasDdot(cublasH, m_Size,
			csr_gpu_r, 1,
			csr_gpu_r, 1,
			&rk_plus1_norm);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);

		if (rk_plus1_norm <= m_MaxInSolve)
		{
			cudaFree(temp);
			break;
		}

		//calculate betak
		double beta = rk_plus1_norm / rk_norm;

		//calculate dk+1
		double* new_temp = NULL;
		cuda_state1 = cudaMalloc((void**)&new_temp, sizeof(double) * m_Size);;
		assert(cuda_state1 == cudaSuccess);

		cuda_state1 = cudaMemset(new_temp,0,sizeof(double) * m_Size);
		assert(cuda_state1 == cudaSuccess);


		status_cub = cublasDaxpy(cublasH, m_Size,
			&beta,
			csr_gpu_d, 1,
			new_temp, 1);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);

		double one = 1;
		status_cub = cublasDaxpy(cublasH, m_Size,
			&one,
			csr_gpu_r, 1,
			new_temp, 1);
		assert(CUBLAS_STATUS_SUCCESS == status_cub);
		cudaMemcpy(csr_gpu_d, new_temp, m_Size * sizeof(double), cudaMemcpyDeviceToDevice);
		cudaFree(temp);
		cudaFree(new_temp);
		//cudaFree(old_r);
	}




	//get result and end
	double max = 0;
	status_cub = cublasDnrm2(cublasH, m_Size, csr_gpu_dx, 1, &max);
	assert(CUBLAS_STATUS_SUCCESS == status_cub);

	cudaMemcpy(m_DX_R, csr_gpu_dx, m_Size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(csr_gpu_value);
	cudaFree(csr_gpu_x);
	cudaFree(csr_gpu_y);
	cudaFree(csr_gpu_b);
	cudaFree(csr_gpu_offset);
	cudaFree(csr_gpu_dx);
	cudaFree(csr_gpu_r);
	cudaFree(csr_gpu_d);

	return max;
}