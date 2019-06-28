/*
filename: TDefine.hpp
description: the hpp that defines points, meshes
date: 2/15/2019
*/



#ifndef TDEFINE_HPP
#define TDEFINE_HPP

#pragma warning(disable:4819)
#pragma warning(disable:4244)
#pragma warning(disable:4267)

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
#include <limits>

#define SIZE 110000

class TPoint
{
public:
	Eigen::Vector3f m_Place;
	Eigen::Vector3i m_RGB;	//0:R 1:G 2:B
	std::vector<int> m_SonMesh;
	std::vector<int> m_PointLink;
	std::vector<int> m_GroupList;
	std::vector<Eigen::Vector3i> m_ColorList;//commensuate with grouplist
	int m_Number;

	TPoint(double x, double y, double z, int num)
	{
		m_Place << x, y, z;
		m_Number = num;
	}

	void m_AddSon(int son_num)
	{
		m_SonMesh.push_back(son_num);
	}
};


class TMesh
{
public:
	Eigen::Vector3i m_Vertex;
	int m_Mosaic;			// mosaic
	Eigen::Vector3i m_RGB;	//0:R 1:G 2:B
	int m_Belong;			// the group in which the point belongs to
	int m_Number;
	Eigen::Vector3f m_Normal;
	std::vector<int> m_MeshLink;

	TMesh(int a, int b, int c, int num)
	{
		m_Vertex << a, b, c;
		m_Number = num;
		m_RGB.fill(0);
		m_Belong = -1;
		m_Mosaic = num % 20;
	}

	
};


class TGroup
{
public:
	int m_Number, m_Mosaic;
	std::vector<int> m_MeshList;
	std::vector<int> m_PointList;
	TGroup(int num, int mosaic)
	{
		m_Number = num;
		m_Mosaic = mosaic;
	}

	void m_AddMesh(int num)
	{
		m_MeshList.push_back(num);
	}

};

class TPicture
{
public:
	Eigen::Matrix<float, 3, 3> m_IntMatrix;
	Eigen::Matrix<float, 3, 4> m_ExtMatrix;
	Eigen::Matrix<float, 3, 3> m_SpinningMatrix;
	Eigen::Matrix<float, 3, 1> m_MovingMatrix;
	Eigen::Matrix<float, 3, 4> m_CameraMatrix;
	Eigen::Vector3f m_UpDirection;
	Eigen::Vector3f m_SeeingDirection;
	Eigen::Vector3f m_View;
	cv::Mat m_Image;
	int m_Number;
	TPicture(int num)
	{
		m_Number = num;
	}
};

#endif