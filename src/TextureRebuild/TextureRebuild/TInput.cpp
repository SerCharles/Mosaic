/*
filename: TInput.cpp
description: used in inputing .ply files (for those with only triangle meshes), matrix files and pictures
also includes .ply output
date: 2/15/2019
*/

//system requirements
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

//cv requirements
#include <opencv2/opencv.hpp>

//algebra requirements
#include <Eigen/Eigen>
#include <limits>

//personal requirements
#include "TDefine.hpp"
#include "TInput.h"
#include "TMath.h"

#define SIZE 110000

vector<TPoint> g_PointList;
vector<TMesh> g_MeshList;
vector<TPicture> g_PictureList;
vector<TGroup> g_GroupList;
Eigen::Vector3f g_Center;


bool JudgeRepeat(int new_num, int num_in_list)
{
	for (int i = 0; i < g_PointList[num_in_list].m_PointLink.size(); i++)
	{
		if (new_num == g_PointList[num_in_list].m_PointLink[i]) return 1;
	}
	return 0;
}


void g_PlyInput(char filename[], std::vector<TPoint>* point_list, std::vector<TMesh>* mesh_list)
{
	int num_point = 0, num_mesh = 0;

	ifstream file(filename);

	g_Center(0) = 0; g_Center(1) = 0; g_Center(2) = 0;


	for (int i = 1; i <= 3; i++)
	{
		char cc[1000];
		file.getline(cc, 1000);
	}
	char c1[100];
	file.getline(c1, 100);
	for (int i = 0; i < strlen(c1); i++)
	{
		if (c1[i] >= '0' && c1[i] <= '9')
		{
			num_point *= 10;
			num_point += (c1[i] - '0');
		}
	}

	for (int i = 1; i <= 3; i++)
	{
		char cc[1000];
		file.getline(cc, 1000);
	}

	char c2[100];
	file.getline(c2, 100);
	for (int i = 0; i < strlen(c2); i++)
	{
		if (c2[i] >= '0' && c2[i] <= '9')
		{
			num_mesh *= 10;
			num_mesh += (c2[i] - '0');
		}
	}

	for (int i = 1; i <= 2; i++)
	{
		char cc[1000];
		file.getline(cc, 1000);
	}

	for (int i = 0; i < num_point; i++)
	{
		double x = 0, y = 0, z = 0;
		file >> x >> y >> z;
		TPoint new_point(x, y, z, i);
		point_list->push_back(new_point);
		g_Center(0) += x;
		g_Center(1) += y;
		g_Center(2) += z;
	}
	g_Center /= num_point;

	for (int i = 0; i < num_mesh; i++)
	{
		int k = 0, a = 0, b = 0, c = 0;
		file >> k >> a >> b >> c;
		TMesh new_mesh(a, b, c, i);
		new_mesh.m_Normal = g_CalculateNormal((*point_list)[a].m_Place, (*point_list)[b].m_Place, (*point_list)[c].m_Place);
		mesh_list->push_back(new_mesh);
		(*point_list)[a].m_AddSon(i);
		(*point_list)[b].m_AddSon(i);
		(*point_list)[c].m_AddSon(i);
		
	}

	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_SonMesh.size(); j++)
		{
			int first_mesh_num = g_PointList[i].m_SonMesh[j];
			for (int k = 0; k < g_PointList[i].m_SonMesh.size(); k++)
			{
				if (k == j)
				{
					continue;
				}
				int second_mesh_num = g_PointList[i].m_SonMesh[k];
				int start_node = -1, end_node = -1;
				g_GetCommonPoints(start_node, end_node, g_MeshList[first_mesh_num], g_MeshList[second_mesh_num]);
				if (start_node != -1 && end_node != -1)
				{
					bool whether_good_1 = 1, whether_good_2 = 1;
					for (int ii = 0; ii < g_MeshList[first_mesh_num].m_MeshLink.size(); ii++)
					{
						if (second_mesh_num == g_MeshList[first_mesh_num].m_MeshLink[ii])
						{
							whether_good_1 = 0;
							break;
						}
					}
					for (int ii = 0; ii < g_MeshList[second_mesh_num].m_MeshLink.size(); ii++)
					{
						if (first_mesh_num == g_MeshList[second_mesh_num].m_MeshLink[ii])
						{
							whether_good_2 = 0;
							break;
						}
					}
					if(whether_good_1 == 1) g_MeshList[first_mesh_num].m_MeshLink.push_back(second_mesh_num);
					if (whether_good_2 == 1) g_MeshList[second_mesh_num].m_MeshLink.push_back(first_mesh_num);
				}

			}
		}
	}

	for (int i = 0; i < g_MeshList.size(); i++)
	{
		int point_0 = g_MeshList[i].m_Vertex(0);
		int point_1 = g_MeshList[i].m_Vertex(1);
		int point_2 = g_MeshList[i].m_Vertex(2);
		if (!JudgeRepeat(point_1, point_0)) g_PointList[point_0].m_PointLink.push_back(point_1);
		if (!JudgeRepeat(point_2, point_0)) g_PointList[point_0].m_PointLink.push_back(point_2);
		if (!JudgeRepeat(point_0, point_1)) g_PointList[point_1].m_PointLink.push_back(point_0);
		if (!JudgeRepeat(point_2, point_1)) g_PointList[point_1].m_PointLink.push_back(point_2);
		if (!JudgeRepeat(point_0, point_2)) g_PointList[point_2].m_PointLink.push_back(point_0);
		if (!JudgeRepeat(point_1, point_2)) g_PointList[point_2].m_PointLink.push_back(point_1);
	}


}

void g_MatrixInput(char filename[], std::vector<TPicture>* picture_list)
{
	std::ifstream fin(filename);
	int num;
	while (fin >> num)
	{
		TPicture new_picture(num);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				fin >> new_picture.m_IntMatrix(i, j);

		double temp;
		fin >> temp;
		fin >> temp;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 4; j++)
				fin >> new_picture.m_ExtMatrix(i, j);

		new_picture.m_CameraMatrix = new_picture.m_IntMatrix *  new_picture.m_ExtMatrix;
		new_picture.m_SpinningMatrix = new_picture.m_ExtMatrix.block(0, 0, 3, 3);
		new_picture.m_MovingMatrix = new_picture.m_ExtMatrix.block(0, 3, 3, 1);
		new_picture.m_View = new_picture.m_SpinningMatrix.inverse() * new_picture.m_MovingMatrix;
		Eigen::Vector3f x_direction(1, 0, 0);
		Eigen::Vector3f y_direction(0, 1, 0);
		Eigen::Vector3f z_direction(0, 0, 1);
		new_picture.m_SeeingDirection = new_picture.m_SpinningMatrix.inverse() * z_direction;
		new_picture.m_UpDirection = - new_picture.m_SpinningMatrix.inverse() * y_direction;
		picture_list->push_back(new_picture);
	}
}

void g_PictureInput(const char* pDir, const char* pPrefix, const char* pSuffix, std::vector<TPicture>* picture_list)
{
	int fileCount = picture_list->size();
	std::string fileName(pDir);
	fileName += '/';
	fileName += pPrefix;
	for (int i = 0; i < fileCount; i++)
	{
		std::cout << fileName + std::to_string(i) + pSuffix << std::endl;
		(*picture_list)[i].m_Image = cv::imread(fileName + std::to_string(i) + pSuffix);
	}
}

void g_RenewRGB()
{
	for (int i = 0; i < g_PointList.size(); i++)
	{
		Eigen::Vector3i color_all(0,0,0);
		for (int j = 0; j < g_PointList[i].m_SonMesh.size(); j++)
		{
			int find_group = g_MeshList[g_PointList[i].m_SonMesh[j]].m_Belong;
			int k = 0;
			for (k = 0; k < g_PointList[i].m_GroupList.size(); k++)
			{
				if (g_PointList[i].m_GroupList[k] == find_group) break;
			}
			color_all += g_PointList[i].m_ColorList[k];
		}
		g_PointList[i].m_RGB = color_all / g_PointList[i].m_SonMesh.size();
	}
}

void g_PlyOutput(char filename[])
{
	ofstream file(filename);
	file << "ply" << endl;
	file << "format ascii 1.0" << endl;
	file << "element vertex "<<g_PointList.size() << endl;
	file << "property float32 x" << endl;
	file << "property float32 y" << endl;
	file << "property float32 z" << endl;
	file << "property uint8 red" << endl;
	file << "property uint8 green" << endl;
	file << "property uint8 blue" << endl;
	file << "element face " << g_MeshList.size() << endl;
	file << "property list uint8 int32 vertex_index" << endl;
	file << "end_header"<<endl;
	for (int i = 0; i < g_PointList.size(); i++)
	{
		file << g_PointList[i].m_Place(0) << " " << g_PointList[i].m_Place(1) << " " << g_PointList[i].m_Place(2) << " ";
		file << g_PointList[i].m_RGB(0) << " " << g_PointList[i].m_RGB(1) << " " << g_PointList[i].m_RGB(2) << " " << endl;
	}
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		file << "3 " << g_MeshList[i].m_Vertex(0) << " " << g_MeshList[i].m_Vertex(1) << " " << g_MeshList[i].m_Vertex(2) << " " << endl;
	}
}
