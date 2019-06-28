/*
filename: TMath.cpp
description: used in Mathematical calculation, such as calculating normal, calculating weight
date: 2/16/2019
*/


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
#include <Eigen/SVD>  



//personal requirements
#include "TDefine.hpp"
#include "TMath.h"
#include "TJudgeHidden.h"

#define D_X 0.05
#define lambda 0.004;
#define MAX_NUM 114514


Eigen::Vector3f g_CalculateNormal(Eigen::Vector3f x, Eigen::Vector3f y, Eigen::Vector3f z)
{
	Eigen::Vector3f a = y - x;
	Eigen::Vector3f b = z - x;
	Eigen::Vector3f answer;
	answer = a.cross(b);
	return answer;
}

double g_CalculateCosine(const Eigen::Vector3f x,const Eigen::Vector3f y)
{
	double a = x.transpose() * y;
	double b = x.norm() * y.norm();
	return a / b;
}


//calculate self value 
double g_GetSelfWeight(const TMesh& the_mesh, const TPicture& the_picture, const std::vector<TPoint>& g_PointList)
{
	if (g_WhetherCanSee[the_picture.m_Number][the_mesh.m_Vertex(0)] == 0 ||
		g_WhetherCanSee[the_picture.m_Number][the_mesh.m_Vertex(1)] == 0 ||
		g_WhetherCanSee[the_picture.m_Number][the_mesh.m_Vertex(2)] == 0 )
	{
		return MAX_NUM;
	}
	//calculate
	Eigen::Vector3f middle = (g_PointList[the_mesh.m_Vertex(0)].m_Place +
		g_PointList[the_mesh.m_Vertex(1)].m_Place + g_PointList[the_mesh.m_Vertex(2)].m_Place) / 3;
	double cosine = g_CalculateCosine(the_mesh.m_Normal, the_picture.m_SeeingDirection);
	double answer = 1 - cosine * cosine;
	return answer;
}

double g_GetDistanceRGB(int r1, int r2, int g1, int g2, int b1,int b2)
{
	double rmean = (r1 + r2) / 2;
	double ans = 0;
	double r = r1 - r2, g = g1 - g2, b = b1 - b2;
	//ans = sqrt((2 + rmean / 256)*(r*r) + 4 * g*g + (2 + 255 / 256 - r / 256)*b*b);
	ans = sqrt(r*r + g * g + b * b);
	return ans;
}


//calculate mutualvalue
double g_GetMutualWeight(const TMesh& mesh_first, const TPicture& picture_first,
	const TMesh& mesh_second, const TPicture& picture_second, const std::vector<TPoint>& g_PointList)
{
	double answer = 0;
	int start_node = -1, end_node = -1;
	g_GetCommonPoints(start_node, end_node, mesh_first, mesh_second);
	if (start_node == -1 || end_node == -1) return -1;
	Eigen::Vector3f start_place = g_PointList[start_node].m_Place;
	Eigen::Vector3f end_place = g_PointList[end_node].m_Place;
	double length = (end_place - start_place).norm();
	double covered = 0;
	while (covered <= length)
	{
		Eigen::Vector3f current_place = start_place + (end_place - start_place) * covered;

		//calculate RGB eucilidian distance
		Eigen::Vector3f color_first, color_second;
		int place_first_x = 0, place_first_y = 0, place_second_x = 0, place_second_y = 0;
		Eigen::Vector4f multiply_place(current_place(0), current_place(1), current_place(2), 1);

		//get project place
		Eigen::Vector3f project_first = picture_first.m_IntMatrix * picture_first.m_ExtMatrix
			* multiply_place;
		double z1 = project_first(2);
		project_first /= z1;

		Eigen::Vector3f project_second = picture_second.m_IntMatrix * picture_second.m_ExtMatrix
			* multiply_place;
		double z2 = project_second(2);
		project_second /= z2;

		place_first_x = int(project_first(0)); place_first_y = int(project_first(1));
		place_second_x = int(project_second(0)); place_second_y = int(project_second(1));

		//judge whether valid
		if (place_first_x < 0 || place_first_x > picture_first.m_Image.cols ||
			place_first_y < 0 || place_first_y > picture_first.m_Image.rows ||
			place_second_x < 0 || place_second_x > picture_second.m_Image.cols ||
			place_second_y < 0 || place_second_y > picture_first.m_Image.rows)
		{
			return -1;
		}

		cv::Vec3b color_first_cv = picture_first.m_Image.at<cv::Vec3b>(place_first_y, place_first_x);
		cv::Vec3b color_second_cv = picture_second.m_Image.at<cv::Vec3b>(place_second_y, place_second_x);

		color_first << color_first_cv[2], color_first_cv[1], color_first_cv[0];
		color_second << color_second_cv[2], color_second_cv[1], color_second_cv[0];

		double distance = 0;
		/*distance = g_GetDistanceRGB(color_first(0), color_second(0), color_first(1),
			color_second(1), color_first(2), color_second(2)) * D_X;*/
		distance = (color_first - color_second).norm()  * D_X * lambda;
		
		answer += distance;
		covered += D_X;
	}
	return answer;
}


void g_GetCommonPoints(int& start_node, int& end_node, const TMesh& mesh_first, const TMesh& mesh_second)
{
	for (int i = 0; i <= 2; i++)
	{
		for (int j = 0; j <= 2; j++)
		{
			if (mesh_first.m_Vertex[i] == mesh_second.m_Vertex[j])
			{
				if (start_node == -1) start_node = i;
				else end_node = i;

				break;
			}
		}
	}
}

double g_GetTotalWeight(int mosaic_list[], const std::vector<TPoint>& g_PointList, 
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList)
{
	double total_energy = 0;
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		double self = g_GetSelfWeight(g_MeshList[i], g_PictureList[mosaic_list[i]], g_PointList);
		total_energy += self;
		for (int j = 0; j < g_MeshList[i].m_MeshLink.size(); j++)
		{
			int next = g_MeshList[i].m_MeshLink[j];
			
			if (i < next && mosaic_list[i] != mosaic_list[next])
			{
				double mutual = g_GetMutualWeight(g_MeshList[i], g_PictureList[mosaic_list[i]],
					g_MeshList[next], g_PictureList[mosaic_list[next]], g_PointList);
				total_energy += mutual;
			}
			
		}
	}
	return total_energy;
}

void g_GetColor(std::vector<TPoint>& g_PointList,const std::vector<TMesh>& g_MeshList,
	const std::vector<TPicture>& g_PictureList, const std::vector<TGroup>& g_GroupList)
{
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_GroupList.size(); j++)
		{
			int picture = g_GroupList[g_PointList[i].m_GroupList[j]].m_Mosaic;
			Eigen::Vector4f multiply_place(g_PointList[i].m_Place(0), g_PointList[i].m_Place(1), g_PointList[i].m_Place(2), 1);
			Eigen::Vector3f project = g_PictureList[picture].m_IntMatrix *
				g_PictureList[picture].m_ExtMatrix * multiply_place;
			int place_x = project(0) / project(2);
			int place_y = project(1) / project(2);
			cv::Vec3b color_cv = g_PictureList[picture].m_Image.at<cv::Vec3b>(place_y, place_x);
			Eigen::Vector3i color;
			color << color_cv[2], color_cv[1], color_cv[0];
			g_PointList[i].m_ColorList.push_back(color);
		}
	}
	/*for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_SonMesh.size(); j++)
		{
			int picture = g_MeshList[j].m_Mosaic;
			Eigen::Vector4f multiply_place(g_PointList[i].m_Place(0), g_PointList[i].m_Place(1), g_PointList[i].m_Place(2), 1);
			Eigen::Vector3f project = g_PictureList[picture].m_IntMatrix *
				g_PictureList[picture].m_ExtMatrix * multiply_place;
			int place_x = project(0) / project(2);
			int place_y = project(1) / project(2);
			cv::Vec3b color_cv = g_PictureList[picture].m_Image.at<cv::Vec3b>(place_y, place_x);
			Eigen::Vector3i color;
			color << color_cv[2], color_cv[1], color_cv[0];
			g_PointList[i].m_RGB += color;
		}
		g_PointList[i].m_RGB /= g_PointList[i].m_SonMesh.size();
	}*/
}

