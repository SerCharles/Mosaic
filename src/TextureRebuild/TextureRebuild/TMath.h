/*
filename: TMath.h
description: used in Mathematical calculation, such as calculating normal, calculating weight
date: 2/16/2019
*/

#ifndef TMATH_H
#define TMATH_H





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

//personal requirements
#include "TDefine.hpp"

extern Eigen::Vector3f g_CalculateNormal(Eigen::Vector3f x, Eigen::Vector3f y, Eigen::Vector3f z);

extern double g_CalculateCosine(const Eigen::Vector3f x, const Eigen::Vector3f y);

extern double g_GetSelfWeight(const TMesh& the_mesh, const TPicture& the_picture, const std::vector<TPoint>& g_PointList)
;

extern double g_GetMutualWeight(const TMesh& mesh_first, const TPicture& picture_first,
	const TMesh& mesh_second, const TPicture& picture_second, const std::vector<TPoint>& g_PointList);

extern void g_GetCommonPoints(int& start_node, int& end_node, const TMesh& mesh_first, const TMesh& mesh_second);

extern double g_GetTotalWeight(int mosaic_list[], const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList);

extern void g_GetColor(std::vector<TPoint>& g_PointList, const std::vector<TMesh>& g_MeshList,
	const std::vector<TPicture>& g_PictureList, const std::vector<TGroup>& g_GroupList);


#endif