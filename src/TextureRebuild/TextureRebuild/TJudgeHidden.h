/*
filename: TJudgeHidden.h
description: use opengl to judge whether a mesh is hidden
date: 2/23/2019
*/

#ifndef TJUDGEHIDDEN_H
#define TJUDGEHIDDEN_H





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

#define SIZE 110000

//ÅÐ¶ÏÊÇ·ñ¿É¼û
extern void g_JudgeHidden(const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList, int j);
	

extern void g_DrawPicture(const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList);

extern bool g_WhetherCanSee[20][SIZE+ 10];

//solve the problem of several points in the front but cannot see
extern void g_RenewHidden(const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList, int j);
#endif