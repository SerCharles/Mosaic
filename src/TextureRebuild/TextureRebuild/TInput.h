/*
filename: TInput.h
description: used in inputing .ply files (for those with only triangle meshes), matrix files and pictures
also includes .ply output
date: 2/15/2019
*/

#ifndef TINPUT_H
#define TINPUT_H

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

extern std::vector<TPoint> g_PointList;
extern std::vector<TMesh> g_MeshList;
extern std::vector<TPicture> g_PictureList;
extern std::vector<TGroup> g_GroupList;
extern Eigen::Vector3f g_Center;

//input ply
extern void g_PlyInput(char filename[], std::vector<TPoint>* point_list, std::vector<TMesh>* mesh_list);

//input camera matrix
extern void g_MatrixInput(char filename[], std::vector<TPicture>* point_list);

//input pictures
extern void g_PictureInput(const char* pDir, const char* pPrefix, const char* pSuffix, std::vector<TPicture>* picture_list);

//output ply
extern void g_PlyOutput(char filename[]);

extern void g_RenewRGB();

#endif