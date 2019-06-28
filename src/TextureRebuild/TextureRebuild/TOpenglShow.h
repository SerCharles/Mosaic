/*
filename: TOpenglShow.h
description: use opengl to show the texture of the model
date: 2/28/2019
*/

#ifndef TOPENGLSHOW_H
#define TOPENGLSHOW_H





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


extern unsigned int g_TextureID[20];

extern int g_Phi, g_Theta;
extern double g_Radius;
extern bool g_MouseDown;
extern int g_MouseX, g_MouseY;

extern void g_GetProjection(int point_num, int picture_num, double& x, double& y);

extern void g_LoadTexture();

extern void g_DrawTexture();

extern void g_ControlClick(int button, int state, int x, int y);

extern void g_HandleMouse(int x, int y);

#endif