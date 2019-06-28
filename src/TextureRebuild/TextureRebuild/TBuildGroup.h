/*
filename: TBuildGroup.h
description: using tarjan to build group 
date: 2/21/2019
*/

#ifndef TBUILDGROUP_H
#define TBUILDGROUP_H





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

//tarjian
extern void g_BuildGroup();
extern void g_LinkGroup();

//brute force
extern void g_BuildGroup(int pp);
#endif