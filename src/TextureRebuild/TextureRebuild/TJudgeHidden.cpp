/*
filename: TJudgeHidden.h
description: use opengl to judge whether a mesh is hidden
date: 2/23/2019
*/




//system requirements
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <windows.h>
using namespace std;

//cv requirements
#include <opencv2/opencv.hpp>

//algebra requirements
#include <Eigen/Eigen>
#include <limits>

//opengl
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <gl/glut.h>
#include <GL/GLAux.h>

#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glu32.lib")
#pragma comment(lib,"glut32.lib")


//personal requirements
#include "TDefine.hpp"
#include "TJudgeHidden.h"

#define D_Y 0.1
#define SIZE 110000

#define WIDTH 1024
#define HEIGHT 768
#define MAX_DIST 5
#define THRESHOLD 0.8
/*
void GetColor(int num, int total, int& r, int& g, int& b)
{
	r = ((num + 1) / (25 * 25)) * 10;
	g = (((num + 1) - (r / 10) * (25 * 25)) / 25) * 10;
	b = ((num + 1) % 25) * 10;
	
}

void GetNum(int& num, int total, int r, int g, int b)
{
	b /= 10;
	g /= 10;
	r /= 10;
	num = r * 25 * 25 + g * 25 + b;
	num -= 1;
}*/

void GetColor(int num, int total, int& r, int& g, int& b)
{
	r = ((num + 1) / (64 * 64)) * 4;
	g = (((num + 1) - (r / 4) * (64 * 64)) / 64) * 4;
	b = ((num + 1) % 64) * 4;

}

void GetNum(int& num, int total, int r, int g, int b)
{
	b /= 4;
	g /= 4;
	r /= 4;
	num = r * 64 * 64 + g * 64 + b;
	num -= 1;
}

bool g_WhetherCanSee[20][SIZE + 10] = { 0 };
unsigned int g_TextureID[20] = { 0 };


/*
void g_DrawPicture(const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList)
{
	//fstream file("cansee.txt", ios::out|ios::app);
	//draw

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);



	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		//glLoadIdentity();
															//glScalef(0.01, 0.01, 0.01);

															//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (int i = 0; i < g_MeshList.size(); i++)
	{

		int c1 = 0, c2 = 0;
		if (i % 4 == 0) c1 = 0; c2 = 0;
		if (i % 4 == 1) c1 = 0; c2 = 1;
		if (i % 4 == 2) c1 = 1; c2 = 0;
		if (i % 4 == 3) c1 = 1; c2 = 1;
		glPolygonOffset(-1.0f * c1, -1.0f * c2);

		int sequence_0 = 255;
		int sequence_1 = 255;
		int sequence_2 = 255;

		GetColor(i, g_MeshList.size(), sequence_0, sequence_1, sequence_2);
		GLubyte color_0 = GLubyte(sequence_0);
		GLubyte color_1 = GLubyte(sequence_1);
		GLubyte color_2 = GLubyte(sequence_2);

		glBegin(GL_TRIANGLES);

		Eigen::Vector3f point_0 = g_PointList[g_MeshList[i].m_Vertex(0)].m_Place;
		glColor3ub(color_0, color_1, color_2);
		glVertex3f(point_0(0), point_0(1), point_0(2));

		Eigen::Vector3f point_1 = g_PointList[g_MeshList[i].m_Vertex(1)].m_Place;
		//glColor3ub(color_0, color_1, color_2);
		glVertex3f(point_1(0), point_1(1), point_1(2));

		Eigen::Vector3f point_2 = g_PointList[g_MeshList[i].m_Vertex(2)].m_Place;
		//glColor3ub(color_0, color_1, color_2);
		glVertex3f(point_2(0), point_2(1), point_2(2));
		glEnd();
	}


}

unsigned char color_red[650000] = { 0 };
unsigned char color_green[650000] = { 0 };
unsigned char color_blue[650000] = { 0 };

void g_JudgeHidden(const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList, int j)
{
	memset(color_red, 0, sizeof(color_red));
	memset(color_green, 0, sizeof(color_green));
	memset(color_blue, 0, sizeof(color_blue));

	GLint eReadFormat;
	glGetIntegerv(GL_IMPLEMENTATION_COLOR_READ_FORMAT, &eReadFormat);
	glReadPixels(0, 0, 800, 800, GL_RED, GL_UNSIGNED_BYTE, color_red);
	glReadPixels(0, 0, 800, 800, GL_GREEN, GL_UNSIGNED_BYTE, color_green);
	glReadPixels(0, 0, 800, 800, GL_BLUE, GL_UNSIGNED_BYTE, color_blue);

	for (int i = 0; i <= 640000; i++)
	{
		int answer = 0;
		GetNum(answer, g_MeshList.size(), int(color_red[i]), int(color_green[i]), int(color_blue[i]));

		if (answer >= 0 && answer < g_MeshList.size())
		{
			g_WhetherCanSee[j][answer] = 1;
		}
	}
}


*/



void g_DrawPicture(const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList)
{
	//fstream file("cansee.txt", ios::out|ios::app);
		//draw

		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LEQUAL);
		glEnable(GL_POLYGON_OFFSET_FILL);
		//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);



		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		//glLoadIdentity();
		//glScalef(0.01, 0.01, 0.01);

		//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		for (int i = 0; i < g_MeshList.size(); i++)
		{

			
			int sequence_0 = 255;
			int sequence_1 = 255;
			int sequence_2 = 255;

			//GetColor(i, g_MeshList.size(), sequence_0, sequence_1, sequence_2);
			GLubyte color_0 = GLubyte(sequence_0);
			GLubyte color_1 = GLubyte(sequence_1);
			GLubyte color_2 = GLubyte(sequence_2);
			
			glBegin(GL_TRIANGLES);
			
			Eigen::Vector3f point_0 = g_PointList[g_MeshList[i].m_Vertex(0)].m_Place;
			glColor3ub(color_0, color_1, color_2);
			glVertex3f(point_0(0), point_0(1), point_0(2));

			Eigen::Vector3f point_1 = g_PointList[g_MeshList[i].m_Vertex(1)].m_Place;
			//glColor3ub(color_0, color_1, color_2);
			glVertex3f(point_1(0), point_1(1), point_1(2));

			Eigen::Vector3f point_2 = g_PointList[g_MeshList[i].m_Vertex(2)].m_Place;
			//glColor3ub(color_0, color_1, color_2);
			glVertex3f(point_2(0), point_2(1), point_2(2));
			glEnd();
		}

		//glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);

		for (int i = 0; i < g_PointList.size(); i++)
		{
			/*int c1 = 0, c2 = 0;
			if (i % 4 == 0) c1 = 0; c2 = 0;
			if (i % 4 == 1) c1 = 0; c2 = 1;
			if (i % 4 == 2) c1 = 1; c2 = 0;
			if (i % 4 == 3) c1 = 1; c2 = 1;
			glPolygonOffset(-1.0f * c1, -1.0f * c2);*/

			int sequence_0 = 0;
			int sequence_1 = 0;
			int sequence_2 = 0;
			GetColor(i, g_PointList.size(), sequence_0, sequence_1, sequence_2);

			GLubyte color_0 = GLubyte(sequence_0);
			GLubyte color_1 = GLubyte(sequence_1);
			GLubyte color_2 = GLubyte(sequence_2);

			Eigen::Vector3f point_0 = g_PointList[i].m_Place;
			glColor3ub(color_0, color_1, color_2);
			glTranslatef(point_0(0), point_0(1), point_0(2));

			int c1 = 0, c2 = 0;
			if (i % 4 == 0) c1 = 0; c2 = 0;
			if (i % 4 == 1) c1 = 0; c2 = 1;
			if (i % 4 == 2) c1 = 1; c2 = 0;
			if (i % 4 == 3) c1 = 1; c2 = 1;
			glPolygonOffset(-1.0f * c1, -1.0f * c2);

			glutSolidCube(D_Y);
			glTranslatef(-point_0(0), -point_0(1), -point_0(2));

		}
}

unsigned char color_red[650000] = { 0 };
unsigned char color_green[650000] = { 0 };
unsigned char color_blue[650000] = { 0 };

void g_JudgeHidden(const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList, int j)
{
	memset(color_red, 0, sizeof(color_red));
	memset(color_green, 0, sizeof(color_green));
	memset(color_blue, 0, sizeof(color_blue));

	GLint eReadFormat;
	glGetIntegerv(GL_IMPLEMENTATION_COLOR_READ_FORMAT, &eReadFormat);
	glReadPixels(0, 0, 800, 800, GL_RED, GL_UNSIGNED_BYTE, color_red);
	glReadPixels(0, 0, 800, 800, GL_GREEN, GL_UNSIGNED_BYTE, color_green);
	glReadPixels(0, 0, 800, 800, GL_BLUE, GL_UNSIGNED_BYTE, color_blue);

	for (int i = 0; i <= 640000; i++)
	{
		int answer = 0;
		GetNum(answer, g_MeshList.size(), int(color_red[i]), int(color_green[i]), int(color_blue[i]));

		if (answer >= 0 && answer < g_PointList.size())
		{
			g_WhetherCanSee[j][answer] = 1;
		}

	}


}

void g_RenewHidden(const std::vector<TPoint>& g_PointList,
	const std::vector<TMesh>& g_MeshList, const std::vector<TPicture>& g_PictureList, int j)
{
	for (int i = 0; i < g_PointList.size(); i++)
	{
		vector<int> Adjacent_List;
		int head = 0, tail = 0, step = 0;
		Adjacent_List.push_back(i);
		while (step < MAX_DIST)
		{
			step++;
			for (int ii = head; ii <= tail; ii++)
			{
				int current = Adjacent_List[ii];
				for (int jj = 0; jj < g_PointList[current].m_PointLink.size(); jj++)
				{
					int nova = g_PointList[current].m_PointLink[jj];
					bool whether_repeat = 0;
					
					for (int kk = 0; kk < Adjacent_List.size(); kk++)
					{
						if (Adjacent_List[kk] == nova)
						{
							whether_repeat = 1;
							break;
						}
					}

					if (whether_repeat == 0)
					{
						Adjacent_List.push_back(nova);
					}
				}
			}
			head = tail + 1;
			tail = Adjacent_List.size() - 1;
		}

		int see_num = 0;
		for (int ii = 0; ii < Adjacent_List.size(); ii++)
		{
			if (g_WhetherCanSee[j][Adjacent_List[ii]] == 1) see_num++;
		}
		if (double(see_num) / double(Adjacent_List.size()) >= THRESHOLD)
		{
			g_WhetherCanSee[j][i] = 1;
		}
		else if (double(see_num) / double(Adjacent_List.size()) <= 1.00 - THRESHOLD)
		{
			g_WhetherCanSee[j][i] = 0;
		}
	}
}
