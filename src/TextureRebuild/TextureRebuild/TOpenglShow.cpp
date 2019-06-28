/*
filename: TOpenglShow.cpp
description: use opengl to show the texture of the model
date: 2/28/2019
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
#include "TOpenglShow.h"
#include "TInput.h"

#define D_Y 0.1
#define SIZE 110000

#define WIDTH 1024
#define HEIGHT 768

#define PI 3.1415926
int g_Phi = 0, g_Theta = 0;
double g_Radius = 0.5;
bool g_MouseDown = 0;
int g_MouseX = 0, g_MouseY = 0;



void g_LoadTexture()
{
	for (int i = 0; i < g_PictureList.size(); i++)
	{
		GLint width = g_PictureList[i].m_Image.cols;
		GLint height = g_PictureList[i].m_Image.rows;

		int bytes = width * height * 3;

		//allocate space for GL picture
		GLubyte* picture_opengl = new GLubyte[bytes];

		memcpy(picture_opengl, g_PictureList[i].m_Image.data,
			sizeof(unsigned char) * bytes);




		//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		//glDisable(GL_LIGHTING);
		glEnable(GL_TEXTURE_2D);

		glGenTextures(1, &g_TextureID[i]);
		glBindTexture(GL_TEXTURE_2D, g_TextureID[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0,
			GL_BGR_EXT, GL_UNSIGNED_BYTE, picture_opengl);
		delete(picture_opengl);
	}
}

void g_GetProjection(int point_num, int picture_num, double& x, double& y)
{
	Eigen::Vector4f multiply_place(g_PointList[point_num].m_Place(0),
		g_PointList[point_num].m_Place(1), g_PointList[point_num].m_Place(2), 1);
	Eigen::Vector3f project = g_PictureList[picture_num].m_IntMatrix *
		g_PictureList[picture_num].m_ExtMatrix * multiply_place;
	int place_x = project(0) / project(2);
	int place_y = project(1) / project(2);
	x = place_x;
	y = place_y;
}


void g_DrawTexture()
{

	glScalef(0.10, 0.10, 0.03);
	glTranslatef(-g_Center(0), -g_Center(1), -g_Center(2));
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	for (int i = 0; i < g_MeshList.size(); i++)
	{
		int mosaic = 0;
		mosaic = g_MeshList[i].m_Mosaic;
		double place_x_0 = 500, place_y_0 = 200;
		double place_x_1 = 500, place_y_1 = 200;
		double place_x_2 = 500, place_y_2 = 200;

		g_GetProjection(g_MeshList[i].m_Vertex(0), mosaic, place_x_0, place_y_0);
		g_GetProjection(g_MeshList[i].m_Vertex(1), mosaic, place_x_1, place_y_1);
		g_GetProjection(g_MeshList[i].m_Vertex(2), mosaic, place_x_2, place_y_2);

		place_y_0 /= HEIGHT; place_y_1 /= HEIGHT; place_y_2 /= HEIGHT;
		place_x_0 /= WIDTH; place_x_1 /= WIDTH; place_x_2 /= WIDTH;

		glColor3f(1, 1, 1);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, g_TextureID[mosaic]);

		glBegin(GL_TRIANGLES);

		//all coordinate(vertex or texture) must be between -1 and 1
		Eigen::Vector3f point_0 = g_PointList[g_MeshList[i].m_Vertex(0)].m_Place;
		glTexCoord2f(place_x_0, place_y_0);
		glVertex3f(point_0(0), point_0(1), point_0(2));

		Eigen::Vector3f point_1 = g_PointList[g_MeshList[i].m_Vertex(1)].m_Place;
		glTexCoord2f(place_x_1, place_y_1);
		glVertex3f(point_1(0), point_1(1), point_1(2));

		Eigen::Vector3f point_2 = g_PointList[g_MeshList[i].m_Vertex(2)].m_Place;
		glTexCoord2f(place_x_2, place_y_2);
		glVertex3f(point_2(0), point_2(1), point_2(2));
		glEnd();
		glDisable(GL_TEXTURE_2D);
		//glutSwapBuffers();
	}
	glFlush();
}

void g_SetCamera()
{
	double place_x = 0, place_y = 0, place_z = 0;
	double up_x = 0, up_y = 0, up_z = 0;
	double theta = g_Theta * PI / 180;
	double phi = g_Phi * PI / 180;
	place_x = g_Radius * cos(theta) * cos(phi);
	place_z = g_Radius * cos(theta) * sin(phi);
	place_y = g_Radius * sin(theta);
	if (abs(sin(theta)) <= 0.05)
	{
		up_x = 0;
		up_y = - 0.1;
		up_z = 0;
	}
	else
	{
		double nova = g_Radius / sin(theta);
		if (nova <= 0)
		{
			up_x = -place_x;
			up_z = -place_z;
			up_y = nova - place_y;
			up_x /= 15; up_y /= 15; up_z /= 15;
		}
		else
		{
			up_x = place_x;
			up_z = place_z;
			up_y = 0 - nova + place_y;
			up_x /= 15; up_y /= 15; up_z /= 15;
		}
	}
	glLoadIdentity();
	gluLookAt(place_x, place_y, place_z, 0, 0, 0, up_x, up_y, up_z);

	g_DrawTexture();

}


void g_HandleMouse(int x, int y)
{
	if (g_MouseDown == 1)
	{
		int dx = x - g_MouseX;
		int dy = y - g_MouseY;
		g_Phi += (dx / 2);
		g_Theta += (dy / 2);
		g_Phi = (g_Phi + 360) % 360;

		if (g_Theta >= 90) g_Theta = 89;
		else if (g_Theta <= -90) g_Theta = -89;

		g_MouseX += dx;
		g_MouseY += dy;
		g_SetCamera();
	}
}

void g_ControlClick(int button, int state, int x, int y)
{
	g_MouseX = x;
	g_MouseY = y;

	if (button == GLUT_LEFT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			g_MouseDown = true;
		}
		else
		{
			g_MouseDown = false;
		}
	}
	else
	{
		g_MouseDown = 0;
	}
}

