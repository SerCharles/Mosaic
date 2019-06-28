/*
filename:main.cpp
description:the main entrance of the file
date: 2/15/2019
*/

#pragma warning(disable:4819)
#pragma warning(disable:4244)
#pragma warning(disable:4267)

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
#include <GL/GL.h>
#include <GL/GLU.h>
#include <gl/glut.h>
#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glu32.lib")
#pragma comment(lib,"glut32.lib")

//cuda
#include <cusolverSp.h>
#include <cuda_runtime_api.h>


//personal requirements
#include "TDefine.hpp"
#include "TInput.h"
#include "TMath.h"
#include "TExpansion.h"
#include "TBuildGroup.h"
#include "TSolver.cuh"
#include "TJudgeHidden.h"
#include "TOpenglShow.h"

#define MAX_ITERATION_TIME 3



TMeshGraph* the_graph = NULL;
TSolver* the_solver = NULL;
char filename_mesh[30] = "../../input/mesh30000.ply";
char filename_matrix[30] = "../../calibParamsI.txt";
char filename_out[30] = "mesh_with_texture.ply";

//�ж��Ƿ�ɼ���OPENGL�����������expansion��solve��������
void mydisplay()
{
	//g_LoadTexture();
	
	for (int i = 0; i < g_PictureList.size(); i++)
	{
		//set
		glLoadIdentity();

		glScalef(0.10, 0.10, 0.03);
		Eigen::Vector3f camera_place = -g_PictureList[i].m_View;
		Eigen::Vector3f seeing_direction = - g_PictureList[i].m_SeeingDirection;
		Eigen::Vector3f seeing_place = camera_place + seeing_direction;
		Eigen::Vector3f up_direction = g_PictureList[i].m_UpDirection;
		gluLookAt(camera_place(0), camera_place(1), camera_place(2),
			seeing_place(0), seeing_place(1), seeing_place(2),
			up_direction(0), up_direction(1), up_direction(2));

		g_DrawPicture(g_PointList, g_MeshList, g_PictureList);

		glFlush();
		//Sleep(500);
		g_JudgeHidden(g_PointList, g_MeshList, g_PictureList, i);
		//Sleep(500);
	}

	for (int i = 1; i <= MAX_ITERATION_TIME; i++)
	{
		for (int j = 0; j < g_PictureList.size(); j++)
		{
			g_RenewHidden(g_PointList, g_MeshList, g_PictureList, j);
		}
	}


	
	//test whether_can_see
	for (int i = 0; i < g_PictureList.size(); i++)
	{
		string pre = "test/picture";
		string c; 
		c = to_string(i);
		pre += c;
		string d = ".ply";
		pre += d;
		for (int j = 0; j < g_PointList.size(); j++)
		{
			Eigen::Vector3i red(255, 0, 0);
			Eigen::Vector3i blue(0, 0, 255);
		

			if (g_WhetherCanSee[i][j] == 1)
			{
				g_PointList[j].m_RGB = red;
			}
			else
			{
				g_PointList[j].m_RGB = blue;
			}
		}
		char filename[30] = { 0 };
		for (int kk = 0; kk < pre.size(); kk++) filename[kk] = pre[kk];

		g_PlyOutput(filename);
	}
	
	

	
	//test finding_whether_can_see
	ofstream bork("how_much_can_see.txt");
	int num[20] = { 0 };
	int same[10] = { 0 };
	int special_6 = 0, special_7 = 0, special_8 = 0, special_9 = 0;;
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		for (int j = 0; j < g_PictureList.size(); j++)
		{
			if (g_WhetherCanSee[j][i] == 1)
			{
				num[j] ++;
			}
			if (j < 10 && g_WhetherCanSee[j][i] == 1 && g_WhetherCanSee[j + 10][i] == 1)
			{
				same[j]++;
			}
		}
		if (g_WhetherCanSee[2][i] == 1 && g_WhetherCanSee[6][i] == 1) special_6 ++;
		if (g_WhetherCanSee[2][i] == 1 && g_WhetherCanSee[7][i] == 1) special_7 ++;
		if (g_WhetherCanSee[2][i] == 1 && g_WhetherCanSee[8][i] == 1) special_8 ++;
		if (g_WhetherCanSee[2][i] == 1 && g_WhetherCanSee[9][i] == 1) special_9 ++;

	}

	for (int i = 0; i < g_PictureList.size(); i++)
	{
		bork << "picture " << i << ", can see " << num[i] << " points" << endl;
	}
	for (int i = 0; i < 10; i++)
	{
		bork << "picture " << i <<" and "<<i+10<< ", can see the same " << same[i] << " points" << endl;
	}
	bork << "picture " << 2 << " and " << 6 << ", can see the same " << special_6 << " points" << endl;
	bork << "picture " << 2 << " and " << 7 << ", can see the same " << special_7 << " points" << endl;
	bork << "picture " << 2 << " and " << 8 << ", can see the same " << special_8 << " points" << endl;
	bork << "picture " << 2 << " and " << 9 << ", can see the same " << special_9 << " points" << endl;
	


	//expansion
	the_graph = new TMeshGraph();
	the_graph->m_MainExpansion();
	

	
	
	//test
	int num_mesh[20] = { 0 };
	ofstream kebab("how_much_used_for_picture.txt");
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		num_mesh[g_MeshList[i].m_Mosaic] ++;
	}
	for (int i = 0; i < g_PictureList.size(); i++)
	{
		kebab << "picture " << i << " is used for " << num_mesh[i] << " meshes" << endl;
	}
	ofstream out("mosaic_result.txt");
	for (int i = 0; i < g_MeshList.size(); i++) out << g_MeshList[i].m_Mosaic << endl;
	
	
	
	//group
	g_BuildGroup(1);
	g_LinkGroup();
	g_GetColor(g_PointList, g_MeshList, g_PictureList, g_GroupList);
	

	
	//seam leveling
	//the_solver = new TSolver();
	//the_solver->m_GaussNewton();
	//the_solver->m_RenewColor();
	
	
	
	//output
	g_RenewRGB();
	g_PlyOutput(filename_out);
	
	exit(0);
}

int main(int argc, char * argv[])
{
	//the_solver = new TSolver(1);
	//��ȡģ��
	g_PlyInput(filename_mesh, &g_PointList, &g_MeshList);

	// ��ȡ������������
	g_MatrixInput(filename_matrix,&g_PictureList);

	// ��ȡͶӰͼƬ
	g_PictureInput("../../wd_data", "WD2_", "_00020.bmp", &g_PictureList);
	
	/*
	//load mosaic
	cout << "Input 0: without alpha-expansion" << endl;
	cout << "Input 1: with alpha-expansion" << endl;
	bool num = 0;
	cin >> num;

	
	ifstream mosaic("out.txt");
	ifstream mosaic_rough("out_rough.txt");

	if (num == 1)
	{
		for (int i = 0; i < g_MeshList.size(); i++) mosaic >> g_MeshList[i].m_Mosaic;
	}
	else
	{
		for (int i = 0; i < g_MeshList.size(); i++) mosaic_rough >> g_MeshList[i].m_Mosaic;
	}
	*/
	glutInit(&argc, argv);//ʹ��glut����Ҫ���г�ʼ��
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);//�趨������ʾģʽ����ɫģ�ͺͻ��棬������RGB��ɫģ�ͺ͵�����
	glutInitWindowPosition(300, 0);//�趨���ڵĳ�ʼλ�ã���Ļ���Ͻ�Ϊԭ�㣬��λΪ����
	glutInitWindowSize(800, 800);//�趨���ڵĴ�С
	glutCreateWindow("TestWhetherCanSee");//����һ�����ڣ������Ǵ��ڱ�����
	glutDisplayFunc(&mydisplay);//��myDisplayָ��Ϊ��ǰ���ڵ���ʾ���ݺ���
	glutMouseFunc(&g_ControlClick);
	glutMotionFunc(&g_HandleMouse);
	glutMainLoop();//ʹ���ڿ������������ʹ��ʾ�ص�������ʼ����

	return 0;
}


