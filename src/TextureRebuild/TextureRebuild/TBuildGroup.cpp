/*
filename: TBuildGroup.cpp
description: using tarjan to build group
date: 2/21/2019
*/





//system requirements
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;


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
#include "TBuildGroup.h"

#define SIZE 110000
int find_time[SIZE + 10] = { 0 }, low_time[SIZE + 10] = { 0 };
bool visit[SIZE + 10] = { 0 };
vector<int> stack;
int current_time = 0;


void Tarjan(int current_point)//����ڼ������ڴ����ݹ���ǵ㡣
{

	current_time++;
	find_time[current_point] = current_time;// �½���ĳ�ʼ����
	low_time[current_point] = current_time;// �½���ĳ�ʼ����

	stack.push_back(current_point);
	visit[current_point] = 1;//��ʾ��ջ��
	for (int i = 0; i < g_MeshList[current_point].m_MeshLink.size(); i++)
	{
		int next = g_MeshList[current_point].m_MeshLink[i];
		if (g_MeshList[current_point].m_Mosaic == g_MeshList[next].m_Mosaic)
		{
			if (find_time[next] == 0)
			{//���û���ʹ�
				Tarjan(next);//���½������죬��ʼ�ݹ�
				low_time[current_point] = min(low_time[current_point], low_time[next]);//�ݹ�������Ƚ�˭��˭�Ķ��ӣ����ף��������Ķ�Ӧ��ϵ���漰��ǿ��ͨ����������С�������顣
			}
			else if (visit[next])
			{  //������ʹ������һ���ջ�
				low_time[current_point] = min(low_time[current_point], find_time[next]);//�Ƚ�˭��˭�Ķ��ӣ����ס��������Ӷ�Ӧ��ϵ
			}
		}
	}
	
	if (low_time[current_point] == find_time[current_point]) //����������ǿ��ͨ�������������С����
	{
		TGroup new_group(g_GroupList.size(), g_MeshList[current_point].m_Mosaic);
		int top = -1;
		while (current_point != top && !stack.empty())
		{
			top = stack[stack.size() - 1];
			stack.pop_back();
			visit[top] = 0;
			g_MeshList[top].m_Belong = new_group.m_Number;
			new_group.m_AddMesh(top);
		}
		g_GroupList.push_back(new_group);
	}
}

void g_BuildGroup()
{
	stack.clear();
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		if (!find_time[i])
		{
			if (i >= 30000)
			{
				int frog = 1;
			}
			Tarjan(i);//�������û�з��ʹ����ʹӴ˵㿪ʼ����ֹͼû����
		}
	}
}



void g_BuildGroup(int pp)
{
	int last_unvisit = 0;
	while (last_unvisit < g_MeshList.size())
	{
		TGroup new_group(g_GroupList.size(), g_MeshList[last_unvisit].m_Mosaic);
		queue<int> que;
		que.push(last_unvisit);
		visit[last_unvisit] = 1;
		g_MeshList[last_unvisit].m_Belong = g_GroupList.size();
		new_group.m_MeshList.push_back(last_unvisit);
		while (!que.empty())
		{
			int top = que.front();
			que.pop();
			for (int i = 0; i < g_MeshList[top].m_MeshLink.size(); i++)
			{
				int next = g_MeshList[top].m_MeshLink[i];
				if (g_MeshList[top].m_Mosaic == g_MeshList[next].m_Mosaic && visit[next] == 0)
				{
					que.push(next);
					visit[next] = 1;
					g_MeshList[next].m_Belong = g_GroupList.size();
					new_group.m_MeshList.push_back(next);
				}
			}
		}
		g_GroupList.push_back(new_group);
		while (visit[last_unvisit] != 0)
		{
			last_unvisit++;
		}
	}
}

bool JudgeRepeatPP(int point, int group)
{
	for (int i = 0; i < g_PointList[point].m_GroupList.size(); i++)
	{
		if (group == g_PointList[point].m_GroupList[i])
		{
			return 1;
		}
	}
	return 0;
}

bool JudgeRepeatReverse(int point, int group)
{
	for (int i = 0; i < g_GroupList[group].m_PointList.size(); i++)
	{
		if (point == g_GroupList[group].m_PointList[i])
		{
			return 1;
		}
	}
	return 0;
}

void g_LinkGroup()
{
	for (int i = 0; i < g_PointList.size(); i++)
	{
		for (int j = 0; j < g_PointList[i].m_SonMesh.size(); j++)
		{
			int new_group = g_MeshList[g_PointList[i].m_SonMesh[j]].m_Belong;
			if (JudgeRepeatPP(i, new_group) == 0)
			{
				g_PointList[i].m_GroupList.push_back(new_group);
			}
			if (JudgeRepeatReverse(i, new_group) == 0)
			{
				g_GroupList[new_group].m_PointList.push_back(i);
			}
		}
	}
}