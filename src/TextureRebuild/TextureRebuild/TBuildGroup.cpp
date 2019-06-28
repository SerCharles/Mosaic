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


void Tarjan(int current_point)//代表第几个点在处理。递归的是点。
{

	current_time++;
	find_time[current_point] = current_time;// 新进点的初始化。
	low_time[current_point] = current_time;// 新进点的初始化。

	stack.push_back(current_point);
	visit[current_point] = 1;//表示在栈里
	for (int i = 0; i < g_MeshList[current_point].m_MeshLink.size(); i++)
	{
		int next = g_MeshList[current_point].m_MeshLink[i];
		if (g_MeshList[current_point].m_Mosaic == g_MeshList[next].m_Mosaic)
		{
			if (find_time[next] == 0)
			{//如果没访问过
				Tarjan(next);//往下进行延伸，开始递归
				low_time[current_point] = min(low_time[current_point], low_time[next]);//递归出来，比较谁是谁的儿子／父亲，就是树的对应关系，涉及到强连通分量子树最小根的事情。
			}
			else if (visit[next])
			{  //如果访问过，并且还在栈里。
				low_time[current_point] = min(low_time[current_point], find_time[next]);//比较谁是谁的儿子／父亲。就是链接对应关系
			}
		}
	}
	
	if (low_time[current_point] == find_time[current_point]) //发现是整个强连通分量子树里的最小根。
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
			Tarjan(i);//当这个点没有访问过，就从此点开始。防止图没走完
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