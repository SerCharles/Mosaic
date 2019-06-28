/*
filename: TExpansion.cpp
description: used in doing alpha-expension
date: 2/17/2019
*/
#define MAX_NUM 1145141919810

//system requirements
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <string>
#include <set>
using namespace std;

//cv requirements
#include <opencv2/opencv.hpp>

//algebra requirements
#include <Eigen/Eigen>
#include <limits>

//personal requirements
#include "TDefine.hpp"
#include "TInput.h"
#include "TMath.h"
#include "TExpansion.h"

#define SIZE 110000
#define MAX_NUM 1145141919810
#define MAX 114514

void TMeshGraph::m_MainExpansion()
{
	int mosaic_start[SIZE + 10] = { 0 };
	min_energy = g_GetTotalWeight(mosaic_start, g_PointList, g_MeshList, g_PictureList);
	while (1)
	{
		bool whether_success = 0;
		int mosaic[SIZE+10] = { 0 };
		//main cycle
		for (int i = 0; i < g_PictureList.size(); i++)
		{
			m_MeshGraph.clear();
			m_CurrentPicture = i;
			m_BuildGraph(1);
			m_NetWorkFlow(1);

			for (int j = 0; j < m_SourceNum; j++)
			{
				int current = m_MeshGraph[m_SourceNum].m_Next[j];
				if (m_MeshGraph[current].m_Type == 1 &&
					m_MeshGraph[m_SourceNum].m_Flow[j] == m_MeshGraph[m_SourceNum].m_Capacity[j])
				{
					mosaic[current] = i;
				}
				else
				{
					mosaic[current] = g_MeshList[current].m_Mosaic;
				}
			}
			double new_value = g_GetTotalWeight(mosaic, g_PointList, g_MeshList, g_PictureList);
			if (new_value < min_energy)
			{
				min_energy = new_value;
				for (int j = 0; j < g_MeshList.size(); j++)
				{
					g_MeshList[j].m_Mosaic = mosaic[j];
				}
				whether_success = 1;
			}
		}
		if (whether_success == 0) break;
	}
	m_Purify();
}

void TMeshGraph::m_BuildGraph()
{
	
	//build mesh points
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		TMeshNode new_node(i, 1, i);
		m_MeshGraph.push_back(new_node);
	}

	//build source and end
	m_SourceNum = m_MeshGraph.size();
	m_EndNum = m_SourceNum + 1;
	TMeshNode source(m_SourceNum, 0, -1);
	TMeshNode end(m_EndNum, 0, -1);
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		double self = g_GetSelfWeight(g_MeshList[i], g_PictureList[m_CurrentPicture], g_PointList);
		source.m_Next.push_back(i);
		source.m_Flow.push_back(0);
		source.m_Capacity.push_back(self);
		m_MeshGraph[i].m_Next.push_back(m_SourceNum);
		m_MeshGraph[i].m_Flow.push_back(self);
		m_MeshGraph[i].m_Capacity.push_back(self);


	}
	for (int i = 0; i < m_SourceNum; i++)
	{
		m_MeshGraph[i].m_Next.push_back(m_EndNum);
		m_MeshGraph[i].m_Flow.push_back(0);

		double self = MAX;
		double self_weight = g_GetSelfWeight(g_MeshList[i], g_PictureList[g_MeshList[i].m_Mosaic], g_PointList);
		if (g_MeshList[i].m_Mosaic != m_CurrentPicture)
		{
			self = self_weight;
		}
		m_MeshGraph[i].m_Capacity.push_back(self);
		end.m_Next.push_back(i);
		end.m_Capacity.push_back(self);
		end.m_Flow.push_back(self);
	}
	m_MeshGraph.push_back(source);
	m_MeshGraph.push_back(end);


	//build appendix points
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		int first_node = i;
		for (int j = 0; j < g_MeshList[i].m_MeshLink.size(); j++)
		{
			int second_node = g_MeshList[i].m_MeshLink[j];
			if (first_node >= second_node)
			{
				continue;
			}
			else if (g_MeshList[first_node].m_Mosaic == g_MeshList[second_node].m_Mosaic)
			{
				// no appendix nodes;
				m_MeshGraph[second_node].m_Next.push_back(first_node);

				m_MeshGraph[first_node].m_Next.push_back(second_node);
				m_MeshGraph[first_node].m_Flow.push_back(0);

				double mutual = 0;
				mutual = g_GetMutualWeight(g_MeshList[first_node], g_PictureList[g_MeshList[first_node].m_Mosaic]
					, g_MeshList[second_node], g_PictureList[m_CurrentPicture],g_PointList);
				if (mutual == -1) mutual = MAX;

				m_MeshGraph[first_node].m_Capacity.push_back(mutual);

				m_MeshGraph[second_node].m_Capacity.push_back(mutual);
				m_MeshGraph[second_node].m_Flow.push_back(mutual);

			}
			else
			{
				//have appendix node
				TMeshNode new_node(m_MeshGraph.size(), 0, -1);
				new_node.m_Next.push_back(first_node);

				m_MeshGraph[first_node].m_Next.push_back(m_MeshGraph.size());
				m_MeshGraph[first_node].m_Flow.push_back(0);
				double mutual = 0;

				mutual = g_GetMutualWeight(g_MeshList[first_node], g_PictureList[g_MeshList[first_node].m_Mosaic]
					, g_MeshList[second_node], g_PictureList[m_CurrentPicture], g_PointList);
				if (mutual == -1) mutual = MAX;



				m_MeshGraph[first_node].m_Capacity.push_back(mutual);
				new_node.m_Flow.push_back(mutual);
				new_node.m_Capacity.push_back(mutual);

				m_MeshGraph[second_node].m_Next.push_back(m_MeshGraph.size());
				new_node.m_Next.push_back(second_node);
				new_node.m_Flow.push_back(0);

				mutual = g_GetMutualWeight(g_MeshList[first_node], g_PictureList[m_CurrentPicture]
					, g_MeshList[second_node], g_PictureList[g_MeshList[second_node].m_Mosaic], g_PointList);
				if (mutual == -1) mutual = MAX;

				new_node.m_Capacity.push_back(mutual);

				m_MeshGraph[second_node].m_Capacity.push_back(mutual);
				m_MeshGraph[second_node].m_Flow.push_back(mutual);

				new_node.m_Next.push_back(m_EndNum);
				new_node.m_Flow.push_back(0);

				mutual = g_GetMutualWeight(g_MeshList[first_node], g_PictureList[g_MeshList[first_node].m_Mosaic]
					, g_MeshList[second_node], g_PictureList[g_MeshList[second_node].m_Mosaic], g_PointList);
				if (mutual == -1) mutual = MAX;



				new_node.m_Capacity.push_back(mutual);

				m_MeshGraph[m_EndNum].m_Next.push_back(m_MeshGraph.size());
				m_MeshGraph[m_EndNum].m_Capacity.push_back(mutual);
				m_MeshGraph[m_EndNum].m_Flow.push_back(mutual);


				m_MeshGraph.push_back(new_node);
			}
		}
	}
}

void TMeshGraph::m_Purify()
{
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		bool whether_same = 1;
		int mosaic = g_MeshList[g_MeshList[i].m_MeshLink[0]].m_Mosaic;
		for (int j = 0; j < g_MeshList[i].m_MeshLink.size(); j++)
		{
			int next = g_MeshList[i].m_MeshLink[j];
			if (g_MeshList[next].m_Mosaic != mosaic)
			{
				whether_same = 0;
				break;
			}
		}
		if (whether_same == 1)
		{
			g_MeshList[i].m_Mosaic = mosaic;
		}
	}
}

struct flow_member
{
	bool m_Type; //positive or negative
	int m_PreviousNumber;
	int m_PreviousNextNumber;
	double m_Value;
	bool m_Visit;
	void m_Set(bool type, int pre, int pre_next, double value) 
	{
		m_Type = type; 
		m_PreviousNumber = pre;
		m_PreviousNextNumber = pre_next; 
		m_Value = value;
		
	}
	flow_member()
	{
		m_Type = 1;
		m_PreviousNumber = -1;
		m_PreviousNextNumber = 0;
		m_Value = 0;
		m_Visit = 0;
	}
	void m_Clear()
	{
		m_Type = 1;
		m_PreviousNumber = -1;
		m_PreviousNextNumber = 0;
		m_Value = 0;
		m_Visit = 0;
	}
};

int depth[3*SIZE] = { 0 };
flow_member temp[3*SIZE];

double TMeshGraph::m_DFS(int u, double dist)
{
	if (u == m_EndNum) return dist;
	for (int i = 0; i < m_MeshGraph[u].m_Next.size(); i++)
	{
		int next = m_MeshGraph[u].m_Next[i];
		if (depth[next] == depth[u] + 1 && depth[u] > 0
			&& m_MeshGraph[u].m_Capacity[i] > m_MeshGraph[u].m_Flow[i])
		{
			double new_dist = m_DFS(next, min(dist,
				m_MeshGraph[u].m_Capacity[i] - m_MeshGraph[u].m_Flow[i]));
			if (new_dist > 0)
			{
				m_MeshGraph[u].m_Flow[i] += new_dist;
				int reverse = -1;
				for (int j = 0; j < m_MeshGraph[next].m_Next.size(); j++)
				{
					if (m_MeshGraph[next].m_Next[j] == u)
					{

						reverse = j;
						break;
					}
				}
				m_MeshGraph[next].m_Flow[reverse] -= new_dist;
				return new_dist;
			}
		}
	}
	return 0;
}

bool TMeshGraph::m_BFS()
{
	//广搜更新层次图层次
	queue<int> que;
	for (int i = 0; i < m_MeshGraph.size(); i++) depth[i] = 0;
	depth[m_SourceNum] = 1; que.push(m_SourceNum);
	while (!que.empty())
	{
		int top = que.front();
		que.pop();
		for (int i = 0; i < m_MeshGraph[top].m_Next.size(); i++)
		{
			int next = m_MeshGraph[top].m_Next[i];
			if (m_MeshGraph[top].m_Capacity[i] > m_MeshGraph[top].m_Flow[i] &&
				depth[next] == 0)
			{
				depth[next] = depth[top] + 1;
				que.push(next);
			}
		}
	}
	if (depth[m_EndNum] > 0) return 1;
	else return 0;
}

//dinic
double TMeshGraph::m_NetWorkFlow()
{
	double ans = 0;
	while (m_BFS())
	{
		while (double d = m_DFS(m_SourceNum, MAX_NUM))
		{
			ans += d;
		}
	}
	return ans;
}




void TMeshGraph::m_BuildGraph(int a)
{
	//build mesh points
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		TMeshNode new_node(i, 1, i);
		m_MeshGraph.push_back(new_node);
	}

	//build source and end
	m_SourceNum = m_MeshGraph.size();
	m_EndNum = m_SourceNum + 1;
	TMeshNode source(m_SourceNum, 0, -1);
	TMeshNode end(m_EndNum, 0, -1);
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		double self = g_GetSelfWeight(g_MeshList[i], g_PictureList[m_CurrentPicture], g_PointList);
		source.m_Next.push_back(i);
		source.m_Flow.push_back(0);
		source.m_Capacity.push_back(self);
		m_MeshGraph[i].m_Previous.push_back(m_SourceNum);
		m_MeshGraph[i].m_PlaceInPrevious.push_back(source.m_Next.size() - 1);

	}
	for (int i = 0; i < m_SourceNum; i++)
	{
		m_MeshGraph[i].m_Next.push_back(m_EndNum);
		m_MeshGraph[i].m_Flow.push_back(0);

		double self = MAX;
		double self_weight = g_GetSelfWeight(g_MeshList[i], g_PictureList[g_MeshList[i].m_Mosaic], g_PointList);
		if (g_MeshList[i].m_Mosaic != m_CurrentPicture)
		{
			self = self_weight;
		}
		m_MeshGraph[i].m_Capacity.push_back(self);
		end.m_Previous.push_back(i);
		end.m_PlaceInPrevious.push_back(m_MeshGraph[i].m_Next.size() - 1);
	}
	m_MeshGraph.push_back(source);
	m_MeshGraph.push_back(end);


	//build appendix points
	for (int i = 0; i < g_MeshList.size(); i++)
	{
		int first_node = i;
		for (int j = 0; j < g_MeshList[i].m_MeshLink.size(); j++)
		{
			int second_node = g_MeshList[i].m_MeshLink[j];
			if (first_node >= second_node)
			{
				continue;
			}
			else if (g_MeshList[first_node].m_Mosaic == g_MeshList[second_node].m_Mosaic)
			{
				// no appendix nodes;
				m_MeshGraph[first_node].m_Next.push_back(second_node);
				m_MeshGraph[first_node].m_Flow.push_back(0);

				double mutual = 0;
				mutual = g_GetMutualWeight(g_MeshList[first_node], g_PictureList[g_MeshList[first_node].m_Mosaic]
					, g_MeshList[second_node], g_PictureList[m_CurrentPicture], g_PointList);
				if (mutual == -1) mutual = MAX;

				m_MeshGraph[first_node].m_Capacity.push_back(mutual);

				m_MeshGraph[second_node].m_Previous.push_back(first_node);
				m_MeshGraph[second_node].m_PlaceInPrevious.push_back(m_MeshGraph[first_node].m_Next.size() - 1);

			}
			else
			{
				//have appendix node
				TMeshNode new_node(m_MeshGraph.size(), 0, -1);

				m_MeshGraph[first_node].m_Next.push_back(m_MeshGraph.size());
				m_MeshGraph[first_node].m_Flow.push_back(0);
				double mutual = 0;

				mutual = g_GetMutualWeight(g_MeshList[first_node], g_PictureList[g_MeshList[first_node].m_Mosaic]
					, g_MeshList[second_node], g_PictureList[m_CurrentPicture], g_PointList);
				if (mutual == -1) mutual = MAX;



				m_MeshGraph[first_node].m_Capacity.push_back(mutual);
				new_node.m_Previous.push_back(first_node);
				new_node.m_PlaceInPrevious.push_back(m_MeshGraph[first_node].m_Next.size() - 1);

				new_node.m_Next.push_back(second_node);
				new_node.m_Flow.push_back(0);

				mutual = g_GetMutualWeight(g_MeshList[first_node], g_PictureList[m_CurrentPicture]
					, g_MeshList[second_node], g_PictureList[g_MeshList[second_node].m_Mosaic], g_PointList);
				if (mutual == -1) mutual = MAX;

				new_node.m_Capacity.push_back(mutual);

				m_MeshGraph[second_node].m_Previous.push_back(m_MeshGraph.size());
				m_MeshGraph[second_node].m_PlaceInPrevious.push_back(new_node.m_Next.size()-1);

				new_node.m_Next.push_back(m_EndNum);
				new_node.m_Flow.push_back(0);

				mutual = g_GetMutualWeight(g_MeshList[first_node], g_PictureList[g_MeshList[first_node].m_Mosaic]
					, g_MeshList[second_node], g_PictureList[g_MeshList[second_node].m_Mosaic], g_PointList);
				if (mutual == -1) mutual = MAX;



				new_node.m_Capacity.push_back(mutual);

				m_MeshGraph[m_EndNum].m_Previous.push_back(m_MeshGraph.size());
				m_MeshGraph[m_EndNum].m_PlaceInPrevious.push_back(new_node.m_Next.size() - 1);


				m_MeshGraph.push_back(new_node);
			}
		}
	}
}


vector<int> A,O;
vector<int> tree_information; //1:s -1:t 0:no
vector<int> previous; //-1:no
queue<int> que_s;
vector<int> route;
vector<int> next_list;
int ans = 0;

//pami04
void TMeshGraph::m_InitializeGraph()
{
	A.clear(); O.clear();
	tree_information.clear();
	previous.clear();
	while (!que_s.empty()) que_s.pop();
	route.clear();
	next_list.clear();
	ans = 0;

	A.push_back(m_SourceNum);
	A.push_back(m_EndNum);
	for (int i = 0; i < m_MeshGraph.size(); i++)
	{
		tree_information.push_back(0);
		previous.push_back(-1);
	}
	tree_information[m_SourceNum] = 1;
	tree_information[m_EndNum] = -1;
}


void TMeshGraph::m_GrowRoute(int& a, int& b)
{
	while (!A.empty())
	{
		int current = A[A.size() - 1];
		if (tree_information[current] == 1)
		{
			for (int i = 0; i < m_MeshGraph[current].m_Next.size(); i++)
			{

				int next = m_MeshGraph[current].m_Next[i];
				double cap = m_MeshGraph[current].m_Capacity[i] - m_MeshGraph[current].m_Flow[i];
				if (cap > 0)
				{
					if (tree_information[next] == 0)
					{
						tree_information[next] = tree_information[current];
						previous[next] = current;
						A.push_back(next);
					}
					else if (tree_information[next] != tree_information[current])
					{
						a = current;
						b = next;
						return;
					}
				}
			}
		}
		else if (tree_information[current] == -1)
		{
			for (int i = 0; i < m_MeshGraph[current].m_Previous.size(); i++)
			{

				int next = m_MeshGraph[current].m_Previous[i];
				double cap = m_MeshGraph[next].m_Capacity[m_MeshGraph[current].m_PlaceInPrevious[i]]
					- m_MeshGraph[next].m_Flow[m_MeshGraph[current].m_PlaceInPrevious[i]];
				if (cap > 0)
				{
					if (tree_information[next] == 0)
					{
						tree_information[next] = tree_information[current];
						previous[next] = current;
						A.push_back(next);
					}
					else if (tree_information[next] != tree_information[current])
					{
						a = current;
						b = next;
						return;
					}
				}
			}
		}
		for (int i = 0; i < A.size(); i++)
		{
			if (A[i] == current)
			{
				A.erase(A.begin() + i);
				break;
			}
		}
	}
	a = -1;
	b = -1;
	return;
}




void TMeshGraph::m_AugmentRoute(int a, int b)
{
	route.clear();
	next_list.clear();
	if (tree_information[a] == -1)
	{
		int temp = a; a = b; b = temp;
	}
	//a:s b::t


	//handle a:s tree
	while (previous[a] >= 0)
	{
		que_s.push(a);
		a = previous[a];
	}
	//if (a == m_SourceNum)
	{
		que_s.push(a);
	}
	/*else
	{
	return;
	}*/

	while (!que_s.empty())
	{
		int temp = que_s.back();
		que_s.pop();
		route.push_back(temp);
	}
	while (previous[b] >= 0)
	{
		route.push_back(b);
		b = previous[b];
	}
	//if (b == m_EndNum)
	{
		route.push_back(b);
	}
	/*else
	{
	return;
	}*/

	double min_increment = MAX_NUM;
	for (int i = 0; i < route.size() - 1; i++)
	{
		for (int j = 0; j < m_MeshGraph[route[i]].m_Next.size(); j++)
		{
			int potential_next = m_MeshGraph[route[i]].m_Next[j];
			if (potential_next == route[i + 1])
			{
				double increment = m_MeshGraph[route[i]].m_Capacity[j]
					- m_MeshGraph[route[i]].m_Flow[j];
				if (increment < min_increment) min_increment = increment;

				next_list.push_back(j);
				break;
			}
		}
	}
	ans += min_increment;
	for (int i = 0; i < route.size() - 1; i++)
	{
		m_MeshGraph[route[i]].m_Flow[next_list[i]] += min_increment;
		if (m_MeshGraph[route[i]].m_Flow[next_list[i]] == m_MeshGraph[route[i]].m_Capacity[next_list[i]])
		{
			if (tree_information[route[i]] == 1 && tree_information[route[i + 1]] == 1)
			{
				previous[route[i + 1]] = -1;
				O.push_back(route[i + 1]);
			}
			else if (tree_information[route[i]] == -1 && tree_information[route[i + 1]] == -1)
			{
				previous[route[i]] = -1;
				O.push_back(route[i]);

			}
		}
	}
}

void TMeshGraph::m_AdoptOrphan()
{
	while (!O.empty())
	{
		int current = O[O.size() - 1];
		O.pop_back();
		bool whether_valid = 0;
		if (tree_information[current] == 1)
		{
			for (int i = 0; i < m_MeshGraph[current].m_Previous.size(); i++)
			{
				int nova = m_MeshGraph[current].m_Previous[i];
				double remain = m_MeshGraph[nova].m_Capacity[m_MeshGraph[current].m_PlaceInPrevious[i]]
					- m_MeshGraph[nova].m_Flow[m_MeshGraph[current].m_PlaceInPrevious[i]];
				if (remain > 0 && tree_information[nova] == 1)
				{
					int judge_nova = nova;
					while (previous[judge_nova] >= 0)
					{
						judge_nova = previous[judge_nova];
					}
					if (judge_nova == m_SourceNum || judge_nova == m_EndNum)
					{
						//valid
						previous[current] = nova;
						whether_valid = 1;
						break;
					}

				}
			}
			if (whether_valid == 0)
			{
				//invalid
				for (int i = 0; i < m_MeshGraph[current].m_Previous.size(); i++)
				{
					int nova = m_MeshGraph[current].m_Previous[i];
					double remain = m_MeshGraph[nova].m_Capacity[m_MeshGraph[current].m_PlaceInPrevious[i]]
						- m_MeshGraph[nova].m_Flow[m_MeshGraph[current].m_PlaceInPrevious[i]];
					if (tree_information[nova] == 1)
					{
						if (remain > 0)
						{
							A.push_back(nova);
						}
						if (previous[nova] == current)
						{
							previous[nova] = -1;
							O.push_back(nova);
						}
					}
				}

				for (int i = 0; i < m_MeshGraph[current].m_Next.size(); i++)
				{
					int nova = m_MeshGraph[current].m_Next[i];
					if (tree_information[nova] == 1)
					{
						if (previous[nova] == current)
						{
							previous[nova] = -1;
							O.push_back(nova);
						}
					}
				}

				tree_information[current] = 0;
				for (int i = 0; i < A.size(); i++)
				{
					if (A[i] == current)
					{
						A.erase(A.begin() + i);
						break;
					}
				}
			}
		}


		else if (tree_information[current] == -1)
		{
			for (int i = 0; i < m_MeshGraph[current].m_Next.size(); i++)
			{
				int nova = m_MeshGraph[current].m_Next[i];
				double remain = m_MeshGraph[current].m_Capacity[i]
					- m_MeshGraph[current].m_Flow[i];
				if (remain > 0 && tree_information[nova] == -1)
				{
					int judge_nova = nova;
					while (previous[judge_nova] >= 0)
					{
						judge_nova = previous[judge_nova];
					}
					if (judge_nova == m_SourceNum || judge_nova == m_EndNum)
					{
						//valid
						previous[current] = nova;
						whether_valid = 1;
						break;
					}

				}
			}
			if (whether_valid == 0)
			{
				//invalid
				for (int i = 0; i < m_MeshGraph[current].m_Next.size(); i++)
				{
					int nova = m_MeshGraph[current].m_Next[i];
					double remain = m_MeshGraph[current].m_Capacity[i]
						- m_MeshGraph[current].m_Flow[i];
					if (tree_information[nova] == -1)
					{
						if (remain > 0)
						{
							A.push_back(nova);
						}
						if (previous[nova] == current)
						{
							previous[nova] = -1;
							O.push_back(nova);
						}
					}
				}
				for (int i = 0; i < m_MeshGraph[current].m_Previous.size(); i++)
				{
					int nova = m_MeshGraph[current].m_Previous[i];
					if (tree_information[nova] == -1)
					{
						if (previous[nova] == current)
						{
							previous[nova] = -1;
							O.push_back(nova);
						}
					}
				}
				tree_information[current] = 0;
				for (int i = 0; i < A.size(); i++)
				{
					if (A[i] == current)
					{
						A.erase(A.begin() + i);
						break;
					}
				}
			}
		}
	}
}

double TMeshGraph::m_NetWorkFlow(int pp)
{
	m_InitializeGraph();
	while (1)
	{
		int a = -1, b = -1;
		m_GrowRoute(a, b);
		if (a == -1) break;
		m_AugmentRoute(a, b);
		m_AdoptOrphan();
	}
	return ans;

}