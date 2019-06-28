/*
filename: TExpansion.h
description: used in doing alpha-expension
date: 2/17/2019
*/

#ifndef TEXPENSION_H
#define TEXPENSION_H



//system requirements
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

#define MAX_NUM 1145141919810


//cv requirements
#include <opencv2/opencv.hpp>

//algebra requirements
#include <Eigen/Eigen>
#include <limits>

//personal requirements
#include "TDefine.hpp"
#include "TInput.h"
#include "TMath.h"


struct tree
{
	vector<int> list;
	int start_tag_in_graph;
	void clear()
	{
		list.clear();
	}
	void push_back(int a)
	{
		list.push_back(a);
	}
};



class TMeshNode
{
public:
	int m_Number;
	bool m_Type;//1:mesh  0:others
	int m_MeshNumber; //non-mesh:-1
	vector<int> m_Next;
	vector<double> m_Capacity;
	vector<double> m_Flow;
	vector<int> m_Previous;
	vector<int> m_PlaceInPrevious;
public:
	TMeshNode(int num, bool type, int mesh_num)
	{
		m_Number = num;
		m_Type = type;
		if (type == 0)
		{
			m_MeshNumber = -1;
		}
		else
		{
			m_MeshNumber = mesh_num;
		}
	}
};


class TMeshGraph
{
public:
	vector<TMeshNode> m_MeshGraph; 
	int m_SourceNum, m_EndNum;
	int m_CurrentPicture;
	double min_energy{ MAX_NUM };
public:
	void m_Set()
	{
		int n, m, s, t;
		cin >> n >> m >> s >> t;
		m_SourceNum = s - 1; m_EndNum = t - 1;
		for (int i = 0; i < n; i++)
		{
			TMeshNode new_node(i, 1, i);
			m_MeshGraph.push_back(new_node);
		}
		//TMeshNode source(m_SourceNum, 0, 0);
		//TMeshNode end(m_EndNum, 0, 0);
		//m_MeshGraph.push_back(source);
		//m_MeshGraph.push_back(end);
		for (int i = 1; i <= m; i++)
		{
			int begin, end;
			double value;
			cin >> begin >> end >> value;
			m_MeshGraph[begin - 1].m_Next.push_back(end - 1);
			m_MeshGraph[begin - 1].m_Flow.push_back(0);
			m_MeshGraph[begin - 1].m_Capacity.push_back(value);

			m_MeshGraph[end - 1].m_Previous.push_back(begin - 1);
			m_MeshGraph[end - 1].m_PlaceInPrevious.push_back(m_MeshGraph[begin - 1].m_Next.size()-1);
		}
	}



	void m_BuildGraph();

	//建图
	void m_BuildGraph(int a);

	//主函数
	void m_MainExpansion();


	void m_Purify();
	bool m_BFS(); //return 1:can use flow
	double m_DFS(int u, double dist);

	//dinic
	double m_NetWorkFlow();
	
	//demi 04，实际用的网络流
	double m_NetWorkFlow(int a);
	void m_InitializeGraph();

	//if has route:return the last search two nodes
	//else :tag = -1
	void m_GrowRoute(int& a, int& b);

	//augment the route
	void m_AugmentRoute(int a, int b);

	//adopt orphans
	void m_AdoptOrphan();
};




#endif