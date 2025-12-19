#ifndef MOTIFMODULARITY_GRAPH_H
#define MOTIFMODULARITY_GRAPH_H

#include <iostream>
#include<cstring>
#include<vector>
#include<fstream>
#include <sstream>
#include<map>
#include<set>
using namespace std;
class graph {
public:
    string graphFile;
    string vertexFile;
    string edgeFile;
    vector<vector<int>>targetGraph;
    vector<set<int>>nbs;
    vector<int>vertexType;
    vector<int>edgeType;
    map<pair<int,int>,int> edges;
    map<int,pair<int,int>> edges2;
    int vertexNum;
    int edgeNum;
    vector<string>vertexName;
    graph();
    graph(string graphFile,string vertexFile,string edgeFile);
    int vertexNumReader(string vertexFile);
    int edgeNumReader(string edgeFile);
    vector<vector<int>>readGraph(string graphFile);
    vector<int>readVertexType(string vertexFile);
    vector<int>readVertexTypeWithLabel(string vertexFile);
    vector<int>readEdgeType(string edgeFile);
};


#endif
