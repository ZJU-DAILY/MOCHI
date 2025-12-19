//
// Created by 周于涵 on 2024-01-04.
//

#ifndef MOCHI_UTIL_H
#define MOCHI_UTIL_H
#define fi first
#define se second
#define mp make_pair

#include "../Result/result.h"
#include<cstdio>
#include<queue>
#include<cstring>
#include<functional>
#include <random>
#include <iterator>
#include<stack>
class Util {
public:
    map<int,graph> generateExpMotif(graph *target);
    void experimentRsl(result *rsl, result *initialRsl,graph *target,set<int> *query);
    double motifDensity(result *rsl);
    double motifSimilarity(result *rsl);
    double density(result *rsl,graph *target);
    double minMaxMotifJarccard(result *rsl,result *initialRsl);
    int BFS(map<int,set<int>> *graph, int start,int maxID);
    int graphDiameter(map<int,set<int>> *graph);//compute diameter
    int queryDistance(map<int,set<int>> *graph,set<int>*query);//计算到查询点的距离
    double computeJaccard(set<int>s1,set<int>s2);
    double avgMotifJaccardSim(result *rsl);
    double avgJaccardSim(result *rsl, graph *target);
    map<int,double> avgEuclideanSim(result *rsl,graph *target);
    void obtainSchema(graph *target);//遍历图，获取图的schema
    double motifRatio(result *rsl,result *initialRsl,graph *target);
    double motifDensityRate(result *rsl,result *initialRsl, graph *target);
    double pairMotifSimilarity(result *rsl);
    result initialRsl;
    map<int,set<int>>schema;
    map<pair<int,int>,double>schemaEdge;//记录存在连边的两种类型节点的实际边数量
    map<pair<int,int>,int>schemaEdgeType;//记录schema每条边的类型
    set<int>readQuery(string queryFile);//读取查询点文件
    set<int>generateTopQuery(int size,string dataFile,string motifFile,result *matchRsl,graph *target,graph *motif);
    graph generateMotif(int motifSize,int edgeSize,graph *target);//直接在原图上随机游走生成motif
    double avgMotifComSim(result *rsl,result *initialRsl);
    graph generateMotifByLabel(int motifSize,int edgeSize,graph *target,set<int>labelSet);
    graph generateMotifByType(int motifSize,int edgeSize,graph *target, int type);
    graph *target;
};


#endif
