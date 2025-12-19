#ifndef MOTIFMODULARITY_RESULT_H
#define MOTIFMODULARITY_RESULT_H


#include "../graph/graph.h"
#include<set>
class result {
public:
    set<int> nodeGraph;//当前图内的节点
    map<int,set<int>> motifGraph;//以motif为单位的图
    map<int,set<int>>conNodeGraph;//满足motif-connected的连通图
    map<pair<int,int>,int>conEdgeWeight;//存储节点之间边的权重
    map<int,set<int>>originGraph;
    vector<vector<int>>matched_subgraph;
//    map<pair<int,int>,int> motifEdges;
//    map<int,pair<int,int>> motifEdges2;
//    vector<map<int, int>> motifToVertex;//vector下标为motifID，<key,value>,key为motif对应的节点，value为匹配的节点
    vector<map<int, int>> vertexToMotif;//vector下标为motifID，<key,value>,key为匹配的节点，value为motif对应的节点
    map<int,set<int>>vertexToMotifID;//存储每个节点所对应的motifID
    double modularity;//存储当前图下的modularity
    int motifNum;//存储motif的数量
    map<int,int> vertexDis;//每个节点到查询点的motif距离
    double motifDensity;
    double density;
    double motifDiameter;
    double diameter;
    double avgJaccardSim;
    map<int,double>avgEuclideanSim;
    vector<int>nodeGraphVis;
    set<int>motifSet;//当前图内的motif
    bool outLimitFlag=false;//标记运行是否超时
    double queryMotifDistance;//到查询点的motif距离
    double queryDistance;//到查询点的距离
    double matchTime;
    double executeTime;
    double totalTime;
    double motifInstance;
    double avgMotifJaccardSim;
    double avgMotifComSim;
    double avgMotifRatio;
    double motifDensityRate;
    double communitySize;
    double minMaxMotifSim;
    double motifConductance;
    double avgMotifSimilarity;
    double avgPairMotifSimilarity;
    double conductance;
    double expansion;
    double motifExpansion;
    map<int,int>volH;//HIN中不同类型节点volume
    int wg;
};


#endif //MOTIFMODULARITY_RESULT_H
