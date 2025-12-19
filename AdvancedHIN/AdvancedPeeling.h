#ifndef MOTIFMODULARITY_ADVANCEDPEELING_H
#define MOTIFMODULARITY_ADVANCEDPEELING_H

#include<set>
#include<math.h>
#include<algorithm>
#include<queue>
#include <chrono>
#include"../graph/graph.h"
#include "../Result/result.h"
#include "../Config/config.h"
#include "../Match/match.h"
class AdvancedPeeling {
public:
    result *search(set<int> inputQuery, result *matchResult, graph *inputMotif, graph *inputTarget,double inputY,int mod,int optimize);
    void modularityCount(int mod);
    void checkConnectivity(map<int,set<int>>*conNodeGraph,set<int>*nodeList);
    bool deleteNode(int ID, set<int> *peelingNode, set<int> *peelingMotif,map<pair<int,int>,int>*oldEdgeWeight,set<pair<int,int>>*peelingConEdge);
    void returnBaseBack(set<int>*peelingNode,set<int>*peelingMotif,map<pair<int,int>,int>*oldEdgeWeight,set<pair<int,int>>*peelingConEdge);
    void greedyPeelingPlus();
    void computeContribution(int mod);
    void motifConnectedNodeGraph();
    void initialMM();
    graph *motif;
    graph *target;
    double y;
    int wg;
    int m;
    map<int,set<int>>classifyMotif;
    map<int,long int>mmVolC;
    map<int,long int>preMMVolC;
    map<int,long int>mmVolH;
    map<int,int> volC;
    map<int,int> preVolC;
    map<int,int> volH;
    long double MM;
    int pren;
    int prewc;
    int n;
    int wc;
    result *initialRsl;
    result nowRsl;
    result rsl;
    set<int> query;
    map<int,set<int>>dominateList;
    map<int,set<pair<int,int>>>vertexDominate;
    map<int,set<pair<int,int>>>oneHopEdge;
    void computeCutNode();
    void tarjan(int now,int root,int fa);
    int mod=1;
    int optimize=0;
    struct gain
    {
        int id;
        long double nodeGain;
        gain(int a,double b) : id(a), nodeGain(b) {}
        bool operator<(const gain& other) const {
            return nodeGain > other.nodeGain;
        }
    };
    vector<gain>contribution;
    set<int>cutNodes;
    int dfs_clock;
    map<int,int>dfn;
    map<int,int>low;
};
#endif
