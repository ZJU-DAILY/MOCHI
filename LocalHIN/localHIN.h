#ifndef MOTIFMODULARITY_LOCALHIN_H
#define MOTIFMODULARITY_LOCALHIN_H
#define INF std::numeric_limits<int>::max()
#include<set>
#include<math.h>
#include<algorithm>
#include<queue>
#include <chrono>
#include"../graph/graph.h"
#include "../Result/result.h"
#include "../Config/config.h"
#include "../Match/match.h"
class localHIN {
public:
    set<int> search(set<int> inputQuery, result *matchResult, graph *inputMotif, graph *inputTarget,double inputY,int mod,int gainMod,int inputModularityMod);
    void MDMcount();
    void returnNowRsl(result *rsl);
    bool deleteNodeWithDistance(int ID);
    bool deleteLayerWithDistance();
    void computeMotifDistance();
    void sortByMotifDistance();
    void distancePeeling();
    void layerPeeling();
    void motifConnectedNodeGraph();
    double computeNodeGain(int node);
    vector<int>peelingOrder;
    graph *motif;
    graph *target;
    double y;
    int wg;
    map<int,int>volH;
    int m;
    double MM;
    map<int,int> preVolC;
    int pren;
    int prewc;
    map<int,int> volC;
    map<int,int> x;
    int n;
    int wc;
    result *initialRsl;
    result nowRsl;
    result rsl;
    set<int> query;
    map<int, set<int>> motifDominateList;
    map<int,set<int>> conDominateList;
    map<int,set<pair<int,int>>>vertexDominate;
    map<int, double> nodeGain;
    struct cmp {
        bool operator()(const pair<double, int>& a, const pair<double, int>& b) const {
            if (a.first != b.first)
                return a.first > b.first;
            return a.second < b.second;
        }
    };
    set<pair<double,int>,cmp>peelingListGain;
    set<int>peelingList;
    int mod=1;
    int layer=-1;
    map<int, set<int>, greater<int>> sortedVertexDis;
    map<int,int> vertexDis;
    int maxWeight;
    bool generateNewQuery();
    int nodeGainMod;
    int modularityMod;
};


#endif
