#ifndef MOTIFMODULARITY_MATCH_H
#define MOTIFMODULARITY_MATCH_H
#include<vector>
#include<map>
#include<set>
#include<utility>
#include <chrono>
#include <queue>
#include <random>
#include <unordered_set>
#include <algorithm>
#include"../graph/graph.h"
#include "../Result/result.h"

class match {
public:
    vector<map<int,int>>result1;
    vector<map<int,int>>result2;
    result currentRsl;
    set<int>newNodes;
    vector<int> enumerated;
    graph *motif;
    graph *target;
    bool isGlobal;
    result *globalMatch(graph *inputMotif, graph *inputTarget);
    result *localMatch(set<int>query,graph *inputMotif, graph *inputTarget);
    set<int>getR2(int nxtu, map<int, int> M1, map<int, int> M2, map<int,int> mcm);
    void matchByDFS(map<int, int> M1, map<int, int> M2, map<int,int> mcm,set<int>extendableNodes);
    map<int, set<int>>* symmetricDetect(graph *motif);
    graph dataGraph;
    graph queryGraph;
    void isIsomorphism(map<int, int> M1, map<int, int> M2,set<int>extendableNodes);
    set<int>getCandidate(int nxtu, map<int, int> M1, map<int, int> M2);
    bool symmetricFlag;
    map<int, set<int>> symmetricNodes;
    void generateMotifGraph(result *rsl);
    int generateConNodeGraph(result *rsl);
    std::chrono::system_clock::time_point match_start;
    std::chrono::system_clock::time_point match_end;
    chrono::duration<double> match_gap;
};


#endif //MOTIFMODULARITY_MATCH_H
