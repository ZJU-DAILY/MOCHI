#ifndef MOTIFMODULARITY_NAIVEORIGINPEELING_H
#define MOTIFMODULARITY_NAIVEORIGINPEELING_H

#include <chrono>
#include<set>
#include<math.h>
#include<algorithm>
#include<queue>
#include"../graph/graph.h"
#include "../Result/result.h"
#include "../Config/config.h"
#include "../Match/match.h"
//without using motif graph
class NaiveOriginPeeling {
public:
    result* search(set<int> inputQuery, result *matchResult, graph *inputMotif, graph *inputTarget,double inputY);
    void MMcount();
    set<int> checkConnectivity();
    bool deleteNode(int ID,set<int>*peelingNode,set<int>*peelingMotif);
    void returnBaseBack(set<int>*peelingNode,set<int>*peelingMotif);
    void greedyPeelingPlus();
    void findConnectComp();
    graph *motif;
    graph *target;
    double y;
    int wg;
    int m;
    double MM;
    int pren;
    int prewc;
    map<int,int> volC;
    map<int,int> preVolC;
    map<int,int> volH;
    int n;
    int wc;
    result *initialRsl;
    result nowRsl;
    result rsl;
    set<int> query;
    set<int>queryMotif;
    struct gain
    {
        int id;
        double nodeGain;
        gain(int a,double b) : id(a), nodeGain(b) {}
        bool operator<(const gain& other) const {
            return nodeGain > other.nodeGain;
        }
    };
    vector<gain>contribution;
    set<int>cutMotifs;
    set<int>cutNodes;
    int dfs_clock;
    map<int,int>dfn;
    map<int,int>low;
};


#endif
