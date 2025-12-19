#ifndef MOTIFMODULARITY_CONFIG_H
#define MOTIFMODULARITY_CONFIG_H
#include <string>
#include <set>
#include <vector>

using namespace std;
class config {
public:
    string dataFile[1] = {"example"};
    string motifFile="motif";
    //default parameters
    double y=0.5;
    int matchTimeLimit=10000;
    int timeLimit=10000;
    int tryLimit=100;
    int nbLimit=10000;
    int motifNumLimit_LB = 0;
    long long motifNumLimit=1000000;
    int diameterTime=200;
    string path="../data/";
    string graphFile = "graph.txt";
    string vertexFile = "vertex.txt";
    string edgeFile = "edge.txt";
    string motifGraphFile = "graph.txt";
    string motifVertexFile = "vertex.txt";
    string motifEdgeFile = "edge.txt";
};


#endif //MOTIFMODULARITY_CONFIG_H
