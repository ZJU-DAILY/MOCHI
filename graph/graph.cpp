
#include <chrono>
#include "graph.h"

graph::graph() {}

graph::graph(string graphFile, string vertexFile, string edgeFile) {
    graph::graphFile = graphFile;
    graph::vertexFile = vertexFile;
    graph::edgeFile = edgeFile;
    graph::vertexNum = graph::vertexNumReader(vertexFile);
    graph::edgeNum = graph::edgeNumReader(edgeFile);
    graph::vertexType = graph::readVertexTypeWithLabel(vertexFile);
    graph::edgeType = graph::readEdgeType(edgeFile);
    graph::targetGraph = graph::readGraph(graphFile);
}

int graph::vertexNumReader(string vertexFile)
{
    ifstream inf;
    inf.open(vertexFile);
    string s;
    vertexNum = 0;
    while (getline(inf, s)) {
        vertexNum++;
    }
    inf.close();
    return vertexNum;
}

int graph::edgeNumReader(string edgeFile)
{
    ifstream inf;
    inf.open(edgeFile);
    string s;
    edgeNum = 0;
    while (getline(inf, s)) {
        edgeNum++;
    }
    inf.close();
    return edgeNum;
}

vector<vector<int>> graph::readGraph(string graphFile) {
    ifstream inf;
    inf.open(graphFile);
    string s;
    vector<vector<int>> graph(vertexNum);
    vector<set<int>> neighbors(vertexNum);
    while (getline(inf, s)) {
        int vertexId = -1;
        vector<int> vertex;
        set<int> nb;
        stringstream ss(s);
        int num;
        int i = 0;
        while (ss >> num) {
            if (vertexId == -1)
                vertexId = num;
            else {
                if (i % 2 == 0)
                {
                    edges[pair<int, int>(vertexId, vertex[vertex.size() - 1])] = num;
                    edges2[num] = pair<int, int>(vertexId, vertex[vertex.size() - 1]);
                } else
                    nb.insert(num);
                vertex.push_back(num);
            }
            i++;
        }
        graph[vertexId] = vertex;
        neighbors[vertexId] = nb;
    }
    inf.close();
    nbs = neighbors;
    return graph;
}

vector<int> graph::readVertexTypeWithLabel(string vertexFile) {
    ifstream inf;
    inf.open(vertexFile);
    string s;
    vector<int> vertexType(vertexNum);
    vector<string> vertexName(vertexNum);
    int line = 0;
    while (getline(inf, s)) {
        line++;
        int vertexId = -1;
        istringstream iss(s);
        string token;
        int count = 0;
        while (std::getline(iss, token, ';')) {
            if (count == 0) {
                istringstream tp(token);
                string token2;
                while (std::getline(tp, token2, ' ')) {
                    int num = stoi(token2);
                    if (vertexId == -1)
                        vertexId = num;
                    else {
                        vertexType[vertexId] = num;
                    }
                }
                count++;
            } else if (count == 1) {
                vertexName[vertexId] = token;
                count++;
            }
        }
    }
    inf.close();
    graph::vertexName = vertexName;
    return vertexType;
}
vector<int> graph::readVertexType(string vertexFile) {
    ifstream inf;
    inf.open(vertexFile);
    string s;
    vector<int> vertexType(vertexNum);
    while (getline(inf, s)) {
        int vertexId = -1;
        stringstream ss(s);
        int num;
        while (ss >> num) {
            if (vertexId == -1)
                vertexId = num;
            else {
                vertexType[vertexId] = num;
            }
        }
    }
    inf.close();
    return vertexType;
}

vector<int> graph::readEdgeType(string edgeFile) {
    ifstream inf;
    inf.open(edgeFile);
    string s;
    vector<int> edgeType(edgeNum);
    while (getline(inf, s)) {
        int edgeId = -1;
        stringstream ss(s);
        int num;
        while (ss >> num) {
            if (edgeId == -1)
                edgeId = num;
            else {
                edgeType[edgeId] = num;
            }
        }
    }
    inf.close();
    return edgeType;
}

