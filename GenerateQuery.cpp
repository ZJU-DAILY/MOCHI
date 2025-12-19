#include <sys/stat.h>
#include <filesystem>
#include <iostream>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include "Config/config.h"
#include "graph/graph.h"
#include "Util/Util.h"
#include "Match/match.h"
#include "rapidMatch.h"

static string dataFile;
static string motifFile;
static graph *target;
static graph *motif;
static int motifSize;
static int edgeSize = motifSize * (motifSize - 1) / 4 + 1;//至少要包含的边的数量
static int motifNum = 50;
static int querySize = 10;
static vector<int> queryNum = {1, 2, 4, 6, 8};
static map<pair<int, int>, int> S;
static map<int, set<int>> Smap;

void transformMotif() {
    S.clear();
    Smap.clear();
    for (int i = 0; i < motif->nbs.size(); i++) {
        set<int> nbs = motif->nbs[i];
        int typeu = motif->vertexType[i];
        if (Smap.count(typeu) == 0) {
            for (int nb: nbs) {
                int typev = motif->vertexType[nb];
                if (S.count(pair<int, int>(typeu, typev)) == 0) {
                    S[pair<int, int>(typeu, typev)] = 1;
                    Smap[typeu].insert(typev);
                } else
                    S[pair<int, int>(typeu, typev)]++;
            }
        } else//可能motif中同一类型节点存在多个，取约束最大的情况
        {
            map<int, int> tmpS;//临时存储当前节点的约束
            for (int nb: nbs) {
                int typev = motif->vertexType[nb];
                if (tmpS.count(typev) == 0)
                    tmpS[typev] = 1;
                else
                    tmpS[typev]++;
            }
            for (const auto &entry: tmpS)//如果当前节点的约束比之前同类型节点的约束更加严格，覆盖为更严格的约束
            {
                int typev = entry.first;
                int num = entry.second;
                if (S.count(pair<int, int>(typeu, typev)) == 0) {
                    S[pair<int, int>(typeu, typev)] = num;
                    Smap[typeu].insert(typev);
                } else {
                    if (S[pair<int, int>(typeu, typev)] < num)
                        S[pair<int, int>(typeu, typev)] = num;
                }
            }
        }
    }
}

bool checkQualified(int v, set<int> *nbs) {
    if (Smap.count(target->vertexType[v]) == 0)
        return false;//如果节点的类型没有出现在motif中，表明该节点不能成为候选节点
    map<pair<int, int>, int> tmpVdegree;
    for (int nb: *nbs)//遍历节点v的邻居nbs，判断是否满足S
    {
        int typeNb = target->vertexType[nb];
        if (tmpVdegree.count(pair<int, int>(v, typeNb)) == 0)
            tmpVdegree[pair<int, int>(v, typeNb)] = 1;
        else
            tmpVdegree[pair<int, int>(v, typeNb)]++;
    }
    int typev = target->vertexType[v];
    for (int typeu: Smap[typev]) {
        int Snum = S[pair<int, int>(typev, typeu)];
        int vNum = 0;
        if (tmpVdegree.count(pair<int, int>(v, typeu)))
            vNum = tmpVdegree[pair<int, int>(v, typeu)];
        if (vNum < Snum)
            return false;
    }
    return true;
}

void transformMatchGraph(graph *dataGraph, string path, int mod)
{
    ofstream outFile(path);
    if (mod == 1)
        outFile << "t " + to_string(dataGraph->vertexNum) + " " + to_string(dataGraph->edgeNum) << endl;
    else if (mod == 0)
        outFile << "t " + to_string(dataGraph->vertexNum) + " " + to_string(dataGraph->edgeNum / 2) << endl;
    for (int i = 0; i < dataGraph->vertexType.size(); i++)
        outFile << "v " + to_string(i) + " " + to_string(dataGraph->vertexType[i]) + " " +
                   to_string(dataGraph->nbs[i].size()) << endl;
    for (int i = 0; i < dataGraph->nbs.size(); i++) {
        for (int j: dataGraph->nbs[i]) {
            if (i < j)
                outFile << "e " + to_string(i) + " " + to_string(j) << endl;
        }
    }
    outFile.close();
}

static void generateGeneralMotif(){
    config config;
    for (int index = 0; index < sizeof(config.dataFile) / sizeof(config.dataFile[0]); index++)//遍历不同的数据集
    {

        string limitSize = "RandomMotif";
        dataFile = config.dataFile[index];
        string graphFile = config.path + dataFile + "/"  + config.graphFile;
        string vertexFile = config.path + dataFile + "/"  + config.vertexFile;
        string edgeFile = config.path + dataFile + "/"  + config.edgeFile;
        auto start_time = std::chrono::high_resolution_clock::now();//记录起始时间
        graph buildGraph(graphFile, vertexFile, edgeFile);
        target = &buildGraph;
        auto end_time = chrono::high_resolution_clock::now();//记录结束时间
        chrono::duration<double> elapsed_time = end_time - start_time;
        cout << "data reading time:" + to_string(elapsed_time.count()) + "s" << endl;
        string matchGraphPath2 = config.path + dataFile + "/matchGraph.txt";
        //随机生成motif
        srand(time(0));
        Util util;
        util.obtainSchema(target);
//        string randomMotifFile = "../data/" + dataFile + "/motifs";
        string randomMotifFile = "../data/" + dataFile + "/" + limitSize;
        string randomMotifCommand;
        randomMotifCommand = "mkdir -p " + randomMotifFile;
        std::mutex mtx;
        mtx.lock();
//        system(randomMotifCommand.c_str());
        std::filesystem::create_directories(randomMotifFile);
        mtx.unlock();
        for (motifSize = 4; motifSize <= 4; motifSize++) {
            int id = 1;
            vector<int> sameNumber;
            while (id <= motifNum) {
                //随机生成边的数量
                cout << "id: " + to_string(id) << endl;
                int lowerBd;
                int upperBd;
                lowerBd = 0;
                upperBd = motifSize * (motifSize - 1) / 2;//平均边数量的上界
                edgeSize = rand() % (upperBd - lowerBd + 1) + lowerBd;
                //生成motif
                graph mtf = util.generateMotif(motifSize, edgeSize, target);
                if (mtf.vertexType.empty())//当找不到motif时重新生成随机数
                    continue;
                motifFile = limitSize + "/" + to_string(motifSize) + "_" + to_string(id);
                //创建motif文件夹
                string motifPath = "../data/" + dataFile + "/" + motifFile;
                string command;
                command = "mkdir -p " + motifPath;
                mtx.lock();
                std::filesystem::create_directories(motifPath);
                mtx.unlock();
                //写入graph文件
                string motifGraphPath = "../data/" + dataFile + "/" + motifFile + "/graph.txt";
                ofstream motifGraphFile(motifGraphPath, std::ios::out | std::ios::trunc);
                for (int i = 0; i < mtf.targetGraph.size(); i++) {
                    motifGraphFile << to_string(i);
                    for (int nb: mtf.targetGraph[i])
                        motifGraphFile << " " + to_string(nb);
                    motifGraphFile << endl;
                }
                motifGraphFile.close();
                //写入vertex文件
                string motifVertexPath = "../data/" + dataFile + "/" + motifFile + "/vertex.txt";
                ofstream motifVertexFile(motifVertexPath, std::ios::out | std::ios::trunc);
                for (int i = 0; i < mtf.vertexType.size(); i++) {
                    motifVertexFile << to_string(i) + " " + to_string(mtf.vertexType[i]) << endl;
                }
                motifVertexFile.close();
                //写入edge文件
                string motifEdgePath = "../data/" + dataFile + "/" + motifFile + "/edge.txt";
                ofstream motifEdgeFile(motifEdgePath, std::ios::out | std::ios::trunc);
                for (int i = 0; i < mtf.edgeType.size(); i++) {
                    motifEdgeFile << to_string(i) + " " + to_string(mtf.edgeType[i]) << endl;
                }
                motifEdgeFile.close();

                cout << "begin global match" << endl;
                graph buildMotifGraph(motifGraphPath, motifVertexPath, motifEdgePath);
                motif = &buildMotifGraph;
                string queryGraphPath = config.path + dataFile + "/" + motifFile + "/queryGraph.txt";
                transformMatchGraph(motif, queryGraphPath, 0);
                start_time = std::chrono::high_resolution_clock::now();
                match match;
                map<int, set<int>> *symmetricNodes;
                symmetricNodes = match.symmetricDetect(motif);

                rapidMatch rapidMatch1;
                result *matchRsl = new result;
                bool judge1 = rapidMatch1.globalMatch(queryGraphPath, matchGraphPath2, matchRsl, symmetricNodes);
                end_time = std::chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed_time = end_time - start_time;
                cout << "global match time: " + to_string(elapsed_time.count()) + "s" << endl;
                double matchTime = elapsed_time.count() - rapidMatch1.loadGraphTime;
                int matchNum = matchRsl->vertexToMotif.size();
                int vertexNum = matchRsl->vertexToMotifID.size();
                if (!judge1)
                    continue;
                cout << "end global match" << endl;
                //生成motif连通图
                start_time = std::chrono::high_resolution_clock::now();
                match.generateConNodeGraph(matchRsl);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed_time = end_time - start_time;
                double generateConTime = elapsed_time.count();
                cout << "generate conNodeGraph time: " + to_string(elapsed_time.count()) + "s" << endl;
                //记录全图motif总数以及匹配的实例数量
                cout << "begin write info.txt" << endl;
                string infoPath = "../data/" + dataFile + "/" + motifFile + "/info.txt";
                ofstream infoFile(infoPath, std::ios::out | std::ios::trunc);
                infoFile << "global match time: " + to_string(matchTime) + " s" << endl;
                infoFile << "generate conNodeGraph time: " + to_string(generateConTime) + " s" << endl;
                infoFile << "match num: " + to_string(matchNum) << endl;
                infoFile << "vertex with embedding num: " + to_string(vertexNum) << endl;
                infoFile.close();
                cout << "end write info.txt" << endl;
                //生成大小分别为1，2，4，6，8的查询点集合
                for (int i: queryNum) {
                    //创建查询点集合文件夹
                    string queryFile = "../data/" + dataFile + "/" + motifFile + "/" + to_string(i);
                    string queryFileCommand;
                    queryFileCommand = "mkdir -p " + queryFile;
                    mtx.lock();
                    std::filesystem::create_directories(queryFile);
                    mtx.unlock();
                }
                set<int> anchorNode;
                bool flag = true;
                transformMotif();
                for (int i = 1; i <= querySize; i++)//随机挑选查询点
                {
                    if (anchorNode.size() == matchRsl->conNodeGraph.size()) {
                        flag = false;
                        break;
                    }
                    int index = rand() % matchRsl->conNodeGraph.size();
                    auto it = matchRsl->conNodeGraph.begin();
                    advance(it, index);
                    int q = (*it).first;//以该节点作为起始点，找到连通图
                    if (anchorNode.count(q))//避免多次采到同一个锚点
                    {
                        i--;
                        continue;
                    }
                    if (!checkQualified(q, &target->nbs[q])) {
                        anchorNode.insert(q);
                        i--;
                        continue;
                    }
                    anchorNode.insert(q);
                    struct point {
                        int nodeID;
                        int dis;
                    };
                    map<int, int> vertexDis;
                    map<int, set<int>> disReverse;
                    point pt;
                    pt.nodeID = q;
                    pt.dis = 0;
                    queue<point> queueList;
                    vector<int> visNode(target->vertexNum, -1);
                    queueList.push(pt);//选取第一个查询点加入队列
                    visNode[q] = 1;
                    while (!queueList.empty()) {
                        point nd = queueList.front();
                        int nowNode = nd.nodeID;
                        int dis = nd.dis;
                        vertexDis[nowNode] = dis;
                        disReverse[dis].insert(nowNode);
                        queueList.pop();
                        for (int nxtNode: (matchRsl->conNodeGraph)[nowNode]) {
                            if (visNode[nxtNode] == -1) {
                                visNode[nxtNode] = 1;
                                point nxtNd;
                                nxtNd.nodeID = nxtNode;
                                nxtNd.dis = dis + 1;
                                queueList.push(nxtNd);
                            }
                        }
                    }
                    cout << "query set size: " + to_string(i) << endl;
                    set<int> nodeList;
                    int dis = 0;
                    while (nodeList.size() <= 8 && dis <= 2 && !disReverse[dis].empty())//最多取两跳内的节点
                    {
                        for (int node: disReverse[dis]) {
                            if (checkQualified(node, &target->nbs[node]))
                                nodeList.insert(node);
                        }
                        dis++;
                    }
                    if (nodeList.size() < 8)//重新找一个锚点
                    {
                        i--;
                        continue;
                    }
                    //生成的查询点嵌套
                    set<int> nodeSet;
                    for (int j: queryNum) {
                        int s = nodeSet.size();
                        while (s < j) {
                            int nodeIndex = rand() % nodeList.size();
                            auto nodeIt = nodeList.begin();
                            advance(nodeIt, nodeIndex);
                            int id = *nodeIt;
                            if (nodeSet.count(id) != 0)
                                continue;
                            nodeSet.insert(id);
                            s++;
                        }
                        string path =
                                "../data/" + dataFile + "/" + motifFile + "/" + to_string(j) + "/" + to_string(i) +
                                ".txt";
                        ofstream inputQuery(path);
                        for (int q: nodeSet)
                            inputQuery << q << endl;
                        inputQuery.close();
                    }
                }
                if (!flag) {
                    cout << "can't find enough queries" << endl;
                    continue;
                }
                sameNumber.push_back(matchNum);
                //在每个motif文件夹下创建result文件夹
                string resultFile = "../data/" + dataFile + "/" + motifFile + "/result";
                string resultCommand;
                resultCommand = "mkdir -p " + resultFile;
                mtx.lock();
                std::filesystem::create_directories(resultFile);
                mtx.unlock();
                id++;
            }
            int maxSame = 0;
            for (int i = 0; i < sameNumber.size(); i++) {
                int same = 0;
                for (int j = 0; j < sameNumber.size(); j++) {
                    if (sameNumber[i] == sameNumber[j] && i != j)
                        same++;
                }
                if (same > maxSame)
                    maxSame = same;
            }
            string infoPath = "../data/" + dataFile + "/" + limitSize + "/sameInfo.txt";
            ofstream infoFile(infoPath, std::ios::out | std::ios::trunc);
            infoFile << "same motif num: " + to_string(maxSame) << endl;
            infoFile.close();
        }
    }
}

int main() {
    generateGeneralMotif();
    return 0;
}
