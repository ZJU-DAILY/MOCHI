
#include "Util.h"
#include "../Match/match.h"
#include "../Config/config.h"
#define INF 0x3f3f3f3f
using namespace std;


double Util::motifDensity(result *rsl) {
    double motifDensity;
    int motifNum = rsl->motifSet.size();
    int nodeNum = rsl->nodeGraph.size();
    motifDensity = (1.0) * motifNum / nodeNum;
    return motifDensity;
}


int Util::BFS(map<int, set<int>> *graph, int start,int maxID) {
    struct nodePoint
    {
        int id;
        int dist;
    };
    queue<nodePoint> q; // 存储节点及其距离
    vector<int> visited(maxID+1,-1); // 记录节点是否已经被访问过
    q.push({start, 0});
    visited[start]=1;
    int max_distance = 0;
    while (!q.empty()) {
        int node = q.front().id;
        int distance = q.front().dist;
        q.pop();
        max_distance = distance;
        // 遍历当前节点的相邻节点
        for (auto neighbor: graph->at(node)) {
            if (visited[neighbor]==-1) {
                q.push({neighbor, distance + 1});
                visited[neighbor]=1;
            }
        }
    }
    return max_distance;
}

// 计算图的直径
int Util::graphDiameter(map<int, set<int>> *graph) {
    int diameter = 0;
    config config;
    auto start_time = std::chrono::high_resolution_clock::now();//记录起始时间
    int maxID=-1;
    for(const auto&entry:*graph)
    {
        if(entry.first>maxID)
            maxID=entry.first;
    }
    if(maxID==-1)
        return -1;
    int count=0;
    for (const auto &entry: *graph) {
        count++;
        auto end_time = std::chrono::high_resolution_clock::now();//记录起始时间
        chrono::duration<double> elapsed_time = end_time - start_time;
        if(elapsed_time.count()>config.diameterTime)
            return diameter;
        diameter = max(diameter, BFS(graph, entry.first,maxID));
        cout<<to_string(count)+"/"+ to_string(graph->size())<<endl;
    }
    return diameter;
}


double Util::motifSimilarity(result *rsl)//两两节点间的motif交集/motif并集
{
    double sumSim=0.0;
    for(int u:rsl->nodeGraph)
    {
        for(int v:rsl->nodeGraph)
        {
            if(u==v and !rsl->vertexToMotifID[u].empty())
            {
                sumSim += 1.0;
                continue;
            }
            if(u==v and rsl->vertexToMotifID[u].empty())
            {
                sumSim += 0;
                continue;
            }
           if(u!=v and rsl->conEdgeWeight.count(pair<int,int>(u,v))==0)
           {
               sumSim += 0;
               continue;
           }
            int wg=rsl->conEdgeWeight[pair<int,int>(u,v)];
            int du=rsl->vertexToMotifID[u].size();
            int dv=rsl->vertexToMotifID[v].size();
            sumSim+=(1.0)*wg/(du+dv-wg);
        }
    }
    double sim=sumSim/(rsl->nodeGraph.size()*rsl->nodeGraph.size());
    return sim;
}

double Util::pairMotifSimilarity(result *rsl)//存在公共motif节点的平均相似度
{
    double sumSim=0.0;
    double pairNum=0.0;
    for(int u:rsl->nodeGraph)
    {
        for(int v:rsl->nodeGraph)
        {
            if(u!=v and rsl->conEdgeWeight.count(pair<int,int>(u,v)))
            {
                int wg=rsl->conEdgeWeight[pair<int,int>(u,v)];
                int du=rsl->vertexToMotifID[u].size();
                int dv=rsl->vertexToMotifID[v].size();
                sumSim+=(1.0)*wg/(du+dv-wg);
                pairNum++;
            }
        }
    }
    double sim;
    if(pairNum == 0)
        sim = 0.0;
    else
        sim=sumSim/(pairNum);
    return sim;
}

double Util::avgJaccardSim(result *rsl, graph *target) {
    double avgJaccardSim;
    map<int,map<int,int>>vertexTypeCount;
    //STEP1:compute multi-set of vertex
    map<int, set<int>> nodeLabels;
    for (int node: rsl->nodeGraph) {
        int tp=target->vertexType[node];
        set<int> nbs = target->nbs[node];
        for (int nb: nbs) {
            int nbtp = target->vertexType[nb];
            if (rsl->nodeGraph.count(nb))
            {
                if(vertexTypeCount.count(node)==0)
                {
                    map<int,int>newMap;
                    vertexTypeCount[node] = newMap;
                }
                if(vertexTypeCount[node].count(nbtp)==0)
                    vertexTypeCount[node][nbtp] = 1;
                else
                    vertexTypeCount[node][nbtp]++;
            }
        }
    }
    //STEP2: compute similarity
    double totalJaccardSim = 0.0;
    int pairCount = 0;
    vector<int> nodes;
    for (const auto& pair : vertexTypeCount) {
        nodes.push_back(pair.first);
    }
    for (size_t i = 0; i < nodes.size(); ++i) {
        for (size_t j = i + 1; j < nodes.size(); ++j) {
            int u = nodes[i];
            int v = nodes[j];
            if(target->vertexType[u]!=target->vertexType[v])
                continue;

            const auto &uCounts = vertexTypeCount[u];
            const auto &vCounts = vertexTypeCount[v];

            int intersectionSum = 0;
            int unionSum = 0;

            for (const auto &[type, uCount]: uCounts) {
                int vCount = vCounts.count(type) ? vCounts.at(type) : 0;
                intersectionSum += min(uCount, vCount);
                unionSum += max(uCount, vCount);
            }

            for (const auto &[type, vCount]: vCounts) {
                if (uCounts.count(type) == 0) {
                    unionSum += vCount;
                }
            }

            double jaccardSim = (unionSum > 0) ? (double) intersectionSum / unionSum : 0.0;

            totalJaccardSim += jaccardSim;
            ++pairCount;
        }
    }

    return (pairCount > 0) ? totalJaccardSim / pairCount : 0.0;
}

void Util::experimentRsl(result *rsl, result *initialRsl,graph *inputTarget,set<int> *query) {
    auto start_time = std::chrono::high_resolution_clock::now();
    target=inputTarget;
    rsl->modularity=rsl->modularity;
    rsl->motifDensity = motifDensity(rsl);
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_time = end_time - start_time;
    cout << "compute motif density time:" + to_string(elapsed_time.count()) + "s" << endl;

    start_time = std::chrono::high_resolution_clock::now();
    rsl->avgMotifSimilarity = motifSimilarity(rsl);
    end_time = chrono::high_resolution_clock::now();
    elapsed_time = end_time - start_time;
    cout << "compute avgMotifSimilarity:" + to_string(elapsed_time.count()) + "s" << endl;

    start_time = std::chrono::high_resolution_clock::now();
    rsl->avgPairMotifSimilarity = pairMotifSimilarity(rsl);
    end_time = chrono::high_resolution_clock::now();
    elapsed_time = end_time - start_time;
    cout << "compute avgPairMotifSimilarity:" + to_string(elapsed_time.count()) + "s" << endl;

    start_time = std::chrono::high_resolution_clock::now();
    rsl->avgJaccardSim = avgJaccardSim(rsl, target);
    end_time = chrono::high_resolution_clock::now();
    elapsed_time = end_time - start_time;
    cout << "compute avgJaccardSim time:" + to_string(elapsed_time.count()) + "s" << endl;

    start_time = std::chrono::high_resolution_clock::now();
    map<int, set<int>> rslGraph;
    for (int u: rsl->nodeGraph) {
        set<int> nbs = target->nbs[u];
        for (int v: nbs) {
            if (rsl->nodeGraph.count(v) > 0)
                rslGraph[u].insert(v);
        }
    }
    rsl->diameter = graphDiameter(&rslGraph);
    end_time = chrono::high_resolution_clock::now();
    elapsed_time = end_time - start_time;
    cout << "compute diameter time:" + to_string(elapsed_time.count()) + "s" << endl;
}

void Util::obtainSchema(graph *target) {
    for (int i = 0; i < target->nbs.size(); i++) {
        set<int> nb = target->nbs[i];
        int typeu = target->vertexType[i];
        for (int node: nb) {
            int typev = target->vertexType[node];
            if (schemaEdge.count(pair<int, int>(typeu, typev)) == 0) {
                schema[typeu].insert(typev);
                schema[typev].insert(typeu);
                int eid1 = target->edges[pair<int, int>(i, node)];
                int eid2 = target->edges[pair<int, int>(node, i)];
                schemaEdgeType[pair<int, int>(typeu, typev)] = target->edgeType[eid1];
                schemaEdgeType[pair<int, int>(typev, typeu)] = target->edgeType[eid2];
                schemaEdge[pair<int, int>(typeu, typev)] = 1;
                schemaEdge[pair<int, int>(typev, typeu)] = 1;
            } else {
                schemaEdge[pair<int, int>(typeu, typev)] = schemaEdge[pair<int, int>(typeu, typev)] + 1;
                schemaEdge[pair<int, int>(typev, typeu)] = schemaEdge[pair<int, int>(typev, typeu)] + 1;
            }
        }
    }
}

set<int>Util::readQuery(string queryFile)
{
    set<int>query;
    ifstream inf;
    inf.open(queryFile);
    string s;
    while (getline(inf, s)) {
        if(s.length()!=0)
        {
            int num= stoi(s);
            query.insert(num);
        }
    }
    inf.close();
    return query;
}
graph Util::generateMotif(int motifSize,int edgeSize,graph *target)//size要求>=2
{
//    obtainSchema(target);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, target->edgeNum-1);
    bool flag=true;
    graph motif;
    int tryTime=0;//记录尝试寻找motif的次数
    config config;
    while(flag)
    {
        tryTime++;
        cout<<"tryTime:"+ to_string(tryTime)<<endl;
        if(tryTime>config.tryLimit)//当尝试次数超出上限时返回空的motif
        {
            graph newMotif;
            return newMotif;
        }
        flag=false;
        graph newMotif;
        motif=newMotif;
        vector<int>vis;
        for(int i=0;i<target->edgeNum;i++)
            vis.push_back(-1);
        set<int>extendableEdges;
        int currentEdge = dis(gen);
        extendableEdges.insert(currentEdge);
        vis[currentEdge]=1;
        int u = target->edges2[currentEdge].first;
        int v = target->edges2[currentEdge].second;
        int symmetricEdge = target->edges[pair<int, int>(v, u)];//找到对称的那条边
        vis[symmetricEdge]=1;
        srand(time(0));
        map<int,int>oldToNew;//建立旧点和新点之间的映射
        while(motif.vertexType.size()<motifSize||motif.edgeType.size()<edgeSize*2)//由于一次只会选择一条可扩展的边，因此不会出现边上的两个节点都不在当前图内
        {
            int candidateSize = extendableEdges.size();
//            cout<<"extendableEdge size: "+ to_string(candidateSize)+"  vertexSize: "+ to_string(motif.vertexType.size())+"  edgeSize: "+
//                                                                                                                         to_string(motif.edgeType.size())<<endl;
            auto it=extendableEdges.begin();
            int index=rand() % candidateSize;
            advance(it,index);
            currentEdge = *it;//随机选择一条边进行连接
            it=extendableEdges.erase(it);
            int u = target->edges2[currentEdge].first;
            int v = target->edges2[currentEdge].second;
            if(motif.vertexType.size()==motifSize)//当节点数量达到要求时，接下来只能选择当前图内的边进行连接
            {
                if(oldToNew.count(u)==0||oldToNew.count(v)==0)//如果节点数量达到要求但当前边有一个节点不在图内，表明会引入新的点，跳过
                {
                    if(extendableEdges.empty())//若当前图没有满足的连边时break，重新选点
                    {
                        flag=true;
                        break;
                    }
                    else
                        continue;
                }
            }
            bool flag1;
            bool flag2;
            int oldEid1 = currentEdge;
            int oldEid2 = target->edges[pair<int, int>(v, u)];//找到对称的那条边
            int newu;
            int newv;
            if (oldToNew.count(u)) {
                newu = oldToNew[u];
                flag1=true;
            } else {
                flag1=false;
                newu = motif.vertexType.size();
                motif.vertexType.push_back(target->vertexType[u]);
                oldToNew[u] = newu;
                vector<int>vct;
                motif.targetGraph.push_back(vct);
                set<int>t;
                motif.nbs.push_back(t);
            }
            if (oldToNew.count(v))
            {
                newv = oldToNew[v];
                flag2= true;
            }
            else {
                flag2=false;
                newv = motif.vertexType.size();
                motif.vertexType.push_back(target->vertexType[v]);
                oldToNew[v] = newv;
                vector<int>vct;
                motif.targetGraph.push_back(vct);
                set<int>t;
                motif.nbs.push_back(t);
            }
            int newEid1=motif.edgeType.size();
            motif.edgeType.push_back(target->edgeType[oldEid1]);
            int newEid2=motif.edgeType.size();
            motif.edgeType.push_back(target->edgeType[oldEid2]);
            int tpu=motif.vertexType[newu];
            int tpv=motif.vertexType[newv];
            int tpE1=schemaEdgeType[pair<int,int>(tpu,tpv)];
            motif.targetGraph[newu].push_back(newv);
            if(motif.nbs[newu].count(newv))
                cout<<""<<endl;
            motif.nbs[newu].insert(newv);
            motif.targetGraph[newu].push_back(newEid1);
            if(tpE1!=motif.edgeType[newEid1])
                cout<<""<<endl;
            motif.targetGraph[newv].push_back(newu);
            if(motif.nbs[newv].count(newu))
                cout<<""<<endl;
            motif.nbs[newv].insert(newu);
            motif.targetGraph[newv].push_back(newEid2);
            int tpE2=schemaEdgeType[pair<int,int>(tpv,tpu)];
            if(tpE2!=motif.edgeType[newEid2])
                cout<<""<<endl;
            if(!flag1&&target->nbs[u].size()<config.nbLimit)
            {
                for (int node: target->nbs[u]) {
                    int eid1 = target->edges[pair<int, int>(u, node)];
                    int eid2 = target->edges[pair<int, int>(node, u)];
                    //将未访问过的边加入集合
                    if (vis[eid1] == -1) {
                        vis[eid1] = 1;
                        vis[eid2] = 1;
                        extendableEdges.insert(eid1);//选一个方向的边加入数组即可
                    }
                }
            }
            if(!flag2&&target->nbs[v].size()<config.nbLimit)
            {
                for (int node: target->nbs[v]) {
                    int eid1 = target->edges[pair<int, int>(v, node)];
                    int eid2 = target->edges[pair<int, int>(node, v)];
                    //将未访问过的边加入集合
                    if (vis[eid1] == -1) {
                        vis[eid1] = 1;
                        vis[eid2] = 1;
                        extendableEdges.insert(eid1);//选一个方向的边加入数组即可
                    }
                }
            }
            if (extendableEdges.empty())//重新选择起始边
            {
                flag=true;
                break;
            }
        }
    }
    return motif;
}