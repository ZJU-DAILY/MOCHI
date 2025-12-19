
#include "AdvancedPeeling.h"

void AdvancedPeeling::modularityCount(int mod)
{
    if (n == 0)
        MM = -1;
    else if (mod == 1) {
        double num1 = (wc * 1.0) / n;
        double num2=(1.0)*wg;
        for(int i=0;i<motif->vertexType.size();i++)
        {
            int tp=motif->vertexType[i];
            num2=num2*(1.0*volC[tp]/volH[tp]);
        }
        num2=num2/n;
        MM = (1 - y) * num1 - y * num2;
    } else if (mod == 2) {
        double num1 = (wc * 1.0) / wg;
        double num2=(1.0)*wg;
        for(int i=0;i<motif->vertexType.size();i++)
        {
            int tp=motif->vertexType[i];
            num2=num2*(1.0*volC[tp]/volH[tp]);
        }
        num2=num2/wg;
        MM = (1 - y) * num1 - y * num2;
    } else if (mod == 3) {
        long double num1 = (wc * 1.0) / wg;
        long double num2 = 1.0;
        for(int i=0;i<motif->vertexType.size();i++)
        {
            num2=num2*(1.0*mmVolC[i]/mmVolH[i]);
        }
        MM = (1 - y) * num1 - y * num2;
    }
}

void AdvancedPeeling::checkConnectivity(map<int, set<int>> *conNodeGraph, set<int> *nodeList) {
    queue<int> queueList;
    vector<int> visNode(target->vertexNum + 1, -1);
    int q = *query.begin();
    queueList.push(q);
    visNode[q] = 1;
    int count = 0;
    while (!queueList.empty()) {
        int nowNode = queueList.front();
        queueList.pop();
        count++;
        for (int nxtNode: (*conNodeGraph)[nowNode]) {
            if (visNode[nxtNode] == -1) {
                visNode[nxtNode] = 1;
                queueList.push(nxtNode);
            }
        }
    }
    if (count != nowRsl.nodeGraph.size()) {
        for (int node: nowRsl.nodeGraph)
        {
            if (visNode[node] == -1)
                nodeList->insert(node);
        }
    }
}

bool
AdvancedPeeling::deleteNode(int ID, set<int> *peelingNode, set<int> *peelingMotif, map<pair<int, int>, int> *oldEdgeWeight,
                        set<pair<int, int>> *peelingConEdge) {
    set<int> influenceMotif = nowRsl.vertexToMotifID[ID];
    map<int, int> vertexInitialDegree;
    for (int imID: influenceMotif) {
        peelingMotif->insert(imID);
        for (const auto &entry: initialRsl->vertexToMotif[imID]) {
            int vID = entry.first;
            if (vID != ID && vertexInitialDegree.count(vID) == 0)
                vertexInitialDegree[vID] = nowRsl.vertexToMotifID[vID].size();
            nowRsl.vertexToMotifID[vID].erase(imID);
            if (nowRsl.vertexToMotifID[vID].empty()) {
                if (query.count(vID) > 0)
                {
                    if(optimize!=2&&optimize!=3)
                    {
                        dominateList[ID].insert(vID);
                    }
                    return false;
                }
                if(optimize!=4)
                {
                    if (influenceMotif.size() == vertexInitialDegree[vID] && vID != ID)
                    {
                        if(mod==3)
                        {
                            int nd=vID;
                            for(int idx:classifyMotif[target->vertexType[nd]])
                            {
                                mmVolC[idx]=mmVolC[idx]-pow(target->nbs[nd].size(),motif->nbs[idx].size());
                            }
                        }
                        nowRsl.nodeGraph.erase(vID);
                        peelingNode->insert(vID);
                        n--;
                        volC[target->vertexType[vID]] -= initialRsl->vertexToMotifID[vID].size();
                        nowRsl.vertexToMotifID.erase(vID);
                    } else if (influenceMotif.size() != vertexInitialDegree[vID] &&
                               vID !=
                               ID)
                    {
                        if(optimize!=2&&optimize!=3) {
                            dominateList[ID].insert(vID);
                        }
                        return false;
                    }
                }
                else if(vID !=
                        ID)
                {
                    if(optimize!=2&&optimize!=3) {
                        dominateList[ID].insert(vID);
                    }
                    return false;
                }
            }
        }
        wc--;
    }
    set<pair<int, int>> outEdge;
    for (const auto &entry: vertexDominate[ID])
    {
        int u = entry.first;
        int v = entry.second;
        if (nowRsl.conEdgeWeight.count(pair<int, int>(u, v)))
        {
            int num = nowRsl.conEdgeWeight[pair<int, int>(u, v)];
            (*oldEdgeWeight)[pair<int, int>(u, v)] = num;
            nowRsl.conEdgeWeight.erase(pair<int, int>(u, v));
            peelingConEdge->insert(pair<int, int>(u, v));
            nowRsl.conNodeGraph[u].erase(v);
            if (nowRsl.conNodeGraph[u].empty())
                nowRsl.conNodeGraph.erase(u);
        } else
            outEdge.insert(pair<int, int>(u, v));
    }
    for (const auto &entry: outEdge)
        vertexDominate[ID].erase(entry);
    set<int> oneHopNb;
    for (int nb: nowRsl.conNodeGraph[ID])
    {
        oneHopNb.insert(nb);
        int num = nowRsl.conEdgeWeight[pair<int, int>(ID, nb)];
        (*oldEdgeWeight)[pair<int, int>(ID, nb)] = num;
        (*oldEdgeWeight)[pair<int, int>(nb, ID)] = num;
        nowRsl.conEdgeWeight.erase(pair<int, int>(ID, nb));
        nowRsl.conEdgeWeight.erase(pair<int, int>(nb, ID));
        peelingConEdge->insert(pair<int, int>(ID, nb));
        peelingConEdge->insert(pair<int, int>(nb, ID));
        nowRsl.conNodeGraph[ID].erase(nb);
        if (nowRsl.conNodeGraph[ID].empty())
            nowRsl.conNodeGraph.erase(ID);
        nowRsl.conNodeGraph[nb].erase(ID);
        if (nowRsl.conNodeGraph[nb].empty())
            nowRsl.conNodeGraph.erase(nb);
    }
    if (oneHopEdge.count(ID)) {
        set<pair<int, int>> outHopEdge;
        for (const auto &entry: oneHopEdge[ID]) {
            int now = entry.first;
            int nxt = entry.second;
            if (nowRsl.conEdgeWeight.count(pair<int, int>(now, nxt))) {
                int num = nowRsl.conEdgeWeight[pair<int, int>(now, nxt)];
                (*oldEdgeWeight)[pair<int, int>(now, nxt)] = num;
                (*oldEdgeWeight)[pair<int, int>(nxt, now)] = num;
                num = num - influenceMotif.size();
                if (num <= 0) {
                    set<int> intersect;
                    set_intersection(nowRsl.vertexToMotifID[now].begin(), nowRsl.vertexToMotifID[now].end(),
                                     nowRsl.vertexToMotifID[nxt].begin(), nowRsl.vertexToMotifID[nxt].end(),
                                     inserter(intersect, intersect.begin()));
                    if (!intersect.empty()) {
                        num = intersect.size();
                        if (intersect.size() > (*oldEdgeWeight)[pair<int, int>(now, nxt)]) {
                            (*oldEdgeWeight)[pair<int, int>(now, nxt)] = num;
                            (*oldEdgeWeight)[pair<int, int>(nxt, now)] = num;
                        }
                        nowRsl.conEdgeWeight[pair<int, int>(now, nxt)] = num;
                        nowRsl.conEdgeWeight[pair<int, int>(nxt, now)] = num;
                    } else {
                        vertexDominate[ID].insert(pair<int, int>(now, nxt));
                        vertexDominate[ID].insert(pair<int, int>(nxt, now));
                        nowRsl.conEdgeWeight.erase(pair<int, int>(now, nxt));
                        nowRsl.conEdgeWeight.erase(pair<int, int>(nxt, now));
                        peelingConEdge->insert(pair<int, int>(now, nxt));
                        peelingConEdge->insert(pair<int, int>(nxt, now));
                        nowRsl.conNodeGraph[now].erase(nxt);
                        if (nowRsl.conNodeGraph[now].empty())
                            nowRsl.conNodeGraph.erase(now);
                        nowRsl.conNodeGraph[nxt].erase(now);
                        if (nowRsl.conNodeGraph[nxt].empty())
                            nowRsl.conNodeGraph.erase(nxt);
                    }
                } else {
                    nowRsl.conEdgeWeight[pair<int, int>(now, nxt)] = num;
                    nowRsl.conEdgeWeight[pair<int, int>(nxt, now)] = num;
                }
            } else
                outHopEdge.insert(pair<int, int>(now, nxt));
        }
        for (const auto &entry: outHopEdge)
            oneHopEdge[ID].erase(entry);
    } else {
        for (int now: oneHopNb) {
            set<int> nodeSet = nowRsl.conNodeGraph[now];
            for (int nxt: nodeSet) {
                if (oneHopNb.count(nxt) == 0)
                    continue;
                if (now < nxt)
                {
                    oneHopEdge[ID].insert(pair<int, int>(now, nxt));
                    int num = nowRsl.conEdgeWeight[pair<int, int>(now, nxt)];
                    (*oldEdgeWeight)[pair<int, int>(now, nxt)] = num;
                    (*oldEdgeWeight)[pair<int, int>(nxt, now)] = num;
                    num = num - influenceMotif.size();
                    if(optimize!=1&&optimize!=3)
                    {
                        if (num <= 0) {
                            set<int> intersect;
                            set_intersection(nowRsl.vertexToMotifID[now].begin(), nowRsl.vertexToMotifID[now].end(),
                                             nowRsl.vertexToMotifID[nxt].begin(), nowRsl.vertexToMotifID[nxt].end(),
                                             inserter(intersect, intersect.begin()));
                            if (!intersect.empty()) {
                                num = intersect.size();
                                if (intersect.size() > (*oldEdgeWeight)[pair<int, int>(now, nxt)]) {
                                    (*oldEdgeWeight)[pair<int, int>(now, nxt)] = num;
                                    (*oldEdgeWeight)[pair<int, int>(nxt, now)] = num;
                                }
                                nowRsl.conEdgeWeight[pair<int, int>(now, nxt)] = num;
                                nowRsl.conEdgeWeight[pair<int, int>(nxt, now)] = num;
                            } else {
                                vertexDominate[ID].insert(pair<int, int>(now, nxt));
                                vertexDominate[ID].insert(pair<int, int>(nxt, now));
                                nowRsl.conEdgeWeight.erase(pair<int, int>(now, nxt));
                                nowRsl.conEdgeWeight.erase(pair<int, int>(nxt, now));
                                peelingConEdge->insert(pair<int, int>(now, nxt));
                                peelingConEdge->insert(pair<int, int>(nxt, now));
                                nowRsl.conNodeGraph[now].erase(nxt);
                                if (nowRsl.conNodeGraph[now].empty())
                                    nowRsl.conNodeGraph.erase(now);
                                nowRsl.conNodeGraph[nxt].erase(now);
                                if (nowRsl.conNodeGraph[nxt].empty())
                                    nowRsl.conNodeGraph.erase(nxt);
                            }
                        } else {
                            nowRsl.conEdgeWeight[pair<int, int>(now, nxt)] = num;
                            nowRsl.conEdgeWeight[pair<int, int>(nxt, now)] = num;
                        }
                    }
                    else
                    {
                        set<int> intersect;
                        set_intersection(nowRsl.vertexToMotifID[now].begin(), nowRsl.vertexToMotifID[now].end(),
                                         nowRsl.vertexToMotifID[nxt].begin(), nowRsl.vertexToMotifID[nxt].end(),
                                         inserter(intersect, intersect.begin()));
                        if (!intersect.empty()) {
                            num = intersect.size();
                            if (intersect.size() > (*oldEdgeWeight)[pair<int, int>(now, nxt)]) {
                                (*oldEdgeWeight)[pair<int, int>(now, nxt)] = num;
                                (*oldEdgeWeight)[pair<int, int>(nxt, now)] = num;
                            }
                            nowRsl.conEdgeWeight[pair<int, int>(now, nxt)] = num;
                            nowRsl.conEdgeWeight[pair<int, int>(nxt, now)] = num;
                        } else {
                            vertexDominate[ID].insert(pair<int, int>(now, nxt));
                            vertexDominate[ID].insert(pair<int, int>(nxt, now));
                            nowRsl.conEdgeWeight.erase(pair<int, int>(now, nxt));
                            nowRsl.conEdgeWeight.erase(pair<int, int>(nxt, now));
                            peelingConEdge->insert(pair<int, int>(now, nxt));
                            peelingConEdge->insert(pair<int, int>(nxt, now));
                            nowRsl.conNodeGraph[now].erase(nxt);
                            if (nowRsl.conNodeGraph[now].empty())
                                nowRsl.conNodeGraph.erase(now);
                            nowRsl.conNodeGraph[nxt].erase(now);
                            if (nowRsl.conNodeGraph[nxt].empty())
                                nowRsl.conNodeGraph.erase(nxt);
                        }
                    }
                }
            }
        }
    }
    if(mod==3)
    {
        int nd=ID;
        for(int idx:classifyMotif[target->vertexType[nd]])
        {
            mmVolC[idx]=mmVolC[idx]-pow(target->nbs[nd].size(),motif->nbs[idx].size());
        }
    }
    nowRsl.nodeGraph.erase(ID);
    peelingNode->insert(ID);
    n--;
    volC[target->vertexType[ID]] -= initialRsl->vertexToMotifID[ID].size();
    nowRsl.vertexToMotifID.erase(ID);
    return true;
}

void AdvancedPeeling::returnBaseBack(set<int> *peelingNode, set<int> *peelingMotif, map<pair<int, int>, int> *oldEdgeWeight,
                                 set<pair<int, int>> *peelingConEdge) {
    volC=preVolC;
    n = pren;
    wc = prewc;
    if(mod==3)
        mmVolC=preMMVolC;
    for (const auto &entry: *oldEdgeWeight) {
        int weight = entry.second;
        nowRsl.conEdgeWeight[entry.first] = weight;
    }
    for (pair<int, int> entry: *peelingConEdge)
        nowRsl.conNodeGraph[entry.first].insert(entry.second);
    for (int node: *peelingNode)
    {
        nowRsl.nodeGraph.insert(node);
    }
    for (int motif: *peelingMotif) {
        for (const auto &entry: initialRsl->vertexToMotif[motif])
        {
            int node = entry.first;
            nowRsl.vertexToMotifID[node].insert(motif);
        }
    }
}

void AdvancedPeeling::computeContribution(int mod)
{
    contribution.clear();
    for (int node: nowRsl.nodeGraph) {
        double num1 = 0;
        double num2 = 0;
        if (mod == 1 || mod == 2) {
            int nown = nowRsl.nodeGraph.size() - 1;
            int nowwc = wc - nowRsl.vertexToMotifID[node].size();
            int nodeTp = target->vertexType[node];
            int volCNode = volC[nodeTp]-initialRsl->vertexToMotifID[node].size();
            if (mod == 1) {
                num1 = (nowwc * 1.0) / nown;
                num2=(1.0)*wg;
                for(int i=0;i<motif->vertexType.size();i++)
                {
                    int tp=motif->vertexType[i];
                    if(tp!=nodeTp)
                        num2=num2*(1.0*volC[tp]/volH[tp]);
                    else
                        num2=num2*(1.0*volCNode/volH[tp]);
                }
                num2=num2/nown;
            } else if (mod == 2) {
                num1 = (nowwc * 1.0) / wg;
                num2=(1.0)*wg;
                for(int i=0;i<motif->vertexType.size();i++)
                {
                    int tp=motif->vertexType[i];
                    if(tp!=nodeTp)
                        num2=num2*(1.0*volC[tp]/volH[tp]);
                    else
                        num2=num2*(1.0*volCNode/volH[tp]);
                }
                num2=num2/wg;
            }
        } else if (mod == 3) {
            int nowwc = wc - nowRsl.vertexToMotifID[node].size();
            num2=1.0;
            int k=target->nbs[node].size();
            for(int i=0;i<motif->vertexType.size();i++)
            {
                if(motif->vertexType[i]==target->vertexType[node])
                {
                    int delta=pow(k,motif->nbs[i].size());
                    num2=num2*(1.0*(mmVolC[i]-delta)/mmVolH[i]);
                }
                else
                {
                    num2=num2*(1.0*mmVolC[i]/mmVolH[i]);
                }
            }
            num1 = (nowwc * 1.0) / wg;
        }
        double ctr = (1 - y) * num1 - y * num2;
        contribution.push_back(gain(node, ctr));
    }
    std::sort(contribution.begin(), contribution.end());
}

void AdvancedPeeling::tarjan(int now, int root, int fa) {
    dfn[now] = low[now] = ++dfs_clock;
    int child = 0;
    for (int to: nowRsl.conNodeGraph[now]) {
        if (!dfn[to])
        {
            child++;
            tarjan(to, root, now);
            low[now] = min(low[now], low[to]);
            if (low[to] >= dfn[now] && now != root) cutNodes.insert(now);
        } else if (to != fa) low[now] = min(low[now], dfn[to]);
    }
    if (child >= 2 && now == root) cutNodes.insert(now);
}

void AdvancedPeeling::computeCutNode() {
    cutNodes.clear();
    dfn.clear();
    low.clear();
    dfs_clock = 0;
    tarjan(*query.begin(), *query.begin(), -1);
}

void AdvancedPeeling::greedyPeelingPlus() {
    bool peelingFlag = true;
    auto start_time = std::chrono::high_resolution_clock::now();
    config config;
    while (peelingFlag) {
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
        if (elapsed_time.count() > config.timeLimit) {
            result newRsl;
            rsl = newRsl;
            break;
        }
        peelingFlag = false;
        preVolC = volC;
        if(mod==3)
            preMMVolC=mmVolC;
        pren = n;
        prewc = wc;
        auto beginPLtime = std::chrono::high_resolution_clock::now();
        computeContribution(mod);
        for (gain entry: contribution) {
            int ID = entry.id;
            if (query.count(ID) > 0)
                continue;
            if(optimize!=2&&optimize!=3) {
                bool flag = false;
                if (!dominateList[ID].empty()) {
                    for (auto it = dominateList[ID].begin(); it != dominateList[ID].end();) {
                        int node = *it;
                        int num1 = nowRsl.vertexToMotifID[ID].size();
                        int num2 = nowRsl.vertexToMotifID[node].size();
                        if (nowRsl.nodeGraph.count(node) > 0 &&
                            (num1 != num2 || (num1 == num2 && nowRsl.conNodeGraph[ID].count(node) == 0)))
                        {
                            flag = true;
                            break;
                        } else
                            it = dominateList[ID].erase(it);
                    }
                    if (flag)
                        continue;
                }
            }
            bool judge;
            set<int> peelingNode;
            set<int> peelingMotif;
            map<pair<int, int>, int> oldEdgeWeight;
            set<pair<int, int>> peelingConEdge;
            judge = deleteNode(ID, &peelingNode, &peelingMotif, &oldEdgeWeight, &peelingConEdge);
            if (!judge) {
                returnBaseBack(&peelingNode, &peelingMotif, &oldEdgeWeight, &peelingConEdge);
                continue;
            }
            set<int> nodeList;
            checkConnectivity(&nowRsl.conNodeGraph, &nodeList);
            if (!nodeList.empty()) {
                if(optimize!=2&&optimize!=3)
                {
                    for (int node: nodeList)
                    {
                        dominateList[ID].insert(node);
                    }
                }
                returnBaseBack(&peelingNode, &peelingMotif, &oldEdgeWeight, &peelingConEdge);
                continue;
            }
            modularityCount(mod);
            nowRsl.modularity = MM;
            if(nowRsl.nodeGraph.size()%100==0)
                cout<<"left node num:" + to_string(nowRsl.nodeGraph.size())<<endl;
            if (nowRsl.modularity > rsl.modularity) {
                rsl.nodeGraph = nowRsl.nodeGraph;
                rsl.modularity = nowRsl.modularity;
                if(nowRsl.nodeGraph.size()%100==0)
                    cout << "optimal modularity:" + to_string(rsl.modularity) << endl;
            }
            peelingFlag = true;
            auto endPLtime = std::chrono::high_resolution_clock::now();
            chrono::duration<double> gap_time = endPLtime - beginPLtime;
            int nodeNum = nowRsl.nodeGraph.size();
            break;
        }
    }
}

void AdvancedPeeling::motifConnectedNodeGraph() {
    queue<int> queueList;
    vector<int> visNode(target->vertexNum + 1, -1);
    int q = *query.begin();
    queueList.push(q);
    nowRsl.nodeGraph.insert(q);
    visNode[q] = 1;
    while (!queueList.empty()) {
        int nowNode = queueList.front();
        queueList.pop();
        nowRsl.vertexToMotifID[nowNode] = initialRsl->vertexToMotifID[nowNode];
        nowRsl.motifSet.insert(initialRsl->vertexToMotifID[nowNode].begin(),
                               initialRsl->vertexToMotifID[nowNode].end());
        nowRsl.conNodeGraph[nowNode].insert(initialRsl->conNodeGraph[nowNode].begin(),
                                            initialRsl->conNodeGraph[nowNode].end());
        int tp=target->vertexType[nowNode];
        int volCNum=initialRsl->vertexToMotifID[nowNode].size();
        if(volC.count(tp)==0)
            volC[tp]=volCNum;
        else
            volC[tp]=volC[tp]+volCNum;
        for (int nxtNode: (initialRsl->conNodeGraph)[nowNode]) {
            int num = initialRsl->conEdgeWeight[pair<int, int>(nowNode, nxtNode)];
            nowRsl.conEdgeWeight[pair<int, int>(nowNode, nxtNode)] = num;
            if (visNode[nxtNode] == -1) {
                visNode[nxtNode] = 1;
                nowRsl.nodeGraph.insert(nxtNode);
                queueList.push(nxtNode);
            }
        }
    }
}

void AdvancedPeeling::initialMM() {
    for(int i=0;i<motif->vertexType.size();i++)
        classifyMotif[motif->vertexType[i]].insert(i);
    for(int i=0;i<target->vertexType.size();i++)
    {
        int tp=target->vertexType[i];
        if(classifyMotif.count(tp)==0)
            continue;
        int k=target->nbs[i].size();
        for(int u:classifyMotif[tp])
        {
            int num=pow(k,motif->nbs[u].size());
            if(mmVolH.count(u)==0)
                mmVolH[u]=num;
            else
                mmVolH[u]=mmVolH[u]+num;
            if(nowRsl.nodeGraph.count(i))
            {
                if(mmVolC.count(u)==0)
                    mmVolC[u]=num;
                else
                    mmVolC[u]=mmVolC[u]+num;
            }
        }
    }
}

result *AdvancedPeeling::search(set<int> inputQuery, result *matchResult, graph *inputMotif, graph *inputTarget,
                            double inputY, int inputMod,int inputOptimize) {
    optimize=inputOptimize;
    mod = inputMod;
    query = inputQuery;
    y = inputY;
    motif = inputMotif;
    target = inputTarget;
    initialRsl = matchResult;
    wg = initialRsl->vertexToMotif.size();
    m = motif->targetGraph.size();
    for(const auto&entry:initialRsl->vertexToMotifID)
    {
        int tp=target->vertexType[entry.first];
        int num=initialRsl->vertexToMotifID[entry.first].size();
        if(volH.count(tp)==0)
            volH[tp]=num;
        else
            volH[tp]=volH[tp]+num;
    }
    match match;
    match.generateConNodeGraph(initialRsl);
    motifConnectedNodeGraph();
    rsl.nodeGraph = nowRsl.nodeGraph;
    cout<<"connected node graph: "+ to_string(rsl.nodeGraph.size())<<endl;
    for (int q: query) {
        if (nowRsl.nodeGraph.count(q) == 0) {
            cout << "query nodes are not connected" << endl;
            result newRsl;
            rsl = newRsl;
            return &rsl;
        }
    }
    wc = nowRsl.motifSet.size();
    prewc = wc;
    n = nowRsl.nodeGraph.size();
    pren = n;
    preVolC=volC;
    if (mod == 3)
    {
        initialMM();
        preMMVolC=mmVolC;
    }
    modularityCount(mod);
    nowRsl.modularity = MM;
    rsl.nodeGraph = nowRsl.nodeGraph;
    rsl.modularity = nowRsl.modularity;
    greedyPeelingPlus();
    return &rsl;
}