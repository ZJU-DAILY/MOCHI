#include "localHIN.h"

void localHIN::MDMcount() {
    if (n == 0)
        MM = -1;
    else if(modularityMod == 0){
        double num1 = (wc * 1.0) / n;
        double num2 = (1.0) * wg;
        for (int i = 0; i < motif->vertexType.size(); i++) {
            int tp = motif->vertexType[i];
            num2 = num2 * (1.0 * volC[tp] / volH[tp]);
        }
        num2 = num2 / n;
        MM = (1 - y) * num1 - y * num2;
    }
    else if(modularityMod == 1)
    {
        double num1 = (wc * 1.0) / wg;
        double num2 = (1.0) * wg;
        for (int i = 0; i < motif->vertexType.size(); i++) {
            int tp = motif->vertexType[i];
            num2 = num2 * (1.0 * volC[tp] / volH[tp]);
        }
        num2 = num2 / wg;
        MM = (1 - y) * num1 - y * num2;
    }
}

void localHIN::computeMotifDistance()
{
    struct point {
        int nodeID;
        int dis;
    };
    for (int q: query) {
        queue<point> queueList;
        vector<int> visNode(target->vertexNum, -1);
        vertexDis[q] = 0;
        visNode[q] = 1;
        struct point nowpt;
        nowpt.nodeID = q;
        nowpt.dis = 0;
        queueList.push(nowpt);
        while (!queueList.empty()) {
            struct point pt = queueList.front();
            int nowNode = pt.nodeID;
            int dis = pt.dis;
            queueList.pop();
            for (int nxtNode: initialRsl->conNodeGraph[nowNode]) {
                if (visNode[nxtNode] == -1) {
                    struct point nxtpt;
                    nxtpt.nodeID = nxtNode;
                    nxtpt.dis = dis + 1;
                    queueList.push(nxtpt);
                    visNode[nxtNode] = 1;
                    if (vertexDis.count(nxtNode) > 0) {
                        if (dis + 1 < vertexDis[nxtNode] &&
                            query.count(nxtNode) == 0)
                            vertexDis[nxtNode] = dis + 1;
                    } else
                        vertexDis[nxtNode] = dis + 1;
                }
            }
        }
    }
}

void localHIN::sortByMotifDistance() {
    for (const auto &entry: vertexDis) {
        int vID = entry.first;
        int d = entry.second;
        sortedVertexDis[d].insert(vID);
    }
}

bool localHIN::deleteNodeWithDistance(int ID) {
    peelingOrder.push_back(ID);
    set<int> influenceMotif = nowRsl.vertexToMotifID[ID];
    set<int> influenceVertex;
    for (int imID: influenceMotif) {
        for (const auto &entry: initialRsl->vertexToMotif[imID]) {
            int vID = entry.first;
            nowRsl.vertexToMotifID[vID].erase(imID);
            if (nowRsl.vertexToMotifID[vID].empty()) {
                if (query.count(vID) > 0)
                    return false;
                if (vID != ID)
                {
                    nowRsl.nodeGraph.erase(vID);
                    peelingOrder.push_back(vID);
                    peelingList.erase(vID);
                    n--;
                    volC[target->vertexType[vID]] -= initialRsl->vertexToMotifID[vID].size();
                    nowRsl.vertexToMotifID.erase(vID);
                    peelingListGain.erase({nodeGain[vID],vID});
                    nodeGain.erase(vID);
                }
            } else if (vID != ID)
                influenceVertex.insert(vID);
        }
        wc--;
    }

    for (int node: influenceVertex) {
        if (nowRsl.nodeGraph.count(node) > 0) {
            peelingListGain.erase({nodeGain[node],node});
            double num = computeNodeGain(node);
            nodeGain[node] = num;
            if(peelingList.count(node))
                peelingListGain.insert({num,node});
        }
    }
    nowRsl.nodeGraph.erase(ID);
    peelingList.erase(ID);
    n--;
    volC[target->vertexType[ID]] -= initialRsl->vertexToMotifID[ID].size();
    nowRsl.vertexToMotifID.erase(ID);
    nodeGain.erase(ID);
    return true;
}

double localHIN::computeNodeGain(int node) {
    double num;
    if (nodeGainMod==1)
    {
        int tp=target->vertexType[node];
        int idx = x[tp];
        int volV = volH[tp];
        int dv = initialRsl->vertexToMotifID[node].size();
        int dvC = nowRsl.vertexToMotifID[node].size();
        double frac = (1.0) * dv / volV;
        num = (1.0) * pow(frac, idx) / dvC;
        nodeGain[node] = num;
    }
    else if (nodeGainMod==2)
    {
        int tp = target->vertexType[node];
        double idx = (1.0)*x[tp]/m;
        int dv = initialRsl->vertexToMotifID[node].size();
        int dvC = nowRsl.vertexToMotifID[node].size();
        double frac = (1.0) * dv / dvC;
        num = (1.0) * pow(frac, idx);
    }
    else if(nodeGainMod==3)
    {
        int dv = initialRsl->vertexToMotifID[node].size();
        int dvC = nowRsl.vertexToMotifID[node].size();
        num = (1.0) * dv / dvC;
    }
    return num;
}

bool localHIN::deleteLayerWithDistance() {
    set<int> influenceMotif;
    for (int node: peelingList) {
        nowRsl.nodeGraph.erase(node);
        n--;
        volC[target->vertexType[node]] -= initialRsl->vertexToMotifID[node].size();
        influenceMotif.insert(nowRsl.vertexToMotifID[node].begin(), nowRsl.vertexToMotifID[node].end());
    }
    for (int imID: influenceMotif) {
        for (const auto &entry: initialRsl->vertexToMotif[imID]) {
            int vID = entry.first;
            nowRsl.vertexToMotifID[vID].erase(imID);
            if (nowRsl.vertexToMotifID[vID].empty()) {
                if (query.count(vID) > 0)
                    return false;
                else
                    nowRsl.vertexToMotifID.erase(vID);
            }
        }
        wc--;
        nowRsl.motifSet.erase(imID);
    }
    peelingList.clear();
    return true;
}

void localHIN::distancePeeling() {
    map<int, set<int>, greater<int>>::iterator entry;
    if (mod == 0) {
        computeMotifDistance();
        sortByMotifDistance();
        entry = sortedVertexDis.begin();
    } else
        entry = sortedVertexDis.find(layer);
    auto start_time = chrono::high_resolution_clock::now();
    while (entry != sortedVertexDis.end()) {
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
        config config;
        if (elapsed_time.count() > config.timeLimit) {
            result newRsl;
            rsl = newRsl;
            break;
        }
        if (entry->first == 0)
            break;
        peelingList.clear();
        peelingList = entry->second;
        bool flag = true;
        //sort the vertices in the same layer
        peelingListGain.clear();
        for(int v:peelingList)
        {
            peelingListGain.insert({nodeGain[v],v});
        }
        while (!peelingList.empty()) {
            int maxNode = -1;
            //select the best vertex in the layer
            auto it = peelingListGain.begin();
            maxNode = it->second;
            peelingListGain.erase(it);
            flag = deleteNodeWithDistance(maxNode);
            if (!flag)
                break;
            MDMcount();
            nowRsl.modularity = MM;
            if(nowRsl.nodeGraph.size()%100==0)
                cout<<"left node num:" + to_string(nowRsl.nodeGraph.size())<<endl;
            if (nowRsl.modularity > rsl.modularity) {
                auto start_time3 = chrono::high_resolution_clock::now();
                rsl.nodeGraph = nowRsl.nodeGraph;
                rsl.modularity = nowRsl.modularity;
//                cout << "optimal modularity updated:" + to_string(rsl.modularity) << endl;
                auto end_time3 = chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed_time3 = end_time3 - start_time3;
//                cout << "update rsl time: " + to_string(elapsed_time3.count()) + "s" << endl;
            }
        }
        if (!flag)
            break;
        entry++;
    }
}

void localHIN::layerPeeling() {
    computeMotifDistance();
    sortByMotifDistance();
    auto start_time = chrono::high_resolution_clock::now();
    for (const auto &entry: sortedVertexDis) {
        if (layer == -1)
            layer = entry.first;
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
        config config;
        if (elapsed_time.count() > config.timeLimit) {
            result newRsl;
            rsl = newRsl;
            break;
        }
        if (entry.first == 0)
            break;
        peelingList = entry.second;
        bool flag = deleteLayerWithDistance();
        if (!flag)
            return;
        else {
            MDMcount();
            nowRsl.modularity = MM;
            if (nowRsl.modularity > rsl.modularity) {
                auto start_time3 = chrono::high_resolution_clock::now();
                prewc = wc;
                pren = n;
                preVolC = volC;
                layer = entry.first - 1;
                rsl.nodeGraph = nowRsl.nodeGraph;
                rsl.modularity = nowRsl.modularity;
                rsl.motifSet = nowRsl.motifSet;
                cout << "optimal modularity updated:" + to_string(rsl.modularity) << endl;
                auto end_time3 = chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed_time3 = end_time3 - start_time3;
                cout << "update rsl time: " + to_string(elapsed_time3.count()) + "s" << endl;
            }
        }
    }
}

void localHIN::returnNowRsl(result *newRsl)
{
    for (int mtf: newRsl->motifSet) {
        for (const auto &entry1: initialRsl->vertexToMotif[mtf])
            newRsl->vertexToMotifID[entry1.first].insert(mtf);
    }
}

void localHIN::motifConnectedNodeGraph() {
    queue<int> queueList;
    vector<int> visNode(target->vertexNum, -1);
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
        int tp = target->vertexType[nowNode];
        int num = initialRsl->vertexToMotifID[nowNode].size();
        if (volC.count(tp) == 0)
            volC[tp] = num;
        else
            volC[tp] = volC[tp] + num;
        for (int nxtNode: (initialRsl->conNodeGraph)[nowNode]) {
            if (visNode[nxtNode] == -1) {
                visNode[nxtNode] = 1;
                nowRsl.nodeGraph.insert(nxtNode);
                queueList.push(nxtNode);
            }
        }
    }
}

bool localHIN::generateNewQuery() {
    map<pair<int, int>, set<int>> path;
    map<pair<int, int>, int> edgeWeight;
    for (int q: query) {
        auto start_time = std::chrono::high_resolution_clock::now();
        cout << "q:" + to_string(q) << endl;
        vector<int> dist;
        vector<bool> visited;
        map<int, int> prev;
        struct nodePoint {
            int id;
            int dist;

            bool operator<(const nodePoint &other) const {
                return dist > other.dist;
            }
        };
        priority_queue<nodePoint> priorityQueue;
        for (int i = 0; i < target->vertexNum; i++) {
            dist.push_back(INF);
            visited.push_back(false);
        }
        dist[q] = 0;
        prev[q] = -1;
        struct nodePoint elm;
        elm.id = q;
        elm.dist = dist[q];
        priorityQueue.push(elm);
        set<int> visQuery = query;
        visQuery.erase(q);
        for (int node: query)
        {
            if (edgeWeight.count(pair<int, int>(q, node)))
                visQuery.erase(node);
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
        cout << "initial time: " + to_string(elapsed_time.count()) << endl;
        start_time = std::chrono::high_resolution_clock::now();
        while (!visQuery.empty() && !priorityQueue.empty()) {
            nodePoint nd = priorityQueue.top();
            if (nd.dist > dist[nd.id] || visited[nd.id])
            {
                priorityQueue.pop();
                continue;
            }
            int u = nd.id;
            priorityQueue.pop();
            visited[u] = true;
            if (visQuery.count(u))
            {
                visQuery.erase(u);
                edgeWeight[pair<int, int>(q, u)] = dist[u];
                edgeWeight[pair<int, int>(u, q)] = dist[u];
                if (visQuery.empty())
                    break;
            }
            for (int v: initialRsl->conNodeGraph[u]) {
//                int weight = maxWeight - initialRsl->conEdgeWeight[pair<int, int>(u, v)]+1;
                int weight = 1;
                if (!visited[v] && dist[u] != INF && dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    prev[v] = u;
                    priorityQueue.push(nodePoint{v, dist[v]});
                }
            }
        }
        end_time = std::chrono::high_resolution_clock::now();
        elapsed_time = end_time - start_time;
        cout << "compute shortest path time: " + to_string(elapsed_time.count()) << endl;
        start_time = std::chrono::high_resolution_clock::now();
        for (int node: query) {
            if (node == q)
                continue;
            if (path.count(pair<int, int>(q, node)))
                continue;
            if (edgeWeight.count(pair<int, int>(q, node)) == 0)
            {
                cout << "query nodes are disconnected: " + to_string(q) + "——" + to_string(node) << endl;
                return false;
            }
            int index = node;
            while (index != -1) {
                path[pair<int, int>(q, node)].insert(index);
                path[pair<int, int>(node, q)].insert(index);
                index = prev[index];
            }
        }
        end_time = std::chrono::high_resolution_clock::now();
        elapsed_time = end_time - start_time;
        cout << "store path time: " + to_string(elapsed_time.count()) << endl;
    }
    map<int, int> oldToNew;
    map<int, int> newToOld;
    int graph[query.size()][query.size()];
    int s = 0;
    for (int q: query) {
        oldToNew[q] = s;
        newToOld[s] = q;
        s++;
    }
    for (int i = 0; i < query.size(); i++) {
        for (int j = 0; j < query.size(); j++) {
            if (edgeWeight.count(pair<int, int>(newToOld[i], newToOld[j])))
                graph[i][j] = edgeWeight[pair<int, int>(newToOld[i], newToOld[j])];
            else
                graph[i][j] = 0;
        }
    }
    cout << "finish construct complete graph" << endl;
    int V = query.size();
    int parent[V];
    int key[V];
    int mstSet[V];
    for (int i = 0; i < V; i++)
        key[i] = INF, mstSet[i] = -1;
    key[0] = 0;
    parent[0] = -1;
    for (int count = 0; count < V - 1; count++) {
        int min = INF, u;
        for (int v = 0; v < V; v++)
            if (mstSet[v] == -1 && key[v] < min)
                min = key[v], u = v;
        mstSet[u] = 1;
        for (int v = 0; v < V; v++)
            if (graph[u][v] && mstSet[v] == -1 && graph[u][v] < key[v])
                parent[v] = u, key[v] = graph[u][v];
    }
    set<int> embeddings;
    set<int> newQuery = query;
    for (int i = 0; i < V; i++) {
        int u = i;
        int v = parent[i];
        if (v == -1)
            continue;
        newQuery.insert(path[pair<int, int>(newToOld[u], newToOld[v])].begin(),
                        path[pair<int, int>(newToOld[u], newToOld[v])].end());
    }
    for (int node: newQuery)
        embeddings.insert(initialRsl->vertexToMotifID[node].begin(), initialRsl->vertexToMotifID[node].end());
    set<int> cs = newQuery;
    int initialNum = newQuery.size();
    set<int> csReverse;
    while (!cs.empty())
    {
        int maxNum = -1;
        int maxCommonNum = -1;
        int maxID = -1;
        auto it = embeddings.begin();
        while (it != embeddings.end()) {
            if (cs.empty())
                break;
            int num = 0;
            int commonNum = 0;
            int motifID = *it;
            for (const auto &entry: initialRsl->vertexToMotif[motifID]) {
                if (cs.count(entry.first))
                    num++;
                if (csReverse.count(entry.first))
                    commonNum++;
            }
            if (num == motif->vertexNum)
            {
                for (const auto &entry: initialRsl->vertexToMotif[motifID]) {
                    cs.erase(entry.first);
                    csReverse.insert(entry.first);
                }
                embeddings.erase(it);
                maxNum = -1;
                maxCommonNum = -1;
                maxID = -1;
                it = embeddings.begin();
                continue;
            } else if (cs.size() == initialNum ||
                       (cs.size() != initialNum && commonNum >= 1))
            {
                if (num > maxNum ||
                    (num == maxNum && commonNum > maxCommonNum))
                {
                    maxNum = num;
                    maxID = motifID;
                    maxCommonNum = commonNum;
                }
            }
            ++it;
        }
        if (!cs.empty()) {
            for (const auto &entry: initialRsl->vertexToMotif[maxID]) {
                if (cs.count(entry.first)) {
                    cs.erase(entry.first);
                    csReverse.insert(entry.first);
                } else {
                    newQuery.insert(entry.first);
                    csReverse.insert(entry.first);
                }
            }
        }
    }
    query = newQuery;
    return true;
}

set<int> localHIN::search(set<int> inputQuery, result *matchResult, graph *inputMotif, graph *inputTarget,
                          double inputY, int inputMod,int gainMod,int inputModularityMod) {
    modularityMod = inputModularityMod;//0 denotes MDM, 1 denotes GMM
    nodeGainMod=gainMod;
    mod = inputMod;
    y = inputY;
    query = inputQuery;
    motif = inputMotif;
    target = inputTarget;
    initialRsl = matchResult;
    wg = initialRsl->vertexToMotif.size();
    m = motif->targetGraph.size();
    for (const auto &entry: initialRsl->vertexToMotifID) {
        int tp = target->vertexType[entry.first];
        int num = initialRsl->vertexToMotifID[entry.first].size();
        if (volH.count(tp) == 0)
            volH[tp] = num;
        else
            volH[tp] = volH[tp] + num;
    }
    for (int tp: motif->vertexType) {
        if (x.count(tp) == 0)
            x[tp] = 1;
        else
            x[tp]++;
    }
    match match;
    auto start_time = chrono::high_resolution_clock::now();
    if(!initialRsl->conNodeGraph.empty())
        cout<<"conNodeGraph not empty!"<<endl;
    maxWeight = match.generateConNodeGraph(initialRsl);
    motifConnectedNodeGraph();
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_time = end_time - start_time;
    cout << "find connected component time: " + to_string(elapsed_time.count()) + "s" << endl;
    bool flag = true;
    for (int q: query) {
        if (nowRsl.nodeGraph.count(q) == 0 || initialRsl->vertexToMotifID[q].size() == 0) {
            flag = false;
            break;
        }
    }
    if (!flag) {
        cout << "query nodes are disconnected or query nodes don't have embedding" << endl;
        result newRsl;
        rsl = newRsl;
        return rsl.nodeGraph;
    }
    if (query.size() > 1) {
        bool flag = generateNewQuery();
        if (!flag) {
            cout << "query nodes are disconnected" << endl;
            set<int> emptySet;
            return emptySet;
        }
    }
    else{
        int q = *query.begin();
        int mtf = *nowRsl.vertexToMotifID[q].begin();
        for(const auto&entry:initialRsl->vertexToMotif[mtf])
            query.insert(entry.first);
    }
    wc = nowRsl.motifSet.size();
    prewc = wc;
    n = nowRsl.nodeGraph.size();
    cout << "connected vertex num: " + to_string(n) << endl;
    cout << "connected motif num: " + to_string(wc) << endl;
    pren = n;
    preVolC = volC;
    MDMcount();
    nowRsl.modularity = MM;
    rsl.nodeGraph = nowRsl.nodeGraph;
    rsl.modularity = nowRsl.modularity;
    rsl.motifSet = nowRsl.motifSet;
    if (mod == 0) {
        for (int node: nowRsl.nodeGraph) {
            double num= computeNodeGain(node);
            nodeGain[node] = num;
        }
        distancePeeling();
    } else if (mod == 1) {
        auto start_time = chrono::high_resolution_clock::now();
        layerPeeling();
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
        cout << "layerPeelingTime: " + to_string(elapsed_time.count()) + "s" << endl;
        if (!rsl.motifSet.empty())
        {
            result newRsl;
            nowRsl = newRsl;
            nowRsl.nodeGraph = rsl.nodeGraph;
            nowRsl.modularity = rsl.modularity;
            nowRsl.motifSet = rsl.motifSet;
            start_time = chrono::high_resolution_clock::now();
            returnNowRsl(&nowRsl);
            end_time = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_time = end_time - start_time;
            cout << "ReConstructTime: " + to_string(elapsed_time.count()) + "s" << endl;
            cout << "left node: " + to_string(rsl.nodeGraph.size()) << endl;
            wc = prewc;
            n = pren;
            volC = preVolC;
            for (int node: nowRsl.nodeGraph) {
                double num = computeNodeGain(node);
                nodeGain[node] = num;
            }
            start_time = chrono::high_resolution_clock::now();
            distancePeeling();
            end_time = chrono::high_resolution_clock::now();
            elapsed_time = end_time - start_time;
            cout << "distancePeelingTime: " + to_string(elapsed_time.count()) + "s" << endl;
        }
    }
    return rsl.nodeGraph;
}
