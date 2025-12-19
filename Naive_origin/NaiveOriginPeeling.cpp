#include "NaiveOriginPeeling.h"

void NaiveOriginPeeling::MMcount() {
    if (n == 0)
        MM = -1;
    else {
        double num1 = (wc * 1.0) / n;
        double num2=(1.0)*wg;
        for(int i=0;i<motif->vertexType.size();i++)
        {
            int tp=motif->vertexType[i];
            num2=num2*(1.0*volC[tp]/volH[tp]);
        }
        num2=num2/n;
        MM = (1 - y) * num1 - y * num2;
    }
}
set<int> NaiveOriginPeeling::checkConnectivity() {
    set<int>nodeSet;
    queue<int> queueList;
    vector<int> visNode(target->vertexNum+1,-1);
    vector<int> visMotif(initialRsl->vertexToMotif.size()+1, -1);
    // find a motif instance of query, if not, return all the graph
    int motifID = -1;
    if (!query.empty()) {
        int firstQ = *query.begin();
        if (nowRsl.vertexToMotifID.count(firstQ) and !nowRsl.vertexToMotifID[firstQ].empty())
            motifID = *nowRsl.vertexToMotifID[firstQ].begin();
    }
    if(motifID==-1)
    {
        nodeSet = nowRsl.nodeGraph;
        return nodeSet;
    }
    // collect vertice with BFS
    int q = *query.begin();
    queueList.push(q);
    visNode[q] = 1;
    nodeSet.insert(q);
    while (!queueList.empty()) {
        int nowNode = queueList.front();
        queueList.pop();
        for (int mtf:nowRsl.vertexToMotifID[nowNode]) {
            if(visMotif[mtf]==1)
                continue;
            visMotif[mtf] = 1;
            for(const auto &entry:initialRsl->vertexToMotif[mtf])
            {
                if(visNode[entry.first]==-1)
                {
                    nodeSet.insert(entry.first);
                    queueList.push(entry.first);
                    visNode[entry.first]=1;
                }
            }
        }
    }
    // filter vertices
    set<int> result;
    set_difference(nowRsl.nodeGraph.begin(), nowRsl.nodeGraph.end(),
                        nodeSet.begin(), nodeSet.end(),
                        std::inserter(result, result.begin()));
    return result;
}

bool NaiveOriginPeeling::deleteNode(int ID, set<int> *peelingNode, set<int> *peelingMotif) {
    set<int> influenceMotif = nowRsl.vertexToMotifID[ID];
    set<int> influenceVertex;
    for (int imID: influenceMotif) {
        peelingMotif->insert(imID);
        for (const auto &entry: initialRsl->vertexToMotif[imID]) {
            int vID = entry.first;
            nowRsl.vertexToMotifID[vID].erase(imID);
            if (nowRsl.vertexToMotifID[vID].empty()&&vID!=ID) {
                if(query.count(vID))
                    return false;
                nowRsl.nodeGraph.erase(vID);
                peelingNode->insert(vID);
                n--;
                volC[target->vertexType[vID]] -= initialRsl->vertexToMotifID[vID].size();
                nowRsl.vertexToMotifID.erase(vID);
            }
        }
        wc--;
    }
    nowRsl.nodeGraph.erase(ID);
    peelingNode->insert(ID);
    n--;
    volC[target->vertexType[ID]] -= initialRsl->vertexToMotifID[ID].size();
    nowRsl.vertexToMotifID.erase(ID);
    return true;
}

void NaiveOriginPeeling::returnBaseBack(set<int> *peelingNode, set<int> *peelingMotif) {
    volC = preVolC;
    n = pren;
    wc = prewc;
    for (int node: *peelingNode)
        nowRsl.nodeGraph.insert(node);
    for (int motif: *peelingMotif) {
        for (const auto &entry: initialRsl->vertexToMotif[motif])
        {
            int node = entry.first;
            nowRsl.vertexToMotifID[node].insert(motif);
        }
    }
}

void NaiveOriginPeeling::greedyPeelingPlus() {
    bool peelingFlag = true;
    auto start_time = std::chrono::high_resolution_clock::now();
    config config;
    double sumTime=0.0;
    int nodeCount=0;
    int nodeNum = nowRsl.nodeGraph.size();
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
        pren = n;
        prewc = wc;
        auto beginPLtime=std::chrono::high_resolution_clock::now();

        // compute gain of each vertex
        set<int>nodeSet = nowRsl.nodeGraph;
        int maxGainID = -1;
        double maxGain = -1;
        int nowID = 0;
        cout<<"round: "<<to_string(nodeCount+1)<<endl;
        for(int nd:nodeSet)
        {
            auto start=std::chrono::high_resolution_clock::now();
//            cout<<to_string(nodeCount+1)+"."+ to_string(nowID+1)<<endl;
            nowID++;
            if(query.count(nd))
                continue;
            set<int> peelingNode;
            set<int> peelingMotif;

            // remove vertex
            bool flag = deleteNode(nd, &peelingNode, &peelingMotif);
            if(flag)
            {
                // check connectivity and return the disconnected vertex
                set<int>vset = checkConnectivity();
                bool valid = true;
                for(int q:query)
                {
                    if(vset.count(q))
                    {
                        valid = false;
                        break;
                    }
                }
                if(valid)
                {
                    for(int node:vset)
                        deleteNode(node, &peelingNode, &peelingMotif);
                    MMcount();
                    if (MM>maxGain)
                    {
                        maxGain = MM;
                        maxGainID = nd;
                    }
                }
            }
            returnBaseBack(&peelingNode, &peelingMotif);
            auto end = std::chrono::high_resolution_clock::now();
            chrono::duration<double> gap_time = end - start;
        }

        // remove best vertex
        if(maxGainID!=-1)
        {
            set<int> peelingNode;
            set<int> peelingMotif;
            deleteNode(maxGainID, &peelingNode, &peelingMotif);
            set<int>vset = checkConnectivity();
            for(int node:vset)
                deleteNode(node, &peelingNode, &peelingMotif);
            MMcount();
            nowRsl.modularity = MM;
            cout << "modularity:" + to_string(nowRsl.modularity)<<endl;
            cout << "left nodeGraph size:" + to_string(nowRsl.nodeGraph.size()) << endl;
            if (nowRsl.modularity > rsl.modularity) {
                rsl.nodeGraph = nowRsl.nodeGraph;
                rsl.modularity = nowRsl.modularity;
                cout << "optimal modularity:" + to_string(rsl.modularity) << endl;
            }
            peelingFlag=true;
            auto endPLtime = std::chrono::high_resolution_clock::now();
            chrono::duration<double> gap_time = endPLtime - beginPLtime;
            sumTime=sumTime+gap_time.count();
            if(sumTime>config.timeLimit)
            {
                cout<<"out of time"<<endl;
                result newRsl;
                rsl = newRsl;
                return;
            }
            nodeCount++;
        }
    }
}
void NaiveOriginPeeling::findConnectComp()
{
    queue<int> queueList;
    vector<int> visMotif(initialRsl->vertexToMotif.size()+1, -1);
    vector<int>visNode(target->vertexNum+1,-1);
    int q = *query.begin();
    queueList.push(q);
    nowRsl.nodeGraph.insert(q);
    while(!queueList.empty())
    {   int nowNode = queueList.front();
        queueList.pop();
        for(int mtf:initialRsl->vertexToMotifID[nowNode])
        {
            if(visMotif[mtf]==1)
                continue;
            visMotif[mtf]=1;
            nowRsl.motifSet.insert(mtf);
            for(const auto&entry:initialRsl->vertexToMotif[mtf])
            {
                int nodeID=entry.first;
                if(visNode[nodeID]==-1)
                {
                    nowRsl.nodeGraph.insert(nodeID);
                    visNode[nodeID]=1;
                    queueList.push(nodeID);
                }
                nowRsl.vertexToMotifID[nodeID].insert(mtf);
            }
        }
    }
}

result*
NaiveOriginPeeling::search(set<int> inputQuery, result *matchResult, graph *inputMotif, graph *inputTarget, double inputY) {
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
    findConnectComp();
    for(int node:nowRsl.nodeGraph)
    {
        int tp=target->vertexType[node];
        int num=initialRsl->vertexToMotifID[node].size();
        if(volC.count(tp)==0)
            volC[tp]=num;
        else
            volC[tp]=volC[tp]+num;
    }
    wc = nowRsl.motifSet.size();
    prewc = wc;
    n = nowRsl.nodeGraph.size();
    pren = n;
    preVolC = volC;
    MMcount();
    nowRsl.modularity = MM;
    rsl.nodeGraph = nowRsl.nodeGraph;
    rsl.modularity = nowRsl.modularity;
    greedyPeelingPlus();
    return &rsl;
}
