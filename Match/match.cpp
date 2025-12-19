
#include "match.h"
#include "../Config/config.h"

set<int> match::
getR2(int nxtu, map<int, int> M1, map<int, int> M2, map<int,int> mcm)
{
    set<int> nowintersectVec;
    bool flag = false;
    for (const auto &entry: M1) {
        int uc = entry.first;
        set<int> nxtintersectVec;
        if (motif->edges.count(
                pair<int, int>(nxtu, uc))) {
            flag = true;
            int vc = M1[uc];
            for (int nxtv:target->nbs[vc]) {
                if (!isGlobal) {
                    if (enumerated[nxtv] == 1)
                        continue;
                }
                if (target->vertexType[nxtv] == motif->vertexType[nxtu] && M2.count(nxtv) == 0)
                {
                    if(mcm.count(nxtu))
                    {
                        if(nxtv>mcm[nxtu])
                            nxtintersectVec.insert(nxtv);
                    }
                    else
                        nxtintersectVec.insert(nxtv);
                }
            }
            if (nxtintersectVec.empty())
            {
                nowintersectVec.clear();
                return nowintersectVec;
            }
            if (nowintersectVec.empty())
                nowintersectVec = nxtintersectVec;
            else {
                set<int>intersect;
                set_intersection(nowintersectVec.begin(),nowintersectVec.end(),nxtintersectVec.begin(),nxtintersectVec.end(),
                                 inserter(intersect,intersect.begin()));
                if (intersect.empty())
                    return intersect;
                else
                    nowintersectVec=intersect;
            }
        }
    }
    if (!flag) {
        for (int i = 0; i < target->vertexType.size(); i++) {
            if (target->vertexType[i] == motif->vertexType[nxtu]&&enumerated[i]==-1)
                nowintersectVec.insert(i);
        }
    }
    return nowintersectVec;
}

void match::matchByDFS(map<int, int> M1, map<int, int> M2, map<int,int> mcm,set<int>extendableNodes) {
    if(currentRsl.outLimitFlag)
        return;
    config config;
    if(match_gap.count()>config.timeLimit||currentRsl.vertexToMotif.size()>config.motifNumLimit)
    {
        currentRsl.outLimitFlag=true;
        return;
    }
    else
    {
        match_end=chrono::high_resolution_clock::now();
        match_gap=match_end-match_start;
    }
    if (M1.size() == motif->targetGraph.size()) {
        currentRsl.vertexToMotif.push_back(M2);
        int mtfID = currentRsl.vertexToMotif.size() - 1;
        cout << "find motif: " + to_string(currentRsl.vertexToMotif.size()) << endl;
        for (const auto &entry: M1) {
            int vID = entry.second;
            cout<<to_string(vID)+"-";
            currentRsl.vertexToMotifID[vID].insert(mtfID);
            if (currentRsl.nodeGraphVis[entry.second]== -1)
            {
                currentRsl.nodeGraphVis[entry.second]=1;
                currentRsl.nodeGraph.insert(entry.second);
                newNodes.insert(entry.second);
            }
        }
        cout<<endl;
        return;
    }
    int nxtu=*extendableNodes.begin();
    set<int>nxtExtendableNodes=extendableNodes;
    nxtExtendableNodes.erase(nxtu);
    for(int node:motif->nbs[nxtu])
    {
        if(M1.count(node)==0)
            nxtExtendableNodes.insert(node);
    }
    set<int> R2 = getR2(nxtu, M1, M2, mcm);
    for (int nxtv: R2) {
        map<int, int> nxtM1 = M1;
        map<int, int> nxtM2 = M2;
        map<int,int> nxtMcm = mcm;
        nxtM1[nxtu] = nxtv;
        nxtM2[nxtv] = nxtu;
        if(symmetricNodes.count(nxtu))
        {
            for(int node:symmetricNodes[nxtu])
                nxtMcm[node]=nxtv;
            nxtMcm[nxtu]=nxtv;
        }
        matchByDFS(nxtM1, nxtM2, nxtMcm,nxtExtendableNodes);
    }
}
void match::generateMotifGraph(result *rsl)
{
    int maxDegree=-1;
    int sum=0;
    int count=0;
    for(int motifID=0;motifID<rsl->vertexToMotif.size();motifID++)
    {
        const auto &motifs=rsl->vertexToMotif[motifID];
        auto start_time = std::chrono::high_resolution_clock::now();
        count++;
        vector<int>visMotif(rsl->vertexToMotif.size()+1,-1);
        visMotif[motifID]=1;
        for(const auto&entry:motifs)
        {
            for(int mtf:rsl->vertexToMotifID[entry.first])
            {
                if(visMotif[mtf]==-1)
                {
                    rsl->motifGraph[motifID].insert(mtf);
                    visMotif[mtf]=1;
                }
            }
        }
        int num=rsl->motifGraph[motifID].size();
        if(num>maxDegree)
            maxDegree=rsl->motifGraph[motifID].size();
        sum+=rsl->motifGraph[motifID].size();
        if(count%100==0)
        {
            cout<<"motif degree: "+ to_string(num)<<endl;
            auto end_time = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_time = end_time - start_time;
            cout<<"finish motif graph: "+ to_string(motifID)+" "+to_string(elapsed_time.count()*100)+"s"<<endl;
            int num=rsl->vertexToMotif.size();
            config config;
            if(elapsed_time.count()*num>config.timeLimit/20)
            {
                cout<<"out of time"<<endl;
                rsl->motifGraph.clear();
                break;
            }
        }
    }
    double avgDegree=(1.0)*sum/(2*rsl->vertexToMotif.size());
    cout<<"maxDegree: "+ to_string(maxDegree)+" avgDegree: "+ to_string(avgDegree)<<endl;
}
int match::generateConNodeGraph(result *rsl) {
    int edgeNum = 0;
    int maxWeight = -1;
    for (const auto &motifs: rsl->vertexToMotif) {
        auto start_time = std::chrono::high_resolution_clock::now();
        for (const auto &entry1: motifs) {
            for (const auto &entry2: motifs) {
                if (entry1.first != entry2.first) {
                    if (rsl->conEdgeWeight.count(pair<int, int>(entry1.first, entry2.first)) == 0) {
                        rsl->conNodeGraph[entry1.first].insert(entry2.first);
                        rsl->conEdgeWeight[pair<int, int>(entry1.first, entry2.first)] = 1;
                        if (maxWeight < 1)
                            maxWeight = 1;
                        edgeNum++;
                    } else {
                        rsl->conEdgeWeight[pair<int, int>(entry1.first, entry2.first)]++;
                        int num = rsl->conEdgeWeight[pair<int, int>(entry1.first, entry2.first)];
                        if (num > maxWeight)
                            maxWeight = num;
                    }
                }
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed_time = end_time - start_time;
    }
    int vertexNum = rsl->nodeGraph.size();
    edgeNum = edgeNum / 2;
    if(vertexNum == 0)
    {
        cout << "no embedding" <<endl;
        return -1;
    }
    double density = edgeNum / vertexNum;
    cout << "vertexNum: " + to_string(vertexNum) << endl;
    cout << "edgeNum: " + to_string(edgeNum) << endl;
    cout << "density: " + to_string(density) << endl;
    return maxWeight;
}
result* match::globalMatch(graph *inputMotif, graph *inputTarget) {
    match_start=std::chrono::high_resolution_clock::now();
    match_end=std::chrono::high_resolution_clock::now();
    match_gap=match_end-match_start;
    auto start_time = std::chrono::high_resolution_clock::now();
    config config;
    motif = inputMotif;
    target = inputTarget;
    isGlobal = true;
    for(int i=0;i<target->vertexNum;i++)
        currentRsl.nodeGraphVis.push_back(-1);
    symmetricDetect(motif);
    int uc = 0;
    int typeuc = motif->vertexType[uc];
    int tag[motif->vertexType.size()];
    memset(tag,-1,sizeof(tag));
    for (int i = 0; i < target->vertexType.size(); i++) {
        if (target->vertexType[i] == typeuc)
        {
            int vc = i;
            map<int, int> M1;
            map<int, int> M2;
            map<int,int> mcm;
            for(const auto&entry:symmetricNodes)
                mcm[entry.first]=-1;
            M1[uc] = vc;
            M2[vc] = uc;
            if(symmetricNodes.count(uc))
            {
                for(int node:symmetricNodes[uc])
                    mcm[node]=vc;
                mcm[uc]=vc;
            }
            set<int>extendableNodes;
            for(int node:motif->nbs[uc])
                extendableNodes.insert(node);
            matchByDFS(M1, M2, mcm,extendableNodes);
            if(currentRsl.outLimitFlag)
            {
                result newRsl;
                currentRsl=newRsl;
                currentRsl.outLimitFlag=true;
                return &currentRsl;
            }
        }
    }
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_time = end_time - start_time;
    cout<<"find match graph time: "+to_string(elapsed_time.count())+"s"<<endl;
    return &currentRsl;
}

result* match::localMatch(set<int> query, graph *inputMotif, graph *inputTarget) {
    match_start=std::chrono::high_resolution_clock::now();
    match_end=std::chrono::high_resolution_clock::now();
    match_gap=match_end-match_start;
    config config;
    motif = inputMotif;
    target = inputTarget;
    isGlobal = false;
    for(int i=0;i<target->vertexNum;i++)
        currentRsl.nodeGraphVis.push_back(-1);
    symmetricDetect(motif);
    vector<int> emrt(target->vertexNum, -1);
    for(int i=0;i<target->vertexNum;i++)
        enumerated.push_back(-1);
    queue<int> qlist;
    vector<int> vis(target->vertexNum, -1);
    int qbegin = *query.begin();
    qlist.push(qbegin);
    vis[qbegin] = 1;
    while (!qlist.empty()) {
        cout << "qlist size:" + to_string(qlist.size()) + "--------------------" << endl;
        int q = qlist.front();
        qlist.pop();
        int typeq = target->vertexType[q];
        int uq = -1;
        newNodes.clear();
        int tag[motif->vertexType.size()];
        memset(tag,-1,sizeof(tag));
        for (int i = 0; i < motif->vertexType.size(); i++) {
            if (motif->vertexType[i] == typeq&&tag[i]==-1) {
                uq = i;
                if(symmetricNodes.count(i))
                {
                    for(int node:symmetricNodes[i])
                        tag[node]=1;
                    tag[i]=1;
                }
                map<int, int> M1;
                map<int, int> M2;
                map<int,int> mcm;
                for(const auto&entry:symmetricNodes)
                    mcm[entry.first]=-1;
                M1[uq] = q;
                M2[q] = uq;
                set<int>extendableNodes;
                for(int node:motif->nbs[uq])
                    extendableNodes.insert(node);
                matchByDFS(M1, M2, mcm,extendableNodes);
                if(currentRsl.outLimitFlag)
                {
                    cout <<"out limtit"<<endl;
                    result newRsl;
                    currentRsl=newRsl;
                    currentRsl.outLimitFlag=true;
                    return &currentRsl;
                }
            }
        }
        if (uq == -1 && query.count(uq)) {
            cout << "motif does not have type of query node,the query node is: " + to_string(q) << endl;
            result newRsl;
            currentRsl = newRsl;
            return &currentRsl;
        }
        enumerated[q] = 1;
        for (int node: newNodes) {
            if (vis[node] == -1)
            {
                qlist.push(node);
                vis[node] = 1;
            }
        }
    }

    for (int q: query) {
        if (currentRsl.nodeGraph.count(q) == 0) {
            cout << "query node is not motif-connected:" + to_string(q) << endl;
            result newRsl;
            currentRsl = newRsl;
            return &currentRsl;
        }
    }
    return &currentRsl;
}

map<int, set<int>>* match::symmetricDetect(graph *inputMotif) {
    motif=inputMotif;
    map<int, vector<int>> classNodes;
    for (int i = 0; i < motif->vertexType.size(); i++) {
        int tp = motif->vertexType[i];
        classNodes[tp].push_back(i);
    }
    for (const auto &entry: classNodes) {
        if (entry.second.size() > 1) {
            vector<int> nodes = entry.second;
            for (int i = 0; i < nodes.size(); i++) {
                for (int j = i + 1; j < nodes.size(); j++) {
                    int u = nodes[i];
                    int v = nodes[j];
                    queryGraph = *motif;
                    dataGraph = *motif;
                    queryGraph.vertexType[u] = -1;
                    queryGraph.vertexType[v] = -2;
                    dataGraph.vertexType[u] = -2;
                    dataGraph.vertexType[v] = -1;
                    symmetricFlag = false;
                    map<int, int> M1;
                    map<int, int> M2;
                    int uc=0;
                    for(int vc=0;vc<dataGraph.vertexType.size();vc++)
                    {
                        if(dataGraph.vertexType[vc]!=queryGraph.vertexType[uc])
                            continue;
                        M1[uc]=vc;
                        M2[vc]=uc;
                        set<int>extendableNodes;
                        for(int node:queryGraph.nbs[uc])
                            extendableNodes.insert(node);
                        isIsomorphism(M1, M2,extendableNodes);
                        if (symmetricFlag) {
                            symmetricNodes[u].insert(v);
                            symmetricNodes[v].insert(u);
                            break;
                        }
                    }
                }
            }
        }
    }
    return &symmetricNodes;
}

set<int> match::
getCandidate(int nxtu, map<int, int> M1, map<int, int> M2) {
    set<int> nowintersectVec;
    bool flag = false;
    for (const auto &entry: M1) {
        int uc = entry.first;
        set<int> nxtintersectVec;
        if (queryGraph.edges.count(
                pair<int, int>(nxtu, uc))) {
            flag = true;
            int vc = M1[uc];
            for (int nxtv:dataGraph.nbs[vc]) {
                if (dataGraph.vertexType[nxtv] == queryGraph.vertexType[nxtu] && M2.count(nxtv) == 0)
                {
                    nxtintersectVec.insert(nxtv);
                }
            }
            if (nxtintersectVec.size() == 0)
            {
                nowintersectVec.clear();
                return nowintersectVec;
            }
            if (nowintersectVec.size() == 0)
                nowintersectVec = nxtintersectVec;
            else {
                set<int>intersect;
                set_intersection(nowintersectVec.begin(),nowintersectVec.end(),nxtintersectVec.begin(),nxtintersectVec.end(),
                                 inserter(intersect,intersect.begin()));
                if (intersect.size() == 0)
                    return intersect;
                else
                    nowintersectVec=intersect;
            }
        }
    }
    if (!flag) {
        for (int i = 0; i < dataGraph.vertexType.size(); i++) {
            if (dataGraph.vertexType[i] == queryGraph.vertexType[nxtu])
                nowintersectVec.insert(i);
        }
    }
    return nowintersectVec;
}

void match::isIsomorphism(map<int, int> M1, map<int, int> M2,set<int>extendableNodes) {
    if (symmetricFlag)
        return;
    if (M1.size() == queryGraph.targetGraph.size()) {
        symmetricFlag = true;
        return;
    }
    int nxtu=*extendableNodes.begin();
    set<int>nxtExtendableNodes=extendableNodes;
    nxtExtendableNodes.erase(nxtu);
    for(int node:motif->nbs[nxtu])
    {
        if(M1.count(node)==0)
            nxtExtendableNodes.insert(node);
    }
    set<int> R2 = getCandidate(nxtu, M1, M2);
    for (int nxtv: R2) {
        map<int, int> nxtM1 = M1;
        map<int, int> nxtM2 = M2;
        nxtM1[nxtu] = nxtv;
        nxtM2[nxtv] = nxtu;
        isIsomorphism(nxtM1, nxtM2,nxtExtendableNodes);
    }
}