#include <iostream>
#include <chrono>
#include <malloc.h>
#include "Config/config.h"
#include "graph/graph.h"
#include "Util/Util.h"
#include "Match/match.h"
#include "RapidMatch/rapidMatch.h"
#include "AdvancedHIN/AdvancedPeeling.h"
namespace fs = std::filesystem;

//compute the statistic information from discovered community
static result constructGraphDetail(const set<int> &nodeGraph, result *matchResult, graph *target, graph *motif) {
    set<int> cutmotif;
    result graphRsl;
    for (int node: nodeGraph) {
        graphRsl.nodeGraph.insert(node);
        set<int> mtfs = matchResult->vertexToMotifID[node];
        for (int mtf: mtfs) {
            bool flag = true;
            for (const auto &entry: matchResult->vertexToMotif[mtf]) {
                int v = entry.first;
                if (nodeGraph.count(v) == 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                graphRsl.motifSet.insert(mtf);
                graphRsl.vertexToMotifID[node].insert(mtf);
            } else
                cutmotif.insert(mtf);
        }
    }

    //compute MW-HIN
    for (int motifID: graphRsl.motifSet) {
        map<int, int> motifs = matchResult->vertexToMotif[motifID];
        for (const auto &entry1: motifs) {
            for (const auto &entry2: motifs) {
                if (entry1.first != entry2.first) {
                    if (graphRsl.conEdgeWeight.count(pair<int, int>(entry1.first, entry2.first)) == 0) {
                        graphRsl.conNodeGraph[entry1.first].insert(entry2.first);
                        graphRsl.conEdgeWeight[pair<int, int>(entry1.first, entry2.first)] = 1;
                    } else
                        graphRsl.conEdgeWeight[pair<int, int>(entry1.first, entry2.first)]++;
                }
            }
        }
    }

    //compute MDM
    int m = motif->vertexNum;
    map<int, int> volC;
    map<int, int> volH;
    int n = nodeGraph.size();
    int wc = graphRsl.motifSet.size();
    int wg = matchResult->vertexToMotif.size();
    for (const auto &entry: matchResult->vertexToMotifID) {
        int tp = target->vertexType[entry.first];
        int num = matchResult->vertexToMotifID[entry.first].size();
        if (volH.count(tp) == 0)
            volH[tp] = num;
        else
            volH[tp] = volH[tp] + num;
        if (graphRsl.nodeGraph.count(entry.first)) {
            if (volC.count(tp) == 0)
                volC[tp] = num;
            else
                volC[tp] = volC[tp] + num;
        }
    }
    double num1 = (wc * 1.0) / n;
    double num2 = (1.0) * wg;
    for (int tp: motif->vertexType) {
        num2 = num2 * (1.0 * volC[tp] / volH[tp]);
    }
    num2 = num2 / n;
    double MM = num1 - num2;
    graphRsl.modularity = MM;

    //compute motif conductance
    int all_volc = 0;
    int all_volc_r = 0;
    for (const auto &entry: volC)
        all_volc += entry.second;
    for (const auto &entry: volH)
        all_volc_r += entry.second;
    all_volc_r = all_volc_r - all_volc;
    int cut = cutmotif.size();

    //compute conductance
    int edge_volc = 0;
    int edge_volc_r = 0;
    int edge_cut = 0;
    for (int node: nodeGraph) {
        int degree = target->nbs[node].size();
        edge_volc += degree;
        for (int nb: target->nbs[node]) {
            if (nodeGraph.count(nb) == 0)
                edge_cut++;
        }
    }
    edge_volc_r = target->edgeNum - edge_volc;
    graphRsl.communitySize = nodeGraph.size();
    return graphRsl;
}

int main() {
    config config;
    int motifNum = 20;
    int queryNum = 1;
    int motifSize = 4;
    int qsize = 1;
    for (const auto &dataFile: config.dataFile) {
        //STEP1: read graph
        string graphFile = config.path + dataFile + "/" + config.graphFile;
        string vertexFile = config.path + dataFile + "/" + config.vertexFile;
        string edgeFile = config.path + dataFile + "/" + config.edgeFile;
        auto start_time = chrono::high_resolution_clock::now();//记录起始时间
        graph buildGraph(graphFile, vertexFile, edgeFile);
        graph *target = &buildGraph;
        auto end_time = chrono::high_resolution_clock::now();//记录结束时间
        chrono::duration<double> elapsed_time = end_time - start_time;
        cout << "data reading time:" + to_string(elapsed_time.count()) + "s" << endl;
        int mid = 0;
        while (mid < motifNum) {
            mid++;
            //STEP2: read motif
            string motifRoot = config.path + dataFile + "/motifs/" + to_string(motifSize) + "_" +
                               to_string(mid);
            string motifGraphFile = motifRoot + "/" + config.graphFile;
            string motifVertexFile = motifRoot + "/" + config.vertexFile;
            string motifEdgeFile = motifRoot + "/" + config.edgeFile;
            graph buildMotifGraph(motifGraphFile, motifVertexFile, motifEdgeFile);
            graph *motif = &buildMotifGraph;

            //STEP3:match motif
            string matchGraphPath = config.path + dataFile + "/matchGraph.txt";
            string queryGraphPath = motifRoot + "/queryGraph.txt";
            start_time = chrono::high_resolution_clock::now();
            match globalMatch;
            map<int, set<int>> *symmetricNodes = globalMatch.symmetricDetect(motif);
            rapidMatch rapidMatch;
            result *matchResult = new result;
            rapidMatch.globalMatch(queryGraphPath, matchGraphPath, matchResult, symmetricNodes);
            cout << "end global match" << endl;
            end_time = chrono::high_resolution_clock::now();
            elapsed_time = end_time - start_time;
            double matchTime = elapsed_time.count() - rapidMatch.loadGraphTime;

            //STEP4:read query
            int qid = 0;
            while (qid < queryNum) {
                qid++;
                string queryPath = motifRoot + "/" + to_string(qsize) + "/" + to_string(qid) + ".txt";
                Util util;
                set<int> query = util.readQuery(queryPath);

                //STEP5: call algorithm to find community
                set<int> community;
                matchResult->motifGraph.clear();
                matchResult->conNodeGraph.clear();
                matchResult->conEdgeWeight.clear();
                double y = 0.5;
                start_time = chrono::high_resolution_clock::now();
                AdvancedPeeling advancedPeeling;
                community = advancedPeeling.search(query, matchResult, motif, target, y,
                                                   1, 0)->nodeGraph;
                end_time = chrono::high_resolution_clock::now();
                elapsed_time = end_time - start_time;
                double peelingTime = elapsed_time.count();

                if (community.empty())
                    continue;

                //STEP6: compute metric
                Util utils;
                result rsl = constructGraphDetail(community, matchResult, target, motif);
                rsl.matchTime = matchTime;
                rsl.executeTime = peelingTime;
                rsl.totalTime = matchTime + peelingTime;
                cout<<"running time: "<<rsl.totalTime<<endl;
                cout<<"MDM: "<<rsl.modularity<<endl;
                cout<<"community: "<<community<<endl;
            }
        }
    }
    return 0;
}
