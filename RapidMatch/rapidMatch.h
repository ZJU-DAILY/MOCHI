#ifndef MOTIFMODULARITY_RAPIDMATCH_H
#define MOTIFMODULARITY_RAPIDMATCH_H

#include "../Result/result.h"
#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include<map>
#include <computesetintersection.h>
#include "matching/matchingcommand.h"
#include "graph/graph.h"
#include "relation/catalog.h"
#include "matching/execution_tree_generator.h"
#include "pretty_print.h"
#include "global_variables.h"
#include "matching/preprocessor.h"
#include "matching/encoder.h"
#include "matching/query_plan_generator.h"
#include"../Config/config.h"
class rapidMatch {
public:
    Graph data_graph=new Graph(true);
    void loadGraph(string data_graph_file);
    bool globalMatchForGenerateQuery(string query_graph_file,result *matchRsl,map<int, set<int>> *inputSymmetricNodes);
    bool globalMatch(string query_graph_file,string data_graph_file, result *matchRsl,map<int, set<int>> *inputSymmetricNodes);
    bool globalMatchByMod(string query_graph_file,string data_graph_file,result *matchRsl,map<int, set<int>> *inputSymmetricNodes,int mod);
    void transformEmbeddingForMMC(uint32_t *matching_order, uint32_t order_length,
                                              uint32_t *embeddings, uint64_t embedding_count,result *matchRsl);
    bool globalMatchForMMC(string query_graph_file,string data_graph_file,result *matchRsl,map<int, set<int>> *inputSymmetricNodes);
    map<int, set<int>> *symmetricNodes;
    void transformEmbedding(uint32_t *matching_order, uint32_t order_length,
            uint32_t *embeddings, uint64_t embedding_count,result *matchRsl);
    double loadGraphTime=0.0;
};


#endif //MOTIFMODULARITY_RAPIDMATCH_H
