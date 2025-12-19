# MOCHI: Motif-based Community Search over Large Heterogeneous Information Networks
-----------------------------------------------------------------------------------------------------------------
## Introduction
This repository contains the source code for **MOCHI: Motif-based Community Search over Large Heterogeneous Information Networks**. In this paper, we investigate the problem of motif-based community search over heterogeneous information networks (MOCHI). We devise a novel motif density modularity (MDM) to assess the motif cohesiveness within communities and propose three algorithms to solve the MOCHI problem.

## Algorithms
We implement all algorithms in C++. All experiments are performed on a Linux machine with 2.2 GHz CPU and 128 GB memory. The following files contain the source code of these algorithms:

1.`Naive_origin/NaiveOriginPeeling.cpp`: The basic algorithm for MOCHI.

2.`AdvancedHIN/AdvancedPeeling.cpp`: The MW-HIN-based algorithm for MOCHI.

3.`LocalHIN/localHIN.cpp`: The motif-distance-based algorithm for MOCHI.

We use **RapidMatch** [1] to find all the motif instances, which can be found [here](https://github.com/RapidsAtHKUST/RapidMatch).

*[1] Shixuan Sun, Xibo Sun, Yulin Che, Qiong Luo, and Bingsheng He. Rapidmatch: A holistic approach to subgraph query processing. Proceedings of the VLDB Endowment, 14(2):176–188, 2020.*

## Input Format
An example format of the input data for MOCHI is shown in the folder `data/example`.

1. `data/example`: An example of the HIN.

2. `data/example/motifs`: Examples of the motif.

3. `data/example/motifs/i_j/p/q.txt`: The q-th query with size p for j-th motif with size i.

4. `vertex.txt`: Each line starts with the vertex id, followed by the vertex type.

5. `edge.txt`: Each line starts with the edge id, followed by the edge type.

6. `graph.txt`: Each line starts with the vertex id, followed by a list of neighbor vertex id and edge id.

7. `queryGraph.txt/matchGraph.txt`: Input format of the graph in RapidMatch. It starts with `t |V| |E|`, where V is the vertex set and E is the edge set. The vertex and an edge are formatted as `v {vertex id} {vertex type} {vertex degree}` and `e {vertex id} {vertex id}`, respectively. The vertex id should start from 0 and the range should be limited within `[0,|V|-1]`.

## Usage
The running examples of three algorithms can be found at `BasicRunner.cpp`, `MWRunner.cpp`, and `MDRunner.cpp`. Please compile the project with the `CMakeLists.txt` and run `BasicRunner`, `MWRunner`, and `MDRunner` to execute the basic algorithm, MW algorithm, and MD algorithm, respectively. You can set the configuration in `Config/config.h`. The algorithm for generating motifs and query vertices can be found at `GenerateQuery.cpp`.

## Requirements
- GCC 11.2.1
- CMake 3.29.2
