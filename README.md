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

## Datasets
We use five real-world HINs in our experiments, including [CIDeRplus](https://mips.helmholtz-muenchen.de/CIDeRplus), [DBLP](https://github.com/Jhy1993/HAN/tree/master), [TMDB](https://www.kaggle.com/tmdb/tmdb-movie-metadata), [Freebase](https://freebase-easy.cs.uni-freiburg.de/dump/), [DBpedia](https://www.dbpedia.org/).
The processed data can be found [here](https://drive.google.com/file/d/1i_HRNss7Sw-OxlY3mtYARn53eJ1cni3c/view?usp=sharing).
The statistics and references of these datasets can be found in our paper.

## Input Format
An example format of the input data for MOCHI is shown in the folder `data/example`.

1. `data/example`: An example of the HIN.

2. `data/example/motifs`: Examples of the motif.

3. `data/example/motifs/i_j/p/q.txt`: The q-th query with size p for j-th motif with size i.

4. `vertex.txt`: Each line starts with the vertex id, followed by the vertex type.

5. `edge.txt`: Each line starts with the edge id, followed by the edge type.

6. `graph.txt`: Each line starts with the vertex id, followed by a list of neighbor vertex id and edge id.

7. `queryGraph.txt/matchGraph.txt`: Input format of the graph in RapidMatch. It starts with `t |V| |E|`, where V is the vertex set and E is the edge set. A vertex and an edge are formatted as `v {vertex id} {vertex type} {vertex degree}` and `e {vertex id} {vertex id}` respectively. The vertex id should start from 0 and the range should be limited within `[0,|V|-1]`.

## Usage
The running examples of three algorithms can be found at `BasicRunner.cpp`, `MWRunner.cpp`, and `MDRunner.cpp`. Please compile the project with the `CMakeLists.txt` and run `BasicRunner`, `MWRunner`, and `MDRunner`, respectively, to execute the basic algorithm, MW algorithm, and MD algorithm. You can set the dataset and path in `Config/config.h`.
The algorithm for generating motifs and query vertices can be found at `GenerateQuery.cpp`. The generated queries will be stored at `data/{dataset}/RandomMotif`. Please move them to the folder `data/{dataset}/motifs` when executing queries.

## Requirements
- GCC 11.2.1
- CMake 3.29.2
