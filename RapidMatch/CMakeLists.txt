set(RAPIDMATCH_SRC
        rapidMatch.h
        rapidMatch.cpp
        matching/matchingcommand.cpp
        matching/matchingcommand.h
        matching/execution_tree.h
        matching/execution_tree_node.h
        matching/execution_tree_node.cpp
        matching/execution_tree_generator.cpp
        matching/execution_tree_generator.h
        matching/execution_tree.cpp
        matching/preprocessor.cpp
        matching/preprocessor.h
        matching/encoder.cpp
        matching/encoder.h
        matching/query_plan_generator.cpp
        matching/query_plan_generator.h
        graph/graph.cpp
        utility/graphoperations.cpp
        utility/commandparser.cpp
        utility/commandparser.h
        utility/computesetintersection.h
        utility/computesetintersection.cpp
        utility/QFilter.h
        utility/global_variables.h
        utility/global_variables.cpp
        utility/han/intersection_algos.cpp
        utility/han/utils/util.cpp
        utility/han/intersection_algos.hpp
        utility/hashjoin/hash_helpers.h
        utility/hashjoin/nop_join.h
        utility/relation/flat_relation.h
        utility/hashjoin/nop_join.cpp
        utility/leapfrogtriejoin/leapfrogtriejoin.cpp
        utility/leapfrogtriejoin/leapfrogtriejoin.h
        utility/relation/trie_relation.h
        utility/primitive/search.cpp
        utility/primitive/search.h
        utility/relation/encoded_trie_relation.h
        utility/relation/hash_relation.h
        utility/relation/catalog.h
        utility/relation/catalog.cpp
        utility/relation/edge_relation.h
        utility/primitive/scan.cpp
        utility/primitive/scan.h
        utility/primitive/semi_join.cpp
        utility/primitive/semi_join.h
        utility/primitive/projection.cpp
        utility/primitive/projection.h
        utility/primitive/nlf_filter.cpp
        utility/primitive/nlf_filter.h
        utility/nucleus_decomposition/nd.h
        utility/nucleus_decomposition/bucket.h
        utility/nucleus_decomposition/bucket.cpp
        utility/nucleus_decomposition/hierarchy.cpp
        utility/nucleus_decomposition/k3-4.cpp
        utility/nucleus_decomposition/kcore.cpp
        utility/nucleus_decomposition/ktruss.cpp
        utility/nucleus_decomposition/util.cpp
        utility/nucleus_decomposition/nd_interface.cpp
        utility/nucleus_decomposition/nd_interface.h
)
add_library(rapid SHARED
        ${RAPIDMATCH_SRC})
