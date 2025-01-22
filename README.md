# Mix & Match: Subgraph Matching for Absolute Coverage
The NP-hard problem ofsubgraph matching calls to return all matchings of a smaller query graph within a larger data graph. The problem is fundamental in graph analysis and query answering, as identifying the location of all matchings facilitates the understanding
and analysis of the larger graph. Nevertheless, existing subgraph
matching methods return results from one location of the graph
before moving to another location, while the total results may be
in the order of billions or even trillions; under these circumstances,
existing methods may only present a portion of the results within
reasonable time or space, which is not representative of the totality
of results. This predicament leads to a biased representation of
the data graph. In this paper, we study the problem of diversity in
subgraph matching and propose an algorithm that quickly returns
results that are representative of the whole data graph.
## **Code**
We build our algorithm using SIGMOD'2024 paper [A Comprehensive Survey and Experimental Study of Subgraph Matching: Trends, Unbiasedness and Interaction](https://dl.acm.org/doi/pdf/10.1145/3639315) 
We kept all the functionalities of the framework and for more specific details we refer to [Code](https://github.com/RapidsAtHKUST/SubgraphMatching).
## General Information
We alter the supported algorighms (enumeration method) to return coverage results.

## Supported Algorithms
|Algorithm|Description|Execution code
|:--------:|:------------:|:------------:
|M&M | Mix & Match | MM
|M&M-I | Mix & Match with Initialization | MMI
|DSQL | Diversified Subgraph Query | DSQL
|DSQL+ | Optimized Diversified Subgraph Query | DSQLP
|VEQ | Versatile Equivalences | VEQ
|LFTJ | the set-intersection based local candidates computation | LFTJ
|RM | Rapid-Match | RM
## Supported Algorithms-TOPK
|Algorithm|Description|Execution code
|:--------:|:------------:|:------------:
|M&M | Mix & Match | MMK
|M&M-I | Mix & Match with Initialization | MMIK
|DSQL | Diversified Subgraph Query | DSQLK
|DSQL+ | Optimized Diversified Subgraph Query | DSQLPK
|LFTJ | the set-intersection based local candidates computation | LFTJK
## Global Exploration
Only for M&M algorithm
|Algorithm|Description|Execution code
|:--------:|:------------:|:------------:
|- | none | 0
|U | Uncovered | 1
|mN | Minimum covered Neighbors | 2
|MN | Maximum uncovered Neighbors (𝑀𝑁) | 3
|HN | Hybrid | 4
## build
Within the `vlabel` directory, create a build directory and compile the source code.
```zsh
mkdir build & cd build
cmake ..
make
```

## Execute

Execute the binary with the following command.
./SubgraphMatching.out -dataset dblp -qsize 16 -qnumber 1 -qnumberL 1 -qprop G -order CFL -filter VEQ -engine MM  -SF DIV-BASIC -time 1 -FairT 0 -symmetry 1
```zsh
./SubgraphMatching.out -dataset data_graphs -qsize query size
-filter method_of_filtering_candidate_vertices -order method_of_ordering_query_vertices -engine method_of_enumerating_partial_results -FairT global diversity method, -symmetry authomorphisms -SF save file name -time Time limit
```

For detailed parameter settings, see `matching/matchingcommand.h`.

**Example 1**: MM experiment setup for 1 sec time limit.

```zsh
./SubgraphMatching.out -d ../../dataset/dblp/data_graph/dblp.graph -q ../../dataset/dblp/query_graph/query_G_32_1.graph -filter VEQ -order CFL -engine MM -num -1 -symmetry 1 -FairT 2 -time 1 -SF Coverage
```
**Example 2**: VEQ experiment setup.

```zsh
./SubgraphMatching.out -d ../../dataset/dblp/data_graph/dblp.graph -q ../../dataset/dblp/query_graph/query_G_32_1.graph -filter VEQ -order VEQ -engine VEQ -num -1 -symmetry 1 -FairT 2 -time 1 -SF Coverage
```


**Example 3**: LFTJ experiment setup.

```zsh
./SubgraphMatching.out -d ../../dataset/dblp/data_graph/dblp.graph -q ../../dataset/dblp/query_graph/query_G_32_1.graph -filter VEQ -order CFL -engine LFTJ -num -1 -symmetry 1 -time 1 -SF Coverage
```

**Example 4**: RM experiment setup.

```zsh
./SubgraphMatching.out -d ../../dataset/dblp/data_graph/dblp.graph -q ../../dataset/dblp/query_graph/query_G_32_1.graph -filter RM -order RM -engine RM -num -1 -symmetry 1 -time 1 -SF Coverage
```

> [1] Zhijie Zhang, Yujie Lu, Weiguo Zheng, and Xuemin Lin. 2024. A Comprehensive Survey and Experimental Study of Subgraph Matching: Trends, Unbiasedness, and Interaction. Proc. ACM Manag. Data 2, 1, Article 60 (February 2024), 29 pages.
> 
> [2] Shixuan Sun and Qiong Luo. 2020. In-Memory Subgraph Matching: An In-depth Study. In Proceedings of the 2020 ACM SIGMOD International Conference on Management of Data (SIGMOD '20). Association for Computing Machinery, New York, NY, USA, 1083–1098.
