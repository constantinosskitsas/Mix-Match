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
We build our algorithm using SIGMOD'2024 paper A Comprehensive Survey and Experimental Study of Subgraph Matching: Trends, Unbiasedness by Dr. Shixuan 
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

```zsh
./SubgraphMatching.out -d data_graphs -q query_graphs
-filter method_of_filtering_candidate_vertices -order method_of_ordering_query_vertices -engine method_of_enumerating_partial_results -num number_of_embeddings,
```

For detailed parameter settings, see `matching/matchingcommand.h`.

**Example 1**: This approach uses the filtering and ordering methods of GraphQL to generate candidate vertex sets and determine the matching order. Results are enumerated using the set-intersection-based local candidate computation method.


```zsh
./SubgraphMatching.out -d ../../dataset/dblp/data_graph/dblp.graph -q ../../dataset/dblp/query_graph/query_G_32_1.graph -filter VEQ -order VEQ -engine MM -num -1 -symmetry 1 -FairT 2 -time 1 -SF saveM
```

**Example2**: (Use the filtering method CaLiG to generate the candidate vertex sets, ordering method of RI to generate the matching order,  KSS engine to enumerate the results with automorphic graphs detection):

```zsh
./SubgraphMatching.out -dataset dblp -qsize 16 -qnumber 1 -qnumberL 10 -qprop G -SF magkas -FairT 2 -filter GQL -order GQL -engine LFTJ -num 10 -symmetry 1 -time 1
```


```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter CaLiG -order RI -engine KSS -num MAX -symmetry 1
```


## Input

The input format of data graph and query graph is fixed which is the same as the previous study [2]. Each graph file starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted as 'v VertexID LabelId Degree' and 'e VertexId VertexId (EdgeLabel)?' respectively. The vertex id must start from 0 and the range is [0,N - 1]. The following is an input sample. You can also find sample data sets and query sets under the test folder. For single graph the format of edges must be the same.

```
t <vertex_number> <edge number>
v <vertex_id> <vertex_label> <vertex_degree>
e <source> <target> (<edge_label>)?
```

For example, the input graph without edge label could be:

```zsh
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

A sample input graph with edge labels is as follows, where the third column of the edge line is the edge label (e.g. the edge label of the first edge is 0).

```zsh
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1 0
e 0 2 1
e 1 2 2
e 1 3 3
e 2 4 4
e 3 4 5
```


## Techniques Supported

For each technique, we support 10 methods, as shown in the following table.

| -filter | -order | -engine |
| :-----: | :----: | :-----: |
|   LDF   |  QSI   |   QSI   |
|   NLF   |  VF3   |   VF3   |
|   GQL   |  GQL   |   GQL   |
|   TSO   |  TSO   | EXPLORE |
|   CFL   |  CFL   |  LFTJ   |
|  DPiso  | DPiso  |  DPiso  |
|  CECI   |  CECI  |  CECI   |
|   VEQ   |   RI   |   VEQ   |
|   RM    |   RM   |   RM    |
|  CaLiG  | VF2PP  |   KSS   |

## Configuration

In addition to setting the filtering, ordering and enumeration methods, you can also configure the set intersection algorithms and optimization techniques by defining macros in 'configuration/config.h'. Please see the comments in 'configuration/config.h' for more details. We briefly introduce these methods in the following.


|         Macro         |                         Description                          |
| :-------------------: | :----------------------------------------------------------: |
|       HYBRID 0        | a hybrid method handling the cardinality skew by integrating the merge-based method with the galloping-based method |
|       HYBRID 1        |               the merge-based set intersection               |
|         SI 0          |          Accelerate the set intersection with AVX2           |
|         SI 1          |         Accelerate the set intersection with AVX512          |
|         SI 2          |                   Scalar set intersection                    |
|    ENABLE_QFLITER     | the [QFilter](https://dl.acm.org/doi/10.1145/3183713.3196924) set intersection algorithm |
|  ENABLE_FAILING_SET   |         the failing set pruning technique in DP-iso          |
| ENABLE_EQUIVALENT_SET |      enable the equivalent set pruning technique in VEQ      |
|    ELABELED_GRAPH     |                  enable edge labeled graph                  |

> [1] Zhijie Zhang, Yujie Lu, Weiguo Zheng, and Xuemin Lin. 2024. A Comprehensive Survey and Experimental Study of Subgraph Matching: Trends, Unbiasedness, and Interaction. Proc. ACM Manag. Data 2, 1, Article 60 (February 2024), 29 pages.
> 
> [2] Shixuan Sun and Qiong Luo. 2020. In-Memory Subgraph Matching: An In-depth Study. In Proceedings of the 2020 ACM SIGMOD International Conference on Management of Data (SIGMOD '20). Association for Computing Machinery, New York, NY, USA, 1083–1098.
