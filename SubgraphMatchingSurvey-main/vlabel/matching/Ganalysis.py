import networkx as nx
from collections import Counter, defaultdict
import numpy as np
Dataset="patents"
QuerN = "../../../../../FairSM/dataset/"
#QuerN1 = "../../../../../Pilos-Subgraph_Matching/dataset/"
QuerN1 = "../../../../Pilos-Subgraph_Matching/dataset/"
QuerN=QuerN1
Data_graph=QuerN+Dataset+"/data_graph/"+Dataset+".graph"

QS="16"
QNR="1"
Query=QuerN1+Dataset+"/query_graph/query_"+"G"+"_"+QS+"_"+QNR+".graph"
print(Query)
def generate_query_paths(base_path, dataset, qnr):
    query_types = ["cycle", "flat_narrow", "path", "tall", "tall_wide", "wide_flat"]
    return [
    f"../dataset/{dataset}/query_graph/{qtype}_{qs}.txt"
    for qtype in query_types
    for qs in range(20)
    ]

def read_undirected_graph(filepath):
    G = nx.Graph()
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
        for line in lines[1:]:  # Skip the first line (t ...)
            parts = line.strip().split()
            if not parts:
                continue

            if parts[0] == 'v':
                node_id = int(parts[1])
                x = int(parts[2])
                y = int(parts[3])
                G.add_node(node_id, x=x, y=y)

            elif parts[0] == 'e':
                u = int(parts[1])
                v = int(parts[2])
                G.add_edge(u, v)

    return G
def classify_query_graph(G):
    n = G.number_of_nodes()
    m = G.number_of_edges()
    NC=True
    for node in G.nodes():
        if G.degree(node) != 2:
            NC= False    
    if nx.is_connected(G) and nx.cycle_basis(G) and NC:
       
        if n == m:  # exactly one cycle
            return 'cycle'
    elif all(deg <= 2 for _, deg in G.degree()):
        return 'path'
    elif nx.is_tree(G):
        depth = nx.diameter(G)
        avg_branching = sum(dict(G.degree()).values()) / n
        if depth >= 3 and avg_branching <= 2:
            return 'tall_tree'
        elif depth <= 2 and avg_branching > 2:
            return 'wide_flat_tree'
        elif depth >= 3 and avg_branching > 2:
            return 'tall_wide_tree'
        elif depth <= 2 and avg_branching <= 2:
            return 'flat_narrow_tree'
        else:
            return 'unclassified_tree'
    
    else:
        depth = nx.diameter(G)
        avg_branching = sum(dict(G.degree()).values()) / n
        if depth >= 3 and avg_branching <= 2:
            return 'tall_graph'
        elif depth <= 2 and avg_branching > 2:
            return 'wide_flat_graph'
        elif depth >= 3 and avg_branching > 2:
            return 'tall_wide_graph'
        elif depth <= 2 and avg_branching <= 2:
            return 'flat_narrow_graph'
        else:
            return 'unclassified__graph'
        return 'graph'



def degree_statistics(G):
    degrees = [deg for _, deg in G.degree()]
    stats = {
        'num_nodes': G.number_of_nodes(),
        'num_edges': G.number_of_edges(),
        'avg_degree': np.mean(degrees),
        'min_degree': np.min(degrees),
        'max_degree': np.max(degrees),
        'std_degree': np.std(degrees),
        #'degree_distribution': dict(Counter(degrees))
    }
    return stats
def analyze_all_degrees(graph_list):
    all_stats = []
    for i, G in enumerate(graph_list):
        stats = degree_statistics(G)
        stats['graph_index'] = i
        all_stats.append(stats)
    return all_stats
def compute_graph_stats(G):
    stats = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'avg_degree': sum(dict(G.degree()).values()) / G.number_of_nodes(),
        'has_cycle': not nx.is_tree(G),
        'diameter': nx.diameter(G) if nx.is_connected(G) else None,
        'tree_depth': (
            nx.diameter(nx.minimum_spanning_tree(G))
            if nx.is_connected(G) else None
        )
    }
    return stats
def analyze_graph_structure(G):
    print("===== Graph Structural Properties =====")
    
    # Basic info
    #print("Type:", "Directed" if G.is_directed() else "Undirected")
    #print("Number of nodes (order):", G.number_of_nodes())
    #print("Number of edges (size):", G.number_of_edges())

    # Degree info
    degrees = dict(G.degree())
    #print("Degrees:", degrees)
    avg_degree = sum(degrees.values()) / len(degrees)
    print("Average Degree:", avg_degree)

    # Connectivity
    if G.is_directed():
        is_connected = nx.is_strongly_connected(G)
        n_components = nx.number_strongly_connected_components(G)
    else:
        is_connected = nx.is_connected(G)
        n_components = nx.number_connected_components(G)

    #print("Is connected:", is_connected)
    #print("Number of connected components:", n_components)

    # Path-based metrics (if connected)
    if is_connected:
        if G.is_directed():
            largest_cc = max(nx.strongly_connected_components(G), key=len)
        else:
            largest_cc = max(nx.connected_components(G), key=len)
        subgraph = G.subgraph(largest_cc)
        print("Diameter:", nx.diameter(subgraph))
        print("Radius:", nx.radius(subgraph))
        print("Average shortest path length:", nx.average_shortest_path_length(subgraph))

    # Clustering
    #if not G.is_directed():
        #print("Clustering Coefficient (per node):", nx.clustering(G))
        #print("Average Clustering Coefficient:", nx.average_clustering(G))

    # Density
    print("Graph Density:", nx.density(G))

    # Degree distribution
    #degree_sequence = [d for n, d in G.degree()]
    #print("Degree Distribution:", degree_sequence)

    print("=======================================\n")
def analyze_graphs_structureGS(G_set):
    #print("===== Aggregate Graph Structural Properties =====")

    total_nodes = 0
    total_edges = 0
    total_degree = 0
    total_density = 0
    count_connected = 0

    diameters = []
    radii = []
    avg_shortest_paths = []

    for G1 in G_set:
        G=read_undirected_graph(G1)
        degrees = dict(G.degree())
        total_nodes += G.number_of_nodes()
        total_edges += G.number_of_edges()
        total_degree += sum(degrees.values())
        total_density += nx.density(G)

        # Check connectivity and collect path metrics
        if G.is_directed():
            is_connected = nx.is_strongly_connected(G)
            if is_connected:
                largest_cc = max(nx.strongly_connected_components(G), key=len)
                subgraph = G.subgraph(largest_cc)
                diameters.append(nx.diameter(subgraph))
                radii.append(nx.radius(subgraph))
                avg_shortest_paths.append(nx.average_shortest_path_length(subgraph))
                count_connected += 1
        else:
            is_connected = nx.is_connected(G)
            if is_connected:
                largest_cc = max(nx.connected_components(G), key=len)
                subgraph = G.subgraph(largest_cc)
                diameters.append(nx.diameter(subgraph))
                radii.append(nx.radius(subgraph))
                avg_shortest_paths.append(nx.average_shortest_path_length(subgraph))
                count_connected += 1

    num_graphs = len(G_set)
    avg_degree = total_degree / total_nodes if total_nodes > 0 else 0
    avg_density = total_density / num_graphs if num_graphs > 0 else 0
    avg_diameter = sum(diameters) / len(diameters) if diameters else None
    avg_radius = sum(radii) / len(radii) if radii else None
    avg_shortest_path = sum(avg_shortest_paths) / len(avg_shortest_paths) if avg_shortest_paths else None

    #print(f"Number of graphs: {num_graphs}")
    #print(f"Total nodes(AVG): {total_nodes/num_graphs}")
    #print(f"Total edges(avg): {total_edges/num_graphs}")
    print(f"Average Degree: {avg_degree:.4f}")
    print(f"Graph Density (avg): {avg_density:.4f}")
    #print(f"Connected Graphs: {count_connected}/{num_graphs}")
    print(f"Average Diameter (connected only): {avg_diameter}")
    print(f"Average Radius (connected only): {avg_radius}")
    print(f"Average Shortest Path Length (connected only): {avg_shortest_path}")
    print("==============================================\n")
def analyze_query_workload(query_graphs):
    type_counter = Counter()
    type_stats = defaultdict(list)

    for G1 in query_graphs:
        G=read_undirected_graph(G1)
        gtype = classify_query_graph(G)
        if (gtype=="flat_narrow_tree"):
            print(G1)
        type_counter[gtype] += 1
        type_stats[gtype].append(compute_graph_stats(G))

    # Aggregate stats
    aggregated = {}
    for gtype, stats_list in type_stats.items():
        if not stats_list:
            continue
        keys = stats_list[0].keys()
        avg_stats = {
            key: sum(s[key] for s in stats_list if s[key] is not None) / len(stats_list)
            for key in keys
        }
        aggregated[gtype] = {
            'count': type_counter[gtype],
            #'avg_stats': avg_stats
        }

    return aggregated
print(read_undirected_graph(Query))
basepath=QuerN1
#print(degree_statistics(read_undirected_graph(Data_graph)))
paths = generate_query_paths(basepath, Dataset, QS)
stats = analyze_query_workload(paths)
simplified = {k: v['count'] for k, v in stats.items()}
print(simplified)

analyze_graphs_structureGS(paths)