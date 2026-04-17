import networkx as nx

#in this function, the network graph is build from the interaction data we fetched from the STRING database
def create_graph(data): 
    G = nx.Graph()
    for _, row in data.iterrows():
        G.add_edge(row["preferredName_A"], 
            row["preferredName_B"], 
            weight = row["tscore"])
    return G

#a few key properties are then calculated to give us an insight of the network structure
def calc_degree_centrality(network): 
    degree = nx.degree_centrality(network)
    return degree

#calculate the betweenness centrality that represents the number of shortest paths that pass through a node
def calc_betweenness_centr(network): 
    betweenness = nx.betweenness_centrality(network)
    return betweenness
#the clustering coefficient is a measurement of the degree to which nodes in a graph tend to cluster together
def calc_clustering_coefficient(network): 
    clustering = nx.clustering(network)
    return clustering

#this will define the hub-proteins, which are the proteins with the highest degree centrality, so we use the degree centrality to identify these 
def find_hub_proteins(degree, num_hubs = 10):
   hub_proteins = sorted(degree.items(), key = lambda x: x[1], reverse = True) 
   return hub_proteins[:num_hubs]

#to be able to reuse the values for visualization, we create a dictionary with the values obtained from the calculations
def get_network_properties(graph, degree, betweenness, clustering, hub_proteins):
    network_properties = {
      "node_count": graph.number_of_nodes(),
      "edge_count": graph.number_of_edges(),
      "degree": degree,
      "betweenness": betweenness,
      "clustering": clustering,
      "hub_proteins": hub_proteins
    }
    return network_properties

#then we print a summary of the network properties to summarize the network properties
def network_summary(graph, degree, betweenness, clustering, hub_proteins): 
   print(f"\nNetwork graph for the gene of interest:")
   print(f"\nNumber of nodes(proteins): {graph.number_of_nodes()}")
   print(f"\nNumber of edges(connections): {graph.number_of_edges()}")
   print(f"\nCalculated degree centrality:")
   for protein, score in sorted(degree.items(), key = lambda x: x[1], reverse = True):
      print(f"{protein}: {score:.3f}")
   print(f"\nCalculated betweenness centrality:")
   for protein, score in sorted(betweenness.items(), key = lambda x: x[1], reverse = True):
      print(f"{protein:15}: {score:.3f}")
   print(f"\nCalculated clustering coefficient:")
   for protein, score in sorted(clustering.items(), key = lambda x: x[1], reverse = True):
      print(f"{protein:15}: {score:.3f}")
   print(f"\nHub proteins:")
   for i, (protein, score) in enumerate(hub_proteins, 1):
      print(f"{i}. {protein:15}: {score:.3f}")