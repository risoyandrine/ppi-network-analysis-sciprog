#make functions to calculate the key properties
#starting off with calculating degree centrality
# we define the nodes as the proteins and the edges between them as their connection

import networkx as nx

def create_graph(data): 
    G = nx.Graph()
    for _, row in data.iterrows():
        G.add_edge(row["preferredName_A"], 
            row["preferredName_B"], 
            weight = row["tscore"])
    return G

# we calculate the degree centrality, this represents the number of edges connected to a node
def calc_degree_centrality(network): 
    degree = nx.degree_centrality(network)
    return degree

# we calculate the betweenness centrality that represents the number of shortest paths that pass through a node
def calc_betweenness_centr(network): 
    betweenness = nx.betweenness_centrality(network)
    return betweenness
#the clustering coefficient is a measurement of the degree to which nodes in a graph tend to cluster together
def calc_clustering_coefficient(network): 
    clustering = nx.clustering(network)
    return clustering

#we also define the hub-proteins, that is the proteins with the highest degree centrality hence we use the degree centrality to find them
#since degree centrality is a measurement of the number of edges connected to a node: connections/total possible connections - 1, the hub proteins are those with the most connections
def find_hub_proteins(degree, num_hubs = 10):
   hub_proteins = sorted(degree.items(), key = lambda x: x[1], reverse = True) # we sort the degree centrality in descending order
   return hub_proteins[:num_hubs] # define the number of hubs to return

#to reuse the values for visualization, we create a dictionary with the values obtained from the calculations
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

#then we print a summary of the network properties
def network_summary(graph, degree, betweenness, clustering, hub_proteins): 
   print(f"Network graph for the gene of interest:")
   print(f"Number of nodes(proteins): {graph.number_of_nodes()}")
   print(f"Number of edges(connections): {graph.number_of_edges()}")
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
   for i, (protein, score) in enumerate(hub_proteins):
      print(f"{i}. {protein:15}: {score:.3f}")