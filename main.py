import argparse
from pathlib import Path

from enrichment import go_enrichment, go_summary
from fetch_data import fetch_string_data, fetch_gene_id
from network import (calc_betweenness_centr, calc_clustering_coefficient, calc_degree_centrality, create_graph, find_hub_proteins, get_network_properties, network_summary)
from visualization import plot_GOenrich, plot_network

# we parse terminal argument
arg = argparse.ArgumentParser()

arg.add_argument("--gene", type=str, required=True)
arg.add_argument("--threshold", type=float, required=True)
arg.add_argument("--species", type=int, default=9606)
arg.add_argument("--network_type", type=str, default="functional")
arg.add_argument("--limit", type=int, default=10)
arg.add_argument("--num_hubs", type=int, default=10)
arg.add_argument("--fdr_threshold", type=float, default=0.05)
arg.add_argument("--output", type=str, default="ppi_network", help="Output filename for the network visualization")
arg.add_argument("--go_output", type=str, default="go_enrichment.png", help="Output filename for the GO enrichment plot")

# we define a function to parse the arguments
def parse_arguments():
    args = arg.parse_args()
    return (args.gene, args.threshold, args.species, args.network_type, args.limit, args.num_hubs, args.fdr_threshold, args.output, args.go_output)

# we define the main function, this is where we call the other functions in the correct order and save the outputs to the desired files
def main():
    (gene, threshold, species, network_type, limit, num_hubs, fdr_threshold, output, go_output) = parse_arguments()
    data = fetch_string_data(gene, threshold, species, network_type, limit)

    if data is None or data.empty:
        print("No interaction data was returned from STRING.")
        return

    graph = create_graph(data)
    degree = calc_degree_centrality(graph)
    betweenness = calc_betweenness_centr(graph)
    clustering = calc_clustering_coefficient(graph)
    hub_proteins = find_hub_proteins(degree, num_hubs=num_hubs)
    network_properties = get_network_properties(graph, degree, betweenness, clustering, hub_proteins)

    hub_names = [protein for protein, score in hub_proteins]
    hub_gene_ids = fetch_gene_id(hub_names, species)

    if hub_gene_ids is not None and not hub_gene_ids.empty:
        print("\nProtein to Gene ID mapping:")
        hub_gene_ids["annotation"] = hub_gene_ids["annotation"].str.split(";").str[0]
        print(hub_gene_ids[["stringId", "preferredName", "annotation"]].to_string(index=False))
    

    network_summary(graph, degree, betweenness, clustering, hub_proteins)

    output_path = Path(output)
    saved_paths = plot_network(graph, network_properties, output_path)
    print("\nSaved network visualizations:")
    for saved_path in saved_paths:
        print(saved_path)

    go_data = go_enrichment(hub_proteins, species, fdr_threshold)
    if go_data is None:
        print("\nNo GO enrichment plot was generated.")
        return

    print()
    go_summary(go_data)
    go_output_path = Path(go_output)
    saved_go_path = plot_GOenrich(go_data, go_output_path)
    print(f"Saved GO enrichment plot to {saved_go_path}")

#we run the main function if the script is run directly
if __name__ == "__main__":
    main()
