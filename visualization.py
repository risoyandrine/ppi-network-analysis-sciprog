# we import the necessary libraries
from math import log10
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
import networkx as nx

# because pyvis is not a standard library, we use a try-except block
try:
    from pyvis.network import Network
except ImportError:
    Network = None

#setting the color palette for the network graph and the GO enrichment plot

network_palette = {
    "hub": "lightcoral",
    "node": "thistle",
    "edge": "lightgrey",
    "label": "saddlebrown",
    "background": "oldlace",
}

go_palette = {
    "Process": "steelblue",
    "Function": "powderblue",
    "Component": "deepskyblue",
}

#we define the function to plot the network, we do it in two different ways: static and interactive
#the static plot is a png image, the interactive plot is an html file

def plot_network(graph, network_properties, output_path):
    outputs = []
    static_path = plot_network_static(graph, network_properties, output_path) # define the static plot path
    outputs.append(static_path)

    if Network is not None:
        interactive_path = plot_network_interactive(graph, network_properties, output_path) # define the interactive plot path
        outputs.append(interactive_path)

    return outputs

#we define the interactive plot function, we use pyvis to create an interactive network graph, if pyvis is not installed it will not be created

def plot_network_interactive(graph, network_properties, output_path):
    net = Network(height="1200px", width="100%", bgcolor="oldlace", font_color="saddlebrown")
    net.from_nx(graph)
# we use the data from the network_properties dictionary from the network_analysis.py file to create the interactive plot
    degree = network_properties["degree"]
    betweenness = network_properties["betweenness"]
    clustering = network_properties["clustering"]
    hub_proteins = network_properties["hub_proteins"]

# we define the hub names and hub scores from the hub_proteins list
    hub_names = {protein for protein, score in hub_proteins}
    hub_scores = {protein: score for protein, score in hub_proteins}

# we iterate over the nodes in the network graph
    for node in net.nodes:
        protein = node["id"]
        degree_score = degree.get(protein, 0) # get the degree score for the current node, with 0 as default value
        betweenness_score = betweenness.get(protein, 0) # get the betweenness score for the current node, with 0 as default value
        clustering_score = clustering.get(protein, 0) # get the clustering score for the current node, with 0 as default value

        node["size"] = 18 + degree_score * 32 # set the size of the node based on the degree score
        node["font"] = {"size": 18, "face": "Arial", "color": network_palette["label"]}
        node["title"] = (
            f"{protein}"
            f"\nDegree centrality: {degree_score:.3f}" #
            f"\nBetweenness centrality: {betweenness_score:.3f}"
            f"\nClustering coefficient: {clustering_score:.3f}"
        )

        if protein in hub_names:
            score = hub_scores[protein]
            node["color"] = network_palette["hub"]
            node["borderWidth"] = 3
            node["title"] += f"\nHub score: {score:.3f}"
        else:
            node["color"] = network_palette["node"]
# we iterate over the edges in the network graph
    for edge in net.edges:
        weight = edge.get("weight") or edge.get("width") or edge.get("value") or 1
        edge["value"] = max(weight, 1)
        edge["color"] = network_palette["edge"]
        edge["title"] = f"Interaction score: {weight:.3f}"
# we use the barnes_hut algorithm to layout the network graph
    net.barnes_hut(gravity=-22000, central_gravity=0.2, spring_length=140, spring_strength=0.03, damping=0.95)
    output_path = Path(output_path).with_suffix(".html")
    net.write_html(str(output_path))
    
# pyvis template has a bug with double headers, so we inject our header manually:
    with open(output_path, "r", encoding="utf-8") as f:
        html_str = f.read()
    html_str = html_str.replace("<h1></h1>", "<h1>Protein-Protein Interaction Network</h1>", 1)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_str)
        
    return output_path

# we define the static plot function
def plot_network_static(graph, network_properties, output_path):
    degree = network_properties["degree"]
    hub_proteins = network_properties["hub_proteins"]
    hub_names = {protein for protein, _ in hub_proteins}

    fig, ax = plt.subplots(figsize=(11, 8.5))
    fig.patch.set_facecolor(network_palette["background"])
    ax.set_facecolor(network_palette["background"])

    pos = nx.spring_layout(graph, k=2, iterations=200, seed=42) # define the layout of the network graph
    node_sizes = [500 + degree.get(node, 0) * 2800 for node in graph.nodes()] # set the size of the nodes based on the degree score
    node_colors = [network_palette["hub"] if node in hub_names else network_palette["node"] for node in graph.nodes()]


    nx.draw_networkx_edges(graph, pos=pos, ax=ax, edge_color=network_palette["edge"], width=1.4, alpha=0.8)
    nx.draw_networkx_nodes(graph, pos=pos, ax=ax, node_size=node_sizes, node_color=node_colors, linewidths=1.3, edgecolors="white")
    nx.draw_networkx_labels(graph, pos=pos, ax=ax, font_size=9, font_weight="bold", font_color=network_palette["label"])

    ax.set_title("Protein-Protein Interaction Network", fontsize=15, fontweight="bold", color="saddlebrown", pad=16)
    hub_label = ", ".join(protein for protein, _ in hub_proteins[:5]) or "None"
    summary_text = f"Nodes: {network_properties['node_count']}\nEdges: {network_properties['edge_count']}\nTop hub proteins: {hub_label}"
    ax.text(0.02, 0.02, summary_text, transform=ax.transAxes, fontsize=9, color="saddlebrown", va="bottom", bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "edgecolor": "lightgrey"})

    ax.axis("off")
    plt.tight_layout()

# we save the plot as a png file
    output_path = Path(output_path).with_suffix(".png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return output_path

# we define the function to plot the GO enrichment
# we us matplotlib to create a bar chart of the GO enrichment results to avoid any more dependencies

def plot_GOenrich(go_data, output_png, top_n=10):
    # we check if the go_data is empty
    if go_data is None or go_data.empty:
        raise ValueError("go_data must contain at least one enrichment result")
    # we take the top n results
    plot_data = go_data.nsmallest(top_n, "fdr").copy()
    plot_data["score"] = plot_data["fdr"].map(lambda value: -log10(max(value, 1e-300)))
    plot_data = plot_data.sort_values("score")
    plot_data["color"] = plot_data["category"].map(go_palette).fillna("grey")
    plot_data["short_label"] = plot_data.apply(lambda row: format_enrichment_label(row["description"], row["category"]), axis=1)

    # we define the figure and the axes

    fig_height = max(5.5, 1.0 * len(plot_data))
    fig, ax = plt.subplots(figsize=(11.5, fig_height))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.barh(plot_data["short_label"], plot_data["score"], color=plot_data["color"], edgecolor="white", linewidth=1.0)
    ax.set_xlabel("-log10(FDR)", fontsize=11, color="saddlebrown")
    ax.set_ylabel("Enriched terms", fontsize=11, color="saddlebrown")
    ax.set_title("Enrichment Results for Hub Proteins", fontsize=15, fontweight="bold", color="saddlebrown", pad=16)
    ax.grid(axis="x", linestyle="--", linewidth=0.7, color="lightgrey", alpha=0.9)
    ax.set_axisbelow(True)

# we remove the top and right spines
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("lightgrey")
    ax.spines["bottom"].set_color("lightgrey")
    ax.tick_params(axis="both", labelsize=9, colors="saddlebrown")

    for index, value in enumerate(plot_data["score"]):
        ax.text(value + 0.05, index, f"{value:.2f}", va="center", fontsize=8.5, color="saddlebrown")
    
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color, label=category) for category, color in go_palette.items() if category in set(plot_data["category"])]
    if legend_handles:
        ax.legend(handles=legend_handles, title="Category", frameon=False, loc="lower right")

    plt.tight_layout()
# we save the plot as a png file
    output_path = Path(output_png)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return output_path

# we define the function to format the enrichment label
def format_enrichment_label(description, category, max_length=58):
    compact_description = " ".join(str(description).split())
    if len(compact_description) > max_length:
        compact_description = compact_description[: max_length - 3].rstrip() + "..."
    return fill(f"{compact_description} [{category}]", width=42)