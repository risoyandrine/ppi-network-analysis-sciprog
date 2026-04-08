# Protein-Protein Interaction (PPI) Network Analysis

This is my Python project for analyzing protein-protein interaction networks (Project 7). The program takes a gene name, fetches interaction data from the STRING database, and builds a network graph. It also runs a Gene Ontology (GO) enrichment analysis to see what biological processes and functions the proteins in the network are involved in.

## Project Structure
The code is split into several modular files so it is easy to read and extend:
- `main.py`: The main script you run from the terminal.
- `fetch_data.py`: Connects to the STRING API to download PPI data and GO enrichment data.
- `network.py`: Builds the network graph and calculates network properties (degree centrality, betweenness centrality, clustering coefficients) to find important "hub" proteins.
- `enrichment.py`: Cleans and filters the enrichment data so we only look at actual GO terms (Component, Process, and Function). It also supports custom background sets if needed.
- `visualization.py`: Creates the final network images and the bar chart.

## Requirements
You need to install a few python libraries before running the code. Run this in your terminal:
```bash
pip install networkx matplotlib pyvis pandas requests
```

## How to use the program
You run the program from the terminal using `main.py`. You have to provide a gene name and a score threshold for the network connections.

**Example run:**
```bash
python main.py --gene TP53 --threshold 400
```
*(In this example, we use the gene TP53, which is highly connected to cancer biology, and set the interaction score limit to 400. You can also add `--limit` to change how many proteins to get, or `--species` if you are not looking at human data). Below is the arguments that is used in the program:*

## Arguments
| Argument | Required | Default | Description |
|---|---|---|---|
| --gene | Yes | - | Gene symbol (e.g. TP53) |
| --threshold | Yes | - | Interaction score 0-1000 |
| --species | No | 9606 | NCBI taxonomy ID |
| --limit | No | 10 | Max interactions |
| --num_hubs | No | 10 | Number of hub proteins |
| --fdr_threshold | No | 0.05 | FDR cutoff for GO |
| --output | No | ppi_network | Output filename for network plots |
| --go_output | No | go_enrichment.png | Output filename for GO plot |

## Output files
When the script is finished, it automatically generates three visual files in your folder:
1. `ppi_network.png`: A static image of the network graph. Hub proteins are colored differently.
2. `ppi_network.html`: An interactive webpage where you can drag the proteins around. If you hold your mouse over a node, you can see its calculated centralities and hub score!
3. `go_enrichment.png`: A bar chart showing the top 10 most significant Gene Ontology terms for our network, so we know what the network actually DOES biologically (like apoptosis or hypoxia).

---
**Author:** Andrine Risøy
**Date:** April 2026