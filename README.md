# Protein-Protein Interaction (PPI) Network Analysis

This is my Python project for analyzing protein-protein interaction networks (Project 7). The program takes one or more gene names, fetches interaction data from the STRING database, and builds a network graph. In addition, it runs a Gene Ontology (GO) enrichment analysis to identify what biological processes and functions the proteins in the network may be involved in.

The structure of the project includes splitting the code into different files as you see below:
- main.py: This is the main script you run from the terminal.
- fetch_data.py: Here it connects to the STRING API to download PPI data and GO enrichment data.
- network.py: This scrip then builds the network graph and calculates a variety of network properties including degree centrality, betweenness centrality, and clustering coefficients to identify important "hub" proteins.
- enrichment.py: and filters the enrichment data so we only look at actual GO terms (Component, Process, and Function). It also supports custom background sets if needed.
- visualization.py: Creates the final network images and the bar chart, this also includes an interactive network graph where you can visualize the interactions.

## Requirements:
To run this project you will need to install a few python libraries. You can install these by running the line of code below:
```bash
pip install networkx matplotlib pyvis pandas requests
```

## How to run the program:
You may run the program from the terminal using `main.py`. Here you have to provide one or more gene names and a score threshold for the network connections. Just to provide an example below, the gene TP53 is used and the threshold is set to 400. 

**Example with a single gene:**
```bash
python main.py --gene TP53 --threshold 400
```

**Example with multiple genes:**
```bash
python main.py --gene TP53 BRCA1 EGFR --threshold 400
```

## Arguments
Some other arguments in this program are optional and can be changed as wanted by the user to fit their individual goals and needs:

| Argument | Required | Default | Description |
|---|---|---|---|
| --gene | Yes | - | One or more gene symbols (for example: TP53 BRCA1 EGFR) |
| --threshold | Yes | - | Interaction score (ranging from 0-1000) |
| --species | No | 9606 | NCBI taxonomy ID |
| --limit | No | 10 | The maximum number of interactions |
| --num_hubs | No | 10 | The number of hub proteins |
| --fdr_threshold | No | 0.05 | FDR cutoff for GO enrichment |
| --output | No | ppi_network | Output filename for the produced network plots |
| --go_output | No | go_enrichment.png | Output filename for the produced GO plot |

## Output files
When the script is finished, it automatically generates three visual files in your folder:
1. ppi_network.png: Here you will see a static image of the network graph. Hub proteins are colored differently.
2. ppi_network.html: This is an interactive webpage where you can physically drag the proteins around. In addition, you may also hold your mouse over a node to see its calculated centralities and hub score.
3. go_enrichment.png: This plot is a bar chart showing the top 10 most significant GO terms for the network.

---
**Author:** Andrine Risøy
**Date:** April 2026