#next task is to perform an go enrichment analysis of the hub proteins that was found in the network analysis

from fetch_data import fetch_go_enrich

# this function takes the raw dataframe and checks that it is there and not empty before it filters the data based on the false discovery rate threshold. If no significant terms are found it will return none
def go_enrichment(hub_proteins, species, fdr_threshold, background=None):
    go_data = fetch_go_enrich(hub_proteins, species, background)
    if go_data is None:
        print("Error: Could not fetch GO enrichment data from STRING database")
        return None
    go_data = go_data[go_data["fdr"] < fdr_threshold]
    go_data = go_data[go_data["category"].isin(["Process", "Function", "Component"])]
    if go_data.empty:
        print("No GO enrichment data found for the given hub proteins and FDR threshold")
        return None
    return go_data

#to give us a more structured picture of what the enrichment analysis found, it will create a summary of the data
def go_summary(go_data):
    print("GO Enrichment Analysis:\n")
    print(f"Number of significant GO terms found: {len(go_data)}\n")
    for category in ["Process", "Function", "Component"]:
        subset = go_data[go_data["category"] == category]
        if subset.empty:
            continue
        print(f"Biological {category}:\n")
        for _, row in subset.head(10).iterrows():
            print(f"{row['description']} ({row['term']}): FDR = {row['fdr']:.4f}")
        print()

