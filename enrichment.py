#we want to perform an go enrichment analysis of the hub proteins obtained

from fetch_data import fetch_go_enrich

def go_enrichment(hub_proteins, species, fdr_threshold, background=None):
    #get the raw dataframe 
    go_data = fetch_go_enrich(hub_proteins, species, background)
    #if none, print error and return
    if go_data is None:
        print("Error: Could not fetch GO enrichment data from STRING database")
        return None
    #filter rows based on the fdr threshold
    go_data = go_data[go_data["fdr"] < fdr_threshold]
    #filter to only include actual GO terms (Process, Function, Component)
    go_data = go_data[go_data["category"].isin(["Process", "Function", "Component"])]
    #return the filtered data
    if go_data.empty:
        print("No GO enrichment data found for the given hub proteins and FDR threshold")
        return None
    return go_data

#to make the data more readable, we create a summary of the GO enrichment data
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

