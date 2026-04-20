import requests 
import pandas as pd

#the function below will give us the interaction data from the STRING database, the parameters follows the STRING documentation 
def fetch_string_data(gene, threshold, species, network_type, limit):
    url = "https://string-db.org/api/json/network"

    params = {
        "identifiers": "\r".join(gene),
        "required_score": threshold,
        "species": species, 
        "network_type": network_type,
        "limit": limit, 
    }

    response = requests.get(url, params=params)

    if response.status_code == 200:
        print("Data fetched successfully!\n") 
        data = response.json()
        return pd.DataFrame(data)
    else: 
        print(f"Error: Could not fetch data from STRING database: {response.status_code}")
        return None

#to later perform the enrichment analysis, the hub proteins will be converted to their gene identifiers
def fetch_gene_id(protein_list , species):
    url = "https://string-db.org/api/json/get_string_ids"

    params = {
        "identifiers": "\r".join(protein_list),
        "species": species,
    }
    response = requests.get(url, params=params)

    if response.status_code == 200: 
        data = response.json()
        return pd.DataFrame(data)
    else: 
        print(f"Error: Could not fetch data Gene IDs: {response.status_code}")
        return None

#this function will give the enrichment analysis by sending the previously identified hub proteins to the STRING database
def fetch_go_enrich(hub_proteins, species, background=None):
    url = "https://string-db.org/api/json/enrichment"

    proteins = [proteins for proteins, score in hub_proteins]

    params = {
        "identifiers": "\r".join(proteins), 
        "species": species,
    }
    
    if background is not None:
        params["background_string_identifiers"] = "\r".join(background)

    response = requests.post(url, data=params) 

    if response.status_code == 200: 
        data = response.json()
        return pd.DataFrame(data)
    else: 
        print(f"Error: Could not fetch data from STRING database: {response.status_code}")
        return None
