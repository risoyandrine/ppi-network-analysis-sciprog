import requests 
import pandas as pd

def fetch_string_data(gene, threshold, species, network_type, limit):
    url = "https://string-db.org/api/json/network"
    # based on the STRING database API documentation I chose this URL, returns JSON data 

    # we define the parameters needed for the API request based on the STRING documentation
    params = {
        "identifiers": gene,
        "required_score": threshold,
        "species": species, 
        "network_type": network_type,
        "limit": limit, 
    }

    # we make the request to the STRING database API
    response = requests.get(url, params=params)

    if response.status_code == 200:
        print("Data fetched successfully!\n") 
        data = response.json()
        return pd.DataFrame(data)
    else: 
        print(f"Error: Could not fetch data from STRING database: {response.status_code}")
        return None

def fetch_go_enrich(hub_proteins, species, background=None):
    url = "https://string-db.org/api/json/enrichment"
    # based on the STRING database API documentation I chose this URL, returns JSON data 

    proteins = [proteins for proteins, score in hub_proteins]

    # we define the parameters needed for the API request based on the STRING documentation
    params = {
        "identifiers": "\r".join(proteins),  # STRING expects newline-separated identifiers
        "species": species,
    }
    # we will only include background if provided (custom background set for enrichment)
    if background is not None:
        params["background_string_identifiers"] = "\r".join(background)
    # we make the request to the STRING database API
    response = requests.post(url, data=params)  # STRING enrichment endpoint requires POST

    if response.status_code == 200: 
        data = response.json()
        return pd.DataFrame(data)
    else: 
        print(f"Error: Could not fetch data from STRING database: {response.status_code}")
        return None
