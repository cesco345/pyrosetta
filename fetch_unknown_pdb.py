# Cell 2: Function to fetch PDB info
import requests
import os 

def fetch_pdb_info(search_term):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                "operator": "exact_match",
                "value": search_term
            }
        },
        "return_type": "entry"
    }
    
    response = requests.post(url, json=query)
    response.raise_for_status()
    
    data = response.json()
    if data['total_count'] > 0:
        return data['result_set'][0]['identifier']
    else:
        return None

# Cell 3: Function to download PDB file
def download_pdb(pdb_id, output_dir='.'):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    response.raise_for_status()
    
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    with open(output_path, 'wb') as f:
        f.write(response.content)
    
    return output_path
# Cell 4: Fetch and download PDB
search_term = "AraC"  # You can change this to the specific gene name if needed
try:
    pdb_id = fetch_pdb_info(search_term)
    if pdb_id:
        pdb_path = download_pdb(pdb_id)
        print(f"Downloaded PDB file: {pdb_path}")
    else:
        print(f"No PDB found for {search_term}")
except requests.exceptions.RequestException as e:
    print(f"Error fetching or downloading PDB: {e}")

