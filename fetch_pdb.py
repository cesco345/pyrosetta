import os
import requests

# Create a directory to store the PDB files
pdb_dir = './pdb_files'
os.makedirs(pdb_dir, exist_ok=True)

# List of PDB IDs
pdb_ids = ['5nkq', '5A4D', '5NKX']  # Add the correct PDB IDs here

# Function to fetch and download a PDB file
def fetch_pdb(pdb_id):
    # Construct the URL
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    
    # File path to save the PDB file
    pdb_file_path = os.path.join(pdb_dir, f'{pdb_id}.pdb')
    
    # Fetch and download the PDB file
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Write the PDB file
        with open(pdb_file_path, 'w') as f:
            f.write(response.text)
        print(f'{pdb_id}.pdb downloaded successfully to {pdb_file_path}')
    else:
        print(f'Failed to download {pdb_id}.pdb, status code: {response.status_code}')

# Loop through each PDB ID and download the file
for pdb_id in pdb_ids:
    fetch_pdb(pdb_id)

