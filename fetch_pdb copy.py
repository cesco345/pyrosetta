import requests
from pyrosetta import pose_from_pdb

def fetch_pdb(pdb_id, output_file=None):
    """
    Fetches a PDB file from the RCSB Protein Data Bank (PDB) and saves it locally.
    
    Parameters:
    - pdb_id: The ID of the PDB file to fetch (e.g., '5tj3').
    - output_file: The local file path to save the PDB file. If not specified, saves as '<pdb_id>.pdb'.
    """
    pdb_url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    
    # Fetch the PDB file from RCSB PDB
    response = requests.get(pdb_url)
    
    if response.status_code == 200:
        # Determine the output file name
        if output_file is None:
            output_file = f'{pdb_id}.pdb'
        
        # Write the PDB content to the local file
        with open(output_file, 'w') as file:
            file.write(response.text)
        
        print(f"PDB file '{pdb_id}' successfully saved as '{output_file}'.")
        return output_file
    else:
        print(f"Failed to fetch PDB file for '{pdb_id}'. Status code: {response.status_code}")

# Example usage: Fetch the PDB file for '5tj3' and save it as '5tj3.pdb'
fetch_pdb('1abc')
