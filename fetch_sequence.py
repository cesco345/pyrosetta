import os
import requests

def fetch_pdb_sequence(pdb_code, output_dir='./sequences'):
    """Fetches the sequence of a protein from the RCSB PDB and saves it in both FASTA and text formats."""
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # RCSB PDB sequence URL (FASTA format)
    url = f"https://www.rcsb.org/fasta/entry/{pdb_code}"

    # Try to fetch the sequence
    try:
        response = requests.get(url)
        if response.status_code == 200:
            # Save as FASTA file
            fasta_file = os.path.join(output_dir, f"{pdb_code}.fasta")
            with open(fasta_file, 'w') as f:
                f.write(response.text)
            print(f"Sequence saved to {fasta_file}")

            # Also save as plain text
            text_file = os.path.join(output_dir, f"{pdb_code}_sequence.txt")
            with open(text_file, 'w') as f:
                f.write(response.text)
            print(f"Sequence saved to {text_file}")
        else:
            print(f"Failed to fetch sequence for {pdb_code}. HTTP status code: {response.status_code}")
    
    except Exception as e:
        print(f"Error fetching sequence for {pdb_code}: {e}")

# Fetch and save the sequence for 5NKQ
fetch_pdb_sequence("5NKQ")

