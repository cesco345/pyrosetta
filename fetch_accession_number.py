import requests

def fetch_uniprot_accession(pdb_id):
  """Fetches UniProt accession number for a given PDB ID.

  Args:
      pdb_id (str): The PDB ID to query.

  Returns:
      str or None: The UniProt accession number if found, otherwise None.
  """

  url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'
  response = requests.get(url)
  response.raise_for_status() 

  data = response.json()

  # Extract UniProt accession from entity_poly
  for entity in data.get('entity_poly', []):
    for chain in entity.get('pdbx_strand_id', '').split(','):
      for ref in entity.get('entity_poly_seq_one_letter_code_can', {}).get(chain, {}).get('reference_database_identifiers', []):
        if ref.get('database_name') == 'UniProt':
          return ref.get('database_accession')

  return None  # Return None if no UniProt accession is found

if __name__ == '__main__':
  pdb_id = '5a4d'
  uniprot_accession = fetch_uniprot_accession(pdb_id)

  if uniprot_accession:
    print(f'UniProt accession for PDB ID {pdb_id}: {uniprot_accession}')
  else:
    print(f'No UniProt accession found for PDB ID {pdb_id}')
