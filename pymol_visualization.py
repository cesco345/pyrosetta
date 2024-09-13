import pyrosetta
import numpy as np
import requests
import os

from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.simple_moves import SmallMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.protocols.relax import FastRelax

# Function to fetch PDB info
def fetch_pdb_info(search_term, apo=False):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": search_term
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True
        }
    }
    
    response = requests.post(url, json=query)
    response.raise_for_status()
    data = response.json()
    
    if data['total_count'] > 0:
        if apo:
            apo_structures = [entry['identifier'] for entry in data['result_set'] 
                              if not any('has_ligands' in group and group['has_ligands'] == 'Y' 
                                         for group in entry.get('group_list', []))]
            return apo_structures[0] if apo_structures else None
        else:
            return data['result_set'][0]['identifier']
    return None

# Function to download PDB file
def download_pdb(pdb_id, output_dir='.'):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    response.raise_for_status()
    
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    with open(output_path, 'wb') as f:
        f.write(response.content)
    
    return output_path

# Fetch and download PDB
search_term = "AraC"
pdb_id = fetch_pdb_info(search_term)
if pdb_id:
    pdb_path = download_pdb(pdb_id)
    print(f"Downloaded PDB file: {pdb_path}")
else:
    raise ValueError("PDB not found")

# Initialize PyRosetta
pyrosetta.init()

# Load and prepare the PDB
def prepare_pdb(pdb_file, output_file):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    filtered_lines = [line for line in lines if line.startswith('ATOM') and line[21] == 'A']
    
    with open(output_file, 'w') as f:
        f.writelines(filtered_lines)

    return output_file

prepared_pdb = prepare_pdb(pdb_path, 'prepared_arac.pdb')

# Load the prepared PDB into PyRosetta
pose = pyrosetta.pose_from_pdb(prepared_pdb)

# Set up MoveMap for N-terminal residues
movemap = MoveMap()
for i in range(1, 16):  # First 15 residues
    movemap.set_bb(i, True)
    movemap.set_chi(i, True)

# Set up scoring function and movers
scorefxn = pyrosetta.get_fa_scorefxn()
small_mover = SmallMover(movemap, 1.0, 1)
min_mover = MinMover()
min_mover.movemap(movemap)
min_mover.score_function(scorefxn)

# Set up FastRelax
relax = FastRelax()
relax.set_scorefxn(scorefxn)
relax.set_movemap(movemap)

# Perform ab initio search
def perform_abinito_search(pose, num_cycles=10):
    for _ in range(num_cycles):
        small_mover.apply(pose)
        min_mover.apply(pose)
    relax.apply(pose)
    return pose

final_pose = perform_abinito_search(pose.clone())

# Fetch and download apo form PDB
apo_pdb_id = fetch_pdb_info("AraC", apo=True)
if apo_pdb_id:
    apo_pdb_path = download_pdb(apo_pdb_id)
    print(f"Downloaded apo PDB file: {apo_pdb_path}")
else:
    raise ValueError("Apo PDB not found")

# Save the final structure
final_pose.dump_pdb('final_arac_structure.pdb')

# Generate PyMOL script
def generate_pymol_script(initial_pdb, predicted_pdb, apo_pdb):
    script = f"""
# Load structures
load {initial_pdb}, initial
load {predicted_pdb}, predicted
load {apo_pdb}, apo

# Set different colors for each structure
color cyan, initial
color magenta, predicted
color yellow, apo

# Show structures as cartoons
show cartoon

# Align structures
align predicted, initial
align apo, initial

# Focus on the N-terminal region (first 15 residues)
select nterm, resi 1-15
zoom nterm

# Label structures
label name CA and initial and resi 1, "Initial"
label name CA and predicted and resi 1, "Predicted"
label name CA and apo and resi 1, "Apo"

# Set nice view
set_view (\\
     0.9641,   -0.2655,   -0.0167,\\
     0.2569,    0.9475,   -0.1915,\\
     0.0673,    0.1775,    0.9814,\\
     0.0000,    0.0000,  -83.6658,\\
    36.2634,   21.9765,   22.6551,\\
   -55.3509,  222.6826,  -20.0000 )

# Save session
save arac_comparison.pse
"""
    with open("visualize_arac.pml", "w") as f:
        f.write(script)
    print("PyMOL script generated: visualize_arac.pml")

# Generate the PyMOL script
generate_pymol_script('prepared_arac.pdb', 'final_arac_structure.pdb', apo_pdb_path)

print("Analysis complete. You can now open PyMOL and run the script 'visualize_arac.pml' to view the structures.")
