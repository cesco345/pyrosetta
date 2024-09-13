# Cell 1: Import necessary modules and initialize PyRosetta
import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import SmallMover, ShearMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.moves import MonteCarlo
import time
from datetime import timedelta

# Initialize PyRosetta
pyrosetta.init()

# Cell 2: Load the pose and create score function
start_time = time.time()

# Load the pose
pose = pose_from_pdb("./pdb_files/5nkq.pdb")
print(f"Pose loaded. Time taken: {timedelta(seconds=int(time.time() - start_time))}")

# Create score function
scorefxn = get_score_function()
print(f"Score function created. Total time: {timedelta(seconds=int(time.time() - start_time))}")

# Cell 3: Create movers
small_mover = SmallMover()
shear_mover = ShearMover()
min_mover = MinMover()
min_mover.score_function(scorefxn)
mc = MonteCarlo(pose, scorefxn, 1.0)
print(f"Movers created. Total time: {timedelta(seconds=int(time.time() - start_time))}")

# Cell 4: Define the folding algorithm function
def folding_algorithm(pose, small_mover, shear_mover, min_mover, mc, magnitude, start_time):
    small_mover.angle_max('H', magnitude)
    small_mover.angle_max('E', magnitude)
    small_mover.angle_max('L', magnitude)
    shear_mover.angle_max('H', magnitude)
    shear_mover.angle_max('E', magnitude)
    shear_mover.angle_max('L', magnitude)
    
    for i in range(100):
        iter_start_time = time.time()
        
        # Five small moves
        for _ in range(5):
            small_mover.apply(pose)
        
        # Minimize
        min_mover.apply(pose)
        
        # Five shear moves
        for _ in range(5):
            shear_mover.apply(pose)
        
        # Minimize
        min_mover.apply(pose)
        
        # Monte Carlo Boltzmann step
        mc.boltzmann(pose)
        
        iter_time = time.time() - iter_start_time
        total_time = time.time() - start_time
        print(f"Iteration {i+1}/100 complete. Iteration time: {timedelta(seconds=int(iter_time))}, Total time: {timedelta(seconds=int(total_time))}")
    
    return pose

# Cell 5: Run the main folding algorithm
magnitudes = [25, 20, 15, 10, 5]
for idx, magnitude in enumerate(magnitudes):
    magnitude_start_time = time.time()
    print(f"Running folding algorithm with magnitude {magnitude} ({idx+1}/5)")
    pose = folding_algorithm(pose, small_mover, shear_mover, min_mover, mc, magnitude, start_time)
    magnitude_time = time.time() - magnitude_start_time
    total_time = time.time() - start_time
    print(f"Magnitude {magnitude} complete. Time taken: {timedelta(seconds=int(magnitude_time))}, Total time: {timedelta(seconds=int(total_time))}")
    print(f"Score after magnitude {magnitude}: {scorefxn(pose)}")

# Cell 6: Save the final pose
pose.dump_pdb("./pdb_files/folded_5nkq.pdb")
total_time = time.time() - start_time
print(f"Folding complete. Final structure saved as 'folded_5nkq.pdb'. Total time: {timedelta(seconds=int(total_time))}")

# Cell 7: Visualize the initial and final structures (optional)
# Note: This cell requires py3Dmol to be installed
import py3Dmol
from IPython.display import display

def visualize_structures(initial_pdb, final_pdb):
    view = py3Dmol.view(width=800, height=400)
    
    view.addModel(open(initial_pdb, 'r').read(), 'pdb')
    view.setStyle({'model': 0}, {'cartoon': {'color': 'lightgray'}})
    
    view.addModel(open(final_pdb, 'r').read(), 'pdb')
    view.setStyle({'model': 1}, {'cartoon': {'color': 'skyblue'}})
    
    view.zoomTo()
    view.render()
    display(view)
    
    print("Visualization Guide:")
    print("- Light gray: Initial structure")
    print("- Sky blue: Final folded structure")

visualize_structures("./pdb_files/5nkq.pdb", "./pdb_files/folded_5nkq.pdb")
