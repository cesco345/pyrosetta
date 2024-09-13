
# Load structures
load prepared_arac.pdb, initial
load final_arac_structure.pdb, predicted
load ./1XJA.pdb, apo

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
set_view (\
     0.9641,   -0.2655,   -0.0167,\
     0.2569,    0.9475,   -0.1915,\
     0.0673,    0.1775,    0.9814,\
     0.0000,    0.0000,  -83.6658,\
    36.2634,   21.9765,   22.6551,\
   -55.3509,  222.6826,  -20.0000 )

# Save images
ray 1024, 768
png initial_structure.png
disable initial
ray 1024, 768
png predicted_structure.png
disable predicted
enable initial
ray 1024, 768
png apo_structure.png

# Quit PyMOL
quit
