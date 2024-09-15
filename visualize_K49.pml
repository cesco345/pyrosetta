
# Load structures
load original_K49.pdb, original
load all_trans_K49.pdb, all_trans

# Color the structures
color white, original
color cyan, all_trans

# Show as cartoon with stick sidechains
show cartoon
show sticks, resi 49

# Focus on K49
zoom resi 49

# Label the residue
label resi 49 and name CA, resn+resi

# Show clashes
distance clash, all_trans and resi 49, all_trans and not resi 49, 3.0
color red, clash

# Align structures
align all_trans, original

# Set view
set_view (\
     0.9641,   -0.2653,   -0.0194,\
     0.2655,    0.9639,    0.0172,\
     0.0150,   -0.0211,    0.9997,\
     0.0000,    0.0000, -197.7769,\
   -10.9400,   -0.2485,   67.1899,\
   155.9365,  239.6174,  -20.0000 )

# Save session
save K49_comparison.pse
    