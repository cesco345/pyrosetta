# Load structures
load original.pdb, original
load optimized.pdb, optimized

# Set view and style
set_view (\
     0.9641,   -0.2653,   -0.0194,\
     0.2655,    0.9639,    0.0172,\
     0.0150,   -0.0211,    0.9997,\
     0.0000,    0.0000, -197.7769,\
   -10.9400,   -0.2485,   67.1899,\
   155.9365,  239.6174,  -20.0000 )
set cartoon_fancy_helices, 1
set cartoon_transparency, 0.5

# Color schemes
color skyblue, original
color lightpink, optimized

# Show K49 and interacting residues
select k49_orig, original and resi 49
select k49_opt, optimized and resi 49
select interact_orig, byres (k49_orig around 5)
select interact_opt, byres (k49_opt around 5)

# Display interactions
show sticks, k49_orig or interact_orig or k49_opt or interact_opt
show cartoon, k49_orig or interact_orig or k49_opt or interact_opt

# Label residues
label (k49_orig and name CA), '%s%s' % (resn, resi)
label (interact_orig and name CA), '%s%s' % (resn, resi)

# Show non-covalent interactions
distance hbonds_orig, k49_orig, interact_orig, 3.5
distance hbonds_opt, k49_opt, interact_opt, 3.5

# Customize appearance
set label_size, 10
set dash_gap, 0.25
set dash_color, yellow

# Center and zoom
zoom resi 49

# Save session
save k49_comparison.pse