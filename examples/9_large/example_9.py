from CageCavityCalc import cavity

cav = cavity()
cav.read_file("cage.pdb")



# This works:
cav.distance_threshold_for_90_deg_angle =30
# I've tried to do it automaticlly, that user does not have to chose, but this does not work:
#radius = cav.calculate_window()
#cav.distance_threshold_for_90_deg_angle = radius*1.4+2

volume = cav.calculate_volume()
print(volume)
cav.print_to_file("cage_cavity.pdb")
cav.print_to_pymol("cage_cavity.pml")

