from CageCavityCalc import cavity
cav = cavity()
#cav.read_file("cage_Pd2L4.xyz")
cav.read_file("cage_Pd12L24.xyz")
#cav.read_file("C16_1831431_spartan_exported.pdb")
cav.grid_spacing = 1.2

window_radius = cav.calculate_window()
print("Cavity_window radius= ", window_radius, " A")

cav.distance_threshold_for_90_deg_angle = window_radius * 2
if cav.distance_threshold_for_90_deg_angle < 5:
    cav.distance_threshold_for_90_deg_angle = 5

volume = cav.calculate_volume()
cav.print_to_file("cage_cavity.pdb")
cav.print_windows_to_file("cage_windows.pdb")
#cav.print_to_pymol("cage_cavity.pml")
print("Cavity_volume= ", volume, " A3")

windows_ellipsoid = cav.windows_ellipsoid_calculation()
print("Fitted windows ellipsoids: ")
print(windows_ellipsoid)
print("Finished!")
