from CageCavityCalc.CageCavityCalc import cavity

cav = cavity()
cav.read_file("large.pdb")
cav.distance_threshold_for_90_deg_angle =30
cav.grid_spacing = 4

cav.calculate_volume()
cav.calculate_esp(metal_name="Pd", metal_charge=2)
#cav.calculate_esp()
cav.print_to_pymol("cage2.pml", "esp")

