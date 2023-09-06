import MDAnalysis
syst = MDAnalysis.Universe("cage.gro")

from CageCavityCalc.CageCavityCalc import cavity
cav = cavity()
cav.read_mdanalysis(syst)
volume = cav.calculate_volume()
print("Volume=", volume)
cav.print_to_file("cage_cavity.gro")

