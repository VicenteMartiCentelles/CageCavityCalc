import MDAnalysis
syst = MDAnalysis.Universe("cage.gro")

from CageCavityCalc.CageCavityCalc import cavity
cav = cavity()
cav.read_mdanalysis(syst)
cav.calculate_volume()
cav.calculate_hydrophobicity()
cav.print_to_file("cage_cavity.pdb", 'h')

