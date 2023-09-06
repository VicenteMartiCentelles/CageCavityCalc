from CageCavityCalc.CageCavityCalc import cavity

cav = cavity()
'''
#Uppsss... forgot to read the structure!

volume = cav.calculate_volume()
print(volume)
cav.print_to_file("cage_cavity.pdb")
'''


#Uppsss... forgot to calculate cavity
cav.read_file("cage.gro")
#volume = cav.calculate_volume()
cav.calculate_hydrophobicity()




