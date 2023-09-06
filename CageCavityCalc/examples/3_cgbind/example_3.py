from cgbind import Linker, Cage
linker = Linker(smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')
cage = Cage(linker, metal='Pd')

from CageCavityCalc.CageCavityCalc import cavity
cav = cavity()
cav.read_cgbind(cage)
cav.calculate_volume()
cav.print_to_file("cage_cavity.pdb")

