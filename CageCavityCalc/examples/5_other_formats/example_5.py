
from CageCavityCalc.CageCavityCalc import cavity
cav = cavity()
cav.read_file("cage.gro")
cav.calculate_volume()
cav.print_to_file("cage_cavity.gro")

