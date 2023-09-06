
from CageCavityCalc.CageCavityCalc import cavity
cav = cavity()

import MDAnalysis
syst = MDAnalysis.Universe("short.gro", "short.xtc")

volume = []
for ts in syst.trajectory:
    cav.read_mdanalysis(syst)
    volume.append(cav.calculate_volume())

print(volume)

