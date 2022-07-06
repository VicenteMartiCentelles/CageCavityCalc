# CageCavityCalc
Automated calculation of cavity in molecular cages

Several option of working exist:

From python:
```
cavity_calc = CageCavCalc()
cavity_calc.read_file("cage.pdb")
cavity_calc.volume()
```
It also alows to visuale it by creating .pse/pdb file with dummy atoms:
```
#TODO
```
You can also calculate hydrophobicity index and visualise it:
```
#TODO
```

Also cgbind:
```
from cgbind import Linker, Cage
linker = Linker(smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')
cage = Cage(linker, metal='Pd')

from main import cavity
cav = cavity()
cav.read_cgbind(cage)
cav.calculate_volume()
cav.print_to_file("cage_cavity.pdb")
```

And lastly from MDAnalysis:
```
import MDAnalysis
syst = MDAnalysis.Universe("cage.gro")

from main import cavity
cav = cavity()
cav.read_mdanalysis(syst)
cav.calculate_volume()
cav.print_to_file("cage_cavity.pdb")
```

By MDAnalysis it is possible to calculate trajectory:
```commandline
from main import cavity
cav = cavity()

import MDAnalysis
syst = MDAnalysis.Universe("short.gro", "short.xtc")

volume = []
for ts in syst.trajectory:
    cav.read_mdanalysis(syst)
    volume.append(cav.calculate_volume())

print(volume)
```

From bash:
```
python ../../main.py -f cage.pdb -o cage_cavity.pdb
```