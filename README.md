# CageCavityCalc
Automated calculation of cavity in molecular cages

![Alt text](CageCavityCalc/pic/principle.png "Principle")



## Installation:
```
conda install rdkit --channel conda-forge
pip install CageCavityCalc
```

### Quick start

Using as a python module:
```
cavity_calc = CageCavCalc()
cavity_calc.read_file("cage.pdb")
cavity_calc.volume()
```
The cavity can be saved into pdb file (and also create pymol *.pml to facilitate cavity visualization) 
```
cavity_calc.print_to_file("cage_cavity.pdb")
cavity_calc.print_to_pymol("cage_cavity.pml")
```
![Alt text](CageCavityCalc/pic/cavity.png "Principle")

You can also calculate hydrophobicity index and visualise it on the cavity:
```
cavity_calc.calculate_hydrophobicity()
cavity_calc.print_to_pymol("cage_cavity_hydrophobicity.pml")
```
![Alt text](CageCavityCalc/pic/hydrophobicity.png "Principle")


From bash:
```
python -m CageCavityCalc -f cage.pdb -o cage_cavity.pdb
```

### Additional support

Reading cage class from cgbind:
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

Reading MDAnalysis universe:
 
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

### Other
To make the calculation loud use CAV_LOG_LEVEL environmental variable

```commandline
export CAV_LOG_LEVEL=INFO
```
