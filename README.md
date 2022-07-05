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
```
You can also calculate hydrophobicity index and visualise it:
```
```

Also cgbind:
```
import cbind
cage =  

cavity_calc = CageCavCalc()
cavity_calc.read_cgbind(cage)
cavity_calc.volume()
```

And lastly from MDAnalysis:
```
import MDAnalysis
syst = MDAnalysis.Universe('cage.gro')
cavity_calc = CageCavCalc()
cavity_calc.read_mdanalysis(syst)
cavity_calc.volume()
```

From bash:
```
CageCavCalc -f cage.mol2 -hydrophobicity -pymol
```