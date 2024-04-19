# CageCavityCalc
Automated calculation of cavity in molecular cages

![Alt text](CageCavityCalc/pic/principle.png "Principle")



## Installation:

The installation of C3 requires the following steps. The software is compatible with Linux, Windows, and Mac.

First it is required to install Miniconda3 (Python 3.7 or later,), that can be obtained from https://docs.conda.io/en/latest/miniconda.html. Then in the “Anaconda Prompt”, the command line version of C3 is installed using pip: “pip install CageCavityCalc”, this will install the required dependencies.

For performing ESP and hydrophobicity calculation, C3 requires OpenBabel that can be installed using “conda install -c conda-forge openbabel” (or “conda config --add channels conda-forge” followed by “conda install openbabel”). For file loading in multiple formats, including files for molecular dynamics cavity analysis, requires installing MDAnalysis using “conda install -c conda-forge mdanalysis” (or “conda config --add channels conda-forge” followed by “conda install mdanalysis”). 

As with any program, to run CageCavityCalc from the command line it is needed to either add its installation folder to the system path or to execute the CageCavityCalc.py file directly from the folder.

To install the C3 PyMol plugin, the open-source version of PyMol must be installed in the “Anaconda Prompt” the command using “conda install -c conda-forge pymol-open-source” for Windows, Mac, and Linux. Alternatively, it can be installed from https://www.cgohlke.com/ for Windows. It is also required to install the following dependencies: “pip install pyqt5 qtpy” (in some cases it may require uninstall pyqt5 with “pip uninstall pyqt5” followed by “pip install pyqt5 qtpy”). For Spanish computers, for running the plugin it is required to change the regional settings of the computer to use points as a decimal separator instead of commas.
Once PyMol is installed, in PyMol the plugin is installed from: Plugin > Plugin Manager > Install New Plugin. Choose “Install from local file” and locate the __init__.py file in the pymol_plugin folder of C3 typically located in C:\Users\UserName\miniconda3\Lib\site-packages\CageCavityCalc\pymol_plugin. To use the plugin, the user just needs to open a cage in PyMol, then go to Plugin and click on CageCavityCalc to run the plugin.  Then, the computed cavity is displayed in PyMol. The user can select the computed property to display by just clicking on the right panel of the generated cavity objects.

Summary of the required installation actions and commands 

Download and install Miniconda3 (https://docs.conda.io/en/latest/miniconda.html)

Open the “Anaconda Prompt” and execute the following commands:

```
pip install CageCavityCalc
conda install -c conda-forge pymol-open-source
pip uninstall pyqt5
pip install pyqt5 qtpy
conda install -c conda-forge openbabel
conda install -c conda-forge mdanalysis
pymol

Then, install the PyMol the plugin: Plugin > Plugin Manager > Install New Plugin.
Choose “Install from local file” and locate the __init__.py file in the pymol_plugin folder of C3
(typically located in C:\Users\UserName\miniconda3\Lib\site-packages\CageCavityCalc\pymol_plugin).
```

### Quick start


The module can be used from the command line or from a python file by loading the CageCavityCalc module. For example, to use C3 from the command line the user needs to executee in the console the following commands: $python CageCavityCalc.py -f cage.pdb -o cage_cavity.pdb -gr 1.5. This order will load the cage.pdb file containing the cage chemical structure and the cavity of the cage will be calculated using a grid spacing of 1.5 Å. Additional arguments can be used as described in Table S1, allowing specifying the distance threshold used to calculate 90º angle, the use of the clustering algorithm to remove noisy cavity points that does not belogin to the main cavity, calculation of hydrophobicity specifying the method and distance function, calculation of hydrophobicity, save a PyMol pml file, or print additional information of the calculations in the terminal.
```
 Arguments that can be used int the C3 Python module though the command line.
-f	Input file (*pdb, *mol2, ...)
-o	Output file (*pdb, *mol2, ...)
-gr X	Grid spacing resolution (Angstroms). Default 1.0
-d90a X	Automatic distance threshold to calculate 90 deg angle as X times window radius. Default 2.0
-d90m X	Manual distance threshold to calculate 90 deg angle in Å
-cluster false, size or dist	Remove cavity noise by dbscan clustering (size or dist)
-hydrophobicity or -hydro	Calculate hydrophobicity
-method Ghose or Crippen	Method to calculate the hydrophobicity
-distfun Audry, Fauchere, Fauchere2, or OnlyValues:	Method to calculate the hydrophobicity
-esp	Calculate the electrostatic potential
-pymol	Create PyMol pml file
-info	Print log INFO on the terminal
```

To use C3 as in a Python script, it is required to load the module and perform. Then, the user needs to specify the name of the file to be loaded, in the example below it is loaded the cage.pdb file, then the cavity is computed using a grid spacing of 1.0 Å and a distance threshold for the 90-degree calculation of 2.0 times the window size. Note that this code uses the same implementation of the distance threshold for the 90-degree of the PyMol plugin, if this is smaller than 5 Å it is set to 5 Å to ensure a correct calculation. The cavity is be saved into .pdb file and also a PyMol *.pml file to facilitate cavity visualization in PyMol. The calculated properties can also be saved, the example below shows how toa calculate the hydrophobicity and save a .pdb file with the hydrophobicity values stored in the B-factor of the .pdb file and as PyMol *.pml file to facilitate cavity visualization with the hydrophobicity in PyMol.

```
cage_name = "cage"
grid_spacing = 1.0
distance_threshold_for_90_deg_angle = 2.0

cav = cavity()
cav.read_file(cage_name+".pdb")
window_radius = cav.calculate_window()
cav.distance_threshold_for_90_deg_angle = window_radius * distance_threshold_for_90_deg_angle
if cav.distance_threshold_for_90_deg_angle < 5:
    cav.distance_threshold_for_90_deg_angle = 5
cav.grid_spacing = float(grid_spacing)
cav.dummy_atom_radii = float(grid_spacing)
volume = cav.calculate_volume()
cav.print_to_file("cage_cavity.pdb")
cav.print_to_pymol("cage_cavity.pml")
cav.calculate_hydrophobicity()
cav.print_to_file("cage_cavity_hydrophobicity.pdb")
cav.print_to_pymol("cage_cavity_hydrophobicity.pml", 'h')

print("Cavity_volume= ", volume, " A3")
```

Example of cavity visulaization in PyMol using teh saved *.pml file.


![Alt text](CageCavityCalc/pic/cavity.png "Principle")

Example of cavity visulaization with hydrophobicity in PyMol using teh saved *.pml file.


![Alt text](CageCavityCalc/pic/hydrophobicity.png "Principle")


The C3 Python module is integrated into PyMol in a plugin. The plugin is integrated into the software through a user interface allowing the selection of the different parameters for the cavity calculation. First, the user needs to initiate PyMol by typing “pymol” in the Anaconda Prompt. Then, in the PyMol interface the user needs to load the desired cage file using File > Open and select the “cage.pdb” file. Then, to initiate the C3 plugin, the user needs access to Plugin > CageCavityCalc. Once all the options are selected, the user needs to click on “Calculate volume” to initiate the calculation of the cavity and all the selected properties. Once the computation is finished, the computed cavity and the cavity with the properties are displayed in PyMol. The PyMol plugin enables the storage of all computed properties in the same PyMol session file, allowing the user to select which one to display and to save PDB files of each property. To save the session file, the user needs to access to File > Save Session As. The user can select the computed property to display by just clicking on the right panel of the generated cavity objects (see Figure 9 in the manuscript). To obtain a good quality image of the cage and the cavity, the user needs to type “ray” in the PyMol command line, then the obtained image can be saved by using File > Export Image As > PNG. 

![Alt text](CageCavityCalc/pic/C3_PyMol_Plugin.png "C3_PyMol_Plugin")

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
