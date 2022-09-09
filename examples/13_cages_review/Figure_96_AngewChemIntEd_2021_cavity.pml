load Figure_96_AngewChemIntEd_2021_cavity.pdb
extract cavity, resname CV
alter name D, vdw=2.0
show_as surface, cavity
