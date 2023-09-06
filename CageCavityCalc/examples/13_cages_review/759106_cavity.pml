load 759106_cavity.pdb
extract cavity, resname CV
alter name D, vdw=1.2
show_as surface, cavity
