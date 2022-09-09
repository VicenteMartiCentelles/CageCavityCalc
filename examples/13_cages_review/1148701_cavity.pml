load 1148701_cavity.pdb
extract cavity, resname CV
alter name D, vdw=0.5
show_as surface, cavity
