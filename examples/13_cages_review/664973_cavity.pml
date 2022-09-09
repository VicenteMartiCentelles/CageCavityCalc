load 664973_cavity.pdb
extract cavity, resname CV
alter name D, vdw=1.0
show_as surface, cavity
