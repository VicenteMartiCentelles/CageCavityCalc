load 646911_cavity.pdb
extract cavity, resname CV
alter name D, vdw=1.5
show_as surface, cavity
