load Figure_73_ChemEurJ_2020_26_1558_cavity.pdb
extract cavity, resname CV
alter name D, vdw=1.2
show_as surface, cavity
