load cage_cavity_hydrophobicity.pdb
extract cavity, resname CV
alter name D, vdw=1
show_as surface, cavity
spectrum b, blue_white_red,cavity, minimum=0.004322, maximum=0.122192
ramp_new 'ramp', cavity, [0.004322,0.040764,0.122192], ['blue','white','red']
recolor
