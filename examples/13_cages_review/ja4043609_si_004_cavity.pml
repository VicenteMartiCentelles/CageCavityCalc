load ja4043609_si_004_cavity.pdb
extract cavity, resname CV
alter name D, vdw=3.0
show_as surface, cavity
spectrum b, blue_white_red,cavity, minimum=-0.048650, maximum=0.107869
ramp_new 'ramp', cavity, [-0.048650,0.018682,0.107869], ['blue','white','red']
recolor
