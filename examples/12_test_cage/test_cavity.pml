load test_cavity.pdb
extract cavity, resname CV
alter name D, vdw=0.7
show_as surface, cavity
spectrum b, blue_white_red,cavity, minimum=0.000000, maximum=0.781007
ramp_new 'ramp', cavity, [0.000000,0.508777,0.781007], ['blue','white','red']
recolor
