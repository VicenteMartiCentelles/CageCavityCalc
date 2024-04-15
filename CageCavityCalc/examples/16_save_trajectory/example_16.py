# This is the same example 6 but saves trajectory
from CageCavityCalc.CageCavityCalc import cavity
import MDAnalysis


syst = MDAnalysis.Universe("short.gro", "short.xtc")
volume = []
max_grid = 0

for idx, ts in enumerate(syst.trajectory):
    cav = cavity()
    cav.read_mdanalysis(syst.atoms)
    volume.append(cav.calculate_volume())
    cav.print_to_file(f"cage_cavity_{idx:}.pdb")
    if len(cav.dummy_atoms_positions) > max_grid:
        max_grid = len(cav.dummy_atoms_positions)

atom_max = max_grid + cav.n_atoms

# Save as a trajectory:
FileTraj = MDAnalysis.Writer("traj.xtc")
for idx, ts in enumerate(syst.trajectory):
    syst = MDAnalysis.Universe(f"cage_cavity_{idx:}.pdb")
    n_missing = atom_max - len(syst.atoms)

    if n_missing != 0:  # save at the same position as existing grid
        positions = np.array([syst.atoms.select_atoms("name D and resname CV")[0].position] * n_missing)
    else:  # if no cavity, save it at the 0,0,0
        positions = np.zeros((n_missing, 3))

    if n_missing > 0:
        sol = MDAnalysis.Universe.empty(n_missing, trajectory=True)
        sol.add_TopologyAttr('name', ['D'] * n_missing)
        sol.add_TopologyAttr('type', ['D'] * n_missing)
        sol.add_TopologyAttr('resname', ['CV'])
        sol.atoms.positions = positions
        FileTraj.write(MDAnalysis.Merge(syst.atoms, sol.atoms).atoms)
    else:
        FileTraj.write(syst.atoms)

FileTraj.close()
