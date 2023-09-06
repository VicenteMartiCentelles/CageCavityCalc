from CageCavityCalc import cavity
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", default=None, help="input file (*pdb, *mol2, ...)")
    parser.add_argument("-o", default="cage_cavity.pdb", help="output file (*pdb, *mol2, ...)")
    parser.add_argument("-hydrophobicity", "-hydro", default=False, action='store_true',
                        help="Calculate hydrophobicity")
    parser.add_argument("-esp", default=False, action='store_true', help="Calculate ESP")
    parser.add_argument("-method", default="Ghose", help="Method to calculate the hydrophobicity: Ghose or Crippen")
    parser.add_argument("-pymol", default=False, action='store_true', help="create pymol script")
    parser.add_argument("-distfun", default="Fauchere",
                        help="Method to calculate the hydrophobicity: Audry, Fauchere, Fauchere2, OnlyValues")
    parser.add_argument("-gr", default=1.0, help="Grid spacing resolution (Angstroms)")
    parser.add_argument("-info", default=False, action='store_true', help="Print log INFO on the terminal")
    parser.add_argument("-cluster", default="false", help="Remove cavity noise by dbscan clustering (size or dist)")
    parser.add_argument("-d90a", default=3.0,
                        help="Automatic distance threshold to calculate 90 deg angle as X times window radious")
    parser.add_argument("-d90m", default=None, help="Manual distance threshold to calculate 90 deg angle in Angstroms")
    return parser.parse_args()


def main():
    args = get_args()
    if args.info == True:
        logger.setLevel('INFO')
    if args.f is None:
        print("input file (-f) is required")
        exit(0)
    # Set the output filename as the input filename + "_cavity.pdb"
    if args.o == "cage_cavity.pdb":
        args.o = args.f[:args.f.find('.')] + "_cavity.pdb"
    # Initialze the cavity
    cav = cavity()
    # Set the grid spacing resolution from the args input. Default = 1 Angstrom
    if args.gr:
        cav.grid_spacing = float(args.gr)
        cav.dummy_atom_radii = float(args.gr)
    # Set the
    if args.cluster:
        cav.clustering_to_remove_cavity_noise = args.cluster
    # Read the input file
    cav.read_file(args.f)
    # Set the distance_threshold_for_90_deg_angle as 3 times the window radius by default or read from input line
    window_radius = cav.calculate_window()

    cav.distance_threshold_for_90_deg_angle = window_radius * float(args.d90a)
    if args.d90m:
        cav.distance_threshold_for_90_deg_angle = float(args.d90m)

    if cav.distance_threshold_for_90_deg_angle < 5:
        cav.distance_threshold_for_90_deg_angle = 5
    logger.info(f"Distance threshold for 90 deg angle = {cav.distance_threshold_for_90_deg_angle:.2f}")
    volume = cav.calculate_volume()
    print("Cage cavity volume = ", volume, " A3")

    if args.esp == True:
        cav.calculate_esp()  # this is problematic if there is metal
        cav.print_to_file(args.o, 'esp')
        if args.pymol == True:
            pymol_filename = args.o[:args.o.find('.')] + ".pml"
            cav.print_to_pymol(pymol_filename, 'esp')
    elif args.hydrophobicity == True:
        cav.hydrophMethod = args.method
        cav.distance_function = args.distfun  # Audry, Fauchere, Fauchere2, OnlyValues
        cavity_hydrophobicity_values = cav.calculate_hydrophobicity()
        average_cavity_hydrophobicity = np.mean(cavity_hydrophobicity_values)
        print(f"Hydrophobicity method and distance function = {cav.hydrophMethod}, {cav.distance_function}")
        print(f"Average cavity hydrophobicity = {average_cavity_hydrophobicity:.5f} A^-3")
        print(f"Total cavity hydrophobicity = {average_cavity_hydrophobicity * volume:.5f}")
        mlp_pos = [i for i in cavity_hydrophobicity_values if i > 0]
        mlp_neg = [i for i in cavity_hydrophobicity_values if i < 0]
        lipophilic_index = sum(mlp_pos) / (sum(mlp_pos) - sum(mlp_neg))
        print(f"Lipophilic_index (LI) = {lipophilic_index:.3f}")

        cav.print_to_file(args.o, 'h')
        if args.pymol == True:
            pymol_filename = args.o[:args.o.find('.')] + ".pml"
            cav.print_to_pymol(pymol_filename, 'h')
    else:
        cav.print_to_file(args.o)
        if args.pymol == True:
            pymol_filename = args.o[:args.o.find('.')] + ".pml"
            cav.print_to_pymol(pymol_filename)



if __name__ == '__main__':
    main()


