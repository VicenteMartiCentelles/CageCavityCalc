from CageCavityCalc import cavity
import argparse
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", default=None, help="input file (*pdb, *mol2, ...)")
    parser.add_argument("-o", default="cage_cavity.pdb", help="output file (*pdb, *mol2, ...)")
    parser.add_argument("-hydrophobicity", "-hydro", default=False, action='store_true',
                        help="Calculate hydrophobicity")
    parser.add_argument("-esp", default=False, action='store_true', help="Calculate ESP")
    parser.add_argument("-metal", default="", help="Metal used in the ESP")
    parser.add_argument("-metal_charge", default=0, help="Charge of the metal used in the ESP")
    parser.add_argument("-charge_method", default="eem", help="Charge method used in the ESP: eem, mmff94, gasteiger, qeq, qtpie, eem2015ha, eem2015hm, eem2015hn, eem2015ba, eem2015bm, eem2015bn")
    parser.add_argument("-method", default="Ghose", help="Method to calculate the hydrophobicity: Ghose or Crippen")
    parser.add_argument("-pymol", default=False, action='store_true', help="create pymol script")
    parser.add_argument("-distfun", default="Fauchere",
                        help="Method to calculate the hydrophobicity: Audry, Fauchere, Fauchere2, OnlyValues")
    parser.add_argument("-gr", default=1.0, help="Grid spacing resolution (Angstroms)")
    parser.add_argument("-info", default=False, action='store_true', help="Print log INFO on the terminal")
    parser.add_argument("-cluster", default="false", help="Remove cavity noise by dbscan clustering (size or dist)")
    parser.add_argument("-d90a", default=2.0,
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
    # Set the cluster
    if args.cluster:
        cav.clustering_to_remove_cavity_noise = args.cluster
    # Read the input file
    cav.read_file(args.f)
    
    window_radius = cav.calculate_window()

    cav.distance_threshold_for_90_deg_angle = window_radius * float(args.d90a)
    if cav.distance_threshold_for_90_deg_angle < 5:
        cav.distance_threshold_for_90_deg_angle = 5

    if args.d90m:
        cav.distance_threshold_for_90_deg_angle = float(args.d90m)

    print(f"Distance threshold for 90 deg angle = {cav.distance_threshold_for_90_deg_angle:.2f}")
    volume = cav.calculate_volume()
    
    #save the PDB file with the computed cavity
    cav.print_to_file(args.o)

    if args.pymol == True:
        cav.print_to_pymol(args.o[:args.o.find('.')] + "_cavity.pml")
    
    print("Cage cavity volume = ", volume, " A3")

    if args.esp == True:
        if (args.metal and args.metal_charge):
            print(f"Metal = {args.metal}, Charge = {args.metal_charge}")
            cav.calculate_esp(metal_name=args.metal, metal_charge=int(args.metal_charge))
            if (args.charge_method):
                cav.calculate_esp(metal_name=args.metal, metal_charge=int(args.metal_charge), method=args.charge_method)
        else:
            cav.calculate_esp()  # this is problematic if there is metal
            if (args.charge_method):
                cav.calculate_esp(method=args.charge_method)
        cav.print_to_file(args.o[:args.o.find('.')] + "_esp.pdb", 'esp')
        if args.pymol == True:
            pymol_filename = args.o[:args.o.find('.')] + "_esp.pml"
            cav.print_to_pymol(pymol_filename, 'esp')
    if args.hydrophobicity == True:
        cav.hydrophMethod = args.method
        cav.distance_function = args.distfun  # Audry, Fauchere, Fauchere2, OnlyValues
        cavity_hydrophobicity_values = cav.calculate_hydrophobicity()
        average_cavity_hydrophobicity = np.mean(cavity_hydrophobicity_values)
        print(f"Hydrophobicity method and distance function = {cav.hydrophMethod}, {cav.distance_function}")
        print(f"Average cavity hydrophobicity = {average_cavity_hydrophobicity:.5f} A^-3")
        print(f"Total cavity hydrophobicity = {average_cavity_hydrophobicity * volume:.5f}")
        mlp_pos = [i for i in cavity_hydrophobicity_values if i > 0]
        mlp_neg = [i for i in cavity_hydrophobicity_values if i < 0]
        hydrophobic_index = sum(mlp_pos) / (sum(mlp_pos) - sum(mlp_neg))
        print(f"Hydrophobic_index (HI) = {hydrophobic_index:.3f}")

        cav.print_to_file(args.o[:args.o.find('.')] + "_hydrophob.pdb", 'h')
        if args.pymol == True:
            pymol_filename = args.o[:args.o.find('.')] + "_hydrophob.pml"
            cav.print_to_pymol(pymol_filename, 'h')



if __name__ == '__main__':
    main()


