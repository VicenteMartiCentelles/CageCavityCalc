from pymol import stored


import sys
sys.path.insert(0, '/u/fd/chem1540/github/CageCavityCalc/')

from CageCavityCalc import CageCavityCalc


#from CageCavityCalc import cavity


def simple():

    cav = CageCavityCalc.cavity()

    # this array will be used to hold the coordinates.  It
    # has access to PyMOL objects and, we have access to it.
    

    # let's just get the alpha carbons, so make the
    # selection just for them
    #userSelection = userSelection + " and n. CA"

    # iterate over state 1, or the userSelection -- this just means
    # for each item in the selection do what the next parameter says.
    # And, that is to append the (x,y,z) coordinates to the stored.alphaCarbon
    # array.
    stored.pos = []
    stored.names = []
    print("welcome")
    cmd.iterate_state(1, "all", "stored.pos.append([x,y,z])")
    cmd.iterate_state(1, "all", "stored.names.append(name)")
    cav.read_pos_name_array(stored.pos, stored.names)
    volume = cav.calculate_volume()
    print("Volume of the cavity", volume)


cmd.extend("xxx", simple)
