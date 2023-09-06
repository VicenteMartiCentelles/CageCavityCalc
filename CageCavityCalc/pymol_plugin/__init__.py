


from __future__ import absolute_import
from __future__ import print_function


#import sys
#sys.path.insert(0, '/u/fd/chem1540/github/CageCavityCalc/')

#from CageCavityCalc import cavity
#UI edited with "Qt Designer"

import os


def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('CageCavCalc', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

#    if dialog is None:  # If the dialos is not created each time, the list of molecules is not updatd
#        dialog = make_dialog()
    dialog = make_dialog() #creating the dialog each time that the plugin is run, the list of new molecules is updated without requiring to close pymol

    dialog.show()


def make_dialog():
    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt

    # create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'cavecav_widget.ui')
    form = loadUi(uifile, dialog)

    
    #form.selection_name.clear()
    form.selection_name.addItems(cmd.get_object_list('all'))

    # callback for the "Ray" button
    def run():
        # get form data
        #height = form.input_height.value()
        #width = form.input_width.value()
        #dpi = form.input_dpi.value()
        #filename = form.input_filename.text()
        #units = form.input_units.currentText()

        # calculate dots per centimeter or inch
        #if units == 'cm':
        #    dots_per_unit = dpi * 2.54
        #else:
        #    dots_per_unit = dpi

        # convert image size to pixels
        #width *= dots_per_unit
        #height *= dots_per_unit

        # render the image
        #if filename:
        #    cmd.png(filename, width, height, dpi=dpi, ray=1, quiet=0)
        #else:
        #    cmd.ray(width, height, quiet=0)
        #    print('No filename selected, only rendering on display')
        selection = form.selection_name.currentText()
        grid_size = float(form.grid_edit.text())
        hydro = form.hydro_check.isChecked()
        aro = form.aro_check.isChecked()
        sas = form.sas_check.isChecked()
        esp = form.esp_check.isChecked()
        method = form.hydro_method.currentText()
        dist = form.hydro_dist.currentText()
        metal_user_charge = form.metalESPcheckBox.isChecked()
        metal_name = form.metal_name_text.text()
        metal_charge = int(form.metal_charge_text.text())
        chargeModel = form.chargeModelComboBox.currentText()

        if form.volume_largest_check.isChecked() == True:
            cluster_text = "size"
        if form.volume_close_to_center_check.isChecked() == True:
            cluster_text = "dist"        
        if form.volume_all_check.isChecked() == True:
            cluster_text = "false"
        
        dist_90_a = form.dist_90_auto.isChecked()
        dist_90_a_v = float(form.dist_90_auto_value.text()) 
        dist_90_m = form.dist_90_manual.isChecked()
        dist_90_m_v = float(form.dist_90_manual_value.text())

        #form.calculate_volume.setText("Calculating... it might take a while (~1 min)")
        
        print(f"Running CageCavCalc with selection={selection}, grid_size={grid_size:}, hydro={hydro:}, aro={aro:}, sas={sas:},esp={esp:}, method={method:}, dist={dist:}, cluster={cluster_text:}, metal_name={metal_name:}, metal_charge={metal_charge:}, metal_user_charge={metal_user_charge:}, chargeModel={chargeModel:}")

        #dialog.close()

        show_cavity_in_pymol(selection, grid_size, hydro, aro, sas, esp, method, dist, cluster_text, dist_90_a, dist_90_a_v, dist_90_m, dist_90_m_v, metal_name, metal_charge, metal_user_charge, chargeModel)




#    def change_slider():
#        slider = int(form.grid_slider.value())+1
#        grid_size = slider*0.2
#        form.grid_edit.setText(f"{grid_size:.1f}")

#    def change_edit():
#        slider = int(float(form.grid_edit.text())/0.2)-1
#        if slider <0:
#           slider=0
#        elif slider>24:
#           slider=24
#        form.grid_slider.setValue(slider)

    # hook up button callbacks
    form.calculate_volume.clicked.connect(run)
    #form.button_browse.clicked.connect(browse_filename)
    #form.button_close.clicked.connect(dialog.close)
    #form.grid_slider.valueChanged.connect(change_slider)
    #form.grid_edit.editingFinished.connect(change_edit)

    return dialog


def show_cavity_in_pymol(selection, grid_size=1, hydro=True, aro=False, sas=False, esp=False, method='Ghose', dist="Fauchere", cluster="false", dist_90_a=True, dist_90_a_v=2.0, dist_90_m=False, dist_90_m_v=5.0, metal_name="Pd", metal_charge=2,metal_user_charge=False, chargeModel="eem"):
    from CageCavityCalc import CageCavityCalc
    from pymol import cmd
    from pymol import stored
    import chempy

    #print(
    #    f"Running CageCavCalc with grid_size={grid_size:}, hydro={hydro:}, aro={aro:}, sas={sas:}, method={method:}, dist={dist:}, cluster={cluster_text:}")
    cav = CageCavityCalc.cavity()
    cav.dummy_atom_radii = grid_size
    cav.grid_spacing = grid_size
    
    cav.clustering_to_remove_cavity_noise = cluster

    # this array will be used to hold the coordinates.  It
    # has access to PyMOL objects and, we have access to it.

    # let's just get the alpha carbons, so make the
    # selection just for them
    # userSelection = userSelection + " and n. CA"

    # iterate over state 1, or the userSelection -- this just means
    # for each item in the selection do what the next parameter says.
    # And, that is to append the (x,y,z) coordinates to the stored.alphaCarbon
    # array.
    stored.pos = []
    stored.names = []

    cmd.iterate_state(1, selection, "stored.pos.append([x,y,z])")  # TODO this should not be all but selected from list
    cmd.iterate_state(1, selection, "stored.names.append(name)")
    cav.read_pos_name_array(stored.pos, stored.names)
    #Set the distance_threshold_for_90_deg_angle as 3 times the window radius
    if dist_90_a == True:
        dist_90_text = "Auto"
        window_radius = cav.calculate_window()
        cav.distance_threshold_for_90_deg_angle = window_radius*dist_90_a_v
        if cav.distance_threshold_for_90_deg_angle < 5:
            cav.distance_threshold_for_90_deg_angle = 5
        print(f"Distance threshold for 90 deg angle = {cav.distance_threshold_for_90_deg_angle:.2f} A")
    if dist_90_m == True:
        dist_90_text = "Manual"
        cav.distance_threshold_for_90_deg_angle = dist_90_m_v
    
    dist_90_text = dist_90_text+str(round(cav.distance_threshold_for_90_deg_angle,1))
    volume = cav.calculate_volume()
    print("Volume of the cavity=", round(volume,2), "A^3")
    
    cavity_text = str(round(volume,2))+"A3_"+str(grid_size)+"A_"+str(selection)+"_"+dist_90_text
    index_to_propery = {0: "aromaticity", 1: "solvent_accessibility", 2: "hydrophobicity", 3: "electrostatics"}
    

    if hydro or aro or sas:
        cavity_hydrophobicity_values = cav.calculate_hydrophobicity()
        average_cavity_hydrophobicity = sum(cavity_hydrophobicity_values)/len(cavity_hydrophobicity_values)
        print(f"Average cavity hydrophobicity = {average_cavity_hydrophobicity:.5f} A^-3")
        print(f"Total cavity hydrophobicity = {average_cavity_hydrophobicity*volume:.5f}")
        mlp_pos = [i for i in cavity_hydrophobicity_values if i > 0]
        mlp_neg = [i for i in cavity_hydrophobicity_values if i < 0]
        lipophilic_index = sum(mlp_pos)/(sum(mlp_pos)-sum(mlp_neg))
        print(f"Lipophilic_index (LI) = {lipophilic_index:.3f}")
    if esp:
        if metal_user_charge:
            cav.calculate_esp(metal_name=metal_name, metal_charge=metal_charge, method=chargeModel) # charges of metals nee to be introduced by the user
        else:
            cav.calculate_esp(method=chargeModel)

    for idx, condition in enumerate([aro, sas, hydro, esp]):
        if condition:
            property_name = index_to_propery[idx]
            property_values = list(cav.get_property_values(property_name))
            model = cmd.get_model(None)

            #print(len(property_values), "==", len(cav.dummy_atoms_positions))

            for property_value, position in zip(property_values[cav.n_atoms:], cav.dummy_atoms_positions):
                # print(property_value, position)
                atom = chempy.Atom()
                atom.name = 'D'
                atom.coord = position
                atom.b = property_value
                model.add_atom(atom)

            cmd.load_model(model, property_name+"_"+str(selection)+"_"+cavity_text)
            cmd.show_as("surface", selection=property_name+"_"+str(selection)+"_"+cavity_text)

            # print(property_name, property_values)
            cmd.spectrum("b", selection=property_name+"_"+str(selection)+"_"+cavity_text, palette="blue_white_red", minimum=min(property_values),
                         maximum=max(property_values))
            cmd.ramp_new("ramp_" + property_name+"_"+str(selection)+"_"+cavity_text, property_name+"_"+str(selection)+"_"+cavity_text,
                         [min(property_values), (min(property_values) + max(property_values)) / 2,
                          max(property_values)], ["blue", "white", "red"])
            cmd.recolor()

    # just volume will be last
    model = cmd.get_model(None)
    for position in cav.dummy_atoms_positions:
        atom = chempy.Atom()
        atom.name = 'D'
        atom.coord = position
        model.add_atom(atom)
    cmd.load_model(model, "cavity_"+cavity_text)
    cmd.alter('name D', 'vdw="'+str(grid_size)+'"')
    cmd.show_as("surface", selection="cavity_"+cavity_text)
    
    cmd.clip("atoms", 5, "All")
    cmd.orient(selection)
    cmd.zoom(selection)
