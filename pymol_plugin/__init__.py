'''
PyMOL Demo Plugin

The plugin resembles the old "Rendering Plugin" from Michael Lerner, which
was written with Tkinter instead of PyQt.

(c) Schrodinger, Inc.

License: BSD-2-Clause
'''

from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

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

    if dialog is None:
        dialog = make_dialog()

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
        cluster = form.largest_check.isChecked()
        
        print(f"Running CageCavCalc with selection={selection}, grid_size={grid_size:}, hydro={hydro:}, aro={aro:}, sas={sas:},esp={esp:}, method={method:}, dist={dist:}, cluster={cluster:}")

        dialog.close()


    def change_slider():
        slider = int(form.grid_slider.value())+1
        grid_size = slider*0.2
        form.grid_edit.setText(f"{grid_size:.1f}")

    def change_edit():
        slider = int(float(form.grid_edit.text())/0.2)-1
        if slider <0:
           slider=0
        elif slider>24:
           slider=24
        form.grid_slider.setValue(slider)


    # hook up button callbacks
    form.calculate_volume.clicked.connect(run)
    #form.button_browse.clicked.connect(browse_filename)
    #form.button_close.clicked.connect(dialog.close)
    form.grid_slider.valueChanged.connect(change_slider)
    form.grid_edit.editingFinished.connect(change_edit)

    return dialog
