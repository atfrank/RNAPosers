from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('RNAPosers Plugin', run_plugin_gui)

# global reference to avoid garbage collection of our dialog
dialog = None

def run_plugin_gui():

    global dialog

    if dialog is None:
        # create a new (empty) Window
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    import os
    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi, getSaveFileNameWithExt
    # from pymol.Qt.QtWidgets.QFileDialog import getOpenFileName, getSaveFileName

    dialog = QtWidgets.QDialog()
    # filename of our UI file
    uifile = os.path.join(os.path.dirname(__file__), 'rnaposers.ui')

    # load the UI file into our dialog
    global form
    form = loadUi(uifile, dialog)
    print("Dialog created")


    def eta_converter(eta_input):
        if eta_input == '2A, 4A and 8A':
            return '248'
        elif eta_input == '2A and 4A':
            return '24'
        else:
            return '2'

    def run():
        import os
        # get form data
        # get form data
        rmsd = form.rmsd.currentText()
        eta = eta_converter(form.eta.currentText())
        pdb = form.pdb_filename.text()
        dcd = form.dcd_filename.text()
        mol2 = form.mol2_filename.text()
        score = form.output_filename.text()

        # some debugging feedback
        print('User', rmsd, eta, pdb, dcd, score)

        # TODO: DO SOMETHING WITH FORM DATA
        cmd.load(pdb)
        cmd.load_traj(dcd, stop=10)
        rnaposers_cmd = " ".join(["./run.sh", pdb, mol2, dcd, rmsd, eta, "test/feature", score])
        os.system(rnaposers_cmd)

    def set_saveas_path():
        filename = QtWidgets.QFileDialog.getSaveFileName()[0]
        if filename:
            form.output_filename.setText(filename)


    def make_set_path(form, name):
        def set_path():
            from pymol.Qt import QtWidgets
            set_command = "".join(["form.", name, "_filename.setText(QtWidgets.QFileDialog.getOpenFileName()[0])"])
            print(set_command)
            eval(set_command)
        return set_path

    def set_dcd_path():
        form.dcd_path.setText(QtWidgets.QFileDialog.getOpenFileName()[0])

    for object in ["pdb", "dcd", "mol2"]:
        set_command = "".join(["form.button_", object, ".clicked.connect(make_set_path(form, '", object, "'))"])
        print(set_command)
        eval(set_command)

    form.button_saveas.clicked.connect(set_saveas_path)
    form.button_run.clicked.connect(run)
    form.button_close.clicked.connect(dialog.close)

    return dialog
