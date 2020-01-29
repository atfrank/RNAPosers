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
    print("[RNAPosers Debugging] Dialog created")


    def eta_converter(eta_input):
        if eta_input == '2A, 4A and 8A':
            return '248'
        elif eta_input == '2A and 4A':
            return '24'
        else:
            return '2'

    def run():
        import os
        # get form data and initialize parameters
        rmsd = form.rmsd.currentText()
        eta = eta_converter(form.eta.currentText())
        pdb = form.pdb_filename.text()
        dcd = form.dcd_filename.text()
        start_frame = 1
        try:
            stop_frame = int(form.stop_frame.text())
        except:
            stop_frame = -1
        mol2 = form.mol2_filename.text()
        score = form.output_filename.text()
        complex_name = "complex"

        # some debugging feedback
        print('[RNAPosers Debugging] Parameters:', rmsd, eta, pdb, dcd, stop_frame, score)

        cmd.delete(complex_name)
        cmd.load(pdb, complex_name)
        cmd.load_traj(dcd, complex_name, state=1, stop=stop_frame)
        rnaposers_cmd = " ".join(["./run.sh", pdb, mol2, dcd, rmsd, eta, "test/feature", score, str(stop_frame)])
        os.system(rnaposers_cmd)
        from reorder_traj import reorder_traj
        reorder_traj(complex_name, score)

    def set_saveas_path():
        filename = QtWidgets.QFileDialog.getSaveFileName()[0]
        if filename:
            form.output_filename.setText(filename)


    def make_set_path(form, name):
        def set_path():
            from pymol.Qt import QtWidgets
            set_command = "".join(["form.", name, "_filename.setText(QtWidgets.QFileDialog.getOpenFileName()[0])"])
            eval(set_command)
        return set_path

    def set_dcd_path():
        form.dcd_path.setText(QtWidgets.QFileDialog.getOpenFileName()[0])

    for object in ["pdb", "dcd", "mol2"]:
        set_command = "".join(["form.button_", object, ".clicked.connect(make_set_path(form, '", object, "'))"])
        eval(set_command)

    form.button_saveas.clicked.connect(set_saveas_path)
    form.button_run.clicked.connect(run)
    form.button_close.clicked.connect(dialog.close)

    return dialog
