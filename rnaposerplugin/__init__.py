from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

import os
# entry point to PyMOL's API
from pymol import cmd

# pymol.Qt provides the PyQt5 interface, but may support PyQt4
# and/or PySide as well
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi, getSaveFileNameWithExt

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
        # TODO: FILL DIALOG WITH WIDGETS HERE

    dialog.show()


def make_dialog():
    dialog = QtWidgets.QDialog()
    # filename of our UI file
    uifile = os.path.join(os.path.dirname(__file__), 'test.ui')

    # load the UI file into our dialog
    form = loadUi(uifile, dialog)
    print("Dialog created")

    def run():
        # get form data
        receptor = form.receptor_path.text()
        poses = form.poses_path.text()

        # some debugging feedback
        print('User', receptor, poses)

        # TODO: DO SOMETHING WITH FORM DATA
        cmd.load(receptor)
        cmd.load(poses)

    def set_poses_path():
        form.poses_path.setText(QtWidgets.QFileDialog.getOpenFileName()[0])

    def set_receptor_path():
        form.receptor_path.setText(QtWidgets.QFileDialog.getOpenFileName()[0])

    form.browse_receptor.clicked.connect(set_receptor_path)
    form.browse_poses.clicked.connect(set_poses_path)
    form.button_run.clicked.connect(run)

    return dialog
