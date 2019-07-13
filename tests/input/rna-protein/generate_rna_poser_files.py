import sys
import inspect
from glob import glob
from pymol import cmd
import numpy
import os
import subprocess
import psico.fullinit
import pymol.plugins
pymol.plugins.autoload = {}
pymol.plugins.preferences = {'instantsave': True, 'verbose': False}
pymol.plugins.set_startup_path( [u'~/local_software/share/pymol/plugins/'] , False)

def pred_rna_poser_3dRPC(pdbid = "1B7F", nposes = 200, renumber = False, home="/home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/"):
    # Initialize 
    cmd.reinitialize()
    cmd.bg_color( "white" )
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("antialias", 2)
    cmd.set('ray_shadows','on')
    cmd.set("ray_trace_mode", 1)
    cmd.set("opaque_background", "off")
    cmd.set("sphere_scale", 0.8)
    cmd.set("transparency", 0.3)
    cmd.set( "ray_trace_gain", 0)
    cmd.cd(home)
    # load poses
    for pose in range(1, nposes+1):
        try:
            cmd.load("%s/3drpc-complex%s.pdb"%(pdbid, pose), pdbid)
        except:
            continue

    cmd.select("rna", "resn U+A+C+G+RU+RA+RC+RG+ADE+URA+CYT+GUA")
    cmd.select("protein", "not rna")
    #cmd.intra_fit("%s"%pdbid)
    
    # rename RNA residues
    cmd.alter("resn C+RC", "resn = 'CYT'")
    cmd.alter("resn A+RA", "resn = 'ADE'")
    cmd.alter("resn U+RU", "resn = 'URA'")
    cmd.alter("resn G+RG", "resn = 'GUA'")
    cmd.alter("protein", "chain = 'A'")
    cmd.alter("rna", "chain = 'B'")
    cmd.reset()
    # rename and renumber protein
    cmd.alter("protein", "resn='UNK'")
    cmd.save (filename = "%s_complex.pdb"%pdbid, selection = pdbid, state = 1)
    if renumber: cmd.alter("protein", "resi = 999")
    # save new reference pdb file and corresponding trajectory
    cmd.save (filename = "%s_complex.pdb"%pdbid, selection = pdbid, state = 1)
    psico.exporting.save_traj(filename = "%s_complex.dcd"%pdbid, selection = pdbid)
    # extract protein and save pdb and sd
    cmd.extract("pro", "protein")
    cmd.show("sticks", "pro")
    cmd.center(pdbid)
    cmd.save (filename = "%s_protein.sd"%pdbid, selection = "pro", state = 1)
    cmd.save (filename = "%s_protein.pdb"%pdbid, selection = "pro", state = 1)
    # generate mol2 files
    cmdline1 = "/export/apps/CentOS7/babel/2.4.1/bin/babel -isd %s_protein.sd -omol2 %s_protein_from_sd.mol2"%(pdbid, pdbid)
    cmdline2 = "/export/apps/CentOS7/babel/2.4.1/bin/babel -ipdb %s_protein.pdb -omol2 %s_protein_from_pdb.mol2"%(pdbid, pdbid)
    p = subprocess.Popen(cmdline1,stdin=None,stdout=None, shell=True)
    p.wait()
    p = subprocess.Popen(cmdline2,stdin=None,stdout=None, shell=True)
    p.wait()
    print("Done with: %s: frame: %s"%(pdbid, cmd.count_frames("(%s)"%pdbid)))

# dataset
# /export/apps/pymol-1.8.6/pymol/pymol -cr generate_rna_poser_files.py
rnas="1B7F 1DFU 1JBS 1P6V 1WPU 1WSU 2ASB 2BH2 2QUX 3BX2".split()
for rna in rnas:
    pred_rna_poser_3dRPC(pdbid = rna, nposes = 100, renumber = False, home="/home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/")

