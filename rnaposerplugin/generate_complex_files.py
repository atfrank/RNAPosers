import sys
import inspect
from glob import glob
from pymol import cmd
import psico.fullinit
import os

def generate_complexes(dir, receptor="receptor.mol2", poses="poses.sd", complex = "complex.dcd", ref_complex = "complex.pdb", ref_ligand_mol2 = "reference_ligand.mol2"):
    def fix_names(obj, is_ligand = False):
        # utility functions to fix names
        if is_ligand:
            cmd.alter("%s"%(obj), "resn = 'UNK'")
            cmd.alter("%s"%(obj), "chain = 'Z'")
            cmd.alter("%s"%(obj), "resi = '1'")
        else:
            cmd.alter("%s"%(obj), "type = 'ATOM'")
            cmd.alter("resn rC+C+RC and %s"%(obj), "resn = 'CYT'")
            cmd.alter("resn rA+A+RA and %s"%(obj), "resn = 'ADE'")
            cmd.alter("resn rU+U+RU and %s"%(obj), "resn = 'URA'")
            cmd.alter("resn rG+G+RG and %s"%(obj), "resn = 'GUA'")

    def split_decoys(poses, dir):
        # utility function to split poses
        for a in range(1, 1+cmd.count_states("(%s)"%poses)):
            cmd.frame(a)
            cmd.save("%s/poses/poses_%i.pdb"%(dir, a))
    
    # utility function to generate complexes
    os.system("rm -rf %s/complex"%(dir))
    os.system("rm -rf %s/poses"%(dir))    
    os.system("mkdir -p %s/complex"%(dir))
    os.system("mkdir -p %s/poses"%(dir))
    
    # split posese
    cmd.delete("*")
    cmd.load(poses, "poses")
    fix_names("poses", is_ligand = True)
    split_decoys("poses", dir)
    poses=range(1, 1+cmd.count_states("(poses)"))
    cmd.delete("poses")
    
    # combine poses + receptor
    for a in poses:
        cmd.delete("*")
        cmd.load(receptor)
        cmd.load("%s/poses/poses_%i.pdb"%(dir, a))
        cmd.create("test", "all")
        fix_names("test", is_ligand = False)
        cmd.save("%s/complex/complex_%i.pdb"%(dir, a),"test")
    
    # save dcd
    cmd.delete("*")
    for a in poses:
        cmd.load("%s/complex/complex_%i.pdb"%(dir, a), "complex")
    psico.exporting.save_traj(filename = "%s/%s"%(dir, complex), selection = "complex")

    # generate reference files
    os.system("cp %s/complex/complex_1.pdb %s/%s"%(dir, dir, ref_complex))
    
    # XX
    cmd.delete("*")
    cmd.load("%s/complex/complex_1.pdb"%(dir), "complex")
    cmd.extract("lig", "resn UNK")
    cmd.save("%s/complex/lig.pdb"%(dir), "lig")
    
    os.system("obabel -ipdb %s/complex/lig.pdb  -omol2 -O %s/%s"%(dir, dir, ref_ligand_mol2))
    
    