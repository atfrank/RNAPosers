import os
import sys
import tempfile
import argparse
from generate_complex_files import generate_complexes
from pymol import cmd

def main():
    # get form data and initialize parameters
    receptor = sys.argv[1] # receptor
    poses = sys.argv[2] # ligand poses
    score = sys.argv[3] # output file path
    rmsd = sys.argv[4]
    eta = sys.argv[5]
    start_frame = 1
    try:
        stop_frame = int(sys.argv[6])
    except:
        stop_frame = -1

    complex_name = "complex"

    # some debugging feedback
    print('[RNAPosers Debugging] Parameters:', receptor, poses, score, rmsd, eta, stop_frame)

    # redefine pdb and dcd

    with tempfile.TemporaryDirectory() as tmpDir:
        mol2 = tmpDir + "/lig.mol2"
        pdb = tmpDir + "/complex.pdb"
        dcd = tmpDir + "/complexes.dcd"
        generate_complexes(receptor, poses, dcd, pdb, mol2)
        cmd.delete(complex_name)
        cmd.load(pdb, complex_name)
        cmd.load_traj(dcd, complex_name, state=1, stop=stop_frame)
        featureFile = tmpDir + "/features"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = os.environ.get('RNAPOSERS_PATH')
        rnaposers_cmd = " ".join(["bash", dir_path + "/src/rna_poser.sh", pdb, mol2, dcd, rmsd, eta, featureFile, score, str(stop_frame)])
        print('[RNAPosers Debugging]',rnaposers_cmd)
        os.system(rnaposers_cmd)
    return 0

if __name__ == '__main__':
    main()
