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
    featureFile = sys.argv[3] # feature file
    start_frame = 1
    try:
        stop_frame = int(sys.argv[4])
    except:
        stop_frame = -1

    complex_name = "complex"

    # some debugging feedback
    print('[RNAPosers Debugging] Parameters:', receptor, poses, stop_frame)

    # redefine pdb and dcd
    with tempfile.TemporaryDirectory() as tmpDir:
        mol2 = tmpDir + "/lig.mol2"
        pdb = tmpDir + "/complex.pdb"
        dcd = tmpDir + "/complexes.dcd"
        generate_complexes(receptor, poses, dcd, pdb, mol2)
        cmd.delete(complex_name)
        cmd.load(pdb, complex_name)
        cmd.load_traj(dcd, complex_name, state=1, stop=stop_frame)
        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = os.environ.get('RNAPOSERS_PATH')
        rnaposers_cmd = " ".join(["bash", dir_path + "/ssMD/rna_poser_ssMD.sh", pdb, mol2, dcd, featureFile, str(stop_frame)])
        print('[RNAPosers Debugging]',rnaposers_cmd)
        os.system(rnaposers_cmd)
    return 0

if __name__ == '__main__':
    main()
