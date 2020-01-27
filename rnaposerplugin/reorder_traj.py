
from pymol import cmd
import pandas as pd

def sele_exists(sele):
    sess = cmd.get_session()
    for name in sess["names"]:
        if name:
            if sele == name[0]:
                return 1
    return 0
cmd.extend('sele_exists', sele_exists)

def combine_frame(dcd_new, dcd_original, frame=1):
    if not sele_exists(dcd_new):
        cmd.create(dcd_new, dcd_original, frame, 1)
    else:
        cmd.create(dcd_new, dcd_original, frame, -1)
cmd.extend('combine_frame', combine_frame)

def get_order_frame(score_file):
    score_table = pd.read_csv(score_file, sep=' ', header = None, usecols=[2], names=['score'])
    score_table.index += 1
    return score_table.sort_values('score', ascending=False).index.values
cmd.extend('get_order_frame', get_order_frame)

def reorder_traj(dcd, score_file):
    ordered_frames = get_order_frame(score_file)
    print("[RNAPosers Debugging] New order of frames (score descending):", ordered_frames)
    dcd_new = "ordered_traj"
    cmd.delete(dcd_new)
    for frame in ordered_frames:
        combine_frame(dcd_new, dcd, frame=frame)
    # cmd.save_traj("ordered_traj.dcd", dcd_new)

cmd.extend('reorder_traj', reorder_traj)
