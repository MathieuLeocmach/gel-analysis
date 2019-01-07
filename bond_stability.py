import sys, os.path, argparse
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='How many time steps is each bond existing.')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--dtmax', type=int, help='History length.', default=10)
    args = parser.parse_args()
    trname = args.trname
    dtmax =  args.dtmax
    print(trname)
    x = xp.Experiment(trname)
    pro = ProgressBar(x.size)
    history = []
    for t, name in x.enum(ext='bonds'):
        history.append(set(
            (a,b) 
            for a,b in np.sort(x.p2tr(t)[
                np.loadtxt(name, dtype=int)
            ], 1)
            ))
        if t<dtmax: continue
        bonds = history[-1]
        if len(bonds)==0:
            np.savetxt(name, np.zeros((0,3), int), fmt='%d')
        else:
            stability = {bond:1 for bond in bonds}
            for his in history[:-1]:
                for bond in bonds & his:
                    stability[bond] += 1
            np.savetxt(
                x.get_format_string("_stability%d"%dtmax, ext='bonds')%t, 
                np.array([(a,b, st) for (a,b), st in stability.items()]), 
                fmt='%d')
        del history[0]
        pro.animate(t)
