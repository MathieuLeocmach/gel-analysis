import sys, os.path, argparse
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='How many times each bond breaks, existence time.')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--dtmax', type=int, help='History length.', default=10)
    args = parser.parse_args()
    trname = args.trname
    dtmax =  args.dtmax
    print(trname)
    x = xp.Experiment(trname)
    pro = ProgressBar(x.size)
    history = []
    broken = []
    for t, name in x.enum(ext='bonds'):
        history.append(set(
            (a,b) 
            for a,b in np.sort(x.p2tr(t)[
                np.loadtxt(name, dtype=int)
            ], 1)
            ))
        if t < 1: continue
        broken.append(set(
            (a,b) 
            for a,b in np.loadtxt(x.get_format_string('_broken', 'bonds')%t)
            ))
        if t<dtmax: continue
        bonds = history[-1]
        if len(bonds)==0:
            np.savetxt(name, np.zeros((0,4), int), fmt='%d')
        else:
            appearance = dict()
            for ago, his in enumerate(history[::-1]):
                for bond in bonds & his:
                    appearance[bond] = ago
            nbbreaks = {bond:0 for bond in bonds}
            for h,g in zip(history[:-2], history[1:-1]):
                for bond in (h&bonds) - g:
                    nbbreaks[bond] += 1
            np.savetxt(
                x.get_format_string("_rate%d"%dtmax, ext='bonds')%t, 
                np.array([(a,b, ap, nbbreaks[(a,b)]) for (a,b), ap in appearance.items()]), 
                fmt='%d')
        del history[0]
        pro.animate(t)
