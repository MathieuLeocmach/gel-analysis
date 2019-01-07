import sys, os.path, argparse
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='What are the stable bonds.')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--dtmax', type=int, help='How long a bond must exist to be labeled stable.', default=10)
    args = parser.parse_args()
    trname = args.trname
    dtmax =  args.dtmax
    print(trname)
    x = xp.Experiment(trname)
    pro = ProgressBar(x.size)
    history = []
    for t, name in x.enum("_stable%d"%dtmax,ext='bonds'):
        history.append(set(
            (a,b) 
            for a,b in np.sort(x.p2tr(t)[
                np.loadtxt(x.get_format_string(ext='bonds')%t, dtype=int)
            ], 1)
            ))
        if t<dtmax: continue
        br = history.pop(0).intersection(*history)
        if len(br)==0:
                np.savetxt(name, np.zeros((0,2), int), fmt='%d')
        else:
            np.savetxt(name, np.array([(a,b) for a,b in br]), fmt='%d')
        pro.animate(t)
        
