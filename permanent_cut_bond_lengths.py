from matplotlib.pyplot import *
import sys, os.path, argparse
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='What are the broken bonds and their on-graph distances. -1 indicates no path on graph. 2 a common neighbour.')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--dtmax', type=int, help='How long a bond must stay broken.', default=10)
    parser.add_argument('--noBrokenBonds', help='Whether to recompute instantaneously broken bonds', action='store_true')
    args = parser.parse_args()
    trname = args.trname
    dtmax =  args.dtmax
    print(trname)
    x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
    if not args.noBrokenBonds:
        print('broken bonds')
        pro = ProgressBar(x.size-1)
        for t, name in x.enum("_broken", "bonds"):
            if t==0: continue
            np.savetxt(name, np.array([(a,b) for a,b in x.broken_bonds(t-1,t)]), fmt='%d')
            pro.animate(t)
    print('bonds that stay broken for at least %d time steps'%dtmax)
    brokens = []
    pro = ProgressBar(x.size-1)
    for t, name in x.enum("_broken", "bonds"):
        if t==0: continue
        br = np.atleast_2d(np.loadtxt(name, dtype=int))
        if len(br) == 0 or len(br.ravel()) == 0:
            brokens.append(set())
        elif len(br.ravel()) == 2:
            brokens.append(set([(br.ravel()[0], br.ravel()[1])]))
        else:
            brokens.append(set((a,b) for a,b in br))
        if t==1: continue
        p2tr = x.p2tr(t)
        bonds = set((a,b) for a,b in np.sort(p2tr[np.loadtxt(x.get_format_string(ext='bonds')%t, dtype=int)], 1))
        for br in brokens[:-1]:
            br -= bonds
        if t>dtmax:
            br = brokens.pop(0)
            if len(br)==0:
                np.savetxt(
                    x.get_format_string("_broken_perm%d"%dtmax, "bonds")%(t-dtmax),
                    np.zeros((0,2), int),
                    fmt='%d')
            else:
                np.savetxt(
                    x.get_format_string("_broken_perm%d"%dtmax, "bonds")%(t-dtmax), 
                    np.array([(a,b) for a,b in br]), 
                    fmt='%d')
        pro.animate(t)
    
    print('Lengths of those bonds')
    pro = ProgressBar(x.size-dtmax)
    for t, name in x.enum("_broken_perm%d"%dtmax, "bonds"):
        if t==0 or t>=x.size-dtmax: continue
        p2tr = x.p2tr(t)
        g = nx.Graph()
        g.add_nodes_from(p2tr)
        g.add_edges_from(p2tr[np.loadtxt(x.get_format_string(ext="bonds")%t, dtype=int)])
        broken = np.atleast_2d(np.loadtxt(name, dtype=int))
        lengths = np.full(len(broken), -1, dtype=int)
        if len(broken) != 0 and len(broken.ravel()) != 0:
            for i, (a,b) in enumerate(broken):
                try:
                    lengths[i] = nx.shortest_path_length(g, a,b)
                except nx.NetworkXNoPath:
                    pass
        np.savetxt(
            x.get_format_string("_broken_perm%d"%dtmax, "length")%t, 
            lengths,
            fmt='%d')
        pro.animate(t)
