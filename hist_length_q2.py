#from matplotlib.pyplot import *
import sys, os.path, argparse, warnings
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
    
if __name__ == '__main__':
    warnings.filterwarnings("error", "(.*)loadtxt: Empty input file:(.*)", UserWarning)
    parser = argparse.ArgumentParser(description='Bin every bond according to its q2 value and its on-graph distances at next time, whether no path on graph, not broken, a common neighbour or a longer distance.')
    parser.add_argument('trname', help='Trajectory file name')
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
    pro = ProgressBar(x.size)
    q2bins = np.linspace(0,1,100)
    for t, name in x.enum("_q2length", "npy"):
        pro.animate(t)
        p2tr1 = x.p2tr(t)
        existin1 = np.zeros(x.nb_trajs, bool)
        existin1[p2tr1] = True
        #associate bonds at t0 with their value of q2
        try:
            bonds = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%t, int))
            q2s = np.loadtxt(x.get_format_string(ext='q2bonds')%t)
            bond2q2_1 = dict(
                ((a,b), q2) 
                for (a,b), q2 in zip(np.sort(p2tr1[bonds], 1), q2s)
                )
        except UserWarning:
            bond2q2_1 = dict()
        if t>0:
            #graph of the bonds between trajectories at t1
            g = nx.Graph()
            g.add_nodes_from(p2tr1)
            g.add_edges_from(bond2q2_1.keys())
            #which trajectories exist at both times step
            existinboth = existin0 & existin1
            #populate arrays before binning
            q2s0 = []
            Ls = []
            for (a,b), q2 in bond2q2_0.items():
                if existinboth[a] & existinboth[b]:
                    q2s0.append(q2)
                    try:
                        L = nx.shortest_path_length(g, a,b)
                    except nx.NetworkXNoPath:
                        L = -1
                    Ls.append(L)
            hist = np.histogram2d(q2s0, Ls, bins=(q2bins,[-1,1,2,3,1000]))[0].astype(int)
            np.save(name, hist)
        existin0 = existin1
        bond2q2_0 = bond2q2_1
        
