from matplotlib.pyplot import *
import sys, os.path
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

def shortest_length(g, trs):
    lengths = np.zeros(len(trs), int)
    for i, (a,b) in enumerate(trs):
        if not g.has_node(a) or not g.has_node(b): continue
        try:
            lengths[i] = nx.shortest_path_length(g, a,b)
        except nx.NetworkXNoPath:
            lengths[i] = -1
    return lengths
    
def broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1):
    #bonds (between trajectories) existing at t but no more at t+dt
    # = broken bonds + lost trajectories
    trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])

    #graph of the bonds between trajectories at t+dt
    g = nx.Graph()
    g.add_nodes_from(p2tr1)
    g.add_edges_from(p2tr1[bonds1])

    return shortest_length(g, trbonds)
    
if __name__ == '__main__':
    for trname in sys.argv[1:]:
        print(trname)
        x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
        maxlength = 0
        pro = ProgressBar(x.size)

        for t, name in x.enum(ext='bonds'):
            bonds1 = np.loadtxt(name, dtype=int)
            p2tr1 = x.p2tr(t)
            if t>0:
                lengths = broken_bonds_lenghts(bonds0, bonds1, p2tr0, p2tr1)
                np.savetxt(x.get_format_string('_broken', ext='length')%t, lengths, fmt='%d')
                if len(lengths)>0:
                    maxlength = max([maxlength, lengths.max()])
            bonds0 = bonds1
            p2tr0 = p2tr1
            pro.animate(t)
    
        histlength = np.sum([
            np.histogram(np.loadtxt(name, dtype=int), bins=np.arange(-1, maxlength+1))[0]
            for t, name in x.enum('_broken', ext='length')
            if t>0], axis=0)
        np.savetxt(
            os.path.join(x.path, 'broken_length.hist'),
            np.column_stack((np.arange(-1, maxlength), histlength)),
            fmt='%d')


#plot(np.arange(-1, lengths.max()), np.histogram(lengths, bins=np.arange(-1, lengths.max()+1))[0])
