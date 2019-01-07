from matplotlib.pyplot import *
import sys, os.path, warnings
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

def shortest_length_values(g, trs, tr2v):
    lengths = np.zeros([len(trs), 3])
    for i, (a,b) in enumerate(trs):
        if not g.has_node(a) or not g.has_node(b): continue
        try:
            lengths[i,0] = nx.shortest_path_length(g, a,b)
        except nx.NetworkXNoPath:
            lengths[i,0] = -1
        lengths[i,1] = tr2v[a]
        lengths[i,2] = tr2v[b]
    return lengths
    
def broken_bonds_lenghts_values(bonds0, bonds1, p2tr0, p2tr1, field=None):
    #bonds (between trajectories) existing at t but no more at t+dt
    # = broken bonds + lost trajectories
    trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])

    #graph of the bonds between trajectories at t+dt
    g = nx.Graph()
    g.add_nodes_from(p2tr1) #([(p, {'v':v}) for p,v in zip(p2tr0, field)])
    g.add_edges_from(p2tr1[bonds1])
    
    #values of the field mapped on trajectories existing at t
    tr2v = dict()
    for tr, v in zip(p2tr0, field):
        tr2v[tr] = v

    return shortest_length_values(g, trbonds, tr2v)
    
if __name__ == '__main__':
    warnings.filterwarnings("error", "(.*)loadtxt: Empty input file:(.*)", UserWarning)
    for trname in sys.argv[1:]:
        print trname
        x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
        maxlength = 0
        pro = ProgressBar(x.size)

        for t, name in x.enum(ext='bonds'):
            try:
                bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
            except UserWarning:
                bonds1 = np.zeros([0,2], int)
            p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
            if t>0:
                q2 = np.loadtxt(x.get_format_string(ext='q2')%(t-1))
                lengths = broken_bonds_lenghts_values(bonds0, bonds1, p2tr0, p2tr1, q2)
                np.savetxt(x.get_format_string('_broken_q2', ext='length')%t, lengths, fmt='%g')
                if len(lengths)>0:
                    maxlength = max([maxlength, lengths.max()])
            bonds0 = bonds1
            p2tr0 = p2tr1
            pro.animate(t)
    
        histvalues = np.zeros(maxlength+1)
        histlength = np.zeros(maxlength+1, int)
        for t, name in x.enum('_broken_q2', ext='length'):
            if t==0: continue
            try:
                lengths = np.atleast_2d(np.loadtxt(name))
            except UserWarning:
                continue
            #length distribution of broken bonds whoes particles 
            #have non zero values of the field
            histlength += np.histogram(
                lengths[:,0], 
                weights=np.sum(lengths[:,1:]>0, -1), 
                bins=np.arange(-1, maxlength+1)
                )[0]
            #length distribution weighted by the values of the field
            histvalues += np.histogram(
                lengths[:,0], 
                weights=lengths[:,1:].sum(-1), 
                bins=np.arange(-1, maxlength+1)
                )[0]
        np.savetxt(
            os.path.join(x.path, 'broken_length_q2.hist'),
            np.column_stack((
                np.arange(-1, maxlength), 
                histlength, 
                histvalues,
                histvalues/np.maximum(1, histlength)
                )),
            fmt='%g')
        #plot(np.arange(-1, maxlength), histvalues/np.maximum(1, histlength))

#plot(np.arange(-1, lengths.max()), np.histogram(lengths, bins=np.arange(-1, lengths.max()+1))[0])
