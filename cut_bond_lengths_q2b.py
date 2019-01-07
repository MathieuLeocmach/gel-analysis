from matplotlib.pyplot import *
import sys, os.path, warnings
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids import boo

def shortest_length_values(g, trs, tr2v):
    lengths = np.zeros([len(trs), 2])
    for i, (a,b) in enumerate(trs):
        if not g.has_node(a) or not g.has_node(b): continue
        try:
            lengths[i,0] = nx.shortest_path_length(g, a,b)
        except nx.NetworkXNoPath:
            lengths[i,0] = -1
        if boo.ql(np.array([tr2v[a]])) == 0:
            lengths[i,1] = boo.ql(np.array([tr2v[b]]))
        elif boo.ql(np.array([tr2v[b]])) == 0:
            lengths[i,1] = boo.ql(np.array([tr2v[a]]))
        else:
            lengths[i,1] = boo.ql(np.array([(tr2v[a] + tr2v[b])/2]))
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
    
def bonds2ngbs(bonds, N):
    """Returns an array of neighbours from bond data"""
    ngbs = -np.ones([N, np.histogram(bonds, bins=np.arange(N+1))[0].max()], int)
    if bonds.shape[-1]>0:
        for a,b in bonds:
            ngbs[a, np.where(ngbs[a]==-1)[0][0]] = b
            ngbs[b, np.where(ngbs[b]==-1)[0][0]] = a
    return ngbs
    
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
                if bonds1.shape[-1]==0:
                    bonds1 = np.zeros([0, 2], int)
            except UserWarning:
                bonds1 = np.zeros([0,2], int)
            p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
            if t>0:
                pos = np.loadtxt(x.get_format_string()%(t-1), skiprows=2)
                inside = np.min((pos-pos.min(0)>14) & (pos.max()-pos>14), -1)
                ngbs = bonds2ngbs(bonds0, len(pos))
                q2m = boo.weave_qlm(pos, ngbs, inside, l=2)
                #q2 = np.loadtxt(x.get_format_string(ext='q2')%(t-1))
                lengths = broken_bonds_lenghts_values(bonds0, bonds1, p2tr0, p2tr1, q2m)
                np.savetxt(x.get_format_string('_broken_q2b', ext='length')%t, lengths, fmt='%g')
                if len(lengths)>0:
                    maxlength = max([maxlength, lengths[:,0].astype(int).max()])
            bonds0 = bonds1
            p2tr0 = p2tr1
            pro.animate(t)
    
        histvalues = np.zeros(maxlength+1)
        histlength = np.zeros(maxlength+1, int)
        hq2blen = np.zeros([maxlength+1, 99], np.int32)
        for t, name in x.enum('_broken_q2b', ext='length'):
            if t==0: continue
            try:
                lengths = np.atleast_2d(np.loadtxt(name))
            except UserWarning:
                continue
            #length distribution of broken bonds whoes particles 
            #have non zero values of the field
            histlength += np.histogram(
                lengths[:,0], 
                weights=lengths[:,1]>0, 
                bins=np.arange(-1, maxlength+1)
                )[0]
            #length distribution weighted by the values of the field
            histvalues += np.histogram(
                lengths[:,0], 
                weights=lengths[:,1], 
                bins=np.arange(-1, maxlength+1)
                )[0]
            #2d histogram length vs q2
            hq2blen += np.histogram2d(
                lengths[:,0], lengths[:,1], 
                bins=(np.arange(-1, maxlength+1), np.linspace(0,1,100))
                )[0]
        np.savetxt(
            os.path.join(x.path, 'broken_length_q2b.hist'),
            np.column_stack((
                np.arange(-1, maxlength), 
                histlength, 
                histvalues,
                histvalues/np.maximum(1, histlength)
                )),
            fmt='%g', header='Lambda histlength histvalues proba\n')
        #plot(np.arange(-1, maxlength), histvalues/np.maximum(1, histlength))
        np.save(os.path.join(x.path, 'broken_length_q2b.npy'), hq2blen)
        

#plot(np.arange(-1, lengths.max()), np.histogram(lengths, bins=np.arange(-1, lengths.max()+1))[0])
