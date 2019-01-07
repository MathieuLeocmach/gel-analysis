import sys, os.path, argparse
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
    
def directed_extant(X, bonds):
    """For each particlem the maximum spatial extant of a all directed paths going through it.
    
    X : spatial coordinate of the nodes along the direction of interest
    bonds : bonds between nodes
    """
    #create a directed graph where the edges are directed along X
    gD = nx.DiGraph()
    gD.add_nodes_from(np.arange(len(X)))
    gD.add_edges_from(
        b[s]
        for b,s in zip(bonds, np.argsort(X[bonds], axis=-1))
        )
    #depth first search
    Ls = np.zeros(len(X), int)
    #iterate on nodes sorted by X
    for v in np.argsort(X):
        if Ls[v]>0: continue
        descendants = np.array(list(nx.dfs_preorder_nodes(gD, source=v)), int)
        L = np.max(X[descendants]-X[v])
        Ls[descendants] = np.maximum(Ls[descendants], L)
    return Ls
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Spatial extant of the largest isostatic cluster directed along X, Y and Z. For each particle the maximal spatial extent of any directed path going through it.')
    
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--tmax', type=int, default=None, help='When to stop')
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname)
    
    if args.tmax is None or args.tmax > x.size:
        args.tmax = x.size

    cluL = np.zeros((args.tmax, 3), int)
    M = np.zeros(3)
    m = np.ones(3)*1024
    
    pro = ProgressBar(len(cluL.ravel()))
    for t, name in x.enum(ext='bonds'):
        if t >= len(cluL):
            break
        pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
        m = np.minimum(m, pos.min(0))
        M = np.maximum(M, pos.max(0))
        bonds = np.loadtxt(name, int)
        nngb = np.histogram(bonds.ravel(), np.arange(len(pos)+1))[0]
        #restrict to bonds involving at two isostatic particles
        bonds = bonds[np.min(nngb[bonds] > 5, axis=-1)]
        for dim in range(pos.shape[1]):
            Ls = directed_extant(pos[:,dim], bonds)
            np.save(x.get_format_string('_perco_isostatic_%sdirected'%('XYZ'[dim]), 'npy')%t, Ls)
            cluL[t, dim] = Ls.max()
            pro.animate(t*3 + dim)
    np.save(os.path.join(x.path, 'perco_isostatic_directed.npy'), cluL)
    np.savetxt(
        os.path.join(x.path, '%s_perco_isostatic_directed.txt'%(os.path.splitext(x.trajfile)[0])), 
        np.column_stack((np.arange(len(cluL)), cluL.max(-1)/(M-m).max())), 
        header=' '.join(['tXYZ'])
        )
