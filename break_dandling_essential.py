from matplotlib.pyplot import *
import sys, os.path
import numpy as np
import networkx as nx
from colloids import vtk, experiment as xp
from colloids.progressbar import ProgressBar
import argparse

    
def breaking_bonds_p0(bonds0, bonds1, p2tr0, p2tr1):
    """Label bonds at t0 between trajectories no more bounded t1.
    
    bonds0, bonds1 are respectively the bonds at t0 and 1 in terms of position
    p2tr0, p2tr1 are respectively the position to trajectory relationship at t0 and t1
    """
    #reversed indexes trajectory to position
    tr2p0 = dict((tr,p) for p, tr in enumerate(p2tr0))
    tr2p1 = dict((tr,p) for p, tr in enumerate(p2tr1))
    #bonds (between trajectories) existing at t but no more at t+dt
    trbonds = set([(a,b) for a,b in np.sort(p2tr0[bonds0], 1)]) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])
    
    breaking0 = set([
            tuple(sorted([tr2p0[tr] for tr in trbond]))
            for trbond in trbonds
            #filter out bonds involving at least a particle that do not exist at t1
            if trbond[0] in tr2p1 and trbond[1] in tr2p1
            ])
    isbreaking = np.zeros(len(bonds0), bool)
    for i, (a,b) in enumerate(bonds0):
        isbreaking[i] = (a,b) in breaking0
    return isbreaking
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Number of bonds, broken bonds, dandling bonds, breaking dandling bonds, essential bonds, breaking essential bonds.')
    parser.add_argument("trajfile", help="path to the file.traj")
    args = parser.parse_args()
    trname = args.trajfile
    print(trname)
    x = xp.Experiment(trname)
    maxlength = 0
    pro = ProgressBar(x.size)
    with open(os.path.join(x.path, x.head+'_dandling_essential.broken'), 'w') as f:
        f.write('#bonds broken dandling broken_dandling essential broken_essential\n')
        for t, name in x.enum(ext='bonds'):
            try:
                bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
                pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
                bondlength1 = np.sqrt((np.diff(pos1[bonds1], axis=1)[:,0]**2).sum(-1))
                bonds1 = bonds1[bondlength1 >= 9]
                bondlength1 = bondlength1[bondlength1 >= 9]
            except UserWarning:
                bonds1 = np.zeros([0,2], int)
                bondlength = np.zeros([0,1], float)
            p2tr1 = x.p2tr(t)
            if t>0:
                #graph, neighbourhood
                g0 = nx.Graph()
                g0.add_nodes_from(range(len(p2tr0)))
                g0.add_edges_from(bonds0)
                nngb = np.histogram(bonds0, np.arange(len(p2tr0)+1))[0]
                
                #which trajectories exists at t1?
                exist1 = np.zeros(len(x.trajs), bool)
                exist1[p2tr1] = True
                #select bonds at t0 where both particles exist at t1
                existb0 = exist1[p2tr0[bonds0]].min(-1)
                bonds0ex1 = bonds0[existb0]
                #further conditionds are restrained to the bonds where 
                #both particles exists at t1
                
                #select breaking bonds
                isbreaking = breaking_bonds_p0(bonds0ex1, bonds1, p2tr0, p2tr1)
                np.savetxt(
                    x.get_format_string('_breaking', 'lengths')%t,
                    bondlength0[existb0][isbreaking]
                    )
                
                #what are the dandling bonds at t0? At least one particle is a dead end
                isdandling = (nngb[bonds0ex1]==1).max(-1)
                
                #what are the essential bonds at t0? No common neighbour, not dandling.
                isessential = np.logical_not(isdandling)
                for i, (a,b) in enumerate(bonds0ex1):
                    isessential[i] &= sum(1 for ngb in nx.common_neighbors(g0, a, b))==0
                
                #to disk
                f.write(' '.join(
                    '%d'%n for n in [
                        len(bonds0ex1), isbreaking.sum(),
                        isdandling.sum(), np.logical_and(isdandling, isbreaking).sum(), 
                        isessential.sum(), np.logical_and(isessential, isbreaking).sum()]
                        )+'\n')
            
            bonds0 = bonds1
            p2tr0 = p2tr1
            bondlength0 = bondlength1
            pro.animate(t)

