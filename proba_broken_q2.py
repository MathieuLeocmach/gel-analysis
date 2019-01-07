import sys, os.path, warnings
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids import boo

def bonds2ngbs(bonds, N):
    """Returns an array of neighbours from bond data"""
    ngbs = -np.ones([N, np.histogram(bonds, bins=np.arange(N+1))[0].max()], int)
    if bonds.shape[-1]>0:
        for a,b in bonds:
            ngbs[a, np.where(ngbs[a]==-1)[0][0]] = b
            ngbs[b, np.where(ngbs[b]==-1)[0][0]] = a
    return ngbs

q2bins = np.linspace(0,1,100)

if __name__ == '__main__':
    warnings.filterwarnings("error", "(.*)loadtxt: Empty input file:(.*)", UserWarning)
    for trname in sys.argv[1:]:
        print trname
        x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
        maxlength = 0
        pro = ProgressBar(x.size)
        
        N_bonds = np.zeros([x.size-1, len(q2bins)-1], int)
        N_broken = np.zeros([x.size-1, len(q2bins)-1], int)

        for t, name in x.enum(ext='bonds'):
            try:
                bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
                if bonds1.shape[-1]==0:
                    bonds1 = np.zeros([0, 2], int)
            except UserWarning:
                bonds1 = np.zeros([0,2], int)
            p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
            if t>0:
                #the bonds0 have been loaded at previous time
                #compute the q2m for each particle
                pos = np.loadtxt(x.get_format_string()%(t-1), skiprows=2)
                inside = np.min((pos-pos.min(0)>14) & (pos.max()-pos>14), -1)
                ngbs = bonds2ngbs(bonds0, len(pos))
                q2m = boo.weave_qlm(pos, ngbs, inside, l=2)
                #trajectories that exist at t but not at t+dt
                terminal_tr = set(p2tr0) - set(p2tr1)
                #positions than do not belong to terminal trajectories
                notterminal = np.ones(len(pos), bool)
                for p, tr in enumerate(p2tr0):
                    notterminal[p] = tr not in terminal_tr
                #the bonds to keep are the one with both members inside and not terminal
                good_p = inside & notterminal
                goodbonds = bonds0[good_p[bonds0].min(-1)]
                #compute q2 for each good bond by averaging the q2m of the members
                q2 = boo.ql(q2m[goodbonds].mean(1))
                #count function of q2
                N_bonds[t-1] = np.histogram(q2, bins=q2bins)[0]
                #convert the goodbonds into trbonds0 while keeping the link with q2
                trbonds0 = {
                    tuple(np.sort(p2tr0[[a,b]]).tolist()): v
                    for (a,b), v in zip(goodbonds, q2)
                    }
                #broken bonds
                broken = set(trbonds0.keys()) - set([(a,b) for a,b in np.sort(p2tr1[bonds1], axis=1)])
                #graph of the bonds between trajectories at t+dt
                g = nx.Graph()
                g.add_nodes_from(p2tr1)
                g.add_edges_from(p2tr1[bonds1])
                #The new distance on graph should be larger than 2 (further than neighbour of neighbour)
                realbroken = []
                for a,b in broken:
                    try:
                        if nx.shortest_path_length(g, a,b)>2:
                            realbroken.append((a,b))
                    except nx.NetworkXNoPath:
                        #if there is no way to go from a to b, add it also
                        realbroken.append((a,b))
                #count function of q2
                N_broken[t-1] = np.histogram([trbonds0[b] for b in realbroken], bins=q2bins)[0]
            #prepare next step
            bonds0 = bonds1
            p2tr0 = p2tr1
            pro.animate(t)
        np.save(
            os.path.join(x.path, x.trajfile[:-5]+'_hist_bonds_broken_q2.npy'), 
            np.dstack((N_bonds, N_broken))
            )
