import sys, os.path, argparse
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Size of the largest cluster of particles having at least N neighbours, with N from 1 to 7. Also corresponding cluster size distributions')
    thresholds = np.arange(7)
    
    parser.add_argument('trname', help='Trajectory file name')
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname)

    cluL = np.zeros((x.size, len(thresholds)), int)
    M = np.zeros(3)
    m = np.ones(3)*1024
    
    pro = ProgressBar(len(cluL.ravel()))
    for t, name in x.enum(ext='bonds'):
        pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
        m = np.minimum(m, pos.min(0))
        M = np.maximum(M, pos.max(0))
        bonds = np.loadtxt(name, int)
        nngb = np.histogram(bonds.ravel(), np.arange(len(pos)+1))[0]
        g = nx.Graph()
        g.add_nodes_from(np.arange(len(pos)))
        g.add_edges_from(bonds)
        for i, thrsq in enumerate(thresholds):
            if np.sum(nngb>thrsq)==0:
                np.savetxt(
                    x.get_format_string('_Nngb%d'%(thrsq+1),'sRg')%t,
                    np.zeros((0,3)),
                    header='s Rg Rg_err\n'
                )
                continue
            gn = g.subgraph(np.where(nngb>thrsq)[0])
            cc = [list(clu) for clu in nx.connected_components(gn)]
            L = [pos[clu].ptp(0).max() for clu in cc]
            if len(L)>0:
                cluL[t,i] = max(L)
            #cluster size distribution
            ss = np.zeros(len(cc), int)
            Rgs = np.zeros(len(cc))
            for ic, clu in enumerate(cc):
                ss[ic] = len(clu)
                Rgs[ic] = np.sqrt(np.var(pos[clu], axis=0).sum())
            nbs, bins = np.histogram(ss, bins=np.arange(ss.max())+1)
            sumRg = np.histogram(ss, bins, weights=Rgs)[0]
            sumRg2 = np.histogram(ss, bins, weights=Rgs**2)[0]
            meanRg = sumRg[nbs>0]/nbs[nbs>0]
            stdRg = np.sqrt(sumRg2[nbs>0]/nbs[nbs>0] - meanRg**2)
            #count both the error on values and on number of clusters at a given s
            errRg = stdRg + 5/nbs[nbs>0]
            np.savetxt(
                x.get_format_string('_Nngb%d'%(thrsq+1),'sRg')%t,
                np.column_stack((
                        bins[:-1][nbs>0], 
                        meanRg, 
                        errRg
                        )),
                header='s Rg Rg_err\n'
                )
            pro.animate(t*len(thresholds) + i)
    np.save(os.path.join(x.path, 'perco_Nngb.npy'), cluL)
    np.savetxt(
        os.path.join(x.path, '%s_perco_Nngb.txt'%(os.path.splitext(x.trajfile)[0])), 
        np.column_stack((np.arange(x.size), cluL/(M-m).max())), 
        header=' '.join(['t']+['N_C%d'%i for i in thresholds+1])
        )
