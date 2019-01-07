import sys, os.path, argparse
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Time evolution of the average coordination number, size of the largest cluster divided by largest box size, effective volume fraction.')
    
    parser.add_argument('trname', help='Trajectory file name')
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname)

    #generate filename_c2p_t???.csv if not existing
    if not os.path.isfile(x.get_format_string('_c2p', 'csv')%0):
        nbparts = x.get_nb()
        bondpattern = x.get_format_string(ext='bonds')
        pro = ProgressBar(x.size)
        for t, name in x.enum('_c2p', 'csv'):
            pro.animate(t)
            g = nx.Graph()
            g.add_nodes_from(np.arange(nbparts[t]))
            bonds = np.atleast_2d(np.loadtxt(bondpattern%t, dtype=int))
            g.add_edges_from(bonds)
            with open(name, 'w') as f:
                for cluster in nx.connected_components(g):
                    f.write(','.join('%d'%p for p in cluster)+'\n')
    #time evolution of the average coordination number
    ncs = x.get_nb_bonds()[:x.size]*2./x.get_nb()[:x.size]
    #time evolution of the maximum extent of the largest culster
    maxptp = []
    #compute the bounding box in the same loop
    m = np.ones(3)*1024
    M = np.zeros(3)
    #compute the sum of radii of gyrations cubed
    Rg3 = []
    pro = ProgressBar(x.size)
    for t, name in x.enum('_c2p', 'csv'):
        pro.animate(t)
        parts = np.loadtxt(x.get_format_string()%t, skiprows=2)
        m = np.minimum(m, parts.min(0))
        M = np.maximum(M, parts.max(0))
        maxptp.append(max([
            parts[[int(e) for e in line.split(',')]].ptp(0).max() 
            for line in open(name)
            ]))
        Rg3.append(np.power(x.rdf_radius()+np.sqrt([
            parts[[int(e) for e in line.split(',')]].var(0).sum(-1)
            for line in open(name)
            ]), 3).sum())
    #normalize the maximum extent of the largest culster by the accessible extent
    maxptp = np.array(maxptp)/(M-m).max()
    #time evolution of the effective volume fraction
    phieff = np.array(Rg3) * 4/3. * np.pi / np.prod(M-m)
    #save the result
    np.savetxt(
        os.path.join(x.path, x.head+'.nc_ptp'), 
        np.column_stack((ncs, maxptp, phieff)),
        header='N_C perco phieff\n'
        )
