from matplotlib.pyplot import *
import sys, os.path
import numpy as np
import networkx as nx
from colloids import vtk, experiment as xp
from colloids.progressbar import ProgressBar


    
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
    
def get_drift(pos0, pos1, p2tr0, p2tr1):
    tr2p0 = dict((tr, p) for p, tr in enumerate(p2tr0))
    tr2p1 = dict((tr, p) for p, tr in enumerate(p2tr1))
    both = np.intersect1d(p2tr0, p2tr1)
    drift = pos1[[tr2p1[tr] for tr in both]] - pos0[[tr2p0[tr] for tr in both]]
    return drift.mean(0)
    
    
if __name__ == '__main__':
    for trname in sys.argv[1:]:
        print(trname)
        x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
        maxlength = 0
        pro = ProgressBar(x.size)
        drift = np.zeros(3)
        #visualise bonds about to break between consecutive time steps
        for t, name in x.enum(ext='bonds'):
            try:
                bonds1 = np.atleast_2d(np.loadtxt(name, dtype=int))
            except UserWarning:
                bonds1 = np.zeros([0,2], int)
            p2tr1 = np.loadtxt(x.get_format_string(ext='p2tr')%t, dtype=int)
            pos1 = np.loadtxt(x.get_format_string()%t, skiprows=2)
            if t>0:
                #export all positions
                v = vtk.Polydata()
                v.points = pos0 - drift
                #export broken paths
                v.bonds = bonds0
                #label broken bonds
                isbreaking = breaking_bonds_p0(bonds0, bonds1, p2tr0, p2tr1)
                v.bondsScalars = [('breaking', isbreaking)]
                #label particles by how many broken bonds they participate in
                nbbreaking = np.zeros(len(pos0), int)
                nbbreaking[bonds0[isbreaking]] += 1
                v.scalars = [
                    ('nbbreaking', nbbreaking),
                    ('Nngb', np.histogram(bonds0.ravel(), np.arange(len(pos0)+1))[0])
                    ]
                #to disk
                v.save(x.get_format_string('_breakingbonds', 'vtk')%t)
                #update drift
                drift += get_drift(pos0, pos1, p2tr0, p2tr1)
            bonds0 = bonds1
            p2tr0 = p2tr1
            pos0 = pos1
            pro.animate(t)

