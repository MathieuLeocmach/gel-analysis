import sys, os.path, warnings
import numpy as np
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids import particles

if __name__ == '__main__':
    for trname in sys.argv[1:]:
        print(trname)
        x = xp.Experiment(trname)
        pro = ProgressBar(x.size)
        
        #bounding box
        m = []
        M = []
        for t, name in x.enum():
            pos = np.loadtxt(name, skiprows=2)
            m.append(pos.min(0))
            M.append(pos.max(0))
            pro.animate(t)
        m = np.min(m,0)
        M = np.max(M, 0)
        
        #Compute structure factor at each time
        Ss = np.zeros([x.size, 128])
        Ss6 = np.zeros_like(Ss)
        for t, name in x.enum():
            pos = np.loadtxt(name, skiprows=2) - m
            Ss[t] = particles.structure_factor(pos, 128, M-m, maxNvec=300)
            bonds = np.loadtxt(x.get_format_string(ext='bonds')%t, int)
            nngb = np.histogram(bonds.ravel(), np.arange(len(pos)+1))[0]
            isostatic = nngb>5
            if np.any(isostatic):
                Ss6[t] = particles.structure_factor(pos[isostatic], 128, M-m, maxNvec=300)
            pro.animate(t)
        #save in 0th value the value of q1, so that the q axis is q1*np.arange(len(S))
        Ss[:,0] = 2*np.pi/np.max(M-m)
        Ss6[:,0] = Ss[:,0]
        np.save(
            os.path.join(x.path, os.path.splitext(x.trajfile)[0]+'_structure_factor.npy'), 
            Ss
            )
        np.save(
            os.path.join(x.path, os.path.splitext(x.trajfile)[0]+'_structure_factor_NC6.npy'), 
            Ss6
            )
        #use the last Ss (a priori the most developped) to get the first and second mimimum
        fm = 1 + np.argmin(Ss[-1,1:])
        fM = fm + np.argmax(Ss[-1,fm:])
        sm = fM + np.argmin(Ss[-1,fM:])
        qmax = (Ss[:,fm:sm]*np.arange(fm,sm)).sum(1)/Ss[:,fm:sm].sum(1)
        peak_value = Ss[:,fm:sm].sum(1)
        np.savetxt(
            os.path.join(x.path, os.path.splitext(x.trajfile)[0]+'.qmax'), 
            np.column_stack((qmax * Ss[0,0], peak_value))
            )
        peak_value6 = Ss6[:,fm:sm].sum(1)
        nonzero = peak_value6>0
        qmax6 = np.zeros(len(nonzero))
        qmax6[nonzero] = (Ss6[nonzero,fm:sm]*np.arange(fm,sm)).sum(1)/Ss6[nonzero,fm:sm].sum(1)
        np.savetxt(
            os.path.join(x.path, os.path.splitext(x.trajfile)[0]+'.qmaxNC6'), 
            np.column_stack((qmax6 * Ss6[0,0], peak_value6))
            )
        
