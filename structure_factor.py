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
        m = np.min(m,0)
        M = np.max(M, 0)
        
        #Compute structure factor at each time
        Ss = np.zeros([x.size, 128])
        for t, name in x.enum():
            Ss[t] = particles.structure_factor(np.loadtxt(name, skiprows=2)-m, 128, M-m, maxNvec=300)
            pro.animate(t)
        #save in 0th value the value of q1, so that the q axis is q1*np.arange(len(S))
        Ss[:,0] = 2*np.pi/np.max(M-m)
        np.save(
            os.path.join(x.path, os.path.splitext(x.trajfile)[0]+'_structure_factor.npy'), 
            Ss
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
        
