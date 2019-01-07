import sys, os.path, warnings
import numpy as np
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids import boo
from colloids.particles import bonds2ngbs

if __name__ == '__main__':
    warnings.filterwarnings("error", "(.*)loadtxt: Empty input file:(.*)", UserWarning)
    for trname in sys.argv[1:]:
        print trname
        x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
        pro = ProgressBar(x.size)
        hq2b = np.zeros([99], int)
        for t, name in x.enum():
            pos = np.loadtxt(name, skiprows=2)
            inside = np.min((pos-pos.min(0)>14) & (pos.max()-pos>14), -1)
            try:
                bonds = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%t, int))
                ngbs = bonds2ngbs(bonds, len(pos))
            except UserWarning:
                bonds = np.zeros([0,2], int)
                ngbs = np.zeros([len(pos),1], int)
            q2m = boo.weave_qlm(pos, ngbs, inside, l=2)
            q2 = boo.ql(q2m[bonds].sum(1)/np.maximum(1, inside[bonds].sum(1))[:,None])
            hq2b += np.histogram(q2, bins= np.linspace(0,1,100))[0]
            pro.animate(t)
        
        np.savetxt(os.path.join(x.path, 'q2b.hist'), np.column_stack((np.linspace(0,1,100)[:-1], hq2b)))
