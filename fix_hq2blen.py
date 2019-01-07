import sys, os.path, warnings
import numpy as np
import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids import boo

if __name__ == '__main__':
    warnings.filterwarnings("error", "(.*)loadtxt: Empty input file:(.*)", UserWarning)
    for trname in sys.argv[1:]:
        print trname
        x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
        maxlength = np.loadtxt(os.path.join(x.path, 'broken_length_q2b.hist')).shape[0]
        hq2blen = np.zeros([maxlength+1, 99], np.int32)
        pro = ProgressBar(x.size)
        for t, name in x.enum('_broken_q2b', ext='length'):
            if t==0: continue
            try:
                lengths = np.atleast_2d(np.loadtxt(name))
                if lengths.shape[-1]>0:
                    hq2blen += np.histogram2d(lengths[:,0], lengths[:,1], bins=(np.arange(-1, maxlength+1), np.linspace(0,1,100)))[0]
            except UserWarning: continue
            pro.animate(t)
        np.save(os.path.join(x.path, 'broken_length_q2b.npy'), hq2blen)
