from matplotlib import pyplot as plt
import sys, os.path, argparse, warnings
import numpy as np
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot and export the breaking probability function of q2 depending on the on-graph distances at next time, whether no path on graph, a common neighbour or a longer distance.')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--tmin', type=int, help='When to start averaging.', default=1)
    parser.add_argument('--tmax', help='When to stop averaging.', default=None)
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
    if args.tmax is None:
        args.tmax = x.size
    else:
        args.tmax = int(args.tmax)
    pro = ProgressBar(x.size)
    q2bins = np.linspace(0,1,100)
    hist = np.zeros((len(q2bins)-1, 4), int)
    for t, name in x.enum("_q2length", "npy"):
        if t<args.tmin or t >= args.tmax:
            continue 
        hist += np.load(name).astype(int)
    htot = np.maximum(1, hist.sum(-1))
    plt.figure('probabroken_q2_type')
    plt.clf()
    for l, h in zip(['common neighbour', 'long path', 'no path'], hist.T[[2,3,0],1:-1]):
        plt.plot(q2bins[1:-2][h>10], (h/htot[1:-1])[h>10], label=l)
        np.savetxt(
            os.path.join(
                x.path, 
                x.head+'_probabroken_q2_%s_from%03d_to%03d.txt'%(
                    '_'.join(l.split()),
                    args.tmin, args.tmax
                    )
                ),
            np.column_stack((q2bins[1:-2], h/htot[1:-1]))[h>10],
            header='q2 proba\n'
            )
    plt.xlabel('$q_2$')
    plt.ylabel(r'$P_\mathrm{break}$')
    plt.yscale('log')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(
        x.path, 
        x.head+'_probabroken_q2_type_from%03d_to%03d.pdf'%(args.tmin, args.tmax)
        ))
