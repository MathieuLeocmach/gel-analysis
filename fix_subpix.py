import sys, os.path, argparse, shutil
import numpy as np
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix subpixel resolution. Fixed data are copied preserving the directory tree from the current directory to a new root.')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--newroot', help='Where to root the new directory tree.', default='../0_Data_fixed/')
    parser.add_argument('--coef', help='Correction factor', default=0.215, type=float)
    args = parser.parse_args()
    trname = args.trname
    n = os.path.splitext(trname)[0]
    print(trname)
    x = xp.Experiment(trname)
    
    newdir = os.path.join(args.newroot, os.path.split(n)[0])
    try:
        os.makedirs(newdir)
    except FileExistsError:
        pass
    
    shutil.copy(n+'.traj', os.path.join(args.newroot, n+'.traj'))
    shutil.move(n+'.rdf', os.path.join(args.newroot, n+'.rdf'))
    x2 = xp.Experiment(os.path.join(args.newroot, n+'.traj'))
    
    pro = ProgressBar(x.size)
    for t,name in x.enum():
        with open(x2.get_format_string()%t, 'wb') as f:
            for i, line in enumerate(open(name, 'rb')):
                if i ==2: break
                f.write(line)
            pos = np.loadtxt(name, skiprows=2); pos[:,-1] /= 1.03709
            xe = pos - np.round(pos)
            xe = xe*np.abs(xe)/args.coef**2*0.5
            pos = np.round(pos)+xe; pos[:,-1] *= 1.03709
            np.savetxt(f, pos)
        pro.animate(t)

