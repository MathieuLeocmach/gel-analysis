import sys, os.path, warnings, argparse
import numpy as np
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids import vtk, boo

if __name__ == '__main__':
    warnings.filterwarnings("error", "(.*)loadtxt: Empty input file:(.*)", UserWarning)
    parser = argparse.ArgumentParser(description='Compute bond orientational order fro l=2, particle centered and coarse-grained on bonds.')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--vtk', help='Whether to export VTK files showing the results(slow)', action='store_true')
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname) #'150A_Ageing_1303.traj')
    if args.vtk:
        print('export VTK')
    pro = ProgressBar(x.size)
    for t, name in x.enum():
        pos = np.loadtxt(name, skiprows=2)
        inside = np.min((pos-pos.min(0)>14) & (pos.max()-pos>14), -1)
        #np.min((pos>14) & ([256,256,128]-pos>14), -1)
        try:
            bonds = np.atleast_2d(np.loadtxt(x.get_format_string(ext='bonds')%t, int))
        except UserWarning:
            bonds = np.zeros([0,2], int)
        q2m = boo.bonds2qlm(pos, bonds, l=2)
        #set q2m of particles on the outside to 0
        q2m[np.logical_not(inside)] = 0
        q2 = boo.ql(q2m)
        np.savetxt(x.get_format_string(ext='q2')%t, q2, fmt='%.2f')
        q2b = boo.ql(q2m[bonds].sum(1)/np.maximum(1, inside[bonds].sum(1))[:,None])
        np.savetxt(x.get_format_string(ext='q2bonds')%t, q2b, fmt='%.2f')
        if args.vtk:
            v = vtk.Polydata()
            v.points = pos
            v.bonds = bonds
            v.scalars = [('q2', q2), ('inside', inside)]
            v.bondsScalars = [('q2', q2b)]
            v.save(x.get_format_string('_q2', ext='vtk')%t)
        pro.animate(t)
