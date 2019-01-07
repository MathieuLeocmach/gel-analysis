import sys, os.path, argparse
import numpy as np
#import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids.povray import *
import subprocess, os
from matplotlib.pyplot import imread, imsave

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export POVray renderings of slices of configurations showing particles colored by isostaticity.')
        
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--zmin', help='Z lower value', default=100.)
    parser.add_argument('--zmax', help='Z higher value', default=125.)
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname)
    #prepare file name pattern to export in current directory
    exportname = os.path.split(x.get_format_string('_isostatic2_slice', 'pov'))[1]
        
    pro = ProgressBar(x.size)
    for t, name in x.enum(ext='bonds'):
        pro.animate(t)
        pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
        #define slice
        good = (pos[:,-1] > args.zmin) & (pos[:,-1] < args.zmax)
        #number of neighbours per particles
        bonds = np.loadtxt(name, dtype=int)
        nngb = np.zeros(len(pos), dtype=int)
        np.add.at(nngb, bonds.ravel(), 1)
        #export to POVray format
        povname = exportname%t
        f = File(povname, "colors.inc", "header_slices.inc")
        isostatic = Union(*[
            Sphere( (x,y,z), 5)
            for x,y,z in pos[(nngb >= 6) & good]
            ]+[Texture(Pigment(color=(0.694,  0.145,  0.623)))])
        f.write(isostatic)
        nonisostatic = Union(*[
                Sphere( (x,y,z), 3)
                for x,y,z in pos[(nngb < 6) & good]
                ]+[Texture(Pigment(color=(1.,  1,  1)))])
        f.write(nonisostatic)
        bondPOV = Union(*[
            Cylinder(tuple(a.tolist()), tuple(b.tolist()), 1)
            for a, b in pos[bonds[good[bonds].min(-1)]]
            ]+[Texture(Pigment(color=(1.,  0.455,  0.156)))])
        f.write(bondPOV)
        f.file.flush()
        #rendering
        jpgname = os.path.splitext(povname)[0]+'.jpg'
        pngname =  os.path.splitext(povname)[0]+'.png'
        subprocess.call(['povray', '-O%s'%pngname, '+W512', '-D00', '+H512' , povname])
        imsave(jpgname, imread(pngname), format='jpg')
        os.remove(pngname)
        f.file.flush()
        #break
#ffmpeg -r 25 -f image2 -s 512x512 -i 163A_1340_percolation_isostatic_slice_t%03d.jpg -vcodec libx264 -crf 25  -pix_fmt yuv420p -vf "drawtext=fontfile=Arial.ttf: fontsize=32 : text=%{eif\\\:n-71\\\\:d\\\\:3}: x=(w-tw)/2: y=(2*lh): fontcolor=black: box=1" isostatic.mp4
