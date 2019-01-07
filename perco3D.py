import sys, os.path, argparse
import numpy as np
#import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids.povray import *
import subprocess, os
from matplotlib.pyplot import imread, imsave

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export POVray renderings of configurations showing particles colored by cluster size.')
    thresholds = np.arange(1, 7)
    
    parser.add_argument('trname', help='Trajectory file name')
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname)
    
    print('count maximum number of clusters')
    maxsize = max(
        max(len(line[:-1].split(',')) for line in open(name))
        for t, name in x.enum('_c2p', ext='csv')
        )-1
    print(maxsize)
    
    pro = ProgressBar(x.size)
    for t, name in x.enum('_c2p', ext='csv'):
        pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
        povname = x.get_format_string('_clusersize', 'pov')%t
        jpgname = x.get_format_string('_clusersize', 'jpg')%t
        pngname =  x.get_format_string('_clusersize', 'png')%t
        f = File(povname,"colors.inc", "perco.inc")
        for line in open(name):
            c = np.array(list(map(int, line[:-1].split(','))))
            if len(c)>1:
                cluster = Union(*[
                    Sphere( (x,y,z), 5)
                    for x,y,z in pos[c]
                    ]+[Texture(Pigment(color='COLORSCALE(%g)'%((len(c)-1)/maxsize)))])
            else:
                cluster = [
                    Sphere(
                        (x,y,z), 5,
                        Texture(Pigment(color='COLORSCALE(0)'))
                        )
                    for x,y,z in pos[c]
                    ][0]
            f.write(cluster)
        f.file.flush()
        
        subprocess.call(['povray', '-O%s'%pngname, '+W256', '+H256' , povname])
        imsave(jpgname, imread(pngname))
        os.remove(pngname)
