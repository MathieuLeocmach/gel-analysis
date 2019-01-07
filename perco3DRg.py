import sys, os.path, argparse
import numpy as np
#import networkx as nx
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids.povray import *
import subprocess, os
from matplotlib.pyplot import imread, imsave

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export POVray renderings of configurations showing particles colored by cluster radius of gyration.')
    
    parser.add_argument('trname', help='Trajectory file name')
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname)
    
    exportname = os.path.split(x.get_format_string('_p2Rg', ext='pov'))[1]
    
    print('find maximum cluster Rg')
    maxsize = max(
        np.loadtxt(name).max()
        for t, name in x.enum('_p2rg', ext='csv')
        )-1
    print(maxsize)
    
    pro = ProgressBar(x.size)
    for t, name in x.enum('_p2rg', ext='csv'):
        pos = np.loadtxt(x.get_format_string()%t, skiprows=2)
        povname = exportname%t
        jpgname = os.path.splitext(povname)[0]+'.jpg'
        pngname = os.path.splitext(povname)[0]+'.png'
        f = File(povname,"colors.inc", "perco.inc")
        for ((X,Y,Z), c) in zip(pos, np.loadtxt(name)/maxsize):
            f.write(
                Sphere(
                    (X,Y,Z), 5,
                    Texture(Pigment(color='COLORSCALE(%g)'%c))
                    )
            )
        f.file.flush()
        
        subprocess.call(['povray', '-O%s'%pngname, '+W512', '+H512' , povname])
        imsave(jpgname, imread(pngname), format='jpg')
        os.remove(pngname)
