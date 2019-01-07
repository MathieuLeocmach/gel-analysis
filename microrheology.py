import sys, os.path, argparse
import numpy as np
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids import microrheology as mr
from scipy import constants
from multiprocessing import Pool

def get_t0s(size, start=0):
    """Logarythmically space initial times"""
    T = size - start - 1
    t0s = [10**i *j -1 for i in range(int(np.log10(T)+1)) for j in range(1,10) if 10**i *j < T]
    return np.array(t0s)

def get_t0s_avg(size, start=0, avg=10, base=10):
    """Logarythmically space initial times with time average"""
    T = size - start - 1 -avg/2
    t0s = np.array([base**i *j -1 for i in range(int(np.log(avg)/np.log(base)), int(np.log(T)/np.log(base)+1)) for j in range(1,base)])
    return t0s[(t0s<T) & (t0s>avg/2)]
    
def G_1P(tx, T=35, time_step=10., px2um=0.28, radius=2.75/2, width=0.7, alldt=True):
    """Compute pulsation-dependent shear moduli from 1P MSD at different t0. 
    Possible to use all dt to help smooth the fit."""
    t0s = get_t0s(tx.size)
    G1Ps = np.zeros((len(t0s),T), np.complex128)
    dts = np.arange(T)+1
    for it0, t0 in enumerate(t0s):
        if alldt:
            msd = np.sum((tx.positions[t0+1:] - tx.positions[t0])**2, -1).mean(-1)
        else:
            msd = np.sum((tx.positions[t0+1:][:T] - tx.positions[t0])**2, -1).mean(-1)
        omega, G = mr.msd2G(np.arange(1, len(msd)+1)*time_step, msd*px2um**2, radius, constants.C2K(28), width=width)
        G1Ps[it0,:min(T, len(G))] = G[:T]
    np.save(os.path.join(tx.xp.path, tx.xp.head)+'_moduli_1P.npy', G1Ps)
    omega = 1/(dts*time_step)
    return omega, G1Ps
    
def G_1P_avg(tx, T=35, avg=10, time_step=10., px2um=0.28, radius=2.75/2, width=0.7, alldt=True):
    """Compute pulsation-dependent shear moduli from 1P MSD at different t0. 
    Possible to use all dt to help smooth the fit. tx must have drift removed."""
    havg = max((int(avg/2),1))
    t0s = get_t0s_avg(tx.size, start=0, avg=avg, base=avg)
    G1Ps = np.zeros((len(t0s),T), np.complex128)
    dts = np.arange(T)+1
    for it0, t0 in enumerate(t0s):
        if alldt:
            msd = np.zeros(tx.size-t0-havg-2)
            for t00 in range(t0-havg, t0+havg+2):
                msd += np.sum((tx.positions[t00+1:][:len(msd)] - tx.positions[t00])**2, -1).mean(-1)
            msd /= 2*havg+1
        else:
            msd = np.zeros(T)
            for t00 in range(t0-havg, t0+havg+2):
                msd += np.sum((tx.positions[t00+1:][:T] - tx.positions[t00])**2, -1).mean(-1)
            msd /= 2*havg+1
        omega, G = mr.msd2G(np.arange(1, len(msd)+1)*time_step, msd*px2um**2, radius, constants.C2K(28), width=width)
        G1Ps[it0,:min(T, len(G))] = G[:T]
    np.save(os.path.join(tx.xp.path, tx.xp.head)+'_moduli_1P_avg%d.npy'%avg, G1Ps)
    omega = 1/(dts*time_step)
    return omega, G1Ps
    
    
def generate_2PBrownian(x, start=0):
    """Generate the Brownian components for logarythmically spaced initial times and all time lags"""
    tx = xp.Txp(x, start=start)

    m = tx.positions.min(1).min(0)
    M = tx.positions.max(1).max(0)
    L = (M-m).min()/4

    t0s = get_t0s(tx.size)
    insides = [
        np.where((tx.positions[t0] - m > L).min(-1) & (M - tx.positions[t0] > L).min(-1))[0] 
        for t0 in t0s
    ]
    tx.remove_drift()

    pro0 = ProgressBar(len(t0s)-1)
    for it0, t0 in enumerate(t0s):
        pro0.animate(it0)
        T = tx.positions.shape[0] - t0
        inside = insides[it0]
        #dts = [10**i *j for i in range(int(np.log10(T)+1)) for j in range(1,10) if 10**i *j < T]
        dts = np.arange(1, T)
        bres, nb, rbins = tx.get_Brownian_corr(t0, dts, 10, L, 10, inside=inside)
        np.save(
            os.path.join(tx.xp.path, tx.xp.head)+'_2msd_0t_%03d.npy'%t0,
            bres
            )
        np.save(
            os.path.join(tx.xp.path, tx.xp.head)+'_2nb_0t_%03d.npy'%t0,
            nb
            )
    return t0s, rbins, tx
    
def inner_2PBrownian(params):
    """Helper function for multithreaded"""
    tx, t0, inside, rmin, rmax, nbins = params
    T = tx.positions.shape[0] - t0
    #dts = [10**i *j for i in range(int(np.log10(T)+1)) for j in range(1,10) if 10**i *j < T]
    dts = np.arange(1, T)
    bres, nb, rbins = tx.get_Brownian_corr(t0, dts, rmin, rmax, nbins, inside=inside)
    np.save(
        os.path.join(tx.xp.path, tx.xp.head)+'_2msd_0t_%03d.npy'%t0,
        bres
        )
    np.save(
        os.path.join(tx.xp.path, tx.xp.head)+'_2nb_0t_%03d.npy'%t0,
        nb
        )
        
def generate_2PBrownian_multi(x, start=0, size=None, nbins=10, rmin=10., workers=2, chunksize=2):
    """Generate the Brownian components for logarythmically spaced initial times and all time lags"""
    tx = xp.Txp(x, start=start, size=size)

    m = tx.positions.min(1).min(0)
    M = tx.positions.max(1).max(0)
    L = (M-m).min()/4

    t0s = get_t0s(tx.size)
    insides = [
        np.where((tx.positions[t0] - m > L).min(-1) & (M - tx.positions[t0] > L).min(-1))[0] 
        for t0 in t0s
    ]
    tx.remove_drift()

    pro0 = ProgressBar(len(t0s)-1)
    with Pool(workers) as p:
        nn = len(t0s)
        for i, _ in enumerate(p.imap_unordered(
                inner_2PBrownian, 
                zip([tx]*nn, t0s, insides, [rmin]*nn, [L]*nn, [nbins]*nn),
                chunksize
            )):
            pro0.animate(i)
        
    lrmin, lrmax = np.log([rmin, L])
    lrbinsize = (lrmax - lrmin) / nbins
    lrbins = np.arange(nbins+1) * lrbinsize
    rbins = np.exp(lrbins+lrmin)
    
    return t0s, rbins, tx
    
def generate_2PBrownian_multi_avg(x, avg=10, start=0, size=None, nbins=10, rmin=10., workers=2, chunksize=2):
    """Generate the Brownian components for logarythmically spaced initial times and all time lags. Time average"""
    tx = xp.Txp(x, start=start, size=size)

    m = tx.positions.min(1).min(0)
    M = tx.positions.max(1).max(0)
    L = (M-m).min()/4

    t0s = get_t0s_avg(tx.size, start=0, avg=avg, base=avg)
    havg = max((1, int(avg/2)))
    
    #initial times on both sides of t0
    t0sw = np.unique(np.sort(np.ravel(t0s[:,None] + np.arange(-havg, havg+2)[None,:])))
    
    insides = [
        np.where((tx.positions[t0] - m > L).min(-1) & (M - tx.positions[t0] > L).min(-1))[0] 
        for t0 in t0sw
    ]
    tx.remove_drift()

    pro0 = ProgressBar(len(t0sw)-1)
    with Pool(workers) as p:
        nn = len(t0sw)
        for i, _ in enumerate(p.imap_unordered(
                inner_2PBrownian, 
                zip([tx]*nn, t0sw, insides, [rmin]*nn, [L]*nn, [nbins]*nn),
                chunksize
            )):
            pro0.animate(i)
    
    #sum the contributions on both sides of t0
    for t0 in t0s:
        nb = np.sum([
            np.load(os.path.join(tx.xp.path, tx.xp.head)+'_2nb_0t_%03d.npy'%t00)
            for t00 in range(t0-havg, t0+havg+2)
        ], axis=0)
        bress = [
            np.load(os.path.join(tx.xp.path, tx.xp.head)+'_2msd_0t_%03d.npy'%t00)
            for t00 in range(t0-havg, t0+havg+2)
        ]
        T = min(b.shape[1] for b in bress)
        bres = np.sum([b[:, :T] for b in bress], axis=0)
        np.save(
            os.path.join(tx.xp.path, tx.xp.head)+'_2msd_0t_%03d_avg%d.npy'%(t0, avg),
            bres
            )
        np.save(
            os.path.join(tx.xp.path, tx.xp.head)+'_2nb_0t_%03d_avg%d.npy'%(t0, avg),
            nb
            )
        
    lrmin, lrmax = np.log([rmin, L])
    lrbinsize = (lrmax - lrmin) / nbins
    lrbins = np.arange(nbins+1) * lrbinsize
    rbins = np.exp(lrbins+lrmin)
    
    return t0s, rbins, tx
    
    
def G_2P(tx, T=35, time_step=10., px2um=0.28, width_min=1., width_max=2.1, alldt=True):
    """Compute pulsation-dependent shear moduli from 2P MSD at different t0. 
    Possible to use all dt to help smooth the fit."""
    t0s = get_t0s(tx.size)
    G2Ps = np.zeros((len(t0s),T), np.complex128)
    dts = np.arange(T)+1
    a = tx.xp.rdf_radius()
    rbins = np.loadtxt(os.path.join(tx.xp.path, tx.xp.head)+'.rbins')
    midradii = np.sqrt(rbins[:-1]*rbins[1:])[2:9]
    for it0, t0 in enumerate(t0s):
        width = width_min + (width_max-width_min) * (it0 / len(t0s))
        bres = np.load(os.path.join(tx.xp.path, tx.xp.head)+'_2msd_0t_%03d.npy'%t0)
        nb = np.load(os.path.join(tx.xp.path, tx.xp.head)+'_2nb_0t_%03d.npy'%t0)
        msds, emsds = tx.MSD_twopoint(bres[...,2:9], nb[2:9])
        #weighted average along r
        msd = 2*(midradii * msds[0] * emsds[0]).sum(-1) / emsds[0].sum(-1)  / a
        if not alldt:
            msd = msd[:T]
        omega2, G = mr.msd2G(
            np.arange(1, len(msd)+1)*time_step, 
            msd* px2um**2, 
            a*px2um, constants.C2K(28), 
            width=width
        )
        G2Ps[it0,:min(T, len(G))] = G[:T]
    np.save(os.path.join(tx.xp.path, tx.xp.head)+'_moduli_2P.npy', G2Ps)
    omega = 1/(dts*time_step)
    return omega, G2Ps
    
    
def G_2P_avg(x, avg=10, start=0, size=None, T=35, time_step=10., px2um=0.28, width_min=1., width_max=2.1, alldt=True):
    """Compute pulsation-dependent shear moduli from time-averaged 2P MSD at different t0. 
    Possible to use all dt to help smooth the fit."""
    if size==None:
        size = x.size
    t0s = get_t0s_avg(size, start, avg, base=avg)
    G2Ps = np.zeros((len(t0s),T), np.complex128)
    dts = np.arange(T)+1
    a = x.rdf_radius()
    rbins = np.loadtxt(os.path.join(x.path, x.head)+'.rbins')
    midradii = np.sqrt(rbins[:-1]*rbins[1:])[2:9]
    for it0, t0 in enumerate(t0s):
        width = width_min + (width_max-width_min) * (it0 / len(t0s))
        bres = np.load(os.path.join(x.path, x.head)+'_2msd_0t_%03d_avg%d.npy'%(t0, avg))
        nb = np.load(os.path.join(x.path, x.head)+'_2nb_0t_%03d_avg%d.npy'%(t0, avg))
        msds, emsds = tx.MSD_twopoint(bres[...,2:9], nb[2:9])
        #weighted average along r
        msd = 2*(midradii * msds[0] * emsds[0]).sum(-1) / emsds[0].sum(-1)  / a
        if not alldt:
            msd = msd[:T]
        omega2, G = mr.msd2G(
            np.arange(1, len(msd)+1)*time_step, 
            msd* px2um**2, 
            a*px2um, constants.C2K(28), 
            width=width
        )
        G2Ps[it0,:min(T, len(G))] = G[:T]
    np.save(os.path.join(x.path, x.head)+'_moduli_2P_avg%d.npy'%avg, G2Ps)
    omega = 1/(dts*time_step)
    return omega, G2Ps

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export 1 and 2-particle microrheology')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--start', help='At which time to start computing', default=0, type=int)
    parser.add_argument('--size', help='Number of time steps after start (to remove late time drift).')
    parser.add_argument('--time_step', help='How many Brownian times per time steps', default=1., type=float)
    parser.add_argument('--avg', help='Do we average and by how many time steps?', type=int)
    parser.add_argument('--px2um', help='How many microns per pixel', default=0.28, type=float)
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname)
    start = args.start
    time_step = args.time_step
    if args.size:
        size = args.size
    else:
        size = x.size - start
    
    t0s, rbins, tx = generate_2PBrownian_multi(
        x, start=start, size=size, workers=10, chunksize=1
        )
    np.savetxt(
        os.path.join(x.path, x.head)+'.t0s',
        t0s, fmt='%d'
    )
    np.savetxt(
        os.path.join(x.path, x.head)+'.rbins',
        rbins,
    )
    T = tx.positions.shape[0] - t0s[-1]-1
    omega, G1Ps = G_1P(
        tx, T=T, time_step=10.*time_step, px2um=args.px2um, 
        radius=2.75/2, width=0.7, alldt=True
        )
    omega, G2Ps = G_2P(
        tx, T=T, time_step=10.*time_step, px2um=args.px2um, 
        width_min=1., width_max=2.1, alldt=True
        )
    omega, G1Psavg = G_1P_avg(
        tx, T=T, avg=10, time_step=10.*time_step, px2um=args.px2um, 
        radius=2.75/2, width=1., alldt=True
        )
    if args.avg:
        generate_2PBrownian_multi_avg(
            x, avg=args.avg, start=start, size=size, 
            workers=10, chunksize=1
            )
        G_2P_avg(
            x, size=size, avg=args.avg, T=T, time_step=time_step, px2um=args.px2um,
            width_min=1., width_max=2.1, alldt=True
            )
